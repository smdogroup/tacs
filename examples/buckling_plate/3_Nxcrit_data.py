"""
Sean Engelstad
Feb 2024, GT SMDO Lab
Goal is to generate data for pure uniaxial compression failure in the x-direction
For now just simply supported only..
"""

from tacs import buckling_surrogate
import numpy as np
import pandas as pd
import os, niceplots
import matplotlib.pyplot as plt

from mpi4py import MPI

comm = MPI.COMM_WORLD

# model inputs
N = 10000

# start of computations
data_dict_list = []

# make a folder for the pandas data
cpath = os.path.dirname(__file__)
data_folder = os.path.join(cpath, "data")
if not os.path.exists(data_folder) and comm.rank == 0:
    os.mkdir(data_folder)

# TODO : change this to a Monte carlo simulation over uniform scales or uniform log scales
# TODO : use the affine transformation of 1_run_analysis.py to make lambda values more consistent.

# arrays to check the values of

# clear the csv file
_clear_data = True

csv_file = os.path.join(data_folder, "Nxcrit.csv")
if _clear_data:
    if os.path.exists(csv_file) and comm.rank == 0:
        os.remove(csv_file)

# maybe it would be better to use a materials database and select a random material then
# and use material classmethods instead

# TODO : rewrite this by randomly choose an isotropic / composite material
# randomly choose a b and h pair until AR and affine_AR in range
# this way a lot of random slenderness values will be used
# iterate by sweeping AR through the data

inner_ct = 0

for foo in range(100):  # until has generated this many samples
    # randomly generate the material
    materials = buckling_surrogate.FlatPlateAnalysis.get_materials()
    material = np.random.choice(np.array(materials))
    ply_angle = np.random.uniform(0.0, 90.0)

    # random geometry, min thickness so that K,G matrices have good norm
    log_slenderness = np.random.uniform(np.log(5.0), np.log(200.0))
    slenderness = np.exp(log_slenderness)
    h = 0.1
    b = h * slenderness  # verified that different b values don't influence non-dim buckling load

    for aspect_ratio in np.linspace(0.2, 5.0, 40):
        a = aspect_ratio * b

        # create the flat plate analysis

        # make the flat plate
        flat_plate: buckling_surrogate.FlatPlateAnalysis = material(
            comm,
            bdf_file="plate.bdf",
            a=a,
            b=b,
            h=h,
        )

        # make sure the affine aspect ratio is in a reasonable range
        _accepted = 0.05 <= flat_plate.affine_aspect_ratio <= 20.0

        if not (_accepted):
            continue  # go to next iteration

        # select number of elements
        # in order to preserve element AR based on overall AR
        min_elem = 30
        max_elem = 150
        if aspect_ratio > 1.0:
            nx = np.min([int(aspect_ratio * min_elem), max_elem])
            ny = min_elem
        else:  # AR < 1.0
            ny = np.min([int(min_elem / aspect_ratio), max_elem])
            nx = min_elem

        _run_buckling = True

        if _run_buckling:
            load_scale = 0.5

            flat_plate.generate_bdf(
                nx=nx,  # my earlier mistake was the #elements was not copied from above!!
                ny=ny,
                exx=flat_plate.affine_exx
                * load_scale,  # scale down to make sure in pre-buckling
                eyy=0.0,
                exy=0.0,
                clamped=False,
            )

            # avg_stresses = flat_plate.run_static_analysis(write_soln=True)
            # if comm.rank == 0:
            #    print(f"avg stresses = {avg_stresses}")

            # Sx0 = avg_stresses[0]
            # Sy0 = avg_stresses[1]
            # Sxy0 = avg_stresses[2]

            tacs_eigvals, errors = flat_plate.run_buckling_analysis(
                sigma=10.0 / load_scale, num_eig=12, write_soln=False
            )

            kx_0 = tacs_eigvals[0] * load_scale
            error_0 = errors[0]

        else:  # just do a model parameter check
            kx_0 = 1.0  # for model parameter check

        reasonable_kx0 = 0.0 < kx_0 < 100.0

        if abs(error_0) < 1e-10 and reasonable_kx0:
            # record the model parameters
            data_dict = {
                # model parameter section
                "Dstar": [flat_plate.Dstar],
                "a0/b0": [flat_plate.affine_aspect_ratio],
                "a/b" : [flat_plate.aspect_ratio],
                "b/h": [flat_plate.slenderness],
                "kx_0": [kx_0],
                "error": [error_0],
                # other parameter section
                "material": [flat_plate.material_name],
                "ply_angle": [flat_plate.ply_angle],
                "nx": [nx],
                "ny": [ny],
            }

            # "s_xx" : [Sx0],
            # "s_yy" : [Sy0],
            # "s_xy" : [Sxy0],

            data_dict_list += [data_dict]

            # write to the csv file
            # convert to a pandas dataframe and save it in the data folder
            inner_ct += 1 # started from 0, so first time is 1
            if comm.rank == 0:
                df = pd.DataFrame(data_dict)
                if inner_ct == 1 and not (os.path.exists(csv_file)):
                    df.to_csv(csv_file, mode="w", index=False, header=True)
                else:
                    df.to_csv(csv_file, mode="a", index=False, header=False)

            # MPI COMM Barrier in case running with multiple procs
            comm.Barrier()

        else:
            continue

# report the percentage of good models out of the Monte Carlo simulation
# frac_good_models = inner * 1.0 / ct
# print(f"fraction of good models: {accepted_ct} / {ct} = {frac_good_models}")

# check the model parameter ranges were covered
# _check_params = False
# plt.style.use(niceplots.get_style())
# if _check_params and comm.rank == 0:
#     plt.figure("Dstar-AR")
#     plt.plot(
#         [_["Dstar"] for _ in data_dict_list], [_["a0/b0"] for _ in data_dict_list], "ko"
#     )
#     plt.yscale("log")
#     plt.xlabel("Dstar")
#     plt.ylabel("AR")
#     plt.show()
#     plt.close("Dstar-AR")

#     plt.figure("Dstar-slenderR")
#     plt.plot(
#         [_["Dstar"] for _ in data_dict_list], [_["b/h"] for _ in data_dict_list], "ko"
#     )
#     plt.yscale("log")
#     plt.xlabel("Dstar")
#     plt.ylabel("slenderR")
#     plt.show()
#     plt.close("Dstar-slenderR")

# # plot aspect ratio versus kx0
# plt.plot([_["a0/b0"] for _ in data_dict_list], [_["kx_0"] for _ in data_dict_list], "o")
# plt.savefig("data/Dstar-kx0.png", dpi=400)
