"""
Sean Engelstad
Feb 2024, GT SMDO Lab
Goal is to generate data for pure uniaxial compression failure in the x-direction
For now just simply supported only..

gen_mc_data.py : generate monte carlo training data for each of the surrogate models
"""

from tacs import buckling_surrogate
import numpy as np
import pandas as pd
import os

from mpi4py import MPI

comm = MPI.COMM_WORLD

# MODEL INPUTS SECTION
# --------------------------------------------
# --------------------------------------------

# number of random samples (traverse all AR for each of them)
N = 100

# select the load style and BC (4 total combinations)
# need to generate all 4 combinations of data to finish this
loading = "Nx"  # "Nx", "Nxy"
BC = "clamped"  # "SS", "clamped"

# END OF MODEL INPUTS SECTION
# --------------------------------------------
# --------------------------------------------

csv_filename = loading + "crit_" + BC + ".csv"

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
_clear_data = False

csv_file = os.path.join(data_folder, csv_filename)
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

for foo in range(N):  # until has generated this many samples
    # randomly generate the material
    materials = buckling_surrogate.FlatPlateAnalysis.get_materials()
    material = np.random.choice(np.array(materials))
    ply_angle = np.random.uniform(0.0, 90.0)

    # random geometry, min thickness so that K,G matrices have good norm
    log_slenderness = np.random.uniform(np.log(5.0), np.log(200.0))
    slenderness = np.exp(log_slenderness)
    h = 1.0
    b = (
        h * slenderness
    )  # verified that different b values don't influence non-dim buckling load

    fail_ct = 0

    for aspect_ratio in np.linspace(0.2, 5.0, 40):
        a = aspect_ratio * b

        if fail_ct > 5:
            break  # exit out of this loop as the solver is failing here..
            # may want to print out to a file when it does this (so the user can know this happened)

        # create the flat plate analysis

        # make the flat plate
        flat_plate: buckling_surrogate.FlatPlateAnalysis = material(
            comm,
            bdf_file="plate.bdf",
            a=a,
            b=b,
            h=h,
            ply_angle=ply_angle,
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
            load_scale = 1.0

            if loading == "Nx":
                exx = flat_plate.affine_exx * load_scale
                exy = 0.0
            elif loading == "Nxy":
                exx = 0.0
                exy = flat_plate.affine_exy * load_scale

            clamped = BC == "clamped"

            flat_plate.generate_bdf(
                nx=nx,  # my earlier mistake was the #elements was not copied from above!!
                ny=ny,
                exx=exx,
                eyy=0.0,
                exy=exy,
                clamped=clamped,
            )

            # avg_stresses = flat_plate.run_static_analysis(write_soln=True)
            # if comm.rank == 0:
            #    print(f"avg stresses = {avg_stresses}")

            # Sx0 = avg_stresses[0]
            # Sy0 = avg_stresses[1]
            # Sxy0 = avg_stresses[2]

            tacs_eigvals, errors = flat_plate.run_buckling_analysis(
                sigma=5.0 / load_scale, num_eig=12, write_soln=False
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
                "a/b": [flat_plate.aspect_ratio],
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
            inner_ct += 1  # started from 0, so first time is 1
            if comm.rank == 0:
                df = pd.DataFrame(data_dict)
                if inner_ct == 1 and not (os.path.exists(csv_file)):
                    df.to_csv(csv_file, mode="w", index=False, header=True)
                else:
                    df.to_csv(csv_file, mode="a", index=False, header=False)

            # MPI COMM Barrier in case running with multiple procs
            comm.Barrier()

        else:
            fail_ct += 1
            continue
