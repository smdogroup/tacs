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
_clear_data = False

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

accepted_ct = 0
ct = 0
while accepted_ct < N:  # until has generated this many samples
    # randomly generate the material constants
    # TODO : change this to randomly selecting a material from a material database
    log_Emax = np.log10(200e9)
    log_E11 = np.random.uniform(9, log_Emax)
    E11 = 10**log_E11
    log_E22 = np.random.uniform(9, log_Emax)
    E22 = 10**log_E22
    log_G12 = np.random.uniform(9, log_Emax)
    G12 = 10**log_G12
    nu12 = np.random.uniform(
        0.1, 0.4
    )  # was [-0.95, 0.45] before which meant a lot of the data was wacky (non-realistic with nu12<0)

    # randomly generate the plate sizes
    # TODO : change this to selecting affine_AR => computing AR
    # select slenderness, scale b such that h > 0.01 (since below that is usually bad for eigenvalue solver)
    log_a = np.random.uniform(-1, 1)
    a = 10**log_a
    log_b = np.random.uniform(-1, 1)
    b = 10**log_b
    log_h = np.random.uniform(-3, -2)
    h = 10**log_h

    # make the flat plate
    flat_plate = buckling_surrogate.FlatPlateAnalysis(
        comm,
        bdf_file="plate.bdf",
        a=a,
        b=b,
        h=h,
        E11=E11,
        nu12=nu12,
        E22=E22,  # set to None if isotropic
        G12=G12,  # set to None if isotropic
    )

    # get main model parameters from the flat plate class
    Dstar = flat_plate.Dstar
    a0_b0 = flat_plate.affine_aspect_ratio
    log_a0_b0 = np.log10(a0_b0)
    slenderR = flat_plate.slenderness
    a_b = flat_plate.aspect_ratio
    AR = a_b
    log_slender = np.log10(slenderR)

    # check the model parameter ranges are reasonable
    valid_Dstar = 0 <= Dstar <= 1.0
    valid_a_b = 0.05 <= a_b <= 10.0
    valid_a0_b0 = 0.05 <= a0_b0 <= 20.0
    valid_bh = 5 <= slenderR <= 100  # maybe I should allow higher slenderness ratios?

    # skip this random model if model parameters are outside ranges
    ct += 1
    if valid_Dstar and valid_a0_b0 and valid_a_b and valid_bh:
        pass
    else:
        _print_out_of_bounds = False
        if _print_out_of_bounds:
            if not (valid_Dstar):
                print(f"invalid Dstar = {Dstar}")
            if not (valid_a_b):
                print(f"invalid a/b = {a_b}")
            if not (valid_a0_b0):
                print(f"invalid a0/b0 = {a0_b0}")
            if not (valid_bh):
                print(f"invalid b/h = {slenderR}")
        continue

    # if model parameters were in range then we can now run the buckling analysis

    # select number of elements
    min_elem = 30
    max_elem = 120
    if AR > 1.0:
        nx = np.min([int(AR * min_elem), max_elem])
        ny = min_elem
    else:  # AR < 1.0
        ny = np.min([int(min_elem / AR), max_elem])
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

    if 1.0 < kx_0 < 100.0 and abs(error_0) < 1e-10:
        accepted_ct += 1

        # log the model parameters
        data_dict = {
            # model parameter section
            "Dstar": [Dstar],
            "a0/b0": [a0_b0],
            # "a/b" : [AR],
            "b/h": [slenderR],
            "kx_0": [kx_0],
            "error": [error_0],
            # "s_xx" : [Sx0],
            # "s_yy" : [Sy0],
            # "s_xy" : [Sxy0],
            # other parameter section
            "E11": [E11],
            "E22": [E22],
            "nu12": [nu12],
            "G12": [G12],
            "a": [a],
            "b": [b],
            "h": [h],
            "nx": [nx],
            "ny": [ny],
        }

        data_dict_list += [data_dict]

        # write to the csv file
        # convert to a pandas dataframe and save it in the data folder
        if comm.rank == 0:
            df = pd.DataFrame(data_dict)
            if accepted_ct == 1 and not (os.path.exists(csv_file)):
                df.to_csv(csv_file, mode="w", index=False)
            else:
                df.to_csv(csv_file, mode="a", index=False, header=False)

        # MPI COMM Barrier in case running with multiple procs
        comm.Barrier()

    else:
        continue

# report the percentage of good models out of the Monte Carlo simulation
frac_good_models = accepted_ct * 1.0 / ct
print(f"fraction of good models: {accepted_ct} / {ct} = {frac_good_models}")

# check the model parameter ranges were covered
_check_params = False
plt.style.use(niceplots.get_style())
if _check_params and comm.rank == 0:
    plt.figure("Dstar-AR")
    plt.plot(
        [_["Dstar"] for _ in data_dict_list], [_["a0/b0"] for _ in data_dict_list], "ko"
    )
    plt.yscale("log")
    plt.xlabel("Dstar")
    plt.ylabel("AR")
    plt.show()
    plt.close("Dstar-AR")

    plt.figure("Dstar-slenderR")
    plt.plot(
        [_["Dstar"] for _ in data_dict_list], [_["b/h"] for _ in data_dict_list], "ko"
    )
    plt.yscale("log")
    plt.xlabel("Dstar")
    plt.ylabel("slenderR")
    plt.show()
    plt.close("Dstar-slenderR")

# plot aspect ratio versus kx0
plt.plot([_["a0/b0"] for _ in data_dict_list], [_["kx_0"] for _ in data_dict_list], "o")
plt.savefig("data/Dstar-kx0.png", dpi=400)
