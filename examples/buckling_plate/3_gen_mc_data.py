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
import argparse
from mpi4py import MPI

comm = MPI.COMM_WORLD

# parse the arguments
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--load', type=str)
parent_parser.add_argument('--BC', type=str)

args = parent_parser.parse_args()

assert args.load in ["Nx", "Nxy", "axial", "shear"]
assert args.BC in ["SS", "CL"]

if args.load in ["Nx", "axial"]:
    loading = "Nx"
else:
    loading = "Nxy"
BC = args.BC

# MODEL INPUTS SECTION
# --------------------------------------------
# --------------------------------------------

# number of random samples (traverse all AR for each of them)
N = 100

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

# run the nominal plate for mode tracking
nominal_plate = buckling_surrogate.FlatPlateAnalysis.hexcelIM7(
    comm=comm,
    bdf_file="plate.bdf",
    a=1.0,
    b=1.0,
    h=0.01,
    ply_angle=0,
)

nominal_plate.generate_bdf(
    nx=30,
    ny=30,
    exx=nominal_plate.affine_exx if loading == "Nx" else 0.0,
    eyy=0.0,
    exy=nominal_plate.affine_exy if loading == "Nxy" else 0.0,
    clamped=BC == "CL",
)
nom_eigvals, _ = nominal_plate.run_buckling_analysis(
    sigma=5.0, num_eig=20, write_soln=False
)

for foo in range(N):  # until has generated this many samples
    # randomly generate the material
    materials = buckling_surrogate.FlatPlateAnalysis.get_materials()
    material = np.random.choice(np.array(materials))
    ply_angle = np.random.uniform(0.0, 90.0)

    # random geometry, min thickness so that K,G matrices have good norm
    log_slenderness = np.random.uniform(np.log(10.0), np.log(200.0))
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
        new_plate: buckling_surrogate.FlatPlateAnalysis = material(
            comm,
            bdf_file="plate.bdf",
            a=a,
            b=b,
            h=h,
            ply_angle=ply_angle,
        )

        # make sure the affine aspect ratio is in a reasonable range
        _accepted = 0.05 <= new_plate.affine_aspect_ratio <= 20.0

        if not (_accepted):
            continue  # go to next iteration

        # select number of elements
        # in order to preserve element AR based on overall AR
        _nelems = 1000
        AR_g1 = aspect_ratio if aspect_ratio > 1 else 1.0/aspect_ratio
        min_elem = int(np.sqrt(_nelems / AR_g1))
        max_elem = int(min_elem * AR_g1)
        if aspect_ratio > 1.0:
            nx = max_elem
            ny = min_elem
        else:  # AR < 1.0
            ny = max_elem
            nx = max(min_elem, 25)

        _run_buckling = True

        if _run_buckling:
            if loading == "Nx":
                exx = new_plate.affine_exx
                exy = 0.0
            elif loading == "Nxy":
                exx = 0.0
                exy = new_plate.affine_exy

            clamped = BC == "CL"

            new_plate.generate_bdf(
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

            new_eigvals, errors = new_plate.run_buckling_analysis(
                sigma=5.0, num_eig=40, write_soln=False
            )

            # min eigenvalue
            kmin = new_eigvals[0]
            error_0 = errors[0]

            if loading == "Nxy": # can have negative lowest eigenvalue (+- for each)
                kmin = np.abs(kmin)

        else:  # just do a model parameter check
            kmin = 1.0  # for model parameter check

        reasonable_min = 0.0 < kmin < 100.0

        if abs(error_0) < 1e-10 and reasonable_min and kmin:
            # perform the mode tracking
            tracked_eigvals, _ = buckling_surrogate.FlatPlateAnalysis.mac_permutation(
                nominal_plate, new_plate, num_modes=20
            )

            # record the model parameters
            data_dict = {
                # model parameter section
                "Dstar": [new_plate.Dstar],
                "a0/b0": [new_plate.affine_aspect_ratio],
                "a/b": [new_plate.aspect_ratio],
                "b/h": [new_plate.slenderness],
                "kmin": [np.real(kmin)],
                "error": [np.real(error_0)],
                # other parameter section
                "material": [new_plate.material_name],
                "ply_angle": [new_plate.ply_angle],
                "nx": [nx],
                "ny": [ny],
            }

            # add the tracked eigenvalues
            for imode, tracked_eigval in enumerate(tracked_eigvals):
                data_dict[f"k_{imode+1}"] = (
                    np.real(tracked_eigval) if tracked_eigval else None
                )

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
