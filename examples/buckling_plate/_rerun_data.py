"""
Sean Engelstad
Feb 2024, GT SMDO Lab
Goal is to rerun certain data points if the mesh wasn't converged

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

# select the load style and BC (4 total combinations)
# need to generate all 4 combinations of data to finish this
loading = "Nx"  # "Nx", "Nxy"
BC = "SS"  # "SS", "CL"

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

inner_ct = 0

# read the original csv file
csv_file = os.path.join(data_folder, csv_filename)
df = pd.read_csv(csv_file)
_Dstar = df["Dstar"].to_numpy()
_AR = df["a/b"].to_numpy()
_SR = df["b/h"].to_numpy()
_materials = df["material"].to_numpy()
_ply_angle = df["ply_angle"].to_numpy()
_nx = df["nx"].to_numpy()
_ny = df["ny"].to_numpy()

# short test here
# df.loc[0,"kmin"] = 2.0
# val = df["kmin"][0]
# print(f"df row 0 = {val}")

# exit()

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

N = _materials.shape[0]
#N = 1
for row in range(N):
    # current criterion
    if _AR[row] < 1.0:
        pass # redo this data point (low AR)
    else:
        continue

    print(f"\n\n Re-Running row {row} with AR = {_AR[row]}")

    # randomly generate the material
    material = buckling_surrogate.FlatPlateAnalysis.get_material_from_str(_materials[row])
    ply_angle = _ply_angle[row]

    # random geometry, min thickness so that K,G matrices have good norm
    slenderness = _SR[row]
    h = 1.0
    b = (
        h * slenderness
    )  # verified that different b values don't influence non-dim buckling load

    fail_ct = 0
    aspect_ratio = _AR[row]
    a = aspect_ratio * b

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

    else:  # just do a model parameter check
        kmin = 1.0  # for model parameter check

    reasonable_min = 0.0 < kmin < 100.0

    if abs(error_0) < 1e-10 and reasonable_min and kmin:
        # perform the mode tracking
        tracked_eigvals, _ = buckling_surrogate.FlatPlateAnalysis.mac_permutation(
            nominal_plate, new_plate, num_modes=20
        )

        # update the model parameters in the original dataframe
        df.loc[row,"kmin"] = np.real(kmin)
        df.loc[row,"error"] = np.real(error_0)
        df.loc[row,"nx"] = nx
        df.loc[row,"ny"] = ny
        df.loc[row,"Dstar"] = new_plate.Dstar
        df.loc[row,"a0/b0"] = new_plate.affine_aspect_ratio
        df.loc[row,"a/b"] = new_plate.aspect_ratio
        df.loc[row,"material"] = new_plate.material_name
        df.loc[row,"ply_angle"] = new_plate.ply_angle

        # add the tracked eigenvalues
        for imode, tracked_eigval in enumerate(tracked_eigvals):
            df.loc[row,f"k_{imode+1}"] = np.real(tracked_eigval) if tracked_eigval else None

        # rewrite the data frame in there
        df.to_csv(csv_file, index=False, mode="w")