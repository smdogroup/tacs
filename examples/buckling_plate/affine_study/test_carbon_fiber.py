"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from tacs import buckling_surrogate
import numpy as np
from mpi4py import MPI
import matplotlib.pyplot as plt
import pandas as pd
import niceplots, os
from mpi4py import MPI

comm = MPI.COMM_WORLD

# vary aspect ratio and for D* = 0.352 for carbon fiber
# show that the equations match the affine transformation from
# "Generic Buckling Curves For Specially Orthotropic Rectangular Plates" by Brunelle

AR_list = []
kx0_FEA_list = []  # FEA
kx0_CF_list = []  # closed-form

# clear the csv file
csv_file = "axial-SS2.csv"
if os.path.exists(csv_file) and comm.rank == 0:
    os.remove(csv_file)

ct = 0
for AR in np.linspace(0.5, 5, 40):
    AR_list += [AR]

    # carbon fiber material with D* = 0.352

    b = 5.0  # try different b values and make sure curve still valid
    slenderness = 100.0  # 33.3

    flat_plate = buckling_surrogate.FlatPlateAnalysis(
        comm=comm,
        bdf_file="plate.bdf",
        a=AR * b,
        b=b,
        h=b / slenderness,
        E11=135e9,
        nu12=0.3,
        E22=10e9,  # set to None if isotropic
        G12=5e9,  # set to None if isotropic
    )

    # select number of elements
    min_elem = 30
    max_elem = 150
    if AR > 1.0:
        nx = np.min([int(AR * min_elem), max_elem])
        ny = min_elem
    else:  # AR < 1.0
        ny = np.min([int(min_elem / AR), max_elem])
        nx = min_elem

    flat_plate.generate_bdf(
        nx=nx,
        ny=ny,
        exx=flat_plate.affine_exx,
        eyy=0.0,
        exy=0.0,
        clamped=False,
    )

    # avg_stresses = flat_plate.run_static_analysis(write_soln=True)

    tacs_eigvals, _ = flat_plate.run_buckling_analysis(
        sigma=5.0, num_eig=6, write_soln=False
    )

    affine_AR = flat_plate.affine_aspect_ratio

    # expect to get ~2.9
    # since k0-2D* = (m a_0/b_0)^2 + (b_0/a_0/m)^2
    # in affine space and D*=1 and k0-2D* = 2.5 in Brunelle paper (buckling-analysis section)
    # "Generic Buckling Curves For Specially Orthotropic Rectangular Plates"
    # but only works in thin plate limit (very thin)

    # the non-dimensional buckling coefficient
    kx0_FEA_list += [tacs_eigvals[0]]

    # compute kx0 from the closed-form solution of the affine transform paper
    # loop over each mode and find the min non-dimensional buckling coefficient
    kx0_modes = []
    for m in range(1, 6):
        kx0_modes += [
            2.0 * flat_plate.Dstar + (m / affine_AR) ** 2 + (affine_AR / m) ** 2
        ]
    kx0_exact = min(kx0_modes)
    kx0_CF_list += [kx0_exact]

    data_dict = {
        "AR": AR_list[-1:],
        "kx0_FEA": kx0_FEA_list[-1:],
        "kx0_CF": kx0_CF_list[-1:],
    }

    # write the data to a file
    if comm.rank == 0:
        df = pd.DataFrame(data_dict)
        if ct == 0:
            df.to_csv(csv_file, mode="w", index=False)
        else:
            df.to_csv(csv_file, mode="a", index=False, header=False)
    comm.Barrier()
    ct += 1

# plot the Closed-Form versus the FEA for affine equations
if comm.rank == 0:
    plt.style.use(niceplots.get_style())
    plt.figure("affine")
    plt.plot(AR_list, kx0_CF_list, label="closed-form")
    plt.plot(AR_list, kx0_FEA_list, "o", label="FEA")
    plt.xlabel(r"$a_0/b_0$")
    plt.ylabel(r"$k_{x_0}$")
    plt.legend()
    plt.savefig("affine-AR-kx0.png", dpi=400)
