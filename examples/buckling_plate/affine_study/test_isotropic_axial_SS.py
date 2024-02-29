import numpy as np
from tacs import buckling_surrogate
import matplotlib.pyplot as plt
import niceplots, pandas, os
from mpi4py import MPI

comm = MPI.COMM_WORLD

# vary aspect ratio and for D* = 1
# show that the equations match the affine transformation from
# "Generic Buckling Curves For Specially Orthotropic Rectangular Plates" by Brunelle

AR_list = []
kx0_FEA_list = []  # FEA
kx0_CF_list = []  # closed-form

# clear the csv file
csv_file = "axial-SS1.csv"
if os.path.exists(csv_file) and comm.rank == 0:
    os.remove(csv_file)

ct = 0
for AR in np.linspace(0.5, 5, 40):
    AR_list += [AR]

    # since material here is isotropic, D11=D22
    # and a0/b0 = a/b = AR
    affine_AR = AR
    slenderness = 100.0
    b = 100.0  # 100.0 originally

    flat_plate = buckling_surrogate.FlatPlateAnalysis(
        comm=comm,
        bdf_file="plate.bdf",
        a=AR * b,
        b=b,
        h=b / slenderness,  # very slender => so near thin plate limit
        E11=70e9,
        nu12=0.33,
        E22=None,  # set to None if isotropic
        G12=None,  # set to None if isotropic
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

    # expect to get ~4.5
    # since k0-2D* = (m a_0/b_0)^2 + (b_0/a_0/m)^2
    # in affine space and D*=1 and k0-2D* = 2.5 in Brunelle paper (buckling-analysis section)
    # "Generic Buckling Curves For Specially Orthotropic Rectangular Plates"
    # but only works in thin plate limit (very thin)

    # the non-dimensional buckling coefficient
    kx0_FEA_list += [np.real(tacs_eigvals[0])]

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
        df = pandas.DataFrame(data_dict)
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
