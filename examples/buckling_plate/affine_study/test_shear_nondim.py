import numpy as np
from tacs import buckling_surrogate
import matplotlib.pyplot as plt
import niceplots, pandas, os
from mpi4py import MPI

comm = MPI.COMM_WORLD

# goal here is to vary the plate size and show that the shear non-dimensionalization is correct
# "Generic Buckling Curves For Specially Orthotropic Rectangular Plates" by Brunelle

AR_list = []
kx0_FEA_list = []  # FEA
kx0_CF_list = []  # closed-form

# clear the csv file
csv_file = "shear-ND.csv"
if os.path.exists(csv_file) and comm.rank == 0:
    os.remove(csv_file)

AR_list = []
SR_list = []
b_list = []
kmin_list = []

ct = 0
for AR in [0.7, 1.0, 1.3]:

    # since material here is isotropic, D11=D22
    # and a0/b0 = a/b = AR
    SR = 100.0
    b = 100.0  # 10.0
    for b in np.linspace(10.0, 100.0, 10):

        flat_plate = buckling_surrogate.FlatPlateAnalysis(
            comm=comm,
            bdf_file="plate.bdf",
            a=AR * b,
            b=b,
            h=b / SR,  # very slender => so near thin plate limit
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
            exx=0.0,
            eyy=0.0,
            exy=flat_plate.affine_exy,
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
        kmin = np.abs(np.real(tacs_eigvals[0]))
        

        data_dict = {
            "AR": [AR],
            "SR" : [SR],
            "kmin" : [kmin],
            "b" : b
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
