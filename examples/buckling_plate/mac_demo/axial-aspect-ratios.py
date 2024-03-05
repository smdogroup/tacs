import numpy as np
from tacs import buckling_surrogate
import matplotlib.pyplot as plt
import niceplots, pandas, os
from mpi4py import MPI

comm = MPI.COMM_WORLD

ct = 0
flat_plates = []
# first solve uniaxial SS then shear SS and do modal comparison


# axial plate 1
# --------------------------------------------------
axial_plate1 = buckling_surrogate.FlatPlateAnalysis(
    comm=comm,
    bdf_file="plate.bdf",
    a=1.0,
    b=1.0,
    h=0.01,  # slender near thin plate limit
    E11=70e9,
    nu12=0.33,
    plate_name="axial1",
)

axial_plate1.generate_bdf(
    nx=30,
    ny=30,
    exx=axial_plate1.affine_exx,
    eyy=0.0,
    exy=0.0,
    clamped=False,
)

# avg_stresses = axial_plate.run_static_analysis(write_soln=True)

tacs_eigvals, _ = axial_plate1.run_buckling_analysis(
    sigma=5.0, num_eig=8, write_soln=True
)

# sheared plate
# --------------------------------------------------
axial_plate2 = buckling_surrogate.FlatPlateAnalysis(
    comm=comm,
    bdf_file="plate.bdf",
    a=2.0,
    b=1.0,
    h=0.01,  # slender near thin plate limit
    E11=70e9,
    nu12=0.33,
    plate_name="axial2",
)

axial_plate2.generate_bdf(
    nx=50,
    ny=25,
    exx=axial_plate2.affine_exx,
    eyy=0.0,
    exy=0.0,
    clamped=False,
)

# avg_stresses = shear_plate.run_static_analysis(write_soln=True)

tacs_eigvals, _ = axial_plate2.run_buckling_analysis(
    sigma=5.0, num_eig=16, write_soln=True
)

# perform modal assurance criterion between these two buckling analyses
# ---------------------------------------------------------------------
buckling_surrogate.FlatPlateAnalysis.mac_permutation(
    axial_plate1, axial_plate2, num_modes=6
)
