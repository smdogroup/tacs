"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from tacs import buckling_surrogate
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD

flat_plate = buckling_surrogate.FlatPlateAnalysis(
    comm=comm,
    bdf_file="plate.bdf",
    a=2.44010,
    b=0.329,
    h=0.00533,
    E11=82.14e9,
    nu12=0.1487,
    E22=16.656e9,  # set to None if isotropic
    G12=3.189e9,  # set to None if isotropic
)

flat_plate.generate_bdf(
    nx=30,
    ny=30,
    exx=0.0,
    eyy=0.0,
    exy=flat_plate.affine_exy,
    clamped=False,
)

# avg_stresses = flat_plate.run_static_analysis(write_soln=True)

tacs_eigvals, errors = flat_plate.run_buckling_analysis(
    sigma=10.0, num_eig=12, write_soln=True
)

# in previous monte carlo iteration, this trial got lam1 = 10.54
# whereas here we get lam1 = 2.56 with 1e-14 error in the eigenvalue solver
# seems like the error was bad on that solution? Need a way from python to check error in solution..

print(f"tacs eigvals = {tacs_eigvals}")
print(f"errors = {errors}")
