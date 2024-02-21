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
    a=1.0,
    b=1.0,
    h=0.005,
    E11=70e9,
    nu12=0.33,
    E22=None,  # set to None if isotropic
    G12=None,  # set to None if isotropic
)

flat_plate.generate_bdf(
    nx=50,
    ny=20,
    exx=flat_plate.affine_exx,
    eyy=0.0,
    exy=0.0,
    clamped=False,
)

avg_stresses = flat_plate.run_static_analysis(write_soln=True)

tacs_eigvals = flat_plate.run_buckling_analysis(sigma=10.0, num_eig=12, write_soln=True)

# expect to get ~4.5
# since k0-2D* = (m a_0/b_0)^2 + (b_0/a_0/m)^2
# in affine space and D*=1 and k0-2D* = 2.5 in Brunelle paper (buckling-analysis section)
# "Generic Buckling Curves For Specially Orthotropic Rectangular Plates"
# but only works in thin plate limit (very thin)

print(f"tacs eigvals = {tacs_eigvals}")