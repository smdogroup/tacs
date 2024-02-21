"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from tacs import buckling_surrogate
import numpy as np

flat_plate = buckling_surrogate.FlatPlateAnalysis(
    bdf_file="plate.bdf",
    a=1.0,
    b=0.2,
    h=0.01,
    E11=70e9,
    nu12=0.33,
    E22=70e9,  # set to None if isotropic
    G12=20e9,  # set to None if isotropic
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

print(f"tacs eigvals = {tacs_eigvals}")
