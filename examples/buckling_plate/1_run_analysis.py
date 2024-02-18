"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from _generate_plate import generate_plate
from _static_analysis import run_static_analysis
from _buckling_analysis import run_buckling_analysis
import numpy as np

generate_plate(Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.001, eyy=0.0, exy=0.0, clamped=False)
# run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)

tacs_eigvals = run_buckling_analysis(
    thickness=0.07, E=70e9, nu=0.33, sigma=30.0, num_eig=12, write_soln=True
)

print(f"tacs eigvals = {tacs_eigvals}")
