"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from .generate_plate import generate_plate
from .static_analysis import run_static_analysis
from .buckling_analysis import run_buckling_analysis

# 3 main verification cases for the buckling analysis
# 1, 2, 3
case = 1
run_static = False

# case 1 - simply supported, uniaxial compression buckling
if case == 1:
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.001, eyy=0.0, exy=0.0, clamped=False
    )
    run_buckling_analysis(
        thickness=0.07, E=70e9, nu=0.33, sigma=30.0, num_eig=12, write_soln=True
    )

# case 2 - 
if case == 2:
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.001, eyy=0.0, exy=0.0, clamped=False
    )
    # run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)
    run_buckling_analysis(
        thickness=0.07, E=70e9, nu=0.33, sigma=30.0, num_eig=12, write_soln=True
    )

# run the static analysis for debugging the resulting displacement field
if run_static:
    run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)