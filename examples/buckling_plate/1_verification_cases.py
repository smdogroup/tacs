"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from _generate_plate import generate_plate
from _static_analysis import run_static_analysis
from _buckling_analysis import run_buckling_analysis
import numpy as np

# 3 main verification cases for the buckling analysis
# 1, 2, 3
case = 2
run_static = False

# case 1 - simply supported, uniaxial compression buckling
if case == 1:
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.001, eyy=0.0, exy=0.0, clamped=False
    )
    tacs_eigvals = run_buckling_analysis(
        thickness=0.07, E=70e9, nu=0.33, sigma=30.0, num_eig=12, write_soln=True
    )

    # eigenvalues from Abaqus for comparison
    abaqus_eigvals = np.array([36.083, 38.000, 51.634, 72.896, 96.711, 113.94])
    rel_error = (tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

    print(
        f"\n\nVerification of case 1 - uniaxial compression,\n\tsimply supported plate buckling modes\n"
    )
    for i in range(6):
        print(f"mode {i+1} eigenvalues:")
        print(
            f"\ttacs = {tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
        )

# case 2 - pure shear, clamped plate
if case == 2:
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.0, eyy=0.0, exy=0.001, clamped=True
    )
    # run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)
    tacs_eigvals = run_buckling_analysis(
        thickness=0.07, E=70e9, nu=0.33, sigma=30.0, num_eig=12, write_soln=True
    )

    # every other shear mode has negative eigenvalue (since reverse shear load will still cause failure by sym)
    pos_tacs_eigvals = tacs_eigvals[::2]

    # eigenvalues from Abaqus for comparison
    abaqus_eigvals = np.array([111.79, 115.45, 169.71, 181.02, 236.06, 242.07])
    rel_error = (pos_tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

    print(f"\n\nVerification of case 2 - pure shear,\n\tclamped plate buckling modes\n")
    for i in range(6):
        print(f"mode {i+1} eigenvalues:")
        print(
            f"\ttacs = {pos_tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
        )

# case 3 - mixed shear and uniaxial compression case, clamped plate
if case == 3:
    # TODO : this case doesn't quite work yet.
    # note that positive exx indicates compression (our convention here)
    generate_plate(
        Lx=1.0, Ly=0.7, nx=30, ny=20, exx=0.001, eyy=0.0, exy=0.001, clamped=True
    )
    # run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)
    tacs_eigvals = run_buckling_analysis(
        thickness=0.07, E=70e9, nu=0.33, sigma=40.0, num_eig=12, write_soln=True
    )

    # every other shear mode has negative eigenvalue (since reverse shear load will still cause failure by sym)
    pos_tacs_eigvals = tacs_eigvals[::2]

    # eigenvalues from Abaqus for comparison
    #abaqus_eigvals = np.array([111.79, 115.45, 169.71, 181.02, 236.06, 242.07])
    rel_error = (pos_tacs_eigvals[:6] - abaqus_eigvals) / abaqus_eigvals

    print(f"\n\nVerification of case 3 - mixed compression + shear,\n\tclamped plate buckling modes\n")
    for i in range(6):
        print(f"mode {i+1} eigenvalues:")
        print(
            f"\ttacs = {pos_tacs_eigvals[i]:.4f}, abaqus = {abaqus_eigvals[i]:.4f}, rel err = {rel_error[i]:.4f}"
        )

# run the static analysis for debugging the resulting displacement field
if run_static:
    run_static_analysis(thickness=0.07, E=70e9, nu=0.33, write_soln=True)
