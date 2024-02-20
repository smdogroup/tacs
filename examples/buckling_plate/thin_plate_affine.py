"""
Sean Engelstad
Feb 2024, GT SMDO Lab
Goal is to generate data for pure uniaxial compression failure in the x-direction
For now just simply supported only..
"""

from _generate_plate import generate_plate
from _static_analysis import run_static_analysis
from _buckling_analysis import run_buckling_analysis
import numpy as np
import pandas as pd
import os

"""
Let's consider the following parameters in exploring the data..
size: Lx, Ly, h (thickness)
material: E, nu, or D = E*h^3/12/(1-nu)
others: aspect ratio Lx/Ly
    slenderness ratio Lx/h, or Ly/h related to shear strain energy
    compute lambda * Nx = Nxcrit as the output (Nxcrit)
"""

# make a data dict for storing pandas data (later converted to dataframe)
keys = ["Lx", "Ly", "thick", "AR", "SR", "D", "Nx", "mu_x", "Nxcrit", "Nxcrit*"]
data_dict = {key: [] for key in keys}

# make a folder for the pandas data
cpath = os.path.dirname(__file__)
data_folder = os.path.join(cpath, "data")
if not os.path.exists(data_folder):
    os.mkdir(data_folder)

# fixed material parameters for rn
# using Aluminum alloy
E = 70e9
nu = 0.33

nx = 30
ny = 20

# try several different cases with different h, Lx, Ly, etc. to demonstrate thin-plate affine transform

# use monte carlo simulation to randomly sample

# for AR in np.linspace(0.5, 10.0, 100):
for AR in [1.0]:
    # in affine space
    exx_affine = 1e-3
    b0 = 1.0
    a0 = b0 * AR
    h0 = (12 * (1 - nu**2) / E) ** (1.0 / 3)

    generate_plate(Lx=a0, Ly=b0, nx=nx, ny=ny, exx=exx_affine, clamped=False)

    avg_stresses_affine = run_static_analysis(
        thickness=h0, E=E, nu=nu, write_soln=False
    )
    Sx = avg_stresses_affine[0]
    Nx0 = -Sx * h0  # positive is compressive

    tacs_eigvals_affine = run_buckling_analysis(
        thickness=h0, E=E, nu=nu, sigma=0.01, num_eig=6, write_soln=False
    )

    # not in affine space
    h = 0.01
    D = E * h**3 / 12 / (1 - nu**2)
    D4 = D**0.25
    a = D4 * a0
    b = D4 * b0
    SR = b / h
    exx = D**0.5 * exx_affine
    print(f"slenderness ratio = {SR}")

    generate_plate(Lx=a0, Ly=b0, nx=nx, ny=ny, exx=exx, clamped=False)

    avg_stresses = run_static_analysis(thickness=h, E=E, nu=nu, write_soln=False)
    Sx = avg_stresses[0]
    Nx = -Sx * h  # positive is compressive

    tacs_eigvals = run_buckling_analysis(
        thickness=h, E=E, nu=nu, sigma=0.01, num_eig=6, write_soln=False
    )

    print(f"affine tacs eigenvalues = {tacs_eigvals_affine}")
    print(f"regular tacs eigenvalues = {tacs_eigvals}")
