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
keys = ["Lx", "Ly", "thick", "AR", "SR", "D", "Nx", "Nxcrit"]
data_dict = {key:[] for key in keys}


# make a folder for the pandas data
cpath = os.path.dirname(__file__)
data_folder = os.path.join(cpath, "data")
if not os.path.exists(data_folder):
    os.mkdir(data_folder)

# fixed material parameters for rn
# using Aluminum alloy
E = 70e9; nu = 0.33

n_Lx = 20
n_AR = 20
n_thick = 10
N = n_Lx * n_AR * n_thick

ct = 1
for log_Lx in np.linspace(-2.0, 2.0, n_Lx):
    Lx = np.power(10, log_Lx)
    data_dict["Lx"] += [Lx]

    for log_AR in np.linspace(-1.5, 1.5, n_AR):
        aspect_ratio = np.power(10, log_AR)
        Ly = Lx / aspect_ratio
        data_dict["AR"] += [aspect_ratio]
        data_dict["Ly"] += [Ly]

        # compute nx, ny based on aspect ratio
        if aspect_ratio > 1.0:
            ny = 20
            nx = int(ny * aspect_ratio)
            if nx > 40:
                nx = 40
        else: # aspect ratio < 1.0
            nx = 20
            ny = int(nx / aspect_ratio)
            if ny > 40:
                ny = 40

        for log_thick in np.linspace(-2.0, 1.0, n_thick):
            print(f"\n\nData run {ct}/{N}")
            ct += 1

            thick = np.power(10, log_thick)
            data_dict["thick"] += [thick]
            slenderness_ratio = Lx / thick
            data_dict["SR"] += [slenderness_ratio]
            # bending stiffness
            D = E*thick**3 / 12.0 / (1 - nu**2)
            data_dict["D"] += [D]

            # make the plate and run each of the analyses
            generate_plate(
                Lx=Lx, Ly=Ly, nx=nx, ny=ny, exx=0.001, clamped=False
            )

            avg_stresses = run_static_analysis(thickness=thick, E=E, nu=nu, write_soln=False)
            Sx = avg_stresses[0]
            Nx = Sx * thick
            data_dict["Nx"] += [Nx]

            tacs_eigvals = run_buckling_analysis(
                thickness=thick, E=E, nu=nu, sigma=30.0, num_eig=6, write_soln=False
            )

            # compute Nxcrit (lowest eigenvalue)
            lam1 = tacs_eigvals[0]
            Nxcrit = Nx * lam1
            data_dict["Nxcrit"] += [Nxcrit]

# convert to a pandas dataframe and save it in the data folder
df = pd.DataFrame(data_dict)
csv_file = os.path.join(data_folder, "Nxcrit.csv")
df.to_csv(csv_file)