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

n_Lx = 2
n_AR = 2
n_thick = 2
N = n_Lx * n_AR * n_thick

# TODO : investigate why there are so many negative eigenvalue
# solutions for the uniaxial compression case ? maybe the guess was too far off and error high?
# may need to add an interface to check the errors in the eigenvalue analysis
ct = 1
for log_Lx in np.linspace(-1.0, 2.0, n_Lx):
    Lx = np.power(10, log_Lx)

    for log_AR in np.linspace(-1.5, 1.5, n_AR):
        aspect_ratio = np.power(10, log_AR)
        Ly = Lx / aspect_ratio

        # compute nx, ny based on aspect ratio
        if aspect_ratio > 1.0:
            ny = 20
            nx = int(ny * aspect_ratio)
            if nx > 40:
                nx = 40
        else:  # aspect ratio < 1.0
            nx = 20
            ny = int(nx / aspect_ratio)
            if ny > 40:
                ny = 40

        for log_thick in np.linspace(-2.0, 1.0, n_thick):
            print(f"\n\nData run {ct}/{N}")
            ct += 1

            thick = np.power(10, log_thick)
            slenderness_ratio = Lx / thick
            # bending stiffness
            D = E * thick**3 / 12.0 / (1 - nu**2)

            # make the plate and run each of the analyses
            exx = 1e-3
            #exx = 1e-3 * (thick/0.07)**2 / Lx**2
            generate_plate(Lx=Lx, Ly=Ly, nx=nx, ny=ny, exx=exx, clamped=False)

            avg_stresses = run_static_analysis(
                thickness=thick, E=E, nu=nu, write_soln=False
            )
            Sx = avg_stresses[0]
            Nx = -Sx * thick  # positive is compressive

            tacs_eigvals = run_buckling_analysis(
                thickness=thick, E=E, nu=nu, sigma=30.0, num_eig=6, write_soln=False
            )

            # get only positive eigenvalues
            pos_eigvals = [eigval for eigval in tacs_eigvals if eigval > 0.0]

            # compute Nxcrit (lowest eigenvalue)
            lam1 = tacs_eigvals[0]
            mu_x = 1.0/lam1
            Nxcrit = Nx * lam1

            # non-dimensionalize the critical load
            # model as function of non-dimensional parameters AR, SR, etc.
            Nxcrit_star = Nxcrit / D * Lx**2

            # make sure the eigenvalue is in a reasonable range (otherwise analysis failed)
            # TODO : potentially change the eigenvalue solver => or guess better
            # and change this failure condition to use errors
            if 0.01 <= lam1 <= 500.0:  # only log data if analysis succeeded
                data_dict["Lx"] += [Lx]
                data_dict["AR"] += [aspect_ratio]
                data_dict["Ly"] += [Ly]
                data_dict["thick"] += [thick]
                data_dict["SR"] += [slenderness_ratio]
                data_dict["D"] += [D]
                data_dict["mu_x"] += [mu_x]
                data_dict["Nx"] += [Nx]
                data_dict["Nxcrit"] += [Nxcrit]
                data_dict["Nxcrit*"] += [Nxcrit_star]

# convert to a pandas dataframe and save it in the data folder
df = pd.DataFrame(data_dict)
csv_file = os.path.join(data_folder, "Nxcrit.csv")
df.to_csv(csv_file)
