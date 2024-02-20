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
keys = ["Lx", "Ly", "thick", "AR", "SR", "D", "Nx", "Nxcrit", "Nxcrit*"]
data_dict = {key: [] for key in keys}

# make a folder for the pandas data
cpath = os.path.dirname(__file__)
data_folder = os.path.join(cpath, "data")
if not os.path.exists(data_folder):
    os.mkdir(data_folder)

# TODO : change this to a Monte carlo simulation over uniform scales or uniform log scales
# TODO : use the affine transformation of 1_run_analysis.py to make lambda values more consistent.

for i in range(1000):
    pass

# convert to a pandas dataframe and save it in the data folder
df = pd.DataFrame(data_dict)
csv_file = os.path.join(data_folder, "Nxcrit.csv")
df.to_csv(csv_file)
