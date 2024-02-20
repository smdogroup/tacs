"""
Sean Engelstad, Feb 2024
GT SMDO Lab
"""
from _generate_plate import generate_plate
from _static_analysis import run_static_analysis
from _buckling_analysis import run_buckling_analysis
import numpy as np

# model inputs
E = 70e9
nu = 0.33
a = 1.0
b = 0.2
nx = 30
ny = 20
exx = 0.001
h = 0.01

# if isotropic it's just the following
E11 = E
D = E * h**3 / 12.0 / (1 - nu**2)
D11 = D
D22 = D
E22 = E
G12 = E/2.0/(1+nu)

# affine transformation to compute k_{x0} buckling coefficients
# with this transformation => output lambda = k_{x0}
exx_T = np.pi**2 * np.sqrt(D11 * D22) / b**2 / h / E11

generate_plate(Lx=a, Ly=b, nx=nx, ny=ny, exx=exx_T, eyy=0.0, exy=0.0, clamped=False)
# run_static_analysis(thickness=h, E=E, nu=nu, write_soln=True)

tacs_eigvals = run_buckling_analysis(
    thickness=h, E11=E11, nu12=nu, E22=E22, G12=G12, sigma=10.0, num_eig=12, write_soln=True
)

print(f"tacs eigvals = {tacs_eigvals}")
