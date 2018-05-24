# A demonstration of basic functions of the Python interface for TACS: loading a
# mesh, creating elements, evaluating functions, solution, and output
from __future__ import print_function

# Import necessary libraries
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

# Load structural mesh from BDF file
tacs_comm = MPI.COMM_WORLD
struct_mesh = TACS.MeshLoader(tacs_comm)
struct_mesh.scanBDFFile("CRM_box_2nd.bdf")

# Set constitutive properties
rho = 2500.0 # density, kg/m^3
E = 70e9 # elastic modulus, Pa
nu = 0.3 # poisson's ratio
kcorr = 5.0 / 6.0 # shear correction factor
ys = 350e6 # yield stress, Pa
min_thickness = 0.002
max_thickness = 0.20
thickness = 0.02

# Loop over components, creating stiffness and element object for each
num_components = struct_mesh.getNumComponents()
for i in range(num_components):
    descriptor = struct_mesh.getElementDescript(i)
    stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness, i,
                                 min_thickness, max_thickness)
    element = None
    if descriptor in ["CQUAD", "CQUADR", "CQUAD4"]:
        element = elements.MITCShell(2, stiff, component_num=i)
    struct_mesh.setElement(i, element)

# Create tacs assembler object from mesh loader
assembler = struct_mesh.createTACS(6)

# Solve the eigenvalue problem
M = assembler.createFEMat()
K = assembler.createFEMat()

pc = TACS.Pc(K)
subspace = 100
restarts = 2
gmres = TACS.KSM(K, pc, subspace, restarts)

# Guess for the lowest natural frequency
sigma_hz = 1.0
sigma = 2.0*np.pi*sigma_hz

# Create the frequency analysis object
num_eigs = 5
freq = TACS.FrequencyAnalysis(assembler, sigma, M, K, gmres, 
                              num_eigs=num_eigs, eig_tol=1e-12)

# Solve the frequency analysis problem
freq.solve()

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS |
        TACS.ToFH5.EXTRAS)
f5 = TACS.ToFH5(assembler, TACS.PY_SHELL, flag)

vec = assembler.createVec()
for i in range(num_eigs):
    freq.extractEigenvector(i, vec)
    assembler.setVariables(vec)
    f5.writeToFile('ucrm_mode%d.f5'%(i))
