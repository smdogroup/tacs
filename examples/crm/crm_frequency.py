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
rho = 2500.0  # density, kg/m^3
E = 70e9  # elastic modulus, Pa
nu = 0.3  # poisson's ratio
kcorr = 5.0 / 6.0  # shear correction factor
ys = 350e6  # yield stress, Pa
min_thickness = 0.002
max_thickness = 0.20
thickness = 0.02

# Loop over components, creating stiffness and element object for each
num_components = struct_mesh.getNumComponents()
for i in range(num_components):
    descriptor = struct_mesh.getElementDescript(i)
    # Setup (isotropic) property and constitutive objects
    prop = constitutive.MaterialProperties(rho=rho, E=E, nu=nu, ys=ys)
    # Set one thickness dv for every component
    stiff = constitutive.IsoShellConstitutive(
        prop, t=thickness, tMin=min_thickness, tMax=max_thickness, tNum=i
    )

    element = None
    transform = None
    if descriptor in ["CQUAD", "CQUADR", "CQUAD4"]:
        element = elements.Quad4Shell(transform, stiff)
    struct_mesh.setElement(i, element)

# Create tacs assembler object from mesh loader
assembler = struct_mesh.createTACS(6)

# Solve the eigenvalue problem
M = assembler.createSchurMat()
K = assembler.createSchurMat()

pc = TACS.Pc(K)

# Assemble and factor the stiffness/Jacobian matrix. Factor the
# Jacobian and solve the linear system for the displacements
alpha = 1.0
beta = 0.0
gamma = 0.0
assembler.assembleJacobian(alpha, beta, gamma, None, K)
pc.factor()  # LU factorization of stiffness matrix

subspace = 100
restarts = 2
gmres = TACS.KSM(K, pc, subspace, restarts)

# Guess for the lowest natural frequency
sigma_hz = 1.0
sigma = 2.0 * np.pi * sigma_hz

# Create the frequency analysis object
num_eigs = 5
freq = TACS.FrequencyAnalysis(
    assembler, sigma, M, K, gmres, num_eigs=num_eigs, eig_tol=1e-12
)

# Solve the frequency analysis problem
freq.solve(print_level=2)

# Output for visualization
flag = (
    TACS.OUTPUT_CONNECTIVITY
    | TACS.OUTPUT_NODES
    | TACS.OUTPUT_DISPLACEMENTS
    | TACS.OUTPUT_STRAINS
    | TACS.OUTPUT_STRESSES
    | TACS.OUTPUT_EXTRAS
)
f5 = TACS.ToFH5(assembler, TACS.BEAM_OR_SHELL_ELEMENT, flag)

vec = assembler.createVec()
for i in range(num_eigs):
    freq.extractEigenvector(i, vec)
    assembler.setVariables(vec)
    f5.writeToFile("ucrm_mode%d.f5" % (i))
