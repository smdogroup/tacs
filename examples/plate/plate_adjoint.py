# Import necessary libraries
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions, solver

# TACS Communicator
comm = MPI.COMM_WORLD

# Input file
bdfFileName = "plate.bdf"

# Load structural mesh from bdf file
mesh = TACS.MeshLoader(comm)
mesh.scanBDFFile(bdfFileName)

# Properties for the structure
rho           = 2500.0  # density, kg/m^3
E             = 70.0e9  # elastic modulus, Pa
nu            = 0.3     # poisson's ratio
kcorr         = 5.0/6.0 # shear correction factor
ys            = 350.0e6 # yield stress, Pa
min_thickness = 0.01    # minimum thickness of elements in m
max_thickness = 0.10    # maximum thickness of elemetns in m  
thickness     = 0.05    # currrent thickness of elements in m

# Properties for dynamics
gravity = elements.GibbsVector(np.array([0.0, 0.0, -9.81]))
v0      = elements.GibbsVector(np.array([0.25, 0.25, 0.25]))
w0      = elements.GibbsVector(np.array([0., 0., 0.]))

# Loop over components, creating stiffness and element object for each
num_components = mesh.getNumComponents()
for i in xrange(num_components):
    descriptor = mesh.getElementDescript(i)
    stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness, i,
                                 min_thickness, max_thickness)
    element = None
    if descriptor in ["CQUAD"]:
        element = elements.MITC(stiff, gravity, v0, w0)        
    mesh.setElement(i, element)

# Create tacs assembler object
tacs = mesh.createTACS(8)

# Create the function list
funcs = []
funcs.append(functions.ksfailure(tacs, 100.0))
funcs.append(functions.compliance(tacs))

# Setup design variables and allocate space for the function values
# and their derivatives wrt design variables
x = np.array([0.03])
fvals = np.zeros(2)
dfdx = np.zeros(2)

fvals_fd = np.zeros(2)
dfdx_fd = np.zeros(2)

# Parameters for time integration and adjoint solve
tinit             = 0.0
tfinal            = 0.02
num_steps_per_sec = 1000

# Test all BDF Integrators for these functions
for bdf_order in [1,2,3]:
    # BDF Integrator
    bdf = solver.BDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, bdf_order)
    bdf.setPrintLevel(0)
    bdf.setJacAssemblyFreq(1)
    bdf.setFunction(funcs)

    bdf.getFuncGrad(1, x, fvals, dfdx)
    bdf.getFDFuncGrad(1, x, fvals_fd, dfdx_fd, 1.0e-6)

    print("Error for BDF order ", bdf_order)
    print(fvals-fvals_fd, dfdx-dfdx_fd)

# Test all DIRK Integrators for these functions
for num_stages in [1,2,3]:
    # DIRK Integrator
    dirk = solver.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, num_stages)
    dirk.setPrintLevel(1)
    dirk.setJacAssemblyFreq(1)
    dirk.setFunction(funcs)

    dirk.getFuncGrad(1, x, fvals, dfdx)
    dirk.getFDFuncGrad(1, x, fvals_fd, dfdx_fd, 1.0e-6)

    print("Error for DIRK order ", num_stages + 1)
    print(fvals-fvals_fd, dfdx-dfdx_fd)
