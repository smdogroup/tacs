# A demonstration of basic functions of the Python interface for TACS: loading a
# mesh, creating elements, evaluating functions, solution, and output

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
for i in xrange(num_components):
    descriptor = struct_mesh.getElementDescript(i)
    stiff = constitutive.isoFSDT(rho, E, nu, kcorr, ys, thickness, i,
                                 min_thickness, max_thickness)
    element = None
    if descriptor in ["CQUAD", "CQUADR", "CQUAD4"]:
        element = elements.MITCShell(2, stiff, component_num=i)
    struct_mesh.setElement(i, element)

# Create tacs assembler object from mesh loader
tacs = struct_mesh.createTACS(6)

# Set up functions to evaluate (just computing mass of structure as example)
mass = functions.StructuralMass(tacs)
funclist = [mass]

# Evaluate functions
fvals = tacs.evalFunctions(funclist)
massval = fvals[0]
if tacs_comm.rank == 0:
    print "Mass: ", massval

x = np.zeros(num_components, TACS.dtype)
tacs.getDesignVars(x)

# Create matrix and distributed vectors
res = tacs.createVec()
forces = tacs.createVec()
ans = tacs.createVec()
mat = tacs.createFEMat()

# Create preconditioner
pc = TACS.Pc(mat)

# Assemble the Jacobian and factor
alpha = 1.0
beta = 0.0
gamma = 0.0
tacs.zeroVariables()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()

# Get numpy array from distributed vector object and write in loads
force_array = forces.getArray() 
force_array[2::6] += 100.0 # uniform load in z direction
forces.applyBCs()

# Solve the linear system
pc.applyFactor(forces, ans)
tacs.setVariables(ans)

# Evaluate the function
ksWeight = 100.0
funcs = [functions.KSFailure(tacs, ksWeight)]
fvals1 = tacs.evalFunctions(funcs)

# Evaluate the derivative
fdvSens = np.zeros(x.shape, TACS.dtype)
product = np.zeros(x.shape, TACS.dtype)
tacs.evalSVSens(funcs[0], res)
pc.applyFactor(res, ans)
tacs.evalDVSens(funcs[0], fdvSens)
tacs.evalAdjointResProduct(ans, product)
fdvSens = fdvSens - product

# Evaluate the result
result = np.sum(fdvSens).real

# Set the complex step
dh = 1e-6
if TACS.dtype is np.complex:
    dh = 1e-30
    x = x + 1j*dh
else:
    x = x + dh

# Set the design variables
tacs.setDesignVars(x)

# Compute the perturbed solution
tacs.zeroVariables()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
pc.applyFactor(forces, ans)
tacs.setVariables(ans)

# Set the new complex step
fvals2 = tacs.evalFunctions(funcs)

if TACS.dtype is np.complex:
    fd = fvals2.imag/dh
else:
    fd = (fvals2 - fvals1)/dh

if tacs_comm.rank == 0:
    print 'FD:     ', fd[0]
    print 'Result: ', result

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)
f5.writeToFile('ucrm.f5')
