# A demonstration of basic functions of the Python interface for TACS: loading a
# mesh, creating elements, evaluating functions, solution, and output

# Import necessary libraries
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions
import time

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

# Set up functions to evaluate
mass = functions.StructuralMass(tacs)
ksWeight = 50.0
ksfailure = functions.KSFailure(tacs, ksWeight)
funclist = [ksfailure]

# Evaluate functions
fvals = tacs.evalFunctions(funclist)
fval = fvals[0]
if tacs_comm.rank == 0:
    print "function eval: ", fval

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
print "Assembling and factoring linear system..."
tic = time.clock()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = time.clock()
print "Time to assemble and factor: ", toc - tic

# Get numpy array from distributed vector object and write in loads
force_array = forces.getArray() 
force_array[2::6] += 100.0 # uniform load in z direction
forces.applyBCs()

# Solve the linear system
print "Solving linear system..."
tic = time.clock()
pc.applyFactor(forces, ans)
toc = time.clock()
print "Time to solve: ", toc - tic
tacs.setVariables(ans)

# Evaluate the function
ksWeight = 100.0
funcs = [functions.KSFailure(tacs, ksWeight)]
fvals1 = tacs.evalFunctions(funcs)

# Create distributed vector for node locations and fill it
struct_X_vec = tacs.createNodeVec()
tacs.getNodes(struct_X_vec)
struct_X = struct_X_vec.getArray().astype(TACS.dtype)
struct_nnodes = len(struct_X)/3

# Evaluate the derivative w.r.t. node locations
dvsens_vec = tacs.createNodeVec()
tacs.evalXptSens(funcs[0], dvsens_vec)
dvsens = dvsens_vec.getArray()

svsens_vec = tacs.createVec()
tacs.evalSVSens(funcs[0], svsens_vec)
svsens_vec.scale(-1.0)

print "Solving for adjoint variables..."
adj_vars_vec = tacs.createVec()
tic = time.clock()
pc.applyFactor(svsens_vec, adj_vars_vec)
toc = time.clock()
print "Time to solve: ", toc - tic

product_vec = tacs.createNodeVec()
tacs.evalAdjointResXptSensProduct(adj_vars_vec, product_vec)
product = product_vec.getArray()
deriv = dvsens + product

# Perturb and set the node locations
struct_X_pert = struct_X
i_pert = 0
dh = 1.0e-6

if TACS.dtype is np.complex:
    dh = 1.0e-30
    struct_X_pert[i_pert] += 1j*dh
else:
    struct_X_pert[i_pert] += dh
    
struct_X[:] = struct_X_pert[:]
struct_X_vec.applyBCs()
tacs.setNodes(struct_X_vec)

# Compute the perturbed solution
tacs.zeroVariables()
print "Assembling and factoring perturbed linear system..."
tic = time.clock()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = time.clock()
print "Time to assemble and factor: ", toc - tic

print "Solving perturbed linear system..."
tic = time.clock()
pc.applyFactor(forces, ans)
toc = time.clock()
print "Time to solve: ", toc - tic
tacs.setVariables(ans)

# Evaluate the function for perturbed solution
fvals2 = tacs.evalFunctions(funcs)

if TACS.dtype is np.complex:
    fd = fvals2.imag/dh
else:
    fd = (fvals2 - fvals1)/dh

if tacs_comm.rank == 0:
    print 'FD:     ', fd[0]
    print 'Result: ', deriv[i_pert].real

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)
f5.writeToFile('ucrm.f5')
