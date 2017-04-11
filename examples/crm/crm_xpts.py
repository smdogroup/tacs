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

# Create distributed vector for node locations and fill it
struct_X_vec = tacs.createNodeVec()
tacs.getNodes(struct_X_vec)
struct_X = struct_X_vec.getArray()
struct_nnodes = len(struct_X)/3

# Set up functions to evaluate
mass = functions.StructuralMass(tacs)
ksWeight = 50.0
ksfailure = functions.KSFailure(tacs, ksWeight)
# funcs = [mass]
funcs = [ksfailure]

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
tic = time.clock()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = time.clock()

if tacs_comm.rank == 0:
    print "Assembled and factored linear system..."
    print "Time to assemble and factor: ", toc - tic
    print

# Get numpy array from distributed vector object and write in loads
force_array = forces.getArray() 
force_array[2::6] += 100.0 # uniform load in z direction
tacs.applyBCs(forces)

# Solve the linear system
tic = time.clock()
pc.applyFactor(forces, ans)
toc = time.clock()
tacs.setVariables(ans)

if tacs_comm.rank == 0:
    print "Solved linear system..."
    print "Time to solve: ", toc - tic
    print

# Evaluate the function
fvals1 = tacs.evalFunctions(funcs)
fval = fvals1[0]
if tacs_comm.rank == 0:
    print "unperturbed function eval: ", fval
    print

# Evaluate the derivative w.r.t. node locations
dvsens_vec = tacs.createNodeVec()
tacs.evalXptSens(funcs[0], dvsens_vec)
dvsens = dvsens_vec.getArray()

svsens_vec = tacs.createVec()
tacs.evalSVSens(funcs[0], svsens_vec)
svsens_vec.scale(-1.0)

adj_vars_vec = tacs.createVec()
tic = time.clock()
pc.applyFactor(svsens_vec, adj_vars_vec)
toc = time.clock()

if tacs_comm.rank == 0:
    print "Solved for adjoint variables..."
    print "Time to solve: ", toc - tic
    print

product_vec = tacs.createNodeVec()
tacs.evalAdjointResXptSensProduct(adj_vars_vec, product_vec)
product = product_vec.getArray()
print dvsens
print product
deriv = dvsens + product

# Perturb and set the node locations
pert_vec = tacs.createNodeVec()
pert_vec.setRand()

if TACS.dtype is np.complex:
    dh = 1.0e-30
    struct_X_vec.axpy(dh*1j, pert_vec)
else:
    dh = 1.0e-9
    struct_X_vec.axpy(dh, pert_vec)

tacs.testElement(0, 2)

pert = pert_vec.getArray().astype(TACS.dtype)
tacs.setNodes(struct_X_vec)

# Compute the perturbed solution
tacs.zeroVariables()
tic = time.clock()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = time.clock()

if tacs_comm.rank == 0:
    print "Assembled and factored perturbed linear system..."
    print "Time to assemble and factor: ", toc - tic
    print

tic = time.clock()
pc.applyFactor(forces, ans)
toc = time.clock()
tacs.setVariables(ans)

if tacs_comm.rank == 0:
    print "Solved perturbed linear system..."
    print "Time to solve: ", toc - tic
    print

# Evaluate the function for perturbed solution
fvals2 = tacs.evalFunctions(funcs)
fval2 = fvals2[0]
if tacs_comm.rank == 0:
    print "perturbed function eval: ", fval2
    print

if TACS.dtype is np.complex:
    fd = fvals2.imag/dh
else:
    fd = (fvals2 - fvals1)/dh

# Print comparison of derivative evaluations
if tacs_comm.rank == 0:
    print 'FD:     ', fd[0]
    print 'Result: ', np.dot(deriv, pert).real

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)
f5.writeToFile('ucrm.f5')
