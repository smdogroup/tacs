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

# Create distributed vector for node locations and fill it
X = tacs.createNodeVec()
tacs.getNodes(X)

# Set up functions to evaluate
mass = functions.StructuralMass(tacs)
ksWeight = 50.0
ksfailure = functions.KSFailure(tacs, ksWeight)
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
tic = MPI.Wtime()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = MPI.Wtime()

if tacs_comm.rank == 0:
    print "Assembled and factored linear system..."
    print "Time to assemble and factor: ", toc - tic
    print

# Get numpy array from distributed vector object and write in loads
force_array = forces.getArray() 
force_array[2::6] += 100.0 # uniform load in z direction
tacs.applyBCs(forces)

# Solve the linear system
tic = MPI.Wtime()
pc.applyFactor(forces, ans)
toc = MPI.Wtime()
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
fXptSens = tacs.createNodeVec()
tacs.evalXptSens(funcs[0], fXptSens)

svSens = tacs.createVec()
tacs.evalSVSens(funcs[0], svSens)
svSens.scale(-1.0)

adjoint = tacs.createVec()
pc.applyFactor(svSens, adjoint)

product = tacs.createNodeVec()
tacs.evalAdjointResXptSensProduct(adjoint, product)
fXptSens.axpy(1.0, product)

# Perturb and set the node locations
pert = tacs.createNodeVec()
pert.setRand()

if TACS.dtype is np.complex:
    dh = 1.0e-30
    X.axpy(dh*1j, pert)
else:
    dh = 1.0e-6
    X.axpy(dh, pert)

tacs.setNodes(X)

# Compute the perturbed solution
tacs.zeroVariables()
tic = MPI.Wtime()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
toc = MPI.Wtime()

if tacs_comm.rank == 0:
    print "Assembled and factored perturbed linear system..."
    print "Time to assemble and factor: ", toc - tic
    print

pc.applyFactor(forces, ans)
tacs.setVariables(ans)

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
result = fXptSens.dot(pert).real
if tacs_comm.rank == 0:
    print 'FD:      ', fd[0]
    print 'Result:  ', result
    print 'Rel err: ', (result - fd[0])/fd[0]

# # Output for visualization 
# flag = (TACS.ToFH5.NODES |
#         TACS.ToFH5.DISPLACEMENTS |
#         TACS.ToFH5.STRAINS)
# f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)
# f5.writeToFile('ucrm.f5')
