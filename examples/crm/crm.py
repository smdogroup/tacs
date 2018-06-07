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
tacs = struct_mesh.createTACS(6)

# Create the KS Function
ksWeight = 100.0
funcs = [functions.KSFailure(tacs, ksWeight)]
# funcs = [functions.StructuralMass(tacs)]
# funcs = [functions.Compliance(tacs)]

# Get the design variable values
x = np.zeros(num_components, TACS.dtype)
tacs.getDesignVars(x)

# Get the node locations
X = tacs.createNodeVec()
tacs.getNodes(X)
tacs.setNodes(X)

# Create the forces
forces = tacs.createVec()
force_array = forces.getArray() 
force_array[2::6] += 100.0 # uniform load in z direction
tacs.applyBCs(forces)

# Set up and solve the analysis problem
res = tacs.createVec()
ans = tacs.createVec()
u = tacs.createVec()
mat = tacs.createFEMat()
pc = TACS.Pc(mat)
subspace = 100
restarts = 2
gmres = TACS.KSM(mat, pc, subspace, restarts)

# Assemble the Jacobian and factor
alpha = 1.0
beta = 0.0
gamma = 0.0
tacs.zeroVariables()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()

# Solve the linear system
gmres.solve(forces, ans)
tacs.setVariables(ans)

# Evaluate the function
fvals1 = tacs.evalFunctions(funcs)

# Solve for the adjoint variables
adjoint = tacs.createVec()
tacs.evalSVSens(funcs[0], res)
gmres.solve(res, adjoint)

# Compute the total derivative w.r.t. material design variables
fdvSens = np.zeros(x.shape, TACS.dtype)
product = np.zeros(x.shape, TACS.dtype)
tacs.evalDVSens(funcs[0], fdvSens)
tacs.evalAdjointResProduct(adjoint, product)
fdvSens = fdvSens - product

# Create a random direction along which to perturb the nodes
pert = tacs.createNodeVec()
X_array = X.getArray()
pert_array = pert.getArray()
pert_array[0::3] = X_array[1::3]
pert_array[1::3] = X_array[0::3]
pert_array[2::3] = X_array[2::3]

# Compute the total derivative w.r.t. nodal locations
fXptSens = tacs.createNodeVec()
product = tacs.createNodeVec()
tacs.evalXptSens(funcs[0], fXptSens)
tacs.evalAdjointResXptSensProduct(adjoint, product)
fXptSens.axpy(-1.0, product)

# Set the complex step
xpert = np.random.uniform(size=x.shape)
xpert = tacs_comm.bcast(xpert, root=0)
if TACS.dtype is np.complex:
    dh = 1e-30
    xnew = x + 1j*dh*xpert
else:
    dh = 1e-6
    xnew = x + dh*xpert

# Set the design variables
tacs.setDesignVars(xnew)

# Compute the perturbed solution
tacs.zeroVariables()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
gmres.solve(forces, u)
tacs.setVariables(u)

# Evaluate the function for perturbed solution
fvals2 = tacs.evalFunctions(funcs)

if TACS.dtype is np.complex:
    fd = fvals2.imag/dh
else:
    fd = (fvals2 - fvals1)/dh

result = np.dot(xpert, fdvSens)
if tacs_comm.rank == 0:
    print('FD:      ', fd[0])
    print('Result:  ', result)
    print('Rel err: ', (result - fd[0])/result)

# Reset the old variable values
tacs.setDesignVars(x)

if TACS.dtype is np.complex:
    dh = 1e-30
    X.axpy(dh*1j, pert)
else:
    dh = 1e-6
    X.axpy(dh, pert)

# Set the perturbed node locations
tacs.setNodes(X)

# Compute the perturbed solution
tacs.zeroVariables()
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()
gmres.solve(forces, u)
tacs.setVariables(u)

# Evaluate the function again
fvals2 = tacs.evalFunctions(funcs)

if TACS.dtype is np.complex:
    fd = fvals2.imag/dh
else:
    fd = (fvals2 - fvals1)/dh

# Compute the projected derivative
result = pert.dot(fXptSens)

if tacs_comm.rank == 0:
    print('FD:      ', fd[0])
    print('Result:  ', result)
    print('Rel err: ', (result - fd[0])/result)

# Output for visualization 
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS |
        TACS.ToFH5.EXTRAS)
f5 = TACS.ToFH5(tacs, TACS.PY_SHELL, flag)
f5.writeToFile('ucrm.f5')
