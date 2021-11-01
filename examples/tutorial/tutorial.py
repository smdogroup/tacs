'''
This example shows how to perform static structural analysis
for a 2D plate.
'''

from tacs import TACS, elements, constitutive, functions
from mpi4py import MPI
import numpy as np

'''
Create TACS Assembler
'''
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size
nx = 20  # number of elements in x direction
ny = 20  # number of elements in y direction
Lx = 1.0
Ly = 1.0
varsPerNode = 2
# varsPerNode = 3  # if we use LinearThermoelasticity2D element type
nodesPerProc = int((nx+1)*(ny+1)/size)
elemsPerProc = int(nx*ny/size)
numOwnedNodes = int(nodesPerProc)
numElements = int(elemsPerProc)
numDependentNodes = 0

# Adjust for the last processor
if (rank == size-1):
    numOwnedNodes = (nx+1)*(ny+1) - nodesPerProc*(size-1)
    numElements = nx*ny - elemsPerProc*(size-1)

assembler = TACS.Assembler.create(comm, varsPerNode,
                                  numOwnedNodes, numElements,
                                  numDependentNodes)

'''
Setup geometry, mesh, material and element
'''

# Set up partition
firstElem = rank*elemsPerProc
firstNode = rank*nodesPerProc
lastElem = (rank+1)*elemsPerProc
lastNode = (rank+1)*nodesPerProc
if (rank == size-1):
    lastElem = nx*ny
    lastNode = (nx+1)*(ny+1)

# Populate connectivity
ptr = np.zeros(numElements + 1, dtype=np.int32)
conn = np.zeros(4*numElements, dtype=np.int32)
ptr[0] = 0
k = 0
for elem in range(firstElem, lastElem):
    i = elem % nx
    j = elem // nx
    conn[4*k] =   i   + j*(nx+1)
    conn[4*k+1] = i+1 + j*(nx+1)
    conn[4*k+2] = i   + (j+1)*(nx+1)
    conn[4*k+3] = i+1 + (j+1)*(nx+1)
    ptr[k+1] = 4*(k+1)
    k += 1

# Set the connectivity
assembler.setElementConnectivity(ptr, conn)

# Create the isotropic material class
props = constitutive.MaterialProperties(rho=2700.0, E=70e3, nu=0.3, ys=270.0)

# Create basis, constitutive, element, etc
linear_basis = elements.LinearQuadBasis()
stiff = constitutive.PlaneStressConstitutive(props)
elements_list = []
for elem in range(firstElem, lastElem):
    stiff = constitutive.PlaneStressConstitutive(props, 1.0, elem)
    model = elements.LinearElasticity2D(stiff);
    # model = elements.LinearThermoelasticity2D(stiff);
    elements_list.append(elements.Element2D(model, linear_basis))

# Set elements into the mesh
assembler.setElements(elements_list)

# Set boundary conditions
for i in range(0, nx + 1):
    # Here nodal indexing is global
    nodes = np.array([i, i + (nx+1)*ny, i*(nx+1), (i+1)*(nx+1) - 1], dtype=np.int32)
    dof = np.array([0, 1, 2], dtype=np.int32)
    values = np.array([0.0, 0.0, 0.0])
    assembler.addBCs(nodes, dof, values)

# Done adding elements
assembler.initialize()

# Create the node location vector
X = assembler.createNodeVec()
Xpts = X.getArray()

# Get nodal locations
k = 0
for node in range(firstNode, lastNode):
    i = node % (nx + 1)
    j = node // (nx + 1)
    Xpts[k] = i*Lx/nx
    Xpts[k+1] = j*Ly/ny
    k += 3

assembler.reorderVec(X)  # Might not needed since we don't reorder the matrix
assembler.setNodes(X)

'''
Solve the static analysis
'''

# Create the finite element matrix
kmat = assembler.createSchurMat()

# Create the Schur preconditioner
pc = TACS.Pc(kmat)

# Allocate space for the vectors
force = assembler.createVec()
res   = assembler.createVec()
ans   = assembler.createVec()
tmp   = assembler.createVec()

# Set force
force_vals = force.getArray()
force_vals[::3] = 1.0
assembler.setBCs(force)

# Assemble the Jacobian for the governing equation
alpha = 1.0
beta  = 0.0
gamma = 0.0
assembler.assembleJacobian(alpha, beta, gamma, res, kmat)

# Factor the preconditioner
pc.factor()

# Create the solver and solve the problem
gmres_iters = 100;
nrestart = 2
is_flexible = 1
gmres = TACS.KSM(kmat, pc, gmres_iters, nrestart, is_flexible)
gmres.setMonitor(comm, freq=1)
gmres.solve(force, ans)  # ans = K\f

# Check residual: res = K*ans - f
kmat.mult(ans, tmp)
tmp.axpy(-1.0, force)
norm = tmp.norm()
if (rank == 0):
    print("|ku - f|: {:15.5e}".format(norm))  

assembler.setVariables(ans)

# Output f5 for visualization
write_flag = TACS.OUTPUT_CONNECTIVITY | \
             TACS.OUTPUT_NODES |  \
             TACS.OUTPUT_DISPLACEMENTS | \
             TACS.OUTPUT_STRAINS | \
             TACS.OUTPUT_STRESSES | \
             TACS.OUTPUT_EXTRAS
etype = TACS.PLANE_STRESS_ELEMENT
f5 = TACS.ToFH5(assembler, etype, write_flag)
f5.writeToFile("pytutorial.f5")

