# Import numpy
import numpy as np

# Import MPI
from mpi4py import MPI

# Import TACS and assorted repositories
from tacs import TACS, elements, constitutive, functions

# Set the MPI communicator
comm = MPI.COMM_WORLD

# Create the stiffness object
props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
stiff = constitutive.PlaneStressConstitutive(props)

# Set up the basis function
model = elements.LinearThermoelasticity2D(stiff)
basis = elements.QuadraticTriangleBasis()
elem = elements.Element2D(model, basis)

# Allocate the TACSCreator object
varsPerNode = model.getVarsPerNode()
creator = TACS.Creator(comm, varsPerNode)

if comm.rank == 0:
    # Create the elements
    nx = 25
    ny = 25

    # Set the nodes
    nnodes = (2 * nx + 1) * (2 * ny + 1)
    nelems = 2 * nx * ny
    nodes = np.arange(nnodes).reshape((2 * nx + 1, 2 * ny + 1))

    conn = []
    for j in range(ny):
        for i in range(nx):
            # Append the first set of nodes
            conn.append(
                [
                    nodes[2 * i, 2 * j],
                    nodes[2 * i + 2, 2 * j],
                    nodes[2 * i + 2, 2 * j + 2],
                    nodes[2 * i + 1, 2 * j],
                    nodes[2 * i + 2, 2 * j + 1],
                    nodes[2 * i + 1, 2 * j + 1],
                ]
            )

            # Append the second set of nodes
            conn.append(
                [
                    nodes[2 * i, 2 * j + 2],
                    nodes[2 * i, 2 * j],
                    nodes[2 * i + 2, 2 * j + 2],
                    nodes[2 * i, 2 * j + 1],
                    nodes[2 * i + 1, 2 * j + 1],
                    nodes[2 * i + 1, 2 * j + 2],
                ]
            )

    # Set the node pointers
    conn = np.array(conn, dtype=np.intc).flatten()
    ptr = np.arange(0, 6 * nelems + 1, 6, dtype=np.intc)
    elem_ids = np.zeros(nelems, dtype=np.intc)
    creator.setGlobalConnectivity(nnodes, ptr, conn, elem_ids)

    # Set up the boundary conditions
    bcnodes = np.array(nodes[0, :], dtype=np.intc)

    # Set the boundary condition variables
    nbcs = 2 * bcnodes.shape[0]
    bcvars = np.zeros(nbcs, dtype=np.intc)
    bcvars[:nbcs:2] = 0
    bcvars[1:nbcs:2] = 1

    # Set the boundary condition pointers
    bcptr = np.arange(0, nbcs + 1, 2, dtype=np.intc)
    creator.setBoundaryConditions(bcnodes, bcvars, bcptr)

    # Set the node locations
    Xpts = np.zeros(3 * nnodes)
    x = np.linspace(0, 10, 2 * nx + 1)
    y = np.linspace(0, 10, 2 * nx + 1)
    for j in range(2 * ny + 1):
        for i in range(2 * nx + 1):
            Xpts[3 * nodes[i, j]] = x[i]
            Xpts[3 * nodes[i, j] + 1] = y[j]

    # Set the node locations
    creator.setNodes(Xpts)

# Set the elements
elements = [elem]
creator.setElements(elements)

# Create the tacs assembler object
assembler = creator.createTACS()

res = assembler.createVec()
ans = assembler.createVec()
mat = assembler.createSchurMat()

# Create the preconditioner for the corresponding matrix
pc = TACS.Pc(mat)

# Assemble the Jacobian
alpha = 1.0
beta = 0.0
gamma = 0.0
assembler.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()

res.getArray()[:] = 1.0
assembler.applyBCs(res)
pc.applyFactor(res, ans)
ans.scale(-1.0)

assembler.setVariables(ans)

# Create the function list
funcs = []

# Create the KS function
ksweight = 100.0
for i in range(1):
    funcs.append(functions.KSFailure(assembler, ksWeight=ksweight))

func_vals = assembler.evalFunctions(funcs)

# Set the element flag
flag = (
    TACS.OUTPUT_CONNECTIVITY
    | TACS.OUTPUT_NODES
    | TACS.OUTPUT_DISPLACEMENTS
    | TACS.OUTPUT_STRAINS
)
f5 = TACS.ToFH5(assembler, TACS.PLANE_STRESS_ELEMENT, flag)
f5.writeToFile("triangle_test.f5")
