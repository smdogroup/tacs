import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, constitutive, functions

dtype = TACS.dtype

Lx = 10.0
Ly = 1.0

# running loads
Nx = 1.0

nx = 10
ny = 10

# Set the MPI communicator
comm = MPI.COMM_WORLD

# Create the stiffness object
props = constitutive.MaterialProperties(rho=2570.0, E=70e9, nu=0.3, ys=350e6)
stiff = constitutive.PlaneStressConstitutive(props)

# Set up the basis function
model = elements.LinearElasticity2D(stiff)
basis = elements.LinearQuadBasis()
elem = elements.Element2D(model, basis)

# Allocate the TACSCreator object
varsPerNode = model.getVarsPerNode()
creator = TACS.Creator(comm, varsPerNode)

if comm.rank == 0:
    num_elems = nx * ny
    num_nodes = (nx + 1) * (ny + 1)

    x = np.linspace(0, Lx, nx + 1, dtype)
    y = np.linspace(0, Ly, ny + 1, dtype)
    xyz = np.zeros([ny + 1, nx + 1, 3], dtype)
    xyz[:, :, 0], xyz[:, :, 1] = np.meshgrid(x, y)

    node_ids = np.arange(num_nodes).reshape(ny + 1, nx + 1)

    conn = []
    for j in range(nx):
        for i in range(ny):
            conn.append([node_ids[i, j],
                         node_ids[i + 1, j],
                         node_ids[i, j + 1],
                         node_ids[i + 1, j + 1]])

    conn = np.array(conn, dtype=np.intc).flatten()
    ptr = np.arange(0, 4 * num_elems + 1, 4, dtype=np.intc)
    comp_ids = np.zeros(num_elems, dtype=np.intc)

    creator.setGlobalConnectivity(num_nodes, ptr, conn, comp_ids)

    # Set up the boundary conditions
    bcnodes1 = np.array(node_ids[:, 0], dtype=np.intc)
    bcvars1 = np.zeros(ny + 1, dtype=np.intc)

    bcnodes2 = np.array(node_ids[0, :], dtype=np.intc)
    bcvars2 = np.ones(nx + 1, dtype=np.intc)

    bcnodes = np.append(bcnodes1, bcnodes2)
    bcvars = np.append(bcvars1, bcvars2)

    # Set the boundary condition pointers
    bcptr = np.arange(len(bcnodes) + 1, dtype=np.intc)
    creator.setBoundaryConditions(bcnodes, bcptr, bcvars)

    # Set the node locations
    creator.setNodes(xyz.flatten())

# Set the elements
elements = [elem]
creator.setElements(elements)

# Create the tacs assembler object
assembler = creator.createTACS()

res = assembler.createVec()
ans = assembler.createVec()
mat = assembler.createSchurMat()

local_xpts = assembler.createNodeVec()
assembler.getNodes(local_xpts)
local_xyz = local_xpts.getArray().copy()
local_num_nodes = len(local_xyz) // 3
local_xyz = local_xyz.reshape(local_num_nodes, 3)

# Create the preconditioner for the corresponding matrix
pc = TACS.Pc(mat)

# Assemble the Jacobian
alpha = 1.0
beta = 0.0
gamma = 0.0
assembler.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()

res_array = res.getArray().reshape(num_nodes, varsPerNode)

# Apply distributed force in the x-direction on right side of plate
res_array[local_xyz[:,0] == Lx, 0] = (Nx * Ly) / (ny + 1)

assembler.applyBCs(res)
pc.applyFactor(res, ans)
ans.scale(-1.0)

assembler.setVariables(ans)

# Create the function list
funcs = []

# Create the KS function
ksweight = 100.0
for i in range(1):
    funcs.append(functions.KSFailure(assembler, ksweight))

func_vals = assembler.evalFunctions(funcs)

# Set the element flag
flag = (TACS.OUTPUT_CONNECTIVITY |
        TACS.OUTPUT_NODES |
        TACS.OUTPUT_DISPLACEMENTS |
        TACS.OUTPUT_STRAINS |
        TACS.OUTPUT_STRESSES)
f5 = TACS.ToFH5(assembler, TACS.PLANE_STRESS_ELEMENT, flag)
f5.writeToFile('d2_test.f5')
