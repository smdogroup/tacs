# Import numpy
import numpy as np

# Import MPI
from mpi4py import MPI

# Import TACS and assorted repositories
from tacs import TACS, elements, constitutive, functions

class PS(constitutive.pyPlaneStress):
    def __init__(self, rho, E, nu, ys):
        self.rho = rho
        self.nu = nu
        self.D = E/(1.0 - nu**2)
        self.G = 0.5*E/(1.0 + nu)
        self.ys = ys
        return

    def calculateStress(self, pt, e):
        '''Compute the stress at the point'''

        # Compute the stresses
        s = np.zeros(e.shape)
        s[0] = self.D*(e[0] + self.nu*e[1])
        s[1] = self.D*(e[1] + self.nu*e[0])
        s[2] = self.G*e[2]
        return s

    def getPointwiseMass(self, pt):
        '''Return the pointwise mass'''
        return self.rho

    def failure(self, pt, e):
        '''Evaluate the stress'''

        # Compute the stresses
        s = np.zeros(e.shape)
        s[0] = self.D*(e[0] + self.nu*e[1])
        s[1] = self.D*(e[1] + self.nu*e[0])
        s[2] = self.G*e[2]

        # Compute the valure of the failure function
        fval = (s[0]**2 + s[1]**2 - s[0]*s[1] + 3*s[2]**2)/self.ys**2 

        return fval

# Allocate the TACS creator
comm = MPI.COMM_WORLD
creator = TACS.Creator(comm, 2)

# Create the stiffness object
rho = 2570.0
E = 70e9
nu = 0.3
ys = 350e6

# stiff = constitutive.PlaneStress(rho, E, nu)
stiff = PS(rho, E, nu, ys) 

# Create the elements
elem_order = 2
elem = elements.PlaneTri6(stiff)
    
if comm.rank == 0:
    # Create the elements
    nx = 25
    ny = 25
    
    # Set the nodes
    nnodes = (2*nx+1)*(2*ny+1)
    nelems = 2*nx*ny
    nodes = np.arange(nnodes).reshape((2*nx+1, 2*ny+1))
    
    conn = []
    for j in range(ny):
        for i in range(nx):
            # Append the first set of nodes
            conn.append([nodes[2*i, 2*j],
                         nodes[2*i+2, 2*j],
                         nodes[2*i+2, 2*j+2],
                         nodes[2*i+1, 2*j],
                         nodes[2*i+2, 2*j+1],
                         nodes[2*i+1, 2*j+1]])
            
            # Append the second set of nodes
            conn.append([nodes[2*i, 2*j+2],
                         nodes[2*i, 2*j],
                         nodes[2*i+2, 2*j+2],
                         nodes[2*i, 2*j+1],
                         nodes[2*i+1, 2*j+1],
                         nodes[2*i+1, 2*j+2]])

    # Set the node pointers
    conn = np.array(conn, dtype=np.intc).flatten()
    ptr = np.arange(0, 6*nelems+1, 6, dtype=np.intc)
    elem_ids = np.zeros(nelems, dtype=np.intc)
    creator.setGlobalConnectivity(nnodes, ptr, conn, elem_ids)

    # Set up the boundary conditions
    bcnodes = np.array(nodes[0,:], dtype=np.intc)

    # Set the boundary condition variables
    nbcs = 2*bcnodes.shape[0]
    bcvars = np.zeros(nbcs, dtype=np.intc)
    bcvars[:nbcs:2] = 0
    bcvars[1:nbcs:2] = 1

    # Set the boundary condition pointers
    bcptr = np.arange(0, nbcs+1, 2, dtype=np.intc)
    creator.setBoundaryConditions(bcnodes, bcvars, bcptr)

    # Set the node locations
    Xpts = np.zeros(3*nnodes)
    x = np.linspace(0, 10, 2*nx+1)
    y = np.linspace(0, 10, 2*nx+1)
    for j in range(2*ny+1):
        for i in range(2*nx+1):
            Xpts[3*nodes[i,j]] = x[i]
            Xpts[3*nodes[i,j]+1] = y[j]
            
    # Set the node locations
    creator.setNodes(Xpts)

# Set the elements
elements = [ elem ]
creator.setElements(elements)

# Create the tacs assembler object
tacs = creator.createTACS()

res = tacs.createVec()
ans = tacs.createVec()
mat = tacs.createFEMat()

# Create the preconditioner for the corresponding matrix
pc = TACS.Pc(mat)

# Assemble the Jacobian
alpha = 1.0
beta = 0.0
gamma = 0.0
tacs.assembleJacobian(alpha, beta, gamma, res, mat)
pc.factor()

res.setRand(1.0, 1.0)
tacs.applyBCs(res)
pc.applyFactor(res, ans)
ans.scale(-1.0)

tacs.setVariables(ans)

tacs.setSimulationTime(0.15)

# Create the function list
funcs = []

# Create the KS function
ksweight = 100.0
for i in range(1):
    funcs.append(functions.KSFailure(tacs, ksweight))

func_vals = tacs.evalFunctions(funcs)

# Set the element flag
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(tacs, TACS.PY_PLANE_STRESS, flag)
f5.writeToFile('triangle_test.f5')
