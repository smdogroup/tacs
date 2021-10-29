from tacs import TACS, elements, constitutive, functions
from mpi4py import MPI
import numpy as np

'''
Create TACS Assembler
'''
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size
nx = 5  # number of elements in x direction
ny = 5  # number of elements in y direction
Lx = 1.0
Ly = 1.0
varsPerNode = 3
nodesPerProc = (nx+1)*(ny+1)/size
elemsPerProc = nx*ny/size
numOwnerNodes = nodesPerProc
numElements = elemsPerProc
numDependentNodes = 0

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
int firstElem = rank*elemsPerProc
int firstNode = rank*nodesPerProc
int lastElem = (rank+1)*elemsPerProc
int lastNode = (rank+1)*nodesPerProc
if (rank == size-1):
    lastElem = nx*ny
    lastNode = (nx+1)*(ny+1)

# Populate connectivity
ptr = np.zeros(numElements + 1)
conn = np.zeros(4*numElements)
ptr[0] = 0
k = 0
for elem in range(firstElem, lastElem):
    i = elem % nx
    j = elem // nx
    conn[4*k] = i + j*(nx+1)
    conn[4*k+1] = i+1 + j*(nx+1)
    conn[4*k+2] = i + (j+1)*(nx+1)
    conn[4*k+3] = i+1 + (j+1)*(nx+1)
    ptr[k+1] = 4*(k+1)
    k += 1

# Set the connectivity
assembler.setElementConnectivity(ptr, conn)

# Create the isotropic material class
props = constitutive.MaterialProperties(rho=2700.0, E=70e3, nu=0.3, ys=270.0)

# Create basis, constitutive, element, etc
linear_basis = element.LinearQuadBasis()
stiff = constitutive.PlaneStressConstitutive(props)
elements = []
k = 0
for elem in range(firstElem, lastElem):
    stiff = constitutive.PlaneStressConstitutive(props, 1.0, elem)
    



