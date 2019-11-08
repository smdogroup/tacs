# from pyfuntofem.model import *
#from pyfuntofem.driver import *
#from pyfuntofem.su2_interface import SU2Interface
from mpi4py import MPI
import numpy as np
import sys
# Import SU2
import pysu2
from tacs import TACS, elements, constitutive
from funtofem import TransferScheme

np.set_printoptions(threshold=sys.maxsize)

comm = MPI.COMM_WORLD

# Create the constitutvie propertes and model
props = constitutive.MaterialProperties()
con = constitutive.PlaneStressConstitutive(props)
heat = elements.HeatConduction2D(con)

# Create the basis class
quad_basis = elements.LinearQuadBasis()

# Create the element
element = elements.Element2D(heat, quad_basis)

# Load in the mesh
mesh = TACS.MeshLoader(comm)
mesh.scanBDFFile('plate_offset.bdf')

# Set the element
mesh.setElement(0, element)

# Create the assembler object
varsPerNode = heat.getVarsPerNode()
assembler = mesh.createTACS(varsPerNode)

# Create su2
ndim = 2
SU2_CFD_ConfigFile = 'aero_plate_NS.cfg'
su2 = pysu2.CFluidDriver(SU2_CFD_ConfigFile, 1, ndim, False, comm)
tags = su2.GetAllBoundaryMarkersTag()
plate_index = tags.index('plate')
X = []

# Initialze MELD
meld = TransferScheme.pyMELDThermal(MPI.COMM_SELF, MPI.COMM_SELF, 0, MPI.COMM_SELF, 0, 1, 10, 0.5) #axis of symmetry, num nearest, beta

for i in range(su2.GetNumberVertices(plate_index)):
    X.extend([su2.GetVertexCoordX(plate_index, i),
              su2.GetVertexCoordY(plate_index, i),
              su2.GetVertexCoordZ(plate_index, i)])

Xpts = assembler.createNodeVec()
assembler.getNodes(Xpts)
Xpts_array = Xpts.getArray()

# find structural nodes on surface of plate
# i.e. structural nodes that are actually exposed to flow
# for this case that's any node with y coord == 0
# for more complex examples tbd
plate_surface = []
for i in range(len(Xpts_array) // 3):
    if Xpts_array[3*i+1] == 0.0:
        plate_surface.extend(Xpts_array[3*i:3*i+3])

plate_surface = np.array(plate_surface)
X = np.array(X)
print(plate_surface)
print(X)

# give meld the aero nodes and plate surface nodes
# does this have to be done every iteration?
meld.setStructNodes(plate_surface)
meld.setAeroNodes(X)
meld.initialize()
    
# Create the vectors/matrices
res = assembler.createVec()
ans = assembler.createVec()
mat = assembler.createSchurMat()
pc = TACS.Pc(mat)

# Assemble the heat conduction matrix
assembler.assembleJacobian(1.0, 0.0, 0.0, res, mat)
pc.factor()
gmres = TACS.KSM(mat, pc, 20)

# Get the boundary temperature
npts = su2.GetNumberVertices(plate_index)
normal_flux = np.zeros(npts)
temp_ret = np.zeros(npts)
theta = np.zeros(npts)
itr_count = 0

# 20 iterations
for iteration in range(20):
    # run su2
    for itr in range(itr_count, itr_count + 10):
        su2.PreprocessExtIter(itr)
        su2.Run()
        su2.Monitor(itr)
        su2.Output(itr)
    itr_count += 10

    # get the normal heat flux from su2
    # this needs to be normalized by area!
    for i in range(npts):
        normal_flux[i] = su2.GetVertexNormalHeatFlux(plate_index, i)
        
    res.zeroEntries()
    res_array = res.getArray()

    # set flux into TACS
    meld.transferFlux(normal_flux, res_array)
    assembler.setBCs(res)

    # solve thermal problem
    gmres.solve(res, ans)
    assembler.setVariables(ans)

    # get temps
    ans_array = ans.getArray()
    meld.transferTemp(ans_array, theta)
    for i in range(npts):
        su2.SetVertexTemperature(plate_index, i, theta[i])
        
# Set the element flag
flag = (TACS.OUTPUT_CONNECTIVITY |
        TACS.OUTPUT_NODES |
        TACS.OUTPUT_DISPLACEMENTS |
        TACS.OUTPUT_STRAINS)
f5 = TACS.ToFH5(assembler, TACS.SCALAR_2D_ELEMENT, flag)
f5.writeToFile('tacs_output_meldtest.f5')
