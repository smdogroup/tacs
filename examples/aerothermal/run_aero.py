# from pyfuntofem.model import *
#from pyfuntofem.driver import *
#from pyfuntofem.su2_interface import SU2Interface
from mpi4py import MPI
import numpy as np
# Import SU2
import pysu2
from tacs import TACS, elements, constitutive
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
mesh.scanBDFFile('plate.bdf')
# Set the element
mesh.setElement(0, element)
# Create the assembler object
varsPerNode = heat.getVarsPerNode()
assembler = mesh.createTACS(varsPerNode)
ndim = 2
SU2_CFD_ConfigFile = 'aero_plate_NS.cfg'
su2 = pysu2.CFluidDriver(SU2_CFD_ConfigFile, 1, ndim, False, comm)
tags = su2.GetAllBoundaryMarkersTag()
plate_index = tags.index('plate')
X = []
for i in range(su2.GetNumberVertices(plate_index)):
    X.extend([su2.GetVertexCoordX(plate_index, i),
              su2.GetVertexCoordY(plate_index, i),
              su2.GetVertexCoordZ(plate_index, i)])
Xpts = assembler.createNodeVec()
assembler.getNodes(Xpts)
Xpts_array = Xpts.getArray()
# Find the mapping from the thermal to the structural nodes
mapping = []
for i in range(0, len(X), 3):
    closest = 0
    min_dist = np.dot(X[i:i+3] - Xpts_array[:3],
                      X[i:i+3] - Xpts_array[:3])
    for j in range(0, len(Xpts_array), 3):
        dist = np.dot(X[i:i+3] - Xpts_array[j:j+3],
                      X[i:i+3] - Xpts_array[j:j+3])
        if dist < min_dist:
            min_dist = dist
            closest = j
    mapping.append(closest//3)
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
for iteration in range(20):
    for itr in range(itr_count, itr_count + 10):
        su2.PreprocessExtIter(itr)
        su2.Run()
        su2.Monitor(itr)
        su2.Output(itr)
    itr_count += 10
    for i in range(npts):
        normal_flux[i] = su2.GetVertexNormalHeatFlux(plate_index, i)
    res.zeroEntries()
    res_array = res.getArray()
    for i, j in enumerate(mapping):
        res_array[j] = normal_flux[i]
        assembler.setBCs(res)
    gmres.solve(res, ans)
    assembler.setVariables(ans)
    ans_array = ans.getArray()
    for i in range(npts):
        theta[i] = ans_array[j]
        su2.SetVertexTemperature(plate_index, i, theta[i])
# Set the element flag
flag = (TACS.OUTPUT_CONNECTIVITY |
        TACS.OUTPUT_NODES |
        TACS.OUTPUT_DISPLACEMENTS |
        TACS.OUTPUT_STRAINS)
f5 = TACS.ToFH5(assembler, TACS.SCALAR_2D_ELEMENT, flag)
f5.writeToFile('tacs_output.f5')
