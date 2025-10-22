import numpy as np
from mpi4py import MPI

from tacs import TACS
from tacs import elements
from tacs import constitutive
from tacs import functions

# number of beam elements
nA = 4
nB = 8
nC = 4

# create the gravity vector
gravity = elements.GibbsVector(0.0, 0.0, -9.81)

# Set the points B, C and D
ptB = elements.GibbsVector(0.00, 0.12, 0.00)
ptC = elements.GibbsVector(0.24, 0.12, 0.00)
ptD = elements.GibbsVector(0.24, 0.00, 0.00)

# Create the revolute direction for B and D
revDirA = elements.GibbsVector(0.0, 0.0, 1.0)
revDirB = elements.GibbsVector(0.0, 0.0, 1.0)
theta = 5.0 * np.pi / 180.0
revDirC = elements.GibbsVector(np.sin(theta), 0.0, np.cos(theta))
revDirD = elements.GibbsVector(0.0, 0.0, 1.0)

# Create the revolute driver and constraints
omega = -0.6  # rad/s
fixed_point = 1
not_fixed = 0
revDriverA = elements.RevoluteDriver(revDirA, omega)
revB = elements.RevoluteConstraint(ptB, revDirB, not_fixed)
revC = elements.RevoluteConstraint(ptC, revDirC, not_fixed)
revD = elements.RevoluteConstraint(ptD, revDirD, fixed_point)

# Set the reference axes for each beam
axis_A = np.array([-1.0, 0.0, 0.0])
axis_B = np.array([0.0, 1.0, 0.0])
axis_C = np.array([1.0, 0.0, 0.0])

# Create the stiffness objects for each beam
rhoA_A = 1.997
rhoIy_A = 42.60e-6
rhoIz_A = 42.60e-6
rhoIyz_A = 0.0
EA_A = 52.99e6
GJ_A = 733.5
EIy_A = 1131.0
EIz_A = 1131.0
kGAy_A = 16.88e6
kGAz_A = 16.88e6
stiffA = constitutive.TimoshenkoConstitutive(
    rhoA_A, rhoIy_A, rhoIz_A, rhoIyz_A, EA_A, GJ_A, EIy_A, EIz_A, kGAy_A, kGAz_A, axis_A
)
stiffB = constitutive.TimoshenkoConstitutive(
    rhoA_A, rhoIy_A, rhoIz_A, rhoIyz_A, EA_A, GJ_A, EIy_A, EIz_A, kGAy_A, kGAz_A, axis_B
)
rhoA_B = 0.4992
rhoIy_B = 2.662e-6
rhoIz_B = 2.662e-6
rhoIyz_B = 0.0
EA_B = 13.25e6
GJ_B = 45.84
EIy_B = 70.66
EIz_B = 70.66
kGAy_B = 4.220e6
kGAz_B = 4.220e6
stiffC = constitutive.TimoshenkoConstitutive(
    rhoA_B, rhoIy_B, rhoIz_B, rhoIyz_B, EA_B, GJ_B, EIy_B, EIz_B, kGAy_B, kGAz_B, axis_C
)

# Create beams
beamA = elements.MITCBeam(stiffA, gravity)
beamB = elements.MITCBeam(stiffB, gravity)
beamC = elements.MITCBeam(stiffC, gravity)

# Create TACS
ndofpernode = 8
nnodes = (2 * nA + 1) + (2 * nB + 1) + (2 * nC + 1) + 4
nelems = nA + nB + nC + 4
assembler = TACS.Assembler.create(MPI.COMM_WORLD, ndofpernode, nnodes, nelems)

# create element nodes
X = assembler.createNodeVec()
xpts = X.getArray()
nodesA = []
nodesB = []
nodesC = []
n = 0
for i in range(2 * nA + 1):
    nodesA.append(n)
    xpts[3 * n + 1] = 0.12 * i / (2.0 * nA)
    n = n + 1
for i in range(2 * nB + 1):
    nodesB.append(n)
    xpts[3 * n + 0] = 0.24 * i / (2.0 * nB)
    xpts[3 * n + 1] = 0.12
    n = n + 1
for i in range(2 * nC + 1):
    nodesC.append(n)
    xpts[3 * n + 0] = 0.24
    xpts[3 * n + 1] = 0.12 * (1.0 - i / (2.0 * nC))
    n = n + 1

# nodes, elements, connectivities and pointers into connectivities
elems = []
ptr = []
conn = []
enum = 0
ptr.append(enum)
for i in range(nA):
    conn.append(nodesA[2 * i])
    conn.append(nodesA[2 * i + 1])
    conn.append(nodesA[2 * i + 2])
    elems.append(beamA)
    ptr.append(ptr[enum] + 3)
    enum = enum + 1

for i in range(nB):
    conn.append(nodesB[2 * i])
    conn.append(nodesB[2 * i + 1])
    conn.append(nodesB[2 * i + 2])
    elems.append(beamB)
    ptr.append(ptr[enum] + 3)
    enum = enum + 1

for i in range(nC):
    conn.append(nodesC[2 * i])
    conn.append(nodesC[2 * i + 1])
    conn.append(nodesC[2 * i + 2])
    elems.append(beamC)
    ptr.append(ptr[enum] + 3)
    enum = enum + 1

# Revolute driver
conn.append(nodesA[0])
conn.append(nnodes - 4)
elems.append(revDriverA)
ptr.append(ptr[enum] + 2)
enum = enum + 1

# revolute constraint
conn.append(nodesA[2 * nA])
conn.append(nodesB[0])
conn.append(nnodes - 3)
elems.append(revB)
ptr.append(ptr[enum] + 3)
enum = enum + 1

# revolute constraint
conn.append(nodesC[0])
conn.append(nodesB[2 * nB])
conn.append(nnodes - 2)
elems.append(revC)
ptr.append(ptr[enum] + 3)
enum = enum + 1

# revolute constraint
conn.append(nodesC[2 * nC])
conn.append(nnodes - 1)
elems.append(revD)
ptr.append(ptr[enum] + 2)
enum = enum + 1

# finalize the creation of assembler
ptr = np.array(ptr, dtype=np.intc)
conn = np.array(conn, dtype=np.intc)

assembler.setElementConnectivity(ptr, conn)
assembler.setElements(elems)
assembler.initialize()

assembler.setNodes(X)

# Create the functions to evaluate on the dynamic system
fmass = functions.StructuralMass(assembler)

ksRho = 100.0
ffail = functions.KSFailure(assembler, ksRho)

funcs = []
funcs.append(fmass)
funcs.append(ffail)


# create Integrator
t0 = 0.0
tf = 12.0
num_steps = 1200
order = 2
integrator = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
integrator.setPrintLevel(1)
integrator.setFunctions(funcs)
integrator.integrate()

fvals = integrator.evalFunctions(funcs)
print("mass=", fvals[0], "fail=", fvals[1])
