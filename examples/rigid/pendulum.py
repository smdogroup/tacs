# Python script to simulate the dynamics of pendulum
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements

# TACS communicator
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print("Num procs :", size, "Rank :", rank)

#=========================================#
# Create vectors for dynamics of pendulum
#=========================================#

print('>> Creating Body A...')

# Setup reference frame for body A
rOA    = elements.GibbsVector(np.array([0., 0., 0.]))      # base point of body A from the global frame
rIA    = elements.GibbsVector(np.array([1., 0., 0.]))      # first coordinate direction wrt the base point
rJA    = elements.GibbsVector(np.array([0., 1., 0.]))      # second coordinate direction wrt the base point

frameA = elements.RefFrame(rOA, rIA, rJA)                  # reference frame attached to body A

# Setup mass properties for body A
massA = 4.0                                                # mass [kg]
cA    = np.array([0.0, 0.0, 0.0])                          # first moment of mass [kg.m]
JA    = np.array([1.0,  0.0, 0.0,                         # second moment of mass [kg.m2]
                        2.0, 0.0,
                             0.75])

# Setup dynamics for body A
g     = elements.GibbsVector(np.array([0.0, 0.0, -9.81]))  # acceleration due to gravity [kg.m.s-2]
r0A   = elements.GibbsVector(np.array([0., 0., 0.]))       # initial position wrt to global frame [m]
v0A   = elements.GibbsVector(np.array([0., 0.1, 0.]))      # initial velocity wrt to global frame [m/s]
w0A   = elements.GibbsVector(np.array([0., 0., 0.2]))      # initial angular velocity wrt to global frame [rad/s]

# Create the body
bodyA = elements.RigidBody(frameA, massA, cA, JA, g, r0A, v0A, w0A)
 
print('>> Creating Body B...')

# Setup reference frame for body B
rOB    = elements.GibbsVector(np.array([0., 0., 0.]))      # base point of body A from the global frame
rIB    = elements.GibbsVector(np.array([1., 0., 0.]))      # first coordinate direction wrt the base point
rJB    = elements.GibbsVector(np.array([0., 1., 0.]))      # second coordinate direction wrt the base point

frameB = elements.RefFrame(rOB, rIB, rJB)                  # reference frame attached to body A

# Setup mass properties for body B
massB = 4.0                                                # mass [kg]
cB    = np.array([0.0, 0.0, 0.0])                          # first moment of mass [kg.m]
JB    = np.array([1.0,  0.0,  0.0,                         # second moment of mass [kg.m2]
                        2.0,  0.0,
                             0.75])

# Setup dynamics for body B
g     = elements.GibbsVector(np.array([0.0, 0.0, -9.81]))   # acceleration due to gravity [kg.m.s-2]
r0B   = elements.GibbsVector(np.array([5., 0., 5.]))        # initial position wrt to global frame [m]
v0B   = elements.GibbsVector(np.array([0., 0.1, 0.]))       # initial velocity wrt to global frame [m/s]
w0B   = elements.GibbsVector(np.array([0., 0., 0.2]))       # initial angular velocity wrt to global frame [rad/s]

# Create the body
bodyB = elements.RigidBody(frameB, massB, cB, JB, g, r0B, v0B, w0B)

################################################
# Create a revolute joint between bodies A and B
################################################

xA  = elements.GibbsVector(np.array([5.0, 0.0, 0.0])) # r(a->joint)
xB  = elements.GibbsVector(np.array([0.0, 0.0, 5.0])) # r(b->joint)
eA  = elements.GibbsVector(np.array([1.0, 0.0, 0.0])) # ?
eB1 = elements.GibbsVector(np.array([1.0, 0.0, 0.0])) # ?
eB2 = elements.GibbsVector(np.array([1.0, 0.0, 0.0])) # ?

spjoint = elements.SphericalConstraint(xA, xB)
# rejoint = elements.RevoluteConstraint(xA, xB, eA, eB1, eB2)

#===================================#
# Create a test element for body A
#===================================#

print(">> Testing element...")

test = elements.TestElem(bodyA)
test.setPrintLevel(0)
test.testResidual()
for i in range(10):
    test.testJacobian(i)

#############################
# Create an instance of TACS
#############################

print(">> Creating TACS...")
num_nodes           = 3
num_owned_nodes     = 3
vars_per_node       = 8
num_elems           = 3
max_csr             = 2

num_dependent_nodes = None

tacs = TACS.Assembler.getInstance(comm, num_owned_nodes,
                                  vars_per_node, num_elems,
                                  num_nodes, max_csr)
                          
# Setup nodes and their connectivities
node = np.arange(num_nodes, dtype=np.intc)
conn = np.arange(num_nodes, dtype=np.intc)

tacs.addNodes(node, conn)

tacs.addElement(bodyA, conn[0:], 1)
tacs.addElement(bodyB, conn[1:], 1)
tacs.addElement(spjoint, conn[0:], 3)

tacs.finalize()

#################################
# Parameters for time integration
#################################

print(">> Setting up time integration")

tinit             = 0.0
tfinal            = 0.2
num_steps_per_sec = 10

# Create an FH5 object for teceplot output
flag = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS)
f5 = TACS.ToFH5(tacs, TACS.PY_RIGID, flag)

# Create functions of interest for adjoint

# Run the problem on the integrators
print(">> BDF Integration")
# BDF Integrator
for bdf_order in [1,2,3]:
    bdf = TACS.BDFIntegrator(tacs, tinit, tfinal,
                             num_steps_per_sec,
                             bdf_order)
    bdf.setPrintLevel(1)
    #bdf.setUseLapack(1)
    bdf.setJacAssemblyFreq(1)
    bdf.configureOutput(f5, 1, "pendulum%04d.f5")
    bdf.integrate()

## print(">> DIRK Integration")
## # DIRK Integrator
## for stages in [1,2,3]:
##     dirk = TACS.DIRKIntegrator(tacs, tinit, tfinal,
##                                num_steps_per_sec,
##                                stages)
##     dirk.setPrintLevel(1)
##     #dirk.setUseLapack(1)
##     dirk.setJacAssemblyFreq(1)
##     dirk.integrate()
