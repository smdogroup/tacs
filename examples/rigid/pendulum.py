# Python script to simulate the dynamics of pendulum
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements, solver

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
rOA    = elements.GibbsVector(np.array([0., 0., 0.]))       # base point of body A from the global frame
rIA    = elements.GibbsVector(np.array([1., 0., 0.]))       # first coordinate direction wrt the base point
rJA    = elements.GibbsVector(np.array([0., 1., 0.]))       # second coordinate direction wrt the base point
frameA = elements.RefFrame(rOA, rIA, rJA)                   # reference frame attached to body A

# Setup mass properties for body A
mass = 4.0                                                # mass [kg]
c    = np.array([0.5, 0.2, -0.1])                         # first moment of mass [kg.m]
J    = np.array([1.0, -0.2, -0.25,                        # second moment of mass [kg.m2]
                 2.0, 0.1,
                 0.75])

# Setup dynamics for body A
g    = elements.GibbsVector(np.array([0.0, 0.0, -9.81]))  # acceleration due to gravity [kg.m.s-2]
r0   = elements.GibbsVector(np.array([0., 0., 0.]))       # initial position wrt to global frame [m]
v0   = elements.GibbsVector(np.array([0., 0.1, 0.]))      # initial velocity wrt to global frame [m/s]
w0   = elements.GibbsVector(np.array([0., 0., 0.2]))      # initial angular velocity wrt to global frame [rad/s]

# Create the body
bodyA = elements.RigidBody(frameA, mass, c, J, g, r0, v0, w0)

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
num_nodes           = 1
num_owned_nodes     = 1
vars_per_node       = 8
num_elems           = 1
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

tacs.finalize()

#################################
# Parameters for time integration
#################################

print(">> Setting up time integration")

tinit             = 0.0
tfinal            = 0.1
num_steps_per_sec = 10

# Run the problem on the integrators
for bdf_order in [1,2,3]:
    # BDF Integrator
    bdf = solver.BDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec, bdf_order)
    bdf.setPrintLevel(1)
    bdf.setJacAssemblyFreq(1)
    bdf.integrate()

for stages in [1,2,3]:
    # BDF Integrator
    bdf = solver.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, stages)
    bdf.setPrintLevel(1)
    bdf.setJacAssemblyFreq(1)
    bdf.integrate()
    
#    bdf.writeSolutionToF5()

## # Set the element flag
flag = (TACS.ToFH5.NODES |
        TACS.ToFH5.DISPLACEMENTS |
        TACS.ToFH5.STRAINS)
f5 = TACS.ToFH5(tacs, TACS.PY_PLANE_STRESS, flag)
f5.writeToFile('triangle_test.f5')
