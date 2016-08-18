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
JA    = np.array([1.0,  0.0, 0.0,                          # second moment of mass [kg.m2]
                        2.0, 0.0,
                             0.75])

# Setup dynamics for body A
g     = elements.GibbsVector(np.array([0.0, 0.0, -10.0]))  # acceleration due to gravity [kg.m.s-2]
rinitA   = elements.GibbsVector(np.array([0., 0., 0.]))    # initial position wrt to global frame [m]
vinitA   = elements.GibbsVector(np.array([1., 1., 1.]))    # initial velocity wrt to global frame [m/s]
winitA   = elements.GibbsVector(np.array([.1, 10., 0.]))   # initial angular velocity wrt to global frame [rad/s]

# Create the body
bodyA = elements.RigidBody(frameA,
                           massA, cA, JA,
                           g, rinitA, vinitA, winitA)
 
#===================================#
# Create a test element for body A
#===================================#

print(">> Testing element...")

for dh in [1.0e-5]:
    test = elements.TestElem(bodyA)
    test.setStepSize(dh)
    test.setPrintLevel(2)
    test.testResidual()
    test.testJacobian()
    
#############################
# Create an instance of TACS
#############################

print(">> Creating TACS...")
num_nodes           = 1
num_owned_nodes     = 1
vars_per_node       = 8 # should equal num_displacements
num_elems           = 1
max_csr             = 2

num_dependent_nodes = None

tacs = TACS.Assembler.getInstance(comm, num_owned_nodes,
                                  vars_per_node, num_elems,
                                  num_nodes, max_csr)
                          
# Setup nodes and their connectivities
node = np.arange(num_nodes, dtype=np.intc) # 0, 1, 2
conn = np.arange(num_nodes, dtype=np.intc) # 0, 1, 2

tacs.addNodes(node, conn)
tacs.addElement(bodyA, conn[0:], 1)

tacs.finalize()

#################################
# Parameters for time integration
#################################

print(">> Setting up time integration")

tinit             = 0.0
tfinal            = 1.0

num_steps_per_sec = 10

# Create an FH5 object for teceplot output
flag = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS)
f5 = TACS.ToFH5(tacs, TACS.PY_RIGID, flag)

# Create functions of interest for adjoint
print(">> NBG Integration")
# NBG Integrator
nbg = TACS.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, 2)
nbg.setPrintLevel(2)
nbg.setJacAssemblyFreq(1)
nbg.configureOutput(f5, 10, "freefall/solution%04d.f5")
nbg.integrate()
