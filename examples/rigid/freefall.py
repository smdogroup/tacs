# Python script to simulate the dynamics of pendulum
import numpy as np
from mpi4py import MPI
from tacs import TACS, elements

# TACS communicator
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print("Num procs :", size, "Rank :", rank)

# =========================================#
# Create vectors for dynamics of pendulum
# =========================================#

print(">> Creating Body A...")

# Setup reference frame for body A
rOA = np.array([0.0, 0.0, 0.0])  # base point of body A from the global frame
rIA = np.array([1.0, 0.0, 0.0])  # first coordinate direction wrt the base point
rJA = np.array([0.0, 1.0, 0.0])  # second coordinate direction wrt the base point

# Inertial properties
massA = 4.0  # mass [kg]
cA = np.array([0.0, 0.0, 0.0])  # first moment of mass [kg.m]
JA = np.array([1.0, 0.0, 0.0, 1.0, 0.0, 1.0])  # second moment of mass [kg.m2]

# Setup dynamics for body A
rinitA = np.array([0.0, 0.0, 0.0])  # initial position wrt to global frame [m]
vinitA = np.array([0.0, 0.0, 0.0])  # initial velocity wrt to global frame [m/s]
winitA = np.array(
    [0.0, 1.0, 0.0]
)  # initial angular velocity wrt to global frame [rad/s]
g = np.array([0.0, 0.0, 0.0])  # acceleration due to gravity [kg.m.s-2]

# Create the body
bodyA = elements.RigidBody(rOA, rIA, rJA, massA, cA, JA, rinitA, vinitA, winitA, g)
## bodyA.setStepSize(dh)
## bodyA.setPrintLevel(2)
## bodyA.testResidual()
## bodyA.testJacobian()

#############################
# Create an instance of TACS
#############################

print(">> Creating TACS...")
num_nodes = 1
num_owned_nodes = 1
vars_per_node = 8  # should equal num_displacements
num_elems = 1
num_dependent_nodes = 0

tacs = TACS.Assembler.create(
    comm, vars_per_node, num_owned_nodes, num_elems, num_dependent_nodes
)

# Setup nodes and their connectivities
node = np.arange(num_nodes, dtype=np.intc)  # 0, 1, 2
conn = np.arange(num_nodes, dtype=np.intc)  # 0, 1, 2
ptr = np.array([0, 1], dtype=np.intc)

elemList = [bodyA]

tacs.setElements(elemList)
tacs.setElementConnectivity(conn, ptr)
tacs.initialize()

#################################
# Parameters for time integration
#################################

print(">> Setting up time integration")

tinit = 0.0
tfinal = 10.0

num_steps_per_sec = 25

# Create an FH5 object for teceplot output
# flag = (TACS.ToFH5.NODES | TACS.ToFH5.DISPLACEMENTS)
# f5 = TACS.ToFH5(tacs, TACS.PY_RIGID, flag)

# Create functions of interest for adjoint
print(">> NBG Integration")
# NBG Integrator
nbg = TACS.DIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec, 2)
nbg.setPrintLevel(2)
nbg.setJacAssemblyFreq(1)
nbg.setOutputFrequency(1)
nbg.integrate()
