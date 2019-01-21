from mpi4py import MPI
from tacs import TACS, elements
import numpy as np
from spring_mass_damper import SpringMassDamper

# Create a SMD element
num_nodes = 1
num_disps = 1
m = 0.5
c = 6.0
k = 18.0
spr = SpringMassDamper(num_nodes, num_disps, m, c, k)

# Create TACS using the element
assembler = TACS.Assembler.create(MPI.COMM_WORLD, 1, 1, 1)
conn = np.array([0], dtype=np.intc)
ptr = np.array([0, 1], dtype=np.intc)
assembler.setElementConnectivity(conn, ptr)
assembler.setElements([spr])
assembler.initialize()

# Time marching setup
tinit  = 0.0
tfinal = 1000.0

# Create integrators for implicit time marching of system
sizes = [1250, 2500, 5000, 10000]
for nsteps in sizes:
    # BDF solution
    bdf_orders = [1,2,3]
    for order in bdf_orders:
        bdf = TACS.BDFIntegrator(assembler, tinit, tfinal, nsteps, order)
        bdf.setPrintLevel(0)
        bdf.integrate()
        bdf.writeRawSolution('smd-bdf' + str(order) + '-' + str(nsteps) + '.dat', 0)
        
    # DIRK solution
    dirk_orders = [2,3,4]
    for order in dirk_orders:
        dirk = TACS.DIRKIntegrator(assembler, tinit, tfinal, nsteps, order-1)
        dirk.setPrintLevel(0)
        dirk.integrate()
        dirk.writeRawSolution('smd-dirk' + str(order) + '-' + str(nsteps) + '.dat',0)
