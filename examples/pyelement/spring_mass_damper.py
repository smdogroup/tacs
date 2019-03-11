import sys

for p in sys.path:
    print(p)

from mpi4py import MPI
from tacs import TACS, elements
import numpy as np

# Define an element in TACS using the pyElement feature
class SpringMassDamper(elements.pyElement):
    def __init__(self, num_nodes, num_disps, m, c, k):
        super(SpringMassDamper, self).__init__(num_disps, num_nodes)
        self.m = m
        self.c = c
        self.k = k    

    def getInitConditions(self, v, dv, ddv, xpts):
        '''Define the initial conditions'''
        v[0] = 1.0 
        dv[0] = 0.0
        ddv[0] = 0.0

        return

    def addResidual(self, time, res, X, v, dv, ddv):
        '''Add the residual of the governing equations'''
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]

        return    

    def addJacobian(self, time, J, alpha, beta, gamma, X, v, dv, ddv):
        '''Add the Jacobian of the governing equations'''
        J[0] += alpha*self.k + beta*self.c + gamma*self.m

        return

# Create instance of user-defined element
num_nodes = 1
num_disps = 1
m = 1.0
c = 0.5 #0.5
k = 5.0 #+1j*1e-30
spr = SpringMassDamper(num_nodes, num_disps, m, c, k)     

# Add user-defined element to TACS
comm = MPI.COMM_WORLD
assembler = TACS.Assembler.create(comm, 1, 1, 1)

conn = np.array([0], dtype=np.intc)
ptr = np.array([0, 1], dtype=np.intc)

assembler.setElementConnectivity(conn, ptr)
assembler.setElements([spr])
assembler.initialize()

# Create instance of integrator
t0 = 0.0
dt = 0.1
num_steps = 100 
tf = num_steps*dt

# order = 2 #BDF integrator order 
# bdf = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
# bdf.setRelTol(1e-9)
# bdf.setAbsTol(1e-12)

# stages = 3
# dirk = TACS.DIRKIntegrator(assembler, t0, tf, num_steps, stages)
# dirk.setRelTol(1e-9)
# dirk.setAbsTol(1e-12)

stages = 6
esdirk = TACS.ESDIRKIntegrator(assembler, t0, tf, num_steps, stages)
esdirk.setRelTol(1e-9)
esdirk.setAbsTol(1e-12)

# Integrate governing equations
# bdf.integrate()
# bdf.iterate(0)
# dirk.iterate(0)
esdirk.iterate(0)
for step in range(1,num_steps+1):
    # bdf.iterate(step)
    # dirk.iterate(step)
    # esdirk.iterate(step)
    for stage in range(stages):
        # dirk.iterateStage(step, stage)
        esdirk.iterateStage(step, stage)

# Get states at final time step
# _, uvec, _, _ = bdf.getStates(num_steps)
# _, uvec, _, _ = dirk.getStates(num_steps)
_, uvec, _, _ = esdirk.getStates(num_steps)

# Get states at specified stage of final time step
# _, qS, qdotS, qddotS = dirk.getStageStates(num_steps,1)

# Get values of states
u = uvec.getArray()
# stage_disp = qS.getArray()
# stage_vel = qdotS.getArray()
# stage_acc = qddotS.getArray()

# Print states
print "u_f = ", u
# print "qS = ", stage_disp
# print "qdotS = ", stage_vel
# print "qddotS = ", stage_acc
# print "df/dx, approx = ", u.imag/1e-30
# dirk.say_hello()

# Write out solution
# bdf.writeRawSolution('spring.dat', 0)
# dirk.writeRawSolution('spring.dat', 0)
# esdirk.writeRawSolution('spring.dat', 0)

#####
# Check error with exact solution at end time
wn = np.sqrt(k/m)
xi = c/(2*m*wn)
x0 = 1
xdot0 = 0
wd = wn*np.sqrt(1-xi**2)
A = np.sqrt(x0**2+((xdot0+xi*wn*x0)/wd)**2)
phi_d = np.arctan(x0*wd/(xdot0+xi*wn*x0))
x_exact = A*np.exp(-xi*wn*tf)*np.sin(wd*tf+phi_d)
error = np.abs((u-x_exact)/x_exact)
print "error = ", error