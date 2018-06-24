'''
Demonstrates how to solve the adjoint equations with user-defined element

The objective function f(x) is the value of the displacement at the final time.
The design variable is simply the stiffness k of the spring-mass-damper system.
We will compute the gradient df/dk by solving the adjoint equations for the
user-defined element, and then check the gradient via finite differences.

'''
from mpi4py import MPI
from tacs import TACS, elements, functions
import numpy as np

np.set_printoptions(precision=15)

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

        return

    def addResidual(self, time, res, X, v, dv, ddv):
        '''Add the residual of the governing equations'''
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]

        return    

    def addJacobian(self, time, J, alpha, beta, gamma, X, v, dv, ddv):
        '''Add the Jacobian of the governing equations'''
        J[0] += alpha*self.k + beta*self.c + gamma*self.m

        return

    def addAdjResProduct(self, time, scale, dvSens, psi, X, v, dv, ddv):
        '''
        Add the derivative of the product of the adjoint variables and the
        residuals to the design variable sensitivities
        '''
        dvSens += scale*psi[0]*v[0]

    def updateStiffness(self, k):
        self.k = k

        return

# Create instance of user-defined element
num_nodes = 1
num_disps = 1
m = 1.0
c = 0.05
k = 5.0
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
dt = 0.01
num_steps = 2000
tf = num_steps*dt
order = 2
bdf = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
bdf.setPrintLevel(0) # turn off printing

# Integrate governing equations
bdf.iterate(0)
for step in range(1,num_steps+1):
    bdf.iterate(step)
_, uvec, _, _ = bdf.getStates(num_steps)
u = uvec.getArray().copy()
bdf.writeRawSolution('spring.dat', 0)

# Specify the number of design variables and the function to the integrator
# (use Structural Mass as a dummy function)
num_dvs = 1
funclist = [functions.StructuralMass(assembler)]
bdf.setFunctions(funclist, num_dvs)

# Solve the adjoint equations and compute the gradient
dfdu_vec = assembler.createVec()
for step in range(num_steps, -1, -1):
    bdf.initAdjoint(step)
    dfdu = dfdu_vec.getArray()
    if step == num_steps:
        dfdu[0] = -1.0
    else:
        dfdu[0] = 0.0
    bdf.iterateAdjoint(step, [dfdu_vec])
    bdf.postAdjoint(step)

dfdx = np.array([0.0], dtype=TACS.dtype)
bdf.getGradient(dfdx)

# Check the gradient by finite difference or complex step
if TACS.dtype == complex:
    h = 1.0e-30
    spr.updateStiffness(k + 1j*h)
else:
    h = 1.0e-6
    spr.updateStiffness(k + h)

bdf.iterate(0)
for step in range(1,num_steps+1):
    bdf.iterate(step)
bdf.writeRawSolution('spring.dat', 0)
_, upos_vec, _, _ = bdf.getStates(num_steps)
upos = upos_vec.getArray().copy()

if TACS.dtype == complex:
    approx = upos.imag/h
else:
    spr.updateStiffness(k - h)
    bdf.iterate(0)
    for step in range(1,num_steps+1):
        bdf.iterate(step)
    bdf.writeRawSolution('spring.dat', 0)
    _, uneg_vec, _, _ = bdf.getStates(num_steps)
    uneg = uneg_vec.getArray().copy()

    approx = 0.5*(upos - uneg)/h

rel_error = (dfdx.real - approx)/approx

print "f = ", u[0].real
print "df/dx =         ", dfdx[0].real
print "df/dx, approx = ", approx[0]
print "rel. error =    ", rel_error[0]
