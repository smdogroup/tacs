from mpi4py import MPI
from tacs import TACS, elements
import numpy as np

# Define an element in TACS using the pyElement feature
class SpringMassDamper(elements.pyElement):
    def __init__(self, num_nodes, num_disps, m, c, k, dt):
        super(SpringMassDamper, self).__init__(num_disp, num_nodes)
        self.m = m
        self.c = c
        self.k = k    
        self.dt = dt

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

    def addAdjResProduct(self, time, scale, dvSens, psi, 
                         Xpts, vars, dvars, ddvars):
        '''
        Add the derivative of the product of the adjoint variables and the
        residuals to the design variable sensitivities
        '''

# Create instance of user-defined element
num_nodes = 1
num_disps = 1
m = 1.0
c = 1.0
k = 5.0
dt = 0.1
spr = SpringMassDamper(num_nodes, num_disps, m, c, k, dt)     

# Create TACS and add user-defined element
comm = MPI.COMM_WORLD
assembler = TACS.Assembler.create(comm, 1, 1, 1)

conn = np.array([0], dtype=np.intc)
ptr = np.array([0, 1], dtype=np.intc)

assembler.setElementConnectivity(conn, ptr)
assembler.setElements([spr])
assembler.initialize()

# Create instance of integrator and integrate governing equations
bdf = TACS.BDFIntegrator(assembler, 0.0, 100.0, 1000, 2)
bdf.iterate(0)

fvec = assembler.createVec()
for step in range(1,1001):
    bdf.iterate(step, fvec)
bdf.writeRawSolution('spring.dat', 0)
