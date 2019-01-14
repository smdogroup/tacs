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
        v[0] = -0.5
        dv[0] = 1.0

        return

    def addResidual(self, time, res, X, v, dv, ddv):
        '''Add the residual of the governing equations'''
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]

        return    

    def addJacobian(self, time, J, alpha, beta, gamma, X, v, dv, ddv):
        '''Add the Jacobian of the governing equations'''
        J[0] += alpha*self.k + beta*self.c + gamma*self.m

        return

if __name__ == '__main__':    
    # Create instance of user-defined element
    num_nodes = 1
    num_disps = 1
    m = 1.0
    c = 0.5
    k = 5.0+1j*1e-30
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
    num_steps = 1000
    tf = num_steps*dt
    order = 2
    bdf = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
    
    # Integrate governing equations
    #bdf.integrate()
    bdf.iterate(0)
    for step in range(1,num_steps+1):
        bdf.iterate(step)
    
    _, uvec, _, _ = bdf.getStates(num_steps)
    u = uvec.getArray()
    print "f = ", u
    print "df/dx, approx = ", u.imag/1e-30
    
    # Write out solution
    bdf.writeRawSolution('spring.dat', 0)
