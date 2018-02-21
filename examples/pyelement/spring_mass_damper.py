from mpi4py import MPI
from tacs import TACS, elements
import numpy as np

class SpringMassDamper(elements.pyElement):
    def __init__(self, num_disp, num_nodes, m, c, k):
        super(SpringMassDamper, self).__init__(num_disp, num_nodes)
        self.m = m
        self.c = c
        self.k = k    

    def getInitConditions(self, v, dv, ddv, xpts):
        v[0] = 1.0

        return

    def addResidual(self, time, res, X, v, dv, ddv):
        res[0] += self.m*ddv[0] + self.c*dv[0] + self.k*v[0]

        return    

    def addJacobian(self, time, J, alpha, beta, gamma, X, v, dv, ddv):
        J[0] += alpha*self.k + beta*self.c + gamma*self.m

        return

spr = SpringMassDamper(1, 1, 1.0, 0.5, 5.0)     

comm = MPI.COMM_WORLD
assembler = TACS.Assembler.create(comm, 1, 1, 1)

conn = np.array([0], dtype=np.intc)
ptr = np.array([0, 1], dtype=np.intc)

assembler.setElementConnectivity(conn, ptr)
assembler.setElements([spr])
assembler.initialize()
bdf = TACS.BDFIntegrator(assembler, 0.0, 100.0, 1000, 2)

bdf.integrate()
bdf.writeRawSolution('spring.dat', 0)
