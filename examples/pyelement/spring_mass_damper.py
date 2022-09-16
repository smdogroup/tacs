from mpi4py import MPI
from tacs import TACS, elements
import numpy as np
import sys
import matplotlib.pylab as plt

# Define an element in TACS using the pyElement feature
class SpringMassDamper(elements.pyElement):
    def __init__(self, num_disps, num_nodes, m, c, k):
        self.m = m
        self.c = c
        self.k = k

    def getMultiplierIndex(self):
        return -1

    def getInitConditions(self, index, X, v, dv, ddv):
        """Define the initial conditions"""
        v[0] = -0.5
        dv[0] = 1.0

        return

    def addResidual(self, index, time, X, v, dv, ddv, res):
        """Add the residual of the governing equations"""
        res[0] += self.m * ddv[0] + self.c * dv[0] + self.k * v[0]

        return

    def addJacobian(self, index, time, alpha, beta, gamma, X, v, dv, ddv, res, mat):
        """Add the Jacobian of the governing equations"""
        res[0] += self.m * ddv[0] + self.c * dv[0] + self.k * v[0]
        mat[0] += alpha * self.k + beta * self.c + gamma * self.m

        return


def createAssembler(m=1.0, c=0.5, k=5.0):
    num_disps = 1
    num_nodes = 1
    spr = SpringMassDamper(num_disps, num_nodes, m, c, k)

    # Add user-defined element to TACS
    comm = MPI.COMM_WORLD
    assembler = TACS.Assembler.create(comm, 1, 1, 1)

    ptr = np.array([0, 1], dtype=np.intc)
    conn = np.array([0], dtype=np.intc)
    assembler.setElementConnectivity(ptr, conn)
    assembler.setElements([spr])
    assembler.initialize()

    return assembler


if __name__ == "__main__":
    # Create instance of user-defined element
    m = 1.0
    c = 0.5
    k = 5.0
    tf = 10.0
    assembler = createAssembler(m=m, c=c, k=k)

    # Compute the exact displacement at the final time
    wn = np.sqrt(k / m)
    xi = c / (2 * m * wn)
    x0 = -0.5
    xdot0 = 1.0
    wd = wn * np.sqrt(1.0 - xi**2)
    A = np.sqrt(x0**2 + ((xdot0 + xi * wn * x0) / wd) ** 2)
    phi_d = np.arctan(x0 * wd / (xdot0 + xi * wn * x0))
    x_exact = A * np.exp(-xi * wn * tf) * np.sin(wd * tf + phi_d)

    all_error = []
    all_dt = []
    for order in [1, 2, 3]:
        all_error.append([])
        all_dt.append([])
        for num_steps in [100, 500, 1000, 5000, 10000, 50000]:
            t0 = 0.0
            if "bdf" in sys.argv:
                integrator = TACS.BDFIntegrator(assembler, t0, tf, num_steps, order)
            else:
                integrator = TACS.DIRKIntegrator(assembler, t0, tf, num_steps, order)

            # Integrate forward in time
            integrator.integrate()

            # Get the solution at the final time
            t, vec, dvec, ddvec = integrator.getStates(num_steps)

            sol = vec.getArray()[0]
            if TACS.dtype is complex:
                all_error[-1].append(np.fabs((sol.real - x_exact.real) / x_exact.real))
            else:
                all_error[-1].append(np.fabs((sol - x_exact) / x_exact))
            all_dt[-1].append(tf / num_steps)

    plt.figure()
    for err, dt in zip(all_error, all_dt):
        plt.loglog(dt, err)
    plt.show()
