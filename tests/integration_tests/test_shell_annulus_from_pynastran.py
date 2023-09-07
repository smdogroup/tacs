import numpy as np
from pyNastran.bdf.bdf import BDF

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, functions

"""
Create annulus out of shell elements using pyNASTRAN's BDF class.
We then use this class object to initialize pyTACS rather than reading in a file.
Apply a load at the spider-webbed RBE2 centered inside the annulus and test KSFailure, StructuralMass, 
Moment of Inertia, and Compliance functions and sensitivities.
"""

# Geometric properties
Ri = 0.5
Ro = 1.0
z = 0.0
t = 0.05
# Material properties
E = 70e9
nu = 0.3
ys = 270e6
rho = 2700.0

# number of elements in r/theta
nr = 5
ntheta = 10
# Center node ID
center_nid = 1000

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "center_moment_Izz": 174.1393226319456,
        "center_moment_compliance": 2047161.6725086374,
        "center_moment_ks_vmfailure": 15.179071361242388,
        "center_moment_mass": 297.56628397306457,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 5e-6
            self.atol = 1e-8
            self.dh = 1e-10
        else:
            self.rtol = 1e-3
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        bdf_info = self.createBDFObject()

        fea_assembler = pytacs.pyTACS(bdf_info, comm)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Apply a static load case with a x moment at the center
        problem = fea_assembler.createStaticProblem("center_moment")
        problem.addLoadToNodes(
            center_nid, [0.0, 0.0, 0.0, 1e7, 0.0, 0.0], nastranOrdering=True
        )

        problem.addFunction("mass", functions.StructuralMass)
        problem.addFunction(
            "Izz", functions.MomentOfInertia, direction1=[0, 0, 1], direction2=[0, 0, 1]
        )
        problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
        problem.addFunction("compliance", functions.Compliance)
        # Set convergence tol to be tight
        problem.setOption("L2Convergence", 1e-20)
        problem.setOption("L2ConvergenceRel", 1e-20)

        return [problem], fea_assembler

    def createBDFObject(self):
        # Initialize our BDF object
        bdfInfo = BDF()
        # Create a cylindrical coord system
        bdfInfo.add_cord2c(1, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 1.0])
        # Discretize in r and theta
        r_array = np.linspace(Ri, Ro, nr)
        theta_array = np.linspace(0.0, 360.0, ntheta, endpoint=False)

        # Loop through and add each node
        curID = 1
        nodeIDs = {}
        for r in r_array:
            for theta in theta_array:
                bdfInfo.add_grid(curID, [r, theta, z], 1)
                nodeIDs[r, theta] = curID
                curID += 1

        # Loop through and add each element
        curID = 1
        for i in range(nr - 1):
            for j in range(ntheta):
                r0 = r_array[i]
                t0 = theta_array[j - 1]
                r1 = r_array[i + 1]
                t1 = theta_array[j]
                conn = [
                    nodeIDs[r0, t0],
                    nodeIDs[r0, t1],
                    nodeIDs[r1, t1],
                    nodeIDs[r1, t0],
                ]
                bdfInfo.add_cquad4(curID, 1, conn)
                curID += 1

        # Simply support the outer edge
        for theta in theta_array:
            bdfInfo.add_spc(0, nodeIDs[Ro, theta], "123", 0.0)

        # Add node at center for applying load
        bdfInfo.add_grid(center_nid, [0.0, 0.0, 0.0], 0)

        # Connect center node with RBE to inner edge
        rbe_dep_nodes = []
        for theta in theta_array:
            rbe_dep_nodes.append(nodeIDs[Ri, theta])
        bdfInfo.add_rbe2(curID, center_nid, "123456", rbe_dep_nodes)

        # Define material card for Aluminum
        bdfInfo.add_mat1(1, E, None, nu, rho, St=ys)
        # Define shell thickness
        bdfInfo.add_pshell(1, 1, t)

        return bdfInfo
