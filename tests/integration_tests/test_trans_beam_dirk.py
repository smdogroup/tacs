import numpy as np
import os
from tacs import pytacs, functions
from pytacs_analysis_base_test import PyTACSTestCase

"""
The test case features a cantilevered beam with a time-varying out-of-plane tip load. 
Varying the two trials are run on each case with varying time step refinement.
The time integration is performed using 2nd-order DIRK in order to verify its 
convergence behavior.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/test_beam.bdf")


f_mag = 10.0

# Integration order for DIRK solver
DIRK_order = 2

# KS function weight
ksweight = 100.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "load_coarse_ks_disp": 4.011613743975663,
        "load_coarse_mass": 0.27000000000000035,
        "load_fine_ks_disp": 5.113070671424432,
        "load_fine_mass": 0.27000000000000396,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """
        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-3
            self.dh = 1e-50
        else:
            self.rtol = 1e-1
            self.atol = 1e-4
            self.dh = 1e-5

        # Instantiate FEA Assembler
        fea_assembler = pytacs.pyTACS(bdf_file, comm)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # set transient problem options
        transientOptions = {"timeIntegrator": "DIRK", "integrationOrder": DIRK_order}

        # get some problem info
        n_vpn = fea_assembler.getVarsPerNode()

        # Create coarse load-specified transient problem
        coarse_prob = fea_assembler.createTransientProblem(
            name="load_coarse",
            tInit=0.0,
            tFinal=1.0,
            numSteps=8,
            options=transientOptions,
        )
        # Create fine load-specified transient problem
        fine_prob = fea_assembler.createTransientProblem(
            name="load_fine",
            tInit=0.0,
            tFinal=1.0,
            numSteps=32,
            options=transientOptions,
        )
        load_probs = [coarse_prob, fine_prob]

        for prob in load_probs:
            forces = np.zeros(n_vpn)
            ns = prob.getNumTimeSteps()
            for k in range(ns + 1):
                t_array = prob.getTimeStages(k)
                for s, t in enumerate(t_array):
                    f = f_mag * t**5
                    forces[2] = f  # applied to z-direction
                    prob.addLoadToNodes(
                        timeStep=k,
                        timeStage=s,
                        nodeIDs=21,
                        F=forces,
                        nastranOrdering=True,
                    )

        for problem in load_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                direction=[0.0, 0.0, 100.0],
                ftype="discrete",
            )

        return load_probs, fea_assembler
