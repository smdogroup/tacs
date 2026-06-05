import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, functions

"""
Tests an unbalanced laminate shell model with the following layup: [0, 45, 30]s.
The laminate information is read in from a PCOMP card in the BDF file.
Two load cases are tested: an in-plane tension and out-of-plane shear.
tests KSDisplacement, KSFailure, StructuralMass, and Compliance functions
and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/comp_plate.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "Tension_compliance": 15637.920804938394,
        "Tension_ks_TsaiWufailure": 8.367136933085968,
        "Tension_mass": 1.1624999999999999,
        "Tension_x_disp": 0.08821099737382809,
        "Tension_y_disp": 0.027801038142742256,
        "Tension_z_disp": -1.926360906255938e-16,
        "VertShear_compliance": 0.00020961674296104168,
        "VertShear_ks_TsaiWufailure": 0.0013092603830388541,
        "VertShear_mass": 1.1624999999999999,
        "VertShear_x_disp": 7.95843722870956e-25,
        "VertShear_y_disp": 9.870262529474192e-23,
        "VertShear_z_disp": 0.010668323131955818,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-6
            self.atol = 1e-6
            self.dh = 1e-50
        else:
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        # Set convergence to be tight for test
        for problem in tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "ks_TsaiWufailure", functions.KSFailure, ksWeight=ksweight
            )
            problem.addFunction(
                "x_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[10.0, 0.0, 0.0],
            )
            problem.addFunction(
                "y_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 10.0, 0.0],
            )
            problem.addFunction(
                "z_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[0.0, 0.0, 10.0],
            )

        return tacs_probs, fea_assembler
