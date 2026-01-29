import os
import unittest

from pytacs_analysis_base_test import PyTACSTestCase

from tacs import TACS, functions, pytacs

"""
9 noded beam model 1 meter long in x direction.
The cross-section is a hollow tube with the following properties:
    d = 0.1
    t = 0.01
We apply apply various tip loads test KSDisplacement, KSFailure,
StructuralMass, and Compliance functions and sensitivities.

This test uses automatic pytacs initialization from BDF file with PBARL cards.
"""

TACS_IS_COMPLEX = TACS.dtype == complex

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/offset_bar.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "load_set_000_cg_x": 0.49999999999999994,
        "load_set_000_cg_y": -0.05,
        "load_set_000_cg_z": 0.05,
        "load_set_000_compliance": 286783.9285716017,
        "load_set_000_ks_vmfailure": 15.833518185315949,
        "load_set_000_mass": 27.000000000000004,
        "load_set_000_x_disp": 0.34227078271035377,
        "load_set_000_y_disp": -0.20790980211605073,
        "load_set_000_z_disp": 2.4826507934089523,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        Uses automatic initialization from BDF file with PBARL cards.
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

        # Set up constitutive objects and elements automatically from BDF
        # This will read material properties and PBARL cards from the BDF file
        fea_assembler.initialize()

        # Read in forces from BDF and create tacs struct problems
        tacs_probs = fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        tacs_probs = tacs_probs.values()
        # Set convergence to be tight for test
        for problem in tacs_probs:
            print(problem.name)
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
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
            problem.addFunction(
                "cg_x", functions.CenterOfMass, direction=[1.0, 0.0, 0.0]
            )
            problem.addFunction(
                "cg_y", functions.CenterOfMass, direction=[0.0, 1.0, 0.0]
            )
            problem.addFunction(
                "cg_z", functions.CenterOfMass, direction=[0.0, 0.0, 1.0]
            )

        return tacs_probs, fea_assembler

    # We have to skip these tests in complex mode because the beam
    # element uses complex step to approximate the Jacobian and this
    # leads to issues with complex stepping the sensitivities.
    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_dv_sensitivities")
    def test_total_dv_sensitivities(self):
        super().test_total_dv_sensitivities()

    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_xpt_sensitivities")
    def test_total_xpt_sensitivities(self):
        super().test_total_xpt_sensitivities()
