import os
import unittest

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions, TACS

"""
6 noded beam model 1 meter long in x direction.
The cross-section is a solid rectangle (defined by PBARL) with the 
following properties:
    w = 0.1
    t = 0.05
We apply two load cases: a distributed gravity and distributed traction case.
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia, 
and Compliance functions and sensitivities.
"""

TACS_IS_COMPLEX = TACS.dtype == complex

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model_pbarl_rect.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "gravity_compliance": 1941.637970569529,
        "gravity_ks_vmfailure": 215.18186143357514,
        "gravity_mass": 13.500000000000004,
        "gravity_x_disp": -0.3750404568908798,
        "gravity_y_disp": 656.1378755065329,
        "gravity_z_disp": 275.45079293281793,
        "gravity_I_xx": 0.0,
        "gravity_I_xy": 0.0,
        "gravity_I_xz": 0.0,
        "gravity_I_yy": 1.13625,
        "gravity_I_yz": 0.0,
        "gravity_I_zz": 1.1278125,
        "traction_compliance": 4.325878285719202,
        "traction_ks_vmfailure": 9.71446798485217,
        "traction_mass": 13.500000000000004,
        "traction_x_disp": 0.009518432861763253,
        "traction_y_disp": -0.7122094288145717,
        "traction_z_disp": 12.02223266609341,
        "traction_I_xx": 0.0,
        "traction_I_xy": 0.0,
        "traction_I_xz": 0.0,
        "traction_I_yy": 1.13625,
        "traction_I_yz": 0.0,
        "traction_I_zz": 1.1278125,
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

        grav_prob = fea_assembler.createStaticProblem("gravity")
        grav_prob.addInertialLoad([-10.0, 3.0, 5.0])

        trac_prob = fea_assembler.createStaticProblem("traction")
        trac_prob.addTractionToComponents([0], [1.0, -2.0, 3.0])

        tacs_probs = [grav_prob, trac_prob]

        # Set convergence to be tight for test
        for problem in tacs_probs:
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
                "I_xx",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[1.0, 0.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xy",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_xz",
                functions.MomentOfInertia,
                direction1=[1.0, 0.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yy",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 1.0, 0.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_yz",
                functions.MomentOfInertia,
                direction1=[0.0, 1.0, 0.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
            )
            problem.addFunction(
                "I_zz",
                functions.MomentOfInertia,
                direction1=[0.0, 0.0, 1.0],
                direction2=[0.0, 0.0, 1.0],
                aboutCM=True,
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
