import os
import unittest

import numpy as np

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, elements, constitutive, functions, TACS

"""
6 noded beam model 1 meter long in x direction.
The cross-section is an I beam (defined by PBARL) with the following properties:

    h_web = 0.1
    t_web = 0.05
    t_flange = 0.01
    w_flange = 0.05

We apply two load cases: a distributed gravity and distributed traction case.
We apply apply various tip loads test KSDisplacement, StructuralMass, MomentOfInertia, 
and Compliance functions and sensitivities.
"""

TACS_IS_COMPLEX = TACS.dtype == complex

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/beam_model_pbarl_Ibeam.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "gravity_I_xx": 0.0,
        "gravity_I_xy": 0.0,
        "gravity_I_xz": 0.0,
        "gravity_I_yy": 0.4055805000000001,
        "gravity_I_yz": 0.0,
        "gravity_I_zz": 0.4116532500000001,
        "gravity_compliance": 2.0369377594102223,
        "gravity_mass": 4.859999999999999,
        "gravity_x_disp": -0.0012712022309512598,
        "gravity_y_disp": 0.04906289294167329,
        "gravity_z_disp": 1.6803054202705177,
        "traction_I_xx": 0.0,
        "traction_I_xy": 0.0,
        "traction_I_xz": 0.0,
        "traction_I_yy": 0.4055805000000001,
        "traction_I_yz": 0.0,
        "traction_I_zz": 0.4116532500000001,
        "traction_compliance": 0.03118108006818735,
        "traction_mass": 4.859999999999999,
        "traction_x_disp": 2.619117597131193e-05,
        "traction_y_disp": -0.005861111298040609,
        "traction_z_disp": 0.13015006699801424,
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

    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_dv_sensitivities")
    def test_total_dv_sensitivities(self):
        super().test_total_dv_sensitivities()

    @unittest.skipIf(TACS_IS_COMPLEX, "Skipping test_total_xpt_sensitivities")
    def test_total_xpt_sensitivities(self):
        super().test_total_xpt_sensitivities()
