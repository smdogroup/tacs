import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
Cylindrical beam constructed from shell elements. The beam is cantilevered at
one end and loaded at the other using an RBE3.

test StructuralMass and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/cylinder.bdf")


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "Axial_compliance": 0.004606292990644062,
        "Axial_mass": 208.58647120309135,
        "Shear-Bending_compliance": 0.32883253964394005,
        "Shear-Bending_mass": 208.58647120309135,
        "Moment-Bending_compliance": 0.03744304905762492,
        "Moment-Bending_mass": 208.58647120309135,
        "Torsion_compliance": 0.05100479047600847,
        "Torsion_mass": 208.58647120309135,
        "Combined_compliance": 0.4218866721681868,
        "Combined_mass": 208.58647120309135,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 1e-8
            self.atol = 1e-8
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

        # Add Functions
        for problem in tacs_probs:
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("compliance", functions.Compliance)

        return tacs_probs, fea_assembler
