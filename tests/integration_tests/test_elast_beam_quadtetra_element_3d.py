import numpy as np
import os
from tacs import pytacs, TACS, elements, constitutive, functions, problems
from pytacs_analysis_base_test import PyTACSTestCase

"""
10m x 1m x 1m Square beam constructed from tetrahedral elements. The beam is cantilevered at
one end and loaded using a pressure on the opposite face. This leads to a pure axial loading on the model.
A second case is run where the beam is hung under a uniform gravity load.

test StructuralMass and Compliance functions and sensitivities
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/solid_beam.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "AxialPressure_compliance": 14226.813534942703,
        "AxialPressure_ks_disp": 0.3088436763967383,
        "AxialPressure_mass": 27000.00000000002,
        "Gravity_compliance": 33000.844246090586,
        "Gravity_ks_disp": 0.3669703989454973,
        "Gravity_mass": 27000.00000000002,
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
            self.rtol = 5e-6
            self.atol = 1e-2
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
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[-100.0, -100.0, -100.0],
            )

        return tacs_probs, fea_assembler
