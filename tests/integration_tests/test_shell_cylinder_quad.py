import os

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, functions

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
        "Axial_compliance": 0.0046138506943161585,
        "Axial_mass": 208.58647120309152,
        "Shear-Bending_compliance": 0.3504119590966936,
        "Shear-Bending_mass": 208.58647120309152,
        "Moment-Bending_compliance": 0.03927297140294969,
        "Moment-Bending_mass": 208.58647120309152,
        "Torsion_compliance": 0.0527402966536922,
        "Torsion_mass": 208.58647120309152,
        "Combined_compliance": 0.447039077847727,
        "Combined_mass": 208.58647120309152,
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
            problem.setOption("useMonitor", True)

        return tacs_probs, fea_assembler
