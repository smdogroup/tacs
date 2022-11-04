import numpy as np
import os
from tacs import pytacs, elements, constitutive, functions
from pytacs_analysis_base_test import PyTACSTestCase

"""
This model tests an I-beam made up of two different element types:
    - Quad shell elements are used to model the web
    - Axial rod elements are used to model the caps
The properties of the cross-section are given below:
    Acaps = 0.05
    tweb = 0.1
We apply apply various tip loads test KSDisplacement, KSFailure, 
StructuralMass, and Compliance functions and sensitivities.

NOTE: Design variables are read in from BDF file.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/I_beam.bdf")

ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "load_set_000_compliance": 100430.19372545928,
        "load_set_000_ks_vmfailure": 0.4486594909084844,
        "load_set_000_mass": 5400.0,
        "load_set_000_x_disp": 0.3558707399683029,
        "load_set_000_y_disp": 1.2614949438778236,
        "load_set_000_z_disp": 0.34011973816621555,
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

        return tacs_probs, fea_assembler
