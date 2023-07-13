import os
import tempfile

from pytacs_analysis_base_test import PyTACSTestCase
from tacs import pytacs, functions

"""
This case tests pyTACS's `writeBDF` method for shell elements.
We first instantiate pyTACS and structural problems from a provided BDF, as usual.
Using that instance of pyTACS we then export a BDF file using the `writeBDF` method,
This should be identical to the original BDF file. We then instantiate a second 
instance of pyTACS from this secondary BDF file and continue the tests with the secondary model.
This test ensures that `writeBDF` generates a consistent file based on the pyTACS model.
The test results should be identical to those in test_beam_bend_coupling.py
 
We test KSDisplacement, StructuralMass, MomentOfInertia, and Compliance functions and sensitivities.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
orig_bdf_file = os.path.join(base_dir, "./input_files/beam_model.bdf")

ksweight = 10.0

from test_beam_bend_coupling import ProblemTest as PT


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = PT.FUNC_REFS

    def setup_tacs_problems(self, comm):
        """
        Setup mesh and pytacs object for problem we will be testing.
        """

        # Create a temporary file on the root proc to hold our intermediate bdf
        if comm.rank == 0:
            temp_file = tempfile.NamedTemporaryFile(suffix=".bdf")
            new_bdf_file = temp_file.name
        else:
            new_bdf_file = None

        # Broadcast the temp file name to other procs
        new_bdf_file = comm.bcast(new_bdf_file, root=0)

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
        orig_fea_assembler = pytacs.pyTACS(orig_bdf_file, comm)
        # Set up constitutive objects and elements
        orig_fea_assembler.initialize()
        # Read in forces from BDF and create tacs struct problems
        orig_tacs_probs = orig_fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        orig_tacs_probs = orig_tacs_probs.values()
        # Set convergence to be tight for test
        for problem in orig_tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)
        # We have to solve each problem before writing bdf file
        for problem in orig_tacs_probs:
            problem.solve()
        # Write out BDF file
        orig_fea_assembler.writeBDF(new_bdf_file, orig_tacs_probs)

        # Create a new pytacs instance from BDF file we wrote
        new_fea_assembler = pytacs.pyTACS(new_bdf_file, comm)
        # Set up constitutive objects and elements
        new_fea_assembler.initialize()
        # Read in forces from BDF and create tacs struct problems
        new_tacs_probs = new_fea_assembler.createTACSProbsFromBDF()
        # Convert from dict to list
        new_tacs_probs = new_tacs_probs.values()

        # Set convergence to be tight for test
        for problem in new_tacs_probs:
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        # Add Functions
        for problem in new_tacs_probs:
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

        return new_tacs_probs, new_fea_assembler
