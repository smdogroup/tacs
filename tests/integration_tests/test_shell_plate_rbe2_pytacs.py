import os
from tacs import pytacs, functions
from pytacs_analysis_base_test import PyTACSTestCase

r"""
Load a  structural model featuring two separate cantilevered plates connected by an RBE2 element.
Apply a load at the RBE2 center node and test KSFailure, StructuralMass, 
and Compliance functions and sensitivities. This is similar to test_shell_plate_rbe2.py, 
except it is run through the pytacs interface.
-----------        ----------- 
|          |\    /|          |
|          | \  / |          |
| Plate 1  |__\/__| Plate 2  |
|          |  /\  |          |
|          | /  \ |          |
|          |/    \|          |
------------       -----------
"""
base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/rbe_test.bdf")

# KS function weight
ksweight = 10.0


class ProblemTest(PyTACSTestCase.PyTACSTest):
    N_PROCS = 2  # this is how many MPI processes to use for this TestCase.

    FUNC_REFS = {
        "load_set_001_compliance": 424230.1227524751,
        "load_set_001_ks_disp": 0.5735091934240975,
        "load_set_001_ks_vmfailure": 0.6843335686556009,
        "load_set_001_mass": 108000.0,
    }

    def setup_tacs_problems(self, comm):
        """
        Setup pytacs object for problems we will be testing.
        """

        # Overwrite default check values
        if self.dtype == complex:
            self.rtol = 5e-6
            self.atol = 1e-8
            self.dh = 1e-9
        else:
            self.rtol = 2e-1
            self.atol = 1e-3
            self.dh = 1e-6

        # Instantiate FEA Assembler
        struct_options = {}

        fea_assembler = pytacs.pyTACS(bdf_file, comm, options=struct_options)

        # Set up constitutive objects and elements
        fea_assembler.initialize()

        # Create case 1 static problem from bdf file
        problems = fea_assembler.createTACSProbsFromBDF()

        for problem in problems.values():
            problem.addFunction("mass", functions.StructuralMass)
            problem.addFunction("ks_vmfailure", functions.KSFailure, ksWeight=ksweight)
            problem.addFunction("compliance", functions.Compliance)
            problem.addFunction(
                "ks_disp",
                functions.KSDisplacement,
                ksWeight=ksweight,
                direction=[1.0, 1.0, 1.0],
            )
            # Set convergence tol to be tight
            problem.setOption("L2Convergence", 1e-20)
            problem.setOption("L2ConvergenceRel", 1e-20)

        return problems.values(), fea_assembler
