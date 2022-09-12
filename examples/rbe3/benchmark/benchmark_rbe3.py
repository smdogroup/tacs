"""
This script is used to regression test the example against historical values.
"""

import numpy as np
import unittest
import sys
import os

# Set the path to the example script we're testing
example_path = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(example_path)

# Reference values for eval functions
FUNC_REF = {
    "Axial_ks_disp": 0.45196571868787827,
    "Axial_mass": 211.71704428436388,
    "Moment-Bending_ks_disp": 1.0494950501687743,
    "Moment-Bending_mass": 211.71704428436388,
    "Shear-Bending_ks_disp": 0.3461665014945353,
    "Shear-Bending_mass": 211.71704428436388,
    "Torsion_ks_disp": 0.300370460972629,
    "Torsion_mass": 211.71704428436388,
}


class ExampleBenchmark(unittest.TestCase):

    N_PROCS = 8  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        # Import the example to automatically run the script
        import analysis

        self.example = analysis

    def benchmark_funcs(self):
        """
        Test the example eval functions against reference values
        """
        func_dict = self.example.funcs

        # Test functions values against historical values
        for func_name in func_dict:
            with self.subTest(function=func_name):
                np.testing.assert_allclose(
                    func_dict[func_name], FUNC_REF[func_name], rtol=1e-6, atol=1e-6
                )
