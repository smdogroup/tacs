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
    "TipForce_Compliance": 2.7437177731573215,
    "TipForce_KSFailure": 0.05892727757722518,
    "TipForce_MaxZDisp": 2.9603338686052387,
    "TipMoment_Compliance": 3.701101650340912,
    "TipMoment_KSFailure": 0.05893060077376249,
    "TipMoment_MaxZDisp": 3.0515204191877943,
}


class ExampleBenchmark(unittest.TestCase):

    N_PROCS = 1  # this is how many MPI processes to use for this TestCase.

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
