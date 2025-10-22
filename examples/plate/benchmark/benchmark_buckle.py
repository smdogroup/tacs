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
    "buckle_eigsb.0": 20.082400947865022,
    "buckle_eigsb.1": 43.01440496545738,
    "buckle_eigsb.2": 80.31233170890538,
    "buckle_eigsb.3": 84.69751385025641,
    "buckle_eigsb.4": 86.82024838978569,
}


class ExampleBenchmark(unittest.TestCase):
    N_PROCS = 8  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        # Import the example to automatically run the script
        import plate_buckle

        self.example = plate_buckle

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
