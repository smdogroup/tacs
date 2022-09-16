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
FUNC_REF = [1.070354064387641]


class ExampleBenchmark(unittest.TestCase):

    N_PROCS = 8  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        # Import the example to automatically run the script
        import crm

        self.example = crm

    def benchmark_funcs(self):
        """
        Test the example eval functions against reference values
        """
        fvals = self.example.fvals1

        # Test functions values against historical values
        np.testing.assert_allclose(fvals, FUNC_REF, rtol=1e-6, atol=1e-6)
