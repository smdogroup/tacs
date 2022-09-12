"""
This script is used to regression test the example against historical values.
"""

import numpy as np
import unittest
import sys
import os

# Import the example to automatically run the script
import beam_opt

# Set the path to the example script we're testing
example_path = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(example_path)

# Reference value for optimal objective
OBJ_REF = 1.5536810027727979


class ExampleBenchmark(unittest.TestCase):

    N_PROCS = 3  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        self.example = beam_opt

    def benchmark_opt(self):
        """
        Test the example optimized objective against reference values
        """
        obj = self.example.m_opt

        # Test objective values against historical values
        np.testing.assert_allclose(obj, OBJ_REF, rtol=1e-3, atol=1e-3)
