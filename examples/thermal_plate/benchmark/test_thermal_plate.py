'''
This script is used to regression test the example against historical values.
'''

import numpy as np
import unittest
import sys
import os

# Set the path to the example script we're testing
example_path = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(example_path)

# Reference values for eval functions
FUNC_REF = {'SteadyState_avg_temp': 69.88016093991516,
            'SteadyState_ks_temp': 98.74014374789108,
            'SteadyState_mass': 39.20272476980967,
            'Transient_avg_temp': 397.96589286752857,
            'Transient_ks_temp': 97.83004117447352,
            'Transient_mass': 392.0272476980816}

class ExampleTest(unittest.TestCase):

    N_PROCS = 8  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        # Import the example to automatically run the script
        import analysis
        self.example = analysis

    def test_funcs(self):
        """
        Test the example eval functions against reference values
        """
        func_dict = self.example.funcs

        # Test functions values against historical values
        for func_name in func_dict:
            with self.subTest(function=func_name):
                np.testing.assert_allclose(func_dict[func_name], FUNC_REF[func_name], rtol=1e-6, atol=1e-6)