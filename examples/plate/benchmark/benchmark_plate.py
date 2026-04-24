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
    "dynamic_modes_eigsm.0": 83438.6990942435,
    "dynamic_modes_eigsm.1": 350167.23587092856,
    "dynamic_modes_eigsm.2": 350167.23587098264,
    "dynamic_modes_eigsm.3": 761399.6241353224,
    "dynamic_modes_eigsm.4": 1143887.0924026824,
    "dynamic_modes_eigsm.5": 1154770.1662673794,
    "dynamic_modes_eigsm.6": 1789335.2577894062,
    "dynamic_modes_eigsm.7": 1789335.2577927846,
    "dynamic_modes_eigsm.8": 2993279.554315713,
    "dynamic_modes_eigsm.9": 3190753.194385857,
    "point_force_ks_vmfailure": 1.4546105780138643,
    "point_force_mass": 12.499999999999943,
    "pressure_ks_vmfailure": 0.352202299866338,
    "pressure_mass": 125.00000000003863,
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
