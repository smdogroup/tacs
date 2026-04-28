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
    "AdjCon_L_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "AdjCon_SPAR_panelThicknessAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerHeightAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerThicknessAdj": np.zeros(42),
    "AdjCon_U_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "DVCon_flangeThicknessMax": -0.091 * np.ones(111),
    "DVCon_stiffenerAspectMax": 0.02 * np.ones(111),
    "DVCon_stiffenerAspectMin": -0.13 * np.ones(111),
    "PanelLengthCon_PanelLength": np.zeros(111),
    "StructAnalysis_compliance": np.float64(267676.5038258135),
    "StructAnalysis_l_skin_ksFailure": np.float64(1.3711772904083315),
    "StructAnalysis_l_skin_mass": np.float64(369.8464611079735),
    "StructAnalysis_mass": np.float64(1118.9051830205872),
    "StructAnalysis_rib_ksFailure": np.float64(0.5057170301961897),
    "StructAnalysis_rib_mass": np.float64(240.18061831648674),
    "StructAnalysis_spar_ksFailure": np.float64(1.3748571539948513),
    "StructAnalysis_spar_mass": np.float64(140.02329928258945),
    "StructAnalysis_u_skin_ksFailure": np.float64(2.2588423985957515),
    "StructAnalysis_u_skin_mass": np.float64(368.8548043135387),
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
        for func_name in FUNC_REF:
            with self.subTest(function=func_name):
                np.testing.assert_allclose(
                    func_dict[func_name], FUNC_REF[func_name], rtol=1e-6, atol=1e-6
                )
