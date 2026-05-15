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
    "StructAnalysis_compliance": np.float64(268273.6603263481),
    "StructAnalysis_l_skin_ksFailure": np.float64(1.4000486346076613),
    "StructAnalysis_l_skin_mass": np.float64(369.8464611079751),
    "StructAnalysis_mass": np.float64(1118.9051830205822),
    "StructAnalysis_rib_ksFailure": np.float64(0.48515871504376035),
    "StructAnalysis_rib_mass": np.float64(240.1806183164843),
    "StructAnalysis_spar_ksFailure": np.float64(1.4114538225383104),
    "StructAnalysis_spar_mass": np.float64(140.0232992825864),
    "StructAnalysis_u_skin_ksFailure": np.float64(2.27386728217599),
    "StructAnalysis_u_skin_mass": np.float64(368.85480431353795),
    "tipZDisp": np.float64(1.5973577336742282),
    "tipTwist": np.float64(3.940918119494375),
    "AdjCon_L_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "AdjCon_SPAR_panelThicknessAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerHeightAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerThicknessAdj": np.zeros(42),
    "AdjCon_U_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "DVCon_flangeThicknessMax": -0.0915 * np.ones(111),
    "DVCon_stiffenerAspectMax": 0.02 * np.ones(111),
    "DVCon_stiffenerAspectMin": -0.13 * np.ones(111),
    "PanelLengthCon_PanelLength": np.zeros(111),
}


class ExampleBenchmark(unittest.TestCase):
    N_PROCS = 8  # this is how many MPI processes to use for this TestCase.

    def setUp(self):
        # Import the example module and run analysis directly.
        import analysis

        self.funcs, _ = analysis.run_analysis()

    def benchmark_funcs(self):
        """
        Test the example eval functions against reference values
        """
        func_dict = self.funcs

        # Test functions values against historical values
        for func_name in FUNC_REF:
            with self.subTest(function=func_name):
                np.testing.assert_allclose(
                    func_dict[func_name], FUNC_REF[func_name], rtol=1e-6, atol=1e-6
                )


if __name__ == "__main__":
    unittest.main(defaultTest="ExampleBenchmark.benchmark_funcs")
