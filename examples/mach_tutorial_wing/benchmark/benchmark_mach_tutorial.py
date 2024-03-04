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
    "StructAnalysis_compliance": 155444.88138037617,
    "StructAnalysis_l_skin_ksFailure": 0.6688691183266002,
    "StructAnalysis_l_skin_mass": 1056.704174594214,
    "StructAnalysis_mass": 3196.8719514874,
    "StructAnalysis_rib_ksFailure": 0.25350286472639316,
    "StructAnalysis_rib_mass": 686.2303380471093,
    "StructAnalysis_spar_ksFailure": 0.6336086764475523,
    "StructAnalysis_spar_mass": 400.06656937882445,
    "StructAnalysis_u_skin_ksFailure": 0.46502733574612126,
    "StructAnalysis_u_skin_mass": 1053.8708694672512,
    "AdjCon_L_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerPitchAdj": np.zeros(21),
    "AdjCon_L_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "AdjCon_SPAR_panelThicknessAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerHeightAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerPitchAdj": np.zeros(42),
    "AdjCon_SPAR_stiffenerThicknessAdj": np.zeros(42),
    "AdjCon_U_SKIN_panelThicknessAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerHeightAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerPitchAdj": np.zeros(21),
    "AdjCon_U_SKIN_stiffenerThicknessAdj": np.zeros(21),
    "DVCon_stiffSpacingMin": 0.08 * np.ones(111),
    "DVCon_stiffenerAspectMax": -0.15 * np.ones(111),
    "DVCon_stiffenerAspectMin": 0.01 * np.ones(111),
    "DVCon_thickDiffLimit": np.zeros(111),
    "PanelLengthCon_PanelLength": np.array(
        [
            2.06290983e-01,
            1.64727032e-01,
            1.11972828e-01,
            4.52157850e-02,
            3.01025244e-02,
            1.49892651e-02,
            -1.23995700e-04,
            -1.52372550e-02,
            -3.03505139e-02,
            -4.54637746e-02,
            -6.05770331e-02,
            -7.56902922e-02,
            -9.08035515e-02,
            -1.05916809e-01,
            -1.21030068e-01,
            -1.36143326e-01,
            -1.51256582e-01,
            -1.66369840e-01,
            -1.81483095e-01,
            -1.96596350e-01,
            -2.11709604e-01,
            -2.26822855e-01,
            -2.41936108e-01,
            2.26290983e-01,
            1.84727032e-01,
            1.31972828e-01,
            6.52157850e-02,
            5.01025244e-02,
            3.49892651e-02,
            1.98760043e-02,
            4.76274500e-03,
            -1.03505139e-02,
            -2.54637746e-02,
            -4.05770331e-02,
            -5.56902922e-02,
            -7.08035515e-02,
            -8.59168091e-02,
            -1.01030068e-01,
            -1.16143326e-01,
            -1.31256582e-01,
            -1.46369840e-01,
            -1.61483095e-01,
            -1.76596350e-01,
            -1.91709604e-01,
            -2.06822855e-01,
            -4.80049559e-02,
            -1.04612092e-02,
            2.52563166e-02,
            2.52563166e-02,
            1.15523772e-02,
            -2.15156100e-03,
            -1.58554993e-02,
            -2.95594381e-02,
            -4.32633758e-02,
            -5.69673146e-02,
            -7.06712519e-02,
            -8.43751890e-02,
            -9.80791272e-02,
            -1.11783064e-01,
            -1.25487000e-01,
            -1.39190937e-01,
            -1.52894872e-01,
            -1.66598808e-01,
            -1.80302742e-01,
            -1.94006675e-01,
            -2.07710608e-01,
            -2.21414538e-01,
            -1.16529874e-01,
            -1.16796551e-01,
            -1.16293732e-01,
            -1.16638202e-01,
            -1.15940931e-01,
            -1.16308447e-01,
            5.22872841e-02,
            5.21318043e-02,
            5.22872765e-02,
            5.21318003e-02,
            5.22872662e-02,
            5.21317938e-02,
            5.22872561e-02,
            5.21317878e-02,
            5.22872452e-02,
            5.21317814e-02,
            5.22872344e-02,
            5.21317755e-02,
            5.22872203e-02,
            5.21317668e-02,
            5.22872061e-02,
            5.21317584e-02,
            5.22871905e-02,
            5.21317492e-02,
            5.22871741e-02,
            5.21317400e-02,
            5.22871538e-02,
            5.21317276e-02,
            5.22871321e-02,
            5.21317148e-02,
            5.22871098e-02,
            5.21317025e-02,
            5.22870747e-02,
            5.21316787e-02,
            5.22870535e-02,
            5.21316706e-02,
            5.22870060e-02,
            5.21316382e-02,
            5.22869736e-02,
            5.21316234e-02,
            5.22869117e-02,
            5.21315825e-02,
            5.22868606e-02,
            5.21315567e-02,
        ]
    ),
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
