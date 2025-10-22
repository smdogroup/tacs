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
    "StructAnalysis_compliance": 181654.50493884285,
    "StructAnalysis_l_skin_ksFailure": 0.6725899627597012,
    "StructAnalysis_l_skin_mass": 845.363339675365,
    "StructAnalysis_mass": 2557.4975611899276,
    "StructAnalysis_rib_ksFailure": 0.25239263216123975,
    "StructAnalysis_rib_mass": 548.9842704376754,
    "StructAnalysis_spar_ksFailure": 0.6259609697539693,
    "StructAnalysis_spar_mass": 320.05325550306407,
    "StructAnalysis_u_skin_ksFailure": 1.9721825626648475,
    "StructAnalysis_u_skin_mass": 843.0966955738023,
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
    "DVCon_flangeThicknessMax": -0.294 * np.ones(111),
    "DVCon_stiffSpacingMin": 0.03 * np.ones(111),
    "DVCon_stiffenerAspectMax": -0.13 * np.ones(111),
    "DVCon_stiffenerAspectMin": 0.02 * np.ones(111),
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
            -2.35228248e-01,
            -2.35294003e-01,
            -2.35033747e-01,
            -2.35263347e-01,
            -2.34708238e-01,
            -2.35050844e-01,
            -3.27127159e-02,
            -3.28681957e-02,
            -3.27127235e-02,
            -3.28681997e-02,
            -3.27127338e-02,
            -3.28682062e-02,
            -3.27127439e-02,
            -3.28682122e-02,
            -3.27127548e-02,
            -3.28682186e-02,
            -3.27127656e-02,
            -3.28682245e-02,
            -3.27127797e-02,
            -3.28682332e-02,
            -3.27127939e-02,
            -3.28682416e-02,
            -3.27128095e-02,
            -3.28682508e-02,
            -3.27128259e-02,
            -3.28682600e-02,
            -3.27128462e-02,
            -3.28682724e-02,
            -3.27128679e-02,
            -3.28682852e-02,
            -3.27128902e-02,
            -3.28682975e-02,
            -3.27129253e-02,
            -3.28683213e-02,
            -3.27129465e-02,
            -3.28683294e-02,
            -3.27129940e-02,
            -3.28683618e-02,
            -3.27130264e-02,
            -3.28683766e-02,
            -3.27130883e-02,
            -3.28684175e-02,
            -3.27131394e-02,
            -3.28684433e-02,
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
