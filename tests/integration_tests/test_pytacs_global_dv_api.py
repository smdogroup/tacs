import os
import unittest

import numpy as np

from tacs import pytacs
from tacs.utilities import Error

try:
    import openmdao.api as om
    from tacs.mphys.dv import TacsDVComp
except ImportError:
    om = None

"""
Unit tests for the pyTACS addGlobalDV/assignMassDV API, covering both the
scalar and array-valued design variable paths (dv numbering, bounds/scale
broadcasting, error handling, and the mphys mass dv split).
No problems are solved.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/point_mass.bdf")


class GlobalDVAPITest(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.fea_assembler = pytacs.pyTACS(bdf_file)

    def test_scalar_dv(self):
        dvNum = self.fea_assembler.addGlobalDV("alpha", 1.0, lower=0.0, upper=2.0)
        self.assertIsInstance(dvNum, int)
        self.assertEqual(dvNum, 0)
        dv_dict = self.fea_assembler.getGlobalDVs()["alpha"]
        self.assertIsInstance(dv_dict["num"], int)
        self.assertEqual(self.fea_assembler.getGlobalDVNums(), [0])
        self.assertEqual(self.fea_assembler.getTotalNumGlobalDVs(), 1)
        self.assertEqual(len(self.fea_assembler.scaleList), 1)

    def test_array_dv(self):
        dvNums = self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        self.assertIsInstance(dvNums, np.ndarray)
        np.testing.assert_array_equal(dvNums, [0, 1, 2])
        self.assertEqual(self.fea_assembler.getGlobalDVNums(), [0, 1, 2])
        self.assertEqual(self.fea_assembler.getTotalNumGlobalDVs(), 3)
        self.assertEqual(len(self.fea_assembler.scaleList), 3)
        # Next dv should be numbered after the array
        nextNum = self.fea_assembler.addGlobalDV("alpha", 1.0)
        self.assertEqual(nextNum, 3)

    def test_array_dv_returned_nums_are_a_copy(self):
        dvNums = self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        dvNums[0] = 99
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["num"], [0, 1, 2])

    def test_array_dv_scalar_bounds_broadcast(self):
        self.fea_assembler.addGlobalDV(
            "g", [0.0, 0.0, -9.81], lower=-20.0, upper=20.0, scale=0.1
        )
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["lowerBound"], [-20.0] * 3)
        np.testing.assert_array_equal(dv_dict["upperBound"], [20.0] * 3)
        self.assertEqual(self.fea_assembler.scaleList, [0.1] * 3)

    def test_array_dv_array_bounds(self):
        self.fea_assembler.addGlobalDV(
            "g",
            [0.0, 0.0, -9.81],
            lower=[-1.0, -2.0, -20.0],
            upper=[1.0, 2.0, 0.0],
            scale=[1.0, 1.0, 0.1],
        )
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        np.testing.assert_array_equal(dv_dict["lowerBound"], [-1.0, -2.0, -20.0])
        np.testing.assert_array_equal(dv_dict["upperBound"], [1.0, 2.0, 0.0])
        self.assertEqual(self.fea_assembler.scaleList, [1.0, 1.0, 0.1])

    def test_array_dv_unbounded(self):
        self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81])
        dv_dict = self.fea_assembler.getGlobalDVs()["g"]
        self.assertIsNone(dv_dict["lowerBound"])
        self.assertIsNone(dv_dict["upperBound"])

    def test_empty_array_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [])

    def test_multidim_array_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", np.zeros((2, 3)))

    def test_bound_length_mismatch_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81], lower=[-1.0, -2.0])

    def test_scale_length_mismatch_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("g", [0.0, 0.0, -9.81], scale=[1.0, 2.0])

    def test_array_bound_with_scalar_value_raises(self):
        with self.assertRaises(Error):
            self.fea_assembler.addGlobalDV("alpha", 1.0, lower=[0.0])

    def test_assign_mass_dv_array_entry(self):
        dvNums = self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=1)
        dv_dict = self.fea_assembler.getGlobalDVs()["point_dvs"]
        np.testing.assert_array_equal(dv_dict["isMassDV"], [False, True])
        # The mass element should be wired to the dv num/value of entry 1
        self.assertEqual(self.fea_assembler.massDVs[1]["mNum"], dvNums[1])
        self.assertEqual(self.fea_assembler.massDVs[1]["m"], 30.0)

    def test_assign_mass_dv_array_entry_bounds(self):
        self.fea_assembler.addGlobalDV(
            "point_dvs", [20.0, 30.0], lower=[1.0, 2.0], upper=[100.0, 200.0]
        )
        self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=0)
        self.assertEqual(self.fea_assembler.massDVs[1]["mlb"], 1.0)
        self.assertEqual(self.fea_assembler.massDVs[1]["mub"], 100.0)

    def test_assign_mass_dv_array_requires_index(self):
        self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_dvs", 1, "m")

    def test_assign_mass_dv_array_index_out_of_range(self):
        self.fea_assembler.addGlobalDV("point_dvs", [20.0, 30.0])
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_dvs", 1, "m", index=2)

    def test_assign_mass_dv_scalar_rejects_index(self):
        self.fea_assembler.addGlobalDV("point_mass", 20.0)
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_mass", 1, "m", index=0)

    def test_assign_mass_dv_missing_key_rejects_index(self):
        with self.assertRaises(Error):
            self.fea_assembler.assignMassDV("point_mass", 1, "m", index=0)

    def test_assign_mass_dv_scalar_still_works(self):
        self.fea_assembler.addGlobalDV("point_mass", 20.0)
        self.fea_assembler.assignMassDV("point_mass", 1, "m")
        dv_dict = self.fea_assembler.getGlobalDVs()["point_mass"]
        self.assertTrue(dv_dict["isMassDV"])


@unittest.skipIf(om is None, "openmdao is not installed")
class TacsDVCompArrayMassTest(unittest.TestCase):
    N_PROCS = 1

    def test_array_mass_dv_split(self):
        fea_assembler = pytacs.pyTACS(bdf_file)
        fea_assembler.addGlobalDV("point_dvs", [20.0, 9.81])
        fea_assembler.assignMassDV("point_dvs", 1, "m", index=0)
        fea_assembler.initialize()

        vals = np.real(fea_assembler.getOrigDesignVars()).astype(float)
        prob = om.Problem()
        prob.model.add_subsystem(
            "dv_comp",
            TacsDVComp(
                fea_assembler=fea_assembler,
                initial_dv_vals=vals,
                separate_mass_dvs=True,
            ),
        )
        prob.setup()
        # Only the mass-assigned entry should be split out of the struct vector
        np.testing.assert_array_equal(prob.get_val("dv_comp.dv_mass_point_dvs"), [20.0])
        np.testing.assert_array_equal(prob.get_val("dv_comp.dv_struct"), [9.81])
        # Full dv vector should be reassembled in the original tacs ordering
        prob.set_val("dv_comp.dv_mass_point_dvs", 25.0)
        prob.run_model()
        np.testing.assert_array_equal(prob.get_val("dv_comp.tacs_dvs"), [25.0, 9.81])


if __name__ == "__main__":
    unittest.main()
