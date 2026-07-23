import os
import unittest
from unittest import mock

from mpi4py import MPI

from tacs import pytacs
from tacs.utilities import Error as TACSError

"""
Unit-style tests for the bounds handling of pyTACS.assignMassDV.

These check that lower/upper bounds passed to (or defaulted by) assignMassDV are
stored on the corresponding global DV entry, as reported by getGlobalDVs, and that
invalid bounds/property names are rejected. No FE analysis is run.
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/rigid_point_mass.bdf")


class AssignMassDVTest(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.fea_assembler = pytacs.pyTACS(bdf_file, MPI.COMM_WORLD)

    def test_diagonal_property_defaults_lower_bound_to_zero(self):
        """A diagonal moment of inertia defaults to a zero lower, unbounded upper."""
        self.fea_assembler.assignMassDV("moi", 1, "I11")
        dv = self.fea_assembler.getGlobalDVs()["moi"]
        self.assertEqual(dv["lowerBound"], 0.0)
        self.assertIsNone(dv["upperBound"])

    def test_product_property_defaults_to_unbounded(self):
        """A product of inertia may be either sign, so both bounds default to None."""
        self.fea_assembler.assignMassDV("poi", 1, "I12")
        dv = self.fea_assembler.getGlobalDVs()["poi"]
        self.assertIsNone(dv["lowerBound"])
        self.assertIsNone(dv["upperBound"])

    def test_explicit_bounds_are_stored(self):
        """Explicit bounds are stored verbatim, including a negative lower on a product."""
        self.fea_assembler.assignMassDV("m1", 1, "m", lower=1.0, upper=5.0)
        self.fea_assembler.assignMassDV("poi", 1, "I12", lower=-5.0, upper=5.0)
        globalDVs = self.fea_assembler.getGlobalDVs()
        self.assertEqual(globalDVs["m1"]["lowerBound"], 1.0)
        self.assertEqual(globalDVs["m1"]["upperBound"], 5.0)
        self.assertEqual(globalDVs["poi"]["lowerBound"], -5.0)
        self.assertEqual(globalDVs["poi"]["upperBound"], 5.0)

    def test_negative_lower_on_non_negative_property_raises(self):
        """Mass and diagonal moments of inertia cannot have a negative lower bound."""
        for dvName in ("m", "I11", "I22", "I33"):
            with self.assertRaises(TACSError):
                self.fea_assembler.assignMassDV(f"bad_{dvName}", 1, dvName, lower=-1.0)

    def test_unrecognized_property_name_raises(self):
        """An unknown mass property name is rejected early."""
        with self.assertRaises(TACSError):
            self.fea_assembler.assignMassDV("typo", 1, "I99")

    def test_lower_greater_than_upper_raises(self):
        """An inverted (empty) bound range is rejected when both bounds are concrete."""
        with self.assertRaises(TACSError):
            self.fea_assembler.assignMassDV("m1", 1, "m", lower=5.0, upper=1.0)

    def test_scale_stored_on_global_dv(self):
        """The scale factor is stored on the global DV entry, defaulting to 1.0."""
        self.fea_assembler.assignMassDV("scaled", 1, "m", scale=3.0)
        self.fea_assembler.assignMassDV("default", 2, "m")
        globalDVs = self.fea_assembler.getGlobalDVs()
        self.assertEqual(globalDVs["scaled"]["scale"], 3.0)
        self.assertEqual(globalDVs["default"]["scale"], 1.0)

    def test_scale_survives_initialize(self):
        """A mass DV scale must not be wiped by the default-callback initialize path."""
        self.fea_assembler.assignMassDV("scaled", 1, "m", scale=4.0)
        dvNum = self.fea_assembler.getGlobalDVs()["scaled"]["num"]
        self.fea_assembler.initialize()
        self.assertEqual(self.fea_assembler.scaleList[dvNum], 4.0)

    def test_unspecified_fields_preserve_prior_addGlobalDV_values(self):
        """Omitted lower/upper/scale must not clobber values set by a prior addGlobalDV."""
        self.fea_assembler.addGlobalDV("pre", 10.0, lower=2.0, upper=8.0, scale=2.0)
        self.fea_assembler.assignMassDV("pre", 1, "m")
        dv = self.fea_assembler.getGlobalDVs()["pre"]
        self.assertEqual(dv["lowerBound"], 2.0)
        self.assertEqual(dv["upperBound"], 8.0)
        self.assertEqual(dv["scale"], 2.0)

    def test_overriding_prior_bounds_warns_and_applies(self):
        """Overriding a bound/scale set by a prior addGlobalDV warns but still applies."""
        self.fea_assembler.addGlobalDV("pre", 10.0, lower=2.0, scale=2.0)
        with mock.patch.object(self.fea_assembler, "_TACSWarning") as mockWarn:
            self.fea_assembler.assignMassDV("pre", 1, "m", lower=0.0, scale=5.0)
        self.assertGreaterEqual(mockWarn.call_count, 1)
        dv = self.fea_assembler.getGlobalDVs()["pre"]
        self.assertEqual(dv["lowerBound"], 0.0)
        self.assertEqual(dv["scale"], 5.0)

    def test_setting_previously_unset_fields_does_not_warn(self):
        """Setting bounds on a DV that only had a value (no bounds) must not warn."""
        self.fea_assembler.addGlobalDV("pre", 10.0)
        with mock.patch.object(self.fea_assembler, "_TACSWarning") as mockWarn:
            self.fea_assembler.assignMassDV("pre", 1, "m", lower=0.0, upper=100.0)
        mockWarn.assert_not_called()


if __name__ == "__main__":
    unittest.main()
