"""
Tests that KSAggregationType.KS_DISCRETE_AVERAGE raises ValueError when passed
to KSTemperature or KSDisplacement (constructor and setKSAggregationType).
"""

import os
import unittest

import numpy as np
from mpi4py import MPI
from tacs import constitutive, elements, functions, pytacs

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "input_files/plate.bdf")


def _make_assembler():
    """Return a minimal TACS Assembler for use in function constructor tests."""

    def elem_callback(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
        props = constitutive.MaterialProperties(rho=2700.0, E=70e9, nu=0.3, ys=270e6)
        con = constitutive.IsoShellConstitutive(props, t=0.01, tNum=dvNum)
        transform = elements.ShellNaturalTransform()
        elem = elements.Quad4Shell(transform, con)
        return elem

    fea = pytacs.pyTACS(bdf_file, comm=MPI.COMM_WORLD)
    fea.initialize(elem_callback)
    return fea.assembler


class TestDiscreteAverageRejection(unittest.TestCase):
    """Verify DISCRETE_AVERAGE is rejected by KSTemperature and KSDisplacement."""

    @classmethod
    def setUpClass(cls):
        cls.assembler = _make_assembler()

    def test_constructor_rejects_discrete_average(self):
        for cls in (functions.KSTemperature, functions.KSDisplacement):
            with self.subTest(cls=cls.__name__):
                with self.assertRaises(ValueError):
                    cls(
                        self.assembler,
                        ksAggregationType=functions.KSAggregationType.KS_DISCRETE_AVERAGE,
                    )

    def test_setter_rejects_discrete_average(self):
        for cls in (functions.KSTemperature, functions.KSDisplacement):
            with self.subTest(cls=cls.__name__):
                func = cls(self.assembler)
                with self.assertRaises(ValueError):
                    func.setKSAggregationType(
                        functions.KSAggregationType.KS_DISCRETE_AVERAGE
                    )

    def test_ks_failure_accepts_discrete_average(self):
        """DISCRETE_AVERAGE must not be rejected for KSFailure."""
        func = functions.KSFailure(
            self.assembler,
            ksAggregationType=functions.KSAggregationType.KS_DISCRETE_AVERAGE,
        )
        func.setKSAggregationType(functions.KSAggregationType.KS_DISCRETE_AVERAGE)


class TestDeprecatedFtypeKwarg(unittest.TestCase):
    """Verify backward-compat behaviour of the deprecated ftype= string kwarg."""

    @classmethod
    def setUpClass(cls):
        cls.assembler = _make_assembler()

    # ------------------------------------------------------------------
    # DeprecationWarning cases
    # ------------------------------------------------------------------

    def test_ftype_warns(self):
        cases = [
            (functions.KSTemperature, "discrete"),
            (functions.KSFailure, "discrete-average"),
            (functions.KSDisplacement, "continuous"),
        ]
        for cls, ftype in cases:
            with self.subTest(cls=cls.__name__, ftype=ftype):
                with self.assertWarns(DeprecationWarning):
                    func = cls(self.assembler, ftype=ftype)
                self.assertIsNotNone(func)

    # ------------------------------------------------------------------
    # ValueError for mixed old + new kwargs
    # ------------------------------------------------------------------

    def test_mixed_kwargs_raises(self):
        for cls in (
            functions.KSTemperature,
            functions.KSFailure,
            functions.KSDisplacement,
        ):
            with self.subTest(cls=cls.__name__):
                with self.assertRaises(ValueError):
                    cls(
                        self.assembler,
                        ftype="continuous",
                        ksAggregationType=functions.KSAggregationType.KS_CONTINUOUS,
                    )

    # ------------------------------------------------------------------
    # ValueError for unrecognised or disallowed ftype strings
    # ------------------------------------------------------------------

    def test_unknown_ftype_raises_value_error(self):
        """An unrecognised ftype string must raise ValueError."""
        with self.assertRaises(ValueError):
            functions.KSTemperature(self.assembler, ftype="bogus")

    def test_ftype_discrete_average_raises_for_temperature(self):
        """'discrete-average' is not in the KSTemperature ftype map and must raise ValueError."""
        with self.assertRaises(ValueError):
            functions.KSTemperature(self.assembler, ftype="discrete-average")


if __name__ == "__main__":
    unittest.main()
