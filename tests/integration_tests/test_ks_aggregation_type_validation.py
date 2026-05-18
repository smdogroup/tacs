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

    def test_ks_temperature_constructor_rejects_discrete_average(self):
        with self.assertRaises(ValueError):
            functions.KSTemperature(
                self.assembler,
                ksAggregationType=functions.KSAggregationType.KS_DISCRETE_AVERAGE,
            )

    def test_ks_temperature_setter_rejects_discrete_average(self):
        func = functions.KSTemperature(self.assembler)
        with self.assertRaises(ValueError):
            func.setKSAggregationType(functions.KSAggregationType.KS_DISCRETE_AVERAGE)

    def test_ks_displacement_constructor_rejects_discrete_average(self):
        with self.assertRaises(ValueError):
            functions.KSDisplacement(
                self.assembler,
                ksAggregationType=functions.KSAggregationType.KS_DISCRETE_AVERAGE,
            )

    def test_ks_displacement_setter_rejects_discrete_average(self):
        func = functions.KSDisplacement(self.assembler)
        with self.assertRaises(ValueError):
            func.setKSAggregationType(functions.KSAggregationType.KS_DISCRETE_AVERAGE)

    def test_ks_failure_accepts_discrete_average(self):
        """DISCRETE_AVERAGE must not be rejected for KSFailure."""
        func = functions.KSFailure(
            self.assembler,
            ksAggregationType=functions.KSAggregationType.KS_DISCRETE_AVERAGE,
        )
        func.setKSAggregationType(functions.KSAggregationType.KS_DISCRETE_AVERAGE)


if __name__ == "__main__":
    unittest.main()
