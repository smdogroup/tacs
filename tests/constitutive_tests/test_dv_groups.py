"""
Tests for the design variable group API on TACS constitutive objects.

Each constitutive class that defines design variables must report its DVs as named "groups"
(one per logical quantity) whose names and value types match the keyword arguments of the
class's Python constructor. For each class, these tests verify:

1. The group names, sizes, scalar/array types, values, and DV numbers match the values the
   object was constructed with.
2. The cross-API invariant: concatenating the active entries of every group, in group order,
   reproduces exactly the output of the pre-existing getDesignVars/getDesignVarNums methods.
"""

import unittest

import numpy as np

from tacs import TACS, constitutive


def makeIsoMaterial():
    return constitutive.MaterialProperties(
        rho=2700.0,
        specific_heat=921.096,
        E=70e3,
        nu=0.3,
        ys=270.0,
        alpha=24.0e-6,
        kappa=230.0,
    )


def makeOrthoMaterial():
    return constitutive.MaterialProperties(
        rho=1550.0,
        specific_heat=921.096,
        E1=54e3,
        E2=18e3,
        nu12=0.25,
        G12=9e3,
        G13=9e3,
        G23=9e3,
        Xt=2410.0,
        Xc=1040.0,
        Yt=73.0,
        Yc=173.0,
        S12=71.0,
        alpha=24.0e-6,
        kappa=230.0,
    )


class DVGroupTestCase(unittest.TestCase):
    """Base class providing the shared design-variable-group assertions."""

    dtype = TACS.dtype

    def assertGroupsConsistent(self, con, expectedValues, expectedNums):
        """Check the DV group API of ``con`` against expected values/nums and the legacy DV API.

        Parameters
        ----------
        con : tacs.TACS.Constitutive
            The constitutive object to check.
        expectedValues : dict[str, float or np.ndarray]
            Expected group values, in group order. A scalar entry indicates a scalar group.
        expectedNums : dict[str, int or np.ndarray]
            Expected DV numbers for each group, with the same scalar/array shapes as
            expectedValues.
        """
        groups = con.getDesignVarGroups()
        nums = con.getDesignVarGroupDVNums()

        # Names must match in order, since group order mirrors the legacy DV array order
        self.assertEqual(list(groups.keys()), list(expectedValues.keys()))
        self.assertEqual(list(nums.keys()), list(expectedValues.keys()))
        self.assertEqual(con.getNumDesignVarGroups(), len(expectedValues))

        for name, expVal in expectedValues.items():
            expNum = expectedNums[name]
            if np.ndim(expVal) == 0:
                # Scalar group: value and num must be plain scalars, not arrays
                self.assertNotIsInstance(groups[name], np.ndarray, msg=name)
                self.assertNotIsInstance(nums[name], np.ndarray, msg=name)
                np.testing.assert_allclose(
                    groups[name], expVal, rtol=1e-12, err_msg=name
                )
                self.assertEqual(nums[name], expNum, msg=name)
            else:
                self.assertIsInstance(groups[name], np.ndarray, msg=name)
                self.assertIsInstance(nums[name], np.ndarray, msg=name)
                np.testing.assert_allclose(
                    groups[name], expVal, rtol=1e-12, err_msg=name
                )
                np.testing.assert_array_equal(nums[name], expNum, err_msg=name)

        # Cross-API invariant: the active group entries concatenated in group order must
        # reproduce the legacy getDesignVars/getDesignVarNums output exactly
        flatNums = []
        flatVals = []
        for name in groups:
            groupNums = np.atleast_1d(nums[name])
            groupVals = np.atleast_1d(groups[name])
            for ii in range(len(groupNums)):
                if groupNums[ii] >= 0:
                    flatNums.append(groupNums[ii])
                    flatVals.append(groupVals[ii])
        np.testing.assert_array_equal(
            np.array(flatNums, dtype=np.intc), con.getDesignVarNums()
        )
        np.testing.assert_allclose(
            np.array(flatVals, dtype=self.dtype), con.getDesignVars(), rtol=1e-12
        )


class TestConstitutiveWithoutDVs(DVGroupTestCase):
    """A constitutive class with no DVs must report zero groups via the base-class defaults."""

    def test_no_dv_groups(self):
        con = constitutive.BasicBeamConstitutive(
            makeIsoMaterial(), A=1.0, J=2.0, Iy=3.0, Iz=4.0
        )
        self.assertEqual(con.getNumDesignVarGroups(), 0)
        self.assertEqual(con.getDesignVarGroups(), {})
        self.assertEqual(con.getDesignVarGroupDVNums(), {})
        self.assertGroupsConsistent(con, {}, {})


class TestIsoShellDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        con = constitutive.IsoShellConstitutive(makeIsoMaterial(), t=0.1, tNum=0)
        self.assertGroupsConsistent(con, {"t": 0.1}, {"t": 0})

    def test_dv_groups_inactive(self):
        con = constitutive.IsoShellConstitutive(makeIsoMaterial(), t=0.1, tNum=-1)
        self.assertGroupsConsistent(con, {"t": 0.1}, {"t": -1})


class TestPlaneStressDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        con = constitutive.PlaneStressConstitutive(makeIsoMaterial(), t=1.0, tNum=0)
        self.assertGroupsConsistent(con, {"t": 1.0}, {"t": 0})


class TestSolidDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        con = constitutive.SolidConstitutive(makeIsoMaterial(), t=1.0, tNum=0)
        self.assertGroupsConsistent(con, {"t": 1.0}, {"t": 0})


class TestPhaseChangeMaterialDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        solidProps = constitutive.MaterialProperties(
            rho=1.0, kappa=2.0, specific_heat=1.0
        )
        liquidProps = constitutive.MaterialProperties(
            rho=0.95, kappa=1.0, specific_heat=1.1
        )
        con = constitutive.PhaseChangeMaterialConstitutive(
            solidProps, liquidProps, lh=10.0, Tm=0.0, t=0.1, tNum=0
        )
        self.assertGroupsConsistent(con, {"t": 0.1}, {"t": 0})


class TestSmearedCompositeShellDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        ply = constitutive.OrthotropicPly(0.1, makeOrthoMaterial())
        plyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(self.dtype)
        plyFracs = np.array([0.5, 0.25, 0.25], dtype=self.dtype)
        plyFracNums = np.array([1, 2, 3], dtype=np.intc)
        con = constitutive.SmearedCompositeShellConstitutive(
            [ply] * 3, 0.1, plyAngles, plyFracs, 0, plyFracNums
        )
        self.assertGroupsConsistent(
            con,
            {"thickness": 0.1, "ply_fractions": plyFracs},
            {"thickness": 0, "ply_fractions": plyFracNums},
        )

    def test_dv_groups_partially_active(self):
        ply = constitutive.OrthotropicPly(0.1, makeOrthoMaterial())
        plyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(self.dtype)
        plyFracs = np.array([0.5, 0.25, 0.25], dtype=self.dtype)
        plyFracNums = np.array([0, -1, 1], dtype=np.intc)
        con = constitutive.SmearedCompositeShellConstitutive(
            [ply] * 3, 0.1, plyAngles, plyFracs, -1, plyFracNums
        )
        self.assertGroupsConsistent(
            con,
            {"thickness": 0.1, "ply_fractions": plyFracs},
            {"thickness": -1, "ply_fractions": plyFracNums},
        )


class TestLamParamFullShellDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        ply = constitutive.OrthotropicPly(0.1, makeOrthoMaterial())
        lpNums = np.array([1, 2, 3, 4, 5, 6], dtype=np.intc)
        con = constitutive.LamParamFullShellConstitutive(ply, 0.1, 0, 0.05, 0.2, lpNums)
        # The lamination parameter values cannot be set through the constructor;
        # they always start at 0 and are restored via setLaminationParameters
        self.assertGroupsConsistent(
            con,
            {"t": 0.1, "lp": np.zeros(6, dtype=self.dtype)},
            {"t": 0, "lp": lpNums},
        )


class TestLamParamSmearedShellDVGroups(DVGroupTestCase):
    def test_dv_groups(self):
        ply = constitutive.OrthotropicPly(0.1, makeOrthoMaterial())
        con = constitutive.LamParamSmearedShellConstitutive(
            ply,
            t=0.1,
            t_num=0,
            min_t=0.05,
            max_t=0.2,
            f0=0.25,
            f45=0.5,
            f90=0.25,
            f0_num=1,
            f45_num=2,
            f90_num=3,
            min_f0=0.1,
            min_f45=0.1,
            min_f90=0.1,
            W1=0.6,
            W3=0.6,
            W1_num=4,
            W3_num=5,
        )
        self.assertGroupsConsistent(
            con,
            {"t": 0.1, "f0": 0.25, "f45": 0.5, "f90": 0.25, "W1": 0.6, "W3": 0.6},
            {"t": 0, "f0": 1, "f45": 2, "f90": 3, "W1": 4, "W3": 5},
        )


if __name__ == "__main__":
    unittest.main()
