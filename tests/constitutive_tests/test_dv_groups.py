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


if __name__ == "__main__":
    unittest.main()
