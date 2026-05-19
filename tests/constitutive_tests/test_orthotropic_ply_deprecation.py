"""Tests for OrthotropicPly failure-criterion API: new enum kwarg, deprecated
boolean kwargs, and the validation errors raised for bad combinations."""

import unittest
import warnings

from tacs import constitutive


def _make_props():
    return constitutive.MaterialProperties(
        rho=1600.0, E1=150e9, E2=9e9, nu12=0.3, G12=7e9, G13=7e9, G23=5e9
    )


class TestOrthotropicPlyDeprecation(unittest.TestCase):
    def setUp(self):
        self.plyThickness = 0.005
        self.props = _make_props()

    # ------------------------------------------------------------------
    # DeprecationWarning cases (single deprecated bool kwarg)
    # ------------------------------------------------------------------

    def test_tsai_wu_criterion_warns(self):
        with self.assertWarns(DeprecationWarning) as cm:
            ply = constitutive.OrthotropicPly(
                self.plyThickness, self.props, tsai_wu_criterion=True
            )
        self.assertIsNotNone(ply)

    def test_max_strain_criterion_warns(self):
        with self.assertWarns(DeprecationWarning) as cm:
            ply = constitutive.OrthotropicPly(
                self.plyThickness, self.props, max_strain_criterion=True
            )
        self.assertIsNotNone(ply)

    # ------------------------------------------------------------------
    # ValueError: two deprecated bools at once
    # ------------------------------------------------------------------

    def test_two_deprecated_bools_raises(self):
        """Specifying two legacy bool kwargs must raise ValueError after the warning."""
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(ValueError, "Only one failure criterion"):
                constitutive.OrthotropicPly(
                    self.plyThickness,
                    self.props,
                    tsai_wu_criterion=True,
                    max_strain_criterion=True,
                )

    # ------------------------------------------------------------------
    # ValueError: mixing new enum kwarg with a deprecated bool kwarg
    # ------------------------------------------------------------------

    def test_mixed_new_and_deprecated_raises(self):
        """Using failure_criterion= alongside a legacy bool kwarg must raise ValueError."""
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(ValueError, "Cannot specify both"):
                constitutive.OrthotropicPly(
                    self.plyThickness,
                    self.props,
                    failure_criterion=constitutive.CompositeFailureCriterion.TSAI_WU,
                    tsai_wu_criterion=True,
                )

    # ------------------------------------------------------------------
    # New API: no warning emitted
    # ------------------------------------------------------------------

    def test_new_api_no_warning(self):
        """Using failure_criterion= enum kwarg must not emit any DeprecationWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ply = constitutive.OrthotropicPly(
                self.plyThickness,
                self.props,
                failure_criterion=constitutive.CompositeFailureCriterion.TSAI_WU,
            )
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        self.assertListEqual(
            deprecation_warnings,
            [],
            msg="No DeprecationWarning should be emitted when using the new enum kwarg.",
        )
        self.assertIsNotNone(ply)


if __name__ == "__main__":
    unittest.main()
