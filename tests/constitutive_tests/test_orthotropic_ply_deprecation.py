"""Tests for OrthotropicPly failure-criterion API: new enum kwarg, deprecated
boolean kwargs, and the validation errors raised for bad combinations.
Also covers the plyThickness/ply_thickness kwarg rename and the deprecated
setUseMaxStrainCriterion / setUseTsaiWuCriterion / setUseModifiedTsaiWuCriterion
setters.
"""

import unittest
import warnings

from tacs import constitutive
from tacs._testing import fails_at_version

_DEPRECATION_REMOVAL_VERSION = "3.14"


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

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_tsai_wu_criterion_warns(self):
        with self.assertWarns(DeprecationWarning):
            ply = constitutive.OrthotropicPly(
                self.plyThickness, self.props, tsai_wu_criterion=True
            )
        self.assertIsNotNone(ply)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_max_strain_criterion_warns(self):
        with self.assertWarns(DeprecationWarning):
            ply = constitutive.OrthotropicPly(
                self.plyThickness, self.props, max_strain_criterion=True
            )
        self.assertIsNotNone(ply)

    # ------------------------------------------------------------------
    # ValueError: two deprecated bools at once
    # ------------------------------------------------------------------

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
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

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_mixed_new_and_deprecated_raises(self):
        """Using failureCriterion= alongside a legacy bool kwarg must raise ValueError."""
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(ValueError, "Cannot specify both"):
                constitutive.OrthotropicPly(
                    self.plyThickness,
                    self.props,
                    failureCriterion=constitutive.CompositeFailureCriterion.TSAI_WU,
                    tsai_wu_criterion=True,
                )

    # ------------------------------------------------------------------
    # New API: no warning emitted
    # ------------------------------------------------------------------

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_new_api_no_warning(self):
        """Using failureCriterion= enum kwarg must not emit any DeprecationWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ply = constitutive.OrthotropicPly(
                self.plyThickness,
                self.props,
                failureCriterion=constitutive.CompositeFailureCriterion.TSAI_WU,
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


class TestOrthotropicPlyKwargRename(unittest.TestCase):
    """Tests for the plyThickness / ply_thickness kwarg rename."""

    def setUp(self):
        self.props = _make_props()

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_ply_thickness_kwarg_warns(self):
        """Passing ply_thickness= as a keyword must emit DeprecationWarning."""
        with self.assertWarns(DeprecationWarning):
            ply = constitutive.OrthotropicPly(ply_thickness=0.005, props=self.props)
        self.assertIsNotNone(ply)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_ply_thickness_and_ply_thickness_conflict_raises(self):
        """Passing both plyThickness= and ply_thickness= must raise ValueError."""
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(ValueError, "Cannot specify both"):
                constitutive.OrthotropicPly(
                    plyThickness=0.005,
                    props=self.props,
                    ply_thickness=0.005,
                )

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_ply_thickness_kwarg_produces_equivalent_ply(self):
        """A ply built with ply_thickness= kwarg must equal one built with plyThickness positionally."""
        with self.assertWarns(DeprecationWarning):
            ply_old = constitutive.OrthotropicPly(ply_thickness=0.005, props=self.props)
        ply_new = constitutive.OrthotropicPly(0.005, self.props)
        self.assertIs(
            type(ply_old.getMaterialProperties()),
            type(ply_new.getMaterialProperties()),
        )

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_failure_criterion_camel_case_no_warning(self):
        """Using failureCriterion= (camelCase) must not emit any DeprecationWarning."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ply = constitutive.OrthotropicPly(
                0.005,
                self.props,
                failureCriterion=constitutive.CompositeFailureCriterion.MAX_STRAIN,
            )
        deprecation_warnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        self.assertListEqual(deprecation_warnings, [])
        self.assertIsNotNone(ply)


class TestOrthotropicPlyDeprecatedSetters(unittest.TestCase):
    """Tests for the deprecated setUse*Criterion setters restored as shims."""

    def setUp(self):
        self.props = _make_props()
        self.ply = constitutive.OrthotropicPly(0.005, self.props)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_set_use_max_strain_criterion_warns(self):
        with self.assertWarns(DeprecationWarning):
            self.ply.setUseMaxStrainCriterion()

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_set_use_tsai_wu_criterion_warns(self):
        with self.assertWarns(DeprecationWarning):
            self.ply.setUseTsaiWuCriterion()

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_set_use_modified_tsai_wu_criterion_warns(self):
        with self.assertWarns(DeprecationWarning):
            self.ply.setUseModifiedTsaiWuCriterion()


if __name__ == "__main__":
    unittest.main()
