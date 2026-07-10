"""Tests for the constitutive derived-output API (getDerivedOutputs)."""

import unittest

import numpy as np

from tacs import TACS, constitutive

PANEL_LENGTH = 0.5
STIFFENER_PITCH = 0.2
PANEL_THICK = 0.02
STIFFENER_HEIGHT = 0.05
STIFFENER_THICK = 0.008
FLANGE_FRACTION = 0.8
PANEL_WIDTH = 0.3


def makePly():
    prop = constitutive.MaterialProperties(
        rho=1550.0,
        E1=54e9,
        E2=18e9,
        nu12=0.25,
        G12=9e9,
        G13=9e9,
        G23=9e9,
        Xt=2410.0e6,
        Xc=1040.0e6,
        Yt=73.0e6,
        Yc=173.0e6,
        S12=71.0e6,
    )
    return constitutive.OrthotropicPly(1e-3, prop)


def makeLayup():
    panelPlyAngles = np.deg2rad([0.0, 45.0, 90.0]).astype(TACS.dtype)
    panelPlyFracs = np.array([0.5, 0.3, 0.2], dtype=TACS.dtype)
    stiffenerPlyAngles = np.deg2rad([0.0, 60.0]).astype(TACS.dtype)
    stiffenerPlyFracs = np.array([0.7, 0.3], dtype=TACS.dtype)
    return panelPlyAngles, panelPlyFracs, stiffenerPlyAngles, stiffenerPlyFracs


def makeBladeCon():
    ply = makePly()
    panelPlyAngles, panelPlyFracs, stiffenerPlyAngles, stiffenerPlyFracs = makeLayup()
    return constitutive.BladeStiffenedShellConstitutive(
        ply,
        ply,
        PANEL_LENGTH,
        STIFFENER_PITCH,
        PANEL_THICK,
        panelPlyAngles,
        panelPlyFracs,
        STIFFENER_HEIGHT,
        STIFFENER_THICK,
        stiffenerPlyAngles,
        stiffenerPlyFracs,
        flangeFraction=FLANGE_FRACTION,
    )


def makeGPBladeCon():
    ply = makePly()
    panelPlyAngles, panelPlyFracs, stiffenerPlyAngles, stiffenerPlyFracs = makeLayup()
    return constitutive.GPBladeStiffenedShellConstitutive(
        ply,
        ply,
        PANEL_LENGTH,
        STIFFENER_PITCH,
        PANEL_THICK,
        panelPlyAngles,
        panelPlyFracs,
        STIFFENER_HEIGHT,
        STIFFENER_THICK,
        stiffenerPlyAngles,
        stiffenerPlyFracs,
        PANEL_WIDTH,
    )


class DerivedOutputTest(unittest.TestCase):
    def test_base_class_has_no_derived_outputs(self):
        """Assert that an IsoShellConstitutive, which declares no derived outputs, returns an empty dict from getDerivedOutputs."""
        prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324e6)
        con = constitutive.IsoShellConstitutive(prop, t=0.01)
        self.assertEqual(con.getDerivedOutputs(), {})

    def test_blade_derived_outputs(self):
        """
        Given a BladeStiffenedShellConstitutive,
        when getDerivedOutputs is called,
        then it returns effectiveThickness and effectiveBendingThickness, with
        effectiveThickness matching the analytic panel-plus-smeared-stiffener
        formula and effectiveBendingThickness exceeding the bare panel thickness.
        """
        con = makeBladeCon()
        outputs = con.getDerivedOutputs()
        self.assertEqual(
            list(outputs.keys()), ["effectiveThickness", "effectiveBendingThickness"]
        )
        # effectiveThickness = panelThick + stiffenerArea / stiffenerPitch,
        # with stiffenerArea = (1 + flangeFraction) * height * thickness
        stiffenerArea = (1.0 + FLANGE_FRACTION) * STIFFENER_HEIGHT * STIFFENER_THICK
        expected = PANEL_THICK + stiffenerArea / STIFFENER_PITCH
        np.testing.assert_allclose(
            np.real(outputs["effectiveThickness"]), expected, rtol=1e-12
        )
        # The stiffened shell must be bending-stiffer than the bare panel
        self.assertGreater(np.real(outputs["effectiveBendingThickness"]), PANEL_THICK)

    def test_gp_blade_derived_outputs(self):
        """
        Given a GPBladeStiffenedShellConstitutive,
        when getDerivedOutputs is called,
        then it returns the full named set of geometric and stiffness-ratio
        outputs in order, with stiffenerAspectRatio matching the height/thickness
        ratio and every returned value finite.
        """
        con = makeGPBladeCon()
        outputs = con.getDerivedOutputs()
        self.assertEqual(
            list(outputs.keys()),
            [
                "effectiveThickness",
                "effectiveBendingThickness",
                "stiffenerAspectRatio",
                "stiffenerAreaRatio",
                "affineAspectRatio",
                "laminateIsotropy",
                "stiffenerStiffnessRatio",
                "transverseShearParameter",
            ],
        )
        np.testing.assert_allclose(
            np.real(outputs["stiffenerAspectRatio"]),
            STIFFENER_HEIGHT / STIFFENER_THICK,
            rtol=1e-12,
        )
        for name, value in outputs.items():
            self.assertTrue(np.isfinite(np.real(value)), msg=name)


if __name__ == "__main__":
    unittest.main()
