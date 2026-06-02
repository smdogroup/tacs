import unittest

import numpy as np

from tacs import TACS, constitutive, elements


class ConstitutiveTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-12
        else:
            self.dh = 1e-8
            self.rtol = 1e-4
        self.dtype = TACS.dtype

        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)

        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)
        self.dvs = np.array([1.0, 0.1], dtype=self.dtype)

        # Create the isotropic material
        self.rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        self.props = constitutive.MaterialProperties(
            rho=self.rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            alpha=cte,
            kappa=kappa,
        )
        self.t = 0.1
        self.d = 1.0
        dNum = 0
        tNum = 1

        # Create stiffness (need class)
        self.con = constitutive.IsoTubeBeamConstitutive(
            self.props, t=self.t, tNum=tNum, d=self.d, dNum=dNum
        )

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_constitutive_density(self):
        # Test density dv sensitivity
        fail = constitutive.TestConstitutiveDensity(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_specific_heat(self):
        # Test specific heat dv sensitivity
        fail = constitutive.TestConstitutiveSpecificHeat(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_heat_flux(self):
        # Test heat flux dv sensitivity
        fail = constitutive.TestConstitutiveHeatFlux(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_stress(self):
        # Test stress dv sensitivity
        fail = constitutive.TestConstitutiveStress(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_thermal_strain(self):
        # Test thermal strain dv sensitivity
        fail = constitutive.TestConstitutiveThermalStrain(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_failure(self):
        # Test failure dv sensitivity
        fail = constitutive.TestConstitutiveFailure(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dvs,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_constitutive_failure_strain_sens(self):
        # Test failure dv sensitivity
        fail = constitutive.TestConstitutiveFailureStrainSens(
            self.con,
            self.elem_index,
            self.pt,
            self.x,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_cross_section_area_and_inertia(self):
        """Verify density and mass moments against the analytic equations.

        The outer diameter is d0 = inner_diameter + 2 * wall_thickness (diameter arithmetic).
        This test catches the bug where d0 was computed as inner + wall (off by one wall
        thickness), since the self-consistency DV-sensitivity tests cannot detect that error —
        the buggy value formula and its sensitivities are consistently wrong so they agree
        with each other even while both are incorrect.
        """
        # Analytic cross-section properties: d0 = d + 2*t (diameter arithmetic)
        d0 = self.d + 2.0 * self.t
        d1 = self.d
        A = np.pi * (d0**2 - d1**2) / 4.0
        Ia = np.pi * (d0**4 - d1**4) / 64.0

        density = self.con.evalDensity(self.elem_index, self.pt, self.x)
        np.testing.assert_allclose(density, self.rho * A, rtol=1e-10)

        moments = self.con.evalMassMoments(self.elem_index, self.pt, self.x)
        np.testing.assert_allclose(moments[0], self.rho * A, rtol=1e-10)
        np.testing.assert_allclose(moments[1], 0.0, atol=1e-14)
        np.testing.assert_allclose(moments[2], 0.0, atol=1e-14)
        np.testing.assert_allclose(moments[3], self.rho * Ia, rtol=1e-10)
        np.testing.assert_allclose(moments[4], self.rho * Ia, rtol=1e-10)
        np.testing.assert_allclose(moments[5], 0.0, atol=1e-14)
