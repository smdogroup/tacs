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
        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        self.props = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            cte=cte,
            kappa=kappa,
        )
        t = 0.1
        d = 1.0
        dNum = 0
        tNum = 1

        # Create stiffness (need class)
        self.con = constitutive.IsoTubeBeamConstitutive(
            self.props, t=t, tNum=tNum, d=d, dNum=dNum
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
