import unittest

import numpy as np

from tacs import TACS, constitutive, elements


class ConstitutiveTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-11
        else:
            self.dh = 1e-8
            self.rtol = 1e-6
        self.dtype = TACS.dtype

        self.atol = np.clip(1e-5 * self.rtol, 1e-8, 1e-14)
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)
        self.dvs = np.array([1.0], dtype=self.dtype)

        # Create the isotropic material
        rho = 1.0  # Density kg/m^3
        kappa = 1.0  # Thermal conductivity W/(m⋅K)
        cp = 1.0  # Specific heat J/(kg⋅K)
        lh = 10.0  # Latent heat J/kg
        mt = 0.0  # Melting temperature (relative) K

        solid_prop = constitutive.MaterialProperties(
            rho=rho, kappa=2.0 * kappa, specific_heat=cp
        )
        liquid_prop = constitutive.MaterialProperties(
            rho=0.95 * rho, kappa=kappa, specific_heat=1.1 * cp
        )

        # Set one thickness value for every component
        self.con = constitutive.PhaseChangeMaterialConstitutive(
            solid_prop, liquid_prop, lh=lh, mt=mt, t=0.1, tNum=-1
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
