import unittest

import numpy as np

from tacs import TACS, constitutive, elements


class ConstitutiveTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-9
        else:
            self.dh = 1e-8
            self.rtol = 1e-4
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)

        nplies = 3
        thickness = 0.1

        # Create the isotropic layup
        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        iso_prop = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E=E,
            nu=nu,
            ys=ys,
            cte=cte,
            kappa=kappa,
        )
        iso_ply = constitutive.OrthotropicPly(thickness, iso_prop)
        iso_layup = [iso_ply] * nplies

        # Create the orthotropic layup
        rho = 1550.0
        specific_heat = 921.096
        E1 = 54e3
        E2 = 18e3
        nu12 = 0.25
        G12 = 9e3
        G13 = 9e3
        Xt = 2410.0
        Xc = 1040.0
        Yt = 73.0
        Yc = 173.0
        S12 = 71.0
        cte = 24.0e-6
        kappa = 230.0
        ortho_prop = constitutive.MaterialProperties(
            rho=rho,
            specific_heat=specific_heat,
            E1=E1,
            E2=E2,
            nu12=nu12,
            G12=G12,
            G13=G13,
            G23=G13,
            Xt=Xt,
            Xc=Xc,
            Yt=Yt,
            Yc=Yc,
            S12=S12,
            cte=cte,
            kappa=kappa,
        )
        ortho_ply = constitutive.OrthotropicPly(thickness, ortho_prop)
        ortho_layup = [ortho_ply] * nplies

        self.layup_list = [iso_layup, ortho_layup]
        self.thickness = thickness
        self.thickness_dv_num = 0
        self.ply_angles = np.deg2rad([0.0, -45.0, 90.0]).astype(self.dtype)
        self.ply_fractions = np.array([0.333, 0.333, 0.333], dtype=self.dtype)
        self.ply_fraction_dv_nums = np.array([1, 2, 3], dtype=np.intc)
        self.dvs = np.array([0.1, 0.333, 0.333, 0.333], dtype=self.dtype)

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_constitutive_density(self):
        # Test density dv sensitivity
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveDensity(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveSpecificHeat(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveHeatFlux(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveStress(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveThermalStrain(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveFailure(
                    con,
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
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.SmearedCompositeShellConstitutive(
                    layup,
                    self.thickness,
                    self.ply_angles,
                    self.ply_fractions,
                    self.thickness_dv_num,
                    self.ply_fraction_dv_nums,
                )
                fail = constitutive.TestConstitutiveFailureStrainSens(
                    con,
                    self.elem_index,
                    self.pt,
                    self.x,
                    self.dh,
                    self.print_level,
                    self.atol,
                    self.rtol,
                )
                self.assertFalse(fail)
