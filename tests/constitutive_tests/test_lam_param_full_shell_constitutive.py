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
            self.rtol = 1e-3
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)

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
            alpha=cte,
            kappa=kappa,
        )
        iso_ply = constitutive.OrthotropicPly(thickness, iso_prop)

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
            alpha=cte,
            kappa=kappa,
        )
        ortho_ply = constitutive.OrthotropicPly(thickness, ortho_prop)

        self.ply_list = [iso_ply, ortho_ply]
        self.thickness = thickness
        self.thickness_dv_num = 0
        self.thinkness_lb = 0.0001
        self.thinkness_ub = 1.0
        self.lpNums = np.arange(0, 6, dtype=np.intc) + 1

        lpDVs = 0.6 * np.ones(6)
        self.dvs = np.concatenate(([self.thickness], lpDVs), dtype=self.dtype)

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def get_con(self, ply):
        con = constitutive.LamParamFullShellConstitutive(
            ply,
            self.thickness,
            self.thickness_dv_num,
            self.thinkness_lb,
            self.thinkness_ub,
            self.lpNums,
        )
        return con

    def test_constitutive_density(self):
        # Test density dv sensitivity
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
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
