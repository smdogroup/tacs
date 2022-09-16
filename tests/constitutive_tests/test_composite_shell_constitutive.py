from tacs import TACS, constitutive, elements
import numpy as np
import unittest

DEG2RAD = np.pi / 180.0


class ConstitutiveTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-9
        else:
            self.dh = 1e-6
            self.rtol = 1e-1
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)
        # This constituitive class has no dvs
        self.dvs = np.array([], dtype=self.dtype)

        nplies = 3
        ply_thickness = 0.1

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
        iso_ply = constitutive.OrthotropicPly(ply_thickness, iso_prop)
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
        ortho_ply = constitutive.OrthotropicPly(ply_thickness, ortho_prop)
        ortho_layup = [ortho_ply] * nplies

        self.layup_list = [iso_layup, ortho_layup]
        self.ply_thicknesses = np.array([ply_thickness] * nplies, dtype=self.dtype)
        self.ply_angles = np.array([0.0, -45.0, 90.0], dtype=self.dtype) * DEG2RAD

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_constitutive_density(self):
        # Test density dv sensitivity
        for layup in self.layup_list:
            with self.subTest(layup=layup):
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
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
                con = constitutive.CompositeShellConstitutive(
                    layup, self.ply_thicknesses, self.ply_angles
                )
                fail = constitutive.TestConstitutiveFailureStrainSens(
                    con,
                    self.elem_index,
                    self.pt,
                    self.x,
                    self.dh,
                    self.print_level,
                    self.rtol,
                    self.atol,
                )
                self.assertFalse(fail)
