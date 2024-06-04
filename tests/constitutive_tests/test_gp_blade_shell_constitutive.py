import unittest

import numpy as np

from tacs import TACS, constitutive, elements

DEG2RAD = np.pi / 180.0
np.random.seed(1342342)


class GPConstitutiveMLTest(unittest.TestCase):
    _my_debug = False

    def setUp(self):

        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-9
        else:
            self.dh = 1e-8
            self.rtol = 1e-3
        self.dtype = TACS.dtype

        # The failure value returned by the model is an aggregate of multiple
        # possible failure modes, we will run the failure tests multiple
        # times to try and cover strain states that cause each failure mode
        self.numFailureTests = 10

        # Basically, only check relative tolerance
        self.atol = self.rtol
        self.print_level = 2 if self._my_debug else 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)

        self.panelLength = 2.0
        self.panelLengthNum = 0
        self.stiffenerPitch = 0.2
        self.stiffenerPitchNum = 1
        self.stiffenerHeight = 0.075
        self.stiffenerHeightNum = 2
        self.stiffenerThickness = 1e-2
        self.stiffenerThicknessNum = 3
        self.panelThickness = 1.5e-2
        self.panelThicknessNum = 4

        self.numPanelPlies = 3
        self.panelPlyAngles = np.array([0.0, 45.0, 90.0], dtype=self.dtype) * DEG2RAD
        self.panelPlyFracs = np.random.random(self.numPanelPlies).astype(self.dtype)
        self.panelPlyFracs /= np.sum(self.panelPlyFracs)  # Make sure ply Fracs sum to 1
        self.panelPlyFracNums = np.arange(5, 5 + self.numPanelPlies, dtype=np.intc)

        self.numStiffenerPlies = 2
        self.stiffenerPlyAngles = np.array([0.0, 60.0], dtype=self.dtype) * DEG2RAD
        self.stiffenerPlyFracs = np.random.random(self.numStiffenerPlies).astype(
            self.dtype
        )
        self.stiffenerPlyFracs /= np.sum(
            self.stiffenerPlyFracs
        )  # Make sure ply Fracs sum to 1
        self.stiffenerPlyFracNums = np.arange(
            5 + self.numPanelPlies,
            5 + self.numPanelPlies + self.numStiffenerPlies,
            dtype=np.intc,
        )

        self.panelWidth = 1.0
        self.panelWidthNum = 5 + self.numPanelPlies + self.numStiffenerPlies

        self.dvs = (
            [
                self.panelLength,
                self.stiffenerPitch,
                self.stiffenerHeight,
                self.stiffenerThickness,
                self.panelThickness,
            ]
            + list(self.panelPlyFracs)
            + list(self.stiffenerPlyFracs)
            + [self.panelWidth]
        )
        self.dvs = np.array(self.dvs, dtype=self.dtype)

        self.flangeFraction = 0.8
        self.kcorr = 5.0 / 6.0

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
        iso_ply = constitutive.OrthotropicPly(1e-3, iso_prop)

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
        ortho_ply = constitutive.OrthotropicPly(1e-3, ortho_prop)

        self.ply_list = [iso_ply, ortho_ply]
        # self.ply_list = [ortho_ply]

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

        # construct the optional ML models
        n_train = 4

        n_param = constitutive.AxialGP.n_param
        self.axialGP = constitutive.AxialGP(
            n_train,
            Xtrain=np.random.rand(n_param * n_train).astype(self.dtype),
            alpha=np.random.rand(n_train).astype(self.dtype),
        )

        n_param = constitutive.ShearGP.n_param
        self.shearGP = constitutive.ShearGP(
            n_train,
            Xtrain=np.random.rand(n_param * n_train).astype(self.dtype),
            alpha=np.random.rand(n_train).astype(self.dtype),
        )

        #n_param = constitutive.CripplingGP.n_param
        #self.cripplingGP = constitutive.CripplingGP(
        #    n_train,
        #    Xtrain=np.random.rand(n_param * n_train).astype(self.dtype),
        #    alpha=np.random.rand(n_train).astype(self.dtype),
        #)
        self.cripplingGP = None

        self.panelGPs = constitutive.PanelGPs(
            self.axialGP,
            self.shearGP,
            self.cripplingGP,
            saveData=False
        )

    def get_con(self, ply):
        con = constitutive.GPBladeStiffenedShellConstitutive(
            ply,
            ply,
            self.panelLength,
            self.stiffenerPitch,
            self.panelThickness,
            self.panelPlyAngles,
            self.panelPlyFracs,
            self.stiffenerHeight,
            self.stiffenerThickness,
            self.stiffenerPlyAngles,
            self.stiffenerPlyFracs,
            self.panelWidth,
            self.kcorr,
            self.flangeFraction,
            self.panelLengthNum,
            self.stiffenerPitchNum,
            self.panelThicknessNum,
            self.panelPlyFracNums,
            self.stiffenerHeightNum,
            self.stiffenerThicknessNum,
            self.stiffenerPlyFracNums,
            self.panelWidthNum,
            self.panelGPs
        )
        # Set the KS weight really low so that all failure modes make a
        # significant contribution to the failure function derivatives
        con.setKSWeight(0.1)
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

    # def test_constitutive_specific_heat(self):
    #     # Test specific heat dv sensitivity
    #     for ply in self.ply_list:
    #         with self.subTest(ply=ply):
    #             con = self.get_con(ply)
    #             fail = constitutive.TestConstitutiveSpecificHeat(
    #                 con,
    #                 self.elem_index,
    #                 self.pt,
    #                 self.x,
    #                 self.dvs,
    #                 self.dh,
    #                 self.print_level,
    #                 self.atol,
    #                 self.rtol,
    #             )
    #             self.assertFalse(fail)

    # def test_constitutive_heat_flux(self):
    #     # Test heat flux dv sensitivity
    #     for ply in self.ply_list:
    #         with self.subTest(ply=ply):
    #             con = self.get_con(ply)
    #             fail = constitutive.TestConstitutiveHeatFlux(
    #                 con,
    #                 self.elem_index,
    #                 self.pt,
    #                 self.x,
    #                 self.dvs,
    #                 self.dh,
    #                 self.print_level,
    #                 self.atol,
    #                 self.rtol,
    #             )
    #             self.assertFalse(fail)

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

    # def test_constitutive_thermal_strain(self):
    #     # Test thermal strain dv sensitivity
    #     for ply in self.ply_list:
    #         with self.subTest(ply=ply):
    #             con = self.get_con(ply)
    #             fail = constitutive.TestConstitutiveThermalStrain(
    #                 con,
    #                 self.elem_index,
    #                 self.pt,
    #                 self.x,
    #                 self.dvs,
    #                 self.dh,
    #                 self.print_level,
    #                 self.atol,
    #                 self.rtol,
    #             )
    #             self.assertFalse(fail)

    def test_constitutive_failure(self):
        # Test failure dv sensitivity
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
                fail = False
                for _ in range(self.numFailureTests):
                    fail = fail or constitutive.TestConstitutiveFailure(
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
                fail = False
                for _ in range(self.numFailureTests):
                    fail = fail or constitutive.TestConstitutiveFailureStrainSens(
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

    def test_constitutive_internal_tests(self):
        """test all the internal or intermediate tests in C++"""
        for ply in self.ply_list:
            with self.subTest(ply=ply):
                con = self.get_con(ply)
                if TACS.dtype is complex:
                    dh = self.dh * 1j
                else:
                    dh = self.dh
                relError = con.test_all_derivative_tests(
                    dh, self.print_level
                )  # self.print_level
                fail = abs(relError.real) > self.rtol
                self.assertFalse(fail)
        if self.print_level != 0:
            print("All internal tests pass!")


class GPConstitutiveCFTest(GPConstitutiveMLTest):
    """
    second test with closed-form only, not ML models
    """

    def setUp(self):
        super(GPConstitutiveCFTest, self).setUp()

        # change the ML models to None so the constitutive class
        # uses closed-form solutions for buckling now.
        self.axialGP = None
        self.shearGP = None
        self.cripplingGP = None
        self.panelGPs = constitutive.PanelGPs(
            self.axialGP,
            self.shearGP,
            self.cripplingGP,
        )


if __name__ == "__main__":
    # current status of these tests
    # case == "CF" [runs internal tests] => all pass
    # case == "ML" [runs internal tests] => all pass
    # case == "all" [runs internal + external tests] => failure criterion fails.

    import argparse

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--case", type=str, default="full")
    args = parent_parser.parse_args()

    assert args.case in ["CF", "ML", "failStrain", "failDV", "full"]

    if args.case == "CF":
        tester = GPConstitutiveCFTest()
        tester._my_debug = True
        tester.setUp()
        tester.test_constitutive_internal_tests()
    elif args.case == "ML":
        tester = GPConstitutiveMLTest()
        tester._my_debug = True
        tester.setUp()
        tester.test_constitutive_internal_tests()
    elif args.case == "failStrain":
        # shouldn't matter which one of these I test as long as internal tests pass
        #tester = GPConstitutiveCFTest()
        tester = GPConstitutiveMLTest()
        tester._my_debug = True
        tester.setUp()
        tester.test_constitutive_failure_strain_sens()
    elif args.case == "failDV":
        # shouldn't matter which one of these I test as long as internal tests pass
        #tester = GPConstitutiveCFTest()
        tester = GPConstitutiveMLTest()
        tester._my_debug = True
        tester.setUp()
        tester.test_constitutive_failure()
    elif args.case == "full":
        unittest.main()
