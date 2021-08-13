from tacs import TACS, constitutive
import numpy as np
import unittest


class ConstitutiveTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
        else:
            self.dh = 1e-6
        self.dtype = TACS.dtype

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.rtol = 5e-5
        self.print_level = 0

        # Set element index
        self.elem_index = 0

        # Set the variable arrays
        self.x = np.ones(3, dtype=self.dtype)
        self.pt = np.zeros(3)
        self.dvs = np.array([1.0], dtype=self.dtype)

        # Create the isotropic material
        rho = 2700.0
        specific_heat = 921.096
        E = 70e3
        nu = 0.3
        ys = 270.0
        cte = 24.0e-6
        kappa = 230.0
        self.props = constitutive.MaterialProperties(rho=rho, specific_heat=specific_heat,
                                                     E=E, nu=nu, ys=ys, cte=cte, kappa=kappa)

        # Create stiffness (need class)
        self.con = constitutive.PlaneStressConstitutive(self.props, t=1.0, tNum=0)

    def testConstitutiveDensity(self):
        # Test density dv sensitivity
        fail = constitutive.TestConstitutiveDensity(self.con, self.elem_index, self.pt, self.x,
                                                    self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveSpecificHeat(self):
        # Test specific heat dv sensitivity
        fail = constitutive.TestConstitutiveSpecificHeat(self.con, self.elem_index, self.pt, self.x,
                                                         self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveHeatFlux(self):
        # Test heat flux dv sensitivity
        fail = constitutive.TestConstitutiveHeatFlux(self.con, self.elem_index, self.pt, self.x,
                                                     self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveStress(self):
        # Test stress dv sensitivity
        fail = constitutive.TestConstitutiveStress(self.con, self.elem_index, self.pt, self.x,
                                                   self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveThermalStrain(self):
        # Test thermal strain dv sensitivity
        fail = constitutive.TestConstitutiveThermalStrain(self.con, self.elem_index, self.pt, self.x,
                                                          self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveFailure(self):
        # Test failure dv sensitivity
        fail = constitutive.TestConstitutiveFailure(self.con, self.elem_index, self.pt, self.x,
                                                    self.dvs, self.dh, self.print_level, self.atol, self.rtol)
        self.assertFalse(fail)

    def testConstitutiveFailureStrainSens(self):
            # Test failure dv sensitivity
            fail = constitutive.TestConstitutiveFailureStrainSens(self.con, self.elem_index, self.pt, self.x,
                                                                  self.dh, self.print_level, self.atol, self.rtol)
            self.assertFalse(fail)
