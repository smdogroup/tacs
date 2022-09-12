from tacs import TACS, elements, constitutive
import numpy as np
import unittest


class ModelTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-50
            self.rtol = 1e-11
        else:
            self.dh = 1e-6
            self.rtol = 1e-2

        # Basically, only check relative tolerance
        self.atol = 1e99
        self.print_level = 0

        # Set element index
        self.elem_index = 0
        # Set the simulation time
        self.time = 0.0

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

        # Create stiffness (need class)
        con = constitutive.PlaneStressConstitutive(self.props, t=1.0, tNum=0)

        # Create the model for 2D
        self.model = elements.LinearElasticity2D(con)

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_element_model_jacobian(self):
        fail = elements.TestElementModelJacobian(
            self.model,
            self.elem_index,
            self.time,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)

    def test_element_model_adj_xpt_sens_product(self):
        fail = elements.TestElementModelAdjXptSensProduct(
            self.model,
            self.elem_index,
            self.time,
            self.dh,
            self.print_level,
            self.atol,
            self.rtol,
        )
        self.assertFalse(fail)
