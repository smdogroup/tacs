import unittest

import numpy as np

from tacs import TACS, elements, constitutive


class ModelTest(unittest.TestCase):
    def setUp(self):
        # fd/cs step size
        if TACS.dtype is complex:
            self.dh = 1e-200
            self.rtol = 1e-11
        else:
            self.dh = 1e-6
            self.rtol = 1e-2

        self.atol = np.clip(1e-5 * self.rtol, 1e-14, 1e-8)
        self.print_level = 0

        # Set element index
        self.elem_index = 0
        # Set the simulation time
        self.time = 0.0

        # Create the isotropic material
        rho = 1.0  # Density kg/m^3
        kappa = 1.0  # Thermal conductivity W/(m⋅K)
        cp = 1.0  # Specific heat J/(kg⋅K)
        lh = 10.0  # Latent heat J/kg
        Tm = 0.0  # Melting temperature (relative) K

        solid_prop = constitutive.MaterialProperties(
            rho=rho, kappa=2.0 * kappa, specific_heat=cp
        )
        liquid_prop = constitutive.MaterialProperties(
            rho=0.95 * rho, kappa=0.1 * kappa, specific_heat=1.1 * cp
        )

        # Set one thickness value for every component
        con = constitutive.PhaseChangeMaterialConstitutive(
            solid_prop, liquid_prop, lh=lh, Tm=Tm, t=0.1, tNum=-1
        )

        # Create the model for 2D
        self.model = elements.PCMHeatConduction2D(con)

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
