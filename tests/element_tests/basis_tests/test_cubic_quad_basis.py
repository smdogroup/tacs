import unittest

import numpy as np

from tacs import TACS, elements


class BasisTest(unittest.TestCase):
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

        # Create the basis functions for 2D
        self.basis = elements.CubicQuadBasis()

        # Seed random number generator in tacs for consistent test results
        elements.SeedRandomGenerator(0)

    def test_element_basis_functions(self):
        fail = elements.TestElementBasisFunctions(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)

    def test_element_basis_face_normals(self):
        fail = elements.TestElementBasisFaceNormals(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)

    def test_element_basis_jacobian_transform(self):
        fail = elements.TestElementBasisJacobianTransform(
            self.basis, self.dh, self.print_level, self.atol, self.rtol
        )
        self.assertFalse(fail)
