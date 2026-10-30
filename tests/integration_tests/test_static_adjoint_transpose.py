import os
import unittest

import numpy as np

from tacs import constitutive, elements, pytacs

"""
StaticProblem.solveAdjoint must return the solution of the transposed system,
K^T phi = rhs, on every row, the boundary-condition rows included: MPhys's
coupled adjoint (LinearBlockGS) checks exactly that residual through
addJacVecProduct(transpose=True), so a residual left on the bc rows stalls it.
The plate is clamped along one edge (three SPC nodes).
"""

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "./input_files/debug_plate.bdf")


def element_callback(
    dv_num, comp_id, comp_descript, elem_descripts, special_dvs, **kwargs
):
    prop = constitutive.MaterialProperties(rho=2780.0, E=73.1e9, nu=0.33, ys=324.0e6)
    con = constitutive.IsoShellConstitutive(prop, t=0.01, tNum=-1)
    transform = elements.ShellNaturalTransform()
    return elements.Quad4Shell(transform, con)


class StaticAdjointTransposeTest(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        fea_assembler = pytacs.pyTACS(bdf_file, options={"printLevel": 0})
        fea_assembler.initialize(element_callback)
        self.problem = fea_assembler.createStaticProblem("adjoint")
        mask = fea_assembler.assembler.createVec()
        mask.getArray()[:] = 1.0
        fea_assembler.assembler.applyBCs(mask)
        self.bc_rows = mask.getArray() == 0.0
        self.n = self.problem.getVariables().size
        self.rhs = np.random.default_rng(0).standard_normal(self.n)

    def test_solve_adjoint_satisfies_transposed_system_on_every_row(self):
        phi = np.zeros(self.n)
        self.problem.solveAdjoint(self.rhs, phi)
        product = np.zeros(self.n)
        self.problem.addJacVecProduct(phi, product, transpose=True)

        self.assertTrue(self.bc_rows.any())
        scale = np.linalg.norm(self.rhs)
        np.testing.assert_allclose(
            product[~self.bc_rows], self.rhs[~self.bc_rows], rtol=0.0, atol=1e-8 * scale
        )
        np.testing.assert_allclose(
            product[self.bc_rows], self.rhs[self.bc_rows], rtol=0.0, atol=1e-8 * scale
        )


if __name__ == "__main__":
    unittest.main()
