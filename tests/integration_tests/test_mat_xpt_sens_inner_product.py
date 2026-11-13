import os
import unittest

import numpy as np
from mpi4py import MPI

from tacs import TACS, constitutive, elements, pytacs

"""
Test the Python binding of TACSAssembler::addMatXptSensInnerProduct.

The node sensitivity of psi^T A(X) phi is compared against central finite
differences of the assembled inner product, for the stiffness and mass
matrices of a small isotropic shell plate. This checks the binding (argument
order, node-vector layout, accumulation), not the element's own derivative.
"""

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BDF_FILE = os.path.join(BASE_DIR, "./input_files/debug_plate.bdf")

# Finite-difference step on the nodes and the agreement it allows
FD_STEP = 1e-6
FD_RTOL = 1e-5


def elemCallBack(dvNum, compID, compDescript, elemDescripts, globalDVs, **kwargs):
    """Isotropic aluminium shell with one thickness variable per component."""
    prop = constitutive.MaterialProperties(rho=2800.0, E=70e9, nu=0.3, ys=400e6)
    con = constitutive.IsoShellConstitutive(prop, t=0.01, tNum=dvNum)
    transform = elements.ShellRefAxisTransform(np.array([0.0, 1.0, 0.0]))
    return elements.Quad4Shell(transform, con)


class MatXptSensInnerProductTest(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.FEAAssembler = pytacs.pyTACS(BDF_FILE, MPI.COMM_WORLD)
        self.FEAAssembler.initialize(elemCallBack)
        problem = self.FEAAssembler.createStaticProblem("plate")
        # Set up the assembler's variables
        problem._updateAssemblerVars()
        self.assembler = self.FEAAssembler.assembler

        rng = np.random.default_rng(0)
        self.psi = self._randomVec(rng)
        self.phi = self._randomVec(rng)

        self.X0 = self.assembler.createNodeVec()
        self.assembler.getNodes(self.X0)
        self.direction = self.assembler.createNodeVec()
        self.direction.getArray()[:] = rng.uniform(
            -1.0, 1.0, self.direction.getArray().size
        )

    def _randomVec(self, rng):
        vec = self.assembler.createVec()
        vec.getArray()[:] = rng.uniform(-1.0, 1.0, vec.getArray().size)
        self.assembler.applyBCs(vec)
        return vec

    def _innerProduct(self, matType, step):
        """psi^T A(X0 + step * direction) phi, reassembled at the moved nodes."""
        X = self.assembler.createNodeVec()
        X.copyValues(self.X0)
        X.axpy(step, self.direction)
        self.assembler.setNodes(X)
        mat = self.assembler.createSchurMat()
        self.assembler.assembleMatType(matType, mat, TACS.NORMAL)
        product = self.assembler.createVec()
        mat.mult(self.phi, product)
        self.assembler.setNodes(self.X0)
        return self.psi.dot(product)

    def _sens(self, matType, calls=1):
        dfdX = self.assembler.createNodeVec()
        for _ in range(calls):
            self.assembler.addMatXptSensInnerProduct(
                1.0, matType, self.psi, self.phi, dfdX
            )
        dfdX.beginSetValues()
        dfdX.endSetValues()
        return dfdX

    def _checkFiniteDifference(self, matType):
        dfdX = self._sens(matType)
        fd = (
            self._innerProduct(matType, FD_STEP) - self._innerProduct(matType, -FD_STEP)
        ) / (2 * FD_STEP)
        np.testing.assert_allclose(
            np.real(dfdX.dot(self.direction)), np.real(fd), rtol=FD_RTOL
        )

    def _checkAccumulates(self, matType):
        once = self._sens(matType).getArray().copy()
        twice = self._sens(matType, calls=2).getArray().copy()
        np.testing.assert_allclose(twice, 2.0 * once, rtol=1e-14)

    def test_stiffness_matchesFiniteDifference(self):
        self._checkFiniteDifference(TACS.STIFFNESS_MATRIX)

    def test_mass_matchesFiniteDifference(self):
        self._checkFiniteDifference(TACS.MASS_MATRIX)

    def test_stiffness_accumulates(self):
        self._checkAccumulates(TACS.STIFFNESS_MATRIX)

    def test_mass_accumulates(self):
        self._checkAccumulates(TACS.MASS_MATRIX)


if __name__ == "__main__":
    unittest.main()
