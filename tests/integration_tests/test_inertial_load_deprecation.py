"""
Tests for the deprecated ``inertiaVector`` argument of ``addInertialLoad``, which was
renamed to ``accelVec``.
"""

import os
import unittest
import warnings

import numpy as np
from mpi4py import MPI
from tacs import constitutive, elements, pytacs
from tacs._testing import fails_at_version

_DEPRECATION_REMOVAL_VERSION = "3.15"

base_dir = os.path.dirname(os.path.abspath(__file__))
bdf_file = os.path.join(base_dir, "input_files/plate.bdf")

ACCEL_VEC = np.array([0.0, 0.0, -9.81])


def _elemCallBack(dvNum, compID, compDescript, elemDescripts, specialDVs, **kwargs):
    """Return a simple isotropic shell element for every component."""
    props = constitutive.MaterialProperties(rho=2700.0, E=70e9, nu=0.3, ys=270e6)
    con = constitutive.IsoShellConstitutive(props, t=0.01, tNum=dvNum)
    transform = elements.ShellNaturalTransform()
    return elements.Quad4Shell(transform, con)


class TestInertiaVectorDeprecation(unittest.TestCase):
    """Verify backward-compat behaviour of the deprecated ``inertiaVector`` argument."""

    @classmethod
    def setUpClass(cls):
        cls.fea = pytacs.pyTACS(
            bdf_file, comm=MPI.COMM_WORLD, options={"printLevel": 0}
        )
        cls.fea.initialize(_elemCallBack)

    def _solveStatic(self, **inertialLoadKwargs):
        """Solve a static problem with a single inertial load and return the state vector."""
        problem = self.fea.createStaticProblem("grav")
        problem.addInertialLoad(**inertialLoadKwargs)
        problem.solve()
        return problem.getVariables()

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_static_inertia_vector_warns_and_matches_accel_vec(self):
        """
        Given a static problem,
        when addInertialLoad is called with the deprecated inertiaVector keyword,
        then a DeprecationWarning is emitted and the solution matches one built with accelVec.
        """
        with self.assertWarns(DeprecationWarning):
            oldStates = self._solveStatic(inertiaVector=ACCEL_VEC)
        newStates = self._solveStatic(accelVec=ACCEL_VEC)
        np.testing.assert_allclose(oldStates, newStates)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_transient_inertia_vector_warns(self):
        """Assert that TransientProblem.addInertialLoad emits a DeprecationWarning when given inertiaVector."""
        problem = self.fea.createTransientProblem(
            "grav", tInit=0.0, tFinal=1.0, numSteps=2
        )
        with self.assertWarns(DeprecationWarning):
            problem.addInertialLoad(0, inertiaVector=ACCEL_VEC)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_both_accel_vec_and_inertia_vector_raises(self):
        """Assert that passing both accelVec and the deprecated inertiaVector raises a ValueError."""
        problem = self.fea.createStaticProblem("grav")
        with self.assertWarns(DeprecationWarning):
            with self.assertRaisesRegex(ValueError, "Cannot specify both"):
                problem.addInertialLoad(accelVec=ACCEL_VEC, inertiaVector=ACCEL_VEC)

    @fails_at_version(_DEPRECATION_REMOVAL_VERSION)
    def test_accel_vec_no_warning(self):
        """Assert that addInertialLoad emits no DeprecationWarning when the acceleration is passed positionally or as accelVec."""
        problem = self.fea.createStaticProblem("grav")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            problem.addInertialLoad(ACCEL_VEC)
            problem.addInertialLoad(accelVec=ACCEL_VEC)
        deprecationWarnings = [
            w for w in caught if issubclass(w.category, DeprecationWarning)
        ]
        self.assertListEqual(deprecationWarnings, [])


if __name__ == "__main__":
    unittest.main()
