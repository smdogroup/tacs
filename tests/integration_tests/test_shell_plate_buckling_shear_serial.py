import unittest

# Import the module rather than the class so testflo does not also collect and
# rerun the parallel ProblemTest from this file.
import test_shell_plate_buckling_shear as shear

"""
Serial run of the shear-loaded plate buckling case in test_shell_plate_buckling_shear.py.

Shell buckling has a large null space of the geometric stiffness (in-plane
modes), whose image under the buckling shift-invert transform is the pole
theta = 1. Under complex step, a Ritz value whose real part lands exactly on 1.0
converts to an eigenvalue whose real part is exactly sigma, which used to sort
ahead of every physical eigenvalue and corrupt the sensitivity checks.

Whether that happens depends on round-off, which changes with the partition. On
the same hardware the serial and two-processor runs of this case have each
triggered it, so running both gives the regression coverage that either alone
does not.
"""


class ProblemTest(shear.ProblemTest):
    N_PROCS = 1  # this is how many MPI processes to use for this TestCase.


if __name__ == "__main__":
    unittest.main()
