import os
import unittest

import numpy as np
import scipy as sp
from mpi4py import MPI

from tacs import TACS, pytacs

"""
Test matrix-vector multiplication operations for TACSSchurMat and TACSParallelMat.

Uses the I_beam.bdf mesh to create a realistic stiffness matrix, then compares
TACS mult and multTranspose operations against scipy reference computed from
a pre-saved matrix file.

The reference matrix can be generated using the following code snippet but only if TACS is run in serial:
```
# Extract scipy sparse matrix (first block is the full matrix in serial)
# The scipy matrix returned by the ParallelMat matches the ordering of the
# results computed by TACS so we use that to generate the reference results
scipyParMat = sp.sparse.coo_array(parallel_mat.getMat()[0])
scipyParMat.eliminate_zeros()
sp.io.mmwrite("I_beam_stiffness_mat.mtx", scipyParMat, precision=16)
```
"""

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BDF_FILE = os.path.join(BASE_DIR, "./input_files/I_beam.bdf")
MAT_FILE = os.path.join(BASE_DIR, "./input_files/I_beam_stiffness_mat.mtx")


class MatrixOperationsTest(unittest.TestCase):
    """
    Test TACS matrix-vector operations against scipy reference.

    Strategy:
    1. Load pre-computed scipy reference matrix from file on all ranks
    2. Set up parallel TACS problem and assemble matrices
    3. Compare TACS mult/multTranspose against scipy @ x
    """

    N_PROCS = 1

    def setUp(self):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.rank
        self.dtype = TACS.dtype

        # Tolerances for comparison
        self.rtol = 1e-10
        self.atol = 1e-10

        # --- Step 1: Load scipy reference matrix from file ---
        self.K_scipy = sp.io.mmread(MAT_FILE).tocsr()

        # --- Step 2: Set up parallel TACS problem ---
        fea_assembler = pytacs.pyTACS(BDF_FILE, self.comm)
        fea_assembler.initialize()

        # Create a static problem (needed to set up vars)
        static_probs = fea_assembler.createTACSProbsFromBDF()
        problem = list(static_probs.values())[0]

        # Set up assembler vars
        problem._updateAssemblerVars()

        self.assembler = fea_assembler.assembler

        # Create both matrix types
        self.schur_mat = self.assembler.createSchurMat(TACS.NATURAL_ORDER)
        self.parallel_mat = self.assembler.createMat()

        # Assemble Jacobian into both matrices
        self.assembler.assembleJacobian(1.0, 0.0, 0.0, None, self.schur_mat)
        self.assembler.assembleJacobian(1.0, 0.0, 0.0, None, self.parallel_mat)

        # Get owner range for this rank (node indices, need to multiply by varsPerNode)
        # getOwnerRange returns array of size (num_procs+1) with node ranges for all procs
        owner_range = self.assembler.getOwnerRange()
        self.node_range = (owner_range[self.rank], owner_range[self.rank + 1])
        self.vars_per_node = self.assembler.getVarsPerNode()

        # Generate deterministic random input vector (same on all ranks)
        vec_size = self.K_scipy.shape[0]
        np.random.seed(42)
        self.x_global = np.random.rand(vec_size).astype(self.dtype)

        # Create TACS vector and set local portion
        self.x = self.assembler.createVec()
        low, high = self.node_range
        low_idx = low * self.vars_per_node
        high_idx = high * self.vars_per_node
        self.x.getArray()[:] = self.x_global[low_idx:high_idx]

    def _compareLocalResult(self, tacs_vec, global_ref, test_name):
        """Compare local portion of TACS vector against global reference."""
        low, high = self.node_range
        low_idx = low * self.vars_per_node
        high_idx = high * self.vars_per_node
        local_result = tacs_vec.getArray().copy()
        local_ref = global_ref[low_idx:high_idx]

        np.testing.assert_allclose(
            local_result,
            local_ref,
            rtol=self.rtol,
            atol=self.atol,
            err_msg=f"{test_name} failed on rank {self.rank}",
        )

    def templateTest(self, mat, testName, transpose=False):
        # Compute scipy reference result
        if transpose:
            y_ref = (self.K_scipy.T @ self.x_global).astype(self.dtype)
        else:
            y_ref = (self.K_scipy @ self.x_global).astype(self.dtype)

        # Compute TACS result
        y = self.assembler.createVec()
        if transpose:
            mat.multTranspose(self.x, y)
        else:
            mat.mult(self.x, y)

        self._compareLocalResult(y, y_ref, testName)

    def test_schur_mat_mult(self):
        """Test TACSSchurMat.mult against scipy matrix-vector product."""
        self.templateTest(self.schur_mat, "SchurMat.mult", transpose=False)

    def test_schur_mat_mult_transpose(self):
        """Test TACSSchurMat.multTranspose against scipy transpose product."""
        self.templateTest(self.schur_mat, "SchurMat.multTranspose", transpose=True)

    def test_parallel_mat_mult(self):
        """Test TACSParallelMat.mult against scipy matrix-vector product."""
        self.templateTest(self.parallel_mat, "ParallelMat.mult", transpose=False)

    def test_parallel_mat_mult_transpose(self):
        """Test TACSParallelMat.multTranspose against scipy transpose product."""
        self.templateTest(self.parallel_mat, "ParallelMat.multTranspose", transpose=True)


if __name__ == "__main__":
    unittest.main()
