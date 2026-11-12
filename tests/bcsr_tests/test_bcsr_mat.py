"""
Correctness tests of the BCSRMat block sparse matrix kernels (through ``TACS.BCSRMatrix``).

Every operation is checked for block sizes 1 through 9 against dense numpy references or
small pure-Python reference implementations that do not depend on the TACS kernels.
"""

from __future__ import annotations

import unittest

import numpy as np

from tacs import TACS

RTOL = 1e-10
SEED = 20261111
IS_COMPLEX = np.issubdtype(TACS.dtype, np.complexfloating)
FULL_LEV_FILL = 1000
FILL_ESTIMATE = 10.0

def makeRandomArray(
    rng: np.random.Generator, shape: int | tuple[int, ...]
) -> np.ndarray:
    """Return a random array of TACS.dtype (with a random imaginary part in complex mode)."""
    values = rng.standard_normal(shape).astype(TACS.dtype)
    if IS_COMPLEX:
        values = values + 1j * rng.standard_normal(shape)
    return values


def makeMask(
    rng: np.random.Generator, numRows: int, numCols: int, density: float = 0.4
) -> np.ndarray:
    """Return a random boolean block pattern whose row 2 is empty."""
    mask = rng.random((numRows, numCols)) < density
    mask[2, :] = False
    return mask


def makeSquareMask(
    rng: np.random.Generator, numRows: int, density: float = 0.35
) -> np.ndarray:
    """Return a random square block pattern with a full diagonal and a diagonal-only row 1."""
    mask = rng.random((numRows, numRows)) < density
    mask[1, :] = False
    np.fill_diagonal(mask, True)
    return mask


def expandBlocks(mask: np.ndarray, blocks: np.ndarray, bsize: int) -> np.ndarray:
    """Assemble the dense matrix from blocks stored in row-major pattern order."""
    dense = np.zeros((mask.shape[0] * bsize, mask.shape[1] * bsize), dtype=TACS.dtype)
    for index, (ii, jj) in enumerate(zip(*np.nonzero(mask), strict=True)):
        dense[ii * bsize : (ii + 1) * bsize, jj * bsize : (jj + 1) * bsize] = blocks[
            index
        ]
    return dense


def makeMatrix(
    rng: np.random.Generator,
    bsize: int,
    mask: np.ndarray,
    isDominant: bool = False,
    numThreads: int = 1,
) -> tuple[TACS.BCSRMatrix, np.ndarray]:
    """Create a BCSRMatrix with random values on `mask` and return it with its dense equivalent."""
    blocks = makeRandomArray(rng, (int(mask.sum()), bsize, bsize))
    if isDominant:
        shift = np.abs(expandBlocks(mask, blocks, bsize)).sum(axis=1).max() + 1.0
        isDiag = np.nonzero(mask)[0] == np.nonzero(mask)[1]
        blocks[isDiag] += shift * np.eye(bsize)
    rowp = np.concatenate(([0], np.cumsum(mask.sum(axis=1))))
    cols = np.nonzero(mask)[1]
    mat = TACS.BCSRMatrix(
        bsize, rowp, cols, ncols=mask.shape[1], A=blocks, num_threads=numThreads
    )
    return mat, expandBlocks(mask, blocks, bsize)


def getDense(mat: TACS.BCSRMatrix) -> np.ndarray:
    """Return a dense copy of the values currently stored in `mat`."""
    return np.array(mat.getArrays().toarray())


def getMask(mat: TACS.BCSRMatrix) -> np.ndarray:
    """Return the boolean block pattern of `mat`."""
    bsr = mat.getArrays()
    mask = np.zeros((mat.getRowDim(), mat.getColDim()), dtype=bool)
    for ii in range(mat.getRowDim()):
        mask[ii, bsr.indices[bsr.indptr[ii] : bsr.indptr[ii + 1]]] = True
    return mask


def expandMask(mask: np.ndarray, bsize: int) -> np.ndarray:
    """Expand a block pattern to a scalar pattern."""
    return np.kron(mask, np.ones((bsize, bsize), dtype=bool))


def getBlock(dense: np.ndarray, bsize: int, ii: int, jj: int) -> np.ndarray:
    """Return a view of block (ii, jj) of a dense matrix."""
    return dense[ii * bsize : (ii + 1) * bsize, jj * bsize : (jj + 1) * bsize]


def computeReferenceILU(dense: np.ndarray, mask: np.ndarray, bsize: int) -> np.ndarray:
    """
    Block ILU restricted to `mask`, returned in the BCSRMat storage convention.

    Strictly-lower blocks hold L, diagonal blocks the inverse of U's diagonal, upper blocks U.
    """
    factor = dense.copy()
    numRows = mask.shape[0]
    for ii in range(numRows):
        for kk in range(ii):
            if mask[ii, kk]:
                lower = getBlock(factor, bsize, ii, kk) @ np.linalg.inv(
                    getBlock(factor, bsize, kk, kk)
                )
                for jj in range(kk + 1, numRows):
                    if mask[ii, jj] and mask[kk, jj]:
                        getBlock(factor, bsize, ii, jj)[:] -= lower @ getBlock(
                            factor, bsize, kk, jj
                        )
                getBlock(factor, bsize, ii, kk)[:] = lower
    for ii in range(numRows):
        getBlock(factor, bsize, ii, ii)[:] = np.linalg.inv(
            getBlock(factor, bsize, ii, ii)
        )
    return factor


def splitStoredFactor(stored: np.ndarray, bsize: int) -> tuple[np.ndarray, np.ndarray]:
    """Return dense (L, U) from a factor in the BCSRMat storage convention."""
    numRows = stored.shape[0] // bsize
    isLower = expandMask(np.tril(np.ones((numRows, numRows), dtype=bool), -1), bsize)
    isUpper = expandMask(np.triu(np.ones((numRows, numRows), dtype=bool), 1), bsize)
    lower = np.where(isLower, stored, 0.0) + np.eye(stored.shape[0])
    upper = np.where(isUpper, stored, 0.0)
    for ii in range(numRows):
        getBlock(upper, bsize, ii, ii)[:] = np.linalg.inv(
            getBlock(stored, bsize, ii, ii)
        )
    return lower, upper


def computeReferenceLowerFactor(
    lower: np.ndarray,
    lowerMask: np.ndarray,
    denseE: np.ndarray,
    maskE: np.ndarray,
    bsize: int,
):
    """Return L^{-1} E with every update restricted to the pattern of E."""
    result = denseE.copy()
    for ii in range(lowerMask.shape[0]):
        for kk in range(ii):
            if lowerMask[ii, kk]:
                for jj in np.nonzero(maskE[ii] & maskE[kk])[0]:
                    getBlock(result, bsize, ii, jj)[:] -= getBlock(
                        lower, bsize, ii, kk
                    ) @ getBlock(result, bsize, kk, jj)
    return result


def computeReferenceUpperFactor(
    upper: np.ndarray,
    upperMask: np.ndarray,
    denseF: np.ndarray,
    maskF: np.ndarray,
    bsize: int,
):
    """Return F U^{-1} with every update restricted to the pattern of F."""
    result = denseF.copy()
    for ii in range(maskF.shape[0]):
        for jj in np.nonzero(maskF[ii])[0]:
            scaled = getBlock(result, bsize, ii, jj) @ np.linalg.inv(
                getBlock(upper, bsize, jj, jj)
            )
            for kk in np.nonzero(maskF[ii] & upperMask[jj])[0]:
                if kk > jj:
                    getBlock(result, bsize, ii, kk)[:] -= scaled @ getBlock(
                        upper, bsize, jj, kk
                    )
            getBlock(result, bsize, ii, jj)[:] = scaled
    return result


def applyReferenceSOR(
    dense: np.ndarray,
    bsize: int,
    start: int,
    end: int,
    omega: float,
    rhs: np.ndarray,
    x: np.ndarray,
    diagShift: np.ndarray | None = None,
    denseB: np.ndarray | None = None,
    varOffset: int = 0,
    xext: np.ndarray | None = None,
) -> None:
    """One SOR sweep on `x` in place: rows start..end-1 forwards, or start-1..end backwards if start > end."""
    rows = range(start, end) if start < end else range(start - 1, end - 1, -1)
    for ii in rows:
        rowSlice = slice(ii * bsize, (ii + 1) * bsize)
        diagBlock = dense[rowSlice, rowSlice]
        residual = rhs[rowSlice] - dense[rowSlice] @ x + diagBlock @ x[rowSlice]
        if denseB is not None and ii >= varOffset:
            rowB = ii - varOffset
            residual -= denseB[rowB * bsize : (rowB + 1) * bsize] @ xext
        if diagShift is not None:
            diagBlock = diagBlock + np.diag(diagShift[rowSlice])
        x[rowSlice] = (1.0 - omega) * x[rowSlice] + omega * np.linalg.solve(
            diagBlock, residual
        )


class TestBCSRMatKernels(unittest.TestCase):
    """Serial kernels for every block size from 1 to 9 (7 and 9 use the generic code)."""

    BLOCK_SIZES = tuple(range(1, 10))
    NUM_THREADS = 1
    NUM_ROWS = 9
    NUM_COLS = 7  # column count of the rectangular matrices

    def assertClose(self, actual: np.ndarray, expected: np.ndarray) -> None:
        """Assert agreement to a relative tolerance scaled by the largest expected magnitude."""
        scale = max(1.0, float(np.abs(expected).max(initial=0.0)))
        np.testing.assert_allclose(actual, expected, rtol=RTOL, atol=RTOL * scale)

    def makeRng(self, bsize: int) -> np.random.Generator:
        """Return the seeded generator for one block size."""
        return np.random.default_rng(SEED + bsize)

    def makeRect(
        self, rng: np.random.Generator, bsize: int, numRows: int, numCols: int
    ):
        """Return a random rectangular (matrix, dense, mask) triple with an empty row."""
        mask = makeMask(rng, numRows, numCols)
        mat, dense = makeMatrix(rng, bsize, mask, numThreads=self.NUM_THREADS)
        return mat, dense, mask

    def makeFactored(self, rng: np.random.Generator, bsize: int, levFill: int):
        """Return (factored ILU(levFill) matrix, original dense matrix, original matrix)."""
        mask = makeSquareMask(rng, self.NUM_ROWS)
        mat, dense = makeMatrix(
            rng, bsize, mask, isDominant=True, numThreads=self.NUM_THREADS
        )
        factor = TACS.BCSRMatrix.createILU(mat, levFill, FILL_ESTIMATE)
        factor.copyValues(mat)
        factor.factor()
        return factor, dense, mat

    def test_mult_and_mult_add(self):
        """
        Given a random rectangular matrix with an empty row,
        when mult and multAdd with z aliased to y are called,
        then the results equal the dense A x and A x + z.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mat, dense, _ = self.makeRect(rng, bsize, self.NUM_ROWS, self.NUM_COLS)
                x = makeRandomArray(rng, self.NUM_COLS * bsize)
                z = makeRandomArray(rng, self.NUM_ROWS * bsize)
                y = makeRandomArray(rng, self.NUM_ROWS * bsize)
                mat.mult(x, y)
                self.assertClose(y, dense @ x)
                zCopy = z.copy()
                mat.multAdd(x, z, z)
                self.assertClose(z, dense @ x + zCopy)

    def test_mult_add_distinct_arrays(self):
        """
        Given a random rectangular matrix with an empty row,
        when multAdd(x, z, y) is called with z and y distinct arrays,
        then y equals the dense A x + z and z is unchanged.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mat, dense, _ = self.makeRect(rng, bsize, self.NUM_ROWS, self.NUM_COLS)
                x = makeRandomArray(rng, self.NUM_COLS * bsize)
                z = makeRandomArray(rng, self.NUM_ROWS * bsize)
                y = makeRandomArray(rng, self.NUM_ROWS * bsize)
                zCopy = z.copy()
                mat.multAdd(x, z, y)
                self.assertClose(y, dense @ x + zCopy)
                self.assertClose(z, zCopy)

    def test_mult_transpose_and_add(self):
        """
        Given a random rectangular matrix with an empty row,
        when multTranspose and multTransposeAdd (z distinct from y and z aliased to y) are called,
        then the results equal the dense A^T x and A^T x + z.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mat, dense, _ = self.makeRect(rng, bsize, self.NUM_ROWS, self.NUM_COLS)
                x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                z = makeRandomArray(rng, self.NUM_COLS * bsize)
                y = makeRandomArray(rng, self.NUM_COLS * bsize)
                mat.multTranspose(x, y)
                self.assertClose(y, dense.T @ x)
                zCopy = z.copy()
                mat.multTransposeAdd(x, z, y)
                self.assertClose(y, dense.T @ x + zCopy)
                mat.multTransposeAdd(x, z, z)
                self.assertClose(z, dense.T @ x + zCopy)

    def test_factor_full_fill(self):
        """
        Given a block-diagonally dominant matrix copied into its full-fill ILU pattern,
        when it is factored and applyFactor is called in two-vector and in-place form,
        then the stored L and U multiply back to A and both solves equal numpy.linalg.solve.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                factor, dense, _ = self.makeFactored(rng, bsize, FULL_LEV_FILL)
                lower, upper = splitStoredFactor(getDense(factor), bsize)
                self.assertClose(lower @ upper, dense)
                rhs = makeRandomArray(rng, self.NUM_ROWS * bsize)
                expected = np.linalg.solve(dense, rhs)
                y = np.zeros_like(rhs)
                rhsCopy = rhs.copy()
                factor.applyFactor(rhs, y)
                self.assertClose(y, expected)
                self.assertClose(rhs, rhsCopy)
                factor.applyFactor(rhs)
                self.assertClose(rhs, expected)

    def test_factor_incomplete_fill(self):
        """
        Given a block-diagonally dominant matrix copied into its ILU(0) or ILU(1) pattern,
        when it is factored,
        then the stored blocks equal a reference block ILU restricted to the same pattern.
        """
        for bsize in self.BLOCK_SIZES:
            for levFill in (0, 1):
                with self.subTest(bsize=bsize, levFill=levFill):
                    rng = self.makeRng(bsize)
                    factor, dense, _ = self.makeFactored(rng, bsize, levFill)
                    self.assertClose(
                        getDense(factor),
                        computeReferenceILU(dense, getMask(factor), bsize),
                    )

    def test_apply_lower_and_upper(self):
        """
        Given an ILU(0)-factored matrix,
        when applyLower and applyUpper are called with separate and with aliased input/output arrays,
        then the results equal dense solves with the L and U built from the stored factor.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                factor, _, _ = self.makeFactored(rng, bsize, 0)
                lower, upper = splitStoredFactor(getDense(factor), bsize)
                for method, triangle in (
                    (factor.applyLower, lower),
                    (factor.applyUpper, upper),
                ):
                    x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                    expected = np.linalg.solve(triangle, x)
                    y = np.zeros_like(x)
                    method(x, y)
                    self.assertClose(y, expected)
                    method(x, x)
                    self.assertClose(x, expected)

    def test_apply_partial_lower_and_upper(self):
        """
        Given an ILU(1)-factored matrix and several row offsets,
        when applyPartialLower and applyPartialUpper are applied to a vector of the trailing rows,
        then the vector equals the dense solve with the trailing diagonal sub-block of L or U.
        """
        for bsize in self.BLOCK_SIZES:
            for varOffset in (0, 3, self.NUM_ROWS - 1):
                with self.subTest(bsize=bsize, varOffset=varOffset):
                    rng = self.makeRng(bsize)
                    factor, _, _ = self.makeFactored(rng, bsize, 1)
                    lower, upper = splitStoredFactor(getDense(factor), bsize)
                    offset = varOffset * bsize
                    for method, triangle in (
                        (factor.applyPartialLower, lower),
                        (factor.applyPartialUpper, upper),
                    ):
                        x = makeRandomArray(rng, (self.NUM_ROWS - varOffset) * bsize)
                        expected = np.linalg.solve(triangle[offset:, offset:], x)
                        method(x, varOffset)
                        self.assertClose(x, expected)

    def test_apply_factor_schur(self):
        """
        Given an ILU(1)-factored matrix and a row offset v,
        when applyFactorSchur is applied to a full-length vector x,
        then x[:v] equals U[:v, :v]^{-1} (x[:v] - U[:v, v:] x[v:]) and x[v:] is unchanged.
        """
        for bsize in self.BLOCK_SIZES:
            for varOffset in (1, 4, self.NUM_ROWS):
                with self.subTest(bsize=bsize, varOffset=varOffset):
                    rng = self.makeRng(bsize)
                    factor, _, _ = self.makeFactored(rng, bsize, 1)
                    _, upper = splitStoredFactor(getDense(factor), bsize)
                    offset = varOffset * bsize
                    x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                    expected = x.copy()
                    expected[:offset] = np.linalg.solve(
                        upper[:offset, :offset],
                        x[:offset] - upper[:offset, offset:] @ x[offset:],
                    )
                    factor.applyFactorSchur(x, varOffset)
                    self.assertClose(x, expected)

    def test_apply_lower_factor(self):
        """
        Given an ILU(0)-factored matrix and a random rectangular matrix E with an empty row,
        when applyLowerFactor(E) is called,
        then E equals the reference L^{-1} E with updates restricted to E's pattern.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                factor, _, _ = self.makeFactored(rng, bsize, 0)
                matE, denseE, maskE = self.makeRect(
                    rng, bsize, self.NUM_ROWS, self.NUM_COLS
                )
                lower, _ = splitStoredFactor(getDense(factor), bsize)
                expected = computeReferenceLowerFactor(
                    lower, getMask(factor), denseE, maskE, bsize
                )
                factor.applyLowerFactor(matE)
                self.assertClose(getDense(matE), expected)

    def test_apply_upper_factor(self):
        """
        Given an ILU(0)-factored matrix and a random rectangular matrix F with an empty row,
        when applyUpperFactor(F) is called,
        then F equals the reference F U^{-1} with updates restricted to F's pattern.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                factor, _, _ = self.makeFactored(rng, bsize, 0)
                matF, denseF, maskF = self.makeRect(
                    rng, bsize, self.NUM_COLS, self.NUM_ROWS
                )
                _, upper = splitStoredFactor(getDense(factor), bsize)
                expected = computeReferenceUpperFactor(
                    upper, getMask(factor), denseF, maskF, bsize
                )
                factor.applyUpperFactor(matF)
                self.assertClose(getDense(matF), expected)

    def test_mat_mult_add(self):
        """
        Given rectangular A and B and a matrix C with its own random pattern and initial values,
        when C.matMultAdd(alpha, A, B) is called for alpha = 1, -1 and 0.37,
        then C equals its initial values plus alpha*A*B restricted to C's pattern.
        """
        numInner = self.NUM_COLS + 1
        for bsize in self.BLOCK_SIZES:
            for alpha in (1.0, -1.0, 0.37):
                with self.subTest(bsize=bsize, alpha=alpha):
                    rng = self.makeRng(bsize)
                    matA, denseA, _ = self.makeRect(rng, bsize, self.NUM_ROWS, numInner)
                    matB, denseB, _ = self.makeRect(rng, bsize, numInner, self.NUM_COLS)
                    matC, denseC, maskC = self.makeRect(
                        rng, bsize, self.NUM_ROWS, self.NUM_COLS
                    )
                    matC.matMultAdd(alpha, matA, matB)
                    expected = denseC + alpha * np.where(
                        expandMask(maskC, bsize), denseA @ denseB, 0.0
                    )
                    self.assertClose(getDense(matC), expected)

    def test_mat_mult_normal(self):
        """
        Given a rectangular matrix B and the normal-equations pattern created from it,
        when matMultNormal(s, B) is called on a matrix holding stale values,
        then the matrix equals the dense B^T diag(s) B.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                matB, denseB, _ = self.makeRect(
                    rng, bsize, self.NUM_ROWS, self.NUM_COLS
                )
                normal = TACS.BCSRMatrix.createNormal(matB, FILL_ESTIMATE)
                normal.setValues(np.ones_like(normal.getArrays().data))
                scale = makeRandomArray(rng, self.NUM_ROWS * bsize)
                normal.matMultNormal(scale, matB)
                self.assertClose(getDense(normal), denseB.T @ np.diag(scale) @ denseB)

    def test_apply_sor_full_sweeps(self):
        """
        Given a dominant matrix whose diagonal has been factored with and without a diagonal shift,
        when applySOR runs two forward sweeps,
        then x equals two sweeps of the reference block SOR.
        """
        omega = 1.2
        for bsize in self.BLOCK_SIZES:
            for doShift in (False, True):
                with self.subTest(bsize=bsize, doShift=doShift):
                    rng = self.makeRng(bsize)
                    mat, dense = makeMatrix(
                        rng,
                        bsize,
                        makeSquareMask(rng, self.NUM_ROWS),
                        True,
                        self.NUM_THREADS,
                    )
                    diagShift = (
                        np.abs(makeRandomArray(rng, self.NUM_ROWS * bsize)).astype(
                            TACS.dtype
                        )
                        if doShift
                        else None
                    )
                    mat.factorDiag(diagShift)
                    rhs = makeRandomArray(rng, self.NUM_ROWS * bsize)
                    x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                    expected = x.copy()
                    for _ in range(2):
                        applyReferenceSOR(
                            dense,
                            bsize,
                            0,
                            self.NUM_ROWS,
                            omega,
                            rhs,
                            expected,
                            diagShift,
                        )
                    mat.applySOR(rhs, x, omega, 2)
                    self.assertClose(x, expected)

    def test_apply_sor_range(self):
        """
        Given a dominant matrix with a factored diagonal, with and without an external coupling matrix B,
        when applySORRange sweeps a forward interval (start < end) and a reverse interval (start > end),
        then x equals the reference sweep over rows start..end-1 or start-1..end respectively.
        """
        omega = 0.8
        varOffset = 3
        intervals = (
            (0, self.NUM_ROWS),
            (2, self.NUM_ROWS - 1),
            (self.NUM_ROWS, 0),
            (self.NUM_ROWS - 1, 2),
        )
        for bsize in self.BLOCK_SIZES:
            for start, end in intervals:
                for doUseB in (False, True):
                    with self.subTest(bsize=bsize, start=start, end=end, doUseB=doUseB):
                        rng = self.makeRng(bsize)
                        mat, dense = makeMatrix(
                            rng,
                            bsize,
                            makeSquareMask(rng, self.NUM_ROWS),
                            True,
                            self.NUM_THREADS,
                        )
                        mat.factorDiag()
                        matB, denseB, _ = self.makeRect(
                            rng, bsize, self.NUM_ROWS - varOffset, self.NUM_COLS
                        )
                        xext = makeRandomArray(rng, self.NUM_COLS * bsize)
                        rhs = makeRandomArray(rng, self.NUM_ROWS * bsize)
                        x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                        expected = x.copy()
                        if doUseB:
                            applyReferenceSOR(
                                dense,
                                bsize,
                                start,
                                end,
                                omega,
                                rhs,
                                expected,
                                None,
                                denseB,
                                varOffset,
                                xext,
                            )
                            mat.applySORRange(
                                matB, start, end, varOffset, omega, rhs, xext, x
                            )
                        else:
                            applyReferenceSOR(
                                dense, bsize, start, end, omega, rhs, expected
                            )
                            mat.applySORRange(None, start, end, 0, omega, rhs, None, x)
                        self.assertClose(x, expected)

    def test_schur_complement_sequence(self):
        """
        Given a dominant matrix partitioned into [[B, E], [F, C]] and the full-fill Schur patterns,
        when Bpc is factored and applyLowerFactor(Epc), applyUpperFactor(Fpc), Sc.matMultAdd(-1, Fpc, Epc) are run,
        then Sc equals the dense Schur complement C - F B^{-1} E.
        """
        numB = self.NUM_ROWS // 2
        offset = None
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mat, dense = makeMatrix(
                    rng,
                    bsize,
                    makeSquareMask(rng, self.NUM_ROWS),
                    True,
                    self.NUM_THREADS,
                )
                offset = numB * bsize
                parts = mat.partition(numB)
                self.assertClose(getDense(parts[0]), dense[:offset, :offset])
                self.assertClose(getDense(parts[1]), dense[:offset, offset:])
                self.assertClose(getDense(parts[2]), dense[offset:, :offset])
                self.assertClose(getDense(parts[3]), dense[offset:, offset:])
                factors = TACS.BCSRMatrix.createSchur(
                    *parts, FULL_LEV_FILL, FILL_ESTIMATE, True
                )
                for factorMat, part in zip(factors, parts, strict=True):
                    factorMat.copyValues(part)
                matBpc, matEpc, matFpc, matSc = factors
                matBpc.factor()
                matBpc.applyLowerFactor(matEpc)
                matBpc.applyUpperFactor(matFpc)
                matSc.matMultAdd(-1.0, matFpc, matEpc)
                expected = dense[offset:, offset:] - dense[
                    offset:, :offset
                ] @ np.linalg.solve(dense[:offset, :offset], dense[:offset, offset:])
                self.assertClose(getDense(matSc), expected)

    def test_mat_mult_pattern(self):
        """
        Given S, A and B with random patterns,
        when the symbolic C = S + A*B constructor is used,
        then C is zero-valued and its pattern is the union of S's pattern and the structural product of A and B.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                matA, _, maskA = self.makeRect(rng, bsize, self.NUM_ROWS, self.NUM_COLS)
                matB, _, maskB = self.makeRect(rng, bsize, self.NUM_COLS, self.NUM_ROWS)
                matS, _, maskS = self.makeRect(rng, bsize, self.NUM_ROWS, self.NUM_ROWS)
                matC = TACS.BCSRMatrix.createMatMultPattern(
                    matS, matA, matB, FILL_ESTIMATE
                )
                expected = maskS | ((maskA.astype(int) @ maskB.astype(int)) > 0)
                np.testing.assert_array_equal(getMask(matC), expected)
                self.assertEqual(np.abs(getDense(matC)).max(), 0.0)

    def test_value_operations(self):
        """
        Given two random matrices X and Y on the same pattern,
        when copyValues, scale, axpy, setValues, addDiag and zeroEntries are applied,
        then the stored values equal the corresponding dense expressions and isEqual tracks equality.
        """
        alpha = 0.37
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mask = makeSquareMask(rng, self.NUM_ROWS)
                matX, denseX = makeMatrix(rng, bsize, mask, numThreads=self.NUM_THREADS)
                matY, denseY = makeMatrix(rng, bsize, mask, numThreads=self.NUM_THREADS)
                self.assertEqual(
                    (matX.getBlockSize(), matX.getRowDim(), matX.getColDim()),
                    (bsize, self.NUM_ROWS, self.NUM_ROWS),
                )
                self.assertFalse(matY.isEqual(matX))
                matY.axpy(alpha, matX)
                self.assertClose(getDense(matY), denseY + alpha * denseX)
                matY.copyValues(matX)
                self.assertTrue(matY.isEqual(matX))
                matY.scale(alpha)
                self.assertClose(getDense(matY), alpha * denseX)
                matY.setValues(matX.getArrays().data)
                self.assertClose(getDense(matY), denseX)
                matY.factorDiag()  # sets up the diagonal index that addDiag needs
                matY.addDiag(alpha)
                self.assertClose(
                    getDense(matY), denseX + alpha * np.eye(self.NUM_ROWS * bsize)
                )
                matY.zeroEntries()
                self.assertEqual(np.abs(getDense(matY)).max(), 0.0)

    def test_axpby(self):
        """
        Given two random matrices X and Y on the same pattern,
        when Y.axpby(alpha, beta, X) is called,
        then Y equals the dense alpha*X + beta*Y.
        """
        alpha = 0.37
        beta = -1.3
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                mask = makeSquareMask(rng, self.NUM_ROWS)
                matX, denseX = makeMatrix(rng, bsize, mask, numThreads=self.NUM_THREADS)
                matY, denseY = makeMatrix(rng, bsize, mask, numThreads=self.NUM_THREADS)
                matY.axpby(alpha, beta, matX)
                self.assertClose(getDense(matY), alpha * denseX + beta * denseY)

    def test_axpby_larger_pattern(self):
        """
        Given a matrix Y whose pattern contains the pattern of X plus extra blocks,
        when Y.axpby(alpha, beta, X) is called,
        then Y equals the dense alpha*X + beta*Y, so the extra blocks are only scaled by beta.
        """
        alpha = 0.37
        beta = -1.3
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                rng = self.makeRng(bsize)
                maskY = makeSquareMask(rng, self.NUM_ROWS)
                maskX = maskY & (rng.random(maskY.shape) < 0.5)
                matX, denseX = makeMatrix(rng, bsize, maskX, numThreads=self.NUM_THREADS)
                matY, denseY = makeMatrix(rng, bsize, maskY, numThreads=self.NUM_THREADS)
                matY.axpby(alpha, beta, matX)
                self.assertClose(getDense(matY), alpha * denseX + beta * denseY)

    def test_generic_matches_block_implementation(self):
        """
        Given a full-fill factorization and rectangular operands computed with the block-specific kernels,
        when the same mult, multTranspose, factor, applyFactor and matMultAdd are repeated after initGenericImpl,
        then the generic results agree with the block-specific ones.
        """
        for bsize in self.BLOCK_SIZES:
            with self.subTest(bsize=bsize):
                results = []
                for doUseGeneric in (False, True):
                    rng = self.makeRng(bsize)
                    mask = makeSquareMask(rng, self.NUM_ROWS)
                    mat, _ = makeMatrix(
                        rng, bsize, mask, isDominant=True, numThreads=self.NUM_THREADS
                    )
                    factor = TACS.BCSRMatrix.createILU(mat, 1, FILL_ESTIMATE)
                    product = TACS.BCSRMatrix.createMatMultPattern(
                        mat, mat, mat, FILL_ESTIMATE
                    )
                    if doUseGeneric:
                        for matrix in (mat, factor, product):
                            matrix.initGenericImpl()
                    x = makeRandomArray(rng, self.NUM_ROWS * bsize)
                    y = np.zeros_like(x)
                    yTrans = np.zeros_like(x)
                    mat.mult(x, y)
                    mat.multTranspose(x, yTrans)
                    factor.copyValues(mat)
                    factor.factor()
                    solution = np.zeros_like(x)
                    factor.applyFactor(x, solution)
                    product.matMultAdd(0.37, mat, mat)
                    results.append(
                        (y, yTrans, getDense(factor), solution, getDense(product))
                    )
                for blockResult, genericResult in zip(*results, strict=True):
                    self.assertClose(genericResult, blockResult)


class TestBCSRMatThreadedKernels(TestBCSRMatKernels):
    """
    The same checks with two threads for block sizes 6 and 8, the only sizes with threaded kernels.

    BCSRMat uses the pthread kernels (mult/multAdd, factor, applyFactor, matMultAdd, applyLowerFactor,
    applyUpperFactor) whenever a threaded kernel exists and the thread count exceeds one; rows are handed
    out in groups of 16 (matrix-vector) or 4 (matrix-matrix), so 40 rows give each scheduler several jobs.
    """

    BLOCK_SIZES = (6, 8)
    NUM_THREADS = 2
    NUM_ROWS = 40
    NUM_COLS = 37


class TestBCSRMatrixValidation(unittest.TestCase):
    """Argument validation of the Python wrapper."""

    def test_wrong_vector_length_raises(self):
        """Assert that vector arguments of the wrong length raise ValueError instead of reaching the C++ kernels."""
        rng = np.random.default_rng(SEED)
        bsize = 3
        numRows = 5
        mat, _ = makeMatrix(rng, bsize, makeSquareMask(rng, numRows), isDominant=True)
        good = np.zeros(numRows * bsize, dtype=TACS.dtype)
        bad = np.zeros(numRows * bsize - 1, dtype=TACS.dtype)
        other = np.zeros(numRows * bsize, dtype=TACS.dtype)
        mat.factorDiag()
        calls = (
            lambda: mat.mult(bad, good),
            lambda: mat.mult(good, bad),
            lambda: mat.multAdd(good, bad, other),
            lambda: mat.multTranspose(bad, good),
            lambda: mat.multTransposeAdd(good, other, bad),
            lambda: mat.applyFactor(bad),
            lambda: mat.applyFactor(good, bad),
            lambda: mat.applyLower(good, bad),
            lambda: mat.applyUpper(bad, good),
            lambda: mat.applyPartialLower(good, 1),
            lambda: mat.applyPartialUpper(good, 1),
            lambda: mat.applyFactorSchur(bad, 1),
            lambda: mat.factorDiag(bad),
            lambda: mat.applySOR(bad, good, 1.0, 1),
            lambda: mat.applySORRange(None, 0, numRows + 1, 0, 1.0, good, None, other),
            lambda: mat.matMultNormal(bad, mat),
            lambda: mat.mult(good, good),
        )
        for index, call in enumerate(calls):
            with self.subTest(call=index):
                self.assertRaises(ValueError, call)

    def test_invalid_pattern_raises(self):
        """Assert that unsorted or out-of-range column indices and inconsistent row pointers raise ValueError."""
        self.assertRaises(ValueError, TACS.BCSRMatrix, 2, [0, 2], [1, 0])
        self.assertRaises(ValueError, TACS.BCSRMatrix, 2, [0, 1], [3], 2)
        self.assertRaises(ValueError, TACS.BCSRMatrix, 2, [0, 2], [0])
        self.assertRaises(
            ValueError, TACS.BCSRMatrix, 2, [0, 1], [0], -1, np.zeros((1, 3, 3))
        )

    def test_incompatible_matrices_raise(self):
        """Assert that matrix operands with mismatched dimensions or block sizes raise ValueError."""
        rng = np.random.default_rng(SEED)
        square, _ = makeMatrix(rng, 2, makeSquareMask(rng, 5), isDominant=True)
        rect, _ = makeMatrix(rng, 2, makeMask(rng, 4, 6))
        otherBlock, _ = makeMatrix(rng, 3, makeSquareMask(rng, 5))
        missingDiag = TACS.BCSRMatrix(2, [0, 1, 2], [1, 0])
        calls = (
            lambda: square.matMultAdd(1.0, rect, rect),
            lambda: square.matMultAdd(1.0, otherBlock, otherBlock),
            lambda: square.applyLowerFactor(rect),
            lambda: square.applyUpperFactor(rect),
            lambda: square.axpy(1.0, rect),
            lambda: square.copyValues(otherBlock),
            lambda: rect.factor(),
            lambda: missingDiag.factor(),
            lambda: missingDiag.factorDiag(),
            lambda: TACS.BCSRMatrix.createILU(rect),
            lambda: square.partition(6),
        )
        for index, call in enumerate(calls):
            with self.subTest(call=index):
                self.assertRaises(ValueError, call)


if __name__ == "__main__":
    unittest.main()
