/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_BCSR_MAT_BLOCK_IMPLEMENTATION_H
#define TACS_BCSR_MAT_BLOCK_IMPLEMENTATION_H

/*
  BCSRMat kernels, serial and threaded, written once for any block size.

  Every kernel is a template on the block size N, so that each loop
  over a block has a compile-time trip count. The kernels are built
  from the small block operations at the top of this file, which are
  written so that compilers can vectorize them or keep their operands
  in registers without any architecture-specific code:

  - The block-vector products accumulate into a local array that the
    caller stores once per block row. Accumulating through the output
    pointer instead forces a store and a reload for every term, because
    the compiler must assume that the output aliases the inputs.

  - In real mode and for block sizes of TACS_BCSR_MIN_AXPY_BLOCK_SIZE or
    more, the block-block products are "axpy" loops along contiguous
    rows of their operands, with __restrict pointers.

  Everywhere else (the block-vector products, the smaller blocks and
  complex mode) the products are summed in the same order as the
  hand-unrolled kernels that this file replaces: all the terms of an
  entry are summed, then added or subtracted.

  All blocks are N x N, stored row-major.
*/

#include <stdio.h>   // fprintf, stderr
#include <string.h>  // memcpy

#include <atomic>  // std::atomic

#include "BCSRMatImpl.h"  // BCSRMatData, BMatComputeInverse, TACS_BCSR_INLINE
#include "TACSObject.h"   // TacsScalar, TacsAddFlops, TACS_USE_COMPLEX

/*
  The smallest block size for which the block-block products use the
  "axpy" loops. Smaller blocks, and all blocks in complex mode, sum the
  terms of each entry first ("dot" form): the axpy form is what compilers
  vectorize well for the larger blocks, but it was several times slower
  than the dot form for block sizes 3 and 4 (GCC 15, aarch64), and the
  dot form keeps the complex results unchanged.
*/
#ifdef TACS_USE_COMPLEX
#define TACS_BCSR_MIN_AXPY_BLOCK_SIZE 1000000
#else
#define TACS_BCSR_MIN_AXPY_BLOCK_SIZE 5
#endif

/*
  Tuning parameters of the threaded kernels (used with OpenMP only).

  Work is measured in blocks times the cube of the block size. The
  factorization loop uses threads when the average block below the
  diagonal has a row of at least TACS_BCSR_MIN_BLOCK_WORK to combine with
  its own. Dynamically scheduled loops hand out their rows in
  TACS_BCSR_CHUNKS_PER_THREAD chunks per thread. The values were chosen
  for block size 6 with GCC 15 on aarch64.
*/
#define TACS_BCSR_MIN_BLOCK_WORK 500L
#define TACS_BCSR_CHUNKS_PER_THREAD 8

// ==============================================================================
// Block operations
// ==============================================================================

/**
 * @brief Compute one entry of a block-vector product
 *
 * @tparam N Block size
 * @param aRow One row of a block (N entries)
 * @param x Vector of N entries
 * @return The dot product of aRow and x, summed from the first term to the
 * last
 */
template <int N>
static TACS_BCSR_INLINE TacsScalar blockRowDot(const TacsScalar *aRow,
                                               const TacsScalar *x) {
  TacsScalar sum = aRow[0] * x[0];
  for (int j = 1; j < N; j++) {
    sum += aRow[j] * x[j];
  }
  return sum;
}

/**
 * @brief Compute y += a * x for one block
 *
 * @tparam N Block size
 * @param a Block entries
 * @param x Vector of N entries
 * @param y Local accumulator of N entries, which must not alias a or x
 */
template <int N>
static TACS_BCSR_INLINE void blockMultAdd(const TacsScalar *a,
                                          const TacsScalar *x,
                                          TacsScalar *__restrict y) {
  for (int i = 0; i < N; i++) {
    y[i] += blockRowDot<N>(&a[N * i], x);
  }
}

/**
 * @brief Compute y -= a * x for one block
 *
 * @tparam N Block size
 * @param a Block entries
 * @param x Vector of N entries
 * @param y Local accumulator of N entries, which must not alias a or x
 */
template <int N>
static TACS_BCSR_INLINE void blockMultSub(const TacsScalar *a,
                                          const TacsScalar *x,
                                          TacsScalar *__restrict y) {
  for (int i = 0; i < N; i++) {
    y[i] -= blockRowDot<N>(&a[N * i], x);
  }
}

/**
 * @brief Compute y = a * x for one block
 *
 * @tparam N Block size
 * @param a Block entries
 * @param x Vector of N entries, which must not alias y
 * @param y Vector of N entries
 */
template <int N>
static TACS_BCSR_INLINE void blockMult(const TacsScalar *a,
                                       const TacsScalar *__restrict x,
                                       TacsScalar *__restrict y) {
  for (int i = 0; i < N; i++) {
    y[i] = blockRowDot<N>(&a[N * i], x);
  }
}

/**
 * @brief Compute one entry of a block-block product
 *
 * @tparam N Block size
 * @param dRow One row of the left block (N entries)
 * @param bCol Pointer to the first entry of a column of the right block
 * @return The product of dRow and the column, summed from the first term to
 * the last
 */
template <int N>
static TACS_BCSR_INLINE TacsScalar blockRowColDot(const TacsScalar *dRow,
                                                  const TacsScalar *bCol) {
  TacsScalar sum = dRow[0] * bCol[0];
  for (int k = 1; k < N; k++) {
    sum += dRow[k] * bCol[N * k];
  }
  return sum;
}

/**
 * @brief Compute c = a * b for blocks
 *
 * @tparam N Block size
 * @param a Left block
 * @param b Right block
 * @param c Result, which must not alias a or b
 */
template <int N>
static TACS_BCSR_INLINE void blockMatMult(const TacsScalar *a,
                                          const TacsScalar *b,
                                          TacsScalar *__restrict c) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      c[N * i + j] = blockRowColDot<N>(&a[N * i], &b[j]);
    }
  }
}

/**
 * @brief Compute a -= d * b for blocks that do not overlap
 *
 * This is the inner operation of the factorization kernels and of the
 * Schur complement product.
 *
 * @tparam N Block size
 * @param d Left block
 * @param b Right block
 * @param a Block to update
 */
template <int N>
static TACS_BCSR_INLINE void blockMatMultSub(const TacsScalar *__restrict d,
                                             const TacsScalar *__restrict b,
                                             TacsScalar *__restrict a) {
  if (N >= TACS_BCSR_MIN_AXPY_BLOCK_SIZE) {
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < N; k++) {
        const TacsScalar dik = d[N * i + k];
        for (int j = 0; j < N; j++) {
          a[N * i + j] -= dik * b[N * k + j];
        }
      }
    }
  } else {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        a[N * i + j] -= blockRowColDot<N>(&d[N * i], &b[j]);
      }
    }
  }
}

/**
 * @brief Compute a += alpha * d * b for blocks that do not overlap
 *
 * @tparam N Block size
 * @param alpha Scalar multiple
 * @param d Left block
 * @param b Right block
 * @param a Block to update
 */
template <int N>
static TACS_BCSR_INLINE void blockMatMultAddScaled(
    const double alpha, const TacsScalar *__restrict d,
    const TacsScalar *__restrict b, TacsScalar *__restrict a) {
  if (N >= TACS_BCSR_MIN_AXPY_BLOCK_SIZE) {
    for (int i = 0; i < N; i++) {
      for (int k = 0; k < N; k++) {
        const TacsScalar dik = alpha * d[N * i + k];
        for (int j = 0; j < N; j++) {
          a[N * i + j] += dik * b[N * k + j];
        }
      }
    }
  } else {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        a[N * i + j] += alpha * blockRowColDot<N>(&d[N * i], &b[j]);
      }
    }
  }
}

/**
 * @brief Update the blocks of a target row that match a source row
 *
 * For every block p in [p, pEnd) of the source row whose column also
 * appears in [k, kEnd) of the target row, compute
 * target[k] -= d * source[p]. Both rows have sorted columns.
 *
 * @tparam N Block size
 * @param d Multiplier block, which must not overlap the target blocks
 * @param sourceCols Column indices of the source matrix
 * @param source Block entries of the source matrix
 * @param p First block of the source row
 * @param pEnd One past the last block of the source row
 * @param targetCols Column indices of the target matrix
 * @param target Block entries of the target matrix
 * @param k First block of the target row
 * @param kEnd One past the last block of the target row
 * @return The number of blocks updated
 */
template <int N>
static TACS_BCSR_INLINE int blockRowUpdate(
    const TacsScalar *d, const int *sourceCols, const TacsScalar *source, int p,
    const int pEnd, const int *targetCols, TacsScalar *target, int k,
    const int kEnd) {
  int numUpdates = 0;
  for (; (p < pEnd) && (k < kEnd); p++) {
    // Determine where the two rows have the same elements
    while (k < kEnd && targetCols[k] < sourceCols[p]) {
      k++;
    }

    if (k < kEnd && targetCols[k] == sourceCols[p]) {
      blockMatMultSub<N>(d, &source[N * N * p], &target[N * N * k]);
      numUpdates++;
    }
  }
  return numUpdates;
}

// ==============================================================================
// Matrix-vector products
// ==============================================================================

/**
 * @brief Compute the matrix-vector product: y = A * x
 *
 * @tparam N Block size
 * @param data BCSR matrix data
 * @param x Input vector
 * @param y Output vector, which must not alias x
 */
template <int N>
void BCSRBlockMatVecMult(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;
  const int numThreads = data->num_threads;

#pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int i = 0; i < nrows; i++) {
    TacsScalar sums[N];
    for (int n = 0; n < N; n++) {
      sums[n] = 0.0;
    }

    const TacsScalar *a = &A[N * N * rowp[i]];
    const int end = rowp[i + 1];
    for (int k = rowp[i]; k < end; k++) {
      blockMultAdd<N>(a, &x[N * cols[k]], sums);
      a += N * N;
    }

    for (int n = 0; n < N; n++) {
      y[N * i + n] = sums[n];
    }
  }

  TacsAddFlops(2 * N * N * rowp[nrows]);
}

/**
 * @brief Compute the matrix-vector product plus addition: y = A * x + z
 *
 * This is the convention of BCSRMat::multAdd, which passes its output
 * vector as y.
 *
 * @tparam N Block size
 * @param data BCSR matrix data
 * @param x Input vector, which must not alias y
 * @param y Output vector
 * @param z Vector to add, which may be the same array as y
 */
template <int N>
void BCSRBlockMatVecMultAdd(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                            TacsScalar *z) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;
  const int numThreads = data->num_threads;

#pragma omp parallel for num_threads(numThreads) schedule(static)
  for (int i = 0; i < nrows; i++) {
    TacsScalar sums[N];
    for (int n = 0; n < N; n++) {
      sums[n] = z[N * i + n];
    }

    const TacsScalar *a = &A[N * N * rowp[i]];
    const int end = rowp[i + 1];
    for (int k = rowp[i]; k < end; k++) {
      blockMultAdd<N>(a, &x[N * cols[k]], sums);
      a += N * N;
    }

    for (int n = 0; n < N; n++) {
      y[N * i + n] = sums[n];
    }
  }

  TacsAddFlops(2 * N * N * rowp[nrows]);
}

// ==============================================================================
// Loop over the blocks below the diagonal
// ==============================================================================

/**
 * @brief Decide whether the loop over the blocks below the diagonal should
 * use its threaded form
 *
 * The threaded form costs some bookkeeping for every block below the
 * diagonal, which only pays off when these blocks carry enough work on
 * average. The work of a block is measured by the size of the row of its
 * column in workData.
 *
 * @param data BCSR matrix data with its diagonal index
 * @param workData Matrix whose rows measure the work
 * @return Whether to use the threaded loop
 */
static inline bool BCSRIsLowerLoopThreaded(const BCSRMatData *data,
                                           const BCSRMatData *workData) {
  if (data->num_threads <= 1) {
    return false;
  }

  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const int *workRowp = workData->rowp;

  long numBlocks = 0;
  long workBlocks = 0;
  for (int i = 0; i < data->nrows; i++) {
    numBlocks += diag[i] - rowp[i];
    for (int j = rowp[i]; j < diag[i]; j++) {
      workBlocks += workRowp[cols[j] + 1] - workRowp[cols[j]];
    }
  }

  const long blockWork = data->bsize * data->bsize * data->bsize;
  return blockWork * workBlocks >= TACS_BCSR_MIN_BLOCK_WORK * numBlocks;
}

/**
 * @brief Process the blocks below the diagonal of every row on several
 * threads, in an order in which a block is processed after the row of its
 * column is finished
 *
 * This is the threaded loop of the factorization kernels. For each row i,
 * entryFunc(i, j) is called for its blocks j = rowp[i], ..., diag[i] - 1
 * in that order, and then finishFunc(i), after which row i is finished.
 * entryFunc(i, j) is only called once row cols[j] is finished. Each row
 * therefore sees the same sequence of calls as in the serial loop of the
 * kernels, and the results do not depend on the number of threads.
 *
 * The threads take whole rows in increasing order from a shared counter.
 * Before each block, a thread waits until the row of the block's column is
 * finished, so later rows eliminate their first blocks while earlier rows
 * are still being finished. A thread only waits for rows that were taken
 * before its own, and the earliest unfinished row never waits, so the loop
 * cannot deadlock. Waiting threads spin, which assumes no more threads than
 * cores.
 *
 * It is kept out of line, and the kernels call their per-block functions
 * directly in their serial loops: with the threaded loop inlined, or with
 * the serial loop going through the lambdas, the serial factorization was
 * 5-7% slower (GCC 15, aarch64).
 *
 * @param data BCSR matrix data with its diagonal index
 * @param entryFunc Function of the row and the index of a block in it;
 * it may write to row i only, and read finished rows
 * @param finishFunc Function of the row
 */
template <class EntryFunc, class FinishFunc>
static TACS_BCSR_NOINLINE void BCSRForEachLowerBlockThreaded(
    BCSRMatData *data, const EntryFunc &entryFunc,
    const FinishFunc &finishFunc) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  std::atomic<bool> *isFinished = new std::atomic<bool>[nrows];
  for (int i = 0; i < nrows; i++) {
    isFinished[i].store(false, std::memory_order_relaxed);
  }
  std::atomic<int> nextRow(0);

#pragma omp parallel num_threads(data->num_threads)
  for (int i = nextRow.fetch_add(1); i < nrows; i = nextRow.fetch_add(1)) {
    for (int j = rowp[i]; j < diag[i]; j++) {
      // The acquire load makes the finished row's values visible
      while (!isFinished[cols[j]].load(std::memory_order_acquire)) {
      }
      entryFunc(i, j);
    }
    finishFunc(i);
    isFinished[i].store(true, std::memory_order_release);
  }

  delete[] isFinished;
}

// ==============================================================================
// Triangular solves with the factored matrix
// ==============================================================================

/**
 * @brief Finish one block row of an upper-triangular solve
 *
 * Subtract the products of the strictly upper blocks of row i with x from
 * rhs, multiply by the stored inverse of the diagonal block, and store the
 * result in xRow.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param i Block row
 * @param x Solution vector, indexed as x[N * cols[k] - offset]
 * @param offset Offset subtracted from the index into x
 * @param rhs Right-hand side for this row (N entries)
 * @param xRow Output for this row (N entries); may be the same array as rhs
 */
template <int N>
static TACS_BCSR_INLINE void blockUpperSolveRow(
    const BCSRMatData *data, const int i, const TacsScalar *x, const int offset,
    const TacsScalar *rhs, TacsScalar *xRow) {
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  TacsScalar acc[N];
  for (int n = 0; n < N; n++) {
    acc[n] = rhs[n];
  }

  const int end = data->rowp[i + 1];
  for (int k = data->diag[i] + 1; k < end; k++) {
    blockMultSub<N>(&A[N * N * k], &x[N * cols[k] - offset], acc);
  }

  // Apply the inverse of the diagonal block
  blockMult<N>(&A[N * N * data->diag[i]], acc, xRow);

  TacsAddFlops(2 * N * N * (end - data->diag[i]));
}

/**
 * @brief Compute one block row of a lower-triangular solve
 *
 * Subtract from rhs the products with y of the strictly lower blocks of row
 * i, and store the result in yRow.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param i Block row
 * @param y Solution vector
 * @param rhs Right-hand side for this row (N entries)
 * @param yRow Output for this row (N entries); may be the same array as rhs
 */
template <int N>
static TACS_BCSR_INLINE void blockLowerSolveRow(const BCSRMatData *data,
                                                const int i,
                                                const TacsScalar *y,
                                                const TacsScalar *rhs,
                                                TacsScalar *yRow) {
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  TacsScalar acc[N];
  for (int n = 0; n < N; n++) {
    acc[n] = rhs[n];
  }

  const int end = data->diag[i];
  for (int k = data->rowp[i]; k < end; k++) {
    blockMultSub<N>(&A[N * N * k], &y[N * cols[k]], acc);
  }

  for (int n = 0; n < N; n++) {
    yRow[n] = acc[n];
  }

  TacsAddFlops(2 * N * N * (end - data->rowp[i]));
}

/**
 * @brief Apply the lower factorization: y = L^{-1} x
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param x Right-hand side
 * @param y Solution, which may be the same array as x
 */
template <int N>
void BCSRBlockMatApplyLower(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  for (int i = 0; i < nrows; i++) {
    blockLowerSolveRow<N>(data, i, y, &x[N * i], &y[N * i]);
  }
}

/**
 * @brief Apply the upper factorization: y = U^{-1} x
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param x Right-hand side
 * @param y Solution, which may be the same array as x
 */
template <int N>
void BCSRBlockMatApplyUpper(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  for (int i = data->nrows - 1; i >= 0; i--) {
    blockUpperSolveRow<N>(data, i, y, 0, &x[N * i], &y[N * i]);
  }
}

/**
 * @brief Apply a portion of the lower factorization: x = L^{-1} x
 *
 * Only the rows and columns from var_offset onwards take part; x holds the
 * entries of those rows only.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param x Vector for rows var_offset to nrows, modified in place
 * @param var_offset First row of the partial factorization
 */
template <int N>
void BCSRBlockMatApplyPartialLower(BCSRMatData *data, TacsScalar *x,
                                   int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;
  const int off = N * var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    TacsScalar *xRow = &x[N * i - off];
    TacsScalar acc[N];
    for (int n = 0; n < N; n++) {
      acc[n] = xRow[n];
    }

    const int end = diag[i];
    int k = rowp[i];
    while (k < end && cols[k] < var_offset) {
      k++;
    }

    for (; k < end; k++) {
      blockMultSub<N>(&A[N * N * k], &x[N * cols[k] - off], acc);
    }

    for (int n = 0; n < N; n++) {
      xRow[n] = acc[n];
    }
  }
}

/**
 * @brief Apply a portion of the upper factorization: x = U^{-1} x
 *
 * Only the rows and columns from var_offset onwards take part; x holds the
 * entries of those rows only.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param x Vector for rows var_offset to nrows, modified in place
 * @param var_offset First row of the partial factorization
 */
template <int N>
void BCSRBlockMatApplyPartialUpper(BCSRMatData *data, TacsScalar *x,
                                   int var_offset) {
  const int off = N * var_offset;

  for (int i = data->nrows - 1; i >= var_offset; i--) {
    TacsScalar *xRow = &x[N * i - off];
    blockUpperSolveRow<N>(data, i, x, off, xRow, xRow);
  }
}

/**
 * @brief Apply the upper factor for the approximate Schur preconditioner
 *
 * Given x = [f; y], compute, in the first var_offset rows of x,
 *
 * x = U_b^{-1} (L_b^{-1} f - (L_b^{-1} E) y)
 *
 * where the factored matrix has the form
 *
 * A = [ L_b          0   ][ U_b  L_b^{-1} E ]
 *     [ F U_b^{-1}   L_s ][ 0    U_s        ]
 *
 * and the first var_offset rows of x already hold L_b^{-1} f.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param x Vector over all the rows, modified in its first var_offset rows
 * @param var_offset Number of rows in the leading block
 */
template <int N>
void BCSRBlockMatApplyFactorSchur(BCSRMatData *data, TacsScalar *x,
                                  int var_offset) {
  for (int i = var_offset - 1; i >= 0; i--) {
    blockUpperSolveRow<N>(data, i, x, 0, &x[N * i], &x[N * i]);
  }
}

// ==============================================================================
// SOR
// ==============================================================================

/**
 * @brief Apply one SOR update to block row i of the system A * x = b
 *
 * @tparam N Block size
 * @param Adata BCSR matrix data
 * @param Bdata Optional off-processor matrix data (may be NULL)
 * @param i Block row to update
 * @param var_offset First row that has entries in Bdata
 * @param Adiag Inverses of the diagonal blocks of A
 * @param omega Relaxation factor
 * @param b Right-hand side
 * @param xext External variables multiplied by Bdata
 * @param x Solution, updated in row i
 */
template <int N>
static TACS_BCSR_INLINE void blockSORRow(
    const BCSRMatData *Adata, const BCSRMatData *Bdata, const int i,
    const int var_offset, const TacsScalar *Adiag, const TacsScalar omega,
    const TacsScalar *b, const TacsScalar *xext, TacsScalar *x) {
  // Copy the right-hand-side to the temporary vector for this row
  TacsScalar t[N];
  for (int n = 0; n < N; n++) {
    t[n] = b[N * i + n];
  }

  // Compute t <- b_i - A_{ij} * x_{j} for j != i
  const int *Acols = Adata->cols;
  const TacsScalar *a = &Adata->A[N * N * Adata->rowp[i]];
  const int end = Adata->rowp[i + 1];
  for (int k = Adata->rowp[i]; k < end; k++) {
    const int j = Acols[k];
    if (i != j) {
      blockMultSub<N>(a, &x[N * j], t);
    }
    a += N * N;
  }

  // Subtract the product with the external variables
  if (Bdata != NULL && i >= var_offset) {
    const int row = i - var_offset;
    const int *Bcols = Bdata->cols;
    a = &Bdata->A[N * N * Bdata->rowp[row]];
    const int bEnd = Bdata->rowp[row + 1];
    for (int k = Bdata->rowp[row]; k < bEnd; k++) {
      blockMultSub<N>(a, &xext[N * Bcols[k]], t);
      a += N * N;
    }
  }

  // x[i] = (1.0 - omega) * x[i] + omega * D^{-1} t
  const TacsScalar *d = &Adiag[N * N * i];
  for (int n = 0; n < N; n++) {
    x[N * i + n] =
        (1.0 - omega) * x[N * i + n] + omega * blockRowDot<N>(&d[N * n], t);
  }
}

/**
 * @brief Apply a sweep of SOR to the system A * x = b
 *
 * The rows are visited in the order start, start + 1, ..., end - 1 when
 * start < end, and start - 1, start - 2, ..., end otherwise.
 *
 * @tparam N Block size
 * @param Adata BCSR matrix data
 * @param Bdata Optional off-processor matrix data (may be NULL)
 * @param start Start of the row range
 * @param end End of the row range
 * @param var_offset First row that has entries in Bdata
 * @param Adiag Inverses of the diagonal blocks of A
 * @param omega Relaxation factor
 * @param b Right-hand side
 * @param xext External variables multiplied by Bdata
 * @param x Solution, updated in place
 */
template <int N>
void BCSRBlockMatApplySOR(BCSRMatData *Adata, BCSRMatData *Bdata,
                          const int start, const int end, const int var_offset,
                          const TacsScalar *Adiag, const TacsScalar omega,
                          const TacsScalar *b, const TacsScalar *xext,
                          TacsScalar *x) {
  if (start < end) {
    for (int i = start; i < end; i++) {
      blockSORRow<N>(Adata, Bdata, i, var_offset, Adiag, omega, b, xext, x);
    }
  } else {
    for (int i = start - 1; i >= end; i--) {
      blockSORRow<N>(Adata, Bdata, i, var_offset, Adiag, omega, b, xext, x);
    }
  }
}

// ==============================================================================
// Factorization
// ==============================================================================

/**
 * @brief Eliminate one block below the diagonal in the ILU factorization
 *
 * For the block at index j of row `row`, in column cj < row, form the
 * multiplier L(row, cj) = A(row, cj) * U(cj, cj)^{-1}, store it, and
 * subtract its products with the blocks of row cj past the diagonal from
 * row `row`.
 *
 * @tparam N Block size
 * @param data BCSR matrix data, whose row cj is already factored
 * @param row Block row to eliminate
 * @param j Index into cols of the block to eliminate
 */
template <int N>
static TACS_BCSR_INLINE void blockFactorEntry(BCSRMatData *data, const int row,
                                              const int j) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  TacsScalar *A = data->A;
  const int cj = cols[j];

  // Form the multiplier d = A[j] * A[diag[cj]] in a local, so that no
  // store to the matrix can alias it, and write it back to the row
  TacsScalar d[N * N];
  blockMatMult<N>(&A[N * N * j], &A[N * N * diag[cj]], d);
  memcpy(&A[N * N * j], d, sizeof(d));

  // A[k] = A[k] - d * A[p] for the blocks past the diagonal of row cj
  const int nz = blockRowUpdate<N>(d, cols, A, diag[cj] + 1, rowp[cj + 1], cols,
                                   A, j + 1, rowp[row + 1]);

  TacsAddFlops(2 * N * N * N * nz + (2 * N - 1) * N * N);
  (void)nz;
}

/**
 * @brief Replace the diagonal block of a row by its inverse
 *
 * @tparam N Block size
 * @param data BCSR matrix data
 * @param row Block row whose blocks below the diagonal are eliminated
 */
template <int N>
static TACS_BCSR_INLINE void blockFactorDiag(BCSRMatData *data, const int row) {
  TacsScalar D[N * N];
  TacsScalar *a = &data->A[N * N * data->diag[row]];
  memcpy(D, a, sizeof(D));

  int ipiv[N];
  const int info = BMatComputeInverse(a, D, ipiv, N);
  if (info > 0) {
    fprintf(stderr,
            "Error during factorization of diagonal %d in \
block row %d \n",
            row + 1, info);
  }
}

/**
 * @brief Perform an ILU factorization of the matrix in place
 *
 * The existing non-zero pattern is used. On exit, the strictly lower
 * blocks hold L (with unit diagonal blocks), the diagonal blocks hold the
 * inverses of the diagonal blocks of U, and the strictly upper blocks hold
 * U.
 *
 * @tparam N Block size
 * @param data BCSR matrix data, overwritten by its factorization
 */
template <int N>
void BCSRBlockMatFactor(BCSRMatData *data) {
  const int nrows = data->nrows;

  for (int i = 0; i < nrows; i++) {
    if (data->diag[i] < 0) {
      fprintf(stderr, "Error in factorization: no diagonal entry for row %d",
              i);
      return;
    }
  }

  if (BCSRIsLowerLoopThreaded(data, data)) {
    BCSRForEachLowerBlockThreaded(
        data,
        [=](const int i, const int j) { blockFactorEntry<N>(data, i, j); },
        [=](const int i) { blockFactorDiag<N>(data, i); });
  } else {
    const int *rowp = data->rowp;
    const int *diag = data->diag;
    for (int i = 0; i < nrows; i++) {
      for (int j = rowp[i]; j < diag[i]; j++) {
        blockFactorEntry<N>(data, i, j);
      }
      blockFactorDiag<N>(data, i);
    }
  }

  // Add flops from the diagonal inversion
  TacsAddFlops(1.333333 * N * N * N * nrows);
}

/**
 * @brief Apply one block of L to a row in E <- L^{-1} E
 *
 * For the strictly lower block at index j of row `row` of the factored
 * matrix, in column cj, subtract its products with row cj of E from row
 * `row` of E.
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param Edata Matrix to update, whose row cj is already complete
 * @param row Block row to update
 * @param j Index into the cols of data of the block of L
 */
template <int N>
static TACS_BCSR_INLINE void blockFactorLowerEntry(const BCSRMatData *data,
                                                   BCSRMatData *Edata,
                                                   const int row, const int j) {
  const int *erowp = Edata->rowp;
  const int *ecols = Edata->cols;
  TacsScalar *E = Edata->A;
  const int cj = data->cols[j];

  // E[row, :] = E[row, :] - A[j] * E[cj, :]
  const int nz =
      blockRowUpdate<N>(&data->A[N * N * j], ecols, E, erowp[cj], erowp[cj + 1],
                        ecols, E, erowp[row], erowp[row + 1]);

  TacsAddFlops(2 * N * N * N * nz);
  (void)nz;
}

/**
 * @brief Compute E <- L^{-1} E, restricted to the non-zero pattern of E
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param Edata Matrix to update, with as many rows as data
 */
template <int N>
void BCSRBlockMatFactorLower(BCSRMatData *data, BCSRMatData *Edata) {
  if (BCSRIsLowerLoopThreaded(data, Edata)) {
    BCSRForEachLowerBlockThreaded(
        data,
        [=](const int i, const int j) {
          blockFactorLowerEntry<N>(data, Edata, i, j);
        },
        [](const int) {});
  } else {
    const int *rowp = data->rowp;
    const int *diag = data->diag;
    for (int i = 0; i < data->nrows; i++) {
      for (int j = rowp[i]; j < diag[i]; j++) {
        blockFactorLowerEntry<N>(data, Edata, i, j);
      }
    }
  }
}

/**
 * @brief Compute one block row of F <- F U^{-1}
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param Fdata Matrix to update
 * @param row Block row of F to update
 */
template <int N>
static TACS_BCSR_INLINE void blockFactorUpperRow(const BCSRMatData *data,
                                                 BCSRMatData *Fdata,
                                                 const int row) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  const int *fcols = Fdata->cols;
  TacsScalar *F = Fdata->A;

  const int rowEnd = Fdata->rowp[row + 1];
  for (int j = Fdata->rowp[row]; j < rowEnd; j++) {
    const int cj = fcols[j];

    // Form the multiplier d = F[j] * A[diag[cj]] in a local and write it
    // back to the row
    TacsScalar d[N * N];
    blockMatMult<N>(&F[N * N * j], &A[N * N * diag[cj]], d);
    memcpy(&F[N * N * j], d, sizeof(d));

    // F[k] = F[k] - d * A[p] for the blocks past the diagonal of row cj
    const int nz = blockRowUpdate<N>(d, cols, A, diag[cj] + 1, rowp[cj + 1],
                                     fcols, F, j + 1, rowEnd);

    TacsAddFlops(2 * N * N * N * nz + (2 * N - 1) * N * N);
    (void)nz;
  }
}

/**
 * @brief Compute F <- F U^{-1}, restricted to the non-zero pattern of F
 *
 * @tparam N Block size
 * @param data Factored BCSR matrix data
 * @param Fdata Matrix to update, with as many columns as data
 */
template <int N>
void BCSRBlockMatFactorUpper(BCSRMatData *data, BCSRMatData *Fdata) {
  const int nrows_f = Fdata->nrows;
  const int numThreads = data->num_threads;
  const int chunk = 1 + nrows_f / (TACS_BCSR_CHUNKS_PER_THREAD * numThreads);

#pragma omp parallel for num_threads(numThreads) schedule(dynamic, chunk)
  for (int i = 0; i < nrows_f; i++) {
    blockFactorUpperRow<N>(data, Fdata, i);
  }
}

// ==============================================================================
// Matrix-matrix products
// ==============================================================================

/**
 * @brief Compute one block row of C += alpha * A * B
 *
 * @tparam N Block size
 * @param alpha Scalar multiple
 * @param Adata Left matrix
 * @param Bdata Right matrix
 * @param Cdata Matrix to update, which must not be A or B
 * @param i Block row of C to update
 */
template <int N>
static TACS_BCSR_INLINE void blockMatMultAddRow(const double alpha,
                                                const BCSRMatData *Adata,
                                                const BCSRMatData *Bdata,
                                                BCSRMatData *Cdata,
                                                const int i) {
  const int *acols = Adata->cols;
  const TacsScalar *A = Adata->A;

  const int *browp = Bdata->rowp;
  const int *bcols = Bdata->cols;
  const TacsScalar *B = Bdata->A;

  const int *ccols = Cdata->cols;
  TacsScalar *C = Cdata->A;
  const int cStart = Cdata->rowp[i];
  const int cEnd = Cdata->rowp[i + 1];

  // C_{ik} += alpha * A_{ij} B_{jk}
  const int jpEnd = Adata->rowp[i + 1];
  for (int jp = Adata->rowp[i]; jp < jpEnd; jp++) {
    const int j = acols[jp];
    const TacsScalar *a = &A[N * N * jp];

    if (alpha == -1.0) {
      const int nz = blockRowUpdate<N>(a, bcols, B, browp[j], browp[j + 1],
                                       ccols, C, cStart, cEnd);
      TacsAddFlops(2 * N * N * N * nz);
      (void)nz;
    } else {
      int cp = cStart;
      const int kpEnd = browp[j + 1];
      for (int kp = browp[j]; (kp < kpEnd) && (cp < cEnd); kp++) {
        while ((cp < cEnd) && (ccols[cp] < bcols[kp])) {
          cp++;
        }

        if ((cp < cEnd) && (bcols[kp] == ccols[cp])) {
          blockMatMultAddScaled<N>(alpha, a, &B[N * N * kp], &C[N * N * cp]);
          TacsAddFlops(2 * N * N * N);
        }
      }
    }
  }
}

/**
 * @brief Compute C += alpha * A * B, restricted to the non-zero pattern of C
 *
 * @tparam N Block size
 * @param alpha Scalar multiple
 * @param Adata Left matrix
 * @param Bdata Right matrix
 * @param Cdata Matrix to update, which must not be A or B
 */
template <int N>
void BCSRBlockMatMatMultAdd(double alpha, BCSRMatData *Adata,
                            BCSRMatData *Bdata, BCSRMatData *Cdata) {
  const int nrows_a = Adata->nrows;
  const int numThreads = Cdata->num_threads;
  const int chunk = 1 + nrows_a / (TACS_BCSR_CHUNKS_PER_THREAD * numThreads);

#pragma omp parallel for num_threads(numThreads) schedule(dynamic, chunk)
  for (int i = 0; i < nrows_a; i++) {
    blockMatMultAddRow<N>(alpha, Adata, Bdata, Cdata, i);
  }
}

/**
 * @brief Compute the scaled normal equations: A += B^{T} * S * B
 *
 * S is the diagonal matrix of the entries of scale, and the update is
 * restricted to the non-zero pattern of A.
 *
 * @tparam N Block size
 * @param Adata Matrix to update
 * @param scale Diagonal scaling, one entry per row of B
 * @param Bdata Matrix whose normal equations are formed
 */
template <int N>
void BCSRBlockMatMatMultNormal(BCSRMatData *Adata, TacsScalar *scale,
                               BCSRMatData *Bdata) {
  const int nrows_a = Adata->nrows;
  const int *arowp = Adata->rowp;
  const int *acols = Adata->cols;
  TacsScalar *A = Adata->A;

  const int nrows_b = Bdata->nrows;
  const int *browp = Bdata->rowp;
  const int *bcols = Bdata->cols;
  const TacsScalar *B = Bdata->A;

  // The next unvisited block in each row of B
  int *kptr = new int[nrows_b];
  memcpy(kptr, browp, nrows_b * sizeof(int));

  // A_{ij} += B_{ki}^{T} * s_{k} * B_{kj}
  for (int i = 0; i < nrows_a; i++) {
    // Scan through column i of the matrix B
    for (int k = 0; k < nrows_b; k++) {
      if ((kptr[k] < browp[k + 1]) && (bcols[kptr[k]] == i)) {
        const TacsScalar *bi = &B[N * N * kptr[k]];
        const TacsScalar *s = &scale[N * k];
        kptr[k]++;

        // Locate j such that (k, j) is in B and (i, j) is in A
        int jpb = browp[k];
        const int jpbEnd = browp[k + 1];
        const int jpaEnd = arowp[i + 1];
        for (int jpa = arowp[i]; (jpa < jpaEnd) && (jpb < jpbEnd); jpa++) {
          while ((jpb < jpbEnd) && (bcols[jpb] < acols[jpa])) {
            jpb++;
          }

          if ((jpb < jpbEnd) && (acols[jpa] == bcols[jpb])) {
            const TacsScalar *bj = &B[N * N * jpb];
            TacsScalar *a = &A[N * N * jpa];

            // a_{nm} += s_{l} * b_{ln} * b_{lm}
            for (int n = 0; n < N; n++) {
              for (int l = 0; l < N; l++) {
                const TacsScalar weight = s[l] * bi[N * l + n];
                for (int m = 0; m < N; m++) {
                  a[N * n + m] += weight * bj[N * l + m];
                }
              }
            }
          }
        }
      }
    }
  }

  delete[] kptr;
}

#endif  // TACS_BCSR_MAT_BLOCK_IMPLEMENTATION_H
