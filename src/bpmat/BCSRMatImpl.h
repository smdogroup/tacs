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

#ifndef TACS_BCSR_MAT_IMPLEMENTATION_H
#define TACS_BCSR_MAT_IMPLEMENTATION_H

/*
  Implementations of the various block-specific operations
*/
#include <pthread.h>

#include <math.h>  // fabs

#include "TACSObject.h"

class BCSRMatData : public TACSObject {
 public:
  BCSRMatData(int _bsize, int _nrows, int _ncols);
  ~BCSRMatData();

  // After initialization these values should never change!
  int nrows;  // The number of rows
  int ncols;  // The number of columns
  int bsize;  // The matrix block-size
  int *diag;  // A pointer to the diagonal entries
  int *rowp;  // A pointer into cols for the start of each row
  int *cols;  // The column number of each variable

  // Information about the size of groups to use
  int matvec_group_size;  // The size of groups for mat-vec operations
  int matmat_group_size;  // The size of groups for mat-mat operations

  // The storage space for each block - this can change
  TacsScalar *A;  // The vector of elements of each block
};

class BCSRMatThread : public TACSObject {
 public:
  BCSRMatThread(BCSRMatData *_mat);
  ~BCSRMatThread();

  // The scheduling functions
  // ------------------------
  // Matrix-matrix and matrix-vector product scheduler
  void init_mat_mult_sched();
  void mat_mult_sched_job(const int group_size, int *row);
  void mat_mult_sched_job_size(const int group_size, int *row, const int nrows);

  // Scheduler for L^{-1} applications
  void init_apply_lower_sched();
  void apply_lower_sched_job(const int group_size, int *index, int *irow,
                             int *jstart, int *jend);
  void apply_lower_mark_completed(const int group_size, int index, int irow,
                                  int jstart, int jend);

  // Scheduler for U^{-1} appilcations
  void init_apply_upper_sched();
  void apply_upper_sched_job(const int group_size, int *index, int *irow,
                             int *jstart, int *jend);
  void apply_upper_mark_completed(const int group_size, int index, int irow,
                                  int jstart, int jend);

  pthread_t threads[TACSThreadInfo::TACS_MAX_NUM_THREADS];

  // The input/output when dealing with vectors
  TacsScalar *input, *output;

  // BCSRMatData to be operated on by the threaded-implementation
  BCSRMatData *mat;

  // Data required for matrix-matrix products
  double alpha;
  BCSRMatData *Amat;
  BCSRMatData *Bmat;

  // The number of completed rows in the implementation
  int num_completed_rows;

  // Compute the row index size
  int *assigned_row_index;   // The index of the fully assigned rows
  int *completed_row_index;  // The indices of the full assigned columns

  // The threaded implementation
  pthread_mutex_t mutex;
  pthread_cond_t cond;
};

#if defined(__GNUC__)
#define TACS_BCSR_INLINE inline __attribute__((always_inline))
#else
#define TACS_BCSR_INLINE inline
#endif

/*!
  Compute the inverse of an n x n row-major matrix with partial pivoting.

  This is defined in the header and always inlined so that, when n is a
  compile-time constant, its loops have constant trip counts.

  Ainv == the inverse (output)
  A == the matrix, overwritten by its LU factorization
  ipiv == work array for the n pivot indices
  n == the matrix dimension

  Returns 0 on success, or k + 1 if the k-th pivot is zero.
*/
TACS_BCSR_INLINE int BMatComputeInverse(TacsScalar *Ainv, TacsScalar *A, int *ipiv, int n) {
  int fail = 0;

  for (int k = 0; k < n - 1; k++) {
    int nk = n * k;

    // Find the maximum value and use it as the pivot
    int r = k;
    double maxv = fabs(TacsRealPart(A[nk + k]));
    for (int j = k + 1; j < n; j++) {
      double t = fabs(TacsRealPart(A[n * j + k]));
      if (t > maxv) {
        maxv = t;
        r = j;
      }
    }

    ipiv[k] = r;

    // If a swap is required, swap the rows
    if (r != k) {
      int nr = n * r;
      for (int j = 0; j < n; j++) {
        TacsScalar t = A[nk + j];
        A[nk + j] = A[nr + j];
        A[nr + j] = t;
      }
    }

    if (A[nk + k] == 0.0) {
      fail = k + 1;
      return fail;
    }

    for (int i = k + 1; i < n; i++) {
      A[n * i + k] = A[n * i + k] / A[nk + k];
    }

    for (int i = k + 1; i < n; i++) {
      int ni = n * i;
      for (int j = k + 1; j < n; j++) {
        A[ni + j] -= A[ni + k] * A[nk + j];
      }
    }
  }

  // Now, compute the matrix-inverse
  for (int k = 0; k < n; k++) {
    int ip = k;
    for (int i = 0; i < n - 1; i++) {
      if (ip == ipiv[i]) {
        ip = i;
      } else if (ip == i) {
        ip = ipiv[i];
      }
    }

    for (int i = 0; i < ip; i++) {
      Ainv[n * i + k] = 0.0;
    }

    Ainv[n * ip + k] = 1.0;

    for (int i = ip + 1; i < n; i++) {
      int ni = n * i;
      Ainv[ni + k] = 0.0;
      for (int j = ip; j < i; j++) {
        Ainv[ni + k] -= A[ni + j] * Ainv[n * j + k];
      }
    }

    for (int i = n - 1; i >= 0; i--) {
      int ni = n * i;
      for (int j = i + 1; j < n; j++) {
        Ainv[ni + k] -= A[ni + j] * Ainv[n * j + k];
      }
      Ainv[ni + k] = Ainv[ni + k] / A[ni + i];
    }
  }

  return fail;
}

/*
  The generic implementation. These will be slow.
  Better to implement faster versions for a block size if you need it!
*/

void BCSRMatVecMult(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                       TacsScalar *z);
void BCSRMatVecMultTranspose(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultTransposeAdd(BCSRMatData *A, TacsScalar *inVec,
                                TacsScalar *addVec, TacsScalar *outVec);
void BCSRMatFactor(BCSRMatData *A);
void BCSRMatApplyLower(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd(double alpha, BCSRMatData *A, BCSRMatData *B,
                       BCSRMatData *C);
void BCSRMatFactorLower(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper(BCSRMatData *A, BCSRMatData *F);
void BCSRMatMatMultNormal(BCSRMatData *A, TacsScalar *s, BCSRMatData *B);
void BCSRMatApplySOR(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                     const int end, const int var_offset,
                     const TacsScalar *Adiag, const TacsScalar omega,
                     const TacsScalar *b, const TacsScalar *xext,
                     TacsScalar *x);

// ==============================================================================
// Block size templated implementations
// ==============================================================================
/**
 * @brief Simple transpose matrix-vector multiplication for small square
 * matrices
 *
 * Computes y += A^T * x, where A is a square matrix of size N x N.
 *
 * @tparam N square matrix size
 * @param A Matrix entries in row-major order
 * @param x Vector to multiply
 * @param y Vector to store result
 */
template <int N>
void addTransposeMatVec(const TacsScalar *const A, const TacsScalar *const x,
                        TacsScalar *const y) {
  for (int ii = 0; ii < N; ++ii) {
    for (int jj = 0; jj < N; ++jj) {
      y[ii] += A[jj * N + ii] * x[jj];
    }
  }
}

/**
 * @brief Compute the matrix-vector product: outVec += A^{T} * inVec
 *
 * @tparam blockSize Size of the blocks in the BCSR matrix
 * @param data BCSR matrix data
 * @param inVec Input vector
 * @param outVec Output vector
 */
template <int blockSize>
void BCSRBlockMatVecMultTranspose(BCSRMatData *data, TacsScalar *inVec,
                                  TacsScalar *outVec) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int b2 = blockSize * blockSize;

  const TacsScalar *a = data->A;
  for (int rowInd = 0; rowInd < nrows; rowInd++) {
    const int rowDataStart = rowp[rowInd];
    const int rowDataEnd = rowp[rowInd + 1];
    for (int k = rowDataStart; k < rowDataEnd; k++) {
      const int outRowInd = blockSize * cols[k];

      addTransposeMatVec<blockSize>(a, inVec, &outVec[outRowInd]);
      a += b2;
    }

    inVec += blockSize;
  }
}

/**
 * @brief Compute the matrix-vector product plus addition: outVec = A^{T} *
 * inVec + addVec
 *
 * @tparam blockSize Size of the blocks in the BCSR matrix
 * @param data BCSR matrix data
 * @param inVec Input vector
 * @param addVec Addition vector
 * @param outVec Output vector
 */
template <int blockSize>
void BCSRBlockMatVecMultTransposeAdd(BCSRMatData *data, TacsScalar *inVec,
                                     TacsScalar *addVec, TacsScalar *outVec) {
  // Copy addVec to outVec if they are different vectors
  const int ncols = data->ncols;
  if (outVec != addVec) {
    memcpy(outVec, addVec, blockSize * ncols * sizeof(TacsScalar));
  }
  BCSRBlockMatVecMultTranspose<blockSize>(data, inVec, outVec);
}

/*
  The block-specific code for block sizes 1 to 8 is in BCSRMatBlockImpl.h,
  as templates on the block size.
*/

#endif
