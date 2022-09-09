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

#include "BCSRMatImpl.h"
#include "tacslapack.h"

/*!
  Code for the general block-size case
*/

/*
  Implementation of the BCSRMatData class
*/
BCSRMatData::BCSRMatData(int _bsize, int _nrows, int _ncols) {
  bsize = _bsize;
  nrows = _nrows;
  ncols = _ncols;

  diag = NULL;
  rowp = NULL;
  cols = NULL;
  A = NULL;

  // The sizes of the groups of procs
  matvec_group_size = 1;
  matmat_group_size = 1;
}

BCSRMatData::~BCSRMatData() {
  if (diag) {
    delete[] diag;
  }
  if (rowp) {
    delete[] rowp;
  }
  if (cols) {
    delete[] cols;
  }
  if (A) {
    delete[] A;
  }
}

/*
  The implementation of the BCSRMatThread class
*/
BCSRMatThread::BCSRMatThread(BCSRMatData *_mat) {
  mat = _mat;
  mat->incref();

  Amat = NULL;
  Bmat = NULL;
  input = NULL;
  output = NULL;

  int nrows = mat->nrows;

  assigned_row_index = new int[nrows];
  completed_row_index = new int[nrows];

  num_completed_rows = 0;
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init(&cond, NULL);
}

BCSRMatThread::~BCSRMatThread() {
  mat->decref();

  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond);

  if (assigned_row_index) {
    delete[] assigned_row_index;
  }
  if (completed_row_index) {
    delete[] completed_row_index;
  }
}

/*
  Schedule the jobs for matrix-vector and matrix-matrix products with
  a given group size
*/
void BCSRMatThread::init_mat_mult_sched() { num_completed_rows = 0; }

void BCSRMatThread::mat_mult_sched_job(const int group_size, int *row) {
  const int nrows = mat->nrows;

  pthread_mutex_lock(&mutex);

  if (num_completed_rows < nrows) {
    *row = num_completed_rows;
    num_completed_rows += group_size;
  } else {
    *row = -1;
  }

  pthread_mutex_unlock(&mutex);
}

void BCSRMatThread::mat_mult_sched_job_size(const int group_size, int *row,
                                            const int nrows) {
  pthread_mutex_lock(&mutex);

  if (num_completed_rows < nrows) {
    *row = num_completed_rows;
    num_completed_rows += group_size;
  } else {
    *row = -1;
  }

  pthread_mutex_unlock(&mutex);
}

/*
  Initialize the data for the lower scheduler
*/
void BCSRMatThread::init_apply_lower_sched() {
  num_completed_rows = 0;
  memset(assigned_row_index, 0, mat->nrows * sizeof(int));
  memset(completed_row_index, 0, mat->nrows * sizeof(int));
}

/*
  Determine the next job to do based on what rows have been completed
  and what tasks are yet to be performed.

  Note that a row cannot be assigned unless all processes are
  completed on that row.
*/
void BCSRMatThread::apply_lower_sched_job(const int group_size, int *index,
                                          int *irow, int *jstart, int *jend) {
  // Obtain the mutex lock now
  pthread_mutex_lock(&mutex);

  *index = -1;
  *irow = -1;
  *jstart = -1;
  *jend = -1;

  const int nrows = mat->nrows;

  while (num_completed_rows < nrows && *index < 0) {
    int ii = num_completed_rows / group_size;
    int i = num_completed_rows;

    for (; i < nrows; i += group_size, ii++) {
      if ((completed_row_index[ii] == assigned_row_index[ii]) &&
          (completed_row_index[ii] == i)) {
        // A digaonal is available - assign it to a thread for completion
        *index = ii;
        *irow = i;
        *jstart = completed_row_index[ii];
        *jend = completed_row_index[ii] + group_size;

        assigned_row_index[ii] = completed_row_index[ii] + group_size;
        break;
      } else if ((completed_row_index[ii] < i + group_size) &&
                 (completed_row_index[ii] < num_completed_rows) &&
                 (completed_row_index[ii] == assigned_row_index[ii])) {
        // Assign the remainder of the avaiable row for completion by a thread
        *index = ii;
        *irow = i;
        *jstart = completed_row_index[ii];
        *jend = num_completed_rows;

        assigned_row_index[ii] = num_completed_rows;
        break;
      }
    }

    if (num_completed_rows < nrows && *index < 0) {
      pthread_cond_wait(&cond, &mutex);
    }
  }

  pthread_mutex_unlock(&mutex);
}

/*
  Once the operations on a row are completed, call this function
  to mark that the operation is completed final values
*/
void BCSRMatThread::apply_lower_mark_completed(const int group_size, int index,
                                               int irow, int jstart, int jend) {
  // Obtain the mutex lock now
  pthread_mutex_lock(&mutex);

  if (irow >= 0) {
    completed_row_index[index] = jend;
    if (jend > num_completed_rows) {
      num_completed_rows = jend;
    }
  }

  pthread_cond_broadcast(&cond);
  pthread_mutex_unlock(&mutex);
}

/*
  Initialize the data for the upper scheduler
*/
void BCSRMatThread::init_apply_upper_sched() {
  num_completed_rows = 0;
  memset(assigned_row_index, 0, mat->nrows * sizeof(int));
  memset(completed_row_index, 0, mat->nrows * sizeof(int));
}

/*
  Determine the next job to do based on what rows have been completed
  and what tasks are yet to be performed.

  Note that a row cannot be assigned unless all processes are
  completed on that row.
*/
void BCSRMatThread::apply_upper_sched_job(const int group_size, int *index,
                                          int *irow, int *jstart, int *jend) {
  // Obtain the mutex lock now
  pthread_mutex_lock(&mutex);

  *index = -1;
  *irow = -1;
  *jstart = -1;
  *jend = -1;

  const int nrows = mat->nrows;

  while (num_completed_rows < nrows && *index < 0) {
    int ii = num_completed_rows / group_size;
    int i = num_completed_rows;

    for (; i < nrows; i += group_size, ii++) {
      if ((completed_row_index[ii] == assigned_row_index[ii]) &&
          (completed_row_index[ii] == i)) {
        // A digaonal is available - assign it to a thread for completion
        *index = ii;
        *irow = nrows - i;
        *jend = nrows - completed_row_index[ii];
        *jstart = nrows - (completed_row_index[ii] + group_size);

        assigned_row_index[ii] = completed_row_index[ii] + group_size;
        break;
      } else if ((completed_row_index[ii] < i + group_size) &&
                 (completed_row_index[ii] < num_completed_rows) &&
                 (completed_row_index[ii] == assigned_row_index[ii])) {
        // Assign the remainder of the avaiable row for completion by a thread
        *index = ii;
        *irow = nrows - i;
        *jstart = nrows - num_completed_rows;
        *jend = nrows - completed_row_index[ii];

        assigned_row_index[ii] = num_completed_rows;
        break;
      }
    }

    if (num_completed_rows < nrows && *index < 0) {
      pthread_cond_wait(&cond, &mutex);
    }
  }

  pthread_mutex_unlock(&mutex);
}

/*
  Once the operations on a row are completed, call this function
  to mark that the operation is completed final values
*/
void BCSRMatThread::apply_upper_mark_completed(const int group_size, int index,
                                               int irow, int jstart, int jend) {
  // Obtain the mutex lock now
  pthread_mutex_lock(&mutex);
  const int nrows = mat->nrows;

  if (irow >= 0) {
    completed_row_index[index] = nrows - jstart;

    if (nrows - jstart > num_completed_rows) {
      num_completed_rows = nrows - jstart;
    }
  }

  pthread_cond_broadcast(&cond);
  pthread_mutex_unlock(&mutex);
}

/*!
  Compute the inverse of a matrix.

  A == A row-major ordered matrix of (bsize)x(bsize)
  w == A work array of size (bsize*bsize)
*/
int BMatComputeInverse(TacsScalar *Ainv, TacsScalar *A, int *ipiv, int n) {
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

/*!
  Compute the matrix-vector product: y = A * x
*/

void BCSRMatVecMult(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  int b2 = bsize * bsize;

  int one = 1;
  TacsScalar alpha = 1.0;
  TacsScalar beta = 1.0;

  TacsScalar *a = data->A;
  for (int i = 0; i < nrows; i++) {
    memset(y, 0, bsize * sizeof(TacsScalar));

    int end = rowp[i + 1];
    for (int k = rowp[i]; k < end; k++) {
      int bj = bsize * cols[k];

      BLASgemv("T", &bsize, &bsize, &alpha, a, &bsize, &x[bj], &one, &beta, y,
               &one);
      a += b2;
    }

    y += bsize;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/

void BCSRMatVecMultAdd(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                       TacsScalar *z) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  int one = 1;
  TacsScalar alpha = 1.0;
  TacsScalar beta = 1.0;

  for (int i = 0; i < nrows; i++) {
    if (z != y) {
      memcpy(z, y, bsize * sizeof(TacsScalar));
    }

    int end = rowp[i + 1];
    int j = rowp[i];
    TacsScalar *a = &data->A[b2 * j];
    for (; j < end; j++) {
      int bj = bsize * cols[j];

      BLASgemv("T", &bsize, &bsize, &alpha, a, &bsize, &x[bj], &one, &beta, y,
               &one);
      a += b2;
    }

    y += bsize;
    z += bsize;
  }
}

/*!
  Compute the matrix-vector product: y = A^{T} * x
*/
void BCSRMatVecMultTranspose(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  int one = 1;
  TacsScalar alpha = 1.0;
  TacsScalar beta = 1.0;

  TacsScalar *a = data->A;
  for (int i = 0; i < nrows; i++) {
    int end = rowp[i + 1];
    for (int k = rowp[i]; k < end; k++) {
      int bj = bsize * cols[k];

      BLASgemv("N", &bsize, &bsize, &alpha, a, &bsize, x, &one, &beta, &y[bj],
               &one);
      a += b2;
    }

    x += bsize;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  int one = 1;
  TacsScalar alpha = -1.0;
  TacsScalar beta = 1.0;

  TacsScalar *yy = y;
  TacsScalar *xx = x;

  for (int i = 0; i < nrows; i++) {
    if (x != y) {
      memcpy(yy, xx, bsize * sizeof(TacsScalar));
    }

    int end = diag[i];
    int j = rowp[i];
    TacsScalar *a = &data->A[b2 * j];
    for (; j < end; j++) {
      int bj = bsize * cols[j];

      BLASgemv("T", &bsize, &bsize, &alpha, a, &bsize, &y[bj], &one, &beta, yy,
               &one);
      a += b2;
    }

    yy += bsize;
    xx += bsize;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  int one = 1;
  TacsScalar alpha = -1.0;
  TacsScalar beta = 1.0;
  TacsScalar zero = 0.0;

  TacsScalar *ty = new TacsScalar[bsize];
  TacsScalar *yy = &y[bsize * (nrows - 1)];
  x = &x[bsize * (nrows - 1)];

  for (int i = nrows - 1; i >= 0; i--) {
    memcpy(ty, x, bsize * sizeof(TacsScalar));

    int end = rowp[i + 1];
    int j = diag[i] + 1;
    TacsScalar *a = &data->A[b2 * j];
    TacsScalar *adiag = a;
    adiag -= b2;

    for (int j = diag[i] + 1; j < end; j++) {
      int bj = bsize * cols[j];

      BLASgemv("T", &bsize, &bsize, &alpha, a, &bsize, &y[bj], &one, &beta, ty,
               &one);
      a += b2;
    }

    BLASgemv("T", &bsize, &bsize, &beta, adiag, &bsize, ty, &one, &zero, yy,
             &one);

    x -= bsize;
    yy -= bsize;
  }

  delete[] ty;
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyPartialLower(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  int off = bsize * var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    int bi = bsize * i - off;
    int k = rowp[i];
    while (cols[k] < var_offset) {
      k++;
    }

    int end = diag[i];
    for (; k < end; k++) {
      int bj = bsize * cols[k] - off;
      TacsScalar *a = &data->A[b2 * k];

      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          x[bi + m] -= a[bm + n] * x[bj + n];
        }
      }
    }
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyPartialUpper(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *ty = new TacsScalar[bsize];
  int off = bsize * var_offset;

  for (int i = nrows - 1; i >= var_offset; i--) {
    int bi = bsize * i - off;
    for (int m = 0; m < bsize; m++) {
      ty[m] = x[bi + m];
    }

    int end = rowp[i + 1];
    for (int k = diag[i] + 1; k < end; k++) {
      int bj = bsize * cols[k] - off;
      TacsScalar *a = &data->A[b2 * k];

      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          ty[m] -= a[bm + n] * x[bj + n];
        }
      }
    }

    // Apply the inverse on the diagonal
    TacsScalar *adiag = &data->A[b2 * diag[i]];
    for (int m = 0; m < bsize; m++) {
      int bm = bsize * m;
      x[bi + m] = 0.0;
      for (int n = 0; n < bsize; n++) {
        x[bi + m] += adiag[bm + n] * ty[n];
      }
    }
  }

  delete[] ty;
}

/*!
  Special function for the approximate Schur preconditioner.

  Given the input vector x = f, and y

  Compute x = U_b^{-1} ( L_b^{-1} f - (L_b^{-1} E) y )

  The matrix is factorized into the following form:
  A     = [ L_b          0   ][ U_b  L_b^{-1} E ]
  .       [ F U_b^{-1}   L_s ][ 0    U_s        ]

  where the division is set by the variable var_offset.
*/
void BCSRMatApplyFactorSchur(BCSRMatData *data, TacsScalar *x, int var_offset) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *tx = new TacsScalar[bsize];

  // Compute x = U_b^{-1} ( x - (L_b^{-1} E) y )
  for (int i = var_offset - 1; i >= 0; i--) {
    int bi = bsize * i;

    for (int n = 0; n < bsize; n++) {
      tx[n] = x[bi + n];
    }

    int end = rowp[i + 1];
    for (int k = diag[i] + 1; k < end; k++) {
      // x[i] = x[i] - A[j] *x[ cols[j] ];
      int j = cols[k];
      int bj = bsize * j;

      TacsScalar *a = &data->A[b2 * k];
      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          tx[m] -= a[bm + n] * x[bj + n];
        }
      }
    }

    // Apply the inverse on the diagonal
    TacsScalar *adiag = &data->A[b2 * diag[i]];
    for (int m = 0; m < bsize; m++) {
      int bm = bsize * m;
      x[bi + m] = 0.0;
      for (int n = 0; n < bsize; n++) {
        x[bi + m] += adiag[bm + n] * tx[n];
      }
    }
  }

  delete[] tx;
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMatApplySOR(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                     const int end, const int var_offset,
                     const TacsScalar *Adiag, const TacsScalar omega,
                     const TacsScalar *b, const TacsScalar *xext,
                     TacsScalar *x) {
  const int *Arowp = Adata->rowp;
  const int *Acols = Adata->cols;

  const int *Browp = NULL;
  const int *Bcols = NULL;
  if (Bdata) {
    Browp = Bdata->rowp;
    Bcols = Bdata->cols;
  }

  int bsize = Adata->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *tx = new TacsScalar[bsize];

  if (start < end) {
    for (int i = start; i < end; i++) {
      int bi = bsize * i;

      // Copy the right-hand-side to the temporary vector
      // for this row
      for (int n = 0; n < bsize; n++) {
        tx[n] = b[bi + n];
      }

      // Set the pointer to the beginning of the current
      // row
      TacsScalar *a = &Adata->A[b2 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        int bj = bsize * j;

        if (i != j) {
          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * x[bj + n];
            }
          }
        }

        // Increment the block pointer by bsize^2
        a += b2;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[36 * Browp[row]];
        end = Browp[row + 1];

        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[bsize * j];

          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * y[n];
            }
          }

          // Increment the block pointer by bsize^2
          a += b2;
        }
      }

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      for (int n = 0; n < bsize; n++) {
        x[bi + n] = (1.0 - omega) * x[bi + n];
      }

      // Apply the diagonal inverse and add the result to
      // the matrix
      const TacsScalar *adiag = &Adiag[b2 * i];
      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          x[bi + m] += omega * adiag[bm + n] * tx[n];
        }
      }
    }
  } else {
    // Go through the matrix with the forward ordering
    for (int i = start; i < end; i++) {
      int bi = bsize * i;

      // Copy the right-hand-side to the temporary vector
      // for this row
      for (int n = 0; n < bsize; n++) {
        tx[n] = b[bi + n];
      }

      // Set the pointer to the beginning of the current
      // row
      TacsScalar *a = &Adata->A[b2 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        int bj = bsize * j;

        if (i != j) {
          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * x[bj + n];
            }
          }
        }

        // Increment the block pointer by bsize^2
        a += b2;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[36 * Browp[row]];
        end = Browp[row + 1];

        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[bsize * j];

          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * y[n];
            }
          }

          // Increment the block pointer by bsize^2
          a += b2;
        }
      }

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      for (int n = 0; n < bsize; n++) {
        x[bi + n] = (1.0 - omega) * x[bi + n];
      }

      // Apply the diagonal inverse and add the result to
      // the matrix
      const TacsScalar *adiag = &Adiag[b2 * i];
      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          x[bi + m] += omega * adiag[bm + n] * tx[n];
        }
      }
    }
  }

  delete[] tx;
}

/*!
  Apply a given number of steps of Symmetric-SOR to the system A*x = b.
*/
void BCSRMatApplySSOR(BCSRMatData *data, TacsScalar *Adiag, TacsScalar omega,
                      int iters, TacsScalar *b, TacsScalar *x) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *tx = new TacsScalar[bsize];

  for (int iter = 0; iter < iters; iter++) {
    // Apply a forward sweep
    for (int i = 0; i < nrows; i++) {
      int bi = bsize * i;

      // Copy over the right-hand-side to the temporary vector tx
      for (int n = 0; n < bsize; n++) {
        tx[n] = b[bi + n];
      }

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = rowp[i + 1];
      for (int k = rowp[i]; k < end; k++) {
        int j = cols[k];
        int bj = bsize * j;

        if (i != j) {
          TacsScalar *a = &data->A[b2 * k];
          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * x[bj + n];
            }
          }
        }
      }

      // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
      for (int n = 0; n < bsize; n++) {
        x[bi + n] = (1.0 - omega) * x[bi + n];
      }

      // Apply the inverse on the diagonal
      TacsScalar *adiag = &Adiag[b2 * i];
      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          x[bi + m] += omega * adiag[bm + n] * tx[n];
        }
      }
    }

    // Apply the backward sweep
    for (int i = nrows - 1; i >= 0; i--) {
      int bi = bsize * i;

      for (int n = 0; n < bsize; n++) {
        tx[n] = b[bi + n];
      }

      int end = rowp[i + 1];
      for (int k = rowp[i]; k < end; k++) {
        int j = cols[k];
        int bj = bsize * j;

        if (i != j) {
          TacsScalar *a = &data->A[b2 * k];
          for (int m = 0; m < bsize; m++) {
            int bm = bsize * m;
            for (int n = 0; n < bsize; n++) {
              tx[m] -= a[bm + n] * x[bj + n];
            }
          }
        }
      }

      // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
      for (int n = 0; n < bsize; n++) {
        x[bi + n] = (1.0 - omega) * x[bi + n];
      }

      // Apply the inverse on the diagonal
      TacsScalar *adiag = &Adiag[b2 * i];
      for (int m = 0; m < bsize; m++) {
        int bm = bsize * m;
        for (int n = 0; n < bsize; n++) {
          x[bi + m] += omega * adiag[bm + n] * tx[n];
        }
      }
    }
  }

  delete[] tx;
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMatApplySOR(BCSRMatData *data, TacsScalar *Adiag, const int *pairs,
                     int npairs, TacsScalar omega, int iters, TacsScalar *b,
                     TacsScalar *x) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *tx = new TacsScalar[2 * bsize];

  for (int iter = 0; iter < iters; iter++) {
    const TacsScalar *D = Adiag;

    for (int i = 0, p = 0; i < nrows; i++) {
      if (p < npairs && i == pairs[p]) {
        int bi = bsize * i;

        // Copy the right-hand-side to the temporary vector
        // for this row
        for (int n = 0; n < 2 * bsize; n++) {
          tx[n] = b[bi + n];
        }

        // Set the pointer to the beginning of the current
        // row
        const TacsScalar *a = &data->A[b2 * rowp[i]];

        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }

          // Increment the block pointer by bsize^2
          a += b2;
        }

        end = rowp[i + 2];
        for (int k = rowp[i + 1]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m + bsize] -= a[bm + n] * x[bj + n];
              }
            }
          }

          // Increment the block pointer by bsize^2
          a += b2;
        }

        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        for (int n = 0; n < 2 * bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the diagonal inverse and add the result to
        // the matrix
        for (int m = 0; m < 2 * bsize; m++) {
          int bm = 2 * bsize * m;
          for (int n = 0; n < 2 * bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }

        D += 4 * b2;
        i++;
        p++;
      } else {
        int bi = bsize * i;

        // Copy the right-hand-side to the temporary vector
        // for this row
        for (int n = 0; n < bsize; n++) {
          tx[n] = b[bi + n];
        }

        // Set the pointer to the beginning of the current
        // row
        TacsScalar *a = &data->A[b2 * rowp[i]];

        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (i != j) {
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }

          // Increment the block pointer by bsize^2
          a += b2;
        }

        // Compute the first term in the update:
        // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
        for (int n = 0; n < bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the diagonal inverse and add the result to
        // the matrix
        for (int m = 0; m < bsize; m++) {
          int bm = bsize * m;
          for (int n = 0; n < bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }

        D += b2;
      }
    }
  }

  delete[] tx;
}

/*!
  Apply a given number of steps of Symmetric-SOR to the system A*x = b.
*/
void BCSRMatApplySSOR(BCSRMatData *data, TacsScalar *Adiag, const int *pairs,
                      int npairs, TacsScalar omega, int iters, TacsScalar *b,
                      TacsScalar *x) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  int bsize = data->bsize;
  const int b2 = bsize * bsize;

  TacsScalar *tx = new TacsScalar[2 * bsize];

  for (int iter = 0; iter < iters; iter++) {
    const TacsScalar *D = Adiag;

    // Apply a forward sweep
    for (int i = 0, p = 0; i < nrows; i++) {
      if (p < npairs && i == pairs[p]) {
        int bi = bsize * i;

        // Copy over the right-hand-side to the temporary vector tx
        for (int n = 0; n < 2 * bsize; n++) {
          tx[n] = b[bi + n];
        }

        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        end = rowp[i + 2];
        for (int k = rowp[i + 1]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[bsize + m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
        for (int n = 0; n < 2 * bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the inverse on the diagonal
        for (int m = 0; m < 2 * bsize; m++) {
          int bm = 2 * bsize * m;
          for (int n = 0; n < 2 * bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }

        D += 4 * b2;
        i++;
        p++;
      } else {
        int bi = bsize * i;

        // Copy over the right-hand-side to the temporary vector tx
        for (int n = 0; n < bsize; n++) {
          tx[n] = b[bi + n];
        }

        // Scan through the row and compute the result:
        // tx <- b_i - A_{ij}*x_{j} for j != i
        int end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
        for (int n = 0; n < bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the inverse on the diagonal
        for (int m = 0; m < bsize; m++) {
          int bm = bsize * m;
          for (int n = 0; n < bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }
      }
    }

    // Apply the backward sweep
    for (int i = nrows - 1, p = npairs - 1; i >= 0; i--) {
      if (p >= 0 && i - 1 == pairs[p]) {
        int bi = bsize * (i - 1);

        for (int n = 0; n < 2 * bsize; n++) {
          tx[n] = b[bi + n];
        }

        int end = rowp[i];
        for (int k = rowp[i - 1]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i - 1 && j != i) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (j != i && j != i + 1) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[bsize + m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
        for (int n = 0; n < 2 * bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the inverse on the diagonal
        D -= 4 * b2;
        for (int m = 0; m < 2 * bsize; m++) {
          int bm = 2 * bsize * m;
          for (int n = 0; n < 2 * bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }

        i--;
        p--;
      } else {
        int bi = bsize * i;

        for (int n = 0; n < bsize; n++) {
          tx[n] = b[bi + n];
        }

        int end = rowp[i + 1];
        for (int k = rowp[i]; k < end; k++) {
          int j = cols[k];
          int bj = bsize * j;

          if (i != j) {
            TacsScalar *a = &data->A[b2 * k];
            for (int m = 0; m < bsize; m++) {
              int bm = bsize * m;
              for (int n = 0; n < bsize; n++) {
                tx[m] -= a[bm + n] * x[bj + n];
              }
            }
          }
        }

        // x[i] = (1.0 - omega)*x[i] + D^{-1} delta x
        for (int n = 0; n < bsize; n++) {
          x[bi + n] = (1.0 - omega) * x[bi + n];
        }

        // Apply the inverse on the diagonal
        D -= b2;
        for (int m = 0; m < bsize; m++) {
          int bm = bsize * m;
          for (int n = 0; n < bsize; n++) {
            x[bi + m] += omega * D[bm + n] * tx[n];
          }
        }
      }
    }
  }

  delete[] tx;
}

/*!
  Perform a matrix-matrix multiplication
*/
void BCSRMatMatMultAdd(double alpha, BCSRMatData *Adata, BCSRMatData *Bdata,
                       BCSRMatData *Cdata) {
  // Retrieve the data required from the matrix
  const int nrows_a = Adata->nrows;
  const int *arowp = Adata->rowp;
  const int *acols = Adata->cols;
  const TacsScalar *A = Adata->A;

  const int *browp = Bdata->rowp;
  const int *bcols = Bdata->cols;
  const TacsScalar *B = Bdata->A;

  // The matrix being written to
  const int *crowp = Cdata->rowp;
  const int *ccols = Cdata->cols;
  TacsScalar *C = Cdata->A;

  const int bsize = Adata->bsize;
  const int b2 = bsize * bsize;

  // C_{ik} = A_{ij} B_{jk}
  for (int i = 0; i < nrows_a; i++) {
    for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
      int j = acols[jp];
      const TacsScalar *a = &A[b2 * jp];

      int kp = browp[j];
      int kp_end = browp[j + 1];

      int cp = crowp[i];
      int cp_end = crowp[i + 1];

      for (; kp < kp_end; kp++) {
        while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
          cp++;
        }
        if (cp >= cp_end) {
          break;
        }

        if (bcols[kp] == ccols[cp]) {
          const TacsScalar *b = &B[b2 * kp];
          TacsScalar *c = &C[b2 * cp];

          // Compute the matrix-matrix multiplication
          for (int n = 0; n < bsize; n++) {
            int nb = n * bsize;
            for (int m = 0; m < bsize; m++) {
              for (int l = 0; l < bsize; l++) {
                c[nb + m] += alpha * a[nb + l] * b[l * bsize + m];
              }
            }
          }
        }
      }
    }
  }
}

/*!
  Compute the scaled normal equations:

  A = B^{T}*s*B
*/
void BCSRMatMatMultNormal(BCSRMatData *Adata, TacsScalar *scale,
                          BCSRMatData *Bdata) {
  // Retrieve the data required from the matrix
  const int nrows_a = Adata->nrows;
  const int *arowp = Adata->rowp;
  const int *acols = Adata->cols;
  TacsScalar *A = Adata->A;

  const int nrows_b = Bdata->nrows;
  const int *browp = Bdata->rowp;
  const int *bcols = Bdata->cols;
  const TacsScalar *B = Bdata->A;

  const int bsize = Adata->bsize;
  const int b2 = bsize * bsize;

  int *kptr = new int[nrows_b];
  memcpy(kptr, browp, nrows_b * sizeof(int));

  // A_{ij} = B_{ki}*s{k}*B_{kj}
  for (int i = 0; i < nrows_a; i++) {
    // Scan through column i of the matrix B_{*i}
    for (int k = 0; k < nrows_b; k++) {
      if ((kptr[k] < browp[k + 1]) && (bcols[kptr[k]] == i)) {
        const TacsScalar *bi = &B[b2 * kptr[k]];
        kptr[k]++;

        // Add the non-zero pattern to
        // for j = 1, n do
        //    A_{ij} += B_{ki}*s_{k}*B_{kj}
        int jpa = arowp[i];
        int jpa_end = arowp[i + 1];

        int jpb = browp[k];
        int jpb_end = browp[k + 1];

        // Locate j such that (k,j) in nz(B_{kj}) and (i,j) in nz(A_{ij})
        for (; jpa < jpa_end; jpa++) {
          while ((jpb < jpb_end) && (bcols[jpb] < acols[jpa])) {
            jpb++;
          }
          if (jpb >= jpb_end) {
            break;
          }

          if (acols[jpa] == bcols[jpb]) {
            const TacsScalar *bj = &B[b2 * jpb];
            const TacsScalar *s = &scale[bsize * k];
            TacsScalar *a = &A[b2 * jpa];

            // a_{nm} += s_{l}*b_{ln}*b_{lm}
            for (int n = 0, nb = 0; n < bsize; n++, nb += bsize) {
              for (int l = 0, lb = 0; l < bsize; l++, lb += bsize) {
                for (int m = 0; m < bsize; m++) {
                  a[nb + m] += s[l] * bi[lb + n] * bj[lb + m];
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
