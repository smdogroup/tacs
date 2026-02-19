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

/*
  Code for the 2x2 block case
*/

/*!
  Compute the matrix-vector product: y = A * x
*/

void BCSRMatVecMult2(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  for (int i = 0; i < nrows; i++) {
    y[0] = 0.0;
    y[1] = 0.0;

    int end = rowp[i + 1];
    int k = rowp[i];
    const TacsScalar *a = &A[4 * k];
    for (; k < end; k++) {
      int bj = 2 * cols[k];

      y[0] += a[0] * x[bj] + a[1] * x[bj + 1];
      y[1] += a[2] * x[bj] + a[3] * x[bj + 1];
      a += 4;
    }
    y += 2;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/

void BCSRMatVecMultAdd2(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  for (int i = 0; i < nrows; i++) {
    z[0] = y[0];
    z[1] = y[1];

    int end = rowp[i + 1];
    int k = rowp[i];
    const TacsScalar *a = &A[4 * k];
    for (; k < end; k++) {
      int bj = 2 * cols[k];

      z[0] += a[0] * x[bj] + a[1] * x[bj + 1];
      z[1] += a[2] * x[bj] + a[3] * x[bj + 1];
      a += 4;
    }
    y += 2;
    z += 2;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/

void BCSRMatApplyLower2(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *yy = y;

  for (int i = 0; i < nrows; i++) {
    yy[0] = x[0];
    yy[1] = x[1];

    int end = diag[i];
    int j = rowp[i];
    const TacsScalar *a = &A[4 * j];
    for (; j < end; j++) {
      int bj = 2 * cols[j];

      yy[0] -= a[0] * y[bj] + a[1] * y[bj + 1];
      yy[1] -= a[2] * y[bj] + a[3] * y[bj + 1];
      a += 4;
    }
    x += 2;
    yy += 2;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/

void BCSRMatApplyUpper2(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y1, y2;
  TacsScalar *yy = &y[2 * (nrows - 1)];
  x = &x[2 * (nrows - 1)];

  for (int i = nrows - 1; i >= 0; i--) {
    y1 = x[0];
    y2 = x[1];

    int end = rowp[i + 1];
    int j = diag[i];
    const TacsScalar *a = &A[4 * j];
    TacsScalar a11, a12, a21, a22;
    a11 = a[0];
    a12 = a[1];
    a21 = a[2];
    a22 = a[3];
    j++;
    a += 4;

    for (; j < end; j++) {
      int bj = 2 * cols[j];

      y1 -= a[0] * y[bj] + a[1] * y[bj + 1];
      y2 -= a[2] * y[bj] + a[3] * y[bj + 1];
      a += 4;
    }

    // Apply the inverse on the diagonal
    yy[0] = a11 * y1 + a12 * y2;
    yy[1] = a21 * y1 + a22 * y2;

    yy -= 2;
    x -= 2;
  }
}

/*!
  Apply the lower factorization x = L^{-1} x
*/

void BCSRMatApplyPartialLower2(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *xx = &x[2];
  int off = 2 * var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    int end = diag[i];
    int j = rowp[i];
    while (cols[j] < var_offset) j++;

    const TacsScalar *a = &A[4 * j];
    for (; j < end; j++) {
      int bj = 2 * cols[j] - off;

      xx[0] -= a[0] * x[bj] + a[1] * x[bj + 1];
      xx[1] -= a[2] * x[bj] + a[3] * x[bj + 1];
      a += 4;
    }
    xx += 2;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/

void BCSRMatApplyPartialUpper2(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y1, y2;
  int off = 2 * var_offset;

  TacsScalar *xx = &x[2 * (nrows - var_offset - 1)];

  for (int i = nrows - 1; i >= var_offset; i--) {
    y1 = xx[0];
    y2 = xx[1];

    int end = rowp[i + 1];
    int j = diag[i];
    const TacsScalar *adiag = &A[4 * j];
    const TacsScalar *a = adiag;
    j++;
    a += 4;

    for (; j < end; j++) {
      int bj = 2 * cols[j] - off;

      y1 -= a[0] * x[bj] + a[1] * x[bj + 1];
      y2 -= a[2] * x[bj] + a[3] * x[bj + 1];
      a += 4;
    }

    xx[0] = adiag[0] * y1 + adiag[1] * y2;
    xx[1] = adiag[2] * y1 + adiag[3] * y2;
    xx -= 2;
  }
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

void BCSRMatApplyFactorSchur2(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y1, y2;
  TacsScalar *xx = &x[2 * (var_offset - 1)];

  for (int i = var_offset - 1; i >= 0; i--) {
    y1 = xx[0];
    y2 = xx[1];

    int j = diag[i];
    int end = rowp[i + 1];
    const TacsScalar *adiag = &A[4 * j];
    const TacsScalar *a = adiag;
    a += 4;
    j++;

    for (; j < end; j++) {
      int bj = 2 * cols[j];
      y1 -= a[0] * x[bj] + a[1] * x[bj + 1];
      y2 -= a[2] * x[bj] + a[3] * x[bj + 1];
      a += 4;
    }

    // Apply the inverse on the diagonal
    xx[0] = adiag[0] * y1 + adiag[1] * y2;
    xx[1] = adiag[2] * y1 + adiag[3] * y2;
    xx -= 2;
  }
}

/*!
  Perform a matrix-matrix multiplication
*/

void BCSRMatMatMultAdd2(double alpha, BCSRMatData *Adata, BCSRMatData *Bdata,
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

  if (alpha == 1.0) {
    // C_{ik} = A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[4 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[4 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[4 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 4;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] += a[0] * b[0] + a[1] * b[2];
            c[1] += a[0] * b[1] + a[1] * b[3];

            c[2] += a[2] * b[0] + a[3] * b[2];
            c[3] += a[2] * b[1] + a[3] * b[3];
          }
          b += 4;
        }
      }
    }
  } else if (alpha == -1.0) {
    // C_{ik} = A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[4 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[4 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[4 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 4;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] -= a[0] * b[0] + a[1] * b[2];
            c[1] -= a[0] * b[1] + a[1] * b[3];

            c[2] -= a[2] * b[0] + a[3] * b[2];
            c[3] -= a[2] * b[1] + a[3] * b[3];
          }
          b += 4;
        }
      }
    }
  } else {
    // C_{ik} = A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[4 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[4 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[4 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 4;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] += alpha * (a[0] * b[0] + a[1] * b[2]);
            c[1] += alpha * (a[0] * b[1] + a[1] * b[3]);

            c[2] += alpha * (a[2] * b[0] + a[3] * b[2]);
            c[3] += alpha * (a[2] * b[1] + a[3] * b[3]);
          }
          b += 4;
        }
      }
    }
  }
}

/*!
  Apply a step of SOR to the system A*x = b.
*/
void BCSRMatApplySOR2(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
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

  // Store temporary data for each row
  TacsScalar t1, t2;

  if (start < end) {
    // Go through the matrix with the forward ordering
    for (int i = start; i < end; i++) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[2 * i];
      t2 = b[2 * i + 1];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[4 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        TacsScalar *y = &x[2 * j];

        if (i != j) {
          t1 -= a[0] * y[0] + a[1] * y[1];
          t2 -= a[2] * y[0] + a[3] * y[1];
        }

        // Increment the block pointer by bsize^2
        a += 4;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[4 * Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[2 * j];

          t1 -= a[0] * y[0] + a[1] * y[1];
          t2 -= a[2] * y[0] + a[3] * y[1];
          a += 4;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[4 * i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[2 * i] = (1.0 - omega) * x[2 * i] + omega * (d[0] * t1 + d[1] * t2);
      x[2 * i + 1] =
          (1.0 - omega) * x[2 * i + 1] + omega * (d[2] * t1 + d[3] * t2);
    }
  } else {
    // Go through the matrix with the forward ordering
    for (int i = start - 1; i >= end; i--) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[2 * i];
      t2 = b[2 * i + 1];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[4 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        TacsScalar *y = &x[2 * j];

        if (i != j) {
          t1 -= a[0] * y[0] + a[1] * y[1];
          t2 -= a[2] * y[0] + a[3] * y[1];
        }

        // Increment the block pointer by bsize^2
        a += 4;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[4 * Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[2 * j];

          t1 -= a[0] * y[0] + a[1] * y[1];
          t2 -= a[2] * y[0] + a[3] * y[1];
          a += 4;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[4 * i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[2 * i] = (1.0 - omega) * x[2 * i] + omega * (d[0] * t1 + d[1] * t2);
      x[2 * i + 1] =
          (1.0 - omega) * x[2 * i + 1] + omega * (d[2] * t1 + d[3] * t2);
    }
  }
}
