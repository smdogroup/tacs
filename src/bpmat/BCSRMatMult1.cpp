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
  Code for the 1x1 block case
*/

/*!
  Compute the matrix-vector product: y = A * x
*/
void BCSRMatVecMult1(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  for (int i = 0; i < nrows; i++) {
    y[0] = 0.0;

    int end = rowp[i + 1];
    int k = rowp[i];
    const TacsScalar *a = &A[k];
    for (; k < end; k++) {
      int bj = cols[k];

      y[0] += a[0] * x[bj];
      a++;
    }
    y++;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/
void BCSRMatVecMultAdd1(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  for (int i = 0; i < nrows; i++) {
    z[0] = y[0];

    int end = rowp[i + 1];
    int k = rowp[i];
    const TacsScalar *a = &A[k];
    for (; k < end; k++) {
      int bj = cols[k];

      z[0] += a[0] * x[bj];
      a++;
    }
    y++;
    z++;
  }
}

/*!
  Compute y = A^{T} * x
*/
void BCSRMatVecMultTranspose1(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  for (int i = 0; i < nrows; i++) {
    int end = rowp[i + 1];
    int k = rowp[i];
    const TacsScalar *a = &A[k];
    for (; k < end; k++) {
      int bj = cols[k];

      y[bj] += a[0] * x[0];
      a++;
    }
    x++;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower1(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *yy = y;

  for (int i = 0; i < nrows; i++) {
    yy[0] = x[0];

    int end = diag[i];
    int j = rowp[i];
    const TacsScalar *a = &A[j];
    for (; j < end; j++) {
      int bj = cols[j];

      yy[0] -= a[0] * y[bj];
      a++;
    }
    x++;
    yy++;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper1(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *yy = &y[(nrows - 1)];
  x = &x[(nrows - 1)];

  for (int i = nrows - 1; i >= 0; i--) {
    TacsScalar y1 = x[0];

    int end = rowp[i + 1];
    int j = diag[i];
    const TacsScalar *a = &A[j];
    TacsScalar a11 = a[0];
    j++;
    a++;

    for (; j < end; j++) {
      int bj = cols[j];

      y1 -= a[0] * y[bj];
      a++;
    }

    // Apply the inverse on the diagonal
    yy[0] = a11 * y1;

    yy--;
    x--;
  }
}

/*!
  Apply the lower factorization x = L^{-1} x
*/
void BCSRMatApplyPartialLower1(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *xx = &x[1];
  int off = var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    int end = diag[i];
    int j = rowp[i];
    while (cols[j] < var_offset) j++;

    const TacsScalar *a = &A[j];
    for (; j < end; j++) {
      int bj = cols[j] - off;

      xx[0] -= a[0] * x[bj];
      a++;
    }
    xx++;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyPartialUpper1(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  int off = var_offset;

  TacsScalar *xx = &x[(nrows - var_offset - 1)];

  for (int i = nrows - 1; i >= var_offset; i--) {
    TacsScalar y1 = xx[0];

    int end = rowp[i + 1];
    int j = diag[i];
    const TacsScalar *adiag = &A[j];
    const TacsScalar *a = adiag;
    j++;
    a++;

    for (; j < end; j++) {
      int bj = cols[j] - off;

      y1 -= a[0] * x[bj];
      a++;
    }

    xx[0] = adiag[0] * y1;
    xx--;
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
void BCSRMatApplyFactorSchur1(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *xx = &x[(var_offset - 1)];

  for (int i = var_offset - 1; i >= 0; i--) {
    TacsScalar y1 = xx[0];

    int j = diag[i];
    int end = rowp[i + 1];
    const TacsScalar *adiag = &A[j];
    const TacsScalar *a = adiag;
    a++;
    j++;

    for (; j < end; j++) {
      int bj = cols[j];
      y1 -= a[0] * x[bj];
      a++;
    }

    // Apply the inverse on the diagonal
    xx[0] = adiag[0] * y1;
    xx--;
  }
}

/*!
  Perform a matrix-matrix multiplication
*/
void BCSRMatMatMultAdd1(double alpha, BCSRMatData *Adata, BCSRMatData *Bdata,
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
        const TacsScalar *a = &A[jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c++;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] += a[0] * b[0];
          }
          b++;
        }
      }
    }
  } else if (alpha == -1.0) {
    // C_{ik} = A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c++;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] -= a[0] * b[0];
          }
          b++;
        }
      }
    }
  } else {
    // C_{ik} = A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c++;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            c[0] += alpha * (a[0] * b[0]);
          }
          b++;
        }
      }
    }
  }
}

/*!
  Compute the scaled normal equations:

  A = B^{T}*s*B
*/
void BCSRMatMatMultNormal1(BCSRMatData *Adata, TacsScalar *scale,
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

  int *kptr = new int[nrows_b];
  memcpy(kptr, browp, nrows_b * sizeof(int));

  // A_{ij} = B_{ki}*s{k}*B_{kj}
  for (int i = 0; i < nrows_a; i++) {
    // Scan through column i of the matrix B_{*i}
    for (int k = 0; k < nrows_b; k++) {
      if ((kptr[k] < browp[k + 1]) && (bcols[kptr[k]] == i)) {
        TacsScalar bi = B[kptr[k]];
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

          // If the column indices are equal, compute the inner product
          if (acols[jpa] == bcols[jpb]) {
            A[jpa] += scale[k] * bi * B[jpb];
          }
        }
      }
    }
  }

  delete[] kptr;
}

/*!
  Apply a step of SOR to the system A*x = b.
*/
void BCSRMatApplySOR1(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
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
  TacsScalar t1;

  if (start < end) {
    // Go through the matrix with the forward ordering
    for (int i = start; i < end; i++) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[i];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        if (i != j) {
          t1 -= a[0] * x[j];
        }

        // Increment the block pointer by bsize^2
        a++;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          t1 -= a[0] * xext[j];
          a++;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[i] = (1.0 - omega) * x[i] + omega * d[0] * t1;
    }
  } else {
    // Go through the matrix with the forward ordering
    for (int i = start - 1; i >= end; i--) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[i];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        if (i != j) {
          t1 -= a[0] * x[j];
        }

        // Increment the block pointer by bsize^2
        a++;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          t1 -= a[0] * xext[j];
          a++;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[i] = (1.0 - omega) * x[i] + omega * d[0] * t1;
    }
  }
}
