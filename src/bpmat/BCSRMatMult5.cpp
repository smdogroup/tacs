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
  The following file contains the specific implementation for block
  size = 5.
*/

/*!
  Compute the matrix-vector product: y = A * x
*/

void BCSRMatVecMult5(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *a = data->A;

  for (int i = 0; i < nrows; i++) {
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;
    y[4] = 0.0;

    int end = rowp[i + 1];
    int k = rowp[i];
    for (; k < end; k++) {
      int j = 5 * cols[k];

      y[0] += a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
              a[3] * x[j + 3] + a[4] * x[j + 4];
      y[1] += a[5] * x[j] + a[6] * x[j + 1] + a[7] * x[j + 2] +
              a[8] * x[j + 3] + a[9] * x[j + 4];
      y[2] += a[10] * x[j] + a[11] * x[j + 1] + a[12] * x[j + 2] +
              a[13] * x[j + 3] + a[14] * x[j + 4];
      y[3] += a[15] * x[j] + a[16] * x[j + 1] + a[17] * x[j + 2] +
              a[18] * x[j + 3] + a[19] * x[j + 4];
      y[4] += a[20] * x[j] + a[21] * x[j + 1] + a[22] * x[j + 2] +
              a[23] * x[j + 3] + a[24] * x[j + 4];
      a += 25;
    }
    y += 5;
  }
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/

void BCSRMatVecMultAdd5(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *a = data->A;

  for (int i = 0; i < nrows; i++) {
    y[0] = z[0];
    y[1] = z[1];
    y[2] = z[2];
    y[3] = z[3];
    y[4] = z[4];

    int end = rowp[i + 1];
    int k = rowp[i];
    for (; k < end; k++) {
      int j = 5 * cols[k];

      y[0] += a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
              a[3] * x[j + 3] + a[4] * x[j + 4];
      y[1] += a[5] * x[j] + a[6] * x[j + 1] + a[7] * x[j + 2] +
              a[8] * x[j + 3] + a[9] * x[j + 4];
      y[2] += a[10] * x[j] + a[11] * x[j + 1] + a[12] * x[j + 2] +
              a[13] * x[j + 3] + a[14] * x[j + 4];
      y[3] += a[15] * x[j] + a[16] * x[j + 1] + a[17] * x[j + 2] +
              a[18] * x[j + 3] + a[19] * x[j + 4];
      y[4] += a[20] * x[j] + a[21] * x[j + 1] + a[22] * x[j + 2] +
              a[23] * x[j + 3] + a[24] * x[j + 4];
      a += 25;
    }
    y += 5;
    z += 5;
  }
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower5(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *z = y;

  for (int i = 0; i < nrows; i++) {
    z[0] = x[0];
    z[1] = x[1];
    z[2] = x[2];
    z[3] = x[3];
    z[4] = x[4];

    int end = diag[i];
    int k = rowp[i];
    const TacsScalar *a = &A[25 * k];
    for (; k < end; k++) {
      int j = 5 * cols[k];

      z[0] -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
              a[3] * y[j + 3] + a[4] * y[j + 4];
      z[1] -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] +
              a[8] * y[j + 3] + a[9] * y[j + 4];
      z[2] -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
              a[13] * y[j + 3] + a[14] * y[j + 4];
      z[3] -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
              a[18] * y[j + 3] + a[19] * y[j + 4];
      z[4] -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
              a[23] * y[j + 3] + a[24] * y[j + 4];
      a += 25;
    }
    z += 5;
    x += 5;
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper5(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y0, y1, y2, y3, y4;

  x = &x[5 * (nrows - 1)];
  for (int i = nrows - 1; i >= 0; i--) {
    y0 = x[0];
    y1 = x[1];
    y2 = x[2];
    y3 = x[3];
    y4 = x[4];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    const TacsScalar *a = &A[25 * k];
    for (; k < end; k++) {
      int j = 5 * cols[k];

      y0 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] + a[3] * y[j + 3] +
            a[4] * y[j + 4];
      y1 -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] + a[8] * y[j + 3] +
            a[9] * y[j + 4];
      y2 -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
            a[13] * y[j + 3] + a[14] * y[j + 4];
      y3 -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
            a[18] * y[j + 3] + a[19] * y[j + 4];
      y4 -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
            a[23] * y[j + 3] + a[24] * y[j + 4];
      a += 25;
    }

    int bi = 5 * i;
    a = &A[25 * diag[i]];
    y[bi] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4;
    y[bi + 1] = a[5] * y0 + a[6] * y1 + a[7] * y2 + a[8] * y3 + a[9] * y4;
    y[bi + 2] = a[10] * y0 + a[11] * y1 + a[12] * y2 + a[13] * y3 + a[14] * y4;
    y[bi + 3] = a[15] * y0 + a[16] * y1 + a[17] * y2 + a[18] * y3 + a[19] * y4;
    y[bi + 4] = a[20] * y0 + a[21] * y1 + a[22] * y2 + a[23] * y3 + a[24] * y4;

    x -= 5;
  }
}

/*!
  Apply a portion of the lower factorization x = L^{-1} x
*/

void BCSRMatApplyPartialLower5(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *xx = &x[5];
  int off = 5 * var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    int end = diag[i];
    int k = rowp[i];
    while (cols[k] < var_offset) k++;

    const TacsScalar *a = &A[25 * k];
    for (; k < end; k++) {
      int j = 5 * cols[k] - off;

      xx[0] -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
               a[3] * x[j + 3] + a[4] * x[j + 4];
      xx[1] -= a[5] * x[j] + a[6] * x[j + 1] + a[7] * x[j + 2] +
               a[8] * x[j + 3] + a[9] * x[j + 4];
      xx[2] -= a[10] * x[j] + a[11] * x[j + 1] + a[12] * x[j + 2] +
               a[13] * x[j + 3] + a[14] * x[j + 4];
      xx[3] -= a[15] * x[j] + a[16] * x[j + 1] + a[17] * x[j + 2] +
               a[18] * x[j + 3] + a[19] * x[j + 4];
      xx[4] -= a[20] * x[j] + a[21] * x[j + 1] + a[22] * x[j + 2] +
               a[23] * x[j + 3] + a[24] * x[j + 4];
      a += 25;
    }
    xx += 5;
  }
}

/*!
  Apply a portion of he upper factorization x = U^{-1} x
*/

void BCSRMatApplyPartialUpper5(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y0, y1, y2, y3, y4;
  TacsScalar *xx = &x[5 * (nrows - var_offset - 1)];
  int off = 5 * var_offset;

  for (int i = nrows - 1; i >= var_offset; i--) {
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];
    y4 = xx[4];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    const TacsScalar *a = &A[25 * k];
    for (; k < end; k++) {
      int j = 5 * cols[k] - off;

      y0 -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] + a[3] * x[j + 3] +
            a[4] * x[j + 4];
      y1 -= a[5] * x[j] + a[6] * x[j + 1] + a[7] * x[j + 2] + a[8] * x[j + 3] +
            a[9] * x[j + 4];
      y2 -= a[10] * x[j] + a[11] * x[j + 1] + a[12] * x[j + 2] +
            a[13] * x[j + 3] + a[14] * x[j + 4];
      y3 -= a[15] * x[j] + a[16] * x[j + 1] + a[17] * x[j + 2] +
            a[18] * x[j + 3] + a[19] * x[j + 4];
      y4 -= a[20] * x[j] + a[21] * x[j + 1] + a[22] * x[j + 2] +
            a[23] * x[j + 3] + a[24] * x[j + 4];
      a += 25;
    }

    a = &A[25 * diag[i]];
    xx[0] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4;
    xx[1] = a[5] * y0 + a[6] * y1 + a[7] * y2 + a[8] * y3 + a[9] * y4;
    xx[2] = a[10] * y0 + a[11] * y1 + a[12] * y2 + a[13] * y3 + a[14] * y4;
    xx[3] = a[15] * y0 + a[16] * y1 + a[17] * y2 + a[18] * y3 + a[19] * y4;
    xx[4] = a[20] * y0 + a[21] * y1 + a[22] * y2 + a[23] * y3 + a[24] * y4;
    xx -= 5;
  }
}

/*!
  Function for the approximate Schur preconditioner
*/

void BCSRMatApplyFactorSchur5(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y0, y1, y2, y3, y4;
  TacsScalar *xx = &x[5 * (var_offset - 1)];

  for (int i = var_offset - 1; i >= 0; i--) {
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];
    y4 = xx[4];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    const TacsScalar *a = &A[25 * k];
    for (; k < end; k++) {
      int j = 5 * cols[k];

      y0 -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] + a[3] * x[j + 3] +
            a[4] * x[j + 4];
      y1 -= a[5] * x[j] + a[6] * x[j + 1] + a[7] * x[j + 2] + a[8] * x[j + 3] +
            a[9] * x[j + 4];
      y2 -= a[10] * x[j] + a[11] * x[j + 1] + a[12] * x[j + 2] +
            a[13] * x[j + 3] + a[14] * x[j + 4];
      y3 -= a[15] * x[j] + a[16] * x[j + 1] + a[17] * x[j + 2] +
            a[18] * x[j + 3] + a[19] * x[j + 4];
      y4 -= a[20] * x[j] + a[21] * x[j + 1] + a[22] * x[j + 2] +
            a[23] * x[j + 3] + a[24] * x[j + 4];
      a += 25;
    }

    a = &A[25 * diag[i]];
    xx[0] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4;
    xx[1] = a[5] * y0 + a[6] * y1 + a[7] * y2 + a[8] * y3 + a[9] * y4;
    xx[2] = a[10] * y0 + a[11] * y1 + a[12] * y2 + a[13] * y3 + a[14] * y4;
    xx[3] = a[15] * y0 + a[16] * y1 + a[17] * y2 + a[18] * y3 + a[19] * y4;
    xx[4] = a[20] * y0 + a[21] * y1 + a[22] * y2 + a[23] * y3 + a[24] * y4;
    xx -= 5;
  }
}

/*!
  Perform a matrix-matrix multiplication
*/

void BCSRMatMatMultAdd5(double alpha, BCSRMatData *Adata, BCSRMatData *Bdata,
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
    // C_{ik} += A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[25 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[25 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[25 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 25;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4;

            b0 = b[0];
            b1 = b[5];
            b2 = b[10];
            b3 = b[15];
            b4 = b[20];
            c[0] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[5] += a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[10] +=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[15] +=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[20] +=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[1];
            b1 = b[6];
            b2 = b[11];
            b3 = b[16];
            b4 = b[21];
            c[1] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[6] += a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[11] +=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[16] +=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[21] +=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[2];
            b1 = b[7];
            b2 = b[12];
            b3 = b[17];
            b4 = b[22];
            c[2] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[7] += a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[12] +=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[17] +=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[22] +=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[3];
            b1 = b[8];
            b2 = b[13];
            b3 = b[18];
            b4 = b[23];
            c[3] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[8] += a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[13] +=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[18] +=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[23] +=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[4];
            b1 = b[9];
            b2 = b[14];
            b3 = b[19];
            b4 = b[24];
            c[4] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[9] += a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[14] +=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[19] +=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[20] +=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;
          }
          b += 25;
        }
      }
    }
  } else if (alpha == -1.0) {
    // C_{ik} -= A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[25 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[25 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[25 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 25;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4;

            b0 = b[0];
            b1 = b[5];
            b2 = b[10];
            b3 = b[15];
            b4 = b[20];
            c[0] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[5] -= a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[10] -=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[15] -=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[20] -=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[1];
            b1 = b[6];
            b2 = b[11];
            b3 = b[16];
            b4 = b[21];
            c[1] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[6] -= a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[11] -=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[16] -=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[21] -=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[2];
            b1 = b[7];
            b2 = b[12];
            b3 = b[17];
            b4 = b[22];
            c[2] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[7] -= a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[12] -=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[17] -=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[22] -=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[3];
            b1 = b[8];
            b2 = b[13];
            b3 = b[18];
            b4 = b[23];
            c[3] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[8] -= a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[13] -=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[18] -=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[23] -=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

            b0 = b[4];
            b1 = b[9];
            b2 = b[14];
            b3 = b[19];
            b4 = b[24];
            c[4] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
            c[9] -= a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
            c[14] -=
                a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
            c[19] -=
                a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
            c[20] -=
                a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;
          }

          b += 25;
        }
      }
    }
  } else {
    // C_{ik} += alpha*A_{ij} B_{jk}
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[25 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[25 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[25 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 25;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4;

            b0 = b[0];
            b1 = b[5];
            b2 = b[10];
            b3 = b[15];
            b4 = b[20];
            c[0] += alpha *
                    (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4);
            c[5] += alpha *
                    (a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4);
            c[10] += alpha * (a[10] * b0 + a[11] * b1 + a[12] * b2 +
                              a[13] * b3 + a[14] * b4);
            c[15] += alpha * (a[15] * b0 + a[16] * b1 + a[17] * b2 +
                              a[18] * b3 + a[19] * b4);
            c[20] += alpha * (a[20] * b0 + a[21] * b1 + a[22] * b2 +
                              a[23] * b3 + a[24] * b4);

            b0 = b[1];
            b1 = b[6];
            b2 = b[11];
            b3 = b[16];
            b4 = b[21];
            c[1] += alpha *
                    (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4);
            c[6] += alpha *
                    (a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4);
            c[11] += alpha * (a[10] * b0 + a[11] * b1 + a[12] * b2 +
                              a[13] * b3 + a[14] * b4);
            c[16] += alpha * (a[15] * b0 + a[16] * b1 + a[17] * b2 +
                              a[18] * b3 + a[19] * b4);
            c[21] += alpha * (a[20] * b0 + a[21] * b1 + a[22] * b2 +
                              a[23] * b3 + a[24] * b4);

            b0 = b[2];
            b1 = b[7];
            b2 = b[12];
            b3 = b[17];
            b4 = b[22];
            c[2] += alpha *
                    (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4);
            c[7] += alpha *
                    (a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4);
            c[12] += alpha * (a[10] * b0 + a[11] * b1 + a[12] * b2 +
                              a[13] * b3 + a[14] * b4);
            c[17] += alpha * (a[15] * b0 + a[16] * b1 + a[17] * b2 +
                              a[18] * b3 + a[19] * b4);
            c[22] += alpha * (a[20] * b0 + a[21] * b1 + a[22] * b2 +
                              a[23] * b3 + a[24] * b4);

            b0 = b[3];
            b1 = b[8];
            b2 = b[13];
            b3 = b[18];
            b4 = b[23];
            c[3] += alpha *
                    (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4);
            c[8] += alpha *
                    (a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4);
            c[13] += alpha * (a[10] * b0 + a[11] * b1 + a[12] * b2 +
                              a[13] * b3 + a[14] * b4);
            c[18] += alpha * (a[15] * b0 + a[16] * b1 + a[17] * b2 +
                              a[18] * b3 + a[19] * b4);
            c[23] += alpha * (a[20] * b0 + a[21] * b1 + a[22] * b2 +
                              a[23] * b3 + a[24] * b4);

            b0 = b[4];
            b1 = b[9];
            b2 = b[14];
            b3 = b[19];
            b4 = b[24];
            c[4] += alpha *
                    (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4);
            c[9] += alpha *
                    (a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4);
            c[14] += alpha * (a[10] * b0 + a[11] * b1 + a[12] * b2 +
                              a[13] * b3 + a[14] * b4);
            c[19] += alpha * (a[15] * b0 + a[16] * b1 + a[17] * b2 +
                              a[18] * b3 + a[19] * b4);
            c[20] += alpha * (a[20] * b0 + a[21] * b1 + a[22] * b2 +
                              a[23] * b3 + a[24] * b4);
          }

          b += 25;
        }
      }
    }
  }
}

/*!
  Apply a step of SOR to the system A*x = b - B*xext
*/
void BCSRMatApplySOR5(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
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
  TacsScalar t1, t2, t3, t4, t5;

  if (start < end) {
    // Go through the matrix with the forward ordering
    for (int i = start; i < end; i++) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[5 * i];
      t2 = b[5 * i + 1];
      t3 = b[5 * i + 2];
      t4 = b[5 * i + 3];
      t5 = b[5 * i + 4];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[25 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        TacsScalar *y = &x[5 * j];

        if (i != j) {
          t1 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                a[3] * y[j + 3] + a[4] * y[j + 4];
          t2 -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] +
                a[8] * y[j + 3] + a[9] * y[j + 4];
          t3 -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
                a[13] * y[j + 3] + a[14] * y[j + 4];
          t4 -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
                a[18] * y[j + 3] + a[19] * y[j + 4];
          t5 -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
                a[23] * y[j + 3] + a[24] * y[j + 4];
        }

        // Increment the block pointer by bsize^2
        a += 25;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[25 * Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[5 * j];

          t1 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                a[3] * y[j + 3] + a[4] * y[j + 4];
          t2 -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] +
                a[8] * y[j + 3] + a[9] * y[j + 4];
          t3 -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
                a[13] * y[j + 3] + a[14] * y[j + 4];
          t4 -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
                a[18] * y[j + 3] + a[19] * y[j + 4];
          t5 -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
                a[23] * y[j + 3] + a[24] * y[j + 4];
          a += 25;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[25 * i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[5 * i] =
          (1.0 - omega) * x[5 * i] +
          omega * (d[0] * t1 + d[1] * t2 + d[2] * t3 + d[3] * t4 + d[4] * t5);
      x[5 * i + 1] =
          (1.0 - omega) * x[5 * i + 1] +
          omega * (d[5] * t1 + d[6] * t2 + d[7] * t3 + d[8] * t4 + d[9] * t5);
      x[5 * i + 2] = (1.0 - omega) * x[5 * i + 2] +
                     omega * (d[10] * t1 + d[11] * t2 + d[12] * t3 +
                              d[13] * t4 + d[14] * t5);
      x[5 * i + 3] = (1.0 - omega) * x[5 * i + 3] +
                     omega * (d[15] * t1 + d[16] * t2 + d[17] * t3 +
                              d[18] * t4 + d[19] * t5);
      x[5 * i + 4] = (1.0 - omega) * x[5 * i + 4] +
                     omega * (d[20] * t1 + d[21] * t2 + d[22] * t3 +
                              d[23] * t4 + d[24] * t5);
    }
  } else {
    // Go through the matrix with the forward ordering
    for (int i = start - 1; i >= end; i--) {
      // Copy the right-hand-side to the temporary vector for this row
      t1 = b[5 * i];
      t2 = b[5 * i + 1];
      t3 = b[5 * i + 2];
      t4 = b[5 * i + 3];
      t5 = b[5 * i + 4];

      // Set the pointer to the beginning of the current row
      const TacsScalar *a = &Adata->A[25 * Arowp[i]];

      // Scan through the row and compute the result:
      // tx <- b_i - A_{ij}*x_{j} for j != i
      int end = Arowp[i + 1];
      for (int k = Arowp[i]; k < end; k++) {
        int j = Acols[k];
        TacsScalar *y = &x[5 * j];

        if (i != j) {
          t1 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                a[3] * y[j + 3] + a[4] * y[j + 4];
          t2 -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] +
                a[8] * y[j + 3] + a[9] * y[j + 4];
          t3 -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
                a[13] * y[j + 3] + a[14] * y[j + 4];
          t4 -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
                a[18] * y[j + 3] + a[19] * y[j + 4];
          t5 -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
                a[23] * y[j + 3] + a[24] * y[j + 4];
        }

        // Increment the block pointer by bsize^2
        a += 25;
      }

      if (Bdata && i >= var_offset) {
        const int row = i - var_offset;

        // Set the pointer to the row in B
        a = &Bdata->A[25 * Browp[row]];
        end = Browp[row + 1];
        for (int k = Browp[row]; k < end; k++) {
          int j = Bcols[k];
          const TacsScalar *y = &xext[5 * j];

          t1 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                a[3] * y[j + 3] + a[4] * y[j + 4];
          t2 -= a[5] * y[j] + a[6] * y[j + 1] + a[7] * y[j + 2] +
                a[8] * y[j + 3] + a[9] * y[j + 4];
          t3 -= a[10] * y[j] + a[11] * y[j + 1] + a[12] * y[j + 2] +
                a[13] * y[j + 3] + a[14] * y[j + 4];
          t4 -= a[15] * y[j] + a[16] * y[j + 1] + a[17] * y[j + 2] +
                a[18] * y[j + 3] + a[19] * y[j + 4];
          t5 -= a[20] * y[j] + a[21] * y[j + 1] + a[22] * y[j + 2] +
                a[23] * y[j + 3] + a[24] * y[j + 4];
          a += 25;
        }
      }

      // Set a pointer to the inverse of the diagonal
      const TacsScalar *d = &Adiag[25 * i];

      // Compute the first term in the update:
      // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
      x[5 * i] =
          (1.0 - omega) * x[5 * i] +
          omega * (d[0] * t1 + d[1] * t2 + d[2] * t3 + d[3] * t4 + d[4] * t5);
      x[5 * i + 1] =
          (1.0 - omega) * x[5 * i + 1] +
          omega * (d[5] * t1 + d[6] * t2 + d[7] * t3 + d[8] * t4 + d[9] * t5);
      x[5 * i + 2] = (1.0 - omega) * x[5 * i + 2] +
                     omega * (d[10] * t1 + d[11] * t2 + d[12] * t3 +
                              d[13] * t4 + d[14] * t5);
      x[5 * i + 3] = (1.0 - omega) * x[5 * i + 3] +
                     omega * (d[15] * t1 + d[16] * t2 + d[17] * t3 +
                              d[18] * t4 + d[19] * t5);
      x[5 * i + 4] = (1.0 - omega) * x[5 * i + 4] +
                     omega * (d[20] * t1 + d[21] * t2 + d[22] * t3 +
                              d[23] * t4 + d[24] * t5);
    }
  }
}
