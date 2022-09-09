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
  Block size = 5 code
*/

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern.  The entries are over-written, all operations are
  performed in place.
*/

void BCSRMatFactor5(BCSRMatData* data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  TacsScalar* A = data->A;

  TacsScalar d00, d01, d02, d03, d04;
  TacsScalar d10, d11, d12, d13, d14;
  TacsScalar d20, d21, d22, d23, d24;
  TacsScalar d30, d31, d32, d33, d34;
  TacsScalar d40, d41, d42, d43, d44;

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr, "Error in factorization: no diagonal entry for row %d",
              i);
      return;
    }

    // Scan from the first entry in the current row, towards the diagonal
    int row_end = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];
      TacsScalar* a = &A[25 * j];
      TacsScalar* b = &A[25 * diag[cj]];

      // Multiply d = A[j] * A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4;

      b0 = b[0];
      b1 = b[5];
      b2 = b[10];
      b3 = b[15];
      b4 = b[20];
      d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d10 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d20 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d30 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d40 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[1];
      b1 = b[6];
      b2 = b[11];
      b3 = b[16];
      b4 = b[21];
      d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d11 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d21 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d31 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d41 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[2];
      b1 = b[7];
      b2 = b[12];
      b3 = b[17];
      b4 = b[22];
      d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d12 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d22 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d32 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d42 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[3];
      b1 = b[8];
      b2 = b[13];
      b3 = b[18];
      b4 = b[23];
      d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d13 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d23 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d33 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d43 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[4];
      b1 = b[9];
      b2 = b[14];
      b3 = b[19];
      b4 = b[24];
      d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d14 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d24 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d34 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d44 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;
      a = &A[25 * k];
      b = &A[25 * p];

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          k++;
          a += 25;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[5];
          b2 = b[10];
          b3 = b[15];
          b4 = b[20];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[5] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[10] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[15] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[20] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[1];
          b1 = b[6];
          b2 = b[11];
          b3 = b[16];
          b4 = b[21];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[6] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[11] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[16] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[21] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[2];
          b1 = b[7];
          b2 = b[12];
          b3 = b[17];
          b4 = b[22];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[7] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[12] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[17] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[22] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[3];
          b1 = b[8];
          b2 = b[13];
          b3 = b[18];
          b4 = b[23];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[13] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[18] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[23] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[4];
          b1 = b[9];
          b2 = b[14];
          b3 = b[19];
          b4 = b[24];
          a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[14] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[19] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[24] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;
        }

        b += 25;
      }

      // Copy the matrix back into the row
      a = &A[25 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d10;
      a[6] = d11;
      a[7] = d12;
      a[8] = d13;
      a[9] = d14;
      a[10] = d20;
      a[11] = d21;
      a[12] = d22;
      a[13] = d23;
      a[14] = d24;
      a[15] = d30;
      a[16] = d31;
      a[17] = d32;
      a[18] = d33;
      a[19] = d34;
      a[20] = d40;
      a[21] = d41;
      a[22] = d42;
      a[23] = d43;
      a[24] = d44;
    }

    // Invert the diagonal portion of the matrix
    TacsScalar D[25];
    TacsScalar* a = &A[25 * diag[i]];
    D[0] = a[0];
    D[1] = a[1];
    D[2] = a[2];
    D[3] = a[3];
    D[4] = a[4];
    D[5] = a[5];
    D[6] = a[6];
    D[7] = a[7];
    D[8] = a[8];
    D[9] = a[9];
    D[10] = a[10];
    D[11] = a[11];
    D[12] = a[12];
    D[13] = a[13];
    D[14] = a[14];
    D[15] = a[15];
    D[16] = a[16];
    D[17] = a[17];
    D[18] = a[18];
    D[19] = a[19];
    D[20] = a[20];
    D[21] = a[21];
    D[22] = a[22];
    D[23] = a[23];
    D[24] = a[24];

    int ipiv[6];
    int info = BMatComputeInverse(a, D, ipiv, 5);

    if (info > 0) {
      fprintf(stderr,
              "Error during factorization of diagonal %d in \
block row %d \n",
              i + 1, info);
    }
  }
}

/*!
  Compute x = L_{B}^{-1} E
*/

void BCSRMatFactorLower5(BCSRMatData* data, BCSRMatData* Edata) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  const TacsScalar* A = data->A;

  // Retrieve the data required from the matrix
  const int* erowp = Edata->rowp;
  const int* ecols = Edata->cols;
  TacsScalar* E = Edata->A;

  for (int i = 0; i < nrows; i++) {
    // Scan from the first entry in the current row, towards the diagonal
    int j_end = diag[i];

    for (int j = rowp[i]; j < j_end; j++) {
      int cj = cols[j];
      const TacsScalar* d = &A[25 * j];

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar* a = &E[25 * k];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar* b = &E[25 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 25;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          TacsScalar b0, b1, b2, b3, b4;
          b0 = b[0];
          b1 = b[5];
          b2 = b[10];
          b3 = b[15];
          b4 = b[20];
          a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4;
          a[5] -= d[5] * b0 + d[6] * b1 + d[7] * b2 + d[8] * b3 + d[9] * b4;
          a[10] -=
              d[10] * b0 + d[11] * b1 + d[12] * b2 + d[13] * b3 + d[14] * b4;
          a[15] -=
              d[15] * b0 + d[16] * b1 + d[17] * b2 + d[18] * b3 + d[19] * b4;
          a[20] -=
              d[20] * b0 + d[21] * b1 + d[22] * b2 + d[23] * b3 + d[24] * b4;

          b0 = b[1];
          b1 = b[6];
          b2 = b[11];
          b3 = b[16];
          b4 = b[21];
          a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4;
          a[6] -= d[5] * b0 + d[6] * b1 + d[7] * b2 + d[8] * b3 + d[9] * b4;
          a[11] -=
              d[10] * b0 + d[11] * b1 + d[12] * b2 + d[13] * b3 + d[14] * b4;
          a[16] -=
              d[15] * b0 + d[16] * b1 + d[17] * b2 + d[18] * b3 + d[19] * b4;
          a[21] -=
              d[20] * b0 + d[21] * b1 + d[22] * b2 + d[23] * b3 + d[24] * b4;

          b0 = b[2];
          b1 = b[7];
          b2 = b[12];
          b3 = b[17];
          b4 = b[22];
          a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4;
          a[7] -= d[5] * b0 + d[6] * b1 + d[7] * b2 + d[8] * b3 + d[9] * b4;
          a[12] -=
              d[10] * b0 + d[11] * b1 + d[12] * b2 + d[13] * b3 + d[14] * b4;
          a[17] -=
              d[15] * b0 + d[16] * b1 + d[17] * b2 + d[18] * b3 + d[19] * b4;
          a[22] -=
              d[20] * b0 + d[21] * b1 + d[22] * b2 + d[23] * b3 + d[24] * b4;

          b0 = b[3];
          b1 = b[8];
          b2 = b[13];
          b3 = b[18];
          b4 = b[23];
          a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4;
          a[8] -= d[5] * b0 + d[6] * b1 + d[7] * b2 + d[8] * b3 + d[9] * b4;
          a[13] -=
              d[10] * b0 + d[11] * b1 + d[12] * b2 + d[13] * b3 + d[14] * b4;
          a[18] -=
              d[15] * b0 + d[16] * b1 + d[17] * b2 + d[18] * b3 + d[19] * b4;
          a[23] -=
              d[20] * b0 + d[21] * b1 + d[22] * b2 + d[23] * b3 + d[24] * b4;

          b0 = b[4];
          b1 = b[9];
          b2 = b[14];
          b3 = b[19];
          b4 = b[24];
          a[4] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4;
          a[9] -= d[5] * b0 + d[6] * b1 + d[7] * b2 + d[8] * b3 + d[9] * b4;
          a[14] -=
              d[10] * b0 + d[11] * b1 + d[12] * b2 + d[13] * b3 + d[14] * b4;
          a[19] -=
              d[15] * b0 + d[16] * b1 + d[17] * b2 + d[18] * b3 + d[19] * b4;
          a[24] -=
              d[20] * b0 + d[21] * b1 + d[22] * b2 + d[23] * b3 + d[24] * b4;
        }
        b += 25;
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/

void BCSRMatFactorUpper5(BCSRMatData* data, BCSRMatData* Fdata) {
  // Retrieve the data required from the matrix
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  const TacsScalar* A = data->A;

  // Retrieve the data required from the matrix
  const int nrows_f = Fdata->nrows;
  const int* frowp = Fdata->rowp;
  const int* fcols = Fdata->cols;
  TacsScalar* F = Fdata->A;

  TacsScalar d00, d01, d02, d03, d04;
  TacsScalar d10, d11, d12, d13, d14;
  TacsScalar d20, d21, d22, d23, d24;
  TacsScalar d30, d31, d32, d33, d34;
  TacsScalar d40, d41, d42, d43, d44;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar* a = &F[25 * j];
      const TacsScalar* b = &A[25 * diag[cj]];

      // Multiply d = F[j] * A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4;

      b0 = b[0];
      b1 = b[5];
      b2 = b[10];
      b3 = b[15];
      b4 = b[20];
      d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d10 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d20 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d30 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d40 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[1];
      b1 = b[6];
      b2 = b[11];
      b3 = b[16];
      b4 = b[21];
      d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d11 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d21 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d31 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d41 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[2];
      b1 = b[7];
      b2 = b[12];
      b3 = b[17];
      b4 = b[22];
      d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d12 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d22 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d32 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d42 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[3];
      b1 = b[8];
      b2 = b[13];
      b3 = b[18];
      b4 = b[23];
      d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d13 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d23 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d33 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d43 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      b0 = b[4];
      b1 = b[9];
      b2 = b[14];
      b3 = b[19];
      b4 = b[24];
      d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4;
      d14 = a[5] * b0 + a[6] * b1 + a[7] * b2 + a[8] * b3 + a[9] * b4;
      d24 = a[10] * b0 + a[11] * b1 + a[12] * b2 + a[13] * b3 + a[14] * b4;
      d34 = a[15] * b0 + a[16] * b1 + a[17] * b2 + a[18] * b3 + a[19] * b4;
      d44 = a[20] * b0 + a[21] * b1 + a[22] * b2 + a[23] * b3 + a[24] * b4;

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &F[25 * k];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &A[25 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 25;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[5];
          b2 = b[10];
          b3 = b[15];
          b4 = b[20];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[5] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[10] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[15] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[20] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[1];
          b1 = b[6];
          b2 = b[11];
          b3 = b[16];
          b4 = b[21];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[6] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[11] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[16] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[21] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[2];
          b1 = b[7];
          b2 = b[12];
          b3 = b[17];
          b4 = b[22];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[7] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[12] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[17] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[22] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[3];
          b1 = b[8];
          b2 = b[13];
          b3 = b[18];
          b4 = b[23];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[13] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[18] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[23] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;

          b0 = b[4];
          b1 = b[9];
          b2 = b[14];
          b3 = b[19];
          b4 = b[24];
          a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4;
          a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4;
          a[14] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4;
          a[19] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4;
          a[24] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4;
        }
        b += 25;
      }

      // Copy over the matrix
      a = &F[25 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d10;
      a[6] = d11;
      a[7] = d12;
      a[8] = d13;
      a[9] = d14;
      a[10] = d20;
      a[11] = d21;
      a[12] = d22;
      a[13] = d23;
      a[14] = d24;
      a[15] = d30;
      a[16] = d31;
      a[17] = d32;
      a[18] = d33;
      a[19] = d34;
      a[20] = d40;
      a[21] = d41;
      a[22] = d42;
      a[23] = d43;
      a[24] = d44;
    }
  }
}
