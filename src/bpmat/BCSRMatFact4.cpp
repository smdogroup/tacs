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
  Block size = 4 code
*/

/*!
  Perform an ILU factorization in place for the block size = 4.
  The entries are over-written, all operations are performed in place.
*/

void BCSRMatFactor4(BCSRMatData* data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  TacsScalar* A = data->A;

  TacsScalar d00, d01, d02, d03;
  TacsScalar d10, d11, d12, d13;
  TacsScalar d20, d21, d22, d23;
  TacsScalar d30, d31, d32, d33;

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr,
              "Error in factorization: no diagonal entry for \
row %d \n",
              i);
      return;
    }
    // Scan from the first entry in the current row, towards the diagonal
    int row_end = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];
      TacsScalar* a = &A[16 * j];
      TacsScalar* b = &A[16 * diag[cj]];

      d00 = a[0] * b[0] + a[1] * b[4] + a[2] * b[8] + a[3] * b[12];
      d10 = a[4] * b[0] + a[5] * b[4] + a[6] * b[8] + a[7] * b[12];
      d20 = a[8] * b[0] + a[9] * b[4] + a[10] * b[8] + a[11] * b[12];
      d30 = a[12] * b[0] + a[13] * b[4] + a[14] * b[8] + a[15] * b[12];

      d01 = a[0] * b[1] + a[1] * b[5] + a[2] * b[9] + a[3] * b[13];
      d11 = a[4] * b[1] + a[5] * b[5] + a[6] * b[9] + a[7] * b[13];
      d21 = a[8] * b[1] + a[9] * b[5] + a[10] * b[9] + a[11] * b[13];
      d31 = a[12] * b[1] + a[13] * b[5] + a[14] * b[9] + a[15] * b[13];

      d02 = a[0] * b[2] + a[1] * b[6] + a[2] * b[10] + a[3] * b[14];
      d12 = a[4] * b[2] + a[5] * b[6] + a[6] * b[10] + a[7] * b[14];
      d22 = a[8] * b[2] + a[9] * b[6] + a[10] * b[10] + a[11] * b[14];
      d32 = a[12] * b[2] + a[13] * b[6] + a[14] * b[10] + a[15] * b[14];

      d03 = a[0] * b[3] + a[1] * b[7] + a[2] * b[11] + a[3] * b[15];
      d13 = a[4] * b[3] + a[5] * b[7] + a[6] * b[11] + a[7] * b[15];
      d23 = a[8] * b[3] + a[9] * b[7] + a[10] * b[11] + a[11] * b[15];
      d33 = a[12] * b[3] + a[13] * b[7] + a[14] * b[11] + a[15] * b[15];

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;
      a = &A[16 * k];
      b = &A[16 * p];

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          k++;
          a += 16;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          TacsScalar b0, b1, b2, b3;
          b0 = b[0];
          b1 = b[4];
          b2 = b[8];
          b3 = b[12];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[4] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[8] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[12] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[1];
          b1 = b[5];
          b2 = b[9];
          b3 = b[13];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[5] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[9] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[13] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[2];
          b1 = b[6];
          b2 = b[10];
          b3 = b[14];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[6] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[10] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[14] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[3];
          b1 = b[7];
          b2 = b[11];
          b3 = b[15];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[7] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[11] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[15] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;
        }

        b += 16;
      }

      // Copy the matrix back into the row
      a = &A[16 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d10;
      a[5] = d11;
      a[6] = d12;
      a[7] = d13;
      a[8] = d20;
      a[9] = d21;
      a[10] = d22;
      a[11] = d23;
      a[12] = d30;
      a[13] = d31;
      a[14] = d32;
      a[15] = d33;
    }

    // Invert the diagonal portion of the matrix
    TacsScalar D[16];
    TacsScalar* a = &A[16 * diag[i]];
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

    int ipiv[4];
    int info = BMatComputeInverse(a, D, ipiv, 4);

    if (info > 0) {
      fprintf(stderr,
              "Error during factorization of diagonal block %d \
in row %d \n",
              i + 1, info);
    }
  }
}

/*!
  Compute x = L_{B}^{-1} E
*/

void BCSRMatFactorLower4(BCSRMatData* data, BCSRMatData* Edata) {
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
      const TacsScalar* d = &A[16 * j];

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar* a = &E[16 * k];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar* b = &E[16 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 16;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          TacsScalar b0, b1, b2, b3;
          b0 = b[0];
          b1 = b[4];
          b2 = b[8];
          b3 = b[12];
          a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3;
          a[4] -= d[4] * b0 + d[5] * b1 + d[6] * b2 + d[7] * b3;
          a[8] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3;
          a[12] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3;

          b0 = b[1];
          b1 = b[5];
          b2 = b[9];
          b3 = b[13];
          a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3;
          a[5] -= d[4] * b0 + d[5] * b1 + d[6] * b2 + d[7] * b3;
          a[9] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3;
          a[13] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3;

          b0 = b[2];
          b1 = b[6];
          b2 = b[10];
          b3 = b[14];
          a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3;
          a[6] -= d[4] * b0 + d[5] * b1 + d[6] * b2 + d[7] * b3;
          a[10] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3;
          a[14] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3;

          b0 = b[3];
          b1 = b[7];
          b2 = b[11];
          b3 = b[15];
          a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3;
          a[7] -= d[4] * b0 + d[5] * b1 + d[6] * b2 + d[7] * b3;
          a[11] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3;
          a[15] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3;
        }
        b += 16;
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/

void BCSRMatFactorUpper4(BCSRMatData* data, BCSRMatData* Fdata) {
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

  TacsScalar d00, d01, d02, d03;
  TacsScalar d10, d11, d12, d13;
  TacsScalar d20, d21, d22, d23;
  TacsScalar d30, d31, d32, d33;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar* a = &F[16 * j];
      const TacsScalar* b = &A[16 * diag[cj]];

      // Multiply d = F[j] * A[diag[cj]]
      TacsScalar b0, b1, b2, b3;

      b0 = b[0];
      b1 = b[4];
      b2 = b[8];
      b3 = b[12];
      d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3;
      d10 = a[4] * b0 + a[5] * b1 + a[6] * b2 + a[7] * b3;
      d20 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3;
      d30 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3;

      b0 = b[1];
      b1 = b[5];
      b2 = b[9];
      b3 = b[13];
      d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3;
      d11 = a[4] * b0 + a[5] * b1 + a[6] * b2 + a[7] * b3;
      d21 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3;
      d31 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3;

      b0 = b[2];
      b1 = b[6];
      b2 = b[10];
      b3 = b[14];
      d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3;
      d12 = a[4] * b0 + a[5] * b1 + a[6] * b2 + a[7] * b3;
      d22 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3;
      d32 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3;

      b0 = b[3];
      b1 = b[7];
      b2 = b[11];
      b3 = b[15];
      d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3;
      d13 = a[4] * b0 + a[5] * b1 + a[6] * b2 + a[7] * b3;
      d23 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3;
      d33 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3;

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &F[16 * k];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &A[16 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 16;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[4];
          b2 = b[8];
          b3 = b[12];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[4] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[8] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[12] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[1];
          b1 = b[5];
          b2 = b[9];
          b3 = b[13];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[5] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[9] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[13] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[2];
          b1 = b[6];
          b2 = b[10];
          b3 = b[14];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[6] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[10] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[14] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;

          b0 = b[3];
          b1 = b[7];
          b2 = b[11];
          b3 = b[15];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3;
          a[7] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3;
          a[11] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3;
          a[15] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3;
        }
        b += 16;
      }

      // Copy over the matrix
      a = &F[16 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d10;
      a[5] = d11;
      a[6] = d12;
      a[7] = d13;
      a[8] = d20;
      a[9] = d21;
      a[10] = d22;
      a[11] = d23;
      a[12] = d30;
      a[13] = d31;
      a[14] = d32;
      a[15] = d33;
    }
  }
}
