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
  Block size = 3 code
*/

/*!
  Perform an ILU factorization in place for the block size = 3.
  The entries are over-written, all operations are performed in place.
*/

void BCSRMatFactor3(BCSRMatData* data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  TacsScalar* A = data->A;

  TacsScalar d00, d01, d02;
  TacsScalar d10, d11, d12;
  TacsScalar d20, d21, d22;

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
      TacsScalar* a = &A[9 * j];
      TacsScalar* b = &A[9 * diag[cj]];

      d00 = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
      d10 = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
      d20 = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];

      d01 = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
      d11 = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
      d21 = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];

      d02 = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];
      d12 = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];
      d22 = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;
      a = &A[9 * k];
      b = &A[9 * p];

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          k++;
          a += 9;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          a[0] -= d00 * b[0] + d01 * b[3] + d02 * b[6];
          a[3] -= d10 * b[0] + d11 * b[3] + d12 * b[6];
          a[6] -= d20 * b[0] + d21 * b[3] + d22 * b[6];

          a[1] -= d00 * b[1] + d01 * b[4] + d02 * b[7];
          a[4] -= d10 * b[1] + d11 * b[4] + d12 * b[7];
          a[7] -= d20 * b[1] + d21 * b[4] + d22 * b[7];

          a[2] -= d00 * b[2] + d01 * b[5] + d02 * b[8];
          a[5] -= d10 * b[2] + d11 * b[5] + d12 * b[8];
          a[8] -= d20 * b[2] + d21 * b[5] + d22 * b[8];
        }

        b += 9;
      }

      // Copy the matrix back into the row
      a = &A[9 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d10;
      a[4] = d11;
      a[5] = d12;
      a[6] = d20;
      a[7] = d21;
      a[8] = d22;
    }

    // Invert the diagonal portion of the matrix
    TacsScalar D[9];
    TacsScalar* a = &A[9 * diag[i]];
    D[0] = a[0];
    D[1] = a[1];
    D[2] = a[2];
    D[3] = a[3];
    D[4] = a[4];
    D[5] = a[5];
    D[6] = a[6];
    D[7] = a[7];
    D[8] = a[8];

    int ipiv[3];
    int info = BMatComputeInverse(a, D, ipiv, 3);

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

void BCSRMatFactorLower3(BCSRMatData* data, BCSRMatData* Edata) {
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
      const TacsScalar* d = &A[9 * j];

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar* a = &E[9 * k];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar* b = &E[9 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 9;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          a[0] -= d[0] * b[0] + d[1] * b[3] + d[2] * b[6];
          a[3] -= d[3] * b[0] + d[4] * b[3] + d[5] * b[6];
          a[6] -= d[6] * b[0] + d[7] * b[3] + d[8] * b[6];

          a[1] -= d[0] * b[1] + d[1] * b[4] + d[2] * b[7];
          a[4] -= d[3] * b[1] + d[4] * b[4] + d[5] * b[7];
          a[7] -= d[6] * b[1] + d[7] * b[4] + d[8] * b[7];

          a[2] -= d[0] * b[2] + d[1] * b[5] + d[2] * b[8];
          a[5] -= d[3] * b[2] + d[4] * b[5] + d[5] * b[8];
          a[8] -= d[6] * b[2] + d[7] * b[5] + d[8] * b[8];
        }
        b += 9;
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/

void BCSRMatFactorUpper3(BCSRMatData* data, BCSRMatData* Fdata) {
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

  TacsScalar d00, d01, d02;
  TacsScalar d10, d11, d12;
  TacsScalar d20, d21, d22;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar* a = &F[9 * j];
      const TacsScalar* b = &A[9 * diag[cj]];

      // Multiply d = F[j] * A[diag[cj]]
      d00 = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
      d10 = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
      d20 = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];

      d01 = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
      d11 = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
      d21 = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];

      d02 = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];
      d12 = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];
      d22 = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &F[9 * k];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &A[9 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 9;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          a[0] -= d00 * b[0] + d01 * b[3] + d02 * b[6];
          a[3] -= d10 * b[0] + d11 * b[3] + d12 * b[6];
          a[6] -= d20 * b[0] + d21 * b[3] + d22 * b[6];

          a[1] -= d00 * b[1] + d01 * b[4] + d02 * b[7];
          a[4] -= d10 * b[1] + d11 * b[4] + d12 * b[7];
          a[7] -= d20 * b[1] + d21 * b[4] + d22 * b[7];

          a[2] -= d00 * b[2] + d01 * b[5] + d02 * b[8];
          a[5] -= d10 * b[2] + d11 * b[5] + d12 * b[8];
          a[8] -= d20 * b[2] + d21 * b[5] + d22 * b[8];
        }
        b += 9;
      }

      // Copy over the matrix
      a = &F[9 * j];
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d10;
      a[4] = d11;
      a[5] = d12;
      a[6] = d20;
      a[7] = d21;
      a[8] = d22;
    }
  }
}
