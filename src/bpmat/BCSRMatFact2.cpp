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
  Block size = 2 implementation.
*/

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern. The entries are over-written and all operations
  are performed in place. This is for an arbitrary block size.
*/

void BCSRMatFactor2(BCSRMatData *data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  TacsScalar *A = data->A;

  TacsScalar d11, d12;
  TacsScalar d21, d22;

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr,
              "Error in factorization: no diagonal entry "
              "for row %d\n",
              i);
      return;
    }

    // Scan from the first entry in the current row, towards the diagonal
    int row_end = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];
      TacsScalar *a = &A[4 * j];
      TacsScalar *b = &A[4 * diag[cj]];

      // D = A[j] * A[diag[cj]]
      d11 = a[0] * b[0] + a[1] * b[2];
      d21 = a[2] * b[0] + a[3] * b[2];

      d12 = a[0] * b[1] + a[1] * b[3];
      d22 = a[2] * b[1] + a[3] * b[3];

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;

      a = &A[4 * k];
      b = &A[4 * p];

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          a += 4;
          k++;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          a[0] -= d11 * b[0] + d12 * b[2];
          a[1] -= d11 * b[1] + d12 * b[3];

          a[2] -= d21 * b[0] + d22 * b[2];
          a[3] -= d21 * b[1] + d22 * b[3];
        }
        b += 4;
      }

      // Copy over the D matrix
      a = &A[4 * j];
      a[0] = d11;
      a[1] = d12;
      a[2] = d21;
      a[3] = d22;
    }

    // Invert the diagonal matrix component -- Invert( &A[4*diag[i] )
    TacsScalar *a = &A[4 * diag[i]];
    d11 = a[0];
    d12 = a[1];
    d21 = a[2];
    d22 = a[3];

    TacsScalar det = (d11 * d22 - d12 * d21);
    if (det == 0.0) {
      fprintf(stderr, "Failure in factorization of row %d\n", i);
    }
    det = 1.0 / det;

    a[0] = det * d22;
    a[1] = -det * d12;
    a[2] = -det * d21;
    a[3] = det * d11;
  }
}

/*!
  Compute x = L_{B}^{-1} E
*/

void BCSRMatFactorLower2(BCSRMatData *data, BCSRMatData *Edata) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  // Retrieve the data required from the matrix
  const int *erowp = Edata->rowp;
  const int *ecols = Edata->cols;
  TacsScalar *E = Edata->A;

  for (int i = 0; i < nrows; i++) {
    // Scan from the first entry in the current row, towards the diagonal
    int j_end = diag[i];

    for (int j = rowp[i]; j < j_end; j++) {
      int cj = cols[j];
      const TacsScalar *d = &A[4 * j];

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar *a = &E[4 * k];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar *b = &E[4 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 4;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          a[0] -= d[0] * b[0] + d[1] * b[2];
          a[1] -= d[0] * b[1] + d[1] * b[3];

          a[2] -= d[2] * b[0] + d[3] * b[2];
          a[3] -= d[2] * b[1] + d[3] * b[3];
        }
        b += 4;
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/

void BCSRMatFactorUpper2(BCSRMatData *data, BCSRMatData *Fdata) {
  // Retrieve the data required from the matrix
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  // Retrieve the data required from the matrix
  const int nrows_f = Fdata->nrows;
  const int *frowp = Fdata->rowp;
  const int *fcols = Fdata->cols;
  TacsScalar *F = Fdata->A;

  TacsScalar d11, d12;
  TacsScalar d21, d22;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar *a = &F[4 * j];
      const TacsScalar *b = &A[4 * diag[cj]];

      // D = A[j] * A[diag[cj]]
      d11 = a[0] * b[0] + a[1] * b[2];
      d21 = a[2] * b[0] + a[3] * b[2];

      d12 = a[0] * b[1] + a[1] * b[3];
      d22 = a[2] * b[1] + a[3] * b[3];

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &F[4 * k];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &A[4 * p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 4;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          a[0] -= d11 * b[0] + d12 * b[2];
          a[1] -= d11 * b[1] + d12 * b[3];

          a[2] -= d21 * b[0] + d22 * b[2];
          a[3] -= d21 * b[1] + d22 * b[3];
        }
        b += 4;
      }

      // Copy over the matrix
      a = &F[4 * j];
      a[0] = d11;
      a[1] = d12;
      a[2] = d21;
      a[3] = d22;
    }
  }
}
