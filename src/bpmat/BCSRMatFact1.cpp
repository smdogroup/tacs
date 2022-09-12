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
  Block size = 1 implementation.
*/

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern. The entries are over-written and all operations
  are performed in place. This is for an arbitrary block size.
*/
void BCSRMatFactor1(BCSRMatData* data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  TacsScalar* A = data->A;

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr,
              "Error in factorization: no diagonal entry \
for row %d\n",
              i);
      return;
    }

    // Scan from the first entry in the current row, towards the diagonal
    int row_end = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];

      // D = A[j] * A[diag[cj]]
      TacsScalar d11 = A[j] * A[diag[cj]];

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;

      TacsScalar* a = &A[k];
      TacsScalar* b = &A[p];

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          a++;
          k++;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          a[0] -= d11 * b[0];
        }
        b++;
      }

      // Copy over the D matrix
      A[j] = d11;
    }

    // Invert the diagonal matrix component
    if (A[diag[i]] != 0.0) {
      A[diag[i]] = 1.0 / A[diag[i]];
    } else {
      fprintf(stderr, "Failure in factorization of row %d\n", i);
    }
  }
}

/*!
  Compute x = L_{B}^{-1} E
*/
void BCSRMatFactorLower1(BCSRMatData* data, BCSRMatData* Edata) {
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
      const TacsScalar* d = &A[j];

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar* a = &E[k];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar* b = &E[p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a++;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          a[0] -= d[0] * b[0];
        }
        b++;
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/
void BCSRMatFactorUpper1(BCSRMatData* data, BCSRMatData* Fdata) {
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

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar* a = &F[j];
      const TacsScalar* b = &A[diag[cj]];

      // D = A[j] * A[diag[cj]]
      TacsScalar d11 = a[0] * b[0];

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &F[k];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &A[p];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a++;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          a[0] -= d11 * b[0];
        }
        b++;
      }

      // Copy over the matrix
      F[j] = d11;
    }
  }
}
