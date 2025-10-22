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
  Implementation of general block code
*/

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern. The entries are over-written and all operations
  are performed in place. This is for an arbitrary block size.
*/

void BCSRMatFactor(BCSRMatData* data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;
  TacsScalar* A = data->A;

  TacsScalar* D = new TacsScalar[bsize * bsize];
  int* ipiv = new int[bsize];

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr, "Error in factorization: no diagonal entry for row %d\n",
              i);
      return;
    }

    // Scan from the first entry in the current row, towards the diagonal
    int row_end = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];
      TacsScalar* a = &A[b2 * j];
      TacsScalar* b = &A[b2 * diag[cj]];

      // D = A[cj] * A[diag[cj]]
      for (int n = 0; n < bsize; n++) {
        for (int m = 0; m < bsize; m++) {
          D[n * bsize + m] = 0.0;
          for (int l = 0; l < bsize; l++) {
            D[n * bsize + m] += a[n * bsize + l] * b[l * bsize + m];
          }
        }
      }

      // Scan through the remainder of the row
      int k = j + 1;

      // The final entry for row: cols[j]
      int end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (int p = diag[cj] + 1; (p < end) && (k < row_end); p++) {
        // Determine where the two rows have the same elements
        while (k < row_end && cols[k] < cols[p]) {
          k++;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < row_end && cols[k] == cols[p]) {
          a = &A[b2 * k];
          b = &A[b2 * p];

          for (int n = 0; n < bsize; n++) {
            int nb = n * bsize;
            for (int m = 0; m < bsize; m++) {
              for (int l = 0; l < bsize; l++) {
                a[nb + m] -= D[nb + l] * b[l * bsize + m];
              }
            }
          }
        }
      }

      // Copy over the matrix
      a = &A[b2 * j];
      for (int n = 0; n < b2; n++) {
        a[n] = D[n];
      }
    }

    // Invert the diagonal matrix component -- Invert( &A[b2*diag[i] )
    TacsScalar* a = &A[b2 * diag[i]];
    for (int n = 0; n < b2; n++) {
      D[n] = a[n];
    }

    int fail = BMatComputeInverse(a, D, ipiv, bsize);
    if (fail != 0) {
      fprintf(stderr, "Failure in factorization of row %d, block row %d\n", i,
              fail);
    }
  }

  delete[] D;
  delete[] ipiv;
}

/*!
  Compute x = L_{B}^{-1} E
*/
void BCSRMatFactorLower(BCSRMatData* data, BCSRMatData* Edata) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  const int bsize = data->bsize;
  const TacsScalar* A = data->A;
  const int b2 = bsize * bsize;

  // Retrieve the data required from the matrix
  const int* erowp = Edata->rowp;
  const int* ecols = Edata->cols;
  TacsScalar* E = Edata->A;

  for (int i = 0; i < nrows; i++) {
    // Scan from the first entry in the current row, towards the diagonal
    int j_end = diag[i];

    for (int j = rowp[i]; j < j_end; j++) {
      int cj = cols[j];
      const TacsScalar* d = &A[b2 * j];

      int k = erowp[i];
      int k_end = erowp[i + 1];

      int p = erowp[cj];
      int p_end = erowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          TacsScalar* a = &E[b2 * k];
          TacsScalar* b = &E[b2 * p];

          for (int n = 0; n < bsize; n++) {
            int nb = n * bsize;
            for (int m = 0; m < bsize; m++) {
              for (int l = 0; l < bsize; l++) {
                a[nb + m] -= d[nb + l] * b[l * bsize + m];
              }
            }
          }
        }
      }
    }
  }
}

/*!
  Compute x = F U_{B}^{-1}
*/
void BCSRMatFactorUpper(BCSRMatData* data, BCSRMatData* Fdata) {
  // Retrieve the data required from the matrix
  const int* rowp = data->rowp;
  const int* cols = data->cols;
  const int* diag = data->diag;
  const int bsize = data->bsize;
  const TacsScalar* A = data->A;
  const int b2 = bsize * bsize;

  // Retrieve the data required from the matrix
  const int nrows_f = Fdata->nrows;
  const int* frowp = Fdata->rowp;
  const int* fcols = Fdata->cols;
  TacsScalar* F = Fdata->A;

  TacsScalar* D = new TacsScalar[bsize * bsize];

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar* a = &F[b2 * j];
      const TacsScalar* b = &A[b2 * diag[cj]];

      // D = A[cj] * A[diag[cj]]
      for (int n = 0; n < bsize; n++) {
        for (int m = 0; m < bsize; m++) {
          D[n * bsize + m] = 0.0;
          for (int l = 0; l < bsize; l++) {
            D[n * bsize + m] += a[n * bsize + l] * b[l * bsize + m];
          }
        }
      }

      int k = j + 1;
      int k_end = frowp[i + 1];

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          a = &F[b2 * k];
          b = &A[b2 * p];

          for (int n = 0; n < bsize; n++) {
            int nb = n * bsize;
            for (int m = 0; m < bsize; m++) {
              for (int l = 0; l < bsize; l++) {
                a[nb + m] -= D[nb + l] * b[l * bsize + m];
              }
            }
          }
        }
      }

      // Copy over the matrix
      a = &F[b2 * j];
      for (int n = 0; n < b2; n++) {
        a[n] = D[n];
      }
    }
  }

  delete[] D;
}
