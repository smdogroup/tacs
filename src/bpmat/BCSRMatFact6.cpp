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
  Block size = 6 code
*/

/*
  Factor the matrix using multiple threads.
*/
void *BCSRMatFactor6_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);

  const int nrows = tdata->mat->nrows;
  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  TacsScalar *A = tdata->mat->A;
  const int group_size = 1;

  TacsScalar d00, d01, d02, d03, d04, d05;
  TacsScalar d10, d11, d12, d13, d14, d15;
  TacsScalar d20, d21, d22, d23, d24, d25;
  TacsScalar d30, d31, d32, d33, d34, d35;
  TacsScalar d40, d41, d42, d43, d44, d45;
  TacsScalar d50, d51, d52, d53, d54, d55;

  while (tdata->num_completed_rows < nrows) {
    int index, row, low, high;
    tdata->apply_lower_sched_job(group_size, &index, &row, &low, &high);

    if (row >= 0) {
      // variable = row
      if (diag[row] < 0) {
        fprintf(stderr, "Error in factorization: no diagonal entry for row %d",
                row);
        pthread_exit(NULL);
        return NULL;
      }

      // Scan from the first entry in the row towards the diagonal
      int kend = rowp[row + 1];
      int jp = rowp[row];
      while (jp < kend && cols[jp] < low) {
        jp++;
      }

      // for j in [low, high)
      for (; (cols[jp] < high) && (cols[jp] < row); jp++) {
        int j = cols[jp];
        TacsScalar *a = &A[36 * jp];
        TacsScalar *b = &A[36 * diag[j]];

        // Multiply d = A[j] *A[diag[cj]]
        TacsScalar b0, b1, b2, b3, b4, b5;

        b0 = b[0];
        b1 = b[6];
        b2 = b[12];
        b3 = b[18];
        b4 = b[24];
        b5 = b[30];
        d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d10 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d20 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d30 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d40 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d50 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[1];
        b1 = b[7];
        b2 = b[13];
        b3 = b[19];
        b4 = b[25];
        b5 = b[31];
        d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d11 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d21 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d31 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d41 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d51 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[2];
        b1 = b[8];
        b2 = b[14];
        b3 = b[20];
        b4 = b[26];
        b5 = b[32];
        d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d12 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d22 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d32 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d42 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d52 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[3];
        b1 = b[9];
        b2 = b[15];
        b3 = b[21];
        b4 = b[27];
        b5 = b[33];
        d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d13 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d23 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d33 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d43 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d53 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[4];
        b1 = b[10];
        b2 = b[16];
        b3 = b[22];
        b4 = b[28];
        b5 = b[34];
        d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d14 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d24 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d34 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d44 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d54 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[5];
        b1 = b[11];
        b2 = b[17];
        b3 = b[23];
        b4 = b[29];
        b5 = b[35];
        d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d15 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d25 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d35 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d45 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d55 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        // Scan through the remainder of the row
        int k = jp + 1;
        int p = diag[j] + 1;
        a = &A[36 * k];
        b = &A[36 * p];

        // The final entry for row: cols[j]
        int pend = rowp[j + 1];

        // Now, scan through row cj starting at the first entry past the
        // diagonal
        for (; (p < pend) && (k < kend); p++) {
          // Determine where the two rows have the same elements
          while (k < kend && cols[k] < cols[p]) {
            k++;
            a += 36;
          }

          // A[k] = A[k] - A[j] * A[p]
          if (k < kend && cols[k] == cols[p]) {
            b0 = b[0];
            b1 = b[6];
            b2 = b[12];
            b3 = b[18];
            b4 = b[24];
            b5 = b[30];
            a[0] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[6] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[12] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[18] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[24] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[30] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[1];
            b1 = b[7];
            b2 = b[13];
            b3 = b[19];
            b4 = b[25];
            b5 = b[31];
            a[1] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[7] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[13] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[19] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[25] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[31] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[2];
            b1 = b[8];
            b2 = b[14];
            b3 = b[20];
            b4 = b[26];
            b5 = b[32];
            a[2] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[8] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[14] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[20] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[26] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[32] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[3];
            b1 = b[9];
            b2 = b[15];
            b3 = b[21];
            b4 = b[27];
            b5 = b[33];
            a[3] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[9] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[15] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[21] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[27] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[33] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[4];
            b1 = b[10];
            b2 = b[16];
            b3 = b[22];
            b4 = b[28];
            b5 = b[34];
            a[4] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[10] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[16] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[22] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[28] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[34] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[5];
            b1 = b[11];
            b2 = b[17];
            b3 = b[23];
            b4 = b[29];
            b5 = b[35];
            a[5] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[11] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[17] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[23] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[29] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[35] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;
          }

          b += 36;
        }

        // Copy the matrix back into the row
        a = &A[36 * jp];
        a[0] = d00;
        a[1] = d01;
        a[2] = d02;
        a[3] = d03;
        a[4] = d04;
        a[5] = d05;
        a[6] = d10;
        a[7] = d11;
        a[8] = d12;
        a[9] = d13;
        a[10] = d14;
        a[11] = d15;
        a[12] = d20;
        a[13] = d21;
        a[14] = d22;
        a[15] = d23;
        a[16] = d24;
        a[17] = d25;
        a[18] = d30;
        a[19] = d31;
        a[20] = d32;
        a[21] = d33;
        a[22] = d34;
        a[23] = d35;
        a[24] = d40;
        a[25] = d41;
        a[26] = d42;
        a[27] = d43;
        a[28] = d44;
        a[29] = d45;
        a[30] = d50;
        a[31] = d51;
        a[32] = d52;
        a[33] = d53;
        a[34] = d54;
        a[35] = d55;
      }

      if (high - 1 == row) {
        // Invert the diagonal portion of the matrix
        TacsScalar D[36];
        TacsScalar *a = &A[36 * diag[row]];
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
        D[25] = a[25];
        D[26] = a[26];
        D[27] = a[27];
        D[28] = a[28];
        D[29] = a[29];
        D[30] = a[30];
        D[31] = a[31];
        D[32] = a[32];
        D[33] = a[33];
        D[34] = a[34];
        D[35] = a[35];

        int ipiv[6];
        int info = BMatComputeInverse(a, D, ipiv, 6);

        if (info > 0) {
          fprintf(stderr,
                  "Error during factorization of diagonal %d in \
block row %d \n",
                  row + 1, info);
        }
      }

      tdata->apply_lower_mark_completed(group_size, index, row, low, high);
    }
  }

  pthread_exit(NULL);
}

/*!
  Compute x = L_{B}^{-1} E
*/
void *BCSRMatFactorLower6_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);

  const int nrows = tdata->mat->nrows;
  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  const TacsScalar *A = tdata->mat->A;
  const int group_size = 1;

  // Retrieve the data required from the matrix
  const int *erowp = tdata->Amat->rowp;
  const int *ecols = tdata->Amat->cols;
  TacsScalar *E = tdata->Amat->A;

  while (tdata->num_completed_rows < nrows) {
    int index, row, low, high;
    tdata->apply_lower_sched_job(group_size, &index, &row, &low, &high);

    if (row >= 0) {
      // Scan from the first entry in the current row, towards the
      // diagonal entry.
      int j_end = diag[row];

      int jp = rowp[row];
      while (jp < j_end && cols[jp] < low) {
        jp++;
      }

      for (; (cols[jp] < high) && (jp < j_end); jp++) {
        int j = cols[jp];
        const TacsScalar *d = &A[36 * jp];

        int k = erowp[row];
        int k_end = erowp[row + 1];
        TacsScalar *a = &E[36 * k];

        int p = erowp[j];
        int p_end = erowp[j + 1];
        TacsScalar *b = &E[36 * p];

        // Now, scan through row cj starting at the first entry past the
        // diagonal
        for (; (p < p_end) && (k < k_end); p++) {
          // Determine where the two rows have the same elements
          while (k < k_end && ecols[k] < ecols[p]) {
            k++;
            a += 36;
          }

          if (k < k_end && ecols[k] == ecols[p]) {
            TacsScalar b0, b1, b2, b3, b4, b5;
            b0 = b[0];
            b1 = b[6];
            b2 = b[12];
            b3 = b[18];
            b4 = b[24];
            b5 = b[30];
            a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[6] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                    d[11] * b5;
            a[12] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[18] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[24] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[30] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;

            b0 = b[1];
            b1 = b[7];
            b2 = b[13];
            b3 = b[19];
            b4 = b[25];
            b5 = b[31];
            a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[7] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                    d[11] * b5;
            a[13] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[19] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[25] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[31] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;

            b0 = b[2];
            b1 = b[8];
            b2 = b[14];
            b3 = b[20];
            b4 = b[26];
            b5 = b[32];
            a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[8] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                    d[11] * b5;
            a[14] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[20] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[26] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[32] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;

            b0 = b[3];
            b1 = b[9];
            b2 = b[15];
            b3 = b[21];
            b4 = b[27];
            b5 = b[33];
            a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[9] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                    d[11] * b5;
            a[15] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[21] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[27] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[33] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;

            b0 = b[4];
            b1 = b[10];
            b2 = b[16];
            b3 = b[22];
            b4 = b[28];
            b5 = b[34];
            a[4] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[10] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 +
                     d[10] * b4 + d[11] * b5;
            a[16] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[22] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[28] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[34] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;

            b0 = b[5];
            b1 = b[11];
            b2 = b[17];
            b3 = b[23];
            b4 = b[29];
            b5 = b[35];
            a[5] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5;
            a[11] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 +
                     d[10] * b4 + d[11] * b5;
            a[17] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                     d[16] * b4 + d[17] * b5;
            a[23] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                     d[22] * b4 + d[23] * b5;
            a[29] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5;
            a[35] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                     d[34] * b4 + d[35] * b5;
          }
          b += 36;
        }
      }

      tdata->apply_lower_mark_completed(group_size, index, row, low, high);
    }
  }

  pthread_exit(NULL);
}

/*!
  Compute x = F U_{B}^{-1}
*/
void *BCSRMatFactorUpper6_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);

  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  const TacsScalar *A = tdata->mat->A;
  const int group_size = 1;

  // Retrieve the data required from the matrix
  const int nrows_f = tdata->Amat->nrows;
  const int *frowp = tdata->Amat->rowp;
  const int *fcols = tdata->Amat->cols;
  TacsScalar *F = tdata->Amat->A;

  TacsScalar d00, d01, d02, d03, d04, d05;
  TacsScalar d10, d11, d12, d13, d14, d15;
  TacsScalar d20, d21, d22, d23, d24, d25;
  TacsScalar d30, d31, d32, d33, d34, d35;
  TacsScalar d40, d41, d42, d43, d44, d45;
  TacsScalar d50, d51, d52, d53, d54, d55;

  while (tdata->num_completed_rows < nrows_f) {
    // Here, low and high are the low/high values of the rows to add into
    // the current row
    int row;
    tdata->mat_mult_sched_job_size(group_size, &row, nrows_f);

    if (row >= 0) {
      int j_end = frowp[row + 1];
      int jp = frowp[row];

      for (; jp < j_end; jp++) {
        int j = fcols[jp];
        TacsScalar *a = &F[36 * jp];
        const TacsScalar *b = &A[36 * diag[j]];

        // Multiply d = F[j] *A[diag[cj]]
        TacsScalar b0, b1, b2, b3, b4, b5;

        b0 = b[0];
        b1 = b[6];
        b2 = b[12];
        b3 = b[18];
        b4 = b[24];
        b5 = b[30];
        d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d10 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d20 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d30 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d40 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d50 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[1];
        b1 = b[7];
        b2 = b[13];
        b3 = b[19];
        b4 = b[25];
        b5 = b[31];
        d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d11 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d21 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d31 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d41 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d51 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[2];
        b1 = b[8];
        b2 = b[14];
        b3 = b[20];
        b4 = b[26];
        b5 = b[32];
        d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d12 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d22 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d32 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d42 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d52 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[3];
        b1 = b[9];
        b2 = b[15];
        b3 = b[21];
        b4 = b[27];
        b5 = b[33];
        d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d13 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d23 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d33 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d43 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d53 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[4];
        b1 = b[10];
        b2 = b[16];
        b3 = b[22];
        b4 = b[28];
        b5 = b[34];
        d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d14 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d24 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d34 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d44 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d54 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        b0 = b[5];
        b1 = b[11];
        b2 = b[17];
        b3 = b[23];
        b4 = b[29];
        b5 = b[35];
        d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5;
        d15 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
              a[11] * b5;
        d25 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
              a[17] * b5;
        d35 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
              a[23] * b5;
        d45 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5;
        d55 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
              a[35] * b5;

        int k = jp + 1;
        int k_end = frowp[row + 1];
        a = &F[36 * k];

        int p = diag[j] + 1;
        int p_end = rowp[j + 1];
        b = &A[36 * p];

        // Now, scan through row j starting at the first entry past the diagonal
        for (; (p < p_end) && (k < k_end); p++) {
          // Determine where the two rows have the same elements
          while (k < k_end && fcols[k] < cols[p]) {
            k++;
            a += 36;
          }

          if (k < k_end && fcols[k] == cols[p]) {
            b0 = b[0];
            b1 = b[6];
            b2 = b[12];
            b3 = b[18];
            b4 = b[24];
            b5 = b[30];
            a[0] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[6] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[12] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[18] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[24] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[30] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[1];
            b1 = b[7];
            b2 = b[13];
            b3 = b[19];
            b4 = b[25];
            b5 = b[31];
            a[1] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[7] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[13] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[19] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[25] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[31] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[2];
            b1 = b[8];
            b2 = b[14];
            b3 = b[20];
            b4 = b[26];
            b5 = b[32];
            a[2] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[8] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[14] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[20] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[26] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[32] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[3];
            b1 = b[9];
            b2 = b[15];
            b3 = b[21];
            b4 = b[27];
            b5 = b[33];
            a[3] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[9] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[15] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[21] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[27] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[33] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[4];
            b1 = b[10];
            b2 = b[16];
            b3 = b[22];
            b4 = b[28];
            b5 = b[34];
            a[4] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[10] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[16] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[22] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[28] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[34] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

            b0 = b[5];
            b1 = b[11];
            b2 = b[17];
            b3 = b[23];
            b4 = b[29];
            b5 = b[35];
            a[5] -=
                d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
            a[11] -=
                d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
            a[17] -=
                d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
            a[23] -=
                d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
            a[29] -=
                d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
            a[35] -=
                d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;
          }
          b += 36;
        }

        // Copy over the matrix
        a = &F[36 * jp];
        a[0] = d00;
        a[1] = d01;
        a[2] = d02;
        a[3] = d03;
        a[4] = d04;
        a[5] = d05;
        a[6] = d10;
        a[7] = d11;
        a[8] = d12;
        a[9] = d13;
        a[10] = d14;
        a[11] = d15;
        a[12] = d20;
        a[13] = d21;
        a[14] = d22;
        a[15] = d23;
        a[16] = d24;
        a[17] = d25;
        a[18] = d30;
        a[19] = d31;
        a[20] = d32;
        a[21] = d33;
        a[22] = d34;
        a[23] = d35;
        a[24] = d40;
        a[25] = d41;
        a[26] = d42;
        a[27] = d43;
        a[28] = d44;
        a[29] = d45;
        a[30] = d50;
        a[31] = d51;
        a[32] = d52;
        a[33] = d53;
        a[34] = d54;
        a[35] = d55;
      }
    }
  }

  pthread_exit(NULL);
}

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern.  The entries are over-written, all operations are
  performed in place.
*/
void BCSRMatFactor6(BCSRMatData *data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  TacsScalar d00, d01, d02, d03, d04, d05;
  TacsScalar d10, d11, d12, d13, d14, d15;
  TacsScalar d20, d21, d22, d23, d24, d25;
  TacsScalar d30, d31, d32, d33, d34, d35;
  TacsScalar d40, d41, d42, d43, d44, d45;
  TacsScalar d50, d51, d52, d53, d54, d55;

  for (int i = 0; i < nrows; i++) {
    // variable = i
    if (diag[i] < 0) {
      fprintf(stderr, "Error in factorization: no diagonal entry for row %d",
              i);
      return;
    }

    // Scan from the first entry in the current row, towards the diagonal
    int kend = rowp[i + 1];

    for (int j = rowp[i]; cols[j] < i; j++) {
      int cj = cols[j];
      TacsScalar *a = &(data->A[36 * j]);
      TacsScalar *b = &(data->A[36 * diag[cj]]);

      // Multiply d = A[j] * A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4, b5;

      b0 = b[0];
      b1 = b[6];
      b2 = b[12];
      b3 = b[18];
      b4 = b[24];
      b5 = b[30];
      d00 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d10 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d20 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d30 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d40 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d50 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[1];
      b1 = b[7];
      b2 = b[13];
      b3 = b[19];
      b4 = b[25];
      b5 = b[31];
      d01 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d11 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d21 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d31 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d41 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d51 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[2];
      b1 = b[8];
      b2 = b[14];
      b3 = b[20];
      b4 = b[26];
      b5 = b[32];
      d02 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d12 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d22 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d32 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d42 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d52 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[3];
      b1 = b[9];
      b2 = b[15];
      b3 = b[21];
      b4 = b[27];
      b5 = b[33];
      d03 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d13 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d23 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d33 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d43 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d53 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[4];
      b1 = b[10];
      b2 = b[16];
      b3 = b[22];
      b4 = b[28];
      b5 = b[34];
      d04 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d14 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d24 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d34 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d44 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d54 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[5];
      b1 = b[11];
      b2 = b[17];
      b3 = b[23];
      b4 = b[29];
      b5 = b[35];
      d05 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d15 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d25 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d35 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d45 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d55 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;
      a = &(data->A[36 * k]);
      b = &(data->A[36 * p]);

      // The final entry for row: cols[j]
      int pend = rowp[cj + 1];

      // Keep track of the number of block matrix products
      int nz = 0;

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < pend) && (k < kend); p++) {
        // Determine where the two rows have the same elements
        while (k < kend && cols[k] < cols[p]) {
          k++;
          a += 36;
        }

        // A[k] = A[k] - A[j] * A[p]
        if (k < kend && cols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[6];
          b2 = b[12];
          b3 = b[18];
          b4 = b[24];
          b5 = b[30];
          a[0] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[6] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[12] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[18] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[24] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[30] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[1];
          b1 = b[7];
          b2 = b[13];
          b3 = b[19];
          b4 = b[25];
          b5 = b[31];
          a[1] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[7] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[13] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[19] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[25] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[31] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[2];
          b1 = b[8];
          b2 = b[14];
          b3 = b[20];
          b4 = b[26];
          b5 = b[32];
          a[2] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[8] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[14] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[20] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[26] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[32] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[3];
          b1 = b[9];
          b2 = b[15];
          b3 = b[21];
          b4 = b[27];
          b5 = b[33];
          a[3] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[9] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[15] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[21] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[27] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[33] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[4];
          b1 = b[10];
          b2 = b[16];
          b3 = b[22];
          b4 = b[28];
          b5 = b[34];
          a[4] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[10] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[16] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[22] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[28] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[34] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[5];
          b1 = b[11];
          b2 = b[17];
          b3 = b[23];
          b4 = b[29];
          b5 = b[35];
          a[5] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[11] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[17] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[23] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[29] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[35] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          nz++;
        }

        b += 36;
      }

      TacsAddFlops(2 * 36 * 6 * nz + 11 * 36);

      // Copy the matrix back into the row
      a = &(data->A[36 * j]);
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d05;
      a[6] = d10;
      a[7] = d11;
      a[8] = d12;
      a[9] = d13;
      a[10] = d14;
      a[11] = d15;
      a[12] = d20;
      a[13] = d21;
      a[14] = d22;
      a[15] = d23;
      a[16] = d24;
      a[17] = d25;
      a[18] = d30;
      a[19] = d31;
      a[20] = d32;
      a[21] = d33;
      a[22] = d34;
      a[23] = d35;
      a[24] = d40;
      a[25] = d41;
      a[26] = d42;
      a[27] = d43;
      a[28] = d44;
      a[29] = d45;
      a[30] = d50;
      a[31] = d51;
      a[32] = d52;
      a[33] = d53;
      a[34] = d54;
      a[35] = d55;
    }

    // Invert the diagonal portion of the matrix
    TacsScalar D[36];
    TacsScalar *a = &(data->A[36 * diag[i]]);
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
    D[25] = a[25];
    D[26] = a[26];
    D[27] = a[27];
    D[28] = a[28];
    D[29] = a[29];
    D[30] = a[30];
    D[31] = a[31];
    D[32] = a[32];
    D[33] = a[33];
    D[34] = a[34];
    D[35] = a[35];

    int ipiv[6];
    int info = BMatComputeInverse(a, D, ipiv, 6);

    if (info > 0) {
      fprintf(stderr,
              "Error during factorization of diagonal %d in \
block row %d \n",
              i + 1, info);
    }
  }

  // Add flops from the diagonal inversion
  TacsAddFlops(1.333333 * 6 * 6 * 6 * nrows);
}

/*!
  Compute x = L_{B}^{-1} E
*/
void BCSRMatFactorLower6(BCSRMatData *data, BCSRMatData *Edata) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  // Retrieve the data required from the matrix
  const int *erowp = Edata->rowp;
  const int *ecols = Edata->cols;

  // Keep track of the number of block matrix products
  int nz = 0;

  for (int i = 0; i < nrows; i++) {
    // Scan from the first entry in the current row, towards the
    // diagonal entry.
    int j_end = diag[i];

    for (int j = rowp[i]; j < j_end; j++) {
      int cj = cols[j];
      TacsScalar *d = &(data->A[36 * j]);

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar *a = &(Edata->A[36 * k]);

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar *b = &(Edata->A[36 * p]);

      // Now, scan through row cj starting at the first entry past the
      // diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 36;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          TacsScalar b0, b1, b2, b3, b4, b5;
          b0 = b[0];
          b1 = b[6];
          b2 = b[12];
          b3 = b[18];
          b4 = b[24];
          b5 = b[30];
          a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[6] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                  d[11] * b5;
          a[12] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[18] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[24] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[30] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          b0 = b[1];
          b1 = b[7];
          b2 = b[13];
          b3 = b[19];
          b4 = b[25];
          b5 = b[31];
          a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[7] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                  d[11] * b5;
          a[13] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[19] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[25] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[31] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          b0 = b[2];
          b1 = b[8];
          b2 = b[14];
          b3 = b[20];
          b4 = b[26];
          b5 = b[32];
          a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[8] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                  d[11] * b5;
          a[14] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[20] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[26] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[32] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          b0 = b[3];
          b1 = b[9];
          b2 = b[15];
          b3 = b[21];
          b4 = b[27];
          b5 = b[33];
          a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[9] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                  d[11] * b5;
          a[15] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[21] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[27] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[33] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          b0 = b[4];
          b1 = b[10];
          b2 = b[16];
          b3 = b[22];
          b4 = b[28];
          b5 = b[34];
          a[4] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[10] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                   d[11] * b5;
          a[16] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[22] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[28] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[34] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          b0 = b[5];
          b1 = b[11];
          b2 = b[17];
          b3 = b[23];
          b4 = b[29];
          b5 = b[35];
          a[5] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5;
          a[11] -= d[6] * b0 + d[7] * b1 + d[8] * b2 + d[9] * b3 + d[10] * b4 +
                   d[11] * b5;
          a[17] -= d[12] * b0 + d[13] * b1 + d[14] * b2 + d[15] * b3 +
                   d[16] * b4 + d[17] * b5;
          a[23] -= d[18] * b0 + d[19] * b1 + d[20] * b2 + d[21] * b3 +
                   d[22] * b4 + d[23] * b5;
          a[29] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5;
          a[35] -= d[30] * b0 + d[31] * b1 + d[32] * b2 + d[33] * b3 +
                   d[34] * b4 + d[35] * b5;

          nz++;
        }
        b += 36;
      }
    }
  }

  TacsAddFlops(2 * 36 * 6 * nz);
}

/*!
  Compute x = F U_{B}^{-1}
*/
void BCSRMatFactorUpper6(BCSRMatData *data, BCSRMatData *Fdata) {
  // Retrieve the data required from the matrix
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  // Retrieve the data required from the matrix
  const int nrows_f = Fdata->nrows;
  const int *frowp = Fdata->rowp;
  const int *fcols = Fdata->cols;

  TacsScalar d00, d01, d02, d03, d04, d05;
  TacsScalar d10, d11, d12, d13, d14, d15;
  TacsScalar d20, d21, d22, d23, d24, d25;
  TacsScalar d30, d31, d32, d33, d34, d35;
  TacsScalar d40, d41, d42, d43, d44, d45;
  TacsScalar d50, d51, d52, d53, d54, d55;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar *a = &(Fdata->A[36 * j]);
      const TacsScalar *b = &(data->A[36 * diag[cj]]);

      // Multiply d = F[j] * A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4, b5;

      b0 = b[0];
      b1 = b[6];
      b2 = b[12];
      b3 = b[18];
      b4 = b[24];
      b5 = b[30];
      d00 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d10 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d20 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d30 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d40 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d50 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[1];
      b1 = b[7];
      b2 = b[13];
      b3 = b[19];
      b4 = b[25];
      b5 = b[31];
      d01 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d11 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d21 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d31 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d41 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d51 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[2];
      b1 = b[8];
      b2 = b[14];
      b3 = b[20];
      b4 = b[26];
      b5 = b[32];
      d02 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d12 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d22 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d32 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d42 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d52 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[3];
      b1 = b[9];
      b2 = b[15];
      b3 = b[21];
      b4 = b[27];
      b5 = b[33];
      d03 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d13 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d23 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d33 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d43 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d53 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[4];
      b1 = b[10];
      b2 = b[16];
      b3 = b[22];
      b4 = b[28];
      b5 = b[34];
      d04 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d14 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d24 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d34 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d44 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d54 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      b0 = b[5];
      b1 = b[11];
      b2 = b[17];
      b3 = b[23];
      b4 = b[29];
      b5 = b[35];
      d05 =
          a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 + a[5] * b5;
      d15 = a[6] * b0 + a[7] * b1 + a[8] * b2 + a[9] * b3 + a[10] * b4 +
            a[11] * b5;
      d25 = a[12] * b0 + a[13] * b1 + a[14] * b2 + a[15] * b3 + a[16] * b4 +
            a[17] * b5;
      d35 = a[18] * b0 + a[19] * b1 + a[20] * b2 + a[21] * b3 + a[22] * b4 +
            a[23] * b5;
      d45 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5;
      d55 = a[30] * b0 + a[31] * b1 + a[32] * b2 + a[33] * b3 + a[34] * b4 +
            a[35] * b5;

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &(Fdata->A[36 * k]);

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &(data->A[36 * p]);

      // Keep track of the number of block matrix-matrix products
      int nz = 0;

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 36;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[6];
          b2 = b[12];
          b3 = b[18];
          b4 = b[24];
          b5 = b[30];
          a[0] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[6] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[12] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[18] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[24] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[30] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[1];
          b1 = b[7];
          b2 = b[13];
          b3 = b[19];
          b4 = b[25];
          b5 = b[31];
          a[1] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[7] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[13] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[19] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[25] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[31] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[2];
          b1 = b[8];
          b2 = b[14];
          b3 = b[20];
          b4 = b[26];
          b5 = b[32];
          a[2] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[8] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[14] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[20] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[26] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[32] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[3];
          b1 = b[9];
          b2 = b[15];
          b3 = b[21];
          b4 = b[27];
          b5 = b[33];
          a[3] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[9] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[15] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[21] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[27] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[33] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[4];
          b1 = b[10];
          b2 = b[16];
          b3 = b[22];
          b4 = b[28];
          b5 = b[34];
          a[4] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[10] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[16] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[22] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[28] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[34] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          b0 = b[5];
          b1 = b[11];
          b2 = b[17];
          b3 = b[23];
          b4 = b[29];
          b5 = b[35];
          a[5] -=
              d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 + d05 * b5;
          a[11] -=
              d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 + d15 * b5;
          a[17] -=
              d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 + d25 * b5;
          a[23] -=
              d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 + d35 * b5;
          a[29] -=
              d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 + d45 * b5;
          a[35] -=
              d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 + d55 * b5;

          nz++;
        }
        b += 36;
      }

      TacsAddFlops(2 * 36 * 6 * nz + 11 * 36);

      // Copy over the matrix
      a = &(Fdata->A[36 * j]);
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d05;
      a[6] = d10;
      a[7] = d11;
      a[8] = d12;
      a[9] = d13;
      a[10] = d14;
      a[11] = d15;
      a[12] = d20;
      a[13] = d21;
      a[14] = d22;
      a[15] = d23;
      a[16] = d24;
      a[17] = d25;
      a[18] = d30;
      a[19] = d31;
      a[20] = d32;
      a[21] = d33;
      a[22] = d34;
      a[23] = d35;
      a[24] = d40;
      a[25] = d41;
      a[26] = d42;
      a[27] = d43;
      a[28] = d44;
      a[29] = d45;
      a[30] = d50;
      a[31] = d51;
      a[32] = d52;
      a[33] = d53;
      a[34] = d54;
      a[35] = d55;
    }
  }
}
