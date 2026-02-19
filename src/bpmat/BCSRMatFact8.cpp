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
  Block size = 8 code
*/

/*
  Factor the matrix using multiple threads.
*/
void *BCSRMatFactor8_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);

  const int nrows = tdata->mat->nrows;
  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  TacsScalar *A = tdata->mat->A;
  const int group_size = 1;

  TacsScalar d00, d01, d02, d03, d04, d05, d06, d07;
  TacsScalar d10, d11, d12, d13, d14, d15, d16, d17;
  TacsScalar d20, d21, d22, d23, d24, d25, d26, d27;
  TacsScalar d30, d31, d32, d33, d34, d35, d36, d37;
  TacsScalar d40, d41, d42, d43, d44, d45, d46, d47;
  TacsScalar d50, d51, d52, d53, d54, d55, d56, d57;
  TacsScalar d60, d61, d62, d63, d64, d65, d66, d67;
  TacsScalar d70, d71, d72, d73, d74, d75, d76, d77;

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
        TacsScalar *a = &A[64 * jp];
        TacsScalar *b = &A[64 * diag[j]];

        // Multiply d = A[j] *A[diag[cj]]
        TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

        b0 = b[0];
        b1 = b[8];
        b2 = b[16];
        b3 = b[24];
        b4 = b[32];
        b5 = b[40];
        b6 = b[48];
        b7 = b[56];
        d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d10 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d20 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d30 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d40 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d50 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d60 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d70 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[1];
        b1 = b[9];
        b2 = b[17];
        b3 = b[25];
        b4 = b[33];
        b5 = b[41];
        b6 = b[49];
        b7 = b[57];
        d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d11 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d21 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d31 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d41 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d51 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d61 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d71 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[2];
        b1 = b[10];
        b2 = b[18];
        b3 = b[26];
        b4 = b[34];
        b5 = b[42];
        b6 = b[50];
        b7 = b[58];
        d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d12 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d22 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d32 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d42 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d52 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d62 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d72 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[3];
        b1 = b[11];
        b2 = b[19];
        b3 = b[27];
        b4 = b[35];
        b5 = b[43];
        b6 = b[51];
        b7 = b[59];
        d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d13 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d23 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d33 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d43 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d53 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d63 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d73 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[4];
        b1 = b[12];
        b2 = b[20];
        b3 = b[28];
        b4 = b[36];
        b5 = b[44];
        b6 = b[52];
        b7 = b[60];
        d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d14 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d24 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d34 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d44 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d54 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d64 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d74 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[5];
        b1 = b[13];
        b2 = b[21];
        b3 = b[29];
        b4 = b[37];
        b5 = b[45];
        b6 = b[53];
        b7 = b[61];
        d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d15 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d25 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d35 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d45 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d55 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d65 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d75 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[6];
        b1 = b[14];
        b2 = b[22];
        b3 = b[30];
        b4 = b[38];
        b5 = b[46];
        b6 = b[54];
        b7 = b[62];
        d06 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d16 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d26 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d36 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d46 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d56 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d66 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d76 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[7];
        b1 = b[15];
        b2 = b[23];
        b3 = b[31];
        b4 = b[39];
        b5 = b[47];
        b6 = b[55];
        b7 = b[63];
        d07 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d17 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d27 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d37 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d47 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d57 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d67 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d77 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        // Scan through the remainder of the row
        int k = jp + 1;
        int p = diag[j] + 1;
        a = &A[64 * k];
        b = &A[64 * p];

        // The final entry for row: cols[j]
        int pend = rowp[j + 1];

        // Now, scan through row cj starting at the first entry past the
        // diagonal
        for (; (p < pend) && (k < kend); p++) {
          // Determine where the two rows have the same elements
          while (k < kend && cols[k] < cols[p]) {
            k++;
            a += 64;
          }

          // A[k] = A[k] - A[j] * A[p]
          if (k < kend && cols[k] == cols[p]) {
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40], b6 = b[48], b7 = b[56];
            a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                    d15 * b5 + d16 * b6 + d17 * b7;
            a[16] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[24] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[32] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[40] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[48] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[56] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41], b6 = b[49], b7 = b[57];
            a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                    d15 * b5 + d16 * b6 + d17 * b7;
            a[17] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[25] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[33] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[41] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[49] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[57] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42], b6 = b[50], b7 = b[58];
            a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[10] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[18] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[26] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[34] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[42] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[50] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[58] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43], b6 = b[51], b7 = b[59];
            a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[11] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[19] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[27] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[35] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[43] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[51] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[59] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44], b6 = b[52], b7 = b[60];
            a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[12] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[20] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[28] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[36] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[44] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[52] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[60] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45], b6 = b[53], b7 = b[61];
            a[5] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[13] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[21] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[29] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[37] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[45] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[53] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[61] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46], b6 = b[54], b7 = b[62];
            a[6] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[14] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[22] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[30] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[38] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[46] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[54] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[62] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47], b6 = b[55], b7 = b[63];
            a[7] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[15] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[23] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[31] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[39] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[47] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[55] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[63] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;
          }

          b += 64;
        }

        // Copy the matrix back into the row
        a = &A[64 * jp];
        a[0] = d00;
        a[1] = d01;
        a[2] = d02;
        a[3] = d03;
        a[4] = d04;
        a[5] = d05;
        a[6] = d06;
        a[7] = d07;
        a[8] = d10;
        a[9] = d11;
        a[10] = d12;
        a[11] = d13;
        a[12] = d14;
        a[13] = d15;
        a[14] = d16;
        a[15] = d17;
        a[16] = d20;
        a[17] = d21;
        a[18] = d22;
        a[19] = d23;
        a[20] = d24;
        a[21] = d25;
        a[22] = d26;
        a[23] = d27;
        a[24] = d30;
        a[25] = d31;
        a[26] = d32;
        a[27] = d33;
        a[28] = d34;
        a[29] = d35;
        a[30] = d36;
        a[31] = d37;
        a[32] = d40;
        a[33] = d41;
        a[34] = d42;
        a[35] = d43;
        a[36] = d44;
        a[37] = d45;
        a[38] = d46;
        a[39] = d47;
        a[40] = d50;
        a[41] = d51;
        a[42] = d52;
        a[43] = d53;
        a[44] = d54;
        a[45] = d55;
        a[46] = d56;
        a[47] = d57;
        a[48] = d60;
        a[49] = d61;
        a[50] = d62;
        a[51] = d63;
        a[52] = d64;
        a[53] = d65;
        a[54] = d66;
        a[55] = d67;
        a[56] = d70;
        a[57] = d71;
        a[58] = d72;
        a[59] = d73;
        a[60] = d74;
        a[61] = d75;
        a[62] = d76;
        a[63] = d77;
      }

      if (high - 1 == row) {
        // Invert the diagonal portion of the matrix
        TacsScalar D[64];
        TacsScalar *a = &A[64 * diag[row]];
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
        D[36] = a[36];
        D[37] = a[37];
        D[38] = a[38];
        D[39] = a[39];
        D[40] = a[40];
        D[41] = a[41];
        D[42] = a[42];
        D[43] = a[43];
        D[44] = a[44];
        D[45] = a[45];
        D[46] = a[46];
        D[47] = a[47];
        D[48] = a[48];
        D[49] = a[49];
        D[50] = a[50];
        D[51] = a[51];
        D[52] = a[52];
        D[53] = a[53];
        D[54] = a[54];
        D[55] = a[55];
        D[56] = a[56];
        D[57] = a[57];
        D[58] = a[58];
        D[59] = a[59];
        D[60] = a[60];
        D[61] = a[61];
        D[62] = a[62];
        D[63] = a[63];

        int ipiv[8];
        int info = BMatComputeInverse(a, D, ipiv, 8);

        if (info > 0) {
          fprintf(
              stderr,
              "Error during factorization of diagonal %d in block row %d \n",
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
void *BCSRMatFactorLower8_thread(void *t) {
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
        const TacsScalar *d = &A[64 * jp];

        int k = erowp[row];
        int k_end = erowp[row + 1];
        TacsScalar *a = &E[64 * k];

        int p = erowp[j];
        int p_end = erowp[j + 1];
        TacsScalar *b = &E[64 * p];

        // Now, scan through row cj starting at the first entry past the
        // diagonal
        for (; (p < p_end) && (k < k_end); p++) {
          // Determine where the two rows have the same elements
          while (k < k_end && ecols[k] < ecols[p]) {
            k++;
            a += 64;
          }

          if (k < k_end && ecols[k] == ecols[p]) {
            TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40];
            b6 = b[48];
            b7 = b[56];
            a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[8] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                    d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[16] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[24] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[32] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[40] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[48] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[56] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41];
            b6 = b[49];
            b7 = b[57];
            a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[9] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                    d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[17] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[25] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[33] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[41] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[49] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[57] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42];
            b6 = b[50];
            b7 = b[58];
            a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[10] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[18] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[26] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[34] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[42] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[50] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[58] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43];
            b6 = b[51];
            b7 = b[59];
            a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[11] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[19] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[27] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[35] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[43] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[51] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[59] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44];
            b6 = b[52];
            b7 = b[60];
            a[4] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[12] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[20] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[28] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[36] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[44] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[52] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[60] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45];
            b6 = b[53];
            b7 = b[61];
            a[5] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[13] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[21] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[29] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[37] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[45] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[53] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[61] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46];
            b6 = b[54];
            b7 = b[62];
            a[6] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[14] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[22] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[30] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[38] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[46] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[54] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[62] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47];
            b6 = b[55];
            b7 = b[63];
            a[7] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                    d[5] * b5 + d[6] * b6 + d[7] * b7;
            a[15] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                     d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
            a[23] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                     d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
            a[31] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                     d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
            a[39] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                     d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
            a[47] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                     d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
            a[55] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                     d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
            a[63] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                     d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;
          }
          b += 64;
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
void *BCSRMatFactorUpper8_thread(void *t) {
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

  TacsScalar d00, d01, d02, d03, d04, d05, d06, d07;
  TacsScalar d10, d11, d12, d13, d14, d15, d16, d17;
  TacsScalar d20, d21, d22, d23, d24, d25, d26, d27;
  TacsScalar d30, d31, d32, d33, d34, d35, d36, d37;
  TacsScalar d40, d41, d42, d43, d44, d45, d46, d47;
  TacsScalar d50, d51, d52, d53, d54, d55, d56, d57;
  TacsScalar d60, d61, d62, d63, d64, d65, d66, d67;
  TacsScalar d70, d71, d72, d73, d74, d75, d76, d77;

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
        TacsScalar *a = &F[64 * jp];
        const TacsScalar *b = &A[64 * diag[j]];

        // Multiply d = F[j] *A[diag[cj]]
        TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

        b0 = b[0];
        b1 = b[8];
        b2 = b[16];
        b3 = b[24];
        b4 = b[32];
        b5 = b[40];
        b6 = b[48];
        b7 = b[56];
        d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d10 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d20 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d30 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d40 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d50 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d60 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d70 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[1];
        b1 = b[9];
        b2 = b[17];
        b3 = b[25];
        b4 = b[33];
        b5 = b[41];
        b6 = b[49];
        b7 = b[57];
        d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d11 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d21 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d31 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d41 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d51 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d61 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d71 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[2];
        b1 = b[10];
        b2 = b[18];
        b3 = b[26];
        b4 = b[34];
        b5 = b[42];
        b6 = b[50];
        b7 = b[58];
        d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d12 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d22 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d32 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d42 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d52 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d62 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d72 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[3];
        b1 = b[11];
        b2 = b[19];
        b3 = b[27];
        b4 = b[35];
        b5 = b[43];
        b6 = b[51];
        b7 = b[59];
        d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d13 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d23 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d33 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d43 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d53 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d63 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d73 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[4];
        b1 = b[12];
        b2 = b[20];
        b3 = b[28];
        b4 = b[36];
        b5 = b[44];
        b6 = b[52];
        b7 = b[60];
        d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d14 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d24 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d34 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d44 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d54 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d64 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d74 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[5];
        b1 = b[13];
        b2 = b[21];
        b3 = b[29];
        b4 = b[37];
        b5 = b[45];
        b6 = b[53];
        b7 = b[61];
        d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d15 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d25 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d35 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d45 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d55 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d65 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d75 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[6];
        b1 = b[14];
        b2 = b[22];
        b3 = b[30];
        b4 = b[38];
        b5 = b[46];
        b6 = b[54];
        b7 = b[62];
        d06 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d16 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d26 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d36 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d46 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d56 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d66 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d76 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        b0 = b[7];
        b1 = b[15];
        b2 = b[23];
        b3 = b[31];
        b4 = b[39];
        b5 = b[47];
        b6 = b[55];
        b7 = b[63];
        d07 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
              a[5] * b5 + a[6] * b6 + a[7] * b7;
        d17 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
              a[13] * b5 + a[14] * b6 + a[15] * b7;
        d27 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
              a[21] * b5 + a[22] * b6 + a[23] * b7;
        d37 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
              a[29] * b5 + a[30] * b6 + a[31] * b7;
        d47 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
              a[37] * b5 + a[38] * b6 + a[39] * b7;
        d57 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
              a[45] * b5 + a[46] * b6 + a[47] * b7;
        d67 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
              a[53] * b5 + a[54] * b6 + a[55] * b7;
        d77 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
              a[61] * b5 + a[62] * b6 + a[63] * b7;

        int k = jp + 1;
        int k_end = frowp[row + 1];
        a = &F[64 * k];

        int p = diag[j] + 1;
        int p_end = rowp[j + 1];
        b = &A[64 * p];

        // Now, scan through row j starting at the first entry past the diagonal
        for (; (p < p_end) && (k < k_end); p++) {
          // Determine where the two rows have the same elements
          while (k < k_end && fcols[k] < cols[p]) {
            k++;
            a += 64;
          }

          if (k < k_end && fcols[k] == cols[p]) {
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40], b6 = b[48], b7 = b[56];
            a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                    d15 * b5 + d16 * b6 + d17 * b7;
            a[16] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[24] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[32] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[40] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[48] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[56] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41], b6 = b[49], b7 = b[57];
            a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                    d15 * b5 + d16 * b6 + d17 * b7;
            a[17] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[25] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[33] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[41] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[49] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[57] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42], b6 = b[50], b7 = b[58];
            a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[10] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[18] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[26] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[34] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[42] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[50] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[58] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43], b6 = b[51], b7 = b[59];
            a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[11] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[19] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[27] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[35] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[43] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[51] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[59] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44], b6 = b[52], b7 = b[60];
            a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[12] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[20] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[28] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[36] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[44] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[52] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[60] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45], b6 = b[53], b7 = b[61];
            a[5] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[13] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[21] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[29] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[37] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[45] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[53] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[61] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46], b6 = b[54], b7 = b[62];
            a[6] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[14] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[22] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[30] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[38] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[46] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[54] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[62] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;

            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47], b6 = b[55], b7 = b[63];
            a[7] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                    d05 * b5 + d06 * b6 + d07 * b7;
            a[15] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                     d15 * b5 + d16 * b6 + d17 * b7;
            a[23] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                     d25 * b5 + d26 * b6 + d27 * b7;
            a[31] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                     d35 * b5 + d36 * b6 + d37 * b7;
            a[39] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                     d45 * b5 + d46 * b6 + d47 * b7;
            a[47] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                     d55 * b5 + d56 * b6 + d57 * b7;
            a[55] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                     d65 * b5 + d66 * b6 + d67 * b7;
            a[63] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                     d75 * b5 + d76 * b6 + d77 * b7;
          }
          b += 64;
        }

        // Copy over the matrix
        a = &F[64 * jp];
        a[0] = d00;
        a[1] = d01;
        a[2] = d02;
        a[3] = d03;
        a[4] = d04;
        a[5] = d05;
        a[6] = d06;
        a[7] = d07;
        a[8] = d10;
        a[9] = d11;
        a[10] = d12;
        a[11] = d13;
        a[12] = d14;
        a[13] = d15;
        a[14] = d16;
        a[15] = d17;
        a[16] = d20;
        a[17] = d21;
        a[18] = d22;
        a[19] = d23;
        a[20] = d24;
        a[21] = d25;
        a[22] = d26;
        a[23] = d27;
        a[24] = d30;
        a[25] = d31;
        a[26] = d32;
        a[27] = d33;
        a[28] = d34;
        a[29] = d35;
        a[30] = d36;
        a[31] = d37;
        a[32] = d40;
        a[33] = d41;
        a[34] = d42;
        a[35] = d43;
        a[36] = d44;
        a[37] = d45;
        a[38] = d46;
        a[39] = d47;
        a[40] = d50;
        a[41] = d51;
        a[42] = d52;
        a[43] = d53;
        a[44] = d54;
        a[45] = d55;
        a[46] = d56;
        a[47] = d57;
        a[48] = d60;
        a[49] = d61;
        a[50] = d62;
        a[51] = d63;
        a[52] = d64;
        a[53] = d65;
        a[54] = d66;
        a[55] = d67;
        a[56] = d70;
        a[57] = d71;
        a[58] = d72;
        a[59] = d73;
        a[60] = d74;
        a[61] = d75;
        a[62] = d76;
        a[63] = d77;
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
void BCSRMatFactor8(BCSRMatData *data) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  TacsScalar d00, d01, d02, d03, d04, d05, d06, d07;
  TacsScalar d10, d11, d12, d13, d14, d15, d16, d17;
  TacsScalar d20, d21, d22, d23, d24, d25, d26, d27;
  TacsScalar d30, d31, d32, d33, d34, d35, d36, d37;
  TacsScalar d40, d41, d42, d43, d44, d45, d46, d47;
  TacsScalar d50, d51, d52, d53, d54, d55, d56, d57;
  TacsScalar d60, d61, d62, d63, d64, d65, d66, d67;
  TacsScalar d70, d71, d72, d73, d74, d75, d76, d77;

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
      TacsScalar *a = &(data->A[64 * j]);
      TacsScalar *b = &(data->A[64 * diag[cj]]);

      // Multiply d = A[j] *A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

      b0 = b[0];
      b1 = b[8];
      b2 = b[16];
      b3 = b[24];
      b4 = b[32];
      b5 = b[40];
      b6 = b[48];
      b7 = b[56];
      d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d10 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d20 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d30 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d40 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d50 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d60 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d70 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[1];
      b1 = b[9];
      b2 = b[17];
      b3 = b[25];
      b4 = b[33];
      b5 = b[41];
      b6 = b[49];
      b7 = b[57];
      d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d11 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d21 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d31 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d41 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d51 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d61 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d71 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[2];
      b1 = b[10];
      b2 = b[18];
      b3 = b[26];
      b4 = b[34];
      b5 = b[42];
      b6 = b[50];
      b7 = b[58];
      d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d12 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d22 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d32 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d42 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d52 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d62 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d72 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[3];
      b1 = b[11];
      b2 = b[19];
      b3 = b[27];
      b4 = b[35];
      b5 = b[43];
      b6 = b[51];
      b7 = b[59];
      d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d13 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d23 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d33 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d43 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d53 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d63 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d73 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[4];
      b1 = b[12];
      b2 = b[20];
      b3 = b[28];
      b4 = b[36];
      b5 = b[44];
      b6 = b[52];
      b7 = b[60];
      d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d14 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d24 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d34 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d44 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d54 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d64 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d74 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[5];
      b1 = b[13];
      b2 = b[21];
      b3 = b[29];
      b4 = b[37];
      b5 = b[45];
      b6 = b[53];
      b7 = b[61];
      d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d15 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d25 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d35 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d45 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d55 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d65 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d75 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[6];
      b1 = b[14];
      b2 = b[22];
      b3 = b[30];
      b4 = b[38];
      b5 = b[46];
      b6 = b[54];
      b7 = b[62];
      d06 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d16 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d26 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d36 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d46 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d56 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d66 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d76 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[7];
      b1 = b[15];
      b2 = b[23];
      b3 = b[31];
      b4 = b[39];
      b5 = b[47];
      b6 = b[55];
      b7 = b[63];
      d07 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d17 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d27 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d37 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d47 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d57 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d67 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d77 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      // Scan through the remainder of the row
      int k = j + 1;
      int p = diag[cj] + 1;
      a = &(data->A[64 * k]);
      b = &(data->A[64 * p]);

      // The final entry for row: cols[j]
      int pend = rowp[cj + 1];

      // Keep track of the number of block matrix products
      int nz = 0;

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < pend) && (k < kend); p++) {
        // Determine where the two rows have the same elements
        while (k < kend && cols[k] < cols[p]) {
          k++;
          a += 64;
        }

        // A[k] = A[k] - A[j] *A[p]
        if (k < kend && cols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[8];
          b2 = b[16];
          b3 = b[24];
          b4 = b[32];
          b5 = b[40], b6 = b[48], b7 = b[56];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                  d15 * b5 + d16 * b6 + d17 * b7;
          a[16] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[24] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[32] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[40] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[48] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[56] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[1];
          b1 = b[9];
          b2 = b[17];
          b3 = b[25];
          b4 = b[33];
          b5 = b[41], b6 = b[49], b7 = b[57];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                  d15 * b5 + d16 * b6 + d17 * b7;
          a[17] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[25] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[33] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[41] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[49] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[57] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[2];
          b1 = b[10];
          b2 = b[18];
          b3 = b[26];
          b4 = b[34];
          b5 = b[42], b6 = b[50], b7 = b[58];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[10] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[18] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[26] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[34] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[42] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[50] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[58] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[3];
          b1 = b[11];
          b2 = b[19];
          b3 = b[27];
          b4 = b[35];
          b5 = b[43], b6 = b[51], b7 = b[59];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[11] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[19] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[27] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[35] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[43] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[51] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[59] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[4];
          b1 = b[12];
          b2 = b[20];
          b3 = b[28];
          b4 = b[36];
          b5 = b[44], b6 = b[52], b7 = b[60];
          a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[12] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[20] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[28] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[36] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[44] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[52] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[60] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[5];
          b1 = b[13];
          b2 = b[21];
          b3 = b[29];
          b4 = b[37];
          b5 = b[45], b6 = b[53], b7 = b[61];
          a[5] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[13] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[21] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[29] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[37] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[45] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[53] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[61] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[6];
          b1 = b[14];
          b2 = b[22];
          b3 = b[30];
          b4 = b[38];
          b5 = b[46], b6 = b[54], b7 = b[62];
          a[6] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[14] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[22] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[30] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[38] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[46] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[54] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[62] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[7];
          b1 = b[15];
          b2 = b[23];
          b3 = b[31];
          b4 = b[39];
          b5 = b[47], b6 = b[55], b7 = b[63];
          a[7] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[15] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[23] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[31] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[39] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[47] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[55] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[63] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          nz++;
        }

        b += 64;
      }

      TacsAddFlops(2 * 64 * 8 * nz + 15 * 64);

      // Copy the matrix back into the row
      a = &(data->A[64 * j]);
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d05;
      a[6] = d06;
      a[7] = d07;
      a[8] = d10;
      a[9] = d11;
      a[10] = d12;
      a[11] = d13;
      a[12] = d14;
      a[13] = d15;
      a[14] = d16;
      a[15] = d17;
      a[16] = d20;
      a[17] = d21;
      a[18] = d22;
      a[19] = d23;
      a[20] = d24;
      a[21] = d25;
      a[22] = d26;
      a[23] = d27;
      a[24] = d30;
      a[25] = d31;
      a[26] = d32;
      a[27] = d33;
      a[28] = d34;
      a[29] = d35;
      a[30] = d36;
      a[31] = d37;
      a[32] = d40;
      a[33] = d41;
      a[34] = d42;
      a[35] = d43;
      a[36] = d44;
      a[37] = d45;
      a[38] = d46;
      a[39] = d47;
      a[40] = d50;
      a[41] = d51;
      a[42] = d52;
      a[43] = d53;
      a[44] = d54;
      a[45] = d55;
      a[46] = d56;
      a[47] = d57;
      a[48] = d60;
      a[49] = d61;
      a[50] = d62;
      a[51] = d63;
      a[52] = d64;
      a[53] = d65;
      a[54] = d66;
      a[55] = d67;
      a[56] = d70;
      a[57] = d71;
      a[58] = d72;
      a[59] = d73;
      a[60] = d74;
      a[61] = d75;
      a[62] = d76;
      a[63] = d77;
    }

    // Invert the diagonal portion of the matrix
    TacsScalar D[64];
    TacsScalar *a = &(data->A[64 * diag[i]]);

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
    D[36] = a[36];
    D[37] = a[37];
    D[38] = a[38];
    D[39] = a[39];
    D[40] = a[40];
    D[41] = a[41];
    D[42] = a[42];
    D[43] = a[43];
    D[44] = a[44];
    D[45] = a[45];
    D[46] = a[46];
    D[47] = a[47];
    D[48] = a[48];
    D[49] = a[49];
    D[50] = a[50];
    D[51] = a[51];
    D[52] = a[52];
    D[53] = a[53];
    D[54] = a[54];
    D[55] = a[55];
    D[56] = a[56];
    D[57] = a[57];
    D[58] = a[58];
    D[59] = a[59];
    D[60] = a[60];
    D[61] = a[61];
    D[62] = a[62];
    D[63] = a[63];

    int ipiv[8];
    int info = BMatComputeInverse(a, D, ipiv, 8);

    if (info > 0) {
      fprintf(stderr,
              "Error during factorization of diagonal %d in block row %d \n",
              i + 1, info);
    }
  }

  // Add flops from the diagonal inversion
  TacsAddFlops(1.333333 * 8 * 8 * 8 * nrows);
}

/*!
  Compute x = L_{B}^{-1} E
*/
void BCSRMatFactorLower8(BCSRMatData *data, BCSRMatData *Edata) {
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
      TacsScalar *d = &(data->A[64 * j]);

      int k = erowp[i];
      int k_end = erowp[i + 1];
      TacsScalar *a = &(Edata->A[64 * k]);

      int p = erowp[cj];
      int p_end = erowp[cj + 1];
      TacsScalar *b = &(Edata->A[64 * p]);

      // Now, scan through row cj starting at the first entry past the
      // diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && ecols[k] < ecols[p]) {
          k++;
          a += 64;
        }

        if (k < k_end && ecols[k] == ecols[p]) {
          TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;
          b0 = b[0];
          b1 = b[8];
          b2 = b[16];
          b3 = b[24];
          b4 = b[32];
          b5 = b[40];
          b6 = b[48];
          b7 = b[56];
          a[0] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[8] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 + d[12] * b4 +
                  d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[16] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[24] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[32] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[40] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[48] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[56] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[1];
          b1 = b[9];
          b2 = b[17];
          b3 = b[25];
          b4 = b[33];
          b5 = b[41];
          b6 = b[49];
          b7 = b[57];
          a[1] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[9] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 + d[12] * b4 +
                  d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[17] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[25] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[33] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[41] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[49] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[57] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[2];
          b1 = b[10];
          b2 = b[18];
          b3 = b[26];
          b4 = b[34];
          b5 = b[42];
          b6 = b[50];
          b7 = b[58];
          a[2] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[10] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[18] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[26] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[34] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[42] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[50] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[58] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[3];
          b1 = b[11];
          b2 = b[19];
          b3 = b[27];
          b4 = b[35];
          b5 = b[43];
          b6 = b[51];
          b7 = b[59];
          a[3] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[11] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[19] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[27] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[35] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[43] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[51] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[59] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[4];
          b1 = b[12];
          b2 = b[20];
          b3 = b[28];
          b4 = b[36];
          b5 = b[44];
          b6 = b[52];
          b7 = b[60];
          a[4] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[12] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[20] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[28] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[36] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[44] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[52] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[60] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[5];
          b1 = b[13];
          b2 = b[21];
          b3 = b[29];
          b4 = b[37];
          b5 = b[45];
          b6 = b[53];
          b7 = b[61];
          a[5] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[13] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[21] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[29] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[37] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[45] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[53] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[61] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[6];
          b1 = b[14];
          b2 = b[22];
          b3 = b[30];
          b4 = b[38];
          b5 = b[46];
          b6 = b[54];
          b7 = b[62];
          a[6] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[14] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[22] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[30] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[38] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[46] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[54] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[62] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          b0 = b[7];
          b1 = b[15];
          b2 = b[23];
          b3 = b[31];
          b4 = b[39];
          b5 = b[47];
          b6 = b[55];
          b7 = b[63];
          a[7] -= d[0] * b0 + d[1] * b1 + d[2] * b2 + d[3] * b3 + d[4] * b4 +
                  d[5] * b5 + d[6] * b6 + d[7] * b7;
          a[15] -= d[8] * b0 + d[9] * b1 + d[10] * b2 + d[11] * b3 +
                   d[12] * b4 + d[13] * b5 + d[14] * b6 + d[15] * b7;
          a[23] -= d[16] * b0 + d[17] * b1 + d[18] * b2 + d[19] * b3 +
                   d[20] * b4 + d[21] * b5 + d[22] * b6 + d[23] * b7;
          a[31] -= d[24] * b0 + d[25] * b1 + d[26] * b2 + d[27] * b3 +
                   d[28] * b4 + d[29] * b5 + d[30] * b6 + d[31] * b7;
          a[39] -= d[32] * b0 + d[33] * b1 + d[34] * b2 + d[35] * b3 +
                   d[36] * b4 + d[37] * b5 + d[38] * b6 + d[39] * b7;
          a[47] -= d[40] * b0 + d[41] * b1 + d[42] * b2 + d[43] * b3 +
                   d[44] * b4 + d[45] * b5 + d[46] * b6 + d[47] * b7;
          a[55] -= d[48] * b0 + d[49] * b1 + d[50] * b2 + d[51] * b3 +
                   d[52] * b4 + d[53] * b5 + d[54] * b6 + d[55] * b7;
          a[63] -= d[56] * b0 + d[57] * b1 + d[58] * b2 + d[59] * b3 +
                   d[60] * b4 + d[61] * b5 + d[62] * b6 + d[63] * b7;

          nz++;
        }
        b += 64;
      }
    }
  }

  TacsAddFlops(2 * 64 * 8 * nz);
}

/*!
  Compute x = F U_{B}^{-1}
*/
void BCSRMatFactorUpper8(BCSRMatData *data, BCSRMatData *Fdata) {
  // Retrieve the data required from the matrix
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  // Retrieve the data required from the matrix
  const int nrows_f = Fdata->nrows;
  const int *frowp = Fdata->rowp;
  const int *fcols = Fdata->cols;

  TacsScalar d00, d01, d02, d03, d04, d05, d06, d07;
  TacsScalar d10, d11, d12, d13, d14, d15, d16, d17;
  TacsScalar d20, d21, d22, d23, d24, d25, d26, d27;
  TacsScalar d30, d31, d32, d33, d34, d35, d36, d37;
  TacsScalar d40, d41, d42, d43, d44, d45, d46, d47;
  TacsScalar d50, d51, d52, d53, d54, d55, d56, d57;
  TacsScalar d60, d61, d62, d63, d64, d65, d66, d67;
  TacsScalar d70, d71, d72, d73, d74, d75, d76, d77;

  for (int i = 0; i < nrows_f; i++) {
    int j_end = frowp[i + 1];

    for (int j = frowp[i]; j < j_end; j++) {
      int cj = fcols[j];
      TacsScalar *a = &(Fdata->A[64 * j]);
      const TacsScalar *b = &(data->A[64 * diag[cj]]);

      // Multiply d = F[j]*A[diag[cj]]
      TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

      b0 = b[0];
      b1 = b[8];
      b2 = b[16];
      b3 = b[24];
      b4 = b[32];
      b5 = b[40];
      b6 = b[48];
      b7 = b[56];
      d00 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d10 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d20 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d30 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d40 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d50 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d60 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d70 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[1];
      b1 = b[9];
      b2 = b[17];
      b3 = b[25];
      b4 = b[33];
      b5 = b[41];
      b6 = b[49];
      b7 = b[57];
      d01 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d11 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d21 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d31 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d41 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d51 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d61 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d71 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[2];
      b1 = b[10];
      b2 = b[18];
      b3 = b[26];
      b4 = b[34];
      b5 = b[42];
      b6 = b[50];
      b7 = b[58];
      d02 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d12 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d22 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d32 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d42 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d52 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d62 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d72 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[3];
      b1 = b[11];
      b2 = b[19];
      b3 = b[27];
      b4 = b[35];
      b5 = b[43];
      b6 = b[51];
      b7 = b[59];
      d03 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d13 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d23 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d33 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d43 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d53 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d63 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d73 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[4];
      b1 = b[12];
      b2 = b[20];
      b3 = b[28];
      b4 = b[36];
      b5 = b[44];
      b6 = b[52];
      b7 = b[60];
      d04 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d14 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d24 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d34 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d44 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d54 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d64 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d74 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[5];
      b1 = b[13];
      b2 = b[21];
      b3 = b[29];
      b4 = b[37];
      b5 = b[45];
      b6 = b[53];
      b7 = b[61];
      d05 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d15 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d25 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d35 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d45 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d55 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d65 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d75 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[6];
      b1 = b[14];
      b2 = b[22];
      b3 = b[30];
      b4 = b[38];
      b5 = b[46];
      b6 = b[54];
      b7 = b[62];
      d06 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d16 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d26 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d36 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d46 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d56 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d66 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d76 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      b0 = b[7];
      b1 = b[15];
      b2 = b[23];
      b3 = b[31];
      b4 = b[39];
      b5 = b[47];
      b6 = b[55];
      b7 = b[63];
      d07 = a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
            a[5] * b5 + a[6] * b6 + a[7] * b7;
      d17 = a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 + a[12] * b4 +
            a[13] * b5 + a[14] * b6 + a[15] * b7;
      d27 = a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 + a[20] * b4 +
            a[21] * b5 + a[22] * b6 + a[23] * b7;
      d37 = a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 + a[28] * b4 +
            a[29] * b5 + a[30] * b6 + a[31] * b7;
      d47 = a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 + a[36] * b4 +
            a[37] * b5 + a[38] * b6 + a[39] * b7;
      d57 = a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 + a[44] * b4 +
            a[45] * b5 + a[46] * b6 + a[47] * b7;
      d67 = a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 + a[52] * b4 +
            a[53] * b5 + a[54] * b6 + a[55] * b7;
      d77 = a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 + a[60] * b4 +
            a[61] * b5 + a[62] * b6 + a[63] * b7;

      int k = j + 1;
      int k_end = frowp[i + 1];
      a = &(Fdata->A[64 * k]);

      int p = diag[cj] + 1;
      int p_end = rowp[cj + 1];
      b = &(data->A[64 * p]);

      // Keep track of the number of block matrix-matrix products
      int nz = 0;

      // Now, scan through row cj starting at the first entry past the diagonal
      for (; (p < p_end) && (k < k_end); p++) {
        // Determine where the two rows have the same elements
        while (k < k_end && fcols[k] < cols[p]) {
          k++;
          a += 64;
        }

        if (k < k_end && fcols[k] == cols[p]) {
          b0 = b[0];
          b1 = b[8];
          b2 = b[16];
          b3 = b[24];
          b4 = b[32];
          b5 = b[40], b6 = b[48], b7 = b[56];
          a[0] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[8] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                  d15 * b5 + d16 * b6 + d17 * b7;
          a[16] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[24] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[32] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[40] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[48] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[56] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[1];
          b1 = b[9];
          b2 = b[17];
          b3 = b[25];
          b4 = b[33];
          b5 = b[41], b6 = b[49], b7 = b[57];
          a[1] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[9] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                  d15 * b5 + d16 * b6 + d17 * b7;
          a[17] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[25] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[33] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[41] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[49] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[57] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[2];
          b1 = b[10];
          b2 = b[18];
          b3 = b[26];
          b4 = b[34];
          b5 = b[42], b6 = b[50], b7 = b[58];
          a[2] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[10] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[18] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[26] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[34] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[42] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[50] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[58] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[3];
          b1 = b[11];
          b2 = b[19];
          b3 = b[27];
          b4 = b[35];
          b5 = b[43], b6 = b[51], b7 = b[59];
          a[3] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[11] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[19] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[27] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[35] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[43] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[51] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[59] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[4];
          b1 = b[12];
          b2 = b[20];
          b3 = b[28];
          b4 = b[36];
          b5 = b[44], b6 = b[52], b7 = b[60];
          a[4] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[12] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[20] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[28] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[36] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[44] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[52] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[60] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[5];
          b1 = b[13];
          b2 = b[21];
          b3 = b[29];
          b4 = b[37];
          b5 = b[45], b6 = b[53], b7 = b[61];
          a[5] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[13] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[21] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[29] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[37] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[45] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[53] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[61] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[6];
          b1 = b[14];
          b2 = b[22];
          b3 = b[30];
          b4 = b[38];
          b5 = b[46], b6 = b[54], b7 = b[62];
          a[6] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[14] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[22] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[30] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[38] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[46] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[54] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[62] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          b0 = b[7];
          b1 = b[15];
          b2 = b[23];
          b3 = b[31];
          b4 = b[39];
          b5 = b[47], b6 = b[55], b7 = b[63];
          a[7] -= d00 * b0 + d01 * b1 + d02 * b2 + d03 * b3 + d04 * b4 +
                  d05 * b5 + d06 * b6 + d07 * b7;
          a[15] -= d10 * b0 + d11 * b1 + d12 * b2 + d13 * b3 + d14 * b4 +
                   d15 * b5 + d16 * b6 + d17 * b7;
          a[23] -= d20 * b0 + d21 * b1 + d22 * b2 + d23 * b3 + d24 * b4 +
                   d25 * b5 + d26 * b6 + d27 * b7;
          a[31] -= d30 * b0 + d31 * b1 + d32 * b2 + d33 * b3 + d34 * b4 +
                   d35 * b5 + d36 * b6 + d37 * b7;
          a[39] -= d40 * b0 + d41 * b1 + d42 * b2 + d43 * b3 + d44 * b4 +
                   d45 * b5 + d46 * b6 + d47 * b7;
          a[47] -= d50 * b0 + d51 * b1 + d52 * b2 + d53 * b3 + d54 * b4 +
                   d55 * b5 + d56 * b6 + d57 * b7;
          a[55] -= d60 * b0 + d61 * b1 + d62 * b2 + d63 * b3 + d64 * b4 +
                   d65 * b5 + d66 * b6 + d67 * b7;
          a[63] -= d70 * b0 + d71 * b1 + d72 * b2 + d73 * b3 + d74 * b4 +
                   d75 * b5 + d76 * b6 + d77 * b7;

          nz++;
        }
        b += 64;
      }

      TacsAddFlops(2 * 64 * 8 * nz + 15 * 64);

      // Copy over the matrix
      a = &(Fdata->A[64 * j]);
      a[0] = d00;
      a[1] = d01;
      a[2] = d02;
      a[3] = d03;
      a[4] = d04;
      a[5] = d05;
      a[6] = d06;
      a[7] = d07;
      a[8] = d10;
      a[9] = d11;
      a[10] = d12;
      a[11] = d13;
      a[12] = d14;
      a[13] = d15;
      a[14] = d16;
      a[15] = d17;
      a[16] = d20;
      a[17] = d21;
      a[18] = d22;
      a[19] = d23;
      a[20] = d24;
      a[21] = d25;
      a[22] = d26;
      a[23] = d27;
      a[24] = d30;
      a[25] = d31;
      a[26] = d32;
      a[27] = d33;
      a[28] = d34;
      a[29] = d35;
      a[30] = d36;
      a[31] = d37;
      a[32] = d40;
      a[33] = d41;
      a[34] = d42;
      a[35] = d43;
      a[36] = d44;
      a[37] = d45;
      a[38] = d46;
      a[39] = d47;
      a[40] = d50;
      a[41] = d51;
      a[42] = d52;
      a[43] = d53;
      a[44] = d54;
      a[45] = d55;
      a[46] = d56;
      a[47] = d57;
      a[48] = d60;
      a[49] = d61;
      a[50] = d62;
      a[51] = d63;
      a[52] = d64;
      a[53] = d65;
      a[54] = d66;
      a[55] = d67;
      a[56] = d70;
      a[57] = d71;
      a[58] = d72;
      a[59] = d73;
      a[60] = d74;
      a[61] = d75;
      a[62] = d76;
      a[63] = d77;
    }
  }
}
