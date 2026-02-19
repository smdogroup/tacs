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
  size = 8.
*/

/*
  Threaded implementation of matrix-multiplication
*/
void *BCSRMatVecMultAdd8_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);
  const int nrows = tdata->mat->nrows;

  // Get the input/output vectors
  const TacsScalar *x = tdata->input;

  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const TacsScalar *A = tdata->mat->A;
  const int group_size = tdata->mat->matvec_group_size;

  while (tdata->num_completed_rows < nrows) {
    int row = -1;
    tdata->mat_mult_sched_job(group_size, &row);

    if (row >= 0) {
      TacsScalar *y = &tdata->output[8 * row];
      int k = rowp[row];
      for (int ii = row; ii < nrows && (ii < row + group_size); ii++) {
        int end = rowp[ii + 1];
        const TacsScalar *a = &A[64 * k];

        for (; k < end; k++) {
          int j = 8 * cols[k];

          y[0] += a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
                  a[3] * x[j + 3] + a[4] * x[j + 4] + a[5] * x[j + 5] +
                  a[6] * x[j + 6] + a[7] * x[j + 7];
          y[1] += a[8] * x[j] + a[9] * x[j + 1] + a[10] * x[j + 2] +
                  a[11] * x[j + 3] + a[12] * x[j + 4] + a[13] * x[j + 5] +
                  a[14] * x[j + 6] + a[15] * x[j + 7];
          y[2] += a[16] * x[j] + a[17] * x[j + 1] + a[18] * x[j + 2] +
                  a[19] * x[j + 3] + a[20] * x[j + 4] + a[21] * x[j + 5] +
                  a[22] * x[j + 6] + a[23] * x[j + 7];
          y[3] += a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
                  a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5] +
                  a[30] * x[j + 6] + a[31] * x[j + 7];
          y[4] += a[32] * x[j] + a[33] * x[j + 1] + a[34] * x[j + 2] +
                  a[35] * x[j + 3] + a[36] * x[j + 4] + a[37] * x[j + 5] +
                  a[38] * x[j + 6] + a[39] * x[j + 7];
          y[5] += a[40] * x[j] + a[41] * x[j + 1] + a[42] * x[j + 2] +
                  a[43] * x[j + 3] + a[44] * x[j + 4] + a[45] * x[j + 5] +
                  a[46] * x[j + 6] + a[47] * x[j + 7];
          y[6] += a[48] * x[j] + a[49] * x[j + 1] + a[50] * x[j + 2] +
                  a[51] * x[j + 3] + a[52] * x[j + 4] + a[53] * x[j + 5] +
                  a[54] * x[j + 6] + a[55] * x[j + 7];
          y[7] += a[56] * x[j] + a[57] * x[j + 1] + a[58] * x[j + 2] +
                  a[59] * x[j + 3] + a[60] * x[j + 4] + a[61] * x[j + 5] +
                  a[62] * x[j + 6] + a[63] * x[j + 7];

          a += 64;
        }

        y += 8;
      }
    }
  }

  pthread_exit(NULL);
}

/*!
  Compute the matrix-vector product: y = A * x
*/
void BCSRMatVecMult8(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  TacsScalar *a = data->A;

  for (int i = 0; i < nrows; i++) {
    int end = rowp[i + 1];

    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = 0.0;
    y[4] = 0.0;
    y[5] = 0.0;
    y[6] = 0.0;
    y[7] = 0.0;

    for (int k = rowp[i]; k < end; k++) {
      int j = 8 * cols[k];

      y[0] += a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
              a[3] * x[j + 3] + a[4] * x[j + 4] + a[5] * x[j + 5] +
              a[6] * x[j + 6] + a[7] * x[j + 7];
      y[1] += a[8] * x[j] + a[9] * x[j + 1] + a[10] * x[j + 2] +
              a[11] * x[j + 3] + a[12] * x[j + 4] + a[13] * x[j + 5] +
              a[14] * x[j + 6] + a[15] * x[j + 7];
      y[2] += a[16] * x[j] + a[17] * x[j + 1] + a[18] * x[j + 2] +
              a[19] * x[j + 3] + a[20] * x[j + 4] + a[21] * x[j + 5] +
              a[22] * x[j + 6] + a[23] * x[j + 7];
      y[3] += a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
              a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5] +
              a[30] * x[j + 6] + a[31] * x[j + 7];
      y[4] += a[32] * x[j] + a[33] * x[j + 1] + a[34] * x[j + 2] +
              a[35] * x[j + 3] + a[36] * x[j + 4] + a[37] * x[j + 5] +
              a[38] * x[j + 6] + a[39] * x[j + 7];
      y[5] += a[40] * x[j] + a[41] * x[j + 1] + a[42] * x[j + 2] +
              a[43] * x[j + 3] + a[44] * x[j + 4] + a[45] * x[j + 5] +
              a[46] * x[j + 6] + a[47] * x[j + 7];
      y[6] += a[48] * x[j] + a[49] * x[j + 1] + a[50] * x[j + 2] +
              a[51] * x[j + 3] + a[52] * x[j + 4] + a[53] * x[j + 5] +
              a[54] * x[j + 6] + a[55] * x[j + 7];
      y[7] += a[56] * x[j] + a[57] * x[j + 1] + a[58] * x[j + 2] +
              a[59] * x[j + 3] + a[60] * x[j + 4] + a[61] * x[j + 5] +
              a[62] * x[j + 6] + a[63] * x[j + 7];

      a += 64;
    }

    y += 8;
  }

  TacsAddFlops(2 * 64 * rowp[nrows]);
}

/*!
  Compute the matrix vector product plus addition: z = A * x + y
*/

void BCSRMatVecMultAdd8(BCSRMatData *data, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z) {
  // Retrieve the data required from the matrix
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  TacsScalar *a = data->A;

  for (int i = 0; i < nrows; i++) {
    int end = rowp[i + 1];

    z[0] = y[0];
    z[1] = y[1];
    z[2] = y[2];
    z[3] = y[3];
    z[4] = y[4];
    z[5] = y[5];
    z[6] = y[6];
    z[7] = y[7];

    for (int k = rowp[i]; k < end; k++) {
      int j = 8 * cols[k];

      z[0] += a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
              a[3] * x[j + 3] + a[4] * x[j + 4] + a[5] * x[j + 5] +
              a[6] * x[j + 6] + a[7] * x[j + 7];
      z[1] += a[8] * x[j] + a[9] * x[j + 1] + a[10] * x[j + 2] +
              a[11] * x[j + 3] + a[12] * x[j + 4] + a[13] * x[j + 5] +
              a[14] * x[j + 6] + a[15] * x[j + 7];
      z[2] += a[16] * x[j] + a[17] * x[j + 1] + a[18] * x[j + 2] +
              a[19] * x[j + 3] + a[20] * x[j + 4] + a[21] * x[j + 5] +
              a[22] * x[j + 6] + a[23] * x[j + 7];
      z[3] += a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
              a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5] +
              a[30] * x[j + 6] + a[31] * x[j + 7];
      z[4] += a[32] * x[j] + a[33] * x[j + 1] + a[34] * x[j + 2] +
              a[35] * x[j + 3] + a[36] * x[j + 4] + a[37] * x[j + 5] +
              a[38] * x[j + 6] + a[39] * x[j + 7];
      z[5] += a[40] * x[j] + a[41] * x[j + 1] + a[42] * x[j + 2] +
              a[43] * x[j + 3] + a[44] * x[j + 4] + a[45] * x[j + 5] +
              a[46] * x[j + 6] + a[47] * x[j + 7];
      z[6] += a[48] * x[j] + a[49] * x[j + 1] + a[50] * x[j + 2] +
              a[51] * x[j + 3] + a[52] * x[j + 4] + a[53] * x[j + 5] +
              a[54] * x[j + 6] + a[55] * x[j + 7];
      z[7] += a[56] * x[j] + a[57] * x[j + 1] + a[58] * x[j + 2] +
              a[59] * x[j + 3] + a[60] * x[j + 4] + a[61] * x[j + 5] +
              a[62] * x[j + 6] + a[63] * x[j + 7];
      a += 64;
    }

    y += 8;
    z += 8;
  }

  TacsAddFlops(2 * 64 * rowp[nrows]);
}

/*
  Apply the lower-triangular matrix over a column range of a row
  obtained from a scheduling function.

  Compute:
  y <- L[irow, jstart:jend]^{-1} y

  where y is the output array.

  Note that the scheduler ensures that no two threads are operating on
  the same data simultaneously.
*/
void *BCSRMatApplyLower8_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);
  const int nrows = tdata->mat->nrows;
  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  const TacsScalar *A = tdata->mat->A;
  const int group_size = tdata->mat->matvec_group_size;

  TacsScalar *y = tdata->output;

  while (tdata->num_completed_rows < nrows) {
    int index, irow, jstart, jend;
    tdata->apply_lower_sched_job(group_size, &index, &irow, &jstart, &jend);

    if (irow >= 0) {
      TacsScalar *z = &y[8 * irow];

      for (int i = irow; (i < nrows) && (i < irow + group_size); i++) {
        int end = diag[i];
        int k = rowp[i];
        while ((k < end) && (cols[k] < jstart)) {
          k++;
        }

        const TacsScalar *a = &A[64 * k];
        for (; (k < end) && (cols[k] < jend); k++) {
          int j = 6 * cols[k];

          z[0] -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                  a[3] * y[j + 3] + a[4] * y[j + 4] + a[5] * y[j + 5] +
                  a[6] * y[j + 6] + a[7] * y[j + 7];
          z[1] -= a[8] * y[j] + a[9] * y[j + 1] + a[10] * y[j + 2] +
                  a[11] * y[j + 3] + a[12] * y[j + 4] + a[13] * y[j + 5] +
                  a[14] * y[j + 6] + a[15] * y[j + 7];
          z[2] -= a[16] * y[j] + a[17] * y[j + 1] + a[18] * y[j + 2] +
                  a[19] * y[j + 3] + a[20] * y[j + 4] + a[21] * y[j + 5] +
                  a[22] * y[j + 6] + a[23] * y[j + 7];
          z[3] -= a[24] * y[j] + a[25] * y[j + 1] + a[26] * y[j + 2] +
                  a[27] * y[j + 3] + a[28] * y[j + 4] + a[29] * y[j + 5] +
                  a[30] * y[j + 6] + a[31] * y[j + 7];
          z[4] -= a[32] * y[j] + a[33] * y[j + 1] + a[34] * y[j + 2] +
                  a[35] * y[j + 3] + a[36] * y[j + 4] + a[37] * y[j + 5] +
                  a[38] * y[j + 6] + a[39] * y[j + 7];
          z[5] -= a[40] * y[j] + a[41] * y[j + 1] + a[42] * y[j + 2] +
                  a[43] * y[j + 3] + a[44] * y[j + 4] + a[45] * y[j + 5] +
                  a[46] * y[j + 6] + a[47] * y[j + 7];
          z[6] -= a[48] * y[j] + a[49] * y[j + 1] + a[50] * y[j + 2] +
                  a[51] * y[j + 3] + a[52] * y[j + 4] + a[53] * y[j + 5] +
                  a[54] * y[j + 6] + a[55] * y[j + 7];
          z[7] -= a[56] * y[j] + a[57] * y[j + 1] + a[58] * y[j + 2] +
                  a[59] * y[j + 3] + a[60] * y[j + 4] + a[61] * y[j + 5] +
                  a[62] * y[j + 6] + a[63] * y[j + 7];
          a += 64;
        }

        z += 8;
      }

      tdata->apply_lower_mark_completed(group_size, index, irow, jstart, jend);
    }
  }

  pthread_exit(NULL);
}

/*
  Apply the upper-triangular matrix over a range of columns dictated
  by a scheduler.

  Compute:
  y <- U[irow, jstart:jend]^{-1} y
*/
void *BCSRMatApplyUpper8_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);
  const int nrows = tdata->mat->nrows;
  const int *rowp = tdata->mat->rowp;
  const int *cols = tdata->mat->cols;
  const int *diag = tdata->mat->diag;
  const TacsScalar *A = tdata->mat->A;
  const int group_size = tdata->mat->matvec_group_size;

  TacsScalar *y = tdata->output;

  while (tdata->num_completed_rows < nrows) {
    int index, irow, jstart, jend;
    tdata->apply_upper_sched_job(group_size, &index, &irow, &jstart, &jend);

    if (irow >= 0) {
      for (int i = irow - 1; (i >= 0) && (i >= irow - group_size); i--) {
        int start = diag[i] + 1;
        int end = rowp[i + 1];

        int k = end - 1;
        while ((k >= start) && (cols[k] >= jend)) {
          k--;
        }

        TacsScalar *z = &y[8 * i];
        const TacsScalar *a = &A[64 * k];
        for (; (k >= start) && (cols[k] >= jstart); k--) {
          int j = 8 * cols[k];

          z[0] -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
                  a[3] * y[j + 3] + a[4] * y[j + 4] + a[5] * y[j + 5] +
                  a[6] * y[j + 6] + a[7] * y[j + 7];
          z[1] -= a[8] * y[j] + a[9] * y[j + 1] + a[10] * y[j + 2] +
                  a[11] * y[j + 3] + a[12] * y[j + 4] + a[13] * y[j + 5] +
                  a[14] * y[j + 6] + a[15] * y[j + 7];
          z[2] -= a[16] * y[j] + a[17] * y[j + 1] + a[18] * y[j + 2] +
                  a[19] * y[j + 3] + a[20] * y[j + 4] + a[21] * y[j + 5] +
                  a[22] * y[j + 6] + a[23] * y[j + 7];
          z[3] -= a[24] * y[j] + a[25] * y[j + 1] + a[26] * y[j + 2] +
                  a[27] * y[j + 3] + a[28] * y[j + 4] + a[29] * y[j + 5] +
                  a[30] * y[j + 6] + a[31] * y[j + 7];
          z[4] -= a[32] * y[j] + a[33] * y[j + 1] + a[34] * y[j + 2] +
                  a[35] * y[j + 3] + a[36] * y[j + 4] + a[37] * y[j + 5] +
                  a[38] * y[j + 6] + a[39] * y[j + 7];
          z[5] -= a[40] * y[j] + a[41] * y[j + 1] + a[42] * y[j + 2] +
                  a[43] * y[j + 3] + a[44] * y[j + 4] + a[45] * y[j + 5] +
                  a[46] * y[j + 6] + a[47] * y[j + 7];
          z[6] -= a[48] * y[j] + a[49] * y[j + 1] + a[50] * y[j + 2] +
                  a[51] * y[j + 3] + a[52] * y[j + 4] + a[53] * y[j + 5] +
                  a[54] * y[j + 6] + a[55] * y[j + 7];
          z[7] -= a[56] * y[j] + a[57] * y[j + 1] + a[58] * y[j + 2] +
                  a[59] * y[j + 3] + a[60] * y[j + 4] + a[61] * y[j + 5] +
                  a[62] * y[j + 6] + a[63] * y[j + 7];
          a -= 64;
        }

        if (irow == jstart + group_size) {
          TacsScalar y0 = z[0], y1 = z[1], y2 = z[2], y3 = z[3];
          TacsScalar y4 = z[4], y5 = z[5], y6 = z[6], y7 = z[7];

          a = &A[64 * (start - 1)];
          z[0] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4 +
                 a[5] * y5 + a[6] * y6 + a[7] * y7;
          z[1] = a[8] * y0 + a[9] * y1 + a[10] * y2 + a[11] * y3 + a[12] * y4 +
                 a[13] * y5 + a[14] * y6 + a[15] * y7;
          z[2] = a[16] * y0 + a[17] * y1 + a[18] * y2 + a[19] * y3 +
                 a[20] * y4 + a[21] * y5 + a[22] * y6 + a[23] * y7;
          z[3] = a[24] * y0 + a[25] * y1 + a[26] * y2 + a[27] * y3 +
                 a[28] * y4 + a[29] * y5 + a[30] * y6 + a[31] * y7;
          z[4] = a[32] * y0 + a[33] * y1 + a[34] * y2 + a[35] * y3 +
                 a[36] * y4 + a[37] * y5 + a[38] * y6 + a[39] * y7;
          z[5] = a[40] * y0 + a[41] * y1 + a[42] * y2 + a[43] * y3 +
                 a[44] * y4 + a[45] * y5 + a[46] * y6 + a[47] * y7;
          z[6] = a[48] * y0 + a[49] * y1 + a[50] * y2 + a[51] * y3 +
                 a[52] * y4 + a[53] * y5 + a[54] * y6 + a[55] * y7;
          z[7] = a[56] * y0 + a[57] * y1 + a[58] * y2 + a[59] * y3 +
                 a[60] * y4 + a[61] * y5 + a[62] * y6 + a[63] * y7;
        }
      }

      tdata->apply_upper_mark_completed(group_size, index, irow, jstart, jend);
    }
  }

  pthread_exit(NULL);
}

/*!
  Apply the lower factorization y = L^{-1} x
*/
void BCSRMatApplyLower8(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *diag = data->diag;
  const int *rowp = data->rowp;
  const int *cols = data->cols;

  TacsScalar *z = y;
  for (int i = 0; i < nrows; i++) {
    z[0] = x[0];
    z[1] = x[1];
    z[2] = x[2];
    z[3] = x[3];
    z[4] = x[4];
    z[5] = x[5];
    z[6] = x[6];
    z[7] = x[7];

    int end = diag[i];
    int k = rowp[i];
    const TacsScalar *a = &(data->A[64 * k]);
    for (; k < end; k++) {
      int j = 8 * cols[k];

      z[0] -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] +
              a[3] * y[j + 3] + a[4] * y[j + 4] + a[5] * y[j + 5] +
              a[6] * y[j + 6] + a[7] * y[j + 7];
      z[1] -= a[8] * y[j] + a[9] * y[j + 1] + a[10] * y[j + 2] +
              a[11] * y[j + 3] + a[12] * y[j + 4] + a[13] * y[j + 5] +
              a[14] * y[j + 6] + a[15] * y[j + 7];
      z[2] -= a[16] * y[j] + a[17] * y[j + 1] + a[18] * y[j + 2] +
              a[19] * y[j + 3] + a[20] * y[j + 4] + a[21] * y[j + 5] +
              a[22] * y[j + 6] + a[23] * y[j + 7];
      z[3] -= a[24] * y[j] + a[25] * y[j + 1] + a[26] * y[j + 2] +
              a[27] * y[j + 3] + a[28] * y[j + 4] + a[29] * y[j + 5] +
              a[30] * y[j + 6] + a[31] * y[j + 7];
      z[4] -= a[32] * y[j] + a[33] * y[j + 1] + a[34] * y[j + 2] +
              a[35] * y[j + 3] + a[36] * y[j + 4] + a[37] * y[j + 5] +
              a[38] * y[j + 6] + a[39] * y[j + 7];
      z[5] -= a[40] * y[j] + a[41] * y[j + 1] + a[42] * y[j + 2] +
              a[43] * y[j + 3] + a[44] * y[j + 4] + a[45] * y[j + 5] +
              a[46] * y[j + 6] + a[47] * y[j + 7];
      z[6] -= a[48] * y[j] + a[49] * y[j + 1] + a[50] * y[j + 2] +
              a[51] * y[j + 3] + a[52] * y[j + 4] + a[53] * y[j + 5] +
              a[54] * y[j + 6] + a[55] * y[j + 7];
      z[7] -= a[56] * y[j] + a[57] * y[j + 1] + a[58] * y[j + 2] +
              a[59] * y[j + 3] + a[60] * y[j + 4] + a[61] * y[j + 5] +
              a[62] * y[j + 6] + a[63] * y[j + 7];
      a += 64;
    }

    z += 8;
    x += 8;
    TacsAddFlops(2 * 64 * nz);
  }
}

/*!
  Apply the upper factorization y = U^{-1} x
*/
void BCSRMatApplyUpper8(BCSRMatData *data, TacsScalar *x, TacsScalar *y) {
  const int nrows = data->nrows;
  const int *diag = data->diag;
  const int *rowp = data->rowp;
  const int *cols = data->cols;

  x = &y[8 * (nrows - 1)];
  for (int i = nrows - 1; i >= 0; i--) {
    TacsScalar y0 = x[0], y1 = x[1], y2 = x[2];
    TacsScalar y3 = x[3], y4 = x[4], y5 = x[5];
    TacsScalar y6 = x[6], y7 = x[7];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    TacsScalar *a = &(data->A[64 * k]);
    for (; k < end; k++) {
      int j = 8 * cols[k];

      y0 -= a[0] * y[j] + a[1] * y[j + 1] + a[2] * y[j + 2] + a[3] * y[j + 3] +
            a[4] * y[j + 4] + a[5] * y[j + 5] + a[6] * y[j + 6] +
            a[7] * y[j + 7];
      y1 -= a[8] * y[j] + a[9] * y[j + 1] + a[10] * y[j + 2] +
            a[11] * y[j + 3] + a[12] * y[j + 4] + a[13] * y[j + 5] +
            a[14] * y[j + 6] + a[15] * y[j + 7];
      y2 -= a[16] * y[j] + a[17] * y[j + 1] + a[18] * y[j + 2] +
            a[19] * y[j + 3] + a[20] * y[j + 4] + a[21] * y[j + 5] +
            a[22] * y[j + 6] + a[23] * y[j + 7];
      y3 -= a[24] * y[j] + a[25] * y[j + 1] + a[26] * y[j + 2] +
            a[27] * y[j + 3] + a[28] * y[j + 4] + a[29] * y[j + 5] +
            a[30] * y[j + 6] + a[31] * y[j + 7];
      y4 -= a[32] * y[j] + a[33] * y[j + 1] + a[34] * y[j + 2] +
            a[35] * y[j + 3] + a[36] * y[j + 4] + a[37] * y[j + 5] +
            a[38] * y[j + 6] + a[39] * y[j + 7];
      y5 -= a[40] * y[j] + a[41] * y[j + 1] + a[42] * y[j + 2] +
            a[43] * y[j + 3] + a[44] * y[j + 4] + a[45] * y[j + 5] +
            a[46] * y[j + 6] + a[47] * y[j + 7];
      y6 -= a[48] * y[j] + a[49] * y[j + 1] + a[50] * y[j + 2] +
            a[51] * y[j + 3] + a[52] * y[j + 4] + a[53] * y[j + 5] +
            a[54] * y[j + 6] + a[55] * y[j + 7];
      y7 -= a[56] * y[j] + a[57] * y[j + 1] + a[58] * y[j + 2] +
            a[59] * y[j + 3] + a[60] * y[j + 4] + a[61] * y[j + 5] +
            a[62] * y[j + 6] + a[63] * y[j + 7];

      a += 64;
    }

    int bi = 8 * i;
    a = &(data->A[64 * diag[i]]);

    y[bi] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4 +
            a[5] * y5 + a[6] * y6 + a[7] * y7;
    y[bi + 1] = a[8] * y0 + a[9] * y1 + a[10] * y2 + a[11] * y3 + a[12] * y4 +
                a[13] * y5 + a[14] * y6 + a[15] * y7;
    y[bi + 2] = a[16] * y0 + a[17] * y1 + a[18] * y2 + a[19] * y3 + a[20] * y4 +
                a[21] * y5 + a[22] * y6 + a[23] * y7;
    y[bi + 3] = a[24] * y0 + a[25] * y1 + a[26] * y2 + a[27] * y3 + a[28] * y4 +
                a[29] * y5 + a[30] * y6 + a[31] * y7;
    y[bi + 4] = a[32] * y0 + a[33] * y1 + a[34] * y2 + a[35] * y3 + a[36] * y4 +
                a[37] * y5 + a[38] * y6 + a[39] * y7;
    y[bi + 5] = a[40] * y0 + a[41] * y1 + a[42] * y2 + a[43] * y3 + a[44] * y4 +
                a[45] * y5 + a[46] * y6 + a[47] * y7;
    y[bi + 6] = a[48] * y0 + a[49] * y1 + a[50] * y2 + a[51] * y3 + a[52] * y4 +
                a[53] * y5 + a[54] * y6 + a[55] * y7;
    y[bi + 7] = a[56] * y0 + a[57] * y1 + a[58] * y2 + a[59] * y3 + a[60] * y4 +
                a[61] * y5 + a[62] * y6 + a[63] * y7;

    x -= 8;
    TacsAddFlops(2 * 64 * nz + 120);
  }
}

/*!
  Apply a portion of the lower factorization x = L^{-1} x
*/

void BCSRMatApplyPartialLower8(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar *xx = &x[6];
  int off = 6 * var_offset;

  for (int i = var_offset + 1; i < nrows; i++) {
    int end = diag[i];
    int k = rowp[i];
    while (cols[k] < var_offset) k++;

    const TacsScalar *a = &A[36 * k];
    for (; k < end; k++) {
      int j = 6 * cols[k] - off;

      xx[0] -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] +
               a[3] * x[j + 3] + a[4] * x[j + 4] + a[5] * x[j + 5];
      xx[1] -= a[6] * x[j] + a[7] * x[j + 1] + a[8] * x[j + 2] +
               a[9] * x[j + 3] + a[10] * x[j + 4] + a[11] * x[j + 5];
      xx[2] -= a[12] * x[j] + a[13] * x[j + 1] + a[14] * x[j + 2] +
               a[15] * x[j + 3] + a[16] * x[j + 4] + a[17] * x[j + 5];
      xx[3] -= a[18] * x[j] + a[19] * x[j + 1] + a[20] * x[j + 2] +
               a[21] * x[j + 3] + a[22] * x[j + 4] + a[23] * x[j + 5];
      xx[4] -= a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
               a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5];
      xx[5] -= a[30] * x[j] + a[31] * x[j + 1] + a[32] * x[j + 2] +
               a[33] * x[j + 3] + a[34] * x[j + 4] + a[35] * x[j + 5];
      a += 36;
    }

    xx += 6;
    TacsAddFlops(2 * 36 * nz);
  }
}

/*!
  Apply a portion of he upper factorization x = U^{-1} x
*/

void BCSRMatApplyPartialUpper8(BCSRMatData *data, TacsScalar *x,
                               int var_offset) {
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y0, y1, y2, y3, y4, y5, y6, y7;
  TacsScalar *xx = &x[8 * (nrows - var_offset - 1)];
  int off = 8 * var_offset;

  for (int i = nrows - 1; i >= var_offset; i--) {
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];
    y4 = xx[4];
    y5 = xx[5];
    y6 = xx[6];
    y7 = xx[7];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    const TacsScalar *a = &A[64 * k];
    for (; k < end; k++) {
      int j = 8 * cols[k] - off;

      y0 -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] + a[3] * x[j + 3] +
            a[4] * x[j + 4] + a[5] * x[j + 5] + a[6] * x[j + 6] +
            a[7] * x[j + 7];
      y1 -= a[8] * x[j] + a[9] * x[j + 1] + a[10] * x[j + 2] +
            a[11] * x[j + 3] + a[12] * x[j + 4] + a[13] * x[j + 5] +
            a[14] * x[j + 6] + a[15] * x[j + 7];
      y2 -= a[16] * x[j] + a[17] * x[j + 1] + a[18] * x[j + 2] +
            a[19] * x[j + 3] + a[20] * x[j + 4] + a[21] * x[j + 5] +
            a[22] * x[j + 6] + a[23] * x[j + 7];
      y3 -= a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
            a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5] +
            a[30] * x[j + 6] + a[31] * x[j + 7];
      y4 -= a[32] * x[j] + a[33] * x[j + 1] + a[34] * x[j + 2] +
            a[35] * x[j + 3] + a[36] * x[j + 4] + a[37] * x[j + 5] +
            a[38] * x[j + 6] + a[39] * x[j + 7];
      y5 -= a[40] * x[j] + a[41] * x[j + 1] + a[42] * x[j + 2] +
            a[43] * x[j + 3] + a[44] * x[j + 4] + a[45] * x[j + 5] +
            a[46] * x[j + 6] + a[47] * x[j + 7];
      y6 -= a[48] * x[j] + a[49] * x[j + 1] + a[50] * x[j + 2] +
            a[51] * x[j + 3] + a[52] * x[j + 4] + a[53] * x[j + 5] +
            a[54] * x[j + 6] + a[55] * x[j + 7];
      y7 -= a[56] * x[j] + a[57] * x[j + 1] + a[58] * x[j + 2] +
            a[59] * x[j + 3] + a[60] * x[j + 4] + a[61] * x[j + 5] +
            a[62] * x[j + 6] + a[63] * x[j + 7];
      a += 64;
    }

    a = &A[64 * diag[i]];
    xx[0] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4 +
            a[5] * y5 + a[6] * y6 + a[7] * y7;
    xx[1] = a[8] * y0 + a[9] * y1 + a[10] * y2 + a[11] * y3 + a[12] * y4 +
            a[13] * y5 + a[14] * y6 + a[15] * y7;
    xx[2] = a[16] * y0 + a[17] * y1 + a[18] * y2 + a[19] * y3 + a[20] * y4 +
            a[21] * y5 + a[22] * y6 + a[23] * y7;
    xx[3] = a[24] * y0 + a[25] * y1 + a[26] * y2 + a[27] * y3 + a[28] * y4 +
            a[29] * y5 + a[30] * y6 + a[31] * y7;
    xx[4] = a[32] * y0 + a[33] * y1 + a[34] * y2 + a[35] * y3 + a[36] * y4 +
            a[37] * y5 + a[38] * y6 + a[39] * y7;
    xx[5] = a[40] * y0 + a[41] * y1 + a[42] * y2 + a[43] * y3 + a[44] * y4 +
            a[45] * y5 + a[46] * y6 + a[47] * y7;
    xx[6] = a[48] * y0 + a[49] * y1 + a[50] * y2 + a[51] * y3 + a[52] * y4 +
            a[53] * y5 + a[54] * y6 + a[55] * y7;
    xx[7] = a[56] * y0 + a[57] * y1 + a[58] * y2 + a[59] * y3 + a[60] * y4 +
            a[61] * y5 + a[62] * y6 + a[63] * y7;
    xx -= 8;
    TacsAddFlops(2 * 64 * nz + 120);
  }
}

/*!
  Function for the approximate Schur preconditioner
*/

void BCSRMatApplyFactorSchur8(BCSRMatData *data, TacsScalar *x,
                              int var_offset) {
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  TacsScalar y0, y1, y2, y3, y4, y5, y6, y7;
  TacsScalar *xx = &x[6 * (var_offset - 1)];

  for (int i = var_offset - 1; i >= 0; i--) {
    y0 = xx[0];
    y1 = xx[1];
    y2 = xx[2];
    y3 = xx[3];
    y4 = xx[4];
    y5 = xx[5];
    y6 = xx[6];
    y7 = xx[7];

    int end = rowp[i + 1];
    int k = diag[i] + 1;
    const TacsScalar *a = &A[64 * k];
    for (; k < end; k++) {
      int j = 8 * cols[k];

      y0 -= a[0] * x[j] + a[1] * x[j + 1] + a[2] * x[j + 2] + a[3] * x[j + 3] +
            a[4] * x[j + 4] + a[5] * x[j + 5] + a[6] * x[j + 6] +
            a[7] * x[j + 7];
      y1 -= a[8] * x[j] + a[9] * x[j + 1] + a[10] * x[j + 2] +
            a[11] * x[j + 3] + a[12] * x[j + 4] + a[13] * x[j + 5] +
            a[14] * x[j + 6] + a[15] * x[j + 7];
      y2 -= a[16] * x[j] + a[17] * x[j + 1] + a[18] * x[j + 2] +
            a[19] * x[j + 3] + a[20] * x[j + 4] + a[21] * x[j + 5] +
            a[22] * x[j + 6] + a[23] * x[j + 7];
      y3 -= a[24] * x[j] + a[25] * x[j + 1] + a[26] * x[j + 2] +
            a[27] * x[j + 3] + a[28] * x[j + 4] + a[29] * x[j + 5] +
            a[30] * x[j + 6] + a[31] * x[j + 7];
      y4 -= a[32] * x[j] + a[33] * x[j + 1] + a[34] * x[j + 2] +
            a[35] * x[j + 3] + a[36] * x[j + 4] + a[37] * x[j + 5] +
            a[38] * x[j + 6] + a[39] * x[j + 7];
      y5 -= a[40] * x[j] + a[41] * x[j + 1] + a[42] * x[j + 2] +
            a[43] * x[j + 3] + a[44] * x[j + 4] + a[45] * x[j + 5] +
            a[46] * x[j + 6] + a[47] * x[j + 7];
      y6 -= a[48] * x[j] + a[49] * x[j + 1] + a[50] * x[j + 2] +
            a[51] * x[j + 3] + a[52] * x[j + 4] + a[53] * x[j + 5] +
            a[54] * x[j + 6] + a[55] * x[j + 7];
      y7 -= a[56] * x[j] + a[57] * x[j + 1] + a[58] * x[j + 2] +
            a[59] * x[j + 3] + a[60] * x[j + 4] + a[61] * x[j + 5] +
            a[62] * x[j + 6] + a[63] * x[j + 7];
      a += 64;
    }

    a = &A[64 * diag[i]];
    xx[0] = a[0] * y0 + a[1] * y1 + a[2] * y2 + a[3] * y3 + a[4] * y4 +
            a[5] * y5 + a[6] * y6 + a[7] * y7;
    xx[1] = a[8] * y0 + a[9] * y1 + a[10] * y2 + a[11] * y3 + a[12] * y4 +
            a[13] * y5 + a[14] * y6 + a[15] * y7;
    xx[2] = a[16] * y0 + a[17] * y1 + a[18] * y2 + a[19] * y3 + a[20] * y4 +
            a[21] * y5 + a[22] * y6 + a[23] * y7;
    xx[3] = a[24] * y0 + a[25] * y1 + a[26] * y2 + a[27] * y3 + a[28] * y4 +
            a[29] * y5 + a[30] * y6 + a[31] * y7;
    xx[4] = a[32] * y0 + a[33] * y1 + a[34] * y2 + a[35] * y3 + a[36] * y4 +
            a[37] * y5 + a[38] * y6 + a[39] * y7;
    xx[5] = a[40] * y0 + a[41] * y1 + a[42] * y2 + a[43] * y3 + a[44] * y4 +
            a[45] * y5 + a[46] * y6 + a[47] * y7;
    xx[6] = a[48] * y0 + a[49] * y1 + a[50] * y2 + a[51] * y3 + a[52] * y4 +
            a[53] * y5 + a[54] * y6 + a[55] * y7;
    xx[7] = a[56] * y0 + a[57] * y1 + a[58] * y2 + a[59] * y3 + a[60] * y4 +
            a[61] * y5 + a[62] * y6 + a[63] * y7;
    xx -= 8;
    TacsAddFlops(2 * 64 * nz + 120);
  }
}

/*
  Perform the matrix-matrix multiplication in parallel using pthreads
*/
void *BCSRMatMatMultAdd8_thread(void *t) {
  BCSRMatThread *tdata = static_cast<BCSRMatThread *>(t);
  const double alpha = tdata->alpha;

  // Retrieve the data required from the matrix
  const int nrows_a = tdata->Amat->nrows;
  const int *arowp = tdata->Amat->rowp;
  const int *acols = tdata->Amat->cols;
  const TacsScalar *A = tdata->Amat->A;
  const int group_size = tdata->mat->matmat_group_size;

  const int *browp = tdata->Bmat->rowp;
  const int *bcols = tdata->Bmat->cols;
  const TacsScalar *B = tdata->Bmat->A;

  // The matrix being written to
  const int nrows_c = tdata->mat->nrows;
  const int *crowp = tdata->mat->rowp;
  const int *ccols = tdata->mat->cols;
  TacsScalar *C = tdata->mat->A;

  while (tdata->num_completed_rows < nrows_c) {
    int row = -1;
    tdata->mat_mult_sched_job(group_size, &row);

    if (row < 0) {
      break;
    }

    if (alpha == 1.0) {
      // C_{ik} = A_{ij} B_{jk}
      for (int i = row; (i < nrows_a) && (i < row + group_size); i++) {
        for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
          int j = acols[jp];
          const TacsScalar *a = &A[64 * jp];

          int kp = browp[j];
          int kp_end = browp[j + 1];
          const TacsScalar *b = &B[64 * kp];

          int cp = crowp[i];
          int cp_end = crowp[i + 1];
          TacsScalar *c = &C[64 * cp];

          for (; kp < kp_end; kp++) {
            while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
              cp++;
              c += 64;
            }
            if (cp >= cp_end) {
              break;
            }

            if (bcols[kp] == ccols[cp]) {
              // Compute the matrix-matrix multiplication
              TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

              // first
              b0 = b[0];
              b1 = b[8];
              b2 = b[16];
              b3 = b[24];
              b4 = b[32];
              b5 = b[40];
              b6 = b[48];
              b7 = b[56];
              c[0] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[8] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[16] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[24] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[32] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[40] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[48] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[56] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // second
              b0 = b[1];
              b1 = b[9];
              b2 = b[17];
              b3 = b[25];
              b4 = b[33];
              b5 = b[41];
              b6 = b[49];
              b7 = b[57];
              c[1] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[9] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[17] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[25] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[33] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[41] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[49] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[57] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // third
              b0 = b[2];
              b1 = b[10];
              b2 = b[18];
              b3 = b[26];
              b4 = b[34];
              b5 = b[42];
              b6 = b[50];
              b7 = b[58];
              c[2] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[10] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[18] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[26] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[34] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[42] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[50] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[58] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // four
              b0 = b[3];
              b1 = b[11];
              b2 = b[19];
              b3 = b[27];
              b4 = b[35];
              b5 = b[43];
              b6 = b[51];
              b7 = b[59];
              c[3] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[11] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[19] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[27] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[35] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[43] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[51] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[59] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // five
              b0 = b[4];
              b1 = b[12];
              b2 = b[20];
              b3 = b[28];
              b4 = b[36];
              b5 = b[44];
              b6 = b[52];
              b7 = b[60];
              c[4] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[12] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[20] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[28] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[36] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[44] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[52] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[60] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // sixth
              b0 = b[5];
              b1 = b[13];
              b2 = b[21];
              b3 = b[29];
              b4 = b[37];
              b5 = b[45];
              b6 = b[53];
              b7 = b[61];
              c[5] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[13] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[21] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[29] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[37] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[45] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[53] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[61] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // seventh
              b0 = b[6];
              b1 = b[14];
              b2 = b[22];
              b3 = b[30];
              b4 = b[38];
              b5 = b[46];
              b6 = b[54];
              b7 = b[62];
              c[6] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[14] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[22] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[30] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[38] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[46] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[54] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[62] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // eight
              b0 = b[7];
              b1 = b[15];
              b2 = b[23];
              b3 = b[31];
              b4 = b[39];
              b5 = b[47];
              b6 = b[55];
              b7 = b[63];
              c[7] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[15] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[23] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[31] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[39] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[47] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[55] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[63] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;
            }
            b += 64;
          }
        }
      }
    } else if (alpha == -1.0) {
      // C_{ik} = A_{ij} B_{jk}
      for (int i = row; (i < nrows_a) && (i < row + group_size); i++) {
        for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
          int j = acols[jp];
          const TacsScalar *a = &A[64 * jp];

          int kp = browp[j];
          int kp_end = browp[j + 1];
          const TacsScalar *b = &B[64 * kp];

          int cp = crowp[i];
          int cp_end = crowp[i + 1];
          TacsScalar *c = &C[64 * cp];

          for (; kp < kp_end; kp++) {
            while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
              cp++;
              c += 64;
            }
            if (cp >= cp_end) {
              break;
            }

            if (bcols[kp] == ccols[cp]) {
              // Compute the matrix-matrix multiplication
              TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

              // first
              b0 = b[0];
              b1 = b[8];
              b2 = b[16];
              b3 = b[24];
              b4 = b[32];
              b5 = b[40];
              b6 = b[48];
              b7 = b[56];
              c[0] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[8] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[16] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[24] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[32] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[40] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[48] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[56] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // second
              b0 = b[1];
              b1 = b[9];
              b2 = b[17];
              b3 = b[25];
              b4 = b[33];
              b5 = b[41];
              b6 = b[49];
              b7 = b[57];
              c[1] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[9] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[17] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[25] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[33] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[41] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[49] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[57] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // third
              b0 = b[2];
              b1 = b[10];
              b2 = b[18];
              b3 = b[26];
              b4 = b[34];
              b5 = b[42];
              b6 = b[50];
              b7 = b[58];
              c[2] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[10] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[18] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[26] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[34] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[42] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[50] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[58] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // four
              b0 = b[3];
              b1 = b[11];
              b2 = b[19];
              b3 = b[27];
              b4 = b[35];
              b5 = b[43];
              b6 = b[51];
              b7 = b[59];
              c[3] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[11] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[19] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[27] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[35] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[43] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[51] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[59] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // five
              b0 = b[4];
              b1 = b[12];
              b2 = b[20];
              b3 = b[28];
              b4 = b[36];
              b5 = b[44];
              b6 = b[52];
              b7 = b[60];
              c[4] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[12] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[20] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[28] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[36] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[44] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[52] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[60] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // sixth
              b0 = b[5];
              b1 = b[13];
              b2 = b[21];
              b3 = b[29];
              b4 = b[37];
              b5 = b[45];
              b6 = b[53];
              b7 = b[61];
              c[5] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[13] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[21] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[29] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[37] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[45] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[53] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[61] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // seventh
              b0 = b[6];
              b1 = b[14];
              b2 = b[22];
              b3 = b[30];
              b4 = b[38];
              b5 = b[46];
              b6 = b[54];
              b7 = b[62];
              c[6] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[14] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[22] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[30] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[38] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[46] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[54] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[62] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // eight
              b0 = b[7];
              b1 = b[15];
              b2 = b[23];
              b3 = b[31];
              b4 = b[39];
              b5 = b[47];
              b6 = b[55];
              b7 = b[63];
              c[7] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[15] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[23] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[31] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[39] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[47] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[55] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[63] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;
            }

            b += 64;
          }
        }
      }
    } else {
      // C_{ik} = A_{ij} B_{jk}
      for (int i = row; (i < nrows_a) && (i < row + group_size); i++) {
        for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
          int j = acols[jp];
          const TacsScalar *a = &A[36 * jp];

          int kp = browp[j];
          int kp_end = browp[j + 1];
          const TacsScalar *b = &B[36 * kp];

          int cp = crowp[i];
          int cp_end = crowp[i + 1];
          TacsScalar *c = &C[36 * cp];

          for (; kp < kp_end; kp++) {
            while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
              cp++;
              c += 36;
            }
            if (cp >= cp_end) {
              break;
            }

            if (bcols[kp] == ccols[cp]) {
              // Compute the matrix-matrix multiplication
              TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

              // first
              b0 = b[0];
              b1 = b[8];
              b2 = b[16];
              b3 = b[24];
              b4 = b[32];
              b5 = b[40];
              b6 = b[48];
              b7 = b[56];
              c[0] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[8] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[16] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[24] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[32] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[40] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[48] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[56] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // second
              b0 = b[1];
              b1 = b[9];
              b2 = b[17];
              b3 = b[25];
              b4 = b[33];
              b5 = b[41];
              b6 = b[49];
              b7 = b[57];
              c[1] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[9] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                      a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[17] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[25] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[33] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[41] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[49] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[57] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // third
              b0 = b[2];
              b1 = b[10];
              b2 = b[18];
              b3 = b[26];
              b4 = b[34];
              b5 = b[42];
              b6 = b[50];
              b7 = b[58];
              c[2] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[10] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[18] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[26] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[34] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[42] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[50] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[58] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // four
              b0 = b[3];
              b1 = b[11];
              b2 = b[19];
              b3 = b[27];
              b4 = b[35];
              b5 = b[43];
              b6 = b[51];
              b7 = b[59];
              c[3] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[11] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[19] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[27] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[35] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[43] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[51] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[59] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // five
              b0 = b[4];
              b1 = b[12];
              b2 = b[20];
              b3 = b[28];
              b4 = b[36];
              b5 = b[44];
              b6 = b[52];
              b7 = b[60];
              c[4] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[12] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[20] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[28] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[36] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[44] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[52] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[60] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // sixth
              b0 = b[5];
              b1 = b[13];
              b2 = b[21];
              b3 = b[29];
              b4 = b[37];
              b5 = b[45];
              b6 = b[53];
              b7 = b[61];
              c[5] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[13] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[21] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[29] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[37] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[45] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[53] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[61] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // seventh
              b0 = b[6];
              b1 = b[14];
              b2 = b[22];
              b3 = b[30];
              b4 = b[38];
              b5 = b[46];
              b6 = b[54];
              b7 = b[62];
              c[6] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[14] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[22] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[30] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[38] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[46] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[54] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[62] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

              // eight
              b0 = b[7];
              b1 = b[15];
              b2 = b[23];
              b3 = b[31];
              b4 = b[39];
              b5 = b[47];
              b6 = b[55];
              b7 = b[63];
              c[7] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                      a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7;
              c[15] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                       a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
              c[23] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                       a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
              c[31] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                       a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
              c[39] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                       a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
              c[47] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                       a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
              c[55] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                       a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
              c[63] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                       a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;
            }

            b += 64;
          }
        }
      }
    }
  }

  pthread_exit(NULL);
}

/*!
  Perform a matrix-matrix multiplication
*/
void BCSRMatMatMultAdd8(double alpha, BCSRMatData *Adata, BCSRMatData *Bdata,
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
    int nz = 0;
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[64 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[64 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[64 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 64;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

            // first
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40];
            b6 = b[48];
            b7 = b[56];
            c[0] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[8] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                    a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[16] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[24] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[32] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[40] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[48] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[56] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // second
            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41];
            b6 = b[49];
            b7 = b[57];
            c[1] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[9] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                    a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[17] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[25] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[33] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[41] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[49] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[57] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // third
            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42];
            b6 = b[50];
            b7 = b[58];
            c[2] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[10] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[18] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[26] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[34] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[42] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[50] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[58] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // four
            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43];
            b6 = b[51];
            b7 = b[59];
            c[3] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[11] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[19] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[27] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[35] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[43] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[51] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[59] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // five
            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44];
            b6 = b[52];
            b7 = b[60];
            c[4] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[12] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[20] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[28] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[36] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[44] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[52] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[60] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // sixth
            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45];
            b6 = b[53];
            b7 = b[61];
            c[5] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[13] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[21] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[29] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[37] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[45] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[53] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[61] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // seventh
            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46];
            b6 = b[54];
            b7 = b[62];
            c[6] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[14] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[22] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[30] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[38] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[46] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[54] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[62] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // eight
            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47];
            b6 = b[55];
            b7 = b[63];
            c[7] += a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[15] += a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[23] += a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[31] += a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[39] += a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[47] += a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[55] += a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[63] += a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            nz++;
          }
          b += 64;
        }
      }
    }
    TacsAddFlops(2 * 64 * 8 * nz);
  } else if (alpha == -1.0) {
    // C_{ik} = A_{ij} B_{jk}
    int nz = 0;
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[64 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[64 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[64 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 64;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

            // first
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40];
            b6 = b[48];
            b7 = b[56];
            c[0] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[8] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                    a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[16] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[24] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[32] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[40] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[48] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[56] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // second
            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41];
            b6 = b[49];
            b7 = b[57];
            c[1] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[9] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                    a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[17] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[25] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[33] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[41] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[49] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[57] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // third
            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42];
            b6 = b[50];
            b7 = b[58];
            c[2] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[10] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[18] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[26] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[34] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[42] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[50] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[58] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // four
            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43];
            b6 = b[51];
            b7 = b[59];
            c[3] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[11] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[19] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[27] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[35] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[43] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[51] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[59] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // five
            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44];
            b6 = b[52];
            b7 = b[60];
            c[4] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[12] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[20] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[28] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[36] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[44] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[52] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[60] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // sixth
            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45];
            b6 = b[53];
            b7 = b[61];
            c[5] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[13] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[21] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[29] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[37] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[45] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[53] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[61] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // seventh
            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46];
            b6 = b[54];
            b7 = b[62];
            c[6] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[14] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[22] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[30] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[38] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[46] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[54] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[62] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            // eight
            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47];
            b6 = b[55];
            b7 = b[63];
            c[7] -= a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 + a[4] * b4 +
                    a[5] * b5 + a[6] * b6 + a[7] * b7;
            c[15] -= a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                     a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7;
            c[23] -= a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                     a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7;
            c[31] -= a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                     a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7;
            c[39] -= a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                     a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7;
            c[47] -= a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                     a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7;
            c[55] -= a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                     a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7;
            c[63] -= a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                     a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7;

            nz++;
          }

          b += 64;
        }
      }
    }
    TacsAddFlops(2 * 64 * 8 * nz);
  } else {
    // C_{ik} = A_{ij} B_{jk}
    int nz = 0;
    for (int i = 0; i < nrows_a; i++) {
      for (int jp = arowp[i]; jp < arowp[i + 1]; jp++) {
        int j = acols[jp];
        const TacsScalar *a = &A[64 * jp];

        int kp = browp[j];
        int kp_end = browp[j + 1];
        const TacsScalar *b = &B[64 * kp];

        int cp = crowp[i];
        int cp_end = crowp[i + 1];
        TacsScalar *c = &C[64 * cp];

        for (; kp < kp_end; kp++) {
          while ((cp < cp_end) && (ccols[cp] < bcols[kp])) {
            cp++;
            c += 64;
          }
          if (cp >= cp_end) {
            break;
          }

          if (bcols[kp] == ccols[cp]) {
            // Compute the matrix-matrix multiplication
            TacsScalar b0, b1, b2, b3, b4, b5, b6, b7;

            // first
            b0 = b[0];
            b1 = b[8];
            b2 = b[16];
            b3 = b[24];
            b4 = b[32];
            b5 = b[40];
            b6 = b[48];
            b7 = b[56];
            c[0] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[8] += alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                             a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[16] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[24] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[32] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[40] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[48] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[56] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // second
            b0 = b[1];
            b1 = b[9];
            b2 = b[17];
            b3 = b[25];
            b4 = b[33];
            b5 = b[41];
            b6 = b[49];
            b7 = b[57];
            c[1] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[9] += alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                             a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[17] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[25] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[33] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[41] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[49] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[57] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // third
            b0 = b[2];
            b1 = b[10];
            b2 = b[18];
            b3 = b[26];
            b4 = b[34];
            b5 = b[42];
            b6 = b[50];
            b7 = b[58];
            c[2] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[10] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[18] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[26] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[34] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[42] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[50] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[58] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // four
            b0 = b[3];
            b1 = b[11];
            b2 = b[19];
            b3 = b[27];
            b4 = b[35];
            b5 = b[43];
            b6 = b[51];
            b7 = b[59];
            c[3] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[11] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[19] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[27] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[35] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[43] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[51] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[59] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // five
            b0 = b[4];
            b1 = b[12];
            b2 = b[20];
            b3 = b[28];
            b4 = b[36];
            b5 = b[44];
            b6 = b[52];
            b7 = b[60];
            c[4] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[12] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[20] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[28] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[36] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[44] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[52] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[60] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // sixth
            b0 = b[5];
            b1 = b[13];
            b2 = b[21];
            b3 = b[29];
            b4 = b[37];
            b5 = b[45];
            b6 = b[53];
            b7 = b[61];
            c[5] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[13] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[21] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[29] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[37] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[45] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[53] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[61] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // seventh
            b0 = b[6];
            b1 = b[14];
            b2 = b[22];
            b3 = b[30];
            b4 = b[38];
            b5 = b[46];
            b6 = b[54];
            b7 = b[62];
            c[6] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[14] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[22] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[30] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[38] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[46] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[54] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[62] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);

            // eight
            b0 = b[7];
            b1 = b[15];
            b2 = b[23];
            b3 = b[31];
            b4 = b[39];
            b5 = b[47];
            b6 = b[55];
            b7 = b[63];
            c[7] += alpha * (a[0] * b0 + a[1] * b1 + a[2] * b2 + a[3] * b3 +
                             a[4] * b4 + a[5] * b5 + a[6] * b6 + a[7] * b7);
            c[15] +=
                alpha * (a[8] * b0 + a[9] * b1 + a[10] * b2 + a[11] * b3 +
                         a[12] * b4 + a[13] * b5 + a[14] * b6 + a[15] * b7);
            c[23] +=
                alpha * (a[16] * b0 + a[17] * b1 + a[18] * b2 + a[19] * b3 +
                         a[20] * b4 + a[21] * b5 + a[22] * b6 + a[23] * b7);
            c[31] +=
                alpha * (a[24] * b0 + a[25] * b1 + a[26] * b2 + a[27] * b3 +
                         a[28] * b4 + a[29] * b5 + a[30] * b6 + a[31] * b7);
            c[39] +=
                alpha * (a[32] * b0 + a[33] * b1 + a[34] * b2 + a[35] * b3 +
                         a[36] * b4 + a[37] * b5 + a[38] * b6 + a[39] * b7);
            c[47] +=
                alpha * (a[40] * b0 + a[41] * b1 + a[42] * b2 + a[43] * b3 +
                         a[44] * b4 + a[45] * b5 + a[46] * b6 + a[47] * b7);
            c[55] +=
                alpha * (a[48] * b0 + a[49] * b1 + a[50] * b2 + a[51] * b3 +
                         a[52] * b4 + a[53] * b5 + a[54] * b6 + a[55] * b7);
            c[63] +=
                alpha * (a[56] * b0 + a[57] * b1 + a[58] * b2 + a[59] * b3 +
                         a[60] * b4 + a[61] * b5 + a[62] * b6 + a[63] * b7);
            nz++;
          }
          b += 64;
        }
      }
    }

    TacsAddFlops((2 * 64 * 8 + 64) * nz);
  }
}

/*!
  Apply an iteration of SOR to the system A*x = b.
*/
void BCSRMatApplySOR8(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
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
  TacsScalar t1, t2, t3, t4, t5, t6, t7, t8;

  // Go through the matrix with the forward ordering
  for (int i = start; i < end; i++) {
    // Copy the right-hand-side to the temporary vector for this row
    t1 = b[8 * i];
    t2 = b[8 * i + 1];
    t3 = b[8 * i + 2];
    t4 = b[8 * i + 3];
    t5 = b[8 * i + 4];
    t6 = b[8 * i + 5];
    t7 = b[8 * i + 6];
    t8 = b[8 * i + 7];

    // Set the pointer to the beginning of the current row
    const TacsScalar *a = &Adata->A[64 * Arowp[i]];

    // Scan through the row and compute the result:
    // tx <- b_i - A_{ij}*x_{j} for j != i
    int end = Arowp[i + 1];
    for (int k = Arowp[i]; k < end; k++) {
      int j = Acols[k];
      TacsScalar *y = &x[8 * j];

      if (i != j) {
        t1 -= a[0] * y[0] + a[1] * y[1] + a[2] * y[2] + a[3] * y[3] +
              a[4] * y[4] + a[5] * y[5] + a[6] * y[6] + a[7] * y[7];
        t2 -= a[8] * y[0] + a[9] * y[1] + a[10] * y[2] + a[11] * y[3] +
              a[12] * y[4] + a[13] * y[5] + a[14] * y[6] + a[15] * y[7];
        t3 -= a[16] * y[0] + a[17] * y[1] + a[18] * y[2] + a[19] * y[3] +
              a[20] * y[4] + a[21] * y[5] + a[22] * y[6] + a[23] * y[7];
        t4 -= a[24] * y[0] + a[25] * y[1] + a[26] * y[2] + a[27] * y[3] +
              a[28] * y[4] + a[29] * y[5] + a[30] * y[6] + a[31] * y[7];
        t5 -= a[32] * y[0] + a[33] * y[1] + a[34] * y[2] + a[35] * y[3] +
              a[36] * y[4] + a[37] * y[5] + a[38] * y[6] + a[39] * y[7];
        t6 -= a[40] * y[0] + a[41] * y[1] + a[42] * y[2] + a[43] * y[3] +
              a[44] * y[4] + a[45] * y[5] + a[46] * y[6] + a[47] * y[7];
        t7 -= a[48] * y[0] + a[49] * y[1] + a[50] * y[2] + a[51] * y[3] +
              a[52] * y[4] + a[53] * y[5] + a[54] * y[6] + a[55] * y[7];
        t8 -= a[56] * y[0] + a[57] * y[1] + a[58] * y[2] + a[59] * y[3] +
              a[60] * y[4] + a[61] * y[5] + a[62] * y[6] + a[63] * y[7];
      }

      // Increment the block pointer by bsize^2
      a += 64;
    }

    if (Bdata && i >= var_offset) {
      const int row = i - var_offset;

      // Set the pointer to the row in B
      a = &Bdata->A[64 * Browp[row]];
      end = Browp[row + 1];
      for (int k = Browp[row]; k < end; k++) {
        int j = Bcols[k];
        const TacsScalar *y = &xext[8 * j];

        t1 -= a[0] * y[0] + a[1] * y[1] + a[2] * y[2] + a[3] * y[3] +
              a[4] * y[4] + a[5] * y[5] + a[6] * y[6] + a[7] * y[7];
        t2 -= a[8] * y[0] + a[9] * y[1] + a[10] * y[2] + a[11] * y[3] +
              a[12] * y[4] + a[13] * y[5] + a[14] * y[6] + a[15] * y[7];
        t3 -= a[16] * y[0] + a[17] * y[1] + a[18] * y[2] + a[19] * y[3] +
              a[20] * y[4] + a[21] * y[5] + a[22] * y[6] + a[23] * y[7];
        t4 -= a[24] * y[0] + a[25] * y[1] + a[26] * y[2] + a[27] * y[3] +
              a[28] * y[4] + a[29] * y[5] + a[30] * y[6] + a[31] * y[7];
        t5 -= a[32] * y[0] + a[33] * y[1] + a[34] * y[2] + a[35] * y[3] +
              a[36] * y[4] + a[37] * y[5] + a[38] * y[6] + a[39] * y[7];
        t6 -= a[40] * y[0] + a[41] * y[1] + a[42] * y[2] + a[43] * y[3] +
              a[44] * y[4] + a[45] * y[5] + a[46] * y[6] + a[47] * y[7];
        t7 -= a[48] * y[0] + a[49] * y[1] + a[50] * y[2] + a[51] * y[3] +
              a[52] * y[4] + a[53] * y[5] + a[54] * y[6] + a[55] * y[7];
        t8 -= a[56] * y[0] + a[57] * y[1] + a[58] * y[2] + a[59] * y[3] +
              a[60] * y[4] + a[61] * y[5] + a[62] * y[6] + a[63] * y[7];

        a += 64;
      }
    }

    // Set a pointer to the inverse of the diagonal
    const TacsScalar *d = &Adiag[64 * i];

    // Compute the first term in the update:
    // x[i] = (1.0 - omega)*x[i] + omega*D^{-1}tx
    x[8 * i] = (1.0 - omega) * x[8 * i] +
               omega * (d[0] * t1 + d[1] * t2 + d[2] * t3 + d[3] * t4 +
                        d[4] * t5 + d[5] * t6 + d[6] * t7 + d[7] * t8);
    x[8 * i + 1] = (1.0 - omega) * x[8 * i + 1] +
                   omega * (d[8] * t1 + d[9] * t2 + d[10] * t3 + d[11] * t4 +
                            d[12] * t5 + d[13] * t6 + d[14] * t7 + d[15] * t8);
    x[8 * i + 2] = (1.0 - omega) * x[8 * i + 2] +
                   omega * (d[16] * t1 + d[17] * t2 + d[18] * t3 + d[19] * t4 +
                            d[20] * t5 + d[21] * t6 + d[22] * t7 + d[23] * t8);
    x[8 * i + 3] = (1.0 - omega) * x[8 * i + 3] +
                   omega * (d[24] * t1 + d[25] * t2 + d[26] * t3 + d[27] * t4 +
                            d[28] * t5 + d[29] * t6 + d[30] * t7 + d[31] * t8);
    x[8 * i + 4] = (1.0 - omega) * x[8 * i + 4] +
                   omega * (d[32] * t1 + d[33] * t2 + d[34] * t3 + d[35] * t4 +
                            d[36] * t5 + d[37] * t6 + d[38] * t7 + d[39] * t8);
    x[8 * i + 5] = (1.0 - omega) * x[8 * i + 5] +
                   omega * (d[40] * t1 + d[41] * t2 + d[42] * t3 + d[43] * t4 +
                            d[44] * t5 + d[45] * t6 + d[46] * t7 + d[47] * t8);
    x[8 * i + 6] = (1.0 - omega) * x[8 * i + 6] +
                   omega * (d[48] * t1 + d[49] * t2 + d[50] * t3 + d[51] * t4 +
                            d[52] * t5 + d[53] * t6 + d[54] * t7 + d[55] * t8);
    x[8 * i + 7] = (1.0 - omega) * x[8 * i + 7] +
                   omega * (d[56] * t1 + d[57] * t2 + d[58] * t3 + d[59] * t4 +
                            d[60] * t5 + d[61] * t6 + d[62] * t7 + d[63] * t8);
  }
}
