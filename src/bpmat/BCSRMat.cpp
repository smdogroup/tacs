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

#include "BCSRMat.h"

#include <stdio.h>

#include "BCSRMatImpl.h"
#include "TacsUtilities.h"
#include "tacslapack.h"

/*
  BCSR matrix implementation
*/

/*
  Merge two uniquely sorted arrays with levels associated with them.

  - All elements in a/alevs remain in a/alevs -- although they may be in
  different locations after the call.
  - Elements in b are only merged with a if blevs[j]+add_lev <= lev_fill
  - If a[i] == b[j] then the new level value will be min(alevs[i],blevs[j])
  - To complicate matters, b/blevs must be constant!

  1. Count the number of elements in b that will be discarded
  2. Count the number of duplicates
  3. Scan from the end of the whole array to the beginning
*/
static int mergeArraysWithLevels(int *a, int *alevs, int na, const int *b,
                                 const int *blevs, int nb, int add_lev,
                                 int lev_fill) {
  int ndup = 0;
  int ndiscard = 0;
  for (int i = 0; i < nb; i++) {
    if (blevs[i] + add_lev > lev_fill) {
      ndiscard++;
    }
  }

  int j = 0, i = 0;
  for (; i < na; i++) {
    while ((j < nb) && ((blevs[j] + add_lev > lev_fill) || (b[j] < a[i]))) {
      j++;
    }
    if (j >= nb) {
      break;
    }
    if (a[i] == b[j]) {
      ndup++;
    }
  }

  int len = na + nb - ndup - ndiscard;  // End of the array
  int end = len - 1;

  j = nb - 1;
  i = na - 1;
  while (i >= 0 && j >= 0) {
    while (j >= 0 && (blevs[j] + add_lev > lev_fill)) {
      j--;
    }

    if (j < 0) {
      break;
    }

    if (a[i] > b[j]) {
      a[end] = a[i];
      alevs[end] = alevs[i];
      end--, i--;
    } else if (b[j] > a[i]) {
      a[end] = b[j];
      alevs[end] = blevs[j] + add_lev;
      end--, j--;
    } else {  // b[j] == a[i]
      a[end] = a[i];
      // Take the minimum of the levels
      alevs[end] =
          (alevs[i] < (blevs[j] + add_lev) ? alevs[i] : (blevs[j] + add_lev));
      end--, j--, i--;
    }
  }

  // Only need to copy over remaining elements from b - if any
  if (end >= 0) {
    while (j >= 0) {
      while (j >= 0 && (blevs[j] + add_lev > lev_fill)) {
        j--;
      }

      if (j < 0) {
        break;
      }

      a[end] = b[j];
      alevs[end] = blevs[j] + add_lev;
      end--, j--;
    }
  }

  return len;
}

/*!
  The BCSR set up function with a set ILU fill level.

  Given the non-zero pattern of BCSRMat bmat, determine the non-zero
  ILU factorization.

  Input:
  matv    == the BCSR matrix that serves as the initial non-zero pattern
  levFill == the level of fill for the new matrix
  fill    == the approximate size of the new matrix

  1. Allocate space based on the parameter fill for the new non-zero pattern.
  This is approximate and need not be exact.

  2. Go through and perform a symbolic factorization to determine the
  new non-zero pattern.

  3. Allocate space for the non-zero pattern that is determined.
*/
BCSRMat::BCSRMat(MPI_Comm _comm, BCSRMat *mat, int levFill, double fill,
                 const char *fname) {
  comm = _comm;
  thread_info = mat->thread_info;
  thread_info->incref();
  tdata = NULL;
  Adiag = NULL;

  if (fill < 1.0) {
    fprintf(stderr, "BCSRMat(): fill must be greater than 1.0\n");
    fill = 2.0;
  }

  // Set the row map/bptr
  data = new BCSRMatData(mat->data->bsize, mat->data->nrows, mat->data->ncols);
  data->incref();
  initBlockImpl();

  int *levs;
  computeILUk(mat, levFill, fill, &levs);

  int length = data->rowp[data->nrows];
  int bsize = data->bsize;
  length *= bsize * bsize;

  data->A = new TacsScalar[length];
  memset(data->A, 0, length * sizeof(TacsScalar));

  // Go through and print out the nz-pattern of the matrix
  if (fname) {
    FILE *fp = fopen(fname, "w");
    if (fp) {
      fprintf(fp, "VARIABLES = \"i\",\"j\"\n");

      for (int k = 0; k < levFill; k++) {
        int header_flag = 1;

        for (int i = 0; i < data->nrows; i++) {
          for (int j = data->rowp[i]; j < data->rowp[i + 1]; j++) {
            if (k == levs[j]) {
              if (header_flag) {
                fprintf(fp, "ZONE T = \"Lev %d\"\n", k);
                header_flag = 0;
              }
              fprintf(fp, "%d %d\n", i, data->cols[j]);
            }
          }
        }
      }

      fclose(fp);
    }
  }

  delete[] levs;
}

/*!
  Create the BCSRMatrix from the rowp/cols arrays. This matrix may
  be rectangular. Note that the matrix retains ownership of the arrays
  rowp/cols by stealing the pointers and setting them to NULL.

  input:
  comm:         the communicator reference for this matrix
  thread_info:  the POSIX threads class for threading
  bsize:        the block size of the elements in the matrix
  nrows:        the number of rows in the matrix
  ncols:        the number of columns in the matrix
  rowp:         the CSR row pointer
  cols:         the column indices
  A:            the matrix values
*/
BCSRMat::BCSRMat(MPI_Comm _comm, TACSThreadInfo *_thread_info, int bsize,
                 int nrows, int ncols, int **_rowp, int **_cols,
                 TacsScalar **_A) {
  comm = _comm;
  thread_info = _thread_info;
  thread_info->incref();
  tdata = NULL;
  Adiag = NULL;

  data = new BCSRMatData(bsize, nrows, ncols);
  data->incref();
  initBlockImpl();

  // Take the pointer from the input
  data->rowp = *_rowp;
  data->cols = *_cols;
  *_rowp = NULL;
  *_cols = NULL;

  if (_A) {
    data->A = *_A;
    *_A = NULL;
  } else {
    // Find the size of the array
    int length = bsize * bsize * data->rowp[nrows];
    data->A = new TacsScalar[length];
    memset(data->A, 0, length * sizeof(TacsScalar));
  }
}

/*!
  This function performs a few calculations.

  Given the matrices with the following set up:
  [ B, E ]
  [ F    ]

  The lower diagonal entry is computed elsewhere
  with the Schur complement approach

  Compute the following:
  1. ILU factorization of L_{B} U_{B} = B + R
  2. The factor Epc = L_{B}^{-1} * E + Re
  3. The factor Fpc = F * U_{B}^{-1} + Rf
  4. Form the Schur complement system S = C - Fpc*Epc
*/
BCSRMat::BCSRMat(MPI_Comm _comm, BCSRMat *Bmat, BCSRMat *Emat, BCSRMat *Fmat,
                 BCSRMat *Cmat, int levFill, double fill, BCSRMat **Epc,
                 BCSRMat **Fpc, BCSRMat **Smat, int use_full_schur) {
  comm = _comm;
  thread_info = Bmat->thread_info;
  thread_info->incref();
  tdata = NULL;
  Adiag = NULL;

  // Check that the dimensions of the matrices match
  if (Bmat->data->nrows != Emat->data->nrows ||
      Bmat->data->ncols != Fmat->data->ncols) {
    fprintf(stderr, "BCSRMat error: Matrix dimensions do not agree\n");
    return;
  }

  // Set the row map/bptr
  data =
      new BCSRMatData(Bmat->data->bsize, Bmat->data->nrows, Bmat->data->ncols);
  data->incref();
  initBlockImpl();

  int *levs;
  computeILUk(Bmat, levFill, fill, &levs);

  if (use_full_schur) {
    // Don't use the levels, compute the full Schur complement
    *Epc = computeILUkEpc(Emat, levs, levFill, fill, NULL);
    *Fpc = computeILUkFpc(Fmat, levs, levFill, fill, NULL);
    *Smat = new BCSRMat(comm, Cmat, NULL, *Fpc, NULL, *Epc, levFill, fill);
  } else {
    // Use the levels from the factors Epc and Fpc in the Schur complement
    int *elevs, *flevs;
    *Epc = computeILUkEpc(Emat, levs, levFill, fill, &elevs);
    *Fpc = computeILUkFpc(Fmat, levs, levFill, fill, &flevs);
    *Smat = new BCSRMat(comm, Cmat, flevs, *Fpc, elevs, *Epc, levFill, fill);

    delete[] elevs;
    delete[] flevs;
  }

  delete[] levs;

  int bsize = data->bsize;
  int length = data->rowp[data->nrows];
  length *= bsize * bsize;
  data->A = new TacsScalar[length];
  memset(data->A, 0, length * sizeof(TacsScalar));
}

/*!
  Compute the non-zero pattern of the matrix C formed from the
  matrix multiplication of A and B. These matrices may be rectangular
  but the dimensions must agree - of course.

  'fill' is the expected fill in based on
  fill = nnz(this)/(nnz(A) + nnz(B) + nnz(S))

  If alevs and blevs are supplied then use them in conjunction with
  levFill to limit the number of non-zeros in the new matrix. If
  alevs and blevs are not supplied, ie. alevs = NULL, blevs = NULL,
  find all non-zero entries. In this case, levFill is not significant
*/
BCSRMat::BCSRMat(MPI_Comm _comm, BCSRMat *smat, int *alevs, BCSRMat *amat,
                 int *blevs, BCSRMat *bmat, int levFill, double fill) {
  comm = _comm;
  thread_info = smat->thread_info;
  thread_info->incref();
  tdata = NULL;
  Adiag = NULL;

  // Check that the block sizes are the same
  if (amat->data->bsize != bmat->data->bsize) {
    fprintf(stderr,
            "BCSRMat symbolic multiplication error: Matrix block "
            "sizes must be equal\n");
    return;
  }
  if (amat->data->nrows != bmat->data->ncols) {
    fprintf(stderr,
            "BCSRMat symbolic multiplication error: Matrix sizes "
            "must agree ncol(A) != nrow(B)\n");
    return;
  }

  data =
      new BCSRMatData(amat->data->bsize, amat->data->nrows, bmat->data->ncols);
  data->incref();
  initBlockImpl();

  const int nrows = data->nrows;

  // Allocate the row and column data
  int init_size = (amat->data->rowp[amat->data->nrows] +
                   bmat->data->rowp[bmat->data->nrows] +
                   smat->data->rowp[smat->data->nrows]);
  int max_size = (int)(fill * init_size);

  int *rowp = new int[nrows + 1];
  int *cols = new int[max_size];

  // C_{ik} = A_{ij}*B_{jk}
  int *tcols = new int[nrows];

  int nnz = 0;
  rowp[0] = 0;

  if (alevs && blevs) {
    int *tlevs = new int[nrows];
    memset(tlevs, 0, nrows * sizeof(int));

    for (int i = 0; i < amat->data->nrows; i++) {
      int num_cols = 0;  // The size of the temporary cols array

      // Include all the column indices from this row of S
      for (int jp = smat->data->rowp[i]; jp < smat->data->rowp[i + 1]; jp++) {
        tcols[num_cols] = smat->data->cols[jp];
        tlevs[num_cols] = 0;
        num_cols++;
      }

      for (int jp = amat->data->rowp[i]; jp < amat->data->rowp[i + 1]; jp++) {
        int j = amat->data->cols[jp];
        int alev = alevs[jp];

        // Merge the two arrays into cols
        int brpj = bmat->data->rowp[j];
        int brsize = bmat->data->rowp[j + 1] - brpj;
        num_cols = mergeArraysWithLevels(
            tcols, tlevs, num_cols, &(bmat->data->cols[brpj]), &blevs[brpj],
            brsize, alev + 1, levFill);
      }

      // Check if adding this row will exceed the allowable size
      if (nnz + num_cols > max_size) {
        max_size = 2 * max_size + num_cols;
        TacsExtendArray(&cols, nnz, max_size);
      }

      for (int k = 0; k < num_cols; k++, nnz++) {
        cols[nnz] = tcols[k];
      }
      rowp[i + 1] = nnz;
    }
  } else {
    // Compute the non-zero structure of the resulting matrix
    // nnz(S + A*B) one row at a time
    for (int i = 0; i < amat->data->nrows; i++) {
      int num_cols = 0;  // The size of the temporary cols array

      // Include all the column indices from this row of S
      for (int jp = smat->data->rowp[i]; jp < smat->data->rowp[i + 1]; jp++) {
        tcols[num_cols] = smat->data->cols[jp];
        num_cols++;
      }

      // Add the non-zero pattern to this matrix from each row of B
      // for each column of A. Merge the sorted arrays into a single
      // array - that will be at most length num cols in B (= num cols
      // in S). This produces a single sorted array
      for (int jp = amat->data->rowp[i]; jp < amat->data->rowp[i + 1]; jp++) {
        int j = amat->data->cols[jp];

        // Merge the two arrays into cols
        int brpj = bmat->data->rowp[j];
        int brsize = bmat->data->rowp[j + 1] - brpj;
        num_cols = TacsMergeSortedArrays(num_cols, tcols, brsize,
                                         &(bmat->data->cols[brpj]));
      }

      // Check if adding this row will exceed the allowable size
      if (nnz + num_cols > max_size) {
        max_size = 2 * max_size + num_cols;
        TacsExtendArray(&cols, nnz, max_size);
      }

      for (int k = 0; k < num_cols; k++, nnz++) {
        cols[nnz] = tcols[k];
      }
      rowp[i + 1] = nnz;
    }
  }

  // Clip the cols array to the correct size
  if (max_size > nnz) {
    TacsExtendArray(&cols, nnz, nnz);
  }

  delete[] tcols;

  data->rowp = rowp;
  data->cols = cols;

  int bsize = data->bsize;
  int length = rowp[nrows];
  length *= bsize * bsize;

  data->A = new TacsScalar[length];
  memset(data->A, 0, length * sizeof(TacsScalar));
}

/*
  This constructor assembles the BCSR matrix based on the non-zero
  pattern of the matrix computed as follows:

  A = B^{T}*B

  This represents the normal equations. The BCSRMat class has code
  for computing the scaled normal equations:

  A = B^{T}*S*B

  where S is a diagonal matrix. This is useful for some optimization
  methods. The non-zero pattern is computed using a variant of the
  following algorithm:

  Compute the non-zero pattern of A_{ij} = B_{ki}*B_{kj}

  for i = 1, n
  .  for k = 1, m
  .     for j = 1, n
  .        if (k, i) in nz(B) and (k, j) in nz(B)
  .        then add (i,j) to nz(A)

  Note that this computes the non-zero pattern of A one row at a time.
  This is more useful for the CSR-type format which is stored by row.

  input:
  comm:  the communicator
  B:     the matrix for
  fill:  the expected fill in
*/
BCSRMat::BCSRMat(MPI_Comm _comm, BCSRMat *B, double fill) {
  comm = _comm;
  thread_info = B->thread_info;
  thread_info->incref();
  tdata = NULL;
  Adiag = NULL;

  data = new BCSRMatData(B->data->bsize, B->data->ncols, B->data->ncols);
  data->incref();
  initBlockImpl();

  const int nrows = data->nrows;  // = B->data->ncols
  const int nrows_b = B->data->nrows;

  // Compute the non-zero pattern for each row
  int mat_size = B->data->rowp[nrows];
  int max_size = (int)(fill * mat_size);  // The maximum size - for now

  int *cols = new int[max_size];
  int *rowp = new int[nrows + 1];
  int *diag = new int[nrows];

  int *new_nz = new int[nrows];
  int *kptr = new int[nrows_b];
  rowp[0] = 0;
  memcpy(kptr, B->data->rowp, nrows_b * sizeof(int));

  // Go through each row and check the
  for (int i = 0; i < nrows; i++) {
    // Store the non-zero entries in this array. Indicate
    // non-zero values with a 1.
    memset(new_nz, 0, nrows * sizeof(int));

    // If (k, i) is in nz(B), then add all the terms in the
    // row (k, j) in nz(B)
    for (int k = 0; k < nrows_b; k++) {
      // Find the entries in column i
      if ((kptr[k] < B->data->rowp[k + 1]) && (B->data->cols[kptr[k]] == i)) {
        kptr[k]++;

        // Add the non-zero pattern from
        for (int jp = B->data->rowp[k]; jp < B->data->rowp[k + 1]; jp++) {
          new_nz[B->data->cols[jp]] = 1;
        }
      }
    }

    // Ensure that the diagonal entry exists
    new_nz[i] = 1;

    // Scan through all of the added elements, and add them to the
    // data structure if there is still room.
    int nz = 0;
    for (int k = 0; k < nrows; k++) {
      if (new_nz[k] == 1) {
        new_nz[nz] = k;

        // Set the diagonal entry
        if (k == i) {
          diag[i] = rowp[i] + nz;
        }

        // Increment the number of non-zeros in this row
        nz++;
      }
    }

    // Extend the array if required
    if (rowp[i] + nz > max_size) {
      max_size = 2 * max_size + nz;
      TacsExtendArray(&cols, rowp[i], max_size);
    }

    // Set the non-zeros from this row
    rowp[i + 1] = rowp[i] + nz;
    memcpy(&cols[rowp[i]], new_nz, nz * sizeof(int));
  }

  // Clip the cols array to the correct size
  int nnz = rowp[nrows];
  if (max_size > nnz) {
    TacsExtendArray(&cols, nnz, nnz);
  }

  delete[] new_nz;
  delete[] kptr;

  // Set the cols/rowp arrays into the data structure
  data->rowp = rowp;
  data->cols = cols;

  int bsize = data->bsize;
  int length = rowp[nrows];
  length *= bsize * bsize;

  data->A = new TacsScalar[length];
  memset(data->A, 0, length * sizeof(TacsScalar));
}

BCSRMat::~BCSRMat() {
  data->decref();
  thread_info->decref();
  if (tdata) {
    tdata->decref();
  }
  if (Adiag) {
    delete[] Adiag;
  }
}

MPI_Comm BCSRMat::getMPIComm() { return comm; }

TACSThreadInfo *BCSRMat::getThreadInfo() { return thread_info; }

/*!
  Compute the ILU(levFill) preconditioner

  fill == The expected degree of fill-in after the computation
  levs == The level set of the entry
*/
void BCSRMat::computeILUk(BCSRMat *mat, int levFill, double fill, int **_levs) {
  int nrows = mat->data->nrows;  // Record the number of rows/columns
  int ncols = mat->data->ncols;

  // Number of non-zeros in the original matrix
  int mat_size = mat->data->rowp[nrows];
  int size = 0;
  int max_size = (int)(fill * mat_size);  // The maximum size - for now

  int *cols = new int[max_size];
  int *levs = new int[max_size];  // The level of fill of an entry
  int *rowp = new int[nrows + 1];
  int *diag = new int[nrows];

  // Fill in the first entries
  rowp[0] = 0;

  // Allocate space for the temporary row info
  int *rlevs = new int[ncols];
  int *rcols = new int[ncols];

  for (int i = 0; i < nrows; i++) {
    int nr = 0;  // Number of entries in the current row

    // Add the matrix elements to the current row of the matrix.
    // These new elements are sorted.
    int diag_flag = 0;
    for (int j = mat->data->rowp[i]; j < mat->data->rowp[i + 1]; j++) {
      if (mat->data->cols[j] == i) {
        diag_flag = 1;
      }
      rcols[nr] = mat->data->cols[j];
      rlevs[nr] = 0;
      nr++;
    }

    // No diagonal element associated with row i, add one!
    if (!diag_flag) {
      nr = TacsMergeSortedArrays(nr, rcols, 1, &i);
    }

    // Now, perform the symbolic factorization -- this generates new entries
    int j = 0;
    for (; rcols[j] < i; j++) {  // For entries in this row, before the diagonal
      int clev = rlevs[j];       // the level of fill for this entry

      int p = j + 1;                   // The index into rcols
      int k_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

      // Start with the first entry after the diagonal in row, cols[j]
      // k is the index into cols for row cols[j]
      for (int k = diag[rcols[j]] + 1; k < k_end; k++) {
        // Increment p to an entry where we may have cols[k] == rcols[p]
        while (p < nr && rcols[p] < cols[k]) {
          p++;
        }

        // The element already exists, check if it has a lower level of fill
        // and update the fill level if necessary
        if (p < nr && rcols[p] == cols[k]) {
          if (rlevs[p] > (clev + levs[k] + 1)) {
            rlevs[p] = clev + levs[k] + 1;
          }
        } else if ((clev + levs[k] + 1) <= levFill) {
          // The element does not exist but should since the level of
          // fill is low enough. Insert the new entry into the list,
          // but keep the list sorted
          for (int n = nr; n > p; n--) {
            rlevs[n] = rlevs[n - 1];
            rcols[n] = rcols[n - 1];
          }

          rlevs[p] = clev + levs[k] + 1;
          rcols[p] = cols[k];
          nr++;
        }
      }
    }

    // Check if the size will be exceeded by adding the new elements
    if (size + nr > max_size) {
      int mat_ext = (int)((fill - 1.0) * mat_size);
      if (nr > mat_ext) {
        mat_ext = nr;
      }
      max_size = max_size + mat_ext;
      TacsExtendArray(&cols, size, max_size);
      TacsExtendArray(&levs, size, max_size);
    }

    // Now, put the new entries into the cols/levs arrays
    for (int k = 0; k < nr; k++) {
      cols[size] = rcols[k];
      levs[size] = rlevs[k];
      size++;
    }

    rowp[i + 1] = size;
    diag[i] = j + rowp[i];
  }

  // Clip the cols array to the correct size
  if (max_size > size) {
    TacsExtendArray(&cols, size, size);
  }

  if (mat->data->rowp[nrows] > 0) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    printf(
        "[%d] BCSRMat: ILU(%d) Input fill ratio %4.2f, actual "
        "fill ratio: %4.2f, nnz(ILU) = %d\n",
        rank, levFill, fill, (1.0 * rowp[nrows]) / mat->data->rowp[nrows],
        rowp[nrows]);
  }

  delete[] rcols;
  delete[] rlevs;

  // Store the rowp/cols and diag arrays
  data->rowp = rowp;
  data->cols = cols;
  data->diag = diag;

  *_levs = levs;
}

/*!
  Compute the ILU fill in for the off-diagonal matrix E.
  [ L_{B} U_{B} |  L_{B}^{-1} E ]
*/
BCSRMat *BCSRMat::computeILUkEpc(BCSRMat *Emat, const int *levs, int levFill,
                                 double fill, int **_elevs) {
  // Retrieve the matrix data
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  // Compute compute fill in for the off-diagonal matrices
  // First compute the extra fill-in for the E-matrix
  int *rlevs = new int[Emat->data->ncols];
  int *rcols = new int[Emat->data->ncols];

  // Set up data for keeping track of the column entries
  int size = 0;

  // Number of non-zeros in the original matrix
  int mat_size = Emat->data->rowp[Emat->data->nrows];
  int max_size = (int)(fill * mat_size);  // The maximum size - for now

  int *ecols = new int[max_size];
  int *elevs = new int[max_size];  // The level of fill of an entry
  int *erowp = new int[Emat->data->nrows + 1];
  erowp[0] = 0;

  for (int i = 0; i < Emat->data->nrows; i++) {
    // Copy the columns from E to the temporary arrays
    int nr = 0;
    for (int j = Emat->data->rowp[i]; j < Emat->data->rowp[i + 1]; j++) {
      rcols[nr] = Emat->data->cols[j];
      rlevs[nr] = 0;
      nr++;
    }

    // Compute the symbolic factorization, generating new entries when necessary
    // Scan from the first to the diagonal entry of B
    for (int jp = rowp[i]; jp < diag[i]; jp++) {
      int clev = levs[jp];

      int p = 0;  // Start from the first entry of row i in Emat
      int j = cols[jp];

      for (int kp = erowp[j]; kp < erowp[j + 1]; kp++) {
        while (p < nr && rcols[p] < ecols[kp]) {
          p++;
        }

        if (p < nr && rcols[p] == ecols[kp]) {
          if (rlevs[p] > (clev + elevs[kp] + 1)) {
            rlevs[p] = clev + elevs[kp] + 1;
          }
        } else if ((clev + elevs[kp] + 1) <= levFill) {
          // The element does not exist but should since the level of
          // fill is low enough. Insert the new entry into the list,
          // but keep the list sorted
          for (int n = nr; n > p; n--) {
            rlevs[n] = rlevs[n - 1];
            rcols[n] = rcols[n - 1];
          }

          rlevs[p] = clev + elevs[kp] + 1;
          rcols[p] = ecols[kp];
          nr++;
        }
      }
    }

    // Check if the size will be exceeded by adding the new elements
    if (size + nr > max_size) {
      int mat_ext = (int)((fill - 1.0) * mat_size);
      if (nr > mat_ext) {
        mat_ext = nr;
      }
      max_size = max_size + mat_ext;
      TacsExtendArray(&ecols, size, max_size);
      TacsExtendArray(&elevs, size, max_size);
    }

    // Now, put the new entries into the cols/levs arrays
    for (int j = 0; j < nr; j++) {
      ecols[size] = rcols[j];
      elevs[size] = rlevs[j];
      size++;
    }

    erowp[i + 1] = size;
  }

  if (_elevs) {
    *_elevs = elevs;
  } else {
    delete[] elevs;
  }
  delete[] rlevs;
  delete[] rcols;

  // Clip the ecols array to the correct length
  if (max_size > size) {
    TacsExtendArray(&ecols, size, size);
  }

  return new BCSRMat(comm, thread_info, Emat->data->bsize, Emat->data->nrows,
                     Emat->data->ncols, &erowp, &ecols);
}

BCSRMat *BCSRMat::computeILUkFpc(BCSRMat *Fmat, const int *levs, int levFill,
                                 double fill, int **_flevs) {
  // Retrieve the matrix data
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int *diag = data->diag;

  // Compute compute fill in for the off-diagonal matrices
  // First compute the extra fill-in for the E-matrix
  int *rlevs = new int[Fmat->data->ncols];
  int *rcols = new int[Fmat->data->ncols];

  // Set up data for keeping track of the column entries
  int size = 0;

  // Number of non-zeros in the original matrix
  int mat_size = Fmat->data->rowp[Fmat->data->nrows];
  int max_size = (int)(fill * mat_size);  // The maximum size - for now

  int *fcols = new int[max_size];
  int *flevs = new int[max_size];  // The level of fill of an entry
  int *frowp = new int[Fmat->data->nrows + 1];
  frowp[0] = 0;

  // Compute extra fill-in for the F-matrix
  for (int i = 0; i < Fmat->data->nrows; i++) {
    int nr = 0;

    // Copy the columns from F to the temporary array
    for (int j = Fmat->data->rowp[i]; j < Fmat->data->rowp[i + 1]; j++) {
      rcols[nr] = Fmat->data->cols[j];
      rlevs[nr] = 0;
      nr++;
    }

    // Scan through row i -- this is row N+i of the global matrix
    int j = 0;
    for (; j < nr; j++) {   // Scan through all row elements
      int clev = rlevs[j];  // the level of fill for this entry

      int p = j + 1;                   // The index into rcols
      int k_end = rowp[rcols[j] + 1];  // the end of row number cols[j]

      // Start with the first entry after the diagonal in row, cols[j]
      // k is the index into cols for row cols[j]
      for (int k = diag[rcols[j]] + 1; k < k_end; k++) {
        // Increment p to an entry where we may have cols[k] == rcols[p]
        while (p < nr && rcols[p] < cols[k]) {
          p++;
        }

        // The element already exists, check if it has a lower level of fill
        // and update the fill level if necessary
        if (p < nr && rcols[p] == cols[k]) {
          if (rlevs[p] > (clev + levs[k] + 1)) {
            rlevs[p] = clev + levs[k] + 1;
          }
        } else if ((clev + levs[k] + 1) <= levFill) {
          // The element does not exist but should since the level of
          // fill is low enough. Insert the new entry into the list,
          // but keep the list sorted
          for (int n = nr; n > p; n--) {
            rlevs[n] = rlevs[n - 1];
            rcols[n] = rcols[n - 1];
          }

          rlevs[p] = clev + levs[k] + 1;
          rcols[p] = cols[k];
          nr++;
        }
      }
    }

    // Check if the size will be exceeded by adding the new elements
    if (size + nr > max_size) {
      int mat_ext = (int)((fill - 1.0) * mat_size);
      if (nr > mat_ext) {
        mat_ext = nr;
      }
      max_size = max_size + mat_ext;
      TacsExtendArray(&fcols, size, max_size);
      TacsExtendArray(&flevs, size, max_size);
    }

    // Now, put the new entries into the cols/levs arrays
    for (int j = 0; j < nr; j++) {
      fcols[size] = rcols[j];
      flevs[size] = rlevs[j];
      size++;
    }

    frowp[i + 1] = size;
  }

  if (_flevs) {
    *_flevs = flevs;
  } else {
    delete[] flevs;
  }
  delete[] rlevs;
  delete[] rcols;

  // Clip the fcols array to the correct length
  if (max_size > size) {
    TacsExtendArray(&fcols, size, size);
  }

  return new BCSRMat(comm, thread_info, Fmat->data->bsize, Fmat->data->nrows,
                     Fmat->data->ncols, &frowp, &fcols);
}

/*
  Compute the location of the diagonal entry for each row
*/
void BCSRMat::setUpDiag() {
  if (!data->diag) {
    data->diag = new int[data->nrows];
  }

  for (int i = 0; i < data->nrows; i++) {
    int row_size = data->rowp[i + 1] - data->rowp[i];

    // Figure out the location of the diagonal entry
    int *item = TacsSearchArray(i, row_size, &data->cols[data->rowp[i]]);
    if (item == NULL) {
      data->diag[i] = -1;  // No diagonal entry
    } else {
      data->diag[i] = item - data->cols;
    }
  }
}

/*
  Initialize the generic implementation of each of the low-level
  routines
*/
void BCSRMat::initGenericImpl() {
  // The serial versions
  bfactor = BCSRMatFactor;
  applylower = BCSRMatApplyLower;
  applyupper = BCSRMatApplyUpper;
  bmult = BCSRMatVecMult;
  bmultadd = BCSRMatVecMultAdd;
  bmulttrans = BCSRMatVecMultTranspose;
  bmatmult = BCSRMatMatMultAdd;
  bfactorlower = BCSRMatFactorLower;
  bfactorupper = BCSRMatFactorUpper;
  bmatmatmultnormal = BCSRMatMatMultNormal;
  applypartiallower = BCSRMatApplyPartialLower;
  applypartialupper = BCSRMatApplyPartialUpper;
  applyschur = BCSRMatApplyFactorSchur;
  applysor = BCSRMatApplySOR;

  // No default threaded versions
  bmultadd_thread = NULL;
  bfactor_thread = NULL;
  applylower_thread = NULL;
  applyupper_thread = NULL;
  bmatmult_thread = NULL;
  bfactorlower_thread = NULL;
  bfactorupper_thread = NULL;
}

/*
  Initialize the block-specific implementations of each low-level
  routine
*/
void BCSRMat::initBlockImpl() {
  initGenericImpl();

  data->matvec_group_size = 32;
  data->matmat_group_size = 4;

  switch (data->bsize) {
    case 1:
      bfactor = BCSRMatFactor1;
      applylower = BCSRMatApplyLower1;
      applyupper = BCSRMatApplyUpper1;
      bmult = BCSRMatVecMult1;
      bmultadd = BCSRMatVecMultAdd1;
      bmulttrans = BCSRMatVecMultTranspose1;
      bmatmult = BCSRMatMatMultAdd1;
      bfactorlower = BCSRMatFactorLower1;
      bfactorupper = BCSRMatFactorUpper1;
      applypartiallower = BCSRMatApplyPartialLower1;
      applypartialupper = BCSRMatApplyPartialUpper1;
      applyschur = BCSRMatApplyFactorSchur1;
      bmatmatmultnormal = BCSRMatMatMultNormal1;
      applysor = BCSRMatApplySOR1;
      break;
    case 2:
      bfactor = BCSRMatFactor2;
      applylower = BCSRMatApplyLower2;
      applyupper = BCSRMatApplyUpper2;
      bmult = BCSRMatVecMult2;
      bmultadd = BCSRMatVecMultAdd2;
      bmatmult = BCSRMatMatMultAdd2;
      bfactorlower = BCSRMatFactorLower2;
      bfactorupper = BCSRMatFactorUpper2;
      applypartiallower = BCSRMatApplyPartialLower2;
      applypartialupper = BCSRMatApplyPartialUpper2;
      applyschur = BCSRMatApplyFactorSchur2;
      applysor = BCSRMatApplySOR2;
      break;
    case 3:
      bfactor = BCSRMatFactor3;
      applylower = BCSRMatApplyLower3;
      applyupper = BCSRMatApplyUpper3;
      bmult = BCSRMatVecMult3;
      bmultadd = BCSRMatVecMultAdd3;
      bmatmult = BCSRMatMatMultAdd3;
      bfactorlower = BCSRMatFactorLower3;
      bfactorupper = BCSRMatFactorUpper3;
      applypartiallower = BCSRMatApplyPartialLower3;
      applypartialupper = BCSRMatApplyPartialUpper3;
      applyschur = BCSRMatApplyFactorSchur3;
      applysor = BCSRMatApplySOR3;
      break;
    case 4:
      bfactor = BCSRMatFactor4;
      applylower = BCSRMatApplyLower4;
      applyupper = BCSRMatApplyUpper4;
      bmult = BCSRMatVecMult4;
      bmultadd = BCSRMatVecMultAdd4;
      bmatmult = BCSRMatMatMultAdd4;
      bfactorlower = BCSRMatFactorLower4;
      bfactorupper = BCSRMatFactorUpper4;
      applypartiallower = BCSRMatApplyPartialLower4;
      applypartialupper = BCSRMatApplyPartialUpper4;
      applyschur = BCSRMatApplyFactorSchur4;
      applysor = BCSRMatApplySOR4;
      break;
    case 5:
      bfactor = BCSRMatFactor5;
      applylower = BCSRMatApplyLower5;
      applyupper = BCSRMatApplyUpper5;
      bmult = BCSRMatVecMult5;
      bmultadd = BCSRMatVecMultAdd5;
      bmatmult = BCSRMatMatMultAdd5;
      bfactorlower = BCSRMatFactorLower5;
      bfactorupper = BCSRMatFactorUpper5;
      applypartiallower = BCSRMatApplyPartialLower5;
      applypartialupper = BCSRMatApplyPartialUpper5;
      applyschur = BCSRMatApplyFactorSchur5;
      applysor = BCSRMatApplySOR5;
      break;
    case 6:
      // These are tuning parameters
      data->matvec_group_size = 16;
      data->matmat_group_size = 4;

      // The serial versions
      bfactor = BCSRMatFactor6;
      applylower = BCSRMatApplyLower6;
      applyupper = BCSRMatApplyUpper6;
      bmult = BCSRMatVecMult6;
      bmultadd = BCSRMatVecMultAdd6;
      bmatmult = BCSRMatMatMultAdd6;
      bfactorlower = BCSRMatFactorLower6;
      bfactorupper = BCSRMatFactorUpper6;
      applypartiallower = BCSRMatApplyPartialLower6;
      applypartialupper = BCSRMatApplyPartialUpper6;
      applyschur = BCSRMatApplyFactorSchur6;
      applysor = BCSRMatApplySOR6;

      // The threaded versions
      bmultadd_thread = BCSRMatVecMultAdd6_thread;
      bfactor_thread = BCSRMatFactor6_thread;
      applylower_thread = BCSRMatApplyLower6_thread;
      applyupper_thread = BCSRMatApplyUpper6_thread;
      bmatmult_thread = BCSRMatMatMultAdd6_thread;
      bfactorlower_thread = BCSRMatFactorLower6_thread;
      bfactorupper_thread = BCSRMatFactorUpper6_thread;
      break;
    case 8:
      // These are tuning parameters
      data->matvec_group_size = 16;
      data->matmat_group_size = 4;

      // The serial versions
      bfactor = BCSRMatFactor8;
      applylower = BCSRMatApplyLower8;
      applyupper = BCSRMatApplyUpper8;
      bmult = BCSRMatVecMult8;
      bmultadd = BCSRMatVecMultAdd8;
      bmatmult = BCSRMatMatMultAdd8;
      bfactorlower = BCSRMatFactorLower8;
      bfactorupper = BCSRMatFactorUpper8;
      applypartiallower = BCSRMatApplyPartialLower8;
      applypartialupper = BCSRMatApplyPartialUpper8;
      applyschur = BCSRMatApplyFactorSchur8;
      applysor = BCSRMatApplySOR8;

      // The threaded versions
      bmultadd_thread = BCSRMatVecMultAdd8_thread;
      bfactor_thread = BCSRMatFactor8_thread;
      applylower_thread = BCSRMatApplyLower8_thread;
      applyupper_thread = BCSRMatApplyUpper8_thread;
      bmatmult_thread = BCSRMatMatMultAdd8_thread;
      bfactorlower_thread = BCSRMatFactorLower8_thread;
      bfactorupper_thread = BCSRMatFactorUpper8_thread;
      break;
    default:
      break;
  }
}

// Functions related to solving the system of equations
// ----------------------------------------------------

/*!
  Perform an ILU factorization of the matrix using the existing
  non-zero pattern. The entries are over-written, all operations are
  performed in place.
*/
void BCSRMat::factor() {
  if (!data->diag) {
    setUpDiag();
  }

  if (bfactor_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    tdata->init_apply_lower_sched();

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bfactor_thread, (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }
  } else {
    bfactor(data);
  }
}

/*!
  Compute y = A*x
*/
void BCSRMat::mult(TacsScalar *xvec, TacsScalar *yvec) {
  if (bmultadd_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    tdata->init_mat_mult_sched();
    tdata->input = xvec;
    tdata->output = yvec;
    memset(yvec, 0, data->bsize * data->nrows * sizeof(TacsScalar));

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bmultadd_thread, (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }

    tdata->input = tdata->output = NULL;
  } else {
    bmult(data, xvec, yvec);
  }
}

/*!
  Compute y = A*x + z
*/
void BCSRMat::multAdd(TacsScalar *xvec, TacsScalar *zvec, TacsScalar *yvec) {
  if (bmultadd_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    // Assign the input/output info
    tdata->init_mat_mult_sched();
    tdata->input = xvec;
    tdata->output = yvec;

    if (zvec != yvec) {
      memcpy(yvec, zvec, data->bsize * data->nrows * sizeof(TacsScalar));
    }

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bmultadd_thread, (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }

    tdata->input = tdata->output = NULL;
  } else {
    bmultadd(data, xvec, yvec, zvec);
  }
}

/*!
  Compute y = A^{T}*x
*/
void BCSRMat::multTranspose(TacsScalar *xvec, TacsScalar *yvec) {
  memset(yvec, 0, data->bsize * data->ncols * sizeof(TacsScalar));
  bmulttrans(data, xvec, yvec);
}

/*!
  Apply the ILU factorization of the matrix to the input vector

  Apply (LU)^{-1} x = y -- only for a factored matrix

  y = U^{-1} L^{-1} x
*/
void BCSRMat::applyFactor(TacsScalar *xvec, TacsScalar *yvec) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyFactor error: matrix not factored\n");
  } else {
    if (applylower_thread && applyupper_thread &&
        thread_info->getNumThreads() > 1) {
      // If not allocated, allocate the threaded data
      if (!tdata) {
        tdata = new BCSRMatThread(data);
        tdata->incref();
      }

      if (yvec != xvec) {
        memcpy(yvec, xvec, data->bsize * data->nrows * sizeof(TacsScalar));
      }
      tdata->output = yvec;

      // Create the joinable attribute
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

      // Apply L^{-1}
      tdata->init_apply_lower_sched();

      // Go through and run the threads
      for (long k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_create(&tdata->threads[k], &attr, applylower_thread,
                       (void *)tdata);
      }

      for (int k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_join(tdata->threads[k], NULL);
      }

      // Apply U^{-1}
      tdata->init_apply_upper_sched();

      // Go through and run the threads
      for (long k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_create(&tdata->threads[k], &attr, applyupper_thread,
                       (void *)tdata);
      }

      for (int k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_join(tdata->threads[k], NULL);
      }

      // Destroy the attribute and join all the threads
      pthread_attr_destroy(&attr);
    } else {
      applylower(data, xvec, yvec);
      applyupper(data, yvec, yvec);
    }
  }
}

/*!
  Apply the ILU factorization of the matrix to the input vector

  Apply (LU)^{-1} x = x -- only for a factored matrix

  x = U^{-1} L^{-1} x
*/
void BCSRMat::applyFactor(TacsScalar *xvec) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyFactor error: matrix not factored\n");
  } else {
    if (applylower_thread && applyupper_thread &&
        thread_info->getNumThreads() > 1) {
      // If not allocated, allocate the threaded data
      if (!tdata) {
        tdata = new BCSRMatThread(data);
        tdata->incref();
      }

      tdata->output = xvec;

      // Create the joinable attribute
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

      // Apply L^{-1}
      tdata->init_apply_lower_sched();

      // Go through and run the threads
      for (long k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_create(&tdata->threads[k], &attr, applylower_thread,
                       (void *)tdata);
      }

      for (int k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_join(tdata->threads[k], NULL);
      }

      // Apply U^{-1}
      tdata->init_apply_upper_sched();

      // Go through and run the threads
      for (long k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_create(&tdata->threads[k], &attr, applyupper_thread,
                       (void *)tdata);
      }

      for (int k = 0; k < thread_info->getNumThreads(); k++) {
        pthread_join(tdata->threads[k], NULL);
      }

      // Destroy the attribute and join all the threads
      pthread_attr_destroy(&attr);
    } else {
      applylower(data, xvec, xvec);
      applyupper(data, xvec, xvec);
    }
  }
}

/*!
  Apply only the upper portion of the ILU factorization

  y = U^{-1} x
*/
void BCSRMat::applyUpper(TacsScalar *xvec, TacsScalar *yvec) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyUpper error: matrix not factored\n");
  } else {
    applyupper(data, xvec, yvec);
  }
}

/*!
  Apply only the upper portion of the ILU factorization

  y = L^{-1} x
*/
void BCSRMat::applyLower(TacsScalar *xvec, TacsScalar *yvec) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyLower error: matrix not factored\n");
  } else {
    applylower(data, xvec, yvec);
  }
}

/*!
  Apply only a part of L^{-1} to the input vector.

  Take the vector input. Assume that it is only of length length nrows
  - var_offset ie. x = x_full[var_offset:nrows] and its block
  structure is set up the same way as the full vector.  Apply the part
  of L^{-1} from var_offset to nrows, with ordering starting at zero.
*/
void BCSRMat::applyPartialLower(TacsScalar *xvec, int var_offset) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyPartialLower error: matrix not factored\n");
  } else {
    applypartiallower(data, xvec, var_offset);
  }
}

/*!
  Apply only a part of U^{-1} to the input vector.

  Take the vector input. Assume that it is only of length nrows -
  var_offset ie. x = x_full[var_offset:nrows] and its block structure
  is set up the same way as the full vector.  Apply the part of U^{-1}
  from var_offset to nrows, with ordering starting at zero.
*/
void BCSRMat::applyPartialUpper(TacsScalar *xvec, int var_offset) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyPartialUpper error: matrix not factored\n");
  } else {
    applypartialupper(data, xvec, var_offset);
  }
}

/*!
  Special function for the approximate Schur preconditioner.

  Given the input vector x = f, and y

  Compute x = U_b^{-1} ( L_b^{-1} f - (L_b^{-1} E) y )

  The matrix is factorized into the following form:
  A     = [ L_b          0   ][ U_b  L_b^{-1} E ]
  .       [ F U_b^{-1}   L_s ][ 0    U_s        ]

  where the division is set by the variable var_offset.
*/
void BCSRMat::applyFactorSchur(TacsScalar *x, int var_offset) {
  if (!data->diag) {
    fprintf(stderr, "BCSRMat applyFactorSchur error: matrix not factored\n");
  } else {
    applyschur(data, x, var_offset);
  }
}

/*!
  Copy the diagonal entries to a set of diagonal matrices.  Factor
  these matrices and store the result.
*/
void BCSRMat::factorDiag(const TacsScalar *diag) {
  if (!data->diag) {
    setUpDiag();
  }

  const int bsize = data->bsize;
  const int b2 = bsize * bsize;
  const int nrows = data->nrows;

  int *ipiv = new int[2 * bsize];
  TacsScalar *D = new TacsScalar[4 * b2];

  if (!Adiag) {
    Adiag = new TacsScalar[nrows * b2];
  }

  for (int i = 0; i < nrows; i++) {
    int d = b2 * data->diag[i];
    memcpy(D, &(data->A[d]), b2 * sizeof(TacsScalar));
    if (diag) {
      for (int j = 0; j < bsize; j++) {
        D[j * (bsize + 1)] += diag[bsize * i + j];
      }
    }
    BMatComputeInverse(&Adiag[b2 * i], D, ipiv, bsize);
  }

  delete[] D;
  delete[] ipiv;
}

/*!
  Apply a given number of steps of SOR to the system A*x = b.
*/
void BCSRMat::applySOR(TacsScalar *b, TacsScalar *x, TacsScalar omega,
                       int iters) {
  if (Adiag) {
    for (int i = 0; i < iters; i++) {
      applysor(data, NULL, 0, data->nrows, 0, Adiag, omega, b, NULL, x);
    }
  } else {
    fprintf(stderr, "Cannot apply SOR: diagonal has not been factored\n");
  }
}

/*!
  Apply SOR to the matrix over the given interval
*/
void BCSRMat::applySOR(BCSRMat *B, int start, int end, int var_offset,
                       TacsScalar omega, const TacsScalar *b,
                       const TacsScalar *xext, TacsScalar *x) {
  if (Adiag) {
    if (B) {
      applysor(data, B->data, start, end, var_offset, Adiag, omega, b, xext, x);
    } else {
      applysor(data, NULL, start, end, var_offset, Adiag, omega, b, xext, x);
    }
  } else {
    fprintf(stderr, "Cannot apply SOR: diagonal has not been factored\n");
  }
}

/*!
  Compute the matrix-matrix product.

  This matrix must be created with a call to BCSRMat( A, B, fill ) for this
  matrix to have the correct non-zero pattern.
*/
void BCSRMat::matMultAdd(double alpha, BCSRMat *amat, BCSRMat *bmat) {
  // Check that the sizes work
  if (data->bsize != amat->data->bsize || data->bsize != bmat->data->bsize) {
    fprintf(stderr,
            "BCSRMat error: cannot multiply matrices with "
            "different block sizes\n");
  }
  if (data->nrows != amat->data->nrows ||
      amat->data->ncols != bmat->data->nrows ||
      bmat->data->ncols != data->ncols) {
    fprintf(stderr,
            "BCSRMat error: cannot multiply matrices, "
            "incorrect dimensions:\n");
    fprintf(stderr, " dim(this) = %d,%d\n dim(A) = %d,%d\n dim(B) = %d,%d\n",
            data->nrows, data->ncols, amat->data->nrows, amat->data->ncols,
            bmat->data->nrows, bmat->data->ncols);
  }

  if (bmatmult_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    // Set the data for the matrix multiplication
    tdata->init_mat_mult_sched();
    tdata->Amat = amat->data;
    tdata->Bmat = bmat->data;
    tdata->alpha = alpha;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bmatmult_thread, (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }

    tdata->alpha = 0.0;
    tdata->Amat = tdata->Bmat = NULL;
  } else {
    bmatmult(alpha, amat->data, bmat->data, data);
  }
}

/*!
  Compute L^{-1} E
*/
void BCSRMat::applyLowerFactor(BCSRMat *emat) {
  if (!data->diag) {
    fprintf(stderr,
            "BCSRMat error: cannot use applyLowerFactor "
            "with an un-factored matrix\n");
    return;
  }
  if (data->nrows != emat->data->nrows) {
    fprintf(stderr,
            "BCSRMat error: matrices are not the "
            "correction dimensions\n");
  }

  if (bfactorlower_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    tdata->Amat = emat->data;
    tdata->init_apply_lower_sched();

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bfactorlower_thread,
                     (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }

    tdata->Amat = NULL;
  } else {
    bfactorlower(data, emat->data);
  }
}

/*!
  Compute F U^{-1}
*/
void BCSRMat::applyUpperFactor(BCSRMat *fmat) {
  if (!data->diag) {
    fprintf(stderr,
            "BCSRMat error: cannot use applyUpperFactor with "
            "an un-factored matrix\n");
    return;
  }
  if (data->nrows != fmat->data->ncols) {
    fprintf(stderr,
            "BCSRMat error: matrices are not the "
            "correction dimensions\n");
  }

  if (bfactorupper_thread && thread_info->getNumThreads() > 1) {
    // If not allocated, allocate the threaded data
    if (!tdata) {
      tdata = new BCSRMatThread(data);
      tdata->incref();
    }

    tdata->Amat = fmat->data;
    tdata->init_mat_mult_sched();

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Go through and run the threads
    for (long k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&tdata->threads[k], &attr, bfactorupper_thread,
                     (void *)tdata);
    }

    // Destroy the attribute and join all the threads
    pthread_attr_destroy(&attr);
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(tdata->threads[k], NULL);
    }

    tdata->Amat = NULL;
  } else {
    bfactorupper(data, fmat->data);
  }
}

/*
  Compute A = B^{T}*S*B where S = diag{s}

  Note that there is not yet a threaded implementation of
  this code.
*/
void BCSRMat::matMultNormal(TacsScalar *s, BCSRMat *bmat) {
  if (data->nrows != bmat->data->ncols) {
    fprintf(stderr,
            "BCSRMat error: matMultNormal matrices are not the "
            "correction dimensions\n");
  }

  zeroEntries();
  bmatmatmultnormal(data, s, bmat->data);
}

/*!
  Zero all entries of the matrix
*/
void BCSRMat::zeroEntries() {
  int bsize = data->bsize;
  int length = data->rowp[data->nrows];
  length *= bsize * bsize;

  memset(data->A, 0, length * sizeof(TacsScalar));
}

/*!
  Scale all the entries in the matrix by a factor
*/
void BCSRMat::scale(TacsScalar alpha) {
  const int bsize = data->bsize;
  int length = data->rowp[data->nrows];
  length *= bsize * bsize;

  int one = 1;
  BLASscal(&length, &alpha, data->A, &one);
}

/*!
  Add values into the matrix row.  The values may only be added
  into parts of the matrix that have existing non-zero pattern. Trying
  to insert values outside these entries will generate an error
  message.

  Input:
  row:       the value of the row index
  ncol:      the number of columns to add
  col:       the column indices (before being offset)
  nca:       the number of columns in a
  avals:     the values
*/
void BCSRMat::addRowValues(int row, int ncol, const int *col, int nca,
                           const TacsScalar *avals) {
  if (ncol <= 0) {
    return;
  }

  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  if (row >= 0 && row < nrows) {
    int row_size = rowp[row + 1] - rowp[row];
    const int *col_array = &cols[rowp[row]];

    for (int j = 0; j < ncol; j++) {
      int c = col[j];
      if (c < 0) {  // Skip entries that are negative
        continue;
      } else if (c < ncols) {
        int *item = TacsSearchArray(c, row_size, col_array);

        if (item == NULL) {
          fprintf(stderr, "BCSRMat error: no entry for (%d,%d)\n", row, c);
        } else {
          // Place the values into the array
          int cp = item - cols;
          int bj = bsize * j;
          TacsScalar *a = &(data->A[b2 * cp]);

          for (int jj = 0; jj < bsize; jj++) {
            int njj = nca * jj;
            int bjj = bsize * jj;
            for (int ii = 0; ii < bsize; ii++) {
              a[ii + bjj] += avals[ii + bj + njj];
            }
          }
        }
      } else {
        fprintf(stderr, "BCSRMat error: column %d out of range [0,%d)\n", c,
                ncols);
      }
    }
  } else {
    fprintf(stderr, "BCSRMat error: row %d out of range [0,%d)\n", row, nrows);
  }
}

/*!
  Add weighted values into the matrix row. Values can only be added to
  locations that have an existing non-zero pattern. Trying to add
  values to locations without a non-zero entry produces an error
  message. The entries in avals are multiplied both by the scalar
  alpha and the values stored in weights.

  NORMAL:

  | A11 | A12 | A13 |
  | A21 | A22 | A23 |
  | A31 | A32 | A33 |

  TRANSPOSE:

  | A11' | A21' | A31' |
  | A12' | A22' | A32' |
  | A13' | A23' | A33' |

  input:
  alpha    the value of the scalar multiplier
  row      the row index
  ncol     the number of columns
  col      the column indices
  weights  the number of weights; len(weights) = ncol
  nca      number of columns in the matrix avals
  avals    the values to add
  matOr:   the matrix orientation
*/
void BCSRMat::addRowWeightValues(TacsScalar alpha, int row, int nwrows,
                                 const int *wrowp, const int *wcols,
                                 const TacsScalar *weights, int nca,
                                 const TacsScalar *avals,
                                 MatrixOrientation matOr) {
  if (nwrows <= 0 || alpha == 0.0) {
    return;
  }

  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  if (row >= 0 && row < nrows) {
    int row_size = rowp[row + 1] - rowp[row];
    const int *col_array = &cols[rowp[row]];

    for (int i = 0; i < nwrows; i++) {
      for (int jp = wrowp[i]; jp < wrowp[i + 1]; jp++) {
        int c = wcols[jp];
        TacsScalar aw = alpha * weights[jp];

        // Skip entries that are negative or have a
        // zero weight value
        if (c < 0 || aw == 0.0) {
          continue;
        } else if (c < ncols) {
          int *item = TacsSearchArray(c, row_size, col_array);

          if (item == NULL) {
            fprintf(stderr, "BCSRMat error: no entry for (%d,%d)\n", row, c);
          } else {
            // Place the values into the array
            int cp = item - cols;
            TacsScalar *A = &(data->A[b2 * cp]);

            if (matOr == TACS_MAT_NORMAL) {
              // Set the offset in the dense input row to the current
              // block matrix entry added to the BCSRMat
              const int offset = bsize * i;

              // Loop over the rows of the input matrix
              for (int ii = 0; ii < bsize; ii++) {
                const TacsScalar *arow = &avals[nca * ii + offset];

                // Add the row entries to the block row
                TacsScalar *Arow = &A[ii * bsize];
                for (int jj = 0; jj < bsize; jj++) {
                  Arow[jj] += aw * arow[jj];
                }
              }
            } else {
              // Set the offset in the dense input row to the current
              // block matrix entry added to the BCSRMat
              const int offset = nca * bsize * i;

              // Loop over the rows of the input matrix
              for (int ii = 0; ii < bsize; ii++) {
                const TacsScalar *acol = &avals[ii + offset];

                // Add the row entries to the block row
                TacsScalar *Arow = &A[ii * bsize];
                for (int jj = 0; jj < bsize; jj++) {
                  Arow[jj] += aw * acol[nca * jj];
                }
              }
            }
          }
        } else {
          fprintf(stderr, "BCSRMat error: column %d out of range [0,%d)\n", c,
                  ncols);
        }
      }
    }
  } else {
    fprintf(stderr, "BCSRMat error: row %d out of range [0,%d)\n", row, nrows);
  }
}

/*!
  Add values into the matrix row.  The values may only be added into
  parts of the matrix that have existing non-zero pattern. Trying to
  insert values outside these entries will generate an error message.

  Input:
  row:        the value of the row index
  ncol:       the number of columns to add
  col:        the column indices (before being offset)
  avals:      the values
*/
void BCSRMat::addBlockRowValues(int row, int ncol, const int *col,
                                const TacsScalar *avals) {
  if (ncol <= 0) {
    return;
  }

  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  if (row >= 0 && row < nrows) {
    int row_size = rowp[row + 1] - rowp[row];
    const int *col_array = &cols[rowp[row]];

    for (int j = 0; j < ncol; j++) {
      int c = col[j];

      if (c < 0) {  // Skip entries that are negative
        continue;
      } else if (c < ncols) {
        int *item = TacsSearchArray(c, row_size, col_array);

        if (item == NULL) {
          fprintf(stderr, "BCSRMat error: no entry for (%d,%d)\n", row, c);
        } else {
          // Add values into the array
          int cp = item - cols;
          TacsScalar *a = &(data->A[b2 * cp]);

          for (int k = 0; k < b2; k++) {
            a[k] += avals[b2 * j + k];
          }
        }
      } else {
        fprintf(stderr,
                "BCSRMat error: column %d "
                "out of range [0,%d)\n",
                c, ncols);
      }
    }
  } else {
    fprintf(stderr, "BCSRMat error: row %d out of range [0,%d)\n", row, nrows);
  }
}

/*!
  Zero the row values. Possibly set the diagonal elements to unity.

  row:      the row of the matrix
  vars:     an integer containing binary variables
  ident:    flag to indicate whether to set the diagonal to 1
*/
void BCSRMat::zeroRow(int row, int vars, int ident) {
  if (row >= 0 && row < data->nrows) {
    const int *rowp = data->rowp;
    const int *cols = data->cols;
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;

    int end = rowp[row + 1];
    for (int j = rowp[row]; j < end; j++) {
      TacsScalar *a = &(data->A[b2 * j]);
      for (int ii = 0; ii < bsize; ii++) {
        if (vars & (1 << ii)) {
          for (int jj = 0; jj < bsize; jj++) {
            a[bsize * ii + jj] = 0.0;
          }
        }
      }
      if (ident && row == cols[j]) {
        for (int ii = 0; ii < bsize; ii++) {
          if (vars & (1 << ii)) {
            a[(bsize + 1) * ii] = 1.0;
          }
        }
      }
    }
  }
}

/*!
  Zero the column values
*/
void BCSRMat::zeroColumns(int num_zero_cols, const int *zero_cols,
                          const int *zero_vars, int ident) {
  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;
  const int *diag = data->diag;

  // Count up the number of columns
  int *colp = new int[ncols + 1];
  memset(colp, 0, (ncols + 1) * sizeof(int));

  // Set the entries of colp to -1 and only those entries
  // with the correct index = -1
  for (int i = 0; i < ncols; i++) {
    colp[i] = -1;
  }
  for (int i = 0; i < num_zero_cols; i++) {
    int col = zero_cols[i];
    if (col >= 0 && col < ncols) {
      colp[col] = 0;
    }
  }

  for (int jp = 0; jp < rowp[nrows]; jp++) {
    int col = cols[jp];
    if (colp[col] >= 0) {
      colp[col]++;
    }
  }

  // Adjust the indices
  int index = 0;
  for (int i = 0; i < ncols; i++) {
    int temp = colp[i];
    colp[i] = index;
    if (temp >= 0) {
      index += temp;
    }
  }
  colp[ncols] = index;

  // Allocate space for the rows
  int *ptr = new int[colp[ncols]];

  for (int jp = 0; jp < rowp[nrows]; jp++) {
    int col = cols[jp];
    if (colp[col + 1] - colp[col] > 0) {
      ptr[colp[col]] = jp;
      colp[col]++;
    }
  }

  // Reset the columns
  for (int i = ncols; i > 0; i--) {
    colp[i] = colp[i - 1];
  }
  colp[0] = 0;

  // Now, zero the columns
  for (int i = 0; i < num_zero_cols; i++) {
    int col = zero_cols[i];
    if (col >= 0 && col < ncols) {
      for (int j = colp[col]; j < colp[col + 1]; j++) {
        int jp = ptr[j];

        TacsScalar *a = &(data->A[b2 * jp]);
        for (int jj = 0; jj < bsize; jj++) {
          if (zero_vars[i] & (1 << jj)) {
            for (int ii = 0; ii < bsize; ii++) {
              a[bsize * ii + jj] = 0.0;
            }
          }
        }
      }
    }
  }

  if (ident && diag) {
    for (int i = 0; i < num_zero_cols; i++) {
      int col = zero_cols[i];
      if (col >= 0 && col < ncols) {
        int jp = diag[col];

        TacsScalar *a = &(data->A[b2 * jp]);
        for (int jj = 0; jj < bsize; jj++) {
          if (zero_vars[i] & (1 << jj)) {
            for (int ii = 0; ii < bsize; ii++) {
              a[bsize * ii + jj] = 0.0;
            }
          }
        }
      }
    }
  }

  // Free the allocated space
  delete[] colp;
  delete[] ptr;
}

/*!
  Partition the existing matrix into four sub-matrices.
  This routine should be avoided if possible however,
  it's useful for testing purposes.

  This partitions the matrix into:

  A = [ B, E ]
  .   [ F, C ]
*/
void BCSRMat::partition(int nrows_p, BCSRMat **Bmat, BCSRMat **Emat,
                        BCSRMat **Fmat, BCSRMat **Cmat) {
  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int *diag = data->diag;
  const TacsScalar *A = data->A;

  int b2 = bsize * bsize;
  *Bmat = NULL;
  *Emat = NULL;
  *Fmat = NULL;
  *Cmat = NULL;

  int nrows_b = nrows_p;
  int ncols_b = nrows_p;
  int nrows_c = nrows - nrows_p;
  int ncols_c = ncols - nrows_p;

  // Create the B and E non-zero pattern
  int *browp = new int[nrows_b + 1];
  int *erowp = new int[nrows_b + 1];

  int nb = 0, ne = 0;
  browp[0] = 0;
  erowp[0] = 0;

  for (int i = 0; i < nrows_b; i++) {
    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] < nrows_b) {
        nb++;
      } else {
        ne++;
      }
    }
    browp[i + 1] = nb;
    erowp[i + 1] = ne;
  }

  int *bcols = new int[nb];
  int *ecols = new int[ne];

  nb = ne = 0;
  for (int i = 0; i < nrows_b; i++) {
    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] < nrows_b) {
        bcols[nb] = cols[j];
        nb++;
      } else {
        ecols[ne] = cols[j] - nrows_b;
        ne++;
      }
    }
  }

  *Bmat =
      new BCSRMat(comm, thread_info, bsize, nrows_b, ncols_b, &browp, &bcols);
  *Emat =
      new BCSRMat(comm, thread_info, bsize, nrows_b, ncols_c, &erowp, &ecols);
  if (diag) {
    (*Bmat)->setUpDiag();
  }

  // Copy over the values to the matrix
  {
    const int *browp, *bcols;
    const int *erowp, *ecols;
    int bs;
    TacsScalar *B, *E;
    (*Bmat)->getArrays(&bs, &nrows_b, &ncols_b, &browp, &bcols, &B);
    (*Emat)->getArrays(&bs, &nrows_b, &ncols_c, &erowp, &ecols, &E);

    nb = 0, ne = 0;
    for (int i = 0; i < nrows_b; i++) {
      for (int j = rowp[i]; j < rowp[i + 1]; j++) {
        if (cols[j] < nrows_b) {
          for (int k = 0; k < b2; k++, nb++) {
            B[nb] = A[b2 * j + k];
          }
        } else {
          for (int k = 0; k < b2; k++, ne++) {
            E[ne] = A[b2 * j + k];
          }
        }
      }
    }
  }

  // Create the B and E non-zero pattern
  int *frowp = new int[nrows_b + 1];
  int *crowp = new int[nrows_b + 1];

  int nf = 0, nc = 0;
  frowp[0] = 0;
  crowp[0] = 0;

  for (int i = nrows_b; i < nrows; i++) {
    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] < nrows_b) {
        nf++;
      } else {
        nc++;
      }
    }
    frowp[i + 1 - nrows_b] = nf;
    crowp[i + 1 - nrows_b] = nc;
  }

  int *fcols = new int[nf];
  int *ccols = new int[nc];

  nf = nc = 0;
  for (int i = nrows_b; i < nrows; i++) {
    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] < nrows_b) {
        fcols[nf] = cols[j];
        nf++;
      } else {
        ccols[nc] = cols[j] - nrows_b;
        nc++;
      }
    }
  }

  *Fmat =
      new BCSRMat(comm, thread_info, bsize, nrows_c, ncols_b, &frowp, &fcols);
  *Cmat =
      new BCSRMat(comm, thread_info, bsize, nrows_c, ncols_c, &crowp, &ccols);
  if (diag) {
    (*Cmat)->setUpDiag();
  }

  // Copy over the values to the matrix
  {
    const int *frowp, *fcols;
    const int *crowp, *ccols;
    int bs;
    TacsScalar *F, *C;
    (*Fmat)->getArrays(&bs, &nrows_c, &ncols_b, &frowp, &fcols, &F);
    (*Cmat)->getArrays(&bs, &nrows_c, &ncols_c, &crowp, &ccols, &C);

    nf = 0, nc = 0;
    for (int i = nrows_b; i < nrows; i++) {
      for (int j = rowp[i]; j < rowp[i + 1]; j++) {
        if (cols[j] < nrows_b) {
          for (int k = 0; k < b2; k++, nf++) {
            F[nf] = A[b2 * j + k];
          }
        } else {
          for (int k = 0; k < b2; k++, nc++) {
            C[nc] = A[b2 * j + k];
          }
        }
      }
    }
  }
}

/*!
  Retrieve the underlying array representation of the BCSRMatrix
*/
void BCSRMat::getArrays(int *_bsize, int *_nrows, int *_ncols,
                        const int **_rowp, const int **_cols,
                        TacsScalar **Avals) {
  if (_bsize) {
    *_bsize = data->bsize;
  }
  if (_nrows) {
    *_nrows = data->nrows;
  }
  if (_ncols) {
    *_ncols = data->ncols;
  }
  if (_rowp) {
    *_rowp = data->rowp;
  }
  if (_cols) {
    *_cols = data->cols;
  }
  if (Avals) {
    *Avals = data->A;
  }
}

/*
  Get the matrix in a dense column-major format appropriate for LAPACK
*/
void BCSRMat::getDenseColumnMajor(TacsScalar *D) {
  const int bsize = data->bsize;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const TacsScalar *A = data->A;

  // Zero the array first and fill the non-zeros alone
  memset(D, 0, bsize * bsize * nrows * nrows * sizeof(TacsScalar));

  // The leading dimension of the dense column matrix
  const int ldd = bsize * nrows;

  // Loop over the block rows
  for (int ib = 0; ib < nrows; ib++) {
    // Loop over the block columns
    for (int jp = rowp[ib]; jp < rowp[ib + 1]; jp++) {
      int jb = cols[jp];
      const TacsScalar *a = &A[bsize * bsize * jp];

      // Now iterate over the block indices
      for (int ii = 0; ii < bsize; ii++) {
        int i = ii + bsize * ib;
        for (int jj = 0; jj < bsize; jj++) {
          int j = jj + bsize * jb;
          D[i + ldd * j] = a[jj + bsize * ii];
        }
      }
    }
  }
}

/*!
  Copy values from the given matrix into this matrix.

  Scan through each row of the matrix, copying entries.
*/
void BCSRMat::copyValues(BCSRMat *mat) {
  if (mat->data->nrows != data->nrows || mat->data->ncols != data->ncols ||
      data->bsize != mat->data->bsize) {
    fprintf(stderr,
            "BCSRMat: Matrices are not the same size "
            "cannot copy values\n");
    return;
  }

  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  this->zeroEntries();

  for (int i = 0; i < nrows; i++) {
    int p = rowp[i];
    int end = rowp[i + 1];
    int mat_end = mat->data->rowp[i + 1];

    for (int j = mat->data->rowp[i]; (j < mat_end) && (p < end); j++) {
      while (cols[p] < mat->data->cols[j] && p < end) {
        p++;
      }

      // Copy the matrix if the two entries are equal
      // and the size of both block--matrices is the same
      if (p < end) {
        if (cols[p] == mat->data->cols[j]) {
          int n = b2 * p;
          int m = b2 * j;

          int last = b2 * (p + 1);
          for (; n < last; n++, m++) {
            data->A[n] = mat->data->A[m];
          }
        } else {
          fprintf(stderr, "BCSRMat: NZ-pattern error cannot copy values\n");
        }
      }
    }
  }
}

/*!
  Add the entries of one matrix times a scalar value into this matrix.

  Compute y <- y + alpha * x

  The non-zero pattern of the two matrices need not match exactly -
  the matrix being copied into may be larger than the original matrix
*/

void BCSRMat::axpy(TacsScalar alpha, BCSRMat *mat) {
  if (mat->data->nrows != data->nrows || mat->data->ncols != data->ncols ||
      mat->data->bsize != data->bsize) {
    fprintf(stderr,
            "BCSRMat: Matrices are not the same "
            "size cannot apply axpy\n");
    return;
  }

  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  for (int i = 0; i < nrows; i++) {
    int p = rowp[i];
    int end = rowp[i + 1];
    int mat_end = mat->data->rowp[i + 1];

    for (int j = mat->data->rowp[i]; (j < mat_end) && (p < end); j++) {
      while (cols[p] < mat->data->cols[j] && p < end) {
        p++;
      }

      // Copy the matrix if the two entries are equal
      // and the size of both block--matrices is the same
      if (p < end) {
        if (cols[p] == mat->data->cols[j]) {
          int n = b2 * p;
          int m = b2 * j;

          int last = b2 * (p + 1);
          for (; n < last; n++, m++) {
            data->A[n] += alpha * mat->data->A[m];
          }
        } else {
          fprintf(stderr, "BCSRMat: NZ-pattern error cannot apply axpy\n");
        }
      }
    }
  }
}

/*!
  Scale the values in this matrix and then add the entries of
  one matrix times a scalar value into this one.

  Compute y <- alpha * x + beta * y

  The non-zero pattern of the two matrices need not be exactly the
  same.  The y matrix (this matrix) - can have more non-zero entries,
  that are simply scaled by beta since there are no corresponding
  entries in x. However, if there are more entries in x, no new space
  will be allocated in y for these additional entries. The non-zero
  patterns are static.
*/
void BCSRMat::axpby(TacsScalar alpha, TacsScalar beta, BCSRMat *mat) {
  if (mat->data->nrows != data->nrows || mat->data->ncols != data->ncols) {
    fprintf(stderr,
            "BCSRMat: Matrices are not the same "
            "size cannot apply axpby\n");
  }

  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;
  const int b2 = bsize * bsize;

  for (int i = 0; i < nrows; i++) {
    int p = rowp[i];
    int end = rowp[i + 1];
    int mat_end = mat->data->rowp[i + 1];

    for (int j = mat->data->rowp[i]; (j < mat_end) && (p < end); j++) {
      while (cols[p] < mat->data->cols[j] && p < end) {
        int last = b2 * cols[p + 1];
        for (int n = b2 * cols[p]; n < last; n++) {
          data->A[n] = beta * data->A[n];
        }

        p++;
      }

      // Copy the matrix if the two entries are equal
      // and the size of both block--matrices is the same
      if (p < end) {
        if (cols[p] == mat->data->cols[j]) {
          int n = b2 * p;
          int m = b2 * j;

          int last = b2 * (p + 1);
          for (; n < last; n++, m++) {
            data->A[n] = beta * data->A[n] + alpha * mat->data->A[m];
          }
        } else {
          fprintf(stderr, "BCSRMat: NZ-pattern error cannot apply axpy\n");
        }
      }
    }
  }
}

/*!
  Add a value to the diagonal entries of the matrix
*/
void BCSRMat::addDiag(TacsScalar alpha) {
  if (data->diag) {
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;
    for (int i = 0; i < data->nrows; i++) {
      TacsScalar *a = &(data->A[b2 * data->diag[i]]);
      for (int k = 0; k < bsize; k++) {
        a[(bsize + 1) * k] += alpha;
      }
    }
  } else {
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;
    for (int i = 0; i < data->nrows; i++) {
      for (int jp = data->rowp[i]; jp < data->rowp[i + 1]; jp++) {
        int j = data->cols[jp];
        if (i == j) {
          TacsScalar *a = &(data->A[b2 * jp]);
          for (int k = 0; k < bsize; k++) {
            a[(bsize + 1) * k] += alpha;
          }
        }
      }
    }
  }
}

/*!
  Add an array of values to the diagonal entries of the matrix
*/
void BCSRMat::addDiag(TacsScalar *alpha) {
  if (data->diag) {
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;
    for (int i = 0; i < data->nrows; i++) {
      TacsScalar *d = &alpha[bsize * i];
      TacsScalar *a = &(data->A[b2 * data->diag[i]]);
      for (int k = 0; k < bsize; k++) {
        a[(bsize + 1) * k] += d[k];
      }
    }
  } else {
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;
    for (int i = 0; i < data->nrows; i++) {
      TacsScalar *d = &alpha[bsize * i];
      for (int jp = data->rowp[i]; jp < data->rowp[i + 1]; jp++) {
        int j = data->cols[jp];
        if (i == j) {
          TacsScalar *a = &(data->A[b2 * jp]);
          for (int k = 0; k < bsize; k++) {
            a[(bsize + 1) * k] += d[k];
          }
        }
      }
    }
  }
}

/*!
  Compute the bandwidth of the matrix.

  Bandwidth = max b_{i}

  where b_{i} = max |cols[j] - i|
*/
void BCSRMat::getNumUpperLowerDiagonals(int *_bl, int *_bu) {
  int bl = 0;
  int bu = 0;

  for (int i = 0; i < data->nrows; i++) {
    if (data->rowp[i + 1] - data->rowp[i] > 0) {
      int ibl = i - data->cols[data->rowp[i]] - 1;
      int ibu = data->cols[data->rowp[i + 1] - 1] - i - 1;

      if (ibl > bl) {
        bl = ibl;
      }
      if (ibu > bu) {
        bu = ibu;
      }
    }
  }

  *_bl = bl * data->bsize;
  *_bu = bu * data->bsize;
}

/*
  Extract the stored matrix in a banded matrix storage format
  compatible with LAPACK.

  if symm_flag = False, supply the full matrix:

  matrix entry i,j is stored in A[bu+i-j+(bu+bl+1)*j]

  if symm_flag = True, supply the upper matrix only:

  matrix entry i,j is stored in A[bu+i-j+(bu+1)*j]

  input:
  size:      the total size of the array
  symm_flag: the symmetry flag - use a symmetric banded form

  in/out:
  A:         the array that will contain the required entries

  returns: the matrix band size
*/
void BCSRMat::getBandedMatrix(TacsScalar *A, int size, int symm_flag) {
  // Compute the matrix bandwidth
  int bl, bu;
  getNumUpperLowerDiagonals(&bl, &bu);

  if (symm_flag) {
    // Check if the size of the array A is sufficiently large
    int mat_size = data->nrows * data->bsize * (bu + 1);

    if (mat_size >= size) {
      memset(A, 0, mat_size * sizeof(TacsScalar));

      // Set the block size
      int bsize = data->bsize;

      // Scan through the entries in the matrix
      for (int i = 0; i < data->nrows; i++) {
        for (int jp = data->rowp[i]; jp < data->rowp[i + 1]; jp++) {
          int j = data->cols[jp];
          if (j > i) {
            // Extract the block into the matrix
            for (int ii = 0; ii < bsize; ii++) {
              for (int jj = 0; jj < bsize; jj++) {
                int index = bu + (i * bsize + ii) + bu * (j * bsize + jj);
                A[index] = data->A[bsize * bsize * jp + ii + bsize * jj];
              }
            }
          } else if (j == i) {
            // Extract the block into the matrix
            for (int ii = 0; ii < bsize; ii++) {
              for (int jj = ii; jj < bsize; jj++) {
                int index = bu + (i * bsize + ii) + bu * (j * bsize + jj);
                A[index] = data->A[bsize * bsize * jp + ii + bsize * jj];
              }
            }
          }
        }
      }
    }
  } else {
    // Check if the size of the array A is sufficiently large
    int mat_size = data->nrows * data->bsize * (bl + bu + 1);

    if (mat_size >= size) {
      memset(A, 0, mat_size * sizeof(TacsScalar));

      // Set the block size
      int bsize = data->bsize;

      // Scan through the entries in the matrix
      for (int i = 0; i < data->nrows; i++) {
        for (int jp = data->rowp[i]; jp < data->rowp[i + 1]; jp++) {
          int j = data->cols[jp];

          // Extract the block into the matrix
          for (int ii = 0; ii < bsize; ii++) {
            for (int jj = 0; jj < bsize; jj++) {
              int index = bu + (i * bsize + ii) + (bu + bl) * (j * bsize + jj);
              A[index] = data->A[bsize * bsize * jp + ii + bsize * jj];
            }
          }
        }
      }
    }
  }
}

/*!
  Test if the matrices are equal to a prescirbed value.

  1. Test the non-zero pattern to see if they are the same.
  2. Test matrix-multiplication against a randomly generated vector
*/
int BCSRMat::isEqual(BCSRMat *mat, double tol) {
  if (mat->data->nrows != data->nrows || mat->data->ncols != data->ncols) {
    printf("Matrices do not have the same dimensions\n");
    return 0;
  }

  const int ncols = data->ncols;
  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;
  const int bsize = data->bsize;

  for (int i = 0; i < nrows; i++) {
    if (rowp[i + 1] != mat->data->rowp[i + 1]) {
      printf("Row pointers do not match for entry %3d\n", i + 1);

      printf("A row %d \n", i);
      for (int j = rowp[i]; j < rowp[i + 1]; j++) {
        printf("%d ", cols[j]);
      }
      printf("\nB row %d\n", i);
      for (int j = mat->data->rowp[i]; j < mat->data->rowp[i + 1]; j++) {
        printf("%d ", mat->data->cols[j]);
      }
      printf("\n");

      return 0;
    }

    for (int j = rowp[i]; j < rowp[i + 1]; j++) {
      if (cols[j] != mat->data->cols[j]) {
        printf("Column arrays do not match in entry (%3d,%3d)\n ", i, cols[j]);
        return 0;
      }
    }
  }

  printf("Non-zero pattern matches\n");

  srand(time(NULL));

  TacsScalar *a = new TacsScalar[bsize * nrows];
  TacsScalar *b = new TacsScalar[bsize * nrows];
  TacsScalar *x = new TacsScalar[bsize * ncols];

  for (int i = 0; i < bsize * ncols; i++) {
    x[i] = (2.0 * rand()) / (RAND_MAX - 1.0) - 1.0;
  }

  mat->mult(x, a);
  mult(x, b);

  TacsScalar norm = 0.0;
  for (int i = 0; i < bsize * nrows; i++) {
    norm += (a[i] - b[i]) * (a[i] - b[i]);
  }

  delete[] a;
  delete[] b;
  delete[] x;

  norm = sqrt(norm);
  if (TacsRealPart(norm) < tol) {
    printf("Matrices are essentially equal |A x - B x| = %15.5e \n",
           TacsRealPart(norm));
    return 1;
  }

  printf("Matrices do not match |A x - B x| = %15.5e \n", TacsRealPart(norm));

  return 0;
}

/*!
  Test the implementation of the Schur-complement approaches

  - Take the current matrix A and split it,
  A = [ B, E ]
  .   [ F, C ]

  - Compute the approximate factorization of A and split that according
  Apc = [ Bpc, Epc ]
  .     [ Fpc, Cpc ]

  - Form the Schur complement with the following approach
  1. Construct Bsc, Esc, Fsc from the matrix factorizations
  2. Compute Csc = C - Fsc * Bsc^{-1} Esc

  - Compare to results computed with the full matrix factorization
  technique
*/
void BCSRMat::testSchur(int nrows_p, int lev, double fill, double tol) {
  BCSRMat *B, *E, *F, *C;
  partition(nrows_p, &B, &E, &F, &C);
  B->incref();
  E->incref();
  F->incref();
  C->incref();

  printf("dim(B) = (%d,%d) \n", B->data->nrows, B->data->ncols);
  printf("dim(E) = (%d,%d) \n", E->data->nrows, E->data->ncols);
  printf("dim(F) = (%d,%d) \n", F->data->nrows, F->data->ncols);
  printf("dim(C) = (%d,%d) \n", C->data->nrows, C->data->ncols);

  BCSRMat *Apc = new BCSRMat(comm, this, lev, fill);
  Apc->incref();
  Apc->copyValues(this);
  Apc->factor();

  BCSRMat *Bpc, *Epc, *Fpc, *Cpc;
  Apc->partition(nrows_p, &Bpc, &Epc, &Fpc, &Cpc);
  Bpc->incref();
  Epc->incref();
  Fpc->incref();
  Cpc->incref();

  BCSRMat *Bsc, *Esc, *Fsc;
  BCSRMat *S;
  int use_full_schur = 1;
  Bsc =
      new BCSRMat(comm, B, E, F, C, lev, fill, &Esc, &Fsc, &S, use_full_schur);
  Bsc->incref();
  Esc->incref();
  Fsc->incref();
  S->incref();

  Bsc->copyValues(B);
  Bsc->factor();

  Esc->copyValues(E);
  Fsc->copyValues(F);

  Bsc->applyLowerFactor(Esc);
  Bsc->applyUpperFactor(Fsc);

  printf("Comparing Epc and Esc\n");
  Epc->isEqual(Esc, tol);

  printf("Comparing Fpc and Fsc\n");
  Fpc->isEqual(Fsc, tol);

  // Test the Schur complement itself
  S->copyValues(C);
  S->matMultAdd(-1.0, Fsc, Esc);

  printf("Comparing Schur complement matrices\n");
  S->factor();
  S->isEqual(Cpc);

  S->decref();
  Bsc->decref();
  Esc->decref();
  Fsc->decref();
  Bpc->decref();
  Epc->decref();
  Fpc->decref();
  Cpc->decref();
  Apc->decref();
  B->decref();
  E->decref();
  F->decref();
  C->decref();
}

/*!
  Print the matrix to a file that can be read in later.

  The format of the file is as follows:

  nrows, ncols, bsize // (matrix size/block sizes)
  rowp                // the pointer to each new matrix row
  cols                // the column indices
  A                   // the matrix entries
*/
void BCSRMat::printMat(const char *fname) {
  // Print the matrix in a somewhat human readable format
  FILE *fp = fopen(fname, "w");

  if (fp) {
    const int nrows = data->nrows;
    const int ncols = data->ncols;
    const int *rowp = data->rowp;
    const int *cols = data->cols;
    const int bsize = data->bsize;
    const int b2 = bsize * bsize;

    // Print in the following format:
    // number or rows, number of columns, block size
    fprintf(fp, "%d %d %d\n", nrows, ncols, bsize);

    // Print out the row-pointer array
    for (int i = 0; i < nrows + 1; i++) {
      fprintf(fp, "%d\n", rowp[i]);
    }

    // Print out the cols array
    for (int jp = 0; jp < rowp[nrows]; jp++) {
      fprintf(fp, "%d\n", cols[jp]);
    }

    // Print out the block data
    for (int jp = 0; jp < rowp[nrows]; jp++) {
      for (int ii = 0; ii < b2; ii += bsize) {
        for (int jj = 0; jj < bsize; jj++) {
          fprintf(fp, "%.16e ", TacsRealPart(data->A[b2 * jp + jj + ii]));
        }
        fprintf(fp, "\n");
      }
    }

    fclose(fp);
  }
}

void BCSRMat::printNzPattern(const char *fname) {
  // Print the matrix in a somewhat human readable format
  FILE *fp = fopen(fname, "w");

  const int nrows = data->nrows;
  const int *rowp = data->rowp;
  const int *cols = data->cols;

  if (fp) {
    fprintf(fp, "VARIABLES = \"i\",\"j\"\n");
    for (int i = 0; i < nrows; i++) {
      for (int j = rowp[i]; j < rowp[i + 1]; j++) {
        fprintf(fp, "%d %d\n", i, cols[j]);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
  }
}
