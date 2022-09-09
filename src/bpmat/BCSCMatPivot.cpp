/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "BCSCMatPivot.h"

#include "TacsUtilities.h"
#include "tacslapack.h"

/*
  Matrix classes based on compressed-sparse-column formats with
  partial pivoting
*/

/*
  The following code is used within the matrix factorization
  routines. This is block-specific code that is designed to run faster
  than non-block-specific code.

  This computes the following:

  for ( int i = 0; i < n; i++ )
  .  for ( int j = 0; j < m; j++ )
  .     X[xdim*index[k] + j] += Y[ydim*k + j]

  input:
  xdim:  the dimension of the x-array
  Y:     the y column
  ydim:  the dimension of the y-array
  index: row indices into the destination column
  n:     the number of rows to transfer
  m:     the number of entries per row

  output:
  X:     the values are added to this array
*/
static inline void addScatterColumn(TacsScalar *X, const int xdim,
                                    const TacsScalar *Y, const int ydim,
                                    const int *index, const int n,
                                    const int m) {
  if (m == 1) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      Y += ydim;
      index++;
    }
  } else if (m == 2) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      Y += ydim;
      index++;
    }
  } else if (m == 3) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      Y += ydim;
      index++;
    }
  } else if (m == 4) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      Y += ydim;
      index++;
    }
  } else if (m == 5) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      Y += ydim;
      index++;
    }
  } else if (m == 6) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      Y += ydim;
      index++;
    }
  } else if (m == 7) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      Y += ydim;
      index++;
    }
  } else if (m == 8) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      Y += ydim;
      index++;
    }
  } else if (m == 9) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      Y += ydim;
      index++;
    }
  } else if (m == 10) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      Y += ydim;
      index++;
    }
  } else if (m == 11) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      Y += ydim;
      index++;
    }
  } else if (m == 12) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      x[11] += Y[11];
      Y += ydim;
      index++;
    }
  } else if (m == 13) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      x[11] += Y[11];
      x[12] += Y[12];
      Y += ydim;
      index++;
    }
  } else if (m == 14) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      x[11] += Y[11];
      x[12] += Y[12];
      x[13] += Y[13];
      Y += ydim;
      index++;
    }
  } else if (m == 15) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      x[11] += Y[11];
      x[12] += Y[12];
      x[13] += Y[13];
      x[14] += Y[14];
      Y += ydim;
      index++;
    }
  } else if (m == 16) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] += Y[0];
      x[1] += Y[1];
      x[2] += Y[2];
      x[3] += Y[3];
      x[4] += Y[4];
      x[5] += Y[5];
      x[6] += Y[6];
      x[7] += Y[7];
      x[8] += Y[8];
      x[9] += Y[9];
      x[10] += Y[10];
      x[11] += Y[11];
      x[12] += Y[12];
      x[13] += Y[13];
      x[14] += Y[14];
      x[15] += Y[15];
      Y += ydim;
      index++;
    }
  } else {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      for (int j = 0; j < m; j++) {
        x[j] += Y[j];
      }
      Y += ydim;
      index++;
    }
  }
}

/*
  The following code is used within the matrix factorization
  routines. This is block-specific code that is designed to run faster
  than non-block-specific code.

  This computes the following:

  for ( int i = 0; i < n; i++ )
  .  for ( int j = 0; j < m; j++ )
  .     X[xsize*index[k] + j] -= Y[ydim*k + j]

  input:
  xdim:  the dimension of the x-array
  Y:     the y column
  ydim:  the dimension of the y-array
  index: row indices into the destination column
  n:     the number of rows to transfer
  m:     the number of entries per row

  output:
  X:     the values are added to this array
*/
static inline void subtractScatterColumn(TacsScalar *X, const int xdim,
                                         const TacsScalar *Y, const int ydim,
                                         const int *index, const int n,
                                         const int m) {
  if (m == 1) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      Y += ydim;
      index++;
    }
  } else if (m == 2) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      Y += ydim;
      index++;
    }
  } else if (m == 3) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      Y += ydim;
      index++;
    }
  } else if (m == 4) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      Y += ydim;
      index++;
    }
  } else if (m == 5) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      Y += ydim;
      index++;
    }
  } else if (m == 6) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      Y += ydim;
      index++;
    }
  } else if (m == 7) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      Y += ydim;
      index++;
    }
  } else if (m == 8) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      Y += ydim;
      index++;
    }
  } else if (m == 9) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      Y += ydim;
      index++;
    }
  } else if (m == 10) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      Y += ydim;
      index++;
    }
  } else if (m == 11) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      Y += ydim;
      index++;
    }
  } else if (m == 12) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      x[11] -= Y[11];
      Y += ydim;
      index++;
    }
  } else if (m == 13) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      x[11] -= Y[11];
      x[12] -= Y[12];
      Y += ydim;
      index++;
    }
  } else if (m == 14) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      x[11] -= Y[11];
      x[12] -= Y[12];
      x[13] -= Y[13];
      Y += ydim;
      index++;
    }
  } else if (m == 15) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      x[11] -= Y[11];
      x[12] -= Y[12];
      x[13] -= Y[13];
      x[14] -= Y[14];
      Y += ydim;
      index++;
    }
  } else if (m == 16) {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      x[0] -= Y[0];
      x[1] -= Y[1];
      x[2] -= Y[2];
      x[3] -= Y[3];
      x[4] -= Y[4];
      x[5] -= Y[5];
      x[6] -= Y[6];
      x[7] -= Y[7];
      x[8] -= Y[8];
      x[9] -= Y[9];
      x[10] -= Y[10];
      x[11] -= Y[11];
      x[12] -= Y[12];
      x[13] -= Y[13];
      x[14] -= Y[14];
      x[15] -= Y[15];
      Y += ydim;
      index++;
    }
  } else {
    for (int i = 0; i < n; i++) {
      TacsScalar *x = &X[xdim * index[0]];
      for (int j = 0; j < m; j++) {
        x[j] -= Y[j];
      }
      Y += ydim;
      index++;
    }
  }
}

/*
  Create the block-based compressed sparse column matrix for
  storage. This constructor takes in a local communicator and the
  non-zero pattern of the matrix. This call steals the ownership of
  the colp/rows arrays. You should not deallocate/modify/use the
  colp/rows/A arrays after you've called this function.

  input:
  comm:   the MPI communicator
  bsize:  the block size used for the matrix
  nrows:  the number of block rows in the matrix
  ncols:  the number of block columns in the matrix
  colp:   the column pointer into the block rows array
  rows:   the block row index for each block matrix
  A:      the (optional) values of the entries in the matrix
*/
BCSCMat::BCSCMat(MPI_Comm _comm, int bsize, int _nrows, int _ncols, int **_colp,
                 int **_rows, TacsScalar **_A) {
  comm = _comm;
  rows_sorted = 0;

  // Copy the data from the matrix
  nrows = _nrows;
  ncols = _ncols;
  nblock_cols = ncols / bsize;
  colp = *_colp;
  rows = *_rows;

  // Create the bptr array
  max_block_size = bsize;
  bptr = new int[nblock_cols + 1];
  for (int i = 0; i < nblock_cols + 1; i++) {
    bptr[i] = i * bsize;
  }

  // Create the aptr array
  aptr = new int[nblock_cols + 1];

  // Set the offsets into the blocked array
  aptr[0] = 0;
  for (int i = 0; i < nblock_cols; i++) {
    int cdim = bptr[i + 1] = bptr[i];
    aptr[i + 1] = aptr[i] + cdim * (colp[i + 1] - colp[i]);
  }

  // Set the input pointers to NULL
  _colp = _rows = NULL;

  // Set the buffer size in doubles
  temp_array_size = 64000;
  temp_array = NULL;

  if (_A) {
    A = *_A;
    *_A = NULL;
  } else {
    // Allocate space for the matrix data
    int mat_size = bsize * colp[ncols];
    A = new TacsScalar[mat_size];
    memset(A, 0, mat_size * sizeof(TacsScalar));
  }
}

/*
  Create a block-based compressed sparse column matrix for storage.
  This constructor takes in a block CSR-type format and converts the
  internal storage into a CSC-type format. Note that the input arrays
  are not modified.

  input:
  comm:            the MPI communicator for this object
  bsize:           the block size for the matrix
  num_block_rows:  the number of block rows
  num_block_cols:  the number of block columns
  block_rowp:      the pointer to the beginning of each row
  block_cols:      the index of each block column in the matrix
*/
BCSCMat::BCSCMat(MPI_Comm _comm, int bsize, int num_block_rows,
                 int num_block_cols, const int *block_rowp,
                 const int *block_cols) {
  comm = _comm;
  rows_sorted = 1;  // Yes, the row indices will be sorted
  nrows = bsize * num_block_rows;
  ncols = bsize * num_block_cols;
  nblock_cols = num_block_cols;

  // Create the bptr array
  max_block_size = bsize;
  bptr = new int[nblock_cols + 1];
  for (int i = 0; i < nblock_cols + 1; i++) {
    bptr[i] = i * bsize;
  }

  // Create the CSC format from CSR
  initCSCfromCSR(num_block_rows, num_block_cols, block_rowp, block_cols);

  // Create the aptr array
  aptr = new int[nblock_cols + 1];

  // Set the offsets into the blocked array
  aptr[0] = 0;
  for (int i = 0; i < nblock_cols; i++) {
    int cdim = bptr[i + 1] - bptr[i];
    aptr[i + 1] = aptr[i] + cdim * (colp[i + 1] - colp[i]);
  }

  // Set the buffer size in doubles
  temp_array_size = 64000;
  temp_array = NULL;

  // Allocate space for the matrix data
  int mat_size = aptr[nblock_cols];
  A = new TacsScalar[mat_size];
  memset(A, 0, mat_size * sizeof(TacsScalar));
}

/*
  Create a block-based compressed sparse column matrix for storage.

  This constructor takes in a CSR-type data format with block sizes
  and converts the the internal storage into a CSC-type format.  Note
  that the input arrays are not modified.

  input:
  comm:            the MPI communicator for this object
  bsize:           the input column block sizes
  num_block_rows:  the number of block rows
  num_block_cols:  the number of block columns
  block_rowp:      the pointer to the beginning of each row
  block_cols:      the index of each block column in the matrix
*/
BCSCMat::BCSCMat(MPI_Comm _comm, const int bsize[], int num_block_rows,
                 int num_block_cols, const int *block_rowp,
                 const int *block_cols) {
  comm = _comm;
  rows_sorted = 1;  // The row indices will be sorted

  // Create the bptr array. Note that if the
  max_block_size = 0;
  int bptr_size = num_block_cols;
  if (num_block_rows > bptr_size) {
    bptr_size = num_block_rows;
  }

  bptr = new int[bptr_size + 1];
  bptr[0] = 0;
  for (int i = 0; i < bptr_size; i++) {
    bptr[i + 1] = bptr[i] + bsize[i];
    if (bsize[i] > max_block_size) {
      max_block_size = bsize[i];
    }
  }

  // Find the blocked row sizes
  nblock_cols = num_block_cols;
  nrows = bptr[num_block_rows];
  ncols = bptr[num_block_cols];

  // Create the CSC format from CSR
  initCSCfromCSR(num_block_rows, num_block_cols, block_rowp, block_cols);

  // Create the aptr array
  aptr = new int[nblock_cols + 1];

  // Set the offsets into the blocked array
  aptr[0] = 0;
  for (int i = 0; i < nblock_cols; i++) {
    int cdim = bptr[i + 1] - bptr[i];
    aptr[i + 1] = aptr[i] + cdim * (colp[i + 1] - colp[i]);
  }

  // Set the buffer size in doubles
  temp_array_size = 64000;
  temp_array = NULL;

  // Allocate space for the matrix data
  int mat_size = aptr[nblock_cols];
  A = new TacsScalar[mat_size];
  memset(A, 0, mat_size * sizeof(TacsScalar));
}

/*
  Create a block-based compressed sparse column matrix for storage.

  This constructor takes in a CSR-type data format with block sizes
  and converts the the internal storage into a CSC-type format.  Note
  that the input arrays are not modified.

  input:
  comm:            the MPI communicator for this object
  bsize:           the input column block sizes
  num_block_rows:  the number of block rows
  num_block_cols:  the number of block columns
  block_colp:      the pointer to the beginning of each column
  block_rows:      the index of each block index in the matrix
  rows_sorted:     flag to indicate whether the row indices are sorted
*/
BCSCMat::BCSCMat(MPI_Comm _comm, const int bsize[], int num_block_rows,
                 int num_block_cols, const int *block_colp,
                 const int *block_rows, int _rows_sorted) {
  comm = _comm;
  rows_sorted = _rows_sorted;

  // Create the bptr array. Note that if the
  max_block_size = 0;
  int bptr_size = num_block_cols;
  if (num_block_rows > bptr_size) {
    bptr_size = num_block_rows;
  }

  bptr = new int[bptr_size + 1];
  bptr[0] = 0;
  for (int i = 0; i < bptr_size; i++) {
    bptr[i + 1] = bptr[i] + bsize[i];
    if (bsize[i] > max_block_size) {
      max_block_size = bsize[i];
    }
  }

  // Find the blocked row sizes
  nblock_cols = num_block_cols;
  nrows = bptr[num_block_rows];
  ncols = bptr[num_block_cols];

  // Exapand each row to find the dimension of the array
  nblock_cols = num_block_cols;

  // Allocate space for the new arrays
  colp = new int[num_block_cols + 1];
  memset(colp, 0, (num_block_cols + 1) * sizeof(int));

  // Count up the full block column pointer array, expanding
  // the count for each block
  colp[0] = 0;
  for (int j = 0; j < num_block_cols; j++) {
    colp[j + 1] = colp[j];
    for (int ip = block_colp[j]; ip < block_colp[j + 1]; ip++) {
      int i = block_rows[ip];
      int dim = bptr[i + 1] - bptr[i];
      colp[j + 1] += dim;
    }
  }

  // Allocate space for the row indices
  rows = new int[colp[num_block_cols]];

  // Set the 'rows' array by scanning through the inupt CSR data
  for (int j = 0; j < num_block_cols; j++) {
    int irow = colp[j];
    for (int ip = block_colp[j]; ip < block_colp[j + 1]; ip++) {
      int i = block_rows[ip];
      int rdim = bptr[i + 1] - bptr[i];
      for (int ii = 0; ii < rdim; ii++, irow++) {
        rows[irow] = bptr[i] + ii;
      }
    }
  }

  // Create the aptr array
  aptr = new int[nblock_cols + 1];

  // Set the offsets into the blocked array
  aptr[0] = 0;
  for (int i = 0; i < nblock_cols; i++) {
    int cdim = bptr[i + 1] - bptr[i];
    aptr[i + 1] = aptr[i] + cdim * (colp[i + 1] - colp[i]);
  }

  // Set the buffer size in doubles
  temp_array_size = 64000;
  temp_array = NULL;

  // Allocate space for the matrix data
  int mat_size = aptr[nblock_cols];
  A = new TacsScalar[mat_size];
  memset(A, 0, mat_size * sizeof(TacsScalar));
}

/*
  Create the matrix from an input binary file

  This file should be generated by calling writeToFile() on another
  BCSCMat matrix object. This can be used to archive the matrix for
  later testing purposes.

  input:
  comm:      the communicator
  filename:  the binary file
*/
BCSCMat::BCSCMat(MPI_Comm _comm, const char *filename) {
  comm = _comm;
  FILE *fp = fopen(filename, "rb");

  // NULL out everything before we start to read things in
  bptr = NULL;
  aptr = NULL;
  colp = NULL;
  rows = NULL;
  A = NULL;

  // Write out all that data that is needed to reconstruct the matrix
  int values[5];
  if (fread(values, sizeof(int), 5, fp) != 5) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
    return;
  }
  rows_sorted = values[0];
  max_block_size = values[1];
  nrows = values[2];
  ncols = values[3];
  nblock_cols = values[4];

  // Next read in the pointer arrays
  bptr = new int[nblock_cols + 1];
  aptr = new int[nblock_cols + 1];
  colp = new int[nblock_cols + 1];
  unsigned int size = nblock_cols + 1;
  if (fread(bptr, sizeof(int), nblock_cols + 1, fp) != size) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
    return;
  }
  if (fread(aptr, sizeof(int), nblock_cols + 1, fp) != size) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
    return;
  }
  if (fread(colp, sizeof(int), nblock_cols + 1, fp) != size) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
    return;
  }

  // Create the rows array and read it in
  size = colp[nblock_cols];
  rows = new int[colp[nblock_cols]];
  if (fread(rows, sizeof(int), colp[nblock_cols], fp) != size) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
  }

  size = aptr[nblock_cols];
  A = new TacsScalar[aptr[nblock_cols]];
  if (fread(A, sizeof(TacsScalar), aptr[nblock_cols], fp) != size) {
    fprintf(stderr, "BCSCMat: Error reading binary file. Incorrect format?\n");
  }

  fclose(fp);

  // Set the buffer size in doubles
  temp_array_size = 64000;
  temp_array = NULL;
}

/*
  Copy the values and non-zero pattern of the input matrix into the
  current matrix.

  input:
  mat: the matrix to be copied into *this matrix
*/
BCSCMat::BCSCMat(BCSCMat *mat) {
  comm = mat->comm;
  rows_sorted = mat->rows_sorted;
  max_block_size = mat->max_block_size;
  nrows = mat->nrows;
  ncols = mat->ncols;
  nblock_cols = mat->nblock_cols;
  bptr = new int[nblock_cols + 1];
  aptr = new int[nblock_cols + 1];
  colp = new int[nblock_cols + 1];
  rows = new int[mat->colp[nblock_cols]];
  A = new TacsScalar[mat->aptr[nblock_cols]];

  // Set the buffer size
  temp_array_size = mat->temp_array_size;
  temp_array = NULL;

  memcpy(bptr, mat->bptr, (nblock_cols + 1) * sizeof(int));
  memcpy(aptr, mat->aptr, (nblock_cols + 1) * sizeof(int));
  memcpy(colp, mat->colp, (nblock_cols + 1) * sizeof(int));
  memcpy(rows, mat->rows, colp[nblock_cols] * sizeof(int));
  memcpy(A, mat->A, aptr[nblock_cols] * sizeof(TacsScalar));
}

/*
  Deallocate the data stored in the BCSC matrix
*/
BCSCMat::~BCSCMat() {
  if (aptr) {
    delete[] aptr;
  }
  if (bptr) {
    delete[] bptr;
  }
  if (colp) {
    delete[] colp;
  }
  if (rows) {
    delete[] rows;
  }
  if (A) {
    delete[] A;
  }
  if (temp_array) {
    delete[] temp_array;
  }
}

/*
  Write the matrix to a binary file that can be read in at a later
  time.

  This is useful for saving/debugging a matrix.

  input:
  filename:  the name of the file to be created
*/
void BCSCMat::writeToFile(const char *filename) {
  FILE *fp = fopen(filename, "wb");

  // Write out all that data that is needed to reconstruct the matrix
  fwrite(&rows_sorted, sizeof(int), 1, fp);
  fwrite(&max_block_size, sizeof(int), 1, fp);
  fwrite(&nrows, sizeof(int), 1, fp);
  fwrite(&ncols, sizeof(int), 1, fp);
  fwrite(&nblock_cols, sizeof(int), 1, fp);
  fwrite(bptr, sizeof(int), nblock_cols + 1, fp);
  fwrite(aptr, sizeof(int), nblock_cols + 1, fp);
  fwrite(colp, sizeof(int), nblock_cols + 1, fp);
  fwrite(rows, sizeof(int), colp[nblock_cols], fp);
  fwrite(A, sizeof(TacsScalar), aptr[nblock_cols], fp);

  fclose(fp);
}

/*
  Create the internal CSC data structures from the given (row and
  column) blocked CSR data

  This code assumes that the bptr object has already been allocated
  and filled in accordingly. Note that in this case, the bptr array
  must be for both rows/columns regardless of whether the matrix is
  square or rectangular.

  input:
  num_block_rows:  the number of block rows in the matrix
  num_block_cols:  the number of block columns
  block_rowp:      CSR-type pointer into each block row of the matrix
  block_cols:      the block column indicies
*/
void BCSCMat::initCSCfromCSR(int num_block_rows, int num_block_cols,
                             const int *block_rowp, const int *block_cols) {
  // Set the number of blocked columns in the matrix
  nblock_cols = num_block_cols;

  // Allocate space for the new arrays
  colp = new int[num_block_cols + 1];
  memset(colp, 0, (num_block_cols + 1) * sizeof(int));

  // Count up the number of references to each column
  colp[0] = 0;
  for (int i = 0; i < num_block_rows; i++) {
    int dim = bptr[i + 1] - bptr[i];
    for (int jp = block_rowp[i]; jp < block_rowp[i + 1]; jp++) {
      int j = block_cols[jp];
      colp[j + 1] += dim;
    }
  }

  // Set the colp array so that it points to the first entry of each column
  for (int i = 1; i <= num_block_cols; i++) {
    colp[i] += colp[i - 1];
  }

  // Allocate space for the row indices
  rows = new int[colp[num_block_cols]];

  // Set the 'rows' array by scanning through the inupt CSR data
  for (int i = 0; i < num_block_rows; i++) {
    int rdim = bptr[i + 1] - bptr[i];
    for (int j = block_rowp[i]; j < block_rowp[i + 1]; j++) {
      int block_column = block_cols[j];
      for (int ii = 0; ii < rdim; ii++) {
        rows[colp[block_column]] = bptr[i] + ii;
        colp[block_column]++;
      }
    }
  }

  // Reset the column pointer array
  int count = 0;
  for (int i = 0; i < num_block_cols; i++) {
    int next = colp[i];
    colp[i] = count;
    count = next;
  }
  colp[num_block_cols] = count;
}

/*
  Get the MPI communicator
*/
MPI_Comm BCSCMat::getMPIComm() { return comm; }

/*
  Zero all the entries in the matrix
*/
void BCSCMat::zeroEntries() {
  memset(A, 0, aptr[nblock_cols] * sizeof(TacsScalar));
}

/*
  Add the values into the matrix corresponding to the blocked entries
  in the matrix.

  The values are added to the specified block matrix entry. Note that
  the input matrix must be the right block size and this is not
  checked in advance.

  For efficiency, you should sort the rows before calling this
  function repeatedly or construct the matrix using the block
  compressed sparse row constructor since this sorts the row indices
  in each column automatically.

  input:
  num_rows:     the number of blocked rows
  brows:        the block matrix row
  num_cols:     the number of blocked columns
  bcols:        the block matrix column
  mat:          the row matrix stored in row-major order
  lda:          the leading row-dimension of the matrix
  is_tranpose:  are we adding the transpose of the block?
*/
void BCSCMat::addMatBlockValues(int num_rows, const int brows[], int num_cols,
                                const int bcols[], const TacsScalar *mat,
                                int ldmat, int is_transpose) {
  // There is probably a more inteligent way to do this, but this
  // code just breaks everything down by case: If adding the transpose
  // and by whether the rows are sorted or not.
  if (is_transpose) {
    if (rows_sorted) {
      // Add values using a bisection search for each column
      for (int i = 0, ioffset = 0; i < num_rows; i++) {
        // Get the block column for this portion of the matrix
        int col = brows[i];
        if (col < 0) {
          continue;
        }
        TacsScalar *column = &A[aptr[col]];

        // Determine the number of columns in the block
        int cdim = bptr[col + 1] - bptr[col];

        // Determine the number of non-zero block column entries in the
        // matrix to search them
        int col_size = colp[col + 1] - colp[col];
        const int *col_rows = &rows[colp[col]];

        // Loop over all the block rows in this column
        for (int j = 0, joffset = 0; j < num_cols; j++) {
          // Find the block row index
          int row = bcols[j];
          if (row < 0) {
            continue;
          }

          // Find the block row dimension
          int rdim = bptr[row + 1] - bptr[row];

          for (int ii = 0; ii < rdim; ii++) {
            int r = bptr[row] + ii;

            // Search the column to obtain the row index
            int *row_pos = TacsSearchArray(r, col_size, col_rows);
            if (row_pos) {
              int loc = row_pos - col_rows;
              TacsScalar *a = &column[cdim * loc];

              // Add the values into the matrix
              for (int jj = 0; jj < cdim; jj++) {
                a[jj] += mat[ldmat * (ioffset + jj) + (joffset + ii)];
              }
            } else {
              fprintf(stderr,
                      "BCSCMat: Entry (%d, %d) not in non-zero pattern\n", row,
                      col);
            }
          }

          // Increment the row pointer within the block matrix
          joffset += rdim;
        }

        // Increment the column pointer within the block matrix
        ioffset += cdim;
      }
    } else {
      // Add values using a bisection search for each column
      for (int i = 0, ioffset = 0; i < num_rows; i++) {
        // Get the block column for this portion of the matrix
        int col = brows[i];
        if (col < 0) {
          continue;
        }
        TacsScalar *column = &A[aptr[col]];

        // Determine the number of columns in the block
        int cdim = bptr[col + 1] - bptr[col];

        // Determine the number of non-zero block column entries in the
        // matrix to search them
        int col_size = colp[col + 1] - colp[col];
        const int *col_rows = &rows[colp[col]];

        // Loop over all the block rows in this column
        for (int j = 0, joffset = 0; j < num_cols; j++) {
          // Find the block row index
          int row = bcols[j];
          if (row < 0) {
            continue;
          }

          // Find the block row dimension
          int rdim = bptr[row + 1] - bptr[row];

          for (int ii = 0; ii < rdim; ii++) {
            int r = bptr[row] + ii;

            // Perform a (slow) linear search in the column
            const int *row_pos = col_rows;
            for (int k = 0; k < col_size; k++) {
              if (row_pos[0] == r) {
                break;
              }
              row_pos++;
            }

            if (row_pos[0] == r) {
              int loc = row_pos - col_rows;
              TacsScalar *a = &column[cdim * loc];

              // Add the values into the matrix
              for (int jj = 0; jj < cdim; jj++) {
                a[jj] += mat[ldmat * (ioffset + jj) + (joffset + ii)];
              }
            } else {
              fprintf(stderr,
                      "BCSCMat: Entry (%d, %d) not in non-zero pattern\n", row,
                      col);
            }
          }

          // Increment the row pointer within the block matrix
          joffset += rdim;
        }

        // Increment the column pointer within the block matrix
        ioffset += cdim;
      }
    }
  } else {
    if (rows_sorted) {
      // Add values using a bisection search for each column
      for (int j = 0, joffset = 0; j < num_cols; j++) {
        // Get the block column for this portion of the matrix
        int col = bcols[j];
        if (col < 0) {
          continue;
        }
        TacsScalar *column = &A[aptr[col]];

        // Determine the number of columns in the block
        int cdim = bptr[col + 1] - bptr[col];

        // Determine the number of non-zero block column entries in the
        // matrix to search them
        int col_size = colp[col + 1] - colp[col];
        const int *col_rows = &rows[colp[col]];

        // Loop over all the block rows in this column
        for (int i = 0, ioffset = 0; i < num_rows; i++) {
          // Find the block row index
          int row = brows[i];
          if (row < 0) {
            continue;
          }

          // Find the block row dimension
          int rdim = bptr[row + 1] - bptr[row];

          for (int ii = 0; ii < rdim; ii++) {
            int r = bptr[row] + ii;

            // Search the column to obtain the row index
            int *row_pos = TacsSearchArray(r, col_size, col_rows);
            if (row_pos) {
              int loc = row_pos - col_rows;
              TacsScalar *a = &column[cdim * loc];

              // Add the values into the matrix
              for (int jj = 0; jj < cdim; jj++) {
                a[jj] += mat[ldmat * (ioffset + ii) + (joffset + jj)];
              }
            } else {
              fprintf(stderr,
                      "BCSCMat: Entry (%d, %d) not in non-zero pattern\n", row,
                      col);
            }
          }

          // Increment the row pointer within the block matrix
          ioffset += rdim;
        }

        // Increment the column pointer within the block matrix
        joffset += cdim;
      }
    } else {
      // Add values using a linear search for each column
      for (int j = 0, joffset = 0; j < num_cols; j++) {
        // Get the block column for this portion of the matrix
        int col = bcols[j];
        if (col < 0) {
          continue;
        }
        TacsScalar *column = &A[aptr[col]];

        // Determine the size of the column and the column row indices
        int col_size = colp[col + 1] - colp[col];
        const int *col_rows = &rows[colp[col]];

        // Compute the column dimension
        int cdim = bptr[col + 1] - bptr[col];

        // Loop over all the block rows in this column
        for (int i = 0, ioffset = 0; i < num_rows; i++) {
          int row = brows[i];
          if (row < 0) {
            continue;
          }

          // Find the block row dimension
          int rdim = bptr[row + 1] - bptr[row];

          // Now find all the rows in this block
          for (int ii = 0; ii < rdim; ii++) {
            int r = bptr[row] + ii;

            // Perform a (slow) linear search in the column
            const int *row_pos = col_rows;
            for (int k = 0; k < col_size; k++) {
              if (row_pos[0] == r) {
                break;
              }
              row_pos++;
            }

            // If the row position exists, add the values into the matrix
            if (row_pos[0] == r) {
              int loc = row_pos - col_rows;
              TacsScalar *a = &column[cdim * loc];

              // Add the values into the matrix
              for (int jj = 0; jj < cdim; jj++) {
                a[jj] += mat[ldmat * (ioffset + ii) + (joffset + jj)];
              }
            } else {
              fprintf(stderr,
                      "BCSCMat: Entry (%d, %d) not in non-zero pattern\n", row,
                      col);
            }
          }

          // Increment the offset in the row dimension
          ioffset += rdim;
        }

        // Increment the offset in the column dimension
        joffset += cdim;
      }
    }
  }
}

/*
  Compute the matrix-vector product:

  Y = A*X

  where Y and X are dense matrices with row-major ordering. Note that
  typically the column dimension of Y and X will be small, e.g. O(10).

  input:
  X:           the dense input data of size (ncols, vec_bsize)
  vec_bsize:   the column dimension of the input/output matrices

  output:
  Y:           the dense output data of size (nrows, bsize)
*/
void BCSCMat::mult(TacsScalar *X, TacsScalar *Y, int vec_bsize) {
  memset(Y, 0, nrows * vec_bsize * sizeof(TacsScalar));
  multAdd(X, Y, Y, vec_bsize);
}

/*
  Compute the matrix-vector product:

  Y = Z + A*X

  where Y, Z and X are dense matrices with row-major ordering. Note
  that typically the column dimension of Y and X will be small,
  e.g. O(10). Y and Z may be the same, in which case no copying is
  performed.

  input:
  X:         the dense input data of size (ncols, vec_bsize)
  Z:         a dense input vector of size (ncols, vec_bsize)
  vec_bsize: the column dimension of the input/output matrices

  output:
  Y:         the dense output data of size (nrows, vec_bsize)
*/
void BCSCMat::multAdd(TacsScalar *X, TacsScalar *Z, TacsScalar *Y,
                      int vec_bsize) {
  if (Y != Z) {
    memcpy(Y, Z, nrows * vec_bsize * sizeof(TacsScalar));
  }

  // Allocate temporary storage for the matrix-matrix products if it
  // hasn't been allocated already
  if (!temp_array) {
    temp_array = new TacsScalar[temp_array_size];
  }

  // Determine the column block size to use
  int col_block_size = temp_array_size / vec_bsize;

  for (int j = 0; j < nblock_cols; j++) {
    // Determine the column dimenion
    int cdim = bptr[j + 1] - bptr[j];

    // Determine the starting location of the storage in this row
    TacsScalar *col = &A[aptr[j]];

    // Scan through the entire column
    int ip_end = colp[j + 1];
    for (int ip = colp[j]; ip < ip_end; ip += col_block_size) {
      // Compute the size of the column that will be computed in the
      // next set of operations: min(column_block_size, ip_end - ip)
      int col_size = col_block_size;
      if ((ip_end - ip) < col_block_size) {
        col_size = (ip_end - ip);
      }

      // We want to compute block = A*X but since the matrices are
      // stored in row-major ordering (rather than fortran col-major
      // order, we compute: block^{T} = A^{T}*X^{T}. There is no BLAS
      // routine to do this, instead we compute the equivalent:
      // block = alpha*X*A + beta*block
      TacsScalar alpha = 1.0, beta = 0.0;
      BLASgemm("N", "N", &vec_bsize, &col_size, &cdim, &alpha,
               &X[bptr[j] * vec_bsize], &vec_bsize, col, &cdim, &beta,
               temp_array, &vec_bsize);

      // Update the position of the column array
      col += cdim * col_size;

      // Now add the contributions from the termporary block to
      // the dense output matrix - this is a sparse AXPY operation
      // Y[rows[ip:col_size], :] += tblock[:col_size, :]
      addScatterColumn(Y, vec_bsize, temp_array, vec_bsize, &rows[ip], col_size,
                       vec_bsize);
    }
  }
}

/*
  Get the values from the matrix

  output:
  nrows:       the number of rows in the matrix
  ncols:       the number of columns in the matrix
  nblock_cols: the number of block columns in the matrix
  bptr:        pointer to the beginning/ends of each column block
  aptr:        pointer to the beginning of each column in A
  colp:        the pointer to each new column
  rows:        the row indices of the entries
  A:           the matrix values
*/
void BCSCMat::getArrays(int *_nrows, int *_ncols, int *_nblock_cols,
                        const int **_bptr, const int **_aptr, const int **_colp,
                        const int **_rows, TacsScalar **_A) {
  if (_nrows) {
    *_nrows = nrows;
  }
  if (_ncols) {
    *_ncols = ncols;
  }
  if (_nblock_cols) {
    *_nblock_cols = nblock_cols;
  }
  if (_bptr) {
    *_bptr = bptr;
  }
  if (_aptr) {
    *_aptr = aptr;
  }
  if (_colp) {
    *_colp = colp;
  }
  if (_rows) {
    *_rows = rows;
  }
  if (_A) {
    *_A = A;
  }
}

/*
  Retrieve the maximum block size of any given sub-block in the entire
  matrix and return it.
*/
int BCSCMat::getMaxBlockSize() { return max_block_size; }

/*
  This class performs sparse matrix factorization with blocking based
  on the nodal structure of the matrix. This blocking enables more use
  of BLAS level-3 routines that have better floating point
  performance.  The BCSCMatPivot data structure uses partial
  pivoting. Note that the storage format is arranged by groups of
  columns, but is not blocked by row - enabling easier implementation
  of pivoting but sacrificing additional performance, perhaps.
*/
BCSCMatPivot::BCSCMatPivot(BCSCMat *_mat) {
  mat = _mat;
  mat->incref();

  // Set the default fill-in
  fill = 5.0;

  // Retrieve the data from the BCSRMat data structure
  int mat_nrows, mat_ncols, mat_nblock_cols;
  mat->getArrays(&mat_nrows, &mat_ncols, &mat_nblock_cols, &bptr, NULL, NULL,
                 NULL, NULL);

  // Retrieve the maximum block size in the matrix
  max_block_size = mat->getMaxBlockSize();

  // Set the sizes of the arrays
  nrows = mat_nrows;
  ncols = mat_ncols;
  nblock_cols = mat_nblock_cols;

  // Allocate the temporary buffer
  temp_array_size = 64000;
  temp_array =
      new TacsScalar[max_block_size * max_block_size + temp_array_size];

  // Allocate the temp column
  temp_column = new TacsScalar[nrows * max_block_size];

  // Allocate the integer space required for factorization
  int size = 3 * nblock_cols + 2 * nrows;
  temp_iarray = new int[size];

  // Allocate the permutation array
  perm = new int[nrows];

  // Allocate storage for the nodal information
  var_to_node = new int[ncols];
  node_to_vars_ptr = new int[nblock_cols + 1];
  node_to_vars = new int[ncols];

  // Allocate storage for the LU factorization
  lu_aptr = new int[nblock_cols + 1];
  lu_diag_ptr = new int[nblock_cols];
  lu_col_ptr = new int[nblock_cols + 1];

  // Null-out the initial arrays
  max_lu_size = 0;
  LU = 0;

  max_lu_rows_size = 0;
  lu_rows = NULL;
}

/*
  The destructor for the matrix
*/
BCSCMatPivot::~BCSCMatPivot() {
  mat->decref();

  // If space has been allocated, deallocate it now
  if (temp_array) {
    delete[] temp_array;
  }
  if (temp_iarray) {
    delete[] temp_iarray;
  }
  if (temp_column) {
    delete[] temp_column;
  }
  if (perm) {
    delete[] perm;
  }
  if (var_to_node) {
    delete[] var_to_node;
  }
  if (node_to_vars_ptr) {
    delete[] node_to_vars_ptr;
  }
  if (node_to_vars) {
    delete[] node_to_vars;
  }
  if (LU) {
    delete[] LU;
  }
  if (lu_diag_ptr) {
    delete[] lu_diag_ptr;
  }
  if (lu_col_ptr) {
    delete[] lu_col_ptr;
  }
  if (lu_rows) {
    delete[] lu_rows;
  }
}

/*
  Determine the topological ordering of the nodes required for the
  application of L^{-1} to a sparse right-hand-side. Here, 'iteration'
  is used to create unique vertex labels that are required for the
  depth-first search. Calling this function twice with the same value
  for iteration will not work.

  This code works by finding the reach(nnz(rhs), L[:iteration,
  :iteration]^{T}), the reachable set of nodes from the vertices of
  the DAG L^{T}. These are the non-zeros of L^{-1}*rhs, and a
  topological ordering is found by performing a reverse post-order of
  the vertices found using a depth-first search. (Note that this
  function returns the post-order, which must be scanned in reverse to
  obtain the reverse post-order.)

  input:
  iteration:      used to generate labels - must be unique between calls
  rhs_rows:       the non-zero rows in the initial right-hand-side
  num_rhs_rows:   the number of rows in the right-hand-side

  output:
  topo_order:      the post-order of the nodes
  num_nodes:       the number of nodes labeled by this function

  temporary storage:
  node_stack:      the stack of nodes
  node_labels:     the labels for the nodes
*/
void BCSCMatPivot::computeTopologicalOrder(int iteration, const int *rhs_rows,
                                           const int num_rhs_rows,
                                           int *topo_order, int *_num_nodes,
                                           int *node_stack, int *node_labels) {
  // Prepare for the depth first search to determine the
  // non-zero structure of the next column vector
  // Set the vertex labels - no reset necessary
  int vertex_label_visit = 2 * iteration + 1;
  int vertex_label_done = 2 * (iteration + 1);

  // Set the number of vertices
  int num_nodes = 0;

  // Begin the depth-first search from every node in
  // the column by adding all the nodes to a stack.
  for (int j = 0; j < num_rhs_rows; j++) {
    // Check if this node exists and whether it has
    // been visited before. If it hasn't,
    // do a depth-first search starting from this node
    int root_node = var_to_node[rhs_rows[j]];
    if (root_node >= 0 && node_labels[root_node] != vertex_label_done) {
      // Set the initial stack-size as 1 - set the node
      // to be labeled as visited
      int stack_size = 1;
      node_stack[0] = root_node;
      node_stack[1] = lu_diag_ptr[root_node] + 1;
      node_labels[root_node] = vertex_label_visit;

      // While the stack has entries continue ...
      while (stack_size > 0) {
        // View the top entry in the stack
        int node = node_stack[2 * stack_size - 2];
        int jp = node_stack[2 * stack_size - 1];

        // Search through all the variables in the current column
        // below the diagonal
        for (; jp < lu_col_ptr[node + 1]; jp++) {
          int row = lu_rows[jp];
          int next_node = var_to_node[row];

          // Check that the node exists, and that it has not
          // yet been visited or created
          if (next_node >= 0 && node_labels[next_node] != vertex_label_visit &&
              node_labels[next_node] != vertex_label_done) {
            // We can start searching further down the list next time
            // we come back to this node
            node_stack[2 * stack_size - 1] = jp;

            // Set the next node in the stack
            node_labels[next_node] = vertex_label_visit;
            node_stack[2 * stack_size] = next_node;
            node_stack[2 * stack_size + 1] = lu_diag_ptr[next_node] + 1;
            stack_size++;
            break;
          }
        }

        // We've exhausted all possible nodes in the list
        if (jp == lu_col_ptr[node + 1]) {
          // Label the node as done, and update the vertex ordering
          node_labels[node] = vertex_label_done;
          topo_order[num_nodes] = node;
          num_nodes++;
          stack_size--;
        }
      }
    }
  }

  *_num_nodes = num_nodes;
}

/*
  Determine the non-zero pattern of the right hand side, given the
  nodes that will contribute an update. The updating nodes should be
  found using the 'computeTopologicalOrder' found above. This function
  obtains the non-zero pattern by determining the union of rhs and
  L[:iteration, node], for all nodes in the updating set (stored in
  the array update_nodes).

  input:
  iteration:      the current iteration count
  rhs_rows:       the rows in the non-zero pattern of the RHS
  num_rhs_rows:   the number of rows in rhs_rows
  update_nodes:   the updating nodes (the actual order doesn't matter)
  num_nodes:      the number of nodes in the updating set

  temporary storage:
  var_labels:     an array that stores the variable labels

  output:
  rhs_vars:       the non-zero variables
  num_vars:       the number of non-zero variables
*/
void BCSCMatPivot::computeColNzPattern(int iteration, const int *rhs_rows,
                                       const int num_rhs_rows,
                                       const int *update_nodes,
                                       const int num_nodes, int *var_labels,
                                       int *rhs_vars, int *_num_vars) {
  // Set the vertex label
  int vertex_label_done = 2 * (iteration + 1);

  // Keep track of the ordered variables
  int num_vars = 0;

  // Add all variables from the right hand side that are either
  // not finished (with node < 0)
  for (int j = 0; j < num_rhs_rows; j++) {
    int row = rhs_rows[j];
    int node = var_to_node[row];

    // Check that the node has not been created yet and if it
    // has not, add the row to the variable list
    if (node < 0) {
      var_labels[row] = vertex_label_done;
      rhs_vars[0] = row;
      rhs_vars++;
      num_vars++;
    }
  }

  // Add all variables from the updating nodes
  for (int i = 0; i < num_nodes; i++) {
    int node = update_nodes[i];

    // Scan through the node's column of updating rows
    for (int jp = lu_diag_ptr[node]; jp < lu_col_ptr[node + 1]; jp++) {
      int row = lu_rows[jp];
      int node = var_to_node[row];

      // Add the variable if it is required
      if (node < 0 && var_labels[row] != vertex_label_done) {
        var_labels[row] = vertex_label_done;
        rhs_vars[0] = row;
        rhs_vars++;
        num_vars++;
      }
    }
  }

  *_num_vars = num_vars;
}

/*
  Compute the update to the columns from the specified node
  columns. This code first computes the application of L on the
  non-zero entries of the column as follows:

  spa[r:s, :] <- L[r:s, r:s]^{-1}spa[r:s, :]

  Next, it computes the update to the remainder of the columns
  in f as follows:

  f[s:, :] <- f[s:, :] - L[s:, r:s]*spa[r:s, :]

  Note that the range r:s depends on the non-zero structure of the
  right-hand side of f. It should not neccessarily be the full size of
  the node, but instead should be the first non-zero entry in f.

  Note that the sparse accumulator is stored in row-major order. This
  facilitates access to the rows of the super-node.

  input:
  node:            the node to use
  node_dim:        the column dimension of the current node
  spa:             the sparse accumulator array of size nrows*spa_dim
  spa_dim:         the column dimension of the sparse accumulator array
  temp_block:      the temporary diagonal block
  temp_block_size: the size of the temp_block_array >= spa_dim*cdim
  temp_cols:       temporary storage for columns
  temp_cols_size   size of the temporary storage
*/
void BCSCMatPivot::applyNodeUpdate(int node, int node_dim, TacsScalar *spa,
                                   int spa_dim, TacsScalar *temp_block,
                                   int temp_block_size, TacsScalar *temp_cols,
                                   int temp_cols_size) {
  // Determine the location of the diagonal of the L matrix
  int loffset =
      lu_aptr[node] + node_dim * (lu_diag_ptr[node] - lu_col_ptr[node]);
  TacsScalar *L = &LU[loffset];

  // Compute the column size
  int max_col_size = temp_cols_size / spa_dim;

  // Retrieve the starting location in the sparse accumulator object
  int node_ptr = node_to_vars_ptr[node];

  // Fill in the block array with the right hand side
  const int *var = &node_to_vars[node_ptr];
  TacsScalar *tblock = &temp_block[0];

  for (int i = 0; i < node_dim; i++) {
    // Retrieve the data from the sparse accumulator
    TacsScalar *tspa = &spa[spa_dim * var[0]];
    for (int j = 0; j < spa_dim; j++) {
      tblock[0] = tspa[0];
      tblock++;
      tspa++;
    }
    var++;
  }

  // Apply the inverse of L[r:s, r:s]^{-1}*spa[r:s, :] without partial
  // pivoting or numerical checks.  Note that LU is stored in row-major
  // ordering therefore "L" is in the "upper" portion of the matrix.
  // This solves the equation: x^{T}*L^{T} = f^{T}, however, since
  // we're using row-major ordering this is equivalent to solving x*L =
  // f, (where L is stored as an upper-triangular matrix.
  TacsScalar alpha = 1.0;
  BLAStrsm("R", "U", "N", "U", &spa_dim, &node_dim, &alpha, L, &node_dim,
           temp_block, &spa_dim);

  // Place the result back into the sparse accumulator
  var = &node_to_vars[node_ptr];
  tblock = &temp_block[0];

  for (int i = 0; i < node_dim; i++) {
    // Retrieve the data from the sparse accumulator
    TacsScalar *tspa = &spa[spa_dim * var[0]];
    for (int j = 0; j < spa_dim; j++) {
      tspa[0] = tblock[0];
      tblock++;
      tspa++;
    }
    var++;
  }

  // Update the location of L
  L += node_dim * node_dim;

  // Compute the product of temp_block with the column L.  This stores
  // the array in the temporary vector, then applies the updates to the
  // rows in f
  int j_end = lu_col_ptr[node + 1];
  for (int j = lu_diag_ptr[node] + node_dim; j < j_end; j += max_col_size) {
    // Compute the size of the updating block
    int col_size = max_col_size;
    if (j_end - j < max_col_size) {
      col_size = j_end - j;
    }

    // Compute the matrix-matrix product:
    // x = L[s:, r:s]*y where y = L[r:s, r:s]^{-1}*spa[r:s, :]
    TacsScalar alpha = 1.0, beta = 0.0;
    BLASgemm("N", "N", &spa_dim, &col_size, &node_dim, &alpha, temp_block,
             &spa_dim, L, &node_dim, &beta, temp_cols, &spa_dim);

    // Update the location of the pointer to the lower block
    L += node_dim * col_size;

    // Add the result of the matrix-vector product back into the sparse
    // accumulator array
    subtractScatterColumn(spa, spa_dim, temp_cols, spa_dim, &lu_rows[j],
                          col_size, spa_dim);
  }
}

/*
  Compute the update to the columns from the upper portion of the LU
  factorization.  This code first computes the application of U on the
  non-zero entries of the column as follows:

  spa[r:s, :] <- U[r:s, r:s]^{-1}spa[r:s, :]

  Next, it computes the update to the remainder of the columns
  in f as follows:

  spa[:r, :] <- spa[:r, :] - U[:r, r:s]*spa[r:s, :]

  Note that the range r:s depends on the non-zero structure of the
  right-hand side of f. It should not neccessarily be the full size of
  the node, but instead should be the first non-zero entry in f.

  Note that the sparse accumulator is stored in row-major order. This
  facilitates access to the rows of the node.

  input:
  snode:           the super node to use
  spa:             the sparse accumulator array of size nrows*spa_width
  temp_block:      the temporary diagonal block
  temp_block_size: the size of the temp_block_array >= spa_width*bsize
  temp_cols:       temporary storage for columns
  temp_cols_size   size of the temporary storage
*/
void BCSCMatPivot::applyNodeUpperUpdate(int node, int node_dim, TacsScalar *spa,
                                        int spa_dim, TacsScalar *temp_block,
                                        int temp_block_size,
                                        TacsScalar *temp_cols,
                                        int temp_cols_size) {
  // Determine the location after the last entry of the U matrix -
  // this is the diagonal block
  int loffset =
      lu_aptr[node] + node_dim * (lu_diag_ptr[node] - lu_col_ptr[node]);
  TacsScalar *U = &LU[loffset];

  // Compute the column size
  int max_col_size = temp_cols_size / spa_dim;

  // Retrieve the starting location in the sparse accumulator object
  int node_ptr = node_to_vars_ptr[node];

  // Fill in the block array with the right hand side
  const int *var = &node_to_vars[node_ptr];
  TacsScalar *tblock = temp_block;

  for (int i = 0; i < node_dim; i++) {
    // Retrieve the data from the sparse accumulator
    TacsScalar *tspa = &spa[spa_dim * var[0]];
    for (int j = 0; j < spa_dim; j++) {
      tblock[0] = tspa[0];
      tblock++;
      tspa++;
    }
    var++;
  }

  // Apply the inverse of U[r:s, r:s]^{-1}*spa[r:s, :] without partial
  // pivoting or numerical checks.  Note that LU is stored in row-major
  // ordering therefore "U" is in the apparent "lower" portion of the
  // matrix from LAPACK's perspective (column-major order). This code
  // solves the equation: U^{T}*x^{T} = f^{T}, however, since
  // we're using row-major ordering this is equivalent to solving x*U = f
  TacsScalar alpha = 1.0;
  BLAStrsm("R", "L", "N", "N", &spa_dim, &node_dim, &alpha, U, &node_dim,
           temp_block, &spa_dim);

  // Place the result back into the sparse accumulator
  var = &node_to_vars[node_ptr];
  tblock = &temp_block[0];

  for (int i = 0; i < node_dim; i++) {
    // Retrieve the data from the sparse accumulator
    TacsScalar *tspa = &spa[spa_dim * var[0]];
    for (int j = 0; j < spa_dim; j++) {
      tspa[0] = tblock[0];
      tblock++;
      tspa++;
    }
    var++;
  }

  // Compute the product of temp_block with the column U.  This stores
  // the array in the temporary vector, then applies the updates to the
  // rows in f
  int j_end = lu_diag_ptr[node];

  // Reset the location of the U pointer
  U = &LU[lu_aptr[node]];
  for (int j = lu_col_ptr[node]; j < j_end; j += max_col_size) {
    // Compute the size of the updating block
    int col_size = max_col_size;
    if (j_end - j < max_col_size) {
      col_size = j_end - j;
    }

    // Compute the matrix-matrix product:
    // x = U[s:, r:s]*y where y = U[r:s, r:s]^{-1}*f[r:s, :]
    TacsScalar alpha = 1.0, beta = 0.0;
    BLASgemm("N", "N", &spa_dim, &col_size, &node_dim, &alpha, temp_block,
             &spa_dim, U, &node_dim, &beta, temp_cols, &spa_dim);

    // Update the location of the pointer to the U block
    U += node_dim * col_size;

    // Add the result of the matrix-vector product back into the
    // sparse accumulator array
    subtractScatterColumn(spa, spa_dim, temp_cols, spa_dim, &lu_rows[j],
                          col_size, spa_dim);
  }
}

/*
  Factor the node and perform any pivoting operations that are
  required. Note that at this point, all columns within this node have
  the same non-zero structure, this greatly simplifies the
  update/pivoting method.

  At the moment, this code is implemented without BLAS, but several
  operations in this function could be converted to use level-2 and
  level-1 BLAS operations.

  This modifies the following data:
  perm:             the permutation sequence
  var_to_nodes:     for each new variable added to the node
  node_to_vars:     for the current node, update
  node_to_vars_ptr: pointer to the node_to_vars variable

  input:
  node:         the node number that is being factored
  col:          pointer to the panel values being factored
  node_dim:     width of the panel
  rows:         row indices of the panel
  num_rows:     number of rows in the panel
  diag_index:   index of the first diagonal entry of the left-most column
*/
void BCSCMatPivot::factorNode(int node, TacsScalar *col, int node_dim,
                              int *rows, int num_rows, int diag_index) {
  int init_diag_index = diag_index;

  // Update the data structure for this newly constructed node
  node_to_vars_ptr[node + 1] = node_to_vars_ptr[node];

  for (int j = 0; j < node_dim; diag_index++, j++) {
    // Apply the update from the previous rows
    if (j > 0) {
      // Set pointers to the portion of the column corresponding to
      // the lower diagonal matrix
      TacsScalar *L = &col[node_dim * init_diag_index];
      TacsScalar *U = &col[node_dim * init_diag_index + j];

      // Solve the problem L[init_diag:diag, init_diag:diag]^{-1}*U
      BLAStrsv("U", "T", "U", &j, L, &node_dim, U, &node_dim);

      // Update the remaining portion of the panel
      for (int i = diag_index; i < num_rows; i++) {
        for (int k = 0; k < j; k++) {
          col[node_dim * i + j] -= col[node_dim * i + k] * U[k * node_dim];
        }
      }
    }

    // Determine the maximum value in the current column
    int max_row_index = diag_index;
    TacsScalar *c = &col[node_dim * diag_index + j];
    TacsScalar diag_entry = c[0];
    TacsScalar max_entry = c[0];
    c += node_dim;

    for (int i = diag_index + 1; i < num_rows; i++) {
      double abs_val = fabs(TacsRealPart(c[0]));
      if (abs_val > fabs(TacsRealPart(max_entry))) {
        max_row_index = i;
        max_entry = c[0];
      }
      c += node_dim;
    }

    if (max_entry == 0.0) {
      fprintf(stderr, "BCSCMatPivot: Error, zero pivot in block column %d\n",
              node);
    } else if (max_entry != max_entry) {
      fprintf(stderr, "BCSCMatPivot: Error, NaN pivot in block column %d\n",
              node);
    }

    // Record the pivot variable
    int pivot = rows[diag_index];
    TacsScalar entry = diag_entry;

    // Pivot only if required by the pivot tolerance
    if (fabs(TacsRealPart(max_entry)) > fabs(TacsRealPart(diag_entry))) {
      // Set the pivot and entry values
      pivot = rows[max_row_index];
      entry = max_entry;

      // Swap rows diag and max_row_index
      if (diag_index != max_row_index) {
        // Swap the index itself
        int temp_row = rows[max_row_index];
        rows[max_row_index] = rows[diag_index];
        rows[diag_index] = temp_row;

        // Swap the values in the rows diag/max_row_index
        for (int i = 0; i < node_dim; i++) {
          TacsScalar temp = col[node_dim * max_row_index + i];
          col[node_dim * max_row_index + i] = col[node_dim * diag_index + i];
          col[node_dim * diag_index + i] = temp;
        }
      }
    }

    // Compute the panel-level factorization
    c = &col[node_dim * (diag_index + 1) + j];
    for (int i = diag_index + 1; i < num_rows; i++) {
      c[0] = c[0] / entry;
      c += node_dim;
    }

    // Record the pivot sequence
    perm[pivot] = node_to_vars_ptr[node + 1];

    // Record the node->variables array
    node_to_vars[node_to_vars_ptr[node + 1]] = pivot;
    node_to_vars_ptr[node + 1]++;

    // Record the var->node information
    var_to_node[pivot] = node;
  }
}

/*
  Factor the matrix using a node-based technique with node-column
  updates. This call computes P*A = L*U, a factorization of the matrix
  using partial pivoting. This code uses an estimate of the predicted
  fill in. The code attempts to extend the arrays if the esimated fill
  is exceeded. If it cannot allocate additional space, it returns a
  negative actual fill: -1.0, otherwise it returns the actual amount of
  fill used in the computation.

  The algorithm to compute P*A = L*U is as follows

  for j in range(0, ncols):
  .   f = A[:, j:j+nstep] # Sparse copy
  .   # symbolic factorization: determine the non-zero pattern
  .   # of the following components of the matrix:
  .
  .   nz(L[:j, :j]^{-1} f[:j]) and
  .   nz(f[:j] + L[:j, j:]*f[j:])

  .   for each supernode r:s in topological order:
  .       f[r:s] = L[r:s, r:s]^{-1}*f[r:s]
  .       f[s:] -= L[s:, r:s]*f[r:s]

  .   # complete sparse factorization on the panel
  .   for k in range(nstep):
  .       jpiv = argmax(f[j+k:, j+k])   # compute pivot row
  .       f[jpiv] <-> f[j+k]            # swap max row
  .       perm[jpiv] = j+k              # record the pivot
  .       U[:j+k+1] = f[:j+k+1]         # copy over the column of U
  .       L[j+k+1:] = f[j+k+1:]/f[j+k]  # copy over the column of L

  input:
  mat:  the input BCSCMat matrix to be factored
  fill: the esimated fill level nnz(LU)/nnz(A) > 1

  returns: the actual fill in from the fully factored matrix
*/
double BCSCMatPivot::factor(double _fill) {
  // Retrieve the data from the BCSRMat data structure
  int mat_nrows, mat_ncols, mat_nblock_cols;
  const int *mat_bptr, *mat_aptr;
  const int *mat_colp, *mat_rows;
  TacsScalar *A;
  mat->getArrays(&mat_nrows, &mat_ncols, &mat_nblock_cols, &mat_bptr, &mat_aptr,
                 &mat_colp, &mat_rows, &A);

  // For now, hard code this estimate
  int max_lu_size = fill * mat_aptr[mat_nblock_cols];

  // First, if LU/lu_rows have not been allocated, allocate them
  if (!LU) {
    // Update the fill based on the input fill-in
    if (_fill > 1.0) {
      fill = _fill;
    }

    // Set the lu size and allocate the array
    max_lu_size = fill * mat_aptr[nblock_cols];
    LU = new TacsScalar[max_lu_size];

    max_lu_rows_size = fill * mat_colp[nblock_cols];
    lu_rows = new int[max_lu_rows_size];
  }

  // Set the var -> node data structure
  for (int i = 0; i < ncols; i++) {
    var_to_node[i] = -1;
  }

  // Initialize the node -> var data structure
  memset(node_to_vars_ptr, 0, (nblock_cols + 1) * sizeof(int));
  for (int i = 0; i < ncols; i++) {
    node_to_vars[i] = -1;
  }

  // Set the initial permutation
  for (int i = 0; i < nrows; i++) {
    perm[i] = -1;
  }

  // Initialize the temporary storage
  TacsScalar *temp_block = &temp_array[0];
  TacsScalar *temp_cols = &temp_array[max_block_size * max_block_size];

  // Initialize the sparse accumulator object
  TacsScalar *column = temp_column;
  memset(column, 0, nrows * max_block_size * sizeof(TacsScalar));

  // Allocate temporary arrays for data used during the factorization
  int *node_stack = &temp_iarray[0];
  int *topo_order = &temp_iarray[2 * nblock_cols];
  int *node_labels = &temp_iarray[2 * nblock_cols + nrows];
  int *var_labels = &temp_iarray[3 * nblock_cols + nrows];

  // Set the labels
  memset(node_labels, 0, nblock_cols * sizeof(int));
  memset(var_labels, 0, nrows * sizeof(int));

  // Set the initial positions within the array
  lu_aptr[0] = 0;
  lu_col_ptr[0] = 0;
  lu_diag_ptr[0] = 0;

  for (int node = 0; node < nblock_cols; node++) {
    // Copy the values from the next group of columns to the matrix
    int node_dim = bptr[node + 1] - bptr[node];
    int col_size = mat_colp[node + 1] - mat_colp[node];

    for (int ip = mat_colp[node]; ip < mat_colp[node + 1];
         ip++, A += node_dim) {
      // Copy over the values to the temporary column
      memcpy(&column[node_dim * mat_rows[ip]], A,
             node_dim * sizeof(TacsScalar));
    }

    // Compute the topological post-order of the nodes
    int num_nodes = 0;
    computeTopologicalOrder(node, &mat_rows[mat_colp[node]], col_size,
                            topo_order, &num_nodes, node_stack, node_labels);

    // Compute the non-zero pattern of the
    int num_vars = 0;
    computeColNzPattern(node, &mat_rows[mat_colp[node]], col_size, topo_order,
                        num_nodes, var_labels, &topo_order[num_nodes],
                        &num_vars);

    // Compute the additional space required for the non-zero pattern
    int nnz_rows = num_vars;
    for (int i = 0; i < num_nodes; i++) {
      int node = topo_order[i];
      nnz_rows += node_to_vars_ptr[node + 1] - node_to_vars_ptr[node];
    }

    // Check if this new column will exceed the available space
    if (node_dim * nnz_rows + lu_aptr[node] > max_lu_size) {
      // Determine an approximate new max size
      max_lu_size = node_dim * nnz_rows + 2 * max_lu_size;

      // Allocate a new LU array and free the old one
      TacsScalar *old_LU = LU;
      LU = new (std::nothrow) TacsScalar[max_lu_size];
      if (LU) {
        memcpy(LU, old_LU, lu_aptr[node] * sizeof(TacsScalar));
      }
      delete[] old_LU;

      // If we can't allocate new space, return a fail condition
      // and free the temporary variables
      if (!LU) {
        return -1.0;
      }
    }

    // Check if this new column will exceed the available colp space
    if (nnz_rows + lu_col_ptr[node] > max_lu_rows_size) {
      max_lu_rows_size = nnz_rows + 2 * max_lu_rows_size;

      // Allocate a new lu_rows array and free the old one
      int *old_rows = lu_rows;
      lu_rows = new (std::nothrow) int[max_lu_rows_size];
      if (lu_rows) {
        memcpy(lu_rows, old_rows, lu_col_ptr[node] * sizeof(int));
      }
      delete[] old_rows;

      // If we can't allocate new space, return a fail condition
      // and free the temporary variables
      if (!lu_rows) {
        return -1.0;
      }
    }

    // Apply the numerical updates to the spa array
    for (int i = num_nodes - 1; i >= 0; i--) {
      int node = topo_order[i];
      int cdim = bptr[node + 1] - bptr[node];
      applyNodeUpdate(node, cdim, column, node_dim, temp_block,
                      node_dim * node_dim, temp_cols, temp_array_size);
    }

    // Extract the matrix values from the the temporary sparse
    // accumulator and set the non-zero pattern in the LU data
    // Keep track of the element of the matrix added last
    lu_col_ptr[node + 1] = lu_col_ptr[node] + nnz_rows;

    // Update the pointer to the next column
    lu_aptr[node + 1] = lu_aptr[node] + node_dim * nnz_rows;

    // Set the pointer into the LU vector
    TacsScalar *lu = &LU[lu_aptr[node]];
    int *lu_row = &lu_rows[lu_col_ptr[node]];

    // Scan through all the new data in topological order
    for (int i = num_nodes - 1; i >= 0; i--) {
      int node = topo_order[i];
      for (int k = node_to_vars_ptr[node]; k < node_to_vars_ptr[node + 1];
           k++) {
        int row = node_to_vars[k];

        // Copy the data from the row
        memcpy(lu, &column[node_dim * row], node_dim * sizeof(TacsScalar));
        lu += node_dim;

        // Zero the corresponding row in the dense column
        memset(&column[node_dim * row], 0, node_dim * sizeof(TacsScalar));

        // Set the row index
        lu_row[0] = row;
        lu_row++;
      }
    }

    // Set the pointer to the diagonal entry of the LU factorization
    lu_diag_ptr[node] = lu_row - lu_rows;

    for (int i = num_nodes; i < num_vars + num_nodes; i++) {
      int row = topo_order[i];

      // Copy the data from the row
      memcpy(lu, &column[node_dim * row], node_dim * sizeof(TacsScalar));
      lu += node_dim;

      // Zero the corresponding row in 'column'
      memset(&column[node_dim * row], 0, node_dim * sizeof(TacsScalar));

      // Set the row index and increment the pointer
      lu_row[0] = row;
      lu_row++;
    }

    // Now all the updates have been applied to the temporary column
    // vector and the non-zero entries are labeled in the correct order.
    // Next, apply the factorization to the column.
    int col_ptr = lu_col_ptr[node];
    int num_rows = lu_col_ptr[node + 1] - col_ptr;
    int diag_index = lu_diag_ptr[node] - col_ptr;
    factorNode(node, &LU[lu_aptr[node]], node_dim, &lu_rows[col_ptr], num_rows,
               diag_index);
  }

  return (1.0 * lu_col_ptr[nblock_cols]) / mat_colp[nblock_cols];
}

/*
  Apply the pivot sequence and the lower/upper factorization to the
  input vector to obtain X <- U^{-1}*L^{-1}*P*X.

  input:
  X:         the dense matrix of size (nrows, bsize)
  vec_size:  the column dimension of the input vector

  output:
  X:        the dense matrix of size (nrows, bsize)
*/
void BCSCMatPivot::applyFactor(TacsScalar *X, int vec_bsize) {
  if (vec_bsize == 1) {
    memcpy(temp_column, X, nrows * sizeof(TacsScalar));

    // Apply the lower/upper parts of the factorization
    applyLower(temp_column, 1, temp_array, temp_array_size);
    applyUpper(temp_column, 1, temp_array, temp_array_size);

    for (int i = 0; i < nrows; i++) {
      X[perm[i]] = temp_column[i];
    }
  } else {
    // Transform from column-major order to row-major order
    // for each block
    int num_rhs_cols = vec_bsize;

    while (num_rhs_cols > 0) {
      int num_cols = num_rhs_cols;
      if (num_cols > max_block_size) {
        num_cols = max_block_size;
      }

      // Keep a pointer to where were currently copying
      // out the values
      TacsScalar *xcopy = X;

      // Copy the tranpose values to temp_array
      for (int j = 0; j < num_cols; j++) {
        for (int i = 0; i < nrows; i++) {
          temp_column[j + i * num_cols] = xcopy[0];
          xcopy++;
        }
      }

      // Apply the lower/upper parts of the factorization
      applyLower(temp_column, num_cols, temp_array, temp_array_size);
      applyUpper(temp_column, num_cols, temp_array, temp_array_size);

      // Transpose the answer and copy the result back into the array
      TacsScalar *c = temp_column;
      for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < vec_bsize; j++) {
          X[nrows * j + perm[i]] = c[0];
          c++;
        }
      }

      X += num_cols * nrows;
      num_rhs_cols -= num_cols;
    }
  }
}

/*
  Apply the lower factorization to a dense column.

  Note that the matrix B is stored in row-major ordering.

  input:
  B:        the dense matrix of size (nrows, bsize)

  output:
  B = L^{-1}*B
*/
void BCSCMatPivot::applyLower(TacsScalar *B, int vec_bsize, TacsScalar *temp,
                              int temp_size) {
  TacsScalar *temp_block = &temp[0];
  TacsScalar *temp_cols = &temp[vec_bsize * max_block_size];

  int temp_block_size = vec_bsize * max_block_size;
  int temp_cols_size = temp_size - temp_block_size;

  for (int node = 0; node < nblock_cols; node++) {
    int node_dim = bptr[node + 1] - bptr[node];
    // Compute L^{-1}*B
    applyNodeUpdate(node, node_dim, B, vec_bsize, temp_block, temp_block_size,
                    temp_cols, temp_cols_size);
  }
}

/*
  Apply the upper factorization to a dense column

  Note that the matrix B is stored in row-major ordering.

  input:
  B;        the dense matrix of size (nrows, bsize)

  output:
  B = U^{-1}*B
*/
void BCSCMatPivot::applyUpper(TacsScalar *B, int vec_bsize, TacsScalar *temp,
                              int temp_size) {
  TacsScalar *temp_block = &temp[0];
  TacsScalar *temp_cols = &temp[vec_bsize * max_block_size];

  int temp_block_size = vec_bsize * max_block_size;
  int temp_cols_size = temp_size - temp_block_size;

  for (int node = nblock_cols - 1; node >= 0; node--) {
    int node_dim = bptr[node + 1] - bptr[node];
    // Compute U^{-1}*B
    applyNodeUpperUpdate(node, node_dim, B, vec_bsize, temp_block,
                         temp_block_size, temp_cols, temp_cols_size);
  }
}
