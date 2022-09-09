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

#include "TACSMatDistribute.h"

#include "TacsUtilities.h"

/*
  Given the non-zero pattern of a local matrix, the global indices
  of the entries and the variable map defining the parallel layout of
  the data, construct a distributed global matrix.

  This prodecure requires a considerable amount of indexing working.
  The goal is to determine two parts of the non-zero matrix pattern.
  On each processor, the non-zero pattern is stored as follows:

  A_i*x_i + B_ext*x_ext = b_i,

  where x_i are the local vector components, and x_ext are all
  external vector components stored on other processors. This
  constructor takes as input a variable layout (rmap), a non-zero
  pattern on each processor, represented by a CSR data structure with
  rowp/cols, and an index map that represents the indexing between the
  local CSR data structure and the global ordering of variables.

  This information is used to construct the non-zero pattern of A_i
  and B_ext, as well as initialize the communication paths required to
  set the values of A_i and B_ext. These communication paths assume
  that the components are set with the same non-zero pattern used to
  construct the DistMat, or a sub-set of it.

  input:
  rmap:          variable mapping assignment to processors
  num_nodes:     the number of nodes on this processor
  rowp:          the pointer to the cols for the non-zero pattern
  cols:          the column indices for the non-zero pattern
  bindex:        the global block indices of the CSR data structure

  The underlying PMat matrix is assembled as follows:

  1. The ownership intervals of the sorted ext_row_vars
  2. The non-local row numbers are sent to non-local procs
  3. The non-local column numbers are sent to non-local procs
  4. Temporary storage arrays are initialized
  5. The CSR data structure is split into local and non-local parts
  6. The non-local CSR data structure is sent to non-local procs
  7. The CSR data structures from all procs are merged
  8. The persistent data communication paths are initialized
*/
TACSMatDistribute::TACSMatDistribute(TACSThreadInfo *thread_info,
                                     TACSNodeMap *_row_map, int _bsize,
                                     int numNodes, const int *rowp,
                                     const int *cols, TACSBVecIndices *rowIndex,
                                     BCSRMat **Aloc, BCSRMat **Bext,
                                     TACSBVecDistribute **colMap) {
  // Copy over the row map
  row_map = _row_map;
  row_map->incref();

  // Get the rank/size of this communicator
  int mpiRank, mpiSize;
  comm = row_map->getMPIComm();
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);

  // Copy over the block size for the matrix
  bsize = _bsize;

  // Get the external variables
  const int *ext_vars;

  if (rowIndex->getIndices(&ext_vars) != numNodes) {
    fprintf(stderr,
            "[%d] TACSMatDistribute error: number of indices "
            "provided must be equal to the number of rows\n",
            mpiRank);
    return;
  }

  // Get the owner range for the number of variables owned by
  // each node
  const int *ownerRange;
  row_map->getOwnerRange(&ownerRange);

  // Get the non-zero pattern of the input matrix
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank + 1];

  // Count up the number of equations that must be sent to other
  // processors. Get only the variables associated with those ones.
  num_ext_rows = 0;
  for (int k = 0; k < numNodes; k++) {
    if (ext_vars[k] < lower || ext_vars[k] >= upper) {
      num_ext_rows++;
    }
  }

  ext_rows = new int[num_ext_rows];
  num_ext_rows = 0;
  for (int k = 0; k < numNodes; k++) {
    if (ext_vars[k] < lower || ext_vars[k] >= upper) {
      ext_rows[num_ext_rows] = ext_vars[k];
      num_ext_rows++;
    }
  }

  // Check that the row indices provided are unique and sort them
  // in preparation
  int num_ext_row_init = num_ext_rows;
  num_ext_rows = TacsUniqueSort(num_ext_rows, ext_rows);
  if (num_ext_row_init != num_ext_rows) {
    fprintf(stderr, "[%d] TACSMatDistrubte error: ext_vars are not unique\n",
            mpiRank);
  }

  // Create the off-process CSR data structure that will be sent to
  // other processes. First, calculate ext_rowp
  ext_rowp = new int[num_ext_rows + 1];
  memset(ext_rowp, 0, (num_ext_rows + 1) * sizeof(int));

  for (int i = 0; i < numNodes; i++) {
    if (ext_vars[i] < lower || ext_vars[i] >= upper) {
      int *item = TacsSearchArray(ext_vars[i], num_ext_rows, ext_rows);
      int index = item - ext_rows;
      ext_rowp[index + 1] = rowp[i + 1] - rowp[i];
    }
  }

  for (int i = 0; i < num_ext_rows; i++) {
    ext_rowp[i + 1] += ext_rowp[i];
  }

  // Next, calculate ext_cols. Find only the external rows and
  // convert them to the global numbering scheme
  ext_cols = new int[ext_rowp[num_ext_rows]];

  for (int i = 0; i < numNodes; i++) {
    if (ext_vars[i] < lower || ext_vars[i] >= upper) {
      int *item = TacsSearchArray(ext_vars[i], num_ext_rows, ext_rows);
      int index = item - ext_rows;

      for (int j = ext_rowp[index], jj = rowp[i]; j < ext_rowp[index + 1];
           j++, jj++) {
        ext_cols[j] = ext_vars[cols[jj]];
      }

      int size = ext_rowp[index + 1] - ext_rowp[index];
      if (size != TacsUniqueSort(size, &ext_cols[ext_rowp[index]])) {
        fprintf(stderr, "[%d] TACSMatDistribute error: Array is not unique\n",
                mpiRank);
      }
    }
  }

  // Match the intervals of the external variables to be sent to other
  // processes
  int *ext_ptr = new int[mpiSize + 1];
  TacsMatchIntervals(mpiSize, ownerRange, num_ext_rows, ext_rows, ext_ptr);

  // Count up the processors that will be sending information
  num_ext_procs = 0;
  for (int k = 0; k < mpiSize; k++) {
    int num_rows = ext_ptr[k + 1] - ext_ptr[k];
    if (num_rows > 0) {
      num_ext_procs++;
    }
  }

  // Find the external processors
  ext_procs = new int[num_ext_procs];
  ext_count = new int[num_ext_procs];
  ext_row_ptr = new int[num_ext_procs + 1];
  ext_row_ptr[0] = 0;
  for (int k = 0, count = 0; k < mpiSize; k++) {
    int num_rows = ext_ptr[k + 1] - ext_ptr[k];
    if (num_rows > 0) {
      ext_procs[count] = k;
      ext_count[count] = num_rows;
      ext_row_ptr[count + 1] = ext_row_ptr[count] + num_rows;
      count++;
    }
  }

  // Adjust the pointer array so that it is an array of counts
  for (int k = 0, offset = 0; k < mpiSize; k++) {
    ext_ptr[k] = ext_ptr[k + 1] - offset;
    offset = ext_ptr[k + 1];
  }

  // Allocate space to store the number of in-coming entries
  int *in_full = new int[mpiSize];
  MPI_Alltoall(ext_ptr, 1, MPI_INT, in_full, 1, MPI_INT, comm);
  delete[] ext_ptr;

  // Count up the number of processors contributing entries to
  // this procesor
  num_in_procs = 0;
  for (int k = 0; k < mpiSize; k++) {
    if (in_full[k] > 0) {
      num_in_procs++;
    }
  }

  // Allocate space to store which processors are sendin
  in_procs = new int[num_in_procs];
  in_count = new int[num_in_procs];
  in_row_ptr = new int[num_in_procs + 1];
  in_row_ptr[0] = 0;
  for (int k = 0, count = 0; k < mpiSize; k++) {
    if (in_full[k] > 0) {
      in_procs[count] = k;
      in_count[count] = in_full[k];
      in_row_ptr[count + 1] = in_row_ptr[count] + in_count[count];
      count++;
    }
  }

  // Prepare to receive the equation numbers from the other processes
  in_row_ptr = new int[num_in_procs + 1];
  in_row_ptr[0] = 0;
  for (int k = 0; k < num_in_procs; k++) {
    in_row_ptr[k + 1] = in_row_ptr[k] + in_count[k];
  }

  // Allocate space for the integers
  num_in_rows = in_row_ptr[num_in_procs];
  in_rows = new int[num_in_rows];

  // Allocate space for the statuses
  in_requests = new MPI_Request[num_in_procs];
  ext_requests = new MPI_Request[num_ext_procs];

  // Post the recvs
  for (int k = 0, offset = 0; k < num_in_procs; k++) {
    int count = in_count[k];
    int source = in_procs[k];
    int tag = 1;
    MPI_Irecv(&in_rows[offset], count, MPI_INT, source, tag, comm,
              &in_requests[k]);
    offset += count;
  }

  // Post the sends
  for (int k = 0, offset = 0; k < num_ext_procs; k++) {
    int count = ext_count[k];
    int dest = ext_procs[k];
    int tag = 1;
    MPI_Isend(&ext_rows[offset], count, MPI_INT, dest, tag, comm,
              &ext_requests[k]);
    offset += count;
  }

  // Wait until everything completes
  if (num_ext_procs > 0) {
    MPI_Waitall(num_ext_procs, ext_requests, MPI_STATUSES_IGNORE);
  }
  if (num_in_procs > 0) {
    MPI_Waitall(num_in_procs, in_requests, MPI_STATUSES_IGNORE);
  }

  // Prepare to pass information from ext_rowp to in_rowp
  int *ext_row_var_count = new int[num_ext_rows];
  for (int k = 0; k < num_ext_rows; k++) {
    ext_row_var_count[k] = ext_rowp[k + 1] - ext_rowp[k];
  }

  in_rowp = new int[num_in_rows + 1];
  in_rowp[0] = 0;

  // Post the recvs
  for (int k = 0, offset = 1; k < num_in_procs; k++) {
    int count = in_count[k];
    int source = in_procs[k];
    int tag = 2;
    MPI_Irecv(&in_rowp[offset], count, MPI_INT, source, tag, comm,
              &in_requests[k]);
    offset += count;
  }

  // Post the sends
  for (int k = 0, offset = 0; k < num_ext_procs; k++) {
    int count = ext_count[k];
    int dest = ext_procs[k];
    int tag = 2;
    MPI_Isend(&ext_row_var_count[offset], count, MPI_INT, dest, tag, comm,
              &ext_requests[k]);
    offset += count;
  }

  if (num_ext_procs > 0) {
    MPI_Waitall(num_ext_procs, ext_requests, MPI_STATUSES_IGNORE);
  }
  if (num_in_procs > 0) {
    MPI_Waitall(num_in_procs, in_requests, MPI_STATUSES_IGNORE);
  }

  // Convert the counts to a pointer offset
  for (int k = 0; k < num_in_rows; k++) {
    in_rowp[k + 1] += in_rowp[k];
  }

  delete[] ext_row_var_count;

  // Pass ext_cols to in_cols
  in_cols = new int[in_rowp[num_in_rows]];

  // Post the recvs
  for (int k = 0, offset = 0, buff_offset = 0; k < num_in_procs; k++) {
    int count = in_rowp[offset + in_count[k]] - in_rowp[offset];
    int source = in_procs[k];
    int tag = 3;
    MPI_Irecv(&in_cols[buff_offset], count, MPI_INT, source, tag, comm,
              &in_requests[k]);
    offset += in_count[k];
    buff_offset += count;
  }

  // Post the sends
  for (int k = 0, offset = 0, buff_offset = 0; k < num_ext_procs; k++) {
    int count = ext_rowp[offset + ext_count[k]] - ext_rowp[offset];
    int dest = ext_procs[k];
    int tag = 3;
    MPI_Isend(&ext_cols[buff_offset], count, MPI_INT, dest, tag, comm,
              &ext_requests[k]);
    offset += ext_count[k];
    buff_offset += count;
  }

  if (num_ext_procs > 0) {
    MPI_Waitall(num_ext_procs, ext_requests, MPI_STATUSES_IGNORE);
  }
  if (num_in_procs > 0) {
    MPI_Waitall(num_in_procs, in_requests, MPI_STATUSES_IGNORE);
  }

  // estimate the non-zero entries per row
  int nz_per_row = 10;
  if (numNodes > 0) {
    nz_per_row = (int)(rowp[numNodes] / numNodes + 1);
  }

  int max_col_vars_size = num_ext_rows;
  for (int i = 0; i < num_in_rows; i++) {
    for (int j = in_rowp[i]; j < in_rowp[i + 1]; j++) {
      if (in_cols[j] < ownerRange[mpiRank] ||
          in_cols[j] >= ownerRange[mpiRank + 1]) {
        max_col_vars_size++;
      }
    }
  }

  int col_vars_size = 0;
  int *temp_col_vars = new int[max_col_vars_size];

  // Get contributions from ext_rows
  for (int i = 0; i < num_ext_rows; i++) {
    temp_col_vars[col_vars_size] = ext_rows[i];
    col_vars_size++;
  }

  for (int i = 0; i < num_in_rows; i++) {
    for (int j = in_rowp[i]; j < in_rowp[i + 1]; j++) {
      if (in_cols[j] < ownerRange[mpiRank] ||
          in_cols[j] >= ownerRange[mpiRank + 1]) {
        temp_col_vars[col_vars_size] = in_cols[j];
        col_vars_size++;
      }
    }
  }

  // Uniquely sort the array
  col_vars_size = TacsUniqueSort(col_vars_size, temp_col_vars);

  int *col_vars = new int[col_vars_size];
  memcpy(col_vars, temp_col_vars, col_vars_size * sizeof(int));
  delete[] temp_col_vars;

  TACSBVecIndices *col_indices = new TACSBVecIndices(&col_vars, col_vars_size);
  *colMap = new TACSBVecDistribute(row_map, col_indices);
  col_map_size = col_indices->getIndices(&col_map_vars);

  // Assemble everything into on and off-diagonal parts
  int *Arowp, *Acols;  // The diagonal entries
  int np, *Browp, *Bcols;

  computeLocalCSR(numNodes, ext_vars, rowp, cols, ownerRange[mpiRank],
                  ownerRange[mpiRank + 1], nz_per_row, &Arowp, &Acols, &np,
                  &Browp, &Bcols);

  // Allocate the local/external matrices
  int n = ownerRange[mpiRank + 1] - ownerRange[mpiRank];
  *Aloc = new BCSRMat(comm, thread_info, bsize, n, n, &Arowp, &Acols);
  *Bext = new BCSRMat(comm, thread_info, bsize, n - np, col_vars_size, &Browp,
                      &Bcols);

  // Allocate space for in-coming matrix elements
  ext_A = new TacsScalar[bsize * bsize * ext_rowp[num_ext_rows]];
  in_A = new TacsScalar[bsize * bsize * in_rowp[num_in_rows]];
  zeroEntries();
}

TACSMatDistribute::~TACSMatDistribute() {
  row_map->decref();
  delete[] ext_procs;
  delete[] ext_count;
  delete[] ext_row_ptr;
  delete[] ext_rows;
  delete[] ext_rowp;
  delete[] ext_cols;
  delete[] ext_A;
  delete[] ext_requests;

  delete[] in_procs;
  delete[] in_count;
  delete[] in_row_ptr;
  delete[] in_rows;
  delete[] in_rowp;
  delete[] in_cols;
  delete[] in_A;
  delete[] in_requests;
}

void TACSMatDistribute::zeroEntries() {
  int ext_len = bsize * bsize * ext_rowp[num_ext_rows];
  memset(ext_A, 0, ext_len * sizeof(TacsScalar));

  int in_len = bsize * bsize * in_rowp[num_in_rows];
  memset(in_A, 0, in_len * sizeof(TacsScalar));
}

/*!
  Find the non-zero pattern for the A = [B, E; F, C] and Bext
  matrices such that:

  [ B, E ][ x ]
  [ F, C ][ y ] + [ Bext ][ y_ext ] = 0
*/
void TACSMatDistribute::computeLocalCSR(int numNodes, const int *ext_vars,
                                        const int *rowp, const int *cols,
                                        int lower, int upper, int nz_per_row,
                                        int **_Arowp, int **_Acols, int *_np,
                                        int **_Browp, int **_Bcols) {
  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);

  // For each row, find the non-zero pattern
  int A_max_row_size = 2 * nz_per_row;
  int *A_row_vars = new int[A_max_row_size];

  int B_max_row_size = 2 * nz_per_row;
  int *B_row_vars = new int[B_max_row_size];

  int A_rows = upper - lower;
  int *A_rowp = new int[A_rows + 1];
  int *B_rowp = new int[A_rows + 1];

  for (int i = 0; i < A_rows + 1; i++) {
    A_rowp[i] = 0;
    B_rowp[i] = 0;
  }

  // Set up a temporary array to keep track of the variables in ext_vars
  int *var_map = new int[A_rows];
  for (int i = 0; i < A_rows; i++) {
    var_map[i] = -1;
  }

  for (int i = 0; i < numNodes; i++) {
    int var = ext_vars[i];
    if (var >= lower && var < upper) {
      var_map[var - lower] = i;
    }
  }

  // Size the A_rowp/B_rowp arrays
  for (int i = 0; i < A_rows; i++) {  // For each variable
    int A_row_size = 0;
    int B_row_size = 0;

    int ei = var_map[i];  // Index into the external variable map
    int var = lower + i;

    if (ei >= 0) {
      // Add variables in this range to the row as well
      int start = rowp[ei];
      int end = rowp[ei + 1];

      // Add rowp[ei], to A_row_vars, B_row_vars
      if (A_row_size + end - start > A_max_row_size) {
        A_max_row_size = A_max_row_size + end - start;
        TacsExtendArray(&A_row_vars, A_row_size, A_max_row_size);
      }

      if (B_row_size + end - start > B_max_row_size) {
        B_max_row_size = B_max_row_size + end - start;
        TacsExtendArray(&B_row_vars, B_row_size, B_max_row_size);
      }

      for (int j = start; j < end; j++) {
        int col_var = ext_vars[cols[j]];
        if (col_var >= lower && col_var < upper) {
          A_row_vars[A_row_size] = col_var;
          A_row_size++;
        } else {
          B_row_vars[B_row_size] = col_var;
          B_row_size++;
        }
      }
    }

    // Merge the off-processor contributions to the rows of A/B
    for (int k = 0; k < num_in_procs; k++) {
      // Try to find the variable in the k-th input - these are sorted
      // locally between in_rows[in_row_ptr[k]:in_row_ptr[k+1]]
      int count = in_row_ptr[k + 1] - in_row_ptr[k];
      int *item = TacsSearchArray(var, count, &in_rows[in_row_ptr[k]]);

      if (item) {
        int row = item - &in_rows[in_row_ptr[k]];
        row += in_row_ptr[k];

        // Add variables in this range to the row as well
        int start = in_rowp[row];
        int end = in_rowp[row + 1];

        if (A_row_size + end - start > A_max_row_size) {
          A_max_row_size = A_max_row_size + end - start;
          TacsExtendArray(&A_row_vars, A_row_size, A_max_row_size);
        }

        if (B_row_size + end - start > B_max_row_size) {
          B_max_row_size = B_max_row_size + end - start;
          TacsExtendArray(&B_row_vars, B_row_size, B_max_row_size);
        }

        for (int j = start; j < end; j++) {
          int col_var = in_cols[j];
          if (col_var >= lower && col_var < upper) {
            A_row_vars[A_row_size] = col_var;
            A_row_size++;
          } else {
            B_row_vars[B_row_size] = col_var;
            B_row_size++;
          }
        }
      }
    }

    // Sort the entries and remove duplicates
    A_row_size = TacsUniqueSort(A_row_size, A_row_vars);
    A_rowp[var - lower + 1] = A_row_size;

    B_row_size = TacsUniqueSort(B_row_size, B_row_vars);
    B_rowp[var - lower + 1] = B_row_size;
  }

  // Now, set up A_rowp/B_rowp
  int np = -1;
  A_rowp[0] = 0;
  B_rowp[0] = 0;

  for (int i = 0; i < A_rows; i++) {
    A_rowp[i + 1] = A_rowp[i + 1] + A_rowp[i];
    B_rowp[i + 1] = B_rowp[i + 1] + B_rowp[i];
    if (B_rowp[i + 1] > 0 && np == -1) {
      np = i;
    }
  }
  if (np == -1) {
    np = A_rows;
  }

  int nc = A_rows - np;
  if (np > 0) {
    int *temp = new int[nc + 1];
    for (int i = 0; i < nc + 1; i++) {
      temp[i] = B_rowp[i + np];
    }
    delete[] B_rowp;
    B_rowp = temp;
  }

  // Now, go through and build up A_cols/B_cols
  int *A_cols = new int[A_rowp[A_rows]];
  int *B_cols = new int[B_rowp[nc]];

  // Size the A_rowp/B_rowp arrays
  for (int i = 0; i < A_rows; i++) {  // For each variable
    int A_row_size = 0;
    int B_row_size = 0;

    int ei = var_map[i];  // Index into the external variable map
    int var = lower + i;

    if (ei >= 0) {
      // Add variables in this range to the row as well
      int start = rowp[ei];
      int end = rowp[ei + 1];

      for (int j = start; j < end; j++) {
        int col_var = ext_vars[cols[j]];
        if (col_var >= lower && col_var < upper) {
          A_row_vars[A_row_size] = col_var - lower;
          A_row_size++;
        } else {
          B_row_vars[B_row_size] = col_var;
          B_row_size++;
        }
      }
    }

    // Merge the off-processor contributions to the rows of A/B
    for (int k = 0; k < num_in_procs; k++) {
      // Try to find the variable in the k-th input - these are sorted
      // locally between in_rows[in_row_ptr[k]:in_row_ptr[k+1]]
      int count = in_row_ptr[k + 1] - in_row_ptr[k];
      int *item = TacsSearchArray(var, count, &in_rows[in_row_ptr[k]]);

      if (item) {
        int row = item - &in_rows[in_row_ptr[k]];
        row += in_row_ptr[k];

        // Add variables in this range to the row as well
        int start = in_rowp[row];
        int end = in_rowp[row + 1];

        for (int j = start; j < end; j++) {
          int col_var = in_cols[j];
          if (col_var >= lower && col_var < upper) {
            A_row_vars[A_row_size] = col_var - lower;
            A_row_size++;
          } else {
            B_row_vars[B_row_size] = col_var;
            B_row_size++;
          }
        }
      }
    }

    // Sort the entries and remove duplicates
    A_row_size = TacsUniqueSort(A_row_size, A_row_vars);
    for (int j = A_rowp[var - lower], k = 0; k < A_row_size; j++, k++) {
      A_cols[j] = A_row_vars[k];
    }

    // Convert the global indices into the local ordering
    B_row_size = TacsUniqueSort(B_row_size, B_row_vars);
    if (var - lower >= np) {
      for (int k = 0; k < B_row_size; k++) {
        int *item = TacsSearchArray(B_row_vars[k], col_map_size, col_map_vars);

        if (!item) {
          fprintf(stderr, "[%d] Error: variable %d not in column map\n",
                  mpiRank, B_row_vars[k]);
        } else {
          B_row_vars[k] = item - col_map_vars;
        }
      }

      int index = var - lower - np;
      for (int j = B_rowp[index], k = 0; k < B_row_size; j++, k++) {
        B_cols[j] = B_row_vars[k];
      }
    }
  }

  delete[] A_row_vars;
  delete[] B_row_vars;
  delete[] var_map;

  *_Arowp = A_rowp;
  *_Acols = A_cols;
  *_np = np;
  *_Browp = B_rowp;
  *_Bcols = B_cols;
}

/*!
  Add values to the matrix.

  If the indices are owned on this processor add them here.  If not,
  add them to arrays that will be passed to other processors during
  the final assembly process.

  If the row or column index passed in is negative, then the entry is
  skipped.

  input:
  nrow:     number of block rows in the dense matrix
  row:      block row numbers in the global matrix
  ncol:     number of block columns in the dense matrix
  col:      block column numbers in the global matrix
  nv, mv:   number of true rows/columns in the dense matrix
  values:   the dense matrix values
*/
void TACSMatDistribute::addValues(TACSParallelMat *mat, int nrow,
                                  const int *row, int ncol, const int *col,
                                  int nv, int mv, const TacsScalar *values) {
  // Get the block matrices
  BCSRMat *Aloc, *Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  // Get the number of local variables and number of coupling
  // variables
  int N, Nc;
  mat->getRowMap(NULL, &N, &Nc);
  int Np = N - Nc;

  int bsize = Aloc->getBlockSize();
  int b2 = bsize * bsize;

  // Set up storage for the variable numbers
  int array[256];
  int *temp = NULL, *acols = NULL, *bcols = NULL;

  // If the number of columns is greater than 128,
  // allocate a new array that can store that many
  if (ncol > 128) {
    temp = new int[2 * ncol];
    acols = &temp[0];
    bcols = &temp[ncol];
  } else {
    acols = &array[0];
    bcols = &array[ncol];
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  const int *ownerRange;
  row_map->getOwnerRange(&ownerRange);

  // The lower/upper variable ranges
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank + 1];

  // Determine what is the diagonal block and what is off-diagonal
  int nb = 0;
  for (int i = 0; i < ncol; i++) {
    int c = col[i];
    acols[i] = -1;
    bcols[i] = -1;

    if (c >= lower && c < upper) {
      // The column is in A
      acols[i] = c - lower;
    } else if (c >= 0) {
      // The column is in the B off-processor part
      int *item = TacsSearchArray(c, col_map_size, col_map_vars);
      if (item) {
        nb++;
        bcols[i] = item - col_map_vars;
      }
    }

    // If we were supposed to find something and didn't,
    // print an error
    if (c >= 0 && acols[i] == -1 && bcols[i] == -1) {
      fprintf(stderr,
              "[%d] TACSMatDistribute: Could not find "
              "a match for column %d\n",
              mpiRank, c);
    }
  }

  for (int i = 0; i < nrow; i++) {
    int r = row[i];

    // Check if the row is on this processor
    if (r >= lower && r < upper) {
      r = r - lower;
      // Add values to the diagonal
      Aloc->addRowValues(r, ncol, acols, mv, &values[mv * i * bsize]);

      // Add values to the off-diagonal
      r = r - Np;
      if (r >= 0 && r < Nc) {
        Bext->addRowValues(r, ncol, bcols, mv, &values[mv * i * bsize]);
      } else if (nb > 0) {
        fprintf(stderr,
                "[%d] TACSMatDistribute error: some values "
                "were not added\n",
                mpiRank);
      }
    } else if (r >= 0) {
      // The row is not on this processor, search for it in the
      // off-processor part that we'll have to communicate later
      // Locate the row within ext_rows
      int *item = TacsSearchArray(r, num_ext_rows, ext_rows);
      if (item) {
        int r_ext = item - ext_rows;

        // Find the values within ext_cols
        for (int j = 0; j < ncol; j++) {
          int c = col[j];

          if (c >= 0) {
            int start = ext_rowp[r_ext];
            int size = ext_rowp[r_ext + 1] - start;
            item = TacsSearchArray(c, size, &ext_cols[start]);

            if (item) {
              TacsScalar *a = &ext_A[b2 * (item - ext_cols)];

              for (int ii = 0; ii < bsize; ii++) {
                for (int jj = 0; jj < bsize; jj++) {
                  a[ii * bsize + jj] +=
                      values[mv * (ii + i * bsize) + (jj + j * bsize)];
                }
              }
            } else {
              fprintf(stderr,
                      "[%d] TACSMatDistribute error: could not "
                      "find col (%d,%d) r_ext = %d\n",
                      mpiRank, r, c, r_ext);
            }
          }
        }
      } else {
        fprintf(stderr,
                "[%d] TACSMatDistribute error: could not "
                "find row %d\n",
                mpiRank, r);
      }
    }
  }
  if (temp) {
    delete[] temp;
  }
}

/*
  Add a weighted sum of the dense input matrix.

  This function adds an inner product of a weighting matrix with a
  dense matrix to the DistMat matrix. The weight matrix is a sparse,
  low-dimensional matrix given in a CSR-type format. The code takes
  this representation of W and adds the terms:

  self <- self + W^{T}*Kv*W

  to the values in the matrix. This code can be used to add the effect
  of dependent nodes to the matrix.

  input:
  nvars:    the number of block variables in the dense input matrix
  varp:     pointer to the start of each weight matrix; len(varp) = nvars+1
  vars:     variable numbers in the matrix; len(vars) = varp[nvars]
  weights:  weighting values for each variable; len(weights) = len(vars)
  nv, nw:   the number of rows and number of columns in the values matrix
  values:   the dense input matrix
*/
void TACSMatDistribute::addWeightValues(TACSParallelMat *mat, int nvars,
                                        const int *varp, const int *vars,
                                        const TacsScalar *weights, int nv,
                                        int mv, const TacsScalar *values,
                                        MatrixOrientation matOr) {
  // Get the block matrices
  BCSRMat *Aloc, *Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  // Get the number of local variables and number of coupling
  // variables
  int N, Nc;
  mat->getRowMap(NULL, &N, &Nc);
  int Np = N - Nc;

  int bsize = Aloc->getBlockSize();
  int b2 = bsize * bsize;

  // The number of variables we'll have to convert = row dim(W^{T})
  int n = varp[nvars];

  // Set up storage for the variable numbers
  int array[256];
  int *temp = NULL, *avars = NULL, *bvars = NULL;

  // If the number of columns is greater than 128,
  // allocate a new array that can store that many
  if (n > 128) {
    temp = new int[2 * n];
    avars = &temp[0];
    bvars = &temp[n];
  } else {
    avars = &array[0];
    bvars = &array[n];
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  const int *ownerRange;
  row_map->getOwnerRange(&ownerRange);

  // The lower/upper variable ranges
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank + 1];

  // Determine what is the diagonal block and what is off-diagonal
  int nb = 0;
  for (int i = 0; i < n; i++) {
    int c = vars[i];
    avars[i] = -1;
    bvars[i] = -1;

    if (c >= lower && c < upper) {
      // The column is in A
      avars[i] = c - lower;
    } else if (c >= 0) {
      // The column is in the B off-processor part
      int *item = TacsSearchArray(c, col_map_size, col_map_vars);
      if (item) {
        nb++;
        bvars[i] = item - col_map_vars;
      }
    }

    // If we were supposed to find something and didn't,
    // print an error
    if (c >= 0 && avars[i] == -1 && bvars[i] == -1) {
      fprintf(stderr,
              "[%d] TACSMatDistribute error: Could not find "
              "a match for column %d\n",
              mpiRank, c);
    }
  }

  // Set the increment along the row or column of the matrix depending
  // on whether we are adding the original matrix or its transpose
  int incr = mv;
  if (matOr == TACS_MAT_TRANSPOSE) {
    incr = 1;
  }

  for (int i = 0; i < nvars; i++) {
    for (int ip = varp[i]; ip < varp[i + 1]; ip++) {
      // Check if the row is on this processor
      if (avars[ip] >= 0) {
        // Add values to the diagonal
        Aloc->addRowWeightValues(weights[ip], avars[ip], nvars, varp, avars,
                                 weights, mv, &values[incr * i * bsize], matOr);

        // Add values to the off-diagonal
        int r = avars[ip] - Np;
        if (r >= 0 && r < Nc) {
          Bext->addRowWeightValues(weights[ip], r, nvars, varp, bvars, weights,
                                   mv, &values[incr * i * bsize], matOr);
        } else if (nb > 0) {
          fprintf(stderr,
                  "[%d] TACSMatDistribute error, some "
                  "values were not added\n",
                  mpiRank);
        }
      } else if (bvars[ip] >= 0) {
        // The row is not on this processor, search for it in the
        // off-processor part that we'll have to communicate later
        // Locate the row within ext_rows
        int *item = TacsSearchArray(vars[ip], num_ext_rows, ext_rows);
        if (item) {
          int r_ext = item - ext_rows;

          // Find the values within ext_cols
          for (int j = 0; j < nvars; j++) {
            for (int jp = varp[j]; jp < varp[j + 1]; jp++) {
              int c = vars[jp];
              TacsScalar aw = weights[ip] * weights[jp];

              if (c >= 0 && aw != 0.0) {
                int start = ext_rowp[r_ext];
                int size = ext_rowp[r_ext + 1] - start;
                item = TacsSearchArray(c, size, &ext_cols[start]);

                if (item) {
                  TacsScalar *a = &ext_A[b2 * (item - ext_cols)];

                  if (matOr == TACS_MAT_NORMAL) {
                    for (int ii = 0; ii < bsize; ii++) {
                      for (int jj = 0; jj < bsize; jj++) {
                        a[ii * bsize + jj] +=
                            aw *
                            values[mv * (ii + i * bsize) + (jj + j * bsize)];
                      }
                    }
                  } else {
                    for (int ii = 0; ii < bsize; ii++) {
                      for (int jj = 0; jj < bsize; jj++) {
                        a[ii * bsize + jj] +=
                            aw *
                            values[mv * (jj + j * bsize) + (ii + i * bsize)];
                      }
                    }
                  }
                } else {
                  fprintf(stderr,
                          "[%d] TACSMatDistribute error: could not "
                          "find col (%d,%d) r_ext = %d\n",
                          mpiRank, vars[jp], c, r_ext);
                }
              }
            }
          }
        } else {
          fprintf(stderr,
                  "[%d] TACSMatDistribute error: could not "
                  "find row %d\n",
                  mpiRank, vars[ip]);
        }
      }
    }
  }

  if (temp) {
    delete[] temp;
  }
}

/*!
  Given a non-zero pattern, pass in the values for the array
*/
void TACSMatDistribute::setValues(TACSParallelMat *mat, int nvars,
                                  const int *ext_vars, const int *rowp,
                                  const int *cols, TacsScalar *avals) {
  // Get the block matrices
  BCSRMat *Aloc, *Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  // Get the number of local variables and number of coupling
  // variables
  int N, Nc;
  mat->getRowMap(NULL, &N, &Nc);
  int Np = N - Nc;

  int bsize = Aloc->getBlockSize();
  int b2 = bsize * bsize;

  // Determine the maximum size of the array
  int max_row = 0;
  for (int i = 0; i < nvars; i++) {
    int size = rowp[i + 1] - rowp[i];
    if (size > max_row) {
      max_row = size;
    }
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  const int *ownerRange;
  row_map->getOwnerRange(&ownerRange);

  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank + 1];

  int *acols = new int[2 * max_row];
  int *bcols = &acols[max_row];

  for (int i = 0; i < nvars; i++) {
    int row = ext_vars[i];

    // Check if this row is in the local or external block
    int nb = 0;
    if (row >= lower && row < upper) {
      row = row - lower;

      // Convert the cols
      int start = rowp[i];
      int end = rowp[i + 1];
      for (int j = rowp[i], k = 0; j < end; j++, k++) {
        int col = ext_vars[cols[j]];
        acols[k] = -1;
        bcols[k] = -1;

        if (col >= lower && col < upper) {
          acols[k] = col - lower;
        } else {
          int *item = TacsSearchArray(col, col_map_size, col_map_vars);
          bcols[k] = item - col_map_vars;
          nb++;
        }
      }

      Aloc->addBlockRowValues(row, end - start, acols, &avals[b2 * start]);

      if (nb > 0) {
        row = row - Np;
        if (row >= 0 && row < Nc) {
          Bext->addBlockRowValues(row, end - start, bcols, &avals[b2 * start]);
        } else {
          fprintf(stderr, "[%d] DistMat error: could not find row %d\n",
                  mpiRank, row);
        }
      }
    } else {
      int *item = TacsSearchArray(row, num_ext_rows, ext_rows);

      if (item) {
        int r_ext = item - ext_rows;

        int end = rowp[i + 1];
        for (int j = rowp[i], k = 0; j < end; j++, k++) {
          int c = cols[j];
          if (c >= 0 && c < nvars) {
            int col = ext_vars[c];

            int ext_start = ext_rowp[r_ext];
            int ext_size = ext_rowp[r_ext + 1] - ext_start;
            item = TacsSearchArray(col, ext_size, &ext_cols[ext_start]);

            if (item) {
              TacsScalar *a = &ext_A[b2 * (item - ext_cols)];
              memcpy(a, &avals[b2 * j], b2 * sizeof(TacsScalar));
            } else {
              fprintf(stderr,
                      "[%d] DistMat error: could not find col "
                      "(%d,%d) r_ext = %d \n",
                      mpiRank, row, col, r_ext);
            }
          } else {
            fprintf(stderr,
                    "[%d] DistMat error: local column out of "
                    "range 0 <= %d < %d\n",
                    mpiRank, c, nvars);
          }
        }
      } else {
        fprintf(stderr, "[%d] DistMat error: could not find row %d\n", mpiRank,
                row);
      }
    }
  }

  delete[] acols;
}

/*
   Initiate the communication of the off-process matrix entries
*/
void TACSMatDistribute::beginAssembly(TACSParallelMat *mat) {
  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  // Get the block size squared
  const int b2 = bsize * bsize;

  // Post the recvs
  for (int k = 0, offset = 0, buff_offset = 0; k < num_in_procs; k++) {
    int count = in_rowp[offset + in_count[k]] - in_rowp[offset];
    count *= b2;

    int source = in_procs[k];
    int tag = 5;
    MPI_Irecv(&in_A[buff_offset], count, TACS_MPI_TYPE, source, tag, comm,
              &in_requests[k]);
    offset += in_count[k];
    buff_offset += count;
  }

  // Post the sends
  for (int k = 0, offset = 0, buff_offset = 0; k < num_ext_procs; k++) {
    int count = ext_rowp[offset + ext_count[k]] - ext_rowp[offset];
    count *= b2;

    int dest = ext_procs[k];
    int tag = 5;
    MPI_Isend(&ext_A[buff_offset], count, TACS_MPI_TYPE, dest, tag, comm,
              &ext_requests[k]);
    offset += ext_count[k];
    buff_offset += count;
  }
}

/*
  Finish the communication of the off-process matrix entries.
  Once the communication is completed, add the off-processor
  entries to the matrix.
*/
void TACSMatDistribute::endAssembly(TACSParallelMat *mat) {
  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  // Get the map between the global-external variables and the local
  // variables (for Bext)
  BCSRMat *Aloc, *Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  // Get the number of local variables and number of coupling
  // variables
  int N, Nc;
  mat->getRowMap(NULL, &N, &Nc);
  int Np = N - Nc;

  // Get the block size squared
  const int b2 = bsize * bsize;

  // Find the owner range for mapping variables
  const int *ownerRange;
  row_map->getOwnerRange(&ownerRange);
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank + 1];

  for (int i = 0; i < num_in_procs; i++) {
    // Get the recv that just completed
    int index;
    MPI_Status status;
    int ierr = MPI_Waitany(num_in_procs, in_requests, &index, &status);

    // Check whether the recv was successful
    if (ierr != MPI_SUCCESS) {
      int len;
      char err_str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(ierr, err_str, &len);
      fprintf(stderr,
              "[%d] TACSMatDistribute::endAssembly MPI_Waitany "
              "error\n%s\n",
              mpiRank, err_str);
    }

    // Identify which group of rows were just recv'd from another
    // processor
    for (int j = in_row_ptr[index]; j < in_row_ptr[index + 1]; j++) {
      // Find the local index of the row
      int row = in_rows[j] - lower;

      // Find the local row index
      for (int k = in_rowp[j]; k < in_rowp[j + 1]; k++) {
        TacsScalar *a = &in_A[b2 * k];

        // Set the column indices
        int col = in_cols[k];
        if (col >= lower && col < upper) {
          // Get the local column index
          col = col - lower;
          Aloc->addBlockRowValues(row, 1, &col, a);
        } else {
          // Use the map from the global column index back to the
          // processor
          int *item = TacsSearchArray(col, col_map_size, col_map_vars);
          if (item) {
            int c = item - col_map_vars;
            Bext->addBlockRowValues(row - Np, 1, &c, a);
          }
        }
      }
    }
  }

  // Wait for all the sending requests
  if (num_ext_procs > 0) {
    MPI_Waitall(num_ext_procs, ext_requests, MPI_STATUSES_IGNORE);
  }
}
