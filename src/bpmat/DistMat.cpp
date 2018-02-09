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

#include "DistMat.h"
#include "FElibrary.h"
#include "MatUtils.h"

/*
  TACSDistMat implementation
*/

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
  thread_info:   object that stores pthreads info
  rmap:          variable mapping assignment to processors
  num_ext_vars:  the size of the block matrix on this processor
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
TACSDistMat::TACSDistMat( TACSThreadInfo *thread_info, 
                          TACSVarMap *map, int bsize, int num_ext_vars, 
                          const int *rowp, const int *cols, 
                          TACSBVecIndices *bindex ){
  comm = map->getMPIComm();

  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);

  // Get the external variables
  const int *ext_vars;
  if (bindex->getIndices(&ext_vars) != num_ext_vars){
    fprintf(stderr, "[%d] DistMat error: number of indices provided must be \
equal to the number of rows\n", mpiRank);
    return;
  }
  
  const int *ownerRange;
  map->getOwnerRange(&ownerRange);

  // Get the non-zero pattern of the input matrix
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank+1];

  // Count up the number of equations that must be sent to other
  // processors. 
  // Get only the variables associated with those ones.
  next_rows = 0;
  for ( int k = 0; k < num_ext_vars; k++ ){
    if (ext_vars[k] < lower || ext_vars[k] >= upper){
      next_rows++;
    }
  }

  ext_rows = new int[ next_rows ];
  next_rows = 0;
  for ( int k = 0; k < num_ext_vars; k++ ){
    if (ext_vars[k] < lower || ext_vars[k] >= upper){
      ext_rows[next_rows] = ext_vars[k];
      next_rows++;
    }
  }
  
  int next_row_init = next_rows;
  next_rows = FElibrary::uniqueSort(ext_rows, next_rows);
  if (next_row_init != next_rows){
    fprintf(stderr, "[%d] DistMat error, ext_vars are not unique\n", 
            mpiRank);
  }

  // Create the off-process CSR data structure that will be sent to 
  // other processes. First, calculate ext_rowp
  ext_rowp = new int[ next_rows+1 ];
  memset(ext_rowp, 0, (next_rows+1)*sizeof(int));

  for ( int i = 0; i < num_ext_vars; i++ ){
    if (ext_vars[i] < lower || ext_vars[i] >= upper){
      int *item = (int*)bsearch(&ext_vars[i], ext_rows, next_rows, 
                                sizeof(int), FElibrary::comparator);
      int index = item - ext_rows;
      ext_rowp[index+1] = rowp[i+1] - rowp[i];    
    }
  }

  for ( int i = 0; i < next_rows; i++ ){
    ext_rowp[i+1] += ext_rowp[i];
  }

  // Next, calculate ext_cols. Find only the external rows and 
  // convert them to the global numbering scheme
  ext_cols = new int[ ext_rowp[next_rows] ];

  for ( int i = 0; i < num_ext_vars; i++ ){
    if (ext_vars[i] < lower || ext_vars[i] >= upper){
      int *item = (int*)bsearch(&ext_vars[i], ext_rows, next_rows, 
                                sizeof(int), FElibrary::comparator);     
      int index = item - ext_rows;

      for ( int j = ext_rowp[index], jj = rowp[i]; 
            j < ext_rowp[index+1]; j++, jj++ ){
        ext_cols[j] = ext_vars[cols[jj]];
      }

      int size = ext_rowp[index+1] - ext_rowp[index];
      if (size != FElibrary::uniqueSort(&ext_cols[ext_rowp[index]], size)){
        fprintf(stderr, "[%d] DistMat error, array is not unique\n", 
                mpiRank);
      }
    }
  }
  
  // Match the intervals of the external variables to be 
  // sent to other processes
  ext_row_ptr = new int[ mpiSize+1 ];
  ext_row_count = new int[ mpiSize ];

  FElibrary::matchIntervals(mpiSize, ownerRange, next_rows, 
    ext_rows, ext_row_ptr);

  for ( int k = 0; k < mpiSize; k++ ){
    ext_row_count[k] = ext_row_ptr[k+1] - ext_row_ptr[k];    
  }

  in_row_count = new int[ mpiSize ];
  in_row_ptr= new int[ mpiSize+1 ];

  MPI_Alltoall(ext_row_count, 1, MPI_INT, in_row_count, 1, MPI_INT, comm);

  // Prepare to receive the equation numbers from the other processes
  in_row_ptr[0] = 0;
  for ( int k = 0; k < mpiSize; k++ ){
    in_row_ptr[k+1] = in_row_ptr[k] + in_row_count[k];
  }
  
  in_rows = new int[ in_row_ptr[mpiSize] ];

  MPI_Alltoallv(ext_rows, ext_row_count, ext_row_ptr, MPI_INT, 
                in_rows, in_row_count, in_row_ptr, MPI_INT, comm);
  
  // Prepare to pass information from ext_rowp to in_rowp
  int *ext_row_var_count = new int[ next_rows ];
  for ( int k = 0; k < next_rows; k++ ){
    ext_row_var_count[k] = ext_rowp[k+1] - ext_rowp[k];
  }

  nin_rows = in_row_ptr[mpiSize];
  in_rowp = new int[ nin_rows+1 ];
  in_rowp[0] = 0;

  MPI_Alltoallv(ext_row_var_count, ext_row_count, ext_row_ptr, MPI_INT,
                &in_rowp[1], in_row_count, in_row_ptr, MPI_INT, comm);

  for ( int k = 0; k < nin_rows; k++ ){
    in_rowp[k+1] += in_rowp[k];
  }

  delete [] ext_row_var_count;

  // Pass ext_cols to in_cols
  in_cols = new int[ in_rowp[nin_rows] ];

  int *ext_cols_ptr = new int[ mpiSize ];
  int *ext_cols_count = new int[ mpiSize ];
  int *in_cols_ptr = new int[ mpiSize ];
  int *in_cols_count = new int[ mpiSize ];

  for ( int i = 0; i < mpiSize; i++ ){
    ext_cols_ptr[i] = ext_rowp[ ext_row_ptr[i] ];    
    ext_cols_count[i] = 
      ext_rowp[ ext_row_ptr[i+1] ] - ext_rowp[ ext_row_ptr[i] ];
    in_cols_ptr[i] = in_rowp[ in_row_ptr[i] ];
    in_cols_count[i] = 
      in_rowp[ in_row_ptr[i+1] ] - in_rowp[ in_row_ptr[i] ]; 
  }

  // Send the column numbers to the other processors
  MPI_Alltoallv(ext_cols, ext_cols_count, ext_cols_ptr, MPI_INT, 
                in_cols, in_cols_count, in_cols_ptr, MPI_INT, comm);

  delete [] ext_cols_count;
  delete [] ext_cols_ptr;
  delete [] in_cols_count;
  delete [] in_cols_ptr;

  // estimate the non-zero entries per row
  int nz_per_row = 10;
  if (num_ext_vars > 0){
    nz_per_row = (int)(rowp[num_ext_vars]/num_ext_vars + 1);
  }
  
  // Now, set up the column map using in_cols  
  int max_col_vars_size = 5*nz_per_row + next_rows;
  int col_vars_size = 0;
  int *col_vars = new int[ max_col_vars_size ];

  // Get contributions from ext_rows
  for ( int i = 0; i < next_rows; i++ ){
    col_vars[col_vars_size] = ext_rows[i];
    col_vars_size++;
  }

  for ( int i = 0; i < nin_rows; i++ ){
    if (col_vars_size + in_rowp[i+1] - in_rowp[i] > max_col_vars_size){
      max_col_vars_size = 2.0*max_col_vars_size;
      matutils::ExtendArray(&col_vars, col_vars_size, max_col_vars_size);
    }    

    for ( int j = in_rowp[i]; j < in_rowp[i+1]; j++ ){
      if (in_cols[j] <  ownerRange[mpiRank] ||
          in_cols[j] >= ownerRange[mpiRank+1]){
        col_vars[col_vars_size] = in_cols[j];
        col_vars_size++;
      }
    }

    col_vars_size = FElibrary::uniqueSort(col_vars, col_vars_size);
  }

  TACSBVecIndices *col_indices = new TACSBVecIndices(&col_vars, col_vars_size);
  TACSBVecDistribute *col_map = new TACSBVecDistribute(map, col_indices);
  col_map_size = col_indices->getIndices(&col_map_vars);

  // Assemble everything into on and off-diagonal parts
  int *Arowp, *Acols; // The diagonal entries
  int np, *Browp, *Bcols;

  setUpLocalExtCSR(num_ext_vars, ext_vars, rowp, cols,
                   ownerRange[mpiRank], ownerRange[mpiRank+1],
                   nz_per_row, &Arowp, &Acols,
                   &np, &Browp, &Bcols);

  // Allocate space for in-coming matrix elements
  ext_A = new TacsScalar[ bsize*bsize*ext_rowp[next_rows] ];
  in_A = new TacsScalar[ bsize*bsize*in_rowp[nin_rows] ];

  int len = bsize*bsize*ext_rowp[next_rows];
  for ( int k = 0; k < len; k++ ){
    ext_A[k] = 0.0;
  }

  // Allocate the local/external matrices  
  int n = ownerRange[mpiRank+1] - ownerRange[mpiRank];
  int m = n;

  BCSRMat *aloc = new BCSRMat(comm, thread_info, 
                              bsize, n, m, &Arowp, &Acols);
  BCSRMat *bext = new BCSRMat(comm, thread_info,
                              bsize, n-np, col_vars_size, &Browp, &Bcols);
  
  // Finally, initialize PMat
  init(map, aloc, bext, col_map);

  // Set up the presistent communication amongst processors
  initPersistent();
}

TACSDistMat::~TACSDistMat(){
  delete [] ext_rows;
  delete [] ext_rowp;
  delete [] ext_cols;
  delete [] ext_row_ptr;
  delete [] ext_row_count;
  
  delete [] in_rows;
  delete [] in_rowp;
  delete [] in_cols;
  delete [] in_row_ptr;
  delete [] in_row_count;

  delete [] in_A;
  delete [] ext_A;
  
  if (nsends > 0){
    delete [] sends;
    delete [] send_proc;
    delete [] send_status;
  }
  if (nreceives > 0){
    delete [] receives;
    delete [] receive_proc;
    delete [] receive_status;
  }
}

/*!
  Find the non-zero pattern for the A = [B, E; F, C] and Bext 
  matrices such that:

  [ B, E ][ x ]
  [ F, C ][ y ] + [ Bext ][ y_ext ] = 0  
*/
void TACSDistMat::setUpLocalExtCSR( int num_ext_vars, const int *ext_vars, 
                                    const int *rowp, const int *cols,
                                    int lower, int upper,
                                    int nz_per_row,
                                    int **_Arowp, int **_Acols,
                                    int *_np, int **_Browp, int **_Bcols ){
  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);

  // For each row, find the non-zero pattern
  int A_max_row_size = 2 *nz_per_row;
  int *A_row_vars = new int[ A_max_row_size ];

  int B_max_row_size = 2 *nz_per_row;
  int *B_row_vars = new int[ B_max_row_size ];

  int A_rows = upper - lower;
  int *A_rowp = new int[ A_rows+1 ];
  int *B_rowp = new int[ A_rows+1 ];

  for ( int i = 0; i < A_rows+1; i++ ){
    A_rowp[i] = 0;
    B_rowp[i] = 0;
  }

  // Set up a temporary array to keep track of the variables in ext_vars
  int *var_map = new int[ A_rows ];
  for ( int i = 0; i < A_rows; i++ ){
    var_map[i] = -1;
  }
  
  for ( int i = 0; i < num_ext_vars; i++ ){
    int var = ext_vars[i];
    if (var >= lower && var < upper){
      var_map[var - lower] = i;
    }
  }

  // Size the A_rowp/B_rowp arrays
  for ( int i = 0; i < A_rows; i++ ){ // For each variable
    int A_row_size = 0;
    int B_row_size = 0;

    int ei  = var_map[i]; // Index into the external variable map
    int var = lower + i;

    if (ei >= 0){
      // Add variables in this range to the row as well
      int start = rowp[ei];
      int end   = rowp[ei+1];
      
      // Add rowp[ei], to A_row_vars, B_row_vars
      if (A_row_size + end-start > A_max_row_size){
        A_max_row_size = A_max_row_size + end-start;
        matutils::ExtendArray(&A_row_vars, A_row_size, A_max_row_size);
      }
      
      if (B_row_size + end-start > B_max_row_size){
        B_max_row_size = B_max_row_size + end-start;
        matutils::ExtendArray(&B_row_vars, B_row_size, B_max_row_size);
      }

      for ( int j = start; j < end; j++ ){
        int col_var = ext_vars[cols[j]];
        if (col_var >= lower && col_var < upper){
          A_row_vars[A_row_size] = col_var;
          A_row_size++;
        }
        else {
          B_row_vars[B_row_size] = col_var;
          B_row_size++;
        }
      }
    }

    // Merge the off-processor contributions to the rows of A/B
    for ( int k = 0; k < mpiSize; k++ ){
      if (in_row_count[k] > 0){
        // Try to find the variable in the k-th input - these are sorted 
        // locally between in_rows[in_row_ptr[k]:in_row_ptr[k+1]]
        int *item = (int*)bsearch(&var, &in_rows[in_row_ptr[k]], 
                                  in_row_count[k], sizeof(int), 
                                  FElibrary::comparator);

        if (item){
          int row = item - &in_rows[in_row_ptr[k]];
          row += in_row_ptr[k];

          // Add variables in this range to the row as well
          int start = in_rowp[row];
          int end   = in_rowp[row+1];

          if (A_row_size + end-start > A_max_row_size){
            A_max_row_size = A_max_row_size + end-start;
            matutils::ExtendArray(&A_row_vars, A_row_size, A_max_row_size);
          }

          if (B_row_size + end-start > B_max_row_size){
            B_max_row_size = B_max_row_size + end-start;
            matutils::ExtendArray(&B_row_vars, B_row_size, B_max_row_size);
          }
          
          for ( int j = start; j < end; j++ ){
            int col_var = in_cols[j];
            if (col_var >= lower && col_var < upper){
              A_row_vars[A_row_size] = col_var;
              A_row_size++;
            }
            else {
              B_row_vars[B_row_size] = col_var;
              B_row_size++;
            }
          }
        }
      }
    }

    // Sort the entries and remove duplicates
    A_row_size = FElibrary::uniqueSort(A_row_vars, A_row_size);
    A_rowp[var - lower+1] = A_row_size;

    B_row_size = FElibrary::uniqueSort(B_row_vars, B_row_size);
    B_rowp[var - lower+1] = B_row_size;
  }

  // Now, set up A_rowp/B_rowp
  int np = -1;
  A_rowp[0] = 0;
  B_rowp[0] = 0;

  for ( int i = 0; i < A_rows; i++ ){
    A_rowp[i+1] = A_rowp[i+1] + A_rowp[i];
    B_rowp[i+1] = B_rowp[i+1] + B_rowp[i];
    if (B_rowp[i+1] > 0 && np == -1){
      np = i;
    }
  }
  if (np == -1){ 
    np = A_rows; 
  }

  int nc = A_rows - np;
  if (np > 0){
    int *temp = new int[ nc+1 ];
    for ( int i = 0; i < nc+1; i++ ){
      temp[i] = B_rowp[i+np];
    }
    delete [] B_rowp;
    B_rowp = temp;
  }

  // Now, go through and build up A_cols/B_cols
  int *A_cols = new int[ A_rowp[A_rows] ];
  int *B_cols = new int[ B_rowp[nc] ];
  
  // Size the A_rowp/B_rowp arrays
  for ( int i = 0; i < A_rows; i++ ){ // For each variable
    int A_row_size = 0;
    int B_row_size = 0;

    int ei  = var_map[i]; // Index into the external variable map
    int var = lower + i;

    if (ei >= 0){
      // Add variables in this range to the row as well
      int start = rowp[ei];
      int end   = rowp[ei+1];

      for ( int j = start; j < end; j++ ){
        int col_var = ext_vars[cols[j]];
        if (col_var >= lower && col_var < upper){
          A_row_vars[A_row_size] = col_var - lower;
          A_row_size++;
        }
        else {
          B_row_vars[B_row_size] = col_var;
          B_row_size++;
        }
      }
    }

    // Merge the off-processor contributions to the rows of A/B
    for ( int k = 0; k < mpiSize; k++ ){
      if (in_row_count[k] > 0){
        // Try to find the variable in the k-th input - these are sorted 
        // locally between in_rows[in_row_ptr[k]:in_row_ptr[k+1]]
        int *item = (int*)bsearch(&var, &in_rows[in_row_ptr[k]], 
                                  in_row_count[k], sizeof(int), 
                                  FElibrary::comparator);

        if (item){
          int row = item - &in_rows[in_row_ptr[k]];
          row += in_row_ptr[k];

          // Add variables in this range to the row as well
          int start = in_rowp[row];
          int end   = in_rowp[row+1];
          
          for ( int j = start; j < end; j++ ){
            int col_var = in_cols[j];
            if (col_var >= lower && col_var < upper){
              A_row_vars[A_row_size] = col_var - lower;
              A_row_size++;
            }
            else {
              B_row_vars[B_row_size] = col_var;
              B_row_size++;
            }
          }
        }
      }
    }

    // Sort the entries and remove duplicates
    A_row_size = FElibrary::uniqueSort(A_row_vars, A_row_size);
    for ( int j = A_rowp[var - lower], k = 0; k < A_row_size; j++, k++ ){
      A_cols[j] = A_row_vars[k];      
    }

    // Convert the global indices into the local ordering
    B_row_size = FElibrary::uniqueSort(B_row_vars, B_row_size);
    if (var - lower >= np){
      for ( int k = 0; k < B_row_size; k++ ){
        int *item = (int*)bsearch(&B_row_vars[k], col_map_vars, 
                                  col_map_size, sizeof(int), 
                                  FElibrary::comparator);
        
        if (!item){
          fprintf(stderr, "[%d] Error: variable %d not in column map\n", 
                  mpiRank, B_row_vars[k]);
        }
        else {
          B_row_vars[k] = item - col_map_vars;
        }
      }

      int index = var - lower - np;
      for ( int j = B_rowp[index], k = 0; k < B_row_size; j++, k++ ){
        B_cols[j] = B_row_vars[k];      
      }
    }
  }

  delete [] A_row_vars;
  delete [] B_row_vars;
  delete [] var_map;

  *_Arowp = A_rowp;
  *_Acols = A_cols;
  *_np = np;
  *_Browp = B_rowp;
  *_Bcols = B_cols;
}

/*!  
  Initialize the persistent communication for the transfer of the
  matrix values that are set on this processor, but are destined for
  other processors.
*/
void TACSDistMat::initPersistent(){   
  int bsize = Aloc->getBlockSize();
  
  int mpiRank, mpiSize;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);

  // Count up the number of non-zero send/receives
  nreceives = 0;
  nsends = 0;

  for ( int i = 0; i < mpiSize; i++ ){
    if (ext_row_count[i] > 0){
      nsends++;
    }
    if (in_row_count[i] > 0){
      nreceives++;
    }
  }

  send_proc = NULL;
  receive_proc = NULL;
  send_status = NULL;
  receive_status = NULL;

  if (nsends > 0){
    sends = new MPI_Request[ nsends ];
    send_proc = new int[ nsends ];      
    send_status = new MPI_Status[ nsends ];      
  }
  else {
    sends = NULL;
    send_proc = NULL;
    send_status = NULL;
  }
  if (nreceives > 0){
    receives = new MPI_Request[ nreceives ];
    receive_proc = new int[ nreceives ]; 
    receive_status = new MPI_Status[ nreceives ];
  }
  else {
    receives = NULL;
    receive_proc = NULL;
    receive_status = NULL;
  }

  // Set up the sends
  int n = 0;
  int tag = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    if (ext_row_count[i] > 0){
      int var = ext_row_ptr[i];
      int count = bsize*bsize*(ext_rowp[ext_row_ptr[i+1]] - ext_rowp[var]);
      MPI_Send_init(&ext_A[bsize*bsize*ext_rowp[var]], count, TACS_MPI_TYPE, 
                    i, tag, comm, &sends[n]);
      send_proc[n] = i;
      n++;
    }
  }

  // Set up the receives
  n = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    if (in_row_count[i] > 0){
      int var = in_row_ptr[i];
      int count = bsize*bsize*(in_rowp[in_row_ptr[i+1]] - in_rowp[var]);
      MPI_Recv_init(&in_A[bsize*bsize*in_rowp[var]], count, TACS_MPI_TYPE, 
                    i, tag, comm, &receives[n]);
      receive_proc[n] = i;
      n++;
    }
  }
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
void TACSDistMat::addValues( int nrow, const int *row, 
                             int ncol, const int *col,
                             int nv, int mv, const TacsScalar *values ){
  int bsize = Aloc->getBlockSize();
  int b2 = bsize*bsize;

  // Set up storage for the variable numbers
  int array[256];
  int *temp = NULL, *acols = NULL, *bcols = NULL;
    
  // If the number of columns is greater than 128,
  // allocate a new array that can store that many  
  if (ncol > 128){
    temp = new int[ 2*ncol ];
    acols = &temp[0];
    bcols = &temp[ncol];
  }
  else {
    acols = &array[0];
    bcols = &array[ncol];
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);
  
  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  // The lower/upper variable ranges
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank+1];

  // Determine what is the diagonal block and what is off-diagonal
  int nb = 0;
  for ( int i = 0; i < ncol; i++ ){
    int c = col[i];
    acols[i] = -1;
    bcols[i] = -1;

    if (c >= lower && c < upper){
      // The column is in A
      acols[i] = c - lower;
    }
    else if (c >= 0){
      // The column is in the B off-processor part 
      int *item = (int*)bsearch(&c, col_map_vars, col_map_size, 
                                sizeof(int), FElibrary::comparator);
      if (item){
        nb++;
        bcols[i] = item - col_map_vars;
      }
    }

    // If we were supposed to find something and didn't,
    // print an error
    if (c >= 0 && 
        acols[i] == -1 &&
        bcols[i] == -1){
      fprintf(stderr, "[%d] Could not find a match for column %d\n",
              mpiRank, c);
    }
  }
  
  for ( int i = 0; i < nrow; i++ ){
    int r = row[i];

    // Check if the row is on this processor
    if (r >= lower && r < upper){
      r = r - lower;
      // Add values to the diagonal
      Aloc->addRowValues(r, ncol, acols, mv, &values[mv*i*bsize]);
      
      // Add values to the off-diagonal
      r = r - Np;
      if (r >= 0 && r < Nc){
        Bext->addRowValues(r, ncol, bcols, mv, &values[mv*i*bsize]);
      }
      else if (nb > 0){
        fprintf(stderr, "[%d] DistMat error, some values were not added\n",
                mpiRank);
      }
    }
    else if (r >= 0){
      // The row is not on this processor, search for it in the
      // off-processor part that we'll have to communicate later
      // Locate the row within ext_rows
      int *item = (int*)bsearch(&r, ext_rows, next_rows, 
                                sizeof(int), FElibrary::comparator);
      if (item){
        int r_ext = item - ext_rows;

        // Find the values within ext_cols
        for ( int j = 0; j < ncol; j++ ){
          int c = col[j];
          
          if (c >= 0){
            int start = ext_rowp[r_ext];
            int size = ext_rowp[r_ext+1] - start;
            item = (int*)bsearch(&c, &ext_cols[start], size, 
                                 sizeof(int), FElibrary::comparator);

            if (item){
              TacsScalar *a = &ext_A[b2*(item - ext_cols)];
              
              for ( int ii = 0; ii < bsize; ii++ ){
                for ( int jj = 0; jj < bsize; jj++ ){
                  a[ii*bsize + jj] += 
                    values[mv*(ii + i*bsize) + (jj + j*bsize)];
                }
              }
            }
            else {
              fprintf(stderr, "[%d] DistMat error: could not find col \
(%d,%d) r_ext = %d\n", mpiRank, r, c, r_ext);
            }
          }
        }
      }
      else {
        fprintf(stderr, "[%d] DistMat error: could not find row %d\n", 
                mpiRank, r);
      }
    }
  }
  if (temp){ delete [] temp; }
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
void TACSDistMat::addWeightValues( int nvars, const int *varp, const int *vars,
                                   const TacsScalar *weights,
                                   int nv, int mv, const TacsScalar *values,
                                   MatrixOrientation matOr ){
  int bsize = Aloc->getBlockSize();
  int b2 = bsize*bsize;

  // The number of variables we'll have to convert = row dim(W^{T})
  int n = varp[nvars];

  // Set up storage for the variable numbers
  int array[256];
  int *temp = NULL, *avars = NULL, *bvars = NULL;
    
  // If the number of columns is greater than 128,
  // allocate a new array that can store that many  
  if (n > 128){
    temp = new int[ 2*n ];
    avars = &temp[0];
    bvars = &temp[n];
  }
  else {
    avars = &array[0];
    bvars = &array[n];
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);
  
  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  // The lower/upper variable ranges
  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank+1];

  // Determine what is the diagonal block and what is off-diagonal
  int nb = 0;
  for ( int i = 0; i < n; i++ ){
    int c = vars[i];
    avars[i] = -1;
    bvars[i] = -1;

    if (c >= lower && c < upper){
      // The column is in A
      avars[i] = c - lower;
    }
    else if (c >= 0){
      // The column is in the B off-processor part 
      int *item = (int*)bsearch(&c, col_map_vars, col_map_size, 
                                sizeof(int), FElibrary::comparator);
      if (item){
        nb++;
        bvars[i] = item - col_map_vars;
      }
    }

    // If we were supposed to find something and didn't,
    // print an error
    if (c >= 0 && 
        avars[i] == -1 &&
        bvars[i] == -1){
      fprintf(stderr, "[%d] Could not find a match for column %d\n",
              mpiRank, c);
    }
  }
  
  // Set the increment along the row or column of the matrix depending
  // on whether we are adding the original matrix or its transpose
  int incr = mv;
  if (matOr == TRANSPOSE){
    incr = 1;
  }

  for ( int i = 0; i < nvars; i++ ){
    for ( int ip = varp[i]; ip < varp[i+1]; ip++ ){
      // Check if the row is on this processor
      if (avars[ip] >= 0){
        // Add values to the diagonal
        Aloc->addRowWeightValues(weights[ip], avars[ip], 
                                 nvars, varp, avars, weights, 
                                 mv, &values[incr*i*bsize], matOr);
        
        // Add values to the off-diagonal
        int r = avars[ip] - Np;
        if (r >= 0 && r < Nc){
          Bext->addRowWeightValues(weights[ip], r, 
                                   nvars, varp, bvars, weights, 
                                   mv, &values[incr*i*bsize], matOr);
        }
        else if (nb > 0){
          fprintf(stderr, "[%d] DistMat error, some values were not added\n",
                  mpiRank);
        }
      }
      else if (bvars[ip] >= 0){
        // The row is not on this processor, search for it in the
        // off-processor part that we'll have to communicate later
        // Locate the row within ext_rows
        int *item = (int*)bsearch(&vars[ip], ext_rows, next_rows, 
                                  sizeof(int), FElibrary::comparator);
        if (item){
          int r_ext = item - ext_rows;
          
          // Find the values within ext_cols
          for ( int j = 0; j < nvars; j++ ){
            for ( int jp = varp[j]; jp < varp[j+1]; jp++ ){
              int c = vars[jp];
              TacsScalar aw = weights[ip]*weights[jp];
            
              if (c >= 0 && aw != 0.0){
                int start = ext_rowp[r_ext];
                int size = ext_rowp[r_ext+1] - start;
                item = (int*)bsearch(&c, &ext_cols[start], size, 
                                     sizeof(int), FElibrary::comparator);
                
                if (item){
                  TacsScalar *a = &ext_A[b2*(item - ext_cols)];
                  
                  if (matOr == NORMAL){
                    for ( int ii = 0; ii < bsize; ii++ ){
                      for ( int jj = 0; jj < bsize; jj++ ){
                        a[ii*bsize + jj] += 
                          aw*values[mv*(ii + i*bsize) + (jj + j*bsize)];
                      }
                    }
                  }
                  else {
                    for ( int ii = 0; ii < bsize; ii++ ){
                      for ( int jj = 0; jj < bsize; jj++ ){
                        a[ii*bsize + jj] += 
                          aw*values[mv*(jj + j*bsize) + (ii + i*bsize)];
                      }
                    }
                  }
                }
                else {
                  fprintf(stderr, "[%d] DistMat error: could not find col \
(%d,%d) r_ext = %d\n", mpiRank, vars[jp], c, r_ext);
                }
              }
            }
          }
        }
        else {
          fprintf(stderr, "[%d] DistMat error: could not find row %d\n", 
                  mpiRank, vars[ip]);
        }
      }
    }
  }

  if (temp){ delete [] temp; }
}

/*!
  Given a non-zero pattern, pass in the values for the array
*/
void TACSDistMat::setValues( int nvars, const int *ext_vars,
                             const int *rowp, const int *cols, 
                             TacsScalar *avals ){
  zeroEntries();

  int bsize = Aloc->getBlockSize();
  int b2 = bsize*bsize;

  // Determine the maximum size of the array
  int max_row = 0;
  for ( int i = 0; i < nvars; i++ ){
    int size = rowp[i+1] - rowp[i];
    if (size > max_row){
      max_row = size;
    }
  }

  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);
  
  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank+1];
  
  int *acols = new int[ 2*max_row ];
  int *bcols = &acols[max_row];

  for ( int i = 0; i < nvars; i++ ){
    int row = ext_vars[i];
    
    // Check if this row is in the local or external block
    int nb = 0;
    if (row >= lower && row < upper){
      row = row - lower;

      // Convert the cols
      int start = rowp[i];
      int end   = rowp[i+1];
      for ( int j = rowp[i], k = 0; j < end; j++, k++ ){
        int col = ext_vars[cols[j]];
        acols[k] = -1;
        bcols[k] = -1;

        if (col >= lower && col < upper){
          acols[k] = col - lower;
        }
        else {
          int *item = (int*)bsearch(&col, col_map_vars, col_map_size, 
                                    sizeof(int), FElibrary::comparator);
          bcols[k] = item - col_map_vars;
          nb++;
        }
      }

      Aloc->addBlockRowValues(row, end-start, acols, &avals[b2*start]);

      if (nb > 0){
        row = row - Np;
        if (row >= 0 && row < Nc){
          Bext->addBlockRowValues(row, end-start, bcols, &avals[b2*start]);
        }
        else {
          fprintf(stderr, "[%d] DistMat error: could not find row %d\n", 
                  mpiRank, row);
        }
      }
    }
    else {
      int *item = (int*)bsearch(&row, ext_rows, next_rows, 
                                sizeof(int), FElibrary::comparator);

      if (item){
        int r_ext = item - ext_rows;

        int end   = rowp[i+1];
        for ( int j = rowp[i], k = 0; j < end; j++, k++ ){
          int c = cols[j];
          if (c >= 0 && c < nvars){
            int col = ext_vars[c];

            int ext_start = ext_rowp[r_ext];
            int ext_size = ext_rowp[r_ext+1] - ext_start;
            item = (int*)bsearch(&col, &ext_cols[ext_start], ext_size, 
                                 sizeof(int), FElibrary::comparator);
            
            if (item){
              TacsScalar *a = &ext_A[b2*(item - ext_cols)];    
              memcpy(a, &avals[b2*j], b2*sizeof(TacsScalar));
            }
            else {
              fprintf(stderr, "[%d] DistMat error: could not find col \
(%d,%d) r_ext = %d \n", mpiRank, row, col, r_ext);
            }
          }
          else {
            fprintf(stderr, "[%d] DistMat error: local column out of \
range 0 <= %d < %d\n", mpiRank, c, nvars);
          }
        }
      }
      else {
        fprintf(stderr, "[%d] DistMat error: could not find row %d\n", 
                mpiRank, row);
      }
    }
  }

  delete [] acols;
}

void TACSDistMat::zeroEntries(){
  int bsize = Aloc->getBlockSize();
  Aloc->zeroEntries();
  Bext->zeroEntries();

  int len = bsize*bsize*ext_rowp[next_rows];
  memset(ext_A, 0, len*sizeof(TacsScalar));
}

/* 
   Initiate the communication of the off-process matrix entries 
*/
void TACSDistMat::beginAssembly(){
  if (nreceives > 0){
    int ierr = MPI_Startall(nreceives, receives);
    if (ierr != MPI_SUCCESS){
      int mpiRank;
      MPI_Comm_rank(comm, &mpiRank);
      int len;
      char err_str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(ierr, err_str, &len);
      fprintf(stderr, 
              "[%d] DistMat::beginAssembly MPI startall receives error\n%s\n", 
              mpiRank, err_str);
    }      
  }
  if (nsends > 0){
    int ierr = MPI_Startall(nsends, sends);
    if (ierr != MPI_SUCCESS){ 
      int mpiRank;
      MPI_Comm_rank(comm, &mpiRank);
      int len;
      char err_str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(ierr, err_str, &len);
      fprintf(stderr, 
              "[%d] DistMat::beginAssembly MPI startall sends error\n%s\n", 
              mpiRank, err_str);
    }      
  }
}

/*
  Finish the communication of the off-process matrix entries. 
  Once the communication is completed, add the off-processor
  entries to the matrix.
*/
void TACSDistMat::endAssembly(){
  // Get the map between the global-external variables 
  // and the local variables (for Bext)
  int bsize = Aloc->getBlockSize();
  int b2 = bsize*bsize;
  
  int mpiRank;
  MPI_Comm_rank(comm, &mpiRank);

  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  int lower = ownerRange[mpiRank];
  int upper = ownerRange[mpiRank+1];

  for ( int i = 0; i < nreceives; i++ ){
    int index;
    MPI_Status status;

    int ierr = MPI_Waitany(nreceives, receives, &index, &status);
    if (ierr != MPI_SUCCESS){
      int len;
      char err_str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(ierr, err_str, &len);
      fprintf(stderr, "[%d] DistMat::endAssembly MPI waitany error\n%s\n",
              mpiRank, err_str);
    }
    int n = receive_proc[index];

    for ( int j = in_row_ptr[n]; j < in_row_ptr[n+1]; j++ ){
      int row = in_rows[j] - lower;

      for ( int k = in_rowp[j]; k < in_rowp[j+1]; k++ ){
        TacsScalar *a = &in_A[b2*k];

        int col = in_cols[k];
        if (col >= lower && col < upper){
          col = col - lower;
          Aloc->addBlockRowValues(row, 1, &col, a);
        }
        else {
          int *item = (int*)bsearch(&col, col_map_vars, col_map_size, 
                                    sizeof(int), FElibrary::comparator);
          if (item){
            int c = item - col_map_vars;
            Bext->addBlockRowValues(row - Np, 1, &c, a);
          }
        }
      }      
    }
  }

  // Wait for all the sending requests
  if (nsends > 0){
    MPI_Waitall(nsends, sends, MPI_STATUSES_IGNORE);
  }
}
