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

#include "BVecInterp.h"
#include "FElibrary.h"
#include "MatUtils.h"

/*
  BVecInterp: Interpolate with constant weights between two vectors. 
*/

/*
  These are the definitions for the generic and block-specific
  matrix-vector operations required in the BVecInterp class.
*/
void BVecInterpMultAddGen( int bsize, int nrows, 
                           const int *rowp, const int *cols,
                           const TacsScalar *weights,
                           const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAddGen( int bsize, int nrows, 
                                    const int *rowp, 
                                    const int *cols,
                                    const TacsScalar *weights,
                                    const TacsScalar *x, 
                                    TacsScalar *y );

void BVecInterpMultAdd1( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd1( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );

void BVecInterpMultAdd2( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd2( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );

void BVecInterpMultAdd3( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd3( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );
void BVecInterpMultAdd4( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd4( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );
void BVecInterpMultAdd5( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd5( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );

void BVecInterpMultAdd6( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *weights,
                         const TacsScalar *x, TacsScalar *y );
void BVecInterpMultTransposeAdd6( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *weights,
                                  const TacsScalar *x, 
                                  TacsScalar *y );

/*
  This object represents a matrix that interpolates between
  vectors of different sizes. Each interpolation is performed block-wise
  in the sense that the same weights are applied to each component 
  of the block.

  This code works in the following manner: 

  1. The input/output maps of the operator are defined
  2. Each processor adds an interpolation between the input and output
  3. A call to finalize() is made to finish initialization
  4. Calls can be made to mult, multAdd, multTranspose, and multTransposeAdd
*/
TACSBVecInterp::TACSBVecInterp( TACSVarMap *_inMap,
                                TACSVarMap *_outMap,
                                int _bsize ){
  inTacs = NULL;
  outTacs = NULL;
  init(_inMap, _outMap, _bsize);
}

/*
  Initialize TACSBVecInterp with the TACSAssembler objects themselves.
  
  This is required when the TACSAssembler objects reorder the variables
*/
TACSBVecInterp::TACSBVecInterp( TACSAssembler *_inTacs,
                                TACSAssembler *_outTacs ){
  inTacs = _inTacs;
  inTacs->incref();

  outTacs = _outTacs;
  outTacs->incref();
  init(inTacs->getVarMap(), outTacs->getVarMap(), 
       inTacs->getVarsPerNode());
}

/*
  Initialize the underlying TACSAssembler classes
*/
void TACSBVecInterp::init( TACSVarMap *_inMap,
                           TACSVarMap *_outMap,
                           int _bsize ){
  inMap = _inMap;
  inMap->incref();

  outMap = _outMap;
  outMap->incref();
  
  // The block size to be used 
  bsize = _bsize;

  // Make sure that the MPI communicators from the input/output
  // vectors are the same
  int result;
  MPI_Comm_compare(outMap->getMPIComm(), 
                   inMap->getMPIComm(), &result);
  
  if (!(result == MPI_IDENT || result == MPI_CONGRUENT)){
    fprintf(stderr, "Error in TACSBVecInterp: MPI groups are not idential \
or congruent. Cannot form interpolant.\n");
    return;
  }

  // Record the communicator for later usage
  comm = inMap->getMPIComm();

  // Initialize the on- and off-processor temporary storage that 
  // will be used to store the interpolation weights
  N = outMap->getDim();
  M = inMap->getDim();

  on_size = 0;
  max_on_size = N;
  max_on_weights = 27*max_on_size;
  on_nums = new int[ max_on_size ];
  on_rowp = new int[ max_on_size+1 ];
  on_vars = new int[ max_on_weights ];
  on_weights = new TacsScalar[ max_on_weights ];
  on_rowp[0] = 0;

  off_size = 0;
  max_off_size = int(0.1*N);
  if (max_off_size < 100){
    max_off_size = 100;
  }
  max_off_weights = 27*max_off_size;
  off_nums = new int[ max_off_size ];
  off_rowp = new int[ max_off_size+1 ];
  off_vars = new int[ max_off_weights ];
  off_weights = new TacsScalar[ max_off_weights ];
  off_rowp[0] = 0;

  // For now, set the interpolation data to null
  vecDist = NULL;
  ctx = NULL;
  rowp = NULL;
  cols = NULL;
  weights = NULL;

  ext_rowp = NULL;
  ext_cols = NULL;
  ext_weights = NULL;
  x_ext = NULL;

  // NULL the transpose weight vector
  transpose_weights = NULL;

  // Initialize the implementation
  multadd = BVecInterpMultAddGen;
  multtransadd = BVecInterpMultTransposeAddGen;

  // Initialize the block-specific implementations
  switch (bsize) {
  case 1:
    multadd = BVecInterpMultAdd1;
    multtransadd = BVecInterpMultTransposeAdd1;
    break;
  case 2:
    multadd = BVecInterpMultAdd2;
    multtransadd = BVecInterpMultTransposeAdd2;
    break;
  case 3:
    multadd = BVecInterpMultAdd3;
    multtransadd = BVecInterpMultTransposeAdd3;
    break;
  case 4:
    multadd = BVecInterpMultAdd4;
    multtransadd = BVecInterpMultTransposeAdd4;
    break;
  case 5:
    multadd = BVecInterpMultAdd5;
    multtransadd = BVecInterpMultTransposeAdd5;
    break;
  case 6:
    multadd = BVecInterpMultAdd6;
    multtransadd = BVecInterpMultTransposeAdd6;
    break;
  default:
    break;
  }
}

/*
  Delete the interpolation object, and free all the objects allocated
  internally.  
*/
TACSBVecInterp::~TACSBVecInterp(){
  inMap->decref();
  outMap->decref();
  if (inTacs){ inTacs->decref(); }
  if (outTacs){ outTacs->decref(); }

  // Deallocate data that may have been freed in finalize()
  if (on_nums){ delete [] on_nums; }
  if (on_rowp){ delete [] on_rowp; }
  if (on_vars){ delete [] on_vars; }
  if (on_weights){ delete [] on_weights; }

  if (off_nums){ delete [] off_nums; }
  if (off_rowp){ delete [] off_rowp; }
  if (off_vars){ delete [] off_vars; }
  if (off_weights){ delete [] off_weights; }

  // Deallocate data that may have been allocated in initialize()
  if (vecDist){ vecDist->decref(); }
  if (ctx){ ctx->decref(); }
  if (rowp){ delete [] rowp; }
  if (cols){ delete [] cols; }
  if (weights){ delete [] weights; }
  if (x_ext){ delete [] x_ext; }
  if (ext_rowp){ delete [] ext_rowp; }
  if (ext_cols){ delete [] ext_cols; }
  if (ext_weights){ delete [] ext_weights; }
  if (transpose_weights){ delete [] transpose_weights; }
}

/*
  Add an interopolation between an output variable, and a series of
  input variables. Variables can be added from anywhere to anywhere,
  but it is more efficient if variables are primarily added on the
  processors to which they belong.
  
  This should be called for every variable in the
  interpolation/extrapolation.

  input:
  num:      the interpolation variable number
  weights:  the interpolation weights
  vars:     the interpolating variable numbers
  size:     the number of points used in the interpolation
*/
void TACSBVecInterp::addInterp( int num, TacsScalar w[], 
                                int vars[], int size ){
  // Get the MPI rank
  int mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);

  // Get the ownership range
  const int *outOwnerRange;
  outMap->getOwnerRange(&outOwnerRange);
  
  if (num >= outOwnerRange[mpi_rank] &&
      num < outOwnerRange[mpi_rank+1]){
    // This code stores the values locally until the finalize call is
    // made. First, check if the current size of the allocated memory
    // is sufficient to store the interpolation. If not, allocate
    // larger arrays.
    if (on_size >= max_on_size){
      // Double the current size of the array and copy the old 
      // values into the newly allocated part.
      max_on_size += max_on_size;
      int *tmp_nums = on_nums;
      int *tmp_rowp = on_rowp;
      on_nums = new int[ max_on_size ];
      on_rowp = new int[ max_on_size+1 ];

      memcpy(on_nums, tmp_nums, on_size*sizeof(int));
      memcpy(on_rowp, tmp_rowp, (on_size+1)*sizeof(int));
      delete [] tmp_nums;
      delete [] tmp_rowp;
    }
    if (on_rowp[on_size] + size > max_on_weights){
      // Increase the size of the working pointer array, and
      // copy over the old values into the new array
      max_on_weights += size + max_on_weights;
      int *tmp_vars = on_vars;
      TacsScalar *tmp_weights = on_weights;
      on_vars = new int[ max_on_weights ];
      on_weights = new TacsScalar[ max_on_weights ];

      memcpy(on_vars, tmp_vars, on_rowp[on_size]*sizeof(int));
      memcpy(on_weights, tmp_weights, 
             on_rowp[on_size]*sizeof(TacsScalar));
      delete [] tmp_vars;
      delete [] tmp_weights;
    }
    
    // Store the values temporarily 
    on_nums[on_size] = num;
    on_rowp[on_size+1] = on_rowp[on_size] + size;
    for ( int i = 0, j = on_rowp[on_size]; i < size; i++, j++ ){
      on_vars[j] = vars[i];
      on_weights[j] = w[i];
    }
    on_size++;
  }
  else {
    // Add the values to the off-processor part. This will store
    // the values locally until the finalize call is made. First,
    // check if the current size of the allocated memory is sufficient
    // to store the interpolation. If not, allocate larger arrays.
    if (off_size >= max_off_size){
      // Double the current size of the array and copy the old 
      // values into the newly allocated part.
      max_off_size += max_off_size;
      int *tmp_nums = off_nums;
      int *tmp_rowp = off_rowp;
      off_nums = new int[ max_off_size ];
      off_rowp = new int[ max_off_size+1 ];

      memcpy(off_nums, tmp_nums, off_size*sizeof(int));
      memcpy(off_rowp, tmp_rowp, (off_size+1)*sizeof(int));
      delete [] tmp_nums;
      delete [] tmp_rowp;
    }
    if (off_rowp[off_size] + size > max_off_weights){
      // Increase the size of the working pointer array, and
      // copy over the old values into the new array
      max_off_weights += size + max_off_weights;
      int *tmp_vars = off_vars;
      TacsScalar *tmp_weights = off_weights;
      off_vars = new int[ max_off_weights ];
      off_weights = new TacsScalar[ max_off_weights ];

      memcpy(off_vars, tmp_vars, off_rowp[off_size]*sizeof(int));
      memcpy(off_weights, tmp_weights, 
             off_rowp[off_size]*sizeof(TacsScalar));
      delete [] tmp_vars;
      delete [] tmp_weights;
    }
    
    // Store the values temporarily 
    off_nums[off_size] = num;
    off_rowp[off_size+1] = off_rowp[off_size] + size;
    for ( int i = 0, j = off_rowp[off_size]; i < size; i++, j++ ){
      off_vars[j] = vars[i];
      off_weights[j] = w[i];
    }
    off_size++;
  }
}

/*
  Initialize the interpolation and set up the internal data structures
  so that the object can be used for interpolation/extrapolation.
  This code is collective on all processors in the communicator.

  This finalize call performs a few critical tasks: 

  1. All interpolation weights are passed to the processors which own
  them.  

  2. The interpolation is divided into two parts: the local part and
  the external part. The local part only acts on local variables,
  while the external part acts on any external variable.

  3. All weights are normalized by the sum of the weights in each row.

  After these operations, the TACSBVecInterp class can be utilized to
  perform restriction or prolongation operations.
*/
void TACSBVecInterp::initialize(){
  // Retrieve the MPI comm size
  int mpi_size, mpi_rank;
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);

  // Get the ownership range
  const int *outOwnerRange, *inOwnerRange;
  outMap->getOwnerRange(&outOwnerRange);
  inMap->getOwnerRange(&inOwnerRange);

  // Find the number of contributions that need to be sent to other
  // processors. Count the number of times each off-processor equation
  // is referenced.
  int *tmp_count = new int[ mpi_size ];
  int *tmp_ptr = new int[ mpi_size+1 ];
  int *tmp_weights_count = new int[ mpi_size ];
  int *tmp_weights_ptr = new int[ mpi_size+1 ];

  memset(tmp_count, 0, mpi_size*sizeof(int));
  memset(tmp_weights_count, 0, mpi_size*sizeof(int));

  for ( int i = 0; i < off_size; i++ ){
    int index = FElibrary::findInterval(off_nums[i], 
                                        outOwnerRange, mpi_size+1);
    tmp_count[index]++;
    tmp_weights_count[index] += off_rowp[i+1] - off_rowp[i];
  }

  // Set a pointer array into the newly allocated block of memory that
  // will store the arranged interpolants
  tmp_ptr[0] = 0;
  tmp_weights_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    tmp_ptr[i+1] = tmp_ptr[i] + tmp_count[i];
    tmp_weights_ptr[i+1] = tmp_weights_ptr[i] + tmp_weights_count[i];
  }

  // Allocate memory for the temporary data storage. Note that we
  // store the number of weights per row, rather than a pointer to
  // data (like CSR format) since this will be easier to transfer
  // between processors since the relative offset won't change.
  int *tmp_nums = new int[ off_size ];
  int *tmp_weights_per_row = new int[ off_size ];
  int *tmp_vars = new int[ off_rowp[off_size] ];
  TacsScalar *tmp_weights = new TacsScalar[ off_rowp[off_size] ];

  // Go through the entire list of interpolant contributions and 
  // copy over the data to the required off-diagonal spot in the
  // temporary data arrays.
  for ( int i = 0; i < off_size; i++ ){
    int index = FElibrary::findInterval(off_nums[i], 
                                        outOwnerRange, mpi_size+1);

    if (!(index < 0 || index >= mpi_size)){
      // Copy the values over to a temporary array
      tmp_nums[tmp_ptr[index]] = off_nums[i];
      
      // Copy the weight and variable information to the new array
      int j = tmp_weights_ptr[index];
      for ( int k = off_rowp[i]; k < off_rowp[i+1]; k++, j++ ){
        tmp_weights[j] = off_weights[k];
        tmp_vars[j] = off_vars[k];
      }
      
      // Record the number of weights stored for this rows
      tmp_weights_per_row[tmp_ptr[index]] = off_rowp[i+1] - off_rowp[i];

      // Increment the pointers to where the next data will be added
      tmp_ptr[index]++;
      tmp_weights_ptr[index] += off_rowp[i+1] - off_rowp[i];
    }
  }

  // Reset the pointer arrays so that they again point to the
  // beginning of each block of data to be sent to the other
  // processors
  tmp_ptr[0] = 0;
  tmp_weights_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    tmp_ptr[i+1] = tmp_ptr[i] + tmp_count[i];
    tmp_weights_ptr[i+1] = tmp_weights_ptr[i] + tmp_weights_count[i];
  }

  // Delete the actual values of the off-diagonal contributions - they
  // are no-longer required since we have copied all the information
  // to the temporary arrays.
  delete [] off_nums;
  delete [] off_rowp;
  delete [] off_vars;
  delete [] off_weights;

  // Send the number of out-going interpolation values to the
  // recieving processors.
  int *in_count = new int[ mpi_size ];
  int *in_ptr = new int[ mpi_size+1 ];
  int *in_weights_count = new int[ mpi_size ];
  int *in_weights_ptr = new int[ mpi_size+1 ];

  // Send/recieve the in-coming counts
  MPI_Alltoall(tmp_count, 1, MPI_INT, in_count, 1, MPI_INT, comm);

  // Set a pointer for the incoming num/count data
  in_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    in_ptr[i+1] = in_ptr[i] + in_count[i];
  }

  // Send/recieve the in-coming weight counts
  MPI_Alltoall(tmp_weights_count, 1, MPI_INT, 
               in_weights_count, 1, MPI_INT, comm);

  // Set a pointer for the incoming weights/variables data
  in_weights_ptr[0] = 0;
  for ( int i = 0; i < mpi_size; i++ ){
    in_weights_ptr[i+1] = in_weights_ptr[i] + in_weights_count[i];
  }

  // Allocate the required space for all the incoming data
  int *in_nums = new int[ in_ptr[mpi_size] ];
  int *in_weights_per_row = new int[ in_ptr[mpi_size] ];
  int *in_vars = new int[ in_weights_ptr[mpi_size] ];
  TacsScalar *in_weights = new TacsScalar[ in_weights_ptr[mpi_size] ];

  // Send and recieve all the data destined for this processor
  // Send the variable numbers
  MPI_Alltoallv(tmp_nums, tmp_count, tmp_ptr, MPI_INT,
                in_nums, in_count, in_ptr, MPI_INT, comm);
  
  // Send the number of variables per row
  MPI_Alltoallv(tmp_weights_per_row, tmp_count, tmp_ptr, MPI_INT,
                in_weights_per_row, in_count, in_ptr, MPI_INT, comm);

  // Send the variables for each row
  MPI_Alltoallv(tmp_vars, tmp_weights_count, tmp_weights_ptr, MPI_INT,
                in_vars, in_weights_count, in_weights_ptr, MPI_INT, comm);

  // Send the weights for each variable
  MPI_Alltoallv(tmp_weights, tmp_weights_count, tmp_weights_ptr, TACS_MPI_TYPE,
                in_weights, in_weights_count, in_weights_ptr, 
                TACS_MPI_TYPE, comm);

  // Delete the temporary data - this is no longer required since we
  // now have the final data on the actual processors
  delete [] tmp_count;
  delete [] tmp_ptr;
  delete [] tmp_weights_count;
  delete [] tmp_weights_ptr;
  delete [] tmp_nums;
  delete [] tmp_weights_per_row;
  delete [] tmp_vars;
  delete [] tmp_weights;

  // Delete the pointers that indicated where the data came from since
  // this information is not relevant, but record how many total entries
  // were sent to this processor.
  int in_size = in_ptr[mpi_size];
  delete [] in_count;
  delete [] in_ptr;
  delete [] in_weights_count;
  delete [] in_weights_ptr;

  // Now all the data required to assemble the internal data
  // structures are on the local processors. Now, we just assemble the
  // on and off-processor parts, but all assignments (i.e. all the
  // output) is local.

  // Allocate space for the internal and external portions of the
  // matrix
  rowp = new int[ N+1 ];
  memset(rowp, 0, (N+1)*sizeof(int));

  ext_rowp = new int[ N+1 ];
  memset(ext_rowp, 0, (N+1)*sizeof(int));  

  // Count up the contributions to the CSR data structures from the
  // on-processor data
  for ( int i = 0; i < on_size; i++ ){
    // Compute the on-processor variable number
    int num = on_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
      
      // Count up the number of internal and external variables
      // in this row. Note that local/external variables are based
      // on the input variable range!
      int num_local = 0, num_ext = 0;
      
      for ( int j = on_rowp[i]; j < on_rowp[i+1]; j++ ){
        if (on_vars[j] >= inOwnerRange[mpi_rank] &&
            on_vars[j] < inOwnerRange[mpi_rank+1]){
          num_local++;
        }
        else {
          num_ext++;
        }
      }
        
      rowp[num+1] += num_local;
      ext_rowp[num+1] += num_ext;
    }
    else {
      // Print an error message. This should never happen. 
      fprintf(stderr, "Error, local interpolation variable out of range\n");
    }
  }

  // Count up the contributions to the CSR data structures from the
  // off-processor data now stored in the input arrays
  for ( int i = 0, k = 0; i < in_size; i++ ){
    // Compute the on-processor variable number
    int num = in_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
      
      // Count up the number of internal and external variables
      // in this row. Note that local/external variables are based
      // on the input variable range!
      int num_local = 0, num_ext = 0;
      
      for ( int j = 0; j < in_weights_per_row[i]; j++, k++ ){
        if (in_vars[k] >= inOwnerRange[mpi_rank] &&
            in_vars[k] < inOwnerRange[mpi_rank+1]){
          num_local++;
        }
        else {
          num_ext++;
        }
      }
        
      rowp[num+1] += num_local;
      ext_rowp[num+1] += num_ext;
    }
    else {
      // Print an error message. This should never happen.
      fprintf(stderr, "Error, local interpolation variable out of range\n");
    
      // Increment the pointer to the input array
      k += in_weights_per_row[i];
    }
  }

  // Increment the rowp pointers to the local and external data so
  // that we have not just the counts, but the offsets into the global
  // variable arrays
  rowp[0] = 0;
  ext_rowp[0] = 0;
  for ( int i = 0; i < N; i++ ){
    rowp[i+1] += rowp[i];
    ext_rowp[i+1] += ext_rowp[i];
  }

  // Now, compute the size of the overall arrays required to store all
  // of this data
  cols = new int[ rowp[N] ];
  ext_cols = new int[ ext_rowp[N] ];

  // Add the contributions to the CSR data structure from the
  // on-processor data. Note that this adjusts the rowp and ext_rowp
  // data. This must be restored afterwards.
  for ( int i = 0; i < on_size; i++ ){
    // Compute the on-processor variable number
    int num = on_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
      
      for ( int j = on_rowp[i]; j < on_rowp[i+1]; j++ ){
        if (on_vars[j] >= inOwnerRange[mpi_rank] &&
            on_vars[j] < inOwnerRange[mpi_rank+1]){
          int index = on_vars[j];
          if (inTacs){
            inTacs->reorderNodes(&index, 1);
          }
          cols[rowp[num]] = index;
          rowp[num]++;
        }
        else {
          ext_cols[ext_rowp[num]] = on_vars[j];
          ext_rowp[num]++;
        }
      }
    }
  }

  // Add the entries into the CSR data structure from the 
  // off-processor data
  for ( int i = 0, k = 0; i < in_size; i++ ){
    // Compute the on-processor variable number
    int num = in_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
            
      for ( int j = 0; j < in_weights_per_row[i]; j++, k++ ){
        if (in_vars[k] >= inOwnerRange[mpi_rank] &&
            in_vars[k] < inOwnerRange[mpi_rank+1]){
          int index = in_vars[k];
          if (inTacs){
            inTacs->reorderNodes(&index, 1);
          }
          cols[rowp[num]] = index;
          rowp[num]++;
        }
        else {
          ext_cols[ext_rowp[num]] = in_vars[k];
          ext_rowp[num]++;
        }
      }
    }
    else {
      // Increment the pointer to the input array
      k += in_weights_per_row[i];
    }
  }

  // Adjust the rowp/ext_rowp array so that they again point to their
  // proper positions
  for ( int i = N; i > 0; i-- ){
    rowp[i] = rowp[i-1];
    ext_rowp[i] = ext_rowp[i-1];
  }
  rowp[0] = 0;
  ext_rowp[0] = 0;

  // Sort and uniquify the CSR data structures so that we don't have
  // duplicate entries
  int nodiag = 0; // Don't remove the diagonal from the matrix
  matutils::SortAndUniquifyCSR(N, rowp, cols, nodiag);
  matutils::SortAndUniquifyCSR(N, ext_rowp, ext_cols, nodiag);

  // Allocate space for the weights. Initialize the weight values
  // to zero
  weights = new TacsScalar[ rowp[N] ];
  memset(weights, 0, rowp[N]*sizeof(TacsScalar));

  ext_weights = new TacsScalar[ ext_rowp[N] ];
  memset(ext_weights, 0, ext_rowp[N]*sizeof(TacsScalar));  

  // Add the weight values themselves to the CSR data structure
  for ( int i = 0; i < on_size; i++ ){
    // Compute the on-processor variable number
    int num = on_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
      
      for ( int j = on_rowp[i]; j < on_rowp[i+1]; j++ ){
        if (on_vars[j] >= inOwnerRange[mpi_rank] &&
            on_vars[j] < inOwnerRange[mpi_rank+1]){
          int index = on_vars[j];
          if (inTacs){
            inTacs->reorderNodes(&index, 1);
          }
          int size = rowp[num+1] - rowp[num];
          int *item = (int*)bsearch(&index, &cols[rowp[num]], size,  
                                    sizeof(int), FElibrary::comparator);
          if (item){
            weights[item - cols] += on_weights[j]; 
          }
        }
        else {
          int size = ext_rowp[num+1] - ext_rowp[num];
          int *item = (int*)bsearch(&on_vars[j], 
                                    &ext_cols[ext_rowp[num]], size,  
                                    sizeof(int), FElibrary::comparator);
          if (item){
            ext_weights[item - ext_cols] += on_weights[j]; 
          }
        }
      }
    }
  }

  // Add the entries into the CSR data structure from the
  // off-processor data
  for ( int i = 0, k = 0; i < in_size; i++ ){
    // Compute the on-processor variable number
    int num = in_nums[i];
    
    if (num >= outOwnerRange[mpi_rank] &&
        num < outOwnerRange[mpi_rank+1]){
      // Adjust the range of the output variable to the local index
      if (outTacs){
        outTacs->reorderNodes(&num, 1);
      }
      num = num - outOwnerRange[mpi_rank];
            
      for ( int j = 0; j < in_weights_per_row[i]; j++, k++ ){
        if (in_vars[k] >= inOwnerRange[mpi_rank] &&
            in_vars[k] < inOwnerRange[mpi_rank+1]){
          int index = in_vars[k];
          if (inTacs){
            inTacs->reorderNodes(&index, 1);
          }
          int size = rowp[num+1] - rowp[num];
          int *item = (int*)bsearch(&index, &cols[rowp[num]], size,  
                                    sizeof(int), FElibrary::comparator);
          if (item){
            weights[item - cols] += in_weights[k]; 
          }
        }
        else {
          int size = ext_rowp[num+1] - ext_rowp[num];
          int *item = (int*)bsearch(&in_vars[k], 
                                    &ext_cols[ext_rowp[num]], size,  
                                    sizeof(int), FElibrary::comparator);
          if (item){
            ext_weights[item - ext_cols] += in_weights[k]; 
          }
        }
      }
    }
    else {
      // Increment the pointer to the input array
      k += in_weights_per_row[i];
    }
  }

  // Delete the incoming data
  delete [] in_nums;
  delete [] in_weights_per_row;
  delete [] in_vars;
  delete [] in_weights;
  
  // Delete the on-processor data
  delete [] on_nums;
  delete [] on_rowp;
  delete [] on_vars;
  delete [] on_weights;

  // Adjust both the internal and external CSR data structures to
  // reflect the internal ordering. We order the local on-processor
  // components of the input map based on shifting the variable input
  // from inOwnerRange[mpi_rank] to start at zero. The external data
  // ordering is based on the ext_vars array.
  for ( int j = 0; j < rowp[N]; j++ ){
    cols[j] = cols[j] - inOwnerRange[mpi_rank];
  }

  // Allocate space for the external variable numbers
  int *ext_vars = new int[ ext_rowp[N] ];
  memcpy(ext_vars, ext_cols, ext_rowp[N]*sizeof(int));
  num_ext_vars = FElibrary::uniqueSort(ext_vars, ext_rowp[N]);

  // Check if the version of TACS has been reordered. If so, 
  if (inTacs && inTacs->isReordered()){
    // Match the intervals for the external node numbers
    int *ext_ptr = new int[ mpi_size+1 ];
    int *ext_count = new int[ mpi_size ];
    FElibrary::matchIntervals(mpi_size, inOwnerRange,
                              num_ext_vars, ext_vars, ext_ptr);

    // Send the nodes owned by other processors the information. First
    // count up how many will go to each process.
    for ( int i = 0; i < mpi_size; i++ ){
      ext_count[i] = ext_ptr[i+1] - ext_ptr[i];
    }

    int *recv_count = new int[ mpi_size ];
    int *recv_ptr = new int[ mpi_size+1 ];
    MPI_Alltoall(ext_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    // Now prepare to send the node numbers to the other processors
    recv_ptr[0] = 0;
    for ( int i = 0; i < mpi_size; i++ ){
      recv_ptr[i+1] = recv_ptr[i] + recv_count[i];
    }

    // Number of nodes that will be received from other procs
    int *recv_vars = new int[ recv_ptr[mpi_size] ];
    MPI_Alltoallv(ext_vars, ext_count, ext_ptr, MPI_INT, 
                  recv_vars, recv_count, recv_ptr, MPI_INT, comm);

    // Reorder the variables that are local to this processor
    inTacs->reorderNodes(recv_vars, recv_ptr[mpi_size]);  

    // Send the new variables back to the original processors in the
    // same order as we recvd them
    int *new_ext_vars = new int[ num_ext_vars ];
    MPI_Alltoallv(recv_vars, recv_count, recv_ptr, MPI_INT,
                  new_ext_vars, ext_count, ext_ptr, MPI_INT, comm);

    // Adjust the ordering of the external variables
    for ( int j = 0; j < ext_rowp[N]; j++ ){
      int *item = (int*)bsearch(&ext_cols[j], ext_vars, num_ext_vars,
                                sizeof(int), FElibrary::comparator);
      ext_cols[j] = new_ext_vars[item - ext_vars];
    }
    
    // Free the old ext_vars array - this is no longer required
    delete [] ext_vars;

    // Set the new ext vars array and sort it
    ext_vars = new_ext_vars;
    FElibrary::uniqueSort(ext_vars, num_ext_vars);

    // Free the recv_vars array - it is no longer required
    delete [] recv_vars;
    delete [] ext_ptr;
    delete [] ext_count;
    delete [] recv_count;
    delete [] recv_ptr;
  }

  // Adjust the external ordering
  for ( int j = 0; j < ext_rowp[N]; j++ ){
    int *item = (int*)bsearch(&ext_cols[j], ext_vars, num_ext_vars,
                              sizeof(int), FElibrary::comparator);
    ext_cols[j] = item - ext_vars;
  }

  // ext_rowp/ext_cols/ext_weights now contains the off-diagonal
  // components this object now can distribute them
  TACSBVecIndices *bindex = new TACSBVecIndices(&ext_vars, num_ext_vars);
  vecDist = new TACSBVecDistribute(inMap, bindex); 
  vecDist->incref();

  // Create a vect distribution context
  ctx = vecDist->createCtx(bsize);
  ctx->incref();

  // Allocate memroy for the in-coming data
  x_ext = new TacsScalar[ bsize*num_ext_vars ];

  // Set all the pointers to null
  on_nums = NULL;
  on_rowp = NULL;
  on_vars = NULL;
  on_weights = NULL;

  off_nums = NULL;
  off_rowp = NULL;
  off_vars = NULL;
  off_weights = NULL;

  // Normalize the weights across both the internal/external mappings
  for ( int i = 0; i < N; i++ ){
    TacsScalar w = 0.0;

    // Compute the sum of the internal and external weights
    for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
      w += weights[j];
    }
    for ( int j = ext_rowp[i]; j < ext_rowp[i+1]; j++ ){
      w += ext_weights[j];
    }

    if (w != 0.0){
      // Compute the inverse of the weights
      w = 1.0/w;

      // Normalize the weights
      for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
        weights[j] *= w;
      }
      for ( int j = ext_rowp[i]; j < ext_rowp[i+1]; j++ ){
        ext_weights[j] *= w;
      }
    }
  }

  // Allocate the transpose weights
  transpose_weights = new TacsScalar[ M ];
  memset(transpose_weights, 0, M*sizeof(TacsScalar));

  // Zero the external contribution and begin
  memset(x_ext, 0, num_ext_vars*sizeof(TacsScalar));
  for ( int i = 0; i < N; i++ ){
    for ( int jp = ext_rowp[i]; jp < ext_rowp[i+1]; jp++ ){
      x_ext[ext_cols[jp]] += ext_weights[jp];
    }
  }

  // Allocate the context used for communicating the weights
  TACSBVecDistCtx *ctx_weights = vecDist->createCtx(1);
  ctx_weights->incref();
  
  // Begin transfering the column sums 
  vecDist->beginReverse(ctx_weights, x_ext, transpose_weights, TACS_ADD_VALUES);
  
  for ( int i = 0; i < N; i++ ){
    for ( int jp = rowp[i]; jp < rowp[i+1]; jp++ ){
      transpose_weights[cols[jp]] += weights[jp];
    }
  }

  // End transfering the weights
  vecDist->endReverse(ctx_weights, x_ext, transpose_weights, TACS_ADD_VALUES);
  ctx_weights->decref();

  // Compute the inverse of the weights
  for ( int i = 0; i < M; i++ ){
    if (transpose_weights[i] != 0.0){
      transpose_weights[i] = 1.0/transpose_weights[i];
    }
  }
  
  // Zero the external vector
  memset(x_ext, 0, bsize*num_ext_vars*sizeof(TacsScalar));
}

/*
  Perform the interpolation from inVec to outVec:

  Interp*inVec -> outVec

  input:
  inVec:  the input vector

  output:
  outVec: the interpolated output vector
*/
void TACSBVecInterp::mult( TACSBVec *inVec, TACSBVec *outVec ){
  if (!vecDist){
    fprintf(stderr, 
            "Must call initialize() before using TACSBVecInterp object\n");
    return;
  }

  // Zero the entries
  outVec->zeroEntries();

  // Get the internal arrays for in/out
  TacsScalar *in, *out;
  inVec->getArray(&in);
  outVec->getArray(&out);

  // Initialize the communication from off-processor components to the
  // on-processor values
  vecDist->beginForward(ctx, in, x_ext);

  // Multiply the on-processor part
  multadd(bsize, N, rowp, cols, weights, in, out);

  // Finish the off-processor communication
  vecDist->endForward(ctx, in, x_ext);

  // Multiply the off-processor part
  multadd(bsize, N, ext_rowp, ext_cols, ext_weights, x_ext, out);
}

/*
  Perform the interpolation from inVec to outVec:

  addVec + Interp*inVec -> outVec

  input:
  inVec:  the input vector
  addVec: the vector to add to the output

  output:
  outVec: the interpolated output vector
*/
void TACSBVecInterp::multAdd( TACSBVec *inVec, TACSBVec *addVec, 
                              TACSBVec *outVec ){
  if (!vecDist){
    fprintf(stderr, 
            "Must call initialize() before using TACSBVecInterp object\n");
    return;
  }

  // If the two vectors are not the same, copy the values to
  // the outVec
  if (outVec != addVec){
    outVec->copyValues(addVec);
  }

  // Get the internal arrays for in/out
  TacsScalar *in, *out;
  inVec->getArray(&in);
  outVec->getArray(&out);

  // Initialize the communication from off-processor components to the
  // on-processor values
  vecDist->beginForward(ctx, in, x_ext);

  // Multiply the on-processor part
  multadd(bsize, N, rowp, cols, weights, in, out);

  // Finish the off-processo communication
  vecDist->endForward(ctx, in, x_ext);

  // Multiply the off-processor part
  multadd(bsize, N, ext_rowp, ext_cols, ext_weights, x_ext, out);
}

/*
  Perform the interpolation from inVec to outVec:

  Interp*inVec -> outVec

  input:
  inVec:  the input vector

  output:
  outVec: the interpolated output vector
*/
void TACSBVecInterp::multTranspose( TACSBVec *inVec, TACSBVec *outVec ){
  if (!vecDist){
    fprintf(stderr, 
            "Must call initialize() before using TACSBVecInterp object\n");
    return;
  }

  // Zero the entries of the output array
  outVec->zeroEntries();

  // Get the local arrays
  TacsScalar *in, *out;
  inVec->getArray(&in);
  outVec->getArray(&out);

  // Zero the off-processor contribution
  memset(x_ext, 0, bsize*num_ext_vars*sizeof(TacsScalar));

  // Multiply the off-processor part first
  multtransadd(bsize, N, ext_rowp, ext_cols, ext_weights, in, x_ext);

  // Initialize communication to the off-processor part
  vecDist->beginReverse(ctx, x_ext, out, TACS_ADD_VALUES);

  // Multiply the on-processor part
  multtransadd(bsize, N, rowp, cols, weights, in, out);
  
  // Finalize the communication to the off-processor part
  vecDist->endReverse(ctx, x_ext, out, TACS_ADD_VALUES);
}

/*
  Perform the interpolation from inVec to outVec:

  addVec + Interp*inVec -> outVec

  input:
  inVec:  the input vector
  addVec: the vector to add to the output

  output:
  outVec: the interpolated output vector
*/
void TACSBVecInterp::multTransposeAdd( TACSBVec *inVec, TACSBVec *addVec, 
                                       TACSBVec *outVec ){
  if (!vecDist){
    fprintf(stderr, 
            "Must call initialize() before using TACSBVecInterp object\n");
    return;
  }

  // If the add and output vectors are different, copy the values
  // over to the new entry
  if (outVec != addVec){
    outVec->copyValues(addVec);
  }

  // Get the local arrays
  TacsScalar *in, *out;
  inVec->getArray(&in);
  outVec->getArray(&out);

  // Zero the off-processor contribution
  memset(x_ext, 0, bsize*num_ext_vars*sizeof(TacsScalar));

  // Multiply the off-processor part first
  multtransadd(bsize, N, ext_rowp, ext_cols, ext_weights, in, x_ext);

  // Initialize communication to the off-processor part
  vecDist->beginReverse(ctx, x_ext, out, TACS_ADD_VALUES);

  // Multiply the on-processor part
  multtransadd(bsize, N, rowp, cols, weights, in, out);
  
  // Finalize the communication to the off-processor part
  vecDist->endReverse(ctx, x_ext, out, TACS_ADD_VALUES);
}

/*
  Perform the weighted transpose interpolation

  outVec <- Interp*inVec

  input:
  inVec:  the input vector

  output:
  outVec: the interpolated output vector
*/
void TACSBVecInterp::multWeightTranspose( TACSBVec *inVec, 
                                          TACSBVec *outVec ){
  // Perform the transpose operation
  multTranspose(inVec, outVec);

  // Retrieve the array from the output vector
  TacsScalar *out;
  outVec->getArray(&out);

  // Normalize each component of the output vector by the weights
  for ( int i = 0; i < M; i++ ){
    for ( int k = 0; k < bsize; k++, out++ ){
      out[0] *= transpose_weights[i];
    }
  }
}

/*
  Print the weights to the specified file name
*/
void TACSBVecInterp::printInterp( const char *filename ){
  FILE *fp = fopen(filename, "w");
  if (fp){
    fprintf(fp, "TACSBVecInterp \n");
    for ( int i = 0; i < N; i++ ){
      fprintf(fp, "Row: %d\n", i);
      
      for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
        fprintf(fp, "(%d,%f) ", cols[j], TacsRealPart(weights[j]));
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
}

/*
  The following are the block-specific and generic code for the
  matrix-vector multiplications required within the BVecInterp class.
  The idea is that these will run faster than a generic implementation
  of the matrix-vector multiplication. These will be most important
  for very large meshes where the number of interpolations requried
  (say within a multigrid algorithm) is considerably higher.  
*/

/*
  Compute a matrix-vector product for generic bsize
*/
void BVecInterpMultAddGen( int bsize, int nrows, 
                           const int *rowp, const int *cols,
                           const TacsScalar *w,
                           const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      for ( int k = 0; k < bsize; k++ ){
        y[bsize*i+k] += w[0]*x[bsize*cols[j]+k];
      }
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for generic bsize
*/
void BVecInterpMultTransposeAddGen( int bsize, int nrows, 
                                    const int *rowp, const int *cols,
                                    const TacsScalar *w,
                                    const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      for ( int k = 0; k < bsize; k++ ){
        y[bsize*cols[j]+k] += w[0]*x[bsize*i+k];
      }
      w++;
    }
  }
}

/*
  Compute a matrix-vector product for bsize = 1
*/
void BVecInterpMultAdd1( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      y[i] += w[0]*x[cols[j]];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 1
*/
void BVecInterpMultTransposeAdd1( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[cols[j]] += w[0]*x[i];
      w++;
    }
  }
}

/*
  Compute a matrix-vector product for bsize = 2
*/
void BVecInterpMultAdd2( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[2*i]   += w[0]*x[2*cols[j]];
      y[2*i+1] += w[0]*x[2*cols[j]+1];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 2
*/
void BVecInterpMultTransposeAdd2( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[2*cols[j]]   += w[0]*x[2*i];
      y[2*cols[j]+1] += w[0]*x[2*i+1];
      w++;
    }
  }
}

/*
  Compute a matrix-vector product for bsize = 3
*/
void BVecInterpMultAdd3( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      y[3*i]   += w[0]*x[3*cols[j]];
      y[3*i+1] += w[0]*x[3*cols[j]+1];
      y[3*i+2] += w[0]*x[3*cols[j]+2];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 3
*/
void BVecInterpMultTransposeAdd3( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[3*cols[j]]   += w[0]*x[3*i];
      y[3*cols[j]+1] += w[0]*x[3*i+1];
      y[3*cols[j]+2] += w[0]*x[3*i+2];
      w++;
    }
  }
}
/*
  Compute a matrix-vector product for bsize = 4
*/
void BVecInterpMultAdd4( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      y[4*i]   += w[0]*x[4*cols[j]];
      y[4*i+1] += w[0]*x[4*cols[j]+1];
      y[4*i+2] += w[0]*x[4*cols[j]+2];
      y[4*i+3] += w[0]*x[4*cols[j]+3];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 4
*/
void BVecInterpMultTransposeAdd4( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[4*cols[j]]   += w[0]*x[4*i];
      y[4*cols[j]+1] += w[0]*x[4*i+1];
      y[4*cols[j]+2] += w[0]*x[4*i+2];
      y[4*cols[j]+3] += w[0]*x[4*i+3];
      w++;
    }
  }
}
/*
  Compute a matrix-vector product for bsize = 5
*/
void BVecInterpMultAdd5( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      y[5*i]   += w[0]*x[5*cols[j]];
      y[5*i+1] += w[0]*x[5*cols[j]+1];
      y[5*i+2] += w[0]*x[5*cols[j]+2];
      y[5*i+3] += w[0]*x[5*cols[j]+3];
      y[5*i+4] += w[0]*x[5*cols[j]+4];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 5
*/
void BVecInterpMultTransposeAdd5( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[5*cols[j]]   += w[0]*x[5*i];
      y[5*cols[j]+1] += w[0]*x[5*i+1];
      y[5*cols[j]+2] += w[0]*x[5*i+2];
      y[5*cols[j]+3] += w[0]*x[5*i+3];
      y[5*cols[j]+4] += w[0]*x[5*i+4];
      w++;
    }
  }
}

/*
  Compute a matrix-vector product for bsize = 5
*/
void BVecInterpMultAdd6( int bsize, int nrows, 
                         const int *rowp, const int *cols,
                         const TacsScalar *w,
                         const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];
    
    for (; j < end; j++ ){
      y[6*i]   += w[0]*x[6*cols[j]];
      y[6*i+1] += w[0]*x[6*cols[j]+1];
      y[6*i+2] += w[0]*x[6*cols[j]+2];
      y[6*i+3] += w[0]*x[6*cols[j]+3];
      y[6*i+4] += w[0]*x[6*cols[j]+4];
      y[6*i+5] += w[0]*x[6*cols[j]+5];
      w++;
    }
  }
}

/*
  Compute the matrix-vector transpose product for bsize = 6
*/
void BVecInterpMultTransposeAdd6( int bsize, int nrows, 
                                  const int *rowp, const int *cols,
                                  const TacsScalar *w,
                                  const TacsScalar *x, TacsScalar *y ){
  for ( int i = 0; i < nrows; i++ ){
    int j = rowp[i];
    int end = rowp[i+1];

    for (; j < end; j++ ){
      y[6*cols[j]]   += w[0]*x[6*i];
      y[6*cols[j]+1] += w[0]*x[6*i+1];
      y[6*cols[j]+2] += w[0]*x[6*i+2];
      y[6*cols[j]+3] += w[0]*x[6*i+3];
      y[6*cols[j]+4] += w[0]*x[6*i+4];
      y[6*cols[j]+5] += w[0]*x[6*i+5];
      w++;
    }
  }
}
