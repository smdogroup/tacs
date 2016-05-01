#include "BVecDist.h"
#include "FElibrary.h"

/*
  Distribute/collect block-vector code

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*!
  BVecIndices class definitions
*/

/*
  The static array arg_sort_list and the function compare_arg_sort are
  designed to be used to sort an array of integers 'such that:

  arg_sort_list[array[i]] 

  is sorted in non-decreasing order. This is useful for reverse
  look-up where given a value, you might like to find the index that
  contains that values - if it exists.

  Care should be exercised when using this function. You must always
  assign the pointer arg_sort_list before using it, because other
  functions may use it as well.
*/

static const int * arg_sort_list = NULL;

static int compare_arg_sort( const void * a, const void * b ){
  // return (*(int*)a - *(int*)b)
  return arg_sort_list[*(int*)a] - arg_sort_list[*(int*)b];
}

/*
  Store an array of indices. 

  The constructor takes a pointer to an array of indices and takes
  ownership of the array. The BVecIndices object is responsible for
  freeing the memory. The pointer is set to NULL, but can be
  accessed with the getIndices function.
*/
BVecIndices::BVecIndices( int ** _indices, int _nindices ){
  indices = *_indices; 
  *_indices = NULL;
  index_args = NULL;
  nindices = _nindices;
  
  // Check if the indices are sorted and unique
  issorted = 1;
  for ( int i = 0; i < nindices-1; i++ ){
    if (indices[i+1] <= indices[i]){
      issorted = 0;
      break;
    }
    if (indices[i] < 0){
      fprintf(stderr, "BVecIndices: replacing negative index at \
entry %d with 0\n", i);
      indices[i] = 0;
    }
  }
}

BVecIndices::~BVecIndices(){
  if (indices){ delete [] indices; }
  if (index_args){ delete [] index_args; }
}

/*
  Retrieve the number of indices stored by this object, and their
  values.
*/
int BVecIndices::getNumIndices(){ return nindices; }

int BVecIndices::getIndices( int ** _indices ){
  *_indices = indices;
  return nindices;
}

/*
  Am I sorted? Or just a random array.
*/
int BVecIndices::isSorted(){ return issorted; }

/*
  Set up the inverse of the index set - if it isn't already.

  This function allocates and sets an array index_args[k] such that
  indices[index_args[k]] is sorted in ascending order.  This
  faciliates performing inverse look-ups. For instance, finding
  the index k such that indices[k] = var.  
*/
void BVecIndices::setUpInverse(){ 
  if (!index_args){ 
    index_args = new int[ nindices+1 ]; 
    for ( int k = 0; k < nindices; k++ ){ 
      index_args[k] = k; 
    }

    // Arg-sort it
    arg_sort_list = indices;
    qsort(index_args, nindices, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;
  }
}

/*
  Find the index such that indices[index] = var.

  If var is not found, or the inverse look-up has not been set up,
  return -1. Otherwise return the index >= 0.

  This function should run in approximately O(log(N)) where N are
  the number of indices.
*/
int BVecIndices::findIndex( int var ){
  const int n_linear_search = 10; // Start linear search when list is this long

  if (index_args){
    // Unforunately bsearch cannot be used in this situation, 
    // because var is the value, while the arg-sort implemented
    // above uses indirect addressing. This implements a binary 
    // search to achieve moderately okay times, I think.

    // If the number of indices is smaller than a certain value, 
    // just perform a linear search
    if (nindices < n_linear_search){
      for ( int k = 0; k < nindices; k++ ){
	int item = indices[k];
	if (item == var){
	  return k;
	}
      }
    }
    else {
      // Binary search an array to find k such that indices[k] = var,
      // where the array indices[index_args[k]] is sorted in ascending
      // order
      int high = nindices-1;
      int low = 0;
      int high_val = indices[index_args[high]];
      int low_val = indices[index_args[low]];

      // Check if the index is at the end points
      if (var == low_val){
	return index_args[low];
      }
      else if (var < low_val){
	return -1;
      }

      if (var == high_val){
	return index_args[high];
      }
      else if (var > high_val){
	return -1;
      }

      int mid = low + (int)((high - low)/2);
      
      // While there are indices left in the list
      while (low != mid){
	int mid_val = indices[index_args[mid]];
	if (mid_val == var){
	  return index_args[mid];
	}	

	if (var < mid_val){
	  high = mid;
	  high_val = mid_val;
	}
	else {
	  low = mid;
	  low_val = mid_val;
	}

	mid = low + (int)((high - low)/2);
      }      

      return -1;
    }
  }
  
  return -1;  
}

/*!
  The following are block-specific implementations of
  code required to copy over values from one array to another
*/
void VecDistGetVars( int bsize, int nvars, int * vars, int lower,
		     TacsScalar * x, TacsScalar * y, 
		     BVecDistribute::OpType op );
void VecDistSetVars( int bsize, int nvars, int * vars, int lower,
		     TacsScalar * x, TacsScalar * y, 
		     BVecDistribute::OpType op );

void VecDistGetVars1( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );
void VecDistSetVars1( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );

void VecDistGetVars2( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );
void VecDistSetVars2( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );

void VecDistGetVars3( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );
void VecDistSetVars3( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );

void VecDistGetVars5( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );
void VecDistSetVars5( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );

void VecDistGetVars6( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );
void VecDistSetVars6( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op );

/*! 
  Distribute vector components to other processors.

  This class is used to pass variables between processes 
  during parallel matrix-vector products.

  The user must provide the global indices of the required 
  variables. This class performs the following collection operation:

  for i = 1:nvars:
  local[i] = vec[ vars[i] ]

  where vars[i] are possibly non-local variables.

  Additionally, a reverse communication is also permitted where
  the following operation takes place,

  for i = 1:nvars:
  vec[ vars[i] ] += local[i]

  This operation is useful for assembling the residual equations
  for the finite--element method.

  The set up proceeds as follows:
  -------------------------------

  1. The ownership range of the input vector is determined.

  2. The break-points in the local vector are determined such that 
  for n in proc:
  for i in [buffRange[n],buffRange[n+1])
  varRange[n] < vars[i] < varRange[n+1] 

  All variables eminating from [buffRange[n],buffRange[n+1])
  For each process from which variables will be extracted, 
  in the array proc.

  3. Communicate the required number of variables between processes.
  (Create appropriate buffers on each recieving process of appropriate size)
  Communicate the block sizes back to the requesting processes. 
  (The buffer on the local/requesting process will be an offset into
  index into the local array supplied by the user).
*/
BVecDistribute::BVecDistribute( VarMap * _rmap, BVecIndices * _bindex ){
  rmap = _rmap;
  rmap->incref();

  // Retrieve the ownership data
  comm = rmap->getMPIComm();
  rmap->getOwnerRange(&ownerRange, &mpiRank, &mpiSize);

  // Set the default implementation
  bgetvars = VecDistGetVars;
  bsetvars = VecDistSetVars;

  // The external variables
  bindex = _bindex;  
  bindex->incref();
  int * vars;
  int nvars = bindex->getIndices(&vars);

  // Check whether the indices are sorted
  sorted_flag = bindex->isSorted();

  if (sorted_flag){
    // If the variables are sorted, do not allocate an array
    ext_vars = vars;
    next_vars = nvars;
    nvars_unsorted = 0;
    ext_unsorted_index = NULL;
    ext_sorted_vals = NULL;
  }
  else {
    // The following are only required when ext_vars is not initially sorted
    // input_vars[i] = ext_vars[ ext_unsorted_index[i] ]
    nvars_unsorted = nvars;

    ext_vars = new int[ nvars ];
    ext_unsorted_index = new int[ nvars_unsorted ];

    memcpy(ext_vars, vars, nvars_unsorted*sizeof(int));
    memcpy(ext_unsorted_index, vars, nvars_unsorted*sizeof(int));

    // Uniquely sort the array
    next_vars = FElibrary::uniqueSort( ext_vars, nvars );

    // For each value, go through and find the matching index
    for ( int i = 0; i < nvars_unsorted; i++ ){
      int * item = (int*) bsearch(&ext_unsorted_index[i], ext_vars, next_vars, 
                                  sizeof(int), FElibrary::comparator);

      ext_unsorted_index[i] = item - ext_vars;
    }
  }

  // Count up the number of processors that variables will be accquired from
  // Number of variables required by mpiRank from processor i
  int * ext_count = new int[ mpiSize ]; 
  ext_ptr  = new int[ mpiSize+1 ];
  ext_ptr[0] = 0;

  // Number of variables requested from mpiRank by processor i
  int * req_count = new int[ mpiSize ]; 
  req_ptr = new int[ mpiSize+1 ];
  req_ptr[0] = 0;

  FElibrary::matchIntervals(mpiSize, ownerRange, 
                            next_vars, ext_vars, ext_ptr);

  for ( int i = 0; i < mpiSize; i++ ){
    ext_count[i] = ext_ptr[i+1] - ext_ptr[i];
    if (i == mpiRank){
      ext_count[i] = 0;
    }
  }
  
  // Transmit the number of external variables to the other processes
  MPI_Alltoall(ext_count, 1, MPI_INT, req_count, 1, MPI_INT, comm);

  // Set up the index into the array that will be filled 
  // with the requesting variables
  req_ptr[0] = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    req_ptr[i+1] = req_ptr[i] + req_count[i];
  }

  req_vars = new int[ req_ptr[mpiSize] ];

  // Transmit the external variables to the owners
  MPI_Alltoallv(ext_vars, ext_count, ext_ptr, MPI_INT, 
		req_vars, req_count, req_ptr, MPI_INT, comm);

  delete [] ext_count;
  delete [] req_count;

  // Set the external buffer size and the request val size to zero
  extval_size = 0;
  reqval_size = 0;
  reqvals = NULL;
  ext_sorted_vals = NULL;  

  // Many processors will likely have no variables associated with them.
  // Create a data structure so that they are not looped over.
  int ne = 0, nr = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    if (i != mpiRank && ext_ptr[i+1] - ext_ptr[i] > 0){
      ne++;
    }
    if (req_ptr[i+1] - req_ptr[i] > 0){
      nr++;
    }
  }
  
  n_ext_proc = ne;
  n_req_proc = nr;
  ext_proc = new int[ ne ];
  req_proc = new int[ nr ];
  ne = nr = 0;

  for ( int i = 0; i < mpiSize; i++ ){
    if (i != mpiRank && ext_ptr[i+1] - ext_ptr[i] > 0){
      ext_proc[ne] = i;  ne++;
    }
    if (req_ptr[i+1] - req_ptr[i] > 0){
      req_proc[nr] = i; nr++;
    }
  }

  // Allocate space for the send/recevies
  if (n_req_proc > 0){
    sends = new MPI_Request[ n_req_proc ];
    send_status = new MPI_Status[ n_req_proc ];
  }
  else {
    sends = NULL;
    send_status = NULL;
  }

  if (n_ext_proc > 0){
    receives = new MPI_Request[ n_ext_proc ];
    receive_status = new MPI_Status[ n_ext_proc ];
  }
  else {
    receives = NULL;
    receive_status = NULL;
  }
}

/*
  Free the memory associated with the block vec distribute object
*/
BVecDistribute::~BVecDistribute(){
  delete [] ext_ptr;
  delete [] ext_proc;
  delete [] req_ptr;
  delete [] req_vars;
  delete [] req_proc;
  delete [] reqvals;

  if (!sorted_flag){
    delete [] ext_vars;
    delete [] ext_unsorted_index;
    delete [] ext_sorted_vals;
  }
  
  bindex->decref();
  rmap->decref();

  // Now, clean up the sends/receives
  if (n_req_proc > 0){
    delete [] sends;
    delete [] send_status;
  }
  
  if (n_ext_proc > 0){
    delete [] receives;
    delete [] receive_status;
  }
}

int BVecDistribute::getDim(){ return bindex->getNumIndices(); }

MPI_Comm BVecDistribute::getMPIComm(){ return comm; }

BVecIndices * BVecDistribute::getBVecIndices(){ return bindex; }

/*!
  Initiate the distribution of values from vec to the local array.
  
  The local array is copied to the ext_vars array. Both local and
  vec may be used immediately following the call to begin.

  for i in nvars:
  .   local[i] = vec[vars[i]]  

  Providing a non-zero varoffset modifies the result as follows:

  for i in nvars:
  .  local[i] = vec[vars[i] - varoffset]

  Note that for unsorted input the output is:
  local[i] = vec[ext_unsorted_index[i]]
*/
void BVecDistribute::beginForward( BVec *vec, TacsScalar *local, 
				   int varoffset ){
  // Copy the values into a request array -- this will be sorted over
  // each interval, but is not sorted between intervals
  TacsScalar *x;
  vec->getArray(&x);

  // Get the input vector block size - the local array is assumed to
  // match this block size
  int bsize = vec->getBlockSize();

  // Initialize the transfer operators
  initImpl(bsize);
    
  // Check if the request array is large enough to recieve the
  // incoming data. If not allocate a buffer that is large enough
  if (reqval_size < bsize*req_ptr[mpiSize]){
    reqval_size = bsize*req_ptr[mpiSize];
    if (reqvals){ delete [] reqvals; }
    reqvals = new TacsScalar[ reqval_size ];
  }

  // If the index array is not sorted, then allocate an array that is
  // large enough to store the sorted external values
  if (!sorted_flag){
    if (extval_size < bsize*ext_ptr[mpiSize]){
      extval_size = bsize*ext_ptr[mpiSize];
      if (ext_sorted_vals){ delete [] ext_sorted_vals; }
      ext_sorted_vals = new TacsScalar[ extval_size ];
    }
  }

  // Set the lower offset
  int lower = bsize*ownerRange[mpiRank];
  if (varoffset > 0){
    lower = lower + bsize*varoffset;
  }

  int tag = 0;

  // Do the off-process part
  bgetvars(bsize, req_ptr[mpiSize], req_vars, lower,
           x, reqvals, INSERT);

  for ( int i = 0; i < n_req_proc; i++ ){
    // Initiate the sends and receives
    int n = req_proc[i];
    int start = bsize*req_ptr[n];
    int size = bsize*req_ptr[n+1] - start;
    MPI_Isend(&reqvals[start], size, TACS_MPI_TYPE, n, tag, comm, &sends[i]);
  }

  if (sorted_flag){
    // First, copy over the local values - if any
    int start = ext_ptr[mpiRank];
    int size = ext_ptr[mpiRank+1] - start;

    bgetvars(bsize, size, &ext_vars[start], lower,
             x, &local[bsize*start], INSERT);

    // If the receiving array is sorted, it can be placed directly into local
    for ( int i = 0; i < n_ext_proc; i++ ){
      int n = ext_proc[i];
      if ( n != mpiRank ){
	int start = bsize*ext_ptr[n];
	int size = bsize*ext_ptr[n+1] - start;

	MPI_Irecv(&local[start], size, TACS_MPI_TYPE, n, tag, 
                  comm, &receives[i]);
      }
    }
  }
  else {
    // First, copy over the local values - if any
    int start = ext_ptr[mpiRank];
    int size = ext_ptr[mpiRank+1] - start;
    bgetvars(bsize, size, &ext_vars[start], lower,
             x, &ext_sorted_vals[bsize*start], INSERT);

    // If the receiving array is not sorted, the data must 
    // first be placed in a receiving array
    for ( int i = 0; i < n_ext_proc; i++ ){
      int n = ext_proc[i];
      if (n != mpiRank){
	int start = bsize*ext_ptr[n];
	int size = bsize*ext_ptr[n+1] - start;

	MPI_Irecv(&ext_sorted_vals[start], size, TACS_MPI_TYPE, n, tag, 
                  comm, &receives[i]);
      }
    }
  }
}

/*
  Finish the forward transfer of the data to the local vector
*/
void BVecDistribute::endForward( BVec *vec, TacsScalar *local ){
  // Get the input vector block size - the local array is assumed to
  // match this block size
  int bsize = vec->getBlockSize();

  // Finalize the sends and receives 
  MPI_Waitall(n_req_proc, sends, send_status);
  MPI_Waitall(n_ext_proc, receives, receive_status);  

  if (!sorted_flag){
    // Copy over the values from the sorted to the unsorted array
    bgetvars(bsize, nvars_unsorted, ext_unsorted_index, 0,
             ext_sorted_vals, local, INSERT);
  }
}

/*!
  Initiate the distribution of values from vec to the local array.
  
  The local array is copied to the ext_vars array. Both local and
  vec may be used immediately following the call to begin.

  for i in nvars:
  Vec[vars[i]] += local[i]
*/
void BVecDistribute::beginReverse( TacsScalar * local, BVec * vec, 
				   enum OpType op ){
  // Copy the values into a request array -- this will be sorted over 
  // each interval, but is not sorted between intervals
  TacsScalar * x;
  vec->getArray(&x);

  // Get the input vector block size - the local array is assumed to
  // match this block size
  int bsize = vec->getBlockSize();

  // Initialize the transfer operators
  initImpl(bsize);
    
  // Check if the request array is large enough to recieve the
  // incoming data. If not allocate a buffer that is large enough
  if (reqval_size < bsize*req_ptr[mpiSize]){
    reqval_size = bsize*req_ptr[mpiSize];
    if (reqvals){ delete [] reqvals; }
    reqvals = new TacsScalar[ reqval_size ];
  }

  // If the index array is not sorted, then allocate an array that is
  // large enough to store the sorted external values
  if (!sorted_flag){
    if (extval_size < bsize*ext_ptr[mpiSize]){
      extval_size = bsize*ext_ptr[mpiSize];
      if (ext_sorted_vals){ delete [] ext_sorted_vals; }
      ext_sorted_vals = new TacsScalar[ extval_size ];
    }

    // Copy over the values from the sorted to the unsorted array
    int len = bsize*ext_ptr[mpiSize];
    memset(ext_sorted_vals, 0, len*sizeof(TacsScalar));
    bsetvars(bsize, nvars_unsorted, ext_unsorted_index, 0,
             local, ext_sorted_vals, op);
    local = ext_sorted_vals;
  }

  int lower = bsize*ownerRange[mpiRank];

  // First do the on-process part
  int start = ext_ptr[mpiRank];
  int size = ext_ptr[mpiRank+1] - start;
  bsetvars(bsize, size, &ext_vars[start], lower,
           &local[bsize*start], x, op);
    
  // Note that the receives/send MPI_Requests have been reversed --
  // this is because they are allocated for the forward operation -
  // the sizes must be reversed for the reverse operation done here
  int tag = 1;
  for ( int i = 0; i < n_ext_proc; i++ ){
    int n = ext_proc[i];
    start = bsize*ext_ptr[n];
    size = bsize*ext_ptr[n+1] - start;
    MPI_Isend(&local[start], size, TACS_MPI_TYPE, n, tag, comm, &receives[i]);
  }
  
  for ( int i = 0; i < n_req_proc; i++ ){
    int n = req_proc[i];
    start = bsize*req_ptr[n];
    size = bsize*req_ptr[n+1] - start;
    MPI_Irecv(&reqvals[start], size, TACS_MPI_TYPE, n, tag, comm, &sends[i]);
  }  
}

/*
  End the reverse transfer of information from the local array to the
  vector using the supplied operation.
*/
void BVecDistribute::endReverse( TacsScalar * local, BVec * vec, 
				 enum OpType op ){
  // Get the input vector block size - the local array is assumed to
  // match this block size
  int bsize = vec->getBlockSize();

  // Finalize the sends and receives 
  MPI_Waitall(n_req_proc, sends, send_status);
  MPI_Waitall(n_ext_proc, receives, receive_status);  
 
  // Now go back through and add these values into their positions
  TacsScalar *x;
  vec->getArray(&x);

  int lower = bsize*ownerRange[mpiRank];

  bsetvars(bsize, req_ptr[mpiSize], req_vars, lower,
           reqvals, x, op);
}

const char * BVecDistribute::TACSObjectName(){
  return name;
}

const char * BVecDistribute::name = "BVecDistribute";

void BVecDistribute::initImpl( int bsize ){
  bgetvars = VecDistGetVars;
  bsetvars = VecDistSetVars;

  switch (bsize) {
  case 1:
    bgetvars = VecDistGetVars1;
    bsetvars = VecDistSetVars1;
    break;
  case 2:
    bgetvars = VecDistGetVars2;
    bsetvars = VecDistSetVars2;
    break;
  case 3:
    bgetvars = VecDistGetVars3;
    bsetvars = VecDistSetVars3;
    break;
  case 5:
    bgetvars = VecDistGetVars5;
    bsetvars = VecDistSetVars5;
    break;
  case 6:
    bgetvars = VecDistGetVars6;
    bsetvars = VecDistSetVars6;
    break;
  default:
    break;
  }
}

/*
  Block-specific implementations that should run slightly faster
*/

/*
  Copy variables distributed arbitrarily in one array to
  a second array where they are ordered sequentially.

  y[i] op x[var[i]]
*/
void VecDistGetVars( int bsize, int nvars, int * vars, int lower,
		     TacsScalar * x, TacsScalar * y, 
		     BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = bsize*vars[i] - lower;
      for ( int k = 0; k < bsize; k++ ){
	y[k] = x[v + k];
      }
 
      y += bsize;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = bsize * vars[i] - lower;
      for ( int k = 0; k < bsize; k++ ){
	y[k] += x[v + k];
      }          
       y += bsize;
     }    
  }
}

/*
  Copy variables from an ordered array to an array where
  they are distributed arbitrarily
  
  y[var[i]] op x[i]
*/
void VecDistSetVars( int bsize, int nvars, int * vars, int lower,
		     TacsScalar * x, TacsScalar * y, 
		     BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = bsize*vars[i] - lower;
      for ( int k = 0; k < bsize; k++ ){
	y[v + k] = x[k];
      }          
 
      x += bsize;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = bsize*vars[i] - lower;
      for ( int k = 0; k < bsize; k++ ){
	y[v + k] += x[k];
      }          
 
      x += bsize;
    }    
  }
}  

// ---------------
// Block size == 1
// ---------------
void VecDistGetVars1( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = vars[i] - lower;
      *y = x[v];
      y++;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = vars[i] - lower;
      *y += x[v];                
      y++;
    }    
  }
}

void VecDistSetVars1( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = vars[i] - lower;
      y[v] = *x;
      x++;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = vars[i] - lower;
      y[v] += *x;
      x++;
    }    
  }
}  

// ---------------
// Block size == 2
// ---------------
void VecDistGetVars2( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 2*vars[i] - lower;
      y[0] = x[v];
      y[1] = x[v+1];
      y += 2;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 2*vars[i] - lower;
      y[0] += x[v];      
      y[1] += x[v+1];      
      y += 2;
    }    
  }
}

void VecDistSetVars2( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 2*vars[i] - lower;
      y[v  ] = x[0];
      y[v+1] = x[1];
      x += 2;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 2*vars[i] - lower;
      y[v  ] += x[0];
      y[v+1] += x[1];
      x += 2;
    }    
  }
}  

// ---------------
// Block size == 3
// ---------------
void VecDistGetVars3( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 3*vars[i] - lower;
      y[0] = x[v];
      y[1] = x[v+1];
      y[2] = x[v+2];
      y += 3;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 3*vars[i] - lower;
      y[0] += x[v];      
      y[1] += x[v+1];      
      y[2] += x[v+2];      
      y += 3;
    }    
  }
}

void VecDistSetVars3( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 3*vars[i] - lower;
      y[v  ] = x[0];
      y[v+1] = x[1];
      y[v+2] = x[2];
      x += 3;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 3*vars[i] - lower;
      y[v  ] += x[0];
      y[v+1] += x[1];
      y[v+2] += x[2];
      x += 3;
    }    
  }
}  

// ---------------
// Block size == 5
// ---------------
void VecDistGetVars5( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 5*vars[i] - lower;
      y[0] = x[v];
      y[1] = x[v+1];
      y[2] = x[v+2];
      y[3] = x[v+3];
      y[4] = x[v+4];
      y += 5;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 5*vars[i] - lower;
      y[0] += x[v];      
      y[1] += x[v+1];      
      y[2] += x[v+2];      
      y[3] += x[v+3];      
      y[4] += x[v+4];      
      y += 5;
    }    
  }
}

void VecDistSetVars5( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 5*vars[i] - lower;
      y[v  ] = x[0];
      y[v+1] = x[1];
      y[v+2] = x[2];
      y[v+3] = x[3];
      y[v+4] = x[4];
      x += 5;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 5*vars[i] - lower;
      y[v  ] += x[0];
      y[v+1] += x[1];
      y[v+2] += x[2];
      y[v+3] += x[3];
      y[v+4] += x[4];
      x += 5;
    }    
  }
}  

// ---------------
// Block size == 6
// ---------------
void VecDistGetVars6( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 6*vars[i] - lower;
      y[0] = x[v];
      y[1] = x[v+1];
      y[2] = x[v+2];
      y[3] = x[v+3];
      y[4] = x[v+4];
      y[5] = x[v+5];
      y += 6;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 6*vars[i] - lower;
      y[0] += x[v];      
      y[1] += x[v+1];      
      y[2] += x[v+2];      
      y[3] += x[v+3];      
      y[4] += x[v+4];      
      y[5] += x[v+5];      
      y += 6;
    }    
  }
}

void VecDistSetVars6( int bsize, int nvars, int * vars, int lower,
		      TacsScalar * x, TacsScalar * y, 
		      BVecDistribute::OpType op ){
  if (op == BVecDistribute::INSERT){
    for ( int i = 0; i < nvars; i++ ){
      int v = 6*vars[i] - lower;
      y[v  ] = x[0];
      y[v+1] = x[1];
      y[v+2] = x[2];
      y[v+3] = x[3];
      y[v+4] = x[4];
      y[v+5] = x[5];
      x += 6;
    }
  }
  else { // Add
    for ( int i = 0; i < nvars; i++ ){
      int v = 6*vars[i] - lower;
      y[v  ] += x[0];
      y[v+1] += x[1];
      y[v+2] += x[2];
      y[v+3] += x[3];
      y[v+4] += x[4];
      y[v+5] += x[5];
      x += 6;
    }    
  }
}  
