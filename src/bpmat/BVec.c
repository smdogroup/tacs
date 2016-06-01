#include "BVec.h"
#include "FElibrary.h"
#include "tacslapack.h"

/*
  Code for the block-vector basis class

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*!
  VarMap
  
  Defines the variable map from the parallel distribution of variables
  to each process

  input:
  comm:  this object is defined over all processors in this comm
  N:     the number of nodes for this processor
*/
VarMap::VarMap( MPI_Comm _comm, int _N ){
  comm = _comm;
  MPI_Comm_rank(comm, &mpiRank);
  MPI_Comm_size(comm, &mpiSize);
  
  N = _N;
  // The ownership ranges for all processes
  ownerRange = new int[ mpiSize+1 ];
  memset(ownerRange, 0, (mpiSize+1)*sizeof(int));

  // Get the number of variables 
  ownerRange[0] = 0;
  MPI_Allgather(&N, 1, MPI_INT, &ownerRange[1], 1, MPI_INT, comm);
  
  for ( int i = 0; i < mpiSize; i++ ){    
    ownerRange[i+1] += ownerRange[i];
  }
}

VarMap::~VarMap(){
  delete [] ownerRange;
}

/*!
  BCMap class 

  Defines the Dirichlet boundary conditions for the vector
  and matrix classes.

  input:
  num_bcs:  an estimate of the number of boundary conditions
*/
BCMap::BCMap( int num_bcs ){
  num_bcs = (num_bcs >= 0 ? num_bcs : 0);
  max_size = num_bcs;
  // Usually, there are 6 or fewer dof per node
  max_var_ptr_size = 6*(num_bcs+1); 

  // Set the increment to be equal to the number of bcs set
  bc_increment = max_size+1;
  var_ptr_increment = max_var_ptr_size;

  nbcs = 0;
  global = new int[ max_size ];
  local = new int[ max_size ];
  var_ptr = new int[ max_size+1 ];
  vars = new int[ max_var_ptr_size ];  
  values = new TacsScalar[ max_var_ptr_size ];  
  var_ptr[0] = 0;
}

/*
  Delete all boundary condition information that this
  object allocated
*/
BCMap::~BCMap(){
  delete [] local;
  delete [] global;
  delete [] var_ptr;
  delete [] vars;  
  delete [] values;
}

/*
  Add a Dirichlet boundary condition with the specified local
  variable, global variable and the local Dirichlet BC number/value
  pair. Note, if no values are specified, they are assumed to be
  zero for each variable.

  input:
  local_var:  the local node number
  global_var: the global node number
  bc_nums:    the variable number to apply the BC
  bc_vals:    the value to apply
  nvals:      the number of values to apply at this node
*/
void BCMap::addBC( int local_var, int global_var, 
		   const int bc_nums[], const TacsScalar bc_vals[], 
		   int _nvals ){
  // If the number of boundary conditions exceeds the available
  // space, allocate more space and copy over the arrays
  if (nbcs+1 >= max_size){
    max_size = max_size + bc_increment;
    int * temp_local = new int[ max_size ];
    int * temp_global = new int[ max_size ];
    int * temp_ptr = new int[ max_size+1 ];

    for ( int i = 0; i < nbcs; i++ ){
      temp_local[i] = local[i];
      temp_global[i] = global[i];
      temp_ptr[i] = var_ptr[i];
    }
    temp_ptr[nbcs] = var_ptr[nbcs];

    delete [] local;
    delete [] global;
    delete [] var_ptr;
    local = temp_local;
    global = temp_global;
    var_ptr = temp_ptr;
  }

  // Set the new variable information
  local[nbcs] = local_var;
  global[nbcs] = global_var;
  var_ptr[nbcs+1] = var_ptr[nbcs] + _nvals;

  // If the number of boundary condition values/vars exceeds
  // the available space, allocate more space and copy over the
  // values from the old array
  if (var_ptr[nbcs+1]+1 >= max_var_ptr_size){
    max_var_ptr_size = var_ptr[nbcs+1] + var_ptr_increment;

    int * temp_vars = new int[ max_var_ptr_size ];
    TacsScalar * temp_vals = new TacsScalar[ max_var_ptr_size ];
    for ( int i = 0; i < var_ptr[nbcs]; i++ ){
      temp_vars[i] = vars[i];
      temp_vals[i] = values[i];
    }

    delete [] vars;
    delete [] values;
    vars = temp_vars;
    values = temp_vals;
  }

  if (bc_vals){
    // If values are specified, copy over the values
    for ( int i = var_ptr[nbcs], k = 0; i < var_ptr[nbcs+1]; i++, k++ ){
      vars[i] = bc_nums[k];
      values[i] = bc_vals[k];
    }  
  }
  else {
    // If no values are specified, assume the values are zero
    for ( int i = var_ptr[nbcs], k = 0; i < var_ptr[nbcs+1]; i++, k++ ){
      vars[i] = bc_nums[k];
      values[i] = 0.0;
    }  
  }
  nbcs++;
}

/*
  Retrieve the boundary conditions that have been set locally
  within this object

  output:
  local_vars:  the local node numbers
  global_vars: the global node numbers
  var_ptr:     the pointer into the list of nodal vars/values
  vars:        the variable number to apply the BC
  values:      the value to apply
*/
int BCMap::getBCs( const int ** _local, const int ** _global,
		   const int ** _var_ptr, 
		   const int ** _vars, const TacsScalar ** _values ){
  *_local = local;
  *_global = global;
  *_var_ptr = var_ptr;
  *_vars = vars;
  *_values = values;

  return nbcs;
}

/*!
  Create a block-based parallel vector

  input:
  rmap:   the variable->processor map for the unknowns
  bcs:    the boundary conditions associated with this vector
*/
BVec::BVec( VarMap * _rmap, int _bsize, BCMap * _bcs ){
  rmap = _rmap;
  rmap->incref();

  // Copy the boundary conditions
  bcs = _bcs;
  if (bcs){ bcs->incref(); }

  bsize = _bsize;
  comm = rmap->getMPIComm();
  rmap->getOwnerRange(&ownerRange, &mpiRank, &mpiSize);
  size = bsize*rmap->getDim();

  displaced = NULL;
  x = new TacsScalar[ size ];
  memset(x, 0, size*sizeof(TacsScalar));
}

/*
  Create the block-based parallel vector without a VarMap object
  or boundary conditions - this is required for some parallel 
  matrix objects, or other situtations in which neither object
  exists.

  input:
  comm:  the communicator that this object is distributed over
  size:  the number of local entries stored by this vector
*/
BVec::BVec( MPI_Comm _comm, int _size, int _bsize ){
  bsize = _bsize;
  size = _size;
  comm = _comm;
  rmap = NULL;
  bcs = NULL;
  ownerRange = NULL;
  displaced = NULL;
  x = new TacsScalar[ size ];
  memset(x, 0, size*sizeof(TacsScalar));
}

BVec::~BVec(){
  if (rmap){ rmap->decref(); }
  if (bcs){  bcs->decref();  }

  if (x){ delete [] x; }
  if (displaced){ delete [] displaced; }
}

/*
  Get the local size of the vector on this processor
*/
void BVec::getSize( int * _size ){
  *_size = size;
}

/*
  Compute the norm of the vector
*/
TacsScalar BVec::norm(){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::scale Error, no array\n", mpiRank);
    return (0.0);
  }

  // Compute the norm for each processor
  TacsScalar res, sum;
#if defined(TACS_USE_COMPLEX) 
  res = TacsScalar(0.0);
  int i = 0;
  int rem = size%4;
  TacsScalar * y = x;
  for ( ; i < rem; i++ ){
    res += y[0]*y[0];
    y++;
  }
  
  for ( ; i < size; i += 4 ){
    res += y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3];
    y += 4;
  }
#else
  int one = 1;
  res = BLASnrm2(&size, x, &one);
  res *= res;
#endif
  TacsAddFlops(2*size);

  MPI_Allreduce(&res, &sum, 1, TACS_MPI_TYPE, MPI_SUM, comm);

  return sqrt(sum);
}

/*
  Scale the vector by a scalar
*/
void BVec::scale( TacsScalar alpha ){        
  if (!x){ 
    fprintf(stderr, "[%d] BVec::scale Error, no array\n", mpiRank);
    return;
  }

  int one = 1;
  BLASscal(&size, &alpha, x, &one);

  TacsAddFlops(size);
}

/*
  Compute the dot product of two vectors
*/
TacsScalar BVec::dot( TACSVec * tvec ){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::dot Error, no array\n", mpiRank);
    return TacsScalar(0.0);
  }

  TacsScalar sum = 0.0;
  BVec * vec = dynamic_cast<BVec*>(tvec);
  if (vec){
    if (vec->size != size){
      fprintf(stderr, "[%d] BVec::dot Error, the sizes must be \
the same\n", mpiRank);
      return TacsScalar(0.0);
    }
    
    TacsScalar res;
#if defined(TACS_USE_COMPLEX) 
    res = TacsScalar(0.0);
    int i = 0;
    int rem = size%4;
    TacsScalar * y = x;
    TacsScalar * z = vec->x;
    for ( ; i < rem; i++ ){
      res += y[0]*z[0];
      y++; z++;
    }
    
    for ( ; i < size; i += 4 ){
      res += y[0]*z[0] + y[1]*z[1] + y[2]*z[2] + y[3]*z[3];
      y += 4;
      z += 4;
    }
#else
    int one = 1;
    res = BLASdot(&size, x, &one, vec->x, &one);
#endif
    MPI_Allreduce(&res, &sum, 1, TACS_MPI_TYPE, MPI_SUM, comm);
  }
  else {
    fprintf(stderr, "BVec type error: Input must be BVec\n");
  }

  TacsAddFlops(2*size);

  return sum;
}

/*
  Compute multiple dot products. This is more efficient for parallel
  computations since there are fewer gather operations for the same
  number of dot products.
*/
void BVec::mdot( TACSVec ** tvec, TacsScalar * ans, int nvecs ){
  if (!x){ 
    fprintf( stderr, "[%d] BVec::dot Error, no array\n", mpiRank );
    return;
  }

  for ( int k = 0; k < nvecs; k++ ){
    ans[k] = 0.0;

    BVec * vec = dynamic_cast<BVec*>(tvec[k]);
    if (vec){
      if (vec->size != size){
        fprintf(stderr, "[%d] BVec::dot Error, the sizes must be \
the same\n", mpiRank);    
        continue;
      }

#if defined(TACS_USE_COMPLEX) 
      TacsScalar res = TacsScalar(0.0);
      int i = 0;
      int rem = size % 4;
      TacsScalar * y = x;
      TacsScalar * z = vec->x;
      for ( ; i < rem; i++ ){
        res += y[0]*z[0];
        y++; z++;
      }
      
      for ( ; i < size; i += 4 ){
        res += y[0]*z[0] + y[1]*z[1] + y[2]*z[2] + y[3]*z[3];
        y += 4;
        z += 4;
      }
      
      ans[k] = res;
#else    
      int one = 1;
      ans[k] = BLASdot(&size, x, &one, vec->x, &one);
#endif
    }
    else {
      fprintf(stderr, "BVec type error: Input must be BVec\n");
    }
  }

  TacsAddFlops(2*nvecs*size);
  
  MPI_Allreduce(MPI_IN_PLACE, ans, nvecs, TACS_MPI_TYPE, MPI_SUM, comm);
}

/*
  Compute y = alpha*x + y
*/
void BVec::axpy( TacsScalar alpha, TACSVec * tvec ){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::axpy Error, no array\n", mpiRank);
    return;
  }
  
  BVec * vec = dynamic_cast<BVec*>(tvec);

  if (vec){
    if (vec->size != size){
      fprintf(stderr, "[%d] BVec::axpy Error, the sizes must be \
the same\n", mpiRank);
      return;
    }
    
    int one = 1;
    BLASaxpy(&size, &alpha, vec->x, &one, x, &one);
  }
  else {
    fprintf(stderr, "BVec type error: Input must be BVec\n");
  }

  TacsAddFlops(2*size);
}

/*
  Compute x <- alpha * vec + beta * x
*/
void BVec::axpby( TacsScalar alpha, TacsScalar beta, TACSVec * tvec ){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::axpby Error no array\n", mpiRank);
    return;
  }

  BVec * vec = dynamic_cast<BVec*>(tvec);

  if (vec){
    if (vec->size != size){
      fprintf(stderr, "[%d] BVec::axpby Error sizes must be \
the same\n", mpiRank);
      return;
    }

    int i = 0;
    int rem = size % 4;
    
    TacsScalar * y = x;
    TacsScalar * z = vec->x;
    
    for ( ; i < rem; i++ ){
      y[0] = beta*y[0] + alpha*z[0];
      y++; z++;
    }
    
    for ( ; i < size; i += 4 ){
      y[0] = beta*y[0] + alpha*z[0];
      y[1] = beta*y[1] + alpha*z[1];
      y[2] = beta*y[2] + alpha*z[2];
      y[3] = beta*y[3] + alpha*z[3];
      
      y += 4;
      z += 4;
    }
  }
  else {
    fprintf(stderr, "BVec type error: Input must be BVec\n");
  }

  TacsAddFlops(3*size);
}

/*
  Copy the values x <- vec->x
*/
void BVec::copyValues( TACSVec * tvec ){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::scale Error, no array\n", mpiRank);
    return;
  }
  
  BVec * vec = dynamic_cast<BVec*>(tvec);
  if (vec){
    if (vec->size != size){
      fprintf(stderr, "[%d] BVec::copyValues error, sizes must be \
the same\n", mpiRank);
      return;
    }
    
    int one = 1;
    BLAScopy(&size, vec->x, &one, x, &one);
  }
  else {
    fprintf(stderr, "BVec type error: Input must be BVec\n");
  }
}

/*
  Zero all the entries in the vector
*/
void BVec::zeroEntries(){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::zeroEntries Error, no array\n", mpiRank);
    return;
  }

  memset(x, 0, size*sizeof(TacsScalar));
}

/*
  Set all the entries in the vector to val
*/
void BVec::set( TacsScalar val ){
  if (!x){ 
    fprintf(stderr, "[%d] BVec::set Error, no array\n", mpiRank);
    return;
  }

  int i = 0;
  int rem = size%4;
  TacsScalar * y = x;

  for ( ; i < rem; i++ ){
    y[0] = val;
    y++;
  }

  for ( ; i < size; i += 4 ){
    y[0] = y[1] = y[2] = y[3] = val;
    y += 4;
  } 
}

/*
  Initialize the random value generator 
*/
void BVec::initRand(){
  unsigned int t = time(NULL);
  MPI_Bcast(&t, 1, MPI_INT, 0, comm);
  srand(t);
}

/*
  Set all the values in the vector using a uniform pseudo-random
  distribution over an interval between lower/upper.
*/  
void BVec::setRand( double lower, double upper ){
  if (!x ){ 
    fprintf(stderr, "[%d] BVec::setRand Error, no array\n", mpiRank);
    return;
  }

  if (!ownerRange){
    for ( int i = 0; i < size; i++ ){
      x[i] = lower + ((upper - lower)*rand())/(1.0*RAND_MAX);
    }
  }
  else {
    // Generate random values for each processor sequentially.
    // This should result in the same number of calls on all
    // processors.
    for ( int k = 0; k < mpiSize; k++ ){
      if (k != mpiRank){
	int end = bsize*ownerRange[k+1];
	for ( int i = bsize*ownerRange[k]; i < end; i++ ){
	  rand(); 
	}
      }
      else {
	for ( int i = 0; i < size; i++ ){
	  x[i] = lower + ((upper - lower)*rand())/(1.0*RAND_MAX);
	}
      }    
    }
  }
}

/*
  Set an array into this vector
*/
void BVec::setArray(TacsScalar * _x){
  x = _x;
}

/*
  Place an array into the place of the local storage for this vector
*/
void BVec::placeArray( TacsScalar * _x ){
  if (!displaced){
    displaced = x;
    x = _x;
  }
  else {
    fprintf(stderr, "[%d] Bvec::placeArray Error, already displaced \
array\n", mpiRank);
  }
}

/*
  Restore the array displaced by the call to placeArray() with the
  original arrayremove the 'placed' array with the original array 
*/
void BVec::restoreArray(){
  if (displaced){
    x = displaced;
    displaced = NULL;
  }
  else {
    fprintf(stderr, "[%d] Bvec::restoreArray Error, no array to restore\n", 
            mpiRank);
  }
}

/*
  Retrieve the locally stored values from the array
*/
int BVec::getArray( TacsScalar ** array ){
  *array = x;
  return size;
}

/*
  Apply the Dirichlet boundary conditions to the vector
*/
void BVec::applyBCs(){
  // apply the boundary conditions
  if (bcs && x){
    const int *local, *global, *var_ptr, *vars;
    const TacsScalar *values;
    int nbcs = bcs->getBCs(&local, &global, &var_ptr, &vars, &values);

    for ( int i = 0; i < nbcs; i++ ){
      if ( global[i] >= ownerRange[mpiRank] &&
	   global[i] <  ownerRange[mpiRank+1] ){
	int var = bsize*(global[i] - ownerRange[mpiRank]);
	
	for ( int k = var_ptr[i]; k < var_ptr[i+1]; k++ ){ 
	  // Scan through the rows to be zeroed
	  x[var + vars[k]] = TacsScalar(0.0);      
	}
      }
    }
  }
}

/*
  Copy values from vec to this. This can be done without parallel 
  communication.
*/
void BVec::copySeqValues( BVec * vec ){
  if (!x){ 
    fprintf( stderr, "[%d] BVec::copySeqValues Error, no array\n", mpiRank );
    return;
  }
  else if (vec->size != ownerRange[mpiSize]){
    fprintf( stderr, "[%d] BVec::copySeqValues Error, vectors must be the \
same size\n", mpiRank );
    return;
  }

  int end   = ownerRange[mpiRank+1];
  int start = ownerRange[mpiRank];
  for ( int i = start; i < end; i++ ){
    x[i - start] = vec->x[i];
  }
}
  
/*!
  Set values from this into vec. This must be done with MPI_Allgather. 
*/
void BVec::setValuesSeq( BVec * vec ){
  if (!x){ 
    fprintf( stderr, "[%d] BVec::setValuesSeq Error, no array\n", mpiRank );
    return;
  }
  else if (vec->size != ownerRange[mpiSize]){
    fprintf( stderr, "[%d] BVec::setValuesSeq: error, vectors must be \
the same size\n", mpiRank );
    return;
  }

  int * count = new int[ mpiSize ];
  int * ptr   = new int[ mpiSize ];
  for ( int i = 0; i < mpiSize; i++ ){
    ptr[i] = ownerRange[i];
    count[i] = ownerRange[i+1] - ownerRange[i];
  }

  MPI_Allgatherv(x, size, TACS_MPI_TYPE, 
                 vec->x, count, ptr, TACS_MPI_TYPE, comm);

  delete [] ptr;
  delete [] count;
}

const char * BVec::TACSObjectName(){
  return vecName;
}

/*!
  Write the values to a file.

  This uses MPI file I/O. The filenames must be the same on all
  processors. The format is independent of the number of processors.

  The file format is as follows:
  int                       The length of the vector
  len * sizeof(TacsScalar)  The vector entries
*/
int BVec::writeToFile( const char * filename ){
  char * fname = new char[ strlen(filename)+1 ];
  strcpy(fname, filename);

  int * range = new int[ mpiSize+1 ];
  range[0] = 0;
  MPI_Allgather(&size, 1, MPI_INT, &range[1], 1, MPI_INT, comm);
  for ( int i = 0; i < mpiSize; i++ ){
    range[i+1] += range[i];
  }

  int fail = 0;
  MPI_File fp = NULL;
  MPI_File_open(comm, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, 
                MPI_INFO_NULL, &fp);

  if (fp){
    if ( mpiRank == 0 ){
      MPI_File_write(fp, &range[mpiSize], 1, MPI_INT, MPI_STATUS_IGNORE);
    }

    char datarep[] = "native";
    MPI_File_set_view(fp, sizeof(int), TACS_MPI_TYPE, TACS_MPI_TYPE, 
                      datarep, MPI_INFO_NULL);
    MPI_File_write_at_all(fp, range[mpiRank], x, size, TACS_MPI_TYPE, 
                          MPI_STATUS_IGNORE);
    MPI_File_close(&fp);
  }
  else {
    fail = 1;
  }

  delete [] range;
  delete [] fname;

  return fail;
}

/*!
  Read values from a binary data file.

  The size of this vector must be the size of the vector
  originally stored in the file otherwise nothing is read in.

  The file format is as follows:
  int                       The length of the vector
  len * sizeof(TacsScalar)  The vector entries
*/
  
int BVec::readFromFile( const char * filename ){
  char * fname = new char[ strlen(filename)+1 ];
  strcpy(fname, filename);

  int * range = new int[ mpiSize+1 ];
  range[0] = 0;
  MPI_Allgather(&size, 1, MPI_INT, &range[1], 1, MPI_INT, comm);
  for ( int i = 0; i < mpiSize; i++ ){
    range[i+1] += range[i];
  }

  int fail = 0;
  MPI_File fp = NULL;
  MPI_File_open(comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);

  if (fp){
    int len = 0;
    if ( mpiRank == 0 ){
      MPI_File_read(fp, &len, 1, MPI_INT, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&len, 1, MPI_INT, 0, comm);
    if (len != range[mpiSize]){
      fprintf( stderr, "[%d] Cannot read BVec from file, incorrect size \
 %d != %d \n", mpiRank, range[mpiSize], len );
      memset(x, 0, size*sizeof(TacsScalar));
    }

    char datarep[] = "native";
    MPI_File_set_view(fp, sizeof(int), TACS_MPI_TYPE, TACS_MPI_TYPE, 
                      datarep, MPI_INFO_NULL);
    MPI_File_read_at_all(fp, range[mpiRank], x, size, TACS_MPI_TYPE, 
                         MPI_STATUS_IGNORE);
    MPI_File_close(&fp);
  }
  else {
    fail = 1;
  }

  delete [] range;
  delete [] fname;

  return fail;
}

const char * BVec::vecName = "BVec";
