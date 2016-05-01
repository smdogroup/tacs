#ifndef TACS_BVEC_DISTRIBUTE_H
#define TACS_BVEC_DISTRIBUTE_H

/*
  Distribute/gather elements of a BVec

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "BVec.h"

/*
  A class containing a pointer to an array of indices.  These indices
  may be either sorted in asscending order and unique or completely
  arbitrary. These cases are detected automatically by the index
  class. The indices may not be negative.  Negative indices are
  replaced by 0, so be warned!!

  This object takes ownership of the array that is passed in.
*/
class BVecIndices : public TACSObject {
 public:
  BVecIndices( int **_indices, int _nindices );
  ~BVecIndices();

  // Retrieve information about the indices
  // --------------------------------------
  int getNumIndices();
  int getIndices( int **_indices );
  int isSorted();

  // Set up/use an arg-sorted array to find the reverse
  // map to find k such that indices[k] = var 
  // --------------------------------------------------
  void setUpInverse();
  int findIndex( int index ); // Find the index k such that indices[k] = index

 private:
  int *indices;
  int nindices;
  int issorted;
  int *index_args;
};

/*!
  Distribute vector components to other processors and collect 
  contributions from other processors.

  This class is used to pass external interface variables between
  processes during parallel matrix-vector products.

  This class performs the following operation:

  for i = 1:nvars:
  local[i] = vec[vars[i]]

  where vars[i] are possibly non-local variables.

  Additionally, a reverse communication is also permitted where the
  following operation takes place,

  for i = 1:nvars:
  vec[vars[i]] += local[i]

  This operation is useful for assembling the residual equations
  within the finite--element method.
*/
class BVecDistribute : public TACSObject {
 public:
  enum OpType { INSERT, ADD };
  BVecDistribute( VarMap *rmap, BVecIndices *bindex );
  ~BVecDistribute();

  // Get the size of the local array
  // All arrays passed must be at least this size
  // --------------------------------------------
  int getDim();
  BVecIndices *getBVecIndices();

  // Transfer the data to the array provided
  // ---------------------------------------
  void beginForward( BVec *vec, TacsScalar *local, 
                     int var_offset = 0 ); 
  void endForward( BVec *vec, TacsScalar *local );

  // Add or insert data back into the vector
  // ---------------------------------------
  void beginReverse( TacsScalar *local, BVec *vec, enum OpType op=ADD ); 
  void endReverse( TacsScalar *local, BVec *vec, enum OpType op=ADD );    

  MPI_Comm getMPIComm();
  const char *TACSObjectName();

 private:
  // Block-specific implementation pointers
  // --------------------------------------
  void initImpl( int bsize );
  void (*bgetvars)( int bsize, int nvars, int * vars, int lower,
		    TacsScalar * x, TacsScalar * y, 
		    BVecDistribute::OpType op );
  void (*bsetvars)( int bsize, int nvars, int * vars, int lower,
		    TacsScalar * x, TacsScalar * y, 
		    BVecDistribute::OpType op );

  // Data defining the distribution of the variables
  VarMap * rmap;

  // The communicator and the MPI data
  MPI_Comm comm;
  int mpiRank, mpiSize;
  const int * ownerRange;

  // Object containing the indices of the external variables
  // -------------------------------------------------------
  BVecIndices * bindex;

  // Data for dealing with an array that is not sorted or unique
  // -----------------------------------------------------------
  int sorted_flag;
  int nvars_unsorted;
  int *ext_unsorted_index;
  TacsScalar *ext_sorted_vals;

  // Data for collecting external variables
  // --------------------------------------
  int next_vars;
  int *ext_ptr;  // Displacements into the local external array
  int *ext_vars; // External variables that are requested by this process
  int extval_size;

  // Data for the requested values
  int *req_ptr;  // Displacement into requested array
  int *req_vars; // Variables that have been requested
  int reqval_size;
  TacsScalar *reqvals;

  // Sending data
  int n_req_proc; // Processes to send non-zero mesages to
  int *req_proc; 
  MPI_Request *sends;
  MPI_Status *send_status;

  // Receiving data
  int n_ext_proc; // Externall processes to expect non-zero receives from
  int *ext_proc;

  MPI_Request *receives;
  MPI_Status *receive_status;

  static const char *name;
};

#endif
