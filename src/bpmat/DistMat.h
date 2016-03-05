#ifndef TACS_DISTRIBUTED_MAT_H
#define TACS_DISTRIBUTED_MAT_H

/*
  Distributed matrix implementation

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "PMat.h"

/*
  Distributed matrix based on a parallel distribution of rows
  (equations).  Inherits from PMat.
*/
class DistMat : public PMat {
 public:
  DistMat( TACSThreadInfo * thread_info, VarMap * _rmap, 
	   int next_vars, const int * rowp, const int * cols, 
	   BVecIndices * bindex,
	   BCMap * _bcs = NULL );
  ~DistMat();
    
  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();
  void addValues( int nrow, const int * row, 
		  int ncol, const int * col,  
		  int nv, int mv, const TacsScalar * values ); 
  void addWeightValues( int nvars, const int *varp, const int *vars,
			const TacsScalar *weights,
			int nv, int mv, const TacsScalar *values );
  void beginAssembly();
  void endAssembly();
  
  // Set values into the matrix from the local BCSRMat
  // -------------------------------------------------
  void setValues( int nvars, const int * ext_vars,
		  const int * rowp, const int * cols, TacsScalar * avals );

 private:
  // Set up the local/external CSR data structure
  void setUpLocalExtCSR( int next_vars, int * ext_vars, 
			 const int * rowp, const int * cols,
			 int lower, int upper,
			 int nz_per_row,
			 int ** _Arowp, int ** _Acols,
			 int * _Np, int ** _Browp, int ** _Bcols );

  // Initialize the persistent communication requests 
  void initPersistent();
   
  // Variables for the column map
  // ----------------------------
  int col_map_size;
  int * col_map_vars;
 
  // Information about assembling the non-zero pattern
  // -------------------------------------------------
  MPI_Comm comm;
  const int * ownerRange;
  int mpiSize, mpiRank;
  
  // Data destined for other processes
  // ---------------------------------
  int next_rows;
  int *ext_rows;
  int *ext_row_ptr, *ext_row_count;
  int *ext_rowp;
  int *ext_cols;

  TacsScalar * ext_A; // Pointer to the external data
  
  // Data received from other processes
  // ----------------------------------
  int nin_rows;
  int *in_rows;
  int *in_row_ptr, *in_row_count;
  int *in_rowp;
  int *in_cols;

  TacsScalar * in_A; // Pointer to incoming data

  // Information for the persistent communication set up
  // ---------------------------------------------------
  int nsends, nreceives;
  MPI_Request *sends, *receives;
  MPI_Status *send_status, *receive_status;
  int *send_proc, *receive_proc;
};

#endif
