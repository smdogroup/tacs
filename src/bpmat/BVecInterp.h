#ifndef BVEC_INTERP_H
#define BVEC_INTERP_H

/*
  Interpolate between two BVecs

  Copyright (c) 2014 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

#include "BVec.h"
#include "BVecDist.h"

/*
  BVecInterp: Interpolate with constant weights between two vectors
  that are not the same size. This applies the same weights to all
  components of the blocks in each BVec object.

  This interpolation class is limited in the following sense:

  1. The output from the forward operation (mult/multAdd) must only
  produce results on this processor.

  2. The input from the forward operation (mult/multAdd) can be from
  any processor. 

  The operation of the BVecInterp object is essentially equivalent to
  a matrix-vector product with a non-rectangular matrix. However, the
  implementation is more efficient than the block-matrix methods
  implemented elsewhere in TACS, since each block would be a bsize x
  bsize identity matrix.

  The BVecInterp operator performs the following multiplication:

  1. y <- A*x
  2. y <- A*x + z

  or the following transpose operations:

  3. y <- A^{T}*x
  4. y <- A^{T}*x + z

  This class is used extensively in the TACS implementation of
  multigrid.
*/
class TACSBVecInterp : public TACSObject {
 public:
  TACSBVecInterp( TACSVarMap *_inMap, TACSVarMap *_outMap, int _bsize );
  ~TACSBVecInterp();

  // Add components of the interpolation
  // -----------------------------------
  void addInterp( int vNum, TacsScalar weights[], int inNums[], int size );
  void initialize();

  // Perform the foward interpolation
  // --------------------------------
  void mult( TACSBVec *in, TACSBVec *out );
  void multAdd( TACSBVec *in, TACSBVec *add, TACSBVec *out );

  // Perform the transpose interpolation
  // -----------------------------------
  void multTranspose( TACSBVec *in, TACSBVec *out );
  void multTransposeAdd( TACSBVec *in, TACSBVec *add, TACSBVec *out );

  // Perform the weighted transpose interpolation
  // --------------------------------------------
  void multWeightTranspose( TACSBVec *in, TACSBVec *out );

  // Print the interpolation
  // -----------------------
  void printInterp( const char *filename );

 private:
  // The MPI communicator
  MPI_Comm comm;

  void (*multadd)( int bsize, int nrows, 
		   const int *rowp, const int *cols,
		   const TacsScalar *weights,
		   const TacsScalar *x, TacsScalar *y );
  void (*multtransadd)( int bsize, int nrows, 
			const int *rowp, const int *cols,
			const TacsScalar *weights,
			const TacsScalar *x, TacsScalar *y );

  // The on and off-processor parts of the interpolation
  // These are dynamically expanded if they are not large enough
  int max_on_size, max_on_weights;
  int on_size, *on_nums, *on_rowp, *on_vars;
  TacsScalar *on_weights;

  int max_off_size, max_off_weights;
  int off_size, *off_nums, *off_rowp, *off_vars;
  TacsScalar *off_weights;

  // The local weight contributions
  int *rowp, *cols;
  TacsScalar *weights;

  // The external weight contributions
  int *ext_rowp, *ext_cols;
  TacsScalar *ext_weights;

  int num_ext_vars; // The number of external variables
  TacsScalar *x_ext; // Variable values from other processors

  // The number of local rows from outMap
  int N, M, bsize;
  TACSVarMap *inMap, *outMap;

  // Store the transpose weights
  TacsScalar *transpose_weights;

  // The object responsible for fetching/distributing the 
  // external variables
  TACSBVecDistribute *vecDist;
  TACSBVecDistCtx *ctx;
};

#endif
