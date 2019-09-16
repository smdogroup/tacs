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

#ifndef TACS_DISTRIBUTED_MAT_H
#define TACS_DISTRIBUTED_MAT_H

#include "TACSBVecDistribute.h"
#include "TACSParallelMat.h"

class TACSMatDistCtx {
 public:



  // Functions for setting values in the matrix
  // ------------------------------------------
  void addValues( int nrow, const int *row,
                  int ncol, const int *col,
                  int nv, int mv, const TacsScalar *values );
  void addWeightValues( int nvars, const int *varp, const int *vars,
                        const TacsScalar *weights,
                        int nv, int mv, const TacsScalar *values,
                        MatrixOrientation matOr=NORMAL );

  // Information for the persistent communication set up
  // ---------------------------------------------------
  int nsends, nreceives;
  MPI_Request *sends, *receives;
  MPI_Status *send_status, *receive_status;
  int *send_proc, *receive_proc;



  // Pointer to the external data
  TacsScalar *ext_A;

  // Pointer to incoming data
  TacsScalar *in_A;
};


/*
  Distribute components of a matrix from a local CSR format to a global
  distributed format compatible with TACSParallelMat.

  This distribution may be performed in essentially two different ways:

  First, the values may be set
  on the desired

*/
class TACSMatDistribute : public TACSObject {
 public:
  TACSMatDistribute( TACSNodeMap *rmap, int num_nodes,
                     const int *rowp, const int *cols,
                     TACSBVecIndices *bindex );
  ~TACSMatDistribute();

  void beginAssembly();
  void endAssembly();

  // Set values into the matrix from the local BCSRMat
  // -------------------------------------------------
  void setValues( int nvars, const int *ext_vars,
                  const int *rowp, const int *cols, TacsScalar *avals );

 private:
  // Set up the local/external CSR data structure
  void setUpLocalExtCSR( int next_vars, const int *ext_vars,
                         const int *rowp, const int *cols,
                         int lower, int upper,
                         int nz_per_row,
                         int ** _Arowp, int ** _Acols,
                         int *_Np, int ** _Browp, int ** _Bcols );

  // Initialize the persistent communication requests
  void initPersistent();

  // Variables for the column map
  // ----------------------------
  int col_map_size;
  const int *col_map_vars;

  // Information about assembling the non-zero pattern
  // -------------------------------------------------
  MPI_Comm comm;

  // Data destined for other processes
  // ---------------------------------
  int next_rows;
  int *ext_rows;
  int *ext_row_ptr, *ext_row_count;
  int *ext_rowp;
  int *ext_cols;

  // Data received from other processes
  // ----------------------------------
  int nin_rows;
  int *in_rows;
  int *in_row_ptr, *in_row_count;
  int *in_rowp;
  int *in_cols;
};

#endif // TACS_DISTRIBUTED_MAT_H
