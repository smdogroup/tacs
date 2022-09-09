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

#ifndef TACS_MAT_DISTRIBUTE_H
#define TACS_MAT_DISTRIBUTE_H

class TACSParallelMat;

#include "BCSRMat.h"
#include "TACSBVecDistribute.h"
#include "TACSParallelMat.h"

/*
  Distribute components of a matrix from a local CSR format to a global
  distributed format compatible with TACSParallelMat.
*/
class TACSMatDistribute : public TACSObject {
 public:
  TACSMatDistribute(TACSThreadInfo *thread_info, TACSNodeMap *rmap, int bsize,
                    int num_nodes, const int *rowp, const int *cols,
                    TACSBVecIndices *bindex, BCSRMat **_Aloc, BCSRMat **_Bext,
                    TACSBVecDistribute **_colMap);
  ~TACSMatDistribute();

  // Zero the entries in the temporary storage
  // -----------------------------------------
  void zeroEntries();

  // Add values to the entries in the matrix
  // ---------------------------------------
  void addValues(TACSParallelMat *mat, int nrow, const int *row, int ncol,
                 const int *col, int nv, int mv, const TacsScalar *values);
  void addWeightValues(TACSParallelMat *mat, int nvars, const int *varp,
                       const int *vars, const TacsScalar *weights, int nv,
                       int mv, const TacsScalar *values,
                       MatrixOrientation matOr = TACS_MAT_NORMAL);

  // Set values into the matrix from the local BCSRMat
  // -------------------------------------------------
  void setValues(TACSParallelMat *mat, int nvars, const int *ext_vars,
                 const int *rowp, const int *cols, TacsScalar *avals);

  // Begin/end distributing values into the matrix
  // ---------------------------------------------
  void beginAssembly(TACSParallelMat *mat);
  void endAssembly(TACSParallelMat *mat);

 private:
  // Set up the local/external CSR data structure
  void computeLocalCSR(int num_ext_vars, const int *ext_vars, const int *rowp,
                       const int *cols, int lower, int upper, int nz_per_row,
                       int **_Arowp, int **_Acols, int *_Np, int **_Browp,
                       int **_Bcols);

  // Variables for the column map
  // ----------------------------
  int col_map_size;
  const int *col_map_vars;
  TACSNodeMap *row_map;

  // Information about assembling the non-zero pattern
  // -------------------------------------------------
  MPI_Comm comm;

  // Block size for this matrix
  int bsize;

  // Data destined for other processes
  // ---------------------------------
  int num_ext_procs;          // Number of processors that will be sent data
  int num_ext_rows;           // Total number of rows that are send data
  int *ext_procs;             // External proc numbers
  int *ext_count;             // Number of rows sent to each proc
  int *ext_row_ptr;           // Pointer from proc into the ext_rows array
  int *ext_rows;              // Row indices
  int *ext_rowp;              // Pointer into the rows
  int *ext_cols;              // Global column indices
  TacsScalar *ext_A;          // Pointer to the data accumulated on this proc
  MPI_Request *ext_requests;  // Requests for sending info

  // Data received from other processes
  // ----------------------------------
  int num_in_procs;  // Number of processors that give contributions
  int num_in_rows;   // Total number of rows given by other processors
  int *in_procs;     // Processor numbers that give info (num_in_procs)
  int *in_count;     // Count of rows from other processors (num_in_procs)
  int *in_row_ptr;   // Offset from proc index to location in rows
  int *in_rows;      // Row numbers for each row (num_in_rows)
  int *in_rowp;      // Pointer into the column numbers (num_in_rows)
  int *in_cols;      // Global column indices
  MPI_Request *in_requests;  // Requests for recving data
  TacsScalar *in_A;
};

#endif  // TACS_MAT_DISTRIBUTE_H
