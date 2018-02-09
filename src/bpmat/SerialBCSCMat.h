/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#ifndef TACS_SERIAL_BCSC_MATRIX_H
#define TACS_SERIAL_BCSC_MATRIX_H

/*
  The following class wraps a serial BCSC matrix
*/
#include "BCSCMatPivot.h"
#include "BVec.h"
#include "KSM.h"

class SerialBCSCMat : public TACSMat {
 public:
  SerialBCSCMat( TACSVarMap *_rmap, int bsize, 
                 int num_block_rows, int num_block_cols,
                 const int *block_rowp, const int *block_cols );
  ~SerialBCSCMat();

  // Set entries into the matrix
  void zeroEntries();
  void addValues( int nrow, const int *row, 
                  int ncol, const int *col,
                  int nv, int mv, const TacsScalar *values );
  void addWeightValues( int nvars, const int *varp, const int *vars,
                        const TacsScalar *weights,
                        int nv, int mv, const TacsScalar *values,
                        MatrixOrientation matOr=NORMAL );
  void applyBCs( TACSBcMap *bcmap );

  // Create vectors
  // --------------
  TACSVec *createVec();

  // Operations required for solving problems
  // ----------------------------------------
  void mult( TACSVec *tx, TACSVec *ty );

  // Retrieve the pointer to the underlying matrix
  // ---------------------------------------------
  BCSCMat *getBCSCMat();

 private:
  // The non-zero information associated with the CSR data
  int bsize, nrows;
  int *rowp, *cols;

  // The variable
  TACSVarMap *rmap;

  // The serial CSC matrix itself
  BCSCMat *mat;
};

class SerialBCSCPc : public TACSPc {
 public:
  SerialBCSCPc( SerialBCSCMat *_mat );
  ~SerialBCSCPc();

  // Factor the matrix and apply the factorization
  // ---------------------------------------------
  void factor();
  void applyFactor( TACSVec *txvec, TACSVec *tyvec );
  void getMat( TACSMat **_mat );

 private:
  SerialBCSCMat *mat;
  double fill;
  BCSCMatPivot *pivot;
};

#endif // TACS_SERIAL_BCSC_MATRIX_H
