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
                 const int *block_rowp, const int *block_cols,
                 TACSBcMap *_bcs );
  ~SerialBCSCMat();

  // Set entries into the matrix
  void zeroEntries();
  void addValues( int nrow, const int *row, 
                  int ncol, const int *col,
                  int nv, int mv, const TacsScalar *values );
  void applyBCs();
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
  TACSBcMap *bcs;

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

 private:
  double fill;
  BCSCMatPivot *pivot;
};

#endif // TACS_SERIAL_BCSC_MATRIX_H
