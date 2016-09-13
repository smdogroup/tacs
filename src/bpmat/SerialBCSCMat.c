#include "SerialBCSCMat.h"


SerialBCSCMat::SerialBCSCMat( TACSVarMap *_rmap, int bsize, 
                              int num_block_rows, 
                              int num_block_cols,
                              const int *block_rowp,
                              const int *block_cols,
                              TACSBcMap *_bcs ){
  rmap = _rmap;
  rmap->incref();
  mat = new BCSCMat(rmap->getMPIComm(), bsize,
                    num_block_rows, num_block_cols,
                    block_rowp, block_cols);
  mat->incref();
  bcs = _bcs;
  if (bcs){ bcs->incref(); } 
}

SerialBCSCMat::~SerialBCSCMat(){
  rmap->decref();
  mat->decref();
  bcs->decref();
}

// Set entries into the matrix
void SerialBCSCMat::zeroEntries(){
  mat->zeroEntries();
}
void SerialBCSCMat::addValues( int nrow, const int *row, 
                               int ncol, const int *col,
                               int nv, int mv, const TacsScalar *values ){
  mat->addMatBlockValues(nrow, row, ncol, col, values, nv);
}

void SerialBCSCMat::applyBCs(){
    
}

// Create vectors/retrieve sizes
// -----------------------------
TACSVec *SerialBCSCMat::createVec(){
  return new TACSBVec(rmap, mat->getMaxBlockSize(), bcs);
}

// Operations required for solving problems
// ----------------------------------------
void SerialBCSCMat::mult( TACSVec *tx, TACSVec *ty ){
  // Dynamic cast to TACSBVec
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(tx);
  yvec = dynamic_cast<TACSBVec*>(ty);

  if (xvec && yvec){
    // Get the entries in the array
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);
      
    mat->mult(x, y, 1);
  }
}

BCSCMat *SerialBCSCMat::getBCSCMat(){
  return mat;
}


SerialBCSCPc::SerialBCSCPc( SerialBCSCMat *mat ){
  fill = 10.0;
  pivot = new BCSCMatPivot(mat->getBCSCMat());
  pivot->incref();
}

SerialBCSCPc::~SerialBCSCPc(){
  pivot->decref();
}

void SerialBCSCPc::factor(){
  // Factor the matrix
  double new_fill = pivot->factor(fill);
  
  // Update the matrix fill estimate
  if (new_fill > fill){
    fill = new_fill;
  }
  else {
    fill += 0.25*(new_fill - fill);
  }
}

void SerialBCSCPc::applyFactor( TACSVec *txvec, TACSVec *tyvec ){
  // Covert to TACSBVec objects
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    // Get the arrays
    TacsScalar *x, *y;
    int size = xvec->getArray(&x);
    yvec->getArray(&y);

    // Copy over the array
    memcpy(y, x, size*sizeof(TacsScalar));
      
    // Apply the factor
    pivot->applyFactor(y, 1);
  }
}
