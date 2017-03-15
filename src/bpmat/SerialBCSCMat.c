#include "SerialBCSCMat.h"


/*
  Compare integers for binary searches
*/
int compare_integers( const void *a, const void *b ){
  return (*(int*)a - *(int*)b);
}

/*
  Create a serial CSC matrix compatible with TACS

  This class serves as a simple wrapper for the BCSCMat class.  The
  main difficulty in using this class is that it uses a
  column-oriented storage format which makes applying boundary
  conditions difficult (since this involves zeroing a row of the
  matrix.) To circumvent this problem, we store the non-zero CSR data
  for the matrix as well, enabling quick access to the columns for
  each row. Note that the row indices in the CSC matrix are guaranteed
  to be sorted.

  input:
  rmap:             the variable->processor map for TACS unknowns
  bsize:            the block size (varsPerNode)
  num_block_rows:   number of block rows
  num_block_cols:   number of block columns
  block_rowp:       the block CSR non-zero data
  block_cols:       the columns in each row
  bcs:              the boundary conditions for the problem
*/
SerialBCSCMat::SerialBCSCMat( TACSVarMap *_rmap, int _bsize, 
                              int num_block_rows, 
                              int num_block_cols,
                              const int *block_rowp,
                              const int *block_cols ){
  rmap = _rmap;
  rmap->incref();
  mat = new BCSCMat(rmap->getMPIComm(), _bsize,
                    num_block_rows, num_block_cols,
                    block_rowp, block_cols);
  mat->incref();

  // Keep a copy of the CSR data to make applying boundary conditions
  // more efficient
  bsize = _bsize;
  nrows = num_block_rows;
  rowp = new int[ nrows+1 ];
  memcpy(rowp, block_rowp, (nrows+1)*sizeof(int));
  
  cols = new int[ rowp[nrows] ];
  memcpy(cols, block_cols, rowp[nrows]*sizeof(int));
}

/*
  Free the SerialBCSC matrix data
*/
SerialBCSCMat::~SerialBCSCMat(){
  // Free the CSR matrix data
  delete [] rowp;
  delete [] cols;

  // Decrease ref counts to the underlying data
  rmap->decref();
  mat->decref();
}

/*
  Zero all the entries in the matrix
*/
void SerialBCSCMat::zeroEntries(){
  mat->zeroEntries();
}

/*
  Add values to the matrix

  input:
  nrow:    the number of rows
  row:     the rows
  ncol:    the number of cols
  cols:    the columns
  nv:      the row-dimension of the dense values matrix
  mv:      the column dimension of the dense vaules matrix
  values:  the values to add to the matrix
*/
void SerialBCSCMat::addValues( int nrow, const int *row, 
                               int ncol, const int *col,
                               int nv, int mv, const TacsScalar *values ){
  mat->addMatBlockValues(nrow, row, ncol, col, values, nv);
}

/*
  Apply the Dirichlet boundary conditions to the matrix

  This code zeros the rows corresponding to the Dirichlet boundary
  conditions in the finite-element problem. This code uses the CSR
  matrix data to simplify searching for the columns with the rows
  corresponding to each boundary condition.
*/
void SerialBCSCMat::applyBCs( TACSBcMap *bcmap ){
  // Set up data so that we can quickly zero rows associated with the
  // boundary conditions in the column-oriented storage format. This
  // relies on the boundary conditions remaining fixed after they are
  // set/initialized. This is guaranteed by the TACSAssembler object.
  int ncols;
  const int *bptr, *aptr, *colp, *rows;
  TacsScalar *A;
  mat->getArrays(NULL, NULL, &ncols,
                 &bptr, &aptr, &colp, &rows, &A);

  const int *nodes, *vars;
  int nbcs = bcmap->getBCs(&nodes, &vars, NULL);
  for ( int i = 0; i < nbcs; i++ ){
    int node = nodes[i];

    // Search through the columns 
    for ( int jp = rowp[node]; jp < rowp[node+1]; jp++ ){
      int col = cols[jp];

      // Jump to column jp and search for the variables
      int col_size = colp[col+1] - colp[col];
      const int *col_rows = &rows[colp[col]];

      // Loop over each of the variables
      for ( int j = 0; j < bsize; j++ ){
        if (vars[i] & (1 << j)){
          int row = bsize*node + j;

          // Search for the row
          int *item = (int*)bsearch(&row, col_rows, col_size,
                                    sizeof(int), compare_integers);
          if (item){
            TacsScalar *a = &A[aptr[col] + bsize*(item - col_rows)];

            // Zero the row
            memset(a, 0, bsize*sizeof(TacsScalar));

            // Set the diagonal to the identity matrix
            if (node == col){
              a[j] = 1.0;
            }
          }
        }
      }
    }
  }
}

/*
  Create the TACSBVec object that matches with this matrix
*/
TACSVec *SerialBCSCMat::createVec(){
  return new TACSBVec(rmap, mat->getMaxBlockSize());
}

/*
  Perform a matrix multiplication
*/
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

/*
  Retrieve the underlying BCSCMat object
*/
BCSCMat *SerialBCSCMat::getBCSCMat(){
  return mat;
}

/*
  Create the preconditioner (a direct solve) associated with the
  SerialBCSCMat class

  This code uses a direct factorization and applyFactor can be used to
  solve the problem.
*/
SerialBCSCPc::SerialBCSCPc( SerialBCSCMat *mat ){
  fill = 10.0;
  pivot = new BCSCMatPivot(mat->getBCSCMat());
  pivot->incref();
}

/*
  Free the data associated with the preconditioner
*/
SerialBCSCPc::~SerialBCSCPc(){
  pivot->decref();
}

/*
  Factor the matrix using the BCSCMatPivot code and update the
  estimate for the fill-in experienced during the factorization.
*/
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

/*
  Apply the factorization to the input vector and store the result in
  the output vector without over-writing the right-hand-side.
*/
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
