#include <stdio.h>
#include "PMat.h"
#include "FElibrary.h"

/*
  Parallel matrix implementation

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

/*!
  The set up for the parallel block-CSR matrix.

  The parallel matrix is split into two parts that are identified
  in the initialization. The diagonal matrix and the off-diagonal 
  matrix. The off-diagonal matrix corresponds to the coupling 
  terms to the external-interface unknowns. The internal-interface
  unknowns must be ordered last on each process. External-interface
  unknowns can only be coupled to other interface-unknowns (either
  external or internal). Thus the global matrix can be represented as
  follows

  A_i = [ B_i, F_i ; G_i, C_i ]
  u_i = [ x_i, y_i ]^{T}

  On each process the unknowns are divided into internal-variables x_i, 
  and internal-iterface variables y_i.
  
  Each domain is coupled to other domains only through the interface
  variables y_i.

  A_i u_i + P * E_{ij} y_j = b_i

  where P = [ 0, I_{size(y_i)} ]^{T}

  The matrix structure outlined above can be exploited to achieve 
  efficient and effective parallel preconditioning.
*/
TACSPMat::TACSPMat( TACSVarMap *_rmap,
                    BCSRMat *_Aloc, BCSRMat *_Bext,
                    TACSBVecDistribute *_ext_dist ){
  init(_rmap, _Aloc, _Bext, _ext_dist);
}

TACSPMat::TACSPMat(){
  rmap = NULL;
  Aloc = NULL;
  Bext = NULL;
  ext_dist = NULL;
  ctx = NULL;
  x_ext = NULL;
  N = 0; Nc = 0; Np = 0; bsize = 0;
}

/*
  Initialize the PMat object
*/
void TACSPMat::init( TACSVarMap *_rmap,
                     BCSRMat *_Aloc, BCSRMat *_Bext,
                     TACSBVecDistribute *_ext_dist ){
  // Set the variable map and the local matrix components
  rmap = _rmap;
  Aloc = _Aloc;
  Bext = _Bext;
  rmap->incref();
  Aloc->incref();
  Bext->incref();
  
  // No external column map
  ext_dist = NULL;
  x_ext = NULL;

  N = Aloc->getRowDim();
  if (N != Aloc->getColDim()){
    fprintf(stderr, "PMat error: Block-diagonal matrix must be square\n");
    return;
  }

  Nc = Bext->getRowDim();
  Np = N-Nc;
  if (Nc > N){    
    fprintf(stderr, "PMat error: Block-diagonal matrix must be square\n");
    return;
  }

  int rank;
  MPI_Comm_rank(rmap->getMPIComm(), &rank);
  printf("[%d] PMat diagnostics: N = %d, Nc = %d\n", rank, N, Nc);

  // Copy the distribution array vector
  ext_dist = _ext_dist;
  ext_dist->incref();
  if (Bext->getColDim() != ext_dist->getDim()){
    fprintf(stderr, "PMat error: Dimensions of external variables and \
external block matrix do not match\n");
    return;
  }

  bsize = Aloc->getBlockSize();
  if (Bext->getBlockSize() != bsize){
    fprintf(stderr, "Block sizes do not match\n");
    return;
  }

  // Create a context for distributing the non-local unknowns
  ctx = ext_dist->createCtx(bsize);
  ctx->incref();

  // Allocate the external array
  int len = bsize*ext_dist->getDim();
  x_ext = new TacsScalar[ len ];
  memset(x_ext, 0, len*sizeof(TacsScalar));
  ext_offset = bsize*Np;
}

TACSPMat::~TACSPMat(){
  if (rmap){ rmap->decref(); }
  if (Aloc){ Aloc->decref(); }
  if (Bext){ Bext->decref(); }
  if (ext_dist){ ext_dist->decref(); }
  if (ctx){ ctx->decref(); }
  if (x_ext){ delete [] x_ext; }
}

/*!
  Determine the local dimensions of the matrix - the diagonal part
*/
void TACSPMat::getSize( int *_nr, int *_nc ){
  *_nr = N*bsize;
  *_nc = N*bsize;
}

/*!
  Zero all matrix-entries
*/
void TACSPMat::zeroEntries(){
  Aloc->zeroEntries();
  Bext->zeroEntries();
}

/*!
  Copy the values from the another matrix
*/
void TACSPMat::copyValues( TACSMat *mat ){
  TACSPMat *pmat = dynamic_cast<TACSPMat*>(mat);
  if (pmat){  
    Aloc->copyValues(pmat->Aloc);
    Bext->copyValues(pmat->Bext);
  }
  else {
    fprintf(stderr, "Cannot copy matrices of different types\n");
  }
}

/*!
  Scale the entries in the other matrices by a given scalar
*/
void TACSPMat::scale( TacsScalar alpha ){
  Aloc->scale(alpha);
  Bext->scale(alpha);
}

/*!
  Compute y <- y + alpha * x
*/
void TACSPMat::axpy( TacsScalar alpha, TACSMat *mat ){
  TACSPMat *pmat = dynamic_cast<TACSPMat*>(mat);
  if (pmat){
    Aloc->axpy(alpha, pmat->Aloc);
    Bext->axpy(alpha, pmat->Bext);
  }
  else {
    fprintf(stderr, "Cannot apply axpy to matrices of different types\n");
  }
}

/*!
  Compute y <- alpha * x + beta * y
*/
void TACSPMat::axpby( TacsScalar alpha, TacsScalar beta, TACSMat *mat ){
  TACSPMat *pmat = dynamic_cast<TACSPMat*>(mat);
  if (pmat){
    Aloc->axpby(alpha, beta, pmat->Aloc);
    Bext->axpby(alpha, beta, pmat->Bext);
  }
  else {
    fprintf(stderr, "Cannot apply axpby to matrices of different types\n");
  }
}

/*
  Add a scalar to the diagonal
*/
void TACSPMat::addDiag( TacsScalar alpha ){
  Aloc->addDiag(alpha);
}

/*!
  Matrix multiplication
*/
void TACSPMat::mult( TACSVec *txvec, TACSVec *tyvec ){
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);

    ext_dist->beginForward(ctx, x, x_ext);
    Aloc->mult(x, y);
    ext_dist->endForward(ctx, x, x_ext);
    Bext->multAdd(x_ext, &y[ext_offset], &y[ext_offset]);    
  }
  else {
    fprintf(stderr, "PMat type error: Input/output must be TACSBVec\n");
  }
}

/*!
  Access the underlying matrices
*/
void TACSPMat::getBCSRMat( BCSRMat ** A, BCSRMat ** B ){
  if (A){ *A = Aloc; }
  if (B){ *B = Bext; }
}

void TACSPMat::getRowMap( int *_bs, int *_N, int *_Nc ){
  if (_bs){ *_bs = bsize; }
  if (_Nc){ *_Nc = Nc; }
  if (_N){ *_N = N; }
}

void TACSPMat::getColMap( int *_bs, int *_M ){
  if (_bs){ *_bs = bsize; }
  if (_M){ *_M = N; }
}

void TACSPMat::getExtColMap( TACSBVecDistribute ** ext_map ){
  if (ext_map){ *ext_map = ext_dist; }
}

/*!
  Apply the boundary conditions

  This code applies the boundary conditions supplied to the matrix
*/
void TACSPMat::applyBCs( TACSBcMap *bcmap ){
  // Get the MPI rank and ownership range
  int mpi_rank;
  const int *ownerRange;
  MPI_Comm_rank(rmap->getMPIComm(), &mpi_rank);
  rmap->getOwnerRange(&ownerRange);

  // apply the boundary conditions
  const int *nodes, *vars;
  const TacsScalar *values;
  int nbcs = bcmap->getBCs(&nodes, &vars, &values);

  // Get the matrix values
  for ( int i = 0; i < nbcs; i++){
    // Find block i and zero out the variables associated with it
    if (nodes[i] >= ownerRange[mpi_rank] &&
        nodes[i] < ownerRange[mpi_rank+1]){
      int bvar  = nodes[i] - ownerRange[mpi_rank];
      int ident = 1; // Replace the diagonal with the identity matrix
      Aloc->zeroRow(bvar, vars[i], ident);

      // Now, check if the variable will be
      // in the off diagonal block (potentially)
      bvar = bvar - (N-Nc);
      if (bvar >= 0){
        ident = 0;
        Bext->zeroRow(bvar, vars[i], ident);
      }
    }
  }      
}

/*
  Create a vector for the matrix
*/
TACSVec *TACSPMat::createVec(){
  return new TACSBVec(rmap, Aloc->getBlockSize());
}

/*!
  Print the matrix non-zero pattern to the screen.
*/
void TACSPMat::printNzPattern( const char *fileName ){
  int mpi_rank;
  const int *ownerRange;
  MPI_Comm_rank(rmap->getMPIComm(), &mpi_rank);
  rmap->getOwnerRange(&ownerRange);

  // Get the sizes of the Aloc and Bext matrices
  int b, Na, Ma;
  const int *rowp, *cols;
  TacsScalar *Avals;
  Aloc->getArrays(&b, &Na, &Ma, &rowp, &cols, &Avals);
  
  int Nb, Mb;
  const int *browp, *bcols;
  TacsScalar *Bvals;
  Bext->getArrays(&b, &Nb, &Mb, &browp, &bcols, &Bvals);

  // Get the map between the global-external 
  // variables and the local variables (for Bext)
  TACSBVecIndices *bindex = ext_dist->getIndices();
  const int *col_vars;
  bindex->getIndices(&col_vars);

  FILE *fp = fopen(fileName, "w");
  if (fp){
    fprintf(fp, "VARIABLES = \"i\", \"j\" \nZONE T = \"Diagonal block %d\"\n", 
            mpi_rank);

    // Print out the diagonal components
    for ( int i = 0; i < Na; i++ ){
      for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
        fprintf(fp, "%d %d\n", i + ownerRange[mpi_rank], 
                cols[j] + ownerRange[mpi_rank]);
      }
    }
    
    if (browp[Nb] > 0){
      fprintf(fp, "ZONE T = \"Off-diagonal block %d\"\n", mpi_rank);
      // Print out the off-diagonal components
      for ( int i = 0; i < Nb; i++ ){
        for ( int j = browp[i]; j < browp[i+1]; j++ ){
          fprintf(fp, "%d %d\n", i + N-Nc + ownerRange[mpi_rank], 
                  col_vars[bcols[j]]);
        }
      }
    }
    
    fclose(fp);
  }
}

const char *TACSPMat::TACSObjectName(){
  return matName;
}

const char *TACSPMat::matName = "TACSPMat";

/*
  Build a simple SOR or Symmetric-SOR preconditioner for the matrix
*/
TACSGaussSeidel::TACSGaussSeidel( TACSPMat *_mat, int _zero_guess, 
                                  TacsScalar _omega, int _iters, 
                                  int _symmetric ){
  mat = _mat;
  mat->incref();

  // Get the on- and off-diagonal components of the matrix
  mat->getBCSRMat(&Aloc, &Bext);
  Aloc->incref();
  Bext->incref();

  // Create a vector to store temporary data for the relaxation
  TACSVec *tbvec = mat->createVec();
  bvec = dynamic_cast<TACSBVec*>(tbvec);
  if (bvec){
    bvec->incref();
  }
  else {
    fprintf(stderr, 
            "TACSGaussSeidel error: Input/output must be TACSBVec\n");
  }

  // Get the number of variables in the row map
  int bsize, N, Nc;
  mat->getRowMap(&bsize, &N, &Nc);

  // Compute the offset to the off-processor terms
  ext_offset = bsize*(N-Nc);

  // Get the external column map - a VecDistribute object
  mat->getExtColMap(&ext_dist);
  ext_dist->incref();
  ctx = ext_dist->createCtx(bsize);
  ctx->incref();

  // Compute the size of the external components
  int ysize = bsize*ext_dist->getDim();
  yext = new TacsScalar[ ysize ];  
  
  // Store the relaxation options 
  zero_guess = _zero_guess;
  omega = _omega;
  iters = _iters;
  symmetric = _symmetric;
}
  
/*
  Free the SOR preconditioner
*/
TACSGaussSeidel::~TACSGaussSeidel(){
  mat->decref();
  Aloc->decref();
  Bext->decref();
  ext_dist->decref();
  ctx->decref();
  delete [] yext;
  if (bvec){ bvec->decref(); }
}

/*
  Factor the diagonal of the matrix
*/
void TACSGaussSeidel::factor(){
  Aloc->factorDiag();
}

/*!
  Apply the preconditioner to the input vector
  Apply SOR to the system A y = x

  The SOR is applied by computing
  Aloc * y = x - Bext * yext

  Then applying the matrix-smoothing for the system of equations,

  Aloc * y = b

  where b = x - Bext * yext
*/
void TACSGaussSeidel::applyFactor( TACSVec *txvec, TACSVec *tyvec ){
  // Covert to TACSBVec objects
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    // Retrieve the values from the input vectors
    TacsScalar *x, *y, *b;
    yvec->getArray(&y);
    xvec->getArray(&x);
    bvec->getArray(&b);
    
    /*
    // Start smoothing the
    top_ext->beginForward(top_ctx, x, xtop);
    
    // Start smoothing the first group of interior nodes
    Aloc->applySOR(NULL, start, end, incr, ext_offset,
                   omega, b, NULL, x);

    top_ext->endForward(top_ctx, x, xtop);
    // Insert the top nodes into their correct locations


    // Smooth the top nodes
    

    // Send the top node values back to the destination




    // Finish smoothing the first group of interior nodes    
    Aloc->applySOR(NULL, start, end, incr, ext_offset,
                   omega, b, NULL, x);
    */

    if (zero_guess){
      yvec->zeroEntries();      
      Aloc->applySOR(x, y, omega, iters);
    }
    else {
      // Begin sending the external-interface values
      ext_dist->beginForward(ctx, y, yext);
      
      // Zero entries in the local vector
      bvec->zeroEntries();
      
      // Finish sending the external-interface unknowns
      ext_dist->endForward(ctx, y, yext);
      
      // Compute b[ext_offset] = Bext*yext
      Bext->mult(yext, &b[ext_offset]); 

      // Compute b = xvec - Bext*yext
      bvec->axpby(1.0, -1.0, xvec);         
      
      Aloc->applySOR(b, y, omega, iters);
    }  

    /*
    // Get the number of variables in the row map
    int bsize, N, Nc;
    mat->getRowMap(&bsize, &N, &Nc);

    int start = 0;
    int end = N;
    int incr = 1;
    Aloc->applySOR(NULL, start, end, incr, 0,
                   omega, x, NULL, y);
    
    start = N-1;
    end = -1;
    incr = -1;
    Aloc->applySOR(NULL, start, end, incr, 0,
                   omega, x, NULL, y);
    */
  }
  else {
    fprintf(stderr, 
            "TACSGaussSeidel type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrix
*/
void TACSGaussSeidel::getMat( TACSMat **_mat ){
  *_mat = mat;
}

/*!
  Build the additive Schwarz preconditioner 
*/
TACSAdditiveSchwarz::TACSAdditiveSchwarz( TACSPMat *_mat, 
                                          int levFill, double fill ){
  mat = _mat;
  mat->incref();

  // Get the underlying matrices in the distributed matrix class
  BCSRMat *B;
  mat->getBCSRMat(&Aloc, &B);
  Aloc->incref();

  // Form the preconditioner matrix for the on-processor (block-diagonal)
  // components of the matrix. Incref the pointer to the matrix
  Apc = new BCSRMat(mat->getMPIComm(), Aloc, levFill, fill);
  Apc->incref();

  alpha = 0.0; // Diagonal scalar to be added to the preconditioner
}

/*
  Free the memory from the additive Schwarz preconditioner
*/
TACSAdditiveSchwarz::~TACSAdditiveSchwarz(){
  mat->decref();
  Aloc->decref();
  Apc->decref();
}

/*
  Add the diagonal shift factor to the preconditioner. The shift
  defaults to zero
*/
void TACSAdditiveSchwarz::setDiagShift( TacsScalar _alpha ){
  alpha = _alpha;
}

/*
  Factor the preconditioner by copying the values from the 
  block-diagonal matrix and then factoring the copy.
*/
void TACSAdditiveSchwarz::factor(){
  Apc->copyValues(Aloc);
  if (alpha != 0.0){ 
    Apc->addDiag(alpha);
  }
  Apc->factor();
}

/*!
  Apply the preconditioner to the input vector

  For the additive Schwarz method that simply involves apply the ILU 
  factorization of the diagonal to the input vector:

  y = U^{-1} L^{-1} x
*/
void TACSAdditiveSchwarz::applyFactor( TACSVec *txvec, TACSVec *tyvec ){
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){  
    // Apply the ILU factorization to a vector
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);
    
    Apc->applyFactor(x, y);
  }
  else {
    fprintf(stderr, 
            "TACSAdditiveSchwarz type error: Input/output must be TACSBVec\n");
  }
}

/*!
  Apply the preconditioner to the input vector

  For the additive Schwarz method that simply involves apply the ILU 
  factorization of the diagonal to the input vector:

  y = U^{-1} L^{-1} y
*/
void TACSAdditiveSchwarz::applyFactor( TACSVec *txvec ){
  TACSBVec *xvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);

  if (xvec){  
    // Apply the ILU factorization to a vector
    // This is the default Additive-Scharwz method
    TacsScalar *x;
    xvec->getArray(&x);
    
    Apc->applyFactor(x);
  }
  else {
    fprintf(stderr, 
            "TACSAdditiveSchwarz type error: Input/output must be TACSBVec\n");
  }
}

/*!
  Retrieve the underlying matrix
*/
void TACSAdditiveSchwarz::getMat( TACSMat **_mat ){
  *_mat = mat;
}

/*!
  The approximate Schur preconditioner class.
*/
TACSApproximateSchur::TACSApproximateSchur( TACSPMat *_mat, 
                                            int levFill, double fill, 
                                            int inner_gmres_iters, 
                                            double inner_rtol, 
                                            double inner_atol ){
  mat = _mat;
  mat->incref();

  // Copy the diagonal matrix
  BCSRMat *Bext;
  mat->getBCSRMat(&Aloc, &Bext);
  Aloc->incref();

  int rank, size;
  MPI_Comm_rank(mat->getMPIComm(), &rank);
  MPI_Comm_size(mat->getMPIComm(), &size);

  Apc = new BCSRMat(mat->getMPIComm(), Aloc, levFill, fill);
  Apc->incref();
  alpha = 0.0;

  // Check if we're dealing with a serial case here...
  gsmat = NULL;
  rvec = wvec = NULL;
  inner_ksm = NULL;
  
  if (size > 1){
    gsmat = new TACSGlobalSchurMat(mat, Apc);
    gsmat->incref();

    TACSVec *trvec = gsmat->createVec();
    TACSVec *twvec = gsmat->createVec();
    rvec = dynamic_cast<TACSBVec*>(trvec);
    wvec = dynamic_cast<TACSBVec*>(twvec);

    // The code relies on these vectors being TACSBVecs
    if (rvec && twvec){
      rvec->incref();
      wvec->incref();
    }
    else {
      fprintf(stderr, 
              "TACSApproximateSchur error: Input/output must be TACSBVec\n");
    }
    
    int nrestart = 0;
    inner_ksm = new GMRES(gsmat, inner_gmres_iters, nrestart);
    inner_ksm->incref();
    inner_ksm->setTolerances(inner_rtol, inner_atol);

    int b, N, Nc;    
    mat->getRowMap(&b, &N, &Nc);
    
    var_offset = N-Nc;
    start = b*(N-Nc);
    end   = b*N; 
  }
}

/*
  Free the data associated with the approximate Schur preconditioner
*/
TACSApproximateSchur::~TACSApproximateSchur(){
  Aloc->decref();
  Apc->decref();
  mat->decref();
  if (gsmat){ gsmat->decref(); }
  if (rvec){ rvec->decref(); }
  if (wvec){ wvec->decref(); }
  if (inner_ksm){ inner_ksm->decref(); }
}

void TACSApproximateSchur::setDiagShift( TacsScalar _alpha ){
  alpha = _alpha;
}

void TACSApproximateSchur::setMonitor( KSMPrint *ksm_print ){
  if (inner_ksm){
    inner_ksm->setMonitor(ksm_print);
  }
}       

/*
  Factor preconditioner based on the values in the matrix.
*/
void TACSApproximateSchur::factor(){
  Apc->copyValues(Aloc);
  if (alpha != 0.0){
    Apc->addDiag(alpha);
  }
  Apc->factor();
}

/*
  Print the non-zero pattern to a file
*/
void TACSApproximateSchur::printNzPattern( const char *fileName ){
  // Get the sizes of the Aloc and Bext matrices
  int b, Na, Ma;
  const int *rowp;
  const int *cols;
  TacsScalar *Avals;
  Apc->getArrays(&b, &Na, &Ma,
                 &rowp, &cols, &Avals);

  BCSRMat *Aloc, *Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  int Nb, Mb;
  const int *browp;
  const int *bcols;
  TacsScalar *Bvals;
  Bext->getArrays(&b, &Nb, &Mb,
                  &browp, &bcols, &Bvals);

  // Get the map between the global-external 
  // variables and the local variables (for Bext)
  TACSBVecDistribute *ext_dist;
  mat->getExtColMap(&ext_dist);
  TACSBVecIndices *bindex = ext_dist->getIndices();
  const int *col_vars;
  bindex->getIndices(&col_vars);

  TACSVarMap *rmap = mat->getRowMap();
  int mpi_rank;
  MPI_Comm_rank(rmap->getMPIComm(), &mpi_rank);
  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  FILE *fp = fopen(fileName, "w");
  fprintf(fp, "VARIABLES = \"i\", \"j\" \nZONE T = \"Diagonal block %d\"\n", 
          mpi_rank);

  // Print out the diagonal components
  for ( int i = 0; i < Na; i++ ){
    for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
      fprintf(fp, "%d %d\n", i + ownerRange[mpi_rank], 
              cols[j] + ownerRange[mpi_rank]);
    }
  }

  if (browp[Nb] > 0){
    fprintf(fp, "ZONE T = \"Off-diagonal block %d\"\n", mpi_rank);
    // Print out the off-diagonal components
    for ( int i = 0; i < Nb; i++ ){
      for ( int j = browp[i]; j < browp[i+1]; j++ ){
        fprintf(fp, "%d %d\n", i + var_offset + ownerRange[mpi_rank], 
                col_vars[bcols[j]]);
      }
    }
  }

  fclose(fp);
}

/*
  For the approximate Schur method the following calculations are
  caried out:
   
  Compute y = U^{-1} L^{-1} x
  Determine the interface unknowns v_{0} = P * y (the last Nc equations)
  Normalize v_{0} = v_{0}/|| v_{0} ||
  for j = 1, m
  .   Exchange interface unknowns. Obtain v_{j,ext}.
  .   Apply approximate Schur complement. 
  .   v_{j+1} = S^{-1} Bext * v_{j,ext} + v_{j}
  .   For k = 1, j
  .      h_{k,j} = dot(v_{j+1}, v_{j})
  .      v_{j+1} = v_{j} - h_{k,j} v_{k]
  .   h_{j+1,j} = || v_{j+1} ||
  .   v_{j+1} = v_{j+1}/h_{j+1,j}
  Compute z_{m} = argmin || \beta e_{1} - H_{m} z ||_{2}
  V_{m} = [ v_{1}, v_{2}, .. v_{m} ]
  compute w = V_{m} z_{m}
  Exchange interface variables w
  Compute: x <- x - [0,w]
  Compute: y <- U^{-1} L^{-1} x

  -----------------------------------------

  Application of the Schur preconditioner:

  Perform the factorization:
  
  A_{i} = [ L_b          0   ][ U_b  L_b^{-1} E ]
  .       [ F U_b^{-1}   L_s ][ 0    U_s        ]

  Find an approximate solution to:

  [ B  E ][ x_i ]                     [ f_i ]
  [ F  C ][ y_i ] + [ sum F_j y_j ] = [ g_i ]

  Compute the modified RHS:

  g_i' = U_s^{-1} L_s^{-1} ( g_i - F B^{-1} f_i)
  .    = U_s^{-1} L_s^{-1} ( g_i - F U_b^{-1} L_b^{-1} f_i)

  Solve for the interface unknowns (with GMRES):

  y_i + sum  U_s^{-1} L_s^{-1} F_j y_j = g_i'

  Compute the interior unknowns:

  x_i = U_b^{-1} L_b^{-1} ( f_i - E * y_i)  
*/
void TACSApproximateSchur::applyFactor( TACSVec * txvec, TACSVec * tyvec ){
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    // Apply the ILU factorization to a vector
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);
    
    if (inner_ksm){
      // y = L_{B}^{-1} x
      // Compute the modified RHS g' = U_s^{-1} L_s^{-1} (g_i - F B^{-1} f_i)
      Apc->applyLower(x, y);
      Apc->applyPartialUpper(&y[start], var_offset);
      
      TacsScalar *r, *w;
      rvec->getArray(&r);
      wvec->getArray(&w);
      memcpy(r, &y[start], (end-start)*sizeof(TacsScalar));
      
      // Solve the global Schur system: S * w = r
      inner_ksm->solve(rvec, wvec);
      memcpy(&y[start], w, (end-start)*sizeof(TacsScalar));
            
      // Compute the interior unknowns from the interface values
      // x_i = U_b^{-1} L_b^{-1} (f_i - E y_i)  
      // x_i = U_b^{-1} (L_b^{-1} f_i - L_b^{-1} E y_i)  
      Apc->applyFactorSchur(y, var_offset);
    }
    else {
      // y = U^{-1} L^{-1} x
      Apc->applyFactor(x, y);
    }
  }
  else {
    fprintf(stderr, 
            "ApproximateSchur type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrix
*/
void TACSApproximateSchur::getMat( TACSMat **_mat ){
  *_mat = mat;
}

/*!
  The block-Jacobi-preconditioned approximate global Schur matrix. 

  This matrix is used within the ApproximateSchur preconditioner.
*/
TACSGlobalSchurMat::TACSGlobalSchurMat( TACSPMat *mat, BCSRMat *_Apc ){
  Apc = _Apc;
  Apc->incref();

  BCSRMat * Aloc;
  mat->getBCSRMat(&Aloc, &Bext);
  Bext->incref();

  int bsize, N, Nc;
  mat->getRowMap(&bsize, &N, &Nc);

  varoffset = N-Nc;
  nvars = bsize*Nc;
  rmap = new TACSVarMap(mat->getMPIComm(), Nc);

  mat->getExtColMap(&ext_dist);
  ext_dist->incref();
  ctx = ext_dist->createCtx(bsize);
  ctx->incref();

  int xsize = bsize*ext_dist->getDim();
  x_ext = new TacsScalar[ xsize ];  
}

TACSGlobalSchurMat::~TACSGlobalSchurMat(){
  Apc->decref();
  Bext->decref(); 
  ext_dist->decref();
  ctx->decref();
  delete [] x_ext;
}

void TACSGlobalSchurMat::getSize( int *_nr, int *_nc ){
  // Get the local dimensions of the matrix
  *_nr = nvars;
  *_nc = nvars;
}

/*
  Compute y <- A * x 
*/
void TACSGlobalSchurMat::mult( TACSVec *txvec, TACSVec *tyvec ){
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);  
    
    // Begin sending the external-interface values
    ext_dist->beginForward(ctx, x, x_ext, varoffset);
    
    // Finish sending the external-interface unknowns
    ext_dist->endForward(ctx, x, x_ext, varoffset);
    Bext->mult(x_ext, y);
    
    // Apply L^{-1}
    Apc->applyPartialLower(y, varoffset);
    
    // Apply U^{-1}
    Apc->applyPartialUpper(y, varoffset);
    
    // Finish the matrix-vector product
    yvec->axpy(1.0, xvec);
  }
  else {
    fprintf(stderr, 
            "GlobalSchurMat type error: Input/output must be TACSBVec\n");
  }
}

/*
  Compute y <- Bext * xext
*/
void TACSGlobalSchurMat::multOffDiag( TACSBVec *xvec, TACSBVec *yvec ){  
  TacsScalar *x, *y;
  xvec->getArray(&x);
  yvec->getArray(&y);  

  // Begin sending the external-interface values
  ext_dist->beginForward(ctx, x, x_ext, varoffset);
  
  // Finish sending the external-interface unknowns
  ext_dist->endForward(ctx, x, x_ext, varoffset);
  Bext->mult(x_ext, y);
}

// Return a new Vector
TACSVec *TACSGlobalSchurMat::createVec(){
  return new TACSBVec(rmap, Apc->getBlockSize());
}
