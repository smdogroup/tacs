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
PMat::PMat( VarMap * _rmap,
	    BCSRMat * _Aloc, BCSRMat * _Bext,
	    BVecDistribute * _col_map,
	    BCMap * _bcs ){
  init(_rmap, _Aloc, _Bext, _col_map, _bcs);
}

PMat::PMat(){
  rmap = NULL;
  Aloc = NULL;
  Bext = NULL;
  bcs  = NULL;
  col_map = NULL;
  x_ext = NULL;
  N = 0; Nc = 0; Np = 0; bsize = 0;
}

void PMat::init( VarMap * _rmap,
		 BCSRMat * _Aloc, BCSRMat * _Bext,
		 BVecDistribute * _col_map,
		 BCMap * _bcs ){

  rmap = _rmap;  rmap->incref();  
  Aloc = _Aloc;  Aloc->incref();
  Bext = _Bext;  Bext->incref();
  bcs  = _bcs;   if (bcs){ bcs->incref(); }
  col_map = NULL;
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

  col_map = _col_map;
  col_map->incref();
  if (Bext->getColDim() != col_map->getDim()){
    fprintf(stderr, "PMat error: Dimensions of external variables and \
external block matrix do not match\n");
    return;
  }

  bsize = Aloc->getBlockSize();
  if (Bext->getBlockSize() != bsize ||
       col_map->getBlockSize() != bsize){
    fprintf(stderr, "Block sizes do not match\n");
    return;
  }

  int len = col_map->getSize();
  x_ext = new TacsScalar[ len ];
  memset(x_ext,0,len*sizeof(TacsScalar));
  ext_offset = bsize * Np;
}

PMat::~PMat(){
  if (rmap){ rmap->decref(); }
  if (Aloc){ Aloc->decref(); }
  if (Bext){ Bext->decref(); }
  if (bcs){ bcs->decref(); }
  if (col_map){ col_map->decref(); }
  if (x_ext){ delete [] x_ext; }
}

/*!
  Determine the local dimensions of the matrix - the diagonal part
*/
void PMat::getSize( int * _nr, int * _nc ){
  *_nr = N*bsize;
  *_nc = N*bsize;
}

/*!
  Zero all matrix-entries
*/
void PMat::zeroEntries(){
  Aloc->zeroEntries();
  Bext->zeroEntries();
}

/*!
  Copy the values from the another matrix
*/
void PMat::copyValues( TACSMat * mat ){
  PMat * pmat = dynamic_cast<PMat*>(mat);
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
void PMat::scale( TacsScalar alpha ){
  Aloc->scale(alpha);
  Bext->scale(alpha);
}

/*!
  Compute y <- y + alpha * x
*/
void PMat::axpy( TacsScalar alpha, TACSMat * mat ){
  PMat * pmat = dynamic_cast<PMat*>(mat);
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
void PMat::axpby( TacsScalar alpha, TacsScalar beta, TACSMat * mat ){
  PMat * pmat = dynamic_cast<PMat*>(mat);
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
void PMat::addDiag( TacsScalar alpha ){
  Aloc->addDiag(alpha);
}

// Functions required for solving linear systems
// ---------------------------------------------

/*!
  Matrix multiplication
*/
void PMat::mult( TACSVec * txvec, TACSVec * tyvec ){
  BVec *xvec, *yvec;

  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);

    col_map->beginForward(xvec, x_ext);
    Aloc->mult(x, y);
    col_map->endForward(xvec, x_ext);
    Bext->multAdd(x_ext, &y[ext_offset], &y[ext_offset]);    
  }
  else {
    fprintf(stderr, "PMat type error: Input/output must be BVec\n");
  }
}

/*!
  Access the underlying matrices
*/
void PMat::getBCSRMat( BCSRMat ** A, BCSRMat ** B ){
  *A = Aloc;
  *B = Bext;
}

void PMat::getRowMap( int * _bs, int * _N, int * _Nc ){
  *_bs = bsize;
  *_Nc = Nc;
  *_N = N;
}

void PMat::getColMap( int * _bs, int * _M ){
  *_bs = bsize;
  *_M = N;
}

void PMat::getExtColMap( BVecDistribute ** ext_map ){
  *ext_map = col_map;
}

/*!
  Apply the boundary conditions

  For the serial case, this simply involves zeroing the appropriate rows. 
*/
void PMat::applyBCs(){
  if (bcs){
    const int * ownerRange;
    int mpiRank, mpiSize;
    rmap->getOwnerRange(&ownerRange, &mpiRank, &mpiSize);

    // apply the boundary conditions
    const int *local, *global, *var_ptr, *vars;
    const TacsScalar *values;
    int nbcs = bcs->getBCs(&local, &global, &var_ptr, &vars, &values);

    // Get the matrix values
    for ( int i = 0; i < nbcs; i++){
      // Find block i and zero out the variables associated with it
      if (global[i] >= ownerRange[mpiRank] &&
	   global[i] <  ownerRange[mpiRank+1]){

	int bvar  = global[i] - ownerRange[mpiRank];
	int start = var_ptr[i];
	int nvars = var_ptr[i+1] - start;

	int ident = 1; // Replace the diagonal with the identity matrix
	Aloc->zeroRow(bvar, nvars, &vars[start], ident);

	// Now, check if the variable will be
	// in the off diagonal block (potentially)
	bvar = bvar - (N-Nc);
	if (bvar >= 0){
	  ident = 0;
	  Bext->zeroRow(bvar, nvars, &vars[start], ident);
	}
      }
    }      
  }
}

TACSVec * PMat::createVec(){
  return new BVec(rmap, bcs);
}

/*!
  Print the matrix non-zero pattern to the screen.
*/
void PMat::printNzPattern( const char * fileName ){
  const int * ownerRange;
  int mpiRank, mpiSize;
  rmap->getOwnerRange(&ownerRange, &mpiRank, &mpiSize);

  // Get the sizes of the Aloc and Bext matrices
  int b, Na, Ma;
  const int * rowp;
  const int * cols;
  TacsScalar * Avals;
  Aloc->getArrays(&b, &Na, &Ma,
                  &rowp, &cols, &Avals);
  
  int Nb, Mb;
  const int * browp;
  const int * bcols;
  TacsScalar * Bvals;
  Bext->getArrays(&b, &Nb, &Mb,
                  &browp, &bcols, &Bvals);

  // Get the map between the global-external 
  // variables and the local variables (for Bext)
  int * col_vars;
  BVecIndices * bindex = col_map->getBVecIndices();
  bindex->getIndices(&col_vars);

  FILE * fp = fopen(fileName, "w");
  if (fp){
    fprintf(fp, "VARIABLES = \"i\", \"j\" \nZONE T = \"Diagonal block %d\"\n", 
            mpiRank);

    // Print out the diagonal components
    for ( int i = 0; i < Na; i++ ){
      for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
        fprintf(fp, "%d %d\n", i + ownerRange[mpiRank], 
                cols[j] + ownerRange[mpiRank]);
      }
    }
    
    if (browp[Nb] > 0){
      fprintf(fp, "ZONE T = \"Off-diagonal block %d\"\n", mpiRank);
      // Print out the off-diagonal components
      for ( int i = 0; i < Nb; i++ ){
        for ( int j = browp[i]; j < browp[i+1]; j++ ){
          fprintf(fp, "%d %d\n", i + N-Nc + ownerRange[mpiRank], 
                  col_vars[bcols[j]]);
        }
      }
    }
    
    fclose(fp);
  }
}


const char * PMat::TACSObjectName(){
  return matName;
}

const char * PMat::matName = "PMat";

/*!
  Build a simple SOR or Symmetric-SOR preconditioner for A
*/
PSOR::PSOR( PMat * mat, int _zero_guess, TacsScalar _omega, int _iters, 
	    int _isSymmetric ){
  mat->getBCSRMat(&Aloc, &Bext);
  Aloc->incref();
  Bext->incref();

  TACSVec * tbvec = mat->createVec();
  bvec = dynamic_cast<BVec*>(tbvec);
  if (bvec){
    bvec->incref();
  }
  else {
    fprintf(stderr, "PSOR type error: Input/output must be BVec\n");
  }

  // Get the external column map - a VecDistribute object
  mat->getExtColMap(&col_map);
  col_map->incref();

  // Get the number of variables in the row map
  int b, N, Nc;
  mat->getRowMap(&b, &N, &Nc);

  ext_offset = b*(N-Nc);
  
  int ysize = col_map->getSize();
  yext = new TacsScalar[ ysize ];  

  zero_guess = _zero_guess;
  omega = _omega;
  iters = _iters;
  isSymmetric = _isSymmetric;
}

PSOR::~PSOR(){
  Aloc->decref();
  Bext->decref();
  if (bvec){ bvec->decref(); }
  col_map->decref();

  delete [] yext;
}

/*
  Factor the diagonal of the matrix
*/
void PSOR::factor(){
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
void PSOR::applyFactor( TACSVec * txvec, TACSVec * tyvec ){
  BVec *xvec, *yvec;

  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

  if (xvec && yvec){
    // Apply the ILU factorization to a vector
    TacsScalar *x, *y, *b;

    yvec->getArray(&y);
    xvec->getArray(&x);
    bvec->getArray(&b);
    
    if (zero_guess){
      yvec->zeroEntries();
      
      if (isSymmetric){
        Aloc->applySSOR(x, y, omega, iters);
      }
      else {
        Aloc->applySOR(x, y, omega, iters);
      }
    }
    else {
      // Begin sending the external-interface values
      col_map->beginForward(yvec, yext);
      
      bvec->zeroEntries();
      
      // Finish sending the external-interface unknowns
      col_map->endForward(yvec, yext);
      
      Bext->mult(yext, &b[ext_offset]);     // b = Bext * yext
      bvec->axpby(1.0, -1.0, xvec);         // b = xvec - Bext * yext
      
      if (isSymmetric){
        Aloc->applySSOR(b, y, omega, iters);
      }
      else {
        Aloc->applySOR(b, y, omega, iters);
      }
    }
  }
  else {
    fprintf(stderr, "PSOR type error: Input/output must be BVec\n");
  }
}

/*!
  Build the additive Schwarz preconditioner 
*/
AdditiveSchwarz::AdditiveSchwarz( PMat * mat, int levFill, double fill ){
  BCSRMat * B;

  mat->getBCSRMat(&Aloc, &B);
  Aloc->incref();

  Apc = new BCSRMat(mat->getMPIComm(), Aloc, levFill, fill);
  Apc->incref();

  alpha = 0.0; // Diagonal scalar to be added to the preconditioner
}

AdditiveSchwarz::~AdditiveSchwarz(){
  Aloc->decref();
  Apc->decref();
}

void AdditiveSchwarz::setDiagShift( TacsScalar _alpha ){
  alpha = _alpha;
}

void AdditiveSchwarz::factor(){
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
void AdditiveSchwarz::applyFactor( TACSVec * txvec, TACSVec * tyvec ){
  BVec *xvec, *yvec;
  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

  if (xvec && yvec){  
    // Apply the ILU factorization to a vector
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);
    
    Apc->applyFactor(x, y);
  }
  else {
    fprintf(stderr, "AdditiveSchwarz type error: Input/output must be BVec\n");
  }
}

/*!
  Apply the preconditioner to the input vector

  For the additive Schwarz method that simply involves apply the ILU 
  factorization of the diagonal to the input vector:

  y = U^{-1} L^{-1} y
*/
void AdditiveSchwarz::applyFactor( TACSVec * txvec ){
  BVec *xvec;
  xvec = dynamic_cast<BVec*>(txvec);

  if (xvec){  
    // Apply the ILU factorization to a vector
    // This is the default Additive-Scharwz method
    TacsScalar * x;
    xvec->getArray(&x);
    
    Apc->applyFactor(x);
  }
  else {
    fprintf(stderr, "AdditiveSchwarz type error: Input/output must be BVec\n");
  }
}

/*!
  The approximate Schur preconditioner class.
*/
ApproximateSchur::ApproximateSchur( PMat * _mat, int levFill, double fill, 
				    int inner_gmres_iters, double inner_rtol, 
				    double inner_atol ){
  mat = _mat;
  mat->incref();

  BCSRMat * Bext;
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
    gsmat = new GlobalSchurMat(mat, Apc);
    gsmat->incref();

    TACSVec * trvec = gsmat->createVec();
    TACSVec * twvec = gsmat->createVec();
    rvec = dynamic_cast<BVec*>(trvec);
    wvec = dynamic_cast<BVec*>(twvec);

    // The code relies on these vectors being BVecs
    if (rvec && twvec){
      rvec->incref();
      wvec->incref();
    }
    else {
      fprintf(stderr, "ApproximateSchur type error: Input/output must be BVec\n");
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

ApproximateSchur::~ApproximateSchur(){
  Aloc->decref();
  Apc->decref();
  mat->decref();

  if (gsmat){ gsmat->decref(); }
  if (rvec){ rvec->decref(); }
  if (wvec){ wvec->decref(); }
  if (inner_ksm){ inner_ksm->decref(); }
}

void ApproximateSchur::setDiagShift( TacsScalar _alpha ){
  alpha = _alpha;
}

void ApproximateSchur::setMonitor( KSMPrint * ksm_print ){
  if (inner_ksm){
    inner_ksm->setMonitor(ksm_print);
  }
} 	

/*
  Factor preconditioner based on the values in the matrix.
*/
void ApproximateSchur::factor(){
  Apc->copyValues(Aloc);
  if (alpha != 0.0){
    Apc->addDiag(alpha);
  }
  Apc->factor();
}

void ApproximateSchur::printNzPattern( const char * fileName ){
 // Get the sizes of the Aloc and Bext matrices
  int b, Na, Ma;
  const int * rowp;
  const int * cols;
  TacsScalar * Avals;
  Apc->getArrays(&b, &Na, &Ma,
                 &rowp, &cols, &Avals);

  BCSRMat * Aloc, * Bext;
  mat->getBCSRMat(&Aloc, &Bext);

  int Nb, Mb;
  const int * browp;
  const int * bcols;
  TacsScalar * Bvals;
  Bext->getArrays(&b, &Nb, &Mb,
                  &browp, &bcols, &Bvals);

  // Get the map between the global-external 
  // variables and the local variables (for Bext)
  int * col_vars;
  BVecDistribute * col_map;
  mat->getExtColMap(&col_map);
  BVecIndices * bindex = col_map->getBVecIndices();
  bindex->getIndices(&col_vars);

  VarMap * rmap = mat->getRowMap();
  int mpiRank, s;
  const int * ownerRange;
  rmap->getOwnerRange(&ownerRange, &mpiRank, &s);

  FILE * fp = fopen(fileName, "w");
  fprintf(fp, "VARIABLES = \"i\", \"j\" \nZONE T = \"Diagonal block %d\"\n", 
          mpiRank);

  // Print out the diagonal components
  for ( int i = 0; i < Na; i++ ){
    for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
      fprintf(fp, "%d %d\n", i + ownerRange[mpiRank], 
              cols[j] + ownerRange[mpiRank]);
    }
  }

  if (browp[Nb] > 0){
    fprintf(fp, "ZONE T = \"Off-diagonal block %d\"\n", mpiRank);
    // Print out the off-diagonal components
    for ( int i = 0; i < Nb; i++ ){
      for ( int j = browp[i]; j < browp[i+1]; j++ ){
	fprintf(fp, "%d %d\n", i + var_offset + ownerRange[mpiRank], 
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
void ApproximateSchur::applyFactor( TACSVec * txvec, TACSVec * tyvec ){
  BVec *xvec, *yvec;
  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

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
      
      TacsScalar * r;
      rvec->getArray(&r);
      memcpy(r, &y[start], (end-start)*sizeof(TacsScalar));
      
      // Solve the global Schur system: S * w = r
      wvec->placeArray(&y[start]);
      inner_ksm->solve(rvec, wvec);
      wvec->restoreArray();
      
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
    fprintf(stderr, "ApproximateSchur type error: Input/output must be BVec\n");
  }
}

/*!
  The block-Jacobi-preconditioned approximate global Schur matrix. 

  This matrix is used within the ApproximateSchur preconditioner.
*/
GlobalSchurMat::GlobalSchurMat( PMat * mat, BCSRMat * _Apc ){
  Apc = _Apc;
  Apc->incref();

  BCSRMat * Aloc;
  mat->getBCSRMat(&Aloc, &Bext);
  Bext->incref();

  mat->getExtColMap(&col_map);
  col_map->incref();

  int bsize, N, Nc;
  mat->getRowMap(&bsize, &N, &Nc);

  varoffset = N-Nc;
  nvars = bsize*Nc;
  rmap = new VarMap(mat->getMPIComm(), Nc, bsize);

  int xsize = col_map->getSize();
  x_ext = new TacsScalar[ xsize ];  
}

GlobalSchurMat::~GlobalSchurMat(){
  Apc->decref();
  Bext->decref(); 
  col_map->decref();

  delete [] x_ext;
}

void GlobalSchurMat::getSize( int * _nr, int * _nc ){
  // Get the local dimensions of the matrix
  *_nr = nvars;
  *_nc = nvars;
}

// Compute y <- A * x 
void GlobalSchurMat::mult( TACSVec * txvec, TACSVec * tyvec ){
  BVec *xvec, *yvec;

  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);  
    
    // Begin sending the external-interface values
    col_map->beginForward(xvec, x_ext, varoffset);
    
    // Finish sending the external-interface unknowns
    col_map->endForward(xvec, x_ext);
    Bext->mult(x_ext, y);
    
    // Apply L^{-1}
    Apc->applyPartialLower(y, varoffset);
    
    // Apply U^{-1}
    Apc->applyPartialUpper(y, varoffset);
    
    // Finish the matrix-vector product
    yvec->axpy(1.0, xvec);
  }
  else {
    fprintf(stderr, "GlobalSchurMat type error: Input/output must be BVec\n");
  }
}

// Compute y <- Bext * xext 
void GlobalSchurMat::multOffDiag( BVec * xvec, BVec * yvec ){  
  TacsScalar *x, *y;
  xvec->getArray(&x);
  yvec->getArray(&y);  

  // Begin sending the external-interface values
  col_map->beginForward(xvec, x_ext, varoffset);

  // Finish sending the external-interface unknowns
  col_map->endForward(xvec, x_ext);
  Bext->mult(x_ext, y);
}

// Return a new Vector
TACSVec * GlobalSchurMat::createVec(){
  return new BVec(rmap);
}
