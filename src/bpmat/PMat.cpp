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

#include <stdio.h>
#include "PMat.h"
#include "FElibrary.h"
#include "tacslapack.h"

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

/*
  The protected constructor that does not take any arguments.
*/
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
  Matrix multiplication
*/
void TACSPMat::multTranspose( TACSVec *txvec, TACSVec *tyvec ){
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);

    Bext->multTranspose(&x[ext_offset], x_ext);
    Aloc->multTranspose(x, y);

    ext_dist->beginReverse(ctx, x_ext, y, TACS_ADD_VALUES);
    ext_dist->endReverse(ctx, x_ext, y, TACS_ADD_VALUES);
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

/*
  Get the local number of rows/coupling rows in the matrix
*/
void TACSPMat::getRowMap( int *_bs, int *_N, int *_Nc ){
  if (_bs){ *_bs = bsize; }
  if (_Nc){ *_Nc = Nc; }
  if (_N){ *_N = N; }
}

/*
  Get the number of columns in the matrix
*/
void TACSPMat::getColMap( int *_bs, int *_M ){
  if (_bs){ *_bs = bsize; }
  if (_M){ *_M = N; }
}

/*
  Get the distribute object to distribute values to other processors
*/
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
                                  int _symmetric,
                                  int _use_l1_gauss_seidel ){
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
  use_l1_gauss_seidel = _use_l1_gauss_seidel;
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
  TacsScalar *diag = NULL;

  if (use_l1_gauss_seidel){
    // Get the block and off-block matrices
    BCSRMat *A, *B;
    mat->getBCSRMat(&A, &B);

    // Get the number of variables in the row
    int bsize, N, Nc;
    mat->getRowMap(&bsize, &N, &Nc);
    const int var_offset = N - Nc;
    const int b2 = bsize*bsize;

    // Allocate space for the diagonal entries
    diag = new TacsScalar[ bsize*N ];
    memset(diag, 0, bsize*N*sizeof(TacsScalar));

    // Sum up the entries in the off diagonal part of the matrix
    BCSRMatData *Bdata = B->getMatData();

    // Compute the sum of the absolute value of the contributions from
    // the off-processor contributions to the matrix
    for ( int i = var_offset; i < N; i++ ){
      int ib = i - var_offset;

      // Set the pointer into the additional diagonal entries
      TacsScalar *d = &diag[bsize*i];

      // Add the contribution from the local part of the matrix
      for ( int jp = Bdata->rowp[ib]; jp < Bdata->rowp[ib+1]; jp++ ){
        // Set the pointer
        const TacsScalar *a = &Bdata->A[b2*jp];

        // For each row, add up the absolute sum
        for ( int ii = 0; ii < bsize; ii++ ){
          for ( int jj = 0; jj < bsize; jj++ ){
            d[ii] += fabs(TacsRealPart(a[0]));
            a++;
          }
        }
      }
    }

    // Loop through the diagonal entries of the A-matrix and ensure that
    // they are positive
    BCSRMatData *Adata = A->getMatData();

    // Loop over the diagoal entries for the off-diagonal part
    for ( int i = var_offset; i < N; i++ ){
      int jp = 0;
      if (Adata->diag){
        jp = Adata->diag[i];
      }
      else {
        for ( jp = Adata->rowp[i]; jp < Adata->rowp[i+1]; jp++ ){
          if (Adata->cols[jp] == i){
            break;
          }
        }
      }

      // Check if any of the diagonal entries are actually negative
      const TacsScalar *a = &Adata->A[b2*jp];
      for ( int ii = 0; ii < bsize; ii++ ){
        if (TacsRealPart(a[ii*(bsize+1)]) < 0.0){
          diag[bsize*i + ii] *= -1.0;
        }
      }
    }
  }

  // Factor the diagonal
  Aloc->factorDiag(diag);

  if (diag){
    delete [] diag;
  }
}

/*!
  Apply the preconditioner to the input vector.

  This code applies multiple SOR steps to the system A*y = x.
  On each processor, the system of equations is

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

    // Get the number of variables in the row map
    int bsize, N, Nc;
    mat->getRowMap(&bsize, &N, &Nc);

    // Set the start/end values for the A-matrix
    const int start = 0;
    const int end = N - Nc;
    const int offset = N - Nc;

    // Set the start/end values for the B-matrix
    const int bstart = end;
    const int bend = N;

    if (symmetric){
      if (zero_guess){
        yvec->zeroEntries();
        Aloc->applySOR(x, y, omega, 1);
      }
      else {
        // Begin sending the external-interface values
        ext_dist->beginForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(NULL, start, end, offset, omega,
                       x, NULL, y);

        // Finish sending the external-interface unknowns
        ext_dist->endForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(Bext, bstart, bend, offset, omega,
                       x, yext, y);
      }

      // Reverse the smoother
      Aloc->applySOR(Bext, bend, bstart, offset, omega,
                     x, yext, y);

      ext_dist->beginForward(ctx, y, yext);

      // Apply the smoother to the local part of the matrix
      Aloc->applySOR(NULL, end, start, offset, omega,
                     x, NULL, y);

      // Finish sending the external-interface unknowns
      ext_dist->endForward(ctx, y, yext);

      for ( int i = 1; i < iters; i++ ){
        // Begin sending the external-interface values
        ext_dist->beginForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(NULL, start, end, offset, omega,
                       x, NULL, y);

        // Finish sending the external-interface unknowns
        ext_dist->endForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(Bext, bstart, bend, offset, omega,
                       x, yext, y);

        // Reverse the smoother
        Aloc->applySOR(Bext, bend, bstart, offset, omega,
                       x, yext, y);

        ext_dist->beginForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(NULL, end, start, offset, omega,
                       x, NULL, y);

        // Finish sending the external-interface unknowns
        ext_dist->endForward(ctx, y, yext);
      }
    }
    else {
      if (zero_guess){
        yvec->zeroEntries();
        Aloc->applySOR(x, y, omega, 1);
      }
      else {
        // Begin sending the external-interface values
        ext_dist->beginForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(NULL, start, end, offset, omega,
                       x, NULL, y);

        // Finish sending the external-interface unknowns
        ext_dist->endForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(Bext, bstart, bend, offset, omega,
                       x, yext, y);
      }

      for ( int i = 1; i < iters; i++ ){
        // Begin sending the external-interface values
        ext_dist->beginForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(NULL, start, end, offset, omega,
                       x, NULL, y);

        // Finish sending the external-interface unknowns
        ext_dist->endForward(ctx, y, yext);

        // Apply the smoother to the local part of the matrix
        Aloc->applySOR(Bext, bstart, bend, offset, omega,
                       x, yext, y);
      }
    }
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

/*
  Create the Chebyshev smoother object
*/
TACSChebyshevSmoother::TACSChebyshevSmoother( TACSPMat *_mat, int _degree,
                                              double _lower_factor,
                                              double _upper_factor,
                                              int _iters ){
  mat = _mat;
  mat->incref();

  // Create the vectors that are needed in the computation
  TACSVec *tt = mat->createVec();
  TACSVec *ht = mat->createVec();
  TACSVec *rest = mat->createVec();

  // Convert the vectors to TACSBVecs
  t = dynamic_cast<TACSBVec*>(tt);
  h = dynamic_cast<TACSBVec*>(ht);
  res = dynamic_cast<TACSBVec*>(rest);
  t->incref();
  h->incref();
  res->incref();

  // Set the lower/upper factors
  lower_factor = _lower_factor;
  upper_factor = _upper_factor;

  // Set the numger of iterations
  iters = _iters;

  // Compute the order
  degree = _degree;
  alpha = beta = 0.0;
  r = new double[ degree ];
  c = new double[ degree+1 ];

  // Set default values
  for ( int i = 0; i < degree; i++ ){
    r[i] = 0.0;
  }
  for ( int i = 0; i < degree+1; i++ ){
    c[i] = 1.0;
  }
}

/*
  Free the memory associated with the Chebyshev smoother
*/
TACSChebyshevSmoother::~TACSChebyshevSmoother(){
  mat->decref();
  delete [] r;
  delete [] c;
  h->decref();
  t->decref();
  res->decref();
}

/*
  Factor the smoother.

  This involves first estimating the spectrum of the matrix and
  then applying the smoother.
*/
void TACSChebyshevSmoother::factor(){
  // Compute the approximate maximum eigenvalue of the matrix
  alpha = 0.0;
  beta = 0.0;

  double rho = gershgorin();
  // double rho = arnoldi(15);

  // Compute the factor
  alpha = lower_factor*rho;
  beta = upper_factor*rho;

  // Compute the roots of the characteristic polynomial
  for ( int k = 0; k < degree; k++ ){
    r[k] = cos(M_PI*(0.5 + k)/degree);
  }

  // Scale the roots from the interval [-1, 1]
  // to the correct interval [alpha, beta]
  for ( int k = 0; k < degree; k++ ){
    r[k] = 0.5*(beta - alpha)*(r[k] + 1.0) + alpha;
  }

  // Compute the coefficients of the polynomial
  memset(c, 0, (degree+1)*sizeof(double));
  c[0] = 1.0;
  for ( int j = 0; j < degree; j++ ){
    for ( int k = j; k >= 0; k-- ){
      c[k+1] = c[k+1] - r[j]*c[k];
    }
  }

  // Now scale the coefficients so that q(0) = 1.0
  // Note that q(A) = 1 - p(A)*A
  for ( int k = 0; k < degree; k++ ){
    c[k] = c[k]/c[degree];
  }
  c[degree] = 1.0;
}

/*
  Apply the Chebyshev smoother to the preconditioner
*/
void TACSChebyshevSmoother::applyFactor( TACSVec *tx,
                                         TACSVec *ty ){
  // Covert to TACSBVec objects
  TACSBVec *x, *y;
  x = dynamic_cast<TACSBVec*>(tx);
  y = dynamic_cast<TACSBVec*>(ty);

  if (x && y){
    for ( int i = 0; i < iters; i++ ){
      // Compute the initial residual
      // res := x - A*y
      res->copyValues(x);
      mat->mult(y, t);
      res->axpy(-1.0, t);

      // h = c[0]*res
      h->copyValues(res);
      h->scale(-c[0]);

      // Compute the polynomial smoother
      for ( int j = 0; j < degree-1; j++ ){
        // t = A*h
        mat->mult(h, t);

        // h = c[j+1]*res + t = c[j+1]*res + A*h
        h->copyValues(res);
        h->scale(-c[j+1]);
        h->axpy(1.0, t);
      }

      // y <- y + h
      y->axpy(1.0, h);
    }
  }
}

/*
  Retrieve the matrix associated with the smoother
*/
void TACSChebyshevSmoother::getMat( TACSMat **_mat ){
  *_mat = mat;
}

/*
  Estimate the spectral radius using Gershgorin disks
*/
double TACSChebyshevSmoother::gershgorin(){
  double rho = 0.0;

  // Get the matrices
  BCSRMat *A, *B;
  mat->getBCSRMat(&A, &B);
  BCSRMatData *Adata = A->getMatData();
  BCSRMatData *Bdata = B->getMatData();

  // Compute the local sizes
  const int bsize = Adata->bsize;
  const int b2 = bsize*bsize;

  // Get the number of variables in the row map
  int N, Nc;
  mat->getRowMap(NULL, &N, &Nc);
  int var_offset = N - Nc;

  // Allocate space to store the
  double *eig = new double[ bsize ];

  // Use Gershgorin theorem to estimate the max eigenvalue
  for ( int i = 0; i < Adata->nrows; i++ ){
    memset(eig, 0, bsize*sizeof(double));

    // Add the contribution from the local part of the matrix
    for ( int jp = Adata->rowp[i]; jp < Adata->rowp[i+1]; jp++ ){
      int j = Adata->cols[jp];

      // Set the pointer
      const TacsScalar *a = &Adata->A[b2*jp];

      // This is the diagonal
      if (i == j){
        // For each row, compute the max eigenvalue
        for ( int ii = 0; ii < bsize; ii++ ){
          for ( int jj = 0; jj < bsize; jj++ ){
            if (ii == jj){
              eig[ii] += TacsRealPart(a[0]);
            }
            else {
              eig[ii] += fabs(TacsRealPart(a[0]));
            }
            a++;
          }
        }
      }
      else {
        // For each row, add up the absolute sum from each entry
        for ( int ii = 0; ii < bsize; ii++ ){
          for ( int jj = 0; jj < bsize; jj++ ){
            eig[ii] += fabs(TacsRealPart(a[0]));
            a++;
          }
        }
      }
    }

    if (i >= var_offset){
      const int ib = i - var_offset;

      // Add the contribution from the local part of the matrix
      for ( int jp = Bdata->rowp[ib]; jp < Bdata->rowp[ib+1]; jp++ ){
        // Set the pointer
        const TacsScalar *a = &Bdata->A[b2*jp];

        // For each row, add up the absolute sum
        for ( int ii = 0; ii < bsize; ii++ ){
          for ( int jj = 0; jj < bsize; jj++ ){
            eig[ii] += fabs(TacsRealPart(a[0]));
            a++;
          }
        }
      }
    }

    // Update beta
    for ( int ii = 0; ii < bsize; ii++ ){
      if (eig[ii] > rho){
        rho = eig[ii];
      }
    }
  }

  delete [] eig;

  // Perform an Allreduce across all processors
  double temp = 0.0;
  MPI_Allreduce(&rho, &temp, 1, MPI_DOUBLE, MPI_MAX, mat->getMPIComm());

  return temp;
}

/*
  Estimate the spectral radius using Arnoldi
*/
double TACSChebyshevSmoother::arnoldi( int size ){
  double *H = new double[ size*(size+1) ];
  memset(H, 0, size*(size+1)*sizeof(double));

  // Allocate space for the vectors
  TACSVec **W = new TACSVec*[ size+1 ];

  // Create an initial random vector
  W[0] = mat->createVec();
  W[0]->incref();
  W[0]->setRand(-1.0, 1.0);
  W[0]->scale(1.0/W[0]->norm());

  // Apply the boundary conditions
  for ( int i = 0; i < size; i++ ){
    W[i+1] = mat->createVec();
    W[i+1]->incref();

    // Multiply by the matrix to get the next vector
    mat->mult(W[i], W[i+1]);

    // Orthogonalize against the existing subspace
    for ( int j = 0; j <= i; j++ ){
      int index = j + i*(size+1);
      H[index] = TacsRealPart(W[i+1]->dot(W[j]));
      W[i+1]->axpy(-H[index], W[j]);
    }

    // Add the term to the matrix
    int index = i+1 + i*(size+1);
    H[index] = TacsRealPart(W[i+1]->norm());
    W[i+1]->scale(1.0/H[index]);
  }

  // Allocate space for the real/complex eigenvalue
  double *eigreal = new double[ size ];
  double *eigimag = new double[ size ];

  // Compute the eigenspectrum of the Hessenberg matrix
  int lwork = 4*size;
  double *work = new double[ lwork ];
  int ldv = 1;
  int ldh = size+1;
  int info = 0;
  LAPACKdgeev("N", "N", &size, H, &ldh, eigreal, eigimag,
              NULL, &ldv, NULL, &ldv, work, &lwork, &info);

  // Find the maximum absolute eigenvalue
  double rho = 0.0;
  for ( int i = 0; i < size; i++ ){
    double val = sqrt(eigreal[i]*eigreal[i] + eigimag[i]*eigimag[i]);
    if (val > rho){
      rho = val;
    }
  }

  for ( int i = 0; i < size+1; i++ ){
    W[i]->decref();
  }

  delete [] eigreal;
  delete [] eigimag;
  delete [] work;
  delete [] H;
  delete [] W;

  return rho;
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

/*
  Add a diagonal contribution to the preconditioner
*/
void TACSApproximateSchur::setDiagShift( TacsScalar _alpha ){
  alpha = _alpha;
}

/*
  Set a monitor for the inner Krylov method
*/
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

/*
  Free the information associated with the global Schur complement matrix
*/
TACSGlobalSchurMat::~TACSGlobalSchurMat(){
  Apc->decref();
  Bext->decref();
  ext_dist->decref();
  ctx->decref();
  delete [] x_ext;
}

/*
  Get the size of the number of local rows/columns
*/
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

/*
  Return a new Vector
*/
TACSVec *TACSGlobalSchurMat::createVec(){
  return new TACSBVec(rmap, Apc->getBlockSize());
}
