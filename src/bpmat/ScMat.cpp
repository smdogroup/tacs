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

#include "ScMat.h"
#include "FElibrary.h"
#include "tacslapack.h"

/*
  Schur-complement based matrix
*/

/*!
  A matrix for Schur-complement based preconditioning.

  This matrix stores entries in the form:

  [ B, E ][ x ]   [ f ]
  [ F, C ][ y ] = [ g ]

  where B, E, F and C are block matrices stored in a BCSR format.
  This matrix is re-ordered from the global numbering scheme using the
  BVecDistribute objects b_map and c_map that map the distributed
  global vector to the rows of B and C respectively.

  input:
  rmap:  the variable map for all variables
  B:     the diagonal block BCSR matrices
  E, F:  the off-diagonal block BCSR matrices
  C:     the diagonal block corresponding to the local Schur complement
  b_map: the global variables corresponding to the B-variables
  c_map: the global variables corresponding to the C-variables
*/
ScMat::ScMat( TACSVarMap *_rmap,
              BCSRMat *_B, BCSRMat *_E, BCSRMat *_F, BCSRMat *_C,
              TACSBVecDistribute *_b_map,
              TACSBVecDistribute *_c_map ){
  init(_rmap, _B, _E, _F, _C, _b_map, _c_map);
}

/*
  Create a ScMat class without any data - set everything to NULL
*/
ScMat::ScMat(){
  rmap = NULL;
  B = NULL;
  E = NULL;
  F = NULL;
  C = NULL;
  b_map = NULL;
  c_map = NULL;
  b_ctx = NULL;
  c_ctx = NULL;
  xlocal = ylocal = NULL;
}

/*
  Initialize the data for the ScMat class

  input:
  rmap:  the variable map for all variables
  B:     the diagonal block BCSR matrices
  E, F:  the off-diagonal block BCSR matrices
  C:     the diagonal block corresponding to the local Schur complement
  b_map: the global variables corresponding to the B-variables
  c_map: the global variables corresponding to the C-variables
*/
void ScMat::init( TACSVarMap *_rmap,
                  BCSRMat *_B, BCSRMat *_E, BCSRMat *_F, BCSRMat *_C,
                  TACSBVecDistribute *_b_map,
                  TACSBVecDistribute *_c_map ){
  rmap = _rmap;
  rmap->incref();

  B = _B;  B->incref();
  E = _E;  E->incref();
  F = _F;  F->incref();
  C = _C;  C->incref();

  b_map = _b_map;
  b_map->incref();
  c_map = _c_map;
  c_map->incref();

  // Check that the block dimensions work out
  int bs = B->getBlockSize();
  if (bs != E->getBlockSize() ||
      bs != F->getBlockSize() ||
      bs != C->getBlockSize()){
    fprintf(stderr, "ScMat error: block sizes do not match\n");
    return;
  }

  // Check that things are the correct dimensions
  if (B->getRowDim() != E->getRowDim()){
    fprintf(stderr, "ScMat error: B, E row dimensions do not match\n");
    return;
  }
  if (F->getRowDim() != C->getRowDim()){
    fprintf(stderr, "ScMat error: F, C row dimensions do not match\n");
    return;
  }

  if (B->getColDim() != F->getColDim()){
    fprintf(stderr, "ScMat error: B, F column dimensions do not match\n");
    return;
  }
  if (E->getColDim() != C->getColDim()){
    fprintf(stderr, "ScMat error: E, C column dimensions do not match\n");
    return;
  }

  if (B->getColDim() != b_map->getDim()){
    fprintf(stderr, "ScMat error: b_map dimensions do not \
match dimensions of B\n");
    return;
  }

  if (C->getColDim() != c_map->getDim()){
    fprintf(stderr, "ScMat error: c_map dimensions do not \
match dimensions of C\n");
    return;
  }

  // Allocate the memory for the preconditioning operations
  local_size = bs*(b_map->getDim() + c_map->getDim());
  local_offset = bs*b_map->getDim();

  b_ctx = b_map->createCtx(bs);
  b_ctx->incref();
  c_ctx = c_map->createCtx(bs);
  c_ctx->incref();

  xlocal = new TacsScalar[local_size];
  ylocal = new TacsScalar[local_size];
  memset(xlocal, 0, local_size*sizeof(TacsScalar));
  memset(ylocal, 0, local_size*sizeof(TacsScalar));
}

/*
  The destructor for the ScMat
*/
ScMat::~ScMat(){
  if (rmap){ rmap->decref(); }
  if (B){ B->decref(); }
  if (E){ E->decref(); }
  if (F){ F->decref(); }
  if (C){ C->decref(); }
  if (b_map){ b_map->decref(); }
  if (c_map){ c_map->decref(); }
  if (b_ctx){ b_ctx->decref(); }
  if (c_ctx){ c_ctx->decref(); }
  if (xlocal){ delete [] xlocal; }
  if (ylocal){ delete [] ylocal; }
}

/*
  Get the row/column dimension of the matrix

  output:
  nr:  the row dimension
  nc:  the column dimension
*/
void ScMat::getSize( int *nr, int *nc ){
  int bs = B->getBlockSize();
  *nr = bs*rmap->getDim();
  *nc = bs*rmap->getDim();
}

/*
  Set the matrix to zero
*/
void ScMat::zeroEntries(){
  B->zeroEntries();
  E->zeroEntries();
  F->zeroEntries();
  C->zeroEntries();
}

/*!
  Copy the values from the another matrix

  input:
  mat:  the matrix to copy the values from
*/
void ScMat::copyValues( TACSMat *mat ){
  // Safely down-cast the matrix to an ScMat - returns NULL
  // if this is not possible
  ScMat *smat = dynamic_cast<ScMat*>(mat);
  if (smat){
    B->copyValues(smat->B);
    E->copyValues(smat->E);
    F->copyValues(smat->F);
    C->copyValues(smat->C);
  }
  else {
    fprintf(stderr, "Cannot copy matrices of different types\n");
  }
}

/*!
  Scale the entries in the other matrices by a given scalar

  input:
  alpha:  Scale the matrix by alpha: A <- alpha*A
*/
void ScMat::scale( TacsScalar alpha ){
  B->scale(alpha);
  E->scale(alpha);
  F->scale(alpha);
  C->scale(alpha);
}

/*!
  Compute y <- y + alpha * x
*/
void ScMat::axpy( TacsScalar alpha, TACSMat *mat ){
  ScMat *smat = dynamic_cast<ScMat*>(mat);
  if (smat){
    B->axpy(alpha, smat->B);
    E->axpy(alpha, smat->E);
    F->axpy(alpha, smat->F);
    C->axpy(alpha, smat->C);
  }
  else {
    fprintf(stderr, "Cannot apply axpy to matrices of different types\n");
  }
}

/*!
  Compute y <- alpha * x + beta * y
*/
void ScMat::axpby( TacsScalar alpha, TacsScalar beta, TACSMat *mat ){
  ScMat *smat = dynamic_cast<ScMat*>(mat);
  if (smat){
    B->axpby(alpha, beta, smat->B);
    E->axpby(alpha, beta, smat->E);
    F->axpby(alpha, beta, smat->F);
    C->axpby(alpha, beta, smat->C);
  }
  else {
    fprintf(stderr, "Cannot apply axpby to matrices of different types\n");
  }
}

/*
  Add a scalar to the diagonal components
*/
void ScMat::addDiag( TacsScalar alpha ){
  B->addDiag(alpha);
  C->addDiag(alpha);
}

/*!
  Matrix multiplication

  y = A * x

  where A is stored in the block form:

  A =
  [ B, E ]
  [ F, C ]

  All inter-process communication contributions are contained
  in the sub-block C.

  Procedure:

  Initiate the distribution of the off-process contributions.  These
  are entirely contained within c_map. Next, perform the on-process
  reordering in b_map, and compute the contributions: B*xlocal,
  E*xlocal. Finalize communications with c_map, and compute the C
  mat-vec product. Initiate reverse communication and compute the
  local product with E. Perform the on-process reverse ordering and
  complete the reverse c_map communication.
*/
void ScMat::mult( TACSVec *txvec, TACSVec *tyvec ){
  tyvec->zeroEntries();

  // Safely down-cast the TACSVec vectors to TACSBVecs
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec*>(txvec);
  yvec = dynamic_cast<TACSBVec*>(tyvec);

  if (xvec && yvec){
    TacsScalar *x, *y;
    xvec->getArray(&x);
    yvec->getArray(&y);

    // First begin the communication of x to the local values
    b_map->beginForward(b_ctx, x, xlocal);
    c_map->beginForward(c_ctx, x, &xlocal[local_offset]);
    b_map->endForward(b_ctx, x, xlocal);

    // Perform the matrix-vector multiplication
    B->mult(xlocal, ylocal);
    F->mult(xlocal, &ylocal[local_offset]);
    c_map->endForward(c_ctx, x, &xlocal[local_offset]);

    C->multAdd(&xlocal[local_offset],
               &ylocal[local_offset], &ylocal[local_offset]);

    // Start sending the values back to y
    c_map->beginReverse(c_ctx, &ylocal[local_offset], y, TACS_ADD_VALUES);
    E->multAdd(&xlocal[local_offset], ylocal, ylocal);

    // Finish transmitting the values back to y
    b_map->beginReverse(b_ctx, ylocal, y, TACS_INSERT_VALUES);
    c_map->endReverse(c_ctx, &ylocal[local_offset], y, TACS_ADD_VALUES);
    b_map->endReverse(b_ctx, ylocal, y, TACS_INSERT_VALUES);
  }
  else {
    fprintf(stderr, "ScMat type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrices

  output:
  B, E, F, C: the matrices in the ScMat class
*/
void ScMat::getBCSRMat( BCSRMat **_B, BCSRMat **_E,
                        BCSRMat **_F, BCSRMat **_C ){
  if (_B){ *_B = B; }
  if (_E){ *_E = E; }
  if (_F){ *_F = F; }
  if (_C){ *_C = C; }
}

/*!
  The Global Schur preconditioner. Some assembly required.

  The global Schur complement is formed by first assembling the Schur
  complement contributions on each of the processors.  These
  contributions are then added to a global matrix that is assembled on
  either all processors, or (if use_root = 1) only the root process.

  levFill and fill are the level of fill and fill-in on the global
  Schur complement

  input:
  smat:    the ScMat matrix for the preconditioner
  levFill: the level of fill to use
  fill:    the expected/best estimate of the fill-in factor
  reorder: flag to indicate whether to re-order the global Schur complement
*/
PcScMat::PcScMat( ScMat *_mat, int levFill, double fill,
                  int reorder_schur_complement ){
  mat = _mat;
  mat->incref();

  mat->getBCSRMat(&B, &E, &F, &C);
  B->incref();
  E->incref();
  F->incref();
  C->incref();

  monitor_factor = 0;
  monitor_back_solve = 0;

  // By default use the less-memory intensive option
  use_pdmat_alltoall = 0;

  // Perform the symbolic factorization of the [ B, E; F, C ] matrix
  int use_full_schur = 1; // Use the exact F * B^{-1} * E

  // Get the block size
  int bsize = B->getBlockSize();

  b_map = mat->getLocalMap();
  c_map = mat->getSchurMap();
  b_map->incref();
  c_map->incref();

  // Create the contexts
  c_ctx = c_map->createCtx(bsize);
  b_ctx = b_map->createCtx(bsize);
  c_ctx->incref();
  b_ctx->incref();

  // Symbolically calculate Sc = C - F * B^{-1} * E
  TACSVarMap *rmap = mat->getVarMap();
  MPI_Comm comm = rmap->getMPIComm();
  Bpc = new BCSRMat(comm, B, E, F, C, levFill, fill,
                    &Epc, &Fpc, &Sc, use_full_schur);
  Bpc->incref();
  Epc->incref();
  Fpc->incref();
  Sc->incref();

  // Determine the ordering for the global Schur variables
  // -----------------------------------------------------
  int root = 0;

  int rank, size;
  MPI_Comm_rank(rmap->getMPIComm(), &rank);
  MPI_Comm_size(rmap->getMPIComm(), &size);

  const int *ownerRange;
  rmap->getOwnerRange(&ownerRange);

  // Compute the max proc grid size
  int max_grid_size = size;
  for ( int i = size; i > 0; i-- ){
    if (ownerRange[i] - ownerRange[i-1] > 0){
      max_grid_size = i;
      break;
    }
  }

  // Get the indices of the global variables
  TACSBVecIndices *c_map_indx = c_map->getIndices();

  // Get the number of global Schur variables
  const int *schur_vars;
  num_local_schur_vars = c_map_indx->getIndices(&schur_vars);

  // Gather the size of the indices on each process
  int *schur_count = NULL;
  int *schur_ptr = NULL;
  if (rank == root){
    schur_count = new int[ size ];
    schur_ptr = new int[ size ];
  }

  MPI_Gather(&num_local_schur_vars, 1, MPI_INT,
             schur_count, 1, MPI_INT, root, comm);

  // schur_root is first stored as the list of tacs variable numbers
  // for the schur complement contribution from each processor. After
  // processing on the root node, schur_root is an array of these same
  // variables but now with the order obtained from uniquely sorting
  // the list of the original schur_root variables.
  int *schur_root = NULL;
  int num_schur_root = 0;
  if (rank == root){
    for ( int i = 0; i < size; i++ ){
      schur_ptr[i] = num_schur_root;
      num_schur_root += schur_count[i];
    }
    schur_root = new int[ num_schur_root ];
  }

  // Gather all the global Schur variables to the root process
  MPI_Gatherv((void*)schur_vars, num_local_schur_vars, MPI_INT,
              schur_root, schur_count, schur_ptr,
              MPI_INT, root, comm);

  // Duplicate the global list and uniquify the result
  int *unique_schur = NULL;
  int num_unique_schur = 0;

  if (rank == root){
    unique_schur = new int[ num_schur_root ];
    memcpy(unique_schur, schur_root, num_schur_root*sizeof(int));
    num_unique_schur = FElibrary::uniqueSort(unique_schur, num_schur_root);

    // For each global Schur variable, now assign an output - the index
    // into the unique list of global Schur variables
    for ( int i = 0; i < num_schur_root; i++ ){
      int *item = (int*)bsearch(&schur_root[i], unique_schur,
                                num_unique_schur, sizeof(int),
                                FElibrary::comparator);
      schur_root[i] = item - unique_schur;
    }
  }

  // Broadcast the global size of the unique Schur complement matrix
  MPI_Bcast(&num_unique_schur, 1, MPI_INT, root, comm);

  // Now, pass the unique Schur indices back to the original procs
  // This is an initial ordering of the variables,
  local_schur_vars = new int[ num_local_schur_vars ];
  MPI_Scatterv(schur_root, schur_count, schur_ptr, MPI_INT,
               local_schur_vars, num_local_schur_vars, MPI_INT, root, comm);

  // Retrieve the non-zero pattern from the local Schur complement
  // and pass it into the block-cyclic matrix
  int bs, num_schur_vars;
  const int *rowp, *cols;
  TacsScalar *vals;
  Sc->getArrays(&bs, &num_schur_vars, &num_schur_vars,
                &rowp, &cols, &vals);

  // Determine the size of the global Schur complement
  int M = num_unique_schur, N = num_unique_schur;

  // Determine the number of blocks to use per block-cylic block
  int csr_blocks_per_block = 1;
  if (bsize < 36){
    csr_blocks_per_block = 36/bsize;
  }

  // Create the global block-cyclic Schur complement matrix
  pdmat = new PDMat(comm, M, N, bsize, local_schur_vars,
                    num_schur_vars, rowp, cols, csr_blocks_per_block,
                    reorder_schur_complement, max_grid_size);
  pdmat->incref();

  // Get the information about the reordering/blocks from the matrix
  int nrows, ncols;
  const int *bptr, *xbptr, *perm, *iperm, *orig_bptr;
  pdmat->getBlockPointers(&nrows, &ncols, &bptr, &xbptr,
                          &perm, &iperm, &orig_bptr);
  if (!orig_bptr){
    orig_bptr = bptr;
  }

  // Set the number of local variables that are defined
  int local_var_count = xbptr[nrows]/bsize;

  // Allocate space for the new pointers
  int *tacs_schur_count = NULL;
  int *tacs_schur_ptr = NULL;
  if (rank == root){
    tacs_schur_count = new int[ size ];
    tacs_schur_ptr = new int[ size ];
  }

  // Gather the local variable count on this processor
  MPI_Gather(&local_var_count, 1, MPI_INT,
             tacs_schur_count, 1, MPI_INT, root, comm);

  // Now, reorder the variables in the Schur complement
  int *local_tacs_schur_vars = NULL;
  if (rank == root){
    // Compute the offset into the pdmat order
    int num_schur = 0;
    for ( int i = 0; i < size; i++ ){
      tacs_schur_ptr[i] = num_schur;
      num_schur += tacs_schur_count[i];
    }

    // Find out where to place variable i from the unique list of
    // local schur variables
    local_tacs_schur_vars = new int[ num_schur ];
    for ( int i = 0, j = 0; (i < nrows) && (j < num_unique_schur); i++ ){
      while (j < num_unique_schur &&
             bsize*j >= orig_bptr[i] &&
             bsize*j < orig_bptr[i+1]){
        // Get the re-ordered block
        int block = i;
        if (iperm){ block = iperm[i]; }
        int owner = pdmat->get_block_owner(block, block);

        // Count up the number of local blocks to offset.  This is a
        // double loop which could be avoided in future. This might be
        // a bottle neck for very large cases, but just at set up.
        int index = (bsize*j - orig_bptr[i])/bsize + tacs_schur_ptr[owner];
        for ( int k = 0; k < block; k++ ){
          if (owner == pdmat->get_block_owner(k, k)){
            index += (bptr[k+1] - bptr[k])/bsize;
          }
        }

        // Set the new value of the schur index
        local_tacs_schur_vars[index] = unique_schur[j];
        unique_schur[j] = index;

        // Increment the index
        j++;
      }
    }

    // For each global Schur variable, now assign an output - the index
    // into the unique list of global Schur variables
    for ( int i = 0; i < num_schur_root; i++ ){
      schur_root[i] = unique_schur[schur_root[i]];
    }

    // Free the original set of unique schur variables
    delete [] unique_schur;
  }

  // Send unique_schur back to the owning processes
  int *tacs_schur_vars = new int[ local_var_count ];
  MPI_Scatterv(local_tacs_schur_vars, tacs_schur_count,
               tacs_schur_ptr, MPI_INT,
               tacs_schur_vars, local_var_count, MPI_INT, root, comm);

  int *local_schur = new int[ num_local_schur_vars ];
  MPI_Scatterv(schur_root, schur_count, schur_ptr, MPI_INT,
               local_schur, num_local_schur_vars, MPI_INT, root, comm);

  // Free memory not required anymore
  if (rank == root){
    delete [] tacs_schur_count;
    delete [] tacs_schur_ptr;
    delete [] schur_count;
    delete [] schur_ptr;
    delete [] schur_root;
    delete [] local_tacs_schur_vars;
  }

  // Set up information required for the global Schur complement matrix
  // Set the variable map
  schur_map = new TACSVarMap(comm, local_var_count);
  schur_map->incref();

  // Create the index set for the new Schur complement variables
  TACSBVecIndices *schur_index =
    new TACSBVecIndices(&local_schur, num_local_schur_vars);

  // Create the Schur complement variable distribution object
  schur_dist = new TACSBVecDistribute(schur_map, schur_index);
  schur_dist->incref();

  // Create the Schur complement context
  schur_ctx = schur_dist->createCtx(bsize);
  schur_ctx->incref();

  TACSBVecIndices *tacs_schur_index =
    new TACSBVecIndices(&tacs_schur_vars, local_var_count);

  // Create the index set for the global Schur complement variables
  tacs_schur_dist = new TACSBVecDistribute(mat->getVarMap(),
                                           tacs_schur_index);
  tacs_schur_dist->incref();

  // Create the context for the Schur complement in the global indices
  tacs_schur_ctx = tacs_schur_dist->createCtx(bsize);
  tacs_schur_ctx->incref();

  // Allocate space for local storage of vectors
  int xsize = bsize*b_map->getDim();
  int ysize = bsize*c_map->getDim();
  xlocal = new TacsScalar[ xsize ];
  yinterface = new TacsScalar[ ysize ];

  memset(xlocal, 0, xsize*sizeof(TacsScalar));
  memset(yinterface, 0, ysize*sizeof(TacsScalar));

  // Allocate the Schur complement vectors
  yschur = new TACSBVec(schur_map, bsize);
  gschur = new TACSBVec(schur_map, bsize);
  yschur->incref();
  gschur->incref();
}

/*
  Destructor for the PcScMat preconditioner object
*/
PcScMat::~PcScMat(){
  // Decrease reference counts to the matrices
  mat->decref();
  B->decref();
  E->decref();
  F->decref();
  C->decref();

  // Decrease reference counts for the corresponding preconditioner
  // objects associated with the matrices
  Bpc->decref();
  Epc->decref();
  Fpc->decref();
  Sc->decref();

  // Decrease the reference count to the variable maps
  c_map->decref();
  b_map->decref();
  schur_map->decref();
  schur_dist->decref();
  tacs_schur_dist->decref();

  // Free the contexts
  c_ctx->decref();
  b_ctx->decref();
  schur_ctx->decref();
  tacs_schur_ctx->decref();

  // Free the global Schur complement matrix
  pdmat->decref();

  // Deallocate the local vectors
  gschur->decref();
  yschur->decref();
  delete [] xlocal;
  delete [] yinterface;
}

/*
  Test three different methods for calculating the Schur complement
  contributions.

  This test ensures that the methods are consistent for the given
  matrix. Note that this test can only be performed after the
  preconditioner has been factored since the matrix Sc is not
  populated until this time.
*/
void PcScMat::testSchurComplement( TACSVec *tin, TACSVec *tout ){
  TACSBVec *invec, *outvec;
  invec = dynamic_cast<TACSBVec*>(tin);
  outvec = dynamic_cast<TACSBVec*>(tout);

  if (invec && outvec){
    // Test two methods of computing the effect of the Schur complement
    outvec->zeroEntries();

    // Get the array
    TacsScalar *in, *out;
    invec->getArray(&in);
    outvec->getArray(&out);

    // Allocate a temporary array to store c-entries
    int bsize = B->getBlockSize();
    int c_size = bsize*c_map->getDim();
    TacsScalar *temp = new TacsScalar[c_size];

    // Comput the schur complement product
    c_map->beginForward(c_ctx, in, yinterface);
    c_map->endForward(c_ctx, in, yinterface);
    Sc->mult(yinterface, temp);
    c_map->beginReverse(c_ctx, temp, out, TACS_ADD_VALUES);
    c_map->endReverse(c_ctx, temp, out, TACS_ADD_VALUES);

    delete [] temp;

    // Compute the schur complement product a second way
    TacsScalar *y;
    yschur->getArray(&y);
    schur_dist->beginReverse(schur_ctx, yinterface,
                             y, TACS_INSERT_VALUES);
    schur_dist->endReverse(schur_ctx, yinterface,
                           y, TACS_INSERT_VALUES);

    TacsScalar *g = NULL;
    yschur->getArray(&y);
    pdmat->mult(y, g);

    TacsScalar outnorm = outvec->norm();
    TacsScalar gnorm = gschur->norm();

    // Now collect the elements directly using the tacs_schur_dist object
    tacs_schur_dist->beginForward(tacs_schur_ctx, in, y);
    tacs_schur_dist->endForward(tacs_schur_ctx, in, y);

    pdmat->mult(y, g);
    TacsScalar gnorm2 = gschur->norm();

    int rank;
    MPI_Comm_rank(c_map->getMPIComm(), &rank);

    if (rank == 0){
      printf("Schur complement consistency test: \n");
      printf("|Full matrix|     = %25.15e \n", TacsRealPart(outnorm));
      printf("|Local to Schur|  = %25.15e \n", TacsRealPart(gnorm));
      printf("|Global to Schur| = %25.15e \n", TacsRealPart(gnorm2));
    }
  }
  else {
    fprintf(stderr, "PcScMat type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrices

  output:
  B, E, F, C: the matrices in the ScMat class
*/
void PcScMat::getBCSRMat( BCSRMat **_Bpc, BCSRMat **_Epc,
                          BCSRMat **_Fpc, BCSRMat **_Sc ){
  if (_Bpc){ *_Bpc = Bpc; }
  if (_Epc){ *_Epc = Epc; }
  if (_Fpc){ *_Fpc = Fpc; }
  if (_Sc){ *_Sc = Sc; }
}

/*
  Set the flag that prints out the factorization time

  input:
  flag: the flag value for the factor-time monitor
*/
void PcScMat::setMonitorFactorFlag( int flag ){
  monitor_factor = flag;
}

/*
  Set the flag that prints out the back solve time

  input:
  flag:  the flag value for the back-solve monitor
*/
void PcScMat::setMonitorBackSolveFlag( int flag ){
  monitor_back_solve = flag;
}

/*
  Set the flag that controls which matrix assembly code to use.

  This flag controls whether the Alltoall version of the PDMat matrix
  assembly is used. When true, this uses a faster, but more
  memory-intensive version of the code. Be aware that this can cause
  the code to run out of memory, but can be very beneficial in terms
  of CPU time.

  input:
  flag:  the flag value to use for the Alltoall flag
*/
void PcScMat::setAlltoallAssemblyFlag( int flag ){
  use_pdmat_alltoall = flag;
}

/*
  Factor the Schur-complement based preconditioner

  This proceeds as follows:

  1. Copy the values of the original matrix in the diagonal blocks B
  and C and the off-diagonal blocks E and F.  These are copied into
  Bpc, Sc and Epc and Fpc respectively.

  2. Factor the diagonal block Bpc = Lb*Ub, and apply it to the
  off-diagonals:

  Epc <- Lb^{-1} Epc
  Fpc <- Fpc Ub^{-1}

  3. Compute the Schur complement contribution on all blocks:
  Sc <- (Sc - Fpc * Epc)

  4. Assemble the contributions of Sc on all processors into the
  global Schur complement matrix scmat.

  Factor the preconditioner for this matrix (pc).
*/
void PcScMat::factor(){
  // Set the time variables
  double diag_factor_time = 0.0;
  double schur_complement_time = 0.0;
  double global_schur_assembly = 0.0;
  double global_schur_time = 0.0;

  if (monitor_factor){
    diag_factor_time = -MPI_Wtime();
  }

  // Copy the diagonal matrix B and factor it
  Bpc->copyValues(B);
  Bpc->factor();

  if (monitor_factor){
    diag_factor_time += MPI_Wtime();
    schur_complement_time = -MPI_Wtime();
  }

  // Copy C, E and F and apply B to obtain L^{-1}*E and F*U^{-1}
  Sc->copyValues(C);
  Epc->copyValues(E);
  Fpc->copyValues(F);
  Bpc->applyLowerFactor(Epc);
  Bpc->applyUpperFactor(Fpc);

  // Compute the Schur complement matrix Sc
  Sc->matMultAdd(-1.0, Fpc, Epc);

  if (monitor_factor){
    schur_complement_time += MPI_Wtime();
    global_schur_assembly = -MPI_Wtime();
  }

  // Assemble the global Schur complement system into PDMat
  // First, zero the Schur complement matrix
  pdmat->zeroEntries();

  // Retrieve the local arrays for the local Schur complement
  int bsize, mlocal, nlocal;
  const int *rowp, *cols;
  TacsScalar *scvals;
  Sc->getArrays(&bsize, &mlocal, &nlocal,
                &rowp, &cols, &scvals);

  // Add the values into the global Schur complement matrix
  // using either the alltoall approach or a sequential add values
  // approach that uses less memory
  if (use_pdmat_alltoall){
    pdmat->addAlltoallValues(bsize, mlocal, local_schur_vars,
                             rowp, cols, scvals);
  }
  else {
    pdmat->addAllValues(bsize, mlocal, local_schur_vars,
                        rowp, cols, scvals);
  }

  if (monitor_factor){
    global_schur_assembly += MPI_Wtime();
    global_schur_time = -MPI_Wtime();
  }

  // Factor the global Schur complement
  pdmat->factor();

  if (monitor_factor){
    global_schur_time += MPI_Wtime();

    int rank;
    MPI_Comm_rank(b_map->getMPIComm(), &rank);
    printf("[%d] Diagonal factor time:  %8.4f\n", rank, diag_factor_time);
    printf("[%d] Local Schur time:      %8.4f\n", rank, schur_complement_time);
    printf("[%d] Global Schur assembly: %8.4f\n", rank, global_schur_assembly);
    printf("[%d] Global Schur time:     %8.4f\n", rank, global_schur_time);
  }
}

/*!
  Apply the preconditioner to the input vector

  [ B,  E ][ x ]   [ f ]
  [ F,  C ][ y ] = [ g ]

  The following relationships hold:

  x = B^{-1}(f - E y)
  F B^{-1} (f - E y) + C y = g

  Collecting these, results in the Schur complement system:

  ( C - F B^{-1} E) y = g - F * B^{-1} f

  Solve this equation for the interface unknowns, then solve,

  x = B^{-1} (f - E y)

  for the remaining unknowns.

  ------------------------------------
  The individual matrices are given as follows,

  Bpc->applyFactor(x, y) -> y = Ub^{-1} Lb^{-1} x

  [ B, E ]   [ Lb       , 0  ][ Ub,  Lb^{-1} E ]
  [ F, C ] = [ F Ub^{-1}, Lc ][  0,  Uc        ]

  Epc = Lb^{-1} E
  Fpc = F Ub^{-1}

  1. Compute x = L^{-1} f
  2. Compute g' = g - F U^{-1} x = g - Fpc * x
  3. Solve approximately (C - F B^{-1} E) y = g'
  4. Compute x <- U^{-1} (x - L^{-1} E * y) = U^{-1} (x - Epc * y)
*/
void PcScMat::applyFactor( TACSVec *tin, TACSVec *tout ){
  // First, perform a safe down-cast from TACSVec to BVec
  TACSBVec *invec, *outvec;
  invec = dynamic_cast<TACSBVec*>(tin);
  outvec = dynamic_cast<TACSBVec*>(tout);

  if (invec && outvec){
    // Set the variables for the back-solve monitor
    double local_time = 0.0;
    double schur_time = 0.0;

    if (monitor_back_solve){
      local_time = -MPI_Wtime();
    }

    // Get the input and output arrays
    TacsScalar *in, *out;
    invec->getArray(&in);
    outvec->getArray(&out);

    // Pass g to the global Schur complement
    TacsScalar *g = NULL;
    gschur->getArray(&g);
    tacs_schur_dist->beginForward(tacs_schur_ctx, in, g);

    // Re-order the local variables into xlocal. Note that this is a
    // local-only reordering and does not require communication.
    b_map->beginForward(b_ctx, in, xlocal);
    b_map->endForward(b_ctx, in, xlocal);

    // xlocal = L^{-1} f
    Bpc->applyLower(xlocal, xlocal);

    // yinterface = F U^{-1} xlocal = F U^{-1} L^{-1} f
    Fpc->mult(xlocal, yinterface);

    // Transmit yinterface to the Schur complement system
    // Pass F U^{-1} L^{-1} f to yschur
    yschur->zeroEntries();
    TacsScalar *y = NULL;
    yschur->getArray(&y);
    schur_dist->beginReverse(schur_ctx, yinterface, y, TACS_ADD_VALUES);

    // Finish accepting g from the global input vector
    tacs_schur_dist->endForward(tacs_schur_ctx, in, g);
    schur_dist->endReverse(schur_ctx, yinterface, y, TACS_ADD_VALUES);

    // Compute the right hand side: g - F U^{-1} L^{-1} f
    gschur->axpy(-1.0, yschur);

    if (monitor_back_solve){
      schur_time = -MPI_Wtime();
    }

    // Apply the global Schur complement factorization to the right
    // hand side
    pdmat->applyFactor(g);

    if (monitor_back_solve){
      schur_time += MPI_Wtime();
      local_time -= schur_time;
    }

    // copy the values from gschur to yschur
    yschur->copyValues(gschur);

    // Pass yschur solution back to the local variables
    schur_dist->beginForward(schur_ctx, y, yinterface);

    // Pass yschur also to the global variables
    tacs_schur_dist->beginReverse(tacs_schur_ctx, y, out,
                                  TACS_INSERT_VALUES);

    // Compute yinterface = yinterface - L^{-1} E y
    // Note: scale y, by -1 first
    schur_dist->endForward(schur_ctx, y, yinterface);
    int one = 1;
    int len = Bpc->getBlockSize()*c_map->getDim();
    TacsScalar alpha = -1.0;
    BLASscal(&len, &alpha, yinterface, &one);

    // Compute xlocal = xlocal - L^{-1} E * yinterface
    Epc->multAdd(yinterface, xlocal, xlocal);

    // Compute xlocal = U^{-1} xlocal
    Bpc->applyUpper(xlocal, xlocal);

    b_map->beginReverse(b_ctx, xlocal, out, TACS_INSERT_VALUES);
    b_map->endReverse(b_ctx, xlocal, out, TACS_INSERT_VALUES);

    // Finish inserting the Schur complement values
    tacs_schur_dist->endReverse(tacs_schur_ctx, y, out,
                                TACS_INSERT_VALUES);

    if (monitor_back_solve){
      local_time += MPI_Wtime();

      int rank;
      MPI_Comm_rank(b_map->getMPIComm(), &rank);
      printf("[%d] Local back-solve time:  %8.4f\n", rank, local_time);
      printf("[%d] Schur back-solve time:  %8.4f\n", rank, schur_time);
    }
  }
  else {
    fprintf(stderr, "PcScMat type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrix
*/
void PcScMat::getMat( TACSMat **_mat ){
  *_mat = mat;
}
