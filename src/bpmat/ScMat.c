#include "ScMat.h"
#include "FElibrary.h"
#include "tacslapack.h"

/*
  Schur-complement based matrix

  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
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
ScMat::ScMat( VarMap * _rmap,
	      BCSRMat * _B, BCSRMat * _E, BCSRMat * _F, BCSRMat * _C,
	      BVecDistribute * _b_map, 
	      BVecDistribute * _c_map ){
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
void ScMat::init( VarMap * _rmap,
		  BCSRMat * _B, BCSRMat * _E, BCSRMat * _F, BCSRMat * _C,
                  BVecDistribute * _b_map, 
                  BVecDistribute * _c_map ){
  rmap = _rmap;
  rmap->incref();

  B = _B;  B->incref();
  E = _E;  E->incref();
  F = _F;  F->incref();
  C = _C;  C->incref();

  b_map = _b_map;
  c_map = _c_map;

  b_map->incref();
  c_map->incref();

  // Check that the block dimensions work out 
  int bs = B->getBlockSize();
  if (bs != E->getBlockSize() ||
      bs != F->getBlockSize() ||
      bs != C->getBlockSize()){
    fprintf(stderr, "ScMat error: block sizes do not match\n" );
    return;
  }

  // Check that things are the correct dimensions
  if (B->getRowDim() != E->getRowDim()){
    fprintf(stderr, "ScMat error: B, E row dimensions do not match\n" );
    return;
  }
  if (F->getRowDim() != C->getRowDim()){
    fprintf(stderr, "ScMat error: F, C row dimensions do not match\n" );
    return;
  }

  if (B->getColDim() != F->getColDim()){
    fprintf(stderr, "ScMat error: B, F column dimensions do not match\n" );
    return;
  }
  if (E->getColDim() != C->getColDim()){
    fprintf(stderr, "ScMat error: E, C column dimensions do not match\n" );
    return;
  }

  if (B->getColDim() != b_map->getDim()){
    fprintf(stderr, "ScMat error: b_map dimensions do not \
match dimensions of B\n" );
    return;
  }

  if (C->getColDim() != c_map->getDim()){
    fprintf(stderr, "ScMat error: c_map dimensions do not \
match dimensions of C\n" );
    return;
  }
  
  // Allocate the memory for the preconditioning operations
  local_size = bs*(b_map->getDim() + c_map->getDim());
  local_offset = bs*b_map->getDim();

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
  if (xlocal){ delete [] xlocal; }
  if (ylocal){ delete [] ylocal; }
}

/*
  Get the row/column dimension of the matrix

  output:
  nr:  the row dimension
  nc:  the column dimension
*/
void ScMat::getSize( int * nr, int * nc ){
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
void ScMat::copyValues( TACSMat * mat ){
  // Safely down-cast the matrix to an ScMat - returns NULL
  // if this is not possible
  ScMat * smat = dynamic_cast<ScMat*>(mat);
  if (smat ){
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
void ScMat::axpy( TacsScalar alpha, TACSMat * mat ){
  ScMat * smat = dynamic_cast<ScMat*>(mat);
  if (smat ){
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
void ScMat::axpby( TacsScalar alpha, TacsScalar beta, TACSMat * mat ){
  ScMat * smat = dynamic_cast<ScMat*>(mat);
  if (smat ){
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
void ScMat::mult( TACSVec * txvec, TACSVec * tyvec ){
  tyvec->zeroEntries();

  // Safely down-cast the TACSVec vectors to BVecs
  BVec *xvec, *yvec;
  xvec = dynamic_cast<BVec*>(txvec);
  yvec = dynamic_cast<BVec*>(tyvec);

  if (xvec && yvec){
    // First begin the communication of x to the local values
    c_map->beginForward(xvec, &xlocal[local_offset]);
    b_map->beginForward(xvec, xlocal);
    b_map->endForward(xvec, xlocal);
    
    // Perform the matrix-vector multiplication 
    B->mult(xlocal, ylocal);
    F->mult(xlocal, &ylocal[local_offset]);
    
    c_map->endForward(xvec, &xlocal[local_offset]);
    
    C->multAdd(&xlocal[local_offset], 
               &ylocal[local_offset], &ylocal[local_offset]);
    
    // Start sending the values back to y
    c_map->beginReverse(&ylocal[local_offset], yvec, BVecDistribute::ADD);
    
    E->multAdd(&xlocal[local_offset], ylocal, ylocal);

    // Finish transmitting the values back to y
    b_map->beginReverse(ylocal, yvec, BVecDistribute::INSERT);
    c_map->endReverse(&ylocal[local_offset], yvec, BVecDistribute::ADD);
    b_map->endReverse(ylocal, yvec, BVecDistribute::INSERT);
  }
  else {
    fprintf(stderr, "ScMat type error: Input/output must be BVec\n");
  }
}

/*
  Retrieve the underlying matrices 

  output:
  B, E, F, C: the matrices in the ScMat class
*/
void ScMat::getBCSRMat( BCSRMat ** _B, BCSRMat ** _E,
			BCSRMat ** _F, BCSRMat ** _C ){
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
PcScMat::PcScMat( ScMat * smat, int levFill, double fill,
                  int reorder_schur_complement ){
  smat->getBCSRMat(&B, &E, &F, &C);
  B->incref(); 
  E->incref();
  F->incref(); 
  C->incref();

  monitor_factor = 0;
  monitor_back_solve = 0;

  use_pdmat_alltoall = 0; // By default use the less-memory intensive option 
  // (but this may be slower!)

  // Perform the symbolic factorization of the [ B, E; F, C ] matrix
  int use_full_schur = 1; // Use the exact F * B^{-1} * E 
  // int use_full_schur = 0; // Use an approximation - discarding elements

  b_map = smat->getLocalMap();
  b_map->incref();

  c_map = smat->getSchurMap();
  c_map->incref();

  // Symbolically calculate Sc = C - F * B^{-1} * E
  VarMap * rmap = smat->getVarMap();
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

  const int * ownerRange;
  int rank, size;
  rmap->getOwnerRange(&ownerRange, &rank, &size);

  // Get the indices of the global variables
  BVecIndices * c_map_indx = c_map->getBVecIndices();
  int * global_schur_vars;
  int nglobal_schur = c_map_indx->getIndices(&global_schur_vars);
  int * new_global_schur_vars = new int[ nglobal_schur ];

  // Gather the size of the indices on each process
  int * global_schur_count = NULL;
  int * global_schur_ptr = NULL;
  if (rank == root){
    global_schur_count = new int[ size ];
    global_schur_ptr = new int[ size ];
  }

  MPI_Gather(&nglobal_schur, 1, MPI_INT, 
             global_schur_count, 1, MPI_INT, root, comm);

  // Set up the ptr into the array
  int * global_schur_root = NULL;
  int nglobal_schur_root = 0;
  if (rank == root){
    for ( int i = 0; i < size; i++ ){
      global_schur_ptr[i] = nglobal_schur_root;
      nglobal_schur_root += global_schur_count[i];
    }
    global_schur_root = new int[ nglobal_schur_root ];
  }

  // Gather all the global Schur variables to the root process
  MPI_Gatherv(global_schur_vars, nglobal_schur, MPI_INT,
              global_schur_root, global_schur_count, global_schur_ptr,
              MPI_INT, root, comm);

  // Duplicate the global list and uniquify the result
  int * unique_schur = NULL;
  int nunique_schur = 0;

  if (rank == root){
    unique_schur = new int[ nglobal_schur_root ];
    memcpy(unique_schur, global_schur_root, nglobal_schur_root*sizeof(int));
    nunique_schur = FElibrary::uniqueSort(unique_schur, nglobal_schur_root);

    // For each global Schur variable, now assign an output - the index
    // into the unique list of global Schur variables
    for ( int i = 0; i < nglobal_schur_root; i++ ){
      int * item = (int*)bsearch(&global_schur_root[i], unique_schur, 
                                 nunique_schur, sizeof(int),
                                 FElibrary::comparator);
      global_schur_root[i] = item - unique_schur; 
    }
  }
  
  // Broadcast the global size of the unique Schur complement matrix
  MPI_Bcast(&nunique_schur, 1, MPI_INT, root, comm);

  // Now, pass the unique Schur indices back to the original procs
  // This is an initial ordering of the variables, 
  MPI_Scatterv(global_schur_root, global_schur_count, global_schur_ptr, MPI_INT,
               new_global_schur_vars, nglobal_schur, MPI_INT, root, comm);

  int nremain = nunique_schur % size;
  int nlocal_schur = (int)((nunique_schur - nremain)/size);
  
  if (rank < nremain){
    nlocal_schur += 1;  
  }

  // unique_schur are the tacs variables corresponding to
  // the local ordering, stored on the root processor. Gather the 
  // number of expected variables on each proc to the root.
  // Reuse the global_schur_count and global_schur_ptr from above
  MPI_Gather(&nlocal_schur, 1, MPI_INT, 
             global_schur_count, 1, MPI_INT, root, comm);

  if (rank == root){
    nglobal_schur_root = 0;
    for ( int i = 0; i < size; i++ ){
      global_schur_ptr[i] = nglobal_schur_root;
      nglobal_schur_root += global_schur_count[i];
    }
  }

  // Send unique_schur back to the owning processes
  int * tacs_schur_vars = new int[nlocal_schur];
  MPI_Scatterv(unique_schur, global_schur_count, global_schur_ptr, MPI_INT,
               tacs_schur_vars, nlocal_schur, MPI_INT, root, comm);

  // Free memory not required anymore
  if (rank == root){
    delete [] global_schur_count;
    delete [] global_schur_ptr;
    delete [] global_schur_root;
    delete [] unique_schur;
  }

  // Set up information required for the global Schur complement matrix
  // Set the variable map
  schur_map = new VarMap(comm, nlocal_schur);
  schur_map->incref();

  // Create the index set for the new Schur complement variables
  BVecIndices * schur_index =
    new BVecIndices(&new_global_schur_vars, nglobal_schur);
  
  // Create the Schur complement matrix
  schur_dist = new BVecDistribute(schur_map, schur_index);
  schur_dist->incref();

  BVecIndices * tacs_schur_index = 
    new BVecIndices(&tacs_schur_vars, nlocal_schur);

  // Create the index set for the global Schur complement variables
  tacs_schur_dist = new BVecDistribute(smat->getVarMap(), 
                                       tacs_schur_index);
  tacs_schur_dist->incref();

  // Retrieve the non-zero pattern from the local Schur complement
  int bs, nschur_vars;
  const int *rowp, *cols;
  TacsScalar * vals;
  Sc->getArrays(&bs, &nschur_vars, &nschur_vars,
                &rowp, &cols, &vals);

  // Determine the size of the global Schur complement
  int M = nunique_schur;
  int N = nunique_schur;
  schur_index->getIndices(&new_global_schur_vars);

  // Determine the number of blocks to use per block-cylic block
  int bsize = B->getBlockSize();
  int csr_blocks_per_block = 36/bsize;

  // Create the global block-cyclic Schur complement matrix
  pdmat = new PDMat(comm, M, N, bsize, new_global_schur_vars, 
                    nschur_vars, rowp, cols, csr_blocks_per_block,
                    reorder_schur_complement);
  pdmat->incref();

  // Allocate space for local storage of vectors
  int xsize = bsize*b_map->getDim();
  int ysize = bsize*c_map->getDim();
  xlocal = new TacsScalar[ xsize ];
  yinterface = new TacsScalar[ ysize ];

  memset(xlocal, 0, xsize*sizeof(TacsScalar));
  memset(yinterface, 0, ysize*sizeof(TacsScalar));

  // Allocate the Schur complement vectors 
  yschur = new BVec(schur_map, bsize);
  gschur = new BVec(schur_map, bsize);
  yschur->incref();
  gschur->incref();
}

/*
  Destructor for the PcScMat preconditioner object
*/
PcScMat::~PcScMat(){
  // Decrease reference counts to the matrices
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
void PcScMat::testSchurComplement( TACSVec * tin, TACSVec * tout ){
  BVec *in, *out;
  in = dynamic_cast<BVec*>(tin);
  out = dynamic_cast<BVec*>(tout);

  if (in && out){
    // Test two methods of computing the effect of the Schur complement
    out->zeroEntries();
    
    int bsize = B->getBlockSize();
    int c_size = bsize*c_map->getDim();
    TacsScalar * temp = new TacsScalar[c_size];
    c_map->beginForward(in, yinterface);
    c_map->endForward(in, yinterface);
    
    Sc->mult(yinterface, temp);
    
    c_map->beginReverse(temp, out, BVecDistribute::ADD);
    c_map->endReverse(temp, out, BVecDistribute::ADD);
    
    delete [] temp;
    
    schur_dist->beginReverse(yinterface, yschur, BVecDistribute::INSERT);
    schur_dist->endReverse(yinterface, yschur, BVecDistribute::INSERT);
    
    TacsScalar *y, *g;
    int y_size = yschur->getArray(&y);
    gschur->getArray(&g);
    pdmat->mult(y_size, y, g);
    
    TacsScalar outnorm = out->norm();
    TacsScalar gnorm = gschur->norm();
    
    // Now collect the elements directly using the tacs_schur_dist object
    yschur->getArray(&y);
    tacs_schur_dist->beginForward(in, y);
    tacs_schur_dist->endForward(in, y);
    
    pdmat->mult(y_size, y, g);
    TacsScalar gnorm2 = gschur->norm(); 
    
    int rank;
    MPI_Comm_rank(c_map->getMPIComm(), &rank);
    
    if (rank == 0){
      printf("Schur complement consistency test: \n");
      printf("|Full matrix|     = %25.15e \n", RealPart(outnorm));
      printf("|Local to Schur|  = %25.15e \n", RealPart(gnorm));
      printf("|Global to Schur| = %25.15e \n", RealPart(gnorm2));
    }
  }
  else {
    fprintf(stderr, "PcScMat type error: Input/output must be BVec\n");
  }
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
  TacsScalar * scvals;
  Sc->getArrays(&bsize, &mlocal, &nlocal,
		&rowp, &cols, &scvals);

  // Retrieve the indices for the local Schur complement
  BVecIndices * schur_index = schur_dist->getBVecIndices();
  int * sc_vars;
  schur_index->getIndices(&sc_vars);

  // Add the values into the global Schur complement matrix
  // using either the alltoall approach or a sequential add values 
  // approach that uses less memory
  if (use_pdmat_alltoall){
    pdmat->addAlltoallValues(bsize, mlocal, sc_vars,
			     rowp, cols, scvals);
  }
  else {
    pdmat->addAllValues(bsize, mlocal, sc_vars,
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

  x = B^{-1} (f - E y)
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
void PcScMat::applyFactor( TACSVec * tin, TACSVec * tout ){
  // First, perform a safe down-cast from TACSVec to BVec
  BVec *in, *out;
  in = dynamic_cast<BVec*>(tin);
  out = dynamic_cast<BVec*>(tout);

  if (in && out){
    // Set the variables for the back-solve monitor
    double local_time = 0.0;
    double schur_time = 0.0;
    
    if (monitor_back_solve){
      local_time = -MPI_Wtime();
    }
    
    // Pass g to the global Schur complement
    TacsScalar * g;
    int gschur_size = gschur->getArray(&g);
    tacs_schur_dist->beginForward(in, g);

    // Re-order the local variables into xlocal
    b_map->beginForward(in, xlocal);
    b_map->endForward(in, xlocal);

    // xlocal = L^{-1} f 
    Bpc->applyLower(xlocal, xlocal);

    // yinterface = F U^{-1} xlocal = F U^{-1} L^{-1} f
    Fpc->mult(xlocal, yinterface);

    // Transmit yinterface to the Schur complement system
    // Pass F U^{-1} L^{-1} f to yschur
    yschur->zeroEntries();
    schur_dist->beginReverse(yinterface, yschur, BVecDistribute::ADD);

    // Finish accepting g from the global input vector
    tacs_schur_dist->endForward(in, g);
    schur_dist->endReverse(yinterface, yschur, BVecDistribute::ADD); 
    
    // Compute the right hand side: g - F U^{-1} L^{-1} f
    gschur->axpy(-1.0, yschur);

    if (monitor_back_solve){
      schur_time = -MPI_Wtime();
    }

    // Apply the global Schur complement factorization to the right
    // hand side 
    pdmat->applyFactor(gschur_size, g);

    if (monitor_back_solve){
      schur_time += MPI_Wtime();
      local_time -= schur_time;
    }
    
    // copy the values from gschur to yschur
    yschur->copyValues(gschur);

    // Pass yschur solution back to the local variables
    schur_dist->beginForward(yschur, yinterface);
    
    // Pass yschur also to the global variables
    TacsScalar * y;
    yschur->getArray(&y);
    tacs_schur_dist->beginReverse(y, out, BVecDistribute::INSERT);
    
    // Compute yinterface = yinterface - L^{-1} E y
    // Note: scale y, by -1 first 
    schur_dist->endForward(yschur, yinterface);
    int one = 1;
    int len = Bpc->getBlockSize()*c_map->getDim();
    TacsScalar alpha = -1.0; 
    BLASscal(&len, &alpha, yinterface, &one);
    
    // Compute xlocal = xlocal - L^{-1} E * yinterface
    Epc->multAdd(yinterface, xlocal, xlocal);
    
    // Compute xlocal = U^{-1} xlocal
    Bpc->applyUpper(xlocal, xlocal);
    
    b_map->beginReverse(xlocal, out, BVecDistribute::INSERT);
    tacs_schur_dist->endReverse(y, out, BVecDistribute::INSERT);
    b_map->endReverse(xlocal, out, BVecDistribute::INSERT);    
    
    if (monitor_back_solve){
      local_time += MPI_Wtime();
      
      int rank;
      MPI_Comm_rank(b_map->getMPIComm(), &rank);
      printf("[%d] Local back-solve time:  %8.4f\n", rank, local_time);
      printf("[%d] Schur back-solve time:  %8.4f\n", rank, schur_time);
    }
  }
  else {
    fprintf(stderr, "PcScMat type error: Input/output must be BVec\n");
  }
}
