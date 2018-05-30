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

#include "FEMat.h"
#include "FElibrary.h"

/*
  Finite-element matrix implementation
*/

/*!
  This constructor produces the data required for the ScMat matrix
  from the data provided by TACS.

  The matrix is stored using the following block structure:

  [ B, E ]
  [ F, C ]

  The blocks are partitioned such that the variables requried for the
  B-matrix are only local - they are not referenced on other
  processors. All coupling between processors occurs in the matrix
  C. No matrix assembly is required to perform matrix-vector products
  because the effect of assembly is handled by adding contributions
  from the C components.

  This matrix is especially well suited for almost dense
  preconditioning. In this case, the factorization can be split into
  two stages. First, a factorization on each process to determine the
  local Schur complement contributions: S = (C - F B^{-1} E). Then, an
  assembly of the Schur complement contributions from each processor
  must be done, in parallel. Followed by some factorization process
  for the global Schur complement that represents the interface
  problem.

  Input:

  rmap: The variable map that defines the row-distribution of the
  global matrix. This is required for matrix-vector products.

  nlocal_vars, rowp, cols: The CSR non-zero structure of all local
  variables

  b_local_indices: The local indicies of the B-matrix
  b_map: The map from the global variables to the local indices

  c_local_indices: The local indices of the C-matrix
  c_map: The map from the global

  Note that b_local_indices/c_local_indices are not required after
  initialization of the FEMat class.

  Procedure:

  The non-zero pattern of the matrices B, E, F and C are
  determined. The index sets b_local_indices and c_local_indices are
  defined such that:

  int *bi, *ci;
  b_local_indices->getIndices(&bi);
  c_local_indices->getIndices(&ci);

  The diagonal blocks are formed such that:
  B[i,j] = A[bi[i], bi[j]]
  C[i,j] = A[ci[i], ci[j]]

  The off-diagonal components are formed such that:
  E[i,j] = A[bi[i], ci[j]]
  F[i,j] = A[ci[i], bi[j]]

  However, the matrix A is provided in a CSR format. This complicates
  the code somewhat.
*/
FEMat::FEMat( TACSThreadInfo *thread_info, TACSVarMap *_rmap,
              int bsize, int nlocal_vars,
              const int *rowp, const int *cols,
              TACSBVecIndices *b_local_indices, TACSBVecDistribute *_b_map,
              TACSBVecIndices *c_local_indices, TACSBVecDistribute *_c_map ){
  // Get the block size
  int rank;
  MPI_Comm_rank(_rmap->getMPIComm(), &rank);

  // bi : i in B -> bi[i] in A
  const int *bi, *ci;
  b_local_indices->getIndices(&bi);
  c_local_indices->getIndices(&ci);
  int Nb = b_local_indices->getNumIndices();
  int Nc = c_local_indices->getNumIndices();

  if (nlocal_vars != Nb + Nc){
    fprintf(stderr, "[%d] FEMat error, CSR data is incorrect size\n",
            rank);
  }

  // Compute the inverse mapping
  // binv : i in A -> binv[i] in B
  int *binv = new int[ nlocal_vars ];
  int *cinv = new int[ nlocal_vars ];

  for ( int i = 0; i < nlocal_vars; i++ ){
    binv[i] = cinv[i] = -1;
  }
  for ( int i = 0; i < Nb; i++ ){
    binv[bi[i]] = i;
  }
  for ( int i = 0; i < Nc; i++ ){
    cinv[ci[i]] = i;
  }

  // Count up the size of the b matrix
  int *browp = new int[ Nb+1 ];
  int *erowp = new int[ Nb+1 ];
  int *frowp = new int[ Nc+1 ];
  int *crowp = new int[ Nc+1 ];
  memset(browp, 0, (Nb+1)*sizeof(int));
  memset(erowp, 0, (Nb+1)*sizeof(int));
  memset(frowp, 0, (Nc+1)*sizeof(int));
  memset(crowp, 0, (Nc+1)*sizeof(int));

  // Count up the size of the different matrices
  for ( int i = 0; i < Nb; i++ ){
    int row = bi[i];

    // Add the variables in the row to either B or E
    for ( int jp = rowp[row]; jp < rowp[row+1]; jp++ ){
      if (binv[cols[jp]] >= 0){
        browp[i+1]++;
      }
      else if (cinv[cols[jp]] >= 0){
        erowp[i+1]++;
      }
      else {
        fprintf(stderr,
                "[%d] FEMat error, C/B indices not complete\n",
                rank);
      }
    }
  }

  for ( int i = 0; i < Nc; i++ ){
    int row = ci[i];

    // Add the variables in the row to either B or E
    for ( int jp = rowp[row]; jp < rowp[row+1]; jp++ ){
      if (binv[cols[jp]] >= 0){
        frowp[i+1]++;
      }
      else if (cinv[cols[jp]] >= 0){
        crowp[i+1]++;
      }
      else {
        fprintf(stderr,
                "[%d] FEMat error, C/B indices not complete\n",
                rank);
      }
    }
  }

  // Now, add up all the indices
  for ( int i = 0; i < Nb; i++ ){
    browp[i+1] = browp[i+1] + browp[i];
    erowp[i+1] = erowp[i+1] + erowp[i];
  }

  for ( int i = 0; i < Nc; i++ ){
    frowp[i+1] = frowp[i+1] + frowp[i];
    crowp[i+1] = crowp[i+1] + crowp[i];
  }

  // Now, prepare to add in all the indices.
  // This modifies the pointers *rowp, these
  // will be adjusted back after the computation
  int *bcols = new int[ browp[Nb] ];
  int *ecols = new int[ erowp[Nb] ];
  int *fcols = new int[ frowp[Nc] ];
  int *ccols = new int[ crowp[Nc] ];

  // Count up the size of the different matrices
  for ( int i = 0; i < Nb; i++ ){
    int row = bi[i];

    // Add the variables in the row to either B or E
    for ( int jp = rowp[row]; jp < rowp[row+1]; jp++ ){
      if (binv[cols[jp]] >= 0){
        bcols[browp[i]] = binv[cols[jp]];
        browp[i]++;
      }
      else if (cinv[cols[jp]] >= 0){
        ecols[erowp[i]] = cinv[cols[jp]];
        erowp[i]++;
      }
      else {
        fprintf(stderr,
                "[%d] FEMat error, C/B indices not complete\n",
                rank);
      }
    }
  }

  for ( int i = 0; i < Nc; i++ ){
    int row = ci[i];

    // Add the variables in the row to either B or E
    for ( int jp = rowp[row]; jp < rowp[row+1]; jp++ ){
      if (binv[cols[jp]] >= 0){
        fcols[frowp[i]] = binv[cols[jp]];
        frowp[i]++;
      }
      else if (cinv[cols[jp]] >= 0){
        ccols[crowp[i]] = cinv[cols[jp]];
        crowp[i]++;
      }
      else {
        fprintf(stderr,
                "[%d] FEMat error, C/B indices not complete\n",
                rank);
      }
    }
  }

  delete [] binv;
  delete [] cinv;

  // Adjust the pointers to the correct values
  for ( int i = 0, bnext = 0, enext = 0; i <= Nb; i++ ){
    int tb = browp[i];
    browp[i] = bnext;
    bnext = tb;

    int te = erowp[i];
    erowp[i] = enext;
    enext = te;
  }

  for ( int i = 0, fnext = 0, cnext = 0; i <= Nc; i++ ){
    int tf = frowp[i];
    frowp[i] = fnext;
    fnext = tf;

    int tc = crowp[i];
    crowp[i] = cnext;
    cnext = tc;
  }

  // Sort the rows
  for ( int i = 0; i < Nb; i++ ){
    int nb = browp[i+1] - browp[i];
    if (nb != FElibrary::uniqueSort(&bcols[browp[i]], nb)){
      fprintf(stderr, "FEMat error, B input nz-pattern not unique\n");
    }

    int ne = erowp[i+1] - erowp[i];
    if (ne != FElibrary::uniqueSort(&ecols[erowp[i]], ne)){
      fprintf(stderr, "FEMat error, E input nz-pattern not unique\n");
    }
  }

  for ( int i = 0; i < Nc; i++ ){
    int nf = frowp[i+1] - frowp[i];
    if (nf != FElibrary::uniqueSort(&fcols[frowp[i]], nf) ){
      fprintf(stderr, "FEMat error, F input nz-pattern not unique\n");
    }

    int nc = crowp[i+1] - crowp[i];
    if (nc != FElibrary::uniqueSort(&ccols[crowp[i]], nc)){
      fprintf(stderr, "FEMat error, C input nz-pattern not unique\n");
    }
  }

  // Create the block matrices
  MPI_Comm comm = _rmap->getMPIComm();
  BCSRMat *_B = new BCSRMat(comm, thread_info, bsize, Nb, Nb, &browp, &bcols);
  BCSRMat *_E = new BCSRMat(comm, thread_info, bsize, Nb, Nc, &erowp, &ecols);
  BCSRMat *_F = new BCSRMat(comm, thread_info, bsize, Nc, Nb, &frowp, &fcols);
  BCSRMat *_C = new BCSRMat(comm, thread_info, bsize, Nc, Nc, &crowp, &ccols);

  // Insure that the inverse look-up has been allocated
  (_b_map->getIndices())->setUpInverse();
  (_c_map->getIndices())->setUpInverse();

  // Initialize the underlying class
  init(_rmap, _B, _E, _F, _C, _b_map, _c_map);
}

FEMat::~FEMat(){}

/*!
  Add the values into the appropriate block matrix.

  This performs a few steps. First, determine the indices for the
  appropriate matrix. If no entries are found for the given non-zero
  entries, print an error.

  Temporary arrays are allocated for the indices. If the number of
  rows/columns in the dense sub-matrix are less than 128, then use the
  statically allocated array, otherwise dynamically allocate temporary
  arrays to store the indices.
*/
void FEMat::addValues( int nrow, const int *row,
                       int ncol, const int *col,
                       int nv, int mv, const TacsScalar *values ){
  int bsize = B->getBlockSize();
  TACSBVecIndices *bindx = b_map->getIndices();
  TACSBVecIndices *cindx = c_map->getIndices();

  // Set up storage for the values of the local variable numbers
  int array[256];
  int *temp = NULL, *bcols = NULL, *ccols = NULL;

  if (ncol > 128){
    temp = new int[ 2*ncol ];
    bcols = &temp[0];
    ccols = &temp[ncol];
  }
  else {
    bcols = &array[0];
    ccols = &array[ncol];
  }

  // Flag to indicate whether we found any items in the c-index set
  int cflag = 0;

  // Convert the columns first so that when we loop over the
  // rows, we can add the entire row in a single shot
  for ( int i = 0; i < ncol; i++ ){
    int c = col[i];
    bcols[i] = -1;
    ccols[i] = -1;

    if (c >= 0){
      // First look for the column in the c-index set
      bcols[i] = bindx->findIndex(c);

      // If it wasn't there, search in the b-index set
      if (bcols[i] < 0){
        cflag = 1;
        ccols[i] = cindx->findIndex(c);
      }

      // If neither turned up anything, print an error
      if ((ccols[i] < 0) &&
          (bcols[i] < 0)){
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global column variable %d not found\n",
                rank, c);
      }
    }
  }

  // Now, loop over the rows and add each one to the corresponding
  // B and F matrices and possibly F or C matrix if cflag is true
  for ( int i = 0; i < nrow; i++ ){
    // Check if the row is on this processor
    int r = row[i];

    if (r >= 0){
      int br = 0, cr = 0;
      if ((br = bindx->findIndex(r)) >= 0){
        B->addRowValues(br, ncol, bcols, mv, &values[mv*i*bsize]);
        if (cflag){
          E->addRowValues(br, ncol, ccols, mv, &values[mv*i*bsize]);
        }
      }
      else if ((cr = cindx->findIndex(r)) >= 0){
        F->addRowValues(cr, ncol, bcols, mv, &values[mv*i*bsize]);
        if (cflag){
          C->addRowValues(cr, ncol, ccols, mv, &values[mv*i*bsize]);
        }
      }
      else {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global row variable %d not found\n",
                rank, r);
      }
    }
  }

  if (temp){ delete [] temp; }
}

/*
  Add a weighted sum of the dense input matrix.

  This function adds an inner product of a weighting matrix with a
  dense matrix to the FEmat matrix. The weight matrix is a sparse,
  low-dimensional matrix given in a CSR-type format. The code takes
  this representation of W and adds the terms:

  self <- self + W^{T}*Kv*W

  to the values in the matrix. This code can be used to add the effect
  of dependent nodes to the matrix.

  input:
  nvars:    the number of block variables in the dense input matrix
  varp:     pointer to the start of each weight row; len(varp) = nvars+1
  vars:     variable numbers in the matrix; len(vars) = varp[nvars]
  weights:  weighting values for each variable; len(weights) = len(vars)
  nv, nw:   the number of rows and number of columns in the values matrix
  values:   the dense input matrix
*/
void FEMat::addWeightValues( int nvars, const int *varp, const int *vars,
                             const TacsScalar *weights,
                             int nv, int mv, const TacsScalar *values,
                             MatrixOrientation matOr ){
  // The block size for the matrix
  int bsize = B->getBlockSize();

  // Get the index sets corresponding to the B/C matrices
  TACSBVecIndices *bindx = b_map->getIndices();
  TACSBVecIndices *cindx = c_map->getIndices();

  // The number of variables we'll have to convert = row dim(W^{T})
  int n = varp[nvars];

  // Set the default array for the local indices - the array
  // will be expanded when it is not large enough
  int array[256];
  int *temp = NULL, *bvars = NULL, *cvars = NULL;

  if (n > 128){
    temp = new int[ 2*n ];
    bvars = &temp[0];
    cvars = &temp[n];
  }
  else {
    bvars = &array[0];
    cvars = &array[n];
  }

  // Flag to indicate whether we found any items in the c-index set
  int cflag = 0;

  // Convert the columns first so that when we loop over the
  // rows, we can add the entire row in a single shot
  for ( int i = 0; i < n; i++ ){
    int c = vars[i];
    bvars[i] = -1;
    cvars[i] = -1;

    if (c >= 0){
      // First look for the column in the c-index set
      bvars[i] = bindx->findIndex(c);

      // If it wasn't there, search in the b-index set
      if (bvars[i] < 0){
        cflag = 1;
        cvars[i] = cindx->findIndex(c);
      }

      // If neither turned up anything, print an error
      if ((cvars[i] < 0) &&
          (bvars[i] < 0)){
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global column variable %d not found\n",
                rank, c);
      }
    }
  }

  // Set the increment along the row or column of the matrix depending
  // on whether we are adding the original matrix or its transpose
  int incr = mv;
  if (matOr == TRANSPOSE){
    incr = 1;
  }

  // Now, loop over the rows and add each one to the corresponding
  // B nd F matrices and possibly F or C matrix if cflag is true
  for ( int i = 0; i < nvars; i++ ){
    for ( int j = varp[i]; j < varp[i+1]; j++ ){
      // Check where to place this value
      if (bvars[j] >= 0){
        B->addRowWeightValues(weights[j], bvars[j],
                              nvars, varp, bvars, weights,
                              mv, &values[incr*i*bsize], matOr);
        if (cflag){
          E->addRowWeightValues(weights[j], bvars[j],
                                nvars, varp, cvars, weights,
                                mv, &values[incr*i*bsize], matOr);
        }
      }
      else if (cvars[j] >= 0){
        F->addRowWeightValues(weights[j], cvars[j],
                              nvars, varp, bvars, weights,
                              mv, &values[incr*i*bsize], matOr);
        if (cflag){
          C->addRowWeightValues(weights[j], cvars[j],
                                nvars, varp, cvars, weights,
                                mv, &values[incr*i*bsize], matOr);
        }
      }
      else {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global row variable %d not found\n",
                rank, vars[j]);
      }
    }
  }

  if (temp){ delete [] temp; }

}

/*
  Create a vector that is compatible with this matrix
*/
TACSVec *FEMat::createVec(){
  return new TACSBVec(rmap, B->getBlockSize());
}

/*
  Apply the Dirichlet boundary conditions by zeroing the rows
  associated with the boundary conditions and setting corresponding
  the diagonal entries to 1.
*/
void FEMat::applyBCs( TACSBcMap *bcmap ){
  if (bcmap){
    TACSBVecIndices *bindx = b_map->getIndices();
    TACSBVecIndices *cindx = c_map->getIndices();

    int mpiRank;
    MPI_Comm_rank(rmap->getMPIComm(), &mpiRank);

    const int *ownerRange;
    rmap->getOwnerRange(&ownerRange);

    // apply the boundary conditions
    const int *nodes, *vars;
    const TacsScalar *values;
    int nbcs = bcmap->getBCs(&nodes, &vars, &values);

    // Get the matrix values
    for ( int i = 0; i < nbcs; i++ ){
      // Find block i and zero out the variables associated with it
      int br = bindx->findIndex(nodes[i]), cr = 0;
      if (br >= 0){
        int ident = 1; // Replace the diagonal with the identity matrix
        B->zeroRow(br, vars[i], ident);
        ident = 0;
        E->zeroRow(br, vars[i], ident);
      }
      else if ((cr = cindx->findIndex(nodes[i])) >= 0){
        int cident = 0;
        int fident = 0;
        if (nodes[i] >= ownerRange[mpiRank] &&
            nodes[i] < ownerRange[mpiRank+1]){
          // Only use the identity here if it is locally owned
          cident = 1;
        }
        F->zeroRow(cr, vars[i], fident);
        C->zeroRow(cr, vars[i], cident);
      }
    }
  }
}
