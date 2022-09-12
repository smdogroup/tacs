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

#include "TACSSchurMat.h"

#include "TacsUtilities.h"
#include "tacslapack.h"

/*
  Schur-complement based matrix
*/

/**
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
TACSSchurMat::TACSSchurMat(TACSNodeMap *_rmap, BCSRMat *_B, BCSRMat *_E,
                           BCSRMat *_F, BCSRMat *_C, TACSBVecDistribute *_b_map,
                           TACSBVecDistribute *_c_map) {
  init(_rmap, _B, _E, _F, _C, _b_map, _c_map);
}

/**
  This constructor produces the data required for the TACSSchurMat
  object from data provided by TACS.

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
TACSSchurMat::TACSSchurMat(TACSThreadInfo *thread_info, TACSNodeMap *_rmap,
                           int bsize, int nlocal_vars, const int *rowp,
                           const int *cols, TACSBVecIndices *b_local_indices,
                           TACSBVecDistribute *_b_map,
                           TACSBVecIndices *c_local_indices,
                           TACSBVecDistribute *_c_map) {
  // Get the block size
  int rank;
  MPI_Comm_rank(_rmap->getMPIComm(), &rank);

  // bi : i in B -> bi[i] in A
  const int *bi, *ci;
  b_local_indices->getIndices(&bi);
  c_local_indices->getIndices(&ci);
  int Nb = b_local_indices->getNumIndices();
  int Nc = c_local_indices->getNumIndices();

  if (nlocal_vars != Nb + Nc) {
    fprintf(stderr, "[%d] FEMat error, CSR data is incorrect size\n", rank);
  }

  // Compute the inverse mapping
  // binv : i in A -> binv[i] in B
  int *binv = new int[nlocal_vars];
  int *cinv = new int[nlocal_vars];

  for (int i = 0; i < nlocal_vars; i++) {
    binv[i] = cinv[i] = -1;
  }
  for (int i = 0; i < Nb; i++) {
    binv[bi[i]] = i;
  }
  for (int i = 0; i < Nc; i++) {
    cinv[ci[i]] = i;
  }

  // Count up the size of the b matrix
  int *browp = new int[Nb + 1];
  int *erowp = new int[Nb + 1];
  int *frowp = new int[Nc + 1];
  int *crowp = new int[Nc + 1];
  memset(browp, 0, (Nb + 1) * sizeof(int));
  memset(erowp, 0, (Nb + 1) * sizeof(int));
  memset(frowp, 0, (Nc + 1) * sizeof(int));
  memset(crowp, 0, (Nc + 1) * sizeof(int));

  // Count up the size of the different matrices
  for (int i = 0; i < Nb; i++) {
    int row = bi[i];

    // Add the variables in the row to either B or E
    for (int jp = rowp[row]; jp < rowp[row + 1]; jp++) {
      if (binv[cols[jp]] >= 0) {
        browp[i + 1]++;
      } else if (cinv[cols[jp]] >= 0) {
        erowp[i + 1]++;
      } else {
        fprintf(stderr, "[%d] FEMat error, C/B indices not complete\n", rank);
      }
    }
  }

  for (int i = 0; i < Nc; i++) {
    int row = ci[i];

    // Add the variables in the row to either B or E
    for (int jp = rowp[row]; jp < rowp[row + 1]; jp++) {
      if (binv[cols[jp]] >= 0) {
        frowp[i + 1]++;
      } else if (cinv[cols[jp]] >= 0) {
        crowp[i + 1]++;
      } else {
        fprintf(stderr, "[%d] FEMat error, C/B indices not complete\n", rank);
      }
    }
  }

  // Now, add up all the indices
  for (int i = 0; i < Nb; i++) {
    browp[i + 1] = browp[i + 1] + browp[i];
    erowp[i + 1] = erowp[i + 1] + erowp[i];
  }

  for (int i = 0; i < Nc; i++) {
    frowp[i + 1] = frowp[i + 1] + frowp[i];
    crowp[i + 1] = crowp[i + 1] + crowp[i];
  }

  // Now, prepare to add in all the indices.
  // This modifies the pointers *rowp, these
  // will be adjusted back after the computation
  int *bcols = new int[browp[Nb]];
  int *ecols = new int[erowp[Nb]];
  int *fcols = new int[frowp[Nc]];
  int *ccols = new int[crowp[Nc]];

  // Count up the size of the different matrices
  for (int i = 0; i < Nb; i++) {
    int row = bi[i];

    // Add the variables in the row to either B or E
    for (int jp = rowp[row]; jp < rowp[row + 1]; jp++) {
      if (binv[cols[jp]] >= 0) {
        bcols[browp[i]] = binv[cols[jp]];
        browp[i]++;
      } else if (cinv[cols[jp]] >= 0) {
        ecols[erowp[i]] = cinv[cols[jp]];
        erowp[i]++;
      } else {
        fprintf(stderr, "[%d] FEMat error, C/B indices not complete\n", rank);
      }
    }
  }

  for (int i = 0; i < Nc; i++) {
    int row = ci[i];

    // Add the variables in the row to either B or E
    for (int jp = rowp[row]; jp < rowp[row + 1]; jp++) {
      if (binv[cols[jp]] >= 0) {
        fcols[frowp[i]] = binv[cols[jp]];
        frowp[i]++;
      } else if (cinv[cols[jp]] >= 0) {
        ccols[crowp[i]] = cinv[cols[jp]];
        crowp[i]++;
      } else {
        fprintf(stderr, "[%d] FEMat error, C/B indices not complete\n", rank);
      }
    }
  }

  delete[] binv;
  delete[] cinv;

  // Adjust the pointers to the correct values
  for (int i = 0, bnext = 0, enext = 0; i <= Nb; i++) {
    int tb = browp[i];
    browp[i] = bnext;
    bnext = tb;

    int te = erowp[i];
    erowp[i] = enext;
    enext = te;
  }

  for (int i = 0, fnext = 0, cnext = 0; i <= Nc; i++) {
    int tf = frowp[i];
    frowp[i] = fnext;
    fnext = tf;

    int tc = crowp[i];
    crowp[i] = cnext;
    cnext = tc;
  }

  // Sort the rows
  for (int i = 0; i < Nb; i++) {
    int nb = browp[i + 1] - browp[i];
    if (nb != TacsUniqueSort(nb, &bcols[browp[i]])) {
      fprintf(stderr, "FEMat error, B input nz-pattern not unique\n");
    }

    int ne = erowp[i + 1] - erowp[i];
    if (ne != TacsUniqueSort(ne, &ecols[erowp[i]])) {
      fprintf(stderr, "FEMat error, E input nz-pattern not unique\n");
    }
  }

  for (int i = 0; i < Nc; i++) {
    int nf = frowp[i + 1] - frowp[i];
    if (nf != TacsUniqueSort(nf, &fcols[frowp[i]])) {
      fprintf(stderr, "FEMat error, F input nz-pattern not unique\n");
    }

    int nc = crowp[i + 1] - crowp[i];
    if (nc != TacsUniqueSort(nc, &ccols[crowp[i]])) {
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

/*
  Initialize the data for the TACSSchurMat class

  input:
  rmap:  the variable map for all variables
  B:     the diagonal block BCSR matrices
  E, F:  the off-diagonal block BCSR matrices
  C:     the diagonal block corresponding to the local Schur complement
  b_map: the global variables corresponding to the B-variables
  c_map: the global variables corresponding to the C-variables
*/
void TACSSchurMat::init(TACSNodeMap *_rmap, BCSRMat *_B, BCSRMat *_E,
                        BCSRMat *_F, BCSRMat *_C, TACSBVecDistribute *_b_map,
                        TACSBVecDistribute *_c_map) {
  rmap = _rmap;
  rmap->incref();

  B = _B;
  B->incref();
  E = _E;
  E->incref();
  F = _F;
  F->incref();
  C = _C;
  C->incref();

  b_map = _b_map;
  b_map->incref();
  c_map = _c_map;
  c_map->incref();

  // Check that the block dimensions work out
  int bs = B->getBlockSize();
  if (bs != E->getBlockSize() || bs != F->getBlockSize() ||
      bs != C->getBlockSize()) {
    fprintf(stderr, "TACSSchurMat error: block sizes do not match\n");
    return;
  }

  // Check that things are the correct dimensions
  if (B->getRowDim() != E->getRowDim()) {
    fprintf(stderr, "TACSSchurMat error: B, E row dimensions do not match\n");
    return;
  }
  if (F->getRowDim() != C->getRowDim()) {
    fprintf(stderr, "TACSSchurMat error: F, C row dimensions do not match\n");
    return;
  }

  if (B->getColDim() != F->getColDim()) {
    fprintf(stderr,
            "TACSSchurMat error: B, F column dimensions do not match\n");
    return;
  }
  if (E->getColDim() != C->getColDim()) {
    fprintf(stderr,
            "TACSSchurMat error: E, C column dimensions do not match\n");
    return;
  }

  if (B->getColDim() != b_map->getNumNodes()) {
    fprintf(stderr,
            "TACSSchurMat error: b_map dimensions do not "
            "match dimensions of B\n");
    return;
  }

  if (C->getColDim() != c_map->getNumNodes()) {
    fprintf(stderr,
            "TACSSchurMat error: c_map dimensions do not "
            "match dimensions of C\n");
    return;
  }

  // Allocate the memory for the preconditioning operations
  local_size = bs * (b_map->getNumNodes() + c_map->getNumNodes());
  local_offset = bs * b_map->getNumNodes();

  b_ctx = b_map->createCtx(bs);
  b_ctx->incref();
  c_ctx = c_map->createCtx(bs);
  c_ctx->incref();

  xlocal = new TacsScalar[local_size];
  ylocal = new TacsScalar[local_size];
  memset(xlocal, 0, local_size * sizeof(TacsScalar));
  memset(ylocal, 0, local_size * sizeof(TacsScalar));
}

/*
  The destructor for the TACSSchurMat
*/
TACSSchurMat::~TACSSchurMat() {
  if (rmap) {
    rmap->decref();
  }
  if (B) {
    B->decref();
  }
  if (E) {
    E->decref();
  }
  if (F) {
    F->decref();
  }
  if (C) {
    C->decref();
  }
  if (b_map) {
    b_map->decref();
  }
  if (c_map) {
    c_map->decref();
  }
  if (b_ctx) {
    b_ctx->decref();
  }
  if (c_ctx) {
    c_ctx->decref();
  }
  if (xlocal) {
    delete[] xlocal;
  }
  if (ylocal) {
    delete[] ylocal;
  }
}

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
void TACSSchurMat::addValues(int nrow, const int *row, int ncol, const int *col,
                             int nv, int mv, const TacsScalar *values) {
  int bsize = B->getBlockSize();
  TACSBVecIndices *bindx = b_map->getIndices();
  TACSBVecIndices *cindx = c_map->getIndices();

  // Set up storage for the values of the local variable numbers
  int array[256];
  int *temp = NULL, *bcols = NULL, *ccols = NULL;

  if (ncol > 128) {
    temp = new int[2 * ncol];
    bcols = &temp[0];
    ccols = &temp[ncol];
  } else {
    bcols = &array[0];
    ccols = &array[ncol];
  }

  // Flag to indicate whether we found any items in the c-index set
  int cflag = 0;

  // Convert the columns first so that when we loop over the
  // rows, we can add the entire row in a single shot
  for (int i = 0; i < ncol; i++) {
    int c = col[i];
    bcols[i] = -1;
    ccols[i] = -1;

    if (c >= 0) {
      // First look for the column in the c-index set
      bcols[i] = bindx->findIndex(c);

      // If it wasn't there, search in the b-index set
      if (bcols[i] < 0) {
        cflag = 1;
        ccols[i] = cindx->findIndex(c);
      }

      // If neither turned up anything, print an error
      if ((ccols[i] < 0) && (bcols[i] < 0)) {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global column variable %d not found\n", rank, c);
      }
    }
  }

  // Now, loop over the rows and add each one to the corresponding
  // B and F matrices and possibly F or C matrix if cflag is true
  for (int i = 0; i < nrow; i++) {
    // Check if the row is on this processor
    int r = row[i];

    if (r >= 0) {
      int br = 0, cr = 0;
      if ((br = bindx->findIndex(r)) >= 0) {
        B->addRowValues(br, ncol, bcols, mv, &values[mv * i * bsize]);
        if (cflag) {
          E->addRowValues(br, ncol, ccols, mv, &values[mv * i * bsize]);
        }
      } else if ((cr = cindx->findIndex(r)) >= 0) {
        F->addRowValues(cr, ncol, bcols, mv, &values[mv * i * bsize]);
        if (cflag) {
          C->addRowValues(cr, ncol, ccols, mv, &values[mv * i * bsize]);
        }
      } else {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global row variable %d not found\n", rank, r);
      }
    }
  }

  if (temp) {
    delete[] temp;
  }
}

/*
  Add a weighted sum of the dense input matrix.

  This function adds an inner product of a weighting matrix with a
  dense matrix to the  matrix. The weight matrix is a sparse,
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
void TACSSchurMat::addWeightValues(int nvars, const int *varp, const int *vars,
                                   const TacsScalar *weights, int nv, int mv,
                                   const TacsScalar *values,
                                   MatrixOrientation matOr) {
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

  if (n > 128) {
    temp = new int[2 * n];
    bvars = &temp[0];
    cvars = &temp[n];
  } else {
    bvars = &array[0];
    cvars = &array[n];
  }

  // Flag to indicate whether we found any items in the c-index set
  int cflag = 0;

  // Convert the columns first so that when we loop over the
  // rows, we can add the entire row in a single shot
  for (int i = 0; i < n; i++) {
    int c = vars[i];
    bvars[i] = -1;
    cvars[i] = -1;

    if (c >= 0) {
      // First look for the column in the c-index set
      bvars[i] = bindx->findIndex(c);

      // If it wasn't there, search in the b-index set
      if (bvars[i] < 0) {
        cflag = 1;
        cvars[i] = cindx->findIndex(c);
      }

      // If neither turned up anything, print an error
      if ((cvars[i] < 0) && (bvars[i] < 0)) {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global column variable %d not found\n", rank, c);
      }
    }
  }

  // Set the increment along the row or column of the matrix depending
  // on whether we are adding the original matrix or its transpose
  int incr = mv;
  if (matOr == TACS_MAT_TRANSPOSE) {
    incr = 1;
  }

  // Now, loop over the rows and add each one to the corresponding
  // B nd F matrices and possibly F or C matrix if cflag is true
  for (int i = 0; i < nvars; i++) {
    for (int j = varp[i]; j < varp[i + 1]; j++) {
      // Check where to place this value
      if (bvars[j] >= 0) {
        B->addRowWeightValues(weights[j], bvars[j], nvars, varp, bvars, weights,
                              mv, &values[incr * i * bsize], matOr);
        if (cflag) {
          E->addRowWeightValues(weights[j], bvars[j], nvars, varp, cvars,
                                weights, mv, &values[incr * i * bsize], matOr);
        }
      } else if (cvars[j] >= 0) {
        F->addRowWeightValues(weights[j], cvars[j], nvars, varp, bvars, weights,
                              mv, &values[incr * i * bsize], matOr);
        if (cflag) {
          C->addRowWeightValues(weights[j], cvars[j], nvars, varp, cvars,
                                weights, mv, &values[incr * i * bsize], matOr);
        }
      } else {
        int rank;
        MPI_Comm_rank(b_map->getMPIComm(), &rank);
        fprintf(stderr, "[%d] Global row variable %d not found\n", rank,
                vars[j]);
      }
    }
  }

  if (temp) {
    delete[] temp;
  }
}

/*
  Create a vector that is compatible with this matrix
*/
TACSVec *TACSSchurMat::createVec() {
  return new TACSBVec(rmap, B->getBlockSize());
}

/*
  Apply the Dirichlet boundary conditions by zeroing the rows
  associated with the boundary conditions and setting corresponding
  the diagonal entries to 1.
*/
void TACSSchurMat::applyBCs(TACSBcMap *bcmap) {
  if (bcmap) {
    TACSBVecIndices *bindx = b_map->getIndices();
    TACSBVecIndices *cindx = c_map->getIndices();

    int mpiRank;
    MPI_Comm_rank(rmap->getMPIComm(), &mpiRank);

    const int *ownerRange;
    rmap->getOwnerRange(&ownerRange);

    // apply the boundary conditions
    const int *nodes, *vars;
    int nbcs = bcmap->getBCs(&nodes, &vars, NULL);

    // Get the matrix values
    for (int i = 0; i < nbcs; i++) {
      // Find block i and zero out the variables associated with it
      int br = bindx->findIndex(nodes[i]), cr = 0;
      if (br >= 0) {
        int ident = 1;  // Replace the diagonal with the identity matrix
        B->zeroRow(br, vars[i], ident);
        ident = 0;
        E->zeroRow(br, vars[i], ident);
      } else if ((cr = cindx->findIndex(nodes[i])) >= 0) {
        int cident = 0;
        int fident = 0;
        if (nodes[i] >= ownerRange[mpiRank] &&
            nodes[i] < ownerRange[mpiRank + 1]) {
          // Only use the identity here if it is locally owned
          cident = 1;
        }
        F->zeroRow(cr, vars[i], fident);
        C->zeroRow(cr, vars[i], cident);
      }
    }
  }
}

/*!
  Apply the transpose of the Dirichlet boundary conditions (zero-ing)
  columns and setting the diagonal entries to unity.
*/
void TACSSchurMat::applyTransposeBCs(TACSBcMap *bcmap) {
  if (bcmap) {
    /*
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
    */
  }
}

/*
  Get the row/column dimension of the matrix

  output:
  nr:  the row dimension
  nc:  the column dimension
*/
void TACSSchurMat::getSize(int *nr, int *nc) {
  int bs = B->getBlockSize();
  *nr = bs * rmap->getNumNodes();
  *nc = bs * rmap->getNumNodes();
}

/*
  Set the matrix to zero
*/
void TACSSchurMat::zeroEntries() {
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
void TACSSchurMat::copyValues(TACSMat *mat) {
  // Safely down-cast the matrix to an TACSSchurMat - returns NULL
  // if this is not possible
  TACSSchurMat *smat = dynamic_cast<TACSSchurMat *>(mat);
  if (smat) {
    B->copyValues(smat->B);
    E->copyValues(smat->E);
    F->copyValues(smat->F);
    C->copyValues(smat->C);
  } else {
    fprintf(stderr, "Cannot copy matrices of different types\n");
  }
}

/*!
  Scale the entries in the other matrices by a given scalar

  input:
  alpha:  Scale the matrix by alpha: A <- alpha*A
*/
void TACSSchurMat::scale(TacsScalar alpha) {
  B->scale(alpha);
  E->scale(alpha);
  F->scale(alpha);
  C->scale(alpha);
}

/*!
  Compute y <- y + alpha * x
*/
void TACSSchurMat::axpy(TacsScalar alpha, TACSMat *mat) {
  TACSSchurMat *smat = dynamic_cast<TACSSchurMat *>(mat);
  if (smat) {
    B->axpy(alpha, smat->B);
    E->axpy(alpha, smat->E);
    F->axpy(alpha, smat->F);
    C->axpy(alpha, smat->C);
  } else {
    fprintf(stderr, "Cannot apply axpy to matrices of different types\n");
  }
}

/*!
  Compute y <- alpha * x + beta * y
*/
void TACSSchurMat::axpby(TacsScalar alpha, TacsScalar beta, TACSMat *mat) {
  TACSSchurMat *smat = dynamic_cast<TACSSchurMat *>(mat);
  if (smat) {
    B->axpby(alpha, beta, smat->B);
    E->axpby(alpha, beta, smat->E);
    F->axpby(alpha, beta, smat->F);
    C->axpby(alpha, beta, smat->C);
  } else {
    fprintf(stderr, "Cannot apply axpby to matrices of different types\n");
  }
}

/*
  Add a scalar to the diagonal components
*/
void TACSSchurMat::addDiag(TacsScalar alpha) {
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
void TACSSchurMat::mult(TACSVec *txvec, TACSVec *tyvec) {
  tyvec->zeroEntries();

  // Safely down-cast the TACSVec vectors to TACSBVecs
  TACSBVec *xvec, *yvec;
  xvec = dynamic_cast<TACSBVec *>(txvec);
  yvec = dynamic_cast<TACSBVec *>(tyvec);

  if (xvec && yvec) {
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

    C->multAdd(&xlocal[local_offset], &ylocal[local_offset],
               &ylocal[local_offset]);

    // Start sending the values back to y
    c_map->beginReverse(c_ctx, &ylocal[local_offset], y, TACS_ADD_VALUES);
    E->multAdd(&xlocal[local_offset], ylocal, ylocal);

    // Finish transmitting the values back to y
    b_map->beginReverse(b_ctx, ylocal, y, TACS_INSERT_VALUES);
    c_map->endReverse(c_ctx, &ylocal[local_offset], y, TACS_ADD_VALUES);
    b_map->endReverse(b_ctx, ylocal, y, TACS_INSERT_VALUES);
  } else {
    fprintf(stderr, "TACSSchurMat type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrices

  output:
  B, E, F, C: the matrices in the TACSSchurMat class
*/
void TACSSchurMat::getBCSRMat(BCSRMat **_B, BCSRMat **_E, BCSRMat **_F,
                              BCSRMat **_C) {
  if (_B) {
    *_B = B;
  }
  if (_E) {
    *_E = E;
  }
  if (_F) {
    *_F = F;
  }
  if (_C) {
    *_C = C;
  }
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
  smat:    the TACSSchurMat matrix for the preconditioner
  levFill: the level of fill to use
  fill:    the expected/best estimate of the fill-in factor
  reorder: flag to indicate whether to re-order the global Schur complement
*/
TACSSchurPc::TACSSchurPc(TACSSchurMat *_mat, int levFill, double fill,
                         int reorder_schur_complement) {
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
  use_cyclic_alltoall = 0;

  // Perform the symbolic factorization of the [ B, E; F, C ] matrix
  int use_full_schur = 1;  // Use the exact F * B^{-1} * E

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
  TACSNodeMap *rmap = mat->getNodeMap();
  MPI_Comm comm = rmap->getMPIComm();
  Bpc = new BCSRMat(comm, B, E, F, C, levFill, fill, &Epc, &Fpc, &Sc,
                    use_full_schur);
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
  for (int i = size; i > 0; i--) {
    if (ownerRange[i] - ownerRange[i - 1] > 0) {
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
  if (rank == root) {
    schur_count = new int[size];
    schur_ptr = new int[size];
  }

  MPI_Gather(&num_local_schur_vars, 1, MPI_INT, schur_count, 1, MPI_INT, root,
             comm);

  // schur_root is first stored as the list of tacs variable numbers
  // for the schur complement contribution from each processor. After
  // processing on the root node, schur_root is an array of these same
  // variables but now with the order obtained from uniquely sorting
  // the list of the original schur_root variables.
  int *schur_root = NULL;
  int num_schur_root = 0;
  if (rank == root) {
    for (int i = 0; i < size; i++) {
      schur_ptr[i] = num_schur_root;
      num_schur_root += schur_count[i];
    }
    schur_root = new int[num_schur_root];
  }

  // Gather all the global Schur variables to the root process
  MPI_Gatherv((void *)schur_vars, num_local_schur_vars, MPI_INT, schur_root,
              schur_count, schur_ptr, MPI_INT, root, comm);

  // Duplicate the global list and uniquify the result
  int *unique_schur = NULL;
  int num_unique_schur = 0;

  if (rank == root) {
    unique_schur = new int[num_schur_root];
    memcpy(unique_schur, schur_root, num_schur_root * sizeof(int));
    num_unique_schur = TacsUniqueSort(num_schur_root, unique_schur);

    // For each global Schur variable, now assign an output - the index
    // into the unique list of global Schur variables
    for (int i = 0; i < num_schur_root; i++) {
      int *item =
          TacsSearchArray(schur_root[i], num_unique_schur, unique_schur);
      schur_root[i] = item - unique_schur;
    }
  }

  // Broadcast the global size of the unique Schur complement matrix
  MPI_Bcast(&num_unique_schur, 1, MPI_INT, root, comm);

  // Now, pass the unique Schur indices back to the original procs
  // This is an initial ordering of the variables,
  local_schur_vars = new int[num_local_schur_vars];
  MPI_Scatterv(schur_root, schur_count, schur_ptr, MPI_INT, local_schur_vars,
               num_local_schur_vars, MPI_INT, root, comm);

  // Retrieve the non-zero pattern from the local Schur complement
  // and pass it into the block-cyclic matrix
  int bs, num_schur_vars;
  const int *rowp, *cols;
  TacsScalar *vals;
  Sc->getArrays(&bs, &num_schur_vars, &num_schur_vars, &rowp, &cols, &vals);

  // Determine the size of the global Schur complement
  int M = num_unique_schur, N = num_unique_schur;

  // Determine the number of blocks to use per block-cylic block
  int csr_blocks_per_block = 1;
  if (bsize < 36) {
    csr_blocks_per_block = 36 / bsize;
  }

  int *schur_cols = new int[rowp[num_schur_vars]];
  for (int i = 0; i < rowp[num_schur_vars]; i++) {
    schur_cols[i] = local_schur_vars[cols[i]];
  }

  // Create the global block-cyclic Schur complement matrix
  bcyclic = new TACSBlockCyclicMat(
      comm, M, N, bsize, local_schur_vars, num_schur_vars, rowp, schur_cols,
      csr_blocks_per_block, reorder_schur_complement, max_grid_size);
  bcyclic->incref();
  delete[] schur_cols;

  // Get the information about the reordering/blocks from the matrix
  int nrows, ncols;
  const int *bptr, *xbptr, *perm, *iperm, *orig_bptr;
  bcyclic->getBlockPointers(&nrows, &ncols, &bptr, &xbptr, &perm, &iperm,
                            &orig_bptr);
  if (!orig_bptr) {
    orig_bptr = bptr;
  }

  // Set the number of local variables that are defined
  int local_var_count = xbptr[nrows] / bsize;

  // Allocate space for the new pointers
  int *tacs_schur_count = NULL;
  int *tacs_schur_ptr = NULL;
  if (rank == root) {
    tacs_schur_count = new int[size];
    tacs_schur_ptr = new int[size];
  }

  // Gather the local variable count on this processor
  MPI_Gather(&local_var_count, 1, MPI_INT, tacs_schur_count, 1, MPI_INT, root,
             comm);

  // Now, reorder the variables in the Schur complement
  int *local_tacs_schur_vars = NULL;
  if (rank == root) {
    // Compute the offset into the block order
    int num_schur = 0;
    for (int i = 0; i < size; i++) {
      tacs_schur_ptr[i] = num_schur;
      num_schur += tacs_schur_count[i];
    }

    // Find out where to place variable i from the unique list of
    // local schur variables
    local_tacs_schur_vars = new int[num_schur];
    for (int i = 0, j = 0; (i < nrows) && (j < num_unique_schur); i++) {
      while (j < num_unique_schur && bsize * j >= orig_bptr[i] &&
             bsize * j < orig_bptr[i + 1]) {
        // Get the re-ordered block
        int block = i;
        if (iperm) {
          block = iperm[i];
        }
        int owner = bcyclic->get_block_owner(block, block);

        // Count up the number of local blocks to offset.  This is a
        // double loop which could be avoided in future. This might be
        // a bottle neck for very large cases, but just at set up.
        int index = (bsize * j - orig_bptr[i]) / bsize + tacs_schur_ptr[owner];
        for (int k = 0; k < block; k++) {
          if (owner == bcyclic->get_block_owner(k, k)) {
            index += (bptr[k + 1] - bptr[k]) / bsize;
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
    for (int i = 0; i < num_schur_root; i++) {
      schur_root[i] = unique_schur[schur_root[i]];
    }

    // Free the original set of unique schur variables
    delete[] unique_schur;
  }

  // Send unique_schur back to the owning processes
  int *tacs_schur_vars = new int[local_var_count];
  MPI_Scatterv(local_tacs_schur_vars, tacs_schur_count, tacs_schur_ptr, MPI_INT,
               tacs_schur_vars, local_var_count, MPI_INT, root, comm);

  int *local_schur = new int[num_local_schur_vars];
  MPI_Scatterv(schur_root, schur_count, schur_ptr, MPI_INT, local_schur,
               num_local_schur_vars, MPI_INT, root, comm);

  // Free memory not required anymore
  if (rank == root) {
    delete[] tacs_schur_count;
    delete[] tacs_schur_ptr;
    delete[] schur_count;
    delete[] schur_ptr;
    delete[] schur_root;
    delete[] local_tacs_schur_vars;
  }

  // Set up information required for the global Schur complement matrix
  // Set the variable map
  schur_map = new TACSNodeMap(comm, local_var_count);
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
  tacs_schur_dist = new TACSBVecDistribute(mat->getNodeMap(), tacs_schur_index);
  tacs_schur_dist->incref();

  // Create the context for the Schur complement in the global indices
  tacs_schur_ctx = tacs_schur_dist->createCtx(bsize);
  tacs_schur_ctx->incref();

  // Allocate space for local storage of vectors
  int xsize = bsize * b_map->getNumNodes();
  int ysize = bsize * c_map->getNumNodes();
  xlocal = new TacsScalar[xsize];
  yinterface = new TacsScalar[ysize];

  memset(xlocal, 0, xsize * sizeof(TacsScalar));
  memset(yinterface, 0, ysize * sizeof(TacsScalar));

  // Allocate the Schur complement vectors
  yschur = new TACSBVec(schur_map, bsize);
  gschur = new TACSBVec(schur_map, bsize);
  yschur->incref();
  gschur->incref();
}

/*
  Destructor for the TACSSchurPc preconditioner object
*/
TACSSchurPc::~TACSSchurPc() {
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
  bcyclic->decref();

  // Deallocate the local vectors
  gschur->decref();
  yschur->decref();
  delete[] xlocal;
  delete[] yinterface;
}

/*
  Test three different methods for calculating the Schur complement
  contributions.

  This test ensures that the methods are consistent for the given
  matrix. Note that this test can only be performed after the
  preconditioner has been factored since the matrix Sc is not
  populated until this time.
*/
void TACSSchurPc::testSchurComplement(TACSVec *tin, TACSVec *tout) {
  TACSBVec *invec, *outvec;
  invec = dynamic_cast<TACSBVec *>(tin);
  outvec = dynamic_cast<TACSBVec *>(tout);

  if (invec && outvec) {
    // Test two methods of computing the effect of the Schur complement
    outvec->zeroEntries();

    // Get the array
    TacsScalar *in, *out;
    invec->getArray(&in);
    outvec->getArray(&out);

    // Allocate a temporary array to store c-entries
    int bsize = B->getBlockSize();
    int c_size = bsize * c_map->getNumNodes();
    TacsScalar *temp = new TacsScalar[c_size];

    // Comput the schur complement product
    c_map->beginForward(c_ctx, in, yinterface);
    c_map->endForward(c_ctx, in, yinterface);
    Sc->mult(yinterface, temp);
    c_map->beginReverse(c_ctx, temp, out, TACS_ADD_VALUES);
    c_map->endReverse(c_ctx, temp, out, TACS_ADD_VALUES);

    delete[] temp;

    // Compute the schur complement product a second way
    TacsScalar *y;
    yschur->getArray(&y);
    schur_dist->beginReverse(schur_ctx, yinterface, y, TACS_INSERT_VALUES);
    schur_dist->endReverse(schur_ctx, yinterface, y, TACS_INSERT_VALUES);

    TacsScalar *g = NULL;
    yschur->getArray(&y);
    bcyclic->mult(y, g);

    TacsScalar outnorm = outvec->norm();
    TacsScalar gnorm = gschur->norm();

    // Now collect the elements directly using the tacs_schur_dist object
    tacs_schur_dist->beginForward(tacs_schur_ctx, in, y);
    tacs_schur_dist->endForward(tacs_schur_ctx, in, y);

    bcyclic->mult(y, g);
    TacsScalar gnorm2 = gschur->norm();

    int rank;
    MPI_Comm_rank(c_map->getMPIComm(), &rank);

    if (rank == 0) {
      printf("Schur complement consistency test: \n");
      printf("|Full matrix|     = %25.15e \n", TacsRealPart(outnorm));
      printf("|Local to Schur|  = %25.15e \n", TacsRealPart(gnorm));
      printf("|Global to Schur| = %25.15e \n", TacsRealPart(gnorm2));
    }
  } else {
    fprintf(stderr, "TACSSchurPc type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrices

  output:
  B, E, F, C: the matrices in the TACSSchurMat class
*/
void TACSSchurPc::getBCSRMat(BCSRMat **_Bpc, BCSRMat **_Epc, BCSRMat **_Fpc,
                             BCSRMat **_Sc) {
  if (_Bpc) {
    *_Bpc = Bpc;
  }
  if (_Epc) {
    *_Epc = Epc;
  }
  if (_Fpc) {
    *_Fpc = Fpc;
  }
  if (_Sc) {
    *_Sc = Sc;
  }
}

/*
  Set the flag that prints out the factorization time

  input:
  flag: the flag value for the factor-time monitor
*/
void TACSSchurPc::setMonitorFactorFlag(int flag) { monitor_factor = flag; }

/*
  Set the flag that prints out the back solve time

  input:
  flag:  the flag value for the back-solve monitor
*/
void TACSSchurPc::setMonitorBackSolveFlag(int flag) {
  monitor_back_solve = flag;
}

/*
  Set the flag that controls which matrix assembly code to use.

  This flag controls whether the Alltoall version of the block matrix
  assembly is used. When true, this uses a faster, but more
  memory-intensive version of the code. Be aware that this can cause
  the code to run out of memory, but can be very beneficial in terms
  of CPU time.

  input:
  flag:  the flag value to use for the Alltoall flag
*/
void TACSSchurPc::setAlltoallAssemblyFlag(int flag) {
  use_cyclic_alltoall = flag;
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
void TACSSchurPc::factor() {
  // Set the time variables
  double diag_factor_time = 0.0;
  double schur_complement_time = 0.0;
  double global_schur_assembly = 0.0;
  double global_schur_time = 0.0;

  if (monitor_factor) {
    diag_factor_time = -MPI_Wtime();
  }

  // Copy the diagonal matrix B and factor it
  Bpc->copyValues(B);
  Bpc->factor();

  if (monitor_factor) {
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

  if (monitor_factor) {
    schur_complement_time += MPI_Wtime();
    global_schur_assembly = -MPI_Wtime();
  }

  // Assemble the global Schur complement system into block matrix.
  // First, zero the Schur complement matrix
  bcyclic->zeroEntries();

  // Retrieve the local arrays for the local Schur complement
  int bsize, mlocal, nlocal;
  const int *rowp, *cols;
  TacsScalar *scvals;
  Sc->getArrays(&bsize, &mlocal, &nlocal, &rowp, &cols, &scvals);

  int *schur_cols = new int[rowp[mlocal]];
  for (int i = 0; i < rowp[mlocal]; i++) {
    schur_cols[i] = local_schur_vars[cols[i]];
  }

  // Add the values into the global Schur complement matrix
  // using either the alltoall approach or a sequential add values
  // approach that uses less memory
  if (use_cyclic_alltoall) {
    bcyclic->addAlltoallValues(bsize, mlocal, local_schur_vars, rowp,
                               schur_cols, scvals);
  } else {
    bcyclic->addAllValues(bsize, mlocal, local_schur_vars, rowp, schur_cols,
                          scvals);
  }
  delete[] schur_cols;

  if (monitor_factor) {
    global_schur_assembly += MPI_Wtime();
    global_schur_time = -MPI_Wtime();
  }

  // Factor the global Schur complement
  bcyclic->factor();

  if (monitor_factor) {
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
void TACSSchurPc::applyFactor(TACSVec *tin, TACSVec *tout) {
  // First, perform a safe down-cast from TACSVec to BVec
  TACSBVec *invec, *outvec;
  invec = dynamic_cast<TACSBVec *>(tin);
  outvec = dynamic_cast<TACSBVec *>(tout);

  if (invec && outvec) {
    // Set the variables for the back-solve monitor
    double local_time = 0.0;
    double schur_time = 0.0;

    if (monitor_back_solve) {
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

    if (monitor_back_solve) {
      schur_time = -MPI_Wtime();
    }

    // Apply the global Schur complement factorization to the right
    // hand side
    bcyclic->applyFactor(g);

    if (monitor_back_solve) {
      schur_time += MPI_Wtime();
      local_time -= schur_time;
    }

    // copy the values from gschur to yschur
    yschur->copyValues(gschur);

    // Pass yschur solution back to the local variables
    schur_dist->beginForward(schur_ctx, y, yinterface);

    // Pass yschur also to the global variables
    tacs_schur_dist->beginReverse(tacs_schur_ctx, y, out, TACS_INSERT_VALUES);

    // Compute yinterface = yinterface - L^{-1} E y
    // Note: scale y, by -1 first
    schur_dist->endForward(schur_ctx, y, yinterface);
    int one = 1;
    int len = Bpc->getBlockSize() * c_map->getNumNodes();
    TacsScalar alpha = -1.0;
    BLASscal(&len, &alpha, yinterface, &one);

    // Compute xlocal = xlocal - L^{-1} E * yinterface
    Epc->multAdd(yinterface, xlocal, xlocal);

    // Compute xlocal = U^{-1} xlocal
    Bpc->applyUpper(xlocal, xlocal);

    b_map->beginReverse(b_ctx, xlocal, out, TACS_INSERT_VALUES);
    b_map->endReverse(b_ctx, xlocal, out, TACS_INSERT_VALUES);

    // Finish inserting the Schur complement values
    tacs_schur_dist->endReverse(tacs_schur_ctx, y, out, TACS_INSERT_VALUES);

    if (monitor_back_solve) {
      local_time += MPI_Wtime();

      int rank;
      MPI_Comm_rank(b_map->getMPIComm(), &rank);
      printf("[%d] Local back-solve time:  %8.4f\n", rank, local_time);
      printf("[%d] Schur back-solve time:  %8.4f\n", rank, schur_time);
    }
  } else {
    fprintf(stderr, "TACSSchurPc type error: Input/output must be TACSBVec\n");
  }
}

/*
  Retrieve the underlying matrix
*/
void TACSSchurPc::getMat(TACSMat **_mat) { *_mat = mat; }
