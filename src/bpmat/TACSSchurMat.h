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

#ifndef TACS_SCHUR_MATRIX_H
#define TACS_SCHUR_MATRIX_H

/*
  Matrix and preconditioner based upon the use of the Schur-complement
*/

#include "BCSRMat.h"
#include "TACSBVec.h"
#include "TACSBVecDistribute.h"
#include "TACSBlockCyclicMat.h"

/*!
  A class for a distributed finite-element matrix.

  This matrix can be used in conjunction with a "preconditioner" that
  is a full factorization that uses static pivoting. This uses a Schur
  complement approach, that is sometimes called sub-structuring, where
  we are limited to one sub-structure per processor. The
  "sub-structure" is the entire domain of the finite-elements owned by
  this processor. The advantage of this approach is that the size of
  the interface domain that forms the global Schur complement, only
  depends on the size of the interface, and not on the element
  order. This reduces the number of degrees of freedom in the Schur
  complement. This matrix stores the equations in the following form
  on each processor:

  [ B, E ][ x ]   [ f ]
  [ F, C ][ y ] = [ g ]

  where x are the local unknowns and y is some portion of a global
  vector. Note that the y variables are non-unique between processors,
  but x is unique. Therefore, we can eliminate the local variables and
  form an equation for all y globally across the finite-element mesh.
  Each processor contributes its own Schur complement:

  C - F*B^{-1}*E

  to the global Schur complement system.
*/
class TACSSchurMat : public TACSMat {
 public:
  TACSSchurMat(TACSNodeMap *_rmap, BCSRMat *_B, BCSRMat *_E, BCSRMat *_F,
               BCSRMat *_C, TACSBVecDistribute *_b_map,
               TACSBVecDistribute *_c_map);
  TACSSchurMat(TACSThreadInfo *thread_info, TACSNodeMap *_rmap, int bsize,
               int nlocal_vars, const int *rowp, const int *cols,
               TACSBVecIndices *b_local_indices, TACSBVecDistribute *_b_map,
               TACSBVecIndices *c_local_indices, TACSBVecDistribute *_c_map);
  ~TACSSchurMat();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void addValues(int nrow, const int *row, int ncol, const int *col, int nv,
                 int mv, const TacsScalar *values);
  void addWeightValues(int nvars, const int *varp, const int *vars,
                       const TacsScalar *weights, int nv, int mv,
                       const TacsScalar *values,
                       MatrixOrientation matOr = TACS_MAT_NORMAL);
  void applyBCs(TACSBcMap *bcmap);
  void applyTransposeBCs(TACSBcMap *bcmap);
  TACSVec *createVec();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();

  // Functions required for solving linear systems
  // ---------------------------------------------
  void getSize(int *_nr, int *_nc);
  void mult(TACSVec *x, TACSVec *y);
  void copyValues(TACSMat *mat);
  void scale(TacsScalar alpha);
  void axpy(TacsScalar alpha, TACSMat *x);
  void axpby(TacsScalar alpha, TacsScalar beta, TACSMat *x);
  void addDiag(TacsScalar alpha);

  // Get the underlying representation for ScMat
  // -------------------------------------------
  void getBCSRMat(BCSRMat **_B, BCSRMat **_E, BCSRMat **_F, BCSRMat **_C);

  TACSBVecDistribute *getLocalMap() { return b_map; }
  TACSBVecDistribute *getSchurMap() { return c_map; }
  TACSNodeMap *getNodeMap() { return rmap; }

 protected:
  TACSSchurMat();
  void init(TACSNodeMap *_rmap, BCSRMat *_B, BCSRMat *_E, BCSRMat *_F,
            BCSRMat *_C, TACSBVecDistribute *_b_map,
            TACSBVecDistribute *_c_map);

  BCSRMat *B, *E, *F, *C;     // The local matrices
  TACSNodeMap *rmap;          // Global row map
  TACSBVecDistribute *b_map;  // Collect variables for B/E
  TACSBVecDistribute *c_map;  // Collect variables for C/F
  TACSBVecDistCtx *b_ctx, *c_ctx;

 private:
  // The local sizes/offsets
  int local_size, local_offset;

  // The local variables used for mat-vec products
  TacsScalar *xlocal, *ylocal;

  static const char *matName;
};

/*!
  The global Schur complement preconditioner.

  This preconditioner can be used to form a full factorization with
  static pivoting for the ScMat class defined above. This uses the
  decomposition of the matrix into local and global contributions to
  form a local Schur complement on each processor as follows:

  C - F*U_{B}^{-1}*L_{B}^{-1}*E

  where B = L_{B}*U_{B} is the LU decomposition of B (with static
  pivoting.) This preconditioner forms the full, exact Schur
  complement from all processors and computes the full factorization
  using a block-cyclic approach with the PDMat class. This is an
  efficient sparse parallel block-cyclic factorization code that
  employs static pivoting. The Schur complement is significantly more
  filled-in than the original matrices, so this factorization often
  performs extremely well in parallel, and is much smaller than the
  full matrix.
*/
class TACSSchurPc : public TACSPc {
 public:
  TACSSchurPc(TACSSchurMat *_mat, int levFill, double fill,
              int reorder_schur_complement);
  ~TACSSchurPc();

  // Functions associated with the factorization
  // -------------------------------------------
  void factor();
  void applyFactor(TACSVec *xvec, TACSVec *yvec);
  void getMat(TACSMat **_mat);
  void testSchurComplement(TACSVec *in, TACSVec *out);

  // Monitor the factorization time on each process
  // ----------------------------------------------
  void setMonitorFactorFlag(int flag);
  void setMonitorBackSolveFlag(int flag);

  // Set the type of matrix assembly to use
  // --------------------------------------
  void setAlltoallAssemblyFlag(int flag);

  // Get the underlying precondition representation
  // ----------------------------------------------
  void getBCSRMat(BCSRMat **_Bpc, BCSRMat **_Epc, BCSRMat **_Fpc,
                  BCSRMat **_Sc);

 private:
  TACSSchurMat *mat;
  BCSRMat *B, *E, *F, *C;    // The block matrices
  BCSRMat *Bpc, *Epc, *Fpc;  // The diagonal contributions
  BCSRMat *Sc;               // Schur complement from this proc.

  TACSBVecDistribute *b_map;  // The map for the local entries
  TACSBVecDistribute *c_map;  // The map for the Schur complement
  TACSBVecDistCtx *b_ctx;
  TACSBVecDistCtx *c_ctx;

  int monitor_factor;      // Monitor the factorization time
  int monitor_back_solve;  // Monitor the back-solves

  // The sparse block cyclic matrix
  TACSBlockCyclicMat *bcyclic;  // This stores the Schur complement

  // This set of indices defines the ordering passed into the
  // parallel block-cyclic matrix factorization
  int num_local_schur_vars;
  int *local_schur_vars;

  // There are two Schur maps here:
  // These objects defines a mapping between the local C-variables
  // (in xlocal/ylocal) to the global Schur complement system (scmat).
  TACSNodeMap *schur_map;          // The variable map associated with Sc
  TACSBVecDistribute *schur_dist;  // Map that distributes the Schur complement
  TACSBVecDistCtx *schur_ctx;      // The context for the distribution object
  int use_cyclic_alltoall;  // Use the Alltoall version for matrix assembly

  // This object defines a mapping between the variables in the
  // global vectors (from ScMat - in/out in applyFactor) and the
  // global Schur complement system (scmat).
  TACSBVecDistribute *tacs_schur_dist;
  TACSBVecDistCtx *tacs_schur_ctx;

  TacsScalar *xlocal;         // The local variables
  TacsScalar *yinterface;     // The interface variables
  TACSBVec *gschur, *yschur;  // The Schur complement vectors
};

#endif  // TACS_SCHUR_MATRIX_H
