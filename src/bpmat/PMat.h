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

#ifndef TACS_PARALLEL_MATRIX_H
#define TACS_PARALLEL_MATRIX_H

#include "BVec.h"
#include "BVecDist.h"
#include "BCSRMat.h"
#include "KSM.h"

/*!
  Parallel matrix class with Block Compressed-Sparse row format.

  Each processor contains the following information:

  1. The local matrix that acts on the local variables (A_loc)
  2. The external-coupling matrix (B_ext)
  3. A map that takes the global vector x and maps it to the local external
  vector components: x_ext

  Matrix-vector products are done row-wise:

  y_loc = A_loc * x_loc + B_ext * x_ext

  Matrix-vector products
  ----------------------
  1. Begin scatter operation from global to local-external unknowns
  2. Perform local matrix-vector product
  3. End scatter operation 
  4. Perform local, exteral matrix-vector product  
*/
class TACSPMat : public TACSMat {
 public:
  TACSPMat( TACSVarMap *rmap,
            BCSRMat *_Aloc, BCSRMat *_Bext,
            TACSBVecDistribute *_col_map );
  ~TACSPMat();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();                // Zero the matrix values 
  void applyBCs( TACSBcMap *bcmap ); // Apply the boundary conditions

  // Functions required for solving linear systems
  // ---------------------------------------------
  void getSize( int *_nr, int *_nc );          // Get the local dimensions
  void mult( TACSVec *x, TACSVec *y );         // y <- A*x
  void multTranspose( TACSVec *x, TACSVec *y ); // y <- A^{T}*x
  TACSVec *createVec();                        // Create a vector
  void copyValues( TACSMat *mat );             // Copy matrix entries
  void scale( TacsScalar alpha );              // Scale the matrix
  void axpy( TacsScalar alpha, TACSMat *mat ); // Compute y <- y + alpha*x
  void axpby( TacsScalar alpha, 
              TacsScalar beta, TACSMat *mat ); // Compute y <- alpha*x + beta*y
  void addDiag( TacsScalar alpha );

  // Other miscelaneous functions
  // ----------------------------
  void getBCSRMat( BCSRMat ** A, BCSRMat ** B );  // Access the underlying mats
  void getRowMap( int *bs, int *_N, int *_Nc );
  void getColMap( int *bs, int *_M );
  TACSVarMap *getRowMap(){ return rmap; }
  void getExtColMap( TACSBVecDistribute **ext_map ); // Access the column map  
  void printNzPattern( const char *fileName ); // Print the non-zero pattern

  const char *TACSObjectName();
  MPI_Comm getMPIComm(){ return rmap->getMPIComm(); }

 protected:
  TACSPMat();

  // Common initialization routine
  void init( TACSVarMap *_rmap, 
             BCSRMat *_Aloc, BCSRMat *_Bext,
             TACSBVecDistribute *_col_map );

  // Local entries for the matrix
  BCSRMat *Aloc, *Bext;
    
  // Map the local entries into the global data
  TACSVarMap *rmap;
  TACSBVecDistribute *ext_dist;
  TACSBVecDistCtx *ctx;

  // Sizes/dimensions of the matrix
  int bsize; // The block size
  int N; // Number of local rows
  int Nc; // Number of equations that are coupled to other processors
  int Np; // The number of local-only equations Np + Nc = N

 private: 
  // External values - used for matrix-vector products
  TacsScalar *x_ext; 
  int ext_offset;

  static const char *matName;
};

/*
  Parallel Gauss--Seidel/SOR

  This uses repeated applications of block Gauss--Seidel. Off-processor 
  updates are delayed, effectively making this a hybrid Jacobi Gauss-Seidel
  method.
*/
class TACSGaussSeidel : public TACSPc {
 public:
  TACSGaussSeidel( TACSPMat *_mat, int _zero_guess, 
                   TacsScalar _omega, int _iters, int _symmetric,
                   int _use_l1_gauss_seidel=0 );
  ~TACSGaussSeidel();

  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void getMat( TACSMat **_mat );

 private:
  // Parallel matrix pointer
  TACSPMat *mat;

  // Information about how to handle the smoother
  int iters; // The number of iterations to apply
  TacsScalar omega; // The over/under relaxation factor
  int zero_guess; // Zero the initial guess
  int symmetric; // Apply the symmetric variant 
  int use_l1_gauss_seidel; // Apply the L1 variant of Gauss-Seidel

  // Pointers to the local/external matrix
  BCSRMat *Aloc, *Bext;
  TACSBVec *bvec;

  // Parallel data for the PSOR object
  TACSBVecDistribute *ext_dist;
  TACSBVecDistCtx *ctx;
  int ext_offset;
  TacsScalar *yext;
};

/*
  Chebyshev Smoother
*/
class TACSChebyshevSmoother : public TACSPc {
 public:
  TACSChebyshevSmoother( TACSPMat *_mat, int _degree,
                         double _lower_factor=1.0/30.0,
                         double _upper_factor=1.1, int _iters=1 );
  ~TACSChebyshevSmoother();

  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void getMat( TACSMat **_mat );

 private:
  // Estimate the spectral radius using Gershgorin method
  double gershgorin();

  // Estimate the spectral radius using Arnoldi
  double arnoldi( int size );

  // Parallel matrix pointer
  TACSPMat *mat;

  // The factor to apply to the largest eigenvalue
  double lower_factor;
  double upper_factor;

  // Set the values for the upper/lower eigenvalues
  double alpha, beta;

  // The number of iterations to apply
  int iters;

  // The degree of the polynomial, the roots and coefficients
  int degree;
  double *r, *c;

  // Temporary vectors
  TACSBVec *res, *t, *h;
};

/*
  Additive Schwarz Method (ASM)

  Set up involves factoring the diagonal portion of the matrix.  Apply
  the local preconditioner to the local components of the residual.
*/
class TACSAdditiveSchwarz : public TACSPc {
 public:
  TACSAdditiveSchwarz( TACSPMat *mat, int levFill, double fill );
  ~TACSAdditiveSchwarz();

  void setDiagShift( TacsScalar _alpha );
  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void applyFactor( TACSVec *yvec );
  void getMat( TACSMat **_mat );

 private:
  TACSPMat *mat;
  BCSRMat *Aloc;
  TacsScalar alpha;
  BCSRMat *Apc;
};

/*
  Approximate Schur

  Set up involves factoring the diagonal portion of the matrix

  AS: Apply preconditioner
  1. Restrict to the interface unknowns
  2. Solve a GMRES-accelerated, Jacobi-preconditioned problem for the
  interface unknowns 
  3. Determine the solution at the internal interface unknowns
*/
class TACSGlobalSchurMat : public TACSMat {
 public:
  TACSGlobalSchurMat( TACSPMat *mat, BCSRMat *Apc );
  ~TACSGlobalSchurMat();

  // Functions used to solve the linear system
  // -----------------------------------------
  void getSize( int *_nr, int *_nc );
  void mult( TACSVec *x, TACSVec *y );
  TACSVec *createVec();

  // Multiply y <- Bext *x_ext 
  void multOffDiag( TACSBVec *x, TACSBVec *y );

 private:
  TACSVarMap *rmap; // The variable map for the interface variables
  BCSRMat *Apc; // The factored diagonal matrix 
  BCSRMat *Bext; // The off-diagonal part
  TACSBVecDistribute *ext_dist; 
  TACSBVecDistCtx *ctx;

  int nvars, varoffset;
  TacsScalar *x_ext;
};

/*!
  The approximate Schur preconditioner

  This preconditioner solves an approximate Schur complement system at
  each iteration. The Schur complement system is formed globally in
  parallel across all processors.
*/
class TACSApproximateSchur : public TACSPc {
 public:
  TACSApproximateSchur( TACSPMat *mat, int levFill, double fill, 
                        int inner_gmres_iters, double inner_rtol=1e-3, 
                        double inner_atol=1e-30 );
  ~TACSApproximateSchur();

  void setDiagShift( TacsScalar _alpha );
  void setMonitor( KSMPrint *ksm_print );
  void factor();
  void applyFactor( TACSVec *xvec, TACSVec *yvec );
  void getMat( TACSMat **_mat );
  void printNzPattern( const char *fileName );

 private:
  TACSBVec *rvec, *wvec;
  TACSPMat *mat;
  BCSRMat *Aloc, *Apc;
  TacsScalar alpha;
 
  // Offsets into the array
  int start, end, var_offset;

  // Global Schur matrix and its associated KSM object
  TACSGlobalSchurMat *gsmat;
  TACSKsm *inner_ksm;
};

#endif // TACS_PARALLEL_MATRIX_H
