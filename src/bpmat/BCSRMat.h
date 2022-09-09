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

#ifndef TACS_BCSR_MATRIX_H
#define TACS_BCSR_MATRIX_H

#include "BCSRMatImpl.h"
#include "KSM.h"
#include "TACSObject.h"
#include "pthread.h"

/*!
  Block Compressed Sparse Row format matrix.

  This matrix requires that the blocks be all the same size.
  Different, optimized code is used for performance-critical
  operations for each block size.
*/
class BCSRMat : public TACSObject {
 public:
  // Factor the matrix - construct the non--zero pattern from an ILU(p)
  BCSRMat(MPI_Comm _comm, BCSRMat *mat, int levFill, double fill,
          const char *fname = NULL);

  // Create a (possibly rectangular matrix)
  BCSRMat(MPI_Comm _comm, TACSThreadInfo *_thread_info, int _bsize, int _nrows,
          int _ncols, int **_rowp, int **_cols, TacsScalar **_A = NULL);

  // Perform the symbolic computation C = S + A*B
  BCSRMat(MPI_Comm _comm, BCSRMat *Smat, int *alevs, BCSRMat *Amat, int *blevs,
          BCSRMat *Bmat, int levFill, double fill);

  // Perform a symbolic computation of the Schur complement
  BCSRMat(MPI_Comm _comm, BCSRMat *Bmat, BCSRMat *Emat, BCSRMat *Fmat,
          BCSRMat *Cmat, int levFill, double fill, BCSRMat **Epc, BCSRMat **Fpc,
          BCSRMat **Sc, int use_full_schur);

  // Perform a symbolic computation of A = B^{T}*B
  BCSRMat(MPI_Comm comm, BCSRMat *Bmat, double fill);

  // Free the matrix
  ~BCSRMat();

  // Access data
  // -----------
  MPI_Comm getMPIComm();
  TACSThreadInfo *getThreadInfo();

  // Functions for setting values in the matrix
  // ------------------------------------------
  void zeroEntries();
  void addRowValues(int row, int ncol, const int *col, int nca,
                    const TacsScalar *avals);
  void addRowWeightValues(TacsScalar alpha, int row, int nwrows,
                          const int *wrowp, const int *wcols,
                          const TacsScalar *weights, int nca,
                          const TacsScalar *avals,
                          MatrixOrientation matOr = TACS_MAT_NORMAL);
  void addBlockRowValues(int row, int ncol, const int *col,
                         const TacsScalar *a);
  void zeroRow(int row, int vars, int ident = 0);
  void zeroColumns(int num_zero_cols, const int *zero_cols,
                   const int *zero_vars, int ident = 0);
  void getArrays(int *_bsize, int *nrows, int *ncols, const int **rowp,
                 const int **cols, TacsScalar **Avals);

  void copyValues(BCSRMat *mat);
  void scale(TacsScalar alpha);
  void axpy(TacsScalar alpha, BCSRMat *x);
  void axpby(TacsScalar alpha, TacsScalar beta, BCSRMat *x);
  void addDiag(TacsScalar alpha);
  void addDiag(TacsScalar *alpha);

  // Functions related to solving the system of equations
  // ----------------------------------------------------
  void factor();
  void mult(TacsScalar *xvec, TacsScalar *yvec);
  void multAdd(TacsScalar *xvec, TacsScalar *zvec, TacsScalar *yvec);
  void multTranspose(TacsScalar *xvec, TacsScalar *yvec);
  void applyFactor(TacsScalar *xvec, TacsScalar *yvec);
  void applyFactor(TacsScalar *xvec);
  void applyUpper(TacsScalar *xvec, TacsScalar *yvec);
  void applyLower(TacsScalar *xvec, TacsScalar *yvec);
  void applyPartialLower(TacsScalar *xvec, int var_offset);
  void applyPartialUpper(TacsScalar *xvec, int var_offset);
  void applyFactorSchur(TacsScalar *x, int var_offset);
  void setDiagPairs(const int *_pairs, int _npairs);
  void factorDiag(const TacsScalar *diag = NULL);
  void applySOR(TacsScalar *x, TacsScalar *y, TacsScalar omega, int iters);
  void applySOR(BCSRMat *B, int start, int end, int var_offset,
                TacsScalar omega, const TacsScalar *b, const TacsScalar *xext,
                TacsScalar *x);

  void matMultAdd(double alpha, BCSRMat *amat, BCSRMat *bmat);
  void applyLowerFactor(BCSRMat *emat);
  void applyUpperFactor(BCSRMat *fmat);

  // Compute the normal equations A = B^{T}*S*B where S = diag{s}
  // ------------------------------------------------------------
  void matMultNormal(TacsScalar *s, BCSRMat *bmat);

  // Get the matrix dimensions
  // -------------------------
  int getBlockSize() { return data->bsize; }
  int getRowDim() { return data->nrows; }
  int getColDim() { return data->ncols; }
  BCSRMatData *getMatData() { return data; }

  // Extract the matrix in a  LAPACK format
  // --------------------------------------
  void getNumUpperLowerDiagonals(int *_bl, int *_bu);
  void getDenseColumnMajor(TacsScalar *A);
  void getBandedMatrix(TacsScalar *A, int size, int symm_flag);

  // Other miscelaneous functions
  // ----------------------------
  int isEqual(BCSRMat *mat, double tol = 1e-12);
  void partition(int Np, BCSRMat **Bmat, BCSRMat **Emat, BCSRMat **Fmat,
                 BCSRMat **Cmat);
  void testSchur(int Np, int lev, double fill, double tol = 1e-12);
  void printMat(const char *fname);
  void printNzPattern(const char *fname);

  // Initialize either the generic or block-specific implementations
  // ---------------------------------------------------------------
  void initGenericImpl();
  void initBlockImpl();

 private:
  void setUpDiag();  // Set up the diagonal entry pointer 'diag'
  void computeILUk(BCSRMat *mat, int levFill, double fill, int **_levs);
  BCSRMat *computeILUkEpc(BCSRMat *EMat, const int *levs, int levFill,
                          double fill, int **_elevs);
  BCSRMat *computeILUkFpc(BCSRMat *EMat, const int *levs, int levFill,
                          double fill, int **_flevs);

  // (potential) block-specific implementations for each operation
  void (*bmult)(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
  void (*bmultadd)(BCSRMatData *A, TacsScalar *x, TacsScalar *y, TacsScalar *z);
  void (*bmulttrans)(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
  void (*bfactor)(BCSRMatData *A);
  void (*applylower)(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
  void (*applyupper)(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

  void (*applypartialupper)(BCSRMatData *A, TacsScalar *x, int var_offset);
  void (*applypartiallower)(BCSRMatData *A, TacsScalar *x, int var_offset);
  void (*applyschur)(BCSRMatData *A, TacsScalar *x, int var_offset);
  void (*applysor)(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                   const int end, const int var_offset, const TacsScalar *Adiag,
                   const TacsScalar omega, const TacsScalar *b,
                   const TacsScalar *xext, TacsScalar *x);

  void (*bmatmult)(double alpha, BCSRMatData *A, BCSRMatData *B,
                   BCSRMatData *C);
  void (*bfactorlower)(BCSRMatData *A, BCSRMatData *E);
  void (*bfactorupper)(BCSRMatData *A, BCSRMatData *F);
  void (*bmatmatmultnormal)(BCSRMatData *A, TacsScalar *s, BCSRMatData *B);

  // The thread-specific implementations
  void *(*bmultadd_thread)(void *);
  void *(*bfactor_thread)(void *);
  void *(*applylower_thread)(void *);
  void *(*applyupper_thread)(void *);
  void *(*bmatmult_thread)(void *);
  void *(*bfactorlower_thread)(void *);
  void *(*bfactorupper_thread)(void *);

  // The MPI communicator - this is useful for output
  MPI_Comm comm;

  // Information about the threaded execution
  TACSThreadInfo *thread_info;

  // The matrix BCSR data
  BCSRMatData *data;

  // The threaded BCSR matrix data
  BCSRMatThread *tdata;

  // Storage space for the factored diagonal entries
  TacsScalar *Adiag;
  int npairs;
  int *pairs;
};

#endif  // TACS_BCSR_MATRIX_H
