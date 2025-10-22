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

#ifndef TACS_BCSR_MAT_IMPLEMENTATION_H
#define TACS_BCSR_MAT_IMPLEMENTATION_H

/*
  Implementations of the various block-specific operations
*/
#include "TACSObject.h"

class BCSRMatData : public TACSObject {
 public:
  BCSRMatData(int _bsize, int _nrows, int _ncols);
  ~BCSRMatData();

  // After initialization these values should never change!
  int nrows;  // The number of rows
  int ncols;  // The number of columns
  int bsize;  // The matrix block-size
  int *diag;  // A pointer to the diagonal entries
  int *rowp;  // A pointer into cols for the start of each row
  int *cols;  // The column number of each variable

  // Information about the size of groups to use
  int matvec_group_size;  // The size of groups for mat-vec operations
  int matmat_group_size;  // The size of groups for mat-mat operations

  // The storage space for each block - this can change
  TacsScalar *A;  // The vector of elements of each block
};

class BCSRMatThread : public TACSObject {
 public:
  BCSRMatThread(BCSRMatData *_mat);
  ~BCSRMatThread();

  // The scheduling functions
  // ------------------------
  // Matrix-matrix and matrix-vector product scheduler
  void init_mat_mult_sched();
  void mat_mult_sched_job(const int group_size, int *row);
  void mat_mult_sched_job_size(const int group_size, int *row, const int nrows);

  // Scheduler for L^{-1} applications
  void init_apply_lower_sched();
  void apply_lower_sched_job(const int group_size, int *index, int *irow,
                             int *jstart, int *jend);
  void apply_lower_mark_completed(const int group_size, int index, int irow,
                                  int jstart, int jend);

  // Scheduler for U^{-1} appilcations
  void init_apply_upper_sched();
  void apply_upper_sched_job(const int group_size, int *index, int *irow,
                             int *jstart, int *jend);
  void apply_upper_mark_completed(const int group_size, int index, int irow,
                                  int jstart, int jend);

  pthread_t threads[TACSThreadInfo::TACS_MAX_NUM_THREADS];

  // The input/output when dealing with vectors
  TacsScalar *input, *output;

  // BCSRMatData to be operated on by the threaded-implementation
  BCSRMatData *mat;

  // Data required for matrix-matrix products
  double alpha;
  BCSRMatData *Amat;
  BCSRMatData *Bmat;

  // The number of completed rows in the implementation
  int num_completed_rows;

  // Compute the row index size
  int *assigned_row_index;   // The index of the fully assigned rows
  int *completed_row_index;  // The indices of the full assigned columns

  // The threaded implementation
  pthread_mutex_t mutex;
  pthread_cond_t cond;
};

int BMatComputeInverse(TacsScalar *Ainv, TacsScalar *A, int *ipiv, int n);

/*
  The generic implementation. These will be slow.
  Better to implement faster versions for a block size if you need it!
*/

void BCSRMatVecMult(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                       TacsScalar *z);
void BCSRMatVecMultTranspose(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor(BCSRMatData *A);
void BCSRMatApplyLower(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd(double alpha, BCSRMatData *A, BCSRMatData *B,
                       BCSRMatData *C);
void BCSRMatFactorLower(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper(BCSRMatData *A, BCSRMatData *F);
void BCSRMatMatMultNormal(BCSRMatData *A, TacsScalar *s, BCSRMatData *B);
void BCSRMatApplySOR(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                     const int end, const int var_offset,
                     const TacsScalar *Adiag, const TacsScalar omega,
                     const TacsScalar *b, const TacsScalar *xext,
                     TacsScalar *x);

/*
  These are the definitions for the block-specific code.
  These should provide better performance.
*/

// The bsize == 1 code
void BCSRMatVecMult1(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd1(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose1(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor1(BCSRMatData *A);
void BCSRMatApplyLower1(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper1(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower1(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper1(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur1(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd1(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower1(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper1(BCSRMatData *A, BCSRMatData *F);
void BCSRMatMatMultNormal1(BCSRMatData *A, TacsScalar *s, BCSRMatData *B);
void BCSRMatApplySOR1(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The bsize == 2 code
void BCSRMatVecMult2(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd2(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose2(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor2(BCSRMatData *A);
void BCSRMatApplyLower2(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper2(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower2(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper2(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur2(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd2(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower2(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper2(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR2(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The bsize == 3 code
void BCSRMatVecMult3(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd3(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose3(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor3(BCSRMatData *A);
void BCSRMatApplyLower3(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper3(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower3(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper3(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur3(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd3(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower3(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper3(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR3(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The bsize == 4 code
void BCSRMatVecMult4(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd4(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose4(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor4(BCSRMatData *A);
void BCSRMatApplyLower4(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper4(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower4(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper4(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur4(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd4(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower4(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper4(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR4(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The bsize == 5 code
void BCSRMatVecMult5(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd5(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose5(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor5(BCSRMatData *A);
void BCSRMatApplyLower5(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper5(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower5(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper5(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur5(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd5(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower5(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper5(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR5(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The bsize == 6 code
void BCSRMatVecMult6(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd6(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose6(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor6(BCSRMatData *A);
void BCSRMatApplyLower6(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper6(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower6(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper6(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur6(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd6(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower6(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper6(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR6(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The threaded implementations for bsize = 6
void *BCSRMatVecMultAdd6_thread(void *t);
void *BCSRMatFactor6_thread(void *t);
void *BCSRMatFactorLower6_thread(void *t);
void *BCSRMatFactorUpper6_thread(void *t);

void *BCSRMatApplyLower6_thread(void *t);
void *BCSRMatApplyUpper6_thread(void *t);

void *BCSRMatMatMultAdd6_thread(void *t);

// The bsize == 8 code
void BCSRMatVecMult8(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatVecMultAdd8(BCSRMatData *A, TacsScalar *x, TacsScalar *y,
                        TacsScalar *z);
void BCSRMatVecMultTranspose8(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatFactor8(BCSRMatData *A);
void BCSRMatApplyLower8(BCSRMatData *A, TacsScalar *x, TacsScalar *y);
void BCSRMatApplyUpper8(BCSRMatData *A, TacsScalar *x, TacsScalar *y);

void BCSRMatApplyPartialLower8(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyPartialUpper8(BCSRMatData *A, TacsScalar *x, int var_offset);
void BCSRMatApplyFactorSchur8(BCSRMatData *A, TacsScalar *x, int var_offset);

void BCSRMatMatMultAdd8(double alpha, BCSRMatData *A, BCSRMatData *B,
                        BCSRMatData *C);
void BCSRMatFactorLower8(BCSRMatData *A, BCSRMatData *E);
void BCSRMatFactorUpper8(BCSRMatData *A, BCSRMatData *F);
void BCSRMatApplySOR8(BCSRMatData *Adata, BCSRMatData *Bdata, const int start,
                      const int end, const int var_offset,
                      const TacsScalar *Adiag, const TacsScalar omega,
                      const TacsScalar *b, const TacsScalar *xext,
                      TacsScalar *x);

// The threaded implementations for bsize = 8
void *BCSRMatVecMultAdd8_thread(void *t);
void *BCSRMatFactor8_thread(void *t);
void *BCSRMatFactorLower8_thread(void *t);
void *BCSRMatFactorUpper8_thread(void *t);

void *BCSRMatApplyLower8_thread(void *t);
void *BCSRMatApplyUpper8_thread(void *t);

void *BCSRMatMatMultAdd8_thread(void *t);

#endif
