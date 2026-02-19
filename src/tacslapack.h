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

#ifndef TACS_LAPACK_H_
#define TACS_LAPACK_H_

/*
  This file contains the definitions of several LAPACK/BLAS functions.
*/

#include "TACSObject.h"

#define LAPACKsyevd dsyevd_
#define LAPACKstevr dstevr_
#define LAPACKdgeev dgeev_
#define LAPACKdspev dspev_
#define LAPACKdspgv dspgv_
#define LAPACKdsbev dsbev_
#define LAPACKdsbgv dsbgv_
#define LAPACKdgelss dgelss_
#define LAPACKdsbgvx dsbgvx_
#define LAPACKdsyev dsyev_
#define LAPACKdpbsv dpbsv_
#define LAPACKzggev zggev_
#define LAPACKdsygvd dsygvd_
#define LAPACKdggev dggev_

#ifdef TACS_USE_COMPLEX
#define BLASdot zdotu_
#define BLASaxpy zaxpy_
#define BLASscal zscal_
#define BLASnrm2 dznrm2_
#define BLAScopy zcopy_
#define BLASgemv zgemv_
#define BLASgbmv zgbmv_
#define BLASsbmv zsbmv_
#define BLASgemm zgemm_
#define BLAStrsm ztrsm_
#define BLAStrsv ztrsv_
#define BLAStbsv ztbsv_
#define LAPACKpbtrf zpbtrf_
#define LAPACKstev zstev_
#define LAPACKgesv zgesv_
#define LAPACKspsv zspsv_
#define LAPACKgetrf zgetrf_
#define LAPACKgetrs zgetrs_
#define LAPACKgetri zgetri_
#else
#define BLASdot ddot_
#define BLASaxpy daxpy_
#define BLASscal dscal_
#define BLASnrm2 dnrm2_
#define BLAScopy dcopy_
#define BLASgemv dgemv_
#define BLASsbmv dsbmv_
#define BLASgbmv dgbmv_
#define BLASgemm dgemm_
#define BLAStrsm dtrsm_
#define BLAStrsv dtrsv_
#define BLAStbsv dtbsv_
#define LAPACKstev dstev_
#define LAPACKpbtrf dpbtrf_
#define LAPACKgesv dgesv_
#define LAPACKspsv dspsv_
#define LAPACKgetrf dgetrf_
#define LAPACKgetrs dgetrs_
#define LAPACKgetri dgetri_
#endif

extern "C" {
// Level 1 BLAS routines
extern TacsScalar BLASdot(int *n, TacsScalar *x, int *incx, TacsScalar *y,
                          int *incy);
extern double BLASnrm2(int *n, TacsScalar *x, int *incx);
extern void BLASaxpy(int *n, TacsScalar *a, TacsScalar *x, int *incx,
                     TacsScalar *y, int *incy);
extern void BLASscal(int *n, TacsScalar *a, TacsScalar *x, int *incx);
extern void BLAScopy(int *n, TacsScalar *x, int *incx, TacsScalar *y,
                     int *incy);  // Copy x into y

// Level 2 BLAS routines
// y = alpha * A * x + beta * y, for a general matrix
extern void BLASgemv(const char *c, int *m, int *n, TacsScalar *alpha,
                     TacsScalar *a, int *lda, TacsScalar *x, int *incx,
                     TacsScalar *beta, TacsScalar *y, int *incy);

// y = alpha * A * x + beta * y, for a symmetric, banded matrix
extern void BLASsbmv(const char *uplo, int *n, int *kb, TacsScalar *alpha,
                     TacsScalar *a, int *lda, TacsScalar *x, int *incx,
                     TacsScalar *beta, TacsScalar *y, int *incy);

// y = alpha * A * x + beta * y, for a general banded matrix
extern void BLASgbmv(const char *trans, int *m, int *n, int *kl, int *ku,
                     TacsScalar *alpha, TacsScalar *a, int *lda, TacsScalar *x,
                     int *incx, TacsScalar *beta, TacsScalar *y, int *incy);

// Solve A*x = b,   or   A**T*x = b in general format
extern void BLAStrsv(const char *uplo, const char *trans, const char *diag,
                     int *n, TacsScalar *a, int *lda, TacsScalar *x, int *incx);

// Solve A*x = b,   or   A**T*x = b in banded format
extern void BLAStbsv(const char *uplo, const char *trans, const char *diag,
                     int *n, int *kd, TacsScalar *a, int *lda, TacsScalar *x,
                     int *incx);

// Level 3 BLAS routines
// C := alpha*op( A )*op( B ) + beta*C,
extern void BLASgemm(const char *ta, const char *tb, int *m, int *n, int *k,
                     TacsScalar *alpha, TacsScalar *a, int *lda, TacsScalar *b,
                     int *ldb, TacsScalar *beta, TacsScalar *c, int *ldc);

// Solve op( A )*X = alpha*B,   or   X*op( A ) = alpha*B
extern void BLAStrsm(const char *side, const char *uplo, const char *transa,
                     const char *diag, int *m, int *n, TacsScalar *alpha,
                     TacsScalar *a, int *lda, TacsScalar *b, int *ldb);

// This routine solves a system of equations with a given number
// of right hand sides
extern void LAPACKgesv(int *n, int *nrhs, TacsScalar *a, int *lda, int *ipiv,
                       TacsScalar *b, int *ldb, int *info);

// Solve a symmetric, indefinite system of linear equations
extern void LAPACKspsv(const char *uplo, int *n, int *nrhs, TacsScalar *a,
                       int *ipiv, TacsScalar *b, int *ldb, int *info);

// Compute an LU factorization of a matrix
extern void LAPACKgetrf(int *m, int *n, TacsScalar *a, int *lda, int *ipiv,
                        int *info);

// Compute a Cholesky factorization of a banded Hermitian positive
// definite matrix
extern void LAPACKpbtrf(const char *uplo, int *n, int *kd, TacsScalar *ab,
                        int *ldab, int *info);

// Given the triangular factorization of the matrix A, compute the inverse
extern void LAPACKgetri(int *n, TacsScalar *a, int *lda, int *ipiv,
                        TacsScalar *work, int *lwork, int *info);

// This routine solves a system of equations with a factored matrix
extern void LAPACKgetrs(const char *c, int *n, int *nrhs, TacsScalar *a,
                        int *lda, int *ipiv, TacsScalar *b, int *ldb,
                        int *info);

// Compute the eigenvalues of a symmetric matrix
extern void LAPACKdsyev(const char *JOBZ, const char *UPLO, int *N, double *A,
                        int *LDA, double *W, double *WORK, int *LWORK,
                        int *INFO);

// Compute the eigenvalues of a symmetric matrix
extern void LAPACKsyevd(const char *jobz, const char *uplo, int *N, double *A,
                        int *lda, double *w, double *work, int *lwork,
                        int *iwork, int *liwork, int *info);

// DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
// of a real generalized symmetric-definite eigenproblem, of the form
// A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
extern void LAPACKdsygvd(int *ITYPE, const char *JOBZ, const char *UPLO, int *N,
                         double *A, int *LDA, double *B, int *LDB, double *W,
                         double *WORK, int *LWORK, int *IWORK, int *LIWORK,
                         int *INFO);

extern void LAPACKdggev(const char *jobvl, const char *jobvr, int *N, double *A,
                        int *lda, double *B, int *ldb, double *alphar,
                        double *alphai, double *beta, double *vl, int *ldvl,
                        double *vr, int *ldvr, double *work, int *lwork,
                        int *info);

// Compute selected eigenvalues of a symmetric tridiagonal matrix
extern void LAPACKstevr(const char *jobz, const char *range, int *n, double *d,
                        double *e, double *vl, double *vu, int *il, int *iu,
                        double *abstol, int *m, double *w, double *z, int *ldz,
                        int *isuppz, double *work, int *lwork, int *iwork,
                        int *liwork, int *info);

// Compute the eigenvalues of a non-symmetric eigenvalue problem
extern void LAPACKdgeev(const char *JOBVL, const char *JOBVR, int *N, double *A,
                        int *LDA, double *WR, double *WI, double *VL, int *LDVL,
                        double *VR, int *LDVR, double *WORK, int *LWORK,
                        int *INFO);

// Compute the eigenvalues of a packed symmetric matrix stored in a packed
extern void LAPACKdspev(const char *JOBZ, const char *UPLO, int *N, double *AP,
                        double *W, double *Z, int *LDZ, double *WORK,
                        int *INFO);

// Solve a symmetric banded positive definite system of equations
extern void LAPACKdpbsv(const char *UPLO, int *N, int *KD, int *NRHS,
                        double *AB, int *LDAB, double *B, int *LDB, int *INFO);

// Compute the eigenvalues of a generalized eigenvalue problem
// with the matrices stored in packed formats
extern void LAPACKdspgv(int *ITYPE, const char *JOBZ, const char *UPLO, int *N,
                        double *AP, double *BP, double *W, double *Z, int *LDZ,
                        double *WORK, int *INFO);

// Solve the eigenvalue problem for a symmetric matrix stored with a
// banded storage format
extern void LAPACKdsbev(const char *JOBZ, const char *UPLO, int *N, int *KD,
                        double *AB, int *LDAB, double *W, double *Z, int *LDZ,
                        double *WORK, int *info);

// Solve the generalized eigenvalue problem for a symmetric matrix
// stored with a banded storage format
extern void LAPACKdsbgv(const char *JOBZ, const char *UPLO, int *N, int *KA,
                        int *KB, double *AB, int *LDAB, double *BB, int *LDBB,
                        double *W, double *Z, int *LDZ, double *WORK,
                        int *info);

// Compute selected eigenvalues and eigenvectors for a symmetric matrix
// storied in banded storage
extern void LAPACKdsbgvx(const char *JOBZ, const char *RANGE, const char *UPLO,
                         int *N, int *KA, int *KB, double *AB, int *LDAB,
                         double *B, int *LDBB, double *Q, int *LDQ, double *VL,
                         double *VU, int *IL, int *IU, double *ABSTOL, int *M,
                         double *W, double *Z, int *LDZ, double *work,
                         int *iwork, int *ifail, int *info);

// Solve an over or underdetermined system of equations
extern void LAPACKdgelss(int *m, int *n, int *nrhs, double *a, int *lda,
                         double *b, int *ldb, double *s, double *rcond,
                         int *rank, double *work, int *lwork, int *info);

// Compute the eigenvalues and optionally the eigenvectors of a
// symmetric, tridiagonal system
extern void LAPACKstev(const char *jobz, int *n, TacsScalar *d, TacsScalar *e,
                       TacsScalar *z, int *ldz, TacsScalar *work, int *info);

// Define the lapack complex type used here
struct LAPACK_cplx_double {
  LAPACK_cplx_double() {
    real = 0.0;
    cplx = 0.0;
  }
  LAPACK_cplx_double(double _real, double _cplx) {
    real = _real;
    cplx = _cplx;
  }
  double real;
  double cplx;
};

extern void LAPACKzggev(const char *JOBVL, const char *JOBVR, int *N,
                        LAPACK_cplx_double *A, int *LDA, LAPACK_cplx_double *B,
                        int *LDB, LAPACK_cplx_double *ALPHA,
                        LAPACK_cplx_double *BETA, LAPACK_cplx_double *VL,
                        int *LDVL, LAPACK_cplx_double *VR, int *LDVR,
                        LAPACK_cplx_double *WORK, int *LWORK, double *RWORK,
                        int *INFO);
}

#endif
