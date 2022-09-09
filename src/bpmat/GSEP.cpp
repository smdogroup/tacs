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

#include "GSEP.h"

#include "tacslapack.h"

/*
  The following file contains the definitions for two main types of classes:

  EPOperator - Defines the eigen problem

  This is a abstract base class that defines an eigenvalue problem with a
  number of properties. It defines an operator which allows a matrix to
  be modified with a spectral transformation. The shifted/inverted eigenvalue
  is extracted from the operator. The class also provides a dot product that
  is used in the construction of the orthogonal basis of the Krylov subspace.
  In the case of symmetric generalized eigenvalue problems where B is symmetric
  this dot product may be replaced with the inner product <x,y> = x^{T} B y.

  SEP - Iterative eigenvalue solver for symmetric eigenvalue problems

  This class generates an orthogonal basis for a Krylov subpsace using
  either Lanczos or full Gram-Schmidt orthogonalization. The Gram-Schmidt
  option is to be prefered. The additional cost is out-weighted by the better
  numerical properties of the full-orthogonalization. The initial subspace
  vector is generated randomly. All subsequent vectors are generated through
  matrix-vector products with the EPOperator class defined above. This operator
  may involve a shifted/inverted eigen problem and so may require the solution
  of a linear system of equations.
*/

/*
  A regular implementation of the operator class - use the
  matrix directly as the operator. This is only efficient for
  determining the eigenvalues with large relative separation. For
  structural problems these are typically the largest eigenvalues
  - the ones that are "uninteresting". This can, however, be used
  to determine the condition number of a system of equations.
*/
EPRegular::EPRegular(TACSMat *_mat) {
  mat = _mat;
  mat->incref();
}

EPRegular::~EPRegular() { mat->decref(); }

TACSVec *EPRegular::createVec() { return mat->createVec(); }

void EPRegular::mult(TACSVec *x, TACSVec *y) { return mat->mult(x, y); }

/*
  Shift and invert operator. This requires that the matrix be shifted
  by sigma away from origin. The KSM object should represent a
  solution to the linear system of equations

  y = (A - sigma I)^{-1} x

  for some value of sigma. This method is intendend to provide a way
  of extracting the eigenvalues close to sigma. The shift and invert
  strategy modifies the eigenvalues so that the modified eigenproblem
  has eigenvalues of,

  mu = 1.0/( lambda - sigma )

  This modification improves the separation of the eigenvalues near
  sigma, allowing them to converge faster. The eigenvectors of the
  modified problem are the same as the eigenvectors of the original
  problem.
*/
EPShiftInvert::EPShiftInvert(TacsScalar _sigma, TACSKsm *_ksm) {
  sigma = _sigma;
  ksm = _ksm;
  ksm->incref();
}

EPShiftInvert::~EPShiftInvert() { ksm->decref(); }

TACSVec *EPShiftInvert::createVec() { return ksm->createVec(); }

void EPShiftInvert::mult(TACSVec *x, TACSVec *y) { return ksm->solve(x, y); }

// The eigenvalues are computed as mu = 1.0/( eig - sigma )
// eig = 1.0/mu + sigma
TacsScalar EPShiftInvert::convertEigenvalue(TacsScalar value) {
  return (1.0 / value + sigma);
}

/*
  A shift and invert operator for genearlized eigenvalue problems.
  This operator transforms the initial generalized eigenvalue problem,

  A x = lambda B x

  into the following shifted and inverted problem,

  mu x = ( A - sigma B )^{-1} B x
  mu = 1.0/( lambda - sigma )

  The eigenvectors of the modified problem are the same as the
  eigenvectors of the initial problem.  Symmetry of the eigenvalue
  problem is maintained by utilizing the inner product <x,y> = x^{T} B
  y.  The Krylov basis is thus B-conjugate.

  This modification to the Lanczos algorithm is often called B-Lanczos
  in the literature.
*/
EPGeneralizedShiftInvert::EPGeneralizedShiftInvert(TacsScalar _sigma,
                                                   TACSKsm *_ksm,
                                                   TACSMat *_inner) {
  sigma = _sigma;
  ksm = _ksm;
  ksm->incref();
  inner = _inner;
  inner->incref();
  temp = inner->createVec();
  temp->incref();
}

EPGeneralizedShiftInvert::~EPGeneralizedShiftInvert() {
  ksm->decref();
  inner->decref();
  temp->decref();
}

/*
  Set the shift and invert value
*/
void EPGeneralizedShiftInvert::setSigma(TacsScalar _sigma) { sigma = _sigma; }

/*
  Create a vector associated with the generalized eigenvalue problem
*/
TACSVec *EPGeneralizedShiftInvert::createVec() { return ksm->createVec(); }

/*
  Compute y = ( A - sigma B )^{-1}*inner x
*/
void EPGeneralizedShiftInvert::mult(TACSVec *x, TACSVec *y) {
  inner->mult(x, temp);
  return ksm->solve(temp, y);
}

/*
  Compute <x,y> = x^{T} inner y
*/
TacsScalar EPGeneralizedShiftInvert::dot(TACSVec *x, TACSVec *y) {
  inner->mult(y, temp);
  return temp->dot(x);
}

/*
  Compute ||B*x|| - this is used to compute the eigenvalue error
*/
TacsScalar EPGeneralizedShiftInvert::errorNorm(TACSVec *x) {
  inner->mult(x, temp);
  return temp->norm();
}

/*
  Convert the shifted eigenvalue to the actual eigenproblem
*/
TacsScalar EPGeneralizedShiftInvert::convertEigenvalue(TacsScalar value) {
  return (1.0 / value + sigma);
}

/*
  Solve the generalized buckling eigenvalue problem

  Ax = - lambda Bx

  for, x, lambda. A shift and invert strategy similar to
  the one employed above is used. Here instead, A-inner
  products are used and the modified eigenvalue problem is,

  (A + sigma B)^{-1} A x = lambda/(lambda - sigma) x
  mu = lambda/(lambda + sigma)

  The original eigenvalues may be obtained using,

  lambda = - mu * sigma/(1 - mu)
*/
EPBucklingShiftInvert::EPBucklingShiftInvert(TacsScalar _sigma, TACSKsm *_ksm,
                                             TACSMat *_inner) {
  sigma = _sigma;
  ksm = _ksm;
  ksm->incref();
  inner = _inner;
  inner->incref();
  temp = inner->createVec();
  temp->incref();
}

EPBucklingShiftInvert::~EPBucklingShiftInvert() {
  ksm->decref();
  inner->decref();
  temp->decref();
}

void EPBucklingShiftInvert::setSigma(TacsScalar _sigma) { sigma = _sigma; }

TACSVec *EPBucklingShiftInvert::createVec() { return ksm->createVec(); }

// Compute y = ( A - sigma B )^{-1} *  x
void EPBucklingShiftInvert::mult(TACSVec *x, TACSVec *y) {
  inner->mult(x, temp);
  return ksm->solve(temp, y);
}

// Compute <x,y> = x^{T} inner y
TacsScalar EPBucklingShiftInvert::dot(TACSVec *x, TACSVec *y) {
  inner->mult(y, temp);
  return temp->dot(x);
}

// Compute || B * x || - this is used to compute the eigenvalue error
TacsScalar EPBucklingShiftInvert::errorNorm(TACSVec *x) {
  inner->mult(x, temp);
  return temp->norm();
}

TacsScalar EPBucklingShiftInvert::convertEigenvalue(TacsScalar value) {
  return -(value * sigma) / (1.0 - value);
}

/*
  Compute the eigenvalues and eigenvectors of a symmetric tridiagonal
  matrix.

  input:
  n:        the order of the matrix
  diag:     the diagonal entries
  upper:    the upper/lower entries of the matrix

  output:
  eigs:     the eigenvalues computed using LAPACK
  eigvecs:  the eigenvectors computed using LAPACK
*/
static void ComputeEigsTriDiag(int n, TacsScalar *_diag, TacsScalar *_upper,
                               TacsScalar *_eigs, TacsScalar *_eigvecs) {
  // The input arguments required for LAPACK
  const char *jobz = "V";
  const char *range = "A";

  // Specify the range of eigenvalues to use - not used in this
  // case. We compute all the eigenvalues in the spectrum
  double vl = 0.0, vu = 0.0;
  int il = 0, iu = 0;

  // The tolerance used in the solution - LAPACK makes this machine
  // precision internally
  double abstol = 0.0;

  // Convert the data to real data in the case of complex arithmetic.
  // Duplicate the data because LAPACK over-writes the matrix on exit.
  double *diag = new double[n];
  double *upper = new double[n - 1];
  for (int i = 0; i < n - 1; i++) {
    diag[i] = TacsRealPart(_diag[i]);
    upper[i] = TacsRealPart(_upper[i]);
  }
  diag[n - 1] = TacsRealPart(_diag[n - 1]);

  // The total number of eigenvalues found on output
  int m = 0;

  // Output and work arrays required by dstevr
  int *isuppz = new int[2 * n];
  int lwork = 20 * n;
  double *work = new double[lwork];
  int liwork = 10 * n;
  int *iwork = new int[liwork];
  int ldz = n;
  int info = -1;

#ifdef TACS_USE_COMPLEX
  /*
    Treat the imaginary part as a perturbation about the given point.
    In the result, replace the imaginary part with the sensitivity to
    this perturbation. This only makes sense in the context of the
    complex step method. It should only be used in this context!
  */
  double *eigs = new double[n];
  double *eigvecs = new double[n * n];

  LAPACKstevr(jobz, range, &n, diag, upper, &vl, &vu, &il, &iu, &abstol, &m,
              eigs, eigvecs, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  if (info != 0) {
    fprintf(stderr, "Error encountered in LAPACK function dstevr\n");
  }

  /*
    The imaginary part of the eigenvalue problem can be determined as
    follows:
    A*x[n*i:(n+1)*i] = eigs[i]*x[n*i:(n+1)*i]

    The sensitivity of the eigenvalues are:
    dA/ds*x + A*dx/ds = deig/ds*x + eig*dx/ds

    Premultiplying by x^{T} yields,
    x^{T}dA/ds*x = x^{T}x*deig/ds
  */

  // Cycle over all the eigenvalues
  for (int i = 0; i < n; i++) {
    double sens = 0.0, dot = 0.0;
    for (int j = 0; j < n; j++) {  // Cycle over matrix entries
      double ans = TacsImagPart(_diag[j]) * eigvecs[n * i + j];
      if (j > 0) {
        ans += TacsImagPart(_upper[j - 1]) * eigvecs[n * i + j - 1];
      }
      if (j < n - 1) {
        ans += TacsImagPart(_upper[j]) * eigvecs[n * i + j + 1];
      }

      dot += eigvecs[n * i + j] * eigvecs[n * i + j];
      sens += ans * eigvecs[n * i + j];
    }

    _eigs[i] = TacsScalar(eigs[i], sens / dot);
  }

  // Discard the eigen sensitivity
  for (int i = 0; i < n * n; i++) {
    _eigvecs[i] = TacsScalar(eigvecs[i], 0.0);
  }

  delete[] eigs;
  delete[] eigvecs;
#else
  // This is the real part of the code that computes the eigenvalues and
  // eigenvectors of the symmetric tridiagonal system
  double *eigs = _eigs;
  double *eigvecs = _eigvecs;

  LAPACKstevr(jobz, range, &n, diag, upper, &vl, &vu, &il, &iu, &abstol, &m,
              eigs, eigvecs, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);

  if (info != 0) {
    fprintf(stderr, "Error encountered in LAPACK function dstevr\n");
  }
#endif  // TACS_USE_COMPLEX

  delete[] isuppz;
  delete[] work;
  delete[] iwork;
  delete[] diag;
  delete[] upper;
}

/*
  Create the symmetric eigenvalue problem solver

  This object uses a Lanczos method to reduce the symmetric eigenvalue
  problem to a small-rank symmetric tridiagonal problem that can be
  solved efficiently using LAPACK routines.

  The symmetric eigenproblem could be derived from either a regular or
  generalized eigenproblem so the inner product is generalized to
  accomodate different options.

  input:
  Op:          the eigenvalue problem operator
  max_iters:   the maximum number of iterations before giving up
  ortho_type:  the type of orthogonalization to use FULL or LOCAL
*/
SEP::SEP(EPOperator *_Op, int _max_iters, OrthoType _ortho_type,
         TACSBcMap *_bcs) {
  // Store the pointer to the eigenproblem operator
  Op = _Op;
  Op->incref();

  // Set the boundary conditions
  bcs = _bcs;
  if (bcs) {
    bcs->incref();
  }

  // Store the information about the subspace vectors
  max_iters = _max_iters;
  ortho_type = _ortho_type;
  Q = new TACSVec *[max_iters + 1];

  // The coefficients of the Lanczos tridiagonal system
  Alpha = new TacsScalar[max_iters];
  Beta = new TacsScalar[max_iters];

  // The eigenvalues and eigenvectors of the tridiagonal system
  eigs = new TacsScalar[max_iters];
  eigvecs = new TacsScalar[max_iters * max_iters];

  // Permutation of the order of the eigenvalues
  perm = new int[max_iters];

  // Default values for convergence tests/sorting of the eigenvalues
  tol = 1e-12;
  spectrum = SMALLEST;
  neigvals = 4;
  niters = -1;

  // Create the vectors required for the Lanczos subspace
  for (int i = 0; i < max_iters + 1; i++) {
    Q[i] = Op->createVec();
    Q[i]->incref();
  }
}

/*
  Deallocate all the information stored in the eigenproblem
*/
SEP::~SEP() {
  Op->decref();
  for (int i = 0; i < max_iters + 1; i++) {
    Q[i]->decref();
  }
  delete[] Q;

  if (bcs) {
    bcs->decref();
  }

  delete[] Alpha;
  delete[] Beta;
  delete[] eigs;
  delete[] eigvecs;
  delete[] perm;
}

/*
  Set the orthogonalization type
*/
void SEP::setOrthoType(enum OrthoType _ortho_type) { ortho_type = _ortho_type; }

/*
  Set the tolerances to use, the desired spectrum, and the number of
  eigenvalues that are requested in the solve
*/
void SEP::setTolerances(double _tol, enum EigenSpectrum _spectrum,
                        int _neigvals) {
  tol = _tol;
  spectrum = _spectrum;
  neigvals = _neigvals;
}

/*
  Reset the underlying operator
*/
void SEP::setOperator(EPOperator *_Op) {
  _Op->incref();
  if (Op) {
    Op->decref();
  }
  Op = _Op;
}

/*
  Solve the eigenvalue problem using the Lanczos method with full or
  local orthogonalization.

  This method uses a spectral shift approach to solve a generalized
  eigenvalue problem with a Lanczos method. The method generates a
  series or orthonormal vectors with respect to a given inner product.
*/
void SEP::solve(KSMPrint *ksm_print, KSMPrint *ksm_file) {
  // Select the initial vector randomly
  Q[0]->setRand();
  if (bcs) {
    Q[0]->applyBCs(bcs);
  }

  // Normalize the first vector
  TacsScalar norm = sqrt(Op->dot(Q[0], Q[0]));
  Q[0]->scale(1.0 / norm);

  if (ortho_type == LOCAL) {
    // Only local orthogonalization is utilized. This code does not
    // orthogonalize the vector against previous vectors.
    int i = 0;
    for (; i < max_iters; i++) {
      // First compute U = A*Q[i] - Q[i-1]*Beta[i]
      Op->mult(Q[i], Q[i + 1]);
      if (bcs) {
        Q[i + 1]->applyBCs(bcs);
      }

      if (i > 0) {
        Q[i + 1]->axpy(-Beta[i - 1], Q[i - 1]);
      }

      // Next, compute the new alpha term
      Alpha[i] = Op->dot(Q[i], Q[i + 1]);
      Q[i + 1]->axpy(-Alpha[i], Q[i]);

      // Compute the new beta term in the Lanczos sequence
      Beta[i] = sqrt(Op->dot(Q[i + 1], Q[i + 1]));
      Q[i + 1]->scale(1.0 / Beta[i]);

      // Check if the desired eigenvalues have converged
      if (checkConverged(Alpha, Beta, i + 1)) {
        niters = i + 1;
        break;
      }
    }

    // Readjust the max number of iterations
    if (i == max_iters) {
      niters = max_iters;
    }
  } else {
    // Perform a full orthogonalization using modified Gram-Schmidt
    int i = 0;
    for (; i < max_iters; i++) {
      // Compute the new vector using the provided operator
      Op->mult(Q[i], Q[i + 1]);
      if (bcs) {
        Q[i + 1]->applyBCs(bcs);
      }

      for (int j = i; j >= 0; j--) {
        TacsScalar h = Op->dot(Q[i + 1], Q[j]);
        Q[i + 1]->axpy(-h, Q[j]);

        // Store the diagonal term (and discard all other terms which
        // will only be non-zero due to numerical issues
        if (j == i) {
          Alpha[i] = h;
        }
      }

      // Evalute the sub-digonal
      Beta[i] = sqrt(Op->dot(Q[i + 1], Q[i + 1]));
      Q[i + 1]->scale(1.0 / Beta[i]);

      // Check if the desired eigenvalues have converged
      if (checkConverged(Alpha, Beta, i + 1)) {
        niters = i + 1;
        break;
      }
    }

    // Readjust the number of iterations if the max iterations
    // has been reached.
    if (i == max_iters) {
      niters = max_iters;
    }
  }

  // Compute the norm of the last vector in the inner product
  TacsScalar er = Op->errorNorm(Q[niters]);

  // Print out a summary of the eigenvalues and errors
  if (ksm_print) {
    char line[256];
    sprintf(line, "%3s %18s %18s %10s\n", " ", "eigenvalue", "shift-invert eig",
            "error");
    ksm_print->print(line);

    for (int i = 0; i < niters; i++) {
      int index = perm[i];
      char line[256];
      sprintf(line, "%3d %18.10e %18.10e %10.3e\n", i,
              TacsRealPart(Op->convertEigenvalue(eigs[index])),
              TacsRealPart(eigs[index]),
              fabs(TacsRealPart(Beta[niters - 1] *
                                eigvecs[index * niters + (niters - 1)] * er)));
      ksm_print->print(line);
    }
  }
  // Print the iteration count to file
  if (ksm_file) {
    char line[256];
    sprintf(line, "%2d\n", niters);
    ksm_file->print(line);
  }
}

/*!
  Extract the n-th eigenvalue from the probelm.
*/
TacsScalar SEP::extractEigenvalue(int n, TacsScalar *error) {
  if (n < 0 || n >= niters) {
    fprintf(stderr, "Eigenvalue out of range\n");
    *error = -1.0;
    return 0.0;
  }

  n = perm[n];

  TacsScalar er = Op->errorNorm(Q[niters]);
  *error = fabs(
      TacsRealPart(Beta[niters - 1] * eigvecs[n * niters + (niters - 1)] * er));

  return Op->convertEigenvalue(eigs[n]);
}

/*!
  Extract the n-th eigenvector from the matrix.

  Compute the eigenvector and store it in ans.  Return the eigenvalue
  and the error computed as error = || A x_n - lambda_n x_n ||_2 via a
  pointer.
*/
TacsScalar SEP::extractEigenvector(int n, TACSVec *ans, TacsScalar *error) {
  if (n < 0 || n >= niters) {
    fprintf(stderr, "Eigenvector out of range\n");
    *error = -1.0;
    ans->zeroEntries();
    return 0.0;
  }

  n = perm[n];

  ans->zeroEntries();
  for (int i = 0; i < niters; i++) {
    ans->axpy(eigvecs[n * niters + i], Q[i]);
  }

  TacsScalar er = Op->errorNorm(Q[niters]);
  *error = fabs(Beta[niters - 1] * eigvecs[n * niters + (niters - 1)] * er);

  return Op->convertEigenvalue(eigs[n]);
}

/*!
  Find the Frobenius norm of the matrix:

  I - Q^{T} Q

  which should be zero if the vectors are completely orthogonal
*/
TacsScalar SEP::checkOrthogonality() {
  // Check to see if all the vectors are orthogonal
  TacsScalar norm = TacsScalar(0.0);

  for (int i = 0; i < niters; i++) {
    for (int j = 0; j < i; j++) {
      TacsScalar aij = Op->dot(Q[i], Q[j]);
      norm += 2.0 * aij * aij;
    }

    TacsScalar aii = 1.0 - Op->dot(Q[i], Q[i]);
    norm += aii * aii;
  }

  return sqrt(norm);
}

/*!
  Print the dot product of the Krylov subspace to the screen
*/
void SEP::printOrthogonality() {
  for (int i = 0; i < niters; i++) {
    for (int j = 0; j < i; j++) {
      TacsScalar aij = Op->dot(Q[i], Q[j]);
      printf("%8.1e ", TacsRealPart(aij));
    }

    TacsScalar aii = Op->dot(Q[i], Q[i]);
    printf("%8.1e \n", TacsRealPart(aii));
  }
}

/*!
  Sort the values based on the actual eigenvalues not their
  transformed values.

  This uses insertion sort. The eigenvalues are left in their same
  ordering - this is required so that the corresponding eigenvectors
  match. The permutation array is altered so that the array indexed by
  p[i] => eigs[ p[i] ] is sorted in ascending order.
*/
void SEP::sortEigenvalues(TacsScalar *values, int neigs, int *p) {
  // The default ordering
  for (int i = 0; i < neigs; i++) {
    p[i] = i;
  }

  // Set the flags based on the sorting criteria: Use absolute value
  // if we only care about the magnitude and sort ascending if we only
  // care about the smallest values
  int use_abs =
      (spectrum == SMALLEST_MAGNITUDE || spectrum == LARGEST_MAGNITUDE);
  int sort_ascending = (spectrum == SMALLEST_MAGNITUDE || spectrum == SMALLEST);

  // Sort the array using insertion sort
  for (int i = 0; i < neigs; i++) {
    // Convert the transformed eigenvalue into the correct
    // range
    TacsScalar eig_new = Op->convertEigenvalue(values[p[i]]);

    // Take the absolute value of the eigenvalue
    if (use_abs) {
      if (TacsRealPart(eig_new) < 0.0) {
        eig_new *= -1.0;
      }
    }

    // Find out where we should place this within the sorted
    // array
    int j = i - 1;
    for (; j >= 0; j--) {
      // Convert the j-th eigenvalue
      TacsScalar eig_j = Op->convertEigenvalue(values[p[j]]);
      if (use_abs) {
        if (TacsRealPart(eig_j) < 0.0) {
          eig_j *= -1.0;
        }
      }

      // If this is the right place, place the new index here
      if (sort_ascending && TacsRealPart(eig_new) > TacsRealPart(eig_j)) {
        break;
      } else if (!sort_ascending &&
                 TacsRealPart(eig_new) < TacsRealPart(eig_j)) {
        break;
      } else {
        // This is not the right place, just move the indices back
        p[j + 1] = p[j];
      }
    }

    // Place the index in the sorted array
    p[j + 1] = i;
  }
}

/*!
  Perform a convergence test.

  Check whether the eigenvalues at the appropriate end of the spectrum
  have converged based on either a relative or absolute tolerance.xx
*/
int SEP::checkConverged(TacsScalar *A, TacsScalar *B, int n) {
  // If we haven't iterated n-times, we can quit immediately
  if (n < neigvals) {
    return 0;
  }

  // Compute the eigenvalues and eigenvectors of the symmetric
  // tridiagonal matrix whose coefficients are stored in A/B
  // which are the diagonal and upper diagonal of the matrix
  TacsScalar beta = B[n - 1];
  ComputeEigsTriDiag(n, A, B, eigs, eigvecs);

  // Find the permutation which sorts the matrix in the desired
  // order
  sortEigenvalues(eigs, n, perm);

  // Check for convergence of each of the desired eigenvalues
  int is_converged = 1;
  for (int k = 0; k < neigvals; k++) {
    // Read off the index
    int index = perm[k];
    TacsScalar er = Op->errorNorm(Q[n]);

    // Read out the predicted error for the eigenvector
    TacsScalar eig_err = fabs(beta * eigvecs[index * n + (n - 1)] * er);
    if (TacsRealPart(eig_err) > tol) {
      is_converged = 0;
      break;
    }
  }

  return is_converged;
}
