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

#include "FElibrary.h"

/*
  Various useful helper-functions for finite-element calculations
*/

TACS_BEGIN_NAMESPACE(FElibrary)

/*
  Test the derivative of the b-spline basis

  input:
  k: the order of the basis to test
*/
void bspline_basis_test(int k) {
  const int max_k = 8;

  if (k > max_k) {
    k = max_k;
  }

  double dh = 1e-6;
  double work[2 * max_k + max_k * max_k];
  double Na_forward[max_k * (max_k + 1)], Na_reverse[max_k * (max_k + 1)];
  double Na[max_k * (max_k + 1)];
  int ideriv = k;

  int n = 3 * max_k - k;

  double Tu[3 * max_k];
  for (int i = 0; i < k; i++) {
    Tu[i] = 0.0;
    Tu[n + i] = 1.0;
  }

  // n - k + 2 points
  for (int i = 0; i < (n - k + 2); i++) {
    double v = (1.0 * i) / (n - k + 1);
    Tu[i + k - 1] = v;
  }

  double u = 0.8425;
  int idx = bspline_interval(u, Tu, n, k);

  printf("Testing basis functions of order: %d\n", k);

  bspline_basis(Na_forward, idx, u, Tu, k, work);
  for (int i = 0; i < k; i++) {
    printf("N[%d] = %15.5e\n", i, Na_forward[i]);
  }

  // Compute the derivatives of the b-spline basis
  bspline_basis_derivative(Na, idx, u, ideriv, Tu, k, work);

  bspline_basis(Na_forward, idx, u + dh, Tu, k, work);
  bspline_basis(Na_reverse, idx, u - dh, Tu, k, work);
  for (int i = 0; i < k; i++) {
    double fd = 0.5 * (Na_forward[i] - Na_reverse[i]) / dh;
    printf("N^(%d)_[%d,%d]: %15.5e  FD: %15.5e\n", 1, k - 1 - i, k - 1,
           Na[i + k], fd);
  }

  for (int j = 1; j < ideriv; j++) {
    bspline_basis_derivative(Na_forward, idx, u + dh, ideriv, Tu, k, work);
    bspline_basis_derivative(Na_reverse, idx, u - dh, ideriv, Tu, k, work);

    for (int i = 0; i < k; i++) {
      double fd = 0.5 * (Na_forward[i + j * k] - Na_reverse[i + j * k]) / dh;
      printf("N^(%d)_[%d,%d]: %15.5e  FD: %15.5e\n", j + 1, k - 1 - i, k - 1,
             Na[i + (j + 1) * k], fd);
    }
  }
}

/*
  Find the interval for the computing the basis function

  input:
  u: the parametric location
  T: the knot locations
  n: the number of control points
  k: the order of the b-spline
*/
int bspline_interval(double u, const double *T, int n, int k) {
  if (u >= T[n]) {
    return n - 1;
  } else if (u < T[k - 1]) {
    return k - 1;
  } else {
    // Use a binary search to locate the interval required
    int low = k - 1;
    int high = n;
    int mid = low + (high - low) / 2;

    while (u < T[mid] || u >= T[mid + 1]) {
      if (u < T[mid]) {
        high = mid;
      } else {
        low = mid;
      }

      mid = low + (high - low) / 2;
    }

    return mid;
  }
}

/*
  Evaluate the basis functions and optionally their derivatives

  output:
  N:  the basis functions - an array of size k

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1])
*/
void bspline_basis(double *N, const int idx, const double u, const double *Tu,
                   const int ku, double *work) {
  N[0] = 1.0;

  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - Tu[i+1 - j]
  // and right[j] = Tu[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];

  for (int j = 1; j < ku; j++) {
    left[j] = u - Tu[idx + 1 - j];
    right[j] = Tu[idx + j] - u;

    N[j] = 0.0;
    for (int i = 0; i < j; i++) {
      double temp = N[i] / (right[i + 1] + left[j - i]);
      N[i] = N[j] + right[i + 1] * temp;
      N[j] = left[j - i] * temp;
    }
  }
}

/*
  Compute the derivative of the b-spline basis

  output:
  N:  the basis functions - an array of size k

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k + k*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1])
*/
void bspline_basis_derivative(double *N, const int idx, const double u,
                              int ideriv, const double *Tu, const int ku,
                              double *work) {
  if (ideriv >= ku) {
    // Set all higher-order derivatives to zero
    memset(&N[ku * ku], 0, (ideriv + 1 - ku) * ku * sizeof(double));

    // set the maximum derivative to ideriv
    ideriv = ku - 1;
  }

  // Set the pointers for the temporary work arrays
  // Note that left[j] = u - Tu[i+1 - j]
  // and right[j] = Tu[i+j] - u
  double *left = &work[0];
  double *right = &work[ku];

  // ndu is a matrix of total size ku*ku
  // The basis functions are stored in the upper triangular part of the matrix
  // such that N_{idx-j, i} = ndu[i + j*ku] if i >= j is the basis function.
  // The knot differences are stored in the lower portion of the matrix such
  // that ndu[i + j*ku] = u_{idx+i+1} - u_{idx+j-i}
  double *ndu = &work[2 * ku];
  ndu[0] = 1.0;

  // Compute all orders of the basis functions
  for (int j = 1; j < ku; j++) {
    left[j] = u - Tu[idx + 1 - j];
    right[j] = Tu[idx + j] - u;

    double njj = 0.0;
    for (int i = 0; i < j; i++) {
      // Compute the difference Tu[idx+i+1] - Tu[idx+j-i]
      ndu[i + j * ku] = right[i + 1] + left[j - i];
      double temp = ndu[j - 1 + i * ku] / ndu[i + j * ku];

      // Compute the basis function
      ndu[j + i * ku] = njj + right[i + 1] * temp;
      njj = left[j - i] * temp;
    }

    // Store the basis function
    ndu[j * (ku + 1)] = njj;
  }

  // Set the basis functions
  for (int i = 0; i < ku; i++) {
    N[i] = ndu[(ku - 1) + i * ku];
  }

  // Set the temporary arrays for the a-coefficients
  double *a0 = &work[0];
  double *a1 = &work[ku];

  for (int i = 0; i < ku; i++) {
    a0[0] = 1.0;

    for (int k = 1; k <= ideriv; k++) {
      double d = 0.0;

      // Compute the first of the a-terms
      // a_{k,0} = a_{k-1,0}/(u_{i+ku-k} - u_{i})
      if (i >= k) {
        a1[0] = a0[0] / ndu[(i - k) + (ku - k) * ku];
        d = a1[0] * ndu[(ku - k - 1) + (i - k) * ku];
      }

      int jstart = k - i;
      if (i >= k - 1) {
        jstart = 1;
      }

      int jend = ku - i - 1;
      if (i <= ku - k) {
        jend = k - 1;
      }

      for (int j = jstart; j <= jend; j++) {
        // Compute a_{k,j} = (a_{k-1,j} - a_{k-1,j-1})/(u_{i+ku+j-k} - u_{i+j})
        a1[j] = (a0[j] - a0[j - 1]) / ndu[(i - k + j) + (ku - k) * ku];
        d += a1[j] * ndu[(ku - k - 1) + (i - k + j) * ku];
      }

      // Compute the term
      // a_{k,k} = -a_{k-1}/(u_{i+ku} - u_{i+k})
      if (i <= ku - k - 1) {
        a1[k] = -a0[k - 1] / ndu[i + (ku - k) * ku];
        d += a1[k] * ndu[(ku - k - 1) + i * ku];
      }

      // Set the basis function
      N[i + k * ku] = d;

      // Swap the rows of a
      double *t = a0;
      a0 = a1;
      a1 = t;
    }
  }

  // Multiply the basis by the factorial term
  int r = ku - 1;
  for (int k = 1; k <= ideriv; k++) {
    for (int j = 0; j < ku; j++) {
      N[j + k * ku] *= r;
    }
    r *= (ku - 1 - k);
  }
}

/*
  Evaluate the one-dimensional b-spline

  input:
  u:     the parametric location for the spline evaluation
  idu:   the order of the derivative to use
  T:     the knot vector of length n + k
  n:     the number of knots
  k:     the order of the spline to evaluate
  coef:  the spline coefficients
  work:  a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = 3*ku
  otherwise: len = (idu+3)*ku + ku*ku

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline1d(const double u, const int idu, const double *Tu,
                     const int nu, const int ku, const TacsScalar *coef,
                     double *work) {
  double *Nu = work;

  // Compute the knot span
  int intu = bspline_interval(u, Tu, nu, ku);

  // Evaluate the basis functions
  if (idu > 0) {
    bspline_basis_derivative(Nu, intu, u, idu, Tu, ku, &work[(idu + 1) * ku]);
  } else {
    bspline_basis(Nu, intu, u, Tu, ku, &work[ku]);
  }

  // Set the interval to the initial control point
  intu = intu - ku + 1;

  // Evaluate the interpolant
  TacsScalar fval = 0.0;
  for (int i = 0; i < ku; i++) {
    fval += Nu[i + idu * ku] * coef[intu + i];
  }

  return fval;
}

/*
  Evaluate a two-dimensional tensor product b-spline

  input:
  u, v:       the parametric location for the spline evaluation
  idu, idv:   the order of the derivative to use
  Tu, Tv:     the knot vector of length n + k
  nu, nv:     the number of knots
  ku, kv:     the order of the spline to evaluate
  coef:       the spline coefficients
  work:       a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = ku + kv + 2*max(ku, kv)
  otherwise: len = (idu+1)*ku + (idv+1)*kv + max(2*ku + ku**2, 2*kv + kv**2)

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline2d(const double u, const double v, const int idu,
                     const int idv, const double *Tu, const double *Tv,
                     const int nu, const int nv, const int ku, const int kv,
                     const TacsScalar *coef, double *work) {
  // The basis functions
  double *Nu = &work[0];
  double *Nv = &work[(idu + 1) * ku];

  // Compute the knot intervals
  int intu = bspline_interval(u, Tu, nu, ku);
  int intv = bspline_interval(v, Tv, nv, kv);

  // Evaluate the basis functions
  if (idu > 0) {
    bspline_basis_derivative(Nu, intu, u, idu, Tu, ku, &work[(idu + 1) * ku]);
  } else {
    bspline_basis(Nu, intu, u, Tu, ku, &work[ku]);
  }

  // Evaluate the basis functions
  if (idv > 0) {
    bspline_basis_derivative(Nv, intv, v, idv, Tv, kv,
                             &work[(idu + 1) * ku + (idv + 1) * kv]);
  } else {
    bspline_basis(Nv, intv, v, Tv, kv, &work[(idu + 1) * ku + kv]);
  }

  // Set the interval to the initial control point
  intu = intu - ku + 1;
  intv = intv - kv + 1;

  TacsScalar fval = 0.0;
  for (int j = 0; j < kv; j++) {
    for (int i = 0; i < ku; i++) {
      fval += Nu[i + idu * ku] * Nv[j + idv * kv] *
              coef[(intu + i) + (intv + j) * nu];
    }
  }

  return fval;
}

/*
  Evaluate a three-dimensional tensor product b-spline

  input:
  u, v, w:        the parametric location for the spline evaluation
  idu, idv, idw:  the order of the derivative to use
  Tu, Tv, Tw:     the knot vector of length n + k
  nu, nv, nw:     the number of knots
  ku, kv, kw:     the order of the spline to evaluate
  coef:           the spline coefficients
  work:           a working array for temporary storage

  the work array must be of size:
  if idu == 0: len = ku + kv + kw + 2*max(ku, kv, kw)
  otherwise:
  len = (idu+1)*ku + (idv+1)*kv + (idw+1)*kw +
  max(2*ku + ku**2, 2*kv + kv**2, 2*kw + kw**2)

  returns:
  the value of the interpolant (or its derivative) at u
*/
TacsScalar bspline3d(const double u, const double v, const double w,
                     const int idu, const int idv, const int idw,
                     const double *Tu, const double *Tv, const double *Tw,
                     const int nu, const int nv, const int nw, const int ku,
                     const int kv, const int kw, const TacsScalar *coef,
                     double *work) {
  // The basis functions
  double *Nu = &work[0];
  double *Nv = &work[(idu + 1) * ku];
  double *Nw = &work[(idu + 1) * ku + (idv + 1) * kv];

  // Compute the knot intervals
  int intu = bspline_interval(u, Tu, nu, ku);
  int intv = bspline_interval(v, Tv, nv, kv);
  int intw = bspline_interval(w, Tw, nw, kw);

  // Evaluate the basis functions
  if (idu > 0) {
    bspline_basis_derivative(Nu, intu, u, idu, Tu, ku, &work[(idu + 1) * ku]);
  } else {
    bspline_basis(Nu, intu, u, Tu, ku, &work[ku]);
  }

  // Evaluate the basis functions
  if (idv > 0) {
    bspline_basis_derivative(Nv, intv, v, idv, Tv, kv,
                             &work[(idu + 1) * ku + (idv + 1) * kv]);
  } else {
    bspline_basis(Nv, intv, v, Tv, kv, &work[(idu + 1) * ku + kv]);
  }

  // Evaluate the basis functions
  if (idw > 0) {
    bspline_basis_derivative(
        Nw, intw, w, idw, Tw, kw,
        &work[(idu + 1) * ku + (idv + 1) * kv + (idw + 1) * kw]);
  } else {
    bspline_basis(Nw, intw, w, Tw, kw,
                  &work[(idu + 1) * ku + (idv + 1) * kv + (idw + 1) * kw]);
  }

  // Set the interval to the initial control point
  intu = intu - ku + 1;
  intv = intv - kv + 1;
  intw = intw - kw + 1;

  TacsScalar fval = 0.0;
  for (int k = 0; k < kw; k++) {
    for (int j = 0; j < kv; j++) {
      for (int i = 0; i < ku; i++) {
        fval += (Nu[i + idu * ku] * Nv[j + idv * kv] * Nw[k + idw * kw] *
                 coef[(intu + i) + (intv + j) * nu + (intw + k) * nu * nv]);
      }
    }
  }

  return fval;
}

/*
  Compute a cubic Hermite polynomial interpolant

  input:
  a: the parametric location on the interval [-1, 1]

  output:
  N: the shape functions
  Na: the derivative of the shape functions
  Naa: the second derivative of the shape functions
*/
void cubicHP(double N[], double Na[], double Naa[], double a) {
  // Shape functions
  N[0] = 0.5 + (-0.75 + 0.25 * a * a) * a;
  N[1] = 0.25 + (-0.25 + (-0.25 + 0.25 * a) * a) * a;

  N[2] = 0.5 + (0.75 - 0.25 * a * a) * a;
  N[3] = -0.25 + (-0.25 + (0.25 + 0.25 * a) * a) * a;

  // First derivatives
  Na[0] = -0.75 + 0.75 * a * a;
  Na[1] = -0.25 - 0.5 * a + 0.75 * a * a;

  Na[2] = 0.75 - 0.75 * a * a;
  Na[3] = -0.25 + 0.5 * a + 0.75 * a * a;

  // Second derivatives
  Naa[0] = 1.5 * a;
  Naa[1] = -0.5 + 1.5 * a;

  Naa[2] = -1.5 * a;
  Naa[3] = 0.5 + 1.5 * a;
}

/*
  A quintic Hermite polynomial interpolant

  input:
  a: the parametric location on the interval [-1, 1]

  output:
  N: the shape functions
  Na: the derivative of the shape functions
  Naa: the second derivative of the shape functions
*/
void quinticHP(double N[], double Na[], double Naa[], double a) {
  // Shape functions
  N[0] = (1.0 + (-1.25 + (-0.5 + 0.75 * a) * a) * a) * a * a;
  N[1] = 0.25 * ((1.0 + (-1.0 + (-1.0 + a) * a) * a) * a * a);

  N[2] = 1.0 + (-2.0 + a * a) * a * a;
  N[3] = (1.0 + (-2.0 + a * a) * a * a) * a;

  N[4] = (1.0 + (1.25 + (-0.5 - 0.75 * a) * a) * a) * a * a;
  N[5] = 0.25 * ((-1.0 + (-1.0 + (1.0 + a) * a) * a) * a * a);

  // First shape function derivatives
  Na[0] = (2.0 + (-3.75 + (-2.0 + 3.75 * a) * a) * a) * a;
  Na[1] = (0.5 + (-0.75 + (-1.0 + 1.25 * a) * a) * a) * a;

  Na[2] = (-4.0 + 4.0 * a * a) * a;
  Na[3] = (1.0 + (-6.0 + 5.0 * a * a) * a * a);

  Na[4] = (2.0 + (3.75 + (-2.0 - 3.75 * a) * a) * a) * a;
  Na[5] = (-0.5 + (-0.75 + (1.0 + 1.25 * a) * a) * a) * a;

  // Second derivatives
  Naa[0] = 2.0 + (-7.5 + (-6.0 + 15.0 * a) * a) * a;
  Naa[1] = 0.5 + (-1.5 + (-3.0 + 5.0 * a) * a) * a;

  Naa[2] = -4.0 + 12.0 * a * a;
  Naa[3] = (-12.0 + 20.0 * a * a) * a;

  Naa[4] = 2.0 + (7.5 + (-6.0 - 15.0 * a) * a) * a;
  Naa[5] = -0.5 + (-1.5 + (3.0 + 5.0 * a) * a) * a;
}

TACS_END_NAMESPACE
