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

#ifndef TACS_FE_LIBRARY_H
#define TACS_FE_LIBRARY_H

/*!
  FElibrary.h contains many important functions and data that are
  repeatedly used in the formation of various finite element stiffness
  matricies. The intent of this code is to provide functions for
  calculating shape functions, Gauss points and functions for
  calculating the Jacobians for element stiffness matricies.
*/

#include <math.h>
#include <stdlib.h>

#include "TACSObject.h"

/*
  The following data defines the quadrature rules that can be used in
  the elements. The default is to use tensor-product Gauss quadrature
  rules, however, more sophisticated methods can be used.

  Currently, the options are Gauss quadrature or Lobatto (or
  Gauss-Lobatto) quadrature schemes that include the end points of the
  interval.
*/
enum QuadratureType { GAUSS_QUADRATURE, LOBATTO_QUADRATURE };

TACS_BEGIN_NAMESPACE(FElibrary)

/*!
  Solve the quadratic equation and return the positive and negative
  roots.

  The code returns the roots of the equation:

  a*x^2 + b*x + c = 0

  This code avoids truncation error by testing the sign of the
  b coefficient and using the corresponding expression with the
  least susceptibility to truncation error.
*/
template <class ScalarType>
void solveQERoots(ScalarType* r1, ScalarType* r2, ScalarType a, ScalarType b,
                  ScalarType c) {
  ScalarType discrim = b * b - 4.0 * a * c;
  if (TacsRealPart(discrim) < 0.0) {
    *r1 = *r2 = 0.0;
    return;
  }

  if (TacsRealPart(a) == 0.0) {
    if (TacsRealPart(b) == 0.0) {
      *r1 = *r2 = 0.0;
      return;
    }

    // Solve b*x + c = 0
    *r1 = -c / b;
    *r2 = 0.0;
    return;
  }

  // Depending on the sign of b, use different expressions to
  // avoid truncation error
  discrim = sqrt(discrim);

  if (TacsRealPart(b) > 0.0) {
    *r1 = -(b + discrim) / (2.0 * a);
    *r2 = c / ((*r1) * a);
  } else {  // b < 0.0
    *r1 = -(b - discrim) / (2.0 * a);
    *r2 = c / ((*r1) * a);
  }
}

/*
  Test the derivative of the b-spline basis

  input:
  k: the order of the basis to test
*/
void bspline_basis_test(int k);

/*
  Find the interval for the computing the basis function

  input:
  u: the parametric location
  T: the knot locations
  n: the number of control points
  k: the order of the b-spline
*/
int bspline_interval(double u, const double* T, int n, int k);

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
void bspline_basis(double* N, const int idx, const double u, const double* Tu,
                   const int ku, double* work);

/*
  Compute the derivative of the b-spline basis

  output:
  N:  the basis functions - an array of size k*(ideriv+1)

  input:
  idx:  the index of the interval for u
  u:    the parametric position
  Tu:   the knot vector
  ku:   the order of the basis functions
  work: a temporary array of size 2*k + k*k

  u is on the idx-th knot span such that u is in the interval
  u \in [Tu[idx], Tu[idx+1])
*/
void bspline_basis_derivative(double* N, const int idx, const double u,
                              int ideriv, const double* Tu, const int ku,
                              double* work);

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
TacsScalar bspline1d(const double u, const int idu, const double* Tu,
                     const int nu, const int ku, const TacsScalar* coef,
                     double* work);

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
                     const int idv, const double* Tu, const double* Tv,
                     const int nu, const int nv, const int ku, const int kv,
                     const TacsScalar* coef, double* work);

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
                     const double* Tu, const double* Tv, const double* Tw,
                     const int nu, const int nv, const int nw, const int ku,
                     const int kv, const int kw, const TacsScalar* coef,
                     double* work);

/*
  C1 functions for one dimensional problems
*/
void cubicHP(double N[], double Na[], double Naa[], double a);
void quinticHP(double N[], double Na[], double Naa[], double a);

TACS_END_NAMESPACE

#endif
