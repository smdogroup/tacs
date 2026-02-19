/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_BERNSTEIN_INTERPOLATION_H
#define TACS_BERNSTEIN_INTERPOLATION_H

/*
  Evaluate the Bernstein shape functions at the given parametric point

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  idx:    the index of the interval for u
  knots:  the interpolation knots in parameter space
  work:   a temporary array of size 2*k

  output:
  N:      the values of the shape functions at u
*/
inline void TacsBernsteinShapeFunctions(const int order, const double u,
                                        double *N) {
  double u1 = 0.5 * (1.0 - u);
  double u2 = 0.5 * (u + 1.0);

  N[0] = 1.0;
  for (int j = 1; j < order; j++) {
    double s = 0.0;
    for (int k = 0; k < j; k++) {
      double t = N[k];
      N[k] = s + u1 * t;
      s = u2 * t;
    }
    N[j] = s;
  }
}

/*
  Evaluate the shape functions and the derivative of the shape functions
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
*/
inline void TacsBernsteinShapeFuncDerivative(const int order, const double u,
                                             double *N, double *Nd) {
  double u1 = 0.5 * (1.0 - u);
  double u2 = 0.5 * (u + 1.0);

  // Compute the basis for the order-1 bernstein polynomial
  N[0] = 1.0;
  for (int j = 1; j < order - 1; j++) {
    double s = 0.0;
    for (int k = 0; k < j; k++) {
      double t = N[k];
      N[k] = s + u1 * t;
      s = u2 * t;
    }
    N[j] = s;
  }

  // Add the contributions to the derivative
  for (int j = 0; j < order; j++) {
    Nd[j] = 0.0;
    if (j > 0) {
      Nd[j] += 0.5 * (order - 1) * N[j - 1];
    }
    if (j < order - 1) {
      Nd[j] -= 0.5 * (order - 1) * N[j];
    }
  }

  // Now compute the full order basis
  TacsBernsteinShapeFunctions(order, u, N);
}

/*
  Evaluate the shape functions and their first and second derivatives
  with respect to the parameter coordinate

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
  Nd:     the derivative of the shape functions at u
  Ndd:    the second derivative of the shape functions at u
*/
inline void TacsBernsteinShapeFuncSecondDerivative(const int order,
                                                   const double u, double *N,
                                                   double *Nd, double *Ndd) {
  double u1 = 0.5 * (1.0 - u);
  double u2 = 0.5 * (u + 1.0);

  // Compute the basis for the order-1 bernstein polynomial
  N[0] = 1.0;
  for (int j = 1; j < order - 2; j++) {
    double s = 0.0;
    for (int k = 0; k < j; k++) {
      double t = N[k];
      N[k] = s + u1 * t;
      s = u2 * t;
    }
    N[j] = s;
  }

  // Add the contributions to the derivative
  for (int j = 0; j < order - 1; j++) {
    Nd[j] = 0.0;
    if (j > 0) {
      Nd[j] += 0.5 * (order - 2) * N[j - 1];
    }
    if (j < order - 2) {
      Nd[j] -= 0.5 * (order - 2) * N[j];
    }
  }

  for (int j = 0; j < order; j++) {
    Ndd[j] = 0.0;
    if (j > 0) {
      Ndd[j] += 0.5 * (order - 1) * Nd[j - 1];
    }
    if (j < order - 1) {
      Ndd[j] -= 0.5 * (order - 1) * Nd[j];
    }
  }

  TacsBernsteinShapeFuncDerivative(order, u, N, Nd);
}

#endif  // TACS_BERNSTEIN_INTERPOLATION_H
