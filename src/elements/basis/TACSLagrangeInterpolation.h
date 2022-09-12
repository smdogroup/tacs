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

#ifndef TACS_LAGRANGE_INTERPOLATION_H
#define TACS_LAGRANGE_INTERPOLATION_H

static const double TacsGaussLobattoPoints2[] = {-1.0, 1.0};

static const double TacsGaussLobattoPoints3[] = {-1.0, 0.0, 1.0};

static const double TacsGaussLobattoPoints4[] = {-1.0, -0.5, 0.5, 1.0};

static const double TacsGaussLobattoPoints5[] = {-1.0, -0.7071067811865475, 0.0,
                                                 0.7071067811865475, 1.0};

static const double TacsGaussLobattoPoints6[] = {-1.0,
                                                 -0.8090169943749475,
                                                 -0.30901699437494745,
                                                 0.30901699437494745,
                                                 0.8090169943749475,
                                                 1.0};

/*
  Evaluate the shape functions at the given parametric point

  input:
  order:  the order of the polynomial and number of knots
  u:      the parametric coordinate
  knots:  the interpolation knots in parameter space

  output:
  N:      the values of the shape functions at u
*/
inline void TacsLagrangeShapeFunctions(const int order, const double u,
                                       const double *knots, double *N) {
  // Loop over the shape functions
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (u - knots[j]) * d;
      }
    }
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
inline void TacsLagrangeShapeFuncDerivative(const int order, const double u,
                                            const double *knots, double *N,
                                            double *Nd) {
  // Loop over the shape function knot locations
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    Nd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double d = 1.0 / (knots[i] - knots[j]);
        N[i] *= (u - knots[j]) * d;

        // Now add up the contribution to the derivative
        for (int k = 0; k < order; k++) {
          if (k != i && k != j) {
            d *= (u - knots[k]) / (knots[i] - knots[k]);
          }
        }

        // Add the derivative contribution
        Nd[i] += d;
      }
    }
  }
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
inline void TacsLagrangeShapeFuncSecondDerivative(const int order,
                                                  const double u,
                                                  const double *knots,
                                                  double *N, double *Nd,
                                                  double *Ndd) {
  // Loop over the shape function control points
  for (int i = 0; i < order; i++) {
    N[i] = 1.0;
    Nd[i] = 0.0;
    Ndd[i] = 0.0;

    // Loop over each point again, except for the current control
    // point, adding the contribution to the shape function
    for (int j = 0; j < order; j++) {
      if (i != j) {
        double tj = 1.0 / (knots[i] - knots[j]);
        double dj = tj;
        N[i] = N[i] * (u - knots[j]) * dj;

        // Loop over all the knots again to determine the
        // contribution to the derivative of the shape function
        for (int k = 0; k < order; k++) {
          if (k != i && k != j) {
            double dk = 1.0 / (knots[i] - knots[k]);
            dj *= (u - knots[k]) * dk;
            dk *= tj;

            for (int m = 0; m < order; m++) {
              if (m != i && m != j && m != k) {
                dk *= (u - knots[m]) / (knots[i] - knots[m]);
              }
            }
            Ndd[i] += dk;
          }
        }
        Nd[i] += dj;
      }
    }
  }
}

#endif  // TACS_LAGRANGE_INTERPOLATION_H
