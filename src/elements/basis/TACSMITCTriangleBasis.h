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

#ifndef TACS_MIXED_INTERPOLATION_BASIS_H
#define TACS_MIXED_INTERPOLATION_BASIS_H

#include "TACSElementBasis.h"

class TACSMixedQuadraticPlateTriangleBasis : public TACSElementBasis {
 public:
  /**
    Get the number of tying fields associated with the element

    @return The number of tying fields defined by this basis
  */
  virtual int getNumTyingFields() { return 2; }

  /**
    Get the number of tying point locations for the specified field

    @param The tying field index
    @return The number of tying points within the element
  */
  virtual int getNumTyingPoints(int field) {
    if (field == 0 || field == 1) {
      return 3;
    }
    return 0;
  }

  /**
    Get the parametric location of the given tying point

    @param field The tying field index
    @param ty The tying point index
    @param pt The parametric location of the tying point
  */
  virtual void getTyingPoint(int field, int ty, double pt[]) {
    if (field == 0) {
      if (ty == 0) {
      } else if (ty == 1) {
      } else if (ty == 2) {
      }
    } else if (field == 1) {
    }
  }

  /**
    Evaluate the tying field interpolation at the specified quadrature point

    @param n The element quadrature point index
    @param pt The parametric location of the tying point
    @param qty The quantity values at the tying point locations
    @param Uty The quantity field values
  */
  void getTyingFieldValues(int n, const double pt[], const TacsScalar qty[],
                           TacsScalar Uty[]) {
    // Evaluate the tying field interpolation
    double N23[3], N13[3];

    Uty[0] = N23[0] * qty[0] + N23[1] * qty[1] + N23[2] * qty[2];
    Uty[1] = N13[0] * qty[3] + N13[1] * qty[4] + N13[2] * qty[5];
  }

  /*
    Add contributions to the derivative of each tying field quantity
  */
  void addTyingFieldSVSens(int n, const double pt[], const TacsScalar scale,
                           const TacsScalar dUty[], TacsScalar dqty[]) {
    TacsScalar s1 = scale * dUty[0];
    dqty[0] += s1 * N23[0];
    dqty[1] += s1 * N23[1];
    dqty[2] += s1 * N23[2];

    TacsScalar s2 = scale * dUty[1];
    dqty[0] += s1 * N13[0];
    dqty[1] += s1 * N13[1];
    dqty[2] += s1 * N13[2];
  }
};

#endif  // TACS_MIXED_INTERPOLATION_BASIS_H
