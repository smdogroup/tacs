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

#ifndef TACS_MITC_BASIS_H
#define TACS_MITC_BASIS_H

#include "TACSElementBasis.h"

class TACSMITCBasis : public TACSElementBasis {
 public:
  //! Set the maximum number of tying field values
  static const int MAX_NUM_TYING_FIELDS = 5;

  //! Set the maximum number of points in any one tying field
  static const int MAX_NUM_TYING_POINTS = 16;

  //! Set the total number of tying field points
  static const int MAX_TOTAL_TYING_POINTS =
      MAX_NUM_TYING_FIELDS * MAX_NUM_TYING_POINTS;

  /**
    Get the number of tying fields associated with the element

    @return The number of tying fields defined by this basis
  */
  virtual int getNumTyingFields() = 0;

  /**
    Get the number of tying point locations for the specified field

    @param The tying field index
    @return The number of tying points within the element
  */
  virtual int getNumTyingPoints(int field) = 0;

  /**
    Get the parametric location of the given tying point

    @param field The tying field index
    @param ty The tying point index
    @param pt The parametric location of the tying point
  */
  virtual void getTyingPoint(int field, int ty, double pt[]) = 0;

  /**
    Get the field gradient of the solution fields at the specified tying point

    This computes the gradient in the computational space at the specified
    tying point. Note that the values of the field variables are computed but
    not their time derivatives.

    @param field The tying field index
    @param ty The tying field point index
    @param pt The parametric location of the tying point
    @param Xpts The node locations
    @param vars_per_node The number of variables at each node
    @param vars The node locations
    @param Xd The derivative of the position w.r.t. computational parameters
    @param U The values of the field variables
    @param Ud The derivative of the field variables w.r.t. computational
    parameters
  */
  void getTyingPointGradient(int field, int ty, const double pt[],
                             const TacsScalar Xpts[], const int vars_per_node,
                             const TacsScalar vars[], TacsScalar Xd[],
                             TacsScalar U[], TacsScalar Ud[]);

  /**
    Interpolate the tying quantities in the parametric space to obtain
    the values at the specified quadrature point.

    @param n The element quadrature point index
    @param pt The parametric location of the quadrature point
    @param qty The quantity values at the tying point locations
    @param Uty The tying field values at the quadrature point
  */
  virtual void interpTyingField(int n, const double pt[],
                                const TacsScalar qty[], TacsScalar Uty[]);

  /**
    Add the transpose of the interpolation to the tying field values

    @param n The element quadrature point index
    @param pt The parametric location of the quadrature point
    @param Uty The derivatives w.r.t. the tying field values
    @param rty The residual at the tying points
  */
  virtual void addInterpTyingFieldTranspose(int n, const double pt[],
                                            const TacsScalar Uty[],
                                            TacsScalar rty[]);

  /**
    Add the weak form of the governing equations to the residual

    @param n The quadrautre point index
    @param pt The quadrature point value
    @param weight The quadrature weight
    @param J The Jacobian coordinate transformation
    @param vars_per_node The number of variables per node
    @param DUt The coefficients of the temporal part of the weak form
    @param DUx The coefficients of the spatial part of the weak form
    @param res The residual
    @param resty The residual for the tying points
  */
  void addWeakResidual(int n, const double pt[], TacsScalar weight,
                       const TacsScalar J[], const int vars_per_node,
                       TacsScalar DUt[], TacsScalar DUx[], TacsScalar res[],
                       TacsScalar resty[]);

  /*
    Add the residual from a tying point

  */
  void addTyingResidual(int field, int ty, const double pt[],
                        const int vars_per_node, const TacsScalar DU[],
                        const TacsScalar DUd[], TacsScalar res[]);

  /**
    Compute the tying field interpolation.

    Note that the number of shape functions is equal to the quantity
    returned by getNumTyingPoints(field)

    @param field The tying field index
    @param pt The parametric location of the quadrature point
    @param N The shape functions for the tying field
  */
  virtual void computeTyingBasis(int field, const double pt[], double N[]) = 0;
};

#endif  // TACS_MITC_BASIS_H
