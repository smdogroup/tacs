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

class TACSMixedInterpElementBasis : public TACSElementBasis {
public:
  //! Set the maximum number of tying field values
  static const int MAX_NUM_TYING_FIELDS = 5;

  //! Set the maximum number of points in any one tying field
  static const int MAX_NUM_TYING_POINTS = 16;

  //! Set the total number of tying field points
  static const int MAX_TOTAL_TYING_POINTS = MAX_NUM_TYING_FIELDS*MAX_NUM_TYING_POINTS;

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
  virtual int getNumTyingPoints( int field ) = 0;

  /**
     Get the parametric location of the given tying point

     @param field The tying field index
     @param ty The tying point index
     @param pt The parametric location of the tying point
  */
  virtual void getTyingPoint( int field, int ty, double pt[] ) = 0;

  /**
     Compute the tying field interpolation.

     Note that the number of shape functions is equal to the quantity
     returned by getNumTyingPoints(field)

     @param field The tying field index
     @param n The quadrature point index (Note: not the tying point index)
     @param pt The parametric location of the quadrature point
     @param N The shape functions for the tying field
  */
  virtual void computeTyingFieldBasis( int field,
                                       int n,
                                       const double pt[],
                                       double N[] ) = 0;

  /**
     Evaluate the tying field interpolation at the specified quadrature point

     @param n The element quadrature point index
     @param pt The parametric location of the tying point
     @param qty The quantity values at the tying point locations
     @param U The quantity field values
  */
  void getTyingFieldValues( int n, const double pt[],
                            const TacsScalar tyvars[],
                            TacsScalar U[] );

  /**
     Add the weak residual to the

     @param n The element quadrature point index
     @param pt The parametric location of the tying point
     @param weight The weight applied to the residual
     @param vars_per_node Number of variables at each node
     @param DUty The residual coefficient for each field
     @param rty The residual of the tying strain
  */
  void addMixedWeakFormResidual( int n,
                                 const double pt[],
                                 TacsScalar weight,
                                 const int vars_per_node,
                                 const TacsScalar DUty[],
                                 TacsScalar rty[] );

private:

  static void addMixedWeakFormResidual( const int num_ty_fields,
                                        const int num_ty_points[],
                                        const double Nty[],
                                        TacsScalar weight,
                                        const int vars_per_node,
                                        const TacsScalar DUty[],
                                        TacsScalar rty[] );
};

#endif // TACS_MIXED_INTERPOLATION_BASIS_H
