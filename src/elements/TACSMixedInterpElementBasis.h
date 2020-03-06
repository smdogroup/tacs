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

#ifndef TACS_MIXED_INTERPOLATION_BASIS_H
#define TACS_MIXED_INTERPOLATION_BASIS_H

#include "TACSElementBasis.h"

class TACSMixedInterpElementBasis : public TACSElementBasis {
 public:

  /**
    Get the number of tying point locations

    @return The number of tying points within the element
  */
  virtual int getNumTyingPoints() = 0;

  /**
    Get the number of tying point field values.

    Note that this must match the expected number of tying point values in the element
  */
  virtual int getNumTyingFieldValues() = 0;

  /**
    Get the parametric location of the given tying point

    @param ty The tying point index
    @param pt The parametric location of the tying point
  */
  virtual void getTyingPoint( int ty, double pt[] ) = 0;

  /**
    Evaluate the tying field interpolation at the specified quadrature point

    @param ty The quadrature point index
    @param pt The parametric location of the tying point
    @param qty The quantity values at the tying point locations
    @param Uty The quantity field values
  */
  virtual void getTyingFieldValue( int n, const double pt[],
                                   const TacsScalar tyvars[],
                                   TacsScalar Uty[] );
};

#endif // TACS_MIXED_INTERPOLATION_BASIS_H
