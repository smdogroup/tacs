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

#ifndef TACS_MIXED_INTERP_ELEMENT_MODEL_H
#define TACS_MIXED_INTERP_ELEMENT_MODEL_H

#include "TACSElementModel.h"

/**
  TACSMixedInterpElementModel defines a physical model class independent of a
  finite element basis
*/
class TACSMixedInterpElementModel : public TACSElementModel {
 public:
  /**
    Get the number of tying field values required by this element

    @return The number of tying field values
  */
  virtual int getNumTyingFields() = 0;

  /**
    Evaluate the given field value at the specified tying point.

    @param field The field index
    @param ty The tying point index
    @param pt The parametric location of the tying point
    @param Xd The derivative of the spatial coordinates w.r.t. computational parameters
    @param U The values of the displacements
    @param Ud The derivative of the displacements w.r.t. computational parameters
    @return The tying field quantity at this point
   */
  virtual TacsScalar evalTyingQuantity( int field, int ty,
                                        const double pt[],
                                        const TacsScalar Xd[],
                                        const TacsScalar U[],
                                        const TacsScalar Ud[] ) = 0;

};

#endif // TACS_MIXED_INTERP_ELEMENT_MODEL_H
