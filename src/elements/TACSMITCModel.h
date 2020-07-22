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

#ifndef TACS_MITC_MODEL_H
#define TACS_MITC_MODEL_H

#include "TACSElementModel.h"

/**
  TACSMITCModel defines a physical model class independent of a
  finite element basis
*/
class TACSMITCModel : public TACSElementModel {
 public:
  /**
    Get the number of tying field values required by this element

    @return The number of tying field values
  */
  virtual int getNumTyingFields() = 0;

  /**
    Evaluate the given field value at the specified tying point.

    @param field The tying field index
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

  /**
    Evaluate the derivative of the tying point quantity with respect to the state
    variables at the specified tying point.

    This code comptues the

    @param field The tying field index
    @param ty The tying point index
    @param pt The parametric location of the tying point
    @param dfdqty The input derivative seed w.r.t. the tying point
    @param Xd The derivative of the spatial coordinates w.r.t. computational parameters
    @param U The values of the field variables at the tying point
    @param Ud The derivative of the field variables w.r.t. computational parameters
    @param DU The derivative w.r.t. the tying field variables
    @param DUd The derivative w.r.t. the spatial derivatives of the field
  */
  virtual void evalTyingQuantitySVSens( int field, int ty, const double pt[],
                                        const TacsScalar dfdqty,
                                        const TacsScalar Xd[],
                                        const TacsScalar U[], const TacsScalar Ud[],
                                        TacsScalar DU[], TacsScalar DUd[] ) = 0;
};

#endif // TACS_MITC_MODEL_H
