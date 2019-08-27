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

#ifndef TACS_ELEMENT_2D_H
#define TACS_ELEMENT_2D_H

#include "TACSElement.h"
#include "TACSElementModel.h"
#include "TACSElementBasis.h"

class TACSElement2D : public TACSElement {
 public:
  static const int MAX_VARS_PER_NODE = 8;

  TACSElement2D( TACSElementModel *_model, TACSElementBasis *_basis );
  ~TACSElement2D();

  // Get the layout properties of the element
  int getVarsPerNode();
  int getNumNodes();
  ElementLayout getLayoutType();

  /**
    Add the residual to the provided vector

    @param time The simulation time
    @param Xpts The element node locations
    @param vars The element state variable values
    @param dvars The time derivative of the element state variables
    @param ddvars The second time derivative of the element state variables
    @param res The element residual output
  */
  void addResidual( int elemIndex, double time, const TacsScalar *Xpts,
                    const TacsScalar *vars, const TacsScalar *dvars,
                    const TacsScalar *ddvars, TacsScalar *res );

  /**
    Add the residual and Jacobians to the arrays

    @param time The simulation time
    @param Xpts The element node locations
    @param vars The element state variable values
    @param dvars The time derivative of the element state variables
    @param ddvars The second time derivative of the element state variables
    @param res The element residual output
    @param mat The element Jacobian output
  */
  void addJacobian( int elemIndex, double time,
                    double alpha, double beta, double gamma,
                    const TacsScalar *Xpts, const TacsScalar *vars,
                    const TacsScalar *dvars, const TacsScalar *ddvars,
                    TacsScalar *res, TacsScalar *mat );

  // Functions for the adjoint
  /*
  void addAdjResProduct( double time, double scale, const TacsScalar psi[],
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        int dvLen, TacsScalar dvSens[] ){}
  void addAdjResXptProduct( double time, double scale, const TacsScalar psi[],
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            TacsScalar fXptSens[] ){}
  */
 private:
  TACSElementModel *model;
  TACSElementBasis *basis;
};

#endif // TACS_ELEMENT_2D_H
