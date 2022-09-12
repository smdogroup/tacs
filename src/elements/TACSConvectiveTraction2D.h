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

#ifndef TACS_CONVECTIVE_TRACTION_2D_H
#define TACS_CONVECTIVE_TRACTION_2D_H

#include "TACSElement2D.h"

class TACSConvectiveTraction2D : public TACSElement {
 public:
  TACSConvectiveTraction2D(int _varsPerNode, int _faceIndex, int _fieldIndex,
                           TacsScalar alpha, TacsScalar _refValue,
                           TACSElementBasis *_basis);
  ~TACSConvectiveTraction2D();

  // Get the layout properties of the element
  int getVarsPerNode();
  int getNumNodes();
  ElementLayout getLayoutType();
  TACSElementBasis *getElementBasis();
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double tangent[]);

  /**
    Add the residual to the provided vector
  */
  void addResidual(int elemIndex, double time, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res);

  /**
    Add the residual and Jacobians to the arrays
  */
  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res, TacsScalar *mat);

  /**
    Add the derivative of the product of the adjoint variables w.r.t.
    the material design variables
  */
  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dvSens[]);

  /**
    Add the derivative of the product of the adjoint variables and the
    residuals with respect to the node locations
  */
  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

 private:
  int varsPerNode, faceIndex, fieldIndex;
  TacsScalar alpha, refValue;
  TACSElementBasis *basis;
};

#endif  // TACS_CONVECTIVE_TRACTION_2D_H
