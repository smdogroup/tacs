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

#ifndef TACS_PRESSURE_2D_H
#define TACS_PRESSURE_2D_H

#include "TACSElement2D.h"

class TACSPressure2D : public TACSElement {
 public:
  TACSPressure2D(int _varsPerNode, int _faceIndex, TACSElementBasis *_basis,
                 TacsScalar p);
  ~TACSPressure2D();

  // Get the layout properties of the element
  const char *getObjectName();
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

 private:
  int varsPerNode, faceIndex;
  TACSElementBasis *basis;
  TacsScalar p;
};

#endif  // TACS_PRESSURE_2D_H
