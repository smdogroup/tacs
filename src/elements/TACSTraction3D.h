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

#ifndef TACS_TRACTION_3D_H
#define TACS_TRACTION_3D_H

#include "TACSElement3D.h"

class TACSTraction3D : public TACSElement {
 public:
  TACSTraction3D(int _varsPerNode, int _faceIndex, TACSElementBasis *_basis,
                 const TacsScalar trac[], int tractionCoordinateComponent = 1);
  TACSTraction3D(int _varsPerNode, int _faceIndex, TACSElementBasis *_basis,
                 void (*_getTractionComponents)(int, int, double,
                                                const TacsScalar *,
                                                const TacsScalar *,
                                                TacsScalar *));
  ~TACSTraction3D();

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
  int tractionCoordinateComponent;
  TacsScalar trac[3 * TACSElement3D::MAX_VARS_PER_NODE];
  void (*getTractionComponents)(int, int, double, const TacsScalar *,
                                const TacsScalar *, TacsScalar *);
};

#endif  // TACS_TRACTION_3D_H
