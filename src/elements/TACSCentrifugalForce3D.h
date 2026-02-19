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

#ifndef TACS_CENTRIFUGAL_FORCE_3D_H
#define TACS_CENTRIFUGAL_FORCE_3D_H

#include "TACSConstitutive.h"
#include "TACSElement3D.h"

class TACSCentrifugalForce3D : public TACSElement {
 public:
  TACSCentrifugalForce3D(int _varsPerNode, TACSConstitutive *_con,
                         TACSElementBasis *_basis, const TacsScalar _omegaVec[],
                         const TacsScalar _rotCenter[]);
  ~TACSCentrifugalForce3D();

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
      Retrieve the global design variable numbers associated with this element
    */
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  /**
    Set the element design variables from the design vector
  */
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  /**
    Get the element design variables values
  */
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  /**
    Get the lower and upper bounds for the design variable values
  */
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

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
  int varsPerNode;
  TACSConstitutive *con;
  TACSElementBasis *basis;
  TacsScalar omegaVec[3], rotCenter[3];
};

#endif  // TACS_CENTRIFUGAL_FORCE_3D_H
