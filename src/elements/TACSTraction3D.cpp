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

#include "TACSTraction3D.h"

#include "TACSElementAlgebra.h"

TACSTraction3D::TACSTraction3D(int _varsPerNode, int _faceIndex,
                               TACSElementBasis *_basis,
                               const TacsScalar _trac[],
                               int _tractionCoordinateComponent) {
  varsPerNode = _varsPerNode;
  faceIndex = _faceIndex;
  basis = _basis;
  basis->incref();
  tractionCoordinateComponent = _tractionCoordinateComponent;
  getTractionComponents = NULL;
  if (tractionCoordinateComponent) {
    memcpy(trac, _trac, varsPerNode * sizeof(TacsScalar));
  } else {
    memcpy(trac, _trac, 3 * varsPerNode * sizeof(TacsScalar));
  }
}

TACSTraction3D::TACSTraction3D(
    int _varsPerNode, int _faceIndex, TACSElementBasis *_basis,
    void (*_getTractionComponents)(int, int, double, const TacsScalar *,
                                   const TacsScalar *, TacsScalar *)) {
  varsPerNode = _varsPerNode;
  faceIndex = _faceIndex;
  basis = _basis;
  basis->incref();
  tractionCoordinateComponent = 0;
  getTractionComponents = _getTractionComponents;
}

TACSTraction3D::~TACSTraction3D() { basis->decref(); }

const char *TACSTraction3D::getObjectName() { return "TACSTraction3D"; }

// Get the layout properties of the element
int TACSTraction3D::getVarsPerNode() { return varsPerNode; }

int TACSTraction3D::getNumNodes() { return basis->getNumNodes(); }

ElementLayout TACSTraction3D::getLayoutType() { return basis->getLayoutType(); }

TACSElementBasis *TACSTraction3D::getElementBasis() { return basis; }

int TACSTraction3D::getNumQuadraturePoints() {
  return basis->getNumQuadraturePoints();
}

double TACSTraction3D::getQuadratureWeight(int n) {
  return basis->getQuadratureWeight(n);
}

double TACSTraction3D::getQuadraturePoint(int n, double pt[]) {
  return basis->getQuadraturePoint(n, pt);
}

int TACSTraction3D::getNumElementFaces() { return basis->getNumElementFaces(); }

int TACSTraction3D::getNumFaceQuadraturePoints(int face) {
  return basis->getNumFaceQuadraturePoints(face);
}

double TACSTraction3D::getFaceQuadraturePoint(int face, int n, double pt[],
                                              double tangent[]) {
  return basis->getFaceQuadraturePoint(face, n, pt, tangent);
}

/*
  Add the residual to the provided vector
*/
void TACSTraction3D::addResidual(int elemIndex, double time,
                                 const TacsScalar *Xpts, const TacsScalar *vars,
                                 const TacsScalar *dvars,
                                 const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumFaceQuadraturePoints(faceIndex);

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3], tangent[6];
    double weight = basis->getFaceQuadraturePoint(faceIndex, n, pt, tangent);

    // Get the face normal
    TacsScalar X[3], Xd[9], normal[3];
    TacsScalar area = basis->getFaceNormal(faceIndex, n, Xpts, X, Xd, normal);

    // Compute the inverse of the transformation
    TacsScalar J[9];
    inv3x3(Xd, J);

    // Multiply the weight by the quadrature point
    area *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement3D::MAX_VARS_PER_NODE];
    TacsScalar DUx[3 * TACSElement3D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 3 * varsPerNode * sizeof(TacsScalar));

    if (getTractionComponents) {
      getTractionComponents(elemIndex, time, faceIndex, X, normal, trac);

      // Set the coefficients for the traction
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -trac[k];
      }
    } else if (tractionCoordinateComponent) {
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -trac[k];
      }
    } else {
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -vec3Dot(&trac[3 * k], normal);
      }
    }

    // Add the weak form of the residual at this point
    basis->addWeakResidual(n, pt, area, J, varsPerNode, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSTraction3D::addJacobian(int elemIndex, double time, TacsScalar alpha,
                                 TacsScalar beta, TacsScalar gamma,
                                 const TacsScalar *Xpts, const TacsScalar *vars,
                                 const TacsScalar *dvars,
                                 const TacsScalar *ddvars, TacsScalar *res,
                                 TacsScalar *mat) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumFaceQuadraturePoints(faceIndex);

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3], tangent[6];
    double weight = basis->getFaceQuadraturePoint(faceIndex, n, pt, tangent);

    // Get the face normal
    TacsScalar X[3], Xd[9], normal[3];
    TacsScalar area = basis->getFaceNormal(faceIndex, n, Xpts, X, Xd, normal);

    // Compute the inverse of the transformation
    TacsScalar J[9];
    inv3x3(Xd, J);

    // Multiply the weight by the quadrature point
    area *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement3D::MAX_VARS_PER_NODE];
    TacsScalar DUx[3 * TACSElement3D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 3 * varsPerNode * sizeof(TacsScalar));

    if (getTractionComponents) {
      getTractionComponents(elemIndex, time, faceIndex, X, normal, trac);

      // Set the coefficients for the traction
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -trac[k];
      }
    } else if (tractionCoordinateComponent) {
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -trac[k];
      }
    } else {
      for (int k = 0; k < varsPerNode; k++) {
        DUt[3 * k] = -vec3Dot(&trac[3 * k], normal);
      }
    }

    // Add the weak form of the residual at this point
    if (res) {
      basis->addWeakResidual(n, pt, area, J, varsPerNode, DUt, DUx, res);
    }
  }
}
