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

#include <algorithm>

#include "TACSElementAlgebra.h"

TACSTraction3D::TACSTraction3D(int _varsPerNode, int _faceIndex,
                               TACSElementBasis *_basis,
                               const TacsScalar _trac[],
                               int _tractionCoordinateComponent,
                               const int *_tracDVNums) {
  varsPerNode = _varsPerNode;
  faceIndex = _faceIndex;
  basis = _basis;
  basis->incref();
  tractionCoordinateComponent = _tractionCoordinateComponent;
  getTractionComponents = NULL;
  if (tractionCoordinateComponent) {
    std::copy_n(_trac, varsPerNode, trac);
  } else {
    std::copy_n(_trac, 3 * varsPerNode, trac);
  }
  initTracDVNums(_tracDVNums);
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
  initTracDVNums(NULL);
}

void TACSTraction3D::initTracDVNums(const int *_tracDVNums) {
  if (_tracDVNums) {
    std::copy_n(_tracDVNums, varsPerNode, tracDVNums);
  } else {
    for (int i = 0; i < varsPerNode; i++) {
      tracDVNums[i] = -1;
    }
  }
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

int TACSTraction3D::getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
  int num = 0;
  for (int i = 0; i < varsPerNode; i++) {
    if (tracDVNums[i] >= 0) {
      if (dvNums && num < dvLen) {
        dvNums[num] = tracDVNums[i];
      }
      num++;
    }
  }
  return num;
}

int TACSTraction3D::setDesignVars(int elemIndex, int dvLen,
                                  const TacsScalar dvs[]) {
  int num = 0;
  for (int i = 0; i < varsPerNode; i++) {
    if (tracDVNums[i] >= 0) {
      if (num < dvLen) {
        trac[i] = dvs[num];
      }
      num++;
    }
  }
  return num;
}

int TACSTraction3D::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  int num = 0;
  for (int i = 0; i < varsPerNode; i++) {
    if (tracDVNums[i] >= 0) {
      if (dvs && num < dvLen) {
        dvs[num] = trac[i];
      }
      num++;
    }
  }
  return num;
}

int TACSTraction3D::getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                                      TacsScalar ub[]) {
  int num = 0;
  for (int i = 0; i < varsPerNode; i++) {
    if (tracDVNums[i] >= 0) {
      if (num < dvLen) {
        lb[num] = -TACS_LARGE_DV_BOUND;
        ub[num] = TACS_LARGE_DV_BOUND;
      }
      num++;
    }
  }
  return num;
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

/*
  Add the derivative of the adjoint-residual product with respect to the
  traction design variables

  The residual is the sum over the face quadrature points of
  area * N^T DUt, so psi^T R = sum area * sum_k DUt[3k] * Psi[k], where Psi
  is the adjoint interpolated to the quadrature point. Design variable i sets
  trac[i]. With coordinate components, DUt[3k] = -trac[k], so the derivative
  is -area * Psi[i]. Otherwise trac stores one 3-vector per field and
  DUt[3k] = -dot(trac[3k:3k+3], normal), so the derivative is
  -area * normal[i % 3] * Psi[i / 3].
*/
void TACSTraction3D::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // Tractions from a user-supplied callback are never design variables
  if (getTractionComponents) {
    return;
  }

  // Compute the number of quadrature points
  const int nquad = basis->getNumFaceQuadraturePoints(faceIndex);

  // Loop over each quadrature point and add the derivative contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3], tangent[6];
    double weight = basis->getFaceQuadraturePoint(faceIndex, n, pt, tangent);

    // Get the face normal
    TacsScalar X[3], Xd[9], normal[3];
    TacsScalar area = basis->getFaceNormal(faceIndex, n, Xpts, X, Xd, normal);

    // Multiply the weight by the quadrature point
    area *= weight;

    // Interpolate the adjoint to the quadrature point
    TacsScalar Psi[TACSElement3D::MAX_VARS_PER_NODE];
    basis->interpFields(n, pt, varsPerNode, psi, 1, Psi);

    int index = 0;
    for (int i = 0; i < varsPerNode; i++) {
      if (tracDVNums[i] >= 0) {
        if (index < dvLen) {
          if (tractionCoordinateComponent) {
            dfdx[index] -= scale * area * Psi[i];
          } else {
            dfdx[index] -= scale * area * normal[i % 3] * Psi[i / 3];
          }
        }
        index++;
      }
    }
  }
}
