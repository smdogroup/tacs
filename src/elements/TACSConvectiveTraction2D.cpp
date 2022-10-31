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

#include "TACSConvectiveTraction2D.h"

#include "TACSElementAlgebra.h"

TACSConvectiveTraction2D::TACSConvectiveTraction2D(
    int _varsPerNode, int _faceIndex, int _fieldIndex, TacsScalar _alpha,
    TacsScalar _refValue, TACSElementBasis *_basis) {
  varsPerNode = _varsPerNode;
  faceIndex = _faceIndex;
  fieldIndex = _fieldIndex;
  if (fieldIndex < 0 || fieldIndex >= varsPerNode) {
    fieldIndex = 0;
  }
  alpha = _alpha;
  refValue = _refValue;
  basis = _basis;
  basis->incref();
}

TACSConvectiveTraction2D::~TACSConvectiveTraction2D() { basis->decref(); }

// Get the layout properties of the element
int TACSConvectiveTraction2D::getVarsPerNode() { return varsPerNode; }

int TACSConvectiveTraction2D::getNumNodes() { return basis->getNumNodes(); }

ElementLayout TACSConvectiveTraction2D::getLayoutType() {
  return basis->getLayoutType();
}

TACSElementBasis *TACSConvectiveTraction2D::getElementBasis() { return basis; }

int TACSConvectiveTraction2D::getNumQuadraturePoints() {
  return basis->getNumQuadraturePoints();
}

double TACSConvectiveTraction2D::getQuadratureWeight(int n) {
  return basis->getQuadratureWeight(n);
}

double TACSConvectiveTraction2D::getQuadraturePoint(int n, double pt[]) {
  return basis->getQuadraturePoint(n, pt);
}

int TACSConvectiveTraction2D::getNumElementFaces() {
  return basis->getNumElementFaces();
}

int TACSConvectiveTraction2D::getNumFaceQuadraturePoints(int face) {
  return basis->getNumFaceQuadraturePoints(face);
}

double TACSConvectiveTraction2D::getFaceQuadraturePoint(int face, int n,
                                                        double pt[],
                                                        double tangent[]) {
  return basis->getFaceQuadraturePoint(face, n, pt, tangent);
}

/*
  Add the residual to the provided vector
*/
void TACSConvectiveTraction2D::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumFaceQuadraturePoints(faceIndex);

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3], tangent[6];
    double weight = basis->getFaceQuadraturePoint(faceIndex, n, pt, tangent);

    // Get the face normal
    TacsScalar X[3], Xd[4], normal[3];
    TacsScalar area = basis->getFaceNormal(faceIndex, n, Xpts, X, Xd, normal);

    // Get the field value
    TacsScalar U[TACSElement2D::MAX_VARS_PER_NODE];
    basis->getFieldValues(-1, pt, Xpts, varsPerNode, vars, X, U);

    // Compute the inverse of the transformation
    TacsScalar J[4];
    inv2x2(Xd, J);

    // Multiply the weight by the quadrature point
    area *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement2D::MAX_VARS_PER_NODE];
    TacsScalar DUx[2 * TACSElement2D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 2 * varsPerNode * sizeof(TacsScalar));

    // Compute the contribution to the right-hand-side
    DUt[3 * fieldIndex] = -alpha * (U[fieldIndex] - refValue);

    // Add the weak form of the residual at this point
    basis->addWeakResidual(n, pt, area, J, varsPerNode, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSConvectiveTraction2D::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *mat) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumFaceQuadraturePoints(faceIndex);

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3], tangent[6];
    double weight = basis->getFaceQuadraturePoint(faceIndex, n, pt, tangent);

    // Get the face normal
    TacsScalar X[3], Xd[4], normal[3];
    TacsScalar area = basis->getFaceNormal(faceIndex, n, Xpts, X, Xd, normal);

    // Get the field value
    TacsScalar U[TACSElement2D::MAX_VARS_PER_NODE];
    basis->getFieldValues(-1, pt, Xpts, varsPerNode, vars, X, U);

    // Compute the inverse of the transformation
    TacsScalar J[4];
    inv2x2(Xd, J);

    // Multiply the weight by the quadrature point
    area *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement2D::MAX_VARS_PER_NODE];
    TacsScalar DUx[2 * TACSElement2D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 2 * varsPerNode * sizeof(TacsScalar));

    // Compute the contribution to the right-hand-side
    DUt[3 * fieldIndex] = -alpha * (U[fieldIndex] - refValue);

    // Add the weak form of the residual at this point
    if (res) {
      basis->addWeakResidual(n, pt, area, J, varsPerNode, DUt, DUx, res);
    }

    // Add the weak form of the residual at this point
    int Jac_nnz = 1;
    TacsScalar Jac = -alpha;
    int Jac_pairs[2] = {3 * fieldIndex, 3 * fieldIndex};

    // Add the weak form of the residual at this point
    basis->scaleWeakMatrix(area, alpha, beta, gamma, Jac_nnz, Jac_pairs, &Jac);
    basis->addWeakMatrix(n, pt, J, varsPerNode, Jac_nnz, Jac_pairs, &Jac, mat);
  }
}

void TACSConvectiveTraction2D::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {}

void TACSConvectiveTraction2D::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdXpts[]) {}
