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

#include "TACSCentrifugalForce3D.h"

#include "TACSElementAlgebra.h"

TACSCentrifugalForce3D::TACSCentrifugalForce3D(int _varsPerNode,
                                               TACSConstitutive *_con,
                                               TACSElementBasis *_basis,
                                               const TacsScalar _omegaVec[],
                                               const TacsScalar _rotCenter[]) {
  varsPerNode = _varsPerNode;
  con = _con;
  con->incref();
  basis = _basis;
  basis->incref();
  memcpy(omegaVec, _omegaVec, 3 * sizeof(TacsScalar));
  memcpy(rotCenter, _rotCenter, 3 * sizeof(TacsScalar));
}

TACSCentrifugalForce3D::~TACSCentrifugalForce3D() {
  con->decref();
  basis->decref();
}

const char *TACSCentrifugalForce3D::getObjectName() {
  return "TACSCentrifugalForce3D";
}

// Get the layout properties of the element
int TACSCentrifugalForce3D::getVarsPerNode() { return varsPerNode; }

int TACSCentrifugalForce3D::getNumNodes() { return basis->getNumNodes(); }

ElementLayout TACSCentrifugalForce3D::getLayoutType() {
  return basis->getLayoutType();
}

TACSElementBasis *TACSCentrifugalForce3D::getElementBasis() { return basis; }

int TACSCentrifugalForce3D::getNumQuadraturePoints() {
  return basis->getNumQuadraturePoints();
}

double TACSCentrifugalForce3D::getQuadratureWeight(int n) {
  return basis->getQuadratureWeight(n);
}

double TACSCentrifugalForce3D::getQuadraturePoint(int n, double pt[]) {
  return basis->getQuadraturePoint(n, pt);
}

int TACSCentrifugalForce3D::getNumElementFaces() {
  return basis->getNumElementFaces();
}

int TACSCentrifugalForce3D::getNumFaceQuadraturePoints(int face) {
  return basis->getNumFaceQuadraturePoints(face);
}

double TACSCentrifugalForce3D::getFaceQuadraturePoint(int face, int n,
                                                      double pt[],
                                                      double tangent[]) {
  return basis->getFaceQuadraturePoint(face, n, pt, tangent);
}

int TACSCentrifugalForce3D::getDesignVarNums(int elemIndex, int dvLen,
                                             int dvNums[]) {
  return con->getDesignVarNums(elemIndex, dvLen, dvNums);
}

int TACSCentrifugalForce3D::setDesignVars(int elemIndex, int dvLen,
                                          const TacsScalar dvs[]) {
  return con->setDesignVars(elemIndex, dvLen, dvs);
}

int TACSCentrifugalForce3D::getDesignVars(int elemIndex, int dvLen,
                                          TacsScalar dvs[]) {
  return con->getDesignVars(elemIndex, dvLen, dvs);
}

int TACSCentrifugalForce3D::getDesignVarRange(int elemIndex, int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[]) {
  return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Add the residual to the provided vector
*/
void TACSCentrifugalForce3D::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the face normal
    TacsScalar X[3], Xd[9], J[9];
    basis->interpFields(n, pt, 3, Xpts, 1, X);
    TacsScalar detXd = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

    // Multiply the quadrature weight by the quadrature point
    TacsScalar volume = weight * detXd;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement3D::MAX_VARS_PER_NODE];
    TacsScalar DUx[3 * TACSElement3D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 3 * varsPerNode * sizeof(TacsScalar));

    // Get the element density
    TacsScalar density = con->evalDensity(elemIndex, pt, X);

    TacsScalar r[3], wxr[3], ac[3];

    // Create vector pointing from rotation center to element gpt
    r[0] = X[0] - rotCenter[0];
    r[1] = X[1] - rotCenter[1];
    r[2] = X[2] - rotCenter[2];

    // Compute omega x r
    crossProduct(omegaVec, r, wxr);

    // Compute centrifugal acceleration
    crossProduct(omegaVec, wxr, ac);

    for (int k = 0; k < 3; k++) {
      DUt[3 * k] = density * ac[k];
    }

    // Add the weak form of the residual at this point
    basis->addWeakResidual(n, pt, volume, J, varsPerNode, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSCentrifugalForce3D::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *mat) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the face normal
    TacsScalar X[3], Xd[9], J[9];
    basis->interpFields(n, pt, 3, Xpts, 1, X);
    TacsScalar detXd = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

    // Multiply the weight by the quadrature point
    TacsScalar volume = weight * detXd;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement3D::MAX_VARS_PER_NODE];
    TacsScalar DUx[3 * TACSElement3D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 3 * varsPerNode * sizeof(TacsScalar));

    // Get the element density
    TacsScalar density = con->evalDensity(elemIndex, pt, X);

    TacsScalar r[3], wxr[3], ac[3];

    // Create vector pointing from rotation center to element gpt
    r[0] = X[0] - rotCenter[0];
    r[1] = X[1] - rotCenter[1];
    r[2] = X[2] - rotCenter[2];

    // Compute omega x r
    crossProduct(omegaVec, r, wxr);

    // Compute centrifugal acceleration
    crossProduct(omegaVec, wxr, ac);

    for (int k = 0; k < 3; k++) {
      DUt[3 * k] = density * ac[k];
    }

    // Add the weak form of the residual at this point
    if (res) {
      basis->addWeakResidual(n, pt, volume, J, varsPerNode, DUt, DUx, res);
    }
  }
}
