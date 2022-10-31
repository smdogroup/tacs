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

#include "TACSInertialForce2D.h"

#include "TACSElementAlgebra.h"

TACSInertialForce2D::TACSInertialForce2D(int _varsPerNode,
                                         TACSConstitutive *_con,
                                         TACSElementBasis *_basis,
                                         const TacsScalar _inertiaVec[]) {
  varsPerNode = _varsPerNode;
  con = _con;
  con->incref();
  basis = _basis;
  basis->incref();
  memcpy(inertiaVec, _inertiaVec, 2 * sizeof(TacsScalar));
}

TACSInertialForce2D::~TACSInertialForce2D() {
  con->decref();
  basis->decref();
}

const char *TACSInertialForce2D::getObjectName() {
  return "TACSInertialForce2D";
}

// Get the layout properties of the element
int TACSInertialForce2D::getVarsPerNode() { return varsPerNode; }

int TACSInertialForce2D::getNumNodes() { return basis->getNumNodes(); }

ElementLayout TACSInertialForce2D::getLayoutType() {
  return basis->getLayoutType();
}

TACSElementBasis *TACSInertialForce2D::getElementBasis() { return basis; }

int TACSInertialForce2D::getNumQuadraturePoints() {
  return basis->getNumQuadraturePoints();
}

double TACSInertialForce2D::getQuadratureWeight(int n) {
  return basis->getQuadratureWeight(n);
}

double TACSInertialForce2D::getQuadraturePoint(int n, double pt[]) {
  return basis->getQuadraturePoint(n, pt);
}

int TACSInertialForce2D::getNumElementFaces() {
  return basis->getNumElementFaces();
}

int TACSInertialForce2D::getNumFaceQuadraturePoints(int face) {
  return basis->getNumFaceQuadraturePoints(face);
}

double TACSInertialForce2D::getFaceQuadraturePoint(int face, int n, double pt[],
                                                   double tangent[]) {
  return basis->getFaceQuadraturePoint(face, n, pt, tangent);
}

int TACSInertialForce2D::getDesignVarNums(int elemIndex, int dvLen,
                                          int dvNums[]) {
  return con->getDesignVarNums(elemIndex, dvLen, dvNums);
}

int TACSInertialForce2D::setDesignVars(int elemIndex, int dvLen,
                                       const TacsScalar dvs[]) {
  return con->setDesignVars(elemIndex, dvLen, dvs);
}

int TACSInertialForce2D::getDesignVars(int elemIndex, int dvLen,
                                       TacsScalar dvs[]) {
  return con->getDesignVars(elemIndex, dvLen, dvs);
}

int TACSInertialForce2D::getDesignVarRange(int elemIndex, int dvLen,
                                           TacsScalar lb[], TacsScalar ub[]) {
  return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Add the residual to the provided vector
*/
void TACSInertialForce2D::addResidual(
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
    TacsScalar Xd[4], J[4];
    TacsScalar detXd = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

    // Multiply the quadrature weight by the quadrature point
    TacsScalar volume = weight * detXd;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement2D::MAX_VARS_PER_NODE];
    TacsScalar DUx[2 * TACSElement2D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 2 * varsPerNode * sizeof(TacsScalar));

    // Get the element density
    TacsScalar density = con->evalDensity(elemIndex, pt, Xpts);

    for (int k = 0; k < 2; k++) {
      DUt[3 * k] = -density * inertiaVec[k];
    }

    // Add the weak form of the residual at this point
    basis->addWeakResidual(n, pt, volume, J, varsPerNode, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSInertialForce2D::addJacobian(int elemIndex, double time,
                                      TacsScalar alpha, TacsScalar beta,
                                      TacsScalar gamma, const TacsScalar *Xpts,
                                      const TacsScalar *vars,
                                      const TacsScalar *dvars,
                                      const TacsScalar *ddvars, TacsScalar *res,
                                      TacsScalar *mat) {
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the residual contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the face normal
    TacsScalar Xd[4], J[4];
    TacsScalar detXd = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

    // Multiply the quadrature weight by the quadrature point
    TacsScalar volume = weight * detXd;

    // Evaluate the weak form of the model
    TacsScalar DUt[3 * TACSElement2D::MAX_VARS_PER_NODE];
    TacsScalar DUx[2 * TACSElement2D::MAX_VARS_PER_NODE];
    memset(DUt, 0, 3 * varsPerNode * sizeof(TacsScalar));
    memset(DUx, 0, 2 * varsPerNode * sizeof(TacsScalar));

    // Get the element density
    TacsScalar density = con->evalDensity(elemIndex, pt, Xpts);

    for (int k = 0; k < 2; k++) {
      DUt[3 * k] = -density * inertiaVec[k];
    }

    // Add the weak form of the residual at this point
    if (res) {
      basis->addWeakResidual(n, pt, volume, J, varsPerNode, DUt, DUx, res);
    }
  }
}
