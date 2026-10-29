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

#include <algorithm>

#include "TACSElementAlgebra.h"

TACSCentrifugalForce3D::TACSCentrifugalForce3D(
    int _varsPerNode, TACSConstitutive *_con, TACSElementBasis *_basis,
    const TacsScalar _omegaVec[], const TacsScalar _rotCenter[],
    const int *_omegaDVNums, const int *_rotCenterDVNums) {
  varsPerNode = _varsPerNode;
  con = _con;
  con->incref();
  basis = _basis;
  basis->incref();
  std::copy_n(_omegaVec, 3, omegaVec);
  std::copy_n(_rotCenter, 3, rotCenter);
  if (_omegaDVNums) {
    std::copy_n(_omegaDVNums, 3, omegaDVNums);
  } else {
    omegaDVNums[0] = omegaDVNums[1] = omegaDVNums[2] = -1;
  }
  if (_rotCenterDVNums) {
    std::copy_n(_rotCenterDVNums, 3, rotCenterDVNums);
  } else {
    rotCenterDVNums[0] = rotCenterDVNums[1] = rotCenterDVNums[2] = -1;
  }
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
  int num = con->getDesignVarNums(elemIndex, dvLen, dvNums);
  for (int i = 0; i < 3; i++) {
    if (omegaDVNums[i] >= 0) {
      if (dvNums && num < dvLen) {
        dvNums[num] = omegaDVNums[i];
      }
      num++;
    }
  }
  for (int i = 0; i < 3; i++) {
    if (rotCenterDVNums[i] >= 0) {
      if (dvNums && num < dvLen) {
        dvNums[num] = rotCenterDVNums[i];
      }
      num++;
    }
  }
  return num;
}

int TACSCentrifugalForce3D::setDesignVars(int elemIndex, int dvLen,
                                          const TacsScalar dvs[]) {
  int num = con->setDesignVars(elemIndex, dvLen, dvs);
  for (int i = 0; i < 3; i++) {
    if (omegaDVNums[i] >= 0) {
      if (num < dvLen) {
        omegaVec[i] = dvs[num];
      }
      num++;
    }
  }
  for (int i = 0; i < 3; i++) {
    if (rotCenterDVNums[i] >= 0) {
      if (num < dvLen) {
        rotCenter[i] = dvs[num];
      }
      num++;
    }
  }
  return num;
}

int TACSCentrifugalForce3D::getDesignVars(int elemIndex, int dvLen,
                                          TacsScalar dvs[]) {
  int num = con->getDesignVars(elemIndex, dvLen, dvs);
  for (int i = 0; i < 3; i++) {
    if (omegaDVNums[i] >= 0) {
      if (dvs && num < dvLen) {
        dvs[num] = omegaVec[i];
      }
      num++;
    }
  }
  for (int i = 0; i < 3; i++) {
    if (rotCenterDVNums[i] >= 0) {
      if (dvs && num < dvLen) {
        dvs[num] = rotCenter[i];
      }
      num++;
    }
  }
  return num;
}

int TACSCentrifugalForce3D::getDesignVarRange(int elemIndex, int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[]) {
  int num = con->getDesignVarRange(elemIndex, dvLen, lb, ub);
  for (int i = 0; i < 3; i++) {
    if (omegaDVNums[i] >= 0) {
      if (num < dvLen) {
        lb[num] = -TACS_LARGE_DV_BOUND;
        ub[num] = TACS_LARGE_DV_BOUND;
      }
      num++;
    }
  }
  for (int i = 0; i < 3; i++) {
    if (rotCenterDVNums[i] >= 0) {
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
    TacsScalar X[3], Xd[9], J[9], U[3];
    basis->interpFields(n, pt, 3, Xpts, 1, X);
    basis->interpFields(n, pt, 3, vars, 1, U);
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

    // Create vector pointing from rotation center to element gpt, including
    // deformations
    r[0] = X[0] - rotCenter[0];  // + U[0];
    r[1] = X[1] - rotCenter[1];  // + U[1];
    r[2] = X[2] - rotCenter[2];  // + U[2];

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
    TacsScalar X[3], Xd[9], J[9], U[3];
    basis->interpFields(n, pt, 3, Xpts, 1, X);
    basis->interpFields(n, pt, 3, vars, 1, U);
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

    // Create vector pointing from rotation center to element gpt, including
    // deformations
    r[0] = X[0] - rotCenter[0];  // + U[0];
    r[1] = X[1] - rotCenter[1];  // + U[1];
    r[2] = X[2] - rotCenter[2];  // + U[2];

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

/*
  Add the derivative of the adjoint-residual product psi^T R with respect to
  the element design variables

  The residual at each quadrature point is R_k = volume * density * ac[k] * N,
  where ac = omega x (omega x r) and r = X - rotCenter, so psi^T R =
  sum(volume * dot(ac, b)) with b = density * Psi and Psi the adjoint
  interpolated to the quadrature point. The derivative with respect to the
  constitutive design variables enters through the density. The derivatives
  of dot(ac, b) with respect to omega and r come from centrifugalAccelDotSens,
  and d(r)/d(rotCenter) = -I.
*/
void TACSCentrifugalForce3D::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // The omega and rotCenter design variables follow the constitutive ones
  const int numConDVs = con->getDesignVarNums(elemIndex, dvLen, NULL);

  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the product contribution
  for (int n = 0; n < nquad; n++) {
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the quadrature point location and the Jacobian transformation
    TacsScalar X[3], Xd[9], J[9];
    basis->interpFields(n, pt, 3, Xpts, 1, X);
    TacsScalar detXd = basis->getJacobianTransform(n, pt, Xpts, Xd, J);

    // Multiply the quadrature weight by the quadrature point
    TacsScalar volume = weight * detXd;

    // Interpolate the adjoint variables to the quadrature point
    TacsScalar Psi[TACSElement3D::MAX_VARS_PER_NODE];
    basis->interpFields(n, pt, varsPerNode, psi, 1, Psi);

    // Compute the centrifugal acceleration, ac = omega x (omega x r)
    TacsScalar r[3], wxr[3], ac[3];
    r[0] = X[0] - rotCenter[0];
    r[1] = X[1] - rotCenter[1];
    r[2] = X[2] - rotCenter[2];
    crossProduct(omegaVec, r, wxr);
    crossProduct(omegaVec, wxr, ac);

    // Add the derivative with respect to the density
    con->addDensityDVSens(elemIndex, scale * volume * vec3Dot(ac, Psi), pt, X,
                          dvLen, dfdx);

    // Compute the derivatives of dot(ac, b) with respect to omega and r
    TacsScalar density = con->evalDensity(elemIndex, pt, X);
    TacsScalar b[3], dOmega[3], dr[3];
    b[0] = density * Psi[0];
    b[1] = density * Psi[1];
    b[2] = density * Psi[2];
    centrifugalAccelDotSens(omegaVec, r, b, dOmega, dr);

    int index = numConDVs;
    for (int i = 0; i < 3; i++) {
      if (omegaDVNums[i] >= 0) {
        if (index < dvLen) {
          dfdx[index] += scale * volume * dOmega[i];
        }
        index++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (index < dvLen) {
          dfdx[index] -= scale * volume * dr[i];
        }
        index++;
      }
    }
  }
}
