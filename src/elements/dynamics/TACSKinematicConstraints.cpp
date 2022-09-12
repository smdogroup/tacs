/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSKinematicConstraints.h"

#include "TACSElementAlgebra.h"
#include "TACSElementQuaternion.h"
#include "TACSGaussQuadrature.h"
#include "TACSRigidBody.h"

/*
  Construct a spherical constraint with the two bodies involved and a
  position vector measured from the global frame to the point where
  the spherical joint is located.
*/
TACSSphericalConstraint::TACSSphericalConstraint(
    TACSRigidBody *_bodyA, TACSRigidBody *_bodyB, TACSGibbsVector *_point,
    TACSGibbsVector *_axis, TACSTranslationConType _con_type) {
  // Copy over the arguments
  inertial_fixed_point = 0;
  bodyA = _bodyA;
  bodyA->incref();
  bodyB = _bodyB;
  bodyB->incref();
  point = _point;
  point->incref();
  axis = _axis;
  con_type = _con_type;
  if (axis) {
    axis->incref();
  } else {
    con_type = COINCIDENT;
  }
}

/*
  Construct a spherical constraint with a body involved and a position
  vector measured from the global frame to the point where the
  spherical joint is located.
*/
TACSSphericalConstraint::TACSSphericalConstraint(
    TACSRigidBody *_bodyA, TACSGibbsVector *_point, TACSGibbsVector *_axis,
    TACSTranslationConType _con_type) {
  // Copy over the arguments
  inertial_fixed_point = 1;
  bodyA = _bodyA;
  bodyA->incref();
  bodyB = NULL;
  point = _point;
  point->incref();
  axis = _axis;
  con_type = _con_type;
  if (axis) {
    axis->incref();
  } else {
    con_type = COINCIDENT;
  }
}

/*
  Destructor for spherical constraint
*/
TACSSphericalConstraint::~TACSSphericalConstraint() {
  if (bodyA) {
    bodyA->decref();
  }
  if (bodyB) {
    bodyB->decref();
  }
  point->decref();
  if (axis) {
    axis->decref();
  }
}

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSSphericalConstraint::getNumNodes() {
  if (inertial_fixed_point) {
    return 2;
  } else {
    return 3;
  }
}

/*
  Returns the multiplier index
*/
int TACSSphericalConstraint::getMultiplierIndex() {
  if (inertial_fixed_point) {
    return 1;
  } else {
    return 2;
  }
}

const char *TACSSphericalConstraint::elem_name = "TACSSphericalConstraint";

/*
  Compute the residual of the governing equations
*/
void TACSSphericalConstraint::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Retrieve the spherical joint location in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the directions that are orthogonal to the axis
  TacsScalar dir1[3], dir2[3];
  dir1[0] = dir1[1] = dir1[2] = 0.0;
  dir2[0] = dir2[1] = dir2[2] = 0.0;

  // If the constraint is colinear, then compute the directions dir1
  // and dir2 based on the axis directions.
  if (con_type == COLINEAR) {
    // Get the axis vector
    const TacsScalar *a;
    axis->getVector(&a);

    // Compute the absolute values of the
    double a0 = fabs(TacsRealPart(a[0]));
    double a1 = fabs(TacsRealPart(a[1]));
    double a2 = fabs(TacsRealPart(a[2]));

    // Compute a vector that is orthogonal to the axis
    TacsScalar e[3];
    e[0] = e[1] = e[2] = 0.0;
    if ((a0 <= a1) && (a0 <= a2)) {
      e[0] = 1.0;
    } else if ((a1 <= a0) && (a1 <= a2)) {
      e[1] = 1.0;
    } else {
      e[2] = 1.0;
    }

    // Compute dir1 = (axis x e) and dir2 = (dir1 x axis)
    crossProduct(1.0, a, e, dir1);
    crossProduct(1.0, dir1, a, dir2);
  }

  // Get the initial position for bodyA and bodyB - if they exist
  const TacsScalar *rA = NULL, *rB = NULL;
  if (bodyA) {
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  } else {
    rA = &Xpts[0];
  }
  if (bodyB) {
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  } else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set pointers to the residual of each body
  TacsScalar *resA = &res[0];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar *uB = &vars[8];
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // The residual for the constraint equations
  TacsScalar *resC = NULL;

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not.
  if (inertial_fixed_point) {
    resC = &res[8];
    lam = &vars[8];
  } else {
    resC = &res[16];
    lam = &vars[16];
  }

  // Compute the rotation matrix for the first body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Compute the distance between the point
  // r = rA + uA + CA^{T}*xA - pt = 0
  TacsScalar r[3];
  mat3x3MultTrans(CA, xA, r);
  vec3Axpy(1.0, uA, r);
  vec3Axpy(1.0, rA, r);

  if (inertial_fixed_point) {
    // Finish computing the constraint
    vec3Axpy(-1.0, pt, r);

    // Add the residual
    if (con_type == COINCIDENT) {
      resC[0] += r[0];
      resC[1] += r[1];
      resC[2] += r[2];

      // Add the terms for body A
      vec3Axpy(1.0, lam, &resA[0]);
      addEMatTransProduct(1.0, xA, lam, etaA, epsA, &resA[3], &resA[4]);
    } else if (con_type == COLINEAR) {
      // resC = dot(rA + uA + CA^{T}*xA - pt, dir) = 0
      resC[0] += vec3Dot(r, dir1);
      resC[1] += vec3Dot(r, dir2);
      resC[2] += lam[2];

      // Add the terms for body A
      resA[0] += lam[0] * dir1[0] + lam[1] * dir2[0];
      resA[1] += lam[0] * dir1[1] + lam[1] * dir2[1];
      resA[2] += lam[0] * dir1[2] + lam[1] * dir2[2];

      // Add the terms from the transpose derivative of the constraint
      // times the multiplier
      addEMatTransProduct(lam[0], xA, dir1, etaA, epsA, &resA[3], &resA[4]);
      addEMatTransProduct(lam[1], xA, dir2, etaA, epsA, &resA[3], &resA[4]);
    } else {  // con_type == COPLANAR
      // Get the axis vector
      const TacsScalar *a;
      axis->getVector(&a);

      // resC = dot(rA + uA + CA^{T}*xA - pt, axis) = 0
      resC[0] += vec3Dot(r, a);
      resC[1] += lam[1];
      resC[2] += lam[2];

      // Add the terms from the transpose derivative of the constraint
      // times the multiplier
      vec3Axpy(lam[0], a, &resA[0]);
      addEMatTransProduct(lam[0], xA, a, etaA, epsA, &resA[3], &resA[4]);
    }
  } else {
    // Set the pointer into the residual
    TacsScalar *resB = &res[8];

    // Compute the rotation matrix for bodyB
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);

    // Compute t = CB^{T}*xB + uB
    TacsScalar t[3];
    mat3x3MultTrans(CB, xB, t);
    vec3Axpy(1.0, uB, t);
    vec3Axpy(1.0, rB, t);

    // r = rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB = 0
    vec3Axpy(-1.0, t, r);

    if (con_type == COINCIDENT) {
      // Set the point so that it is coincident
      resC[0] += r[0];
      resC[1] += r[1];
      resC[2] += r[2];

      // Add the terms for body A
      vec3Axpy(1.0, lam, &resA[0]);
      addEMatTransProduct(1.0, xA, lam, etaA, epsA, &resA[3], &resA[4]);

      // Add the terms for body B
      vec3Axpy(-1.0, lam, &resB[0]);
      addEMatTransProduct(-1.0, xB, lam, etaB, epsB, &resB[3], &resB[4]);
    } else if (con_type == COLINEAR) {
      // Pre-multiply to get the local directions
      TacsScalar d1[3], d2[3];
      mat3x3MultTrans(CA, dir1, d1);
      mat3x3MultTrans(CA, dir2, d2);

      // resC = dot(rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB, dir) = 0
      resC[0] += vec3Dot(r, d1);
      resC[1] += vec3Dot(r, d2);
      resC[2] += lam[2];

      // Add the terms for body A
      resA[0] += lam[0] * d1[0] + lam[1] * d2[0];
      resA[1] += lam[0] * d1[1] + lam[1] * d2[1];
      resA[2] += lam[0] * d1[2] + lam[1] * d2[2];

      // Add the terms for body B
      resB[0] -= lam[0] * d1[0] + lam[1] * d2[0];
      resB[1] -= lam[0] * d1[1] + lam[1] * d2[1];
      resB[2] -= lam[0] * d1[2] + lam[1] * d2[2];

      // Add the terms from the transpose derivative of the constraint
      // times the multiplier
      addEMatTransProduct(lam[0], xA, d1, etaA, epsA, &resA[3], &resA[4]);
      addEMatTransProduct(lam[1], xA, d2, etaA, epsA, &resA[3], &resA[4]);

      addEMatTransProduct(-lam[0], xA, d1, etaA, epsA, &resB[3], &resB[4]);
      addEMatTransProduct(-lam[1], xA, d2, etaA, epsA, &resB[3], &resB[4]);
    } else {  // con_type == COPLANAR
      // Get the axis vector
      const TacsScalar *a;
      axis->getVector(&a);

      // Get the vector in the A-axis
      TacsScalar aA[3];
      mat3x3MultTrans(CA, a, aA);

      // resC = dot(rA + uA + CA^{T}*xA - pt, axis) = 0
      resC[0] += vec3Dot(r, aA);
      resC[1] += lam[1];
      resC[2] += lam[2];

      // Add the terms from the transpose derivative of the constraint
      // times the multiplier
      vec3Axpy(lam[0], aA, &resA[0]);
      addEMatTransProduct(lam[0], xA, aA, etaA, epsA, &resA[3], &resA[4]);
    }
  }

  // Add the dummy constraints for the remaining Lagrange multiplier
  // variables
  for (int i = 3; i < 8; i++) {
    resC[i] += lam[i];
  }
}

/*
  Compute the Jacobian of the residuals of the governing equations
*/
void TACSSphericalConstraint::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *J) {
  // Add the residual
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  // Retrieve the spherical joint location in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the directions that are orthogonal to the axis
  TacsScalar dir1[3], dir2[3];
  dir1[0] = dir1[1] = dir1[2] = 0.0;
  dir2[0] = dir2[1] = dir2[2] = 0.0;

  // If the constraint is colinear, then compute the directions dir1
  // and dir2 based on the axis directions.
  if (con_type == COLINEAR) {
    // Get the axis vector
    const TacsScalar *a;
    axis->getVector(&a);

    // Compute the absolute values of the
    double a0 = fabs(TacsRealPart(a[0]));
    double a1 = fabs(TacsRealPart(a[1]));
    double a2 = fabs(TacsRealPart(a[2]));

    // Compute a vector that is orthogonal to the axis
    TacsScalar e[3];
    e[0] = e[1] = e[2] = 0.0;
    if ((a0 <= a1) && (a0 <= a2)) {
      e[0] = 1.0;
    } else if ((a1 <= a0) && (a1 <= a2)) {
      e[1] = 1.0;
    } else {
      e[2] = 1.0;
    }

    // Compute dir1 = (axis x e) and dir2 = (dir1 x axis)
    crossProduct(1.0, a, e, dir1);
    crossProduct(1.0, dir1, a, dir2);
  }

  // Get the initial position for bodyA and bodyB - if they exist
  const TacsScalar *rA = NULL, *rB = NULL;
  if (bodyA) {
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  } else {
    rA = &Xpts[0];
  }
  if (bodyB) {
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  } else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // Set the Lagrange multipliers for the constraint
  const TacsScalar *lam = NULL;

  // Get the number of variables
  const int nvars = getNumVariables();

  // Set the offset to the Lagrange multipliers
  int offset = 0;
  if (inertial_fixed_point) {
    offset = 8;
    lam = &vars[8];
  } else {
    offset = 16;
    lam = &vars[16];
  }

  if (inertial_fixed_point) {
    if (con_type == COINCIDENT) {
      // Add the identity matricies to the Jacobian
      addBlockIdent(alpha, &J[offset], nvars);
      addBlockIdent(alpha, &J[offset * nvars], nvars);

      // Add the second derivative terms
      addBlockDMatTransDeriv(alpha, lam, xA, &J[3 * (nvars + 1)], nvars);

      // Add the term from the derivative w.r.t. lambda
      addBlockEMatTrans(alpha, etaA, epsA, xA, &J[3 * nvars + offset], nvars);

      // Add the term from the derivative of the constraint
      addBlockEMat(alpha, etaA, epsA, xA, &J[offset * nvars + 3], nvars);
    } else if (con_type == COLINEAR) {
      // resC = dot(rA, uA + CA^{T}*xA - pt, dir) = 0
      addVecMat(alpha, dir1, &J[offset], nvars);
      addVecMat(alpha, dir2, &J[offset + 1], nvars);

      // Add the derivative of to the matrix dir1^{T}*CA^{T}*xA

      J[(offset + 2) * (nvars + 1)] += alpha;
    } else {  // con_type == COPLANAR
    }
  } else {
    // Add the terms required for body B if it is defined
    if (con_type == COINCIDENT) {
      // Add the identity matricies to the Jacobian
      addBlockIdent(alpha, &J[offset], nvars);
      addBlockIdent(alpha, &J[offset * nvars], nvars);

      // Add the second derivative terms
      addBlockDMatTransDeriv(alpha, lam, xA, &J[3 * (nvars + 1)], nvars);

      // Add the term from the derivative w.r.t. lambda
      addBlockEMatTrans(alpha, etaA, epsA, xA, &J[3 * nvars + offset], nvars);

      // Add the term from the derivative of the constraint
      addBlockEMat(alpha, etaA, epsA, xA, &J[offset * nvars + 3], nvars);

      // Add the identity matrices to the Jacobian
      addBlockIdent(-alpha, &J[8 * nvars + offset], nvars);
      addBlockIdent(-alpha, &J[offset * nvars + 8], nvars);

      // Add the second derivative terms
      addBlockDMatTransDeriv(-alpha, lam, xB, &J[11 * (nvars + 1)], nvars);

      // Add the terms from the derivatives w.r.t. lambdas
      addBlockEMatTrans(-alpha, etaB, epsB, xB, &J[11 * nvars + offset], nvars);

      // Add the terms from the derivatives of the constraint
      addBlockEMat(-alpha, etaB, epsB, xB, &J[offset * nvars + 11], nvars);
    } else if (con_type == COLINEAR) {
    } else {  // con_type == COPLANAR
    }
  }

  // Add the Jacobian entries for the dummy constraints
  for (int i = offset + 3; i < nvars; i++) {
    J[(nvars + 1) * i] += alpha;
  }
}

/*
  Get the design variable numbers
*/
int TACSSphericalConstraint::getDesignVarNums(int elemIndex, int dvLen,
                                              int dvNums[]) {
  return point->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the design variable values
*/
int TACSSphericalConstraint::setDesignVars(int elemIndex, int dvLen,
                                           const TacsScalar dvs[]) {
  return point->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable values associated with the joint location
*/
int TACSSphericalConstraint::getDesignVars(int elemIndex, int dvLen,
                                           TacsScalar dvs[]) {
  return point->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable values associated with the joint location
*/
int TACSSphericalConstraint::getDesignVarRange(int elemIndex, int dvLen,
                                               TacsScalar lb[],
                                               TacsScalar ub[]) {
  return point->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Constructor for revolute contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA:  pointer to bodyA
  bodyB:  pointer to bodyB
  point:  the position of the joint from the global reference point
  rev:    the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint(TACSRigidBody *_bodyA,
                                               TACSRigidBody *_bodyB,
                                               TACSGibbsVector *_point,
                                               TACSGibbsVector *_eAVec,
                                               int _inertial_rev_axis) {
  // Copy over the input arguments
  inertial_fixed_point = 0;
  bodyA = _bodyA;
  bodyA->incref();
  bodyB = _bodyB;
  bodyB->incref();
  point = _point;
  point->incref();
  eAVec = _eAVec;
  eAVec->incref();
  inertial_rev_axis = _inertial_rev_axis;

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Constructor for revolute contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA:  pointer to bodyA
  point:  the position of the joint from the global reference point
  rev:    the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint(TACSRigidBody *_bodyA,
                                               TACSGibbsVector *_point,
                                               TACSGibbsVector *_eAVec) {
  // Copy over the input arguments
  inertial_fixed_point = 1;
  inertial_rev_axis = 1;
  bodyA = _bodyA;
  bodyA->incref();
  bodyB = NULL;
  point = _point;
  point->incref();
  eAVec = _eAVec;
  eAVec->incref();

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Create a revolute constraint joining two nodes together directly -
  this can be used to connect two flexible bodies

  input:
  inertial_fixed_point: Is the reference point fixed?
  point:   the position of the joint from the global reference point
  rev:     the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint(int _inertial_fixed_point,
                                               TACSGibbsVector *_point,
                                               TACSGibbsVector *_eAVec,
                                               int _inertial_rev_axis) {
  // Copy over the input arguments
  inertial_fixed_point = _inertial_fixed_point;
  bodyA = bodyB = NULL;
  point = _point;
  point->incref();
  eAVec = _eAVec;
  eAVec->incref();
  inertial_rev_axis = _inertial_rev_axis;

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Destuctor for the revolute constraint
*/
TACSRevoluteConstraint::~TACSRevoluteConstraint() {
  if (bodyA) {
    bodyA->decref();
  }
  if (bodyB) {
    bodyB->decref();
  }
  point->decref();
  eAVec->decref();
  eVec->decref();
  eB1Vec->decref();
  eB2Vec->decref();
}

void TACSRevoluteConstraint::setRevoluteAxis(TACSGibbsVector *_eAVec) {
  eAVec = _eAVec;
  int init_vector = 1;
  updatePoints(init_vector);
}

const char *TACSRevoluteConstraint::elem_name = "TACSRevoluteConstraint";

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSRevoluteConstraint::getNumNodes() {
  if (inertial_fixed_point) {
    return 2;
  } else {
    return 3;
  }
}

/*
  Get the local multiplier node number
*/
int TACSRevoluteConstraint::getMultiplierIndex() {
  if (inertial_fixed_point) {
    return 1;
  } else {
    return 2;
  }
}

/*
  Read the data from the given initial point vectors/locations and
  re-compute the internal data that is requied to evaluate the
  kinematic constraints.
*/
void TACSRevoluteConstraint::updatePoints(int init_vector) {
  // Find the minimum absolute component of eAVec along any coordinate
  // direction. Set the vector components of e along this direction
  // to maximize orthogonality among the coordinate directions. For
  // the purpose of optimization, this direction is fixed at
  // initialization.
  // Retrieve the revolute direction in global frame
  const TacsScalar *rev;
  eAVec->getVector(&rev);

  TacsScalar e[3];
  if (init_vector) {
    e[0] = e[1] = e[2] = 0.0;
    if ((fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[1]))) &&
        (fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[2])))) {
      e[0] = 1.0;
    } else if ((fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[0]))) &&
               (fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[2])))) {
      e[1] = 1.0;
    } else {
      e[2] = 1.0;
    }
    eVec = new TACSGibbsVector(e);
    eVec->incref();
  } else {
    const TacsScalar *etmp;
    eVec->getVector(&etmp);
    memcpy(e, etmp, 3 * sizeof(TacsScalar));
  }

  // Compute/recompute the eB1 and eB2 directions based on e
  TacsScalar eB1[3], eB2[3];
  crossProduct(1.0, rev, e, eB2);
  if (eB2Vec) {
    eB2Vec->decref();
  }
  eB2Vec = new TACSGibbsVector(eB2);
  eB2Vec->incref();

  crossProduct(1.0, eB2, rev, eB1);
  if (eB1Vec) {
    eB1Vec->decref();
  }
  eB1Vec = new TACSGibbsVector(eB1);
  eB1Vec->incref();
}

/*
  Compute the residual of the governing equations
*/
void TACSRevoluteConstraint::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Get the initial position vectors - depending on whether we're
  // taking the initial position from the node locations
  const TacsScalar *rA, *rB;
  if (bodyA) {
    // Get the initial position for bodyA
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  } else {
    rA = &Xpts[0];
  }
  if (bodyB) {
    // Get the initial position for bodyA
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  } else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set pointers to the residual of each body
  TacsScalar *resA = &res[0];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar *uB = &vars[8];
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // The residual for the constraint equations
  TacsScalar *resC = NULL;

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not.
  if (inertial_fixed_point) {
    resC = &res[8];
    lam = &vars[8];
  } else {
    resC = &res[16];
    lam = &vars[16];
  }

  // Retrieve the pointers to eA, eB1, eB2
  const TacsScalar *eA, *eB1, *eB2;
  eAVec->getVector(&eA);
  eB1Vec->getVector(&eB1);
  eB2Vec->getVector(&eB2);

  // Compute the rotation matrices
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Evaluate the constraint
  // resC = rA + uA + CA^{T}*xA - pt = 0 or
  // resC = rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB = 0
  mat3x3MultTransAdd(CA, xA, resC);
  vec3Axpy(1.0, uA, resC);
  vec3Axpy(1.0, rA, resC);

  // Add the terms for body A
  vec3Axpy(1.0, lam, &resA[0]);
  addEMatTransProduct(1.0, xA, lam, etaA, epsA, &resA[3], &resA[4]);

  // Next, handle the revolute constraint
  if (inertial_fixed_point) {
    // Subtract pt vec if bodyB is not present
    vec3Axpy(-1.0, pt, resC);

    // Add the revolute direction constraint
    TacsScalar tA[3];
    mat3x3MultTrans(CA, eA, tA);

    // Compute the contributions to the first revolute constraint
    resC[3] += vec3Dot(tA, eB1);

    // Add the derivative (d(CA^{T}*eA)/dqA)*eB1 = E(eA)*tB
    addEMatTransProduct(lam[3], eA, eB1, etaA, epsA, &resA[3], &resA[4]);

    // Compute the contributions to the second revolute constraint
    resC[4] += vec3Dot(tA, eB2);

    // Add the derivative (d(CA^{T}*eA)/dqA)*eB2 = E(eA)*eB2
    addEMatTransProduct(lam[4], eA, eB2, etaA, epsA, &resA[3], &resA[4]);

    // Add the dummy constraints
    for (int i = 5; i < 8; i++) {
      resC[i] += lam[i];
    }
  } else {
    // Set the residual for bodyB
    TacsScalar *resB = &res[8];

    // Compute the rotation matrix for bodyB
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);

    // Compute t = rB + uB + CB^{T}*xB
    TacsScalar t[3];
    mat3x3MultTrans(CB, xB, t);
    vec3Axpy(1.0, uB, t);
    vec3Axpy(1.0, rB, t);

    // Complete the evaluation of the constraint
    vec3Axpy(-1.0, t, resC);

    // Add the terms for body B
    vec3Axpy(-1.0, lam, &resB[0]);
    addEMatTransProduct(-1.0, xB, lam, etaB, epsB, &resB[3], &resB[4]);

    if (inertial_rev_axis) {
      // In this case, eA, eB1 and eB2 are all in the inertial
      // reference frame, and fixed in this frame. This code
      // adds a revolute constraint relative to the reference
      // frame for both body A and body B
      TacsScalar tA[3];
      mat3x3MultTrans(CA, eA, tA);

      // Compute the contributions to the first revolute constraint
      resC[3] += vec3Dot(tA, eB1);
      addEMatTransProduct(lam[3], eA, eB1, etaA, epsA, &resA[3], &resA[4]);

      // Compute the contributions to the second revolute constraint
      resC[4] += vec3Dot(tA, eB2);
      addEMatTransProduct(lam[4], eA, eB2, etaA, epsA, &resA[3], &resA[4]);

      // Compute tB = CB^{T}*eA
      TacsScalar tB[3];
      mat3x3MultTrans(CB, eA, tB);

      // Compute the contributions to the first revolute constraint
      resC[5] += vec3Dot(tB, eB1);
      addEMatTransProduct(lam[5], eA, eB1, etaB, epsB, &resB[3], &resB[4]);

      // Compute the contributions to the second revolute constraint
      resC[6] += vec3Dot(tB, eB2);
      addEMatTransProduct(lam[6], eA, eB2, etaB, epsB, &resB[3], &resB[4]);

      // Add the dummy constraints
      resC[7] += lam[7];
    } else {
      // Add the revolute direction constraint
      TacsScalar tA[3], tB1[3], tB2[3];
      mat3x3MultTrans(CA, eA, tA);
      mat3x3MultTrans(CB, eB1, tB1);
      mat3x3MultTrans(CB, eB2, tB2);

      // Compute the contributions to the first revolute constraint
      resC[3] += vec3Dot(tA, tB1);

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
      addEMatTransProduct(lam[3], eA, tB1, etaA, epsA, &resA[3], &resA[4]);
      // Add the derivative d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA
      addEMatTransProduct(lam[3], eB1, tA, etaB, epsB, &resB[3], &resB[4]);

      // Compute the contributions to the second revolute constraint
      resC[4] += vec3Dot(tA, tB2);

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
      addEMatTransProduct(lam[4], eA, tB2, etaA, epsA, &resA[3], &resA[4]);
      // Add the derivative d(CB^{T}*eB2)/dqB)*tA = E(eB2)*tA
      addEMatTransProduct(lam[4], eB2, tA, etaB, epsB, &resB[3], &resB[4]);

      // Add the dummy constraints
      for (int i = 5; i < 8; i++) {
        resC[i] += lam[i];
      }
    }
  }
}

/*
  Compute the Jacobian of the residuals of the governing equations
*/
void TACSRevoluteConstraint::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *J) {
  // Add the residual
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Get the number of variables
  const int nvars = getNumVariables();

  // Get the initial position vectors - depending on whether we're
  // taking the initial position from the node locations
  const TacsScalar *rA, *rB;
  if (bodyA) {
    // Get the initial position for bodyA
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  } else {
    rA = &Xpts[0];
  }
  if (bodyB) {
    // Get the initial position for bodyA
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  } else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not. Also, set the offset to the Lagrange multipliers
  int offset = 0;
  if (inertial_fixed_point) {
    offset = 8;
    lam = &vars[8];
  } else {
    offset = 16;
    lam = &vars[16];
  }

  // Retrieve the pointers to eA, eB1, eB2
  const TacsScalar *eA, *eB1, *eB2;
  eAVec->getVector(&eA);
  eB1Vec->getVector(&eB1);
  eB2Vec->getVector(&eB2);

  // Compute the rotation matrix
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Add the identity matricies to the Jacobian
  addBlockIdent(alpha, &J[offset], nvars);
  addBlockIdent(alpha, &J[offset * nvars], nvars);

  // Add the second derivative terms
  addBlockDMatTransDeriv(alpha, lam, xA, &J[3 * (nvars + 1)], nvars);

  // Add the term from the derivative w.r.t. lambda
  addBlockEMatTrans(alpha, etaA, epsA, xA, &J[3 * nvars + offset], nvars);

  // Add the term from the derivative of the constraint
  addBlockEMat(alpha, etaA, epsA, xA, &J[offset * nvars + 3], nvars);

  if (inertial_fixed_point) {
    // Add the revolute direction constraint
    TacsScalar tA[3];
    mat3x3MultTrans(CA, eA, tA);

    // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
    TacsScalar gA[4];
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, eB1, etaA, epsA, &gA[0], &gA[1]);

    // Add the constraint terms to the matrix
    for (int i = 0; i < 4; i++) {
      J[nvars * (i + 3) + offset + 3] += gA[i];
      J[(offset + 3) * nvars + i + 3] += gA[i];
    }

    // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, eB2, etaA, epsA, &gA[0], &gA[1]);

    for (int i = 0; i < 4; i++) {
      J[nvars * (i + 3) + offset + 4] += gA[i];
      J[(offset + 4) * nvars + i + 3] += gA[i];
    }

    // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
    addBlockDMatTransDeriv(alpha * lam[3], eB1, eA, &J[3 * (nvars + 1)], nvars);

    // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
    addBlockDMatTransDeriv(alpha * lam[4], eB2, eA, &J[3 * (nvars + 1)], nvars);

    // Add the Jacobian entries for the dummy constraints
    for (int i = offset + 5; i < nvars; i++) {
      J[(nvars + 1) * i] += alpha;
    }
  } else {
    // Compute the rotation matrix
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);

    // Add the block identities
    addBlockIdent(-alpha, &J[8 * nvars + offset], nvars);
    addBlockIdent(-alpha, &J[offset * nvars + 8], nvars);

    // Add the second derivative terms
    addBlockDMatTransDeriv(-alpha, lam, xB, &J[11 * (nvars + 1)], nvars);

    // Add the terms from the derivatives w.r.t. lambdas
    addBlockEMatTrans(-alpha, etaB, epsB, xB, &J[11 * nvars + offset], nvars);

    // Add the terms from the derivatives of the constraint
    addBlockEMat(-alpha, etaB, epsB, xB, &J[offset * nvars + 11], nvars);

    if (inertial_rev_axis) {
      // Add the revolute direction constraint
      TacsScalar tA[3];
      mat3x3MultTrans(CA, eA, tA);

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
      TacsScalar gA[4];
      memset(gA, 0, sizeof(gA));
      addEMatTransProduct(alpha, eA, eB1, etaA, epsA, &gA[0], &gA[1]);

      // Add the constraint terms to the matrix
      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 3) + offset + 3] += gA[i];
        J[(offset + 3) * nvars + i + 3] += gA[i];
      }

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
      memset(gA, 0, sizeof(gA));
      addEMatTransProduct(alpha, eA, eB2, etaA, epsA, &gA[0], &gA[1]);

      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 3) + offset + 4] += gA[i];
        J[(offset + 4) * nvars + i + 3] += gA[i];
      }

      // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
      addBlockDMatTransDeriv(alpha * lam[3], eB1, eA, &J[3 * (nvars + 1)],
                             nvars);

      // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
      addBlockDMatTransDeriv(alpha * lam[4], eB2, eA, &J[3 * (nvars + 1)],
                             nvars);

      // Add the revolute direction constraint
      TacsScalar tB[3];
      mat3x3MultTrans(CB, eA, tB);

      // Add the derivative (d(CB^{T}*eA)/dqB)*tB1 = E(eA)*tB
      TacsScalar gB[4];
      memset(gB, 0, sizeof(gB));
      addEMatTransProduct(alpha, eA, eB1, etaB, epsB, &gB[0], &gB[1]);

      // Add the constraint terms to the matrix
      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 11) + offset + 5] += gB[i];
        J[(offset + 5) * nvars + i + 11] += gB[i];
      }

      // Add the derivative (d(CB^{T}*eA)/dqB)*tB2 = E(eA)*tB2
      memset(gB, 0, sizeof(gB));
      addEMatTransProduct(alpha, eA, eB2, etaB, epsB, &gB[0], &gB[1]);

      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 11) + offset + 6] += gB[i];
        J[(offset + 6) * nvars + i + 11] += gB[i];
      }

      // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
      addBlockDMatTransDeriv(alpha * lam[5], eB1, eA, &J[11 * (nvars + 1)],
                             nvars);

      // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
      addBlockDMatTransDeriv(alpha * lam[6], eB2, eA, &J[11 * (nvars + 1)],
                             nvars);

      // Add the Jacobian entries for the dummy constraints
      J[(nvars + 1) * (nvars - 1)] += alpha;
    } else {
      // Add the revolute direction constraint
      TacsScalar tA[3], tB1[3], tB2[3];
      mat3x3MultTrans(CA, eA, tA);
      mat3x3MultTrans(CB, eB1, tB1);
      mat3x3MultTrans(CB, eB2, tB2);

      TacsScalar gA[4], gB[4];

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
      memset(gA, 0, sizeof(gA));
      addEMatTransProduct(alpha, eA, tB1, etaA, epsA, &gA[0], &gA[1]);
      // Add the derivative d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA
      memset(gB, 0, sizeof(gB));
      addEMatTransProduct(alpha, eB1, tA, etaB, epsB, &gB[0], &gB[1]);

      // Add the constraint terms to the matrix
      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 3) + 19] += gA[i];
        J[nvars * (i + 11) + 19] += gB[i];
        J[19 * nvars + i + 3] += gA[i];
        J[19 * nvars + i + 11] += gB[i];
      }

      // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
      memset(gA, 0, sizeof(gA));
      addEMatTransProduct(alpha, eA, tB2, etaA, epsA, &gA[0], &gA[1]);
      // Add the derivative d(CB^{T}*eB2)/dqB)*tA = E(eB2)*tA
      memset(gB, 0, sizeof(gB));
      addEMatTransProduct(alpha, eB2, tA, etaB, epsB, &gB[0], &gB[1]);

      // Add the constraint terms to the matrix
      for (int i = 0; i < 4; i++) {
        J[nvars * (i + 3) + 20] += gA[i];
        J[nvars * (i + 11) + 20] += gB[i];
        J[20 * nvars + i + 3] += gA[i];
        J[20 * nvars + i + 11] += gB[i];
      }

      // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
      addBlockDMatTransDeriv(alpha * lam[3], tB1, eA, &J[3 * (nvars + 1)],
                             nvars);
      addBlockDMatTransDeriv(alpha * lam[3], tA, eB1, &J[11 * (nvars + 1)],
                             nvars);

      // Add the contributions from the off-diagonal blocks
      TacsScalar EA[12], EB[12];
      computeEMat(etaA, epsA, eA, EA);
      computeEMat(etaB, epsB, eB1, EB);

      addBlock3x4Product(alpha * lam[3], EA, EB, &J[3 * nvars + 11], nvars);
      addBlock3x4Product(alpha * lam[3], EB, EA, &J[11 * nvars + 3], nvars);

      // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
      addBlockDMatTransDeriv(alpha * lam[4], tB2, eA, &J[3 * (nvars + 1)],
                             nvars);
      addBlockDMatTransDeriv(alpha * lam[4], tA, eB2, &J[11 * (nvars + 1)],
                             nvars);

      // Add the contributions from the off-diagonal blocks for the second
      // constraint
      computeEMat(etaB, epsB, eB2, EB);

      addBlock3x4Product(alpha * lam[4], EA, EB, &J[3 * nvars + 11], nvars);
      addBlock3x4Product(alpha * lam[4], EB, EA, &J[11 * nvars + 3], nvars);

      // Add the Jacobian entries for the dummy constraints
      for (int i = offset + 5; i < nvars; i++) {
        J[(nvars + 1) * i] += alpha;
      }
    }
  }
}

/*
  Get the design variable numbers
*/
int TACSRevoluteConstraint::getDesignVarNums(int elemIndex, int dvLen,
                                             int dvNums[]) {
  return point->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the design variable values
*/
int TACSRevoluteConstraint::setDesignVars(int elemIndex, int dvLen,
                                          const TacsScalar dvs[]) {
  int ndvs = point->setDesignVars(elemIndex, dvLen, dvs);
  updatePoints();
  return ndvs;
}

/*
  Get the design variable values associated with the joint location
*/
int TACSRevoluteConstraint::getDesignVars(int elemIndex, int dvLen,
                                          TacsScalar dvs[]) {
  return point->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable values associated with the joint location
*/
int TACSRevoluteConstraint::getDesignVarRange(int elemIndex, int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[]) {
  return point->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Two-point rigid link constraint
*/
TACSRigidLink::TACSRigidLink(TACSRigidBody *_bodyA) {
  bodyA = _bodyA;
  bodyA->incref();
}

TACSRigidLink::~TACSRigidLink() { bodyA->decref(); }

const char *TACSRigidLink::elem_name = "TACSRigidLink";

/*
  Return the number of displacements
*/
int TACSRigidLink::getVarsPerNode() { return 8; }

/*
  Return the number of nodes
*/
int TACSRigidLink::getNumNodes() { return 3; }

/*
  Return the element name
*/
const char *TACSRigidLink::getObjectName() { return elem_name; }

const TacsScalar TACSRigidLink::delta = 20.0;

/*
  Compute the residual of the governing equations
*/
void TACSRigidLink::addResidual(int elemIndex, double time,
                                const TacsScalar *Xpts, const TacsScalar *vars,
                                const TacsScalar *dvars,
                                const TacsScalar *ddvars, TacsScalar *res) {
  // Set pointers to the residual
  TacsScalar *resA = &res[0];
  TacsScalar *resB = &res[8];
  TacsScalar *resC = &res[16];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for point B
  const TacsScalar *uB = &vars[8];
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // Set the pointer to the multipliers
  const TacsScalar *lam = &vars[16];

  // Retrieve the initial position of body A
  TACSGibbsVector *xAVec = bodyA->getInitPosition();
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Read from the node locations, the initial position of body B
  const TacsScalar *xB = &Xpts[3];

  // Compute the rotation matrix for body A
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Compute the distance between body A and the point B in the
  // initial configuration
  TacsScalar t[3];
  t[0] = xA[0] - xB[0];
  t[1] = xA[1] - xB[1];
  t[2] = xA[2] - xB[2];

  // Add the residual
  // resC = uB - uA + (xB - xA) + CA^{T}*(xA - xB)
  // resC = uB - uA + t + CA^{T}*(xA - xB)
  resC[0] += uB[0] - uA[0];
  resC[1] += uB[1] - uA[1];
  resC[2] += uB[2] - uA[2];

  vec3Axpy(-1.0, t, resC);
  mat3x3MultTransAdd(CA, t, resC);

  // Add the residual for the quaternions
  resC[3] += lam[3];
  resC[4] += epsB[0] - epsA[0];
  resC[5] += epsB[1] - epsA[1];
  resC[6] += epsB[2] - epsA[2];

  // Add the dummy constraint for the remaining multiplier
  resC[7] += lam[7];

  // Add the terms from the first constraint
  vec3Axpy(-1.0, &lam[0], &resA[0]);
  addEMatTransProduct(-1.0, t, &lam[0], etaA, epsA, &resA[3], &resA[4]);
  vec3Axpy(1.0, &lam[0], &resB[0]);

  // Add the terms from the second constraint
  vec3Axpy(-1.0, &lam[4], &resA[4]);
  vec3Axpy(1.0, &lam[4], &resB[4]);

  // Add the scalar part to the residual of the linked body
  resA[7] += delta * (etaB - etaA);
  resB[7] += delta * (etaB - etaA);
  resB[3] += delta * (vars[7] + vars[15]);
  resA[3] -= delta * (vars[7] + vars[15]);
}

/*
  Compute the Jacobian of the governing equations
*/
void TACSRigidLink::addJacobian(int elemIndex, double time, TacsScalar alpha,
                                TacsScalar beta, TacsScalar gamma,
                                const TacsScalar *Xpts, const TacsScalar *vars,
                                const TacsScalar *dvars,
                                const TacsScalar *ddvars, TacsScalar *res,
                                TacsScalar *J) {
  // Add the residual
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the Lagrange multiplier variables
  const TacsScalar *lam = &vars[16];

  // The number of variables in the Jacobian matrix
  const int nvars = 3 * 8;

  // Retrieve the initial position of body A
  TACSGibbsVector *xAVec = bodyA->getInitPosition();
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Read from the node locations, the initial position of body B
  const TacsScalar *xB = &Xpts[3];

  // Compute the rotation matrix for body A
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Compute the distance between body A and the point B in the
  // initial configuration
  TacsScalar t[3];
  t[0] = xA[0] - xB[0];
  t[1] = xA[1] - xB[1];
  t[2] = xA[2] - xB[2];

  // Derivatives of the position constraint
  addBlockIdent(-alpha, &J[16 * nvars], nvars);
  addBlockIdent(alpha, &J[16 * nvars + 8], nvars);
  addBlockEMat(alpha, etaA, epsA, t, &J[16 * nvars + 3], nvars);

  // Derivatives of the quaternion constraint
  addBlockIdent(alpha, &J[20 * nvars + 12], nvars);
  addBlockIdent(-alpha, &J[20 * nvars + 4], nvars);

  // Add the Jacobian contribution from the dummy constraint
  J[nvars * nvars - 1] += alpha;

  // Add the contributions from the derivative of resA
  addBlockIdent(-alpha, &J[16], nvars);
  addBlockIdent(alpha, &J[8 * nvars + 16], nvars);
  addBlockDMatTransDeriv(-alpha, lam, t, &J[3 * nvars + 3], nvars);
  addBlockEMatTrans(-alpha, etaA, epsA, t, &J[3 * nvars + 16], nvars);

  // Add the derivatives of the quaternion constraint w.r.t. lam[3]
  J[19 * nvars + 19] += alpha;

  // Add the remaining quaternion constraint derivatives w.r.t. lam[4:]
  addBlockIdent(-alpha, &J[4 * nvars + 20], nvars);
  addBlockIdent(alpha, &J[12 * nvars + 20], nvars);

  // Add the contributions to the residual of the linked body
  J[7 * nvars + 11] += delta * alpha;
  J[11 * nvars + 7] += delta * alpha;
  J[7 * nvars + 3] -= delta * alpha;
  J[3 * nvars + 7] -= delta * alpha;

  J[15 * nvars + 11] += delta * alpha;
  J[11 * nvars + 15] += delta * alpha;
  J[15 * nvars + 3] -= delta * alpha;
  J[3 * nvars + 15] -= delta * alpha;
}

TACSRevoluteDriver::TACSRevoluteDriver(TACSGibbsVector *rev,
                                       TacsScalar _omega) {
  revVec = rev;
  revVec->incref();
  omega = _omega;
}

TACSRevoluteDriver::~TACSRevoluteDriver() { revVec->decref(); }

int TACSRevoluteDriver::getVarsPerNode() { return 8; }

int TACSRevoluteDriver::getNumNodes() { return 2; }

const char *TACSRevoluteDriver::getObjectName() { return "TACSRevoluteDriver"; }

const TacsScalar TACSRevoluteDriver::delta = 10.0;

void TACSRevoluteDriver::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the angle of rotation based on the time
  TacsScalar theta = omega * time;

  // Compute the sin of the half angle
  TacsScalar s = sin(0.5 * theta);

  // Extract the components of the vector
  const TacsScalar *rev;
  revVec->getVector(&rev);

  // Set the multipliers
  const TacsScalar *lam = &vars[8];
  TacsScalar *resC = &res[8];

  // Scale the coefficient
  s *= 1.0 / sqrt(rev[0] * rev[0] + rev[1] * rev[1] + rev[2] * rev[2]);

  // Set the displacements equal to zero
  resC[0] += vars[0];
  resC[1] += vars[1];
  resC[2] += vars[2];

  // Set the quaternion parameters equal to their specified values
  resC[4] += vars[4] - rev[0] * s;
  resC[5] += vars[5] - rev[1] * s;
  resC[6] += vars[6] - rev[2] * s;

  // Add dummy constraints for the remaining multipliers
  resC[3] += lam[3];
  resC[7] += lam[7];

  // Add the multipliers (forces) to the equations of motion
  res[0] += lam[0];
  res[1] += lam[1];
  res[2] += lam[2];

  res[4] += lam[4];
  res[5] += lam[5];
  res[6] += lam[6];

  res[7] += delta * (vars[3] - cos(0.5 * theta));
  res[3] += delta * vars[7];
}

void TACSRevoluteDriver::addJacobian(int elemIndex, double time,
                                     TacsScalar alpha, TacsScalar beta,
                                     TacsScalar gamma, const TacsScalar *Xpts,
                                     const TacsScalar *vars,
                                     const TacsScalar *dvars,
                                     const TacsScalar *ddvars, TacsScalar *res,
                                     TacsScalar *J) {
  // Add the residual to the governing equations
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  const int nvars = 16;

  // Add the terms from the displacement
  J[8 * nvars] += alpha;
  J[9 * nvars + 1] += alpha;
  J[10 * nvars + 2] += alpha;

  // Add terms from the quaternion constraint
  J[12 * nvars + 4] += alpha;
  J[13 * nvars + 5] += alpha;
  J[14 * nvars + 6] += alpha;

  // Add the dummy terms
  J[11 * nvars + 11] += alpha;
  J[15 * nvars + 15] += alpha;

  // Add the multiplier/rigid-body coupling terms
  J[8] += alpha;
  J[nvars + 9] += alpha;
  J[2 * nvars + 10] += alpha;

  J[4 * nvars + 12] += alpha;
  J[5 * nvars + 13] += alpha;
  J[6 * nvars + 14] += alpha;

  // Add the extra constraint
  J[7 * nvars + 3] += delta * alpha;
  J[3 * nvars + 7] += delta * alpha;
}

/*
  Copy the data and increment the reference counts
*/
TACSAverageConstraint::TACSAverageConstraint(TACSRigidBody *_bodyA,
                                             TACSGibbsVector *_point,
                                             TACSRefFrame *_refFrame,
                                             int _moment_flag) {
  moment_flag = _moment_flag;
  bodyA = _bodyA;
  bodyA->incref();
  point = _point;
  point->incref();
  refFrame = _refFrame;
  refFrame->incref();
}

/*
  Free the data
*/
TACSAverageConstraint::~TACSAverageConstraint() {
  bodyA->decref();
  point->decref();
  refFrame->decref();
}

const char *TACSAverageConstraint::elem_name = "TACSAverageConstraint";

/*
  The number of DOF per node
*/
int TACSAverageConstraint::getVarsPerNode() { return 8; }

/*
  The number of nodes for the element
*/
int TACSAverageConstraint::getNumNodes() { return 5; }

/*
  The element name
*/
const char *TACSAverageConstraint::getObjectName() { return elem_name; }

/*
  Add the residual of the governing equations
*/
void TACSAverageConstraint::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Get the initial position for bodyA  from global origin
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Retrieve the reference point location
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the reference point location relative to the intial body
  // point location
  TacsScalar bref[3];
  bref[0] = pt[0] - rA[0];
  bref[1] = pt[1] - rA[1];
  bref[2] = pt[2] - rA[2];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Compute the rotation matrix for the body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // The Lagrange multipliers and constraint pointers
  TacsScalar *con = &res[8 * 4];
  const TacsScalar *lam = &vars[8 * 4];

  // Get the quadrature points and weights
  const int numGauss = 3;
  const double *gaussPts = TacsGaussQuadPts3;
  const double *gaussWts = TacsGaussQuadWts3;

  // Perform the numerical quadrature
  for (int i = 0; i < numGauss; i++) {
    // Get the quadrature point
    const double xi = gaussPts[i];

    // Compute the shape functions
    double N[3];
    N[0] = -0.5 * xi * (1.0 - xi);
    N[1] = (1.0 - xi) * (1.0 + xi);
    N[2] = 0.5 * (1.0 + xi) * xi;

    double Na[3];
    Na[0] = -0.5 + xi;
    Na[1] = -2.0 * xi;
    Na[2] = 0.5 + xi;

    // Compute the position and displacement vector
    TacsScalar Xp[3];
    Xp[0] = N[0] * Xpts[3] + N[1] * Xpts[6] + N[2] * Xpts[9];
    Xp[1] = N[0] * Xpts[4] + N[1] * Xpts[7] + N[2] * Xpts[10];
    Xp[2] = N[0] * Xpts[5] + N[1] * Xpts[8] + N[2] * Xpts[11];

    TacsScalar up[3];
    up[0] = N[0] * vars[8] + N[1] * vars[16] + N[2] * vars[24];
    up[1] = N[0] * vars[9] + N[1] * vars[17] + N[2] * vars[25];
    up[2] = N[0] * vars[10] + N[1] * vars[18] + N[2] * vars[26];

    // Compute the derivative of the position vector along the length
    // of the edge
    TacsScalar Xa[3];
    Xa[0] = Na[0] * Xpts[3] + Na[1] * Xpts[6] + Na[2] * Xpts[9];
    Xa[1] = Na[0] * Xpts[4] + Na[1] * Xpts[7] + Na[2] * Xpts[10];
    Xa[2] = Na[0] * Xpts[5] + Na[1] * Xpts[8] + Na[2] * Xpts[11];

    // Compute the position vector on the element surface relative to
    // the initial reference point
    TacsScalar Xref[3];
    Xref[0] = Xp[0] - rA[0] - bref[0];
    Xref[1] = Xp[1] - rA[1] - bref[1];
    Xref[2] = Xp[2] - rA[2] - bref[2];

    // Evaluate the displacement in the body-fixed coordinate frame:
    TacsScalar atmp[3];
    atmp[0] = Xp[0] + up[0] - uA[0] - rA[0];
    atmp[1] = Xp[1] + up[1] - uA[1] - rA[1];
    atmp[2] = Xp[2] + up[2] - uA[2] - rA[2];

    // uref = CA*(xp + up - uA - rA) - bref - Xref
    TacsScalar uref[3];
    mat3x3Mult(CA, atmp, uref);
    uref[0] = uref[0] - bref[0] - Xref[0];
    uref[1] = uref[1] - bref[1] - Xref[1];
    uref[2] = uref[2] - bref[2] - Xref[2];

    // Compute the quadrature weight for this point
    TacsScalar h = sqrt(vec3Dot(Xa, Xa)) * gaussWts[i];

    // Get the reference axis associated with
    const TacsScalar *Cref;
    refFrame->getRotation(&Cref);

    // Compute the displacements in the local frame
    TacsScalar x[3], u[3];
    mat3x3Mult(Cref, uref, u);

    // Add the integration along the coordinate directions
    con[0] += h * u[0];
    con[1] += h * u[1];
    con[2] += h * u[2];

    // Add the multipliers times the derivative of the constraints
    // w.r.t. the state variables
    for (int j = 0; j < 3; j++) {
      // Iterate over each displacement component and add the
      // contributions from the displacement degrees of freedom.
      for (int k = 0; k < 3; k++) {
        TacsScalar d[3], t[3];
        d[0] = N[j] * CA[k];
        d[1] = N[j] * CA[3 + k];
        d[2] = N[j] * CA[6 + k];
        mat3x3Mult(Cref, d, t);

        // Add the derivative from the elastic degrees of freedom
        res[8 * (j + 1) + k] +=
            h * (t[0] * lam[0] + t[1] * lam[1] + t[2] * lam[2]);

        // Add the contribution from the rigid degrees of freedom
        res[k] -= h * (t[0] * lam[0] + t[1] * lam[1] + t[2] * lam[2]);
      }
    }

    // Add the contributions to the derivative w.r.t. the quaternion
    // parameterization. Compute the transpose of the derivative of
    // (h*lam^{T}*Cref*Cbi*uref)
    TacsScalar s[3];
    mat3x3MultTrans(Cref, &lam[0], s);
    addDMatTransProduct(h, atmp, s, etaA, epsA, &res[3], &res[4]);

    if (moment_flag) {
      // Evaluate the position along the first and second directions
      // in the body-fixed coordinate system
      mat3x3Mult(Cref, Xref, x);

      if (moment_flag & X_MOMENT) {
        con[3] += h * (x[1] * u[2] - x[2] * u[1]);
      }
      if (moment_flag & Y_MOMENT) {
        con[4] += h * x[1] * u[0];
      }
      if (moment_flag & Z_MOMENT) {
        con[5] += h * x[2] * u[0];
      }

      // Add the multipliers times the derivative of the constraints
      // w.r.t. the state variables
      for (int j = 0; j < 3; j++) {
        // Iterate over each displacement component and add the
        // contributions from the displacement degrees of freedom.
        for (int k = 0; k < 3; k++) {
          TacsScalar d[3], t[3];
          d[0] = N[j] * CA[k];
          d[1] = N[j] * CA[3 + k];
          d[2] = N[j] * CA[6 + k];
          mat3x3Mult(Cref, d, t);

          // Add the derivative from the elastic degrees of freedom
          TacsScalar r = 0.0;
          if (moment_flag & X_MOMENT) {
            r += (x[1] * t[2] - x[2] * t[1]) * lam[3];
          }
          if (moment_flag & Y_MOMENT) {
            r += x[1] * t[0] * lam[4];
          }
          if (moment_flag & Z_MOMENT) {
            r += x[2] * t[0] * lam[5];
          }

          // Add the terms to the matrix
          res[8 * (j + 1) + k] += h * r;
          res[k] -= h * r;
        }
      }

      // Compute t = B^{T}*lam
      TacsScalar t[3];
      t[0] = t[1] = t[2] = 0.0;
      if (moment_flag & X_MOMENT) {
        t[1] = -x[2] * lam[3];
        t[2] = x[1] * lam[3];
      }
      if (moment_flag & Y_MOMENT) {
        t[0] = x[1] * lam[4];
      }
      if (moment_flag & Z_MOMENT) {
        t[0] += x[2] * lam[5];
      }
      mat3x3MultTrans(Cref, t, s);
      addDMatTransProduct(h, atmp, s, etaA, epsA, &res[3], &res[4]);

      // Set dummy constraints
      if (!(moment_flag & X_MOMENT)) {
        con[3] += lam[3];
      }
      if (!(moment_flag & Y_MOMENT)) {
        con[4] += lam[4];
      }
      if (!(moment_flag & Z_MOMENT)) {
        con[5] += lam[5];
      }
      con[6] += lam[6];
      con[7] += lam[7];
    } else {
      con[3] += lam[3];
      con[4] += lam[4];
      con[5] += lam[5];
      con[6] += lam[6];
      con[7] += lam[7];
    }
  }
}

void TACSAverageConstraint::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *J) {
  addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  const int nvars = 5 * 8;

  // Get the initial position for bodyA  from global origin
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Retrieve the reference point location
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the reference point location relative to the intial body
  // point location
  TacsScalar bref[3];
  bref[0] = pt[0] - rA[0];
  bref[1] = pt[1] - rA[1];
  bref[2] = pt[2] - rA[2];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Compute the rotation matrix for the body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // The Lagrange multipliers and constraint pointers
  const TacsScalar *lam = &vars[8 * 4];

  // Get the quadrature points and weights
  const int numGauss = 3;
  const double *gaussPts = TacsGaussQuadPts3;
  const double *gaussWts = TacsGaussQuadWts3;

  // Perform the numerical quadrature
  for (int i = 0; i < numGauss; i++) {
    // Get the quadrature point
    const double xi = gaussPts[i];

    // Compute the shape functions
    double N[3];
    N[0] = -0.5 * xi * (1.0 - xi);
    N[1] = (1.0 - xi) * (1.0 + xi);
    N[2] = 0.5 * (1.0 + xi) * xi;

    double Na[3];
    Na[0] = -0.5 + xi;
    Na[1] = -2.0 * xi;
    Na[2] = 0.5 + xi;

    // Compute the position and displacement vector
    TacsScalar Xp[3];
    Xp[0] = N[0] * Xpts[3] + N[1] * Xpts[6] + N[2] * Xpts[9];
    Xp[1] = N[0] * Xpts[4] + N[1] * Xpts[7] + N[2] * Xpts[10];
    Xp[2] = N[0] * Xpts[5] + N[1] * Xpts[8] + N[2] * Xpts[11];

    TacsScalar up[3];
    up[0] = N[0] * vars[8] + N[1] * vars[16] + N[2] * vars[24];
    up[1] = N[0] * vars[9] + N[1] * vars[17] + N[2] * vars[25];
    up[2] = N[0] * vars[10] + N[1] * vars[18] + N[2] * vars[26];

    // Compute the derivative of the position vector along the length
    // of the edge
    TacsScalar Xa[3];
    Xa[0] = Na[0] * Xpts[3] + Na[1] * Xpts[6] + Na[2] * Xpts[9];
    Xa[1] = Na[0] * Xpts[4] + Na[1] * Xpts[7] + Na[2] * Xpts[10];
    Xa[2] = Na[0] * Xpts[5] + Na[1] * Xpts[8] + Na[2] * Xpts[11];

    // Compute the position vector on the element surface relative to
    // the initial reference point
    TacsScalar Xref[3];
    Xref[0] = Xp[0] - rA[0] - bref[0];
    Xref[1] = Xp[1] - rA[1] - bref[1];
    Xref[2] = Xp[2] - rA[2] - bref[2];

    // Evaluate the displacement in the body-fixed coordinate frame:
    TacsScalar atmp[3];
    atmp[0] = Xp[0] + up[0] - uA[0] - rA[0];
    atmp[1] = Xp[1] + up[1] - uA[1] - rA[1];
    atmp[2] = Xp[2] + up[2] - uA[2] - rA[2];

    // uref = CA*(xp + up - uA - rA) - bref - Xref
    TacsScalar uref[3];
    mat3x3Mult(CA, atmp, uref);
    uref[0] = uref[0] - bref[0] - Xref[0];
    uref[1] = uref[1] - bref[1] - Xref[1];
    uref[2] = uref[2] - bref[2] - Xref[2];

    // Compute the quadrature weight for this point
    TacsScalar h = alpha * sqrt(vec3Dot(Xa, Xa)) * gaussWts[i];

    // Get the reference axis associated with
    const TacsScalar *Cref;
    refFrame->getRotation(&Cref);

    // Compute the displacements in the local frame
    TacsScalar u[3];
    mat3x3Mult(Cref, uref, u);

    // Evaluate the position along the first and second directions
    // in the body-fixed coordinate system
    TacsScalar x[3];
    mat3x3Mult(Cref, Xref, x);

    // Compute s = Cref*lambda
    TacsScalar s0[3], s1[3];
    s0[0] = s0[1] = s0[2] = 0.0;
    if (moment_flag & X_MOMENT) {
      s0[1] = -x[2] * lam[3];
      s0[2] = x[1] * lam[3];
    }
    if (moment_flag & Y_MOMENT) {
      s0[0] = x[1] * lam[4];
    }
    if (moment_flag & Z_MOMENT) {
      s0[0] += x[2] * lam[5];
    }
    mat3x3MultTrans(Cref, s0, s1);
    mat3x3MultTrans(Cref, &lam[0], s0);

    // Add the multipliers times the derivative of the constraints
    // w.r.t. the state variables
    for (int j = 0; j < 3; j++) {
      // Iterate over each displacement component and add the
      // contributions from the displacement degrees of freedom.
      for (int k = 0; k < 3; k++) {
        TacsScalar d[3], t[3];
        d[0] = N[j] * CA[k];
        d[1] = N[j] * CA[3 + k];
        d[2] = N[j] * CA[6 + k];
        mat3x3Mult(Cref, d, t);

        // Add the derivative from the elastic degrees of freedom
        J[(8 * (j + 1) + k) * nvars + 32] += h * t[0];
        J[(8 * (j + 1) + k) * nvars + 33] += h * t[1];
        J[(8 * (j + 1) + k) * nvars + 34] += h * t[2];

        J[32 * nvars + (8 * (j + 1) + k)] += h * t[0];
        J[33 * nvars + (8 * (j + 1) + k)] += h * t[1];
        J[34 * nvars + (8 * (j + 1) + k)] += h * t[2];

        // Add the contribution from the rigid degrees of freedom
        J[k * nvars + 32] -= h * t[0];
        J[k * nvars + 33] -= h * t[1];
        J[k * nvars + 34] -= h * t[2];

        J[32 * nvars + k] -= h * t[0];
        J[33 * nvars + k] -= h * t[1];
        J[34 * nvars + k] -= h * t[2];
      }

      // Add the term from the transpose of the derivative
      // of the constraints times the multipliers
      addBlockEMat(h * N[j], etaA, epsA, s0, &J[(8 * (j + 1)) * nvars + 3],
                   nvars);
      addBlockEMatTrans(h * N[j], etaA, epsA, s0, &J[3 * nvars + 8 * (j + 1)],
                        nvars);
    }

    // Add the contributions to the derivative w.r.t. the quaternion
    // parameterization. Compute the transpose of the derivative of
    // (h*lam^{T}*Cref*Cbi*uref)

    // Add the derivative w.r.t. the multipliers
    addBlockDMatTrans(h, etaA, epsA, atmp, &J[3 * nvars + 4 * 8], nvars);

    // Add the derivative w.r.t. the
    addBlockEMat(-h, etaA, epsA, s0, &J[3], nvars);

    // Add the derivative of D(atmp)^{T}*s w.r.t. qA
    addBlockDMatTransDeriv(h, atmp, s0, &J[3 * nvars + 3], nvars);

    // Add the derivative of D(atmp)^{T}*s w.r.t. uA
    addBlockEMatTrans(-h, etaA, epsA, s0, &J[3 * nvars], nvars);

    // Add the derivative of the constraint w.r.t. the quaternions
    addBlockDMat(h, etaA, epsA, atmp, &J[32 * nvars + 3], nvars);

    if (moment_flag) {
      // Add the multipliers times the derivative of the constraints
      // w.r.t. the state variables
      for (int j = 0; j < 3; j++) {
        // Iterate over each displacement component and add the
        // contributions from the displacement degrees of freedom.
        for (int k = 0; k < 3; k++) {
          TacsScalar d[3], t[3];
          d[0] = N[j] * CA[k];
          d[1] = N[j] * CA[3 + k];
          d[2] = N[j] * CA[6 + k];
          mat3x3Mult(Cref, d, t);

          // Add the derivative from the elastic degrees of freedom
          if (moment_flag & X_MOMENT) {
            J[(8 * (j + 1) + k) * nvars + 35] +=
                h * (x[1] * t[2] - x[2] * t[1]);
            J[35 * nvars + (8 * (j + 1) + k)] +=
                h * (x[1] * t[2] - x[2] * t[1]);
            J[k * nvars + 35] -= h * (x[1] * t[2] - x[2] * t[1]);
            J[35 * nvars + k] -= h * (x[1] * t[2] - x[2] * t[1]);
          }
          if (moment_flag & Y_MOMENT) {
            J[(8 * (j + 1) + k) * nvars + 36] += h * x[1] * t[0];
            J[36 * nvars + (8 * (j + 1) + k)] += h * x[1] * t[0];
            J[k * nvars + 36] -= h * x[1] * t[0];
            J[36 * nvars + k] -= h * x[1] * t[0];
          }
          if (moment_flag & Z_MOMENT) {
            J[(8 * (j + 1) + k) * nvars + 37] += h * x[2] * t[0];
            J[37 * nvars + (8 * (j + 1) + k)] += h * x[2] * t[0];
            J[k * nvars + 37] -= h * x[2] * t[0];
            J[37 * nvars + k] -= h * x[2] * t[0];
          }
        }

        // Add the term from the transpose of the derivative
        // of the constraints times the multipliers
        addBlockEMat(h * N[j], etaA, epsA, s1, &J[(8 * (j + 1)) * nvars + 3],
                     nvars);
        addBlockEMatTrans(h * N[j], etaA, epsA, s1, &J[3 * nvars + 8 * (j + 1)],
                          nvars);
      }

      // Add the derivative w.r.t. uA
      addBlockEMatTrans(-h, etaA, epsA, s1, &J[3 * nvars], nvars);

      // Add the derivative contribution to uA
      addBlockDMatTransDeriv(h, atmp, s1, &J[3 * nvars + 3], nvars);

      // Add the derivative w.r.t. the displacement variables
      addBlockEMat(-h, etaA, epsA, s1, &J[3], nvars);

      // Add the derivative w.r.t. the multipliers
      TacsScalar D[12];
      computeDMat(etaA, epsA, atmp, D);
      for (int k = 0; k < 4; k++) {
        if (moment_flag & X_MOMENT) {
          J[(3 + k) * nvars + 35] += h * (x[1] * D[8 + k] - x[2] * D[4 + k]);
          J[35 * nvars + 3 + k] += h * (x[1] * D[8 + k] - x[2] * D[4 + k]);
        }
        if (moment_flag & Y_MOMENT) {
          J[(3 + k) * nvars + 36] += h * x[1] * D[k];
          J[36 * nvars + 3 + k] += h * x[1] * D[k];
        }
        if (moment_flag & Z_MOMENT) {
          J[(3 + k) * nvars + 37] += h * x[2] * D[k];
          J[37 * nvars + 3 + k] += h * x[2] * D[k];
        }
      }

      if (!(moment_flag & X_MOMENT)) {
        J[(nvars - 5) * (nvars + 1)] += alpha;
      }
      if (!(moment_flag & Y_MOMENT)) {
        J[(nvars - 4) * (nvars + 1)] += alpha;
      }
      if (!(moment_flag & Z_MOMENT)) {
        J[(nvars - 3) * (nvars + 1)] += alpha;
      }
      J[(nvars - 2) * (nvars + 1)] += alpha;
      J[(nvars - 1) * (nvars + 1)] += alpha;
    } else {
      for (int k = 3; k < 8; k++) {
        J[(8 * 4 + k) * (nvars + 1)] += alpha;
      }
    }
  }
}

/*
  Construct a fixed constraint with the two bodies involved and a
  position vector measured from the global frame to the point where
  the fixed joint is located.
*/
TACSFixedConstraint::TACSFixedConstraint(TACSRigidBody *_body,
                                         TACSGibbsVector *_point) {
  // Copy over the arguments
  body = _body;
  body->incref();
  point = _point;
  point->incref();

  // Set class variables to NULL
  xVec = NULL;

  updatePoints();
}

/*
  Destructor for fixed constraint
*/
TACSFixedConstraint::~TACSFixedConstraint() {
  body->decref();
  point->decref();
  xVec->decref();
}

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSFixedConstraint::getNumNodes() { return 2; }

const char *TACSFixedConstraint::elem_name = "TACSFixedConstraint";

/*
  Update the local data. This takes the initial reference points for
  the body, given as rA and the initial point for
  the connection between the A and B bodies, "point pt", and computes
  the vectors xA and xB in the global frame.
*/
void TACSFixedConstraint::updatePoints() {
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Fetch the positions of body in global frame
  TACSGibbsVector *rAVec = body->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Determine the position of the joint from bodyA in the global frame
  // xAVec = point - rAVec
  TacsScalar xA[3];
  for (int i = 0; i < 3; i++) {
    xA[i] = pt[i] - rA[i];
  }
  if (xVec) {
    xVec->decref();
  }
  xVec = new TACSGibbsVector(xA);
  xVec->incref();
}

/*
  Compute the residual of the governing equations
*/
void TACSFixedConstraint::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Set pointers to the residual
  TacsScalar *resA = &res[0];
  TacsScalar *resC = &res[8];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar *epsA = &vars[4];

  // Set the pointer to the multipliers
  const TacsScalar *lam = &vars[8];

  // Add the constraint residual equations
  resC[0] += uA[0];
  resC[1] += uA[1];
  resC[2] += uA[2];
  resC[3] += lam[3];
  resC[4] += epsA[0];
  resC[5] += epsA[1];
  resC[6] += epsA[2];

  // Add the dummy constraint for the remaining multiplier
  resC[7] += lam[7];

  // Add the constraint reaction forces/moments to the body residual
  vec3Axpy(1.0, &lam[0], &resA[0]);
  vec3Axpy(1.0, &lam[4], &resA[4]);
}

/*
  Get the design variable numbers
*/
int TACSFixedConstraint::getDesignVarNums(int elemIndex, int dvLen,
                                          int dvNums[]) {
  return point->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the design variable values
*/
int TACSFixedConstraint::setDesignVars(int elemIndex, int dvLen,
                                       const TacsScalar dvs[]) {
  int ndvs = point->setDesignVars(elemIndex, dvLen, dvs);
  updatePoints();
  return ndvs;
}

/*
  Get the design variable values associated with the joint location
*/
int TACSFixedConstraint::getDesignVars(int elemIndex, int dvLen,
                                       TacsScalar dvs[]) {
  return point->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable values associated with the joint location
*/
int TACSFixedConstraint::getDesignVarRange(int elemIndex, int dvLen,
                                           TacsScalar lb[], TacsScalar ub[]) {
  return point->getDesignVarRange(elemIndex, dvLen, lb, ub);
}
