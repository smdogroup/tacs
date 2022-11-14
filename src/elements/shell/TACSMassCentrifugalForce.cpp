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

#include "TACSMassCentrifugalForce.h"

#include "TACSElementAlgebra.h"

TACSMassCentrifugalForce::TACSMassCentrifugalForce(
    TACSGeneralMassConstitutive *_con, const TacsScalar _omegaVec[],
    const TacsScalar _rotCenter[]) {
  con = _con;
  con->incref();
  memcpy(omegaVec, _omegaVec, 3 * sizeof(TacsScalar));
  memcpy(rotCenter, _rotCenter, 3 * sizeof(TacsScalar));
}

TACSMassCentrifugalForce::~TACSMassCentrifugalForce() { con->decref(); }

const char *TACSMassCentrifugalForce::getObjectName() {
  return "TACSMassCentrifugalForce";
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSMassCentrifugalForce::getVarsPerNode() { return NUM_DISPS; }

int TACSMassCentrifugalForce::getNumNodes() { return NUM_NODES; }

/*
  Add the residual to the provided vector
*/
void TACSMassCentrifugalForce::addResidual(
    int elemIndex, double time, const TacsScalar *X, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  double pt[3] = {0.0, 0.0, 0.0};

  TacsScalar r[3], wxr[3];
  TacsScalar ac[NUM_DISPS], f[NUM_DISPS];
  memset(ac, 0, NUM_DISPS * sizeof(TacsScalar));

  // Create vector pointing from rotation center to element gpt
  r[0] = X[0] - rotCenter[0];
  r[1] = X[1] - rotCenter[1];
  r[2] = X[2] - rotCenter[2];

  // Compute omega x r
  crossProduct(omegaVec, r, wxr);

  // Compute centrifugal acceleration
  crossProduct(omegaVec, wxr, ac);

  con->evalInertia(elemIndex, pt, X, ac, f);

  for (int i = 0; i < NUM_DISPS; i++) {
    res[i] += f[i];
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSMassCentrifugalForce::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *X, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *mat) {
  if (res) {
    addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
  }
}
