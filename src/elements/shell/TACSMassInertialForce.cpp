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

#include "TACSMassInertialForce.h"

TACSMassInertialForce::TACSMassInertialForce(TACSGeneralMassConstitutive *_con,
                                             const TacsScalar _inertiaVec[],
                                             const int *_inertiaVecDVNums) {
  con = _con;
  con->incref();
  memset(inertiaVec, 0, NUM_DISPS * sizeof(TacsScalar));
  memcpy(inertiaVec, _inertiaVec, 3 * sizeof(TacsScalar));
  if (_inertiaVecDVNums) {
    memcpy(inertiaVecDVNums, _inertiaVecDVNums, 3 * sizeof(int));
  } else {
    inertiaVecDVNums[0] = inertiaVecDVNums[1] = inertiaVecDVNums[2] = -1;
  }
}

TACSMassInertialForce::~TACSMassInertialForce() { con->decref(); }

const char *TACSMassInertialForce::getObjectName() {
  return "TACSMassInertialForce";
}

/*
  Retrieve the numbers of displacements, nodes, stress and variables
*/
int TACSMassInertialForce::getVarsPerNode() { return NUM_DISPS; }

int TACSMassInertialForce::getNumNodes() { return NUM_NODES; }

/*
  Add the residual to the provided vector
*/
void TACSMassInertialForce::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  TacsScalar f[NUM_DISPS];
  double pt[3] = {0.0, 0.0, 0.0};
  con->evalInertia(elemIndex, pt, Xpts, inertiaVec, f);
  for (int i = 0; i < NUM_DISPS; i++) {
    res[i] -= f[i];
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSMassInertialForce::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res,
    TacsScalar *mat) {
  if (res) {
    addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  }
}

void TACSMassInertialForce::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  double pt[3] = {0.0, 0.0, 0.0};
  // Add the product of the derivative of the inertial force w.r.t. constitutive
  // DVs
  con->addInertiaDVSens(elemIndex, -scale, pt, Xpts, inertiaVec, psi, dvLen,
                        dfdx);
  // Add sensitivity w.r.t. inertia vector DVs.
  // res[i] -= f[i] where f = M * inertiaVec, so d(res)/d(inertiaVec[j]) =
  // -M[:,j]. The adjoint product is: -scale * psi^T * M[:,j] = -scale * (M *
  // psi)[j], which equals -scale * evalInertia(psi)[j] by symmetry of M.
  TacsScalar g[NUM_DISPS];
  con->evalInertia(elemIndex, pt, Xpts, psi, g);
  for (int i = 0; i < 3; i++) {
    if (inertiaVecDVNums[i] >= 0 && inertiaVecDVNums[i] < dvLen) {
      dfdx[inertiaVecDVNums[i]] += -scale * g[i];
    }
  }
}
