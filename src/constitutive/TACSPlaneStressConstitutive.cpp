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

#include "TACSPlaneStressConstitutive.h"

const char* TACSPlaneStressConstitutive::psName = "TACSPlaneStressConstitutive";

const char* TACSPlaneStressConstitutive::getObjectName() { return psName; }

/*
  PlaneStressStiffness member function definitions
*/
TACSPlaneStressConstitutive::TACSPlaneStressConstitutive(
    TACSMaterialProperties* props, TacsScalar _t, int _tNum, TacsScalar _tlb,
    TacsScalar _tub) {
  properties = props;
  if (properties) {
    properties->incref();
  }
  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
}

TACSPlaneStressConstitutive::~TACSPlaneStressConstitutive() {
  if (properties) {
    properties->decref();
  }
}

int TACSPlaneStressConstitutive::getNumStresses() { return NUM_STRESSES; }

// Retrieve the global design variable numbers
int TACSPlaneStressConstitutive::getDesignVarNums(int elemIndex, int dvLen,
                                                  int dvNums[]) {
  if (tNum >= 0) {
    if (dvNums && dvLen >= 1) {
      dvNums[0] = tNum;
    }
    return 1;
  }
  return 0;
}

// Set the element design variable from the design vector
int TACSPlaneStressConstitutive::setDesignVars(int elemIndex, int dvLen,
                                               const TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSPlaneStressConstitutive::getDesignVars(int elemIndex, int dvLen,
                                               TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSPlaneStressConstitutive::getDesignVarRange(int elemIndex, int dvLen,
                                                   TacsScalar lb[],
                                                   TacsScalar ub[]) {
  if (tNum >= 0 && dvLen >= 1) {
    if (lb) {
      lb[0] = tlb;
    }
    if (ub) {
      ub[0] = tub;
    }
    return 1;
  }
  return 0;
}

// Evaluate the material density
TacsScalar TACSPlaneStressConstitutive::evalDensity(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[]) {
  if (properties) {
    return t * properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSPlaneStressConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    dfdx[0] += scale * properties->getDensity();
  }
}

// Evaluate the specific heat
TacsScalar TACSPlaneStressConstitutive::evalSpecificHeat(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[]) {
  if (properties) {
    return properties->getSpecificHeat();
  }
  return 0.0;
}

// Evaluate the stresss
void TACSPlaneStressConstitutive::evalStress(int elemIndex, const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar e[],
                                             TacsScalar s[]) {
  TacsScalar C[6];
  if (properties) {
    properties->evalTangentStiffness2D(C);

    s[0] = t * (C[0] * e[0] + C[1] * e[1] + C[2] * e[2]);
    s[1] = t * (C[1] * e[0] + C[3] * e[1] + C[4] * e[2]);
    s[2] = t * (C[2] * e[0] + C[4] * e[1] + C[5] * e[2]);
  } else {
    s[0] = s[1] = s[2] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSPlaneStressConstitutive::evalTangentStiffness(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       TacsScalar C[]) {
  if (properties) {
    properties->evalTangentStiffness2D(C);
    C[0] *= t;
    C[1] *= t;
    C[2] *= t;
    C[3] *= t;
    C[4] *= t;
    C[5] *= t;
  } else {
    C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
  }
}

// Add the derivative of the stress w.r.t. design variables
void TACSPlaneStressConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0] * e[0] + C[1] * e[1] + C[2] * e[2]);
    s[1] = (C[1] * e[0] + C[3] * e[1] + C[4] * e[2]);
    s[2] = (C[2] * e[0] + C[4] * e[1] + C[5] * e[2]);

    // Compute the derivative w.r.t. the design vector
    dfdx[0] += scale * (s[0] * psi[0] + s[1] * psi[1] + s[2] * psi[2]);
  }
}

// Evaluate the thermal strain
void TACSPlaneStressConstitutive::evalThermalStrain(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    TacsScalar theta,
                                                    TacsScalar e[]) {
  if (properties) {
    properties->evalThermalStrain2D(e);
    e[0] *= theta;
    e[1] *= theta;
    e[2] *= theta;
  } else {
    e[0] = e[1] = e[2] = 0.0;
  }
}

// Evaluate the heat flux, given the thermal gradient
void TACSPlaneStressConstitutive::evalHeatFlux(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar grad[],
                                               TacsScalar flux[]) {
  if (properties) {
    TacsScalar Kc[3];
    properties->evalTangentHeatFlux2D(Kc);
    flux[0] = t * (Kc[0] * grad[0] + Kc[1] * grad[1]);
    flux[1] = t * (Kc[1] * grad[0] + Kc[2] * grad[1]);
  }
}

// Evaluate the tangent of the heat flux
void TACSPlaneStressConstitutive::evalTangentHeatFlux(int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      TacsScalar Kc[]) {
  if (properties) {
    properties->evalTangentHeatFlux2D(Kc);
    Kc[0] *= t;
    Kc[1] *= t;
    Kc[2] *= t;
  }
}

// Add the derivative of the heat flux
void TACSPlaneStressConstitutive::addHeatFluxDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar grad[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    TacsScalar Kc[3];
    properties->evalTangentHeatFlux2D(Kc);
    dfdx[0] += scale * (psi[0] * (Kc[0] * grad[0] + Kc[1] * grad[1]) +
                        psi[1] * (Kc[1] * grad[0] + Kc[2] * grad[1]));
  }
}

// Evaluate the material failure index
TacsScalar TACSPlaneStressConstitutive::evalFailure(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar e[]) {
  if (properties) {
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0] * e[0] + C[1] * e[1] + C[2] * e[2]);
    s[1] = (C[1] * e[0] + C[3] * e[1] + C[4] * e[2]);
    s[2] = (C[2] * e[0] + C[4] * e[1] + C[5] * e[2]);

    return properties->vonMisesFailure2D(s);
  }
  return 0.0;
}

// Evaluate the derivative of the failure criteria w.r.t. strain
TacsScalar TACSPlaneStressConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar dfde[]) {
  if (properties) {
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0] * e[0] + C[1] * e[1] + C[2] * e[2]);
    s[1] = (C[1] * e[0] + C[3] * e[1] + C[4] * e[2]);
    s[2] = (C[2] * e[0] + C[4] * e[1] + C[5] * e[2]);

    TacsScalar sens[3];
    TacsScalar fail = properties->vonMisesFailure2DStressSens(s, sens);

    dfde[0] = (C[0] * sens[0] + C[1] * sens[1] + C[2] * sens[2]);
    dfde[1] = (C[1] * sens[0] + C[3] * sens[1] + C[4] * sens[2]);
    dfde[2] = (C[2] * sens[0] + C[4] * sens[1] + C[5] * sens[2]);

    return fail;
  }
  return 0.0;
}
