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

#include "TACSSolidConstitutive.h"

const char* TACSSolidConstitutive::sName = "TACSSolidConstitutive";

const char* TACSSolidConstitutive::getObjectName() { return sName; }

/*
  Get the material properties object
*/
TACSMaterialProperties* TACSSolidConstitutive::getMaterialProperties() {
  return properties;
}

/*
  SolidStiffness member function definitions
*/
TACSSolidConstitutive::TACSSolidConstitutive(TACSMaterialProperties* props,
                                             TacsScalar _t, int _tNum,
                                             TacsScalar _tlb, TacsScalar _tub) {
  properties = props;
  if (properties) {
    properties->incref();
  }
  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
}

TACSSolidConstitutive::~TACSSolidConstitutive() {
  if (properties) {
    properties->decref();
  }
}

int TACSSolidConstitutive::getNumStresses() { return NUM_STRESSES; }

// Retrieve the global design variable numbers
int TACSSolidConstitutive::getDesignVarNums(int elemIndex, int dvLen,
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
int TACSSolidConstitutive::setDesignVars(int elemIndex, int dvLen,
                                         const TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSSolidConstitutive::getDesignVars(int elemIndex, int dvLen,
                                         TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSSolidConstitutive::getDesignVarRange(int elemIndex, int dvLen,
                                             TacsScalar lb[], TacsScalar ub[]) {
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
TacsScalar TACSSolidConstitutive::evalDensity(int elemIndex, const double pt[],
                                              const TacsScalar X[]) {
  if (properties) {
    return t * properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSSolidConstitutive::addDensityDVSens(int elemIndex, TacsScalar scale,
                                             const double pt[],
                                             const TacsScalar X[], int dvLen,
                                             TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    dfdx[0] += scale * properties->getDensity();
  }
}

// Evaluate the specific heat
TacsScalar TACSSolidConstitutive::evalSpecificHeat(int elemIndex,
                                                   const double pt[],
                                                   const TacsScalar X[]) {
  if (properties) {
    return t * properties->getSpecificHeat();
  }
  return 0.0;
}

// Add the derivative of the specific heat
void TACSSolidConstitutive::addSpecificHeatDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    dfdx[0] += scale * properties->getSpecificHeat();
  }
}

// Evaluate the stresss
void TACSSolidConstitutive::evalStress(int elemIndex, const double pt[],
                                       const TacsScalar X[],
                                       const TacsScalar e[], TacsScalar s[]) {
  TacsScalar C[21];
  if (properties) {
    properties->evalTangentStiffness3D(C);
    s[0] = t * (C[0] * e[0] + C[1] * e[1] + C[2] * e[2] + C[3] * e[3] +
                C[4] * e[4] + C[5] * e[5]);
    s[1] = t * (C[1] * e[0] + C[6] * e[1] + C[7] * e[2] + C[8] * e[3] +
                C[9] * e[4] + C[10] * e[5]);
    s[2] = t * (C[2] * e[0] + C[7] * e[1] + C[11] * e[2] + C[12] * e[3] +
                C[13] * e[4] + C[14] * e[5]);
    s[3] = t * (C[3] * e[0] + C[8] * e[1] + C[12] * e[2] + C[15] * e[3] +
                C[16] * e[4] + C[17] * e[5]);
    s[4] = t * (C[4] * e[0] + C[9] * e[1] + C[13] * e[2] + C[16] * e[3] +
                C[18] * e[4] + C[19] * e[5]);
    s[5] = t * (C[5] * e[0] + C[10] * e[1] + C[14] * e[2] + C[17] * e[3] +
                C[19] * e[4] + C[20] * e[5]);
  } else {
    s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSSolidConstitutive::evalTangentStiffness(int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 TacsScalar C[]) {
  if (properties) {
    properties->evalTangentStiffness3D(C);
    for (int i = 0; i < 21; i++) {
      C[i] *= t;
    }
  } else {
    memset(C, 0, 21 * sizeof(TacsScalar));
  }
}

// Add the derivative of the stress w.r.t. design variables
void TACSSolidConstitutive::addStressDVSens(int elemIndex, TacsScalar scale,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar e[],
                                            const TacsScalar psi[], int dvLen,
                                            TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    TacsScalar C[21];
    properties->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2] + C[3] * e[3] + C[4] * e[4] +
           C[5] * e[5];
    s[1] = C[1] * e[0] + C[6] * e[1] + C[7] * e[2] + C[8] * e[3] + C[9] * e[4] +
           C[10] * e[5];
    s[2] = C[2] * e[0] + C[7] * e[1] + C[11] * e[2] + C[12] * e[3] +
           C[13] * e[4] + C[14] * e[5];
    s[3] = C[3] * e[0] + C[8] * e[1] + C[12] * e[2] + C[15] * e[3] +
           C[16] * e[4] + C[17] * e[5];
    s[4] = C[4] * e[0] + C[9] * e[1] + C[13] * e[2] + C[16] * e[3] +
           C[18] * e[4] + C[19] * e[5];
    s[5] = C[5] * e[0] + C[10] * e[1] + C[14] * e[2] + C[17] * e[3] +
           C[19] * e[4] + C[20] * e[5];

    // Compute the derivative w.r.t. the design vector
    dfdx[0] += scale * (s[0] * psi[0] + s[1] * psi[1] + s[2] * psi[2] +
                        s[3] * psi[3] + s[4] * psi[4] + s[5] * psi[5]);
  }
}

// Evaluate the thermal strain
void TACSSolidConstitutive::evalThermalStrain(int elemIndex, const double pt[],
                                              const TacsScalar X[],
                                              TacsScalar theta,
                                              TacsScalar e[]) {
  if (properties) {
    properties->evalThermalStrain3D(e);
    e[0] *= theta;
    e[1] *= theta;
    e[2] *= theta;
    e[3] *= theta;
    e[4] *= theta;
    e[5] *= theta;
  } else {
    e[0] = e[1] = e[2] = e[3] = e[4] = e[5] = 0.0;
  }
}

// Evaluate the heat flux, given the thermal gradient
void TACSSolidConstitutive::evalHeatFlux(int elemIndex, const double pt[],
                                         const TacsScalar X[],
                                         const TacsScalar grad[],
                                         TacsScalar flux[]) {
  if (properties) {
    TacsScalar C[6];
    properties->evalTangentHeatFlux3D(C);
    flux[0] = t * (C[0] * grad[0] + C[1] * grad[1] + C[2] * grad[2]);
    flux[1] = t * (C[1] * grad[0] + C[3] * grad[1] + C[4] * grad[2]);
    flux[2] = t * (C[2] * grad[0] + C[4] * grad[1] + C[5] * grad[2]);
  }
}

// Evaluate the tangent of the heat flux
void TACSSolidConstitutive::evalTangentHeatFlux(int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                TacsScalar C[]) {
  if (properties) {
    properties->evalTangentHeatFlux3D(C);
    for (int i = 0; i < 6; i++) {
      C[i] *= t;
    }
  }
}

// Add the derivative of the heat flux
void TACSSolidConstitutive::addHeatFluxDVSens(int elemIndex, TacsScalar scale,
                                              const double pt[],
                                              const TacsScalar X[],
                                              const TacsScalar grad[],
                                              const TacsScalar psi[], int dvLen,
                                              TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    TacsScalar C[6];
    properties->evalTangentHeatFlux3D(C);
    dfdx[0] +=
        scale * (psi[0] * (C[0] * grad[0] + C[1] * grad[1] + C[2] * grad[2]) +
                 psi[1] * (C[1] * grad[0] + C[3] * grad[1] + C[4] * grad[2]) +
                 psi[2] * (C[2] * grad[0] + C[4] * grad[1] + C[5] * grad[2]));
  }
}

// Evaluate the material failure index
TacsScalar TACSSolidConstitutive::evalFailure(int elemIndex, const double pt[],
                                              const TacsScalar X[],
                                              const TacsScalar e[]) {
  if (properties) {
    TacsScalar C[21];
    properties->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2] + C[3] * e[3] + C[4] * e[4] +
           C[5] * e[5];
    s[1] = C[1] * e[0] + C[6] * e[1] + C[7] * e[2] + C[8] * e[3] + C[9] * e[4] +
           C[10] * e[5];
    s[2] = C[2] * e[0] + C[7] * e[1] + C[11] * e[2] + C[12] * e[3] +
           C[13] * e[4] + C[14] * e[5];
    s[3] = C[3] * e[0] + C[8] * e[1] + C[12] * e[2] + C[15] * e[3] +
           C[16] * e[4] + C[17] * e[5];
    s[4] = C[4] * e[0] + C[9] * e[1] + C[13] * e[2] + C[16] * e[3] +
           C[18] * e[4] + C[19] * e[5];
    s[5] = C[5] * e[0] + C[10] * e[1] + C[14] * e[2] + C[17] * e[3] +
           C[19] * e[4] + C[20] * e[5];

    return properties->vonMisesFailure3D(s);
  }
  return 0.0;
}

// Evaluate the derivative of the failure criteria w.r.t. strain
TacsScalar TACSSolidConstitutive::evalFailureStrainSens(int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        const TacsScalar e[],
                                                        TacsScalar dfde[]) {
  if (properties) {
    TacsScalar C[21];
    properties->evalTangentStiffness3D(C);

    TacsScalar s[6];
    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2] + C[3] * e[3] + C[4] * e[4] +
           C[5] * e[5];
    s[1] = C[1] * e[0] + C[6] * e[1] + C[7] * e[2] + C[8] * e[3] + C[9] * e[4] +
           C[10] * e[5];
    s[2] = C[2] * e[0] + C[7] * e[1] + C[11] * e[2] + C[12] * e[3] +
           C[13] * e[4] + C[14] * e[5];
    s[3] = C[3] * e[0] + C[8] * e[1] + C[12] * e[2] + C[15] * e[3] +
           C[16] * e[4] + C[17] * e[5];
    s[4] = C[4] * e[0] + C[9] * e[1] + C[13] * e[2] + C[16] * e[3] +
           C[18] * e[4] + C[19] * e[5];
    s[5] = C[5] * e[0] + C[10] * e[1] + C[14] * e[2] + C[17] * e[3] +
           C[19] * e[4] + C[20] * e[5];

    TacsScalar sens[6];
    TacsScalar fail = properties->vonMisesFailure3DStressSens(s, sens);

    dfde[0] = C[0] * sens[0] + C[1] * sens[1] + C[2] * sens[2] +
              C[3] * sens[3] + C[4] * sens[4] + C[5] * sens[5];
    dfde[1] = C[1] * sens[0] + C[6] * sens[1] + C[7] * sens[2] +
              C[8] * sens[3] + C[9] * sens[4] + C[10] * sens[5];
    dfde[2] = C[2] * sens[0] + C[7] * sens[1] + C[11] * sens[2] +
              C[12] * sens[3] + C[13] * sens[4] + C[14] * sens[5];
    dfde[3] = C[3] * sens[0] + C[8] * sens[1] + C[12] * sens[2] +
              C[15] * sens[3] + C[16] * sens[4] + C[17] * sens[5];
    dfde[4] = C[4] * sens[0] + C[9] * sens[1] + C[13] * sens[2] +
              C[16] * sens[3] + C[18] * sens[4] + C[19] * sens[5];
    dfde[5] = C[5] * sens[0] + C[10] * sens[1] + C[14] * sens[2] +
              C[17] * sens[3] + C[19] * sens[4] + C[20] * sens[5];

    return fail;
  }
  return 0.0;
}
