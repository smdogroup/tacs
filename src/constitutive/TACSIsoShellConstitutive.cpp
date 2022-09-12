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

#include "TACSIsoShellConstitutive.h"

#include "TACSElementAlgebra.h"

const char *TACSIsoShellConstitutive::constName = "TACSIsoShellConstitutive";

/*
  Create the shell constitutive object
*/
TACSIsoShellConstitutive::TACSIsoShellConstitutive(
    TACSMaterialProperties *props, TacsScalar _t, int _tNum, TacsScalar _tlb,
    TacsScalar _tub) {
  properties = props;
  if (properties) {
    properties->incref();
  }

  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
  kcorr = 5.0 / 6.0;
  ksWeight = 100.0;
}

TACSIsoShellConstitutive::~TACSIsoShellConstitutive() {
  if (properties) {
    properties->decref();
  }
}

// Retrieve the global design variable numbers
int TACSIsoShellConstitutive::getDesignVarNums(int elemIndex, int dvLen,
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
int TACSIsoShellConstitutive::setDesignVars(int elemIndex, int dvLen,
                                            const TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSIsoShellConstitutive::getDesignVars(int elemIndex, int dvLen,
                                            TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSIsoShellConstitutive::getDesignVarRange(int elemIndex, int dvLen,
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
TacsScalar TACSIsoShellConstitutive::evalDensity(int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[]) {
  if (properties) {
    return t * properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSIsoShellConstitutive::addDensityDVSens(int elemIndex, TacsScalar scale,
                                                const double pt[],
                                                const TacsScalar X[], int dvLen,
                                                TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    dfdx[0] += scale * properties->getDensity();
  }
}

// Evaluate the mass moments
void TACSIsoShellConstitutive::evalMassMoments(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               TacsScalar moments[]) {
  if (properties) {
    TacsScalar rho = properties->getDensity();
    moments[0] = rho * t;
    moments[1] = 0.0;
    moments[2] = rho * t * t * t / 12.0;
  }
}

// Add the sensitivity of the mass moments
void TACSIsoShellConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    TacsScalar rho = properties->getDensity();
    dfdx[0] += rho * (scale[0] + 0.25 * t * t * scale[2]);
  }
}

// Evaluate the specific heat
TacsScalar TACSIsoShellConstitutive::evalSpecificHeat(int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[]) {
  if (properties) {
    return properties->getSpecificHeat();
  }
  return 0.0;
}

// Evaluate the stresss
void TACSIsoShellConstitutive::evalStress(int elemIndex, const double pt[],
                                          const TacsScalar X[],
                                          const TacsScalar e[],
                                          TacsScalar s[]) {
  if (properties) {
    TacsScalar A[6], B[6], D[6], As[3], drill;

    // Compute the tangent stiffness matrix
    properties->evalTangentStiffness2D(A);

    // The bending-stretch coupling matrix is zero in this case
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

    // Scale the in-plane matrix and bending stiffness
    // matrix by the appropriate quantities
    TacsScalar I = t * t * t / 12.0;
    for (int i = 0; i < 6; i++) {
      D[i] = I * A[i];
      A[i] *= t;
    }

    // Set the through-thickness shear stiffness
    As[0] = As[2] = (5.0 / 6.0) * A[5];
    As[1] = 0.0;

    drill = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);

    // Evaluate the stress
    computeStress(A, B, D, As, drill, e, s);
  } else {
    s[0] = s[1] = s[2] = 0.0;
    s[3] = s[4] = s[5] = 0.0;
    s[6] = s[7] = s[8] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSIsoShellConstitutive::evalTangentStiffness(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    TacsScalar C[]) {
  if (properties) {
    TacsScalar *A = &C[0];
    TacsScalar *B = &C[6];
    TacsScalar *D = &C[12];
    TacsScalar *As = &C[18];

    // Compute the tangent stiffness matrix
    properties->evalTangentStiffness2D(A);

    // The bending-stretch coupling matrix is zero in this case
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

    // Scale the in-plane matrix and bending stiffness
    // matrix by the appropriate quantities
    TacsScalar I = t * t * t / 12.0;
    for (int i = 0; i < 6; i++) {
      D[i] = I * A[i];
      A[i] *= t;
    }

    // Set the through-thickness shear stiffness
    As[0] = As[2] = (5.0 / 6.0) * A[5];
    As[1] = 0.0;

    C[21] = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
  } else {
    memset(C, 0, 22 * sizeof(TacsScalar));
  }
}

// Add the contribution
void TACSIsoShellConstitutive::addStressDVSens(int elemIndex, TacsScalar scale,
                                               const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar e[],
                                               const TacsScalar psi[],
                                               int dvLen, TacsScalar dfdx[]) {
  if (properties && tNum >= 0) {
    // Compute the tangent stiffness matrix
    TacsScalar A[6];
    properties->evalTangentStiffness2D(A);

    TacsScalar dI = t * t / 4.0;

    dfdx[0] += scale * (mat3x3SymmInner(A, &psi[0], &e[0]) +
                        dI * mat3x3SymmInner(A, &psi[3], &e[3]) +
                        (5.0 / 6.0) * A[5] *
                            (psi[6] * e[6] + psi[7] * e[7] +
                             DRILLING_REGULARIZATION * psi[8] * e[8]));
  }
}

// Calculate the point-wise failure criteria
TacsScalar TACSIsoShellConstitutive::evalFailure(int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar e[]) {
  if (properties) {
    TacsScalar et[3], eb[3];
    TacsScalar ht = 0.5 * t;

    et[0] = e[0] + ht * e[3];
    et[1] = e[1] + ht * e[4];
    et[2] = e[2] + ht * e[5];

    eb[0] = e[0] - ht * e[3];
    eb[1] = e[1] - ht * e[4];
    eb[2] = e[2] - ht * e[5];

    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar st[3], sb[3];
    mat3x3SymmMult(C, et, st);
    mat3x3SymmMult(C, eb, sb);

    TacsScalar top = properties->vonMisesFailure2D(st);
    TacsScalar bottom = properties->vonMisesFailure2D(sb);

    TacsScalar ksMax;

    if (TacsRealPart(top) > TacsRealPart(bottom)) {
      ksMax = top;
    } else {
      ksMax = bottom;
    }

    // Use a ks approximation for the max value
    TacsScalar ksSum =
        exp(ksWeight * (top - ksMax)) + exp(ksWeight * (bottom - ksMax));
    TacsScalar ksVal = ksMax + log(ksSum) / ksWeight;
    return ksVal;
  }

  return 0.0;
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSIsoShellConstitutive::evalFailureStrainSens(int elemIndex,
                                                           const double pt[],
                                                           const TacsScalar X[],
                                                           const TacsScalar e[],
                                                           TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = 0.0;
  sens[3] = sens[4] = sens[5] = 0.0;
  sens[6] = sens[7] = sens[8] = 0.0;

  if (properties) {
    TacsScalar et[3], eb[3];
    TacsScalar ht = 0.5 * t;

    et[0] = e[0] + ht * e[3];
    et[1] = e[1] + ht * e[4];
    et[2] = e[2] + ht * e[5];

    eb[0] = e[0] - ht * e[3];
    eb[1] = e[1] - ht * e[4];
    eb[2] = e[2] - ht * e[5];

    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar st[3], sb[3];
    mat3x3SymmMult(C, et, st);
    mat3x3SymmMult(C, eb, sb);

    TacsScalar top = properties->vonMisesFailure2D(st);
    TacsScalar bottom = properties->vonMisesFailure2D(sb);

    TacsScalar ksMax;
    if (TacsRealPart(top) > TacsRealPart(bottom)) {
      ksMax = top;
    } else {
      ksMax = bottom;
    }

    // Use a ks approximation for the max value
    TacsScalar ksSum =
        exp(ksWeight * (top - ksMax)) + exp(ksWeight * (bottom - ksMax));
    TacsScalar ksVal = ksMax + log(ksSum) / ksWeight;

    TacsScalar psi[3], phi[3];

    // Contribution from plate top
    properties->vonMisesFailure2DStressSens(st, psi);
    mat3x3SymmMult(C, psi, phi);
    TacsScalar ksFactor = exp(ksWeight * (top - ksMax)) / ksSum;

    sens[0] = ksFactor * phi[0];
    sens[1] = ksFactor * phi[1];
    sens[2] = ksFactor * phi[2];

    sens[3] = ksFactor * ht * phi[0];
    sens[4] = ksFactor * ht * phi[1];
    sens[5] = ksFactor * ht * phi[2];

    // Contribution from plate bottom
    properties->vonMisesFailure2DStressSens(sb, psi);
    mat3x3SymmMult(C, psi, phi);
    ksFactor = exp(ksWeight * (bottom - ksMax)) / ksSum;

    sens[0] += ksFactor * phi[0];
    sens[1] += ksFactor * phi[1];
    sens[2] += ksFactor * phi[2];

    sens[3] -= ksFactor * ht * phi[0];
    sens[4] -= ksFactor * ht * phi[1];
    sens[5] -= ksFactor * ht * phi[2];

    sens[6] = sens[7] = sens[8] = 0.0;

    return ksVal;
  }

  return 0.0;
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSIsoShellConstitutive::addFailureDVSens(int elemIndex, TacsScalar scale,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar e[], int dvLen,
                                                TacsScalar dfdx[]) {
  if (properties && tNum >= 0 && dvLen >= 1) {
    TacsScalar et[3], eb[3];
    TacsScalar ht = 0.5 * t;

    et[0] = e[0] + ht * e[3];
    et[1] = e[1] + ht * e[4];
    et[2] = e[2] + ht * e[5];

    eb[0] = e[0] - ht * e[3];
    eb[1] = e[1] - ht * e[4];
    eb[2] = e[2] - ht * e[5];

    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar st[3], sb[3];
    mat3x3SymmMult(C, et, st);
    mat3x3SymmMult(C, eb, sb);

    TacsScalar top = properties->vonMisesFailure2D(st);
    TacsScalar bottom = properties->vonMisesFailure2D(sb);

    TacsScalar ksMax;
    if (TacsRealPart(top) > TacsRealPart(bottom)) {
      ksMax = top;
    } else {
      ksMax = bottom;
    }
    // Use a ks approximation for the max value
    TacsScalar ksSum =
        exp(ksWeight * (top - ksMax)) + exp(ksWeight * (bottom - ksMax));

    TacsScalar psi[3], phi[3];

    // Contribution from plate top
    properties->vonMisesFailure2DStressSens(st, psi);
    mat3x3SymmMult(C, psi, phi);
    TacsScalar ksFactor = exp(ksWeight * (top - ksMax)) / ksSum;
    dfdx[0] += ksFactor * 0.5 * scale *
               (phi[0] * e[3] + phi[1] * e[4] + phi[2] * e[5]);

    // Contribution from plate bottom
    properties->vonMisesFailure2DStressSens(sb, psi);
    mat3x3SymmMult(C, psi, phi);
    ksFactor = exp(ksWeight * (bottom - ksMax)) / ksSum;
    dfdx[0] -= ksFactor * 0.5 * scale *
               (phi[0] * e[3] + phi[1] * e[4] + phi[2] * e[5]);
  }
}

// Evaluate the thermal strain
void TACSIsoShellConstitutive::evalThermalStrain(int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 TacsScalar theta,
                                                 TacsScalar e[]) {
  if (properties) {
    properties->evalThermalStrain2D(e);
    e[0] *= theta;
    e[1] *= theta;
    e[2] *= theta;

    e[3] = e[4] = e[5] = 0.0;
    e[6] = e[7] = e[8] = 0.0;
  } else {
    e[0] = e[1] = e[2] = 0.0;
    e[3] = e[4] = e[5] = 0.0;
    e[6] = e[7] = e[8] = 0.0;
  }
}

// Evaluate the heat flux, given the thermal gradient
void TACSIsoShellConstitutive::evalHeatFlux(int elemIndex, const double pt[],
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
void TACSIsoShellConstitutive::evalTangentHeatFlux(int elemIndex,
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
void TACSIsoShellConstitutive::addHeatFluxDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar grad[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  if (properties && tNum >= 0 && dvLen >= 1) {
    TacsScalar Kc[3];
    properties->evalTangentHeatFlux2D(Kc);

    dfdx[0] += scale * (psi[0] * (Kc[0] * grad[0] + Kc[1] * grad[1]) +
                        psi[1] * (Kc[1] * grad[0] + Kc[2] * grad[1]));
  }
}

/*
  Return the constitutive name
*/
const char *TACSIsoShellConstitutive::getObjectName() { return constName; }

TacsScalar TACSIsoShellConstitutive::evalDesignFieldValue(int elemIndex,
                                                          const double pt[],
                                                          const TacsScalar X[],
                                                          int index) {
  if (index == 0) {
    return t;
  } else {
    return 0.0;
  }
}