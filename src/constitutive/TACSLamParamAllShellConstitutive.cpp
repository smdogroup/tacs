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

#include "TACSLamParamAllShellConstitutive.h"

#include "TACSElementAlgebra.h"

const char* TACSLamParamAllShellConstitutive::constName = "TACSLamParamAllShellConstitutive";

TACSLamParamAllShellConstitutive::TACSLamParamAllShellConstitutive(TACSOrthotropicPly* _orthoPly,
                                 TacsScalar _t,
                                 int _tNum,
                                 TacsScalar _tlb,
                                 TacsScalar _tub,
                                 int _lpNums[],
                                 TacsScalar _ksWeight) {
  orthoPly = _orthoPly;
  orthoPly->incref();

  t = _t;
  tNum = _tNum;

  numDesignVars = 0;
  if (tNum >= 0) {
    numDesignVars++;
  }

  tlb = _tlb;
  tub = _tub;

  for (int k = 0; k < 6; k++) {
    lp[k] = 0.0;
    lpNums[k] = _lpNums[k];
    if (lpNums[k] >= 0) {
      numDesignVars++;
    }
  }

  ksWeight = _ksWeight;

  kcorr = 5.0 / 6.0;

  numFailAngles = 4;

  // Calculate the invariant properties of the laminate
  TacsScalar Q11, Q12, Q22, Q44, Q55, Q66;
  orthoPly->getLaminateStiffness(&Q11, &Q12, &Q22, &Q44, &Q55, &Q66);

  // The invariant properties for the bending and in-plane stiffness
  U1 = (3.0 * Q11 + 3.0 * Q22 + 2.0 * Q12 + 4.0 * Q66) / 8.0;
  U2 = (Q11 - Q22) / 2.0;
  U3 = (Q11 + Q22 - 2.0 * Q12 - 4.0 * Q66) / 8.0;
  U4 = (Q11 + Q22 + 6.0 * Q12 - 4.0 * Q66) / 8.0;
  U5 = (Q11 + Q22 - 2.0 * Q12 + 4.0 * Q66) / 8.0;

  // The invariant coefficients for the shear coefficients
  U6 = (Q44 + Q55) / 2.0;
  U7 = (Q44 - Q55) / 2.0;
}

TACSLamParamAllShellConstitutive::~TACSLamParamAllShellConstitutive() { orthoPly->decref(); }

// Set the lamination parameter values directly
void TACSLamParamAllShellConstitutive::setLaminationParameters(TacsScalar _lp[]) {
  for (int k = 0; k < 6; k++) {
    lp[k] = _lp[k];
  }
}

// Set the number of failure angles to set
void TACSLamParamAllShellConstitutive::setNumFailAngles(int _numFailAngles) {
  if (_numFailAngles >= 4 && _numFailAngles <= MAX_NUM_FAIL_ANGLES) {
    numFailAngles = _numFailAngles;
  }
}


int TACSLamParamAllShellConstitutive::getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
  int index = 0;
  if (dvNums && dvLen >= numDesignVars) {
    if (tNum >= 0) {
      dvNums[index] = tNum;
      index++;
    }

    for (int k = 0; k < 6; k++) {
      if (lpNums[k] >= 0) {
        dvNums[index] = lpNums[k];
        index++;
      }
    }
  }
  return numDesignVars; // NOTE: Could probably return index instead or compute index automatically in constructor
}

int TACSLamParamAllShellConstitutive::setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
  int index = 0;
  if (tNum >= 0) {
    t = dvs[index];
    index++;
  }
  for (int k = 0; k < 6; k++) {
    if (lpNums[k] >= 0) {
      lp[k] = dvs[index];
      index++;
    }
  }
  return numDesignVars;
}

int TACSLamParamAllShellConstitutive::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  int index = 0;
  if (tNum >= 0) {
    dvs[index] = t;
    index++;
  }
  for (int k = 0; k < 6; k++) {
    if (lpNums[k] >= 0) {
      dvs[index] = lp[k];
      index++;
    }
  }
  return numDesignVars;
}

int TACSLamParamAllShellConstitutive::getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[], TacsScalar ub[]) {
  int index = 0;
  if (tNum >= 0) {
    lb[tNum] = tlb;
    ub[tNum] = tub;
    index++;
  }
  for (int k = 0; k < 6; k++) {
    if (lpNums[k] >= 0) {
      lb[index] = -1.0;
      ub[index] = 1.0;
      index++;
    }
  }
  return numDesignVars;
}

/*!
  Check the determinant of the symmetric 3x3 matrix
  stored in the following format:

  [a[0], a[1], a[2]]
  [a[1], a[3], a[4]]
  [a[2], a[4], a[5]]

  The determinant is:

  a[0]*(a[3]*a[5] - a[4]*a[4]) -
  a[1]*(a[1]*a[5] - a[2]*a[4]) +
  a[2]*(a[1]*a[4] - a[2]*a[3])
*/

int TACSLamParamAllShellConstitutive::checkDeterminant(const TacsScalar a[]) {
  TacsScalar d =
      (a[0] * (a[3] * a[5] - a[4] * a[4]) - a[1] * (a[1] * a[5] - a[2] * a[4]) +
       a[2] * (a[1] * a[4] - a[2] * a[3]));

  if (TacsRealPart(d) <= 0.0) {
    return 0;
  }

  d = a[0] * a[3] - a[1] * a[1];
  if (TacsRealPart(d) <= 0.0) {
    return 0;
  }

  d = a[3] * a[5] - a[4] * a[4];
  if (TacsRealPart(d) <= 0.0) {
    return 0;
  }

  return 1;
}

// Evaluate the mass per unit area
TacsScalar TACSLamParamAllShellConstitutive::evalDensity(int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[]) {
  return t * orthoPly->getDensity();
}

// Add the derivative of the density w.r.t. the design variables
void TACSLamParamAllShellConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  if (tNum >= 0) {
    dfdx[0] += scale * orthoPly->getDensity();
  }
}

// Evaluate the mass moments
void TACSLamParamAllShellConstitutive::evalMassMoments(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    TacsScalar moments[]) {
  TacsScalar rho = orthoPly->getDensity();
  moments[0] = rho * t;
  moments[1] = 0.0;
  moments[2] = rho * t * t * t / 12.0;
}

// Add the sensitivity of the mass moments
void TACSLamParamAllShellConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  if (tNum >= 0) {
    TacsScalar rho = orthoPly->getDensity();
    dfdx[0] += rho * (scale[0] + 0.25 * t * t * scale[2]);
  }
}

// Evaluate the specific heat
TacsScalar TACSLamParamAllShellConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

void TACSLamParamAllShellConstitutive::getStiffness(TacsScalar A[], TacsScalar B[],
                                                 TacsScalar D[],
                                                 TacsScalar As[],
                                                 TacsScalar *drill) {
  // Calculate the in-plane stiffness using the lamination
  // parameters

  TacsScalar V1 = lp[0];  // cos(2*theta)
  TacsScalar V3 = lp[1];  // cos(4*theta)

  A[0] = t * (U1 + U2 * V1 + U3 * V3);
  A[1] = t * (U4 - U3 * V3);
  A[2] = 0.0;
  A[3] = t * (U1 - U2 * V1 + U3 * V3);
  A[4] = 0.0;
  A[5] = t * (U5 - U3 * V3);

  As[0] = kcorr * t * (U6 + U7 * V1);
  As[1] = 0.0;
  As[2] = kcorr * t * (U6 - U7 * V1);

  B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

  TacsScalar W1 = lp[2];  // z**2 * cos(2*theta)
  TacsScalar W2 = lp[3];  // z**2 * sin(2*theta)
  TacsScalar W3 = lp[4];  // z**2 * cos(4*theta)
  TacsScalar W4 = lp[5];  // z**2 * sin(4*theta)

  // Calculate the bending stiffness using the lamination
  // parameters
  TacsScalar p = t * t * t / 12.0;
  D[0] = p * (U1 + U2 * W1 + U3 * W3);
  D[1] = p * (U4 - U3 * W3);
  D[2] = p * (0.5 * U2 * W2 + U3 * W4);
  D[3] = p * (U1 - U2 * W1 + U3 * W3);
  D[4] = p * (0.5 * U2 * W2 - U3 * W4);
  D[5] = p * (U5 - U3 * W3);

  if (!checkDeterminant(A)) {
    fprintf(stderr, "TACSLamParamAllShellConstitutive: Error, A has negative eigenvalues\n");
  }

  if (!checkDeterminant(D)) {
    fprintf(stderr, "TACSLamParamAllShellConstitutive: Error, D has negative eigenvalues\n");
    fprintf(stderr, "W = [%12.6e, %12.6e, %12.6e %12.6e]\n", TacsRealPart(W1), TacsRealPart(W2), TacsRealPart(W3), TacsRealPart(W4));
  }

  *drill = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
}

// Evaluate the stress
void TACSLamParamAllShellConstitutive::evalStress(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar e[],
                                               TacsScalar s[]) {
  TacsScalar A[6], B[6], D[6], As[3], drill;
  getStiffness(A, B, D, As, &drill);
  TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);
}

// Evaluate the derivative of the product of the stress with a vector psi.
void TACSLamParamAllShellConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar sA[6], sD[6];

  if (tNum >= 0) {
    // Calculate the in-plane stiffness using the lamination
    // parameters
    sA[0] = (U1 + U2 * lp[0] + U3 * lp[1]);
    sA[1] = (U4 - U3 * lp[1]);
    sA[2] = 0.0;
    sA[3] = (U1 - U2 * lp[0] + U3 * lp[1]);
    sA[4] = 0.0;
    sA[5] = (U5 - U3 * lp[1]);
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);

    TacsScalar sAs0 = kcorr * (U6 + U7 * lp[0]);
    TacsScalar sAs2 = kcorr * (U6 - U7 * lp[0]);
    dfdx[0] +=
        scale * (sAs0 * psi[6] * e[6] + sAs2 * psi[7] * e[7] +
                 0.5 * DRILLING_REGULARIZATION * (sAs0 + sAs2) * psi[8] * e[8]);

    // Calculate the bending stiffness using the lamination
    // parameters
    TacsScalar p = t * t / 4.0;
    sD[0] = p * (U1 + U2 * lp[2] + U3 * lp[4]);
    sD[1] = p * (U4 - U3 * lp[4]);
    sD[2] = p * (0.5 * U2 * lp[3] + U3 * lp[5]);
    sD[3] = p * (U1 - U2 * lp[2] + U3 * lp[4]);
    sD[4] = p * (0.5 * U2 * lp[3] - U3 * lp[5]);
    sD[5] = p * (U5 - U3 * lp[4]);
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
  if (lpNums[0] >= 0) {
    sA[0] = t * U2;
    sA[1] = sA[2] = 0.0;
    sA[3] = -t * U2;
    sA[4] = sA[5] = 0.0;
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);

    TacsScalar sAs0 = kcorr * t * U7;
    TacsScalar sAs2 = -kcorr * t * U7;
    dfdx[0] +=
        scale * (sAs0 * psi[6] * e[6] + sAs2 * psi[7] * e[7] +
                 0.5 * DRILLING_REGULARIZATION * (sAs0 + sAs2) * psi[8] * e[8]);
    dfdx++;
  }
  if (lpNums[1] >= 0) {
    sA[0] = t * U3;
    sA[1] = -t * U3;
    sA[2] = 0.0;
    sA[3] = t * U3;
    sA[4] = 0.0;
    sA[5] = -t * U3;
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);
    dfdx++;
  }
  if (lpNums[2] >= 0) {
    TacsScalar p = t * t * t / 12.0;
    sD[0] = p * U2;
    sD[1] = sD[2] = 0.0;
    sD[3] = -p * U2;
    sD[4] = sD[5] = 0.0;
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
  if (lpNums[3] >= 0) {
    TacsScalar p = t * t * t / 12.0;
    sD[0] = sD[1] = 0.0;
    sD[2] = 0.5 * p * U2;
    sD[3] = 0.0;
    sD[4] = 0.5 * p * U2;
    sD[5] = 0.0;
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
  if (lpNums[4] >= 0) {
    TacsScalar p = t * t * t / 12.0;
    sD[0] = p * U3;
    sD[1] = -p * U3;
    sD[2] = 0.0;
    sD[3] = p * U3;
    sD[4] = 0.0;
    sD[5] = -p * U3;
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
  if (lpNums[5] >= 0) {
    TacsScalar p = t * t * t / 12.0;
    sD[0] = sD[1] = 0.0;
    sD[2] = p * U3;
    sD[3] = 0.0;
    sD[4] = -p * U3;
    sD[5] = 0.0;
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
}

// Evaluate the tangent stiffness
void TACSLamParamAllShellConstitutive::evalTangentStiffness(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[],
                                                         TacsScalar C[]) {
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];
  TacsScalar *drill = &C[21];
  getStiffness(A, B, D, As, drill);
}


/*!
  Compute the failure loads for each of the ply angles in the
  range

  [-90, -90+d, ... 90-d]

  where d = 180/numFailAngles.

  Note that the calculations are performed using radians.
*/
void TACSLamParamAllShellConstitutive::computeFailure(const TacsScalar strain[],
                                                   TacsScalar fvals[],
                                                   TacsScalar *_max) {
  TacsScalar max = 0.0;
  for (int k = 0; k < numFailAngles; k++) {
    TacsScalar angle = -0.5 * M_PI + k * M_PI / numFailAngles;

    // Compute the strain and failure criteria at the upper surface
    TacsScalar e[3];
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    fvals[2 * k] = orthoPly->failure(angle, e);
    if (TacsRealPart(fvals[2 * k]) > TacsRealPart(max)) {
      max = fvals[2 * k];
    }

    // Compute the failure criteria at the lower surface
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    fvals[2 * k + 1] = orthoPly->failure(angle, e);
    if (TacsRealPart(fvals[2 * k + 1]) > TacsRealPart(max)) {
      max = fvals[2 * k + 1];
    }
  }

  *_max = max;
}

/*!
  Compute the derivative of the failure load w.r.t. the strain and
  accumulate the weighted sensitivity into the array 'sens'
*/
void TACSLamParamAllShellConstitutive::computeFailureStrainSens(
    const TacsScalar strain[], const TacsScalar weights[], TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = 0.0;
  sens[3] = sens[4] = sens[5] = 0.0;
  sens[6] = sens[7] = sens[8] = 0.0;

  for (int k = 0; k < numFailAngles; k++) {
    TacsScalar angle = -0.5 * M_PI + k * M_PI / numFailAngles;

    // Compute the strain and failure criteria at the upper surface
    TacsScalar e[3], e_sens[3];
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    sens[0] += weights[2 * k] * e_sens[0];
    sens[1] += weights[2 * k] * e_sens[1];
    sens[2] += weights[2 * k] * e_sens[2];

    sens[3] += 0.5 * t * weights[2 * k] * e_sens[0];
    sens[4] += 0.5 * t * weights[2 * k] * e_sens[1];
    sens[5] += 0.5 * t * weights[2 * k] * e_sens[2];

    // Compute the strain at the lower surface
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    sens[0] += weights[2 * k + 1] * e_sens[0];
    sens[1] += weights[2 * k + 1] * e_sens[1];
    sens[2] += weights[2 * k + 1] * e_sens[2];

    sens[3] -= 0.5 * t * weights[2 * k + 1] * e_sens[0];
    sens[4] -= 0.5 * t * weights[2 * k + 1] * e_sens[1];
    sens[5] -= 0.5 * t * weights[2 * k + 1] * e_sens[2];
  }
}

/*!
  Compute the derivative of the failure load w.r.t. the plate
  thickness
*/
TacsScalar TACSLamParamAllShellConstitutive::computeFailureDVSens(const TacsScalar strain[], const TacsScalar weights[]) {
  TacsScalar failSens = 0.0;

  for (int k = 0; k < numFailAngles; k++) {
    TacsScalar angle = -0.5 * M_PI + k * M_PI / numFailAngles;

    TacsScalar e[3], e_sens[3];
    // Compute the strain and failure criteria at the upper surface
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    failSens += 0.5 * weights[2 * k] * (strain[3] * e_sens[0] + strain[4] * e_sens[1] + strain[5] * e_sens[2]);

    // Compute the strain at the lower surface
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    failSens -= 0.5 * weights[2 * k + 1] * (strain[3] * e_sens[0] + strain[4] * e_sens[1] + strain[5] * e_sens[2]);
  }

  return failSens;
}

/*!
  Compute the failure load for a series of ply angles and take the
  approximate maximum using the KS function.
*/
TacsScalar TACSLamParamAllShellConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[]) {
  TacsScalar fvals[2 * MAX_NUM_FAIL_ANGLES];
  TacsScalar max, ks_sum = 0.0;

  computeFailure(strain, fvals, &max);
  for (int k = 0; k < 2 * numFailAngles; k++) {
    ks_sum += exp(ksWeight * (fvals[k] - max));
  }

  return max + log(ks_sum) / ksWeight;
}

/*!
  Compute the derivative of the failure load w.r.t. the strain
  values
*/
TacsScalar TACSLamParamAllShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], TacsScalar sens[]) {
  TacsScalar fvals[2 * MAX_NUM_FAIL_ANGLES], weights[2 * MAX_NUM_FAIL_ANGLES];
  TacsScalar max, ks_sum = 0.0;

  computeFailure(strain, fvals, &max);
  for (int k = 0; k < 2 * numFailAngles; k++) {
    ks_sum += exp(ksWeight * (fvals[k] - max));
  }

  for (int k = 0; k < 2 * numFailAngles; k++) {
    weights[k] = exp(ksWeight * (fvals[k] - max)) / ks_sum;
  }

  computeFailureStrainSens(strain, weights, sens);

  return max + log(ks_sum) / ksWeight;
}

/*!
  Functions to determine the derivative of the failure
  load w.r.t. the design variables
*/
void TACSLamParamAllShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  if (tNum >= 0) {
    TacsScalar fvals[2 * MAX_NUM_FAIL_ANGLES], weights[2 * MAX_NUM_FAIL_ANGLES];
    TacsScalar max, ks_sum = 0.0;

    computeFailure(strain, fvals, &max);
    for (int k = 0; k < 2 * numFailAngles; k++) {
      ks_sum += exp(ksWeight * (fvals[k] - max));
    }

    for (int k = 0; k < 2 * numFailAngles; k++) {
      weights[k] = exp(ksWeight * (fvals[k] - max)) / ks_sum;
    }

    dfdx[index] += scale *computeFailureDVSens(strain, weights);
    index++;
  }
}

// Retrieve the design variable for plotting purposes
// --------------------------------------------------
TacsScalar TACSLamParamAllShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return t;
  }
  else if (index <= 6) {
    // Return the bending parameters
    return lp[index + 1];
  }

  return 0.0;
}
