/*
  This file is part of the package TACS.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include "TACSLamParamShellConstitutive.h"

#include "TACSElementAlgebra.h"

const char *TACSLamParamShellConstitutive::constName =
    "TACSLamParamShellConstitutive";

TACSLamParamShellConstitutive::TACSLamParamShellConstitutive(
    TACSOrthotropicPly *_orthoPly, TacsScalar _t, int _t_num, TacsScalar _min_t,
    TacsScalar _max_t, TacsScalar _f0, TacsScalar _f45, TacsScalar _f90,
    int _f0_num, int _f45_num, int _f90_num, TacsScalar _min_f0,
    TacsScalar _min_f45, TacsScalar _min_f90, TacsScalar _W1, TacsScalar _W3,
    int _W1_num, int _W3_num, TacsScalar _ksWeight, TacsScalar _epsilon) {
  orthoPly = _orthoPly;
  orthoPly->incref();

  t = _t;
  tNum = _t_num;

  numDesignVars = 0;
  if (tNum >= 0) {
    numDesignVars++;
  }

  tlb = _min_t;
  tub = _max_t;

  // Record the fraction dv numbers
  n0 = _f0_num;
  n45 = _f45_num;
  n90 = _f90_num;

  // Store the ply fractions at 0, 45, 90
  f0 = _f0;
  f45 = _f45;
  f90 = _f90;

  // Minimum bound on the ply fraction
  min_f0 = _min_f0;
  min_f45 = _min_f45;
  min_f90 = _min_f90;

  // Adjust it to something reasonable
  if (TacsRealPart(min_f0) < 0.0) {
    min_f0 = 0.0;
  }
  if (TacsRealPart(min_f45) < 0.0) {
    min_f45 = 0.0;
  }
  if (TacsRealPart(min_f90) < 0.0) {
    min_f90 = 0.0;
  }

  if ((n0 >= 0 || n45 >= 0 || n90 >= 0) && (n0 < 0 || n45 < 0 || n90 < 0)) {
    fprintf(stderr,
            "TACSLamParamShellConstitutive: Either all ply fractions "
            "must be dvs, or none can be\n");
  }
  if (n0 >= 0) {
    numDesignVars += 3;
  }

  // Record the lamination parameters
  nW1 = _W1_num;
  nW3 = _W3_num;

  W1 = _W1;
  W3 = _W3;

  if ((nW1 >= 0 || nW3 >= 0) && (nW1 < 0 || nW3 < 0)) {
    fprintf(stderr,
            "TACSLamParamShellConstitutive: Either all lamination "
            "parameters must be dvs, or none can be\n");
  }
  if (nW1 >= 0) {
    numDesignVars += 2;
  }

  ksWeight = _ksWeight;
  epsilon = _epsilon;
  if (TacsRealPart(epsilon) < 0.0) {
    epsilon = 0.0;
  }

  kcorr = 5.0 / 6.0;

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

TACSLamParamShellConstitutive::~TACSLamParamShellConstitutive() {
  orthoPly->decref();
}

int TACSLamParamShellConstitutive::getDesignVarNums(int elemIndex, int dvLen,
                                                    int dvNums[]) {
  if (dvNums && dvLen >= numDesignVars) {
    int i = 0;
    if (tNum >= 0) {
      dvNums[i] = tNum;
      i++;
    }
    if (n0 >= 0) {
      dvNums[i] = n0;
      i++;
    }
    if (n45 >= 0) {
      dvNums[i] = n45;
      i++;
    }
    if (n90 >= 0) {
      dvNums[i] = n90;
      i++;
    }
    if (nW1 >= 0) {
      dvNums[i] = nW1;
      i++;
    }
    if (nW3 >= 0) {
      dvNums[i] = nW3;
      i++;
    }
  }
  return numDesignVars;
}

int TACSLamParamShellConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                 const TacsScalar dvs[]) {
  int i = 0;
  if (tNum >= 0) {
    t = dvs[i];
    i++;
  }
  if (n0 >= 0) {
    f0 = dvs[i];
    i++;
  }
  if (n45 >= 0) {
    f45 = dvs[i];
    i++;
  }
  if (n90 >= 0) {
    f90 = dvs[i];
    i++;
  }
  if (nW1 >= 0) {
    W1 = dvs[i];
    i++;
  }
  if (nW3 >= 0) {
    W3 = dvs[i];
    i++;
  }
  return numDesignVars;
}

int TACSLamParamShellConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                 TacsScalar dvs[]) {
  int i = 0;
  if (tNum >= 0) {
    dvs[i] = t;
    i++;
  }
  if (n0 >= 0) {
    dvs[i] = f0;
    i++;
  }
  if (n45 >= 0) {
    dvs[i] = f45;
    i++;
  }
  if (n90 >= 0) {
    dvs[i] = f90;
    i++;
  }
  if (nW1 >= 0) {
    dvs[i] = W1;
    i++;
  }
  if (nW3 >= 0) {
    dvs[i] = W3;
    i++;
  }
  return numDesignVars;
}

int TACSLamParamShellConstitutive::getDesignVarRange(int elemIndex, int dvLen,
                                                     TacsScalar lb[],
                                                     TacsScalar ub[]) {
  int i = 0;
  if (tNum >= 0) {
    lb[i] = tlb;
    ub[i] = tub;
    i++;
  }
  if (n0 >= 0) {
    lb[i] = min_f0;
    ub[i] = 1.0;
    i++;
  }
  if (n45 >= 0) {
    lb[i] = min_f45;
    ub[i] = 1.0;
    i++;
  }
  if (n90 >= 0) {
    lb[i] = min_f90;
    ub[i] = 1.0;
    i++;
  }
  if (nW1 >= 0) {
    lb[i] = -1.0;
    ub[i] = 1.0;
    i++;
  }
  if (nW3 >= 0) {
    lb[i] = -1.0;
    ub[i] = 1.0;
    i++;
  }
  return numDesignVars;
}

/*
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
int TACSLamParamShellConstitutive::checkDeterminant(const TacsScalar a[]) {
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
TacsScalar TACSLamParamShellConstitutive::evalDensity(int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[]) {
  return t * orthoPly->getDensity();
}

// Add the derivative of the density w.r.t. the design variables
void TACSLamParamShellConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  if (tNum >= 0) {
    dfdx[0] += scale * orthoPly->getDensity();
  }
}

// Evaluate the mass moments
void TACSLamParamShellConstitutive::evalMassMoments(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    TacsScalar moments[]) {
  TacsScalar rho = orthoPly->getDensity();
  moments[0] = rho * t;
  moments[1] = 0.0;
  moments[2] = rho * t * t * t / 12.0;
}

// Add the sensitivity of the mass moments
void TACSLamParamShellConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  if (tNum >= 0) {
    TacsScalar rho = orthoPly->getDensity();
    dfdx[0] += rho * (scale[0] + 0.25 * t * t * scale[2]);
  }
}

// Evaluate the specific heat
TacsScalar TACSLamParamShellConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

// Get the stiffness values
void TACSLamParamShellConstitutive::getStiffness(TacsScalar A[], TacsScalar B[],
                                                 TacsScalar D[],
                                                 TacsScalar As[],
                                                 TacsScalar *drill) {
  // Calculate the in-plane stiffness using the lamination
  // parameters
  TacsScalar V1 = f0 - f90;
  TacsScalar V3 = f0 + f90 - f45;

  A[0] = t * (U1 + U2 * V1 + U3 * V3);
  A[1] = t * (U4 - U3 * V3);
  A[3] = t * (U1 - U2 * V1 + U3 * V3);
  A[5] = t * (U5 - U3 * V3);
  A[2] = A[4] = 0.0;

  As[0] = kcorr * t * (U6 + U7 * V1);
  As[2] = kcorr * t * (U6 - U7 * V1);
  As[1] = 0.0;

  // The bending-stretching coupling matrix is zero
  B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

  // Calculate the bending stiffness using the lamination
  // parameters
  TacsScalar p = t * t * t / 12.0;
  D[0] = p * (U1 + U2 * W1 + U3 * W3);
  D[1] = p * (U4 - U3 * W3);
  D[3] = p * (U1 - U2 * W1 + U3 * W3);
  D[5] = p * (U5 - U3 * W3);
  D[4] = D[2] = 0.0;

  if (!checkDeterminant(A)) {
    fprintf(
        stderr,
        "TACSLamParamShellConstitutive: Error, A has negative eigenvalues\n");
  }

  if (!checkDeterminant(D)) {
    fprintf(
        stderr,
        "TACSLamParamShellConstitutive: Error, D has negative eigenvalues\n");
    fprintf(stderr, "W = [%12.6e, %12.6e, %12.6e %12.6e]\n", TacsRealPart(W1),
            0.0, TacsRealPart(W3), 0.0);
  }

  *drill = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
}

// Evaluate the stress
void TACSLamParamShellConstitutive::evalStress(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar e[],
                                               TacsScalar s[]) {
  TacsScalar A[6], B[6], D[6], As[3], drill;
  getStiffness(A, B, D, As, &drill);
  TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);
}

// Evaluate the tangent stiffness
void TACSLamParamShellConstitutive::evalTangentStiffness(int elemIndex,
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

// Evaluate the derivative of the product of the stress with a vector
void TACSLamParamShellConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar sA[6], sD[6];

  if (tNum >= 0) {
    TacsScalar V1 = f0 - f90;
    TacsScalar V3 = f0 + f90 - f45;
    sA[0] = (U1 + U2 * V1 + U3 * V3);
    sA[1] = (U4 - U3 * V3);
    sA[2] = 0.0;
    sA[3] = (U1 - U2 * V1 + U3 * V3);
    sA[4] = 0.0;
    sA[5] = (U5 - U3 * V3);
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);

    TacsScalar sAs0 = kcorr * (U6 + U7 * V1);
    TacsScalar sAs2 = kcorr * (U6 - U7 * V1);
    dfdx[0] +=
        scale * (sAs0 * psi[6] * e[6] + sAs2 * psi[7] * e[7] +
                 0.5 * DRILLING_REGULARIZATION * (sAs0 + sAs2) * psi[8] * e[8]);

    // Calculate the bending stiffness using the lamination
    // parameters
    TacsScalar p = t * t / 4.0;
    sD[0] = p * (U1 + U2 * W1 + U3 * W3);
    sD[1] = p * (U4 - U3 * W3);
    sD[2] = 0.0;
    sD[3] = p * (U1 - U2 * W1 + U3 * W3);
    sD[4] = 0.0;
    sD[5] = p * (U5 - U3 * W3);
    dfdx[0] += scale * mat3x3SymmInner(sD, &psi[3], &e[3]);
    dfdx++;
  }
  if (n0 >= 0) {
    sA[0] = t * (U2 + U3);
    sA[1] = -t * U3;
    sA[2] = 0.0;
    sA[3] = t * (-U2 + U3);
    sA[4] = 0.0;
    sA[5] = -t * U3;
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);

    TacsScalar sAs0 = t * kcorr * U7;
    TacsScalar sAs2 = -t * kcorr * U7;
    dfdx[0] +=
        scale * (sAs0 * psi[6] * e[6] + sAs2 * psi[7] * e[7] +
                 0.5 * DRILLING_REGULARIZATION * (sAs0 + sAs2) * psi[8] * e[8]);
    dfdx++;
  }
  if (n45 >= 0) {
    sA[0] = -t * U3;
    sA[1] = t * U3;
    sA[2] = 0.0;
    sA[3] = -t * U3;
    sA[4] = 0.0;
    sA[5] = t * U3;

    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);
    dfdx++;
  }
  if (n90 >= 0) {
    sA[0] = t * (-U2 + U3);
    sA[1] = -t * U3;
    sA[2] = 0.0;
    sA[3] = t * (U2 + U3);
    sA[4] = 0.0;
    sA[5] = -t * U3;
    dfdx[0] += scale * mat3x3SymmInner(sA, &psi[0], &e[0]);

    TacsScalar sAs0 = -t * kcorr * U7;
    TacsScalar sAs2 = t * kcorr * U7;
    dfdx[0] +=
        scale * (sAs0 * psi[6] * e[6] + sAs2 * psi[7] * e[7] +
                 0.5 * DRILLING_REGULARIZATION * (sAs0 + sAs2) * psi[8] * e[8]);
    dfdx++;
  }
  if (nW1 >= 0) {
    TacsScalar p = t * t * t / 12.0;
    sD[0] = p * U2;
    sD[3] = -p * U2;
    dfdx[0] += scale * (sD[0] * psi[3] * e[3] + sD[3] * psi[4] * e[4]);
    dfdx++;
  }
  if (nW3 >= 0) {
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
}

/*!
  Compute the failure loads for each of the ply angles in the
  range

  [-90, -90+d, ... 90-d]

  where d = 180/numFailAngles.

  Note that the calculations are performed using radians.
*/
void TACSLamParamShellConstitutive::computeFailure(const TacsScalar strain[],
                                                   TacsScalar fvals[],
                                                   TacsScalar *_max) {
  TacsScalar max = 0.0;
  for (int k = 0; k < NUM_FAIL_ANGLES; k++) {
    TacsScalar angle = 0.0, factor = 1.0;
    if (k == 0) {
      angle = 0.0;
      factor = f0 / (f0 + epsilon);
    } else if (k == 1) {
      angle = 0.25 * M_PI;
      factor = f45 / (f45 + epsilon);
    } else if (k == 2) {
      angle = -0.25 * M_PI;
      factor = f45 / (f45 + epsilon);
    } else {
      angle = 0.5 * M_PI;
      factor = f90 / (f90 + epsilon);
    }

    // Compute the strain and failure criteria at the upper surface
    TacsScalar e[3];
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    fvals[2 * k] = factor * orthoPly->failure(angle, e);
    if (TacsRealPart(fvals[2 * k]) > TacsRealPart(max)) {
      max = fvals[2 * k];
    }

    // Compute the failure criteria at the lower surface
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    fvals[2 * k + 1] = factor * orthoPly->failure(angle, e);
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
void TACSLamParamShellConstitutive::computeFailureStrainSens(
    const TacsScalar strain[], const TacsScalar weights[], TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = sens[3] = 0.0;
  sens[4] = sens[5] = sens[6] = sens[7] = 0.0;
  sens[8] = 0.0;

  for (int k = 0; k < NUM_FAIL_ANGLES; k++) {
    TacsScalar angle = 0.0, factor = 1.0;
    if (k == 0) {
      angle = 0.0;
      factor = f0 / (f0 + epsilon);
    } else if (k == 1) {
      angle = 0.25 * M_PI;
      factor = f45 / (f45 + epsilon);
    } else if (k == 2) {
      angle = -0.25 * M_PI;
      factor = f45 / (f45 + epsilon);
    } else {
      angle = 0.5 * M_PI;
      factor = f90 / (f90 + epsilon);
    }

    // Compute the strain and failure criteria at the upper surface
    TacsScalar e[3], e_sens[3];
    e[0] = strain[0] + 0.5 * t * strain[3];
    e[1] = strain[1] + 0.5 * t * strain[4];
    e[2] = strain[2] + 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    sens[0] += factor * weights[2 * k] * e_sens[0];
    sens[1] += factor * weights[2 * k] * e_sens[1];
    sens[2] += factor * weights[2 * k] * e_sens[2];

    sens[3] += 0.5 * factor * t * weights[2 * k] * e_sens[0];
    sens[4] += 0.5 * factor * t * weights[2 * k] * e_sens[1];
    sens[5] += 0.5 * factor * t * weights[2 * k] * e_sens[2];

    // Compute the strain at the lower sufrace
    e[0] = strain[0] - 0.5 * t * strain[3];
    e[1] = strain[1] - 0.5 * t * strain[4];
    e[2] = strain[2] - 0.5 * t * strain[5];

    orthoPly->failureStrainSens(angle, e, e_sens);

    sens[0] += factor * weights[2 * k + 1] * e_sens[0];
    sens[1] += factor * weights[2 * k + 1] * e_sens[1];
    sens[2] += factor * weights[2 * k + 1] * e_sens[2];

    sens[3] -= 0.5 * factor * t * weights[2 * k + 1] * e_sens[0];
    sens[4] -= 0.5 * factor * t * weights[2 * k + 1] * e_sens[1];
    sens[5] -= 0.5 * factor * t * weights[2 * k + 1] * e_sens[2];
  }
}

/**
  Compute the failure load for a series of ply angles and take the
  approximate maximum using the KS function.
*/
TacsScalar TACSLamParamShellConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[]) {
  TacsScalar fvals[2 * NUM_FAIL_ANGLES];
  TacsScalar max, ks_sum = 0.0;

  computeFailure(strain, fvals, &max);
  for (int k = 0; k < 2 * NUM_FAIL_ANGLES; k++) {
    ks_sum += exp(ksWeight * (fvals[k] - max));
  }

  return max + log(ks_sum) / ksWeight;
}

/**
  Compute the derivative of the failure load w.r.t. the strain
  values
*/
TacsScalar TACSLamParamShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], TacsScalar sens[]) {
  TacsScalar fvals[2 * NUM_FAIL_ANGLES], weights[2 * NUM_FAIL_ANGLES];
  TacsScalar max, ks_sum = 0.0;

  computeFailure(strain, fvals, &max);

  for (int k = 0; k < 2 * NUM_FAIL_ANGLES; k++) {
    weights[k] = exp(ksWeight * (fvals[k] - max));
    ks_sum += weights[k];
  }

  TacsScalar inv_ks_sum = 1.0 / ks_sum;
  for (int k = 0; k < 2 * NUM_FAIL_ANGLES; k++) {
    weights[k] *= inv_ks_sum;
  }

  computeFailureStrainSens(strain, weights, sens);

  return max + log(ks_sum) / ksWeight;
}

/*!
  Functions to determine the derivative of the failure
  load w.r.t. the design variables
*/
void TACSLamParamShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar fvals[2 * NUM_FAIL_ANGLES], weights[2 * NUM_FAIL_ANGLES];
  TacsScalar max, ks_sum = 0.0;

  computeFailure(strain, fvals, &max);

  for (int k = 0; k < 2 * NUM_FAIL_ANGLES; k++) {
    weights[k] = exp(ksWeight * (fvals[k] - max));
    ks_sum += weights[k];
  }

  TacsScalar inv_ks_sum = 1.0 / ks_sum;
  for (int k = 0; k < 2 * NUM_FAIL_ANGLES; k++) {
    weights[k] *= inv_ks_sum;
  }

  int index = 0;
  if (tNum >= 0) {
    for (int k = 0; k < NUM_FAIL_ANGLES; k++) {
      TacsScalar angle = 0.0, factor = scale;
      if (k == 0) {
        angle = 0.0;
        factor *= f0 / (f0 + epsilon);
      } else if (k == 1) {
        angle = 0.25 * M_PI;
        factor *= f45 / (f45 + epsilon);
      } else if (k == 2) {
        angle = -0.25 * M_PI;
        factor *= f45 / (f45 + epsilon);
      } else {
        angle = 0.5 * M_PI;
        factor *= f90 / (f90 + epsilon);
      }

      // Compute the strain and failure criteria at the upper surface
      TacsScalar e[3], e_sens[3];
      e[0] = strain[0] + 0.5 * t * strain[3];
      e[1] = strain[1] + 0.5 * t * strain[4];
      e[2] = strain[2] + 0.5 * t * strain[5];
      orthoPly->failureStrainSens(angle, e, e_sens);

      dfdx[index] += 0.5 * factor * weights[2 * k] *
                     (strain[3] * e_sens[0] + strain[4] * e_sens[1] +
                      strain[5] * e_sens[2]);

      // Compute the strain at the lower sufrace
      e[0] = strain[0] - 0.5 * t * strain[3];
      e[1] = strain[1] - 0.5 * t * strain[4];
      e[2] = strain[2] - 0.5 * t * strain[5];
      orthoPly->failureStrainSens(angle, e, e_sens);

      dfdx[index] -= 0.5 * factor * weights[2 * k + 1] *
                     (strain[3] * e_sens[0] + strain[4] * e_sens[1] +
                      strain[5] * e_sens[2]);
    }

    index++;
  }
  if (n0 >= 0) {  // If one of the ply fractions is a variable, they all are
    for (int k = 0; k < NUM_FAIL_ANGLES; k++) {
      TacsScalar angle = 0.0, dfactor = 0.0;
      if (k == 0) {
        angle = 0.0;
        TacsScalar d = 1.0 / (f0 + epsilon);
        dfactor = scale * epsilon * d * d;
      } else if (k == 1) {
        angle = 0.25 * M_PI;
        TacsScalar d = 1.0 / (f45 + epsilon);
        dfactor = scale * epsilon * d * d;
      } else if (k == 2) {
        angle = -0.25 * M_PI;
        TacsScalar d = 1.0 / (f45 + epsilon);
        dfactor = scale * epsilon * d * d;
      } else {
        angle = 0.5 * M_PI;
        TacsScalar d = 1.0 / (f90 + epsilon);
        dfactor = scale * epsilon * d * d;
      }

      // Compute the strain and failure criteria at the upper surface
      TacsScalar e[3];
      e[0] = strain[0] + 0.5 * t * strain[3];
      e[1] = strain[1] + 0.5 * t * strain[4];
      e[2] = strain[2] + 0.5 * t * strain[5];
      dfdx[index] += dfactor * weights[2 * k] * orthoPly->failure(angle, e);

      // Compute the failure criteria at the lower surface
      e[0] = strain[0] - 0.5 * t * strain[3];
      e[1] = strain[1] - 0.5 * t * strain[4];
      e[2] = strain[2] - 0.5 * t * strain[5];
      dfdx[index] += dfactor * weights[2 * k + 1] * orthoPly->failure(angle, e);

      index++;
    }
  }
}

// Retrieve the design variable for plotting purposes
TacsScalar TACSLamParamShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return t;
  } else if (index == 1) {
    return f0;
  } else if (index == 2) {
    return f45;
  } else if (index == 3) {
    return f90;
  } else if (index == 4) {
    return W1;
  } else if (index == 5) {
    return W3;
  }

  return 0.0;
}
