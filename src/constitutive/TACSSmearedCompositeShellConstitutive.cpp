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

#include "TACSSmearedCompositeShellConstitutive.h"

#include "TACSElementAlgebra.h"
#include "TacsUtilities.h"

const char *TACSSmearedCompositeShellConstitutive::constName =
    "TACSSmearedCompositeShellConstitutive";

/*
  Create the shell constitutive
*/
TACSSmearedCompositeShellConstitutive::TACSSmearedCompositeShellConstitutive(
    int _num_plies, TACSOrthotropicPly **_ply_props, TacsScalar _thickness,
    const TacsScalar *_ply_angles, const TacsScalar *_ply_fractions,
    int _thickness_dv_num, const int *_ply_fraction_dv_nums,
    TacsScalar _thickness_lb, TacsScalar _thickness_ub,
    const TacsScalar *_ply_fraction_lb, const TacsScalar *_ply_fraction_ub) {
  num_plies = _num_plies;

  thickness = _thickness;
  thickness_dv_num = _thickness_dv_num;
  thickness_lb = _thickness_lb;
  thickness_ub = _thickness_ub;

  kcorr = 5.0 / 6.0;

  ply_fractions = new TacsScalar[num_plies];
  ply_fraction_dv_nums = new int[num_plies];
  ply_angles = new TacsScalar[num_plies];
  ply_fraction_lb = new TacsScalar[num_plies];
  ply_fraction_ub = new TacsScalar[num_plies];
  ply_props = new TACSOrthotropicPly *[num_plies];

  for (int i = 0; i < num_plies; i++) {
    ply_props[i] = _ply_props[i];
    ply_props[i]->incref();

    ply_fractions[i] = _ply_fractions[i];
    ply_angles[i] = _ply_angles[i];

    if (_ply_fraction_dv_nums) {
      ply_fraction_dv_nums[i] = _ply_fraction_dv_nums[i];
    } else {
      ply_fraction_dv_nums[i] = -1;
    }

    if (_ply_fraction_lb) {
      ply_fraction_lb[i] = _ply_fraction_lb[i];
    } else {
      ply_fraction_lb[i] = 0.0;
    }

    if (_ply_fraction_ub) {
      ply_fraction_ub[i] = _ply_fraction_ub[i];
    } else {
      ply_fraction_ub[i] = 1.0;
    }
  }

  ks_weight = 100.0;
  nfvals = 2 * num_plies;
  fvals = new TacsScalar[nfvals];
  dks_vals = new TacsScalar[nfvals];
  dfvals = new TacsScalar *[nfvals];
  for (int i = 0; i < nfvals; i++) {
    dfvals[i] = new TacsScalar[NUM_STRESSES];
  }
}

TACSSmearedCompositeShellConstitutive::
    ~TACSSmearedCompositeShellConstitutive() {
  for (int i = 0; i < num_plies; i++) {
    ply_props[i]->decref();
  }
  delete[] ply_props;
  delete[] ply_angles;
  delete[] ply_fractions;
  delete[] ply_fraction_dv_nums;
  delete[] ply_fraction_lb;
  delete[] ply_fraction_ub;
  delete[] fvals;
  delete[] dks_vals;
  for (int i = 0; i < nfvals; i++) {
    delete[] dfvals[i];
  }
  delete[] dfvals;
}

int TACSSmearedCompositeShellConstitutive::getDesignVarNums(int elemIndex,
                                                            int dvLen,
                                                            int dvNums[]) {
  int index = 0;
  if (thickness_dv_num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = thickness_dv_num;
    }
    index++;
  }
  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      if (dvNums && dvLen > index) {
        dvNums[index] = ply_fraction_dv_nums[i];
      }
      index++;
    }
  }
  return index;
}

int TACSSmearedCompositeShellConstitutive::setDesignVars(
    int elemIndex, int dvLen, const TacsScalar dvs[]) {
  int index = 0;
  if (thickness_dv_num >= 0) {
    thickness = dvs[index];
    index++;
  }
  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      ply_fractions[i] = dvs[index];
      index++;
    }
  }
  return index;
}

int TACSSmearedCompositeShellConstitutive::getDesignVars(int elemIndex,
                                                         int dvLen,
                                                         TacsScalar dvs[]) {
  int index = 0;
  if (thickness_dv_num >= 0) {
    dvs[index] = thickness;
    index++;
  }
  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      dvs[index] = ply_fractions[i];
      index++;
    }
  }
  return index;
}

int TACSSmearedCompositeShellConstitutive::getDesignVarRange(int elemIndex,
                                                             int dvLen,
                                                             TacsScalar lb[],
                                                             TacsScalar ub[]) {
  int index = 0;
  if (thickness_dv_num >= 0) {
    lb[index] = thickness_lb;
    ub[index] = thickness_ub;
    index++;
  }
  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      lb[index] = ply_fraction_lb[i];
      ub[index] = ply_fraction_ub[i];
      index++;
    }
  }
  return index;
}

// Evaluate the material density
TacsScalar TACSSmearedCompositeShellConstitutive::evalDensity(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  // Compute the thickness-weighted density across all plies
  TacsScalar rhot = 0.0;
  for (int i = 0; i < num_plies; i++) {
    TacsScalar rho_ply = ply_props[i]->getDensity();
    TacsScalar t_ply = thickness * ply_fractions[i];
    rhot += rho_ply * t_ply;
  }
  return rhot;
}

void TACSSmearedCompositeShellConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  int index = 0;

  if (thickness_dv_num >= 0) {
    for (int i = 0; i < num_plies; i++) {
      TacsScalar rho_ply = ply_props[i]->getDensity();
      dfdx[index] += scale * rho_ply * ply_fractions[i];
    }
    index++;
  }

  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      TacsScalar rho_ply = ply_props[i]->getDensity();
      dfdx[index] += scale * rho_ply * thickness;
      index++;
    }
  }
}

TacsScalar TACSSmearedCompositeShellConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return thickness;
  } else {
    for (int i = 0; i < num_plies; i++) {
      if (index == i + 1) {
        return ply_fractions[i];
      }
    }
  }
  return 0.0;
}

// Evaluate the mass moments
void TACSSmearedCompositeShellConstitutive::evalMassMoments(
    int elemIndex, const double pt[], const TacsScalar X[],
    TacsScalar moments[]) {
  moments[0] = 0.0;
  moments[1] = 0.0;
  moments[2] = 0.0;

  // Compute the contribution to the mass moment from each layer
  for (int i = 0; i < num_plies; i++) {
    TacsScalar rho_ply = ply_props[i]->getDensity();

    moments[0] += thickness * rho_ply * ply_fractions[i];
    moments[2] +=
        thickness * thickness * thickness * rho_ply * ply_fractions[i] / 12.0;
  }
}

// Add the sensitivity of the mass moments
void TACSSmearedCompositeShellConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  if (thickness_dv_num >= 0) {
    for (int i = 0; i < num_plies; i++) {
      TacsScalar rho_ply = ply_props[i]->getDensity();

      dfdx[index] +=
          rho_ply * (scale[0] + 0.25 * thickness * thickness * scale[2]);
    }
    index++;
  }
  for (int i = 0; i < num_plies; i++) {
    if (ply_fraction_dv_nums[i] >= 0) {
      TacsScalar rho_ply = ply_props[i]->getDensity();

      dfdx[index] +=
          rho_ply * (scale[0] * thickness +
                     thickness * thickness * thickness * scale[2] / 12.0);
      index++;
    }
  }
}

// Evaluate the specific heat
TacsScalar TACSSmearedCompositeShellConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

// Evaluate the FSDT stiffness matrices
TacsScalar TACSSmearedCompositeShellConstitutive::evalFSDTStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar A[],
    TacsScalar B[], TacsScalar D[], TacsScalar As[]) {
  TacsScalar drill;

  // Zero the stiffness matrices
  for (int k = 0; k < 6; k++) {
    A[k] = B[k] = D[k] = 0.0;
  }

  for (int k = 0; k < 3; k++) {
    As[k] = 0.0;
  }

  // Compute the contribution to the stiffness from each layer
  for (int k = 0; k < num_plies; k++) {
    TacsScalar Qbar[6], Abar[3];
    ply_props[k]->calculateQbar(ply_angles[k], Qbar);
    ply_props[k]->calculateAbar(ply_angles[k], Abar);

    TacsScalar a = ply_fractions[k] * thickness;
    TacsScalar d =
        1.0 / 12.0 * ply_fractions[k] * (thickness * thickness * thickness);

    for (int i = 0; i < 6; i++) {
      A[i] += a * Qbar[i];
      D[i] += d * Qbar[i];
    }

    for (int i = 0; i < 3; i++) {
      As[i] += kcorr * a * Abar[i];
    }
  }

  drill = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);

  return drill;
}

// Evaluate the stress
void TACSSmearedCompositeShellConstitutive::evalStress(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       const TacsScalar e[],
                                                       TacsScalar s[]) {
  TacsScalar A[6], B[6], D[6], As[3], drill;
  drill = evalFSDTStiffness(elemIndex, pt, X, A, B, D, As);

  // Evaluate the stress
  TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);
}

// Add the contribution
void TACSSmearedCompositeShellConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar dA[6] = {0.0}, dD[6] = {0.0}, dAs[3] = {0.0}, ddrill;

  int index = 0;
  if (thickness_dv_num >= 0) {
    for (int k = 0; k < num_plies; k++) {
      TacsScalar Qbar[6], Abar[3];
      ply_props[k]->calculateQbar(ply_angles[k], Qbar);
      ply_props[k]->calculateAbar(ply_angles[k], Abar);

      TacsScalar da = ply_fractions[k];
      TacsScalar dd = 0.25 * ply_fractions[k] * (thickness * thickness);

      for (int i = 0; i < 6; i++) {
        dA[i] += da * Qbar[i];
        dD[i] += dd * Qbar[i];
      }

      for (int i = 0; i < 3; i++) {
        dAs[i] += kcorr * da * Abar[i];
      }
    }

    ddrill = 0.5 * DRILLING_REGULARIZATION * (dAs[0] + dAs[2]);

    dfdx[index] +=
        scale * (mat3x3SymmInner(dA, &psi[0], &e[0]) +
                 mat3x3SymmInner(dD, &psi[3], &e[3]) +
                 mat2x2SymmInner(dAs, &psi[6], &e[6]) + ddrill * psi[8] * e[8]);
    index++;
  }

  for (int k = 0; k < num_plies; k++) {
    if (ply_fraction_dv_nums[k] >= 0) {
      TacsScalar Qbar[6], Abar[3];
      ply_props[k]->calculateQbar(ply_angles[k], Qbar);
      ply_props[k]->calculateAbar(ply_angles[k], Abar);

      TacsScalar da = thickness;
      TacsScalar dd = 1.0 / 12.0 * (thickness * thickness * thickness);

      for (int i = 0; i < 6; i++) {
        dA[i] = da * Qbar[i];
        dD[i] = dd * Qbar[i];
      }

      for (int i = 0; i < 3; i++) {
        dAs[i] = kcorr * da * Abar[i];
      }

      ddrill = 0.5 * DRILLING_REGULARIZATION * (dAs[0] + dAs[2]);

      dfdx[index] += scale * (mat3x3SymmInner(dA, &psi[0], &e[0]) +
                              mat3x3SymmInner(dD, &psi[3], &e[3]) +
                              mat2x2SymmInner(dAs, &psi[6], &e[6]) +
                              ddrill * psi[8] * e[8]);
      index++;
    }
  }
}

/*
  Compute the most critical failure criteria for the laminate
*/
TacsScalar TACSSmearedCompositeShellConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[]) {
  evalPlyTopBottomFailure(strain, fvals);
  TacsScalar ks_val = ksAggregation(fvals, nfvals, ks_weight);

  return ks_val;
}

/*
  Compute the most critical failure criteria for the laminate
*/
void TACSSmearedCompositeShellConstitutive::evalPlyTopBottomFailure(
    const TacsScalar strain[], TacsScalar fvals[]) {
  // Compute the total thickness of the laminate
  TacsScalar tb = -0.5 * thickness;
  TacsScalar tt = 0.5 * thickness;

  for (int i = 0; i < num_plies; i++) {
    TacsScalar lamStrain[3] = {0.0};

    getLaminaStrain(strain, tb, lamStrain);
    fvals[2 * i] = ply_props[i]->failure(ply_angles[i], lamStrain);

    getLaminaStrain(strain, tt, lamStrain);
    fvals[2 * i + 1] = ply_props[i]->failure(ply_angles[i], lamStrain);
  }
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSSmearedCompositeShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = 0.0;
  sens[3] = sens[4] = sens[5] = 0.0;
  sens[6] = sens[7] = sens[8] = 0.0;

  // Compute the total thickness of the laminate
  TacsScalar tb = -0.5 * thickness;
  TacsScalar tt = 0.5 * thickness;

  evalPlyTopBottomFailure(strain, fvals);

  for (int i = 0; i < num_plies; i++) {
    TacsScalar lamStrain[3], phi[3];

    // plate bottom stress
    getLaminaStrain(strain, tb, lamStrain);
    TacsScalar fval_bot =
        ply_props[i]->failureStrainSens(ply_angles[i], lamStrain, phi);
    dfvals[2 * i + 0][0] = phi[0];
    dfvals[2 * i + 0][1] = phi[1];
    dfvals[2 * i + 0][2] = phi[2];
    dfvals[2 * i + 0][3] = tb * phi[0];
    dfvals[2 * i + 0][4] = tb * phi[1];
    dfvals[2 * i + 0][5] = tb * phi[2];
    dfvals[2 * i + 0][6] = dfvals[2 * i + 0][7] = dfvals[2 * i + 0][8] = 0.0;

    // plate top stress
    getLaminaStrain(strain, tt, lamStrain);
    TacsScalar fval_top =
        ply_props[i]->failureStrainSens(ply_angles[i], lamStrain, phi);
    dfvals[2 * i + 1][0] = phi[0];
    dfvals[2 * i + 1][1] = phi[1];
    dfvals[2 * i + 1][2] = phi[2];
    dfvals[2 * i + 1][3] = tt * phi[0];
    dfvals[2 * i + 1][4] = tt * phi[1];
    dfvals[2 * i + 1][5] = tt * phi[2];
    dfvals[2 * i + 1][6] = dfvals[2 * i + 1][7] = dfvals[2 * i + 1][8] = 0.0;
  }

  TacsScalar ks_val = ksAggregationSensProduct(fvals, nfvals, NUM_STRESSES,
                                               ks_weight, dfvals, sens);

  return ks_val;
}

// Add the derivative of the failure criteria w.r.t. the design variables
void TACSSmearedCompositeShellConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  if (thickness_dv_num >= 0 && dvLen >= 1) {
    // Compute the total thickness of the laminate
    TacsScalar tb = -0.5 * thickness;
    TacsScalar tt = 0.5 * thickness;

    evalPlyTopBottomFailure(e, fvals);
    TacsScalar ks_val = ksAggregationSens(fvals, nfvals, ks_weight, dks_vals);

    for (int i = 0; i < num_plies; i++) {
      TacsScalar lamStrain[3], phi[3];

      // plate bottom stress
      getLaminaStrain(e, tb, lamStrain);
      TacsScalar fval_bot =
          ply_props[i]->failureStrainSens(ply_angles[i], lamStrain, phi);
      dfdx[index] += -0.5 * scale * dks_vals[2 * i + 0] *
                     (phi[0] * e[3] + phi[1] * e[4] + phi[2] * e[5]);

      // plate top stress
      getLaminaStrain(e, tt, lamStrain);
      TacsScalar fval_top =
          ply_props[i]->failureStrainSens(ply_angles[i], lamStrain, phi);
      dfdx[index] += 0.5 * scale * dks_vals[2 * i + 1] *
                     (phi[0] * e[3] + phi[1] * e[4] + phi[2] * e[5]);
    }

    index++;
  }
}

/*
  Get the strain in a single ply
*/
void TACSSmearedCompositeShellConstitutive::getLaminaStrain(
    const TacsScalar rmStrain[], TacsScalar tp, TacsScalar strain[]) {
  strain[0] = rmStrain[0] + tp * rmStrain[3];
  strain[1] = rmStrain[1] + tp * rmStrain[4];
  strain[2] = rmStrain[2] + tp * rmStrain[5];
}

// Evaluate the tangent stiffness
void TACSSmearedCompositeShellConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];

  C[21] = evalFSDTStiffness(elemIndex, pt, X, A, B, D, As);
}

// Evaluate the thermal strain
void TACSSmearedCompositeShellConstitutive::evalThermalStrain(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar theta,
    TacsScalar e[]) {
  e[0] = e[1] = e[2] = 0.0;
  e[3] = e[4] = e[5] = 0.0;
  e[6] = e[7] = e[8] = 0.0;
}

// Evaluate the heat flux, given the thermal gradient
void TACSSmearedCompositeShellConstitutive::evalHeatFlux(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar grad[], TacsScalar flux[]) {
  flux[0] = 0.0;
  flux[1] = 0.0;
}

// Evaluate the tangent of the heat flux
void TACSSmearedCompositeShellConstitutive::evalTangentHeatFlux(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar Kc[]) {
  Kc[0] = Kc[1] = Kc[2] = 0.0;
}

/*
  Return the constitutive name
*/
const char *TACSSmearedCompositeShellConstitutive::getObjectName() {
  return constName;
}

/*
  Get laminate thicknesses
*/
TacsScalar TACSSmearedCompositeShellConstitutive::getLaminateThickness() {
  return thickness;
}

/*
  Get ply angles
*/
void TACSSmearedCompositeShellConstitutive::getPlyAngles(
    TacsScalar *_ply_angles) {
  for (int i = 0; i < num_plies; i++) {
    _ply_angles[i] = ply_angles[i];
  }
}

/*
  Get ply angles
*/
void TACSSmearedCompositeShellConstitutive::getPlyFractions(
    TacsScalar *_ply_fractions) {
  for (int i = 0; i < num_plies; i++) {
    _ply_fractions[i] = ply_fractions[i];
  }
}
