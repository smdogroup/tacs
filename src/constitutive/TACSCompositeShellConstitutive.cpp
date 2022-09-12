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

#include "TACSCompositeShellConstitutive.h"

const char *TACSCompositeShellConstitutive::constName =
    "TACSCompositeShellConstitutive";

/*
  Create the shell constitutive
*/
TACSCompositeShellConstitutive::TACSCompositeShellConstitutive(
    int _num_plies, TACSOrthotropicPly **_ply_props,
    const TacsScalar *_ply_thickness, const TacsScalar *_ply_angles,
    TacsScalar _kcorr) {
  num_plies = _num_plies;
  ply_thickness = new TacsScalar[num_plies];
  ply_angles = new TacsScalar[num_plies];
  ply_props = new TACSOrthotropicPly *[num_plies];

  for (int i = 0; i < num_plies; i++) {
    ply_props[i] = _ply_props[i];
    ply_props[i]->incref();

    ply_thickness[i] = _ply_thickness[i];
    ply_angles[i] = _ply_angles[i];
  }

  kcorr = _kcorr;
}

TACSCompositeShellConstitutive::~TACSCompositeShellConstitutive() {
  for (int i = 0; i < num_plies; i++) {
    ply_props[i]->decref();
  }
  delete[] ply_props;
  delete[] ply_thickness;
  delete[] ply_angles;
}

// Evaluate the material density
TacsScalar TACSCompositeShellConstitutive::evalDensity(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[]) {
  // Compute the thickness-weighted density across all plies
  TacsScalar rhot = 0.0;
  for (int i = 0; i < num_plies; i++) {
    TacsScalar rho_ply = ply_props[i]->getDensity();
    TacsScalar t_ply = ply_thickness[i];
    rhot += rho_ply * t_ply;
  }
  return rhot;
}

// Evaluate the mass moments
void TACSCompositeShellConstitutive::evalMassMoments(int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     TacsScalar moments[]) {
  moments[0] = 0.0;
  moments[1] = 0.0;
  moments[2] = 0.0;

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += ply_thickness[i];
  }

  // Compute the contribution to the mass moment from each layer
  TacsScalar t0 = -0.5 * t;
  for (int i = 0; i < num_plies; i++) {
    TacsScalar rho_ply = ply_props[i]->getDensity();
    TacsScalar t1 = t0 + ply_thickness[i];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5 * (t1 * t1 - t0 * t0);
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    moments[0] += a * rho_ply;
    moments[1] += b * rho_ply;
    moments[2] += d * rho_ply;

    // Update the position of the bottom interface
    t0 = t1;
  }
}

// Evaluate the specific heat
TacsScalar TACSCompositeShellConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

// Evaluate the stresss
void TACSCompositeShellConstitutive::evalStress(int elemIndex,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar e[],
                                                TacsScalar s[]) {
  TacsScalar A[6], B[6], D[6], As[3], drill;

  // Zero the stiffness matrices
  for (int k = 0; k < 6; k++) {
    A[k] = B[k] = D[k] = 0.0;
  }

  for (int k = 0; k < 3; k++) {
    As[k] = 0.0;
  }

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += ply_thickness[i];
  }

  // Compute the contribution to the stiffness from each layer
  TacsScalar t0 = -0.5 * t;
  for (int k = 0; k < num_plies; k++) {
    TacsScalar Qbar[6], Abar[3];
    ply_props[k]->calculateQbar(ply_angles[k], Qbar);
    ply_props[k]->calculateAbar(ply_angles[k], Abar);

    TacsScalar t1 = t0 + ply_thickness[k];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5 * (t1 * t1 - t0 * t0);
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    for (int i = 0; i < 6; i++) {
      A[i] += a * Qbar[i];
      B[i] += b * Qbar[i];
      D[i] += d * Qbar[i];
    }

    for (int i = 0; i < 3; i++) {
      As[i] += kcorr * a * Abar[i];
    }

    // Update the position of the bottom interface
    t0 = t1;
  }

  drill = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);

  // Evaluate the stress
  TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);
}

/*
  Compute the most critical failure criteria for the laminate
*/
TacsScalar TACSCompositeShellConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[]) {
  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += ply_thickness[i];
  }
  TacsScalar t0 = -0.5 * t;

  // Keep track of the maximum failure criterion
  TacsScalar max = 0.0;

  for (int i = 0; i < num_plies; i++) {
    TacsScalar lamStrain[3];
    TacsScalar tp = t0 + 0.5 * ply_thickness[i];
    getLaminaStrain(lamStrain, strain, tp);
    TacsScalar fval = ply_props[i]->failure(ply_angles[i], lamStrain);

    if (TacsRealPart(fval) > TacsRealPart(max)) {
      max = fval;
    }
    t0 += ply_thickness[i];
  }

  return max;
}

// Evaluate the derivative of the failure criteria w.r.t. the strain
TacsScalar TACSCompositeShellConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar strain[], TacsScalar sens[]) {
  sens[0] = sens[1] = sens[2] = 0.0;
  sens[3] = sens[4] = sens[5] = 0.0;
  sens[6] = sens[7] = sens[8] = 0.0;

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += ply_thickness[i];
  }
  TacsScalar t0 = -0.5 * t;

  // Keep track of the maximum failure criterion
  TacsScalar max = 0.0;

  for (int i = 0; i < num_plies; i++) {
    TacsScalar lamStrain[3], phi[3];
    TacsScalar tp = t0 + 0.5 * ply_thickness[i];
    getLaminaStrain(lamStrain, strain, tp);
    TacsScalar fval = ply_props[i]->failure(ply_angles[i], lamStrain);
    ply_props[i]->failureStrainSens(ply_angles[i], lamStrain, phi);

    if (TacsRealPart(fval) > TacsRealPart(max)) {
      max = fval;
      sens[0] = phi[0];
      sens[1] = phi[1];
      sens[2] = phi[2];
      sens[3] = tp * phi[0];
      sens[4] = tp * phi[1];
      sens[5] = tp * phi[2];
    }
    t0 += ply_thickness[i];
  }

  return max;
}

/*
  Get the strain in a single ply
*/
void TACSCompositeShellConstitutive::getLaminaStrain(
    TacsScalar strain[], const TacsScalar rmStrain[], TacsScalar tp) {
  strain[0] = rmStrain[0] + tp * rmStrain[3];
  strain[1] = rmStrain[1] + tp * rmStrain[4];
  strain[2] = rmStrain[2] + tp * rmStrain[5];
}

// Evaluate the tangent stiffness
void TACSCompositeShellConstitutive::evalTangentStiffness(int elemIndex,
                                                          const double pt[],
                                                          const TacsScalar X[],
                                                          TacsScalar C[]) {
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];

  // Zero the stiffness matrices
  for (int k = 0; k < 6; k++) {
    A[k] = B[k] = D[k] = 0.0;
  }

  for (int k = 0; k < 3; k++) {
    As[k] = 0.0;
  }

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += ply_thickness[i];
  }

  // Compute the contribution to the stiffness from each layer
  TacsScalar t0 = -0.5 * t;
  for (int k = 0; k < num_plies; k++) {
    TacsScalar Qbar[6], Abar[3];
    ply_props[k]->calculateQbar(ply_angles[k], Qbar);
    ply_props[k]->calculateAbar(ply_angles[k], Abar);

    TacsScalar t1 = t0 + ply_thickness[k];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5 * (t1 * t1 - t0 * t0);
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    for (int i = 0; i < 6; i++) {
      A[i] += a * Qbar[i];
      B[i] += b * Qbar[i];
      D[i] += d * Qbar[i];
    }

    for (int i = 0; i < 3; i++) {
      As[i] += kcorr * a * Abar[i];
    }

    // Update the position of the bottom interface
    t0 = t1;
  }

  C[21] = 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
}

// Evaluate the thermal strain
void TACSCompositeShellConstitutive::evalThermalStrain(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       TacsScalar theta,
                                                       TacsScalar e[]) {
  e[0] = e[1] = e[2] = 0.0;
  e[3] = e[4] = e[5] = 0.0;
  e[6] = e[7] = e[8] = 0.0;
}

// Evaluate the heat flux, given the thermal gradient
void TACSCompositeShellConstitutive::evalHeatFlux(int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar grad[],
                                                  TacsScalar flux[]) {
  flux[0] = 0.0;
  flux[1] = 0.0;
}

// Evaluate the tangent of the heat flux
void TACSCompositeShellConstitutive::evalTangentHeatFlux(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[],
                                                         TacsScalar Kc[]) {
  Kc[0] = Kc[1] = Kc[2] = 0.0;
}

/*
  Return the constitutive name
*/
const char *TACSCompositeShellConstitutive::getObjectName() {
  return constName;
}
