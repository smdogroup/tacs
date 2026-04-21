#include "TACSCompositeTubeBeamConstitutive.h"

#include <math.h>
#include <string.h>

TACSCompositeTubeBeamConstitutive::TACSCompositeTubeBeamConstitutive(
    TacsScalar E11_in, TacsScalar E22, TacsScalar G12, TacsScalar nu12,
    TacsScalar rho_in, TacsScalar X_c_in, TacsScalar X_t_in,
    const TacsScalar *layup_angles, int n_plies,
    TacsScalar inner_init, TacsScalar wall_init,
    int inner_dv, int wall_dv,
    TacsScalar inner_lb, TacsScalar inner_ub,
    TacsScalar wall_lb, TacsScalar wall_ub,
    TacsScalar buckle_length,
    TacsScalar buckle_length_factor) {
  // Store ply-level material properties
  E11 = E11_in;
  rho = rho_in;
  X_c = X_c_in;
  X_t = X_t_in;

  // Compute CLT-smeared effective moduli (mean Qbar over equal-thickness plies)
  TacsScalar nu21 = nu12 * E22 / E11_in;
  TacsScalar denom = 1.0 - nu12 * nu21;
  TacsScalar Q11 = E11_in / denom;
  TacsScalar Q22 = E22 / denom;
  TacsScalar Q12 = nu12 * E22 / denom;
  TacsScalar Q66 = G12;

  TacsScalar sum_Qbar11 = 0.0;
  TacsScalar sum_Qbar12 = 0.0;
  TacsScalar sum_Qbar66 = 0.0;
  for (int k = 0; k < n_plies; k++) {
    TacsScalar theta = layup_angles[k];
    TacsScalar m = cos(theta);
    TacsScalar n = sin(theta);
    TacsScalar m2 = m * m;
    TacsScalar n2 = n * n;
    TacsScalar m4 = m2 * m2;
    TacsScalar n4 = n2 * n2;
    TacsScalar m2n2 = m2 * n2;

    TacsScalar Qbar11 = m4 * Q11 + 2.0 * m2n2 * (Q12 + 2.0 * Q66) + n4 * Q22;
    TacsScalar Qbar12 = m2n2 * (Q11 + Q22 - 4.0 * Q66) + (m4 + n4) * Q12;
    TacsScalar Qbar66 =
        (m2 - n2) * (m2 - n2) * Q66 + m2n2 * (Q11 + Q22 - 2.0 * Q12);

    sum_Qbar11 += Qbar11;
    sum_Qbar12 += Qbar12;
    sum_Qbar66 += Qbar66;
  }
  E_eff = sum_Qbar11 / n_plies;
  G_eff = sum_Qbar66 / n_plies;
  TacsScalar A12_mean = sum_Qbar12 / n_plies;
  nu_eff = A12_mean / E_eff;

  // Geometry DVs
  inner = inner_init;
  wall = wall_init;
  innerDV = inner_dv;
  wallDV = wall_dv;
  innerLb = inner_lb;
  innerUb = inner_ub;
  wallLb = wall_lb;
  wallUb = wall_ub;

  // Buckling parameters
  buckleLength = buckle_length;
  buckleLengthFactor = buckle_length_factor;

  ks_weight = 100.0;
}

TACSCompositeTubeBeamConstitutive::~TACSCompositeTubeBeamConstitutive() {}

int TACSCompositeTubeBeamConstitutive::getDesignVarNums(int elemIndex,
                                                        int dvLen,
                                                        int dvNums[]) {
  int index = 0;
  if (innerDV >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = innerDV;
    }
    index++;
  }
  if (wallDV >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = wallDV;
    }
    index++;
  }
  return index;
}

int TACSCompositeTubeBeamConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                     const TacsScalar dvs[]) {
  int index = 0;
  if (innerDV >= 0) {
    inner = dvs[index];
    index++;
  }
  if (wallDV >= 0) {
    wall = dvs[index];
    index++;
  }
  return index;
}

int TACSCompositeTubeBeamConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                     TacsScalar dvs[]) {
  int index = 0;
  if (innerDV >= 0) {
    dvs[index] = inner;
    index++;
  }
  if (wallDV >= 0) {
    dvs[index] = wall;
    index++;
  }
  return index;
}

int TACSCompositeTubeBeamConstitutive::getDesignVarRange(int elemIndex,
                                                         int dvLen,
                                                         TacsScalar lb[],
                                                         TacsScalar ub[]) {
  int index = 0;
  if (innerDV >= 0) {
    lb[index] = innerLb;
    ub[index] = innerUb;
    index++;
  }
  if (wallDV >= 0) {
    lb[index] = wallLb;
    ub[index] = wallUb;
    index++;
  }
  return index;
}

TacsScalar TACSCompositeTubeBeamConstitutive::evalDensity(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  return rho * A;
}

void TACSCompositeTubeBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  TacsScalar d0 = inner + wall;
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    dfdx[index] += scale * rho * dA;
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    dfdx[index] += scale * rho * dA;
    index++;
  }
}

void TACSCompositeTubeBeamConstitutive::evalMassMoments(
    int elemIndex, const double pt[], const TacsScalar X[],
    TacsScalar moments[]) {
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  moments[0] = rho * A;
  moments[1] = 0.0;
  moments[2] = 0.0;
  moments[3] = rho * Ia;
  moments[4] = rho * Ia;
  moments[5] = 0.0;
}

void TACSCompositeTubeBeamConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;
    dfdx[index] += rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0) / 16.0;
    dfdx[index] += rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
}

TacsScalar TACSCompositeTubeBeamConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  return 0.0;
}

void TACSCompositeTubeBeamConstitutive::evalStress(int elemIndex,
                                                   const double pt[],
                                                   const TacsScalar X[],
                                                   const TacsScalar e[],
                                                   TacsScalar s[]) {
  TacsScalar kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  s[0] = E_eff * A * e[0];
  s[1] = 2.0 * G_eff * Ia * e[1];
  s[2] = E_eff * Ia * e[2];
  s[3] = E_eff * Ia * e[3];
  s[4] = kcorr * G_eff * A * e[4];
  s[5] = kcorr * G_eff * A * e[5];
}

void TACSCompositeTubeBeamConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  TacsScalar kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  C[0] = E_eff * A;
  C[6] = 2.0 * G_eff * Ia;
  C[11] = E_eff * Ia;
  C[15] = E_eff * Ia;
  C[18] = kcorr * G_eff * A;
  C[20] = kcorr * G_eff * A;
}

void TACSCompositeTubeBeamConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;
    dfdx[index] +=
        scale * (E_eff * dA * e[0] * psi[0] + 2.0 * G_eff * dIa * e[1] * psi[1] +
                 E_eff * dIa * e[2] * psi[2] + E_eff * dIa * e[3] * psi[3] +
                 kcorr * G_eff * dA * e[4] * psi[4] +
                 kcorr * G_eff * dA * e[5] * psi[5]);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0) / 16.0;
    dfdx[index] +=
        scale * (E_eff * dA * e[0] * psi[0] + 2.0 * G_eff * dIa * e[1] * psi[1] +
                 E_eff * dIa * e[2] * psi[2] + E_eff * dIa * e[3] * psi[3] +
                 kcorr * G_eff * dA * e[4] * psi[4] +
                 kcorr * G_eff * dA * e[5] * psi[5]);
    index++;
  }
}

/*
  Failure modes:
    0: fiber compression  F_1c = -(E11 * eps_1) / X_c  (positive = critical)
    1: fiber tension      F_1t =  (E11 * eps_1) / X_t  (skipped when X_t <= 0)
    2: Euler buckling     F_eu  = Nx / Pcr              (skipped when Kb == 0)

  eps_1 = e[0] + r0*(e[2] + e[3]) — axial strain at outer fibre.
  Uses the same r0 convention as TACSIsoTubeBeamConstitutive: r0 = 0.5*inner + wall.
*/
TacsScalar TACSCompositeTubeBeamConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[]) {
  TacsScalar r0 = 0.5 * inner + wall;
  TacsScalar eps1 = e[0] + r0 * (e[2] + e[3]);

  TacsScalar fail_checks[3];
  int count = 0;
  TacsScalar max_fail = -1e20;

  // Fiber compression
  fail_checks[count] = -(E11 * eps1) / X_c;
  if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
    max_fail = fail_checks[count];
  }
  count++;

  // Fiber tension (optional)
  if (TacsRealPart(X_t) > 0.0) {
    fail_checks[count] = (E11 * eps1) / X_t;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;
  }

  // Euler column buckling
  if (TacsRealPart(buckleLengthFactor) != 0.0) {
    TacsScalar d0 = inner + wall;
    TacsScalar d1 = inner;
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar Leff = buckleLengthFactor * buckleLength;
    TacsScalar Nx = -A * e[0];
    TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    fail_checks[count] = Nx / Pcr;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;
  }

  TacsScalar ks_sum = 0.0;
  for (int i = 0; i < count; i++) {
    ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
  }
  return max_fail + log(ks_sum) / ks_weight;
}

TacsScalar TACSCompositeTubeBeamConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  TacsScalar r0 = 0.5 * inner + wall;
  TacsScalar eps1 = e[0] + r0 * (e[2] + e[3]);

  TacsScalar fail_checks[3];
  TacsScalar fail_checks_sens[3][6];
  int count = 0;
  TacsScalar max_fail = -1e20;

  memset(sens, 0, 6 * sizeof(TacsScalar));
  memset(fail_checks_sens, 0, 3 * 6 * sizeof(TacsScalar));

  // Fiber compression: F_1c = -(E11/X_c) * eps1
  //   d(F_1c)/d(e[0]) = -E11/X_c
  //   d(F_1c)/d(e[2]) = d(F_1c)/d(e[3]) = -E11*r0/X_c
  TacsScalar dF1c_de0 = -E11 / X_c;
  TacsScalar dF1c_de23 = -E11 * r0 / X_c;
  fail_checks[count] = -(E11 * eps1) / X_c;
  fail_checks_sens[count][0] = dF1c_de0;
  fail_checks_sens[count][2] = dF1c_de23;
  fail_checks_sens[count][3] = dF1c_de23;
  if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
    max_fail = fail_checks[count];
  }
  count++;

  // Fiber tension (optional)
  if (TacsRealPart(X_t) > 0.0) {
    TacsScalar dF1t_de0 = E11 / X_t;
    TacsScalar dF1t_de23 = E11 * r0 / X_t;
    fail_checks[count] = (E11 * eps1) / X_t;
    fail_checks_sens[count][0] = dF1t_de0;
    fail_checks_sens[count][2] = dF1t_de23;
    fail_checks_sens[count][3] = dF1t_de23;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;
  }

  // Euler column buckling — only e[0] contributes
  if (TacsRealPart(buckleLengthFactor) != 0.0) {
    TacsScalar d0 = inner + wall;
    TacsScalar d1 = inner;
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar Leff = buckleLengthFactor * buckleLength;
    TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    TacsScalar Nx = -A * e[0];
    fail_checks[count] = Nx / Pcr;
    fail_checks_sens[count][0] = -A / Pcr;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;
  }

  TacsScalar ks_sum = 0.0;
  for (int ii = 0; ii < count; ii++) {
    TacsScalar val = exp(ks_weight * (fail_checks[ii] - max_fail));
    ks_sum += val;
    for (int kk = 0; kk < NUM_STRESSES; kk++) {
      sens[kk] += val * fail_checks_sens[ii][kk];
    }
  }
  for (int kk = 0; kk < NUM_STRESSES; kk++) {
    sens[kk] /= ks_sum;
  }
  return max_fail + log(ks_sum) / ks_weight;
}

void TACSCompositeTubeBeamConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], int dvLen, TacsScalar dfdx[]) {
  int dvNums[2];
  dvNums[0] = innerDV;
  dvNums[1] = wallDV;
  int index = 0;

  for (int dv_index = 0; dv_index < 2; dv_index++) {
    int dvNum = dvNums[dv_index];
    if (dvNum < 0) {
      continue;
    }

    TacsScalar dinner = (dvNum == innerDV) ? 1.0 : 0.0;
    TacsScalar dwall = (dvNum == wallDV) ? 1.0 : 0.0;
    TacsScalar r0 = 0.5 * inner + wall;
    TacsScalar dr0 = 0.5 * dinner + dwall;
    TacsScalar eps1 = e[0] + r0 * (e[2] + e[3]);
    TacsScalar deps1 = dr0 * (e[2] + e[3]);

    TacsScalar fail_checks[3];
    TacsScalar fail_checks_dsens[3];
    int count = 0;
    TacsScalar max_fail = -1e20;

    // Fiber compression: d(F_1c)/d(dv) = -(E11/X_c) * deps1
    fail_checks[count] = -(E11 * eps1) / X_c;
    fail_checks_dsens[count] = -(E11 / X_c) * deps1;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;

    // Fiber tension (optional)
    if (TacsRealPart(X_t) > 0.0) {
      fail_checks[count] = (E11 * eps1) / X_t;
      fail_checks_dsens[count] = (E11 / X_t) * deps1;
      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }
      count++;
    }

    // Euler column buckling
    if (TacsRealPart(buckleLengthFactor) != 0.0) {
      TacsScalar d0 = inner + wall;
      TacsScalar d1 = inner;
      TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
      TacsScalar Ia =
          M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
      TacsScalar dA = M_PI * (wall * dinner + d0 * dwall) / 2.0;
      TacsScalar dIa = M_PI / 16.0 * ((d0 * d0 * d0 - d1 * d1 * d1) * dinner +
                                       d0 * d0 * d0 * dwall);
      TacsScalar Leff = buckleLengthFactor * buckleLength;
      TacsScalar Nx = -A * e[0];
      TacsScalar dNx = -dA * e[0];
      TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
      TacsScalar dPcr = M_PI * M_PI * dIa / (Leff * Leff);
      fail_checks[count] = Nx / Pcr;
      fail_checks_dsens[count] = dNx / Pcr - Nx * dPcr / (Pcr * Pcr);
      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }
      count++;
    }

    TacsScalar ks_sum = 0.0;
    for (int i = 0; i < count; i++) {
      ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
    }
    for (int i = 0; i < count; i++) {
      dfdx[index] += scale * exp(ks_weight * (fail_checks[i] - max_fail)) *
                     fail_checks_dsens[i] / ks_sum;
    }
    index++;
  }
}

TacsScalar TACSCompositeTubeBeamConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return inner;
  } else if (index == 1) {
    return wall;
  }
  return 0.0;
}

const char *TACSCompositeTubeBeamConstitutive::constName =
    "TACSCompositeTubeBeamConstitutive";

const char *TACSCompositeTubeBeamConstitutive::getObjectName() {
  return constName;
}
