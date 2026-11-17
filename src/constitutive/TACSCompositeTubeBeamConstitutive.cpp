#include "TACSCompositeTubeBeamConstitutive.h"

#include <math.h>
#include <string.h>

TACSCompositeTubeBeamConstitutive::TACSCompositeTubeBeamConstitutive(
    TacsScalar E11_in, TacsScalar E22, TacsScalar G12, TacsScalar nu12,
    TacsScalar rho_in, TacsScalar X_c_in, TacsScalar X_t_in,
    const TacsScalar *layup_angles, int n_plies, TacsScalar inner_init,
    TacsScalar wall_init, int inner_dv, int wall_dv, TacsScalar inner_lb,
    TacsScalar inner_ub, TacsScalar wall_lb, TacsScalar wall_ub,
    TacsScalar buckle_length, TacsScalar buckle_length_factor, int x_dv,
    TacsScalar p_penalty_in, TacsScalar k_floor_in, TacsScalar eps_m_in) {
  // Store ply-level material properties
  E11 = E11_in;
  rho = rho_in;
  X_c = X_c_in;
  X_t = X_t_in;

  // Compute CLT-smeared effective moduli (mean Qbar over equal-thickness plies)
  TacsScalar const nu21 = nu12 * E22 / E11_in;
  TacsScalar const denom = 1.0 - nu12 * nu21;
  TacsScalar const Q11 = E11_in / denom;
  TacsScalar const Q22 = E22 / denom;
  TacsScalar const Q12 = nu12 * E22 / denom;
  TacsScalar const Q66 = G12;

  TacsScalar sum_Qbar11 = 0.0;
  TacsScalar sum_Qbar12 = 0.0;
  TacsScalar sum_Qbar66 = 0.0;
  for (int k = 0; k < n_plies; k++) {
    TacsScalar const theta = layup_angles[k];
    TacsScalar const m = cos(theta);
    TacsScalar const n = sin(theta);
    TacsScalar const m2 = m * m;
    TacsScalar const n2 = n * n;
    TacsScalar const m4 = m2 * m2;
    TacsScalar const n4 = n2 * n2;
    TacsScalar const m2n2 = m2 * n2;

    TacsScalar const Qbar11 =
        m4 * Q11 + 2.0 * m2n2 * (Q12 + 2.0 * Q66) + n4 * Q22;
    TacsScalar const Qbar12 = m2n2 * (Q11 + Q22 - 4.0 * Q66) + (m4 + n4) * Q12;
    TacsScalar const Qbar66 =
        (m2 - n2) * (m2 - n2) * Q66 + m2n2 * (Q11 + Q22 - 2.0 * Q12);

    sum_Qbar11 += Qbar11;
    sum_Qbar12 += Qbar12;
    sum_Qbar66 += Qbar66;
  }
  E_eff = sum_Qbar11 / (double)n_plies;
  G_eff = sum_Qbar66 / (double)n_plies;
  TacsScalar const A12_mean = sum_Qbar12 / (double)n_plies;
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

  xDV = x_dv;
  x_val = 1.0;
  p_penalty = p_penalty_in;
  k_floor = k_floor_in;
  eps_m = eps_m_in;
}

TACSCompositeTubeBeamConstitutive::~TACSCompositeTubeBeamConstitutive() {}

TacsScalar TACSCompositeTubeBeamConstitutive::getMaterialDensity() const {
  if (xDV < 0) {
    return rho;
  }
  return (eps_m + (1.0 - eps_m) * x_val) * rho;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getMaterialDensityDVSens() const {
  if (xDV < 0) {
    return 0.0;
  }
  return (1.0 - eps_m) * rho;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getTopologyValue() const {
  if (xDV < 0) {
    return 1.0;
  }
  return x_val;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getTopologyValueDVSens() const {
  if (xDV < 0) {
    return 0.0;
  }
  return 1.0;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getYoungsModulus() const {
  if (xDV < 0) {
    return E_eff;
  }
  TacsScalar const stiffness_scale =
      k_floor + (1.0 - k_floor) * pow(x_val, p_penalty);
  return stiffness_scale * E_eff;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getYoungsModulusDVSens() const {
  if (xDV < 0) {
    return 0.0;
  }
  TacsScalar const dstiffness_scale =
      (1.0 - k_floor) * p_penalty * pow(x_val, p_penalty - 1.0);
  return dstiffness_scale * E_eff;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getShearModulus() const {
  if (xDV < 0) {
    return G_eff;
  }
  TacsScalar const stiffness_scale =
      k_floor + (1.0 - k_floor) * pow(x_val, p_penalty);
  return stiffness_scale * G_eff;
}

TacsScalar TACSCompositeTubeBeamConstitutive::getShearModulusDVSens() const {
  if (xDV < 0) {
    return 0.0;
  }
  TacsScalar const dstiffness_scale =
      (1.0 - k_floor) * p_penalty * pow(x_val, p_penalty - 1.0);
  return dstiffness_scale * G_eff;
}

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
  if (xDV >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = xDV;
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
  if (xDV >= 0) {
    x_val = dvs[index];
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
  if (xDV >= 0) {
    dvs[index] = x_val;
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
  if (xDV >= 0) {
    lb[index] = 0.0;
    ub[index] = 1.0;
    index++;
  }
  return index;
}

TacsScalar TACSCompositeTubeBeamConstitutive::evalDensity(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  return getMaterialDensity() * A;
}

void TACSCompositeTubeBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const density = getMaterialDensity();
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar const dA = M_PI * wall;
    dfdx[index] += scale * density * dA;
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar const dA = M_PI * d0;
    dfdx[index] += scale * density * dA;
    index++;
  }
  if (xDV >= 0) {
    TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    dfdx[index] += scale * getMaterialDensityDVSens() * A;
    index++;
  }
}

void TACSCompositeTubeBeamConstitutive::evalMassMoments(int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        TacsScalar moments[]) {
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar const Ia =
      M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
  TacsScalar const density = getMaterialDensity();

  moments[0] = density * A;
  moments[1] = 0.0;
  moments[2] = 0.0;
  moments[3] = density * Ia;
  moments[4] = density * Ia;
  moments[5] = 0.0;
}

void TACSCompositeTubeBeamConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const density = getMaterialDensity();
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar const dA = M_PI * wall;
    TacsScalar const dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;
    dfdx[index] += density * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar const dA = M_PI * d0;
    TacsScalar const dIa = M_PI * (d0 * d0 * d0) / 8.0;
    dfdx[index] += density * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (xDV >= 0) {
    TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar const Ia =
        M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    dfdx[index] += getMaterialDensityDVSens() *
                   (scale[0] * A + scale[3] * Ia + scale[4] * Ia);
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
  TacsScalar const Ep = getYoungsModulus();
  TacsScalar const Gp = getShearModulus();
  TacsScalar const kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar const Ia =
      M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  s[0] = Ep * A * e[0];
  s[1] = 2.0 * Gp * Ia * e[1];
  s[2] = Ep * Ia * e[2];
  s[3] = Ep * Ia * e[3];
  s[4] = kcorr * Gp * A * e[4];
  s[5] = kcorr * Gp * A * e[5];
}

void TACSCompositeTubeBeamConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  TacsScalar const Ep = getYoungsModulus();
  TacsScalar const Gp = getShearModulus();
  TacsScalar const kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar const Ia =
      M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  C[0] = Ep * A;
  C[6] = 2.0 * Gp * Ia;
  C[11] = Ep * Ia;
  C[15] = Ep * Ia;
  C[18] = kcorr * Gp * A;
  C[20] = kcorr * Gp * A;
}

void TACSCompositeTubeBeamConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar const Ep = getYoungsModulus();
  TacsScalar const Gp = getShearModulus();
  TacsScalar const kcorr = 2.0 * (1.0 + nu_eff) / (4.0 + 3.0 * nu_eff);
  TacsScalar const d0 = inner + 2.0 * wall;
  TacsScalar const d1 = inner;
  int index = 0;
  if (innerDV >= 0) {
    TacsScalar const dA = M_PI * wall;
    TacsScalar const dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;
    dfdx[index] +=
        scale *
        (Ep * dA * e[0] * psi[0] + 2.0 * Gp * dIa * e[1] * psi[1] +
         Ep * dIa * e[2] * psi[2] + Ep * dIa * e[3] * psi[3] +
         kcorr * Gp * dA * e[4] * psi[4] + kcorr * Gp * dA * e[5] * psi[5]);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar const dA = M_PI * d0;
    TacsScalar const dIa = M_PI * (d0 * d0 * d0) / 8.0;
    dfdx[index] +=
        scale *
        (Ep * dA * e[0] * psi[0] + 2.0 * Gp * dIa * e[1] * psi[1] +
         Ep * dIa * e[2] * psi[2] + Ep * dIa * e[3] * psi[3] +
         kcorr * Gp * dA * e[4] * psi[4] + kcorr * Gp * dA * e[5] * psi[5]);
    index++;
  }
  if (xDV >= 0) {
    TacsScalar const dEp = getYoungsModulusDVSens();
    TacsScalar const dGp = getShearModulusDVSens();
    TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar const Ia =
        M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    dfdx[index] +=
        scale *
        (dEp * A * e[0] * psi[0] + 2.0 * dGp * Ia * e[1] * psi[1] +
         dEp * Ia * e[2] * psi[2] + dEp * Ia * e[3] * psi[3] +
         kcorr * dGp * A * e[4] * psi[4] + kcorr * dGp * A * e[5] * psi[5]);
    index++;
  }
}

/*
  Failure modes:
    0: fiber compression  F_1c = -(E11 * eps_1) / X_c  (positive = critical)
    1: fiber tension      F_1t =  (E11 * eps_1) / X_t  (skipped when X_t <= 0)
    2: Euler buckling     F_eu  = Nx / Pcr              (skipped when Kb == 0)

  eps_1 = e[0] + r0*(e[2] + e[3]) — axial strain at outer fibre.
  Uses the same r0 convention as TACSIsoTubeBeamConstitutive: r0 = 0.5*inner +
  wall.
*/
TacsScalar TACSCompositeTubeBeamConstitutive::evalFailure(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[]) {
  TacsScalar const r0 = 0.5 * inner + wall;
  TacsScalar const eps1 = e[0] + r0 * (e[2] + e[3]);

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
    TacsScalar const d0 = inner + 2.0 * wall;
    TacsScalar const d1 = inner;
    TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar const Ia =
        M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar const Leff = buckleLengthFactor * buckleLength;
    TacsScalar const Nx = -A * e[0];
    TacsScalar const Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    TacsScalar const buckle_ratio = Nx / Pcr;
    fail_checks[count] = getTopologyValue() * buckle_ratio;
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
  TacsScalar const r0 = 0.5 * inner + wall;
  TacsScalar const eps1 = e[0] + r0 * (e[2] + e[3]);

  TacsScalar fail_checks[3];
  TacsScalar fail_checks_sens[3][6];
  int count = 0;
  TacsScalar max_fail = -1e20;

  memset(sens, 0, 6 * sizeof(TacsScalar));
  memset(fail_checks_sens, 0, 3 * 6 * sizeof(TacsScalar));

  // Fiber compression: F_1c = -(E11/X_c) * eps1
  //   d(F_1c)/d(e[0]) = -E11/X_c
  //   d(F_1c)/d(e[2]) = d(F_1c)/d(e[3]) = -E11*r0/X_c
  TacsScalar const dF1c_de0 = -E11 / X_c;
  TacsScalar const dF1c_de23 = -E11 * r0 / X_c;
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
    TacsScalar const dF1t_de0 = E11 / X_t;
    TacsScalar const dF1t_de23 = E11 * r0 / X_t;
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
    TacsScalar const d0 = inner + 2.0 * wall;
    TacsScalar const d1 = inner;
    TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar const Ia =
        M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar const Leff = buckleLengthFactor * buckleLength;
    TacsScalar const Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    TacsScalar const Nx = -A * e[0];
    TacsScalar const buckle_ratio = Nx / Pcr;
    TacsScalar const x_eff = getTopologyValue();
    fail_checks[count] = x_eff * buckle_ratio;
    fail_checks_sens[count][0] = x_eff * (-A / Pcr);
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }
    count++;
  }

  TacsScalar ks_sum = 0.0;
  for (int ii = 0; ii < count; ii++) {
    TacsScalar const val = exp(ks_weight * (fail_checks[ii] - max_fail));
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
  int dvNums[3];
  dvNums[0] = innerDV;
  dvNums[1] = wallDV;
  dvNums[2] = xDV;
  int index = 0;

  for (int dv_index = 0; dv_index < 3; dv_index++) {
    int const dvNum = dvNums[dv_index];
    if (dvNum < 0) {
      continue;
    }

    TacsScalar const dinner = (dvNum == innerDV) ? 1.0 : 0.0;
    TacsScalar const dwall = (dvNum == wallDV) ? 1.0 : 0.0;
    TacsScalar const r0 = 0.5 * inner + wall;
    TacsScalar const dr0 = 0.5 * dinner + dwall;
    TacsScalar const eps1 = e[0] + r0 * (e[2] + e[3]);
    TacsScalar const deps1 = dr0 * (e[2] + e[3]);

    TacsScalar fail_checks[3];
    TacsScalar fail_checks_dsens[3];
    int count = 0;
    TacsScalar max_fail = -1e20;

    // Fiber compression: d(F_1c)/d(dv) = -(E11/X_c) * deps1
    // (zero for xDV since fiber failure doesn't depend directly on x)
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
      TacsScalar const d0 = inner + 2.0 * wall;
      TacsScalar const d1 = inner;
      TacsScalar const A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
      TacsScalar const Ia =
          M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
      TacsScalar const dA = M_PI * (wall * dinner + d0 * dwall);
      TacsScalar const dIa =
          M_PI / 16.0 *
          ((d0 * d0 * d0 - d1 * d1 * d1) * dinner + 2.0 * d0 * d0 * d0 * dwall);
      TacsScalar const Leff = buckleLengthFactor * buckleLength;
      TacsScalar const Nx = -A * e[0];
      TacsScalar const dNx = -dA * e[0];
      TacsScalar const Pcr = M_PI * M_PI * Ia / (Leff * Leff);
      TacsScalar const dPcr = M_PI * M_PI * dIa / (Leff * Leff);
      TacsScalar const buckle_ratio = Nx / Pcr;
      TacsScalar const dbuckle_ratio = dNx / Pcr - Nx * dPcr / (Pcr * Pcr);
      TacsScalar const x_eff = getTopologyValue();

      fail_checks[count] = x_eff * buckle_ratio;
      if (dvNum == xDV) {
        fail_checks_dsens[count] = getTopologyValueDVSens() * buckle_ratio;
      } else {
        fail_checks_dsens[count] = x_eff * dbuckle_ratio;
      }
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
  } else if (index == 2) {
    return x_val;
  }
  return 0.0;
}

const char *TACSCompositeTubeBeamConstitutive::constName =
    "TACSCompositeTubeBeamConstitutive";

const char *TACSCompositeTubeBeamConstitutive::getObjectName() {
  return constName;
}
