#include "TACSIsoTubeBeamConstitutive.h"

namespace {
// Finite SIMP stiffness floor for numerical regularization so near-void
// members do not create near-singular stiffness matrices. This improves
// conditioning for eigenvalue solves and parallel runs.
inline TacsScalar simp_stiffness_scale(int xDV, TacsScalar x_val,
                                       TacsScalar p_penalty,
                                       TacsScalar k_floor) {
  if (xDV < 0) {
    return 1.0;
  }
  return k_floor + (1.0 - k_floor) * pow(x_val, p_penalty);
}

inline TacsScalar simp_stiffness_scale_derivative(int xDV, TacsScalar x_val,
                                                  TacsScalar p_penalty,
                                                  TacsScalar k_floor) {
  if (xDV < 0) {
    return 0.0;
  }
  return (1.0 - k_floor) * p_penalty * pow(x_val, p_penalty - 1.0);
}
}  // namespace

TACSIsoTubeBeamConstitutive::TACSIsoTubeBeamConstitutive(
    TACSMaterialProperties *properties, TacsScalar inner_init,
    TacsScalar wall_init, int inner_dv, int wall_dv, TacsScalar inner_lb,
    TacsScalar inner_ub, TacsScalar wall_lb, TacsScalar wall_ub,
    TacsScalar buckle_length, int buckle_length_dv,
    TacsScalar buckle_length_factor,
    int x_dv, TacsScalar p_penalty_in, TacsScalar k_floor_in,
    TacsScalar eps_m_in) {
  props = properties;
  props->incref();

  inner = inner_init;
  wall = wall_init;
  innerDV = inner_dv;
  wallDV = wall_dv;
  innerLb = inner_lb;
  innerUb = inner_ub;
  wallLb = wall_lb;
  wallUb = wall_ub;
  buckleLength = buckle_length;
  buckleLengthDV = buckle_length_dv;
  buckleLengthFactor = buckle_length_factor;
  buckleLengthLb = 1e-20;
  buckleLengthUb = 1e20;
  ks_weight = 100.0;
  xDV = x_dv;
  x_val = 1.0;
  p_penalty = p_penalty_in;
  k_floor = k_floor_in;
  eps_m = eps_m_in;
}

TACSIsoTubeBeamConstitutive::~TACSIsoTubeBeamConstitutive() { props->decref(); }

int TACSIsoTubeBeamConstitutive::getDesignVarNums(int elemIndex, int dvLen,
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
  if (buckleLengthDV >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = buckleLengthDV;
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

int TACSIsoTubeBeamConstitutive::setDesignVars(int elemIndex, int dvLen,
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
  if (buckleLengthDV >= 0) {
    buckleLength = dvs[index];
    index++;
  }
  if (xDV >= 0) {
    x_val = dvs[index];
    index++;
  }
  return index;
}

int TACSIsoTubeBeamConstitutive::getDesignVars(int elemIndex, int dvLen,
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
  if (buckleLengthDV >= 0) {
    dvs[index] = buckleLength;
    index++;
  }
  if (xDV >= 0) {
    dvs[index] = x_val;
    index++;
  }
  return index;
}

int TACSIsoTubeBeamConstitutive::getDesignVarRange(int elemIndex, int dvLen,
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
  if (buckleLengthDV >= 0) {
    lb[index] = buckleLengthLb;
    ub[index] = buckleLengthUb;
    index++;
  }
  if (xDV >= 0) {
    lb[index] = 0.0;
    ub[index] = 1.0;
    index++;
  }
  return index;
}

void TACSIsoTubeBeamConstitutive::evalMassMoments(int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  TacsScalar moments[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  TacsScalar rho_scale = (xDV >= 0) ? (eps_m + (1.0 - eps_m) * x_val) : 1.0;
  moments[0] = rho_scale * rho * A;
  moments[1] = 0.0;
  moments[2] = 0.0;
  moments[3] = rho_scale * rho * Ia;
  moments[4] = rho_scale * rho * Ia;
  moments[5] = 0.0;
}

void TACSIsoTubeBeamConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar rho_scale = (xDV >= 0) ? (eps_m + (1.0 - eps_m) * x_val) : 1.0;

  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;

    dfdx[index] += rho_scale * rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0) / 16.0;

    dfdx[index] += rho_scale * rho * (scale[0] * dA + scale[3] * dIa + scale[4] * dIa);
    index++;
  }
  if (xDV >= 0) {
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    dfdx[index] += (1.0 - eps_m) * rho *
                   (scale[0] * A + scale[3] * Ia + scale[4] * Ia);
    index++;
  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalSpecificHeat(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[]) {
  if (props) {
    return props->getSpecificHeat();
  }
  return 0.0;
}

TacsScalar TACSIsoTubeBeamConstitutive::evalDensity(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[]) {
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar rho = props->getDensity();
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar rho_scale = (xDV >= 0) ? (eps_m + (1.0 - eps_m) * x_val) : 1.0;

  return rho_scale * rho * A;
}

void TACSIsoTubeBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar rho = props->getDensity();
  TacsScalar rho_scale = (xDV >= 0) ? (eps_m + (1.0 - eps_m) * x_val) : 1.0;

  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    dfdx[index] += scale * rho_scale * rho * dA;
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    dfdx[index] += scale * rho_scale * rho * dA;
    index++;
  }
  if (xDV >= 0) {
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    dfdx[index] += scale * (1.0 - eps_m) * rho * A;
    index++;
  }
}

void TACSIsoTubeBeamConstitutive::evalStress(int elemIndex, const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar e[],
                                             TacsScalar s[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar stiffness_scale =
      simp_stiffness_scale(xDV, x_val, p_penalty, k_floor);
  E *= stiffness_scale;
  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  s[0] = E * A * e[0];
  s[1] = 2.0 * G * Ia * e[1];
  s[2] = E * Ia * e[2];
  s[3] = E * Ia * e[3];
  s[4] = kcorr * G * A * e[4];
  s[5] = kcorr * G * A * e[5];
}

void TACSIsoTubeBeamConstitutive::evalTangentStiffness(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       TacsScalar C[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar stiffness_scale =
      simp_stiffness_scale(xDV, x_val, p_penalty, k_floor);
  E *= stiffness_scale;
  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;
  TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
  TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  C[0] = E * A;
  C[6] = 2.0 * G * Ia;
  C[11] = E * Ia;
  C[15] = E * Ia;
  C[18] = kcorr * G * A;
  C[20] = kcorr * G * A;
}

void TACSIsoTubeBeamConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar stiffness_scale =
      simp_stiffness_scale(xDV, x_val, p_penalty, k_floor);
  TacsScalar Ep = E * stiffness_scale;
  TacsScalar Gp = 0.5 * Ep / (1.0 + nu);
  TacsScalar kcorr = 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;

  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;

    dfdx[index] +=
        scale *
        (Ep * dA * e[0] * psi[0] + 2.0 * Gp * dIa * e[1] * psi[1] +
         Ep * dIa * e[2] * psi[2] + Ep * dIa * e[3] * psi[3] +
         kcorr * Gp * dA * e[4] * psi[4] + kcorr * Gp * dA * e[5] * psi[5]);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0) / 16.0;

    dfdx[index] +=
        scale *
        (Ep * dA * e[0] * psi[0] + 2.0 * Gp * dIa * e[1] * psi[1] +
         Ep * dIa * e[2] * psi[2] + Ep * dIa * e[3] * psi[3] +
         kcorr * Gp * dA * e[4] * psi[4] + kcorr * Gp * dA * e[5] * psi[5]);
    index++;
  }
  // Skip buckleLengthDV (stiffness doesn't depend on buckle length)
  if (buckleLengthDV >= 0) {
    index++;
  }
  if (xDV >= 0) {
    TacsScalar dstiffness_scale =
        simp_stiffness_scale_derivative(xDV, x_val, p_penalty, k_floor);
    TacsScalar dEp = E * dstiffness_scale;
    TacsScalar dGp = 0.5 * dEp / (1.0 + nu);
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;

    dfdx[index] +=
        scale *
        (dEp * A * e[0] * psi[0] + 2.0 * dGp * Ia * e[1] * psi[1] +
         dEp * Ia * e[2] * psi[2] + dEp * Ia * e[3] * psi[3] +
         kcorr * dGp * A * e[4] * psi[4] + kcorr * dGp * A * e[5] * psi[5]);
    index++;
  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailure(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar e[]) {
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;
  TacsScalar fail_checks[2];
  TacsScalar max_fail = -1e20, ks_sum = 0.0;
  int count = 1;

  // Von Mises failure at outer fibre
  e0[0] = e[0] + r0 * e[2] + r0 * e[3];
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];
  props->evalStress3D(e0, s0);
  fail_checks[0] = props->vonMisesFailure3D(s0);
  if (TacsRealPart(fail_checks[0]) > TacsRealPart(max_fail)) {
    max_fail = fail_checks[0];
  }

  // Euler column buckling (pin-pin, E cancels in Nx/Pcr)
  if (TacsRealPart(buckleLengthFactor) != 0.0) {
    TacsScalar d0 = inner + wall;
    TacsScalar d1 = inner;
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar Leff = buckleLengthFactor * buckleLength;
    TacsScalar Nx = -A * e[0];
    TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    TacsScalar buckle_ratio = Nx / Pcr;
    fail_checks[1] = (xDV >= 0) ? x_val * buckle_ratio : buckle_ratio;
    if (TacsRealPart(fail_checks[1]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[1];
    }
    count = 2;
  }

  for (int i = 0; i < count; i++) {
    ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
  }
  return max_fail + log(ks_sum) / ks_weight;
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  TacsScalar e0[6], s0[6], s0d[6], e0d[6];
  TacsScalar r0 = 0.5 * inner + wall;
  TacsScalar fail_checks[2];
  TacsScalar fail_checks_sens[2][6];
  TacsScalar max_fail = -1e20, ks_sum = 0.0;
  int count = 1;

  memset(sens, 0, 6 * sizeof(TacsScalar));
  memset(fail_checks_sens, 0, 2 * 6 * sizeof(TacsScalar));

  // Von Mises failure at outer fibre
  e0[0] = e[0] + r0 * e[2] + r0 * e[3];
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];
  props->evalStress3D(e0, s0);
  fail_checks[0] = props->vonMisesFailure3DStressSens(s0, s0d);
  props->evalStress3D(s0d, e0d);
  fail_checks_sens[0][0] = e0d[0];
  fail_checks_sens[0][2] = r0 * e0d[0];
  fail_checks_sens[0][3] = r0 * e0d[0];
  fail_checks_sens[0][1] = r0 * e0d[3];
  fail_checks_sens[0][5] = e0d[4];
  fail_checks_sens[0][4] = e0d[5];
  if (TacsRealPart(fail_checks[0]) > TacsRealPart(max_fail)) {
    max_fail = fail_checks[0];
  }

  // Euler column buckling strain sensitivity: only e[0] (axial) contributes
  if (TacsRealPart(buckleLengthFactor) != 0.0) {
    TacsScalar d0 = inner + wall;
    TacsScalar d1 = inner;
    TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
    TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
    TacsScalar Leff = buckleLengthFactor * buckleLength;
    TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
    TacsScalar Nx = -A * e[0];
    TacsScalar buckle_ratio = Nx / Pcr;
    TacsScalar x_relax = (xDV >= 0) ? x_val : 1.0;
    fail_checks[1] = x_relax * buckle_ratio;
    fail_checks_sens[1][0] = x_relax * (-A / Pcr);
    if (TacsRealPart(fail_checks[1]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[1];
    }
    count = 2;
  }

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

void TACSIsoTubeBeamConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], int dvLen, TacsScalar dfdx[]) {
  int dvNums[4];
  dvNums[0] = innerDV;
  dvNums[1] = wallDV;
  dvNums[2] = buckleLengthDV;
  dvNums[3] = xDV;
  int index = 0;

  TacsScalar x_relax = (xDV >= 0) ? x_val : 1.0;

  for (int dv_index = 0; dv_index < 4; dv_index++) {
    int dvNum = dvNums[dv_index];
    if (dvNum < 0) {
      continue;
    }

    TacsScalar e0[6], s0[6], s0d[6], e0d[6];
    TacsScalar fail_checks[2];
    TacsScalar fail_checks_sens[2];
    TacsScalar max_fail = -1e20, ks_sum = 0.0;
    int count = 1;

    TacsScalar dinner = (dvNum == innerDV) ? 1.0 : 0.0;
    TacsScalar dwall = (dvNum == wallDV) ? 1.0 : 0.0;
    TacsScalar r0 = 0.5 * inner + wall;
    TacsScalar dr0 = 0.5 * dinner + dwall;

    // Von Mises DV sensitivity (zero for x and buckleLengthDV)
    e0[0] = e[0] + r0 * e[2] + r0 * e[3];
    e0[1] = 0.0;
    e0[2] = 0.0;
    e0[3] = r0 * e[1];
    e0[4] = e[5];
    e0[5] = e[4];
    props->evalStress3D(e0, s0);
    fail_checks[0] = props->vonMisesFailure3DStressSens(s0, s0d);
    props->evalStress3D(s0d, e0d);
    fail_checks_sens[0] =
        e0d[0] * dr0 * (e[2] + e[3]) + e0d[3] * dr0 * e[1];
    if (TacsRealPart(fail_checks[0]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[0];
    }

    // Euler buckling DV sensitivity
    if (TacsRealPart(buckleLengthFactor) != 0.0) {
      TacsScalar d0 = inner + wall;
      TacsScalar d1 = inner;
      TacsScalar A = M_PI * ((d0 * d0) - (d1 * d1)) / 4.0;
      TacsScalar Ia = M_PI * ((d0 * d0 * d0 * d0) - (d1 * d1 * d1 * d1)) / 64.0;
      TacsScalar dA = M_PI * (wall * dinner + d0 * dwall) / 2.0;
      TacsScalar dIa =
          M_PI / 16.0 * ((d0 * d0 * d0 - d1 * d1 * d1) * dinner +
                         d0 * d0 * d0 * dwall);
      TacsScalar dL = (dvNum == buckleLengthDV) ? 1.0 : 0.0;
      TacsScalar Leff = buckleLengthFactor * buckleLength;
      TacsScalar dLeff = buckleLengthFactor * dL;
      TacsScalar Nx = -A * e[0];
      TacsScalar dNx = -dA * e[0];
      TacsScalar Pcr = M_PI * M_PI * Ia / (Leff * Leff);
      TacsScalar dPcr =
          M_PI * M_PI *
          (dIa / (Leff * Leff) - 2.0 * Ia * dLeff / (Leff * Leff * Leff));
      TacsScalar buckle_ratio = Nx / Pcr;
      TacsScalar dbuckle_ratio = dNx / Pcr - Nx * dPcr / (Pcr * Pcr);

      fail_checks[1] = x_relax * buckle_ratio;
      if (dvNum == xDV) {
        fail_checks_sens[1] = buckle_ratio;
      } else {
        fail_checks_sens[1] = x_relax * dbuckle_ratio;
      }
      if (TacsRealPart(fail_checks[1]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[1];
      }
      count = 2;
    }

    for (int i = 0; i < count; i++) {
      ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
    }
    for (int i = 0; i < count; i++) {
      dfdx[index] += scale * exp(ks_weight * (fail_checks[i] - max_fail)) *
                     fail_checks_sens[i] / ks_sum;
    }
    index++;
  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return inner;
  } else if (index == 1) {
    return wall;
  } else if (index == 2) {
    return buckleLength;
  } else if (index == 3) {
    return x_val;
  }
  return 0.0;
}

const char *TACSIsoTubeBeamConstitutive::constName =
    "TACSIsoTubeBeamConstitutive";

/*
  Return the constitutive name
*/
const char *TACSIsoTubeBeamConstitutive::getObjectName() { return constName; }
