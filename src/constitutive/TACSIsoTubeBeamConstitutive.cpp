#include "TACSIsoTubeBeamConstitutive.h"

TACSIsoTubeBeamConstitutive::TACSIsoTubeBeamConstitutive(
    TACSMaterialProperties* properties, TacsScalar inner_init,
    TacsScalar wall_init, int inner_dv, int wall_dv, TacsScalar inner_lb,
    TacsScalar inner_ub, TacsScalar wall_lb, TacsScalar wall_ub) {
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

  moments[0] = rho * A;
  moments[1] = 0.0;
  moments[2] = 0.0;
  moments[3] = rho * Ia;
  moments[4] = rho * Ia;
  moments[5] = 0.0;
}

void TACSIsoTubeBeamConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar rho = props->getDensity();
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

  return rho * A;
}

void TACSIsoTubeBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  TacsScalar d0 = inner + wall;
  TacsScalar rho = props->getDensity();

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

void TACSIsoTubeBeamConstitutive::evalStress(int elemIndex, const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar e[],
                                             TacsScalar s[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

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

  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 2.0 * (1.0 + nu) / (4.0 + 3.0 * nu);
  TacsScalar d0 = inner + wall;
  TacsScalar d1 = inner;

  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dA = M_PI * wall / 2.0;
    TacsScalar dIa = M_PI * ((d0 * d0 * d0) - (d1 * d1 * d1)) / 16.0;

    dfdx[index] +=
        scale *
        (E * dA * e[0] * psi[0] + 2.0 * G * dIa * e[1] * psi[1] +
         E * dIa * e[2] * psi[2] + E * dIa * e[3] * psi[3] +
         kcorr * G * dA * e[4] * psi[4] + kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dA = M_PI * d0 / 2.0;
    TacsScalar dIa = M_PI * (d0 * d0 * d0) / 16.0;

    dfdx[index] +=
        scale *
        (E * dA * e[0] * psi[0] + 2.0 * G * dIa * e[1] * psi[1] +
         E * dIa * e[2] * psi[2] + E * dIa * e[3] * psi[3] +
         kcorr * G * dA * e[4] * psi[4] + kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailure(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar e[]) {
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3];  // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  return props->vonMisesFailure3D(s0);
}

TacsScalar TACSIsoTubeBeamConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3];  // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  TacsScalar s0d[6];
  TacsScalar fail = props->vonMisesFailure3DStressSens(s0, s0d);

  TacsScalar e0d[6];
  props->evalStress3D(s0d, e0d);

  sens[0] = e0d[0];
  sens[2] = r0 * e0d[0];
  sens[3] = r0 * e0d[0];
  sens[1] = r0 * e0d[3];
  sens[5] = e0d[4];
  sens[4] = e0d[5];

  return fail;
}

void TACSIsoTubeBeamConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], int dvLen, TacsScalar dfdx[]) {
  // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
  TacsScalar e0[6], s0[6];
  TacsScalar r0 = 0.5 * inner + wall;

  e0[0] = e[0] + r0 * e[2] + r0 * e[3];  // ex
  e0[1] = 0.0;
  e0[2] = 0.0;
  e0[3] = r0 * e[1];
  e0[4] = e[5];
  e0[5] = e[4];

  // Compute the stress
  props->evalStress3D(e0, s0);

  // Compute the von Mises stress
  TacsScalar s0d[6];
  props->vonMisesFailure3DStressSens(s0, s0d);

  TacsScalar e0d[6];
  props->evalStress3D(s0d, e0d);

  int index = 0;
  if (innerDV >= 0) {
    TacsScalar dr0 = 0.5;
    dfdx[index] += scale * (e0d[0] * dr0 * (e[2] + e[3]) + e0d[3] * dr0 * e[1]);
    index++;
  }
  if (wallDV >= 0) {
    TacsScalar dr0 = 1.0;
    dfdx[index] += scale * (e0d[0] * dr0 * (e[2] + e[3]) + e0d[3] * dr0 * e[1]);

    index++;
  }
}

TacsScalar TACSIsoTubeBeamConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return inner;
  } else if (index == 1) {
    return wall;
  }
  return 0.0;
}

const char* TACSIsoTubeBeamConstitutive::constName =
    "TACSIsoTubeBeamConstitutive";

/*
  Return the constitutive name
*/
const char* TACSIsoTubeBeamConstitutive::getObjectName() { return constName; }
