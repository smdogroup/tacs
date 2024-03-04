#include "TACSIsoRectangleBeamConstitutive.h"

TACSIsoRectangleBeamConstitutive::TACSIsoRectangleBeamConstitutive(
    TACSMaterialProperties* properties, TacsScalar _width,
    TacsScalar _thickness, int _width_num, int _thickness_num,
    TacsScalar _lb_width, TacsScalar _ub_width, TacsScalar _lb_thickness,
    TacsScalar _ub_thickness, TacsScalar _w_offset, TacsScalar _t_offset) {
  props = properties;
  props->incref();

  thickness = _thickness;
  width = _width;
  thickness_num = _thickness_num;
  width_num = _width_num;
  lb_thickness = _lb_thickness;
  ub_thickness = _ub_thickness;
  lb_width = _lb_width;
  ub_width = _ub_width;
  t_offset = _t_offset;
  w_offset = _w_offset;

  ks_weight = 100.0;
}

TACSIsoRectangleBeamConstitutive::~TACSIsoRectangleBeamConstitutive() {
  props->decref();
}

int TACSIsoRectangleBeamConstitutive::getDesignVarNums(int elemIndex, int dvLen,
                                                       int dvNums[]) {
  int index = 0;
  if (width_num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = width_num;
    }
    index++;
  }
  if (thickness_num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = thickness_num;
    }
    index++;
  }
  return index;
}

int TACSIsoRectangleBeamConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                    const TacsScalar dvs[]) {
  int index = 0;
  if (width_num >= 0) {
    width = dvs[index];
    index++;
  }
  if (thickness_num >= 0) {
    thickness = dvs[index];
    index++;
  }
  return index;
}

int TACSIsoRectangleBeamConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                    TacsScalar dvs[]) {
  int index = 0;
  if (width_num >= 0) {
    dvs[index] = width;
    index++;
  }
  if (thickness_num >= 0) {
    dvs[index] = thickness;
    index++;
  }
  return index;
}

int TACSIsoRectangleBeamConstitutive::getDesignVarRange(int elemIndex,
                                                        int dvLen,
                                                        TacsScalar lb[],
                                                        TacsScalar ub[]) {
  int index = 0;
  if (width_num >= 0) {
    lb[index] = lb_width;
    ub[index] = ub_width;
    index++;
  }
  if (thickness_num >= 0) {
    lb[index] = lb_thickness;
    ub[index] = ub_thickness;
    index++;
  }
  return index;
}

void TACSIsoRectangleBeamConstitutive::evalMassMoments(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       TacsScalar moments[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar A = thickness * width;
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar Iy = 1.0 / 12.0 * thickness * w3 + delta_z * delta_z * A;
  TacsScalar Iz = 1.0 / 12.0 * width * t3 + delta_y * delta_y * A;
  TacsScalar Iyz = -delta_y * delta_z * A;

  moments[0] = rho * A;
  moments[1] = rho * delta_y * A;  // centroid offset y?
  moments[2] = rho * delta_z * A;  // centroid offset z?
  moments[3] = rho * Iz;
  moments[4] = rho * Iy;
  moments[5] = -rho * Iyz;  // -rho * Iyz
}

void TACSIsoRectangleBeamConstitutive::addMassMomentsDVSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar scale[], int dvLen, TacsScalar dfdx[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar A = thickness * width;
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar Iy = 1.0 / 12.0 * thickness * w3 + delta_z * delta_z * A;
  TacsScalar Iz = 1.0 / 12.0 * width * t3 + delta_y * delta_y * A;
  TacsScalar Iyz = -delta_y * delta_z * A;

  int index = 0;
  if (width_num >= 0) {
    TacsScalar dA = A / width;
    TacsScalar dAz = delta_z * dA + w_offset * A;
    TacsScalar dAy = delta_y * dA;
    TacsScalar dIy = 0.25 * thickness * width * width +
                     2.0 * w_offset * delta_z * A + delta_z * delta_z * dA;
    TacsScalar dIz = 1.0 / 12.0 * t3 + delta_y * delta_y * dA;
    TacsScalar dIyz = -delta_y * w_offset * A - delta_y * delta_z * dA;

    dfdx[index] += rho * (scale[0] * dA + scale[1] * dAy + scale[2] * dAz +
                          scale[3] * dIz + scale[4] * dIy - scale[5] * dIyz);
    index++;
  }
  if (thickness_num >= 0) {
    TacsScalar dA = A / thickness;
    TacsScalar dAy = delta_y * dA + t_offset * A;
    TacsScalar dAz = delta_z * dA;
    TacsScalar dIy = 1.0 / 12.0 * w3 + delta_z * delta_z * dA;
    TacsScalar dIz = 0.25 * thickness * thickness * width +
                     2.0 * t_offset * delta_y * A + delta_y * delta_y * dA;
    TacsScalar dIyz = -delta_z * t_offset * A - delta_y * delta_z * dA;

    dfdx[index] += rho * (scale[0] * dA + scale[1] * dAy + scale[2] * dAz +
                          scale[3] * dIz + scale[4] * dIy - scale[5] * dIyz);
    index++;
  }
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[]) {
  if (props) {
    return props->getSpecificHeat();
  }
  return 0.0;
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalDensity(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar A = thickness * width;

  return rho * A;
}

void TACSIsoRectangleBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  TacsScalar rho = props->getDensity();

  if (width_num >= 0) {
    TacsScalar dA = thickness;
    dfdx[index] += scale * rho * dA;
    index++;
  }
  if (thickness_num >= 0) {
    TacsScalar dA = width;
    dfdx[index] += scale * rho * dA;
    index++;
  }
}

void TACSIsoRectangleBeamConstitutive::evalStress(int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar e[],
                                                  TacsScalar s[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar A = thickness * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar Iy = 1.0 / 12.0 * thickness * w3 + delta_z * delta_z * A;
  TacsScalar Iz = 1.0 / 12.0 * width * t3 + delta_y * delta_y * A;
  TacsScalar Iyz = -delta_y * delta_z * A;
  // Torsion constant for rectangle
  TacsScalar a = width / 2.0;
  TacsScalar b = thickness / 2.0;
  TacsScalar b3 = b * b * b;
  TacsScalar b4 = b * b3;
  TacsScalar a4 = a * a * a * a;
  TacsScalar J = a * b3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0));

  s[0] = E * A * (e[0] - delta_y * e[2] - delta_z * e[3]);
  s[1] = G * J * e[1];
  s[2] = E * (Iz * e[2] - delta_y * A * e[0] - Iyz * e[3]);
  s[3] = E * (Iy * e[3] - delta_z * A * e[0] - Iyz * e[2]);
  s[4] = kcorr * G * A * e[4];
  s[5] = kcorr * G * A * e[5];
}

void TACSIsoRectangleBeamConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar A = thickness * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar Iy = 1.0 / 12.0 * thickness * w3 + delta_z * delta_z * A;
  TacsScalar Iz = 1.0 / 12.0 * width * t3 + delta_y * delta_y * A;
  TacsScalar Iyz = -delta_y * delta_z * A;
  // Torsion constant for rectangle
  TacsScalar a = width / 2.0;
  TacsScalar b = thickness / 2.0;
  TacsScalar b3 = b * b * b;
  TacsScalar b4 = b * b3;
  TacsScalar a4 = a * a * a * a;
  TacsScalar J = a * b3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0));

  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  C[0] = E * A;
  C[2] = -E * delta_y * A;
  C[3] = -E * delta_z * A;
  C[6] = G * J;
  C[11] = E * Iz;
  C[12] = -E * Iyz;
  C[15] = E * Iy;
  C[18] = kcorr * G * A;
  C[20] = kcorr * G * A;
}

void TACSIsoRectangleBeamConstitutive::addStressDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[]) {
  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);

  TacsScalar G = 0.5 * E / (1.0 + nu);
  TacsScalar kcorr = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);
  TacsScalar A = thickness * width;
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar a = width / 2.0;
  TacsScalar b = thickness / 2.0;
  TacsScalar b3 = b * b * b;
  TacsScalar b4 = b * b3;
  TacsScalar a4 = a * a * a * a;

  int index = 0;
  if (width_num >= 0) {
    TacsScalar ddelta_z = w_offset;
    TacsScalar dA = thickness;
    TacsScalar dIz =
        1.0 / 12.0 * thickness * thickness * thickness + delta_y * delta_y * dA;
    TacsScalar dIy = width * width * thickness / 4.0 +
                     2.0 * ddelta_z * delta_z * A + delta_z * delta_z * dA;
    TacsScalar dIyz = -delta_y * ddelta_z * A + -delta_y * delta_z * dA;
    TacsScalar da = 0.5;
    TacsScalar da4 = 4.0 * a * a * a * da;
    TacsScalar dJ =
        da * b3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0)) +
        a * b3 * (3.36 * b / a * da / a * (1.0 - b4 / a4 / 12.0)) +
        a * b3 * (3.36 * b / a * (-b4 / a4 / 12.0 * da4 / a4));

    dfdx[index] +=
        scale *
        (E *
             (dA * e[0] - delta_y * dA * e[2] -
              (delta_z * dA + ddelta_z * A) * e[3]) *
             psi[0] +
         G * dJ * e[1] * psi[1] +
         E * (dIz * e[2] - dIyz * e[3] - delta_y * dA * e[0]) * psi[2] +
         E * (dIy * e[3] - dIyz * e[2] - (delta_z * dA + ddelta_z * A) * e[0]) *
             psi[3] +
         kcorr * G * dA * e[4] * psi[4] + kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
  if (thickness_num >= 0) {
    TacsScalar ddelta_y = t_offset;
    TacsScalar dA = width;
    TacsScalar dIz = width * thickness * thickness / 4.0 +
                     2.0 * ddelta_y * delta_y * A + delta_y * delta_y * dA;
    TacsScalar dIy = width * width * width / 12.0 + delta_z * delta_z * dA;
    TacsScalar dIyz = -ddelta_y * delta_z * A + -delta_y * delta_z * dA;
    TacsScalar db = 0.5;
    TacsScalar db3 = 3.0 * b * b * db;
    TacsScalar db4 = 4.0 * b3 * db;
    TacsScalar dJ =
        a * db3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0)) +
        a * b3 * (-3.36 * db / a * (1.0 - b4 / a4 / 12.0)) +
        a * b3 * (-3.36 * b / a * (-db4 / a4 / 12.0));

    dfdx[index] +=
        scale *
        (E *
             (dA * e[0] - (delta_y * dA + ddelta_y * A) * e[2] -
              delta_z * dA * e[3]) *
             psi[0] +
         G * dJ * e[1] * psi[1] +
         E * (dIz * e[2] - dIyz * e[3] - (delta_y * dA + ddelta_y * A) * e[0]) *
             psi[2] +
         E * (dIy * e[3] - dIyz * e[2] - delta_z * dA * e[0]) * psi[3] +
         kcorr * G * dA * e[4] * psi[4] + kcorr * G * dA * e[5] * psi[5]);
    index++;
  }
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalFailure(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[],
                                                         const TacsScalar e[]) {
  // Check the cross-section for failure at the four corners
  TacsScalar e0[6], s0[6];    // 3D stress, NOT beam stresses
  TacsScalar fail_checks[4];  // fail value a four corners
  TacsScalar max_fail = -1e20, ks_sum = 0.0;
  TacsScalar y_lim[] = {-(0.5 + t_offset) * thickness,
                        (0.5 - t_offset) * thickness};
  TacsScalar z_lim[] = {-(0.5 + w_offset) * width, (0.5 - w_offset) * width};

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++, count++) {
      // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
      e0[0] = e[0] + y_lim[i] * e[2] + z_lim[j] * e[3];  // ex
      e0[1] = 0.0;
      e0[2] = 0.0;
      e0[3] = 0.0;
      e0[4] = e[5] + (y_lim[i] + t_offset * thickness) * e[1];
      e0[5] = e[4] - (z_lim[j] + w_offset * width) * e[1];

      // Compute the stress
      props->evalStress3D(e0, s0);

      // Compute the von Mises stress
      fail_checks[count] = props->vonMisesFailure3D(s0);

      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }
    }
  }

  for (int i = 0; i < count; i++) {
    ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
  }

  return max_fail + log(ks_sum) / ks_weight;
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  // Check the cross-section for failure at the four corners
  TacsScalar e0[6], s0[6];
  TacsScalar dfde0[6], dfds0[6];
  TacsScalar fail_checks[4];
  TacsScalar fail_checks_sens[4][6];
  TacsScalar max_fail = -1e20, ks_sum = 0.0;
  TacsScalar y_lim[] = {-(0.5 + t_offset) * thickness,
                        (0.5 - t_offset) * thickness};
  TacsScalar z_lim[] = {-(0.5 + w_offset) * width, (0.5 - w_offset) * width};

  memset(sens, 0, 6 * sizeof(TacsScalar));

  int count = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++, count++) {
      // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
      e0[0] = e[0] + y_lim[i] * e[2] + z_lim[j] * e[3];  // ex
      e0[1] = 0.0;
      e0[2] = 0.0;
      e0[3] = 0.0;
      e0[4] = e[5] + (y_lim[i] + t_offset * thickness) * e[1];
      e0[5] = e[4] - (z_lim[j] + w_offset * width) * e[1];

      // Compute the stress
      props->evalStress3D(e0, s0);

      // Compute the von Mises stress
      fail_checks[count] = props->vonMisesFailure3DStressSens(s0, dfds0);

      props->evalStress3D(dfds0, dfde0);
      fail_checks_sens[count][0] = dfde0[0];
      fail_checks_sens[count][2] = y_lim[i] * dfde0[0];
      fail_checks_sens[count][3] = z_lim[j] * dfde0[0];
      fail_checks_sens[count][1] =
          (y_lim[i] + t_offset * thickness) * dfde0[4] -
          (z_lim[j] + w_offset * width) * dfde0[5];
      fail_checks_sens[count][5] = dfde0[4];
      fail_checks_sens[count][4] = dfde0[5];

      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }
    }
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

void TACSIsoRectangleBeamConstitutive::addFailureDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar e[], int dvLen, TacsScalar dfdx[]) {
  int dvNums[2];
  dvNums[0] = width_num;
  dvNums[1] = thickness_num;
  int index = 0;
  for (int dv_index = 0; dv_index < 2; dv_index++) {
    int dvNum = dvNums[dv_index];
    // Check the cross-section for failure at the four corners
    TacsScalar e0[6], s0[6];
    TacsScalar e0d[6], s0d[6];
    TacsScalar fail_checks[4];
    TacsScalar fail_checks_sens[4];
    TacsScalar max_fail = -1e20, ks_sum = 0.0;
    TacsScalar y_lim[] = {-(0.5 + t_offset) * thickness,
                          (0.5 - t_offset) * thickness};
    TacsScalar z_lim[] = {-(0.5 + w_offset) * width, (0.5 - w_offset) * width};
    TacsScalar dwidth = 0.0, dthickness = 0.0;

    if (dvNum == width_num) {
      dwidth = 1.0;
    }
    if (dvNum == thickness_num) {
      dthickness = 1.0;
    }

    TacsScalar dy_lim[] = {-(0.5 + t_offset) * dthickness,
                           (0.5 - t_offset) * dthickness};
    TacsScalar dz_lim[] = {-(0.5 + w_offset) * dwidth,
                           (0.5 - w_offset) * dwidth};

    int count = 0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++, count++) {
        // Compute the combined strain state e0 = [ex, ey, ez, gyz, gxz, gxy]
        e0[0] = e[0] + y_lim[i] * e[2] + z_lim[j] * e[3];  // ex
        e0[1] = 0.0;
        e0[2] = 0.0;
        e0[3] = 0.0;
        e0[4] = e[5] + (y_lim[i] + t_offset * thickness) * e[1];
        e0[5] = e[4] - (z_lim[j] + w_offset * width) * e[1];

        // Compute the stress
        props->evalStress3D(e0, s0);

        // Compute the von Mises stress
        fail_checks[count] = props->vonMisesFailure3DStressSens(s0, s0d);
        props->evalStress3D(s0d, e0d);

        fail_checks_sens[count] =
            e0d[0] * (dy_lim[i] * e[2] + dz_lim[j] * e[3]) +
            (e0d[4] * (dy_lim[i] + t_offset * dthickness) -
             e0d[5] * (dz_lim[j] + w_offset * dwidth)) *
                e[1];

        if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
          max_fail = fail_checks[count];
        }
      }
    }

    for (int i = 0; i < count; i++) {
      ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
    }

    if (dvNum >= 0) {
      for (int i = 0; i < count; i++) {
        dfdx[index] += scale * exp(ks_weight * (fail_checks[i] - max_fail)) *
                       fail_checks_sens[i] / ks_sum;
      }
      index++;
    }
  }
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return width;
  } else if (index == 1) {
    return thickness;
  }
  return 0.0;
}

const char* TACSIsoRectangleBeamConstitutive::constName =
    "TACSIsoRectangleBeamConstitutive";

/*
  Return the constitutive name
*/
const char* TACSIsoRectangleBeamConstitutive::getObjectName() {
  return constName;
}
