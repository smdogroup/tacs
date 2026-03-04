#include "TACSIsoRectangleBeamConstitutive.h"

constexpr TacsScalar TACSIsoRectangleBeamConstitutive::eps;

TACSIsoRectangleBeamConstitutive::TACSIsoRectangleBeamConstitutive(
    TACSMaterialProperties *properties, TacsScalar _width,
    TacsScalar _thickness, TacsScalar _buckle_length, int _width_num,
    int _thickness_num, int _buckle_length_num, TacsScalar _lb_width,
    TacsScalar _ub_width, TacsScalar _lb_thickness, TacsScalar _ub_thickness,
    TacsScalar _w_offset, TacsScalar _t_offset,
    TacsScalar _buckle_length_factor) {
  props = properties;
  props->incref();

  thickness = _thickness;
  width = _width;
  buckle_length = _buckle_length;
  thickness_num = _thickness_num;
  width_num = _width_num;
  buckle_length_num = _buckle_length_num;
  lb_thickness = _lb_thickness;
  ub_thickness = _ub_thickness;
  lb_width = _lb_width;
  ub_width = _ub_width;
  lb_buckle_length = 1e-20;
  ub_buckle_length = 1e20;
  t_offset = _t_offset;
  w_offset = _w_offset;
  buckle_length_factor = _buckle_length_factor;

  props->getIsotropicProperties(&E, &nu);

  G = 0.5 * E / (1.0 + nu);
  kcorr = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);

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
  if (buckle_length_num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = buckle_length_num;
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
  if (buckle_length_num >= 0) {
    buckle_length = dvs[index];
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
  if (buckle_length_num >= 0) {
    dvs[index] = buckle_length;
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
  if (buckle_length_num >= 0) {
    lb[index] = lb_buckle_length;
    ub[index] = ub_buckle_length;
    index++;
  }
  return index;
}

void TACSIsoRectangleBeamConstitutive::evalMassMoments(int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[],
                                                       TacsScalar moments[]) {
  TacsScalar rho = props->getDensity();
  TacsScalar A = evalArea();
  TacsScalar I[3];
  evalMomentsOfInertia(I);
  TacsScalar Iy = I[0], Iz = I[1], Iyz = I[2];
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;

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
  TacsScalar A = evalArea();
  TacsScalar I[3];
  evalMomentsOfInertia(I);
  TacsScalar Iy = I[0], Iz = I[1], Iyz = I[2];
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;

  int index = 0;
  if (width_num >= 0) {
    TacsScalar dA = evalAreaSens(width_num);
    TacsScalar dAz = delta_z * dA + w_offset * A;
    TacsScalar dAy = delta_y * dA;
    TacsScalar sI[3];
    evalMomentsOfInertiaSens(width_num, sI);
    TacsScalar dIy = sI[0], dIz = sI[1], dIyz = sI[2];

    dfdx[index] += rho * (scale[0] * dA + scale[1] * dAy + scale[2] * dAz +
                          scale[3] * dIz + scale[4] * dIy - scale[5] * dIyz);
    index++;
  }
  if (thickness_num >= 0) {
    TacsScalar dA = evalAreaSens(thickness_num);
    TacsScalar dAy = delta_y * dA + t_offset * A;
    TacsScalar dAz = delta_z * dA;
    TacsScalar sI[3];
    evalMomentsOfInertiaSens(thickness_num, sI);
    TacsScalar dIy = sI[0], dIz = sI[1], dIyz = sI[2];

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
  TacsScalar A = evalArea();

  return rho * A;
}

void TACSIsoRectangleBeamConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[]) {
  int index = 0;
  TacsScalar rho = props->getDensity();

  if (width_num >= 0) {
    TacsScalar dA = evalAreaSens(width_num);
    dfdx[index] += scale * rho * dA;
    index++;
  }
  if (thickness_num >= 0) {
    TacsScalar dA = evalAreaSens(thickness_num);
    dfdx[index] += scale * rho * dA;
    index++;
  }
}

void TACSIsoRectangleBeamConstitutive::evalStress(int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar e[],
                                                  TacsScalar s[]) {
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar A = evalArea();
  TacsScalar I[3];
  evalMomentsOfInertia(I);
  TacsScalar Iy = I[0], Iz = I[1], Iyz = I[2];
  TacsScalar J = evalTorsionalConstant();

  s[0] = E * A * (e[0] - delta_y * e[2] - delta_z * e[3]);
  s[1] = G * J * e[1];
  s[2] = E * (Iz * e[2] - delta_y * A * e[0] - Iyz * e[3]);
  s[3] = E * (Iy * e[3] - delta_z * A * e[0] - Iyz * e[2]);
  s[4] = kcorr * G * A * e[4];
  s[5] = kcorr * G * A * e[5];
}

void TACSIsoRectangleBeamConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar A = evalArea();
  TacsScalar I[3];
  evalMomentsOfInertia(I);
  TacsScalar Iy = I[0], Iz = I[1], Iyz = I[2];
  TacsScalar J = evalTorsionalConstant();

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
  TacsScalar A = evalArea();
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;

  int index = 0;
  if (width_num >= 0) {
    TacsScalar ddelta_z = w_offset;
    TacsScalar dA = evalAreaSens(width_num);
    TacsScalar sI[3];
    evalMomentsOfInertiaSens(width_num, sI);
    TacsScalar dIy = sI[0], dIz = sI[1], dIyz = sI[2];
    TacsScalar dJ = evalTorsionalConstantSens(width_num);

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
    TacsScalar dA = evalAreaSens(thickness_num);
    TacsScalar sI[3];
    evalMomentsOfInertiaSens(thickness_num, sI);
    TacsScalar dIy = sI[0], dIz = sI[1], dIyz = sI[2];
    TacsScalar dJ = evalTorsionalConstantSens(thickness_num);

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
  TacsScalar fail_checks[6];  // fail value a four corners + 2 buckling checks
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

  // Check the buckling failure criteria
  if (buckle_length_factor != 0.0) {
    TacsScalar delta_y = t_offset * thickness;
    TacsScalar delta_z = w_offset * width;
    TacsScalar A = evalArea();
    TacsScalar t3 = thickness * thickness * thickness;
    TacsScalar w3 = width * width * width;
    // NOTE: Iy and Iz are taken about centroid, not offset
    TacsScalar Iy = 1.0 / 12.0 * thickness * w3;
    TacsScalar Iz = 1.0 / 12.0 * width * t3;
    TacsScalar Leff = buckle_length_factor * buckle_length;
    // We don't include young modulus here since it cancels out in the ratio
    // later
    TacsScalar Nx = -A * (e[0] - delta_y * e[2] - delta_z * e[3]);

    TacsScalar Pcr1 = M_PI * M_PI * Iy / (Leff * Leff);
    fail_checks[count] = Nx / Pcr1;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }

    count += 1;

    TacsScalar Pcr2 = M_PI * M_PI * Iz / (Leff * Leff);
    fail_checks[count] = Nx / Pcr2;
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }

    count += 1;
  }

  for (int i = 0; i < count; i++) {
    ks_sum += exp(ks_weight * (fail_checks[i] - max_fail));
  }

  return max_fail + log(ks_sum) / ks_weight;
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalFailureStrainSens(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar e[], TacsScalar sens[]) {
  // Check the cross-section for failure at the four corners + 2 buckling checks
  TacsScalar e0[6], s0[6];
  TacsScalar dfde0[6], dfds0[6];
  TacsScalar fail_checks[6];
  TacsScalar fail_checks_sens[6][6];
  TacsScalar max_fail = -1e20, ks_sum = 0.0;
  TacsScalar y_lim[] = {-(0.5 + t_offset) * thickness,
                        (0.5 - t_offset) * thickness};
  TacsScalar z_lim[] = {-(0.5 + w_offset) * width, (0.5 - w_offset) * width};

  memset(sens, 0, 6 * sizeof(TacsScalar));
  memset(fail_checks_sens, 0, 6 * 6 * sizeof(TacsScalar));

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

  // Check the buckling failure criteria
  if (TacsRealPart(buckle_length_factor) != 0.0) {
    TacsScalar delta_y = t_offset * thickness;
    TacsScalar delta_z = w_offset * width;
    TacsScalar A = evalArea();
    TacsScalar t3 = thickness * thickness * thickness;
    TacsScalar w3 = width * width * width;
    TacsScalar Iy = 1.0 / 12.0 * thickness * w3;
    TacsScalar Iz = 1.0 / 12.0 * width * t3;
    TacsScalar Leff = buckle_length_factor * buckle_length;
    TacsScalar Nx = -A * (e[0] - delta_y * e[2] - delta_z * e[3]);

    TacsScalar Pcr1 = M_PI * M_PI * Iy / (Leff * Leff);
    fail_checks[count] = Nx / Pcr1;
    fail_checks_sens[count][0] = -A / (M_PI * M_PI * Iy / (Leff * Leff));
    fail_checks_sens[count][2] =
        delta_y * A / (M_PI * M_PI * Iy / (Leff * Leff));
    fail_checks_sens[count][3] =
        delta_z * A / (M_PI * M_PI * Iy / (Leff * Leff));
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }

    count += 1;

    TacsScalar Pcr2 = M_PI * M_PI * Iz / (Leff * Leff);
    fail_checks[count] = Nx / Pcr2;
    fail_checks_sens[count][0] = -A / (M_PI * M_PI * Iz / (Leff * Leff));
    fail_checks_sens[count][2] =
        delta_y * A / (M_PI * M_PI * Iz / (Leff * Leff));
    fail_checks_sens[count][3] =
        delta_z * A / (M_PI * M_PI * Iz / (Leff * Leff));
    if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
      max_fail = fail_checks[count];
    }

    count += 1;
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
  int dvNums[3];
  dvNums[0] = width_num;
  dvNums[1] = thickness_num;
  dvNums[2] = buckle_length_num;
  int index = 0;
  for (int dv_index = 0; dv_index < 3; dv_index++) {
    int dvNum = dvNums[dv_index];

    // If the design variable is not being optimized, skip it
    if (dvNum < 0) {
      continue;
    }

    // Check the cross-section for failure at the four corners + 2 buckling
    // checks
    TacsScalar e0[6], s0[6];
    TacsScalar e0d[6], s0d[6];
    TacsScalar fail_checks[6];
    TacsScalar fail_checks_sens[6];
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

    // Check the buckling failure criteria
    if (TacsRealPart(buckle_length_factor) != 0.0) {
      TacsScalar dL = 0.0;
      if (dvNum == buckle_length_num) {
        dL = 1.0;
      }
      TacsScalar delta_y = t_offset * thickness;
      TacsScalar ddelta_y = t_offset * dthickness;
      TacsScalar delta_z = w_offset * width;
      TacsScalar ddelta_z = w_offset * dwidth;
      TacsScalar A = evalArea();
      TacsScalar dA = evalAreaSens(dvNum);
      TacsScalar t3 = thickness * thickness * thickness;
      TacsScalar w3 = width * width * width;
      TacsScalar Iy = 1.0 / 12.0 * thickness * w3;
      TacsScalar Iz = 1.0 / 12.0 * width * t3;
      TacsScalar dIy = 0.0, dIz = 0.0;
      if (dvNum == width_num) {
        dIy = 0.25 * thickness * width * width;
        dIz = 1.0 / 12.0 * t3;
      } else if (dvNum == thickness_num) {
        dIy = 1.0 / 12.0 * w3;
        dIz = 0.25 * thickness * thickness * width;
      }
      TacsScalar Leff = buckle_length_factor * buckle_length;
      TacsScalar dLeff = buckle_length_factor * dL;
      TacsScalar Nx = -A * (e[0] - delta_y * e[2] - delta_z * e[3]);
      TacsScalar dNx = A * (ddelta_y * e[2] + ddelta_z * e[3]) -
                       dA * (e[0] - delta_y * e[2] - delta_z * e[3]);

      TacsScalar Pcr1 = M_PI * M_PI * Iy / (Leff * Leff);
      TacsScalar dPcr1 =
          M_PI * M_PI *
          (dIy / (Leff * Leff) - 2.0 * Iy * dLeff / (Leff * Leff * Leff));
      fail_checks[count] = Nx / Pcr1;
      fail_checks_sens[count] = dNx / Pcr1 - Nx * dPcr1 / (Pcr1 * Pcr1);
      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }

      count += 1;

      TacsScalar Pcr2 = M_PI * M_PI * Iz / (Leff * Leff);
      TacsScalar dPcr2 =
          M_PI * M_PI *
          (dIz / (Leff * Leff) - 2.0 * Iz * dLeff / (Leff * Leff * Leff));
      fail_checks[count] = Nx / Pcr2;
      fail_checks_sens[count] = dNx / Pcr2 - Nx * dPcr2 / (Pcr2 * Pcr2);
      if (TacsRealPart(fail_checks[count]) > TacsRealPart(max_fail)) {
        max_fail = fail_checks[count];
      }

      count += 1;
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

TacsScalar TACSIsoRectangleBeamConstitutive::evalAreaSens(int dvNum) {
  TacsScalar dA = 0.0;
  if (dvNum == width_num) {
    dA = thickness;
  }
  if (dvNum == thickness_num) {
    dA = width;
  }
  return dA;
}

void TACSIsoRectangleBeamConstitutive::evalMomentsOfInertia(
    TacsScalar moments[]) {
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar A = evalArea();
  TacsScalar Iy = 1.0 / 12.0 * thickness * w3 + delta_z * delta_z * A;
  TacsScalar Iz = 1.0 / 12.0 * width * t3 + delta_y * delta_y * A;
  TacsScalar Iyz = -delta_y * delta_z * thickness * width;
  moments[0] = Iy;
  moments[1] = Iz;
  moments[2] = Iyz;
}

void TACSIsoRectangleBeamConstitutive::evalMomentsOfInertiaSens(
    int dvNum, TacsScalar momentsSens[]) {
  TacsScalar delta_y = t_offset * thickness;
  TacsScalar delta_z = w_offset * width;
  TacsScalar A = evalArea();
  TacsScalar t3 = thickness * thickness * thickness;
  TacsScalar w3 = width * width * width;
  TacsScalar dA = evalAreaSens(dvNum);
  momentsSens[0] = 0.0;
  momentsSens[1] = 0.0;
  momentsSens[2] = 0.0;
  if (dvNum == width_num) {
    momentsSens[0] = 0.25 * thickness * width * width +
                     2.0 * w_offset * delta_z * A + delta_z * delta_z * dA;
    momentsSens[1] = 1.0 / 12.0 * t3 + delta_y * delta_y * dA;
    momentsSens[2] = -delta_y * w_offset * A - delta_y * delta_z * dA;
  } else if (dvNum == thickness_num) {
    momentsSens[0] = 1.0 / 12.0 * w3 + delta_z * delta_z * dA;
    momentsSens[1] = 0.25 * thickness * thickness * width +
                     2.0 * t_offset * delta_y * A + delta_y * delta_y * dA;
    momentsSens[2] = -delta_z * t_offset * A - delta_y * delta_z * dA;
  }
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalTorsionalConstant() {
  // Torsion constant for rectangle
  // continuous approx of max(thickness, width), min(thickness, width)
  TacsScalar tpw = thickness + width;
  TacsScalar tmw = thickness - width;
  TacsScalar max_dim = (tpw + std::sqrt(tmw * tmw + eps)) / 2.0;
  TacsScalar min_dim = (tpw - std::sqrt(tmw * tmw + eps)) / 2.0;
  TacsScalar a = max_dim / 2.0;
  TacsScalar b = min_dim / 2.0;
  TacsScalar b3 = b * b * b;
  TacsScalar b4 = b * b3;
  TacsScalar a4 = a * a * a * a;
  TacsScalar J = a * b3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0));
  return J;
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalTorsionalConstantSens(
    int dvNum) {
  TacsScalar tpw = thickness + width;
  TacsScalar tmw = thickness - width;
  TacsScalar max_dim = (tpw + std::sqrt(tmw * tmw + eps)) / 2.0;
  TacsScalar min_dim = (tpw - std::sqrt(tmw * tmw + eps)) / 2.0;
  TacsScalar a = max_dim / 2.0;
  TacsScalar b = min_dim / 2.0;
  TacsScalar b3 = b * b * b;
  TacsScalar b4 = b * b3;
  TacsScalar a4 = a * a * a * a;
  TacsScalar dmax_dim = 0.0;
  TacsScalar dmin_dim = 0.0;
  if (dvNum == width_num) {
    dmax_dim = (1.0 - tmw / std::sqrt(tmw * tmw + eps)) / 2.0;
    dmin_dim = (1.0 + tmw / std::sqrt(tmw * tmw + eps)) / 2.0;
  } else if (dvNum == thickness_num) {
    dmax_dim = (1.0 + tmw / std::sqrt(tmw * tmw + eps)) / 2.0;
    dmin_dim = (1.0 - tmw / std::sqrt(tmw * tmw + eps)) / 2.0;
  }
  TacsScalar da = dmax_dim / 2.0;
  TacsScalar db = dmin_dim / 2.0;
  TacsScalar da4 = 4.0 * a * a * a * da;
  TacsScalar db3 = 3.0 * b * b * db;
  TacsScalar db4 = 4.0 * b3 * db;
  TacsScalar dJ =
      da * b3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0)) +
      a * b3 * (3.36 * b / a * da / a * (1.0 - b4 / a4 / 12.0)) +
      a * b3 * (3.36 * b / a * (-b4 / a4 / 12.0 * da4 / a4)) +
      a * db3 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - b4 / a4 / 12.0)) +
      a * b3 * (-3.36 * db / a * (1.0 - b4 / a4 / 12.0)) +
      a * b3 * (-3.36 * b / a * (-db4 / a4 / 12.0));
  return dJ;
}

TacsScalar TACSIsoRectangleBeamConstitutive::evalDesignFieldValue(
    int elemIndex, const double pt[], const TacsScalar X[], int index) {
  if (index == 0) {
    return width;
  } else if (index == 1) {
    return thickness;
  } else if (index == 2) {
    return buckle_length;
  }
  return 0.0;
}

const char *TACSIsoRectangleBeamConstitutive::constName =
    "TACSIsoRectangleBeamConstitutive";

/*
  Return the constitutive name
*/
const char *TACSIsoRectangleBeamConstitutive::getObjectName() {
  return constName;
}
