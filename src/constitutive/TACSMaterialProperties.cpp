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

#include "TACSMaterialProperties.h"

const static double LINEAR_STRESS_CUTOFF = 1e-15;
const static double HUGE_FAILURE_LOAD = 1e20;

TACSMaterialProperties::TACSMaterialProperties(
    TacsScalar _rho, TacsScalar _specific_heat, TacsScalar _E, TacsScalar _nu,
    TacsScalar _ys, TacsScalar _alpha, TacsScalar _kappa) {
  rho = _rho;
  specific_heat = _specific_heat;
  E = _E;
  nu = _nu;
  G = 0.5 * E / (1.0 + nu);
  ys = _ys;
  alpha = _alpha;
  kappa = _kappa;

  // Set the orthotropic material properties from the
  // isotropic properties.
  E1 = E;
  E2 = E;
  E3 = E;
  G12 = G;
  G13 = G;
  G23 = G;
  nu12 = nu;
  nu13 = nu;
  nu23 = nu;

  // Set the failure properties
  T1 = C1 = ys;
  T2 = C2 = ys;
  T3 = C3 = ys;
  S12 = S13 = S23 = sqrt(3.0) * ys;

  // Set the coefficients of thermal expansion
  alpha1 = alpha;
  alpha2 = alpha;
  alpha3 = alpha;

  // Set the coefficients of
  kappa1 = kappa;

  mat_type = TACS_ISOTROPIC_MATERIAL;
}

TACSMaterialProperties::TACSMaterialProperties(
    TacsScalar _rho, TacsScalar _specific_heat, TacsScalar _E1, TacsScalar _E2,
    TacsScalar _E3, TacsScalar _nu12, TacsScalar _nu13, TacsScalar _nu23,
    TacsScalar _G12, TacsScalar _G13, TacsScalar _G23, TacsScalar _T1,
    TacsScalar _C1, TacsScalar _T2, TacsScalar _C2, TacsScalar _T3,
    TacsScalar _C3, TacsScalar _S12, TacsScalar _S13, TacsScalar _S23,
    TacsScalar _alpha1, TacsScalar _alpha2, TacsScalar _alpha3,
    TacsScalar _kappa1, TacsScalar _kappa2, TacsScalar _kappa3) {
  mat_type = TACS_ANISOTROPIC_MATERIAL;

  rho = _rho;
  specific_heat = _specific_heat;

  // Set the isotropic properties to 0
  E = nu = G = ys = alpha = 0.0;

  // Record the orthotropic properties
  E1 = _E1;
  E2 = _E2;
  E3 = _E3;
  G12 = _G12;
  G13 = _G13;
  G23 = _G23;
  nu12 = _nu12;
  nu13 = _nu13;
  nu23 = _nu23;

  T1 = _T1;
  C1 = _C1;
  T2 = _T2;
  C2 = _C2;
  T3 = _T3;
  C3 = _C3;
  S12 = _S12;
  S13 = _S13;
  S23 = _S23;

  // Ensure that the compressive strength values are negative
  if (TacsRealPart(C1) < 0.0) {
    C1 *= -1.0;
  }
  if (TacsRealPart(C2) < 0.0) {
    C2 *= -1.0;
  }
  if (TacsRealPart(C3) < 0.0) {
    C3 *= -1.0;
  }

  alpha1 = _alpha1;
  alpha2 = _alpha2;
  alpha3 = _alpha3;
}

// Get the material type
MaterialType TACSMaterialProperties::getMaterialType() { return mat_type; }

// Extract the density property values
TacsScalar TACSMaterialProperties::getDensity() { return rho; }

TacsScalar TACSMaterialProperties::getSpecificHeat() { return specific_heat; }

// Set the density property values
void TACSMaterialProperties::setDensity(TacsScalar _rho) { rho = _rho; }

void TACSMaterialProperties::setSpecificHeat(TacsScalar _specific_heat) {
  specific_heat = _specific_heat;
}

// Extract the coefficients
void TACSMaterialProperties::getIsotropicProperties(TacsScalar *_E,
                                                    TacsScalar *_nu) {
  if (_E) {
    *_E = E;
  }
  if (_nu) {
    *_nu = nu;
  }
}

void TACSMaterialProperties::getOrthotropicProperties(
    TacsScalar *_E1, TacsScalar *_E2, TacsScalar *_E3, TacsScalar *_nu12,
    TacsScalar *_nu13, TacsScalar *_nu23, TacsScalar *_G12, TacsScalar *_G13,
    TacsScalar *_G23) {
  if (_E1) {
    *_E1 = E1;
  }
  if (_E2) {
    *_E2 = E2;
  }
  if (_E3) {
    *_E3 = E3;
  }
  if (_nu12) {
    *_nu12 = nu12;
  }
  if (_nu13) {
    *_nu13 = nu13;
  }
  if (_nu23) {
    *_nu23 = nu23;
  }
  if (_G12) {
    *_G12 = G12;
  }
  if (_G13) {
    *_G13 = G13;
  }
  if (_G23) {
    *_G23 = G23;
  }
}

void TACSMaterialProperties::getStrengthProperties(
    TacsScalar *_T1, TacsScalar *_C1, TacsScalar *_T2, TacsScalar *_C2,
    TacsScalar *_T3, TacsScalar *_C3, TacsScalar *_S12, TacsScalar *_S13,
    TacsScalar *_S23) {
  if (_T1) {
    *_T1 = T1;
  }
  if (_C1) {
    *_C1 = C1;
  }
  if (_T2) {
    *_T2 = T2;
  }
  if (_C2) {
    *_C2 = C2;
  }
  if (_T3) {
    *_T3 = T3;
  }
  if (_C3) {
    *_C3 = C3;
  }
  if (_S12) {
    *_S12 = S12;
  }
  if (_S13) {
    *_S13 = S13;
  }
  if (_S23) {
    *_S23 = S23;
  }
}

void TACSMaterialProperties::getCoefThermalExpansion(TacsScalar *_a1,
                                                     TacsScalar *_a2,
                                                     TacsScalar *_a3) {
  if (_a1) {
    *_a1 = alpha1;
  }
  if (_a2) {
    *_a2 = alpha2;
  }
  if (_a3) {
    *_a3 = alpha3;
  }
}

void TACSMaterialProperties::getThermalConductivity(TacsScalar *_k1,
                                                    TacsScalar *_k2,
                                                    TacsScalar *_k3) {
  if (_k1) {
    *_k1 = kappa1;
  }
  if (_k2) {
    *_k2 = kappa2;
  }
  if (_k3) {
    *_k3 = kappa3;
  }
}

void TACSMaterialProperties::evalTangentStiffness3D(TacsScalar C[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    TacsScalar D = E / ((1.0 + nu) * (1.0 - 2.0 * nu));

    C[0] = (1.0 - nu) * D;
    C[1] = nu * D;
    C[2] = nu * D;
    C[3] = C[4] = C[5] = 0.0;

    C[6] = (1.0 - nu) * D;
    C[7] = nu * D;
    C[8] = C[9] = C[10] = 0.0;

    C[11] = (1.0 - nu) * D;
    C[12] = C[13] = C[14] = 0.0;

    C[15] = G;
    C[16] = C[17] = 0.0;

    C[18] = G;
    C[19] = 0.0;

    C[20] = G;
  } else {
    TacsScalar nu21 = (E2 / E1) * nu12;
    TacsScalar nu31 = (E3 / E1) * nu13;
    TacsScalar nu32 = (E3 / E2) * nu23;

    TacsScalar D = 1.0 / (1.0 - nu12 * nu21 - nu13 * nu31 - nu23 * nu32 -
                          2.0 * nu21 * nu32 * nu13);

    C[0] = (1.0 - nu23 * nu32) * E1 * D;
    C[1] = (nu21 + nu31 * nu23) * E1 * D;
    C[2] = (nu31 + nu21 * nu32) * E1 * D;
    C[3] = C[4] = C[5] = 0.0;

    C[6] = (1.0 - nu13 * nu31) * E2 * D;
    C[7] = (nu32 + nu12 * nu31) * E2 * D;
    C[8] = C[9] = C[10] = 0.0;

    C[11] = (1.0 - nu12 * nu21) * E3 * D;
    C[12] = C[13] = C[14] = 0.0;

    C[15] = G23;
    C[16] = C[17] = 0.0;

    C[18] = G13;
    C[19] = 0.0;

    C[20] = G12;
  }
}

void TACSMaterialProperties::evalTangentStiffness2D(TacsScalar C[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    TacsScalar D = E / (1.0 - nu * nu);
    C[0] = D;
    C[1] = nu * D;
    C[2] = 0.0;
    C[3] = D;
    C[4] = 0.0;
    C[5] = G;
  } else {
    TacsScalar nu21 = nu12 * E2 / E1;
    C[0] = E1 / (1.0 - nu12 * nu21);
    C[1] = nu12 * E2 / (1.0 - nu12 * nu21);
    C[2] = 0.0;
    C[3] = E2 / (1.0 - nu12 * nu21);
    C[4] = 0.0;
    C[5] = G12;
  }
}

void TACSMaterialProperties::evalTangentHeatFlux3D(TacsScalar C[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    C[0] = C[3] = C[5] = kappa;
    C[1] = C[2] = C[4] = 0.0;
  } else {
    C[0] = kappa1;
    C[3] = kappa2;
    C[5] = kappa2;
    C[1] = C[2] = C[4] = 0.0;
  }
}

void TACSMaterialProperties::evalTangentHeatFlux2D(TacsScalar C[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    C[0] = C[2] = kappa;
    C[1] = 0.0;
  } else {
    C[0] = kappa1;
    C[2] = kappa2;
    C[1] = 0.0;
  }
}

void TACSMaterialProperties::evalThermalStrain3D(TacsScalar e[]) {
  e[0] = alpha1;
  e[1] = alpha2;
  e[2] = alpha3;
  e[3] = e[4] = e[5] = 0.0;
}

void TACSMaterialProperties::evalThermalStrain2D(TacsScalar e[]) {
  e[0] = alpha1;
  e[1] = alpha2;
  e[2] = 0.0;
}

void TACSMaterialProperties::evalStress2D(const TacsScalar e[],
                                          TacsScalar s[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    TacsScalar D = E / (1.0 - nu * nu);

    s[0] = D * (e[0] + nu * e[1]);
    s[1] = D * (e[1] + nu * e[0]);
    s[2] = G * e[2];
  } else {
    TacsScalar C[6];
    evalTangentStiffness2D(C);
    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2];
    s[1] = C[1] * e[0] + C[3] * e[1] + C[4] * e[2];
    s[2] = C[2] * e[0] + C[4] * e[1] + C[5] * e[2];
  }
}

void TACSMaterialProperties::evalStress3D(const TacsScalar e[],
                                          TacsScalar s[]) {
  if (mat_type == TACS_ISOTROPIC_MATERIAL) {
    TacsScalar D = E / ((1.0 + nu) * (1.0 - 2.0 * nu));

    s[0] = D * ((1.0 - nu) * e[0] + nu * e[1] + nu * e[2]);
    s[1] = D * ((1.0 - nu) * e[1] + nu * e[0] + nu * e[2]);
    s[2] = D * ((1.0 - nu) * e[2] + nu * e[0] + nu * e[1]);
    s[3] = G * e[3];
    s[4] = G * e[4];
    s[5] = G * e[5];
  } else {
    TacsScalar C[21];
    evalTangentStiffness3D(C);

    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2] + C[3] * e[3] + C[4] * e[4] +
           C[5] * e[5];
    s[1] = C[1] * e[0] + C[6] * e[1] + C[7] * e[2] + C[8] * e[3] + C[9] * e[4] +
           C[10] * e[5];
    s[2] = C[2] * e[0] + C[7] * e[1] + C[11] * e[2] + C[12] * e[3] +
           C[13] * e[4] + C[14] * e[5];
    s[3] = C[3] * e[0] + C[8] * e[1] + C[12] * e[2] + C[15] * e[3] +
           C[16] * e[4] + C[17] * e[5];
    s[4] = C[4] * e[0] + C[9] * e[1] + C[13] * e[2] + C[16] * e[3] +
           C[18] * e[4] + C[19] * e[5];
    s[5] = C[5] * e[0] + C[10] * e[1] + C[14] * e[2] + C[17] * e[3] +
           C[19] * e[4] + C[20] * e[5];
  }
}

/*!
  The von Mises failure criteria is the following:

  (sx - sy)^2 + (sx - sz)^2 + (sy - sz)^2 +
  6.0 * ( sxy^2 + sxz^2 + syz^2 ) = 2.0 * ys^2

  sx  = s[0]
  sy  = s[1]
  sz  = s[2]
  syz = s[3]
  sxz = s[4]
  sxy = s[5]
*/
TacsScalar TACSMaterialProperties::vonMisesFailure3D(const TacsScalar s[]) {
  TacsScalar fail =
      sqrt(0.5 *
           ((s[0] - s[1]) * (s[0] - s[1]) + (s[0] - s[2]) * (s[0] - s[2]) +
            (s[1] - s[2]) * (s[1] - s[2]) +
            6.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))) /
      ys;
  return fail;
}

TacsScalar TACSMaterialProperties::vonMisesFailure3DStressSens(
    const TacsScalar s[], TacsScalar sens[]) {
  TacsScalar fail = sqrt(
      0.5 * ((s[0] - s[1]) * (s[0] - s[1]) + (s[0] - s[2]) * (s[0] - s[2]) +
             (s[1] - s[2]) * (s[1] - s[2]) +
             6.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5])));

  if (fail != 0.0) {
    TacsScalar fact = 0.5 / (ys * fail);
    sens[0] = fact * (2.0 * s[0] - s[1] - s[2]);
    sens[1] = fact * (2.0 * s[1] - s[0] - s[2]);
    sens[2] = fact * (2.0 * s[2] - s[0] - s[1]);
    sens[3] = 6.0 * fact * s[3];
    sens[4] = 6.0 * fact * s[4];
    sens[5] = 6.0 * fact * s[5];
  } else {
    sens[0] = sens[1] = sens[2] = 0.0;
    sens[3] = sens[4] = sens[5] = 0.0;
  }

  return fail;
}

/*
  Compute the von Mises failure criteria:

  (sx^2 + sy^2 - sx*sy + 3.0*sxy^2)/ys^2 < 1.0

  sx  = s[0]
  sy  = s[1]
  sxy = s[2]
*/
TacsScalar TACSMaterialProperties::vonMisesFailure2D(const TacsScalar s[]) {
  TacsScalar fail =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]) / ys;
  return fail;
}

TacsScalar TACSMaterialProperties::vonMisesFailure2DStressSens(
    const TacsScalar s[], TacsScalar sens[]) {
  TacsScalar fail =
      sqrt(s[0] * s[0] + s[1] * s[1] - s[0] * s[1] + 3.0 * s[2] * s[2]);

  if (fail != 0.0) {
    sens[0] = (s[0] - 0.5 * s[1]) / (fail * ys);
    sens[1] = (s[1] - 0.5 * s[0]) / (fail * ys);
    sens[2] = (3.0 * s[2]) / (fail * ys);
  } else {
    sens[0] = sens[1] = sens[2] = 0.0;
  }
  fail = fail / ys;

  return fail;
}

/*
  The OrthoPly material class.

  This uses a Tsai-Wu tensor failure criteria.

  Inputs are:
  E1, E2, nu12, G12 = In-plane stiffness properties
  G23, G13 = Out of plane shear stiffness

  Xt, Xc = Tensile and compressive fiber-aligned failure loads
  Yt, Yc = Tensile and compressive off-axis failure loads
  S12 = In plane shear failure load
  C = Interaction strength such that sigma_1 = sigma_2 = C
*/
TACSOrthotropicPly::TACSOrthotropicPly(TacsScalar _plyThickness,
                                       TACSMaterialProperties *_properties) {
  plyThickness = _plyThickness;
  properties = _properties;
  properties->incref();

  rho = properties->getDensity();

  properties->getOrthotropicProperties(&E1, &E2, NULL, &nu12, NULL, NULL, &G12,
                                       &G13, &G23);
  nu21 = nu12 * E2 / E1;
  Q11 = E1 / (1.0 - nu12 * nu21);
  Q22 = E2 / (1.0 - nu12 * nu21);
  Q12 = nu12 * E2 / (1.0 - nu12 * nu21);

  Q44 = G23;
  Q55 = G13;
  Q66 = G12;

  // Definitions from Jones, Mechanics of Composite Materials
  C12 = (Q11 + Q22 - 4.0 * Q66);
  C16 = (Q11 - Q12 - 2.0 * Q66);
  C26 = (Q12 - Q22 + 2.0 * Q66);
  C66 = (Q11 + Q22 - 2.0 * Q12 - 2.0 * Q66);

  // Record the failure data. If the coefficients
  // are negative, make them positive
  properties->getStrengthProperties(&Xt, &Xc, &Yt, &Yc, NULL, NULL, &S12, NULL,
                                    NULL);
  C = 0.0;

  // Determine the strain strength values based on the
  // supplied stress strength values.
  eXt = Xt / E1;
  eXc = Xc / E1;
  eYt = Yt / E2;
  eYc = Yc / E2;
  eS12 = S12 / G12;

  // By default, use Tsai-Wu
  useTsaiWuCriterion = 1;

  // Compute the coefficients for the Tsai-Wu failure criteria
  F1 = (Xc - Xt) / (Xt * Xc);
  F2 = (Yc - Yt) / (Yt * Yc);
  F11 = 1.0 / (Xt * Xc);
  F22 = 1.0 / (Yt * Yc);
  F66 = 1.0 / (S12 * S12);
  if (TacsRealPart(C) != 0.0) {
    F12 = 0.5 * (1.0 - (F1 + F2) * C - (F11 + F22) * C * C) / (C * C);
  } else {
    F12 = 0.0;  // Assume no interaction
  }

  // Check the stability criterion
  if (TacsRealPart(F12 * F12) >= TacsRealPart(F11 * F22)) {
    fprintf(stderr,
            "TACSOrthotropicPly: Value of C = %e results in "
            "non-physical F12 = %e. Setting F12 = 0.\n",
            TacsRealPart(C), TacsRealPart(F12));
    fprintf(stderr,
            "TACSOrthotropicPly: Tsai-Wu coefficients: F11: "
            "%e, F22: %e, F66: %e, F1: %e, F2: %e\n",
            TacsRealPart(F11), TacsRealPart(F22), TacsRealPart(F66),
            TacsRealPart(F1), TacsRealPart(F2));
    F12 = 0.0;
  }

  // Set the default value of the KS penalty
  ksWeight = 100.0;
}

const char *TACSOrthotropicPly::name = "TACSOrthotropicPly";

/*
  Set the KS penalty factor
*/
void TACSOrthotropicPly::setKSWeight(TacsScalar _ksWeight) {
  ksWeight = _ksWeight;
}

/*
  Set the failure criteria to use
*/
void TACSOrthotropicPly::setUseMaxStrainCriterion() { useTsaiWuCriterion = 0; }

void TACSOrthotropicPly::setUseTsaiWuCriterion() { useTsaiWuCriterion = 1; }

/*
  Get the density of the material
*/
TacsScalar TACSOrthotropicPly::getDensity() { return rho; }

/*
  Get the thickness of a single laminae
*/
TacsScalar TACSOrthotropicPly::getPlyThickness() { return plyThickness; }

/*
  Get the stiffness constants
*/
void TACSOrthotropicPly::getStiffness(TacsScalar *_E1, TacsScalar *_E2,
                                      TacsScalar *_nu12, TacsScalar *_G12,
                                      TacsScalar *_G23, TacsScalar *_G13) {
  *_E1 = E1;
  *_E2 = E2;
  *_nu12 = nu12;
  *_G12 = G12;
  *_G23 = G23;
  *_G13 = G13;
}

/*
  Get the lamination stiffness objects
*/
void TACSOrthotropicPly::getLaminateStiffness(
    TacsScalar *_Q11, TacsScalar *_Q12, TacsScalar *_Q22, TacsScalar *_Q44,
    TacsScalar *_Q55, TacsScalar *_Q66) {
  *_Q11 = Q11;
  *_Q12 = Q12;
  *_Q22 = Q22;
  *_Q44 = Q44;
  *_Q55 = Q55;
  *_Q66 = Q66;
}

/*
  Get the strength parameters
*/
void TACSOrthotropicPly::getStrength(TacsScalar *_Xt, TacsScalar *_Xc,
                                     TacsScalar *_Yt, TacsScalar *_Yc,
                                     TacsScalar *_S12, TacsScalar *_C) {
  *_Xt = Xt;
  *_Xc = Xc;
  *_Yt = Yt;
  *_Yc = Yc;
  *_S12 = S12;
  *_C = C;
}

/*
  Get the strength parameters for strain
*/
void TACSOrthotropicPly::getStrainStrength(TacsScalar *_eXt, TacsScalar *_eXc,
                                           TacsScalar *_eYt, TacsScalar *_eYc,
                                           TacsScalar *_eS12) {
  *_eXt = eXt;
  *_eXc = eXc;
  *_eYt = eYt;
  *_eYc = eYc;
  *_eS12 = eS12;
}

/*
  Get the Tsai-Wu failure constants
*/
void TACSOrthotropicPly::getTsaiWu(TacsScalar *_F1, TacsScalar *_F2,
                                   TacsScalar *_F11, TacsScalar *_F12,
                                   TacsScalar *_F22, TacsScalar *_F66) {
  *_F1 = F1;
  *_F2 = F2;
  *_F11 = F11;
  *_F12 = F12;
  *_F22 = F22;
  *_F66 = F66;
}

/*
  Retrieve the stiffness invariants for the laminate
*/
void TACSOrthotropicPly::getLaminateInvariants(TacsScalar *U1, TacsScalar *U2,
                                               TacsScalar *U3, TacsScalar *U4,
                                               TacsScalar *U5, TacsScalar *U6) {
  *U1 = 0.125 * (3.0 * Q11 + 3.0 * Q22 + 2.0 * Q12 + 4.0 * Q66);
  *U2 = 0.5 * (Q11 - Q22);
  *U3 = 0.125 * (Q11 + Q22 - 2.0 * Q12 - 4.0 * Q66);
  *U4 = 0.125 * (Q11 + Q22 + 6.0 * Q12 - 4.0 * Q66);
  *U5 = 0.5 * (Q44 + Q55);
  *U6 = 0.5 * (Q44 - Q55);
}

/*!
  Calculate the in-plane stresses in the ply-coordinate axis
*/
void TACSOrthotropicPly::getPlyStress(const TacsScalar strain[],
                                      TacsScalar stress[]) {
  stress[0] = Q11 * strain[0] + Q12 * strain[1];
  stress[1] = Q12 * strain[0] + Q22 * strain[1];
  stress[2] = Q66 * strain[2];
}

/*
  Compute the rotation
*/
void TACSOrthotropicPly::calculateAbar(TacsScalar angle, TacsScalar Abar[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  Abar[0] = cos2 * Q44 + sin2 * Q55;
  Abar[1] = cos1 * sin1 * (Q55 - Q44);
  Abar[2] = sin2 * Q44 + cos2 * Q55;
}

/*
  Calculate the derivative of the Abar matrix w.r.t. the angle
*/
void TACSOrthotropicPly::calculateAbarAngleSens(TacsScalar angle,
                                                TacsScalar Abar[]) {
  TacsScalar cos_2 = cos(2.0 * angle);
  TacsScalar sin_2 = sin(2.0 * angle);

  Abar[0] = -sin_2 * Q44 + sin_2 * Q55;
  Abar[1] = cos_2 * (Q55 - Q44);
  Abar[2] = sin_2 * Q44 - sin_2 * Q55;
}

/*
  Taken from Jones, Mechanics of composite materials pg. 51
*/
void TACSOrthotropicPly::calculateQbar(TacsScalar angle, TacsScalar Qbar[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  TacsScalar cos4 = cos2 * cos2;
  TacsScalar sin4 = sin2 * sin2;

  // First row of the matrix
  Qbar[0] = Q11 * cos4 + 2.0 * (Q12 + 2.0 * Q66) * sin2 * cos2 + Q22 * sin4;
  Qbar[1] = C12 * sin2 * cos2 + Q12 * (sin4 + cos4);
  Qbar[2] = C16 * sin1 * cos2 * cos1 + C26 * sin2 * sin1 * cos1;

  // Second row of the matrix
  Qbar[3] = Q11 * sin4 + 2.0 * (Q12 + 2.0 * Q66) * sin2 * cos2 + Q22 * cos4;
  Qbar[4] = C16 * sin2 * sin1 * cos1 + C26 * sin1 * cos2 * cos1;

  // Last row of the matrix
  Qbar[5] = C66 * sin2 * cos2 + Q66 * (sin4 + cos4);
}

/*
  Calculate the sensitivity of the Qbar matrix w.r.t. the angle
*/
void TACSOrthotropicPly::calculateQbarAngleSens(TacsScalar angle,
                                                TacsScalar Qbar[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  TacsScalar s_cos2 = -2.0 * cos1 * sin1;
  TacsScalar s_sin2 = 2.0 * sin1 * cos1;

  TacsScalar s_cos4 = -4.0 * cos2 * cos1 * sin1;
  TacsScalar s_sin4 = 4.0 * sin2 * sin1 * cos1;

  TacsScalar s_sin2cos2 = (sin2 * s_cos2 + s_sin2 * cos2);

  TacsScalar s_sin3cos1 = 3.0 * sin2 * cos2 - sin2 * sin2;
  TacsScalar s_cos3sin1 = -3.0 * cos2 * sin2 + cos2 * cos2;

  // First row of the matrix
  Qbar[0] = Q11 * s_cos4 + 2.0 * (Q12 + 2.0 * Q66) * s_sin2cos2 + Q22 * s_sin4;
  Qbar[1] = C12 * s_sin2cos2 + Q12 * (s_sin4 + s_cos4);
  Qbar[2] = C16 * s_cos3sin1 + C26 * s_sin3cos1;

  // Second row of the matrix
  Qbar[3] = Q11 * s_sin4 + 2.0 * (Q12 + 2.0 * Q66) * s_sin2cos2 + Q22 * s_cos4;
  Qbar[4] = C16 * s_sin3cos1 + C26 * s_cos3sin1;

  // Last row of the matrix
  Qbar[5] = C66 * s_sin2cos2 + Q66 * (s_sin4 + s_cos4);
}

void TACSOrthotropicPly::calculateStress(TacsScalar angle,
                                         const TacsScalar strain[],
                                         TacsScalar stress[]) {
  TacsScalar strn[3], strs[3];
  transformStrainGlobal2Ply(angle, strain, strn);
  getPlyStress(strn, strs);
  transformStressPly2Global(angle, strs, stress);
}

/*
  Calculate the failure criteria

  fail (the return value) is such that

  fail <  1.0 ==> Okay
  fail >= 1.0 ==> Material failure

  One of two failure criteria are used:

  1. The Tsai-Wu failure criterion:

  F(s) =
  F1*s[0] + F2*s[1]
  F11*s[0]**2 + F22*s[1]**2 +
  2.0*F12*s[0]*s[1] + F66*s[2]**2 <= 1.0

  2. The maximum strain failure criteria:

  KS(e_{i}/e_max^{+/-}, ksWeight) <= 1.0

  where e_max^{+/-} are the postive and negative strains at failure
  for component i. This ensures that the ratio of the maximum strain
  over the strain at failure is below 1.
*/
TacsScalar TACSOrthotropicPly::failure(TacsScalar angle,
                                       const TacsScalar strain[]) {
  TacsScalar e[3];  // Ply strain
  transformStrainGlobal2Ply(angle, strain, e);

  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion) {
    TacsScalar s[3];  // Ply stress
    getPlyStress(e, s);

    fail = (F11 * s[0] * s[0] + F22 * s[1] * s[1] + 2.0 * F12 * s[0] * s[1] +
            F66 * s[2] * s[2] + F1 * s[0] + F2 * s[1]);
  } else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6];
    f[0] = e[0] / eXt;
    f[1] = -e[0] / eXc;
    f[2] = e[1] / eYt;
    f[3] = -e[1] / eYc;
    f[4] = e[2] / eS12;
    f[5] = -e[2] / eS12;

    TacsScalar max = f[0];
    for (int k = 1; k < 6; k++) {
      if (TacsRealPart(f[k]) > TacsRealPart(max)) {
        max = f[k];
      }
    }

    TacsScalar ksSum = 0.0;
    for (int k = 0; k < 6; k++) {
      ksSum += exp(ksWeight * (f[k] - max));
    }
    fail = max + log(ksSum) / ksWeight;
  }

  return fail;
}

TacsScalar TACSOrthotropicPly::failureStrainSens(TacsScalar angle,
                                                 const TacsScalar strain[],
                                                 TacsScalar sens[]) {
  TacsScalar e[3];  // Ply strain
  transformStrainGlobal2Ply(angle, strain, e);

  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion) {
    TacsScalar s[3];  // Ply stress
    getPlyStress(e, s);

    fail = (F11 * s[0] * s[0] + F22 * s[1] * s[1] + 2.0 * F12 * s[0] * s[1] +
            F66 * s[2] * s[2] + F1 * s[0] + F2 * s[1]);

    sens[0] = F1 + 2.0 * F11 * s[0] + 2.0 * F12 * s[1];
    sens[1] = F2 + 2.0 * F22 * s[1] + 2.0 * F12 * s[0];
    sens[2] = 2.0 * F66 * s[2];

    TacsScalar sSens[3];
    getPlyStress(sens, sSens);

    transformStressPly2Global(angle, sSens, sens);
  } else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6], fexp[6];
    f[0] = e[0] / eXt;
    f[1] = -e[0] / eXc;
    f[2] = e[1] / eYt;
    f[3] = -e[1] / eYc;
    f[4] = e[2] / eS12;
    f[5] = -e[2] / eS12;

    TacsScalar max = f[0];
    for (int k = 1; k < 6; k++) {
      if (TacsRealPart(f[k]) > TacsRealPart(max)) {
        max = f[k];
      }
    }

    TacsScalar ksSum = 0.0;
    for (int k = 0; k < 6; k++) {
      fexp[k] = exp(ksWeight * (f[k] - max));
      ksSum += fexp[k];
    }

    fail = max + log(ksSum) / ksWeight;

    TacsScalar sSens[3];
    sSens[0] = (fexp[0] * eXc - fexp[1] * eXt) / (ksSum * eXt * eXc);
    sSens[1] = (fexp[2] * eYc - fexp[3] * eYt) / (ksSum * eYt * eYc);
    sSens[2] = (fexp[4] - fexp[5]) / (ksSum * eS12);

    transformStressPly2Global(angle, sSens, sens);
  }

  return fail;
}

TacsScalar TACSOrthotropicPly::failureAngleSens(TacsScalar angle,
                                                const TacsScalar strain[],
                                                TacsScalar *failSens) {
  TacsScalar e[3], se[3];  // The ply strain
  transformStrainGlobal2Ply(angle, strain, e);

  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion) {
    TacsScalar s[3];  // The ply stress
    getPlyStress(e, s);

    fail = (F11 * s[0] * s[0] + F22 * s[1] * s[1] + 2.0 * F12 * s[0] * s[1] +
            F66 * s[2] * s[2] + F1 * s[0] + F2 * s[1]);

    // Compute the sensitivity of the transformation
    transformStrainGlobal2PlyAngleSens(angle, strain, se);

    // Compute the sensitivity of the stress
    TacsScalar ss[3];
    getPlyStress(se, ss);

    *failSens =
        (2.0 * (F11 * s[0] * ss[0] + F22 * s[1] * ss[1] +
                F12 * (ss[0] * s[1] + s[0] * ss[1]) + F66 * s[2] * ss[2]) +
         F1 * ss[0] + F2 * ss[1]);
  } else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6], fs[6];
    f[0] = e[0] / eXt;
    f[1] = -e[0] / eXc;
    f[2] = e[1] / eYt;
    f[3] = -e[1] / eYc;
    f[4] = e[2] / eS12;
    f[5] = -e[2] / eS12;

    // Compute the sensitivity of the transformation
    transformStrainGlobal2PlyAngleSens(angle, strain, se);
    fs[0] = se[0] / eXt;
    fs[1] = -se[0] / eXc;
    fs[2] = se[1] / eYt;
    fs[3] = -se[1] / eYc;
    fs[4] = se[2] / eS12;
    fs[5] = -se[2] / eS12;

    TacsScalar max = f[0];
    for (int k = 1; k < 6; k++) {
      if (TacsRealPart(f[k]) > TacsRealPart(max)) {
        max = f[k];
      }
    }

    TacsScalar ksSum = 0.0;
    TacsScalar fSens = 0.0;
    for (int k = 0; k < 6; k++) {
      TacsScalar fexp = exp(ksWeight * (f[k] - max));
      fSens += fexp * fs[k];
      ksSum += fexp;
    }

    fail = max + log(ksSum) / ksWeight;

    *failSens = fSens / ksSum;
  }

  return fail;
}

/*
  Calculate the failure load based on constant and linear components
  of the strain. The components of the strain are in the laminate
  coordinates (global coordinates) but just the in-plane
  contributions. The angle specified in radians is angle that the ply
  is oriented at - to get to the lamina coordinate system. See Jones -
  Mechanics of Composite Materials for details.

  returns:
  the multiple of lstrain at which failure occurs

  input:
  angle = [radians] the orientation angle of the lamina
  cstrain = the constant in-plane strain components
  lstrain = the linear in-plane strain components

  This works by computing F(cstrain, P*lstrain) = 0 for some P. This
  turns out to be a quadratic in P.
*/
TacsScalar TACSOrthotropicPly::calculateFailLoad(TacsScalar angle,
                                                 const TacsScalar cstrain[],
                                                 const TacsScalar lstrain[]) {
  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3];
  transformStrainGlobal2Ply(angle, cstrain, cstn);
  transformStrainGlobal2Ply(angle, lstrain, lstn);

  TacsScalar cstr[3], lstr[3];
  getPlyStress(cstn, cstr);
  getPlyStress(lstn, lstr);

  TacsScalar c = (F1 * cstr[0] + F2 * cstr[1] + F11 * cstr[0] * cstr[0] +
                  F22 * cstr[1] * cstr[1] + 2.0 * F12 * cstr[0] * cstr[1] +
                  F66 * cstr[2] * cstr[2]) -
                 1.0;
  TacsScalar b = (F1 * lstr[0] + F2 * lstr[1] + 2.0 * F11 * cstr[0] * lstr[0] +
                  2.0 * F22 * cstr[1] * lstr[1] +
                  2.0 * F12 * (cstr[0] * lstr[1] + cstr[1] * lstr[0]) +
                  2.0 * F66 * cstr[2] * lstr[2]);
  TacsScalar a = (F11 * lstr[0] * lstr[0] + F22 * lstr[1] * lstr[1] +
                  2.0 * F12 * lstr[0] * lstr[1] + F66 * lstr[2] * lstr[2]);
  TacsScalar pos = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF * TacsRealPart(F11 + F22)) {
    pos = HUGE_FAILURE_LOAD;
  } else if (TacsRealPart(c) >= 0.0) {
    pos = 0.0;
  } else {
    TacsScalar discrim = b * b - 4.0 * a * c;
    if (TacsRealPart(discrim) < 0.0) {
      pos = 0.0;
    } else {
      discrim = sqrt(discrim);

      if (TacsRealPart(b) >= 0.0) {
        pos = 2.0 * c / (b + discrim);
      } else {  // b < 0.0
        pos = 0.5 * (-b + discrim) / a;
      }
    }
  }

  return pos;
}

/*
  Determine the sensitivity of the failure load calculation above to
  the strain components cstrain and lstrain. Here the arguments angle
  [radians], and the in-plane strain components cstrain[] and
  lstrain[] are the same as above. The output cSens[], lSens[] are
  the sensitivity of the failure load w.r.t. to the constant and linear
  strain components.
*/
TacsScalar TACSOrthotropicPly::calculateFailLoadStrainSens(
    TacsScalar angle, const TacsScalar cstrain[], const TacsScalar lstrain[],
    TacsScalar cSens[], TacsScalar lSens[]) {
  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3];
  transformStrainGlobal2Ply(angle, cstrain, cstn);
  transformStrainGlobal2Ply(angle, lstrain, lstn);

  // The constant and linearly varying components of stress - in the ply axis
  TacsScalar cstr[3], lstr[3];
  getPlyStress(cstn, cstr);
  getPlyStress(lstn, lstr);

  TacsScalar c = (F1 * cstr[0] + F2 * cstr[1] + F11 * cstr[0] * cstr[0] +
                  F22 * cstr[1] * cstr[1] + 2.0 * F12 * cstr[0] * cstr[1] +
                  F66 * cstr[2] * cstr[2]) -
                 1.0;
  TacsScalar b = (F1 * lstr[0] + F2 * lstr[1] + 2.0 * F11 * cstr[0] * lstr[0] +
                  2.0 * F22 * cstr[1] * lstr[1] +
                  2.0 * F12 * (cstr[0] * lstr[1] + cstr[1] * lstr[0]) +
                  2.0 * F66 * cstr[2] * lstr[2]);
  TacsScalar a = (F11 * lstr[0] * lstr[0] + F22 * lstr[1] * lstr[1] +
                  2.0 * F12 * lstr[0] * lstr[1] + F66 * lstr[2] * lstr[2]);

  TacsScalar pos = 0.0;
  TacsScalar pa = 0.0, pb = 0.0, pc = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF * TacsRealPart(F11 + F22)) {
    pos = HUGE_FAILURE_LOAD;
  } else if (TacsRealPart(c) >= 0.0) {
    pos = 0.0;
  } else {
    TacsScalar discrim = b * b - 4.0 * a * c;
    if (TacsRealPart(discrim) < 0.0) {
      pos = 0.0;
    } else {
      discrim = sqrt(discrim);

      if (TacsRealPart(b) >= 0.0) {
        pos = 2.0 * c / (b + discrim);
      } else {
        pos = 0.5 * (-b + discrim) / a;
      }

      pa = -(c + pos * discrim) / (discrim * a);
      pb = -pos / discrim;
      pc = -1.0 / discrim;
    }
  }

  // Now, multiply the derivative of dp/da, dp/db, dp/dc by
  // da/dcstr, db/dcstr,

  cSens[0] = (pc * (F1 + 2.0 * F11 * cstr[0] + 2.0 * F12 * cstr[1]) +
              pb * (2.0 * F11 * lstr[0] + 2.0 * F12 * lstr[1]));
  cSens[1] = (pc * (F2 + 2.0 * F12 * cstr[0] + 2.0 * F22 * cstr[1]) +
              pb * (2.0 * F12 * lstr[0] + 2.0 * F22 * lstr[1]));
  cSens[2] = (pc * (2.0 * F66 * cstr[2]) + pb * (2.0 * F66 * lstr[2]));

  lSens[0] = (pb * (F1 + 2.0 * F11 * cstr[0] + 2.0 * F12 * cstr[1]) +
              pa * (2.0 * F11 * lstr[0] + 2.0 * F12 * lstr[1]));
  lSens[1] = (pb * (F2 + 2.0 * F12 * cstr[1] + 2.0 * F22 * cstr[1]) +
              pa * (2.0 * F12 * lstr[0] + 2.0 * F22 * lstr[1]));
  lSens[2] = (pb * (2.0 * F66 * cstr[2]) + pa * (2.0 * F66 * lstr[2]));

  TacsScalar cstrSens[3], lstrSens[3];
  getPlyStress(cSens, cstrSens);
  getPlyStress(lSens, lstrSens);

  transformStressPly2Global(angle, cstrSens, cSens);
  transformStressPly2Global(angle, lstrSens, lSens);

  return pos;
}

/*
  Determine the sensitivity of the failure load calculation to to ply
  angle. Here, the arguments are the same. The return value posSens is
  the sensitivity of the failure load w.r.t. the ply angle.
*/
TacsScalar TACSOrthotropicPly::calculateFailLoadAngleSens(
    TacsScalar angle, const TacsScalar cstrain[], const TacsScalar lstrain[],
    TacsScalar *posSens) {
  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3];
  transformStrainGlobal2Ply(angle, cstrain, cstn);
  transformStrainGlobal2Ply(angle, lstrain, lstn);

  // The constant and linearly varying components of stress - in the ply axis
  TacsScalar cstr[3], lstr[3];
  getPlyStress(cstn, cstr);
  getPlyStress(lstn, lstr);

  // Now, determine the sensitivity of the transformation
  transformStrainGlobal2PlyAngleSens(angle, cstrain, cstn);
  transformStrainGlobal2PlyAngleSens(angle, lstrain, lstn);

  // The constant and linearly varying sensitivites w.r.t. the strain
  TacsScalar scstr[3], slstr[3];
  getPlyStress(cstn, scstr);
  getPlyStress(lstn, slstr);

  TacsScalar c = (F1 * cstr[0] + F2 * cstr[1] + F11 * cstr[0] * cstr[0] +
                  F22 * cstr[1] * cstr[1] + 2.0 * F12 * cstr[0] * cstr[1] +
                  F66 * cstr[2] * cstr[2]) -
                 1.0;
  TacsScalar b = (F1 * lstr[0] + F2 * lstr[1] + 2.0 * F11 * cstr[0] * lstr[0] +
                  2.0 * F22 * cstr[1] * lstr[1] +
                  2.0 * F12 * (cstr[0] * lstr[1] + cstr[1] * lstr[0]) +
                  2.0 * F66 * cstr[2] * lstr[2]);
  TacsScalar a = (F11 * lstr[0] * lstr[0] + F22 * lstr[1] * lstr[1] +
                  2.0 * F12 * lstr[0] * lstr[1] + F66 * lstr[2] * lstr[2]);

  *posSens = 0.0;
  TacsScalar pos = 0.0;
  TacsScalar pa = 0.0, pb = 0.0, pc = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF * TacsRealPart(F11 + F22)) {
    pos = HUGE_FAILURE_LOAD;
  } else if (TacsRealPart(c) >= 0.0) {
    pos = 0.0;
  } else {
    TacsScalar discrim = b * b - 4.0 * a * c;
    if (TacsRealPart(discrim) < 0.0) {
      pos = 0.0;
    } else {
      discrim = sqrt(discrim);

      if (TacsRealPart(b) >= 0.0) {
        pos = 2.0 * c / (b + discrim);
      } else {
        pos = 0.5 * (-b + discrim) / a;
      }

      pa = -(c + pos * discrim) / (discrim * a);
      pb = -pos / discrim;
      pc = -1.0 / discrim;
    }
  }

  TacsScalar sc =
      (F1 * scstr[0] + F2 * scstr[1] + 2.0 * F11 * cstr[0] * scstr[0] +
       2.0 * F22 * cstr[1] * scstr[1] +
       2.0 * F12 * (cstr[0] * scstr[1] + scstr[0] * cstr[1]) +
       2.0 * F66 * cstr[2] * scstr[2]);
  TacsScalar sb = (F1 * slstr[0] + F2 * slstr[1] +
                   2.0 * F11 * (cstr[0] * slstr[0] + scstr[0] * lstr[0]) +
                   2.0 * F22 * (cstr[1] * slstr[1] + scstr[1] * lstr[1]) +
                   2.0 * F12 *
                       (cstr[0] * slstr[1] + scstr[1] * lstr[0] +
                        scstr[0] * lstr[1] + cstr[1] * slstr[0]) +
                   2.0 * F66 * (cstr[2] * slstr[2] + scstr[2] * lstr[2]));
  TacsScalar sa =
      (2.0 * F11 * lstr[0] * slstr[0] + 2.0 * F22 * lstr[1] * slstr[1] +
       2.0 * F12 * (lstr[0] * slstr[1] + slstr[0] * lstr[1]) +
       2.0 * F66 * lstr[2] * slstr[2]);

  *posSens = pa * sa + pb * sb + pc * sc;
  return pos;
}

/*!
  The following function tests the accuracy of the implementation
  of the sensitivities for the failure calculation and the sensitivity
  of the transformation equations
*/
void TACSOrthotropicPly::testFailSens(double dh, TacsScalar angle) {
  // Select various strain components ...
  TacsScalar strain[3];

  printf("\nTesting failure sensitivity for angle = %5.2f \n",
         TacsRealPart(angle));

  for (int k = 0; k < 3; k++) {
    strain[k] = -1.0;

    // Calculate the failure load
    TacsScalar p = failure(angle, strain);
    printf("Failure criteria = %15.8e \n", TacsRealPart(p));

    // Calculate the sensitivity of the failure load
    TacsScalar sens[3];
    failureStrainSens(angle, strain, sens);

    // Compare the result to a finite-difference calculation
    for (int j = 0; j < 3; j++) {
      // Test the constant component
      TacsScalar val = strain[j];
      strain[j] = val + dh;
      TacsScalar p1 = failure(angle, strain);

      strain[j] = val - dh;
      TacsScalar p2 = failure(angle, strain);
      strain[j] = val;

      TacsScalar fd = 0.5 * (p1 - p2) / dh;
      printf("sens[%d] FD: %15.8e An: %15.8e Error: %10.3e \n", j,
             TacsRealPart(fd), TacsRealPart(sens[j]),
             TacsRealPart((fd - sens[j]) / sens[j]));
    }

    // Calculate the sensitivity w.r.t. the angle
    TacsScalar paSens;
    failureAngleSens(angle, strain, &paSens);
    TacsScalar p1 = failure(angle + dh, strain);
    TacsScalar p2 = failure(angle - dh, strain);
    TacsScalar fd = 0.5 * (p1 - p2) / dh;

    printf("Angle sensitivity FD: %15.8e An: %15.8e Error: %10.3e \n",
           TacsRealPart(fd), TacsRealPart(paSens),
           TacsRealPart((fd - paSens) / paSens));
  }
}

/*
  Print the stiffness and failure properties
*/
void TACSOrthotropicPly::printProperties() {
  printf("\nStiffness properties \n");
  printf("E1   = %15.5e \n", TacsRealPart(E1));
  printf("E2   = %15.5e \n", TacsRealPart(E2));
  printf("nu12 = %15.5e \n", TacsRealPart(nu12));
  printf("nu21 = %15.5e \n", TacsRealPart(nu21));
  printf("G12  = %15.5e \n", TacsRealPart(G12));
  printf("G23  = %15.5e \n", TacsRealPart(G23));
  printf("G13  = %15.5e \n", TacsRealPart(G13));

  printf("\nFailure Properties \n");
  printf("Xt   = %15.5e \n", TacsRealPart(Xt));
  printf("Xc   = %15.5e \n", TacsRealPart(Xc));
  printf("Yt   = %15.5e \n", TacsRealPart(Yt));
  printf("Yc   = %15.5e \n", TacsRealPart(Yc));
  printf("S12  = %15.5e \n", TacsRealPart(S12));
  printf("C    = %15.5e \n", TacsRealPart(C));

  printf("\nStrain Failure Properties \n");
  printf("eXt  = %15.5e \n", TacsRealPart(eXt));
  printf("eXc  = %15.5e \n", TacsRealPart(eXc));
  printf("eYt  = %15.5e \n", TacsRealPart(eYt));
  printf("eYc  = %15.5e \n", TacsRealPart(eYc));
  printf("eS12 = %15.5e \n", TacsRealPart(eS12));

  printf("\nTsai-Wu tensor coefficients\n");
  printf("F1   = %15.5e \n", TacsRealPart(F1));
  printf("F2   = %15.5e \n", TacsRealPart(F2));
  printf("F11  = %15.5e \n", TacsRealPart(F11));
  printf("F12  = %15.5e \n", TacsRealPart(F12));
  printf("F22  = %15.5e \n", TacsRealPart(F22));
  printf("F66  = %15.5e \n", TacsRealPart(F66));
}

// First, the stress transformations
// Transform stress from the global to the local frame
void TACSOrthotropicPly::transformStressGlobal2Ply(TacsScalar angle,
                                                   const TacsScalar global[],
                                                   TacsScalar plyStress[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  plyStress[0] =
      cos2 * global[0] + sin2 * global[1] + 2.0 * sin1 * cos1 * global[2];
  plyStress[1] =
      sin2 * global[0] + cos2 * global[1] - 2.0 * sin1 * cos1 * global[2];
  plyStress[2] =
      sin1 * cos1 * (-global[0] + global[1]) + (cos2 - sin2) * global[2];
}

// Transform stress from the ply frame to the global frame
void TACSOrthotropicPly::transformStressPly2Global(TacsScalar angle,
                                                   const TacsScalar plyStress[],
                                                   TacsScalar global[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  global[0] = cos2 * plyStress[0] + sin2 * plyStress[1] -
              2.0 * sin1 * cos1 * plyStress[2];
  global[1] = sin2 * plyStress[0] + cos2 * plyStress[1] +
              2.0 * sin1 * cos1 * plyStress[2];
  global[2] = sin1 * cos1 * (plyStress[0] - plyStress[1]) +
              (cos2 - sin2) * plyStress[2];
}

// The sensitivity of the transformation of
void TACSOrthotropicPly::transformStressGlobal2PlyAngleSens(
    TacsScalar angle, const TacsScalar global[], TacsScalar plyStress[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = -2.0 * cos1 * sin1;           // sensitivity of cos**2
  TacsScalar s_sin2 = 2.0 * sin1 * cos1;            // sensitivity of sin**2
  TacsScalar s_sincos = cos1 * cos1 - sin1 * sin1;  // sensitivity of cos*sin

  plyStress[0] =
      s_cos2 * global[0] + s_sin2 * global[1] + 2.0 * s_sincos * global[2];
  plyStress[1] =
      s_sin2 * global[0] + s_cos2 * global[1] - 2.0 * s_sincos * global[2];
  plyStress[2] =
      s_sincos * (global[1] - global[0]) + (s_cos2 - s_sin2) * global[2];
}

// Transform stress from the ply frame to the global frame
void TACSOrthotropicPly::transformStressPly2GlobalAngleSens(
    TacsScalar angle, const TacsScalar plyStress[], TacsScalar global[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = -2.0 * cos1 * sin1;           // sensitivity of cos**2
  TacsScalar s_sin2 = 2.0 * sin1 * cos1;            // sensitivity of sin**2
  TacsScalar s_sincos = cos1 * cos1 - sin1 * sin1;  // sensitivity of cos*sin

  global[0] = s_cos2 * plyStress[0] + s_sin2 * plyStress[1] -
              2.0 * s_sincos * plyStress[2];
  global[1] = s_sin2 * plyStress[0] + s_cos2 * plyStress[1] +
              2.0 * s_sincos * plyStress[2];
  global[2] = s_sincos * (plyStress[0] - plyStress[1]) +
              (s_cos2 - s_sin2) * plyStress[2];
}

// Next, the strain transformations
// Transform strain from the global to the local frame
void TACSOrthotropicPly::transformStrainGlobal2Ply(TacsScalar angle,
                                                   const TacsScalar global[],
                                                   TacsScalar plyStrain[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  plyStrain[0] = cos2 * global[0] + sin2 * global[1] + sin1 * cos1 * global[2];
  plyStrain[1] = sin2 * global[0] + cos2 * global[1] - sin1 * cos1 * global[2];
  plyStrain[2] =
      2.0 * sin1 * cos1 * (global[1] - global[0]) + (cos2 - sin2) * global[2];
}

// Transform stress from the ply frame to the global frame
void TACSOrthotropicPly::transformStrainPly2Global(TacsScalar angle,
                                                   const TacsScalar plyStrain[],
                                                   TacsScalar global[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1 * cos1;
  TacsScalar sin2 = sin1 * sin1;

  global[0] =
      cos2 * plyStrain[0] + sin2 * plyStrain[1] - sin1 * cos1 * plyStrain[2];
  global[1] =
      sin2 * plyStrain[0] + cos2 * plyStrain[1] + sin1 * cos1 * plyStrain[2];
  global[2] = 2.0 * sin1 * cos1 * (plyStrain[0] - plyStrain[1]) +
              (cos2 - sin2) * plyStrain[2];
}

// The sensitivity of the transformation to the
void TACSOrthotropicPly::transformStrainGlobal2PlyAngleSens(
    TacsScalar angle, const TacsScalar global[], TacsScalar plyStrain[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = -2.0 * cos1 * sin1;           // sensitivity of cos**2
  TacsScalar s_sin2 = 2.0 * sin1 * cos1;            // sensitivity of sin**2
  TacsScalar s_sincos = cos1 * cos1 - sin1 * sin1;  // sensitivity of cos*sin

  plyStrain[0] = s_cos2 * global[0] + s_sin2 * global[1] + s_sincos * global[2];
  plyStrain[1] = s_sin2 * global[0] + s_cos2 * global[1] - s_sincos * global[2];
  plyStrain[2] =
      2.0 * s_sincos * (-global[0] + global[1]) + (s_cos2 - s_sin2) * global[2];
}

// Transform stress from the ply frame to the global frame
void TACSOrthotropicPly::transformStrainPly2GlobalAngleSens(
    TacsScalar angle, const TacsScalar plyStrain[], TacsScalar global[]) {
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = -2.0 * cos1 * sin1;           // sensitivity of cos**2
  TacsScalar s_sin2 = 2.0 * sin1 * cos1;            // sensitivity of sin**2
  TacsScalar s_sincos = cos1 * cos1 - sin1 * sin1;  // sensitivity of cos*sin

  global[0] =
      s_cos2 * plyStrain[0] + s_sin2 * plyStrain[1] - s_sincos * plyStrain[2];
  global[1] =
      s_sin2 * plyStrain[0] + s_cos2 * plyStrain[1] + s_sincos * plyStrain[2];
  global[2] = 2.0 * s_sincos * (plyStrain[0] - plyStrain[1]) +
              (s_cos2 - s_sin2) * plyStrain[2];
}

const char *TACSOrthotropicPly::getObjectName() { return name; }
