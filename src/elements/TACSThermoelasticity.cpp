/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSThermoelasticity.h"

TACSLinearThermoelasticity2D::TACSLinearThermoelasticity2D(
    TACSPlaneStressConstitutive *_stiff, ElementStrainType _strain_type,
    int _steady_state_flag) {
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
  steady_state_flag = _steady_state_flag;
}

TACSLinearThermoelasticity2D::~TACSLinearThermoelasticity2D() {
  stiff->decref();
}

// 0;   1;    2;   3;   4; 5;   6;    7;   8;   9;10;  11;   12;  13;  14
// u; u,t; u,tt; u,x; u,y; v; v,t; v,tt; v,x; v,y; T; T,t; T,tt; T,x; T,y

const int TACSLinearThermoelasticity2D::linear_Jac_pairs[] = {
    2, 2, 7, 7, 11, 11, 3, 3,  3, 4,  3,  8,  3,  9,  3,  10, 4,  3,
    4, 4, 4, 8, 4,  9,  4, 10, 8, 3,  8,  4,  8,  8,  8,  9,  8,  10,
    9, 3, 9, 4, 9,  8,  9, 9,  9, 10, 13, 13, 13, 14, 14, 13, 14, 14};

int TACSLinearThermoelasticity2D::getNumParameters() { return 2; }

int TACSLinearThermoelasticity2D::getVarsPerNode() { return 3; }

int TACSLinearThermoelasticity2D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearThermoelasticity2D::getDesignVarNums(int elemIndex, int dvLen,
                                                   int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearThermoelasticity2D::setDesignVars(int elemIndex, int dvLen,
                                                const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearThermoelasticity2D::getDesignVars(int elemIndex, int dvLen,
                                                TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearThermoelasticity2D::getDesignVarRange(int elemIndex, int dvLen,
                                                    TacsScalar lb[],
                                                    TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearThermoelasticity2D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;

  DUt[3] = 0.0;
  DUt[4] = 0.0;

  if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
    DUt[2] = 0.0;
    DUt[5] = 0.0;
  } else {
    DUt[2] = rho * Ut[2];
    DUt[5] = rho * Ut[5];
  }

  DUt[6] = 0.0;
  DUt[8] = 0.0;

  if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
    DUt[7] = 0.0;
  } else {
    DUt[7] = c * rho * Ut[7];
  }

  // Compute the thermal strain components
  TacsScalar theta = Ut[6];  // The temperature value
  TacsScalar et[3];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[3] - et[1];
    e[2] = Ux[1] + Ux[2] - et[2];
  } else {
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
  }

  // Evaluate the components of the stress
  TacsScalar s[3];
  stiff->evalStress(elemIndex, pt, X, e, s);
  DUx[0] = s[0];
  DUx[1] = s[2];

  DUx[2] = s[2];
  DUx[3] = s[1];

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2], flux[2];
  grad[0] = Ux[4];
  grad[1] = Ux[5];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[4] = flux[0];
  DUx[5] = flux[1];
}

void TACSLinearThermoelasticity2D::getWeakMatrixNonzeros(
    ElementMatrixType matType, int elemIndex, int *Jac_nnz,
    const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 27;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSLinearThermoelasticity2D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density and specific heat
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
    TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = 0.0;

    DUt[3] = 0.0;
    DUt[4] = 0.0;

    if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
      DUt[2] = 0.0;
      DUt[5] = 0.0;
    } else {
      DUt[2] = rho * Ut[2];
      DUt[5] = rho * Ut[5];
    }

    DUt[6] = 0.0;
    DUt[8] = 0.0;

    if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
      DUt[7] = 0.0;
    } else {
      DUt[7] = c * rho * Ut[7];
    }

    // Compute the thermal strain components
    TacsScalar theta = Ut[6];  // The temperature value
    TacsScalar et[3];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[3] - et[1];
      e[2] = Ux[1] + Ux[2] - et[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
    }

    // Evaluate the components of the stress
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);
    DUx[0] = s[0];
    DUx[1] = s[2];

    DUx[2] = s[2];
    DUx[3] = s[1];

    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[2], flux[2];
    grad[0] = Ux[4];
    grad[1] = Ux[5];

    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
    DUx[4] = flux[0];
    DUx[5] = flux[1];

    // Set the time-dependent terms
    if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
      Jac[0] = 0.0;
      Jac[1] = 0.0;
    } else {
      Jac[0] = rho;
      Jac[1] = rho;
    }

    if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
      Jac[2] = 0.0;
    } else {
      Jac[2] = c * rho;
    }

    // Compute the unit strain
    TacsScalar C[6], Kc[3];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);
    stiff->evalTangentStiffness(elemIndex, pt, X, C);
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    // Add the terms for linear thermoelasticity
    if (strain_type == TACS_LINEAR_STRAIN) {
      // s[0] = C[0]*(u,x - theta*et[0]) + C[1]*(v,y - theta*et[1])
      //      + C[2]*(u,y + v,x - theta*et[2])
      // s[1] = C[1]*(u,x) + C[3]*(v,y) + C[4]*(u,y + v,x)
      // s[2] = C[2]*(u,x - theta*et[0]) + C[4]*(v,y - theta*et[1])
      //      + C[5]*(u,y + v,x - theta*et[2])

      // i == 3 (s[0])
      Jac[3] = C[0];                                           // j == 3
      Jac[4] = C[2];                                           // j == 4
      Jac[5] = C[2];                                           // j == 8
      Jac[6] = C[1];                                           // j == 9
      Jac[7] = -(C[0] * et[0] + C[1] * et[1] + C[2] * et[2]);  // j == 10

      // i == 4 (s[2])
      Jac[8] = C[2];                                            // j == 3
      Jac[9] = C[5];                                            // j == 4
      Jac[10] = C[5];                                           // j == 8
      Jac[11] = C[4];                                           // j == 9
      Jac[12] = -(C[2] * et[0] + C[4] * et[1] + C[5] * et[2]);  // j == 10

      // i == 8 (s[2])
      Jac[13] = C[2];                                           // j == 3
      Jac[14] = C[5];                                           // j == 4
      Jac[15] = C[5];                                           // j == 8
      Jac[16] = C[4];                                           // j == 9
      Jac[17] = -(C[2] * et[0] + C[4] * et[1] + C[5] * et[2]);  // j == 10

      // i == 9 (s[1])
      Jac[18] = C[1];                                           // j == 3
      Jac[19] = C[4];                                           // j == 4
      Jac[20] = C[4];                                           // j == 8
      Jac[21] = C[3];                                           // j == 9
      Jac[22] = -(C[1] * et[0] + C[3] * et[1] + C[4] * et[2]);  // j == 10
    }

    // i == 13
    Jac[23] = Kc[0];  // j == 13
    Jac[24] = Kc[1];  // j == 14

    // i == 14
    Jac[25] = Kc[1];  // j == 13
    Jac[26] = Kc[2];  // j == 14
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSLinearThermoelasticity2D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar *dfdx) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  TacsScalar rho_coef =
      scale * (Ut[2] * Psi[0] + Ut[5] * Psi[1] + c * Ut[7] * Psi[2]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, dfdx);

  TacsScalar c_coef = scale * rho * Ut[7] * Psi[2];
  stiff->addSpecificHeatDVSens(elemIndex, c_coef, pt, X, dvLen, dfdx);

  // Compute the thermal strain components
  TacsScalar theta = Ut[6];  // The temperature value
  TacsScalar et[3];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[3] - et[1];
    e[2] = Ux[1] + Ux[2] - et[2];
  } else {
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
  }

  TacsScalar phi[3];
  phi[0] = Psix[0];
  phi[1] = Psix[3];
  phi[2] = Psix[1] + Psix[2];
  stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, dfdx);

  // Evaluate the components of the stress using the components of the
  // phi vector
  TacsScalar s[3];
  stiff->evalStress(elemIndex, pt, X, phi, s);
  stiff->addThermalStrainDVSens(elemIndex, pt, X, -theta * scale, s, dvLen,
                                dfdx);

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2];
  grad[0] = Ux[4];
  grad[1] = Ux[5];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  TacsScalar psi[2];
  psi[0] = Psix[4];
  psi[1] = Psix[5];
  stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, grad, psi, dvLen, dfdx);
}

void TACSLinearThermoelasticity2D::evalWeakAdjXptSensProduct(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
    TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
    TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;
  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;
  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = 0.0;
  dfdPsix[3] = dfdPsix[4] = dfdPsix[5] = 0.0;

  // Compute the material density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  // Compute the thermal strain components
  TacsScalar theta = Ut[6];  // The temperature value
  TacsScalar et[3];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[3] - et[1];
    e[2] = Ux[1] + Ux[2] - et[2];
  } else {
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
  }

  TacsScalar phi[3];
  phi[0] = Psix[0];
  phi[1] = Psix[3];
  phi[2] = Psix[1] + Psix[2];

  TacsScalar t1[3], t2[3];
  stiff->evalStress(elemIndex, pt, X, e, t1);
  stiff->evalStress(elemIndex, pt, X, phi, t2);

  *product = (rho * (Psi[0] * Ut[2] + Psi[1] * Ut[5] + c * Psi[2] * Ut[7]) +
              t2[0] * e[0] + t2[1] * e[1] + t2[2] * e[2]);

  if (strain_type == TACS_LINEAR_STRAIN) {
    dfdPsix[0] = t1[0];
    dfdPsix[3] = t1[1];
    dfdPsix[1] = t1[2];
    dfdPsix[2] = t1[2];

    dfdUx[0] = t2[0];
    dfdUx[3] = t2[1];
    dfdUx[1] = t2[2];
    dfdUx[2] = t2[2];
  }

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2];
  grad[0] = Ux[4];
  grad[1] = Ux[5];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  phi[0] = Psix[4];
  phi[1] = Psix[5];
  stiff->evalHeatFlux(elemIndex, pt, X, grad, t1);
  stiff->evalHeatFlux(elemIndex, pt, X, phi, t2);

  dfdPsix[4] = t1[0];
  dfdPsix[5] = t1[1];

  dfdUx[4] = t2[0];
  dfdUx[5] = t2[1];

  *product += t2[0] * grad[0] + t2[1] * grad[1];
}
/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearThermoelasticity2D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
    if (quantity) {
      // Compute the thermal strain components
      TacsScalar theta = Ut[6];  // The temperature value
      TacsScalar et[3];
      stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

      TacsScalar e[3];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0] - et[0];
        e[1] = Ux[3] - et[1];
        e[2] = Ux[1] + Ux[2] - et[2];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
        e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
        e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
      }

      *quantity = stiff->evalFailure(elemIndex, pt, X, e);
    }

    return 1;
  } else if (quantityType == TACS_HEAT_FLUX) {
    if (quantity) {
      // Compute the thermal flux from the thermal gradient
      TacsScalar grad[2], flux[2];
      grad[0] = Ux[4];
      grad[1] = Ux[5];

      // Add the flux components to the heat transfer portion
      // of the governing equations
      stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
      quantity[0] = flux[0];
      quantity[1] = flux[1];
    }

    return 2;
  } else if (quantityType == TACS_TEMPERATURE) {
    if (quantity) {
      *quantity = Ut[6];
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }

    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      // Compute the thermal strain components
      TacsScalar theta = Ut[6];  // The temperature value
      TacsScalar et[3];
      stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

      TacsScalar e[3];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0] - et[0];
        e[1] = Ux[3] - et[1];
        e[2] = Ux[1] + Ux[2] - et[2];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
        e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
        e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
      }

      // Evaluate the sress
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Evaluate the strain energy density
      *quantity = (e[0] * s[0] + e[1] * s[1] + e[2] * s[2]);
    }

    return 1;
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      TacsScalar e[3];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0];
        e[1] = Ux[3];
        e[2] = Ux[1] + Ux[2];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
        e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
        e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);
      }

      // Evaluate the sress
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Evaluate the strain energy density
      *quantity = (e[0] * s[0] + e[1] * s[1] + e[2] * s[2]);
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      // Return displacement components
      quantity[0] = Ut[0];
      quantity[1] = Ut[3];
    }

    return 2;
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    if (quantity) {
      TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
      quantity[0] = density * X[0];
      quantity[1] = density * X[1];
    }

    return 2;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
      quantity[0] = density * X[1] * X[1];   // Ixx
      quantity[1] = -density * X[1] * X[0];  // Ixy
      quantity[2] = density * X[0] * X[0];   // Iyy
    }

    return 3;
  }

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSLinearThermoelasticity2D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[6];  // The temperature value
    TacsScalar et[3];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[3] - et[1];
      e[2] = Ux[1] + Ux[2] - et[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
    }

    stiff->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, e, dvLen, dfdx);

    TacsScalar dfde[3];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, dfde);
    stiff->addThermalStrainDVSens(elemIndex, pt, X, -scale * dfdq[0] * theta,
                                  dfde, dvLen, dfdx);
  } else if (quantityType == TACS_HEAT_FLUX) {
    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[2];
    grad[0] = Ux[4];
    grad[1] = Ux[5];

    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, grad, dfdq, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[6];  // The temperature value
    TacsScalar et[3];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[3] - et[1];
      e[2] = Ux[1] + Ux[2] - et[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - et[0];
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - et[1];
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - et[2];
    }

    // Add the contributions to the derivative from the strain
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);
    stiff->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);
    stiff->addThermalStrainDVSens(
        elemIndex, pt, X, -2.0 * scale * dfdq[0] * theta, s, dvLen, dfdx);
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[3];
      e[2] = Ux[1] + Ux[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);
    }

    // Add the contributions to the derivative from the strain
    stiff->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar dfdmass = 0.0;
    dfdmass += scale * dfdq[0] * X[0];
    dfdmass += scale * dfdq[1] * X[1];
    stiff->addDensityDVSens(elemIndex, dfdmass, pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar dfdmass = 0.0;
    dfdmass += scale * dfdq[0] * X[1] * X[1];
    dfdmass += -scale * dfdq[1] * X[1] * X[0];
    dfdmass += scale * dfdq[2] * X[0] * X[0];
    stiff->addDensityDVSens(elemIndex, dfdmass, pt, X, dvLen, dfdx);
  }
}

/*
   Evaluate the derivatives of the point-wise quantity of interest
   with respect to X, Ut and Ux.
*/
void TACSLinearThermoelasticity2D::evalPointQuantitySens(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
    TacsScalar dfdUx[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = 0.0;
  dfdXd[2] = dfdXd[3] = 0.0;

  dfdUt[0] = dfdUt[1] = dfdUt[2] = 0.0;
  dfdUt[3] = dfdUt[4] = dfdUt[5] = 0.0;
  dfdUt[6] = dfdUt[7] = dfdUt[8] = 0.0;

  dfdUx[0] = dfdUx[1] = 0.0;
  dfdUx[2] = dfdUx[3] = 0.0;
  dfdUx[4] = dfdUx[5] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[6];  // The temperature value
    TacsScalar et[3];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);

    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - theta * et[0];
      e[1] = Ux[3] - theta * et[1];
      e[2] = Ux[1] + Ux[2] - theta * et[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - theta * et[0];
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - theta * et[1];
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - theta * et[2];
    }

    TacsScalar sens[3];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = dfdq[0] * sens[0];
      dfdUx[3] = dfdq[0] * sens[1];

      dfdUx[1] = dfdq[0] * sens[2];
      dfdUx[2] = dfdq[0] * sens[2];

      dfdUt[6] =
          -dfdq[0] * (sens[0] * et[0] + sens[1] * et[1] + sens[2] * et[2]);
    }
  } else if (quantityType == TACS_TEMPERATURE) {
    dfdUt[6] = dfdq[0];
  } else if (quantityType == TACS_HEAT_FLUX) {
    // flux = Kc*grad
    TacsScalar Kc[3];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    dfdUx[4] = dfdq[0] * Kc[0] + dfdq[1] * Kc[1];
    dfdUx[5] = dfdq[0] * Kc[1] + dfdq[1] * Kc[2];
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[6];  // The temperature value
    TacsScalar et[3];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);

    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - theta * et[0];
      e[1] = Ux[3] - theta * et[1];
      e[2] = Ux[1] + Ux[2] - theta * et[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]) - theta * et[0];
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]) - theta * et[1];
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]) - theta * et[2];
    }

    // Add the contributions to the derivative from the strain
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = 2.0 * dfdq[0] * s[0];
      dfdUx[3] = 2.0 * dfdq[0] * s[1];

      dfdUx[1] = 2.0 * dfdq[0] * s[2];
      dfdUx[2] = 2.0 * dfdq[0] * s[2];

      dfdUt[6] = -2.0 * dfdq[0] * (s[0] * et[0] + s[1] * et[1] + s[2] * et[2]);
    }
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[3];
      e[2] = Ux[1] + Ux[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);
    }

    // Add the contributions to the derivative from the strain
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = 2.0 * dfdq[0] * s[0];
      dfdUx[3] = 2.0 * dfdq[0] * s[1];

      dfdUx[1] = 2.0 * dfdq[0] * s[2];
      dfdUx[2] = 2.0 * dfdq[0] * s[2];
    }
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    dfdUt[0] = dfdq[0];
    dfdUt[3] = dfdq[1];
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
    dfdX[0] = density * dfdq[0];
    dfdX[1] = density * dfdq[1];
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
    dfdX[0] = density * (2.0 * dfdq[2] * X[0] - dfdq[1] * X[1]);
    dfdX[1] = density * (2.0 * dfdq[0] * X[1] - dfdq[1] * X[0]);
  }
}

/*
  Get the data for visualization at a given point
*/
void TACSLinearThermoelasticity2D::getOutputData(
    int elemIndex, const double time, ElementType etype, int write_flag,
    const double pt[], const TacsScalar X[], const TacsScalar Ut[],
    const TacsScalar Ux[], int ld_data, TacsScalar *data) {
  if (etype == TACS_PLANE_STRESS_ELEMENT) {
    if (write_flag & TACS_OUTPUT_NODES) {
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      data[0] = Ut[0];
      data[1] = Ut[3];
      data += 2;
    }

    TacsScalar e[3];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[3];
      e[2] = Ux[1] + Ux[2];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);
    }

    if (write_flag & TACS_OUTPUT_STRAINS) {
      data[0] = e[0];
      data[1] = e[1];
      data[2] = e[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_STRESSES) {
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);
      data[0] = s[0];
      data[1] = s[1];
      data[2] = s[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = stiff->evalFailure(elemIndex, pt, X, e);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}

TACSLinearThermoelasticity3D::TACSLinearThermoelasticity3D(
    TACSSolidConstitutive *_stiff, ElementStrainType _strain_type,
    int _steady_state_flag) {
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
  steady_state_flag = _steady_state_flag;
}

TACSLinearThermoelasticity3D::~TACSLinearThermoelasticity3D() {
  stiff->decref();
}

// 0;   1;    2;   3;   4;   5;
// u; u,t; u,tt; u,x; u,y; u,z;

// 6;   7;    8;   9;  10;  11;
// v; v,t; v,tt; v,x; v,y; v,z;

// 12;  13;   14;  15;  16;  17;
//  w; w,t; w,tt; w,x; w,y; w,z;

// 18;  19;   20;  21;  22;  23;
//  T; T,t; T,tt; T,x; T,y; T,z;

const int TACSLinearThermoelasticity3D::linear_Jac_pairs[] = {
    2,  2,  8,  8,  14, 14, 19, 19, 3,  3,  3,  4,  3,  5,  3,  9,  3,  10, 3,
    11, 3,  15, 3,  16, 3,  17, 3,  18, 4,  3,  4,  4,  4,  5,  4,  9,  4,  10,
    4,  11, 4,  15, 4,  16, 4,  17, 4,  18, 5,  3,  5,  4,  5,  5,  5,  9,  5,
    10, 5,  11, 5,  15, 5,  16, 5,  17, 5,  18, 9,  3,  9,  4,  9,  5,  9,  9,
    9,  10, 9,  11, 9,  15, 9,  16, 9,  17, 9,  18, 10, 3,  10, 4,  10, 5,  10,
    9,  10, 10, 10, 11, 10, 15, 10, 16, 10, 17, 10, 18, 11, 3,  11, 4,  11, 5,
    11, 9,  11, 10, 11, 11, 11, 15, 11, 16, 11, 17, 11, 18, 15, 3,  15, 4,  15,
    5,  15, 9,  15, 10, 15, 11, 15, 15, 15, 16, 15, 17, 15, 18, 16, 3,  16, 4,
    16, 5,  16, 9,  16, 10, 16, 11, 16, 15, 16, 16, 16, 17, 16, 18, 17, 3,  17,
    4,  17, 5,  17, 9,  17, 10, 17, 11, 17, 15, 17, 16, 17, 17, 17, 18, 21, 21,
    21, 22, 21, 23, 22, 21, 22, 22, 22, 23, 23, 21, 23, 22, 23, 23};

int TACSLinearThermoelasticity3D::getNumParameters() { return 3; }

int TACSLinearThermoelasticity3D::getVarsPerNode() { return 4; }

int TACSLinearThermoelasticity3D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearThermoelasticity3D::getDesignVarNums(int elemIndex, int dvLen,
                                                   int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearThermoelasticity3D::setDesignVars(int elemIndex, int dvLen,
                                                const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearThermoelasticity3D::getDesignVars(int elemIndex, int dvLen,
                                                TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearThermoelasticity3D::getDesignVarRange(int elemIndex, int dvLen,
                                                    TacsScalar lb[],
                                                    TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearThermoelasticity3D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;

  DUt[3] = 0.0;
  DUt[4] = 0.0;

  DUt[6] = 0.0;
  DUt[7] = 0.0;

  if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
    DUt[2] = 0.0;
    DUt[5] = 0.0;
    DUt[8] = 0.0;
  } else {
    DUt[2] = rho * Ut[2];
    DUt[5] = rho * Ut[5];
    DUt[8] = rho * Ut[8];
  }

  DUt[9] = 0.0;
  DUt[11] = 0.0;

  if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
    DUt[10] = 0.0;
  } else {
    DUt[10] = c * rho * Ut[10];
  }

  // Compute the thermal strain components
  TacsScalar theta = Ut[9];  // The temperature value
  TacsScalar et[6];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[6];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[4] - et[1];
    e[2] = Ux[8] - et[2];

    e[3] = Ux[5] + Ux[7] - et[3];
    e[4] = Ux[2] + Ux[6] - et[4];
    e[5] = Ux[1] + Ux[3] - et[5];
  } else {
    e[0] =
        Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
    e[1] =
        Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
    e[2] =
        Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

    e[3] =
        Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) - et[3];
    e[4] =
        Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) - et[4];
    e[5] =
        Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) - et[5];
  }

  // Evaluate the stress
  TacsScalar s[6];
  stiff->evalStress(elemIndex, pt, X, e, s);

  DUx[0] = s[0];
  DUx[1] = s[5];
  DUx[2] = s[4];

  DUx[3] = s[5];
  DUx[4] = s[1];
  DUx[5] = s[3];

  DUx[6] = s[4];
  DUx[7] = s[3];
  DUx[8] = s[2];

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[3], flux[3];
  grad[0] = Ux[9];
  grad[1] = Ux[10];
  grad[2] = Ux[11];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[9] = flux[0];
  DUx[10] = flux[1];
  DUx[11] = flux[2];
}

void TACSLinearThermoelasticity3D::getWeakMatrixNonzeros(
    ElementMatrixType matType, int elemIndex, int *Jac_nnz,
    const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 103;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSLinearThermoelasticity3D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density and specific heat
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
    TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = 0.0;

    DUt[3] = 0.0;
    DUt[4] = 0.0;

    DUt[6] = 0.0;
    DUt[7] = 0.0;

    if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
      DUt[2] = 0.0;
      DUt[5] = 0.0;
      DUt[8] = 0.0;
    } else {
      DUt[2] = rho * Ut[2];
      DUt[5] = rho * Ut[5];
      DUt[8] = rho * Ut[8];
    }

    DUt[9] = 0.0;
    DUt[11] = 0.0;

    if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
      DUt[10] = 0.0;
    } else {
      DUt[10] = c * rho * Ut[10];
    }

    // Compute the thermal strain components
    TacsScalar theta = Ut[9];  // The temperature value
    TacsScalar et[6];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[4] - et[1];
      e[2] = Ux[8] - et[2];

      e[3] = Ux[5] + Ux[7] - et[3];
      e[4] = Ux[2] + Ux[6] - et[4];
      e[5] = Ux[1] + Ux[3] - et[5];
    } else {
      e[0] =
          Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
      e[1] =
          Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
      e[2] =
          Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
             et[3];
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
             et[4];
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
             et[5];
    }

    // Evaluate the stress
    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    DUx[0] = s[0];
    DUx[1] = s[5];
    DUx[2] = s[4];

    DUx[3] = s[5];
    DUx[4] = s[1];
    DUx[5] = s[3];

    DUx[6] = s[4];
    DUx[7] = s[3];
    DUx[8] = s[2];

    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[3], flux[3];
    grad[0] = Ux[9];
    grad[1] = Ux[10];
    grad[2] = Ux[11];

    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
    DUx[9] = flux[0];
    DUx[10] = flux[1];
    DUx[11] = flux[2];

    if (steady_state_flag & TACS_STEADY_STATE_MECHANICAL) {
      Jac[0] = 0.0;
      Jac[1] = 0.0;
      Jac[2] = 0.0;
    } else {
      Jac[0] = rho;
      Jac[1] = rho;
      Jac[2] = rho;
    }

    if (steady_state_flag & TACS_STEADY_STATE_THERMAL) {
      Jac[3] = 0.0;
    } else {
      Jac[3] = c * rho;
    }

    // Compute the unit strain
    TacsScalar C[21], Kc[6];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);
    stiff->evalTangentStiffness(elemIndex, pt, X, C);
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    // Compute the thermal stress
    s[0] = C[0] * et[0] + C[1] * et[1] + C[2] * et[2] + C[3] * et[3] +
           C[4] * et[4] + C[5] * et[5];
    s[1] = C[1] * et[0] + C[6] * et[1] + C[7] * et[2] + C[8] * et[3] +
           C[9] * et[4] + C[10] * et[5];
    s[2] = C[2] * et[0] + C[7] * et[1] + C[11] * et[2] + C[12] * et[3] +
           C[13] * et[4] + C[14] * et[5];
    s[3] = C[3] * et[0] + C[8] * et[1] + C[12] * et[2] + C[15] * et[3] +
           C[16] * et[4] + C[17] * et[5];
    s[4] = C[4] * et[0] + C[9] * et[1] + C[13] * et[2] + C[16] * et[3] +
           C[18] * et[4] + C[19] * et[5];
    s[5] = C[5] * et[0] + C[10] * et[1] + C[14] * et[2] + C[17] * et[3] +
           C[19] * et[4] + C[20] * et[5];

    // Add the terms for linear thermoelasticity
    if (strain_type == TACS_LINEAR_STRAIN) {
      // s[0]
      Jac[4] = C[0];
      Jac[5] = C[5];
      Jac[6] = C[4];
      Jac[7] = C[5];
      Jac[8] = C[1];
      Jac[9] = C[3];
      Jac[10] = C[4];
      Jac[11] = C[3];
      Jac[12] = C[2];
      Jac[13] = -s[0];

      // s[5]
      Jac[14] = C[5];
      Jac[15] = C[20];
      Jac[16] = C[19];
      Jac[17] = C[20];
      Jac[18] = C[10];
      Jac[19] = C[17];
      Jac[20] = C[19];
      Jac[21] = C[17];
      Jac[22] = C[14];
      Jac[23] = -s[5];

      // s[4]
      Jac[24] = C[4];
      Jac[25] = C[19];
      Jac[26] = C[18];
      Jac[27] = C[19];
      Jac[28] = C[9];
      Jac[29] = C[16];
      Jac[30] = C[18];
      Jac[31] = C[16];
      Jac[32] = C[13];
      Jac[33] = -s[4];

      // s[5]
      Jac[34] = C[5];
      Jac[35] = C[20];
      Jac[36] = C[19];
      Jac[37] = C[20];
      Jac[38] = C[10];
      Jac[39] = C[17];
      Jac[40] = C[19];
      Jac[41] = C[17];
      Jac[42] = C[14];
      Jac[43] = -s[5];

      // s[1]
      Jac[44] = C[1];
      Jac[45] = C[10];
      Jac[46] = C[9];
      Jac[47] = C[10];
      Jac[48] = C[6];
      Jac[49] = C[8];
      Jac[50] = C[9];
      Jac[51] = C[8];
      Jac[52] = C[7];
      Jac[53] = -s[1];

      // s[3]
      Jac[54] = C[3];
      Jac[55] = C[17];
      Jac[56] = C[16];
      Jac[57] = C[17];
      Jac[58] = C[8];
      Jac[59] = C[15];
      Jac[60] = C[16];
      Jac[61] = C[15];
      Jac[62] = C[12];
      Jac[63] = -s[3];

      // s[4]
      Jac[64] = C[4];
      Jac[65] = C[19];
      Jac[66] = C[18];
      Jac[67] = C[19];
      Jac[68] = C[9];
      Jac[69] = C[16];
      Jac[70] = C[18];
      Jac[71] = C[16];
      Jac[72] = C[13];
      Jac[73] = -s[4];

      // s[3]
      Jac[74] = C[3];
      Jac[75] = C[17];
      Jac[76] = C[16];
      Jac[77] = C[17];
      Jac[78] = C[8];
      Jac[79] = C[15];
      Jac[80] = C[16];
      Jac[81] = C[15];
      Jac[82] = C[12];
      Jac[83] = -s[3];

      // s[2]
      Jac[84] = C[2];
      Jac[85] = C[14];
      Jac[86] = C[13];
      Jac[87] = C[14];
      Jac[88] = C[7];
      Jac[89] = C[12];
      Jac[90] = C[13];
      Jac[91] = C[12];
      Jac[92] = C[11];
      Jac[93] = -s[2];
    }

    Jac[94] = Kc[0];
    Jac[95] = Kc[1];
    Jac[96] = Kc[2];

    Jac[97] = Kc[1];
    Jac[98] = Kc[3];
    Jac[99] = Kc[4];

    Jac[100] = Kc[2];
    Jac[101] = Kc[4];
    Jac[102] = Kc[5];
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSLinearThermoelasticity3D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar dfdx[]) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  TacsScalar rho_coef = scale * (Ut[2] * Psi[0] + Ut[5] * Psi[1] +
                                 Ut[8] * Psi[2] + c * Ut[10] * Psi[3]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, dfdx);

  TacsScalar c_coef = scale * rho * Ut[10] * Psi[3];
  stiff->addSpecificHeatDVSens(elemIndex, c_coef, pt, X, dvLen, dfdx);

  // Compute the thermal strain components
  TacsScalar theta = Ut[9];  // The temperature value
  TacsScalar et[6];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[6];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[4] - et[1];
    e[2] = Ux[8] - et[2];

    e[3] = Ux[5] + Ux[7] - et[3];
    e[4] = Ux[2] + Ux[6] - et[4];
    e[5] = Ux[1] + Ux[3] - et[5];
  } else {
    e[0] =
        Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
    e[1] =
        Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
    e[2] =
        Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

    e[3] =
        Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) - et[3];
    e[4] =
        Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) - et[4];
    e[5] =
        Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) - et[5];
  }

  TacsScalar phi[6];
  phi[0] = Psix[0];
  phi[1] = Psix[4];
  phi[2] = Psix[8];
  phi[3] = Psix[5] + Psix[7];
  phi[4] = Psix[2] + Psix[6];
  phi[5] = Psix[1] + Psix[3];
  stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, dfdx);

  // Evaluate the components of the stress using the components of the
  // phi vector
  TacsScalar s[6];
  stiff->evalStress(elemIndex, pt, X, phi, s);
  stiff->addThermalStrainDVSens(elemIndex, pt, X, -theta * scale, s, dvLen,
                                dfdx);

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[3];
  grad[0] = Ux[9];
  grad[1] = Ux[10];
  grad[2] = Ux[11];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  TacsScalar psi[3];
  psi[0] = Psix[9];
  psi[1] = Psix[10];
  psi[2] = Psix[11];
  stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, grad, psi, dvLen, dfdx);
}

void TACSLinearThermoelasticity3D::evalWeakAdjXptSensProduct(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
    TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
    TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;
  dfdXd[6] = dfdXd[7] = dfdXd[8] = 0.0;

  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;
  dfdUx[9] = dfdUx[10] = dfdUx[11] = 0.0;

  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = 0.0;
  dfdPsix[3] = dfdPsix[4] = dfdPsix[5] = 0.0;
  dfdPsix[6] = dfdPsix[7] = dfdPsix[8] = 0.0;
  dfdPsix[9] = dfdPsix[10] = dfdPsix[11] = 0.0;

  // Compute the material density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  // Compute the thermal strain components
  TacsScalar theta = Ut[9];  // The temperature value
  TacsScalar et[6];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[6];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0] - et[0];
    e[1] = Ux[4] - et[1];
    e[2] = Ux[8] - et[2];

    e[3] = Ux[5] + Ux[7] - et[3];
    e[4] = Ux[2] + Ux[6] - et[4];
    e[5] = Ux[1] + Ux[3] - et[5];
  } else {
    e[0] =
        Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
    e[1] =
        Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
    e[2] =
        Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

    e[3] =
        Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) - et[3];
    e[4] =
        Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) - et[4];
    e[5] =
        Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) - et[5];
  }

  TacsScalar phi[6];
  phi[0] = Psix[0];
  phi[1] = Psix[4];
  phi[2] = Psix[8];
  phi[3] = Psix[5] + Psix[7];
  phi[4] = Psix[2] + Psix[6];
  phi[5] = Psix[1] + Psix[3];

  TacsScalar t1[6], t2[6];
  stiff->evalStress(elemIndex, pt, X, e, t1);
  stiff->evalStress(elemIndex, pt, X, phi, t2);

  *product = (rho * (Psi[0] * Ut[2] + Psi[1] * Ut[5] + Psi[2] * Ut[8] +
                     c * Psi[3] * Ut[10]) +
              (t2[0] * e[0] + t2[1] * e[1] + t2[2] * e[2] + t2[3] * e[3] +
               t2[4] * e[4] + t2[5] * e[5]));

  if (strain_type == TACS_LINEAR_STRAIN) {
    dfdUx[0] = t2[0];
    dfdUx[4] = t2[1];
    dfdUx[8] = t2[2];
    dfdUx[5] = t2[3];
    dfdUx[7] = t2[3];
    dfdUx[2] = t2[4];
    dfdUx[6] = t2[4];
    dfdUx[1] = t2[5];
    dfdUx[3] = t2[5];

    dfdPsix[0] = t1[0];
    dfdPsix[4] = t1[1];
    dfdPsix[8] = t1[2];
    dfdPsix[5] = t1[3];
    dfdPsix[7] = t1[3];
    dfdPsix[2] = t1[4];
    dfdPsix[6] = t1[4];
    dfdPsix[1] = t1[5];
    dfdPsix[3] = t1[5];
  }

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[3];
  grad[0] = Ux[9];
  grad[1] = Ux[10];
  grad[2] = Ux[11];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  phi[0] = Psix[9];
  phi[1] = Psix[10];
  phi[2] = Psix[11];
  stiff->evalHeatFlux(elemIndex, pt, X, grad, t1);
  stiff->evalHeatFlux(elemIndex, pt, X, phi, t2);

  dfdPsix[9] = t1[0];
  dfdPsix[10] = t1[1];
  dfdPsix[11] = t1[2];

  dfdUx[9] = t2[0];
  dfdUx[10] = t2[1];
  dfdUx[11] = t2[2];

  *product += t2[0] * grad[0] + t2[1] * grad[1] + t2[2] * grad[2];
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearThermoelasticity3D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
    if (quantity) {
      // Compute the thermal strain components
      TacsScalar theta = Ut[9];  // The temperature value
      TacsScalar et[6];
      stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

      // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
      TacsScalar e[6];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0] - et[0];
        e[1] = Ux[4] - et[1];
        e[2] = Ux[8] - et[2];

        e[3] = Ux[5] + Ux[7] - et[3];
        e[4] = Ux[2] + Ux[6] - et[4];
        e[5] = Ux[1] + Ux[3] - et[5];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) -
               et[0];
        e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) -
               et[1];
        e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) -
               et[2];

        e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
               et[3];
        e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
               et[4];
        e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
               et[5];
      }

      *quantity = stiff->evalFailure(elemIndex, pt, X, e);
    }

    return 1;
  } else if (quantityType == TACS_HEAT_FLUX) {
    if (quantity) {
      // Compute the thermal flux from the thermal gradient
      TacsScalar grad[3], flux[3];
      grad[0] = Ux[9];
      grad[1] = Ux[10];
      grad[2] = Ux[11];

      // Add the flux components to the heat transfer portion
      // of the governing equations
      stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
      quantity[0] = flux[0];
      quantity[1] = flux[1];
      quantity[2] = flux[2];
    }

    return 3;
  } else if (quantityType == TACS_TEMPERATURE) {
    if (quantity) {
      *quantity = Ut[9];
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }

    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      // Compute the thermal strain components
      TacsScalar theta = Ut[9];  // The temperature value
      TacsScalar et[6];
      stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

      // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
      TacsScalar e[6];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0] - et[0];
        e[1] = Ux[4] - et[1];
        e[2] = Ux[8] - et[2];

        e[3] = Ux[5] + Ux[7] - et[3];
        e[4] = Ux[2] + Ux[6] - et[4];
        e[5] = Ux[1] + Ux[3] - et[5];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) -
               et[0];
        e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) -
               et[1];
        e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) -
               et[2];

        e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
               et[3];
        e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
               et[4];
        e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
               et[5];
      }

      TacsScalar s[6];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Evaluate the strain energy density
      *quantity = (e[0] * s[0] + e[1] * s[1] + e[2] * s[2] + e[3] * s[3] +
                   e[4] * s[4] + e[5] * s[5]);
    }

    return 1;
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      // Compute the total strain e = 0.5*(u,x + u,x^{T})
      TacsScalar e[6];
      if (strain_type == TACS_LINEAR_STRAIN) {
        e[0] = Ux[0];
        e[1] = Ux[4];
        e[2] = Ux[8];

        e[3] = Ux[5] + Ux[7];
        e[4] = Ux[2] + Ux[6];
        e[5] = Ux[1] + Ux[3];
      } else {
        e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
        e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
        e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

        e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
        e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
        e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
      }

      TacsScalar s[6];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Evaluate the strain energy density
      *quantity = (e[0] * s[0] + e[1] * s[1] + e[2] * s[2] + e[3] * s[3] +
                   e[4] * s[4] + e[5] * s[5]);
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      // Return displacement components
      quantity[0] = Ut[0];
      quantity[1] = Ut[3];
      quantity[2] = Ut[6];
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    if (quantity) {
      TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
      quantity[0] = density * X[0];
      quantity[1] = density * X[1];
      quantity[2] = density * X[2];
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
      quantity[0] = density * (X[1] * X[1] + X[2] * X[2]);  // Ixx
      quantity[1] = -density * X[1] * X[0];                 // Ixy
      quantity[2] = -density * X[2] * X[0];                 // Ixz
      quantity[3] = density * (X[0] * X[0] + X[2] * X[2]);  // Iyy
      quantity[4] = -density * X[1] * X[2];                 // Iyz
      quantity[5] = density * (X[1] * X[1] + X[0] * X[0]);  // Izz
    }

    return 6;
  }

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSLinearThermoelasticity3D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[9];  // The temperature value
    TacsScalar et[6];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[4] - et[1];
      e[2] = Ux[8] - et[2];

      e[3] = Ux[5] + Ux[7] - et[3];
      e[4] = Ux[2] + Ux[6] - et[4];
      e[5] = Ux[1] + Ux[3] - et[5];
    } else {
      e[0] =
          Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
      e[1] =
          Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
      e[2] =
          Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
             et[3];
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
             et[4];
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
             et[5];
    }

    stiff->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, e, dvLen, dfdx);

    TacsScalar dfde[6];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, dfde);
    stiff->addThermalStrainDVSens(elemIndex, pt, X, -scale * dfdq[0] * theta,
                                  dfde, dvLen, dfdx);
  } else if (quantityType == TACS_HEAT_FLUX) {
    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[3];
    grad[0] = Ux[9];
    grad[1] = Ux[10];
    grad[2] = Ux[11];

    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, grad, dfdq, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[9];  // The temperature value
    TacsScalar et[6];
    stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - et[0];
      e[1] = Ux[4] - et[1];
      e[2] = Ux[8] - et[2];

      e[3] = Ux[5] + Ux[7] - et[3];
      e[4] = Ux[2] + Ux[6] - et[4];
      e[5] = Ux[1] + Ux[3] - et[5];
    } else {
      e[0] =
          Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) - et[0];
      e[1] =
          Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) - et[1];
      e[2] =
          Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) - et[2];

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
             et[3];
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
             et[4];
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
             et[5];
    }

    // Add the contributions to the derivative from the strain
    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);
    stiff->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);

    // Add the result from the derivative of the thermal strain
    stiff->addThermalStrainDVSens(
        elemIndex, pt, X, -2.0 * scale * dfdq[0] * theta, s, dvLen, dfdx);
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    // Compute the total strain e = 0.5*(u,x + u,x^{T})
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[4];
      e[2] = Ux[8];

      e[3] = Ux[5] + Ux[7];
      e[4] = Ux[2] + Ux[6];
      e[5] = Ux[1] + Ux[3];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
    }

    stiff->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar dfdmass = 0.0;
    dfdmass += scale * dfdq[0] * X[0];
    dfdmass += scale * dfdq[1] * X[1];
    dfdmass += scale * dfdq[2] * X[2];
    stiff->addDensityDVSens(elemIndex, dfdmass, pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar dfdmass = 0.0;
    dfdmass += scale * dfdq[0] * (X[1] * X[1] + X[2] * X[2]);
    dfdmass += -scale * dfdq[1] * X[1] * X[0];
    dfdmass += -scale * dfdq[2] * X[2] * X[0];
    dfdmass += scale * dfdq[3] * (X[0] * X[0] + X[2] * X[2]);
    dfdmass += -scale * dfdq[4] * X[1] * X[2];
    dfdmass += scale * dfdq[5] * (X[1] * X[1] + X[0] * X[0]);
    stiff->addDensityDVSens(elemIndex, dfdmass, pt, X, dvLen, dfdx);
  }
}

/*
   Evaluate the derivatives of the point-wise quantity of interest
   with respect to X, Ut and Ux.
*/
void TACSLinearThermoelasticity3D::evalPointQuantitySens(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
    TacsScalar dfdUx[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;
  dfdXd[6] = dfdXd[7] = dfdXd[8] = 0.0;

  dfdUt[0] = dfdUt[1] = dfdUt[2] = 0.0;
  dfdUt[3] = dfdUt[4] = dfdUt[5] = 0.0;
  dfdUt[6] = dfdUt[7] = dfdUt[8] = 0.0;
  dfdUt[9] = dfdUt[10] = dfdUt[11] = 0.0;

  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;
  dfdUx[9] = dfdUx[10] = dfdUx[11] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[9];  // The temperature value
    TacsScalar et[6];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - theta * et[0];
      e[1] = Ux[4] - theta * et[1];
      e[2] = Ux[8] - theta * et[2];

      e[3] = Ux[5] + Ux[7] - theta * et[3];
      e[4] = Ux[2] + Ux[6] - theta * et[4];
      e[5] = Ux[1] + Ux[3] - theta * et[5];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) -
             theta * et[0];
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) -
             theta * et[1];
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) -
             theta * et[2];

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
             theta * et[3];
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
             theta * et[4];
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
             theta * et[5];
    }

    TacsScalar sens[6];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = dfdq[0] * sens[0];
      dfdUx[4] = dfdq[0] * sens[1];
      dfdUx[8] = dfdq[0] * sens[2];

      dfdUx[5] = dfdq[0] * sens[3];
      dfdUx[7] = dfdq[0] * sens[3];

      dfdUx[2] = dfdq[0] * sens[4];
      dfdUx[6] = dfdq[0] * sens[4];

      dfdUx[1] = dfdq[0] * sens[5];
      dfdUx[3] = dfdq[0] * sens[5];

      dfdUt[9] =
          -dfdq[0] * (sens[0] * et[0] + sens[1] * et[1] + sens[2] * et[2] +
                      sens[3] * et[3] + sens[4] * et[4] + sens[5] * et[5]);
    }
  } else if (quantityType == TACS_TEMPERATURE) {
    dfdUt[9] = dfdq[0];
  } else if (quantityType == TACS_HEAT_FLUX) {
    // flux = Kc*grad
    TacsScalar Kc[6];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    dfdUx[9] = dfdq[0] * Kc[0] + dfdq[1] * Kc[1] + dfdq[2] * Kc[2];
    dfdUx[10] = dfdq[0] * Kc[1] + dfdq[1] * Kc[3] + dfdq[2] * Kc[4];
    dfdUx[11] = dfdq[0] * Kc[2] + dfdq[1] * Kc[4] + dfdq[2] * Kc[5];
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the thermal strain components
    TacsScalar theta = Ut[9];  // The temperature value
    TacsScalar et[6];
    stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);

    // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0] - theta * et[0];
      e[1] = Ux[4] - theta * et[1];
      e[2] = Ux[8] - theta * et[2];

      e[3] = Ux[5] + Ux[7] - theta * et[3];
      e[4] = Ux[2] + Ux[6] - theta * et[4];
      e[5] = Ux[1] + Ux[3] - theta * et[5];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]) -
             theta * et[0];
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]) -
             theta * et[1];
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]) -
             theta * et[2];

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]) -
             theta * et[3];
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]) -
             theta * et[4];
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]) -
             theta * et[5];
    }

    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = 2.0 * dfdq[0] * s[0];
      dfdUx[4] = 2.0 * dfdq[0] * s[1];
      dfdUx[8] = 2.0 * dfdq[0] * s[2];

      dfdUx[5] = 2.0 * dfdq[0] * s[3];
      dfdUx[7] = 2.0 * dfdq[0] * s[3];

      dfdUx[2] = 2.0 * dfdq[0] * s[4];
      dfdUx[6] = 2.0 * dfdq[0] * s[4];

      dfdUx[1] = 2.0 * dfdq[0] * s[5];
      dfdUx[3] = 2.0 * dfdq[0] * s[5];

      dfdUt[9] = -2.0 * dfdq[0] *
                 (s[0] * et[0] + s[1] * et[1] + s[2] * et[2] + s[3] * et[3] +
                  s[4] * et[4] + s[5] * et[5]);
    }
  } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
    // Compute the total strain e = 0.5*(u,x + u,x^{T})
    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[4];
      e[2] = Ux[8];

      e[3] = Ux[5] + Ux[7];
      e[4] = Ux[2] + Ux[6];
      e[5] = Ux[1] + Ux[3];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
    }

    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = 2.0 * dfdq[0] * s[0];
      dfdUx[4] = 2.0 * dfdq[0] * s[1];
      dfdUx[8] = 2.0 * dfdq[0] * s[2];

      dfdUx[5] = 2.0 * dfdq[0] * s[3];
      dfdUx[7] = 2.0 * dfdq[0] * s[3];

      dfdUx[2] = 2.0 * dfdq[0] * s[4];
      dfdUx[6] = 2.0 * dfdq[0] * s[4];

      dfdUx[1] = 2.0 * dfdq[0] * s[5];
      dfdUx[3] = 2.0 * dfdq[0] * s[5];
    }
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    dfdUt[0] = dfdq[0];
    dfdUt[3] = dfdq[1];
    dfdUt[6] = dfdq[2];
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
    dfdX[0] = density * dfdq[0];
    dfdX[1] = density * dfdq[1];
    dfdX[2] = density * dfdq[2];
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar density = stiff->evalDensity(elemIndex, pt, X);
    dfdX[0] = density * (2.0 * X[0] * (dfdq[3] + dfdq[5]) - dfdq[1] * X[1] -
                         dfdq[2] * X[2]);
    dfdX[1] = density * (2.0 * X[1] * (dfdq[0] + dfdq[5]) - dfdq[1] * X[0] -
                         dfdq[4] * X[2]);
    dfdX[2] = density * (2.0 * X[2] * (dfdq[0] + dfdq[3]) - dfdq[2] * X[0] -
                         dfdq[4] * X[1]);
  }
}

/*
  Get the data for visualization at a given point
*/
void TACSLinearThermoelasticity3D::getOutputData(
    int elemIndex, const double time, ElementType etype, int write_flag,
    const double pt[], const TacsScalar X[], const TacsScalar Ut[],
    const TacsScalar Ux[], int ld_data, TacsScalar *data) {
  if (etype == TACS_SOLID_ELEMENT) {
    if (write_flag & TACS_OUTPUT_NODES) {
      // doesn't this depend whether it's linear/quadratic/etc?
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      data[0] = Ut[0];
      data[1] = Ut[3];
      data[3] = Ut[6];
      data += 3;
    }

    TacsScalar e[6];
    if (strain_type == TACS_LINEAR_STRAIN) {
      e[0] = Ux[0];
      e[1] = Ux[4];
      e[2] = Ux[8];

      e[3] = Ux[5] + Ux[7];
      e[4] = Ux[2] + Ux[6];
      e[5] = Ux[1] + Ux[3];
    } else {
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
    }

    if (write_flag & TACS_OUTPUT_STRAINS) {
      data[0] = e[0];
      data[1] = e[1];
      data[2] = e[2];
      data[3] = e[3];
      data[4] = e[4];
      data[5] = e[5];
      data += 6;
    }
    if (write_flag & TACS_OUTPUT_STRESSES) {
      TacsScalar s[6];
      stiff->evalStress(elemIndex, pt, X, e, s);
      data[0] = s[0];
      data[1] = s[1];
      data[2] = s[2];
      data[3] = s[3];
      data[4] = s[4];
      data[5] = s[5];
      data += 6;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = stiff->evalFailure(elemIndex, pt, X, e);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}
