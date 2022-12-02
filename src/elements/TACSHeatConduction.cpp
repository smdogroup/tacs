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

#include "TACSHeatConduction.h"

TACSHeatConduction2D::TACSHeatConduction2D(
    TACSPlaneStressConstitutive *_stiff) {
  stiff = _stiff;
  stiff->incref();
}

TACSHeatConduction2D::~TACSHeatConduction2D() { stiff->decref(); }

// 0;   1;    2;   3;   4;
// T; T,t; T,tt; T,x; T,y

const int TACSHeatConduction2D::linear_Jac_pairs[] = {1, 1, 3, 3, 3,
                                                      4, 4, 3, 4, 4};

int TACSHeatConduction2D::getNumParameters() { return 2; }

int TACSHeatConduction2D::getVarsPerNode() { return 1; }

int TACSHeatConduction2D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSHeatConduction2D::getDesignVarNums(int elemIndex, int dvLen,
                                           int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSHeatConduction2D::setDesignVars(int elemIndex, int dvLen,
                                        const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSHeatConduction2D::getDesignVars(int elemIndex, int dvLen,
                                        TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSHeatConduction2D::getDesignVarRange(int elemIndex, int dvLen,
                                            TacsScalar lb[], TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSHeatConduction2D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = c * rho * Ut[1];
  DUt[2] = 0.0;

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2], flux[2];
  grad[0] = Ux[0];
  grad[1] = Ux[1];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[0] = flux[0];
  DUx[1] = flux[1];
}

void TACSHeatConduction2D::getWeakMatrixNonzeros(ElementMatrixType matType,
                                                 int elemIndex, int *Jac_nnz,
                                                 const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 5;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSHeatConduction2D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density and specific heat
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
    TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = c * rho * Ut[1];
    DUt[2] = 0.0;

    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[2], flux[2];
    grad[0] = Ux[0];
    grad[1] = Ux[1];

    stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
    DUx[0] = flux[0];
    DUx[1] = flux[1];

    // Set the time-dependent terms
    Jac[0] = c * rho;

    // Compute the unit strain
    TacsScalar Kc[3];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    Jac[1] = Kc[0];
    Jac[2] = Kc[1];
    Jac[3] = Kc[1];
    Jac[4] = Kc[2];
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSHeatConduction2D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar *dfdx) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  TacsScalar rho_coef = scale * (c * Ut[1] * Psi[0]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, dfdx);

  TacsScalar c_coef = scale * rho * Ut[1] * Psi[0];
  stiff->addSpecificHeatDVSens(elemIndex, c_coef, pt, X, dvLen, dfdx);

  stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, Ux, Psix, dvLen, dfdx);
}

void TACSHeatConduction2D::evalWeakAdjXptSensProduct(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
    TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
    TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;

  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  stiff->evalHeatFlux(elemIndex, pt, X, Ux, dfdPsix);
  stiff->evalHeatFlux(elemIndex, pt, X, Psix, dfdUx);

  *product = c * rho * Psi[0] * Ut[1] + dfdUx[0] * Ux[0] + dfdUx[1] * Ux[1];
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSHeatConduction2D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_HEAT_FLUX) {
    if (quantity) {
      stiff->evalHeatFlux(elemIndex, pt, X, Ux, quantity);
    }

    return 2;
  } else if (quantityType == TACS_TEMPERATURE) {
    if (quantity) {
      *quantity = Ut[0];
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }

    return 1;
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
void TACSHeatConduction2D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_HEAT_FLUX) {
    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, Ux, dfdq, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
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
void TACSHeatConduction2D::evalPointQuantitySens(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
    TacsScalar dfdUx[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = 0.0;
  dfdXd[2] = dfdXd[3] = 0.0;

  dfdUt[0] = dfdUt[1] = dfdUt[2] = 0.0;

  dfdUx[0] = dfdUx[1] = 0.0;

  if (quantityType == TACS_TEMPERATURE) {
    dfdUt[0] = dfdq[0];
  } else if (quantityType == TACS_HEAT_FLUX) {
    // flux = Kc*grad
    TacsScalar Kc[3];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    dfdUx[0] = dfdq[0] * Kc[0] + dfdq[1] * Kc[1];
    dfdUx[1] = dfdq[0] * Kc[1] + dfdq[1] * Kc[2];
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
void TACSHeatConduction2D::getOutputData(
    int elemIndex, const double time, ElementType etype, int write_flag,
    const double pt[], const TacsScalar X[], const TacsScalar Ut[],
    const TacsScalar Ux[], int ld_data, TacsScalar *data) {
  if (etype == TACS_SCALAR_2D_ELEMENT) {
    if (write_flag & TACS_OUTPUT_NODES) {
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      data[0] = Ut[0];
      data += 1;
    }

    TacsScalar grad[2];
    grad[0] = Ux[0];
    grad[1] = Ux[1];

    if (write_flag & TACS_OUTPUT_STRAINS) {
      data[0] = grad[0];
      data[1] = grad[1];
      data += 2;
    }
    if (write_flag & TACS_OUTPUT_STRESSES) {
      TacsScalar flux[2];
      stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
      data[0] = flux[0];
      data[1] = flux[1];
      data += 2;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = stiff->evalDensity(elemIndex, pt, X);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}

TACSHeatConduction3D::TACSHeatConduction3D(TACSSolidConstitutive *_stiff) {
  stiff = _stiff;
  stiff->incref();
}

TACSHeatConduction3D::~TACSHeatConduction3D() { stiff->decref(); }

// 0;   1;    2;   3;   4;   5;
// T; T,t; T,tt; T,x; T,y; T,z;

const int TACSHeatConduction3D::linear_Jac_pairs[] = {
    1, 1, 3, 3, 3, 4, 3, 5, 4, 3, 4, 4, 4, 5, 5, 3, 5, 4, 5, 5};

int TACSHeatConduction3D::getNumParameters() { return 3; }

int TACSHeatConduction3D::getVarsPerNode() { return 1; }

int TACSHeatConduction3D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSHeatConduction3D::getDesignVarNums(int elemIndex, int dvLen,
                                           int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSHeatConduction3D::setDesignVars(int elemIndex, int dvLen,
                                        const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSHeatConduction3D::getDesignVars(int elemIndex, int dvLen,
                                        TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSHeatConduction3D::getDesignVarRange(int elemIndex, int dvLen,
                                            TacsScalar lb[], TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSHeatConduction3D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = c * rho * Ut[1];
  DUt[2] = 0.0;

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[3], flux[3];
  grad[0] = Ux[0];
  grad[1] = Ux[1];
  grad[2] = Ux[2];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[0] = flux[0];
  DUx[1] = flux[1];
  DUx[2] = flux[2];
}

void TACSHeatConduction3D::getWeakMatrixNonzeros(ElementMatrixType matType,
                                                 int elemIndex, int *Jac_nnz,
                                                 const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 10;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSHeatConduction3D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density and specific heat
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
    TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = c * rho * Ut[1];
    DUt[2] = 0.0;

    // Compute the thermal flux from the thermal gradient
    TacsScalar grad[3], flux[3];
    grad[0] = Ux[0];
    grad[1] = Ux[1];
    grad[2] = Ux[2];

    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
    DUx[0] = flux[0];
    DUx[1] = flux[1];
    DUx[2] = flux[2];

    // Set the time-dependent terms
    Jac[0] = c * rho;

    // Compute the unit strain
    TacsScalar Kc[6];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    Jac[1] = Kc[0];
    Jac[2] = Kc[1];
    Jac[3] = Kc[2];

    Jac[4] = Kc[1];
    Jac[5] = Kc[3];
    Jac[6] = Kc[4];

    Jac[7] = Kc[2];
    Jac[8] = Kc[4];
    Jac[9] = Kc[5];
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSHeatConduction3D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar dfdx[]) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  TacsScalar rho_coef = scale * (c * Ut[1] * Psi[0]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, dfdx);

  TacsScalar c_coef = scale * rho * Ut[1] * Psi[0];
  stiff->addSpecificHeatDVSens(elemIndex, c_coef, pt, X, dvLen, dfdx);

  stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, Ux, Psix, dvLen, dfdx);
}

void TACSHeatConduction3D::evalWeakAdjXptSensProduct(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
    TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
    TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;
  dfdXd[6] = dfdXd[7] = dfdXd[8] = 0.0;

  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  stiff->evalHeatFlux(elemIndex, pt, X, Ux, dfdPsix);
  stiff->evalHeatFlux(elemIndex, pt, X, Psix, dfdUx);

  *product = c * rho * Psi[0] * Ut[1] + dfdUx[0] * Ux[0] + dfdUx[1] * Ux[1] +
             dfdUx[2] * Ux[2];
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSHeatConduction3D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_HEAT_FLUX) {
    if (quantity) {
      // Add the flux components to the heat transfer portion
      // of the governing equations
      stiff->evalHeatFlux(elemIndex, pt, X, Ux, quantity);
    }

    return 3;
  } else if (quantityType == TACS_TEMPERATURE) {
    if (quantity) {
      *quantity = Ut[0];
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }

    return 1;
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
void TACSHeatConduction3D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_HEAT_FLUX) {
    // Add the flux components to the heat transfer portion
    // of the governing equations
    stiff->addHeatFluxDVSens(elemIndex, scale, pt, X, Ux, dfdq, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
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
void TACSHeatConduction3D::evalPointQuantitySens(
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

  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;

  if (quantityType == TACS_TEMPERATURE) {
    dfdUt[0] = dfdq[0];
  } else if (quantityType == TACS_HEAT_FLUX) {
    // flux = Kc*grad
    TacsScalar Kc[6];
    stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

    dfdUx[0] = dfdq[0] * Kc[0] + dfdq[1] * Kc[1] + dfdq[2] * Kc[2];
    dfdUx[1] = dfdq[0] * Kc[1] + dfdq[1] * Kc[3] + dfdq[2] * Kc[4];
    dfdUx[2] = dfdq[0] * Kc[2] + dfdq[1] * Kc[4] + dfdq[2] * Kc[5];
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
void TACSHeatConduction3D::getOutputData(
    int elemIndex, const double time, ElementType etype, int write_flag,
    const double pt[], const TacsScalar X[], const TacsScalar Ut[],
    const TacsScalar Ux[], int ld_data, TacsScalar *data) {
  if (etype == TACS_SCALAR_3D_ELEMENT) {
    if (write_flag & TACS_OUTPUT_NODES) {
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      data[0] = Ut[0];
      data += 1;
    }

    TacsScalar grad[3];
    grad[0] = Ux[0];
    grad[1] = Ux[1];
    grad[2] = Ux[2];

    if (write_flag & TACS_OUTPUT_STRAINS) {
      data[0] = grad[0];
      data[1] = grad[1];
      data[2] = grad[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_STRESSES) {
      TacsScalar flux[3];
      stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
      data[0] = flux[0];
      data[1] = flux[1];
      data[2] = flux[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = stiff->evalDensity(elemIndex, pt, X);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}
