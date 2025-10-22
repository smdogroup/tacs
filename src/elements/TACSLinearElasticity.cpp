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

#include "TACSLinearElasticity.h"

TACSLinearElasticity2D::TACSLinearElasticity2D(
    TACSPlaneStressConstitutive *_stiff, ElementStrainType _strain_type) {
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
}

TACSLinearElasticity2D::~TACSLinearElasticity2D() { stiff->decref(); }

// 0;   1;    2;   3;   4; 5;   6;    7;   8;   9;
// u; u,t; u,tt; u,x; u,y; v; v,t; v,tt; v,x; v,y;

const int TACSLinearElasticity2D::linear_Jac_pairs[] = {
    2, 2, 7, 7, 3, 3, 3, 4, 3, 8, 3, 9, 4, 3, 4, 4, 4, 8,
    4, 9, 8, 3, 8, 4, 8, 8, 8, 9, 9, 3, 9, 4, 9, 8, 9, 9};

int TACSLinearElasticity2D::getNumParameters() { return 2; }

int TACSLinearElasticity2D::getVarsPerNode() { return 2; }

int TACSLinearElasticity2D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearElasticity2D::getDesignVarNums(int elemIndex, int dvLen,
                                             int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearElasticity2D::setDesignVars(int elemIndex, int dvLen,
                                          const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearElasticity2D::getDesignVars(int elemIndex, int dvLen,
                                          TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearElasticity2D::getDesignVarRange(int elemIndex, int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearElasticity2D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho * Ut[2];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho * Ut[5];

  if (strain_type == TACS_LINEAR_STRAIN) {
    TacsScalar e[3];
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];

    // Evaluate the stress
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    DUx[0] = s[0];  // u,x
    DUx[1] = s[2];  // u,y

    DUx[2] = s[2];  // v,x
    DUx[3] = s[1];  // v,y
  } else {
    TacsScalar e[3];
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);

    // Evaluate the stress
    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    // Set the coefficients for the weak form
    // Coef = S + Ux*S
    DUx[0] = Ux[1] * s[2] + s[0] * (Ux[0] + 1.0);
    DUx[1] = Ux[1] * s[1] + s[2] * (Ux[0] + 1.0);
    DUx[2] = Ux[2] * s[0] + s[2] * (Ux[3] + 1.0);
    DUx[3] = Ux[2] * s[2] + s[1] * (Ux[3] + 1.0);
  }
}

/*
  Add the design variable derivative of the product of the adjoint
  vector with the weak form of the residual
*/
void TACSLinearElasticity2D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar *fdvSens) {
  // Evaluate the density
  TacsScalar rho_coef = scale * (Ut[2] * Psi[0] + Ut[5] * Psi[1]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, fdvSens);

  if (strain_type == TACS_LINEAR_STRAIN) {
    TacsScalar e[3];
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];

    TacsScalar phi[3];
    phi[0] = Psix[0];
    phi[1] = Psix[3];
    phi[2] = Psix[1] + Psix[2];
    stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, fdvSens);
  } else {
    TacsScalar e[3];
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);

    TacsScalar phi[3];
    phi[0] = Psix[0] * (Ux[0] + 1.0) + Psix[2] * Ux[2];
    phi[1] = Psix[1] * Ux[1] + Psix[3] * (Ux[3] + 1.0);
    phi[2] = Psix[0] * Ux[1] + Psix[1] * (Ux[0] + 1.0) +
             Psix[2] * (Ux[3] + 1.0) + Psix[3] * Ux[2];

    stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, fdvSens);
  }
}

void TACSLinearElasticity2D::evalWeakAdjXptSensProduct(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
    TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
    TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;
  dfdUx[0] = dfdUx[1] = dfdUx[2] = dfdUx[3] = 0.0;
  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = dfdPsix[3] = 0.0;
  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;

  TacsScalar e[3], phi[3];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];

    phi[0] = Psix[0];
    phi[1] = Psix[3];
    phi[2] = Psix[1] + Psix[2];
  } else {
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
    e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
    e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);

    phi[0] = Psix[0] * (Ux[0] + 1.0) + Psix[2] * Ux[2];
    phi[1] = Psix[1] * Ux[1] + Psix[3] * (Ux[3] + 1.0);
    phi[2] = Psix[0] * Ux[1] + Psix[1] * (Ux[0] + 1.0) +
             Psix[2] * (Ux[3] + 1.0) + Psix[3] * Ux[2];
  }

  // Compute the material density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  TacsScalar s[3], a[3];
  stiff->evalStress(elemIndex, pt, X, e, s);
  stiff->evalStress(elemIndex, pt, X, phi, a);

  *product = (rho * (Psi[0] * Ut[2] + Psi[1] * Ut[5]) + a[0] * e[0] +
              a[1] * e[1] + a[2] * e[2]);

  if (strain_type == TACS_LINEAR_STRAIN) {
    dfdPsix[0] = s[0];
    dfdPsix[3] = s[1];
    dfdPsix[1] = s[2];
    dfdPsix[2] = s[2];

    dfdUx[0] = a[0];
    dfdUx[3] = a[1];
    dfdUx[1] = a[2];
    dfdUx[2] = a[2];
  } else {
    dfdPsix[0] = Ux[1] * s[2] + s[0] * (Ux[0] + 1.0);
    dfdPsix[1] = Ux[1] * s[1] + s[2] * (Ux[0] + 1.0);
    dfdPsix[2] = Ux[2] * s[0] + s[2] * (Ux[3] + 1.0);
    dfdPsix[3] = Ux[2] * s[2] + s[1] * (Ux[3] + 1.0);

    dfdUx[0] =
        Psix[0] * s[0] + Psix[1] * s[2] + Ux[1] * a[2] + a[0] * (Ux[0] + 1.0);
    dfdUx[1] =
        Psix[0] * s[2] + Psix[1] * s[1] + Ux[1] * a[1] + a[2] * (Ux[0] + 1.0);
    dfdUx[2] =
        Psix[2] * s[0] + Psix[3] * s[2] + Ux[2] * a[0] + a[2] * (Ux[3] + 1.0);
    dfdUx[3] =
        Psix[2] * s[2] + Psix[3] * s[1] + Ux[2] * a[2] + a[1] * (Ux[3] + 1.0);
  }
}

void TACSLinearElasticity2D::getWeakMatrixNonzeros(ElementMatrixType matType,
                                                   int elemIndex, int *Jac_nnz,
                                                   const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 18;
    *Jac_pairs = linear_Jac_pairs;
  } else if (matType == TACS_MASS_MATRIX) {
    *Jac_nnz = 2;
    *Jac_pairs = linear_Jac_pairs;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    *Jac_nnz = 16;
    *Jac_pairs = &linear_Jac_pairs[4];
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    *Jac_nnz = 16;
    *Jac_pairs = &linear_Jac_pairs[4];
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSLinearElasticity2D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = 0.0;
    DUt[2] = rho * Ut[2];  // u,tt

    DUt[3] = 0.0;
    DUt[4] = 0.0;
    DUt[5] = rho * Ut[5];  // v,tt

    // Set the acceleration terms
    Jac[0] = rho;
    Jac[1] = rho;

    if (strain_type == TACS_LINEAR_STRAIN) {
      TacsScalar e[3];
      e[0] = Ux[0];
      e[1] = Ux[3];
      e[2] = Ux[1] + Ux[2];

      // Evaluate the stress
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);

      DUx[0] = s[0];  // u,x
      DUx[1] = s[2];  // u,y

      DUx[2] = s[2];  // v,x
      DUx[3] = s[1];  // v,y

      TacsScalar C[6];
      stiff->evalTangentStiffness(elemIndex, pt, X, C);

      // Index:         3            9            4     8
      // s[0] = C[0]*(u,x) + C[1]*(v,y) + C[2]*(u,y + v,x)
      // s[1] = C[1]*(u,x) + C[3]*(v,y) + C[4]*(u,y + v,x)
      // s[2] = C[2]*(u,x) + C[4]*(v,y) + C[5]*(u,y + v,x)

      // i == 1 (s[0])
      Jac[2] = C[0];  // j == 3
      Jac[3] = C[2];  // j == 4
      Jac[4] = C[2];  // j == 8
      Jac[5] = C[1];  // j == 9

      // i == 2 (s[2])
      Jac[6] = C[2];  // j == 3
      Jac[7] = C[5];  // j == 4
      Jac[8] = C[5];  // j == 8
      Jac[9] = C[4];  // j == 9

      // i == 4 (s[2])
      Jac[10] = C[2];  // j == 3
      Jac[11] = C[5];  // j == 4
      Jac[12] = C[5];  // j == 8
      Jac[13] = C[4];  // j == 9

      // i == 5 (s[1])
      Jac[14] = C[1];  // j == 3
      Jac[15] = C[4];  // j == 4
      Jac[16] = C[4];  // j == 8
      Jac[17] = C[3];  // j == 9
    } else {
      TacsScalar e[3];
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[2] * Ux[2]);
      e[1] = Ux[3] + 0.5 * (Ux[1] * Ux[1] + Ux[3] * Ux[3]);
      e[2] = Ux[1] + Ux[2] + (Ux[0] * Ux[1] + Ux[2] * Ux[3]);

      // Evaluate the stress
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Set the coefficients for the weak form
      // Coef = S + Ux*S
      DUx[0] = Ux[1] * s[2] + s[0] * (Ux[0] + 1.0);
      DUx[1] = Ux[1] * s[1] + s[2] * (Ux[0] + 1.0);
      DUx[2] = Ux[2] * s[0] + s[2] * (Ux[3] + 1.0);
      DUx[3] = Ux[2] * s[2] + s[1] * (Ux[3] + 1.0);

      TacsScalar C[6];
      stiff->evalTangentStiffness(elemIndex, pt, X, C);

      // Compute the derivative of the stress w.r.t. the components of Ux
      TacsScalar dsdUx[3 * 4];
      dsdUx[0] = C[0] * (Ux[0] + 1.0) + C[2] * Ux[1];
      dsdUx[1] = C[1] * Ux[1] + C[2] * (Ux[0] + 1.0);
      dsdUx[2] = C[0] * Ux[2] + C[2] * (Ux[3] + 1.0);
      dsdUx[3] = C[1] * (Ux[3] + 1.0) + C[2] * Ux[2];

      dsdUx[4] = C[1] * (Ux[0] + 1.0) + C[4] * Ux[1];
      dsdUx[5] = C[3] * Ux[1] + C[4] * (Ux[0] + 1.0);
      dsdUx[6] = C[1] * Ux[2] + C[4] * (Ux[3] + 1.0);
      dsdUx[7] = C[3] * (Ux[3] + 1.0) + C[4] * Ux[2];

      dsdUx[8] = C[2] * (Ux[0] + 1.0) + C[5] * Ux[1];
      dsdUx[9] = C[4] * Ux[1] + C[5] * (Ux[0] + 1.0);
      dsdUx[10] = C[2] * Ux[2] + C[5] * (Ux[3] + 1.0);
      dsdUx[11] = C[4] * (Ux[3] + 1.0) + C[5] * Ux[2];

      Jac[2] = Ux[1] * dsdUx[8] + dsdUx[0] * (Ux[0] + 1.0) + s[0];
      Jac[3] = Ux[1] * dsdUx[9] + dsdUx[1] * (Ux[0] + 1.0) + s[2];
      Jac[4] = Ux[1] * dsdUx[10] + dsdUx[2] * (Ux[0] + 1.0);
      Jac[5] = Ux[1] * dsdUx[11] + dsdUx[3] * (Ux[0] + 1.0);

      Jac[6] = Ux[1] * dsdUx[4] + dsdUx[8] * (Ux[0] + 1.0) + s[2];
      Jac[7] = Ux[1] * dsdUx[5] + dsdUx[9] * (Ux[0] + 1.0) + s[1];
      Jac[8] = Ux[1] * dsdUx[6] + dsdUx[10] * (Ux[0] + 1.0);
      Jac[9] = Ux[1] * dsdUx[7] + dsdUx[11] * (Ux[0] + 1.0);

      Jac[10] = Ux[2] * dsdUx[0] + dsdUx[8] * (Ux[3] + 1.0);
      Jac[11] = Ux[2] * dsdUx[1] + dsdUx[9] * (Ux[3] + 1.0);
      Jac[12] = Ux[2] * dsdUx[2] + dsdUx[10] * (Ux[3] + 1.0) + s[0];
      Jac[13] = Ux[2] * dsdUx[3] + dsdUx[11] * (Ux[3] + 1.0) + s[2];

      Jac[14] = Ux[2] * dsdUx[8] + dsdUx[4] * (Ux[3] + 1.0);
      Jac[15] = Ux[2] * dsdUx[9] + dsdUx[5] * (Ux[3] + 1.0);
      Jac[16] = Ux[2] * dsdUx[10] + dsdUx[6] * (Ux[3] + 1.0) + s[2];
      Jac[17] = Ux[2] * dsdUx[11] + dsdUx[7] * (Ux[3] + 1.0) + s[1];
    }
  } else if (matType == TACS_MASS_MATRIX) {
    // Evaluate the density
    TacsScalar rho = stiff->evalMassMatrixDensity(elemIndex, pt, X);

    // Set the acceleration terms
    Jac[0] = rho;
    Jac[1] = rho;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    TacsScalar C[6];
    stiff->evalTangentStiffness(elemIndex, pt, X, C);

    Jac[0] = C[0];  // j == 3
    Jac[1] = C[2];  // j == 4
    Jac[2] = C[2];  // j == 8
    Jac[3] = C[1];  // j == 9

    // i == 2 (s[2])
    Jac[4] = C[2];  // j == 3
    Jac[5] = C[5];  // j == 4
    Jac[6] = C[5];  // j == 8
    Jac[7] = C[4];  // j == 9

    // i == 4 (s[2])
    Jac[8] = C[2];   // j == 3
    Jac[9] = C[5];   // j == 4
    Jac[10] = C[5];  // j == 8
    Jac[11] = C[4];  // j == 9

    // i == 5 (s[1])
    Jac[12] = C[1];  // j == 3
    Jac[13] = C[4];  // j == 4
    Jac[14] = C[4];  // j == 8
    Jac[15] = C[3];  // j == 9
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    // Compute the tangent stiffness
    TacsScalar C[6];
    stiff->evalGeometricTangentStiffness(elemIndex, pt, X, C);

    // Compute the tangent strain and stress
    TacsScalar e[3], s[3];
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];

    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2];
    s[1] = C[1] * e[0] + C[3] * e[1] + C[4] * e[2];
    s[2] = C[2] * e[0] + C[4] * e[1] + C[5] * e[2];

    Jac[0] = s[0];
    Jac[1] = s[2];
    Jac[2] = 0.0;
    Jac[3] = 0.0;

    Jac[4] = s[2];
    Jac[5] = s[1];
    Jac[6] = 0.0;
    Jac[7] = 0.0;

    Jac[8] = 0.0;
    Jac[9] = 0.0;
    Jac[10] = s[0];
    Jac[11] = s[2];

    Jac[12] = 0.0;
    Jac[13] = 0.0;
    Jac[14] = s[2];
    Jac[15] = s[1];
  }
}

void TACSLinearElasticity2D::addWeakMatDVSens(
    ElementMatrixType matType, int elemIndex, const double time,
    TacsScalar scale, int n, const double pt[], const TacsScalar X[],
    const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
    const TacsScalar Psi[], const TacsScalar Psix[], const TacsScalar Phi[],
    const TacsScalar Phix[], int dvLen, TacsScalar *dfdx) {
  if (matType == TACS_MASS_MATRIX) {
    TacsScalar rho_coef = Psi[0] * Phi[0] + Psi[3] * Phi[3];

    stiff->addMassMatrixDensityDVSens(elemIndex, scale * rho_coef, pt, X, dvLen,
                                      dfdx);
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    TacsScalar ePsi[3];
    ePsi[0] = Psix[0];
    ePsi[1] = Psix[3];
    ePsi[2] = Psix[1] + Psix[2];

    TacsScalar ePhi[3];
    ePhi[0] = Phix[0];
    ePhi[1] = Phix[3];
    ePhi[2] = Phix[1] + Phix[2];

    stiff->addStressDVSens(elemIndex, scale, pt, X, ePsi, ePhi, dvLen, dfdx);
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    // Compute the tangent strain and stress
    TacsScalar e[3];
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];

    TacsScalar psi[3];
    psi[0] = Phix[0] * Psix[0] + Phix[2] * Psix[2];
    psi[1] = Phix[1] * Psix[1] + Phix[3] * Psix[3];
    psi[2] = Phix[0] * Psix[1] + Phix[1] * Psix[0] + Phix[2] * Psix[3] +
             Phix[3] * Psix[2];

    stiff->addGeometricTangentStressDVSens(elemIndex, scale, pt, X, e, psi,
                                           dvLen, dfdx);
  }
}

void TACSLinearElasticity2D::evalWeakMatSVSens(
    ElementMatrixType matType, int elemIndex, const double time,
    TacsScalar scale, int n, const double pt[], const TacsScalar X[],
    const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
    const TacsScalar Psi[], const TacsScalar Psix[], const TacsScalar Phi[],
    const TacsScalar Phix[], TacsScalar dfdU[], TacsScalar dfdUx[]) {
  dfdU[0] = dfdU[1] = 0.0;
  dfdU[2] = dfdU[3] = 0.0;
  dfdU[4] = dfdU[5] = 0.0;
  dfdUx[0] = dfdUx[1] = 0.0;
  dfdUx[2] = dfdUx[3] = 0.0;

  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TacsScalar e[3];
    e[0] = Phix[0] * Psix[0] + Phix[2] * Psix[2];
    e[1] = Phix[1] * Psix[1] + Phix[3] * Psix[3];
    e[2] = Phix[0] * Psix[1] + Phix[1] * Psix[0] + Phix[2] * Psix[3] +
           Phix[3] * Psix[2];

    // Compute the tangent stiffness
    TacsScalar C[6];
    stiff->evalGeometricTangentStiffness(elemIndex, pt, X, C);

    // Compute the tangent strain and stress
    TacsScalar s[3];
    s[0] = C[0] * e[0] + C[1] * e[1] + C[2] * e[2];
    s[1] = C[1] * e[0] + C[3] * e[1] + C[4] * e[2];
    s[2] = C[2] * e[0] + C[4] * e[1] + C[5] * e[2];

    dfdUx[0] = scale * s[0];
    dfdUx[3] = scale * s[1];

    dfdUx[1] = scale * s[2];
    dfdUx[2] = scale * s[2];
  }
}
/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearElasticity2D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
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

      *quantity = stiff->evalFailure(elemIndex, pt, X, e);
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }
    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
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
void TACSLinearElasticity2D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX) {
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

    stiff->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, e, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
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

    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);
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
void TACSLinearElasticity2D::evalPointQuantitySens(
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

  dfdUx[0] = dfdUx[1] = 0.0;
  dfdUx[2] = dfdUx[3] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX) {
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

    TacsScalar sens[3];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = dfdq[0] * sens[0];
      dfdUx[3] = dfdq[0] * sens[1];

      dfdUx[1] = dfdq[0] * sens[2];
      dfdUx[2] = dfdq[0] * sens[2];
    } else {
      dfdUx[0] = dfdq[0] * (sens[0] * (1.0 + Ux[0]) + sens[2] * Ux[1]);
      dfdUx[1] = dfdq[0] * (sens[1] * Ux[1] + sens[2] * (1.0 + Ux[0]));

      dfdUx[2] = dfdq[0] * (sens[0] * Ux[2] + sens[2] * (1.0 + Ux[3]));
      dfdUx[3] = dfdq[0] * (sens[1] * (1.0 + Ux[3]) + sens[2] * Ux[2]);
    }
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
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

    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN) {
      dfdUx[0] = 2.0 * dfdq[0] * s[0];
      dfdUx[3] = 2.0 * dfdq[0] * s[1];

      dfdUx[1] = 2.0 * dfdq[0] * s[2];
      dfdUx[2] = 2.0 * dfdq[0] * s[2];
    } else {
      dfdUx[0] = 2.0 * dfdq[0] * (s[0] * (1.0 + Ux[0]) + s[2] * Ux[1]);
      dfdUx[1] = 2.0 * dfdq[0] * (s[1] * Ux[1] + s[2] * (1.0 + Ux[0]));

      dfdUx[2] = 2.0 * dfdq[0] * (s[0] * Ux[2] + s[2] * (1.0 + Ux[3]));
      dfdUx[3] = 2.0 * dfdq[0] * (s[1] * (1.0 + Ux[3]) + s[2] * Ux[2]);
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
void TACSLinearElasticity2D::getOutputData(
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

TACSLinearElasticity3D::TACSLinearElasticity3D(TACSSolidConstitutive *_stiff,
                                               ElementStrainType _strain_type) {
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
}

TACSLinearElasticity3D::~TACSLinearElasticity3D() { stiff->decref(); }

// 0;   1;    2;   3;   4;   5;
// u; u,t; u,tt; u,x; u,y; u,z;

// 6;   7;    8;   9;  10;  11;
// v; v,t; v,tt; v,x; v,y; v,z;

// 12;  13;   14;  15;  16;  17;
//  w; w,t; w,tt; w,x; w,y; w,z;

const int TACSLinearElasticity3D::linear_Jac_pairs[] = {
    2,  2,  8,  8,  14, 14, 3,  3,  3,  4,  3,  5,  3,  9,  3,  10, 3,  11, 3,
    15, 3,  16, 3,  17, 4,  3,  4,  4,  4,  5,  4,  9,  4,  10, 4,  11, 4,  15,
    4,  16, 4,  17, 5,  3,  5,  4,  5,  5,  5,  9,  5,  10, 5,  11, 5,  15, 5,
    16, 5,  17, 9,  3,  9,  4,  9,  5,  9,  9,  9,  10, 9,  11, 9,  15, 9,  16,
    9,  17, 10, 3,  10, 4,  10, 5,  10, 9,  10, 10, 10, 11, 10, 15, 10, 16, 10,
    17, 11, 3,  11, 4,  11, 5,  11, 9,  11, 10, 11, 11, 11, 15, 11, 16, 11, 17,
    15, 3,  15, 4,  15, 5,  15, 9,  15, 10, 15, 11, 15, 15, 15, 16, 15, 17, 16,
    3,  16, 4,  16, 5,  16, 9,  16, 10, 16, 11, 16, 15, 16, 16, 16, 17, 17, 3,
    17, 4,  17, 5,  17, 9,  17, 10, 17, 11, 17, 15, 17, 16, 17, 17};

int TACSLinearElasticity3D::getNumParameters() { return 3; }

int TACSLinearElasticity3D::getVarsPerNode() { return 3; }

int TACSLinearElasticity3D::getDesignVarsPerNode() {
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearElasticity3D::getDesignVarNums(int elemIndex, int dvLen,
                                             int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearElasticity3D::setDesignVars(int elemIndex, int dvLen,
                                          const TacsScalar dvs[]) {
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearElasticity3D::getDesignVars(int elemIndex, int dvLen,
                                          TacsScalar dvs[]) {
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearElasticity3D::getDesignVarRange(int elemIndex, int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[]) {
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearElasticity3D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho * Ut[2];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho * Ut[5];

  DUt[6] = 0.0;
  DUt[7] = 0.0;
  DUt[8] = rho * Ut[8];

  if (strain_type == TACS_LINEAR_STRAIN) {
    TacsScalar e[6];
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];

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
  } else {
    TacsScalar e[6];
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
    e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
    e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

    e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
    e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
    e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);

    // Evaluate the stress
    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    // Set the coefficients for the weak form
    // Coef = (I + Ux)*S
    DUx[0] = Ux[1] * s[5] + Ux[2] * s[4] + s[0] * (Ux[0] + 1.0);
    DUx[1] = Ux[1] * s[1] + Ux[2] * s[3] + s[5] * (Ux[0] + 1.0);
    DUx[2] = Ux[1] * s[3] + Ux[2] * s[2] + s[4] * (Ux[0] + 1.0);

    DUx[3] = Ux[3] * s[0] + Ux[5] * s[4] + s[5] * (Ux[4] + 1.0);
    DUx[4] = Ux[3] * s[5] + Ux[5] * s[3] + s[1] * (Ux[4] + 1.0);
    DUx[5] = Ux[3] * s[4] + Ux[5] * s[2] + s[3] * (Ux[4] + 1.0);

    DUx[6] = Ux[6] * s[0] + Ux[7] * s[5] + s[4] * (Ux[8] + 1.0);
    DUx[7] = Ux[6] * s[5] + Ux[7] * s[1] + s[3] * (Ux[8] + 1.0);
    DUx[8] = Ux[6] * s[4] + Ux[7] * s[3] + s[2] * (Ux[8] + 1.0);
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSLinearElasticity3D::addWeakAdjProduct(
    int elemIndex, const double time, TacsScalar scale, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar Psi[],
    const TacsScalar Psix[], int dvLen, TacsScalar *fdvSens) {
  // Evaluate the density
  TacsScalar rho_coef =
      scale * (Ut[2] * Psi[0] + Ut[5] * Psi[1] + Ut[8] * Psi[2]);
  stiff->addDensityDVSens(elemIndex, rho_coef, pt, X, dvLen, fdvSens);

  if (strain_type == TACS_LINEAR_STRAIN) {
    TacsScalar e[6];
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];

    TacsScalar phi[6];
    phi[0] = Psix[0];
    phi[1] = Psix[4];
    phi[2] = Psix[8];
    phi[3] = Psix[5] + Psix[7];
    phi[4] = Psix[2] + Psix[6];
    phi[5] = Psix[1] + Psix[3];

    stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, fdvSens);
  } else {
    TacsScalar e[6];
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
    e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
    e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

    e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
    e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
    e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);

    // Compute the contribution to the adjoint
    TacsScalar phi[6];
    phi[0] = Psix[0] * (Ux[0] + 1.0) + Psix[3] * Ux[3] + Psix[6] * Ux[6];
    phi[1] = Psix[1] * Ux[1] + Psix[4] * (Ux[4] + 1.0) + Psix[7] * Ux[7];
    phi[2] = Psix[2] * Ux[2] + Psix[5] * Ux[5] + Psix[8] * (Ux[8] + 1.0);

    phi[3] = Psix[1] * Ux[2] + Psix[2] * Ux[1] + Psix[4] * Ux[5] +
             Psix[5] * (Ux[4] + 1.0) + Psix[7] * (Ux[8] + 1.0) +
             Psix[8] * Ux[7];
    phi[4] = Psix[0] * Ux[2] + Psix[2] * (Ux[0] + 1.0) + Psix[3] * Ux[5] +
             Psix[5] * Ux[3] + Psix[6] * (Ux[8] + 1.0) + Psix[8] * Ux[6];
    phi[5] = Psix[0] * Ux[1] + Psix[1] * (Ux[0] + 1.0) +
             Psix[3] * (Ux[4] + 1.0) + Psix[4] * Ux[3] + Psix[6] * Ux[7] +
             Psix[7] * Ux[6];

    stiff->addStressDVSens(elemIndex, scale, pt, X, e, phi, dvLen, fdvSens);
  }
}

void TACSLinearElasticity3D::evalWeakAdjXptSensProduct(
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

  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = 0.0;
  dfdPsix[3] = dfdPsix[4] = dfdPsix[5] = 0.0;
  dfdPsix[6] = dfdPsix[7] = dfdPsix[8] = 0.0;

  TacsScalar e[6], phi[6];
  if (strain_type == TACS_LINEAR_STRAIN) {
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];
    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];

    phi[0] = Psix[0];
    phi[1] = Psix[4];
    phi[2] = Psix[8];

    phi[3] = Psix[5] + Psix[7];
    phi[4] = Psix[2] + Psix[6];
    phi[5] = Psix[1] + Psix[3];
  } else {
    e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
    e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
    e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

    e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
    e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
    e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);

    phi[0] = Psix[0] * (Ux[0] + 1.0) + Psix[3] * Ux[3] + Psix[6] * Ux[6];
    phi[1] = Psix[1] * Ux[1] + Psix[4] * (Ux[4] + 1.0) + Psix[7] * Ux[7];
    phi[2] = Psix[2] * Ux[2] + Psix[5] * Ux[5] + Psix[8] * (Ux[8] + 1.0);

    phi[3] = Psix[1] * Ux[2] + Psix[2] * Ux[1] + Psix[4] * Ux[5] +
             Psix[5] * (Ux[4] + 1.0) + Psix[7] * (Ux[8] + 1.0) +
             Psix[8] * Ux[7];
    phi[4] = Psix[0] * Ux[2] + Psix[2] * (Ux[0] + 1.0) + Psix[3] * Ux[5] +
             Psix[5] * Ux[3] + Psix[6] * (Ux[8] + 1.0) + Psix[8] * Ux[6];
    phi[5] = Psix[0] * Ux[1] + Psix[1] * (Ux[0] + 1.0) +
             Psix[3] * (Ux[4] + 1.0) + Psix[4] * Ux[3] + Psix[6] * Ux[7] +
             Psix[7] * Ux[6];
  }

  // Compute the material density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  TacsScalar s[6], a[6];
  stiff->evalStress(elemIndex, pt, X, e, s);
  stiff->evalStress(elemIndex, pt, X, phi, a);

  *product =
      (rho * (Psi[0] * Ut[2] + Psi[1] * Ut[5] + Psi[2] * Ut[8]) + a[0] * e[0] +
       a[1] * e[1] + a[2] * e[2] + a[3] * e[3] + a[4] * e[4] + a[5] * e[5]);

  if (strain_type == TACS_LINEAR_STRAIN) {
    dfdUx[0] = a[0];
    dfdUx[4] = a[1];
    dfdUx[8] = a[2];
    dfdUx[5] = a[3];
    dfdUx[7] = a[3];
    dfdUx[2] = a[4];
    dfdUx[6] = a[4];
    dfdUx[1] = a[5];
    dfdUx[3] = a[5];

    dfdPsix[0] = s[0];
    dfdPsix[4] = s[1];
    dfdPsix[8] = s[2];
    dfdPsix[5] = s[3];
    dfdPsix[7] = s[3];
    dfdPsix[2] = s[4];
    dfdPsix[6] = s[4];
    dfdPsix[1] = s[5];
    dfdPsix[3] = s[5];
  } else {
    dfdPsix[0] = Ux[1] * s[5] + Ux[2] * s[4] + s[0] * (Ux[0] + 1.0);
    dfdPsix[1] = Ux[1] * s[1] + Ux[2] * s[3] + s[5] * (Ux[0] + 1.0);
    dfdPsix[2] = Ux[1] * s[3] + Ux[2] * s[2] + s[4] * (Ux[0] + 1.0);

    dfdPsix[3] = Ux[3] * s[0] + Ux[5] * s[4] + s[5] * (Ux[4] + 1.0);
    dfdPsix[4] = Ux[3] * s[5] + Ux[5] * s[3] + s[1] * (Ux[4] + 1.0);
    dfdPsix[5] = Ux[3] * s[4] + Ux[5] * s[2] + s[3] * (Ux[4] + 1.0);

    dfdPsix[6] = Ux[6] * s[0] + Ux[7] * s[5] + s[4] * (Ux[8] + 1.0);
    dfdPsix[7] = Ux[6] * s[5] + Ux[7] * s[1] + s[3] * (Ux[8] + 1.0);
    dfdPsix[8] = Ux[6] * s[4] + Ux[7] * s[3] + s[2] * (Ux[8] + 1.0);

    dfdUx[0] = Psix[0] * s[0] + Psix[1] * s[5] + Psix[2] * s[4] + Ux[1] * a[5] +
               Ux[2] * a[4] + a[0] * (Ux[0] + 1.0);
    dfdUx[1] = Psix[0] * s[5] + Psix[1] * s[1] + Psix[2] * s[3] + Ux[1] * a[1] +
               Ux[2] * a[3] + a[5] * (Ux[0] + 1.0);
    dfdUx[2] = Psix[0] * s[4] + Psix[1] * s[3] + Psix[2] * s[2] + Ux[1] * a[3] +
               Ux[2] * a[2] + a[4] * (Ux[0] + 1.0);

    dfdUx[3] = Psix[3] * s[0] + Psix[4] * s[5] + Psix[5] * s[4] + Ux[3] * a[0] +
               Ux[5] * a[4] + a[5] * (Ux[4] + 1.0);
    dfdUx[4] = Psix[3] * s[5] + Psix[4] * s[1] + Psix[5] * s[3] + Ux[3] * a[5] +
               Ux[5] * a[3] + a[1] * (Ux[4] + 1.0);
    dfdUx[5] = Psix[3] * s[4] + Psix[4] * s[3] + Psix[5] * s[2] + Ux[3] * a[4] +
               Ux[5] * a[2] + a[3] * (Ux[4] + 1.0);

    dfdUx[6] = Psix[6] * s[0] + Psix[7] * s[5] + Psix[8] * s[4] + Ux[6] * a[0] +
               Ux[7] * a[5] + a[4] * (Ux[8] + 1.0);
    dfdUx[7] = Psix[6] * s[5] + Psix[7] * s[1] + Psix[8] * s[3] + Ux[6] * a[5] +
               Ux[7] * a[1] + a[3] * (Ux[8] + 1.0);
    dfdUx[8] = Psix[6] * s[4] + Psix[7] * s[3] + Psix[8] * s[2] + Ux[6] * a[4] +
               Ux[7] * a[3] + a[2] * (Ux[8] + 1.0);
  }
}

void TACSLinearElasticity3D::getWeakMatrixNonzeros(ElementMatrixType matType,
                                                   int elemIndex, int *Jac_nnz,
                                                   const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 84;
    *Jac_pairs = linear_Jac_pairs;
  } else if (matType == TACS_MASS_MATRIX) {
    *Jac_nnz = 3;
    *Jac_pairs = linear_Jac_pairs;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    *Jac_nnz = 81;
    *Jac_pairs = &linear_Jac_pairs[6];
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    *Jac_nnz = 81;
    *Jac_pairs = &linear_Jac_pairs[6];
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSLinearElasticity3D::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    // Evaluate the density
    TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

    DUt[0] = 0.0;
    DUt[1] = 0.0;
    DUt[2] = rho * Ut[2];

    DUt[3] = 0.0;
    DUt[4] = 0.0;
    DUt[5] = rho * Ut[5];

    DUt[6] = 0.0;
    DUt[7] = 0.0;
    DUt[8] = rho * Ut[8];

    // Acceleration terms
    Jac[0] = rho;
    Jac[1] = rho;
    Jac[2] = rho;

    if (strain_type == TACS_LINEAR_STRAIN) {
      TacsScalar e[6];
      e[0] = Ux[0];
      e[1] = Ux[4];
      e[2] = Ux[8];

      e[3] = Ux[5] + Ux[7];
      e[4] = Ux[2] + Ux[6];
      e[5] = Ux[1] + Ux[3];

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

      TacsScalar C[21];
      stiff->evalTangentStiffness(elemIndex, pt, X, C);

      // u,x u,y u,z v,x v,y v,z w,x w,y w,z
      // e[0] = Ux[0]; e[1] = Ux[4]; e[2] = Ux[8];
      // e[3] = Ux[5] + Ux[7]; e[4] = Ux[2] + Ux[6]; e[5] = Ux[1] + Ux[3];

      // s =  [s0 s5 s4]
      //      [s5 s1 s3]
      //      [s4 s3 s2]

      // Index:
      // s[0] = C[0]*(u,x) + C[1]*(v,y) + C[2]*(w,z) + C[3]*(v,z + w,y)
      //                   + C[4]*(u,z + w,x) + C[5]*(u,y + v,x)
      // s[1] = C[1]*(u,x) + C[6]*(v,y) + C[7]*(w,z) + C[8]*(v,z + w,y)
      //                   + C[9]*(u,z + w,x) + C[10]*(u,y + v,x)
      // s[2] = C[2]*(u,x) + C[7]*(v,y) + C[11]*(w,z) + C[12]*(v,z + w,y)
      //                   + C[13]*(u,z + w,x) + C[14]*(u,y + v,x)
      // s[3] = C[3]*(u,x) + C[8]*(v,y) + C[12]*(w,z) + C[15]*(v,z + w,y)
      //                   + C[16]*(u,z + w,x) + C[17]*(u,y + v,x)
      // s[4] = C[4]*(u,x) + C[9]*(v,y) + C[13]*(w,z) + C[16]*(v,z + w,y)
      //                   + C[18]*(u,z + w,x) + C[19]*(u,y + v,x)
      // s[5] = C[5]*(u,x) + C[10]*(v,y) + C[14]*(w,z) + C[17]*(v,z + w,y)
      //                   + C[19]*(u,z + w,x) + C[20]*(u,y + v,x)

      // 0;   1;    2;   3;   4;   5;
      // u; u,t; u,tt; u,x; u,y; u,z;

      // 6;   7;    8;   9;  10;  11;
      // v; v,t; v,tt; v,x; v,y; v,z;

      // 12;  13;   14;  15;  16;  17;
      //  w; w,t; w,tt; w,x; w,y; w,z;

      // s[0]
      Jac[3] = C[0];   // u,x 3
      Jac[4] = C[5];   // u,y 4
      Jac[5] = C[4];   // u,z 5
      Jac[6] = C[5];   // v,x 9
      Jac[7] = C[1];   // v,y 10
      Jac[8] = C[3];   // v,z 11
      Jac[9] = C[4];   // w,x 15
      Jac[10] = C[3];  // w,y 16
      Jac[11] = C[2];  // w,z 17

      // s[5]
      Jac[12] = C[5];
      Jac[13] = C[20];
      Jac[14] = C[19];
      Jac[15] = C[20];
      Jac[16] = C[10];
      Jac[17] = C[17];
      Jac[18] = C[19];
      Jac[19] = C[17];
      Jac[20] = C[14];

      // s[4]
      Jac[21] = C[4];
      Jac[22] = C[19];
      Jac[23] = C[18];
      Jac[24] = C[19];
      Jac[25] = C[9];
      Jac[26] = C[16];
      Jac[27] = C[18];
      Jac[28] = C[16];
      Jac[29] = C[13];

      // s[5]
      Jac[30] = C[5];
      Jac[31] = C[20];
      Jac[32] = C[19];
      Jac[33] = C[20];
      Jac[34] = C[10];
      Jac[35] = C[17];
      Jac[36] = C[19];
      Jac[37] = C[17];
      Jac[38] = C[14];

      // s[1]
      Jac[39] = C[1];
      Jac[40] = C[10];
      Jac[41] = C[9];
      Jac[42] = C[10];
      Jac[43] = C[6];
      Jac[44] = C[8];
      Jac[45] = C[9];
      Jac[46] = C[8];
      Jac[47] = C[7];

      // s[3]
      Jac[48] = C[3];
      Jac[49] = C[17];
      Jac[50] = C[16];
      Jac[51] = C[17];
      Jac[52] = C[8];
      Jac[53] = C[15];
      Jac[54] = C[16];
      Jac[55] = C[15];
      Jac[56] = C[12];

      // s[4]
      Jac[57] = C[4];
      Jac[58] = C[19];
      Jac[59] = C[18];
      Jac[60] = C[19];
      Jac[61] = C[9];
      Jac[62] = C[16];
      Jac[63] = C[18];
      Jac[64] = C[16];
      Jac[65] = C[13];

      // s[3]
      Jac[66] = C[3];
      Jac[67] = C[17];
      Jac[68] = C[16];
      Jac[69] = C[17];
      Jac[70] = C[8];
      Jac[71] = C[15];
      Jac[72] = C[16];
      Jac[73] = C[15];
      Jac[74] = C[12];

      // s[2]
      Jac[75] = C[2];
      Jac[76] = C[14];
      Jac[77] = C[13];
      Jac[78] = C[14];
      Jac[79] = C[7];
      Jac[80] = C[12];
      Jac[81] = C[13];
      Jac[82] = C[12];
      Jac[83] = C[11];
    } else {
      TacsScalar e[6];
      e[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
      e[1] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
      e[2] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);

      e[3] = Ux[5] + Ux[7] + (Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
      e[4] = Ux[2] + Ux[6] + (Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
      e[5] = Ux[1] + Ux[3] + (Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);

      // Evaluate the stress
      TacsScalar s[6];
      stiff->evalStress(elemIndex, pt, X, e, s);

      // Set the coefficients for the weak form
      // Coef = (I + Ux)*S
      DUx[0] = Ux[1] * s[5] + Ux[2] * s[4] + s[0] * (Ux[0] + 1.0);
      DUx[1] = Ux[1] * s[1] + Ux[2] * s[3] + s[5] * (Ux[0] + 1.0);
      DUx[2] = Ux[1] * s[3] + Ux[2] * s[2] + s[4] * (Ux[0] + 1.0);

      DUx[3] = Ux[3] * s[0] + Ux[5] * s[4] + s[5] * (Ux[4] + 1.0);
      DUx[4] = Ux[3] * s[5] + Ux[5] * s[3] + s[1] * (Ux[4] + 1.0);
      DUx[5] = Ux[3] * s[4] + Ux[5] * s[2] + s[3] * (Ux[4] + 1.0);

      DUx[6] = Ux[6] * s[0] + Ux[7] * s[5] + s[4] * (Ux[8] + 1.0);
      DUx[7] = Ux[6] * s[5] + Ux[7] * s[1] + s[3] * (Ux[8] + 1.0);
      DUx[8] = Ux[6] * s[4] + Ux[7] * s[3] + s[2] * (Ux[8] + 1.0);

      TacsScalar C[21];
      stiff->evalTangentStiffness(elemIndex, pt, X, C);

      // Compute the derivative of each stress component w.r.t. Ux
      TacsScalar dsdUx[6 * 9];
      dsdUx[0] = C[0] * (Ux[0] + 1.0) + C[4] * Ux[2] + C[5] * Ux[1];
      dsdUx[1] = C[1] * Ux[1] + C[3] * Ux[2] + C[5] * (Ux[0] + 1.0);
      dsdUx[2] = C[2] * Ux[2] + C[3] * Ux[1] + C[4] * (Ux[0] + 1.0);
      dsdUx[3] = C[0] * Ux[3] + C[4] * Ux[5] + C[5] * (Ux[4] + 1.0);
      dsdUx[4] = C[1] * (Ux[4] + 1.0) + C[3] * Ux[5] + C[5] * Ux[3];
      dsdUx[5] = C[2] * Ux[5] + C[3] * (Ux[4] + 1.0) + C[4] * Ux[3];
      dsdUx[6] = C[0] * Ux[6] + C[4] * (Ux[8] + 1.0) + C[5] * Ux[7];
      dsdUx[7] = C[1] * Ux[7] + C[3] * (Ux[8] + 1.0) + C[5] * Ux[6];
      dsdUx[8] = C[2] * (Ux[8] + 1.0) + C[3] * Ux[7] + C[4] * Ux[6];

      dsdUx[9] = C[10] * Ux[1] + C[1] * (Ux[0] + 1.0) + C[9] * Ux[2];
      dsdUx[10] = C[10] * (Ux[0] + 1.0) + C[6] * Ux[1] + C[8] * Ux[2];
      dsdUx[11] = C[7] * Ux[2] + C[8] * Ux[1] + C[9] * (Ux[0] + 1.0);
      dsdUx[12] = C[10] * (Ux[4] + 1.0) + C[1] * Ux[3] + C[9] * Ux[5];
      dsdUx[13] = C[10] * Ux[3] + C[6] * (Ux[4] + 1.0) + C[8] * Ux[5];
      dsdUx[14] = C[7] * Ux[5] + C[8] * (Ux[4] + 1.0) + C[9] * Ux[3];
      dsdUx[15] = C[10] * Ux[7] + C[1] * Ux[6] + C[9] * (Ux[8] + 1.0);
      dsdUx[16] = C[10] * Ux[6] + C[6] * Ux[7] + C[8] * (Ux[8] + 1.0);
      dsdUx[17] = C[7] * (Ux[8] + 1.0) + C[8] * Ux[7] + C[9] * Ux[6];

      dsdUx[18] = C[13] * Ux[2] + C[14] * Ux[1] + C[2] * (Ux[0] + 1.0);
      dsdUx[19] = C[12] * Ux[2] + C[14] * (Ux[0] + 1.0) + C[7] * Ux[1];
      dsdUx[20] = C[11] * Ux[2] + C[12] * Ux[1] + C[13] * (Ux[0] + 1.0);
      dsdUx[21] = C[13] * Ux[5] + C[14] * (Ux[4] + 1.0) + C[2] * Ux[3];
      dsdUx[22] = C[12] * Ux[5] + C[14] * Ux[3] + C[7] * (Ux[4] + 1.0);
      dsdUx[23] = C[11] * Ux[5] + C[12] * (Ux[4] + 1.0) + C[13] * Ux[3];
      dsdUx[24] = C[13] * (Ux[8] + 1.0) + C[14] * Ux[7] + C[2] * Ux[6];
      dsdUx[25] = C[12] * (Ux[8] + 1.0) + C[14] * Ux[6] + C[7] * Ux[7];
      dsdUx[26] = C[11] * (Ux[8] + 1.0) + C[12] * Ux[7] + C[13] * Ux[6];

      dsdUx[27] = C[16] * Ux[2] + C[17] * Ux[1] + C[3] * (Ux[0] + 1.0);
      dsdUx[28] = C[15] * Ux[2] + C[17] * (Ux[0] + 1.0) + C[8] * Ux[1];
      dsdUx[29] = C[12] * Ux[2] + C[15] * Ux[1] + C[16] * (Ux[0] + 1.0);
      dsdUx[30] = C[16] * Ux[5] + C[17] * (Ux[4] + 1.0) + C[3] * Ux[3];
      dsdUx[31] = C[15] * Ux[5] + C[17] * Ux[3] + C[8] * (Ux[4] + 1.0);
      dsdUx[32] = C[12] * Ux[5] + C[15] * (Ux[4] + 1.0) + C[16] * Ux[3];
      dsdUx[33] = C[16] * (Ux[8] + 1.0) + C[17] * Ux[7] + C[3] * Ux[6];
      dsdUx[34] = C[15] * (Ux[8] + 1.0) + C[17] * Ux[6] + C[8] * Ux[7];
      dsdUx[35] = C[12] * (Ux[8] + 1.0) + C[15] * Ux[7] + C[16] * Ux[6];

      dsdUx[36] = C[18] * Ux[2] + C[19] * Ux[1] + C[4] * (Ux[0] + 1.0);
      dsdUx[37] = C[16] * Ux[2] + C[19] * (Ux[0] + 1.0) + C[9] * Ux[1];
      dsdUx[38] = C[13] * Ux[2] + C[16] * Ux[1] + C[18] * (Ux[0] + 1.0);
      dsdUx[39] = C[18] * Ux[5] + C[19] * (Ux[4] + 1.0) + C[4] * Ux[3];
      dsdUx[40] = C[16] * Ux[5] + C[19] * Ux[3] + C[9] * (Ux[4] + 1.0);
      dsdUx[41] = C[13] * Ux[5] + C[16] * (Ux[4] + 1.0) + C[18] * Ux[3];
      dsdUx[42] = C[18] * (Ux[8] + 1.0) + C[19] * Ux[7] + C[4] * Ux[6];
      dsdUx[43] = C[16] * (Ux[8] + 1.0) + C[19] * Ux[6] + C[9] * Ux[7];
      dsdUx[44] = C[13] * (Ux[8] + 1.0) + C[16] * Ux[7] + C[18] * Ux[6];

      dsdUx[45] = C[19] * Ux[2] + C[20] * Ux[1] + C[5] * (Ux[0] + 1.0);
      dsdUx[46] = C[10] * Ux[1] + C[17] * Ux[2] + C[20] * (Ux[0] + 1.0);
      dsdUx[47] = C[14] * Ux[2] + C[17] * Ux[1] + C[19] * (Ux[0] + 1.0);
      dsdUx[48] = C[19] * Ux[5] + C[20] * (Ux[4] + 1.0) + C[5] * Ux[3];
      dsdUx[49] = C[10] * (Ux[4] + 1.0) + C[17] * Ux[5] + C[20] * Ux[3];
      dsdUx[50] = C[14] * Ux[5] + C[17] * (Ux[4] + 1.0) + C[19] * Ux[3];
      dsdUx[51] = C[19] * (Ux[8] + 1.0) + C[20] * Ux[7] + C[5] * Ux[6];
      dsdUx[52] = C[10] * Ux[7] + C[17] * (Ux[8] + 1.0) + C[20] * Ux[6];
      dsdUx[53] = C[14] * (Ux[8] + 1.0) + C[17] * Ux[7] + C[19] * Ux[6];

      // Set the Jacobian coefficients
      Jac[3] = Ux[1] * dsdUx[45] + Ux[2] * dsdUx[36] +
               dsdUx[0] * (Ux[0] + 1.0) + s[0];
      Jac[4] = Ux[1] * dsdUx[46] + Ux[2] * dsdUx[37] +
               dsdUx[1] * (Ux[0] + 1.0) + s[5];
      Jac[5] = Ux[1] * dsdUx[47] + Ux[2] * dsdUx[38] +
               dsdUx[2] * (Ux[0] + 1.0) + s[4];
      Jac[6] = Ux[1] * dsdUx[48] + Ux[2] * dsdUx[39] + dsdUx[3] * (Ux[0] + 1.0);
      Jac[7] = Ux[1] * dsdUx[49] + Ux[2] * dsdUx[40] + dsdUx[4] * (Ux[0] + 1.0);
      Jac[8] = Ux[1] * dsdUx[50] + Ux[2] * dsdUx[41] + dsdUx[5] * (Ux[0] + 1.0);
      Jac[9] = Ux[1] * dsdUx[51] + Ux[2] * dsdUx[42] + dsdUx[6] * (Ux[0] + 1.0);
      Jac[10] =
          Ux[1] * dsdUx[52] + Ux[2] * dsdUx[43] + dsdUx[7] * (Ux[0] + 1.0);
      Jac[11] =
          Ux[1] * dsdUx[53] + Ux[2] * dsdUx[44] + dsdUx[8] * (Ux[0] + 1.0);

      Jac[12] = Ux[1] * dsdUx[9] + Ux[2] * dsdUx[27] +
                dsdUx[45] * (Ux[0] + 1.0) + s[5];
      Jac[13] = Ux[1] * dsdUx[10] + Ux[2] * dsdUx[28] +
                dsdUx[46] * (Ux[0] + 1.0) + s[1];
      Jac[14] = Ux[1] * dsdUx[11] + Ux[2] * dsdUx[29] +
                dsdUx[47] * (Ux[0] + 1.0) + s[3];
      Jac[15] =
          Ux[1] * dsdUx[12] + Ux[2] * dsdUx[30] + dsdUx[48] * (Ux[0] + 1.0);
      Jac[16] =
          Ux[1] * dsdUx[13] + Ux[2] * dsdUx[31] + dsdUx[49] * (Ux[0] + 1.0);
      Jac[17] =
          Ux[1] * dsdUx[14] + Ux[2] * dsdUx[32] + dsdUx[50] * (Ux[0] + 1.0);
      Jac[18] =
          Ux[1] * dsdUx[15] + Ux[2] * dsdUx[33] + dsdUx[51] * (Ux[0] + 1.0);
      Jac[19] =
          Ux[1] * dsdUx[16] + Ux[2] * dsdUx[34] + dsdUx[52] * (Ux[0] + 1.0);
      Jac[20] =
          Ux[1] * dsdUx[17] + Ux[2] * dsdUx[35] + dsdUx[53] * (Ux[0] + 1.0);

      Jac[21] = Ux[1] * dsdUx[27] + Ux[2] * dsdUx[18] +
                dsdUx[36] * (Ux[0] + 1.0) + s[4];
      Jac[22] = Ux[1] * dsdUx[28] + Ux[2] * dsdUx[19] +
                dsdUx[37] * (Ux[0] + 1.0) + s[3];
      Jac[23] = Ux[1] * dsdUx[29] + Ux[2] * dsdUx[20] +
                dsdUx[38] * (Ux[0] + 1.0) + s[2];
      Jac[24] =
          Ux[1] * dsdUx[30] + Ux[2] * dsdUx[21] + dsdUx[39] * (Ux[0] + 1.0);
      Jac[25] =
          Ux[1] * dsdUx[31] + Ux[2] * dsdUx[22] + dsdUx[40] * (Ux[0] + 1.0);
      Jac[26] =
          Ux[1] * dsdUx[32] + Ux[2] * dsdUx[23] + dsdUx[41] * (Ux[0] + 1.0);
      Jac[27] =
          Ux[1] * dsdUx[33] + Ux[2] * dsdUx[24] + dsdUx[42] * (Ux[0] + 1.0);
      Jac[28] =
          Ux[1] * dsdUx[34] + Ux[2] * dsdUx[25] + dsdUx[43] * (Ux[0] + 1.0);
      Jac[29] =
          Ux[1] * dsdUx[35] + Ux[2] * dsdUx[26] + dsdUx[44] * (Ux[0] + 1.0);

      Jac[30] =
          Ux[3] * dsdUx[0] + Ux[5] * dsdUx[36] + dsdUx[45] * (Ux[4] + 1.0);
      Jac[31] =
          Ux[3] * dsdUx[1] + Ux[5] * dsdUx[37] + dsdUx[46] * (Ux[4] + 1.0);
      Jac[32] =
          Ux[3] * dsdUx[2] + Ux[5] * dsdUx[38] + dsdUx[47] * (Ux[4] + 1.0);
      Jac[33] = Ux[3] * dsdUx[3] + Ux[5] * dsdUx[39] +
                dsdUx[48] * (Ux[4] + 1.0) + s[0];
      Jac[34] = Ux[3] * dsdUx[4] + Ux[5] * dsdUx[40] +
                dsdUx[49] * (Ux[4] + 1.0) + s[5];
      Jac[35] = Ux[3] * dsdUx[5] + Ux[5] * dsdUx[41] +
                dsdUx[50] * (Ux[4] + 1.0) + s[4];
      Jac[36] =
          Ux[3] * dsdUx[6] + Ux[5] * dsdUx[42] + dsdUx[51] * (Ux[4] + 1.0);
      Jac[37] =
          Ux[3] * dsdUx[7] + Ux[5] * dsdUx[43] + dsdUx[52] * (Ux[4] + 1.0);
      Jac[38] =
          Ux[3] * dsdUx[8] + Ux[5] * dsdUx[44] + dsdUx[53] * (Ux[4] + 1.0);

      Jac[39] =
          Ux[3] * dsdUx[45] + Ux[5] * dsdUx[27] + dsdUx[9] * (Ux[4] + 1.0);
      Jac[40] =
          Ux[3] * dsdUx[46] + Ux[5] * dsdUx[28] + dsdUx[10] * (Ux[4] + 1.0);
      Jac[41] =
          Ux[3] * dsdUx[47] + Ux[5] * dsdUx[29] + dsdUx[11] * (Ux[4] + 1.0);
      Jac[42] = Ux[3] * dsdUx[48] + Ux[5] * dsdUx[30] +
                dsdUx[12] * (Ux[4] + 1.0) + s[5];
      Jac[43] = Ux[3] * dsdUx[49] + Ux[5] * dsdUx[31] +
                dsdUx[13] * (Ux[4] + 1.0) + s[1];
      Jac[44] = Ux[3] * dsdUx[50] + Ux[5] * dsdUx[32] +
                dsdUx[14] * (Ux[4] + 1.0) + s[3];
      Jac[45] =
          Ux[3] * dsdUx[51] + Ux[5] * dsdUx[33] + dsdUx[15] * (Ux[4] + 1.0);
      Jac[46] =
          Ux[3] * dsdUx[52] + Ux[5] * dsdUx[34] + dsdUx[16] * (Ux[4] + 1.0);
      Jac[47] =
          Ux[3] * dsdUx[53] + Ux[5] * dsdUx[35] + dsdUx[17] * (Ux[4] + 1.0);

      Jac[48] =
          Ux[3] * dsdUx[36] + Ux[5] * dsdUx[18] + dsdUx[27] * (Ux[4] + 1.0);
      Jac[49] =
          Ux[3] * dsdUx[37] + Ux[5] * dsdUx[19] + dsdUx[28] * (Ux[4] + 1.0);
      Jac[50] =
          Ux[3] * dsdUx[38] + Ux[5] * dsdUx[20] + dsdUx[29] * (Ux[4] + 1.0);
      Jac[51] = Ux[3] * dsdUx[39] + Ux[5] * dsdUx[21] +
                dsdUx[30] * (Ux[4] + 1.0) + s[4];
      Jac[52] = Ux[3] * dsdUx[40] + Ux[5] * dsdUx[22] +
                dsdUx[31] * (Ux[4] + 1.0) + s[3];
      Jac[53] = Ux[3] * dsdUx[41] + Ux[5] * dsdUx[23] +
                dsdUx[32] * (Ux[4] + 1.0) + s[2];
      Jac[54] =
          Ux[3] * dsdUx[42] + Ux[5] * dsdUx[24] + dsdUx[33] * (Ux[4] + 1.0);
      Jac[55] =
          Ux[3] * dsdUx[43] + Ux[5] * dsdUx[25] + dsdUx[34] * (Ux[4] + 1.0);
      Jac[56] =
          Ux[3] * dsdUx[44] + Ux[5] * dsdUx[26] + dsdUx[35] * (Ux[4] + 1.0);

      Jac[57] =
          Ux[6] * dsdUx[0] + Ux[7] * dsdUx[45] + dsdUx[36] * (Ux[8] + 1.0);
      Jac[58] =
          Ux[6] * dsdUx[1] + Ux[7] * dsdUx[46] + dsdUx[37] * (Ux[8] + 1.0);
      Jac[59] =
          Ux[6] * dsdUx[2] + Ux[7] * dsdUx[47] + dsdUx[38] * (Ux[8] + 1.0);
      Jac[60] =
          Ux[6] * dsdUx[3] + Ux[7] * dsdUx[48] + dsdUx[39] * (Ux[8] + 1.0);
      Jac[61] =
          Ux[6] * dsdUx[4] + Ux[7] * dsdUx[49] + dsdUx[40] * (Ux[8] + 1.0);
      Jac[62] =
          Ux[6] * dsdUx[5] + Ux[7] * dsdUx[50] + dsdUx[41] * (Ux[8] + 1.0);
      Jac[63] = Ux[6] * dsdUx[6] + Ux[7] * dsdUx[51] +
                dsdUx[42] * (Ux[8] + 1.0) + s[0];
      Jac[64] = Ux[6] * dsdUx[7] + Ux[7] * dsdUx[52] +
                dsdUx[43] * (Ux[8] + 1.0) + s[5];
      Jac[65] = Ux[6] * dsdUx[8] + Ux[7] * dsdUx[53] +
                dsdUx[44] * (Ux[8] + 1.0) + s[4];

      Jac[66] =
          Ux[6] * dsdUx[45] + Ux[7] * dsdUx[9] + dsdUx[27] * (Ux[8] + 1.0);
      Jac[67] =
          Ux[6] * dsdUx[46] + Ux[7] * dsdUx[10] + dsdUx[28] * (Ux[8] + 1.0);
      Jac[68] =
          Ux[6] * dsdUx[47] + Ux[7] * dsdUx[11] + dsdUx[29] * (Ux[8] + 1.0);
      Jac[69] =
          Ux[6] * dsdUx[48] + Ux[7] * dsdUx[12] + dsdUx[30] * (Ux[8] + 1.0);
      Jac[70] =
          Ux[6] * dsdUx[49] + Ux[7] * dsdUx[13] + dsdUx[31] * (Ux[8] + 1.0);
      Jac[71] =
          Ux[6] * dsdUx[50] + Ux[7] * dsdUx[14] + dsdUx[32] * (Ux[8] + 1.0);
      Jac[72] = Ux[6] * dsdUx[51] + Ux[7] * dsdUx[15] +
                dsdUx[33] * (Ux[8] + 1.0) + s[5];
      Jac[73] = Ux[6] * dsdUx[52] + Ux[7] * dsdUx[16] +
                dsdUx[34] * (Ux[8] + 1.0) + s[1];
      Jac[74] = Ux[6] * dsdUx[53] + Ux[7] * dsdUx[17] +
                dsdUx[35] * (Ux[8] + 1.0) + s[3];

      Jac[75] =
          Ux[6] * dsdUx[36] + Ux[7] * dsdUx[27] + dsdUx[18] * (Ux[8] + 1.0);
      Jac[76] =
          Ux[6] * dsdUx[37] + Ux[7] * dsdUx[28] + dsdUx[19] * (Ux[8] + 1.0);
      Jac[77] =
          Ux[6] * dsdUx[38] + Ux[7] * dsdUx[29] + dsdUx[20] * (Ux[8] + 1.0);
      Jac[78] =
          Ux[6] * dsdUx[39] + Ux[7] * dsdUx[30] + dsdUx[21] * (Ux[8] + 1.0);
      Jac[79] =
          Ux[6] * dsdUx[40] + Ux[7] * dsdUx[31] + dsdUx[22] * (Ux[8] + 1.0);
      Jac[80] =
          Ux[6] * dsdUx[41] + Ux[7] * dsdUx[32] + dsdUx[23] * (Ux[8] + 1.0);
      Jac[81] = Ux[6] * dsdUx[42] + Ux[7] * dsdUx[33] +
                dsdUx[24] * (Ux[8] + 1.0) + s[4];
      Jac[82] = Ux[6] * dsdUx[43] + Ux[7] * dsdUx[34] +
                dsdUx[25] * (Ux[8] + 1.0) + s[3];
      Jac[83] = Ux[6] * dsdUx[44] + Ux[7] * dsdUx[35] +
                dsdUx[26] * (Ux[8] + 1.0) + s[2];
    }
  } else if (matType == TACS_MASS_MATRIX) {
    // Evaluate the density
    TacsScalar rho = stiff->evalMassMatrixDensity(elemIndex, pt, X);

    // Set the acceleration terms
    Jac[0] = rho;
    Jac[1] = rho;
    Jac[2] = rho;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    TacsScalar C[21];
    stiff->evalTangentStiffness(elemIndex, pt, X, C);

    Jac[0] = C[0];
    Jac[1] = C[5];
    Jac[2] = C[4];
    Jac[3] = C[5];
    Jac[4] = C[1];
    Jac[5] = C[3];
    Jac[6] = C[4];
    Jac[7] = C[3];
    Jac[8] = C[2];

    Jac[9] = C[5];
    Jac[10] = C[20];
    Jac[11] = C[19];
    Jac[12] = C[20];
    Jac[13] = C[10];
    Jac[14] = C[17];
    Jac[15] = C[19];
    Jac[16] = C[17];
    Jac[17] = C[14];

    Jac[18] = C[4];
    Jac[19] = C[19];
    Jac[20] = C[18];
    Jac[21] = C[19];
    Jac[22] = C[9];
    Jac[23] = C[16];
    Jac[24] = C[18];
    Jac[25] = C[16];
    Jac[26] = C[13];

    Jac[27] = C[5];
    Jac[28] = C[20];
    Jac[29] = C[19];
    Jac[30] = C[20];
    Jac[31] = C[10];
    Jac[32] = C[17];
    Jac[33] = C[19];
    Jac[34] = C[17];
    Jac[35] = C[14];

    Jac[36] = C[1];
    Jac[37] = C[10];
    Jac[38] = C[9];
    Jac[39] = C[10];
    Jac[40] = C[6];
    Jac[41] = C[8];
    Jac[42] = C[9];
    Jac[43] = C[8];
    Jac[44] = C[7];

    Jac[45] = C[3];
    Jac[46] = C[17];
    Jac[47] = C[16];
    Jac[48] = C[17];
    Jac[49] = C[8];
    Jac[50] = C[15];
    Jac[51] = C[16];
    Jac[52] = C[15];
    Jac[53] = C[12];

    Jac[54] = C[4];
    Jac[55] = C[19];
    Jac[56] = C[18];
    Jac[57] = C[19];
    Jac[58] = C[9];
    Jac[59] = C[16];
    Jac[60] = C[18];
    Jac[61] = C[16];
    Jac[62] = C[13];

    Jac[63] = C[3];
    Jac[64] = C[17];
    Jac[65] = C[16];
    Jac[66] = C[17];
    Jac[67] = C[8];
    Jac[68] = C[15];
    Jac[69] = C[16];
    Jac[70] = C[15];
    Jac[71] = C[12];

    Jac[72] = C[2];
    Jac[73] = C[14];
    Jac[74] = C[13];
    Jac[75] = C[14];
    Jac[76] = C[7];
    Jac[77] = C[12];
    Jac[78] = C[13];
    Jac[79] = C[12];
    Jac[80] = C[11];
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    // Compute the tangent stiffness
    TacsScalar C[21];
    stiff->evalGeometricTangentStiffness(elemIndex, pt, X, C);

    // Compute the tangent strain and stress
    TacsScalar e[6], s[6];
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];

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

    Jac[0] = s[0];
    Jac[1] = s[5];
    Jac[2] = s[4];
    Jac[3] = 0.0;
    Jac[4] = 0.0;
    Jac[5] = 0.0;
    Jac[6] = 0.0;
    Jac[7] = 0.0;
    Jac[8] = 0.0;

    Jac[9] = s[5];
    Jac[10] = s[1];
    Jac[11] = s[3];
    Jac[12] = 0.0;
    Jac[13] = 0.0;
    Jac[14] = 0.0;
    Jac[15] = 0.0;
    Jac[16] = 0.0;
    Jac[17] = 0.0;

    Jac[18] = s[4];
    Jac[19] = s[3];
    Jac[20] = s[2];
    Jac[21] = 0.0;
    Jac[22] = 0.0;
    Jac[23] = 0.0;
    Jac[24] = 0.0;
    Jac[25] = 0.0;
    Jac[26] = 0.0;

    Jac[27] = 0.0;
    Jac[28] = 0.0;
    Jac[29] = 0.0;
    Jac[30] = s[0];
    Jac[31] = s[5];
    Jac[32] = s[4];
    Jac[33] = 0.0;
    Jac[34] = 0.0;
    Jac[35] = 0.0;

    Jac[36] = 0.0;
    Jac[37] = 0.0;
    Jac[38] = 0.0;
    Jac[39] = s[5];
    Jac[40] = s[1];
    Jac[41] = s[3];
    Jac[42] = 0.0;
    Jac[43] = 0.0;
    Jac[44] = 0.0;

    Jac[45] = 0.0;
    Jac[46] = 0.0;
    Jac[47] = 0.0;
    Jac[48] = s[4];
    Jac[49] = s[3];
    Jac[50] = s[2];
    Jac[51] = 0.0;
    Jac[52] = 0.0;
    Jac[53] = 0.0;

    Jac[54] = 0.0;
    Jac[55] = 0.0;
    Jac[56] = 0.0;
    Jac[57] = 0.0;
    Jac[58] = 0.0;
    Jac[59] = 0.0;
    Jac[60] = s[0];
    Jac[61] = s[5];
    Jac[62] = s[4];

    Jac[63] = 0.0;
    Jac[64] = 0.0;
    Jac[65] = 0.0;
    Jac[66] = 0.0;
    Jac[67] = 0.0;
    Jac[68] = 0.0;
    Jac[69] = s[5];
    Jac[70] = s[1];
    Jac[71] = s[3];

    Jac[72] = 0.0;
    Jac[73] = 0.0;
    Jac[74] = 0.0;
    Jac[75] = 0.0;
    Jac[76] = 0.0;
    Jac[77] = 0.0;
    Jac[78] = s[4];
    Jac[79] = s[3];
    Jac[80] = s[2];
  }
}

void TACSLinearElasticity3D::addWeakMatDVSens(
    ElementMatrixType matType, int elemIndex, const double time,
    TacsScalar scale, int n, const double pt[], const TacsScalar X[],
    const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
    const TacsScalar Psi[], const TacsScalar Psix[], const TacsScalar Phi[],
    const TacsScalar Phix[], int dvLen, TacsScalar *dfdx) {
  if (matType == TACS_MASS_MATRIX) {
    TacsScalar rho_coef = Psi[0] * Phi[0] + Psi[3] * Phi[3] + Psi[6] * Phi[6];

    stiff->addMassMatrixDensityDVSens(elemIndex, scale * rho_coef, pt, X, dvLen,
                                      dfdx);
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    TacsScalar ePsi[6];
    ePsi[0] = Psix[0];
    ePsi[1] = Psix[4];
    ePsi[2] = Psix[8];

    ePsi[3] = Psix[5] + Psix[7];
    ePsi[4] = Psix[2] + Psix[6];
    ePsi[5] = Psix[1] + Psix[3];

    TacsScalar ePhi[6];
    ePhi[0] = Phix[0];
    ePhi[1] = Phix[4];
    ePhi[2] = Phix[8];

    ePhi[3] = Phix[5] + Phix[7];
    ePhi[4] = Phix[2] + Phix[6];
    ePhi[5] = Phix[1] + Phix[3];

    stiff->addStressDVSens(elemIndex, scale, pt, X, ePsi, ePhi, dvLen, dfdx);
  } else if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TacsScalar e[6];
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];

    TacsScalar psi[6];
    psi[0] = Phix[0] * Psix[0] + Phix[3] * Psix[3] + Phix[6] * Psix[6];
    psi[1] = Phix[1] * Psix[1] + Phix[4] * Psix[4] + Phix[7] * Psix[7];
    psi[2] = Phix[2] * Psix[2] + Phix[5] * Psix[5] + Phix[8] * Psix[8];

    psi[3] = Phix[1] * Psix[2] + Phix[2] * Psix[1] + Phix[4] * Psix[5] +
             Phix[5] * Psix[4] + Phix[7] * Psix[8] + Phix[8] * Psix[7];
    psi[4] = Phix[0] * Psix[2] + Phix[2] * Psix[0] + Phix[3] * Psix[5] +
             Phix[5] * Psix[3] + Phix[6] * Psix[8] + Phix[8] * Psix[6];
    psi[5] = Phix[0] * Psix[1] + Phix[1] * Psix[0] + Phix[3] * Psix[4] +
             Phix[4] * Psix[3] + Phix[6] * Psix[7] + Phix[7] * Psix[6];

    stiff->addGeometricTangentStressDVSens(elemIndex, scale, pt, X, e, psi,
                                           dvLen, dfdx);
  }
}

void TACSLinearElasticity3D::evalWeakMatSVSens(
    ElementMatrixType matType, int elemIndex, const double time,
    TacsScalar scale, int n, const double pt[], const TacsScalar X[],
    const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
    const TacsScalar Psi[], const TacsScalar Psix[], const TacsScalar Phi[],
    const TacsScalar Phix[], TacsScalar dfdU[], TacsScalar dfdUx[]) {
  dfdU[0] = dfdU[1] = dfdU[2] = 0.0;
  dfdU[3] = dfdU[4] = dfdU[5] = 0.0;
  dfdU[6] = dfdU[7] = dfdU[8] = 0.0;
  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;

  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TacsScalar e[6];
    e[0] = Phix[0] * Psix[0] + Phix[3] * Psix[3] + Phix[6] * Psix[6];
    e[1] = Phix[1] * Psix[1] + Phix[4] * Psix[4] + Phix[7] * Psix[7];
    e[2] = Phix[2] * Psix[2] + Phix[5] * Psix[5] + Phix[8] * Psix[8];

    e[3] = Phix[1] * Psix[2] + Phix[2] * Psix[1] + Phix[4] * Psix[5] +
           Phix[5] * Psix[4] + Phix[7] * Psix[8] + Phix[8] * Psix[7];
    e[4] = Phix[0] * Psix[2] + Phix[2] * Psix[0] + Phix[3] * Psix[5] +
           Phix[5] * Psix[3] + Phix[6] * Psix[8] + Phix[8] * Psix[6];
    e[5] = Phix[0] * Psix[1] + Phix[1] * Psix[0] + Phix[3] * Psix[4] +
           Phix[4] * Psix[3] + Phix[6] * Psix[7] + Phix[7] * Psix[6];

    // Compute the tangent stiffness
    TacsScalar C[21];
    stiff->evalGeometricTangentStiffness(elemIndex, pt, X, C);

    // Compute the tangent strain and stress
    TacsScalar s[6];
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

    dfdUx[0] = scale * s[0];
    dfdUx[4] = scale * s[1];
    dfdUx[8] = scale * s[2];

    dfdUx[5] = scale * s[3];
    dfdUx[7] = scale * s[3];

    dfdUx[2] = scale * s[4];
    dfdUx[6] = scale * s[4];

    dfdUx[1] = scale * s[5];
    dfdUx[3] = scale * s[5];
  }
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearElasticity3D::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
    if (quantity) {
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

      *quantity = stiff->evalFailure(elemIndex, pt, X, e);
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = stiff->evalDensity(elemIndex, pt, X);
    }

    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
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
void TACSLinearElasticity3D::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX) {
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

    stiff->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, e, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    stiff->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
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
void TACSLinearElasticity3D::evalPointQuantitySens(
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

  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX) {
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
    } else {
      dfdUx[0] = dfdq[0] * (Ux[0] * sens[0] + Ux[1] * sens[5] +
                            Ux[2] * sens[4] + sens[0]);
      dfdUx[1] = dfdq[0] *
                 (Ux[1] * sens[1] + Ux[2] * sens[3] + sens[5] * (Ux[0] + 1.0));
      dfdUx[2] = dfdq[0] *
                 (Ux[1] * sens[3] + Ux[2] * sens[2] + sens[4] * (Ux[0] + 1.0));

      dfdUx[3] = dfdq[0] *
                 (Ux[3] * sens[0] + Ux[5] * sens[4] + sens[5] * (Ux[4] + 1.0));
      dfdUx[4] = dfdq[0] * (Ux[3] * sens[5] + Ux[4] * sens[1] +
                            Ux[5] * sens[3] + sens[1]);
      dfdUx[5] = dfdq[0] *
                 (Ux[3] * sens[4] + Ux[5] * sens[2] + sens[3] * (Ux[4] + 1.0));

      dfdUx[6] = dfdq[0] *
                 (Ux[6] * sens[0] + Ux[7] * sens[5] + sens[4] * (Ux[8] + 1.0));
      dfdUx[7] = dfdq[0] *
                 (Ux[6] * sens[5] + Ux[7] * sens[1] + sens[3] * (Ux[8] + 1.0));
      dfdUx[8] = dfdq[0] * (Ux[6] * sens[4] + Ux[7] * sens[3] +
                            Ux[8] * sens[2] + sens[2]);
    }
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
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
    } else {
      dfdUx[0] =
          2.0 * dfdq[0] * (Ux[0] * s[0] + Ux[1] * s[5] + Ux[2] * s[4] + s[0]);
      dfdUx[1] =
          2.0 * dfdq[0] * (Ux[1] * s[1] + Ux[2] * s[3] + s[5] * (Ux[0] + 1.0));
      dfdUx[2] =
          2.0 * dfdq[0] * (Ux[1] * s[3] + Ux[2] * s[2] + s[4] * (Ux[0] + 1.0));

      dfdUx[3] =
          2.0 * dfdq[0] * (Ux[3] * s[0] + Ux[5] * s[4] + s[5] * (Ux[4] + 1.0));
      dfdUx[4] =
          2.0 * dfdq[0] * (Ux[3] * s[5] + Ux[4] * s[1] + Ux[5] * s[3] + s[1]);
      dfdUx[5] =
          2.0 * dfdq[0] * (Ux[3] * s[4] + Ux[5] * s[2] + s[3] * (Ux[4] + 1.0));

      dfdUx[6] =
          2.0 * dfdq[0] * (Ux[6] * s[0] + Ux[7] * s[5] + s[4] * (Ux[8] + 1.0));
      dfdUx[7] =
          2.0 * dfdq[0] * (Ux[6] * s[5] + Ux[7] * s[1] + s[3] * (Ux[8] + 1.0));
      dfdUx[8] =
          2.0 * dfdq[0] * (Ux[6] * s[4] + Ux[7] * s[3] + Ux[8] * s[2] + s[2]);
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
void TACSLinearElasticity3D::getOutputData(
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
