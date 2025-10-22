/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2020 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSThermoelasticPlateModel.h"

/*
  Evaluate the strain components for the plate.

  These strain components consist of the in-plane strains,
  curvatures and the out-of-plane shear strains.

  These strains are organized in the following manner:

  (exx, eyy, gxy, kxx, kyy, kxy, gyz, gxz)

  The strains are computed based on the assumed displacement expression:

  u(x,y,z) = u + z*roty
  v(x,y,z) = v - z*rotx
  w(x,y,z) = w

  Note the rotations (rotx, roty) are oriented about the x and y axes,
  respectively. Based on this definition, the strain components are given by

  (u,x; v,y; u,y + v,x; roty,x; -rotx,y; roty,y - rotx,x; w,y - rotx; w,x +
  roty)
*/
static inline void evalPlateStrain(const TacsScalar Ut[], const TacsScalar Ux[],
                                   const TacsScalar et[], TacsScalar e[]) {
  //   0,  1,  2   3   4   5      6       7      8      9
  // u,x u,y v,x v,y w,x w,y rotx,x  rotx,y roty,x roty,y
  e[0] = Ux[0] - et[0];          // exx = u,x
  e[1] = Ux[3] - et[1];          // eyy = v,y
  e[2] = Ux[1] + Ux[2] - et[2];  // exy = u,y + v,x

  e[3] = Ux[8] - et[3];          // kxx = roty,x
  e[4] = -Ux[7] - et[4];         // kyy = -rotx,y
  e[5] = Ux[9] - Ux[6] - et[5];  // kxy = roty,y - rotx,x

  e[6] = Ux[5] - Ut[9] - et[6];   // eyz = w,y - rotx
  e[7] = Ux[4] + Ut[12] - et[7];  // exz = w,x + roty

  e[8] = 0.0;
}

static inline void evalPlateStrain(const TacsScalar Ut[], const TacsScalar Ux[],
                                   const TacsScalar theta,
                                   const TacsScalar et[], TacsScalar e[]) {
  //   0,  1,  2   3   4   5      6       7      8      9
  // u,x u,y v,x v,y w,x w,y rotx,x  rotx,y roty,x roty,y
  e[0] = Ux[0] - theta * et[0];          // exx = u,x
  e[1] = Ux[3] - theta * et[1];          // eyy = v,y
  e[2] = Ux[1] + Ux[2] - theta * et[2];  // exy = u,y + v,x

  e[3] = Ux[8] - theta * et[3];          // kxx = roty,x
  e[4] = -Ux[7] - theta * et[4];         // kyy = -rotx,y
  e[5] = Ux[9] - Ux[6] - theta * et[5];  // kxy = roty,y - rotx,x

  e[6] = Ux[5] - Ut[9] - theta * et[6];   // eyz = w,y - rotx
  e[7] = Ux[4] + Ut[12] - theta * et[7];  // exz = w,x + roty

  e[8] = 0.0;
}

TACSThermoelasticPlateModel::TACSThermoelasticPlateModel(
    TACSShellConstitutive *_con) {
  con = _con;
  con->incref();
}

TACSThermoelasticPlateModel::~TACSThermoelasticPlateModel() { con->decref(); }

const int TACSThermoelasticPlateModel::linear_Jac_pairs[] = {
    2,  2,  3,  3,  3,  4,  3,  8,  3,  9,  3,  18, 3,  19, 3,  23, 3,  24, 3,
    25, 4,  3,  4,  4,  4,  8,  4,  9,  4,  18, 4,  19, 4,  23, 4,  24, 4,  25,
    7,  7,  8,  3,  8,  4,  8,  8,  8,  9,  8,  18, 8,  19, 8,  23, 8,  24, 8,
    25, 9,  3,  9,  4,  9,  8,  9,  9,  9,  18, 9,  19, 9,  23, 9,  24, 9,  25,
    12, 12, 13, 13, 13, 14, 13, 15, 13, 20, 13, 25, 14, 13, 14, 14, 14, 15, 14,
    20, 14, 25, 15, 13, 15, 14, 15, 15, 15, 20, 15, 25, 18, 3,  18, 4,  18, 8,
    18, 9,  18, 18, 18, 19, 18, 23, 18, 24, 18, 25, 19, 3,  19, 4,  19, 8,  19,
    9,  19, 18, 19, 19, 19, 23, 19, 24, 19, 25, 20, 13, 20, 14, 20, 15, 20, 20,
    20, 25, 23, 3,  23, 4,  23, 8,  23, 9,  23, 18, 23, 19, 23, 23, 23, 24, 23,
    25, 24, 3,  24, 4,  24, 8,  24, 9,  24, 18, 24, 19, 24, 23, 24, 24, 24, 25,
    26, 26, 28, 28, 28, 29, 29, 28, 29, 29};

int TACSThermoelasticPlateModel::getNumParameters() { return 2; }

int TACSThermoelasticPlateModel::getVarsPerNode() { return 6; }

int TACSThermoelasticPlateModel::getDesignVarsPerNode() {
  return con->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSThermoelasticPlateModel::getDesignVarNums(int elemIndex, int dvLen,
                                                  int dvNums[]) {
  return con->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSThermoelasticPlateModel::setDesignVars(int elemIndex, int dvLen,
                                               const TacsScalar dvs[]) {
  return con->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSThermoelasticPlateModel::getDesignVars(int elemIndex, int dvLen,
                                               TacsScalar dvs[]) {
  return con->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSThermoelasticPlateModel::getDesignVarRange(int elemIndex, int dvLen,
                                                   TacsScalar lb[],
                                                   TacsScalar ub[]) {
  return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSThermoelasticPlateModel::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  // Evaluate the density
  TacsScalar rho = con->evalDensity(elemIndex, pt, X);
  TacsScalar c = con->evalSpecificHeat(elemIndex, pt, X);

  // Compute the thermal strain
  TacsScalar et[9];
  TacsScalar theta = Ut[3 * 5];  // The temperature
  con->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the 8 components of the strain, plus a place-holder for the
  // in-plane rotation (zero always in this case)
  TacsScalar e[9];
  evalPlateStrain(Ut, Ux, et, e);

  // Evaluate the stress
  TacsScalar s[9];
  con->evalStress(elemIndex, pt, X, e, s);

  // Evaluate the heat flux from the thermal gradient
  TacsScalar grad[2], flux[2];
  grad[0] = Ux[10];
  grad[1] = Ux[11];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  con->evalHeatFlux(elemIndex, pt, X, grad, flux);

  // Set pointers to the in-plane, bending and shear resultants
  const TacsScalar *N = &s[0];
  const TacsScalar *M = &s[3];
  const TacsScalar *Q = &s[6];

  // Set the coefficients
  memset(DUt, 0, 5 * 3 * sizeof(TacsScalar));
  DUt[2] = rho * Ut[2];
  DUt[5] = rho * Ut[5];
  DUt[8] = rho * Ut[8];

  DUt[9] = -Q[0];  // rotx
  DUt[12] = Q[1];  // roty
  DUt[16] = c * rho * Ut[16];

  DUx[0] = N[0];  // u,x
  DUx[1] = N[2];  // u,y

  DUx[2] = N[2];  // v,x
  DUx[3] = N[1];  // v,y

  DUx[4] = Q[1];  // w,x
  DUx[5] = Q[0];  // w,y

  DUx[6] = -M[2];  // rotx,x
  DUx[7] = -M[1];  // rotx,y

  DUx[8] = M[0];  // roty,x
  DUx[9] = M[2];  // roty,y

  // Add the components from the heat flux
  DUx[10] = flux[0];
  DUx[11] = flux[1];
}

void TACSThermoelasticPlateModel::getWeakMatrixNonzeros(
    ElementMatrixType matType, int elemIndex, int *Jac_nnz,
    const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 100;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSThermoelasticPlateModel::evalWeakMatrix(
    ElementMatrixType matType, int elemIndex, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar DUt[],
    TacsScalar DUx[], TacsScalar Jac[]) {
  // Evaluate the density
  TacsScalar rho = con->evalDensity(elemIndex, pt, X);
  TacsScalar c = con->evalSpecificHeat(elemIndex, pt, X);

  // Compute the thermal strain
  TacsScalar et[9];
  TacsScalar theta = Ut[3 * 5];  // The temperature
  con->evalThermalStrain(elemIndex, pt, X, 1.0, et);

  // Compute the 8 components of the strain, plus a place-holder for the
  // in-plane rotation (zero always in this case)
  TacsScalar e[9];
  evalPlateStrain(Ut, Ux, theta, et, e);

  // Compute the stiffness matrix
  TacsScalar C[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  con->evalTangentStiffness(elemIndex, pt, X, C);

  // Extract the stiffnesses
  TacsScalar drill;
  const TacsScalar *A, *B, *D, *As;
  con->extractTangentStiffness(C, &A, &B, &D, &As, &drill);

  TacsScalar Kc[3];
  con->evalTangentHeatFlux(elemIndex, pt, X, Kc);

  // Get the stress components
  TacsScalar s[9];
  con->computeStress(A, B, D, As, drill, e, s);

  // Set pointers to the in-plane, bending and shear resultants
  const TacsScalar *N = &s[0];
  const TacsScalar *M = &s[3];
  const TacsScalar *Q = &s[6];

  // Set the coefficients
  memset(DUt, 0, 5 * 3 * sizeof(TacsScalar));
  DUt[2] = rho * Ut[2];
  DUt[5] = rho * Ut[5];
  DUt[8] = rho * Ut[8];

  DUt[9] = -Q[0];  // rotx
  DUt[12] = Q[1];  // roty
  DUt[16] = c * rho * Ut[16];

  DUx[0] = N[0];  // u,x
  DUx[1] = N[2];  // u,y

  DUx[2] = N[2];  // v,x
  DUx[3] = N[1];  // v,y

  DUx[4] = Q[1];  // w,x
  DUx[5] = Q[0];  // w,y

  DUx[6] = -M[2];  // rotx,x
  DUx[7] = -M[1];  // rotx,y

  DUx[8] = M[0];  // roty,x
  DUx[9] = M[2];  // roty,y

  // Add the components from the heat flux
  DUx[10] = Kc[0] * Ux[10] + Kc[1] * Ux[11];
  DUx[11] = Kc[1] * Ux[10] + Kc[2] * Ux[11];

  Jac[0] = rho;
  Jac[1] = A[0];
  Jac[2] = A[2];
  Jac[3] = A[2];
  Jac[4] = A[1];
  Jac[5] = -B[2];
  Jac[6] = -B[1];
  Jac[7] = B[0];
  Jac[8] = B[2];
  Jac[9] = -A[0] * et[0] - A[1] * et[1] - A[2] * et[2] - B[0] * et[3] -
           B[1] * et[4] - B[2] * et[5];
  Jac[10] = A[2];
  Jac[11] = A[5];
  Jac[12] = A[5];
  Jac[13] = A[4];
  Jac[14] = -B[5];
  Jac[15] = -B[4];
  Jac[16] = B[2];
  Jac[17] = B[5];
  Jac[18] = -A[2] * et[0] - A[4] * et[1] - A[5] * et[2] - B[2] * et[3] -
            B[4] * et[4] - B[5] * et[5];
  Jac[19] = rho;
  Jac[20] = A[2];
  Jac[21] = A[5];
  Jac[22] = A[5];
  Jac[23] = A[4];
  Jac[24] = -B[5];
  Jac[25] = -B[4];
  Jac[26] = B[2];
  Jac[27] = B[5];
  Jac[28] = -A[2] * et[0] - A[4] * et[1] - A[5] * et[2] - B[2] * et[3] -
            B[4] * et[4] - B[5] * et[5];
  Jac[29] = A[1];
  Jac[30] = A[4];
  Jac[31] = A[4];
  Jac[32] = A[3];
  Jac[33] = -B[4];
  Jac[34] = -B[3];
  Jac[35] = B[1];
  Jac[36] = B[4];
  Jac[37] = -A[1] * et[0] - A[3] * et[1] - A[4] * et[2] - B[1] * et[3] -
            B[3] * et[4] - B[4] * et[5];
  Jac[38] = rho;
  Jac[39] = As[2];
  Jac[40] = As[1];
  Jac[41] = -As[1];
  Jac[42] = As[2];
  Jac[43] = -As[1] * et[6] - As[2] * et[7];
  Jac[44] = As[1];
  Jac[45] = As[0];
  Jac[46] = -As[0];
  Jac[47] = As[1];
  Jac[48] = -As[0] * et[6] - As[1] * et[7];
  Jac[49] = -As[1];
  Jac[50] = -As[0];
  Jac[51] = As[0];
  Jac[52] = -As[1];
  Jac[53] = As[0] * et[6] + As[1] * et[7];
  Jac[54] = -B[2];
  Jac[55] = -B[5];
  Jac[56] = -B[5];
  Jac[57] = -B[4];
  Jac[58] = D[5];
  Jac[59] = D[4];
  Jac[60] = -D[2];
  Jac[61] = -D[5];
  Jac[62] = B[2] * et[0] + B[4] * et[1] + B[5] * et[2] + D[2] * et[3] +
            D[4] * et[4] + D[5] * et[5];
  Jac[63] = -B[1];
  Jac[64] = -B[4];
  Jac[65] = -B[4];
  Jac[66] = -B[3];
  Jac[67] = D[4];
  Jac[68] = D[3];
  Jac[69] = -D[1];
  Jac[70] = -D[4];
  Jac[71] = B[1] * et[0] + B[3] * et[1] + B[4] * et[2] + D[1] * et[3] +
            D[3] * et[4] + D[4] * et[5];
  Jac[72] = As[2];
  Jac[73] = As[1];
  Jac[74] = -As[1];
  Jac[75] = As[2];
  Jac[76] = -As[1] * et[6] - As[2] * et[7];
  Jac[77] = B[0];
  Jac[78] = B[2];
  Jac[79] = B[2];
  Jac[80] = B[1];
  Jac[81] = -D[2];
  Jac[82] = -D[1];
  Jac[83] = D[0];
  Jac[84] = D[2];
  Jac[85] = -B[0] * et[0] - B[1] * et[1] - B[2] * et[2] - D[0] * et[3] -
            D[1] * et[4] - D[2] * et[5];
  Jac[86] = B[2];
  Jac[87] = B[5];
  Jac[88] = B[5];
  Jac[89] = B[4];
  Jac[90] = -D[5];
  Jac[91] = -D[4];
  Jac[92] = D[2];
  Jac[93] = D[5];
  Jac[94] = -B[2] * et[0] - B[4] * et[1] - B[5] * et[2] - D[2] * et[3] -
            D[4] * et[4] - D[5] * et[5];
  Jac[95] = c * rho;
  Jac[96] = Kc[0];
  Jac[97] = Kc[1];
  Jac[98] = Kc[1];
  Jac[99] = Kc[2];
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSThermoelasticPlateModel::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
    TacsScalar e[9], et[9];
    TacsScalar theta = Ut[3 * 5];  // The temperature
    con->evalThermalStrain(elemIndex, pt, X, theta, et);
    evalPlateStrain(Ut, Ux, et, e);
    *quantity = con->evalFailure(elemIndex, pt, X, e);

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    *quantity = con->evalDensity(elemIndex, pt, X);

    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    TacsScalar e[9], et[9];
    TacsScalar theta = Ut[3 * 5];  // The temperature
    con->evalThermalStrain(elemIndex, pt, X, theta, et);
    evalPlateStrain(Ut, Ux, et, e);

    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    // Evaluate the strain energy density
    *quantity =
        (e[0] * s[0] + e[1] * s[1] + e[2] * s[2] + e[3] * s[3] + e[4] * s[4] +
         e[5] * s[5] + e[6] * s[6] + e[7] * s[7] + e[8] * s[8]);

    return 1;
  }

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSThermoelasticPlateModel::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {}

/*
  Evaluate the derivatives of the point-wise quantity of interest
  with respect to X, Ut and Ux.
*/
void TACSThermoelasticPlateModel::evalPointQuantitySens(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
    TacsScalar dfdUx[]) {}

/*
  Get the data for visualization at a given point
*/
void TACSThermoelasticPlateModel::getOutputData(
    int elemIndex, const double time, ElementType etype, int write_flag,
    const double pt[], const TacsScalar X[], const TacsScalar Ut[],
    const TacsScalar Ux[], int ld_data, TacsScalar *data) {
  if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    if (write_flag & TACS_OUTPUT_NODES) {
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      data[0] = Ut[0];
      data[1] = Ut[3];
      data[2] = Ut[6];
      data[3] = Ut[9];
      data[4] = Ut[12];
      data[5] = 0.0;
      data += 6;
    }

    TacsScalar e[9], et[9];
    TacsScalar theta = Ut[3 * 5];  // The temperature
    con->evalThermalStrain(elemIndex, pt, X, theta, et);
    evalPlateStrain(Ut, Ux, et, e);

    if (write_flag & TACS_OUTPUT_STRAINS) {
      for (int i = 0; i < 9; i++) {
        data[i] = e[i];
      }
      data += 9;
    }
    if (write_flag & TACS_OUTPUT_STRESSES) {
      TacsScalar s[9];
      con->evalStress(elemIndex, pt, X, e, s);
      for (int i = 0; i < 9; i++) {
        data[i] = s[i];
      }
      data += 9;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = con->evalFailure(elemIndex, pt, X, e);
      data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}
