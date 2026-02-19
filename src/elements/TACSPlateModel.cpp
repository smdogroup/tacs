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

#include "TACSPlateModel.h"

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
                                   TacsScalar e[]) {
  //   0,  1,  2   3   4   5      6       7      8      9
  // u,x u,y v,x v,y w,x w,y rotx,x  rotx,y roty,x roty,y
  e[0] = Ux[0];          // exx = u,x
  e[1] = Ux[3];          // eyy = v,y
  e[2] = Ux[1] + Ux[2];  // exy = u,y + v,x

  e[3] = Ux[8];          // kxx = roty,x
  e[4] = -Ux[7];         // kyy = -rotx,y
  e[5] = Ux[9] - Ux[6];  // kxy = roty,y - rotx,x

  e[6] = Ux[5] - Ut[9];   // eyz = w,y - rotx
  e[7] = Ux[4] + Ut[12];  // exz = w,x + roty

  e[8] = 0.0;
}

TACSPlateModel::TACSPlateModel(TACSShellConstitutive *_con) {
  con = _con;
  con->incref();
}

TACSPlateModel::~TACSPlateModel() { con->decref(); }

const int TACSPlateModel::linear_Jac_pairs[] = {
    2,  2,  3,  3,  3,  4,  3,  8,  3,  9,  3,  18, 3,  19, 3,  23, 3,  24, 4,
    3,  4,  4,  4,  8,  4,  9,  4,  18, 4,  19, 4,  23, 4,  24, 7,  7,  8,  3,
    8,  4,  8,  8,  8,  9,  8,  18, 8,  19, 8,  23, 8,  24, 9,  3,  9,  4,  9,
    8,  9,  9,  9,  18, 9,  19, 9,  23, 9,  24, 12, 12, 13, 13, 13, 14, 13, 15,
    13, 20, 14, 13, 14, 14, 14, 15, 14, 20, 15, 13, 15, 14, 15, 15, 15, 20, 18,
    3,  18, 4,  18, 8,  18, 9,  18, 18, 18, 19, 18, 23, 18, 24, 19, 3,  19, 4,
    19, 8,  19, 9,  19, 18, 19, 19, 19, 23, 19, 24, 20, 13, 20, 14, 20, 15, 20,
    20, 23, 3,  23, 4,  23, 8,  23, 9,  23, 18, 23, 19, 23, 23, 23, 24, 24, 3,
    24, 4,  24, 8,  24, 9,  24, 18, 24, 19, 24, 23, 24, 24};

int TACSPlateModel::getNumParameters() { return 2; }

int TACSPlateModel::getVarsPerNode() { return 5; }

int TACSPlateModel::getDesignVarsPerNode() {
  return con->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSPlateModel::getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
  return con->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSPlateModel::setDesignVars(int elemIndex, int dvLen,
                                  const TacsScalar dvs[]) {
  return con->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSPlateModel::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  return con->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSPlateModel::getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                                      TacsScalar ub[]) {
  return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSPlateModel::evalWeakIntegrand(int elemIndex, const double time, int n,
                                       const double pt[], const TacsScalar X[],
                                       const TacsScalar Xd[],
                                       const TacsScalar Ut[],
                                       const TacsScalar Ux[], TacsScalar DUt[],
                                       TacsScalar DUx[]) {
  // Evaluate the density
  TacsScalar rho = con->evalDensity(elemIndex, pt, X);

  // Compute the 8 components of the strain, plus a place-holder for the
  // in-plane rotation (zero always in this case)
  TacsScalar e[9];
  evalPlateStrain(Ut, Ux, e);

  // Evaluate the stress
  TacsScalar s[9];
  con->evalStress(elemIndex, pt, X, e, s);

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
}

void TACSPlateModel::getWeakMatrixNonzeros(ElementMatrixType matType,
                                           int elemIndex, int *Jac_nnz,
                                           const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 83;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

void TACSPlateModel::evalWeakMatrix(ElementMatrixType matType, int elemIndex,
                                    const double time, int n, const double pt[],
                                    const TacsScalar X[], const TacsScalar Xd[],
                                    const TacsScalar Ut[],
                                    const TacsScalar Ux[], TacsScalar DUt[],
                                    TacsScalar DUx[], TacsScalar Jac[]) {
  // Evaluate the density
  TacsScalar rho = con->evalDensity(elemIndex, pt, X);

  // Compute the 8 components of the strain, plus a place-holder for the
  // in-plane rotation (zero always in this case)
  TacsScalar e[9];
  evalPlateStrain(Ut, Ux, e);

  // Compute the stiffness matrix
  TacsScalar C[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  con->evalTangentStiffness(elemIndex, pt, X, C);

  // Extract the stiffnesses
  TacsScalar drill;
  const TacsScalar *A, *B, *D, *As;
  con->extractTangentStiffness(C, &A, &B, &D, &As, &drill);

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

  Jac[0] = rho;
  Jac[1] = A[0];
  Jac[2] = A[2];
  Jac[3] = A[2];
  Jac[4] = A[1];
  Jac[5] = -B[2];
  Jac[6] = -B[1];
  Jac[7] = B[0];
  Jac[8] = B[2];
  Jac[9] = A[2];
  Jac[10] = A[5];
  Jac[11] = A[5];
  Jac[12] = A[4];
  Jac[13] = -B[5];
  Jac[14] = -B[4];
  Jac[15] = B[2];
  Jac[16] = B[5];

  Jac[17] = rho;
  Jac[18] = A[2];
  Jac[19] = A[5];
  Jac[20] = A[5];
  Jac[21] = A[4];
  Jac[22] = -B[5];
  Jac[23] = -B[4];
  Jac[24] = B[2];
  Jac[25] = B[5];
  Jac[26] = A[1];
  Jac[27] = A[4];
  Jac[28] = A[4];
  Jac[29] = A[3];
  Jac[30] = -B[4];
  Jac[31] = -B[3];
  Jac[32] = B[1];
  Jac[33] = B[4];

  Jac[34] = rho;
  Jac[35] = As[2];
  Jac[36] = As[1];
  Jac[37] = -As[1];
  Jac[38] = As[2];
  Jac[39] = As[1];
  Jac[40] = As[0];
  Jac[41] = -As[0];
  Jac[42] = As[1];
  Jac[43] = -As[1];
  Jac[44] = -As[0];
  Jac[45] = As[0];
  Jac[46] = -As[1];
  Jac[47] = -B[2];
  Jac[48] = -B[5];
  Jac[49] = -B[5];
  Jac[50] = -B[4];
  Jac[51] = D[5];
  Jac[52] = D[4];
  Jac[53] = -D[2];
  Jac[54] = -D[5];
  Jac[55] = -B[1];
  Jac[56] = -B[4];
  Jac[57] = -B[4];
  Jac[58] = -B[3];
  Jac[59] = D[4];
  Jac[60] = D[3];
  Jac[61] = -D[1];
  Jac[62] = -D[4];
  Jac[63] = As[2];
  Jac[64] = As[1];
  Jac[65] = -As[1];
  Jac[66] = As[2];
  Jac[67] = B[0];
  Jac[68] = B[2];
  Jac[69] = B[2];
  Jac[70] = B[1];
  Jac[71] = -D[2];
  Jac[72] = -D[1];
  Jac[73] = D[0];
  Jac[74] = D[2];
  Jac[75] = B[2];
  Jac[76] = B[5];
  Jac[77] = B[5];
  Jac[78] = B[4];
  Jac[79] = -D[5];
  Jac[80] = -D[4];
  Jac[81] = D[2];
  Jac[82] = D[5];
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSPlateModel::evalPointQuantity(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], TacsScalar *quantity) {
  if (quantityType == TACS_FAILURE_INDEX) {
    TacsScalar e[9];
    evalPlateStrain(Ut, Ux, e);
    *quantity = con->evalFailure(elemIndex, pt, X, e);

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    *quantity = con->evalDensity(elemIndex, pt, X);

    return 1;
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    TacsScalar e[9];
    evalPlateStrain(Ut, Ux, e);

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
void TACSPlateModel::addPointQuantityDVSens(
    int elemIndex, const int quantityType, const double time, TacsScalar scale,
    int n, const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    int dvLen, TacsScalar dfdx[]) {}

/*
  Evaluate the derivatives of the point-wise quantity of interest
  with respect to X, Ut and Ux.
*/
void TACSPlateModel::evalPointQuantitySens(
    int elemIndex, const int quantityType, const double time, int n,
    const double pt[], const TacsScalar X[], const TacsScalar Xd[],
    const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
    TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
    TacsScalar dfdUx[]) {}

/*
  Get the data for visualization at a given point
*/
void TACSPlateModel::getOutputData(int elemIndex, const double time,
                                   ElementType etype, int write_flag,
                                   const double pt[], const TacsScalar X[],
                                   const TacsScalar Ut[], const TacsScalar Ux[],
                                   int ld_data, TacsScalar *data) {
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

    TacsScalar e[9];
    evalPlateStrain(Ut, Ux, e);

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
