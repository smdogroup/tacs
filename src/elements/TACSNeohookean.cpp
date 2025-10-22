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

#include "TACSNeohookean.h"

#include "TACSElementAlgebra.h"

TACSNeohookean3D::TACSNeohookean3D(TacsScalar _C1, TacsScalar _D1) {
  C1 = _C1;
  D1 = _D1;
}

TACSNeohookean3D::~TACSNeohookean3D() {}

const int TACSNeohookean3D::linear_Jac_pairs[] = {
    3,  3, 3,  4, 3,  5, 3,  9, 3,  10, 3,  11, 3,  15, 3,  16, 3,  17,
    4,  3, 4,  4, 4,  5, 4,  9, 4,  10, 4,  11, 4,  15, 4,  16, 4,  17,
    5,  3, 5,  4, 5,  5, 5,  9, 5,  10, 5,  11, 5,  15, 5,  16, 5,  17,
    9,  3, 9,  4, 9,  5, 9,  9, 9,  10, 9,  11, 9,  15, 9,  16, 9,  17,
    10, 3, 10, 4, 10, 5, 10, 9, 10, 10, 10, 11, 10, 15, 10, 16, 10, 17,
    11, 3, 11, 4, 11, 5, 11, 9, 11, 10, 11, 11, 11, 15, 11, 16, 11, 17,
    15, 3, 15, 4, 15, 5, 15, 9, 15, 10, 15, 11, 15, 15, 15, 16, 15, 17,
    16, 3, 16, 4, 16, 5, 16, 9, 16, 10, 16, 11, 16, 15, 16, 16, 16, 17,
    17, 3, 17, 4, 17, 5, 17, 9, 17, 10, 17, 11, 17, 15, 17, 16, 17, 17};

int TACSNeohookean3D::getNumParameters() { return 3; }

int TACSNeohookean3D::getVarsPerNode() { return 3; }

/**
   Evaluate the coefficients of the weak form integrand
*/
void TACSNeohookean3D::evalWeakIntegrand(
    int elemIndex, const double time, int n, const double pt[],
    const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
    const TacsScalar Ux[], TacsScalar DUt[], TacsScalar DUx[]) {
  TacsScalar F[9];
  F[0] = 1.0 + Ux[0];
  F[1] = Ux[1];
  F[2] = Ux[2];

  F[3] = Ux[3];
  F[4] = 1.0 + Ux[4];
  F[5] = Ux[5];

  F[6] = Ux[6];
  F[7] = Ux[7];
  F[8] = 1.0 + Ux[8];

  // Compute tr(C) = tr(F^{T}*F) = sum_{ij} F_{ij}^2
  TacsScalar I1 =
      (F[0] * F[0] + F[1] * F[1] + F[2] * F[2] + F[3] * F[3] + F[4] * F[4] +
       F[5] * F[5] + F[6] * F[6] + F[7] * F[7] + F[8] * F[8]);
  TacsScalar J = det3x3(F);

  TacsScalar dI1, dJ;
  evalStrainEnergyDeriv(I1, J, &dI1, &dJ);

  memset(DUt, 0, 9 * sizeof(TacsScalar));
  memset(DUx, 0, 9 * sizeof(TacsScalar));

  // Add dU0/dJ*dJ/dUx
  addDet3x3Sens(dJ, F, DUx);

  // Add dU0/dI1*dI1/dUx
  DUx[0] += 2.0 * F[0] * dI1;
  DUx[1] += 2.0 * F[1] * dI1;
  DUx[2] += 2.0 * F[2] * dI1;

  DUx[3] += 2.0 * F[3] * dI1;
  DUx[4] += 2.0 * F[4] * dI1;
  DUx[5] += 2.0 * F[5] * dI1;

  DUx[6] += 2.0 * F[6] * dI1;
  DUx[7] += 2.0 * F[7] * dI1;
  DUx[8] += 2.0 * F[8] * dI1;
}

/**
  Get the non-zero pattern for the matrix
*/
void TACSNeohookean3D::getWeakMatrixNonzeros(ElementMatrixType matType,
                                             int elemIndex, int *Jac_nnz,
                                             const int *Jac_pairs[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    *Jac_nnz = 81;
    *Jac_pairs = linear_Jac_pairs;
  } else {
    *Jac_nnz = 0;
    *Jac_pairs = NULL;
  }
}

/**
  Evaluate the derivatives of the weak form coefficients
*/
void TACSNeohookean3D::evalWeakMatrix(ElementMatrixType matType, int elemIndex,
                                      const double time, int n,
                                      const double pt[], const TacsScalar X[],
                                      const TacsScalar Xd[],
                                      const TacsScalar Ut[],
                                      const TacsScalar Ux[], TacsScalar DUt[],
                                      TacsScalar DUx[], TacsScalar Jac[]) {
  if (matType == TACS_JACOBIAN_MATRIX) {
    TacsScalar F[9];
    F[0] = 1.0 + Ux[0];
    F[1] = Ux[1];
    F[2] = Ux[2];

    F[3] = Ux[3];
    F[4] = 1.0 + Ux[4];
    F[5] = Ux[5];

    F[6] = Ux[6];
    F[7] = Ux[7];
    F[8] = 1.0 + Ux[8];

    // Compute tr(C) = tr(F^{T}*F) = sum_{ij} F_{ij}^2
    TacsScalar I1 =
        (F[0] * F[0] + F[1] * F[1] + F[2] * F[2] + F[3] * F[3] + F[4] * F[4] +
         F[5] * F[5] + F[6] * F[6] + F[7] * F[7] + F[8] * F[8]);
    TacsScalar J = det3x3(F);

    TacsScalar dI1, dJ;
    evalStrainEnergyDeriv(I1, J, &dI1, &dJ);

    memset(DUt, 0, 9 * sizeof(TacsScalar));
    memset(DUx, 0, 9 * sizeof(TacsScalar));

    // Add dU0/dJ*dJ/dUx
    addDet3x3Sens(dJ, F, DUx);

    // Add dU0/dI1*dI1/dUx
    DUx[0] += 2.0 * F[0] * dI1;
    DUx[1] += 2.0 * F[1] * dI1;
    DUx[2] += 2.0 * F[2] * dI1;

    DUx[3] += 2.0 * F[3] * dI1;
    DUx[4] += 2.0 * F[4] * dI1;
    DUx[5] += 2.0 * F[5] * dI1;

    DUx[6] += 2.0 * F[6] * dI1;
    DUx[7] += 2.0 * F[7] * dI1;
    DUx[8] += 2.0 * F[8] * dI1;

    // Add the second derivative contributions from the first derivatives
    // of the strain energy
    det3x32ndSens(dJ, F, Jac);
    for (int i = 0; i < 9; i++) {
      Jac[10 * i] += 2.0 * dI1;
    }

    // Add the other second derivative contributions directly from the
    // strain energy
    TacsScalar d2I1, d2J, dIJ;
    evalStrainEnergy2ndDeriv(I1, J, &d2I1, &d2J, &dIJ);

    if (TacsRealPart(d2I1) != 0.0) {
      for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
          Jac[9 * i + j] += 4.0 * F[i] * F[j] * d2I1;
        }
      }
    }
    if (TacsRealPart(d2J) != 0.0) {
      TacsScalar t[9];
      det3x3Sens(F, t);

      for (int i = 0; i < 9; i++) {
        addDet3x3Sens(d2J * t[i], F, &Jac[9 * i]);
      }
    }
    if (TacsRealPart(dIJ) != 0.0) {
      for (int i = 0; i < 9; i++) {
        addDet3x3Sens(2.0 * dIJ * F[i], F, &Jac[9 * i]);
      }
    }
  }
}

/**
   Get the output for a single node in the mesh
*/
void TACSNeohookean3D::getOutputData(int elemIndex, const double time,
                                     ElementType etype, int write_flag,
                                     const double pt[], const TacsScalar X[],
                                     const TacsScalar Ut[],
                                     const TacsScalar Ux[], int ld_data,
                                     TacsScalar *data) {
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

    TacsScalar F[9];
    F[0] = 1.0 + Ux[0];
    F[1] = Ux[1];
    F[2] = Ux[2];

    F[3] = Ux[3];
    F[4] = 1.0 + Ux[4];
    F[5] = Ux[5];

    F[6] = Ux[6];
    F[7] = Ux[7];
    F[8] = 1.0 + Ux[8];

    TacsScalar C[9];
    mat3x3TransMatMult(F, F, C);

    TacsScalar e[6];
    e[0] = C[0];
    e[1] = C[4];
    e[2] = C[8];
    e[3] = C[5];
    e[4] = C[2];
    e[5] = C[1];

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
      data[0] = e[0];
      data[1] = e[1];
      data[2] = e[2];
      data[3] = e[3];
      data[4] = e[4];
      data[5] = e[5];
      data += 6;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS) {
      data[0] = data[1] = data[2] = data[3] = 0.0;
      data += 4;
    }
  }
}
