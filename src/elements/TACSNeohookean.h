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

#ifndef TACS_NEOHOOKEAN_H
#define TACS_NEOHOOKEAN_H

#include <math.h>

#include "TACSElementModel.h"

class TACSNeohookean3D : public TACSElementModel {
 public:
  TACSNeohookean3D(TacsScalar C1, TacsScalar D1);
  ~TACSNeohookean3D();

  int getNumParameters();
  int getVarsPerNode();

  /**
    Evaluate the coefficients of the weak form integrand
  */
  void evalWeakIntegrand(int elemIndex, const double time, int n,
                         const double pt[], const TacsScalar X[],
                         const TacsScalar Xd[], const TacsScalar Ut[],
                         const TacsScalar Ux[], TacsScalar DUt[],
                         TacsScalar DUx[]);

  /**
    Get the non-zero pattern for the matrix
  */
  void getWeakMatrixNonzeros(ElementMatrixType matType, int elemIndex,
                             int *Jac_nnz, const int *Jac_pairs[]);

  /**
    Evaluate the derivatives of the weak form coefficients
  */
  void evalWeakMatrix(ElementMatrixType matType, int elemIndex,
                      const double time, int n, const double pt[],
                      const TacsScalar X[], const TacsScalar Xd[],
                      const TacsScalar Ut[], const TacsScalar Ux[],
                      TacsScalar DUt[], TacsScalar DUx[], TacsScalar Jac[]);

  /**
    Get the output for a single node in the mesh
  */
  void getOutputData(int elemIndex, const double time, ElementType etype,
                     int write_flag, const double pt[], const TacsScalar X[],
                     const TacsScalar Ut[], const TacsScalar Ux[], int ld_data,
                     TacsScalar *data);

 private:
  // Store the coefficients
  TacsScalar C1, D1;

  // Evaluate the strain energy and its derivative
  TacsScalar evalStrainEnergy(TacsScalar I1, TacsScalar J) {
    return C1 * (I1 - 3.0 - 2.0 * log(J)) + D1 * (J - 1.0) * (J - 1.0);
  }
  void evalStrainEnergyDeriv(TacsScalar I1, TacsScalar J, TacsScalar *dI1,
                             TacsScalar *dJ) {
    *dI1 = C1;
    *dJ = -2.0 * C1 / J + 2.0 * D1 * (J - 1.0);
  }
  void evalStrainEnergy2ndDeriv(TacsScalar I1, TacsScalar J, TacsScalar *dI1,
                                TacsScalar *dJ, TacsScalar *dIJ) {
    *dI1 = 0.0;
    *dJ = 2.0 * C1 / (J * J) + 2.0 * D1;
    *dIJ = 0.0;
  }

  // Constant member data
  static const int linear_Jac_pairs[2 * 81];
};

#endif  // TACS_NEOHOOKEAN
