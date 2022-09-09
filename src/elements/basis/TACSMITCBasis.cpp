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

#include "TACSMITCBasis.h"

#include "TACSElementAlgebra.h"

void TACSMITCBasis::getTyingPointGradient(int field, int ty, const double pt[],
                                          const TacsScalar Xpts[],
                                          const int vars_per_node,
                                          const TacsScalar vars[],
                                          TacsScalar Xd[], TacsScalar U[],
                                          TacsScalar Ud[]) {
  interpFieldsGrad(-1, pt, 3, Xpts, Xd);
  interpFields(-1, pt, vars_per_node, vars, 1, U);
  interpFieldsGrad(-1, pt, vars_per_node, vars, Ud);
}

void TACSMITCBasis::interpTyingField(int n, const double pt[],
                                     const TacsScalar qty[], TacsScalar Uty[]) {
  const int num_ty_fields = getNumTyingFields();
  for (int field = 0; field < num_ty_fields; field++) {
    double N[MAX_NUM_TYING_POINTS];
    computeTyingBasis(field, pt, N);

    const int num_ty_points = getNumTyingPoints(field);
    Uty[field] = 0.0;
    for (int i = 0; i < num_ty_points; i++) {
      Uty[field] += qty[0] * N[i];
      qty++;
    }
  }
}

void TACSMITCBasis::addInterpTyingFieldTranspose(int n, const double pt[],
                                                 const TacsScalar Uty[],
                                                 TacsScalar qty[]) {
  const int num_ty_fields = getNumTyingFields();
  for (int field = 0; field < num_ty_fields; field++) {
    double N[MAX_NUM_TYING_POINTS];
    computeTyingBasis(field, pt, N);

    const int num_ty_points = getNumTyingPoints(field);
    for (int i = 0; i < num_ty_points; i++) {
      qty[0] += N[i] * Uty[field];
      qty++;
    }
  }
}

void TACSMITCBasis::addWeakResidual(int n, const double pt[], TacsScalar weight,
                                    const TacsScalar J[],
                                    const int vars_per_node, TacsScalar DUt[],
                                    TacsScalar DUx[], TacsScalar res[],
                                    TacsScalar rty[]) {
  const int num_params = getNumParameters();

  // Add contributions from DUt
  for (int i = 0; i < 3 * vars_per_node; i++) {
    DUt[i] *= weight;
  }
  addInterpFieldsTranspose(n, pt, 3, &DUt[0], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[1], vars_per_node, res);
  addInterpFieldsTranspose(n, pt, 3, &DUt[2], vars_per_node, res);

  if (num_params == 3) {
    for (int i = 0; i < vars_per_node; i++) {
      TacsScalar dx[3];
      mat3x3Mult(J, &DUx[3 * i], dx);
      DUx[3 * i] = weight * dx[0];
      DUx[3 * i + 1] = weight * dx[1];
      DUx[3 * i + 2] = weight * dx[2];
    }
  } else if (num_params == 2) {
    for (int i = 0; i < vars_per_node; i++) {
      TacsScalar dx[2];
      mat2x2Mult(J, &DUx[2 * i], dx);
      DUx[2 * i] = weight * dx[0];
      DUx[2 * i + 1] = weight * dx[1];
    }
  } else {
    for (int i = 0; i < vars_per_node; i++) {
      DUx[i] *= J[0] * weight;
    }
  }

  addInterpFieldsGradTranspose(n, pt, vars_per_node, DUx, res);

  // Add the contributions to the tying residual
  addInterpTyingFieldTranspose(n, pt, &DUt[3 * vars_per_node], rty);
}

/*
  Add the contributions to the residual from a tying point location
*/
void TACSMITCBasis::addTyingResidual(int field, int ty, const double pt[],
                                     const int vars_per_node,
                                     const TacsScalar DU[],
                                     const TacsScalar DUd[], TacsScalar res[]) {
  addInterpFieldsTranspose(-1, pt, 1, DU, vars_per_node, res);
  addInterpFieldsGradTranspose(-1, pt, vars_per_node, DUd, res);
}