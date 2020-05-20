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

#include "TACSMixedInterpElementBasis.h"

void TACSMixedInterpElementBasis::getTyingFieldValues( int n,
                                                       const double pt[],
                                                       const TacsScalar tyvars[],
                                                       TacsScalar U[] ){
  const int num_ty_fields = getNumTyingFields();
  for ( int field = 0; field < num_ty_fields; field++ ){
    double N[MAX_NUM_TYING_POINTS];
    computeTyingFieldBasis(field, n, pt, N);

    const double *n = N;
    const int num_ty_points = getNumTyingPoints(field);
    U[field] = 0.0;
    for ( int i = 0; i < num_ty_points; i++ ){
      U[field] += tyvars[0]*N[i];
      tyvars++;
    }
  }
}

void TACSMixedInterpElementBasis::addMixedWeakFormResidual( int n,
                                                            const double pt[],
                                                            TacsScalar weight,
                                                            const int vars_per_node,
                                                            const TacsScalar DUty[],
                                                            TacsScalar rty[] ){
  // Compute the shape functions for each tying field
  int num_ty_fields = getNumTyingFields();
  int num_ty_points[MAX_NUM_TYING_FIELDS];

  double Nty[MAX_TOTAL_TYING_POINTS];
  for ( int field = 0, index = 0; field < num_ty_fields; field++ ){
    num_ty_points[field] = getNumTyingPoints(field);
    computeTyingFieldBasis(field, n, pt, &Nty[index]);
    index += num_ty_points[field];
  }

  // Add the residual contributions to both
  addMixedWeakFormResidual(num_ty_fields, num_ty_points,
                           Nty, weight, vars_per_node, DUty, rty);
}

void TACSMixedInterpElementBasis::addMixedWeakFormResidual( const int num_ty_fields,
                                                            const int num_ty_points[],
                                                            const double Nty[],
                                                            TacsScalar weight,
                                                            const TacsScalar DUty[],
                                                            TacsScalar rty[] ){
  // Add contributions to the tying residual
  const TacsScalar *dty = DUt[num_params*vars_per_node];
  for ( int i = 0; i < num_ty_fields; i++ ){
    for ( int j = 0; j < num_ty_points[i]; j++ ){
      rty[0] += weight*Nty[0]*dty[i];
      rty++;
      Nty++;
    }
  }
}
