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

#include "TACSConstitutive.h"

#include <math.h>
#include <stdlib.h>

#include "tacslapack.h"

const char *TACSConstitutive::constName = "TACSConstitutive";

/**
  Return the generic constitutive class name
*/
const char *TACSConstitutive::getObjectName() { return constName; }

/*
  Compute a two-dimensional representation of the failure envelope.
  Store the values in the output x_vals and y_vals.

  input:
  npts:       the number of points to include in the envelope plot
  pt:         the quadrature point to evaluate the failure criteria
  X:          the X coordinate of the point
  x_stress:   the x components of the stress
  y_stress:   the y components of the stress
  x_vals:     the x values of the failure abscissa
  y_vals:     the y values of the failure abscissa
*/
void TACSConstitutive::getFailureEnvelope(
    int npts, int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar x_stress[], const TacsScalar y_stress[],
    TacsScalar x_vals[], TacsScalar y_vals[]) {
  // The failure criteria solution tolerance
  double tol = 1e-8;

  // The max number of Newton iterations to use
  int max_newton_iters = 50;

  int nstress = getNumStresses();

  // Compute an explicit form of the constitutive relationship
  int *ipiv = new int[nstress];
  TacsScalar *e = new TacsScalar[nstress];
  TacsScalar *f_sens = new TacsScalar[nstress];
  TacsScalar *x_strain = new TacsScalar[nstress];
  TacsScalar *y_strain = new TacsScalar[nstress];
  TacsScalar *C = new TacsScalar[nstress * nstress];

  for (int i = 0; i < nstress; i++) {
    memset(e, 0, nstress * sizeof(TacsScalar));
    e[i] = 1.0;
    evalStress(elemIndex, pt, X, e, &C[i * nstress]);
  }

  // Factor the constitutive matrix
  int info;
  LAPACKgetrf(&nstress, &nstress, C, &nstress, ipiv, &info);

  // Copy over the stress values and solve for the strain
  memcpy(x_strain, x_stress, nstress * sizeof(TacsScalar));
  memcpy(y_strain, y_stress, nstress * sizeof(TacsScalar));

  int one = 1;
  LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, x_strain, &nstress,
              &info);
  LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, y_strain, &nstress,
              &info);

  TacsScalar P = 1.0;
  for (int k = 0; k < npts; k++) {
    TacsScalar theta = (2.0 * M_PI * k) / (npts - 1);
    TacsScalar c = cos(theta), s = sin(theta);

    // Solve the equation:
    // failure(pt, P*(c*x_strain + s*y_strain)) - 1.0 = 0
    for (int i = 0; i < max_newton_iters; i++) {
      for (int j = 0; j < nstress; j++) {
        e[j] = P * (c * x_strain[j] + s * y_strain[j]);
      }

      // Compute the failure criterion and the derivative of the
      // failure criterion w.r.t. the load parameter P
      TacsScalar fail, failSens;
      fail = evalFailure(elemIndex, pt, X, e);

      if (fabs(TacsRealPart(fail) - 1.0) < tol) {
        break;
      }

      evalFailureStrainSens(elemIndex, pt, X, e, f_sens);

      failSens = 0.0;
      for (int j = 0; j < nstress; j++) {
        failSens += f_sens[j] * (c * x_strain[j] + s * y_strain[j]);
      }

      // Compute the Newton update to the boundary
      P = P - (fail - 1.0) / failSens;
    }

    // Compute the strain at the final point
    for (int j = 0; j < nstress; j++) {
      e[j] = P * (c * x_strain[j] + s * y_strain[j]);
    }

    // Compute the corresponding stress at the final point
    evalStress(elemIndex, pt, X, e, f_sens);

    // Set the values of the output
    x_vals[k] = P * c;
    y_vals[k] = P * s;
  }

  delete[] ipiv;
  delete[] e;
  delete[] f_sens;
  delete[] x_strain;
  delete[] y_strain;
  delete[] C;
}
