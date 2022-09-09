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

#include "TACSConstitutiveVerification.h"

#include "tacslapack.h"

int TacsTestConstitutiveDensity(TACSConstitutive *con, int elemIndex,
                                const double pt[], const TacsScalar X[],
                                int ndvs, const TacsScalar *dvs, double dh,
                                int test_print_level, double test_fail_atol,
                                double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Sensitivity w.r.t. the density
  TacsScalar rho = con->evalDensity(elemIndex, pt, X);

  // Compute the sensitivity w.r.t. the density
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar scale;
  TacsGenerateRandomArray(&scale, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addDensityDVSens(elemIndex, scale, pt, X, ndvs, dfdx);

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    TacsScalar rho_forward = con->evalDensity(elemIndex, pt, X);
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = scale * TacsImagPart(rho_forward) / dh;
#else
    dfdx_approx[i] = scale * (rho_forward - rho) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the density for %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(rho)/dx", dfdx, dfdx_approx, ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveSpecificHeat(TACSConstitutive *con, int elemIndex,
                                     const double pt[], const TacsScalar X[],
                                     int ndvs, const TacsScalar *dvs, double dh,
                                     int test_print_level,
                                     double test_fail_atol,
                                     double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the sensitivity w.r.t. the density
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar scale;
  TacsGenerateRandomArray(&scale, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addSpecificHeatDVSens(elemIndex, scale, pt, X, ndvs, dfdx);

  // Evaluate the specific heat
  TacsScalar c = con->evalSpecificHeat(elemIndex, pt, X);

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    TacsScalar c_forward = con->evalSpecificHeat(elemIndex, pt, X);
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = scale * TacsImagPart(c_forward) / dh;
#else
    dfdx_approx[i] = scale * (c_forward - c) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the specific heat for %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(c)/dx", dfdx, dfdx_approx, ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveHeatFlux(TACSConstitutive *con, int elemIndex,
                                 const double pt[], const TacsScalar X[],
                                 int ndvs, const TacsScalar *dvs, double dh,
                                 int test_print_level, double test_fail_atol,
                                 double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  TacsScalar flux[3], grad[3], psi[3];
  TacsGenerateRandomArray(grad, 3);
  TacsGenerateRandomArray(psi, 3);

  // Allocate space for the derivatives
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar scale;
  TacsGenerateRandomArray(&scale, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addHeatFluxDVSens(elemIndex, scale, pt, X, grad, psi, ndvs, dfdx);

  flux[0] = flux[1] = flux[2] = 0.0;
  con->evalHeatFlux(elemIndex, pt, X, grad, flux);
  TacsScalar prod = 0.0;
  for (int j = 0; j < 3; j++) {
    prod += psi[j] * flux[j];
  }

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    flux[0] = flux[1] = flux[2] = 0.0;
    con->evalHeatFlux(elemIndex, pt, X, grad, flux);
    TacsScalar prod_forward = 0.0;
    for (int j = 0; j < 3; j++) {
      prod_forward += psi[j] * flux[j];
    }
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = scale * TacsImagPart(prod_forward) / dh;
#else
    dfdx_approx[i] = scale * (prod_forward - prod) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the heat-flux product %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(psi^{T}flux)/dx", dfdx, dfdx_approx,
                             ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveStress(TACSConstitutive *con, int elemIndex,
                               const double pt[], const TacsScalar X[],
                               int ndvs, const TacsScalar *dvs, double dh,
                               int test_print_level, double test_fail_atol,
                               double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  int nstress = con->getNumStresses();
  TacsScalar *s = new TacsScalar[nstress];
  TacsScalar *e = new TacsScalar[nstress];
  TacsScalar *psi = new TacsScalar[nstress];
  TacsGenerateRandomArray(e, nstress);
  TacsGenerateRandomArray(psi, nstress);

  // Allocate space for the derivatives
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar scale;
  TacsGenerateRandomArray(&scale, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addStressDVSens(elemIndex, scale, pt, X, e, psi, ndvs, dfdx);

  con->evalStress(elemIndex, pt, X, e, s);
  TacsScalar prod = 0.0;
  for (int j = 0; j < nstress; j++) {
    prod += psi[j] * s[j];
  }

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    con->evalStress(elemIndex, pt, X, e, s);
    TacsScalar prod_forward = 0.0;
    for (int j = 0; j < nstress; j++) {
      prod_forward += psi[j] * s[j];
    }
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = scale * TacsImagPart(prod_forward) / dh;
#else
    dfdx_approx[i] = scale * (prod_forward - prod) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the stress-adjoint product %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(psi^{T}s)/dx", dfdx, dfdx_approx, ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;
  delete[] s;
  delete[] e;
  delete[] psi;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveThermalStrain(TACSConstitutive *con, int elemIndex,
                                      const double pt[], const TacsScalar X[],
                                      int ndvs, const TacsScalar *dvs,
                                      double dh, int test_print_level,
                                      double test_fail_atol,
                                      double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  int nstress = con->getNumStresses();
  TacsScalar *e = new TacsScalar[nstress];
  TacsScalar *psi = new TacsScalar[nstress];
  TacsGenerateRandomArray(e, nstress);
  TacsGenerateRandomArray(psi, nstress);

  // Allocate space for the derivatives
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar theta;
  TacsGenerateRandomArray(&theta, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addThermalStrainDVSens(elemIndex, pt, X, theta, psi, ndvs, dfdx);

  con->evalThermalStrain(elemIndex, pt, X, theta, e);
  TacsScalar prod = 0.0;
  for (int j = 0; j < nstress; j++) {
    prod += psi[j] * e[j];
  }

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    con->evalThermalStrain(elemIndex, pt, X, theta, e);
    TacsScalar prod_forward = 0.0;
    for (int j = 0; j < nstress; j++) {
      prod_forward += psi[j] * e[j];
    }
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = TacsImagPart(prod_forward) / dh;
#else
    dfdx_approx[i] = (prod_forward - prod) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the thermal strain-adjoint product %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(psi^{T}et)/dx", dfdx, dfdx_approx,
                             ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;
  delete[] e;
  delete[] psi;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveFailure(TACSConstitutive *con, int elemIndex,
                                const double pt[], const TacsScalar X[],
                                int ndvs, const TacsScalar *dvs, double dh,
                                int test_print_level, double test_fail_atol,
                                double test_fail_rtol) {
  con->setDesignVars(elemIndex, ndvs, dvs);

  int nstress = con->getNumStresses();
  TacsScalar *e = new TacsScalar[nstress];
  TacsGenerateRandomArray(e, nstress);

  // Allocate space for the derivatives
  TacsScalar *xtemp = new TacsScalar[ndvs];
  TacsScalar *dfdx = new TacsScalar[ndvs];
  TacsScalar *dfdx_approx = new TacsScalar[ndvs];

  // Copy the design variable values
  memcpy(xtemp, dvs, ndvs * sizeof(TacsScalar));

  // Generate a random array for the scale factor
  TacsScalar scale;
  TacsGenerateRandomArray(&scale, 1);

  // Evaluate the derivatives
  memset(dfdx, 0, ndvs * sizeof(TacsScalar));
  con->addFailureDVSens(elemIndex, scale, pt, X, e, ndvs, dfdx);

  TacsScalar failure = con->evalFailure(elemIndex, pt, X, e);

  // Compute the approximate derivative
  for (int i = 0; i < ndvs; i++) {
    TacsScalar x = xtemp[i];
#ifdef TACS_USE_COMPLEX
    xtemp[i] = x + TacsScalar(0.0, dh);
#else
    xtemp[i] = x + dh;
#endif
    con->setDesignVars(elemIndex, ndvs, xtemp);

    TacsScalar failure_forward = con->evalFailure(elemIndex, pt, X, e);
#ifdef TACS_USE_COMPLEX
    dfdx_approx[i] = scale * TacsImagPart(failure_forward) / dh;
#else
    dfdx_approx[i] = scale * (failure_forward - failure) / dh;
#endif
    xtemp[i] = x;
  }

  // Reset the design variable values
  con->setDesignVars(elemIndex, ndvs, dvs);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfdx, dfdx_approx, ndvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(dfdx, dfdx_approx, ndvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the failure index for object %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(failure)/dx", dfdx, dfdx_approx, ndvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xtemp;
  delete[] dfdx;
  delete[] dfdx_approx;
  delete[] e;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestConstitutiveFailureStrainSens(TACSConstitutive *con, int elemIndex,
                                          const double pt[],
                                          const TacsScalar X[], double dh,
                                          int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol) {
  int nstress = con->getNumStresses();
  TacsScalar *e = new TacsScalar[nstress];
  TacsGenerateRandomArray(e, nstress);

  // Allocate space for the derivatives
  TacsScalar *dfde = new TacsScalar[nstress];
  TacsScalar *dfde_approx = new TacsScalar[nstress];

  // Evaluate the derivatives
  con->evalFailureStrainSens(elemIndex, pt, X, e, dfde);

  TacsScalar failure = con->evalFailure(elemIndex, pt, X, e);

  // Compute the approximate derivative
  for (int i = 0; i < nstress; i++) {
    TacsScalar et = e[i];
#ifdef TACS_USE_COMPLEX
    e[i] = et + TacsScalar(0.0, dh);
#else
    e[i] = et + dh;
#endif

    TacsScalar failure_forward = con->evalFailure(elemIndex, pt, X, e);
#ifdef TACS_USE_COMPLEX
    dfde_approx[i] = TacsImagPart(failure_forward) / dh;
#else
    dfde_approx[i] = (failure_forward - failure) / dh;
#endif
    e[i] = et;
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(dfde, dfde_approx, nstress, &max_err_index);
  double max_rel =
      TacsGetMaxRelError(dfde, dfde_approx, nstress, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the failure index w.r.t. strain for "
            "object %s\n",
            con->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d(failure)/de", dfde, dfde_approx,
                             nstress);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] e;
  delete[] dfde;
  delete[] dfde_approx;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestConstitutive(TACSConstitutive *con, int elemIndex, double dh,
                         int test_print_level, double test_fail_atol,
                         double test_fail_rtol) {
  // Retrieve the number of variables
  int ndvs = con->getDesignVarNums(elemIndex, 0, NULL);
  ndvs *= con->getDesignVarsPerNode();

  // Test the sensitivities w.r.t. the design variables first...
  TacsScalar *dvs = new TacsScalar[ndvs];
  con->getDesignVars(elemIndex, ndvs, dvs);

  srand(time(NULL));

  double pt[3];
  TacsGenerateRandomArray(pt, 3);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing derivatives for the constitutive object %s at "
            "parametric point: \npt[0]  %15.6e\npt[1]  %15.6e\npt[2]"
            "  %15.6e\n",
            con->getObjectName(), pt[0], pt[1], pt[2]);
  }

  TacsScalar X[3];
  TacsGenerateRandomArray(X, 3);

  int flag, fail = 0;
  flag = TacsTestConstitutiveDensity(con, elemIndex, pt, X, ndvs, dvs, dh,
                                     test_print_level, test_fail_atol,
                                     test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveSpecificHeat(con, elemIndex, pt, X, ndvs, dvs, dh,
                                          test_print_level, test_fail_atol,
                                          test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveHeatFlux(con, elemIndex, pt, X, ndvs, dvs, dh,
                                      test_print_level, test_fail_atol,
                                      test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveStress(con, elemIndex, pt, X, ndvs, dvs, dh,
                                    test_print_level, test_fail_atol,
                                    test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveThermalStrain(con, elemIndex, pt, X, ndvs, dvs, dh,
                                           test_print_level, test_fail_atol,
                                           test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveFailure(con, elemIndex, pt, X, ndvs, dvs, dh,
                                     test_print_level, test_fail_atol,
                                     test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestConstitutiveFailureStrainSens(con, elemIndex, pt, X, dh,
                                               test_print_level, test_fail_atol,
                                               test_fail_rtol);
  fail = flag || fail;

  return fail;
}
