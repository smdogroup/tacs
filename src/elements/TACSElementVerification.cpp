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

#include "TACSElementVerification.h"

#include "TACSElementAlgebra.h"
#include "tacslapack.h"

/*
  Seed random number generator for repeatable test results
*/
void TacsSeedRandomGenerator(int seed) { srand(seed); }

/*
  Assign variables randomly to an array. This is useful for
  testing various things.
*/
void TacsGenerateRandomArray(TacsReal *array, int size, TacsReal lower,
                             TacsReal upper) {
  for (int i = 0; i < size; i++) {
    array[i] = (upper - lower) * (rand() / ((double)RAND_MAX + 1)) + lower;
  }
}

void TacsGenerateRandomArray(TacsComplex *array, int size, TacsComplex lower,
                             TacsComplex upper) {
  for (int i = 0; i < size; i++) {
    array[i] = (upper - lower) * (rand() / ((double)RAND_MAX + 1)) + lower;
  }
}

/*
  Find the largest absolute value of the difference between the
  arrays a and b
*/
double TacsGetMaxError(TacsScalar *a, TacsScalar *b, int size, int *max_index) {
  double max_error = 0.0;
  *max_index = -1;

  for (int i = 0; i < size; i++) {
    double er = fabs(TacsRealPart(a[i] - b[i]));
    if (i == 0 || er > max_error) {
      max_error = er;
      *max_index = i;
    }
  }
  return max_error;
}

/*
  Find the maximum relative error between a and b and return the
*/
double TacsGetMaxRelError(TacsScalar *a, TacsScalar *b, int size,
                          int *max_index) {
  double max_error = 0.0;
  *max_index = -1;

  for (int i = 0; i < size; i++) {
    double er = 0.0;
    if (a[i] != 0.0) {
      er = fabs(TacsRealPart((a[i] - b[i]) / a[i]));
    }
    if (i == 0 || er > max_error) {
      max_error = er;
      *max_index = i;
    }
  }
  return max_error;
}

/*
  Print out the values and the relative errors
*/
void TacsPrintErrorComponents(FILE *fp, const char *descript, TacsScalar *a,
                              TacsScalar *b, int size) {
  fprintf(fp, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
          "Analytic", "Approximate", "Rel. Error");
  for (int i = 0; i < size; i++) {
    if (a[i] != 0.0) {
      fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, i,
              TacsRealPart(a[i]), TacsRealPart(b[i]),
              fabs(TacsRealPart((a[i] - b[i]) / a[i])));
    } else {
      fprintf(fp, "%s[%3d] %15.6e %15.6e\n", descript, i, TacsRealPart(a[i]),
              TacsRealPart(b[i]));
    }
  }
}

/*
  Perturb the input variables in the forward sense
*/
void TacsForwardDiffPerturb(TacsScalar *out, int size, const TacsScalar *orig,
                            const TacsScalar *pert, double dh) {
#ifdef TACS_USE_COMPLEX
  for (int k = 0; k < size; k++) {
    out[k] = orig[k] + TacsScalar(0.0, dh) * pert[k];
  }
#else
  for (int k = 0; k < size; k++) {
    out[k] = orig[k] + dh * pert[k];
  }
#endif  // TACS_USE_COMPLEX
}

/*
  Perturb the variables in the backward sense
*/
void TacsBackwardDiffPerturb(TacsScalar *out, int size, const TacsScalar *orig,
                             const TacsScalar *pert, double dh) {
#ifdef TACS_USE_COMPLEX
  for (int k = 0; k < size; k++) {
    out[k] = orig[k];
  }
#else
  for (int k = 0; k < size; k++) {
    out[k] = orig[k] - dh * pert[k];
  }
#endif  // TACS_USE_COMPLEX
}

/*
  Form the forward approximation
*/
void TacsFormDiffApproximate(TacsScalar *forward, const TacsScalar *backward,
                             int size, TacsScalar dh) {
#ifdef TACS_USE_COMPLEX
  for (int k = 0; k < size; k++) {
    forward[k] = TacsImagPart(forward[k]) / dh;
  }
#else
  for (int k = 0; k < size; k++) {
    forward[k] = (forward[k] - backward[k]) / (2.0 * dh);
  }
#endif  // TACS_USE_COMPLEX
}

/*
  The following function tests the consistency of the implementation
  of the residuals and the energy expressions, relying on Lagrange's
  equations.

  NOTE: This function only works when the Lagrange's equations do not
  have any constraints. Element-specific code is needed when Largrange
  multipliers are embedded within the element variables.

  This function uses finite-differences to compute the derivatives
  within Lagrange's equations and compares the result with the
  residual computed using the residual routine.

  Lagrange's equations of motion are given as follows:

  d/dt(dL/d(dot{q})^{T}) - dL/dq^{T} = 0

  This can be evaluated using finite-differencing as follows:

  dL/dqi(q, dq) .= (L(q, dq + h*ei) - L(q, dq - h*ei))/h

  d(f(q, dq))/dt .=
  (f(q + dt*dq, dq + dt*ddq) - f(q - dt*dq, dq - dt*ddq))/dt
*/
int TacsTestElementResidual(TACSElement *element, int elemIndex, double time,
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            double dh, int test_print_level,
                            double test_fail_atol, double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  // Allocate temporary arrays for the computation
  TacsScalar *q = new TacsScalar[nvars];
  TacsScalar *dq = new TacsScalar[nvars];

  // Temporary storage for the residuals
  TacsScalar *res1 = new TacsScalar[nvars];
  TacsScalar *res2 = new TacsScalar[nvars];
  TacsScalar *fd = new TacsScalar[nvars];

  // Compute the values of the variables at (t + dt)
  for (int i = 0; i < nvars; i++) {
    q[i] = vars[i] + dh * dvars[i];
    dq[i] = dvars[i] + dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for (int i = 0; i < nvars; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
#ifdef TACS_USE_COMPLEX
    dq[i] = dqtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);
    res1[i] = TacsImagPart((T1 - P1)) / dh;
#else
    dq[i] = dqtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif
    dq[i] = dqtmp;
  }

  // Compute the values of the variables at (t - dt)
  for (int i = 0; i < nvars; i++) {
    q[i] = vars[i] - dh * dvars[i];
    dq[i] = dvars[i] - dh * ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for (int i = 0; i < nvars; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
#ifdef TACS_USE_COMPLEX
    dq[i] = dqtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);
    res2[i] = TacsImagPart((T1 - P1)) / dh;
#else
    dq[i] = dqtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);

    // Compute and store the approximation
    res2[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif
    dq[i] = dqtmp;
  }

  // Evaluate the finite-difference for the first term in Largrange's
  // equations of motion
  for (int i = 0; i < nvars; i++) {
    fd[i] = 0.5 * (res1[i] - res2[i]) / dh;
  }

  // Reset the values of q and dq at time t
  for (int i = 0; i < nvars; i++) {
    q[i] = vars[i];
    dq[i] = dvars[i];
  }

  // Compute the contribution from dL/dq^{T}
  for (int i = 0; i < nvars; i++) {
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar qtmp = q[i];

#ifdef TACS_USE_COMPLEX
    q[i] = qtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);
    res1[i] = TacsImagPart((T1 - P1)) / dh;
#else
    q[i] = qtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    q[i] = qtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);

    // Compute and store the approximation
    res1[i] = 0.5 * ((T1 - P1) - (T2 - P2)) / dh;
#endif
    q[i] = qtmp;
  }

  // Add the result to the finite-difference result
  for (int i = 0; i < nvars; i++) {
    fd[i] -= res1[i];
  }

  // Evaluate the residual using the code
  memset(res1, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res1);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(res1, fd, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(res1, fd, nvars, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the residual implementation for element %s.\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1) {
    fprintf(stderr, "The difference between the FD and true residual is:\n");
    TacsPrintErrorComponents(stderr, "Res error", res1, fd, nvars);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] q;
  delete[] dq;
  delete[] res1;
  delete[] res2;
  delete[] fd;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  The following function tests the consistency between the
  implementation of the residuals and the implementation of the system
  Jacobian.

  input:
  col:   test only the specified column of the matrix
*/
int TacsTestElementJacobian(TACSElement *element, int elemIndex, double time,
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            int col, double dh, int test_print_level,
                            double test_fail_atol, double test_fail_rtol) {
  // Set the failure flag
  int fail = 0;

  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  TacsScalar *result = new TacsScalar[nvars];
  TacsScalar *temp = new TacsScalar[nvars];
  TacsScalar *pert = new TacsScalar[nvars];

  TacsScalar *q = new TacsScalar[nvars];
  TacsScalar *dq = new TacsScalar[nvars];
  TacsScalar *ddq = new TacsScalar[nvars];
  TacsScalar *res = new TacsScalar[nvars];
  TacsScalar *mat = new TacsScalar[nvars * nvars];

  if (col >= 0 && col < nvars) {
    memset(pert, 0, nvars * sizeof(TacsScalar));
    pert[col] = 1.0;
  } else {
    TacsGenerateRandomArray(pert, nvars);
  }

  // Compute the Jacobian
  double alpha = (1.0 * rand()) / RAND_MAX;
  double beta = (1.0 * rand()) / RAND_MAX;
  double gamma = (1.0 * rand()) / RAND_MAX;

  memset(mat, 0, nvars * nvars * sizeof(TacsScalar));
  element->addJacobian(elemIndex, time, alpha, beta, gamma, Xpts, vars, dvars,
                       ddvars, res, mat);

  // Evaluate the Jacobian
  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("T", &nvars, &nvars, &a, mat, &nvars, pert, &one, &b, result, &one);

  // Perturb the variables in the forward sense
  TacsForwardDiffPerturb(q, nvars, vars, pert, alpha * dh);
  TacsForwardDiffPerturb(dq, nvars, dvars, pert, beta * dh);
  TacsForwardDiffPerturb(ddq, nvars, ddvars, pert, gamma * dh);
  memset(res, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, q, dq, ddq, res);

  // Perturb the variables in the backward sens
  TacsBackwardDiffPerturb(q, nvars, vars, pert, alpha * dh);
  TacsBackwardDiffPerturb(dq, nvars, dvars, pert, beta * dh);
  TacsBackwardDiffPerturb(ddq, nvars, ddvars, pert, gamma * dh);
  memset(temp, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, q, dq, ddq, temp);

  // Form the FD/CS approximate
  TacsFormDiffApproximate(res, temp, nvars, dh);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, res, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, res, nvars, &max_rel_index);

  if (test_print_level > 0) {
    if (col >= 0 && col < nvars) {
      fprintf(stderr,
              "Testing column %d of the Jacobian matrix for element %s.\n", col,
              element->getObjectName());
    } else {
      fprintf(stderr, "Testing the Jacobian matrix for element %s.\n",
              element->getObjectName());
    }
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    if (col >= 0 && col < nvars) {
      fprintf(stderr, "Column %d of the Jacobian matrix is\n", col);
    } else {
      fprintf(stderr,
              "The product of a random vector and the Jacobian matrix is\n");
    }
    TacsPrintErrorComponents(stderr, "K*u", result, res, nvars);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Test residual computed through addJacobian, make sure it's consistent with
  // addResidual
  memset(temp, 0, nvars * sizeof(TacsScalar));
  memset(res, 0, nvars * sizeof(TacsScalar));

  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, temp);
  element->addJacobian(elemIndex, time, alpha, beta, gamma, Xpts, vars, dvars,
                       ddvars, res, mat);

  // Compute the error
  max_err = TacsGetMaxError(temp, res, nvars, &max_err_index);
  max_rel = TacsGetMaxRelError(temp, res, nvars, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(
        stderr,
        "Testing residual consistency of addJacobian method for element %s.\n",
        element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    fprintf(
        stderr,
        "The difference between addResidual and addJacobian residuals is\n");
    TacsPrintErrorComponents(stderr, "res", temp, res, nvars);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  delete[] temp;
  delete[] q;
  delete[] dq;
  delete[] ddq;
  delete[] mat;
  delete[] pert;
  delete[] res;
  delete[] result;

  return fail;
}

/*
  The following function tests the consistency between the
  implementation of the residuals and the implementation of the system
  Jacobian.

  input:
  col:   test only the specified column of the matrix
*/
int TacsTestElementMatFreeJacobian(
    TACSElement *element, int elemIndex, double time, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int col, double dh, int test_print_level,
    double test_fail_atol, double test_fail_rtol) {
  // Set the failure flag
  int fail = 0;

  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  TacsScalar *result = new TacsScalar[nvars];
  TacsScalar *pert = new TacsScalar[nvars];
  TacsScalar *res = new TacsScalar[nvars];
  TacsScalar *mat = new TacsScalar[nvars * nvars];

  if (col >= 0 && col < nvars) {
    memset(pert, 0, nvars * sizeof(TacsScalar));
    pert[col] = 1.0;
  } else {
    TacsGenerateRandomArray(pert, nvars);
  }

  // Compute the Jacobian
  double alpha = (1.0 * rand()) / RAND_MAX;
  double beta = (1.0 * rand()) / RAND_MAX;
  double gamma = (1.0 * rand()) / RAND_MAX;

  memset(mat, 0, nvars * nvars * sizeof(TacsScalar));
  element->addJacobian(elemIndex, time, alpha, beta, gamma, Xpts, vars, dvars,
                       ddvars, res, mat);

  // Evaluate the Jacobian
  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("T", &nvars, &nvars, &a, mat, &nvars, pert, &one, &b, result, &one);

  // Get the size of the matrix
  int dsize, tsize;
  element->getMatVecDataSizes(TACS_JACOBIAN_MATRIX, elemIndex, &dsize, &tsize);

  // Allocate the matrix
  TacsScalar *tarray = new TacsScalar[tsize];
  TacsScalar *data = new TacsScalar[dsize];
  element->getMatVecProductData(TACS_JACOBIAN_MATRIX, elemIndex, time, alpha,
                                beta, gamma, Xpts, vars, dvars, ddvars, data);

  memset(res, 0, nvars * sizeof(TacsScalar));
  element->addMatVecProduct(TACS_JACOBIAN_MATRIX, elemIndex, data, tarray, pert,
                            res);

  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, res, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, res, nvars, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the matrix-free Jacobian for element %s.\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    fprintf(stderr,
            "The product of a random vector and the Jacobian matrix is\n");
    TacsPrintErrorComponents(stderr, "mat-free", result, res, nvars);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] pert;
  delete[] res;
  delete[] result;
  delete[] mat;
  delete[] data;
  delete[] tarray;

  return fail;
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestAdjResProduct(TACSElement *element, int elemIndex, double time,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          int dvLen, const TacsScalar *x, double dh,
                          int test_print_level, double test_fail_atol,
                          double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();
  int dvs_per_node = element->getDesignVarsPerNode();
  int num_dvs = dvs_per_node * dvLen;

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[num_dvs];
  memset(result, 0, num_dvs * sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[nvars];
  TacsGenerateRandomArray(adjoint, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0 * rand() / RAND_MAX;

  element->setDesignVars(elemIndex, dvLen, x);
  element->addAdjResProduct(elemIndex, time, scale, adjoint, Xpts, vars, dvars,
                            ddvars, dvLen, result);

  // Compute the product of the result with a perturbation
  // vector that is equal to perturb = sign(result[k])
  TacsScalar dpdx = 0.0;
  for (int k = 0; k < num_dvs; k++) {
    dpdx += fabs(result[k]);
  }

  // Allocate an array to store the perturbed design variable
  // values
  TacsScalar *xpert = new TacsScalar[num_dvs];
  TacsScalar fd_dpdx = 0.0;

  // Zero the residual
  TacsScalar *res = new TacsScalar[nvars];

#ifdef TACS_USE_COMPLEX
  // Perturb the design variables: xpert = x + dh*sign(result[k])
  for (int k = 0; k < num_dvs; k++) {
    if (TacsRealPart(result[k]) >= 0.0) {
      xpert[k] = x[k] + TacsScalar(0.0, dh);
    } else {
      xpert[k] = x[k] - TacsScalar(0.0, dh);
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  TacsScalar p1 = 0.0;
  memset(res, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  for (int k = 0; k < nvars; k++) {
    p1 += scale * res[k] * adjoint[k];
  }

  fd_dpdx = TacsImagPart(p1) / dh;
#else
  // Perturb the design variables: xpert = x + dh*sign(result[k])
  for (int k = 0; k < num_dvs; k++) {
    if (result[k] >= 0.0) {
      xpert[k] = x[k] + dh;
    } else {
      xpert[k] = x[k] - dh;
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  memset(res, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  TacsScalar p1 = 0.0;
  for (int k = 0; k < nvars; k++) {
    p1 += scale * res[k] * adjoint[k];
  }

  // Pertub the design variables: xpert = x - dh*sign(result[k])
  for (int k = 0; k < num_dvs; k++) {
    if (result[k] >= 0.0) {
      xpert[k] = x[k] - dh;
    } else {
      xpert[k] = x[k] + dh;
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  // Compute the residual again
  memset(res, 0, nvars * sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  TacsScalar p2 = 0.0;
  for (int k = 0; k < nvars; k++) {
    p2 += scale * res[k] * adjoint[k];
  }

  // Compute the finite-difference approximation
  fd_dpdx = 0.5 * (p1 - p2) / dh;
#endif

  // Set the design variable values
  element->setDesignVars(elemIndex, dvLen, x);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(&dpdx, &fd_dpdx, 1, &max_err_index);
  double max_rel = TacsGetMaxRelError(&dpdx, &fd_dpdx, 1, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the adjoint-residual product for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "Adj-Res product", &dpdx, &fd_dpdx, 1);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] adjoint;
  delete[] xpert;
  delete[] res;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestAdjResXptProduct(TACSElement *element, int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[], double dh,
                             int test_print_level, double test_fail_atol,
                             double test_fail_rtol) {
  int nvars = element->getNumNodes() * element->getVarsPerNode();
  int nnodes = element->getNumNodes();

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[3 * nnodes];
  memset(result, 0, 3 * nnodes * sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[nvars];
  TacsGenerateRandomArray(adjoint, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0 * rand() / RAND_MAX;
  element->addAdjResXptProduct(elemIndex, time, scale, adjoint, Xpts, vars,
                               dvars, ddvars, result);

  // Allocate space to store the results
  TacsScalar *fd = new TacsScalar[3 * nnodes];
  TacsScalar *X = new TacsScalar[3 * nnodes];
  TacsScalar *res = new TacsScalar[nvars];

  for (int k = 0; k < 3 * nnodes; k++) {
    // Copy the points
    memcpy(X, Xpts, 3 * nnodes * sizeof(TacsScalar));

    // Perturb the nodes in the forward sense
    TacsScalar one = 1.0;
    TacsForwardDiffPerturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars * sizeof(TacsScalar));
    element->addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
    TacsScalar p1 = 0.0;
    for (int i = 0; i < nvars; i++) {
      p1 += scale * adjoint[i] * res[i];
    }

    // Perturb the nodes in the reverse sense
    TacsBackwardDiffPerturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars * sizeof(TacsScalar));
    element->addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
    TacsScalar p2 = 0.0;
    for (int i = 0; i < nvars; i++) {
      p2 += scale * adjoint[i] * res[i];
    }

    // Form the approximation
    TacsFormDiffApproximate(&p1, &p2, 1, dh);

    // Set the
    fd[k] = p1;
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, 3 * nnodes, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, 3 * nnodes, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the adjoint-residual product for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "Adj-Res Xpt product", result, fd,
                             3 * nnodes);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] adjoint;
  delete[] fd;
  delete[] X;
  delete[] res;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestElementMatDVSens(TACSElement *element, ElementMatrixType matType,
                             int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             int dvLen, const TacsScalar *x, double dh,
                             int test_print_level, double test_fail_atol,
                             double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();
  int dvs_per_node = element->getDesignVarsPerNode();
  int num_dvs = dvs_per_node * dvLen;

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *xcopy = new TacsScalar[num_dvs];
  TacsScalar *result = new TacsScalar[num_dvs];
  TacsScalar *fd = new TacsScalar[num_dvs];
  memset(result, 0, num_dvs * sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *mat = new TacsScalar[nvars * nvars];
  TacsScalar *psi = new TacsScalar[nvars];
  TacsScalar *phi = new TacsScalar[nvars];
  TacsGenerateRandomArray(psi, nvars);
  TacsGenerateRandomArray(phi, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0 * rand() / RAND_MAX;

  element->setDesignVars(elemIndex, dvLen, x);
  memset(result, 0, num_dvs * sizeof(TacsScalar));
  element->addMatDVSensInnerProduct(matType, elemIndex, time, scale, psi, phi,
                                    Xpts, vars, dvLen, result);

  for (int k = 0; k < num_dvs; k++) {
    memcpy(xcopy, x, num_dvs * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    // Perturb the design variables: xpert = x + dh*sign(result[k])
    xcopy[k] = x[k] + TacsScalar(0.0, dh);
    element->setDesignVars(elemIndex, dvLen, xcopy);

    TacsScalar p1 = 0.0;
    element->getMatType(matType, elemIndex, time, Xpts, vars, mat);
    for (int i = 0; i < nvars; i++) {
      for (int j = 0; j < nvars; j++) {
        p1 += scale * mat[i + j * nvars] * psi[i] * phi[j];
      }
    }

    fd[k] = TacsImagPart(p1) / dh;
#else
    xcopy[k] = x[k] + dh;
    element->setDesignVars(elemIndex, dvLen, xcopy);

    TacsScalar p1 = 0.0;
    element->getMatType(matType, elemIndex, time, Xpts, vars, mat);
    for (int i = 0; i < nvars; i++) {
      for (int j = 0; j < nvars; j++) {
        p1 += scale * mat[i + j * nvars] * psi[i] * phi[j];
      }
    }

    xcopy[k] = x[k] - dh;
    element->setDesignVars(elemIndex, dvLen, xcopy);

    TacsScalar p2 = 0.0;
    element->getMatType(matType, elemIndex, time, Xpts, vars, mat);
    for (int i = 0; i < nvars; i++) {
      for (int j = 0; j < nvars; j++) {
        p2 += scale * mat[i + j * nvars] * psi[i] * phi[j];
      }
    }

    fd[k] = 0.5 * (p1 - p2) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Set the design variable values
  element->setDesignVars(elemIndex, dvLen, x);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, num_dvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, num_dvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr,
            "Testing the derivative of the matrix inner product for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "Matrix inner product", result, fd,
                             num_dvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] xcopy;
  delete[] psi;
  delete[] phi;
  delete[] mat;
  delete[] result;
  delete[] fd;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestElementMatXptSens(TACSElement *element, ElementMatrixType elemType,
                              int elemIndex, double time,
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              double dh, int test_print_level,
                              double test_fail_atol, double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  TacsScalar *psi = new TacsScalar[nvars];
  TacsScalar *phi = new TacsScalar[nvars];

  int nnodes = element->getNumNodes();

  TacsScalar *result = new TacsScalar[3 * nnodes];
  TacsScalar *pert = new TacsScalar[3 * nnodes];  // perturbation for nodes
  TacsScalar *Xp = new TacsScalar[3 * nnodes];
  TacsScalar *mat = new TacsScalar[nvars * nvars];

  // Generate a perturbation vector
  TacsGenerateRandomArray(pert, 3 * nnodes);
  TacsGenerateRandomArray(psi, nvars);
  TacsGenerateRandomArray(phi, nvars);
  memset(result, 0, 3 * nnodes * sizeof(TacsScalar));

  double scale = 1.0 * rand() / RAND_MAX;

  // Compute the element matrix
  element->addMatXptSensInnerProduct(elemType, elemIndex, time, scale, psi, phi,
                                     Xpts, vars, result);

  // Perturb the nodes in the forward sense
  TacsForwardDiffPerturb(Xp, 3 * nnodes, Xpts, pert, dh);
  element->getMatType(elemType, elemIndex, time, Xp, vars, mat);
  TacsScalar forward = 0.0;
  for (int i = 0; i < nvars; i++) {
    for (int j = 0; j < nvars; j++) {
      forward += mat[nvars * i + j] * psi[i] * phi[j];
    }
  }

  // Perturb the nodes in the reverse sense
  TacsBackwardDiffPerturb(Xp, 3 * nnodes, Xpts, pert, dh);
  element->getMatType(elemType, elemIndex, time, Xp, vars, mat);
  TacsScalar backward = 0.0;
  for (int i = 0; i < nvars; i++) {
    for (int j = 0; j < nvars; j++) {
      backward += mat[nvars * i + j] * psi[i] * phi[j];
    }
  }

  // Form the FD/CS approximate
  TacsFormDiffApproximate(&forward, &backward, 1, dh);

  // Compute the fd result
  TacsScalar res = 0.0;
  for (int i = 0; i < 3 * nnodes; i++) {
    res += result[i] * pert[i];
  }
  TacsScalar fd = scale * forward;

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(&res, &fd, 1, &max_err_index);
  double max_rel = TacsGetMaxRelError(&res, &fd, 1, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the element matrix type for element %s.\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    fprintf(stderr,
            "The product of a random vector and the stiffness matrix\n");
    TacsPrintErrorComponents(stderr, "Element Xpt product", &res, &fd, 1);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] psi;
  delete[] phi;
  delete[] pert;
  delete[] Xp;
  delete[] mat;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestElementMatSVSens(TACSElement *element, ElementMatrixType elemType,
                             int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             double dh, int test_print_level,
                             double test_fail_atol, double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  TacsScalar *result = new TacsScalar[nvars];
  TacsScalar *psi = new TacsScalar[nvars];
  TacsScalar *phi = new TacsScalar[nvars];

  TacsScalar *pert = new TacsScalar[nvars];  // perturbation for vars
  TacsScalar *q = new TacsScalar[nvars];
  TacsScalar *mat = new TacsScalar[nvars * nvars];

  // Generate a perturbation vector
  TacsGenerateRandomArray(pert, nvars);
  TacsGenerateRandomArray(psi, nvars);
  TacsGenerateRandomArray(phi, nvars);

  // Compute the element matrix
  element->getMatSVSensInnerProduct(elemType, elemIndex, time, psi, phi, Xpts,
                                    vars, result);

  // Perturb the variables in the forward sense
  TacsForwardDiffPerturb(q, nvars, vars, pert, dh);
  element->getMatType(elemType, elemIndex, time, Xpts, q, mat);
  TacsScalar forward = 0.0;
  for (int i = 0; i < nvars; i++) {
    for (int j = 0; j < nvars; j++) {
      forward += mat[nvars * i + j] * psi[i] * phi[j];
    }
  }

  // Perturb the variables in the reverse sense
  TacsBackwardDiffPerturb(q, nvars, vars, pert, dh);
  element->getMatType(elemType, elemIndex, time, Xpts, q, mat);
  TacsScalar backward = 0.0;
  for (int i = 0; i < nvars; i++) {
    for (int j = 0; j < nvars; j++) {
      backward += mat[nvars * i + j] * psi[i] * phi[j];
    }
  }

  // Form the FD/CS approximate
  TacsFormDiffApproximate(&forward, &backward, 1, dh);

  // Compute the fd result
  TacsScalar res = 0.0;
  for (int i = 0; i < nvars; i++) {
    res += result[i] * pert[i];
  }
  TacsScalar fd = forward;

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(&res, &fd, 1, &max_err_index);
  double max_rel = TacsGetMaxRelError(&res, &fd, 1, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the element matrix type for element %s.\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    fprintf(stderr,
            "The product of a random vector and the stiffness matrix\n");
    TacsPrintErrorComponents(stderr, "Element SV product", &res, &fd, 1);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] psi;
  delete[] phi;
  delete[] pert;
  delete[] q;
  delete[] mat;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

int TacsTestElementBasisFunctions(TACSElementBasis *basis, double dh,
                                  int test_print_level, double test_fail_atol,
                                  double test_fail_rtol) {
  // Set the failure parameter
  int fail = 0;

  // Get the number of parameters/number of nodes
  int nparams = basis->getNumParameters();
  int nnodes = basis->getNumNodes();

  // Create an array to store the values of the adjoint-residual
  // product
  double *result = new double[nparams * nnodes];
  memset(result, 0, nparams * nnodes * sizeof(double));

  // Generate a random array of values
  double *fd = new double[nparams * nnodes];

  double pt0[3];
  TacsGenerateRandomArray(pt0, 3);

  double *N0 = new double[nnodes];
  double *N = new double[nnodes];
  basis->computeBasisGradient(pt0, N0, result);

  // Compute the finite-difference
  for (int i = 0; i < nparams; i++) {
    double pt[3];
    memcpy(pt, pt0, 3 * sizeof(double));
    pt[i] = pt[i] + dh;
    basis->computeBasis(pt, N);

    for (int j = 0; j < nnodes; j++) {
      fd[nparams * j + i] = (N[j] - N0[j]) / dh;
    }
  }

#ifndef TACS_USE_COMPLEX
  // Compute the error
  int max_err_index, max_rel_index;
  double max_err =
      TacsGetMaxError(result, fd, nparams * nnodes, &max_err_index);
  double max_rel =
      TacsGetMaxRelError(result, fd, nparams * nnodes, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the basis functions\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dN/dp", result, fd, nparams * nnodes);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
#endif  // not TACS_USE_COMPLEX

  delete[] result;
  delete[] fd;
  delete[] N;
  delete[] N0;

#ifndef TACS_USE_COMPLEX
  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);
#endif  // not TACS_USE_COMPLEX

  return fail;
}

int TacsTestElementBasisFaceNormals(TACSElementBasis *basis, double dh,
                                    int test_print_level, double test_fail_atol,
                                    double test_fail_rtol) {
  int fail = 0;

  // The failure value
  int nfaces = basis->getNumElementFaces();
  int nnodes = basis->getNumNodes();

  // The sensitivity values
  TacsScalar dfdA, dfdX[3], dfdXd[9], dfdn[3];

  // Generate random sensitivity inputs
  TacsGenerateRandomArray(&dfdA, 1);
  TacsGenerateRandomArray(dfdX, 3);
  TacsGenerateRandomArray(dfdXd, 9);
  TacsGenerateRandomArray(dfdn, 3);

  // Other vectors
  TacsScalar *Xpts = new TacsScalar[3 * nnodes];
  TacsScalar *pert = new TacsScalar[3 * nnodes];
  TacsGenerateRandomArray(pert, 3 * nnodes);
  TacsGenerateRandomArray(Xpts, 3 * nnodes);

  TacsScalar *Xpts_pert = new TacsScalar[3 * nnodes];
  TacsScalar *dfdXpts = new TacsScalar[3 * nnodes];

  for (int face = 0; face < nfaces; face++) {
    for (int n = 0; n < basis->getNumFaceQuadraturePoints(face); n++) {
      TacsScalar X[3], Xd[9], normal[3];
      TacsScalar A = basis->getFaceNormal(face, n, Xpts, X, Xd, normal);
      memset(dfdXpts, 0, 3 * nnodes * sizeof(TacsScalar));

      basis->addFaceNormalXptSens(face, n, A, Xd, normal, dfdA, dfdX, dfdXd,
                                  dfdn, dfdXpts);

      TacsScalar proj = 0.0;
      for (int k = 0; k < 3 * nnodes; k++) {
        proj += pert[k] * dfdXpts[k];
      }

      // Perturb the variables in the forward sense
      TacsForwardDiffPerturb(Xpts_pert, 3 * nnodes, Xpts, pert, dh);

      TacsScalar forward_X[3], forward_Xd[9], forward_normal[3];
      memset(forward_X, 0, 3 * sizeof(TacsScalar));
      memset(forward_Xd, 0, 9 * sizeof(TacsScalar));
      memset(forward_normal, 0, 3 * sizeof(TacsScalar));
      TacsScalar forward_A = basis->getFaceNormal(face, n, Xpts_pert, forward_X,
                                                  forward_Xd, forward_normal);

      TacsBackwardDiffPerturb(Xpts_pert, 3 * nnodes, Xpts, pert, dh);

      TacsScalar backward_X[3], backward_Xd[9], backward_normal[3];
      memset(backward_X, 0, 3 * sizeof(TacsScalar));
      memset(backward_Xd, 0, 9 * sizeof(TacsScalar));
      memset(backward_normal, 0, 3 * sizeof(TacsScalar));
      TacsScalar backward_A = basis->getFaceNormal(
          face, n, Xpts_pert, backward_X, backward_Xd, backward_normal);

      // Form the FD/CS approximate
      TacsFormDiffApproximate(&forward_A, &backward_A, 1, dh);
      TacsFormDiffApproximate(forward_X, backward_X, 3, dh);
      TacsFormDiffApproximate(forward_Xd, backward_Xd, 9, dh);
      TacsFormDiffApproximate(forward_normal, backward_normal, 3, dh);

      TacsScalar fd = 0.0;
      fd = dfdA * forward_A;
      for (int i = 0; i < 3; i++) {
        fd += forward_X[i] * dfdX[i] + forward_normal[i] * dfdn[i];
      }
      for (int i = 0; i < 9; i++) {
        fd += forward_Xd[i] * dfdXd[i];
      }

      // Compute the error
      int max_err_index, max_rel_index;
      double max_err = TacsGetMaxError(&proj, &fd, 1, &max_err_index);
      double max_rel = TacsGetMaxRelError(&proj, &fd, 1, &max_rel_index);
      int flag = (max_err > test_fail_atol || max_rel > test_fail_rtol);
      fail = flag || fail;

      if (test_print_level > 0) {
        fprintf(stderr,
                "Testing the face normal derivative on face %d quadrature "
                "point %d\n",
                face, n);
        fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
                max_err_index);
        fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
                max_rel_index);
      }

      // Print the error if required
      if (test_print_level > 1) {
        TacsPrintErrorComponents(stderr, "df/dXpts", &proj, &fd, 1);
      }
      if (test_print_level) {
        fprintf(stderr, "\n");
      }
    }
  }

  delete[] Xpts;
  delete[] pert;
  delete[] Xpts_pert;
  delete[] dfdXpts;

  return fail;
}

int TacsTestElementBasisJacobianTransform(TACSElementBasis *basis, double dh,
                                          int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol) {
  int fail = 0;

  // The failure value
  int nquad = basis->getNumQuadraturePoints();
  int nnodes = basis->getNumNodes();

  // The sensitivity values
  TacsScalar dfddetJ, dfdXd[9], dfdJ[9];

  // Generate random sensitivity inputs
  TacsGenerateRandomArray(&dfddetJ, 1);
  TacsGenerateRandomArray(dfdXd, 9);
  TacsGenerateRandomArray(dfdJ, 9);

  // Other vectors
  TacsScalar *Xpts = new TacsScalar[3 * nnodes];
  TacsScalar *pert = new TacsScalar[3 * nnodes];
  TacsGenerateRandomArray(pert, 3 * nnodes);
  TacsGenerateRandomArray(Xpts, 3 * nnodes);

  TacsScalar *Xpts_pert = new TacsScalar[3 * nnodes];
  TacsScalar *dfdXpts = new TacsScalar[3 * nnodes];

  for (int n = 0; n < nquad; n++) {
    double pt[3];
    basis->getQuadraturePoint(n, pt);

    TacsScalar Xd[9], J[9];
    basis->getJacobianTransform(n, pt, Xpts, Xd, J);
    memset(dfdXpts, 0, 3 * nnodes * sizeof(TacsScalar));

    basis->addJacobianTransformXptSens(n, pt, Xd, J, dfddetJ, dfdXd, dfdJ,
                                       dfdXpts);

    TacsScalar proj = 0.0;
    for (int k = 0; k < 3 * nnodes; k++) {
      proj += pert[k] * dfdXpts[k];
    }

    // Perturb the variables in the forward sense
    TacsForwardDiffPerturb(Xpts_pert, 3 * nnodes, Xpts, pert, dh);

    TacsScalar forward_Xd[9], forward_J[9];
    memset(forward_Xd, 0, 9 * sizeof(TacsScalar));
    memset(forward_J, 0, 9 * sizeof(TacsScalar));
    TacsScalar forward_detJ =
        basis->getJacobianTransform(n, pt, Xpts_pert, forward_Xd, forward_J);

    TacsBackwardDiffPerturb(Xpts_pert, 3 * nnodes, Xpts, pert, dh);

    TacsScalar backward_Xd[9], backward_J[9];
    memset(backward_Xd, 0, 9 * sizeof(TacsScalar));
    memset(backward_J, 0, 9 * sizeof(TacsScalar));
    TacsScalar backward_detJ =
        basis->getJacobianTransform(n, pt, Xpts_pert, backward_Xd, backward_J);

    // Form the FD/CS approximate
    TacsFormDiffApproximate(&forward_detJ, &backward_detJ, 1, dh);
    TacsFormDiffApproximate(forward_Xd, backward_Xd, 9, dh);
    TacsFormDiffApproximate(forward_J, backward_J, 9, dh);

    TacsScalar fd = dfddetJ * forward_detJ;
    for (int i = 0; i < 9; i++) {
      fd += forward_Xd[i] * dfdXd[i] + forward_J[i] * dfdJ[i];
    }

    // Compute the error
    int max_err_index, max_rel_index;
    double max_err = TacsGetMaxError(&proj, &fd, 1, &max_err_index);
    double max_rel = TacsGetMaxRelError(&proj, &fd, 1, &max_rel_index);
    int flag = (max_err > test_fail_atol || max_rel > test_fail_rtol);
    fail = flag || fail;

    if (test_print_level > 0) {
      fprintf(stderr,
              "Testing the Jacobian transformation derivative at quadrature "
              "point %d\n",
              n);
      fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
              max_err_index);
      fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
              max_rel_index);
    }

    // Print the error if required
    if (test_print_level > 1) {
      TacsPrintErrorComponents(stderr, "df/dXpts", &proj, &fd, 1);
    }
    if (test_print_level) {
      fprintf(stderr, "\n");
    }
  }

  delete[] Xpts;
  delete[] pert;
  delete[] Xpts_pert;
  delete[] dfdXpts;

  return fail;
}

/*
  Test the quantity output design variable sensitivities
*/
int TacsTestElementQuantityDVSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh, int test_print_level,
    double test_fail_atol, double test_fail_rtol) {
  // Retrieve the number of variables
  int dvs_per_node = element->getDesignVarsPerNode();
  int ndvs = element->getDesignVarNums(elemIndex, 0, NULL);
  int num_dvs = dvs_per_node * ndvs;

  TacsScalar *x = new TacsScalar[num_dvs];
  element->getDesignVars(elemIndex, num_dvs, x);

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[num_dvs];
  TacsScalar *fd = new TacsScalar[num_dvs];
  memset(result, 0, num_dvs * sizeof(TacsScalar));
  memset(fd, 0, num_dvs * sizeof(TacsScalar));

  // Generate a random array
  TacsScalar f0[9], dfdq[9];
  TacsGenerateRandomArray(dfdq, 9);

  // Evaluate the pointwise quantity of interest
  double pt[3] = {-0.125, 0.383, -0.233};
  TacsScalar detXd;
  int count = element->evalPointQuantity(elemIndex, quantityType, time, -1, pt,
                                         Xpts, vars, dvars, ddvars, &detXd, f0);

  if (count > 0) {
    // Evaluate the derivative
    double scale = 1.0 * rand() / RAND_MAX;
    element->addPointQuantityDVSens(elemIndex, quantityType, time, scale, -1,
                                    pt, Xpts, vars, dvars, ddvars, dfdq,
                                    num_dvs, result);

    for (int i = 0; i < num_dvs; i++) {
      TacsScalar xi = x[i];

#ifdef TACS_USE_COMPLEX
      x[i] = xi + TacsScalar(0.0, dh);
#else
      x[i] = xi + dh;
#endif  // TACS_USE_COMPLEX

      TacsScalar f1[9];
      element->setDesignVars(elemIndex, num_dvs, x);
      element->evalPointQuantity(elemIndex, quantityType, time, -1, pt, Xpts,
                                 vars, dvars, ddvars, &detXd, f1);

      for (int j = 0; j < count; j++) {
#ifdef TACS_USE_COMPLEX
        fd[i] += scale * dfdq[j] * TacsImagPart(f1[j]) / dh;
#else
        fd[i] += scale * dfdq[j] * (f1[j] - f0[j]) / dh;
#endif  // TACS_USE_COMPLEX

        x[i] = xi;
      }
    }
  }

  element->setDesignVars(elemIndex, num_dvs, x);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, num_dvs, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, num_dvs, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the pointwise quantity for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdx", result, fd, num_dvs);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] fd;
  delete[] x;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the quantity output state variable sensitivities
*/
int TacsTestElementQuantitySVSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh, int test_print_level,
    double test_fail_atol, double test_fail_rtol) {
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[nvars];
  TacsScalar *fd = new TacsScalar[nvars];
  TacsScalar *q = new TacsScalar[nvars];
  TacsScalar *qdot = new TacsScalar[nvars];
  TacsScalar *qddot = new TacsScalar[nvars];
  memset(result, 0, nvars * sizeof(TacsScalar));

  // Generate a random array
  TacsScalar f0[9], dfdq[9];
  TacsGenerateRandomArray(dfdq, 9);

  // Evaluate the pointwise quantity of interest
  double pt[3] = {-0.125, 0.383, -0.233};
  TacsScalar detXd;
  int count = element->evalPointQuantity(elemIndex, quantityType, time, -1, pt,
                                         Xpts, vars, dvars, ddvars, &detXd, f0);

  if (count > 0) {
    // Evaluate the derivative
    double scale = 1.0 * rand() / RAND_MAX;
    TacsScalar alpha = 1.0 * rand() / RAND_MAX;
    TacsScalar beta = 1.0 * rand() / RAND_MAX;
    TacsScalar gamma = 1.0 * rand() / RAND_MAX;
    element->addPointQuantitySVSens(
        elemIndex, quantityType, time, scale * alpha, scale * beta,
        scale * gamma, -1, pt, Xpts, vars, dvars, ddvars, dfdq, result);

    for (int i = 0; i < nvars; i++) {
      memcpy(q, vars, nvars * sizeof(TacsScalar));
      memcpy(qdot, dvars, nvars * sizeof(TacsScalar));
      memcpy(qddot, ddvars, nvars * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
      q[i] = vars[i] + alpha * TacsScalar(0.0, dh);
      qdot[i] = dvars[i] + beta * TacsScalar(0.0, dh);
      qddot[i] = ddvars[i] + gamma * TacsScalar(0.0, dh);
#else
      q[i] = vars[i] + alpha * dh;
      qdot[i] = dvars[i] + beta * dh;
      qddot[i] = ddvars[i] + gamma * dh;
#endif  // TACS_USE_COMPLEX

      TacsScalar f1[9];
      element->evalPointQuantity(elemIndex, quantityType, time, -1, pt, Xpts, q,
                                 qdot, qddot, &detXd, f1);

      for (int j = 0; j < count; j++) {
#ifdef TACS_USE_COMPLEX
        fd[i] += scale * dfdq[j] * TacsImagPart(f1[j]) / dh;
#else
        fd[i] += scale * dfdq[j] * (f1[j] - f0[j]) / dh;
#endif  // TACS_USE_COMPLEX
      }
    }
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, nvars, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the pointwise quantity for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdu", result, fd, nvars);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] fd;
  delete[] q;
  delete[] qdot;
  delete[] qddot;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the quantity output element node sensitivities
*/
int TacsTestElementQuantityXptSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh, int test_print_level,
    double test_fail_atol, double test_fail_rtol) {
  int nnodes = element->getNumNodes();

  // Create an array to store the result
  TacsScalar *result = new TacsScalar[3 * nnodes];
  TacsScalar *fd = new TacsScalar[3 * nnodes];
  memset(result, 0, 3 * nnodes * sizeof(TacsScalar));
  memset(fd, 0, 3 * nnodes * sizeof(TacsScalar));

  TacsScalar *Xt = new TacsScalar[3 * nnodes];

  // Generate a random array
  TacsScalar f0[9], dfdq[9];
  TacsGenerateRandomArray(dfdq, 9);

  TacsScalar dfdetXd = 0.0;
  TacsGenerateRandomArray(&dfdetXd, 1);

  // Evaluate the pointwise quantity of interest
  double pt[3] = {-0.125, 0.383, -0.233};
  TacsScalar detXd0;
  int count =
      element->evalPointQuantity(elemIndex, quantityType, time, -1, pt, Xpts,
                                 vars, dvars, ddvars, &detXd0, f0);

  if (count > 0) {
    // Evaluate the derivative
    double scale = 1.0 * rand() / RAND_MAX;
    element->addPointQuantityXptSens(elemIndex, quantityType, time, scale, -1,
                                     pt, Xpts, vars, dvars, ddvars, dfdetXd,
                                     dfdq, result);

    for (int i = 0; i < 3 * nnodes; i++) {
      memcpy(Xt, Xpts, 3 * nnodes * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
      Xt[i] = Xpts[i] + TacsScalar(0.0, dh);
#else
      Xt[i] = Xpts[i] + dh;
#endif  // TACS_USE_COMPLEX

      TacsScalar detXd1, f1[9];
      element->evalPointQuantity(elemIndex, quantityType, time, -1, pt, Xt,
                                 vars, dvars, ddvars, &detXd1, f1);

#ifdef TACS_USE_COMPLEX
      fd[i] += scale * dfdetXd * TacsImagPart(detXd1) / dh;
#else
      fd[i] += scale * dfdetXd * (detXd1 - detXd0) / dh;
#endif  // TACS_USE_COMPLEX

      for (int j = 0; j < count; j++) {
#ifdef TACS_USE_COMPLEX
        fd[i] += scale * dfdq[j] * TacsImagPart(f1[j]) / dh;
#else
        fd[i] += scale * dfdq[j] * (f1[j] - f0[j]) / dh;
#endif  // TACS_USE_COMPLEX
      }
    }
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, 3 * nnodes, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, 3 * nnodes, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative of the pointwise quantity for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdXpts", result, fd, 3 * nnodes);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  delete[] result;
  delete[] fd;
  delete[] Xt;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestElementBasis(TACSElementBasis *basis, double dh,
                         int test_print_level, double test_fail_atol,
                         double test_fail_rtol) {
  int fail = 0;
  int flag = TacsTestElementBasisFunctions(basis, dh, test_print_level,
                                           test_fail_atol, test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestElementBasisFaceNormals(basis, dh, test_print_level,
                                         test_fail_atol, test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestElementBasisJacobianTransform(basis, dh, test_print_level,
                                               test_fail_atol, test_fail_rtol);
  fail = flag || fail;

  return fail;
}

// Set the maximum number of variables at a given node
static const int MAX_VARS_PER_NODE = 10;

int TacsTestElementModelJacobian(TACSElementModel *model, int elemIndex,
                                 const double time, double dh,
                                 int test_print_level, double test_fail_atol,
                                 double test_fail_rtol) {
  int fail = 0;

  const int vars_per_node = model->getVarsPerNode();
  const int num_params = model->getNumParameters();

  // Set the quadrature point
  const int n = 0;
  const double pt[3] = {-0.125, 0.383, -0.233};

  int Jac_nnz;
  const int *Jac_pairs;
  model->getWeakMatrixNonzeros(TACS_JACOBIAN_MATRIX, elemIndex, &Jac_nnz,
                               &Jac_pairs);

  // Fill in this data with random numbers...
  TacsScalar X[3], Xd[9];
  TacsGenerateRandomArray(X, 3);
  TacsGenerateRandomArray(Xd, 9);

  TacsScalar Ut[3 * MAX_VARS_PER_NODE], Ux[3 * MAX_VARS_PER_NODE];
  TacsGenerateRandomArray(Ut, 3 * vars_per_node);
  TacsGenerateRandomArray(Ux, num_params * vars_per_node);

  // Check consistency with the integrad
  TacsScalar pDUt[3 * MAX_VARS_PER_NODE], pDUx[3 * MAX_VARS_PER_NODE];
  model->evalWeakIntegrand(elemIndex, time, n, pt, X, Xd, Ut, Ux, pDUt, pDUx);

  // Evaluate the weak form of the model
  TacsScalar DUt[3 * MAX_VARS_PER_NODE], DUx[3 * MAX_VARS_PER_NODE];
  TacsScalar Jac[36 * MAX_VARS_PER_NODE * MAX_VARS_PER_NODE];
  model->evalWeakMatrix(TACS_JACOBIAN_MATRIX, elemIndex, time, n, pt, X, Xd, Ut,
                        Ux, DUt, DUx, Jac);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err =
      TacsGetMaxError(DUt, pDUt, 3 * vars_per_node, &max_err_index);
  double max_rel =
      TacsGetMaxRelError(DUt, pDUt, 3 * vars_per_node, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the consistency of the weak form for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "DUt", DUt, pDUt, 3 * vars_per_node);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  max_err =
      TacsGetMaxError(DUx, pDUx, num_params * vars_per_node, &max_err_index);
  max_rel =
      TacsGetMaxRelError(DUx, pDUx, num_params * vars_per_node, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the Jacobian of the weak form for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "DUx", DUx, pDUx,
                             num_params * vars_per_node);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Now compute the derivative, component by component. Assume that
  // the Jacobian is supplied in a sparse format (the most common)
  TacsScalar fd[36 * MAX_VARS_PER_NODE * MAX_VARS_PER_NODE];
  for (int k = 0; k < Jac_nnz; k++) {
    TacsScalar pUt[3 * MAX_VARS_PER_NODE], pUx[3 * MAX_VARS_PER_NODE];
    memcpy(pUt, Ut, 3 * vars_per_node * sizeof(TacsScalar));
    memcpy(pUx, Ux, num_params * vars_per_node * sizeof(TacsScalar));

    int ix = Jac_pairs[2 * k];
    int jx = Jac_pairs[2 * k + 1];

    // Check which parameter to perturb
    if (jx % (num_params + 3) < 3) {
      int i = jx / (num_params + 3);
      int j = jx % (num_params + 3);
#ifdef TACS_USE_COMPLEX
      pUt[3 * i + j] = Ut[3 * i + j] + TacsScalar(0.0, dh);
#else
      pUt[3 * i + j] = Ut[3 * i + j] + dh;
#endif  // TACS_USE_COMPLEX
    } else {
      int i = jx / (num_params + 3);
      int j = (jx % (num_params + 3)) - 3;
#ifdef TACS_USE_COMPLEX
      pUx[num_params * i + j] = Ux[num_params * i + j] + TacsScalar(0.0, dh);
#else
      pUx[num_params * i + j] = Ux[num_params * i + j] + dh;
#endif  // TACS_USE_COMPLEX
    }

    model->evalWeakIntegrand(elemIndex, time, n, pt, X, Xd, pUt, pUx, pDUt,
                             pDUx);

    if (ix % (num_params + 3) < 3) {
      int i = ix / (num_params + 3);
      int j = ix % (num_params + 3);
#ifdef TACS_USE_COMPLEX
      fd[k] = TacsImagPart(pDUt[3 * i + j]) / dh;
#else
      fd[k] = (pDUt[3 * i + j] - DUt[3 * i + j]) / dh;
#endif  // TACS_USE_COMPLEX
    } else {
      int i = ix / (num_params + 3);
      int j = (ix % (num_params + 3)) - 3;
#ifdef TACS_USE_COMPLEX
      fd[k] = TacsImagPart(pDUx[num_params * i + j]) / dh;
#else
      fd[k] = (pDUx[num_params * i + j] - DUx[num_params * i + j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  max_err = TacsGetMaxError(Jac, fd, Jac_nnz, &max_err_index);
  max_rel = TacsGetMaxRelError(Jac, fd, Jac_nnz, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the Jacobian of the weak form for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "Jac", Jac, fd, Jac_nnz);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

int TacsTestElementModelAdjXptSensProduct(TACSElementModel *model,
                                          int elemIndex, const double time,
                                          double dh, int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol) {
  int fail = 0;

  const int num_params = model->getNumParameters();
  const int vars_per_node = model->getVarsPerNode();

  // Set the quadrature point
  const int n = 0;
  const double pt[3] = {-0.233, 0.383, -0.125};

  // Fill in this data with random numbers...
  TacsScalar X[3], Xd[9];
  TacsGenerateRandomArray(X, 3);
  TacsGenerateRandomArray(Xd, 9);

  TacsScalar Ut[3 * MAX_VARS_PER_NODE], Ux[3 * MAX_VARS_PER_NODE];
  TacsGenerateRandomArray(Ut, 3 * vars_per_node);
  TacsGenerateRandomArray(Ux, num_params * vars_per_node);

  TacsScalar Psi[MAX_VARS_PER_NODE], Psix[3 * MAX_VARS_PER_NODE];
  TacsGenerateRandomArray(Psi, vars_per_node);
  TacsGenerateRandomArray(Psix, num_params * vars_per_node);

  // Evaluate the weak integrand
  TacsScalar DUt[3 * MAX_VARS_PER_NODE], DUx[3 * MAX_VARS_PER_NODE];
  model->evalWeakIntegrand(elemIndex, time, n, pt, X, Xd, Ut, Ux, DUt, DUx);

  TacsScalar result = 0.0;
  for (int i = 0; i < vars_per_node; i++) {
    result += Psi[i] * (DUt[3 * i] + DUt[3 * i + 1] + DUt[3 * i + 2]);
    for (int j = 0; j < num_params; j++) {
      result += Psix[num_params * i + j] * DUx[num_params * i + j];
    }
  }

  // Evaluate the derivative
  TacsScalar product;
  TacsScalar dfdX[3], dfdXd[9], dfdUx[3 * MAX_VARS_PER_NODE];
  TacsScalar dfdPsix[3 * MAX_VARS_PER_NODE];
  model->evalWeakAdjXptSensProduct(elemIndex, time, n, pt, X, Xd, Ut, Ux, Psi,
                                   Psix, &product, dfdX, dfdXd, dfdUx, dfdPsix);

  int max_err_index, max_rel_index;
  double max_err = 0.0, max_rel = 0.0;
  max_err = TacsGetMaxError(&product, &result, 1, &max_err_index);
  max_rel = TacsGetMaxRelError(&product, &result, 1, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the product term for consistency for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "product", &product, &result, 1);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Perturb X
  TacsScalar fdX[3];
  for (int k = 0; k < 3; k++) {
    TacsScalar pX[3];
    memcpy(pX, X, 3 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    pX[k] = X[k] + TacsScalar(0.0, dh);
#else
    pX[k] = X[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar pDUt[3 * MAX_VARS_PER_NODE], pDUx[3 * MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, n, pt, pX, Xd, Ut, Ux, pDUt,
                             pDUx);

    product = 0.0;
    for (int i = 0; i < vars_per_node; i++) {
      product += Psi[i] * (pDUt[3 * i] + pDUt[3 * i + 1] + pDUt[3 * i + 2]);
      for (int j = 0; j < num_params; j++) {
        product += Psix[num_params * i + j] * pDUx[num_params * i + j];
      }
    }

#ifdef TACS_USE_COMPLEX
    fdX[k] = TacsImagPart(product) / dh;
#else
    fdX[k] = (product - result) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(dfdX, fdX, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dfdX, fdX, 3, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative dfdX for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdX", dfdX, fdX, 3);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Perturb Xd
  TacsScalar fdXd[9];
  for (int k = 0; k < 3 * num_params; k++) {
    TacsScalar pXd[9];
    memcpy(pXd, Xd, 9 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    pXd[k] = Xd[k] + TacsScalar(0.0, dh);
#else
    pXd[k] = Xd[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar pDUt[3 * MAX_VARS_PER_NODE], pDUx[3 * MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, n, pt, X, pXd, Ut, Ux, pDUt,
                             pDUx);

    product = 0.0;
    for (int i = 0; i < vars_per_node; i++) {
      product += Psi[i] * (pDUt[3 * i] + pDUt[3 * i + 1] + pDUt[3 * i + 2]);
      for (int j = 0; j < num_params; j++) {
        product += Psix[num_params * i + j] * pDUx[num_params * i + j];
      }
    }

#ifdef TACS_USE_COMPLEX
    fdXd[k] = TacsImagPart(product) / dh;
#else
    fdXd[k] = (product - result) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(dfdXd, fdXd, 3 * num_params, &max_err_index);
  max_rel = TacsGetMaxRelError(dfdXd, fdXd, 3 * num_params, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative dfdXd for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdXd", dfdXd, fdXd, 3 * num_params);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Perturb Ux
  TacsScalar fdUx[3 * MAX_VARS_PER_NODE];
  for (int k = 0; k < num_params * vars_per_node; k++) {
    TacsScalar pUx[3 * MAX_VARS_PER_NODE];
    memcpy(pUx, Ux, num_params * vars_per_node * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    pUx[k] = Ux[k] + TacsScalar(0.0, dh);
#else
    pUx[k] = Ux[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar pDUt[3 * MAX_VARS_PER_NODE], pDUx[3 * MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, n, pt, X, Xd, Ut, pUx, pDUt,
                             pDUx);

    product = 0.0;
    for (int i = 0; i < vars_per_node; i++) {
      product += Psi[i] * (pDUt[3 * i] + pDUt[3 * i + 1] + pDUt[3 * i + 2]);
      for (int j = 0; j < num_params; j++) {
        product += Psix[num_params * i + j] * pDUx[num_params * i + j];
      }
    }

#ifdef TACS_USE_COMPLEX
    fdUx[k] = TacsImagPart(product) / dh;
#else
    fdUx[k] = (product - result) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err =
      TacsGetMaxError(dfdUx, fdUx, vars_per_node * num_params, &max_err_index);
  max_rel = TacsGetMaxRelError(dfdUx, fdUx, vars_per_node * num_params,
                               &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative dfdUx for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdUx", dfdUx, fdUx,
                             vars_per_node * num_params);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Perturb Psix
  TacsScalar fdPsix[3 * MAX_VARS_PER_NODE];
  for (int k = 0; k < num_params * vars_per_node; k++) {
    TacsScalar pPsix[3 * MAX_VARS_PER_NODE];
    memcpy(pPsix, Psix, num_params * vars_per_node * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    pPsix[k] = Psix[k] + TacsScalar(0.0, dh);
#else
    pPsix[k] = Psix[k] + dh;
#endif  // TACS_USE_COMPLEX

    TacsScalar pDUt[3 * MAX_VARS_PER_NODE], pDUx[3 * MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, n, pt, X, Xd, Ut, Ux, pDUt, pDUx);

    product = 0.0;
    for (int i = 0; i < vars_per_node; i++) {
      product += Psi[i] * (pDUt[3 * i] + pDUt[3 * i + 1] + pDUt[3 * i + 2]);
      for (int j = 0; j < num_params; j++) {
        product += pPsix[num_params * i + j] * pDUx[num_params * i + j];
      }
    }

#ifdef TACS_USE_COMPLEX
    fdPsix[k] = TacsImagPart(product) / dh;
#else
    fdPsix[k] = (product - result) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(dfdPsix, fdPsix, vars_per_node * num_params,
                            &max_err_index);
  max_rel = TacsGetMaxRelError(dfdPsix, fdPsix, vars_per_node * num_params,
                               &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative dfdPsix for model %s\n",
            model->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dfdPsix", dfdPsix, fdPsix,
                             vars_per_node * num_params);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestElementModel(TACSElementModel *model, int elemIndex,
                         const double time, double dh, int test_print_level,
                         double test_fail_atol, double test_fail_rtol) {
  int fail = 0;
  int flag =
      TacsTestElementModelJacobian(model, elemIndex, time, dh, test_print_level,
                                   test_fail_atol, test_fail_rtol);
  fail = flag || fail;

  flag = TacsTestElementModelAdjXptSensProduct(model, elemIndex, time, dh,
                                               test_print_level, test_fail_atol,
                                               test_fail_rtol);
  fail = flag || fail;

  return fail;
}
