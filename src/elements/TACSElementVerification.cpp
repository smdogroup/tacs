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
#include "tacslapack.h"

/*
  Assign variables randomly to an array. This is useful for
  testing various things.
*/
void TacsGenerateRandomArray( TacsScalar *array, int size,
                              TacsScalar lower,
                              TacsScalar upper ){
  for ( int i = 0; i < size; i++ ){
    array[i] = (upper - lower)*(rand()/((double)RAND_MAX+1)) + lower;
  }
}

/*
  Find the largest absolute value of the difference between the
  arrays a and b
*/
double TacsGetMaxError( TacsScalar *a, TacsScalar *b, int size,
                        int *max_index ){
  double max_error = 0.0;
  *max_index = -1;

  for ( int i = 0; i < size; i++ ){
    double er = fabs(TacsRealPart(a[i] - b[i]));
    if (i == 0 || er > max_error){
      max_error = er;
      *max_index = i;
    }
  }
  return max_error;
}

/*
  Find the maximum relative error between a and b and return the
*/
double TacsGetMaxRelError( TacsScalar *a, TacsScalar *b, int size,
                           int *max_index ){
  double max_error = 0.0;
  *max_index = -1;

  for ( int i = 0; i < size; i++ ){
    double er = 0.0;
    if (a[i] != 0.0){
      er = fabs(TacsRealPart((a[i] - b[i])/a[i]));
    }
    if (i == 0 || er > max_error){
      max_error = er;
      *max_index = i;
    }
  }
  return max_error;
}

/*
  Print out the values and the relative errors
*/
void TacsPrintErrorComponents( FILE *fp, const char *descript,
                               TacsScalar *a, TacsScalar *b,
                               int size ){
  fprintf(fp, "%*s[   ] %15s %15s %15s\n",
          (int)strlen(descript), "Val", "Analytic", "Approximate", "Rel. Error");
  for ( int i = 0; i < size; i++ ){
    if (a[i] != 0.0){
      fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n",
              descript, i, TacsRealPart(a[i]), TacsRealPart(b[i]),
              fabs(TacsRealPart((a[i] - b[i])/a[i])));
    }
    else {
      fprintf(fp, "%s[%3d] %15.6e %15.6e\n",
              descript, i, TacsRealPart(a[i]), TacsRealPart(b[i]));
    }
  }
}

/*
  Perturb the input variables in the forward sense
*/
void TacsForwardDiffPerturb( TacsScalar *out, int size,
                             const TacsScalar *orig,
                             const TacsScalar *pert,
                             double dh ){
#ifdef TACS_USE_COMPLEX
  for ( int k = 0; k < size; k++ ){
    out[k] = orig[k] + TacsScalar(0.0, dh)*pert[k];
  }
#else
  for ( int k = 0; k < size; k++ ){
    out[k] = orig[k] + dh*pert[k];
  }
#endif // TACS_USE_COMPLEX
}

/*
  Perturb the variables in the backward sense
*/
void TacsBackwardDiffPerturb( TacsScalar *out, int size,
                              const TacsScalar *orig,
                              const TacsScalar *pert,
                              double dh ){
#ifdef TACS_USE_COMPLEX
  for ( int k = 0; k < size; k++ ){
    out[k] = orig[k];
  }
#else
  for ( int k = 0; k < size; k++ ){
    out[k] = orig[k] - dh*pert[k];
  }
#endif // TACS_USE_COMPLEX
}

/*
  Form the forward approximation
*/
void TacsFormDiffApproximate( TacsScalar *forward,
                              const TacsScalar *backward,
                              int size,
                              TacsScalar dh ){
#ifdef TACS_USE_COMPLEX
  for ( int k = 0; k < size; k++ ){
    forward[k] = TacsImagPart(forward[k])/dh;
  }
#else
  for ( int k = 0; k < size; k++ ){
    forward[k] = (forward[k] - backward[k])/(2.0*dh);
  }
#endif // TACS_USE_COMPLEX
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
int TacsTestElementResidual( TACSElement *element,
                             int elemIndex,
                             double time,
                             const TacsScalar Xpts[],
                             const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[],
                             double dh,
                             int test_print_level,
                             double test_fail_atol,
                             double test_fail_rtol ){
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  // Allocate temporary arrays for the computation
  TacsScalar *q = new TacsScalar[ nvars ];
  TacsScalar *dq = new TacsScalar[ nvars ];

  // Temporary storage for the residuals
  TacsScalar *res1 = new TacsScalar[ nvars ];
  TacsScalar *res2 = new TacsScalar[ nvars ];
  TacsScalar *fd = new TacsScalar[ nvars ];

  // Compute the values of the variables at (t + dt)
  for ( int i = 0; i < nvars; i++ ){
    q[i] = vars[i] + dh*dvars[i];
    dq[i] = dvars[i] + dh*ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for ( int i = 0; i < nvars; i++ ){
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
#ifdef TACS_USE_COMPLEX
    dq[i] = dqtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq&, T1, &P1);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    dq[i] = dqtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);
    res1[i] = 0.5*((T1 - P1) - (T2 - P2))/dh;
#endif
    dq[i] = dqtmp;
  }

  // Compute the values of the variables at (t - dt)
  for ( int i = 0; i < nvars; i++ ){
    q[i] = vars[i] - dh*dvars[i];
    dq[i] = dvars[i] - dh*ddvars[i];
  }

  // Evaluate the derivative w.r.t. dot{q}
  for ( int i = 0; i < nvars; i++ ){
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar dqtmp = dq[i];
#ifdef TACS_USE_COMPLEX
    dq[i] = dqtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    dq[i] = dqtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    dq[i] = dqtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);

    // Compute and store the approximation
    res2[i] = 0.5*((T1 - P1) - (T2 - P2))/dh;
#endif
    dq[i] = dqtmp;
  }

  // Evaluate the finite-difference for the first term in Largrange's
  // equations of motion
  for ( int i = 0; i < nvars; i++ ){
    fd[i] = 0.5*(res1[i] - res2[i])/dh;
  }

  // Reset the values of q and dq at time t
  for ( int i = 0; i < nvars; i++ ){
    q[i] = vars[i];
    dq[i] = dvars[i];
  }

  // Compute the contribution from dL/dq^{T}
  for ( int i = 0; i < nvars; i++ ){
    // Evaluate the finite-difference for component i
    TacsScalar T1, P1, T2, P2;
    TacsScalar qtmp = q[i];

#ifdef TACS_USE_COMPLEX
    q[i] = qtmp + TacsScalar(0.0, dh);
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    q[i] = qtmp + dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T1, &P1);

    q[i] = qtmp - dh;
    element->computeEnergies(elemIndex, time, Xpts, q, dq, &T2, &P2);

    // Compute and store the approximation
    res1[i] = 0.5*((T1 - P1) - (T2 - P2))/dh;
#endif
    q[i] = qtmp;
  }

  // Add the result to the finite-difference result
  for ( int i = 0; i < nvars; i++ ){
    fd[i] -= res1[i];
  }

  // Evaluate the residual using the code
  memset(res1, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res1);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(res1, fd, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(res1, fd, nvars, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr,
            "Testing the residual implementation for element %s.\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr,
            "The difference between the FD and true residual is:\n");
    TacsPrintErrorComponents(stderr, "Res error",
                             res1, fd, nvars);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] q;
  delete [] dq;
  delete [] res1;
  delete [] res2;
  delete [] fd;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  The following function tests the consistency between the
  implementation of the residuals and the implementation of the system
  Jacobian.

  input:
  col:   test only the specified column of the matrix
*/
int TacsTestElementJacobian( TACSElement *element,
                             int elemIndex,
                             double time,
                             const TacsScalar Xpts[],
                             const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[],
                             int col,
                             double dh,
                             int test_print_level,
                             double test_fail_atol,
                             double test_fail_rtol ){
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  TacsScalar *result = new TacsScalar[nvars];
  TacsScalar *temp = new TacsScalar[nvars];
  TacsScalar *pert = new TacsScalar[nvars]; // perturbation for vars

  TacsScalar *q = new TacsScalar[nvars];
  TacsScalar *dq = new TacsScalar[nvars];
  TacsScalar *ddq = new TacsScalar[nvars];
  TacsScalar *res = new TacsScalar[nvars];
  TacsScalar *mat = new TacsScalar[nvars*nvars];

  if (col >= 0 && col < nvars){
    memset(pert, 0, nvars*sizeof(TacsScalar));
    pert[col] = 1.0;
  }
  else {
    TacsGenerateRandomArray(pert, nvars);
  }

  // Compute the Jacobian
  double alpha = (1.0*rand())/RAND_MAX;
  double beta = (1.0*rand())/RAND_MAX;
  double gamma = (1.0*rand())/RAND_MAX;

  memset(mat, 0, nvars*nvars*sizeof(TacsScalar));
  element->addJacobian(elemIndex, time, alpha, beta, gamma,
                       Xpts, vars, dvars, ddvars, res, mat);

  // Evaluate the Jacobian
  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("T", &nvars, &nvars, &a, mat, &nvars,
           pert, &one, &b, result, &one);

  // Perturb the variables in the forward sense
  TacsForwardDiffPerturb(q, nvars, vars, pert, alpha*dh);
  TacsForwardDiffPerturb(dq, nvars, dvars, pert, beta*dh);
  TacsForwardDiffPerturb(ddq, nvars, ddvars, pert, gamma*dh);
  memset(res, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, q, dq, ddq, res);

  // Perturb the variables in the backward sens
  TacsBackwardDiffPerturb(q, nvars, vars, pert, alpha*dh);
  TacsBackwardDiffPerturb(dq, nvars, dvars, pert, beta*dh);
  TacsBackwardDiffPerturb(ddq, nvars, ddvars, pert, gamma*dh);
  memset(temp, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, q, dq, ddq, temp);

  // Form the FD/CS approximate
  TacsFormDiffApproximate(res, temp, nvars, dh);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, res, nvars, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, res, nvars, &max_rel_index);

  if (test_print_level > 0){
    if (col >= 0 && col < nvars){
      fprintf(stderr,
              "Testing column %d of the stiffness matrix for element %s.\n",
              col, element->getObjectName());
    }
    else {
      fprintf(stderr,
              "Testing the stiffness matrix for element %s.\n",
              element->getObjectName());
    }
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    if (col >= 0 && col < nvars){
      fprintf(stderr,
              "Column %d of the stiffness matrix is\n", col);
    }
    else {
      fprintf(stderr,
              "The product of a random vector and the stiffness matrix is\n");
    }
    TacsPrintErrorComponents(stderr, "K*u",
                             result, res, nvars);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] result;
  delete [] temp;
  delete [] pert;
  delete [] q;
  delete [] dq;
  delete [] ddq;
  delete [] res;
  delete [] mat;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestAdjResProduct( TACSElement *element,
                           int elemIndex,
                           int dvLen,
                           const TacsScalar *x,
                           double time,
                           const TacsScalar Xpts[],
                           const TacsScalar vars[],
                           const TacsScalar dvars[],
                           const TacsScalar ddvars[],
                           double dh,
                           int test_print_level,
                           double test_fail_atol,
                           double test_fail_rtol ){
  // Retrieve the number of variables
  int nvars = element->getNumVariables();

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[ dvLen ];
  memset(result, 0, dvLen*sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[ nvars ];
  TacsGenerateRandomArray(adjoint, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0*rand()/RAND_MAX;

  element->addAdjResProduct(elemIndex, time, scale, adjoint,
                            Xpts, vars, dvars, ddvars, dvLen, result);

  // Compute the product of the result with a perturbation
  // vector that is equal to perturb = sign(result[k])
  TacsScalar dpdx = 0.0;
  for ( int k = 0; k < dvLen; k++ ){
    dpdx += fabs(result[k]);
  }

  // Allocate an array to store the perturbed design variable
  // values
  TacsScalar *xpert = new TacsScalar[ dvLen ];
  TacsScalar fd_dpdx = 0.0;

  // Zero the residual
  TacsScalar *res = new TacsScalar[ nvars ];

#ifdef TACS_USE_COMPLEX
  // Perturb the design variables: xpert = x + dh*sign(result[k])
  for ( int k = 0; k < dvLen; k++ ){
    if (TacsRealPart(result[k]) >= 0.0){
      xpert[k] = x[k] + TacsScalar(0.0, dh);
    }
    else {
      xpert[k] = x[k] - TacsScalar(0.0, dh);
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  TacsScalar p1 = 0.0;
  memset(res, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  for ( int k = 0; k < nvars; k++ ){
    p1 += scale*res[k]*adjoint[k];
  }

  fd_dpdx = TacsImagPart(p1)/dh;
#else
  // Perturb the design variables: xpert = x + dh*sign(result[k])
  for ( int k = 0; k < dvLen; k++ ){
    if (result[k] >= 0.0){
      xpert[k] = x[k] + dh;
    }
    else {
      xpert[k] = x[k] - dh;
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  memset(res, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);

  TacsScalar p1 = 0.0;
  for ( int k = 0; k < nvars; k++ ){
    p1 += scale*res[k]*adjoint[k];
  }

  // Pertub the design variables: xpert = x - dh*sign(result[k])
  for ( int k = 0; k < dvLen; k++ ){
    if (result[k] >= 0.0){
      xpert[k] = x[k] - dh;
    }
    else {
      xpert[k] = x[k] + dh;
    }
  }
  element->setDesignVars(elemIndex, dvLen, xpert);

  // Compute the residual again
  memset(res, 0, nvars*sizeof(TacsScalar));
  element->addResidual(elemIndex, time, Xpts, vars, dvars, ddvars, res);
  TacsScalar p2 = 0.0;
  for ( int k = 0; k < nvars; k++ ){
    p2 += scale*res[k]*adjoint[k];
  }

  // Compute the finite-difference approximation
  fd_dpdx = 0.5*(p1 - p2)/dh;
#endif

  // Set the design variable values
  element->setDesignVars(elemIndex, dvLen, x);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(&dpdx, &fd_dpdx, 1, &max_err_index);
  double max_rel = TacsGetMaxRelError(&dpdx, &fd_dpdx, 1,
                                      &max_rel_index);

  test_print_level = 2;
  if (test_print_level > 0){
    fprintf(stderr,
            "Testing the derivative of the adjoint-residual product for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "Adj-Res product",
                             &dpdx, &fd_dpdx, 1);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] result;
  delete [] adjoint;
  delete [] xpert;
  delete [] res;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestAdjResXptProduct( TACSElement *element,
                              int elemIndex,
                              double time,
                              const TacsScalar Xpts[],
                              const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              double dh,
                              int test_print_level,
                              double test_fail_atol,
                              double test_fail_rtol ){
  int nvars = element->getNumNodes()*element->getVarsPerNode();
  int nnodes = element->getNumNodes();

  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[ 3*nnodes ];
  memset(result, 0, 3*nnodes*sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[ nvars ];
  TacsGenerateRandomArray(adjoint, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0*rand()/RAND_MAX;
  element->addAdjResXptProduct(elemIndex, time, scale, adjoint,
                               Xpts, vars, dvars, ddvars, result);

  // Allocate space to store the results
  TacsScalar *fd = new TacsScalar[ 3*nnodes ];
  TacsScalar *X = new TacsScalar[ 3*nnodes ];
  TacsScalar *res = new TacsScalar[ nvars ];

  for ( int k = 0; k < 3*nnodes; k++ ){
    // Copy the points
    memcpy(X, Xpts, 3*nnodes*sizeof(TacsScalar));

    // Perturb the nodes in the forward sense
    TacsScalar one = 1.0;
    TacsForwardDiffPerturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars*sizeof(TacsScalar));
    element->addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
    TacsScalar p1 = 0.0;
    for ( int i = 0; i < nvars; i++ ){
      p1 += scale*adjoint[i]*res[i];
    }

    // Perturb the nodes in the reverse sense
    TacsBackwardDiffPerturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars*sizeof(TacsScalar));
    element->addResidual(elemIndex, time, X, vars, dvars, ddvars, res);
    TacsScalar p2 = 0.0;
    for ( int i = 0; i < nvars; i++ ){
      p2 += scale*adjoint[i]*res[i];
    }

    // Form the approximation
    TacsFormDiffApproximate(&p1, &p2, 1, dh);

    // Set the
    fd[k] = p1;
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, 3*nnodes, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, 3*nnodes, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr,
            "Testing the derivative of the adjoint-residual product for %s\n",
            element->getObjectName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "Adj-Res Xpt product",
                             result, fd, 3*nnodes);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] result;
  delete [] adjoint;
  delete [] fd;
  delete [] X;
  delete [] res;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TacsTestElementBasis( TACSElementBasis *basis,
                          double dh,
                          int test_print_level,
                          double test_fail_atol,
                          double test_fail_rtol ){
  int nparams = basis->getNumParameters();
  int nnodes = basis->getNumNodes();

  // Create an array to store the values of the adjoint-residual
  // product
  double *result = new double[ nparams*nnodes ];
  memset(result, 0, 3*nnodes*sizeof(double));

  // Generate a random array of values
  double *fd = new double[ nparams*nnodes ];

  double pt0[3];
  TacsGenerateRandomArray(pt0, 3);

  double *N0 = new double[ nnodes ];
  double *N = new double[ nnodes ];
  basis->computeBasis(pt0, N0);

  // Compute the finite-difference
  for ( int i = 0; i < nparams; i++ ){
    double pt[3];
    memcpy(pt, pt0, 3*sizeof(double));
    pt[i] = pt[0] + dh;
    basis->computeBasis(pt, N);

    for ( int j = 0; j < nnodes; j++ ){
      fd[nparams*j + i] = (N[j] - N0[j])/dh;
    }
  }

#ifndef TACS_USE_COMPLEX
  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(result, fd, 3*nnodes, &max_err_index);
  double max_rel = TacsGetMaxRelError(result, fd, 3*nnodes, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr,
            "Testing the derivative of the basis functions\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "dN/dp",
                             result, fd, nparams*nnodes);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }
#endif // TACS_USE_COMPLEX

  delete [] result;
  delete [] fd;
  delete [] N;
  delete [] N0;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}
