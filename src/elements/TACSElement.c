#include <stdlib.h>
#include "TACSElement.h"
#include "tacslapack.h"

/*
  The TACSElement and Element Traction class definitions.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

int TACSElement::test_print_level = 0;
double TACSElement::test_step_size = 1e-6;
double TACSElement::test_fail_rtol = 1e-5;
double TACSElement::test_fail_atol = 1e-30;

/*
  Assign variables randomly to an array. This is useful for
  testing various things.
*/
static void generate_random_array( TacsScalar *array, int size, 
                                   TacsScalar lower=-1.0, 
                                   TacsScalar upper=1.0 ){
  for ( int i = 0; i < size; i++ ){
    array[i] = (upper - lower)*(rand()/((double)RAND_MAX+1)) + lower;
  }
}

/*
  Find the largest absolute value of the difference between the
  arrays a and b
*/
static double get_max_error( TacsScalar *a, TacsScalar *b, int size, 
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
static double get_max_rel_error( TacsScalar *a, TacsScalar *b, int size, 
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
static void print_error_components( FILE *fp, const char *descript,
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
static void forward_perturb( TacsScalar *out, int size,
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
static void backward_perturb( TacsScalar *out, int size,
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
static void form_approximate( TacsScalar *forward, 
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
  Set the static testing information
*/
void TACSElement::setFailTolerances( double fail_rtol, double fail_atol ){
  test_fail_rtol = fail_rtol;
  test_fail_atol = fail_atol;
}
 
void TACSElement::setPrintLevel( int flag ){
  test_print_level = flag;
}

void TACSElement::setStepSize( double dh ){
  test_step_size = dh;
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
int TACSElement::testResidual( double time, 
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[] ){
  // Retrieve the number of variables
  int nvars = numVariables();

  // Allocate temporary arrays for the computation
  TacsScalar *q = new TacsScalar[ nvars ];
  TacsScalar *dq = new TacsScalar[ nvars ];

  // Temporary storage for the residuals
  TacsScalar *res1 = new TacsScalar[ nvars ];
  TacsScalar *res2 = new TacsScalar[ nvars ];
  TacsScalar *fd = new TacsScalar[ nvars ];

  // The step length
  double dh = test_step_size;

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
    computeEnergies(time, &T1, &P1, Xpts, q, dq);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    dq[i] = dqtmp + dh;
    computeEnergies(time, &T1, &P1, Xpts, q, dq);

    dq[i] = dqtmp - dh;
    computeEnergies(time, &T2, &P2, Xpts, q, dq);
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
    computeEnergies(time, &T1, &P1, Xpts, q, dq);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    dq[i] = dqtmp + dh;
    computeEnergies(time, &T1, &P1, Xpts, q, dq);

    dq[i] = dqtmp - dh;
    computeEnergies(time, &T2, &P2, Xpts, q, dq);

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
    computeEnergies(time, &T1, &P1, Xpts, q, dq);
    res1[i] = TacsImagPart((T1 - P1))/dh;
#else
    q[i] = qtmp + dh;
    computeEnergies(time, &T1, &P1, Xpts, q, dq);

    q[i] = qtmp - dh;
    computeEnergies(time, &T2, &P2, Xpts, q, dq);

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
  addResidual(time, res1, Xpts, vars, dvars, ddvars);
  
  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(res1, fd, nvars, &max_err_index);
  double max_rel = get_max_rel_error(res1, fd, nvars, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, 
            "Testing the residual implementation for element %s.\n",
            elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }

  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The difference between the FD and true residual is:\n");
    print_error_components(stderr, "Res error",
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
int TACSElement::testJacobian( double time, 
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[], 
                               int col ){
  int nvars = numVariables();
  
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
    generate_random_array(pert, nvars);
  }
  
  // Compute the Jacobian
  double alpha = (1.0*rand())/RAND_MAX;
  double beta = (1.0*rand())/RAND_MAX;
  double gamma = (1.0*rand())/RAND_MAX;
  
  memset(mat, 0, nvars*nvars*sizeof(TacsScalar));
  addJacobian(time, mat, alpha, beta, gamma,
              Xpts, vars, dvars, ddvars);
  
  // Evaluate the Jacobian
  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("T", &nvars, &nvars, &a, mat, &nvars,
           pert, &one, &b, result, &one);
  
  // The step length
  double dh = test_step_size;

  // Perturb the variables in the forward sense
  forward_perturb(q, nvars, vars, pert, alpha*dh);
  forward_perturb(dq, nvars, dvars, pert, beta*dh);
  forward_perturb(ddq, nvars, ddvars, pert, gamma*dh);
  memset(res, 0, nvars*sizeof(TacsScalar));
  addResidual(time, res, Xpts, q, dq, ddq);

  // Perturb the variables in the backward sens
  backward_perturb(q, nvars, vars, pert, alpha*dh);
  backward_perturb(dq, nvars, dvars, pert, beta*dh);
  backward_perturb(ddq, nvars, ddvars, pert, gamma*dh);
  memset(temp, 0, nvars*sizeof(TacsScalar));
  addResidual(time, temp, Xpts, q, dq, ddq);

  // Form the FD/CS approximate
  form_approximate(res, temp, nvars, dh);
  /* for (int i = 0; i < nvars; i++){ */
  /*   printf("res[%d]: %f\n", i,res[i]); */
  /* } */
  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(result, res, nvars, &max_err_index);
  double max_rel = get_max_rel_error(result, res, nvars, &max_rel_index);

  if (test_print_level > 0){
    if (col >= 0 && col < nvars){
      fprintf(stderr, 
              "Testing column %d of the stiffness matrix for element %s.\n",
              col, elementName());
    }
    else {
      fprintf(stderr, 
              "Testing the stiffness matrix for element %s.\n",
              elementName());
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
    print_error_components(stderr, "K*u",
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
  Test the derivative of the strain with respect to the state
  variables

  addPtwiseStrainSVSens adds the derivative of a function of the
  strain with respect to the element state variables.
*/
int TACSElement::testStrainSVSens( const TacsScalar Xpts[],
                                   const TacsScalar vars[] ){
  // Set the parametric point within the element
  double pt[3];
  pt[0] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[1] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[2] = -1.0 + 2.0*rand()/RAND_MAX;

  // Allocate temporary arrays for the strain values
  int nvars = numVariables();
  int nstress = numStresses();
  TacsScalar *strain = new TacsScalar[nstress];
  TacsScalar *elementSens = new TacsScalar[nvars];
  TacsScalar *elementSensApprox = new TacsScalar[nvars];
  TacsScalar *temp = new TacsScalar[nvars];
  TacsScalar *vars_copy = new TacsScalar[nvars];

  // Randomly assign the strainSens
  TacsScalar *strainSens = new TacsScalar[nstress];
  generate_random_array(strainSens, nstress);

  TacsScalar scale = 1.0;
  memset(elementSens, 0, nvars*sizeof(TacsScalar));
  addStrainSVSens(elementSens, pt, scale,
                  strainSens, Xpts, vars);

  memset(elementSensApprox, 0, nvars*sizeof(TacsScalar));
  memset(temp, 0, nvars*sizeof(TacsScalar));

  // The step length
  double dh = test_step_size;

  for ( int k = 0; k < nvars; k++ ){
#ifdef TACS_USE_COMPLEX
    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = vars[k] + TacsScalar(0.0, dh);
    getStrain(strain, pt, Xpts, vars_copy);

    for ( int j = 0; j < nstress; j++ ){
      elementSensApprox[k] += strainSens[j]*TacsImagPart(strain[j])/dh;
    }
#else 
    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = vars[k] + dh;
    getStrain(strain, pt, Xpts, vars_copy);
        
    for ( int j = 0; j < nstress; j++ ){
      elementSensApprox[k] += strainSens[j]*strain[j];
    }

    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = vars[k] - dh;
    getStrain(strain, pt, Xpts, vars_copy);
    
    TacsScalar temp = 0.0;     
    for ( int j = 0; j < nstress; j++ ){
      temp += strainSens[j]*strain[j];
    }

    elementSensApprox[k] = (elementSensApprox[k] - temp)/(2.0*dh);
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(elementSens, elementSensApprox, nvars,
                                 &max_err_index);
  double max_rel = get_max_rel_error(elementSens, elementSensApprox, nvars,
                                     &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, 
	    "Testing the strain sensivity w.r.t. state variables for %s.\n",
	    elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the strain w.r.t. the state variables ons\n");
    print_error_components(stderr, "strainSVSens", 
                           elementSens, elementSensApprox, nvars);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Free the allocated data
  delete [] strain;
  delete [] strainSens;
  delete [] elementSens;
  delete [] elementSensApprox;
  delete [] temp;
  delete [] vars_copy;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}


/*
  Test the derivative of the strain with respect to the state
  variables

  addPtwiseStrainSVSens adds the derivative of a function of the
  strain with respect to the element state variables.
*/
int TACSElement::testStrainXptSens( const TacsScalar Xpts[],
                                    const TacsScalar vars[] ){
  // Set the parametric point within the element
  double pt[3];
  pt[0] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[1] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[2] = -1.0 + 2.0*rand()/RAND_MAX;

  // Allocate temporary arrays for the strain values
  int nnodes = numNodes();
  int nstress = numStresses();
  TacsScalar *strain = new TacsScalar[ nstress ];
  
  // Set a random derivative of the strain
  TacsScalar *strainSens = new TacsScalar[ nstress ];
  generate_random_array(strainSens, nstress);

  // Get the derivative of the strain w.r.t. the nodes
  TacsScalar *deriv = new TacsScalar[ 3*nnodes ];
  memset(deriv, 0, 3*nnodes*sizeof(TacsScalar));

  TacsScalar scale = 1.0*rand()/RAND_MAX;
  addStrainXptSens(deriv, pt, scale, strainSens, Xpts, vars);

  // Allocate an array to store the derivative
  TacsScalar *fd = new TacsScalar[ 3*nnodes ];
  TacsScalar *X = new TacsScalar[ 3*nnodes ];

  // The step length
  double dh = test_step_size;

  // Evaluate the derivative of the strain w.r.t. the node locations
    for ( int k = 0; k < 3*nnodes; k++ ){
    // Copy the points
    memcpy(X, Xpts, 3*nnodes*sizeof(TacsScalar));

    // Perturb the nodes in the forward sense
    TacsScalar one = 1.0;
    forward_perturb(&X[k], 1, &Xpts[k], &one, dh);
    getStrain(strain, pt, X, vars);
    TacsScalar p1 = 0.0;
    for ( int i = 0; i < nstress; i++ ){
      p1 += scale*strainSens[i]*strain[i];
    }
    
    // Perturb the nodes in the reverse sense
    backward_perturb(&X[k], 1, &Xpts[k], &one, dh);
    getStrain(strain, pt, X, vars);
    TacsScalar p2 = 0.0;
    for ( int i = 0; i < nstress; i++ ){
      p2 += scale*strainSens[i]*strain[i];
    }
    
    // Form the approximation
    form_approximate(&p1, &p2, 1, dh);

    // Set the
    fd[k] = p1;
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(deriv, fd, 3*nnodes,
                                 &max_err_index);
  double max_rel = get_max_rel_error(deriv, fd, 3*nnodes,
                                     &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, 
	    "Testing the strain sensivity w.r.t. node locations for %s.\n",
	    elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the strain w.r.t. the state variables ons\n");
    print_error_components(stderr, "strainXptSens", 
                           deriv, fd, 3*nnodes);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Free the allocated data
  delete [] strain;
  delete [] X;
  delete [] fd;
  delete [] strainSens;
  delete [] deriv;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}

/*
  Test the derivative of the inner product of the adjoint vector and
  the residual with respect to material design variables.
*/
int TACSElement::testAdjResProduct( const TacsScalar *x, int dvLen,
                                    double time, const TacsScalar Xpts[],
                                    const TacsScalar vars[], 
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[] ){

  int nvars = numVariables();
  setDesignVars(x, dvLen);
  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[ dvLen ];
  memset(result, 0, dvLen*sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[ nvars ];
  generate_random_array(adjoint, nvars);
  
  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0;
  
  addAdjResProduct(time, scale,
                   result, dvLen, adjoint,
                   Xpts, vars, dvars, ddvars);
  
  // Compute the product of the result with a perturbation 
  // vector that is equal to perturb = sign(result[k])
  TacsScalar dpdx = 0.0;
  for ( int k = 0; k < dvLen; k++ ){
    dpdx += fabs(result[k]);
  }

  // The step length
  double dh = test_step_size;

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
  setDesignVars(xpert, dvLen);

  TacsScalar p1 = 0.0;
  memset(res, 0, nvars*sizeof(TacsScalar));
  addResidual(time, res, Xpts, vars, dvars, ddvars);
  for ( int k = 0; k < nvars; k++ ){
    p1 += res[k]*adjoint[k];
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
  setDesignVars(xpert, dvLen);

  memset(res, 0, nvars*sizeof(TacsScalar));
  addResidual(time, res, Xpts, vars, dvars, ddvars);

  TacsScalar p1 = 0.0;
  for ( int k = 0; k < nvars; k++ ){
    p1 += res[k]*adjoint[k];
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
  setDesignVars(xpert, dvLen);

  // Compute the residual again
  memset(res, 0, nvars*sizeof(TacsScalar));
  addResidual(time, res, Xpts, vars, dvars, ddvars);
  TacsScalar p2 = 0.0;
  for ( int k = 0; k < nvars; k++ ){
    p2 += res[k]*adjoint[k];
  }

  // Compute the finite-difference approximation
  fd_dpdx = 0.5*(p1 - p2)/dh;
#endif

  // Set the design variable values
  setDesignVars(x, dvLen);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(&dpdx, &fd_dpdx, 1, &max_err_index);
  double max_rel = get_max_rel_error(&dpdx, &fd_dpdx, 1,
                                     &max_rel_index);
 
  test_print_level = 2;
  if (test_print_level > 0){
    fprintf(stderr, 
            "Testing the derivative of the adjoint-residual product for %s\n",
            elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    print_error_components(stderr, "Adj-Res product",
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
int TACSElement::testAdjResXptProduct( double time, 
                                       const TacsScalar Xpts[],
                                       const TacsScalar vars[], 
                                       const TacsScalar dvars[],
                                       const TacsScalar ddvars[] ){
  int nnodes = numNodes();
  int nvars = numVariables();
 
  // Create an array to store the values of the adjoint-residual
  // product
  TacsScalar *result = new TacsScalar[ 3*nnodes ];
  memset(result, 0, 3*nnodes*sizeof(TacsScalar));

  // Generate a random array of values
  TacsScalar *adjoint = new TacsScalar[ nvars ];
  generate_random_array(adjoint, nvars);

  // Evaluate the derivative of the adjoint-residual product
  double scale = 1.0*rand()/RAND_MAX;
  addAdjResXptProduct(time, scale, result, adjoint,
                      Xpts, vars, dvars, ddvars);
  
  // The step length
  double dh = test_step_size;

  // Allocate space to store the results
  TacsScalar *fd = new TacsScalar[ 3*nnodes ];
  TacsScalar *X = new TacsScalar[ 3*nnodes ];
  TacsScalar *res = new TacsScalar[ nvars ];

  for ( int k = 0; k < 3*nnodes; k++ ){
    // Copy the points
    memcpy(X, Xpts, 3*nnodes*sizeof(TacsScalar));

    // Perturb the nodes in the forward sense
    TacsScalar one = 1.0;
    forward_perturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars*sizeof(TacsScalar));
    addResidual(time, res, X, vars, dvars, ddvars);
    TacsScalar p1 = 0.0;
    for ( int i = 0; i < nvars; i++ ){
      p1 += scale*adjoint[i]*res[i];
    }
    
    // Perturb the nodes in the reverse sense
    backward_perturb(&X[k], 1, &Xpts[k], &one, dh);
    memset(res, 0, nvars*sizeof(TacsScalar));
    addResidual(time, res, X, vars, dvars, ddvars);
    TacsScalar p2 = 0.0;
    for ( int i = 0; i < nvars; i++ ){
      p2 += scale*adjoint[i]*res[i];
    }
    
    // Form the approximation
    form_approximate(&p1, &p2, 1, dh);

    // Set the
    fd[k] = p1;
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(result, fd, 3*nnodes, &max_err_index);
  double max_rel = get_max_rel_error(result, fd, 3*nnodes, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, 
            "Testing the derivative of the adjoint-residual product for %s\n",
            elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    print_error_components(stderr, "Adj-Res Xpt product",
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
  Test the derivative of the strain w.r.t. the nodal coordinates
*/
int TACSElement::testJacobianXptSens( const TacsScalar Xpts[] ){
  // Set the parametric point within the element
  double pt[3];
  pt[0] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[1] = -1.0 + 2.0*rand()/RAND_MAX;
  pt[2] = -1.0 + 2.0*rand()/RAND_MAX;

  // First, test the derivative w.r.t. the nodal coordinates
  int nnodes = 3*numNodes(); // actually 3 times the number of nodes
  TacsScalar jacSens, jacSensApprox, temp;
  TacsScalar *jacXptSens = new TacsScalar[nnodes];
  TacsScalar *Xpt_pert = new TacsScalar[nnodes];
  TacsScalar *Xpt_copy = new TacsScalar[nnodes];
  
  generate_random_array(Xpt_pert, nnodes);

  // Compute the sensitivity
  getDetJacobianXptSens(jacXptSens, pt, Xpts);
  
  // The step length
  double dh = test_step_size;

  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("N", &one, &nnodes, &a, jacXptSens, &one,
           Xpt_pert, &one, &b, &jacSens, &one);
  
  forward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  jacSensApprox = getDetJacobian(pt, Xpt_copy);

  backward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  temp = getDetJacobian(pt, Xpt_copy);

  form_approximate(&jacSensApprox, &temp, 1, dh);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(&jacSens, &jacSensApprox, 1,
                                 &max_err_index);
  double max_rel = get_max_rel_error(&jacSens, &jacSensApprox, 1,
                                     &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, 
	    "Testing the det. Jacobian sensivity w.r.t. the nodes for %s.\n",
	    elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the det. Jacobian w.r.t. the nodes\n");
    print_error_components(stderr, "jacXptSens", 
                           &jacSens, &jacSensApprox, 1);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] jacXptSens;
  delete [] Xpt_pert;
  delete [] Xpt_copy;

  return (max_err > test_fail_atol || max_rel > test_fail_rtol);
}


