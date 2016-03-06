#include <stdlib.h>
#include "TACSElement.h"
#include "tacslapack.h"

/*
  The TACSElement and Element Traction class definitions.

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  Assign variables randomly to an array. This is useful for
  testing various things.
*/
static void generate_random_array( TacsScalar * array, int size, 
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
static double get_max_error( TacsScalar * a, TacsScalar * b, int size, 
                             int * max_index ){
  double max_error = 0.0;
  *max_index = -1;
  
  for ( int i = 0; i < size; i++ ){
    double er = fabs(RealPart(a[i] - b[i]));
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
static double get_max_rel_error( TacsScalar * a, TacsScalar * b, int size, 
                                 int * max_index ){
  double max_error = 0.0;
  *max_index = -1;
  
  for ( int i = 0; i < size; i++ ){
    double er = 0.0;
    if (a[i] != 0.0){
      er = fabs(RealPart((a[i] - b[i])/a[i]));
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
static void print_error_components( FILE * fp, const char * descript,
                                    TacsScalar * a, TacsScalar * b, 
                                    int size ){
  fprintf(fp, "%*s[   ] %15s %15s %15s\n",
          (int)strlen(descript), "Val", "Analytic", "Approximate", "Rel. Error");
  for ( int i = 0; i < size; i++ ){
    if (a[i] != 0.0){
      fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n", 
              descript, i, RealPart(a[i]), RealPart(b[i]), 
              fabs(RealPart((a[i] - b[i])/a[i])));
    }
    else {
      fprintf(fp, "%s[%3d] %15.6e %15.6e\n", 
              descript, i, RealPart(a[i]), RealPart(b[i]));
    }
  }  
}

/*
  Perturb the input variables in the forward sense
*/
static void forward_perturb( TacsScalar * out, int size,
                             const TacsScalar * orig,
                             const TacsScalar * pert,
                             TacsScalar dh ){
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
static void backward_perturb( TacsScalar * out, int size,
                              const TacsScalar * orig,
                              const TacsScalar * pert,
                              TacsScalar dh ){
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
static void form_approximate( TacsScalar * forward, 
                              const TacsScalar * backward, 
                              int size,
                              TacsScalar dh ){
#ifdef TACS_USE_COMPLEX
  for ( int k = 0; k < size; k++ ){
    forward[k] = ImagPart(forward[k])/dh;
  }
#else
  for ( int k = 0; k < size; k++ ){
    forward[k] = (forward[k] - backward[k])/(2.0*dh);
  }
#endif // TACS_USE_COMPLEX
}

/*
  The following class performs element-level verification of the
  underlying routines to test that they are correct.  This code uses
  finite-difference when compiled in real mode and complex-step when
  compiled in complex mode.

  input:
  element:  the element object to be tested
  Xpts:     the element nodal locations
*/
TestElement::TestElement( TACSElement * _element, 
                          const TacsScalar _Xpts[] ){
  element = _element;
  element->incref();
  
  // Copy over the values of the variables
  // and the nodes.
  int nvars = element->numVariables();
  int nnodes = 3*element->numNodes();

  vars = new TacsScalar[ nvars ];
  dvars = new TacsScalar[ nvars ];
  ddvars = new TacsScalar[ nvars ];
  Xpts = new TacsScalar[ nnodes ];
  memcpy(Xpts, _Xpts, nnodes*sizeof(TacsScalar));

  generate_random_array(vars, nvars);
  generate_random_array(dvars, nvars);
  generate_random_array(ddvars, nvars);
  
  // Default parameter values
  dh = 1e-7;
  print_level = 0;

  // Set the default values for the failure tolerances. 
  // The most important is the relative tolerance as the absolute
  // error can be large due to large magnitudes of the quantities
  // themselves. (This is especially true for randomly-generated
  // values of the state variables.)
  fail_rtol = 1e-5;
  fail_atol = 10.0;
}

TestElement::~TestElement(){
  element->decref();
  delete [] vars;
  delete [] dvars;
  delete [] ddvars;
  delete [] Xpts;
}

/*
  Test the stiffness matrix using the residual function
*/
int TestElement::testStiffnessMat( int col ){
  int nvars = element->numVariables();
  
  TacsScalar * result = new TacsScalar[nvars];
  TacsScalar * temp = new TacsScalar[nvars];
  TacsScalar * vars_pert = new TacsScalar[nvars]; // perturbation for vars
  TacsScalar * vars_copy = new TacsScalar[nvars];
  TacsScalar * res = new TacsScalar[nvars];
  TacsScalar * mat = new TacsScalar[nvars*nvars];

  if (col >= 0 && col < nvars){
    memset(vars_pert, 0, nvars*sizeof(TacsScalar));
    vars_pert[col] = 1.0;
  }
  else {
    generate_random_array(vars_pert, nvars);
  }

  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  element->getJacobian(mat, alpha, beta, gamma,
		       Xpts, vars, dvars, ddvars);

  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("T", &nvars, &nvars, &a, mat, &nvars,
           vars_pert, &one, &b, result, &one);

  forward_perturb(vars_copy, nvars, vars, vars_pert, dh);
  element->getResidual(res, Xpts, vars_copy, dvars, ddvars);

  backward_perturb(vars_copy, nvars, vars, vars_pert, dh);
  element->getResidual(res, Xpts, vars_copy, dvars, ddvars);

  form_approximate(res, temp, nvars, dh);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(result, res, nvars, &max_err_index);
  double max_rel = get_max_rel_error(result, res, nvars, &max_rel_index);

  if (print_level > 0){
    if (col >= 0 && col < nvars){
      fprintf(stderr, 
              "Testing column %d of the stiffness matrix for element %s.\n",
              col, element->elementName());
    }
    else {
      fprintf(stderr, 
              "Testing the stiffness matrix for element %s.\n",
              element->elementName());
    }
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (print_level > 1){
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
  if (print_level){ fprintf(stderr, "\n"); }

  delete [] result;
  delete [] temp;
  delete [] vars_pert;
  delete [] vars_copy;
  delete [] res;
  delete [] mat;

  return (max_err > fail_atol || max_rel > fail_rtol);
}

/*
  Test the sensitivity of the matrix w.r.t. the design variables
*/
/*
int TestElement::testMatDVSens( ElementMatrixType type ){
  // Get the design variables associated with this element
  int ndvs = element->getNumDesignVars();
  int * dv_nums = new int[ ndvs ];
  int dv_index = 0;
  if (!element->getDesignVarNums(dv_nums, &dv_index, ndvs)){
    fprintf(stderr, 
            "The number of design variables defined by %s is inconsistent\n",
            element->elementName());
    return 0;
  }
  
  int max_dv = 0;
  for ( int k = 0; k < ndvs; k++ ){
    if (dv_nums[k]+1 > max_dv){
      max_dv = dv_nums[k]+1;
    }
  }

  TacsScalar * dvs = new TacsScalar[ max_dv ];
  element->getDesignVars(dvs, max_dv);

  if (print_level){
    fprintf(stderr, 
	    "Testing the matrix sensivity w.r.t. the design variables \
for element %s.\n",
	    element->elementName());
  }

  int fail_flag = 0;

  int nvars = element->numVariables();
  TacsScalar * res = new TacsScalar[nvars];
  TacsScalar * matSens = new TacsScalar[nvars*nvars];
  TacsScalar * matSensApprox = new TacsScalar[nvars*nvars];
  TacsScalar * temp = new TacsScalar[nvars*nvars];

  const char * descript = "StiffMatDVSens";
  if (type == MASS_MATRIX){
    descript = "MassMatDVSens";
  }
  else if (type == GEOMETRIC_STIFFNESS_MATRIX){
    descript = "GeoStiffMatDVSens";
  }

  for ( int k = 0; k < dv_index; k++ ){
    // Compute the derivative of the residual here
    element->getMatTypeDVSens(dv_nums[k], type, 1.0, matSens, vars, 
			      Xpts, NORMAL);

    // Test the residual here
    TacsScalar x = dvs[dv_nums[k]];
    dvs[dv_nums[k]] = x + dh;
    element->setDesignVars(dvs, max_dv);
    if (type == STIFFNESS_MATRIX){
      element->getMat(matSensApprox, res, vars, Xpts, NORMAL);
    }
    else {
      element->getMatType(type, 1.0, matSensApprox, vars, Xpts, NORMAL);
    }

    dvs[dv_nums[k]] = x - dh;
    element->setDesignVars(dvs, max_dv);
    if (type == STIFFNESS_MATRIX){
      element->getMat(temp, res, vars, Xpts, NORMAL);
    }
    else {
      element->getMatType(type, 1.0, temp, vars,  Xpts, NORMAL);
    }
    
    dvs[dv_nums[k]] = x;
    element->setDesignVars(dvs, max_dv);

    form_approximate(matSensApprox, temp, nvars*nvars, dh);

    // Compute the error
    int max_err_index, max_rel_index;
    double max_err = get_max_error(matSens, matSensApprox, 
                                   nvars*nvars, &max_err_index);
    double max_rel = get_max_rel_error(matSens, matSensApprox, 
                                       nvars*nvars, &max_rel_index);
   
    if (print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (print_level > 1){
      fprintf(stderr, 
              "The sensitivity of the matrix w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, descript,
                             matSens, matSensApprox, nvars*nvars);
    }
    if (print_level){ fprintf(stderr, "\n"); }

    fail_flag = (fail_flag || (max_err > fail_atol || max_rel > fail_rtol));
  }

  delete [] dv_nums;
  delete [] dvs;
  delete [] matSens;
  delete [] matSensApprox;
  delete [] temp;
  delete [] res;

  return fail_flag;
}
*/
/*
  Test the derivative of the strain with respect to the state
  variables

  addPtwiseStrainSVSens adds the derivative of a function of the
  strain with respect to the element state variables.
*/
int TestElement::testStrainSVSens( const double pt[] ){
  int nvars = element->numVariables();
  int nstress = element->numStresses();

  TacsScalar * strain = new TacsScalar[nstress];
  TacsScalar * elementSens = new TacsScalar[nvars];
  TacsScalar * elementSensApprox = new TacsScalar[nvars];
  TacsScalar * temp = new TacsScalar[nvars];

  TacsScalar * vars_copy = new TacsScalar[nvars];

  // Randomly assign the strainSens
  TacsScalar * strainSens = new TacsScalar[nstress];
  generate_random_array(strainSens, nstress);

  TacsScalar scale = 1.0;
  memset(elementSens, 0, nvars*sizeof(TacsScalar));
  element->addStrainSVSens(elementSens, pt, scale,
			   strainSens, Xpts, vars);

  memset(elementSensApprox, 0, nvars*sizeof(TacsScalar));
  memset(temp, 0, nvars*sizeof(TacsScalar));

  for ( int k = 0; k < nvars; k++ ){
#ifdef TACS_USE_COMPLEX
    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = TacsScalar(vars[k], dh);
    element->getStrain(strain, pt, Xpts, vars_copy);

    for ( int j = 0; j < nstress; j++ ){
      elementSensApprox[k] += strainSens[j]*imag(strain[j])/dh;
    }
#else 
    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = vars[k] + dh;
    element->getStrain(strain, pt, Xpts, vars_copy);
        
    for ( int j = 0; j < nstress; j++ ){
      elementSensApprox[k] += strainSens[j]*strain[j];
    }

    memcpy(vars_copy, vars, nvars*sizeof(TacsScalar));
    vars_copy[k] = vars[k] - dh;
    element->getStrain(strain, pt, Xpts, vars_copy);
    
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

  if (print_level > 0){
    fprintf(stderr, 
	    "Testing the strain sensivity w.r.t. state variables for %s.\n",
	    element->elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the strain w.r.t. the state variables ons\n");
    print_error_components(stderr, "strainSVSens", 
                           elementSens, elementSensApprox, nvars);
  }
  if (print_level){ fprintf(stderr, "\n"); }

  return (max_err > fail_atol || max_rel > fail_rtol);
}

/*
  Test the derivative of the strain w.r.t. the nodal coordinates
*/
int TestElement::testJacobianXptSens( const double pt[] ){
  // First, test the derivative w.r.t. the nodal coordinates
  int nnodes = 3*element->numNodes(); // actually 3 times the number of nodes
  TacsScalar jacSens, jacSensApprox, temp;
  TacsScalar * jacXptSens = new TacsScalar[nnodes];
  TacsScalar * Xpt_pert = new TacsScalar[nnodes];
  TacsScalar * Xpt_copy = new TacsScalar[nnodes];
  
  generate_random_array(Xpt_pert, nnodes);

  // Compute the sensitivity
  element->getJacobianXptSens(jacXptSens, pt, Xpts);
  
  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("N", &one, &nnodes, &a, jacXptSens, &one,
           Xpt_pert, &one, &b, &jacSens, &one);
  
  forward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  jacSensApprox = element->getJacobian(pt, Xpt_copy);

  backward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  temp = element->getJacobian(pt, Xpt_copy);

  form_approximate(&jacSensApprox, &temp, 1, dh);

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(&jacSens, &jacSensApprox, 1,
                                 &max_err_index);
  double max_rel = get_max_rel_error(&jacSens, &jacSensApprox, 1,
                                     &max_rel_index);

  if (print_level > 0){
    fprintf(stderr, 
	    "Testing the det. Jacobian sensivity w.r.t. the nodes for %s.\n",
	    element->elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the det. Jacobian w.r.t. the nodes\n");
    print_error_components(stderr, "jacXptSens", 
                           &jacSens, &jacSensApprox, 1);
  }
  if (print_level){ fprintf(stderr, "\n"); }

  delete [] jacXptSens;
  delete [] Xpt_pert;
  delete [] Xpt_copy;

  return (max_err > fail_atol || max_rel > fail_rtol);
}

/*
  Test the derivative of the strain w.r.t. the nodal coordinates
*/
/*
int TestElement::testStrainXptSens( const double pt[] ){
  // First, test the derivative w.r.t. the nodal coordinates
  int nstress = element->numStresses();
  int nnodes = 3*element->numNodes(); // actually 3 times the number of nodes
  TacsScalar * temp = new TacsScalar[nstress];
  TacsScalar * strainSens = new TacsScalar[nstress];
  TacsScalar * strainSensApprox = new TacsScalar[nnodes];
  TacsScalar * strainXptSens = new TacsScalar[nnodes];
  TacsScalar * Xpt_pert = new TacsScalar[nnodes];
  TacsScalar * Xpt_copy = new TacsScalar[nnodes];
  
  generate_random_array(Xpt_pert, nnodes);

  // Compute the sensitivity
  element->addStrainXptSens(temp, strainXptSens, pt, Xpts, vars);
  element->getStrain(strainSens, pt, Xpts, vars);

  if (print_level > 0){
    int max_err_index, max_rel_index;
    double max_err = get_max_error(strainSens, temp,
                                   nstress, &max_err_index);
    double max_rel = get_max_rel_error(strainSens, temp,
                                       nstress, &max_rel_index);

    fprintf(stderr, 
	    "Testing consistency of the strain calculations for element %s\n",
            element->elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
    fprintf(stderr, "\n");
  }

  int one = 1;
  TacsScalar a = 1.0, b = 0.0;
  BLASgemv("N", &nstress, &nnodes, &a, strainXptSens, &nstress,
           Xpt_pert, &one, &b, strainSens, &one);
  
  forward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  element->getStrain(strainSensApprox, pt, Xpt_copy, vars);
    
  backward_perturb(Xpt_copy, nnodes, Xpts, Xpt_pert, dh);
  element->getStrain(temp, pt, Xpt_copy, vars);
    
  form_approximate(strainSensApprox, temp, nstress, dh);
    
  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(strainSens, strainSensApprox, 
                                 nstress, &max_err_index);
  double max_rel = get_max_rel_error(strainSens, strainSensApprox, 
                                     nstress, &max_rel_index);
    
  if (print_level > 0){
    fprintf(stderr, 
            "Testing the strain sensivity w.r.t. the nodes for element %s.\n",
            element->elementName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }  
  // Print the error if required
  if (print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the strain w.r.t. the nodal locations\n");
    print_error_components(stderr, "strainXptSens", 
                           strainSens, strainSensApprox, nstress);
  }
  if (print_level){ fprintf(stderr, "\n"); }
  
  delete [] strainXptSens;
  delete [] Xpt_pert;
  delete [] strainSens;
  delete [] strainSensApprox;
  delete [] Xpt_copy;
  delete [] temp;

  return (max_err > fail_atol || max_rel > fail_rtol);
}
*/
/*
  Test the implementation of a constitutive class
*/
TestConstitutive::TestConstitutive( TACSConstitutive * _con ){
  con = _con;
  con->incref();

  // Default parameter values
  dh = 0.5e-7;
  print_level = 0;

  // Set the default values for the failure tolerances. 
  fail_rtol = 1e-6;
  fail_atol = 10.0;
}

TestConstitutive::~TestConstitutive(){
  con->decref();
}

/*
  This generates strains from a given stress.

  This results in more realistic stress levels and avoids issues
  realted to very large stress resulting in unrealistic failure
  prediction.  
*/
void TestConstitutive::compute_strain( TacsScalar strain[], 
                                       const double pt[],
                                       const TacsScalar stress[] ){
  int nstress = con->getNumStresses();

  TacsScalar * C = new TacsScalar[ nstress*nstress ];
  memset(C, 0, nstress*nstress*sizeof(TacsScalar));

  for ( int k = 0; k < nstress; k++ ){
    memset(strain, 0, nstress*sizeof(TacsScalar));

    strain[k] = 1.0;
    con->calculateStress(pt, strain, &C[nstress*k]);
  }

  memcpy(strain, stress, nstress*sizeof(TacsScalar));

  // Solve for the strain given the stresses
  int * ipiv = new int[ nstress ];
  int one = 1;
  int info = 0;
  LAPACKgesv(&nstress, &one, C, &nstress, ipiv, 
             strain, &nstress, &info);
  if (info != 0){
    fprintf(stderr, "Problem with the constitutive matrix!\n");
  }

  delete [] C;
  delete [] ipiv;
}

/*
  Test the sensitivity of the strain to a failure 
*/
int TestConstitutive::testFailStrainSens( const double pt[] ){
  int nstress = con->getNumStresses();

  TacsScalar * strain = new TacsScalar[ nstress ];
  TacsScalar * sens = new TacsScalar[ nstress ];
  TacsScalar * strainCopy = new TacsScalar[ nstress ];
  TacsScalar * sensApprox = new TacsScalar[ nstress ];

  // First use sens as a temporary randomly generated array
  // These respresent randomly generated stresses on the interval [-1, 1]
  generate_random_array(sens, nstress, -1.0, 1.0);

  // Now geneate the corresponding strains
  compute_strain(strain, pt, sens);
  con->failureStrainSens(pt, strain, sens);

  // Determine approximations of sens
  for ( int k = 0; k < nstress; k++ ){
    memcpy(strainCopy, strain, nstress*sizeof(TacsScalar));
#ifdef TACS_USE_COMPLEX
    strainCopy[k] += TacsScalar(0.0, dh);
    TacsScalar forward;
    con->failure(pt, strainCopy, &forward);    
    sensApprox[k] = ImagPart(forward)/dh;
#else  
    TacsScalar e = strainCopy[k];
    strainCopy[k] = e + dh;
    TacsScalar forward;
    con->failure(pt, strainCopy, &forward);    

    strainCopy[k] = e - dh;
    TacsScalar backward;
    con->failure(pt, strainCopy, &backward);
    sensApprox[k] = (forward - backward)/(2.0*dh);
#endif
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(sens, sensApprox, 
				 nstress, &max_err_index);
  double max_rel = get_max_rel_error(sens, sensApprox, 
				     nstress, &max_rel_index);

  int fail_flag = (max_err > fail_atol || max_rel > fail_rtol);

  if (print_level > 0){
    fprintf(stderr, 
	    "Testing the failure sensivity for %s.\n",
	    con->constitutiveName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the failure criteria w.r.t. \
the components of strain\n");
    print_error_components(stderr, "failStrainSens", 
                           sens, sensApprox, nstress);
  }
  if (print_level){ fprintf(stderr, "\n"); }

  delete [] sens;
  delete [] strain;
  delete [] strainCopy;
  delete [] sensApprox;

  return fail_flag;
}

/*
int TestConstitutive::testFailDVSens( const double pt[], int ndvs ){ 
  int nstress = con->getNumStresses();

  TacsScalar * strain = new TacsScalar[ nstress ];
  TacsScalar * stress = new TacsScalar[ nstress ];

  // Randomly generate the stresses for failure criteria calculations
  generate_random_array(stress, nstress);

  // Determine the corresponding strains
  compute_strain(strain, pt, stress);
  delete [] stress;

  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int * dv_nums = new int[ ndvs ];
  int dv_index = 0;
  if (!con->getDesignVarNums(dv_nums, &dv_index, ndvs)){
    fprintf(stderr, 
            "The number of design variables defined by %s is inconsistent\n",
            con->constitutiveName());
    return 0;
  }
  
  int max_dv = 0;
  for ( int k = 0; k < ndvs; k++ ){
    if (dv_nums[k]+1 > max_dv){
      max_dv = dv_nums[k]+1;
    }
  }

  TacsScalar * dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (print_level){
    fprintf(stderr, 
	    "Testing the sensivity of the failure criteria w.r.t. the \
design variables for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int fail_flag = 0;

  for ( int k = 0; k < dv_index; k++ ){
    // Compute the derivative of the residual here
    TacsScalar failSens;
    con->failureDVSens(dv_nums[k], pt, strain, &failSens);

    TacsScalar x = dvs[dv_nums[k]];
#ifdef TACS_USE_COMPLEX
    dvs[dv_nums[k]] = x + TacsScalar(0.0, dh);
    con->setDesignVars(dvs, max_dv);
    TacsScalar forward;
    con->failure(pt, strain, &forward);

    TacsScalar failSensApprox = ImagPart(forward)/dh;
#else
    // Compute the finite-difference derivative
    dvs[dv_nums[k]] = x + dh;
    con->setDesignVars(dvs, max_dv);
    TacsScalar forward;
    con->failure(pt, strain, &forward);

    dvs[dv_nums[k]] = x - dh;
    con->setDesignVars(dvs, max_dv);
    TacsScalar backward;
    con->failure(pt, strain, &backward);

    TacsScalar failSensApprox = (forward - backward)/(2.0*dh);
#endif // TACS_USE_COMPLEX

    dvs[dv_nums[k]] = x;
    con->setDesignVars(dvs, max_dv);

    // Compute the error
    int max_err_index, max_rel_index;
    double max_err = get_max_error(&failSens, &failSensApprox, 
                                   1, &max_err_index);
    double max_rel = get_max_rel_error(&failSens, &failSensApprox,
                                       1, &max_rel_index);
   
    if (print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (print_level > 1){
      fprintf(stderr, 
              "The sensitivity failure criteria w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "failDVSens", 
                             &failSens, &failSensApprox, 1);
    }
    if (print_level){ fprintf(stderr, "\n"); }

    fail_flag = (fail_flag || (max_err > fail_atol || max_rel > fail_rtol));
  }

  delete [] strain;
  delete [] dvs;
  delete [] dv_nums;

  return fail_flag;
}
*/
/*
int TestConstitutive::testMassDVSens( const double pt[] ){ 
  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int ndvs = con->getNumDesignVars();
  int * dv_nums = new int[ ndvs ];
  int dv_index = 0;
  if (!con->getDesignVarNums(dv_nums, &dv_index, ndvs)){
    fprintf(stderr, 
            "The number of design variables defined by %s is inconsistent\n",
            con->constitutiveName());
    return 0;
  }
  
  int max_dv = 0;
  for ( int k = 0; k < ndvs; k++ ){
    if (dv_nums[k]+1 > max_dv){
      max_dv = dv_nums[k]+1;
    }
  }

  TacsScalar * dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (print_level){
    fprintf(stderr, 
	    "Testing the sensivity w.r.t. the design variables \
for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int fail_flag = 0;

  for ( int k = 0; k < dv_index; k++ ){
    // Compute the derivative of the residual here
    TacsScalar massSens[6];
    con->pointwiseMassDVSens(dv_nums[k], pt, massSens);

    // Test the residual here
    TacsScalar x = dvs[dv_nums[k]];
#ifdef TACS_USE_COMPLEX
    dvs[dv_nums[k]] = x + TacsScalar(0.0, dh);
    con->setDesignVars(dvs, max_dv);
    TacsScalar forward[6];
    con->pointwiseMass(pt, forward);    
    TacsScalar massSensApprox = ImagPart(forward[0])/dh;
#else
    dvs[dv_nums[k]] = x + dh;
    TacsScalar forward[6];
    con->setDesignVars(dvs, max_dv);
    con->pointwiseMass(pt, forward);

    dvs[dv_nums[k]] = x - dh;
    TacsScalar backward[6];
    con->setDesignVars(dvs, max_dv);
    con->pointwiseMass(pt, backward);

    TacsScalar massSensApprox = (forward[0] - backward[0])/(2.0*dh);
#endif // TACS_USE_COMPLEX
    dvs[dv_nums[k]] = x;
    con->setDesignVars(dvs, max_dv);

    // Compute the error
    int max_err_index, max_rel_index;
    double max_err = get_max_error(massSens, &massSensApprox, 
                                   1, &max_err_index);
    double max_rel = get_max_rel_error(massSens, &massSensApprox,
                                       1, &max_rel_index);
   
    if (print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (print_level > 1){
      fprintf(stderr, 
              "The sensitivity mass w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "massDVSens", 
                             massSens, &massSensApprox, 1);
    }
    if (print_level){ fprintf(stderr, "\n"); }

    fail_flag = (fail_flag || (max_err > fail_atol || max_rel > fail_rtol));
  }

  delete [] dvs;
  delete [] dv_nums;

  return fail_flag;
}
*/
/*
  Test the sensitivity of the strain to a buckling 
*/
int TestConstitutive::testBucklingStrainSens(){
  int nstress = con->getNumStresses();

  TacsScalar * strain = new TacsScalar[ nstress ];
  TacsScalar * sens = new TacsScalar[ nstress ];
  TacsScalar * strainCopy = new TacsScalar[ nstress ];
  TacsScalar * sensApprox = new TacsScalar[ nstress ];

  // First use conSens and linSens as temporary, randomly generated arrays
  // These respresent randomly generated stresses on the interval [-1, 1]
  // memset(conSens, 0, nstress*sizeof(TacsScalar));
  generate_random_array(sens, nstress, -1.0, 1.0);

  // Now geneate the corresponding strains
  double pt[3] = {0.0, 0.0, 0.0};
  compute_strain(strain, pt, sens);
  con->bucklingStrainSens(strain, sens);

  // Determine approximations of sens
  for ( int k = 0; k < nstress; k++ ){
    memcpy(strainCopy, strain, nstress*sizeof(TacsScalar));
    TacsScalar e = strainCopy[k];
    strainCopy[k] = e + dh;
    TacsScalar forward;
    con->buckling(strainCopy, &forward);

    strainCopy[k] = e - dh;
    TacsScalar backward;
    con->buckling(strainCopy, &backward);
    sensApprox[k] = (forward - backward)/(2.0*dh);
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = get_max_error(sens, sensApprox, 
				 nstress, &max_err_index);
  double max_rel = get_max_rel_error(sens, sensApprox, 
				     nstress, &max_rel_index);

  int fail_flag = (max_err > fail_atol || max_rel > fail_rtol);

  if (print_level > 0){
    fprintf(stderr, 
	    "Testing the buckling sensivity for %s.\n",
	    con->constitutiveName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the buckling criteria w.r.t. \
the components of strain\n");
    print_error_components(stderr, "bucklingStrainSens", 
                           sens, sensApprox, nstress);
  }
  if (print_level){ fprintf(stderr, "\n"); }

  delete [] sens;
  delete [] strain;
  delete [] strainCopy;
  delete [] sensApprox;

  return fail_flag;
}

/*
int TestConstitutive::testBucklingDVSens(){ 
  int nstress = con->getNumStresses();

  TacsScalar * strain = new TacsScalar[ nstress ];
  TacsScalar * stress = new TacsScalar[ nstress ];

  // Randomly generate the stresses for buckling criteria calculations
  generate_random_array(stress, nstress);

  // Determine the corresponding strains
  double pt[3] = {0.0, 0.0, 0.0};
  compute_strain(strain, pt, stress);
  delete [] stress;

  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int ndvs = con->getNumDesignVars();
  int * dv_nums = new int[ ndvs ];
  int dv_index = 0;
  if (!con->getDesignVarNums(dv_nums, &dv_index, ndvs)){
    fprintf(stderr, 
            "The number of design variables defined by %s is inconsistent\n",
            con->constitutiveName());
    return 0;
  }
  
  int max_dv = 0;
  for ( int k = 0; k < ndvs; k++ ){
    if (dv_nums[k]+1 > max_dv){
      max_dv = dv_nums[k]+1;
    }
  }

  TacsScalar * dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (print_level){
    fprintf(stderr, 
	    "Testing the sensivity of the buckling criteria w.r.t. the \
design variables for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int fail_flag = 0;

  for ( int k = 0; k < dv_index; k++ ){
    // Compute the derivative of the residual here
    TacsScalar bucklingSens;
    con->bucklingDVSens(dv_nums[k], strain, &bucklingSens);

    // Test the residual here
    TacsScalar x = dvs[dv_nums[k]];
    dvs[dv_nums[k]] = x + dh;
    con->setDesignVars(dvs, max_dv);
    TacsScalar forward;
    con->buckling(strain, &forward);

    dvs[dv_nums[k]] = x - dh;
    con->setDesignVars(dvs, max_dv);
    TacsScalar backward;
    con->buckling(strain, &backward);
    
    dvs[dv_nums[k]] = x;
    con->setDesignVars(dvs, max_dv);

    TacsScalar bucklingSensApprox = (forward - backward)/(2.0*dh);

    // Compute the error
    int max_err_index, max_rel_index;
    double max_err = get_max_error(&bucklingSens, &bucklingSensApprox, 
                                   1, &max_err_index);
    double max_rel = get_max_rel_error(&bucklingSens, &bucklingSensApprox,
                                       1, &max_rel_index);
   
    if (print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (print_level > 1){
      fprintf(stderr, 
              "The sensitivity buckling criteria w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "bucklingDVSens", 
                             &bucklingSens, &bucklingSensApprox, 1);
    }
    if (print_level){ fprintf(stderr, "\n"); }

    fail_flag = (fail_flag || (max_err > fail_atol || max_rel > fail_rtol));
  }

  delete [] strain;
  delete [] dvs;
  delete [] dv_nums;

  return fail_flag;
}
*/
