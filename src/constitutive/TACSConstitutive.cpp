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

#include <stdlib.h>
#include <math.h>
#include "tacslapack.h"
#include "TACSConstitutive.h"
#include "FElibrary.h"

/*
  It is not necessary to override these function in order to implement the class
*/
const char * TACSConstitutive::TACSObjectName(){ 
  return this->constitutiveName(); 
}

/*
  Write a two-dimensional representation of the failure envelope to a
  file

  input:
  file_name:  the file to write
  npts:       the number of points to include in the envelope plot
  gpt:        the Gauss point to evaluate the failure criteria
  x_stress:   the x components of the stress
  y_stress:   the y components of the stress
*/
void TACSConstitutive::writeFailureEnvelope( const char * file_name, int npts,
					     const double gpt[], 
					     const TacsScalar x_stress[], 
					     const TacsScalar y_stress[] ){
  FILE * fp = fopen(file_name, "w");

  // The failure criteria solution tolerance
  double tol = 1e-8;

  // The max number of Newton iterations to use
  int max_newton_iters = 50;

  if (fp){
    int nstress = getNumStresses();

    // Compute an explicit form of the constitutive relationship
    int * ipiv = new int[ nstress ];
    TacsScalar * e = new TacsScalar[ nstress ];
    TacsScalar * f_sens = new TacsScalar[ nstress ];
    TacsScalar * x_strain = new TacsScalar[ nstress ];
    TacsScalar * y_strain = new TacsScalar[ nstress ];
    TacsScalar * C = new TacsScalar[ nstress*nstress ];

    for ( int i = 0; i < nstress; i++ ){
      memset(e, 0, nstress*sizeof(TacsScalar));
      e[i] = 1.0;
      calculateStress(gpt, e, &C[i*nstress]);
    }

    // Factor the constitutive matrix
    int info;
    LAPACKgetrf(&nstress, &nstress, C, &nstress, ipiv, &info);
    
    // Copy over the stress values and solve for the strain
    memcpy(x_strain, x_stress, nstress*sizeof(TacsScalar));
    memcpy(y_strain, y_stress, nstress*sizeof(TacsScalar));

    int one = 1;
    LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, 
		x_strain, &nstress, &info);
    LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, 
		y_strain, &nstress, &info);

    fprintf(fp, "Variables = ");
    for ( int j = 0; j < nstress; j++ ){
      fprintf(fp, "s%d ", j);
    }
    fprintf(fp, "\n");

    TacsScalar P = 1.0;
    for ( int k = 0; k < npts; k++ ){
      TacsScalar theta = (2.0*M_PI*k)/(npts-1);
      TacsScalar c = cos(theta), s = sin(theta);

      // Solve the equation: 
      // failure(gpt, P*(c*x_strain + s*y_strain)) - 1.0 = 0           
      for ( int i = 0; i < max_newton_iters; i++ ){
	for ( int j = 0; j < nstress; j++ ){
	  e[j] = P*(c*x_strain[j] + s*y_strain[j]);
	}

	// Compute the failure criterion and the derivative of the
	// failure criterion w.r.t. the load parameter P
	TacsScalar fail, failSens;
	failure(gpt, e, &fail);
	
	if (fabs(TacsRealPart(fail) - 1.0) < tol){
	  break;
	}

	failureStrainSens(gpt, e, f_sens);

	failSens = 0.0;
	for ( int j = 0; j < nstress; j++ ){
	  failSens += f_sens[j]*(c*x_strain[j] + s*y_strain[j]);
	}

	// Compute the Newton update to the boundary
	P = P - (fail - 1.0)/failSens;
      }

      // Compute the strain at the final point
      for ( int j = 0; j < nstress; j++ ){
	e[j] = P*(c*x_strain[j] + s*y_strain[j]);
      }

      // Compute the corresponding stress at the final point
      calculateStress(gpt, e, f_sens);

      // Print out the result to the file
      for ( int j = 0; j < nstress; j++ ){
	fprintf(fp, "%15.8f ", TacsRealPart(f_sens[j]));
      }
      fprintf(fp, "\n");
    }

    fclose(fp);

    delete [] ipiv;
    delete [] e;
    delete [] f_sens;
    delete [] x_strain;
    delete [] y_strain;
    delete [] C;
  }
}

/*
  Write a two-dimensional representation of the buckling envelope to a
  file

  input:
  file_name:  the file to write
  npts:       the number of points to include in the envelope plot
  gpt:        the Gauss point to evaluate the failure criteria
  x_stress:   the x components of the stress
  y_stress:   the y components of the stress
  theta_min:  start the radial search from this angle
  theta_max:  end the radial search at this angle
*/
void TACSConstitutive::writeBucklingEnvelope( const char * file_name, int npts,
					      const double gpt[], 
					      const TacsScalar x_stress[], 
					      const TacsScalar y_stress[],
					      double theta_min, double theta_max ){
  FILE * fp = fopen(file_name, "w");

  // The failure criteria solution tolerance
  double tol = 1e-8;

  // The max number of Newton iterations to use
  int max_newton_iters = 50;

  if (fp){
    int nstress = getNumStresses();

    // Compute an explicit form of the constitutive relationship
    int * ipiv = new int[ nstress ];
    TacsScalar * e = new TacsScalar[ nstress ];
    TacsScalar * f_sens = new TacsScalar[ nstress ];
    TacsScalar * x_strain = new TacsScalar[ nstress ];
    TacsScalar * y_strain = new TacsScalar[ nstress ];
    TacsScalar * C = new TacsScalar[ nstress*nstress ];

    for ( int i = 0; i < nstress; i++ ){
      memset(e, 0, nstress*sizeof(TacsScalar));
      e[i] = 1.0;
      calculateStress(gpt, e, &C[i*nstress]);
    }

    // Factor the constitutive matrix
    int info;
    LAPACKgetrf(&nstress, &nstress, C, &nstress, ipiv, &info);
    
    // Copy over the stress values and solve for the strain
    memcpy(x_strain, x_stress, nstress*sizeof(TacsScalar));
    memcpy(y_strain, y_stress, nstress*sizeof(TacsScalar));

    int one = 1;
    LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, 
		x_strain, &nstress, &info);
    LAPACKgetrs("N", &nstress, &one, C, &nstress, ipiv, 
		y_strain, &nstress, &info);

    fprintf(fp, "Variables = ");
    for ( int j = 0; j < nstress; j++ ){
      fprintf(fp, "s%d ", j);
    }
    fprintf(fp, "\n");

    // Estimate the smallest postive root of P and use
    // this as a starting point for Newton's method
    TacsScalar P = 1000.0;
    TacsScalar bval0 = 0.0, bval1 = 0.0, dbval1 = 0.0;
    memset(e, 0, nstress*sizeof(TacsScalar));
    buckling(e, &bval0);
    
    TacsScalar c = cos(theta_min), s = sin(theta_min);
    for ( int j = 0; j < nstress; j++ ){
      e[j] = P*(c*x_strain[j] + s*y_strain[j]);
    }

    buckling(e, &bval1);
    bucklingStrainSens(e, f_sens);
    for ( int j = 0; j < nstress; j++ ){
      dbval1 += f_sens[j]*(c*x_strain[j] + s*y_strain[j]);
    }

    // Compute the roots of a quadratic such that:
    // b(0) = bval0, b(P) = bval1, b'(P) = dbval1
    // b(P) = qa*P^2 + qb*P + qc - 1.0 = 0.0
    TacsScalar P1 = 0.0, P2 = 0.0;
    TacsScalar qa = (dbval1*P - (bval1 - bval0))/(P*P);
    TacsScalar qb = 2.0*(bval1 - bval0)/P - dbval1;
    TacsScalar qc = bval0 - 1.0;
    
    // Determine the minimum positive root
    FElibrary::solveQERoots(&P1, &P2, qa, qb, qc);

    // Assign the minimum positive root (if any) to P
    if (TacsRealPart(P1) > 0.0 && TacsRealPart(P2) > 0.0){
      P = P1;
      if (TacsRealPart(P2) < TacsRealPart(P1)){
	P = P2;
      }
    }
    else if (TacsRealPart(P1) > 0.0){
      P = P1;
    }
    else if (TacsRealPart(P2) > 0.0){
      P = P2;
    }
    
    for ( int k = 0; k < npts; k++ ){
      TacsScalar theta = theta_min + (theta_max - theta_min)*k/(npts-1);
      c = cos(theta);
      s = sin(theta);

      // Solve the equation: 
      // buckling(P*(c*x_strain + s*y_strain)) - 1.0 = 0           
      for ( int i = 0; i < max_newton_iters; i++ ){
	for ( int j = 0; j < nstress; j++ ){
	  e[j] = P*(c*x_strain[j] + s*y_strain[j]);
	}

	// Compute the buckling criterion and the derivative of the
	// buckling criterion w.r.t. the load parameter P
	TacsScalar bval, bvalSens;
	buckling(e, &bval);
	
	if (fabs(TacsRealPart(bval) - 1.0) < tol){
	  break;
	}

	bucklingStrainSens(e, f_sens);

	bvalSens = 0.0;
	for ( int j = 0; j < nstress; j++ ){
	  bvalSens += f_sens[j]*(c*x_strain[j] + s*y_strain[j]);
	}

	// Compute the Newton update to the boundary
	P = P - (bval - 1.0)/bvalSens;
      }

      // Compute the strain at the final point
      for ( int j = 0; j < nstress; j++ ){
	e[j] = P*(c*x_strain[j] + s*y_strain[j]);
      }

      // Compute the corresponding stress at the final point
      calculateStress(gpt, e, f_sens);

      // Print out the result to the file
      for ( int j = 0; j < nstress; j++ ){
	fprintf(fp, "%15.8f ", TacsRealPart(f_sens[j]));
      }
      fprintf(fp, "\n");
    }

    fclose(fp);

    delete [] ipiv;
    delete [] e;
    delete [] f_sens;
    delete [] x_strain;
    delete [] y_strain;
    delete [] C;
  }
}


/*
  This generates strains from a given stress.

  This results in more realistic stress levels and avoids issues
  realted to very large stress resulting in unrealistic failure
  prediction.  
*/
/*
void TestConstitutive::compute_strain( TacsScalar strain[], 
                                       const double pt[],
                                       const TacsScalar stress[] ){
  int nstress = con->getNumStresses();

  TacsScalar *C = new TacsScalar[ nstress*nstress ];
  memset(C, 0, nstress*nstress*sizeof(TacsScalar));

  for ( int k = 0; k < nstress; k++ ){
    memset(strain, 0, nstress*sizeof(TacsScalar));

    strain[k] = 1.0;
    con->calculateStress(pt, strain, &C[nstress*k]);
  }

  memcpy(strain, stress, nstress*sizeof(TacsScalar));

  // Solve for the strain given the stresses
  int *ipiv = new int[ nstress ];
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
*/
/*
  Test the sensitivity of the strain to a failure 
*/
/*
int TestConstitutive::testFailStrainSens( const double pt[] ){
  int nstress = con->getNumStresses();

  TacsScalar *strain = new TacsScalar[ nstress ];
  TacsScalar *sens = new TacsScalar[ nstress ];
  TacsScalar *strainCopy = new TacsScalar[ nstress ];
  TacsScalar *sensApprox = new TacsScalar[ nstress ];

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
    sensApprox[k] = TacsImagPart(forward)/dh;
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

  int test_failflag = (max_err > test_failatol || max_rel > test_failrtol);

  if (test_print_level > 0){
    fprintf(stderr, 
	    "Testing the failure sensivity for %s.\n",
	    con->constitutiveName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the failure criteria w.r.t. \
the components of strain\n");
    print_error_components(stderr, "failStrainSens", 
                           sens, sensApprox, nstress);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] sens;
  delete [] strain;
  delete [] strainCopy;
  delete [] sensApprox;

  return test_failflag;
}
*/
/*
int TestConstitutive::testFailDVSens( const double pt[], int ndvs ){ 
  int nstress = con->getNumStresses();

  TacsScalar *strain = new TacsScalar[ nstress ];
  TacsScalar *stress = new TacsScalar[ nstress ];

  // Randomly generate the stresses for failure criteria calculations
  generate_random_array(stress, nstress);

  // Determine the corresponding strains
  compute_strain(strain, pt, stress);
  delete [] stress;

  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int *dv_nums = new int[ ndvs ];
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

  TacsScalar *dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (test_print_level){
    fprintf(stderr, 
	    "Testing the sensivity of the failure criteria w.r.t. the \
design variables for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int test_failflag = 0;

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

    TacsScalar failSensApprox = TacsImagPart(forward)/dh;
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
   
    if (test_print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (test_print_level > 1){
      fprintf(stderr, 
              "The sensitivity failure criteria w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "failDVSens", 
                             &failSens, &failSensApprox, 1);
    }
    if (test_print_level){ fprintf(stderr, "\n"); }

    test_failflag = (test_failflag || (max_err > test_failatol || max_rel > test_failrtol));
  }

  delete [] strain;
  delete [] dvs;
  delete [] dv_nums;

  return test_failflag;
}
*/
/*
int TestConstitutive::testMassDVSens( const double pt[] ){ 
  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int ndvs = con->getNumDesignVars();
  int *dv_nums = new int[ ndvs ];
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

  TacsScalar *dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (test_print_level){
    fprintf(stderr, 
	    "Testing the sensivity w.r.t. the design variables \
for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int test_failflag = 0;

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
    TacsScalar massSensApprox = TacsImagPart(forward[0])/dh;
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
   
    if (test_print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (test_print_level > 1){
      fprintf(stderr, 
              "The sensitivity mass w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "massDVSens", 
                             massSens, &massSensApprox, 1);
    }
    if (test_print_level){ fprintf(stderr, "\n"); }

    test_failflag = (test_failflag || (max_err > test_failatol || max_rel > test_failrtol));
  }

  delete [] dvs;
  delete [] dv_nums;

  return test_failflag;
}
*/
/*
  Test the sensitivity of the strain to a buckling 
*/
/*
int TestConstitutive::testBucklingStrainSens(){
  int nstress = con->getNumStresses();

  TacsScalar *strain = new TacsScalar[ nstress ];
  TacsScalar *sens = new TacsScalar[ nstress ];
  TacsScalar *strainCopy = new TacsScalar[ nstress ];
  TacsScalar *sensApprox = new TacsScalar[ nstress ];

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

  int test_failflag = (max_err > test_failatol || max_rel > test_failrtol);

  if (test_print_level > 0){
    fprintf(stderr, 
	    "Testing the buckling sensivity for %s.\n",
	    con->constitutiveName());
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
	    max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
	    max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    fprintf(stderr, 
            "The sensitivity of the buckling criteria w.r.t. \
the components of strain\n");
    print_error_components(stderr, "bucklingStrainSens", 
                           sens, sensApprox, nstress);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  delete [] sens;
  delete [] strain;
  delete [] strainCopy;
  delete [] sensApprox;

  return test_failflag;
}
*/
/*
int TestConstitutive::testBucklingDVSens(){ 
  int nstress = con->getNumStresses();

  TacsScalar *strain = new TacsScalar[ nstress ];
  TacsScalar *stress = new TacsScalar[ nstress ];

  // Randomly generate the stresses for buckling criteria calculations
  generate_random_array(stress, nstress);

  // Determine the corresponding strains
  double pt[3] = {0.0, 0.0, 0.0};
  compute_strain(strain, pt, stress);
  delete [] stress;

  // Next, determine the design variables to use.
  // Get the design variables associated with this constitutive class.
  int ndvs = con->getNumDesignVars();
  int *dv_nums = new int[ ndvs ];
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

  TacsScalar *dvs = new TacsScalar[ max_dv ];
  con->getDesignVars(dvs, max_dv);

  if (test_print_level){
    fprintf(stderr, 
	    "Testing the sensivity of the buckling criteria w.r.t. the \
design variables for constitutive class %s.\n",
	    con->constitutiveName());
  }

  int test_failflag = 0;

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
   
    if (test_print_level > 0){
      fprintf(stderr, "Max Err dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_err, max_err_index);
      fprintf(stderr, "Max REr dv %3d: %10.4e in component %d.\n",
	      dv_nums[k], max_rel, max_rel_index);
    }
    // Print the error if required
    if (test_print_level > 1){
      fprintf(stderr, 
              "The sensitivity buckling criteria w.r.t. dv %d is \n",
              dv_nums[k]);
      print_error_components(stderr, "bucklingDVSens", 
                             &bucklingSens, &bucklingSensApprox, 1);
    }
    if (test_print_level){ fprintf(stderr, "\n"); }

    test_failflag = (test_failflag || (max_err > test_failatol || max_rel > test_failrtol));
  }

  delete [] strain;
  delete [] dvs;
  delete [] dv_nums;

  return test_failflag;
}
*/
