#include <stdlib.h>
#include <math.h>
#include "tacslapack.h"
#include "TACSConstitutive.h"
#include "FElibrary.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*! 
  The maximum number of stress components. This maximum applies to 
  testing only and can be increased if necessary.
*/

static const int MAX_STRESSES = 20;

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
	
	if (fabs(fail - 1.0) < tol){
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
	fprintf(fp, "%15.8f ", RealPart(f_sens[j]));
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
    if (P1 > 0.0 && P2 > 0.0){
      P = P1;
      if (P2 < P1){
	P = P2;
      }
    }
    else if (P1 > 0.0){
      P = P1;
    }
    else if (P2 > 0.0){
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
	
	if (fabs(bval - 1.0) < tol){
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
	fprintf(fp, "%15.8f ", RealPart(f_sens[j]));
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
