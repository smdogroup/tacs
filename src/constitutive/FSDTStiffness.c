#include "FSDTStiffness.h"
#include "tacslapack.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * FSDTStiffness::constName = "FSDTStiffness";

/*
  Construct the FSDTStiffness object
*/
FSDTStiffness::FSDTStiffness(){
  transform_type = NATURAL;
  axis[0] = axis[1] = axis[2] = 0.0;
  axis[0] = 1.0;
}

/*
  Return the constitutive name
*/
const char * FSDTStiffness::constitutiveName(){ 
  return constName; 
}

/*
  Set the default drilling regularization value
*/
double FSDTStiffness::DRILLING_REGULARIZATION = 10.0;

/*
  Set the drilling stiffness regularization parameter
*/ 
void FSDTStiffness::setDrillingRegularization( double kval ){
  DRILLING_REGULARIZATION = kval;
}

/*
  Return the number of stresses associated with this constitutive
  object 
*/
int FSDTStiffness::getNumStresses(){ 
  return NUM_STRESSES; 
}

/*
  Return the number of mass moments defined by this group of
  objects
*/
int FSDTStiffness::getNumMassMoments(){
  return NUM_MASS_MOMENTS;
}

/*
  Given the parametric point within the element and the value of the 
  strain, compute the stress

  input:
  pt:   the parametric point
  e:    the components of the strain

  output:
  s:    the components of the stress
*/
void FSDTStiffness::calculateStress( const double pt[],
                                     const TacsScalar e[],
                                     TacsScalar s[] ){
  // Compute the stiffness matrix (ignore the in-plane rotation penalty)
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(pt, A, B, D, As);

  // Cmopute the stress
  calculateStress(A, B, D, As, e, s);
}

/*
  Given the parametric point, the comonents of the strain, and a
  vector psi equal in length to the size of the strain, add the
  derivative of the product with the stiffness with respect to the
  design variables to the input vector.

  In other words, compute the derivative:

  fdvSens += alpha*d(psi^{T}*C*e)/dx

  input:
  pt:       the parametric point
  e:        the components of the strain
  alpha:    a scalar multiplier
  psi:      the components of the multiplying vector
  dvLen:    the length of the design variable array

  output:
  fdvSens:  the design variable vector
*/
void FSDTStiffness::addStressDVSens( const double pt[], const TacsScalar e[],
                                     TacsScalar alpha, const TacsScalar psi[], 
                                     TacsScalar fdvSens[], int dvLen ){
  // Scale the input vector psi by the scalar alpha
  TacsScalar p[8];
  if (alpha != 1.0){
    p[0] = alpha*psi[0];
    p[1] = alpha*psi[1];
    p[2] = alpha*psi[2];
    p[3] = alpha*psi[3];
    p[4] = alpha*psi[4];
    p[5] = alpha*psi[5];
    p[6] = alpha*psi[6];
    p[7] = alpha*psi[7];
  }

  // Add the derivative of the product to the design vector
  TacsScalar rotPsi = 0.0;
  addStiffnessDVSens(pt, e, p, rotPsi, fdvSens, dvLen);
}

/*
  Print the stiffness information from this object
*/
void FSDTStiffness::printStiffness( const double pt[] ){
  double zero[3] = {0.0, 0.0, 0.0};

  // Evaluate the stiffness
  TacsScalar A[6], B[6], D[6], As[3];
  if (pt){ 
    getStiffness(pt, A, B, D, As);
  }
  else {
    getStiffness(zero, A, B, D, As);
  }

  // Print out the stiffness matrices
  printf("\nThe A matrix: \n");
  printf("%20.4f %20.4f %20.4f \n", 
         RealPart(A[0]), RealPart(A[1]), RealPart(A[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(A[1]), RealPart(A[3]), RealPart(A[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(A[2]), RealPart(A[4]), RealPart(A[5]));

  printf("\nThe B matrix: \n");
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(B[0]), RealPart(B[1]), RealPart(B[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(B[1]), RealPart(B[3]), RealPart(B[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(B[2]), RealPart(B[4]), RealPart(B[5]));

  printf("\nThe D matrix: \n");
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(D[0]), RealPart(D[1]), RealPart(D[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(D[1]), RealPart(D[3]), RealPart(D[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 RealPart(D[2]), RealPart(D[4]), RealPart(D[5]));

  printf("\nThe Ashear matrix: \n");
  printf("%20.4f %20.4f \n", 
	 RealPart(As[0]), RealPart(As[1]));
  printf( "%20.4f %20.4f \n", 
	 RealPart(As[1]), RealPart(As[2]));
}
