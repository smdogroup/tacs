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

#include "FSDTStiffness.h"
#include "tacslapack.h"

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
         TacsRealPart(A[0]), TacsRealPart(A[1]), TacsRealPart(A[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(A[1]), TacsRealPart(A[3]), TacsRealPart(A[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(A[2]), TacsRealPart(A[4]), TacsRealPart(A[5]));

  printf("\nThe B matrix: \n");
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(B[0]), TacsRealPart(B[1]), TacsRealPart(B[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(B[1]), TacsRealPart(B[3]), TacsRealPart(B[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(B[2]), TacsRealPart(B[4]), TacsRealPart(B[5]));

  printf("\nThe D matrix: \n");
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(D[0]), TacsRealPart(D[1]), TacsRealPart(D[2]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(D[1]), TacsRealPart(D[3]), TacsRealPart(D[4]));
  printf("%20.4f %20.4f %20.4f \n", 
	 TacsRealPart(D[2]), TacsRealPart(D[4]), TacsRealPart(D[5]));

  printf("\nThe Ashear matrix: \n");
  printf("%20.4f %20.4f \n", 
	 TacsRealPart(As[0]), TacsRealPart(As[1]));
  printf( "%20.4f %20.4f \n", 
	 TacsRealPart(As[1]), TacsRealPart(As[2]));
}
