#include "FSDTStiffness.h"
#include "tacslapack.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * FSDTStiffness::constName = "FSDTStiffness";

FSDTStiffness::FSDTStiffness(){
  is_param_variable = 0;
  transform_type = NATURAL;
  axis[0] = axis[1] = axis[2] = 0.0;
  axis[0] = 1.0;
}

const char * FSDTStiffness::constitutiveName() const { 
  return constName; 
}

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
int FSDTStiffness::getNumStresses() const { 
  return NUM_STRESSES; 
}

/*
  Calculate the stress based on the strain
*/
void FSDTStiffness::calculateStress( const double gpt[], 
				     const TacsScalar strain[],
				     TacsScalar stress[] ){
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);
  calcStress(A, B, D, As, strain, stress);
}

/*
  Calculate the derivative of the stress w.r.t. the design index
*/
void FSDTStiffness::calculateStressDVSens( int dvNum, 
                                           const double gpt[], 
					   const TacsScalar strain[],
					   TacsScalar stress[] ){
  TacsScalar sA[6], sB[6], sD[6], sAs[3];
  getStiffnessDVSens(dvNum, gpt, sA, sB, sD, sAs);
  calcStress(sA, sB, sD, sAs, strain, stress);
}

/*
  Compute the derivative of the stiffness matrix w.r.t. the design
  variable 
*/
TacsScalar FSDTStiffness::getStiffnessDVSens( int dvNum, 
                                              const double gpt[],
                                              TacsScalar sA[], 
                                              TacsScalar sB[],
                                              TacsScalar sD[], 
                                              TacsScalar sAs[] ){
  sA[0] = sA[1] = sA[2] = sA[3] = sA[4] = sA[5] = 0.0;
  sB[0] = sB[1] = sB[2] = sB[3] = sB[4] = sB[5] = 0.0;
  sD[0] = sD[1] = sD[2] = sD[3] = sD[4] = sD[5] = 0.0;
  sAs[0] = sAs[1] = sAs[2] = 0.0;

  return 0.0;
}

/*
  Compute the inner product of the derivative of the stress times the
  input strain

  Note that this function might not be used as often as for other
  constitutive classes since we also have the rotational stiffness
  that is not stored in the stiffness matrix itself (although perhaps
  it should be...)
*/
void FSDTStiffness::addStressDVSens( const double pt[], 
				     const TacsScalar strain[], 
				     TacsScalar alpha, const TacsScalar psi[], 
				     TacsScalar dvSens[], int dvLen ){
  TacsScalar sA[6], sB[6], sD[6], sAs[3];

  int ndvs = getNumDesignVars();
  if (ndvs <= 256){
    int dvNums[256];
    ndvs = 0;
    getDesignVarNums(dvNums, &ndvs, 256); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	getStiffnessDVSens(dvNums[k], pt, sA, sB, sD, sAs);

	TacsScalar s[8];
	calcStress(sA, sB, sD, sAs, strain, s);
	
	dvSens[dvNums[k]] += alpha*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
				    s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5] +
				    s[6]*psi[6] + s[7]*psi[7]);
      }
    }
  }
  else {
    int *dvNums = new int[ ndvs ];
    int ndvs = 0;
    ndvs = getDesignVarNums(dvNums, &ndvs, ndvs); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	getStiffnessDVSens(dvNums[k], pt, sA, sB, sD, sAs);

	TacsScalar s[8];
	calcStress(sA, sB, sD, sAs, strain, s);
	
	dvSens[dvNums[k]] += alpha*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
				    s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5] +
				    s[6]*psi[6] + s[7]*psi[7]);
      }
    }

    delete [] dvNums;
  }
}

/*
  Compute the inner product of the derivative of the stress times the
  input strain

  Note that this function might not be used as often as for other
  constitutive classes since we also have the rotational stiffness
  that is not stored in the stiffness matrix itself (although perhaps
  it should be...)
*/
void FSDTStiffness::addStiffnessDVSens( const double pt[], TacsScalar alpha, 
					const TacsScalar psi[], const TacsScalar phi[], 
					const TacsScalar psi_rot, const TacsScalar phi_rot,
					TacsScalar dvSens[], int dvLen ){
  TacsScalar sA[6], sB[6], sD[6], sAs[3];

  int ndvs = getNumDesignVars();
  if (ndvs <= 256){
    int dvNums[256];
    ndvs = 0;
    getDesignVarNums(dvNums, &ndvs, 256); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar krot = getStiffnessDVSens(dvNums[k], pt, sA, sB, sD, sAs);

	TacsScalar s[8];
	calcStress(sA, sB, sD, sAs, phi, s);
	
	dvSens[dvNums[k]] += alpha*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
				    s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5] +
				    s[6]*psi[6] + s[7]*psi[7] + 
				    krot*psi_rot*phi_rot);
      }
    }
  }
  else {
    int *dvNums = new int[ ndvs ];
    int ndvs = 0;
    ndvs = getDesignVarNums(dvNums, &ndvs, ndvs); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar krot = getStiffnessDVSens(dvNums[k], pt, sA, sB, sD, sAs);

	TacsScalar s[8];
	calcStress(sA, sB, sD, sAs, phi, s);
	
	dvSens[dvNums[k]] += alpha*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] + 
				    s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5] +
				    s[6]*psi[6] + s[7]*psi[7] + 
				    krot*psi_rot*phi_rot);
      }
    }

    delete [] dvNums;
  }
}

/*
  Get the derivative of the pointwise mass w.r.t. the given design
  variable 
*/
void FSDTStiffness::pointwiseMassDVSens( int dvNum,
                                         const double gpt[],
                                         TacsScalar massDVSens[] ){
  massDVSens[0] = 0.0;
  massDVSens[1] = 0.0;
  massDVSens[2] = 0.0;
}

/*
  Add the derivative of the mass with respect to the design variables
  to the given vector.
*/
void FSDTStiffness::addPointwiseMassDVSens( const double pt[], 
					    const TacsScalar alpha[],
					    TacsScalar dvSens[], int dvLen ){
  TacsScalar massDVSens[3];

  int ndvs = getNumDesignVars();
  if (ndvs <= 256){
    int dvNums[256];
    ndvs = 0;
    getDesignVarNums(dvNums, &ndvs, 256); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	pointwiseMassDVSens(dvNums[k], pt, massDVSens);
	
	dvSens[dvNums[k]] += (alpha[0]*massDVSens[0] +
			      alpha[1]*massDVSens[1] +
			      alpha[2]*massDVSens[2]);
      }
    }
  }
  else {
    int *dvNums = new int[ ndvs ];
    int ndvs = 0;
    ndvs = getDesignVarNums(dvNums, &ndvs, ndvs); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	pointwiseMassDVSens(dvNums[k], pt, massDVSens);
	
	dvSens[dvNums[k]] += (alpha[0]*massDVSens[0] +
			      alpha[1]*massDVSens[1] +
			      alpha[2]*massDVSens[2]);
      }
    }

    delete [] dvNums;
  }  
}

/*
  Add the product of the failure constraint with the alpha
*/
void FSDTStiffness::addFailureDVSens( const double pt[], 
				      const TacsScalar strain[],
				      TacsScalar alpha, 
				      TacsScalar dvSens[], int dvLen ){
  int ndvs = getNumDesignVars();
  if (ndvs <= 256){
    int dvNums[256];
    ndvs = 0;
    getDesignVarNums(dvNums, &ndvs, 256); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar failSens;
	failureDVSens(dvNums[k], pt, strain, &failSens);
	dvSens[dvNums[k]] += alpha*failSens;
      }
    }
  }
  else {
    int *dvNums = new int[ ndvs ];
    int ndvs = 0;
    ndvs = getDesignVarNums(dvNums, &ndvs, ndvs); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar failSens;
	failureDVSens(dvNums[k], pt, strain, &failSens);
	dvSens[dvNums[k]] += alpha*failSens;
      }
    }

    delete [] dvNums;
  }
}
 
/*
  Add the derivative of buckling criteria with respect to the design
  variables 
*/
void FSDTStiffness::addBucklingDVSens( const TacsScalar strain[], 
				       TacsScalar alpha,
				       TacsScalar dvSens[], int dvLen ){
  int ndvs = getNumDesignVars();
  if (ndvs <= 256){
    int dvNums[256];
    ndvs = 0;
    getDesignVarNums(dvNums, &ndvs, 256); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar bucklingSens;
	bucklingDVSens(dvNums[k], strain, &bucklingSens);
	dvSens[dvNums[k]] += alpha*bucklingSens;
      }
    }
  }
  else {
    int *dvNums = new int[ ndvs ];
    int ndvs = 0;
    ndvs = getDesignVarNums(dvNums, &ndvs, ndvs); 

    for ( int k = 0; k < ndvs; k++ ){
      if (dvNums[k] >= 0 && dvNums[k] < dvLen){
	TacsScalar bucklingSens;
	bucklingDVSens(dvNums[k], strain, &bucklingSens);
	dvSens[dvNums[k]] += alpha*bucklingSens;
      }
    }

    delete [] dvNums;
  }
}

/*
  Print the stiffness information from this object
*/
void FSDTStiffness::printInfo(){
  double gpt[3] = {0.0, 0.0, 0.0};

  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(gpt, A, B, D, As);

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

/*
  Determine the strain that corresponds to the input stress at the
  present point

  input:
  pt: the parametric point to evaluate the constitutive matrix
  s: the stress

  output:
  e: the strain components
*/
void FSDTStiffness::computeStrainFromStress( const double pt[], 
                                             const TacsScalar s[],
                                             TacsScalar e[] ){
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffness(pt, A, B, D, As);

  TacsScalar C[64];
  memset(C, 0, 64*sizeof(TacsScalar));

  // Set the A matrix
  C[0]  = A[0];  C[1]  = A[1];  C[2]  = A[2];
  C[8]  = A[1];  C[9]  = A[3];  C[10] = A[4];
  C[16] = A[2];  C[17] = A[4];  C[18] = A[5];

  // Set the B matrix
  C[3]  = B[0];  C[4]  = B[1];  C[5]  = B[2];
  C[11] = B[1];  C[12] = B[3];  C[13] = B[4];
  C[19] = B[2];  C[20] = B[4];  C[21] = B[5];

  C[24] = B[0];  C[25] = B[1];  C[26] = B[2];
  C[32] = B[1];  C[33] = B[3];  C[34] = B[4];
  C[40] = B[2];  C[41] = B[4];  C[42] = B[5];

  // Set the D matrix
  C[27] = D[0];  C[28] = D[1];  C[29] = D[2];
  C[35] = D[1];  C[36] = D[3];  C[37] = D[4];
  C[43] = D[2];  C[44] = D[4];  C[45] = D[5];

  // Set the As matrix
  C[54] = As[0];  C[55] = As[1];
  C[62] = As[1];  C[63] = As[2];

  int n = 8;
  int one = 1, ipiv[8];
  int info;
  memcpy(e, s, n*sizeof(TacsScalar));
  LAPACKgesv(&n, &one, C, &n, ipiv, e, &n, &info);
}

/*
  Evaluate the derivative of the strain with respect to the design
  variable dvNum.

  input:
  dvNum: the design variable number
  pt: the parametric point

  output:
  es: the sensitivity of the strain components

  The strain is given by:
  e = C^{-1} s

  The derivative is determined as follows:
  dC/dx * e + C * de/dx = 0
  de/ex = - C^{-1} * dC/dx * e
*/
void FSDTStiffness::computeStrainFromStressDVSens( int dvNum, 
                                                   const double pt[], 
                                                   const TacsScalar e[],
                                                   TacsScalar es[] ){
  // Evaluate the derivative of the stiffness w.r.t. the design variable
  TacsScalar A[6], B[6], D[6], As[3];
  getStiffnessDVSens(dvNum, pt, A, B, D, As);

  // Compute es = - dC/dx * e
  calcStress(A, B, D, As, e, es);
  for ( int i = 0; i < 8; i++ ){
    es[i] *= -1;
  }

  // Compute the stiffness
  getStiffness(pt, A, B, D, As);

  TacsScalar C[64];
  memset(C, 0, 64*sizeof(TacsScalar));

  // Set the A matrix
  C[0]  = A[0];  C[1]  = A[1];  C[2]  = A[2];
  C[8]  = A[1];  C[9]  = A[3];  C[10] = A[4];
  C[16] = A[2];  C[17] = A[4];  C[18] = A[5];

  // Set the B matrix
  C[3]  = B[0];  C[4]  = B[1];  C[5]  = B[2];
  C[11] = B[1];  C[12] = B[3];  C[13] = B[4];
  C[19] = B[2];  C[20] = B[4];  C[21] = B[5];

  C[24] = B[0];  C[25] = B[1];  C[26] = B[2];
  C[32] = B[1];  C[33] = B[3];  C[34] = B[4];
  C[40] = B[2];  C[41] = B[4];  C[42] = B[5];

  // Set the D matrix
  C[27] = D[0];  C[28] = D[1];  C[29] = D[2];
  C[35] = D[1];  C[36] = D[3];  C[37] = D[4];
  C[43] = D[2];  C[44] = D[4];  C[45] = D[5];

  // Set the As matrix
  C[54] = As[0];  C[55] = As[1];
  C[62] = As[1];  C[63] = As[2];

  int n = 8;
  int one = 1, ipiv[8];
  int info;

  // Evaluate es = de/dx = - C^{-1} * dC/dx * e
  LAPACKgesv(&n, &one, C, &n, ipiv, es, &n, &info);
}

/*
  Evaluate the failure criterion for the function based 
  on the input value of the stress

  input:
  pt: the parametric point
  s: the stress
  
  returns:
  the failure constraint value
*/
TacsScalar FSDTStiffness::computeFailure( const double pt[], 
                                          const TacsScalar s[] ){
  TacsScalar e[8];
  computeStrainFromStress(pt, s, e);

  TacsScalar fail;
  failure(pt, e, &fail);
  return fail;
}

/*
  Evaluate the buckling envelope constraint based on the
  the input value of the stress

  input:
  pt: the parametric point
  s: the stress
  
  returns:
  the buckling constraint value
*/
TacsScalar FSDTStiffness::computeBuckling( const double pt[], 
                                           const TacsScalar s[] ){
  TacsScalar e[8];
  computeStrainFromStress(pt, s, e);

  TacsScalar fail;
  buckling(e, &fail);
  return fail;
}

/*
  Evaluate the derivative of the failure function based on 
  the input value of the stress

  input:
  pt: the parametric point
  s: the input compoments of the stress
  numDVs: the number of design variables

  output:
  fdvSens: the sensitivity w.r.t. the design variables
*/
void FSDTStiffness::computeFailureDVSens( const double pt[], const TacsScalar s[], 
                                          TacsScalar fdvSens[], int numDVs ){
  TacsScalar e[8];
  computeStrainFromStress(pt, s, e);

  for ( int i = 0; i < numDVs; i++ ){
    if (ownsDesignVar(i)){
      failureDVSens(i, pt, e, &fdvSens[i]);

      // Compute the derivative of the failure function w.r.t. the
      // strain compoments
      TacsScalar fes[8];
      failureStrainSens(pt, e, fes); 

      // Compute the derivative of the strain w.r.t. the design variable
      TacsScalar es[8];
      computeStrainFromStressDVSens(i, pt, e, es);
      fdvSens[i] += (fes[0]*es[0] + fes[1]*es[1] + fes[2]*es[2] +
                     fes[3]*es[3] + fes[4]*es[4] + fes[5]*es[5] +
                     fes[6]*es[6] + fes[7]*es[7]);
    }
  }
}

/*
  Evaluate the derivative of the buckling function w.r.t. the
  design variables

  input:
  pt: the parametric point
  s: the input components of the stress
  numDVs: the number of design variables
 
  ouput: 
  fdvSens: the derivative of the buckling function w.r.t. the
  design variables
*/
void FSDTStiffness::computeBucklingDVSens( const double pt[], const TacsScalar s[],
                                           TacsScalar fdvSens[], int numDVs ){
  TacsScalar e[8];
  computeStrainFromStress(pt, s, e);

  for ( int i = 0; i < numDVs; i++ ){
    if (ownsDesignVar(i)){
      bucklingDVSens(i, e, &fdvSens[i]);

      // Compute the derivative of the buckling function w.r.t. the
      // strain compoments
      TacsScalar fes[8];
      bucklingStrainSens(e, fes); 

      // Compute the derivative of the strain w.r.t. the design variable
      TacsScalar es[8];
      computeStrainFromStressDVSens(i, pt, e, es);
      
      // Add this constribution to the sensitivity
      fdvSens[i] += (fes[0]*es[0] + fes[1]*es[1] + fes[2]*es[2] +
                     fes[3]*es[3] + fes[4]*es[4] + fes[5]*es[5] +
                     fes[6]*es[6] + fes[7]*es[7]);
    }
  }
}
