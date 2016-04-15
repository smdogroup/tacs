#include "RigidBody.h"
#include "TACSElementAlgebra.h"

/*
  Rigid-body dynamics routines for TACS

  Copyright (c) 2015 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  Write the relative error for components of a vector for a
  finite-difference or complex-step study to a given file

  input:
  fp:         the output file 
  descript:   description of the component
  a:          analytic values
  b:          set of finite-difference/complex-step values
  size:       the number of values
  rel_err:    relative error tolerance
*/
void writeErrorComponents( FILE * fp, const char * descript,
			   TacsScalar *a, TacsScalar *fd, 
			   int size, double rel_err ){
  int print_flag = 1;
  for ( int i = 0; i < size; i++ ){
    double rel = 0.0;
    if (a[i] != 0.0){
      rel = fabs(RealPart((a[i] - fd[i])/a[i]));
    }
    else {
      rel = fabs(RealPart((a[i] - fd[i])));
    }

    if (rel > rel_err || a[i] != a[i] || fd[i] != fd[i]){
      if (print_flag){
	fprintf(fp, "%*s[   ] %15s %15s %15s\n",
		(int)strlen(descript), "Val", "Analytic", 
		"Approximate", "Rel. Error");
	print_flag = 0;
      }
      if (a[i] != 0.0){
	fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n", 
		descript, i, RealPart(a[i]), RealPart(fd[i]), 
		fabs(RealPart((a[i] - fd[i])/a[i])));
      }
      else {
	fprintf(fp, "%s[%3d] %15.6e %15.6e\n", 
		descript, i, RealPart(a[i]), RealPart(fd[i]));
      }
    }
  }  
}

/*
  A reference coordinate frame

  This code generates a reference frame from three vectors.  The
  vectors define a primary direction and a secondary direction which
  are used to three orthonormal/right-handed vectors.

  input:
  r0:    the base point of the coordinate frame
  r1:    (r1 - r0)/||r1 - r0|| forms the first direction
  r2:    direction perpendicular to (r1 - r0) used as second direction
*/
TACSRefFrame::TACSRefFrame( TACSGibbsVector *_r0, 
			    TACSGibbsVector *_r1,
			    TACSGibbsVector *_r2 ){
  r0 = _r0;
  r1 = _r1;
  r2 = _r2;
  r0->incref();
  r1->incref();
  r2->incref();
  initialize();
}

/*
  Free the allocated memory
*/
TACSRefFrame::~TACSRefFrame(){
  r0->decref();
  r1->decref();
  r2->decref();
}

/*
  Retrieve the rotation matrix from the global to local coordinates
*/
void TACSRefFrame::getRotation( const TacsScalar **_C ){
  *_C = C;
}

/*
  Update the coordinate frame to reflect any changes to the vectors
  that form the basis of the initial transformation matrix. 

  This should be called once during the initial condition
  calculation. This will enable the modification of reference frames
  based solely on moving points.
*/
void TACSRefFrame::initialize(){
  // Compute the initial position based on v1
  const TacsScalar *x0, *x1, *x2;
  r0->getVector(&x0);
  r1->getVector(&x1);
  r2->getVector(&x2);

  // Compute the first row of the transformation matrix
  TacsScalar d1[3];
  d1[0] = x1[0] - x0[0];
  d1[1] = x1[1] - x0[1];  
  d1[2] = x1[2] - x0[2];

  // Compute the first row of the transformation matrix
  TacsScalar C1norm = sqrt(vecDot(d1, d1));
  TacsScalar invC1norm = 1.0/C1norm;
  C[0] = invC1norm*d1[0];
  C[1] = invC1norm*d1[1];
  C[2] = invC1norm*d1[2];
  
  // Compute the second vector for the off-axis direction
  TacsScalar d2[3];
  d2[0] = x2[0] - x0[0];  
  d2[1] = x2[1] - x0[1];  
  d2[2] = x2[2] - x0[2];

  // Orthogonalize the first and second vectors
  TacsScalar dot = vecDot(&C[0], d2);
  TacsScalar s2[3];
  s2[0] = d2[0] - dot*C[0];
  s2[1] = d2[1] - dot*C[1];
  s2[2] = d2[2] - dot*C[2];

  // Normalize the vector
  TacsScalar C2norm = sqrt(vecDot(s2, s2));
  TacsScalar invC2norm = 1.0/C2norm;
  C[3] = invC2norm*s2[0];
  C[4] = invC2norm*s2[1];
  C[5] = invC2norm*s2[2];

  // Compute the third basis vector
  crossProduct(1.0, &C[0], &C[3], &C[6]);

  // Code to compute the derivative
  // The derivative of C1 w.r.t. d1 = (x1 - x0)
  vecNormDeriv(C1norm, d1, dC1d1);

  // Compute the derivative of C2 w.r.t. d2
  TacsScalar dC2ds2[9];
  vecNormDeriv(C2norm, s2, dC2ds2);

  TacsScalar ds2d2[9];
  ds2d2[0] = 1.0 - C[0]*C[0];
  ds2d2[1] = -C[0]*C[1];
  ds2d2[2] = -C[0]*C[2];

  ds2d2[3] = -C[1]*C[0];
  ds2d2[4] = 1.0 - C[1]*C[1];
  ds2d2[5] = -C[1]*C[2];

  ds2d2[6] = -C[2]*C[0];
  ds2d2[7] = -C[2]*C[1];
  ds2d2[8] = 1.0 - C[2]*C[2];

  // Compute dC2d2 = dC2ds2*ds2d2
  matMatMult(dC2ds2, ds2d2, dC2d2);

  // Compute the derivative dC2d1 = dC2ds2*ds2dC1*dC1d1
  // First find the derivative of s2 w.r.t. C1
  TacsScalar ds2dC1[9];
  ds2dC1[0] = -(dot + C[0]*d2[0]);
  ds2dC1[1] = -C[0]*d2[1];
  ds2dC1[2] = -C[0]*d2[2];

  ds2dC1[3] = -C[1]*d2[0];
  ds2dC1[4] = -(dot + C[1]*d2[1]);
  ds2dC1[5] = -C[1]*d2[2];

  ds2dC1[6] = -C[2]*d2[0];
  ds2dC1[7] = -C[2]*d2[1];
  ds2dC1[8] = -(dot + C[2]*d2[2]);

  // Compute the product to find the derivative of C2 w.r.t. d1
  TacsScalar tmp[9];
  matMatMult(dC2ds2, ds2dC1, tmp);
  matMatMult(tmp, dC1d1, dC2d1);
}

/*
  Set the design variables into the three points which define the
  reference frame.

  input:
  dvs:     the design variable values
  numDVs:  the number of design vars/array length
*/
void TACSRefFrame::setDesignVars( const TacsScalar *dvs, 
				  int numDVs ){
  r0->setDesignVars(dvs, numDVs);
  r1->setDesignVars(dvs, numDVs);
  r2->setDesignVars(dvs, numDVs);
  initialize();
}

/*
  Get the design variable values from the points that define the
  reference frame.

  input:
  numDVs:  the number of design vars/array length

  output:
  dvs:     the design variable values
*/
void TACSRefFrame::getDesignVars( TacsScalar *dvs, int numDVs ){
  r0->getDesignVars(dvs, numDVs);
  r1->getDesignVars(dvs, numDVs);
  r2->getDesignVars(dvs, numDVs);
}

/*
  Add the derivative of the inner product of the rotation matrix with
  two input vectors to an array.

  The code computes the following:

  fdvSens = d(psi^{T}*C*phi)/dx

  The code adds the result to the design variable locations of the r0,
  r1, and r2 points. The derivatives of the rotation matrix with
  respect to the points are computed when the point locations are
  updated. These derivatives are stored in dC1d1, dC2d1, and dC2d2.
  These represent the derivatives of the first and second rows of the
  rotation matrix with respect to the difference between (r1 - r0) and
  (r2 - r0), respectively.

  input:
  numDVs:   the number of design variables in the array
  psi:      the pre-multiplying vector
  phi:      the post-multiplying vector

  output:
  fdvSens:  the derivatives are added to this array
*/
void TACSRefFrame::addRotationAdjResProduct( TacsScalar fdvSens[],
					     int numDVs,
					     const TacsScalar psi[],
					     const TacsScalar phi[] ){
  // Temporary vector
  TacsScalar tmp[3];

  // Components that store the contributions to each point
  TacsScalar dx0[3], dx1[3], dx2[3];

  // Compute the contribution from the derivative of the first row
  matMultTrans(dC1d1, phi, dx1);
  vecScale(psi[0], dx1);

  // Compute the contributions from the derviative of the second row
  matMultTrans(dC2d1, phi, tmp);
  vecAxpy(psi[1], tmp, dx1);

  matMultTrans(dC2d2, phi, dx2);
  vecScale(psi[1], dx2);

  // Compute the contribution from the third row such that c3 = c1 x c2
  crossProduct(-psi[2], &C[0], phi, tmp);
  matMultTransAdd(dC2d1, tmp, dx1);
  matMultTransAdd(dC2d2, tmp, dx2);

  crossProduct(psi[2], &C[3], phi, tmp);
  matMultTransAdd(dC1d1, tmp, dx1);
  
  // Add the derivative to the sensitivity vector
  r1->addPointAdjResProduct(fdvSens, numDVs, 1.0, dx1);
  r2->addPointAdjResProduct(fdvSens, numDVs, 1.0, dx2);

  // Add the result from the base point
  dx0[0] = -(dx1[0] + dx2[0]);
  dx0[1] = -(dx1[1] + dx2[1]);
  dx0[2] = -(dx1[2] + dx2[2]);
  r0->addPointAdjResProduct(fdvSens, numDVs, 1.0, dx0);
}

/*
  Test the implementation of the derivative of the rotation matrix
  
  input:
  numDVs: the number of design variables
  dh:     the step size for the central-difference formula
*/
void TACSRefFrame::testRotation( int numDVs, double dh ){
  // Allocate design variable arrays
  TacsScalar *x = new TacsScalar[ numDVs ];
  TacsScalar *fdvSens = new TacsScalar[ numDVs ];

  // Set random values into the psi/phi vectors
  TacsScalar psi[3], phi[3];
  for ( int k = 0; k < 3; k++ ){
    psi[k] = -0.5 + (1.0*rand())/RAND_MAX;
    phi[k] = -0.5 + (1.0*rand())/RAND_MAX;
  }

  // Zero the sensitivity
  memset(fdvSens, 0, numDVs*sizeof(TacsScalar));

  // Get the variables and make sure that they are consistent
  getDesignVars(x, numDVs);
  setDesignVars(x, numDVs);

  // Compute the derivative and store it in fdvSens
  addRotationAdjResProduct(fdvSens, numDVs, psi, phi);

  // Set the vectors to zero
  TacsScalar dfdx[9], fd[9];
  memset(dfdx, 0, sizeof(dfdx));
  memset(fd, 0, sizeof(fd));

  // Set the vectors
  TACSGibbsVector *vec[3] = {r0, r1, r2};
  for ( int k = 0; k < 3; k++ ){
    // Retrieve the design variables
    const int *dvs;
    vec[k]->getVectorDesignVarNums(&dvs);

    for ( int i = 0; i < 3; i++ ){
      if (dvs[i] >= 0 && dvs[i] < numDVs){
	TacsScalar t[3];
	TacsScalar xtmp = x[dvs[i]];

#ifdef TACS_USE_COMPLEX
	// Evaluate the matrix at x + j*dh
	x[dvs[i]] = xtmp + TacsScalar(0.0, dh);
	setDesignVars(x, numDVs);
	matMult(C, phi, t);
	fd[3*k+i] = ImagPart(vecDot(psi, t))/dh;
#else
	// Evaluate C at (x + dh)
	x[dvs[i]] = xtmp + dh;
	setDesignVars(x, numDVs);
	matMult(C, phi, t);
	TacsScalar f1 = vecDot(psi, t);

	// Evaluate C at (x - dh)
	x[dvs[i]] = xtmp - dh;
	setDesignVars(x, numDVs);
	matMult(C, phi, t);
	TacsScalar f2 = vecDot(psi, t);
	
	fd[3*k+i] = 0.5*(f1 - f2)/dh;
#endif // TACS_USE_COMPLEX
	x[dvs[i]] = xtmp;

	// Track the values
	dfdx[3*k+i] = fdvSens[dvs[i]];
      }
    }
  }

  // Set the design variable values back to their original values
  setDesignVars(x, numDVs);

  writeErrorComponents(stdout, "TACSRefFrame", dfdx, fd, 9);

  // Free the allocated memory
  delete [] x;
  delete [] fdvSens;
}

/*
  The constructor for the rigid body dynamics 

  input:
  rot_param:  the rotation matrix parametrization
  mass:       the mass of the body
  c:          the first moment of inertia
  J:          the symmetric second moment of inertia
  g:          the acceleration due to gravity in the global frame
*/
TACSRigidBody::TACSRigidBody( const TacsScalar _mass, 
                              const TacsScalar _c[], 
                              const TacsScalar _J[], 
                              TACSRefFrame *_CRef,
                              TACSGibbsVector *_rInit, 
                              TACSGibbsVector *_vInit, 
                              TACSGibbsVector *_omegaInit, 
                              TACSGibbsVector *gvec ){
  // Set the initial values of the Euler parameters
  p[0] = 1.0;
  p[1] = p[2] = p[3] = 0.0;

  // Copy over the property values
  mass = _mass;

  // Copy over the first moment of inertia properties
  c[0] = _c[0];
  c[1] = _c[1];
  c[2] = _c[2];

  // Copy over the second moment of inertia properties
  J[0] = _J[0];
  J[1] = _J[1];
  J[2] = _J[2];
  J[3] = _J[3];
  J[4] = _J[4];
  J[5] = _J[5];

  // Copy the accel. due to gravity
  const TacsScalar *_g;
  gvec->getVector(&_g);
  g[0] = _g[0];
  g[1] = _g[1];
  g[2] = _g[2];

  // Copy over the initial reference frame
  CRef = _CRef;

  // Copy over the initial vectors of everything
  rInit = _rInit;
  vInit = _vInit;
  omegaInit = _omegaInit;

  // Increment the reference counts for things
  CRef->incref();
  rInit->incref();
  vInit->incref();
  omegaInit->incref();

  // Set the design variable numbers for the inertial properties
  massDV = -1;
  cDV[0] = cDV[1] = cDV[2] = -1;
  JDV[0] = JDV[1] = JDV[2] = 
    JDV[3] = JDV[4] = JDV[5] = -1;
}

/*
  Decrease the reference counts to everything
*/
TACSRigidBody::~TACSRigidBody(){
  CRef->decref();
  rInit->decref();
  vInit->decref();
  omegaInit->decref();
}

/*
  Set the design variable numbers associated with the inertial
  properties of the body

  input:
  massDV:   design variable number for the mass
  cDV:      design variable numbers for the first moment of mass
  JDV:      design variable numbers for the second moments of mass
*/
void TACSRigidBody::setDesignVarNums( int _massDV, 
                                    const int _cDV[], 
				    const int _JDV[] ){
  massDV = _massDV;
  if (_cDV){
    cDV[0] = _cDV[0];
    cDV[1] = _cDV[1];
    cDV[2] = _cDV[2];
  }
  if (_JDV){
    JDV[0] = _JDV[0];
    JDV[1] = _JDV[1];
    JDV[2] = _JDV[2];
    JDV[3] = _JDV[3];
    JDV[4] = _JDV[4];
    JDV[5] = _JDV[5];
  }
}

/*
  Set the design variable values
*/
void TACSRigidBody::setDesignVars( const TacsScalar dvs[], int numDVs ){
  // Set the mass design variable
  if (massDV >= 0 && massDV < numDVs){
    mass = dvs[massDV];
  }

  // Set the moment of mass variable
  for ( int k = 0; k < 3; k++ ){
    if (cDV[k] >= 0 && cDV[k] < numDVs){
      c[k] = dvs[cDV[k]];
    }
  }
  
  // Set the second moment of mass variables
  for ( int k = 0; k < 6; k++ ){
    if (JDV[k] >= 0 && JDV[k] < numDVs){
      J[k] = dvs[JDV[k]];
    }
  }

  // Set the reference frame design variables
  CRef->setDesignVars(dvs, numDVs);

  // Set the design variable values for the initial vectors
  rInit->setDesignVars(dvs, numDVs);
  vInit->setDesignVars(dvs, numDVs);
  omegaInit->setDesignVars(dvs, numDVs);
}

/*
  Retrieve the design variable values
*/
void TACSRigidBody::getDesignVars( TacsScalar dvs[], int numDVs ){
  // Get the mass design variable
  if (massDV >= 0 && massDV < numDVs){
    dvs[massDV] = mass;
  }

  // Get the moment of mass variables
  for ( int k = 0; k < 3; k++ ){
    if (cDV[k] >= 0 && cDV[k] < numDVs){
      dvs[cDV[k]] = c[k];
    }
  }
  
  // Get the second moment of mass variables
  for ( int k = 0; k < 6; k++ ){
    if (JDV[k] >= 0 && JDV[k] < numDVs){
      dvs[JDV[k]] = J[k];
    }
  }

  // Get the reference frame design variables
  CRef->getDesignVars(dvs, numDVs);

  // Get the design variable values for the initial vectors
  rInit->getDesignVars(dvs, numDVs);
  vInit->getDesignVars(dvs, numDVs);
  omegaInit->getDesignVars(dvs, numDVs);
}

/*
  Retrieve the design variable range
*/
void TACSRigidBody::getDesignVarRange( TacsScalar lb[], 
                                     TacsScalar ub[],
                                     int numDVs ){


}

/*
  Evaluate the rotation matrix (C), the angular velocity matrix (S),
  and their derivatives with respect to the rotational variables as a
  function of the input states.

  This code uses an Euler sequence. In general this should probably
  enable different options to be used.

  input:
  q:   the rigid-body states

  The function calculates and stores the following:
  C:   the rotation matrix
  S:   the angular rate matrix
  Cd:  the derivatives of C
  Sd:  the derivatives of S
*/
void TACSRigidBody::setVariables( const TacsScalar qkin[],
                                  const TacsScalar qdyn[] ){
  // Set the position of the rigid body
  r0[0] = qkin[0];
  r0[1] = qkin[1];
  r0[2] = qkin[2];

  // Extract the velocity
  v0[0] = qdyn[0];
  v0[1] = qdyn[1];
  v0[2] = qdyn[2];

  // Extract the angular velocity
  omega[0] = qdyn[3];
  omega[1] = qdyn[4];
  omega[2] = qdyn[5];

  // Retrieve the initial matrix
  const TacsScalar *C0;
  CRef->getRotation(&C0);

  if (rot_param == EULER_ANGLES){
    // Compute the sin/cos of the Euler angles
    const TacsScalar *theta = &qkin[3];
    TacsScalar c1 = cos(theta[0]);
    TacsScalar s1 = sin(theta[0]);
    TacsScalar c2 = cos(theta[1]);
    TacsScalar s2 = sin(theta[1]);
    TacsScalar c3 = cos(theta[2]);
    TacsScalar s3 = sin(theta[2]);
    
    // Compute the rotation matrix
    Cr[0] = c2*c3;
    Cr[1] = c2*s3;
    Cr[2] = -s2;
    
    Cr[3] = s1*s2*c3 - c1*s3;
    Cr[4] = s1*s2*s3 + c1*c3;
    Cr[5] = s1*c2;
    
    Cr[6] = c1*s2*c3 + s1*s3;
    Cr[7] = c1*s2*s3 - s1*c3;
    Cr[8] = c1*c2;
    matMatMult(Cr, C0, C);
        
    // Derivatives w.r.t. c1/s1
    TacsScalar D[9];
    D[0] = D[1] = D[2] = 0.0;

    D[3] = c1*s2*c3 + s1*s3;
    D[4] = c1*s2*s3 - s1*c3;
    D[5] = c1*c2;
    
    D[6] = -s1*s2*c3 + c1*s3;
    D[7] = -s1*s2*s3 - c1*c3;
    D[8] = -s1*c2;
    matMatMult(D, C0, &Cd[0]);
    
    // Derivatives w.r.t. c2/s2
    D[0] = -s2*c3;
    D[1] = -s2*s3;
    D[2] = -c2;
    
    D[3] = s1*c2*c3;
    D[4] = s1*c2*s3;
    D[5] = -s1*s2;
    
    D[6] = c1*c2*c3;
    D[7] = c1*c2*s3;
    D[8] = -c1*s2;
    matMatMult(D, C0, &Cd[9]);

    // Derivatives w.r.t. c3/s3
    D[0] = -c2*s3;
    D[1] = c2*c3;
    D[2] = 0.0;

    D[3] = -s1*s2*s3 - c1*c3;
    D[4] = s1*s2*c3 - c1*s3;
    D[5] = 0.0;

    D[6] = -c1*s2*s3 + s1*c3;
    D[7] = c1*s2*c3 + s1*s3;
    D[8] = 0.0;
    matMatMult(D, C0, &Cd[18]);

    // Set the derivatives to zero 
    memset(Sd, 0, 3*9*sizeof(TacsScalar));

    // Compute the angular rate matrix
    S[0] = 1.0;
    S[1] = 0.0;
    S[2] = -s2;
    
    S[3] = 0.0;
    S[4] = c1;
    S[5] = s1*c2;
    
    S[6] = 0.0;
    S[7] = -s1;
    S[8] = c1*c2;

    // Derivatives w.r.t. c1/s1
    TacsScalar *sd = Sd;
    sd[4] =-s1;
    sd[5] = c1*c2;
    sd[7] = -c1;
    sd[8] = -s1*c2;
    
    // Derivatives w.r.t c2/s2
    sd += 9;
    sd[2] = -c2;
    sd[5] = -s1*s2;
    sd[8] = -c1*s2;
  }
  else if (rot_param == EULER_PARAMETERS){
    // Compute the sin/cos of the Euler angles
    p[0] = qkin[3];
    p[1] = qkin[4];
    p[2] = qkin[5];
    p[3] = qkin[6];
  
    // Compute the rotation matrix as a function of the
    // Euler parameters
    Cr[0] = 1.0 - 2.0*(p[2]*p[2] + p[3]*p[3]);
    Cr[1] = 2.0*(p[1]*p[2] + p[3]*p[0]);
    Cr[2] = 2.0*(p[1]*p[3] - p[2]*p[0]);

    Cr[3] = 2.0*(p[2]*p[1] - p[3]*p[0]);
    Cr[4] = 1.0 - 2.0*(p[1]*p[1] + p[3]*p[3]);
    Cr[5] = 2.0*(p[2]*p[3] + p[1]*p[0]);

    Cr[6] = 2.0*(p[3]*p[1] + p[2]*p[0]);
    Cr[7] = 2.0*(p[3]*p[2] - p[1]*p[0]);
    Cr[8] = 1.0 - 2.0*(p[1]*p[1] + p[2]*p[2]);
    matMatMult(Cr, C0, C);
    
    // Compute the derivatives of the rotation matrix as 
    // a function of the Euler parameters
    TacsScalar D[9];
    D[0] = 0.0;
    D[1] = 2.0*p[3];
    D[2] =-2.0*p[2];

    D[3] =-2.0*p[3];
    D[4] = 0.0;
    D[5] = 2.0*p[1];

    D[6] = 2.0*p[2];
    D[7] =-2.0*p[1];
    D[8] = 0.0;
    matMatMult(D, C0, &Cd[0]);

    // Derivative w.r.t. p[1]
    D[0] = 0.0;
    D[1] = 2.0*p[2];
    D[2] = 2.0*p[3];

    D[3] = 2.0*p[2];
    D[4] =-4.0*p[1];
    D[5] = 2.0*p[0];

    D[6] = 2.0*p[3];
    D[7] =-2.0*p[0];
    D[8] =-4.0*p[1];
    matMatMult(D, C0, &Cd[9]);

    // Derivative w.r.t. p[2]
    D[0] =-4.0*p[2];
    D[1] = 2.0*p[1];
    D[2] =-2.0*p[0];

    D[3] = 2.0*p[1];
    D[4] = 0.0;
    D[5] = 2.0*p[3];

    D[6] = 2.0*p[0];
    D[7] = 2.0*p[3];
    D[8] =-4.0*p[2];
    matMatMult(D, C0, &Cd[18]);

    // Derivative w.r.t. p[3]
    D[0] =-4.0*p[3];
    D[1] = 2.0*p[0];
    D[2] = 2.0*p[1];

    D[3] =-2.0*p[0];
    D[4] =-4.0*p[3];
    D[5] = 2.0*p[2];

    D[6] = 2.0*p[1];
    D[7] = 2.0*p[2];
    D[8] = 0.0;
    matMatMult(D, C0, &Cd[27]);
  }
}

/*
  Retrieve the position, velocity and angular velocity vectors (all
  stored in the body-fixed frame).
  
  output:
  r0:    the position vector
  v0:    the velocity vector
  omega: the angular velocity vector
*/
void TACSRigidBody::getVariables( const TacsScalar *_r0[], 
				const TacsScalar *_v0[],
				const TacsScalar *_omega[] ){
  if (_r0){ *_r0 = r0; }
  if (_v0){ *_v0 = v0; }
  if (_omega){ *_omega = omega; }
}

/*
  Retrieve the initial values if the kinematic and dynamic variables

  output:
  qkin:  the kinematic variables
  qdyn:  the dynamic variable values
*/
void TACSRigidBody::getInitVariables( TacsScalar qkin[],
				    TacsScalar qdyn[] ){
  // Compute the initial points
  const TacsScalar *r;
  rInit->getVector(&r);
  qkin[0] = r[0];
  qkin[1] = r[1];
  qkin[2] = r[2];

  // Retrieve the initial coordinate transformation
  if (rot_param == EULER_ANGLES){
    qkin[3] = qkin[4] = qkin[5] = 0.0;
  }
  else {
    qkin[3] = 1.0;
    qkin[4] = qkin[5] = qkin[6] = 0.0;
  }

  // Get the initial velocity vector and transform it into the
  // local reference frame
  const TacsScalar *C0;
  CRef->getRotation(&C0);

  const TacsScalar *v;
  vInit->getVector(&v);
  matMult(C0, v, &qdyn[0]);

  // Get the initial angular velocity vector and transformm it
  // into the local reference frame
  const TacsScalar *omeg;
  omegaInit->getVector(&omeg);
  matMult(C0, omeg, &qdyn[3]);
}

/*
  Retrieve the initial variable vectors

  output:
  CRef:       the initial rotation matrix 
  rInit:      the initial position vector
  vInit:      the initial velocity
  omegaInit:  the initial angular velocity
*/
void TACSRigidBody::getInitVariables( TACSRefFrame **_CRef,
				    TACSGibbsVector **_rInit, 
				    TACSGibbsVector **_vInit,
				    TACSGibbsVector **_omegaInit ){
  if (_CRef){ *_CRef = CRef; }
  if (_rInit){ *_rInit = rInit; }
  if (_vInit){ *_vInit = vInit; }
  if (_omegaInit){ *_omegaInit = omegaInit; }
}

/*
  Retrieve the rotation matrix after it has been computed by the
  setRotation() code. No check is performed (or can be performed) to
  ensure that you've made the calls in the correct order, so be
  careful.

  output:
  C:     constant pointer to the rotation matrix
*/
void TACSRigidBody::getRotation( const TacsScalar *_C[] ){
  if (_C){ *_C = C; }
}

/*
  Retrieve the rotation matrix and its derivative after they 
  have been computed by the setRotation() code.

  output:
  C:     constant pointer to the rotation matrix
  C:     constant pointer to the derivatives of the rotation matrix
*/
void TACSRigidBody::getRotationDeriv( const TacsScalar *_C[],
				    const TacsScalar *_Cd[],
				    int *nparam ){
  if (_C){ *_C = C; }
  if (_Cd){ *_Cd = Cd; }
  *nparam = 3;
  if (rot_param == EULER_PARAMETERS){
    *nparam = 4;
  }
}

/*
  Add the derivative of the inner product of the rotation matrix with
  the two given input vectors. 

  The result of the computation can be written:

  fdvSens = d(psi^{T}*C*phi)/dx
  .       = d(psi^{T}*Cr*CRef*phi)/dx
  .       = d(t^{T}*CRef*phi)/dx
 
  This is required for computing the adjoint-residual vector product
  for the adjoint method. 

  input:
  numDVs:   the number of design variables/length of the fdvSens array
  psi:      the first vector
  phi:      the second vector

  output:
  fdvSens:  the array where the derivative is added 
*/
void TACSRigidBody::addRotationAdjResProduct( TacsScalar fdvSens[], 
					    int numDVs,
					    const TacsScalar psi[],
					    const TacsScalar phi[] ){
  TacsScalar t[3];
  matMultTrans(Cr, psi, t);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, t, phi);
}


/*
  Compute the kinematic and potential energies of the rigid body

  The kinetic energy is given by:
  
  T = 0.5*m*v0^{T}*v0 + 0.5*omega^{T}*J*omega
  .   + v0^{T}*C^{T}*omega^{x}*c

  where r0 and v0 are expressed in the global reference frame, and C
  is the rotation matrix from the global to body-fixed frames. The
  potential energy is due to the force of gravity

  U = -m*g^{T}*r0 - g^{T}*C^{T}*c

  input:
  time:   the simulation time
  X:      the reference node location
  vars:   the state vector
  dvars:  the time derivative of the states

  output:
  Te:   the kinematic energy
  Pe:   the potential energy
*/
void TACSRigidBody::computeEnergies( double time,
                                     TacsScalar *Te, 
                                     TacsScalar *Pe,
                                     const TacsScalar X[],
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[] ){  
  // Set the location
  const TaccScalar *r0 = &vars[0];
  const TacsScalar *v0 = &dvars[0]; 

  // Set the pointers to the Euler parameters
  TacsScalar eta = vars[3];
  const TacsScalar *eps = &vars[4];
  TacsScalar deta = dvars[3];
  const TacsScalar *deps = &dvars[4];

  // Compute the rotation matrix
  TacsScalar C[9];
  computeRotationMat(eta, eps, C);

  // Compute the angular velocity from the Euler parameters
  TacsScalar omega[3];
  computeSRateProduct(eta, eps, deta, deps, omega);

  // Add the contribution from the linear motion
  *Te = 0.5*mass*vecDot(v0, v0);
  
  // Add the contribution from the angular velocity
  TacsScalar tmp[3];
  matSymmMult(J, omega, tmp);
  *Te += 0.5*vecDot(omega, tmp);

  // Add the coupled contribution from the angular velocity/rotation
  crossProduct(1.0, c, omega, tmp);
  *Te -= vecDot(v0, tmp);

  // Compute the potential energy
  matMult(C, g, tmp);
  *Pe = -(mass*vecDot(r0, g) + vecDot(c, tmp));
}

/*
  Compute the residual of the governing equations for the rigid body
  motion.

  The equations of motion for the rigid body are derived using
  Lagrange's equations of motion with constraints. Within each body,
  the constraint enforcing the norm on the quaternion is imposed
  directly. The motion equations are written as follows:

  m*ddot{r} - c^{x}*dot{omega} - mass*g = 0
  S^{T}*dot{y} + 2*dot{S}^{T}*y - A^{T}*lambda = 0

  where y = J*omega + c^{x}*v

  Note that the matrix S relates the quaternion rates to the angular
  acceleration such that omega = S*dot{q}. The matrix S is given by:
  
  S(q) = 2[ -eps | (eta*I - eps^{x}) ]

  Note that the matrix S has the property that dot{S} = S(dot{q}).
  The transpose of S is given as follows:

  S^{T} = [      2*eps^{T}      ]
  .       [ 2*(eta*I + eps^{x}) ] 
  
  


*/
void TACSRigidBody::getResidual( double time, 
                                 TacsScalar res[],
                                 const TacsScalar X[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){
  // Zero the residual entries
  memset(res, 0, 8*sizeof(TacsScalar));

  // Set the location and its time derivatives
  const TaccScalar *r0 = &vars[0];
  const TacsScalar *v0 = &dvars[0]; 
  const Tacsscalar *a0 = &ddvars[0];

  // Set the pointers to the Euler parameters and all their time
  // derivatives
  TacsScalar eta = vars[3];
  const TacsScalar *eps = &vars[4];
  TacsScalar deta = dvars[3];
  const TacsScalar *deps = &dvars[4];
  TacsScalar ddeta = ddvars[3];
  const TacsScalar *ddeps = &ddvars[4];

  // Compute the angular velocity from the Euler parameters
  TacsScalar omega[3];
  computeSRateProduct(eta, eps, deta, deps, omega);

  // Compute the angular acceleration
  // domega = S(dot{q})*dot{q} + S(q)*ddot{q}
  TacsScalar domega[3];
  computeSRateProduct(deta, deps, deta, deps, domega);
  addSRateProduct(eta, eps, ddeta, ddeps, domega);

  // Compute the cross product
  TacsScalar cdw[3];
  crossProduct(1.0, c, domega, cdw);
  
  // The dynamics for the reference point
  r[0] = mass*(a0[0] - g[0]) - cdw[0];
  r[1] = mass*(a0[1] - g[1]) - cdw[1];
  r[2] = mass*(a0[2] - g[2]) - cdw[2];
  
  // Compute y = J*omega + c^{x}*v0
  TacsScalar y[3];
  matSymmMult(J, omega, y);
  crossProductAdd(1.0, c, v0, y);

  // Compute dy = J*dot{omega} + c^{x}*dot{v0}
  TacsScalar dy[3];
  matSymmMult(J, domega, dy);
  crossProductAdd(1.0, c, a0, dy);

  // The dynamics for the rotation
  // Add S^{T}*dy
  addSRateTransProduct(1.0, eta, eps, dy, &r[3], &r[4]);

  // Add 2*dot{S}^{T}*y
  addSRateTransProduct(2.0, deta, deps, y, &r[3], &r[4]);


}

/*
  Compute the Jacobian of the governing equations
*/
void TACSRigidBody::getJacobian( double time, TacsScalar J[],
                                 double alpha, double beta, double gamma,
                                 const TacsScalar X[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){


}


/*
  Test the rotation matrix and S matrix to see if they are implemented
  correctly and that their derivatives match.  Use the
  finite-difference or complex-step to check the implementation
  depending on whether the code is compiled in complex mode.

  input:
  q:   temporary state variables
  dh:  finite-difference (complex-step) step size
*/
void TACSRigidBody::testRotation( TacsScalar qkin[], double dh ){
  // Evaluate the derivative
  TacsScalar qdyn[6] = {0.0, 0.0, 0.0, 
			0.0, 0.0, 0.0};
  setVariables(qkin, qdyn);

  // Pointer to the rotation DOF
  TacsScalar *param = &qkin[3];

  // Temporary storage since the setRotation() code over-writes
  // previously computed values
  TacsScalar Ct[9], Cdt[4*9];

  // Copy over the matrices that we're testing
  memcpy(Ct, C, sizeof(Ct));
  memcpy(Cdt, Cd, sizeof(Cdt));

  int nparam = 3;
  if (rot_param == EULER_PARAMETERS){
    nparam = 4;
  }

  // Test the derivatives using finite-difference or complex-step
  for ( int i = 0; i < nparam; i++ ){
    TacsScalar qtmp = param[i];

#ifdef TACS_USE_COMPLEX
    param[i] = qtmp + TacsScalar(0.0, dh);
    setVariables(qkin, qdyn);

    // Compute the derivative using complex step
    for ( int k = 0; k < 9; k++ ){
      C[k] = ImagPart(C[k])/dh;
    }
#else 
    param[i] = qtmp + dh;
    setVariables(qkin, qdyn);

    // Compute the derivative using central difference
    for ( int k = 0; k < 9; k++ ){
      C[k] = (C[k] - Ct[k])/dh;
    }
#endif // TACS_USE_COMPLEX
    // Reset the DOF
    param[i] = qtmp;

    char descript[64];
    sprintf(descript, "dC/dtheta%d", i+1);
    writeErrorComponents(stdout, descript, &Cdt[9*i], C, 9);
  }
}

/*
  Retrieve the inertial properties from the body

  Note that you can input a NULL pointer for unwanted inertial
  properties. This code checks for NULL pointers before assignment.

  input/output:
  mass:   the mass of the body
  c:      the first moment of mass
  J:      the moment of inertias of the body
*/
void TACSRigidBody::getInertial( TacsScalar *_mass, 
			       const TacsScalar *_c[], 
			       const TacsScalar *_J[] ){
  if (_mass){ *_mass = mass; }
  if (_c){ *_c = c; }
  if (_J){ *_J = J; }
}

/*
  Add the derivative of the inertial properties to the given
  derivative array. 
  
  Note that you can supply a NULL pointer for the dcdx and dJdx
  derivatives. 

  input:
  numDVs:    the length of the fdvSens array
  dmdx:      the derivative w.r.t. the mass
  dcdx:      the derivative w.r.t. the first moment of mass
  dJdx:      the derivative w.r.t. the moment of inertia

  output:
  fdvSens:   the derivative array
*/
void TACSRigidBody::addInertialAdjResProduct( TacsScalar fdvSens[],
					    int numDVs, 
					    TacsScalar dmdx, 
					    const TacsScalar dcdx[],
					    const TacsScalar dJdx[] ){
  // Add the derivative from the mass
  if (massDV >= 0 && massDV < numDVs){ 
    fdvSens[massDV] += dmdx; 
  }
  
  // Add the derivative from the first moment of mass
  if (dcdx){
    for ( int k = 0; k < 3; k++ ){
      if (cDV[k] >= 0 && cDV[k] < numDVs){
	fdvSens[cDV[k]] += dcdx[k];
      }
    }
  }
  
  // Add the result from the inertial terms
  if (dJdx){
    for ( int k = 0; k < 6; k++ ){
      if (JDV[k] >= 0 && JDV[k] < numDVs){
	fdvSens[JDV[k]] += dJdx[k];
      }
    }
  }
}

/*
  Compute the kinematic and potential energies of the rigid body

  input:
  q:    the dynamics state vector

  output:
  Te:   the kinematic energy
  Pe:   the potential energy
*/
void TACSRigidBody::getEnergies( TacsScalar *Te,
			       TacsScalar *Pe ){
  TacsScalar tmp[3]; // temporary vector

  // Add the contribution from the linear motion
  *Te = 0.5*mass*vecDot(v0, v0);
  
  // Add the contribution from the angular velocity
  matSymmMult(J, omega, tmp);
  *Te += 0.5*vecDot(omega, tmp);

  // Add the coupled contribution from the angular velocity/rotation
  crossProduct(1.0, c, omega, tmp);
  *Te -= vecDot(v0, tmp);

  // Compute the potential energy
  matMult(C, g, tmp);
  *Pe = -(mass*vecDot(r0, g) + vecDot(c, tmp));
}

/*
  Compute the residual of the governing equations of dynamics using
  the input DOF in q and their derivatives in qdot.

  input:
  q:     the dynamic state variables
  qdot:  the time-derivatives of the dynamic states
  
  output:
  res:   the residual
*/
void TACSRigidBody::getResidual( TacsScalar rkin[], 
			       TacsScalar rdyn[],
			       const TacsScalar qkinDot[], 
			       const TacsScalar qdynDot[] ){
  TacsScalar tmp[3]; // temporary vector

  // Set pointers to the time-derivatives of the DOF 
  const TacsScalar *rdot = &qkinDot[0];
  const TacsScalar *vdot = &qdynDot[0];
  const TacsScalar *omegadot = &qdynDot[3];

  // First kinematic equation: C*dot{r} - v = 0 
  matMult(C, rdot, &rkin[0]);
  vecAxpy(-1.0, v0, &rkin[0]);

  if (rot_param == EULER_ANGLES){
    // Second kinematic equation: S*dot{theta} - omega = 0
    const TacsScalar *thetadot = &qkinDot[3];
    matMult(S, thetadot, &rkin[3]);
    vecAxpy(-1.0, omega, &rkin[3]);
  }
  else if (rot_param == EULER_PARAMETERS){
    // Second kinematic equation:
    const TacsScalar *pdot = &qkinDot[3];
    // omega - 2*p[0]*dot{p[1:]} + 2.0*p[1:]^{x}*dot{p[1:]} 
    //   + dot{p[0]}*p[1:]) = 0 
    crossProduct(2.0, &p[1], &pdot[1], &rkin[3]);
    vecAxpy(1.0, omega, &rkin[3]);
    vecAxpy(-2.0*p[0], &pdot[1], &rkin[3]);
    vecAxpy(2.0*pdot[0], &p[1], &rkin[3]);

    rkin[6] = (p[0]*p[0] + p[1]*p[1] + 
	       p[2]*p[2] + p[3]*p[3] - 1.0);
  }

  // Set the residual to zero
  memset(rdyn, 0, ndyn*sizeof(TacsScalar));
  
  // Add the force balance equation:
  // -------------------------------
  // m*dot{v} - c^x*dot{omega} + 
  //   omega^{x}(m*v_{0} - c^{x}*{omega}) - m*Cbi*g = 0
  vecAxpy(mass, vdot, &rdyn[0]);
  crossProductAdd(-1.0, c, omegadot, &rdyn[0]);
  crossProductAdd(mass, omega, v0, &rdyn[0]);

  // Add: -omega^{x}*c^{x}*omega
  crossProduct(1.0, c, omega, tmp);
  crossProductAdd(-1.0, omega, tmp, &rdyn[0]);

  // Add -m*Cbi*g
  matMult(C, g, tmp);
  vecAxpy(-mass, tmp, &rdyn[0]);

  // Add the moment balance equation:
  // --------------------------------
  // c^{x}*dot{v} + J*dot{omega} + 
  //   c^{x}*omega^{x}*v + omega^{x}*J*omega - c^{x}*Cbi*g = 0
  crossProductAdd(1.0, c, vdot, &rdyn[3]);
  matSymmMultAdd(J, omegadot, &rdyn[3]);

  // Add c^{x}*omega^{x}*v
  crossProduct(1.0, omega, v0, tmp);
  crossProductAdd(1.0, c, tmp, &rdyn[3]);
  
  // Add omega^{x}*J*omega
  matSymmMult(J, omega, tmp);
  crossProductAdd(1.0, omega, tmp, &rdyn[3]);

  // Add - c^{x}*Cbi*g
  matMult(C, g, tmp);
  crossProductAdd(-1.0, c, tmp, &rdyn[3]);
}

/*
  Compute the Jacobian of the governing equations of motion
  with respect to a linear combination of q and qdot.

  input:
  alpha:   the coefficient for the state variables
  beta:    the coefficient for the time-derivative terms
  dkinDot: the time-derivatives of the kinematic states
  dkinDot: the time-derivatives of the dynamic states

  output:
  D11:   the kinematic-kinematic Jacobian matrix
  D12:   the kinematic-dynamic Jacobian matrix
  D21:   the dynamic-kinematic Jacobian matrix
  D22:   the dynamic-dynamic Jacobian matrix
*/
void TACSRigidBody::getJacobian( TacsScalar D11[],
			       TacsScalar D12[],
			       TacsScalar D21[],
			       TacsScalar D22[],
			       double alpha,
			       double beta,
			       const TacsScalar qkinDot[], 
			       const TacsScalar qdynDot[] ){
  TacsScalar tmp[12], tmp2[12]; // temporary matrix

  // Set pointers to the time-derivatives of the DOF 
  const TacsScalar *rdot = &qkinDot[0];
  
  // Zero the components of the Jacobian matrix
  memset(D11, 0, nkin*nkin*sizeof(TacsScalar));
  memset(D12, 0, nkin*ndyn*sizeof(TacsScalar));
  memset(D21, 0, ndyn*nkin*sizeof(TacsScalar));
  memset(D22, 0, ndyn*ndyn*sizeof(TacsScalar));

  int nparam = 3;
  if (rot_param == EULER_PARAMETERS){
    nparam = 4;
  }

  // Derivative of C*dot{r} - v = 0
  // ------------------------------
  // Add beta*C to the diagonal block
  addBlockMat(beta, C, nkin, &D11[0]);
  
  for ( int k = 0; k < nparam; k++ ){
    // Add d(C)/d(param)*dot{r} to the D11 matrix
    matMult(&Cd[9*k], rdot, tmp);
    addVecMat(alpha, tmp, nkin, &D11[3+k]);
 
    // Add the derivative of mass*C*g
    matMult(&Cd[9*k], g, tmp);
    addVecMat(-mass*alpha, tmp, nkin, &D21[3+k]);

    // Derivative w.r.t. theta
    // d(c^{x}*C*g)/d(theta)
    crossProduct(1.0, c, tmp, tmp2);
    int offset = 3*nkin;
    addVecMat(-alpha, tmp2, nkin, &D21[offset+3+k]);
  }

  // Add -I from the v and omega terms
  addBlockIdent(-alpha, ndyn, &D12[0]);

  if (rot_param == EULER_ANGLES){
    const TacsScalar *thetadot = &qkinDot[3];

    // Derivative of S*dot{theta} - omega = 0
    // --------------------------------------
    int offset = 3*nkin;
    addBlockMat(beta, S, nkin, &D11[offset+3]);
    
    // Add d(S)/dtheta*dot{theta} to the block matrix
    matMult(&Sd[0], thetadot, &tmp[0]);
    matMult(&Sd[9], thetadot, &tmp[3]);
    matMult(&Sd[18], thetadot, &tmp[6]);
    addBlockMatTrans(alpha, tmp, nkin, &D11[offset+3]);

    // Add -I from the omega term
    offset = 3*ndyn;
    addBlockIdent(-alpha, nkin, &D12[offset+3]);
  }
  else if (rot_param == EULER_PARAMETERS){
    // Take the derivative w.r.t. p/dot{p}
    // dot{p[1:]} - 0.5*(p[1:]^{x}*omega + p[0]*omega) = 0
    int offset = 3*nkin;
    const TacsScalar *pdot = &qkinDot[3];
    
    addBlockSkew(2.0*beta, &p[1], nkin, &D11[offset+4]);
    addBlockSkew(-2.0*alpha, &pdot[1], nkin, &D11[offset+4]);

    addVecMat(-2.0*alpha, &pdot[1], nkin, &D11[offset+3]);
    addVecMat(2.0*beta, &p[1], nkin, &D11[offset+3]);

    addBlockIdent(-2.0*beta*p[0], nkin, &D11[offset+4]);
    addBlockIdent(2.0*pdot[0]*alpha, nkin, &D11[offset+4]);

    // Take the derivative of the constraint
    // p^{T}*p - 1.0 = 0
    for ( int k = 0; k < 4; k++ ){
      D11[6*nkin + 3+k] = 2.0*p[k]*alpha;
    }
    
    // Take the derivative w.r.t. omega
    offset = 3*ndyn;
    addBlockIdent(alpha, ndyn, &D12[offset+3]);
  }

  // Derivative of the system dynamics
  // ---------------------------------

  // Derivative w.r.t. v
  addBlockIdent(mass*beta, ndyn, &D22[0]);
  addBlockSkew(mass*alpha, omega, ndyn, &D22[0]); 
  
  // Derivative w.r.t. omega
  addBlockSkew(-beta, c, ndyn, &D22[3]); 
  addBlockSkew(-mass*alpha, v0, ndyn, &D22[3]); 
  addBlockSkewSkew(-alpha, omega, c, ndyn, &D22[3]);
  crossProduct(1.0, c, omega, tmp);
  addBlockSkew(alpha, tmp, ndyn, &D22[3]);
  
  // Derivative of the moment dynamics
  // ---------------------------------

  // Derivative w.r.t. v
  int offset = 3*ndyn;
  addBlockSkew(beta, c, ndyn, &D22[offset]);
  addBlockSkewSkew(alpha, c, omega, ndyn, &D22[offset]); 
  
  // Derivative w.r.t. omega
  addBlockSymmMat(beta, J, ndyn, &D22[offset+3]);

  // Add (J*omega)^{x}
  matSymmMult(J, omega, tmp);
  addBlockSkew(-alpha, tmp, ndyn, &D22[offset+3]);
  addBlockSkewSkew(-alpha, c, v0, ndyn, &D22[offset+3]);
  
  // Add omega^{x}*J to the block matrix
  tmp[0] = J[2]*omega[1] - J[1]*omega[2];
  tmp[1] = J[4]*omega[1] - J[3]*omega[2];
  tmp[2] = J[5]*omega[1] - J[4]*omega[2];

  tmp[3] = J[0]*omega[2] - J[2]*omega[0];
  tmp[4] = J[1]*omega[2] - J[4]*omega[0];
  tmp[5] = J[2]*omega[2] - J[5]*omega[0];

  tmp[6] = J[1]*omega[0] - J[0]*omega[1];
  tmp[7] = J[3]*omega[0] - J[1]*omega[1];
  tmp[8] = J[4]*omega[0] - J[2]*omega[1];
  addBlockMat(alpha, tmp, ndyn, &D22[offset+3]);
}

/*
  Add the derivative of the product of the adjoint vector with the
  residual vector to the given design variable vector.

  Note that when the design variables are just inertial properties,
  the derivative of the residual only involves the kinematics.
 
  input:
  numDVs:   the number of design variables/length of the fdvSens array
  kinAdj:   the kinematic adjoint variables
  dynAdj:   the dynamic adjoint variables

  output:
  fdvSens:  the adjoint-residual product
*/
void TACSRigidBody::addAdjResProduct( TacsScalar fdvSens[], 
				    int numDVs,
				    const TacsScalar qkinDot[],
				    const TacsScalar qdynDot[],
				    const TacsScalar kinAdj[],
				    const TacsScalar dynAdj[] ){
  TacsScalar tmp[3]; // temporary vector

  // Set pointers to the time-derivatives of the DOF 
  const TacsScalar *rdot = &qkinDot[0];
  const TacsScalar *vdot = &qdynDot[0];
  const TacsScalar *omegadot = &qdynDot[3];

  // Add the derivative: d(kinAdj^{T}*(Cr*CRef*dot{r} - v0))/dx
  matMultTrans(Cr, &kinAdj[0], tmp);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, tmp, rdot);

  // Add the derivative: -mass*d(dynAd^{T}*Cr*CRef*g)/dx
  matMultTrans(Cr, &dynAdj[0], tmp);
  vecScale(-mass, tmp);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, tmp, g);

  // Add the derivative: -d(dynAdj^{T}*c^{x}*Cr*CRef*g)/dx
  TacsScalar s[3];
  crossProduct(1.0, c, &dynAdj[3], s);
  matMultTrans(Cr, s, tmp);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, tmp, g);

  // The components of the derivative that are computed
  TacsScalar dmdx, dcdx[3], dJdx[6];

  // Add the terms to the derivative of the mass
  dmdx = vecDot(vdot, &dynAdj[0]);

  crossProduct(1.0, omega, v0, tmp);
  dmdx += vecDot(tmp, &dynAdj[0]);

  matMult(C, g, tmp);
  dmdx -= vecDot(tmp, &dynAdj[0]);

  // Add the terms from the derivative of the first moment of mass
  crossProduct(-1.0, omegadot, &dynAdj[0], dcdx);
  
  crossProduct(1.0, omega, &dynAdj[0], tmp);
  crossProductAdd(1.0, omega, tmp, dcdx);

  crossProductAdd(1.0, vdot, &dynAdj[3], dcdx);
  
  crossProduct(1.0, omega, v0, tmp);
  crossProductAdd(-1.0, &dynAdj[3], tmp, dcdx);

  matMult(C, g, tmp);
  crossProductAdd(1.0, &dynAdj[3], tmp, dcdx);

  // Add the contributions to the derivative of the second
  // moments/products of inertia
  crossProduct(1.0, omega, &dynAdj[3], tmp);

  const TacsScalar *adj = &dynAdj[3];
  // First compute the contributions from the diagonal parts
  // of the inertial matrix
  dJdx[0] = adj[0]*omegadot[0] - tmp[0]*omega[0];
  dJdx[3] = adj[1]*omegadot[1] - tmp[1]*omega[1];
  dJdx[5] = adj[2]*omegadot[2] - tmp[2]*omega[2];

  // Compute the contributions from the off-diagonals in J
  dJdx[1] = adj[0]*omegadot[1] + adj[1]*omegadot[0] 
    - (tmp[0]*omega[1] + tmp[1]*omega[0]);
  dJdx[2] = adj[0]*omegadot[2] + adj[2]*omegadot[0] 
    - (tmp[0]*omega[2] + tmp[2]*omega[0]);
  dJdx[4] = adj[1]*omegadot[2] + adj[2]*omegadot[1] 
    - (tmp[1]*omega[2] + tmp[2]*omega[1]);

  // Add the derivative from the mass
  if (massDV >= 0 && massDV < numDVs){
    fdvSens[massDV] += dmdx;
  }

  // Add the derivative from the first moment of mass
  for ( int k = 0; k < 3; k++ ){
    if (cDV[k] >= 0 && cDV[k] < numDVs){
      fdvSens[cDV[k]] += dcdx[k];
    }
  }
  
  // Add the result from the inertial terms
  for ( int k = 0; k < 6; k++ ){
    if (JDV[k] >= 0 && JDV[k] < numDVs){
      fdvSens[JDV[k]] += dJdx[k];
    }
  }
}

/*
  Add the derivative of the product of the given adjoint vectors with
  the initial design variable values to the given array.

  input:
  numDVs:  the number of design variables
  kinAdj:  adjoint variables associated with the kinematic states
  dynAdj:  adjoint variables associated with the dynamic states

  output:
  fdvSens: derivative array containing the result
*/
void TACSRigidBody::addAdjInitVariableProduct( TacsScalar fdvSens[], 
					     int numDVs,
					     const TacsScalar kinAdj[], 
					     const TacsScalar dynAdj[] ){
  // Add the contributions from the initial reference point
  rInit->addPointAdjResProduct(fdvSens, numDVs, 1.0, &kinAdj[0]);

  // Add the contributions that depend on an initial rotation from the
  // global reference frame to the local body-fixed frame
  TacsScalar t[3];
  const TacsScalar *C0;
  CRef->getRotation(&C0);

  // Add the terms from the initial velocity
  const TacsScalar *v;
  vInit->getVector(&v);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, &dynAdj[0], v);
  
  matMultTrans(C0, &dynAdj[0], t);
  vInit->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);

  // Add the terms from the initial angular velocity
  const TacsScalar *omeg;
  omegaInit->getVector(&omeg);
  CRef->addRotationAdjResProduct(fdvSens, numDVs, &dynAdj[3], omeg);
  
  matMultTrans(C0, &dynAdj[3], t);
  omegaInit->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);
}


/*
  The following class implements a spherical joint between two rigid
  bodies.

  The spherical joint does to permit the transmission of a torque
  between the bodies. The boides are fixed at points that are
  specified within each of the body-fixed frames.

  input:
  bodyA:  the pointer to the first body
  bodyB:  the pointer to the second body
  point:  the attachment point in the inertial frame
*/
TACSDynSphericalJoint::TACSDynSphericalJoint( TACSRigidBody *_bodyA,
					      TACSRigidBody *_bodyB,
					      TACSGibbsVector *_point ){
  // Copy over the pointer to the first body object
  bodyA = _bodyA;
  bodyA->incref();

  // Copy the pointer to the second body
  bodyB = _bodyB;
  bodyB->incref();

  // Increase the reference count to the point object
  point = _point;
  point->incref();

  updatePoints();
}

/*
  The following class implements a spherical joint between a fixed
  point and a given body.

  The spherical joint does to permit the transmission of a torque
  between the body and the fixed point. 

  input:
  bodyA:  the pointer to the first body
  point:  the inertial reference point
*/
TACSDynSphericalJoint::TACSDynSphericalJoint( TACSRigidBody *_bodyA,
					      TACSGibbsVector *_point ){
  // Copy over the pointer to the first body object
  bodyA = _bodyA;
  bodyA->incref();

  bodyB = NULL;
  xB[0] = xB[1] = xB[2] = 0.0;

  point = _point;
  point->incref();

  updatePoints();
}

/*
  Free/deallocate the joint
*/
TACSDynSphericalJoint::~TACSDynSphericalJoint(){
  bodyA->decref();
  if (bodyB){ 
    bodyB->decref();
  }
  point->decref();
}

/*
  Update the data used internally to compute the constraints.

  This takes the initial reference points for the bodies A and B,
  given as rA0 and rB0, and the initial point for the connection
  between the A and B bodies, x0, and computes the vectors xA and xB
  in the local body-fixed frame. This makes it easier to specify the
  initial locations of all the points in the inertial (global)
  reference frame.
*/
void TACSDynSphericalJoint::updatePoints(){
  // The body-fixed frames and vectors
  TACSRefFrame *FA0, *FB0;
  TACSGibbsVector *rA0vec, *rB0vec;
  const TacsScalar *CA0, *CB0;
  const TacsScalar *rA0, *rB0, *x0;

  // Retrieve the initial frames/position vectors
  bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);

  // Retrieve the rotation matricies
  FA0->getRotation(&CA0);

  // Retrieve the initial point locations
  rA0vec->getVector(&rA0);
  point->getVector(&x0);

  // Set the position of the constrained point within the first
  // body-fixed frame
  TacsScalar t[3];
  t[0] = x0[0] - rA0[0];
  t[1] = x0[1] - rA0[1];
  t[2] = x0[2] - rA0[2];
  matMult(CA0, t, xA);

  if (bodyB){
    bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);
    FB0->getRotation(&CB0);
    rB0vec->getVector(&rB0);

    // Set the position of the constrained point within the second
    // body-fixed frame
    t[0] = x0[0] - rB0[0];
    t[1] = x0[1] - rB0[1];
    t[2] = x0[2] - rB0[2];
    matMult(CB0, t, xB);
  }
  else {
    // Set the fixed point
    xI[0] = x0[0];
    xI[1] = x0[1];
    xI[2] = x0[2];
  }
}

/*
  Set the design variable values
*/
void TACSDynSphericalJoint::setDesignVars( const TacsScalar dvs[], 
					   int numDVs ){
  bodyA->setDesignVars(dvs, numDVs);
  if (bodyB){
    bodyB->setDesignVars(dvs, numDVs);
  }
  point->setDesignVars(dvs, numDVs);
  updatePoints();
}

/*
  Get the design variable values
*/
void TACSDynSphericalJoint::getDesignVars( TacsScalar dvs[], int numDVs ){
  bodyA->getDesignVars(dvs, numDVs);
  if (bodyB){
    bodyB->getDesignVars(dvs, numDVs);
  }
  point->getDesignVars(dvs, numDVs);
}

/*
  Retrieve the bodies associated with this kinematic constraint
*/
void TACSDynSphericalJoint::getBodies( TACSRigidBody **_bodyA,
				       TACSRigidBody **_bodyB ){
  *_bodyA = bodyA;
  *_bodyB = bodyB;
}

/*
  Compute the residual of the constraint equation.

  The constraint equation is expressed in the global coordinate frame
  as follows:

  CA^{T}*xA + rA - CB^{T}*xB - rB = 0
  g = 0

  Note that the use of the transpose converts the components to the
  global inertial reference frame.

  input:
  q:    the kinematic/dynamic states
  
  output:
  res:  the residual
*/
void TACSDynSphericalJoint::getResidual( TacsScalar res[], 
					 const TacsScalar fr[] ){
  // Temporary vector
  TacsScalar tmp[3];

  // Retrieve the position vector and rotation matrix
  const TacsScalar *rA, *CA;
  bodyA->getVariables(&rA, NULL, NULL);
  bodyA->getRotation(&CA);

  // Form the residual of the first constraint equation
  // CA^{T}*xA + rA - CB^{T}*xB - rB - xI = 0
  matMultTrans(CA, xA, &res[0]);
  vecAxpy(1.0, rA, &res[0]);

  if (bodyB){
    // Retrieve the position vector and rotation matrix
    // for the second body
    const TacsScalar *rB, *CB;
    bodyB->getVariables(&rB, NULL, NULL);
    bodyB->getRotation(&CB);

    matMultTrans(CB, xB, tmp);
    vecAxpy(1.0, rB, tmp);
    vecAxpy(-1.0, tmp, &res[0]);
  }
  else {
    // Add in the fixed position
    vecAxpy(-1.0, xI, &res[0]);
  }

  // Form the constraint for the torques
  res[3] = fr[3];
  res[4] = fr[4];
  res[5] = fr[5];
}

/*
  Add the derivative of the product of the adjoint variables and
  the kinematic constraint residuals to the given input array.

  input:
  numDVs:  the number of design variables
  resAdj:  the adjoint variables
  fr:      the reaction forces/torques

  output:
  fdvSens: the array to which the result is added
*/
void TACSDynSphericalJoint::addAdjResProduct( TacsScalar fdvSens[], 
					      int numDVs,
					      const TacsScalar resAdj[],
					      const TacsScalar fr[] ){
  // The body-fixed frames and vectors
  TACSRefFrame *FA0, *FB0;
  TACSGibbsVector *rA0vec, *rB0vec;
  const TacsScalar *CA0, *CB0, *CA, *CB;
  const TacsScalar *rA0, *rB0, *x0;

  // Retrieve the initial frames/position vectors
  bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);

  // Retrieve the rotation matricies
  bodyA->getRotation(&CA);
  FA0->getRotation(&CA0);

  // Retrieve the initial point locations
  rA0vec->getVector(&rA0);
  point->getVector(&x0);

  // Add the derivative of from the points
  TacsScalar t[3], s[3];
  matMult(CA, &resAdj[0], s); // s = CA*resAdj
  matMultTrans(CA0, s, t); // t = CA0^{T}*CA*resAdj
  point->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);
  rA0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, t);

  // Add the derivative from the initial reference frame position
  t[0] = x0[0] - rA0[0];
  t[1] = x0[1] - rA0[1];
  t[2] = x0[2] - rA0[2];
  FA0->addRotationAdjResProduct(fdvSens, numDVs, s, t);

  // Add the derivative of: resAdj^{T}*CA^{T}*xA = xA^{T}*CA*resAdj
  bodyA->addRotationAdjResProduct(fdvSens, numDVs, xA, &resAdj[0]);

  // If it exists, add the derivative of the rotation matrix from
  // bodyB using the same approach as above
  if (bodyB){
    bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);

    // Retrieve the rotation matricies
    bodyB->getRotation(&CB);
    FB0->getRotation(&CB0);

    // Retrieve the initial point locations
    rB0vec->getVector(&rB0);

    // Add the derivative of from the points
    matMult(CB, &resAdj[0], s); // s = CB*resAdj
    matMultTrans(CB0, s, t); // t = CB0^{T}*CB*resAdj
    point->addPointAdjResProduct(fdvSens, numDVs, -1.0, t);
    rB0vec->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);

    // Add the derivative from the initial reference frame position
    t[0] = -(x0[0] - rB0[0]);
    t[1] = -(x0[1] - rB0[1]);
    t[2] = -(x0[2] - rB0[2]);
    FB0->addRotationAdjResProduct(fdvSens, numDVs, s, t);

    // Add the derivative of: -resAdj^{T}*CB^{T}*xB = xB^{T}*CB*resAdj
    TacsScalar xBneg[3];
    xBneg[0] = -xB[0];
    xBneg[1] = -xB[1];
    xBneg[2] = -xB[2];
    bodyB->addRotationAdjResProduct(fdvSens, numDVs, xBneg, &resAdj[0]);
  }
  else {
    point->addPointAdjResProduct(fdvSens, numDVs, -1.0, &resAdj[0]);
  }
}

/*
  Compute the Jacobian of the constraint equation w.r.t. the dynamic
  variables (e.g. forces/moments)

  input:
  alpha:  the coefficient for the variables 
  fr:     the dynamic variables

  output:
  D:  the Jacobian w.r.t. the dynamic variables
*/
void TACSDynSphericalJoint::getJacobian( TacsScalar D[], 
					 double alpha,
					 const TacsScalar fr[] ){
  // Set this part of the Jacobian to zero
  memset(D, 0, ncon*ncon*sizeof(TacsScalar));

  // Set the entries corresponding to the torques to the identity
  // matrix
  for ( int k = 3; k < ncon; k++ ){ 
    D[(ncon+1)*k] = alpha;
  }
}

/*
  Add terms to the residual of the governing equations to account for
  the reaction forces/torques internal to the multibody dynamics
  system.

  input:
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  input/output:
  rdyn:   the dynamic residuals
*/
void TACSDynSphericalJoint::addBodyResidual( TacsScalar rdyn[],
					     TACSRigidBody *rbody,
					     const TacsScalar fr[] ){
  if (rbody == bodyA){
    // Get the rotation matrix
    const TacsScalar *CA;
    bodyA->getRotation(&CA);

    // Add terms to the body residual
    TacsScalar fA[3];

    // Transform the constraint forces to the local frame
    matMult(CA, &fr[0], fA);
    vecAxpy(-1.0, fA, &rdyn[0]);
    
    // Add the constraint 
    crossProductAdd(-1.0, xA, fA, &rdyn[3]);
  }
  else if (bodyB && rbody == bodyB){
    // Get the rotation matrix
    const TacsScalar *CB;
    bodyB->getRotation(&CB);

    // Add terms to the body residual
    TacsScalar fB[3];

    // Transform the constraint forces to the local frame
    matMult(CB, &fr[0], fB);
    vecAxpy(1.0, fB, &rdyn[0]);
    
    // Add the constraint 
    crossProductAdd(1.0, xB, fB, &rdyn[3]);
  }
}

/*
  Add components to the Jacobian of the dynamics with respect to the
  kinematic variables of the rigid body.

  input:
  alpha:  the coefficient for the variables 
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  input/output:
  D:      the Jacobian of the dynamics w.r.t. the kinematic variables
*/
void TACSDynSphericalJoint::addBodyJacobian( TacsScalar D[],
					     double alpha,
					     TACSRigidBody *rbody,
					     const TacsScalar fr[] ){
  // Add the contributions to the Jacobian
  if (rbody == bodyA){
    // Get the number of kinematic terms
    int nkin = 0;
    bodyA->getNumDof(&nkin, NULL);

    // Get the rotation matrix
    int nparam = 0;
    const TacsScalar *CA, *CAd;
    bodyA->getRotationDeriv(&CA, &CAd, &nparam);

    // Add terms to the body residual
    TacsScalar fA[3], xAxI1[3];

    // Transform the constraint forces to the local frame
    for ( int k = 0; k < nparam; k++ ){
      matMult(&CAd[9*k], &fr[0], fA);
      addVecMat(-alpha, fA, nkin, &D[3+k]); 

      // Add the terms from the cross-product
      crossProduct(1.0, xA, fA, xAxI1);
      addVecMat(-alpha, xAxI1, nkin, &D[3*nkin+3+k]);
    }
  }
  else if (bodyB && rbody == bodyB){
   // Get the number of kinematic terms
    int nkin = 0;
    bodyB->getNumDof(&nkin, NULL);

    // Get the rotation matrix and its derivative
    int nparam = 0;
    const TacsScalar *CB, *CBd;
    bodyB->getRotationDeriv(&CB, &CBd, &nparam);

    // Add terms to the body residual
    TacsScalar fB[3], xBxI1[3];

    // Transform the constraint forces to the local frame
    for ( int k = 0; k < nparam; k++ ){
      matMult(&CBd[9*k], &fr[0], fB);
      addVecMat(alpha, fB, nkin, &D[3+k]); 

      // Add the terms from the cross-product
      crossProduct(1.0, xB, fB, xBxI1);
      addVecMat(alpha, xBxI1, nkin, &D[3*nkin+3+k]);
    }
  }
}

/*
  Compute the off-diagonal contributions to the Jacobian

  input:
  alpha:  the coefficient for the variables 
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  output:
  Dcon:    derivative of kinematic constraint w.r.t. kinematics vars
  Ddyn:    derivative of the dynamic contribution w.r.t. forces/torques
*/
void TACSDynSphericalJoint::getOffDiagJacobian( TacsScalar Dcon[], 
						TacsScalar Ddyn[], 
						double alpha,
						TACSRigidBody *rbody,
						const TacsScalar fr[] ){  
  if (rbody == bodyA){
    // Get the number of kinematic/dynamic state variables
    int nkin, ndyn;
    bodyA->getNumDof(&nkin, &ndyn);

    // Zero the off-diagonal entries
    memset(Dcon, 0, ncon*nkin*sizeof(TacsScalar));
    memset(Ddyn, 0, ndyn*ncon*sizeof(TacsScalar));

    // Get the rotation matrix and its derivative
    int nparam = 0;
    const TacsScalar *CA, *CAd;
    bodyA->getRotationDeriv(&CA, &CAd, &nparam);

    // Temporary matrix
    TacsScalar tmp[9];

    // Add the derivative of r
    addBlockIdent(alpha, nkin, &Dcon[0]);

    // Find the derivative of (CA^{T}*xA) w.r.t. the rotation matrix
    // parametrization
    for ( int k = 0; k < nparam; k++ ){
      matMultTrans(&CAd[9*k], xA, tmp);
      addVecMat(alpha, tmp, nkin, &Dcon[3+k]);
    }

    // Transform the constraint forces to the local frame
    addBlockMat(-alpha, CA, ncon, &Ddyn[0]); 

    // Add the terms from the cross-product
    TacsScalar col[3]; // Column of the rotation matrix

    // Extract the first column of the rotation matrix
    col[0] = CA[0];  col[1] = CA[3];  col[2] = CA[6];
    crossProduct(1.0, xA, col, &tmp[0]);

    // Extract the second column of the rotation matrix
    col[0] = CA[1];  col[1] = CA[4];  col[2] = CA[7];
    crossProduct(1.0, xA, col, &tmp[3]);

    // Extract the third column of the rotation matrix
    col[0] = CA[2];  col[1] = CA[5];  col[2] = CA[8];
    crossProduct(1.0, xA, col, &tmp[6]);
    addBlockMatTrans(-alpha, tmp, ncon, &Ddyn[3*ncon]);
  }
  else if (bodyB && rbody == bodyB){
    // Get the number of kinematic/dynamic state variables
    int nkin, ndyn;
    bodyA->getNumDof(&nkin, &ndyn);

    // Zero the off-diagonal entries
    memset(Dcon, 0, ncon*nkin*sizeof(TacsScalar));
    memset(Ddyn, 0, ndyn*ncon*sizeof(TacsScalar));

    // Get the rotation matrix and its derivative
    int nparam = 0;
    const TacsScalar *CB, *CBd;
    bodyB->getRotationDeriv(&CB, &CBd, &nparam);

    // Temporary matrix
    TacsScalar tmp[9];

    // Add the derivative of r
    addBlockIdent(-alpha, nkin, &Dcon[0]);

    // Find the derivative of (CB^{T}*xB) w.r.t. the rotation matrix
    // parametrization
    for ( int k = 0; k < nparam; k++ ){
      matMultTrans(&CBd[9*k], xB, tmp);
      addVecMat(-alpha, tmp, nkin, &Dcon[3+k]);
    }

    // Transform the constraint forces to the local frame
    addBlockMat(alpha, CB, ncon, &Ddyn[0]); 

    // Add the terms from the cross-product
    TacsScalar col[3]; // Column of the rotation matrix

    // Extract the first column of the rotation matrix
    col[0] = CB[0];  col[1] = CB[3];  col[2] = CB[6];
    crossProduct(1.0, xB, col, &tmp[0]);

    // Extract the second column of the rotation matrix
    col[0] = CB[1];  col[1] = CB[4];  col[2] = CB[7];
    crossProduct(1.0, xB, col, &tmp[3]);

    // Extract the third column of the rotation matrix
    col[0] = CB[2];  col[1] = CB[5];  col[2] = CB[8];
    crossProduct(1.0, xB, col, &tmp[6]);
    addBlockMatTrans(alpha, tmp, ncon, &Ddyn[3*ncon]);
  }
}

/*
  Add the derivative of the product of the adjoint variables and the
  contribution to the dynamic residuals from this constraint to the
  given input array.

  input:
  numDVs:  the number of design variables
  resAdj:  the adjoint variables
  body:    pointer to the body associated with the residuals
  fr:      the reaction forces/torques

  output:
  fdvSens: the array to which the result is added
*/
void TACSDynSphericalJoint::addBodyAdjResProduct( TacsScalar fdvSens[], 
						  int numDVs,
						  const TacsScalar dynAdj[],
						  TACSRigidBody *body,
						  const TacsScalar fr[] ){
  if (body == bodyA){
    // The body-fixed frames and vectors
    TACSRefFrame *FA0;
    TACSGibbsVector *rA0vec;
    const TacsScalar *CA0, *CA;
    const TacsScalar *rA0, *x0;

    // Retrieve the initial frames/position vectors
    bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);
    
    // Retrieve the rotation matricies
    bodyA->getRotation(&CA);
    FA0->getRotation(&CA0);

    // Retrieve the initial point locations
    rA0vec->getVector(&rA0);
    point->getVector(&x0);

    // Take the negative of the dynamic adjoint variables
    TacsScalar t[3];
    t[0] = -dynAdj[0];
    t[1] = -dynAdj[1];
    t[2] = -dynAdj[2];

    // - xA^{x}*CA*f[0] 
    // = - (CA0*(x0 - rA0))^{x}*CA*f[0]
    // = (CA*f[0])^{x}*CA0*(x0 - rA0)
    // psi^{T}*(CA*f[0])^{x}*CA0*(x0 - rA0)
    // = - ((CA*f[0])^{x}*psi)^{T}*CA0*(x0 - rA0)
    // Take the cross-product of the variables
    crossProductAdd(1.0, xA, &dynAdj[3], t);
    
    // Add the derivative to the result
    bodyA->addRotationAdjResProduct(fdvSens, numDVs, t, &fr[0]);

    // Compute the derivative contribution from xA
    TacsScalar fA[3];
    matMult(CA, &fr[0], fA);
    crossProduct(-1.0, fA, &dynAdj[3], t);

    // Add the term t^{T}d(CA0)/dx*s = 
    // ((CA*fr)^{x}*dynAdj)^{T}*dCA0/dx*(x0 - xA0)
    TacsScalar s[3];
    s[0] = x0[0] - rA0[0];
    s[1] = x0[1] - rA0[1];
    s[2] = x0[2] - rA0[2];
    FA0->addRotationAdjResProduct(fdvSens, numDVs, t, s);

    // Compute s = -CA0^{T}*((CA*fr)^{x}*dynAdj)
    matMultTrans(CA0, t, s);

    // Add the contribution from the derivatives of the points
    rA0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, s);
    point->addPointAdjResProduct(fdvSens, numDVs, 1.0, s);
  }
  else if (body == bodyB){
    // The body-fixed frames and vectors
    TACSRefFrame *FB0;
    TACSGibbsVector *rB0vec;
    const TacsScalar *CB0, *CB;
    const TacsScalar *rB0, *x0;

    // Retrieve the initial frames/position vectors
    bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);
    
    // Retrieve the rotation matricies
    bodyB->getRotation(&CB);
    FB0->getRotation(&CB0);

    // Retrieve the initial point locations
    rB0vec->getVector(&rB0);
    point->getVector(&x0);

    // Take the negative of the dynamic adjoint variables
    TacsScalar t[3];
    t[0] = dynAdj[0];
    t[1] = dynAdj[1];
    t[2] = dynAdj[2];

    // xB^{x}*CB*f[0] 
    // = (CB0*(x0 - rB0))^{x}*CB*f[0]
    // = - (CB*f[0])^{x}*CB0*(x0 - rB0)
    // - psi^{T}*(CB*f[0])^{x}*CB0*(x0 - rB0)
    // = ((CB*f[0])^{x}*psi)^{T}*CB0*(x0 - rB0)
    // Take the cross-product of the variables
    crossProductAdd(-1.0, xB, &dynAdj[3], t);
    
    // Add the derivative to the result
    bodyB->addRotationAdjResProduct(fdvSens, numDVs, t, &fr[0]);

    // Compute the derivative contribution from xB
    TacsScalar fB[3];
    matMult(CB, &fr[0], fB);
    crossProduct(1.0, fB, &dynAdj[3], t);

    // Add the term t^{T}d(CB0)/dx*s = 
    // ((CB*fr)^{x}*dynAdj)^{T}*dCB0/dx*(x0 - xB0)
    TacsScalar s[3];
    s[0] = x0[0] - rB0[0];
    s[1] = x0[1] - rB0[1];
    s[2] = x0[2] - rB0[2];
    FB0->addRotationAdjResProduct(fdvSens, numDVs, t, s);

    // Compute s = CB0^{T}*((CB*fr)^{x}*dynAdj)
    matMultTrans(CB0, t, s);

    // Add the contribution from the derivatives of the points
    rB0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, s);
    point->addPointAdjResProduct(fdvSens, numDVs, 1.0, s);
  }
}

/*
  The following class implements a revolute joint between two rigid
  bodies or between a rigid body and a fixed point.

  The revolute joint does to permit the transmission of a torque
  between the bodies along the revolute joint axis. The boides are
  fixed at points that are specified within each of the body-fixed
  frames.

  input:
  bodyA:  the pointer to the first body
  xA:     the position of the constraint point within the first body
  revA:   the axis of revolution in the first body
  bodyB:  the pointer to the second body
  xB:     the position of the constraint point within the second body
  revB:   the axis of revolution in the second body
*/
TACSDynRevoluteJoint::TACSDynRevoluteJoint( TACSRigidBody *_bodyA,
					    TACSRigidBody *_bodyB,
					    TACSGibbsVector *_point,
					    TACSGibbsVector *_rev ){
  // Copy over the pointer to the first body object
  bodyA = _bodyA;
  bodyA->incref();

  // Copy the pointer to the second body
  bodyB = _bodyB;
  bodyB->incref();

  // Get the point
  point = _point;
  point->incref();

  // Get the direction of revolution
  rev = _rev;
  rev->incref();

  // Update the points to define revA/revB and xA/xB etc.
  int init_ek = 1;
  updatePoints(init_ek);
}

/*
  Free/deallocate the joint
*/
TACSDynRevoluteJoint::~TACSDynRevoluteJoint(){
  bodyA->decref();
  bodyB->decref();
  point->decref();
  rev->decref();
}

/*
  Read the data from the given initial point vectors/locations and
  re-compute the internal data that is requied to evaluate the
  kinematic constraints.

  This code computes the required vectors in the body-fixed coordinate
  frames. In particular, the code computes:
  
  xA = CA0*(x0 - rA0)
  xB = CB0*(x0 - rB0)

  where x0 is the attachment point, and rA0 and rB0 are the initial
  points of bodies A and B in the global (inertial) reference frame.
  In addition, the code computes the initial revolute direction in
  the A and B reference frames, given by:

  revA = CA0*rev
  revB = CB0*rev

  where rev is the revolute direction in the global frame.
*/
void TACSDynRevoluteJoint::updatePoints( int init_ek ){
  // The body-fixed frames and vectors
  TACSRefFrame *FA0, *FB0;
  TACSGibbsVector *rA0vec, *rB0vec;
  const TacsScalar *CA0, *CB0;
  const TacsScalar *rA0, *rB0, *x0;

  // Retrieve the initial frames/position vectors
  bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);
  bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);

  // Retrieve the rotation matricies
  FA0->getRotation(&CA0);
  FB0->getRotation(&CB0);

  // Retrieve the initial point locations
  rA0vec->getVector(&rA0);
  rB0vec->getVector(&rB0);
  point->getVector(&x0);

  // Set the position of the constrained point within the first
  // body-fixed frame
  TacsScalar t[3];
  t[0] = x0[0] - rA0[0];
  t[1] = x0[1] - rA0[1];
  t[2] = x0[2] - rA0[2];
  matMult(CA0, t, xA);

  // Get the revolute direction
  const TacsScalar *rv;
  rev->getVector(&rv);

  // Find the axis of revolution in the first body frame
  matMult(CA0, rv, revA);

  // Set the position of the constrained point within the second
  // body-fixed frame
  t[0] = x0[0] - rB0[0];
  t[1] = x0[1] - rB0[1];
  t[2] = x0[2] - rB0[2];
  matMult(CB0, t, xB);

  // Find the axis of revolution in the second body frame
  matMult(CB0, rv, revB);

  // Find the minimum absolute component of revB along any
  // coordinate direction. Set the vector components of
  // ek along this direction to maximize orthogonality among
  // the coordinate directions. For the purpose of optimization,
  // this direction is fixed at initialization.
  if (init_ek){
    ek[0] = ek[1] = ek[2] = 0.0;
    if ((fabs(revB[0]) <= fabs(revB[1])) && 
	(fabs(revB[0]) <= fabs(revB[2]))){
      ek[0] = 1.0;
    }
    else if ((fabs(revB[1]) <= fabs(revB[0])) && 
	     (fabs(revB[1]) <= fabs(revB[2]))){
      ek[1] = 1.0;
    }
    else {
      ek[2] = 1.0;
    }
  }

  // Compute the eB1 and eB2 directions based on ek
  crossProduct(1.0, revB, ek, eB2);
  crossProduct(1.0, eB2, revB, eB1);
}

/*
  Set the design variable values into the object
*/
void TACSDynRevoluteJoint::setDesignVars( const TacsScalar dvs[],
					  int numDVs ){
  bodyA->setDesignVars(dvs, numDVs);
  if (bodyB){
    bodyB->setDesignVars(dvs, numDVs);
  }
  rev->setDesignVars(dvs, numDVs);
  point->setDesignVars(dvs, numDVs);
  updatePoints();
}

/*
  Retrieve the design variable values from the object
*/
void TACSDynRevoluteJoint::getDesignVars( TacsScalar dvs[], int numDVs ){
  bodyA->getDesignVars(dvs, numDVs);
  if (bodyB){
    bodyB->getDesignVars(dvs, numDVs);
  }
  rev->getDesignVars(dvs, numDVs);
  point->getDesignVars(dvs, numDVs);
}

/*
  Retrieve the bodies associated with this kinematic constraint
*/
void TACSDynRevoluteJoint::getBodies( TACSRigidBody **_bodyA,
				      TACSRigidBody **_bodyB ){
  *_bodyA = bodyA;
  *_bodyB = bodyB;
}

/*
  Compute the residual of the constraint equation.

  The constraint equation is expressed in the global coordinate frame
  as follows:

  CA^{T}*xA + rA - CB^{T}*xB - rB = 0
  g = 0

  Note that the use of the transpose converts the components to the
  global inertial reference frame.

  input:
  q:    the kinematic/dynamic states
  
  output:
  res:  the residual
*/
void TACSDynRevoluteJoint::getResidual( TacsScalar res[], 
					const TacsScalar fr[] ){
  // Temporary vector
  TacsScalar tmp[3];

  // Retrieve the position vector and rotation matrix
  const TacsScalar *rA, *CA;
  bodyA->getVariables(&rA, NULL, NULL);
  bodyA->getRotation(&CA);

  // Form the residual of the first constraint equation
  // CA^{T}*xA + rA - CB^{T}*xB - rB - xI = 0
  matMultTrans(CA, xA, &res[0]);
  vecAxpy(1.0, rA, &res[0]);

  // Retrieve the position vector and rotation matrix
  // for the second body
  const TacsScalar *rB, *CB;
  bodyB->getVariables(&rB, NULL, NULL);
  bodyB->getRotation(&CB);
  
  matMultTrans(CB, xB, tmp);
  vecAxpy(1.0, rB, tmp);
  vecAxpy(-1.0, tmp, &res[0]);
  
  // Form the constraint for the torques
  // dot(revA, gr) = revA^{T}*CA*gr
  matMult(CA, &fr[3], tmp);
  res[3] = vecDot(revA, tmp);

  // Compute the dot product of of eB1 and eB2 with revA
  TacsScalar s[3];
  matMultTrans(CA, revA, tmp);
  matMultTrans(CB, eB1, s);
  res[4] = vecDot(tmp, s);

  matMultTrans(CB, eB2, s);
  res[5] = vecDot(tmp, s);
}

/*
  Compute the derivative of the product of the adjoint vector and the
  residual vector of the kinematic constraints.

  The output is added to the given sensitivity array.

  input:
  numDVs:  the number of design variables
  resAdj:  the adjoint variables
  fr:      the reaction forces/torques

  output:
  fdvSens: the array to which the result is added
*/
void TACSDynRevoluteJoint::addAdjResProduct( TacsScalar fdvSens[], 
					     int numDVs,
					     const TacsScalar resAdj[],
					     const TacsScalar fr[] ){
  // The body-fixed frames and vectors
  TACSRefFrame *FA0, *FB0;
  TACSGibbsVector *rA0vec, *rB0vec;
  const TacsScalar *CA0, *CB0, *CA, *CB;
  const TacsScalar *rA0, *rB0, *x0, *rv;

  // Retrieve the initial frames/position vectors
  bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);
  bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);

  // Retrieve the rotation matricies
  bodyA->getRotation(&CA);
  bodyB->getRotation(&CB);
  FA0->getRotation(&CA0);
  FB0->getRotation(&CB0);

  // Retrieve the initial point locations
  rA0vec->getVector(&rA0);
  rB0vec->getVector(&rB0);
  point->getVector(&x0);
  rev->getVector(&rv);

  // Add the derivative of from the points
  TacsScalar t[3], s[3];
  matMult(CA, &resAdj[0], s); // s = CA*resAdj
  matMultTrans(CA0, s, t); // t = CA0^{T}*CA*resAdj
  point->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);
  rA0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, t);

  // Add the derivative from the initial reference frame position
  t[0] = x0[0] - rA0[0];
  t[1] = x0[1] - rA0[1];
  t[2] = x0[2] - rA0[2];
  FA0->addRotationAdjResProduct(fdvSens, numDVs, s, t);

  // Add the derivative of: resAdj^{T}*CA^{T}*xA = xA^{T}*CA*resAdj
  bodyA->addRotationAdjResProduct(fdvSens, numDVs, xA, &resAdj[0]);

  // Add the contribution from the second body 
  // Add the derivative of from the points
  matMult(CB, &resAdj[0], s); // s = CB*resAdj
  matMultTrans(CB0, s, t); // t = CB0^{T}*CB*resAdj
  point->addPointAdjResProduct(fdvSens, numDVs, -1.0, t);
  rB0vec->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);
  
  // Add the derivative from the initial reference frame position
  t[0] = -(x0[0] - rB0[0]);
  t[1] = -(x0[1] - rB0[1]);
  t[2] = -(x0[2] - rB0[2]);
  FB0->addRotationAdjResProduct(fdvSens, numDVs, s, t);
  
  // Add the derivative of: -resAdj^{T}*CB^{T}*xB = xB^{T}*CB*resAdj
  TacsScalar xBneg[3];
  xBneg[0] = -xB[0];
  xBneg[1] = -xB[1];
  xBneg[2] = -xB[2];
  bodyB->addRotationAdjResProduct(fdvSens, numDVs, xBneg, &resAdj[0]);  

  // Add the contributions from the fourth constraint:
  // -------------------------------------------------
  // d(resAdj[3]*revA^{T}*CA*fr[3])/dx
  t[0] = resAdj[3]*fr[3];
  t[1] = resAdj[3]*fr[4];
  t[2] = resAdj[3]*fr[5];
  bodyA->addRotationAdjResProduct(fdvSens, numDVs, revA, t);

  // Add the contributions from the derivative of revA = CA0*rv
  matMult(CA, t, s);
  FA0->addRotationAdjResProduct(fdvSens, numDVs, s, rv);

  // Multiply by C0^{T} to get the sensitivity for rev
  matMultTrans(CA0, s, t);
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);

  // Add the contribution from the fifth constraint:
  // -----------------------------------------------
  // Note that eB1 = ek*dot(rv, rv) - revB*dot(revB, ek)
  matMultTrans(CA, revA, t);
  vecScale(resAdj[4], t);
  bodyB->addRotationAdjResProduct(fdvSens, numDVs, eB1, t);

  // s = CB*CA^{T}*revA
  matMult(CB, t, s);

  // Add the derivative from d(s^{T}*ek*dot(rv, rv))/dx
  TacsScalar d = vecDot(s, ek);
  rev->addPointAdjResProduct(fdvSens, numDVs, 2.0*d, rv);

  // Add the derivative: - d/dx(s^{T}*CB0*rv*(ek^{T}*CB0*rv))
  TacsScalar d1 = -vecDot(s, revB);
  TacsScalar d2 = -vecDot(ek, revB);

  // t = d1*ek
  t[0] = d1*ek[0];
  t[1] = d1*ek[1];
  t[2] = d1*ek[2];

  // Add the derivative: t*d(CB0)/dx*rv
  FB0->addRotationAdjResProduct(fdvSens, numDVs, t, rv);

  // Add the derivative t*CB0*d(rv)/dx
  TacsScalar r[3];
  matMultTrans(CB0, t, r);
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, r);

  // t = d2*s
  t[0] = d2*s[0];
  t[1] = d2*s[1];
  t[2] = d2*s[2];

  // Add the derivative: t*d(CB0)/dx*rv
  FB0->addRotationAdjResProduct(fdvSens, numDVs, t, rv);

  // Add the derivative t*CB0*d(rv)/dx
  matMultTrans(CB0, t, r);
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, r);

  // Compute the contribution to the derivative from the other side
  matMultTrans(CB, eB1, t);
  vecScale(resAdj[4], t);

  // t = CB^{T}*eB1
  // Add the derivative: t^{T}*CA^{T}*revA
  bodyA->addRotationAdjResProduct(fdvSens, numDVs, revA, t);
  
  // s = CA*t
  matMult(CA, t, s);

  // Add the derivative: s^{T}*CA0*rv
  FA0->addRotationAdjResProduct(fdvSens, numDVs, s, rv);

  // t = CA0^{T}*s
  matMultTrans(CA0, s, t);
  
  // Add the derivative: t = CA0^{T}*s
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);

  // Add the contribution from the sixth consraint:
  // ----------------------------------------------
  matMultTrans(CA, revA, t);
  vecScale(resAdj[5], t);
  bodyB->addRotationAdjResProduct(fdvSens, numDVs, eB2, t);

  // s = CB*CA^{T}*revA
  matMult(CB, t, s);

  // t = ek^{x}*s
  crossProduct(1.0, ek, s, t);
  
  // Add the derivative: t*CB0*rv, where t = e^{x}*CB*CA^{T}*revA
  FB0->addRotationAdjResProduct(fdvSens, numDVs, t, rv);

  // Multiply s = CB0^{T}*t
  matMultTrans(CB0, t, s);
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, s);

  // Compute from the other side
  matMultTrans(CB, eB2, t);
  vecScale(resAdj[5], t);

  // t = CB^{T}*eB2
  // Add the derivative: t^{T}*CA^{T}*revA
  bodyA->addRotationAdjResProduct(fdvSens, numDVs, revA, t);
  
  // s = CA*t
  matMult(CA, t, s);

  // Add the derivative: s^{T}*CA0*rv
  FA0->addRotationAdjResProduct(fdvSens, numDVs, s, rv);

  // t = CA0^{T}*s
  matMultTrans(CA0, s, t);
  
  // Add the derivative: t = CA0^{T}*s
  rev->addPointAdjResProduct(fdvSens, numDVs, 1.0, t);
}

/*
  Compute the Jacobian of the constraint equation w.r.t. the dynamic
  variables (e.g. forces/moments)

  input:
  alpha:  the coefficient for the variables 
  fr:     the dynamic variables

  output:
  Dff:  the Jacobian w.r.t. the dynamic variables
*/
void TACSDynRevoluteJoint::getJacobian( TacsScalar D[], 
					double alpha,
					const TacsScalar fr[] ){
  // Set this part of the Jacobian to zero
  memset(D, 0, ncon*ncon*sizeof(TacsScalar));

  // Retrieve the rotation matrix
  TacsScalar tmp[3];
  const TacsScalar *CA;
  bodyA->getRotation(&CA);
  matMultTrans(CA, revA, tmp);

  for ( int k = 0; k < 3; k++ ){
    D[3*(ncon+1)+k] = alpha*tmp[k];
  }
}

/*
  Add terms to the residual of the governing equations to account for
  the reaction forces/torques internal to the multibody dynamics
  system.

  input:
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  input/output:
  rdyn:   the dynamic residuals
*/
void TACSDynRevoluteJoint::addBodyResidual( TacsScalar rdyn[],
					    TACSRigidBody *rbody,
					    const TacsScalar fr[] ){
  if (rbody == bodyA){
    // Get the rotation matrix
    const TacsScalar *CA;
    bodyA->getRotation(&CA);

    // Add terms to the body residual
    TacsScalar fA[3], gA[3];

    // Transform the constraint forces to the local frame
    matMult(CA, &fr[0], fA);
    matMult(CA, &fr[3], gA);
    vecAxpy(-1.0, fA, &rdyn[0]);
    
    // Add the torque
    crossProductAdd(-1.0, xA, fA, &rdyn[3]);
    vecAxpy(-1.0, gA, &rdyn[3]);
  }
  else if (rbody == bodyB){
    // Get the rotation matrix
    const TacsScalar *CB;
    bodyB->getRotation(&CB);

    // Add terms to the body residual
    TacsScalar fB[3], gB[3];

    // Transform the constraint forces to the local frame
    matMult(CB, &fr[0], fB);
    matMult(CB, &fr[3], gB);
    vecAxpy(1.0, fB, &rdyn[0]);
    
    // Add the torque
    crossProductAdd(1.0, xB, fB, &rdyn[3]);
    vecAxpy(1.0, gB, &rdyn[3]);
  }
}

/*
  Add components to the Jacobian of the dynamics with respect to the
  kinematic variables of the rigid body.

  input:
  alpha:  the coefficient for the variables 
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  input/output:
  D:      the Jacobian of the dynamics w.r.t. the kinematic variables
*/
void TACSDynRevoluteJoint::addBodyJacobian( TacsScalar D[],
					    double alpha,
					    TACSRigidBody *rbody,
					    const TacsScalar fr[] ){
  // Add the contributions to the Jacobian
  if (rbody == bodyA){
    // Get the number of kinematic terms
    int nkin = 0;
    bodyA->getNumDof(&nkin, NULL);

    // Get the rotation matrix
    int nparam = 0;
    const TacsScalar *CA, *CAd;
    bodyA->getRotationDeriv(&CA, &CAd, &nparam);

    // Add terms to the body residual
    TacsScalar fA[3], gA[3], xAxI[3];

    for ( int k = 0; k < nparam; k++ ){
      // Transform the constraint forces to the local frame
      matMult(&CAd[9*k], &fr[0], fA);
      addVecMat(-alpha, fA, nkin, &D[3+k]); 

      // Add the contributions from the torque
      matMult(&CAd[9*k], &fr[3], gA);
      addVecMat(-alpha, gA, nkin, &D[3*nkin+3+k]);

      // Add the terms from the cross-product
      crossProduct(1.0, xA, fA, xAxI);
      addVecMat(-alpha, xAxI, nkin, &D[3*nkin+3+k]);
    }
  }
  else if (bodyB && rbody == bodyB){
   // Get the number of kinematic terms
    int nkin = 0;
    bodyB->getNumDof(&nkin, NULL);

    // Get the rotation matrix
    int nparam = 0;
    const TacsScalar *CB, *CBd;
    bodyB->getRotationDeriv(&CB, &CBd, &nparam);

    // Add terms to the body residual
    TacsScalar fB[3], gB[3], xBxI[3];

    for ( int k = 0; k < nparam; k++ ){
      // Transform the constraint forces to the local frame
      matMult(&CBd[9*k], &fr[0], fB);
      addVecMat(alpha, fB, nkin, &D[3+k]); 

      // Add the contributions from the torque
      matMult(&CBd[9*k], &fr[3], gB);
      addVecMat(alpha, gB, nkin, &D[3*nkin+3+k]);

      // Add the terms from the cross-product
      crossProduct(1.0, xB, fB, xBxI);
      addVecMat(alpha, xBxI, nkin, &D[3*nkin+3+k]);
    }
  }
}

/*
  Compute the off-diagonal contributions to the Jacobian

  input:
  alpha:  the coefficient for the variables 
  body:   a pointer to the rigid body
  fr:     the forces/torque variables

  output:
  Dcon:    derivative of kinematic constraint w.r.t. kinematics vars
  Ddyn:    derivative of the dynamic contribution w.r.t. forces/torques
*/
void TACSDynRevoluteJoint::getOffDiagJacobian( TacsScalar Dcon[], 
					       TacsScalar Ddyn[], 
					       double alpha,
					       TACSRigidBody *rbody,
					       const TacsScalar fr[] ){
  if (rbody == bodyA){
    // Get the number of kinematic/dynamic state variables
    int nkin, ndyn;
    bodyA->getNumDof(&nkin, &ndyn);

    // Zero the off-diagonal entries
    memset(Dcon, 0, ncon*nkin*sizeof(TacsScalar));
    memset(Ddyn, 0, ndyn*ncon*sizeof(TacsScalar));

    // Get the rotation matrix and its derivative
    int nparam = 0;
    const TacsScalar *CA, *CAd;
    bodyA->getRotationDeriv(&CA, &CAd, &nparam);
    const TacsScalar *CB;
    bodyB->getRotation(&CB);

    // Temporary matrix
    TacsScalar tmp[9];

    // Add the derivative of r
    addBlockIdent(alpha, nkin, &Dcon[0]);
    
    // Compute the components of the orthogonal vectors for the
    // kinematic constraint in the local frame
    TacsScalar s1[3], s2[3];
    matMultTrans(CB, eB1, s1);
    matMultTrans(CB, eB2, s2);

    for ( int k = 0; k < nparam; k++ ){
      // Find the derivative of (CA^{T}*xA) w.r.t. the rotation matrix
      // parametrization
      matMultTrans(&CAd[9*k], xA, tmp);
      addVecMat(alpha, tmp, nkin, &Dcon[3+k]);

      // Take the derivative of the constraint 
      // dot(revA, gr) = revA^{T}*CA*gr = 0
      matMult(&CAd[9*k], &fr[3], tmp);
      Dcon[3*nkin+3+k] = alpha*vecDot(revA, tmp);

      // Take the derivative of the constraint 
      // dot(revA, eB1) = eB1^{T}*CB*CA^{T}*revA
      matMultTrans(&CAd[9*k], revA, tmp);
      Dcon[4*nkin+3+k] = alpha*vecDot(s1, tmp);

      // Take the derivative of the constraint 
      // dot(revA, eB1) = eB1^{T}*CB*CA^{T}*revA
      matMultTrans(&CAd[9*k], revA, tmp);
      Dcon[5*nkin+3+k] = alpha*vecDot(s2, tmp);
    }

    // Transform the constraint forces to the local frame
    addBlockMat(-alpha, CA, ncon, &Ddyn[0]); 

    // Transform the constraint torques to the local frame
    addBlockMat(-alpha, CA, ncon, &Ddyn[3*ncon + 3]); 

    // Add the terms from the cross-product
    TacsScalar col[3]; // Column of the rotation matrix

    // Extract the first column of the rotation matrix
    col[0] = CA[0];  col[1] = CA[3];  col[2] = CA[6];
    crossProduct(1.0, xA, col, &tmp[0]);

    // Extract the second column of the rotation matrix
    col[0] = CA[1];  col[1] = CA[4];  col[2] = CA[7];
    crossProduct(1.0, xA, col, &tmp[3]);

    // Extract the third column of the rotation matrix
    col[0] = CA[2];  col[1] = CA[5];  col[2] = CA[8];
    crossProduct(1.0, xA, col, &tmp[6]);
    addBlockMatTrans(-alpha, tmp, ncon, &Ddyn[3*ncon]);
  }
  else if (rbody == bodyB){
    // Get the number of kinematic/dynamic state variables
    int nkin, ndyn;
    bodyB->getNumDof(&nkin, &ndyn);

    // Zero the off-diagonal entries
    memset(Dcon, 0, ncon*nkin*sizeof(TacsScalar));
    memset(Ddyn, 0, ndyn*ncon*sizeof(TacsScalar));

    // Get the rotation matrix and its derivative
    int nparam = 0;
    const TacsScalar *CA;
    bodyA->getRotation(&CA);
    const TacsScalar *CB, *CBd;
    bodyB->getRotationDeriv(&CB, &CBd, &nparam);

    // Temporary matrix
    TacsScalar tmp[9];

    // Add the derivative of r
    addBlockIdent(-alpha, nkin, &Dcon[0]);
    
    // Compute the components of revA in the global frame
    TacsScalar s[3];
    matMultTrans(CA, revA, s);

    for ( int k = 0; k < nparam; k++ ){
      // Find the derivative of (CB^{T}*xB) w.r.t. the rotation matrix
      // parametrization
      matMultTrans(&CBd[9*k], xB, tmp);
      addVecMat(-alpha, tmp, nkin, &Dcon[3+k]);

      // Take the derivative of the constraint 
      // dot(revA, eB1) = eB1^{T}*CB*CA^{T}*revA
      matMultTrans(&CBd[9*k], eB1, tmp);
      Dcon[4*nkin+3+k] = alpha*vecDot(s, tmp);

      // Take the derivative of the constraint 
      // dot(revA, eB1) = eB1^{T}*CB*CA^{T}*revA
      matMultTrans(&CBd[9*k], eB2, tmp);
      Dcon[5*nkin+3+k] = alpha*vecDot(s, tmp);
    }

    // Transform the constraint forces to the local frame
    addBlockMat(alpha, CB, ncon, &Ddyn[0]); 

    // Transform the constraint torques to the local frame
    addBlockMat(alpha, CB, ncon, &Ddyn[3*ncon + 3]); 

    // Add the terms from the cross-product
    TacsScalar col[3]; // Column of the rotation matrix

    // Extract the first column of the rotation matrix
    col[0] = CB[0];  col[1] = CB[3];  col[2] = CB[6];
    crossProduct(1.0, xB, col, &tmp[0]);

    // Extract the second column of the rotation matrix
    col[0] = CB[1];  col[1] = CB[4];  col[2] = CB[7];
    crossProduct(1.0, xB, col, &tmp[3]);

    // Extract the third column of the rotation matrix
    col[0] = CB[2];  col[1] = CB[5];  col[2] = CB[8];
    crossProduct(1.0, xB, col, &tmp[6]);
    addBlockMatTrans(alpha, tmp, ncon, &Ddyn[3*ncon]);
  }
}

/*
  Add the derivative of the product of the adjoint variables and the
  contribution to the dynamic residuals from this constraint to the
  given input array.

  input:
  numDVs:  the number of design variables
  resAdj:  the adjoint variables
  body:    pointer to the body associated with the residuals
  fr:      the reaction forces/torques

  output:
  fdvSens: the array to which the result is added
*/
void TACSDynRevoluteJoint::addBodyAdjResProduct( TacsScalar fdvSens[], 
						 int numDVs,
						 const TacsScalar dynAdj[],
						 TACSRigidBody *body,
						 const TacsScalar fr[] ){
  if (body == bodyA){
    // The body-fixed frames and vectors
    TACSRefFrame *FA0;
    TACSGibbsVector *rA0vec;
    const TacsScalar *CA0, *CA;
    const TacsScalar *rA0, *x0;

    // Retrieve the initial frames/position vectors
    bodyA->getInitVariables(&FA0, &rA0vec, NULL, NULL);
    
    // Retrieve the rotation matricies
    bodyA->getRotation(&CA);
    FA0->getRotation(&CA0);

    // Retrieve the initial point locations
    rA0vec->getVector(&rA0);
    point->getVector(&x0);

    // Take the negative of the dynamic adjoint variables
    TacsScalar t[3];
    t[0] = -dynAdj[0];
    t[1] = -dynAdj[1];
    t[2] = -dynAdj[2];

    // - xA^{x}*CA*f[0] 
    // = - (CA0*(x0 - rA0))^{x}*CA*f[0]
    // = (CA*f[0])^{x}*CA0*(x0 - rA0)
    // psi^{T}*(CA*f[0])^{x}*CA0*(x0 - rA0)
    // = - ((CA*f[0])^{x}*psi)^{T}*CA0*(x0 - rA0)
    // Take the cross-product of the variables
    crossProductAdd(1.0, xA, &dynAdj[3], t);
    
    // Add the derivative to the result
    bodyA->addRotationAdjResProduct(fdvSens, numDVs, t, &fr[0]);

    // Add the contribution from the reaction torque
    t[0] = -dynAdj[3];
    t[1] = -dynAdj[4];
    t[2] = -dynAdj[5];
    bodyA->addRotationAdjResProduct(fdvSens, numDVs, t, &fr[3]);

    // Compute the derivative contribution from xA
    TacsScalar fA[3];
    matMult(CA, &fr[0], fA);
    crossProduct(-1.0, fA, &dynAdj[3], t);

    // Add the term t^{T}d(CA0)/dx*s = 
    // ((CA*fr)^{x}*dynAdj)^{T}*dCA0/dx*(x0 - xA0)
    TacsScalar s[3];
    s[0] = x0[0] - rA0[0];
    s[1] = x0[1] - rA0[1];
    s[2] = x0[2] - rA0[2];
    FA0->addRotationAdjResProduct(fdvSens, numDVs, t, s);

    // Compute s = -CA0^{T}*((CA*fr)^{x}*dynAdj)
    matMultTrans(CA0, t, s);

    // Add the contribution from the derivatives of the points
    rA0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, s);
    point->addPointAdjResProduct(fdvSens, numDVs, 1.0, s);
  }
  else if (body == bodyB){
    // The body-fixed frames and vectors
    TACSRefFrame *FB0;
    TACSGibbsVector *rB0vec;
    const TacsScalar *CB0, *CB;
    const TacsScalar *rB0, *x0;

    // Retrieve the initial frames/position vectors
    bodyB->getInitVariables(&FB0, &rB0vec, NULL, NULL);
    
    // Retrieve the rotation matricies
    bodyB->getRotation(&CB);
    FB0->getRotation(&CB0);

    // Retrieve the initial point locations
    rB0vec->getVector(&rB0);
    point->getVector(&x0);

    // Take the negative of the dynamic adjoint variables
    TacsScalar t[3];
    t[0] = dynAdj[0];
    t[1] = dynAdj[1];
    t[2] = dynAdj[2];

    // xB^{x}*CB*f[0] 
    // = (CB0*(x0 - rB0))^{x}*CB*f[0]
    // = - (CB*f[0])^{x}*CB0*(x0 - rB0)
    // - psi^{T}*(CB*f[0])^{x}*CB0*(x0 - rB0)
    // = ((CB*f[0])^{x}*psi)^{T}*CB0*(x0 - rB0)
    // Take the cross-product of the variables
    crossProductAdd(-1.0, xB, &dynAdj[3], t);
    
    // Add the derivative to the result
    bodyB->addRotationAdjResProduct(fdvSens, numDVs, t, &fr[0]);

    // Add the contribution from the torque
    bodyB->addRotationAdjResProduct(fdvSens, numDVs, 
				    &dynAdj[3], &fr[3]);

    // Compute the derivative contribution from xB
    TacsScalar fB[3];
    matMult(CB, &fr[0], fB);
    crossProduct(1.0, fB, &dynAdj[3], t);

    // Add the term t^{T}d(CB0)/dx*s = 
    // ((CB*fr)^{x}*dynAdj)^{T}*dCB0/dx*(x0 - xB0)
    TacsScalar s[3];
    s[0] = x0[0] - rB0[0];
    s[1] = x0[1] - rB0[1];
    s[2] = x0[2] - rB0[2];
    FB0->addRotationAdjResProduct(fdvSens, numDVs, t, s);

    // Compute s = CB0^{T}*((CB*fr)^{x}*dynAdj)
    matMultTrans(CB0, t, s);

    // Add the contribution from the derivatives of the points
    rB0vec->addPointAdjResProduct(fdvSens, numDVs, -1.0, s);
    point->addPointAdjResProduct(fdvSens, numDVs, 1.0, s);
  }
}
