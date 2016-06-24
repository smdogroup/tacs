#include "MITC9.h"
#include "TACSElementAlgebra.h"
#include "FElibrary.h"

/*
  Rigid-body dynamics routines for TACS
  
  Copyright (c) 2015-2016 Graeme Kennedy. All rights reserved. 
*/

/* 
   Return the number of displacements
*/
int MITC9::numDisplacements(){ return NUM_DISPS; }

/*
  Return the number of stresses
*/
int MITC9::numStresses(){ return NUM_STRESSES; }
  
/*
  Return the number of extras
*/
int MITC9::numExtras(){ return NUM_EXTRAS; }

/*
  Return the number of FE nodes
*/
int MITC9::numNodes(){ return NUM_NODES; }

/*
  Return the ElementType
*/
ElementType MITC9::getElementType(){ return SHELL; }

/* 
   Set up the internal static data for the names of the element,
   displacements, stresses, strains and extra variables, respectively.
*/
const char * MITC9::elemName = "MITC9";
  
const char * MITC9::dispNames[] = { "u0", "v0", "w0", 
				    "eta", "rotx", "roty", "rotz", "lam" };
  
const char * MITC9::stressNames[] = { "sx0", "sy0", "sxy0", 
				      "sx1", "sy1", "sxy1", 
				      "syz0", "sxz0" };
  
const char * MITC9::strainNames[] = { "ex0", "ey0", "exy0", 
				      "ex1", "ey1", "exy1", 
				      "eyz0", "exz0" };
  
const char * MITC9::extraNames[] = { "lambda", "buckling",
				     "dv1", "dv2" };

/*
  Returns the elementName
*/
const char * MITC9::elementName(){ 
  return elemName;
}

/*
  Returns the displacement names
*/
const char * MITC9::displacementName( int i ){
  if (i >= 0 && i < NUM_DISPS){
    return dispNames[i];
  }
  return "";
}

/*
  Returns the name of the stresses
*/
const char * MITC9::stressName( int i ){ 
  if (i >= 0 && i < NUM_STRESSES){
    return stressNames[i];
  }
  return "";
}

/*
  Returns the name of the strains
*/
const char * MITC9::strainName( int i ){
  if (i >= 0 && i < NUM_STRESSES){
    return strainNames[i];
  }
  return "";
}

/*
  Returns the extra names
*/
const char * MITC9::extraName( int i ){
  if (i >= 0 && i < NUM_EXTRAS){
    return extraNames[i];
  }
  return "";
}

/*
  Set the design variable values
*/
void MITC9::setDesignVars( const TacsScalar dvs[], int numDVs ){
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Get the design variable values
*/
void MITC9::getDesignVars( TacsScalar dvs[], int numDVs ){
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Get the design variable range
*/
void MITC9::getDesignVarRange( TacsScalar lb[], 
                               TacsScalar ub[], int numDVs ){
  stiff->getDesignVarRange(lb, ub, numDVs);
}

/*
  Evaluate the shape functions of the element given the u/v
  coordinates of the point

  input:
  u:   the first parametric coordinate
  v:   the second parametric coordinate
  
  output:
  N:   the shape functions
*/
static inline void computeShapeFunc( const double u,
				     const double v,
				     double N[] ){
  // Compute the shape functions
  double nu[3];
  nu[0] = -0.5*u*(1.0 - u);
  nu[1] = (1.0 - u)*(1.0 + u);
  nu[2] = 0.5*(1.0 + u)*u;

  double nv[3];
  nv[0] = -0.5*v*(1.0 - v);
  nv[1] = (1.0 - v)*(1.0 + v);
  nv[2] = 0.5*(1.0 + v)*v;

  // Compute the shape functions
  N[0] = nu[0]*nv[0];
  N[1] = nu[1]*nv[0];
  N[2] = nu[2]*nv[0];

  N[3] = nu[0]*nv[1];
  N[4] = nu[1]*nv[1];
  N[5] = nu[2]*nv[1];

  N[6] = nu[0]*nv[2];
  N[7] = nu[1]*nv[2];
  N[8] = nu[2]*nv[2];
}

/*
  Evaluate the derivatives of the shape functions of the element at
  given the u/v coordinates of the point

  input:
  u:    the first parametric coordinate
  v:    the second parametric coordinate
  
  output:
  Na:   the derivative of the shape functions w.r.t. u
  Nb:   the derivative of the shape functions w.r.t. v
*/
static inline void computeShapeFunc( const double u,
				     const double v,
				     double Na[],
				     double Nb[] ){
  // Compute the shape functions and their derivatives
  double nu[3], dnu[3];
  nu[0] = -0.5*u*(1.0 - u);
  nu[1] = (1.0 - u)*(1.0 + u);
  nu[2] = 0.5*(1.0 + u)*u;
  
  dnu[0] = -0.5 + u;
  dnu[1] = -2.0*u;
  dnu[2] = 0.5 + u;

  double nv[3], dnv[3];
  nv[0] = -0.5*v*(1.0 - v);
  nv[1] = (1.0 - v)*(1.0 + v);
  nv[2] = 0.5*(1.0 + v)*v;

  dnv[0] = -0.5 + v;
  dnv[1] = -2.0*v;
  dnv[2] = 0.5 + v;

  // Compute the derivative of the shape funcs w.r.t. u
  Na[0] = dnu[0]*nv[0];
  Na[1] = dnu[1]*nv[0];
  Na[2] = dnu[2]*nv[0];

  Na[3] = dnu[0]*nv[1];
  Na[4] = dnu[1]*nv[1];
  Na[5] = dnu[2]*nv[1];

  Na[6] = dnu[0]*nv[2];
  Na[7] = dnu[1]*nv[2];
  Na[8] = dnu[2]*nv[2];

  // Compute the derivative of the shape funcs w.r.t. v
  Nb[0] = nu[0]*dnv[0];
  Nb[1] = nu[1]*dnv[0];
  Nb[2] = nu[2]*dnv[0];

  Nb[3] = nu[0]*dnv[1];
  Nb[4] = nu[1]*dnv[1];
  Nb[5] = nu[2]*dnv[1];

  Nb[6] = nu[0]*dnv[2];
  Nb[7] = nu[1]*dnv[2];
  Nb[8] = nu[2]*dnv[2];
}

/*
  Compute the inner product of a set of shape functions (or
  their derivative) with a 3*NUM_NODES array of nodes, variables
  etc.

  input:
  Na:   the shape function
  X:    the variables 

  ouput:
  Xa:   the inner product = Na^{T}*X
*/
static inline void innerProduct( const double Na[],
				 const TacsScalar X[],
				 TacsScalar Xa[] ){
  Xa[0] = (Na[0]*X[0] + Na[1]*X[3] + Na[2]*X[6] + 
	   Na[3]*X[9] + Na[4]*X[12] + Na[5]*X[15] +
	   Na[6]*X[18] + Na[7]*X[21] + Na[8]*X[24]);
  Xa[1] = (Na[0]*X[1] + Na[1]*X[4] + Na[2]*X[7] + 
	   Na[3]*X[10] + Na[4]*X[13] + Na[5]*X[16] +
	   Na[6]*X[19] + Na[7]*X[22] + Na[8]*X[25]);
  Xa[2] = (Na[0]*X[2] + Na[1]*X[5] + Na[2]*X[8] + 
	   Na[3]*X[11] + Na[4]*X[14] + Na[5]*X[17] +
	   Na[6]*X[20] + Na[7]*X[23] + Na[8]*X[26]);
}

/*
  Compute the inner product of a set of shape functions (or
  their derivative) with a 8*NUM_NODES array of nodes, variables
  etc.

  input:
  Na:   the shape function
  X:    the variables 

  ouput:
  Xa:   the inner product = Na^{T}*X
*/
static inline void innerProduct8( const double Na[],
				  const TacsScalar X[],
				  TacsScalar Xa[] ){
  Xa[0] = (Na[0]*X[0] + Na[1]*X[8] + Na[2]*X[16] + 
	   Na[3]*X[24] + Na[4]*X[32] + Na[5]*X[40] +
	   Na[6]*X[48] + Na[7]*X[56] + Na[8]*X[64]);
  Xa[1] = (Na[0]*X[1] + Na[1]*X[9] + Na[2]*X[17] + 
	   Na[3]*X[25] + Na[4]*X[33] + Na[5]*X[41] +
	   Na[6]*X[49] + Na[7]*X[57] + Na[8]*X[65]);
  Xa[2] = (Na[0]*X[2] + Na[1]*X[10] + Na[2]*X[18] + 
	   Na[3]*X[26] + Na[4]*X[34] + Na[5]*X[42] +
	   Na[6]*X[50] + Na[7]*X[58] + Na[8]*X[66]);
}

/*
  Compute the normal implied by the interpolation of the normal at
  each node to an arbitrary point in the element. This will not be a
  normal vector, in general, but is required to obtain strain-free
  rigid-body motion.

  input:
  N:    the shape functions for the element
  Xr:   the frames at each node

  output:
  fnormal:   the "frame normal" obtained by interpolating the 
  .          normal at each node
*/
static inline void computeFrameNormal( const double N[],
				       const TacsScalar Xr[],
				       TacsScalar fnormal[] ){
  fnormal[0] = (N[0]*Xr[2] + N[1]*Xr[11] + N[2]*Xr[20] +
		N[3]*Xr[29] + N[4]*Xr[38] + N[5]*Xr[47] +
		N[6]*Xr[56] + N[7]*Xr[65] + N[8]*Xr[74]);
  fnormal[1] = (N[0]*Xr[5] + N[1]*Xr[14] + N[2]*Xr[23] +
		N[3]*Xr[32] + N[4]*Xr[41] + N[5]*Xr[50] +
		N[6]*Xr[59] + N[7]*Xr[68] + N[8]*Xr[77]);
  fnormal[2] = (N[0]*Xr[8] + N[1]*Xr[17] + N[2]*Xr[26] +
		N[3]*Xr[35] + N[4]*Xr[44] + N[5]*Xr[53] +
		N[6]*Xr[62] + N[7]*Xr[71] + N[8]*Xr[80]);
}

/*
  Compute the derivative of the normal direction w.r.t. the parametric
  directions and multiply them

  This code computes the (approximate) through-thickness derivative of
  the inverse of the Jacobian matrix given as follows:

  d([X,r]^{-1})/dz

  The approximation enters through the fact that the normal direction
  is interpolated using the shape functions. This approach provides
  an approximation, but avoids the use of second derivatives. 

  input:
  Na:   the derivative of the shape functions
  Nb:   the derivative of the shape functions
*/
static inline void computeNormalRateMat( const double Na[],
					 const double Nb[],
					 const TacsScalar Xr[], 
					 const TacsScalar Xdinv[],
					 TacsScalar zXdinv[] ){
  // Compute the derivatives of the normal direction along the
  // parametric directions
  TacsScalar dn[9];

  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = (Na[0]*Xr[2] + Na[1]*Xr[11] + Na[2]*Xr[20] +
	   Na[3]*Xr[29] + Na[4]*Xr[38] + Na[5]*Xr[47] +
	   Na[6]*Xr[56] + Na[7]*Xr[65] + Na[8]*Xr[74]);
  dn[3] = (Na[0]*Xr[5] + Na[1]*Xr[14] + Na[2]*Xr[23] +
	   Na[3]*Xr[32] + Na[4]*Xr[41] + Na[5]*Xr[50] +
	   Na[6]*Xr[59] + Na[7]*Xr[68] + Na[8]*Xr[77]);
  dn[6] = (Na[0]*Xr[8] + Na[1]*Xr[17] + Na[2]*Xr[26] +
	   Na[3]*Xr[35] + Na[4]*Xr[44] + Na[5]*Xr[53] +
	   Na[6]*Xr[62] + Na[7]*Xr[71] + Na[8]*Xr[80]);
  
  dn[1] = (Nb[0]*Xr[2] + Nb[1]*Xr[11] + Nb[2]*Xr[20] +
	   Nb[3]*Xr[29] + Nb[4]*Xr[38] + Nb[5]*Xr[47] +
	   Nb[6]*Xr[56] + Nb[7]*Xr[65] + Nb[8]*Xr[74]);
  dn[4] = (Nb[0]*Xr[5] + Nb[1]*Xr[14] + Nb[2]*Xr[23] +
	   Nb[3]*Xr[32] + Nb[4]*Xr[41] + Nb[5]*Xr[50] +
	   Nb[6]*Xr[59] + Nb[7]*Xr[68] + Nb[8]*Xr[77]);
  dn[7] = (Nb[0]*Xr[8] + Nb[1]*Xr[17] + Nb[2]*Xr[26] +
	   Nb[3]*Xr[35] + Nb[4]*Xr[44] + Nb[5]*Xr[53] +
	   Nb[6]*Xr[62] + Nb[7]*Xr[71] + Nb[8]*Xr[80]);

  // Compute zXdinv = -Xdinv*dn*Xdinv
  TacsScalar tmp[9];
  matMatMult(dn, Xdinv, tmp);
  matMatMult(Xdinv, tmp, zXdinv);

  // Scale all the entries in zXdinv by -1
  zXdinv[0] *= -1.0;
  zXdinv[1] *= -1.0;
  zXdinv[2] *= -1.0;
  zXdinv[3] *= -1.0;
  zXdinv[4] *= -1.0;
  zXdinv[5] *= -1.0;
  zXdinv[6] *= -1.0;
  zXdinv[7] *= -1.0;
  zXdinv[8] *= -1.0;
}

/*
  Compute the tying shape functions for the element for the
  out-of-plane shear components.
  
  This uses the quadrature points from a two-point Gauss quadrature
  scheme as the tying points. The tying point shape functions for the
  g13 and g23 shear strains are numbered as follows:

  g13               g23
  +--4-----5--+     +-----------+
  |           |     3     4     5
  |  2  +  3  |     |     +     |
  |           |     0     1     2
  +--0-----1--+     +-----------+

  input:
  u:    the first parametric coordinate
  v:    the second parametric coordinate
  
  output:
  N13:  the shape functions associated with the g13 shear strain
  N23:  the shape functions associated with the g23 shear strain
*/
static inline void computeTyingFunc( const double u,
				     const double v,
				     double N13[], 
				     double N23[] ){
  // The tying point offset
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;
  const double tinv = 1.0/t;
  const double sinv = 1.0/(s*s);

  // Compute the shape functions
  double nu[3], nv[3];
  nu[0] = 0.5*sinv*u*(u - s);
  nu[1] = sinv*(s - u)*(s + u);
  nu[2] = 0.5*sinv*u*(s + u);

  nv[0] = 0.5*sinv*v*(v - s);
  nv[1] = sinv*(s - v)*(s + v);
  nv[2] = 0.5*sinv*v*(s + v);

  // Compute the shape functions for the reduced dimension
  double ntu[2], ntv[2];
  ntu[0] = 0.5*tinv*(t - u);
  ntu[1] = 0.5*tinv*(t + u);

  ntv[0] = 0.5*tinv*(t - v);
  ntv[1] = 0.5*tinv*(t + v);

  // Compute the shape functions for g13
  N13[0] = ntu[0]*nv[0];
  N13[1] = ntu[1]*nv[0];
  N13[2] = ntu[0]*nv[1];
  N13[3] = ntu[1]*nv[1];
  N13[4] = ntu[0]*nv[2];
  N13[5] = ntu[1]*nv[2];
  
  // Compute the shape functions for g23
  N23[0] = nu[0]*ntv[0];
  N23[1] = nu[1]*ntv[0];
  N23[2] = nu[2]*ntv[0];
  N23[3] = nu[0]*ntv[1];
  N23[4] = nu[1]*ntv[1];
  N23[5] = nu[2]*ntv[1];
}

/*
  Given the input vectors that are the derivative of the position
  vector along the u/v directions, and the orthogonal normal vector,
  assemble the results into a reference frame.  Note that this is not
  an orthonormal frame.

  input:
  Xa:      the first in-plane direction
  Xb:      the second in-plane direction
  normal:  the normal direction perpendicular to Xa and Xb
  
  output:
  Xr:      the reference frame associated with the given point
*/
static inline void assembleFrame( const TacsScalar Xa[],
				  const TacsScalar Xb[],
				  const TacsScalar normal[],
				  TacsScalar Xr[] ){
  // Fill in the values of Xr
  Xr[0] = Xa[0];
  Xr[3] = Xa[1];
  Xr[6] = Xa[2];
  
  Xr[1] = Xb[0];
  Xr[4] = Xb[1];
  Xr[7] = Xb[2];
  
  // Add the values for the normal
  Xr[2] = normal[0];
  Xr[5] = normal[1];
  Xr[8] = normal[2];
}

/*
  Compute the 3x4 matrix from the following matrix-matrix product:

  A = J*S = 2(I - n*n^{T})*[ -eps | (eta*I - eps^{x}) ] 

  Note that this is the product of the inertia contribution from the
  director and the angular rate kinematic S matrix.

  input:
  n:     the surface normal vector
  eta:   the quaternion scalar
  eps:   the quaternion vector
  
  output:
  A:     the product of the matrix (I - n*n^{T}) with S
*/
static inline void computeNormalRateProduct( const TacsScalar n[], 
					     const TacsScalar eta,
					     const TacsScalar eps[],
					     TacsScalar A[] ){
  // Compute the rate matrix
  TacsScalar S[12];
  computeSRateMat(eta, eps, S);
  
  // Pre-multiply the S matrix by (I - n*n^{T})
  TacsScalar ns = 0.0;

  // Compute the first column of A
  ns = n[0]*S[0] + n[1]*S[4] + n[2]*S[8];
  A[0] = S[0] - n[0]*ns;
  A[4] = S[4] - n[1]*ns;
  A[8] = S[8] - n[2]*ns;

  // Compute the second column of A
  ns = n[0]*S[1] + n[1]*S[5] + n[2]*S[9];
  A[1] = S[1] - n[0]*ns;
  A[5] = S[5] - n[1]*ns;
  A[9] = S[9] - n[2]*ns;

  // Compute the third column of A
  ns = n[0]*S[2] + n[1]*S[6] + n[2]*S[10];
  A[2] = S[2] - n[0]*ns;
  A[6] = S[6] - n[1]*ns;
  A[10] = S[10] - n[2]*ns;

  // Compute the third column of A
  ns = n[0]*S[3] + n[1]*S[7] + n[2]*S[11];
  A[3] = S[3] - n[0]*ns;
  A[7] = S[7] - n[1]*ns;
  A[11] = S[11] - n[2]*ns;
}

/*
  Constructor for the MITC9 element class

  input:
  stiff:      the stiffness object
  gravity:    the gravity vector
  vInit:      the initial velocity
  omegaInit:  the initial angular velocity
*/
MITC9::MITC9( FSDTStiffness *_stiff, TACSGibbsVector *_gravity,
              TACSGibbsVector *_vInit,
              TACSGibbsVector *_omegaInit ){
  // Set the stiffness 
  stiff = _stiff;
  stiff->incref();

  // Copy over the vectors (if they are defined)
  gravity = _gravity;
  vInit = _vInit;
  omegaInit = _omegaInit;
  if (gravity){ gravity->incref(); }
  if (vInit){ vInit->incref(); }
  if (omegaInit){ omegaInit->incref(); }

  // Get the Gauss quadrature points and weights
  FElibrary::getGaussPtsWts(ORDER, &gaussPts, &gaussWts);
}

MITC9::~MITC9(){
  stiff->decref();
  if (gravity){ gravity->decref(); }
  if (vInit){ vInit->decref(); }
  if (omegaInit){ omegaInit->decref(); }
}

/*
  Retrieve the initial values of the design variables
*/
void MITC9::getInitCondition( TacsScalar vars[], 
			      TacsScalar dvars[],
			      const TacsScalar X[] ){
  memset(vars, 0, 8*NUM_NODES*sizeof(TacsScalar));
  memset(dvars, 0, 8*NUM_NODES*sizeof(TacsScalar));

  // The initial quaternions are eta = 1.0, eps = 0
  for ( int i = 0; i < NUM_NODES; i++ ){
    vars[8*i + 3] = 1.0;
  }

  // If the initial velocity is defined
  if (vInit){
    const TacsScalar *v0;
    vInit->getVector(&v0);
    for ( int i = 0; i < NUM_NODES; i++ ){
      dvars[8*i] = v0[0];
      dvars[8*i+1] = v0[1];
      dvars[8*i+2] = v0[2];
    }
  }
  
  // If the initial angular velocity is defined
  if (omegaInit){
    const TacsScalar *omega;
    omegaInit->getVector(&omega);
    
    for ( int i = 0; i < NUM_NODES; i++ ){
      // dot{u} = v + r^{x}*omega
      crossProductAdd(1.0, omega, &X[3*i], &dvars[8*i]);
      
      // d{eps}/dt = 0.5*omega
      dvars[8*i+4] = 0.5*omega[0];
      dvars[8*i+5] = 0.5*omega[1];
      dvars[8*i+6] = 0.5*omega[2];
    }
  }
}

/*
  The following function evaluates the kinetic energy and potential
  and elastic energies of the element.

  These can be used to verify that the equations of motion are
  implemented correctly, since the element implements a method based
  on Lagrange's equations of motion.
  
  input:
  time:   the simulation time
  vars:   the values of the variables
  dvars:  the time-derivative of the variables

  output:
  Te:     the kinetic energy
  Pe:     the potential energy
*/
void MITC9::computeEnergies( double time,
                             TacsScalar *_Te, 
			     TacsScalar *_Pe,
			     const TacsScalar X[],
			     const TacsScalar vars[],
			     const TacsScalar dvars[] ){
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity){
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];  g[1] = _g[1];  g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity at the nodes
  TacsScalar omega[3*NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);

  // Compute the directors at the nodes
  TacsScalar dir[3*NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Initialize the velocities
  TacsScalar Te = 0.0, Pe = 0.0;

  // Evaluate the kinetic energy of the element
  for ( int j = 0; j < ORDER; j++ ){
    for ( int i = 0; i < ORDER; i++ ){
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);
      
      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9], Xdinv[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i]*gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar v0[3];
      innerProduct8(N, dvars, v0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3];
      innerProduct(N, omega, omeg);

      // Compute the dot product with the normal
      TacsScalar omegn = vecDot(omeg, fn);

      // Add the contributions to the kinetic energy
      Te += 0.5*h*(rho[0]*vecDot(v0, v0) + 
		   rho[1]*(vecDot(omeg, omeg) - omegn*omegn));

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar rot = computeRotPenalty(N, Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9]; 

      // Compute the cross product to find the normal direction
      TacsScalar normal[3];
      crossProduct(1.0, Xa, Xb, normal);
      TacsScalar nrm = sqrt(vecDot(normal, normal));
      vecScale(1.0/nrm, normal);

      // Scale the Xa direction so that it is a unit vector
      nrm = sqrt(vecDot(Xa, Xa));
      vecScale(1.0/nrm, Xa);

      // Compute the second perpendicular direction 
      crossProduct(1.0, normal, Xa, Xb);
   
      // Assemble the transformation matrix
      assembleFrame(Xa, Xb, normal, T);

      // Compute the displacement-based strain
      TacsScalar e[8], s[8];
      evalStrain(e, Ur, dr, Xdinv, zXdinv, T);

      // Add the contribution from the tying strain
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);

      // Compute the stress based on the strain values
      TacsScalar A[6], Bc[6], D[6], As[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);
      stiff->calculateStress(A, Bc, D, As, e, s);

      // Compute the terms for the potential energy due to gravity
      TacsScalar U[3];
      innerProduct8(N, vars, U);

      // Add the product of the stress times strain to the
      // potential energy computation
      Pe += 0.5*h*(strainProduct(s, e) + 
		   kpenalty*rot*rot - 
		   2.0*rho[0]*vecDot(g, U));
    }
  }

  *_Te = Te;
  *_Pe = Pe;
}

/*
  Get the residuals of the equations of motion.

  The equations of motion for this element are derived using
  Lagrange's equations with constraints. Each node imposes a
  constraint that its own quaternions satisfy the required unit norm
  constraint.

  The equations of motion can be divided into two parts: (1) the
  motion of the deformed surface and (2) the motion of the normals of
  the surface (or directors since they are no longer normal to the
  deformed surface during deformation). The linear motion takes the
  form:

  M*ddot{u} + dU/dx - fg = 0

  where M is a mass matrix. The rotational degrees of freedom satisfy
  the following equations of motion:

  S^{T}*J*d{omega} + 2*dot{S}^{T}*J*omega + dU/dq - A^{T}*lamb = 0

  where J = (I - n*n^{T}) is a rotational inertia term and the
  constraints produce the term A^{T}*lamb. 

  input:
  X:      the element nodal coordinates
  vars:   the element degrees of freedom
  dvars:  the first time derivative of the degrees of freedom
  ddvars: the second time derivative of the degrees of freedom

  output:
  res:    the residuals
*/
void MITC9::addResidual( double time,
                         TacsScalar res[],
			 const TacsScalar X[],
			 const TacsScalar vars[],
			 const TacsScalar dvars[],
			 const TacsScalar ddvars[] ){
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity){
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];  g[1] = _g[1];  g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3*NUM_NODES], domega[3*NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3*NUM_NODES], dirdq[12*NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  for ( int j = 0; j < ORDER; j++ ){
    for ( int i = 0; i < ORDER; i++ ){
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];
      
      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar Xdinv[9];
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i]*gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar a0[3];
      innerProduct8(N, ddvars, a0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp*fn[0];
      w[1] = omeg[1] - tmp*fn[1];
      w[2] = omeg[2] - tmp*fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp*fn[0];
      dw[1] = domeg[1] - tmp*fn[1];
      dw[2] = domeg[2] - tmp*fn[2];

      // Add the contribution to the residual
      TacsScalar *r = res;
      const TacsScalar *q = vars, *dq = dvars;
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
	// Add the contributions from the rectilinear velocity
	r[0] += h*N[ii]*rho[0]*a0[0];
	r[1] += h*N[ii]*rho[0]*a0[1];
	r[2] += h*N[ii]*rho[0]*a0[2];

	// Add the contributions from the angular velocity
	// S^{T}*dw + 2*dot{S}^{T}*w
	TacsScalar eta = q[3];
	const TacsScalar *eps = &q[4];
	TacsScalar deta = dq[3];
	const TacsScalar *deps = &dq[4];

	// Add S^{T}*dw
	r[3] -= 2.0*h*N[ii]*rho[1]*vecDot(eps, dw);
	crossProductAdd(2.0*h*N[ii]*rho[1], eps, dw, &r[4]);
	vecAxpy(2.0*h*N[ii]*eta*rho[1], dw, &r[4]);

	// Add 2*dot{S}^{T}*w
	r[3] -= 4.0*h*N[ii]*rho[1]*vecDot(deps, w);
	crossProductAdd(4.0*h*N[ii]*rho[1], deps, w, &r[4]);
	vecAxpy(4.0*h*N[ii]*deta*rho[1], w, &r[4]);

	r += 8;
	q += 8;
	dq += 8;
      }

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar Brot[8*NUM_NODES];
      TacsScalar rot = computeBRotPenalty(Brot, N, Na, Nb,
					  Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9]; 

      // Compute the cross product to find the normal direction
      TacsScalar normal[3];
      crossProduct(1.0, Xa, Xb, normal);
      TacsScalar nrm = sqrt(vecDot(normal, normal));
      vecScale(1.0/nrm, normal);

      // Scale the Xa direction so that it is a unit vector
      nrm = sqrt(vecDot(Xa, Xa));
      vecScale(1.0/nrm, Xa);
      // Compute the second perpendicular direction 
      crossProduct(1.0, normal, Xa, Xb);
      assembleFrame(Xa, Xb, normal, T);

      // Compute the displacement-based strain
      TacsScalar e[8], B[64*NUM_NODES];
      evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

      // Add the contribution from the tying straint
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
      addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

      // Compute the stress based on the strain values
      TacsScalar s[8]; // The stress components
      TacsScalar A[6], Bc[6], D[6], As[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);
      stiff->calculateStress(A, Bc, D, As, e, s);

      // Scale the rotation by the in-plane penalty term
      rot *= kpenalty;

      // Add the contribution to the residual
      r = res;
      const TacsScalar *b = B, *br = Brot;
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
	r[0] += h*(strainProduct(s, &b[0]) + rot*br[0] - rho[0]*N[ii]*g[0]);
	r[1] += h*(strainProduct(s, &b[8]) + rot*br[1] - rho[0]*N[ii]*g[1]);
	r[2] += h*(strainProduct(s, &b[16]) + rot*br[2] - rho[0]*N[ii]*g[2]);
	r[3] += h*(strainProduct(s, &b[24]) + rot*br[3]);
	r[4] += h*(strainProduct(s, &b[32]) + rot*br[4]);
	r[5] += h*(strainProduct(s, &b[40]) + rot*br[5]);
	r[6] += h*(strainProduct(s, &b[48]) + rot*br[6]);
	r[7] += h*(strainProduct(s, &b[56]) + rot*br[7]);

	r += 8;
	b += 64;
	br += 8;
      }
    }
  }

  // Add the constraints from the quaternion parametrization
  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i+3];
    TacsScalar lamb = vars[8*i+7];

    // Enforce the quaternion constraint
    res[8*i+7] += q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] - 1.0;

    // Add the result to the governing equations
    res[8*i+3] += 2.0*q[0]*lamb;
    res[8*i+4] += 2.0*q[1]*lamb;
    res[8*i+5] += 2.0*q[2]*lamb;
    res[8*i+6] += 2.0*q[3]*lamb;
  }
}

/*
  Add the Jacobian for the governing equations

  The following code computes the Jacobian of the equations of motion
  involving a linear combination of the Jacobian w.r.t. the variables
  and their first and second time derivatives as follows:

  J = alpha*dR/dq + beta*dR/d(dot{q}) + gamma*dR/d(ddot{q})

  The alpha, beta, and gamma coefficients are used for the values and
  their first and second time-derivatives respectively.  This function
  is used to assemble the element-contribution to the multibody system
  Jacobian.

  input:
  alpha:   the variable coefficient
  beta:    the first-derivative coefficient
  gamma:   the second-derivative coefficient
  X:       the nodal locations (of the initial configuration)
  vars:    the variable values
  dvars:   the time derivative of the variable values
  ddvars:  the second time derivative of the variable values

  output:
  J:       the Jacobian matrix
*/
void MITC9::addJacobian( double time, TacsScalar J[],
			 double alpha, double beta, double gamma,
			 const TacsScalar X[],
			 const TacsScalar vars[],
			 const TacsScalar dvars[],
			 const TacsScalar ddvars[] ){
  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3*NUM_NODES], domega[3*NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3*NUM_NODES], dirdq[12*NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // The weights that are used for the geometric stiffness from
  // the tying strain
  TacsScalar w13[6], w23[6];
  memset(w13, 0, 6*sizeof(TacsScalar));
  memset(w23, 0, 6*sizeof(TacsScalar));

  for ( int j = 0; j < ORDER; j++ ){
    for ( int i = 0; i < ORDER; i++ ){
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];
 
      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9], Xdinv[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i]*gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // Add the contributions from the linear motion
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
	for ( int jj = 0; jj < NUM_NODES; jj++ ){
	  const TacsScalar scale = gamma*h*N[ii]*N[jj]*rho[0];
	  // Add the contributions from the rectilinear velocity
	  J[8*NUM_NODES*(8*ii) + 8*jj] += scale;
	  J[8*NUM_NODES*(8*ii+1) + 8*jj+1] += scale;
	  J[8*NUM_NODES*(8*ii+2) + 8*jj+2] += scale;
	}
      }

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp*fn[0];
      w[1] = omeg[1] - tmp*fn[1];
      w[2] = omeg[2] - tmp*fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp*fn[0];
      dw[1] = domeg[1] - tmp*fn[1];
      dw[2] = domeg[2] - tmp*fn[2];

      // Add the contributions from the rotational DOF
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
	// Compute and store the product of (I - n*n^{T}) with the
	// angular rate matrices
	TacsScalar JSii[12], JdSii[12];
	computeNormalRateProduct(fn, vars[8*ii+3], 
				 &vars[8*ii+4], JSii);
	computeNormalRateProduct(fn, dvars[8*ii+3], 
				 &dvars[8*ii+4], JdSii);

	// Set the pointer to the Jacobian entries that will
	// be added
	TacsScalar *Jp = &J[(8*NUM_NODES+1)*(8*ii+3)];
	const int ldj = 8*NUM_NODES;

	// Add the diagonal terms
	const TacsScalar dscale = h*N[ii]*rho[1];
	addSRateMatTransDeriv(alpha*dscale, dw, Jp, ldj);
	addSRateMatTransDeriv(2.0*beta*dscale, w, Jp, ldj);

	// Add the result to the Jacobian matrix
	const TacsScalar *q = vars, *dq = dvars, *ddq = ddvars;
	for ( int jj = 0; jj < NUM_NODES; jj++ ){
	  // Set the common scaling factor for all terms
	  const TacsScalar scale = h*N[ii]*N[jj]*rho[1];

	  // Set the pointer to the Jacobian entries that will
	  // be added
	  TacsScalar *Jp = &J[8*NUM_NODES*(8*ii+3) + 8*jj+3];
	  
	  // Compute S = S(q)
	  TacsScalar Sjj[12];
	  computeSRateMat(q[3], &q[4], Sjj);

	  // Compute dot{S} = S(dot{q})
	  TacsScalar dSjj[12];
	  computeSRateMat(dq[3], &dq[4], dSjj);

	  // Compute ddot{S} = S(ddot{q})
	  TacsScalar ddSjj[12];
	  computeSRateMat(ddq[3], &ddq[4], ddSjj);

	  // Add the Jacobian terms from the DOF:
	  // T(dw) - S^{T}*J*S(ddot{q}) - 2*dot{S}^{T}*J*dot{S}
	  addBlock3x4Product(-alpha*scale, JSii, ddSjj, Jp, ldj);
	  addBlock3x4Product(-2.0*alpha*scale, JdSii, dSjj, Jp, ldj);

	  // Add the Jacobian terms from the first time derivatives:
	  // 2*dot{S}^{T}*J*S
	  addBlock3x4Product(2.0*beta*scale, JdSii, Sjj, Jp, ldj);

	  // Add the Jacobian terms from the second time derivatives:
	  // S^{T}*J*S
	  addBlock3x4Product(gamma*scale, JSii, Sjj, Jp, ldj);

	  q += 8;
	  dq += 8;
	  ddq += 8;
	}
      }

      // This code is expensive, and is only required when alpha 
      // is non-zero
      if (alpha != 0.0){
	// Everything left in this function is proportional to
	// h*alpha, so we scale h by alpha
	h *= alpha;

        // Evaluate the stiffness terms at the parametric point
        TacsScalar A[6], Bc[6], D[6], As[3];
        TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);

	// Evaluate the tying strain interpolation
	double N13[6], N23[6];
	computeTyingFunc(u, v, N13, N23);

	// Compute the through-thickness derivative of [X,r]^{-1}
	TacsScalar zXdinv[9];
	computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

	// Compute the derivatives of Ua/Ub along the given directions
	TacsScalar Ur[9], Ua[3], Ub[3], d[3];
	innerProduct8(Na, vars, Ua);
	innerProduct8(Nb, vars, Ub);
	innerProduct(N, dir, d);
	assembleFrame(Ua, Ub, d, Ur);
	
	// Now compute the derivatives of the director along each
	// coordinate direction
	TacsScalar dr[9], da[3], db[3], zero[3];
	innerProduct(Na, dir, da);
	innerProduct(Nb, dir, db);
	zero[0] = zero[1] = zero[2] = 0.0;
	assembleFrame(da, db, zero, dr);

	// Add the in-plane penalty terms
	TacsScalar Brot[8*NUM_NODES];
	TacsScalar rot = computeBRotPenalty(Brot, N, Na, Nb,
					    Xa, Xb, Ua, Ub, vars);

	// Add the contribution from the penalty
	addGRotMat(J, h*kpenalty*rot, N, Na, Nb, 
		   Xa, Xb, Ua, Ub, vars);
	
	// Compute the transformation to the locally-aligned frame
	TacsScalar T[9]; 
	
	// Compute the cross product to find the normal direction
	TacsScalar normal[3];
	crossProduct(1.0, Xa, Xb, normal);
	TacsScalar nrm = sqrt(vecDot(normal, normal));
	vecScale(1.0/nrm, normal);

	// Scale the Xa direction so that it is a unit vector
	nrm = sqrt(vecDot(Xa, Xa));
	vecScale(1.0/nrm, Xa);

	// Compute the second perpendicular direction 
	crossProduct(1.0, normal, Xa, Xb);
	
	// Assemble the transformation matrix
	assembleFrame(Xa, Xb, normal, T);
	
	// Compute the displacement-based strain
	TacsScalar e[8], B[64*NUM_NODES];
	evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

	// Add the contribution from the tying straint
	addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
	addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

        // Compute the stress based on the strain values
        TacsScalar s[8]; // The stress components
        stiff->calculateStress(A, Bc, D, As, e, s);

	// Add to the weights
	addTyingGmatWeights(w13, w23, h, s,
			    N13, N23, Xdinv, T);

	// Add the geometric stiffness terms
	addGmat(J, h, s, N, Na, Nb,
		Ur, dr, Xdinv, zXdinv, T, Xr, dirdq);
	
	// Add the contribution to the residual
	for ( int ii = 0; ii < NUM_NODES; ii++ ){
	  for ( int ik = 0; ik < 7; ik++ ){
	    // Compute the stress from the 8*i + ik component
	    TacsScalar sbii[8];
            stiff->calculateStress(A, Bc, D, As, 
                                   &B[8*(8*ii + ik)], sbii);

	    // Compute the 
	    TacsScalar pr = kpenalty*Brot[8*ii + ik];

	    const TacsScalar *b = B, *br = Brot;
	    TacsScalar *Jp = &J[8*NUM_NODES*(8*ii + ik)];
	    for ( int jj = 0; jj < NUM_NODES; jj++ ){
	      Jp[0] += h*(strainProduct(sbii, &b[0]) + pr*br[0]);
	      Jp[1] += h*(strainProduct(sbii, &b[8]) + pr*br[1]);
	      Jp[2] += h*(strainProduct(sbii, &b[16]) + pr*br[2]);
	      Jp[3] += h*(strainProduct(sbii, &b[24]) + pr*br[3]);
	      Jp[4] += h*(strainProduct(sbii, &b[32]) + pr*br[4]);
	      Jp[5] += h*(strainProduct(sbii, &b[40]) + pr*br[5]);
	      Jp[6] += h*(strainProduct(sbii, &b[48]) + pr*br[6]);
	  
	      Jp += 8;
	      b += 64;
	      br += 8;
	    }
	  }
	}
      }
    }
  }

  // Add the geometric stiffness terms from the tying strain
  addTyingGmat(J, w13, w23, X, Xr, vars, dir, dirdq);

  // Add the constraints from the quaternions
  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i+3];
    const int ldj = 8*NUM_NODES;

    TacsScalar *Jp = &J[(8*NUM_NODES+1)*(8*i + 3)];
    TacsScalar lamb = vars[8*i+7];

    // Add the constraint terms
    Jp[4] += 2.0*alpha*q[0];
    Jp[4+ldj] += 2.0*alpha*q[1];
    Jp[4+2*ldj] += 2.0*alpha*q[2];
    Jp[4+3*ldj] += 2.0*alpha*q[3];

    // Enforce the quaternion constraint
    Jp[4*ldj] += 2.0*alpha*q[0];
    Jp[4*ldj+1] += 2.0*alpha*q[1];
    Jp[4*ldj+2] += 2.0*alpha*q[2];
    Jp[4*ldj+3] += 2.0*alpha*q[3];

    // Add the terms to the diagonal
    Jp[0] += 2.0*alpha*lamb;
    Jp[ldj+1] += 2.0*alpha*lamb;
    Jp[2*(ldj+1)] += 2.0*alpha*lamb;
    Jp[3*(ldj+1)] += 2.0*alpha*lamb;
  }
}

/*
  Add the derivative of the product of the adjoint variables with the
  residuals to the given design variable vector.

  input:
  time:     the simulation time
  scale:    the scalar factor applied to the result
  fdvSens:  the derivative (times scale) is accumulated in this array
  dvLen:    the design array length
  psi:      the adjoint vector
  X:        the node locations
  vars:     the variable values
  dvars:    the time derivatives of the variable values
  ddvars:   the second time derivative of the variable values
*/
void MITC9::addAdjResProduct( double time, double scale,
                              TacsScalar fdvSens[], int dvLen,
                              const TacsScalar psi[],
                              const TacsScalar X[],
                              const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[] ){
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity){
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];  g[1] = _g[1];  g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3*NUM_NODES], domega[3*NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3*NUM_NODES], dirdq[12*NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  for ( int j = 0; j < ORDER; j++ ){
    for ( int i = 0; i < ORDER; i++ ){
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];
      
      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar Xdinv[9];
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= scale*gaussWts[i]*gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar a0[3];
      innerProduct8(N, ddvars, a0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp*fn[0];
      w[1] = omeg[1] - tmp*fn[1];
      w[2] = omeg[2] - tmp*fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp*fn[0];
      dw[1] = domeg[1] - tmp*fn[1];
      dw[2] = domeg[2] - tmp*fn[2];

      // Add the contribution to the residual
      TacsScalar mscale[2] = {0.0, 0.0};
      const TacsScalar *p = psi;
      const TacsScalar *q = vars, *dq = dvars;
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
	// Add the contributions from the rectilinear velocity
        mscale[0] += h*N[ii]*(a0[0]*p[0] + a0[1]*p[1] + a0[2]*p[2]);

	// Add the contributions from the angular velocity
	// S^{T}*dw + 2*dot{S}^{T}*w
	TacsScalar eta = q[3];
	const TacsScalar *eps = &q[4];
	TacsScalar deta = dq[3];
	const TacsScalar *deps = &dq[4];

	// Add p^{T}*S^{T}*dw
        TacsScalar t[3];
        crossProduct(1.0, eps, dw, t);
	mscale[1] -= 2.0*p[3]*h*N[ii]*vecDot(eps, dw);
        mscale[1] += 2.0*h*N[ii]*(eta*vecDot(dw, &p[4]) + vecDot(t, &p[4]));

	// Add p^{T}*2*dot{S}^{T}*w
        crossProduct(1.0, deps, w, t);
	mscale[1] -= 4.0*p[3]*h*N[ii]*vecDot(deps, w);
        mscale[1] += 4.0*h*N[ii]*(deta*vecDot(w, &p[4]) + vecDot(t, &p[4]));

        // Increment the pointers to the variables/multipliers
	q += 8;
	dq += 8;
        p += 8;
      }

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar Brot[8*NUM_NODES];
      TacsScalar rot = computeBRotPenalty(Brot, N, Na, Nb,
					  Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9]; 

      // Compute the cross product to find the normal direction
      TacsScalar normal[3];
      crossProduct(1.0, Xa, Xb, normal);
      TacsScalar nrm = sqrt(vecDot(normal, normal));
      vecScale(1.0/nrm, normal);

      // Scale the Xa direction so that it is a unit vector
      nrm = sqrt(vecDot(Xa, Xa));
      vecScale(1.0/nrm, Xa);
      // Compute the second perpendicular direction 
      crossProduct(1.0, normal, Xa, Xb);
      assembleFrame(Xa, Xb, normal, T);

      // Compute the displacement-based strain
      TacsScalar e[8], B[64*NUM_NODES];
      evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

      // Add the contribution from the tying strain
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
      addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

      // Set the Lagrange multiplier associated with the strain
      TacsScalar epsi[8];
      TacsScalar rotPsi = 0.0;
      epsi[0] = epsi[1] = epsi[2] = epsi[3] =
        epsi[4] = epsi[5] = epsi[6] = epsi[7] = 0.0;

      const TacsScalar *b = B, *br = Brot;
      for ( int ii = 0; ii < 8*NUM_NODES; ii++ ){
        epsi[0] += b[0]*psi[ii];
        epsi[1] += b[1]*psi[ii];
        epsi[2] += b[2]*psi[ii];
        epsi[3] += b[3]*psi[ii];
        epsi[4] += b[4]*psi[ii];
        epsi[5] += b[5]*psi[ii];
        epsi[6] += b[6]*psi[ii];
        epsi[7] += b[7]*psi[ii];
        rotPsi += br[0]*psi[ii];
        b += 8;
        br++;
      }

      // Add the contribution from the gravity load
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
        mscale[0] -= 
          h*N[ii]*(g[0]*psi[8*ii] + g[1]*psi[8*ii+1] + g[2]*psi[8*ii+2]);
      }

      // Scale the psi vector by the determinant of the Jacobian
      // transformation
      for ( int k = 0; k < 8; k++ ){
        epsi[k] *= h;
      }

      // Add the derivative contribution from the mass/area
      stiff->addPointwiseMassDVSens(pt, mscale, fdvSens, dvLen);

      // Add the derivative
      stiff->addStiffnessDVSens(pt, e, epsi, h*rot*rotPsi,
                                fdvSens, dvLen);
    }
  }
}

/*
  Given the nodal degrees of freedom and their time-derivatives,
  compute the angular velocity at each node
*/
void MITC9::computeAngularVelocity( TacsScalar omega[],
				    const TacsScalar vars[],
				    const TacsScalar dvars[] ){
  for ( int i = 0; i < NUM_NODES; i++ ){
    TacsScalar eta = vars[3];
    const TacsScalar *eps = &vars[4];
    TacsScalar deta = dvars[3];
    const TacsScalar *deps = &dvars[4];

    // omega = -2*eps^{x}*deps + 2*eta*deps - eps*deta
    crossProduct(-2.0, eps, deps, omega);
    vecAxpy(2.0*eta, deps, omega);
    vecAxpy(-2.0*deta, eps, omega);

    omega += 3;
    vars += 8;
    dvars += 8;
  }
}

/*
  Given the nodal degrees of freedom and their first and second
  time-derivatives, compute the angular acceleration at the nodes.
*/
void MITC9::computeAngularAccel( TacsScalar domega[],
				 const TacsScalar vars[],
				 const TacsScalar dvars[],
				 const TacsScalar ddvars[] ){
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set pointers to the values 
    TacsScalar eta = vars[3];
    const TacsScalar *eps = &vars[4];
    TacsScalar deta = dvars[3];
    const TacsScalar *deps = &dvars[4];
    TacsScalar ddeta = ddvars[3];
    const TacsScalar *ddeps = &ddvars[4];

    // domega = S(q)*ddot{q}
    crossProduct(-2.0, eps, ddeps, domega);
    vecAxpy(2.0*eta, ddeps, domega);
    vecAxpy(-2.0*ddeta, eps, domega);

    domega += 3;
    vars += 8;
    dvars += 8;
    ddvars += 8;
  }
}

/*
  At each node in the finite-element, compute the derivatives of
  the coordinates directions and assemble a locally-aligned reference
  frame.

  Each locally-aligned reference frame Xr[] consists of a 3x3 matrix
  stored in row-major order where the first direction is aligned
  along the 1-direction and the third direction is normal to the 
  local surface and the second direction is perpendicular to both
  these directions. This can be written as follows:

  Xr = [X,xi1; X,xi2; n]

  input:
  X:    the initial nodal locations
  
  output:
  Xr:   the locally-aligned frames
*/
void MITC9::computeFrames( TacsScalar Xr[],
			   const TacsScalar X[] ){
  for ( int j = 0; j < ORDER; j++ ){
    for ( int i = 0; i < ORDER; i++ ){
      // Find the u/v values at the node locations
      double u = -1.0 + 2.0*i/(ORDER-1.0);
      double v = -1.0 + 2.0*j/(ORDER-1.0);

      // Evaluate the shape functions
      double Na[9], Nb[9];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the cross product
      TacsScalar normal[3];
      crossProduct(1.0, Xa, Xb, normal);
      TacsScalar nrm = sqrt(vecDot(normal, normal));
      vecScale(1.0/nrm, normal);

      // Assemble the frame into the Xr matrix
      assembleFrame(Xa, Xb, normal, Xr);

      // Increment the pointer to the frames
      Xr += 9;
    }
  }
}

/*
  Based on the values of the element variables, determine the director
  values.
  
  The director is the difference between the deformed normal direction
  and the initial normal. The directors are used to determine the
  through-thickness strain distribution in the finite-element model.
  The advantage of the director approach is that it enables a
  geometrically exact displacement representation and also allows for
  shell-shell intersections at 90 degree angles.

  input:
  vars:   the element variable values
  Xr:     the local frames at each node in the mesh

  output:
  d:      the director at each node in the finite-element
*/
void MITC9::computeDirectors( TacsScalar d[],
			      const TacsScalar vars[],
			      const TacsScalar Xr[] ){
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute C = rot - I
    TacsScalar C[9];
    C[0] =-2.0*(q[2]*q[2] + q[3]*q[3]);
    C[1] = 2.0*(q[1]*q[2] + q[3]*q[0]);
    C[2] = 2.0*(q[1]*q[3] - q[2]*q[0]);

    C[3] = 2.0*(q[2]*q[1] - q[3]*q[0]);
    C[4] =-2.0*(q[1]*q[1] + q[3]*q[3]);
    C[5] = 2.0*(q[2]*q[3] + q[1]*q[0]);

    C[6] = 2.0*(q[3]*q[1] + q[2]*q[0]);
    C[7] = 2.0*(q[3]*q[2] - q[1]*q[0]);
    C[8] =-2.0*(q[1]*q[1] + q[2]*q[2]);

    // Compute d = C^{T}*n
    d[0] = C[0]*Xr[2] + C[3]*Xr[5] + C[6]*Xr[8];
    d[1] = C[1]*Xr[2] + C[4]*Xr[5] + C[7]*Xr[8];
    d[2] = C[2]*Xr[2] + C[5]*Xr[5] + C[8]*Xr[8];

    d += 3; // Each director is a 3-vector
    Xr += 9; // Increment over each frame
    vars += 8; // 8 variables per node
  }
}

/*
  Compute the derivative of the director values w.r.t. the rotation
  matrix parametrization. 

  This code computes the derivative of (C^{T} - I)*n w.r.t. the
  quaternions at each node in the finite-element. This is required
  to compute the governing equations of motion for the element.

  input:
  vars:  the element variables
  Xr:    the local frame for each node in the mesh

  output:
  ddq:   the derivative of the directors w.r.t. the quaterions
*/
void MITC9::computeDirectorDeriv( TacsScalar ddq[],
				  const TacsScalar vars[],
				  const TacsScalar Xr[] ){
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute the derivative of the rotation matrix w.r.t.
    // the quaternions
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0*q[3];
    Q[2] =-2.0*q[2];

    Q[3] =-2.0*q[3];
    Q[4] = 0.0;
    Q[5] = 2.0*q[1];

    Q[6] = 2.0*q[2];
    Q[7] =-2.0*q[1];
    Q[8] = 0.0;

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0]*Xr[2] + Q[3]*Xr[5] + Q[6]*Xr[8];
    ddq[1] = Q[1]*Xr[2] + Q[4]*Xr[5] + Q[7]*Xr[8];
    ddq[2] = Q[2]*Xr[2] + Q[5]*Xr[5] + Q[8]*Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0*q[2];
    Q[2] = 2.0*q[3];

    Q[3] = 2.0*q[2];
    Q[4] =-4.0*q[1];
    Q[5] = 2.0*q[0];

    Q[6] = 2.0*q[3];
    Q[7] =-2.0*q[0];
    Q[8] =-4.0*q[1];

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0]*Xr[2] + Q[3]*Xr[5] + Q[6]*Xr[8];
    ddq[1] = Q[1]*Xr[2] + Q[4]*Xr[5] + Q[7]*Xr[8];
    ddq[2] = Q[2]*Xr[2] + Q[5]*Xr[5] + Q[8]*Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[2]
    Q[0] =-4.0*q[2];
    Q[1] = 2.0*q[1];
    Q[2] =-2.0*q[0];

    Q[3] = 2.0*q[1];
    Q[4] = 0.0;
    Q[5] = 2.0*q[3];

    Q[6] = 2.0*q[0];
    Q[7] = 2.0*q[3];
    Q[8] =-4.0*q[2];

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0]*Xr[2] + Q[3]*Xr[5] + Q[6]*Xr[8];
    ddq[1] = Q[1]*Xr[2] + Q[4]*Xr[5] + Q[7]*Xr[8];
    ddq[2] = Q[2]*Xr[2] + Q[5]*Xr[5] + Q[8]*Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[3]
    Q[0] =-4.0*q[3];
    Q[1] = 2.0*q[0];
    Q[2] = 2.0*q[1];

    Q[3] =-2.0*q[0];
    Q[4] =-4.0*q[3];
    Q[5] = 2.0*q[2];

    Q[6] = 2.0*q[1];
    Q[7] = 2.0*q[2];
    Q[8] = 0.0;

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0]*Xr[2] + Q[3]*Xr[5] + Q[6]*Xr[8];
    ddq[1] = Q[1]*Xr[2] + Q[4]*Xr[5] + Q[7]*Xr[8];
    ddq[2] = Q[2]*Xr[2] + Q[5]*Xr[5] + Q[8]*Xr[8];
    ddq += 3;

    Xr += 9; // Increment over each frame
    vars += 8; // 8 variables per node
  }
}  

/*
  Given the derivatives of the displacements and the transformations
  to the local coordinates, evaluate strain.

  The expressions for the strain involve the surface derivatives of
  the displacements along the coordinate directions and the derivative
  of the director along the surface coordinate directions:

  Ur = [ U0,xi; d(xi) ]
  dr = [ d,xi; 0 ]

  The geometry enters through both an surface-based transformation and
  a derivative of the transformation through the thickness of the
  shell as follows:

  Xdinv = [X,r]^{-1}
  zXdinv = -[X,r]^{-1}*[n,xi]*[X,r]^{-1}
  
  input:
  Ur:     the derivative of the displacements w.r.t. all shell coordinates
  dr:     the derivative of the director w.r.t. all shell coordinates
  Xdinv:  the inverse of the Jacobian matrix
  zXdinv: the derivative of the inverse of the Jacobian w.r.t. z

  output:
  e:      the displacement-based components of the strain 
*/
void MITC9::evalStrain( TacsScalar e[], 
			const TacsScalar Ur[], 
			const TacsScalar dr[],
			const TacsScalar Xdinv[],
			const TacsScalar zXdinv[],
			const TacsScalar T[] ){
  TacsScalar U0[9], U1[9], tmp[9];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Compute the in-plane strain
  e[0] = U0[0] + 0.5*(U0[0]*U0[0] + U0[3]*U0[3] + U0[6]*U0[6]);
  e[1] = U0[4] + 0.5*(U0[1]*U0[1] + U0[4]*U0[4] + U0[7]*U0[7]);
  e[2] = U0[1] + U0[3] + (U0[0]*U0[1] + U0[3]*U0[4] + U0[6]*U0[7]);

  // Compute the bending strain
  e[3] = U1[0] + (U0[0]*U1[0] + U0[3]*U1[3] + U0[6]*U1[6]);
  e[4] = U1[4] + (U0[1]*U1[1] + U0[4]*U1[4] + U0[7]*U1[7]);
  e[5] = U1[1] + U1[3] + (U0[0]*U1[1] + U0[3]*U1[4] + U0[6]*U1[7] +
			  U1[0]*U0[1] + U1[3]*U0[4] + U1[6]*U0[7]);
}

/*
  Evaluate the strain and the derivative of the strain w.r.t. the
  element variables. 

*/
void MITC9::evalBmat( TacsScalar e[],
		      TacsScalar B[],
		      const double N[],
		      const double Na[],
		      const double Nb[],
		      const TacsScalar Ur[],
		      const TacsScalar dr[],
		      const TacsScalar Xdinv[],
		      const TacsScalar zXdinv[],
		      const TacsScalar T[],
		      const TacsScalar dirdq[] ){
  TacsScalar U0[9], U1[9], tmp[9];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Compute the in-plane strain
  e[0] = U0[0] + 0.5*(U0[0]*U0[0] + U0[3]*U0[3] + U0[6]*U0[6]);
  e[1] = U0[4] + 0.5*(U0[1]*U0[1] + U0[4]*U0[4] + U0[7]*U0[7]);
  e[2] = U0[1] + U0[3] + (U0[0]*U0[1] + U0[3]*U0[4] + U0[6]*U0[7]);

  // Compute the bending strain
  e[3] = U1[0] + (U0[0]*U1[0] + U0[3]*U1[3] + U0[6]*U1[6]);
  e[4] = U1[4] + (U0[1]*U1[1] + U0[4]*U1[4] + U0[7]*U1[7]);
  e[5] = U1[1] + U1[3] + (U0[0]*U1[1] + U0[3]*U1[4] + U0[6]*U1[7] +
			  U1[0]*U0[1] + U1[3]*U0[4] + U1[6]*U0[7]);

  // Compute the derivatives of the strain w.r.t. the displacement
  // variables. This code takes advantage of the sparsity of the
  // derivatives to simplify the computations
  for ( int i = 0; i < NUM_NODES; i++ ){
    TacsScalar *b = &B[8*(8*i)];

    // Compute the values of the displacements
    TacsScalar dU0[9], dU1[9], ztmp[3];

    // Compute dU0 = d(Ur)/dq_{k} * Xdinv
    // [  0,  0, 0 ][ 0  1  2 ][ T0  T1  T2 ]
    // [ Na, Nb, 0 ][ 3  4  5 ][ T3  T4  T5 ]
    // [  0,  0, 0 ][ 6  7  8 ][ T6  T7  T8 ]
    dU0[0] = Na[i]*Xdinv[0] + Nb[i]*Xdinv[3];
    dU0[1] = Na[i]*Xdinv[1] + Nb[i]*Xdinv[4];
    dU0[2] = Na[i]*Xdinv[2] + Nb[i]*Xdinv[5];
    matMultTrans(T, dU0, tmp);

    dU1[0] = Na[i]*zXdinv[0] + Nb[i]*zXdinv[3];
    dU1[1] = Na[i]*zXdinv[1] + Nb[i]*zXdinv[4];
    dU1[2] = Na[i]*zXdinv[2] + Nb[i]*zXdinv[5];
    matMultTrans(T, dU1, ztmp);

    for ( int k = 0; k < 3; k++ ){
      // dU0 = T[3*k+i]*tmp[3*k+j];
      dU0[0] = T[3*k]*tmp[0];
      dU0[1] = T[3*k]*tmp[1];
      dU0[2] = T[3*k]*tmp[2];
      dU0[3] = T[3*k+1]*tmp[0];
      dU0[4] = T[3*k+1]*tmp[1];
      dU0[5] = T[3*k+1]*tmp[2];
      dU0[6] = T[3*k+2]*tmp[0];
      dU0[7] = T[3*k+2]*tmp[1];
      dU0[8] = T[3*k+2]*tmp[2];

      // dU1 = T[3*k+i]*ztmp[3*k+j];
      dU1[0] = T[3*k]*ztmp[0];
      dU1[1] = T[3*k]*ztmp[1];
      dU1[2] = T[3*k]*ztmp[2];
      dU1[3] = T[3*k+1]*ztmp[0];
      dU1[4] = T[3*k+1]*ztmp[1];
      dU1[5] = T[3*k+1]*ztmp[2];
      dU1[6] = T[3*k+2]*ztmp[0];
      dU1[7] = T[3*k+2]*ztmp[1];
      dU1[8] = T[3*k+2]*ztmp[2];
    
      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0]*dU0[0] + U0[3]*dU0[3] + U0[6]*dU0[6]);
      b[1] = dU0[4] + (U0[1]*dU0[1] + U0[4]*dU0[4] + U0[7]*dU0[7]);
      b[2] = dU0[1] + dU0[3] + (U0[0]*dU0[1] + U0[3]*dU0[4] + U0[6]*dU0[7] +
				dU0[0]*U0[1] + dU0[3]*U0[4] + dU0[6]*U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0]*dU1[0] + U0[3]*dU1[3] + U0[6]*dU1[6] + 
		       dU0[0]*U1[0] + dU0[3]*U1[3] + dU0[6]*U1[6]);
      b[4] = dU1[4] + (U0[1]*dU1[1] + U0[4]*dU1[4] + U0[7]*dU1[7] +
		       dU0[1]*U1[1] + dU0[4]*U1[4] + dU0[7]*U1[7]);
      b[5] = dU1[1] + dU1[3] + (U0[0]*dU1[1] + U0[3]*dU1[4] + U0[6]*dU1[7] +
				U1[0]*dU0[1] + U1[3]*dU0[4] + U1[6]*dU0[7] +
				dU0[0]*U1[1] + dU0[3]*U1[4] + dU0[6]*U1[7] +
				dU1[0]*U0[1] + dU1[3]*U0[4] + dU1[6]*U0[7]);
      b[6] = b[7] = 0.0;
      b += 8;
    }
  }
  
  // Add the contributions from the derivative of the strain w.r.t. the
  // rotation variables. These derivatives only make contributions to the
  // bending strains, not the in-plane strains.
  for ( int i = 0; i < NUM_NODES; i++ ){
    TacsScalar *b = &B[8*(8*i+3)];

    // Compute the values of the displacements
    TacsScalar dU0[9], dU1[9], drdq[9];

    for ( int k = 0; k < 4; k++ ){
      // Compute dU0 = T^{T}*dUr*Xdinv*T
      // T^{T}*[ 0 | 0 | N*dirdq[0] ]*Xdinv*T
      // .     [ 0 | 0 | N*dirdq[1] ]*Xdinv*T
      //       [ 0 | 0 | N*dirdq[2] ]*Xdinv*T
      dU0[0] = N[i]*dirdq[0]*Xdinv[6];
      dU0[1] = N[i]*dirdq[0]*Xdinv[7];
      dU0[2] = N[i]*dirdq[0]*Xdinv[8];
      dU0[3] = N[i]*dirdq[1]*Xdinv[6];
      dU0[4] = N[i]*dirdq[1]*Xdinv[7];
      dU0[5] = N[i]*dirdq[1]*Xdinv[8];
      dU0[6] = N[i]*dirdq[2]*Xdinv[6];
      dU0[7] = N[i]*dirdq[2]*Xdinv[7];
      dU0[8] = N[i]*dirdq[2]*Xdinv[8];
      matMatMult(dU0, T, tmp);
      matTransMatMult(T, tmp, dU0);

      // Compute the derivative for d
      drdq[0] = Na[i]*dirdq[0];
      drdq[3] = Na[i]*dirdq[1];
      drdq[6] = Na[i]*dirdq[2];
      drdq[1] = Nb[i]*dirdq[0];
      drdq[4] = Nb[i]*dirdq[1];
      drdq[7] = Nb[i]*dirdq[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, dU1);

      dU1[0] += N[i]*dirdq[0]*zXdinv[6];
      dU1[1] += N[i]*dirdq[0]*zXdinv[7];
      dU1[2] += N[i]*dirdq[0]*zXdinv[8];
      dU1[3] += N[i]*dirdq[1]*zXdinv[6];
      dU1[4] += N[i]*dirdq[1]*zXdinv[7];
      dU1[5] += N[i]*dirdq[1]*zXdinv[8];
      dU1[6] += N[i]*dirdq[2]*zXdinv[6];
      dU1[7] += N[i]*dirdq[2]*zXdinv[7];
      dU1[8] += N[i]*dirdq[2]*zXdinv[8];
      matMatMult(dU1, T, tmp);
      matTransMatMult(T, tmp, dU1);

      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0]*dU0[0] + U0[3]*dU0[3] + U0[6]*dU0[6]);
      b[1] = dU0[4] + (U0[1]*dU0[1] + U0[4]*dU0[4] + U0[7]*dU0[7]);
      b[2] = dU0[1] + dU0[3] + (U0[0]*dU0[1] + U0[3]*dU0[4] + U0[6]*dU0[7] +
				dU0[0]*U0[1] + dU0[3]*U0[4] + dU0[6]*U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0]*dU1[0] + U0[3]*dU1[3] + U0[6]*dU1[6] + 
		       dU0[0]*U1[0] + dU0[3]*U1[3] + dU0[6]*U1[6]);
      b[4] = dU1[4] + (U0[1]*dU1[1] + U0[4]*dU1[4] + U0[7]*dU1[7] +
		       dU0[1]*U1[1] + dU0[4]*U1[4] + dU0[7]*U1[7]);
      b[5] = dU1[1] + dU1[3] + (U0[0]*dU1[1] + U0[3]*dU1[4] + U0[6]*dU1[7] +
				U1[0]*dU0[1] + U1[3]*dU0[4] + U1[6]*dU0[7] +
				dU0[0]*U1[1] + dU0[3]*U1[4] + dU0[6]*U1[7] +
				dU1[0]*U0[1] + dU1[3]*U0[4] + dU1[6]*U0[7]);
      b[6] = b[7] = 0.0;
      b += 8;
      dirdq += 3;
    }

    // Zero the contribution from the multiplier
    b[0] = b[1] = b[2] = b[3] = b[4] = b[5] = b[6] = b[7] = 0.0;
  }
}

/*
  Add the contribution from the geometric stiffness 
*/
void MITC9::addGmat( TacsScalar J[],
		     const TacsScalar scale,
		     const TacsScalar s[],
		     const double N[],
		     const double Na[],
		     const double Nb[],
		     const TacsScalar Ur[],
		     const TacsScalar dr[],
		     const TacsScalar Xdinv[],
		     const TacsScalar zXdinv[],
		     const TacsScalar T[],
		     const TacsScalar Xr[],
		     const TacsScalar dirdq[] ){
  // The gradient of the displacement field
  TacsScalar U0[9], U1[9], tmp[9];

  // Store the first derivatives w.r.t. the displacements
  // and rotations in the element
  TacsScalar dU0[7*9*NUM_NODES], dU1[7*9*NUM_NODES];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Pre-compute terms that are required for the geometric
  // stiffness matrix computation
  TacsScalar *du0 = dU0, *du1 = dU1; 
  for ( int i = 0; i < NUM_NODES; i++ ){
    TacsScalar t[3], ztmp[3];
    t[0] = Na[i]*Xdinv[0] + Nb[i]*Xdinv[3];
    t[1] = Na[i]*Xdinv[1] + Nb[i]*Xdinv[4];
    t[2] = Na[i]*Xdinv[2] + Nb[i]*Xdinv[5];
    matMultTrans(T, t, tmp);

    t[0] = Na[i]*zXdinv[0] + Nb[i]*zXdinv[3];
    t[1] = Na[i]*zXdinv[1] + Nb[i]*zXdinv[4];
    t[2] = Na[i]*zXdinv[2] + Nb[i]*zXdinv[5];
    matMultTrans(T, t, ztmp);

    for ( int k = 0; k < 3; k++ ){
      // dU0 = T[3*k+i]*tmp[3*k+j];
      du0[0] = T[3*k]*tmp[0];
      du0[1] = T[3*k]*tmp[1];
      du0[2] = T[3*k]*tmp[2];
      du0[3] = T[3*k+1]*tmp[0];
      du0[4] = T[3*k+1]*tmp[1];
      du0[5] = T[3*k+1]*tmp[2];
      du0[6] = T[3*k+2]*tmp[0];
      du0[7] = T[3*k+2]*tmp[1];
      du0[8] = T[3*k+2]*tmp[2];

      // dU1 = T[3*k+i]*ztmp[3*k+j];
      du1[0] = T[3*k]*ztmp[0];
      du1[1] = T[3*k]*ztmp[1];
      du1[2] = T[3*k]*ztmp[2];
      du1[3] = T[3*k+1]*ztmp[0];
      du1[4] = T[3*k+1]*ztmp[1];
      du1[5] = T[3*k+1]*ztmp[2];
      du1[6] = T[3*k+2]*ztmp[0];
      du1[7] = T[3*k+2]*ztmp[1];
      du1[8] = T[3*k+2]*ztmp[2];
      
      du0 += 9;
      du1 += 9;
    }

    for ( int k = 0; k < 4; k++ ){
      // Compute du0 = T^{T}*dur*Xdinv*T
      du0[0] = N[i]*dirdq[0]*Xdinv[6];
      du0[1] = N[i]*dirdq[0]*Xdinv[7];
      du0[2] = N[i]*dirdq[0]*Xdinv[8];
      du0[3] = N[i]*dirdq[1]*Xdinv[6];
      du0[4] = N[i]*dirdq[1]*Xdinv[7];
      du0[5] = N[i]*dirdq[1]*Xdinv[8];
      du0[6] = N[i]*dirdq[2]*Xdinv[6];
      du0[7] = N[i]*dirdq[2]*Xdinv[7];
      du0[8] = N[i]*dirdq[2]*Xdinv[8];
      matMatMult(du0, T, tmp);
      matTransMatMult(T, tmp, du0);

      // Add the contributions from the other 
      // derivatives
      TacsScalar drdq[9];
      drdq[0] = Na[i]*dirdq[0];
      drdq[3] = Na[i]*dirdq[1];
      drdq[6] = Na[i]*dirdq[2];
      drdq[1] = Nb[i]*dirdq[0];
      drdq[4] = Nb[i]*dirdq[1];
      drdq[7] = Nb[i]*dirdq[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, du1);

      du1[0] += N[i]*dirdq[0]*zXdinv[6];
      du1[1] += N[i]*dirdq[0]*zXdinv[7];
      du1[2] += N[i]*dirdq[0]*zXdinv[8];
      du1[3] += N[i]*dirdq[1]*zXdinv[6];
      du1[4] += N[i]*dirdq[1]*zXdinv[7];
      du1[5] += N[i]*dirdq[1]*zXdinv[8];
      du1[6] += N[i]*dirdq[2]*zXdinv[6];
      du1[7] += N[i]*dirdq[2]*zXdinv[7];
      du1[8] += N[i]*dirdq[2]*zXdinv[8];
      matMatMult(du1, T, tmp);
      matTransMatMult(T, tmp, du1);
   
      du0 += 9;
      du1 += 9;
      dirdq += 3;
    }
  }

  // Compute the derivatives of the strain w.r.t. the displacement
  // variables. This code takes advantage of the sparsity of the
  // derivatives to simplify the computations
  for ( int i = 0; i < 7*NUM_NODES; i++ ){
    const TacsScalar *dU0i = &dU0[9*i]; 
    const TacsScalar *dU1i = &dU1[9*i]; 
    for ( int j = i; j < 7*NUM_NODES; j++ ){
      const TacsScalar *dU0j = &dU0[9*j];
      const TacsScalar *dU1j = &dU1[9*j];

      // Compute the real indices
      int ii = 8*(i/7) + (i%7);
      int jj = 8*(j/7) + (j%7);
      int idx = 8*NUM_NODES*ii + jj;
      int sym = 8*NUM_NODES*jj + ii;
    
      // Compute the derivative of the in-plane strain
      TacsScalar b[6];
      b[0] = dU0i[0]*dU0j[0] + dU0i[3]*dU0j[3] + dU0i[6]*dU0j[6];
      b[1] = dU0i[1]*dU0j[1] + dU0i[4]*dU0j[4] + dU0i[7]*dU0j[7];
      b[2] = (dU0i[0]*dU0j[1] + dU0i[3]*dU0j[4] + dU0i[6]*dU0j[7] +
	      dU0j[0]*dU0i[1] + dU0j[3]*dU0i[4] + dU0j[6]*dU0i[7]);
      
      b[3] = (dU0i[0]*dU1j[0] + dU0i[3]*dU1j[3] + dU0i[6]*dU1j[6] + 
	      dU0j[0]*dU1i[0] + dU0j[3]*dU1i[3] + dU0j[6]*dU1i[6]);
      b[4] = (dU0i[1]*dU1j[1] + dU0i[4]*dU1j[4] + dU0i[7]*dU1j[7] +
	      dU0j[1]*dU1i[1] + dU0j[4]*dU1i[4] + dU0j[7]*dU1i[7]);
      b[5] = (dU0i[0]*dU1j[1] + dU0i[3]*dU1j[4] + dU0i[6]*dU1j[7] +
	      dU1i[0]*dU0j[1] + dU1i[3]*dU0j[4] + dU1i[6]*dU0j[7] +
	      dU0j[0]*dU1i[1] + dU0j[3]*dU1i[4] + dU0j[6]*dU1i[7] +
	      dU1j[0]*dU0i[1] + dU1j[3]*dU0i[4] + dU1j[6]*dU0i[7]);
      
      TacsScalar Jadd = 
	scale*(b[0]*s[0] + b[1]*s[1] + b[2]*s[2] +
	       b[3]*s[3] + b[4]*s[4] + b[5]*s[5]);
      
      // Add the values symmetrically
      J[idx] += Jadd;
      if (ii != jj){
	J[sym] += Jadd;
      }
    }
  }

  // Add the contributions from the second derivatives of the quaternions
  for ( int i = 0; i < NUM_NODES; i++ ){
    // d = N[i]*C(q_{i})^{T}*n - these terms only get added along
    // each diagonal in the matrix
    
    // Extract the normal from the frame
    TacsScalar normal[3];
    normal[0] = Xr[2];
    normal[1] = Xr[5];
    normal[2] = Xr[8];
    Xr += 9;
    
    // Compute the second derivative w.r.t. the quaternion
    TacsScalar dCtndq[3*9];
    computeQtr2ndDeriv(normal, dCtndq);
    const TacsScalar *dC = dCtndq;

    // Compute the partials derivatives w.r.t. eta,eps 
    for ( int ii = 0; ii < 9; ii++ ){
      // Compute dU0 = T^{T}*dUr*Xdinv*T
      dU0[0] = N[i]*dC[0]*Xdinv[6];
      dU0[1] = N[i]*dC[0]*Xdinv[7];
      dU0[2] = N[i]*dC[0]*Xdinv[8];
      dU0[3] = N[i]*dC[1]*Xdinv[6];
      dU0[4] = N[i]*dC[1]*Xdinv[7];
      dU0[5] = N[i]*dC[1]*Xdinv[8];
      dU0[6] = N[i]*dC[2]*Xdinv[6];
      dU0[7] = N[i]*dC[2]*Xdinv[7];
      dU0[8] = N[i]*dC[2]*Xdinv[8];
      matMatMult(dU0, T, tmp);
      matTransMatMult(T, tmp, dU0);

      // Add the contributions from the other 
      // derivatives
      TacsScalar drdq[9];
      drdq[0] = Na[i]*dC[0];
      drdq[3] = Na[i]*dC[1];
      drdq[6] = Na[i]*dC[2];
      drdq[1] = Nb[i]*dC[0];
      drdq[4] = Nb[i]*dC[1];
      drdq[7] = Nb[i]*dC[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, dU1);

      dU1[0] += N[i]*dC[0]*zXdinv[6];
      dU1[1] += N[i]*dC[0]*zXdinv[7];
      dU1[2] += N[i]*dC[0]*zXdinv[8];
      dU1[3] += N[i]*dC[1]*zXdinv[6];
      dU1[4] += N[i]*dC[1]*zXdinv[7];
      dU1[5] += N[i]*dC[1]*zXdinv[8];
      dU1[6] += N[i]*dC[2]*zXdinv[6];
      dU1[7] += N[i]*dC[2]*zXdinv[7];
      dU1[8] += N[i]*dC[2]*zXdinv[8];
      matMatMult(dU1, T, tmp);
      matTransMatMult(T, tmp, dU1);

      // Compute the derivative of the in-plane strain
      TacsScalar b[6];
      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0]*dU0[0] + U0[3]*dU0[3] + U0[6]*dU0[6]);
      b[1] = dU0[4] + (U0[1]*dU0[1] + U0[4]*dU0[4] + U0[7]*dU0[7]);
      b[2] = dU0[1] + dU0[3] + (U0[0]*dU0[1] + U0[3]*dU0[4] + U0[6]*dU0[7] +
				dU0[0]*U0[1] + dU0[3]*U0[4] + dU0[6]*U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0]*dU1[0] + U0[3]*dU1[3] + U0[6]*dU1[6] + 
		       dU0[0]*U1[0] + dU0[3]*U1[3] + dU0[6]*U1[6]);
      b[4] = dU1[4] + (U0[1]*dU1[1] + U0[4]*dU1[4] + U0[7]*dU1[7] +
		       dU0[1]*U1[1] + dU0[4]*U1[4] + dU0[7]*U1[7]);
      b[5] = dU1[1] + dU1[3] + (U0[0]*dU1[1] + U0[3]*dU1[4] + U0[6]*dU1[7] +
				U1[0]*dU0[1] + U1[3]*dU0[4] + U1[6]*dU0[7] +
				dU0[0]*U1[1] + dU0[3]*U1[4] + dU0[6]*U1[7] +
				dU1[0]*U0[1] + dU1[3]*U0[4] + dU1[6]*U0[7]);

      TacsScalar Jadd = 
	scale*(b[0]*s[0] + b[1]*s[1] + b[2]*s[2] +
	       b[3]*s[3] + b[4]*s[4] + b[5]*s[5]);

      if (ii < 3){
	int iv = 8*i + 3;
	int jv = 8*i + 4 + ii;
	J[(8*NUM_NODES)*iv + jv] += Jadd;
	J[(8*NUM_NODES)*jv + iv] += Jadd;
      }
      else {
	if (ii == 3){
	  int iv = 8*i + 4;
	  J[(8*NUM_NODES+1)*iv] += Jadd;
	}
	else if (ii == 4){
	  int iv = 8*i + 4, jv = 8*i + 5;
	  J[(8*NUM_NODES)*iv + jv] += Jadd;
	  J[(8*NUM_NODES)*jv + iv] += Jadd;
	}
	else if (ii == 5){
	  int iv = 8*i + 4, jv = 8*i + 6;
	  J[(8*NUM_NODES)*iv + jv] += Jadd;
	  J[(8*NUM_NODES)*jv + iv] += Jadd;
	}
	else if (ii == 6){
	  int iv = 8*i + 5;
	  J[(8*NUM_NODES+1)*iv] += Jadd;
	}
	else if (ii == 7){
	  int iv = 8*i + 5, jv = 8*i + 6;
	  J[(8*NUM_NODES)*iv + jv] += Jadd;
	  J[(8*NUM_NODES)*jv + iv] += Jadd;
	}
	else if (ii == 8){
	  int iv = 8*i + 6;
	  J[(8*NUM_NODES+1)*iv] += Jadd;
	}
      }

      dC += 3;
    }
  }
}

/*
  Compute the value of the tensorial strain at the tying points
  in the element. 

  This code evaluates the tensorial shear strain values at the
  tying points which consist of the 2-point Gauss quadrature
  points in one direction, and the nodal locations along the other
  direction. 

  input:
  X:     the initial values of the nodal coordinates 
  vars:  the values of the variables
  dr:    the director values at every node in the element
  
  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::computeTyingStrain( TacsScalar g13[], 
				TacsScalar g23[],
				const TacsScalar X[],
				const TacsScalar Xr[],
				const TacsScalar vars[],
				const TacsScalar dir[] ){
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t,  t,  -t,   t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0,  s, s};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions 
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3], Xb[3];
    innerProduct(Na, X, Xa);
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    g13[k] = 0.5*(vecDot(Xa, d) + vecDot(Ua, fn) + vecDot(Ua, d));
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0,  s, -s, 0.0, s};
  const double g23_vpts[] = {-t,  -t, -t,  t,   t, t};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions 
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3], Xb[3];
    innerProduct(Na, X, Xa);
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    g23[k] = 0.5*(vecDot(Xb, d) + vecDot(Ub, fn) + vecDot(Ub, d));
  }
}

/*
  Compute the value of the tensorial strain and the derivative of the
  tensorial strain w.r.t. the state variables at the tying points in
  the element.

  This code evaluates the tensorial shear strain values at the tying
  points which consist of the 2-point Gauss quadrature points in one
  direction, and the nodal locations along the other direction.

  input:
  X:     the initial values of the nodal coordinates 
  vars:  the values of the variables
  dr:    the director values at every node in the element
  
  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::computeTyingBmat( TacsScalar g13[], 
			      TacsScalar g23[],
			      TacsScalar B13[],
			      TacsScalar B23[],
			      const TacsScalar X[],
			      const TacsScalar Xr[],
			      const TacsScalar vars[],
			      const TacsScalar dir[],
			      const TacsScalar dirdq[] ){
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t,  t,  -t,   t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0,  s, s};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions 
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    g13[k] = 0.5*(vecDot(Xa, d) + vecDot(Ua, fn) + vecDot(Ua, d));

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xa[0] + Ua[0];
    t[1] = Xa[1] + Ua[1];
    t[2] = Xa[2] + Ua[2];
  
    TacsScalar *b13 = &B13[8*NUM_NODES*k];
    const TacsScalar *dq = dirdq;
    for ( int i = 0; i < NUM_NODES; i++ ){
      // Compute the derivatives w.r.t. U
      b13[0] = 0.5*Na[i]*(fn[0] + d[0]);
      b13[1] = 0.5*Na[i]*(fn[1] + d[1]);
      b13[2] = 0.5*Na[i]*(fn[2] + d[2]);
      
      // Compute the derivatives w.r.t. q
      b13[3] = 0.5*N[i]*(dq[0]*t[0] + dq[1]*t[1] + dq[2]*t[2]);
      b13[4] = 0.5*N[i]*(dq[3]*t[0] + dq[4]*t[1] + dq[5]*t[2]);
      b13[5] = 0.5*N[i]*(dq[6]*t[0] + dq[7]*t[1] + dq[8]*t[2]);
      b13[6] = 0.5*N[i]*(dq[9]*t[0] + dq[10]*t[1] + dq[11]*t[2]);
      b13[7] = 0.0;

      b13 += 8;
      dq += 12;
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0,  s, -s, 0.0, s};
  const double g23_vpts[] = {-t,  -t, -t,  t,   t, t};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions 
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    g23[k] = 0.5*(vecDot(Xb, d) + vecDot(Ub, fn) + vecDot(Ub, d));

    TacsScalar t[3];
    t[0] = Xb[0] + Ub[0];
    t[1] = Xb[1] + Ub[1];
    t[2] = Xb[2] + Ub[2];

    TacsScalar *b23 = &B23[8*NUM_NODES*k];
    const TacsScalar *dq = dirdq;
    for ( int i = 0; i < NUM_NODES; i++ ){
      // Compute the derivatives w.r.t. U
      b23[0] = 0.5*Nb[i]*(fn[0] + d[0]);
      b23[1] = 0.5*Nb[i]*(fn[1] + d[1]);
      b23[2] = 0.5*Nb[i]*(fn[2] + d[2]);
      
      // Compute the derivatives w.r.t. q
      b23[3] = 0.5*N[i]*(dq[0]*t[0] + dq[1]*t[1] + dq[2]*t[2]);
      b23[4] = 0.5*N[i]*(dq[3]*t[0] + dq[4]*t[1] + dq[5]*t[2]);
      b23[5] = 0.5*N[i]*(dq[6]*t[0] + dq[7]*t[1] + dq[8]*t[2]);
      b23[6] = 0.5*N[i]*(dq[9]*t[0] + dq[10]*t[1] + dq[11]*t[2]);
      b23[7] = 0.0;

      b23 += 8;
      dq += 12;
    }
  }
}

/*
  Compute the value of the tensorial strain and the derivative of the
  tensorial strain w.r.t. the state variables at the tying points in
  the element.

  This code evaluates the tensorial shear strain values at the tying
  points which consist of the 2-point Gauss quadrature points in one
  direction, and the nodal locations along the other direction.

  input:
  X:     the initial values of the nodal coordinates 
  vars:  the values of the variables
  dr:    the director values at every node in the element
  
  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::addTyingGmat( TacsScalar J[],
			  const TacsScalar w13[], 
			  const TacsScalar w23[],
			  const TacsScalar X[],
			  const TacsScalar Xr[],
			  const TacsScalar vars[],
			  const TacsScalar dir[],
			  const TacsScalar dirdq[] ){
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  const int iv[] = {3, 3, 3, 4, 4, 4, 5, 5, 6};
  const int jv[] = {4, 5, 6, 4, 5, 6, 5, 6, 6};

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t,  t,  -t,   t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0,  s, s};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions 
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xa[0] + Ua[0];
    t[1] = Xa[1] + Ua[1];
    t[2] = Xa[2] + Ua[2];

    const TacsScalar *xr = Xr;
    for ( int i = 0; i < NUM_NODES; i++ ){
      // Extract the normal from the frame
      TacsScalar normal[3];
      normal[0] = xr[2];
      normal[1] = xr[5];
      normal[2] = xr[8];
      xr += 9;

      // Compute the second derivative of the quaternion
      TacsScalar dCtndq[3*9];
      computeQtr2ndDeriv(normal, dCtndq);

      // Compute the second derivatives from the quaternions alone
      const TacsScalar *dC = dCtndq;
      for ( int ii = 0; ii < 9; ii++ ){
	TacsScalar Jadd = 0.5*w13[k]*N[i]*vecDot(t, dC);
	J[(8*NUM_NODES)*(8*i + iv[ii]) + (8*i + jv[ii])] += Jadd;
	if (iv[ii] != jv[ii]){
	  J[(8*NUM_NODES)*(8*i + jv[ii]) + (8*i + iv[ii])] += Jadd;
	}
	dC += 3;
      }
      
      // Set the derivative of the director w.r.t. q
      const TacsScalar *dq = dirdq;
      for ( int j = 0; j < NUM_NODES; j++ ){
	for ( int ii = 0; ii < 4; ii++ ){ // Loop over quaternions
	  for ( int jj = 0; jj < 3; jj++ ){ // Loop over displacements
	    TacsScalar Jadd = 0.5*w13[k]*N[j]*Na[i]*dq[3*ii + jj]; 
	    J[8*NUM_NODES*(8*i + jj) + 8*j + 3+ii] += Jadd;
	    J[8*NUM_NODES*(8*j + 3+ii) + 8*i + jj] += Jadd;
	  }
	}
	dq += 12;
      }
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0,  s, -s, 0.0, s};
  const double g23_vpts[] = {-t,  -t, -t,  t,   t, t};

  for ( int k = 0; k < 6; k++ ){
    // Evaluate the element shape functions
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    TacsScalar t[3];
    t[0] = Xb[0] + Ub[0];
    t[1] = Xb[1] + Ub[1];
    t[2] = Xb[2] + Ub[2];

    const TacsScalar *xr = Xr;
    for ( int i = 0; i < NUM_NODES; i++ ){
      // Extract the normal from the frame
      TacsScalar normal[3];
      normal[0] = xr[2];
      normal[1] = xr[5];
      normal[2] = xr[8];
      xr += 9;

      // Compute the second derivative of the quaternion
      TacsScalar dCtndq[3*9];
      computeQtr2ndDeriv(normal, dCtndq);

      // Compute the second derivatives from the quaternions alone
      const TacsScalar *dC = dCtndq;
      for ( int ii = 0; ii < 9; ii++ ){
	TacsScalar Jadd = 0.5*w23[k]*N[i]*vecDot(t, dC);
	J[(8*NUM_NODES)*(8*i + iv[ii]) + (8*i + jv[ii])] += Jadd;
	if (iv[ii] != jv[ii]){
	  J[(8*NUM_NODES)*(8*i + jv[ii]) + (8*i + iv[ii])] += Jadd;
	}
	dC += 3;
      }
      
      // Set the derivative of the director w.r.t. q
      const TacsScalar *dq = dirdq;
      for ( int j = 0; j < NUM_NODES; j++ ){
	for ( int ii = 0; ii < 4; ii++ ){ // Loop over quaternions
	  for ( int jj = 0; jj < 3; jj++ ){ // Loop over displacements
	    TacsScalar Jadd = 0.5*w23[k]*N[j]*Nb[i]*dq[3*ii + jj];
	    J[8*NUM_NODES*(8*i + jj) + 8*j + 3+ii] += Jadd;
	    J[8*NUM_NODES*(8*j + 3+ii) + 8*i + jj] += Jadd;
	  }
	}
	dq += 12;
      }
    }
  }
}

/*
  Add the assumed strain distribution interpolated from the tying
  points to the strain at the present point. 

  This code interpolates the strain from the tying points, which
  interpolates the tensorial components of the strain, and adds these
  values to the strain. This involes an interpolation and then a
  transformation. Note that not all components of the transformation
  are required. In general, the transformation takes the form:

  e = A^{T}*g*A

  where g are the tensorial components of the strain and A is the
  transformation matrix (which is not necessarily orthonormal!)

  input:
  N13:  the shape functions for the g13 tensorial strain
  N23:  the shape functions for the g23 tensorial strain
  g13:  the g13 tensorial strain at the tying points
  g23:  the g23 tensorial strain at the tying points
  T:    the transformation matrix

  output:
  e:    the shear strain components are over-written
*/
void MITC9::addTyingStrain( TacsScalar e[],
			    const double N13[],
			    const double N23[],
			    const TacsScalar g13[],
			    const TacsScalar g23[],
			    const TacsScalar Xdinv[],
			    const TacsScalar T[] ){
  // Compute the strain using the assumed strain distribution
  // and the strain evaluated at the tying points
  TacsScalar G13, G23;
  G13 = (g13[0]*N13[0] + g13[1]*N13[1] + g13[2]*N13[2] +
	 g13[3]*N13[3] + g13[4]*N13[4] + g13[5]*N13[5]);
  G23 = (g23[0]*N23[0] + g23[1]*N23[1] + g23[2]*N23[2] +
	 g23[3]*N23[3] + g23[4]*N23[4] + g23[5]*N23[5]);

  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0]*T[0] + Xdinv[1]*T[3] + Xdinv[2]*T[6];
  A12 = Xdinv[0]*T[1] + Xdinv[1]*T[4] + Xdinv[2]*T[7];

  A21 = Xdinv[3]*T[0] + Xdinv[4]*T[3] + Xdinv[5]*T[6];
  A22 = Xdinv[3]*T[1] + Xdinv[4]*T[4] + Xdinv[5]*T[7];

  // Compute and set the final strain values
  // e = 2*A^{T}*G
  e[6] = 2.0*(A12*G13 + A22*G23);
  e[7] = 2.0*(A11*G13 + A21*G23);
}

/*
  Add the tying coefficients required to compute the geometric
  stiffness matrix term
*/
void MITC9::addTyingGmatWeights( TacsScalar w13[],
				 TacsScalar w23[],
				 const TacsScalar scalar,
				 const TacsScalar s[],
				 const double N13[],
				 const double N23[],
				 const TacsScalar Xdinv[],
				 const TacsScalar T[] ){
  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0]*T[0] + Xdinv[1]*T[3] + Xdinv[2]*T[6];
  A12 = Xdinv[0]*T[1] + Xdinv[1]*T[4] + Xdinv[2]*T[7];

  A21 = Xdinv[3]*T[0] + Xdinv[4]*T[3] + Xdinv[5]*T[6];
  A22 = Xdinv[3]*T[1] + Xdinv[4]*T[4] + Xdinv[5]*T[7];

  // Add the contributions to the weights
  for ( int k = 0; k < 6; k++ ){
    w13[k] += 2.0*scalar*(s[6]*A12 + s[7]*A11)*N13[k];
    w23[k] += 2.0*scalar*(s[6]*A22 + s[7]*A21)*N23[k];
  }
}

/*
  Add the contribution from the assumed strain derivatives
  interpolated from the tying points to the derivative of the strain
  at the present point.

  Note that this code is analagous to the addTyingStrain function
  except for the B matrix.

  input:
  N13:  the shape functions for the g13 tensorial strain
  N23:  the shape functions for the g23 tensorial strain
  B13:  the g13 tensorial strain at the tying points
  B23:  the g23 tensorial strain at the tying points
  T:    the transformation matrix

  output:
  B:    the derivative of the strain
*/
void MITC9::addTyingBmat( TacsScalar B[],
			  const double N13[],
			  const double N23[],
			  const TacsScalar B13[],
			  const TacsScalar B23[],
			  const TacsScalar Xdinv[],
			  const TacsScalar T[] ){
  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0]*T[0] + Xdinv[1]*T[3] + Xdinv[2]*T[6];
  A12 = Xdinv[0]*T[1] + Xdinv[1]*T[4] + Xdinv[2]*T[7];

  A21 = Xdinv[3]*T[0] + Xdinv[4]*T[3] + Xdinv[5]*T[6];
  A22 = Xdinv[3]*T[1] + Xdinv[4]*T[4] + Xdinv[5]*T[7];

  const int offset = 8*NUM_NODES;
  for ( int k = 0; k < offset; k++ ){
    // Compute the strain using the assumed strain distribution
    // and the strain evaluated at the tying points
    TacsScalar G13, G23;
    G13 = (B13[0]*N13[0] + B13[offset]*N13[1] +
	   B13[2*offset]*N13[2] + B13[3*offset]*N13[3] + 
	   B13[4*offset]*N13[4] + B13[5*offset]*N13[5]);
    G23 = (B23[0]*N23[0] + B23[offset]*N23[1] +
	   B23[2*offset]*N23[2] + B23[3*offset]*N23[3] + 
	   B23[4*offset]*N23[4] + B23[5*offset]*N23[5]);

    // Compute and set the final strain values
    // e = 2*A^{T}*G
    B[6] = 2.0*(A12*G13 + A22*G23);
    B[7] = 2.0*(A11*G13 + A21*G23);
    
    B += 8;
    B13 += 1;
    B23 += 1;
  }
}

/*
  Evaluate the in-plane rotation penalty term 

  This term is used to add an artificial drilling stiffness to the
  plate to enable the use of full rotations at the element nodes.

  The rotation penalty term is calculated as follows:

  rot = Xa^{T}*C*(Xb + Ub) - Xb^{T}*C*(Xa + Ua)

  In this code, the approximate rotation matrix is interpolated from
  the nodes of the mesh and is computed directly from the directors
*/
TacsScalar MITC9::computeRotPenalty( const double N[],
				     const TacsScalar Xa[],
				     const TacsScalar Xb[],
				     const TacsScalar Ua[],
				     const TacsScalar Ub[],
				     const TacsScalar vars[] ){
  TacsScalar Ci[9];
  Ci[0] = Ci[1] = Ci[2] = 0.0;
  Ci[3] = Ci[4] = Ci[5] = 0.0;
  Ci[6] = Ci[7] = Ci[8] = 0.0;
  
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set the pointer to the quaternions
    const TacsScalar *q = &vars[8*i + 3];

    // Compute the full rotation matrix
    TacsScalar C[9];
    C[0] = 1.0 - 2.0*(q[2]*q[2] + q[3]*q[3]);
    C[1] = 2.0*(q[1]*q[2] + q[3]*q[0]);
    C[2] = 2.0*(q[1]*q[3] - q[2]*q[0]);

    C[3] = 2.0*(q[2]*q[1] - q[3]*q[0]);
    C[4] = 1.0 - 2.0*(q[1]*q[1] + q[3]*q[3]);
    C[5] = 2.0*(q[2]*q[3] + q[1]*q[0]);

    C[6] = 2.0*(q[3]*q[1] + q[2]*q[0]);
    C[7] = 2.0*(q[3]*q[2] - q[1]*q[0]);
    C[8] = 1.0 - 2.0*(q[1]*q[1] + q[2]*q[2]);

    // Interpolate the rotation matrix
    Ci[0] += N[i]*C[0];
    Ci[1] += N[i]*C[1];
    Ci[2] += N[i]*C[2];
    Ci[3] += N[i]*C[3];
    Ci[4] += N[i]*C[4];
    Ci[5] += N[i]*C[5];
    Ci[6] += N[i]*C[6];
    Ci[7] += N[i]*C[7];
    Ci[8] += N[i]*C[8];
  }

  // Compute va = C*(Xa + Ua)
  TacsScalar va[3];
  matMult(Ci, Xa, va);
  matMultAdd(Ci, Ua, va);

  // Compute vb = C*(Xa + Ua)
  TacsScalar vb[3];
  matMult(Ci, Xb, vb);
  matMultAdd(Ci, Ub, vb);

  // Compute the penalty term
  return vecDot(Xa, vb) - vecDot(Xb, va);
}

/*
  Evaluate the derivative of the in-plane rotation penalty term 

  This term is used to add an artificial drilling stiffness to the
  plate to enable the use of full rotations at the element nodes.

  The rotation penalty term is calculated as follows:

  rot = Xa^{T}*C*(Xb + Ub) - Xb^{T}*C*(Xa + Ua)

  In this code, the approximate rotation matrix is interpolated from
  the nodes of the mesh and is computed directly from the directors
*/
TacsScalar MITC9::computeBRotPenalty( TacsScalar brot[],
				      const double N[],
				      const double Na[],
				      const double Nb[],
				      const TacsScalar Xa[],
				      const TacsScalar Xb[],
				      const TacsScalar Ua[],
				      const TacsScalar Ub[],
				      const TacsScalar vars[] ){
  TacsScalar Ci[9];
  Ci[0] = Ci[1] = Ci[2] = 0.0;
  Ci[3] = Ci[4] = Ci[5] = 0.0;
  Ci[6] = Ci[7] = Ci[8] = 0.0;
  
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set the pointer to the quaternions
    const TacsScalar *q = &vars[8*i + 3];

    // Compute the rotation matrix
    TacsScalar C[9];
    C[0] = 1.0 - 2.0*(q[2]*q[2] + q[3]*q[3]);
    C[1] = 2.0*(q[1]*q[2] + q[3]*q[0]);
    C[2] = 2.0*(q[1]*q[3] - q[2]*q[0]);

    C[3] = 2.0*(q[2]*q[1] - q[3]*q[0]);
    C[4] = 1.0 - 2.0*(q[1]*q[1] + q[3]*q[3]);
    C[5] = 2.0*(q[2]*q[3] + q[1]*q[0]);

    C[6] = 2.0*(q[3]*q[1] + q[2]*q[0]);
    C[7] = 2.0*(q[3]*q[2] - q[1]*q[0]);
    C[8] = 1.0 - 2.0*(q[1]*q[1] + q[2]*q[2]);

    // Interpolate the rotation matrix to the nodes
    Ci[0] += N[i]*C[0];
    Ci[1] += N[i]*C[1];
    Ci[2] += N[i]*C[2];
    Ci[3] += N[i]*C[3];
    Ci[4] += N[i]*C[4];
    Ci[5] += N[i]*C[5];
    Ci[6] += N[i]*C[6];
    Ci[7] += N[i]*C[7];
    Ci[8] += N[i]*C[8];
  }

  // Compute va = Ci*(Xa + Ua)
  TacsScalar va[3], ta[3];
  ta[0] = Xa[0] + Ua[0];
  ta[1] = Xa[1] + Ua[1];
  ta[2] = Xa[2] + Ua[2];
  matMult(Ci, ta, va);

  // Compute vb = Ci*(Xb + Ub)
  TacsScalar vb[3], tb[3];
  tb[0] = Xb[0] + Ub[0];
  tb[1] = Xb[1] + Ub[1];
  tb[2] = Xb[2] + Ub[2];
  matMult(Ci, tb, vb);

  // Compute wa = Ci^{T}*Xa and wb = Ci^{T}*Xb
  TacsScalar wa[3], wb[3];
  matMultTrans(Ci, Xa, wa);
  matMultTrans(Ci, Xb, wb);

  // Add the values from the
  TacsScalar *b = brot;
  for ( int i = 0; i < NUM_NODES; i++ ){
    // wa*Nb - wb*Na
    b[0] = Nb[i]*wa[0] - Na[i]*wb[0];
    b[1] = Nb[i]*wa[1] - Na[i]*wb[1];
    b[2] = Nb[i]*wa[2] - Na[i]*wb[2];

    // Add the terms from the derivatives w.r.t. Ci. Note that 
    // each term adds a contribution which takes the following form:
    //  N[i]*(Xa^{T}*D*tb - Xb^{T}*D*ta)
    const TacsScalar *q = &vars[8*i + 3];
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0*q[3];
    Q[2] =-2.0*q[2];

    Q[3] =-2.0*q[3];
    Q[4] = 0.0;
    Q[5] = 2.0*q[1];

    Q[6] = 2.0*q[2];
    Q[7] =-2.0*q[1];
    Q[8] = 0.0;
    b[3] = N[i]*(mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0*q[2];
    Q[2] = 2.0*q[3];

    Q[3] = 2.0*q[2];
    Q[4] =-4.0*q[1];
    Q[5] = 2.0*q[0];

    Q[6] = 2.0*q[3];
    Q[7] =-2.0*q[0];
    Q[8] =-4.0*q[1];
    b[4] = N[i]*(mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[2]
    Q[0] =-4.0*q[2];
    Q[1] = 2.0*q[1];
    Q[2] =-2.0*q[0];

    Q[3] = 2.0*q[1];
    Q[4] = 0.0;
    Q[5] = 2.0*q[3];

    Q[6] = 2.0*q[0];
    Q[7] = 2.0*q[3];
    Q[8] =-4.0*q[2];
    b[5] = N[i]*(mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[3]
    Q[0] =-4.0*q[3];
    Q[1] = 2.0*q[0];
    Q[2] = 2.0*q[1];

    Q[3] =-2.0*q[0];
    Q[4] =-4.0*q[3];
    Q[5] = 2.0*q[2];

    Q[6] = 2.0*q[1];
    Q[7] = 2.0*q[2];
    Q[8] = 0.0;
    b[6] = N[i]*(mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    b[7] = 0.0;

    b += 8;
  }

  // Compute the penalty term
  return  vecDot(Xa, vb) - vecDot(Xb, va);
}

/*
  Add the second derivative of the in-plane penalty term to the
  stiffness matrix.
*/
void MITC9::addGRotMat( TacsScalar J[],
			const TacsScalar scale, 
			const double N[],
			const double Na[],
			const double Nb[],
			const TacsScalar Xa[],
			const TacsScalar Xb[],
			const TacsScalar Ua[],
			const TacsScalar Ub[],
			const TacsScalar vars[] ){
  // Compute va = Ci*(Xa + Ua)
  TacsScalar ta[3];
  ta[0] = Xa[0] + Ua[0];
  ta[1] = Xa[1] + Ua[1];
  ta[2] = Xa[2] + Ua[2];

  // Compute vb = Ci*(Xb + Ub)
  TacsScalar tb[3];
  tb[0] = Xb[0] + Ub[0];
  tb[1] = Xb[1] + Ub[1];
  tb[2] = Xb[2] + Ub[2];

  // Compute the second derivative w.r.t. the quaternion
  TacsScalar dCXa[3*9], dCXb[3*9];
  computeQtr2ndDeriv(Xa, dCXa);
  computeQtr2ndDeriv(Xb, dCXb);

  // Pre-compute terms for the second derivatives
  TacsScalar Wa[12*NUM_NODES], Wb[12*NUM_NODES];
  TacsScalar *wa = Wa, *wb = Wb;

  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i + 3];
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0*q[3];
    Q[2] =-2.0*q[2];

    Q[3] =-2.0*q[3];
    Q[4] = 0.0;
    Q[5] = 2.0*q[1];

    Q[6] = 2.0*q[2];
    Q[7] =-2.0*q[1];
    Q[8] = 0.0;
    matMultTrans(Q, Xa, &wa[0]);
    matMultTrans(Q, Xb, &wb[0]);

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0*q[2];
    Q[2] = 2.0*q[3];

    Q[3] = 2.0*q[2];
    Q[4] =-4.0*q[1];
    Q[5] = 2.0*q[0];

    Q[6] = 2.0*q[3];
    Q[7] =-2.0*q[0];
    Q[8] =-4.0*q[1];
    matMultTrans(Q, Xa, &wa[3]);
    matMultTrans(Q, Xb, &wb[3]);

    // Derivative w.r.t. q[2]
    Q[0] =-4.0*q[2];
    Q[1] = 2.0*q[1];
    Q[2] =-2.0*q[0];

    Q[3] = 2.0*q[1];
    Q[4] = 0.0;
    Q[5] = 2.0*q[3];

    Q[6] = 2.0*q[0];
    Q[7] = 2.0*q[3];
    Q[8] =-4.0*q[2];
    matMultTrans(Q, Xa, &wa[6]);
    matMultTrans(Q, Xb, &wb[6]);

    // Derivative w.r.t. q[3]
    Q[0] =-4.0*q[3];
    Q[1] = 2.0*q[0];
    Q[2] = 2.0*q[1];

    Q[3] =-2.0*q[0];
    Q[4] =-4.0*q[3];
    Q[5] = 2.0*q[2];

    Q[6] = 2.0*q[1];
    Q[7] = 2.0*q[2];
    Q[8] = 0.0;
    matMultTrans(Q, Xa, &wa[9]);
    matMultTrans(Q, Xb, &wb[9]);

    // Increment the pointers
    wa += 12;  wb += 12;
  }

  const int iv[] = {3, 3, 3, 4, 4, 4, 5, 5, 6};
  const int jv[] = {4, 5, 6, 4, 5, 6, 5, 6, 6};

  // Add the contributions from the second derivatives of the quaternions
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Add the second derivative coupling between the drilling
    // rotation and the in-plane displacement
    wa = Wa;  wb = Wb;
    for ( int j = 0; j < NUM_NODES; j++ ){
      for ( int ii = 0; ii < 4; ii++ ){ // Loop over the quaternions
	for ( int jj = 0; jj < 3; jj++ ){ // Loop over the displacements
	  J[8*NUM_NODES*(8*i + jj) + 8*j + 3+ii] += 
	    scale*N[j]*(Nb[i]*wa[3*ii+jj] - Na[i]*wb[3*ii+jj]);

	  J[8*NUM_NODES*(8*j + 3+ii) + 8*i + jj] +=
	    scale*N[j]*(Nb[i]*wa[3*ii+jj] - Na[i]*wb[3*ii+jj]);
	}
      }
      wa += 12;  wb += 12;
    }

    // Set pointers to the second derivatives
    const TacsScalar *dCa = dCXa, *dCb = dCXb;

    // Compute the second derivatives from the quaternions alone
    for ( int ii = 0; ii < 9; ii++ ){
      TacsScalar Jadd = scale*N[i]*(vecDot(tb, dCa) - vecDot(ta, dCb)); 
      J[(8*NUM_NODES)*(8*i + iv[ii]) + (8*i + jv[ii])] += Jadd;
      if (iv[ii] != jv[ii]){
	J[(8*NUM_NODES)*(8*i + jv[ii]) + (8*i + iv[ii])] += Jadd;
      }
      dCa += 3;  dCb += 3;
    }
  }
}

/*
  Get the constitutive object
*/
TACSConstitutive *MITC9::getConstitutive(){
  return stiff;
}

/*
  Return the number of quadrature points
*/
int MITC9::getNumGaussPts(){ 
  return ORDER*ORDER;
}

/*
  Return the quadrature points and weights
*/
double MITC9::getGaussWtsPts( const int num, double pt[] ){
  int m = (int)(num/ORDER);
  int n = num % ORDER;
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];    
  
  return gaussWts[n]*gaussWts[m];
}

/*
  Get the values of the shape functions
*/
void MITC9::getShapeFunctions( const double pt[], double N[] ){
  computeShapeFunc(pt[0], pt[1], N);
}

/*
  Retrieve the determinant of the Jacobian transformation matrix
*/
TacsScalar MITC9::getDetJacobian( const double pt[],
                                  const TacsScalar X[] ){
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);
  
  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  TacsScalar h = inv3x3(Xd, Xdinv);

  return h;
}

/*
  Evaluate the strain at a parametric point within the element
*/
void MITC9::getStrain( TacsScalar e[],
                       const double pt[],
                       const TacsScalar X[],
                       const TacsScalar vars[] ){
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3*NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);
  
  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  TacsScalar h = inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);
  
  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9]; 
  
  // Take the cross product to find the normal direction
  TacsScalar normal[3];
  crossProduct(1.0, Xa, Xb, normal);
  TacsScalar nrm = sqrt(vecDot(normal, normal));
  vecScale(1.0/nrm, normal);

  // Scale the Xa direction so that it is a unit vector
  nrm = sqrt(vecDot(Xa, Xa));
  vecScale(1.0/nrm, Xa);
  
  // Compute the second perpendicular direction 
  crossProduct(1.0, normal, Xa, Xb);
  
  // Assemble the transformation matrix
  assembleFrame(Xa, Xb, normal, T);
  
  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);
  
  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);
  
  // Compute the displacement-based strain
  evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
  
  // Add the contribution from the tying strain
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
}

/*
  Add the derivative of the product of the array esens with the strain
  with respect to the state variables
*/
void MITC9::addStrainSVSens( TacsScalar sens[],
                             const double pt[], 
                             const TacsScalar scale,
                             const TacsScalar esens[], 
                             const TacsScalar X[],
                             const TacsScalar vars[] ){
  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the derivatives of the directors
  TacsScalar dir[3*NUM_NODES], dirdq[12*NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Set the parametric locations
  const double u = pt[0];
  const double v = pt[1];

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  TacsScalar Xdinv[9];
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9]; 

  // Compute the cross product to find the normal direction
  TacsScalar normal[3];
  crossProduct(1.0, Xa, Xb, normal);
  TacsScalar nrm = sqrt(vecDot(normal, normal));
  vecScale(1.0/nrm, normal);
  
  // Scale the Xa direction so that it is a unit vector
  nrm = sqrt(vecDot(Xa, Xa));
  vecScale(1.0/nrm, Xa);
  // Compute the second perpendicular direction 
  crossProduct(1.0, normal, Xa, Xb);
  assembleFrame(Xa, Xb, normal, T);
  
  // Compute the displacement-based strain
  TacsScalar e[8], B[64*NUM_NODES];
  evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
  
  // Add the contribution from the tying straint
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
  addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

  const TacsScalar *b = B;
  for ( int ii = 0; ii < 8*NUM_NODES; ii++ ){
    sens[ii] += 
      scale*(b[0]*esens[0] + b[1]*esens[1] + b[2]*esens[2] + b[3]*esens[3] +
             b[4]*esens[4] + b[5]*esens[5] + b[6]*esens[6] + b[7]*esens[7]);
    b += 8;
  }
}

/*
  Determine the number of nodes and elements for visualization 
  generated by the data in this element. Note that higher-order
  elements are broken down into bi-linear elements for visualization.

  output:
  nelems:  the number of visualization elements
  nnodes:  the number of nodes used for visualization
  ncsr:    the number of entries in a CSR-type data structure used
  to store the connectivity
*/
void MITC9::addOutputCount( int * nelems, 
                            int * nnodes, int * ncsr ){
  *nelems += 4;
  *nnodes += 9;
  *ncsr += 16;
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is 
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below. 

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization
  
  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
void MITC9::getOutputData( unsigned int out_type,
                           double * data, int ld_data,
                           const TacsScalar Xpts[],
                           const TacsScalar vars[] ){
  for ( int m = 0, p = 0; m < 3; m++ ){
    for ( int n = 0; n < 3; n++, p++ ){
      double pt[2];
      pt[0] = -1.0 + 1.0*n;
      pt[1] = -1.0 + 1.0*m;

      TacsScalar strain[8], stress[8];
      getStrain(strain, pt, Xpts, vars);

      int index = 0;
      if (out_type & TACSElement::OUTPUT_NODES){
	for ( int k = 0; k < 3; k++ ){
	  data[index+k] = RealPart(Xpts[3*p+k]);
	}
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
	for ( int k = 0; k < NUM_DISPS; k++ ){
	  data[index+k] = RealPart(vars[NUM_DISPS*p+k]);
	}
        index += NUM_DISPS;
      }
      if (out_type & TACSElement::OUTPUT_STRAINS){
        for ( int k = 0; k < NUM_STRESSES; k++ ){
          data[index+k] = RealPart(strain[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_STRESSES){
        // Evaluate the stiffness at the current point 
        // and then calculate the stress
        stiff->calculateStress(pt, strain, stress);
        
        for ( int k = 0; k < NUM_STRESSES; k++ ){
          data[index+k] = RealPart(stress[k]);
	}
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_EXTRAS){
	// Compute the failure value
	TacsScalar lambda;
	stiff->failure(pt, strain, &lambda);
	data[index] = RealPart(lambda);

	// Compute the buckling constraint value
	TacsScalar bval;
	stiff->buckling(strain, &bval);
	data[index+1] = RealPart(bval);

	data[index+2] = RealPart(stiff->getDVOutputValue(0, pt));
	data[index+3] = RealPart(stiff->getDVOutputValue(1, pt));

        index += NUM_EXTRAS;
      }
      if (out_type & TACSElement::OUTPUT_COORDINATES){
        index += 9;
      }

      data += ld_data;
    }
  }  
}

/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or 
  less global
*/
void MITC9::getOutputConnectivity( int * con, int node ){
  int p = 0;
  for ( int m = 0; m < 2; m++ ){
    for ( int n = 0; n < 2; n++ ){
      con[4*p]   = node + n   + 3*m;
      con[4*p+1] = node + n+1 + 3*m;
      con[4*p+2] = node + n+1 + 3*(m+1);
      con[4*p+3] = node + n   + 3*(m+1);
      p++;
    }
  }
}


/*
  Test the implementation of the strain by comparing against a
  rigid-body displacement and rotation.  
*/
/*
  void MITC9::testStrain( const TacsScalar X[] ){
  TacsScalar vars[8*NUM_NODES];
  memset(vars, 0, sizeof(vars));

  const double u = 0.125;
  const double v = -0.63;

  // Pick a random set of values for the quaternions 
  TacsScalar q[4];
  q[1] = 0.25;
  q[2] = 0.5;
  q[3] = 0.125;
  q[0] = sqrt(1.0 - vecDot(&q[1], &q[1]));

  // Compute the rotation matrix C = rot - I
  TacsScalar C[9];
  C[0] =-2.0*(q[2]*q[2] + q[3]*q[3]);
  C[1] = 2.0*(q[1]*q[2] + q[3]*q[0]);
  C[2] = 2.0*(q[1]*q[3] - q[2]*q[0]);
  
  C[3] = 2.0*(q[2]*q[1] - q[3]*q[0]);
  C[4] =-2.0*(q[1]*q[1] + q[3]*q[3]);
  C[5] = 2.0*(q[2]*q[3] + q[1]*q[0]);
  
  C[6] = 2.0*(q[3]*q[1] + q[2]*q[0]);
  C[7] = 2.0*(q[3]*q[2] - q[1]*q[0]);
  C[8] =-2.0*(q[1]*q[1] + q[2]*q[2]);

  // Set the rigid displacement
  TacsScalar u0[3] = {1.25, -2.5, -4.0};

  // Compute the variables and set the quaternion values
  for ( int k = 0; k < NUM_NODES; k++ ){
  // Compute the displacements
  matMultTrans(C, &X[3*k], &vars[8*k]);
  for ( int i = 0; i < 3; i++ ){
  vars[8*k+i] += u0[i];
  }

  // Copy the values of the quaternions
  memcpy(&vars[8*k+3], q, 4*sizeof(TacsScalar));
  }

  // Compute the strain given the rigid rotation/translation
  // -------------------------------------------------------

  // Compute the reference frames at the nodes
  TacsScalar Xr[9*NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3*NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);
  
  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  TacsScalar h = inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);
  
  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);
  
  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);
  
  // Compute the rotation penalty
  TacsScalar rot = computeRotPenalty(N, Xa, Xb, Ua, Ub, vars);
  printf("rot = %15.5e\n", RealPart(rot));

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9]; 
  
  // Take the cross product to find the normal direction
  TacsScalar normal[3];
  crossProduct(1.0, Xa, Xb, normal);
  TacsScalar nrm = sqrt(vecDot(normal, normal));
  vecScale(1.0/nrm, normal);

  // Scale the Xa direction so that it is a unit vector
  nrm = sqrt(vecDot(Xa, Xa));
  vecScale(1.0/nrm, Xa);
  
  // Compute the second perpendicular direction 
  crossProduct(1.0, normal, Xa, Xb);
  
  // Assemble the transformation matrix
  assembleFrame(Xa, Xb, normal, T);

  // Compute the displacement-based strain
  TacsScalar e[8];
  evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
  
  // Add the contribution from the tying strain
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);

  for ( int k = 0; k < 8; k++ ){
  printf("strain[%d] = %15.5e\n", k, RealPart(e[k]));
  }
  
  // Compute the derivatives of the directors
  TacsScalar dirdq[12*NUM_NODES];
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the B matrix
  TacsScalar B[64*NUM_NODES];
  evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

  // Add the tying strain
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
  addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

  // Compute the derivative of the strain w.r.t.
  double dh = 1e-6;
  TacsScalar fd[8];

  for ( int k = 0; k < 8*NUM_NODES; k++ ){
  TacsScalar vtmp = vars[k];
  vars[k] = vtmp + dh;

  // Compute the directors at the nodes
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Compute the derivatives of Ua/Ub along the given directions
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);
    
  evalStrain(fd, Ur, dr, Xdinv, zXdinv, T);
  
  // Add the contribution from the tying strain
  addTyingStrain(fd, N13, N23, g13, g23, Xdinv, T);

  for ( int i = 0; i < 8; i++ ){
  fd[i] = (fd[i] - e[i])/dh;
  }

  vars[k] = vtmp;
    
  // Write out the error components
  char descript[64];
  sprintf(descript, "B%d", k);
  writeErrorComponents(stdout, descript,
  &B[8*k], fd, 8);
  }
  }
			
*/
