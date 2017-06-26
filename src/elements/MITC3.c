#include "MITC3.h"
#include "TACSElementAlgebra.h"
#include "FElibrary.h"

/*
  Rigid-body dynamics routines for TACS
  
  Copyright (c) 2015-2017 Graeme Kennedy. All rights reserved. 
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
static void writeErrorComponents( FILE *fp, const char *descript,
                                  TacsScalar *a, TacsScalar *fd, 
                                  int size, double rel_err=0.0 ){
  int print_flag = 1;
  for ( int i = 0; i < size; i++ ){
    double rel = 0.0;
    if (a[i] != 0.0){
      rel = fabs(TacsRealPart((a[i] - fd[i])/a[i]));
    }
    else {
      rel = fabs(TacsRealPart((a[i] - fd[i])));
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
                descript, i, TacsRealPart(a[i]), TacsRealPart(fd[i]), 
                fabs(TacsRealPart((a[i] - fd[i])/a[i])));
      }
      else {
        fprintf(fp, "%s[%3d] %15.6e %15.6e\n", 
                descript, i, TacsRealPart(a[i]), TacsRealPart(fd[i]));
      }
    }
  }  
}

/* 
   Return the number of displacements
*/
int MITC3::numDisplacements(){ return NUM_DISPS; }

/*
  Return the number of stresses
*/
int MITC3::numStresses(){ return NUM_STRESSES; }
  
/*
  Return the number of extras
*/
int MITC3::numExtras(){ return NUM_EXTRAS; }

/*
  Return the number of FE nodes
*/
int MITC3::numNodes(){ return NUM_NODES; }

/*
  Return the ElementType
*/
ElementType MITC3::getElementType(){ return TACS_SHELL; }

/* 
   Set up the internal static data for the names of the element,
   displacements, stresses, strains and extra variables, respectively.
*/
const char * MITC3::elemName = "MITC3";
  
const char * MITC3::dispNames[] = { "u0", "v0", "w0", 
                                    "eta", "epsx", "epsy", "epsz", "lam" };
  
const char * MITC3::stressNames[] = { "sx0", "sx1", "sy1", 
                                      "st0", "sxy0", "sxz0" };

const char * MITC3::strainNames[] = { "sx0", "sx1", "sy1", 
                                      "st0", "sxy0", "sxz0" };
  
const char * MITC3::extraNames[] = { "lambda", "buckling",
                                     "dv1", "dv2" };

/*
  Returns the elementName
*/
const char * MITC3::elementName(){ 
  return elemName;
}

/*
  Returns the displacement names
*/
const char * MITC3::displacementName( int i ){
  if (i >= 0 && i < NUM_DISPS){
    return dispNames[i];
  }
  return "";
}

/*
  Returns the name of the stresses
*/
const char * MITC3::stressName( int i ){ 
  if (i >= 0 && i < NUM_STRESSES){
    return stressNames[i];
  }
  return "";
}

/*
  Returns the name of the strains
*/
const char * MITC3::strainName( int i ){
  if (i >= 0 && i < NUM_STRESSES){
    return strainNames[i];
  }
  return "";
}

/*
  Returns the extra names
*/
const char * MITC3::extraName( int i ){
  if (i >= 0 && i < NUM_EXTRAS){
    return extraNames[i];
  }
  return "";
}

/*
  Set the design variable values
*/
void MITC3::setDesignVars( const TacsScalar dvs[], int numDVs ){
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Get the design variable values
*/
void MITC3::getDesignVars( TacsScalar dvs[], int numDVs ){
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Get the design variable range
*/
void MITC3::getDesignVarRange( TacsScalar lb[], 
                               TacsScalar ub[], int numDVs ){
  stiff->getDesignVarRange(lb, ub, numDVs);
}

/*
  Evaluate the shape functions of the element given the u-coordinate
  of the point

  input:
  u:   the parametric coordinate
  
  output:
  N:   the shape functions
  Na:   the derivative of the shape functions w.r.t. u
*/
static inline void computeShapeFunc( const double u,
                                     double N[],
                                     double Na[] ){
  // Compute the shape functions and their derivatives
  N[0] = -0.5*u*(1.0 - u);
  N[1] = (1.0 - u)*(1.0 + u);
  N[2] = 0.5*(1.0 + u)*u;

  Na[0] = -0.5 + u;
  Na[1] = -2.0*u;
  Na[2] = 0.5 + u;
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
  Xa[0] = (Na[0]*X[0] + Na[1]*X[3] + Na[2]*X[6]);
  Xa[1] = (Na[0]*X[1] + Na[1]*X[4] + Na[2]*X[7]);
  Xa[2] = (Na[0]*X[2] + Na[1]*X[5] + Na[2]*X[8]);
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
  Xa[0] = (Na[0]*X[0] + Na[1]*X[8] + Na[2]*X[16]);
  Xa[1] = (Na[0]*X[1] + Na[1]*X[9] + Na[2]*X[17]);
  Xa[2] = (Na[0]*X[2] + Na[1]*X[10] + Na[2]*X[18]);
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
  n1:   first normal to the beam curve
  n2:   second normal to the beam curve
*/
static inline void computeFrameNormals( const double N[],
                                        const TacsScalar Xr[],
                                        TacsScalar n1[],
                                        TacsScalar n2[] ){
  n1[0] = (N[0]*Xr[1] + N[1]*Xr[10] + N[2]*Xr[19]);
  n1[1] = (N[0]*Xr[4] + N[1]*Xr[13] + N[2]*Xr[22]);
  n1[2] = (N[0]*Xr[7] + N[1]*Xr[16] + N[2]*Xr[25]);

  n2[0] = (N[0]*Xr[2] + N[1]*Xr[11] + N[2]*Xr[20]);
  n2[1] = (N[0]*Xr[5] + N[1]*Xr[14] + N[2]*Xr[23]);
  n2[2] = (N[0]*Xr[8] + N[1]*Xr[17] + N[2]*Xr[26]);
}

/*
  Compute the derivative of the frame normal
*/
static inline void addFrameNormalSens( const TacsScalar n1d[],
                                       const TacsScalar n2d[],
                                       const double N[],
                                       TacsScalar Xrd[] ){
  Xrd[1] += N[0]*n1d[0];
  Xrd[4] += N[0]*n1d[1];
  Xrd[7] += N[0]*n1d[2];
  Xrd[10] += N[1]*n1d[0];
  Xrd[13] += N[1]*n1d[1];
  Xrd[16] += N[1]*n1d[2];
  Xrd[19] += N[2]*n1d[0];
  Xrd[22] += N[2]*n1d[1];
  Xrd[25] += N[2]*n1d[2];

  Xrd[2] += N[0]*n2d[0];
  Xrd[5] += N[0]*n2d[1];
  Xrd[8] += N[0]*n2d[2];
  Xrd[11] += N[1]*n2d[0];
  Xrd[14] += N[1]*n2d[1];
  Xrd[17] += N[1]*n2d[2];
  Xrd[20] += N[2]*n2d[0];
  Xrd[23] += N[2]*n2d[1];
  Xrd[26] += N[2]*n2d[2];
}

/*
  Compute the derivative of the normal direction w.r.t. the parametric
  directions

  This code computes the (approximate) through-thickness derivative of
  the inverse of the Jacobian matrix given as follows:

  d([X,r]^{-1})/dz1 and d([X,r]^{-1})/dz2

  The approximation enters through the fact that the normal direction
  is interpolated using the shape functions. This approach provides
  an approximation, but avoids the use of second derivatives. 

  input:
  Na:     the derivative of the shape functions
  Xr:     the local frames at each node
  Xdinv:  the inverse of the Jacobian transformation
  
  output:
  z1Xdinv:  the approximate first derivative w.r.t. z1
  z2Xdinv:  the approximate first derivative w.r.t. z2
*/
static inline void computeNormalRateMat( const double Na[],
                                         const TacsScalar Xr[], 
                                         const TacsScalar Xdinv[],
                                         TacsScalar z1Xdinv[],
                                         TacsScalar z2Xdinv[] ){
  // Compute the derivatives of the normal direction along the
  // parametric directions
  TacsScalar dn[9];
  dn[1] = dn[4] = dn[7] = 0.0;
  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = (Na[0]*Xr[1] + Na[1]*Xr[10] + Na[2]*Xr[19]);
  dn[3] = (Na[0]*Xr[4] + Na[1]*Xr[13] + Na[2]*Xr[22]);
  dn[6] = (Na[0]*Xr[7] + Na[1]*Xr[16] + Na[2]*Xr[25]);
  
  // Compute z1Xdinv = -Xdinv*dn*Xdinv
  TacsScalar tmp[9];
  matMatMult(dn, Xdinv, tmp);
  matMatMult(Xdinv, tmp, z1Xdinv);

  // Scale all the entries in z1Xdinv by -1
  z1Xdinv[0] = -z1Xdinv[0];
  z1Xdinv[1] = -z1Xdinv[1];
  z1Xdinv[2] = -z1Xdinv[2];
  z1Xdinv[3] = -z1Xdinv[3];
  z1Xdinv[4] = -z1Xdinv[4];
  z1Xdinv[5] = -z1Xdinv[5];
  z1Xdinv[6] = -z1Xdinv[6];
  z1Xdinv[7] = -z1Xdinv[7];
  z1Xdinv[8] = -z1Xdinv[8];

  // Compute the derivative w.r.t. z2
  dn[1] = dn[4] = dn[7] = 0.0;
  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = (Na[0]*Xr[2] + Na[1]*Xr[11] + Na[2]*Xr[20]);
  dn[3] = (Na[0]*Xr[5] + Na[1]*Xr[14] + Na[2]*Xr[23]);
  dn[6] = (Na[0]*Xr[8] + Na[1]*Xr[17] + Na[2]*Xr[26]);
  
  // Compute z2Xdinv = -Xdinv*dn*Xdinv
  matMatMult(dn, Xdinv, tmp);
  matMatMult(Xdinv, tmp, z2Xdinv);

  // Scale all the entries in z2Xdinv by -1
  z2Xdinv[0] = -z2Xdinv[0];
  z2Xdinv[1] = -z2Xdinv[1];
  z2Xdinv[2] = -z2Xdinv[2];
  z2Xdinv[3] = -z2Xdinv[3];
  z2Xdinv[4] = -z2Xdinv[4];
  z2Xdinv[5] = -z2Xdinv[5];
  z2Xdinv[6] = -z2Xdinv[6];
  z2Xdinv[7] = -z2Xdinv[7];
  z2Xdinv[8] = -z2Xdinv[8];
}

/*
  Compute the derivative of the matrix zXdinv with respect to the
  inputs. The matrix zXdinv is computed using the following expression
  in index notation:

  zXdinv_{kl} = -Xdinv_{kn}*dn_{nm}*Xdinv_{ml}

  There are two main computations:
  1. The evaluation of dn (the derivative of the normal direction)
  2. The evaluation of zXdinv once dn is known

  The contributions to the sensitivity of dn are:
  
  d(zXdinv_{kl})/d(dn{ij})
  .  = -Xdinv_{kn}*delta_{ni}*delta_{mj}*Xdinv_{ml}
  .  = -Xdinv_{ki}*Xdinv_{jl}

  dnd_{ij} = zXdinvd_{kl}*d(zXdinv_{kl})/d(d_{ij})
  .  = -zXdinvd_{kl}*Xdinv_{ki}*Xdinv_{jl}
  
  dnd = -Xdinv^{T}*zXdinvd*Xdinv^{T}
  
  The second contribution is:

  d(zXdinv_{kl})/d(Xdinv_{ij})
  .  = -Xdinv_{kn}*dn_{nm}*delta_{mi}*delta_{lj} 
  .    -delta_{ki}*delta_{nj}*d_{nm}*Xdinv_{ml}

  Xdinvd_{ij} = zXdinvd_{kl}*d(zXdinv_{kl})/d(Xdinv_{ij})
  .  = -Xdinv_{kn}*dn_{ni}*zXdinvd_{kj}
  .    -zXdinvd_{il}*dn_{jm}*Xdinv_{ml}

  dnd += -(Xdinv*dn)^{T}*Xdinvd -Xdinvd*(dn*Xdinv)^{T}
*/
/*
static inline void addNormalRateMatSens( TacsScalar Xrd[],
                                         TacsScalar Xdinvd[],
                                         const TacsScalar zXdinvd[],
                                         const double Na[],
                                         const double Nb[],
                                         const TacsScalar Xr[], 
                                         const TacsScalar Xdinv[] ){
  // Compute the derivatives of the normal direction along the
  // parametric directions
  TacsScalar dn[9];
  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = (Na[0]*Xr[2] + Na[1]*Xr[11] + Na[2]*Xr[20]);
  dn[3] = (Na[0]*Xr[5] + Na[1]*Xr[14] + Na[2]*Xr[23]);
  dn[6] = (Na[0]*Xr[8] + Na[1]*Xr[17] + Na[2]*Xr[26]);
  
  dn[1] = (Nb[0]*Xr[2] + Nb[1]*Xr[11] + Nb[2]*Xr[20]);
  dn[4] = (Nb[0]*Xr[5] + Nb[1]*Xr[14] + Nb[2]*Xr[23]);
  dn[7] = (Nb[0]*Xr[8] + Nb[1]*Xr[17] + Nb[2]*Xr[26]);
  
  // Compute the contribution
  // -dnd = Xdinv^{T}*zXdinvd*Xdinv^{T}
  TacsScalar dnd[9], t[9];
  matMatTransMult(zXdinvd, Xdinv, t);
  matTransMatMult(Xdinv, t, dnd);

  // Add the contribution to Xdinvd
  // Xdinvd -= (Xdinv*dn)^{T}*Xdinvd + Xdinvd*(dn*Xdinv)^{T}
  dn[0] = -dn[0];
  dn[3] = -dn[3];
  dn[6] = -dn[6];

  dn[1] = -dn[1];
  dn[4] = -dn[4];
  dn[7] = -dn[7];
  matMatMult(Xdinv, dn, t);
  matTransMatMultAdd(t, zXdinvd, Xdinvd);
  matMatMult(dn, Xdinv, t);
  matMatTransMultAdd(zXdinvd, t, Xdinvd);

  // Compute the contribution to Xrd
  for ( int i = 0; i < 3; i++ ){
    Xrd[2+9*i] -= (dnd[0]*Na[i] + dnd[1]*Nb[i]);
    Xrd[5+9*i] -= (dnd[3]*Na[i] + dnd[4]*Nb[i]);
    Xrd[8+9*i] -= (dnd[6]*Na[i] + dnd[7]*Nb[i]);
  }
}
*/

/*
  Compute the tying shape functions for the element for the
  out-of-plane shear components.
  
  This uses the quadrature points from a two-point Gauss quadrature
  scheme as the tying points. The tying point shape functions for the
  g13 and g23 shear strains are numbered as follows:

  g13               g23
  +--0-----1--+     +--0-----1--+

  input:
  u:    the parametric coordinate
  
  output:
  N13:  the shape functions associated with the g13 shear strain
*/
static inline void computeTyingFunc( const double u,
                                     double N13[] ){
  // The tying point offset
  const double t = 0.577350269189626;
  const double tinv = 1.0/t;

  // Compute the shape functions for the reduced dimension
  N13[0] = 0.5*tinv*(t - u);
  N13[1] = 0.5*tinv*(t + u);
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
  Constructor for the MITC3 element class

  input:
  stiff:      the stiffness object
  gravity:    the gravity vector
  vInit:      the initial velocity
  omegaInit:  the initial angular velocity
*/
MITC3::MITC3( TimoshenkoStiffness *_stiff, 
              TACSGibbsVector *_gravity,
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

MITC3::~MITC3(){
  stiff->decref();
  if (gravity){ gravity->decref(); }
  if (vInit){ vInit->decref(); }
  if (omegaInit){ omegaInit->decref(); }
}

/*
  Retrieve the initial values of the design variables
*/
void MITC3::getInitConditions( TacsScalar vars[], 
                               TacsScalar dvars[],
                               TacsScalar ddvars[],
                               const TacsScalar X[] ){
  memset(vars, 0, 8*NUM_NODES*sizeof(TacsScalar));
  memset(dvars, 0, 8*NUM_NODES*sizeof(TacsScalar));
  memset(ddvars, 0, 8*NUM_NODES*sizeof(TacsScalar));

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
      // dot{u} = v + omega^{x}*r
      crossProductAdd(1.0, omega, &X[3*i], &dvars[8*i]);
      
      // d{eps}/dt = 0.5*omega
      dvars[8*i+4] = 0.5*omega[0];
      dvars[8*i+5] = 0.5*omega[1];
      dvars[8*i+6] = 0.5*omega[2];

      // ddot{u} = omega^{x}*omega^{x}*r
      // Note that the second derivative of the quaternions is zero since
      // there is no angular acceleration
      TacsScalar omegar[3];
      crossProduct(1.0, omega, &X[3*i], omegar);
      crossProduct(1.0, omega, omegar, &ddvars[8*i]);
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
void MITC3::computeEnergies( double time,
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
  TacsScalar d1[3*NUM_NODES], d2[3*NUM_NODES];
  computeDirectors(d1, d2, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[2], g23[2];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Initialize the velocities
  TacsScalar Te = 0.0, Pe = 0.0;

  // Evaluate the kinetic energy of the element
  for ( int i = 0; i < ORDER; i++ ){
      // Set the Gauss quadrature points
    const double u = gaussPts[i];
  
    // Evaluate the shape functions
    double N[NUM_NODES], Na[NUM_NODES];
    computeShapeFunc(u, N, Na);

    // Compute the derivative along the shape function direction
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);
      
    // Compute the frame normal
    TacsScalar n1[3], n2[3];
    computeFrameNormals(N, Xr, n1, n2);

    // Evaluate the derivatives in the locally-aligned frame
    TacsScalar Xd[9], Xdinv[9];
    assembleFrame(Xa, n1, n2, Xd);

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
    computeTransform(T, Xa, Xb);
      
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
void MITC3::addResidual( double time,
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

  // Compute the area for this element: this is used to scale the
  // constraint equations within the code
  TacsScalar area = 0.0;

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

      // Compute the area
      area += h;

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
      computeTransform(T, Xa, Xb);

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

        r += 8;
        b += 64;
        br += 8;
      }
    }
  }
  
  // Set the scaling for the constraints
  TacsScalar scale = area;

  // Add the constraints from the quaternion parametrization
  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i+3];
    TacsScalar lamb = vars[8*i+7];

    // Add the result to the governing equations
    res[8*i+3] += 2.0*scale*q[0]*lamb;
    res[8*i+4] += 2.0*scale*q[1]*lamb;
    res[8*i+5] += 2.0*scale*q[2]*lamb;
    res[8*i+6] += 2.0*scale*q[3]*lamb;
    
    // Enforce the quaternion constraint
    res[8*i+7] += scale*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] - 1.0);
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
void MITC3::addJacobian( double time, TacsScalar J[],
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

  // Compute the area of this element: used to evaluate the
  // constraints
  TacsScalar area = 0.0;

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

      // Compute the area
      area += h;

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
        computeTransform(T, Xa, Xb);
        
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

  // Set the scaling for the constraints
  TacsScalar scale = 2.0*alpha*area;

  // Add the constraints from the quaternions
  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i+3];
    const int ldj = 8*NUM_NODES;

    TacsScalar *Jp = &J[(8*NUM_NODES+1)*(8*i + 3)];
    TacsScalar lamb = vars[8*i+7];
    
    // Add the constraint terms
    Jp[4] += scale*q[0];
    Jp[4+ldj] += scale*q[1];
    Jp[4+2*ldj] += scale*q[2];
    Jp[4+3*ldj] += scale*q[3];

    // Enforce the quaternion constraint
    Jp[4*ldj] += scale*q[0];
    Jp[4*ldj+1] += scale*q[1];
    Jp[4*ldj+2] += scale*q[2];
    Jp[4*ldj+3] += scale*q[3];

    // Add the terms to the diagonal
    Jp[0] += scale*lamb;
    Jp[ldj+1] += scale*lamb;
    Jp[2*(ldj+1)] += scale*lamb;
    Jp[3*(ldj+1)] += scale*lamb;
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
void MITC3::addAdjResProduct( double time, double scale,
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
      computeTransform(T, Xa, Xb);

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
  Compute the derivative of the product of the adjoint variables and
  the residual vector.

  This code uses a reverse-mode method, which makes the derivative
  code faster, but a bit more difficult to follow.
  
  input:
  time:      the simulation time
  scale:     scale the derivative vector by this factor
  psi:       the adjoint variables associated with this element
  vars:      the state variables
  dvars:     the first derivative of the state variables
  ddvars:    the second derivative of the state variables

  output:
  fXptSens:  add the derivative to this vector
*/
void MITC3::addAdjResXptProduct( double time, double scale,
                                 TacsScalar fXptSens[],
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

  // The derivatives w.r.t. reference frames
  TacsScalar Xrd[9*NUM_NODES];
  memset(Xrd, 0, 9*NUM_NODES*sizeof(TacsScalar));

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3*NUM_NODES], domega[3*NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3*NUM_NODES], dirdq[12*NUM_NODES];  
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Memory to keep track of the derivatives w.r.t. directors
  TacsScalar dird[3*NUM_NODES], dirdqd[12*NUM_NODES];  
  memset(dird, 0, 3*NUM_NODES*sizeof(TacsScalar));
  memset(dirdqd, 0, 12*NUM_NODES*sizeof(TacsScalar));

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6*8*NUM_NODES], B23[6*8*NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Memory to keep track of the derivative contributions
  TacsScalar g13d[6], g23d[6];
  TacsScalar B13d[6*8*NUM_NODES], B23d[6*8*NUM_NODES];
  memset(g13d, 0, 6*sizeof(TacsScalar));
  memset(g23d, 0, 6*sizeof(TacsScalar));
  memset(B13d, 0, 6*8*NUM_NODES*sizeof(TacsScalar));
  memset(B23d, 0, 6*8*NUM_NODES*sizeof(TacsScalar));

  // Add up the contributions to the derivative of the area
  TacsScalar aread = 0.0;

  // Add the constraints from the quaternion parametrization
  for ( int i = 0; i < NUM_NODES; i++ ){
    const TacsScalar *q = &vars[8*i+3];
    TacsScalar lamb = vars[8*i+7];
    const TacsScalar *p = &psi[8*i+3];
    
    // Enforce the quaternion constraint
    aread += (p[7]*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] - 1.0) +
              2.0*lamb*(psi[0]*q[0] + psi[1]*q[1] + psi[2]*q[2] + psi[3]*q[3]));
  }

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

      // Add the contribution from the derivative of the inner product
      // w.r.t. h (the determinant of the Jacobian matrix)
      TacsScalar hd = aread;

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
      TacsScalar w[3], dw[3];
      TacsScalar tmp1 = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp1*fn[0];
      w[1] = omeg[1] - tmp1*fn[1];
      w[2] = omeg[2] - tmp1*fn[2];

      TacsScalar tmp2 = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp2*fn[0];
      dw[1] = domeg[1] - tmp2*fn[1];
      dw[2] = domeg[2] - tmp2*fn[2];

      // Accumulate the derivative contributions
      TacsScalar wd[3], dwd[3];
      wd[0] = wd[1] = wd[2] = 0.0;
      dwd[0] = dwd[1] = dwd[2] = 0.0;

      // Add the contribution to the residual
      const TacsScalar *p = psi;
      const TacsScalar *q = vars, *dq = dvars;
      for ( int ii = 0; ii < NUM_NODES; ii++ ){
        // Add the contributions from the rectilinear velocity
        hd += N[ii]*rho[0]*(p[0]*a0[0] + p[1]*a0[1] + p[2]*a0[2]);

        // Add the contributions from the angular velocity
        // S^{T}*dw + 2*dot{S}^{T}*w
        TacsScalar eta = q[3];
        const TacsScalar *eps = &q[4];
        TacsScalar deta = dq[3];
        const TacsScalar *deps = &dq[4];

        // Add S^{T}*dw
        TacsScalar scl = 2.0*N[ii]*rho[1];
        TacsScalar t[3];
        crossProduct(1.0, eps, dw, t);
        hd -= scl*p[3]*vecDot(eps, dw);
        hd += scl*vecDot(t, &p[4]);
        hd += scl*eta*vecDot(dw, &p[4]);

        // Add 2*dot{S}^{T}*w
        crossProduct(1.0, deps, w, t);
        hd -= 2.0*scl*p[3]*vecDot(deps, w);
        hd += 2.0*scl*vecDot(t, &p[4]); 
        hd += 2.0*scl*deta*vecDot(w, &p[4]);

        // Add the contributions to wd and dwd
        scl *= h;
        vecAxpy(-scl*p[3], eps, dwd);
        crossProductAdd(scl, &p[4], eps, dwd);
        vecAxpy(scl*eta, &p[4], dwd);

        vecAxpy(-2.0*scl*p[3], eps, wd);
        crossProductAdd(2.0*scl, &p[4], eps, wd);
        vecAxpy(2.0*scl*eta, &p[4], wd);

        p += 8;
        q += 8;
        dq += 8;
      }

      TacsScalar fnd[3];
      TacsScalar tmpw = vecDot(wd, fn);
      TacsScalar tmpdw = vecDot(dwd, fn);
      fnd[0] = -(tmp1*wd[0] + tmp2*dwd[0] + tmpw*omeg[0] + tmpdw*domeg[0]);
      fnd[1] = -(tmp1*wd[1] + tmp2*dwd[1] + tmpw*omeg[1] + tmpdw*domeg[1]);
      fnd[2] = -(tmp1*wd[2] + tmp2*dwd[2] + tmpw*omeg[2] + tmpdw*domeg[2]);

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

      // Add the terms from the potential energy
      TacsScalar Brot[8*NUM_NODES];
      TacsScalar rot = computeBRotPenalty(Brot, N, Na, Nb,
                                          Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9]; 
      computeTransform(T, Xa, Xb);

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
      
      // Compute the product Bpsi = B*psi
      TacsScalar Bpsi[NUM_STRESSES];
      memset(Bpsi, 0, NUM_STRESSES*sizeof(TacsScalar));

      // Compute the product BRpsi = Brot*psi
      TacsScalar BRpsi = 0.0;

      // Set the pointer to the psi variables back to the begining
      p = psi;

      // Compute the products
      const TacsScalar *b = B;
      const TacsScalar *br = Brot;
      for ( int k = 0; k < 8*NUM_NODES; k++ ){
        Bpsi[0] += b[0]*p[0];
        Bpsi[1] += b[1]*p[0];
        Bpsi[2] += b[2]*p[0];
        Bpsi[3] += b[3]*p[0];
        Bpsi[4] += b[4]*p[0];
        Bpsi[5] += b[5]*p[0];
        Bpsi[6] += b[6]*p[0];
        Bpsi[7] += b[7]*p[0];
        BRpsi += br[0]*p[0];
        hd += p[0]*(strainProduct(s, b) + rot*br[0]);        
        p++;
        br++;
        b += NUM_STRESSES;
      }

      // Add the contribution from the gravity vector
      for ( int k = 0; k < NUM_NODES; k++ ){
        hd -= rho[0]*N[k]*(g[0]*psi[8*k] + g[1]*psi[8*k+1] + g[2]*psi[8*k+2]);
      }

      // Scale by the derivative w.r.t. h by the quadrature weight
      hd *= gaussWts[i]*gaussWts[j]*scale;

      // Compute the product spsi = C*B*psi
      TacsScalar spsi[NUM_STRESSES];
      stiff->calculateStress(A, Bc, D, As, Bpsi, spsi);

      // Add the contribution from the derivative of the
      // term h*e^{T}*C*B*psi with respect to the stress     
      TacsScalar Urd[9], drd[9];
      TacsScalar Xdinvd[9], zXdinvd[9];
      TacsScalar Td[9];
      evalStrainSens(Urd, drd, Xdinvd, zXdinvd, Td,
                     h*scale, spsi, Ur, dr, Xdinv, zXdinv, T);
      addTyingStrainSens(g13d, g23d, Xdinvd, Td,
                         h*scale, spsi, N13, N23, g13, g23, Xdinv, T);

      // Compute the quantity hs = h*s;
      TacsScalar hs[NUM_STRESSES];
      for ( int k = 0; k < NUM_STRESSES; k++ ){
        hs[k] = s[k]*h*scale;
      }      

      // Add the contribution from the derivative of the 
      // term h*s^{T}*B*psi w.r.t. B
      addBmatSens(Urd, drd, Xdinvd, zXdinvd, Td, dirdqd,
                  hs, psi, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
      addTyingBmatSens(B13d, B23d, Xdinvd, Td,
                       hs, psi, N13, N23, B13, B23, Xdinv, T);

      // Add the contribution from the director computation
      for ( int k = 0; k < NUM_NODES; k++ ){
        dird[3*k] += Na[k]*drd[0] + Nb[k]*drd[1] + N[k]*Urd[2];
        dird[3*k+1] += Na[k]*drd[3] + Nb[k]*drd[4] + N[k]*Urd[5];
        dird[3*k+2] += Na[k]*drd[6] + Nb[k]*drd[7] + N[k]*Urd[8];
      }

      // Compute the derivatives Xad and Xbd
      TacsScalar Xad[3], Xbd[3];
      computeTransformSens(Xad, Xbd, Td, Xa, Xb);

      // Add the contribution from the rotation
      TacsScalar rotd = scale*h*BRpsi*kpenalty;
      addBRotPenaltySens(Xad, Xbd, rotd, scale*h*rot, psi,
                         N, Na, Nb, Xa, Xb, Ua, Ub, vars);

      // Compute the derivatives Xrd
      addNormalRateMatSens(Xrd, Xdinvd, zXdinvd,
                           Na, Nb, Xr, Xdinv);

      // Compute the derivative of the inverse of the Jacobian 
      TacsScalar Xdd[9];
      inv3x3Sens(Xdd, Xdinvd, Xdinv);

      // Add the contributions from the derivative of the determinant
      TacsScalar hXdd[9];
      det3x3Sens(Xd, hXdd);
      for ( int k = 0; k < 9; k++ ){
        Xdd[k] += hXdd[k]*hd;
      }

      // Extract/add the sensitivities from the frame
      fnd[0] += Xdd[2];  fnd[1] += Xdd[5];  fnd[2] += Xdd[8];
  
      // Add the contributions to Xad and Xbd
      Xad[0] += Xdd[0];  Xad[1] += Xdd[3];  Xad[2] += Xdd[6];
      Xbd[0] += Xdd[1];  Xbd[1] += Xdd[4];  Xbd[2] += Xdd[7];
  
      // // Compute the frame normal
      addFrameNormalSens(fnd, N, Xrd);

      // Add the derivatives the shape function directions
      for ( int k = 0; k < NUM_NODES; k++ ){
        fXptSens[3*k] += Na[k]*Xad[0] + Nb[k]*Xbd[0];
        fXptSens[3*k+1] += Na[k]*Xad[1] + Nb[k]*Xbd[1];
        fXptSens[3*k+2] += Na[k]*Xad[2] + Nb[k]*Xbd[2];
      }
    }
  }

  addComputeTyingBmatSens(fXptSens, Xrd, dird, dirdqd, 
                          g13d, g23d, B13d, B23d, X, Xr, vars, dir, dirdq);
  addDirectorsSens(Xrd, dird, vars);
  addDirectorDerivSens(Xrd, dirdqd, vars);
  addFramesSens(fXptSens, Xrd, X);
}

/*
  Given the nodal degrees of freedom and their time-derivatives,
  compute the angular velocity at each node
*/
void MITC3::computeAngularVelocity( TacsScalar omega[],
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
void MITC3::computeAngularAccel( TacsScalar domega[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Set pointers to the values 
    TacsScalar eta = vars[3];
    const TacsScalar *eps = &vars[4];
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
  At each node in the finite-element, compute the derivatives of the
  coordinate directions and assemble a locally-aligned reference
  frame.

  Each locally-aligned reference frame Xr[] consists of a 3x3 matrix
  stored in row-major order where the first direction is aligned along
  the axial direction, the second direction is normal to the axial
  direction along th


1-direction and the third direction is normal to the local
  surface and the second direction is perpendicular to both these
  directions. This can be written as follows:

  Xr = [X,xi1; X,xi2; n]

  input:
  X:    the initial nodal locations
  
  output:
  Xr:   the locally-aligned frames
*/
void MITC3::computeFrames( TacsScalar Xr[],
                           const TacsScalar X[] ){
  // Get the reference axis associated with the beam
  const TacsScalar *axis;
  stiff->getRefAxis(&axis);  

  for ( int i = 0; i < ORDER; i++ ){
    // Find the u/v values at the node locations
    double u = -1.0 + 2.0*i/(ORDER-1.0);
    
    // Evaluate the shape functions
    double Na[NUM_NODES], Nb[NUM_NODES];
    computeShapeFunc(u, Na, Nb);

    // Compute the derivative along the shape function directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the cross product
    TacsScalar n1[3], n2[3];

    
    crossProduct(1.0, Xa, Xb, normal);
    TacsScalar nrm = sqrt(vecDot(normal, normal));
    vecScale(1.0/nrm, normal);

    // Assemble the frame into the Xr matrix
    assembleFrame(Xa, n1, n2, Xr);
    
    // Increment the pointer to the frames
    Xr += 9;  
  }
}

/*
  Given the non-orthonormal in-plane directions Xa/Xb which are the
  parametric tangents, compute the transformation to the
  locally-aligned frame for the strain.

  input:
  Xa:   tangent to the first parametric direction
  Xb:   tangent to the second parametric direction

  output:
  T:    the transformation matrix
*/
void MITC3::computeTransform( TacsScalar T[],
                              const TacsScalar Xa[],
                              const TacsScalar Xb[] ){      
  // Compute the cross product of Xa/Xb to find the normal direction
  TacsScalar normal[3];
  crossProduct(1.0, Xa, Xb, normal);
  TacsScalar nrm = sqrt(vecDot(normal, normal));
  vecScale(1.0/nrm, normal);

  // The first coordinate direction
  TacsScalar d1[3];
  if (stiff->getTransformType() == FSDTStiffness::NATURAL){
    // Scale the Xa direction so that it is a unit vector
    TacsScalar invXa = 1.0/sqrt(vecDot(Xa, Xa));
    d1[0] = invXa*Xa[0];
    d1[1] = invXa*Xa[1];
    d1[2] = invXa*Xa[2];
  }
  else {
    // Retrieve the reference axis from the constitutive
    // object
    const TacsScalar *axis = stiff->getRefAxis();
    TacsScalar an = vecDot(axis, normal);
    
    // Take the component of the reference axis perpendicular
    // to the surface
    d1[0] = axis[0] - an*normal[0];
    d1[1] = axis[1] - an*normal[1];
    d1[2] = axis[2] - an*normal[2];
    
    // Normalize the new direction
    TacsScalar inv = 1.0/sqrt(vecDot(d1, d1));
    vecScale(inv, d1);
  }

  // Compute the second perpendicular direction 
  TacsScalar d2[3];
  crossProduct(1.0, normal, d1, d2);
  
  // Assemble the transformation matrix
  assembleFrame(d1, d2, normal, T);
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
void MITC3::computeDirectors( TacsScalar d1[],
                              TacsScalar d2[],
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

    // Compute d1 = C^{T}*n1
    d1[0] = C[0]*Xr[1] + C[3]*Xr[4] + C[6]*Xr[7];
    d1[1] = C[1]*Xr[1] + C[4]*Xr[4] + C[7]*Xr[7];
    d1[2] = C[2]*Xr[1] + C[5]*Xr[4] + C[8]*Xr[7];

    // Compute d2 = C^{T}*n2
    d2[0] = C[0]*Xr[2] + C[3]*Xr[5] + C[6]*Xr[8];
    d2[1] = C[1]*Xr[2] + C[4]*Xr[5] + C[7]*Xr[8];
    d2[2] = C[2]*Xr[2] + C[5]*Xr[5] + C[8]*Xr[8];

    d1 += 3; // Each director is a 3-vector
    d2 += 3; // Each director is a 3-vector
    Xr += 9; // Increment over each frame
    vars += 8; // 8 variables per node
  }
}

/*
  Add the derivative of the directors to the local frame
*/
void MITC3::addDirectorsSens( TacsScalar Xrd[],
                              const TacsScalar d1d[],
                              const TacsScalar d2d[],
                              const TacsScalar vars[] ){
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

    Xrd[1] += C[0]*d1d[0] + C[1]*d1d[1] + C[2]*d1d[2];
    Xrd[4] += C[3]*d1d[0] + C[4]*d1d[1] + C[5]*d1d[2];
    Xrd[7] += C[6]*d1d[0] + C[7]*d1d[1] + C[8]*d1d[2];

    Xrd[2] += C[0]*d2d[0] + C[1]*d2d[1] + C[2]*d2d[2];
    Xrd[5] += C[3]*d2d[0] + C[4]*d2d[1] + C[5]*d2d[2];
    Xrd[8] += C[6]*d2d[0] + C[7]*d2d[1] + C[8]*d2d[2];

    d1d += 3; // Each director is a 3-vector
    d2d += 3; // Each director is a 3-vector
    Xrd += 9; // Increment over each frame
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
void MITC3::computeDirectorDeriv( TacsScalar ddq[],
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
void MITC3::evalStrain( TacsScalar e[], 
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

  This function is used in the evaluation of the residuals and in the
  computation of the Jacobian matrix.

  input:
  N, Na, Nb:   shape functions and their derivatives
  Ur:          derivatives of displacements in the parametric direction
  dr:          derivative of the director field
  Xdinv:       Jacobian transformation
  zXdinv:      through-thickness derivative of the Jacobian transformation
  T:           transformation to the local strain axis
  dirdq:       derivative of the director field with respect to states

  output:
  e:           the strain
  B;           the derivative of the strain w.r.t. the variables
*/
void MITC3::evalBmat( TacsScalar e[],
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
void MITC3::addGmat( TacsScalar J[],
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

  // Add the contributions from the second derivatives of the
  // quaternions
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
void MITC3::computeTyingStrain( TacsScalar g13[], 
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
void MITC3::computeTyingBmat( TacsScalar g13[], 
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
void MITC3::addTyingGmat( TacsScalar J[],
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
void MITC3::addTyingStrain( TacsScalar e[],
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
void MITC3::addTyingGmatWeights( TacsScalar w13[],
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
void MITC3::addTyingBmat( TacsScalar B[],
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
  Add the sensitivities from the tying strain components
*/
void MITC3::addTyingBmatSens( TacsScalar B13d[], TacsScalar B23d[],
                              TacsScalar Xdinvd[], TacsScalar Td[],
                              const TacsScalar eSens[],
                              const TacsScalar psi[],
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

  TacsScalar A11d = 0.0, A12d = 0.0;
  TacsScalar A21d = 0.0, A22d = 0.0;

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
    TacsScalar G13d = 2.0*psi[0]*(A12*eSens[6] + A11*eSens[7]);
    TacsScalar G23d = 2.0*psi[0]*(A22*eSens[6] + A21*eSens[7]);
  
    A11d += 2.0*psi[0]*G13*eSens[7];
    A12d += 2.0*psi[0]*G13*eSens[6];
    A21d += 2.0*psi[0]*G23*eSens[7];
    A22d += 2.0*psi[0]*G23*eSens[6];

    for ( int i = 0; i < 6; i++ ){
      B13d[offset*i] += G13d*N13[i];
      B23d[offset*i] += G23d*N23[i];
    }

    B13++;   B23++;
    B13d++;  B23d++;
    psi++;
  }

  // Add the derivative to the transformation terms
  Td[0] += Xdinv[0]*A11d + Xdinv[3]*A21d;
  Td[1] += Xdinv[0]*A12d + Xdinv[3]*A22d;
  Td[3] += Xdinv[1]*A11d + Xdinv[4]*A21d;
  Td[4] += Xdinv[1]*A12d + Xdinv[4]*A22d;
  Td[6] += Xdinv[2]*A11d + Xdinv[5]*A21d;
  Td[7] += Xdinv[2]*A12d + Xdinv[5]*A22d;

  // Add the terms to the Xdinv
  Xdinvd[0] += T[0]*A11d + T[1]*A12d;
  Xdinvd[1] += T[3]*A11d + T[4]*A12d;
  Xdinvd[2] += T[6]*A11d + T[7]*A12d;

  Xdinvd[3] += T[0]*A21d + T[1]*A22d;
  Xdinvd[4] += T[3]*A21d + T[4]*A22d;
  Xdinvd[5] += T[6]*A21d + T[7]*A22d;
}

/*
  Get the constitutive object
*/
TACSConstitutive *MITC3::getConstitutive(){
  return stiff;
}

/*
  Return the number of quadrature points
*/
int MITC3::getNumGaussPts(){ 
  return ORDER*ORDER;
}

/*
  Return the quadrature points and weights
*/
double MITC3::getGaussWtsPts( const int num, double pt[] ){
  int m = (int)(num/ORDER);
  int n = num % ORDER;
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];    
  
  return gaussWts[n]*gaussWts[m];
}

/*
  Get the values of the shape functions
*/
void MITC3::getShapeFunctions( const double pt[], double N[] ){
  double Na[3];
  computeShapeFunc(pt[0], N, Na);
}

/*
  Retrieve the determinant of the Jacobian transformation matrix
*/
TacsScalar MITC3::getDetJacobian( const double pt[],
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
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the determinant of the Jacobian transformation
  return det3x3(Xd);
}

/*
  Evaluate the determinant of the Jacobian transformation w.r.t.
  the node locations
*/
TacsScalar MITC3::getDetJacobianXptSens( TacsScalar *hXptSens, 
                                         const double pt[], 
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
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the determinant of the Jacobian transformation
  TacsScalar h = det3x3(Xd);

  // Compute the derivative of the determinant
  TacsScalar Xdd[9];
  det3x3Sens(Xd, Xdd);

  // Extract/add the sensitivities from the frame
  TacsScalar fnd[3];
  fnd[0] = Xdd[2];  fnd[1] = Xdd[5];  fnd[2] = Xdd[8];
  
  // Add the contributions to Xad and Xbd
  TacsScalar Xad[3], Xbd[3];
  Xad[0] = Xdd[0];  Xad[1] = Xdd[3];  Xad[2] = Xdd[6];
  Xbd[0] = Xdd[1];  Xbd[1] = Xdd[4];  Xbd[2] = Xdd[7];
  
  // // Compute the frame normal
  TacsScalar Xrd[9*NUM_NODES];
  memset(Xrd, 0, 9*NUM_NODES*sizeof(TacsScalar));
  addFrameNormalSens(fnd, N, Xrd);

  // Add the derivatives the shape function directions
  for ( int k = 0; k < NUM_NODES; k++ ){
    hXptSens[3*k] = Na[k]*Xad[0] + Nb[k]*Xbd[0];
    hXptSens[3*k+1] = Na[k]*Xad[1] + Nb[k]*Xbd[1];
    hXptSens[3*k+2] = Na[k]*Xad[2] + Nb[k]*Xbd[2];
  }
  
  addFramesSens(hXptSens, Xrd, X);

  return h;
}

/*
  Evaluate the strain at a parametric point within the element
*/
void MITC3::getStrain( TacsScalar e[],
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
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);
  
  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9];
  computeTransform(T, Xa, Xb);
    
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
void MITC3::addStrainSVSens( TacsScalar sens[],
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
  computeTransform(T, Xa, Xb);
  
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
void MITC3::addOutputCount( int * nelems, 
                            int * nnodes, int * ncsr ){
  *nelems += 2;
  *nnodes += 3;
  *ncsr += 4;
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
void MITC3::getOutputData( unsigned int out_type,
                           double * data, int ld_data,
                           const TacsScalar Xpts[],
                           const TacsScalar vars[] ){
  for ( int n = 0; n < 3; n++, p++ ){
    double pt[2];
    pt[0] = -1.0 + 1.0*n;

    TacsScalar strain[6], stress[6];
    getStrain(strain, pt, Xpts, vars);
      
    int index = 0;
    if (out_type & TACSElement::OUTPUT_NODES){
      for ( int k = 0; k < 3; k++ ){
        data[index+k] = TacsRealPart(Xpts[3*p+k]);
      }
      index += 3;
    }
    if (out_type & TACSElement::OUTPUT_DISPLACEMENTS){
      for ( int k = 0; k < NUM_DISPS; k++ ){
        data[index+k] = TacsRealPart(vars[NUM_DISPS*p+k]);
      }
      index += NUM_DISPS;
    }
    if (out_type & TACSElement::OUTPUT_STRAINS){
      // Add the term due to the potential energy
      for ( int k = 0; k < NUM_STRESSES; k++ ){
        data[index+k] = TacsRealPart(strain[k]);
      }
      index += NUM_STRESSES;
    }
    if (out_type & TACSElement::OUTPUT_STRESSES){
      // Evaluate the stiffness at the current point 
      // and then calculate the stress
      stiff->calculateStress(pt, strain, stress);
      
      for ( int k = 0; k < NUM_STRESSES; k++ ){
        data[index+k] = TacsRealPart(stress[k]);
      }
      index += NUM_STRESSES;
    }
    if (out_type & TACSElement::OUTPUT_EXTRAS){
      // Compute the failure value
      TacsScalar lambda;
      stiff->failure(pt, strain, &lambda);
      data[index] = TacsRealPart(lambda);
      
      // Compute the buckling constraint value
      TacsScalar bval;
      stiff->buckling(strain, &bval);
      data[index+1] = TacsRealPart(bval);
      
      data[index+2] = TacsRealPart(stiff->getDVOutputValue(0, pt));
      data[index+3] = TacsRealPart(stiff->getDVOutputValue(1, pt));
      
      index += NUM_EXTRAS;
    }
    if (out_type & TACSElement::OUTPUT_COORDINATES){
      index += 9;
    }
    
    data += ld_data;    
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
void MITC3::getOutputConnectivity( int * con, int node ){
  con[0] = node;
  con[1] = node+1;
  con[2] = node+1;
  con[3] = node+2;
}
