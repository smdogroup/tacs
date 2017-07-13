#include "RigidBody.h"
#include "KinematicConstraints.h"
#include "TACSElementAlgebra.h"
#include "FElibrary.h"

/*
  Construct a spherical constraint with the two bodies involved and a
  position vector measured from the global frame to the point where
  the spherical joint is located.
*/
TACSSphericalConstraint::TACSSphericalConstraint( TACSRigidBody *_bodyA, 
                                                  TACSRigidBody *_bodyB, 
                                                  TACSGibbsVector *_point ){
  // Copy over the arguments
  bodyA = _bodyA; bodyA->incref();
  bodyB = _bodyB; bodyB->incref();
  point = _point; point->incref();

  // Set class variables to NULL
  xAVec = NULL;
  xBVec = NULL;
  
  updatePoints();
}


/*
  Construct a spherical constraint with a body involved and a position
  vector measured from the global frame to the point where the
  spherical joint is located.
*/
TACSSphericalConstraint::TACSSphericalConstraint( TACSRigidBody *_bodyA, 
                                                  TACSGibbsVector *_point ){
  // Copy over the arguments
  bodyA = _bodyA; bodyA->incref();
  bodyB = NULL;
  point = _point; point->incref();

  // Set class variables to NULL
  xAVec = NULL;
  xBVec = NULL;
  
  updatePoints();
}

/*
  Destructor for spherical constraint
*/
TACSSphericalConstraint::~TACSSphericalConstraint(){
  bodyA->decref();
  if(bodyB){ bodyB->decref(); }
  point->decref();
  xAVec->decref();
  if(xBVec){ xBVec->decref(); }
}

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSSphericalConstraint::numNodes(){
  if(bodyA && bodyB){
    return 3;
  } 
  else {
    return 2;
  }
}

const char *TACSSphericalConstraint::elem_name = "TACSSphericalConstraint";

/*
  Update the local data. This takes the initial reference points for
  the bodies A and B, given as rA and rB, and the initial point for
  the connection between the A and B bodies, "point pt", and computes
  the vectors xA and xB in the global frame.
*/
void TACSSphericalConstraint::updatePoints(){
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Fetch the positions of body in global frame
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Determine the position of the joint from bodyA in the global frame
  // xAVec = point - rAVec
  TacsScalar xA[3];
  for ( int i = 0; i < 3; i++ ){
    xA[i] = pt[i] - rA[i];
  }
  if (xAVec) {
    xAVec->decref();
  }
  xAVec = new TACSGibbsVector(xA);
  xAVec->incref();

  if (bodyB){
    // Fetch the positions of body in global frame
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    const TacsScalar *rB;
    rBVec->getVector(&rB);

    // Determine the position of the joint from bodyB in the global frame
    // xBVec = point - rBVec
    TacsScalar xB[3];
    for ( int i = 0; i < 3; i++ ){
      xB[i] = pt[i] - rB[i];
    }
    if (xBVec) {
      xBVec->decref();
    }
    xBVec = new TACSGibbsVector(xB);
    xBVec->incref();
  }
}

/*
  Compute the kinetic and potential energy within the element
*/
void TACSSphericalConstraint::computeEnergies( double time,
                                               TacsScalar *_Te, 
                                               TacsScalar *_Pe,
                                               const TacsScalar Xpts[],
                                               const TacsScalar vars[],
                                               const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

/*
  Compute the residual of the governing equations
*/
void TACSSphericalConstraint::addResidual( double time, TacsScalar res[],
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[] ){
  // Retrieve the joint location from  global origin
  const TacsScalar *pt;
  point->getVector(&pt);

  // Get the initial position for bodyA  from global origin
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Set pointers to the residual of each body
  TacsScalar *resA = &res[0];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];
  
  // The residual for the constraint equations
  TacsScalar *resC = NULL;

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not.
  if (bodyB){
    resC = &res[16];
    lam = &vars[16];
  }
  else {
    resC = &res[8];
    lam = &vars[8];
  }

  // Compute the rotation matrices for each body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Retrieve the pointers to xAVec and xBVec
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Add the terms for body A
  vecAxpy(1.0, lam, &resA[0]);
  addEMatTransProduct(1.0, xA, lam, etaA, epsA, 
                      &resA[3], &resA[4]);

  // Evaluate the constraint 
  // resC = rA + uA + CA^{T}*xA - pt = 0 or 
  // resC = rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB = 0
  matMultTrans(CA, xA, resC);
  vecAxpy(1.0, uA, resC);
  vecAxpy(1.0, rA, resC);
 
  if (bodyB){
    // Set the residual for bodyB
    TacsScalar *resB = &res[8];

    // Set the variables for body B
    const TacsScalar *uB = &vars[8];
    const TacsScalar etaB = vars[11];
    const TacsScalar *epsB = &vars[12];

    // Compute the rotation matrix for bodyB
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);
    
    const TacsScalar *xB;
    xBVec->getVector(&xB);
    
    // Add the terms for body B
    vecAxpy(-1.0, lam, &resB[0]);
    addEMatTransProduct(-1.0, xB, lam, etaB, epsB, 
                        &resB[3], &resB[4]);

    // Compute t = CB^{T}*xB + uB
    TacsScalar t[3];
    matMultTrans(CB, xB, t);
    vecAxpy(1.0, uB, t);

    // Complete the evaluation of the constraint
    vecAxpy(-1.0, t, resC);

    // Get the initial position for bodyB
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    const TacsScalar *rB;
    rBVec->getVector(&rB);
    vecAxpy(-1.0, rB, resC);

  } else {
    // subtract the joint location if bodyB is not present
    vecAxpy(-1.0, pt, resC);      
  }
  

  // Add the dummy constraints for the remaining Lagrange multiplier
  // variables
  for ( int i = 3; i < 8; i++ ){
    resC[i] += lam[i];
  }
}

/*
  Compute the Jacobian of the residuals of the governing equations
*/
void TACSSphericalConstraint::addJacobian( double time, TacsScalar J[],
                                           double alpha, 
                                           double beta, 
                                           double gamma,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[] ){
  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the Lagrange multipliers for the constraint
  const TacsScalar *lam = NULL;

  // Get the number of variables
  const int nvars = numVariables();
  
  // Set the offset to the Lagrange multipliers
  int offset = 0;
  if (bodyB){
    offset = 16;
    lam = &vars[16];
  }
  else {
    offset = 8;
    lam = &vars[8];
  }

  // Add the identity matricies to the Jacobian
  addBlockIdent(alpha, &J[offset], nvars);
  addBlockIdent(alpha, &J[offset*nvars], nvars);

  // Retrieve the pointers to xAVec
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Add the second derivative terms
  addBlockDMatTransDeriv(alpha, lam, xA, &J[3*(nvars+1)], nvars);

  // Add the term from the derivative w.r.t. lambda
  addBlockEMatTrans(alpha, etaA, epsA, xA, &J[3*nvars + offset], nvars);
  
  // Add the term from the derivative of the constraint
  addBlockEMat(alpha, etaA, epsA, xA, &J[offset*nvars + 3], nvars);

  // Add the terms required for body B if it is defined
  if (bodyB){
    addBlockIdent(-alpha, &J[8*nvars + offset], nvars);
    addBlockIdent(-alpha, &J[offset*nvars + 8], nvars);

    // Set the variables for body B
    const TacsScalar etaB = vars[11];
    const TacsScalar *epsB = &vars[12];

    // Retrieve the pointer to xBVec
    const TacsScalar *xB;
    xBVec->getVector(&xB);

    // Add the second derivative terms
    addBlockDMatTransDeriv(-alpha, lam, xB, &J[11*(nvars+1)], nvars);
    
    // Add the terms from the derivatives w.r.t. lambdas
    addBlockEMatTrans(-alpha, etaB, epsB, xB, &J[11*nvars + offset], nvars);
    
    // Add the terms from the derivatives of the constraint
    addBlockEMat(-alpha, etaB, epsB, xB, &J[offset*nvars + 11], nvars);
  }

  // Add the Jacobian entries for the dummy constraints
  for ( int i = offset+3; i < nvars; i++ ){
    J[(nvars+1)*i] += alpha;
  }
}

/*
  Set the design variable values
*/
void TACSSphericalConstraint::setDesignVars( const TacsScalar dvs[], 
                                             int numDVs ){
  point->setDesignVars(dvs, numDVs);
  updatePoints();
}

/*
  Get the design variable values associated with the joint location
*/
void TACSSphericalConstraint::getDesignVars( TacsScalar dvs[], int numDVs ){
  point->getDesignVars(dvs, numDVs);
}


/*
  Constructor for revolute contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA:  pointer to bodyA
  bodyB:  pointer to bodyB
  point:  the position of the joint from the global reference point
  rev:    the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint( TACSRigidBody *_bodyA, 
                                                TACSRigidBody *_bodyB, 
                                                TACSGibbsVector *_point, 
                                                TACSGibbsVector *_eAVec ){
  // Copy over the input arguments
  fixed_ref_point = 0;
  bodyA = _bodyA; bodyA->incref();
  bodyB = _bodyB; bodyB->incref();
  point = _point; point->incref();
  eAVec = _eAVec; eAVec->incref();

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Constructor for revolute contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA:  pointer to bodyA
  point:  the position of the joint from the global reference point
  rev:    the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint( TACSRigidBody *_bodyA, 
                                                TACSGibbsVector *_point, 
                                                TACSGibbsVector *_eAVec ){
  // Copy over the input arguments
  fixed_ref_point = 1;
  bodyA = _bodyA; bodyA->incref();
  bodyB = NULL;
  point = _point; point->incref();
  eAVec = _eAVec; eAVec->incref();

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Create a revolute constraint joining two nodes together directly - 
  this can be used to connect two flexible bodies

  input:
  fixed_ref_point: Is the reference point fixed?
  point:   the position of the joint from the global reference point
  rev:     the revolute direction in global frame
*/
TACSRevoluteConstraint::TACSRevoluteConstraint( int _fixed_ref_point,
                                                TACSGibbsVector *_point,
                                                TACSGibbsVector *_eAVec ){
  // Copy over the input arguments
  fixed_ref_point = _fixed_ref_point;
  bodyA = bodyB = NULL;
  point = _point; point->incref();
  eAVec = _eAVec; eAVec->incref();

  // Set class variables to NULL
  eB1Vec = eB2Vec = eVec = NULL;

  int init_vector = 1;
  updatePoints(init_vector);
}

/*
  Destuctor for the revolute constraint
*/
TACSRevoluteConstraint::~TACSRevoluteConstraint(){
  if (bodyA){ bodyA->decref(); }
  if (bodyB){ bodyB->decref(); }
  point->decref();
  eAVec->decref();
  eVec->decref();
  eB1Vec->decref();
  eB2Vec->decref();
}

const char *TACSRevoluteConstraint::elem_name = "TACSRevoluteConstraint";

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSRevoluteConstraint::numNodes(){
  if (fixed_ref_point){
    return 2;
  } 
  else {
    return 3;
  }
}

/*
  Read the data from the given initial point vectors/locations and
  re-compute the internal data that is requied to evaluate the
  kinematic constraints.
*/
void TACSRevoluteConstraint::updatePoints( int init_vector ){
  // Find the minimum absolute component of eAVec along any coordinate
  // direction. Set the vector components of e along this direction
  // to maximize orthogonality among the coordinate directions. For
  // the purpose of optimization, this direction is fixed at
  // initialization.
  // Retrieve the revolute direction in global frame
  const TacsScalar *rev;
  eAVec->getVector(&rev);

  TacsScalar e[3];
  if (init_vector){
    e[0] = e[1] = e[2] = 0.0;
    if ((fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[1]))) && 
        (fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[2])))){
      e[0] = 1.0;
    }
    else if ((fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[0]))) && 
             (fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[2])))){
      e[1] = 1.0;
    }
    else {
      e[2] = 1.0;
    }
    eVec = new TACSGibbsVector(e);
    eVec->incref();
  }
  else {
    const TacsScalar *etmp;
    eVec->getVector(&etmp);
    memcpy(e, etmp, 3*sizeof(TacsScalar));
  }

  // Compute/recompute the eB1 and eB2 directions based on e
  TacsScalar eB1[3], eB2[3];
  crossProduct(1.0, rev, e, eB2);
  if (eB2Vec){
    eB2Vec->decref();
  }
  eB2Vec = new TACSGibbsVector(eB2);
  eB2Vec->incref();
  
  crossProduct(1.0, eB2, rev, eB1);
  if (eB1Vec){
    eB1Vec->decref();
  }
  eB1Vec = new TACSGibbsVector(eB1);
  eB1Vec->incref();
}

/*
  Compute the kinetic and potential energy within the element
*/
void TACSRevoluteConstraint::computeEnergies( double time,
                                              TacsScalar *_Te, 
                                              TacsScalar *_Pe,
                                              const TacsScalar Xpts[],
                                              const TacsScalar vars[],
                                              const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

/*
  Compute the residual of the governing equations
*/
void TACSRevoluteConstraint::addResidual( double time, TacsScalar res[],
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[] ){
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Get the initial position vectors - depending on whether we're
  // taking the initial position from the node locations
  const TacsScalar *rA, *rB;
  if (bodyA){
    // Get the initial position for bodyA
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  }
  else {
    rA = &Xpts[0];
  }
  if (bodyB){
    // Get the initial position for bodyA
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  }
  else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set pointers to the residual of each body
  TacsScalar *resA = &res[0];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar *uB = &vars[8];
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // The residual for the constraint equations
  TacsScalar *resC = NULL;

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not.
  if (fixed_ref_point){
    resC = &res[8];
    lam = &vars[8];
  }
  else {
    resC = &res[16];
    lam = &vars[16];
  }

  // Retrieve the pointers to eA, eB1, eB2
  const TacsScalar *eA, *eB1, *eB2;
  eAVec->getVector(&eA);
  eB1Vec->getVector(&eB1);
  eB2Vec->getVector(&eB2);

  // Compute the rotation matrices for each body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Add the terms for body A
  vecAxpy(1.0, lam, &resA[0]);
  addEMatTransProduct(1.0, xA, lam, etaA, epsA, 
                      &resA[3], &resA[4]);

  // Evaluate the constraint 
  // resC = rA + uA + CA^{T}*xA - pt = 0 or 
  // resC = rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB = 0
  matMultTransAdd(CA, xA, resC);
  vecAxpy(1.0, uA, resC);
  vecAxpy(1.0, rA, resC);

  if (!fixed_ref_point){
    // Set the residual for bodyB
    TacsScalar *resB = &res[8];

    // Compute the rotation matrix for bodyB
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);
        
    // Add the terms for body B
    vecAxpy(-1.0, lam, &resB[0]);
    addEMatTransProduct(-1.0, xB, lam, etaB, epsB, 
                        &resB[3], &resB[4]);

    // Compute t = CB^{T}*xB + uB
    TacsScalar t[3];
    matMultTrans(CB, xB, t);
    vecAxpy(1.0, uB, t);

    // Complete the evaluation of the constraint
    vecAxpy(-1.0, t, resC);

    // Get the initial position for bodyB
    vecAxpy(-1.0, rB, resC);

    // Add the revolute direction constraint
    TacsScalar tA[3], tB1[3], tB2[3];
    matMultTrans(CA, eA, tA);
    matMultTrans(CB, eB1, tB1);
    matMultTrans(CB, eB2, tB2);

    // Compute the contributions to the first revolute constraint
    resC[3] += vecDot(tA, tB1);
    
    // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
    addEMatTransProduct(lam[3], eA, tB1, etaA, epsA, 
                        &resA[3], &resA[4]);
    // Add the derivative d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA
    addEMatTransProduct(lam[3], eB1, tA, etaB, epsB,
                        &resB[3], &resB[4]);
    
    // Compute the contributions to the second revolute constraint
    resC[4] += vecDot(tA, tB2);
    
    // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
    addEMatTransProduct(lam[4], eA, tB2, etaA, epsA, 
                        &resA[3], &resA[4]);
    // Add the derivative d(CB^{T}*eB2)/dqB)*tA = E(eB2)*tA
    addEMatTransProduct(lam[4], eB2, tA, etaB, epsB,
                        &resB[3], &resB[4]);
  }
  else {
    // Subtract pt vec if bodyB is not present
    vecAxpy(-1.0, pt, resC);
  
    // Add the revolute direction constraint
    TacsScalar tA[3];
    matMultTrans(CA, eA, tA);

    // Compute the contributions to the first revolute constraint
    resC[3] += vecDot(tA, eB1);
    
    // Add the derivative (d(CA^{T}*eA)/dqA)*eB1 = E(eA)*tB
    addEMatTransProduct(lam[3], eA, eB1, etaA, epsA, 
                        &resA[3], &resA[4]);
    
    // Compute the contributions to the second revolute constraint
    resC[4] += vecDot(tA, eB2);

    // Add the derivative (d(CA^{T}*eA)/dqA)*eB2 = E(eA)*eB2
    addEMatTransProduct(lam[4], eA, eB2, etaA, epsA, 
                        &resA[3], &resA[4]);
  }

  // Add the dummy constraints for the remaining Lagrange multiplier
  // variables
  for ( int i = 5; i < 8; i++ ){
    resC[i] += lam[i];
  }
}

/*
  Compute the Jacobian of the residuals of the governing equations
*/
void TACSRevoluteConstraint::addJacobian( double time, TacsScalar J[],
                                          double alpha, 
                                          double beta, 
                                          double gamma,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[] ){
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Get the number of variables
  const int nvars = numVariables();

  // Get the initial position vectors - depending on whether we're
  // taking the initial position from the node locations
  const TacsScalar *rA, *rB;
  if (bodyA){
    // Get the initial position for bodyA
    TACSGibbsVector *rAVec = bodyA->getInitPosition();
    rAVec->getVector(&rA);
  }
  else {
    rA = &Xpts[0];
  }
  if (bodyB){
    // Get the initial position for bodyA
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    rBVec->getVector(&rB);
  }
  else {
    rB = &Xpts[3];
  }

  // Set the positions relative to the initial location
  TacsScalar xA[3];
  xA[0] = pt[0] - rA[0];
  xA[1] = pt[1] - rA[1];
  xA[2] = pt[2] - rA[2];

  TacsScalar xB[3];
  xB[0] = pt[0] - rB[0];
  xB[1] = pt[1] - rB[1];
  xB[2] = pt[2] - rB[2];

  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for body B
  const TacsScalar etaB = vars[11];
  const TacsScalar *epsB = &vars[12];

  // The Lagrange multipliers
  const TacsScalar *lam = NULL;

  // Set the pointers depending on whether both body A and body B
  // exist or not. Also, set the offset to the Lagrange multipliers
  int offset = 0;
  if (fixed_ref_point){
    offset = 8;
    lam = &vars[8];
  }
  else {
    offset = 16;
    lam = &vars[16];
  }

  // Retrieve the pointers to eA, eB1, eB2
  const TacsScalar *eA, *eB1, *eB2;
  eAVec->getVector(&eA);
  eB1Vec->getVector(&eB1);
  eB2Vec->getVector(&eB2);
  
  // Compute the rotation matrix
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Add the identity matricies to the Jacobian
  addBlockIdent(alpha, &J[offset], nvars);
  addBlockIdent(alpha, &J[offset*nvars], nvars);

  // Add the second derivative terms
  addBlockDMatTransDeriv(alpha, lam, xA, &J[3*(nvars+1)], nvars);

  // Add the term from the derivative w.r.t. lambda
  addBlockEMatTrans(alpha, etaA, epsA, xA, &J[3*nvars + offset], nvars);
  
  // Add the term from the derivative of the constraint
  addBlockEMat(alpha, etaA, epsA, xA, &J[offset*nvars + 3], nvars);

  // Add the terms required for body B if it is defined
  if (!fixed_ref_point){
    // Compute the rotation matrix
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);

    // Add the block identities
    addBlockIdent(-alpha, &J[8*nvars + offset], nvars);
    addBlockIdent(-alpha, &J[offset*nvars + 8], nvars);

    // Add the second derivative terms
    addBlockDMatTransDeriv(-alpha, lam, xB, &J[11*(nvars+1)], nvars);
    
    // Add the terms from the derivatives w.r.t. lambdas
    addBlockEMatTrans(-alpha, etaB, epsB, xB, &J[11*nvars + offset], nvars);
    
    // Add the terms from the derivatives of the constraint
    addBlockEMat(-alpha, etaB, epsB, xB, &J[offset*nvars + 11], nvars);

    // Add the revolute direction constraint
    TacsScalar tA[3], tB1[3], tB2[3];
    matMultTrans(CA, eA, tA);
    matMultTrans(CB, eB1, tB1);
    matMultTrans(CB, eB2, tB2);
    
    TacsScalar gA[4], gB[4];
    
    // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, tB1, etaA, epsA, 
                        &gA[0], &gA[1]);
    // Add the derivative d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA
    memset(gB, 0, sizeof(gB));
    addEMatTransProduct(alpha, eB1, tA, etaB, epsB,
                        &gB[0], &gB[1]);
    
    // Add the constraint terms to the matrix
    for ( int i = 0; i < 4; i++ ){
      J[nvars*(i+3) + 19] += gA[i];
      J[nvars*(i+11) + 19] += gB[i];
      J[19*nvars + i+3] += gA[i];
      J[19*nvars + i+11] += gB[i];
    }
    
    // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, tB2, etaA, epsA, 
                        &gA[0], &gA[1]);
    // Add the derivative d(CB^{T}*eB2)/dqB)*tA = E(eB2)*tA
    memset(gB, 0, sizeof(gB));
    addEMatTransProduct(alpha, eB2, tA, etaB, epsB,
                        &gB[0], &gB[1]);
    
    // Add the constraint terms to the matrix
    for ( int i = 0; i < 4; i++ ){
      J[nvars*(i+3) + 20] += gA[i];
      J[nvars*(i+11) + 20] += gB[i];
      J[20*nvars + i+3] += gA[i];
      J[20*nvars + i+11] += gB[i];
    }
    
    // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
    addBlockDMatTransDeriv(alpha*lam[3], tB1, eA, &J[3*(nvars+1)], nvars);
    addBlockDMatTransDeriv(alpha*lam[3], tA, eB1, &J[11*(nvars+1)], nvars);
    
    // Add the contributions from the off-diagonal blocks
    TacsScalar EA[12], EB[12];
    computeEMat(etaA, epsA, eA, EA);
    computeEMat(etaB, epsB, eB1, EB);
    
    addBlock3x4Product(alpha*lam[3], EA, EB, &J[3*nvars+11], nvars);
    addBlock3x4Product(alpha*lam[3], EB, EA, &J[11*nvars+3], nvars);
    
    // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
    addBlockDMatTransDeriv(alpha*lam[4], tB2, eA, &J[3*(nvars+1)], nvars);
    addBlockDMatTransDeriv(alpha*lam[4], tA, eB2, &J[11*(nvars+1)], nvars);
    
    // Add the contributions from the off-diagonal blocks for the second
    // constraint
    computeEMat(etaB, epsB, eB2, EB);
    
    addBlock3x4Product(alpha*lam[4], EA, EB, &J[3*nvars+11], nvars);
    addBlock3x4Product(alpha*lam[4], EB, EA, &J[11*nvars+3], nvars);
  }
  else {
    // Add the revolute direction constraint
    TacsScalar tA[3];
    matMultTrans(CA, eA, tA);

    // Add the derivative (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB
    TacsScalar gA[4];
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, eB1, etaA, epsA, 
                        &gA[0], &gA[1]);
    
    // Add the constraint terms to the matrix
    for ( int i = 0; i < 4; i++ ){
      J[nvars*(i+3) + offset+3] += gA[i];
      J[(offset+3)*nvars + i+3] += gA[i];
    }

    // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
    memset(gA, 0, sizeof(gA));
    addEMatTransProduct(alpha, eA, eB2, etaA, epsA, 
                        &gA[0], &gA[1]);

    for ( int i = 0; i < 4; i++ ){
      J[nvars*(i+3) + offset+4] += gA[i];
      J[(offset+4)*nvars + i+3] += gA[i];
    }

    // Add the diagonal contributions to the constraint tA^{T}*tB1 = 0
    addBlockDMatTransDeriv(alpha*lam[3], eB1, eA, &J[3*(nvars+1)], nvars);
        
    // Add the diagonal contributions to the constraint tA^{T}*tB2 = 0
    addBlockDMatTransDeriv(alpha*lam[4], eB2, eA, &J[3*(nvars+1)], nvars);
  }

  // Add the Jacobian entries for the dummy constraints
  for ( int i = offset+5; i < nvars; i++ ){
    J[(nvars+1)*i] += alpha;
  }
}

/*
  Set the design variable values
*/
void TACSRevoluteConstraint::setDesignVars( const TacsScalar dvs[], 
                                            int numDVs ){
  point->setDesignVars(dvs, numDVs);
  updatePoints();
}

/*
  Get the design variable values associated with the joint location
*/
void TACSRevoluteConstraint::getDesignVars( TacsScalar dvs[], int numDVs ){
  point->getDesignVars(dvs, numDVs);
}

/*
  Two-point rigid link constraint
*/
TACSRigidLink::TACSRigidLink( TACSRigidBody *_bodyA ){
  bodyA = _bodyA;
  bodyA->incref();
}

TACSRigidLink::~TACSRigidLink(){
  bodyA->decref();
}

const char *TACSRigidLink::elem_name = "TACSRigidLink";

/*
  Return the number of displacements
*/
int TACSRigidLink::numDisplacements(){ 
  return 8; 
}

/*
  Return the number of nodes
*/
int TACSRigidLink::numNodes(){ 
  return 3; 
}

/*
  Return the element name
*/
const char* TACSRigidLink::elementName(){ 
  return elem_name;
}

/*
  Compute the kinetic and potential energy within the element
*/
void TACSRigidLink::computeEnergies( double time,
                                     TacsScalar *_Te, 
                                     TacsScalar *_Pe,
                                     const TacsScalar Xpts[],
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

/*
  Compute the residual of the governing equations
*/
void TACSRigidLink::addResidual( double time, TacsScalar res[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){
  // Set pointers to the residual
  TacsScalar *resA = &res[0];
  TacsScalar *resB = &res[8];
  TacsScalar *resC = &res[16];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the variables for point B
  const TacsScalar *uB = &vars[8];
  const TacsScalar *epsB = &vars[12];

  // Set the pointer to the multipliers
  const TacsScalar *lam = &vars[16];

  // Retrieve the initial position of body A
  TACSGibbsVector *xAVec = bodyA->getInitPosition();
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Read from the node locations, the initial position of body B
  const TacsScalar *xB = &Xpts[3];

  // Compute the rotation matrix for body A
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Compute the distance between body A and the point B in the
  // initial configuration
  TacsScalar t[3];
  t[0] = xA[0] - xB[0];
  t[1] = xA[1] - xB[1];
  t[2] = xA[2] - xB[2];

  // Add the residual 
  // resC = uB - uA + (xB - xA) + CA^{T}*(xA - xB)
  // resC = uB - uA + t + CA^{T}*(xA - xB)
  resC[0] += uB[0] - uA[0];
  resC[1] += uB[1] - uA[1];
  resC[2] += uB[2] - uA[2];
  
  vecAxpy(-1.0, t, resC);
  matMultTransAdd(CA, t, resC);

  // Add the residual for the quaternions
  resC[3] += lam[3];
  resC[4] += epsB[0] - epsA[0];
  resC[5] += epsB[1] - epsA[1];
  resC[6] += epsB[2] - epsA[2];
  
  // Add the dummy constraint for the remaining multiplier 
  resC[7] += lam[7];

  // Add the terms from the first constraint
  vecAxpy(-1.0, &lam[0], &resA[0]);
  addEMatTransProduct(-1.0, t, &lam[0], etaA, epsA,
                      &resA[3], &resA[4]);
  vecAxpy(1.0, &lam[0], &resB[0]);
  
  // Add the terms from the second constraint
  vecAxpy(-1.0, &lam[4], &resA[4]);
  vecAxpy(1.0, &lam[4], &resB[4]);
}

/*
  Compute the Jacobian of the governing equations
*/
void TACSRigidLink::addJacobian( double time, TacsScalar J[],
                                 double alpha, double beta, double gamma,
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){
  // Set the variables for body A
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];

  // Set the Lagrange multiplier variables
  const TacsScalar *lam = &vars[16];

  // The number of variables in the Jacobian matrix
  const int nvars = 3*8;

  // Retrieve the initial position of body A
  TACSGibbsVector *xAVec = bodyA->getInitPosition();
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Read from the node locations, the initial position of body B
  const TacsScalar *xB = &Xpts[3];

  // Compute the rotation matrix for body A
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // Compute the distance between body A and the point B in the
  // initial configuration
  TacsScalar t[3];
  t[0] = xA[0] - xB[0];
  t[1] = xA[1] - xB[1];
  t[2] = xA[2] - xB[2];

  // Derivatives of the position constraint
  addBlockIdent(-alpha, &J[16*nvars], nvars);
  addBlockIdent(alpha, &J[16*nvars+8], nvars);
  addBlockEMat(alpha, etaA, epsA, t, &J[16*nvars + 3], nvars);

  // Derivatives of the quaternion constraint
  addBlockIdent(alpha, &J[20*nvars + 12], nvars);
  addBlockIdent(-alpha, &J[20*nvars + 4], nvars);

  // Add the Jacobian contribution from the dummy constraint
  J[nvars*nvars-1] += alpha; 

  // Add the contributions from the derivative of resA
  addBlockIdent(-alpha, &J[16], nvars);
  addBlockIdent(alpha, &J[8*nvars + 16], nvars);
  addBlockDMatTransDeriv(-alpha, lam, t, &J[3*nvars + 3], nvars);
  addBlockEMatTrans(-alpha, etaA, epsA, t, &J[3*nvars + 16], nvars);

  // Add the derivatives of the quaternion constraint w.r.t. lam[3]
  J[19*nvars + 19] += alpha;
  
  // Add the remaining quaternion constraint derivatives w.r.t. lam[4:]
  addBlockIdent(-alpha, &J[4*nvars + 20], nvars);
  addBlockIdent(alpha, &J[12*nvars + 20], nvars);
}

TACSRevoluteDriver::TACSRevoluteDriver( TACSGibbsVector *orig, 
                                        TACSGibbsVector *rev,
                                        TacsScalar _omega ){
  origVec = orig;
  revVec = rev;
  origVec->incref();
  revVec->incref();
  omega = _omega;
}

TACSRevoluteDriver::~TACSRevoluteDriver(){
  origVec->decref();
  revVec->decref();
}

int TACSRevoluteDriver::numDisplacements(){ 
  return 8; 
}
  
int TACSRevoluteDriver::numNodes(){ 
  return 2; 
}

const char* TACSRevoluteDriver::elementName(){ 
  return "TACSRevoluteDriver"; 
}

void TACSRevoluteDriver::computeEnergies( double time,
                                          TacsScalar *_Te, 
                                          TacsScalar *_Pe,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

void TACSRevoluteDriver::addResidual( double time, TacsScalar res[],
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[],
                                      const TacsScalar ddvars[] ){
  // The displacements at the final time are given
  TacsScalar theta = omega*time;
  TacsScalar s = sin(theta);
  TacsScalar c = cos(theta);

  // Extract the components of the vector
  const TacsScalar *rev, *orig;
  revVec->getVector(&rev);
  origVec->getVector(&orig);

  // Compute the initial position between the 
  TacsScalar d1[3];
  d1[0] = Xpts[0] - orig[0];
  d1[1] = Xpts[1] - orig[1];
  d1[2] = Xpts[2] - orig[2];

  // Compute the component of the vector d1 perpendicular to the
  // revolute direction
  TacsScalar p[3];
  TacsScalar a = vecDot(rev, d1);
  p[0] = d1[0] - a*rev[0];
  p[1] = d1[1] - a*rev[1];
  p[2] = d1[2] - a*rev[2];
    
  // Find the cross product e = rev x d1
  TacsScalar e[3];
  crossProduct(1.0, rev, d1, e);

  // Combine the perpendicular components of the vector
  TacsScalar d2[3];
  d2[0] = a*rev[0] + c*p[0] + s*e[0];
  d2[1] = a*rev[1] + c*p[1] + s*e[1];
  d2[2] = a*rev[2] + c*p[2] + s*e[2];

  // Compute the displacement from the initial point
  TacsScalar u[3];
  u[0] = d2[0] - d1[0];
  u[1] = d2[1] - d1[1];
  u[2] = d2[2] - d1[2];
    
  // Set the pointers to the Lagrange multipliers
  const TacsScalar *lam = &vars[8];

  // Add the multipliers (forces) to the equations of motion
  res[0] += lam[0];
  res[1] += lam[1];
  res[2] += lam[2];

  // Add the kinematic constraint equations
  res[8] += vars[0] - u[0];
  res[9] += vars[1] - u[1];
  res[10] += vars[2] - u[2];

  // Add dummy constraints for the remaining multipliers
  res[11] += lam[3];
  res[12] += lam[4];
  res[13] += lam[5];
  res[14] += lam[6];
  res[15] += lam[7];
}
  
void TACSRevoluteDriver::addJacobian( double time, TacsScalar J[],
                                      double alpha, double beta, double gamma,
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[],
                                      const TacsScalar ddvars[] ){
  // The number of variables in the Jacobian matrix
  const int nvars = 2*8;

  // Add the block from the multipliers
  addBlockIdent(alpha, &J[8], nvars);

  // Add the block from the constraint equations
  addBlockIdent(alpha, &J[8*nvars], nvars);

  // Add the diagonal block from the dummy equations
  for ( int i = 11; i < nvars; i++ ){
    J[(nvars+1)*i] += alpha;
  }
}

/*
  Constructor for cylindrical contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA : pointer to bodyA
  bodyB : pointer to bodyB
  point : the position of the joint from the global reference point
  rev   : the cylindrical direction in global frame
*/
TACSCylindricalConstraint::TACSCylindricalConstraint( TACSRigidBody *_bodyA, 
                                                      TACSRigidBody *_bodyB, 
                                                      TACSGibbsVector *_point, 
                                                      TACSGibbsVector *_eAVec ){
  // Copy over the input arguments
  bodyA = _bodyA; bodyA->incref();
  bodyB = _bodyB; bodyB->incref();
  point = _point; point->incref();
  eAVec = _eAVec; eAVec->incref();

  // Set class variables to NULL
  xAVec  = xBVec = NULL;
  eB1Vec = eB2Vec = eVec = NULL;

  int init_e = 1;
  updatePoints(init_e);
}

/*
  Constructor for cylindrical contraint taking Gibbs vectors as
  inputs. A refers to bodyA and B refers to bodyB.

  input:
  bodyA : pointer to bodyA
  point : the position of the joint from the global reference point
  rev   : the cylindrical direction in global frame
*/
TACSCylindricalConstraint::TACSCylindricalConstraint( TACSRigidBody *_bodyA, 
                                                      TACSGibbsVector *_point, 
                                                      TACSGibbsVector *_eAVec ){
  // Copy over the input arguments
  bodyA = _bodyA; bodyA->incref();
  bodyB = NULL;
  point = _point; point->incref();
  eAVec = _eAVec; eAVec->incref();

  // Set class variables to NULL
  xAVec  = xBVec = NULL;
  eB1Vec = eB2Vec = eVec = NULL;

  int init_e = 1;
  updatePoints(init_e);
}

/*
  Destuctor for the cylindrical constraint
*/
TACSCylindricalConstraint::~TACSCylindricalConstraint(){
  bodyA->decref();
  if(bodyB){ bodyB->decref(); }
  point->decref();

  eAVec->decref();
  eVec->decref();

  xAVec->decref();
  if(xBVec){ xBVec->decref(); }

  eB1Vec->decref();
  eB2Vec->decref();
}

const char *TACSCylindricalConstraint::elem_name = "TACSCylindricalConstraint";

/*
  Returns the number of nodes based on the constraint nature
*/
int TACSCylindricalConstraint::numNodes(){
  if(bodyA && bodyB){
    return 3;
  } 
  else {
    return 2;
  }
}

/*
  Read the data from the given initial point vectors/locations and
  re-compute the internal data that is requied to evaluate the
  kinematic constraints.

  This code computes the required vectors in the body-fixed coordinate
  frames. In particular, the code computes:
  
  xA = pt - rA
  xB = pt - rB

  where pt is the attachment point, and rA and rB are the initial
  points of bodies A and B in the global (inertial) reference frame.
*/
void TACSCylindricalConstraint::updatePoints( int init_e ){
  // Retrieve the coordinates of the joint point in the global frame
  const TacsScalar *pt;
  point->getVector(&pt);

  // Fetch the positions of body in global frame
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Determine the position of the joint from bodyA in the global frame
  // xAVec = point - rAVec
  TacsScalar xA[3];
  for ( int i = 0; i < 3; i++ ){
    xA[i] = pt[i] - rA[i];
  }
  if (xAVec) {
    xAVec->decref();
  }
  xAVec = new TACSGibbsVector(xA);
  xAVec->incref();

  if (bodyB){
    // Fetch the positions of body in global frame
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    const TacsScalar *rB;
    rBVec->getVector(&rB);

    // Determine the position of the joint from bodyB in the global frame
    // xBVec = point - rBVec
    TacsScalar xB[3];
    for ( int i = 0; i < 3; i++ ){
      xB[i] = pt[i] - rB[i];
    }
    if (xBVec) {
      xBVec->decref();
    }
    xBVec = new TACSGibbsVector(xB);
    xBVec->incref();
  }

  // Find the minimum absolute component of eAVec along any coordinate
  // direction. Set the vector components of e along this direction
  // to maximize orthogonality among the coordinate directions. For
  // the purpose of optimization, this direction is fixed at
  // initialization.
  // Retrieve the cylindrical direction in global frame
  const TacsScalar *rev;
  eAVec->getVector(&rev);

  TacsScalar e[3];
  if (init_e){
    e[0] = e[1] = e[2] = 0.0;
    if ((fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[1]))) && 
        (fabs(TacsRealPart(rev[0])) <= fabs(TacsRealPart(rev[2])))){
      e[0] = 1.0;
    }
    else if ((fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[0]))) && 
             (fabs(TacsRealPart(rev[1])) <= fabs(TacsRealPart(rev[2])))){
      e[1] = 1.0;
    }
    else {
      e[2] = 1.0;
    }
    eVec = new TACSGibbsVector(e);
    eVec->incref();
  }
  else {
    const TacsScalar *etmp;
    eVec->getVector(&etmp);
    memcpy(e, etmp, 3*sizeof(TacsScalar));
  }

  // Compute/recompute the eB1 and eB2 directions based on e
  TacsScalar eB1[3], eB2[3];
  crossProduct(1.0, rev, e, eB2);
  if (eB2Vec){
    eB2Vec->decref();
  }
  eB2Vec = new TACSGibbsVector(eB2);
  eB2Vec->incref();
  
  crossProduct(1.0, eB2, rev, eB1);
  if (eB1Vec){
    eB1Vec->decref();
  }
  eB1Vec = new TACSGibbsVector(eB1);
  eB1Vec->incref();
}

/*
  Compute the kinetic and potential energy within the element
*/
void TACSCylindricalConstraint::computeEnergies( double time,
                                                 TacsScalar *_Te, 
                                                 TacsScalar *_Pe,
                                                 const TacsScalar Xpts[],
                                                 const TacsScalar vars[],
                                                 const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

/*
  Compute the residual of the governing equations
*/
void TACSCylindricalConstraint::addResidual( double time, TacsScalar res[],
                                             const TacsScalar Xpts[],
                                             const TacsScalar vars[],
                                             const TacsScalar dvars[],
                                             const TacsScalar ddvars[] ){
  // Retrieve the joint location from  global origin
  const TacsScalar *pt;
  point->getVector(&pt);

  // Set the bodyA pointers
  TacsScalar *resA       = &res[0];
  const TacsScalar *uA   = &vars[0];
  const TacsScalar etaA  = vars[3];
  const TacsScalar *epsA = &vars[4];
  
  // Set the constraint pointers
  TacsScalar *resC = NULL;
  const TacsScalar *lam = NULL;
  if (bodyB){
    resC = &res[16];
    lam = &vars[16];
  }
  else {
    resC = &res[8];
    lam = &vars[8];
  }

  // Retrieve the pointers to eA, eB1, eB2
  const TacsScalar *eA, *eB1, *eB2;
  eAVec->getVector(&eA);
  eB1Vec->getVector(&eB1);
  eB2Vec->getVector(&eB2);

  // Compute the rotation matrices for bodyA
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  TacsScalar tA[3];
  matMultTrans(CA, eA, tA);

  // Get the initial position for bodyA
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Retrieve the pointers to xAVec and xBVec
  const TacsScalar *xA;
  xAVec->getVector(&xA);

  // Evaluate the constraint 
  // resC = rA + uA + CA^{T}*xA - pt = 0 or 
  // resC = rA + uA + CA^{T}*xA - rB - uB - CB^{T}*xB = 0

  TacsScalar dA[3];
  matMultTransAdd(CA, xA, dA);
  vecAxpy(1.0, uA, dA);
  vecAxpy(1.0, rA, dA);
 
  if (bodyB){
    // Set the pointers to bodyB
    TacsScalar *resB = &res[8];
    const TacsScalar *uB = &vars[8];
    const TacsScalar etaB = vars[11];
    const TacsScalar *epsB = &vars[12];

    // Compute the rotation matrix for bodyB
    TacsScalar CB[9];
    computeRotationMat(etaB, epsB, CB);

    // Get the initial position for bodyB
    TACSGibbsVector *rBVec = bodyB->getInitPosition();
    const TacsScalar *rB;
    rBVec->getVector(&rB);

    // Get the joint location from body B    
    const TacsScalar *xB;
    xBVec->getVector(&xB);

    // Position constraint on body B: dB = CB^{T}*xB + uB + rB0
    TacsScalar dB[3];  
    matMultTrans(CB, xB, dB);
    vecAxpy(1.0, uB, dB);
    vecAxpy(1.0, rB, dB);

    // Rotate the joint axes into different frames
    TacsScalar tB1[3], tB2[3];
    matMultTrans(CB, eB1, tB1);
    matMultTrans(CB, eB2, tB2);

    // Complete the evaluation of the constraint 1
    resC[0] += vecDot(tB1, dA);
    resC[1] += vecDot(tB2, dA);
    resC[0] -= vecDot(tB1, dB);
    resC[1] -= vecDot(tB2, dB);

    // Add the contributions of this constraint into body A    
    TacsScalar f[3];
    crossProduct(1.0, eA, lam, f);

    vecAxpy(1.0, f, &resA[0]);
    addEMatTransProduct(1.0, xA, f, etaA, epsA, 
                        &resA[3], &resA[4]);

    // Add the contributions of this constraint into body B
    vecAxpy(-1.0, f, &resB[0]);
    addEMatTransProduct(-1.0, xB, f, etaB, epsB, 
                        &resB[3], &resB[4]);

    // Compute the contributions to the first cylindrical constraint
    resC[3] += vecDot(tA, tB1);

    // Add the derivative of this constraint 
    // (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB to body A
    addEMatTransProduct(lam[3], eA, tB1, etaA, epsA, 
                        &resA[3], &resA[4]);

    // Add the derivative of this constraint 
    // d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA to body B
    addEMatTransProduct(lam[3], eB1, tA, etaB, epsB,
                        &resB[3], &resB[4]);

    // Compute the contributions to the second cylindrical constraint
    resC[4] += vecDot(tA, tB2);

    // Add the derivative (d(CA^{T}*eA)/dqA)*tB2 = E(eA)*tB2
    addEMatTransProduct(lam[4], eA, tB2, etaA, epsA, 
                        &resA[3], &resA[4]);
    // Add the derivative d(CB^{T}*eB2)/dqB)*tA = E(eB2)*tA
    addEMatTransProduct(lam[4], eB2, tA, etaB, epsB,
                        &resB[3], &resB[4]);
  }
  else {
    // This term should take us back to the origin if bodyB is not present
    vecAxpy(-1.0, pt, dA);

    // Complete the evaluation of the constraint 1
    resC[0] += vecDot(eB1, dA);
    resC[1] += vecDot(eB2, dA);

    // Add the contributions of this constraint into body A    
    TacsScalar f[3];
    crossProduct(1.0, tA, lam, f);

    vecAxpy(1.0, f, &resA[0]);
    addEMatTransProduct(1.0, xA, f, etaA, epsA, 
                        &resA[3], &resA[4]);

    // Compute the contributions to the first cylindrical constraint
    resC[3] += vecDot(tA, eB1);

    // Add the derivative (d(CA^{T}*eA)/dqA)*eB1 = E(eA)*tB
    addEMatTransProduct(lam[3], eA, eB1, etaA, epsA, 
                        &resA[3], &resA[4]);

    // Compute the contributions to the second cylindrical constraint
    resC[4] += vecDot(tA, eB2);

    // Add the derivative (d(CA^{T}*eA)/dqA)*eB2 = E(eA)*eB2
    addEMatTransProduct(lam[4], eA, eB2, etaA, epsA, 
                        &resA[3], &resA[4]);
  }

  // Add the dummy constraints for the remaining Lagrange multiplier
  // variables
  resC[2] += lam[2];
  for ( int i = 5; i < 8; i++ ){
    resC[i] += lam[i];
  }
}

/*
  Set the design variable values
*/
void TACSCylindricalConstraint::setDesignVars( const TacsScalar dvs[], 
                                               int numDVs ){
  point->setDesignVars(dvs, numDVs);
  updatePoints();
}

/*
  Get the design variable values associated with the joint location
*/
void TACSCylindricalConstraint::getDesignVars( TacsScalar dvs[], int numDVs ){
  point->getDesignVars(dvs, numDVs);
}


/*
  Copy the data and increment the reference counts
*/
TACSAverageConstraint::TACSAverageConstraint( TACSRigidBody *_bodyA,
                                              TACSGibbsVector *_point,
                                              TACSRefFrame *_refFrame,
                                              int _moment_flag ){
  moment_flag = _moment_flag;
  bodyA = _bodyA;
  bodyA->incref();
  point = _point;
  point->incref();
  refFrame = _refFrame;
  refFrame->incref();
}

/*
  Free the data
*/
TACSAverageConstraint::~TACSAverageConstraint(){
  bodyA->decref();
  point->decref();
  refFrame->decref();
}

const char *TACSAverageConstraint::elem_name = "TACSAverageConstraint";

/*
  The number of DOF per node
*/
int TACSAverageConstraint::numDisplacements(){
  return 8;
}

/*
  The number of nodes for the element
*/
int TACSAverageConstraint::numNodes(){
  return 5;
}

/*
  The element name
*/
const char* TACSAverageConstraint::elementName(){
  return elem_name;
}

/*
  Compute the kinetic and potential energies
*/
void TACSAverageConstraint::computeEnergies( double time,
                                             TacsScalar *_Te, 
                                             TacsScalar *_Pe,
                                             const TacsScalar Xpts[],
                                             const TacsScalar vars[],
                                             const TacsScalar dvars[] ){
  *_Te = 0.0;
  *_Pe = 0.0;
}

/*
  Add the residual of the governing equations
*/
void TACSAverageConstraint::addResidual( double time, 
                                         TacsScalar res[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[] ){
  // Get the initial position for bodyA  from global origin
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Retrieve the reference point location
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the reference point location relative to the intial body
  // point location
  TacsScalar bref[3];
  bref[0] = pt[0] - rA[0];
  bref[1] = pt[1] - rA[1];
  bref[2] = pt[2] - rA[2];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];
  
  // Compute the rotation matrix for the body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // The Lagrange multipliers and constraint pointers
  TacsScalar *con = &res[8*4];
  const TacsScalar *lam = &vars[8*4];

  // Get the quadrature points and weights
  const double *gaussPts, *gaussWts;
  int numGauss = FElibrary::getGaussPtsWts(3, &gaussPts, &gaussWts);

  // Perform the numerical quadrature
  for ( int i = 0; i < numGauss; i++ ){
    // Get the quadrature point
    const double xi = gaussPts[i];
    
    // Compute the shape functions
    double N[3];
    N[0] = -0.5*xi*(1.0 - xi);
    N[1] = (1.0 - xi)*(1.0 + xi);
    N[2] = 0.5*(1.0 + xi)*xi;

    double Na[3];
    Na[0] = -0.5 + xi;
    Na[1] = -2.0*xi;
    Na[2] = 0.5 + xi;

    // Compute the position and displacement vector
    TacsScalar Xp[3];
    Xp[0] = N[0]*Xpts[3] + N[1]*Xpts[6] + N[2]*Xpts[9];
    Xp[1] = N[0]*Xpts[4] + N[1]*Xpts[7] + N[2]*Xpts[10];
    Xp[2] = N[0]*Xpts[5] + N[1]*Xpts[8] + N[2]*Xpts[11];

    TacsScalar up[3];
    up[0] = N[0]*vars[8]  + N[1]*vars[16] + N[2]*vars[24];
    up[1] = N[0]*vars[9]  + N[1]*vars[17] + N[2]*vars[25];
    up[2] = N[0]*vars[10] + N[1]*vars[18] + N[2]*vars[26];

    // Compute the derivative of the position vector along the length
    // of the edge
    TacsScalar Xa[3];
    Xa[0] = Na[0]*Xpts[3] + Na[1]*Xpts[6] + Na[2]*Xpts[9];
    Xa[1] = Na[0]*Xpts[4] + Na[1]*Xpts[7] + Na[2]*Xpts[10];
    Xa[2] = Na[0]*Xpts[5] + Na[1]*Xpts[8] + Na[2]*Xpts[11];

    // Compute the position vector on the element surface relative to
    // the initial reference point
    TacsScalar Xref[3];
    Xref[0] = Xp[0] - rA[0] - bref[0];
    Xref[1] = Xp[1] - rA[1] - bref[1];
    Xref[2] = Xp[2] - rA[2] - bref[2];

    // Evaluate the displacement in the body-fixed coordinate frame:
    TacsScalar atmp[3];
    atmp[0] = Xp[0] + up[0] - uA[0] - rA[0];
    atmp[1] = Xp[1] + up[1] - uA[1] - rA[1];
    atmp[2] = Xp[2] + up[2] - uA[2] - rA[2];

    // uref = CA*(xp + up - uA - rA) - bref - Xref
    TacsScalar uref[3];
    matMult(CA, atmp, uref);
    uref[0] = uref[0] - bref[0] - Xref[0];
    uref[1] = uref[1] - bref[1] - Xref[1];
    uref[2] = uref[2] - bref[2] - Xref[2];

    // Compute the quadrature weight for this point
    TacsScalar h = sqrt(vecDot(Xa, Xa))*gaussWts[i];

    // Get the reference axis associated with 
    const TacsScalar *Cref;
    refFrame->getRotation(&Cref);

    // Compute the displacements in the local frame
    TacsScalar x[3], u[3];
    matMult(Cref, uref, u);
    
    // Add the integration along the coordinate directions
    con[0] += h*u[0];
    con[1] += h*u[1];
    con[2] += h*u[2];

    // Add the multipliers times the derivative of the constraints
    // w.r.t. the state variables
    for ( int j = 0; j < 3; j++ ){
      // Iterate over each displacement component and add the
      // contributions from the displacement degrees of freedom.
      for ( int k = 0; k < 3; k++ ){
        TacsScalar d[3], t[3];
        d[0] = N[j]*CA[k];
        d[1] = N[j]*CA[3+k];
        d[2] = N[j]*CA[6+k];
        matMult(Cref, d, t);
    
        // Add the derivative from the elastic degrees of freedom
        res[8*(j+1)+k] += h*(t[0]*lam[0] + t[1]*lam[1] + t[2]*lam[2]);

        // Add the contribution from the rigid degrees of freedom
        res[k] -= h*(t[0]*lam[0] + t[1]*lam[1] + t[2]*lam[2]);
      }
    }

    // Add the contributions to the derivative w.r.t. the quaternion
    // parameterization. Compute the transpose of the derivative of
    // (h*lam^{T}*Cref*Cbi*uref)
    TacsScalar s[3];
    matMultTrans(Cref, &lam[0], s);
    addDMatTransProduct(h, atmp, s, etaA, epsA, &res[3], &res[4]);

    if (moment_flag){
      // Evaluate the position along the first and second directions
      // in the body-fixed coordinate system
      matMult(Cref, Xref, x);

      if (moment_flag & X_MOMENT){
        con[3] += h*(x[1]*u[2] - x[2]*u[1]);
      }
      if (moment_flag & Y_MOMENT){
        con[4] += h*x[1]*u[0];
      }
      if (moment_flag & Z_MOMENT){
        con[5] += h*x[2]*u[0];
      }

      // Add the multipliers times the derivative of the constraints
      // w.r.t. the state variables
      for ( int j = 0; j < 3; j++ ){
        // Iterate over each displacement component and add the
        // contributions from the displacement degrees of freedom.
        for ( int k = 0; k < 3; k++ ){
          TacsScalar d[3], t[3];
          d[0] = N[j]*CA[k];
          d[1] = N[j]*CA[3+k];
          d[2] = N[j]*CA[6+k];
          matMult(Cref, d, t);

          // Add the derivative from the elastic degrees of freedom
          TacsScalar r = 0.0;
          if (moment_flag & X_MOMENT){
            r += (x[1]*t[2] - x[2]*t[1])*lam[3];
          }
          if (moment_flag & Y_MOMENT){
            r += x[1]*t[0]*lam[4];
          }
          if (moment_flag & Z_MOMENT){
            r += x[2]*t[0]*lam[5];
          }
    
          // Add the terms to the matrix
          res[8*(j+1)+k] += h*r;
          res[k] -= h*r;
        }
      }

      // Compute t = B^{T}*lam
      TacsScalar t[3];
      t[0] = t[1] = t[2] = 0.0;
      if (moment_flag & X_MOMENT){
        t[1] =-x[2]*lam[3];
        t[2] = x[1]*lam[3];
      }
      if (moment_flag & Y_MOMENT){
        t[0] = x[1]*lam[4];
      }
      if (moment_flag & Z_MOMENT){
        t[0] += x[2]*lam[5];
      }
      matMultTrans(Cref, t, s);
      addDMatTransProduct(h, atmp, s, etaA, epsA, &res[3], &res[4]);

      // Set dummy constraints
      if (!(moment_flag & X_MOMENT)){
        con[3] += lam[3];
      }
      if (!(moment_flag & Y_MOMENT)){
        con[4] += lam[4];
      }
      if (!(moment_flag & Z_MOMENT)){
        con[5] += lam[5];
      }
      con[6] += lam[6];
      con[7] += lam[7];
    }
    else {
      con[3] += lam[3];
      con[4] += lam[4];
      con[5] += lam[5];
      con[6] += lam[6];
      con[7] += lam[7];
    }
  }
}

void TACSAverageConstraint::addJacobian( double time, 
                                         TacsScalar J[],
                                         double alpha, 
                                         double beta, 
                                         double gamma,
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[] ){
  const int nvars = 5*8;

  // Get the initial position for bodyA  from global origin
  TACSGibbsVector *rAVec = bodyA->getInitPosition();
  const TacsScalar *rA;
  rAVec->getVector(&rA);

  // Retrieve the reference point location
  const TacsScalar *pt;
  point->getVector(&pt);

  // Compute the reference point location relative to the intial body
  // point location
  TacsScalar bref[3];
  bref[0] = pt[0] - rA[0];
  bref[1] = pt[1] - rA[1];
  bref[2] = pt[2] - rA[2];

  // Set the variables for body A
  const TacsScalar *uA = &vars[0];
  const TacsScalar etaA = vars[3];
  const TacsScalar *epsA = &vars[4];
  
  // Compute the rotation matrix for the body
  TacsScalar CA[9];
  computeRotationMat(etaA, epsA, CA);

  // The Lagrange multipliers and constraint pointers
  const TacsScalar *lam = &vars[8*4];

  // Get the quadrature points and weights
  const double *gaussPts, *gaussWts;
  int numGauss = FElibrary::getGaussPtsWts(3, &gaussPts, &gaussWts);

  // Perform the numerical quadrature
  for ( int i = 0; i < numGauss; i++ ){
    // Get the quadrature point
    const double xi = gaussPts[i];
    
    // Compute the shape functions
    double N[3];
    N[0] = -0.5*xi*(1.0 - xi);
    N[1] = (1.0 - xi)*(1.0 + xi);
    N[2] = 0.5*(1.0 + xi)*xi;

    double Na[3];
    Na[0] = -0.5 + xi;
    Na[1] = -2.0*xi;
    Na[2] = 0.5 + xi;

    // Compute the position and displacement vector
    TacsScalar Xp[3];
    Xp[0] = N[0]*Xpts[3] + N[1]*Xpts[6] + N[2]*Xpts[9];
    Xp[1] = N[0]*Xpts[4] + N[1]*Xpts[7] + N[2]*Xpts[10];
    Xp[2] = N[0]*Xpts[5] + N[1]*Xpts[8] + N[2]*Xpts[11];

    TacsScalar up[3];
    up[0] = N[0]*vars[8]  + N[1]*vars[16] + N[2]*vars[24];
    up[1] = N[0]*vars[9]  + N[1]*vars[17] + N[2]*vars[25];
    up[2] = N[0]*vars[10] + N[1]*vars[18] + N[2]*vars[26];

    // Compute the derivative of the position vector along the length
    // of the edge
    TacsScalar Xa[3];
    Xa[0] = Na[0]*Xpts[3] + Na[1]*Xpts[6] + Na[2]*Xpts[9];
    Xa[1] = Na[0]*Xpts[4] + Na[1]*Xpts[7] + Na[2]*Xpts[10];
    Xa[2] = Na[0]*Xpts[5] + Na[1]*Xpts[8] + Na[2]*Xpts[11];

    // Compute the position vector on the element surface relative to
    // the initial reference point
    TacsScalar Xref[3];
    Xref[0] = Xp[0] - rA[0] - bref[0];
    Xref[1] = Xp[1] - rA[1] - bref[1];
    Xref[2] = Xp[2] - rA[2] - bref[2];

    // Evaluate the displacement in the body-fixed coordinate frame:
    TacsScalar atmp[3];
    atmp[0] = Xp[0] + up[0] - uA[0] - rA[0];
    atmp[1] = Xp[1] + up[1] - uA[1] - rA[1];
    atmp[2] = Xp[2] + up[2] - uA[2] - rA[2];

    // uref = CA*(xp + up - uA - rA) - bref - Xref
    TacsScalar uref[3];
    matMult(CA, atmp, uref);
    uref[0] = uref[0] - bref[0] - Xref[0];
    uref[1] = uref[1] - bref[1] - Xref[1];
    uref[2] = uref[2] - bref[2] - Xref[2];

    // Compute the quadrature weight for this point
    TacsScalar h = alpha*sqrt(vecDot(Xa, Xa))*gaussWts[i];

    // Get the reference axis associated with 
    const TacsScalar *Cref;
    refFrame->getRotation(&Cref);

    // Compute the displacements in the local frame
    TacsScalar u[3];
    matMult(Cref, uref, u);

    // Evaluate the position along the first and second directions
    // in the body-fixed coordinate system
    TacsScalar x[3];
    matMult(Cref, Xref, x);

    // Compute s = Cref*lambda
    TacsScalar s0[3], s1[3];
    s0[0] = s0[1] = s0[2] = 0.0;
    if (moment_flag & X_MOMENT){
      s0[1] =-x[2]*lam[3];
      s0[2] = x[1]*lam[3];      
    }
    if (moment_flag & Y_MOMENT){
      s0[0] = x[1]*lam[4];
    }
    if (moment_flag & Z_MOMENT){
      s0[0] += x[2]*lam[5];
    }
    matMultTrans(Cref, s0, s1);
    matMultTrans(Cref, &lam[0], s0);

    // Add the multipliers times the derivative of the constraints
    // w.r.t. the state variables
    for ( int j = 0; j < 3; j++ ){
      // Iterate over each displacement component and add the
      // contributions from the displacement degrees of freedom.
      for ( int k = 0; k < 3; k++ ){
        TacsScalar d[3], t[3];
        d[0] = N[j]*CA[k];
        d[1] = N[j]*CA[3+k];
        d[2] = N[j]*CA[6+k];
        matMult(Cref, d, t);
    
        // Add the derivative from the elastic degrees of freedom
        J[(8*(j+1)+k)*nvars + 32] += h*t[0];
        J[(8*(j+1)+k)*nvars + 33] += h*t[1];
        J[(8*(j+1)+k)*nvars + 34] += h*t[2];

        J[32*nvars + (8*(j+1)+k)] += h*t[0];
        J[33*nvars + (8*(j+1)+k)] += h*t[1];
        J[34*nvars + (8*(j+1)+k)] += h*t[2];

        // Add the contribution from the rigid degrees of freedom
        J[k*nvars + 32] -= h*t[0];
        J[k*nvars + 33] -= h*t[1];
        J[k*nvars + 34] -= h*t[2];

        J[32*nvars + k] -= h*t[0];
        J[33*nvars + k] -= h*t[1];
        J[34*nvars + k] -= h*t[2];
      }
      
      // Add the term from the transpose of the derivative 
      // of the constraints times the multipliers
      addBlockEMat(h*N[j], etaA, epsA, s0, 
                   &J[(8*(j+1))*nvars + 3], nvars);
      addBlockEMatTrans(h*N[j], etaA, epsA, s0, 
                        &J[3*nvars + 8*(j+1)], nvars);
    }

    // Add the contributions to the derivative w.r.t. the quaternion
    // parameterization. Compute the transpose of the derivative of
    // (h*lam^{T}*Cref*Cbi*uref)

    // Add the derivative w.r.t. the multipliers
    addBlockDMatTrans(h, etaA, epsA, atmp,
                      &J[3*nvars + 4*8], nvars);

    // Add the derivative w.r.t. the
    addBlockEMat(-h, etaA, epsA, s0, &J[3], nvars);

    // Add the derivative of D(atmp)^{T}*s w.r.t. qA
    addBlockDMatTransDeriv(h, atmp, s0, &J[3*nvars + 3], nvars);

    // Add the derivative of D(atmp)^{T}*s w.r.t. uA
    addBlockEMatTrans(-h, etaA, epsA, s0, &J[3*nvars], nvars);

    // Add the derivative of the constraint w.r.t. the quaternions
    addBlockDMat(h, etaA, epsA, atmp, &J[32*nvars + 3], nvars);

    if (moment_flag){
      // Add the multipliers times the derivative of the constraints
      // w.r.t. the state variables
      for ( int j = 0; j < 3; j++ ){
        // Iterate over each displacement component and add the
        // contributions from the displacement degrees of freedom.
        for ( int k = 0; k < 3; k++ ){
          TacsScalar d[3], t[3];
          d[0] = N[j]*CA[k];
          d[1] = N[j]*CA[3+k];
          d[2] = N[j]*CA[6+k];
          matMult(Cref, d, t);
    
          // Add the derivative from the elastic degrees of freedom
          if (moment_flag & X_MOMENT){
            J[(8*(j+1)+k)*nvars + 35] += h*(x[1]*t[2] - x[2]*t[1]);
            J[35*nvars + (8*(j+1)+k)] += h*(x[1]*t[2] - x[2]*t[1]);
            J[k*nvars + 35] -= h*(x[1]*t[2] - x[2]*t[1]);
            J[35*nvars + k] -= h*(x[1]*t[2] - x[2]*t[1]);
          }
          if (moment_flag & Y_MOMENT){
            J[(8*(j+1)+k)*nvars + 36] += h*x[1]*t[0];
            J[36*nvars + (8*(j+1)+k)] += h*x[1]*t[0];
            J[k*nvars + 36] -= h*x[1]*t[0];
            J[36*nvars + k] -= h*x[1]*t[0];
          }
          if (moment_flag & Z_MOMENT){
            J[(8*(j+1)+k)*nvars + 37] += h*x[2]*t[0];
            J[37*nvars + (8*(j+1)+k)] += h*x[2]*t[0];
            J[k*nvars + 37] -= h*x[2]*t[0];
            J[37*nvars + k] -= h*x[2]*t[0];
          }
        }

        // Add the term from the transpose of the derivative 
        // of the constraints times the multipliers
        addBlockEMat(h*N[j], etaA, epsA, s1, 
                     &J[(8*(j+1))*nvars + 3], nvars);
        addBlockEMatTrans(h*N[j], etaA, epsA, s1, 
                          &J[3*nvars + 8*(j+1)], nvars);
      }

      // Add the derivative w.r.t. uA
      addBlockEMatTrans(-h, etaA, epsA, s1, &J[3*nvars], nvars);

      // Add the derivative contribution to uA
      addBlockDMatTransDeriv(h, atmp, s1, &J[3*nvars + 3], nvars);

      // Add the derivative w.r.t. the displacement variables
      addBlockEMat(-h, etaA, epsA, s1, &J[3], nvars);

      // Add the derivative w.r.t. the multipliers
      TacsScalar D[12];
      computeDMat(etaA, epsA, atmp, D);
      for ( int k = 0; k < 4; k++ ){
        if (moment_flag & X_MOMENT){
          J[(3+k)*nvars + 35] += h*(x[1]*D[8+k] - x[2]*D[4+k]);
          J[35*nvars + 3+k] += h*(x[1]*D[8+k] - x[2]*D[4+k]);
        }
        if (moment_flag & Y_MOMENT){
          J[(3+k)*nvars + 36] += h*x[1]*D[k];
          J[36*nvars + 3+k] += h*x[1]*D[k];
        }
        if (moment_flag & Z_MOMENT){
          J[(3+k)*nvars + 37] += h*x[2]*D[k];
          J[37*nvars + 3+k] += h*x[2]*D[k];
        }
      }

      if (!(moment_flag & X_MOMENT)){
        J[(nvars-5)*(nvars+1)] += alpha;
      }
      if (!(moment_flag & Y_MOMENT)){
        J[(nvars-4)*(nvars+1)] += alpha;
      }
      if (!(moment_flag & Z_MOMENT)){
        J[(nvars-3)*(nvars+1)] += alpha;
      }
      J[(nvars-2)*(nvars+1)] += alpha;
      J[(nvars-1)*(nvars+1)] += alpha;
    }
    else {
      for ( int k = 3; k < 8; k++ ){
        J[(8*4 + k)*(nvars+1)] += alpha;
      }
    }
  }
}
