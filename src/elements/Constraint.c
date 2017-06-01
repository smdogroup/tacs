#include "RigidBody.h"
#include "Constraint.h"
#include "TACSElementAlgebra.h"

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

  // Evaluate the constraint that enforces the position rA + CA^{T}*xA + uA
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

    // Add the derivative of this constraint to (d(CA^{T}*eA)/dqA)*tB1 = E(eA)*tB body A
    addEMatTransProduct(lam[3], eA, tB1, etaA, epsA, 
                        &resA[3], &resA[4]);
    // Add the derivative of this constraint to  d(CB^{T}*eB1)/dqB)*tA = E(eB1)*tA body B
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
