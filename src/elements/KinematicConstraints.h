#ifndef TACS_KINEMATIC_CONSTRAINTS_H
#define TACS_KINEMATIC_CONSTRAINTS_H

/*
  Constraints for Rigid-body dynamics in TACS.
  
  Copyright (c) 2015-2017 Graeme Kennedy. All rights reserved.
*/

#include "TACSElement.h"
#include "TACSGibbsVector.h"

/*
  Spherical constraint
*/
class TACSSphericalConstraint : public TACSElement {
 public:
  TACSSphericalConstraint( TACSRigidBody *_bodyA,
                           TACSRigidBody *_bodyB,
                           TACSGibbsVector *_point );
  TACSSphericalConstraint( TACSRigidBody *_bodyA,
                           TACSGibbsVector *_point );
  ~TACSSphericalConstraint();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    if (bodyA && bodyB){
      *multiplier = 2;
    }
    else {
      *multiplier = 1;
    }    
  }

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes();
  const char* elementName(){ return elem_name; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

 private:
  void updatePoints();            // Update the local data
  TACSRigidBody *bodyA, *bodyB;   // The rigid bodies involved in the joint
  TACSGibbsVector *point;         // The point where the joint is located in global frame
  TACSGibbsVector *xAVec, *xBVec; // The positions of joint from each body in global frame
  static const char *elem_name;   // The name of the element
};

/*
  Revolute constraint
*/
class TACSRevoluteConstraint : public TACSElement {
 public:
  TACSRevoluteConstraint( TACSRigidBody *_bodyA,
                          TACSRigidBody *_bodyB,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec );
  TACSRevoluteConstraint( TACSRigidBody *_bodyA,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec );
  ~TACSRevoluteConstraint();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    *multiplier = 2;
  }

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes();
  const char* elementName(){ return elem_name; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
 private:
  // Update the local data
  void updatePoints( int init_e=0 );

  TACSRigidBody *bodyA, *bodyB; // The rigid bodies involved in the joint
  TACSGibbsVector *point; // Point where the joint is located in global frame
  TACSGibbsVector *eAVec; // Revolute direction in global frame
  TACSGibbsVector *xAVec, *xBVec; // Joint ositions from body in global frame
  TACSGibbsVector *eB1Vec, *eB2Vec; 
  TACSGibbsVector *eVec; // The coordinate direction in global frame
  static const char *elem_name; // The name of the element
};

/*
  Spherical constraint
*/
class TACSRigidLink : public TACSElement {
 public:
  TACSRigidLink( TACSRigidBody *_bodyA );
  ~TACSRigidLink();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    *multiplier = 2;
  }

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements();
  int numNodes();
  const char* elementName();

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
 
 private:
  TACSRigidBody *bodyA; // The rigid body
  static const char *elem_name; // The name of the element
};

/*
  Drives the connected points at a specified angular rate about
  the specified revolute direction fixed at the given origin.
*/
class TACSRevoluteDriver : public TACSElement {
 public:
  TACSRevoluteDriver( TACSGibbsVector *orig, 
                      TACSGibbsVector *rev,
                      TacsScalar _omega );
  ~TACSRevoluteDriver();

  int numDisplacements();
  int numNodes();
  const char* elementName();

  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
  void addJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

 private:
  TacsScalar omega;
  TACSGibbsVector *origVec, *revVec;
};

/*
  Cylindrical constraint between rigid bodies
*/
class TACSCylindricalConstraint : public TACSElement {
 public:
  TACSCylindricalConstraint( TACSRigidBody *_bodyA,
                             TACSRigidBody *_bodyB,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  TACSCylindricalConstraint( TACSRigidBody *_bodyA,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  ~TACSCylindricalConstraint();

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes();
  const char* elementName(){ return elem_name; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );
 private:
  // Update the local data 
  void updatePoints( int init_e=0 );

  // The rigid bodies involved in the joint
  TACSRigidBody *bodyA, *bodyB;

  // The point where the joint is located in global frame
  TACSGibbsVector *point;

  // The revolute direction in global frame
  TACSGibbsVector *eAVec;

  // The positions of joint from each body in global frame
  TACSGibbsVector *xAVec, *xBVec;

  // The positions of joint from each body in global fram
  TACSGibbsVector *eB1Vec, *eB2Vec;

  // The coordinate direction in global frame with minimal dot product
  // with eAVec
  TACSGibbsVector *eVec;

  // The name of the element
  static const char *elem_name;
};

#endif // TACS_KINEMATIC_CONSTRAINTS_H
