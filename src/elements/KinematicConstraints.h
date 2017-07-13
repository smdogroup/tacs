#ifndef TACS_KINEMATIC_CONSTRAINTS_H
#define TACS_KINEMATIC_CONSTRAINTS_H

/*
  Constraints for Rigid-body dynamics in TACS.
  
  Copyright (c) 2015-2017 Graeme Kennedy. All rights reserved.
*/

#include "TACSElement.h"
#include "TACSGibbsVector.h"
#include "RigidBody.h"

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
  // Update the local data
  void updatePoints();
        
  // The rigid bodies involved in the joint
  TACSRigidBody *bodyA, *bodyB;

  // The point where the joint is located in global frame
  TACSGibbsVector *point;

  // The positions of joint from each body in global frame
  TACSGibbsVector *xAVec, *xBVec;

  // The name of the element
  static const char *elem_name;
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
  TACSRevoluteConstraint( int _fixed_ref_point,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec );
  ~TACSRevoluteConstraint();

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
  // Update the local data
  void updatePoints( int init_e=0 );

  // Is the reference point fixed or not?
  int fixed_ref_point;

  // The rigid bodies involved in the joint
  TACSRigidBody *bodyA, *bodyB; 

  // Point where the joint is located in global frame
  TACSGibbsVector *point; 

  // Revolute direction in global frame
  TACSGibbsVector *eAVec; 
  
  // Local axis in the B-frame
  TACSGibbsVector *eB1Vec, *eB2Vec; 

  // The coordinate direction in global frame
  TACSGibbsVector *eVec; 

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

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    *multiplier = 1;
  }

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

/*
  The following constraint is designed to connect rigid and flexible
  elements in an average sense. The element enforces that the zeroth
  and first moments of the displacements about a reference point in a
  given frame are zero. Constraint elements must be added along the
  entire interface, and these constraints must share the same rigid
  body and multiplier nodes.  This approach is designed to avoid
  issues of non-physical, perfectly rigid connections between rigid
  and flexible components.

  The absolute position of the reference point in the inertial frame
  is computed as a function of the displacement and rotation of the
  rigid body as follows:

  ref(t = 0) = r0 + bref
  ref(t) = r0 + u0 + Cbi^{T}*bref

  where bref is the vector from the origin of the initial body
  position to the reference pooint.
  
  The position of the flexible point, anywhere along the surface, is
  given as

  pos(t = 0) = X(xi)
  pos(t) = X(xi) + u(xi)

  where xi is a coordinate that runs along the edge of the connecting
  element.

  Initial position vector from the reference point (in the inertial
  frame):
  
  Xref(t = 0) = X(xi) - r0 - bref
  
  Position vector at t is:
  
  Xref(t) = X(xi) + u(xi) - (r0 + u0 + Cbi^{T}*bref)

  The relative displacement of the position in the flexible body is
  given as:

  uref = u(xi) - u0 - (Cbi^{T} - I)*bref

  In the body-fixed coordinate frame, the displacements are:

  u = Cbi*uref 
  . = Cbi*(u(xi) - u0) + (Cbi - I)*bref
  . = Cbi*(u(xi) - u0 + bref) - bref

  We then form a reference frame Cref (from global to local) and
  compute

  ulocal = Cref*u, and xlocal = Cref*Xref

  and perform the integration in the local coordinate axes.
*/
class TACSAverageConstraint : public TACSElement {
 public:
  TACSAverageConstraint( TACSRigidBody *_bodyA,
                         TACSGibbsVector *_point,
                         TACSRefFrame *_refFrame,
                         int _use_moments );
  ~TACSAverageConstraint();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    *multiplier = 4;
  }

  // Get the number of displacements/nodes and the element name
  // ----------------------------------------------------------
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

  // Assemble the residual and the Jacobian
  // --------------------------------------
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
  // Flag to indicate whether to constraint the displacements or the
  // displacements and the moments of the displacement
  int use_moments;

  // The rigid body 
  TACSRigidBody *bodyA;

  // The point in the inertial reference frame
  TACSGibbsVector *point;

  // The reference frame used to define the local coordinate system in
  // the initial configuration (moments are taken about the y-z plane
  // in this reference frame)
  TACSRefFrame *refFrame;  

  // The element name
  static const char *elem_name;
};

#endif // TACS_KINEMATIC_CONSTRAINTS_H
