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

  The absolute position of a point along the flexible connection is
  given by

  X(xi) + U(xi)

  This point on the flexible body is also observed from a frame fixed
  relative to the rigid body, called C. The average displacements and
  displacement moments observed in the C frame are set to zero. The
  position of the C frame within the rigid body frame is given by the
  vector bref, and is fixed in the body frame. At the initial point
  the body-frame and inertial frame are aligned so that bref at the
  initial point is also in the body frame.

  At t=0 the relative position of the point X(xi) on the flexible
  body relative to frame is given as

  Xref = X(xi) - r0 - bref

  Note that r0 is the inital body location. At time t, the flexible
  and rigid body have moved to point X(xi) + U(xi) and r0 + u0,
  respectively. But Xref has also been convected in frame C.  The new
  position of Xref in the inertial frame is given as:

  X' = r0 + u0 + CB^{T}(bref + Xref)
  
  The difference between X(xi) + U(xi) and X' is the displacement
  observed in frame C. This is then
  
  U' = X(xi) + U(xi) - X'
  .  = X(xi) + U(xi) - r0 - u0 - CB^{T}(bref + Xref)

  The position vector in the local frame C is then:

  u = Cref*CB*U'
  . = Cref*(CB*(X(xi) + U(xi) - r0 - u0) - bref + Xref)

  The integration is performed using the local Xref location:

  x = Cref*Xref.
*/
class TACSAverageConstraint : public TACSElement {
 public:
  static const int X_MOMENT = 1;
  static const int Y_MOMENT = 2;
  static const int Z_MOMENT = 4; 

  TACSAverageConstraint( TACSRigidBody *_bodyA,
                         TACSGibbsVector *_point,
                         TACSRefFrame *_refFrame,
                         int _moment_flag );
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
  int moment_flag;

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
