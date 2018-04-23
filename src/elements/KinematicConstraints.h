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

#ifndef TACS_KINEMATIC_CONSTRAINTS_H
#define TACS_KINEMATIC_CONSTRAINTS_H

/*
  Constraints for Rigid-body dynamics in TACS.
*/

#include "TACSElement.h"
#include "TACSGibbsVector.h"
#include "RigidBody.h"

enum TACSTranslationConType { COINCIDENT, COLINEAR, COPLANAR };

/*
  Generic base class for constraint elements
*/
class TACSKinematicConstraint : public TACSElement {
 public:
  TACSKinematicConstraint( TACSGibbsVector *_point,
                           int _dof_flag,
                           TACSRigidBody *_bodyA,
                           TACSRigidBody *_bodyB=NULL );
  ~TACSKinematicConstraint();
  
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

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes() = 0;
  const char* elementName(){ return elem_name; }

 private:
  // The point where the joint is located in global frame
  TACSGibbsVector *point;

  // The rigid bodies involved in the joint
  TACSRigidBody *bodyA, *bodyB;
  
  // The positions of joint from each body in global frame
  TACSGibbsVector *xAVec, *xBVec;
  
  // The name of the element
  static const char *elem_name;
};

/*
  Spherical constraint

  The spherical constraint class imposes no restriction on the
  relative rotation between the two reference frames (either the first
  frame A and the inertial frame or the first frame A and a second
  frame B). 

  The translational degrees of freedom are constrained either to be
  1. Coincident 
  2. Colinear
  3. Coplanar

  The nature of these translational constraints is dictated by the
  arguments passed to the constructor.
*/
class TACSSphericalConstraint : public TACSElement {
 public:
  TACSSphericalConstraint( TACSRigidBody *_bodyA,
                           TACSRigidBody *_bodyB,
                           TACSGibbsVector *_point,
                           TACSGibbsVector *_axis=NULL,
                           TACSTranslationConType con_type=COINCIDENT );
  TACSSphericalConstraint( TACSRigidBody *_bodyA,
                           TACSGibbsVector *_point,
                           TACSGibbsVector *_axis=NULL, 
                           TACSTranslationConType con_type=COINCIDENT );
  ~TACSSphericalConstraint();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier );

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
  // Keep a record of the constraint type
  TACSTranslationConType con_type;

  // Indicates whether the point/constraint is fixed inertially (true)
  // or whether it is fixed within body A (false)
  int inertial_fixed_point;

  // The rigid bodies involved in the joint
  TACSRigidBody *bodyA, *bodyB;

  // The axis used to define either the axis or plane for the
  // constraints
  TACSGibbsVector *axis;

  // The point where the joint is located in global frame
  TACSGibbsVector *point;

  // The name of the element
  static const char *elem_name;
};

/*
  Revolute constraint

  This constraint forces the relative rotation between two bodies 
  (body A and body B) at a point to lie on a single axis. There are 
  several different ways that this constraint can be formualted. Body A
  is treated as the primary body, while body B is treated as a secondary
  body. When TACSRigidBody objects are passed in, the reference point
  location is retrieved from the body object, whereas if no body is 
  passed in, the point is taken from the node location.

  First, the revolute axis may be either fixed to the inertial reference
  frame or fixed/convected with body B's body-aligned frame. In the 
  latter case, the revolute axis changes continuously as a function of 
  the orientation of body B.

  Second, the reference point may be constrained such that it is:
  1. Fixed such that the point is fixed inertially
  2. Constrained so that the two components in frames A and B coincide
  3. Free - unconstrained
*/
class TACSRevoluteConstraint : public TACSElement {
 public:
  TACSRevoluteConstraint( TACSRigidBody *_bodyA,
                          TACSRigidBody *_bodyB,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec,
                          int _inertial_rev_axis=0 );
  TACSRevoluteConstraint( TACSRigidBody *_bodyA,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec );
  TACSRevoluteConstraint( int _fixed_ref_point,
                          TACSGibbsVector *_point,
                          TACSGibbsVector *_eAVec,
                          int _inertial_rev_axis=0 );
  ~TACSRevoluteConstraint();

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier );

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

  // Is the reference axis fixed in body B's body-fixed frame
  // or is it fixed in the inertial refernece frame
  int inertial_rev_axis;

  // Are there two bodies or just one?
  int inertial_fixed_point;

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

  // Parameter used to avoid linearly dependent Jacobian rows
  static const TacsScalar delta;
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

/*
  Drives the attached body sinusoidally along the given direction
*/
class TACSMotionDriver : public TACSElement {
 public:
  TACSMotionDriver( TACSGibbsVector *_dir, 
                    TacsScalar _omega,
                    int _fix_rotations=0 ){
    // Copy over the direction
    dir = _dir;
    dir->incref();
    
    // Copy over the angular rate
    omega = _omega;
    
    // Fix the rotations if requested
    fix_rotations = _fix_rotations;
  }

  ~TACSMotionDriver(){    
    dir->decref();
  }

  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
    *multiplier = 1;
  }

  int numDisplacements(){
    return 8;
  }

  int numNodes(){
    return 2;
  }

  const char* elementName(){
    return "TACSMotionDriver";
  }

  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] ){
    *_Te = 0.0;
    *_Pe = 0.0;
  }

  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Retrieve the direction
    const TacsScalar *d;
    dir->getVector(&d);
    
    // The Lagrange multipliers
    const TacsScalar *lam = &vars[8];
    
    // Specify a sinusoidal motion
    const TacsScalar scale = sin(omega*time);

    // Equations to specify sinusoidal translatory motion
    res[8] += vars[0] - scale*d[0]; 
    res[9] += vars[1] - scale*d[1];
    res[10] += vars[2] - scale*d[2];

    // Add the constraint reaction forces to the body
    res[0] += lam[0];
    res[1] += lam[1];
    res[2] += lam[2];

    if (fix_rotations){
      // Equations to arrest rotational DOF
      res[8+4] += vars[4];
      res[8+5] += vars[5];
      res[8+6] += vars[6];
      
      // Apply reaction moments to arrest rotational dof
      res[4] += lam[4];
      res[5] += lam[5];
      res[6] += lam[6];

      //Add dummy constraint eqns
      res[8+3] += lam[3];
      res[8+7] += lam[7];
    } 
    else {
      // Add the dummy constraints for remaining constraint equations
      for ( int i = 3; i < 8; i++ ){
        res[8+i] += lam[i];
      }
    }
  }
  
 private:
  TacsScalar omega;
  TACSGibbsVector *dir;
  int fix_rotations;
};

/*

  Drives the attached body sinusoidally along the given direction for bodies
  with six degrees of freedom (for linearized rotations)

*/

class TACSLinearizedMotionDriver : public TACSElement {
 public:
  TACSLinearizedMotionDriver( TACSGibbsVector *_dir, 
                    TacsScalar _omega,
                    int _arrest_rotations = 0){
    // Copy over the direction
    dir = _dir;
    dir->incref();
    
    // Copy over the angular rate
    omega = _omega;
    
    // Arrest rotations if requested
    arrest_rotations = _arrest_rotations;
  }

  ~TACSLinearizedMotionDriver(){
    dir->decref();
  }

  void getMultiplierIndex( int *multiplier ){
    *multiplier = 1;
  }

  int numDisplacements(){
    return 6;
  }

  int numNodes(){
    return 2;
  }

  const char* elementName(){
    return "TACSMotionDriver";
  }

  void computeEnergies( double time,
                        TacsScalar *_Te, 
                        TacsScalar *_Pe,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[] ){
    *_Te = 0.0;
    *_Pe = 0.0;
  }

  void addResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] ){
    // Retrieve the direction
    const TacsScalar *d;
    dir->getVector(&d);
    
    // The Lagrange multipliers
    const TacsScalar *lam = &vars[6];
    
    // Specify a sinusoidal motion
    const TacsScalar scale = sin(omega*time);

    // Equations to specify sinusoidal translatory motion
    res[6] += vars[0] - scale*d[0]; 
    res[7] += vars[1] - scale*d[1];
    res[8] += vars[2] - scale*d[2];
    res[9] += vars[3];
    res[10] += vars[4];
    res[11] += vars[5];

    // Add the constraint reaction forces to the body
    res[0] += lam[0];
    res[1] += lam[1];
    res[2] += lam[2];
    res[3] += lam[3];
    res[4] += lam[4];
    res[5] += lam[5];
  }
  
 private:
  TacsScalar omega;
  TACSGibbsVector *dir;
  int arrest_rotations;
};

/*
  Drives the connected points at a specified angular rate about
  the specified revolute direction fixed at the given origin.
*/
class TACSRevoluteDriver : public TACSElement {
 public:
  TACSRevoluteDriver( TACSGibbsVector *rev,
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

  // Parameter used to avoid linearly dependent Jacobian rows
  static const TacsScalar delta;
};



/*
  These will become deprecated...
*/




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
  Prismatic constraint between rigid bodies
*/
class TACSPrismaticConstraint : public TACSElement {
 public:
  TACSPrismaticConstraint( TACSRigidBody *_bodyA,
                             TACSRigidBody *_bodyB,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  TACSPrismaticConstraint( TACSRigidBody *_bodyA,
                             TACSGibbsVector *_point,
                             TACSGibbsVector *_eAVec );
  ~TACSPrismaticConstraint();

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
  Sliding pivot constraint
*/
class TACSSlidingPivotConstraint : public TACSElement {
 public:
  TACSSlidingPivotConstraint( TACSRigidBody *_bodyA,
                              TACSRigidBody *_bodyB,
                              TACSGibbsVector *_point,
                              TACSGibbsVector *_eAVec );
  TACSSlidingPivotConstraint( TACSRigidBody *_bodyA,
                              TACSGibbsVector *_point,
                              TACSGibbsVector *_eAVec );
  ~TACSSlidingPivotConstraint();

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
  Fixed constraint constrains all the degrees of freedom 
  of the attached body
*/
class TACSFixedConstraint : public TACSElement {
 public:
  TACSFixedConstraint( TACSRigidBody *_bodyA,                       
                       TACSGibbsVector *_point );
  ~TACSFixedConstraint();
  
  // Get the multiplier precedent to ensure they are ordered last
  // ------------------------------------------------------------
  void getMultiplierIndex( int *multiplier ){
      *multiplier = 1;
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
  void updatePoints();
        
  // The rigid bodies involved in the joint
  TACSRigidBody *body;

  // The point where the joint is located in global frame
  TACSGibbsVector *point;

  // The positions of joint from each body in global frame
  TACSGibbsVector *xVec;

  // The name of the element
  static const char *elem_name;
};

#endif // TACS_KINEMATIC_CONSTRAINTS_H
