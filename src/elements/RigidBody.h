#ifndef TACS_RIGID_BODY_DYNAMICS_H
#define TACS_RIGID_BODY_DYNAMICS_H

/*
  Rigid-body dynamics routines for TACS

  Copyright (c) 2015-2016 Graeme Kennedy. All rights reserved. 
*/

#include "TACSElement.h"
#include "TACSGibbsVector.h"

/*
  A reference coordinate system that is built to enable a change in
  orientation of coordinate frames as a function of input vectors.

  This class builds a rotation matrix based on 3 input vectors in the
  global reference frame. These vectors represent the base point, r0,
  and two other points r1 and r2 which are used to construct the
  reference frame. The rotation matrix is used to transform from the
  global coordinate from to a local, body-fixed coordinate frame.  In
  the body-frame, (r1 - r0) is aligned along the first local
  coordinate vector and the orthogonal complement along (r2 - r0)
  forms the second basis vector.
*/
class TACSRefFrame : public TACSObject {
 public:
  TACSRefFrame( TACSGibbsVector *_r0, 
		TACSGibbsVector *_r1,
		TACSGibbsVector *_r2 );
  ~TACSRefFrame();
  void initialize();
  void getRotation( const TacsScalar ** _C );
  void setDesignVars( const TacsScalar *dvs, int numDVs );
  void getDesignVars( TacsScalar *dvs, int numDVs );
  void addRotationAdjResProduct( TacsScalar fdvSens[], int numDVs,
				 const TacsScalar psi[],
				 const TacsScalar phi[] );
  void testRotation( int numDVs, double dh );

 private:
  // The rotation matrix associated with this reference frame
  TacsScalar C[9];

  // Store the derivatives of the first and second rows of the
  // rotation matrix with respect to the vectors d1 = r1 - r0
  // and d2 = r2 - r0
  TacsScalar dC1d1[9];
  TacsScalar dC2d1[9], dC2d2[9]; 

  // The basis points for the reference frame
  TACSGibbsVector *r0, *r1, *r2;
};

/*
  Dynamics for a single rigid body.

  This class handles the dynamics of a single free body subject to
  external moments and torques and a constant gravity load (although
  this may change in the future.)

  In this class, all inertial properties (with the exception of the
  body mass of course) are expressed in a body-fixed reference frame.
  The body-fixed frame is defined by an initial reference
  configuration plus a rotation parametrized either using Euler angles
  or Euler parameters/quaternions.

  Kinematic constraints add additional internal reaction forces and
  torques are required to complete the full multibody system. These
  internal reactions are treated within the kinematic constraint class
  defined below.
*/
class TACSRigidBody : public TACSElement {
 public:
  TACSRigidBody( const TacsScalar _mass, 
                 const TacsScalar _c[], 
                 const TacsScalar _J[] );
  ~TACSRigidBody();

  // Set design variables numbers associated with the inertial props.
  // ----------------------------------------------------------------
  void setDesignVarNums( int _massDV, 
                         const int _cDV[], 
			 const int _JDV[] );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes(){ return 1; }

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[], int numDVs );

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitCondition( TacsScalar vars[],
			 TacsScalar dvars[],
			 const TacsScalar X[] );

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
  void getResidual( double time, TacsScalar res[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void getJacobian( double time, TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Test the residual implementation
  // --------------------------------
  void testResidual( double dh );

 private:
  TacsScalar mass; // The mass of the rigid body
  TacsScalar c[3]; // The first moment of inertia
  TacsScalar J[6]; // The second moment of inertia

  // The initial position, velocity, angular velocity
  TACSGibbsVector *rInit, *vInit, *omegaInit;

  // The initial reference frame
  TACSRefFrame *CRef;

  // The design variable numbers associated with the inertial properties
  int massDV, cDV[3], JDV[6];

  TacsScalar g[3]; // The accel. due to gravity in the global ref. frame
};

#endif // TACS_RIGID_BODY_DYNAMICS_H
