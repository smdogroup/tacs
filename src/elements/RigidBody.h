#ifndef TACS_RIGID_BODY_DYNAMICS_H
#define TACS_RIGID_BODY_DYNAMICS_H

/*
  Rigid-body dynamics routines for TACS

  Copyright (c) 2015-2016 Graeme Kennedy. All rights reserved. 
*/

#include "TACSElement.h"
#include "TACSGibbVector.h"

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
	       const TacsScalar _J[],
	       TACSRefFrame *_CInit,
	       TACSGibbsVector *_rInit, 
	       TACSGibbsVector *_vInit, 
	       TACSGibbsVector *_omegaInit, 
	       TACSGibbsVector *_gvec );
  ~TACSRigidBody();

  // Set design variables numbers associated with the inertial props.
  // ----------------------------------------------------------------
  void setDesignVarNums( int _massDV, const int _cDV[], 
			 const int _JDV[] );

  // Return the number of displacements and nodes
  // --------------------------------------------
  int numDisplacements(){ return 8; }
  int numNodes(){ return 1; }

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

  // Set and retrieve design variable values
  // ---------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[], int numDVs );

















  // Set and retrieve the variables
  // ------------------------------
  void setVariables( const TacsScalar qkin[], 
		     const TacsScalar qdyn[] );
  void getVariables( const TacsScalar *_r0[], 
		     const TacsScalar *_v0[],
		     const TacsScalar *_omega[] );
  void getInitVariables( TacsScalar qkin[],
			 TacsScalar qdyn[] );
  void getInitVariables( TACSRefFrame **_CInit,
			 TACSGibbsVector **_rInit, 
			 TACSGibbsVector **_vInit,
			 TACSGibbsVector **_omegaInit );

  // Routines for evaluating the transformation matrix
  // -------------------------------------------------
  void getRotation( const TacsScalar *C[] );
  void getRotationDeriv( const TacsScalar *C[],
			 const TacsScalar *Cd[],
			 int *nparam );
  void addRotationAdjResProduct( TacsScalar fdvSens[], int numDVs,
				 const TacsScalar psi[],
				 const TacsScalar phi[] );
  void testRotation( TacsScalar qkin[], double dh ); 
  
  // Retrieve the interial properties (in the body-fixed frame)
  // ----------------------------------------------------------
  void getInertial( TacsScalar *_mass, 
		    const TacsScalar *_c[], 
		    const TacsScalar *_J[] );
  void addInertialAdjResProduct( TacsScalar fdvSens[], int numDVs, 
				 TacsScalar dmdx, 
				 const TacsScalar dcdx[],
				 const TacsScalar dJdx[] );

  // Compute the kinematic and potential energies
  // --------------------------------------------
  void getEnergies( TacsScalar *Te, TacsScalar *Pe );

  // Add the residuals from the governing equations: R(q, qdot) = 0
  // --------------------------------------------------------------
  void getResidual( TacsScalar rkin[], TacsScalar rdyn[],
		    const TacsScalar qkinDot[], 
		    const TacsScalar qdynDot[] );

  // Add the linear combination of the Jacobians [dR/dq + alpha*dR/dqdot]
  // --------------------------------------------------------------------
  void getJacobian( TacsScalar D11[], TacsScalar D12[],
		    TacsScalar D21[], TacsScalar D22[], 
		    double alpha, double beta,
		    const TacsScalar qkinDot[],
		    const TacsScalar qdynDot[] );

  // Test the Jacobian implementation
  // --------------------------------
  void testJacobian( double alpha, double beta, 
		     const TacsScalar qkin[],
		     const TacsScalar qdyn[],
		     const TacsScalar qkinDot[],
		     const TacsScalar qdynDot[],
		     double dh );

  // Add the derivative of the residual w.r.t. mass/c/J
  // --------------------------------------------------
  void addAdjResProduct( TacsScalar fdvSens[], int numDVs,
			 const TacsScalar qkinDot[],
			 const TacsScalar qdynDot[],
			 const TacsScalar kinAdj[],
			 const TacsScalar dynAdj[] );
  void addAdjInitVariableProduct( TacsScalar fdvSens[], int numDVs,
				  const TacsScalar kinAdj[], 
				  const TacsScalar dynAdj[] );
  void testAdjResProduct( int numDVs, double dh );

 private:
  // The rotation matrix parametrization
  enum RotationParam rot_param;

  int nkin; // The number of kinematic variables
  static const int ndyn = 6; // The number of dynamic DOF

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

  // The position, velocity and angular velocity vectors
  TacsScalar r0[3], v0[3], omega[3];

  // The rotation matrix, the relative rotation matrix and the
  // derivative of the rotation matrix w.r.t. the parametrization
  // Note: C = Cr*CRef
  TacsScalar C[9], Cr[9], Cd[4*9];

  // The Euler parameter values - these are only used for the Euler
  // parameters
  TacsScalar p[4];

  // The Euler angles and the rate matrix - this is only used
  // for Euler angles
  TacsScalar S[9], Sd[3*9];
};

/*
  Base kinematic constraint class.

  All joint classes/boundary conditions are defined using this
  kinematic constraint class. The constraint class uses one or more
  bodies that are linked through kinematic constraints to other bodies
  or boundary conditions.  The kinematic constraints are formulated so
  that the forces and torques are expressed in the global reference
  frame.
*/
class TACSDynKinematicCon : public TACSObject {
 public:
  virtual ~TACSDynKinematicCon(){}

  // Set and retrieve design variable values
  // ---------------------------------------
  virtual void setDesignVars( const TacsScalar dvs[], int numDVs ) = 0;
  virtual void getDesignVars( TacsScalar dvs[], int numDVs ) = 0;

  // Retrieve the number of constraint force/torque dof
  // --------------------------------------------------
  virtual int getNumDof( int *_ncon ) = 0;

  // Get the rigid bodies associated with this constraint
  // ----------------------------------------------------
  virtual void getBodies( TACSDynBody **_bodyA,
			  TACSDynBody **_bodyB ) = 0;

  // Get the residual/Jacobians of the constraint equation
  // -----------------------------------------------------
  virtual void getResidual( TacsScalar rcon[], const TacsScalar fr[] ) = 0;
  virtual void getJacobian( TacsScalar D[], double alpha, 
			    const TacsScalar fr[] ) = 0;
  virtual void addAdjResProduct( TacsScalar fdvSens[], int numDVs,
  				 const TacsScalar resAdj[],
  				 const TacsScalar fr[] ) = 0;

  // Add the residual and Jacobian associated with the body
  // ------------------------------------------------------
  virtual void addBodyResidual( TacsScalar res[],
				TACSDynBody *body,
				const TacsScalar fr[] ) = 0;
  virtual void addBodyJacobian( TacsScalar D[], double alpha,
				TACSDynBody *body,
				const TacsScalar fr[] ) = 0;
  virtual void addBodyAdjResProduct( TacsScalar fdvSens[], int numDVs,
  				     const TacsScalar dynAdj[],
  				     TACSDynBody *body,
  				     const TacsScalar fr[] ) = 0;

  // Compute the off-diagonal Jacobian contributions
  // -----------------------------------------------
  virtual void getOffDiagJacobian( TacsScalar Dcon[], TacsScalar Ddyn[],
				   double alpha,
				   TACSDynBody *body,
				   const TacsScalar fr[] ) = 0;

  // Test the Jacobian of the kinematic constraints
  // ----------------------------------------------
  void testJacobian( double alpha, const TacsScalar fr[], double dh );
  void testAdjResProduct( int numDVs, double dh );
};

/*
  A spherical joint between either a fixed point and a body or two
  bodies.

  The spherical joint does not transmit a torque between the two
  connected bodies (or fixed point). The body is free to rotate in an
  arbitrary direction without restriction.
*/
class TACSDynSphericalJoint : public TACSDynKinematicCon {
 public:
  TACSDynSphericalJoint( TACSDynBody *_body1,
			 TACSDynBody *_body2,
			 TACSGibbsVector *_point );
  TACSDynSphericalJoint( TACSDynBody *_body1, 
			 TACSGibbsVector *_point );
  ~TACSDynSphericalJoint();

  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Retrieve the number of constraint force/torque dof
  // --------------------------------------------------
  int getNumDof( int *_ncon ){
    *_ncon = 6;
  }

  // Get the rigid bodies associated with this constraint
  // ----------------------------------------------------
  void getBodies( TACSDynBody **_bodyA,
		  TACSDynBody **_bodyB );

  // Get the residual/Jacobians of the constraint equation
  // -----------------------------------------------------
  void getResidual( TacsScalar res[], const TacsScalar fr[] );
  void getJacobian( TacsScalar D[], double alpha,
		    const TacsScalar fr[] );
  void addAdjResProduct( TacsScalar fdvSens[], int numDVs,
  			 const TacsScalar resAdj[],
  			 const TacsScalar fr[] );

  // Add the residual and Jacobian associated with the body
  // ------------------------------------------------------
  void addBodyResidual( TacsScalar rdyn[],
			TACSDynBody *body,
			const TacsScalar fr[] );
  void addBodyJacobian( TacsScalar D[], double alpha,
			TACSDynBody *body,
			const TacsScalar fr[] );
  void addBodyAdjResProduct( TacsScalar fdvSens[], int numDVs,
			     const TacsScalar dynAdj[],
			     TACSDynBody *body,
			     const TacsScalar fr[] );

  // Compute the off-diagonal Jacobian contributions
  // -----------------------------------------------
  void getOffDiagJacobian( TacsScalar Dcon[], TacsScalar Ddyn[],
			   double alpha, TACSDynBody *body,
			   const TacsScalar fr[] );

  // Get the points associated with the joint
  // ----------------------------------------
  void getFixedPoints( const TacsScalar *_xA[],
		       const TacsScalar *_xB[],
		       const TacsScalar *_xI[] ){
    *_xA = xA;
    *_xB = xB;
    *_xI = xI;
  }

 private:
  // Update the vector data stored locally in the object
  void updatePoints();

  // The number of constraints = number of forces + torques
  static const int ncon = 6;

  TACSDynBody *bodyA, *bodyB; // The rigid bodies
  TACSGibbsVector *point; // The point where the bodies are fixed

  TacsScalar xI[3]; // Fixed point/fixed separation
  TacsScalar xA[3]; // Position with first body-frame
  TacsScalar xB[3]; // Position with second body-frame
};

/*
  A revolute joint between either two bodies or a body and a fixed
  point. 

  The revolute joint permits a rotation about a specified axis in each
  body without torque. The reaction torque is confined to lie in the
  plane perpendicular to the revolute direction. The kinematic
  constraints for the revolute joint consist of a conincident point
  condition as well as a condition that the revolute direction
  expressed in one body frame remain perpendicular to two mutually
  orthogonal directions defined in the second coordinate frame.

  Given the body-fixed revolute direction revA in bodyA, and the
  revolute direction revB in bodyB, the code first computes the
  directions eB1 and eB2:

  eB2 = crossProduct(revB, e)
  eB1 = crossProduct(eB2, revB)

  where e is the coordinate direction in B with minimal dot product
  with revB. The kinematic condition on the relative motion is then:

  rA + rA1 - rB - rB1 = 0
  dot(revA, eB1) = 0
  dot(revA, eB2) = 0

  with the dynamic constraint that the torque along the revolute
  direction must be zero:
  
  dot(revA, gr) = 0
*/
class TACSDynRevoluteJoint : public TACSDynKinematicCon {
 public:
  TACSDynRevoluteJoint( TACSDynBody *_bodyA, 
			TACSDynBody *_bodyB, 
			TACSGibbsVector *_point,
			TACSGibbsVector *_rev );
  ~TACSDynRevoluteJoint();

  // Set and get the design variable values
  // --------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );

  // Retrieve the number of constraint force/torque dof
  // --------------------------------------------------
  int getNumDof( int *_ncon ){
    *_ncon = ncon;
  }

  // Get the rigid bodies associated with this constraint
  // ----------------------------------------------------
  void getBodies( TACSDynBody **_bodyA,
		  TACSDynBody **_bodyB );

  // Get the residual/Jacobians of the constraint equation
  // -----------------------------------------------------
  void getResidual( TacsScalar res[], const TacsScalar fr[] );
  void getJacobian( TacsScalar D[], double alpha,
		    const TacsScalar fr[] );
  void addAdjResProduct( TacsScalar fdvSens[], int numDVs,
  			 const TacsScalar resAdj[],
  			 const TacsScalar fr[] );

  // Add the residual and Jacobian associated with the body
  // ------------------------------------------------------
  void addBodyResidual( TacsScalar rdyn[],
			TACSDynBody *body,
			const TacsScalar fr[] );
  void addBodyJacobian( TacsScalar D[], double alpha,
			TACSDynBody *body,
			const TacsScalar fr[] );
  void addBodyAdjResProduct( TacsScalar fdvSens[], int numDVs,
			     const TacsScalar dynAdj[],
			     TACSDynBody *body,
			     const TacsScalar fr[] );

  // Compute the off-diagonal Jacobian contributions
  // -----------------------------------------------
  void getOffDiagJacobian( TacsScalar Dcon[], TacsScalar Ddyn[],
			   double alpha, 
			   TACSDynBody *body,
			   const TacsScalar fr[] );

  // Get the points associated with the joint
  // ----------------------------------------
  void getFixedPoints( const TacsScalar *_xA[],
		       const TacsScalar *_xB[] ){
    *_xA = xA;
    *_xB = xB;
  }

  // Get the rotation vectors in the local frames
  // --------------------------------------------
  void getRevoluteAxis( const TacsScalar *_revA[],
			const TacsScalar *_revB[] ){
    *_revA = revA;
    *_revB = revB;
  }

 private:
  // Update the internal points used to evaluate the kinematic constraint
  void updatePoints( int init_ek = 0 );

  // The number of constraints = number of forces + torques
  static const int ncon = 6;

  // The bodies and vectors of interest
  TACSDynBody *bodyA, *bodyB; // The rigid bodies
  TACSGibbsVector *point, *rev; // The fixed point and revolute direction

  // The revolute direction in the first and second bodies
  TacsScalar revA[3], revB[3];

  // Two mutually perpendicular directions to the revolute direction
  // which are also perpendicular to each other such that they form a
  // triad that is right-handed
  TacsScalar eB1[3], eB2[3];
  TacsScalar ek[3]; // Fixed point/fixed separation

  // The fixed point location within the two bodies
  TacsScalar xA[3]; // Position with first body-frame
  TacsScalar xB[3]; // Position with second body-frame
};

#endif // TACS_RIGID_BODY_DYNAMICS_H
