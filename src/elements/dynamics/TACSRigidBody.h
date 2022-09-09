/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_RIGID_BODY_DYNAMICS_H
#define TACS_RIGID_BODY_DYNAMICS_H

/*
  Rigid-body dynamics classes for TACS.

  1. TACSRefFrame
  2. TACSRigidBody
  3. TACSSphericalConstraint
  4. TACSRevoluteConstraint
*/

#include "TACSElement.h"
#include "TACSElementTypes.h"
#include "TACSGibbsVector.h"

/*
  Class for visualizing rigid bodies with different geometries
*/
class TACSRigidBodyViz : public TACSObject {
 public:
  TACSRigidBodyViz(int _npts, int _nelems, TacsScalar *_Xpt, int _conn[],
                   TACSGibbsVector *_vref = NULL);
  TACSRigidBodyViz(TacsScalar Lx, TacsScalar Ly, TacsScalar Lz);
  ~TACSRigidBodyViz();
  void getMesh(int *_npts, int *_nelems, const TacsScalar **_Xpts,
               const int **_conn);

 private:
  void init(int _npts, int _nelems, TacsScalar *_Xpt, int _conn[],
            TACSGibbsVector *_vref = NULL);

  TacsScalar *Xpts;
  int npts, nelems;
  int *conn;
};

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
  TACSRefFrame(TACSGibbsVector *_r0, TACSGibbsVector *_r1,
               TACSGibbsVector *_r2);
  ~TACSRefFrame();

  // Returns the rotation matrix associated with the frame of ref.
  //--------------------------------------------------------------
  void getRotation(const TacsScalar **_C);

  // Getters and setters for DVs
  //----------------------------
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Product of adjoint variables with rotation parametrization
  //-----------------------------------------------------------
  void addRotationAdjResProduct(const TacsScalar psi[], const TacsScalar phi[],
                                int dvLen, TacsScalar dfdx[]);

  // Routine to perform sanity checks on the implementations of the
  // frame of reference
  // --------------------------------------------------------------
  void testRotation(double dh);

 private:
  // Recompute the rotation matrix
  //------------------------------
  void initialize();

  // The rotation matrix associated with this reference frame
  TacsScalar C[9];

  // Store the derivatives of the first and second rows of the
  // rotation matrix with respect to the vectors d1 = r1 - r0
  // and d2 = r2 - r0
  TacsScalar dC1d1[9];
  TacsScalar dC2d1[9], dC2d2[9];

  // The basis points for the reference frame
  TACSGibbsVector *r0, *r1, *r2;
  int r0offset, r1offset, r2offset;  // Offset indices
};

/*
  Dynamics for a single rigid body.

  This class handles the dynamics of a single free body subject to
  external moments and torques and a constant gravity load (although
  this may change in the future.)

  In this class, all inertial properties (with the exception of the
  body mass of course) are expressed initially in a body-fixed
  reference frame, but are later converted into a Newtonian frame.
  The body-fixed frame is defined by an initial reference
  configuration plus a rotation parametrized either using Euler
  parameters/quaternions.

  Kinematic constraints add additional internal reaction forces and
  torques are required to complete the full multibody system. These
  internal reactions are treated within the kinematic constraint class.
*/
class TACSRigidBody : public TACSElement {
 public:
  TACSRigidBody(TACSRefFrame *_CRef, const TacsScalar _mass,
                const TacsScalar _cRef[], const TacsScalar _JRef[],
                TACSGibbsVector *_rInit, TACSGibbsVector *_vInit,
                TACSGibbsVector *_omegaInit, TACSGibbsVector *_gvec);
  ~TACSRigidBody();

  // Set design variables numbers associated with the inertial props.
  // ----------------------------------------------------------------
  void setDesignVarNums(int _massDV, const int _cDV[], const int _JDV[]);

  // Set and retrieve design variable values
  // ---------------------------------------
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Return the number of displacements and nodes
  // --------------------------------------------
  int getVarsPerNode();
  int getNumNodes();
  ElementLayout getLayoutType();

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitConditions(int elemIndex, const TacsScalar X[], TacsScalar vars[],
                         TacsScalar dvars[], TacsScalar ddvars[]);

  // Retrieve the position of the rigid body
  // ---------------------------------------
  TACSGibbsVector *getInitPosition();

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                       const TacsScalar vars[], const TacsScalar dvars[],
                       TacsScalar *_Te, TacsScalar *_Pe);

  // Compute the residual and Jacobian of the governing equations
  // ------------------------------------------------------------
  void addResidual(int elemIndex, double time, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res);
  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res, TacsScalar *mat);

  // Test the residual implementation
  // --------------------------------
  void testResidual(double dh);
  void testJacobian(double dh, TacsScalar alpha, TacsScalar beta,
                    TacsScalar gamma);

  // Functions for post-processing
  // -----------------------------
  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

  // Rigid body visualization
  //-------------------------
  void setVisualization(TACSRigidBodyViz *viz);

 private:
  // Recompute the inertial properties in the global ref. frame
  void updateInertialProperties();

  // The inertial properties in the global reference frame
  TacsScalar mass;  // The mass of the rigid body
  TacsScalar c[3];  // The first moment of inertia
  TacsScalar J[6];  // The second moment of inertia

  // The initial position, velocity, angular velocity
  TACSGibbsVector *rInit, *vInit, *omegaInit;

  // The accel. due to gravity in the global ref. frame
  TACSGibbsVector *gvec;

  // The initial reference frame
  TACSRefFrame *CRef;

  // The inertial properties in the locally-aligned body frame
  // oriented with respect to CRef
  TacsScalar cRef[3];
  TacsScalar JRef[6];

  // The design variable numbers associated with the inertial properties
  int massDV, cDV[3], JDV[6];

  // Visualization object
  TACSRigidBodyViz *viz;

  // The name of the element
  static const char *elem_name;
};

#endif  // TACS_RIGID_BODY_DYNAMICS_H
