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

#ifndef TACS_MITC3_H
#define TACS_MITC3_H

#include "TACSBeamBasis.h"
#include "TACSBeamConstitutive.h"
#include "TACSElement.h"
#include "TACSGibbsVector.h"

/*
  A Mixed-Interpolation of Tensorial Components element for dynamic
  analysis.

  The following class implements a geometrically nonlinear
  finite-element shell for large displacement/rotation problems.  The
  element permits arbitrary rotation/displacement rigid body motion.
  The rotational parametrization is based on the quaternions with an
  added constraint at each node.

  The shell formulation utilizes through-thickness strain/kinematic
  assumptions that are based on first-order shear deformation
  theory. The theory takes into account the nonlinear rotational
  kinematics that are required to obtain strain-free rotation of the
  elements.

  The drilling degree of freedom is handled through the use of a penalty
  term that penalizes the discrepancy between the in-plane rotations
  predicted from nonlinear shell theory and those predicted by stress-
*/
class MITC3 : public TACSElement {
 public:
  static const int ORDER = 3;
  static const int NUM_NODES = ORDER;
  static const int NUM_DISPS = 8;
  static const int NUM_STRESSES = 6;

  MITC3(TACSBeamConstitutive *_stiff, TACSGibbsVector *_gravity = NULL,
        TACSGibbsVector *_vInit = NULL, TACSGibbsVector *_omegaInit = NULL);
  ~MITC3();

  // Return the sizes of the array components
  // ----------------------------------------
  int getVarsPerNode();
  int getNumNodes();

  // Functions to determine the variable names and quantities
  // --------------------------------------------------------
  const char *getObjectName();
  ElementLayout getLayoutType();

  // Get the element basis
  // ---------------------
  TACSElementBasis *getElementBasis() { return &basis; }

  // Functions for handling the design variables
  // -------------------------------------------
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Retrieve the initial values of the state variables
  // --------------------------------------------------
  void getInitConditions(int elemIndex, const TacsScalar X[], TacsScalar vars[],
                         TacsScalar dvars[], TacsScalar ddvars[]);

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                       const TacsScalar vars[], const TacsScalar dvars[],
                       TacsScalar *Te, TacsScalar *Pe);

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[]);

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]);

  // Add the adjoint-residual product
  // --------------------------------
  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dfdx[]);

  // Functions for post-processing
  // -----------------------------
  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

  // Evaluate a point-wise quantity of interest.
  // ------------------------------------------
  int evalPointQuantity(int elemIndex, int quantityType, double time, int n,
                        double pt[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar *quantity);

  // Add the derivative of the point quantity w.r.t. the design variables
  // --------------------------------------------------------------------
  void addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                              TacsScalar scale, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], int dvLen,
                              TacsScalar dfdx[]);

  // Add the derivative of the point quantity w.r.t. the state variables
  // -------------------------------------------------------------------
  void addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                              TacsScalar alpha, TacsScalar beta,
                              TacsScalar gamma, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], TacsScalar dfdu[]);

  // Add the derivative of the point quantity w.r.t. the node locations
  // ------------------------------------------------------------------
  void addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfdq[], TacsScalar dfdXpts[]);

  // Member functions for evaluating global functions of interest
  // ------------------------------------------------------------
  TACSConstitutive *getConstitutive();

  // Get the number of Gauss quadrature points
  // -----------------------------------------
  int getNumGaussPts();
  double getGaussWtsPts(const int num, double pt[]);

  // Get the shape functions from the element
  // ----------------------------------------
  void getShapeFunctions(const double pt[], double N[]);

  // Return the determinant of the Jacobian of the transformation
  // ------------------------------------------------------------
  TacsScalar getDetJacobian(const double pt[], const TacsScalar X[]);

  // Get the strain and the parametric location from the element
  // -----------------------------------------------------------
  void getStrain(const double pt[], const TacsScalar X[],
                 const TacsScalar vars[], TacsScalar e[]);

  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  void addStrainSVSens(const double pt[], const TacsScalar scale,
                       const TacsScalar esens[], const TacsScalar Xpts[],
                       const TacsScalar vars[], TacsScalar sens[]);

  // Test the strain expressions
  void testStrain(const TacsScalar X[]);

 private:
  // Helper functions required for analysis
  void computeAngularVelocity(TacsScalar omega[], const TacsScalar vars[],
                              const TacsScalar dvars[]);

  // Compute the angular acceleration at the nodes
  void computeAngularAccel(TacsScalar domega[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[]);

  // Compute the local frame for strain computations
  TacsScalar computeTransform(TacsScalar T[], const TacsScalar Xa[]);

  // Compute the reference frames at each node of the element
  void computeFrames(TacsScalar Xr[], const TacsScalar X[]);

  // Compute the derivative of the transformation
  void computeTransformSens(TacsScalar Xad[], TacsScalar Xbd[],
                            const TacsScalar Td[], const TacsScalar Xa[],
                            const TacsScalar Xb[]);

  // Compute the directors for each node
  void computeDirectors(TacsScalar d1[], TacsScalar d2[],
                        const TacsScalar vars[], const TacsScalar Xr[]);
  void addDirectorsSens(TacsScalar Xrd[], const TacsScalar d1d[],
                        const TacsScalar d2d[], const TacsScalar vars[]);

  // Compute the derivative of the directors w.r.t. the variables
  void computeDirectorDeriv(TacsScalar d1dq[], TacsScalar d2dq[],
                            const TacsScalar vars[], const TacsScalar Xr[]);

  // Evaluate the strain
  void evalStrain(TacsScalar e[], const TacsScalar Ur[], const TacsScalar d1a[],
                  const TacsScalar d2a[], const TacsScalar Xdinv[],
                  const TacsScalar z1Xdinv[], const TacsScalar z2Xdinv[],
                  const TacsScalar T[]);

  // Evaluate the bmat matrix
  void evalBmat(TacsScalar e[], TacsScalar B[], const double N[],
                const double Na[], const TacsScalar Ur[],
                const TacsScalar d1a[], const TacsScalar d2a[],
                const TacsScalar Xdinv[], const TacsScalar z1Xdinv[],
                const TacsScalar z2Xdinv[], const TacsScalar T[],
                const TacsScalar d1dq[], const TacsScalar d2dq[]);

  // Add the terms from the geometric stiffness matrix
  void addGmat(TacsScalar J[], const TacsScalar scale, const TacsScalar s[],
               const double N[], const double Na[], const TacsScalar Ur[],
               const TacsScalar d1a[], const TacsScalar d2a[],
               const TacsScalar Xdinv[], const TacsScalar z1Xdinv[],
               const TacsScalar z2Xdinv[], const TacsScalar T[],
               const TacsScalar Xr[], const TacsScalar d1dq[],
               const TacsScalar d2dq[]);

  // Compute the tying strain
  void computeTyingStrain(TacsScalar g12[], TacsScalar g13[],
                          const TacsScalar X[], const TacsScalar Xr[],
                          const TacsScalar vars[], const TacsScalar d1[],
                          const TacsScalar d2[]);

  // Compute the tying bmat matrices
  void computeTyingBmat(TacsScalar g12[], TacsScalar g13[], TacsScalar B12[],
                        TacsScalar B13[], const TacsScalar X[],
                        const TacsScalar Xr[], const TacsScalar vars[],
                        const TacsScalar d1[], const TacsScalar d2[],
                        const TacsScalar dir1dq[], const TacsScalar dir2dq[]);

  // Add the tying strain
  void addTyingStrain(TacsScalar e[], const double N12[],
                      const TacsScalar g12[], const TacsScalar g13[]);

  // Add the derivative of the tying strain
  void addTyingBmat(TacsScalar B[], const double N12[], const TacsScalar b12[],
                    const TacsScalar b13[]);

  // Add the second derivative of the strain at the tying points
  void addTyingGmat(TacsScalar J[], const TacsScalar w12[],
                    const TacsScalar w13[], const TacsScalar X[],
                    const TacsScalar Xr[], const TacsScalar vars[],
                    const TacsScalar d1[], const TacsScalar d2[],
                    const TacsScalar d1dq[], const TacsScalar d2dq[]);

  // Compute the inertia tensor
  void computeInertiaTensor(const TacsScalar rho[], const TacsScalar n1[],
                            const TacsScalar n2[], TacsScalar Jr[]);

  // Compute the product of the stress and the strain
  inline TacsScalar strainProduct(const TacsScalar s[], const TacsScalar e[]) {
    return (e[0] * s[0] + e[1] * s[1] + e[2] * s[2] + e[3] * s[3] +
            e[4] * s[4] + e[5] * s[5]);
  }

  // Set the pointers to quadrature points/weights
  const double *gaussPts, *gaussWts;

  // The stiffness object
  TACSBeamConstitutive *stiff;

  // The gravity vector (if any)
  TACSGibbsVector *gravity;

  // Initial velocity/angular velocity
  TACSGibbsVector *vInit, *omegaInit;

  // The element name
  static const char *elemName;

  // The element basis
  TACSQuadraticBeamBasis basis;
};

#endif  // TACS_MITC3_H
