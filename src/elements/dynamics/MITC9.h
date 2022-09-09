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

#ifndef TACS_MITC9_H
#define TACS_MITC9_H

#include "FSDTStiffness.h"
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
class MITC9 : public TACSElement {
 public:
  static const int ORDER = 3;
  static const int NUM_NODES = ORDER * ORDER;
  static const int NUM_DISPS = 8;
  static const int NUM_STRESSES = 8;

  MITC9(FSDTStiffness *_stiff, TACSGibbsVector *_gravity = NULL,
        TACSGibbsVector *_vInit = NULL, TACSGibbsVector *_omegaInit = NULL);
  ~MITC9();

  // Return the sizes of the array components
  // ----------------------------------------
  int getVarsPerNode();
  int getNumNodes();

  // Functions to determine the variable names and quantities
  // --------------------------------------------------------
  const char *getObjectName();
  ElementLayout getLayoutType();

  // Functions for handling the design variables
  // -------------------------------------------
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);
  void setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);
  void getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);
  void getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
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

  // Derivatives for the adjoint equations
  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dvSens[]);
  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

  // Functions for post-processing
  // -----------------------------
  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

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
  TacsScalar getDetJacobianXptSens(TacsScalar hXptSens[], const double pt[],
                                   const TacsScalar X[]);

  // Get the strain and the parametric location from the element
  // -----------------------------------------------------------
  void getStrain(TacsScalar e[], const double pt[], const TacsScalar X[],
                 const TacsScalar vars[]);

  // This function adds the sensitivity of the strain to the state variables
  // -----------------------------------------------------------------------
  void addStrainSVSens(TacsScalar sens[], const double pt[],
                       const TacsScalar scale, const TacsScalar esens[],
                       const TacsScalar Xpts[], const TacsScalar vars[]);

  // This function adds the sensitivity of the strain w.r.t. Xpts
  // ------------------------------------------------------------
  void addStrainXptSens(TacsScalar strainXptSens[], const double pt[],
                        const TacsScalar scale, const TacsScalar strainSens[],
                        const TacsScalar Xpts[], const TacsScalar vars[]);

  // Functions for post-processing
  // -----------------------------
  void addOutputCount(int *nelems, int *nnodes, int *ncsr);
  void getOutputData(unsigned int out_type, double *data, int ld_data,
                     const TacsScalar Xpts[], const TacsScalar vars[]);
  void getOutputConnectivity(int *con, int node);

  // Test the strain implementation
  // ------------------------------
  void testStrain(const TacsScalar X[]);

  // Test member helper functions for geometric derivatives
  // ------------------------------------------------------
  void testXptSens(double dh = 1e-6);

 private:
  // Helper functions required for analysis
  void computeAngularVelocity(TacsScalar omega[], const TacsScalar vars[],
                              const TacsScalar dvars[]);

  // Compute the angular acceleration at the nodes
  void computeAngularAccel(TacsScalar domega[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[]);

  // Compute the reference frames at each node of the element
  void computeFrames(TacsScalar Xr[], const TacsScalar X[]);
  void addFramesSens(TacsScalar Xd[], const TacsScalar Xrd[],
                     const TacsScalar X[]);

  // Compute the local frame for strain computations
  void computeTransform(TacsScalar T[], const TacsScalar Xa[],
                        const TacsScalar Xb[]);

  // Compute the derivative of the transformation
  void computeTransformSens(TacsScalar Xad[], TacsScalar Xbd[],
                            const TacsScalar Td[], const TacsScalar Xa[],
                            const TacsScalar Xb[]);

  // Compute the directors for each node
  void computeDirectors(TacsScalar d[], const TacsScalar vars[],
                        const TacsScalar Xr[]);
  void addDirectorsSens(TacsScalar Xrd[], const TacsScalar dd[],
                        const TacsScalar vars[]);

  // Compute the derivative of the directors w.r.t. the variables
  void computeDirectorDeriv(TacsScalar dirdq[], const TacsScalar vars[],
                            const TacsScalar Xr[]);
  void addDirectorDerivSens(TacsScalar Xrd[], const TacsScalar ddqd[],
                            const TacsScalar vars[]);

  // Evaluate the strain
  void evalStrain(TacsScalar e[], const TacsScalar Ur[], const TacsScalar dr[],
                  const TacsScalar Xdinv[], const TacsScalar zXdinv[],
                  const TacsScalar T[]);

  // Evaluate the derivative of the strain w.r.t. inputs
  void evalStrainSens(TacsScalar Urd[], TacsScalar drd[], TacsScalar Xdinvd[],
                      TacsScalar zXdinvd[], TacsScalar Td[], TacsScalar scale,
                      const TacsScalar eSens[], const TacsScalar Ur[],
                      const TacsScalar dr[], const TacsScalar Xdinv[],
                      const TacsScalar zXdinv[], const TacsScalar T[]);

  // Evaluate the derivative of the strain w.r.t. the element variables
  void evalBmat(TacsScalar e[], TacsScalar B[], const double N[],
                const double Na[], const double Nb[], const TacsScalar Ur[],
                const TacsScalar dr[], const TacsScalar Xdinv[],
                const TacsScalar zXdinv[], const TacsScalar T[],
                const TacsScalar dirdq[]);

  // Evaluate the derivative of the eSens^{T}*Bmat*psi w.r.t. node locations
  void addBmatSens(TacsScalar Urd[], TacsScalar drd[], TacsScalar Xdinvd[],
                   TacsScalar zXdinvd[], TacsScalar Td[], TacsScalar dirdqd[],
                   const TacsScalar eSens[], const TacsScalar psi[],
                   const double N[], const double Na[], const double Nb[],
                   const TacsScalar Ur[], const TacsScalar dr[],
                   const TacsScalar Xdinv[], const TacsScalar zXdinv[],
                   const TacsScalar T[], const TacsScalar dirdq[]);

  // Add the interpolated strain to the strain vector
  void addTyingStrain(TacsScalar e[], const double N13[], const double N23[],
                      const TacsScalar g13[], const TacsScalar g23[],
                      const TacsScalar Xdinv[], const TacsScalar T[]);

  // Add the derivative of the strain w.r.t. the outputs
  void addTyingStrainSens(TacsScalar g13d[], TacsScalar g23d[],
                          TacsScalar Xdinvd[], TacsScalar Td[],
                          TacsScalar scale, const TacsScalar e[],
                          const double N13[], const double N23[],
                          const TacsScalar g13[], const TacsScalar g23[],
                          const TacsScalar Xdinv[], const TacsScalar T[]);

  // Add the contribution from the tying strain to the b-matrix
  void addTyingBmat(TacsScalar B[], const double N13[], const double N23[],
                    const TacsScalar b13[], const TacsScalar b23[],
                    const TacsScalar Xdinv[], const TacsScalar T[]);
  void addTyingBmatSens(TacsScalar B13d[], TacsScalar B23d[],
                        TacsScalar Xdinvd[], TacsScalar Td[],
                        const TacsScalar eSens[], const TacsScalar psi[],
                        const double N13[], const double N23[],
                        const TacsScalar B13[], const TacsScalar B23[],
                        const TacsScalar Xdinv[], const TacsScalar T[]);

  // Compute the shear strain at the tying points
  void computeTyingStrain(TacsScalar g13[], TacsScalar g23[],
                          const TacsScalar X[], const TacsScalar Xr[],
                          const TacsScalar vars[], const TacsScalar dir[]);
  void addComputeTyingStrainSens(TacsScalar Xd[], TacsScalar Xrd[],
                                 TacsScalar dird[], const TacsScalar g13d[],
                                 const TacsScalar g23d[], const TacsScalar X[],
                                 const TacsScalar Xr[], const TacsScalar vars[],
                                 const TacsScalar dir[]);

  // Compute the derivative of the strain at the tying points
  void computeTyingBmat(TacsScalar g13[], TacsScalar g23[], TacsScalar b13[],
                        TacsScalar b23[], const TacsScalar X[],
                        const TacsScalar Xr[], const TacsScalar vars[],
                        const TacsScalar dir[], const TacsScalar dirdq[]);
  void addComputeTyingBmatSens(TacsScalar Xd[], TacsScalar Xrd[],
                               TacsScalar dird[], TacsScalar dirdqd[],
                               const TacsScalar g13d[], const TacsScalar g23d[],
                               const TacsScalar B13d[], const TacsScalar B23d[],
                               const TacsScalar X[], const TacsScalar Xr[],
                               const TacsScalar vars[], const TacsScalar dir[],
                               const TacsScalar dirdq[]);

  // Add the terms from the geometric stiffness matrix
  void addGmat(TacsScalar J[], const TacsScalar scale, const TacsScalar s[],
               const double N[], const double Na[], const double Nb[],
               const TacsScalar Ur[], const TacsScalar dr[],
               const TacsScalar Xdinv[], const TacsScalar zXdinv[],
               const TacsScalar T[], const TacsScalar Xr[],
               const TacsScalar dirdq[]);

  // Add the term from the drilling rotation
  TacsScalar computeRotPenalty(const double N[], const TacsScalar Xa[],
                               const TacsScalar Xb[], const TacsScalar Ua[],
                               const TacsScalar Ub[], const TacsScalar vars[]);

  // Compute the derivative of the rotation term
  TacsScalar computeBRotPenalty(TacsScalar brot[], const double N[],
                                const double Na[], const double Nb[],
                                const TacsScalar Xa[], const TacsScalar Xb[],
                                const TacsScalar Ua[], const TacsScalar Ub[],
                                const TacsScalar vars[]);

  // Add the derivative of the penalty term
  void addBRotPenaltySens(TacsScalar Xad[], TacsScalar Xbd[],
                          const TacsScalar rotd, const TacsScalar scale,
                          const TacsScalar psi[], const double N[],
                          const double Na[], const double Nb[],
                          const TacsScalar Xa[], const TacsScalar Xb[],
                          const TacsScalar Ua[], const TacsScalar Ub[],
                          const TacsScalar vars[]);

  // Add the geometric stiffness term from the rotation
  void addGRotMat(TacsScalar J[], const TacsScalar scale, const double N[],
                  const double Na[], const double Nb[], const TacsScalar Xa[],
                  const TacsScalar Xb[], const TacsScalar Ua[],
                  const TacsScalar Ub[], const TacsScalar vars[]);

  // Add the geometric stiffness matrix from the tying strain
  void addTyingGmat(TacsScalar J[], const TacsScalar w13[],
                    const TacsScalar w23[], const TacsScalar X[],
                    const TacsScalar Xr[], const TacsScalar vars[],
                    const TacsScalar dir[], const TacsScalar dirdq[]);

  // Add to the weights required to compute the geometric stiffness
  void addTyingGmatWeights(TacsScalar w13[], TacsScalar w23[],
                           const TacsScalar scalar, const TacsScalar s[],
                           const double N13[], const double N23[],
                           const TacsScalar Xdinv[], const TacsScalar T[]);

  // Perform tests of the underlying implementation for the Xpt derivatives
  void testInv3x3Sens(double dh);
  void testStrainSens(double dh);
  void testTransformSens(double dh);
  void testNormalRateSens(double dh);
  void testTyingStrainSens(double dh);
  void testBmatSens(double dh);
  void testTyingBmatSens(double dh);
  void testBrotSens(double dh);
  void testFrameSens(double dh);

  // Compute the product of the stress and the strain
  inline TacsScalar strainProduct(const TacsScalar s[], const TacsScalar e[]) {
    return (e[0] * s[0] + e[1] * s[1] + e[2] * s[2] + e[3] * s[3] +
            e[4] * s[4] + e[5] * s[5] + e[6] * s[6] + e[7] * s[7]);
  }

  // Set the pointers to quadrature points/weights
  const double *gaussPts, *gaussWts;

  // Set the drill rotation omega factor
  double drill_inertia_factor;

  // The stiffness object
  FSDTStiffness *stiff;

  // The gravity vector (if any)
  TACSGibbsVector *gravity;

  // Initial velocity/angular velocity
  TACSGibbsVector *vInit, *omegaInit;

  // The names of the displacements, stresses etc.
  static const char *elemName;
  static const char *dispNames[NUM_DISPS];
  static const char *stressNames[NUM_STRESSES];
  static const char *strainNames[NUM_STRESSES];
  static const char *extraNames[NUM_EXTRAS];
};

#endif  // TACS_MITC9_H
