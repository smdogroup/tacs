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


#ifndef TACS_LAM_PARAM_ALL_SHELL_CONSTITUTIVE_H
#define TACS_LAM_PARAM_ALL_SHELL_CONSTITUTIVE_H

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

/*!
  This class implements a lamination parameter based parametrization of
  the shell stiffness and strength properties.

  There are six lamination parameters that define a symmetric balanced
  laminate. These must be combined with appropriate feasibility domain
  for the constraint on the values of the lamination parameters.

  The failure calculations are based on the maximum failure criteria
  taken from a set of set ply angles.

  This implementation is based on the original lpFSDTStiffness class
*/

class TACSLamParamAllShellConstitutive : public TACSShellConstitutive {
 public:
  static const int MAX_NUM_FAIL_ANGLES = 12;

  TACSLamParamAllShellConstitutive(TACSOrthotropicPly* _orthoPly,
                  TacsScalar _t,
                  int _tNum,
                  TacsScalar _tlb,
                  TacsScalar _tub,
                  int _lpNums[],
                  TacsScalar _ksWeight);
  ~TACSLamParamAllShellConstitutive();

  // -------------------------------------
  // Set functions for non-default values
  // --------------------------------------------

  // Set the lamination parameter values directly
  void setLaminationParameters(TacsScalar _lp[]);

  // Set the number of fail angles to test
  // -------------------------------------
  void setNumFailAngles(int _numFailAngles);

  // -------------------------------------
  // Functions for design variable control
  // -------------------------------------

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[], TacsScalar ub[]);



  // -------------------------------------
  // Evaluate mass properties
  // -------------------------------------
  // Evaluate the mass per unit area
  TacsScalar evalDensity(int elemIndex, const double pt[], const TacsScalar X[]);

  // Add the derivative of the density w.r.t. the design variables
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the mass moments
  void evalMassMoments(int elemIndex, const double pt[], const TacsScalar X[],
                       TacsScalar moments[]);

  // Add the sensitivity of the mass moments
  void addMassMomentsDVSens(int elemIndex, const double pt[],
                            const TacsScalar X[], const TacsScalar scale[],
                            int dvLen, TacsScalar dfdx[]);


  // -------------------------------------
  // Evaluate thermal properties
  // -------------------------------------

  // Evaluate the specific heat. Not implemented for this class.
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);


  // -------------------------------------
  // Compute stress/strain/stiffness
  // -------------------------------------
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar e[], TacsScalar s[]);

  // Add the derivative of the product of stress with a vector psi to dfdx
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);


  // -------------------------------------
  // Compute failure criteria
  // -------------------------------------
  // Calculate the point-wise failure criteria
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Evaluate the derivative of the failure criteria w.r.t. the strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

  // Add the derivative of the failure criteria w.r.t. the design variables
  void addFailureDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], const TacsScalar strain[],
                        int dvLen, TacsScalar dfdx[]);

  // -------------------------------------
  // Compute output quantities
  // -------------------------------------
  // Get the object name
  const char* getObjectName() { return constName; }

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

 private:

  // Calculate the failure properties
  void computeFailure(const TacsScalar strain[], TacsScalar fvals[],
                      TacsScalar *_max);
  void computeFailureStrainSens(const TacsScalar strain[],
                                const TacsScalar weights[], TacsScalar sens[]);
  TacsScalar computeFailureDVSens(const TacsScalar strain[], const TacsScalar weights[]);

  // Check that the matrix is positive definite (used for testing)
  int checkDeterminant(const TacsScalar a[]);

  // Get the stiffness matrices based on the current parameter values
  void getStiffness(TacsScalar A[], TacsScalar B[], TacsScalar D[],
                    TacsScalar As[], TacsScalar *drill);

  // The number of angles to check for failure < MAX_NUM_FAIL_ANGLES
  int numFailAngles;

  // The number of design variables
  int numDesignVars;

  // The shear correction factor
  TacsScalar kcorr;

  // The values of the material invariants
  TacsScalar U1, U2, U3, U4, U5;
  TacsScalar U6, U7;  // The invariants for shear

  TACSOrthotropicPly* orthoPly;
  TacsScalar ksWeight;

  // The thickness information
  TacsScalar t;         // The thickness of the laminate
  int tNum;             // The design variable number
  TacsScalar tlb, tub;  // The lower and upper bounds

  // The lamination parameter values
  TacsScalar lp[6];  // The lamination parameters
  int lpNums[6];     // The design variable numbers for the lamination parameters

  static const char* constName;
};

#endif // TACS_LAM_PARAM_ALL_SHELL_CONSTITUTIVE_H
