/*
  This file is part of the package TACS.

  Copyright (C) 2015 Georgia Tech Research Corporation.
  Additional copyright (C) 2015 Graeme Kennedy.
  All rights reserved.

  TMR is licensed under the Apache License, Version 2.0 (the "License");
  you may not use this software except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#ifndef TACS_LAM_PARAM_SHELL_CONSTITUTIVE_H
#define TACS_LAM_PARAM_SHELL_CONSTITUTIVE_H

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

/*
  This class implements a lamination parameter based parametrization of
  the shell stiffness and strength properties.

  The class is restricted to symmetric balanced laminates.
  The in-plane properties are parametrized in terms of laminate fractions
  at the angles 0, +/-45 and 90 degrees. The thickness is treated as a
  continuous design variable. Additional lamination parameters W1 and W3
  are used for the bending stiffness.

  The failure index calculations are based on the maximum strain failure
  criteria imposed at 0, pm 45, 90 ply angles.
*/
class TACSLamParamShellConstitutive : public TACSShellConstitutive {
 public:
  TACSLamParamShellConstitutive(TACSOrthotropicPly *_orthoPly, TacsScalar _t,
                                int _t_num, TacsScalar _min_t,
                                TacsScalar _max_t, TacsScalar _f0,
                                TacsScalar _f45, TacsScalar _f90, int _f0_num,
                                int _f45_num, int _f90_num, TacsScalar _min_f0,
                                TacsScalar _min_f45, TacsScalar _min_f90,
                                TacsScalar _W1, TacsScalar _W3, int _W1_num,
                                int _W3_num, TacsScalar _ksWeight,
                                TacsScalar _epsilon);
  ~TACSLamParamShellConstitutive();

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Evaluate the mass per unit area
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

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

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // Evaluate the stress
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar e[], TacsScalar s[]);

  // Add the derivative of the product of stress with a vector psi to dfdx
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

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

  // Get the object name
  const char *getObjectName() { return constName; }

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

 private:
  static const int NUM_FAIL_ANGLES = 4;

  // Calculate the failure properties
  void computeFailure(const TacsScalar strain[], TacsScalar fvals[],
                      TacsScalar *_max);
  void computeFailureStrainSens(const TacsScalar strain[],
                                const TacsScalar weights[], TacsScalar sens[]);

  // Check that the matrix is positive definite (used for testing)
  int checkDeterminant(const TacsScalar a[]);

  // Get the stiffness matrices based on the current parameter values
  void getStiffness(TacsScalar A[], TacsScalar B[], TacsScalar D[],
                    TacsScalar As[], TacsScalar *drill);

  // The number of design variables
  int numDesignVars;

  // The shear correction factor
  TacsScalar kcorr;

  // The values of the material invariants
  TacsScalar U1, U2, U3, U4, U5;
  TacsScalar U6, U7;  // The invariants for shear

  TACSOrthotropicPly *orthoPly;
  TacsScalar ksWeight, epsilon;

  // The thickness information
  TacsScalar t;         // The thickness of the laminate
  int tNum;             // The design variable number
  TacsScalar tlb, tub;  // The lower and upper bounds

  // The ply fraction information
  int n0, n45, n90;                     // The fraction dv numbers
  TacsScalar f0, f45, f90;              // The fraction values
  TacsScalar min_f0, min_f45, min_f90;  // The lower fraction bounds

  // The lamination parameter values
  int nW1, nW3;       // The design variable numbers
  TacsScalar W1, W3;  // The lamination parameter values

  static const char *constName;
};

#endif  // TACS_LAM_PARAM_SHELL_CONSTITUTIVE_H
