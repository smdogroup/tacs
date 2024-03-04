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

#ifndef TACS_SMEARED_COMPOSITE_SHELL_CONSTITUTIVE_H
#define TACS_SMEARED_COMPOSITE_SHELL_CONSTITUTIVE_H

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

class TACSSmearedCompositeShellConstitutive : public TACSShellConstitutive {
 public:
  TACSSmearedCompositeShellConstitutive(
      int _num_plies, TACSOrthotropicPly **_ply_props, TacsScalar _thickness,
      const TacsScalar *_ply_angles, const TacsScalar *_ply_fractions,
      int _thickness_dv_num = -1, const int *_ply_fraction_dv_nums = NULL,
      TacsScalar _thickness_lb = 0.0, TacsScalar _thickness_ub = 1e20,
      const TacsScalar *_ply_fraction_lb = NULL,
      const TacsScalar *_ply_fraction_ub = NULL, TacsScalar _t_offset = 0.0);
  ~TACSSmearedCompositeShellConstitutive();

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Retrieve the design variable for plotting purposes
  TacsScalar evalDesignFieldValue(int elemIndex, const double pt[],
                                  const TacsScalar X[], int index);

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Add the derivative of the density
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
                  const TacsScalar strain[], TacsScalar stress[]);

  // Add the contribution of the product of stress with a vector psi to dfdx
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate failure
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

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // Evaluate the thermal strain
  void evalThermalStrain(int elemIndex, const double pt[], const TacsScalar X[],
                         TacsScalar theta, TacsScalar strain[]);

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux(int elemIndex, const double pt[], const TacsScalar X[],
                    const TacsScalar grad[], TacsScalar flux[]);

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux(int elemIndex, const double pt[],
                           const TacsScalar X[], TacsScalar C[]);

  // The name of the constitutive object
  const char *getObjectName();

  // Get ply angles, thicknesses, and fractions
  TacsScalar getLaminateThickness();
  void getPlyAngles(TacsScalar *_ply_angles);
  void getPlyFractions(TacsScalar *_ply_fractions);
  TacsScalar getThicknessOffset();

 private:
  // Store information about the design variable
  int num_plies;
  TACSOrthotropicPly **ply_props;
  TacsScalar thickness;
  TacsScalar *ply_fractions, *ply_angles;
  int thickness_dv_num;
  int *ply_fraction_dv_nums;
  TacsScalar thickness_lb, thickness_ub;
  TacsScalar *ply_fraction_lb, *ply_fraction_ub;
  TacsScalar kcorr;
  double ks_weight;
  int nfvals;
  TacsScalar *fvals, *dks_vals;
  TacsScalar **dfvals;
  TacsScalar t_offset;

  // The object name
  static const char *constName;

  TacsScalar evalFSDTStiffness(int elemIndex, const double pt[],
                               const TacsScalar X[], TacsScalar A[],
                               TacsScalar B[], TacsScalar D[], TacsScalar As[]);
  void getLaminaStrain(const TacsScalar rmStrain[], TacsScalar tp,
                       TacsScalar strain[]);
  void evalPlyTopBottomFailure(const TacsScalar strain[], TacsScalar fvals[]);
};

#endif  // TACS_SMEARED_COMPOSITE_SHELL_CONSTITUTIVE_H
