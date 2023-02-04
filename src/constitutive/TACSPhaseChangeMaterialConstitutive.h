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

#ifndef TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H
#define TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H

#include "TACSConstitutive.h"
#include "TACSMaterialProperties.h"

/*
  This is the base class for the phase change material constitutive objects.
*/
class TACSPhaseChangeMaterialConstitutive : public TACSConstitutive {
 public:
  TACSPhaseChangeMaterialConstitutive(TACSMaterialProperties *solid_properties,
                                      TACSMaterialProperties *liquid_properties,
                                      TacsScalar _lh, TacsScalar _Tm,
                                      TacsScalar _dT = 10.0,
                                      TacsScalar _t = 1.0, int _tNum = -1,
                                      TacsScalar _tlb = 0.0,
                                      TacsScalar _tub = 1.0);
  ~TACSPhaseChangeMaterialConstitutive();

  int getNumStresses();

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Compute the phase change coefficient
  TacsScalar evalTransitionCoef(const TacsScalar T);

  TacsScalar evalTransitionCoefSVSens(const TacsScalar T);

  // Evaluate the material's phase
  int evalPhase(const TacsScalar T);

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar u[]);

  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]) {
    TacsScalar u[1];
    u[0] = 0.0;
    return evalDensity(elemIndex, pt, X, u);
  }

  // Add the derivative of the density
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[],
                        const TacsScalar u[]);

  void addDensitySVSens(int elemIndex, const double pt[], const TacsScalar X[],
                        TacsScalar dfdu[], const TacsScalar u[]);

  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[], const TacsScalar u[]);

  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]) {
    TacsScalar u[1];
    u[0] = 0.0;
    return evalSpecificHeat(elemIndex, pt, X, u);
  }

  void addSpecificHeatSVSens(int elemIndex, const double pt[],
                             const TacsScalar X[], TacsScalar dfdu[],
                             const TacsScalar u[]);

  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux(int elemIndex, const double pt[], const TacsScalar X[],
                    const TacsScalar grad[], TacsScalar flux[],
                    const TacsScalar u[]);

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux(int elemIndex, const double pt[],
                           const TacsScalar X[], TacsScalar Kc[],
                           const TacsScalar u[]);

  // Add the derivative of the heat flux
  void addHeatFluxDVSens(int elemIndex, TacsScalar scale, const double pt[],
                         const TacsScalar X[], const TacsScalar grad[],
                         const TacsScalar psi[], int dvLen, TacsScalar dfdx[],
                         const TacsScalar u[]);

  void addKappaSVSens(int elemIndex, const double pt[], const TacsScalar X[],
                      TacsScalar dfdu[], const TacsScalar u[]);

  // Extra info about the constitutive class
  const char *getObjectName();

 protected:
  // Material properties class
  TACSMaterialProperties *solid_properties;
  TACSMaterialProperties *liquid_properties;

 private:
  // Store information about the design variable
  TacsScalar lh, Tm, t, tlb, tub, dT, b;
  int tNum;

  static const char *psName;
};

#endif  // TACS_PHASE_CHANGE_MATERIAL_CONSTITUTIVE_H
