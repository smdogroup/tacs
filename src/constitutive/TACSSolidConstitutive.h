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

#ifndef TACS_SOLID_CONSTITUTIVE_H
#define TACS_SOLID_CONSTITUTIVE_H

#include "TACSConstitutive.h"
#include "TACSMaterialProperties.h"

/*
  This is the base class for the solid constitutive objects.

  All objects performing solid elastic analysis should utilize this class.
*/
class TACSSolidConstitutive : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 6;

  TACSSolidConstitutive(TACSMaterialProperties *properties, TacsScalar _t = 1.0,
                        int _tNum = -1, TacsScalar _tlb = 0.0,
                        TacsScalar _tub = 1.0);
  ~TACSSolidConstitutive();

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

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Add the derivative of the density
  void addDensityDVSens(int elemIndex, TacsScalar scale, const double pt[],
                        const TacsScalar X[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // Add the derivative of the density
  void addSpecificHeatDVSens(int elemIndex, TacsScalar scale, const double pt[],
                             const TacsScalar X[], int dvLen,
                             TacsScalar dfdx[]);

  // Evaluate the stresss
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]);

  // Evaluate the thermal strain
  void evalThermalStrain(int elemIndex, const double pt[], const TacsScalar X[],
                         TacsScalar theta, TacsScalar strain[]);

  // Add the contribution
  void addStressDVSens(int elemIndex, TacsScalar scale, const double pt[],
                       const TacsScalar X[], const TacsScalar strain[],
                       const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux(int elemIndex, const double pt[], const TacsScalar X[],
                    const TacsScalar grad[], TacsScalar flux[]);

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux(int elemIndex, const double pt[],
                           const TacsScalar X[], TacsScalar C[]);

  // Add the derivative of the heat flux
  void addHeatFluxDVSens(int elemIndex, TacsScalar scale, const double pt[],
                         const TacsScalar X[], const TacsScalar grad[],
                         const TacsScalar psi[], int dvLen, TacsScalar dfdx[]);

  // Evaluate the material failure index
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Evaluate the derivative of the failure criteria w.r.t. strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

  // Extra info about the constitutive class
  const char *getObjectName();

  // Return the material properties object
  TACSMaterialProperties *getMaterialProperties();

 protected:
  // Materiial properties class
  TACSMaterialProperties *properties;

 private:
  // Store information about the design variable
  TacsScalar t, tlb, tub;
  int tNum;

  static const char *sName;
};

#endif  // TACS_SOLID_CONSTITUTIVE_H
