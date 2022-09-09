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

#ifndef TACS_COMPOSITE_SHELL_CONSTITUTIVE_H
#define TACS_COMPOSITE_SHELL_CONSTITUTIVE_H

#include "TACSMaterialProperties.h"
#include "TACSShellConstitutive.h"

class TACSCompositeShellConstitutive : public TACSShellConstitutive {
 public:
  TACSCompositeShellConstitutive(int _num_plies,
                                 TACSOrthotropicPly **_ply_props,
                                 const TacsScalar *_ply_thickness,
                                 const TacsScalar *_ply_angles,
                                 TacsScalar _kcorr = 5.0 / 6.0);
  ~TACSCompositeShellConstitutive();

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Evaluate the mass moments
  void evalMassMoments(int elemIndex, const double pt[], const TacsScalar X[],
                       TacsScalar moments[]);

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]);

  // Evaluate the stresss
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  // Evaluate failure
  TacsScalar evalFailure(int elemIndex, const double pt[], const TacsScalar X[],
                         const TacsScalar e[]);

  // Evaluate the derivative of the failure criteria w.r.t. the strain
  TacsScalar evalFailureStrainSens(int elemIndex, const double pt[],
                                   const TacsScalar X[], const TacsScalar e[],
                                   TacsScalar sens[]);

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

 private:
  // Store information about the design variable
  int num_plies;
  TACSOrthotropicPly **ply_props;
  TacsScalar *ply_thickness, *ply_angles;
  TacsScalar kcorr;

  // The object name
  static const char *constName;

  void getLaminaStrain(TacsScalar strain[], const TacsScalar rmStrain[],
                       TacsScalar tp);
};

#endif  // TACS_COMPOSITE_SHELL_CONSTITUTIVE_H
