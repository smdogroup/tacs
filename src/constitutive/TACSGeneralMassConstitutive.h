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

#ifndef TACS_GENERAL_MASS_CONSTITUTIVE_H
#define TACS_GENERAL_MASS_CONSTITUTIVE_H

#include "TACSConstitutive.h"

/*
  This is the base class for the point mass constitutive objects.
  Assumes 6 dofs.

*/
class TACSGeneralMassConstitutive : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 0;

  TACSGeneralMassConstitutive(const TacsScalar _M[]);

  TACSGeneralMassConstitutive();

  int getNumStresses();

  // Given the mass matrix and an acceleration vector, evaluate the inertial
  // forces
  void evalInertia(int elemIndex, const double pt[], const TacsScalar X[],
                   const TacsScalar ddu[], TacsScalar f[]);

  // Evaluate the mass matrix
  void evalMassMatrix(int elemIndex, const double pt[], const TacsScalar X[],
                      TacsScalar C[]);

  // Evaluate the material density
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]);

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]) {
    return 0.0;
  }

  // Evaluate the stresss
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]) {
    return;
  }

  // Evaluate the tangent stiffness
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[]) {
    return;
  }

  // Extra info about the constitutive class
  const char *getObjectName();

 protected:
  // Mass matrix
  TacsScalar M[21];

 private:
  static const char *name;
};

#endif  // TACS_GENERAL_MASS_CONSTITUTIVE_H
