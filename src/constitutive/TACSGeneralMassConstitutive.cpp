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

#include "TACSGeneralMassConstitutive.h"

const char* TACSGeneralMassConstitutive::name = "TACSGeneralMassConstitutive";

const char* TACSGeneralMassConstitutive::getObjectName() { return name; }

/*
  GeneralMassConstitutive member function definitions
*/
TACSGeneralMassConstitutive::TACSGeneralMassConstitutive(
    const TacsScalar _M[]) {
  memcpy(M, _M, 21 * sizeof(TacsScalar));
}

/*
  Default constructor
*/
TACSGeneralMassConstitutive::TACSGeneralMassConstitutive() {
  memset(M, 0, 21 * sizeof(TacsScalar));
}

int TACSGeneralMassConstitutive::getNumStresses() { return NUM_STRESSES; }

/*
  Given the mass matrix and an acceleration vector, evaluate the inertial forces
*/
void TACSGeneralMassConstitutive::evalInertia(int elemIndex, const double pt[],
                                              const TacsScalar X[],
                                              const TacsScalar ddu[],
                                              TacsScalar f[]) {
  TacsScalar M[21];
  evalMassMatrix(elemIndex, pt, X, M);

  f[0] = M[0] * ddu[0] + M[1] * ddu[1] + M[2] * ddu[2] + M[3] * ddu[3] +
         M[4] * ddu[4] + M[5] * ddu[5];
  f[1] = M[1] * ddu[0] + M[6] * ddu[1] + M[7] * ddu[2] + M[8] * ddu[3] +
         M[9] * ddu[4] + M[10] * ddu[5];
  f[2] = M[2] * ddu[0] + M[7] * ddu[1] + M[11] * ddu[2] + M[12] * ddu[3] +
         M[13] * ddu[4] + M[14] * ddu[5];
  f[3] = M[3] * ddu[0] + M[8] * ddu[1] + M[12] * ddu[2] + M[15] * ddu[3] +
         M[16] * ddu[4] + M[17] * ddu[5];
  f[4] = M[4] * ddu[0] + M[9] * ddu[1] + M[13] * ddu[2] + M[16] * ddu[3] +
         M[18] * ddu[4] + M[19] * ddu[5];
  f[5] = M[5] * ddu[0] + M[10] * ddu[1] + M[14] * ddu[2] + M[17] * ddu[3] +
         M[19] * ddu[4] + M[20] * ddu[5];
}

// Evaluate the mass matrix
void TACSGeneralMassConstitutive::evalMassMatrix(int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 TacsScalar C[]) {
  memcpy(C, M, 21 * sizeof(TacsScalar));
}

// Evaluate the material density
TacsScalar TACSGeneralMassConstitutive::evalDensity(int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[]) {
  TacsScalar mass = (M[0] + M[6] + M[11]) / 3.0;
  return mass;
}
