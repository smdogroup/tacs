#include "TACSGeneralSpringConstitutive.h"

/*
  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/

const char* TACSGeneralSpringConstitutive::constName =
    "TACSGeneralSpringConstitutive";
const char* TACSGeneralSpringConstitutive::constitutiveName() {
  return constName;
}

TACSGeneralSpringConstitutive::TACSGeneralSpringConstitutive(TacsScalar _C[]) {
  // Copy stiffness matrix entries
  memcpy(C, _C, 21 * sizeof(TacsScalar));
}

TACSGeneralSpringConstitutive::TACSGeneralSpringConstitutive() {
  // Spring has zero stiffness
  memset(C, 0, 21 * sizeof(TacsScalar));
}

int TACSGeneralSpringConstitutive::getNumStresses() { return NUM_STRESSES; }

void TACSGeneralSpringConstitutive::evalStress(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar strain[],
                                               TacsScalar stress[]) {
  TacsScalar Ct[21];
  evalTangentStiffness(elemIndex, pt, X, Ct);
  calcStress(Ct, strain, stress);
}

void TACSGeneralSpringConstitutive::evalTangentStiffness(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[],
                                                         TacsScalar Ct[]) {
  memcpy(Ct, C, 21 * sizeof(TacsScalar));
}