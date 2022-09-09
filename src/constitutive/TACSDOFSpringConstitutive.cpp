#include "TACSDOFSpringConstitutive.h"

/*
  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/

const char* TACSDOFSpringConstitutive::constName = "TACSDOFSpringConstitutive";
const char* TACSDOFSpringConstitutive::constitutiveName() { return constName; }

TACSDOFSpringConstitutive::TACSDOFSpringConstitutive(TacsScalar k[]) {
  // Set the diagonal terms based on 6 DOF stiffness
  C[0] = k[0];
  C[6] = k[1];
  C[11] = k[2];
  C[15] = k[3];
  C[18] = k[4];
  C[20] = k[5];
}
