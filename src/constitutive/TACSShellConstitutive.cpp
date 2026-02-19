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

#include "TACSShellConstitutive.h"

#include "tacslapack.h"

const char *TACSShellConstitutive::constName = "TACSShellConstitutive";

/*
  Return the constitutive name
*/
const char *TACSShellConstitutive::getObjectName() { return constName; }

/*
  Get the number of stresses
*/
int TACSShellConstitutive::getNumStresses() { return NUM_STRESSES; }

// Extract the tangent stiffness components from the matrix
void TACSShellConstitutive::extractTangentStiffness(
    const TacsScalar *C, const TacsScalar **A, const TacsScalar **B,
    const TacsScalar **D, const TacsScalar **As, TacsScalar *drill) {
  if (A) {
    *A = &C[0];
  }
  if (B) {
    *B = &C[6];
  }
  if (D) {
    *D = &C[12];
  }
  if (As) {
    *As = &C[18];
  }
  if (drill) {
    *drill = C[21];
  }
}

/*
  Set the default drilling regularization value
*/
double TACSShellConstitutive::DRILLING_REGULARIZATION = 10.0;

/*
  Set the drilling stiffness regularization parameter
*/
void TACSShellConstitutive::setDrillingRegularization(double kval) {
  DRILLING_REGULARIZATION = kval;
}
