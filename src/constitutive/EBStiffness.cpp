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

#include "EBStiffness.h"

/*
  Set the stiffness and mass properties of the EBStiffness class.
*/
EBStiffness::EBStiffness(TacsScalar rho, TacsScalar E, TacsScalar G,
                         TacsScalar A, TacsScalar Ix, TacsScalar Iy,
                         TacsScalar J, TacsScalar _ref_dir[3],
                         EBReferenceDirection _ref_type)
    : ref_type(_ref_type) {
  // Save the input parameters
  ref_dir[0] = _ref_dir[0];
  ref_dir[1] = _ref_dir[1];
  ref_dir[2] = _ref_dir[2];

  // Set the values of C
  memset(C, 0, sizeof(C));
  C[0] = E * A;
  C[4] = E * Ix;
  C[7] = E * Iy;
  C[9] = G * J;

  // Set the values of the mass moments
  memset(mass, 0, sizeof(mass));
  mass[0] = rho * A;
  mass[3] = rho * Ix;
  mass[5] = rho * Iy;
}

/*
  Return the constitutive name
*/
const char *EBStiffness::constName = "EBStiffness";

/*
  Get the constitutive name
*/
const char *EBStiffness::constitutiveName() { return constName; }

/*
  Get the stiffness matrix
*/
void EBStiffness::getStiffness(const double pt[], TacsScalar Ct[]) {
  memcpy(Ct, C, sizeof(C));
}

/*
  Get the number of stresses
*/
int EBStiffness::getNumStresses() { return NUM_STRESSES; }

/*
  Compute the stress
*/
void EBStiffness::calculateStress(const double pt[], const TacsScalar strain[],
                                  TacsScalar stress[]) {
  calcStress(C, strain, stress);
}

/*
  Add the derivative of the inner product of the stiffness matrix to
  the design variable vector

  Note that this constitutive object does not define any design
  variables.
*/
void addStressDVSens(const double pt[], const TacsScalar strain[],
                     TacsScalar alpha, const TacsScalar psi[],
                     TacsScalar fdvSens, int numDVs) {}

/*
  Return the number of mass moments
*/
int EBStiffness::getNumMassMoments() { return 6; }

/*
  Retrieve the mass moments at the given parametric point
*/
void EBStiffness::getPointwiseMass(const double pt[], TacsScalar _mass[]) {
  memcpy(_mass, mass, 6 * sizeof(TacsScalar));
}

/*
  Add the derivative of the product of the mass components with the
  given input vector.

  Note this constitutive class does not define any design variables.
*/
void EBStiffness::addPointwiseMassDVSens(const double pt[],
                                         const TacsScalar alpha[],
                                         TacsScalar dvSens[], int dvLen) {}
