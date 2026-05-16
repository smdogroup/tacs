/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSBasicBeamConstitutive.h"

/*
  Constructor for Timoshenko beam theory based constitutive object.

  EA               : axial stiffness
  EI22, EI33, EI23 : bending stiffness
  GJ               : torsional stiffness
  kG22, kG33, kG23 : shearing stiffness
  m00              : mass per unit span
  m11, m22, m33    : moments of inertia such that m11 = m22 + m33
  xm2, xm3         : cross sectional center of mass location
  xc2, xc3         : cross sectional centroid
  xk2, xk3         : cross sectional shear center
  muS              : viscous damping coefficient
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TacsScalar EA, TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
    TacsScalar GJ, TacsScalar kG22, TacsScalar kG33, TacsScalar kG23,
    TacsScalar m00, TacsScalar m11, TacsScalar m22, TacsScalar m33,
    TacsScalar xm2, TacsScalar xm3, TacsScalar xc2, TacsScalar xc3,
    TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  props = NULL;

  populateMats(EA, EI22, EI33, EI23, GJ, kG22, kG33, kG23, m00, m11, m22, m33,
               xm2, xm3, xc2, xc3, xk2, xk3, muS);
}

void TACSBasicBeamConstitutive::populateMats(
    TacsScalar EA, TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
    TacsScalar GJ, TacsScalar kG22, TacsScalar kG33, TacsScalar kG23,
    TacsScalar m00, TacsScalar m11, TacsScalar m22, TacsScalar m33,
    TacsScalar xm2, TacsScalar xm3, TacsScalar xc2, TacsScalar xc3,
    TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  // Set the entries of the stiffness matrix
  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  // row 1 for axial force
  C[0] = EA;
  C[2] = xc3 * EA;
  C[3] = -xc2 * EA;

  // row 2 for twisting moment
  C[6] = GJ + xk2 * xk2 * kG33 + xk3 * xk3 * kG22 + 2.0 * xk2 * xk3 * kG23;
  C[9] = -xk2 * kG23 - xk3 * kG22;
  C[10] = xk2 * kG33 + xk3 * kG23;

  // row 3 for bending moment about axis 2
  C[11] = EI22 + xc3 * xc3 * EA;
  C[12] = -(EI23 + xc2 * xc3 * EA);

  // row 4 for bending moment about axis 3
  C[15] = EI33 + xc2 * xc2 * EA;

  // row 5 for shear 2
  C[18] = kG22;
  C[19] = -kG23;

  // row 6 for shear 3
  C[20] = kG33;

  // Set the entries of the density matrix
  rho[0] = m00;
  rho[1] = xm2 * m00;
  rho[2] = xm3 * m00;
  rho[3] = m11;
  rho[4] = m22;
  rho[5] = m33;  // m00*xm2*xm3;
}

/*
  Set the diagonal components of the stiffness matrix and the mass
  moments of the cross-section.
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TacsScalar rhoA, TacsScalar rhoIy, TacsScalar rhoIz, TacsScalar rhoIyz,
    TacsScalar EA, TacsScalar GJ, TacsScalar EIy, TacsScalar EIz,
    TacsScalar kGAy, TacsScalar kGAz) {
  props = NULL;

  // Set the entries of the stiffness matrix
  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));
  C[0] = EA;
  C[6] = GJ;
  C[11] = EIz;
  C[15] = EIy;
  C[18] = kGAy;
  C[20] = kGAz;

  // Set the entries of the density matrix
  rho[0] = rhoA;
  rho[1] = 0.0;
  rho[2] = 0.0;
  rho[3] = rhoIy;
  rho[4] = rhoIz;
  rho[5] = rhoIyz;
}

/*
  Set the diagonal components of the stiffness matrix and the mass
  moments of the cross-section.
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TACSMaterialProperties *properties, TacsScalar A, TacsScalar J,
    TacsScalar Iy, TacsScalar Iz, TacsScalar Iyz, TacsScalar ky, TacsScalar kz,
    TacsScalar nsm, TacsScalar xm2, TacsScalar xm3, TacsScalar xc2,
    TacsScalar xc3, TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  props = properties;
  props->incref();

  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);
  const TacsScalar G = 0.5 * E / (1.0 + nu);

  const TacsScalar density = props->getDensity();

  const TacsScalar EA = E * A;
  const TacsScalar EI22 = E * Iz;
  const TacsScalar EI33 = E * Iy;
  const TacsScalar EI23 = E * Iyz;
  const TacsScalar GJ = G * J;
  const TacsScalar kG22 = ky * G * A;
  const TacsScalar kG33 = kz * G * A;
  const TacsScalar kG23 = 0.0;
  // NSM sits at (xm2, xm3) — same location as structural CoM.
  // Full parallel-axis contributions are baked into rho[] here so that
  // evalMassMoments can return rho[] directly without runtime NSM logic.
  const TacsScalar m00 = density * A + nsm;
  const TacsScalar m11 = density * Iz + nsm * xm2 * xm2;
  const TacsScalar m22 = density * Iy + nsm * xm3 * xm3;
  const TacsScalar m33 = -density * Iyz + nsm * xm2 * xm3;

  populateMats(EA, EI22, EI33, EI23, GJ, kG22, kG33, kG23, m00, m11, m22, m33,
               xm2, xm3, xc2, xc3, xk2, xk3, muS);
}

TACSBasicBeamConstitutive::~TACSBasicBeamConstitutive() {
  if (props) {
    props->decref();
  }
}

const char *TACSBasicBeamConstitutive::constName = "TACSBasicBeamConstitutive";

/*
  Return the constitutive name
*/
const char *TACSBasicBeamConstitutive::getObjectName() { return constName; }
