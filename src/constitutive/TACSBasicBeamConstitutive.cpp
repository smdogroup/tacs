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
  kGA22, kGA33, kGA23 : shearing stiffness
  m00              : mass per unit span
  m11, m22, m33    : second mass moments of inertia (m11 = m22 + m33)
  xm2, xm3         : combined center of mass (positive geometric y, z)
  xc2, xc3         : structural centroid (xc2 = -z_centroid, xc3 = +y_centroid)
  xk2, xk3         : shear center
  muS              : viscous damping coefficient
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TacsScalar EA, TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
    TacsScalar GJ, TacsScalar kGA22, TacsScalar kGA33, TacsScalar kGA23,
    TacsScalar m00, TacsScalar m11, TacsScalar m22, TacsScalar m33,
    TacsScalar xm2, TacsScalar xm3, TacsScalar xc2, TacsScalar xc3,
    TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  props = NULL;

  populateMats(EA, EI22, EI33, EI23, GJ, kGA22, kGA33, kGA23, m00, xm2 * m00,
               xm3 * m00, m11, m22, m33, xc2, xc3, xk2, xk3, muS);
}

void TACSBasicBeamConstitutive::populateMats(
    TacsScalar EA, TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
    TacsScalar GJ, TacsScalar kGA22, TacsScalar kGA33, TacsScalar kGA23,
    TacsScalar m00, TacsScalar m1, TacsScalar m2, TacsScalar m11,
    TacsScalar m22, TacsScalar m33, TacsScalar xc2, TacsScalar xc3,
    TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  // Set the entries of the stiffness matrix
  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));

  C[0] = EA;
  C[2] = EA * xc2;
  C[3] = EA * xc3;
  C[6] = GJ + xk2 * (kGA23 * xk3 + kGA33 * xk2) +
         xk3 * (kGA22 * xk3 + kGA23 * xk2);
  C[9] = -kGA22 * xk3 - kGA23 * xk2;
  C[10] = kGA23 * xk3 + kGA33 * xk2;
  C[11] = EA * xc2 * xc2 + EI33;
  C[12] = EA * xc2 * xc3 + EI23;
  C[15] = EA * xc3 * xc3 + EI22;
  C[18] = kGA22;
  C[19] = -kGA23;
  C[20] = kGA33;

  // Set the entries of the density matrix
  rho[0] = m00;
  rho[1] = m1;
  rho[2] = m2;
  rho[3] = m11;
  rho[4] = m22;
  rho[5] = m33;
}

/*
  Set the diagonal components of the stiffness matrix and the mass
  moments of the cross-section.
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TacsScalar rhoA, TacsScalar rhoI22, TacsScalar rhoI33, TacsScalar rhoI23,
    TacsScalar EA, TacsScalar GJ, TacsScalar EI22, TacsScalar EI33,
    TacsScalar kGA2, TacsScalar kGA3) {
  props = NULL;

  // Set the entries of the stiffness matrix
  memset(C, 0, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));
  C[0] = EA;
  C[6] = GJ;
  C[11] = EI33;
  C[15] = EI22;
  C[18] = kGA2;
  C[20] = kGA3;

  // Set the entries of the density matrix
  rho[0] = rhoA;
  rho[1] = 0.0;
  rho[2] = 0.0;
  rho[3] = rhoI22;
  rho[4] = rhoI33;
  rho[5] = rhoI23;
}

/*
  Set the diagonal components of the stiffness matrix and the mass
  moments of the cross-section.
*/
TACSBasicBeamConstitutive::TACSBasicBeamConstitutive(
    TACSMaterialProperties *properties, TacsScalar A, TacsScalar J,
    TacsScalar I22, TacsScalar I33, TacsScalar I23, TacsScalar k2,
    TacsScalar k3, TacsScalar nsm, TacsScalar xm2, TacsScalar xm3,
    TacsScalar xnsm2, TacsScalar xnsm3, TacsScalar xc2, TacsScalar xc3,
    TacsScalar xk2, TacsScalar xk3, TacsScalar muS) {
  props = properties;
  props->incref();

  TacsScalar E, nu;
  props->getIsotropicProperties(&E, &nu);
  const TacsScalar G = 0.5 * E / (1.0 + nu);

  const TacsScalar density = props->getDensity();

  const TacsScalar EA = E * A;
  const TacsScalar EI22 = E * I22;
  const TacsScalar EI33 = E * I33;
  const TacsScalar EI23 = E * I23;
  const TacsScalar GJ = G * J;
  const TacsScalar kGA22 = k2 * G * A;
  const TacsScalar kGA33 = k3 * G * A;
  const TacsScalar kGA23 = 0.0;

  // xm2, xm3   : structural CoM (positive geometric y, z)
  // xnsm2, xnsm3 : NSM CoM (positive geometric y, z); default 0 = at reference
  // axis
  const TacsScalar m00 = density * A + nsm;
  // First mass moments: structural + NSM contributions combined
  const TacsScalar m1 = density * A * xm2 + nsm * xnsm2;
  const TacsScalar m2 = density * A * xm3 + nsm * xnsm3;
  // Second mass moments about the reference axis (parallel-axis for each piece)
  const TacsScalar m11 =
      density * I33 + density * A * xm2 * xm2 + nsm * xnsm2 * xnsm2;
  const TacsScalar m22 =
      density * I22 + density * A * xm3 * xm3 + nsm * xnsm3 * xnsm3;
  const TacsScalar m33 =
      -density * I23 + density * A * xm2 * xm3 + nsm * xnsm2 * xnsm3;

  populateMats(EA, EI22, EI33, EI23, GJ, kGA22, kGA33, kGA23, m00, m1, m2, m11,
               m22, m33, xc2, xc3, xk2, xk3, muS);
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
