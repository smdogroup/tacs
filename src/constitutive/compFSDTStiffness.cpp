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

#include "compFSDTStiffness.h"

/*
  Create the stiffness object that consists of a series of
  orthotropic material plies.

  Note that this constitutive class does not define any design
  variables.
*/
compFSDTStiffness::compFSDTStiffness(OrthoPly **_ortho_ply, TacsScalar _kcorr,
                                     TacsScalar *_thickness,
                                     TacsScalar *_ply_angles, int _num_plies) {
  num_plies = _num_plies;
  ortho_ply = new OrthoPly *[num_plies];
  for (int i = 0; i < num_plies; i++) {
    ortho_ply[i] = _ortho_ply[i];
  }

  kcorr = _kcorr;
  thickness = new TacsScalar[num_plies];
  memcpy(thickness, _thickness, num_plies * sizeof(TacsScalar));

  ply_angles = new TacsScalar[num_plies];
  memcpy(ply_angles, _ply_angles, num_plies * sizeof(TacsScalar));
}

/*
  Free/dereference variables that were used by compFSDTStiffness
*/
compFSDTStiffness::~compFSDTStiffness() {
  if (ortho_ply) {
    for (int i = 0; i < num_plies; i++) {
      ortho_ply[i]->decref();
    }
  }
  if (thickness) {
    delete[] thickness;
  }
  if (ply_angles) {
    delete[] ply_angles;
  }
}

/*
  Set/get the name of the material property
*/
const char *compFSDTStiffness::constName = "compFSDTStiffness";

const char *compFSDTStiffness::constitutiveName() { return constName; }

/*
  Compute the laminate stiffness matricies given the specified
  lamination sequence and material types.
*/
TacsScalar compFSDTStiffness::getStiffness(const double gpt[], TacsScalar A[],
                                           TacsScalar B[], TacsScalar D[],
                                           TacsScalar As[]) {
  // Zero the stiffness matrices
  for (int k = 0; k < 6; k++) {
    A[k] = B[k] = D[k] = 0.0;
  }

  for (int k = 0; k < 3; k++) {
    As[k] = 0.0;
  }

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += thickness[i];
  }

  // Compute the contribution to the stiffness from each layer
  TacsScalar t0 = -0.5 * t;
  for (int k = 0; k < num_plies; k++) {
    TacsScalar Qbar[6], Abar[3];
    ortho_ply[k]->calculateQbar(Qbar, ply_angles[k]);
    ortho_ply[k]->calculateAbar(Abar, ply_angles[k]);

    TacsScalar t1 = t0 + thickness[k];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5 * (t1 * t1 - t0 * t0);
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    for (int i = 0; i < 6; i++) {
      A[i] += a * Qbar[i];
      B[i] += b * Qbar[i];
      D[i] += d * Qbar[i];
    }

    for (int i = 0; i < 3; i++) {
      As[i] += kcorr * a * Abar[i];
    }

    // Update the position of the bottom interface
    t0 = t1;
  }

  return 0.5 * DRILLING_REGULARIZATION * (As[0] + As[2]);
}

void compFSDTStiffness::getPointwiseMass(const double gpt[],
                                         TacsScalar mass[]) {
  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += thickness[i];
  }

  // zero the moments of inertia
  mass[0] = mass[1] = mass[2] = 0.0;

  // Compute the mass properties
  TacsScalar t0 = -0.5 * t;
  for (int i = 0; i < num_plies; i++) {
    TacsScalar t1 = t0 + thickness[i];
    TacsScalar d = 1.0 / 3.0 * (t1 * t1 * t1 - t0 * t0 * t0);

    mass[0] += thickness[i] * ortho_ply[i]->getRho();
    mass[2] += d * ortho_ply[i]->getRho() / 12.0;
  }
}

/*
  Compute the most critical failure criteria for the laminate
*/
void compFSDTStiffness::failure(const double pt[], const TacsScalar strain[],
                                TacsScalar *fail) {
  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for (int i = 0; i < num_plies; i++) {
    t += thickness[i];
  }
  TacsScalar t0 = -0.5 * t;

  // Keep track of the maximum failure criterion
  TacsScalar max = 0.0;

  for (int i = 0; i < num_plies; i++) {
    TacsScalar lamStrain[3];
    TacsScalar tp = t0 + 0.5 * thickness[i];
    getLaminaStrain(lamStrain, strain, tp);
    TacsScalar fval = ortho_ply[i]->failure(ply_angles[i], lamStrain);

    if (TacsRealPart(fval) > TacsRealPart(max)) {
      max = fval;
    }
    t0 += thickness[i];
  }

  *fail = max;
}

/*
  Get the strain in a single ply
*/
void compFSDTStiffness::getLaminaStrain(TacsScalar strain[],
                                        const TacsScalar rmStrain[],
                                        TacsScalar tp) {
  strain[0] = rmStrain[0] + tp * rmStrain[3];
  strain[1] = rmStrain[1] + tp * rmStrain[4];
  strain[2] = rmStrain[2] + tp * rmStrain[5];
}
