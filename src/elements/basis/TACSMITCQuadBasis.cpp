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

#include "TACSMITCQuadBasis.h"

#include "TACSGaussQuadrature.h"
#include "TACSLagrangeInterpolation.h"

static void getEdgeTangent(int edge, double t[]) {
  if (edge == 0) {
    // -X edge
    t[0] = 0.0;
    t[1] = -1.0;
  } else if (edge == 1) {
    // +X edge
    t[0] = 0.0;
    t[1] = 1.0;
  } else if (edge == 2) {
    // -Y edge
    t[0] = 1.0;
    t[1] = 0.0;
  } else if (edge == 3) {
    // +Y edge
    t[0] = -1.0;
    t[1] = 0.0;
  }
}

/*
  Linear Quad basis class functions
*/
ElementLayout TACSLinearMITCQuadBasis::getLayoutType() {
  return TACS_QUAD_ELEMENT;
}

void TACSLinearMITCQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 2);
  pt[1] = -1.0 + 2.0 * (n / 2);
}

int TACSLinearMITCQuadBasis::getNumNodes() { return 4; }

int TACSLinearMITCQuadBasis::getNumParameters() { return 2; }

int TACSLinearMITCQuadBasis::getNumQuadraturePoints() { return 4; }

double TACSLinearMITCQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
}

double TACSLinearMITCQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts2[n % 2];
  pt[1] = TacsGaussQuadPts2[n / 2];

  return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
}

int TACSLinearMITCQuadBasis::getNumElementFaces() { return 4; }

int TACSLinearMITCQuadBasis::getNumFaceQuadraturePoints(int face) { return 2; }

double TACSLinearMITCQuadBasis::getFaceQuadraturePoint(int face, int n,
                                                       double pt[],
                                                       double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts2[n];
  } else {
    pt[0] = TacsGaussQuadPts2[n];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getEdgeTangent(face, t);

  return TacsGaussQuadWts2[n];
}

void TACSLinearMITCQuadBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);
}

void TACSLinearMITCQuadBasis::computeBasisGradient(const double pt[],
                                                   double N[], double Nxi[]) {
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

  Nxi[0] = -0.25 * (1.0 - pt[1]);
  Nxi[1] = -0.25 * (1.0 - pt[0]);

  Nxi[2] = 0.25 * (1.0 - pt[1]);
  Nxi[3] = -0.25 * (1.0 + pt[0]);

  Nxi[4] = -0.25 * (1.0 + pt[1]);
  Nxi[5] = 0.25 * (1.0 - pt[0]);

  Nxi[6] = 0.25 * (1.0 + pt[1]);
  Nxi[7] = 0.25 * (1.0 + pt[0]);
}

int TACSLinearMITCQuadBasis::getNumTyingFields() { return 2; }

int TACSLinearMITCQuadBasis::getNumTyingPoints(int field) { return 2; }

void TACSLinearMITCQuadBasis::getTyingPoint(int field, int ty, double pt[]) {
  if (field == 0) {
    if (ty == 0) {
      pt[0] = 0.0;
      pt[1] = -1.0;
    } else if (ty == 1) {
      pt[0] = 0.0;
      pt[1] = 1.0;
    }
  } else if (field == 1) {
    if (ty == 0) {
      pt[0] = -1.0;
      pt[1] = 0.0;
    } else if (ty == 1) {
      pt[0] = 1.0;
      pt[1] = 0.0;
    }
  }
}

void TACSLinearMITCQuadBasis::computeTyingBasis(int field, const double pt[],
                                                double N[]) {
  if (field == 0) {
    N[0] = 0.5 * (1.0 - pt[1]);
    N[1] = 0.5 * (1.0 + pt[1]);
  } else if (field == 1) {
    N[0] = 0.5 * (1.0 - pt[0]);
    N[1] = 0.5 * (1.0 + pt[0]);
  }
}

/*
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticMITCQuadBasis::getLayoutType() {
  return TACS_QUAD_QUADRATIC_ELEMENT;
}

void TACSQuadraticMITCQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 1.0 * (n % 3);
  pt[1] = -1.0 + 1.0 * (n / 3);
}

int TACSQuadraticMITCQuadBasis::getNumNodes() { return 9; }

int TACSQuadraticMITCQuadBasis::getNumParameters() { return 2; }

int TACSQuadraticMITCQuadBasis::getNumQuadraturePoints() { return 9; }

double TACSQuadraticMITCQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

double TACSQuadraticMITCQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[n / 3];

  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

int TACSQuadraticMITCQuadBasis::getNumElementFaces() { return 4; }

int TACSQuadraticMITCQuadBasis::getNumFaceQuadraturePoints(int face) {
  return 3;
}

double TACSQuadraticMITCQuadBasis::getFaceQuadraturePoint(int face, int n,
                                                          double pt[],
                                                          double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts3[n];
  } else {
    pt[0] = TacsGaussQuadPts3[n];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getEdgeTangent(face, t);

  return TacsGaussQuadWts3[n];
}

void TACSQuadraticMITCQuadBasis::computeBasis(const double pt[], double N[]) {
  double na[3];
  na[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  na[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  na[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double nb[3];
  nb[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  nb[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  nb[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  N[0] = na[0] * nb[0];
  N[1] = na[1] * nb[0];
  N[2] = na[2] * nb[0];
  N[3] = na[0] * nb[1];
  N[4] = na[1] * nb[1];
  N[5] = na[2] * nb[1];
  N[6] = na[0] * nb[2];
  N[7] = na[1] * nb[2];
  N[8] = na[2] * nb[2];
}

void TACSQuadraticMITCQuadBasis::computeBasisGradient(const double pt[],
                                                      double N[],
                                                      double Nxi[]) {
  double na[3];
  na[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  na[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  na[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double nb[3];
  nb[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  nb[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  nb[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  N[0] = na[0] * nb[0];
  N[1] = na[1] * nb[0];
  N[2] = na[2] * nb[0];
  N[3] = na[0] * nb[1];
  N[4] = na[1] * nb[1];
  N[5] = na[2] * nb[1];
  N[6] = na[0] * nb[2];
  N[7] = na[1] * nb[2];
  N[8] = na[2] * nb[2];

  double dna[3];
  dna[0] = -0.5 + pt[0];
  dna[1] = -2.0 * pt[0];
  dna[2] = 0.5 + pt[0];

  double dnb[3];
  dnb[0] = -0.5 + pt[1];
  dnb[1] = -2.0 * pt[1];
  dnb[2] = 0.5 + pt[1];

  Nxi[0] = dna[0] * nb[0];
  Nxi[1] = na[0] * dnb[0];
  Nxi[2] = dna[1] * nb[0];
  Nxi[3] = na[1] * dnb[0];
  Nxi[4] = dna[2] * nb[0];
  Nxi[5] = na[2] * dnb[0];
  Nxi[6] = dna[0] * nb[1];
  Nxi[7] = na[0] * dnb[1];
  Nxi[8] = dna[1] * nb[1];
  Nxi[9] = na[1] * dnb[1];
  Nxi[10] = dna[2] * nb[1];
  Nxi[11] = na[2] * dnb[1];
  Nxi[12] = dna[0] * nb[2];
  Nxi[13] = na[0] * dnb[2];
  Nxi[14] = dna[1] * nb[2];
  Nxi[15] = na[1] * dnb[2];
  Nxi[16] = dna[2] * nb[2];
  Nxi[17] = na[2] * dnb[2];
}

int TACSQuadraticMITCQuadBasis::getNumTyingFields() { return 2; }

int TACSQuadraticMITCQuadBasis::getNumTyingPoints(int field) { return 6; }

void TACSQuadraticMITCQuadBasis::getTyingPoint(int field, int ty, double pt[]) {
  if (field == 0) {
    pt[0] = TacsGaussQuadPts2[ty % 2];
    pt[1] = -1.0 + 2.0 * (ty / 2);
  } else if (field == 1) {
    pt[0] = -1.0 + 2.0 * (ty / 2);
    pt[1] = TacsGaussQuadPts2[ty % 2];
  }
}

void TACSQuadraticMITCQuadBasis::computeTyingBasis(int field, const double pt[],
                                                   double N[]) {
  if (field == 0) {
    double n0[2];
    n0[0] = (pt[0] - TacsGaussQuadPts2[1]) /
            (TacsGaussQuadPts2[0] - TacsGaussQuadPts2[1]);
    n0[1] = (pt[0] - TacsGaussQuadPts2[0]) /
            (TacsGaussQuadPts2[1] - TacsGaussQuadPts2[0]);

    double n1[3];
    n1[0] = -0.5 * pt[1] * (1.0 - pt[1]);
    n1[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
    n1[2] = 0.5 * (1.0 + pt[1]) * pt[1];

    N[0] = n0[0] * n1[0];
    N[1] = n0[1] * n1[0];
    N[2] = n0[0] * n1[1];
    N[3] = n0[1] * n1[1];
    N[4] = n0[0] * n1[2];
    N[5] = n0[1] * n1[2];
  } else if (field == 1) {
    double n0[2];
    n0[0] = (pt[1] - TacsGaussQuadPts2[1]) /
            (TacsGaussQuadPts2[0] - TacsGaussQuadPts2[1]);
    n0[1] = (pt[1] - TacsGaussQuadPts2[0]) /
            (TacsGaussQuadPts2[1] - TacsGaussQuadPts2[0]);

    double n1[3];
    n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
    n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
    n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

    N[0] = n0[0] * n1[0];
    N[1] = n0[1] * n1[0];
    N[2] = n0[0] * n1[1];
    N[3] = n0[1] * n1[1];
    N[4] = n0[0] * n1[2];
    N[5] = n0[1] * n1[2];
  }
}