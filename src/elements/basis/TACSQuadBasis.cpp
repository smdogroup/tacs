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

#include "TACSQuadBasis.h"

#include "TACSBasisMacros.h"
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
ElementLayout TACSLinearQuadBasis::getLayoutType() { return TACS_QUAD_ELEMENT; }

void TACSLinearQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 2);
  pt[1] = -1.0 + 2.0 * (n / 2);
}

int TACSLinearQuadBasis::getNumNodes() { return 4; }

int TACSLinearQuadBasis::getNumParameters() { return 2; }

int TACSLinearQuadBasis::getNumQuadraturePoints() { return 4; }

double TACSLinearQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
}

double TACSLinearQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts2[n % 2];
  pt[1] = TacsGaussQuadPts2[n / 2];

  return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
}

int TACSLinearQuadBasis::getNumElementFaces() { return 4; }

int TACSLinearQuadBasis::getNumFaceQuadraturePoints(int face) { return 2; }

double TACSLinearQuadBasis::getFaceQuadraturePoint(int face, int n, double pt[],
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

void TACSLinearQuadBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
  N[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
  N[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
  N[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);
}

void TACSLinearQuadBasis::computeBasisGradient(const double pt[], double N[],
                                               double Nxi[]) {
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

/*
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticQuadBasis::getLayoutType() {
  return TACS_QUAD_QUADRATIC_ELEMENT;
}

void TACSQuadraticQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 1.0 * (n % 3);
  pt[1] = -1.0 + 1.0 * (n / 3);
}

int TACSQuadraticQuadBasis::getNumNodes() { return 9; }

int TACSQuadraticQuadBasis::getNumParameters() { return 2; }

int TACSQuadraticQuadBasis::getNumQuadraturePoints() { return 9; }

double TACSQuadraticQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

double TACSQuadraticQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[n / 3];

  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

int TACSQuadraticQuadBasis::getNumElementFaces() { return 4; }

int TACSQuadraticQuadBasis::getNumFaceQuadraturePoints(int face) { return 3; }

double TACSQuadraticQuadBasis::getFaceQuadraturePoint(int face, int n,
                                                      double pt[], double t[]) {
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

void TACSQuadraticQuadBasis::computeBasis(const double pt[], double N[]) {
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

void TACSQuadraticQuadBasis::computeBasisGradient(const double pt[], double N[],
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

/*
  Cubic Quad basis class functions
*/
ElementLayout TACSCubicQuadBasis::getLayoutType() {
  return TACS_QUAD_CUBIC_ELEMENT;
}

void TACSCubicQuadBasis::getVisPoint(int n, double pt[]) {
  static const double pts[4] = {-1.0, -0.5, 0.5, 1.0};
  pt[0] = pts[n % 4];
  pt[1] = pts[n / 4];
}

int TACSCubicQuadBasis::getNumNodes() { return 16; }

int TACSCubicQuadBasis::getNumParameters() { return 2; }

int TACSCubicQuadBasis::getNumQuadraturePoints() { return 16; }

double TACSCubicQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
}

double TACSCubicQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[n / 4];

  return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
}

int TACSCubicQuadBasis::getNumElementFaces() { return 4; }

int TACSCubicQuadBasis::getNumFaceQuadraturePoints(int face) { return 4; }

double TACSCubicQuadBasis::getFaceQuadraturePoint(int face, int n, double pt[],
                                                  double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts4[n];
  } else {
    pt[0] = TacsGaussQuadPts4[n];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getEdgeTangent(face, t);

  return TacsGaussQuadWts4[n];
}

void TACSCubicQuadBasis::computeBasis(const double pt[], double N[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  N[0] = na[0] * nb[0];
  N[1] = na[1] * nb[0];
  N[2] = na[2] * nb[0];
  N[3] = na[3] * nb[0];
  N[4] = na[0] * nb[1];
  N[5] = na[1] * nb[1];
  N[6] = na[2] * nb[1];
  N[7] = na[3] * nb[1];
  N[8] = na[0] * nb[2];
  N[9] = na[1] * nb[2];
  N[10] = na[2] * nb[2];
  N[11] = na[3] * nb[2];
  N[12] = na[0] * nb[3];
  N[13] = na[1] * nb[3];
  N[14] = na[2] * nb[3];
  N[15] = na[3] * nb[3];
}

void TACSCubicQuadBasis::computeBasisGradient(const double pt[], double N[],
                                              double Nxi[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  N[0] = na[0] * nb[0];
  N[1] = na[1] * nb[0];
  N[2] = na[2] * nb[0];
  N[3] = na[3] * nb[0];
  N[4] = na[0] * nb[1];
  N[5] = na[1] * nb[1];
  N[6] = na[2] * nb[1];
  N[7] = na[3] * nb[1];
  N[8] = na[0] * nb[2];
  N[9] = na[1] * nb[2];
  N[10] = na[2] * nb[2];
  N[11] = na[3] * nb[2];
  N[12] = na[0] * nb[3];
  N[13] = na[1] * nb[3];
  N[14] = na[2] * nb[3];
  N[15] = na[3] * nb[3];

  double dna[4];
  dna[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  dna[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  dna[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  dna[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double dnb[4];
  dnb[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  dnb[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  dnb[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  dnb[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  Nxi[0] = dna[0] * nb[0];
  Nxi[1] = na[0] * dnb[0];

  Nxi[2] = dna[1] * nb[0];
  Nxi[3] = na[1] * dnb[0];

  Nxi[4] = dna[2] * nb[0];
  Nxi[5] = na[2] * dnb[0];

  Nxi[6] = dna[3] * nb[0];
  Nxi[7] = na[3] * dnb[0];

  Nxi[8] = dna[0] * nb[1];
  Nxi[9] = na[0] * dnb[1];

  Nxi[10] = dna[1] * nb[1];
  Nxi[11] = na[1] * dnb[1];

  Nxi[12] = dna[2] * nb[1];
  Nxi[13] = na[2] * dnb[1];

  Nxi[14] = dna[3] * nb[1];
  Nxi[15] = na[3] * dnb[1];

  Nxi[16] = dna[0] * nb[2];
  Nxi[17] = na[0] * dnb[2];

  Nxi[18] = dna[1] * nb[2];
  Nxi[19] = na[1] * dnb[2];

  Nxi[20] = dna[2] * nb[2];
  Nxi[21] = na[2] * dnb[2];

  Nxi[22] = dna[3] * nb[2];
  Nxi[23] = na[3] * dnb[2];

  Nxi[24] = dna[0] * nb[3];
  Nxi[25] = na[0] * dnb[3];

  Nxi[26] = dna[1] * nb[3];
  Nxi[27] = na[1] * dnb[3];

  Nxi[28] = dna[2] * nb[3];
  Nxi[29] = na[2] * dnb[3];

  Nxi[30] = dna[3] * nb[3];
  Nxi[31] = na[3] * dnb[3];
}

void TACSCubicQuadBasis::interpFields(const int n, const double pt[],
                                      const int m, const TacsScalar v[],
                                      const int incr, TacsScalar u[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  for (int i = 0; i < m; i++, u += incr, v++) {
    u[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER4(na, nb, m, v);
  }
}

void TACSCubicQuadBasis::addInterpFieldsTranspose(const int n,
                                                  const double pt[],
                                                  const int incr,
                                                  const TacsScalar u[],
                                                  const int m, TacsScalar v[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  for (int i = 0; i < m; i++, u += incr, v++) {
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER4(na, nb, u[0], temp, m, v);
  }
}

void TACSCubicQuadBasis::interpFieldsGrad(const int n, const double pt[],
                                          const int m, const TacsScalar v[],
                                          TacsScalar g[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double dna[4];
  dna[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  dna[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  dna[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  dna[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double dnb[4];
  dnb[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  dnb[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  dnb[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  dnb[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  for (int i = 0; i < m; i++, g += 2, v++) {
    g[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER4(dna, nb, m, v);
    g[1] = TACS_BASIS_EVAL_TENSOR2D_ORDER4(na, dnb, m, v);
  }
}

void TACSCubicQuadBasis::addInterpFieldsGradTranspose(int n, const double pt[],
                                                      const int m,
                                                      const TacsScalar g[],
                                                      TacsScalar v[]) {
  double na[4];
  na[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  na[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  na[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double nb[4];
  nb[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  nb[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  nb[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double dna[4];
  dna[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  dna[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  dna[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  dna[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double dnb[4];
  dnb[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  dnb[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  dnb[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  dnb[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  for (int i = 0; i < m; i++, g += 2, v++) {
    TacsScalar temp;
    TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER4(dna, nb, g[0], temp, m, v);
    TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER4(na, dnb, g[1], temp, m, v);
  }
}

/*
  Quartic Quad basis class functions
*/
TACSQuarticQuadBasis::TACSQuarticQuadBasis() {
  for (int i = 0; i < 5; i++) {
    TacsLagrangeShapeFuncDerivative(5, TacsGaussQuadPts5[i], cosine_pts,
                                    &Nf[5 * i], &Nfxi[5 * i]);
  }
}

const double TACSQuarticQuadBasis::cosine_pts[5] = {
    -1.0, -0.7071067811865475, 0.0, 0.7071067811865475, 1.0};

ElementLayout TACSQuarticQuadBasis::getLayoutType() {
  return TACS_QUAD_QUARTIC_ELEMENT;
}

void TACSQuarticQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = cosine_pts[n % 5];
  pt[1] = cosine_pts[n / 5];
}

int TACSQuarticQuadBasis::getNumNodes() { return 25; }

int TACSQuarticQuadBasis::getNumParameters() { return 2; }

int TACSQuarticQuadBasis::getNumQuadraturePoints() { return 25; }

double TACSQuarticQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[n / 5];
}

double TACSQuarticQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts5[n % 5];
  pt[1] = TacsGaussQuadPts5[n / 5];

  return TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[n / 5];
}

int TACSQuarticQuadBasis::getNumElementFaces() { return 4; }

int TACSQuarticQuadBasis::getNumFaceQuadraturePoints(int face) { return 5; }

double TACSQuarticQuadBasis::getFaceQuadraturePoint(int face, int n,
                                                    double pt[], double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts5[n];
  } else {
    pt[0] = TacsGaussQuadPts5[n];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getEdgeTangent(face, t);

  return TacsGaussQuadWts5[n];
}

void TACSQuarticQuadBasis::computeBasis(const double pt[], double N[]) {
  double na[5], nb[5];
  TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, nb);

  for (int j = 0; j < 5; j++) {
    for (int i = 0; i < 5; i++) {
      N[i + 5 * j] = na[i] * nb[j];
    }
  }
}

void TACSQuarticQuadBasis::computeBasisGradient(const double pt[], double N[],
                                                double Nxi[]) {
  double na[5], nb[5];
  double dna[5], dnb[5];
  TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, nb, dnb);

  for (int j = 0; j < 5; j++) {
    for (int i = 0; i < 5; i++) {
      N[i + 5 * j] = na[i] * nb[j];
      Nxi[2 * (i + 5 * j)] = dna[i] * nb[j];
      Nxi[2 * (i + 5 * j) + 1] = na[i] * dnb[j];
    }
  }
}

void TACSQuarticQuadBasis::interpFields(const int n, const double pt[],
                                        const int m, const TacsScalar v[],
                                        const int incr, TacsScalar u[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * (n / 5)];

    for (int i = 0; i < m; i++, u += incr, v++) {
      u[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1, n2, m, v);
    }
  } else {
    double n1[5], n2[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);

    for (int i = 0; i < m; i++, u += incr, v++) {
      u[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1, n2, m, v);
    }
  }
}

void TACSQuarticQuadBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * (n / 5)];

    for (int i = 0; i < m; i++, u += incr, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1, n2, u[0], temp, m, v);
    }
  } else {
    double n1[5], n2[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);

    for (int i = 0; i < m; i++, u += incr, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1, n2, u[0], temp, m, v);
    }
  }
}

void TACSQuarticQuadBasis::interpFieldsGrad(const int n, const double pt[],
                                            const int m, const TacsScalar v[],
                                            TacsScalar g[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * (n / 5)];
    const double *n1xi = &Nfxi[5 * (n % 5)];
    const double *n2xi = &Nfxi[5 * (n / 5)];

    for (int i = 0; i < m; i++, g += 2, v++) {
      g[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1xi, n2, m, v);
      g[1] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1, n2xi, m, v);
    }
  } else {
    double n1[5], n2[5], n1xi[5], n2xi[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1xi);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2xi);

    for (int i = 0; i < m; i++, g += 2, v++) {
      g[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1xi, n2, m, v);
      g[1] = TACS_BASIS_EVAL_TENSOR2D_ORDER5(n1, n2xi, m, v);
    }
  }
}

void TACSQuarticQuadBasis::addInterpFieldsGradTranspose(int n,
                                                        const double pt[],
                                                        const int m,
                                                        const TacsScalar g[],
                                                        TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * (n / 5)];
    const double *n1xi = &Nfxi[5 * (n % 5)];
    const double *n2xi = &Nfxi[5 * (n / 5)];

    for (int i = 0; i < m; i++, g += 2, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1xi, n2, g[0], temp, m, v);
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1, n2xi, g[1], temp, m, v);
    }
  } else {
    double n1[5], n2[5], n1xi[5], n2xi[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1xi);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2xi);

    for (int i = 0; i < m; i++, g += 2, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1xi, n2, g[0], temp, m, v);
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER5(n1, n2xi, g[1], temp, m, v);
    }
  }
}

/*
  Quintic Quad basis class functions
*/
TACSQuinticQuadBasis::TACSQuinticQuadBasis() {
  for (int i = 0; i < 6; i++) {
    TacsLagrangeShapeFuncDerivative(6, TacsGaussQuadPts6[i], cosine_pts,
                                    &Nf[6 * i], &Nfxi[6 * i]);
  }
}

const double TACSQuinticQuadBasis::cosine_pts[6] = {-1.0,
                                                    -0.8090169943749475,
                                                    -0.30901699437494745,
                                                    0.30901699437494745,
                                                    0.8090169943749475,
                                                    1.0};

ElementLayout TACSQuinticQuadBasis::getLayoutType() {
  return TACS_QUAD_QUINTIC_ELEMENT;
}

void TACSQuinticQuadBasis::getVisPoint(int n, double pt[]) {
  pt[0] = cosine_pts[n % 6];
  pt[1] = cosine_pts[n / 6];
}

int TACSQuinticQuadBasis::getNumNodes() { return 36; }

int TACSQuinticQuadBasis::getNumParameters() { return 2; }

int TACSQuinticQuadBasis::getNumQuadraturePoints() { return 36; }

double TACSQuinticQuadBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[n / 6];
}

double TACSQuinticQuadBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts6[n % 6];
  pt[1] = TacsGaussQuadPts6[n / 6];

  return TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[n / 6];
}

int TACSQuinticQuadBasis::getNumElementFaces() { return 4; }

int TACSQuinticQuadBasis::getNumFaceQuadraturePoints(int face) { return 6; }

double TACSQuinticQuadBasis::getFaceQuadraturePoint(int face, int n,
                                                    double pt[], double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts6[n];
  } else {
    pt[0] = TacsGaussQuadPts6[n];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getEdgeTangent(face, t);

  return TacsGaussQuadWts6[n];
}

void TACSQuinticQuadBasis::computeBasis(const double pt[], double N[]) {
  double na[6], nb[6];
  TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, nb);

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++) {
      N[i + 6 * j] = na[i] * nb[j];
    }
  }
}

void TACSQuinticQuadBasis::computeBasisGradient(const double pt[], double N[],
                                                double Nxi[]) {
  double na[6], nb[6];
  double dna[6], dnb[6];
  TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, nb, dnb);

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++) {
      N[i + 6 * j] = na[i] * nb[j];
      Nxi[2 * (i + 6 * j)] = dna[i] * nb[j];
      Nxi[2 * (i + 6 * j) + 1] = na[i] * dnb[j];
    }
  }
}

void TACSQuinticQuadBasis::interpFields(const int n, const double pt[],
                                        const int m, const TacsScalar v[],
                                        const int incr, TacsScalar u[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * (n / 6)];

    for (int i = 0; i < m; i++, u += incr, v++) {
      u[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1, n2, m, v);
    }
  } else {
    double n1[6], n2[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);

    for (int i = 0; i < m; i++, u += incr, v++) {
      u[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1, n2, m, v);
    }
  }
}

void TACSQuinticQuadBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * (n / 6)];

    for (int i = 0; i < m; i++, u += incr, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1, n2, u[0], temp, m, v);
    }
  } else {
    double n1[6], n2[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);

    for (int i = 0; i < m; i++, u += incr, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1, n2, u[0], temp, m, v);
    }
  }
}

void TACSQuinticQuadBasis::interpFieldsGrad(const int n, const double pt[],
                                            const int m, const TacsScalar v[],
                                            TacsScalar g[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * (n / 6)];
    const double *n1xi = &Nfxi[6 * (n % 6)];
    const double *n2xi = &Nfxi[6 * (n / 6)];

    for (int i = 0; i < m; i++, g += 2, v++) {
      g[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1xi, n2, m, v);
      g[1] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1, n2xi, m, v);
    }
  } else {
    double n1[6], n2[6], n1xi[6], n2xi[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1xi);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2xi);

    for (int i = 0; i < m; i++, g += 2, v++) {
      g[0] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1xi, n2, m, v);
      g[1] = TACS_BASIS_EVAL_TENSOR2D_ORDER6(n1, n2xi, m, v);
    }
  }
}

void TACSQuinticQuadBasis::addInterpFieldsGradTranspose(int n,
                                                        const double pt[],
                                                        const int m,
                                                        const TacsScalar g[],
                                                        TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * (n / 6)];
    const double *n1xi = &Nfxi[6 * (n % 6)];
    const double *n2xi = &Nfxi[6 * (n / 6)];

    for (int i = 0; i < m; i++, g += 2, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1xi, n2, g[0], temp, m, v);
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1, n2xi, g[1], temp, m, v);
    }
  } else {
    double n1[6], n2[6], n1xi[6], n2xi[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1xi);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2xi);

    for (int i = 0; i < m; i++, g += 2, v++) {
      TacsScalar temp;
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1xi, n2, g[0], temp, m, v);
      TACS_BASIS_TRANSPOSE_TENSOR2D_ORDER6(n1, n2xi, g[1], temp, m, v);
    }
  }
}