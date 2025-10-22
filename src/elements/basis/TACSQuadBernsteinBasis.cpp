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

#include "TACSQuadBernsteinBasis.h"

#include "TACSBernsteinInterpolation.h"
#include "TACSGaussQuadrature.h"

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
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticQuadBernsteinBasis::getLayoutType() {
  return TACS_QUAD_QUADRATIC_ELEMENT;
}

void TACSQuadraticQuadBernsteinBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 1.0 * (n % 3);
  pt[1] = -1.0 + 1.0 * (n / 3);
}

int TACSQuadraticQuadBernsteinBasis::getNumNodes() { return 9; }

int TACSQuadraticQuadBernsteinBasis::getNumParameters() { return 2; }

int TACSQuadraticQuadBernsteinBasis::getNumQuadraturePoints() { return 9; }

double TACSQuadraticQuadBernsteinBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

double TACSQuadraticQuadBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[n / 3];

  return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
}

int TACSQuadraticQuadBernsteinBasis::getNumElementFaces() { return 4; }

int TACSQuadraticQuadBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 3;
}

double TACSQuadraticQuadBernsteinBasis::getFaceQuadraturePoint(int face, int n,
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

void TACSQuadraticQuadBernsteinBasis::computeBasis(const double pt[],
                                                   double N[]) {
  double na[3], nb[3];
  TacsBernsteinShapeFunctions(3, pt[0], na);
  TacsBernsteinShapeFunctions(3, pt[1], nb);

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

void TACSQuadraticQuadBernsteinBasis::computeBasisGradient(const double pt[],
                                                           double N[],
                                                           double Nxi[]) {
  double na[3], nb[3];
  double dna[3], dnb[3];
  TacsBernsteinShapeFuncDerivative(3, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(3, pt[1], nb, dnb);

  N[0] = na[0] * nb[0];
  N[1] = na[1] * nb[0];
  N[2] = na[2] * nb[0];
  N[3] = na[0] * nb[1];
  N[4] = na[1] * nb[1];
  N[5] = na[2] * nb[1];
  N[6] = na[0] * nb[2];
  N[7] = na[1] * nb[2];
  N[8] = na[2] * nb[2];

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
ElementLayout TACSCubicQuadBernsteinBasis::getLayoutType() {
  return TACS_QUAD_CUBIC_ELEMENT;
}

void TACSCubicQuadBernsteinBasis::getVisPoint(int n, double pt[]) {
  static const double pts[4] = {-1.0, -0.5, 0.5, 1.0};
  pt[0] = pts[n % 4];
  pt[1] = pts[n / 4];
}

int TACSCubicQuadBernsteinBasis::getNumNodes() { return 16; }

int TACSCubicQuadBernsteinBasis::getNumParameters() { return 2; }

int TACSCubicQuadBernsteinBasis::getNumQuadraturePoints() { return 16; }

double TACSCubicQuadBernsteinBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
}

double TACSCubicQuadBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[n / 4];

  return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
}

int TACSCubicQuadBernsteinBasis::getNumElementFaces() { return 4; }

int TACSCubicQuadBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 4;
}

double TACSCubicQuadBernsteinBasis::getFaceQuadraturePoint(int face, int n,
                                                           double pt[],
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

void TACSCubicQuadBernsteinBasis::computeBasis(const double pt[], double N[]) {
  double na[4], nb[4];
  TacsBernsteinShapeFunctions(4, pt[0], na);
  TacsBernsteinShapeFunctions(4, pt[1], nb);

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

void TACSCubicQuadBernsteinBasis::computeBasisGradient(const double pt[],
                                                       double N[],
                                                       double Nxi[]) {
  double na[4], nb[4];
  double dna[4], dnb[4];
  TacsBernsteinShapeFuncDerivative(4, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(4, pt[1], nb, dnb);

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

/*
  Quartic Quad basis class functions
*/
ElementLayout TACSQuarticQuadBernsteinBasis::getLayoutType() {
  return TACS_QUAD_QUARTIC_ELEMENT;
}

void TACSQuarticQuadBernsteinBasis::getVisPoint(int n, double pt[]) {
  static const double pts[5] = {-1.0, -0.7071067811865475, 0.0,
                                0.7071067811865475, 1.0};
  pt[0] = pts[n % 5];
  pt[1] = pts[n / 5];
}

int TACSQuarticQuadBernsteinBasis::getNumNodes() { return 25; }

int TACSQuarticQuadBernsteinBasis::getNumParameters() { return 2; }

int TACSQuarticQuadBernsteinBasis::getNumQuadraturePoints() { return 25; }

double TACSQuarticQuadBernsteinBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[n / 5];
}

double TACSQuarticQuadBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts5[n % 5];
  pt[1] = TacsGaussQuadPts5[n / 5];

  return TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[n / 5];
}

int TACSQuarticQuadBernsteinBasis::getNumElementFaces() { return 4; }

int TACSQuarticQuadBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 5;
}

double TACSQuarticQuadBernsteinBasis::getFaceQuadraturePoint(int face, int n,
                                                             double pt[],
                                                             double t[]) {
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

void TACSQuarticQuadBernsteinBasis::computeBasis(const double pt[],
                                                 double N[]) {
  double na[5], nb[5];
  TacsBernsteinShapeFunctions(5, pt[0], na);
  TacsBernsteinShapeFunctions(5, pt[1], nb);

  for (int j = 0; j < 5; j++) {
    for (int i = 0; i < 5; i++) {
      N[i + 5 * j] = na[i] * nb[j];
    }
  }
}

void TACSQuarticQuadBernsteinBasis::computeBasisGradient(const double pt[],
                                                         double N[],
                                                         double Nxi[]) {
  double na[5], nb[5];
  double dna[5], dnb[5];
  TacsBernsteinShapeFuncDerivative(5, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(5, pt[1], nb, dnb);

  for (int j = 0; j < 5; j++) {
    for (int i = 0; i < 5; i++) {
      N[i + 5 * j] = na[i] * nb[j];
      Nxi[2 * (i + 5 * j)] = dna[i] * nb[j];
      Nxi[2 * (i + 5 * j) + 1] = na[i] * dnb[j];
    }
  }
}

/*
  Quintic Quad basis class functions
*/
ElementLayout TACSQuinticQuadBernsteinBasis::getLayoutType() {
  return TACS_QUAD_QUINTIC_ELEMENT;
}

void TACSQuinticQuadBernsteinBasis::getVisPoint(int n, double pt[]) {
  static const double pts[6] = {-1.0,
                                -0.8090169943749475,
                                -0.30901699437494745,
                                0.30901699437494745,
                                0.8090169943749475,
                                1.0};
  pt[0] = pts[n % 6];
  pt[1] = pts[n / 6];
}

int TACSQuinticQuadBernsteinBasis::getNumNodes() { return 36; }

int TACSQuinticQuadBernsteinBasis::getNumParameters() { return 2; }

int TACSQuinticQuadBernsteinBasis::getNumQuadraturePoints() { return 36; }

double TACSQuinticQuadBernsteinBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[n / 6];
}

double TACSQuinticQuadBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts6[n % 6];
  pt[1] = TacsGaussQuadPts6[n / 6];

  return TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[n / 6];
}

int TACSQuinticQuadBernsteinBasis::getNumElementFaces() { return 4; }

int TACSQuinticQuadBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 6;
}

double TACSQuinticQuadBernsteinBasis::getFaceQuadraturePoint(int face, int n,
                                                             double pt[],
                                                             double t[]) {
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

void TACSQuinticQuadBernsteinBasis::computeBasis(const double pt[],
                                                 double N[]) {
  double na[6], nb[6];
  TacsBernsteinShapeFunctions(6, pt[0], na);
  TacsBernsteinShapeFunctions(6, pt[1], nb);

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++) {
      N[i + 6 * j] = na[i] * nb[j];
    }
  }
}

void TACSQuinticQuadBernsteinBasis::computeBasisGradient(const double pt[],
                                                         double N[],
                                                         double Nxi[]) {
  double na[6], nb[6];
  double dna[6], dnb[6];
  TacsBernsteinShapeFuncDerivative(6, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(6, pt[1], nb, dnb);

  for (int j = 0; j < 6; j++) {
    for (int i = 0; i < 6; i++) {
      N[i + 6 * j] = na[i] * nb[j];
      Nxi[2 * (i + 6 * j)] = dna[i] * nb[j];
      Nxi[2 * (i + 6 * j) + 1] = na[i] * dnb[j];
    }
  }
}
