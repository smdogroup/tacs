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

#include "TACSLineBasis.h"

#include "TACSGaussQuadrature.h"

/*
  Linear Quad basis class functions
*/
ElementLayout TACSLinearLineBasis::getLayoutType() { return TACS_LINE_ELEMENT; }

void TACSLinearLineBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * n;
}

int TACSLinearLineBasis::getNumNodes() { return 2; }

int TACSLinearLineBasis::getNumParameters() { return 1; }

int TACSLinearLineBasis::getNumQuadraturePoints() { return 2; }

double TACSLinearLineBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts2[n];
}

double TACSLinearLineBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts2[n];
  return TacsGaussQuadWts2[n];
}

int TACSLinearLineBasis::getNumElementFaces() { return 2; }

int TACSLinearLineBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSLinearLineBasis::getFaceQuadraturePoint(int face, int n, double pt[],
                                                   double t[]) {
  if (face == 0) {
    t[0] = -1.0;
    pt[0] = -1.0;
  } else {
    t[0] = 1.0;
    pt[1] = 1.0;
  }

  return 1.0;
}

void TACSLinearLineBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 0.5 * (1.0 - pt[0]);
  N[1] = 0.5 * (1.0 + pt[0]);
}

void TACSLinearLineBasis::computeBasisGradient(const double pt[], double N[],
                                               double Nxi[]) {
  N[0] = 0.5 * (1.0 - pt[0]);
  N[1] = 0.5 * (1.0 + pt[0]);

  Nxi[0] = -0.5;
  Nxi[1] = 0.5;
}

/*
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticLineBasis::getLayoutType() {
  return TACS_LINE_QUADRATIC_ELEMENT;
}

void TACSQuadraticLineBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 1.0 * n;
}

int TACSQuadraticLineBasis::getNumNodes() { return 3; }

int TACSQuadraticLineBasis::getNumParameters() { return 1; }

int TACSQuadraticLineBasis::getNumQuadraturePoints() { return 3; }

double TACSQuadraticLineBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts3[n];
}

double TACSQuadraticLineBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n];
  return TacsGaussQuadWts3[n];
}

int TACSQuadraticLineBasis::getNumElementFaces() { return 2; }

int TACSQuadraticLineBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSQuadraticLineBasis::getFaceQuadraturePoint(int face, int n,
                                                      double pt[], double t[]) {
  if (face == 0) {
    t[0] = -1.0;
    pt[0] = -1.0;
  } else {
    t[0] = 1.0;
    pt[1] = 1.0;
  }

  return 1.0;
}

void TACSQuadraticLineBasis::computeBasis(const double pt[], double N[]) {
  N[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  N[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  N[2] = 0.5 * (1.0 + pt[0]) * pt[0];
}

void TACSQuadraticLineBasis::computeBasisGradient(const double pt[], double N[],
                                                  double Nxi[]) {
  N[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  N[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  N[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  Nxi[0] = -0.5 + pt[0];
  Nxi[1] = -2.0 * pt[0];
  Nxi[2] = 0.5 + pt[0];
}

/*
  Cubic Quad basis class functions
*/
ElementLayout TACSCubicLineBasis::getLayoutType() {
  return TACS_LINE_CUBIC_ELEMENT;
}

void TACSCubicLineBasis::getVisPoint(int n, double pt[]) {
  static const double pts[4] = {-1.0, -0.5, 0.5, 1.0};
  pt[0] = pts[n];
}

int TACSCubicLineBasis::getNumNodes() { return 4; }

int TACSCubicLineBasis::getNumParameters() { return 1; }

int TACSCubicLineBasis::getNumQuadraturePoints() { return 4; }

double TACSCubicLineBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts4[n];
}

double TACSCubicLineBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n];
  return TacsGaussQuadWts4[n];
}

int TACSCubicLineBasis::getNumElementFaces() { return 2; }

int TACSCubicLineBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSCubicLineBasis::getFaceQuadraturePoint(int face, int n, double pt[],
                                                  double t[]) {
  if (face == 0) {
    t[0] = -1.0;
    pt[0] = -1.0;
  } else {
    t[0] = 1.0;
    pt[1] = 1.0;
  }

  return 1.0;
}

void TACSCubicLineBasis::computeBasis(const double pt[], double N[]) {
  N[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  N[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);
}

void TACSCubicLineBasis::computeBasisGradient(const double pt[], double N[],
                                              double Nxi[]) {
  N[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  N[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  Nxi[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  Nxi[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  Nxi[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  Nxi[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;
}
