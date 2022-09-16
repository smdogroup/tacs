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

#include "TACSBeamBasis.h"

#include "TACSGaussQuadrature.h"

/*
  Linear Quad basis class functions
*/
ElementLayout TACSLinearBeamBasis::getLayoutType() { return TACS_LINE_ELEMENT; }

void TACSLinearBeamBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * n;
}

int TACSLinearBeamBasis::getNumNodes() { return 2; }

int TACSLinearBeamBasis::getNumParameters() { return 1; }

int TACSLinearBeamBasis::getNumQuadraturePoints() { return 2; }

double TACSLinearBeamBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts2[n];
}

double TACSLinearBeamBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts2[n];
  return TacsGaussQuadWts2[n];
}

int TACSLinearBeamBasis::getNumElementFaces() { return 2; }

int TACSLinearBeamBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSLinearBeamBasis::getFaceQuadraturePoint(int face, int n, double pt[],
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

void TACSLinearBeamBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 0.5 * (1.0 - pt[0]);
  N[1] = 0.5 * (1.0 + pt[0]);
}

void TACSLinearBeamBasis::computeBasisGradient(const double pt[], double N[],
                                               double Nxi[]) {
  N[0] = 0.5 * (1.0 - pt[0]);
  N[1] = 0.5 * (1.0 + pt[0]);

  Nxi[0] = -0.5;
  Nxi[1] = 0.5;
}

/*
  Quadratic Quad basis class functions
*/
ElementLayout TACSQuadraticBeamBasis::getLayoutType() {
  return TACS_LINE_QUADRATIC_ELEMENT;
}

void TACSQuadraticBeamBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 1.0 * n;
}

int TACSQuadraticBeamBasis::getNumNodes() { return 3; }

int TACSQuadraticBeamBasis::getNumParameters() { return 1; }

int TACSQuadraticBeamBasis::getNumQuadraturePoints() { return 3; }

double TACSQuadraticBeamBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts3[n];
}

double TACSQuadraticBeamBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n];
  return TacsGaussQuadWts3[n];
}

int TACSQuadraticBeamBasis::getNumElementFaces() { return 2; }

int TACSQuadraticBeamBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSQuadraticBeamBasis::getFaceQuadraturePoint(int face, int n,
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

void TACSQuadraticBeamBasis::computeBasis(const double pt[], double N[]) {
  N[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  N[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  N[2] = 0.5 * (1.0 + pt[0]) * pt[0];
}

void TACSQuadraticBeamBasis::computeBasisGradient(const double pt[], double N[],
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
ElementLayout TACSCubicBeamBasis::getLayoutType() {
  return TACS_LINE_CUBIC_ELEMENT;
}

void TACSCubicBeamBasis::getVisPoint(int n, double pt[]) {
  static const double pts[4] = {-1.0, -0.5, 0.5, 1.0};
  pt[0] = pts[n];
}

int TACSCubicBeamBasis::getNumNodes() { return 4; }

int TACSCubicBeamBasis::getNumParameters() { return 1; }

int TACSCubicBeamBasis::getNumQuadraturePoints() { return 4; }

double TACSCubicBeamBasis::getQuadratureWeight(int n) {
  return TacsGaussQuadWts4[n];
}

double TACSCubicBeamBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n];
  return TacsGaussQuadWts4[n];
}

int TACSCubicBeamBasis::getNumElementFaces() { return 2; }

int TACSCubicBeamBasis::getNumFaceQuadraturePoints(int face) { return 1; }

double TACSCubicBeamBasis::getFaceQuadraturePoint(int face, int n, double pt[],
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

void TACSCubicBeamBasis::computeBasis(const double pt[], double N[]) {
  N[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  N[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  N[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);
}

void TACSCubicBeamBasis::computeBasisGradient(const double pt[], double N[],
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
