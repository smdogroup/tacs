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

#include "TACSTriangularBasis.h"

#include "TACSGaussQuadrature.h"
#include "TACSTriangleQuadrature.h"

static void getFaceTangents(int face, double t[]) {
  if (face == 0) {
    t[0] = 1.0;
    t[1] = 0.0;
  } else if (face == 1) {
    t[0] = -1.0;
    t[1] = 1.0;
  } else if (face == 2) {
    t[0] = 0.0;
    t[1] = 1.0;
  }
}

/*
  Linear Triangle basis class functions
*/
ElementLayout TACSLinearTriangleBasis::getLayoutType() {
  return TACS_TRI_ELEMENT;
}

void TACSLinearTriangleBasis::getVisPoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = pt[1] = 0.0;
  } else if (n == 1) {
    pt[0] = 1.0;
    pt[1] = 0.0;
  } else {
    pt[0] = 0.0;
    pt[1] = 1.0;
  }
}

int TACSLinearTriangleBasis::getNumNodes() { return 3; }

int TACSLinearTriangleBasis::getNumParameters() { return 2; }

int TACSLinearTriangleBasis::getNumQuadraturePoints() { return 3; }

double TACSLinearTriangleBasis::getQuadratureWeight(int n) {
  return TacsTriangleWts3[n];
}

double TACSLinearTriangleBasis::getQuadraturePoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = TacsTrianglePts3[0];
    pt[1] = TacsTrianglePts3[1];
  } else if (n == 1) {
    pt[0] = TacsTrianglePts3[2];
    pt[1] = TacsTrianglePts3[3];
  } else if (n == 2) {
    pt[0] = TacsTrianglePts3[4];
    pt[1] = TacsTrianglePts3[5];
  }

  return TacsTriangleWts3[n];
}

int TACSLinearTriangleBasis::getNumElementFaces() { return 3; }

int TACSLinearTriangleBasis::getNumFaceQuadraturePoints(int face) { return 2; }

double TACSLinearTriangleBasis::getFaceQuadraturePoint(int face, int n,
                                                       double pt[],
                                                       double t[]) {
  getFaceTangents(face, t);

  if (face == 0) {
    pt[0] = 0.5 * TacsGaussQuadPts2[n] + 0.5;
    pt[1] = 0.0;
    return 0.5 * TacsGaussQuadWts2[n];
  } else if (face == 2) {
    pt[0] = 0.0;
    pt[1] = 0.5 * TacsGaussQuadPts2[n] + 0.5;
    return 0.5 * TacsGaussQuadWts2[n];
  } else if (face == 1) {
    pt[0] = 1.0 - (0.5 * TacsGaussQuadPts2[n] + 0.5);
    pt[1] = 0.5 * TacsGaussQuadPts2[n] + 0.5;
    return sqrt(2) * 0.5 * TacsGaussQuadWts2[n];
  }

  return 0.0;
}

void TACSLinearTriangleBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 1 - pt[0] - pt[1];
  N[1] = pt[0];
  N[2] = pt[1];
}

void TACSLinearTriangleBasis::computeBasisGradient(const double pt[],
                                                   double N[], double Nxi[]) {
  N[0] = 1.0 - pt[0] - pt[1];
  N[1] = pt[0];
  N[2] = pt[1];

  Nxi[0] = -1.0;
  Nxi[1] = -1.0;
  Nxi[2] = 1.0;
  Nxi[3] = 0.0;
  Nxi[4] = 0.0;
  Nxi[5] = 1.0;
}

/*
  Quadratic Triangle basis class functions
*/
ElementLayout TACSQuadraticTriangleBasis::getLayoutType() {
  return TACS_TRI_QUADRATIC_ELEMENT;
}

void TACSQuadraticTriangleBasis::getVisPoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = 0.0;
    pt[1] = 0.0;
  } else if (n == 1) {
    pt[0] = 1.0;
    pt[1] = 0.0;
  } else if (n == 2) {
    pt[0] = 0.0;
    pt[1] = 1.0;
  } else if (n == 3) {
    pt[0] = 0.5;
    pt[1] = 0.0;
  } else if (n == 4) {
    pt[0] = 0.5;
    pt[1] = 0.5;
  } else if (n == 5) {
    pt[0] = 0.0;
    pt[1] = 0.5;
  }
}

int TACSQuadraticTriangleBasis::getNumNodes() { return 6; }

int TACSQuadraticTriangleBasis::getNumParameters() { return 2; }

int TACSQuadraticTriangleBasis::getNumQuadraturePoints() { return 4; }

double TACSQuadraticTriangleBasis::getQuadratureWeight(int n) {
  return TacsTriangleWts4[0];
}

double TACSQuadraticTriangleBasis::getQuadraturePoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = TacsTrianglePts4[0];
    pt[1] = TacsTrianglePts4[1];
    return TacsTriangleWts4[0];
  } else if (n == 1) {
    pt[0] = TacsTrianglePts4[2];
    pt[1] = TacsTrianglePts4[3];
    return TacsTriangleWts4[1];
  } else if (n == 2) {
    pt[0] = TacsTrianglePts4[4];
    pt[1] = TacsTrianglePts4[5];
    return TacsTriangleWts4[2];
  } else if (n == 3) {
    pt[0] = TacsTrianglePts4[6];
    pt[1] = TacsTrianglePts4[7];
    return TacsTriangleWts4[3];
  }

  return 0.0;
}

int TACSQuadraticTriangleBasis::getNumElementFaces() { return 3; }

int TACSQuadraticTriangleBasis::getNumFaceQuadraturePoints(int face) {
  return 3;
}

double TACSQuadraticTriangleBasis::getFaceQuadraturePoint(int face, int n,
                                                          double pt[],
                                                          double t[]) {
  getFaceTangents(face, t);

  if (face == 0) {
    pt[0] = 0.5 * TacsGaussQuadPts3[n] + 0.5;
    pt[1] = 0.0;
    return 0.5 * TacsGaussQuadWts3[n];
  } else if (face == 2) {
    pt[0] = 0.0;
    pt[1] = 0.5 * TacsGaussQuadPts3[n] + 0.5;
    return 0.5 * TacsGaussQuadWts3[n];
  } else if (face == 1) {
    pt[0] = 1.0 - (0.5 * TacsGaussQuadPts3[n] + 0.5);
    pt[1] = 0.5 * TacsGaussQuadPts3[n] + 0.5;
    return sqrt(2) * 0.5 * TacsGaussQuadWts3[n];
  }

  return 0.0;
}

void TACSQuadraticTriangleBasis::computeBasis(const double pt[], double N[]) {
  N[0] = (1.0 - pt[0] - pt[1]) * (1.0 - 2.0 * pt[0] - 2.0 * pt[1]);
  N[1] = pt[0] * (2.0 * pt[0] - 1.0);
  N[2] = pt[1] * (2.0 * pt[1] - 1.0);
  N[3] = 4.0 * pt[0] * (1.0 - pt[0] - pt[1]);
  N[4] = 4.0 * pt[0] * pt[1];
  N[5] = 4.0 * pt[1] * (1.0 - pt[0] - pt[1]);
}

void TACSQuadraticTriangleBasis::computeBasisGradient(const double pt[],
                                                      double N[],
                                                      double Nxi[]) {
  N[0] = (1.0 - pt[0] - pt[1]) * (1.0 - 2.0 * pt[0] - 2.0 * pt[1]);
  N[1] = pt[0] * (2.0 * pt[0] - 1.0);
  N[2] = pt[1] * (2.0 * pt[1] - 1.0);
  N[3] = 4.0 * pt[0] * (1.0 - pt[0] - pt[1]);
  N[4] = 4.0 * pt[0] * pt[1];
  N[5] = 4.0 * pt[1] * (1.0 - pt[0] - pt[1]);

  Nxi[0] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
  Nxi[1] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
  Nxi[2] = 4.0 * pt[0] - 1.0;
  Nxi[3] = 0.0;
  Nxi[4] = 0.0;
  Nxi[5] = 4.0 * pt[1] - 1.0;
  Nxi[6] = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
  Nxi[7] = -4.0 * pt[0];
  Nxi[8] = 4.0 * pt[1];
  Nxi[9] = 4.0 * pt[0];
  Nxi[10] = -4.0 * pt[1];
  Nxi[11] = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];
}

/*
  Cubic Triangle basis class functions
*/
ElementLayout TACSCubicTriangleBasis::getLayoutType() {
  return TACS_TRI_CUBIC_ELEMENT;
}

void TACSCubicTriangleBasis::getVisPoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = 0.0;
    pt[1] = 0.0;
  } else if (n == 1) {
    pt[0] = 1.0;
    pt[1] = 0.0;
  } else if (n == 2) {
    pt[0] = 0.0;
    pt[1] = 1.0;
  } else if (n == 3) {
    pt[0] = 0.25;
    pt[1] = 0.0;
  } else if (n == 4) {
    pt[0] = 0.0;
    pt[1] = 0.75;
  } else if (n == 5) {
    pt[0] = 0.75;
    pt[1] = 0.25;
  } else if (n == 6) {
    pt[0] = 0.25;
    pt[1] = 0.75;
  } else if (n == 7) {
    pt[0] = 0.0;
    pt[1] = 0.75;
  } else if (n == 8) {
    pt[0] = 0.0;
    pt[1] = 0.75;
  } else if (n == 9) {
    pt[0] = 0.25;
    pt[1] = 0.25;
  }
}

int TACSCubicTriangleBasis::getNumNodes() { return 10; }

int TACSCubicTriangleBasis::getNumParameters() { return 2; }

int TACSCubicTriangleBasis::getNumQuadraturePoints() { return 6; }

double TACSCubicTriangleBasis::getQuadratureWeight(int n) {
  return TacsTriangleWts6[n];
}

double TACSCubicTriangleBasis::getQuadraturePoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = TacsTrianglePts6[0];
    pt[1] = TacsTrianglePts6[1];
    return TacsTriangleWts6[0];
  } else if (n == 1) {
    pt[0] = TacsTrianglePts6[2];
    pt[1] = TacsTrianglePts6[3];
    return TacsTriangleWts6[0];
  } else if (n == 2) {
    pt[0] = TacsTrianglePts6[4];
    pt[1] = TacsTrianglePts6[5];
    return TacsTriangleWts6[0];
  } else if (n == 3) {
    pt[0] = TacsTrianglePts6[6];
    pt[1] = TacsTrianglePts6[7];
    return TacsTriangleWts6[1];
  } else if (n == 4) {
    pt[0] = TacsTrianglePts6[8];
    pt[1] = TacsTrianglePts6[9];
    return TacsTriangleWts6[1];
  } else if (n == 5) {
    pt[0] = TacsTrianglePts6[10];
    pt[1] = TacsTrianglePts6[11];
    return TacsTriangleWts6[1];
  }
  return 0.0;
}

int TACSCubicTriangleBasis::getNumElementFaces() { return 3; }

int TACSCubicTriangleBasis::getNumFaceQuadraturePoints(int face) { return 4; }

double TACSCubicTriangleBasis::getFaceQuadraturePoint(int face, int n,
                                                      double pt[], double t[]) {
  getFaceTangents(face, t);

  if (face == 0) {
    pt[0] = 0.5 * TacsGaussQuadPts4[n] + 0.5;
    pt[1] = 0.0;
    return 0.5 * TacsGaussQuadWts4[n];
  } else if (face == 2) {
    pt[0] = 0.0;
    pt[1] = 0.5 * TacsGaussQuadPts4[n] + 0.5;
    return 0.5 * TacsGaussQuadWts4[n];
  } else if (face == 1) {
    pt[0] = 1.0 - (0.5 * TacsGaussQuadPts4[n] + 0.5);
    pt[1] = 0.5 * TacsGaussQuadPts4[n] + 0.5;
    return sqrt(2) * 0.5 * TacsGaussQuadWts4[n];
  }

  return 0.0;
}

void TACSCubicTriangleBasis::computeBasis(const double pt[], double N[]) {
  // corner nodes
  N[0] = 0.5 * (3.0 * (1 - pt[0] - pt[1]) - 1) *
         (3.0 * (1 - pt[0] - pt[1]) - 2) * (1 - pt[0] - pt[1]);
  N[1] = 0.5 * (3.0 * pt[0] - 1) * (3.0 * pt[0] - 2) * pt[0];
  N[2] = 0.5 * (3.0 * pt[1] - 1) * (3.0 * pt[1] - 2) * pt[1];

  // mid-nodes
  N[3] =
      9.0 / 2.0 * pt[0] * (1 - pt[0] - pt[1]) * (3.0 * (1 - pt[0] - pt[1]) - 1);
  N[4] = 9.0 / 2.0 * pt[0] * (1 - pt[0] - pt[1]) * (3.0 * pt[0] - 1);
  N[5] = 9.0 / 2.0 * pt[0] * pt[1] * (3.0 * pt[0] - 1);
  N[6] = 9.0 / 2.0 * pt[0] * pt[1] * (3.0 * pt[1] - 1);
  N[7] = 9.0 / 2.0 * pt[1] * (1 - pt[0] - pt[1]) * (3.0 * pt[1] - 1);
  N[8] =
      9.0 / 2.0 * pt[1] * (1 - pt[0] - pt[1]) * (3.0 * (1 - pt[0] - pt[1]) - 1);

  // center node
  N[9] = 27.0 * pt[0] * pt[1] * (1 - pt[0] - pt[1]);
}

void TACSCubicTriangleBasis::computeBasisGradient(const double pt[], double N[],
                                                  double Nxi[]) {
  // corner nodes
  N[0] = 0.5 * (3.0 * (1 - pt[0] - pt[1]) - 1) *
         (3.0 * (1 - pt[0] - pt[1]) - 2) * (1 - pt[0] - pt[1]);
  N[1] = 0.5 * (3.0 * pt[0] - 1) * (3.0 * pt[0] - 2) * pt[0];
  N[2] = 0.5 * (3.0 * pt[1] - 1) * (3.0 * pt[1] - 2) * pt[1];

  // mid-nodes
  N[3] =
      9.0 / 2.0 * pt[0] * (1 - pt[0] - pt[1]) * (3.0 * (1 - pt[0] - pt[1]) - 1);
  N[4] = 9.0 / 2.0 * pt[0] * (1 - pt[0] - pt[1]) * (3.0 * pt[0] - 1);
  N[5] = 9.0 / 2.0 * pt[0] * pt[1] * (3.0 * pt[0] - 1);
  N[6] = 9.0 / 2.0 * pt[0] * pt[1] * (3.0 * pt[1] - 1);
  N[7] = 9.0 / 2.0 * pt[1] * (1 - pt[0] - pt[1]) * (3.0 * pt[1] - 1);
  N[8] =
      9.0 / 2.0 * pt[1] * (1 - pt[0] - pt[1]) * (3.0 * (1 - pt[0] - pt[1]) - 1);

  // center node
  N[9] = 27.0 * pt[0] * pt[1] * (1 - pt[0] - pt[1]);

  // gradient components
  Nxi[0] = -13.5 * pt[0] * pt[0] + pt[0] * (18.0 - 27.0 * pt[1]) -
           13.5 * pt[1] * pt[1] + 18.0 * pt[1] - 5.5;
  Nxi[1] = -13.5 * pt[0] * pt[0] + pt[0] * (18.0 - 27.0 * pt[1]) -
           13.5 * pt[1] * pt[1] + 18.0 * pt[1] - 5.5;
  Nxi[2] = 13.5 * pt[0] * pt[0] - 9.0 * pt[0] + 1.0;
  Nxi[3] = 0.0;
  Nxi[4] = 0.0;
  Nxi[5] = 13.5 * pt[1] * pt[1] - 9.0 * pt[1] + 1.0;
  Nxi[6] = 40.5 * pt[0] * pt[0] + pt[0] * (54.0 * pt[1] - 45.0) +
           13.5 * pt[1] * pt[1] - 22.5 * pt[1] + 9.0;
  Nxi[7] = pt[0] * (27.0 * pt[0] + 27.0 * pt[1] - 22.5);
  Nxi[8] =
      -40.5 * pt[0] * pt[0] + pt[0] * (36.0 - 27.0 * pt[1]) + 4.5 * pt[1] - 4.5;
  Nxi[9] = -13.5 * (pt[0] - 1.0 / 3.0) * pt[0];
  Nxi[10] = pt[1] * (27.0 * pt[0] - 4.5);
  Nxi[11] = 4.5 * pt[0] * (3 * pt[0] - 1.0);
  Nxi[12] = 4.5 * pt[1] * (3 * pt[1] - 1.0);
  Nxi[13] = pt[0] * (27.0 * pt[1] - 4.5);
  Nxi[14] = -13.5 * (pt[1] - 1.0 / 3.0) * pt[1];
  Nxi[15] =
      pt[0] * (4.5 - 27.0 * pt[1]) - 40.5 * pt[1] * pt[1] + 36.0 * pt[1] - 4.5;
  Nxi[16] = pt[1] * (27.0 * pt[0] + 27.0 * pt[0] - 22.5);
  Nxi[17] = 13.5 * pt[0] * pt[0] + pt[0] * (54.0 * pt[1] - 22.5) +
            40.5 * pt[1] * pt[1] - 45.0 * pt[1] + 9.0;
  Nxi[18] = -27.0 * pt[1] * (2 * pt[0] + pt[1] - 1);
  Nxi[19] = -27.0 * pt[0] * (2 * pt[1] + pt[0] - 1);
}
