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

#include "TACSHexaBernsteinBasis.h"

#include "TACSBernsteinInterpolation.h"
#include "TACSGaussQuadrature.h"

static void getFaceTangents(int face, double t[]) {
  if (face == 0) {
    // Z - Y plane make a -X normal direction
    t[0] = 0.0;
    t[1] = 0.0;
    t[2] = 1.0;
    t[3] = 0.0;
    t[4] = 1.0;
    t[5] = 0.0;
  } else if (face == 1) {
    // Y - Z plane makes a X normal direction
    t[0] = 0.0;
    t[1] = 1.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 0.0;
    t[5] = 1.0;
  } else if (face == 2) {
    // X - Z plane makes a -Y normal direction
    t[0] = 1.0;
    t[1] = 0.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 0.0;
    t[5] = 1.0;
  } else if (face == 3) {
    // Z - X plane makes a Y normal direction
    t[0] = 0.0;
    t[1] = 0.0;
    t[2] = 1.0;
    t[3] = 1.0;
    t[4] = 0.0;
    t[5] = 0.0;
  } else if (face == 4) {
    // Y - X plane makes a -Z normal direction
    t[0] = 0.0;
    t[1] = 1.0;
    t[2] = 0.0;
    t[3] = 1.0;
    t[4] = 0.0;
    t[5] = 0.0;
  } else if (face == 5) {
    // X - Y plane makes a Z normal direction
    t[0] = 1.0;
    t[1] = 0.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 0.0;
    t[5] = 1.0;
  }
}

/*
  Quadratic Hexa basis class functions
*/
ElementLayout TACSQuadraticHexaBernsteinBasis::getLayoutType() {
  return TACS_HEXA_QUADRATIC_ELEMENT;
}

void TACSQuadraticHexaBernsteinBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 3);
  pt[1] = -1.0 + 2.0 * ((n % 9) / 3);
  pt[2] = -1.0 + 2.0 * (n / 9);
}

int TACSQuadraticHexaBernsteinBasis::getNumNodes() { return 27; }

int TACSQuadraticHexaBernsteinBasis::getNumParameters() { return 3; }

int TACSQuadraticHexaBernsteinBasis::getNumQuadraturePoints() { return 27; }

double TACSQuadraticHexaBernsteinBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[(n % 9) / 3] *
          TacsGaussQuadWts3[n / 9]);
}

double TACSQuadraticHexaBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[(n % 9) / 3];
  pt[2] = TacsGaussQuadPts3[n / 9];

  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[(n % 9) / 3] *
          TacsGaussQuadWts3[n / 9]);
}

int TACSQuadraticHexaBernsteinBasis::getNumElementFaces() { return 6; }

int TACSQuadraticHexaBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 9;
}

double TACSQuadraticHexaBernsteinBasis::getFaceQuadraturePoint(int face, int n,
                                                               double pt[],
                                                               double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts3[n % 3];
    pt[2] = TacsGaussQuadPts3[n / 3];
  } else if (face / 2 == 1) {
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = -1.0 + 2.0 * (face % 2);
    pt[2] = TacsGaussQuadPts3[n / 3];
  } else {
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = TacsGaussQuadPts3[n / 3];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3]);
}

void TACSQuadraticHexaBernsteinBasis::computeBasis(const double pt[],
                                                   double N[]) {
  double na[3], nb[3], nc[3];
  TacsBernsteinShapeFunctions(3, pt[0], na);
  TacsBernsteinShapeFunctions(3, pt[1], nb);
  TacsBernsteinShapeFunctions(3, pt[2], nc);

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        N++;
      }
    }
  }
}

void TACSQuadraticHexaBernsteinBasis::computeBasisGradient(const double pt[],
                                                           double N[],
                                                           double Nxi[]) {
  double na[3], nb[3], nc[3];
  double dna[3], dnb[3], dnc[3];
  TacsBernsteinShapeFuncDerivative(3, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(3, pt[1], nb, dnb);
  TacsBernsteinShapeFuncDerivative(3, pt[2], nc, dnc);

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        Nxi[0] = dna[i] * nb[j] * nc[k];
        Nxi[1] = na[i] * dnb[j] * nc[k];
        Nxi[2] = na[i] * nb[j] * dnc[k];
        N++;
        Nxi += 3;
      }
    }
  }
}

/*
  Cubic Hexa basis class functions
*/
ElementLayout TACSCubicHexaBernsteinBasis::getLayoutType() {
  return TACS_HEXA_CUBIC_ELEMENT;
}

void TACSCubicHexaBernsteinBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 4);
  pt[1] = -1.0 + 2.0 * ((n % 16) / 4);
  pt[2] = -1.0 + 2.0 * (n / 16);
}

int TACSCubicHexaBernsteinBasis::getNumNodes() { return 64; }

int TACSCubicHexaBernsteinBasis::getNumParameters() { return 3; }

int TACSCubicHexaBernsteinBasis::getNumQuadraturePoints() { return 64; }

double TACSCubicHexaBernsteinBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[(n % 16) / 4] *
          TacsGaussQuadWts4[n / 16]);
}

double TACSCubicHexaBernsteinBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[(n % 16) / 4];
  pt[2] = TacsGaussQuadPts4[n / 16];

  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[(n % 16) / 4] *
          TacsGaussQuadWts4[n / 16]);
}

int TACSCubicHexaBernsteinBasis::getNumElementFaces() { return 4; }

int TACSCubicHexaBernsteinBasis::getNumFaceQuadraturePoints(int face) {
  return 4;
}

double TACSCubicHexaBernsteinBasis::getFaceQuadraturePoint(int face, int n,
                                                           double pt[],
                                                           double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts4[n % 4];
    pt[2] = TacsGaussQuadPts4[n / 4];
  } else if (face / 2 == 1) {
    pt[0] = TacsGaussQuadPts3[n % 4];
    pt[1] = -1.0 + 2.0 * (face % 2);
    pt[2] = TacsGaussQuadPts3[n / 4];
  } else {
    pt[0] = TacsGaussQuadPts3[n % 4];
    pt[1] = TacsGaussQuadPts3[n / 4];
    pt[1] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4]);
}

void TACSCubicHexaBernsteinBasis::computeBasis(const double pt[], double N[]) {
  double na[4], nb[4], nc[4];
  TacsBernsteinShapeFunctions(4, pt[0], na);
  TacsBernsteinShapeFunctions(4, pt[1], nb);
  TacsBernsteinShapeFunctions(4, pt[2], nc);

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        N++;
      }
    }
  }
}

void TACSCubicHexaBernsteinBasis::computeBasisGradient(const double pt[],
                                                       double N[],
                                                       double Nxi[]) {
  double na[4], nb[4], nc[4];
  double dna[4], dnb[4], dnc[4];
  TacsBernsteinShapeFuncDerivative(4, pt[0], na, dna);
  TacsBernsteinShapeFuncDerivative(4, pt[1], nb, dnb);
  TacsBernsteinShapeFuncDerivative(4, pt[2], nc, dnc);

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        N[0] = na[i] * nb[j] * nc[k];
        Nxi[0] = dna[i] * nb[j] * nc[k];
        Nxi[1] = na[i] * dnb[j] * nc[k];
        Nxi[2] = na[i] * nb[j] * dnc[k];
        N++;
        Nxi += 3;
      }
    }
  }
}
