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

#include "TACSHexaBasis.h"

#include "TACSGaussQuadrature.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSTensorProductBasisImpl.h"

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
  Linear Hexa basis class functions
*/
TACSLinearHexaBasis::TACSLinearHexaBasis() {
  for (int i = 0; i < 2; i++) {
    TacsLagrangeShapeFuncDerivative(2, TacsGaussQuadPts2[i],
                                    TacsGaussLobattoPoints2, &Nf[2 * i],
                                    &Nfx[2 * i]);
  }
}

ElementLayout TACSLinearHexaBasis::getLayoutType() { return TACS_HEXA_ELEMENT; }

void TACSLinearHexaBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 2);
  pt[1] = -1.0 + 2.0 * ((n % 4) / 2);
  pt[2] = -1.0 + 2.0 * (n / 4);
}

int TACSLinearHexaBasis::getNumNodes() { return 8; }

int TACSLinearHexaBasis::getNumParameters() { return 3; }

int TACSLinearHexaBasis::getNumQuadraturePoints() { return 8; }

double TACSLinearHexaBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[(n % 4) / 2] *
          TacsGaussQuadWts2[n / 4]);
}

double TACSLinearHexaBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts2[n % 2];
  pt[1] = TacsGaussQuadPts2[(n % 4) / 2];
  pt[2] = TacsGaussQuadPts2[n / 4];

  return (TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[(n % 4) / 2] *
          TacsGaussQuadWts2[n / 4]);
}

int TACSLinearHexaBasis::getNumElementFaces() { return 6; }

int TACSLinearHexaBasis::getNumFaceQuadraturePoints(int face) { return 4; }

double TACSLinearHexaBasis::getFaceQuadraturePoint(int face, int n, double pt[],
                                                   double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts2[n % 2];
    pt[2] = TacsGaussQuadPts2[n / 2];
  } else if (face / 2 == 1) {
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = -1.0 + 2.0 * (face % 2);
    pt[2] = TacsGaussQuadPts2[n / 2];
  } else {
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = TacsGaussQuadPts2[n / 2];
    pt[2] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2]);
}

void TACSLinearHexaBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
  N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
  N[2] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
  N[3] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
  N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
  N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
  N[6] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
  N[7] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
}

void TACSLinearHexaBasis::computeBasisGradient(const double pt[], double N[],
                                               double Nxi[]) {
  N[0] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
  N[1] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 - pt[2]);
  N[2] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
  N[3] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 - pt[2]);
  N[4] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
  N[5] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]) * (1.0 + pt[2]);
  N[6] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);
  N[7] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]) * (1.0 + pt[2]);

  Nxi[0] = -0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
  Nxi[1] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
  Nxi[2] = -0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);

  Nxi[3] = 0.125 * (1.0 - pt[1]) * (1.0 - pt[2]);
  Nxi[4] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
  Nxi[5] = -0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);

  Nxi[6] = -0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
  Nxi[7] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[2]);
  Nxi[8] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

  Nxi[9] = 0.125 * (1.0 + pt[1]) * (1.0 - pt[2]);
  Nxi[10] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[2]);
  Nxi[11] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);

  Nxi[12] = -0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
  Nxi[13] = -0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
  Nxi[14] = 0.125 * (1.0 - pt[0]) * (1.0 - pt[1]);

  Nxi[15] = 0.125 * (1.0 - pt[1]) * (1.0 + pt[2]);
  Nxi[16] = -0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
  Nxi[17] = 0.125 * (1.0 + pt[0]) * (1.0 - pt[1]);

  Nxi[18] = -0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
  Nxi[19] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[2]);
  Nxi[20] = 0.125 * (1.0 - pt[0]) * (1.0 + pt[1]);

  Nxi[21] = 0.125 * (1.0 + pt[1]) * (1.0 + pt[2]);
  Nxi[22] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[2]);
  Nxi[23] = 0.125 * (1.0 + pt[0]) * (1.0 + pt[1]);
}

void TACSLinearHexaBasis::interpFields(const int n, const double pt[],
                                       const int m, const TacsScalar v[],
                                       const int incr, TacsScalar u[]) {
  double n1[2], n2[2], n3[2];
  n1[0] = 0.5 * (1.0 - pt[0]);
  n1[1] = 0.5 * (1.0 + pt[0]);
  n2[0] = 0.5 * (1.0 - pt[1]);
  n2[1] = 0.5 * (1.0 + pt[1]);
  n3[0] = 0.5 * (1.0 - pt[2]);
  n3[1] = 0.5 * (1.0 + pt[2]);

  for (int p = 0; p < m; p++, v++, u += incr) {
    TacsInterpTensor3DInterp2(m, n1, n2, n3, v, &u[0]);
  }
}

void TACSLinearHexaBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  double n1[2], n2[2], n3[2];
  n1[0] = 0.5 * (1.0 - pt[0]);
  n1[1] = 0.5 * (1.0 + pt[0]);
  n2[0] = 0.5 * (1.0 - pt[1]);
  n2[1] = 0.5 * (1.0 + pt[1]);
  n3[0] = 0.5 * (1.0 - pt[2]);
  n3[1] = 0.5 * (1.0 + pt[2]);

  for (int p = 0; p < m; p++, v++, u += incr) {
    TacsAddTransTensor3DInterp2(m, n1, n2, n3, u, v);
  }
}

void TACSLinearHexaBasis::interpFieldsGrad(const int n, const double pt[],
                                           const int m, const TacsScalar v[],
                                           TacsScalar g[]) {
  double n1[2], n2[2], n3[2], n1x[2], n2x[2], n3x[3];
  n1[0] = 0.5 * (1.0 - pt[0]);
  n1[1] = 0.5 * (1.0 + pt[0]);
  n2[0] = 0.5 * (1.0 - pt[1]);
  n2[1] = 0.5 * (1.0 + pt[1]);
  n3[0] = 0.5 * (1.0 - pt[2]);
  n3[1] = 0.5 * (1.0 + pt[2]);
  n1x[0] = -0.5;
  n1x[1] = 0.5;
  n2x[0] = -0.5;
  n2x[1] = 0.5;
  n3x[0] = -0.5;
  n3x[1] = 0.5;

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsGradTensor3DInterp2(m, n1, n2, n3, n1x, n2x, n3x, v, g);
  }
}

void TACSLinearHexaBasis::addInterpFieldsGradTranspose(int n, const double pt[],
                                                       const int m,
                                                       const TacsScalar g[],
                                                       TacsScalar v[]) {
  double n1[2], n2[2], n3[2], n1x[2], n2x[2], n3x[3];
  n1[0] = 0.5 * (1.0 - pt[0]);
  n1[1] = 0.5 * (1.0 + pt[0]);
  n2[0] = 0.5 * (1.0 - pt[1]);
  n2[1] = 0.5 * (1.0 + pt[1]);
  n3[0] = 0.5 * (1.0 - pt[2]);
  n3[1] = 0.5 * (1.0 + pt[2]);
  n1x[0] = -0.5;
  n1x[1] = 0.5;
  n2x[0] = -0.5;
  n2x[1] = 0.5;
  n3x[0] = -0.5;
  n3x[1] = 0.5;

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsAddGradTransTensor3DInterp2(m, n1, n2, n3, n1x, n2x, n3x, g, v);
  }
}

void TACSLinearHexaBasis::interpAllFieldsGrad(const int m, const TacsScalar *v,
                                              TacsScalar *out) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TACSInterpAllTensor3DInterp2(1, 0, Nf, Nfx, v, out);
  } else if (m == 3) {
    TACSInterpAllTensor3DInterp2(3, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp2(3, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp2(3, 2, Nf, Nfx, &v[2], out);
  } else if (m == 4) {
    TACSInterpAllTensor3DInterp2(4, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp2(4, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp2(4, 2, Nf, Nfx, &v[2], out);
    TACSInterpAllTensor3DInterp2(4, 3, Nf, Nfx, &v[3], out);
  } else {
    for (int i = 0; i < m; i++) {
      TACSInterpAllTensor3DInterp2(m, i, Nf, Nfx, &v[i], out);
    }
  }
}

void TACSLinearHexaBasis::addInterpAllFieldsGradTranspose(const int m,
                                                          const TacsScalar *in,
                                                          TacsScalar *v) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TacsAddAllTransTensor3DInterp2(1, 0, Nf, Nfx, in, v);
  } else if (m == 3) {
    TacsAddAllTransTensor3DInterp2(3, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp2(3, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp2(3, 2, Nf, Nfx, in, &v[2]);
  } else if (m == 4) {
    TacsAddAllTransTensor3DInterp2(4, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp2(4, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp2(4, 2, Nf, Nfx, in, &v[2]);
    TacsAddAllTransTensor3DInterp2(4, 3, Nf, Nfx, in, &v[3]);
  } else {
    for (int i = 0; i < m; i++) {
      TacsAddAllTransTensor3DInterp2(m, i, Nf, Nfx, in, &v[i]);
    }
  }
}

/*
  Quadratic Hexa basis class functions
*/
TACSQuadraticHexaBasis::TACSQuadraticHexaBasis() {
  for (int i = 0; i < 3; i++) {
    TacsLagrangeShapeFuncDerivative(3, TacsGaussQuadPts3[i],
                                    TacsGaussLobattoPoints3, &Nf[3 * i],
                                    &Nfx[3 * i]);
  }
}

ElementLayout TACSQuadraticHexaBasis::getLayoutType() {
  return TACS_HEXA_QUADRATIC_ELEMENT;
}

void TACSQuadraticHexaBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 3);
  pt[1] = -1.0 + 2.0 * ((n % 9) / 3);
  pt[2] = -1.0 + 2.0 * (n / 9);
}

int TACSQuadraticHexaBasis::getNumNodes() { return 27; }

int TACSQuadraticHexaBasis::getNumParameters() { return 3; }

int TACSQuadraticHexaBasis::getNumQuadraturePoints() { return 27; }

double TACSQuadraticHexaBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[(n % 9) / 3] *
          TacsGaussQuadWts3[n / 9]);
}

double TACSQuadraticHexaBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts3[n % 3];
  pt[1] = TacsGaussQuadPts3[(n % 9) / 3];
  pt[2] = TacsGaussQuadPts3[n / 9];

  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[(n % 9) / 3] *
          TacsGaussQuadWts3[n / 9]);
}

int TACSQuadraticHexaBasis::getNumElementFaces() { return 6; }

int TACSQuadraticHexaBasis::getNumFaceQuadraturePoints(int face) { return 9; }

double TACSQuadraticHexaBasis::getFaceQuadraturePoint(int face, int n,
                                                      double pt[], double t[]) {
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
    pt[2] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3]);
}

void TACSQuadraticHexaBasis::computeBasis(const double pt[], double N[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        N[0] = n1[i] * n2[j] * n3[k];
        N++;
      }
    }
  }
}

void TACSQuadraticHexaBasis::computeBasisGradient(const double pt[], double N[],
                                                  double Nxi[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  double n1x[3];
  n1x[0] = -0.5 + pt[0];
  n1x[1] = -2.0 * pt[0];
  n1x[2] = 0.5 + pt[0];

  double n2x[3];
  n2x[0] = -0.5 + pt[1];
  n2x[1] = -2.0 * pt[1];
  n2x[2] = 0.5 + pt[1];

  double n3x[3];
  n3x[0] = -0.5 + pt[2];
  n3x[1] = -2.0 * pt[2];
  n3x[2] = 0.5 + pt[2];

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        N[0] = n1[i] * n2[j] * n3[k];
        Nxi[0] = n1x[i] * n2[j] * n3[k];
        Nxi[1] = n1[i] * n2x[j] * n3[k];
        Nxi[2] = n1[i] * n2[j] * n3x[k];
        N++;
        Nxi += 3;
      }
    }
  }
}

void TACSQuadraticHexaBasis::interpFields(const int n, const double pt[],
                                          const int m, const TacsScalar v[],
                                          const int incr, TacsScalar u[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  for (int p = 0; p < m; p++, v++, u += incr) {
    TacsInterpTensor3DInterp3(m, n1, n2, n3, v, u);
  }
}

void TACSQuadraticHexaBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  for (int i = 0; i < m; i++, u += incr, v++) {
    TacsAddTransTensor3DInterp3(m, n1, n2, n3, u, v);
  }
}

void TACSQuadraticHexaBasis::interpFieldsGrad(const int n, const double pt[],
                                              const int m, const TacsScalar v[],
                                              TacsScalar g[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  double n1x[3];
  n1x[0] = -0.5 + pt[0];
  n1x[1] = -2.0 * pt[0];
  n1x[2] = 0.5 + pt[0];

  double n2x[3];
  n2x[0] = -0.5 + pt[1];
  n2x[1] = -2.0 * pt[1];
  n2x[2] = 0.5 + pt[1];

  double n3x[3];
  n3x[0] = -0.5 + pt[2];
  n3x[1] = -2.0 * pt[2];
  n3x[2] = 0.5 + pt[2];

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsGradTensor3DInterp3(m, n1, n2, n3, n1x, n2x, n3x, v, g);
  }
}

void TACSQuadraticHexaBasis::addInterpFieldsGradTranspose(int n,
                                                          const double pt[],
                                                          const int m,
                                                          const TacsScalar g[],
                                                          TacsScalar v[]) {
  double n1[3];
  n1[0] = -0.5 * pt[0] * (1.0 - pt[0]);
  n1[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
  n1[2] = 0.5 * (1.0 + pt[0]) * pt[0];

  double n2[3];
  n2[0] = -0.5 * pt[1] * (1.0 - pt[1]);
  n2[1] = (1.0 - pt[1]) * (1.0 + pt[1]);
  n2[2] = 0.5 * (1.0 + pt[1]) * pt[1];

  double n3[3];
  n3[0] = -0.5 * pt[2] * (1.0 - pt[2]);
  n3[1] = (1.0 - pt[2]) * (1.0 + pt[2]);
  n3[2] = 0.5 * (1.0 + pt[2]) * pt[2];

  double n1x[3];
  n1x[0] = -0.5 + pt[0];
  n1x[1] = -2.0 * pt[0];
  n1x[2] = 0.5 + pt[0];

  double n2x[3];
  n2x[0] = -0.5 + pt[1];
  n2x[1] = -2.0 * pt[1];
  n2x[2] = 0.5 + pt[1];

  double n3x[3];
  n3x[0] = -0.5 + pt[2];
  n3x[1] = -2.0 * pt[2];
  n3x[2] = 0.5 + pt[2];

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsAddGradTransTensor3DInterp3(m, n1, n2, n3, n1x, n2x, n3x, g, v);
  }
}

void TACSQuadraticHexaBasis::interpAllFieldsGrad(const int m,
                                                 const TacsScalar *v,
                                                 TacsScalar *out) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TACSInterpAllTensor3DInterp3(1, 0, Nf, Nfx, v, out);
  } else if (m == 3) {
    TACSInterpAllTensor3DInterp3(3, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp3(3, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp3(3, 2, Nf, Nfx, &v[2], out);
  } else if (m == 4) {
    TACSInterpAllTensor3DInterp3(4, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp3(4, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp3(4, 2, Nf, Nfx, &v[2], out);
    TACSInterpAllTensor3DInterp3(4, 3, Nf, Nfx, &v[3], out);
  } else {
    for (int i = 0; i < m; i++) {
      TACSInterpAllTensor3DInterp3(m, i, Nf, Nfx, &v[i], out);
    }
  }
}

void TACSQuadraticHexaBasis::addInterpAllFieldsGradTranspose(
    const int m, const TacsScalar *in, TacsScalar *v) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TacsAddAllTransTensor3DInterp3(1, 0, Nf, Nfx, in, v);
  } else if (m == 3) {
    TacsAddAllTransTensor3DInterp3(3, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp3(3, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp3(3, 2, Nf, Nfx, in, &v[2]);
  } else if (m == 4) {
    TacsAddAllTransTensor3DInterp3(4, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp3(4, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp3(4, 2, Nf, Nfx, in, &v[2]);
    TacsAddAllTransTensor3DInterp3(4, 3, Nf, Nfx, in, &v[3]);
  } else {
    for (int i = 0; i < m; i++) {
      TacsAddAllTransTensor3DInterp3(m, i, Nf, Nfx, in, &v[i]);
    }
  }
}

/*
  Cubic Hexa basis class functions
*/
TACSCubicHexaBasis::TACSCubicHexaBasis() {
  for (int i = 0; i < 4; i++) {
    TacsLagrangeShapeFuncDerivative(4, TacsGaussQuadPts4[i],
                                    TacsGaussLobattoPoints4, &Nf[4 * i],
                                    &Nfx[4 * i]);
  }
}

ElementLayout TACSCubicHexaBasis::getLayoutType() {
  return TACS_HEXA_CUBIC_ELEMENT;
}

void TACSCubicHexaBasis::getVisPoint(int n, double pt[]) {
  pt[0] = -1.0 + 2.0 * (n % 4);
  pt[1] = -1.0 + 2.0 * ((n % 16) / 4);
  pt[2] = -1.0 + 2.0 * (n / 16);
}

int TACSCubicHexaBasis::getNumNodes() { return 64; }

int TACSCubicHexaBasis::getNumParameters() { return 3; }

int TACSCubicHexaBasis::getNumQuadraturePoints() { return 64; }

double TACSCubicHexaBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[(n % 16) / 4] *
          TacsGaussQuadWts4[n / 16]);
}

double TACSCubicHexaBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts4[n % 4];
  pt[1] = TacsGaussQuadPts4[(n % 16) / 4];
  pt[2] = TacsGaussQuadPts4[n / 16];

  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[(n % 16) / 4] *
          TacsGaussQuadWts4[n / 16]);
}

int TACSCubicHexaBasis::getNumElementFaces() { return 6; }

int TACSCubicHexaBasis::getNumFaceQuadraturePoints(int face) { return 16; }

double TACSCubicHexaBasis::getFaceQuadraturePoint(int face, int n, double pt[],
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
    pt[2] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4]);
}

void TACSCubicHexaBasis::computeBasis(const double pt[], double N[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        N[0] = n1[i] * n2[j] * n3[k];
        N++;
      }
    }
  }
}

void TACSCubicHexaBasis::computeBasisGradient(const double pt[], double N[],
                                              double Nxi[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  double n1x[4];
  n1x[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  n1x[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  n1x[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  n1x[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double n2x[4];
  n2x[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  n2x[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  n2x[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  n2x[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  double n3x[4];
  n3x[0] = -2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] + 1.0 / 6.0;
  n3x[1] = 4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] - 4.0 / 3.0;
  n3x[2] = -4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] + 4.0 / 3.0;
  n3x[3] = 2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] - 1.0 / 6.0;

  for (int k = 0; k < 4; k++) {
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        N[0] = n1[i] * n2[j] * n3[k];
        Nxi[0] = n1x[i] * n2[j] * n3[k];
        Nxi[1] = n1[i] * n2x[j] * n3[k];
        Nxi[2] = n1[i] * n2[j] * n3x[k];
        N++;
        Nxi += 3;
      }
    }
  }
}

void TACSCubicHexaBasis::interpFields(const int n, const double pt[],
                                      const int m, const TacsScalar v[],
                                      const int incr, TacsScalar u[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  for (int p = 0; p < m; p++, v++, u += incr) {
    TacsInterpTensor3DInterp4(m, n1, n2, n3, v, u);
  }
}

void TACSCubicHexaBasis::addInterpFieldsTranspose(const int n,
                                                  const double pt[],
                                                  const int incr,
                                                  const TacsScalar u[],
                                                  const int m, TacsScalar v[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  for (int p = 0; p < m; p++, v++, u += incr) {
    TacsAddTransTensor3DInterp4(m, n1, n2, n3, u, v);
  }
}

void TACSCubicHexaBasis::interpFieldsGrad(const int n, const double pt[],
                                          const int m, const TacsScalar v[],
                                          TacsScalar g[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  double n1x[4];
  n1x[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  n1x[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  n1x[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  n1x[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double n2x[4];
  n2x[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  n2x[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  n2x[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  n2x[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  double n3x[4];
  n3x[0] = -2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] + 1.0 / 6.0;
  n3x[1] = 4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] - 4.0 / 3.0;
  n3x[2] = -4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] + 4.0 / 3.0;
  n3x[3] = 2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] - 1.0 / 6.0;

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsGradTensor3DInterp4(m, n1, n2, n3, n1x, n2x, n3x, v, g);
  }
}

void TACSCubicHexaBasis::addInterpFieldsGradTranspose(int n, const double pt[],
                                                      const int m,
                                                      const TacsScalar g[],
                                                      TacsScalar v[]) {
  double n1[4];
  n1[0] = -(2.0 / 3.0) * (0.5 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[1] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 - pt[0]) * (1.0 - pt[0]);
  n1[2] = (4.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (1.0 - pt[0]);
  n1[3] = -(2.0 / 3.0) * (1.0 + pt[0]) * (0.5 + pt[0]) * (0.5 - pt[0]);

  double n2[4];
  n2[0] = -(2.0 / 3.0) * (0.5 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[1] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 - pt[1]) * (1.0 - pt[1]);
  n2[2] = (4.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (1.0 - pt[1]);
  n2[3] = -(2.0 / 3.0) * (1.0 + pt[1]) * (0.5 + pt[1]) * (0.5 - pt[1]);

  double n3[4];
  n3[0] = -(2.0 / 3.0) * (0.5 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[1] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 - pt[2]) * (1.0 - pt[2]);
  n3[2] = (4.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (1.0 - pt[2]);
  n3[3] = -(2.0 / 3.0) * (1.0 + pt[2]) * (0.5 + pt[2]) * (0.5 - pt[2]);

  double n1x[4];
  n1x[0] = -2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] + 1.0 / 6.0;
  n1x[1] = 4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] - 4.0 / 3.0;
  n1x[2] = -4.0 * pt[0] * pt[0] - (4.0 / 3.0) * pt[0] + 4.0 / 3.0;
  n1x[3] = 2.0 * pt[0] * pt[0] + (4.0 / 3.0) * pt[0] - 1.0 / 6.0;

  double n2x[4];
  n2x[0] = -2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] + 1.0 / 6.0;
  n2x[1] = 4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] - 4.0 / 3.0;
  n2x[2] = -4.0 * pt[1] * pt[1] - (4.0 / 3.0) * pt[1] + 4.0 / 3.0;
  n2x[3] = 2.0 * pt[1] * pt[1] + (4.0 / 3.0) * pt[1] - 1.0 / 6.0;

  double n3x[4];
  n3x[0] = -2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] + 1.0 / 6.0;
  n3x[1] = 4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] - 4.0 / 3.0;
  n3x[2] = -4.0 * pt[2] * pt[2] - (4.0 / 3.0) * pt[2] + 4.0 / 3.0;
  n3x[3] = 2.0 * pt[2] * pt[2] + (4.0 / 3.0) * pt[2] - 1.0 / 6.0;

  for (int p = 0; p < m; p++, v++, g += 3) {
    TacsAddGradTransTensor3DInterp4(m, n1, n2, n3, n1x, n2x, n3x, g, v);
  }
}

void TACSCubicHexaBasis::interpAllFieldsGrad(const int m, const TacsScalar v[],
                                             TacsScalar out[]) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TACSInterpAllTensor3DInterp4(1, 0, Nf, Nfx, v, out);
  } else if (m == 3) {
    TACSInterpAllTensor3DInterp4(3, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp4(3, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp4(3, 2, Nf, Nfx, &v[2], out);
  } else if (m == 4) {
    TACSInterpAllTensor3DInterp4(4, 0, Nf, Nfx, &v[0], out);
    TACSInterpAllTensor3DInterp4(4, 1, Nf, Nfx, &v[1], out);
    TACSInterpAllTensor3DInterp4(4, 2, Nf, Nfx, &v[2], out);
    TACSInterpAllTensor3DInterp4(4, 3, Nf, Nfx, &v[3], out);
  } else {
    for (int i = 0; i < m; i++) {
      TACSInterpAllTensor3DInterp4(m, i, Nf, Nfx, &v[i], out);
    }
  }
}

void TACSCubicHexaBasis::addInterpAllFieldsGradTranspose(const int m,
                                                         const TacsScalar in[],
                                                         TacsScalar v[]) {
  // Try to force the compiler to do inline optimization
  if (m == 1) {
    TacsAddAllTransTensor3DInterp4(1, 0, Nf, Nfx, in, v);
  } else if (m == 3) {
    TacsAddAllTransTensor3DInterp4(3, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp4(3, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp4(3, 2, Nf, Nfx, in, &v[2]);
  } else if (m == 4) {
    TacsAddAllTransTensor3DInterp4(4, 0, Nf, Nfx, in, &v[0]);
    TacsAddAllTransTensor3DInterp4(4, 1, Nf, Nfx, in, &v[1]);
    TacsAddAllTransTensor3DInterp4(4, 2, Nf, Nfx, in, &v[2]);
    TacsAddAllTransTensor3DInterp4(4, 3, Nf, Nfx, in, &v[3]);
  } else {
    for (int i = 0; i < m; i++) {
      TacsAddAllTransTensor3DInterp4(m, i, Nf, Nfx, in, &v[i]);
    }
  }
}

/*
  Quartic hexahedral basis class functions
*/
TACSQuarticHexaBasis::TACSQuarticHexaBasis() {
  for (int i = 0; i < 5; i++) {
    TacsLagrangeShapeFuncDerivative(5, TacsGaussQuadPts5[i], cosine_pts,
                                    &Nf[5 * i], &Nfxi[5 * i]);
  }
}

const double TACSQuarticHexaBasis::cosine_pts[5] = {
    -1.0, -0.7071067811865475, 0.0, 0.7071067811865475, 1.0};

ElementLayout TACSQuarticHexaBasis::getLayoutType() {
  return TACS_HEXA_QUARTIC_ELEMENT;
}

void TACSQuarticHexaBasis::getVisPoint(int n, double pt[]) {
  pt[0] = cosine_pts[(n % 25) % 5];
  pt[1] = cosine_pts[(n / 25) % 5];
  pt[2] = cosine_pts[n / 25];
}

int TACSQuarticHexaBasis::getNumNodes() { return 125; }

int TACSQuarticHexaBasis::getNumParameters() { return 3; }

int TACSQuarticHexaBasis::getNumQuadraturePoints() { return 125; }

double TACSQuarticHexaBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[(n % 25) / 5] *
          TacsGaussQuadWts5[n / 25]);
}

double TACSQuarticHexaBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts5[n % 5];
  pt[1] = TacsGaussQuadPts5[(n % 25) / 5];
  pt[2] = TacsGaussQuadPts5[n / 25];

  return (TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[(n % 25) / 5] *
          TacsGaussQuadWts5[n / 25]);
}

int TACSQuarticHexaBasis::getNumElementFaces() { return 6; }

int TACSQuarticHexaBasis::getNumFaceQuadraturePoints(int face) { return 25; }

double TACSQuarticHexaBasis::getFaceQuadraturePoint(int face, int n,
                                                    double pt[], double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts5[n % 5];
    pt[2] = TacsGaussQuadPts5[n / 5];
  } else if (face / 2 == 1) {
    pt[0] = TacsGaussQuadPts5[n % 5];
    pt[1] = -1.0 + 2.0 * (face % 2);
    pt[2] = TacsGaussQuadPts5[n / 5];
  } else {
    pt[0] = TacsGaussQuadPts5[n % 5];
    pt[1] = TacsGaussQuadPts5[n / 5];
    pt[2] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts5[n % 5] * TacsGaussQuadWts5[n / 5]);
}

void TACSQuarticHexaBasis::computeBasis(const double pt[], double N[]) {
  double na[5], nb[5], nc[5];
  TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, nb);
  TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, nc);

  for (int k = 0; k < 5; k++) {
    for (int j = 0; j < 5; j++) {
      for (int i = 0; i < 5; i++) {
        N[i + 5 * j + 25 * k] = na[i] * nb[j] * nc[k];
      }
    }
  }
}

void TACSQuarticHexaBasis::computeBasisGradient(const double pt[], double N[],
                                                double Nxi[]) {
  double na[5], nb[5], nc[5];
  double dna[5], dnb[5], dnc[5];
  TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, nb, dnb);
  TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, nc, dnc);

  for (int k = 0; k < 5; k++) {
    for (int j = 0; j < 5; j++) {
      for (int i = 0; i < 5; i++) {
        N[i + 5 * j + 25 * k] = na[i] * nb[j] * nc[k];
        Nxi[3 * (i + 5 * j + 25 * k)] = dna[i] * nb[j] * nc[k];
        Nxi[3 * (i + 5 * j + 25 * k) + 1] = na[i] * dnb[j] * nc[k];
        Nxi[3 * (i + 5 * j + 25 * k) + 2] = na[i] * nb[j] * dnc[k];
      }
    }
  }
}

void TACSQuarticHexaBasis::interpFields(const int n, const double pt[],
                                        const int m, const TacsScalar v[],
                                        const int incr, TacsScalar u[]) {
  for (int p = 0; p < m; p++) {
    u[p * incr] = 0.0;
  }

  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * ((n % 25) / 5)];
    const double *n3 = &Nf[5 * (n / 25)];

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
                           n1[3] * v[3 * m] + n1[4] * v[4 * m]);
          u[p * incr] += n23 * t1;
          v++;
        }
        v += 4 * m;
      }
    }
  } else {
    double n1[5], n2[5], n3[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, n3);

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
                           n1[3] * v[3 * m] + n1[4] * v[4 * m]);
          u[p * incr] += n23 * t1;
          v++;
        }
        v += 4 * m;
      }
    }
  }
}

void TACSQuarticHexaBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * ((n % 25) / 5)];
    const double *n3 = &Nf[5 * (n / 25)];

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * u[p * incr];
          v[0] += a * n1[0];
          v[m] += a * n1[1];
          v[2 * m] += a * n1[2];
          v[3 * m] += a * n1[3];
          v[4 * m] += a * n1[4];
          v++;
        }
        v += 4 * m;
      }
    }
  } else {
    double n1[5], n2[5], n3[5];
    TacsLagrangeShapeFunctions(5, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(5, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(5, pt[2], cosine_pts, n3);

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * u[p * incr];
          v[0] += a * n1[0];
          v[m] += a * n1[1];
          v[2 * m] += a * n1[2];
          v[3 * m] += a * n1[3];
          v[4 * m] += a * n1[4];
          v++;
        }
        v += 4 * m;
      }
    }
  }
}

void TACSQuarticHexaBasis::interpFieldsGrad(const int n, const double pt[],
                                            const int m, const TacsScalar v[],
                                            TacsScalar grad[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * ((n % 25) / 5)];
    const double *n3 = &Nf[5 * (n / 25)];
    const double *n1x = &Nfxi[5 * (n % 5)];
    const double *n2x = &Nfxi[5 * ((n % 25) / 5)];
    const double *n3x = &Nfxi[5 * (n / 25)];

    memset(grad, 0, 3 * m * sizeof(TacsScalar));

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          g[0] += n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                         n1x[3] * v[3 * m] + n1x[4] * v[4 * m]);
          TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
                           n1[3] * v[3 * m] + n1[4] * v[4 * m]);
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }
  } else {
    double n1[5], n2[5], n3[5], n1x[5], n2x[5], n3x[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, n3, n3x);

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          g[0] += n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                         n1x[3] * v[3 * m] + n1x[4] * v[4 * m]);
          TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
                           n1[3] * v[3 * m] + n1[4] * v[4 * m]);
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }
  }
}

void TACSQuarticHexaBasis::addInterpFieldsGradTranspose(int n,
                                                        const double pt[],
                                                        const int m,
                                                        const TacsScalar grad[],
                                                        TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[5 * (n % 5)];
    const double *n2 = &Nf[5 * ((n % 25) / 5)];
    const double *n3 = &Nf[5 * (n / 25)];
    const double *n1x = &Nfxi[5 * (n % 5)];
    const double *n2x = &Nfxi[5 * ((n % 25) / 5)];
    const double *n3x = &Nfxi[5 * (n / 25)];

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        const TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          TacsScalar a = n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }
  } else {
    double n1[5], n2[5], n3[5], n1x[5], n2x[5], n3x[5];
    TacsLagrangeShapeFuncDerivative(5, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(5, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(5, pt[2], cosine_pts, n3, n3x);

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        const TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          TacsScalar a = n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }
  }
}

void TACSQuarticHexaBasis::interpAllFieldsGrad(const int m,
                                               const TacsScalar values[],
                                               TacsScalar out[]) {
  if (m == 1) {
    TACSInterpAllTensor3DInterp5VarsPerNode1(Nf, Nfxi, values, out);
  } else if (m == 3) {
    TACSInterpAllTensor3DInterp5VarsPerNode3(Nf, Nfxi, values, out);
  } else if (m == 4) {
    TACSInterpAllTensor3DInterp5VarsPerNode4(Nf, Nfxi, values, out);
  } else {
    TACSInterpAllTensor3DInterp5(m, Nf, Nfxi, values, out);
  }
}

void TACSQuarticHexaBasis::addInterpAllFieldsGradTranspose(
    const int m, const TacsScalar in[], TacsScalar values[]) {
  if (m == 1) {
    TacsAddAllTransTensor3DInterp5VarsPerNode1(Nf, Nfxi, in, values);
  } else if (m == 3) {
    TacsAddAllTransTensor3DInterp5VarsPerNode3(Nf, Nfxi, in, values);
  } else if (m == 4) {
    TacsAddAllTransTensor3DInterp5VarsPerNode4(Nf, Nfxi, in, values);
  } else {
    TacsAddAllTransTensor3DInterp5(m, Nf, Nfxi, in, values);
  }
}

/*
  Quintic Quad basis class functions
*/
TACSQuinticHexaBasis::TACSQuinticHexaBasis() {
  for (int i = 0; i < 6; i++) {
    TacsLagrangeShapeFuncDerivative(6, TacsGaussQuadPts6[i], cosine_pts,
                                    &Nf[6 * i], &Nfxi[6 * i]);
  }
}

const double TACSQuinticHexaBasis::cosine_pts[6] = {-1.0,
                                                    -0.8090169943749475,
                                                    -0.30901699437494745,
                                                    0.30901699437494745,
                                                    0.8090169943749475,
                                                    1.0};

ElementLayout TACSQuinticHexaBasis::getLayoutType() {
  return TACS_HEXA_QUINTIC_ELEMENT;
}

void TACSQuinticHexaBasis::getVisPoint(int n, double pt[]) {
  pt[0] = cosine_pts[n % 6];
  pt[1] = cosine_pts[n / 6];
}

int TACSQuinticHexaBasis::getNumNodes() { return 216; }

int TACSQuinticHexaBasis::getNumParameters() { return 3; }

int TACSQuinticHexaBasis::getNumQuadraturePoints() { return 216; }

double TACSQuinticHexaBasis::getQuadratureWeight(int n) {
  return (TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[(n % 36) / 6] *
          TacsGaussQuadWts6[n / 36]);
}

double TACSQuinticHexaBasis::getQuadraturePoint(int n, double pt[]) {
  pt[0] = TacsGaussQuadPts6[n % 6];
  pt[1] = TacsGaussQuadPts6[(n % 36) / 6];
  pt[2] = TacsGaussQuadPts6[n / 36];

  return (TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[(n % 36) / 6] *
          TacsGaussQuadWts6[n / 36]);
}

int TACSQuinticHexaBasis::getNumElementFaces() { return 6; }

int TACSQuinticHexaBasis::getNumFaceQuadraturePoints(int face) { return 36; }

double TACSQuinticHexaBasis::getFaceQuadraturePoint(int face, int n,
                                                    double pt[], double t[]) {
  if (face / 2 == 0) {
    pt[0] = -1.0 + 2.0 * (face % 2);
    pt[1] = TacsGaussQuadPts6[n % 6];
    pt[2] = TacsGaussQuadPts6[n / 6];
  } else if (face / 2 == 1) {
    pt[0] = TacsGaussQuadPts6[n % 6];
    pt[1] = -1.0 + 2.0 * (face % 2);
    pt[2] = TacsGaussQuadPts6[n / 6];
  } else {
    pt[0] = TacsGaussQuadPts6[n % 6];
    pt[1] = TacsGaussQuadPts6[n / 6];
    pt[2] = -1.0 + 2.0 * (face % 2);
  }

  getFaceTangents(face, t);

  return (TacsGaussQuadWts6[n % 6] * TacsGaussQuadWts6[n / 6]);
}

void TACSQuinticHexaBasis::computeBasis(const double pt[], double N[]) {
  double na[6], nb[6], nc[6];
  TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, na);
  TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, nb);
  TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, nc);

  for (int k = 0; k < 6; k++) {
    for (int j = 0; j < 6; j++) {
      for (int i = 0; i < 6; i++) {
        N[i + 6 * j + 36 * k] = na[i] * nb[j] * nc[k];
      }
    }
  }
}

void TACSQuinticHexaBasis::computeBasisGradient(const double pt[], double N[],
                                                double Nxi[]) {
  double na[6], nb[6], nc[6];
  double dna[6], dnb[6], dnc[6];
  TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, na, dna);
  TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, nb, dnb);
  TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, nc, dnc);

  for (int k = 0; k < 6; k++) {
    for (int j = 0; j < 6; j++) {
      for (int i = 0; i < 6; i++) {
        N[i + 6 * j + 36 * k] = na[i] * nb[j] * nc[k];
        Nxi[3 * (i + 6 * j + 36 * k)] = dna[i] * nb[j] * nc[k];
        Nxi[3 * (i + 6 * j + 36 * k) + 1] = na[i] * dnb[j] * nc[k];
        Nxi[3 * (i + 6 * j + 36 * k) + 2] = na[i] * nb[j] * dnc[k];
      }
    }
  }
}

void TACSQuinticHexaBasis::interpFields(const int n, const double pt[],
                                        const int m, const TacsScalar v[],
                                        const int incr, TacsScalar u[]) {
  for (int p = 0; p < m; p++) {
    u[p * incr] = 0.0;
  }

  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * ((n % 36) / 6)];
    const double *n3 = &Nf[6 * (n / 36)];

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar t1 =
              (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
               n1[3] * v[3 * m] + n1[4] * v[4 * m] + n1[5] * v[5 * m]);
          u[p * incr] += n23 * t1;
          v++;
        }
        v += 5 * m;
      }
    }
  } else {
    double n1[6], n2[6], n3[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, n3);

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar t1 =
              (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
               n1[3] * v[3 * m] + n1[4] * v[4 * m] + n1[5] * v[5 * m]);
          u[p * incr] += n23 * t1;
          v++;
        }
        v += 5 * m;
      }
    }
  }
}

void TACSQuinticHexaBasis::addInterpFieldsTranspose(
    const int n, const double pt[], const int incr, const TacsScalar u[],
    const int m, TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * ((n % 36) / 6)];
    const double *n3 = &Nf[6 * (n / 36)];

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * u[p * incr];
          v[0] += a * n1[0];
          v[m] += a * n1[1];
          v[2 * m] += a * n1[2];
          v[3 * m] += a * n1[3];
          v[4 * m] += a * n1[4];
          v[5 * m] += a * n1[5];
          v++;
        }
        v += 5 * m;
      }
    }
  } else {
    double n1[6], n2[6], n3[6];
    TacsLagrangeShapeFunctions(6, pt[0], cosine_pts, n1);
    TacsLagrangeShapeFunctions(6, pt[1], cosine_pts, n2);
    TacsLagrangeShapeFunctions(6, pt[2], cosine_pts, n3);

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * u[p * incr];
          v[0] += a * n1[0];
          v[m] += a * n1[1];
          v[2 * m] += a * n1[2];
          v[3 * m] += a * n1[3];
          v[4 * m] += a * n1[4];
          v[5 * m] += a * n1[5];
          v++;
        }
        v += 5 * m;
      }
    }
  }
}

void TACSQuinticHexaBasis::interpFieldsGrad(const int n, const double pt[],
                                            const int m, const TacsScalar v[],
                                            TacsScalar grad[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * ((n % 36) / 6)];
    const double *n3 = &Nf[6 * (n / 36)];
    const double *n1x = &Nfxi[6 * (n % 6)];
    const double *n2x = &Nfxi[6 * ((n % 36) / 6)];
    const double *n3x = &Nfxi[6 * (n / 36)];

    memset(grad, 0, 3 * m * sizeof(TacsScalar));

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          g[0] +=
              n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                     n1x[3] * v[3 * m] + n1x[4] * v[4 * m] + n1x[5] * v[5 * m]);
          TacsScalar t1 =
              (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
               n1[3] * v[3 * m] + n1[4] * v[4 * m] + n1[5] * v[5 * m]);
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }
  } else {
    double n1[6], n2[6], n3[6], n1x[6], n2x[6], n3x[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, n3, n3x);

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          g[0] +=
              n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                     n1x[3] * v[3 * m] + n1x[4] * v[4 * m] + n1x[5] * v[5 * m]);
          TacsScalar t1 =
              (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
               n1[3] * v[3 * m] + n1[4] * v[4 * m] + n1[5] * v[5 * m]);
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }
  }
}

void TACSQuinticHexaBasis::addInterpFieldsGradTranspose(int n,
                                                        const double pt[],
                                                        const int m,
                                                        const TacsScalar grad[],
                                                        TacsScalar v[]) {
  if (n >= 0) {
    const double *n1 = &Nf[6 * (n % 6)];
    const double *n2 = &Nf[6 * ((n % 36) / 6)];
    const double *n3 = &Nf[6 * (n / 36)];
    const double *n1x = &Nfxi[6 * (n % 6)];
    const double *n2x = &Nfxi[6 * ((n % 36) / 6)];
    const double *n3x = &Nfxi[6 * (n / 36)];

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        const TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          TacsScalar a = n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          v[5 * m] += a * n1[5] + b * n1x[5];
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }
  } else {
    double n1[6], n2[6], n3[6], n1x[6], n2x[6], n3x[6];
    TacsLagrangeShapeFuncDerivative(6, pt[0], cosine_pts, n1, n1x);
    TacsLagrangeShapeFuncDerivative(6, pt[1], cosine_pts, n2, n2x);
    TacsLagrangeShapeFuncDerivative(6, pt[2], cosine_pts, n3, n3x);

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar n23 = n2[j] * n3[k];
        TacsScalar n2x3 = n2x[j] * n3[k];
        TacsScalar n23x = n2[j] * n3x[k];

        const TacsScalar *g = grad;
        for (int p = 0; p < m; p++) {
          TacsScalar a = n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          v[5 * m] += a * n1[5] + b * n1x[5];
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }
  }
}

void TACSQuinticHexaBasis::interpAllFieldsGrad(const int m,
                                               const TacsScalar values[],
                                               TacsScalar out[]) {
  if (m == 1) {
    TACSInterpAllTensor3DInterp6VarsPerNode1(Nf, Nfxi, values, out);
  } else if (m == 3) {
    TACSInterpAllTensor3DInterp6VarsPerNode3(Nf, Nfxi, values, out);
  } else if (m == 4) {
    TACSInterpAllTensor3DInterp6VarsPerNode4(Nf, Nfxi, values, out);
  } else {
    TACSInterpAllTensor3DInterp6(m, Nf, Nfxi, values, out);
  }
}

void TACSQuinticHexaBasis::addInterpAllFieldsGradTranspose(
    const int m, const TacsScalar in[], TacsScalar values[]) {
  if (m == 1) {
    TacsAddAllTransTensor3DInterp6VarsPerNode1(Nf, Nfxi, in, values);
  } else if (m == 3) {
    TacsAddAllTransTensor3DInterp6VarsPerNode3(Nf, Nfxi, in, values);
  } else if (m == 4) {
    TacsAddAllTransTensor3DInterp6VarsPerNode4(Nf, Nfxi, in, values);
  } else {
    TacsAddAllTransTensor3DInterp6(m, Nf, Nfxi, in, values);
  }
}
