#ifndef TACS_SHELL_ELEMENT_QUADRATURE_H
#define TACS_SHELL_ELEMENT_QUADRATURE_H

#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSTriangleQuadrature.h"

/**
  Defines the quadrature over both the face and quadrature
*/
class TACSQuadLinearQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 4;

  static int getNumParameters() { return 2; }
  static int getNumQuadraturePoints() { return 4; }
  static double getQuadratureWeight(int n) {
    return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
  }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts2[n % 2];
    pt[1] = TacsGaussQuadPts2[n / 2];

    return TacsGaussQuadWts2[n % 2] * TacsGaussQuadWts2[n / 2];
  }
  static int getNumElementFaces() { return 4; }
  static int getNumFaceQuadraturePoints(int face) { return 2; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    if (face / 2 == 0) {
      pt[0] = -1.0 + 2.0 * (face % 2);
      pt[1] = TacsGaussQuadPts2[n];
    } else {
      pt[0] = TacsGaussQuadPts2[n];
      pt[1] = -1.0 + 2.0 * (face % 2);
    }

    if (face == 0) {
      // -X edge
      t[0] = 0.0;
      t[1] = -1.0;
    } else if (face == 1) {
      // +X edge
      t[0] = 0.0;
      t[1] = 1.0;
    } else if (face == 2) {
      // -Y edge
      t[0] = 1.0;
      t[1] = 0.0;
    } else if (face == 3) {
      // +Y edge
      t[0] = -1.0;
      t[1] = 0.0;
    }

    return TacsGaussQuadWts2[n];
  }
};

class TACSQuadQuadraticQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 9;

  static int getNumParameters() { return 2; }
  static int getNumQuadraturePoints() { return 9; }
  static double getQuadratureWeight(int n) {
    return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
  }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts3[n % 3];
    pt[1] = TacsGaussQuadPts3[n / 3];

    return TacsGaussQuadWts3[n % 3] * TacsGaussQuadWts3[n / 3];
  }
  static int getNumElementFaces() { return 4; }
  static int getNumFaceQuadraturePoints(int face) { return 3; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    if (face / 2 == 0) {
      pt[0] = -1.0 + 2.0 * (face % 2);
      pt[1] = TacsGaussQuadPts3[n];
    } else {
      pt[0] = TacsGaussQuadPts3[n];
      pt[1] = -1.0 + 2.0 * (face % 2);
    }

    if (face == 0) {
      // -X edge
      t[0] = 0.0;
      t[1] = -1.0;
    } else if (face == 1) {
      // +X edge
      t[0] = 0.0;
      t[1] = 1.0;
    } else if (face == 2) {
      // -Y edge
      t[0] = 1.0;
      t[1] = 0.0;
    } else if (face == 3) {
      // +Y edge
      t[0] = -1.0;
      t[1] = 0.0;
    }

    return TacsGaussQuadWts3[n];
  }
};

class TACSQuadCubicQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 16;

  static int getNumParameters() { return 2; }
  static int getNumQuadraturePoints() { return 16; }
  static double getQuadratureWeight(int n) {
    return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
  }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts4[n % 4];
    pt[1] = TacsGaussQuadPts4[n / 4];

    return TacsGaussQuadWts4[n % 4] * TacsGaussQuadWts4[n / 4];
  }
  static int getNumElementFaces() { return 4; }
  static int getNumFaceQuadraturePoints(int face) { return 4; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    if (face / 2 == 0) {
      pt[0] = -1.0 + 2.0 * (face % 2);
      pt[1] = TacsGaussQuadPts4[n];
    } else {
      pt[0] = TacsGaussQuadPts4[n];
      pt[1] = -1.0 + 2.0 * (face % 2);
    }

    if (face == 0) {
      // -X edge
      t[0] = 0.0;
      t[1] = -1.0;
    } else if (face == 1) {
      // +X edge
      t[0] = 0.0;
      t[1] = 1.0;
    } else if (face == 2) {
      // -Y edge
      t[0] = 1.0;
      t[1] = 0.0;
    } else if (face == 3) {
      // +Y edge
      t[0] = -1.0;
      t[1] = 0.0;
    }

    return TacsGaussQuadWts4[n];
  }
};

class TACSTriLinearQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 3;

  static int getNumParameters() { return 2; }
  static int getNumQuadraturePoints() { return 3; }
  static double getQuadratureWeight(int n) { return TacsTriangleWts3[n]; }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsTrianglePts3[2 * n];
    pt[1] = TacsTrianglePts3[2 * n + 1];

    return TacsTriangleWts3[n];
  }
  static int getNumElementFaces() { return 3; }
  static int getNumFaceQuadraturePoints(int face) { return 2; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
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
      return 0.5 * sqrt(2.0) * TacsGaussQuadWts2[n];
    }

    return 0.0;
  }
};

class TACSTriQuadraticQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 6;

  static int getNumParameters() { return 2; }
  static int getNumQuadraturePoints() { return 4; }
  static double getQuadratureWeight(int n) { return TacsTriangleWts4[n]; }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsTrianglePts4[2 * n];
    pt[1] = TacsTrianglePts4[2 * n + 1];

    return TacsTriangleWts4[n];
  }
  static int getNumElementFaces() { return 3; }
  static int getNumFaceQuadraturePoints(int face) { return 2; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
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
      return 0.5 * sqrt(2.0) * TacsGaussQuadWts2[n];
    }

    return 0.0;
  }
};

#endif  // TACS_SHELL_ELEMENT_QUADRATURE_H