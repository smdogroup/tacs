#ifndef TACS_BEAM_ELEMENT_QUADRATURE_H
#define TACS_BEAM_ELEMENT_QUADRATURE_H

#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"

/**
  Defines the quadrature for a linear beam
*/
class TACSBeamLinearQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 2;

  static int getNumParameters() { return 1; }
  static int getNumQuadraturePoints() { return 2; }
  static double getQuadratureWeight(int n) { return TacsGaussQuadWts2[n]; }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts2[n];

    return TacsGaussQuadWts2[n];
  }
  static int getNumElementFaces() { return 2; }
  static int getNumFaceQuadraturePoints(int face) { return 1; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    pt[0] = -1.0 + 2.0 * face;
    return 1.0;
  }
};

/**
  Defines the quadrature for a quadratic beam
*/
class TACSBeamQuadraticQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 3;

  static int getNumParameters() { return 1; }
  static int getNumQuadraturePoints() { return 3; }
  static double getQuadratureWeight(int n) { return TacsGaussQuadWts3[n]; }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts3[n];

    return TacsGaussQuadWts3[n];
  }
  static int getNumElementFaces() { return 2; }
  static int getNumFaceQuadraturePoints(int face) { return 1; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    pt[0] = -1.0 + 2.0 * face;
    return 1.0;
  }
};

/**
  Defines the quadrature for a cubic beam
*/
class TACSBeamCubicQuadrature {
 public:
  static const int NUM_QUADRATURE_POINTS = 4;

  static int getNumParameters() { return 1; }
  static int getNumQuadraturePoints() { return 4; }
  static double getQuadratureWeight(int n) { return TacsGaussQuadWts4[n]; }
  static double getQuadraturePoint(int n, double pt[]) {
    pt[0] = TacsGaussQuadPts4[n];

    return TacsGaussQuadWts4[n];
  }
  static int getNumElementFaces() { return 2; }
  static int getNumFaceQuadraturePoints(int face) { return 1; }
  static double getFaceQuadraturePoint(int face, int n, double pt[],
                                       double t[]) {
    pt[0] = -1.0 + 2.0 * face;
    return 1.0;
  }
};

#endif  // TACS_BEAM_ELEMENT_QUADRATURE_H
