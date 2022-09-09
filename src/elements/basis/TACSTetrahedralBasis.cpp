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

#include "TACSTetrahedralBasis.h"

#include "TACSTetrahedronQuadrature.h"
#include "TACSTriangleQuadrature.h"

static void getFaceTangents(int face, double t[]) {
  if (face == 0) {
    // XY plane
    // return X and Y axes
    t[0] = 1.0;
    t[1] = 0.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 1.0;
    t[5] = 0.0;
  } else if (face == 1) {
    // XZ plane
    // return X and Z axes
    t[0] = 1.0;
    t[1] = 0.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 0.0;
    t[5] = 1.0;
  } else if (face == 2) {
    // slanted face
    t[0] = -1.0;
    t[1] = 0.0;
    t[2] = 1.0;
    t[3] = 0.0;
    t[4] = -1.0;
    t[5] = 1.0;
  } else if (face == 3) {
    // YZ plane
    t[0] = 0.0;
    t[1] = 1.0;
    t[2] = 0.0;
    t[3] = 0.0;
    t[4] = 0.0;
    t[5] = 1.0;
  }
}

/*
  Linear Tetrahedral basis class functions
*/
ElementLayout TACSLinearTetrahedralBasis::getLayoutType() {
  return TACS_TETRA_ELEMENT;
}

void TACSLinearTetrahedralBasis::getVisPoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = pt[1] = pt[2] = 0.0;
  } else if (n == 1) {
    pt[0] = 1.0;
    pt[1] = 0.0;
    pt[2] = 0.0;
  } else if (n == 2) {
    pt[0] = 0.0;
    pt[1] = 1.0;
    pt[2] = 0.0;
  } else if (n == 3) {
    pt[0] = 0.0;
    pt[1] = 0.0;
    pt[2] = 1.0;
  }
}

int TACSLinearTetrahedralBasis::getNumNodes() { return 4; }

int TACSLinearTetrahedralBasis::getNumParameters() { return 3; }

int TACSLinearTetrahedralBasis::getNumQuadraturePoints() { return 4; }

double TACSLinearTetrahedralBasis::getQuadratureWeight(int n) {
  if (n == 0 || n == 1 || n == 2 || n == 3) {
    return TacsTetrahedronWts4[0];
  }
  return 0.0;
}

double TACSLinearTetrahedralBasis::getQuadraturePoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = TacsTetrahedronPts4[0];
    pt[1] = TacsTetrahedronPts4[1];
    pt[2] = TacsTetrahedronPts4[2];
    return TacsTetrahedronWts4[0];
  } else if (n == 1) {
    pt[0] = TacsTetrahedronPts4[3];
    pt[1] = TacsTetrahedronPts4[4];
    pt[2] = TacsTetrahedronPts4[5];
    return TacsTetrahedronWts4[0];
  } else if (n == 2) {
    pt[0] = TacsTetrahedronPts4[6];
    pt[1] = TacsTetrahedronPts4[7];
    pt[2] = TacsTetrahedronPts4[8];
    return TacsTetrahedronWts4[0];
  } else if (n == 3) {
    pt[0] = TacsTetrahedronPts4[9];
    pt[1] = TacsTetrahedronPts4[10];
    pt[2] = TacsTetrahedronPts4[11];
    return TacsTetrahedronWts4[0];
  }
  return 0.0;
}

int TACSLinearTetrahedralBasis::getNumElementFaces() { return 4; }

int TACSLinearTetrahedralBasis::getNumFaceQuadraturePoints(int face) {
  return 3;
}

double TACSLinearTetrahedralBasis::getFaceQuadraturePoint(int face, int n,
                                                          double pt[],
                                                          double t[]) {
  getFaceTangents(face, t);

  if (face == 0) {
    // XY plane face
    pt[2] = 0;
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
    return TacsTriangleWts3[0];
  } else if (face == 1) {
    // XZ plane face
    pt[1] = 0;
    if (n == 0) {
      pt[0] = TacsTrianglePts3[0];
      pt[2] = TacsTrianglePts3[1];
    } else if (n == 1) {
      pt[0] = TacsTrianglePts3[2];
      pt[2] = TacsTrianglePts3[3];
    } else if (n == 2) {
      pt[0] = TacsTrianglePts3[4];
      pt[2] = TacsTrianglePts3[5];
    }
    return TacsTriangleWts3[0];
  } else if (face == 2) {
    // points should be the same but multiplied by rt 2
    // weight is already multiplied by 1/2 (area of regular tri)
    // needs to me multiplied by 2 to remove that and then by the area of
    // slanted tri
    if (n == 0) {
      pt[0] = 1.0 / 6.0;
      pt[1] = 1.0 / 6.0;
      pt[2] = 2.0 / 3.0;
    } else if (n == 1) {
      pt[0] = 2.0 / 3.0;
      pt[1] = 1.0 / 6.0;
      pt[2] = 1.0 / 6.0;
    } else if (n == 2) {
      pt[0] = 1.0 / 6.0;
      pt[1] = 2.0 / 3.0;
      pt[2] = 1.0 / 6.0;
    }
    return sqrt(3.0) * TacsTriangleWts3[0];
  } else if (face == 3) {
    // YZ face
    pt[0] = 0;
    if (n == 0) {
      pt[1] = TacsTrianglePts3[0];
      pt[2] = TacsTrianglePts3[1];
    } else if (n == 1) {
      pt[1] = TacsTrianglePts3[2];
      pt[2] = TacsTrianglePts3[3];
    } else if (n == 2) {
      pt[1] = TacsTrianglePts3[4];
      pt[2] = TacsTrianglePts3[5];
    }
    return TacsTriangleWts3[0];
  }
  return 0.0;
}

void TACSLinearTetrahedralBasis::computeBasis(const double pt[], double N[]) {
  N[0] = 1 - pt[0] - pt[1] - pt[2];
  N[1] = pt[0];
  N[2] = pt[1];
  N[3] = pt[2];
}

void TACSLinearTetrahedralBasis::computeBasisGradient(const double pt[],
                                                      double N[],
                                                      double Nxi[]) {
  N[0] = 1 - pt[0] - pt[1] - pt[2];
  N[1] = pt[0];
  N[2] = pt[1];
  N[3] = pt[2];

  Nxi[0] = -1.0;
  Nxi[1] = -1.0;
  Nxi[2] = -1.0;
  Nxi[3] = 1.0;
  Nxi[4] = 0.0;
  Nxi[5] = 0.0;
  Nxi[6] = 0.0;
  Nxi[7] = 1.0;
  Nxi[8] = 0.0;
  Nxi[9] = 0.0;
  Nxi[10] = 0.0;
  Nxi[11] = 1.0;
}

/*
  Quadratic Tetrahedral basis class functions
*/
ElementLayout TACSQuadraticTetrahedralBasis::getLayoutType() {
  return TACS_TETRA_QUADRATIC_ELEMENT;
}

/*
  Get the parametric node locations for visualization

*/
void TACSQuadraticTetrahedralBasis::getVisPoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = pt[1] = pt[2] = 0.0;
  } else if (n == 1) {
    pt[0] = 1.0;
    pt[1] = 0.0;
    pt[2] = 0.0;
  } else if (n == 2) {
    pt[0] = 0.0;
    pt[1] = 1.0;
    pt[2] = 0.0;
  } else if (n == 3) {
    pt[0] = 0.0;
    pt[1] = 0.0;
    pt[2] = 1.0;
  } else if (n == 4) {
    pt[0] = 0.5;
    pt[1] = 0.0;
    pt[2] = 0.0;
  } else if (n == 5) {
    pt[0] = 0.5;
    pt[1] = 0.5;
    pt[2] = 0.0;
  } else if (n == 6) {
    pt[0] = 0.0;
    pt[1] = 0.5;
    pt[2] = 0.0;
  } else if (n == 7) {
    pt[0] = 0.0;
    pt[1] = 0.0;
    pt[2] = 0.5;
  } else if (n == 8) {
    pt[0] = 0.5;
    pt[1] = 0.0;
    pt[2] = 0.5;
  } else if (n == 9) {
    pt[0] = 0.0;
    pt[1] = 0.5;
    pt[2] = 0.5;
  }
}

int TACSQuadraticTetrahedralBasis::getNumNodes() { return 10; }

int TACSQuadraticTetrahedralBasis::getNumParameters() { return 3; }

int TACSQuadraticTetrahedralBasis::getNumQuadraturePoints() { return 5; }

double TACSQuadraticTetrahedralBasis::getQuadratureWeight(int n) {
  if (n == 0) {
    return TacsTetrahedronWts5[0];
  } else if (n == 1 || n == 2 || n == 3 || n == 4) {
    return TacsTetrahedronWts5[1];
  }
  return 0.0;
}

double TACSQuadraticTetrahedralBasis::getQuadraturePoint(int n, double pt[]) {
  if (n == 0) {
    pt[0] = TacsTetrahedronPts5[0];
    pt[1] = TacsTetrahedronPts5[1];
    pt[2] = TacsTetrahedronPts5[2];
    return TacsTetrahedronWts5[0];
  } else if (n == 1) {
    pt[0] = TacsTetrahedronPts5[3];
    pt[1] = TacsTetrahedronPts5[4];
    pt[2] = TacsTetrahedronPts5[5];
    return TacsTetrahedronWts5[1];
  } else if (n == 2) {
    pt[0] = TacsTetrahedronPts5[6];
    pt[1] = TacsTetrahedronPts5[7];
    pt[2] = TacsTetrahedronPts5[8];
    return TacsTetrahedronWts5[1];
  } else if (n == 3) {
    pt[0] = TacsTetrahedronPts5[9];
    pt[1] = TacsTetrahedronPts5[10];
    pt[2] = TacsTetrahedronPts5[11];
    return TacsTetrahedronWts5[1];
  } else if (n == 4) {
    pt[0] = TacsTetrahedronPts5[12];
    pt[1] = TacsTetrahedronPts5[13];
    pt[2] = TacsTetrahedronPts5[14];
    return TacsTetrahedronWts5[1];
  }
  return 0.0;
}

int TACSQuadraticTetrahedralBasis::getNumElementFaces() { return 4; }

int TACSQuadraticTetrahedralBasis::getNumFaceQuadraturePoints(int face) {
  return 4;
}

double TACSQuadraticTetrahedralBasis::getFaceQuadraturePoint(int face, int n,
                                                             double pt[],
                                                             double t[]) {
  getFaceTangents(face, t);

  if (face == 0) {
    // XY plane face
    pt[2] = 0;
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
      return TacsTriangleWts4[1];
    } else if (n == 3) {
      pt[0] = TacsTrianglePts4[6];
      pt[1] = TacsTrianglePts4[7];
      return TacsTriangleWts4[1];
    }
  } else if (face == 1) {
    // XZ plane face
    pt[1] = 0;
    if (n == 0) {
      pt[0] = TacsTrianglePts4[0];
      pt[2] = TacsTrianglePts4[1];
      return TacsTriangleWts4[0];
    } else if (n == 1) {
      pt[0] = TacsTrianglePts4[2];
      pt[2] = TacsTrianglePts4[3];
      return TacsTriangleWts4[1];
    } else if (n == 2) {
      pt[0] = TacsTrianglePts4[4];
      pt[2] = TacsTrianglePts4[5];
      return TacsTriangleWts4[1];
    } else if (n == 3) {
      pt[0] = TacsTrianglePts4[6];
      pt[2] = TacsTrianglePts4[7];
      return TacsTriangleWts4[1];
    }
  } else if (face == 2) {
    // slanted face -- same procedure as linear
    if (n == 0) {
      pt[0] = 1.0 / 3.0;
      pt[1] = 1.0 / 3.0;
      pt[2] = 1.0 / 3.0;
      return sqrt(3.0) * TacsTriangleWts4[0];
    } else if (n == 1) {
      pt[0] = 1.0 / 5.0;
      pt[1] = 1.0 / 5.0;
      pt[2] = 3.0 / 5.0;
      return sqrt(3.0) * TacsTriangleWts4[1];
    } else if (n == 2) {
      pt[0] = 3.0 / 5.0;
      pt[1] = 1.0 / 5.0;
      pt[2] = 1.0 / 5.0;
      return sqrt(3.0) * TacsTriangleWts4[1];
    } else if (n == 3) {
      pt[0] = 1.0 / 5.0;
      pt[1] = 3.0 / 5.0;
      pt[2] = 1.0 / 5.0;
      return sqrt(3.0) * TacsTriangleWts4[1];
    }
  } else if (face == 3) {
    // YZ face
    pt[0] = 0;
    if (n == 0) {
      pt[1] = TacsTrianglePts4[0];
      pt[2] = TacsTrianglePts4[1];
      return TacsTriangleWts4[0];
    } else if (n == 1) {
      pt[1] = TacsTrianglePts4[2];
      pt[2] = TacsTrianglePts4[3];
      return TacsTriangleWts4[1];
    } else if (n == 2) {
      pt[1] = TacsTrianglePts4[4];
      pt[2] = TacsTrianglePts4[5];
      return TacsTriangleWts4[1];
    } else if (n == 3) {
      pt[1] = TacsTrianglePts4[6];
      pt[2] = TacsTrianglePts4[7];
      return TacsTriangleWts4[1];
    }
  }
  return 0.0;
}

void TACSQuadraticTetrahedralBasis::computeBasis(const double pt[],
                                                 double N[]) {
  double l0 = 1 - pt[0] - pt[1] - pt[2];
  double l1 = pt[0];
  double l2 = pt[1];
  double l3 = pt[2];

  // Corner nodes
  N[0] = l0 * (2 * l0 - 1);
  N[1] = l1 * (2 * l1 - 1);
  N[2] = l2 * (2 * l2 - 1);
  N[3] = l3 * (2 * l3 - 1);

  // Mid-side nodes
  N[4] = 4 * l1 * l0;
  N[5] = 4 * l1 * l2;
  N[6] = 4 * l2 * l0;
  // N[7] = 4*l3*l1;
  N[7] = 4 * l3 * l0;
  // N[8] = 4*l2*l3;
  N[8] = 4 * l1 * l3;
  // N[9] = 4*l3*l0;
  N[9] = 4 * l3 * l2;
}

void TACSQuadraticTetrahedralBasis::computeBasisGradient(const double pt[],
                                                         double N[],
                                                         double Nxi[]) {
  double l0 = 1 - pt[0] - pt[1] - pt[2];
  double l1 = pt[0];
  double l2 = pt[1];
  double l3 = pt[2];

  // Corner nodes
  N[0] = l0 * (2 * l0 - 1);
  N[1] = l1 * (2 * l1 - 1);
  N[2] = l2 * (2 * l2 - 1);
  N[3] = l3 * (2 * l3 - 1);

  // Mid-side nodes
  N[4] = 4 * l1 * l0;
  N[5] = 4 * l1 * l2;
  N[6] = 4 * l2 * l0;
  N[7] = 4 * l3 * l0;
  N[8] = 4 * l1 * l3;
  N[9] = 4 * l3 * l2;

  // Corner node derivatives
  Nxi[0] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
  Nxi[1] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
  Nxi[2] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
  Nxi[3] = 4 * pt[0] - 1;
  Nxi[4] = 0;
  Nxi[5] = 0;
  Nxi[6] = 0;
  Nxi[7] = 4 * pt[1] - 1;
  Nxi[8] = 0;
  Nxi[9] = 0;
  Nxi[10] = 0;
  Nxi[11] = 4 * pt[2] - 1;

  // Mid node derivatives
  Nxi[12] = -4 * (2 * pt[0] + pt[1] + pt[2] - 1);
  Nxi[13] = -4 * pt[0];
  Nxi[14] = -4 * pt[0];

  Nxi[15] = 4 * pt[1];
  Nxi[16] = 4 * pt[0];
  Nxi[17] = 0.0;

  Nxi[18] = -4 * pt[1];
  Nxi[19] = -4 * (pt[0] + 2 * pt[1] + pt[2] - 1);
  Nxi[20] = -4 * pt[1];

  Nxi[21] = -4 * pt[2];
  Nxi[22] = -4 * pt[2];
  Nxi[23] = -4 * (pt[0] + pt[1] + 2 * pt[2] - 1);

  Nxi[24] = 4 * pt[2];
  Nxi[25] = 0.0;
  Nxi[26] = 4 * pt[0];

  Nxi[27] = 0.0;
  Nxi[28] = 4 * pt[2];
  Nxi[29] = 4 * pt[1];
}

/*
  Cubic Tetrahedral basis class functions
*/
ElementLayout TACSCubicTetrahedralBasis::getLayoutType() {
  return TACS_TETRA_CUBIC_ELEMENT;
}

void TACSCubicTetrahedralBasis::getVisPoint(int n, double pt[]) {
  pt[0] = pt[1] = pt[2] = 0.0;
}

int TACSCubicTetrahedralBasis::getNumNodes() { return 3; }

int TACSCubicTetrahedralBasis::getNumParameters() { return 2; }

int TACSCubicTetrahedralBasis::getNumQuadraturePoints() { return 1; }

double TACSCubicTetrahedralBasis::getQuadratureWeight(int n) { return 1.0; }

double TACSCubicTetrahedralBasis::getQuadraturePoint(int n, double pt[]) {
  return 1.0;
}

int TACSCubicTetrahedralBasis::getNumElementFaces() { return 3; }

int TACSCubicTetrahedralBasis::getNumFaceQuadraturePoints(int face) {
  return 1;
}

double TACSCubicTetrahedralBasis::getFaceQuadraturePoint(int face, int n,
                                                         double pt[],
                                                         double t[]) {
  return 0.1;
}

void TACSCubicTetrahedralBasis::computeBasis(const double pt[], double N[]) {}

void TACSCubicTetrahedralBasis::computeBasisGradient(const double pt[],
                                                     double N[], double Nxi[]) {
}
