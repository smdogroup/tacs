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

#include "TACSElementMappedBasis.h"

#include "TACSElementAlgebra.h"

TacsScalar TACSElementMappedBasis::getFaceNormal(int face, int n,
                                                 const TacsScalar Xpts[],
                                                 TacsScalar X[],
                                                 TacsScalar Xd[],
                                                 TacsScalar normal[]) {
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();

  double pt[3];
  double tangents[2 * 3];
  getFaceQuadraturePoint(face, n, pt, tangents);

  double N[256], Nxi[3 * 256];
  computeBasisGradient(pt, N, Nxi);

  return computeMappedFaceNormal(num_params, num_nodes, N, Nxi, Xpts, tangents,
                                 X, Xd, normal);
}

TacsScalar TACSElementMappedBasis::computeMappedFaceNormal(
    const int num_params, const int num_nodes, const double N[],
    const double Nxi[], const TacsScalar Xpts[], const double tangents[],
    TacsScalar X[], TacsScalar Xd[], TacsScalar n[]) {
  if (num_params == 2) {
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;

    // Loop over each quadrature point for each basis function
    TacsScalar Xd1[3] = {0.0, 0.0, 0.0};
    TacsScalar Xd2[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < num_nodes; i++) {
      X[0] += N[0] * Xpts[0];
      X[1] += N[0] * Xpts[1];
      X[2] += N[0] * Xpts[2];

      Xd1[0] += Nxi[0] * Xpts[0];
      Xd1[1] += Nxi[0] * Xpts[1];
      Xd1[2] += Nxi[0] * Xpts[2];

      Xd2[0] += Nxi[1] * Xpts[0];
      Xd2[1] += Nxi[1] * Xpts[1];
      Xd2[2] += Nxi[1] * Xpts[2];

      Xpts += 3;
      Nxi += 2;
      N++;
    }

    // Set the values of Xd
    Xd[0] = Xd1[0];
    Xd[2] = Xd1[1];
    Xd[4] = Xd1[2];

    Xd[1] = Xd2[0];
    Xd[3] = Xd2[1];
    Xd[5] = Xd2[2];

    // Compute the surface normal = tangent x refdir
    TacsScalar sn[3];
    crossProduct(Xd1, Xd2, sn);

    // Compute the tangent direction
    TacsScalar t[3];
    t[0] = Xd[0] * tangents[0] + Xd[1] * tangents[1];
    t[1] = Xd[2] * tangents[0] + Xd[3] * tangents[1];
    t[2] = Xd[4] * tangents[0] + Xd[5] * tangents[1];

    // Compute the edge normal n = t x n
    crossProduct(t, sn, n);

    // Compute the norm of the vector
    TacsScalar A = sqrt(vec2Dot(t, t));
    TacsScalar invA = 1.0 / A;

    n[0] = n[0] * invA;
    n[1] = n[1] * invA;
    n[2] = n[2] * invA;

    return A;
  } else if (num_params == 1) {
    // Zero the values of the coordinate and its derivative
    X[0] = X[1] = X[2] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = 0.0;

    // Loop over each quadrature point for each basis function
    for (int i = 0; i < num_nodes; i++) {
      X[0] += N[0] * Xpts[0];
      X[1] += N[0] * Xpts[1];
      X[2] += N[0] * Xpts[2];

      Xd[0] += Nxi[0] * Xpts[0];
      Xd[1] += Nxi[0] * Xpts[1];
      Xd[2] += Nxi[0] * Xpts[2];

      Xpts += 3;
      Nxi++;
      N++;
    }

    TacsScalar A = sqrt(vec2Dot(Xd, Xd));
    TacsScalar invA = 1.0 / A;

    n[0] = invA * tangents[0] * Xd[0];
    n[1] = invA * tangents[0] * Xd[1];
    n[2] = invA * tangents[0] * Xd[2];

    return 1.0;
  }

  return 0.0;
}

void TACSElementMappedBasis::addFaceNormalXptSens(
    int face, int n, const TacsScalar A, const TacsScalar Xd[],
    const TacsScalar normal[], const TacsScalar dfdA, const TacsScalar dfdX[],
    const TacsScalar dfdXd[], const TacsScalar dfdn[], TacsScalar dfdXpts[]) {}

TacsScalar TACSElementMappedBasis::getJacobianTransform(int n,
                                                        const double pt[],
                                                        const TacsScalar Xpts[],
                                                        TacsScalar Xd[],
                                                        TacsScalar J[]) {
  const int num_params = getNumParameters();
  const int num_nodes = getNumNodes();
  double N[256], Nxi[3 * 256];
  computeBasisGradient(pt, N, Nxi);

  return computeMappedJacobianTransform(num_params, num_nodes, Nxi, Xpts, Xd,
                                        J);
}

TacsScalar TACSElementMappedBasis::computeMappedJacobianTransform(
    const int num_params, const int num_nodes, const double Nxi[],
    const TacsScalar Xpts[], TacsScalar Xd[], TacsScalar J[]) {
  if (num_params == 2) {
    // Loop over each quadrature point for each basis function
    TacsScalar Xd1[3] = {0.0, 0.0, 0.0};
    TacsScalar Xd2[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < num_nodes; i++) {
      Xd1[0] += Nxi[0] * Xpts[0];
      Xd1[1] += Nxi[0] * Xpts[1];
      Xd1[2] += Nxi[0] * Xpts[2];

      Xd2[0] += Nxi[1] * Xpts[0];
      Xd2[1] += Nxi[1] * Xpts[1];
      Xd2[2] += Nxi[1] * Xpts[2];

      Xpts += 3;
      Nxi += 2;
    }

    // Set the values of Xd
    Xd[0] = Xd1[0];
    Xd[2] = Xd1[1];
    Xd[4] = Xd1[2];

    Xd[1] = Xd2[0];
    Xd[3] = Xd2[1];
    Xd[5] = Xd2[2];

    // Compute the surface normal = tangent x refdir
    TacsScalar sn[3];
    crossProduct(Xd1, Xd2, sn);

    TacsScalar A = sqrt(vec3Dot(sn, sn));

    return A;
  } else if (num_params == 1) {
    // Zero the values of the coordinate and its derivative
    Xd[0] = Xd[1] = Xd[2] = 0.0;

    // Loop over each quadrature point for each basis function
    for (int i = 0; i < num_nodes; i++) {
      Xd[0] += Nxi[0] * Xpts[0];
      Xd[1] += Nxi[0] * Xpts[1];
      Xd[2] += Nxi[0] * Xpts[2];

      Xpts += 3;
      Nxi++;
    }

    TacsScalar A = sqrt(vec2Dot(Xd, Xd));
    J[0] = 1.0 / A;

    return A;
  }

  return 0.0;
}

void TACSElementMappedBasis::addJacobianTransformXptSens(
    int n, const double pt[], const TacsScalar Xd[], const TacsScalar J[],
    TacsScalar dfddetJ, const TacsScalar dfXd[], const TacsScalar dfdJ[],
    TacsScalar dfdXpts[]) {}
