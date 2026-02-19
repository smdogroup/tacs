/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef PLANE_STRESS_TRACTION_H
#define PLANE_STRESS_TRACTION_H

/*
  The following file contains definitions for plane stress
  surface tractions
*/

/*
  The following class defines a traction for use with a quadrilateral
  element.

  The surface argument is an integer that denotes the element surface
  over which the traction should be applied. The surfaces are indexed
  as follows in the xi/eta parameter space:

  ---- 3 ----
  |         |
  0         1
  |         |
  ---- 2 ----

  Note that the integration within the element occurs over each face.
*/
template <int order>
class PSQuadTraction : public TACSElement {
 public:
  PSQuadTraction(int surf, TacsScalar _tx, TacsScalar _ty) {
    for (int k = 0; k < order; k++) {
      tx[k] = _tx;
      ty[k] = _ty;
    }
    initBaseDir(surf);
    // Set the knot locations
    if (order == 2) {
      knots[0] = -1.0;
      knots[1] = 1.0;
    } else if (order == 3) {
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    } else {
      // Set a co-sine spacing for the knot locations
      for (int k = 0; k < order; k++) {
        knots[k] = -cos(M_PI * k / (order - 1));
      }
    }
  }
  PSQuadTraction(int surf, TacsScalar _tx[], TacsScalar _ty[]) {
    for (int k = 0; k < order; k++) {
      tx[k] = _tx[k];
      ty[k] = _ty[k];
    }
    initBaseDir(surf);
    // Set the knot locations
    if (order == 2) {
      knots[0] = -1.0;
      knots[1] = 1.0;
    } else if (order == 3) {
      knots[0] = -1.0;
      knots[1] = 0.0;
      knots[2] = 1.0;
    } else {
      // Set a co-sine spacing for the knot locations
      for (int k = 0; k < order; k++) {
        knots[k] = -cos(M_PI * k / (order - 1));
      }
    }
  }

  // Get the number of displacements/nodes
  int numDisplacements() { return 2; }
  int numNodes() { return order * order; }

  // Add the residual from the forces
  // --------------------------------
  void addResidual(double time, TacsScalar res[], const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]) {
    // Retrieve the quadrature scheme of the appropriate order
    const double *gaussPts, *gaussWts;
    int numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Integrate over the specified element surface
    for (int n = 0; n < numGauss; n++) {
      double pt[2];
      pt[0] = gaussPts[n] * dir[0] + base[0];
      pt[1] = gaussPts[n] * dir[1] + base[1];

      // Evaluate the Lagrange basis in each direction
      double na[order], nb[order], dna[order], dnb[order];
      FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
      FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);

      // Calcualte the Jacobian at the current point
      const TacsScalar *x = Xpts;
      TacsScalar Xd[4] = {0.0, 0.0, 0.0, 0.0};
      for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
          Xd[0] += x[0] * dna[i] * nb[j];
          Xd[1] += x[0] * na[i] * dnb[j];

          Xd[2] += x[1] * dna[i] * nb[j];
          Xd[3] += x[1] * na[i] * dnb[j];
          x += 3;
        }
      }

      // Compute the derivative along each direction
      TacsScalar dx = Xd[0] * dir[0] + Xd[2] * dir[1];
      TacsScalar dy = Xd[1] * dir[0] + Xd[3] * dir[1];
      TacsScalar hsurf = gaussWts[n] * sqrt(dx * dx + dy * dy);

      // Calculate the traction at the current point
      TacsScalar Tx = 0.0, Ty = 0.0;
      double N[order];
      FElibrary::lagrangeSFKnots(N, gaussPts[n], knots, order);
      for (int i = 0; i < order; i++) {
        Tx += N[i] * tx[i];
        Ty += N[i] * ty[i];
      }

      // Add the result to the element
      TacsScalar *r = res;
      for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
          r[0] -= hsurf * na[i] * nb[j] * Tx;
          r[1] -= hsurf * na[i] * nb[j] * Ty;
          r += 2;
        }
      }
    }
  }

  // Function used for localizing the error to nodes with PU-weights
  // ---------------------------------------------------------------
  void addLocalizedError(double time, TacsScalar err[],
                         const TacsScalar adjoint[], const TacsScalar Xpts[],
                         const TacsScalar vars[]) {
    // Retrieve the quadrature scheme of the appropriate order
    const double *gaussPts, *gaussWts;
    int numGauss = FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Integrate over the specified element surface
    for (int n = 0; n < numGauss; n++) {
      double pt[2];
      pt[0] = gaussPts[n] * dir[0] + base[0];
      pt[1] = gaussPts[n] * dir[1] + base[1];

      // Evaluate the Lagrange basis in each direction
      double na[order], nb[order], dna[order], dnb[order];
      FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
      FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);

      // Calcualte the Jacobian at the current point
      const TacsScalar *x = Xpts;
      TacsScalar Xd[4] = {0.0, 0.0, 0.0, 0.0};
      TacsScalar Ax = 0.0, Ay = 0.0;
      for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
          Xd[0] += x[0] * dna[i] * nb[j];
          Xd[1] += x[0] * na[i] * dnb[j];

          Xd[2] += x[1] * dna[i] * nb[j];
          Xd[3] += x[1] * na[i] * dnb[j];
          x += 3;

          Ax += adjoint[2 * i] * na[i] * nb[j];
          Ay += adjoint[2 * i + 1] * na[i] * nb[j];
        }
      }

      // Compute the derivative along each direction
      TacsScalar dx = Xd[0] * dir[0] + Xd[2] * dir[1];
      TacsScalar dy = Xd[1] * dir[0] + Xd[3] * dir[1];
      TacsScalar hsurf = gaussWts[n] * sqrt(dx * dx + dy * dy);

      // Calculate the traction at the current point
      TacsScalar Tx = 0.0, Ty = 0.0;
      double N[order];
      FElibrary::lagrangeSFKnots(N, gaussPts[n], knots, order);
      for (int i = 0; i < order; i++) {
        Tx += N[i] * tx[i];
        Ty += N[i] * ty[i];
      }

      // Add the contribution to the residual - the minus sign is due
      // to the fact that this is a work term
      TacsScalar product = -hsurf * (Tx * Ax + Ty * Ay);

      // Add the product using the linear partition of unity basis
      // functions
      double Nerr[4];
      Nerr[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
      Nerr[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
      Nerr[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
      Nerr[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

      err[0] += Nerr[0] * product;
      err[order - 1] += Nerr[1] * product;
      err[order * (order - 1)] += Nerr[2] * product;
      err[order * order - 1] += Nerr[3] * product;
    }
  }

 private:
  // Set the base point and direction
  // --------------------------------
  void initBaseDir(int surf) {
    // Determine the base point and integration direction
    if (surf == 0 || surf == 1) {
      dir[0] = 0.0;
      dir[1] = 1.0;
    } else {
      dir[0] = 1.0;
      dir[1] = 0.0;
    }

    // Set the base point: The mid-point of the edge
    base[0] = base[1] = 0.0;
    if (surf == 0 || surf == 1) {
      base[0] = -1.0 + 2.0 * (surf % 2);
    } else {
      base[1] = -1.0 + 2.0 * (surf % 2);
    }
  }

  // The traction information
  TacsScalar tx[order], ty[order];

  // The parametric base point and direction along which the
  // integration will occur.
  double base[2], dir[2];

  double knots[order];
};

#endif  // PLANE_STRESS_TRACTION_H
