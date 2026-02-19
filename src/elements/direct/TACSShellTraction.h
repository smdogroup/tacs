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

#ifndef TACS_SHELL_TRACTION_H
#define TACS_SHELL_TRACTION_H

/*
  Shell element traction class
*/

#include "ShellUtils.h"
#include "TACSElement.h"

/*
  TACSShellTraction class

  This class defines a general shell traction. The traction may apply
  a force either normal to or in-plane of the shell.
*/
template <int order>
class TACSShellTraction : public TACSElement {
 public:
  // Set the number of nodes
  static const int NUM_NODES = order * order;

  // Constructor for the shell element
  // ---------------------------------
  TACSShellTraction(TacsScalar _tx[], TacsScalar _ty[], TacsScalar _tz[]) {
    self = NULL;
    evalf = NULL;
    memcpy(tx, _tx, NUM_NODES * sizeof(TacsScalar));
    memcpy(ty, _ty, NUM_NODES * sizeof(TacsScalar));
    memcpy(tz, _tz, NUM_NODES * sizeof(TacsScalar));
  }
  TACSShellTraction(TacsScalar _tx, TacsScalar _ty, TacsScalar _tz) {
    self = NULL;
    evalf = NULL;
    for (int i = 0; i < NUM_NODES; i++) {
      tx[i] = _tx;
      ty[i] = _ty;
      tz[i] = _tz;
    }
  }
  TACSShellTraction(void *_self,
                    void (*_evalf)(void *, const TacsScalar *, TacsScalar *)) {
    self = _self;
    evalf = _evalf;
  }

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements() { return 6; }
  int numNodes() { return NUM_NODES; }
  int numStresses() { return 0; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *Te, TacsScalar *Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]) {
    *Te = 0.0, *Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual(double time, TacsScalar res[], const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]) {
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for (int m = 0; m < order; m++) {
      for (int n = 0; n < order; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
        Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
        Tensor::normalize3D(&Xd[6]);

        // Compute the determinant of the Jacobian
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;
        if (evalf) {
          TacsScalar X[3] = {0.0, 0.0, 0.0};
          for (int i = 0; i < order * order; i++) {
            X[0] += N[i] * Xpts[3 * i];
            X[1] += N[i] * Xpts[3 * i + 1];
            X[2] += N[i] * Xpts[3 * i + 2];
          }
          TacsScalar t[3];
          evalf(self, X, t);
          Tx = t[0];
          Ty = t[1];
          Tz = t[2];
        } else {
          for (int i = 0; i < NUM_NODES; i++) {
            Tx += tx[i] * N[i];
            Ty += ty[i] * N[i];
            Tz += tz[i] * N[i];
          }
        }

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar *r = res;
        for (int i = 0; i < NUM_NODES; i++) {
          r[0] -= h * Tx * N[i];
          r[1] -= h * Ty * N[i];
          r[2] -= h * Tz * N[i];
          r += 6;
        }
      }
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian(double time, TacsScalar J[], double alpha, double beta,
                   double gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]) {}

  // Function used for localizing the error to nodes with PU-weights
  // ---------------------------------------------------------------
  void addLocalizedError(double time, TacsScalar err[],
                         const TacsScalar adjoint[], const TacsScalar Xpts[],
                         const TacsScalar vars[]) {
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    int npts = FElibrary::getGaussPtsWts(order + 1, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for (int m = 0; m < npts; m++) {
      for (int n = 0; n < npts; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
        Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
        Tensor::normalize3D(&Xd[6]);

        // Compute the determinant of the Jacobian
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;
        if (evalf) {
          TacsScalar X[3] = {0.0, 0.0, 0.0};
          for (int i = 0; i < order * order; i++) {
            X[0] += N[i] * Xpts[3 * i];
            X[1] += N[i] * Xpts[3 * i + 1];
            X[2] += N[i] * Xpts[3 * i + 2];
          }
          TacsScalar t[3];
          evalf(self, X, t);
          Tx = t[0];
          Ty = t[1];
          Tz = t[2];
        } else {
          for (int i = 0; i < NUM_NODES; i++) {
            Tx += tx[i] * N[i];
            Ty += ty[i] * N[i];
            Tz += tz[i] * N[i];
          }
        }

        // Compute the adjoint terms
        TacsScalar Ax = 0.0, Ay = 0.0, Az = 0.0;
        for (int i = 0; i < NUM_NODES; i++) {
          Ax += adjoint[6 * i] * N[i];
          Ay += adjoint[6 * i + 1] * N[i];
          Az += adjoint[6 * i + 2] * N[i];
        }

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar product = -h * (Tx * Ax + Ty * Ay + Tz * Az);

        // Add the product using the linear partition of unity basis
        // functions
        double Nerr[4];
        Nerr[0] = 0.25 * (1.0 - pt[0]) * (1.0 - pt[1]);
        Nerr[1] = 0.25 * (1.0 + pt[0]) * (1.0 - pt[1]);
        Nerr[2] = 0.25 * (1.0 - pt[0]) * (1.0 + pt[1]);
        Nerr[3] = 0.25 * (1.0 + pt[0]) * (1.0 + pt[1]);

        // Add the contributions to the error
        for (int node = 0; node < 4; node++) {
          err[(node % 2) * (order - 1) + (node / 2) * order * (order - 1)] +=
              Nerr[node] * product;
        }
      }
    }
  }

 private:
  TacsScalar tx[NUM_NODES], ty[NUM_NODES], tz[NUM_NODES];

  // The data/function for the right-hand-side
  void *self;
  void (*evalf)(void *, const TacsScalar *, TacsScalar *);
};

/*
  TACSShellPressure class

  This class defines a surface pressure. The pressure is always normal
  to the surface.
*/
template <int order>
class TACSShellPressure : public TACSElement {
 public:
  // Set the number of nodes
  static const int NUM_NODES = order * order;

  // Constructor for the shell element
  // ---------------------------------
  TACSShellPressure(TacsScalar _p[]) {
    memcpy(p, _p, NUM_NODES * sizeof(TacsScalar));
  }
  TACSShellPressure(TacsScalar _p) {
    for (int i = 0; i < NUM_NODES; i++) {
      p[i] = _p;
    }
  }

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements() { return 6; }
  int numNodes() { return NUM_NODES; }
  int numStresses() { return 0; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *Te, TacsScalar *Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]) {
    *Te = 0.0, *Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual(double time, TacsScalar res[], const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]) {
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for (int m = 0; m < order; m++) {
      for (int n = 0; n < order; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
        Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
        Tensor::normalize3D(&Xd[6]);

        // Compute the determinant of the Jacobian
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar P = 0.0;
        for (int i = 0; i < NUM_NODES; i++) {
          P += p[i] * N[i];
        }

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar *r = res;
        for (int i = 0; i < NUM_NODES; i++) {
          r[0] -= h * P * Xd[6] * N[i];
          r[1] -= h * P * Xd[7] * N[i];
          r[2] -= h * P * Xd[8] * N[i];
          r += 6;
        }
      }
    }
  }

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian(double time, TacsScalar J[], double alpha, double beta,
                   double gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[]) {}

  // Function used for localizing the error to nodes with PU-weights
  // ---------------------------------------------------------------
  void addLocalizedError(double time, TacsScalar err[],
                         const TacsScalar adjoint[], const TacsScalar Xpts[],
                         const TacsScalar vars[]) {
    // Get the quadrature points and weights
    const double *gaussPts, *gaussWts;
    FElibrary::getGaussPtsWts(order, &gaussPts, &gaussWts);

    // Add the residual due to the shell traction
    for (int m = 0; m < order; m++) {
      for (int n = 0; n < order; n++) {
        // Set the quadrature point
        double pt[2];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
        TacsScalar X[3], Xd[9];
        shellutils::shell_jacobian(order, X, Xd, N, Na, Nb, pt, Xpts);

        // Determine the normal direction and normalize it
        Tensor::crossProduct3D(&Xd[6], &Xd[0], &Xd[3]);
        Tensor::normalize3D(&Xd[6]);

        // Compute the determinant of the Jacobian
        TacsScalar h = FElibrary::jacobian3d(Xd);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar P = 0.0;
        TacsScalar Ax = 0.0, Ay = 0.0, Az = 0.0;
        for (int i = 0; i < NUM_NODES; i++) {
          P += p[i] * N[i];
          Ax += adjoint[6 * i] * N[i];
          Ay += adjoint[6 * i + 1] * N[i];
          Az += adjoint[6 * i + 2] * N[i];
        }

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = P * Xd[6], Ty = P * Xd[7], Tz = P * Xd[8];

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        TacsScalar product = -h * (Tx * Ax + Ty * Ay + Tz * Az);

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
  }

 private:
  TacsScalar p[NUM_NODES];
};

#endif  // TACS_SHELL_TRACTION
