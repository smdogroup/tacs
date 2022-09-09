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

#ifndef TACS_3D_TRACTION_H
#define TACS_3D_TRACTION_H

#include "FElibrary.h"
#include "TACSElement.h"
#include "TensorToolbox.h"

/*
  The surface traction class for 3D elements
*/
template <int order>
class TACS3DTraction : public TACSElement {
 public:
  static const int NUM_NODES = order * order * order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;

  // Constructor for the shell element
  TACS3DTraction(int _surface, TacsScalar _tx[], TacsScalar _ty[],
                 TacsScalar _tz[]) {
    surface = surface;
    memcpy(tx, _tx, order * order * sizeof(TacsScalar));
    memcpy(ty, _ty, order * order * sizeof(TacsScalar));
    memcpy(tz, _tz, order * order * sizeof(TacsScalar));

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
  TACS3DTraction(int _surface, TacsScalar _tx, TacsScalar _ty, TacsScalar _tz) {
    surface = _surface;
    for (int i = 0; i < order * order; i++) {
      tx[i] = _tx;
      ty[i] = _ty;
      tz[i] = _tz;
    }

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

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements() { return 3; }
  int numNodes() { return NUM_NODES; }
  int numStresses() { return 0; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *Te, TacsScalar *Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]) {
    *Te = 0.0, *Pe = 0.0;
  }

  /*
    Compute the shape functions and their derivatives w.r.t. the
    parametric element location
  */
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[]) {
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        N[0] = na[i] * nb[j];
        Na[0] = dna[i] * nb[j];
        Nb[0] = na[i] * dnb[j];
        N++;
        Na++;
        Nb++;
      }
    }
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
        double pt[3];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[order * order];
        double Na[order * order], Nb[order * order];
        getShapeFunctions(pt, N, Na, Nb);

        TacsScalar Xa[3], Xb[3];
        Xa[0] = Xa[1] = Xa[2] = 0.0;
        Xb[0] = Xb[1] = Xb[2] = 0.0;

        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[j + k * order] * Xpts[3 * node];
              Xa[1] += Na[j + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[j + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[j + k * order] * Xpts[3 * node];
              Xb[1] += Nb[j + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[j + k * order] * Xpts[3 * node + 2];
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + k * order] * Xpts[3 * node];
              Xa[1] += Na[i + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + k * order] * Xpts[3 * node];
              Xb[1] += Nb[i + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + k * order] * Xpts[3 * node + 2];
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + j * order] * Xpts[3 * node];
              Xa[1] += Na[i + j * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + j * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + j * order] * Xpts[3 * node];
              Xb[1] += Nb[i + j * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + j * order] * Xpts[3 * node + 2];
            }
          }
        }

        // Determine the normal direction
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, Xa, Xb);
        TacsScalar h = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                            normal[2] * normal[2]);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;
        for (int i = 0; i < order * order; i++) {
          Tx += tx[i] * N[i];
          Ty += ty[i] * N[i];
          Tz += tz[i] * N[i];
        }

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[j + k * order];
              res[3 * node + 1] -= h * Ty * N[j + k * order];
              res[3 * node + 2] -= h * Tz * N[j + k * order];
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + k * order];
              res[3 * node + 1] -= h * Ty * N[i + k * order];
              res[3 * node + 2] -= h * Tz * N[i + k * order];
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + j * order];
              res[3 * node + 1] -= h * Ty * N[i + j * order];
              res[3 * node + 2] -= h * Tz * N[i + j * order];
            }
          }
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

 private:
  int surface;
  TacsScalar tx[order * order];
  TacsScalar ty[order * order];
  TacsScalar tz[order * order];

  // The knot locations for the basis functions
  double knots[order];
};

/*
  The surface pressure class for 3D elements
*/
template <int order>
class TACS3DPressureTraction : public TACSElement {
 public:
  static const int NUM_NODES = order * order * order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;

  // Constructor for the shell element
  TACS3DPressureTraction(int _surface, TacsScalar _pressure) {
    surface = surface;
    pressure = _pressure;

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

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements() { return 3; }
  int numNodes() { return NUM_NODES; }
  int numStresses() { return 0; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *Te, TacsScalar *Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]) {
    *Te = 0.0, *Pe = 0.0;
  }

  /*
    Compute the shape functions and their derivatives w.r.t. the
    parametric element location
  */
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[]) {
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        N[0] = na[i] * nb[j];
        Na[0] = dna[i] * nb[j];
        Nb[0] = na[i] * dnb[j];
        N++;
        Na++;
        Nb++;
      }
    }
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
        double pt[3];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[order * order];
        double Na[order * order], Nb[order * order];
        getShapeFunctions(pt, N, Na, Nb);

        TacsScalar Xa[3], Xb[3];
        Xa[0] = Xa[1] = Xa[2] = 0.0;
        Xb[0] = Xb[1] = Xb[2] = 0.0;

        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[j + k * order] * Xpts[3 * node];
              Xa[1] += Na[j + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[j + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[j + k * order] * Xpts[3 * node];
              Xb[1] += Nb[j + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[j + k * order] * Xpts[3 * node + 2];
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + k * order] * Xpts[3 * node];
              Xa[1] += Na[i + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + k * order] * Xpts[3 * node];
              Xb[1] += Nb[i + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + k * order] * Xpts[3 * node + 2];
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + j * order] * Xpts[3 * node];
              Xa[1] += Na[i + j * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + j * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + j * order] * Xpts[3 * node];
              Xb[1] += Nb[i + j * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + j * order] * Xpts[3 * node + 2];
            }
          }
        }

        // Determine the normal direction
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, Xa, Xb);
        TacsScalar h = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                            normal[2] * normal[2]);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the
        // quadrature point within the element
        TacsScalar Tx = -pressure * normal[0];
        TacsScalar Ty = -pressure * normal[1];
        TacsScalar Tz = -pressure * normal[2];

        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[j + k * order];
              res[3 * node + 1] -= h * Ty * N[j + k * order];
              res[3 * node + 2] -= h * Tz * N[j + k * order];
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + k * order];
              res[3 * node + 1] -= h * Ty * N[i + k * order];
              res[3 * node + 2] -= h * Tz * N[i + k * order];
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + j * order];
              res[3 * node + 1] -= h * Ty * N[i + j * order];
              res[3 * node + 2] -= h * Tz * N[i + j * order];
            }
          }
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

 private:
  int surface;
  TacsScalar pressure;

  // The knot locations for the basis functions
  double knots[order];
};

template <int order>
class TACS3DBoundingTraction : public TACSElement {
 public:
  static const int NUM_NODES = order * order * order;

  static const int U_NEGATIVE = 0;
  static const int U_POSITIVE = 1;
  static const int V_NEGATIVE = 2;
  static const int V_POSITIVE = 3;
  static const int W_NEGATIVE = 4;
  static const int W_POSITIVE = 5;

  // Constructor for the shell element
  TACS3DBoundingTraction(int _surface, TacsScalar _tx[], TacsScalar _ty[],
                         TacsScalar _tz[], TacsScalar _box[]) {
    surface = surface;
    memcpy(tx, _tx, order * order * sizeof(TacsScalar));
    memcpy(ty, _ty, order * order * sizeof(TacsScalar));
    memcpy(tz, _tz, order * order * sizeof(TacsScalar));

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
    memcpy(box, _box, 6 * sizeof(TacsScalar));
  }
  TACS3DBoundingTraction(int _surface, TacsScalar _tx, TacsScalar _ty,
                         TacsScalar _tz, TacsScalar _box[]) {
    surface = _surface;
    for (int i = 0; i < order * order; i++) {
      tx[i] = _tx;
      ty[i] = _ty;
      tz[i] = _tz;
    }

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
    memcpy(box, _box, 6 * sizeof(TacsScalar));
  }

  // Return the number of displacements, stresses and nodes
  // ------------------------------------------------------
  int numDisplacements() { return 3; }
  int numNodes() { return NUM_NODES; }
  int numStresses() { return 0; }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies(double time, TacsScalar *Te, TacsScalar *Pe,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[]) {
    *Te = 0.0, *Pe = 0.0;
  }

  /*
    Compute the shape functions and their derivatives w.r.t. the
    parametric element location
  */
  void getShapeFunctions(const double pt[], double N[], double Na[],
                         double Nb[]) {
    double na[order], nb[order];
    double dna[order], dnb[order];
    FElibrary::lagrangeSFKnots(na, dna, pt[0], knots, order);
    FElibrary::lagrangeSFKnots(nb, dnb, pt[1], knots, order);
    for (int j = 0; j < order; j++) {
      for (int i = 0; i < order; i++) {
        N[0] = na[i] * nb[j];
        Na[0] = dna[i] * nb[j];
        Nb[0] = na[i] * dnb[j];
        N++;
        Na++;
        Nb++;
      }
    }
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
        double pt[3];
        pt[0] = gaussPts[n];
        pt[1] = gaussPts[m];

        // Compute X, Xd, N, Na and Nb
        double N[order * order];
        double Na[order * order], Nb[order * order];
        getShapeFunctions(pt, N, Na, Nb);

        TacsScalar Xa[3], Xb[3];
        Xa[0] = Xa[1] = Xa[2] = 0.0;
        Xb[0] = Xb[1] = Xb[2] = 0.0;

        int box_x = 0, box_y = 0, box_z = 0;

        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[j + k * order] * Xpts[3 * node];
              Xa[1] += Na[j + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[j + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[j + k * order] * Xpts[3 * node];
              Xb[1] += Nb[j + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[j + k * order] * Xpts[3 * node + 2];

              if ((TacsRealPart(Xpts[3 * node]) <= TacsRealPart(box[3])) &&
                  (TacsRealPart(Xpts[3 * node]) >= TacsRealPart(box[0]))) {
                box_x = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 1]) <= TacsRealPart(box[4])) &&
                  (TacsRealPart(Xpts[3 * node + 1]) >= TacsRealPart(box[1]))) {
                box_y = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 2]) <= TacsRealPart(box[5])) &&
                  (TacsRealPart(Xpts[3 * node + 2]) >= TacsRealPart(box[2]))) {
                box_z = 1;
              }
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + k * order] * Xpts[3 * node];
              Xa[1] += Na[i + k * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + k * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + k * order] * Xpts[3 * node];
              Xb[1] += Nb[i + k * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + k * order] * Xpts[3 * node + 2];

              if ((TacsRealPart(Xpts[3 * node]) <= TacsRealPart(box[3])) &&
                  (TacsRealPart(Xpts[3 * node]) >= TacsRealPart(box[0]))) {
                box_x = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 1]) <= TacsRealPart(box[4])) &&
                  (TacsRealPart(Xpts[3 * node + 1]) >= TacsRealPart(box[1]))) {
                box_y = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 2]) <= TacsRealPart(box[5])) &&
                  (TacsRealPart(Xpts[3 * node + 2]) >= TacsRealPart(box[2]))) {
                box_z = 1;
              }
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;

              Xa[0] += Na[i + j * order] * Xpts[3 * node];
              Xa[1] += Na[i + j * order] * Xpts[3 * node + 1];
              Xa[2] += Na[i + j * order] * Xpts[3 * node + 2];

              Xb[0] += Nb[i + j * order] * Xpts[3 * node];
              Xb[1] += Nb[i + j * order] * Xpts[3 * node + 1];
              Xb[2] += Nb[i + j * order] * Xpts[3 * node + 2];

              if ((TacsRealPart(Xpts[3 * node]) <= TacsRealPart(box[3])) &&
                  (TacsRealPart(Xpts[3 * node]) >= TacsRealPart(box[0]))) {
                box_x = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 1]) <= TacsRealPart(box[4])) &&
                  (TacsRealPart(Xpts[3 * node + 1]) >= TacsRealPart(box[1]))) {
                box_y = 1;
              }
              if ((TacsRealPart(Xpts[3 * node + 2]) <= TacsRealPart(box[5])) &&
                  (TacsRealPart(Xpts[3 * node + 2]) >= TacsRealPart(box[2]))) {
                box_z = 1;
              }
            }
          }
        }

        // Determine the normal direction
        TacsScalar normal[3];
        Tensor::crossProduct3D(normal, Xa, Xb);
        TacsScalar h = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                            normal[2] * normal[2]);
        h *= gaussWts[n] * gaussWts[m];

        // Evaluate the traction force evaluated at the quadrature point within
        // the element if it is within the bounding box
        TacsScalar Tx = 0.0, Ty = 0.0, Tz = 0.0;

        if (box_x && box_y && box_z) {
          for (int i = 0; i < order * order; i++) {
            Tx += tx[i] * N[i];
            Ty += ty[i] * N[i];
            Tz += tz[i] * N[i];
          }
        }
        // Add the contribution to the residual - the minus sign
        // is due to the fact that this is a work term
        if (surface < 2) {
          const int i = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int j = 0; j < order; j++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[j + k * order];
              res[3 * node + 1] -= h * Ty * N[j + k * order];
              res[3 * node + 2] -= h * Tz * N[j + k * order];
            }
          }
        } else if (surface < 4) {
          const int j = (order - 1) * (surface % 2);
          for (int k = 0; k < order; k++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + k * order];
              res[3 * node + 1] -= h * Ty * N[i + k * order];
              res[3 * node + 2] -= h * Tz * N[i + k * order];
            }
          }
        } else {
          const int k = (order - 1) * (surface % 2);
          for (int j = 0; j < order; j++) {
            for (int i = 0; i < order; i++) {
              const int node = i + j * order + k * order * order;
              res[3 * node] -= h * Tx * N[i + j * order];
              res[3 * node + 1] -= h * Ty * N[i + j * order];
              res[3 * node + 2] -= h * Tz * N[i + j * order];
            }
          }
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

 private:
  int surface;
  TacsScalar tx[order * order];
  TacsScalar ty[order * order];
  TacsScalar tz[order * order];

  // The knot locations for the basis functions
  double knots[order];
  TacsScalar box[6];
};

#endif  // TACS_3D_TRACTION_H
