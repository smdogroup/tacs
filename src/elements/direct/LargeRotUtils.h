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

#ifndef LARGE_ROT_UTILS_H
#define LARGE_ROT_UTILS_H

#include "ShellUtils.h"
#include "TACSObject.h"
#include "TensorToolbox.h"

/*
  Common shell related computations for large rotations/moderate
  strains.
*/

TACS_BEGIN_NAMESPACE(largerot)

using namespace shellutils;

/*
  Compute the matrix that is used to evaluate the rate of change
  of displacements through the thickness. This matrix can also
  be used to determine the rotation of the normal in the deformed
  configuration. C is defined as:

  C = Q^{T} - I

  where Q is an orthonormal rotation matrix and I is the identity
  matrix. Here Q is defined using Euler angles. In most cases the
  singularities in Q will not be reached unless extreme deformation
  of the structure occurs.

  Note that due to the formulation, up to third-derivatives of C
  are required!
*/
void compute_rate_matrix(TacsScalar C[], TacsScalar c1, TacsScalar s1,
                         TacsScalar c2, TacsScalar s2, TacsScalar c3,
                         TacsScalar s3);

void compute_rate_matrix(TacsScalar C[], TacsScalar Ct[], TacsScalar c1,
                         TacsScalar s1, TacsScalar c2, TacsScalar s2,
                         TacsScalar c3, TacsScalar s3);

void compute_2nd_rate_matrix(TacsScalar Ctt[], TacsScalar c1, TacsScalar s1,
                             TacsScalar c2, TacsScalar s2, TacsScalar c3,
                             TacsScalar s3);

void compute_3rd_rate_matrix(TacsScalar Cttt[], TacsScalar c1, TacsScalar s1,
                             TacsScalar c2, TacsScalar s2, TacsScalar c3,
                             TacsScalar s3);

/*
  This function tests the implementation of the rate matrix
  and its derivatives. This uses a finite-difference step
  of dh and prints out all matrix components to stderr.
*/
void test_rate_matrix(TacsScalar U[], double dh);

/*
  Compute the in-plane penalty term that penalizes the difference
  between the drilling rotation degree of freedom and the in-plane
  rotations defined by the deformation tensor. This technique is based
  on the polar decomposition of the deformation tensor due to Simo et
  al. - a very clever trick!

  This code computes the in-plane penalty and its derivative at a
  point.
*/
TacsScalar compute_inplane_penalty(TacsScalar drot[], const int num_points,
                                   const TacsScalar Xd[], const TacsScalar Ud[],
                                   const TacsScalar C[], const TacsScalar Ct[],
                                   const double N[], const double Na[],
                                   const double Nb[]);

/*
  This code adds the second derivative of the in-plane penalty term to
  the lower portion of the stiffness matrix. Note that this
  computation only includes terms from the rotation matrices since the
  displacements only appear as a quadratic term handled through the
  mat[i,j] = drot[i]*drot[j] penalty term in the stiffness matrix
  computation.
*/
void add_inplane_penalty(TacsScalar matrix[], const int num_points,
                         const TacsScalar scale, const TacsScalar Xd[],
                         const TacsScalar Ud[], const TacsScalar Ct[],
                         const TacsScalar Ctt[], const double N[],
                         const double Na[], const double Nb[]);

/*
  This code evaluates the strain based on the displacements and
  rotations, the rates of change of the displacements and rotations
  along the parametric coordinate directions as well as the local
  coordinate transformation and information about the normal
  direction.
*/
void large_rot_strain(TacsScalar strain[], const TacsScalar Ux[],
                      const TacsScalar Uxd[], const TacsScalar C[],
                      const TacsScalar Ct[], const TacsScalar t[],
                      const TacsScalar tx[], const TacsScalar ztx[],
                      const TacsScalar n[], const TacsScalar n_xi[],
                      const TacsScalar n_eta[]);

/*
  Evaluate the strain and the derivative of the strain using the
  large-rotation formula. Note that the input derivativess are stored
  in essentially column-major format - they are stacked on one
  another. This is the convention throughout TACS for these
  derivatives!
*/
void large_rot_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                           const TacsScalar Ux[], const TacsScalar Uxd[],
                           const TacsScalar C[], const TacsScalar Ct[],
                           const TacsScalar t[], const TacsScalar dt[],
                           const TacsScalar tx[], const TacsScalar dtx[],
                           const TacsScalar ztx[], const TacsScalar dztx[],
                           const TacsScalar n[], const TacsScalar dn[],
                           const TacsScalar n_xi[], const TacsScalar dn_xi[],
                           const TacsScalar n_eta[], const TacsScalar dn_eta[],
                           const int num_components);

/*
  This code evaluates the derivative of the strain expressions
  w.r.t. the nodal displacements and rotations and returns the result
  in the matrix B. Note that the first and second derivatives of the
  rotation matrix are required.
*/
void large_rot_bmat(TacsScalar B[], const int num_points, const double N[],
                    const double Na[], const double Nb[], const TacsScalar Ux[],
                    const TacsScalar Uxd[], const TacsScalar C[],
                    const TacsScalar Ct[], const TacsScalar Ctt[],
                    const TacsScalar t[], const TacsScalar tx[],
                    const TacsScalar ztx[], const TacsScalar n[],
                    const TacsScalar n_xi[], const TacsScalar n_eta[]);

/*
  Add the derivative of product the strain times a stress vector to the
  residual vector. In other words, compute d(stress^{T}B)/dx.

  This is required for computing the derivative of the element residual
  w.r.t. the nodes.
*/
void add_large_rot_bmat_sens(
    TacsScalar res[], const int num_points, const TacsScalar scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar C[], const TacsScalar Ct[], const TacsScalar Ctt[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], const int num_components);

/*
  Add the contributions from the inner product of the stress with the
  second derivative of the strain. This is more costly than small-angle
  approximations b/c the second derivative is most costly to compute in
  this case. Note that up to third derivatives of the rotation matrix
  are required.
*/
void add_large_rot_stress_bmat(
    TacsScalar matrix[], const int num_points, const TacsScalar scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar C[], const TacsScalar Ct[], const TacsScalar Ctt[],
    const TacsScalar Cttt[], const TacsScalar t[], const TacsScalar tx[],
    const TacsScalar ztx[], const TacsScalar n[], const TacsScalar n_xi[],
    const TacsScalar n_eta[]);

/*
  The following functions are used to compute the bending strain
  contributions to the terms computed above - e.g. strain, B-matrix,
  geometric stiffness, etc. These contributions are required for the
  MITC-based formulation of the large-rotation shell element.
*/

/*
  This code evaluates the strain based on the displacements and
  rotations, the rates of change of the displacements and rotations
  along the parametric coordinate directions as well as the local
  coordinate transformation and information about the normal
  direction. Unlike the function above, only the bending contributions
  to the strain are evaluated.
*/
void large_rot_bend_strain(TacsScalar strain[], const TacsScalar Ux[],
                           const TacsScalar Uxd[], const TacsScalar C[],
                           const TacsScalar Ct[], const TacsScalar t[],
                           const TacsScalar tx[], const TacsScalar ztx[],
                           const TacsScalar n[], const TacsScalar n_xi[],
                           const TacsScalar n_eta[]);

/*
  Evaluate the strain and the derivative of the strain using the
  large-rotation formula. Only the bending contributions are computed.
  Note that the input derivativess are stored in essentially
  column-major format - they are stacked on one another. This is the
  convention throughout TACS for these derivatives!
*/
void large_rot_bend_strain_sens(
    TacsScalar strain[], TacsScalar dstrain[], const TacsScalar Ux[],
    const TacsScalar Uxd[], const TacsScalar C[], const TacsScalar Ct[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], const int num_components);

/*
  This code evaluates the derivative of the strain expressions
  w.r.t. the nodal displacements and rotations and returns the result
  in the matrix B. Only the bending contributions to B are included.
  Note that the first and second derivatives of the rotation matrix
  are required.
*/
void large_rot_bend_bmat(TacsScalar B[], const int num_points, const double N[],
                         const double Na[], const double Nb[],
                         const TacsScalar Ux[], const TacsScalar Uxd[],
                         const TacsScalar C[], const TacsScalar Ct[],
                         const TacsScalar Ctt[], const TacsScalar t[],
                         const TacsScalar tx[], const TacsScalar ztx[],
                         const TacsScalar n[], const TacsScalar n_xi[],
                         const TacsScalar n_eta[]);

/*
  Add the contributions from the inner product of the stress with the
  second derivative of the strain. This is more costly than small-angle
  approximations b/c the second derivative is most costly to compute in
  this case. Note that up to third derivatives of the rotation matrix
  are required.
*/
void add_large_rot_bend_stress_bmat(
    TacsScalar matrix[], const int num_points, const TacsScalar scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar C[], const TacsScalar Ct[], const TacsScalar Ctt[],
    const TacsScalar Cttt[], const TacsScalar t[], const TacsScalar tx[],
    const TacsScalar ztx[], const TacsScalar n[], const TacsScalar n_xi[],
    const TacsScalar n_eta[]);

/*
  The following template-based functions are required for the
  nonlinear implementation of the MITC shell element with large
  rotations. These functions are used to compute the contributions to
  the strain from the strain interpolation.

  Note that this implementation is very similar to the linear
  implementation in ShellUtils.h
*/

/*
  Evaluate the tensorial displacement-based strain at the tying points.

  The tensorial strain is given as follows,

  E_t = U_{x,xi}*X_{,xi}^{T} + X_{,xi}*U_{x,xi}^{T}

  X_{,xi} = Xd^{T}
  U_{x,xi} = [ Ud, r ]

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points

  vars: the displacements and rotations at the nodal locations
  Xpts: the nodal locations for the element

  output:
  g11: the in-plane normal tensorial strain in the 1 direction
  g22: the in-plane normal tensorial strain in the 2 direction
  g12: the in-plane tensorial shear strain

  g23: the out-of-plane tensorial shear strain in the 2-direction
  g13: the out-of-plane tensorial shear strain in the 1-direction
*/
template <int order, int tying_order>
void compute_lr_tying_strain(TacsScalar g11[], TacsScalar g22[],
                             TacsScalar g12[], TacsScalar g23[],
                             TacsScalar g13[], const double knots[],
                             const double pknots[], const TacsScalar vars[],
                             const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], r[3];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
      g23[0] = 0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                      (normal[0] + r[0]) * Ud[1] + (normal[1] + r[1]) * Ud[3] +
                      (normal[2] + r[2]) * Ud[5]);
      g22++;
      g23++;
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                      Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                      Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
      g12++;
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], r[3];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
      g13[0] = 0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                      (normal[0] + r[0]) * Ud[0] + (normal[1] + r[1]) * Ud[2] +
                      (normal[2] + r[2]) * Ud[4]);
      g11++;
      g13++;
    }
  }
}

/*
  Evaluate the strain at the tying points and the sensitivity of the
  strain at the tying points w.r.t. the nodal locations.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  g11, g22, g12, g23, g13: the in-plane tensorial strains and
  out-of-plane tensorial shear strains at the tying points

  dg11, dg22, dg12, dg23, dg13: the sensitivities of the the in-plane
  tensorial strains and out-of-plane tensorial shear strains at the
  tying points w.r.t. the nodal locataions
*/
template <int order, int tying_order>
void compute_lr_tying_strain_sens(
    TacsScalar g11[], TacsScalar g22[], TacsScalar g12[], TacsScalar g23[],
    TacsScalar g13[], TacsScalar dg11[], TacsScalar dg22[], TacsScalar dg12[],
    TacsScalar dg23[], TacsScalar dg13[], const double knots[],
    const double pknots[], const TacsScalar vars[], const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      // Calculate the rotation rate matrix
      TacsScalar C[9], r[3];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      // Evaluate the tensorial strain
      g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
      g23[0] = 0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                      (normal[0] + r[0]) * Ud[1] + (normal[1] + r[1]) * Ud[3] +
                      (normal[2] + r[2]) * Ud[5]);
      g22++;
      g23++;

      TacsScalar Xt[3], Ut[3];
      Xt[0] = C[0] * Xd[3] + C[3] * Xd[4] + C[6] * Xd[5];
      Xt[1] = C[1] * Xd[3] + C[4] * Xd[4] + C[7] * Xd[5];
      Xt[2] = C[2] * Xd[3] + C[5] * Xd[4] + C[8] * Xd[5];

      Ut[0] = C[0] * Ud[1] + C[3] * Ud[3] + C[6] * Ud[5];
      Ut[1] = C[1] * Ud[1] + C[4] * Ud[3] + C[7] * Ud[5];
      Ut[2] = C[2] * Ud[1] + C[5] * Ud[3] + C[8] * Ud[5];

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          TacsScalar normalSens[3];

          if (k == 0) {
            normalSens[0] = 0.0;
            normalSens[1] = Xd[2] * Nb[i] - Na[i] * Xd[5];
            normalSens[2] = Na[i] * Xd[4] - Xd[1] * Nb[i];
          } else if (k == 1) {
            normalSens[0] = Na[i] * Xd[5] - Xd[2] * Nb[i];
            normalSens[1] = 0.0;
            normalSens[2] = Xd[0] * Nb[i] - Na[i] * Xd[3];
          } else {
            normalSens[0] = Xd[1] * Nb[i] - Na[i] * Xd[4];
            normalSens[1] = Na[i] * Xd[3] - Xd[0] * Nb[i];
            normalSens[2] = 0.0;
          }

          TacsScalar normSens =
              (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
               normal[2] * normalSens[2]);

          normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
          normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
          normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

          // Calculate the rotation
          TacsScalar rSens[3];
          Tensor::crossProduct3D(rSens, Urot, normalSens);

          dg22[0] = Ud[2 * k + 1] * Nb[i];
          dg23[0] =
              0.5 * (r[k] * Nb[i] + normalSens[0] * (Ud[1] + Xt[0] + Ut[0]) +
                     normalSens[1] * (Ud[3] + Xt[1] + Ut[1]) +
                     normalSens[2] * (Ud[5] + Xt[2] + Ut[2]));
          dg22++;
          dg23++;
        }
      }
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                      Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                      Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
      g12++;

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          dg12[0] = 0.5 * (Nb[i] * Ud[2 * k] + Na[i] * Ud[2 * k + 1]);
          dg12++;
        }
      }
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      TacsScalar inv_norm = Tensor::normalize3D(normal);
      inv_norm = 1.0 / inv_norm;

      // Calculate the rotation rate matrix
      TacsScalar C[9], r[3];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
      g13[0] = 0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                      (normal[0] + r[0]) * Ud[0] + (normal[1] + r[1]) * Ud[2] +
                      (normal[2] + r[2]) * Ud[4]);
      g11++;
      g13++;

      TacsScalar Xt[3], Ut[3];
      Xt[0] = C[0] * Xd[0] + C[3] * Xd[1] + C[6] * Xd[2];
      Xt[1] = C[1] * Xd[0] + C[4] * Xd[1] + C[7] * Xd[2];
      Xt[2] = C[2] * Xd[0] + C[5] * Xd[1] + C[8] * Xd[2];

      Ut[0] = C[0] * Ud[0] + C[3] * Ud[2] + C[6] * Ud[4];
      Ut[1] = C[1] * Ud[0] + C[4] * Ud[2] + C[7] * Ud[4];
      Ut[2] = C[2] * Ud[0] + C[5] * Ud[2] + C[8] * Ud[4];

      for (int i = 0; i < order * order; i++) {
        for (int k = 0; k < 3; k++) {
          TacsScalar normalSens[3];

          if (k == 0) {
            normalSens[0] = 0.0;
            normalSens[1] = Xd[2] * Nb[i] - Na[i] * Xd[5];
            normalSens[2] = Na[i] * Xd[4] - Xd[1] * Nb[i];
          } else if (k == 1) {
            normalSens[0] = Na[i] * Xd[5] - Xd[2] * Nb[i];
            normalSens[1] = 0.0;
            normalSens[2] = Xd[0] * Nb[i] - Na[i] * Xd[3];
          } else {
            normalSens[0] = Xd[1] * Nb[i] - Na[i] * Xd[4];
            normalSens[1] = Na[i] * Xd[3] - Xd[0] * Nb[i];
            normalSens[2] = 0.0;
          }

          TacsScalar normSens =
              (normal[0] * normalSens[0] + normal[1] * normalSens[1] +
               normal[2] * normalSens[2]);

          normalSens[0] = (normalSens[0] - normSens * normal[0]) * inv_norm;
          normalSens[1] = (normalSens[1] - normSens * normal[1]) * inv_norm;
          normalSens[2] = (normalSens[2] - normSens * normal[2]) * inv_norm;

          dg11[0] = Ud[2 * k] * Na[i];
          dg13[0] =
              0.5 * (Na[i] * r[k] + normalSens[0] * (Ud[0] + Xt[0] + Ut[0]) +
                     normalSens[1] * (Ud[2] + Xt[1] + Ut[1]) +
                     normalSens[2] * (Ud[4] + Xt[2] + Ut[2]));
          dg11++;
          dg13++;
        }
      }
    }
  }
}

/*
  Compute the tying bmat matrices. This corresponds to the derivative
  of the displacement-based strain at each of the tying points.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  b11: the derivative of the in-plane normal tensorial strain in the 1
  direction
  b22: the derivative of the in-plane normal tensorial strain in the 2
  direction
  b12: the derivative of the in-plane tensorial shear strain

  b23: the derivative of the out-of-plane tensorial shear strain in
  the 2-direction
  b13: the derivative of the out-of-plane tensorial shear strain in
  the 1-direction
*/
template <int order, int tying_order>
void compute_lr_tying_bmat(TacsScalar g11[], TacsScalar g22[], TacsScalar g12[],
                           TacsScalar g23[], TacsScalar g13[], TacsScalar b11[],
                           TacsScalar b22[], TacsScalar b12[], TacsScalar b23[],
                           TacsScalar b13[], const double knots[],
                           const double pknots[], const TacsScalar vars[],
                           const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], Ct[3 * 9], r[3], rt[9];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      // Compute the derivative of C w.r.t. each rotation, times the normal
      // vector - store the result in the array rt[]
      for (int k = 0; k < 3; k++) {
        const TacsScalar* Ca = &Ct[9 * k];
        rt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Evaluate the tensorial strain
      g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
      g23[0] = 0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                      (normal[0] + r[0]) * Ud[1] + (normal[1] + r[1]) * Ud[3] +
                      (normal[2] + r[2]) * Ud[5]);
      g22++;
      g23++;

      // Add the contributions to the b-matrix
      for (int i = 0; i < order * order; i++) {
        b22[0] = (Xd[3] + Ud[1]) * Nb[i];
        b22[1] = (Xd[4] + Ud[3]) * Nb[i];
        b22[2] = (Xd[5] + Ud[5]) * Nb[i];
        b22 += 3;

        // Compute the derivative of g23 w.r.t. the displacements and
        // rotations
        b23[0] = 0.5 * (normal[0] + r[0]) * Nb[i];
        b23[1] = 0.5 * (normal[1] + r[1]) * Nb[i];
        b23[2] = 0.5 * (normal[2] + r[2]) * Nb[i];
        b23[3] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[0] + (Xd[4] + Ud[3]) * rt[1] +
                  (Xd[5] + Ud[5]) * rt[2]);
        b23[4] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[3] + (Xd[4] + Ud[3]) * rt[4] +
                  (Xd[5] + Ud[5]) * rt[5]);
        b23[5] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[6] + (Xd[4] + Ud[3]) * rt[7] +
                  (Xd[5] + Ud[5]) * rt[8]);
        b23 += 6;
      }
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                      Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                      Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
      g12++;

      for (int i = 0; i < order * order; i++) {
        b12[0] = 0.5 * ((Xd[3] + Ud[1]) * Na[i] + (Xd[0] + Ud[0]) * Nb[i]);
        b12[1] = 0.5 * ((Xd[4] + Ud[3]) * Na[i] + (Xd[1] + Ud[2]) * Nb[i]);
        b12[2] = 0.5 * ((Xd[5] + Ud[5]) * Na[i] + (Xd[2] + Ud[4]) * Nb[i]);
        b12 += 3;
      }
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], Ct[3 * 9], r[3], rt[9];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      // Compute the derivative of C w.r.t. each rotation, times the normal
      // vector - store the result in the array rt[]
      for (int k = 0; k < 3; k++) {
        const TacsScalar* Ca = &Ct[9 * k];
        rt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Compute the tensorial strain components
      g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
      g13[0] = 0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                      (normal[0] + r[0]) * Ud[0] + (normal[1] + r[1]) * Ud[2] +
                      (normal[2] + r[2]) * Ud[4]);
      g11++;
      g13++;

      // Compute the derivatives of the tensorial strain w.r.t. the
      // displacements and rotations
      for (int i = 0; i < order * order; i++) {
        b11[0] = (Xd[0] + Ud[0]) * Na[i];
        b11[1] = (Xd[1] + Ud[2]) * Na[i];
        b11[2] = (Xd[2] + Ud[4]) * Na[i];
        b11 += 3;

        // Compute the derivative of the shear strain components
        b13[0] = 0.5 * (normal[0] + r[0]) * Na[i];
        b13[1] = 0.5 * (normal[1] + r[1]) * Na[i];
        b13[2] = 0.5 * (normal[2] + r[2]) * Na[i];
        b13[3] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[0] + (Xd[1] + Ud[2]) * rt[1] +
                  (Xd[2] + Ud[4]) * rt[2]);
        b13[4] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[3] + (Xd[1] + Ud[2]) * rt[4] +
                  (Xd[2] + Ud[4]) * rt[5]);
        b13[5] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[6] + (Xd[1] + Ud[2]) * rt[7] +
                  (Xd[2] + Ud[4]) * rt[8]);
        b13 += 6;
      }
    }
  }
}

/*
  Compute the strain, the derivative of the strain and the second
  derivatives of the strain at the tying points in the matrix. This
  corresponds to the derivative of the displacement-based strain at
  each of the tying points.

  input:
  knots: the 'order'th (ie 2, 3 or 4) Gauss-points
  pknots: the 'order-1'th Gauss-points
  Xpts: the nodal locations for the element

  output:
  g11, g22, g12, g23, g13: the tensorial strain at the tying points
  in the quadrilateral element

  b11, b22, b12, b23, b13: the derivatives of the tensorial strain
  at the tying points w.r.t. the nodal displacements/rotations

  n23, n13: the second derivatives of the tensorial strain at the
  tying points w.r.t. the rotationals

  Note that the size of the n23/n13 matrices are:

  3*(order*order*(order*order+1))*NUM_G23, and
  3*(order*order*(order*order+1))*NUM_G13
*/
template <int order, int tying_order>
void compute_lr_tying_nmat(TacsScalar g11[], TacsScalar g22[], TacsScalar g12[],
                           TacsScalar g23[], TacsScalar g13[], TacsScalar b11[],
                           TacsScalar b22[], TacsScalar b12[], TacsScalar b23[],
                           TacsScalar b13[], TacsScalar n23[], TacsScalar n13[],
                           const double knots[], const double pknots[],
                           const TacsScalar vars[], const TacsScalar Xpts[]) {
  double N[order * order], Na[order * order], Nb[order * order];

  // Evaluate g22 and g23
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order; n++) {
      double pt[2];
      pt[0] = knots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], Ct[3 * 9], Ctt[6 * 9], r[3], rt[9], rtt[18];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
      compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      // Compute the derivative of C w.r.t. each rotation, times the normal
      // vector - store the result in the array rt[]
      for (int k = 0; k < 3; k++) {
        const TacsScalar* Ca = &Ct[9 * k];
        rt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Compute the product of the second derivatives w.r.t. the rotations
      // Note that the order
      for (int k = 0; k < 6; k++) {
        const TacsScalar* Ca = &Ctt[9 * k];
        rtt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rtt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rtt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Evaluate the tensorial strain
      g22[0] = (Xd[3] * Ud[1] + Xd[4] * Ud[3] + Xd[5] * Ud[5] +
                0.5 * (Ud[1] * Ud[1] + Ud[3] * Ud[3] + Ud[5] * Ud[5]));
      g23[0] = 0.5 * (Xd[3] * r[0] + Xd[4] * r[1] + Xd[5] * r[2] +
                      (normal[0] + r[0]) * Ud[1] + (normal[1] + r[1]) * Ud[3] +
                      (normal[2] + r[2]) * Ud[5]);
      g22++;
      g23++;

      // Add the contributions to the b-matrix
      for (int i = 0; i < order * order; i++) {
        b22[0] = (Xd[3] + Ud[1]) * Nb[i];
        b22[1] = (Xd[4] + Ud[3]) * Nb[i];
        b22[2] = (Xd[5] + Ud[5]) * Nb[i];
        b22 += 3;

        // Compute the derivative of g23 w.r.t. the displacements and
        // rotations
        b23[0] = 0.5 * (normal[0] + r[0]) * Nb[i];
        b23[1] = 0.5 * (normal[1] + r[1]) * Nb[i];
        b23[2] = 0.5 * (normal[2] + r[2]) * Nb[i];
        b23[3] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[0] + (Xd[4] + Ud[3]) * rt[1] +
                  (Xd[5] + Ud[5]) * rt[2]);
        b23[4] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[3] + (Xd[4] + Ud[3]) * rt[4] +
                  (Xd[5] + Ud[5]) * rt[5]);
        b23[5] = 0.5 * N[i] *
                 ((Xd[3] + Ud[1]) * rt[6] + (Xd[4] + Ud[3]) * rt[7] +
                  (Xd[5] + Ud[5]) * rt[8]);
        b23 += 6;

        // Evaluate the terms that represent the second derivative g23 w.r.t.
        // the rotation degrees of freedom
        for (int j = 0; j <= i; j++) {
          // Compute the second derivative: d(b23[0:3])/d(theta)
          n23[0] = 0.5 * rt[0] * Nb[i] * N[j];
          n23[1] = 0.5 * rt[3] * Nb[i] * N[j];
          n23[2] = 0.5 * rt[6] * Nb[i] * N[j];
          n23[3] = 0.5 * rt[1] * Nb[i] * N[j];
          n23[4] = 0.5 * rt[4] * Nb[i] * N[j];
          n23[5] = 0.5 * rt[7] * Nb[i] * N[j];
          n23[6] = 0.5 * rt[2] * Nb[i] * N[j];
          n23[7] = 0.5 * rt[5] * Nb[i] * N[j];
          n23[8] = 0.5 * rt[8] * Nb[i] * N[j];
          n23 += 9;

          // Compute the second derivative: d(b23[3:6])/du
          n23[0] = 0.5 * rt[0] * N[i] * Nb[j];
          n23[1] = 0.5 * rt[1] * N[i] * Nb[j];
          n23[2] = 0.5 * rt[2] * N[i] * Nb[j];
          n23[3] = 0.5 * rt[3] * N[i] * Nb[j];
          n23[4] = 0.5 * rt[4] * N[i] * Nb[j];
          n23[5] = 0.5 * rt[5] * N[i] * Nb[j];
          n23[6] = 0.5 * rt[6] * N[i] * Nb[j];
          n23[7] = 0.5 * rt[7] * N[i] * Nb[j];
          n23[8] = 0.5 * rt[8] * N[i] * Nb[j];
          n23 += 9;

          // Compute the symmetric derivative of d(b23[3:6])/d(theta)
          n23[0] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[0] + (Xd[4] + Ud[3]) * rtt[1] +
                    (Xd[5] + Ud[5]) * rtt[2]);
          n23[1] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[3] + (Xd[4] + Ud[3]) * rtt[4] +
                    (Xd[5] + Ud[5]) * rtt[5]);
          n23[2] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[6] + (Xd[4] + Ud[3]) * rtt[7] +
                    (Xd[5] + Ud[5]) * rtt[8]);
          n23[3] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[9] + (Xd[4] + Ud[3]) * rtt[10] +
                    (Xd[5] + Ud[5]) * rtt[11]);
          n23[4] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[12] + (Xd[4] + Ud[3]) * rtt[13] +
                    (Xd[5] + Ud[5]) * rtt[14]);
          n23[5] = 0.5 * N[j] * N[i] *
                   ((Xd[3] + Ud[1]) * rtt[15] + (Xd[4] + Ud[3]) * rtt[16] +
                    (Xd[5] + Ud[5]) * rtt[17]);
          n23 += 6;
        }
      }
    }
  }

  // Evaluate g12
  for (int m = 0; m < tying_order - 1; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = pknots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      g12[0] = 0.5 * (Xd[0] * Ud[1] + Xd[1] * Ud[3] + Xd[2] * Ud[5] +
                      Xd[3] * Ud[0] + Xd[4] * Ud[2] + Xd[5] * Ud[4] +
                      Ud[0] * Ud[1] + Ud[2] * Ud[3] + Ud[4] * Ud[5]);
      g12++;

      for (int i = 0; i < order * order; i++) {
        b12[0] = 0.5 * ((Xd[3] + Ud[1]) * Na[i] + (Xd[0] + Ud[0]) * Nb[i]);
        b12[1] = 0.5 * ((Xd[4] + Ud[3]) * Na[i] + (Xd[1] + Ud[2]) * Nb[i]);
        b12[2] = 0.5 * ((Xd[5] + Ud[5]) * Na[i] + (Xd[2] + Ud[4]) * Nb[i]);
        b12 += 3;
      }
    }
  }

  // Evaluate g11 and g13
  for (int m = 0; m < tying_order; m++) {
    for (int n = 0; n < tying_order - 1; n++) {
      double pt[2];
      pt[0] = pknots[n];
      pt[1] = knots[m];

      TacsScalar Urot[3], Ud[6], Xd[6];
      compute_tensorial_components(order, Urot, Ud, Xd, N, Na, Nb, pt, Xpts,
                                   vars);

      // Caculate the normal
      TacsScalar normal[3];
      Tensor::crossProduct3D(normal, &Xd[0], &Xd[3]);
      Tensor::normalize3D(normal);

      // Calculate the rotation rate matrix
      TacsScalar C[9], Ct[3 * 9], Ctt[6 * 9], r[3], rt[9], rtt[18];
      TacsScalar c1 = cos(Urot[0]), s1 = sin(Urot[0]);
      TacsScalar c2 = cos(Urot[1]), s2 = sin(Urot[1]);
      TacsScalar c3 = cos(Urot[2]), s3 = sin(Urot[2]);
      compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
      compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);

      // Evaluate the through-thickness rate of change of displacement
      r[0] = C[0] * normal[0] + C[1] * normal[1] + C[2] * normal[2];
      r[1] = C[3] * normal[0] + C[4] * normal[1] + C[5] * normal[2];
      r[2] = C[6] * normal[0] + C[7] * normal[1] + C[8] * normal[2];

      // Compute the derivative of C w.r.t. each rotation, times the normal
      // vector - store the result in the array rt[]
      for (int k = 0; k < 3; k++) {
        const TacsScalar* Ca = &Ct[9 * k];
        rt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Compute the product of the second derivatives w.r.t. the rotations
      // Note that the order
      for (int k = 0; k < 6; k++) {
        const TacsScalar* Ca = &Ctt[9 * k];
        rtt[3 * k] = Ca[0] * normal[0] + Ca[1] * normal[1] + Ca[2] * normal[2];
        rtt[3 * k + 1] =
            Ca[3] * normal[0] + Ca[4] * normal[1] + Ca[5] * normal[2];
        rtt[3 * k + 2] =
            Ca[6] * normal[0] + Ca[7] * normal[1] + Ca[8] * normal[2];
      }

      // Compute the tensorial strain components
      g11[0] = (Xd[0] * Ud[0] + Xd[1] * Ud[2] + Xd[2] * Ud[4] +
                0.5 * (Ud[0] * Ud[0] + Ud[2] * Ud[2] + Ud[4] * Ud[4]));
      g13[0] = 0.5 * (Xd[0] * r[0] + Xd[1] * r[1] + Xd[2] * r[2] +
                      (normal[0] + r[0]) * Ud[0] + (normal[1] + r[1]) * Ud[2] +
                      (normal[2] + r[2]) * Ud[4]);
      g11++;
      g13++;

      // Compute the derivatives of the tensorial strain w.r.t. the
      // displacements and rotations
      for (int i = 0; i < order * order; i++) {
        b11[0] = (Xd[0] + Ud[0]) * Na[i];
        b11[1] = (Xd[1] + Ud[2]) * Na[i];
        b11[2] = (Xd[2] + Ud[4]) * Na[i];
        b11 += 3;

        // Compute the derivative of the shear strain components
        b13[0] = 0.5 * (normal[0] + r[0]) * Na[i];
        b13[1] = 0.5 * (normal[1] + r[1]) * Na[i];
        b13[2] = 0.5 * (normal[2] + r[2]) * Na[i];
        b13[3] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[0] + (Xd[1] + Ud[2]) * rt[1] +
                  (Xd[2] + Ud[4]) * rt[2]);
        b13[4] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[3] + (Xd[1] + Ud[2]) * rt[4] +
                  (Xd[2] + Ud[4]) * rt[5]);
        b13[5] = 0.5 * N[i] *
                 ((Xd[0] + Ud[0]) * rt[6] + (Xd[1] + Ud[2]) * rt[7] +
                  (Xd[2] + Ud[4]) * rt[8]);
        b13 += 6;

        // Evaluate the terms that represent the second derivative g23 w.r.t.
        // the rotation degrees of freedom
        for (int j = 0; j <= i; j++) {
          // Compute the second derivatives of d(b13[0:3])/(d(theta))
          n13[0] = 0.5 * rt[0] * Na[i] * N[j];
          n13[1] = 0.5 * rt[3] * Na[i] * N[j];
          n13[2] = 0.5 * rt[6] * Na[i] * N[j];
          n13[3] = 0.5 * rt[1] * Na[i] * N[j];
          n13[4] = 0.5 * rt[4] * Na[i] * N[j];
          n13[5] = 0.5 * rt[7] * Na[i] * N[j];
          n13[6] = 0.5 * rt[2] * Na[i] * N[j];
          n13[7] = 0.5 * rt[5] * Na[i] * N[j];
          n13[8] = 0.5 * rt[8] * Na[i] * N[j];
          n13 += 9;

          // Compute the second derivatives of d(b13[3:6])/du
          n13[0] = 0.5 * rt[0] * N[i] * Na[j];
          n13[1] = 0.5 * rt[1] * N[i] * Na[j];
          n13[2] = 0.5 * rt[2] * N[i] * Na[j];
          n13[3] = 0.5 * rt[3] * N[i] * Na[j];
          n13[4] = 0.5 * rt[4] * N[i] * Na[j];
          n13[5] = 0.5 * rt[5] * N[i] * Na[j];
          n13[6] = 0.5 * rt[6] * N[i] * Na[j];
          n13[7] = 0.5 * rt[7] * N[i] * Na[j];
          n13[8] = 0.5 * rt[8] * N[i] * Na[j];
          n13 += 9;

          // Compute the second derivatives of d(b[3:6])/d(theta) - symmetric
          n13[0] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[0] + (Xd[1] + Ud[2]) * rtt[1] +
                    (Xd[2] + Ud[4]) * rtt[2]);
          n13[1] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[3] + (Xd[1] + Ud[2]) * rtt[4] +
                    (Xd[2] + Ud[4]) * rtt[5]);
          n13[2] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[6] + (Xd[1] + Ud[2]) * rtt[7] +
                    (Xd[2] + Ud[4]) * rtt[8]);
          n13[3] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[9] + (Xd[1] + Ud[2]) * rtt[10] +
                    (Xd[2] + Ud[4]) * rtt[11]);
          n13[4] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[12] + (Xd[1] + Ud[2]) * rtt[13] +
                    (Xd[2] + Ud[4]) * rtt[14]);
          n13[5] = 0.5 * N[j] * N[i] *
                   ((Xd[0] + Ud[0]) * rtt[15] + (Xd[1] + Ud[2]) * rtt[16] +
                    (Xd[2] + Ud[4]) * rtt[17]);
          n13 += 6;
        }
      }
    }
  }
}

/*
  Add the second derivatives of the strain times the stress to the
  stiffness matrix.

  input:
  matrix: the stiffness matrix being assembled
  scale: scale factor to muliply the result by
  stress: the stress at the current location

  n13, n23: the second derivatives of the strain at all the tying
  points -> these are large arrays
  N11, N22, N12: the strain interpolation functions interpolations
*/
template <int order, int tying_order>
void add_lr_tying_stress_nmat(TacsScalar matrix[], TacsScalar scale,
                              const TacsScalar stress[], const TacsScalar n13[],
                              const TacsScalar n23[], const TacsScalar tx[],
                              const double N11[], const double N22[],
                              const double N12[], const double knots[],
                              const double pknots[]) {
  // Number of variables in the stiffness matrix
  const int nvars = 6 * order * order;

  // Compute the second derivatives of the in-plane strain
  TacsScalar h11[order * order * (order * order + 1) / 2];
  TacsScalar h22[order * order * (order * order + 1) / 2];
  TacsScalar h12[order * order * (order * order + 1) / 2];
  memset(h11, 0, sizeof(h11));
  memset(h22, 0, sizeof(h22));
  memset(h12, 0, sizeof(h12));
  shellutils::nonlinear_tying_nmat<order, tying_order>(h11, h22, h12, N11, N22,
                                                       N12, knots, pknots);

  // The second derivatives of the shear strain components
  // There are 2*9 + 6 = 24 second derivative components per point
  TacsScalar h13[12 * (order * order * (order * order + 1))];
  TacsScalar h23[12 * (order * order * (order * order + 1))];
  memset(h13, 0, sizeof(h13));
  memset(h23, 0, sizeof(h23));

  // Loop over each tying-point in the element
  for (int i = 0; i < tying_order * (tying_order - 1); i++) {
    TacsScalar* _h13 = h13;
    TacsScalar* _h23 = h23;

    // Compute the interpolation at this point
    for (int j = 0; j < (order * order * (order * order + 1)) / 2; j++) {
      _h13[0] += n13[0] * N11[0];
      _h13[1] += n13[1] * N11[0];
      _h13[2] += n13[2] * N11[0];
      _h13[3] += n13[3] * N11[0];
      _h13[4] += n13[4] * N11[0];
      _h13[5] += n13[5] * N11[0];
      _h13[6] += n13[6] * N11[0];
      _h13[7] += n13[7] * N11[0];
      _h13[8] += n13[8] * N11[0];
      _h13[9] += n13[9] * N11[0];
      _h13[10] += n13[10] * N11[0];
      _h13[11] += n13[11] * N11[0];
      _h13[12] += n13[12] * N11[0];
      _h13[13] += n13[13] * N11[0];
      _h13[14] += n13[14] * N11[0];
      _h13[15] += n13[15] * N11[0];
      _h13[16] += n13[16] * N11[0];
      _h13[17] += n13[17] * N11[0];
      _h13[18] += n13[18] * N11[0];
      _h13[19] += n13[19] * N11[0];
      _h13[20] += n13[20] * N11[0];
      _h13[21] += n13[21] * N11[0];
      _h13[22] += n13[22] * N11[0];
      _h13[23] += n13[23] * N11[0];
      n13 += 24;
      _h13 += 24;

      _h23[0] += n23[0] * N22[0];
      _h23[1] += n23[1] * N22[0];
      _h23[2] += n23[2] * N22[0];
      _h23[3] += n23[3] * N22[0];
      _h23[4] += n23[4] * N22[0];
      _h23[5] += n23[5] * N22[0];
      _h23[6] += n23[6] * N22[0];
      _h23[7] += n23[7] * N22[0];
      _h23[8] += n23[8] * N22[0];
      _h23[9] += n23[9] * N22[0];
      _h23[10] += n23[10] * N22[0];
      _h23[11] += n23[11] * N22[0];
      _h23[12] += n23[12] * N22[0];
      _h23[13] += n23[13] * N22[0];
      _h23[14] += n23[14] * N22[0];
      _h23[15] += n23[15] * N22[0];
      _h23[16] += n23[16] * N22[0];
      _h23[17] += n23[17] * N22[0];
      _h23[18] += n23[18] * N22[0];
      _h23[19] += n23[19] * N22[0];
      _h23[20] += n23[20] * N22[0];
      _h23[21] += n23[21] * N22[0];
      _h23[22] += n23[22] * N22[0];
      _h23[23] += n23[23] * N22[0];
      n23 += 24;
      _h23 += 24;
    }
    N11++;
    N22++;
  }

  // Now add the second derivative contributions to the
  // stiffness matrix - scaled by 'scale'
  TacsScalar* _h11 = h11;
  TacsScalar* _h22 = h22;
  TacsScalar* _h12 = h12;
  TacsScalar* _h13 = h13;
  TacsScalar* _h23 = h23;

  for (int i = 0; i < order * order; i++) {
    const int row = 6 * i;
    for (int j = 0; j <= i; j++) {
      const int col = 6 * j;

      // Evaluate the second derivative of the tensorial strain.
      // Note that this is NOT the engineering strain!
      TacsScalar g[6], s[6];
      g[0] = _h11[0];
      g[1] = _h22[0];
      g[5] = _h12[0];
      g[2] = g[3] = g[4] = 0.0;

      // Transform to the local axis - note that this employs the
      // transformation for stress - since we're using tensorial strain -
      // not engineering strain
      Tensor::transform3DStress(s, g, tx);

      // Evaluate the contribution - note the factor of 2.0 due to the
      // use of tensorial strain
      TacsScalar val =
          scale *
          (stress[0] * s[0] + stress[1] * s[1] +
           2.0 * (stress[2] * s[5] + stress[6] * s[3] + stress[7] * s[4]));

      // Add values to the diagonal of each component
      for (int k = 0; k < 3; k++) {
        matrix[(row + k) * nvars + col + k] += val;
      }

      // Increment the pointers
      _h11++;
      _h22++;
      _h12++;

      // Compute the mixed-partial terms from the shear strain
      for (int ii = 0; ii < 3; ii++) {
        for (int jj = 3; jj < 6; jj++) {
          // Determine the second derivative of the tensorial strain
          TacsScalar g[6], s[6];
          g[0] = g[1] = g[2] = g[5] = 0.0;
          g[3] = _h23[0];
          g[4] = _h13[0];

          // Transform the tensorial strain to the local axis
          // Note that the stress transformation is used since tensorial
          // strain is used (not engineering strain)
          Tensor::transform3DStress(s, g, tx);

          // Detmine the contribution to the stiffness matrix - the
          // factor of 2.0 is due to the use of tensorial strain
          TacsScalar val =
              scale *
              (stress[0] * s[0] + stress[1] * s[1] +
               2.0 * (stress[2] * s[5] + stress[6] * s[3] + stress[7] * s[4]));

          matrix[(row + ii) * nvars + (col + jj)] += val;

          // Increment the pointers to the second derivatives of the
          // strain
          _h23++;
          _h13++;
        }
      }

      // Compute the mixed-partial terms from the shear strain
      for (int ii = 3; ii < 6; ii++) {
        for (int jj = 0; jj < 3; jj++) {
          // Determine the second derivative of the tensorial strain
          TacsScalar g[6], s[6];
          g[0] = g[1] = g[2] = g[5] = 0.0;
          g[3] = _h23[0];
          g[4] = _h13[0];

          // Transform the tensorial strain to the local axis
          // Note that the stress transformation is used since tensorial
          // strain is used (not engineering strain)
          Tensor::transform3DStress(s, g, tx);

          // Detmine the contribution to the stiffness matrix - the
          // factor of 2.0 is due to the use of tensorial strain
          TacsScalar val =
              scale *
              (stress[0] * s[0] + stress[1] * s[1] +
               2.0 * (stress[2] * s[5] + stress[6] * s[3] + stress[7] * s[4]));

          matrix[(row + ii) * nvars + (col + jj)] += val;

          // Increment the pointers to the second derivatives of the
          // strain
          _h23++;
          _h13++;
        }
      }

      // Now, compute the contributions from the shear strain
      for (int k = 0; k < 6; k++) {
        // Determine the second derivative of the tensorial strain
        TacsScalar g[6], s[6];
        g[0] = g[1] = g[2] = g[5] = 0.0;
        g[3] = _h23[0];
        g[4] = _h13[0];

        // Transform the tensorial strain to the local axis
        // Note that the stress transformation is used since tensorial
        // strain is used (not engineering strain)
        Tensor::transform3DStress(s, g, tx);

        // Detmine the contribution to the stiffness matrix - the
        // factor of 2.0 is due to the use of tensorial strain
        TacsScalar val =
            scale *
            (stress[0] * s[0] + stress[1] * s[1] +
             2.0 * (stress[2] * s[5] + stress[6] * s[3] + stress[7] * s[4]));

        if (k < 3) {
          // Add the derivative C,kk along the sub-diagonals
          matrix[(row + 3 + k) * nvars + col + 3 + k] += val;
        } else if (k == 3) {
          // The derivative for C,12
          matrix[(row + 3) * nvars + col + 4] += val;
          matrix[(row + 4) * nvars + col + 3] += val;
        } else if (k == 4) {
          // The derivative for C,13
          matrix[(row + 3) * nvars + col + 5] += val;
          matrix[(row + 5) * nvars + col + 3] += val;
        } else if (k == 5) {
          // The derivative for C,23
          matrix[(row + 4) * nvars + col + 5] += val;
          matrix[(row + 5) * nvars + col + 4] += val;
        }

        // Increment the pointers to the second derivatives of the
        // strain
        _h23++;
        _h13++;
      }
    }
  }
}

TACS_END_NAMESPACE

#endif
