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

#include "ShellUtils.h"

#include "FElibrary.h"

/*
  Common structural shell utilities for small strains
*/

TACS_BEGIN_NAMESPACE(shellutils)

/*
  Compute the derivative of the inverse of the Jacobian
  transformation with respect to the through thickness direction
  zeta. This will be used to determine the variation of the strain
  through the thickness

  zJinv = [
  x,xi   y,xi   z,xi
  x,eta  y,eta  z,eta
  x,zeta y,zeta z,zeta ],zeta =
  - Jinv * dX_{xi}/d zeta * Jinv

  dX_{xi}/d zeta =
  [ n_xi  ]
  [ n_eta ]
  [ 0     ]

  [ 0 1 2 ][ n_xi  ][ 0 1 2 ]
  [ 3 4 5 ][ n_eta ][ 3 4 5 ]
  [ 6 7 8 ][ 0     ][ 6 7 8 ]

  input:
  Jinv: the inverse of Xd
  n_xi, n_eta: the derivatives of the shell normal along the
  natural coordinates xi and eta
*/
static inline void compute_jacobian_zeta(TacsScalar zJinv[],
                                         const TacsScalar Jinv[],
                                         const TacsScalar n_xi[],
                                         const TacsScalar n_eta[]) {
  TacsScalar tp[9];
  tp[0] = -(Jinv[0] * n_xi[0] + Jinv[1] * n_eta[0]);
  tp[3] = -(Jinv[3] * n_xi[0] + Jinv[4] * n_eta[0]);
  tp[6] = -(Jinv[6] * n_xi[0] + Jinv[7] * n_eta[0]);

  tp[1] = -(Jinv[0] * n_xi[1] + Jinv[1] * n_eta[1]);
  tp[4] = -(Jinv[3] * n_xi[1] + Jinv[4] * n_eta[1]);
  tp[7] = -(Jinv[6] * n_xi[1] + Jinv[7] * n_eta[1]);

  tp[2] = -(Jinv[0] * n_xi[2] + Jinv[1] * n_eta[2]);
  tp[5] = -(Jinv[3] * n_xi[2] + Jinv[4] * n_eta[2]);
  tp[8] = -(Jinv[6] * n_xi[2] + Jinv[7] * n_eta[2]);

  zJinv[0] = tp[0] * Jinv[0] + tp[1] * Jinv[3] + tp[2] * Jinv[6];
  zJinv[3] = tp[3] * Jinv[0] + tp[4] * Jinv[3] + tp[5] * Jinv[6];
  zJinv[6] = tp[6] * Jinv[0] + tp[7] * Jinv[3] + tp[8] * Jinv[6];

  zJinv[1] = tp[0] * Jinv[1] + tp[1] * Jinv[4] + tp[2] * Jinv[7];
  zJinv[4] = tp[3] * Jinv[1] + tp[4] * Jinv[4] + tp[5] * Jinv[7];
  zJinv[7] = tp[6] * Jinv[1] + tp[7] * Jinv[4] + tp[8] * Jinv[7];

  zJinv[2] = tp[0] * Jinv[2] + tp[1] * Jinv[5] + tp[2] * Jinv[8];
  zJinv[5] = tp[3] * Jinv[2] + tp[4] * Jinv[5] + tp[5] * Jinv[8];
  zJinv[8] = tp[6] * Jinv[2] + tp[7] * Jinv[5] + tp[8] * Jinv[8];
}

/*
   Compute the sensitivity of zJinv using the sensitivity information

   input:
   Jinv, dJinv: the inverse and derivative of the inverse of Xd
   n_xi, dn_xi: the derivative of the normal along xi and its sensitivity
   n_eta, dn_eta: the derivative of the normal along eta and its sensitivity
*/
static inline void compute_jacobian_zeta_sens(
    TacsScalar dzJinv[], const TacsScalar Jinv[], const TacsScalar dJinv[],
    const TacsScalar n_xi[], const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[]) {
  TacsScalar tp[9], dtp[9];
  tp[0] = -(Jinv[0] * n_xi[0] + Jinv[1] * n_eta[0]);
  tp[3] = -(Jinv[3] * n_xi[0] + Jinv[4] * n_eta[0]);
  tp[6] = -(Jinv[6] * n_xi[0] + Jinv[7] * n_eta[0]);

  tp[1] = -(Jinv[0] * n_xi[1] + Jinv[1] * n_eta[1]);
  tp[4] = -(Jinv[3] * n_xi[1] + Jinv[4] * n_eta[1]);
  tp[7] = -(Jinv[6] * n_xi[1] + Jinv[7] * n_eta[1]);

  tp[2] = -(Jinv[0] * n_xi[2] + Jinv[1] * n_eta[2]);
  tp[5] = -(Jinv[3] * n_xi[2] + Jinv[4] * n_eta[2]);
  tp[8] = -(Jinv[6] * n_xi[2] + Jinv[7] * n_eta[2]);

  dtp[0] = -(dJinv[0] * n_xi[0] + dJinv[1] * n_eta[0] + Jinv[0] * dn_xi[0] +
             Jinv[1] * dn_eta[0]);
  dtp[3] = -(dJinv[3] * n_xi[0] + dJinv[4] * n_eta[0] + Jinv[3] * dn_xi[0] +
             Jinv[4] * dn_eta[0]);
  dtp[6] = -(dJinv[6] * n_xi[0] + dJinv[7] * n_eta[0] + Jinv[6] * dn_xi[0] +
             Jinv[7] * dn_eta[0]);

  dtp[1] = -(dJinv[0] * n_xi[1] + dJinv[1] * n_eta[1] + Jinv[0] * dn_xi[1] +
             Jinv[1] * dn_eta[1]);
  dtp[4] = -(dJinv[3] * n_xi[1] + dJinv[4] * n_eta[1] + Jinv[3] * dn_xi[1] +
             Jinv[4] * dn_eta[1]);
  dtp[7] = -(dJinv[6] * n_xi[1] + dJinv[7] * n_eta[1] + Jinv[6] * dn_xi[1] +
             Jinv[7] * dn_eta[1]);

  dtp[2] = -(dJinv[0] * n_xi[2] + dJinv[1] * n_eta[2] + Jinv[0] * dn_xi[2] +
             Jinv[1] * dn_eta[2]);
  dtp[5] = -(dJinv[3] * n_xi[2] + dJinv[4] * n_eta[2] + Jinv[3] * dn_xi[2] +
             Jinv[4] * dn_eta[2]);
  dtp[8] = -(dJinv[6] * n_xi[2] + dJinv[7] * n_eta[2] + Jinv[6] * dn_xi[2] +
             Jinv[7] * dn_eta[2]);

  dzJinv[0] = (dtp[0] * Jinv[0] + dtp[1] * Jinv[3] + dtp[2] * Jinv[6] +
               tp[0] * dJinv[0] + tp[1] * dJinv[3] + tp[2] * dJinv[6]);
  dzJinv[3] = (dtp[3] * Jinv[0] + dtp[4] * Jinv[3] + dtp[5] * Jinv[6] +
               tp[3] * dJinv[0] + tp[4] * dJinv[3] + tp[5] * dJinv[6]);
  dzJinv[6] = (dtp[6] * Jinv[0] + dtp[7] * Jinv[3] + dtp[8] * Jinv[6] +
               tp[6] * dJinv[0] + tp[7] * dJinv[3] + tp[8] * dJinv[6]);

  dzJinv[1] = (dtp[0] * Jinv[1] + dtp[1] * Jinv[4] + dtp[2] * Jinv[7] +
               tp[0] * dJinv[1] + tp[1] * dJinv[4] + tp[2] * dJinv[7]);
  dzJinv[4] = (dtp[3] * Jinv[1] + dtp[4] * Jinv[4] + dtp[5] * Jinv[7] +
               tp[3] * dJinv[1] + tp[4] * dJinv[4] + tp[5] * dJinv[7]);
  dzJinv[7] = (dtp[6] * Jinv[1] + dtp[7] * Jinv[4] + dtp[8] * Jinv[7] +
               tp[6] * dJinv[1] + tp[7] * dJinv[4] + tp[8] * dJinv[7]);

  dzJinv[2] = (dtp[0] * Jinv[2] + dtp[1] * Jinv[5] + dtp[2] * Jinv[8] +
               tp[0] * dJinv[2] + tp[1] * dJinv[5] + tp[2] * dJinv[8]);
  dzJinv[5] = (dtp[3] * Jinv[2] + dtp[4] * Jinv[5] + dtp[5] * Jinv[8] +
               tp[3] * dJinv[2] + tp[4] * dJinv[5] + tp[5] * dJinv[8]);
  dzJinv[8] = (dtp[6] * Jinv[2] + dtp[7] * Jinv[5] + dtp[8] * Jinv[8] +
               tp[6] * dJinv[2] + tp[7] * dJinv[5] + tp[8] * dJinv[8]);
}

/*
  Compute the product tx = t*Jinv and ztx = t*zJinv

  These are required for computing the transformation between the
  global axis and the shell coordinate system. These quantities are
  eventually used to compute the strain in the local shell coordinate
  system.

  input:
  t:     the transformation
  Jinv:  the inverse of Xd
  zJinv: the derivative of the inverse of Xd in the through-thickness
  direction

  output:
  tx: the transformation between the global and shell-fixed frames
  ztx: the derivative of the transformation between global and
  shell-fixed frames in the through-thickness direction (zeta)
*/
static inline void transform_jacobian_product(TacsScalar tx[], TacsScalar ztx[],
                                              const TacsScalar t[],
                                              const TacsScalar Jinv[],
                                              const TacsScalar zJinv[]) {
  // Compute the transformation between the global and
  // shell-fixed frames
  tx[0] = t[0] * Jinv[0] + t[1] * Jinv[3] + t[2] * Jinv[6];
  tx[3] = t[3] * Jinv[0] + t[4] * Jinv[3] + t[5] * Jinv[6];
  tx[6] = t[6] * Jinv[0] + t[7] * Jinv[3] + t[8] * Jinv[6];

  tx[1] = t[0] * Jinv[1] + t[1] * Jinv[4] + t[2] * Jinv[7];
  tx[4] = t[3] * Jinv[1] + t[4] * Jinv[4] + t[5] * Jinv[7];
  tx[7] = t[6] * Jinv[1] + t[7] * Jinv[4] + t[8] * Jinv[7];

  tx[2] = t[0] * Jinv[2] + t[1] * Jinv[5] + t[2] * Jinv[8];
  tx[5] = t[3] * Jinv[2] + t[4] * Jinv[5] + t[5] * Jinv[8];
  tx[8] = t[6] * Jinv[2] + t[7] * Jinv[5] + t[8] * Jinv[8];

  // Compute the derivative of the transformation between
  // the global and shell-fixed frames
  ztx[0] = t[0] * zJinv[0] + t[1] * zJinv[3] + t[2] * zJinv[6];
  ztx[3] = t[3] * zJinv[0] + t[4] * zJinv[3] + t[5] * zJinv[6];
  ztx[6] = t[6] * zJinv[0] + t[7] * zJinv[3] + t[8] * zJinv[6];

  ztx[1] = t[0] * zJinv[1] + t[1] * zJinv[4] + t[2] * zJinv[7];
  ztx[4] = t[3] * zJinv[1] + t[4] * zJinv[4] + t[5] * zJinv[7];
  ztx[7] = t[6] * zJinv[1] + t[7] * zJinv[4] + t[8] * zJinv[7];

  ztx[2] = t[0] * zJinv[2] + t[1] * zJinv[5] + t[2] * zJinv[8];
  ztx[5] = t[3] * zJinv[2] + t[4] * zJinv[5] + t[5] * zJinv[8];
  ztx[8] = t[6] * zJinv[2] + t[7] * zJinv[5] + t[8] * zJinv[8];
}

/*
  Compute the product tx = t*Jinv and ztx = t*zJinv and the
  derivative of the products w.r.t. the input perturbations
*/
static inline void transform_jacobian_product_sens(
    TacsScalar dtx[], TacsScalar dztx[], const TacsScalar t[],
    const TacsScalar dt[], const TacsScalar Jinv[], const TacsScalar zJinv[],
    const TacsScalar dJinv[], const TacsScalar dzJinv[]) {
  // Compute the derivative
  dtx[0] = dt[0] * Jinv[0] + dt[1] * Jinv[3] + dt[2] * Jinv[6];
  dtx[3] = dt[3] * Jinv[0] + dt[4] * Jinv[3] + dt[5] * Jinv[6];
  dtx[6] = dt[6] * Jinv[0] + dt[7] * Jinv[3] + dt[8] * Jinv[6];

  dtx[1] = dt[0] * Jinv[1] + dt[1] * Jinv[4] + dt[2] * Jinv[7];
  dtx[4] = dt[3] * Jinv[1] + dt[4] * Jinv[4] + dt[5] * Jinv[7];
  dtx[7] = dt[6] * Jinv[1] + dt[7] * Jinv[4] + dt[8] * Jinv[7];

  dtx[2] = dt[0] * Jinv[2] + dt[1] * Jinv[5] + dt[2] * Jinv[8];
  dtx[5] = dt[3] * Jinv[2] + dt[4] * Jinv[5] + dt[5] * Jinv[8];
  dtx[8] = dt[6] * Jinv[2] + dt[7] * Jinv[5] + dt[8] * Jinv[8];

  // Complete the derivative
  dtx[0] += t[0] * dJinv[0] + t[1] * dJinv[3] + t[2] * dJinv[6];
  dtx[3] += t[3] * dJinv[0] + t[4] * dJinv[3] + t[5] * dJinv[6];
  dtx[6] += t[6] * dJinv[0] + t[7] * dJinv[3] + t[8] * dJinv[6];

  dtx[1] += t[0] * dJinv[1] + t[1] * dJinv[4] + t[2] * dJinv[7];
  dtx[4] += t[3] * dJinv[1] + t[4] * dJinv[4] + t[5] * dJinv[7];
  dtx[7] += t[6] * dJinv[1] + t[7] * dJinv[4] + t[8] * dJinv[7];

  dtx[2] += t[0] * dJinv[2] + t[1] * dJinv[5] + t[2] * dJinv[8];
  dtx[5] += t[3] * dJinv[2] + t[4] * dJinv[5] + t[5] * dJinv[8];
  dtx[8] += t[6] * dJinv[2] + t[7] * dJinv[5] + t[8] * dJinv[8];

  // Compute the derivative
  dztx[0] = dt[0] * zJinv[0] + dt[1] * zJinv[3] + dt[2] * zJinv[6];
  dztx[3] = dt[3] * zJinv[0] + dt[4] * zJinv[3] + dt[5] * zJinv[6];
  dztx[6] = dt[6] * zJinv[0] + dt[7] * zJinv[3] + dt[8] * zJinv[6];

  dztx[1] = dt[0] * zJinv[1] + dt[1] * zJinv[4] + dt[2] * zJinv[7];
  dztx[4] = dt[3] * zJinv[1] + dt[4] * zJinv[4] + dt[5] * zJinv[7];
  dztx[7] = dt[6] * zJinv[1] + dt[7] * zJinv[4] + dt[8] * zJinv[7];

  dztx[2] = dt[0] * zJinv[2] + dt[1] * zJinv[5] + dt[2] * zJinv[8];
  dztx[5] = dt[3] * zJinv[2] + dt[4] * zJinv[5] + dt[5] * zJinv[8];
  dztx[8] = dt[6] * zJinv[2] + dt[7] * zJinv[5] + dt[8] * zJinv[8];

  // Complete the derivative
  dztx[0] += t[0] * dzJinv[0] + t[1] * dzJinv[3] + t[2] * dzJinv[6];
  dztx[3] += t[3] * dzJinv[0] + t[4] * dzJinv[3] + t[5] * dzJinv[6];
  dztx[6] += t[6] * dzJinv[0] + t[7] * dzJinv[3] + t[8] * dzJinv[6];

  dztx[1] += t[0] * dzJinv[1] + t[1] * dzJinv[4] + t[2] * dzJinv[7];
  dztx[4] += t[3] * dzJinv[1] + t[4] * dzJinv[4] + t[5] * dzJinv[7];
  dztx[7] += t[6] * dzJinv[1] + t[7] * dzJinv[4] + t[8] * dzJinv[7];

  dztx[2] += t[0] * dzJinv[2] + t[1] * dzJinv[5] + t[2] * dzJinv[8];
  dztx[5] += t[3] * dzJinv[2] + t[4] * dzJinv[5] + t[5] * dzJinv[8];
  dztx[8] += t[6] * dzJinv[2] + t[7] * dzJinv[5] + t[8] * dzJinv[8];
}

/*
  Compute the cross product of two vectors
*/
static inline void cross_product(TacsScalar n[], const TacsScalar A[],
                                 const TacsScalar B[]) {
  n[0] = A[1] * B[2] - A[2] * B[1];
  n[1] = A[2] * B[0] - A[0] * B[2];
  n[2] = A[0] * B[1] - A[1] * B[0];
}

/*
  Compute the derivative of the cross product
*/
static inline void cross_product_sens(TacsScalar dn[], const TacsScalar A[],
                                      const TacsScalar B[],
                                      const TacsScalar dA[],
                                      const TacsScalar dB[]) {
  dn[0] = (dA[1] * B[2] - dA[2] * B[1] + A[1] * dB[2] - A[2] * dB[1]);
  dn[1] = (dA[2] * B[0] - dA[0] * B[2] + A[2] * dB[0] - A[0] * dB[2]);
  dn[2] = (dA[0] * B[1] - dA[1] * B[0] + A[0] * dB[1] - A[1] * dB[0]);
}

/*
  Normalize a vector in place
*/
static inline TacsScalar normalize_vector(TacsScalar n[]) {
  TacsScalar norm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);

  if (norm != 0.0) {
    TacsScalar norm_inv = 1.0 / norm;
    n[0] *= norm_inv;
    n[1] *= norm_inv;
    n[2] *= norm_inv;
  }

  return norm;
}

/*
  Determine the derivative of a normalize operation.

  input/output:
  dn: the sensitivity of n - the un-normalized input

  input:
  n: the vector that was previously normalized
  norm: the norm of the vector n before normalization

  Note that here the input vector n must already be normalized
*/
static inline void normalize_vector_sens(TacsScalar dn[], TacsScalar norm,
                                         const TacsScalar n[]) {
  if (norm != 0.0) {
    TacsScalar norm_inv = 1.0 / norm;
    TacsScalar dnorm = n[0] * dn[0] + n[1] * dn[1] + n[2] * dn[2];

    dn[0] = (dn[0] - n[0] * dnorm) * norm_inv;
    dn[1] = (dn[1] - n[1] * dnorm) * norm_inv;
    dn[2] = (dn[2] - n[2] * dnorm) * norm_inv;
  } else {
    dn[0] = dn[1] = dn[2] = 0.0;
  }
}

/*
  Normalize a vector and copy over the result to a new vector

  input:
  r: the vector to be normalized - not modified!

  output:
  n: the noralized vector

  returns: the norm of r
*/
static inline TacsScalar normalize_vector(TacsScalar n[],
                                          const TacsScalar r[]) {
  TacsScalar norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

  if (norm != 0.0) {
    TacsScalar norm_inv = 1.0 / norm;
    n[0] = norm_inv * r[0];
    n[1] = norm_inv * r[1];
    n[2] = norm_inv * r[2];
  } else {
    n[0] = n[1] = n[2] = 0.0;
  }

  return norm;
}

/*
  Compute the sensitivity of the normalized vector

  input:
  norm: the norm of r
  r: the vector to be normalized
  dr: the sensitivity of the vector to be normalized

  output:
  dn: the sensitivity of the normalized vector
*/
static inline void normalize_vector_sens(TacsScalar dn[], TacsScalar norm,
                                         const TacsScalar r[],
                                         const TacsScalar dr[]) {
  if (norm != 0.0) {
    TacsScalar norm_inv = 1.0 / norm;
    TacsScalar dnorm_inv = -norm_inv * norm_inv * norm_inv *
                           (r[0] * dr[0] + r[1] * dr[1] + r[2] * dr[2]);

    dn[0] = norm_inv * dr[0] + dnorm_inv * r[0];
    dn[1] = norm_inv * dr[1] + dnorm_inv * r[1];
    dn[2] = norm_inv * dr[2] + dnorm_inv * r[2];
  } else {
    dn[0] = dn[1] = dn[2] = 0.0;
  }
}

/*
  Compute the sensitivity of the normalieze_vector_sens - this is a
  bit confusing

  input:
  norm: the norm of r
  r: the vector to be normalized
  dr: the sensitivity of the vector to be normalized
  rs: the perturbation of r
  drs: the perturbation of dr

  output:
  dns: the sensitivity of the pertubation to r - through the inputs rs/drs
*/
static inline void normalize_vector_second_sens(
    TacsScalar dns[], TacsScalar norm, const TacsScalar r[],
    const TacsScalar dr[], const TacsScalar rs[], const TacsScalar drs[]) {
  if (norm != 0.0) {
    TacsScalar norm_inv = 1.0 / norm;

    TacsScalar dnorm_inv = -norm_inv * norm_inv * norm_inv *
                           (r[0] * dr[0] + r[1] * dr[1] + r[2] * dr[2]);

    TacsScalar norm_invs = -norm_inv * norm_inv * norm_inv *
                           (r[0] * rs[0] + r[1] * rs[1] + r[2] * rs[2]);

    TacsScalar dnorm_invs =
        -norm_inv * norm_inv * norm_inv *
            (r[0] * drs[0] + r[1] * drs[1] + r[2] * drs[2] + rs[0] * dr[0] +
             rs[1] * dr[1] + rs[2] * dr[2]) -
        (3.0 * dnorm_inv) * norm_inv * norm_inv *
            (r[0] * rs[0] + r[1] * rs[1] + r[2] * rs[2]);

    dns[0] = (norm_invs * dr[0] + norm_inv * drs[0] + dnorm_invs * r[0] +
              dnorm_inv * rs[0]);

    dns[1] = (norm_invs * dr[1] + norm_inv * drs[1] + dnorm_invs * r[1] +
              dnorm_inv * rs[1]);

    dns[2] = (norm_invs * dr[2] + norm_inv * drs[2] + dnorm_invs * r[2] +
              dnorm_inv * rs[2]);
  } else {
    dns[0] = dns[1] = dns[2] = 0.0;
  }
}

/*
  Compute the transformation that is required to take the derivatives
  w.r.t. the local shell coordinates (xi,eta,zeta) to the local shell
  frame (e1,e2,e3).

  This transformation involves two steps:

  U_{x,x} = U_{,xi} * X_{,xi}^{-1}
  U_{x,e} = U_{,xi} * X_{,xi}^{-1} * T^{T}
  = U_{,xi} * tx^{T}

  tx^{T} = X_{,xi}^{-1} * T^{T}
  tx = T * X_{,xi}^{-T}
  = T * Jinv

  input:
  Xd  == The derivatives of the shell coordinates
  Xdd == The second derivatives of the shell coordinates

  output:
  n     == The shell normal direction
  n_xi  == The derivative of n w.r.t. xi
  n_eta == The derivative of n w.r.t. eta
  tx    == The transformation
  ztx   == The derivative of the transformation w.r.t. zeta

  intermediaries:
  Jinv == X_{,xi}^{-T}, the inverse of the Jacobian
  t    == The transformation from X to e
*/
TacsScalar compute_transform(TacsScalar t[], TacsScalar tx[], TacsScalar ztx[],
                             TacsScalar n[], TacsScalar n_xi[],
                             TacsScalar n_eta[], TacsScalar Xd[],
                             const TacsScalar Xdd[]) {
  // n = Xd[0] x Xd[3]
  // |   i     j      k    |
  // | Xd[0]  Xd[1]  Xd[2] |
  // | Xd[3]  Xd[4]  Xd[5] |
  cross_product(n, &Xd[0], &Xd[3]);
  cross_product_sens(n_xi, &Xd[0], &Xd[3], &Xdd[0], &Xdd[3]);
  cross_product_sens(n_eta, &Xd[0], &Xd[3], &Xdd[3], &Xdd[6]);

  // Normalize the normal vector
  TacsScalar norm = normalize_vector(n);

  // Now take the derivative of the normal vector along
  // the directions xi and eta
  normalize_vector_sens(n_xi, norm, n);
  normalize_vector_sens(n_eta, norm, n);

  Xd[6] = n[0];
  Xd[7] = n[1];
  Xd[8] = n[2];

  // Jinv = Xd^{-1}
  // zJinv = p(Xd^{-1})/p zeta
  TacsScalar Jinv[9], zJinv[9];
  TacsScalar h = FElibrary::jacobian3d(Xd, Jinv);
  compute_jacobian_zeta(zJinv, Jinv, n_xi, n_eta);

  // t[0:3] = Xd[0:3]//||Xd[0:3]||
  t[0] = Xd[0];
  t[1] = Xd[1];
  t[2] = Xd[2];
  normalize_vector(&t[0]);

  // t[3:6] = n x t[0:3]
  cross_product(&t[3], n, t);

  // t[6:9] = n
  t[6] = n[0];
  t[7] = n[1];
  t[8] = n[2];

  transform_jacobian_product(tx, ztx, t, Jinv, zJinv);

  return h;
}

/*
  Compute the transformation and the sensitivity of the transformation
  that is required to take the derivatives w.r.t. the local shell
  coordinates (xi,eta,zeta) to the local shell frame (e1,e2,e3).

  input:
  Xd, Xdd: the first and second derivatives of the shell surface
  Na, Nb: the first derivatives of the shape functions
  Naa, Nab, Nbb: the second derivatives of the shape functions
  num_nodes: the number of nodes

  output:
  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivities of the transformations

  n: the surface normal
  n_xi, n_eta: the derivatives of the surface normal

  dn: the sensitivity of the surface normal
  dn_xi, dn_eta: the sensitivity of the derivatives of the surface normal
*/
TacsScalar compute_transform_sens(
    TacsScalar dh[], TacsScalar t[], TacsScalar dt[], TacsScalar tx[],
    TacsScalar dtx[], TacsScalar ztx[], TacsScalar dztx[], TacsScalar n[],
    TacsScalar dn[], TacsScalar n_xi[], TacsScalar dn_xi[], TacsScalar n_eta[],
    TacsScalar dn_eta[], TacsScalar Xd[], const TacsScalar Xdd[],
    const double Na[], const double Nb[], const double Naa[],
    const double Nab[], const double Nbb[], int num_nodes) {
  // n = Xd[0] x Xd[3]
  // |   i     j      k    |
  // | Xd[0]  Xd[1]  Xd[2] |
  // | Xd[3]  Xd[4]  Xd[5] |
  TacsScalar rn[3], rn_xi[3], rn_eta[3];
  cross_product(rn, &Xd[0], &Xd[3]);
  cross_product_sens(rn_xi, &Xd[0], &Xd[3], &Xdd[0], &Xdd[3]);
  cross_product_sens(rn_eta, &Xd[0], &Xd[3], &Xdd[3], &Xdd[6]);

  // Normalize the normal vector
  TacsScalar norm = normalize_vector(n, rn);

  // Now take the derivative of the normal vector along
  // the directions xi and eta
  normalize_vector_sens(n_xi, norm, rn, rn_xi);
  normalize_vector_sens(n_eta, norm, rn, rn_eta);

  Xd[6] = n[0];
  Xd[7] = n[1];
  Xd[8] = n[2];

  // Jinv = Xd^{-1}
  // zJinv = p(Xd^{-1})/p zeta
  TacsScalar Jinv[9], zJinv[9];
  TacsScalar h = FElibrary::jacobian3d(Xd, Jinv);
  compute_jacobian_zeta(zJinv, Jinv, n_xi, n_eta);

  // t[0:3] = Xd[0:3]//||Xd[0:3]||
  t[0] = Xd[0];
  t[1] = Xd[1];
  t[2] = Xd[2];
  TacsScalar t1norm = normalize_vector(&t[0]);

  // t[3:6] = n x t[0:3]
  cross_product(&t[3], n, t);

  // t[6:9] = n
  t[6] = n[0];
  t[7] = n[1];
  t[8] = n[2];

  transform_jacobian_product(tx, ztx, t, Jinv, zJinv);

  // Now, go through and compute the sensitivities of these calculations
  // to the
  for (int k = 0; k < num_nodes; k++) {
    for (int kk = 0; kk < 3; kk++) {  // Cycle through the dimensions
      // dXd[kk] = Na[k]; dXd[3+kk] = Nb[k];
      // dXdd[kk] = Naa[k]; dXdd[3+kk] = Nab[k]; dXdd[6+kk] = Nbb[k];

      // The un-normalized vectors
      TacsScalar drn[3], drn_xi[3], drn_eta[3];

      // Compute the sensitivity of n:
      // n[0] = Xd[1]*Xd[5] - Xd[2]*Xd[4];
      // n[1] = Xd[2]*Xd[3] - Xd[0]*Xd[5];
      // n[2] = Xd[0]*Xd[4] - Xd[1]*Xd[3];

      if (kk == 0) {
        drn[0] = 0.0;
        drn[1] = Xd[2] * Nb[k] - Na[k] * Xd[5];
        drn[2] = Na[k] * Xd[4] - Xd[1] * Nb[k];
      } else if (kk == 1) {
        drn[0] = Na[k] * Xd[5] - Xd[2] * Nb[k];
        drn[1] = 0.0;
        drn[2] = Xd[0] * Nb[k] - Na[k] * Xd[3];
      } else {
        drn[0] = Xd[1] * Nb[k] - Na[k] * Xd[4];
        drn[1] = Na[k] * Xd[3] - Xd[0] * Nb[k];
        drn[2] = 0.0;
      }

      normalize_vector_sens(dn, norm, rn, drn);

      // Compute the sensitivity of n_xi
      // n_xi[0] = Xd[1]*Xdd[5] + Xdd[1]*Xd[5] - Xd[2]*Xdd[4] - Xdd[2]*Xd[4];
      // n_xi[1] = Xd[2]*Xdd[3] + Xdd[2]*Xd[3] - Xd[0]*Xdd[5] - Xdd[0]*Xd[5];
      // n_xi[2] = Xd[0]*Xdd[4] + Xdd[0]*Xd[4] - Xd[1]*Xdd[3] - Xdd[1]*Xd[3];

      if (kk == 0) {
        drn_xi[0] = 0.0;
        drn_xi[1] =
            Xd[2] * Nab[k] + Xdd[2] * Nb[k] - Na[k] * Xdd[5] - Naa[k] * Xd[5];
        drn_xi[2] =
            Na[k] * Xdd[4] + Naa[k] * Xd[4] - Xd[1] * Nab[k] - Xdd[1] * Nb[k];
      } else if (kk == 1) {
        drn_xi[0] =
            Na[k] * Xdd[5] + Naa[k] * Xd[5] - Xd[2] * Nab[k] - Xdd[2] * Nb[k];
        drn_xi[1] = 0.0;
        drn_xi[2] =
            Xd[0] * Nab[k] + Xdd[0] * Nb[k] - Na[k] * Xdd[3] - Naa[k] * Xd[3];
      } else {
        drn_xi[0] =
            Xd[1] * Nab[k] + Xdd[1] * Nb[k] - Na[k] * Xdd[4] - Naa[k] * Xd[4];
        drn_xi[1] =
            Na[k] * Xdd[3] + Naa[k] * Xd[3] - Xd[0] * Nab[k] - Xdd[0] * Nb[k];
        drn_xi[2] = 0.0;
      }

      // Compute the sensitivity of the normalization
      normalize_vector_second_sens(dn_xi, norm, rn, rn_xi, drn, drn_xi);

      // Compute the sensitivity of n_eta
      // n_eta[0] = Xd[1]*Xdd[8] + Xdd[4]*Xd[5] - Xd[2]*Xdd[7] - Xdd[5]*Xd[4];
      // n_eta[1] = Xd[2]*Xdd[6] + Xdd[5]*Xd[3] - Xd[0]*Xdd[8] - Xdd[3]*Xd[5];
      // n_eta[2] = Xd[0]*Xdd[7] + Xdd[3]*Xd[4] - Xd[1]*Xdd[6] - Xdd[4]*Xd[3];

      if (kk == 0) {
        drn_eta[0] = 0.0;
        drn_eta[1] =
            Xd[2] * Nbb[k] + Xdd[5] * Nb[k] - Na[k] * Xdd[8] - Nab[k] * Xd[5];
        drn_eta[2] =
            Na[k] * Xdd[7] + Nab[k] * Xd[4] - Xd[1] * Nbb[k] - Xdd[4] * Nb[k];
      } else if (kk == 1) {
        drn_eta[0] =
            Na[k] * Xdd[8] + Nab[k] * Xd[5] - Xd[2] * Nbb[k] - Xdd[5] * Nb[k];
        drn_eta[1] = 0.0;
        drn_eta[2] =
            Xd[0] * Nbb[k] + Xdd[3] * Nb[k] - Na[k] * Xdd[6] - Nab[k] * Xd[3];
      } else {
        drn_eta[0] =
            Xd[1] * Nbb[k] + Xdd[4] * Nb[k] - Na[k] * Xdd[7] - Nab[k] * Xd[4];
        drn_eta[1] =
            Na[k] * Xdd[6] + Nab[k] * Xd[3] - Xd[0] * Nbb[k] - Xdd[3] * Nb[k];
        drn_eta[2] = 0.0;
      }

      // Compute the sensitivity of the normalization
      normalize_vector_second_sens(dn_eta, norm, rn, rn_eta, drn, drn_eta);

      // Compute the sensitivity of the inverse operation
      TacsScalar dXd[9];
      dXd[0] = dXd[1] = dXd[2] = 0.0;
      dXd[3] = dXd[4] = dXd[5] = 0.0;
      dXd[kk] = Na[k];
      dXd[3 + kk] = Nb[k];

      dXd[6] = dn[0];
      dXd[7] = dn[1];
      dXd[8] = dn[2];

      TacsScalar dJinv[9], dzJinv[9];
      FElibrary::jacobian3dSens(Xd, dXd, h, Jinv, dJinv, &dh[0]);

      compute_jacobian_zeta_sens(dzJinv, Jinv, dJinv, n_xi, dn_xi, n_eta,
                                 dn_eta);

      // Compute the transformation
      // t[0:3] = Xd[0:3]//||Xd[0:3]||
      TacsScalar scale = -(Xd[kk] * Na[k]) / (t1norm * t1norm * t1norm);
      dt[0] = scale * Xd[0];
      dt[1] = scale * Xd[1];
      dt[2] = scale * Xd[2];

      dt[kk] += Na[k] / t1norm;

      // t[3:6] = n x t[0:3]
      cross_product_sens(&dt[3], n, t, dn, dt);

      // t[6:9] = n
      dt[6] = dn[0];
      dt[7] = dn[1];
      dt[8] = dn[2];

      transform_jacobian_product_sens(dtx, dztx, t, dt, Jinv, zJinv, dJinv,
                                      dzJinv);

      dt += 9, dtx += 9, dztx += 9;
      dn += 3, dn_xi += 3, dn_eta += 3;
      dh++;
    }
  }

  return h;
}

/*
  Compute the transformation that is required to take the derivatives
  w.r.t. the local shell coordinates (xi,eta,zeta) to the local shell
  frame (e1,e2,e3).

  This transformation involves two steps:

  U_{x,x} = U_{,xi} * X_{,xi}^{-1}
  U_{x,e} = U_{,xi} * X_{,xi}^{-1} * T^{T}
  = U_{,xi} * tx^{T}

  tx^{T} = X_{,xi}^{-1} * T^{T}
  tx = T * X_{,xi}^{-T}
  = T * Jinv

  input:
  Xd  == The derivatives of the shell coordinates
  Xdd == The second derivatives of the shell coordinates

  output:
  n     == The shell normal direction
  n_xi  == The derivative of n w.r.t. xi
  n_eta == The derivative of n w.r.t. eta
  tx    == The transformation
  ztx   == The derivative of the transformation w.r.t. zeta

  intermediaries:
  Jinv == X_{,xi}^{-T}, the inverse of the Jacobian
  t    == The transformation from X to e
*/
TacsScalar compute_transform_refaxis(TacsScalar t[], TacsScalar tx[],
                                     TacsScalar ztx[], TacsScalar n[],
                                     TacsScalar n_xi[], TacsScalar n_eta[],
                                     const TacsScalar axis[], TacsScalar Xd[],
                                     const TacsScalar Xdd[]) {
  // n = Xd[0] x Xd[3]
  // |   i     j      k    |
  // | Xd[0]  Xd[1]  Xd[2] |
  // | Xd[3]  Xd[4]  Xd[5] |
  cross_product(n, &Xd[0], &Xd[3]);
  cross_product_sens(n_xi, &Xd[0], &Xd[3], &Xdd[0], &Xdd[3]);
  cross_product_sens(n_eta, &Xd[0], &Xd[3], &Xdd[3], &Xdd[6]);

  // Normalize the normal vector
  TacsScalar norm = normalize_vector(n);

  // Now take the derivative of the normal vector along
  // the directions xi and eta
  normalize_vector_sens(n_xi, norm, n);
  normalize_vector_sens(n_eta, norm, n);

  Xd[6] = n[0];
  Xd[7] = n[1];
  Xd[8] = n[2];

  TacsScalar Jinv[9], zJinv[9];
  TacsScalar h = FElibrary::jacobian3d(Xd, Jinv);
  compute_jacobian_zeta(zJinv, Jinv, n_xi, n_eta);

  // Compute the transformation. Take the component
  // of the reference axis in the plane of the surface.
  TacsScalar an = n[0] * axis[0] + n[1] * axis[1] + n[2] * axis[2];
  t[0] = axis[0] - an * n[0];
  t[1] = axis[1] - an * n[1];
  t[2] = axis[2] - an * n[2];
  // t[0:3] = axis - (n, axis) n
  // t = t[0:3]/||t[0:3]||
  normalize_vector(&t[0]);

  // t[3:6] = n x t[0:3]
  cross_product(&t[3], n, t);

  // t[6:9] = n
  t[6] = n[0];
  t[7] = n[1];
  t[8] = n[2];

  transform_jacobian_product(tx, ztx, t, Jinv, zJinv);

  return h;
}

/*
  Compute the transformation that is required to take the derivatives
  w.r.t. the local shell coordinates (xi,eta,zeta) to the local shell
  frame (e1,e2,e3).
*/
TacsScalar compute_transform_refaxis_sens(
    TacsScalar *dh, TacsScalar t[], TacsScalar dt[], TacsScalar tx[],
    TacsScalar dtx[], TacsScalar ztx[], TacsScalar dztx[], TacsScalar n[],
    TacsScalar dn[], TacsScalar n_xi[], TacsScalar dn_xi[], TacsScalar n_eta[],
    TacsScalar dn_eta[], const TacsScalar axis[], TacsScalar Xd[],
    const TacsScalar Xdd[], const double Na[], const double Nb[],
    const double Naa[], const double Nab[], const double Nbb[], int num_nodes) {
  // n = Xd[0] x Xd[3]
  // |   i     j      k    |
  // | Xd[0]  Xd[1]  Xd[2] |
  // | Xd[3]  Xd[4]  Xd[5] |
  TacsScalar rn[3], rn_xi[3], rn_eta[3];
  cross_product(rn, &Xd[0], &Xd[3]);
  cross_product_sens(rn_xi, &Xd[0], &Xd[3], &Xdd[0], &Xdd[3]);
  cross_product_sens(rn_eta, &Xd[0], &Xd[3], &Xdd[3], &Xdd[6]);

  // Normalize the normal vector
  TacsScalar norm = normalize_vector(n, rn);

  // Now take the derivative of the normal vector along
  // the directions xi and eta
  normalize_vector_sens(n_xi, norm, rn, rn_xi);
  normalize_vector_sens(n_eta, norm, rn, rn_eta);

  Xd[6] = n[0];
  Xd[7] = n[1];
  Xd[8] = n[2];

  TacsScalar Jinv[9], zJinv[9];
  TacsScalar h = FElibrary::jacobian3d(Xd, Jinv);
  compute_jacobian_zeta(zJinv, Jinv, n_xi, n_eta);

  // Compute the transformation. Take the component
  // of the reference axis in the plane of the surface.
  TacsScalar an = n[0] * axis[0] + n[1] * axis[1] + n[2] * axis[2];
  t[0] = axis[0] - an * n[0];
  t[1] = axis[1] - an * n[1];
  t[2] = axis[2] - an * n[2];
  // t[0:3] = axis - (n, axis) n
  // t = t[0:3]/||t[0:3]||
  TacsScalar t1norm = normalize_vector(&t[0]);

  // t[3:6] = n x t[0:3]
  cross_product(&t[3], n, t);

  // t[6:9] = n
  t[6] = n[0];
  t[7] = n[1];
  t[8] = n[2];

  transform_jacobian_product(tx, ztx, t, Jinv, zJinv);

  // Now, go through and compute the sensitivities of these calculations
  // to the
  for (int k = 0; k < num_nodes; k++) {
    for (int kk = 0; kk < 3; kk++) {  // Cycle through the dimensions
      // dXd[kk] = Na[k]; dXd[3+kk] = Nb[k];
      // dXdd[kk] = Naa[k]; dXdd[3+kk] = Nab[k]; dXdd[6+kk] = Nbb[k];

      // The un-normalized vectors
      TacsScalar drn[3], drn_xi[3], drn_eta[3];

      // Compute the sensitivity of n:
      // n[0] = Xd[1]*Xd[5] - Xd[2]*Xd[4];
      // n[1] = Xd[2]*Xd[3] - Xd[0]*Xd[5];
      // n[2] = Xd[0]*Xd[4] - Xd[1]*Xd[3];

      if (kk == 0) {
        drn[0] = 0.0;
        drn[1] = Xd[2] * Nb[k] - Na[k] * Xd[5];
        drn[2] = Na[k] * Xd[4] - Xd[1] * Nb[k];
      } else if (kk == 1) {
        drn[0] = Na[k] * Xd[5] - Xd[2] * Nb[k];
        drn[1] = 0.0;
        drn[2] = Xd[0] * Nb[k] - Na[k] * Xd[3];
      } else if (kk == 2) {
        drn[0] = Xd[1] * Nb[k] - Na[k] * Xd[4];
        drn[1] = Na[k] * Xd[3] - Xd[0] * Nb[k];
        drn[2] = 0.0;
      }

      normalize_vector_sens(dn, norm, rn, drn);

      // Compute the sensitivity of n_xi
      // n_xi[0] = Xd[1]*Xdd[5] + Xdd[1]*Xd[5] - Xd[2]*Xdd[4] - Xdd[2]*Xd[4];
      // n_xi[1] = Xd[2]*Xdd[3] + Xdd[2]*Xd[3] - Xd[0]*Xdd[5] - Xdd[0]*Xd[5];
      // n_xi[2] = Xd[0]*Xdd[4] + Xdd[0]*Xd[4] - Xd[1]*Xdd[3] - Xdd[1]*Xd[3];

      if (kk == 0) {
        drn_xi[0] = 0.0;
        drn_xi[1] =
            Xd[2] * Nab[k] + Xdd[2] * Nb[k] - Na[k] * Xdd[5] - Naa[k] * Xd[5];
        drn_xi[2] =
            Na[k] * Xdd[4] + Naa[k] * Xd[4] - Xd[1] * Nab[k] - Xdd[1] * Nb[k];
      } else if (kk == 1) {
        drn_xi[0] =
            Na[k] * Xdd[5] + Naa[k] * Xd[5] - Xd[2] * Nab[k] - Xdd[2] * Nb[k];
        drn_xi[1] = 0.0;
        drn_xi[2] =
            Xd[0] * Nab[k] + Xdd[0] * Nb[k] - Na[k] * Xdd[3] - Naa[k] * Xd[3];
      } else if (kk == 2) {
        drn_xi[0] =
            Xd[1] * Nab[k] + Xdd[1] * Nb[k] - Na[k] * Xdd[4] - Naa[k] * Xd[4];
        drn_xi[1] =
            Na[k] * Xdd[3] + Naa[k] * Xd[3] - Xd[0] * Nab[k] - Xdd[0] * Nb[k];
        drn_xi[2] = 0.0;
      }

      // Compute the sensitivity of the normalization
      normalize_vector_second_sens(dn_xi, norm, rn, rn_xi, drn, drn_xi);

      // Compute the sensitivity of n_eta
      // n_eta[0] = Xd[1]*Xdd[8] + Xdd[4]*Xd[5] - Xd[2]*Xdd[7] - Xdd[5]*Xd[4];
      // n_eta[1] = Xd[2]*Xdd[6] + Xdd[5]*Xd[3] - Xd[0]*Xdd[8] - Xdd[3]*Xd[5];
      // n_eta[2] = Xd[0]*Xdd[7] + Xdd[3]*Xd[4] - Xd[1]*Xdd[6] - Xdd[4]*Xd[3];

      if (kk == 0) {
        drn_eta[0] = 0.0;
        drn_eta[1] =
            Xd[2] * Nbb[k] + Xdd[5] * Nb[k] - Na[k] * Xdd[8] - Nab[k] * Xd[5];
        drn_eta[2] =
            Na[k] * Xdd[7] + Nab[k] * Xd[4] - Xd[1] * Nbb[k] - Xdd[4] * Nb[k];
      } else if (kk == 1) {
        drn_eta[0] =
            Na[k] * Xdd[8] + Nab[k] * Xd[5] - Xd[2] * Nbb[k] - Xdd[5] * Nb[k];
        drn_eta[1] = 0.0;
        drn_eta[2] =
            Xd[0] * Nbb[k] + Xdd[3] * Nb[k] - Na[k] * Xdd[6] - Nab[k] * Xd[3];
      } else if (kk == 2) {
        drn_eta[0] =
            Xd[1] * Nbb[k] + Xdd[4] * Nb[k] - Na[k] * Xdd[7] - Nab[k] * Xd[4];
        drn_eta[1] =
            Na[k] * Xdd[6] + Nab[k] * Xd[3] - Xd[0] * Nbb[k] - Xdd[3] * Nb[k];
        drn_eta[2] = 0.0;
      }

      // Compute the sensitivity of the normalization
      normalize_vector_second_sens(dn_eta, norm, rn, rn_eta, drn, drn_eta);

      // Compute the sensitivity of the inverse operation
      TacsScalar dXd[9];
      dXd[0] = dXd[1] = dXd[2] = 0.0;
      dXd[3] = dXd[4] = dXd[5] = 0.0;
      dXd[kk] = Na[k];
      dXd[3 + kk] = Nb[k];

      dXd[6] = dn[0];
      dXd[7] = dn[1];
      dXd[8] = dn[2];

      TacsScalar dJinv[9], dzJinv[9];
      FElibrary::jacobian3dSens(Xd, dXd, h, Jinv, dJinv, &dh[0]);

      compute_jacobian_zeta_sens(dzJinv, Jinv, dJinv, n_xi, dn_xi, n_eta,
                                 dn_eta);

      // Compute the transformation and its derivative
      TacsScalar dan = dn[0] * axis[0] + dn[1] * axis[1] + dn[2] * axis[2];
      dt[0] = -an * dn[0] - dan * n[0];
      dt[1] = -an * dn[1] - dan * n[1];
      dt[2] = -an * dn[2] - dan * n[2];

      normalize_vector_sens(&dt[0], t1norm, &t[0]);

      // t[3:6] = n x t[0:3]
      cross_product_sens(&dt[3], n, t, dn, dt);

      // t[6:9] = n
      dt[6] = dn[0];
      dt[7] = dn[1];
      dt[8] = dn[2];

      transform_jacobian_product_sens(dtx, dztx, t, dt, Jinv, zJinv, dJinv,
                                      dzJinv);

      dt += 9, dtx += 9, dztx += 9;
      dn += 3, dn_xi += 3, dn_eta += 3;
      dh++;
    }
  }

  return h;
}

/*
  Transform a 3d vector in place from one reference frame
  to another.
*/

static inline void transform_vector3d(TacsScalar vec[], const TacsScalar t[]) {
  TacsScalar r[3];
  r[0] = vec[0];
  r[1] = vec[1];
  r[2] = vec[2];

  vec[0] = t[0] * r[0] + t[1] * r[1] + t[2] * r[2];
  vec[1] = t[3] * r[0] + t[4] * r[1] + t[5] * r[2];
  vec[2] = t[6] * r[0] + t[7] * r[1] + t[8] * r[2];
}

static inline void transform_vector3d_sens(TacsScalar vec[],
                                           const TacsScalar r[],
                                           const TacsScalar t[],
                                           const TacsScalar dt[]) {
  TacsScalar dr[3];
  dr[0] = vec[0];
  dr[1] = vec[1];
  dr[2] = vec[2];

  vec[0] = (dt[0] * r[0] + dt[1] * r[1] + dt[2] * r[2] + t[0] * dr[0] +
            t[1] * dr[1] + t[2] * dr[2]);
  vec[1] = (dt[3] * r[0] + dt[4] * r[1] + dt[5] * r[2] + t[3] * dr[0] +
            t[4] * dr[1] + t[5] * dr[2]);
  vec[2] = (dt[6] * r[0] + dt[7] * r[1] + dt[8] * r[2] + t[6] * dr[0] +
            t[7] * dr[1] + t[8] * dr[2]);
}

/*
  Transform the displacement gradient from one reference
  frame to another - note that this only transforms the
  components and not the derivatives.
*/
static inline void transform_displ_gradient(TacsScalar Ued[],
                                            const TacsScalar t[],
                                            const TacsScalar Ud[]) {
  // Transform the displacement gradient
  Ued[0] = t[0] * Ud[0] + t[1] * Ud[2] + t[2] * Ud[4];
  Ued[2] = t[3] * Ud[0] + t[4] * Ud[2] + t[5] * Ud[4];
  Ued[4] = t[6] * Ud[0] + t[7] * Ud[2] + t[8] * Ud[4];

  Ued[1] = t[0] * Ud[1] + t[1] * Ud[3] + t[2] * Ud[5];
  Ued[3] = t[3] * Ud[1] + t[4] * Ud[3] + t[5] * Ud[5];
  Ued[5] = t[6] * Ud[1] + t[7] * Ud[3] + t[8] * Ud[5];

  // Transform the rotation gradient
  Ued[6] = t[0] * Ud[6] + t[1] * Ud[8] + t[2] * Ud[10];
  Ued[8] = t[3] * Ud[6] + t[4] * Ud[8] + t[5] * Ud[10];
  Ued[10] = t[6] * Ud[6] + t[7] * Ud[8] + t[8] * Ud[10];

  Ued[7] = t[0] * Ud[7] + t[1] * Ud[9] + t[2] * Ud[11];
  Ued[9] = t[3] * Ud[7] + t[4] * Ud[9] + t[5] * Ud[11];
  Ued[11] = t[6] * Ud[7] + t[7] * Ud[9] + t[8] * Ud[11];
}

static inline void transform_displ_gradient_bmat(TacsScalar Ued[], int ii,
                                                 const TacsScalar t[],
                                                 double Na, double Nb) {
  if (ii < 3) {  // Ud[2*ii] = Na, Ud[2*ii+1] = Nb
    // Transform the displacement gradient
    Ued[0] = Na * t[ii];
    Ued[2] = Na * t[3 + ii];
    Ued[4] = Na * t[6 + ii];

    Ued[1] = Nb * t[ii];
    Ued[3] = Nb * t[3 + ii];
    Ued[5] = Nb * t[6 + ii];

    Ued[6] = Ued[7] = Ued[8] = 0.0;
    Ued[9] = Ued[10] = Ued[11] = 0.0;
  } else {
    // Transform the rotation gradient
    Ued[6] = Na * t[ii];
    Ued[8] = Na * t[3 + ii];
    Ued[10] = Na * t[6 + ii];

    Ued[7] = Nb * t[ii];
    Ued[9] = Nb * t[3 + ii];
    Ued[11] = Nb * t[6 + ii];

    Ued[0] = Ued[1] = Ued[2] = 0.0;
    Ued[3] = Ued[4] = Ued[5] = 0.0;
  }
}

/*
  Compute the cross product of phi X n and return the result.  These
  terms are used to calculate the contribution of the rotational
  components to the shell strain.

  Compute the cross product of the rotation vector phi = [ rotx, roty,
  rotz ] with the normal. Also required are the derivatives of the
  cross product with respect to the shell coordinates xi, eta

  r = phi X n =
  | i j k |
  | 3 4 5 |
  | 0 1 2 |

  r_xi = phi_xi X n =
  | i j k  |
  | 6 8 10 |
  | 0 1 2  |

  r_eta = phi_eta X n =
  | i j k  |
  | 7 9 11 |
  | 0 1 2  |
*/
static inline void normal_rot(TacsScalar r[], TacsScalar r_xi[],
                              TacsScalar r_eta[], const TacsScalar U[],
                              const TacsScalar Ud[], const TacsScalar n[],
                              const TacsScalar n_xi[],
                              const TacsScalar n_eta[]) {
  r[0] = U[4] * n[2] - U[5] * n[1];
  r[1] = U[5] * n[0] - U[3] * n[2];
  r[2] = U[3] * n[1] - U[4] * n[0];

  r_xi[0] = U[4] * n_xi[2] - U[5] * n_xi[1] + Ud[8] * n[2] - Ud[10] * n[1];
  r_xi[1] = U[5] * n_xi[0] - U[3] * n_xi[2] + Ud[10] * n[0] - Ud[6] * n[2];
  r_xi[2] = U[3] * n_xi[1] - U[4] * n_xi[0] + Ud[6] * n[1] - Ud[8] * n[0];

  r_eta[0] = U[4] * n_eta[2] - U[5] * n_eta[1] + Ud[9] * n[2] - Ud[11] * n[1];
  r_eta[1] = U[5] * n_eta[0] - U[3] * n_eta[2] + Ud[11] * n[0] - Ud[7] * n[2];
  r_eta[2] = U[3] * n_eta[1] - U[4] * n_eta[0] + Ud[7] * n[1] - Ud[9] * n[0];
}

/*
  Compute the derivative w.r.t. the state variables of the cross
  product of the rotation variables with the normal. The derivative
  along coordinate lines is also returned.
*/
static inline void normal_rot_bmat(TacsScalar r[], TacsScalar r_xi[],
                                   TacsScalar r_eta[], int ii, double N,
                                   double Na, double Nb, const TacsScalar n[],
                                   const TacsScalar n_xi[],
                                   const TacsScalar n_eta[]) {
  if (ii < 3) {
    r[0] = r[1] = r[2] = 0.0;
    r_xi[0] = r_xi[1] = r_xi[2] = 0.0;
    r_eta[0] = r_eta[1] = r_eta[2] = 0.0;
  } else if (ii == 3) {
    r[0] = 0.0;
    r[1] = -N * n[2];
    r[2] = N * n[1];

    r_xi[0] = 0.0;
    r_xi[1] = -N * n_xi[2] - Na * n[2];
    r_xi[2] = N * n_xi[1] + Na * n[1];

    r_eta[0] = 0.0;
    r_eta[1] = -N * n_eta[2] - Nb * n[2];
    r_eta[2] = N * n_eta[1] + Nb * n[1];
  } else if (ii == 4) {
    r[0] = N * n[2];
    r[1] = 0.0;
    r[2] = -N * n[0];

    r_xi[0] = N * n_xi[2] + Na * n[2];
    r_xi[1] = 0.0;
    r_xi[2] = -N * n_xi[0] - Na * n[0];

    r_eta[0] = N * n_eta[2] + Nb * n[2];
    r_eta[1] = 0.0;
    r_eta[2] = -N * n_eta[0] - Nb * n[0];
  } else {
    r[0] = -N * n[1];
    r[1] = N * n[0];
    r[2] = 0.0;

    r_xi[0] = -N * n_xi[1] - Na * n[1];
    r_xi[1] = N * n_xi[0] + Na * n[0];
    r_xi[2] = 0.0;

    r_eta[0] = -N * n_eta[1] - Nb * n[1];
    r_eta[1] = N * n_eta[0] + Nb * n[0];
    r_eta[2] = 0.0;
  }
}

/*
  Compute the complete strain expression in the local coordinates

  The displacements through the thickness are:

  U(zeta) = u0 + zeta * phi X n

  The derivative of the displacements w.r.t. the local shell
  coordinates at the mid-surface are:

  U_{,xi}|zeta=0 = [ u0,xi | u0,eta | phi X n ]

  The derivative of the displacements through the thickness in the
  zeta direction:

  U_{,xi,zeta}|zeta=0 =
  [ phi,xi X n + phi X n_xi |  phi,eta X n + phi X n_eta | 0 ]

  The displacements in the global coordinate frame are:

  Ux = [ u0, v0, w0, rotx, roty, rotz ]

  The derivative of the displacements w.r.t. the local shell
  coordinates are:

  Uxd = [
  u0,xi u0,eta
  v0,xi v0,eta
  w0,xi w0,eta
  rotx,xi rotx,eta
  roty,xi roty,eta
  rotz,xi rotz,eta ]

  With these formula, the following expression for the displacement
  gradient can be constructed:

  U_{e,e} = t * U_{x,xi} * tx^{T}

  Note that tx^{T} is the product of the inverse of the Jacobian and
  the transformation.

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  n: the shell normal
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta

  output:
  strain: the value of the strain at the current point
  rotz: the difference between the drilling rotation and actual rotation
*/
void linear_strain(TacsScalar strain[], TacsScalar *rotz, const TacsScalar Ux[],
                   const TacsScalar Uxd[], const TacsScalar t[],
                   const TacsScalar tx[], const TacsScalar ztx[],
                   const TacsScalar n[], const TacsScalar n_xi[],
                   const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // rotz = 0.5 * ( v,x - u,y ) - inplane_rot
  *rotz = (0.5 * (Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2] -
                  (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5])) -
           (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // Compute the in-plane components of the tensorial strain
  strain[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];  // e_11
  strain[1] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];  // e_22
  strain[2] = (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5] + Ud[2] * tx[0] +
               Ud[3] * tx[1] + r[1] * tx[2]);  // e_12

  // Compute the bending components of the strain
  strain[3] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
               Ud[1] * ztx[1] + r[0] * ztx[2]);  // k_11
  strain[4] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
               Ud[3] * ztx[4] + r[1] * ztx[5]);  // k_22

  strain[5] =
      (r_xi[0] * tx[3] + r_eta[0] * tx[4] + r_xi[1] * tx[0] + r_eta[1] * tx[1] +
       Ud[0] * ztx[3] + Ud[1] * ztx[4] + r[0] * ztx[5] + Ud[2] * ztx[0] +
       Ud[3] * ztx[1] + r[1] * ztx[2]);  // k_12

  // Compute the shear terms
  strain[6] = (Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8] + Ud[4] * tx[3] +
               Ud[5] * tx[4] + r[2] * tx[5]);  // e_23
  strain[7] = (Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8] + Ud[4] * tx[0] +
               Ud[5] * tx[1] + r[2] * tx[2]);  // e_13
}

/*
  Compute the derivative of the strain w.r.t. all displacement
  variables simultaneously.

  input:
  N: the shape functions

  Na, Nb: the derivatives of the shape functions along the parameter
  directions

  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  n: the shell normal
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta

  output:
  B: the derivative of the strain w.r.t. the nodal displacements

  rotz: the derivative of difference between the drilling rotation and
  actual rotation
*/
void linear_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                 const double N[], const double Na[], const double Nb[],
                 const TacsScalar t[], const TacsScalar tx[],
                 const TacsScalar ztx[], const TacsScalar n[],
                 const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  // For each point
  for (int i = 0; i < num_points; i++) {
    // For each displacement component
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar Ud[12];
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      TacsScalar r[3], r_xi[3], r_eta[3];
      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // rotz = 0.5 * (v,x - u,y )
      rotz[0] = 0.5 * (Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2] -
                       (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5]));
      if (ii >= 3) {
        rotz[0] -= N[i] * n[ii - 3];
      }

      // Compute the in-plane components of the tensorial strain
      B[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];  // e_11
      B[1] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];  // e_22
      B[2] = (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5] + Ud[2] * tx[0] +
              Ud[3] * tx[1] + r[1] * tx[2]);  // e_12

      // Compute the bending components of the strain
      B[3] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
              Ud[1] * ztx[1] + r[0] * ztx[2]);  // k_11
      B[4] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
              Ud[3] * ztx[4] + r[1] * ztx[5]);  // k_22

      B[5] =
          (r_xi[0] * tx[3] + r_eta[0] * tx[4] + r_xi[1] * tx[0] +
           r_eta[1] * tx[1] + Ud[0] * ztx[3] + Ud[1] * ztx[4] + r[0] * ztx[5] +
           Ud[2] * ztx[0] + Ud[3] * ztx[1] + r[1] * ztx[2]);  // k_12

      // Compute the shear terms
      B[6] = (Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8] + Ud[4] * tx[3] +
              Ud[5] * tx[4] + r[2] * tx[5]);  // e_23
      B[7] = (Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8] + Ud[4] * tx[0] +
              Ud[5] * tx[1] + r[2] * tx[2]);  // e_13

      B += 8;
      rotz += 1;
    }
  }
}

/*
  Compute the design variable sensitivity of the strain to a
  perturbation of the nodal locations.

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t: the transformation
  dt: the sensitivity of the local transformation

  tx: the transformation times the Jacobian
  dtx: the sensitivity of the transformation times the Jacobian

  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  dztx: the sensitivity of ztx

  n: the shell normal
  dn: the sensitivity of the shell normal

  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  dn_xi: the sensitivity of n_xi

  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta
  dn_eta: the sensitivity of n_eta

  num_components: the number of sensitivity components to compute

  output:
  strain: the value of the strain at the current point
  dstrain: the sentivity of the strain at the curren point

  rotz: the difference between the drilling rotation and actual rotation
  drotz: the sensitivity of rotz
*/
void linear_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                        TacsScalar *rotz, TacsScalar *drotz,
                        const TacsScalar Ux[], const TacsScalar Uxd[],
                        const TacsScalar t[], const TacsScalar dt[],
                        const TacsScalar tx[], const TacsScalar dtx[],
                        const TacsScalar ztx[], const TacsScalar dztx[],
                        const TacsScalar n[], const TacsScalar dn[],
                        const TacsScalar n_xi[], const TacsScalar dn_xi[],
                        const TacsScalar n_eta[], const TacsScalar dn_eta[],
                        int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface values
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
  TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  // rotz = 0.5 * (v,x - u,y ) - inplane_rot
  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // The in plane strains
  strain[0] = ux;
  strain[1] = vy;
  strain[2] = uy + vx;

  // Compute the bending components of the strain
  strain[3] = ux1;
  strain[4] = vy1;
  strain[5] = uy1 + vx1;

  // Compute the shear terms
  strain[6] = vz + wy;
  strain[7] = uz + wx;

  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                      Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
    TacsScalar suz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                      Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
    TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                      Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
    TacsScalar svz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                      Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

    TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                      Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
    TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                      Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    drotz[0] =
        (0.5 * (svx - suy) - (dn[0] * Ux[3] + dn[1] * Ux[4] + dn[2] * Ux[5]));

    // The in plane strains
    dstrain[0] = sux;
    dstrain[1] = svy;
    dstrain[2] = suy + svx;

    // Compute the bending components of the strain
    dstrain[3] = sux1;
    dstrain[4] = svy1;
    dstrain[5] = suy1 + svx1;

    // Compute the shear terms
    dstrain[6] = svz + swy;
    dstrain[7] = suz + swx;

    dstrain += 8;
    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
    drotz++;
  }
}

/*
  Add components of the sensitivity w.r.t. the nodes of the B matrix
  (the derivative of the strain w.r.t. the nodal coordinates) - to the
  residual matrix 'res'. This is used for constructing the element
  residuals.

  input:
  num_points: number of points (number of nodes)
  stress_scale: the scale value for the stresses
  rot_scale: the value to scale the rotation derivative
  stress: the components of the stress

  N, Na, Nb: the shape functions and their derivatives
  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivity of the transformations

  n, n_xi, n_eta: the normal and the derivative of the normal along
  the xi and eta directions
  dn, dn_xi, dn_eta: the sensitivity of the normal (and its derivatives)
*/
void add_linear_bmat_sens(
    TacsScalar res[], int num_points, const TacsScalar stress_scale,
    const TacsScalar rot_scale, const TacsScalar stress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar t[],
    const TacsScalar dt[], const TacsScalar tx[], const TacsScalar dtx[],
    const TacsScalar ztx[], const TacsScalar dztx[], const TacsScalar n[],
    const TacsScalar dn[], const TacsScalar n_xi[], const TacsScalar dn_xi[],
    const TacsScalar n_eta[], const TacsScalar dn_eta[], int num_components) {
  // For each component in the derivative
  for (int k = 0; k < num_components; k++) {
    // For each point
    for (int i = 0; i < num_points; i++) {
      // For each displacement component
      for (int ii = 0; ii < 6; ii++) {
        TacsScalar Ud[12], dUd[12];
        transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
        transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

        TacsScalar r[3], r_xi[3], r_eta[3];
        TacsScalar dr[3], dr_xi[3], dr_eta[3];
        normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

        normal_rot_bmat(dr, dr_xi, dr_eta, ii, N[i], Na[i], Nb[i], dn, dn_xi,
                        dn_eta);

        transform_vector3d_sens(dr, r, t, dt);
        transform_vector3d_sens(dr_xi, r_xi, t, dt);
        transform_vector3d_sens(dr_eta, r_eta, t, dt);

        transform_vector3d(r, t);
        transform_vector3d(r_xi, t);
        transform_vector3d(r_eta, t);

        // The mid-surface values
        TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                          Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
        TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                          Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
        TacsScalar suz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                          Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

        TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                          Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
        TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                          Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
        TacsScalar svz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                          Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

        TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                          Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
        TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                          Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);

        // The first-derivative values of the displacements
        TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] +
                           dUd[0] * ztx[0] + dUd[1] * ztx[1] + dr[0] * ztx[2] +
                           r_xi[0] * dtx[0] + r_eta[0] * dtx[1] +
                           Ud[0] * dztx[0] + Ud[1] * dztx[1] + r[0] * dztx[2]);
        TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] +
                           dUd[0] * ztx[3] + dUd[1] * ztx[4] + dr[0] * ztx[5] +
                           r_xi[0] * dtx[3] + r_eta[0] * dtx[4] +
                           Ud[0] * dztx[3] + Ud[1] * dztx[4] + r[0] * dztx[5]);

        TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] +
                           dUd[2] * ztx[0] + dUd[3] * ztx[1] + dr[1] * ztx[2] +
                           r_xi[1] * dtx[0] + r_eta[1] * dtx[1] +
                           Ud[2] * dztx[0] + Ud[3] * dztx[1] + r[1] * dztx[2]);
        TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] +
                           dUd[2] * ztx[3] + dUd[3] * ztx[4] + dr[1] * ztx[5] +
                           r_xi[1] * dtx[3] + r_eta[1] * dtx[4] +
                           Ud[2] * dztx[3] + Ud[3] * dztx[4] + r[1] * dztx[5]);

        // Add the contributions from the stress sensitivities
        res[0] +=
            stress_scale *
            (stress[0] * sux + stress[1] * svy + stress[2] * (suy + svx) +
             stress[3] * sux1 + stress[4] * svy1 + stress[5] * (suy1 + svx1) +
             stress[6] * (svz + swy) + stress[7] * (suz + swx));

        // Add the contribution from the rotational sensitivities
        res[0] += rot_scale * 0.5 * (svx - suy);

        if (ii >= 3) {
          res[0] -= rot_scale * N[i] * dn[ii - 3];
        }

        res++;
      }
    }

    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
  }
}

/*
  Compute the nonlinear components of the strain.

  This function uses the same expression for the displacement
  gradient:

  U_{e,e} = t * U_{x,xi} * tx^{T}

  These are then used to construct the nonlinear expressions for the
  strain in the local Cartesian coordinate frame e.

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  n: the shell normal
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta

  output:
  strain: the value of the strain at the current point
  rotz: the difference between the drilling rotation and actual rotation
*/
void nonlinear_strain(TacsScalar strain[], TacsScalar *rotz,
                      const TacsScalar Ux[], const TacsScalar Uxd[],
                      const TacsScalar t[], const TacsScalar tx[],
                      const TacsScalar ztx[], const TacsScalar n[],
                      const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
  TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
  TacsScalar wz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // The in plane strains
  strain[0] = ux + 0.5 * (ux * ux + vx * vx + wx * wx);
  strain[1] = vy + 0.5 * (uy * uy + vy * vy + wy * wy);
  strain[2] = uy + vx + ux * uy + vx * vy + wx * wy;

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);

  // Compute the shear terms
  strain[6] = vz + wy + uz * uy + vz * vy + wz * wy;
  strain[7] = uz + wx + uz * ux + vz * vx + wz * wx;
}

/*
  Compute the derivative of the nonlinear strain expressions w.r.t
  the nodal displacements/rotations.

  input:
  N: the shape functions
  Na, Nb: the derivatives of the shape functions w.r.t. the natural
  coordinates of the shell

  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  n: the shell normal
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta

  output:
  strain: the value of the strain at the current point
  rotz: the difference between the drilling rotation and actual rotation
*/
void nonlinear_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                    const double N[], const double Na[], const double Nb[],
                    const TacsScalar Ux[], const TacsScalar Uxd[],
                    const TacsScalar t[], const TacsScalar tx[],
                    const TacsScalar ztx[], const TacsScalar n[],
                    const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
  TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
  TacsScalar wz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // For each point
  for (int i = 0; i < num_points; i++) {
    // For each displacement component
    for (int ii = 0; ii < 6; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
      TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
      TacsScalar duz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

      TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
      TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
      TacsScalar dvz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

      TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
      TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
      TacsScalar dwz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

      // The first-derivative values of the displacements
      TacsScalar dux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                         Ud[1] * ztx[1] + r[0] * ztx[2]);
      TacsScalar duy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                         Ud[1] * ztx[4] + r[0] * ztx[5]);

      TacsScalar dvx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                         Ud[3] * ztx[1] + r[1] * ztx[2]);
      TacsScalar dvy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                         Ud[3] * ztx[4] + r[1] * ztx[5]);

      TacsScalar dwx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                         Ud[5] * ztx[1] + r[2] * ztx[2]);
      TacsScalar dwy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                         Ud[5] * ztx[4] + r[2] * ztx[5]);

      rotz[0] = 0.5 * (dvx - duy);
      if (ii >= 3) {
        rotz[0] -= N[i] * n[ii - 3];
      }

      // The in plane strains
      B[0] = dux + (ux * dux + vx * dvx + wx * dwx);
      B[1] = dvy + (uy * duy + vy * dvy + wy * dwy);
      B[2] = duy + dvx +
             (dux * uy + dvx * vy + dwx * wy + ux * duy + vx * dvy + wx * dwy);

      // Compute the bending components of the strain
      B[3] = dux1 + (dux * ux1 + dvx * vx1 + dwx * wx1 + ux * dux1 + vx * dvx1 +
                     wx * dwx1);
      B[4] = dvy1 + (duy * uy1 + dvy * vy1 + dwy * wy1 + uy * duy1 + vy * dvy1 +
                     wy * dwy1);
      B[5] = duy1 + dvx1 +
             (dux1 * uy + dvx1 * vy + dwx1 * wy + dux * uy1 + dvx * vy1 +
              dwx * wy1 + ux1 * duy + vx1 * dvy + wx1 * dwy + ux * duy1 +
              vx * dvy1 + wx * dwy1);

      // Compute the shear terms
      B[6] = dvz + dwy +
             (duz * uy + dvz * vy + dwz * wy + uz * duy + vz * dvy + wz * dwy);
      B[7] = duz + dwx +
             (duz * ux + dvz * vx + dwz * wx + uz * dux + vz * dvx + wz * dwx);

      B += 8;
      rotz += 1;
    }
  }
}

/*
  Add the product of the stress with the second derivative of the
  nonlinear strain multiplied by stress provided.

  Note that this only fills in the bottom half of the matrix.

  input:
  N: the shape functions
  Na, Nb: the derivatives of the shape functions w.r.t. the natural
  coordinates of the shell

  t: the transformation
  tx: the transformation times the Jacobian
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  n: the shell normal
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta
  num_points: the number of points (nodes)
  scale: a scaling factor applied to each new value added to the matrix

  output:
  matrix: the matrix of values
*/
void nonlinear_stress_bmat(TacsScalar matrix[], int num_points,
                           TacsScalar scale, const TacsScalar stress[],
                           const double N[], const double Na[],
                           const double Nb[], const TacsScalar t[],
                           const TacsScalar tx[], const TacsScalar ztx[],
                           const TacsScalar n[], const TacsScalar n_xi[],
                           const TacsScalar n_eta[]) {
  if (num_points > MAX_NUM_NODES) {
    return;
  }

  TacsScalar dU0[9 * 6 * MAX_NUM_NODES];
  TacsScalar dU1[6 * 6 * MAX_NUM_NODES];

  TacsScalar *du0 = dU0;
  TacsScalar *du1 = dU1;

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar Ud[12];
      TacsScalar r[3], r_xi[3], r_eta[3];
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // The mid-surface value of the displacement derivatives
      // ux, uy, uz
      du0[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
      du0[1] = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
      du0[2] = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

      // vx, vy, vz
      du0[3] = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
      du0[4] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
      du0[5] = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

      // wx, wy, wz
      du0[6] = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
      du0[7] = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
      du0[8] = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];
      du0 += 9;

      // The first-derivative values of the displacements
      // u1x, u1y
      du1[0] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                Ud[1] * ztx[1] + r[0] * ztx[2]);
      du1[1] = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                Ud[1] * ztx[4] + r[0] * ztx[5]);

      // v1x, v1y
      du1[2] = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                Ud[3] * ztx[1] + r[1] * ztx[2]);
      du1[3] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                Ud[3] * ztx[4] + r[1] * ztx[5]);

      // w1x, w1y
      du1[4] = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                Ud[5] * ztx[1] + r[2] * ztx[2]);
      du1[5] = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                Ud[5] * ztx[4] + r[2] * ztx[5]);
      du1 += 6;
    }
  }

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar *mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Set the i-counter from the beginning of the row
      TacsScalar *di0 = &dU0[9 * (6 * i + ii)];
      TacsScalar *di1 = &dU1[6 * (6 * i + ii)];

      // Start the j-counter from the i-counter
      TacsScalar *dj0 = dU0;
      TacsScalar *dj1 = dU1;

      // For each point
      for (int j = 0; j <= i; j++) {
        int end = (j == i ? ii : 5);

        // For each displacement component
        for (int jj = 0; jj <= end; jj++) {
          TacsScalar B[8];

          // The in plane strains
          B[0] = (di0[0] * dj0[0] + di0[3] * dj0[3] + di0[6] * dj0[6]);
          B[1] = (di0[1] * dj0[1] + di0[4] * dj0[4] + di0[7] * dj0[7]);
          B[2] = (di0[0] * dj0[1] + di0[3] * dj0[4] + di0[6] * dj0[7] +
                  di0[1] * dj0[0] + di0[4] * dj0[3] + di0[7] * dj0[6]);

          // Compute the bending components of the strain
          B[3] = (di0[0] * dj1[0] + di0[3] * dj1[2] + di0[6] * dj1[4] +
                  di1[0] * dj0[0] + di1[2] * dj0[3] + di1[4] * dj0[6]);
          B[4] = (di0[1] * dj1[1] + di0[4] * dj1[3] + di0[7] * dj1[5] +
                  di1[1] * dj0[1] + di1[3] * dj0[4] + di1[5] * dj0[7]);
          B[5] = (di0[0] * dj1[1] + di0[3] * dj1[3] + di0[6] * dj1[5] +
                  di1[0] * dj0[1] + di1[2] * dj0[4] + di1[4] * dj0[7] +
                  dj0[0] * di1[1] + dj0[3] * di1[3] + dj0[6] * di1[5] +
                  dj1[0] * di0[1] + dj1[2] * di0[4] + dj1[4] * di0[7]);

          // Compute the shear terms
          B[6] = (di0[2] * dj0[1] + di0[5] * dj0[4] + di0[8] * dj0[7] +
                  di0[1] * dj0[2] + di0[4] * dj0[5] + di0[7] * dj0[8]);
          B[7] = (di0[2] * dj0[0] + di0[5] * dj0[3] + di0[8] * dj0[6] +
                  di0[0] * dj0[2] + di0[3] * dj0[5] + di0[6] * dj0[8]);

          mat[0] +=
              scale * (B[0] * stress[0] + B[1] * stress[1] + B[2] * stress[2] +
                       B[3] * stress[3] + B[4] * stress[4] + B[5] * stress[5] +
                       B[6] * stress[6] + B[7] * stress[7]);
          mat += 1;

          dj0 += 9;
          dj1 += 6;
        }
      }
    }
  }
}

/*
  Compute the sensitivities of the nonlinear strain expressions

  input:
  N: the shape functions
  Na, Nb: the derivatives of the shape functions w.r.t. the natural
  coordinates of the shell

  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivities of the transformations

  n, n_xi, n_eta: the shell normals and derivatives
  dn, dn_xi, dn_eta: the sensitivities of the shell normals and their
  derivatives

  output:
  strain: the nonlinear strain
  dstrain: the sensitivity of the nonlinear strain
  rotz: the difference between drilling rotation and out-of-plane rotation
  drotz: the sensitivity of rotz
*/
void nonlinear_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                           TacsScalar *rotz, TacsScalar *drotz,
                           const TacsScalar Ux[], const TacsScalar Uxd[],
                           const TacsScalar t[], const TacsScalar dt[],
                           const TacsScalar tx[], const TacsScalar dtx[],
                           const TacsScalar ztx[], const TacsScalar dztx[],
                           const TacsScalar n[], const TacsScalar dn[],
                           const TacsScalar n_xi[], const TacsScalar dn_xi[],
                           const TacsScalar n_eta[], const TacsScalar dn_eta[],
                           int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface values
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
  TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
  TacsScalar wz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // rotz = 0.5 * (v,x - u,y) - inplane_rot
  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // The in plane strains
  strain[0] = ux + 0.5 * (ux * ux + vx * vx + wx * wx);
  strain[1] = vy + 0.5 * (uy * uy + vy * vy + wy * wy);
  strain[2] = uy + vx + ux * uy + vx * vy + wx * wy;

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);

  // Compute the shear terms
  strain[6] = vz + wy + uz * uy + vz * vy + wz * wy;
  strain[7] = uz + wx + uz * ux + vz * vx + wz * wx;

  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                      Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
    TacsScalar suz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                      Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
    TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                      Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
    TacsScalar svz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                      Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

    TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                      Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
    TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                      Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);
    TacsScalar swz = (dUd[4] * tx[6] + dUd[5] * tx[7] + dr[2] * tx[8] +
                      Ud[4] * dtx[6] + Ud[5] * dtx[7] + r[2] * dtx[8]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    TacsScalar swx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dUd[4] * ztx[0] +
                       dUd[5] * ztx[1] + dr[2] * ztx[2] + r_xi[2] * dtx[0] +
                       r_eta[2] * dtx[1] + Ud[4] * dztx[0] + Ud[5] * dztx[1] +
                       r[2] * dztx[2]);
    TacsScalar swy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dUd[4] * ztx[3] +
                       dUd[5] * ztx[4] + dr[2] * ztx[5] + r_xi[2] * dtx[3] +
                       r_eta[2] * dtx[4] + Ud[4] * dztx[3] + Ud[5] * dztx[4] +
                       r[2] * dztx[5]);

    drotz[0] =
        (0.5 * (svx - suy) - (dn[0] * Ux[3] + dn[1] * Ux[4] + dn[2] * Ux[5]));

    // The in plane strains
    dstrain[0] = sux + (ux * sux + vx * svx + wx * swx);
    dstrain[1] = svy + (uy * suy + vy * svy + wy * swy);
    dstrain[2] =
        suy + svx +
        (sux * uy + svx * vy + swx * wy + ux * suy + vx * svy + wx * swy);

    // Compute the bending components of the strain
    dstrain[3] = sux1 + (sux * ux1 + svx * vx1 + swx * wx1 + ux * sux1 +
                         vx * svx1 + wx * swx1);
    dstrain[4] = svy1 + (suy * uy1 + svy * vy1 + swy * wy1 + uy * suy1 +
                         vy * svy1 + wy * swy1);
    dstrain[5] =
        suy1 + svx1 +
        (sux1 * uy + svx1 * vy + swx1 * wy + sux * uy1 + svx * vy1 + swx * wy1 +
         ux1 * suy + vx1 * svy + wx1 * swy + ux * suy1 + vx * svy1 + wx * swy1);

    // Compute the shear terms
    dstrain[6] =
        svz + swy +
        (suz * uy + svz * vy + swz * wy + uz * suy + vz * svy + wz * swy);
    dstrain[7] =
        suz + swx +
        (suz * ux + svz * vx + swz * wx + uz * sux + vz * svx + wz * swx);

    dstrain += 8;
    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
    drotz++;
  }
}

/*
  Compute the sensitivity of the derivative of the strain w.r.t. the
  displacements and rotations (the sensitivity of B-mat), multiply it
  by the stresses and add it to res[].

  input:
  num_points: number of points (nodes)
  stress_scale: scale all the stresses
  rot_scale: scale the rotation contribution
  stress: the values of the stresses

  N, Na, Nb: the shape functions
  Ux, Uxd: the displacements and the derivatives of the displacements

  t, tx, ztx: the transformation
  dt, dtx, dztx: the sensitivity of the transformation

  n, n_xi, n_eta: the normal of the shell
  dn, dn_xi, dn_eta: the derivative of the normal of the shell

  output:
  res: add the derivative of the product of the stress and B to res
*/
void add_nonlinear_bmat_sens(
    TacsScalar res[], int num_points, const TacsScalar stress_scale,
    const TacsScalar rot_scale, const TacsScalar stress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar Ux[],
    const TacsScalar Uxd[], const TacsScalar t[], const TacsScalar dt[],
    const TacsScalar tx[], const TacsScalar dtx[], const TacsScalar ztx[],
    const TacsScalar dztx[], const TacsScalar n[], const TacsScalar dn[],
    const TacsScalar n_xi[], const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rotational contributions
  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
  TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
  TacsScalar wz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // For each component in the derivative
  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(Ud, t, Uxd);
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                      Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
    TacsScalar suz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                      Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
    TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                      Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
    TacsScalar svz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                      Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

    TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                      Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
    TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                      Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);
    TacsScalar swz = (dUd[4] * tx[6] + dUd[5] * tx[7] + dr[2] * tx[8] +
                      Ud[4] * dtx[6] + Ud[5] * dtx[7] + r[2] * dtx[8]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    TacsScalar swx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dUd[4] * ztx[0] +
                       dUd[5] * ztx[1] + dr[2] * ztx[2] + r_xi[2] * dtx[0] +
                       r_eta[2] * dtx[1] + Ud[4] * dztx[0] + Ud[5] * dztx[1] +
                       r[2] * dztx[2]);
    TacsScalar swy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dUd[4] * ztx[3] +
                       dUd[5] * ztx[4] + dr[2] * ztx[5] + r_xi[2] * dtx[3] +
                       r_eta[2] * dtx[4] + Ud[4] * dztx[3] + Ud[5] * dztx[4] +
                       r[2] * dztx[5]);

    // For each point
    for (int i = 0; i < num_points; i++) {
      // For each displacement component
      for (int ii = 0; ii < 6; ii++) {
        transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
        transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

        normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

        normal_rot_bmat(dr, dr_xi, dr_eta, ii, N[i], Na[i], Nb[i], dn, dn_xi,
                        dn_eta);

        transform_vector3d_sens(dr, r, t, dt);
        transform_vector3d_sens(dr_xi, r_xi, t, dt);
        transform_vector3d_sens(dr_eta, r_eta, t, dt);

        transform_vector3d(r, t);
        transform_vector3d(r_xi, t);
        transform_vector3d(r_eta, t);

        // The mid-surface value of the displacement derivatives
        TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
        TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
        TacsScalar duz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

        TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
        TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
        TacsScalar dvz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

        TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
        TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
        TacsScalar dwz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

        // The first-derivative values of the displacements
        TacsScalar dux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                           Ud[1] * ztx[1] + r[0] * ztx[2]);
        TacsScalar duy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                           Ud[1] * ztx[4] + r[0] * ztx[5]);

        TacsScalar dvx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                           Ud[3] * ztx[1] + r[1] * ztx[2]);
        TacsScalar dvy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                           Ud[3] * ztx[4] + r[1] * ztx[5]);

        TacsScalar dwx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                           Ud[5] * ztx[1] + r[2] * ztx[2]);
        TacsScalar dwy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                           Ud[5] * ztx[4] + r[2] * ztx[5]);

        // The mid-surface values
        TacsScalar sdux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                           Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
        TacsScalar sduy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                           Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
        TacsScalar sduz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                           Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

        TacsScalar sdvx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                           Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
        TacsScalar sdvy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                           Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
        TacsScalar sdvz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                           Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

        TacsScalar sdwx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                           Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
        TacsScalar sdwy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                           Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);
        TacsScalar sdwz = (dUd[4] * tx[6] + dUd[5] * tx[7] + dr[2] * tx[8] +
                           Ud[4] * dtx[6] + Ud[5] * dtx[7] + r[2] * dtx[8]);

        // The first-derivative values of the displacements
        TacsScalar sdux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] +
                            dUd[0] * ztx[0] + dUd[1] * ztx[1] + dr[0] * ztx[2] +
                            r_xi[0] * dtx[0] + r_eta[0] * dtx[1] +
                            Ud[0] * dztx[0] + Ud[1] * dztx[1] + r[0] * dztx[2]);
        TacsScalar sduy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] +
                            dUd[0] * ztx[3] + dUd[1] * ztx[4] + dr[0] * ztx[5] +
                            r_xi[0] * dtx[3] + r_eta[0] * dtx[4] +
                            Ud[0] * dztx[3] + Ud[1] * dztx[4] + r[0] * dztx[5]);

        TacsScalar sdvx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] +
                            dUd[2] * ztx[0] + dUd[3] * ztx[1] + dr[1] * ztx[2] +
                            r_xi[1] * dtx[0] + r_eta[1] * dtx[1] +
                            Ud[2] * dztx[0] + Ud[3] * dztx[1] + r[1] * dztx[2]);
        TacsScalar sdvy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] +
                            dUd[2] * ztx[3] + dUd[3] * ztx[4] + dr[1] * ztx[5] +
                            r_xi[1] * dtx[3] + r_eta[1] * dtx[4] +
                            Ud[2] * dztx[3] + Ud[3] * dztx[4] + r[1] * dztx[5]);

        TacsScalar sdwx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] +
                            dUd[4] * ztx[0] + dUd[5] * ztx[1] + dr[2] * ztx[2] +
                            r_xi[2] * dtx[0] + r_eta[2] * dtx[1] +
                            Ud[4] * dztx[0] + Ud[5] * dztx[1] + r[2] * dztx[2]);
        TacsScalar sdwy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] +
                            dUd[4] * ztx[3] + dUd[5] * ztx[4] + dr[2] * ztx[5] +
                            r_xi[2] * dtx[3] + r_eta[2] * dtx[4] +
                            Ud[4] * dztx[3] + Ud[5] * dztx[4] + r[2] * dztx[5]);

        TacsScalar dB[8];

        // The in plane strains
        dB[0] = sdux + (sux * dux + svx * dvx + swx * dwx + ux * sdux +
                        vx * sdvx + wx * sdwx);
        dB[1] = sdvy + (suy * duy + svy * dvy + swy * dwy + uy * sduy +
                        vy * sdvy + wy * sdwy);
        dB[2] = sduy + sdvx +
                (sdux * uy + sdvx * vy + sdwx * wy + sux * duy + svx * dvy +
                 swx * dwy + dux * suy + dvx * svy + dwx * swy + ux * sduy +
                 vx * sdvy + wx * sdwy);

        // Compute the bending components of the strain
        dB[3] = sdux1 + (sdux * ux1 + sdvx * vx1 + sdwx * wx1 + sux * dux1 +
                         svx * dvx1 + swx * dwx1 + dux * sux1 + dvx * svx1 +
                         dwx * swx1 + ux * sdux1 + vx * sdvx1 + wx * sdwx1);
        dB[4] = sdvy1 + (sduy * uy1 + sdvy * vy1 + sdwy * wy1 + suy * duy1 +
                         svy * dvy1 + swy * dwy1 + duy * suy1 + dvy * svy1 +
                         dwy * swy1 + uy * sduy1 + vy * sdvy1 + wy * sdwy1);
        dB[5] =
            sduy1 + sdvx1 +
            (sdux1 * uy + sdvx1 * vy + sdwx1 * wy + sdux * uy1 + sdvx * vy1 +
             sdwx * wy1 + sux1 * duy + svx1 * dvy + swx1 * dwy + sux * duy1 +
             svx * dvy1 + swx * dwy1 + dux1 * suy + dvx1 * svy + dwx1 * swy +
             dux * suy1 + dvx * svy1 + dwx * swy1 + ux1 * sduy + vx1 * sdvy +
             wx1 * sdwy + ux * sduy1 + vx * sdvy1 + wx * sdwy1);

        // Compute the shear terms
        dB[6] = sdvz + sdwy +
                (sduz * uy + sdvz * vy + sdwz * wy + suz * duy + svz * dvy +
                 swz * dwy + duz * suy + dvz * svy + dwz * swy + uz * sduy +
                 vz * sdvy + wz * sdwy);
        dB[7] = sduz + sdwx +
                (sduz * ux + sdvz * vx + sdwz * wx + suz * dux + svz * dvx +
                 swz * dwx + duz * sux + dvz * svx + dwz * swx + uz * sdux +
                 vz * sdvx + wz * sdwx);

        // Add the contributions from the stress sensitivities
        res[0] += stress_scale *
                  (stress[0] * dB[0] + stress[1] * dB[1] + stress[2] * dB[2] +
                   stress[3] * dB[3] + stress[4] * dB[4] + stress[5] * dB[5] +
                   stress[6] * dB[6] + stress[7] * dB[7]);

        // Add the contribution from the rotational sensitivities
        res[0] += rot_scale * 0.5 * (sdvx - sduy);

        if (ii >= 3) {
          res[0] -= rot_scale * N[i] * dn[ii - 3];
        }

        res++;
      }
    }

    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
  }
}

/*
  Add the product of the stress with the second derivative of the
  nonlinear strain multiplied by stress provided.

  Note that this only fills in the bottom half of the matrix.

  input:
  num_points: number of points (nodes)
  scale: scale all the new entries
  stress: the values of the stresses

  N, Na, Nb: the shape functions
  Ux, Uxd: the displacements and the derivatives of the displacements

  t, tx, ztx: the transformation
  dt, dtx, dztx: the sensitivity of the transformation

  n, n_xi, n_eta: the normal of the shell
  dn, dn_xi, dn_eta: the derivative of the normal of the shell

  output:
  matrix: add values to the matrix
*/
void nonlinear_stress_bmat_sens(
    TacsScalar matrix[], int num_points, TacsScalar scale, TacsScalar sscale,
    const TacsScalar stress[], const TacsScalar sstress[], const double N[],
    const double Na[], const double Nb[], const TacsScalar t[],
    const TacsScalar dt[], const TacsScalar tx[], const TacsScalar dtx[],
    const TacsScalar ztx[], const TacsScalar dztx[], const TacsScalar n[],
    const TacsScalar dn[], const TacsScalar n_xi[], const TacsScalar dn_xi[],
    const TacsScalar n_eta[], const TacsScalar dn_eta[]) {
  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar Ud[12], dUd[12];
      TacsScalar r[3], r_xi[3], r_eta[3];
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
      transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);
      normal_rot_bmat(dr, dr_xi, dr_eta, ii, N[i], Na[i], Nb[i], dn, dn_xi,
                      dn_eta);

      transform_vector3d_sens(dr, r, t, dt);
      transform_vector3d_sens(dr_xi, r_xi, t, dt);
      transform_vector3d_sens(dr_eta, r_eta, t, dt);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // The mid-surface value of the displacement derivatives
      TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
      TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
      TacsScalar uz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

      TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
      TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
      TacsScalar vz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

      TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
      TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
      TacsScalar wz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

      // The first-derivative values of the displacements
      TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                        Ud[1] * ztx[1] + r[0] * ztx[2]);
      TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                        Ud[1] * ztx[4] + r[0] * ztx[5]);

      TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                        Ud[3] * ztx[1] + r[1] * ztx[2]);
      TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                        Ud[3] * ztx[4] + r[1] * ztx[5]);

      TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                        Ud[5] * ztx[1] + r[2] * ztx[2]);
      TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                        Ud[5] * ztx[4] + r[2] * ztx[5]);

      // The mid-surface values
      TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                        Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
      TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                        Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
      TacsScalar suz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                        Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

      TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                        Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
      TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                        Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
      TacsScalar svz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                        Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

      TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                        Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
      TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                        Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);
      TacsScalar swz = (dUd[4] * tx[6] + dUd[5] * tx[7] + dr[2] * tx[8] +
                        Ud[4] * dtx[6] + Ud[5] * dtx[7] + r[2] * dtx[8]);

      // The first-derivative values of the displacements
      TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] +
                         dUd[0] * ztx[0] + dUd[1] * ztx[1] + dr[0] * ztx[2] +
                         r_xi[0] * dtx[0] + r_eta[0] * dtx[1] +
                         Ud[0] * dztx[0] + Ud[1] * dztx[1] + r[0] * dztx[2]);
      TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] +
                         dUd[0] * ztx[3] + dUd[1] * ztx[4] + dr[0] * ztx[5] +
                         r_xi[0] * dtx[3] + r_eta[0] * dtx[4] +
                         Ud[0] * dztx[3] + Ud[1] * dztx[4] + r[0] * dztx[5]);

      TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] +
                         dUd[2] * ztx[0] + dUd[3] * ztx[1] + dr[1] * ztx[2] +
                         r_xi[1] * dtx[0] + r_eta[1] * dtx[1] +
                         Ud[2] * dztx[0] + Ud[3] * dztx[1] + r[1] * dztx[2]);
      TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] +
                         dUd[2] * ztx[3] + dUd[3] * ztx[4] + dr[1] * ztx[5] +
                         r_xi[1] * dtx[3] + r_eta[1] * dtx[4] +
                         Ud[2] * dztx[3] + Ud[3] * dztx[4] + r[1] * dztx[5]);

      TacsScalar swx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] +
                         dUd[4] * ztx[0] + dUd[5] * ztx[1] + dr[2] * ztx[2] +
                         r_xi[2] * dtx[0] + r_eta[2] * dtx[1] +
                         Ud[4] * dztx[0] + Ud[5] * dztx[1] + r[2] * dztx[2]);
      TacsScalar swy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] +
                         dUd[4] * ztx[3] + dUd[5] * ztx[4] + dr[2] * ztx[5] +
                         r_xi[2] * dtx[3] + r_eta[2] * dtx[4] +
                         Ud[4] * dztx[3] + Ud[5] * dztx[4] + r[2] * dztx[5]);

      TacsScalar *mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // For each point
      for (int j = 0; j <= i; j++) {
        int end = (j == i ? ii : 5);

        // For each displacement component
        for (int jj = 0; jj <= end; jj++) {
          transform_displ_gradient_bmat(Ud, jj, t, Na[j], Nb[j]);
          transform_displ_gradient_bmat(dUd, jj, dt, Na[j], Nb[j]);
          normal_rot_bmat(r, r_xi, r_eta, jj, N[j], Na[j], Nb[j], n, n_xi,
                          n_eta);

          normal_rot_bmat(dr, dr_xi, dr_eta, jj, N[j], Na[j], Nb[j], dn, dn_xi,
                          dn_eta);

          transform_vector3d_sens(dr, r, t, dt);
          transform_vector3d_sens(dr_xi, r_xi, t, dt);
          transform_vector3d_sens(dr_eta, r_eta, t, dt);

          transform_vector3d(r, t);
          transform_vector3d(r_xi, t);
          transform_vector3d(r_eta, t);

          // The mid-surface value of the displacement derivatives
          TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
          TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
          TacsScalar duz = Ud[0] * tx[6] + Ud[1] * tx[7] + r[0] * tx[8];

          TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
          TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];
          TacsScalar dvz = Ud[2] * tx[6] + Ud[3] * tx[7] + r[1] * tx[8];

          TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
          TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
          TacsScalar dwz = Ud[4] * tx[6] + Ud[5] * tx[7] + r[2] * tx[8];

          // The first-derivative values of the displacements
          TacsScalar dux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] +
                             Ud[0] * ztx[0] + Ud[1] * ztx[1] + r[0] * ztx[2]);
          TacsScalar duy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] +
                             Ud[0] * ztx[3] + Ud[1] * ztx[4] + r[0] * ztx[5]);

          TacsScalar dvx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] +
                             Ud[2] * ztx[0] + Ud[3] * ztx[1] + r[1] * ztx[2]);
          TacsScalar dvy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] +
                             Ud[2] * ztx[3] + Ud[3] * ztx[4] + r[1] * ztx[5]);

          TacsScalar dwx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] +
                             Ud[4] * ztx[0] + Ud[5] * ztx[1] + r[2] * ztx[2]);
          TacsScalar dwy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] +
                             Ud[4] * ztx[3] + Ud[5] * ztx[4] + r[2] * ztx[5]);

          // The mid-surface values
          TacsScalar sdux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                             Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
          TacsScalar sduy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                             Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
          TacsScalar sduz = (dUd[0] * tx[6] + dUd[1] * tx[7] + dr[0] * tx[8] +
                             Ud[0] * dtx[6] + Ud[1] * dtx[7] + r[0] * dtx[8]);

          TacsScalar sdvx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                             Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
          TacsScalar sdvy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                             Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);
          TacsScalar sdvz = (dUd[2] * tx[6] + dUd[3] * tx[7] + dr[1] * tx[8] +
                             Ud[2] * dtx[6] + Ud[3] * dtx[7] + r[1] * dtx[8]);

          TacsScalar sdwx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                             Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
          TacsScalar sdwy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                             Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);
          TacsScalar sdwz = (dUd[4] * tx[6] + dUd[5] * tx[7] + dr[2] * tx[8] +
                             Ud[4] * dtx[6] + Ud[5] * dtx[7] + r[2] * dtx[8]);

          // The first-derivative values of the displacements
          TacsScalar sdux1 =
              (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
               dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
               r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
               r[0] * dztx[2]);
          TacsScalar sduy1 =
              (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
               dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
               r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
               r[0] * dztx[5]);

          TacsScalar sdvx1 =
              (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
               dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
               r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
               r[1] * dztx[2]);
          TacsScalar sdvy1 =
              (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
               dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
               r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
               r[1] * dztx[5]);

          TacsScalar sdwx1 =
              (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dUd[4] * ztx[0] +
               dUd[5] * ztx[1] + dr[2] * ztx[2] + r_xi[2] * dtx[0] +
               r_eta[2] * dtx[1] + Ud[4] * dztx[0] + Ud[5] * dztx[1] +
               r[2] * dztx[2]);
          TacsScalar sdwy1 =
              (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dUd[4] * ztx[3] +
               dUd[5] * ztx[4] + dr[2] * ztx[5] + r_xi[2] * dtx[3] +
               r_eta[2] * dtx[4] + Ud[4] * dztx[3] + Ud[5] * dztx[4] +
               r[2] * dztx[5]);

          TacsScalar B[8];

          // The in plane strains
          B[0] = (ux * dux + vx * dvx + wx * dwx);
          B[1] = (uy * duy + vy * dvy + wy * dwy);
          B[2] =
              (dux * uy + dvx * vy + dwx * wy + ux * duy + vx * dvy + wx * dwy);

          // Compute the bending components of the strain
          B[3] = (dux * ux1 + dvx * vx1 + dwx * wx1 + ux * dux1 + vx * dvx1 +
                  wx * dwx1);
          B[4] = (duy * uy1 + dvy * vy1 + dwy * wy1 + uy * duy1 + vy * dvy1 +
                  wy * dwy1);
          B[5] = (dux1 * uy + dvx1 * vy + dwx1 * wy + dux * uy1 + dvx * vy1 +
                  dwx * wy1 + ux1 * duy + vx1 * dvy + wx1 * dwy + ux * duy1 +
                  vx * dvy1 + wx * dwy1);

          // Compute the shear terms
          B[6] =
              (duz * uy + dvz * vy + dwz * wy + uz * duy + vz * dvy + wz * dwy);
          B[7] =
              (duz * ux + dvz * vx + dwz * wx + uz * dux + vz * dvx + wz * dwx);

          TacsScalar bs =
              (B[0] * stress[0] + B[1] * stress[1] + B[2] * stress[2] +
               B[3] * stress[3] + B[4] * stress[4] + B[5] * stress[5] +
               B[6] * stress[6] + B[7] * stress[7]);

          TacsScalar ssb =
              (B[0] * sstress[0] + B[1] * sstress[1] + B[2] * sstress[2] +
               B[3] * sstress[3] + B[4] * sstress[4] + B[5] * sstress[5] +
               B[6] * sstress[6] + B[7] * sstress[7]);

          // Compute the sensitivity terms ...

          // The in plane strains
          B[0] = (sux * dux + svx * dvx + swx * dwx + ux * sdux + vx * sdvx +
                  wx * sdwx);
          B[1] = (suy * duy + svy * dvy + swy * dwy + uy * sduy + vy * sdvy +
                  wy * sdwy);
          B[2] = (sdux * uy + sdvx * vy + sdwx * wy + sux * duy + svx * dvy +
                  swx * dwy + dux * suy + dvx * svy + dwx * swy + ux * sduy +
                  vx * sdvy + wx * sdwy);

          // Compute the bending components of the strain
          B[3] = (sdux * ux1 + sdvx * vx1 + sdwx * wx1 + sux * dux1 +
                  svx * dvx1 + swx * dwx1 + dux * sux1 + dvx * svx1 +
                  dwx * swx1 + ux * sdux1 + vx * sdvx1 + wx * sdwx1);
          B[4] = (sduy * uy1 + sdvy * vy1 + sdwy * wy1 + suy * duy1 +
                  svy * dvy1 + swy * dwy1 + duy * suy1 + dvy * svy1 +
                  dwy * swy1 + uy * sduy1 + vy * sdvy1 + wy * sdwy1);
          B[5] =
              (sdux1 * uy + sdvx1 * vy + sdwx1 * wy + sdux * uy1 + sdvx * vy1 +
               sdwx * wy1 + sux1 * duy + svx1 * dvy + swx1 * dwy + sux * duy1 +
               svx * dvy1 + swx * dwy1 + dux1 * suy + dvx1 * svy + dwx1 * swy +
               dux * suy1 + dvx * svy1 + dwx * swy1 + ux1 * sduy + vx1 * sdvy +
               wx1 * sdwy + ux * sduy1 + vx * sdvy1 + wx * sdwy1);

          // Compute the shear terms
          B[6] = (sduz * uy + sdvz * vy + sdwz * wy + suz * duy + svz * dvy +
                  swz * dwy + duz * suy + dvz * svy + dwz * swy + uz * sduy +
                  vz * sdvy + wz * sdwy);
          B[7] = (sduz * ux + sdvz * vx + sdwz * wx + suz * dux + svz * dvx +
                  swz * dwx + duz * sux + dvz * svx + dwz * swx + uz * sdux +
                  vz * sdvx + wz * sdwx);

          TacsScalar sbs =
              (B[0] * stress[0] + B[1] * stress[1] + B[2] * stress[2] +
               B[3] * stress[3] + B[4] * stress[4] + B[5] * stress[5] +
               B[6] * stress[6] + B[7] * stress[7]);

          mat[0] += sscale * bs + scale * (ssb + sbs);

          mat += 1;
        }
      }
    }
  }
}

/*
  The following are functions for the MITC variant of the shells. For
  this implementation, only the bending terms are required. All other
  terms are interpolated from the displacement-based values evaluated
  at certain tying points. The selection of these points has a strong
  effect on the modeling capabilities of the element.
*/

/*
  Compute the linear strain without the in plane components
*/
void linear_bend_strain(TacsScalar strain[], TacsScalar *rotz,
                        const TacsScalar Ux[], const TacsScalar Uxd[],
                        const TacsScalar t[], const TacsScalar tx[],
                        const TacsScalar ztx[], const TacsScalar n[],
                        const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // rotz = 0.5 * (v,x - u,y ) - inplane_rot
  *rotz = (0.5 * (Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2] -
                  (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5])) -
           (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  strain[0] = strain[1] = strain[2] = 0.0;
  strain[6] = strain[7] = 0.0;

  // Compute the bending components of the strain
  strain[3] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
               Ud[1] * ztx[1] + r[0] * ztx[2]);  // k_11
  strain[4] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
               Ud[3] * ztx[4] + r[1] * ztx[5]);  // k_22

  strain[5] =
      (r_xi[0] * tx[3] + r_eta[0] * tx[4] + r_xi[1] * tx[0] + r_eta[1] * tx[1] +
       Ud[0] * ztx[3] + Ud[1] * ztx[4] + r[0] * ztx[5] + Ud[2] * ztx[0] +
       Ud[3] * ztx[1] + r[1] * ztx[2]);  // k_12
}

void linear_bend_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                      const double N[], const double Na[], const double Nb[],
                      const TacsScalar t[], const TacsScalar tx[],
                      const TacsScalar ztx[], const TacsScalar n[],
                      const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  // For each point
  for (int i = 0; i < num_points; i++) {
    // For each displacement component
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar Ud[12];
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      TacsScalar r[3], r_xi[3], r_eta[3];
      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // rotz = 0.5 * (v,x - u,y )
      rotz[0] = 0.5 * (Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2] -
                       (Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5]));
      if (ii >= 3) {
        rotz[0] -= N[i] * n[ii - 3];
      }

      // Compute the bending components of the strain
      B[3] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
              Ud[1] * ztx[1] + r[0] * ztx[2]);  // k_11
      B[4] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
              Ud[3] * ztx[4] + r[1] * ztx[5]);  // k_22

      B[5] =
          (r_xi[0] * tx[3] + r_eta[0] * tx[4] + r_xi[1] * tx[0] +
           r_eta[1] * tx[1] + Ud[0] * ztx[3] + Ud[1] * ztx[4] + r[0] * ztx[5] +
           Ud[2] * ztx[0] + Ud[3] * ztx[1] + r[1] * ztx[2]);  // k_12

      B += 8;
      rotz++;
    }
  }
}

/*
  Compute the design variable sensitivity of the strain to a perturbation
  of the position variables.

  input:
  Ux, Uxd: the displacements and the derivatives of the displacements in
  the global reference frame

  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivities of the transformations

  n, n_xi, n_eta: the normal and derivatives of the normals
  dn, dn_xi, dn_eta: the sensitivities of the normals
*/

void linear_bend_strain_sens(TacsScalar strain[], TacsScalar dstrain[],
                             TacsScalar *rotz, TacsScalar *drotz,
                             const TacsScalar Ux[], const TacsScalar Uxd[],
                             const TacsScalar t[], const TacsScalar dt[],
                             const TacsScalar tx[], const TacsScalar dtx[],
                             const TacsScalar ztx[], const TacsScalar dztx[],
                             const TacsScalar n[], const TacsScalar dn[],
                             const TacsScalar n_xi[], const TacsScalar dn_xi[],
                             const TacsScalar n_eta[],
                             const TacsScalar dn_eta[], int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface values
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];
  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  // rotz = 0.5 * (v,x - u,y ) - inplane_rot
  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // Compute the bending components of the strain
  strain[3] = ux1;
  strain[4] = vy1;
  strain[5] = uy1 + vx1;

  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    drotz[0] =
        (0.5 * (svx - suy) - (dn[0] * Ux[3] + dn[1] * Ux[4] + dn[2] * Ux[5]));

    // Compute the bending components of the strain
    dstrain[3] = sux1;
    dstrain[4] = svy1;
    dstrain[5] = suy1 + svx1;

    dstrain += 8;
    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
    drotz++;
  }
}

/*
  Add the components of the sensitivity to the matrix res.

  input:
  num_points: number of points (number of nodes)
  stress_scale: the scale value for the stresses
  rot_scale: the value to scale the rotation derivative
  stress: the components of the stress

  N, Na, Nb: the shape functions and their derivatives
  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivity of the transformations

  n, n_xi, n_eta: the normal and the derivative of the normal along
  the xi and eta directions
  dn, dn_xi, dn_eta: the sensitivity of the normal (and its derivatives)
*/
void add_linear_bend_bmat_sens(
    TacsScalar fXptSens[], const TacsScalar psi[], int num_points,
    const TacsScalar stress_scale, const TacsScalar rot_scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar t[], const TacsScalar dt[],
    const TacsScalar tx[], const TacsScalar dtx[], const TacsScalar ztx[],
    const TacsScalar dztx[], const TacsScalar n[], const TacsScalar dn[],
    const TacsScalar n_xi[], const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components) {
  // For each component in the derivative
  for (int k = 0; k < num_components; k++) {
    // For each point
    for (int i = 0; i < num_points; i++) {
      // For each displacement component
      for (int ii = 0; ii < 6; ii++) {
        TacsScalar Ud[12], dUd[12];
        transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
        transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

        TacsScalar r[3], r_xi[3], r_eta[3];
        TacsScalar dr[3], dr_xi[3], dr_eta[3];
        normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

        normal_rot_bmat(dr, dr_xi, dr_eta, ii, N[i], Na[i], Nb[i], dn, dn_xi,
                        dn_eta);

        transform_vector3d_sens(dr, r, t, dt);
        transform_vector3d_sens(dr_xi, r_xi, t, dt);
        transform_vector3d_sens(dr_eta, r_eta, t, dt);

        transform_vector3d(r, t);
        transform_vector3d(r_xi, t);
        transform_vector3d(r_eta, t);

        // The mid-surface values
        TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                          Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);
        TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                          Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);

        // The first-derivative values of the displacements
        TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] +
                           dUd[0] * ztx[0] + dUd[1] * ztx[1] + dr[0] * ztx[2] +
                           r_xi[0] * dtx[0] + r_eta[0] * dtx[1] +
                           Ud[0] * dztx[0] + Ud[1] * dztx[1] + r[0] * dztx[2]);
        TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] +
                           dUd[0] * ztx[3] + dUd[1] * ztx[4] + dr[0] * ztx[5] +
                           r_xi[0] * dtx[3] + r_eta[0] * dtx[4] +
                           Ud[0] * dztx[3] + Ud[1] * dztx[4] + r[0] * dztx[5]);

        TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] +
                           dUd[2] * ztx[0] + dUd[3] * ztx[1] + dr[1] * ztx[2] +
                           r_xi[1] * dtx[0] + r_eta[1] * dtx[1] +
                           Ud[2] * dztx[0] + Ud[3] * dztx[1] + r[1] * dztx[2]);
        TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] +
                           dUd[2] * ztx[3] + dUd[3] * ztx[4] + dr[1] * ztx[5] +
                           r_xi[1] * dtx[3] + r_eta[1] * dtx[4] +
                           Ud[2] * dztx[3] + Ud[3] * dztx[4] + r[1] * dztx[5]);

        // Add the contributions from the stress sensitivities
        fXptSens[k] +=
            stress_scale * psi[6 * i + ii] *
            (stress[3] * sux1 + stress[4] * svy1 + stress[5] * (suy1 + svx1));

        // Add the contribution from the rotational sensitivities
        fXptSens[k] += 0.5 * rot_scale * psi[6 * i + ii] * (svx - suy);

        if (ii >= 3) {
          fXptSens[k] -= rot_scale * psi[6 * i + ii] * N[i] * dn[ii - 3];
        }
      }
    }

    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
  }
}

/*
  Compute the nonlinear components of the strain in the local coordinate
  system.

  input:
  Ux, Uxd: the displacements in the global coordinate system
  N, Na, Nb: the shape functions and their derivatives
  t, tx, ztx: the transformations
  n, n_xi, n_eta: the normal and the derivative of the normal along the
  xi and eta directions

  output:
  strain: the strain values
  rotz: the difference between the drilling rotation and the normal
  rotation
*/
void nonlinear_bend_strain(TacsScalar strain[], TacsScalar *rotz,
                           const TacsScalar Ux[], const TacsScalar Uxd[],
                           const TacsScalar t[], const TacsScalar tx[],
                           const TacsScalar ztx[], const TacsScalar n[],
                           const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);
}

/*
  Add the nonlinear components corresponding to the bending strain

  input:
  N, Na, Nb: the shape functions and their derivatives
  Ux, Uxd: the displacements in the global coordinate system and their
  derivatives
  t, tx, ztx: the transformations

  output:
  B: the derivative of the strain w.r.t. the nodal displacements and
  rotations
  rotz: the difference between the drilling and normal rotations
*/
void nonlinear_bend_bmat(TacsScalar B[], TacsScalar rotz[], int num_points,
                         const double N[], const double Na[], const double Nb[],
                         const TacsScalar Ux[], const TacsScalar Uxd[],
                         const TacsScalar t[], const TacsScalar tx[],
                         const TacsScalar ztx[], const TacsScalar n[],
                         const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // For each point
  for (int i = 0; i < num_points; i++) {
    // For each displacement component
    for (int ii = 0; ii < 6; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
      TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

      TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
      TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

      TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
      TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

      // The first-derivative values of the displacements
      TacsScalar dux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                         Ud[1] * ztx[1] + r[0] * ztx[2]);
      TacsScalar duy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                         Ud[1] * ztx[4] + r[0] * ztx[5]);

      TacsScalar dvx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                         Ud[3] * ztx[1] + r[1] * ztx[2]);
      TacsScalar dvy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                         Ud[3] * ztx[4] + r[1] * ztx[5]);

      TacsScalar dwx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                         Ud[5] * ztx[1] + r[2] * ztx[2]);
      TacsScalar dwy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                         Ud[5] * ztx[4] + r[2] * ztx[5]);

      rotz[0] = 0.5 * (dvx - duy);
      if (ii >= 3) {
        rotz[0] -= N[i] * n[ii - 3];
      }

      // Compute the bending components of the strain
      B[3] = dux1 + (dux * ux1 + dvx * vx1 + dwx * wx1 + ux * dux1 + vx * dvx1 +
                     wx * dwx1);
      B[4] = dvy1 + (duy * uy1 + dvy * vy1 + dwy * wy1 + uy * duy1 + vy * dvy1 +
                     wy * dwy1);
      B[5] = duy1 + dvx1 +
             (dux1 * uy + dvx1 * vy + dwx1 * wy + dux * uy1 + dvx * vy1 +
              dwx * wy1 + ux1 * duy + vx1 * dvy + wx1 * dwy + ux * duy1 +
              vx * dvy1 + wx * dwy1);

      B += 8;
      rotz++;
    }
  }
}

/*
  Add the product of the stress with the second derivative of the
  nonlinear strain multiplied by stress provided.

  Note that this only fills in the bottom half of the matrix.

  input:
  scale: scale all new entries into the matrix by this value
  stress: the stress at the integration point
  N, Na, Nb: the shape functions and their derivatives
  Ux, Uxd: the displacements in the global coordinate system and their
  derivatives
  t, tx, ztx: the transformations

  output:
  matrix: the matrix values
  num_points: the number of nodes
*/
void nonlinear_bend_stress_bmat(TacsScalar matrix[], int num_points,
                                TacsScalar scale, const TacsScalar stress[],
                                const double N[], const double Na[],
                                const double Nb[], const TacsScalar t[],
                                const TacsScalar tx[], const TacsScalar ztx[],
                                const TacsScalar n[], const TacsScalar n_xi[],
                                const TacsScalar n_eta[]) {
  if (num_points > MAX_NUM_NODES) {
    return;
  }

  TacsScalar dU0[6 * 6 * MAX_NUM_NODES];
  TacsScalar dU1[6 * 6 * MAX_NUM_NODES];

  TacsScalar *du0 = dU0;
  TacsScalar *du1 = dU1;

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      // Compute the derivative of the deformation gradient w.r.t.
      // the ii-th displacement at the i-th node
      TacsScalar Ud[12];
      TacsScalar r[3], r_xi[3], r_eta[3];
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
      normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

      transform_vector3d(r, t);
      transform_vector3d(r_xi, t);
      transform_vector3d(r_eta, t);

      // The mid-surface value of the displacement derivatives
      // ux, uy
      du0[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
      du0[1] = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

      // vx, vy
      du0[2] = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
      du0[3] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

      // wx, wy
      du0[4] = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
      du0[5] = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];
      du0 += 6;

      // The first-derivative values of the displacements
      // u1x, u1y
      du1[0] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                Ud[1] * ztx[1] + r[0] * ztx[2]);
      du1[1] = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                Ud[1] * ztx[4] + r[0] * ztx[5]);

      // v1x, v1y
      du1[2] = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                Ud[3] * ztx[1] + r[1] * ztx[2]);
      du1[3] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                Ud[3] * ztx[4] + r[1] * ztx[5]);

      // w1x, w1y
      du1[4] = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                Ud[5] * ztx[1] + r[2] * ztx[2]);
      du1[5] = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                Ud[5] * ztx[4] + r[2] * ztx[5]);
      du1 += 6;
    }
  }

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar *mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Set the i-counter from the beginning of the row
      TacsScalar *di0 = &dU0[6 * (6 * i + ii)];
      TacsScalar *di1 = &dU1[6 * (6 * i + ii)];

      // Start the j-counter from the i-counter
      TacsScalar *dj0 = dU0;
      TacsScalar *dj1 = dU1;

      // For each point
      for (int j = 0; j <= i; j++) {
        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        // For each displacement component
        for (int jj = 0; jj < end; jj++) {
          TacsScalar B[3];

          // Compute the bending components of the strain
          B[0] = (di0[0] * dj1[0] + di0[2] * dj1[2] + di0[4] * dj1[4] +
                  di1[0] * dj0[0] + di1[2] * dj0[2] + di1[4] * dj0[4]);
          B[1] = (di0[1] * dj1[1] + di0[3] * dj1[3] + di0[5] * dj1[5] +
                  di1[1] * dj0[1] + di1[3] * dj0[3] + di1[5] * dj0[5]);
          B[2] = (di0[0] * dj1[1] + di0[2] * dj1[3] + di0[4] * dj1[5] +
                  di1[0] * dj0[1] + di1[2] * dj0[3] + di1[4] * dj0[5] +
                  dj0[0] * di1[1] + dj0[2] * di1[3] + dj0[4] * di1[5] +
                  dj1[0] * di0[1] + dj1[2] * di0[3] + dj1[4] * di0[5]);

          mat[0] +=
              scale * (B[0] * stress[3] + B[1] * stress[4] + B[2] * stress[5]);
          mat += 1;

          dj0 += 6;
          dj1 += 6;
        }
      }
    }
  }
}

/*
  Compute the inner product of the second derivative of the nonlinear
  bending components of the strain with respect to two vectors.

  input:
  Uxphi, Uxdphi:  the displacements/derivatives in the global coordinate system
  Uxpsi, Uxdpsi:  the displacements/derivatives in the global coordinate system
  t, tx, ztx:     the transformations
  n, n_xi, n_eta: the normal and its derivatives

  output:
  bstrain:  the product with the second derivative of the strain
*/
void inner_nonlinear_bend_bmat(TacsScalar bstrain[], const TacsScalar Uxpsi[],
                               const TacsScalar Uxdpsi[],
                               const TacsScalar Uxphi[],
                               const TacsScalar Uxdphi[], const TacsScalar t[],
                               const TacsScalar tx[], const TacsScalar ztx[],
                               const TacsScalar n[], const TacsScalar n_xi[],
                               const TacsScalar n_eta[]) {
  // The derivatives of the in-plane and rotations of the phi/psi
  // vectors with respect to the local x/y coordinates
  TacsScalar du0psi[6], du1psi[6];
  TacsScalar du0phi[6], du1phi[6];

  // Temporary storage used in this computation
  TacsScalar Ud[12];
  TacsScalar r[3], r_xi[3], r_eta[3];

  // Compute the displacement gradient
  transform_displ_gradient(Ud, t, Uxdpsi);

  // Compute the rotation of the normal vector and its derivative
  // along xi/eta directions
  normal_rot(r, r_xi, r_eta, Uxpsi, Uxdpsi, n, n_xi, n_eta);

  // Transform the rate vectors through the thickness into
  // the local reference frame
  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  // ux, uy
  du0psi[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  du0psi[1] = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  // vx, vy
  du0psi[2] = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  du0psi[3] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  // wx, wy
  du0psi[4] = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  du0psi[5] = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  // u1x, u1y
  du1psi[0] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
               Ud[1] * ztx[1] + r[0] * ztx[2]);
  du1psi[1] = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
               Ud[1] * ztx[4] + r[0] * ztx[5]);

  // v1x, v1y
  du1psi[2] = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
               Ud[3] * ztx[1] + r[1] * ztx[2]);
  du1psi[3] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
               Ud[3] * ztx[4] + r[1] * ztx[5]);

  // w1x, w1y
  du1psi[4] = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
               Ud[5] * ztx[1] + r[2] * ztx[2]);
  du1psi[5] = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
               Ud[5] * ztx[4] + r[2] * ztx[5]);

  // Compute the displacement gradient
  transform_displ_gradient(Ud, t, Uxdphi);

  // Compute the rotation of the normal vector and its derivative
  // along xi/eta directions
  normal_rot(r, r_xi, r_eta, Uxphi, Uxdphi, n, n_xi, n_eta);

  // Transform the rate vectors through the thickness into
  // the local reference frame
  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  // ux, uy
  du0phi[0] = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  du0phi[1] = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  // vx, vy
  du0phi[2] = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  du0phi[3] = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  // wx, wy
  du0phi[4] = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  du0phi[5] = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  // u1x, u1y
  du1phi[0] = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
               Ud[1] * ztx[1] + r[0] * ztx[2]);
  du1phi[1] = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
               Ud[1] * ztx[4] + r[0] * ztx[5]);

  // v1x, v1y
  du1phi[2] = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
               Ud[3] * ztx[1] + r[1] * ztx[2]);
  du1phi[3] = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
               Ud[3] * ztx[4] + r[1] * ztx[5]);

  // w1x, w1y
  du1phi[4] = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
               Ud[5] * ztx[1] + r[2] * ztx[2]);
  du1phi[5] = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
               Ud[5] * ztx[4] + r[2] * ztx[5]);

  // Compute the bending components of the strain
  bstrain[3] =
      (du0phi[0] * du1psi[0] + du0phi[2] * du1psi[2] + du0phi[4] * du1psi[4] +
       du1phi[0] * du0psi[0] + du1phi[2] * du0psi[2] + du1phi[4] * du0psi[4]);
  bstrain[4] =
      (du0phi[1] * du1psi[1] + du0phi[3] * du1psi[3] + du0phi[5] * du1psi[5] +
       du1phi[1] * du0psi[1] + du1phi[3] * du0psi[3] + du1phi[5] * du0psi[5]);
  bstrain[5] =
      (du0phi[0] * du1psi[1] + du0phi[2] * du1psi[3] + du0phi[4] * du1psi[5] +
       du1phi[0] * du0psi[1] + du1phi[2] * du0psi[3] + du1phi[4] * du0psi[5] +
       du0psi[0] * du1phi[1] + du0psi[2] * du1phi[3] + du0psi[4] * du1phi[5] +
       du1psi[0] * du0phi[1] + du1psi[2] * du0phi[3] + du1psi[4] * du0phi[5]);
}

/*
  Compute the sensitivity of the bending strain

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  t: the transformation
  dt: the sensitivity of the local transformation

  tx: the transformation times the Jacobian
  dtx: the sensitivity of the transformation times the Jacobian

  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction
  dztx: the sensitivity of ztx

  n: the shell normal
  dn: the sensitivity of the shell normal

  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  dn_xi: the sensitivity of n_xi

  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta
  dn_eta: the sensitivity of n_eta

  num_components: the number of sensitivity components to compute

  output:
  strain: the value of the strain at the current point
  dstrain: the sentivity of the strain at the curren point

  rotz: the difference between the drilling rotation and actual rotation
  drotz: the sensitivity of rotz
*/
void nonlinear_bend_strain_sens(
    TacsScalar strain[], TacsScalar dstrain[], TacsScalar *rotz,
    TacsScalar *drotz, const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface values
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // rotz = 0.5 * (v,x - u,y ) - inplane_rot
  *rotz = (0.5 * (vx - uy) - (n[0] * Ux[3] + n[1] * Ux[4] + n[2] * Ux[5]));

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);

  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                      Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);

    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
    TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                      Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);

    TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                      Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
    TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                      Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    TacsScalar swx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dUd[4] * ztx[0] +
                       dUd[5] * ztx[1] + dr[2] * ztx[2] + r_xi[2] * dtx[0] +
                       r_eta[2] * dtx[1] + Ud[4] * dztx[0] + Ud[5] * dztx[1] +
                       r[2] * dztx[2]);
    TacsScalar swy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dUd[4] * ztx[3] +
                       dUd[5] * ztx[4] + dr[2] * ztx[5] + r_xi[2] * dtx[3] +
                       r_eta[2] * dtx[4] + Ud[4] * dztx[3] + Ud[5] * dztx[4] +
                       r[2] * dztx[5]);

    drotz[0] =
        (0.5 * (svx - suy) - (dn[0] * Ux[3] + dn[1] * Ux[4] + dn[2] * Ux[5]));

    // Compute the bending components of the strain
    dstrain[3] = sux1 + (sux * ux1 + svx * vx1 + swx * wx1 + ux * sux1 +
                         vx * svx1 + wx * swx1);
    dstrain[4] = svy1 + (suy * uy1 + svy * vy1 + swy * wy1 + uy * suy1 +
                         vy * svy1 + wy * swy1);
    dstrain[5] =
        suy1 + svx1 +
        (sux1 * uy + svx1 * vy + swx1 * wy + sux * uy1 + svx * vy1 + swx * wy1 +
         ux1 * suy + vx1 * svy + wx1 * swy + ux * suy1 + vx * svy1 + wx * swy1);

    dstrain += 8;
    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
    drotz++;
  }
}

/*
  Add components of the sensitivity w.r.t. the nodes of the B matrix
  (the derivative of the strain w.r.t. the nodal coordinates) - to the
  residual matrix 'res'. This is used for constructing the element
  residuals.

  input:
  num_points: number of points (number of nodes)
  stress_scale: the scale value for the stresses
  rot_scale: the value to scale the rotation derivative
  stress: the components of the stress

  N, Na, Nb: the shape functions and their derivatives
  t, tx, ztx: the transformations
  dt, dtx, dztx: the sensitivity of the transformations

  n, n_xi, n_eta: the normal and the derivative of the normal along
  the xi and eta directions
  dn, dn_xi, dn_eta: the sensitivity of the normal (and its derivatives)
*/
void add_nonlinear_bend_bmat_sens(
    TacsScalar fXptSens[], const TacsScalar psi[], int num_points,
    const TacsScalar stress_scale, const TacsScalar rot_scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], int num_components) {
  TacsScalar Ud[12];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rotational contributions
  TacsScalar r[3], r_xi[3], r_eta[3];
  normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

  transform_vector3d(r, t);
  transform_vector3d(r_xi, t);
  transform_vector3d(r_eta, t);

  // The mid-surface value of the displacement derivatives
  TacsScalar ux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
  TacsScalar uy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

  TacsScalar vx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
  TacsScalar vy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

  TacsScalar wx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
  TacsScalar wy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

  // The first-derivative values of the displacements
  TacsScalar ux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                    Ud[1] * ztx[1] + r[0] * ztx[2]);
  TacsScalar uy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                    Ud[1] * ztx[4] + r[0] * ztx[5]);

  TacsScalar vx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                    Ud[3] * ztx[1] + r[1] * ztx[2]);
  TacsScalar vy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                    Ud[3] * ztx[4] + r[1] * ztx[5]);

  TacsScalar wx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                    Ud[5] * ztx[1] + r[2] * ztx[2]);
  TacsScalar wy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                    Ud[5] * ztx[4] + r[2] * ztx[5]);

  // For each component in the derivative
  for (int k = 0; k < num_components; k++) {
    TacsScalar dUd[12];
    transform_displ_gradient(Ud, t, Uxd);
    transform_displ_gradient(dUd, dt, Uxd);

    TacsScalar dr[3], dr_xi[3], dr_eta[3];
    normal_rot(dr, dr_xi, dr_eta, Ux, Uxd, dn, dn_xi, dn_eta);
    normal_rot(r, r_xi, r_eta, Ux, Uxd, n, n_xi, n_eta);

    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface values
    TacsScalar sux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                      Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
    TacsScalar suy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                      Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);

    TacsScalar svx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                      Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
    TacsScalar svy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                      Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);

    TacsScalar swx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                      Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
    TacsScalar swy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                      Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);

    // The first-derivative values of the displacements
    TacsScalar sux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dUd[0] * ztx[0] +
                       dUd[1] * ztx[1] + dr[0] * ztx[2] + r_xi[0] * dtx[0] +
                       r_eta[0] * dtx[1] + Ud[0] * dztx[0] + Ud[1] * dztx[1] +
                       r[0] * dztx[2]);
    TacsScalar suy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dUd[0] * ztx[3] +
                       dUd[1] * ztx[4] + dr[0] * ztx[5] + r_xi[0] * dtx[3] +
                       r_eta[0] * dtx[4] + Ud[0] * dztx[3] + Ud[1] * dztx[4] +
                       r[0] * dztx[5]);

    TacsScalar svx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dUd[2] * ztx[0] +
                       dUd[3] * ztx[1] + dr[1] * ztx[2] + r_xi[1] * dtx[0] +
                       r_eta[1] * dtx[1] + Ud[2] * dztx[0] + Ud[3] * dztx[1] +
                       r[1] * dztx[2]);
    TacsScalar svy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dUd[2] * ztx[3] +
                       dUd[3] * ztx[4] + dr[1] * ztx[5] + r_xi[1] * dtx[3] +
                       r_eta[1] * dtx[4] + Ud[2] * dztx[3] + Ud[3] * dztx[4] +
                       r[1] * dztx[5]);

    TacsScalar swx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dUd[4] * ztx[0] +
                       dUd[5] * ztx[1] + dr[2] * ztx[2] + r_xi[2] * dtx[0] +
                       r_eta[2] * dtx[1] + Ud[4] * dztx[0] + Ud[5] * dztx[1] +
                       r[2] * dztx[2]);
    TacsScalar swy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dUd[4] * ztx[3] +
                       dUd[5] * ztx[4] + dr[2] * ztx[5] + r_xi[2] * dtx[3] +
                       r_eta[2] * dtx[4] + Ud[4] * dztx[3] + Ud[5] * dztx[4] +
                       r[2] * dztx[5]);
    // For each point
    for (int i = 0; i < num_points; i++) {
      // For each displacement component
      for (int ii = 0; ii < 6; ii++) {
        transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
        transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

        normal_rot_bmat(r, r_xi, r_eta, ii, N[i], Na[i], Nb[i], n, n_xi, n_eta);

        normal_rot_bmat(dr, dr_xi, dr_eta, ii, N[i], Na[i], Nb[i], dn, dn_xi,
                        dn_eta);

        transform_vector3d_sens(dr, r, t, dt);
        transform_vector3d_sens(dr_xi, r_xi, t, dt);
        transform_vector3d_sens(dr_eta, r_eta, t, dt);

        transform_vector3d(r, t);
        transform_vector3d(r_xi, t);
        transform_vector3d(r_eta, t);

        // The mid-surface value of the displacement derivatives
        TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1] + r[0] * tx[2];
        TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4] + r[0] * tx[5];

        TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1] + r[1] * tx[2];
        TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4] + r[1] * tx[5];

        TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1] + r[2] * tx[2];
        TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4] + r[2] * tx[5];

        // The first-derivative values of the displacements
        TacsScalar dux1 = (r_xi[0] * tx[0] + r_eta[0] * tx[1] + Ud[0] * ztx[0] +
                           Ud[1] * ztx[1] + r[0] * ztx[2]);
        TacsScalar duy1 = (r_xi[0] * tx[3] + r_eta[0] * tx[4] + Ud[0] * ztx[3] +
                           Ud[1] * ztx[4] + r[0] * ztx[5]);

        TacsScalar dvx1 = (r_xi[1] * tx[0] + r_eta[1] * tx[1] + Ud[2] * ztx[0] +
                           Ud[3] * ztx[1] + r[1] * ztx[2]);
        TacsScalar dvy1 = (r_xi[1] * tx[3] + r_eta[1] * tx[4] + Ud[2] * ztx[3] +
                           Ud[3] * ztx[4] + r[1] * ztx[5]);

        TacsScalar dwx1 = (r_xi[2] * tx[0] + r_eta[2] * tx[1] + Ud[4] * ztx[0] +
                           Ud[5] * ztx[1] + r[2] * ztx[2]);
        TacsScalar dwy1 = (r_xi[2] * tx[3] + r_eta[2] * tx[4] + Ud[4] * ztx[3] +
                           Ud[5] * ztx[4] + r[2] * ztx[5]);

        // The mid-surface values
        TacsScalar sdux = (dUd[0] * tx[0] + dUd[1] * tx[1] + dr[0] * tx[2] +
                           Ud[0] * dtx[0] + Ud[1] * dtx[1] + r[0] * dtx[2]);
        TacsScalar sduy = (dUd[0] * tx[3] + dUd[1] * tx[4] + dr[0] * tx[5] +
                           Ud[0] * dtx[3] + Ud[1] * dtx[4] + r[0] * dtx[5]);

        TacsScalar sdvx = (dUd[2] * tx[0] + dUd[3] * tx[1] + dr[1] * tx[2] +
                           Ud[2] * dtx[0] + Ud[3] * dtx[1] + r[1] * dtx[2]);
        TacsScalar sdvy = (dUd[2] * tx[3] + dUd[3] * tx[4] + dr[1] * tx[5] +
                           Ud[2] * dtx[3] + Ud[3] * dtx[4] + r[1] * dtx[5]);

        TacsScalar sdwx = (dUd[4] * tx[0] + dUd[5] * tx[1] + dr[2] * tx[2] +
                           Ud[4] * dtx[0] + Ud[5] * dtx[1] + r[2] * dtx[2]);
        TacsScalar sdwy = (dUd[4] * tx[3] + dUd[5] * tx[4] + dr[2] * tx[5] +
                           Ud[4] * dtx[3] + Ud[5] * dtx[4] + r[2] * dtx[5]);

        // The first-derivative values of the displacements
        TacsScalar sdux1 = (dr_xi[0] * tx[0] + dr_eta[0] * tx[1] +
                            dUd[0] * ztx[0] + dUd[1] * ztx[1] + dr[0] * ztx[2] +
                            r_xi[0] * dtx[0] + r_eta[0] * dtx[1] +
                            Ud[0] * dztx[0] + Ud[1] * dztx[1] + r[0] * dztx[2]);
        TacsScalar sduy1 = (dr_xi[0] * tx[3] + dr_eta[0] * tx[4] +
                            dUd[0] * ztx[3] + dUd[1] * ztx[4] + dr[0] * ztx[5] +
                            r_xi[0] * dtx[3] + r_eta[0] * dtx[4] +
                            Ud[0] * dztx[3] + Ud[1] * dztx[4] + r[0] * dztx[5]);

        TacsScalar sdvx1 = (dr_xi[1] * tx[0] + dr_eta[1] * tx[1] +
                            dUd[2] * ztx[0] + dUd[3] * ztx[1] + dr[1] * ztx[2] +
                            r_xi[1] * dtx[0] + r_eta[1] * dtx[1] +
                            Ud[2] * dztx[0] + Ud[3] * dztx[1] + r[1] * dztx[2]);
        TacsScalar sdvy1 = (dr_xi[1] * tx[3] + dr_eta[1] * tx[4] +
                            dUd[2] * ztx[3] + dUd[3] * ztx[4] + dr[1] * ztx[5] +
                            r_xi[1] * dtx[3] + r_eta[1] * dtx[4] +
                            Ud[2] * dztx[3] + Ud[3] * dztx[4] + r[1] * dztx[5]);

        TacsScalar sdwx1 = (dr_xi[2] * tx[0] + dr_eta[2] * tx[1] +
                            dUd[4] * ztx[0] + dUd[5] * ztx[1] + dr[2] * ztx[2] +
                            r_xi[2] * dtx[0] + r_eta[2] * dtx[1] +
                            Ud[4] * dztx[0] + Ud[5] * dztx[1] + r[2] * dztx[2]);
        TacsScalar sdwy1 = (dr_xi[2] * tx[3] + dr_eta[2] * tx[4] +
                            dUd[4] * ztx[3] + dUd[5] * ztx[4] + dr[2] * ztx[5] +
                            r_xi[2] * dtx[3] + r_eta[2] * dtx[4] +
                            Ud[4] * dztx[3] + Ud[5] * dztx[4] + r[2] * dztx[5]);

        TacsScalar dB[3];

        // Compute the bending components of the strain
        dB[0] = sdux1 + (sdux * ux1 + sdvx * vx1 + sdwx * wx1 + sux * dux1 +
                         svx * dvx1 + swx * dwx1 + dux * sux1 + dvx * svx1 +
                         dwx * swx1 + ux * sdux1 + vx * sdvx1 + wx * sdwx1);
        dB[1] = sdvy1 + (sduy * uy1 + sdvy * vy1 + sdwy * wy1 + suy * duy1 +
                         svy * dvy1 + swy * dwy1 + duy * suy1 + dvy * svy1 +
                         dwy * swy1 + uy * sduy1 + vy * sdvy1 + wy * sdwy1);
        dB[2] =
            sduy1 + sdvx1 +
            (sdux1 * uy + sdvx1 * vy + sdwx1 * wy + sdux * uy1 + sdvx * vy1 +
             sdwx * wy1 + sux1 * duy + svx1 * dvy + swx1 * dwy + sux * duy1 +
             svx * dvy1 + swx * dwy1 + dux1 * suy + dvx1 * svy + dwx1 * swy +
             dux * suy1 + dvx * svy1 + dwx * swy1 + ux1 * sduy + vx1 * sdvy +
             wx1 * sdwy + ux * sduy1 + vx * sdvy1 + wx * sdwy1);

        // Add the contributions from the stress sensitivities
        fXptSens[k] +=
            stress_scale * psi[6 * i + ii] *
            (stress[3] * dB[0] + stress[4] * dB[1] + stress[5] * dB[2]);

        // Add the contribution from the rotational sensitivities
        fXptSens[k] += 0.5 * rot_scale * psi[6 * i + ii] * (sdvx - sduy);

        if (ii >= 3) {
          fXptSens[k] -= rot_scale * psi[6 * i + ii] * N[i] * dn[ii - 3];
        }
      }
    }

    dt += 9, dtx += 9, dztx += 9;
    dn += 3, dn_xi += 3, dn_eta += 3;
  }
}

/*
  This is a helper function for the large rotation strain utilities

  Given the order and the tensor-product contributions to the strain
  moments, compute the displacements, rotations and rate of change of
  displacements.
*/
void compute_tensorial_components(const int order, TacsScalar Urot[],
                                  TacsScalar Ud[], TacsScalar Xd[], double N[],
                                  double Na[], double Nb[], const double pt[],
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[]) {
  if (order == 2) {
    FElibrary::biLagrangeSF(N, Na, Nb, pt, 2);

    // Compute the rotations
    Urot[0] =
        vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3];
    Urot[1] =
        vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3];
    Urot[2] =
        vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3];

    // Compute the derivatives of the displacement
    Ud[0] =
        vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] + vars[18] * Na[3];
    Ud[2] =
        vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] + vars[19] * Na[3];
    Ud[4] =
        vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] + vars[20] * Na[3];

    // Compute the derivative of X along a
    Xd[0] =
        Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] + Xpts[9] * Na[3];
    Xd[1] =
        Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] + Xpts[10] * Na[3];
    Xd[2] =
        Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] + Xpts[11] * Na[3];

    // Compute the derivative of the displacement
    Ud[1] =
        vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] + vars[18] * Nb[3];
    Ud[3] =
        vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] + vars[19] * Nb[3];
    Ud[5] =
        vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] + vars[20] * Nb[3];

    // Compute the derivative of X along b
    Xd[3] =
        Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] + Xpts[9] * Nb[3];
    Xd[4] =
        Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] + Xpts[10] * Nb[3];
    Xd[5] =
        Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] + Xpts[11] * Nb[3];
  } else if (order == 3) {
    FElibrary::biLagrangeSF(N, Na, Nb, pt, 3);

    // Compute the rotations
    Urot[0] = (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] +
               vars[21] * N[3] + vars[27] * N[4] + vars[33] * N[5] +
               vars[39] * N[6] + vars[45] * N[7] + vars[51] * N[8]);
    Urot[1] = (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] +
               vars[22] * N[3] + vars[28] * N[4] + vars[34] * N[5] +
               vars[40] * N[6] + vars[46] * N[7] + vars[52] * N[8]);
    Urot[2] = (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] +
               vars[23] * N[3] + vars[29] * N[4] + vars[35] * N[5] +
               vars[41] * N[6] + vars[47] * N[7] + vars[53] * N[8]);

    // Compute the derivatives of U along the a-direction
    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8]);

    // Compute the derivative of X along a
    Xd[0] = Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] +
            Xpts[9] * Na[3] + Xpts[12] * Na[4] + Xpts[15] * Na[5] +
            Xpts[18] * Na[6] + Xpts[21] * Na[7] + Xpts[24] * Na[8];
    Xd[1] = Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] +
            Xpts[10] * Na[3] + Xpts[13] * Na[4] + Xpts[16] * Na[5] +
            Xpts[19] * Na[6] + Xpts[22] * Na[7] + Xpts[25] * Na[8];
    Xd[2] = Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] +
            Xpts[11] * Na[3] + Xpts[14] * Na[4] + Xpts[17] * Na[5] +
            Xpts[20] * Na[6] + Xpts[23] * Na[7] + Xpts[26] * Na[8];

    // Compute the derivatives along the b-direction
    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8]);

    // Compute the derivative of X along a
    Xd[3] = Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] +
            Xpts[9] * Nb[3] + Xpts[12] * Nb[4] + Xpts[15] * Nb[5] +
            Xpts[18] * Nb[6] + Xpts[21] * Nb[7] + Xpts[24] * Nb[8];
    Xd[4] = Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] +
            Xpts[10] * Nb[3] + Xpts[13] * Nb[4] + Xpts[16] * Nb[5] +
            Xpts[19] * Nb[6] + Xpts[22] * Nb[7] + Xpts[25] * Nb[8];
    Xd[5] = Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] +
            Xpts[11] * Nb[3] + Xpts[14] * Nb[4] + Xpts[17] * Nb[5] +
            Xpts[20] * Nb[6] + Xpts[23] * Nb[7] + Xpts[26] * Nb[8];
  } else if (order == 4) {
    FElibrary::biLagrangeSF(N, Na, Nb, pt, 4);

    // Compute the shape functions and rotations
    Urot[0] =
        (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3] +
         vars[27] * N[4] + vars[33] * N[5] + vars[39] * N[6] + vars[45] * N[7] +
         vars[51] * N[8] + vars[57] * N[9] + vars[63] * N[10] +
         vars[69] * N[11] + vars[75] * N[12] + vars[81] * N[13] +
         vars[87] * N[14] + vars[93] * N[15]);
    Urot[1] =
        (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3] +
         vars[28] * N[4] + vars[34] * N[5] + vars[40] * N[6] + vars[46] * N[7] +
         vars[52] * N[8] + vars[58] * N[9] + vars[64] * N[10] +
         vars[70] * N[11] + vars[76] * N[12] + vars[82] * N[13] +
         vars[88] * N[14] + vars[94] * N[15]);
    Urot[2] =
        (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3] +
         vars[29] * N[4] + vars[35] * N[5] + vars[41] * N[6] + vars[47] * N[7] +
         vars[53] * N[8] + vars[59] * N[9] + vars[65] * N[10] +
         vars[71] * N[11] + vars[77] * N[12] + vars[83] * N[13] +
         vars[89] * N[14] + vars[95] * N[15]);

    // Compute the derivatives of the shape functions and the
    // derivatives along a
    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8] +
             vars[54] * Na[9] + vars[60] * Na[10] + vars[66] * Na[11] +
             vars[72] * Na[12] + vars[78] * Na[13] + vars[84] * Na[14] +
             vars[90] * Na[15]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8] +
             vars[55] * Na[9] + vars[61] * Na[10] + vars[67] * Na[11] +
             vars[73] * Na[12] + vars[79] * Na[13] + vars[85] * Na[14] +
             vars[91] * Na[15]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8] +
             vars[56] * Na[9] + vars[62] * Na[10] + vars[68] * Na[11] +
             vars[74] * Na[12] + vars[80] * Na[13] + vars[86] * Na[14] +
             vars[92] * Na[15]);

    Xd[0] =
        (Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] + Xpts[9] * Na[3] +
         Xpts[12] * Na[4] + Xpts[15] * Na[5] + Xpts[18] * Na[6] +
         Xpts[21] * Na[7] + Xpts[24] * Na[8] + Xpts[27] * Na[9] +
         Xpts[30] * Na[10] + Xpts[33] * Na[11] + Xpts[36] * Na[12] +
         Xpts[39] * Na[13] + Xpts[42] * Na[14] + Xpts[45] * Na[15]);
    Xd[1] = (Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] +
             Xpts[10] * Na[3] + Xpts[13] * Na[4] + Xpts[16] * Na[5] +
             Xpts[19] * Na[6] + Xpts[22] * Na[7] + Xpts[25] * Na[8] +
             Xpts[28] * Na[9] + Xpts[31] * Na[10] + Xpts[34] * Na[11] +
             Xpts[37] * Na[12] + Xpts[40] * Na[13] + Xpts[43] * Na[14] +
             Xpts[46] * Na[15]);
    Xd[2] = (Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] +
             Xpts[11] * Na[3] + Xpts[14] * Na[4] + Xpts[17] * Na[5] +
             Xpts[20] * Na[6] + Xpts[23] * Na[7] + Xpts[26] * Na[8] +
             Xpts[29] * Na[9] + Xpts[32] * Na[10] + Xpts[35] * Na[11] +
             Xpts[38] * Na[12] + Xpts[41] * Na[13] + Xpts[44] * Na[14] +
             Xpts[47] * Na[15]);

    // Compute the derivatives of the shape functions along the
    // parametric b-direction
    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8] +
             vars[54] * Nb[9] + vars[60] * Nb[10] + vars[66] * Nb[11] +
             vars[72] * Nb[12] + vars[78] * Nb[13] + vars[84] * Nb[14] +
             vars[90] * Nb[15]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8] +
             vars[55] * Nb[9] + vars[61] * Nb[10] + vars[67] * Nb[11] +
             vars[73] * Nb[12] + vars[79] * Nb[13] + vars[85] * Nb[14] +
             vars[91] * Nb[15]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8] +
             vars[56] * Nb[9] + vars[62] * Nb[10] + vars[68] * Nb[11] +
             vars[74] * Nb[12] + vars[80] * Nb[13] + vars[86] * Nb[14] +
             vars[92] * Nb[15]);

    Xd[3] =
        (Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] + Xpts[9] * Nb[3] +
         Xpts[12] * Nb[4] + Xpts[15] * Nb[5] + Xpts[18] * Nb[6] +
         Xpts[21] * Nb[7] + Xpts[24] * Nb[8] + Xpts[27] * Nb[9] +
         Xpts[30] * Nb[10] + Xpts[33] * Nb[11] + Xpts[36] * Nb[12] +
         Xpts[39] * Nb[13] + Xpts[42] * Nb[14] + Xpts[45] * Nb[15]);
    Xd[4] = (Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] +
             Xpts[10] * Nb[3] + Xpts[13] * Nb[4] + Xpts[16] * Nb[5] +
             Xpts[19] * Nb[6] + Xpts[22] * Nb[7] + Xpts[25] * Nb[8] +
             Xpts[28] * Nb[9] + Xpts[31] * Nb[10] + Xpts[34] * Nb[11] +
             Xpts[37] * Nb[12] + Xpts[40] * Nb[13] + Xpts[43] * Nb[14] +
             Xpts[46] * Nb[15]);
    Xd[5] = (Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] +
             Xpts[11] * Nb[3] + Xpts[14] * Nb[4] + Xpts[17] * Nb[5] +
             Xpts[20] * Nb[6] + Xpts[23] * Nb[7] + Xpts[26] * Nb[8] +
             Xpts[29] * Nb[9] + Xpts[32] * Nb[10] + Xpts[35] * Nb[11] +
             Xpts[38] * Nb[12] + Xpts[41] * Nb[13] + Xpts[44] * Nb[14] +
             Xpts[47] * Nb[15]);
  } else {
    FElibrary::biLagrangeSF(N, Na, Nb, pt, order);

    Urot[0] = Urot[1] = Urot[2] = 0.0;
    Ud[0] = Ud[1] = Ud[2] = Ud[3] = Ud[4] = Ud[5] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = Xd[3] = Xd[4] = Xd[5] = 0.0;

    const int order2 = order * order;
    for (int i = 0; i < order2; i++) {
      Urot[0] += vars[3] * N[0];
      Urot[1] += vars[4] * N[0];
      Urot[2] += vars[5] * N[0];

      // First derivatives of A and B
      Ud[0] += vars[0] * Na[0];
      Ud[2] += vars[1] * Na[0];
      Ud[4] += vars[2] * Na[0];
      Xd[0] += Xpts[0] * Na[0];
      Xd[1] += Xpts[1] * Na[0];
      Xd[2] += Xpts[2] * Na[0];

      Ud[1] += vars[0] * Nb[0];
      Ud[3] += vars[1] * Nb[0];
      Ud[5] += vars[2] * Nb[0];
      Xd[3] += Xpts[0] * Nb[0];
      Xd[4] += Xpts[1] * Nb[0];
      Xd[5] += Xpts[2] * Nb[0];

      N++;
      Na++;
      Nb++;
      Xpts += 3;
      vars += 6;
    }
  }
}

/*
  Compute the displacement within the shell

  input:
  num_nodes:   the number of nodes
  vars:        the values of the state variables
  N:           the shape functions

  output:
  U:           the displacements at the parametric point withinn the element
*/
void compute_shell_U(const int num_nodes, TacsScalar U[],
                     const TacsScalar vars[], const double N[]) {
  if (num_nodes == 4) {
    U[0] = vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] + vars[18] * N[3];
    U[1] = vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] + vars[19] * N[3];
    U[2] = vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] + vars[20] * N[3];
    U[3] = vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3];
    U[4] = vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3];
    U[5] = vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3];
  } else if (num_nodes == 9) {
    U[0] = (vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] +
            vars[18] * N[3] + vars[24] * N[4] + vars[30] * N[5] +
            vars[36] * N[6] + vars[42] * N[7] + vars[48] * N[8]);
    U[1] = (vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] +
            vars[19] * N[3] + vars[25] * N[4] + vars[31] * N[5] +
            vars[37] * N[6] + vars[43] * N[7] + vars[49] * N[8]);
    U[2] = (vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] +
            vars[20] * N[3] + vars[26] * N[4] + vars[32] * N[5] +
            vars[38] * N[6] + vars[44] * N[7] + vars[50] * N[8]);
    U[3] = (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] +
            vars[21] * N[3] + vars[27] * N[4] + vars[33] * N[5] +
            vars[39] * N[6] + vars[45] * N[7] + vars[51] * N[8]);
    U[4] = (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] +
            vars[22] * N[3] + vars[28] * N[4] + vars[34] * N[5] +
            vars[40] * N[6] + vars[46] * N[7] + vars[52] * N[8]);
    U[5] = (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] +
            vars[23] * N[3] + vars[29] * N[4] + vars[35] * N[5] +
            vars[41] * N[6] + vars[47] * N[7] + vars[53] * N[8]);
  } else if (num_nodes == 16) {
    U[0] =
        (vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] + vars[18] * N[3] +
         vars[24] * N[4] + vars[30] * N[5] + vars[36] * N[6] + vars[42] * N[7] +
         vars[48] * N[8] + vars[54] * N[9] + vars[60] * N[10] +
         vars[66] * N[11] + vars[72] * N[12] + vars[78] * N[13] +
         vars[84] * N[14] + vars[90] * N[15]);
    U[1] =
        (vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] + vars[19] * N[3] +
         vars[25] * N[4] + vars[31] * N[5] + vars[37] * N[6] + vars[43] * N[7] +
         vars[49] * N[8] + vars[55] * N[9] + vars[61] * N[10] +
         vars[67] * N[11] + vars[73] * N[12] + vars[79] * N[13] +
         vars[85] * N[14] + vars[91] * N[15]);
    U[2] =
        (vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] + vars[20] * N[3] +
         vars[26] * N[4] + vars[32] * N[5] + vars[38] * N[6] + vars[44] * N[7] +
         vars[50] * N[8] + vars[56] * N[9] + vars[62] * N[10] +
         vars[68] * N[11] + vars[74] * N[12] + vars[80] * N[13] +
         vars[86] * N[14] + vars[92] * N[15]);
    U[3] =
        (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3] +
         vars[27] * N[4] + vars[33] * N[5] + vars[39] * N[6] + vars[45] * N[7] +
         vars[51] * N[8] + vars[57] * N[9] + vars[63] * N[10] +
         vars[69] * N[11] + vars[75] * N[12] + vars[81] * N[13] +
         vars[87] * N[14] + vars[93] * N[15]);
    U[4] =
        (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3] +
         vars[28] * N[4] + vars[34] * N[5] + vars[40] * N[6] + vars[46] * N[7] +
         vars[52] * N[8] + vars[58] * N[9] + vars[64] * N[10] +
         vars[70] * N[11] + vars[76] * N[12] + vars[82] * N[13] +
         vars[88] * N[14] + vars[94] * N[15]);
    U[5] =
        (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3] +
         vars[29] * N[4] + vars[35] * N[5] + vars[41] * N[6] + vars[47] * N[7] +
         vars[53] * N[8] + vars[59] * N[9] + vars[65] * N[10] +
         vars[71] * N[11] + vars[77] * N[12] + vars[83] * N[13] +
         vars[89] * N[14] + vars[95] * N[15]);
  } else {
    U[0] = U[1] = U[2] = U[3] = U[4] = U[5] = 0.0;

    const TacsScalar *v = vars;
    for (int i = 0; i < num_nodes; i++) {
      U[0] += v[0] * N[0];
      U[1] += v[1] * N[0];
      U[2] += v[2] * N[0];
      U[3] += v[3] * N[0];
      U[4] += v[4] * N[0];
      U[5] += v[5] * N[0];
      v += 6;
      N++;
    }
  }
}

/*
  Compute the displacement and the derivative of the displacement
  along the parametric directions

  input:
  num_nodes: the number of nodes
  vars: the x,y,z displacements and rotations for each node
  N: the shape functions
  Na: the derivative of the shape functions along the first
  parametric direction
  Nb: the derivative of the shape functions along the second
  parametric direction

  output:
  U: the displacements and rotations at the parametric point
  Ud: the derivative of the displacements along the parametric directions
*/
void compute_shell_Ud(const int num_nodes, TacsScalar U[], TacsScalar Ud[],
                      const TacsScalar vars[], const double N[],
                      const double Na[], const double Nb[]) {
  if (num_nodes == 4) {
    U[0] = vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] + vars[18] * N[3];
    U[1] = vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] + vars[19] * N[3];
    U[2] = vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] + vars[20] * N[3];
    U[3] = vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3];
    U[4] = vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3];
    U[5] = vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3];

    Ud[0] =
        vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] + vars[18] * Na[3];
    Ud[2] =
        vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] + vars[19] * Na[3];
    Ud[4] =
        vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] + vars[20] * Na[3];
    Ud[6] =
        vars[3] * Na[0] + vars[9] * Na[1] + vars[15] * Na[2] + vars[21] * Na[3];
    Ud[8] = vars[4] * Na[0] + vars[10] * Na[1] + vars[16] * Na[2] +
            vars[22] * Na[3];
    Ud[10] = vars[5] * Na[0] + vars[11] * Na[1] + vars[17] * Na[2] +
             vars[23] * Na[3];

    Ud[1] =
        vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] + vars[18] * Nb[3];
    Ud[3] =
        vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] + vars[19] * Nb[3];
    Ud[5] =
        vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] + vars[20] * Nb[3];
    Ud[7] =
        vars[3] * Nb[0] + vars[9] * Nb[1] + vars[15] * Nb[2] + vars[21] * Nb[3];
    Ud[9] = vars[4] * Nb[0] + vars[10] * Nb[1] + vars[16] * Nb[2] +
            vars[22] * Nb[3];
    Ud[11] = vars[5] * Nb[0] + vars[11] * Nb[1] + vars[17] * Nb[2] +
             vars[23] * Nb[3];
  } else if (num_nodes == 9) {
    U[0] = (vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] +
            vars[18] * N[3] + vars[24] * N[4] + vars[30] * N[5] +
            vars[36] * N[6] + vars[42] * N[7] + vars[48] * N[8]);
    U[1] = (vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] +
            vars[19] * N[3] + vars[25] * N[4] + vars[31] * N[5] +
            vars[37] * N[6] + vars[43] * N[7] + vars[49] * N[8]);
    U[2] = (vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] +
            vars[20] * N[3] + vars[26] * N[4] + vars[32] * N[5] +
            vars[38] * N[6] + vars[44] * N[7] + vars[50] * N[8]);
    U[3] = (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] +
            vars[21] * N[3] + vars[27] * N[4] + vars[33] * N[5] +
            vars[39] * N[6] + vars[45] * N[7] + vars[51] * N[8]);
    U[4] = (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] +
            vars[22] * N[3] + vars[28] * N[4] + vars[34] * N[5] +
            vars[40] * N[6] + vars[46] * N[7] + vars[52] * N[8]);
    U[5] = (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] +
            vars[23] * N[3] + vars[29] * N[4] + vars[35] * N[5] +
            vars[41] * N[6] + vars[47] * N[7] + vars[53] * N[8]);

    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8]);
    Ud[6] = (vars[3] * Na[0] + vars[9] * Na[1] + vars[15] * Na[2] +
             vars[21] * Na[3] + vars[27] * Na[4] + vars[33] * Na[5] +
             vars[39] * Na[6] + vars[45] * Na[7] + vars[51] * Na[8]);
    Ud[8] = (vars[4] * Na[0] + vars[10] * Na[1] + vars[16] * Na[2] +
             vars[22] * Na[3] + vars[28] * Na[4] + vars[34] * Na[5] +
             vars[40] * Na[6] + vars[46] * Na[7] + vars[52] * Na[8]);
    Ud[10] = (vars[5] * Na[0] + vars[11] * Na[1] + vars[17] * Na[2] +
              vars[23] * Na[3] + vars[29] * Na[4] + vars[35] * Na[5] +
              vars[41] * Na[6] + vars[47] * Na[7] + vars[53] * Na[8]);

    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8]);
    Ud[7] = (vars[3] * Nb[0] + vars[9] * Nb[1] + vars[15] * Nb[2] +
             vars[21] * Nb[3] + vars[27] * Nb[4] + vars[33] * Nb[5] +
             vars[39] * Nb[6] + vars[45] * Nb[7] + vars[51] * Nb[8]);
    Ud[9] = (vars[4] * Nb[0] + vars[10] * Nb[1] + vars[16] * Nb[2] +
             vars[22] * Nb[3] + vars[28] * Nb[4] + vars[34] * Nb[5] +
             vars[40] * Nb[6] + vars[46] * Nb[7] + vars[52] * Nb[8]);
    Ud[11] = (vars[5] * Nb[0] + vars[11] * Nb[1] + vars[17] * Nb[2] +
              vars[23] * Nb[3] + vars[29] * Nb[4] + vars[35] * Nb[5] +
              vars[41] * Nb[6] + vars[47] * Nb[7] + vars[53] * Nb[8]);
  } else if (num_nodes == 16) {
    U[0] =
        (vars[0] * N[0] + vars[6] * N[1] + vars[12] * N[2] + vars[18] * N[3] +
         vars[24] * N[4] + vars[30] * N[5] + vars[36] * N[6] + vars[42] * N[7] +
         vars[48] * N[8] + vars[54] * N[9] + vars[60] * N[10] +
         vars[66] * N[11] + vars[72] * N[12] + vars[78] * N[13] +
         vars[84] * N[14] + vars[90] * N[15]);
    U[1] =
        (vars[1] * N[0] + vars[7] * N[1] + vars[13] * N[2] + vars[19] * N[3] +
         vars[25] * N[4] + vars[31] * N[5] + vars[37] * N[6] + vars[43] * N[7] +
         vars[49] * N[8] + vars[55] * N[9] + vars[61] * N[10] +
         vars[67] * N[11] + vars[73] * N[12] + vars[79] * N[13] +
         vars[85] * N[14] + vars[91] * N[15]);
    U[2] =
        (vars[2] * N[0] + vars[8] * N[1] + vars[14] * N[2] + vars[20] * N[3] +
         vars[26] * N[4] + vars[32] * N[5] + vars[38] * N[6] + vars[44] * N[7] +
         vars[50] * N[8] + vars[56] * N[9] + vars[62] * N[10] +
         vars[68] * N[11] + vars[74] * N[12] + vars[80] * N[13] +
         vars[86] * N[14] + vars[92] * N[15]);
    U[3] =
        (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3] +
         vars[27] * N[4] + vars[33] * N[5] + vars[39] * N[6] + vars[45] * N[7] +
         vars[51] * N[8] + vars[57] * N[9] + vars[63] * N[10] +
         vars[69] * N[11] + vars[75] * N[12] + vars[81] * N[13] +
         vars[87] * N[14] + vars[93] * N[15]);
    U[4] =
        (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3] +
         vars[28] * N[4] + vars[34] * N[5] + vars[40] * N[6] + vars[46] * N[7] +
         vars[52] * N[8] + vars[58] * N[9] + vars[64] * N[10] +
         vars[70] * N[11] + vars[76] * N[12] + vars[82] * N[13] +
         vars[88] * N[14] + vars[94] * N[15]);
    U[5] =
        (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3] +
         vars[29] * N[4] + vars[35] * N[5] + vars[41] * N[6] + vars[47] * N[7] +
         vars[53] * N[8] + vars[59] * N[9] + vars[65] * N[10] +
         vars[71] * N[11] + vars[77] * N[12] + vars[83] * N[13] +
         vars[89] * N[14] + vars[95] * N[15]);

    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8] +
             vars[54] * Na[9] + vars[60] * Na[10] + vars[66] * Na[11] +
             vars[72] * Na[12] + vars[78] * Na[13] + vars[84] * Na[14] +
             vars[90] * Na[15]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8] +
             vars[55] * Na[9] + vars[61] * Na[10] + vars[67] * Na[11] +
             vars[73] * Na[12] + vars[79] * Na[13] + vars[85] * Na[14] +
             vars[91] * Na[15]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8] +
             vars[56] * Na[9] + vars[62] * Na[10] + vars[68] * Na[11] +
             vars[74] * Na[12] + vars[80] * Na[13] + vars[86] * Na[14] +
             vars[92] * Na[15]);
    Ud[6] = (vars[3] * Na[0] + vars[9] * Na[1] + vars[15] * Na[2] +
             vars[21] * Na[3] + vars[27] * Na[4] + vars[33] * Na[5] +
             vars[39] * Na[6] + vars[45] * Na[7] + vars[51] * Na[8] +
             vars[57] * Na[9] + vars[63] * Na[10] + vars[69] * Na[11] +
             vars[75] * Na[12] + vars[81] * Na[13] + vars[87] * Na[14] +
             vars[93] * Na[15]);
    Ud[8] = (vars[4] * Na[0] + vars[10] * Na[1] + vars[16] * Na[2] +
             vars[22] * Na[3] + vars[28] * Na[4] + vars[34] * Na[5] +
             vars[40] * Na[6] + vars[46] * Na[7] + vars[52] * Na[8] +
             vars[58] * Na[9] + vars[64] * Na[10] + vars[70] * Na[11] +
             vars[76] * Na[12] + vars[82] * Na[13] + vars[88] * Na[14] +
             vars[94] * Na[15]);
    Ud[10] = (vars[5] * Na[0] + vars[11] * Na[1] + vars[17] * Na[2] +
              vars[23] * Na[3] + vars[29] * Na[4] + vars[35] * Na[5] +
              vars[41] * Na[6] + vars[47] * Na[7] + vars[53] * Na[8] +
              vars[59] * Na[9] + vars[65] * Na[10] + vars[71] * Na[11] +
              vars[77] * Na[12] + vars[83] * Na[13] + vars[89] * Na[14] +
              vars[95] * Na[15]);

    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8] +
             vars[54] * Nb[9] + vars[60] * Nb[10] + vars[66] * Nb[11] +
             vars[72] * Nb[12] + vars[78] * Nb[13] + vars[84] * Nb[14] +
             vars[90] * Nb[15]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8] +
             vars[55] * Nb[9] + vars[61] * Nb[10] + vars[67] * Nb[11] +
             vars[73] * Nb[12] + vars[79] * Nb[13] + vars[85] * Nb[14] +
             vars[91] * Nb[15]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8] +
             vars[56] * Nb[9] + vars[62] * Nb[10] + vars[68] * Nb[11] +
             vars[74] * Nb[12] + vars[80] * Nb[13] + vars[86] * Nb[14] +
             vars[92] * Nb[15]);
    Ud[7] = (vars[3] * Nb[0] + vars[9] * Nb[1] + vars[15] * Nb[2] +
             vars[21] * Nb[3] + vars[27] * Nb[4] + vars[33] * Nb[5] +
             vars[39] * Nb[6] + vars[45] * Nb[7] + vars[51] * Nb[8] +
             vars[57] * Nb[9] + vars[63] * Nb[10] + vars[69] * Nb[11] +
             vars[75] * Nb[12] + vars[81] * Nb[13] + vars[87] * Nb[14] +
             vars[93] * Nb[15]);
    Ud[9] = (vars[4] * Nb[0] + vars[10] * Nb[1] + vars[16] * Nb[2] +
             vars[22] * Nb[3] + vars[28] * Nb[4] + vars[34] * Nb[5] +
             vars[40] * Nb[6] + vars[46] * Nb[7] + vars[52] * Nb[8] +
             vars[58] * Nb[9] + vars[64] * Nb[10] + vars[70] * Nb[11] +
             vars[76] * Nb[12] + vars[82] * Nb[13] + vars[88] * Nb[14] +
             vars[94] * Nb[15]);
    Ud[11] = (vars[5] * Nb[0] + vars[11] * Nb[1] + vars[17] * Nb[2] +
              vars[23] * Nb[3] + vars[29] * Nb[4] + vars[35] * Nb[5] +
              vars[41] * Nb[6] + vars[47] * Nb[7] + vars[53] * Nb[8] +
              vars[59] * Nb[9] + vars[65] * Nb[10] + vars[71] * Nb[11] +
              vars[77] * Nb[12] + vars[83] * Nb[13] + vars[89] * Nb[14] +
              vars[95] * Nb[15]);
  } else {
    U[0] = U[1] = U[2] = U[3] = U[4] = U[5] = 0.0;
    Ud[0] = Ud[1] = Ud[2] = Ud[3] = Ud[4] = Ud[5] = 0.0;
    Ud[6] = Ud[7] = Ud[8] = Ud[9] = Ud[10] = Ud[11] = 0.0;

    const TacsScalar *v = vars;
    for (int i = 0; i < num_nodes; i++) {
      U[0] += v[0] * N[0];
      U[1] += v[1] * N[0];
      U[2] += v[2] * N[0];
      U[3] += v[3] * N[0];
      U[4] += v[4] * N[0];
      U[5] += v[5] * N[0];
      v += 6;
      N++;
    }

    v = vars;
    for (int i = 0; i < num_nodes; i++) {
      Ud[0] += v[0] * Na[0];
      Ud[2] += v[1] * Na[0];
      Ud[4] += v[2] * Na[0];
      Ud[6] += v[3] * Na[0];
      Ud[8] += v[4] * Na[0];
      Ud[10] += v[5] * Na[0];
      v += 6;
      Na++;
    }

    v = vars;
    for (int i = 0; i < num_nodes; i++) {
      Ud[1] += v[0] * Nb[0];
      Ud[3] += v[1] * Nb[0];
      Ud[5] += v[2] * Nb[0];
      Ud[7] += v[3] * Nb[0];
      Ud[9] += v[4] * Nb[0];
      Ud[11] += v[5] * Nb[0];

      v += 6;
      Nb++;
    }
  }
}

/*
  Compute only the terms and derivatives required to compute the
  in-plane and out-of-plane shear tensorial strain components for a
  b-spline shell.

  input:
  order:      the order of the b-spline in both u/v <= 8
  tu, tv:     the knot vector in the u/v directions
  intu, intv: the knot intervals in the u/v directions
  pt:         the parametric point to evaluate the shape functions
  Xpts:       the nodal locations
  vars:       the state varaible values

  output:
  Urot: the rotations at the mid-surface
  Ud: the derivatives of the displacements at the mid-surface
  Xd: the derivatives of the nodal locations along the parametric directions
  N, Na, Nb:  the tensor product shape functions and their derivatives
*/
void bspline_tensorial_components(const int order, TacsScalar Urot[],
                                  TacsScalar Ud[], TacsScalar Xd[], double N[],
                                  double Na[], double Nb[], const double pt[],
                                  int intu, int intv, const double tu[],
                                  const double tv[], const TacsScalar Xpts[],
                                  const TacsScalar vars[]) {
  double Nu[2 * 8], Nv[2 * 8], work[2 * 8 + 8 * 8];
  if (order > 8) {
    return;
  }

  // Re-map the Gauss point to the b-spline knot vector
  double u = 0.5 * (tu[intu] * (1.0 - pt[0]) + tu[intu + 1] * (1.0 + pt[0]));
  double v = 0.5 * (tv[intv] * (1.0 - pt[1]) + tv[intv + 1] * (1.0 + pt[1]));

  // Compute the derivative of pt with respect to the Gauss poitns
  double du = 0.5 * (tu[intu + 1] - tu[intu]);
  double dv = 0.5 * (tv[intv + 1] - tv[intv]);

  // Evaluate the b-spline basis functions
  int ideriv = 1;
  FElibrary::bspline_basis_derivative(Nu, intu, u, ideriv, tu, order, work);
  FElibrary::bspline_basis_derivative(Nv, intv, v, ideriv, tv, order, work);

  double *n = N, *na = Na, *nb = Nb;
  for (int j = 0; j < order; j++) {
    for (int i = 0; i < order; i++) {
      n[0] = Nu[i] * Nv[j];
      na[0] = Nu[order + i] * Nv[j] * du;
      nb[0] = Nu[i] * Nv[order + j] * dv;
      n++;
      na++;
      nb++;
    }
  }

  if (order == 2) {
    // Compute the rotations
    Urot[0] =
        vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3];
    Urot[1] =
        vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3];
    Urot[2] =
        vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3];

    // Compute the derivatives of the displacement
    Ud[0] =
        vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] + vars[18] * Na[3];
    Ud[2] =
        vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] + vars[19] * Na[3];
    Ud[4] =
        vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] + vars[20] * Na[3];

    // Compute the derivative of X along a
    Xd[0] =
        Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] + Xpts[9] * Na[3];
    Xd[1] =
        Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] + Xpts[10] * Na[3];
    Xd[2] =
        Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] + Xpts[11] * Na[3];

    // Compute the derivative of the displacement
    Ud[1] =
        vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] + vars[18] * Nb[3];
    Ud[3] =
        vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] + vars[19] * Nb[3];
    Ud[5] =
        vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] + vars[20] * Nb[3];

    // Compute the derivative of X along b
    Xd[3] =
        Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] + Xpts[9] * Nb[3];
    Xd[4] =
        Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] + Xpts[10] * Nb[3];
    Xd[5] =
        Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] + Xpts[11] * Nb[3];
  } else if (order == 3) {
    // Compute the rotations
    Urot[0] = (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] +
               vars[21] * N[3] + vars[27] * N[4] + vars[33] * N[5] +
               vars[39] * N[6] + vars[45] * N[7] + vars[51] * N[8]);
    Urot[1] = (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] +
               vars[22] * N[3] + vars[28] * N[4] + vars[34] * N[5] +
               vars[40] * N[6] + vars[46] * N[7] + vars[52] * N[8]);
    Urot[2] = (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] +
               vars[23] * N[3] + vars[29] * N[4] + vars[35] * N[5] +
               vars[41] * N[6] + vars[47] * N[7] + vars[53] * N[8]);

    // Compute the derivatives of U along the a-direction
    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8]);

    // Compute the derivative of X along a
    Xd[0] = Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] +
            Xpts[9] * Na[3] + Xpts[12] * Na[4] + Xpts[15] * Na[5] +
            Xpts[18] * Na[6] + Xpts[21] * Na[7] + Xpts[24] * Na[8];
    Xd[1] = Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] +
            Xpts[10] * Na[3] + Xpts[13] * Na[4] + Xpts[16] * Na[5] +
            Xpts[19] * Na[6] + Xpts[22] * Na[7] + Xpts[25] * Na[8];
    Xd[2] = Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] +
            Xpts[11] * Na[3] + Xpts[14] * Na[4] + Xpts[17] * Na[5] +
            Xpts[20] * Na[6] + Xpts[23] * Na[7] + Xpts[26] * Na[8];

    // Compute the derivatives along the b-direction
    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8]);

    // Compute the derivative of X along a
    Xd[3] = Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] +
            Xpts[9] * Nb[3] + Xpts[12] * Nb[4] + Xpts[15] * Nb[5] +
            Xpts[18] * Nb[6] + Xpts[21] * Nb[7] + Xpts[24] * Nb[8];
    Xd[4] = Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] +
            Xpts[10] * Nb[3] + Xpts[13] * Nb[4] + Xpts[16] * Nb[5] +
            Xpts[19] * Nb[6] + Xpts[22] * Nb[7] + Xpts[25] * Nb[8];
    Xd[5] = Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] +
            Xpts[11] * Nb[3] + Xpts[14] * Nb[4] + Xpts[17] * Nb[5] +
            Xpts[20] * Nb[6] + Xpts[23] * Nb[7] + Xpts[26] * Nb[8];
  } else if (order == 4) {
    // Compute the shape functions and rotations
    Urot[0] =
        (vars[3] * N[0] + vars[9] * N[1] + vars[15] * N[2] + vars[21] * N[3] +
         vars[27] * N[4] + vars[33] * N[5] + vars[39] * N[6] + vars[45] * N[7] +
         vars[51] * N[8] + vars[57] * N[9] + vars[63] * N[10] +
         vars[69] * N[11] + vars[75] * N[12] + vars[81] * N[13] +
         vars[87] * N[14] + vars[93] * N[15]);
    Urot[1] =
        (vars[4] * N[0] + vars[10] * N[1] + vars[16] * N[2] + vars[22] * N[3] +
         vars[28] * N[4] + vars[34] * N[5] + vars[40] * N[6] + vars[46] * N[7] +
         vars[52] * N[8] + vars[58] * N[9] + vars[64] * N[10] +
         vars[70] * N[11] + vars[76] * N[12] + vars[82] * N[13] +
         vars[88] * N[14] + vars[94] * N[15]);
    Urot[2] =
        (vars[5] * N[0] + vars[11] * N[1] + vars[17] * N[2] + vars[23] * N[3] +
         vars[29] * N[4] + vars[35] * N[5] + vars[41] * N[6] + vars[47] * N[7] +
         vars[53] * N[8] + vars[59] * N[9] + vars[65] * N[10] +
         vars[71] * N[11] + vars[77] * N[12] + vars[83] * N[13] +
         vars[89] * N[14] + vars[95] * N[15]);

    // Compute the derivatives of the shape functions and the
    // derivatives along a
    Ud[0] = (vars[0] * Na[0] + vars[6] * Na[1] + vars[12] * Na[2] +
             vars[18] * Na[3] + vars[24] * Na[4] + vars[30] * Na[5] +
             vars[36] * Na[6] + vars[42] * Na[7] + vars[48] * Na[8] +
             vars[54] * Na[9] + vars[60] * Na[10] + vars[66] * Na[11] +
             vars[72] * Na[12] + vars[78] * Na[13] + vars[84] * Na[14] +
             vars[90] * Na[15]);
    Ud[2] = (vars[1] * Na[0] + vars[7] * Na[1] + vars[13] * Na[2] +
             vars[19] * Na[3] + vars[25] * Na[4] + vars[31] * Na[5] +
             vars[37] * Na[6] + vars[43] * Na[7] + vars[49] * Na[8] +
             vars[55] * Na[9] + vars[61] * Na[10] + vars[67] * Na[11] +
             vars[73] * Na[12] + vars[79] * Na[13] + vars[85] * Na[14] +
             vars[91] * Na[15]);
    Ud[4] = (vars[2] * Na[0] + vars[8] * Na[1] + vars[14] * Na[2] +
             vars[20] * Na[3] + vars[26] * Na[4] + vars[32] * Na[5] +
             vars[38] * Na[6] + vars[44] * Na[7] + vars[50] * Na[8] +
             vars[56] * Na[9] + vars[62] * Na[10] + vars[68] * Na[11] +
             vars[74] * Na[12] + vars[80] * Na[13] + vars[86] * Na[14] +
             vars[92] * Na[15]);

    Xd[0] =
        (Xpts[0] * Na[0] + Xpts[3] * Na[1] + Xpts[6] * Na[2] + Xpts[9] * Na[3] +
         Xpts[12] * Na[4] + Xpts[15] * Na[5] + Xpts[18] * Na[6] +
         Xpts[21] * Na[7] + Xpts[24] * Na[8] + Xpts[27] * Na[9] +
         Xpts[30] * Na[10] + Xpts[33] * Na[11] + Xpts[36] * Na[12] +
         Xpts[39] * Na[13] + Xpts[42] * Na[14] + Xpts[45] * Na[15]);
    Xd[1] = (Xpts[1] * Na[0] + Xpts[4] * Na[1] + Xpts[7] * Na[2] +
             Xpts[10] * Na[3] + Xpts[13] * Na[4] + Xpts[16] * Na[5] +
             Xpts[19] * Na[6] + Xpts[22] * Na[7] + Xpts[25] * Na[8] +
             Xpts[28] * Na[9] + Xpts[31] * Na[10] + Xpts[34] * Na[11] +
             Xpts[37] * Na[12] + Xpts[40] * Na[13] + Xpts[43] * Na[14] +
             Xpts[46] * Na[15]);
    Xd[2] = (Xpts[2] * Na[0] + Xpts[5] * Na[1] + Xpts[8] * Na[2] +
             Xpts[11] * Na[3] + Xpts[14] * Na[4] + Xpts[17] * Na[5] +
             Xpts[20] * Na[6] + Xpts[23] * Na[7] + Xpts[26] * Na[8] +
             Xpts[29] * Na[9] + Xpts[32] * Na[10] + Xpts[35] * Na[11] +
             Xpts[38] * Na[12] + Xpts[41] * Na[13] + Xpts[44] * Na[14] +
             Xpts[47] * Na[15]);

    // Compute the derivatives of the shape functions along the
    // parametric b-direction
    Ud[1] = (vars[0] * Nb[0] + vars[6] * Nb[1] + vars[12] * Nb[2] +
             vars[18] * Nb[3] + vars[24] * Nb[4] + vars[30] * Nb[5] +
             vars[36] * Nb[6] + vars[42] * Nb[7] + vars[48] * Nb[8] +
             vars[54] * Nb[9] + vars[60] * Nb[10] + vars[66] * Nb[11] +
             vars[72] * Nb[12] + vars[78] * Nb[13] + vars[84] * Nb[14] +
             vars[90] * Nb[15]);
    Ud[3] = (vars[1] * Nb[0] + vars[7] * Nb[1] + vars[13] * Nb[2] +
             vars[19] * Nb[3] + vars[25] * Nb[4] + vars[31] * Nb[5] +
             vars[37] * Nb[6] + vars[43] * Nb[7] + vars[49] * Nb[8] +
             vars[55] * Nb[9] + vars[61] * Nb[10] + vars[67] * Nb[11] +
             vars[73] * Nb[12] + vars[79] * Nb[13] + vars[85] * Nb[14] +
             vars[91] * Nb[15]);
    Ud[5] = (vars[2] * Nb[0] + vars[8] * Nb[1] + vars[14] * Nb[2] +
             vars[20] * Nb[3] + vars[26] * Nb[4] + vars[32] * Nb[5] +
             vars[38] * Nb[6] + vars[44] * Nb[7] + vars[50] * Nb[8] +
             vars[56] * Nb[9] + vars[62] * Nb[10] + vars[68] * Nb[11] +
             vars[74] * Nb[12] + vars[80] * Nb[13] + vars[86] * Nb[14] +
             vars[92] * Nb[15]);

    Xd[3] =
        (Xpts[0] * Nb[0] + Xpts[3] * Nb[1] + Xpts[6] * Nb[2] + Xpts[9] * Nb[3] +
         Xpts[12] * Nb[4] + Xpts[15] * Nb[5] + Xpts[18] * Nb[6] +
         Xpts[21] * Nb[7] + Xpts[24] * Nb[8] + Xpts[27] * Nb[9] +
         Xpts[30] * Nb[10] + Xpts[33] * Nb[11] + Xpts[36] * Nb[12] +
         Xpts[39] * Nb[13] + Xpts[42] * Nb[14] + Xpts[45] * Nb[15]);
    Xd[4] = (Xpts[1] * Nb[0] + Xpts[4] * Nb[1] + Xpts[7] * Nb[2] +
             Xpts[10] * Nb[3] + Xpts[13] * Nb[4] + Xpts[16] * Nb[5] +
             Xpts[19] * Nb[6] + Xpts[22] * Nb[7] + Xpts[25] * Nb[8] +
             Xpts[28] * Nb[9] + Xpts[31] * Nb[10] + Xpts[34] * Nb[11] +
             Xpts[37] * Nb[12] + Xpts[40] * Nb[13] + Xpts[43] * Nb[14] +
             Xpts[46] * Nb[15]);
    Xd[5] = (Xpts[2] * Nb[0] + Xpts[5] * Nb[1] + Xpts[8] * Nb[2] +
             Xpts[11] * Nb[3] + Xpts[14] * Nb[4] + Xpts[17] * Nb[5] +
             Xpts[20] * Nb[6] + Xpts[23] * Nb[7] + Xpts[26] * Nb[8] +
             Xpts[29] * Nb[9] + Xpts[32] * Nb[10] + Xpts[35] * Nb[11] +
             Xpts[38] * Nb[12] + Xpts[41] * Nb[13] + Xpts[44] * Nb[14] +
             Xpts[47] * Nb[15]);
  } else {
    Urot[0] = Urot[1] = Urot[2] = 0.0;
    Ud[0] = Ud[1] = Ud[2] = Ud[3] = Ud[4] = Ud[5] = 0.0;
    Xd[0] = Xd[1] = Xd[2] = Xd[3] = Xd[4] = Xd[5] = 0.0;

    const int order2 = order * order;
    for (int i = 0; i < order2; i++) {
      Urot[0] += vars[3] * N[0];
      Urot[1] += vars[4] * N[0];
      Urot[2] += vars[5] * N[0];

      // First derivatives of A and B
      Ud[0] += vars[0] * Na[0];
      Ud[2] += vars[1] * Na[0];
      Ud[4] += vars[2] * Na[0];
      Xd[0] += Xpts[0] * Na[0];
      Xd[1] += Xpts[1] * Na[0];
      Xd[2] += Xpts[2] * Na[0];

      Ud[1] += vars[0] * Nb[0];
      Ud[3] += vars[1] * Nb[0];
      Ud[5] += vars[2] * Nb[0];
      Xd[3] += Xpts[0] * Nb[0];
      Xd[4] += Xpts[1] * Nb[0];
      Xd[5] += Xpts[2] * Nb[0];

      N++;
      Na++;
      Nb++;
      Xpts += 3;
      vars += 6;
    }
  }
}

/*
  Compute the position and first derivatives along the shell surface.
  At the same time, compute the shape functions.

  The shell surface is based on the shape functions and nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd: position and first derivative of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
*/
void shell_jacobian(const int order, TacsScalar X[], TacsScalar Xd[],
                    double N[], double Na[], double Nb[], const double gpt[],
                    const TacsScalar Xpts[]) {
  double na[8], nb[8];
  double dna[8], dnb[8];

  X[0] = X[1] = X[2] = TacsScalar(0.0);
  Xd[0] = Xd[1] = Xd[2] = TacsScalar(0.0);
  Xd[3] = Xd[4] = Xd[5] = TacsScalar(0.0);

  if (order <= 8) {
    if (order == 2) {
      FElibrary::lagrangeSF(na, dna, gpt[0], 2);
      FElibrary::lagrangeSF(nb, dnb, gpt[1], 2);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[0] * nb[1];
      N[3] = na[1] * nb[1];

      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[0] * nb[1];
      Na[3] = dna[1] * nb[1];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[0] * dnb[1];
      Nb[3] = na[1] * dnb[1];
    } else if (order == 3) {
      FElibrary::lagrangeSF(na, dna, gpt[0], 3);
      FElibrary::lagrangeSF(nb, dnb, gpt[1], 3);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[2] * nb[0];
      N[3] = na[0] * nb[1];
      N[4] = na[1] * nb[1];
      N[5] = na[2] * nb[1];
      N[6] = na[0] * nb[2];
      N[7] = na[1] * nb[2];
      N[8] = na[2] * nb[2];

      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[2] * nb[0];
      Na[3] = dna[0] * nb[1];
      Na[4] = dna[1] * nb[1];
      Na[5] = dna[2] * nb[1];
      Na[6] = dna[0] * nb[2];
      Na[7] = dna[1] * nb[2];
      Na[8] = dna[2] * nb[2];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[2] * dnb[0];
      Nb[3] = na[0] * dnb[1];
      Nb[4] = na[1] * dnb[1];
      Nb[5] = na[2] * dnb[1];
      Nb[6] = na[0] * dnb[2];
      Nb[7] = na[1] * dnb[2];
      Nb[8] = na[2] * dnb[2];
    } else if (order == 4) {
      FElibrary::lagrangeSF(na, dna, gpt[0], 4);
      FElibrary::lagrangeSF(nb, dnb, gpt[1], 4);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[2] * nb[0];
      N[3] = na[3] * nb[0];
      N[4] = na[0] * nb[1];
      N[5] = na[1] * nb[1];
      N[6] = na[2] * nb[1];
      N[7] = na[3] * nb[1];
      N[8] = na[0] * nb[2];
      N[9] = na[1] * nb[2];
      N[10] = na[2] * nb[2];
      N[11] = na[3] * nb[2];
      N[12] = na[0] * nb[3];
      N[13] = na[1] * nb[3];
      N[14] = na[2] * nb[3];
      N[15] = na[3] * nb[3];

      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[2] * nb[0];
      Na[3] = dna[3] * nb[0];
      Na[4] = dna[0] * nb[1];
      Na[5] = dna[1] * nb[1];
      Na[6] = dna[2] * nb[1];
      Na[7] = dna[3] * nb[1];
      Na[8] = dna[0] * nb[2];
      Na[9] = dna[1] * nb[2];
      Na[10] = dna[2] * nb[2];
      Na[11] = dna[3] * nb[2];
      Na[12] = dna[0] * nb[3];
      Na[13] = dna[1] * nb[3];
      Na[14] = dna[2] * nb[3];
      Na[15] = dna[3] * nb[3];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[2] * dnb[0];
      Nb[3] = na[3] * dnb[0];
      Nb[4] = na[0] * dnb[1];
      Nb[5] = na[1] * dnb[1];
      Nb[6] = na[2] * dnb[1];
      Nb[7] = na[3] * dnb[1];
      Nb[8] = na[0] * dnb[2];
      Nb[9] = na[1] * dnb[2];
      Nb[10] = na[2] * dnb[2];
      Nb[11] = na[3] * dnb[2];
      Nb[12] = na[0] * dnb[3];
      Nb[13] = na[1] * dnb[3];
      Nb[14] = na[2] * dnb[3];
      Nb[15] = na[3] * dnb[3];
    } else {
      FElibrary::lagrangeSF(na, dna, gpt[0], order);
      FElibrary::lagrangeSF(nb, dnb, gpt[1], order);

      double *Np = N;
      double *Nap = Na;
      double *Nbp = Nb;

      for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
          Np[0] = na[i] * nb[j];
          Nap[0] = dna[i] * nb[j];
          Nbp[0] = na[i] * dnb[j];
          Np++;
          Nap++;
          Nbp++;
        }
      }
    }

    const int order2 = order * order;
    for (int i = 0; i < order2; i++) {
      // Evaluate the derivatives
      X[0] += Xpts[0] * N[0];
      X[1] += Xpts[1] * N[0];
      X[2] += Xpts[2] * N[0];

      // First derivatives of A and B
      Xd[0] += Xpts[0] * Na[0];
      Xd[1] += Xpts[1] * Na[0];
      Xd[2] += Xpts[2] * Na[0];

      Xd[3] += Xpts[0] * Nb[0];
      Xd[4] += Xpts[1] * Nb[0];
      Xd[5] += Xpts[2] * Nb[0];

      N++;
      Na++;
      Nb++;
      Xpts += 3;
    }
  }
}

/*
  Compute the postion as well ass the first and second derivatives
  along the surface of the shell.

  At the same time compute the shape functions, and their derivatives,
  at the current Gauss point gpt. Here Xpts is a vector of the nodal
  locations.

  input:
  gpt: the current Gauss point
  Xpts: the nodal locations for the

  output:
  X, Xd, Xdd: position, first and second derivatives of the surface
  N: the shape functions
  Na, Nb: the first derivatives of the shape functions
  Naa, Nab, Nbb: the second derivatives of the shape functions
*/
void shell_hessian(const int order, TacsScalar X[], TacsScalar Xd[],
                   TacsScalar Xdd[], double N[], double Na[], double Nb[],
                   double Naa[], double Nab[], double Nbb[], const double gpt[],
                   const TacsScalar Xpts[]) {
  double na[8], nb[8];
  double dna[8], dnb[8];
  double ddna[8], ddnb[8];

  X[0] = X[1] = X[2] = TacsScalar(0.0);

  Xd[0] = Xd[1] = Xd[2] = TacsScalar(0.0);
  Xd[3] = Xd[4] = Xd[5] = TacsScalar(0.0);

  Xdd[0] = Xdd[1] = Xdd[2] = TacsScalar(0.0);
  Xdd[3] = Xdd[4] = Xdd[5] = TacsScalar(0.0);
  Xdd[6] = Xdd[7] = Xdd[8] = TacsScalar(0.0);

  if (order <= 8) {
    if (order == 2) {
      FElibrary::lagrangeSF(na, dna, ddna, gpt[0], 2);
      FElibrary::lagrangeSF(nb, dnb, ddnb, gpt[1], 2);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[0] * nb[1];
      N[3] = na[1] * nb[1];

      // First derivatives
      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[0] * nb[1];
      Na[3] = dna[1] * nb[1];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[0] * dnb[1];
      Nb[3] = na[1] * dnb[1];

      // Second derivatives
      Naa[0] = Naa[1] = Naa[2] = Naa[3] = 0.0;

      Nab[0] = dna[0] * dnb[0];
      Nab[1] = dna[1] * dnb[0];
      Nab[2] = dna[0] * dnb[1];
      Nab[3] = dna[1] * dnb[1];

      Nbb[0] = Nbb[1] = Nbb[2] = Nbb[3] = 0.0;
    } else if (order == 3) {
      FElibrary::lagrangeSF(na, dna, ddna, gpt[0], 3);
      FElibrary::lagrangeSF(nb, dnb, ddnb, gpt[1], 3);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[2] * nb[0];
      N[3] = na[0] * nb[1];
      N[4] = na[1] * nb[1];
      N[5] = na[2] * nb[1];
      N[6] = na[0] * nb[2];
      N[7] = na[1] * nb[2];
      N[8] = na[2] * nb[2];

      // First derivatives
      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[2] * nb[0];
      Na[3] = dna[0] * nb[1];
      Na[4] = dna[1] * nb[1];
      Na[5] = dna[2] * nb[1];
      Na[6] = dna[0] * nb[2];
      Na[7] = dna[1] * nb[2];
      Na[8] = dna[2] * nb[2];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[2] * dnb[0];
      Nb[3] = na[0] * dnb[1];
      Nb[4] = na[1] * dnb[1];
      Nb[5] = na[2] * dnb[1];
      Nb[6] = na[0] * dnb[2];
      Nb[7] = na[1] * dnb[2];
      Nb[8] = na[2] * dnb[2];

      // Second derivatives
      Naa[0] = ddna[0] * nb[0];
      Naa[1] = ddna[1] * nb[0];
      Naa[2] = ddna[2] * nb[0];
      Naa[3] = ddna[0] * nb[1];
      Naa[4] = ddna[1] * nb[1];
      Naa[5] = ddna[2] * nb[1];
      Naa[6] = ddna[0] * nb[2];
      Naa[7] = ddna[1] * nb[2];
      Naa[8] = ddna[2] * nb[2];

      Nab[0] = dna[0] * dnb[0];
      Nab[1] = dna[1] * dnb[0];
      Nab[2] = dna[2] * dnb[0];
      Nab[3] = dna[0] * dnb[1];
      Nab[4] = dna[1] * dnb[1];
      Nab[5] = dna[2] * dnb[1];
      Nab[6] = dna[0] * dnb[2];
      Nab[7] = dna[1] * dnb[2];
      Nab[8] = dna[2] * dnb[2];

      Nbb[0] = na[0] * ddnb[0];
      Nbb[1] = na[1] * ddnb[0];
      Nbb[2] = na[2] * ddnb[0];
      Nbb[3] = na[0] * ddnb[1];
      Nbb[4] = na[1] * ddnb[1];
      Nbb[5] = na[2] * ddnb[1];
      Nbb[6] = na[0] * ddnb[2];
      Nbb[7] = na[1] * ddnb[2];
      Nbb[8] = na[2] * ddnb[2];
    } else if (order == 4) {
      FElibrary::lagrangeSF(na, dna, ddna, gpt[0], 4);
      FElibrary::lagrangeSF(nb, dnb, ddnb, gpt[1], 4);

      N[0] = na[0] * nb[0];
      N[1] = na[1] * nb[0];
      N[2] = na[2] * nb[0];
      N[3] = na[3] * nb[0];
      N[4] = na[0] * nb[1];
      N[5] = na[1] * nb[1];
      N[6] = na[2] * nb[1];
      N[7] = na[3] * nb[1];
      N[8] = na[0] * nb[2];
      N[9] = na[1] * nb[2];
      N[10] = na[2] * nb[2];
      N[11] = na[3] * nb[2];
      N[12] = na[0] * nb[3];
      N[13] = na[1] * nb[3];
      N[14] = na[2] * nb[3];
      N[15] = na[3] * nb[3];

      // First derivatives
      Na[0] = dna[0] * nb[0];
      Na[1] = dna[1] * nb[0];
      Na[2] = dna[2] * nb[0];
      Na[3] = dna[3] * nb[0];
      Na[4] = dna[0] * nb[1];
      Na[5] = dna[1] * nb[1];
      Na[6] = dna[2] * nb[1];
      Na[7] = dna[3] * nb[1];
      Na[8] = dna[0] * nb[2];
      Na[9] = dna[1] * nb[2];
      Na[10] = dna[2] * nb[2];
      Na[11] = dna[3] * nb[2];
      Na[12] = dna[0] * nb[3];
      Na[13] = dna[1] * nb[3];
      Na[14] = dna[2] * nb[3];
      Na[15] = dna[3] * nb[3];

      Nb[0] = na[0] * dnb[0];
      Nb[1] = na[1] * dnb[0];
      Nb[2] = na[2] * dnb[0];
      Nb[3] = na[3] * dnb[0];
      Nb[4] = na[0] * dnb[1];
      Nb[5] = na[1] * dnb[1];
      Nb[6] = na[2] * dnb[1];
      Nb[7] = na[3] * dnb[1];
      Nb[8] = na[0] * dnb[2];
      Nb[9] = na[1] * dnb[2];
      Nb[10] = na[2] * dnb[2];
      Nb[11] = na[3] * dnb[2];
      Nb[12] = na[0] * dnb[3];
      Nb[13] = na[1] * dnb[3];
      Nb[14] = na[2] * dnb[3];
      Nb[15] = na[3] * dnb[3];

      // Second derivatives
      Naa[0] = ddna[0] * nb[0];
      Naa[1] = ddna[1] * nb[0];
      Naa[2] = ddna[2] * nb[0];
      Naa[3] = ddna[3] * nb[0];
      Naa[4] = ddna[0] * nb[1];
      Naa[5] = ddna[1] * nb[1];
      Naa[6] = ddna[2] * nb[1];
      Naa[7] = ddna[3] * nb[1];
      Naa[8] = ddna[0] * nb[2];
      Naa[9] = ddna[1] * nb[2];
      Naa[10] = ddna[2] * nb[2];
      Naa[11] = ddna[3] * nb[2];
      Naa[12] = ddna[0] * nb[3];
      Naa[13] = ddna[1] * nb[3];
      Naa[14] = ddna[2] * nb[3];
      Naa[15] = ddna[3] * nb[3];

      Nab[0] = dna[0] * dnb[0];
      Nab[1] = dna[1] * dnb[0];
      Nab[2] = dna[2] * dnb[0];
      Nab[3] = dna[3] * dnb[0];
      Nab[4] = dna[0] * dnb[1];
      Nab[5] = dna[1] * dnb[1];
      Nab[6] = dna[2] * dnb[1];
      Nab[7] = dna[3] * dnb[1];
      Nab[8] = dna[0] * dnb[2];
      Nab[9] = dna[1] * dnb[2];
      Nab[10] = dna[2] * dnb[2];
      Nab[11] = dna[3] * dnb[2];
      Nab[12] = dna[0] * dnb[3];
      Nab[13] = dna[1] * dnb[3];
      Nab[14] = dna[2] * dnb[3];
      Nab[15] = dna[3] * dnb[3];

      Nbb[0] = na[0] * ddnb[0];
      Nbb[1] = na[1] * ddnb[0];
      Nbb[2] = na[2] * ddnb[0];
      Nbb[3] = na[3] * ddnb[0];
      Nbb[4] = na[0] * ddnb[1];
      Nbb[5] = na[1] * ddnb[1];
      Nbb[6] = na[2] * ddnb[1];
      Nbb[7] = na[3] * ddnb[1];
      Nbb[8] = na[0] * ddnb[2];
      Nbb[9] = na[1] * ddnb[2];
      Nbb[10] = na[2] * ddnb[2];
      Nbb[11] = na[3] * ddnb[2];
      Nbb[12] = na[0] * ddnb[3];
      Nbb[13] = na[1] * ddnb[3];
      Nbb[14] = na[2] * ddnb[3];
      Nbb[15] = na[3] * ddnb[3];
    } else {
      FElibrary::lagrangeSF(na, dna, ddna, gpt[0], order);
      FElibrary::lagrangeSF(nb, dnb, ddnb, gpt[1], order);

      double *Np = N;
      double *Nap = Na;
      double *Nbp = Nb;
      double *Naap = Naa;
      double *Nabp = Nab;
      double *Nbbp = Nbb;

      for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
          Np[0] = na[i] * nb[j];
          Nap[0] = dna[i] * nb[j];
          Nbp[0] = na[i] * dnb[j];

          Naap[0] = ddna[i] * nb[j];
          Nabp[0] = dna[i] * dnb[j];
          Nbbp[0] = na[i] * ddnb[j];

          Np++;
          Nap++;
          Nbp++;
          Naap++;
          Nabp++;
          Nbbp++;
        }
      }
    }

    const int order2 = order * order;
    for (int i = 0; i < order2; i++) {
      X[0] += Xpts[0] * N[0];
      X[1] += Xpts[1] * N[0];
      X[2] += Xpts[2] * N[0];

      // First derivatives
      Xd[0] += Xpts[0] * Na[0];
      Xd[1] += Xpts[1] * Na[0];
      Xd[2] += Xpts[2] * Na[0];

      Xd[3] += Xpts[0] * Nb[0];
      Xd[4] += Xpts[1] * Nb[0];
      Xd[5] += Xpts[2] * Nb[0];

      // Second derivatives
      Xdd[0] += Xpts[0] * Naa[0];
      Xdd[1] += Xpts[1] * Naa[0];
      Xdd[2] += Xpts[2] * Naa[0];

      Xdd[3] += Xpts[0] * Nab[0];
      Xdd[4] += Xpts[1] * Nab[0];
      Xdd[5] += Xpts[2] * Nab[0];

      Xdd[6] += Xpts[0] * Nbb[0];
      Xdd[7] += Xpts[1] * Nbb[0];
      Xdd[8] += Xpts[2] * Nbb[0];

      N++;
      Na++;
      Nb++;
      Naa++;
      Nab++;
      Nbb++;
      Xpts += 3;
    }
  }
}

TACS_END_NAMESPACE
