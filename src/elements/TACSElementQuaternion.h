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

#ifndef TACS_ELEMENT_QUATERNION_H
#define TACS_ELEMENT_QUATERNION_H

#include "TACSElementAlgebra.h"
#include "TACSObject.h"

/*
  Given the quaternion parameters, compute the rotation matrix.

  input:
  eta:    the quaternion scalar
  eps:    the quaternion 3-vector

  output:
  C:      the rotation matrix
*/
static inline void computeRotationMat(const TacsScalar eta,
                                      const TacsScalar eps[], TacsScalar C[]) {
  C[0] = 1.0 - 2.0 * (eps[1] * eps[1] + eps[2] * eps[2]);
  C[1] = 2.0 * (eps[0] * eps[1] + eta * eps[2]);
  C[2] = 2.0 * (eps[0] * eps[2] - eta * eps[1]);

  C[3] = 2.0 * (eps[1] * eps[0] - eta * eps[2]);
  C[4] = 1.0 - 2.0 * (eps[0] * eps[0] + eps[2] * eps[2]);
  C[5] = 2.0 * (eps[1] * eps[2] + eta * eps[0]);

  C[6] = 2.0 * (eps[2] * eps[0] + eta * eps[1]);
  C[7] = 2.0 * (eps[2] * eps[1] - eta * eps[0]);
  C[8] = 1.0 - 2.0 * (eps[0] * eps[0] + eps[1] * eps[1]);
}

/*
  Given the quaternion parameters, and their time-derivatives,
  compute the time derivative of the rotation matrix

  input:
  eta:    the quaternion scalar
  eps:    the quaternion 3-vector
  deta:   the time derivative of the quaternion scalar
  deps:   the time derivative of the quaternion 3-vector

  output:
  C:      the rotation matrix
*/
static inline void computeRotationMatDeriv(const TacsScalar eta,
                                           const TacsScalar eps[],
                                           const TacsScalar deta,
                                           const TacsScalar deps[],
                                           TacsScalar C[]) {
  C[0] = -4.0 * (eps[1] * deps[1] + eps[2] * deps[2]);
  C[1] = 2.0 *
         (eps[0] * deps[1] + deps[0] * eps[1] + eta * deps[2] + deta * eps[2]);
  C[2] = 2.0 *
         (eps[0] * deps[2] + deps[0] * eps[2] - eta * deps[1] - deta * eps[1]);

  C[3] = 2.0 *
         (eps[1] * deps[0] + deps[1] * eps[0] - eta * deps[2] - deta * eps[2]);
  C[4] = -4.0 * (eps[0] * deps[0] + eps[2] * deps[2]);
  C[5] = 2.0 *
         (eps[1] * deps[2] + deps[1] * eps[2] + eta * deps[0] + deta * eps[0]);

  C[6] = 2.0 *
         (eps[2] * deps[0] + deps[2] * eps[0] + eta * deps[1] + deta * eps[1]);
  C[7] = 2.0 *
         (eps[2] * deps[1] + deps[2] * eps[1] - eta * deps[0] - deta * eps[0]);
  C[8] = -4.0 * (eps[0] * deps[0] + eps[1] * deps[1]);
}

/*
  Compute the product of the 3x4 rotation rate matrix with the
  given components of the quaternion vector.

  y <- S(eta, eps)*x
  y <- -2*eps*xeta + 2*(eta*I - eps^{x})*xeps

  input:
  eta:   the quaternion scalar
  eps:   the quaternion vector
  xeta:  the x component of the scalar
  xeps:  the x components of the vector

  output:
  y:     3-vector containing the result
*/
static inline void computeSRateProduct(const TacsScalar eta,
                                       const TacsScalar eps[],
                                       const TacsScalar xeta,
                                       const TacsScalar xeps[],
                                       TacsScalar y[]) {
  crossProduct(-2.0, eps, xeps, y);
  vec3Axpy(2.0 * eta, xeps, y);
  vec3Axpy(-2.0 * xeta, eps, y);
}

/*
  Add the product of the 3x4 rotation rate matrix with the given
  components of the quaternion vector.

  y <- y + a*S(eta, eps)*x
  y <- y - 2*a*eps*xeta + 2*a*(eta*I - eps^{x})*xeps

  input:
  eta:   the quaternion scalar
  eps:   the quaternion vector
  xeta:  the x component of the scalar
  xeps:  the x components of the vector

  output:
  y:     3-vector containing the result
*/
static inline void addSRateProduct(const TacsScalar eta, const TacsScalar eps[],
                                   const TacsScalar xeta,
                                   const TacsScalar xeps[], TacsScalar y[]) {
  crossProductAdd(-2.0, eps, xeps, y);
  vec3Axpy(2.0 * eta, xeps, y);
  vec3Axpy(-2.0 * xeta, eps, y);
}

/*
  Add the product of the transpose of the 3x4 rotation rate matrix
  with the given components of x.

  y <- y + a*S(eta, eps)^{T}*x

  [ yeta ] += a*S^{T}*x = [      2*a*eps^{T}*x      ]
  [ yeps ] +=             [ 2*a*(eta*I + eps^{x})*x ]


  input:
  a:     the scalar input
  eta:   the quaternion scalar
  eps:   the quaternion vector
  x:     the 3-vector input

  output:
  yeps:  the scalar component of the output
  yeta:  the 3-vector component of the output
*/
static inline void addSRateTransProduct(const TacsScalar a,
                                        const TacsScalar eta,
                                        const TacsScalar eps[],
                                        const TacsScalar x[], TacsScalar *yeta,
                                        TacsScalar yeps[]) {
  *yeta -= 2.0 * a * vec3Dot(eps, x);
  crossProductAdd(2.0 * a, eps, x, yeps);
  vec3Axpy(2.0 * a * eta, x, yeps);
}

/*
  Compute the 3x4 rotation rate matrix that takes the quaternion rates
  and returns the angular velocity:

  S = 2[ -eps | (eta*I - eps^{x}) ]

  input:
  eta:   the quaternion scalar
  eps:   the quaternion vector

  output:
  S:     the 3x4 rate matrix
*/
static inline void computeSRateMat(const TacsScalar eta, const TacsScalar eps[],
                                   TacsScalar S[]) {
  S[0] = -2.0 * eps[0];
  S[1] = 2.0 * eta;
  S[2] = 2.0 * eps[2];
  S[3] = -2.0 * eps[1];

  S[4] = -2.0 * eps[1];
  S[5] = -2.0 * eps[2];
  S[6] = 2.0 * eta;
  S[7] = 2.0 * eps[0];

  S[8] = -2.0 * eps[2];
  S[9] = 2.0 * eps[1];
  S[10] = -2.0 * eps[0];
  S[11] = 2.0 * eta;
}

/*
  Compute the product of the matrix

  y <- y + a*D(v)^{T}*x

  The matrix D(v) is the derivative of d(C*v)/dq and is given as
  follows:

  D(v) = 2*[ v^{x}*eps | (v^{x}*eps^{x} + eta*v^{x} - 2*eps^{x}*v^{x}) ]


  D(v)^{T}*x = 2*[ -eps^{T}*v^{x}*x ]
  .              [ eps^{x}*v^{x}*x - eta*v^{x}*x - 2*v^{x}*eps^{x}*x ]

  input:
  a:    the scalar input
  v:    the vector input for D(v)
  x:    the multiplier vector
  eta:  the quaternion scalar
  eps:  the quaternion 3-vector

  output:
  yeta:  the quaternion output scalar
  eps:   the quaternion output 3-vector
*/
static inline void addDMatTransProduct(const TacsScalar a, const TacsScalar v[],
                                       const TacsScalar x[],
                                       const TacsScalar eta,
                                       const TacsScalar eps[], TacsScalar *yeta,
                                       TacsScalar yeps[]) {
  // Compute the cross product
  TacsScalar t[3];
  crossProduct(1.0, v, x, t);
  *yeta -= 2.0 * a * vec3Dot(eps, t);

  // Add the term 2*eps^{x}*v^{x}*x - 2*eta* v^{x}*x
  crossProductAdd(2.0 * a, eps, t, yeps);
  crossProductAdd(-2.0 * a * eta, v, x, yeps);

  // Compute the final term -4*v^{x}*eps^{x}*x
  crossProduct(1.0, eps, x, t);
  crossProductAdd(-4.0 * a, v, t, yeps);
}

/*
  Compute the product of the matrix

  y <- y + a*E(v)^{T}*x

  The matrix E(v) is the derivative of d(C^{T}*v)/dq and is given as
  follows:

  E(v) = 2*[ -v^{x}*eps | (v^{x}*eps^{x} - eta*v^{x} - 2*eps^{x}*v^{x}) ]

  input:
  a:    the scalar input
  v:    the vector input for D(v)
  x:    the multiplier vector
  eta:  the quaternion scalar
  eps:  the quaternion 3-vector

  output:
  yeta:  the quaternion output scalar
  eps:   the quaternion output 3-vector
*/
static inline void addEMatTransProduct(const TacsScalar a, const TacsScalar v[],
                                       const TacsScalar x[],
                                       const TacsScalar eta,
                                       const TacsScalar eps[], TacsScalar *yeta,
                                       TacsScalar yeps[]) {
  // Compute the cross product
  TacsScalar t[3];
  crossProduct(1.0, v, x, t);
  *yeta += 2.0 * a * vec3Dot(eps, t);

  // Add the term 2*eps^{x}*v^{x}*x + 2*eta*v^{x}*x
  crossProductAdd(2.0 * a, eps, t, yeps);
  crossProductAdd(2.0 * a * eta, v, x, yeps);

  // Compute the final term -4*v^{x}*eps^{x}*x
  crossProduct(1.0, eps, x, t);
  crossProductAdd(-4.0 * a, v, t, yeps);
}

/*
  Compute the elements of the D(v) matrix

  D(v) = 2*[ v^{x}*eps | (v^{x}*eps^{x} + eta*v^{x} - 2*eps^{x}*v^{x}) ]

  input:
  eta:   the quaternion scalar
  eps:   the quaternion 3-vector
  v:     the input vector

  output:
  D:     the 3x4 derivative matrix
*/
static inline void computeDMat(const TacsScalar eta, const TacsScalar eps[],
                               const TacsScalar v[], TacsScalar D[]) {
  D[0] = 2.0 * (v[1] * eps[2] - v[2] * eps[1]);
  D[1] = 2.0 * (v[1] * eps[1] + v[2] * eps[2]);
  D[2] = 2.0 * (eps[0] * v[1] - 2.0 * v[0] * eps[1] - eta * v[2]);
  D[3] = 2.0 * (eps[0] * v[2] - 2.0 * v[0] * eps[2] + eta * v[1]);

  D[4] = 2.0 * (v[2] * eps[0] - v[0] * eps[2]);
  D[5] = 2.0 * (eps[1] * v[0] - 2.0 * v[1] * eps[0] + eta * v[2]);
  D[6] = 2.0 * (v[0] * eps[0] + v[2] * eps[2]);
  D[7] = 2.0 * (eps[1] * v[2] - 2.0 * v[1] * eps[2] - eta * v[0]);

  D[8] = 2.0 * (v[0] * eps[1] - v[1] * eps[0]);
  D[9] = 2.0 * (eps[2] * v[0] - 2.0 * v[2] * eps[0] - eta * v[1]);
  D[10] = 2.0 * (eps[2] * v[1] - 2.0 * v[2] * eps[1] + eta * v[0]);
  D[11] = 2.0 * (v[0] * eps[0] + v[1] * eps[1]);
}

/*
  Add the elements of the D(v) matrix

  D(v) = 2*[ v^{x}*eps | (v^{x}*eps^{x} + eta*v^{x} - 2*eps^{x}*v^{x}) ]

  input:
  eta:   the quaternion scalar
  eps:   the quaternion 3-vector
  v:     the input vector

  output:
  D:     the 3x4 derivative matrix
*/
static inline void addBlockDMat(const TacsScalar a, const TacsScalar eta,
                                const TacsScalar eps[], const TacsScalar v[],
                                TacsScalar D[], const int ldd) {
  const TacsScalar d = 2.0 * a;
  D[0] += d * (v[1] * eps[2] - v[2] * eps[1]);
  D[1] += d * (v[1] * eps[1] + v[2] * eps[2]);
  D[2] += d * (eps[0] * v[1] - 2.0 * v[0] * eps[1] - eta * v[2]);
  D[3] += d * (eps[0] * v[2] - 2.0 * v[0] * eps[2] + eta * v[1]);
  D += ldd;

  D[0] += d * (v[2] * eps[0] - v[0] * eps[2]);
  D[1] += d * (eps[1] * v[0] - 2.0 * v[1] * eps[0] + eta * v[2]);
  D[2] += d * (v[0] * eps[0] + v[2] * eps[2]);
  D[3] += d * (eps[1] * v[2] - 2.0 * v[1] * eps[2] - eta * v[0]);
  D += ldd;

  D[0] += d * (v[0] * eps[1] - v[1] * eps[0]);
  D[1] += d * (eps[2] * v[0] - 2.0 * v[2] * eps[0] - eta * v[1]);
  D[2] += d * (eps[2] * v[1] - 2.0 * v[2] * eps[1] + eta * v[0]);
  D[3] += d * (v[0] * eps[0] + v[1] * eps[1]);
  D += ldd;
}

/*
  Compute the elements of the D(v) matrix

  D(v) = 2*[ v^{x}*eps | (v^{x}*eps^{x} + eta*v^{x} - 2*eps^{x}*v^{x}) ]

  input:
  eta:   the quaternion scalar
  eps:   the quaternion 3-vector
  v:     the input vector

  output:
  D:     the 3x4 derivative matrix
*/
static inline void addBlockDMatTrans(const TacsScalar a, const TacsScalar eta,
                                     const TacsScalar eps[],
                                     const TacsScalar v[], TacsScalar D[],
                                     const int ldd) {
  const TacsScalar d = 2.0 * a;
  D[0] += d * (v[1] * eps[2] - v[2] * eps[1]);
  D[1] += d * (v[2] * eps[0] - v[0] * eps[2]);
  D[2] += d * (v[0] * eps[1] - v[1] * eps[0]);
  D += ldd;

  D[0] += d * (v[1] * eps[1] + v[2] * eps[2]);
  D[1] += d * (eps[1] * v[0] - 2.0 * v[1] * eps[0] + eta * v[2]);
  D[2] += d * (eps[2] * v[0] - 2.0 * v[2] * eps[0] - eta * v[1]);
  D += ldd;

  D[0] += d * (eps[0] * v[1] - 2.0 * v[0] * eps[1] - eta * v[2]);
  D[1] += d * (v[0] * eps[0] + v[2] * eps[2]);
  D[2] += d * (eps[2] * v[1] - 2.0 * v[2] * eps[1] + eta * v[0]);
  D += ldd;

  D[0] += d * (eps[0] * v[2] - 2.0 * v[0] * eps[2] + eta * v[1]);
  D[1] += d * (eps[1] * v[2] - 2.0 * v[1] * eps[2] - eta * v[0]);
  D[2] += d * (v[0] * eps[0] + v[1] * eps[1]);
  D += ldd;
}

/*
  Compute the elements of the E(v) matrix

  E(v) = 2*[ -v^{x}*eps | (v^{x}*eps^{x} + eta*v^{x} + 2*eps^{x}*v^{x}) ]

  input:
  eta:   the quaternion scalar
  eps:   the quaternion 3-vector
  v:     the input vector

  output:
  E:     the 3x4 derivative matrix
*/
static inline void computeEMat(const TacsScalar eta, const TacsScalar eps[],
                               const TacsScalar v[], TacsScalar E[]) {
  E[0] = 2.0 * (v[2] * eps[1] - v[1] * eps[2]);
  E[1] = 2.0 * (v[1] * eps[1] + v[2] * eps[2]);
  E[2] = 2.0 * (eps[0] * v[1] - 2.0 * v[0] * eps[1] + eta * v[2]);
  E[3] = 2.0 * (eps[0] * v[2] - 2.0 * v[0] * eps[2] - eta * v[1]);

  E[4] = 2.0 * (v[0] * eps[2] - v[2] * eps[0]);
  E[5] = 2.0 * (eps[1] * v[0] - 2.0 * v[1] * eps[0] - eta * v[2]);
  E[6] = 2.0 * (v[0] * eps[0] + v[2] * eps[2]);
  E[7] = 2.0 * (eps[1] * v[2] - 2.0 * v[1] * eps[2] + eta * v[0]);

  E[8] = 2.0 * (v[1] * eps[0] - v[0] * eps[1]);
  E[9] = 2.0 * (eps[2] * v[0] - 2.0 * v[2] * eps[0] + eta * v[1]);
  E[10] = 2.0 * (eps[2] * v[1] - 2.0 * v[2] * eps[1] - eta * v[0]);
  E[11] = 2.0 * (v[0] * eps[0] + v[1] * eps[1]);
}

/*
  Add the derivative of the matrix d(D(v)^{T}*w)/dq to the provided
  block matrix. The derivative takes the form:

  TD(v,w) = 2*[ 0           (w^{x}*v)^{T}                     ]
  .           [ w^{x}*v     v*w^{T} + w*v^{T} - 2*(v^{T}*w)*I ]

  input:
  a:     scalar factor
  v:     the input 3-vector
  w:     the second input 3-vector
  ldd:   the leading row dimension of the block matrix D

  output:
  D:     the resulting block matrix
*/
static inline void addBlockDMatTransDeriv(const TacsScalar a,
                                          const TacsScalar v[],
                                          const TacsScalar w[], TacsScalar D[],
                                          const int ldd) {
  const TacsScalar b = 2.0 * a;
  // Compute the cross product w^{x}*v
  TacsScalar t[3];
  crossProduct(1.0, w, v, t);

  D[1] += b * t[0];
  D[2] += b * t[1];
  D[3] += b * t[2];
  D += ldd;

  D[0] += b * t[0];
  D[1] -= 2.0 * b * (v[1] * w[1] + v[2] * w[2]);
  D[2] += b * (v[0] * w[1] + w[0] * v[1]);
  D[3] += b * (v[0] * w[2] + w[0] * v[2]);
  D += ldd;

  D[0] += b * t[1];
  D[1] += b * (v[1] * w[0] + w[1] * v[0]);
  D[2] -= 2.0 * b * (v[0] * w[0] + v[2] * w[2]);
  D[3] += b * (v[1] * w[2] + w[1] * v[2]);
  D += ldd;

  D[0] += b * t[2];
  D[1] += b * (v[2] * w[0] + w[2] * v[0]);
  D[2] += b * (v[2] * w[1] + w[2] * v[1]);
  D[3] -= 2.0 * b * (v[0] * w[0] + v[1] * w[1]);
}

/*
  Add the 3x4 matrix E(v) to the block matrix

  E(v) = 2*[ -v^{x}*eps | (eps*v^{T} - 2*v*eps^{T} - eta*^{x} + v^{T}*eps*1) ]

  input:
  a:    the scalar value
  eta:   the scalar quaternion
  eps:   the 3-vector quaternion components
  v:     the 3-vector v
  ldd:   the leading dimension (row-dimension) of D

  output:
  D:     the result is added to this matrix
*/
static inline void addBlockEMat(const TacsScalar a, const TacsScalar eta,
                                const TacsScalar eps[], const TacsScalar v[],
                                TacsScalar D[], const int ldd) {
  const TacsScalar b = 2.0 * a;
  D[0] += b * (v[2] * eps[1] - v[1] * eps[2]);
  D[1] += b * (v[1] * eps[1] + v[2] * eps[2]);
  D[2] += b * (eps[0] * v[1] - 2.0 * v[0] * eps[1] + eta * v[2]);
  D[3] += b * (eps[0] * v[2] - 2.0 * v[0] * eps[2] - eta * v[1]);
  D += ldd;

  D[0] += b * (v[0] * eps[2] - v[2] * eps[0]);
  D[1] += b * (eps[1] * v[0] - 2.0 * v[1] * eps[0] - eta * v[2]);
  D[2] += b * (v[0] * eps[0] + v[2] * eps[2]);
  D[3] += b * (eps[1] * v[2] - 2.0 * v[1] * eps[2] + eta * v[0]);
  D += ldd;

  D[0] += b * (v[1] * eps[0] - v[0] * eps[1]);
  D[1] += b * (eps[2] * v[0] - 2.0 * v[2] * eps[0] + eta * v[1]);
  D[2] += b * (eps[2] * v[1] - 2.0 * v[2] * eps[1] - eta * v[0]);
  D[3] += b * (v[0] * eps[0] + v[1] * eps[1]);
}

/*
  Add the transpose of the 3x4 matrix E(v) to the block matrix

  E(v) = 2*[ -v^{x}*eps | (eps*v^{T} - 2*v*eps^{T} - eta*^{x} + v^{T}*eps*1) ]

  input:
  a:    the scalar value
  eta:   the scalar quaternion
  eps:   the 3-vector quaternion components
  v:     the 3-vector v
  ldd:   the leading dimension (row-dimension) of D

  output:
  D:     the result is added to this matrix
*/
static inline void addBlockEMatTrans(const TacsScalar a, const TacsScalar eta,
                                     const TacsScalar eps[],
                                     const TacsScalar v[], TacsScalar D[],
                                     const int ldd) {
  const TacsScalar b = 2.0 * a;
  D[0] += b * (v[2] * eps[1] - v[1] * eps[2]);
  D[1] += b * (v[0] * eps[2] - v[2] * eps[0]);
  D[2] += b * (v[1] * eps[0] - v[0] * eps[1]);
  D += ldd;

  D[0] += b * (v[1] * eps[1] + v[2] * eps[2]);
  D[1] += b * (eps[1] * v[0] - 2.0 * v[1] * eps[0] - eta * v[2]);
  D[2] += b * (eps[2] * v[0] - 2.0 * v[2] * eps[0] + eta * v[1]);
  D += ldd;

  D[0] += b * (eps[0] * v[1] - 2.0 * v[0] * eps[1] + eta * v[2]);
  D[1] += b * (v[0] * eps[0] + v[2] * eps[2]);
  D[2] += b * (eps[2] * v[1] - 2.0 * v[2] * eps[1] - eta * v[0]);
  D += ldd;

  D[0] += b * (eps[0] * v[2] - 2.0 * v[0] * eps[2] - eta * v[1]);
  D[1] += b * (eps[1] * v[2] - 2.0 * v[1] * eps[2] + eta * v[0]);
  D[2] += b * (v[0] * eps[0] + v[1] * eps[1]);
}

/*
  Add the 4x4 matrix from the derivative of the transpose of the
  angular rate S matrix

  d(S^{T}*v)/dq =
  [ 0 | -v^{T} ]
  [ v | -v^{x} ]

  input:
  a:      a scalar multiplier on the contribution
  v:      the input vector
  ldd:    the leading dimension of the Jacobian matrix

  output:
  D:      the Jacobian matrix to which the the contribution is added
*/
static inline void addSRateMatTransDeriv(const TacsScalar a,
                                         const TacsScalar v[], TacsScalar D[],
                                         const int ldd) {
  const TacsScalar b = 2.0 * a;
  D[1] -= b * v[0];
  D[2] -= b * v[1];
  D[3] -= b * v[2];
  D += ldd;

  D[0] += b * v[0];
  D[2] += b * v[2];
  D[3] -= b * v[1];
  D += ldd;

  D[0] += b * v[1];
  D[1] -= b * v[2];
  D[3] += b * v[0];
  D += ldd;

  D[0] += b * v[2];
  D[1] += b * v[1];
  D[2] -= b * v[0];
}

/*
  Compute the product of a 3x3 and 3x4 matrix and add the result to
  the block matrix D

  D += A*B

  input:
  A:    a 3x3 matrix in row-major order
  B:    a 3x4 matrix in row-major order

  output:
  D:    the result is added to this matrix
*/
static inline void addBlock3x3x4Product(const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar D[],
                                        const int ldd) {
  D[0] += A[0] * B[0] + A[1] * B[4] + A[2] * B[8];
  D[1] += A[0] * B[1] + A[1] * B[5] + A[2] * B[9];
  D[2] += A[0] * B[2] + A[1] * B[6] + A[2] * B[10];
  D[3] += A[0] * B[3] + A[1] * B[7] + A[2] * B[11];
  D += ldd;

  D[0] += A[3] * B[0] + A[4] * B[4] + A[5] * B[8];
  D[1] += A[3] * B[1] + A[4] * B[5] + A[5] * B[9];
  D[2] += A[3] * B[2] + A[4] * B[6] + A[5] * B[10];
  D[3] += A[3] * B[3] + A[4] * B[7] + A[5] * B[11];
  D += ldd;

  D[0] += A[6] * B[0] + A[7] * B[4] + A[8] * B[8];
  D[1] += A[6] * B[1] + A[7] * B[5] + A[8] * B[9];
  D[2] += A[6] * B[2] + A[7] * B[6] + A[8] * B[10];
  D[3] += A[6] * B[3] + A[7] * B[7] + A[8] * B[11];
}

/*
  Compute the product of a 3x3 and 3x4 matrix and add the result to
  the block matrix D

  D += A^{T}*B

  input:
  A:    a 3x4 matrix in row-major order
  B:    a 3x4 matrix in row-major order

  output:
  D:    the result is added to this matrix
*/
static inline void addBlock4x3x3Product(const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar D[],
                                        const int ldd) {
  D[0] += A[0] * B[0] + A[4] * B[3] + A[8] * B[6];
  D[1] += A[0] * B[1] + A[4] * B[4] + A[8] * B[7];
  D[2] += A[0] * B[2] + A[4] * B[5] + A[8] * B[8];
  D += ldd;

  D[0] += A[1] * B[0] + A[5] * B[3] + A[9] * B[6];
  D[1] += A[1] * B[1] + A[5] * B[4] + A[9] * B[7];
  D[2] += A[1] * B[2] + A[5] * B[5] + A[9] * B[8];
  D += ldd;

  D[0] += A[2] * B[0] + A[6] * B[3] + A[10] * B[6];
  D[1] += A[2] * B[1] + A[6] * B[4] + A[10] * B[7];
  D[2] += A[2] * B[2] + A[6] * B[5] + A[10] * B[8];
  D += ldd;

  D[0] += A[3] * B[0] + A[7] * B[3] + A[11] * B[6];
  D[1] += A[3] * B[1] + A[7] * B[4] + A[11] * B[7];
  D[2] += A[3] * B[2] + A[7] * B[5] + A[11] * B[8];
}

/*
  Compute: D += a*A^{T}*B where a is a scalar, and A and B are 3x4
  matrices stored in column-major order.

  input:
  a:    the scalar multiple
  A:    3x4 matrix in row-major order
  B:    3x4 matrix in row-major order
  ldd:  the leading row dimension of the Jacobian matrix D

  output:
  D:    the result is added to this matrix D += a*A^{T}*B
*/
static inline void addBlock3x4Product(const TacsScalar a, const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar D[],
                                      const int ldd) {
  D[0] += a * (A[0] * B[0] + A[4] * B[4] + A[8] * B[8]);
  D[1] += a * (A[0] * B[1] + A[4] * B[5] + A[8] * B[9]);
  D[2] += a * (A[0] * B[2] + A[4] * B[6] + A[8] * B[10]);
  D[3] += a * (A[0] * B[3] + A[4] * B[7] + A[8] * B[11]);
  D += ldd;

  D[0] += a * (A[1] * B[0] + A[5] * B[4] + A[9] * B[8]);
  D[1] += a * (A[1] * B[1] + A[5] * B[5] + A[9] * B[9]);
  D[2] += a * (A[1] * B[2] + A[5] * B[6] + A[9] * B[10]);
  D[3] += a * (A[1] * B[3] + A[5] * B[7] + A[9] * B[11]);
  D += ldd;

  D[0] += a * (A[2] * B[0] + A[6] * B[4] + A[10] * B[8]);
  D[1] += a * (A[2] * B[1] + A[6] * B[5] + A[10] * B[9]);
  D[2] += a * (A[2] * B[2] + A[6] * B[6] + A[10] * B[10]);
  D[3] += a * (A[2] * B[3] + A[6] * B[7] + A[10] * B[11]);
  D += ldd;

  D[0] += a * (A[3] * B[0] + A[7] * B[4] + A[11] * B[8]);
  D[1] += a * (A[3] * B[1] + A[7] * B[5] + A[11] * B[9]);
  D[2] += a * (A[3] * B[2] + A[7] * B[6] + A[11] * B[10]);
  D[3] += a * (A[3] * B[3] + A[7] * B[7] + A[11] * B[11]);
}

/*
  Compute the second derivative of the product of the transpose of a
  quaternion-parametrized rotation matrix with a vector, i.e.:

  d^2/dq^2(C^{T}*v)

  where v is a constant vector. Note that the second derivative of the
  rotation matrix is a constant and so this code only depends on the
  input vector v.

  The order of the derivatives is as follows:
  d^2(C^{T}*v)/(d(eta)d(epsilon_{i})

  d^2(C^{T}*v)/(d(epsilon_{i})d(epsilon_{j})

  for (i,j) = (0,0), (0,1), (0,2), (1,1), (1,2), (2,2)

  The result is stored in a 9x3 array in row-major order.

  input:
  v:    the constant vector in the multiplication C^{T}*v

  output:
  dv    the second derivatives of C^{T}*v w.r.t q
*/
static inline void computeQtr2ndDeriv(const TacsScalar v[], TacsScalar dv[]) {
  // Derivatives of eta and eps
  dv[0] = 0.0;
  dv[1] = -2.0 * v[2];
  dv[2] = 2.0 * v[1];
  dv += 3;

  dv[0] = 2.0 * v[2];
  dv[1] = 0.0;
  dv[2] = -2.0 * v[0];
  dv += 3;

  dv[0] = -2.0 * v[1];
  dv[1] = 2.0 * v[0];
  dv[2] = 0.0;
  dv += 3;

  // Second derivatives w.r.t eps
  // C,11
  dv[0] = 0.0;
  dv[1] = -4.0 * v[1];
  dv[2] = -4.0 * v[2];
  dv += 3;

  // C,12
  dv[0] = 2.0 * v[1];
  dv[1] = 2.0 * v[0];
  dv[2] = 0.0;
  dv += 3;

  // C,13
  dv[0] = 2.0 * v[2];
  dv[1] = 0.0;
  dv[2] = 2.0 * v[0];
  dv += 3;

  // C,22
  dv[0] = -4.0 * v[0];
  dv[1] = 0.0;
  dv[2] = -4.0 * v[2];
  dv += 3;

  // C,23
  dv[0] = 0.0;
  dv[1] = 2.0 * v[2];
  dv[2] = 2.0 * v[1];
  dv += 3;

  // C,33
  dv[0] = -4.0 * v[0];
  dv[1] = -4.0 * v[1];
  dv[2] = 0.0;
  dv += 3;
}

#endif  // TACS_ELEMENT_QUATERNION_H
