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

#ifndef TACS_ELEMENT_ALGEBRA_H
#define TACS_ELEMENT_ALGEBRA_H

/*
  A header file with lots of useful linear algebra. Note that this is
  designed to be included in .c/.cpp files directly.
*/

#include <math.h>

#include "TACSObject.h"

/*
  Compute the cross-product

  out = (x cross y)

  input:
  x:    the first input 3-vector
  y:    the second input 3-vector

  output:
  out:  the resulting vector
*/
static inline void crossProduct(const TacsScalar x[], const TacsScalar y[],
                                TacsScalar out[]) {
  out[0] = (x[1] * y[2] - x[2] * y[1]);
  out[1] = (x[2] * y[0] - x[0] * y[2]);
  out[2] = (x[0] * y[1] - x[1] * y[0]);
}

/*
  Find the cross product and its sensitivity with respect to two
  3D input vectors and their sensitivity

  input:
  A:    the first input 3-vector
  B:    the second input 3-vector
  sA:   sensitivity of the first input 3-vector
  sB:   sensitivity of the second input 3-vector

  output:
  out:  the resulting vector
  sout: sensitivity of the resulting vector
*/
static inline void crossProductSens(const TacsScalar A[], const TacsScalar B[],
                                    const TacsScalar sA[],
                                    const TacsScalar sB[], TacsScalar out[],
                                    TacsScalar sout[]) {
  //   (A2 * B3 - A3 * B2 ), ( A3 *B1 - A1 * B3 ), ( A1 * B2 - A2 * B1 )

  out[0] = A[1] * B[2] - A[2] * B[1];
  out[1] = A[2] * B[0] - A[0] * B[2];
  out[2] = A[0] * B[1] - A[1] * B[0];

  sout[0] = (sA[1] * B[2] - sA[2] * B[1] + A[1] * sB[2] - A[2] * sB[1]);
  sout[1] = (sA[2] * B[0] - sA[0] * B[2] + A[2] * sB[0] - A[0] * sB[2]);
  sout[2] = (sA[0] * B[1] - sA[1] * B[0] + A[0] * sB[1] - A[1] * sB[0]);
}

/*
  Compute the cross-product

  out = a*(x cross y)

  input:
  a:    the scalar multiplier
  x:    the first input 3-vector
  y:    the second input 3-vector

  output:
  out:  the resulting vector
*/
static inline void crossProduct(const TacsScalar a, const TacsScalar x[],
                                const TacsScalar y[], TacsScalar out[]) {
  out[0] = a * (x[1] * y[2] - x[2] * y[1]);
  out[1] = a * (x[2] * y[0] - x[0] * y[2]);
  out[2] = a * (x[0] * y[1] - x[1] * y[0]);
}

/*
  Compute the cross-product and add it to the output vector

  out = a*(x cross y)

  input:
  a:    the scalar multiplier
  x:    the first input 3-vector
  y:    the second input 3-vector

  output:
  out:  the resulting vector
*/
static inline void crossProductAdd(const TacsScalar a, const TacsScalar x[],
                                   const TacsScalar y[], TacsScalar out[]) {
  out[0] += a * (x[1] * y[2] - x[2] * y[1]);
  out[1] += a * (x[2] * y[0] - x[0] * y[2]);
  out[2] += a * (x[0] * y[1] - x[1] * y[0]);
}

/*
  Scale the vector by the given scalar

  input
  a:   the scalar
  x:   the vector
*/
static inline void vec3Scale(const TacsScalar a, TacsScalar x[]) {
  x[0] *= a;
  x[1] *= a;
  x[2] *= a;
}

/*
  Compute the dot-product of two vectors

  return = x^{T}*y

  input:
  x:   the first vector
  y:   the second vector

  returns: the dot product
*/
static inline TacsScalar vec3Dot(const TacsScalar x[], const TacsScalar y[]) {
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

/*
  Add the product using an AXPY operation

  y <- y + alpha*x

  input:
  a:    the alpha scalar
  x:    the input vector

  in/out:
  y:    the result
*/
static inline void vec3Axpy(const TacsScalar a, const TacsScalar x[],
                            TacsScalar y[]) {
  y[0] += a * x[0];
  y[1] += a * x[1];
  y[2] += a * x[2];
}

/*
  Scale the vector by the given scalar

  input
  a:   the scalar
  x:   the vector
*/
static inline void vec2Scale(const TacsScalar a, TacsScalar x[]) {
  x[0] *= a;
  x[1] *= a;
}

/*
  Compute the dot-product of two vectors

  return = x^{T}*y

  input:
  x:   the first vector
  y:   the second vector

  returns: the dot product
*/
static inline TacsScalar vec2Dot(const TacsScalar x[], const TacsScalar y[]) {
  return (x[0] * y[0] + x[1] * y[1]);
}

/*
  Add the product using an AXPY operation

  y <- y + alpha*x

  input:
  a:    the alpha scalar
  x:    the input vector

  in/out:
  y:    the result
*/
static inline void vec2Axpy(const TacsScalar a, const TacsScalar x[],
                            TacsScalar y[]) {
  y[0] += a * x[0];
  y[1] += a * x[1];
}

/*
  Compute the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec3x3Outer(const TacsReal a[], const TacsReal b[],
                               TacsReal C[]) {
  C[0] = a[0] * b[0];
  C[1] = a[0] * b[1];
  C[2] = a[0] * b[2];

  C[3] = a[1] * b[0];
  C[4] = a[1] * b[1];
  C[5] = a[1] * b[2];

  C[6] = a[2] * b[0];
  C[7] = a[2] * b[1];
  C[8] = a[2] * b[2];
}

/*
  Add the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec3x3OuterAdd(const TacsScalar alpha, const TacsScalar a[],
                                  const TacsScalar b[], TacsScalar C[]) {
  C[0] += alpha * a[0] * b[0];
  C[1] += alpha * a[0] * b[1];
  C[2] += alpha * a[0] * b[2];

  C[3] += alpha * a[1] * b[0];
  C[4] += alpha * a[1] * b[1];
  C[5] += alpha * a[1] * b[2];

  C[6] += alpha * a[2] * b[0];
  C[7] += alpha * a[2] * b[1];
  C[8] += alpha * a[2] * b[2];
}

/*
  Compute the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec3x3Outer(const TacsComplex a[], const TacsReal b[],
                               TacsComplex C[]) {
  C[0] = a[0] * b[0];
  C[1] = a[0] * b[1];
  C[2] = a[0] * b[2];

  C[3] = a[1] * b[0];
  C[4] = a[1] * b[1];
  C[5] = a[1] * b[2];

  C[6] = a[2] * b[0];
  C[7] = a[2] * b[1];
  C[8] = a[2] * b[2];
}

/*
  Add the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec3x3OuterAdd(const TacsComplex alpha,
                                  const TacsComplex a[], const TacsReal b[],
                                  TacsComplex C[]) {
  C[0] += alpha * a[0] * b[0];
  C[1] += alpha * a[0] * b[1];
  C[2] += alpha * a[0] * b[2];

  C[3] += alpha * a[1] * b[0];
  C[4] += alpha * a[1] * b[1];
  C[5] += alpha * a[1] * b[2];

  C[6] += alpha * a[2] * b[0];
  C[7] += alpha * a[2] * b[1];
  C[8] += alpha * a[2] * b[2];
}

/*
  Compute the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec2x2Outer(const TacsScalar a[], const TacsScalar b[],
                               TacsScalar C[]) {
  C[0] = a[0] * b[0];
  C[1] = a[0] * b[1];

  C[2] = a[1] * b[0];
  C[3] = a[1] * b[1];
}

/*
  Add the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec2x2OuterAdd(const TacsScalar alpha, const TacsScalar a[],
                                  const TacsScalar b[], TacsScalar C[]) {
  C[0] += alpha * a[0] * b[0];
  C[1] += alpha * a[0] * b[1];

  C[2] += alpha * a[1] * b[0];
  C[3] += alpha * a[1] * b[1];
}

/*
  Compute the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec2x2Outer(const TacsComplex a[], const TacsReal b[],
                               TacsComplex C[]) {
  C[0] = a[0] * b[0];
  C[1] = a[0] * b[1];

  C[2] = a[1] * b[0];
  C[3] = a[1] * b[1];
}

/*
  Add the outer product of two vectors

  C <- a*b^{T}

  input:
  a:   the first vector
  b:   the second vector

  output:
  C:   the resulting matrix
*/
static inline void vec2x2OuterAdd(const TacsComplex alpha,
                                  const TacsComplex a[], const TacsReal b[],
                                  TacsComplex C[]) {
  C[0] += alpha * a[0] * b[0];
  C[1] += alpha * a[0] * b[1];

  C[2] += alpha * a[1] * b[0];
  C[3] += alpha * a[1] * b[1];
}

/*
  Compute the derivative of x/||x||_{2} w.r.t. x using x and the norm
  of x.

  This code computes a 3x3 matrix that takes the form:

  d(x/||x||_{2})/dx = (I*||x||^2 + x*x^{T})/||x||^3
*/
static inline void vec3NormDeriv(TacsScalar nrm, const TacsScalar x[],
                                 TacsScalar D[]) {
  TacsScalar s = 1.0 / (nrm * nrm * nrm);
  TacsScalar t = nrm * nrm;

  D[0] = s * (t - x[0] * x[0]);
  D[1] = -s * x[0] * x[1];
  D[2] = -s * x[0] * x[2];

  D[3] = -s * x[1] * x[0];
  D[4] = s * (t - x[1] * x[1]);
  D[5] = -s * x[1] * x[2];

  D[6] = -s * x[2] * x[0];
  D[7] = -s * x[2] * x[1];
  D[8] = s * (t - x[2] * x[2]);
}

/*
  Normalize a 3D vector in place

  input:
  A: the input 3-vector to be normalized

  output:
  A: normalized 3-vector
  Anrm: norm of input vector
*/
static inline TacsScalar vec3Normalize(TacsScalar A[]) {
  TacsScalar Anrm = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  TacsScalar invAnrm = 1.0 / Anrm;

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;
    A[2] = A[2] * invAnrm;
  }

  return Anrm;
}

/*
  Find the sensitivity of the normalized 3D vector

  input:
  A: the input 3-vector to be normalized

  output:
  A: normalized 3-vector
  sA: sensitivity ot the normalized 3-vector
  Anrm: norm of input vector
  sAnrm: sensitivity of the norm of input vector
*/
static inline TacsScalar vec3NormalizeSens(TacsScalar A[], TacsScalar *sAnrm,
                                           TacsScalar sA[]) {
  TacsScalar Anrm = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
  TacsScalar invAnrm = 1.0 / Anrm;
  *sAnrm = (1.0 / Anrm) * (A[0] * sA[0] + A[1] * sA[1] + A[2] * sA[2]);

  if (Anrm != 0.0) {
    A[0] = A[0] * invAnrm;
    A[1] = A[1] * invAnrm;
    A[2] = A[2] * invAnrm;

    sA[0] = (sA[0] - (A[0]) * (*sAnrm)) * invAnrm;
    sA[1] = (sA[1] - (A[1]) * (*sAnrm)) * invAnrm;
    sA[2] = (sA[2] - (A[2]) * (*sAnrm)) * invAnrm;
  } else {
    sA[0] = sA[1] = sA[2] = 0.0;
  }

  return Anrm;
}

/*
  Compute y <- A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void mat3x3Mult(const TacsScalar A[], const TacsScalar x[],
                              TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

/*
  Compute y <- A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void mat3x3Mult(const TacsComplex A[], const TacsReal x[],
                              TacsComplex y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  y[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

/*
  Compute y <- A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  output:
  y:   the resulting vector
*/
static inline void mat2x2Mult(const TacsScalar A[], const TacsScalar x[],
                              TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1];
  y[1] = A[2] * x[0] + A[3] * x[1];
}

/*
  Compute y <- A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  output:
  y:   the resulting vector
*/
static inline void mat2x2Mult(const TacsComplex A[], const TacsReal x[],
                              TacsComplex y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1];
  y[1] = A[2] * x[0] + A[3] * x[1];
}

/*
  Compute y <- A^{T}*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void mat3x3MultTrans(const TacsScalar A[], const TacsScalar x[],
                                   TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  y[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  y[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

/*
  Compute y <- A^{T}*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void mat3x3MultTrans(const TacsComplex A[], const TacsReal x[],
                                   TacsComplex y[]) {
  y[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  y[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  y[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

/*
  Compute y <- A^{T}*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  output:
  y:   the resulting vector
*/
static inline void mat2x2MultTrans(const TacsScalar A[], const TacsScalar x[],
                                   TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[2] * x[1];
  y[1] = A[1] * x[0] + A[3] * x[1];
}

/*
  Compute y <- A^{T}*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  output:
  y:   the resulting vector
*/
static inline void mat2x2MultTrans(const TacsComplex A[], const TacsReal x[],
                                   TacsComplex y[]) {
  y[0] = A[0] * x[0] + A[2] * x[1];
  y[1] = A[1] * x[0] + A[3] * x[1];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat3x3MultAdd(const TacsScalar A[], const TacsScalar x[],
                                 TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] += A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  y[2] += A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat2x2MultAdd(const TacsScalar A[], const TacsScalar x[],
                                 TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[1] * x[1];
  y[1] += A[2] * x[0] + A[3] * x[1];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat3x3MultTransAdd(const TacsScalar A[],
                                      const TacsScalar x[], TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  y[1] += A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  y[2] += A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat2x2MultTransAdd(const TacsScalar A[],
                                      const TacsScalar x[], TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[2] * x[1];
  y[1] += A[1] * x[0] + A[3] * x[1];
}

/*
  Compute the inner product with a 3x3 matrix:

  return:  x^{T}*A*y

  input:
  A:   a 3x3 matrix in row-major order
  x:   a 3-vector
  y:   a 3-vector
*/
static inline TacsScalar mat3x3Inner(const TacsScalar A[], const TacsScalar x[],
                                     const TacsScalar y[]) {
  return (x[0] * (A[0] * y[0] + A[1] * y[1] + A[2] * y[2]) +
          x[1] * (A[3] * y[0] + A[4] * y[1] + A[5] * y[2]) +
          x[2] * (A[6] * y[0] + A[7] * y[1] + A[8] * y[2]));
}

/*
  Compute the inner product with a 2x2 matrix:

  return:  x^{T}*A*y

  input:
  A:   a 2x2 matrix in row-major order
  x:   a 2-vector
  y:   a 2-vector
*/
static inline TacsScalar mat2x2Inner(const TacsScalar A[], const TacsScalar x[],
                                     const TacsScalar y[]) {
  return (x[0] * (A[0] * y[0] + A[1] * y[1]) +
          x[1] * (A[2] * y[0] + A[3] * y[1]));
}

/*
  Compute y^{T}*A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector
  y:   the resulting vector

  returns:  the inner product
*/
static inline TacsScalar mat3x3SymmInner(const TacsScalar A[],
                                         const TacsScalar x[],
                                         const TacsScalar y[]) {
  return (y[0] * (A[0] * x[0] + A[1] * x[1] + A[2] * x[2]) +
          y[1] * (A[1] * x[0] + A[3] * x[1] + A[4] * x[2]) +
          y[2] * (A[2] * x[0] + A[4] * x[1] + A[5] * x[2]));
}

/*
  Compute y^{T}*A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector
  y:   the resulting vector

  returns:  the inner product
*/
static inline TacsScalar mat2x2SymmInner(const TacsScalar A[],
                                         const TacsScalar x[],
                                         const TacsScalar y[]) {
  return (y[0] * (A[0] * x[0] + A[1] * x[1]) +
          y[1] * (A[1] * x[0] + A[2] * x[1]));
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat3x3SymmMult(const TacsScalar A[], const TacsScalar x[],
                                  TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] = A[1] * x[0] + A[3] * x[1] + A[4] * x[2];
  y[2] = A[2] * x[0] + A[4] * x[1] + A[5] * x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat2x2SymmMult(const TacsScalar A[], const TacsScalar x[],
                                  TacsScalar y[]) {
  y[0] = A[0] * x[0] + A[1] * x[1];
  y[1] = A[1] * x[0] + A[2] * x[1];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat3x3SymmMultAdd(const TacsScalar A[], const TacsScalar x[],
                                     TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  y[1] += A[1] * x[0] + A[3] * x[1] + A[4] * x[2];
  y[2] += A[2] * x[0] + A[4] * x[1] + A[5] * x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 2x2 input matrix in row-major order
  x:   the input 2-vector

  in/out:
  y:   the resulting vector
*/
static inline void mat2x2SymmMultAdd(const TacsScalar A[], const TacsScalar x[],
                                     TacsScalar y[]) {
  y[0] += A[0] * x[0] + A[1] * x[1];
  y[1] += A[1] * x[0] + A[2] * x[1];
}

/*
  Compute C = A*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3MatMult(const TacsScalar A[], const TacsScalar B[],
                                 TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
  C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
  C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];

  C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
  C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
  C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];

  C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
  C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
  C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

/*
  Compute C = A*B

  input:
  A:   the first 2x2 input matrix in row-major order
  B:   the second 2x2 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat2x2MatMult(const TacsScalar A[], const TacsScalar B[],
                                 TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[2];
  C[2] = A[2] * B[0] + A[3] * B[2];

  C[1] = A[0] * B[1] + A[1] * B[3];
  C[3] = A[2] * B[1] + A[3] * B[3];
}

/*
  Compute C = A*B^{T}

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3MatTransMult(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
  C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
  C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];

  C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
  C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
  C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];

  C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
  C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
  C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
}

/*
  Compute C = A*B^{T}

  input:
  A:   the first 2x2 input matrix in row-major order
  B:   the second 2x2 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat2x2MatTransMult(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[1];
  C[1] = A[0] * B[2] + A[1] * B[3];

  C[2] = A[2] * B[0] + A[3] * B[1];
  C[3] = A[2] * B[2] + A[3] * B[3];
}

/*
  Compute C = A^{T}*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3TransMatMult(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
  C[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
  C[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];

  C[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
  C[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
  C[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];

  C[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
  C[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
  C[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
}

/*
  Compute C = A^{T}*B

  input:
  A:   the first 2x2 input matrix in row-major order
  B:   the second 2x2 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat2x2TransMatMult(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[2] * B[2];
  C[1] = A[0] * B[1] + A[2] * B[3];

  C[2] = A[1] * B[0] + A[3] * B[2];
  C[3] = A[1] * B[1] + A[3] * B[3];
}

/*
  Compute C += A*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3MatMultAdd(const TacsScalar A[], const TacsScalar B[],
                                    TacsScalar C[]) {
  C[0] += A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
  C[3] += A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
  C[6] += A[6] * B[0] + A[7] * B[3] + A[8] * B[6];

  C[1] += A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
  C[4] += A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
  C[7] += A[6] * B[1] + A[7] * B[4] + A[8] * B[7];

  C[2] += A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
  C[5] += A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
  C[8] += A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

/*
  Compute C =+ A*B

  input:
  A:   the first 2x2 input matrix in row-major order
  B:   the second 2x2 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat2x2MatMultAdd(const TacsScalar A[], const TacsScalar B[],
                                    TacsScalar C[]) {
  C[0] += A[0] * B[0] + A[1] * B[2];
  C[2] += A[2] * B[0] + A[3] * B[2];

  C[1] += A[0] * B[1] + A[1] * B[3];
  C[3] += A[2] * B[1] + A[3] * B[3];
}

/*
  Compute D = A^{T}*B*C
*/
static inline void mat3x3TransMatTransform(const TacsScalar A[],
                                           const TacsScalar B[],
                                           const TacsScalar C[],
                                           TacsScalar D[]) {
  TacsScalar tmp[9];
  mat3x3TransMatMult(A, B, tmp);
  mat3x3MatMult(tmp, C, D);
}

/*
  Compute D += A^{T}*B*C
*/
static inline void mat3x3TransMatTransformAdd(const TacsScalar A[],
                                              const TacsScalar B[],
                                              const TacsScalar C[],
                                              TacsScalar D[]) {
  TacsScalar tmp[9];
  mat3x3TransMatMult(A, B, tmp);
  mat3x3MatMultAdd(tmp, C, D);
}

/*
  Compute C += A^{T}*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3TransMatMultAdd(const TacsScalar A[],
                                         const TacsScalar B[], TacsScalar C[]) {
  C[0] += A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
  C[1] += A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
  C[2] += A[0] * B[2] + A[3] * B[5] + A[6] * B[8];

  C[3] += A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
  C[4] += A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
  C[5] += A[1] * B[2] + A[4] * B[5] + A[7] * B[8];

  C[6] += A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
  C[7] += A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
  C[8] += A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
}

/*
  Compute the transformation A = T*S*T^{T}

  input:
  T:  the 3x3 transformation
  S:  the 3x3 flattened symmetric matrix

  output:
  A:  the 3x3 flattened symmetric matrix
*/
static inline void mat3x3SymmTransform(const TacsScalar T[],
                                       const TacsScalar S[], TacsScalar A[]) {
  // Compute W = S*T^{T}
  // [S[0] S[1] S[2]][T[0] T[3] T[6]]
  // [S[1] S[3] S[4]][T[1] T[4] T[7]]
  // [S[2] S[4] S[5]][T[2] T[5] T[8]]

  TacsScalar W[9];
  W[0] = S[0] * T[0] + S[1] * T[1] + S[2] * T[2];
  W[1] = S[0] * T[3] + S[1] * T[4] + S[2] * T[5];
  W[2] = S[0] * T[6] + S[1] * T[7] + S[2] * T[8];

  W[3] = S[1] * T[0] + S[3] * T[1] + S[4] * T[2];
  W[4] = S[1] * T[3] + S[3] * T[4] + S[4] * T[5];
  W[5] = S[1] * T[6] + S[3] * T[7] + S[4] * T[8];

  W[6] = S[2] * T[0] + S[4] * T[1] + S[5] * T[2];
  W[7] = S[2] * T[3] + S[4] * T[4] + S[5] * T[5];
  W[8] = S[2] * T[6] + S[4] * T[7] + S[5] * T[8];

  // Compute the symmetric part of T*W
  // [T[0] T[1] T[2]][W[0] W[1] W[2]]
  // [T[3] T[4] T[5]][W[3] W[4] W[5]]
  // [T[6] T[6] T[8]][W[6] W[7] W[8]]
  A[0] = T[0] * W[0] + T[1] * W[3] + T[2] * W[6];
  A[1] = T[0] * W[1] + T[1] * W[4] + T[2] * W[7];
  A[2] = T[0] * W[2] + T[1] * W[5] + T[2] * W[8];

  A[3] = T[3] * W[1] + T[4] * W[4] + T[5] * W[7];
  A[4] = T[3] * W[2] + T[4] * W[5] + T[5] * W[8];

  A[5] = T[6] * W[2] + T[7] * W[5] + T[8] * W[8];
}

/*
  Compute the transformation A = T^{T}*S*T

  input:
  T:  the 3x3 transformation
  S:  the 3x3 flattened symmetric matrix

  output:
  A:  the 3x3 flattened symmetric matrix
*/
static inline void mat3x3SymmTransformTranspose(const TacsScalar T[],
                                                const TacsScalar S[],
                                                TacsScalar A[]) {
  // Compute W = S*T
  // [S[0] S[1] S[2]][T[0] T[1] T[2]]
  // [S[1] S[3] S[4]][T[3] T[4] T[5]]
  // [S[2] S[4] S[5]][T[6] T[7] T[8]]

  TacsScalar W[9];
  W[0] = S[0] * T[0] + S[1] * T[3] + S[2] * T[6];
  W[1] = S[0] * T[1] + S[1] * T[4] + S[2] * T[7];
  W[2] = S[0] * T[2] + S[1] * T[5] + S[2] * T[8];

  W[3] = S[1] * T[0] + S[3] * T[3] + S[4] * T[6];
  W[4] = S[1] * T[1] + S[3] * T[4] + S[4] * T[7];
  W[5] = S[1] * T[2] + S[3] * T[5] + S[4] * T[8];

  W[6] = S[2] * T[0] + S[4] * T[3] + S[5] * T[6];
  W[7] = S[2] * T[1] + S[4] * T[4] + S[5] * T[7];
  W[8] = S[2] * T[2] + S[4] * T[5] + S[5] * T[8];

  // Compute the symmetric part of T^{T}*W
  // [T[0] T[3] T[6]][W[0] W[1] W[2]]
  // [T[1] T[4] T[7]][W[3] W[4] W[5]]
  // [T[2] T[5] T[8]][W[6] W[7] W[8]]
  A[0] = T[0] * W[0] + T[3] * W[3] + T[6] * W[6];
  A[1] = T[0] * W[1] + T[3] * W[4] + T[6] * W[7];
  A[2] = T[0] * W[2] + T[3] * W[5] + T[6] * W[8];

  A[3] = T[1] * W[1] + T[4] * W[4] + T[7] * W[7];
  A[4] = T[1] * W[2] + T[4] * W[5] + T[7] * W[8];

  A[5] = T[2] * W[2] + T[5] * W[5] + T[8] * W[8];
}

/*
  Compute the derivative of the transformation A = T*S*T^{T}

  input:
  T:   the 3x3 transformation
  dA:  the derivative w.r.t. the 3x3 flattened symmetric matrix

  output:
  dS:  the derivative w.r.t. the 3x3 flattened symmetric matrix
*/
static inline void mat3x3SymmTransformSens(const TacsScalar T[],
                                           const TacsScalar dA[],
                                           TacsScalar dS[]) {
  TacsScalar dW[9];

  dW[0] = T[0] * dA[0];
  dW[1] = T[0] * dA[1] + T[3] * dA[3];
  dW[2] = T[0] * dA[2] + T[3] * dA[4] + T[6] * dA[5];

  dW[3] = T[1] * dA[0];
  dW[4] = T[1] * dA[1] + T[4] * dA[3];
  dW[5] = T[1] * dA[2] + T[4] * dA[4] + T[7] * dA[5];

  dW[6] = T[2] * dA[0];
  dW[7] = T[2] * dA[1] + T[5] * dA[3];
  dW[8] = T[2] * dA[2] + T[5] * dA[4] + T[8] * dA[5];

  dS[0] = (T[0] * dW[0] + T[3] * dW[1] + T[6] * dW[2]);
  dS[1] = (T[1] * dW[0] + T[4] * dW[1] + T[7] * dW[2] + T[0] * dW[3] +
           T[3] * dW[4] + T[6] * dW[5]);
  dS[2] = (T[2] * dW[0] + T[5] * dW[1] + T[8] * dW[2] + T[0] * dW[6] +
           T[3] * dW[7] + T[6] * dW[8]);

  dS[3] = (T[1] * dW[3] + T[4] * dW[4] + T[7] * dW[5]);
  dS[4] = (T[2] * dW[3] + T[5] * dW[4] + T[8] * dW[5] + T[1] * dW[6] +
           T[4] * dW[7] + T[7] * dW[8]);

  dS[5] = (T[2] * dW[6] + T[5] * dW[7] + T[8] * dW[8]);
}

/*
  Compute the derivative of the transformation A = T^{T}*S*T

  input:
  T:   the 3x3 transformation
  dA:  the derivative w.r.t. the 3x3 flattened symmetric matrix

  output:
  dS:  the derivative w.r.t. the 3x3 flattened symmetric matrix
*/
static inline void mat3x3SymmTransformTransSens(const TacsScalar T[],
                                                const TacsScalar dA[],
                                                TacsScalar dS[]) {
  TacsScalar dW[9];
  dW[0] = T[0] * dA[0];
  dW[1] = T[0] * dA[1] + T[1] * dA[3];
  dW[2] = T[0] * dA[2] + T[1] * dA[4] + T[2] * dA[5];

  dW[3] = T[3] * dA[0];
  dW[4] = T[3] * dA[1] + T[4] * dA[3];
  dW[5] = T[3] * dA[2] + T[4] * dA[4] + T[5] * dA[5];

  dW[6] = T[6] * dA[0];
  dW[7] = T[6] * dA[1] + T[7] * dA[3];
  dW[8] = T[6] * dA[2] + T[7] * dA[4] + T[8] * dA[5];

  dS[0] = (T[0] * dW[0] + T[1] * dW[1] + T[2] * dW[2]);
  dS[1] = (T[3] * dW[0] + T[4] * dW[1] + T[5] * dW[2] + T[0] * dW[3] +
           T[1] * dW[4] + T[2] * dW[5]);
  dS[2] = (T[6] * dW[0] + T[7] * dW[1] + T[8] * dW[2] + T[0] * dW[6] +
           T[1] * dW[7] + T[2] * dW[8]);

  dS[3] = (T[3] * dW[3] + T[4] * dW[4] + T[5] * dW[5]);
  dS[4] = (T[6] * dW[3] + T[7] * dW[4] + T[8] * dW[5] + T[3] * dW[6] +
           T[4] * dW[7] + T[5] * dW[8]);

  dS[5] = (T[6] * dW[6] + T[7] * dW[7] + T[8] * dW[8]);
}

/*
  Compute the second derivative of the transformation

  A[i,j] = T[k,i]*S[i,j]*T[j,k]

  input:
  T:   The 3x3 transformation matrix
  d2A: The second derivative of the 3x3 symmetric matrix

  output:
  d2S: The second derivative of the 3x3 symmetric matrix
*/
static inline void mat3x3SymmTransformTransHessian(const TacsScalar T[],
                                                   const TacsScalar d2A[],
                                                   TacsScalar d2S[]) {
  TacsScalar tmp[36];
  const TacsScalar *dA = d2A;
  TacsScalar *dS = tmp;
  for (int i = 0; i < 6; i++) {
    TacsScalar dW[9];
    dW[0] = T[0] * dA[0];
    dW[1] = T[0] * dA[1] + T[1] * dA[3];
    dW[2] = T[0] * dA[2] + T[1] * dA[4] + T[2] * dA[5];

    dW[3] = T[3] * dA[0];
    dW[4] = T[3] * dA[1] + T[4] * dA[3];
    dW[5] = T[3] * dA[2] + T[4] * dA[4] + T[5] * dA[5];

    dW[6] = T[6] * dA[0];
    dW[7] = T[6] * dA[1] + T[7] * dA[3];
    dW[8] = T[6] * dA[2] + T[7] * dA[4] + T[8] * dA[5];

    dS[0] = (T[0] * dW[0] + T[1] * dW[1] + T[2] * dW[2]);
    dS[1] = (T[3] * dW[0] + T[4] * dW[1] + T[5] * dW[2] + T[0] * dW[3] +
             T[1] * dW[4] + T[2] * dW[5]);
    dS[2] = (T[6] * dW[0] + T[7] * dW[1] + T[8] * dW[2] + T[0] * dW[6] +
             T[1] * dW[7] + T[2] * dW[8]);

    dS[3] = (T[3] * dW[3] + T[4] * dW[4] + T[5] * dW[5]);
    dS[4] = (T[6] * dW[3] + T[7] * dW[4] + T[8] * dW[5] + T[3] * dW[6] +
             T[4] * dW[7] + T[5] * dW[8]);

    dS[5] = (T[6] * dW[6] + T[7] * dW[7] + T[8] * dW[8]);

    dS += 6;
    dA += 6;
  }

  dA = tmp;
  dS = d2S;
  for (int i = 0; i < 6; i++) {
    TacsScalar dW[9];
    dW[0] = T[0] * dA[0];
    dW[1] = T[0] * dA[6] + T[1] * dA[18];
    dW[2] = T[0] * dA[12] + T[1] * dA[24] + T[2] * dA[30];

    dW[3] = T[3] * dA[0];
    dW[4] = T[3] * dA[6] + T[4] * dA[18];
    dW[5] = T[3] * dA[12] + T[4] * dA[24] + T[5] * dA[30];

    dW[6] = T[6] * dA[0];
    dW[7] = T[6] * dA[6] + T[7] * dA[18];
    dW[8] = T[6] * dA[12] + T[7] * dA[24] + T[8] * dA[30];

    dS[0] = (T[0] * dW[0] + T[1] * dW[1] + T[2] * dW[2]);
    dS[6] = (T[3] * dW[0] + T[4] * dW[1] + T[5] * dW[2] + T[0] * dW[3] +
             T[1] * dW[4] + T[2] * dW[5]);
    dS[12] = (T[6] * dW[0] + T[7] * dW[1] + T[8] * dW[2] + T[0] * dW[6] +
              T[1] * dW[7] + T[2] * dW[8]);

    dS[18] = (T[3] * dW[3] + T[4] * dW[4] + T[5] * dW[5]);
    dS[24] = (T[6] * dW[3] + T[7] * dW[4] + T[8] * dW[5] + T[3] * dW[6] +
              T[4] * dW[7] + T[5] * dW[8]);

    dS[30] = (T[6] * dW[6] + T[7] * dW[7] + T[8] * dW[8]);

    dS++;
    dA++;
  }
}

/*
  Compute D = B*c^{x}

  input:
  a:   the first 3-vector
  B    the 3x3 input matrix
  c:   the second 3-vector

  output:
  D:   the output
*/
static inline void mat3x3MatSkewTransform(const TacsScalar B[],
                                          const TacsScalar c[],
                                          TacsScalar D[]) {
  // [B[0]  B[1]  B[2]][0    -c[2]  c[1]]
  // [B[3]  B[4]  B[5]][c[2]   0   -c[0]]
  // [B[6]  B[7]  B[8]][-c[1] c[0]  a[1]]

  D[0] = c[2] * B[1] - c[1] * B[2];
  D[3] = c[2] * B[4] - c[1] * B[5];
  D[6] = c[2] * B[7] - c[1] * B[8];

  D[1] = c[0] * B[2] - c[2] * B[0];
  D[4] = c[0] * B[5] - c[2] * B[3];
  D[7] = c[0] * B[8] - c[2] * B[6];

  D[2] = c[1] * B[0] - c[0] * B[1];
  D[5] = c[1] * B[3] - c[0] * B[4];
  D[8] = c[1] * B[6] - c[0] * B[7];
}

/*
  Compute D = a^{x}*B

  input:
  a:   the first 3-vector
  B    the 3x3 input matrix
  c:   the second 3-vector

  output:
  D:   the output
*/
static inline void mat3x3SkewMatTransform(const TacsScalar a[],
                                          const TacsScalar B[],
                                          TacsScalar D[]) {
  // [0    -a[2]  a[1]][B[0]  B[1]  B[2]]
  // [a[2]   0   -a[0]][B[3]  B[4]  B[5]]
  // [-a[1] a[0]  a[1]][B[6]  B[7]  B[8]]

  D[0] = a[1] * B[6] - a[2] * B[3];
  D[1] = a[1] * B[7] - a[2] * B[4];
  D[2] = a[1] * B[8] - a[2] * B[5];

  D[3] = a[2] * B[0] - a[0] * B[6];
  D[4] = a[2] * B[1] - a[0] * B[7];
  D[5] = a[2] * B[2] - a[0] * B[8];

  D[6] = a[0] * B[3] - a[1] * B[0];
  D[7] = a[0] * B[4] - a[1] * B[1];
  D[8] = a[0] * B[5] - a[1] * B[2];
}

/*
  Compute D = a^{x}*B*c^{x}

  input:
  a:   the first 3-vector
  B    the 3x3 input matrix
  c:   the second 3-vector

  output:
  D:   the output
*/
static inline void mat3x3SkewMatSkewTransform(const TacsScalar a[],
                                              const TacsScalar B[],
                                              const TacsScalar c[],
                                              TacsScalar D[]) {
  // [0    -a[2]  a[1]][B[0]  B[1]  B[2]][0    -c[2]  c[1]]
  // [a[2]   0   -a[0]][B[3]  B[4]  B[5]][c[2]   0   -c[0]]
  // [-a[1] a[0]  a[1]][B[6]  B[7]  B[8]][-c[1] c[0]  a[1]]

  TacsScalar t[9];
  t[0] = a[1] * B[6] - a[2] * B[3];
  t[1] = a[1] * B[7] - a[2] * B[4];
  t[2] = a[1] * B[8] - a[2] * B[5];

  t[3] = a[2] * B[0] - a[0] * B[6];
  t[4] = a[2] * B[1] - a[0] * B[7];
  t[5] = a[2] * B[2] - a[0] * B[8];

  t[6] = a[0] * B[3] - a[1] * B[0];
  t[7] = a[0] * B[4] - a[1] * B[1];
  t[8] = a[0] * B[5] - a[1] * B[2];

  D[0] = c[2] * t[1] - c[1] * t[2];
  D[3] = c[2] * t[4] - c[1] * t[5];
  D[6] = c[2] * t[7] - c[1] * t[8];

  D[1] = c[0] * t[2] - c[2] * t[0];
  D[4] = c[0] * t[5] - c[2] * t[3];
  D[7] = c[0] * t[8] - c[2] * t[6];

  D[2] = c[1] * t[0] - c[0] * t[1];
  D[5] = c[1] * t[3] - c[0] * t[4];
  D[8] = c[1] * t[6] - c[0] * t[7];
}

/*
  Compute C += A^{T}*B

  input:
  A:   the first 2x2 input matrix in row-major order
  B:   the second 2x2 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat2x2TransMatMultAdd(const TacsScalar A[],
                                         const TacsScalar B[], TacsScalar C[]) {
  C[0] += A[0] * B[0] + A[2] * B[2];
  C[2] += A[1] * B[0] + A[3] * B[2];

  C[1] += A[0] * B[1] + A[2] * B[3];
  C[3] += A[1] * B[1] + A[3] * B[3];
}

/*
  Compute C += A*B^{T}

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void mat3x3MatTransMultAdd(const TacsScalar A[],
                                         const TacsScalar B[], TacsScalar C[]) {
  C[0] += A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
  C[3] += A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
  C[6] += A[6] * B[0] + A[7] * B[1] + A[8] * B[2];

  C[1] += A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
  C[4] += A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
  C[7] += A[6] * B[3] + A[7] * B[4] + A[8] * B[5];

  C[2] += A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
  C[5] += A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
  C[8] += A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
}

/*
  Multiply the matrix by

  [ A[0] A[1] A[2] ][ B[0]  B[1]  B[2]  B[3] ]
  [ A[3] A[4] A[5] ][ B[4]  B[5]  B[6]  B[7] ]
  [ A[6] A[7] A[8] ][ B[8]  B[9]  B[10] B[11] ]

  input:
  A:   a symmetric 3x3 matrix
  B:   a 3x4 matrix in row-major order

  output:
  C:   a 3x4 matrix in row-major order
*/
static inline void matMat3x4Mult(const TacsScalar A[], const TacsScalar B[],
                                 TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[4] + A[2] * B[8];
  C[1] = A[0] * B[1] + A[1] * B[5] + A[2] * B[9];
  C[2] = A[0] * B[2] + A[1] * B[6] + A[2] * B[10];
  C[3] = A[0] * B[3] + A[1] * B[7] + A[2] * B[11];

  C[4] = A[3] * B[0] + A[4] * B[4] + A[5] * B[8];
  C[5] = A[3] * B[1] + A[4] * B[5] + A[5] * B[9];
  C[6] = A[3] * B[2] + A[4] * B[6] + A[5] * B[10];
  C[7] = A[3] * B[3] + A[4] * B[7] + A[5] * B[11];

  C[8] = A[6] * B[0] + A[7] * B[4] + A[8] * B[8];
  C[9] = A[6] * B[1] + A[7] * B[5] + A[8] * B[9];
  C[10] = A[6] * B[2] + A[7] * B[6] + A[8] * B[10];
  C[11] = A[6] * B[3] + A[7] * B[7] + A[8] * B[11];
}

/*
  Multiply a symmetric matrix by a 3x4 row-major matrix

  [ A[0] A[1] A[2] ][ B[0]  B[1]  B[2]  B[3] ]
  [ A[1] A[3] A[4] ][ B[4]  B[5]  B[6]  B[7] ]
  [ A[2] A[4] A[5] ][ B[8]  B[9]  B[10] B[11] ]

  input:
  A:   a symmetric 3x3 matrix
  B:   a 3x4 matrix in row-major order

  output:
  C:   a 3x4 matrix in row-major order
*/
static inline void matSymmMat3x4Mult(const TacsScalar A[], const TacsScalar B[],
                                     TacsScalar C[]) {
  C[0] = A[0] * B[0] + A[1] * B[4] + A[2] * B[8];
  C[1] = A[0] * B[1] + A[1] * B[5] + A[2] * B[9];
  C[2] = A[0] * B[2] + A[1] * B[6] + A[2] * B[10];
  C[3] = A[0] * B[3] + A[1] * B[7] + A[2] * B[11];

  C[4] = A[1] * B[0] + A[3] * B[4] + A[4] * B[8];
  C[5] = A[1] * B[1] + A[3] * B[5] + A[4] * B[9];
  C[6] = A[1] * B[2] + A[3] * B[6] + A[4] * B[10];
  C[7] = A[1] * B[3] + A[3] * B[7] + A[4] * B[11];

  C[8] = A[2] * B[0] + A[4] * B[4] + A[5] * B[8];
  C[9] = A[2] * B[1] + A[4] * B[5] + A[5] * B[9];
  C[10] = A[2] * B[2] + A[4] * B[6] + A[5] * B[10];
  C[11] = A[2] * B[3] + A[4] * B[7] + A[5] * B[11];
}

/*
  Set the 3x3 matrix as a skew symmetric matrix

  C = a*b^{x}

  where a is a scalar and b is 3-vector

  input:
  a:    a scalar
  b:    a 3-vector

  output:
  C:    the result is added to this matrix
*/
static inline void setMatSkew(const TacsScalar a, const TacsScalar b[],
                              TacsScalar C[]) {
  C[0] = 0.0;
  C[1] = -a * b[2];
  C[2] = a * b[1];

  C[3] = a * b[2];
  C[4] = 0.0;
  C[5] = -a * b[0];

  C[6] = -a * b[1];
  C[7] = a * b[0];
  C[8] = 0.0;
}

/*
  Add the skew matrix to a 3x3 matrix

  C += a*b^{x}

  where a is a scalar and b is 3-vector

  input:
  a:    a scalar
  b:    a 3-vector

  output:
  C:    the result is added to this matrix
*/
static inline void addMatSkew(const TacsScalar a, const TacsScalar b[],
                              TacsScalar C[]) {
  C[1] -= a * b[2];
  C[2] += a * b[1];

  C[3] += a * b[2];
  C[5] -= a * b[0];

  C[6] -= a * b[1];
  C[7] += a * b[0];
}

/*
  Set the product of two skew matrices into the matrix

  C = a*b^{x}*c^{x} = a*(c*b^{T} - c^{T}*b*I)

  input:
  a:    a scalar
  b:    a 3-vector

  output:
  C:    the result is added to this matrix
*/
static inline void setMatSkewSkew(const TacsScalar a, const TacsScalar b[],
                                  const TacsScalar c[], TacsScalar C[]) {
  C[0] = -a * (c[1] * b[1] + c[2] * b[2]);
  C[1] = a * c[0] * b[1];
  C[2] = a * c[0] * b[2];

  C[3] = a * c[1] * b[0];
  C[4] = -a * (c[0] * b[0] + c[2] * b[2]);
  C[5] = a * c[1] * b[2];

  C[6] = a * c[2] * b[0];
  C[7] = a * c[2] * b[1];
  C[8] = -a * (c[0] * b[0] + c[1] * b[1]);
}

/*
  Add the product of two skew matrices to the matrix

  C += a*b^{x}*c^{x} = a*(c*b^{T} - c^{T}*b*I)

  input:
  a:    a scalar
  b:    a 3-vector

  output:
  C:    the result is added to this matrix
*/
static inline void addMatSkewSkew(const TacsScalar a, const TacsScalar b[],
                                  const TacsScalar c[], TacsScalar C[]) {
  C[0] -= a * (c[1] * b[1] + c[2] * b[2]);
  C[1] += a * c[0] * b[1];
  C[2] += a * c[0] * b[2];

  C[3] += a * c[1] * b[0];
  C[4] -= a * (c[0] * b[0] + c[2] * b[2]);
  C[5] += a * c[1] * b[2];

  C[6] += a * c[2] * b[0];
  C[7] += a * c[2] * b[1];
  C[8] -= a * (c[0] * b[0] + c[1] * b[1]);
}

/*
  Add a 3x3 block-matrix to a larger block matrix

  D[:,:] += a*A

  input:
  a:    the scalar
  A:    the 3x3 block matrix
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D: the block matrix with in row-major order
*/
static inline void addBlockMat(const TacsScalar a, const TacsScalar A[],
                               TacsScalar D[], const int ldd) {
  D[0] += a * A[0];
  D[1] += a * A[1];
  D[2] += a * A[2];

  D += ldd;
  D[0] += a * A[3];
  D[1] += a * A[4];
  D[2] += a * A[5];

  D += ldd;
  D[0] += a * A[6];
  D[1] += a * A[7];
  D[2] += a * A[8];
}

/*
  Add the transpose of a 3x3 block-matrix to a larger block matrix

  D[:,:] += a*A^{T}

  input:
  a:    the scalar
  A:    the 3x3 block matrix
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix with in row-major order
*/
static inline void addBlockMatTrans(const TacsScalar a, const TacsScalar A[],
                                    TacsScalar D[], const int ldd) {
  D[0] += a * A[0];
  D[1] += a * A[3];
  D[2] += a * A[6];

  D += ldd;
  D[0] += a * A[1];
  D[1] += a * A[4];
  D[2] += a * A[7];

  D += ldd;
  D[0] += a * A[2];
  D[1] += a * A[5];
  D[2] += a * A[8];
}

/*
  Add a vector to the matrix

  D[:,:] += a*A

  input:
  a:    the scalar
  A:    the 3-vector to add to the matrix
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix with in row-major order
*/
static inline void addVecMat(const TacsScalar a, const TacsScalar A[],
                             TacsScalar D[], const int ldd) {
  D[0] += a * A[0];
  D += ldd;
  D[0] += a * A[1];
  D += ldd;
  D[0] += a * A[2];
}

/*
  Add a 3x3 block-matrix to a larger block matrix

  D[:,:] += a*A

  input:
  a:    the scalar
  A:    the 3x3 block matrix
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D: the block matrix with in row-major order
*/
static inline void addBlockSymmMat(const TacsScalar a, const TacsScalar A[],
                                   TacsScalar D[], const int ldd) {
  D[0] += a * A[0];
  D[1] += a * A[1];
  D[2] += a * A[2];

  D += ldd;
  D[0] += a * A[1];
  D[1] += a * A[3];
  D[2] += a * A[4];

  D += ldd;
  D[0] += a * A[2];
  D[1] += a * A[4];
  D[2] += a * A[5];
}

/*
  Add a scalar multiple of the identity matrix to the block matrix

  input:
  a:    the scalar
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix in row-major order
*/
static inline void addBlockIdent(const TacsScalar a, TacsScalar D[],
                                 const int ldd) {
  D[0] += a;

  D += ldd;
  D[1] += a;

  D += ldd;
  D[2] += a;
}

/*
  Add a 3x3 scalar skew-symmetric matrix to the D matrix

  input:
  a:    the scalar
  x:    the vector
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix in row-major order
*/
static inline void addBlockSkew(const TacsScalar a, const TacsScalar x[],
                                TacsScalar D[], const int ldd) {
  D[1] -= a * x[2];
  D[2] += a * x[1];

  D += ldd;
  D[0] += a * x[2];
  D[2] -= a * x[0];

  D += ldd;
  D[0] -= a * x[1];
  D[1] += a * x[0];
}

/*
  Add the product of two skew-symm matrices to the D matrix

  D[:,:] += a*x^{x}*y^{x} = a*(y*x^{T} - I*x^{T}*y)

  input:
  a:    the scalar
  x:    the first vector
  y:    the second vectorr
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix in row-major order
*/
static inline void addBlockSkewSkew(const TacsScalar a, const TacsScalar x[],
                                    const TacsScalar y[], TacsScalar D[],
                                    const int ldd) {
  D[0] -= a * (x[1] * y[1] + x[2] * y[2]);
  D[1] += a * y[0] * x[1];
  D[2] += a * y[0] * x[2];

  D += ldd;
  D[0] += a * y[1] * x[0];
  D[1] -= a * (x[0] * y[0] + x[2] * y[2]);
  D[2] += a * y[1] * x[2];

  D += ldd;
  D[0] += a * y[2] * x[0];
  D[1] += a * y[2] * x[1];
  D[2] -= a * (x[0] * y[0] + x[1] * y[1]);
}

/*
  Compute the determinant of a 3x3 matrix

  input:
  A:        a 3x3 matrix in row-major order

  returns:  the determinant of A
*/
static inline TacsScalar det3x3(const TacsScalar A[]) {
  return (A[8] * (A[0] * A[4] - A[3] * A[1]) -
          A[7] * (A[0] * A[5] - A[3] * A[2]) +
          A[6] * (A[1] * A[5] - A[2] * A[4]));
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void det3x3Sens(const TacsScalar A[], TacsScalar Ad[]) {
  Ad[0] = A[8] * A[4] - A[7] * A[5];
  Ad[1] = A[6] * A[5] - A[8] * A[3];
  Ad[2] = A[7] * A[3] - A[6] * A[4];

  Ad[3] = A[7] * A[2] - A[8] * A[1];
  Ad[4] = A[8] * A[0] - A[6] * A[2];
  Ad[5] = A[6] * A[1] - A[7] * A[0];

  Ad[6] = A[1] * A[5] - A[2] * A[4];
  Ad[7] = A[3] * A[2] - A[0] * A[5];
  Ad[8] = A[0] * A[4] - A[3] * A[1];
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void addDet3x3Sens(const TacsScalar s, const TacsScalar A[],
                                 TacsScalar Ad[]) {
  Ad[0] += s * (A[8] * A[4] - A[7] * A[5]);
  Ad[1] += s * (A[6] * A[5] - A[8] * A[3]);
  Ad[2] += s * (A[7] * A[3] - A[6] * A[4]);

  Ad[3] += s * (A[7] * A[2] - A[8] * A[1]);
  Ad[4] += s * (A[8] * A[0] - A[6] * A[2]);
  Ad[5] += s * (A[6] * A[1] - A[7] * A[0]);

  Ad[6] += s * (A[1] * A[5] - A[2] * A[4]);
  Ad[7] += s * (A[3] * A[2] - A[0] * A[5]);
  Ad[8] += s * (A[0] * A[4] - A[3] * A[1]);
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void det3x32ndSens(const TacsScalar s, const TacsScalar A[],
                                 TacsScalar Ad[]) {
  // Ad[0] = s*(A[8]*A[4] - A[7]*A[5]);
  Ad[0] = 0.0;
  Ad[1] = 0.0;
  Ad[2] = 0.0;
  Ad[3] = 0.0;
  Ad[4] = s * A[8];
  Ad[5] = -s * A[7];
  Ad[6] = 0.0;
  Ad[7] = -s * A[5];
  Ad[8] = s * A[4];
  Ad += 9;

  // Ad[1] += s*(A[6]*A[5] - A[8]*A[3]);
  Ad[0] = 0.0;
  Ad[1] = 0.0;
  Ad[2] = 0.0;
  Ad[3] = -s * A[8];
  Ad[4] = 0.0;
  Ad[5] = s * A[6];
  Ad[6] = s * A[5];
  Ad[7] = 0.0;
  Ad[8] = -s * A[3];
  ;
  Ad += 9;

  // Ad[2] += s*(A[7]*A[3] - A[6]*A[4]);
  Ad[0] = 0.0;
  Ad[1] = 0.0;
  Ad[2] = 0.0;
  Ad[3] = s * A[7];
  Ad[4] = -s * A[6];
  Ad[5] = 0.0;
  Ad[6] = -s * A[4];
  Ad[7] = s * A[3];
  Ad[8] = 0.0;
  Ad += 9;

  // Ad[3] += s*(A[7]*A[2] - A[8]*A[1]);
  Ad[0] = 0.0;
  Ad[1] = -s * A[8];
  Ad[2] = s * A[7];
  Ad[3] = 0.0;
  Ad[4] = 0.0;
  Ad[5] = 0.0;
  Ad[6] = 0.0;
  Ad[7] = s * A[2];
  Ad[8] = -s * A[1];
  Ad += 9;

  // Ad[4] += s*(A[8]*A[0] - A[6]*A[2]);
  Ad[0] = s * A[8];
  Ad[1] = 0.0;
  Ad[2] = -s * A[6];
  Ad[3] = 0.0;
  Ad[4] = 0.0;
  Ad[5] = 0.0;
  Ad[6] = -s * A[2];
  Ad[7] = 0.0;
  Ad[8] = s * A[0];
  Ad += 9;

  // Ad[5] += s*(A[6]*A[1] - A[7]*A[0]);
  Ad[0] = -s * A[7];
  Ad[1] = s * A[6];
  Ad[2] = 0.0;
  Ad[3] = 0.0;
  Ad[4] = 0.0;
  Ad[5] = 0.0;
  Ad[6] = s * A[1];
  Ad[7] = -s * A[0];
  Ad[8] = 0.0;
  Ad += 9;

  // Ad[6] += s*(A[1]*A[5] - A[2]*A[4]);
  Ad[0] = 0.0;
  Ad[1] = s * A[5];
  Ad[2] = -s * A[4];
  Ad[3] = 0.0;
  Ad[4] = -s * A[2];
  Ad[5] = s * A[1];
  Ad[6] = 0.0;
  Ad[7] = 0.0;
  Ad[8] = 0.0;
  Ad += 9;

  // Ad[7] += s*(A[3]*A[2] - A[0]*A[5]);
  Ad[0] = -s * A[5];
  Ad[1] = 0.0;
  Ad[2] = s * A[3];
  Ad[3] = s * A[2];
  Ad[4] = 0.0;
  Ad[5] = -s * A[0];
  Ad[6] = 0.0;
  Ad[7] = 0.0;
  Ad[8] = 0.0;
  Ad += 9;

  // Ad[8] += s*(A[0]*A[4] - A[3]*A[1]);
  Ad[0] = s * A[4];
  Ad[1] = -s * A[3];
  Ad[2] = 0.0;
  Ad[3] = -s * A[1];
  Ad[4] = s * A[0];
  Ad[5] = 0.0;
  Ad[6] = 0.0;
  Ad[7] = 0.0;
  Ad[8] = 0.0;
}

/*
  Compute the inverse of a 3x3 matrix

  input:
  A:          a 3x3 matrix in row major order

  output:
  Ainv:       the inverse of the 3x3 matrix

  returns:    the determinant of A
*/
static inline TacsScalar inv3x3(const TacsScalar A[], TacsScalar Ainv[]) {
  TacsScalar det =
      (A[8] * (A[0] * A[4] - A[3] * A[1]) - A[7] * (A[0] * A[5] - A[3] * A[2]) +
       A[6] * (A[1] * A[5] - A[2] * A[4]));
  TacsScalar detinv = 1.0 / det;

  Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * detinv;
  Ainv[1] = -(A[1] * A[8] - A[2] * A[7]) * detinv;
  Ainv[2] = (A[1] * A[5] - A[2] * A[4]) * detinv;

  Ainv[3] = -(A[3] * A[8] - A[5] * A[6]) * detinv;
  Ainv[4] = (A[0] * A[8] - A[2] * A[6]) * detinv;
  Ainv[5] = -(A[0] * A[5] - A[2] * A[3]) * detinv;

  Ainv[6] = (A[3] * A[7] - A[4] * A[6]) * detinv;
  Ainv[7] = -(A[0] * A[7] - A[1] * A[6]) * detinv;
  Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * detinv;

  return det;
}

/*
  Compute the sensitivity of the 3x3 inverse matrix

  input:
  Ainv:   The 3x3 inverse of the matrix
  Ainvd:  The derivative of the 3x3 inverse

  output:
  Ainvd:      derivative of the inverse of the 3x3 matrix
*/
static inline void inv3x3Sens(const TacsScalar Ainv[], const TacsScalar Ainvd[],
                              TacsScalar Ad[]) {
  // d(Ainv_{kl})/d(A_{ij})
  //  = -Ainv_{kn}*delta_{ni}*delta{mj}*Ainv_{ml}
  //  = -Ainv_{ki}*Ainv_{jl}

  // Ad_{ij}
  //  = d(Ainv_{kl})/d(A_{ij})*Ainvd_{kl}
  //  = -Ainv_{ki}*Ainv_{jl}*Ainvd_{kl}

  // Ad = -Ainv^{T}*Ainvd*Ainv^{T}
  TacsScalar t[9];
  mat3x3TransMatMult(Ainv, Ainvd, t);
  mat3x3MatTransMult(t, Ainv, Ad);

  Ad[0] = -Ad[0];
  Ad[1] = -Ad[1];
  Ad[2] = -Ad[2];
  Ad[3] = -Ad[3];
  Ad[4] = -Ad[4];
  Ad[5] = -Ad[5];
  Ad[6] = -Ad[6];
  Ad[7] = -Ad[7];
  Ad[8] = -Ad[8];
}

/*
  Compute the sensitivity of the 3x3 inverse matrix
  input:
  A: a 3x3 matrix in row major order
  sA: the sensitivity of the input matrix
  output:
  Ainv: the inverse of the matrix
  sAinv: the sensitivity of the matrix
  sh: the sensitivity of the determinant of the matrix
  returns the determinant of A
*/
static inline TacsScalar inv3x3Sens(const TacsScalar A[], const TacsScalar sA[],
                                    TacsScalar Ainv[], TacsScalar sAinv[],
                                    TacsScalar *_sh) {
  TacsScalar h =
      (A[8] * (A[0] * A[4] - A[3] * A[1]) - A[7] * (A[0] * A[5] - A[3] * A[2]) +
       A[6] * (A[1] * A[5] - A[2] * A[4]));
  TacsScalar hinv = 1.0 / h;

  TacsScalar sh;
  sh = sA[8] * (A[0] * A[4] - A[3] * A[1]) +
       A[8] * (A[0] * sA[4] + sA[0] * A[4] - A[3] * sA[1] - sA[3] * A[1]) -
       sA[7] * (A[0] * A[5] - A[3] * A[2]) -
       A[7] * (A[0] * sA[5] + sA[0] * A[5] - A[3] * sA[2] - sA[3] * A[2]) +
       sA[6] * (A[1] * A[5] - A[2] * A[4]) +
       A[6] * (A[1] * sA[5] + sA[1] * A[5] - A[2] * sA[4] - sA[2] * A[4]);

  Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * hinv;
  Ainv[1] = -(A[1] * A[8] - A[2] * A[7]) * hinv;
  Ainv[2] = (A[1] * A[5] - A[2] * A[4]) * hinv;

  Ainv[3] = -(A[3] * A[8] - A[5] * A[6]) * hinv;
  Ainv[4] = (A[0] * A[8] - A[2] * A[6]) * hinv;
  Ainv[5] = -(A[0] * A[5] - A[2] * A[3]) * hinv;

  Ainv[6] = (A[3] * A[7] - A[4] * A[6]) * hinv;
  Ainv[7] = -(A[0] * A[7] - A[1] * A[6]) * hinv;
  Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * hinv;

  for (int i = 0; i < 9; i++) {
    sAinv[i] = -Ainv[i] * hinv * sh;
  }

  sAinv[0] +=
      (A[4] * sA[8] + sA[4] * A[8] - A[5] * sA[7] - sA[5] * A[7]) * hinv;
  sAinv[1] +=
      -(A[1] * sA[8] + sA[1] * A[8] - A[2] * sA[7] - sA[2] * A[7]) * hinv;
  sAinv[2] +=
      (A[1] * sA[5] + sA[1] * A[5] - A[2] * sA[4] - sA[2] * A[4]) * hinv;

  sAinv[3] +=
      -(A[3] * sA[8] + sA[3] * A[8] - A[5] * sA[6] - sA[5] * A[6]) * hinv;
  sAinv[4] +=
      (A[0] * sA[8] + sA[0] * A[8] - A[2] * sA[6] - sA[2] * A[6]) * hinv;
  sAinv[5] +=
      -(A[0] * sA[5] + sA[0] * A[5] - A[2] * sA[3] - sA[2] * A[3]) * hinv;

  sAinv[6] +=
      (A[3] * sA[7] + sA[3] * A[7] - A[4] * sA[6] - sA[4] * A[6]) * hinv;
  sAinv[7] +=
      -(A[0] * sA[7] + sA[0] * A[7] - A[1] * sA[6] - sA[1] * A[6]) * hinv;
  sAinv[8] +=
      (A[0] * sA[4] + sA[0] * A[4] - A[1] * sA[3] - sA[1] * A[3]) * hinv;

  *_sh = sh;
  return h;
}

/*
  Compute the determinant of a 3x3 matrix

  input:
  A:        a 3x3 matrix in row-major order

  returns:  the determinant of A
*/
static inline TacsScalar det2x2(const TacsScalar A[]) {
  return (A[0] * A[3] - A[1] * A[2]);
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void det2x2Sens(const TacsScalar A[], TacsScalar Ad[]) {
  Ad[0] = A[3];
  Ad[1] = -A[2];
  Ad[2] = -A[1];
  Ad[3] = A[0];
}

/*
  Compute the inverse of a 3x3 matrix

  input:
  A:          a 3x3 matrix in row major order

  output:
  Ainv:       the inverse of the 3x3 matrix

  returns:    the determinant of A
*/
static inline TacsScalar inv2x2(const TacsScalar A[], TacsScalar Ainv[]) {
  TacsScalar det = A[0] * A[3] - A[1] * A[2];
  TacsScalar detinv = 1.0 / det;

  Ainv[0] = A[3] * detinv;
  Ainv[1] = -A[1] * detinv;

  Ainv[2] = -A[2] * detinv;
  Ainv[3] = A[0] * detinv;

  return det;
}

/*
  Compute the sensitivity of the 3x3 inverse matrix

  input:
  Ainv:   The 3x3 inverse of the matrix
  Ainvd:  The derivative of the 3x3 inverse

  output:
  Ainvd:      derivative of the inverse of the 3x3 matrix
*/
static inline void inv2x2Sens(const TacsScalar Ainv[], const TacsScalar Ainvd[],
                              TacsScalar Ad[]) {
  // d(Ainv_{kl})/d(A_{ij})
  //  = -Ainv_{kn}*delta_{ni}*delta{mj}*Ainv_{ml}
  //  = -Ainv_{ki}*Ainv_{jl}

  // Ad_{ij}
  //  = d(Ainv_{kl})/d(A_{ij})*Ainvd_{kl}
  //  = -Ainv_{ki}*Ainv_{jl}*Ainvd_{kl}

  // Ad = -Ainv^{T}*Ainvd*Ainv^{T}
  TacsScalar t[4];
  mat2x2TransMatMult(Ainv, Ainvd, t);
  mat2x2MatTransMult(t, Ainv, Ad);

  Ad[0] = -Ad[0];
  Ad[1] = -Ad[1];
  Ad[2] = -Ad[2];
  Ad[3] = -Ad[3];
}

#endif  // TACS_ELEMENT_ALGEBRA_H
