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

/*
  Compute y <- A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void matMult( const TacsScalar A[],
                            const TacsScalar x[],
                            TacsScalar y[] ){
  y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
  y[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

/*
  Compute y <- A^{T}*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  output:
  y:   the resulting vector
*/
static inline void matMultTrans( const TacsScalar A[],
                                 const TacsScalar x[],
                                 TacsScalar y[] ){
  y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
  y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
  y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void matMultAdd( const TacsScalar A[],
                               const TacsScalar x[],
                               TacsScalar y[] ){
  y[0] += A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] += A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
  y[2] += A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void matMultTransAdd( const TacsScalar A[],
                                    const TacsScalar x[],
                                    TacsScalar y[] ){
  y[0] += A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
  y[1] += A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
  y[2] += A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
}

/*
  Compute y^{T}*A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector
  y:   the resulting vector

  returns:  the inner product
*/
static inline TacsScalar matSymmInner( const TacsScalar A[],
                                       const TacsScalar x[],
                                       const TacsScalar y[] ){
  return (y[0]*(A[0]*x[0] + A[1]*x[1] + A[2]*x[2]) +
          y[1]*(A[1]*x[0] + A[3]*x[1] + A[4]*x[2]) +
          y[2]*(A[2]*x[0] + A[4]*x[1] + A[5]*x[2]));
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void matSymmMult( const TacsScalar A[],
                                const TacsScalar x[],
                                TacsScalar y[] ){
  y[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] = A[1]*x[0] + A[3]*x[1] + A[4]*x[2];
  y[2] = A[2]*x[0] + A[4]*x[1] + A[5]*x[2];
}

/*
  Compute y <- y + A*x

  input:
  A:   the 3x3 input matrix in row-major order
  x:   the input 3-vector

  in/out:
  y:   the resulting vector
*/
static inline void matSymmMultAdd( const TacsScalar A[],
                                   const TacsScalar x[],
                                   TacsScalar y[] ){
  y[0] += A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
  y[1] += A[1]*x[0] + A[3]*x[1] + A[4]*x[2];
  y[2] += A[2]*x[0] + A[4]*x[1] + A[5]*x[2];
}

/*
  Compute C = A*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matMatMult( const TacsScalar A[],
                               const TacsScalar B[],
                               TacsScalar C[] ){
  C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
  C[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
  C[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];

  C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
  C[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
  C[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];

  C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
  C[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
  C[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

/*
  Compute C = A*B^{T}

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matMatTransMult( const TacsScalar A[],
                                    const TacsScalar B[],
                                    TacsScalar C[] ){
  C[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  C[3] = A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
  C[6] = A[6]*B[0] + A[7]*B[1] + A[8]*B[2];

  C[1] = A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
  C[4] = A[3]*B[3] + A[4]*B[4] + A[5]*B[5];
  C[7] = A[6]*B[3] + A[7]*B[4] + A[8]*B[5];

  C[2] = A[0]*B[6] + A[1]*B[7] + A[2]*B[8];
  C[5] = A[3]*B[6] + A[4]*B[7] + A[5]*B[8];
  C[8] = A[6]*B[6] + A[7]*B[7] + A[8]*B[8];
}

/*
  Compute C = A^{T}*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matTransMatMult( const TacsScalar A[],
                                    const TacsScalar B[],
                                    TacsScalar C[] ){
  C[0] = A[0]*B[0] + A[3]*B[3] + A[6]*B[6];
  C[1] = A[0]*B[1] + A[3]*B[4] + A[6]*B[7];
  C[2] = A[0]*B[2] + A[3]*B[5] + A[6]*B[8];

  C[3] = A[1]*B[0] + A[4]*B[3] + A[7]*B[6];
  C[4] = A[1]*B[1] + A[4]*B[4] + A[7]*B[7];
  C[5] = A[1]*B[2] + A[4]*B[5] + A[7]*B[8];

  C[6] = A[2]*B[0] + A[5]*B[3] + A[8]*B[6];
  C[7] = A[2]*B[1] + A[5]*B[4] + A[8]*B[7];
  C[8] = A[2]*B[2] + A[5]*B[5] + A[8]*B[8];
}

/*
  Compute C += A^{T}*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matTransMatMultAdd( const TacsScalar A[],
                                       const TacsScalar B[],
                                       TacsScalar C[] ){
  C[0] += A[0]*B[0] + A[3]*B[3] + A[6]*B[6];
  C[1] += A[0]*B[1] + A[3]*B[4] + A[6]*B[7];
  C[2] += A[0]*B[2] + A[3]*B[5] + A[6]*B[8];

  C[3] += A[1]*B[0] + A[4]*B[3] + A[7]*B[6];
  C[4] += A[1]*B[1] + A[4]*B[4] + A[7]*B[7];
  C[5] += A[1]*B[2] + A[4]*B[5] + A[7]*B[8];

  C[6] += A[2]*B[0] + A[5]*B[3] + A[8]*B[6];
  C[7] += A[2]*B[1] + A[5]*B[4] + A[8]*B[7];
  C[8] += A[2]*B[2] + A[5]*B[5] + A[8]*B[8];
}

/*
  Compute C += A*B

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matMatMultAdd( const TacsScalar A[],
                                  const TacsScalar B[],
                                  TacsScalar C[] ){
  C[0] += A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
  C[3] += A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
  C[6] += A[6]*B[0] + A[7]*B[3] + A[8]*B[6];

  C[1] += A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
  C[4] += A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
  C[7] += A[6]*B[1] + A[7]*B[4] + A[8]*B[7];

  C[2] += A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
  C[5] += A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
  C[8] += A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

/*
  Compute C += A*B^{T}

  input:
  A:   the first 3x3 input matrix in row-major order
  B:   the second 3x3 input matrix in row-major order

  output:
  C:   the resulting matrix
*/
static inline void matMatTransMultAdd( const TacsScalar A[],
                                       const TacsScalar B[],
                                       TacsScalar C[] ){
  C[0] += A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  C[3] += A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
  C[6] += A[6]*B[0] + A[7]*B[1] + A[8]*B[2];

  C[1] += A[0]*B[3] + A[1]*B[4] + A[2]*B[5];
  C[4] += A[3]*B[3] + A[4]*B[4] + A[5]*B[5];
  C[7] += A[6]*B[3] + A[7]*B[4] + A[8]*B[5];

  C[2] += A[0]*B[6] + A[1]*B[7] + A[2]*B[8];
  C[5] += A[3]*B[6] + A[4]*B[7] + A[5]*B[8];
  C[8] += A[6]*B[6] + A[7]*B[7] + A[8]*B[8];
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
static inline void matMat3x4Mult( const TacsScalar A[],
                                  const TacsScalar B[],
                                  TacsScalar C[] ){
  C[0] = A[0]*B[0] + A[1]*B[4] + A[2]*B[8];
  C[1] = A[0]*B[1] + A[1]*B[5] + A[2]*B[9];
  C[2] = A[0]*B[2] + A[1]*B[6] + A[2]*B[10];
  C[3] = A[0]*B[3] + A[1]*B[7] + A[2]*B[11];

  C[4] = A[3]*B[0] + A[4]*B[4] + A[5]*B[8];
  C[5] = A[3]*B[1] + A[4]*B[5] + A[5]*B[9];
  C[6] = A[3]*B[2] + A[4]*B[6] + A[5]*B[10];
  C[7] = A[3]*B[3] + A[4]*B[7] + A[5]*B[11];

  C[8] = A[6]*B[0] + A[7]*B[4] + A[8]*B[8];
  C[9] = A[6]*B[1] + A[7]*B[5] + A[8]*B[9];
  C[10]= A[6]*B[2] + A[7]*B[6] + A[8]*B[10];
  C[11]= A[6]*B[3] + A[7]*B[7] + A[8]*B[11];
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
static inline void matSymmMat3x4Mult( const TacsScalar A[],
                                      const TacsScalar B[],
                                      TacsScalar C[] ){
  C[0] = A[0]*B[0] + A[1]*B[4] + A[2]*B[8];
  C[1] = A[0]*B[1] + A[1]*B[5] + A[2]*B[9];
  C[2] = A[0]*B[2] + A[1]*B[6] + A[2]*B[10];
  C[3] = A[0]*B[3] + A[1]*B[7] + A[2]*B[11];

  C[4] = A[1]*B[0] + A[3]*B[4] + A[4]*B[8];
  C[5] = A[1]*B[1] + A[3]*B[5] + A[4]*B[9];
  C[6] = A[1]*B[2] + A[3]*B[6] + A[4]*B[10];
  C[7] = A[1]*B[3] + A[3]*B[7] + A[4]*B[11];

  C[8] = A[2]*B[0] + A[4]*B[4] + A[5]*B[8];
  C[9] = A[2]*B[1] + A[4]*B[5] + A[5]*B[9];
  C[10]= A[2]*B[2] + A[4]*B[6] + A[5]*B[10];
  C[11]= A[2]*B[3] + A[4]*B[7] + A[5]*B[11];
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
static inline void crossProduct( const TacsScalar a,
                                 const TacsScalar x[],
                                 const TacsScalar y[],
                                 TacsScalar out[] ){
  out[0] = a*(x[1]*y[2] - x[2]*y[1]);
  out[1] = a*(x[2]*y[0] - x[0]*y[2]);
  out[2] = a*(x[0]*y[1] - x[1]*y[0]);
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
static inline void crossProductAdd( const TacsScalar a,
                                    const TacsScalar x[],
                                    const TacsScalar y[],
                                    TacsScalar out[] ){
  out[0] += a*(x[1]*y[2] - x[2]*y[1]);
  out[1] += a*(x[2]*y[0] - x[0]*y[2]);
  out[2] += a*(x[0]*y[1] - x[1]*y[0]);
}

/*
  Scale the vector by the given scalar

  input
  a:   the scalar
  x:   the vector
*/
static inline void vecScale( const TacsScalar a,
                             TacsScalar x[] ){
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
static inline TacsScalar vecDot( const TacsScalar x[],
                                 const TacsScalar y[] ){
  return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
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
static inline void vecAxpy( const TacsScalar a,
                            const TacsScalar x[],
                            TacsScalar y[] ){
  y[0] += a*x[0];
  y[1] += a*x[1];
  y[2] += a*x[2];
}

/*
  Compute the derivative of x/||x||_{2} w.r.t. x using x and the norm
  of x.

  This code computes a 3x3 matrix that takes the form:

  d(x/||x||_{2})/dx = (I*||x||^2 + x*x^{T})/||x||^3
*/
static inline void vecNormDeriv( TacsScalar nrm,
                                 const TacsScalar x[],
                                 TacsScalar D[] ){
  TacsScalar s = 1.0/(nrm*nrm*nrm);
  TacsScalar t = nrm*nrm;

  D[0] = s*(t - x[0]*x[0]);
  D[1] =-s*x[0]*x[1];
  D[2] =-s*x[0]*x[2];

  D[3] =-s*x[1]*x[0];
  D[4] = s*(t - x[1]*x[1]);
  D[5] =-s*x[1]*x[2];

  D[6] =-s*x[2]*x[0];
  D[7] =-s*x[2]*x[1];
  D[8] = s*(t - x[2]*x[2]);
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
static inline void setMatSkew( const TacsScalar a,
                               const TacsScalar b[],
                               TacsScalar C[] ){
  C[0] = 0.0;
  C[1] = -a*b[2];
  C[2] = a*b[1];

  C[3] = a*b[2];
  C[4] = 0.0;
  C[5] = -a*b[0];

  C[6] = -a*b[1];
  C[7] = a*b[0];
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
static inline void addMatSkew( const TacsScalar a,
                               const TacsScalar b[],
                               TacsScalar C[] ){
  C[1] -= a*b[2];
  C[2] += a*b[1];

  C[3] += a*b[2];
  C[5] -= a*b[0];

  C[6] -= a*b[1];
  C[7] += a*b[0];
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
static inline void setMatSkewSkew( const TacsScalar a,
                                   const TacsScalar b[],
                                   const TacsScalar c[],
                                   TacsScalar C[] ){
  C[0] = -a*(c[1]*b[1] + c[2]*b[2]);
  C[1] = a*c[0]*b[1];
  C[2] = a*c[0]*b[2];

  C[3] = a*c[1]*b[0];
  C[4] = -a*(c[0]*b[0] + c[2]*b[2]);
  C[5] = a*c[1]*b[2];

  C[6] = a*c[2]*b[0];
  C[7] = a*c[2]*b[1];
  C[8] = -a*(c[0]*b[0] + c[1]*b[1]);
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
static inline void addMatSkewSkew( const TacsScalar a,
                                   const TacsScalar b[],
                                   const TacsScalar c[],
                                   TacsScalar C[] ){
  C[0] -= a*(c[1]*b[1] + c[2]*b[2]);
  C[1] += a*c[0]*b[1];
  C[2] += a*c[0]*b[2];

  C[3] += a*c[1]*b[0];
  C[4] -= a*(c[0]*b[0] + c[2]*b[2]);
  C[5] += a*c[1]*b[2];

  C[6] += a*c[2]*b[0];
  C[7] += a*c[2]*b[1];
  C[8] -= a*(c[0]*b[0] + c[1]*b[1]);
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
static inline void addBlockMat( const TacsScalar a,
                                const TacsScalar A[],
                                TacsScalar D[], const int ldd ){
  D[0] += a*A[0];
  D[1] += a*A[1];
  D[2] += a*A[2];

  D += ldd;
  D[0] += a*A[3];
  D[1] += a*A[4];
  D[2] += a*A[5];

  D += ldd;
  D[0] += a*A[6];
  D[1] += a*A[7];
  D[2] += a*A[8];
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
static inline void addBlockMatTrans( const TacsScalar a,
                                     const TacsScalar A[],
                                     TacsScalar D[], const int ldd ){
  D[0] += a*A[0];
  D[1] += a*A[3];
  D[2] += a*A[6];

  D += ldd;
  D[0] += a*A[1];
  D[1] += a*A[4];
  D[2] += a*A[7];

  D += ldd;
  D[0] += a*A[2];
  D[1] += a*A[5];
  D[2] += a*A[8];
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
static inline void addVecMat( const TacsScalar a,
                              const TacsScalar A[],
                              TacsScalar D[],
                              const int ldd ){
  D[0] += a*A[0]; D += ldd;
  D[0] += a*A[1]; D += ldd;
  D[0] += a*A[2];
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
static inline void addBlockSymmMat( const TacsScalar a,
                                    const TacsScalar A[],
                                    TacsScalar D[], const int ldd ){
  D[0] += a*A[0];
  D[1] += a*A[1];
  D[2] += a*A[2];

  D += ldd;
  D[0] += a*A[1];
  D[1] += a*A[3];
  D[2] += a*A[4];

  D += ldd;
  D[0] += a*A[2];
  D[1] += a*A[4];
  D[2] += a*A[5];
}

/*
  Add a scalar multiple of the identity matrix to the block matrix

  input:
  a:    the scalar
  ldd:  the leading row-dimension of the D matrix

  in/out:
  D:    the block matrix in row-major order
*/
static inline void addBlockIdent( const TacsScalar a,
                                  TacsScalar D[], const int ldd ){
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
static inline void addBlockSkew( const TacsScalar a,
                                 const TacsScalar x[],
                                 TacsScalar D[], const int ldd ){
  D[1] -= a*x[2];
  D[2] += a*x[1];

  D += ldd;
  D[0] += a*x[2];
  D[2] -= a*x[0];

  D += ldd;
  D[0] -= a*x[1];
  D[1] += a*x[0];
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
static inline void addBlockSkewSkew( const TacsScalar a,
                                     const TacsScalar x[],
                                     const TacsScalar y[],
                                     TacsScalar D[], const int ldd ){
  D[0] -= a*(x[1]*y[1] + x[2]*y[2]);
  D[1] += a*y[0]*x[1];
  D[2] += a*y[0]*x[2];

  D += ldd;
  D[0] += a*y[1]*x[0];
  D[1] -= a*(x[0]*y[0] + x[2]*y[2]);
  D[2] += a*y[1]*x[2];

  D += ldd;
  D[0] += a*y[2]*x[0];
  D[1] += a*y[2]*x[1];
  D[2] -= a*(x[0]*y[0] + x[1]*y[1]);
}

/*
  Compute the determinant of a 3x3 matrix

  input:
  A:        a 3x3 matrix in row-major order

  returns:  the determinant of A
*/
static inline TacsScalar det3x3( const TacsScalar A[] ){
  return (A[8]*(A[0]*A[4] - A[3]*A[1]) -
          A[7]*(A[0]*A[5] - A[3]*A[2]) +
          A[6]*(A[1]*A[5] - A[2]*A[4]));
}

/*
  Compute the derivative of the determinant with respect to the
  components of A
*/
static inline void det3x3Sens( const TacsScalar A[],
                               TacsScalar Ad[] ){
  Ad[0] = A[8]*A[4] - A[7]*A[5];
  Ad[1] = A[6]*A[5] - A[8]*A[3];
  Ad[2] = A[7]*A[3] - A[6]*A[4];

  Ad[3] = A[7]*A[2] - A[8]*A[1];
  Ad[4] = A[8]*A[0] - A[6]*A[2];
  Ad[5] = A[6]*A[1] - A[7]*A[0];

  Ad[6] = A[1]*A[5] - A[2]*A[4];
  Ad[7] = A[3]*A[2] - A[0]*A[5];
  Ad[8] = A[0]*A[4] - A[3]*A[1];
}

/*
  Compute the inverse of a 3x3 matrix

  input:
  A:          a 3x3 matrix in row major order

  output:
  Ainv:       the inverse of the 3x3 matrix

  returns:    the determinant of A
*/
static inline TacsScalar inv3x3( const TacsScalar A[],
                                 TacsScalar Ainv[] ){
  TacsScalar det = (A[8]*(A[0]*A[4] - A[3]*A[1]) -
                    A[7]*(A[0]*A[5] - A[3]*A[2]) +
                    A[6]*(A[1]*A[5] - A[2]*A[4]));
  TacsScalar detinv = 1.0/det;

  Ainv[0] = (A[4]*A[8] - A[5]*A[7])*detinv;
  Ainv[1] =-(A[1]*A[8] - A[2]*A[7])*detinv;
  Ainv[2] = (A[1]*A[5] - A[2]*A[4])*detinv;

  Ainv[3] =-(A[3]*A[8] - A[5]*A[6])*detinv;
  Ainv[4] = (A[0]*A[8] - A[2]*A[6])*detinv;
  Ainv[5] =-(A[0]*A[5] - A[2]*A[3])*detinv;

  Ainv[6] = (A[3]*A[7] - A[4]*A[6])*detinv;
  Ainv[7] =-(A[0]*A[7] - A[1]*A[6])*detinv;
  Ainv[8] = (A[0]*A[4] - A[1]*A[3])*detinv;

  return det;
}

/*
  Compute the sensitivity of the 3x3 inverse matrix

  input:
  Ad:         a 3x3 matrix of the derivatives of A
  A:          a 3x3 matrix in row major order

  output:
  Ainvd:      derivative of the inverse of the 3x3 matrix
*/
static inline void inv3x3Sens( TacsScalar Ad[],
                               const TacsScalar Ainvd[],
                               const TacsScalar Ainv[] ){
  // d(Ainv_{kl})/d(A_{ij})
  //  = -Ainv_{kn}*delta_{ni}*delta{mj}*Ainv_{ml}
  //  = -Ainv_{ki}*Ainv_{jl}

  // Ad_{ij}
  //  = d(Ainv_{kl})/d(A_{ij})*Ainvd_{kl}
  //  = -Ainv_{ki}*Ainv_{jl}*Ainvd_{kl}

  // Ad = -Ainv^{T}*Ainvd*Ainv^{T}
  TacsScalar t[9];
  matTransMatMult(Ainv, Ainvd, t);
  matMatTransMult(t, Ainv, Ad);

  Ad[0] = -Ad[0];
  Ad[1] = -Ad[1];
  Ad[2] = -Ad[2];
  Ad[3] = -Ad[3];
  Ad[4] = -Ad[4];
  Ad[5] = -Ad[5];
  Ad[6] = -Ad[6];
  Ad[7] = -Ad[7];
  Ad[8] = -Ad[8];
  Ad[9] = -Ad[9];
}

/*
  Compute the inner product with a 3x3 matrix:

  return:  x^{T}*A*y

  input:
  A:   a 3x3 matrix in row-major order
  x:   a 3-vector
  y:   a 3-vector
*/
static inline TacsScalar mat3x3Inner( const TacsScalar A[],
                                      const TacsScalar x[],
                                      const TacsScalar y[] ){
  return (x[0]*(A[0]*y[0] + A[1]*y[1] + A[2]*y[2]) +
          x[1]*(A[3]*y[0] + A[4]*y[1] + A[5]*y[2]) +
          x[2]*(A[6]*y[0] + A[7]*y[1] + A[8]*y[2]));
}

/*
  Given the quaternion parameters, compute the rotation matrix.

  input:
  eta:    the quaternion scalar
  eps:    the quaternion 3-vector

  output:
  C:      the rotation matrix
*/
static inline void computeRotationMat( const TacsScalar eta,
                                       const TacsScalar eps[],
                                       TacsScalar C[] ){
  C[0] = 1.0 - 2.0*(eps[1]*eps[1] + eps[2]*eps[2]);
  C[1] = 2.0*(eps[0]*eps[1] + eta*eps[2]);
  C[2] = 2.0*(eps[0]*eps[2] - eta*eps[1]);

  C[3] = 2.0*(eps[1]*eps[0] - eta*eps[2]);
  C[4] = 1.0 - 2.0*(eps[0]*eps[0] + eps[2]*eps[2]);
  C[5] = 2.0*(eps[1]*eps[2] + eta*eps[0]);

  C[6] = 2.0*(eps[2]*eps[0] + eta*eps[1]);
  C[7] = 2.0*(eps[2]*eps[1] - eta*eps[0]);
  C[8] = 1.0 - 2.0*(eps[0]*eps[0] + eps[1]*eps[1]);
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
static inline void computeRotationMatDeriv( const TacsScalar eta,
                                            const TacsScalar eps[],
                                            const TacsScalar deta,
                                            const TacsScalar deps[],
                                            TacsScalar C[] ){
  C[0] =-4.0*(eps[1]*deps[1] + eps[2]*deps[2]);
  C[1] = 2.0*(eps[0]*deps[1] + deps[0]*eps[1] + eta*deps[2] + deta*eps[2]);
  C[2] = 2.0*(eps[0]*deps[2] + deps[0]*eps[2] - eta*deps[1] - deta*eps[1]);

  C[3] = 2.0*(eps[1]*deps[0] + deps[1]*eps[0] - eta*deps[2] - deta*eps[2]);
  C[4] =-4.0*(eps[0]*deps[0] + eps[2]*deps[2]);
  C[5] = 2.0*(eps[1]*deps[2] + deps[1]*eps[2] + eta*deps[0] + deta*eps[0]);

  C[6] = 2.0*(eps[2]*deps[0] + deps[2]*eps[0] + eta*deps[1] + deta*eps[1]);
  C[7] = 2.0*(eps[2]*deps[1] + deps[2]*eps[1] - eta*deps[0] - deta*eps[0]);
  C[8] =-4.0*(eps[0]*deps[0] + eps[1]*deps[1]);
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
static inline void computeSRateProduct( const TacsScalar eta,
                                        const TacsScalar eps[],
                                        const TacsScalar xeta,
                                        const TacsScalar xeps[],
                                        TacsScalar y[] ){
  crossProduct(-2.0, eps, xeps, y);
  vecAxpy(2.0*eta, xeps, y);
  vecAxpy(-2.0*xeta, eps, y);
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
static inline void addSRateProduct( const TacsScalar eta,
                                    const TacsScalar eps[],
                                    const TacsScalar xeta,
                                    const TacsScalar xeps[],
                                    TacsScalar y[] ){
  crossProductAdd(-2.0, eps, xeps, y);
  vecAxpy(2.0*eta, xeps, y);
  vecAxpy(-2.0*xeta, eps, y);
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
static inline void addSRateTransProduct( const TacsScalar a,
                                         const TacsScalar eta,
                                         const TacsScalar eps[],
                                         const TacsScalar x[],
                                         TacsScalar *yeta,
                                         TacsScalar yeps[] ){
  *yeta -= 2.0*a*vecDot(eps, x);
  crossProductAdd(2.0*a, eps, x, yeps);
  vecAxpy(2.0*a*eta, x, yeps);
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
static inline void computeSRateMat( const TacsScalar eta,
                                    const TacsScalar eps[],
                                    TacsScalar S[] ){
  S[0] = -2.0*eps[0];
  S[1] = 2.0*eta;
  S[2] = 2.0*eps[2];
  S[3] = -2.0*eps[1];

  S[4] = -2.0*eps[1];
  S[5] = -2.0*eps[2];
  S[6] = 2.0*eta;
  S[7] = 2.0*eps[0];

  S[8] = -2.0*eps[2];
  S[9] = 2.0*eps[1];
  S[10] = -2.0*eps[0];
  S[11] = 2.0*eta;
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
static inline void addDMatTransProduct( const TacsScalar a,
                                        const TacsScalar v[],
                                        const TacsScalar x[],
                                        const TacsScalar eta,
                                        const TacsScalar eps[],
                                        TacsScalar *yeta,
                                        TacsScalar yeps[] ){
  // Compute the cross product
  TacsScalar t[3];
  crossProduct(1.0, v, x, t);
  *yeta -= 2.0*a*vecDot(eps, t);

  // Add the term 2*eps^{x}*v^{x}*x - 2*eta* v^{x}*x
  crossProductAdd(2.0*a, eps, t, yeps);
  crossProductAdd(-2.0*a*eta, v, x, yeps);

  // Compute the final term -4*v^{x}*eps^{x}*x
  crossProduct(1.0, eps, x, t);
  crossProductAdd(-4.0*a, v, t, yeps);
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
static inline void addEMatTransProduct( const TacsScalar a,
                                        const TacsScalar v[],
                                        const TacsScalar x[],
                                        const TacsScalar eta,
                                        const TacsScalar eps[],
                                        TacsScalar *yeta,
                                        TacsScalar yeps[] ){
  // Compute the cross product
  TacsScalar t[3];
  crossProduct(1.0, v, x, t);
  *yeta += 2.0*a*vecDot(eps, t);

  // Add the term 2*eps^{x}*v^{x}*x + 2*eta*v^{x}*x
  crossProductAdd(2.0*a, eps, t, yeps);
  crossProductAdd(2.0*a*eta, v, x, yeps);

  // Compute the final term -4*v^{x}*eps^{x}*x
  crossProduct(1.0, eps, x, t);
  crossProductAdd(-4.0*a, v, t, yeps);
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
static inline void computeDMat( const TacsScalar eta,
                                const TacsScalar eps[],
                                const TacsScalar v[],
                                TacsScalar D[] ){
  D[0] = 2.0*(v[1]*eps[2] - v[2]*eps[1]);
  D[1] = 2.0*(v[1]*eps[1] + v[2]*eps[2]);
  D[2] = 2.0*(eps[0]*v[1] - 2.0*v[0]*eps[1] - eta*v[2]);
  D[3] = 2.0*(eps[0]*v[2] - 2.0*v[0]*eps[2] + eta*v[1]);

  D[4] = 2.0*(v[2]*eps[0] - v[0]*eps[2]);
  D[5] = 2.0*(eps[1]*v[0] - 2.0*v[1]*eps[0] + eta*v[2]);
  D[6] = 2.0*(v[0]*eps[0] + v[2]*eps[2]);
  D[7] = 2.0*(eps[1]*v[2] - 2.0*v[1]*eps[2] - eta*v[0]);

  D[8] = 2.0*(v[0]*eps[1] - v[1]*eps[0]);
  D[9] = 2.0*(eps[2]*v[0] - 2.0*v[2]*eps[0] - eta*v[1]);
  D[10]= 2.0*(eps[2]*v[1] - 2.0*v[2]*eps[1] + eta*v[0]);
  D[11]= 2.0*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void addBlockDMat( const TacsScalar a,
                                 const TacsScalar eta,
                                 const TacsScalar eps[],
                                 const TacsScalar v[],
                                 TacsScalar D[],
                                 const int ldd ){
  const TacsScalar d = 2.0*a;
  D[0] += d*(v[1]*eps[2] - v[2]*eps[1]);
  D[1] += d*(v[1]*eps[1] + v[2]*eps[2]);
  D[2] += d*(eps[0]*v[1] - 2.0*v[0]*eps[1] - eta*v[2]);
  D[3] += d*(eps[0]*v[2] - 2.0*v[0]*eps[2] + eta*v[1]);
  D += ldd;

  D[0] += d*(v[2]*eps[0] - v[0]*eps[2]);
  D[1] += d*(eps[1]*v[0] - 2.0*v[1]*eps[0] + eta*v[2]);
  D[2] += d*(v[0]*eps[0] + v[2]*eps[2]);
  D[3] += d*(eps[1]*v[2] - 2.0*v[1]*eps[2] - eta*v[0]);
  D += ldd;

  D[0] += d*(v[0]*eps[1] - v[1]*eps[0]);
  D[1] += d*(eps[2]*v[0] - 2.0*v[2]*eps[0] - eta*v[1]);
  D[2] += d*(eps[2]*v[1] - 2.0*v[2]*eps[1] + eta*v[0]);
  D[3] += d*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void addBlockDMatTrans( const TacsScalar a,
                                      const TacsScalar eta,
                                      const TacsScalar eps[],
                                      const TacsScalar v[],
                                      TacsScalar D[],
                                      const int ldd ){
  const TacsScalar d = 2.0*a;
  D[0] += d*(v[1]*eps[2] - v[2]*eps[1]);
  D[1] += d*(v[2]*eps[0] - v[0]*eps[2]);
  D[2] += d*(v[0]*eps[1] - v[1]*eps[0]);
  D += ldd;

  D[0] += d*(v[1]*eps[1] + v[2]*eps[2]);
  D[1] += d*(eps[1]*v[0] - 2.0*v[1]*eps[0] + eta*v[2]);
  D[2] += d*(eps[2]*v[0] - 2.0*v[2]*eps[0] - eta*v[1]);
  D += ldd;

  D[0] += d*(eps[0]*v[1] - 2.0*v[0]*eps[1] - eta*v[2]);
  D[1] += d*(v[0]*eps[0] + v[2]*eps[2]);
  D[2] += d*(eps[2]*v[1] - 2.0*v[2]*eps[1] + eta*v[0]);
  D += ldd;

  D[0] += d*(eps[0]*v[2] - 2.0*v[0]*eps[2] + eta*v[1]);
  D[1] += d*(eps[1]*v[2] - 2.0*v[1]*eps[2] - eta*v[0]);
  D[2] += d*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void computeEMat( const TacsScalar eta,
                                const TacsScalar eps[],
                                const TacsScalar v[],
                                TacsScalar E[] ){
  E[0] = 2.0*(v[2]*eps[1] - v[1]*eps[2]);
  E[1] = 2.0*(v[1]*eps[1] + v[2]*eps[2]);
  E[2] = 2.0*(eps[0]*v[1] - 2.0*v[0]*eps[1] + eta*v[2]);
  E[3] = 2.0*(eps[0]*v[2] - 2.0*v[0]*eps[2] - eta*v[1]);

  E[4] = 2.0*(v[0]*eps[2] - v[2]*eps[0]);
  E[5] = 2.0*(eps[1]*v[0] - 2.0*v[1]*eps[0] - eta*v[2]);
  E[6] = 2.0*(v[0]*eps[0] + v[2]*eps[2]);
  E[7] = 2.0*(eps[1]*v[2] - 2.0*v[1]*eps[2] + eta*v[0]);

  E[8] = 2.0*(v[1]*eps[0] - v[0]*eps[1]);
  E[9] = 2.0*(eps[2]*v[0] - 2.0*v[2]*eps[0] + eta*v[1]);
  E[10]= 2.0*(eps[2]*v[1] - 2.0*v[2]*eps[1] - eta*v[0]);
  E[11]= 2.0*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void addBlockDMatTransDeriv( const TacsScalar a,
                                           const TacsScalar v[],
                                           const TacsScalar w[],
                                           TacsScalar D[],
                                           const int ldd ){
  const TacsScalar b = 2.0*a;
  // Compute the cross product w^{x}*v
  TacsScalar t[3];
  crossProduct(1.0, w, v, t);

  D[1] += b*t[0];
  D[2] += b*t[1];
  D[3] += b*t[2];
  D += ldd;

  D[0] += b*t[0];
  D[1] -= 2.0*b*(v[1]*w[1] + v[2]*w[2]);
  D[2] += b*(v[0]*w[1] + w[0]*v[1]);
  D[3] += b*(v[0]*w[2] + w[0]*v[2]);
  D += ldd;

  D[0] += b*t[1];
  D[1] += b*(v[1]*w[0] + w[1]*v[0]);
  D[2] -= 2.0*b*(v[0]*w[0] + v[2]*w[2]);
  D[3] += b*(v[1]*w[2] + w[1]*v[2]);
  D += ldd;

  D[0] += b*t[2];
  D[1] += b*(v[2]*w[0] + w[2]*v[0]);
  D[2] += b*(v[2]*w[1] + w[2]*v[1]);
  D[3] -= 2.0*b*(v[0]*w[0] + v[1]*w[1]);
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
static inline void addBlockEMat( const TacsScalar a,
                                 const TacsScalar eta,
                                 const TacsScalar eps[],
                                 const TacsScalar v[],
                                 TacsScalar D[],
                                 const int ldd ){
  const TacsScalar b = 2.0*a;
  D[0] += b*(v[2]*eps[1] - v[1]*eps[2]);
  D[1] += b*(v[1]*eps[1] + v[2]*eps[2]);
  D[2] += b*(eps[0]*v[1] - 2.0*v[0]*eps[1] + eta*v[2]);
  D[3] += b*(eps[0]*v[2] - 2.0*v[0]*eps[2] - eta*v[1]);
  D += ldd;

  D[0] += b*(v[0]*eps[2] - v[2]*eps[0]);
  D[1] += b*(eps[1]*v[0] - 2.0*v[1]*eps[0] - eta*v[2]);
  D[2] += b*(v[0]*eps[0] + v[2]*eps[2]);
  D[3] += b*(eps[1]*v[2] - 2.0*v[1]*eps[2] + eta*v[0]);
  D += ldd;

  D[0] += b*(v[1]*eps[0] - v[0]*eps[1]);
  D[1] += b*(eps[2]*v[0] - 2.0*v[2]*eps[0] + eta*v[1]);
  D[2] += b*(eps[2]*v[1] - 2.0*v[2]*eps[1] - eta*v[0]);
  D[3] += b*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void addBlockEMatTrans( const TacsScalar a,
                                      const TacsScalar eta,
                                      const TacsScalar eps[],
                                      const TacsScalar v[],
                                      TacsScalar D[],
                                      const int ldd ){
  const TacsScalar b = 2.0*a;
  D[0] += b*(v[2]*eps[1] - v[1]*eps[2]);
  D[1] += b*(v[0]*eps[2] - v[2]*eps[0]);
  D[2] += b*(v[1]*eps[0] - v[0]*eps[1]);
  D += ldd;

  D[0] += b*(v[1]*eps[1] + v[2]*eps[2]);
  D[1] += b*(eps[1]*v[0] - 2.0*v[1]*eps[0] - eta*v[2]);
  D[2] += b*(eps[2]*v[0] - 2.0*v[2]*eps[0] + eta*v[1]);
  D += ldd;

  D[0] += b*(eps[0]*v[1] - 2.0*v[0]*eps[1] + eta*v[2]);
  D[1] += b*(v[0]*eps[0] + v[2]*eps[2]);
  D[2] += b*(eps[2]*v[1] - 2.0*v[2]*eps[1] - eta*v[0]);
  D += ldd;

  D[0] += b*(eps[0]*v[2] - 2.0*v[0]*eps[2] - eta*v[1]);
  D[1] += b*(eps[1]*v[2] - 2.0*v[1]*eps[2] + eta*v[0]);
  D[2] += b*(v[0]*eps[0] + v[1]*eps[1]);
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
static inline void addSRateMatTransDeriv( const TacsScalar a,
                                          const TacsScalar v[],
                                          TacsScalar D[],
                                          const int ldd ){
  const TacsScalar b = 2.0*a;
  D[1] -= b*v[0];
  D[2] -= b*v[1];
  D[3] -= b*v[2];
  D += ldd;

  D[0] += b*v[0];
  D[2] += b*v[2];
  D[3] -= b*v[1];
  D += ldd;

  D[0] += b*v[1];
  D[1] -= b*v[2];
  D[3] += b*v[0];
  D += ldd;

  D[0] += b*v[2];
  D[1] += b*v[1];
  D[2] -= b*v[0];
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
static inline void addBlock3x3x4Product( const TacsScalar A[],
                                         const TacsScalar B[],
                                         TacsScalar D[],
                                         const int ldd ){
  D[0] += A[0]*B[0] + A[1]*B[4] + A[2]*B[8];
  D[1] += A[0]*B[1] + A[1]*B[5] + A[2]*B[9];
  D[2] += A[0]*B[2] + A[1]*B[6] + A[2]*B[10];
  D[3] += A[0]*B[3] + A[1]*B[7] + A[2]*B[11];
  D += ldd;

  D[0] += A[3]*B[0] + A[4]*B[4] + A[5]*B[8];
  D[1] += A[3]*B[1] + A[4]*B[5] + A[5]*B[9];
  D[2] += A[3]*B[2] + A[4]*B[6] + A[5]*B[10];
  D[3] += A[3]*B[3] + A[4]*B[7] + A[5]*B[11];
  D += ldd;

  D[0] += A[6]*B[0] + A[7]*B[4] + A[8]*B[8];
  D[1] += A[6]*B[1] + A[7]*B[5] + A[8]*B[9];
  D[2] += A[6]*B[2] + A[7]*B[6] + A[8]*B[10];
  D[3] += A[6]*B[3] + A[7]*B[7] + A[8]*B[11];
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
static inline void addBlock4x3x3Product( const TacsScalar A[],
                                         const TacsScalar B[],
                                         TacsScalar D[],
                                         const int ldd ){
  D[0] += A[0]*B[0] + A[4]*B[3] + A[8]*B[6];
  D[1] += A[0]*B[1] + A[4]*B[4] + A[8]*B[7];
  D[2] += A[0]*B[2] + A[4]*B[5] + A[8]*B[8];
  D += ldd;

  D[0] += A[1]*B[0] + A[5]*B[3] + A[9]*B[6];
  D[1] += A[1]*B[1] + A[5]*B[4] + A[9]*B[7];
  D[2] += A[1]*B[2] + A[5]*B[5] + A[9]*B[8];
  D += ldd;

  D[0] += A[2]*B[0] + A[6]*B[3] + A[10]*B[6];
  D[1] += A[2]*B[1] + A[6]*B[4] + A[10]*B[7];
  D[2] += A[2]*B[2] + A[6]*B[5] + A[10]*B[8];
  D += ldd;

  D[0] += A[3]*B[0] + A[7]*B[3] + A[11]*B[6];
  D[1] += A[3]*B[1] + A[7]*B[4] + A[11]*B[7];
  D[2] += A[3]*B[2] + A[7]*B[5] + A[11]*B[8];
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
static inline void addBlock3x4Product( const TacsScalar a,
                                       const TacsScalar A[],
                                       const TacsScalar B[],
                                       TacsScalar D[],
                                       const int ldd ){
  D[0] += a*(A[0]*B[0] + A[4]*B[4] + A[8]*B[8]);
  D[1] += a*(A[0]*B[1] + A[4]*B[5] + A[8]*B[9]);
  D[2] += a*(A[0]*B[2] + A[4]*B[6] + A[8]*B[10]);
  D[3] += a*(A[0]*B[3] + A[4]*B[7] + A[8]*B[11]);
  D += ldd;

  D[0] += a*(A[1]*B[0] + A[5]*B[4] + A[9]*B[8]);
  D[1] += a*(A[1]*B[1] + A[5]*B[5] + A[9]*B[9]);
  D[2] += a*(A[1]*B[2] + A[5]*B[6] + A[9]*B[10]);
  D[3] += a*(A[1]*B[3] + A[5]*B[7] + A[9]*B[11]);
  D += ldd;

  D[0] += a*(A[2]*B[0] + A[6]*B[4] + A[10]*B[8]);
  D[1] += a*(A[2]*B[1] + A[6]*B[5] + A[10]*B[9]);
  D[2] += a*(A[2]*B[2] + A[6]*B[6] + A[10]*B[10]);
  D[3] += a*(A[2]*B[3] + A[6]*B[7] + A[10]*B[11]);
  D += ldd;

  D[0] += a*(A[3]*B[0] + A[7]*B[4] + A[11]*B[8]);
  D[1] += a*(A[3]*B[1] + A[7]*B[5] + A[11]*B[9]);
  D[2] += a*(A[3]*B[2] + A[7]*B[6] + A[11]*B[10]);
  D[3] += a*(A[3]*B[3] + A[7]*B[7] + A[11]*B[11]);
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
static inline void computeQtr2ndDeriv( const TacsScalar v[],
                                       TacsScalar dv[] ){
  // Derivatives of eta and eps
  dv[0] = 0.0;
  dv[1] = -2.0*v[2];
  dv[2] = 2.0*v[1];
  dv += 3;

  dv[0] = 2.0*v[2];
  dv[1] = 0.0;
  dv[2] = -2.0*v[0];
  dv += 3;

  dv[0] = -2.0*v[1];
  dv[1] = 2.0*v[0];
  dv[2] = 0.0;
  dv += 3;

  // Second derivatives w.r.t eps
  // C,11
  dv[0] = 0.0;
  dv[1] = -4.0*v[1];
  dv[2] = -4.0*v[2];
  dv += 3;

  // C,12
  dv[0] = 2.0*v[1];
  dv[1] = 2.0*v[0];
  dv[2] = 0.0;
  dv += 3;

  // C,13
  dv[0] = 2.0*v[2];
  dv[1] = 0.0;
  dv[2] = 2.0*v[0];
  dv += 3;

  // C,22
  dv[0] = -4.0*v[0];
  dv[1] = 0.0;
  dv[2] = -4.0*v[2];
  dv += 3;

  // C,23
  dv[0] = 0.0;
  dv[1] = 2.0*v[2];
  dv[2] = 2.0*v[1];
  dv += 3;

  // C,33
  dv[0] = -4.0*v[0];
  dv[1] = -4.0*v[1];
  dv[2] = 0.0;
  dv += 3;
}

#endif // TACS_ALGEBRA_H
