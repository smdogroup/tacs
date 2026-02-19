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

#ifndef TACS_TENSOR_PRODUCT_BASIS_IMPL_H
#define TACS_TENSOR_PRODUCT_BASIS_IMPL_H

#include "TACSObject.h"

/*
  Inline functions for tensor-product basis computations
*/
inline void TacsInterpTensor3DInterp2(const int m, const double *n1,
                                      const double *n2, const double *n3,
                                      const TacsScalar *v, TacsScalar *u) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m]);
  TacsScalar t2 = (n1[0] * v[2 * m] + n1[1] * v[3 * m]);
  TacsScalar t3 = (n1[0] * v[4 * m] + n1[1] * v[5 * m]);
  TacsScalar t4 = (n1[0] * v[6 * m] + n1[1] * v[7 * m]);

  u[0] = n3[0] * (n2[0] * t1 + n2[1] * t2) + n3[1] * (n2[0] * t3 + n2[1] * t4);
}

inline void TacsAddTransTensor3DInterp2(const int m, const double *n1,
                                        const double *n2, const double *n3,
                                        const TacsScalar *u, TacsScalar *v) {
  TacsScalar a1 = n2[0] * u[0];
  TacsScalar a2 = n2[1] * u[0];

  TacsScalar t1 = n3[0] * a1;
  TacsScalar t2 = n3[0] * a2;
  TacsScalar t3 = n3[1] * a1;
  TacsScalar t4 = n3[1] * a2;

  v[0] += n1[0] * t1;
  v[m] += n1[1] * t1;

  v[2 * m] += n1[0] * t2;
  v[3 * m] += n1[1] * t2;

  v[4 * m] += n1[0] * t3;
  v[5 * m] += n1[1] * t3;

  v[6 * m] += n1[0] * t4;
  v[7 * m] += n1[1] * t4;
}

inline void TacsGradTensor3DInterp2(const int m, const double *n1,
                                    const double *n2, const double *n3,
                                    const double *n1x, const double *n2x,
                                    const double *n3x, const TacsScalar *v,
                                    TacsScalar *g) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m]);
  TacsScalar t2 = (n1[0] * v[2 * m] + n1[1] * v[3 * m]);
  TacsScalar t3 = (n1[0] * v[4 * m] + n1[1] * v[5 * m]);
  TacsScalar t4 = (n1[0] * v[6 * m] + n1[1] * v[7 * m]);

  TacsScalar t1x = (n1x[0] * v[0] + n1x[1] * v[m]);
  TacsScalar t2x = (n1x[0] * v[2 * m] + n1x[1] * v[3 * m]);
  TacsScalar t3x = (n1x[0] * v[4 * m] + n1x[1] * v[5 * m]);
  TacsScalar t4x = (n1x[0] * v[6 * m] + n1x[1] * v[7 * m]);

  g[0] =
      n3[0] * (n2[0] * t1x + n2[1] * t2x) + n3[1] * (n2[0] * t3x + n2[1] * t4x);
  g[1] =
      n3[0] * (n2x[0] * t1 + n2x[1] * t2) + n3[1] * (n2x[0] * t3 + n2x[1] * t4);
  g[2] =
      n3x[0] * (n2[0] * t1 + n2[1] * t2) + n3x[1] * (n2[0] * t3 + n2[1] * t4);
}

inline void TacsAddGradTransTensor3DInterp2(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1];
  TacsScalar b2 = n2x[1] * g[1];
  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];

  TacsScalar t1 = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2 = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3 = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t4 = n3[1] * b2 + n3x[1] * c2;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];

  TacsScalar t1x = n3[0] * a1;
  TacsScalar t2x = n3[0] * a2;
  TacsScalar t3x = n3[1] * a1;
  TacsScalar t4x = n3[1] * a2;

  v[0] += n1[0] * t1 + n1x[0] * t1x;
  v[m] += n1[1] * t1 + n1x[1] * t1x;
  v[2 * m] += n1[0] * t2 + n1x[0] * t2x;
  v[3 * m] += n1[1] * t2 + n1x[1] * t2x;
  v[4 * m] += n1[0] * t3 + n1x[0] * t3x;
  v[5 * m] += n1[1] * t3 + n1x[1] * t3x;
  v[6 * m] += n1[0] * t4 + n1x[0] * t4x;
  v[7 * m] += n1[1] * t4 + n1x[1] * t4x;
}

inline void TacsInterpGradTensor3DInterp2(const int m, const double *n1,
                                          const double *n2, const double *n3,
                                          const double *n1x, const double *n2x,
                                          const double *n3x,
                                          const TacsScalar *v, TacsScalar *u,
                                          TacsScalar *g) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m]);
  TacsScalar t2 = (n1[0] * v[2 * m] + n1[1] * v[3 * m]);
  TacsScalar t3 = (n1[0] * v[4 * m] + n1[1] * v[5 * m]);
  TacsScalar t4 = (n1[0] * v[6 * m] + n1[1] * v[7 * m]);

  TacsScalar t1x = (n1x[0] * v[0] + n1x[1] * v[m]);
  TacsScalar t2x = (n1x[0] * v[2 * m] + n1x[1] * v[3 * m]);
  TacsScalar t3x = (n1x[0] * v[4 * m] + n1x[1] * v[5 * m]);
  TacsScalar t4x = (n1x[0] * v[6 * m] + n1x[1] * v[7 * m]);

  TacsScalar a1 = (n2[0] * t1 + n2[1] * t2);
  TacsScalar a2 = (n2[0] * t3 + n2[1] * t4);

  u[0] = n3[0] * a1 + n3[1] * a2;
  g[0] =
      n3[0] * (n2[0] * t1x + n2[1] * t2x) + n3[1] * (n2[0] * t3x + n2[1] * t4x);
  g[1] =
      n3[0] * (n2x[0] * t1 + n2x[1] * t2) + n3[1] * (n2x[0] * t3 + n2x[1] * t4);
  g[2] = n3x[0] * a1 + n3x[1] * a2;
}

inline void TacsAddInterpGradTransTensor3DInterp2(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *u, const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1] + n2[0] * u[0];
  TacsScalar b2 = n2x[1] * g[1] + n2[1] * u[0];
  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];

  TacsScalar t1 = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2 = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3 = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t4 = n3[1] * b2 + n3x[1] * c2;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];

  TacsScalar t1x = n3[0] * a1;
  TacsScalar t2x = n3[0] * a2;
  TacsScalar t3x = n3[1] * a1;
  TacsScalar t4x = n3[1] * a2;

  v[0] += n1[0] * t1 + n1x[0] * t1x;
  v[m] += n1[1] * t1 + n1x[1] * t1x;
  v[2 * m] += n1[0] * t2 + n1x[0] * t2x;
  v[3 * m] += n1[1] * t2 + n1x[1] * t2x;
  v[4 * m] += n1[0] * t3 + n1x[0] * t3x;
  v[5 * m] += n1[1] * t3 + n1x[1] * t3x;
  v[6 * m] += n1[0] * t4 + n1x[0] * t4x;
  v[7 * m] += n1[1] * t4 + n1x[1] * t4x;
}

/*
  3D Tensor product functions for p = 1
*/
inline void TACSInterpAllTensor3DInterp2(const int m, const int i,
                                         const double N[], const double Nx[],
                                         const TacsScalar v[],
                                         TacsScalar out[]) {
  const int j = 3 * i + m;
  TacsInterpGradTensor3DInterp2(m, &N[0], &N[0], &N[0], &Nx[0], &Nx[0], &Nx[0],
                                v, &out[i], &out[j]);
  TacsInterpGradTensor3DInterp2(m, &N[2], &N[0], &N[0], &Nx[2], &Nx[0], &Nx[0],
                                v, &out[4 * m + i], &out[4 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[0], &N[2], &N[0], &Nx[0], &Nx[2], &Nx[0],
                                v, &out[8 * m + i], &out[8 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[2], &N[2], &N[0], &Nx[2], &Nx[2], &Nx[0],
                                v, &out[12 * m + i], &out[12 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[0], &N[0], &N[2], &Nx[0], &Nx[0], &Nx[2],
                                v, &out[16 * m + i], &out[16 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[2], &N[0], &N[2], &Nx[2], &Nx[0], &Nx[2],
                                v, &out[20 * m + i], &out[20 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[0], &N[2], &N[2], &Nx[0], &Nx[2], &Nx[2],
                                v, &out[24 * m + i], &out[24 * m + j]);
  TacsInterpGradTensor3DInterp2(m, &N[2], &N[2], &N[2], &Nx[2], &Nx[2], &Nx[2],
                                v, &out[28 * m + i], &out[28 * m + j]);
}

inline void TacsAddAllTransTensor3DInterp2(const int m, const int i,
                                           const double N[], const double Nx[],
                                           const TacsScalar in[],
                                           TacsScalar v[]) {
  const int j = 3 * i + m;
  TacsAddInterpGradTransTensor3DInterp2(m, &N[0], &N[0], &N[0], &Nx[0], &Nx[0],
                                        &Nx[0], &in[i], &in[j], v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[2], &N[0], &N[0], &Nx[2], &Nx[0],
                                        &Nx[0], &in[4 * m + i], &in[4 * m + j],
                                        v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[0], &N[2], &N[0], &Nx[0], &Nx[2],
                                        &Nx[0], &in[8 * m + i], &in[8 * m + j],
                                        v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[2], &N[2], &N[0], &Nx[2], &Nx[2],
                                        &Nx[0], &in[12 * m + i],
                                        &in[12 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[0], &N[0], &N[2], &Nx[0], &Nx[0],
                                        &Nx[2], &in[16 * m + i],
                                        &in[16 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[2], &N[0], &N[2], &Nx[2], &Nx[0],
                                        &Nx[2], &in[20 * m + i],
                                        &in[20 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[0], &N[2], &N[2], &Nx[0], &Nx[2],
                                        &Nx[2], &in[24 * m + i],
                                        &in[24 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp2(m, &N[2], &N[2], &N[2], &Nx[2], &Nx[2],
                                        &Nx[2], &in[28 * m + i],
                                        &in[28 * m + j], v);
}

/*
  3D Tensor product functions for p = 2
*/
inline void TacsInterpTensor3DInterp3(const int m, const double *n1,
                                      const double *n2, const double *n3,
                                      const TacsScalar *v, TacsScalar *u) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m]);
  TacsScalar t2 = (n1[0] * v[3 * m] + n1[1] * v[4 * m] + n1[2] * v[5 * m]);
  TacsScalar t3 = (n1[0] * v[6 * m] + n1[1] * v[7 * m] + n1[2] * v[8 * m]);
  TacsScalar t4 = (n1[0] * v[9 * m] + n1[1] * v[10 * m] + n1[2] * v[11 * m]);
  TacsScalar t5 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m]);
  TacsScalar t6 = (n1[0] * v[15 * m] + n1[1] * v[16 * m] + n1[2] * v[17 * m]);
  TacsScalar t7 = (n1[0] * v[18 * m] + n1[1] * v[19 * m] + n1[2] * v[20 * m]);
  TacsScalar t8 = (n1[0] * v[21 * m] + n1[1] * v[22 * m] + n1[2] * v[23 * m]);
  TacsScalar t9 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m]);

  u[0] = n3[0] * (n2[0] * t1 + n2[1] * t2 + n2[2] * t3) +
         n3[1] * (n2[0] * t4 + n2[1] * t5 + n2[2] * t6) +
         n3[2] * (n2[0] * t7 + n2[1] * t8 + n2[2] * t9);
}

inline void TacsAddTransTensor3DInterp3(const int m, const double *n1,
                                        const double *n2, const double *n3,
                                        const TacsScalar *u, TacsScalar *v) {
  TacsScalar a1 = n2[0] * u[0];
  TacsScalar a2 = n2[1] * u[0];
  TacsScalar a3 = n2[2] * u[0];
  TacsScalar t1 = n3[0] * a1;
  TacsScalar t2 = n3[0] * a2;
  TacsScalar t3 = n3[0] * a3;
  TacsScalar t4 = n3[1] * a1;
  TacsScalar t5 = n3[1] * a2;
  TacsScalar t6 = n3[1] * a3;
  TacsScalar t7 = n3[2] * a1;
  TacsScalar t8 = n3[2] * a2;
  TacsScalar t9 = n3[2] * a3;

  v[0] += n1[0] * t1;
  v[m] += n1[1] * t1;
  v[2 * m] += n1[2] * t1;
  v[3 * m] += n1[0] * t2;
  v[4 * m] += n1[1] * t2;
  v[5 * m] += n1[2] * t2;
  v[6 * m] += n1[0] * t3;
  v[7 * m] += n1[1] * t3;
  v[8 * m] += n1[2] * t3;
  v[9 * m] += n1[0] * t4;
  v[10 * m] += n1[1] * t4;
  v[11 * m] += n1[2] * t4;
  v[12 * m] += n1[0] * t5;
  v[13 * m] += n1[1] * t5;
  v[14 * m] += n1[2] * t5;
  v[15 * m] += n1[0] * t6;
  v[16 * m] += n1[1] * t6;
  v[17 * m] += n1[2] * t6;
  v[18 * m] += n1[0] * t7;
  v[19 * m] += n1[1] * t7;
  v[20 * m] += n1[2] * t7;
  v[21 * m] += n1[0] * t8;
  v[22 * m] += n1[1] * t8;
  v[23 * m] += n1[2] * t8;
  v[24 * m] += n1[0] * t9;
  v[25 * m] += n1[1] * t9;
  v[26 * m] += n1[2] * t9;
}

inline void TacsGradTensor3DInterp3(const int m, const double *n1,
                                    const double *n2, const double *n3,
                                    const double *n1x, const double *n2x,
                                    const double *n3x, const TacsScalar *v,
                                    TacsScalar *g) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m]);
  TacsScalar t2 = (n1[0] * v[3 * m] + n1[1] * v[4 * m] + n1[2] * v[5 * m]);
  TacsScalar t3 = (n1[0] * v[6 * m] + n1[1] * v[7 * m] + n1[2] * v[8 * m]);
  TacsScalar t4 = (n1[0] * v[9 * m] + n1[1] * v[10 * m] + n1[2] * v[11 * m]);
  TacsScalar t5 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m]);
  TacsScalar t6 = (n1[0] * v[15 * m] + n1[1] * v[16 * m] + n1[2] * v[17 * m]);
  TacsScalar t7 = (n1[0] * v[18 * m] + n1[1] * v[19 * m] + n1[2] * v[20 * m]);
  TacsScalar t8 = (n1[0] * v[21 * m] + n1[1] * v[22 * m] + n1[2] * v[23 * m]);
  TacsScalar t9 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m]);

  TacsScalar t1x = (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m]);
  TacsScalar t2x = (n1x[0] * v[3 * m] + n1x[1] * v[4 * m] + n1x[2] * v[5 * m]);
  TacsScalar t3x = (n1x[0] * v[6 * m] + n1x[1] * v[7 * m] + n1x[2] * v[8 * m]);
  TacsScalar t4x =
      (n1x[0] * v[9 * m] + n1x[1] * v[10 * m] + n1x[2] * v[11 * m]);
  TacsScalar t5x =
      (n1x[0] * v[12 * m] + n1x[1] * v[13 * m] + n1x[2] * v[14 * m]);
  TacsScalar t6x =
      (n1x[0] * v[15 * m] + n1x[1] * v[16 * m] + n1x[2] * v[17 * m]);
  TacsScalar t7x =
      (n1x[0] * v[18 * m] + n1x[1] * v[19 * m] + n1x[2] * v[20 * m]);
  TacsScalar t8x =
      (n1x[0] * v[21 * m] + n1x[1] * v[22 * m] + n1x[2] * v[23 * m]);
  TacsScalar t9x =
      (n1x[0] * v[24 * m] + n1x[1] * v[25 * m] + n1x[2] * v[26 * m]);

  g[0] = n3[0] * (n2[0] * t1x + n2[1] * t2x + n2[2] * t3x) +
         n3[1] * (n2[0] * t4x + n2[1] * t5x + n2[2] * t6x) +
         n3[2] * (n2[0] * t7x + n2[1] * t8x + n2[2] * t9x);

  g[1] = n3[0] * (n2x[0] * t1 + n2x[1] * t2 + n2x[2] * t3) +
         n3[1] * (n2x[0] * t4 + n2x[1] * t5 + n2x[2] * t6) +
         n3[2] * (n2x[0] * t7 + n2x[1] * t8 + n2x[2] * t9);

  g[2] = n3x[0] * (n2[0] * t1 + n2[1] * t2 + n2[2] * t3) +
         n3x[1] * (n2[0] * t4 + n2[1] * t5 + n2[2] * t6) +
         n3x[2] * (n2[0] * t7 + n2[1] * t8 + n2[2] * t9);
}

inline void TacsAddGradTransTensor3DInterp3(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1];
  TacsScalar b2 = n2x[1] * g[1];
  TacsScalar b3 = n2x[2] * g[1];

  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];
  TacsScalar c3 = n2[2] * g[2];

  TacsScalar t1 = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2 = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3 = n3[0] * b3 + n3x[0] * c3;
  TacsScalar t4 = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t5 = n3[1] * b2 + n3x[1] * c2;
  TacsScalar t6 = n3[1] * b3 + n3x[1] * c3;
  TacsScalar t7 = n3[2] * b1 + n3x[2] * c1;
  TacsScalar t8 = n3[2] * b2 + n3x[2] * c2;
  TacsScalar t9 = n3[2] * b3 + n3x[2] * c3;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];
  TacsScalar a3 = n2[2] * g[0];

  TacsScalar t1x = n3[0] * a1;
  TacsScalar t2x = n3[0] * a2;
  TacsScalar t3x = n3[0] * a3;
  TacsScalar t4x = n3[1] * a1;
  TacsScalar t5x = n3[1] * a2;
  TacsScalar t6x = n3[1] * a3;
  TacsScalar t7x = n3[2] * a1;
  TacsScalar t8x = n3[2] * a2;
  TacsScalar t9x = n3[2] * a3;

  v[0] += n1[0] * t1 + n1x[0] * t1x;
  v[m] += n1[1] * t1 + n1x[1] * t1x;
  v[2 * m] += n1[2] * t1 + n1x[2] * t1x;
  v[3 * m] += n1[0] * t2 + n1x[0] * t2x;
  v[4 * m] += n1[1] * t2 + n1x[1] * t2x;
  v[5 * m] += n1[2] * t2 + n1x[2] * t2x;
  v[6 * m] += n1[0] * t3 + n1x[0] * t3x;
  v[7 * m] += n1[1] * t3 + n1x[1] * t3x;
  v[8 * m] += n1[2] * t3 + n1x[2] * t3x;
  v[9 * m] += n1[0] * t4 + n1x[0] * t4x;
  v[10 * m] += n1[1] * t4 + n1x[1] * t4x;
  v[11 * m] += n1[2] * t4 + n1x[2] * t4x;
  v[12 * m] += n1[0] * t5 + n1x[0] * t5x;
  v[13 * m] += n1[1] * t5 + n1x[1] * t5x;
  v[14 * m] += n1[2] * t5 + n1x[2] * t5x;
  v[15 * m] += n1[0] * t6 + n1x[0] * t6x;
  v[16 * m] += n1[1] * t6 + n1x[1] * t6x;
  v[17 * m] += n1[2] * t6 + n1x[2] * t6x;
  v[18 * m] += n1[0] * t7 + n1x[0] * t7x;
  v[19 * m] += n1[1] * t7 + n1x[1] * t7x;
  v[20 * m] += n1[2] * t7 + n1x[2] * t7x;
  v[21 * m] += n1[0] * t8 + n1x[0] * t8x;
  v[22 * m] += n1[1] * t8 + n1x[1] * t8x;
  v[23 * m] += n1[2] * t8 + n1x[2] * t8x;
  v[24 * m] += n1[0] * t9 + n1x[0] * t9x;
  v[25 * m] += n1[1] * t9 + n1x[1] * t9x;
  v[26 * m] += n1[2] * t9 + n1x[2] * t9x;
}

inline void TacsInterpGradTensor3DInterp3(const int m, const double *n1,
                                          const double *n2, const double *n3,
                                          const double *n1x, const double *n2x,
                                          const double *n3x,
                                          const TacsScalar *v, TacsScalar *u,
                                          TacsScalar *g) {
  TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m]);
  TacsScalar t2 = (n1[0] * v[3 * m] + n1[1] * v[4 * m] + n1[2] * v[5 * m]);
  TacsScalar t3 = (n1[0] * v[6 * m] + n1[1] * v[7 * m] + n1[2] * v[8 * m]);
  TacsScalar t4 = (n1[0] * v[9 * m] + n1[1] * v[10 * m] + n1[2] * v[11 * m]);
  TacsScalar t5 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m]);
  TacsScalar t6 = (n1[0] * v[15 * m] + n1[1] * v[16 * m] + n1[2] * v[17 * m]);
  TacsScalar t7 = (n1[0] * v[18 * m] + n1[1] * v[19 * m] + n1[2] * v[20 * m]);
  TacsScalar t8 = (n1[0] * v[21 * m] + n1[1] * v[22 * m] + n1[2] * v[23 * m]);
  TacsScalar t9 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m]);

  TacsScalar t1x = (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m]);
  TacsScalar t2x = (n1x[0] * v[3 * m] + n1x[1] * v[4 * m] + n1x[2] * v[5 * m]);
  TacsScalar t3x = (n1x[0] * v[6 * m] + n1x[1] * v[7 * m] + n1x[2] * v[8 * m]);
  TacsScalar t4x =
      (n1x[0] * v[9 * m] + n1x[1] * v[10 * m] + n1x[2] * v[11 * m]);
  TacsScalar t5x =
      (n1x[0] * v[12 * m] + n1x[1] * v[13 * m] + n1x[2] * v[14 * m]);
  TacsScalar t6x =
      (n1x[0] * v[15 * m] + n1x[1] * v[16 * m] + n1x[2] * v[17 * m]);
  TacsScalar t7x =
      (n1x[0] * v[18 * m] + n1x[1] * v[19 * m] + n1x[2] * v[20 * m]);
  TacsScalar t8x =
      (n1x[0] * v[21 * m] + n1x[1] * v[22 * m] + n1x[2] * v[23 * m]);
  TacsScalar t9x =
      (n1x[0] * v[24 * m] + n1x[1] * v[25 * m] + n1x[2] * v[26 * m]);

  TacsScalar a1 = (n2[0] * t1 + n2[1] * t2 + n2[2] * t3);
  TacsScalar a2 = (n2[0] * t4 + n2[1] * t5 + n2[2] * t6);
  TacsScalar a3 = (n2[0] * t7 + n2[1] * t8 + n2[2] * t9);

  u[0] = n3[0] * a1 + n3[1] * a2 + n3[2] * a3;

  g[0] = n3[0] * (n2[0] * t1x + n2[1] * t2x + n2[2] * t3x) +
         n3[1] * (n2[0] * t4x + n2[1] * t5x + n2[2] * t6x) +
         n3[2] * (n2[0] * t7x + n2[1] * t8x + n2[2] * t9x);

  g[1] = n3[0] * (n2x[0] * t1 + n2x[1] * t2 + n2x[2] * t3) +
         n3[1] * (n2x[0] * t4 + n2x[1] * t5 + n2x[2] * t6) +
         n3[2] * (n2x[0] * t7 + n2x[1] * t8 + n2x[2] * t9);

  g[2] = n3x[0] * a1 + n3x[1] * a2 + n3x[2] * a3;
}

inline void TacsAddInterpGradTransTensor3DInterp3(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *u, const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1] + n2[0] * u[0];
  TacsScalar b2 = n2x[1] * g[1] + n2[1] * u[0];
  TacsScalar b3 = n2x[2] * g[1] + n2[2] * u[0];

  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];
  TacsScalar c3 = n2[2] * g[2];

  TacsScalar t1 = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2 = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3 = n3[0] * b3 + n3x[0] * c3;
  TacsScalar t4 = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t5 = n3[1] * b2 + n3x[1] * c2;
  TacsScalar t6 = n3[1] * b3 + n3x[1] * c3;
  TacsScalar t7 = n3[2] * b1 + n3x[2] * c1;
  TacsScalar t8 = n3[2] * b2 + n3x[2] * c2;
  TacsScalar t9 = n3[2] * b3 + n3x[2] * c3;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];
  TacsScalar a3 = n2[2] * g[0];

  TacsScalar t1x = n3[0] * a1;
  TacsScalar t2x = n3[0] * a2;
  TacsScalar t3x = n3[0] * a3;
  TacsScalar t4x = n3[1] * a1;
  TacsScalar t5x = n3[1] * a2;
  TacsScalar t6x = n3[1] * a3;
  TacsScalar t7x = n3[2] * a1;
  TacsScalar t8x = n3[2] * a2;
  TacsScalar t9x = n3[2] * a3;

  v[0] += n1[0] * t1 + n1x[0] * t1x;
  v[m] += n1[1] * t1 + n1x[1] * t1x;
  v[2 * m] += n1[2] * t1 + n1x[2] * t1x;
  v[3 * m] += n1[0] * t2 + n1x[0] * t2x;
  v[4 * m] += n1[1] * t2 + n1x[1] * t2x;
  v[5 * m] += n1[2] * t2 + n1x[2] * t2x;
  v[6 * m] += n1[0] * t3 + n1x[0] * t3x;
  v[7 * m] += n1[1] * t3 + n1x[1] * t3x;
  v[8 * m] += n1[2] * t3 + n1x[2] * t3x;
  v[9 * m] += n1[0] * t4 + n1x[0] * t4x;
  v[10 * m] += n1[1] * t4 + n1x[1] * t4x;
  v[11 * m] += n1[2] * t4 + n1x[2] * t4x;
  v[12 * m] += n1[0] * t5 + n1x[0] * t5x;
  v[13 * m] += n1[1] * t5 + n1x[1] * t5x;
  v[14 * m] += n1[2] * t5 + n1x[2] * t5x;
  v[15 * m] += n1[0] * t6 + n1x[0] * t6x;
  v[16 * m] += n1[1] * t6 + n1x[1] * t6x;
  v[17 * m] += n1[2] * t6 + n1x[2] * t6x;
  v[18 * m] += n1[0] * t7 + n1x[0] * t7x;
  v[19 * m] += n1[1] * t7 + n1x[1] * t7x;
  v[20 * m] += n1[2] * t7 + n1x[2] * t7x;
  v[21 * m] += n1[0] * t8 + n1x[0] * t8x;
  v[22 * m] += n1[1] * t8 + n1x[1] * t8x;
  v[23 * m] += n1[2] * t8 + n1x[2] * t8x;
  v[24 * m] += n1[0] * t9 + n1x[0] * t9x;
  v[25 * m] += n1[1] * t9 + n1x[1] * t9x;
  v[26 * m] += n1[2] * t9 + n1x[2] * t9x;
}

inline void TACSInterpAllTensor3DInterp3(const int m, const int i,
                                         const double N[], const double Nx[],
                                         const TacsScalar v[],
                                         TacsScalar out[]) {
  const int j = 3 * i + m;
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[0], &N[0], &Nx[0], &Nx[0], &Nx[0],
                                v, &out[i], &out[j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[0], &N[0], &Nx[3], &Nx[0], &Nx[0],
                                v, &out[4 * m + i], &out[4 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[0], &N[0], &Nx[6], &Nx[0], &Nx[0],
                                v, &out[8 * m + i], &out[8 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[3], &N[0], &Nx[0], &Nx[3], &Nx[0],
                                v, &out[12 * m + i], &out[12 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[3], &N[0], &Nx[3], &Nx[3], &Nx[0],
                                v, &out[16 * m + i], &out[16 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[3], &N[0], &Nx[6], &Nx[3], &Nx[0],
                                v, &out[20 * m + i], &out[20 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[6], &N[0], &Nx[0], &Nx[6], &Nx[0],
                                v, &out[24 * m + i], &out[24 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[6], &N[0], &Nx[3], &Nx[6], &Nx[0],
                                v, &out[28 * m + i], &out[28 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[6], &N[0], &Nx[6], &Nx[6], &Nx[0],
                                v, &out[32 * m + i], &out[32 * m + j]);

  TacsInterpGradTensor3DInterp3(m, &N[0], &N[0], &N[3], &Nx[0], &Nx[0], &Nx[3],
                                v, &out[36 * m + i], &out[36 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[0], &N[3], &Nx[3], &Nx[0], &Nx[3],
                                v, &out[40 * m + i], &out[40 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[0], &N[3], &Nx[6], &Nx[0], &Nx[3],
                                v, &out[44 * m + i], &out[44 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[3], &N[3], &Nx[0], &Nx[3], &Nx[3],
                                v, &out[48 * m + i], &out[48 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[3], &N[3], &Nx[3], &Nx[3], &Nx[3],
                                v, &out[52 * m + i], &out[52 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[3], &N[3], &Nx[6], &Nx[3], &Nx[3],
                                v, &out[56 * m + i], &out[56 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[6], &N[3], &Nx[0], &Nx[6], &Nx[3],
                                v, &out[60 * m + i], &out[60 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[6], &N[3], &Nx[3], &Nx[6], &Nx[3],
                                v, &out[64 * m + i], &out[64 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[6], &N[3], &Nx[6], &Nx[6], &Nx[3],
                                v, &out[68 * m + i], &out[68 * m + j]);

  TacsInterpGradTensor3DInterp3(m, &N[0], &N[0], &N[6], &Nx[0], &Nx[0], &Nx[6],
                                v, &out[72 * m + i], &out[72 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[0], &N[6], &Nx[3], &Nx[0], &Nx[6],
                                v, &out[76 * m + i], &out[76 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[0], &N[6], &Nx[6], &Nx[0], &Nx[6],
                                v, &out[80 * m + i], &out[80 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[3], &N[6], &Nx[0], &Nx[3], &Nx[6],
                                v, &out[84 * m + i], &out[84 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[3], &N[6], &Nx[3], &Nx[3], &Nx[6],
                                v, &out[88 * m + i], &out[88 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[3], &N[6], &Nx[6], &Nx[3], &Nx[6],
                                v, &out[92 * m + i], &out[92 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[0], &N[6], &N[6], &Nx[0], &Nx[6], &Nx[6],
                                v, &out[96 * m + i], &out[96 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[3], &N[6], &N[6], &Nx[3], &Nx[6], &Nx[6],
                                v, &out[100 * m + i], &out[100 * m + j]);
  TacsInterpGradTensor3DInterp3(m, &N[6], &N[6], &N[6], &Nx[6], &Nx[6], &Nx[6],
                                v, &out[104 * m + i], &out[104 * m + j]);
}

inline void TacsAddAllTransTensor3DInterp3(const int m, const int i,
                                           const double N[], const double Nx[],
                                           const TacsScalar in[],
                                           TacsScalar v[]) {
  const int j = 3 * i + m;
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[0], &N[0], &Nx[0], &Nx[0],
                                        &Nx[0], &in[i], &in[j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[0], &N[0], &Nx[3], &Nx[0],
                                        &Nx[0], &in[4 * m + i], &in[4 * m + j],
                                        v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[0], &N[0], &Nx[6], &Nx[0],
                                        &Nx[0], &in[8 * m + i], &in[8 * m + j],
                                        v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[3], &N[0], &Nx[0], &Nx[3],
                                        &Nx[0], &in[12 * m + i],
                                        &in[12 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[3], &N[0], &Nx[3], &Nx[3],
                                        &Nx[0], &in[16 * m + i],
                                        &in[16 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[3], &N[0], &Nx[6], &Nx[3],
                                        &Nx[0], &in[20 * m + i],
                                        &in[20 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[6], &N[0], &Nx[0], &Nx[6],
                                        &Nx[0], &in[24 * m + i],
                                        &in[24 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[6], &N[0], &Nx[3], &Nx[6],
                                        &Nx[0], &in[28 * m + i],
                                        &in[28 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[6], &N[0], &Nx[6], &Nx[6],
                                        &Nx[0], &in[32 * m + i],
                                        &in[32 * m + j], v);

  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[0], &N[3], &Nx[0], &Nx[0],
                                        &Nx[3], &in[36 * m + i],
                                        &in[36 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[0], &N[3], &Nx[3], &Nx[0],
                                        &Nx[3], &in[40 * m + i],
                                        &in[40 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[0], &N[3], &Nx[6], &Nx[0],
                                        &Nx[3], &in[44 * m + i],
                                        &in[44 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[3], &N[3], &Nx[0], &Nx[3],
                                        &Nx[3], &in[48 * m + i],
                                        &in[48 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[3], &N[3], &Nx[3], &Nx[3],
                                        &Nx[3], &in[52 * m + i],
                                        &in[52 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[3], &N[3], &Nx[6], &Nx[3],
                                        &Nx[3], &in[56 * m + i],
                                        &in[56 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[6], &N[3], &Nx[0], &Nx[6],
                                        &Nx[3], &in[60 * m + i],
                                        &in[60 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[6], &N[3], &Nx[3], &Nx[6],
                                        &Nx[3], &in[64 * m + i],
                                        &in[64 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[6], &N[3], &Nx[6], &Nx[6],
                                        &Nx[3], &in[68 * m + i],
                                        &in[68 * m + j], v);

  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[0], &N[6], &Nx[0], &Nx[0],
                                        &Nx[6], &in[72 * m + i],
                                        &in[72 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[0], &N[6], &Nx[3], &Nx[0],
                                        &Nx[6], &in[76 * m + i],
                                        &in[76 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[0], &N[6], &Nx[6], &Nx[0],
                                        &Nx[6], &in[80 * m + i],
                                        &in[80 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[3], &N[6], &Nx[0], &Nx[3],
                                        &Nx[6], &in[84 * m + i],
                                        &in[84 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[3], &N[6], &Nx[3], &Nx[3],
                                        &Nx[6], &in[88 * m + i],
                                        &in[88 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[3], &N[6], &Nx[6], &Nx[3],
                                        &Nx[6], &in[92 * m + i],
                                        &in[92 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[0], &N[6], &N[6], &Nx[0], &Nx[6],
                                        &Nx[6], &in[96 * m + i],
                                        &in[96 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[3], &N[6], &N[6], &Nx[3], &Nx[6],
                                        &Nx[6], &in[100 * m + i],
                                        &in[100 * m + j], v);
  TacsAddInterpGradTransTensor3DInterp3(m, &N[6], &N[6], &N[6], &Nx[6], &Nx[6],
                                        &Nx[6], &in[104 * m + i],
                                        &in[104 * m + j], v);
}

/*
  3D Tensor product functions for p = 3
*/
inline void TacsInterpTensor3DInterp4(const int m, const double *n1,
                                      const double *n2, const double *n3,
                                      const TacsScalar *v, TacsScalar *u) {
  TacsScalar t1 =
      (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] + n1[3] * v[3 * m]);
  TacsScalar t2 = (n1[0] * v[4 * m] + n1[1] * v[5 * m] + n1[2] * v[6 * m] +
                   n1[3] * v[7 * m]);
  TacsScalar t3 = (n1[0] * v[8 * m] + n1[1] * v[9 * m] + n1[2] * v[10 * m] +
                   n1[3] * v[11 * m]);
  TacsScalar t4 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m] +
                   n1[3] * v[15 * m]);
  TacsScalar t5 = (n1[0] * v[16 * m] + n1[1] * v[17 * m] + n1[2] * v[18 * m] +
                   n1[3] * v[19 * m]);
  TacsScalar t6 = (n1[0] * v[20 * m] + n1[1] * v[21 * m] + n1[2] * v[22 * m] +
                   n1[3] * v[23 * m]);
  TacsScalar t7 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m] +
                   n1[3] * v[27 * m]);
  TacsScalar t8 = (n1[0] * v[28 * m] + n1[1] * v[29 * m] + n1[2] * v[30 * m] +
                   n1[3] * v[31 * m]);
  TacsScalar t9 = (n1[0] * v[32 * m] + n1[1] * v[33 * m] + n1[2] * v[34 * m] +
                   n1[3] * v[35 * m]);
  TacsScalar t10 = (n1[0] * v[36 * m] + n1[1] * v[37 * m] + n1[2] * v[38 * m] +
                    n1[3] * v[39 * m]);
  TacsScalar t11 = (n1[0] * v[40 * m] + n1[1] * v[41 * m] + n1[2] * v[42 * m] +
                    n1[3] * v[43 * m]);
  TacsScalar t12 = (n1[0] * v[44 * m] + n1[1] * v[45 * m] + n1[2] * v[46 * m] +
                    n1[3] * v[47 * m]);
  TacsScalar t13 = (n1[0] * v[48 * m] + n1[1] * v[49 * m] + n1[2] * v[50 * m] +
                    n1[3] * v[51 * m]);
  TacsScalar t14 = (n1[0] * v[52 * m] + n1[1] * v[53 * m] + n1[2] * v[54 * m] +
                    n1[3] * v[55 * m]);
  TacsScalar t15 = (n1[0] * v[56 * m] + n1[1] * v[57 * m] + n1[2] * v[58 * m] +
                    n1[3] * v[59 * m]);
  TacsScalar t16 = (n1[0] * v[60 * m] + n1[1] * v[61 * m] + n1[2] * v[62 * m] +
                    n1[3] * v[63 * m]);

  u[0] = n3[0] * (n2[0] * t1 + n2[1] * t2 + n2[2] * t3 + n2[3] * t4) +
         n3[1] * (n2[0] * t5 + n2[1] * t6 + n2[2] * t7 + n2[3] * t8) +
         n3[2] * (n2[0] * t9 + n2[1] * t10 + n2[2] * t11 + n2[3] * t12) +
         n3[3] * (n2[0] * t13 + n2[1] * t14 + n2[2] * t15 + n2[3] * t16);
}

inline void TacsAddTransTensor3DInterp4(const int m, const double *n1,
                                        const double *n2, const double *n3,
                                        const TacsScalar *u, TacsScalar *v) {
  TacsScalar a1 = n2[0] * u[0];
  TacsScalar a2 = n2[1] * u[0];
  TacsScalar a3 = n2[2] * u[0];
  TacsScalar a4 = n2[3] * u[0];

  TacsScalar t1 = n3[0] * a1;
  TacsScalar t2 = n3[0] * a2;
  TacsScalar t3 = n3[0] * a3;
  TacsScalar t4 = n3[0] * a4;
  TacsScalar t5 = n3[1] * a1;
  TacsScalar t6 = n3[1] * a2;
  TacsScalar t7 = n3[1] * a3;
  TacsScalar t8 = n3[1] * a4;
  TacsScalar t9 = n3[2] * a1;
  TacsScalar t10 = n3[2] * a2;
  TacsScalar t11 = n3[2] * a3;
  TacsScalar t12 = n3[2] * a4;
  TacsScalar t13 = n3[3] * a1;
  TacsScalar t14 = n3[3] * a2;
  TacsScalar t15 = n3[3] * a3;
  TacsScalar t16 = n3[3] * a4;

  v[0] += n1[0] * t1;
  v[m] += n1[1] * t1;
  v[2 * m] += n1[2] * t1;
  v[3 * m] += n1[3] * t1;
  v[4 * m] += n1[0] * t2;
  v[5 * m] += n1[1] * t2;
  v[6 * m] += n1[2] * t2;
  v[7 * m] += n1[3] * t2;
  v[8 * m] += n1[0] * t3;
  v[9 * m] += n1[1] * t3;
  v[10 * m] += n1[2] * t3;
  v[11 * m] += n1[3] * t3;
  v[12 * m] += n1[0] * t4;
  v[13 * m] += n1[1] * t4;
  v[14 * m] += n1[2] * t4;
  v[15 * m] += n1[3] * t4;
  v[16 * m] += n1[0] * t5;
  v[17 * m] += n1[1] * t5;
  v[18 * m] += n1[2] * t5;
  v[19 * m] += n1[3] * t5;
  v[20 * m] += n1[0] * t6;
  v[21 * m] += n1[1] * t6;
  v[22 * m] += n1[2] * t6;
  v[23 * m] += n1[3] * t6;
  v[24 * m] += n1[0] * t7;
  v[25 * m] += n1[1] * t7;
  v[26 * m] += n1[2] * t7;
  v[27 * m] += n1[3] * t7;
  v[28 * m] += n1[0] * t8;
  v[29 * m] += n1[1] * t8;
  v[30 * m] += n1[2] * t8;
  v[31 * m] += n1[3] * t8;
  v[32 * m] += n1[0] * t9;
  v[33 * m] += n1[1] * t9;
  v[34 * m] += n1[2] * t9;
  v[35 * m] += n1[3] * t9;
  v[36 * m] += n1[0] * t10;
  v[37 * m] += n1[1] * t10;
  v[38 * m] += n1[2] * t10;
  v[39 * m] += n1[3] * t10;
  v[40 * m] += n1[0] * t11;
  v[41 * m] += n1[1] * t11;
  v[42 * m] += n1[2] * t11;
  v[43 * m] += n1[3] * t11;
  v[44 * m] += n1[0] * t12;
  v[45 * m] += n1[1] * t12;
  v[46 * m] += n1[2] * t12;
  v[47 * m] += n1[3] * t12;
  v[48 * m] += n1[0] * t13;
  v[49 * m] += n1[1] * t13;
  v[50 * m] += n1[2] * t13;
  v[51 * m] += n1[3] * t13;
  v[52 * m] += n1[0] * t14;
  v[53 * m] += n1[1] * t14;
  v[54 * m] += n1[2] * t14;
  v[55 * m] += n1[3] * t14;
  v[56 * m] += n1[0] * t15;
  v[57 * m] += n1[1] * t15;
  v[58 * m] += n1[2] * t15;
  v[59 * m] += n1[3] * t15;
  v[60 * m] += n1[0] * t16;
  v[61 * m] += n1[1] * t16;
  v[62 * m] += n1[2] * t16;
  v[63 * m] += n1[3] * t16;
}

inline void TacsGradTensor3DInterp4(const int m, const double *n1,
                                    const double *n2, const double *n3,
                                    const double *n1x, const double *n2x,
                                    const double *n3x, const TacsScalar *v,
                                    TacsScalar *g) {
  TacsScalar t1 =
      (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] + n1[3] * v[3 * m]);
  TacsScalar t2 = (n1[0] * v[4 * m] + n1[1] * v[5 * m] + n1[2] * v[6 * m] +
                   n1[3] * v[7 * m]);
  TacsScalar t3 = (n1[0] * v[8 * m] + n1[1] * v[9 * m] + n1[2] * v[10 * m] +
                   n1[3] * v[11 * m]);
  TacsScalar t4 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m] +
                   n1[3] * v[15 * m]);
  TacsScalar t5 = (n1[0] * v[16 * m] + n1[1] * v[17 * m] + n1[2] * v[18 * m] +
                   n1[3] * v[19 * m]);
  TacsScalar t6 = (n1[0] * v[20 * m] + n1[1] * v[21 * m] + n1[2] * v[22 * m] +
                   n1[3] * v[23 * m]);
  TacsScalar t7 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m] +
                   n1[3] * v[27 * m]);
  TacsScalar t8 = (n1[0] * v[28 * m] + n1[1] * v[29 * m] + n1[2] * v[30 * m] +
                   n1[3] * v[31 * m]);
  TacsScalar t9 = (n1[0] * v[32 * m] + n1[1] * v[33 * m] + n1[2] * v[34 * m] +
                   n1[3] * v[35 * m]);
  TacsScalar t10 = (n1[0] * v[36 * m] + n1[1] * v[37 * m] + n1[2] * v[38 * m] +
                    n1[3] * v[39 * m]);
  TacsScalar t11 = (n1[0] * v[40 * m] + n1[1] * v[41 * m] + n1[2] * v[42 * m] +
                    n1[3] * v[43 * m]);
  TacsScalar t12 = (n1[0] * v[44 * m] + n1[1] * v[45 * m] + n1[2] * v[46 * m] +
                    n1[3] * v[47 * m]);
  TacsScalar t13 = (n1[0] * v[48 * m] + n1[1] * v[49 * m] + n1[2] * v[50 * m] +
                    n1[3] * v[51 * m]);
  TacsScalar t14 = (n1[0] * v[52 * m] + n1[1] * v[53 * m] + n1[2] * v[54 * m] +
                    n1[3] * v[55 * m]);
  TacsScalar t15 = (n1[0] * v[56 * m] + n1[1] * v[57 * m] + n1[2] * v[58 * m] +
                    n1[3] * v[59 * m]);
  TacsScalar t16 = (n1[0] * v[60 * m] + n1[1] * v[61 * m] + n1[2] * v[62 * m] +
                    n1[3] * v[63 * m]);

  TacsScalar t1x =
      (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] + n1x[3] * v[3 * m]);
  TacsScalar t2x = (n1x[0] * v[4 * m] + n1x[1] * v[5 * m] + n1x[2] * v[6 * m] +
                    n1x[3] * v[7 * m]);
  TacsScalar t3x = (n1x[0] * v[8 * m] + n1x[1] * v[9 * m] + n1x[2] * v[10 * m] +
                    n1x[3] * v[11 * m]);
  TacsScalar t4x = (n1x[0] * v[12 * m] + n1x[1] * v[13 * m] +
                    n1x[2] * v[14 * m] + n1x[3] * v[15 * m]);
  TacsScalar t5x = (n1x[0] * v[16 * m] + n1x[1] * v[17 * m] +
                    n1x[2] * v[18 * m] + n1x[3] * v[19 * m]);
  TacsScalar t6x = (n1x[0] * v[20 * m] + n1x[1] * v[21 * m] +
                    n1x[2] * v[22 * m] + n1x[3] * v[23 * m]);
  TacsScalar t7x = (n1x[0] * v[24 * m] + n1x[1] * v[25 * m] +
                    n1x[2] * v[26 * m] + n1x[3] * v[27 * m]);
  TacsScalar t8x = (n1x[0] * v[28 * m] + n1x[1] * v[29 * m] +
                    n1x[2] * v[30 * m] + n1x[3] * v[31 * m]);
  TacsScalar t9x = (n1x[0] * v[32 * m] + n1x[1] * v[33 * m] +
                    n1x[2] * v[34 * m] + n1x[3] * v[35 * m]);
  TacsScalar t10x = (n1x[0] * v[36 * m] + n1x[1] * v[37 * m] +
                     n1x[2] * v[38 * m] + n1x[3] * v[39 * m]);
  TacsScalar t11x = (n1x[0] * v[40 * m] + n1x[1] * v[41 * m] +
                     n1x[2] * v[42 * m] + n1x[3] * v[43 * m]);
  TacsScalar t12x = (n1x[0] * v[44 * m] + n1x[1] * v[45 * m] +
                     n1x[2] * v[46 * m] + n1x[3] * v[47 * m]);
  TacsScalar t13x = (n1x[0] * v[48 * m] + n1x[1] * v[49 * m] +
                     n1x[2] * v[50 * m] + n1x[3] * v[51 * m]);
  TacsScalar t14x = (n1x[0] * v[52 * m] + n1x[1] * v[53 * m] +
                     n1x[2] * v[54 * m] + n1x[3] * v[55 * m]);
  TacsScalar t15x = (n1x[0] * v[56 * m] + n1x[1] * v[57 * m] +
                     n1x[2] * v[58 * m] + n1x[3] * v[59 * m]);
  TacsScalar t16x = (n1x[0] * v[60 * m] + n1x[1] * v[61 * m] +
                     n1x[2] * v[62 * m] + n1x[3] * v[63 * m]);

  g[0] = n3[0] * (n2[0] * t1x + n2[1] * t2x + n2[2] * t3x + n2[3] * t4x) +
         n3[1] * (n2[0] * t5x + n2[1] * t6x + n2[2] * t7x + n2[3] * t8x) +
         n3[2] * (n2[0] * t9x + n2[1] * t10x + n2[2] * t11x + n2[3] * t12x) +
         n3[3] * (n2[0] * t13x + n2[1] * t14x + n2[2] * t15x + n2[3] * t16x);

  g[1] = n3[0] * (n2x[0] * t1 + n2x[1] * t2 + n2x[2] * t3 + n2x[3] * t4) +
         n3[1] * (n2x[0] * t5 + n2x[1] * t6 + n2x[2] * t7 + n2x[3] * t8) +
         n3[2] * (n2x[0] * t9 + n2x[1] * t10 + n2x[2] * t11 + n2x[3] * t12) +
         n3[3] * (n2x[0] * t13 + n2x[1] * t14 + n2x[2] * t15 + n2x[3] * t16);

  g[2] = n3x[0] * (n2[0] * t1 + n2[1] * t2 + n2[2] * t3 + n2[3] * t4) +
         n3x[1] * (n2[0] * t5 + n2[1] * t6 + n2[2] * t7 + n2[3] * t8) +
         n3x[2] * (n2[0] * t9 + n2[1] * t10 + n2[2] * t11 + n2[3] * t12) +
         n3x[3] * (n2[0] * t13 + n2[1] * t14 + n2[2] * t15 + n2[3] * t16);
}

inline void TacsAddGradTransTensor3DInterp4(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1];
  TacsScalar b2 = n2x[1] * g[1];
  TacsScalar b3 = n2x[2] * g[1];
  TacsScalar b4 = n2x[3] * g[1];

  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];
  TacsScalar c3 = n2[2] * g[2];
  TacsScalar c4 = n2[3] * g[2];

  TacsScalar t1x = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2x = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3x = n3[0] * b3 + n3x[0] * c3;
  TacsScalar t4x = n3[0] * b4 + n3x[0] * c4;
  TacsScalar t5x = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t6x = n3[1] * b2 + n3x[1] * c2;
  TacsScalar t7x = n3[1] * b3 + n3x[1] * c3;
  TacsScalar t8x = n3[1] * b4 + n3x[1] * c4;
  TacsScalar t9x = n3[2] * b1 + n3x[2] * c1;
  TacsScalar t10x = n3[2] * b2 + n3x[2] * c2;
  TacsScalar t11x = n3[2] * b3 + n3x[2] * c3;
  TacsScalar t12x = n3[2] * b4 + n3x[2] * c4;
  TacsScalar t13x = n3[3] * b1 + n3x[3] * c1;
  TacsScalar t14x = n3[3] * b2 + n3x[3] * c2;
  TacsScalar t15x = n3[3] * b3 + n3x[3] * c3;
  TacsScalar t16x = n3[3] * b4 + n3x[3] * c4;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];
  TacsScalar a3 = n2[2] * g[0];
  TacsScalar a4 = n2[3] * g[0];

  TacsScalar t1 = n3[0] * a1;
  TacsScalar t2 = n3[0] * a2;
  TacsScalar t3 = n3[0] * a3;
  TacsScalar t4 = n3[0] * a4;
  TacsScalar t5 = n3[1] * a1;
  TacsScalar t6 = n3[1] * a2;
  TacsScalar t7 = n3[1] * a3;
  TacsScalar t8 = n3[1] * a4;
  TacsScalar t9 = n3[2] * a1;
  TacsScalar t10 = n3[2] * a2;
  TacsScalar t11 = n3[2] * a3;
  TacsScalar t12 = n3[2] * a4;
  TacsScalar t13 = n3[3] * a1;
  TacsScalar t14 = n3[3] * a2;
  TacsScalar t15 = n3[3] * a3;
  TacsScalar t16 = n3[3] * a4;

  v[0] += n1[0] * t1x + n1x[0] * t1;
  v[m] += n1[1] * t1x + n1x[1] * t1;
  v[2 * m] += n1[2] * t1x + n1x[2] * t1;
  v[3 * m] += n1[3] * t1x + n1x[3] * t1;
  v[4 * m] += n1[0] * t2x + n1x[0] * t2;
  v[5 * m] += n1[1] * t2x + n1x[1] * t2;
  v[6 * m] += n1[2] * t2x + n1x[2] * t2;
  v[7 * m] += n1[3] * t2x + n1x[3] * t2;
  v[8 * m] += n1[0] * t3x + n1x[0] * t3;
  v[9 * m] += n1[1] * t3x + n1x[1] * t3;
  v[10 * m] += n1[2] * t3x + n1x[2] * t3;
  v[11 * m] += n1[3] * t3x + n1x[3] * t3;
  v[12 * m] += n1[0] * t4x + n1x[0] * t4;
  v[13 * m] += n1[1] * t4x + n1x[1] * t4;
  v[14 * m] += n1[2] * t4x + n1x[2] * t4;
  v[15 * m] += n1[3] * t4x + n1x[3] * t4;
  v[16 * m] += n1[0] * t5x + n1x[0] * t5;
  v[17 * m] += n1[1] * t5x + n1x[1] * t5;
  v[18 * m] += n1[2] * t5x + n1x[2] * t5;
  v[19 * m] += n1[3] * t5x + n1x[3] * t5;
  v[20 * m] += n1[0] * t6x + n1x[0] * t6;
  v[21 * m] += n1[1] * t6x + n1x[1] * t6;
  v[22 * m] += n1[2] * t6x + n1x[2] * t6;
  v[23 * m] += n1[3] * t6x + n1x[3] * t6;
  v[24 * m] += n1[0] * t7x + n1x[0] * t7;
  v[25 * m] += n1[1] * t7x + n1x[1] * t7;
  v[26 * m] += n1[2] * t7x + n1x[2] * t7;
  v[27 * m] += n1[3] * t7x + n1x[3] * t7;
  v[28 * m] += n1[0] * t8x + n1x[0] * t8;
  v[29 * m] += n1[1] * t8x + n1x[1] * t8;
  v[30 * m] += n1[2] * t8x + n1x[2] * t8;
  v[31 * m] += n1[3] * t8x + n1x[3] * t8;
  v[32 * m] += n1[0] * t9x + n1x[0] * t9;
  v[33 * m] += n1[1] * t9x + n1x[1] * t9;
  v[34 * m] += n1[2] * t9x + n1x[2] * t9;
  v[35 * m] += n1[3] * t9x + n1x[3] * t9;
  v[36 * m] += n1[0] * t10x + n1x[0] * t10;
  v[37 * m] += n1[1] * t10x + n1x[1] * t10;
  v[38 * m] += n1[2] * t10x + n1x[2] * t10;
  v[39 * m] += n1[3] * t10x + n1x[3] * t10;
  v[40 * m] += n1[0] * t11x + n1x[0] * t11;
  v[41 * m] += n1[1] * t11x + n1x[1] * t11;
  v[42 * m] += n1[2] * t11x + n1x[2] * t11;
  v[43 * m] += n1[3] * t11x + n1x[3] * t11;
  v[44 * m] += n1[0] * t12x + n1x[0] * t12;
  v[45 * m] += n1[1] * t12x + n1x[1] * t12;
  v[46 * m] += n1[2] * t12x + n1x[2] * t12;
  v[47 * m] += n1[3] * t12x + n1x[3] * t12;
  v[48 * m] += n1[0] * t13x + n1x[0] * t13;
  v[49 * m] += n1[1] * t13x + n1x[1] * t13;
  v[50 * m] += n1[2] * t13x + n1x[2] * t13;
  v[51 * m] += n1[3] * t13x + n1x[3] * t13;
  v[52 * m] += n1[0] * t14x + n1x[0] * t14;
  v[53 * m] += n1[1] * t14x + n1x[1] * t14;
  v[54 * m] += n1[2] * t14x + n1x[2] * t14;
  v[55 * m] += n1[3] * t14x + n1x[3] * t14;
  v[56 * m] += n1[0] * t15x + n1x[0] * t15;
  v[57 * m] += n1[1] * t15x + n1x[1] * t15;
  v[58 * m] += n1[2] * t15x + n1x[2] * t15;
  v[59 * m] += n1[3] * t15x + n1x[3] * t15;
  v[60 * m] += n1[0] * t16x + n1x[0] * t16;
  v[61 * m] += n1[1] * t16x + n1x[1] * t16;
  v[62 * m] += n1[2] * t16x + n1x[2] * t16;
  v[63 * m] += n1[3] * t16x + n1x[3] * t16;
}

inline void TacsInterpGradTensor3DInterp4(const int m, const double *n1,
                                          const double *n2, const double *n3,
                                          const double *n1x, const double *n2x,
                                          const double *n3x,
                                          const TacsScalar *v, TacsScalar *u,
                                          TacsScalar *g) {
  TacsScalar t1 =
      (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] + n1[3] * v[3 * m]);
  TacsScalar t2 = (n1[0] * v[4 * m] + n1[1] * v[5 * m] + n1[2] * v[6 * m] +
                   n1[3] * v[7 * m]);
  TacsScalar t3 = (n1[0] * v[8 * m] + n1[1] * v[9 * m] + n1[2] * v[10 * m] +
                   n1[3] * v[11 * m]);
  TacsScalar t4 = (n1[0] * v[12 * m] + n1[1] * v[13 * m] + n1[2] * v[14 * m] +
                   n1[3] * v[15 * m]);
  TacsScalar t5 = (n1[0] * v[16 * m] + n1[1] * v[17 * m] + n1[2] * v[18 * m] +
                   n1[3] * v[19 * m]);
  TacsScalar t6 = (n1[0] * v[20 * m] + n1[1] * v[21 * m] + n1[2] * v[22 * m] +
                   n1[3] * v[23 * m]);
  TacsScalar t7 = (n1[0] * v[24 * m] + n1[1] * v[25 * m] + n1[2] * v[26 * m] +
                   n1[3] * v[27 * m]);
  TacsScalar t8 = (n1[0] * v[28 * m] + n1[1] * v[29 * m] + n1[2] * v[30 * m] +
                   n1[3] * v[31 * m]);
  TacsScalar t9 = (n1[0] * v[32 * m] + n1[1] * v[33 * m] + n1[2] * v[34 * m] +
                   n1[3] * v[35 * m]);
  TacsScalar t10 = (n1[0] * v[36 * m] + n1[1] * v[37 * m] + n1[2] * v[38 * m] +
                    n1[3] * v[39 * m]);
  TacsScalar t11 = (n1[0] * v[40 * m] + n1[1] * v[41 * m] + n1[2] * v[42 * m] +
                    n1[3] * v[43 * m]);
  TacsScalar t12 = (n1[0] * v[44 * m] + n1[1] * v[45 * m] + n1[2] * v[46 * m] +
                    n1[3] * v[47 * m]);
  TacsScalar t13 = (n1[0] * v[48 * m] + n1[1] * v[49 * m] + n1[2] * v[50 * m] +
                    n1[3] * v[51 * m]);
  TacsScalar t14 = (n1[0] * v[52 * m] + n1[1] * v[53 * m] + n1[2] * v[54 * m] +
                    n1[3] * v[55 * m]);
  TacsScalar t15 = (n1[0] * v[56 * m] + n1[1] * v[57 * m] + n1[2] * v[58 * m] +
                    n1[3] * v[59 * m]);
  TacsScalar t16 = (n1[0] * v[60 * m] + n1[1] * v[61 * m] + n1[2] * v[62 * m] +
                    n1[3] * v[63 * m]);

  TacsScalar t1x =
      (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] + n1x[3] * v[3 * m]);
  TacsScalar t2x = (n1x[0] * v[4 * m] + n1x[1] * v[5 * m] + n1x[2] * v[6 * m] +
                    n1x[3] * v[7 * m]);
  TacsScalar t3x = (n1x[0] * v[8 * m] + n1x[1] * v[9 * m] + n1x[2] * v[10 * m] +
                    n1x[3] * v[11 * m]);
  TacsScalar t4x = (n1x[0] * v[12 * m] + n1x[1] * v[13 * m] +
                    n1x[2] * v[14 * m] + n1x[3] * v[15 * m]);
  TacsScalar t5x = (n1x[0] * v[16 * m] + n1x[1] * v[17 * m] +
                    n1x[2] * v[18 * m] + n1x[3] * v[19 * m]);
  TacsScalar t6x = (n1x[0] * v[20 * m] + n1x[1] * v[21 * m] +
                    n1x[2] * v[22 * m] + n1x[3] * v[23 * m]);
  TacsScalar t7x = (n1x[0] * v[24 * m] + n1x[1] * v[25 * m] +
                    n1x[2] * v[26 * m] + n1x[3] * v[27 * m]);
  TacsScalar t8x = (n1x[0] * v[28 * m] + n1x[1] * v[29 * m] +
                    n1x[2] * v[30 * m] + n1x[3] * v[31 * m]);
  TacsScalar t9x = (n1x[0] * v[32 * m] + n1x[1] * v[33 * m] +
                    n1x[2] * v[34 * m] + n1x[3] * v[35 * m]);
  TacsScalar t10x = (n1x[0] * v[36 * m] + n1x[1] * v[37 * m] +
                     n1x[2] * v[38 * m] + n1x[3] * v[39 * m]);
  TacsScalar t11x = (n1x[0] * v[40 * m] + n1x[1] * v[41 * m] +
                     n1x[2] * v[42 * m] + n1x[3] * v[43 * m]);
  TacsScalar t12x = (n1x[0] * v[44 * m] + n1x[1] * v[45 * m] +
                     n1x[2] * v[46 * m] + n1x[3] * v[47 * m]);
  TacsScalar t13x = (n1x[0] * v[48 * m] + n1x[1] * v[49 * m] +
                     n1x[2] * v[50 * m] + n1x[3] * v[51 * m]);
  TacsScalar t14x = (n1x[0] * v[52 * m] + n1x[1] * v[53 * m] +
                     n1x[2] * v[54 * m] + n1x[3] * v[55 * m]);
  TacsScalar t15x = (n1x[0] * v[56 * m] + n1x[1] * v[57 * m] +
                     n1x[2] * v[58 * m] + n1x[3] * v[59 * m]);
  TacsScalar t16x = (n1x[0] * v[60 * m] + n1x[1] * v[61 * m] +
                     n1x[2] * v[62 * m] + n1x[3] * v[63 * m]);

  TacsScalar a1 = (n2[0] * t1 + n2[1] * t2 + n2[2] * t3 + n2[3] * t4);
  TacsScalar a2 = (n2[0] * t5 + n2[1] * t6 + n2[2] * t7 + n2[3] * t8);
  TacsScalar a3 = (n2[0] * t9 + n2[1] * t10 + n2[2] * t11 + n2[3] * t12);
  TacsScalar a4 = (n2[0] * t13 + n2[1] * t14 + n2[2] * t15 + n2[3] * t16);

  u[0] = n3[0] * a1 + n3[1] * a2 + n3[2] * a3 + n3[3] * a4;

  g[0] = n3[0] * (n2[0] * t1x + n2[1] * t2x + n2[2] * t3x + n2[3] * t4x) +
         n3[1] * (n2[0] * t5x + n2[1] * t6x + n2[2] * t7x + n2[3] * t8x) +
         n3[2] * (n2[0] * t9x + n2[1] * t10x + n2[2] * t11x + n2[3] * t12x) +
         n3[3] * (n2[0] * t13x + n2[1] * t14x + n2[2] * t15x + n2[3] * t16x);

  g[1] = n3[0] * (n2x[0] * t1 + n2x[1] * t2 + n2x[2] * t3 + n2x[3] * t4) +
         n3[1] * (n2x[0] * t5 + n2x[1] * t6 + n2x[2] * t7 + n2x[3] * t8) +
         n3[2] * (n2x[0] * t9 + n2x[1] * t10 + n2x[2] * t11 + n2x[3] * t12) +
         n3[3] * (n2x[0] * t13 + n2x[1] * t14 + n2x[2] * t15 + n2x[3] * t16);

  g[2] = n3x[0] * a1 + n3x[1] * a2 + n3x[2] * a3 + n3x[3] * a4;
}

inline void TacsAddInterpGradTransTensor3DInterp4(
    const int m, const double *n1, const double *n2, const double *n3,
    const double *n1x, const double *n2x, const double *n3x,
    const TacsScalar *u, const TacsScalar *g, TacsScalar *v) {
  TacsScalar b1 = n2x[0] * g[1] + n2[0] * u[0];
  TacsScalar b2 = n2x[1] * g[1] + n2[1] * u[0];
  TacsScalar b3 = n2x[2] * g[1] + n2[2] * u[0];
  TacsScalar b4 = n2x[3] * g[1] + n2[3] * u[0];

  TacsScalar c1 = n2[0] * g[2];
  TacsScalar c2 = n2[1] * g[2];
  TacsScalar c3 = n2[2] * g[2];
  TacsScalar c4 = n2[3] * g[2];

  TacsScalar t1x = n3[0] * b1 + n3x[0] * c1;
  TacsScalar t2x = n3[0] * b2 + n3x[0] * c2;
  TacsScalar t3x = n3[0] * b3 + n3x[0] * c3;
  TacsScalar t4x = n3[0] * b4 + n3x[0] * c4;
  TacsScalar t5x = n3[1] * b1 + n3x[1] * c1;
  TacsScalar t6x = n3[1] * b2 + n3x[1] * c2;
  TacsScalar t7x = n3[1] * b3 + n3x[1] * c3;
  TacsScalar t8x = n3[1] * b4 + n3x[1] * c4;
  TacsScalar t9x = n3[2] * b1 + n3x[2] * c1;
  TacsScalar t10x = n3[2] * b2 + n3x[2] * c2;
  TacsScalar t11x = n3[2] * b3 + n3x[2] * c3;
  TacsScalar t12x = n3[2] * b4 + n3x[2] * c4;
  TacsScalar t13x = n3[3] * b1 + n3x[3] * c1;
  TacsScalar t14x = n3[3] * b2 + n3x[3] * c2;
  TacsScalar t15x = n3[3] * b3 + n3x[3] * c3;
  TacsScalar t16x = n3[3] * b4 + n3x[3] * c4;

  TacsScalar a1 = n2[0] * g[0];
  TacsScalar a2 = n2[1] * g[0];
  TacsScalar a3 = n2[2] * g[0];
  TacsScalar a4 = n2[3] * g[0];

  TacsScalar t1 = n3[0] * a1;
  TacsScalar t2 = n3[0] * a2;
  TacsScalar t3 = n3[0] * a3;
  TacsScalar t4 = n3[0] * a4;
  TacsScalar t5 = n3[1] * a1;
  TacsScalar t6 = n3[1] * a2;
  TacsScalar t7 = n3[1] * a3;
  TacsScalar t8 = n3[1] * a4;
  TacsScalar t9 = n3[2] * a1;
  TacsScalar t10 = n3[2] * a2;
  TacsScalar t11 = n3[2] * a3;
  TacsScalar t12 = n3[2] * a4;
  TacsScalar t13 = n3[3] * a1;
  TacsScalar t14 = n3[3] * a2;
  TacsScalar t15 = n3[3] * a3;
  TacsScalar t16 = n3[3] * a4;

  v[0] += n1[0] * t1x + n1x[0] * t1;
  v[m] += n1[1] * t1x + n1x[1] * t1;
  v[2 * m] += n1[2] * t1x + n1x[2] * t1;
  v[3 * m] += n1[3] * t1x + n1x[3] * t1;
  v[4 * m] += n1[0] * t2x + n1x[0] * t2;
  v[5 * m] += n1[1] * t2x + n1x[1] * t2;
  v[6 * m] += n1[2] * t2x + n1x[2] * t2;
  v[7 * m] += n1[3] * t2x + n1x[3] * t2;
  v[8 * m] += n1[0] * t3x + n1x[0] * t3;
  v[9 * m] += n1[1] * t3x + n1x[1] * t3;
  v[10 * m] += n1[2] * t3x + n1x[2] * t3;
  v[11 * m] += n1[3] * t3x + n1x[3] * t3;
  v[12 * m] += n1[0] * t4x + n1x[0] * t4;
  v[13 * m] += n1[1] * t4x + n1x[1] * t4;
  v[14 * m] += n1[2] * t4x + n1x[2] * t4;
  v[15 * m] += n1[3] * t4x + n1x[3] * t4;
  v[16 * m] += n1[0] * t5x + n1x[0] * t5;
  v[17 * m] += n1[1] * t5x + n1x[1] * t5;
  v[18 * m] += n1[2] * t5x + n1x[2] * t5;
  v[19 * m] += n1[3] * t5x + n1x[3] * t5;
  v[20 * m] += n1[0] * t6x + n1x[0] * t6;
  v[21 * m] += n1[1] * t6x + n1x[1] * t6;
  v[22 * m] += n1[2] * t6x + n1x[2] * t6;
  v[23 * m] += n1[3] * t6x + n1x[3] * t6;
  v[24 * m] += n1[0] * t7x + n1x[0] * t7;
  v[25 * m] += n1[1] * t7x + n1x[1] * t7;
  v[26 * m] += n1[2] * t7x + n1x[2] * t7;
  v[27 * m] += n1[3] * t7x + n1x[3] * t7;
  v[28 * m] += n1[0] * t8x + n1x[0] * t8;
  v[29 * m] += n1[1] * t8x + n1x[1] * t8;
  v[30 * m] += n1[2] * t8x + n1x[2] * t8;
  v[31 * m] += n1[3] * t8x + n1x[3] * t8;
  v[32 * m] += n1[0] * t9x + n1x[0] * t9;
  v[33 * m] += n1[1] * t9x + n1x[1] * t9;
  v[34 * m] += n1[2] * t9x + n1x[2] * t9;
  v[35 * m] += n1[3] * t9x + n1x[3] * t9;
  v[36 * m] += n1[0] * t10x + n1x[0] * t10;
  v[37 * m] += n1[1] * t10x + n1x[1] * t10;
  v[38 * m] += n1[2] * t10x + n1x[2] * t10;
  v[39 * m] += n1[3] * t10x + n1x[3] * t10;
  v[40 * m] += n1[0] * t11x + n1x[0] * t11;
  v[41 * m] += n1[1] * t11x + n1x[1] * t11;
  v[42 * m] += n1[2] * t11x + n1x[2] * t11;
  v[43 * m] += n1[3] * t11x + n1x[3] * t11;
  v[44 * m] += n1[0] * t12x + n1x[0] * t12;
  v[45 * m] += n1[1] * t12x + n1x[1] * t12;
  v[46 * m] += n1[2] * t12x + n1x[2] * t12;
  v[47 * m] += n1[3] * t12x + n1x[3] * t12;
  v[48 * m] += n1[0] * t13x + n1x[0] * t13;
  v[49 * m] += n1[1] * t13x + n1x[1] * t13;
  v[50 * m] += n1[2] * t13x + n1x[2] * t13;
  v[51 * m] += n1[3] * t13x + n1x[3] * t13;
  v[52 * m] += n1[0] * t14x + n1x[0] * t14;
  v[53 * m] += n1[1] * t14x + n1x[1] * t14;
  v[54 * m] += n1[2] * t14x + n1x[2] * t14;
  v[55 * m] += n1[3] * t14x + n1x[3] * t14;
  v[56 * m] += n1[0] * t15x + n1x[0] * t15;
  v[57 * m] += n1[1] * t15x + n1x[1] * t15;
  v[58 * m] += n1[2] * t15x + n1x[2] * t15;
  v[59 * m] += n1[3] * t15x + n1x[3] * t15;
  v[60 * m] += n1[0] * t16x + n1x[0] * t16;
  v[61 * m] += n1[1] * t16x + n1x[1] * t16;
  v[62 * m] += n1[2] * t16x + n1x[2] * t16;
  v[63 * m] += n1[3] * t16x + n1x[3] * t16;
}

inline void TACSInterpAllTensor3DInterp4(const int m, const int i,
                                         const double N[], const double Nx[],
                                         const TacsScalar v[],
                                         TacsScalar out[]) {
  const int j = 3 * i + m;
  const double *Nz = N;
  const double *Nzx = Nx;
  for (int n = 0; n < 4; n++) {
    TacsInterpGradTensor3DInterp4(m, &N[0], &N[0], Nz, &Nx[0], &Nx[0], Nzx, v,
                                  &out[i], &out[j]);
    TacsInterpGradTensor3DInterp4(m, &N[4], &N[0], Nz, &Nx[4], &Nx[0], Nzx, v,
                                  &out[4 * m + i], &out[4 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[8], &N[0], Nz, &Nx[8], &Nx[0], Nzx, v,
                                  &out[8 * m + i], &out[8 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[12], &N[0], Nz, &Nx[12], &Nx[0], Nzx, v,
                                  &out[12 * m + i], &out[12 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[0], &N[4], Nz, &Nx[0], &Nx[4], Nzx, v,
                                  &out[16 * m + i], &out[16 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[4], &N[4], Nz, &Nx[4], &Nx[4], Nzx, v,
                                  &out[20 * m + i], &out[20 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[8], &N[4], Nz, &Nx[8], &Nx[4], Nzx, v,
                                  &out[24 * m + i], &out[24 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[12], &N[4], Nz, &Nx[12], &Nx[4], Nzx, v,
                                  &out[28 * m + i], &out[28 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[0], &N[8], Nz, &Nx[0], &Nx[8], Nzx, v,
                                  &out[32 * m + i], &out[32 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[4], &N[8], Nz, &Nx[4], &Nx[8], Nzx, v,
                                  &out[36 * m + i], &out[36 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[8], &N[8], Nz, &Nx[8], &Nx[8], Nzx, v,
                                  &out[40 * m + i], &out[40 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[12], &N[8], Nz, &Nx[12], &Nx[8], Nzx, v,
                                  &out[44 * m + i], &out[44 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[0], &N[12], Nz, &Nx[0], &Nx[12], Nzx, v,
                                  &out[48 * m + i], &out[48 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[4], &N[12], Nz, &Nx[4], &Nx[12], Nzx, v,
                                  &out[52 * m + i], &out[52 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[8], &N[12], Nz, &Nx[8], &Nx[12], Nzx, v,
                                  &out[56 * m + i], &out[56 * m + j]);
    TacsInterpGradTensor3DInterp4(m, &N[12], &N[12], Nz, &Nx[12], &Nx[12], Nzx,
                                  v, &out[60 * m + i], &out[60 * m + j]);

    out += 64 * m;
    Nz += 4;
    Nzx += 4;
  }
}

inline void TacsAddAllTransTensor3DInterp4(const int m, const int i,
                                           const double N[], const double Nx[],
                                           const TacsScalar in[],
                                           TacsScalar v[]) {
  const int j = 3 * i + m;
  const double *Nz = N;
  const double *Nzx = Nx;
  for (int n = 0; n < 4; n++) {
    TacsAddInterpGradTransTensor3DInterp4(m, &N[0], &N[0], Nz, &Nx[0], &Nx[0],
                                          Nzx, &in[i], &in[j], v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[4], &N[0], Nz, &Nx[4], &Nx[0],
                                          Nzx, &in[4 * m + i], &in[4 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[8], &N[0], Nz, &Nx[8], &Nx[0],
                                          Nzx, &in[8 * m + i], &in[8 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[12], &N[0], Nz, &Nx[12], &Nx[0],
                                          Nzx, &in[12 * m + i], &in[12 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[0], &N[4], Nz, &Nx[0], &Nx[4],
                                          Nzx, &in[16 * m + i], &in[16 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[4], &N[4], Nz, &Nx[4], &Nx[4],
                                          Nzx, &in[20 * m + i], &in[20 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[8], &N[4], Nz, &Nx[8], &Nx[4],
                                          Nzx, &in[24 * m + i], &in[24 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[12], &N[4], Nz, &Nx[12], &Nx[4],
                                          Nzx, &in[28 * m + i], &in[28 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[0], &N[8], Nz, &Nx[0], &Nx[8],
                                          Nzx, &in[32 * m + i], &in[32 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[4], &N[8], Nz, &Nx[4], &Nx[8],
                                          Nzx, &in[36 * m + i], &in[36 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[8], &N[8], Nz, &Nx[8], &Nx[8],
                                          Nzx, &in[40 * m + i], &in[40 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[12], &N[8], Nz, &Nx[12], &Nx[8],
                                          Nzx, &in[44 * m + i], &in[44 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[0], &N[12], Nz, &Nx[0], &Nx[12],
                                          Nzx, &in[48 * m + i], &in[48 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[4], &N[12], Nz, &Nx[4], &Nx[12],
                                          Nzx, &in[52 * m + i], &in[52 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[8], &N[12], Nz, &Nx[8], &Nx[12],
                                          Nzx, &in[56 * m + i], &in[56 * m + j],
                                          v);
    TacsAddInterpGradTransTensor3DInterp4(m, &N[12], &N[12], Nz, &Nx[12],
                                          &Nx[12], Nzx, &in[60 * m + i],
                                          &in[60 * m + j], v);

    in += 64 * m;
    Nz += 4;
    Nzx += 4;
  }
}

/*
  3D Tensor product functions for p = 4
*/
void TACSInterpAllTensor3DInterp5(const int m, const double N[],
                                  const double Nx[], const TacsScalar values[],
                                  TacsScalar out[]);
void TACSInterpAllTensor3DInterp5VarsPerNode1(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);
void TACSInterpAllTensor3DInterp5VarsPerNode3(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);
void TACSInterpAllTensor3DInterp5VarsPerNode4(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);

void TacsAddAllTransTensor3DInterp5(const int m, const double N[],
                                    const double Nx[], const TacsScalar in[],
                                    TacsScalar values[]);
void TacsAddAllTransTensor3DInterp5VarsPerNode1(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);
void TacsAddAllTransTensor3DInterp5VarsPerNode3(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);
void TacsAddAllTransTensor3DInterp5VarsPerNode4(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);

/*
  3D Tensor product functions for p = 5
*/
void TACSInterpAllTensor3DInterp6(const int m, const double N[],
                                  const double Nx[], const TacsScalar values[],
                                  TacsScalar out[]);
void TACSInterpAllTensor3DInterp6VarsPerNode1(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);
void TACSInterpAllTensor3DInterp6VarsPerNode3(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);
void TACSInterpAllTensor3DInterp6VarsPerNode4(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]);

void TacsAddAllTransTensor3DInterp6(const int m, const double N[],
                                    const double Nx[], const TacsScalar in[],
                                    TacsScalar values[]);
void TacsAddAllTransTensor3DInterp6VarsPerNode1(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);
void TacsAddAllTransTensor3DInterp6VarsPerNode3(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);
void TacsAddAllTransTensor3DInterp6VarsPerNode4(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]);

#endif  // TACS_TENSOR_PRODUCT_BASIS_IMPL_H
