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

#include "LargeRotUtils.h"

#include "FElibrary.h"

/*
  Common structural shell utilities for small strains
*/

TACS_BEGIN_NAMESPACE(largerot)

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

static inline void transform_vector3d(TacsScalar out[], const TacsScalar in[],
                                      const TacsScalar t[]) {
  out[0] = t[0] * in[0] + t[1] * in[1] + t[2] * in[2];
  out[1] = t[3] * in[0] + t[4] * in[1] + t[5] * in[2];
  out[2] = t[6] * in[0] + t[7] * in[1] + t[8] * in[2];
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
  Transform the displacement gradient from one reference frame to
  another - note that this only transforms the components and not the
  derivatives.
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
}

/*
  Evaluate the derivative of the displacement gradient w.r.t. the nodal
  displacements.
*/
static inline void transform_displ_gradient_bmat(TacsScalar Ued[], int ii,
                                                 const TacsScalar t[],
                                                 double Na, double Nb) {
  if (ii < 3) {
    // Transform the displacement gradient
    Ued[0] = Na * t[ii];
    Ued[2] = Na * t[3 + ii];
    Ued[4] = Na * t[6 + ii];

    Ued[1] = Nb * t[ii];
    Ued[3] = Nb * t[3 + ii];
    Ued[5] = Nb * t[6 + ii];
  } else {
    Ued[0] = 0.0;
    Ued[2] = 0.0;
    Ued[4] = 0.0;

    Ued[1] = 0.0;
    Ued[3] = 0.0;
    Ued[5] = 0.0;
  }
}

/*
  Evaluate the matrix C such that C*n (where n is the normal) is the
  rate of change of displacements through the thickness of the shell.

  The matrix C = Q^{T} - I, where Q is a rotation matrix.

  Also, evaluate the derivatives of the matrix C w.r.t. the parameters
  theta and beta. These are stored in the arrays Ct and Cb.

  The rotation matrix is as follows:

  Q = C_1(theta_1) C_2(theta_2) C_3(theta_3)

  Q^{T} = C_3^{T} C_2^{T} C_1^{T}

  =
  [ c3  -s3    ][ c2     s2 ][ 1          ]
  [ s3   c3    ][     1     ][    c1  -s1 ]
  [          1 ][-s2     c2 ][    s1   c1 ]

  =
  [ c2*c3 |  -s3 |  s2*c3 ][ 1          ]
  [ c2*s3 |   c3 |  s2*s3 ][    c1  -s1 ]
  [ -s2   |   0  |  c2    ][    s1   c1 ]

  =
  [ c2*c3  | -c1*s3 + s1*s2*c3  |  s1*s3 + c1*s2*c3 ]
  [ c2*s3  |  c1*c3 + s1*s2*s3  | -s1*c3 + c1*s2*s3 ]
  [ -s2    |  s1*c2             |  c1*c2            ]
*/
void compute_rate_matrix(TacsScalar C[], TacsScalar c1, TacsScalar s1,
                         TacsScalar c2, TacsScalar s2, TacsScalar c3,
                         TacsScalar s3) {
  // Evaluate C = Q^{T} - I
  C[0] = c2 * c3 - 1.0;
  C[1] = -c1 * s3 + s1 * s2 * c3;
  C[2] = s1 * s3 + c1 * s2 * c3;

  C[3] = c2 * s3;
  C[4] = c1 * c3 + s1 * s2 * s3 - 1.0;
  C[5] = -s1 * c3 + c1 * s2 * s3;

  C[6] = -s2;
  C[7] = s1 * c2;
  C[8] = c1 * c2 - 1.0;
}

/*
  Compute the rate matrix and the first derivatives of the rate matrix
  C. The derivatives are stored as follows:

  C,1 = &Ctt[0]
  C,2 = &Ctt[9]
  C,3 = &Ctt[18]
*/
void compute_rate_matrix(TacsScalar C[], TacsScalar Ct[], TacsScalar c1,
                         TacsScalar s1, TacsScalar c2, TacsScalar s2,
                         TacsScalar c3, TacsScalar s3) {
  // Evaluate C = Q^{T} - I
  C[0] = c2 * c3 - 1.0;
  C[1] = -c1 * s3 + s1 * s2 * c3;
  C[2] = s1 * s3 + c1 * s2 * c3;

  C[3] = c2 * s3;
  C[4] = c1 * c3 + s1 * s2 * s3 - 1.0;
  C[5] = -s1 * c3 + c1 * s2 * s3;

  C[6] = -s2;
  C[7] = s1 * c2;
  C[8] = c1 * c2 - 1.0;

  // Now evaluate the data in Ct
  // Evaluate the derivative of C w.r.t. theta_1
  Ct[0] = 0.0;
  Ct[1] = s1 * s3 + c1 * s2 * c3;
  Ct[2] = c1 * s3 - s1 * s2 * c3;

  Ct[3] = 0.0;
  Ct[4] = -s1 * c3 + c1 * s2 * s3;
  Ct[5] = -c1 * c3 - s1 * s2 * s3;  // - 1.0

  Ct[6] = 0.0;
  Ct[7] = c1 * c2;  // 1.0
  Ct[8] = -s1 * c2;
  Ct += 9;

  // Evaluate the derivative of C w.r.t. theta_2
  Ct[0] = -s2 * c3;
  Ct[1] = s1 * c2 * c3;
  Ct[2] = c1 * c2 * c3;  // 1.0

  Ct[3] = -s2 * s3;
  Ct[4] = s1 * c2 * s3;
  Ct[5] = c1 * c2 * s3;

  Ct[6] = -c2;  // -1.0
  Ct[7] = -s1 * s2;
  Ct[8] = -c1 * s2;
  Ct += 9;

  // Evaluate the derivative of C w.r.t. theta_3
  Ct[0] = -c2 * s3;
  Ct[1] = -c1 * c3 - s1 * s2 * s3;  // -1.0
  Ct[2] = s1 * c3 - c1 * s2 * s3;

  Ct[3] = c2 * c3;  // 1.0
  Ct[4] = -c1 * s3 + s1 * s2 * c3;
  Ct[5] = s1 * s3 + c1 * s2 * c3;

  Ct[6] = 0.0;
  Ct[7] = 0.0;
  Ct[8] = 0.0;
}

/*
  Compute the second derivatives of C.

  These derivatives are required when computing the derivative of the
  strain expression w.r.t. the nodal displacements.

  The second derivatives of the matrix C are stored as follows:

  C,11 = &Ctt[0]
  C,22 = &Ctt[9]
  C,33 = &Ctt[18]
  C,12 = &Ctt[27]
  C,13 = &Ctt[36]
  C,23 = &Ctt[45]
*/
void compute_2nd_rate_matrix(TacsScalar Ctt[], TacsScalar c1, TacsScalar s1,
                             TacsScalar c2, TacsScalar s2, TacsScalar c3,
                             TacsScalar s3) {
  // Now evaluate the data in Ct
  // Evaluate the derivative of C w.r.t. theta_1, theta_1
  Ctt[0] = 0.0;
  Ctt[1] = c1 * s3 - s1 * s2 * c3;
  Ctt[2] = -s1 * s3 - c1 * s2 * c3;

  Ctt[3] = 0.0;
  Ctt[4] = -c1 * c3 - s1 * s2 * s3;
  Ctt[5] = s1 * c3 - c1 * s2 * s3;

  Ctt[6] = 0.0;
  Ctt[7] = -s1 * c2;
  Ctt[8] = -c1 * c2;
  Ctt += 9;

  // Evaluate the derivative of C w.r.t. theta_2, theta_2
  Ctt[0] = -c2 * c3;
  Ctt[1] = -s1 * s2 * c3;
  Ctt[2] = -c1 * s2 * c3;

  Ctt[3] = -c2 * s3;
  Ctt[4] = -s1 * s2 * s3;
  Ctt[5] = -c1 * s2 * s3;

  Ctt[6] = s2;
  Ctt[7] = -s1 * c2;
  Ctt[8] = -c1 * c2;
  Ctt += 9;

  // Evaluate the derivative of C w.r.t. theta_3, theta_3
  Ctt[0] = -c2 * c3;
  Ctt[1] = c1 * s3 - s1 * s2 * c3;
  Ctt[2] = -s1 * s3 - c1 * s2 * c3;

  Ctt[3] = -c2 * s3;
  Ctt[4] = -c1 * c3 - s1 * s2 * s3;
  Ctt[5] = s1 * c3 - c1 * s2 * s3;

  Ctt[6] = 0.0;
  Ctt[7] = 0.0;
  Ctt[8] = 0.0;
  Ctt += 9;

  // Now evaluate the data in Ct
  // Evaluate the derivative of C w.r.t. theta_1, theta_2
  Ctt[0] = 0.0;
  Ctt[1] = c1 * c2 * c3;
  Ctt[2] = -s1 * c2 * c3;

  Ctt[3] = 0.0;
  Ctt[4] = c1 * c2 * s3;
  Ctt[5] = -s1 * c2 * s3;

  Ctt[6] = 0.0;
  Ctt[7] = -c1 * s2;
  Ctt[8] = s1 * s2;
  Ctt += 9;

  // Evaluate the derivative of C w.r.t. theta_1, theta_3
  Ctt[0] = 0.0;
  Ctt[1] = s1 * c3 - c1 * s2 * s3;
  Ctt[2] = c1 * c3 + s1 * s2 * s3;

  Ctt[3] = 0.0;
  Ctt[4] = s1 * s3 + c1 * s2 * c3;
  Ctt[5] = c1 * s3 - s1 * s2 * c3;

  Ctt[6] = 0.0;
  Ctt[7] = 0.0;
  Ctt[8] = 0.0;
  Ctt += 9;

  // Evaluate the derivative of C w.r.t. theta_2, theta_3
  Ctt[0] = s2 * s3;
  Ctt[1] = -s1 * c2 * s3;
  Ctt[2] = -c1 * c2 * s3;

  Ctt[3] = -s2 * c3;
  Ctt[4] = s1 * c2 * c3;
  Ctt[5] = c1 * c2 * c3;

  Ctt[6] = 0.0;
  Ctt[7] = 0.0;
  Ctt[8] = 0.0;
}

/*
  Compute the thrid derivative of the rate matrix, C, where C*n is the
  rate of change of the displacements through the thickness of the shell.

  Note that since the matrix C is expressed in terms of Euler angles,
  certain combinations of the the third derivatives are the negative
  of the first derivatives:

  C,iii = - C,i

  (this applies to single Euler-angle matrices as well).

  The order of the matrices supplied is:

  C,112 = &Cttt[0];
  C,113 = &Cttt[9];
  C,122 = &Cttt[18];
  C,123 = &Cttt[27];
  C,133 = &Cttt[36];
  C,223 = &Cttt[45];
  C,233 = &Cttt[54];
*/
void compute_3rd_rate_matrix(TacsScalar Cttt[], TacsScalar c1, TacsScalar s1,
                             TacsScalar c2, TacsScalar s2, TacsScalar c3,
                             TacsScalar s3) {
  // Compute C,112
  Cttt[0] = 0.0;
  Cttt[1] = -s1 * c2 * c3;
  Cttt[2] = -c1 * c2 * c3;

  Cttt[3] = 0.0;
  Cttt[4] = -s1 * c2 * s3;
  Cttt[5] = -c1 * c2 * s3;

  Cttt[6] = 0.0;
  Cttt[7] = s1 * s2;
  Cttt[8] = c1 * s2;
  Cttt += 9;

  // Compute C,113
  Cttt[0] = 0.0;
  Cttt[1] = c1 * c3 + s1 * s2 * s3;
  Cttt[2] = -s1 * c3 + c1 * s2 * s3;

  Cttt[3] = 0.0;
  Cttt[4] = c1 * s3 - s1 * s2 * c3;
  Cttt[5] = -s1 * s3 - c1 * s2 * c3;

  Cttt[6] = 0.0;
  Cttt[7] = 0.0;
  Cttt[8] = 0.0;
  Cttt += 9;

  // Compute C,122
  Cttt[0] = 0.0;
  Cttt[1] = -c1 * s2 * c3;
  Cttt[2] = s1 * s2 * c3;

  Cttt[3] = 0.0;
  Cttt[4] = -c1 * s2 * s3;
  Cttt[5] = s1 * s2 * s3;

  Cttt[6] = 0.0;
  Cttt[7] = -c1 * c2;
  Cttt[8] = s1 * c2;
  Cttt += 9;

  // Compute C,123
  Cttt[0] = 0.0;
  Cttt[1] = -c1 * c2 * s3;
  Cttt[2] = s1 * c2 * s3;

  Cttt[3] = 0.0;
  Cttt[4] = c1 * c2 * c3;
  Cttt[5] = -s1 * c2 * c3;

  Cttt[6] = 0.0;
  Cttt[7] = 0.0;
  Cttt[8] = 0.0;
  Cttt += 9;

  // Compute C,133
  Cttt[0] = 0.0;
  Cttt[1] = -s1 * s3 - c1 * s2 * c3;
  Cttt[2] = -c1 * s3 + s1 * s2 * c3;

  Cttt[3] = 0.0;
  Cttt[4] = s1 * c3 - c1 * s2 * s3;
  Cttt[5] = c1 * c3 + s1 * s2 * s3;

  Cttt[6] = 0.0;
  Cttt[7] = 0.0;
  Cttt[8] = 0.0;
  Cttt += 9;

  // Compute C,223
  Cttt[0] = c2 * s3;
  Cttt[1] = s1 * s2 * s3;
  Cttt[2] = c1 * s2 * s3;

  Cttt[3] = -c2 * c3;
  Cttt[4] = -s1 * s2 * c3;
  Cttt[5] = -c1 * s2 * c3;

  Cttt[6] = 0.0;
  Cttt[7] = 0.0;
  Cttt[8] = 0.0;
  Cttt += 9;

  // Compute C,233
  Cttt[0] = s2 * c3;
  Cttt[1] = -s1 * c2 * c3;
  Cttt[2] = -c1 * c2 * c3;

  Cttt[3] = s2 * s3;
  Cttt[4] = -s1 * c2 * s3;
  Cttt[5] = -c1 * c2 * s3;

  Cttt[6] = 0.0;
  Cttt[7] = 0.0;
  Cttt[8] = 0.0;
  Cttt += 9;
}

/*
  Test the implementation of the rate matrix

  input:
  U: an array of the rotation variables
  dh: the finite-difference step size
*/
void test_rate_matrix(TacsScalar U[], double dh) {
  TacsScalar Q[9], C[9], Ct[27], Ctt[54], Cttt[63];
  TacsScalar c1, s1, c2, s2, c3, s3;
  c1 = cos(U[0]);
  s1 = sin(U[0]);
  c2 = cos(U[1]);
  s2 = sin(U[1]);
  c3 = cos(U[2]);
  s3 = sin(U[2]);
  compute_rate_matrix(C, Ct, c1, s1, c2, s2, c3, s3);
  compute_2nd_rate_matrix(Ctt, c1, s1, c2, s2, c3, s3);
  compute_3rd_rate_matrix(Cttt, c1, s1, c2, s2, c3, s3);

  fprintf(stderr, "Orthogonality check\n");
  for (int i = 0; i < 9; i++) {
    Q[i] = C[i];
  }
  for (int i = 0; i < 3; i++) {
    Q[4 * i] += 1.0;
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(stderr, "Q[%d]^{T}Q[%d] = %10.3e\n", i, j,
              TacsRealPart(Q[3 * i] * Q[3 * j] + Q[3 * i + 1] * Q[3 * j + 1] +
                           Q[3 * i + 2] * Q[3 * j + 2]));
    }
  }

  for (int k = 0; k < 3; k++) {
    TacsScalar Cd[9], Ctd[27], Cttd[54];
    c1 = cos(U[0]);
    s1 = sin(U[0]);
    c2 = cos(U[1]);
    s2 = sin(U[1]);
    c3 = cos(U[2]);
    s3 = sin(U[2]);

    if (k == 0) {
      c1 = cos(U[0] + dh);
      s1 = sin(U[0] + dh);
    }
    if (k == 1) {
      c2 = cos(U[1] + dh);
      s2 = sin(U[1] + dh);
    }
    if (k == 2) {
      c3 = cos(U[2] + dh);
      s3 = sin(U[2] + dh);
    }

    compute_rate_matrix(Cd, Ctd, c1, s1, c2, s2, c3, s3);
    compute_2nd_rate_matrix(Cttd, c1, s1, c2, s2, c3, s3);

    for (int j = 0; j < 9; j++) {
      Cd[j] = (Cd[j] - C[j]) / dh;
    }

    for (int j = 0; j < 27; j++) {
      Ctd[j] = (Ctd[j] - Ct[j]) / dh;
    }

    for (int j = 0; j < 54; j++) {
      Cttd[j] = (Cttd[j] - Ctt[j]) / dh;
    }

    fprintf(stderr, "Ct%d[   ] %15s %15s %15s\n", k + 1, "Analytic",
            "Approximate", "Rel. Error");
    for (int j = 0; j < 9; j++) {
      if (Cd[j] != 0.0) {
        fprintf(stderr, "Ct%d[%3d] %15.6e %15.6e %15.4e\n", k + 1, j,
                TacsRealPart(Ct[9 * k + j]), TacsRealPart(Cd[j]),
                fabs(TacsRealPart((Ct[9 * k + j] - Cd[j]) / Cd[j])));
      } else {
        fprintf(stderr, "Ct%d[%3d] %15.6e %15.6e\n", k + 1, j,
                TacsRealPart(Ct[9 * k + j]), TacsRealPart(Cd[j]));
      }
    }

    int nt = 2;
    int t_index[2] = {0, 0};
    int tt_index[2] = {0, 0};

    if (k == 0) {
      t_index[0] = 0;
      tt_index[0] = 0;
      t_index[1] = 9;
      tt_index[1] = 27;
    }
    if (k == 1) {
      t_index[0] = 9;
      tt_index[0] = 9;
      t_index[1] = 18;
      tt_index[1] = 45;
    }
    if (k == 2) {
      t_index[0] = 18;
      tt_index[0] = 18;
      t_index[1] = 0;
      tt_index[1] = 36;
    }

    for (int p = 0; p < nt; p++) {
      int t = t_index[p];
      int tt = tt_index[p];

      fprintf(stderr, "Ctt%d[   ] %15s %15s %15s\n", k + 1, "Analytic",
              "Approximate", "Rel. Error");
      for (int j = 0; j < 9; j++) {
        if (Ctd[t + j] != 0.0) {
          fprintf(stderr, "Ctt%d[%3d] %15.6e %15.6e %15.4e\n", k + 1, j,
                  TacsRealPart(Ctt[tt + j]), TacsRealPart(Ctd[t + j]),
                  fabs(TacsRealPart((Ctt[tt + j] - Ctd[t + j]) / Ctd[t + j])));
        } else {
          fprintf(stderr, "Ctt%d[%3d] %15.6e %15.6e\n", k + 1, j,
                  TacsRealPart(Ctt[tt + j]), TacsRealPart(Ctd[t + j]));
        }
      }
    }

    if (k >= 1) {
      int nttt = 0;
      int tt_idx[5];
      int ttt_idx[5];

      if (k == 1) {
        nttt = 2;
        ttt_idx[0] = 0;
        tt_idx[0] = 0;  // C,112
        ttt_idx[1] = 18;
        tt_idx[1] = 27;  // C,122
      }
      if (k == 2) {
        nttt = 5;
        ttt_idx[0] = 9;
        tt_idx[0] = 0;  // C,113
        ttt_idx[1] = 27;
        tt_idx[1] = 27;  // C,123
        ttt_idx[2] = 36;
        tt_idx[2] = 36;  // C,133
        ttt_idx[3] = 45;
        tt_idx[3] = 9;  // C,223
        ttt_idx[4] = 54;
        tt_idx[4] = 45;  // C,233
      }

      for (int p = 0; p < nttt; p++) {
        int tt = tt_idx[p];
        int ttt = ttt_idx[p];

        fprintf(stderr, "Cttt%d[   ] %15s %15s %15s\n", k + 1, "Analytic",
                "Approximate", "Rel. Error");
        for (int j = 0; j < 9; j++) {
          if (Cttd[tt + j] != 0.0) {
            fprintf(stderr, "Cttt%d[%3d] %15.6e %15.6e %15.4e\n", k + 1, j,
                    TacsRealPart(Cttt[ttt + j]), TacsRealPart(Cttd[tt + j]),
                    fabs(TacsRealPart((Cttt[ttt + j] - Cttd[tt + j]) /
                                      Cttd[tt + j])));
          } else {
            fprintf(stderr, "Cttt%d[%3d] %15.6e %15.6e\n", k + 1, j,
                    TacsRealPart(Cttt[ttt + j]), TacsRealPart(Cttd[tt + j]));
          }
        }
      }
    }
  }
}

/*
  Evaluate the in-plane rotation penalty term

  input:
  num_points: the number of nodes
  Xd: the in-plane normal directions in the global frame
  Ud: the displacement in the global ref. frame
  C = Q^{T} - I: the normal-displacement matrix
  Ct: the derivative of C w.r.t. the theta values
  N, Na, Nb: the shape functions

  output:
  drot: the derivative of the penalty term w.r.t the variables

  returns:
  rot: the penalty term

  The rotation penalty term is calculated as follows:

  rot = X,1^{T}*Q*(X,2 + U,2) - X,2^{T}*Q*(X,1 + U,1)

  where the comma notation is used to denote differentiation.

  In python notation, the rot term is computed as follows

  rot = (dot(Xd[0:3], dot(C^{T} + I, Xd[3:6] + Ud[1:6:2])) -
         dot(Xd[3:6], dot(C^{T} + I, Xd[0:3] + Ud[0:6:2])))
*/
TacsScalar compute_inplane_penalty(TacsScalar drot[], const int num_points,
                                   const TacsScalar Xd[], const TacsScalar Ud[],
                                   const TacsScalar C[], const TacsScalar Ct[],
                                   const double N[], const double Na[],
                                   const double Nb[]) {
  TacsScalar v1[3], v2[3];

  // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
  v1[0] = (C[0] + 1.0) * (Xd[0] + Ud[0]) + C[3] * (Xd[1] + Ud[2]) +
          C[6] * (Xd[2] + Ud[4]);
  v1[1] = C[1] * (Xd[0] + Ud[0]) + (C[4] + 1.0) * (Xd[1] + Ud[2]) +
          C[7] * (Xd[2] + Ud[4]);
  v1[2] = C[2] * (Xd[0] + Ud[0]) + C[5] * (Xd[1] + Ud[2]) +
          (C[8] + 1.0) * (Xd[2] + Ud[4]);

  // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
  v2[0] = (C[0] + 1.0) * (Xd[3] + Ud[1]) + C[3] * (Xd[4] + Ud[3]) +
          C[6] * (Xd[5] + Ud[5]);
  v2[1] = C[1] * (Xd[3] + Ud[1]) + (C[4] + 1.0) * (Xd[4] + Ud[3]) +
          C[7] * (Xd[5] + Ud[5]);
  v2[2] = C[2] * (Xd[3] + Ud[1]) + C[5] * (Xd[4] + Ud[3]) +
          (C[8] + 1.0) * (Xd[5] + Ud[5]);

  // Compute the penalty term
  TacsScalar rot = ((Xd[0] * v2[0] + Xd[1] * v2[1] + Xd[2] * v2[2]) -
                    (Xd[3] * v1[0] + Xd[4] * v1[1] + Xd[5] * v1[2]));

  // Compute the penalty terms
  for (int i = 0; i < num_points; i++) {
    // Compute the derivative w.r.t. the in-plane displacements

    for (int j = 0; j < 3; j++) {
      TacsScalar dv1[3], dv2[3];
      // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
      dv1[0] = C[3 * j] * Na[0];
      dv1[1] = C[3 * j + 1] * Na[0];
      dv1[2] = C[3 * j + 2] * Na[0];

      dv1[j] += Na[0];

      // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
      dv2[0] = C[3 * j] * Nb[0];
      dv2[1] = C[3 * j + 1] * Nb[0];
      dv2[2] = C[3 * j + 2] * Nb[0];

      dv2[j] += Nb[0];

      drot[j] = ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                 (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
    }

    for (int j = 0; j < 3; j++) {
      const TacsScalar* Cd = &Ct[9 * j];
      TacsScalar dv1[3], dv2[3];
      // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
      dv1[0] = Cd[0] * (Xd[0] + Ud[0]) + Cd[3] * (Xd[1] + Ud[2]) +
               Cd[6] * (Xd[2] + Ud[4]);
      dv1[1] = Cd[1] * (Xd[0] + Ud[0]) + Cd[4] * (Xd[1] + Ud[2]) +
               Cd[7] * (Xd[2] + Ud[4]);
      dv1[2] = Cd[2] * (Xd[0] + Ud[0]) + Cd[5] * (Xd[1] + Ud[2]) +
               Cd[8] * (Xd[2] + Ud[4]);

      // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
      dv2[0] = Cd[0] * (Xd[3] + Ud[1]) + Cd[3] * (Xd[4] + Ud[3]) +
               Cd[6] * (Xd[5] + Ud[5]);
      dv2[1] = Cd[1] * (Xd[3] + Ud[1]) + Cd[4] * (Xd[4] + Ud[3]) +
               Cd[7] * (Xd[5] + Ud[5]);
      dv2[2] = Cd[2] * (Xd[3] + Ud[1]) + Cd[5] * (Xd[4] + Ud[3]) +
               Cd[8] * (Xd[5] + Ud[5]);

      drot[3 + j] = N[0] * ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                            (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
    }

    N++, Na++, Nb++, drot += 6;
  }

  return rot;
}

/*
  Add the second derivative of the in-plane penalty term to the
  stiffness matrix.

  input:
  scale: scale the second derivative values by this scalar
  num_points: the number of nodes
  Xd: the in-plane normal directions in the global frame
  Ud: the displacement in the global ref. frame
  C = Q^{T} - I: the normal-displacement matrix
  Ct: the derivative of C w.r.t. the theta values
  N: the shape functions

  output:
  matrix: adds additional terms

  The rotation penalty term is calculated as follows:

  rot = X,1^{T}*Q*(X,2 + U,2) - X,2^{T}*Q*(X,1 + U,1)

  where the comma notation is used to denote differentiation.

  In python notation, the rot term is computed as follows

  rot = (dot(Xd[0:3], dot(C^{T} + I, Xd[3:6] + Ud[1:6:2])) -
         dot(Xd[3:6], dot(C^{T} + I, Xd[0:3] + Ud[0:6:2])))
*/
void add_inplane_penalty(TacsScalar matrix[], const int num_points,
                         const TacsScalar scale, const TacsScalar Xd[],
                         const TacsScalar Ud[], const TacsScalar Ct[],
                         const TacsScalar Ctt[], const double N[],
                         const double Na[], const double Nb[]) {
  TacsScalar ddv1[3 * 6], ddv2[3 * 6];

  // Compute the inner products of the second derivatives and the directions
  // (X + U),1  and (X + U),2

  for (int i = 0; i < 6; i++) {
    const TacsScalar* C = &Ctt[9 * i];

    // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
    ddv1[3 * i] = C[0] * (Xd[0] + Ud[0]) + C[3] * (Xd[1] + Ud[2]) +
                  C[6] * (Xd[2] + Ud[4]);
    ddv1[1 + 3 * i] = C[1] * (Xd[0] + Ud[0]) + C[4] * (Xd[1] + Ud[2]) +
                      C[7] * (Xd[2] + Ud[4]);
    ddv1[2 + 3 * i] = C[2] * (Xd[0] + Ud[0]) + C[5] * (Xd[1] + Ud[2]) +
                      C[8] * (Xd[2] + Ud[4]);

    // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
    ddv2[3 * i] = C[0] * (Xd[3] + Ud[1]) + C[3] * (Xd[4] + Ud[3]) +
                  C[6] * (Xd[5] + Ud[5]);
    ddv2[1 + 3 * i] = C[1] * (Xd[3] + Ud[1]) + C[4] * (Xd[4] + Ud[3]) +
                      C[7] * (Xd[5] + Ud[5]);
    ddv2[2 + 3 * i] = C[2] * (Xd[3] + Ud[1]) + C[5] * (Xd[4] + Ud[3]) +
                      C[8] * (Xd[5] + Ud[5]);
  }

  // Compute the penalty terms
  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 3; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Loop from the start of the i-th row to the diagonal
      for (int j = 0; j <= i; j++) {
        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        mat += 3;

        for (int jj = 3; jj < end; jj++) {
          TacsScalar dv1[3], dv2[3];
          const TacsScalar* C = &Ct[9 * (jj - 3)];

          // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
          dv1[0] = C[3 * ii] * Na[i];
          dv1[1] = C[3 * ii + 1] * Na[i];
          dv1[2] = C[3 * ii + 2] * Na[i];

          // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
          dv2[0] = C[3 * ii] * Nb[i];
          dv2[1] = C[3 * ii + 1] * Nb[i];
          dv2[2] = C[3 * ii + 2] * Nb[i];

          mat[0] += scale * N[j] *
                    ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                     (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          mat += 1;
        }
      }
    }

    // Loop over only the rotation terms
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Loop from the start of the i-th row to the diagonal
      for (int j = 0; j <= i; j++) {
        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        for (int jj = 0; (jj < 3) && (jj < end); jj++) {
          TacsScalar dv1[3], dv2[3];
          const TacsScalar* C = &Ct[9 * (ii - 3)];

          // Compute v1 = dot(C^{T} + I, Xd[0:3] + Ud[0:6:2]))
          dv1[0] = C[3 * jj] * Na[j];
          dv1[1] = C[3 * jj + 1] * Na[j];
          dv1[2] = C[3 * jj + 2] * Na[j];

          // Compute v2 = dot(C^{T} + I, Xd[3:6] + Ud[1:6:2]))
          dv2[0] = C[3 * jj] * Nb[j];
          dv2[1] = C[3 * jj + 1] * Nb[j];
          dv2[2] = C[3 * jj + 2] * Nb[j];

          mat[0] += scale * N[i] *
                    ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                     (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          mat++;
        }

        for (int jj = 3; jj < end; jj++) {
          if ((ii == 3 && jj == 4) || (ii == 4 && jj == 3)) {
            // Compute the derivative C,12
            const TacsScalar* dv1 = &ddv1[3 * 3];
            const TacsScalar* dv2 = &ddv2[3 * 3];

            mat[0] += scale * N[i] * N[j] *
                      ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                       (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          } else if ((ii == 3 && jj == 5) || (ii == 5 && jj == 3)) {
            // Use the derivatives C,13
            const TacsScalar* dv1 = &ddv1[3 * 4];
            const TacsScalar* dv2 = &ddv2[3 * 4];

            mat[0] += scale * N[i] * N[j] *
                      ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                       (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          } else if ((ii == 4 && jj == 5) || (ii == 5 && jj == 4)) {
            // Use the derivatives C,23
            const TacsScalar* dv1 = &ddv1[3 * 5];
            const TacsScalar* dv2 = &ddv2[3 * 5];

            mat[0] += scale * N[i] * N[j] *
                      ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                       (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          } else {  // ii == jj
            // Use the diagonal
            const TacsScalar* dv1 = &ddv1[3 * (ii - 3)];
            const TacsScalar* dv2 = &ddv2[3 * (ii - 3)];

            mat[0] += scale * N[i] * N[j] *
                      ((Xd[0] * dv2[0] + Xd[1] * dv2[1] + Xd[2] * dv2[2]) -
                       (Xd[3] * dv1[0] + Xd[4] * dv1[1] + Xd[5] * dv1[2]));
          }

          mat++;
        }
      }
    }
  }
}

/*
  Determine the rate of change of the displacements through the shell,
  r, and the derivative of r, w.r.t. the coordinates xi.

  r    = C*n
  r,xi = C,xi*n + C*n,xi

  C,xi = C,j*n * theta_{j},xi

  output:
  r_xi: the derivative of r w.r.t. xi
  r_eta: the derivative of r w.r.t. eta

  input:
  Uxd: the rate of change of the displacements + rotation variables
  along the coordinate lines of the shell

  n_xi: the derivative of n w.r.t. xi
  n_eta: the derivative of n w.r.t. eta

  Cn: the vector Cn = C*n = (Q^{T} - I)*n
  Cn_xi: the vector Cn_xi = C*n_xi
  Cn_eta: the vector Cn_eta = C*n_eta
  C1n, C2n, C3n: Cjn = Cj*n the derivative of C w.r.t. theta
*/
static inline void compute_normal_rate(
    TacsScalar r_xi[], TacsScalar r_eta[], const TacsScalar Uxd[],
    const TacsScalar Cn_xi[], const TacsScalar Cn_eta[], const TacsScalar C1n[],
    const TacsScalar C2n[], const TacsScalar C3n[]) {
  r_xi[0] = Cn_xi[0] + C1n[0] * Uxd[6] + C2n[0] * Uxd[8] + C3n[0] * Uxd[10];
  r_xi[1] = Cn_xi[1] + C1n[1] * Uxd[6] + C2n[1] * Uxd[8] + C3n[1] * Uxd[10];
  r_xi[2] = Cn_xi[2] + C1n[2] * Uxd[6] + C2n[2] * Uxd[8] + C3n[2] * Uxd[10];

  r_eta[0] = Cn_eta[0] + C1n[0] * Uxd[7] + C2n[0] * Uxd[9] + C3n[0] * Uxd[11];
  r_eta[1] = Cn_eta[1] + C1n[1] * Uxd[7] + C2n[1] * Uxd[9] + C3n[1] * Uxd[11];
  r_eta[2] = Cn_eta[2] + C1n[2] * Uxd[7] + C2n[2] * Uxd[9] + C3n[2] * Uxd[11];
}

/*
  Compute the derivative of the normal rate of change of displacements
  through the thickness w.r.t. one of the thetas

  output:
  r: The derivative of r w.r.t. theta
  r_xi: The derivative of r along xi direction, w.r.t. theta
  r_eta: The derivative of r along eta direction w.r.t. theta

  input:
  Uxd: the rate of change of the displacements + rotation variables
  along the coordinate lines of the shell

  Cn_xi: the derivative of the surface normal along xi
  Cn_eta: the derivative of the surface normal along eta

  Ctn: the derivative of C*n = (Q^{T} - I)*n
  C1tn, C2tn, C3tn: the second derivatives of C w.r.t. theta
*/
static inline void compute_normal_rate_theta(
    TacsScalar dr_xi[], TacsScalar dr_eta[], const TacsScalar Uxd[], double N,
    double Na, double Nb, const TacsScalar Ctn[], const TacsScalar Ctn_xi[],
    const TacsScalar Ctn_eta[], const TacsScalar C1tn[],
    const TacsScalar C2tn[], const TacsScalar C3tn[]) {
  dr_xi[0] = Na * Ctn[0] + N * (Ctn_xi[0] + C1tn[0] * Uxd[6] +
                                C2tn[0] * Uxd[8] + C3tn[0] * Uxd[10]);
  dr_xi[1] = Na * Ctn[1] + N * (Ctn_xi[1] + C1tn[1] * Uxd[6] +
                                C2tn[1] * Uxd[8] + C3tn[1] * Uxd[10]);
  dr_xi[2] = Na * Ctn[2] + N * (Ctn_xi[2] + C1tn[2] * Uxd[6] +
                                C2tn[2] * Uxd[8] + C3tn[2] * Uxd[10]);

  dr_eta[0] = Nb * Ctn[0] + N * (Ctn_eta[0] + C1tn[0] * Uxd[7] +
                                 C2tn[0] * Uxd[9] + C3tn[0] * Uxd[11]);
  dr_eta[1] = Nb * Ctn[1] + N * (Ctn_eta[1] + C1tn[1] * Uxd[7] +
                                 C2tn[1] * Uxd[9] + C3tn[1] * Uxd[11]);
  dr_eta[2] = Nb * Ctn[2] + N * (Ctn_eta[2] + C1tn[2] * Uxd[7] +
                                 C2tn[2] * Uxd[9] + C3tn[2] * Uxd[11]);
}

/*
  Compute the second derivative of the terms r_xi[] and r_eta[]
  w.r.t. the rotations.

  output:
  ddr_xi: The derivative of r along xi direction, w.r.t. theta_i, theta_j
  ddr_eta: The derivative of r along eta direction w.r.t. theta_i, theta_j

  input:
  Uxd: the rate of change of the displacements + rotation variables
  along the coordinate lines of the shell

  Cijn_xi: the derivative of the surface normal along xi
  Cijn_eta: the derivative of the surface normal along eta

  Cij1n, Cij2n, Cij3n: the third derivatives of C w.r.t. theta
*/
static inline void compute_2nd_normal_rate_theta(
    TacsScalar ddr_xi[], TacsScalar ddr_eta[], const TacsScalar Uxd[],
    double Ni, double Nai, double Nbi, double Nj, double Naj, double Nbj,
    const TacsScalar Cijn[], const TacsScalar Cijn_xi[],
    const TacsScalar Cijn_eta[], const TacsScalar C1ijn[],
    const TacsScalar C2ijn[], const TacsScalar C3ijn[]) {
  ddr_xi[0] = (Nai * Nj + Naj * Ni) * Cijn[0] +
              Ni * Nj *
                  (Cijn_xi[0] + C1ijn[0] * Uxd[6] + C2ijn[0] * Uxd[8] +
                   C3ijn[0] * Uxd[10]);
  ddr_xi[1] = (Nai * Nj + Naj * Ni) * Cijn[1] +
              Ni * Nj *
                  (Cijn_xi[1] + C1ijn[1] * Uxd[6] + C2ijn[1] * Uxd[8] +
                   C3ijn[1] * Uxd[10]);
  ddr_xi[2] = (Nai * Nj + Naj * Ni) * Cijn[2] +
              Ni * Nj *
                  (Cijn_xi[2] + C1ijn[2] * Uxd[6] + C2ijn[2] * Uxd[8] +
                   C3ijn[2] * Uxd[10]);

  ddr_eta[0] = (Nbi * Nj + Nbj * Ni) * Cijn[0] +
               Ni * Nj *
                   (Cijn_eta[0] + C1ijn[0] * Uxd[7] + C2ijn[0] * Uxd[9] +
                    C3ijn[0] * Uxd[11]);
  ddr_eta[1] = (Nbi * Nj + Nbj * Ni) * Cijn[1] +
               Ni * Nj *
                   (Cijn_eta[1] + C1ijn[1] * Uxd[7] + C2ijn[1] * Uxd[9] +
                    C3ijn[1] * Uxd[11]);
  ddr_eta[2] = (Nbi * Nj + Nbj * Ni) * Cijn[2] +
               Ni * Nj *
                   (Cijn_eta[2] + C1ijn[2] * Uxd[7] + C2ijn[2] * Uxd[9] +
                    C3ijn[2] * Uxd[11]);
}

/*
  Compute the nonlinear components of the strain using the large angle
  shell formulation.

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
*/
void large_rot_strain(TacsScalar strain[], const TacsScalar Ux[],
                      const TacsScalar Uxd[], const TacsScalar C[],
                      const TacsScalar Ct[], const TacsScalar t[],
                      const TacsScalar tx[], const TacsScalar ztx[],
                      const TacsScalar n[], const TacsScalar n_xi[],
                      const TacsScalar n_eta[]) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Transform the displacement into the local shell coordinates
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
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
  Compute the derivative of the nonlinear components of the strain
  using the large angle shell formulation.

  This function uses the same expression for the displacement
  gradient:

  U_{e,e} = t * U_{x,xi} * tx^{T}

  These are then used to construct the nonlinear expressions for the
  strain in the local Cartesian coordinate frame e.

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell
  C, Ct: the normal rotation rate matrices

  t, dt: the transformation and the derivative
  tx, dtx: the transformation times the Jacobian and its derivative
  ztx, dztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction and its derivative

  n, dn: the shell normal and the derivative of the shell normal
  n_xi, dn_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  and its derivative
  n_eta, dn_eta: the derivative of the shell normal w.r.t. the coordinate line
  eta and its derivative

  num_componnets: the number of components to differentiate

  output:
  strain: the value of the strain at the current point
  dstrain: the value of the derivative of the strain at the current point
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
                           const int num_components) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Transform the displacement into the local shell coordinates
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
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

  // Now go through for each component of the derivative and compute the
  // derivative of the strain w.r.t. the input
  for (int i = 0; i < num_components; i++) {
    // Transform the displacement into the local shell coordinates
    TacsScalar dUd[6];
    transform_displ_gradient(dUd, dt, Uxd);

    // Compute the rate of change of the displacement through
    // the thickness
    TacsScalar dr[3];  // Compute dr = C*dn
    transform_vector3d(dr, dn, C);

    // Compute the product of the derivatives of n with C
    TacsScalar dCn_xi[3], dCn_eta[3];
    transform_vector3d(dCn_xi, dn_xi, C);
    transform_vector3d(dCn_eta, dn_eta, C);

    // Compute the product of the derivatives of C with n
    TacsScalar dC1n[3], dC2n[3], dC3n[3];
    transform_vector3d(dC1n, dn, C1);
    transform_vector3d(dC2n, dn, C2);
    transform_vector3d(dC3n, dn, C3);

    // Recompute r, r_xi and r_eta before transformation
    // This is required for the derivative operations
    transform_vector3d(r, n, C);
    compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

    // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
    TacsScalar dr_xi[3], dr_eta[3];
    compute_normal_rate(dr_xi, dr_eta, Uxd, dCn_xi, dCn_eta, dC1n, dC2n, dC3n);

    // Transform the displacement rate through the thickness into the
    // local shell coordinates
    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface value of the displacement derivatives
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

    dt += 9;
    dtx += 9;
    dztx += 9;
    dn += 3;
    dn_xi += 3;
    dn_eta += 3;
  }
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
*/
void large_rot_bmat(TacsScalar B[], const int num_points, const double N[],
                    const double Na[], const double Nb[], const TacsScalar Ux[],
                    const TacsScalar Uxd[], const TacsScalar C[],
                    const TacsScalar Ct[], const TacsScalar Ctt[],
                    const TacsScalar t[], const TacsScalar tx[],
                    const TacsScalar ztx[], const TacsScalar n[],
                    const TacsScalar n_xi[], const TacsScalar n_eta[]) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Extract the second derivative information
  const TacsScalar* C11 = &Ctt[0];
  const TacsScalar* C22 = &Ctt[9];
  const TacsScalar* C33 = &Ctt[18];
  const TacsScalar* C12 = &Ctt[27];
  const TacsScalar* C13 = &Ctt[36];
  const TacsScalar* C23 = &Ctt[45];

  // Transform the displacement gradient to the local frame
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
  transform_vector3d(r, t);
  transform_vector3d(Cn_xi, t);
  transform_vector3d(Cn_eta, t);
  transform_vector3d(C1n, t);
  transform_vector3d(C2n, t);
  transform_vector3d(C3n, t);

  // Compute the derivatives of C times n,xi and n,eta
  TacsScalar C1n_xi[3], C2n_xi[3], C3n_xi[3];
  TacsScalar C1n_eta[3], C2n_eta[3], C3n_eta[3];
  transform_vector3d(C1n_xi, n_xi, C1);
  transform_vector3d(C2n_xi, n_xi, C2);
  transform_vector3d(C3n_xi, n_xi, C3);
  transform_vector3d(C1n_eta, n_eta, C1);
  transform_vector3d(C2n_eta, n_eta, C2);
  transform_vector3d(C3n_eta, n_eta, C3);

  // Transform the vectors to the local coordinate system
  transform_vector3d(C1n_xi, t);
  transform_vector3d(C2n_xi, t);
  transform_vector3d(C3n_xi, t);
  transform_vector3d(C1n_eta, t);
  transform_vector3d(C2n_eta, t);
  transform_vector3d(C3n_eta, t);

  // Compute the product of the second derivatives of C with n
  TacsScalar C11n[3], C22n[3], C33n[3], C12n[3], C13n[3], C23n[3];
  transform_vector3d(C11n, n, C11);
  transform_vector3d(C22n, n, C22);
  transform_vector3d(C33n, n, C33);
  transform_vector3d(C12n, n, C12);
  transform_vector3d(C13n, n, C13);
  transform_vector3d(C23n, n, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n, t);
  transform_vector3d(C22n, t);
  transform_vector3d(C33n, t);
  transform_vector3d(C12n, t);
  transform_vector3d(C13n, t);
  transform_vector3d(C23n, t);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

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
    // For the displacement components
    for (int ii = 0; ii < 3; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1];
      TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4];
      TacsScalar duz = Ud[0] * tx[6] + Ud[1] * tx[7];

      TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1];
      TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4];
      TacsScalar dvz = Ud[2] * tx[6] + Ud[3] * tx[7];

      TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1];
      TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4];
      TacsScalar dwz = Ud[4] * tx[6] + Ud[5] * tx[7];

      // The first-derivative values of the displacements
      TacsScalar dux1 = Ud[0] * ztx[0] + Ud[1] * ztx[1];
      TacsScalar duy1 = Ud[0] * ztx[3] + Ud[1] * ztx[4];

      TacsScalar dvx1 = Ud[2] * ztx[0] + Ud[3] * ztx[1];
      TacsScalar dvy1 = Ud[2] * ztx[3] + Ud[3] * ztx[4];

      TacsScalar dwx1 = Ud[4] * ztx[0] + Ud[5] * ztx[1];
      TacsScalar dwy1 = Ud[4] * ztx[3] + Ud[5] * ztx[4];

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
    }

    // For each displacement component
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      if (ii == 3) {
        dr[0] = C1n[0] * N[i];
        dr[1] = C1n[1] * N[i];
        dr[2] = C1n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                  C1n_xi, C1n_eta, C11n, C12n, C13n);
      } else if (ii == 4) {
        dr[0] = C2n[0] * N[i];
        dr[1] = C2n[1] * N[i];
        dr[2] = C2n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                  C2n_xi, C2n_eta, C12n, C22n, C23n);
      } else {  // ii == 5
        dr[0] = C3n[0] * N[i];
        dr[1] = C3n[1] * N[i];
        dr[2] = C3n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                  C3n_xi, C3n_eta, C13n, C23n, C33n);
      }

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = dr[0] * tx[2];
      TacsScalar duy = dr[0] * tx[5];
      TacsScalar duz = dr[0] * tx[8];

      TacsScalar dvx = dr[1] * tx[2];
      TacsScalar dvy = dr[1] * tx[5];
      TacsScalar dvz = dr[1] * tx[8];

      TacsScalar dwx = dr[2] * tx[2];
      TacsScalar dwy = dr[2] * tx[5];
      TacsScalar dwz = dr[2] * tx[8];

      // The first-derivative values of the displacements
      TacsScalar dux1 = dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dr[0] * ztx[2];
      TacsScalar duy1 = dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dr[0] * ztx[5];

      TacsScalar dvx1 = dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dr[1] * ztx[2];
      TacsScalar dvy1 = dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dr[1] * ztx[5];

      TacsScalar dwx1 = dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dr[2] * ztx[2];
      TacsScalar dwy1 = dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dr[2] * ztx[5];

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
    }
  }
}

/*
  Compute the derivative w.r.t. the nodes of the derivative of the
  nonlinear strain expressions w.r.t the nodal
  displacements/rotations. I.e. compute d(stress^{T}B)/dx, where x
  are the set of nodes.

  input:
  N: the shape functions
  Na, Nb: the derivatives of the shape functions w.r.t. the natural
  coordinates of the shell

  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell

  C, Ct, Ctt: the rotational rate matrices and their derivatives

  t, dt: the transformation and its derivative
  tx: the transformation times the Jacobian and its derivative
  ztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction and its derivative
  n: the shell normal and its derivative
  n_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  n_eta: the derivative of the shell normal w.r.t. the coordinate line eta

  output:
  res: the matrix to which the product is added
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
    const TacsScalar dn_eta[], const int num_components) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Extract the second derivative information
  const TacsScalar* C11 = &Ctt[0];
  const TacsScalar* C22 = &Ctt[9];
  const TacsScalar* C33 = &Ctt[18];
  const TacsScalar* C12 = &Ctt[27];
  const TacsScalar* C13 = &Ctt[36];
  const TacsScalar* C23 = &Ctt[45];

  // Transform the displacement gradient to the local frame
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
  transform_vector3d(r, t);
  transform_vector3d(Cn_xi, t);
  transform_vector3d(Cn_eta, t);
  transform_vector3d(C1n, t);
  transform_vector3d(C2n, t);
  transform_vector3d(C3n, t);

  // Compute the derivatives of C times n,xi and n,eta
  TacsScalar C1n_xi[3], C2n_xi[3], C3n_xi[3];
  TacsScalar C1n_eta[3], C2n_eta[3], C3n_eta[3];
  transform_vector3d(C1n_xi, n_xi, C1);
  transform_vector3d(C2n_xi, n_xi, C2);
  transform_vector3d(C3n_xi, n_xi, C3);
  transform_vector3d(C1n_eta, n_eta, C1);
  transform_vector3d(C2n_eta, n_eta, C2);
  transform_vector3d(C3n_eta, n_eta, C3);

  // Transform the vectors to the local coordinate system
  transform_vector3d(C1n_xi, t);
  transform_vector3d(C2n_xi, t);
  transform_vector3d(C3n_xi, t);
  transform_vector3d(C1n_eta, t);
  transform_vector3d(C2n_eta, t);
  transform_vector3d(C3n_eta, t);

  // Compute the product of the second derivatives of C with n
  TacsScalar C11n[3], C22n[3], C33n[3], C12n[3], C13n[3], C23n[3];
  transform_vector3d(C11n, n, C11);
  transform_vector3d(C22n, n, C22);
  transform_vector3d(C33n, n, C33);
  transform_vector3d(C12n, n, C12);
  transform_vector3d(C13n, n, C13);
  transform_vector3d(C23n, n, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n, t);
  transform_vector3d(C22n, t);
  transform_vector3d(C33n, t);
  transform_vector3d(C12n, t);
  transform_vector3d(C13n, t);
  transform_vector3d(C23n, t);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

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

  // Store the derivatives of the displacements w.r.t. each node
  TacsScalar dU0[9 * 6 * shellutils::MAX_NUM_NODES];
  TacsScalar dU1[6 * 6 * shellutils::MAX_NUM_NODES];

  // First, compute the derivatives of the displacements w.r.t. the
  // nodal coordinates
  TacsScalar* du0 = dU0;
  TacsScalar* du1 = dU1;

  // For each point
  for (int i = 0; i < num_points; i++) {
    // For the displacement components
    for (int ii = 0; ii < 3; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      // The mid-surface value of the displacement derivatives
      du0[0] = Ud[0] * tx[0] + Ud[1] * tx[1];
      du0[1] = Ud[0] * tx[3] + Ud[1] * tx[4];
      du0[2] = Ud[0] * tx[6] + Ud[1] * tx[7];

      du0[3] = Ud[2] * tx[0] + Ud[3] * tx[1];
      du0[4] = Ud[2] * tx[3] + Ud[3] * tx[4];
      du0[5] = Ud[2] * tx[6] + Ud[3] * tx[7];

      du0[6] = Ud[4] * tx[0] + Ud[5] * tx[1];
      du0[7] = Ud[4] * tx[3] + Ud[5] * tx[4];
      du0[8] = Ud[4] * tx[6] + Ud[5] * tx[7];
      du0 += 9;

      // The first-derivative values of the displacements
      du1[0] = Ud[0] * ztx[0] + Ud[1] * ztx[1];
      du1[1] = Ud[0] * ztx[3] + Ud[1] * ztx[4];

      du1[2] = Ud[2] * ztx[0] + Ud[3] * ztx[1];
      du1[3] = Ud[2] * ztx[3] + Ud[3] * ztx[4];

      du1[4] = Ud[4] * ztx[0] + Ud[5] * ztx[1];
      du1[5] = Ud[4] * ztx[3] + Ud[5] * ztx[4];
      du1 += 6;
    }

    // For each displacement component
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      if (ii == 3) {
        dr[0] = C1n[0] * N[i];
        dr[1] = C1n[1] * N[i];
        dr[2] = C1n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                  C1n_xi, C1n_eta, C11n, C12n, C13n);
      } else if (ii == 4) {
        dr[0] = C2n[0] * N[i];
        dr[1] = C2n[1] * N[i];
        dr[2] = C2n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                  C2n_xi, C2n_eta, C12n, C22n, C23n);
      } else {  // ii == 5
        dr[0] = C3n[0] * N[i];
        dr[1] = C3n[1] * N[i];
        dr[2] = C3n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                  C3n_xi, C3n_eta, C13n, C23n, C33n);
      }

      // The mid-surface value of the displacement derivatives
      du0[0] = dr[0] * tx[2];
      du0[1] = dr[0] * tx[5];
      du0[2] = dr[0] * tx[8];

      du0[3] = dr[1] * tx[2];
      du0[4] = dr[1] * tx[5];
      du0[5] = dr[1] * tx[8];

      du0[6] = dr[2] * tx[2];
      du0[7] = dr[2] * tx[5];
      du0[8] = dr[2] * tx[8];
      du0 += 9;

      // The first-derivative values of the displacements
      du1[0] = dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dr[0] * ztx[2];
      du1[1] = dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dr[0] * ztx[5];

      du1[2] = dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dr[1] * ztx[2];
      du1[3] = dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dr[1] * ztx[5];

      du1[4] = dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dr[2] * ztx[2];
      du1[5] = dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dr[2] * ztx[5];
      du1 += 6;
    }
  }

  for (int j = 0; j < num_components; j++) {
    // Transform the displacement into the local shell coordinates
    TacsScalar dUd[6];
    transform_displ_gradient(Ud, t, Uxd);
    transform_displ_gradient(dUd, dt, Uxd);

    // Compute the rate of change of the displacement through
    // the thickness
    TacsScalar dr[3];  // Compute dr = C*dn
    transform_vector3d(r, n, C);
    transform_vector3d(dr, dn, C);
    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d(r, t);

    // Compute the product of the derivatives of n with C
    TacsScalar dr_xi[3], dr_eta[3];
    TacsScalar dCn_xi[3], dCn_eta[3];
    transform_vector3d(Cn_xi, n_xi, C);
    transform_vector3d(Cn_eta, n_eta, C);
    transform_vector3d(dCn_xi, dn_xi, C);
    transform_vector3d(dCn_eta, dn_eta, C);
    transform_vector3d_sens(dCn_xi, Cn_xi, t, dt);
    transform_vector3d_sens(dCn_eta, Cn_eta, t, dt);
    transform_vector3d(Cn_xi, t);
    transform_vector3d(Cn_eta, t);

    // Compute the product of the derivatives of C with n
    TacsScalar dC1n[3], dC2n[3], dC3n[3];
    transform_vector3d(C1n, n, C1);
    transform_vector3d(C2n, n, C2);
    transform_vector3d(C3n, n, C3);
    transform_vector3d(dC1n, dn, C1);
    transform_vector3d(dC2n, dn, C2);
    transform_vector3d(dC3n, dn, C3);
    transform_vector3d_sens(dC1n, C1n, t, dt);
    transform_vector3d_sens(dC2n, C2n, t, dt);
    transform_vector3d_sens(dC3n, C3n, t, dt);
    transform_vector3d(C1n, t);
    transform_vector3d(C2n, t);
    transform_vector3d(C3n, t);

    // Compute the derivatives of C times n,xi and n,eta
    TacsScalar dC1n_xi[3], dC2n_xi[3], dC3n_xi[3];
    TacsScalar dC1n_eta[3], dC2n_eta[3], dC3n_eta[3];
    transform_vector3d(C1n_xi, n_xi, C1);
    transform_vector3d(C2n_xi, n_xi, C2);
    transform_vector3d(C3n_xi, n_xi, C3);
    transform_vector3d(C1n_eta, n_eta, C1);
    transform_vector3d(C2n_eta, n_eta, C2);
    transform_vector3d(C3n_eta, n_eta, C3);
    transform_vector3d(dC1n_xi, dn_xi, C1);
    transform_vector3d(dC2n_xi, dn_xi, C2);
    transform_vector3d(dC3n_xi, dn_xi, C3);
    transform_vector3d(dC1n_eta, dn_eta, C1);
    transform_vector3d(dC2n_eta, dn_eta, C2);
    transform_vector3d(dC3n_eta, dn_eta, C3);
    transform_vector3d_sens(dC1n_xi, C1n_xi, t, dt);
    transform_vector3d_sens(dC2n_xi, C2n_xi, t, dt);
    transform_vector3d_sens(dC3n_xi, C3n_xi, t, dt);
    transform_vector3d_sens(dC1n_eta, C1n_eta, t, dt);
    transform_vector3d_sens(dC2n_eta, C2n_eta, t, dt);
    transform_vector3d_sens(dC3n_eta, C3n_eta, t, dt);
    transform_vector3d(C1n_xi, t);
    transform_vector3d(C2n_xi, t);
    transform_vector3d(C3n_xi, t);
    transform_vector3d(C1n_eta, t);
    transform_vector3d(C2n_eta, t);
    transform_vector3d(C3n_eta, t);

    // Compute the product of the second derivatives of C with n
    TacsScalar dC11n[3], dC22n[3], dC33n[3], dC12n[3], dC13n[3], dC23n[3];
    transform_vector3d(C11n, n, C11);
    transform_vector3d(C22n, n, C22);
    transform_vector3d(C33n, n, C33);
    transform_vector3d(C12n, n, C12);
    transform_vector3d(C13n, n, C13);
    transform_vector3d(C23n, n, C23);
    transform_vector3d(dC11n, dn, C11);
    transform_vector3d(dC22n, dn, C22);
    transform_vector3d(dC33n, dn, C33);
    transform_vector3d(dC12n, dn, C12);
    transform_vector3d(dC13n, dn, C13);
    transform_vector3d(dC23n, dn, C23);
    transform_vector3d_sens(dC11n, C11n, t, dt);
    transform_vector3d_sens(dC22n, C22n, t, dt);
    transform_vector3d_sens(dC33n, C33n, t, dt);
    transform_vector3d_sens(dC12n, C12n, t, dt);
    transform_vector3d_sens(dC13n, C13n, t, dt);
    transform_vector3d_sens(dC23n, C23n, t, dt);
    transform_vector3d(C11n, t);
    transform_vector3d(C22n, t);
    transform_vector3d(C33n, t);
    transform_vector3d(C12n, t);
    transform_vector3d(C13n, t);
    transform_vector3d(C23n, t);

    // Recompute r, r_xi and r_eta before transformation
    // This is required for the derivative operations
    compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

    // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
    compute_normal_rate(dr_xi, dr_eta, Uxd, dCn_xi, dCn_eta, dC1n, dC2n, dC3n);

    // The mid-surface value of the displacement derivatives
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

    // Set the pointers to the pre-computed derivatives
    const TacsScalar* du0 = dU0;
    const TacsScalar* du1 = dU1;

    // For each point
    for (int i = 0; i < num_points; i++) {
      // For the displacement components
      for (int ii = 0; ii < 3; ii++) {
        transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);
        transform_displ_gradient_bmat(dUd, ii, dt, Na[i], Nb[i]);

        // The mid-surface values
        TacsScalar sdux =
            (dUd[0] * tx[0] + dUd[1] * tx[1] + Ud[0] * dtx[0] + Ud[1] * dtx[1]);
        TacsScalar sduy =
            (dUd[0] * tx[3] + dUd[1] * tx[4] + Ud[0] * dtx[3] + Ud[1] * dtx[4]);
        TacsScalar sduz =
            (dUd[0] * tx[6] + dUd[1] * tx[7] + Ud[0] * dtx[6] + Ud[1] * dtx[7]);

        TacsScalar sdvx =
            (dUd[2] * tx[0] + dUd[3] * tx[1] + Ud[2] * dtx[0] + Ud[3] * dtx[1]);
        TacsScalar sdvy =
            (dUd[2] * tx[3] + dUd[3] * tx[4] + Ud[2] * dtx[3] + Ud[3] * dtx[4]);
        TacsScalar sdvz =
            (dUd[2] * tx[6] + dUd[3] * tx[7] + Ud[2] * dtx[6] + Ud[3] * dtx[7]);

        TacsScalar sdwx =
            (dUd[4] * tx[0] + dUd[5] * tx[1] + Ud[4] * dtx[0] + Ud[5] * dtx[1]);
        TacsScalar sdwy =
            (dUd[4] * tx[3] + dUd[5] * tx[4] + Ud[4] * dtx[3] + Ud[5] * dtx[4]);
        TacsScalar sdwz =
            (dUd[4] * tx[6] + dUd[5] * tx[7] + Ud[4] * dtx[6] + Ud[5] * dtx[7]);

        // The first-derivative values of the displacements
        TacsScalar sdux1 = (dUd[0] * ztx[0] + dUd[1] * ztx[1] +
                            Ud[0] * dztx[0] + Ud[1] * dztx[1]);
        TacsScalar sduy1 = (dUd[0] * ztx[3] + dUd[1] * ztx[4] +
                            Ud[0] * dztx[3] + Ud[1] * dztx[4]);

        TacsScalar sdvx1 = (dUd[2] * ztx[0] + dUd[3] * ztx[1] +
                            Ud[2] * dztx[0] + Ud[3] * dztx[1]);
        TacsScalar sdvy1 = (dUd[2] * ztx[3] + dUd[3] * ztx[4] +
                            Ud[2] * dztx[3] + Ud[3] * dztx[4]);

        TacsScalar sdwx1 = (dUd[4] * ztx[0] + dUd[5] * ztx[1] +
                            Ud[4] * dztx[0] + Ud[5] * dztx[1]);
        TacsScalar sdwy1 = (dUd[4] * ztx[3] + dUd[5] * ztx[4] +
                            Ud[4] * dztx[3] + Ud[5] * dztx[4]);

        TacsScalar dB[8];

        // The in plane strains
        dB[0] = sdux + (sux * du0[0] + svx * du0[3] + swx * du0[6] + ux * sdux +
                        vx * sdvx + wx * sdwx);
        dB[1] = sdvy + (suy * du0[1] + svy * du0[4] + swy * du0[7] + uy * sduy +
                        vy * sdvy + wy * sdwy);
        dB[2] = sduy + sdvx +
                (sdux * uy + sdvx * vy + sdwx * wy + sux * du0[1] +
                 svx * du0[4] + swx * du0[7] + du0[0] * suy + du0[3] * svy +
                 du0[6] * swy + ux * sduy + vx * sdvy + wx * sdwy);

        // Compute the bending components of the strain
        dB[3] = sdux1 +
                (sdux * ux1 + sdvx * vx1 + sdwx * wx1 + sux * du1[0] +
                 svx * du1[3] + swx * du1[6] + du0[0] * sux1 + du0[3] * svx1 +
                 du0[6] * swx1 + ux * sdux1 + vx * sdvx1 + wx * sdwx1);
        dB[4] = sdvy1 +
                (sduy * uy1 + sdvy * vy1 + sdwy * wy1 + suy * du1[1] +
                 svy * du1[4] + swy * du1[7] + du0[1] * suy1 + du0[4] * svy1 +
                 du0[7] * swy1 + uy * sduy1 + vy * sdvy1 + wy * sdwy1);
        dB[5] = sduy1 + sdvx1 +
                (sdux1 * uy + sdvx1 * vy + sdwx1 * wy + sdux * uy1 +
                 sdvx * vy1 + sdwx * wy1 + sux1 * du0[1] + svx1 * du0[4] +
                 swx1 * du0[7] + sux * du1[1] + svx * du1[4] + swx * du1[7] +
                 du1[0] * suy + du1[3] * svy + du1[6] * swy + du0[0] * suy1 +
                 du0[3] * svy1 + du0[6] * swy1 + ux1 * sduy + vx1 * sdvy +
                 wx1 * sdwy + ux * sduy1 + vx * sdvy1 + wx * sdwy1);

        // Compute the shear terms
        dB[6] = sdvz + sdwy +
                (sduz * uy + sdvz * vy + sdwz * wy + suz * du0[1] +
                 svz * du0[4] + swz * du0[7] + du0[2] * suy + du0[5] * svy +
                 du0[8] * swy + uz * sduy + vz * sdvy + wz * sdwy);
        dB[7] = sduz + sdwx +
                (sduz * ux + sdvz * vx + sdwz * wx + suz * du0[0] +
                 svz * du0[3] + swz * du0[7] + du0[2] * sux + du0[5] * svx +
                 du0[8] * swx + uz * sdux + vz * sdvx + wz * sdwx);

        res[0] +=
            scale * (stress[0] * dB[0] + stress[1] * dB[1] + stress[2] * dB[2] +
                     stress[3] * dB[3] + stress[4] * dB[4] + stress[5] * dB[5] +
                     stress[6] * dB[6] + stress[7] * dB[7]);
        res += 1;

        du0 += 9;
        du1 += 6;
      }

      // For each displacement component
      for (int ii = 3; ii < 6; ii++) {
        // TacsScalar dr[3], dr_xi[3], dr_eta[3];
        TacsScalar sdr[3], sdr_xi[3], sdr_eta[3];
        if (ii == 3) {
          dr[0] = C1n[0] * N[i];
          dr[1] = C1n[1] * N[i];
          dr[2] = C1n[2] * N[i];
          compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                    C1n_xi, C1n_eta, C11n, C12n, C13n);

          sdr[0] = dC1n[0] * N[i];
          sdr[1] = dC1n[1] * N[i];
          sdr[2] = dC1n[2] * N[i];
          compute_normal_rate_theta(sdr_xi, sdr_eta, Uxd, N[i], Na[i], Nb[i],
                                    dC1n, dC1n_xi, dC1n_eta, dC11n, dC12n,
                                    dC13n);
        } else if (ii == 4) {
          dr[0] = C2n[0] * N[i];
          dr[1] = C2n[1] * N[i];
          dr[2] = C2n[2] * N[i];
          compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                    C2n_xi, C2n_eta, C12n, C22n, C23n);

          sdr[0] = dC2n[0] * N[i];
          sdr[1] = dC2n[1] * N[i];
          sdr[2] = dC2n[2] * N[i];
          compute_normal_rate_theta(sdr_xi, sdr_eta, Uxd, N[i], Na[i], Nb[i],
                                    dC2n, dC2n_xi, dC2n_eta, dC12n, dC22n,
                                    dC23n);
        } else {  // ii == 5
          dr[0] = C3n[0] * N[i];
          dr[1] = C3n[1] * N[i];
          dr[2] = C3n[2] * N[i];
          compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                    C3n_xi, C3n_eta, C13n, C23n, C33n);

          sdr[0] = dC3n[0] * N[i];
          sdr[1] = dC3n[1] * N[i];
          sdr[2] = dC3n[2] * N[i];
          compute_normal_rate_theta(sdr_xi, sdr_eta, Uxd, N[i], Na[i], Nb[i],
                                    dC3n, dC3n_xi, dC3n_eta, dC13n, dC23n,
                                    dC33n);
        }

        // The mid-surface value of the displacement derivatives
        TacsScalar sdux = sdr[0] * tx[2] + dr[0] * dtx[2];
        TacsScalar sduy = sdr[0] * tx[5] + dr[0] * dtx[5];
        TacsScalar sduz = sdr[0] * tx[8] + dr[0] * dtx[8];

        TacsScalar sdvx = sdr[1] * tx[2] + dr[1] * dtx[2];
        TacsScalar sdvy = sdr[1] * tx[5] + dr[1] * dtx[5];
        TacsScalar sdvz = sdr[1] * tx[8] + dr[1] * dtx[8];

        TacsScalar sdwx = sdr[2] * tx[2] + dr[2] * dtx[2];
        TacsScalar sdwy = sdr[2] * tx[5] + dr[2] * dtx[5];
        TacsScalar sdwz = sdr[2] * tx[8] + dr[2] * dtx[8];

        // The first-derivative values of the displacements
        TacsScalar sdux1 =
            (sdr_xi[0] * tx[0] + sdr_eta[0] * tx[1] + sdr[0] * ztx[2] +
             dr_xi[0] * dtx[0] + dr_eta[0] * dtx[1] + dr[0] * dztx[2]);
        TacsScalar sduy1 =
            (sdr_xi[0] * tx[3] + sdr_eta[0] * tx[4] + sdr[0] * ztx[5] +
             dr_xi[0] * dtx[3] + dr_eta[0] * dtx[4] + dr[0] * dztx[5]);

        TacsScalar sdvx1 =
            (sdr_xi[1] * tx[0] + sdr_eta[1] * tx[1] + sdr[1] * ztx[2] +
             dr_xi[1] * dtx[0] + dr_eta[1] * dtx[1] + dr[1] * dztx[2]);
        TacsScalar sdvy1 =
            (sdr_xi[1] * tx[3] + sdr_eta[1] * tx[4] + sdr[1] * ztx[5] +
             dr_xi[1] * dtx[3] + dr_eta[1] * dtx[4] + dr[1] * dztx[5]);

        TacsScalar sdwx1 =
            (sdr_xi[2] * tx[0] + sdr_eta[2] * tx[1] + sdr[2] * ztx[2] +
             dr_xi[2] * dtx[0] + dr_eta[2] * dtx[1] + dr[2] * dztx[2]);
        TacsScalar sdwy1 =
            (sdr_xi[2] * tx[3] + sdr_eta[2] * tx[4] + sdr[2] * ztx[5] +
             dr_xi[2] * dtx[3] + dr_eta[2] * dtx[4] + dr[2] * dztx[5]);

        TacsScalar dB[8];

        // The in plane strains
        dB[0] = sdux + (sux * du0[0] + svx * du0[3] + swx * du0[6] + ux * sdux +
                        vx * sdvx + wx * sdwx);
        dB[1] = sdvy + (suy * du0[1] + svy * du0[4] + swy * du0[7] + uy * sduy +
                        vy * sdvy + wy * sdwy);
        dB[2] = sduy + sdvx +
                (sdux * uy + sdvx * vy + sdwx * wy + sux * du0[1] +
                 svx * du0[4] + swx * du0[7] + du0[0] * suy + du0[3] * svy +
                 du0[6] * swy + ux * sduy + vx * sdvy + wx * sdwy);

        // Compute the bending components of the strain
        dB[3] = sdux1 +
                (sdux * ux1 + sdvx * vx1 + sdwx * wx1 + sux * du1[0] +
                 svx * du1[3] + swx * du1[6] + du0[0] * sux1 + du0[3] * svx1 +
                 du0[6] * swx1 + ux * sdux1 + vx * sdvx1 + wx * sdwx1);
        dB[4] = sdvy1 +
                (sduy * uy1 + sdvy * vy1 + sdwy * wy1 + suy * du1[1] +
                 svy * du1[4] + swy * du1[7] + du0[1] * suy1 + du0[4] * svy1 +
                 du0[7] * swy1 + uy * sduy1 + vy * sdvy1 + wy * sdwy1);
        dB[5] = sduy1 + sdvx1 +
                (sdux1 * uy + sdvx1 * vy + sdwx1 * wy + sdux * uy1 +
                 sdvx * vy1 + sdwx * wy1 + sux1 * du0[1] + svx1 * du0[4] +
                 swx1 * du0[7] + sux * du1[1] + svx * du1[4] + swx * du1[7] +
                 du1[0] * suy + du1[3] * svy + du1[6] * swy + du0[0] * suy1 +
                 du0[3] * svy1 + du0[6] * swy1 + ux1 * sduy + vx1 * sdvy +
                 wx1 * sdwy + ux * sduy1 + vx * sdvy1 + wx * sdwy1);

        // Compute the shear terms
        dB[6] = sdvz + sdwy +
                (sduz * uy + sdvz * vy + sdwz * wy + suz * du0[1] +
                 svz * du0[4] + swz * du0[7] + du0[2] * suy + du0[5] * svy +
                 du0[8] * swy + uz * sduy + vz * sdvy + wz * sdwy);
        dB[7] = sduz + sdwx +
                (sduz * ux + sdvz * vx + sdwz * wx + suz * du0[0] +
                 svz * du0[3] + swz * du0[7] + du0[2] * sux + du0[5] * svx +
                 du0[8] * swx + uz * sdux + vz * sdvx + wz * sdwx);

        res[0] +=
            scale * (stress[0] * dB[0] + stress[1] * dB[1] + stress[2] * dB[2] +
                     stress[3] * dB[3] + stress[4] * dB[4] + stress[5] * dB[5] +
                     stress[6] * dB[6] + stress[7] * dB[7]);
        res += 1;

        du0 += 9;
        du1 += 6;
      }
    }

    dt += 9;
    dtx += 9;
    dztx += 9;
    dn += 3;
    dn_xi += 3;
    dn_eta += 3;
  }
}

/*
  Add the product of the stress with the second derivative of the nonlinear
  strain. Note that this only fills half the matrix. The contribution is
  symmetric.

  input:
  N, Na, Nb: the shape functions, and the derivatives of the shape
  functions along the parametric coordinates

  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the Ux with respect to the parametric coordinates

  C, Ct, Ctt, Cttt: the rate matrix and its first, second and third
  (yes, third!) derivative with respect to the rotation angles

  t: the transformation to the local coordinate system
  tx: the transformation t times the local Jacobian
  ztx: the derivative of the transformaion with respect to the
  through-thickness direction

  n, n_xi, n_eta: the normal direction and its derivative along the
  parametric directions

  num_points: the number of points (nodes) in the interpolant
  scale: a scaling factor (for instance the local determinant of the
  Jacobian times the Gauss quadrature weight)

  output:
  matrix: the matrix of values
*/
void add_large_rot_stress_bmat(
    TacsScalar matrix[], const int num_points, const TacsScalar scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar C[], const TacsScalar Ct[], const TacsScalar Ctt[],
    const TacsScalar Cttt[], const TacsScalar t[], const TacsScalar tx[],
    const TacsScalar ztx[], const TacsScalar n[], const TacsScalar n_xi[],
    const TacsScalar n_eta[]) {
  if (num_points > shellutils::MAX_NUM_NODES) {
    return;
  }

  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Extract the second derivative information
  const TacsScalar* C11 = &Ctt[0];
  const TacsScalar* C22 = &Ctt[9];
  const TacsScalar* C33 = &Ctt[18];
  const TacsScalar* C12 = &Ctt[27];
  const TacsScalar* C13 = &Ctt[36];
  const TacsScalar* C23 = &Ctt[45];

  // Extract the thrid derivative information
  const TacsScalar* C112 = &Cttt[0];
  const TacsScalar* C113 = &Cttt[9];
  const TacsScalar* C122 = &Cttt[18];
  const TacsScalar* C123 = &Cttt[27];
  const TacsScalar* C133 = &Cttt[36];
  const TacsScalar* C223 = &Cttt[45];
  const TacsScalar* C233 = &Cttt[54];

  // Transform the displacement gradient to the local frame
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);
  transform_vector3d(r, t);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);
  transform_vector3d(Cn_xi, t);
  transform_vector3d(Cn_eta, t);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);
  transform_vector3d(C1n, t);
  transform_vector3d(C2n, t);
  transform_vector3d(C3n, t);

  // Compute the derivatives of C times n,xi and n,eta
  TacsScalar C1n_xi[3], C2n_xi[3], C3n_xi[3];
  TacsScalar C1n_eta[3], C2n_eta[3], C3n_eta[3];
  transform_vector3d(C1n_xi, n_xi, C1);
  transform_vector3d(C2n_xi, n_xi, C2);
  transform_vector3d(C3n_xi, n_xi, C3);
  transform_vector3d(C1n_eta, n_eta, C1);
  transform_vector3d(C2n_eta, n_eta, C2);
  transform_vector3d(C3n_eta, n_eta, C3);

  // Transform the vectors to the local coordinate system
  transform_vector3d(C1n_xi, t);
  transform_vector3d(C2n_xi, t);
  transform_vector3d(C3n_xi, t);
  transform_vector3d(C1n_eta, t);
  transform_vector3d(C2n_eta, t);
  transform_vector3d(C3n_eta, t);

  // Compute the product of the second derivatives of C with n
  TacsScalar C11n[3], C22n[3], C33n[3], C12n[3], C13n[3], C23n[3];
  transform_vector3d(C11n, n, C11);
  transform_vector3d(C22n, n, C22);
  transform_vector3d(C33n, n, C33);
  transform_vector3d(C12n, n, C12);
  transform_vector3d(C13n, n, C13);
  transform_vector3d(C23n, n, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n, t);
  transform_vector3d(C22n, t);
  transform_vector3d(C33n, t);
  transform_vector3d(C12n, t);
  transform_vector3d(C13n, t);
  transform_vector3d(C23n, t);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

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

  // Store the derivatives of the displacements w.r.t. each node
  // This makes the second derivative computation faster
  TacsScalar dU0[9 * 6 * shellutils::MAX_NUM_NODES];
  TacsScalar dU1[6 * 6 * shellutils::MAX_NUM_NODES];

  TacsScalar* du0 = dU0;
  TacsScalar* du1 = dU1;

  // For each point
  for (int i = 0; i < num_points; i++) {
    // For the displacement components
    for (int ii = 0; ii < 3; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      // The mid-surface value of the displacement derivatives
      du0[0] = Ud[0] * tx[0] + Ud[1] * tx[1];
      du0[1] = Ud[0] * tx[3] + Ud[1] * tx[4];
      du0[2] = Ud[0] * tx[6] + Ud[1] * tx[7];

      du0[3] = Ud[2] * tx[0] + Ud[3] * tx[1];
      du0[4] = Ud[2] * tx[3] + Ud[3] * tx[4];
      du0[5] = Ud[2] * tx[6] + Ud[3] * tx[7];

      du0[6] = Ud[4] * tx[0] + Ud[5] * tx[1];
      du0[7] = Ud[4] * tx[3] + Ud[5] * tx[4];
      du0[8] = Ud[4] * tx[6] + Ud[5] * tx[7];
      du0 += 9;

      // The first-derivative values of the displacements
      du1[0] = Ud[0] * ztx[0] + Ud[1] * ztx[1];
      du1[1] = Ud[0] * ztx[3] + Ud[1] * ztx[4];

      du1[2] = Ud[2] * ztx[0] + Ud[3] * ztx[1];
      du1[3] = Ud[2] * ztx[3] + Ud[3] * ztx[4];

      du1[4] = Ud[4] * ztx[0] + Ud[5] * ztx[1];
      du1[5] = Ud[4] * ztx[3] + Ud[5] * ztx[4];
      du1 += 6;
    }

    // For each displacement component
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      if (ii == 3) {
        dr[0] = C1n[0] * N[i];
        dr[1] = C1n[1] * N[i];
        dr[2] = C1n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                  C1n_xi, C1n_eta, C11n, C12n, C13n);
      } else if (ii == 4) {
        dr[0] = C2n[0] * N[i];
        dr[1] = C2n[1] * N[i];
        dr[2] = C2n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                  C2n_xi, C2n_eta, C12n, C22n, C23n);
      } else {  // ii == 5
        dr[0] = C3n[0] * N[i];
        dr[1] = C3n[1] * N[i];
        dr[2] = C3n[2] * N[i];
        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                  C3n_xi, C3n_eta, C13n, C23n, C33n);
      }

      // The mid-surface value of the displacement derivatives
      du0[0] = dr[0] * tx[2];
      du0[1] = dr[0] * tx[5];
      du0[2] = dr[0] * tx[8];

      du0[3] = dr[1] * tx[2];
      du0[4] = dr[1] * tx[5];
      du0[5] = dr[1] * tx[8];

      du0[6] = dr[2] * tx[2];
      du0[7] = dr[2] * tx[5];
      du0[8] = dr[2] * tx[8];
      du0 += 9;

      // The first-derivative values of the displacements
      du1[0] = dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dr[0] * ztx[2];
      du1[1] = dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dr[0] * ztx[5];

      du1[2] = dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dr[1] * ztx[2];
      du1[3] = dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dr[1] * ztx[5];

      du1[4] = dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dr[2] * ztx[2];
      du1[5] = dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dr[2] * ztx[5];
      du1 += 6;
    }
  }

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Set the i-counter from the beginning of the row
      TacsScalar* di0 = &dU0[9 * (6 * i + ii)];
      TacsScalar* di1 = &dU1[6 * (6 * i + ii)];

      // Start the j-counter from the i-counter
      TacsScalar* dj0 = dU0;
      TacsScalar* dj1 = dU1;

      // For each point
      for (int j = 0; j <= i; j++) {
        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        // For each displacement component
        for (int jj = 0; jj < end; jj++) {
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

  // Compute the product of the second derivative of C with n_xi
  TacsScalar C11n_xi[3], C22n_xi[3], C33n_xi[3];
  TacsScalar C12n_xi[3], C13n_xi[3], C23n_xi[3];
  transform_vector3d(C11n_xi, n_xi, C11);
  transform_vector3d(C22n_xi, n_xi, C22);
  transform_vector3d(C33n_xi, n_xi, C33);
  transform_vector3d(C12n_xi, n_xi, C12);
  transform_vector3d(C13n_xi, n_xi, C13);
  transform_vector3d(C23n_xi, n_xi, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n_xi, t);
  transform_vector3d(C22n_xi, t);
  transform_vector3d(C33n_xi, t);
  transform_vector3d(C12n_xi, t);
  transform_vector3d(C13n_xi, t);
  transform_vector3d(C23n_xi, t);

  // Compute the product of the second derivative of C with n_eta
  TacsScalar C11n_eta[3], C22n_eta[3], C33n_eta[3];
  TacsScalar C12n_eta[3], C13n_eta[3], C23n_eta[3];
  transform_vector3d(C11n_eta, n_eta, C11);
  transform_vector3d(C22n_eta, n_eta, C22);
  transform_vector3d(C33n_eta, n_eta, C33);
  transform_vector3d(C12n_eta, n_eta, C12);
  transform_vector3d(C13n_eta, n_eta, C13);
  transform_vector3d(C23n_eta, n_eta, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n_eta, t);
  transform_vector3d(C22n_eta, t);
  transform_vector3d(C33n_eta, t);
  transform_vector3d(C12n_eta, t);
  transform_vector3d(C13n_eta, t);
  transform_vector3d(C23n_eta, t);

  // Compute the inner product of the third derivatives of C with n
  TacsScalar C112n[3], C113n[3], C122n[3], C123n[3], C133n[3], C223n[3],
      C233n[3];
  transform_vector3d(C112n, n, C112);
  transform_vector3d(C113n, n, C113);
  transform_vector3d(C122n, n, C122);
  transform_vector3d(C123n, n, C123);
  transform_vector3d(C133n, n, C133);
  transform_vector3d(C223n, n, C223);
  transform_vector3d(C233n, n, C233);

  // Transform the vectors to the local system
  transform_vector3d(C112n, t);
  transform_vector3d(C113n, t);
  transform_vector3d(C122n, t);
  transform_vector3d(C123n, t);
  transform_vector3d(C133n, t);
  transform_vector3d(C223n, t);
  transform_vector3d(C233n, t);

  // Compute the remaining terms required for the second derivatives
  C1n[0] = -C1n[0];
  C1n[1] = -C1n[1];
  C1n[2] = -C1n[2];

  C2n[0] = -C2n[0];
  C2n[1] = -C2n[1];
  C2n[2] = -C2n[2];

  C3n[0] = -C3n[0];
  C3n[1] = -C3n[1];
  C3n[2] = -C3n[2];

  // Now add the terms that represent the second derivative of the
  // displacement gradient
  for (int i = 0; i < num_points; i++) {
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // For each point
      for (int j = 0; j <= i; j++) {
        // Skip the displacements - only the rotations are involved in
        // this calculation
        mat += 3;

        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        for (int jj = 3; jj < end; jj++) {
          TacsScalar ddr[3], ddr_xi[3], ddr_eta[3];

          if (ii == 3) {
            if (jj == 3) {
              ddr[0] = C11n[0] * N[i] * N[j];
              ddr[1] = C11n[1] * N[i] * N[j];
              ddr[2] = C11n[2] * N[i] * N[j];
              // Compute the derivatives 111, 112, 113
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C11n, C11n_xi, C11n_eta, C1n, C112n, C113n);
            } else if (jj == 4) {
              ddr[0] = C12n[0] * N[i] * N[j];
              ddr[1] = C12n[1] * N[i] * N[j];
              ddr[2] = C12n[2] * N[i] * N[j];
              // Compute the derivatives 112, 122, 123
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C12n, C12n_xi, C12n_eta, C112n, C122n, C123n);
            } else {  // jj == 5
              ddr[0] = C13n[0] * N[i] * N[j];
              ddr[1] = C13n[1] * N[i] * N[j];
              ddr[2] = C13n[2] * N[i] * N[j];
              // Compute the derivatives 113, 123, 133
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C13n, C13n_xi, C13n_eta, C113n, C123n, C133n);
            }
          } else if (ii == 4) {
            if (jj == 3) {
              ddr[0] = C12n[0] * N[i] * N[j];
              ddr[1] = C12n[1] * N[i] * N[j];
              ddr[2] = C12n[2] * N[i] * N[j];
              // Compute the derivatives 112, 122, 123
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C12n, C12n_xi, C12n_eta, C112n, C122n, C123n);
            } else if (jj == 4) {
              ddr[0] = C22n[0] * N[i] * N[j];
              ddr[1] = C22n[1] * N[i] * N[j];
              ddr[2] = C22n[2] * N[i] * N[j];
              // Compute the derivatives 122, 222, 223
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C22n, C22n_xi, C22n_eta, C122n, C2n, C223n);
            } else {  // jj == 5
              ddr[0] = C23n[0] * N[i] * N[j];
              ddr[1] = C23n[1] * N[i] * N[j];
              ddr[2] = C23n[2] * N[i] * N[j];
              // Compute the derivatives 123, 223, 233
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C23n, C23n_xi, C23n_eta, C123n, C223n, C233n);
            }
          } else {
            if (jj == 3) {
              ddr[0] = C13n[0] * N[i] * N[j];
              ddr[1] = C13n[1] * N[i] * N[j];
              ddr[2] = C13n[2] * N[i] * N[j];
              // Compute the derivatives 113, 123, 133
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C13n, C13n_xi, C13n_eta, C113n, C123n, C133n);
            } else if (jj == 4) {
              ddr[0] = C23n[0] * N[i] * N[j];
              ddr[1] = C23n[1] * N[i] * N[j];
              ddr[2] = C23n[2] * N[i] * N[j];
              // Compute the derivatives 23,  123, 223, 233
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C23n, C23n_xi, C23n_eta, C123n, C223n, C233n);
            } else {  // jj == 5
              ddr[0] = C33n[0] * N[i] * N[j];
              ddr[1] = C33n[1] * N[i] * N[j];
              ddr[2] = C33n[2] * N[i] * N[j];
              // Compute the derivatives 33, 133, 233, 333
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C33n, C33n_xi, C33n_eta, C133n, C233n, C3n);
            }
          }

          // The mid-surface value of the displacement derivatives
          TacsScalar ddux = ddr[0] * tx[2];
          TacsScalar dduy = ddr[0] * tx[5];
          TacsScalar dduz = ddr[0] * tx[8];

          TacsScalar ddvx = ddr[1] * tx[2];
          TacsScalar ddvy = ddr[1] * tx[5];
          TacsScalar ddvz = ddr[1] * tx[8];

          TacsScalar ddwx = ddr[2] * tx[2];
          TacsScalar ddwy = ddr[2] * tx[5];
          TacsScalar ddwz = ddr[2] * tx[8];

          // The first-derivative values of the displacements
          TacsScalar ddux1 =
              (ddr_xi[0] * tx[0] + ddr_eta[0] * tx[1] + ddr[0] * ztx[2]);
          TacsScalar dduy1 =
              (ddr_xi[0] * tx[3] + ddr_eta[0] * tx[4] + ddr[0] * ztx[5]);

          TacsScalar ddvx1 =
              (ddr_xi[1] * tx[0] + ddr_eta[1] * tx[1] + ddr[1] * ztx[2]);
          TacsScalar ddvy1 =
              (ddr_xi[1] * tx[3] + ddr_eta[1] * tx[4] + ddr[1] * ztx[5]);

          TacsScalar ddwx1 =
              (ddr_xi[2] * tx[0] + ddr_eta[2] * tx[1] + ddr[2] * ztx[2]);
          TacsScalar ddwy1 =
              (ddr_xi[2] * tx[3] + ddr_eta[2] * tx[4] + ddr[2] * ztx[5]);

          TacsScalar B[8];
          // Compute the in-plane components
          B[0] = ddux + (ux * ddux + vx * ddvx + wx * ddwx);
          B[1] = ddvy + (uy * dduy + vy * ddvy + wy * ddwy);
          B[2] = dduy + ddvx +
                 (ddux * uy + ddvx * vy + ddwx * wy + ux * dduy + vx * ddvy +
                  wx * ddwy);

          // Compute the bending components of the strain
          B[3] = ddux1 + (ddux * ux1 + ddvx * vx1 + ddwx * wx1 + ux * ddux1 +
                          vx * ddvx1 + wx * ddwx1);
          B[4] = ddvy1 + (dduy * uy1 + ddvy * vy1 + ddwy * wy1 + uy * dduy1 +
                          vy * ddvy1 + wy * ddwy1);
          B[5] = dduy1 + ddvx1 +
                 (ddux1 * uy + ddvx1 * vy + ddwx1 * wy + ddux * uy1 +
                  ddvx * vy1 + ddwx * wy1 + ux1 * dduy + vx1 * ddvy +
                  wx1 * ddwy + ux * dduy1 + vx * ddvy1 + wx * ddwy1);

          // Compute the shear terms
          B[6] = ddvz + ddwy +
                 (dduz * uy + ddvz * vy + ddwz * wy + uz * dduy + vz * ddvy +
                  wz * ddwy);
          B[7] = dduz + ddwx +
                 (dduz * ux + ddvz * vx + ddwz * wx + uz * ddux + vz * ddvx +
                  wz * ddwx);

          mat[0] +=
              scale * (B[0] * stress[0] + B[1] * stress[1] + B[2] * stress[2] +
                       B[3] * stress[3] + B[4] * stress[4] + B[5] * stress[5] +
                       B[6] * stress[6] + B[7] * stress[7]);
          mat += 1;
        }
      }
    }
  }
}

/*
  Compute the nonlinear bending components of the strain using the
  large angle shell formulation.

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
*/
void large_rot_bend_strain(TacsScalar strain[], const TacsScalar Ux[],
                           const TacsScalar Uxd[], const TacsScalar C[],
                           const TacsScalar Ct[], const TacsScalar t[],
                           const TacsScalar tx[], const TacsScalar ztx[],
                           const TacsScalar n[], const TacsScalar n_xi[],
                           const TacsScalar n_eta[]) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Transform the displacement into the local shell coordinates
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
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

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);
}

/*
  Compute the derivative of the nonlinear components of the strain
  using the large angle shell formulation.

  This function uses the same expression for the displacement
  gradient:

  U_{e,e} = t * U_{x,xi} * tx^{T}

  These are then used to construct the nonlinear expressions for the
  strain in the local Cartesian coordinate frame e.

  input:
  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the global displacements along the coordinate
  lines of the local shell
  C, Ct: the normal rotation rate matrices

  t, dt: the transformation and the derivative
  tx, dtx: the transformation times the Jacobian and its derivative
  ztx, dztx: the derivative of the transformation times the Jacobian w.r.t. the
  through-thickness direction and its derivative

  n, dn: the shell normal and the derivative of the shell normal
  n_xi, dn_xi: the derivative of the shell normal w.r.t. the coordinate line xi
  and its derivative
  n_eta, dn_eta: the derivative of the shell normal w.r.t. the coordinate line
  eta and its derivative

  num_componnets: the number of components to differentiate

  output:
  strain: the value of the strain at the current point
  dstrain: the value of the derivative of the strain at the current point
*/
void large_rot_bend_strain_sens(
    TacsScalar strain[], TacsScalar dstrain[], const TacsScalar Ux[],
    const TacsScalar Uxd[], const TacsScalar C[], const TacsScalar Ct[],
    const TacsScalar t[], const TacsScalar dt[], const TacsScalar tx[],
    const TacsScalar dtx[], const TacsScalar ztx[], const TacsScalar dztx[],
    const TacsScalar n[], const TacsScalar dn[], const TacsScalar n_xi[],
    const TacsScalar dn_xi[], const TacsScalar n_eta[],
    const TacsScalar dn_eta[], const int num_components) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Transform the displacement into the local shell coordinates
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
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

  // Compute the bending components of the strain
  strain[3] = ux1 + ux * ux1 + vx * vx1 + wx * wx1;
  strain[4] = vy1 + uy * uy1 + vy * vy1 + wy * wy1;
  strain[5] = uy1 + vx1 +
              (ux1 * uy + vx1 * vy + wx1 * wy + ux * uy1 + vx * vy1 + wx * wy1);

  // Now go through for each component of the derivative and compute the
  // derivative of the strain w.r.t. the input
  for (int i = 0; i < num_components; i++) {
    // Transform the displacement into the local shell coordinates
    TacsScalar dUd[6];
    transform_displ_gradient(dUd, dt, Uxd);

    // Compute the rate of change of the displacement through
    // the thickness
    TacsScalar dr[3];  // Compute dr = C*dn
    transform_vector3d(dr, dn, C);

    // Compute the product of the derivatives of n with C
    TacsScalar dCn_xi[3], dCn_eta[3];
    transform_vector3d(dCn_xi, dn_xi, C);
    transform_vector3d(dCn_eta, dn_eta, C);

    // Compute the product of the derivatives of C with n
    TacsScalar dC1n[3], dC2n[3], dC3n[3];
    transform_vector3d(dC1n, dn, C1);
    transform_vector3d(dC2n, dn, C2);
    transform_vector3d(dC3n, dn, C3);

    // Recompute r, r_xi and r_eta before transformation
    // This is required for the derivative operations
    transform_vector3d(r, n, C);
    compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

    // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
    TacsScalar dr_xi[3], dr_eta[3];
    compute_normal_rate(dr_xi, dr_eta, Uxd, dCn_xi, dCn_eta, dC1n, dC2n, dC3n);

    // Transform the displacement rate through the thickness into the
    // local shell coordinates
    transform_vector3d_sens(dr, r, t, dt);
    transform_vector3d_sens(dr_xi, r_xi, t, dt);
    transform_vector3d_sens(dr_eta, r_eta, t, dt);

    transform_vector3d(r, t);
    transform_vector3d(r_xi, t);
    transform_vector3d(r_eta, t);

    // The mid-surface value of the displacement derivatives
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

    dt += 9;
    dtx += 9;
    dztx += 9;
    dn += 3;
    dn_xi += 3;
    dn_eta += 3;
  }
}

/*
  Compute the derivative of the nonlinear strain expressions w.r.t
  the nodal displacements/rotations. Include only the contributions
  to the bending strain.

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
*/
void large_rot_bend_bmat(TacsScalar B[], const int num_points, const double N[],
                         const double Na[], const double Nb[],
                         const TacsScalar Ux[], const TacsScalar Uxd[],
                         const TacsScalar C[], const TacsScalar Ct[],
                         const TacsScalar Ctt[], const TacsScalar t[],
                         const TacsScalar tx[], const TacsScalar ztx[],
                         const TacsScalar n[], const TacsScalar n_xi[],
                         const TacsScalar n_eta[]) {
  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Extract the second derivative information
  const TacsScalar* C11 = &Ctt[0];
  const TacsScalar* C22 = &Ctt[9];
  const TacsScalar* C33 = &Ctt[18];
  const TacsScalar* C12 = &Ctt[27];
  const TacsScalar* C13 = &Ctt[36];
  const TacsScalar* C23 = &Ctt[45];

  // Transform the displacement gradient to the local frame
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);

  // Transform the displacement rate through the thickness into the
  // local shell coordinates
  transform_vector3d(r, t);
  transform_vector3d(Cn_xi, t);
  transform_vector3d(Cn_eta, t);
  transform_vector3d(C1n, t);
  transform_vector3d(C2n, t);
  transform_vector3d(C3n, t);

  // Compute the derivatives of C times n,xi and n,eta
  TacsScalar C1n_xi[3], C2n_xi[3], C3n_xi[3];
  TacsScalar C1n_eta[3], C2n_eta[3], C3n_eta[3];
  transform_vector3d(C1n_xi, n_xi, C1);
  transform_vector3d(C2n_xi, n_xi, C2);
  transform_vector3d(C3n_xi, n_xi, C3);
  transform_vector3d(C1n_eta, n_eta, C1);
  transform_vector3d(C2n_eta, n_eta, C2);
  transform_vector3d(C3n_eta, n_eta, C3);

  // Transform the vectors to the local coordinate system
  transform_vector3d(C1n_xi, t);
  transform_vector3d(C2n_xi, t);
  transform_vector3d(C3n_xi, t);
  transform_vector3d(C1n_eta, t);
  transform_vector3d(C2n_eta, t);
  transform_vector3d(C3n_eta, t);

  // Compute the product of the second derivatives of C with n
  TacsScalar C11n[3], C22n[3], C33n[3], C12n[3], C13n[3], C23n[3];
  transform_vector3d(C11n, n, C11);
  transform_vector3d(C22n, n, C22);
  transform_vector3d(C33n, n, C33);
  transform_vector3d(C12n, n, C12);
  transform_vector3d(C13n, n, C13);
  transform_vector3d(C23n, n, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n, t);
  transform_vector3d(C22n, t);
  transform_vector3d(C33n, t);
  transform_vector3d(C12n, t);
  transform_vector3d(C13n, t);
  transform_vector3d(C23n, t);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

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
    // For the displacement components
    for (int ii = 0; ii < 3; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = Ud[0] * tx[0] + Ud[1] * tx[1];
      TacsScalar duy = Ud[0] * tx[3] + Ud[1] * tx[4];

      TacsScalar dvx = Ud[2] * tx[0] + Ud[3] * tx[1];
      TacsScalar dvy = Ud[2] * tx[3] + Ud[3] * tx[4];

      TacsScalar dwx = Ud[4] * tx[0] + Ud[5] * tx[1];
      TacsScalar dwy = Ud[4] * tx[3] + Ud[5] * tx[4];

      // The first-derivative values of the displacements
      TacsScalar dux1 = Ud[0] * ztx[0] + Ud[1] * ztx[1];
      TacsScalar duy1 = Ud[0] * ztx[3] + Ud[1] * ztx[4];

      TacsScalar dvx1 = Ud[2] * ztx[0] + Ud[3] * ztx[1];
      TacsScalar dvy1 = Ud[2] * ztx[3] + Ud[3] * ztx[4];

      TacsScalar dwx1 = Ud[4] * ztx[0] + Ud[5] * ztx[1];
      TacsScalar dwy1 = Ud[4] * ztx[3] + Ud[5] * ztx[4];

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
    }

    // For each displacement component
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      if (ii == 3) {
        dr[0] = C1n[0] * N[i];
        dr[1] = C1n[1] * N[i];
        dr[2] = C1n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                  C1n_xi, C1n_eta, C11n, C12n, C13n);
      } else if (ii == 4) {
        dr[0] = C2n[0] * N[i];
        dr[1] = C2n[1] * N[i];
        dr[2] = C2n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                  C2n_xi, C2n_eta, C12n, C22n, C23n);
      } else {  // ii == 5
        dr[0] = C3n[0] * N[i];
        dr[1] = C3n[1] * N[i];
        dr[2] = C3n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                  C3n_xi, C3n_eta, C13n, C23n, C33n);
      }

      // The mid-surface value of the displacement derivatives
      TacsScalar dux = dr[0] * tx[2];
      TacsScalar duy = dr[0] * tx[5];

      TacsScalar dvx = dr[1] * tx[2];
      TacsScalar dvy = dr[1] * tx[5];

      TacsScalar dwx = dr[2] * tx[2];
      TacsScalar dwy = dr[2] * tx[5];

      // The first-derivative values of the displacements
      TacsScalar dux1 = dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dr[0] * ztx[2];
      TacsScalar duy1 = dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dr[0] * ztx[5];

      TacsScalar dvx1 = dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dr[1] * ztx[2];
      TacsScalar dvy1 = dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dr[1] * ztx[5];

      TacsScalar dwx1 = dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dr[2] * ztx[2];
      TacsScalar dwy1 = dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dr[2] * ztx[5];

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
    }
  }
}

/*
  Add the product of the stress with the second derivative of the
  nonlinear strain. Note that this only fills half the matrix. The
  contribution is symmetric. This function includes only those
  contributions from the bending strain.

  input:
  N, Na, Nb: the shape functions, and the derivatives of the shape
  functions along the parametric coordinates

  Ux: the values of the displacements in the global coordinate frame
  Uxd: the derivative of the Ux with respect to the parametric coordinates

  C, Ct, Ctt, Cttt: the rate matrix and its first, second and third
  (yes, third!) derivative with respect to the rotation angles

  t: the transformation to the local coordinate system
  tx: the transformation t times the local Jacobian
  ztx: the derivative of the transformaion with respect to the
  through-thickness direction

  n, n_xi, n_eta: the normal direction and its derivative along the
  parametric directions

  num_points: the number of points (nodes) in the interpolant
  scale: a scaling factor (for instance the local determinant of the
  Jacobian times the Gauss quadrature weight)

  output:
  matrix: the matrix of values
*/
void add_large_rot_bend_stress_bmat(
    TacsScalar matrix[], const int num_points, const TacsScalar scale,
    const TacsScalar stress[], const double N[], const double Na[],
    const double Nb[], const TacsScalar Ux[], const TacsScalar Uxd[],
    const TacsScalar C[], const TacsScalar Ct[], const TacsScalar Ctt[],
    const TacsScalar Cttt[], const TacsScalar t[], const TacsScalar tx[],
    const TacsScalar ztx[], const TacsScalar n[], const TacsScalar n_xi[],
    const TacsScalar n_eta[]) {
  if (num_points > shellutils::MAX_NUM_NODES) {
    return;
  }

  // Extract the derivative information
  const TacsScalar* C1 = &Ct[0];
  const TacsScalar* C2 = &Ct[9];
  const TacsScalar* C3 = &Ct[18];

  // Extract the second derivative information
  const TacsScalar* C11 = &Ctt[0];
  const TacsScalar* C22 = &Ctt[9];
  const TacsScalar* C33 = &Ctt[18];
  const TacsScalar* C12 = &Ctt[27];
  const TacsScalar* C13 = &Ctt[36];
  const TacsScalar* C23 = &Ctt[45];

  // Extract the thrid derivative information
  const TacsScalar* C112 = &Cttt[0];
  const TacsScalar* C113 = &Cttt[9];
  const TacsScalar* C122 = &Cttt[18];
  const TacsScalar* C123 = &Cttt[27];
  const TacsScalar* C133 = &Cttt[36];
  const TacsScalar* C223 = &Cttt[45];
  const TacsScalar* C233 = &Cttt[54];

  // Transform the displacement gradient to the local frame
  TacsScalar Ud[6];
  transform_displ_gradient(Ud, t, Uxd);

  // Compute the rate of change of the displacement through
  // the thickness
  TacsScalar r[3];  // Compute r = C*n
  transform_vector3d(r, n, C);
  transform_vector3d(r, t);

  // Compute the product of the derivatives of n with C
  TacsScalar Cn_xi[3], Cn_eta[3];
  transform_vector3d(Cn_xi, n_xi, C);
  transform_vector3d(Cn_eta, n_eta, C);
  transform_vector3d(Cn_xi, t);
  transform_vector3d(Cn_eta, t);

  // Compute the product of the derivatives of C with n
  TacsScalar C1n[3], C2n[3], C3n[3];
  transform_vector3d(C1n, n, C1);
  transform_vector3d(C2n, n, C2);
  transform_vector3d(C3n, n, C3);
  transform_vector3d(C1n, t);
  transform_vector3d(C2n, t);
  transform_vector3d(C3n, t);

  // Compute the derivatives of C times n,xi and n,eta
  TacsScalar C1n_xi[3], C2n_xi[3], C3n_xi[3];
  TacsScalar C1n_eta[3], C2n_eta[3], C3n_eta[3];
  transform_vector3d(C1n_xi, n_xi, C1);
  transform_vector3d(C2n_xi, n_xi, C2);
  transform_vector3d(C3n_xi, n_xi, C3);
  transform_vector3d(C1n_eta, n_eta, C1);
  transform_vector3d(C2n_eta, n_eta, C2);
  transform_vector3d(C3n_eta, n_eta, C3);

  // Transform the vectors to the local coordinate system
  transform_vector3d(C1n_xi, t);
  transform_vector3d(C2n_xi, t);
  transform_vector3d(C3n_xi, t);
  transform_vector3d(C1n_eta, t);
  transform_vector3d(C2n_eta, t);
  transform_vector3d(C3n_eta, t);

  // Compute the product of the second derivatives of C with n
  TacsScalar C11n[3], C22n[3], C33n[3], C12n[3], C13n[3], C23n[3];
  transform_vector3d(C11n, n, C11);
  transform_vector3d(C22n, n, C22);
  transform_vector3d(C33n, n, C33);
  transform_vector3d(C12n, n, C12);
  transform_vector3d(C13n, n, C13);
  transform_vector3d(C23n, n, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n, t);
  transform_vector3d(C22n, t);
  transform_vector3d(C33n, t);
  transform_vector3d(C12n, t);
  transform_vector3d(C13n, t);
  transform_vector3d(C23n, t);

  // Compute r,xi = C*n,xi + C,xi*n and r,eta = C*n,eta + C,eta*n
  TacsScalar r_xi[3], r_eta[3];
  compute_normal_rate(r_xi, r_eta, Uxd, Cn_xi, Cn_eta, C1n, C2n, C3n);

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

  // Store the derivatives of the displacements w.r.t. each node
  // This makes the second derivative computation faster
  TacsScalar dU0[6 * 6 * shellutils::MAX_NUM_NODES];
  TacsScalar dU1[6 * 6 * shellutils::MAX_NUM_NODES];

  TacsScalar* du0 = dU0;
  TacsScalar* du1 = dU1;

  // For each point
  for (int i = 0; i < num_points; i++) {
    // For the displacement components
    for (int ii = 0; ii < 3; ii++) {
      transform_displ_gradient_bmat(Ud, ii, t, Na[i], Nb[i]);

      // The mid-surface value of the displacement derivatives
      du0[0] = Ud[0] * tx[0] + Ud[1] * tx[1];
      du0[1] = Ud[0] * tx[3] + Ud[1] * tx[4];

      du0[2] = Ud[2] * tx[0] + Ud[3] * tx[1];
      du0[3] = Ud[2] * tx[3] + Ud[3] * tx[4];

      du0[4] = Ud[4] * tx[0] + Ud[5] * tx[1];
      du0[5] = Ud[4] * tx[3] + Ud[5] * tx[4];
      du0 += 6;

      // The first-derivative values of the displacements
      du1[0] = Ud[0] * ztx[0] + Ud[1] * ztx[1];
      du1[1] = Ud[0] * ztx[3] + Ud[1] * ztx[4];

      du1[2] = Ud[2] * ztx[0] + Ud[3] * ztx[1];
      du1[3] = Ud[2] * ztx[3] + Ud[3] * ztx[4];

      du1[4] = Ud[4] * ztx[0] + Ud[5] * ztx[1];
      du1[5] = Ud[4] * ztx[3] + Ud[5] * ztx[4];
      du1 += 6;
    }

    // For each displacement component
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar dr[3], dr_xi[3], dr_eta[3];
      if (ii == 3) {
        dr[0] = C1n[0] * N[i];
        dr[1] = C1n[1] * N[i];
        dr[2] = C1n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C1n,
                                  C1n_xi, C1n_eta, C11n, C12n, C13n);
      } else if (ii == 4) {
        dr[0] = C2n[0] * N[i];
        dr[1] = C2n[1] * N[i];
        dr[2] = C2n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C2n,
                                  C2n_xi, C2n_eta, C12n, C22n, C23n);
      } else {  // ii == 5
        dr[0] = C3n[0] * N[i];
        dr[1] = C3n[1] * N[i];
        dr[2] = C3n[2] * N[i];

        compute_normal_rate_theta(dr_xi, dr_eta, Uxd, N[i], Na[i], Nb[i], C3n,
                                  C3n_xi, C3n_eta, C13n, C23n, C33n);
      }

      // The mid-surface value of the displacement derivatives
      du0[0] = dr[0] * tx[2];
      du0[1] = dr[0] * tx[5];

      du0[2] = dr[1] * tx[2];
      du0[3] = dr[1] * tx[5];

      du0[4] = dr[2] * tx[2];
      du0[5] = dr[2] * tx[5];
      du0 += 6;

      // The first-derivative values of the displacements
      du1[0] = dr_xi[0] * tx[0] + dr_eta[0] * tx[1] + dr[0] * ztx[2];
      du1[1] = dr_xi[0] * tx[3] + dr_eta[0] * tx[4] + dr[0] * ztx[5];

      du1[2] = dr_xi[1] * tx[0] + dr_eta[1] * tx[1] + dr[1] * ztx[2];
      du1[3] = dr_xi[1] * tx[3] + dr_eta[1] * tx[4] + dr[1] * ztx[5];

      du1[4] = dr_xi[2] * tx[0] + dr_eta[2] * tx[1] + dr[2] * ztx[2];
      du1[5] = dr_xi[2] * tx[3] + dr_eta[2] * tx[4] + dr[2] * ztx[5];
      du1 += 6;
    }
  }

  for (int i = 0; i < num_points; i++) {
    for (int ii = 0; ii < 6; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // Set the i-counter from the beginning of the row
      TacsScalar* di0 = &dU0[6 * (6 * i + ii)];
      TacsScalar* di1 = &dU1[6 * (6 * i + ii)];

      // Start the j-counter from the i-counter
      TacsScalar* dj0 = dU0;
      TacsScalar* dj1 = dU1;

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

  // Compute the product of the second derivative of C with n_xi
  TacsScalar C11n_xi[3], C22n_xi[3], C33n_xi[3];
  TacsScalar C12n_xi[3], C13n_xi[3], C23n_xi[3];
  transform_vector3d(C11n_xi, n_xi, C11);
  transform_vector3d(C22n_xi, n_xi, C22);
  transform_vector3d(C33n_xi, n_xi, C33);
  transform_vector3d(C12n_xi, n_xi, C12);
  transform_vector3d(C13n_xi, n_xi, C13);
  transform_vector3d(C23n_xi, n_xi, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n_xi, t);
  transform_vector3d(C22n_xi, t);
  transform_vector3d(C33n_xi, t);
  transform_vector3d(C12n_xi, t);
  transform_vector3d(C13n_xi, t);
  transform_vector3d(C23n_xi, t);

  // Compute the product of the second derivative of C with n_eta
  TacsScalar C11n_eta[3], C22n_eta[3], C33n_eta[3];
  TacsScalar C12n_eta[3], C13n_eta[3], C23n_eta[3];
  transform_vector3d(C11n_eta, n_eta, C11);
  transform_vector3d(C22n_eta, n_eta, C22);
  transform_vector3d(C33n_eta, n_eta, C33);
  transform_vector3d(C12n_eta, n_eta, C12);
  transform_vector3d(C13n_eta, n_eta, C13);
  transform_vector3d(C23n_eta, n_eta, C23);

  // Transform the vectors to the local system
  transform_vector3d(C11n_eta, t);
  transform_vector3d(C22n_eta, t);
  transform_vector3d(C33n_eta, t);
  transform_vector3d(C12n_eta, t);
  transform_vector3d(C13n_eta, t);
  transform_vector3d(C23n_eta, t);

  // Compute the inner product of the third derivatives of C with n
  TacsScalar C112n[3], C113n[3], C122n[3], C123n[3], C133n[3], C223n[3],
      C233n[3];
  transform_vector3d(C112n, n, C112);
  transform_vector3d(C113n, n, C113);
  transform_vector3d(C122n, n, C122);
  transform_vector3d(C123n, n, C123);
  transform_vector3d(C133n, n, C133);
  transform_vector3d(C223n, n, C223);
  transform_vector3d(C233n, n, C233);

  // Transform the vectors to the local system
  transform_vector3d(C112n, t);
  transform_vector3d(C113n, t);
  transform_vector3d(C122n, t);
  transform_vector3d(C123n, t);
  transform_vector3d(C133n, t);
  transform_vector3d(C223n, t);
  transform_vector3d(C233n, t);

  // Compute the remaining terms required for the second derivatives
  C1n[0] = -C1n[0];
  C1n[1] = -C1n[1];
  C1n[2] = -C1n[2];

  C2n[0] = -C2n[0];
  C2n[1] = -C2n[1];
  C2n[2] = -C2n[2];

  C3n[0] = -C3n[0];
  C3n[1] = -C3n[1];
  C3n[2] = -C3n[2];

  // Now add the terms that represent the second derivative of the
  // displacement gradient
  for (int i = 0; i < num_points; i++) {
    for (int ii = 3; ii < 6; ii++) {
      TacsScalar* mat = &matrix[(6 * i + ii) * (6 * num_points)];

      // For each point
      for (int j = 0; j <= i; j++) {
        // Skip the displacements - only the rotations are involved in
        // this calculation
        mat += 3;

        int end = 6;
        if (i == j) {
          end = ii + 1;
        }

        for (int jj = 3; jj < end; jj++) {
          TacsScalar ddr[3], ddr_xi[3], ddr_eta[3];

          if (ii == 3) {
            if (jj == 3) {
              ddr[0] = C11n[0] * N[i] * N[j];
              ddr[1] = C11n[1] * N[i] * N[j];
              ddr[2] = C11n[2] * N[i] * N[j];
              // Compute the derivatives 111, 112, 113
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C11n, C11n_xi, C11n_eta, C1n, C112n, C113n);
            } else if (jj == 4) {
              ddr[0] = C12n[0] * N[i] * N[j];
              ddr[1] = C12n[1] * N[i] * N[j];
              ddr[2] = C12n[2] * N[i] * N[j];
              // Compute the derivatives 112, 122, 123
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C12n, C12n_xi, C12n_eta, C112n, C122n, C123n);
            } else {  // jj == 5
              ddr[0] = C13n[0] * N[i] * N[j];
              ddr[1] = C13n[1] * N[i] * N[j];
              ddr[2] = C13n[2] * N[i] * N[j];
              // Compute the derivatives 113, 123, 133
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C13n, C13n_xi, C13n_eta, C113n, C123n, C133n);
            }
          } else if (ii == 4) {
            if (jj == 3) {
              ddr[0] = C12n[0] * N[i] * N[j];
              ddr[1] = C12n[1] * N[i] * N[j];
              ddr[2] = C12n[2] * N[i] * N[j];
              // Compute the derivatives 112, 122, 123
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C12n, C12n_xi, C12n_eta, C112n, C122n, C123n);
            } else if (jj == 4) {
              ddr[0] = C22n[0] * N[i] * N[j];
              ddr[1] = C22n[1] * N[i] * N[j];
              ddr[2] = C22n[2] * N[i] * N[j];
              // Compute the derivatives 122, 222, 223
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C22n, C22n_xi, C22n_eta, C122n, C2n, C223n);
            } else {  // jj == 5
              ddr[0] = C23n[0] * N[i] * N[j];
              ddr[1] = C23n[1] * N[i] * N[j];
              ddr[2] = C23n[2] * N[i] * N[j];
              // Compute the derivatives 123, 223, 233
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C23n, C23n_xi, C23n_eta, C123n, C223n, C233n);
            }
          } else {
            if (jj == 3) {
              ddr[0] = C13n[0] * N[i] * N[j];
              ddr[1] = C13n[1] * N[i] * N[j];
              ddr[2] = C13n[2] * N[i] * N[j];
              // Compute the derivatives 113, 123, 133
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C13n, C13n_xi, C13n_eta, C113n, C123n, C133n);
            } else if (jj == 4) {
              ddr[0] = C23n[0] * N[i] * N[j];
              ddr[1] = C23n[1] * N[i] * N[j];
              ddr[2] = C23n[2] * N[i] * N[j];
              // Compute the derivatives 23,  123, 223, 233
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C23n, C23n_xi, C23n_eta, C123n, C223n, C233n);
            } else {  // jj == 5
              ddr[0] = C33n[0] * N[i] * N[j];
              ddr[1] = C33n[1] * N[i] * N[j];
              ddr[2] = C33n[2] * N[i] * N[j];
              // Compute the derivatives 33, 133, 233, 333
              compute_2nd_normal_rate_theta(
                  ddr_xi, ddr_eta, Uxd, N[i], Na[i], Nb[i], N[j], Na[j], Nb[j],
                  C33n, C33n_xi, C33n_eta, C133n, C233n, C3n);
            }
          }

          // The mid-surface value of the displacement derivatives
          TacsScalar ddux = ddr[0] * tx[2];
          TacsScalar dduy = ddr[0] * tx[5];

          TacsScalar ddvx = ddr[1] * tx[2];
          TacsScalar ddvy = ddr[1] * tx[5];

          TacsScalar ddwx = ddr[2] * tx[2];
          TacsScalar ddwy = ddr[2] * tx[5];

          // The first-derivative values of the displacements
          TacsScalar ddux1 =
              (ddr_xi[0] * tx[0] + ddr_eta[0] * tx[1] + ddr[0] * ztx[2]);
          TacsScalar dduy1 =
              (ddr_xi[0] * tx[3] + ddr_eta[0] * tx[4] + ddr[0] * ztx[5]);

          TacsScalar ddvx1 =
              (ddr_xi[1] * tx[0] + ddr_eta[1] * tx[1] + ddr[1] * ztx[2]);
          TacsScalar ddvy1 =
              (ddr_xi[1] * tx[3] + ddr_eta[1] * tx[4] + ddr[1] * ztx[5]);

          TacsScalar ddwx1 =
              (ddr_xi[2] * tx[0] + ddr_eta[2] * tx[1] + ddr[2] * ztx[2]);
          TacsScalar ddwy1 =
              (ddr_xi[2] * tx[3] + ddr_eta[2] * tx[4] + ddr[2] * ztx[5]);

          TacsScalar B[3];
          // Compute the bending components of the strain
          B[0] = ddux1 + (ddux * ux1 + ddvx * vx1 + ddwx * wx1 + ux * ddux1 +
                          vx * ddvx1 + wx * ddwx1);
          B[1] = ddvy1 + (dduy * uy1 + ddvy * vy1 + ddwy * wy1 + uy * dduy1 +
                          vy * ddvy1 + wy * ddwy1);
          B[2] = dduy1 + ddvx1 +
                 (ddux1 * uy + ddvx1 * vy + ddwx1 * wy + ddux * uy1 +
                  ddvx * vy1 + ddwx * wy1 + ux1 * dduy + vx1 * ddvy +
                  wx1 * ddwy + ux * dduy1 + vx * ddvy1 + wx * ddwy1);

          mat[0] +=
              scale * (B[0] * stress[3] + B[1] * stress[4] + B[2] * stress[5]);
          mat += 1;
        }
      }
    }
  }
}

TACS_END_NAMESPACE
