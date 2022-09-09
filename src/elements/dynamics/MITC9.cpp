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

#include "MITC9.h"

#include "TACSElementAlgebra.h"
#include "TACSElementQuaternion.h"
#include "TACSGaussQuadrature.h"

/*
  Rigid-body dynamics routines for TACS
*/

/*
  Write the relative error for components of a vector for a
  finite-difference or complex-step study to a given file

  input:
  fp:         the output file
  descript:   description of the component
  a:          analytic values
  b:          set of finite-difference/complex-step values
  size:       the number of values
  rel_err:    relative error tolerance
*/
static void writeErrorComponents(FILE *fp, const char *descript, TacsScalar *a,
                                 TacsScalar *fd, int size,
                                 double rel_err = 0.0) {
  int print_flag = 1;
  for (int i = 0; i < size; i++) {
    double rel = 0.0;
    if (a[i] != 0.0) {
      rel = fabs(TacsRealPart((a[i] - fd[i]) / a[i]));
    } else {
      rel = fabs(TacsRealPart((a[i] - fd[i])));
    }

    if (rel > rel_err || a[i] != a[i] || fd[i] != fd[i]) {
      if (print_flag) {
        fprintf(fp, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
                "Analytic", "Approximate", "Rel. Error");
        print_flag = 0;
      }
      if (a[i] != 0.0) {
        fprintf(fp, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, i,
                TacsRealPart(a[i]), TacsRealPart(fd[i]),
                fabs(TacsRealPart((a[i] - fd[i]) / a[i])));
      } else {
        fprintf(fp, "%s[%3d] %15.6e %15.6e\n", descript, i, TacsRealPart(a[i]),
                TacsRealPart(fd[i]));
      }
    }
  }
}

/*
  Evaluate the shape functions of the element given the u/v
  coordinates of the point

  input:
  u:   the first parametric coordinate
  v:   the second parametric coordinate

  output:
  N:   the shape functions
*/
static inline void computeShapeFunc(const double u, const double v,
                                    double N[]) {
  // Compute the shape functions
  double nu[3];
  nu[0] = -0.5 * u * (1.0 - u);
  nu[1] = (1.0 - u) * (1.0 + u);
  nu[2] = 0.5 * (1.0 + u) * u;

  double nv[3];
  nv[0] = -0.5 * v * (1.0 - v);
  nv[1] = (1.0 - v) * (1.0 + v);
  nv[2] = 0.5 * (1.0 + v) * v;

  // Compute the shape functions
  N[0] = nu[0] * nv[0];
  N[1] = nu[1] * nv[0];
  N[2] = nu[2] * nv[0];

  N[3] = nu[0] * nv[1];
  N[4] = nu[1] * nv[1];
  N[5] = nu[2] * nv[1];

  N[6] = nu[0] * nv[2];
  N[7] = nu[1] * nv[2];
  N[8] = nu[2] * nv[2];
}

/*
  Evaluate the derivatives of the shape functions of the element at
  given the u/v coordinates of the point

  input:
  u:    the first parametric coordinate
  v:    the second parametric coordinate

  output:
  Na:   the derivative of the shape functions w.r.t. u
  Nb:   the derivative of the shape functions w.r.t. v
*/
static inline void computeShapeFunc(const double u, const double v, double Na[],
                                    double Nb[]) {
  // Compute the shape functions and their derivatives
  double nu[3], dnu[3];
  nu[0] = -0.5 * u * (1.0 - u);
  nu[1] = (1.0 - u) * (1.0 + u);
  nu[2] = 0.5 * (1.0 + u) * u;

  dnu[0] = -0.5 + u;
  dnu[1] = -2.0 * u;
  dnu[2] = 0.5 + u;

  double nv[3], dnv[3];
  nv[0] = -0.5 * v * (1.0 - v);
  nv[1] = (1.0 - v) * (1.0 + v);
  nv[2] = 0.5 * (1.0 + v) * v;

  dnv[0] = -0.5 + v;
  dnv[1] = -2.0 * v;
  dnv[2] = 0.5 + v;

  // Compute the derivative of the shape funcs w.r.t. u
  Na[0] = dnu[0] * nv[0];
  Na[1] = dnu[1] * nv[0];
  Na[2] = dnu[2] * nv[0];

  Na[3] = dnu[0] * nv[1];
  Na[4] = dnu[1] * nv[1];
  Na[5] = dnu[2] * nv[1];

  Na[6] = dnu[0] * nv[2];
  Na[7] = dnu[1] * nv[2];
  Na[8] = dnu[2] * nv[2];

  // Compute the derivative of the shape funcs w.r.t. v
  Nb[0] = nu[0] * dnv[0];
  Nb[1] = nu[1] * dnv[0];
  Nb[2] = nu[2] * dnv[0];

  Nb[3] = nu[0] * dnv[1];
  Nb[4] = nu[1] * dnv[1];
  Nb[5] = nu[2] * dnv[1];

  Nb[6] = nu[0] * dnv[2];
  Nb[7] = nu[1] * dnv[2];
  Nb[8] = nu[2] * dnv[2];
}

/*
  Compute the inner product of a set of shape functions (or
  their derivative) with a 3*NUM_NODES array of nodes, variables
  etc.

  input:
  Na:   the shape function
  X:    the variables

  ouput:
  Xa:   the inner product = Na^{T}*X
*/
static inline void innerProduct(const double Na[], const TacsScalar X[],
                                TacsScalar Xa[]) {
  Xa[0] = (Na[0] * X[0] + Na[1] * X[3] + Na[2] * X[6] + Na[3] * X[9] +
           Na[4] * X[12] + Na[5] * X[15] + Na[6] * X[18] + Na[7] * X[21] +
           Na[8] * X[24]);
  Xa[1] = (Na[0] * X[1] + Na[1] * X[4] + Na[2] * X[7] + Na[3] * X[10] +
           Na[4] * X[13] + Na[5] * X[16] + Na[6] * X[19] + Na[7] * X[22] +
           Na[8] * X[25]);
  Xa[2] = (Na[0] * X[2] + Na[1] * X[5] + Na[2] * X[8] + Na[3] * X[11] +
           Na[4] * X[14] + Na[5] * X[17] + Na[6] * X[20] + Na[7] * X[23] +
           Na[8] * X[26]);
}

/*
  Compute the inner product of a set of shape functions (or
  their derivative) with a 8*NUM_NODES array of nodes, variables
  etc.

  input:
  Na:   the shape function
  X:    the variables

  ouput:
  Xa:   the inner product = Na^{T}*X
*/
static inline void innerProduct8(const double Na[], const TacsScalar X[],
                                 TacsScalar Xa[]) {
  Xa[0] = (Na[0] * X[0] + Na[1] * X[8] + Na[2] * X[16] + Na[3] * X[24] +
           Na[4] * X[32] + Na[5] * X[40] + Na[6] * X[48] + Na[7] * X[56] +
           Na[8] * X[64]);
  Xa[1] = (Na[0] * X[1] + Na[1] * X[9] + Na[2] * X[17] + Na[3] * X[25] +
           Na[4] * X[33] + Na[5] * X[41] + Na[6] * X[49] + Na[7] * X[57] +
           Na[8] * X[65]);
  Xa[2] = (Na[0] * X[2] + Na[1] * X[10] + Na[2] * X[18] + Na[3] * X[26] +
           Na[4] * X[34] + Na[5] * X[42] + Na[6] * X[50] + Na[7] * X[58] +
           Na[8] * X[66]);
}

/*
  Compute the normal implied by the interpolation of the normal at
  each node to an arbitrary point in the element. This will not be a
  normal vector, in general, but is required to obtain strain-free
  rigid-body motion.

  input:
  N:    the shape functions for the element
  Xr:   the frames at each node

  output:
  fnormal:   the "frame normal" obtained by interpolating the
  .          normal at each node
*/
static inline void computeFrameNormal(const double N[], const TacsScalar Xr[],
                                      TacsScalar fnormal[]) {
  fnormal[0] = (N[0] * Xr[2] + N[1] * Xr[11] + N[2] * Xr[20] + N[3] * Xr[29] +
                N[4] * Xr[38] + N[5] * Xr[47] + N[6] * Xr[56] + N[7] * Xr[65] +
                N[8] * Xr[74]);
  fnormal[1] = (N[0] * Xr[5] + N[1] * Xr[14] + N[2] * Xr[23] + N[3] * Xr[32] +
                N[4] * Xr[41] + N[5] * Xr[50] + N[6] * Xr[59] + N[7] * Xr[68] +
                N[8] * Xr[77]);
  fnormal[2] = (N[0] * Xr[8] + N[1] * Xr[17] + N[2] * Xr[26] + N[3] * Xr[35] +
                N[4] * Xr[44] + N[5] * Xr[53] + N[6] * Xr[62] + N[7] * Xr[71] +
                N[8] * Xr[80]);
}

/*
  Compute the derivative of the frame normal

*/
static inline void addFrameNormalSens(const TacsScalar fnd[], const double N[],
                                      TacsScalar Xrd[]) {
  Xrd[2] += N[0] * fnd[0];
  Xrd[5] += N[0] * fnd[1];
  Xrd[8] += N[0] * fnd[2];
  Xrd[11] += N[1] * fnd[0];
  Xrd[14] += N[1] * fnd[1];
  Xrd[17] += N[1] * fnd[2];
  Xrd[20] += N[2] * fnd[0];
  Xrd[23] += N[2] * fnd[1];
  Xrd[26] += N[2] * fnd[2];
  Xrd[29] += N[3] * fnd[0];
  Xrd[32] += N[3] * fnd[1];
  Xrd[35] += N[3] * fnd[2];
  Xrd[38] += N[4] * fnd[0];
  Xrd[41] += N[4] * fnd[1];
  Xrd[44] += N[4] * fnd[2];
  Xrd[47] += N[5] * fnd[0];
  Xrd[50] += N[5] * fnd[1];
  Xrd[53] += N[5] * fnd[2];
  Xrd[56] += N[6] * fnd[0];
  Xrd[59] += N[6] * fnd[1];
  Xrd[62] += N[6] * fnd[2];
  Xrd[65] += N[7] * fnd[0];
  Xrd[68] += N[7] * fnd[1];
  Xrd[71] += N[7] * fnd[2];
  Xrd[74] += N[8] * fnd[0];
  Xrd[77] += N[8] * fnd[1];
  Xrd[80] += N[8] * fnd[2];
}

/*
  Compute the derivative of the normal direction w.r.t. the parametric
  directions and multiply them

  This code computes the (approximate) through-thickness derivative of
  the inverse of the Jacobian matrix given as follows:

  d([X,r]^{-1})/dz

  The approximation enters through the fact that the normal direction
  is interpolated using the shape functions. This approach provides
  an approximation, but avoids the use of second derivatives.

  input:
  Na:   the derivative of the shape functions
  Nb:   the derivative of the shape functions
*/
static inline void computeNormalRateMat(const double Na[], const double Nb[],
                                        const TacsScalar Xr[],
                                        const TacsScalar Xdinv[],
                                        TacsScalar zXdinv[]) {
  // Compute the derivatives of the normal direction along the
  // parametric directions
  TacsScalar dn[9];
  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = -(Na[0] * Xr[2] + Na[1] * Xr[11] + Na[2] * Xr[20] + Na[3] * Xr[29] +
            Na[4] * Xr[38] + Na[5] * Xr[47] + Na[6] * Xr[56] + Na[7] * Xr[65] +
            Na[8] * Xr[74]);
  dn[3] = -(Na[0] * Xr[5] + Na[1] * Xr[14] + Na[2] * Xr[23] + Na[3] * Xr[32] +
            Na[4] * Xr[41] + Na[5] * Xr[50] + Na[6] * Xr[59] + Na[7] * Xr[68] +
            Na[8] * Xr[77]);
  dn[6] = -(Na[0] * Xr[8] + Na[1] * Xr[17] + Na[2] * Xr[26] + Na[3] * Xr[35] +
            Na[4] * Xr[44] + Na[5] * Xr[53] + Na[6] * Xr[62] + Na[7] * Xr[71] +
            Na[8] * Xr[80]);

  dn[1] = -(Nb[0] * Xr[2] + Nb[1] * Xr[11] + Nb[2] * Xr[20] + Nb[3] * Xr[29] +
            Nb[4] * Xr[38] + Nb[5] * Xr[47] + Nb[6] * Xr[56] + Nb[7] * Xr[65] +
            Nb[8] * Xr[74]);
  dn[4] = -(Nb[0] * Xr[5] + Nb[1] * Xr[14] + Nb[2] * Xr[23] + Nb[3] * Xr[32] +
            Nb[4] * Xr[41] + Nb[5] * Xr[50] + Nb[6] * Xr[59] + Nb[7] * Xr[68] +
            Nb[8] * Xr[77]);
  dn[7] = -(Nb[0] * Xr[8] + Nb[1] * Xr[17] + Nb[2] * Xr[26] + Nb[3] * Xr[35] +
            Nb[4] * Xr[44] + Nb[5] * Xr[53] + Nb[6] * Xr[62] + Nb[7] * Xr[71] +
            Nb[8] * Xr[80]);

  // Compute zXdinv = -Xdinv*dn*Xdinv
  TacsScalar tmp[9];
  matMatMult(dn, Xdinv, tmp);
  matMatMult(Xdinv, tmp, zXdinv);
}

/*
  Compute the derivative of the matrix zXdinv with respect to the
  inputs. The matrix zXdinv is computed using the following expression
  in index notation:

  zXdinv_{kl} = -Xdinv_{kn}*dn_{nm}*Xdinv_{ml}

  There are two main computations:
  1. The evaluation of dn (the derivative of the normal direction)
  2. The evaluation of zXdinv once dn is known

  The contributions to the sensitivity of dn are:

  d(zXdinv_{kl})/d(dn{ij})
  .  = -Xdinv_{kn}*delta_{ni}*delta_{mj}*Xdinv_{ml}
  .  = -Xdinv_{ki}*Xdinv_{jl}

  dnd_{ij} = zXdinvd_{kl}*d(zXdinv_{kl})/d(d_{ij})
  .  = -zXdinvd_{kl}*Xdinv_{ki}*Xdinv_{jl}

  dnd = -Xdinv^{T}*zXdinvd*Xdinv^{T}

  The second contribution is:

  d(zXdinv_{kl})/d(Xdinv_{ij})
  .  = -Xdinv_{kn}*dn_{nm}*delta_{mi}*delta_{lj}
  .    -delta_{ki}*delta_{nj}*d_{nm}*Xdinv_{ml}

  Xdinvd_{ij} = zXdinvd_{kl}*d(zXdinv_{kl})/d(Xdinv_{ij})
  .  = -Xdinv_{kn}*dn_{ni}*zXdinvd_{kj}
  .    -zXdinvd_{il}*dn_{jm}*Xdinv_{ml}

  dnd += -(Xdinv*dn)^{T}*Xdinvd -Xdinvd*(dn*Xdinv)^{T}
*/
static inline void addNormalRateMatSens(TacsScalar Xrd[], TacsScalar Xdinvd[],
                                        const TacsScalar zXdinvd[],
                                        const double Na[], const double Nb[],
                                        const TacsScalar Xr[],
                                        const TacsScalar Xdinv[]) {
  // Compute the derivatives of the normal direction along the
  // parametric directions
  TacsScalar dn[9];
  dn[2] = dn[5] = dn[8] = 0.0;
  dn[0] = (Na[0] * Xr[2] + Na[1] * Xr[11] + Na[2] * Xr[20] + Na[3] * Xr[29] +
           Na[4] * Xr[38] + Na[5] * Xr[47] + Na[6] * Xr[56] + Na[7] * Xr[65] +
           Na[8] * Xr[74]);
  dn[3] = (Na[0] * Xr[5] + Na[1] * Xr[14] + Na[2] * Xr[23] + Na[3] * Xr[32] +
           Na[4] * Xr[41] + Na[5] * Xr[50] + Na[6] * Xr[59] + Na[7] * Xr[68] +
           Na[8] * Xr[77]);
  dn[6] = (Na[0] * Xr[8] + Na[1] * Xr[17] + Na[2] * Xr[26] + Na[3] * Xr[35] +
           Na[4] * Xr[44] + Na[5] * Xr[53] + Na[6] * Xr[62] + Na[7] * Xr[71] +
           Na[8] * Xr[80]);

  dn[1] = (Nb[0] * Xr[2] + Nb[1] * Xr[11] + Nb[2] * Xr[20] + Nb[3] * Xr[29] +
           Nb[4] * Xr[38] + Nb[5] * Xr[47] + Nb[6] * Xr[56] + Nb[7] * Xr[65] +
           Nb[8] * Xr[74]);
  dn[4] = (Nb[0] * Xr[5] + Nb[1] * Xr[14] + Nb[2] * Xr[23] + Nb[3] * Xr[32] +
           Nb[4] * Xr[41] + Nb[5] * Xr[50] + Nb[6] * Xr[59] + Nb[7] * Xr[68] +
           Nb[8] * Xr[77]);
  dn[7] = (Nb[0] * Xr[8] + Nb[1] * Xr[17] + Nb[2] * Xr[26] + Nb[3] * Xr[35] +
           Nb[4] * Xr[44] + Nb[5] * Xr[53] + Nb[6] * Xr[62] + Nb[7] * Xr[71] +
           Nb[8] * Xr[80]);

  // Compute the contribution
  // -dnd = Xdinv^{T}*zXdinvd*Xdinv^{T}
  TacsScalar dnd[9], t[9];
  matMatTransMult(zXdinvd, Xdinv, t);
  matTransMatMult(Xdinv, t, dnd);

  // Add the contribution to Xdinvd
  // Xdinvd -= (Xdinv*dn)^{T}*Xdinvd + Xdinvd*(dn*Xdinv)^{T}
  dn[0] = -dn[0];
  dn[3] = -dn[3];
  dn[6] = -dn[6];

  dn[1] = -dn[1];
  dn[4] = -dn[4];
  dn[7] = -dn[7];
  matMatMult(Xdinv, dn, t);
  matTransMatMultAdd(t, zXdinvd, Xdinvd);
  matMatMult(dn, Xdinv, t);
  matMatTransMultAdd(zXdinvd, t, Xdinvd);

  // Compute the contribution to Xrd
  for (int i = 0; i < 9; i++) {
    Xrd[2 + 9 * i] -= (dnd[0] * Na[i] + dnd[1] * Nb[i]);
    Xrd[5 + 9 * i] -= (dnd[3] * Na[i] + dnd[4] * Nb[i]);
    Xrd[8 + 9 * i] -= (dnd[6] * Na[i] + dnd[7] * Nb[i]);
  }
}

/*
  Compute the tying shape functions for the element for the
  out-of-plane shear components.

  This uses the quadrature points from a two-point Gauss quadrature
  scheme as the tying points. The tying point shape functions for the
  g13 and g23 shear strains are numbered as follows:

  g13               g23
  +--4-----5--+     +-----------+
  |           |     3     4     5
  |  2  +  3  |     |     +     |
  |           |     0     1     2
  +--0-----1--+     +-----------+

  input:
  u:    the first parametric coordinate
  v:    the second parametric coordinate

  output:
  N13:  the shape functions associated with the g13 shear strain
  N23:  the shape functions associated with the g23 shear strain
*/
static inline void computeTyingFunc(const double u, const double v,
                                    double N13[], double N23[]) {
  // The tying point offset
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;
  const double tinv = 1.0 / t;
  const double sinv = 1.0 / (s * s);

  // Compute the shape functions
  double nu[3], nv[3];
  nu[0] = 0.5 * sinv * u * (u - s);
  nu[1] = sinv * (s - u) * (s + u);
  nu[2] = 0.5 * sinv * u * (s + u);

  nv[0] = 0.5 * sinv * v * (v - s);
  nv[1] = sinv * (s - v) * (s + v);
  nv[2] = 0.5 * sinv * v * (s + v);

  // Compute the shape functions for the reduced dimension
  double ntu[2], ntv[2];
  ntu[0] = 0.5 * tinv * (t - u);
  ntu[1] = 0.5 * tinv * (t + u);

  ntv[0] = 0.5 * tinv * (t - v);
  ntv[1] = 0.5 * tinv * (t + v);

  // Compute the shape functions for g13
  N13[0] = ntu[0] * nv[0];
  N13[1] = ntu[1] * nv[0];
  N13[2] = ntu[0] * nv[1];
  N13[3] = ntu[1] * nv[1];
  N13[4] = ntu[0] * nv[2];
  N13[5] = ntu[1] * nv[2];

  // Compute the shape functions for g23
  N23[0] = nu[0] * ntv[0];
  N23[1] = nu[1] * ntv[0];
  N23[2] = nu[2] * ntv[0];
  N23[3] = nu[0] * ntv[1];
  N23[4] = nu[1] * ntv[1];
  N23[5] = nu[2] * ntv[1];
}

/*
  Given the input vectors that are the derivative of the position
  vector along the u/v directions, and the orthogonal normal vector,
  assemble the results into a reference frame.  Note that this is not
  an orthonormal frame.

  input:
  Xa:      the first in-plane direction
  Xb:      the second in-plane direction
  normal:  the normal direction perpendicular to Xa and Xb

  output:
  Xr:      the reference frame associated with the given point
*/
static inline void assembleFrame(const TacsScalar Xa[], const TacsScalar Xb[],
                                 const TacsScalar normal[], TacsScalar Xr[]) {
  // Fill in the values of Xr
  Xr[0] = Xa[0];
  Xr[3] = Xa[1];
  Xr[6] = Xa[2];

  Xr[1] = Xb[0];
  Xr[4] = Xb[1];
  Xr[7] = Xb[2];

  // Add the values for the normal
  Xr[2] = normal[0];
  Xr[5] = normal[1];
  Xr[8] = normal[2];
}

/*
  Compute the 3x4 matrix from the following matrix-matrix product:

  A = J*S = 2(I - n*n^{T})*[ -eps | (eta*I - eps^{x}) ]

  Note that this is the product of the inertia contribution from the
  director and the angular rate kinematic S matrix.

  input:
  n:     the surface normal vector
  eta:   the quaternion scalar
  eps:   the quaternion vector

  output:
  A:     the product of the matrix (I - n*n^{T}) with S
*/
static inline void computeNormalRateProduct(const TacsScalar n[],
                                            const TacsScalar eta,
                                            const TacsScalar eps[],
                                            TacsScalar A[]) {
  // Compute the rate matrix
  TacsScalar S[12];
  computeSRateMat(eta, eps, S);

  // Pre-multiply the S matrix by (I - n*n^{T})
  TacsScalar ns = 0.0;

  // Compute the first column of A
  ns = n[0] * S[0] + n[1] * S[4] + n[2] * S[8];
  A[0] = S[0] - n[0] * ns;
  A[4] = S[4] - n[1] * ns;
  A[8] = S[8] - n[2] * ns;

  // Compute the second column of A
  ns = n[0] * S[1] + n[1] * S[5] + n[2] * S[9];
  A[1] = S[1] - n[0] * ns;
  A[5] = S[5] - n[1] * ns;
  A[9] = S[9] - n[2] * ns;

  // Compute the third column of A
  ns = n[0] * S[2] + n[1] * S[6] + n[2] * S[10];
  A[2] = S[2] - n[0] * ns;
  A[6] = S[6] - n[1] * ns;
  A[10] = S[10] - n[2] * ns;

  // Compute the third column of A
  ns = n[0] * S[3] + n[1] * S[7] + n[2] * S[11];
  A[3] = S[3] - n[0] * ns;
  A[7] = S[7] - n[1] * ns;
  A[11] = S[11] - n[2] * ns;
}

/*
  Constructor for the MITC9 element class

  input:
  stiff:      the stiffness object
  gravity:    the gravity vector
  vInit:      the initial velocity
  omegaInit:  the initial angular velocity
*/
MITC9::MITC9(FSDTStiffness *_stiff, TACSGibbsVector *_gravity,
             TACSGibbsVector *_vInit, TACSGibbsVector *_omegaInit) {
  // Set the stiffness
  stiff = _stiff;
  stiff->incref();

  // Copy over the vectors (if they are defined)
  gravity = _gravity;
  vInit = _vInit;
  omegaInit = _omegaInit;
  if (gravity) {
    gravity->incref();
  }
  if (vInit) {
    vInit->incref();
  }
  if (omegaInit) {
    omegaInit->incref();
  }

  // Get the Gauss quadrature points and weights
  FElibrary::getGaussPtsWts(ORDER, &gaussPts, &gaussWts);
}

MITC9::~MITC9() {
  stiff->decref();
  if (gravity) {
    gravity->decref();
  }
  if (vInit) {
    vInit->decref();
  }
  if (omegaInit) {
    omegaInit->decref();
  }
}

/*
   Return the number of displacements
*/
int MITC9::getVarsPerNode() { return NUM_DISPS; }

/*
  Return the number of FE nodes
*/
int MITC9::getNumNodes() { return NUM_NODES; }

/*
  Return the ElementLayout
*/
ElementLayout MITC9::getLayoutType() { return TACS_QUAD_QUADRATIC_ELEMENT; }

/*
   Set up the internal static data for the names of the element,
   displacements, stresses, strains and extra variables, respectively.
*/
const char *MITC9::elemName = "MITC9";

/*
  Returns the elementName
*/
const char *MITC9::getObjectName() { return elemName; }

/*
  Get the design variable numbers
*/
int MITC9::getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the design variable values
*/
void MITC9::setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
  stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable values
*/
void MITC9::getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
  stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the design variable range
*/
void MITC9::getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                              TacsScalar ub[]) {
  stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Retrieve the initial values of the design variables
*/
void MITC9::getInitConditions(int elemIndex, const TacsScalar X[],
                              TacsScalar vars[], TacsScalar dvars[],
                              TacsScalar ddvars[]) {
  memset(vars, 0, 8 * NUM_NODES * sizeof(TacsScalar));
  memset(dvars, 0, 8 * NUM_NODES * sizeof(TacsScalar));
  memset(ddvars, 0, 8 * NUM_NODES * sizeof(TacsScalar));

  // The initial quaternions are eta = 1.0, eps = 0
  for (int i = 0; i < NUM_NODES; i++) {
    vars[8 * i + 3] = 1.0;
  }

  // If the initial velocity is defined
  if (vInit) {
    const TacsScalar *v0;
    vInit->getVector(&v0);
    for (int i = 0; i < NUM_NODES; i++) {
      dvars[8 * i] = v0[0];
      dvars[8 * i + 1] = v0[1];
      dvars[8 * i + 2] = v0[2];
    }
  }

  // If the initial angular velocity is defined
  if (omegaInit) {
    const TacsScalar *omega;
    omegaInit->getVector(&omega);

    for (int i = 0; i < NUM_NODES; i++) {
      // dot{u} = v + omega^{x}*r
      crossProductAdd(1.0, omega, &X[3 * i], &dvars[8 * i]);

      // d{eps}/dt = 0.5*omega
      dvars[8 * i + 4] = 0.5 * omega[0];
      dvars[8 * i + 5] = 0.5 * omega[1];
      dvars[8 * i + 6] = 0.5 * omega[2];

      // ddot{u} = omega^{x}*omega^{x}*r
      // Note that the second derivative of the quaternions is zero since
      // there is no angular acceleration
      TacsScalar omegar[3];
      crossProduct(1.0, omega, &X[3 * i], omegar);
      crossProduct(1.0, omega, omegar, &ddvars[8 * i]);
    }
  }
}

/*
  The following function evaluates the kinetic energy and potential
  and elastic energies of the element.

  These can be used to verify that the equations of motion are
  implemented correctly, since the element implements a method based
  on Lagrange's equations of motion.

  input:
  time:   the simulation time
  vars:   the values of the variables
  dvars:  the time-derivative of the variables

  output:
  Te:     the kinetic energy
  Pe:     the potential energy
*/
void MITC9::computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                            const TacsScalar vars[], const TacsScalar dvars[],
                            TacsScalar *Te, TacsScalar *Pe) {
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity) {
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity at the nodes
  TacsScalar omega[3 * NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);

  // Compute the directors at the nodes
  TacsScalar dir[3 * NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Initialize the velocities
  TacsScalar Te = 0.0, Pe = 0.0;

  // Evaluate the kinetic energy of the element
  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9], Xdinv[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i] * gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar v0[3];
      innerProduct8(N, dvars, v0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3];
      innerProduct(N, omega, omeg);

      // Compute the dot product with the normal
      TacsScalar omegn = vecDot(omeg, fn);

      // Add the contributions to the kinetic energy
      Te += 0.5 * h *
            (rho[0] * vecDot(v0, v0) +
             rho[1] * (vecDot(omeg, omeg) - omegn * omegn));

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar rot = computeRotPenalty(N, Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9];
      computeTransform(T, Xa, Xb);

      // Compute the displacement-based strain
      TacsScalar e[8], s[8];
      evalStrain(e, Ur, dr, Xdinv, zXdinv, T);

      // Add the contribution from the tying strain
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);

      // Compute the stress based on the strain values
      TacsScalar A[6], Bc[6], D[6], As[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);
      stiff->calculateStress(A, Bc, D, As, e, s);

      // Compute the terms for the potential energy due to gravity
      TacsScalar U[3];
      innerProduct8(N, vars, U);

      // Add the product of the stress times strain to the
      // potential energy computation
      Pe += 0.5 * h *
            (strainProduct(s, e) + kpenalty * rot * rot -
             2.0 * rho[0] * vecDot(g, U));
    }
  }

  *_Te = Te;
  *_Pe = Pe;
}

/*
  Get the residuals of the equations of motion.

  The equations of motion for this element are derived using
  Lagrange's equations with constraints. Each node imposes a
  constraint that its own quaternions satisfy the required unit norm
  constraint.

  The equations of motion can be divided into two parts: (1) the
  motion of the deformed surface and (2) the motion of the normals of
  the surface (or directors since they are no longer normal to the
  deformed surface during deformation). The linear motion takes the
  form:

  M*ddot{u} + dU/dx - fg = 0

  where M is a mass matrix. The rotational degrees of freedom satisfy
  the following equations of motion:

  S^{T}*J*d{omega} + 2*dot{S}^{T}*J*omega + dU/dq - A^{T}*lamb = 0

  where J = (I - n*n^{T}) is a rotational inertia term and the
  constraints produce the term A^{T}*lamb.

  input:
  X:      the element nodal coordinates
  vars:   the element degrees of freedom
  dvars:  the first time derivative of the degrees of freedom
  ddvars: the second time derivative of the degrees of freedom

  output:
  res:    the residuals
*/
void MITC9::addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar res[]) {
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity) {
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3 * NUM_NODES], domega[3 * NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Compute the area for this element: this is used to scale the
  // constraint equations within the code
  TacsScalar area = 0.0;

  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar Xdinv[9];
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i] * gaussWts[j];

      // Compute the area
      area += h;

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar a0[3];
      innerProduct8(N, ddvars, a0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp * fn[0];
      w[1] = omeg[1] - tmp * fn[1];
      w[2] = omeg[2] - tmp * fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp * fn[0];
      dw[1] = domeg[1] - tmp * fn[1];
      dw[2] = domeg[2] - tmp * fn[2];

      // Add the contribution to the residual
      TacsScalar *r = res;
      const TacsScalar *q = vars, *dq = dvars;
      for (int ii = 0; ii < NUM_NODES; ii++) {
        // Add the contributions from the rectilinear velocity
        r[0] += h * N[ii] * rho[0] * a0[0];
        r[1] += h * N[ii] * rho[0] * a0[1];
        r[2] += h * N[ii] * rho[0] * a0[2];

        // Add the contributions from the angular velocity
        // S^{T}*dw + 2*dot{S}^{T}*w
        TacsScalar eta = q[3];
        const TacsScalar *eps = &q[4];
        TacsScalar deta = dq[3];
        const TacsScalar *deps = &dq[4];

        // Add S^{T}*dw
        r[3] -= 2.0 * h * N[ii] * rho[1] * vecDot(eps, dw);
        crossProductAdd(2.0 * h * N[ii] * rho[1], eps, dw, &r[4]);
        vecAxpy(2.0 * h * N[ii] * eta * rho[1], dw, &r[4]);

        // Add 2*dot{S}^{T}*w
        r[3] -= 4.0 * h * N[ii] * rho[1] * vecDot(deps, w);
        crossProductAdd(4.0 * h * N[ii] * rho[1], deps, w, &r[4]);
        vecAxpy(4.0 * h * N[ii] * deta * rho[1], w, &r[4]);

        r += 8;
        q += 8;
        dq += 8;
      }

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar Brot[8 * NUM_NODES];
      TacsScalar rot =
          computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9];
      computeTransform(T, Xa, Xb);

      // Compute the displacement-based strain
      TacsScalar e[8], B[64 * NUM_NODES];
      evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

      // Add the contribution from the tying straint
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
      addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

      // Compute the stress based on the strain values
      TacsScalar s[8];  // The stress components
      TacsScalar A[6], Bc[6], D[6], As[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);
      stiff->calculateStress(A, Bc, D, As, e, s);

      // Scale the rotation by the in-plane penalty term
      rot *= kpenalty;

      // Add the contribution to the residual
      r = res;
      const TacsScalar *b = B, *br = Brot;
      for (int ii = 0; ii < NUM_NODES; ii++) {
        r[0] +=
            h * (strainProduct(s, &b[0]) + rot * br[0] - rho[0] * N[ii] * g[0]);
        r[1] +=
            h * (strainProduct(s, &b[8]) + rot * br[1] - rho[0] * N[ii] * g[1]);
        r[2] += h * (strainProduct(s, &b[16]) + rot * br[2] -
                     rho[0] * N[ii] * g[2]);
        r[3] += h * (strainProduct(s, &b[24]) + rot * br[3]);
        r[4] += h * (strainProduct(s, &b[32]) + rot * br[4]);
        r[5] += h * (strainProduct(s, &b[40]) + rot * br[5]);
        r[6] += h * (strainProduct(s, &b[48]) + rot * br[6]);

        r += 8;
        b += 64;
        br += 8;
      }
    }
  }

  // Set the scaling for the constraints
  TacsScalar scale = area;

  // Add the constraints from the quaternion parametrization
  for (int i = 0; i < NUM_NODES; i++) {
    const TacsScalar *q = &vars[8 * i + 3];
    TacsScalar lamb = vars[8 * i + 7];

    // Add the result to the governing equations
    res[8 * i + 3] += 2.0 * scale * q[0] * lamb;
    res[8 * i + 4] += 2.0 * scale * q[1] * lamb;
    res[8 * i + 5] += 2.0 * scale * q[2] * lamb;
    res[8 * i + 6] += 2.0 * scale * q[3] * lamb;

    // Enforce the quaternion constraint
    res[8 * i + 7] +=
        scale * (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1.0);
  }
}

/*
  Add the Jacobian for the governing equations

  The following code computes the Jacobian of the equations of motion
  involving a linear combination of the Jacobian w.r.t. the variables
  and their first and second time derivatives as follows:

  J = alpha*dR/dq + beta*dR/d(dot{q}) + gamma*dR/d(ddot{q})

  The alpha, beta, and gamma coefficients are used for the values and
  their first and second time-derivatives respectively.  This function
  is used to assemble the element-contribution to the multibody system
  Jacobian.

  input:
  alpha:   the variable coefficient
  beta:    the first-derivative coefficient
  gamma:   the second-derivative coefficient
  X:       the nodal locations (of the initial configuration)
  vars:    the variable values
  dvars:   the time derivative of the variable values
  ddvars:  the second time derivative of the variable values

  output:
  J:       the Jacobian matrix
*/
void MITC9::addJacobian(int elemIndex, double time, TacsScalar alpha,
                        TacsScalar beta, TacsScalar gamma,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TacsScalar res[], TacsScalar mat[]) {
  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3 * NUM_NODES], domega[3 * NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Compute the area of this element: used to evaluate the
  // constraints
  TacsScalar area = 0.0;

  // The weights that are used for the geometric stiffness from
  // the tying strain
  TacsScalar w13[6], w23[6];
  memset(w13, 0, 6 * sizeof(TacsScalar));
  memset(w23, 0, 6 * sizeof(TacsScalar));

  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9], Xdinv[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i] * gaussWts[j];

      // Compute the area
      area += h;

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // Add the contributions from the linear motion
      for (int ii = 0; ii < NUM_NODES; ii++) {
        for (int jj = 0; jj < NUM_NODES; jj++) {
          const TacsScalar scale = gamma * h * N[ii] * N[jj] * rho[0];
          // Add the contributions from the rectilinear velocity
          J[8 * NUM_NODES * (8 * ii) + 8 * jj] += scale;
          J[8 * NUM_NODES * (8 * ii + 1) + 8 * jj + 1] += scale;
          J[8 * NUM_NODES * (8 * ii + 2) + 8 * jj + 2] += scale;
        }
      }

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp * fn[0];
      w[1] = omeg[1] - tmp * fn[1];
      w[2] = omeg[2] - tmp * fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp * fn[0];
      dw[1] = domeg[1] - tmp * fn[1];
      dw[2] = domeg[2] - tmp * fn[2];

      // Add the contributions from the rotational DOF
      for (int ii = 0; ii < NUM_NODES; ii++) {
        // Compute and store the product of (I - n*n^{T}) with the
        // angular rate matrices
        TacsScalar JSii[12], JdSii[12];
        computeNormalRateProduct(fn, vars[8 * ii + 3], &vars[8 * ii + 4], JSii);
        computeNormalRateProduct(fn, dvars[8 * ii + 3], &dvars[8 * ii + 4],
                                 JdSii);

        // Set the pointer to the Jacobian entries that will
        // be added
        TacsScalar *Jp = &J[(8 * NUM_NODES + 1) * (8 * ii + 3)];
        const int ldj = 8 * NUM_NODES;

        // Add the diagonal terms
        const TacsScalar dscale = h * N[ii] * rho[1];
        addSRateMatTransDeriv(alpha * dscale, dw, Jp, ldj);
        addSRateMatTransDeriv(2.0 * beta * dscale, w, Jp, ldj);

        // Add the result to the Jacobian matrix
        const TacsScalar *q = vars, *dq = dvars, *ddq = ddvars;
        for (int jj = 0; jj < NUM_NODES; jj++) {
          // Set the common scaling factor for all terms
          const TacsScalar scale = h * N[ii] * N[jj] * rho[1];

          // Set the pointer to the Jacobian entries that will
          // be added
          TacsScalar *Jp = &J[8 * NUM_NODES * (8 * ii + 3) + 8 * jj + 3];

          // Compute S = S(q)
          TacsScalar Sjj[12];
          computeSRateMat(q[3], &q[4], Sjj);

          // Compute dot{S} = S(dot{q})
          TacsScalar dSjj[12];
          computeSRateMat(dq[3], &dq[4], dSjj);

          // Compute ddot{S} = S(ddot{q})
          TacsScalar ddSjj[12];
          computeSRateMat(ddq[3], &ddq[4], ddSjj);

          // Add the Jacobian terms from the DOF:
          // T(dw) - S^{T}*J*S(ddot{q}) - 2*dot{S}^{T}*J*dot{S}
          addBlock3x4Product(-alpha * scale, JSii, ddSjj, Jp, ldj);
          addBlock3x4Product(-2.0 * alpha * scale, JdSii, dSjj, Jp, ldj);

          // Add the Jacobian terms from the first time derivatives:
          // 2*dot{S}^{T}*J*S
          addBlock3x4Product(2.0 * beta * scale, JdSii, Sjj, Jp, ldj);

          // Add the Jacobian terms from the second time derivatives:
          // S^{T}*J*S
          addBlock3x4Product(gamma * scale, JSii, Sjj, Jp, ldj);

          q += 8;
          dq += 8;
          ddq += 8;
        }
      }

      // This code is expensive, and is only required when alpha
      // is non-zero
      if (alpha != 0.0) {
        // Everything left in this function is proportional to
        // h*alpha, so we scale h by alpha
        h *= alpha;

        // Evaluate the stiffness terms at the parametric point
        TacsScalar A[6], Bc[6], D[6], As[3];
        TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);

        // Evaluate the tying strain interpolation
        double N13[6], N23[6];
        computeTyingFunc(u, v, N13, N23);

        // Compute the through-thickness derivative of [X,r]^{-1}
        TacsScalar zXdinv[9];
        computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

        // Compute the derivatives of Ua/Ub along the given directions
        TacsScalar Ur[9], Ua[3], Ub[3], d[3];
        innerProduct8(Na, vars, Ua);
        innerProduct8(Nb, vars, Ub);
        innerProduct(N, dir, d);
        assembleFrame(Ua, Ub, d, Ur);

        // Now compute the derivatives of the director along each
        // coordinate direction
        TacsScalar dr[9], da[3], db[3], zero[3];
        innerProduct(Na, dir, da);
        innerProduct(Nb, dir, db);
        zero[0] = zero[1] = zero[2] = 0.0;
        assembleFrame(da, db, zero, dr);

        // Add the in-plane penalty terms
        TacsScalar Brot[8 * NUM_NODES];
        TacsScalar rot =
            computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

        // Add the contribution from the penalty
        addGRotMat(J, h * kpenalty * rot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

        // Compute the transformation to the locally-aligned frame
        TacsScalar T[9];
        computeTransform(T, Xa, Xb);

        // Compute the displacement-based strain
        TacsScalar e[8], B[64 * NUM_NODES];
        evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

        // Add the contribution from the tying straint
        addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
        addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

        // Compute the stress based on the strain values
        TacsScalar s[8];  // The stress components
        stiff->calculateStress(A, Bc, D, As, e, s);

        // Add to the weights
        addTyingGmatWeights(w13, w23, h, s, N13, N23, Xdinv, T);

        // Add the geometric stiffness terms
        addGmat(J, h, s, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, Xr, dirdq);

        // Add the contribution to the residual
        for (int ii = 0; ii < NUM_NODES; ii++) {
          for (int ik = 0; ik < 7; ik++) {
            // Compute the stress from the 8*i + ik component
            TacsScalar sbii[8];
            stiff->calculateStress(A, Bc, D, As, &B[8 * (8 * ii + ik)], sbii);

            // Compute the
            TacsScalar pr = kpenalty * Brot[8 * ii + ik];

            const TacsScalar *b = B, *br = Brot;
            TacsScalar *Jp = &J[8 * NUM_NODES * (8 * ii + ik)];
            for (int jj = 0; jj < NUM_NODES; jj++) {
              Jp[0] += h * (strainProduct(sbii, &b[0]) + pr * br[0]);
              Jp[1] += h * (strainProduct(sbii, &b[8]) + pr * br[1]);
              Jp[2] += h * (strainProduct(sbii, &b[16]) + pr * br[2]);
              Jp[3] += h * (strainProduct(sbii, &b[24]) + pr * br[3]);
              Jp[4] += h * (strainProduct(sbii, &b[32]) + pr * br[4]);
              Jp[5] += h * (strainProduct(sbii, &b[40]) + pr * br[5]);
              Jp[6] += h * (strainProduct(sbii, &b[48]) + pr * br[6]);

              Jp += 8;
              b += 64;
              br += 8;
            }
          }
        }
      }
    }
  }

  // Add the geometric stiffness terms from the tying strain
  addTyingGmat(J, w13, w23, X, Xr, vars, dir, dirdq);

  // Set the scaling for the constraints
  TacsScalar scale = 2.0 * alpha * area;

  // Add the constraints from the quaternions
  for (int i = 0; i < NUM_NODES; i++) {
    const TacsScalar *q = &vars[8 * i + 3];
    const int ldj = 8 * NUM_NODES;

    TacsScalar *Jp = &J[(8 * NUM_NODES + 1) * (8 * i + 3)];
    TacsScalar lamb = vars[8 * i + 7];

    // Add the constraint terms
    Jp[4] += scale * q[0];
    Jp[4 + ldj] += scale * q[1];
    Jp[4 + 2 * ldj] += scale * q[2];
    Jp[4 + 3 * ldj] += scale * q[3];

    // Enforce the quaternion constraint
    Jp[4 * ldj] += scale * q[0];
    Jp[4 * ldj + 1] += scale * q[1];
    Jp[4 * ldj + 2] += scale * q[2];
    Jp[4 * ldj + 3] += scale * q[3];

    // Add the terms to the diagonal
    Jp[0] += scale * lamb;
    Jp[ldj + 1] += scale * lamb;
    Jp[2 * (ldj + 1)] += scale * lamb;
    Jp[3 * (ldj + 1)] += scale * lamb;
  }
}

/*
  Add the derivative of the product of the adjoint variables with the
  residuals to the given design variable vector.

  input:
  time:     the simulation time
  scale:    the scalar factor applied to the result
  fdvSens:  the derivative (times scale) is accumulated in this array
  dvLen:    the design array length
  psi:      the adjoint vector
  X:        the node locations
  vars:     the variable values
  dvars:    the time derivatives of the variable values
  ddvars:   the second time derivative of the variable values
*/
void MITC9::addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                             const TacsScalar psi[], const TacsScalar Xpts[],
                             const TacsScalar vars[], const TacsScalar dvars[],
                             const TacsScalar ddvars[], int dvLen,
                             TacsScalar dfdx[]) {
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity) {
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3 * NUM_NODES], domega[3 * NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar Xdinv[9];
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i] * gaussWts[j];

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar a0[3];
      innerProduct8(N, ddvars, a0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar tmp, w[3], dw[3];
      tmp = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp * fn[0];
      w[1] = omeg[1] - tmp * fn[1];
      w[2] = omeg[2] - tmp * fn[2];

      tmp = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp * fn[0];
      dw[1] = domeg[1] - tmp * fn[1];
      dw[2] = domeg[2] - tmp * fn[2];

      // Add the contribution to the residual
      TacsScalar mscale[2] = {0.0, 0.0};
      const TacsScalar *p = psi;
      const TacsScalar *q = vars, *dq = dvars;
      for (int ii = 0; ii < NUM_NODES; ii++) {
        // Add the contributions from the rectilinear velocity
        mscale[0] += h * N[ii] * (a0[0] * p[0] + a0[1] * p[1] + a0[2] * p[2]);

        // Add the contributions from the angular velocity
        // S^{T}*dw + 2*dot{S}^{T}*w
        TacsScalar eta = q[3];
        const TacsScalar *eps = &q[4];
        TacsScalar deta = dq[3];
        const TacsScalar *deps = &dq[4];

        // Add p^{T}*S^{T}*dw
        TacsScalar t[3];
        crossProduct(1.0, eps, dw, t);
        mscale[1] -= 2.0 * p[3] * h * N[ii] * vecDot(eps, dw);
        mscale[1] +=
            2.0 * h * N[ii] * (eta * vecDot(dw, &p[4]) + vecDot(t, &p[4]));

        // Add p^{T}*2*dot{S}^{T}*w
        crossProduct(1.0, deps, w, t);
        mscale[1] -= 4.0 * p[3] * h * N[ii] * vecDot(deps, w);
        mscale[1] +=
            4.0 * h * N[ii] * (deta * vecDot(w, &p[4]) + vecDot(t, &p[4]));

        // Increment the pointers to the variables/multipliers
        q += 8;
        dq += 8;
        p += 8;
      }

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the term due to the potential energy
      TacsScalar Brot[8 * NUM_NODES];
      TacsScalar rot =
          computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9];
      computeTransform(T, Xa, Xb);

      // Compute the displacement-based strain
      TacsScalar e[8], B[64 * NUM_NODES];
      evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

      // Add the contribution from the tying strain
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
      addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

      // Set the Lagrange multiplier associated with the strain
      TacsScalar epsi[8];
      TacsScalar rotPsi = 0.0;
      epsi[0] = epsi[1] = epsi[2] = epsi[3] = epsi[4] = epsi[5] = epsi[6] =
          epsi[7] = 0.0;

      const TacsScalar *b = B, *br = Brot;
      for (int ii = 0; ii < 8 * NUM_NODES; ii++) {
        epsi[0] += b[0] * psi[ii];
        epsi[1] += b[1] * psi[ii];
        epsi[2] += b[2] * psi[ii];
        epsi[3] += b[3] * psi[ii];
        epsi[4] += b[4] * psi[ii];
        epsi[5] += b[5] * psi[ii];
        epsi[6] += b[6] * psi[ii];
        epsi[7] += b[7] * psi[ii];
        rotPsi += br[0] * psi[ii];
        b += 8;
        br++;
      }

      // Add the contribution from the gravity load
      for (int ii = 0; ii < NUM_NODES; ii++) {
        mscale[0] -= h * N[ii] *
                     (g[0] * psi[8 * ii] + g[1] * psi[8 * ii + 1] +
                      g[2] * psi[8 * ii + 2]);
      }
      mscale[0] *= scale;
      mscale[1] *= scale;

      // Scale the psi vector by the determinant of the Jacobian
      // transformation
      for (int k = 0; k < 8; k++) {
        epsi[k] *= h * scale;
      }

      // Add the derivative contribution from the mass/area
      stiff->addMassMomentsDVSens(pt, mscale, dvLen, dfdx);

      // Add the derivative
      stiff->addStressDVSens(elemIndex, pt, Xpt, e, h * scale * rot * rotPsi,
                             epsi, dvLen, dfdx);
    }
  }
}

/*
  Compute the derivative of the product of the adjoint variables and
  the residual vector.

  This code uses a reverse-mode method, which makes the derivative
  code faster, but a bit more difficult to follow.

  input:
  time:      the simulation time
  scale:     scale the derivative vector by this factor
  psi:       the adjoint variables associated with this element
  vars:      the state variables
  dvars:     the first derivative of the state variables
  ddvars:    the second derivative of the state variables

  output:
  fXptSens:  add the derivative to this vector
*/
void MITC9::addAdjResXptProduct(double time, double scale,
                                TacsScalar fXptSens[], const TacsScalar psi[],
                                const TacsScalar X[], const TacsScalar vars[],
                                const TacsScalar dvars[],
                                const TacsScalar ddvars[]) {
  // Set the gravity vector - if one exists
  TacsScalar g[3] = {0.0, 0.0, 0.0};
  if (gravity) {
    const TacsScalar *_g;
    gravity->getVector(&_g);
    g[0] = _g[0];
    g[1] = _g[1];
    g[2] = _g[2];
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // The derivatives w.r.t. reference frames
  TacsScalar Xrd[9 * NUM_NODES];
  memset(Xrd, 0, 9 * NUM_NODES * sizeof(TacsScalar));

  // Compute the angular velocity and acceleration at the nodes
  TacsScalar omega[3 * NUM_NODES], domega[3 * NUM_NODES];
  computeAngularVelocity(omega, vars, dvars);
  computeAngularAccel(domega, vars, dvars, ddvars);

  // Compute the derivatives of the directors
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Memory to keep track of the derivatives w.r.t. directors
  TacsScalar dird[3 * NUM_NODES], dirdqd[12 * NUM_NODES];
  memset(dird, 0, 3 * NUM_NODES * sizeof(TacsScalar));
  memset(dirdqd, 0, 12 * NUM_NODES * sizeof(TacsScalar));

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Memory to keep track of the derivative contributions
  TacsScalar g13d[6], g23d[6];
  TacsScalar B13d[6 * 8 * NUM_NODES], B23d[6 * 8 * NUM_NODES];
  memset(g13d, 0, 6 * sizeof(TacsScalar));
  memset(g23d, 0, 6 * sizeof(TacsScalar));
  memset(B13d, 0, 6 * 8 * NUM_NODES * sizeof(TacsScalar));
  memset(B23d, 0, 6 * 8 * NUM_NODES * sizeof(TacsScalar));

  // Add up the contributions to the derivative of the area
  TacsScalar aread = 0.0;

  // Add the constraints from the quaternion parametrization
  for (int i = 0; i < NUM_NODES; i++) {
    const TacsScalar *q = &vars[8 * i + 3];
    TacsScalar lamb = vars[8 * i + 7];
    const TacsScalar *p = &psi[8 * i + 3];

    // Enforce the quaternion constraint
    aread +=
        (p[7] * (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3] - 1.0) +
         2.0 * lamb *
             (psi[0] * q[0] + psi[1] * q[1] + psi[2] * q[2] + psi[3] * q[3]));
  }

  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Set the Gauss quadrature points
      const double u = gaussPts[i];
      const double v = gaussPts[j];

      // The parametric point required for evaluating
      // stresses/strains within the element
      const double pt[2] = {u, v};

      // Evaluate the shape functions
      double N[NUM_NODES];
      computeShapeFunc(u, v, N);

      // Evaluate the derivatives of the shape functions
      double Na[NUM_NODES], Nb[NUM_NODES];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the frame normal
      TacsScalar fn[3];
      computeFrameNormal(N, Xr, fn);

      // Evaluate the derivatives in the locally-aligned frame
      TacsScalar Xd[9];
      assembleFrame(Xa, Xb, fn, Xd);

      // Compute the derivatives of the shape functions
      TacsScalar Xdinv[9];
      TacsScalar h = inv3x3(Xd, Xdinv);
      h *= gaussWts[i] * gaussWts[j];

      // Add the contribution from the derivative of the inner product
      // w.r.t. h (the determinant of the Jacobian matrix)
      TacsScalar hd = aread;

      // Evaluate the areal mass properties
      TacsScalar rho[2];
      stiff->getPointwiseMass(pt, rho);

      // The following is used to evaluate the kinetic energy
      // ---------------------------------------------------
      // Evaluate the velocity at the quadrature point
      TacsScalar a0[3];
      innerProduct8(N, ddvars, a0);

      // Compute the value of omega at the current point
      TacsScalar omeg[3], domeg[3];
      innerProduct(N, omega, omeg);
      innerProduct(N, domega, domeg);

      // Remove the normal component from angular velocity/accel.
      // w = (I - fn*fn^{T})*omega
      TacsScalar w[3], dw[3];
      TacsScalar tmp1 = vecDot(omeg, fn);
      w[0] = omeg[0] - tmp1 * fn[0];
      w[1] = omeg[1] - tmp1 * fn[1];
      w[2] = omeg[2] - tmp1 * fn[2];

      TacsScalar tmp2 = vecDot(domeg, fn);
      dw[0] = domeg[0] - tmp2 * fn[0];
      dw[1] = domeg[1] - tmp2 * fn[1];
      dw[2] = domeg[2] - tmp2 * fn[2];

      // Accumulate the derivative contributions
      TacsScalar wd[3], dwd[3];
      wd[0] = wd[1] = wd[2] = 0.0;
      dwd[0] = dwd[1] = dwd[2] = 0.0;

      // Add the contribution to the residual
      const TacsScalar *p = psi;
      const TacsScalar *q = vars, *dq = dvars;
      for (int ii = 0; ii < NUM_NODES; ii++) {
        // Add the contributions from the rectilinear velocity
        hd += N[ii] * rho[0] * (p[0] * a0[0] + p[1] * a0[1] + p[2] * a0[2]);

        // Add the contributions from the angular velocity
        // S^{T}*dw + 2*dot{S}^{T}*w
        TacsScalar eta = q[3];
        const TacsScalar *eps = &q[4];
        TacsScalar deta = dq[3];
        const TacsScalar *deps = &dq[4];

        // Add S^{T}*dw
        TacsScalar scl = 2.0 * N[ii] * rho[1];
        TacsScalar t[3];
        crossProduct(1.0, eps, dw, t);
        hd -= scl * p[3] * vecDot(eps, dw);
        hd += scl * vecDot(t, &p[4]);
        hd += scl * eta * vecDot(dw, &p[4]);

        // Add 2*dot{S}^{T}*w
        crossProduct(1.0, deps, w, t);
        hd -= 2.0 * scl * p[3] * vecDot(deps, w);
        hd += 2.0 * scl * vecDot(t, &p[4]);
        hd += 2.0 * scl * deta * vecDot(w, &p[4]);

        // Add the contributions to wd and dwd
        scl *= h * scale;
        vecAxpy(-scl * p[3], eps, dwd);
        crossProductAdd(scl, &p[4], eps, dwd);
        vecAxpy(scl * eta, &p[4], dwd);

        vecAxpy(-2.0 * scl * p[3], deps, wd);
        crossProductAdd(2.0 * scl, &p[4], deps, wd);
        vecAxpy(2.0 * scl * deta, &p[4], wd);

        p += 8;
        q += 8;
        dq += 8;
      }

      TacsScalar fnd[3];
      TacsScalar tmpw = vecDot(wd, fn);
      TacsScalar tmpdw = vecDot(dwd, fn);
      fnd[0] =
          -(tmp1 * wd[0] + tmp2 * dwd[0] + tmpw * omeg[0] + tmpdw * domeg[0]);
      fnd[1] =
          -(tmp1 * wd[1] + tmp2 * dwd[1] + tmpw * omeg[1] + tmpdw * domeg[1]);
      fnd[2] =
          -(tmp1 * wd[2] + tmp2 * dwd[2] + tmpw * omeg[2] + tmpdw * domeg[2]);

      // The following code is used to evaluate the potential energy
      // -----------------------------------------------------------
      // Evaluate the tying strain interpolation
      double N13[6], N23[6];
      computeTyingFunc(u, v, N13, N23);

      // Compute the through-thickness derivative of [X,r]^{-1}
      TacsScalar zXdinv[9];
      computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

      // Compute the derivatives of Ua/Ub along the given directions
      TacsScalar Ur[9], Ua[3], Ub[3], d[3];
      innerProduct8(Na, vars, Ua);
      innerProduct8(Nb, vars, Ub);
      innerProduct(N, dir, d);
      assembleFrame(Ua, Ub, d, Ur);

      // Now compute the derivatives of the director along each
      // coordinate direction
      TacsScalar dr[9], da[3], db[3], zero[3];
      innerProduct(Na, dir, da);
      innerProduct(Nb, dir, db);
      zero[0] = zero[1] = zero[2] = 0.0;
      assembleFrame(da, db, zero, dr);

      // Add the terms from the potential energy
      TacsScalar Brot[8 * NUM_NODES];
      TacsScalar rot =
          computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

      // Compute the transformation to the locally-aligned frame
      TacsScalar T[9];
      computeTransform(T, Xa, Xb);

      // Compute the displacement-based strain
      TacsScalar e[8], B[64 * NUM_NODES];
      evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

      // Add the contribution from the tying straint
      addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
      addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

      // Compute the stress based on the strain values
      TacsScalar s[8];  // The stress components
      TacsScalar A[6], Bc[6], D[6], As[3];
      TacsScalar kpenalty = stiff->getStiffness(pt, A, Bc, D, As);
      stiff->calculateStress(A, Bc, D, As, e, s);

      // Scale the rotation by the in-plane penalty term
      rot *= kpenalty;

      // Compute the product Bpsi = B*psi
      TacsScalar Bpsi[NUM_STRESSES];
      memset(Bpsi, 0, NUM_STRESSES * sizeof(TacsScalar));

      // Compute the product BRpsi = Brot*psi
      TacsScalar BRpsi = 0.0;

      // Set the pointer to the psi variables back to the begining
      p = psi;

      // Compute the products
      const TacsScalar *b = B;
      const TacsScalar *br = Brot;
      for (int k = 0; k < 8 * NUM_NODES; k++) {
        Bpsi[0] += b[0] * p[0];
        Bpsi[1] += b[1] * p[0];
        Bpsi[2] += b[2] * p[0];
        Bpsi[3] += b[3] * p[0];
        Bpsi[4] += b[4] * p[0];
        Bpsi[5] += b[5] * p[0];
        Bpsi[6] += b[6] * p[0];
        Bpsi[7] += b[7] * p[0];
        BRpsi += br[0] * p[0];
        hd += p[0] * (strainProduct(s, b) + rot * br[0]);
        p++;
        br++;
        b += NUM_STRESSES;
      }

      // Add the contribution from the gravity vector
      for (int k = 0; k < NUM_NODES; k++) {
        hd -=
            rho[0] * N[k] *
            (g[0] * psi[8 * k] + g[1] * psi[8 * k + 1] + g[2] * psi[8 * k + 2]);
      }

      // Scale by the derivative w.r.t. h by the quadrature weight
      hd *= gaussWts[i] * gaussWts[j] * scale;

      // Compute the product spsi = C*B*psi
      TacsScalar spsi[NUM_STRESSES];
      stiff->calculateStress(A, Bc, D, As, Bpsi, spsi);

      // Add the contribution from the derivative of the
      // term h*e^{T}*C*B*psi with respect to the stress
      TacsScalar Urd[9], drd[9];
      TacsScalar Xdinvd[9], zXdinvd[9];
      TacsScalar Td[9];
      evalStrainSens(Urd, drd, Xdinvd, zXdinvd, Td, h * scale, spsi, Ur, dr,
                     Xdinv, zXdinv, T);
      addTyingStrainSens(g13d, g23d, Xdinvd, Td, h * scale, spsi, N13, N23, g13,
                         g23, Xdinv, T);

      // Compute the quantity hs = h*s;
      TacsScalar hs[NUM_STRESSES];
      for (int k = 0; k < NUM_STRESSES; k++) {
        hs[k] = s[k] * h * scale;
      }

      // Add the contribution from the derivative of the
      // term h*s^{T}*B*psi w.r.t. B
      addBmatSens(Urd, drd, Xdinvd, zXdinvd, Td, dirdqd, hs, psi, N, Na, Nb, Ur,
                  dr, Xdinv, zXdinv, T, dirdq);
      addTyingBmatSens(B13d, B23d, Xdinvd, Td, hs, psi, N13, N23, B13, B23,
                       Xdinv, T);

      // Add the contribution from the director computation
      for (int k = 0; k < NUM_NODES; k++) {
        dird[3 * k] += Na[k] * drd[0] + Nb[k] * drd[1] + N[k] * Urd[2];
        dird[3 * k + 1] += Na[k] * drd[3] + Nb[k] * drd[4] + N[k] * Urd[5];
        dird[3 * k + 2] += Na[k] * drd[6] + Nb[k] * drd[7] + N[k] * Urd[8];
      }

      // Compute the derivatives Xad and Xbd
      TacsScalar Xad[3], Xbd[3];
      computeTransformSens(Xad, Xbd, Td, Xa, Xb);

      // Add the contribution from the rotation
      TacsScalar rotd = scale * h * BRpsi * kpenalty;
      addBRotPenaltySens(Xad, Xbd, rotd, scale * h * rot, psi, N, Na, Nb, Xa,
                         Xb, Ua, Ub, vars);

      // Compute the derivatives Xrd
      addNormalRateMatSens(Xrd, Xdinvd, zXdinvd, Na, Nb, Xr, Xdinv);

      // Compute the derivative of the inverse of the Jacobian
      TacsScalar Xdd[9];
      inv3x3Sens(Xdd, Xdinvd, Xdinv);

      // Add the contributions from the derivative of the determinant
      TacsScalar hXdd[9];
      det3x3Sens(Xd, hXdd);
      for (int k = 0; k < 9; k++) {
        Xdd[k] += hXdd[k] * hd;
      }

      // Extract/add the sensitivities from the frame
      fnd[0] += Xdd[2];
      fnd[1] += Xdd[5];
      fnd[2] += Xdd[8];

      // Add the contributions to Xad and Xbd
      Xad[0] += Xdd[0];
      Xad[1] += Xdd[3];
      Xad[2] += Xdd[6];
      Xbd[0] += Xdd[1];
      Xbd[1] += Xdd[4];
      Xbd[2] += Xdd[7];

      // // Compute the frame normal
      addFrameNormalSens(fnd, N, Xrd);

      // Add the derivatives the shape function directions
      for (int k = 0; k < NUM_NODES; k++) {
        fXptSens[3 * k] += Na[k] * Xad[0] + Nb[k] * Xbd[0];
        fXptSens[3 * k + 1] += Na[k] * Xad[1] + Nb[k] * Xbd[1];
        fXptSens[3 * k + 2] += Na[k] * Xad[2] + Nb[k] * Xbd[2];
      }
    }
  }

  addComputeTyingBmatSens(fXptSens, Xrd, dird, dirdqd, g13d, g23d, B13d, B23d,
                          X, Xr, vars, dir, dirdq);
  addDirectorsSens(Xrd, dird, vars);
  addDirectorDerivSens(Xrd, dirdqd, vars);
  addFramesSens(fXptSens, Xrd, X);
}

/*
  Given the nodal degrees of freedom and their time-derivatives,
  compute the angular velocity at each node
*/
void MITC9::computeAngularVelocity(TacsScalar omega[], const TacsScalar vars[],
                                   const TacsScalar dvars[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    TacsScalar eta = vars[3];
    const TacsScalar *eps = &vars[4];
    TacsScalar deta = dvars[3];
    const TacsScalar *deps = &dvars[4];

    // omega = -2*eps^{x}*deps + 2*eta*deps - eps*deta
    crossProduct(-2.0, eps, deps, omega);
    vecAxpy(2.0 * eta, deps, omega);
    vecAxpy(-2.0 * deta, eps, omega);

    omega += 3;
    vars += 8;
    dvars += 8;
  }
}

/*
  Given the nodal degrees of freedom and their first and second
  time-derivatives, compute the angular acceleration at the nodes.
*/
void MITC9::computeAngularAccel(TacsScalar domega[], const TacsScalar vars[],
                                const TacsScalar dvars[],
                                const TacsScalar ddvars[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    // Set pointers to the values
    TacsScalar eta = vars[3];
    const TacsScalar *eps = &vars[4];
    TacsScalar ddeta = ddvars[3];
    const TacsScalar *ddeps = &ddvars[4];

    // domega = S(q)*ddot{q}
    crossProduct(-2.0, eps, ddeps, domega);
    vecAxpy(2.0 * eta, ddeps, domega);
    vecAxpy(-2.0 * ddeta, eps, domega);

    domega += 3;
    vars += 8;
    dvars += 8;
    ddvars += 8;
  }
}

/*
  At each node in the finite-element, compute the derivatives of the
  coordinate directions and assemble a locally-aligned reference
  frame.

  Each locally-aligned reference frame Xr[] consists of a 3x3 matrix
  stored in row-major order where the first direction is aligned
  along the 1-direction and the third direction is normal to the local
  surface and the second direction is perpendicular to both these
  directions. This can be written as follows:

  Xr = [X,xi1; X,xi2; n]

  input:
  X:    the initial nodal locations

  output:
  Xr:   the locally-aligned frames
*/
void MITC9::computeFrames(TacsScalar Xr[], const TacsScalar X[]) {
  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Find the u/v values at the node locations
      double u = -1.0 + 2.0 * i / (ORDER - 1.0);
      double v = -1.0 + 2.0 * j / (ORDER - 1.0);

      // Evaluate the shape functions
      double Na[9], Nb[9];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the cross product
      TacsScalar normal[3];
      crossProduct(1.0, Xa, Xb, normal);
      TacsScalar nrm = sqrt(vecDot(normal, normal));
      vecScale(1.0 / nrm, normal);

      // Assemble the frame into the Xr matrix
      assembleFrame(Xa, Xb, normal, Xr);

      // Increment the pointer to the frames
      Xr += 9;
    }
  }
}

/*
  Add the derivative of the frame computation to the output vector

  This code adds the contribution to the derivative of the frame
  computation with respect to the node locations X.
*/
void MITC9::addFramesSens(TacsScalar Xd[], const TacsScalar Xrd[],
                          const TacsScalar X[]) {
  for (int j = 0; j < ORDER; j++) {
    for (int i = 0; i < ORDER; i++) {
      // Find the u/v values at the node locations
      double u = -1.0 + 2.0 * i / (ORDER - 1.0);
      double v = -1.0 + 2.0 * j / (ORDER - 1.0);

      // Evaluate the shape functions
      double Na[9], Nb[9];
      computeShapeFunc(u, v, Na, Nb);

      // Compute the derivative along the shape function
      // directions
      TacsScalar Xa[3], Xb[3];
      innerProduct(Na, X, Xa);
      innerProduct(Nb, X, Xb);

      // Compute the cross product
      TacsScalar n[3];
      crossProduct(1.0, Xa, Xb, n);
      TacsScalar invNrm = 1.0 / sqrt(vecDot(n, n));

      // Compute the normal
      TacsScalar Xad[3], Xbd[3], nhatd[3];
      Xad[0] = Xrd[0];
      Xad[1] = Xrd[3];
      Xad[2] = Xrd[6];
      Xbd[0] = Xrd[1];
      Xbd[1] = Xrd[4];
      Xbd[2] = Xrd[7];
      nhatd[0] = Xrd[2];
      nhatd[1] = Xrd[5];
      nhatd[2] = Xrd[8];

      // Compute the sensitivity contribution from n = nhat/||nhat||
      TacsScalar nd[3];
      TacsScalar tmp = invNrm * invNrm * invNrm;
      nd[0] = tmp * ((n[1] * n[1] + n[2] * n[2]) * nhatd[0] -
                     n[0] * (n[1] * nhatd[1] + n[2] * nhatd[2]));
      nd[1] = tmp * ((n[0] * n[0] + n[2] * n[2]) * nhatd[1] -
                     n[1] * (n[0] * nhatd[0] + n[2] * nhatd[2]));
      nd[2] = tmp * ((n[0] * n[0] + n[1] * n[1]) * nhatd[2] -
                     n[2] * (n[0] * nhatd[0] + n[1] * nhatd[1]));

      // Add the derivative contribution to nd and Xbd
      crossProductAdd(1.0, Xb, nd, Xad);
      crossProductAdd(1.0, nd, Xa, Xbd);

      // Add the contributions to the derivative
      for (int k = 0; k < NUM_NODES; k++) {
        Xd[3 * k] += Na[k] * Xad[0] + Nb[k] * Xbd[0];
        Xd[1 + 3 * k] += Na[k] * Xad[1] + Nb[k] * Xbd[1];
        Xd[2 + 3 * k] += Na[k] * Xad[2] + Nb[k] * Xbd[2];
      }

      // Increment the pointer to the frames
      Xrd += 9;
    }
  }
}

/*
  Given the non-orthonormal in-plane directions Xa/Xb which are the
  parametric tangents, compute the transformation to the
  locally-aligned frame for the strain.

  input:
  Xa:   tangent to the first parametric direction
  Xb:   tangent to the second parametric direction

  output:
  T:    the transformation matrix
*/
void MITC9::computeTransform(TacsScalar T[], const TacsScalar Xa[],
                             const TacsScalar Xb[]) {
  // Compute the cross product of Xa/Xb to find the normal direction
  TacsScalar normal[3];
  crossProduct(1.0, Xa, Xb, normal);
  TacsScalar nrm = sqrt(vecDot(normal, normal));
  vecScale(1.0 / nrm, normal);

  // The first coordinate direction
  TacsScalar d1[3];
  if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
    // Scale the Xa direction so that it is a unit vector
    TacsScalar invXa = 1.0 / sqrt(vecDot(Xa, Xa));
    d1[0] = invXa * Xa[0];
    d1[1] = invXa * Xa[1];
    d1[2] = invXa * Xa[2];
  } else {
    // Retrieve the reference axis from the constitutive
    // object
    const TacsScalar *axis = stiff->getRefAxis();
    TacsScalar an = vecDot(axis, normal);

    // Take the component of the reference axis perpendicular
    // to the surface
    d1[0] = axis[0] - an * normal[0];
    d1[1] = axis[1] - an * normal[1];
    d1[2] = axis[2] - an * normal[2];

    // Normalize the new direction
    TacsScalar inv = 1.0 / sqrt(vecDot(d1, d1));
    vecScale(inv, d1);
  }

  // Compute the second perpendicular direction
  TacsScalar d2[3];
  crossProduct(1.0, normal, d1, d2);

  // Assemble the transformation matrix
  assembleFrame(d1, d2, normal, T);
}

/*
  Compute the derivative of the transformation
*/
void MITC9::computeTransformSens(TacsScalar Xad[], TacsScalar Xbd[],
                                 const TacsScalar Td[], const TacsScalar Xa[],
                                 const TacsScalar Xb[]) {
  // Compute the cross product of Xa/Xb to find the normal direction
  TacsScalar n[3];
  crossProduct(1.0, Xa, Xb, n);

  TacsScalar invNrm = 1.0 / sqrt(vecDot(n, n));
  TacsScalar nhat[3];
  nhat[0] = invNrm * n[0];
  nhat[1] = invNrm * n[1];
  nhat[2] = invNrm * n[2];

  // The first coordinate direction
  if (stiff->getTransformType() == FSDTStiffness::NATURAL) {
    // Scale the Xa direction so that it is a unit vector
    TacsScalar invXa = 1.0 / sqrt(vecDot(Xa, Xa));
    TacsScalar d1[3];
    d1[0] = invXa * Xa[0];
    d1[1] = invXa * Xa[1];
    d1[2] = invXa * Xa[2];

    // Compute the second perpendicular direction
    TacsScalar d2[3];
    crossProduct(1.0, nhat, d1, d2);

    // Extract the seeds from the Td input array
    TacsScalar d1d[3], d2d[3];
    TacsScalar nhatd[3];  // normal seed
    d1d[0] = Td[0];
    d1d[1] = Td[3];
    d1d[2] = Td[6];
    d2d[0] = Td[1];
    d2d[1] = Td[4];
    d2d[2] = Td[7];
    nhatd[0] = Td[2];
    nhatd[1] = Td[5];
    nhatd[2] = Td[8];

    // Add the derivative contribution from d2 to nd and d1d
    crossProductAdd(1.0, d2d, nhat, d1d);
    crossProductAdd(1.0, d1, d2d, nhatd);

    // Compute df/dT*dT/d(d1)*d(d1/d(Xa))
    TacsScalar tmp = invXa * invXa * invXa;
    Xad[0] = tmp * ((Xa[1] * Xa[1] + Xa[2] * Xa[2]) * d1d[0] -
                    Xa[0] * (Xa[1] * d1d[1] + Xa[2] * d1d[2]));
    Xad[1] = tmp * ((Xa[0] * Xa[0] + Xa[2] * Xa[2]) * d1d[1] -
                    Xa[1] * (Xa[0] * d1d[0] + Xa[2] * d1d[2]));
    Xad[2] = tmp * ((Xa[0] * Xa[0] + Xa[1] * Xa[1]) * d1d[2] -
                    Xa[2] * (Xa[0] * d1d[0] + Xa[1] * d1d[1]));

    // Compute n = nhat/||nhat||
    TacsScalar nd[3];
    tmp = invNrm * invNrm * invNrm;
    nd[0] = tmp * ((n[1] * n[1] + n[2] * n[2]) * nhatd[0] -
                   n[0] * (n[1] * nhatd[1] + n[2] * nhatd[2]));
    nd[1] = tmp * ((n[0] * n[0] + n[2] * n[2]) * nhatd[1] -
                   n[1] * (n[0] * nhatd[0] + n[2] * nhatd[2]));
    nd[2] = tmp * ((n[0] * n[0] + n[1] * n[1]) * nhatd[2] -
                   n[2] * (n[0] * nhatd[0] + n[1] * nhatd[1]));

    // Add the derivative contribution to nd and Xbd
    crossProductAdd(1.0, Xb, nd, Xad);
    crossProduct(1.0, nd, Xa, Xbd);
  } else {
    // Retrieve the reference axis from the constitutive
    // object
    const TacsScalar *axis = stiff->getRefAxis();
    TacsScalar an = vecDot(axis, nhat);

    // Take the component of the reference axis perpendicular
    // to the surface
    TacsScalar t[3];
    t[0] = axis[0] - an * nhat[0];
    t[1] = axis[1] - an * nhat[1];
    t[2] = axis[2] - an * nhat[2];

    // Normalize the new direction
    TacsScalar d1[3];
    TacsScalar inv = 1.0 / sqrt(vecDot(t, t));
    d1[0] = inv * t[0];
    d1[1] = inv * t[1];
    d1[2] = inv * t[2];

    // Extract the seeds from the Td input array
    TacsScalar d1d[3], d2d[3];
    TacsScalar nhatd[3];  // normal seed
    d1d[0] = Td[0];
    d1d[1] = Td[3];
    d1d[2] = Td[6];
    d2d[0] = Td[1];
    d2d[1] = Td[4];
    d2d[2] = Td[7];
    nhatd[0] = Td[2];
    nhatd[1] = Td[5];
    nhatd[2] = Td[8];

    // Add the derivative contribution from d2 to nd and d1d
    crossProductAdd(1.0, d2d, nhat, d1d);
    crossProductAdd(1.0, d1, d2d, nhatd);

    // Compute df/dT*dT/d(d1)*d(d1/d(Xa))
    TacsScalar tmp = inv * inv * inv;
    TacsScalar td[3];
    td[0] = tmp * ((t[1] * t[1] + t[2] * t[2]) * d1d[0] -
                   t[0] * (t[1] * d1d[1] + t[2] * d1d[2]));
    td[1] = tmp * ((t[0] * t[0] + t[2] * t[2]) * d1d[1] -
                   t[1] * (t[0] * d1d[0] + t[2] * d1d[2]));
    td[2] = tmp * ((t[0] * t[0] + t[1] * t[1]) * d1d[2] -
                   t[2] * (t[0] * d1d[0] + t[1] * d1d[1]));

    TacsScalar ta = vecDot(td, nhat);
    nhatd[0] -= ta * axis[0] + an * td[0];
    nhatd[1] -= ta * axis[1] + an * td[1];
    nhatd[2] -= ta * axis[2] + an * td[2];

    // Compute n = nhat/||nhat||
    TacsScalar nd[3];
    tmp = invNrm * invNrm * invNrm;
    nd[0] = tmp * ((n[1] * n[1] + n[2] * n[2]) * nhatd[0] -
                   n[0] * (n[1] * nhatd[1] + n[2] * nhatd[2]));
    nd[1] = tmp * ((n[0] * n[0] + n[2] * n[2]) * nhatd[1] -
                   n[1] * (n[0] * nhatd[0] + n[2] * nhatd[2]));
    nd[2] = tmp * ((n[0] * n[0] + n[1] * n[1]) * nhatd[2] -
                   n[2] * (n[0] * nhatd[0] + n[1] * nhatd[1]));

    // Add the derivative contribution to nd and Xbd
    crossProduct(1.0, Xb, nd, Xad);
    crossProduct(1.0, nd, Xa, Xbd);
  }
}

/*
  Based on the values of the element variables, determine the director
  values.

  The director is the difference between the deformed normal direction
  and the initial normal. The directors are used to determine the
  through-thickness strain distribution in the finite-element model.
  The advantage of the director approach is that it enables a
  geometrically exact displacement representation and also allows for
  shell-shell intersections at 90 degree angles.

  input:
  vars:   the element variable values
  Xr:     the local frames at each node in the mesh

  output:
  d:      the director at each node in the finite-element
*/
void MITC9::computeDirectors(TacsScalar d[], const TacsScalar vars[],
                             const TacsScalar Xr[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute C = rot - I
    TacsScalar C[9];
    C[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
    C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
    C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

    C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
    C[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
    C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

    C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
    C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
    C[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

    // Compute d = C^{T}*n
    d[0] = C[0] * Xr[2] + C[3] * Xr[5] + C[6] * Xr[8];
    d[1] = C[1] * Xr[2] + C[4] * Xr[5] + C[7] * Xr[8];
    d[2] = C[2] * Xr[2] + C[5] * Xr[5] + C[8] * Xr[8];

    d += 3;     // Each director is a 3-vector
    Xr += 9;    // Increment over each frame
    vars += 8;  // 8 variables per node
  }
}

/*
  Add the derivative of the directors to the local frame
*/
void MITC9::addDirectorsSens(TacsScalar Xrd[], const TacsScalar dd[],
                             const TacsScalar vars[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute C = rot - I
    TacsScalar C[9];
    C[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
    C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
    C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

    C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
    C[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
    C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

    C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
    C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
    C[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

    Xrd[2] += C[0] * dd[0] + C[1] * dd[1] + C[2] * dd[2];
    Xrd[5] += C[3] * dd[0] + C[4] * dd[1] + C[5] * dd[2];
    Xrd[8] += C[6] * dd[0] + C[7] * dd[1] + C[8] * dd[2];

    dd += 3;    // Each director is a 3-vector
    Xrd += 9;   // Increment over each frame
    vars += 8;  // 8 variables per node
  }
}

/*
  Compute the derivative of the director values w.r.t. the rotation
  matrix parametrization.

  This code computes the derivative of (C^{T} - I)*n w.r.t. the
  quaternions at each node in the finite-element. This is required
  to compute the governing equations of motion for the element.

  input:
  vars:  the element variables
  Xr:    the local frame for each node in the mesh

  output:
  ddq:   the derivative of the directors w.r.t. the quaterions
*/
void MITC9::computeDirectorDeriv(TacsScalar ddq[], const TacsScalar vars[],
                                 const TacsScalar Xr[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute the derivative of the rotation matrix w.r.t.
    // the quaternions
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0 * q[3];
    Q[2] = -2.0 * q[2];

    Q[3] = -2.0 * q[3];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[1];

    Q[6] = 2.0 * q[2];
    Q[7] = -2.0 * q[1];
    Q[8] = 0.0;

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0] * Xr[2] + Q[3] * Xr[5] + Q[6] * Xr[8];
    ddq[1] = Q[1] * Xr[2] + Q[4] * Xr[5] + Q[7] * Xr[8];
    ddq[2] = Q[2] * Xr[2] + Q[5] * Xr[5] + Q[8] * Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0 * q[2];
    Q[2] = 2.0 * q[3];

    Q[3] = 2.0 * q[2];
    Q[4] = -4.0 * q[1];
    Q[5] = 2.0 * q[0];

    Q[6] = 2.0 * q[3];
    Q[7] = -2.0 * q[0];
    Q[8] = -4.0 * q[1];

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0] * Xr[2] + Q[3] * Xr[5] + Q[6] * Xr[8];
    ddq[1] = Q[1] * Xr[2] + Q[4] * Xr[5] + Q[7] * Xr[8];
    ddq[2] = Q[2] * Xr[2] + Q[5] * Xr[5] + Q[8] * Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[2]
    Q[0] = -4.0 * q[2];
    Q[1] = 2.0 * q[1];
    Q[2] = -2.0 * q[0];

    Q[3] = 2.0 * q[1];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[3];

    Q[6] = 2.0 * q[0];
    Q[7] = 2.0 * q[3];
    Q[8] = -4.0 * q[2];

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0] * Xr[2] + Q[3] * Xr[5] + Q[6] * Xr[8];
    ddq[1] = Q[1] * Xr[2] + Q[4] * Xr[5] + Q[7] * Xr[8];
    ddq[2] = Q[2] * Xr[2] + Q[5] * Xr[5] + Q[8] * Xr[8];
    ddq += 3;

    // Derivative w.r.t. q[3]
    Q[0] = -4.0 * q[3];
    Q[1] = 2.0 * q[0];
    Q[2] = 2.0 * q[1];

    Q[3] = -2.0 * q[0];
    Q[4] = -4.0 * q[3];
    Q[5] = 2.0 * q[2];

    Q[6] = 2.0 * q[1];
    Q[7] = 2.0 * q[2];
    Q[8] = 0.0;

    // Compute ddq = D^{T}*n
    ddq[0] = Q[0] * Xr[2] + Q[3] * Xr[5] + Q[6] * Xr[8];
    ddq[1] = Q[1] * Xr[2] + Q[4] * Xr[5] + Q[7] * Xr[8];
    ddq[2] = Q[2] * Xr[2] + Q[5] * Xr[5] + Q[8] * Xr[8];
    ddq += 3;

    Xr += 9;    // Increment over each frame
    vars += 8;  // 8 variables per node
  }
}

/*
  Add the sensitivity with respect to the nodes from the computation
  of the derivative of the director values w.r.t. the rotation matrix
  parametrization.
*/
void MITC9::addDirectorDerivSens(TacsScalar Xrd[], const TacsScalar ddqd[],
                                 const TacsScalar vars[]) {
  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the
    const TacsScalar *q = &vars[3];

    // Compute the derivative of the rotation matrix w.r.t.
    // the quaternions
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0 * q[3];
    Q[2] = -2.0 * q[2];

    Q[3] = -2.0 * q[3];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[1];

    Q[6] = 2.0 * q[2];
    Q[7] = -2.0 * q[1];
    Q[8] = 0.0;

    Xrd[2] += Q[0] * ddqd[0] + Q[1] * ddqd[1] + Q[2] * ddqd[2];
    Xrd[5] += Q[3] * ddqd[0] + Q[4] * ddqd[1] + Q[5] * ddqd[2];
    Xrd[8] += Q[6] * ddqd[0] + Q[7] * ddqd[1] + Q[8] * ddqd[2];
    ddqd += 3;

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0 * q[2];
    Q[2] = 2.0 * q[3];

    Q[3] = 2.0 * q[2];
    Q[4] = -4.0 * q[1];
    Q[5] = 2.0 * q[0];

    Q[6] = 2.0 * q[3];
    Q[7] = -2.0 * q[0];
    Q[8] = -4.0 * q[1];

    Xrd[2] += Q[0] * ddqd[0] + Q[1] * ddqd[1] + Q[2] * ddqd[2];
    Xrd[5] += Q[3] * ddqd[0] + Q[4] * ddqd[1] + Q[5] * ddqd[2];
    Xrd[8] += Q[6] * ddqd[0] + Q[7] * ddqd[1] + Q[8] * ddqd[2];
    ddqd += 3;

    // Derivative w.r.t. q[2]
    Q[0] = -4.0 * q[2];
    Q[1] = 2.0 * q[1];
    Q[2] = -2.0 * q[0];

    Q[3] = 2.0 * q[1];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[3];

    Q[6] = 2.0 * q[0];
    Q[7] = 2.0 * q[3];
    Q[8] = -4.0 * q[2];

    Xrd[2] += Q[0] * ddqd[0] + Q[1] * ddqd[1] + Q[2] * ddqd[2];
    Xrd[5] += Q[3] * ddqd[0] + Q[4] * ddqd[1] + Q[5] * ddqd[2];
    Xrd[8] += Q[6] * ddqd[0] + Q[7] * ddqd[1] + Q[8] * ddqd[2];
    ddqd += 3;

    // Derivative w.r.t. q[3]
    Q[0] = -4.0 * q[3];
    Q[1] = 2.0 * q[0];
    Q[2] = 2.0 * q[1];

    Q[3] = -2.0 * q[0];
    Q[4] = -4.0 * q[3];
    Q[5] = 2.0 * q[2];

    Q[6] = 2.0 * q[1];
    Q[7] = 2.0 * q[2];
    Q[8] = 0.0;

    Xrd[2] += Q[0] * ddqd[0] + Q[1] * ddqd[1] + Q[2] * ddqd[2];
    Xrd[5] += Q[3] * ddqd[0] + Q[4] * ddqd[1] + Q[5] * ddqd[2];
    Xrd[8] += Q[6] * ddqd[0] + Q[7] * ddqd[1] + Q[8] * ddqd[2];
    ddqd += 3;

    Xrd += 9;   // Increment over each frame
    vars += 8;  // 8 variables per node
  }
}

/*
  Given the derivatives of the displacements and the transformations
  to the local coordinates, evaluate strain.

  The expressions for the strain involve the surface derivatives of
  the displacements along the coordinate directions and the derivative
  of the director along the surface coordinate directions:

  Ur = [ U0,xi; d(xi) ]
  dr = [ d,xi; 0 ]

  The geometry enters through both an surface-based transformation and
  a derivative of the transformation through the thickness of the
  shell as follows:

  Xdinv = [X,r]^{-1}
  zXdinv = -[X,r]^{-1}*[n,xi]*[X,r]^{-1}

  input:
  Ur:     the derivative of the displacements w.r.t. all shell coordinates
  dr:     the derivative of the director w.r.t. all shell coordinates
  Xdinv:  the inverse of the Jacobian matrix
  zXdinv: the derivative of the inverse of the Jacobian w.r.t. z

  output:
  e:      the displacement-based components of the strain
*/
void MITC9::evalStrain(TacsScalar e[], const TacsScalar Ur[],
                       const TacsScalar dr[], const TacsScalar Xdinv[],
                       const TacsScalar zXdinv[], const TacsScalar T[]) {
  TacsScalar U0[9], U1[9], tmp[9];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Compute the in-plane strain
  e[0] = U0[0] + 0.5 * (U0[0] * U0[0] + U0[3] * U0[3] + U0[6] * U0[6]);
  e[1] = U0[4] + 0.5 * (U0[1] * U0[1] + U0[4] * U0[4] + U0[7] * U0[7]);
  e[2] = U0[1] + U0[3] + (U0[0] * U0[1] + U0[3] * U0[4] + U0[6] * U0[7]);

  // Compute the bending strain
  e[3] = U1[0] + (U0[0] * U1[0] + U0[3] * U1[3] + U0[6] * U1[6]);
  e[4] = U1[4] + (U0[1] * U1[1] + U0[4] * U1[4] + U0[7] * U1[7]);
  e[5] = U1[1] + U1[3] +
         (U0[0] * U1[1] + U0[3] * U1[4] + U0[6] * U1[7] + U1[0] * U0[1] +
          U1[3] * U0[4] + U1[6] * U0[7]);
}

/*
  Compute the derivative of the product of the strain with an input
  vector with respect to all of the function inputs.

  input:
  scale:   scalar multiple for all components
  eSens:   input sensitivity vector
  Ur:      derivative of displacements w.r.t. shell coordinates
  dr:      derivative of the direction w.r.t. shell coordinates
  Xdinv:   the inverse of the Jacobian matrix
  zXdinv:  the derivative of the inverse of the Jacobian w.r.t. z

  output:
  Urd:     derivative of the function w.r.t. Ur
  drd:     derivative of the function w.r.t. dr
  Xdinvd:  derivative of the function w.r.t. Xdinv
  zXdinvd: derivative of the function w.r.t. Xdinvd
*/
void MITC9::evalStrainSens(TacsScalar Urd[], TacsScalar drd[],
                           TacsScalar Xdinvd[], TacsScalar zXdinvd[],
                           TacsScalar Td[], TacsScalar scale,
                           const TacsScalar eSens[], const TacsScalar Ur[],
                           const TacsScalar dr[], const TacsScalar Xdinv[],
                           const TacsScalar zXdinv[], const TacsScalar T[]) {
  // Compute U0 = T^{T}*Ur*Xdinv*T
  TacsScalar UrXdinv[9], UrXdinvT[9];
  matMatMult(Ur, Xdinv, UrXdinv);
  matMatMult(UrXdinv, T, UrXdinvT);

  TacsScalar U0[9];
  matTransMatMult(T, UrXdinvT, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  TacsScalar tsum[9];
  matMatMult(Ur, zXdinv, tsum);
  matMatMultAdd(dr, Xdinv, tsum);

  TacsScalar U1[9], drXdinvT[9];
  matMatMult(tsum, T, drXdinvT);
  matTransMatMult(T, drXdinvT, U1);

  // Compute the derivative of the strain product w.r.t. U0
  TacsScalar U0d[9];
  U0d[0] = scale * (eSens[0] * (1.0 + U0[0]) + eSens[2] * U0[1] +
                    eSens[3] * U1[0] + eSens[5] * U1[1]);
  U0d[1] = scale * (eSens[2] * (1.0 + U0[0]) + eSens[1] * U0[1] +
                    eSens[4] * U1[1] + eSens[5] * U1[0]);
  U0d[2] = 0.0;
  U0d[3] = scale * (eSens[2] * (1.0 + U0[4]) + eSens[0] * U0[3] +
                    eSens[3] * U1[3] + eSens[5] * U1[4]);
  U0d[4] = scale * (eSens[1] * (1.0 + U0[4]) + eSens[2] * U0[3] +
                    eSens[4] * U1[4] + eSens[5] * U1[3]);
  U0d[5] = 0.0;
  U0d[6] = scale * (eSens[0] * U0[6] + eSens[2] * U0[7] + eSens[3] * U1[6] +
                    eSens[5] * U1[7]);
  U0d[7] = scale * (eSens[1] * U0[7] + eSens[2] * U0[6] + eSens[4] * U1[7] +
                    eSens[5] * U1[6]);
  U0d[8] = 0.0;

  // Compute the derivative with respect to U1
  TacsScalar U1d[9];
  U1d[0] = scale * (eSens[3] * (1.0 + U0[0]) + eSens[5] * U0[1]);
  U1d[1] = scale * (eSens[5] * (1.0 + U0[0]) + eSens[4] * U0[1]);
  U1d[2] = 0.0;

  U1d[3] = scale * (eSens[5] * (1.0 + U0[4]) + eSens[3] * U0[3]);
  U1d[4] = scale * (eSens[4] * (1.0 + U0[4]) + eSens[5] * U0[3]);
  U1d[5] = 0.0;

  U1d[6] = scale * (eSens[3] * U0[6] + eSens[5] * U0[7]);
  U1d[7] = scale * (eSens[4] * U0[7] + eSens[5] * U0[6]);
  U1d[8] = 0.0;

  // Compute the derivative w.r.t. Ur
  TacsScalar t[9];    // Temporary matrix
  TacsScalar t0d[9];  // = T*U0d*T^{T}
  matMatTransMult(U0d, T, Urd);
  matMatMult(T, Urd, t0d);
  matMatTransMult(t0d, Xdinv, Urd);

  TacsScalar t1d[9];  // = T*U1d*T^{T}
  matMatTransMult(U1d, T, t);
  matMatMult(T, t, t1d);
  matMatTransMultAdd(t1d, zXdinv, Urd);

  // Compute the contribution to the derivative w.r.t. dr
  matMatTransMult(t1d, Xdinv, drd);

  // Compute the contribution to the derivative w.r.t. Xdinv
  matTransMatMult(dr, t1d, Xdinvd);
  matTransMatMultAdd(Ur, t0d, Xdinvd);

  // Compute the contribution to the derivative w.r.t. zXdinv
  matTransMatMult(Ur, t1d, zXdinvd);

  // Compute the contribution to the derivative w.r.t. T
  matMatTransMult(UrXdinvT, U0d, Td);  // Td = Ur*Xdinv*T*U0d^{T}
  matTransMatMult(UrXdinv, T, t);
  matMatMultAdd(t, U0d, Td);  // Td += (UrXdinv)^{T}*T*Ud0

  matMatTransMultAdd(drXdinvT, U1d, Td);  // Td = Ur*Xdinv*T*U0d^{T}
  matTransMatMult(tsum, T, t);
  matMatMultAdd(t, U1d, Td);  // Td += (UrXdinv)^{T}*T*Ud0
}

/*
  Evaluate the strain and the derivative of the strain w.r.t. the
  element variables.

  This function is used in the evaluation of the residuals and in the
  computation of the Jacobian matrix.

  input:
  N, Na, Nb:   shape functions and their derivatives
  Ur:          derivatives of displacements in the parametric direction
  dr:          derivative of the director field
  Xdinv:       Jacobian transformation
  zXdinv:      through-thickness derivative of the Jacobian transformation
  T:           transformation to the local strain axis
  dirdq:       derivative of the director field with respect to states

  output:
  e:           the strain
  B;           the derivative of the strain w.r.t. the variables
*/
void MITC9::evalBmat(TacsScalar e[], TacsScalar B[], const double N[],
                     const double Na[], const double Nb[],
                     const TacsScalar Ur[], const TacsScalar dr[],
                     const TacsScalar Xdinv[], const TacsScalar zXdinv[],
                     const TacsScalar T[], const TacsScalar dirdq[]) {
  TacsScalar U0[9], U1[9], tmp[9];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Compute the in-plane strain
  e[0] = U0[0] + 0.5 * (U0[0] * U0[0] + U0[3] * U0[3] + U0[6] * U0[6]);
  e[1] = U0[4] + 0.5 * (U0[1] * U0[1] + U0[4] * U0[4] + U0[7] * U0[7]);
  e[2] = U0[1] + U0[3] + (U0[0] * U0[1] + U0[3] * U0[4] + U0[6] * U0[7]);

  // Compute the bending strain
  e[3] = U1[0] + (U0[0] * U1[0] + U0[3] * U1[3] + U0[6] * U1[6]);
  e[4] = U1[4] + (U0[1] * U1[1] + U0[4] * U1[4] + U0[7] * U1[7]);
  e[5] = U1[1] + U1[3] +
         (U0[0] * U1[1] + U0[3] * U1[4] + U0[6] * U1[7] + U1[0] * U0[1] +
          U1[3] * U0[4] + U1[6] * U0[7]);

  // Compute the derivatives of the strain w.r.t. the displacement
  // variables. This code takes advantage of the sparsity of the
  // derivatives to simplify the computations
  for (int i = 0; i < NUM_NODES; i++) {
    TacsScalar *b = &B[8 * (8 * i)];

    // Compute the values of the displacements
    TacsScalar dU0[9], dU1[9], ztmp[3];

    // Compute dU0 = d(Ur)/dq_{k} * Xdinv
    // [  0,  0, 0 ][ 0  1  2 ][ T0  T1  T2 ]
    // [ Na, Nb, 0 ][ 3  4  5 ][ T3  T4  T5 ]
    // [  0,  0, 0 ][ 6  7  8 ][ T6  T7  T8 ]
    dU0[0] = Na[i] * Xdinv[0] + Nb[i] * Xdinv[3];
    dU0[1] = Na[i] * Xdinv[1] + Nb[i] * Xdinv[4];
    dU0[2] = Na[i] * Xdinv[2] + Nb[i] * Xdinv[5];
    matMultTrans(T, dU0, tmp);

    dU1[0] = Na[i] * zXdinv[0] + Nb[i] * zXdinv[3];
    dU1[1] = Na[i] * zXdinv[1] + Nb[i] * zXdinv[4];
    dU1[2] = Na[i] * zXdinv[2] + Nb[i] * zXdinv[5];
    matMultTrans(T, dU1, ztmp);

    for (int k = 0; k < 3; k++) {
      // dU0 = T[3*k+i]*tmp[3*k+j];
      dU0[0] = T[3 * k] * tmp[0];
      dU0[1] = T[3 * k] * tmp[1];
      dU0[2] = T[3 * k] * tmp[2];
      dU0[3] = T[3 * k + 1] * tmp[0];
      dU0[4] = T[3 * k + 1] * tmp[1];
      dU0[5] = T[3 * k + 1] * tmp[2];
      dU0[6] = T[3 * k + 2] * tmp[0];
      dU0[7] = T[3 * k + 2] * tmp[1];
      dU0[8] = T[3 * k + 2] * tmp[2];

      // dU1 = T[3*k+i]*ztmp[3*k+j];
      dU1[0] = T[3 * k] * ztmp[0];
      dU1[1] = T[3 * k] * ztmp[1];
      dU1[2] = T[3 * k] * ztmp[2];
      dU1[3] = T[3 * k + 1] * ztmp[0];
      dU1[4] = T[3 * k + 1] * ztmp[1];
      dU1[5] = T[3 * k + 1] * ztmp[2];
      dU1[6] = T[3 * k + 2] * ztmp[0];
      dU1[7] = T[3 * k + 2] * ztmp[1];
      dU1[8] = T[3 * k + 2] * ztmp[2];

      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0] * dU0[0] + U0[3] * dU0[3] + U0[6] * dU0[6]);
      b[1] = dU0[4] + (U0[1] * dU0[1] + U0[4] * dU0[4] + U0[7] * dU0[7]);
      b[2] = dU0[1] + dU0[3] +
             (U0[0] * dU0[1] + U0[3] * dU0[4] + U0[6] * dU0[7] +
              dU0[0] * U0[1] + dU0[3] * U0[4] + dU0[6] * U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0] * dU1[0] + U0[3] * dU1[3] + U0[6] * dU1[6] +
                       dU0[0] * U1[0] + dU0[3] * U1[3] + dU0[6] * U1[6]);
      b[4] = dU1[4] + (U0[1] * dU1[1] + U0[4] * dU1[4] + U0[7] * dU1[7] +
                       dU0[1] * U1[1] + dU0[4] * U1[4] + dU0[7] * U1[7]);
      b[5] =
          dU1[1] + dU1[3] +
          (U0[0] * dU1[1] + U0[3] * dU1[4] + U0[6] * dU1[7] + U1[0] * dU0[1] +
           U1[3] * dU0[4] + U1[6] * dU0[7] + dU0[0] * U1[1] + dU0[3] * U1[4] +
           dU0[6] * U1[7] + dU1[0] * U0[1] + dU1[3] * U0[4] + dU1[6] * U0[7]);
      b[6] = b[7] = 0.0;
      b += 8;
    }
  }

  // Add the contributions from the derivative of the strain w.r.t. the
  // rotation variables. These derivatives only make contributions to the
  // bending strains, not the in-plane strains.
  for (int i = 0; i < NUM_NODES; i++) {
    TacsScalar *b = &B[8 * (8 * i + 3)];

    // Compute the values of the displacements
    TacsScalar dU0[9], dU1[9], drdq[9];

    for (int k = 0; k < 4; k++) {
      // Compute dU0 = T^{T}*dUr*Xdinv*T
      // T^{T}*[ 0 | 0 | N*dirdq[0] ]*Xdinv*T
      // .     [ 0 | 0 | N*dirdq[1] ]*Xdinv*T
      //       [ 0 | 0 | N*dirdq[2] ]*Xdinv*T
      dU0[0] = N[i] * dirdq[0] * Xdinv[6];
      dU0[1] = N[i] * dirdq[0] * Xdinv[7];
      dU0[2] = N[i] * dirdq[0] * Xdinv[8];
      dU0[3] = N[i] * dirdq[1] * Xdinv[6];
      dU0[4] = N[i] * dirdq[1] * Xdinv[7];
      dU0[5] = N[i] * dirdq[1] * Xdinv[8];
      dU0[6] = N[i] * dirdq[2] * Xdinv[6];
      dU0[7] = N[i] * dirdq[2] * Xdinv[7];
      dU0[8] = N[i] * dirdq[2] * Xdinv[8];
      matMatMult(dU0, T, tmp);
      matTransMatMult(T, tmp, dU0);

      // Compute the derivative for d
      drdq[0] = Na[i] * dirdq[0];
      drdq[3] = Na[i] * dirdq[1];
      drdq[6] = Na[i] * dirdq[2];
      drdq[1] = Nb[i] * dirdq[0];
      drdq[4] = Nb[i] * dirdq[1];
      drdq[7] = Nb[i] * dirdq[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, dU1);

      dU1[0] += N[i] * dirdq[0] * zXdinv[6];
      dU1[1] += N[i] * dirdq[0] * zXdinv[7];
      dU1[2] += N[i] * dirdq[0] * zXdinv[8];
      dU1[3] += N[i] * dirdq[1] * zXdinv[6];
      dU1[4] += N[i] * dirdq[1] * zXdinv[7];
      dU1[5] += N[i] * dirdq[1] * zXdinv[8];
      dU1[6] += N[i] * dirdq[2] * zXdinv[6];
      dU1[7] += N[i] * dirdq[2] * zXdinv[7];
      dU1[8] += N[i] * dirdq[2] * zXdinv[8];
      matMatMult(dU1, T, tmp);
      matTransMatMult(T, tmp, dU1);

      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0] * dU0[0] + U0[3] * dU0[3] + U0[6] * dU0[6]);
      b[1] = dU0[4] + (U0[1] * dU0[1] + U0[4] * dU0[4] + U0[7] * dU0[7]);
      b[2] = dU0[1] + dU0[3] +
             (U0[0] * dU0[1] + U0[3] * dU0[4] + U0[6] * dU0[7] +
              dU0[0] * U0[1] + dU0[3] * U0[4] + dU0[6] * U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0] * dU1[0] + U0[3] * dU1[3] + U0[6] * dU1[6] +
                       dU0[0] * U1[0] + dU0[3] * U1[3] + dU0[6] * U1[6]);
      b[4] = dU1[4] + (U0[1] * dU1[1] + U0[4] * dU1[4] + U0[7] * dU1[7] +
                       dU0[1] * U1[1] + dU0[4] * U1[4] + dU0[7] * U1[7]);
      b[5] =
          dU1[1] + dU1[3] +
          (U0[0] * dU1[1] + U0[3] * dU1[4] + U0[6] * dU1[7] + U1[0] * dU0[1] +
           U1[3] * dU0[4] + U1[6] * dU0[7] + dU0[0] * U1[1] + dU0[3] * U1[4] +
           dU0[6] * U1[7] + dU1[0] * U0[1] + dU1[3] * U0[4] + dU1[6] * U0[7]);
      b[6] = b[7] = 0.0;
      b += 8;
      dirdq += 3;
    }

    // Zero the contribution from the multiplier
    b[0] = b[1] = b[2] = b[3] = b[4] = b[5] = b[6] = b[7] = 0.0;
  }
}

/*
  Evaluate the derivative of the quantity eSens^{T}*B*psi with respect
  to the finite-element node locations.

  input:
  eSens:   input sensitivity vector
  psi:     input adjoint variables
  N:       shape functions
  Na, Nb:  derivative of the shape functions
  Ur:      derivative of displacements w.r.t. shell coordinates
  dr:      derivative of the direction w.r.t. shell coordinates
  Xdinv:   the inverse of the Jacobian matrix
  zXdinv:  the derivative of the inverse of the Jacobian w.r.t. z
  dirdq:   the derivative of the director field w.r.t. states

  output:
  Urd:     derivative of the function w.r.t. Ur
  drd:     derivative of the function w.r.t. dr
  Xdinvd:  derivative of the function w.r.t. Xdinv
  zXdinvd: derivative of the function w.r.t. Xdinvd
  dirdqd:  derivative of the function w.r.t. dirdq
*/
void MITC9::addBmatSens(TacsScalar Urd[], TacsScalar drd[], TacsScalar Xdinvd[],
                        TacsScalar zXdinvd[], TacsScalar Td[],
                        TacsScalar dirdqd[], const TacsScalar eSens[],
                        const TacsScalar psi[], const double N[],
                        const double Na[], const double Nb[],
                        const TacsScalar Ur[], const TacsScalar dr[],
                        const TacsScalar Xdinv[], const TacsScalar zXdinv[],
                        const TacsScalar T[], const TacsScalar dirdq[]) {
  // Compute U0 = T^{T}*Ur*Xdinv*T
  TacsScalar UrXdinv[9], UrXdinvT[9];
  matMatMult(Ur, Xdinv, UrXdinv);
  matMatMult(UrXdinv, T, UrXdinvT);

  TacsScalar U0[9];
  matTransMatMult(T, UrXdinvT, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  TacsScalar tsum[9];
  matMatMult(Ur, zXdinv, tsum);
  matMatMultAdd(dr, Xdinv, tsum);

  TacsScalar U1[9], drXdinvT[9];
  matMatMult(tsum, T, drXdinvT);
  matTransMatMult(T, drXdinvT, U1);

  // Zero internal terms that will be accumulated
  TacsScalar U0d[9], U1d[9];
  memset(U0d, 0, 9 * sizeof(TacsScalar));
  memset(U1d, 0, 9 * sizeof(TacsScalar));

  // Compute the derivatives of the strain w.r.t. the displacement
  // variables. This code takes advantage of the sparsity of the
  // derivatives to simplify the computations
  for (int i = 0; i < NUM_NODES; i++) {
    // Pull out the psi vector
    const TacsScalar *p = &psi[8 * i];

    // Compute dU0 = d(Ur)/dq_{k} * Xdinv
    // [  0,  0, 0 ][ 0  1  2 ][ T0  T1  T2 ]
    // [ Na, Nb, 0 ][ 3  4  5 ][ T3  T4  T5 ]
    // [  0,  0, 0 ][ 6  7  8 ][ T6  T7  T8 ]
    TacsScalar t0[3], tmp[3], tmpd[3];
    t0[0] = Na[i] * Xdinv[0] + Nb[i] * Xdinv[3];
    t0[1] = Na[i] * Xdinv[1] + Nb[i] * Xdinv[4];
    t0[2] = Na[i] * Xdinv[2] + Nb[i] * Xdinv[5];
    matMultTrans(T, t0, tmp);

    TacsScalar t1[3], ztmp[3], ztmpd[3];
    t1[0] = Na[i] * zXdinv[0] + Nb[i] * zXdinv[3];
    t1[1] = Na[i] * zXdinv[1] + Nb[i] * zXdinv[4];
    t1[2] = Na[i] * zXdinv[2] + Nb[i] * zXdinv[5];
    matMultTrans(T, t1, ztmp);

    // Zero the seeds
    tmpd[0] = tmpd[1] = tmpd[2] = 0.0;
    ztmpd[0] = ztmpd[1] = ztmpd[2] = 0.0;

    for (int k = 0; k < 3; k++) {
      // dU0[3*i+j] = T[3*k+i]*tmp[j];
      TacsScalar dU0[9];
      dU0[0] = T[3 * k] * tmp[0];
      dU0[1] = T[3 * k] * tmp[1];
      dU0[2] = T[3 * k] * tmp[2];
      dU0[3] = T[3 * k + 1] * tmp[0];
      dU0[4] = T[3 * k + 1] * tmp[1];
      dU0[5] = T[3 * k + 1] * tmp[2];
      dU0[6] = T[3 * k + 2] * tmp[0];
      dU0[7] = T[3 * k + 2] * tmp[1];
      dU0[8] = T[3 * k + 2] * tmp[2];

      // dU1[3*i+j] = T[3*k+i]*ztmp[j];
      TacsScalar dU1[9];
      dU1[0] = T[3 * k] * ztmp[0];
      dU1[1] = T[3 * k] * ztmp[1];
      dU1[2] = T[3 * k] * ztmp[2];
      dU1[3] = T[3 * k + 1] * ztmp[0];
      dU1[4] = T[3 * k + 1] * ztmp[1];
      dU1[5] = T[3 * k + 1] * ztmp[2];
      dU1[6] = T[3 * k + 2] * ztmp[0];
      dU1[7] = T[3 * k + 2] * ztmp[1];
      dU1[8] = T[3 * k + 2] * ztmp[2];

      // Compute the derivative of the function w.r.t. dU0
      TacsScalar dU0d[9];
      dU0d[0] = p[0] * (eSens[0] * (1.0 + U0[0]) + eSens[2] * U0[1] +
                        eSens[3] * U1[0] + eSens[5] * U1[1]);
      dU0d[1] = p[0] * (eSens[2] * (1.0 + U0[0]) + eSens[1] * U0[1] +
                        eSens[4] * U1[1] + eSens[5] * U1[0]);
      dU0d[2] = 0.0;
      dU0d[3] = p[0] * (eSens[2] * (1.0 + U0[4]) + eSens[0] * U0[3] +
                        eSens[3] * U1[3] + eSens[5] * U1[4]);
      dU0d[4] = p[0] * (eSens[1] * (1.0 + U0[4]) + eSens[2] * U0[3] +
                        eSens[4] * U1[4] + eSens[5] * U1[3]);
      dU0d[5] = 0.0;
      dU0d[6] = p[0] * (eSens[0] * U0[6] + eSens[2] * U0[7] + eSens[3] * U1[6] +
                        eSens[5] * U1[7]);
      dU0d[7] = p[0] * (eSens[1] * U0[7] + eSens[2] * U0[6] + eSens[4] * U1[7] +
                        eSens[5] * U1[6]);
      dU0d[8] = 0.0;

      // Compute the derivative with respect to dU1
      TacsScalar dU1d[9];
      dU1d[0] = p[0] * (eSens[3] * (1.0 + U0[0]) + eSens[5] * U0[1]);
      dU1d[1] = p[0] * (eSens[5] * (1.0 + U0[0]) + eSens[4] * U0[1]);
      dU1d[2] = 0.0;

      dU1d[3] = p[0] * (eSens[5] * (1.0 + U0[4]) + eSens[3] * U0[3]);
      dU1d[4] = p[0] * (eSens[4] * (1.0 + U0[4]) + eSens[5] * U0[3]);
      dU1d[5] = 0.0;

      dU1d[6] = p[0] * (eSens[3] * U0[6] + eSens[5] * U0[7]);
      dU1d[7] = p[0] * (eSens[4] * U0[7] + eSens[5] * U0[6]);
      dU1d[8] = 0.0;

      // Add the contributions to U0d
      U0d[0] += p[0] * (eSens[0] * dU0[0] + eSens[2] * dU0[1] +
                        eSens[3] * dU1[0] + eSens[5] * dU1[1]);
      U0d[1] += p[0] * (eSens[2] * dU0[0] + eSens[1] * dU0[1] +
                        eSens[4] * dU1[1] + eSens[5] * dU1[0]);
      U0d[3] += p[0] * (eSens[2] * dU0[4] + eSens[0] * dU0[3] +
                        eSens[3] * dU1[3] + eSens[5] * dU1[4]);
      U0d[4] += p[0] * (eSens[1] * dU0[4] + eSens[2] * dU0[3] +
                        eSens[4] * dU1[4] + eSens[5] * dU1[3]);
      U0d[6] += p[0] * (eSens[0] * dU0[6] + eSens[2] * dU0[7] +
                        eSens[3] * dU1[6] + eSens[5] * dU1[7]);
      U0d[7] += p[0] * (eSens[1] * dU0[7] + eSens[2] * dU0[6] +
                        eSens[4] * dU1[7] + eSens[5] * dU1[6]);

      // Add the contributions to U1d
      U1d[0] += p[0] * (eSens[3] * dU0[0] + eSens[5] * dU0[1]);
      U1d[1] += p[0] * (eSens[5] * dU0[0] + eSens[4] * dU0[1]);
      U1d[3] += p[0] * (eSens[5] * dU0[4] + eSens[3] * dU0[3]);
      U1d[4] += p[0] * (eSens[4] * dU0[4] + eSens[5] * dU0[3]);
      U1d[6] += p[0] * (eSens[3] * dU0[6] + eSens[5] * dU0[7]);
      U1d[7] += p[0] * (eSens[4] * dU0[7] + eSens[5] * dU0[6]);

      // dU0[3*i+j] = T[3*k+i]*tmp[j];
      // dU1[3*i+j] = T[3*k+i]*ztmp[j];
      Td[3 * k] += (dU0d[0] * tmp[0] + dU1d[0] * ztmp[0] + dU0d[1] * tmp[1] +
                    dU1d[1] * ztmp[1] + dU0d[2] * tmp[2] + dU1d[2] * ztmp[2]);
      Td[3 * k + 1] +=
          (dU0d[3] * tmp[0] + dU1d[3] * ztmp[0] + dU0d[4] * tmp[1] +
           dU1d[4] * ztmp[1] + dU0d[5] * tmp[2] + dU1d[5] * ztmp[2]);
      Td[3 * k + 2] +=
          (dU0d[6] * tmp[0] + dU1d[6] * ztmp[0] + dU0d[7] * tmp[1] +
           dU1d[7] * ztmp[1] + dU0d[8] * tmp[2] + dU1d[8] * ztmp[2]);

      tmpd[0] += (dU0d[0] * T[3 * k] + dU0d[3] * T[3 * k + 1] +
                  dU0d[6] * T[3 * k + 2]);
      tmpd[1] += (dU0d[1] * T[3 * k] + dU0d[4] * T[3 * k + 1] +
                  dU0d[7] * T[3 * k + 2]);
      tmpd[2] += (dU0d[2] * T[3 * k] + dU0d[5] * T[3 * k + 1] +
                  dU0d[8] * T[3 * k + 2]);

      ztmpd[0] += (dU1d[0] * T[3 * k] + dU1d[3] * T[3 * k + 1] +
                   dU1d[6] * T[3 * k + 2]);
      ztmpd[1] += (dU1d[1] * T[3 * k] + dU1d[4] * T[3 * k + 1] +
                   dU1d[7] * T[3 * k + 2]);
      ztmpd[2] += (dU1d[2] * T[3 * k] + dU1d[5] * T[3 * k + 1] +
                   dU1d[8] * T[3 * k + 2]);

      // Increment the pointer to the input vector
      p++;
    }

    TacsScalar t0d[3];
    // tmp_{i} = T_{ji}*t0_{j}
    matMult(T, tmpd, t0d);
    Xdinvd[0] += t0d[0] * Na[i];
    Xdinvd[1] += t0d[1] * Na[i];
    Xdinvd[2] += t0d[2] * Na[i];

    Xdinvd[3] += t0d[0] * Nb[i];
    Xdinvd[4] += t0d[1] * Nb[i];
    Xdinvd[5] += t0d[2] * Nb[i];

    TacsScalar t1d[3];
    matMult(T, ztmpd, t1d);
    zXdinvd[0] += t1d[0] * Na[i];
    zXdinvd[1] += t1d[1] * Na[i];
    zXdinvd[2] += t1d[2] * Na[i];

    zXdinvd[3] += t1d[0] * Nb[i];
    zXdinvd[4] += t1d[1] * Nb[i];
    zXdinvd[5] += t1d[2] * Nb[i];

    // Add the contribution to the Td matrix
    Td[0] += t0[0] * tmpd[0] + t1[0] * ztmpd[0];
    Td[1] += t0[0] * tmpd[1] + t1[0] * ztmpd[1];
    Td[2] += t0[0] * tmpd[2] + t1[0] * ztmpd[2];

    Td[3] += t0[1] * tmpd[0] + t1[1] * ztmpd[0];
    Td[4] += t0[1] * tmpd[1] + t1[1] * ztmpd[1];
    Td[5] += t0[1] * tmpd[2] + t1[1] * ztmpd[2];

    Td[6] += t0[2] * tmpd[0] + t1[2] * ztmpd[0];
    Td[7] += t0[2] * tmpd[1] + t1[2] * ztmpd[1];
    Td[8] += t0[2] * tmpd[2] + t1[2] * ztmpd[2];
  }

  // Add the contributions from the derivative of the strain w.r.t. the
  // rotation variables. These derivatives only make contributions to the
  // bending strains, not the in-plane strains.
  for (int i = 0; i < NUM_NODES; i++) {
    const TacsScalar *p = &psi[8 * i + 3];

    for (int k = 0; k < 4; k++) {
      // Compute dU0 = T^{T}*dUr*Xdinv*T
      // T^{T}*[ 0 | 0 | N*dirdq[0] ]*Xdinv*T
      // .     [ 0 | 0 | N*dirdq[1] ]*Xdinv*T
      //       [ 0 | 0 | N*dirdq[2] ]*Xdinv*T
      TacsScalar dU0[9], tu0[9], tu0T[9];
      tu0[0] = N[i] * dirdq[0] * Xdinv[6];
      tu0[1] = N[i] * dirdq[0] * Xdinv[7];
      tu0[2] = N[i] * dirdq[0] * Xdinv[8];
      tu0[3] = N[i] * dirdq[1] * Xdinv[6];
      tu0[4] = N[i] * dirdq[1] * Xdinv[7];
      tu0[5] = N[i] * dirdq[1] * Xdinv[8];
      tu0[6] = N[i] * dirdq[2] * Xdinv[6];
      tu0[7] = N[i] * dirdq[2] * Xdinv[7];
      tu0[8] = N[i] * dirdq[2] * Xdinv[8];
      matMatMult(tu0, T, tu0T);
      matTransMatMult(T, tu0T, dU0);

      // Compute the derivative for d
      TacsScalar drdq[9];
      drdq[0] = Na[i] * dirdq[0];
      drdq[3] = Na[i] * dirdq[1];
      drdq[6] = Na[i] * dirdq[2];
      drdq[1] = Nb[i] * dirdq[0];
      drdq[4] = Nb[i] * dirdq[1];
      drdq[7] = Nb[i] * dirdq[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;

      // Compute dU1
      TacsScalar dU1[9], tu1[9], tu1T[9];
      matMatMult(drdq, Xdinv, tu1);
      tu1[0] += N[i] * dirdq[0] * zXdinv[6];
      tu1[1] += N[i] * dirdq[0] * zXdinv[7];
      tu1[2] += N[i] * dirdq[0] * zXdinv[8];
      tu1[3] += N[i] * dirdq[1] * zXdinv[6];
      tu1[4] += N[i] * dirdq[1] * zXdinv[7];
      tu1[5] += N[i] * dirdq[1] * zXdinv[8];
      tu1[6] += N[i] * dirdq[2] * zXdinv[6];
      tu1[7] += N[i] * dirdq[2] * zXdinv[7];
      tu1[8] += N[i] * dirdq[2] * zXdinv[8];
      matMatMult(tu1, T, tu1T);
      matTransMatMult(T, tu1T, dU1);

      // Compute the derivative of the function w.r.t. dU0
      TacsScalar dU0d[9];
      dU0d[0] = p[0] * (eSens[0] * (1.0 + U0[0]) + eSens[2] * U0[1] +
                        eSens[3] * U1[0] + eSens[5] * U1[1]);
      dU0d[1] = p[0] * (eSens[2] * (1.0 + U0[0]) + eSens[1] * U0[1] +
                        eSens[4] * U1[1] + eSens[5] * U1[0]);
      dU0d[2] = 0.0;
      dU0d[3] = p[0] * (eSens[2] * (1.0 + U0[4]) + eSens[0] * U0[3] +
                        eSens[3] * U1[3] + eSens[5] * U1[4]);
      dU0d[4] = p[0] * (eSens[1] * (1.0 + U0[4]) + eSens[2] * U0[3] +
                        eSens[4] * U1[4] + eSens[5] * U1[3]);
      dU0d[5] = 0.0;
      dU0d[6] = p[0] * (eSens[0] * U0[6] + eSens[2] * U0[7] + eSens[3] * U1[6] +
                        eSens[5] * U1[7]);
      dU0d[7] = p[0] * (eSens[1] * U0[7] + eSens[2] * U0[6] + eSens[4] * U1[7] +
                        eSens[5] * U1[6]);
      dU0d[8] = 0.0;

      // Compute the derivative with respect to dU1
      TacsScalar dU1d[9];
      dU1d[0] = p[0] * (eSens[3] * (1.0 + U0[0]) + eSens[5] * U0[1]);
      dU1d[1] = p[0] * (eSens[5] * (1.0 + U0[0]) + eSens[4] * U0[1]);
      dU1d[2] = 0.0;

      dU1d[3] = p[0] * (eSens[5] * (1.0 + U0[4]) + eSens[3] * U0[3]);
      dU1d[4] = p[0] * (eSens[4] * (1.0 + U0[4]) + eSens[5] * U0[3]);
      dU1d[5] = 0.0;

      dU1d[6] = p[0] * (eSens[3] * U0[6] + eSens[5] * U0[7]);
      dU1d[7] = p[0] * (eSens[4] * U0[7] + eSens[5] * U0[6]);
      dU1d[8] = 0.0;

      // Add the contributions to U0d
      U0d[0] += p[0] * (eSens[0] * dU0[0] + eSens[2] * dU0[1] +
                        eSens[3] * dU1[0] + eSens[5] * dU1[1]);
      U0d[1] += p[0] * (eSens[2] * dU0[0] + eSens[1] * dU0[1] +
                        eSens[4] * dU1[1] + eSens[5] * dU1[0]);
      U0d[3] += p[0] * (eSens[2] * dU0[4] + eSens[0] * dU0[3] +
                        eSens[3] * dU1[3] + eSens[5] * dU1[4]);
      U0d[4] += p[0] * (eSens[1] * dU0[4] + eSens[2] * dU0[3] +
                        eSens[4] * dU1[4] + eSens[5] * dU1[3]);
      U0d[6] += p[0] * (eSens[0] * dU0[6] + eSens[2] * dU0[7] +
                        eSens[3] * dU1[6] + eSens[5] * dU1[7]);
      U0d[7] += p[0] * (eSens[1] * dU0[7] + eSens[2] * dU0[6] +
                        eSens[4] * dU1[7] + eSens[5] * dU1[6]);

      // Add the contributions to U1d
      U1d[0] += p[0] * (eSens[3] * dU0[0] + eSens[5] * dU0[1]);
      U1d[1] += p[0] * (eSens[5] * dU0[0] + eSens[4] * dU0[1]);
      U1d[3] += p[0] * (eSens[5] * dU0[4] + eSens[3] * dU0[3]);
      U1d[4] += p[0] * (eSens[4] * dU0[4] + eSens[5] * dU0[3]);
      U1d[6] += p[0] * (eSens[3] * dU0[6] + eSens[5] * dU0[7]);
      U1d[7] += p[0] * (eSens[4] * dU0[7] + eSens[5] * dU0[6]);

      // dU0 = T^{T}*tu*T
      // dU0_{kl} = T_{nk}*tu0_{nm}*T_{ml}
      // Add the contribution from the term: dU0 = T^{T}*tu0*T
      TacsScalar t[9];
      matMatTransMultAdd(tu0T, dU0d, Td);
      matTransMatMult(tu0, T, t);
      matMatMultAdd(t, dU0d, Td);

      // Add the contribution from the term: dU1 = T^{T}*tu1*T
      matMatTransMultAdd(tu1T, dU1d, Td);
      matTransMatMult(tu1, T, t);
      matMatMultAdd(t, dU1d, Td);

      // Add the contribution to zXdinvd
      TacsScalar tu1d[9];
      matMatTransMult(dU1d, T, t);
      matMatMult(T, t, tu1d);
      zXdinvd[6] +=
          N[i] * (dirdq[0] * tu1d[0] + dirdq[1] * tu1d[3] + dirdq[2] * tu1d[6]);
      zXdinvd[7] +=
          N[i] * (dirdq[0] * tu1d[1] + dirdq[1] * tu1d[4] + dirdq[2] * tu1d[7]);
      zXdinvd[8] +=
          N[i] * (dirdq[0] * tu1d[2] + dirdq[1] * tu1d[5] + dirdq[2] * tu1d[8]);

      // Add the contribution to Xdinvd
      TacsScalar tu0d[9];
      matMatTransMult(dU0d, T, t);
      matMatMult(T, t, tu0d);
      Xdinvd[6] +=
          N[i] * (dirdq[0] * tu0d[0] + dirdq[1] * tu0d[3] + dirdq[2] * tu0d[6]);
      Xdinvd[7] +=
          N[i] * (dirdq[0] * tu0d[1] + dirdq[1] * tu0d[4] + dirdq[2] * tu0d[7]);
      Xdinvd[8] +=
          N[i] * (dirdq[0] * tu0d[2] + dirdq[1] * tu0d[5] + dirdq[2] * tu0d[8]);

      // Add the sensitivity from the the multiplication
      // Given that: tu1_{kl} = drdq_{kn}*Xdinv_{nl}
      // Then the derivative is
      // d(tu1_{kl})/d{Xdvin_{ij}} = drdq_{ki}*delta_{lj}
      // Xdinvd_{ij}
      // = df/d{tu1_{kl}}x{drdq_{ki}*delta_{lj}}
      // = tu1d_{kj}*drdq_{ki}
      matTransMatMultAdd(drdq, tu1d, Xdinvd);

      // Add the contributions to dirdqd
      dirdqd[0] +=
          N[i] * (Xdinv[6] * tu0d[0] + Xdinv[7] * tu0d[1] + Xdinv[8] * tu0d[2]);
      dirdqd[1] +=
          N[i] * (Xdinv[6] * tu0d[3] + Xdinv[7] * tu0d[4] + Xdinv[8] * tu0d[5]);
      dirdqd[2] +=
          N[i] * (Xdinv[6] * tu0d[6] + Xdinv[7] * tu0d[7] + Xdinv[8] * tu0d[8]);

      dirdqd[0] += N[i] * (zXdinv[6] * tu1d[0] + zXdinv[7] * tu1d[1] +
                           zXdinv[8] * tu1d[2]);
      dirdqd[1] += N[i] * (zXdinv[6] * tu1d[3] + zXdinv[7] * tu1d[4] +
                           zXdinv[8] * tu1d[5]);
      dirdqd[2] += N[i] * (zXdinv[6] * tu1d[6] + zXdinv[7] * tu1d[7] +
                           zXdinv[8] * tu1d[8]);

      // Add the derivative contribution from drdq
      // tu1_{kl} = drdq_{kn}*Xdinv_{nl}
      // tu1_{kl}/d(drdq_{ij}) = delta_{ik}*Xdinv{jl}
      // dtu1_{kl}*delta_{ik}*Xdinv{jl} = dtu1_{il}*Xdinv_{jl}
      TacsScalar drdqd[9];
      matMatTransMult(tu1d, Xdinv, drdqd);
      dirdqd[0] += Na[i] * drdqd[0] + Nb[i] * drdqd[1];
      dirdqd[1] += Na[i] * drdqd[3] + Nb[i] * drdqd[4];
      dirdqd[2] += Na[i] * drdqd[6] + Nb[i] * drdqd[7];

      // Increment the pointers to the data
      dirdqd += 3;
      dirdq += 3;
      p++;
    }
  }

  // Compute the derivative w.r.t. Ur
  TacsScalar t[9];    // Temporary matrix
  TacsScalar t0d[9];  // = T*U0d*T^{T}
  matMatTransMult(U0d, T, t);
  matMatMult(T, t, t0d);
  matMatTransMultAdd(t0d, Xdinv, Urd);

  TacsScalar t1d[9];  // = T*U1d*T^{T}
  matMatTransMult(U1d, T, t);
  matMatMult(T, t, t1d);
  matMatTransMultAdd(t1d, zXdinv, Urd);

  // Compute the contribution to the derivative w.r.t. dr
  matMatTransMultAdd(t1d, Xdinv, drd);

  // Compute the contribution to the derivative w.r.t. Xdinv
  matTransMatMultAdd(dr, t1d, Xdinvd);
  matTransMatMultAdd(Ur, t0d, Xdinvd);

  // Compute the contribution to the derivative w.r.t. zXdinv
  matTransMatMultAdd(Ur, t1d, zXdinvd);

  // Compute the contribution to the derivative w.r.t. T
  matMatTransMultAdd(UrXdinvT, U0d, Td);  // Td = Ur*Xdinv*T*U0d^{T}
  matTransMatMult(UrXdinv, T, t);
  matMatMultAdd(t, U0d, Td);  // Td += (UrXdinv)^{T}*T*Ud0

  matMatTransMultAdd(drXdinvT, U1d, Td);  // Td = Ur*Xdinv*T*U0d^{T}
  matTransMatMult(tsum, T, t);
  matMatMultAdd(t, U1d, Td);  // Td += (UrXdinv)^{T}*T*Ud0
}

/*
  Add the contribution from the geometric stiffness
*/
void MITC9::addGmat(TacsScalar J[], const TacsScalar scale,
                    const TacsScalar s[], const double N[], const double Na[],
                    const double Nb[], const TacsScalar Ur[],
                    const TacsScalar dr[], const TacsScalar Xdinv[],
                    const TacsScalar zXdinv[], const TacsScalar T[],
                    const TacsScalar Xr[], const TacsScalar dirdq[]) {
  // The gradient of the displacement field
  TacsScalar U0[9], U1[9], tmp[9];

  // Store the first derivatives w.r.t. the displacements
  // and rotations in the element
  TacsScalar dU0[7 * 9 * NUM_NODES], dU1[7 * 9 * NUM_NODES];

  // Compute U0 = T^{T}*Ur*Xdinv*T
  matMatMult(Ur, Xdinv, U0);
  matMatMult(U0, T, tmp);
  matTransMatMult(T, tmp, U0);

  // Compute U1 = T^{T}*(dr*Xdinv + Ur*zXdinv)*T
  matMatMult(Ur, zXdinv, U1);
  matMatMultAdd(dr, Xdinv, U1);
  matMatMult(U1, T, tmp);
  matTransMatMult(T, tmp, U1);

  // Pre-compute terms that are required for the geometric
  // stiffness matrix computation
  TacsScalar *du0 = dU0, *du1 = dU1;
  for (int i = 0; i < NUM_NODES; i++) {
    TacsScalar t[3], ztmp[3];
    t[0] = Na[i] * Xdinv[0] + Nb[i] * Xdinv[3];
    t[1] = Na[i] * Xdinv[1] + Nb[i] * Xdinv[4];
    t[2] = Na[i] * Xdinv[2] + Nb[i] * Xdinv[5];
    matMultTrans(T, t, tmp);

    t[0] = Na[i] * zXdinv[0] + Nb[i] * zXdinv[3];
    t[1] = Na[i] * zXdinv[1] + Nb[i] * zXdinv[4];
    t[2] = Na[i] * zXdinv[2] + Nb[i] * zXdinv[5];
    matMultTrans(T, t, ztmp);

    for (int k = 0; k < 3; k++) {
      // dU0 = T[3*k+i]*tmp[3*k+j];
      du0[0] = T[3 * k] * tmp[0];
      du0[1] = T[3 * k] * tmp[1];
      du0[2] = T[3 * k] * tmp[2];
      du0[3] = T[3 * k + 1] * tmp[0];
      du0[4] = T[3 * k + 1] * tmp[1];
      du0[5] = T[3 * k + 1] * tmp[2];
      du0[6] = T[3 * k + 2] * tmp[0];
      du0[7] = T[3 * k + 2] * tmp[1];
      du0[8] = T[3 * k + 2] * tmp[2];

      // dU1 = T[3*k+i]*ztmp[3*k+j];
      du1[0] = T[3 * k] * ztmp[0];
      du1[1] = T[3 * k] * ztmp[1];
      du1[2] = T[3 * k] * ztmp[2];
      du1[3] = T[3 * k + 1] * ztmp[0];
      du1[4] = T[3 * k + 1] * ztmp[1];
      du1[5] = T[3 * k + 1] * ztmp[2];
      du1[6] = T[3 * k + 2] * ztmp[0];
      du1[7] = T[3 * k + 2] * ztmp[1];
      du1[8] = T[3 * k + 2] * ztmp[2];

      du0 += 9;
      du1 += 9;
    }

    for (int k = 0; k < 4; k++) {
      // Compute du0 = T^{T}*dur*Xdinv*T
      du0[0] = N[i] * dirdq[0] * Xdinv[6];
      du0[1] = N[i] * dirdq[0] * Xdinv[7];
      du0[2] = N[i] * dirdq[0] * Xdinv[8];
      du0[3] = N[i] * dirdq[1] * Xdinv[6];
      du0[4] = N[i] * dirdq[1] * Xdinv[7];
      du0[5] = N[i] * dirdq[1] * Xdinv[8];
      du0[6] = N[i] * dirdq[2] * Xdinv[6];
      du0[7] = N[i] * dirdq[2] * Xdinv[7];
      du0[8] = N[i] * dirdq[2] * Xdinv[8];
      matMatMult(du0, T, tmp);
      matTransMatMult(T, tmp, du0);

      // Add the contributions from the other
      // derivatives
      TacsScalar drdq[9];
      drdq[0] = Na[i] * dirdq[0];
      drdq[3] = Na[i] * dirdq[1];
      drdq[6] = Na[i] * dirdq[2];
      drdq[1] = Nb[i] * dirdq[0];
      drdq[4] = Nb[i] * dirdq[1];
      drdq[7] = Nb[i] * dirdq[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, du1);

      du1[0] += N[i] * dirdq[0] * zXdinv[6];
      du1[1] += N[i] * dirdq[0] * zXdinv[7];
      du1[2] += N[i] * dirdq[0] * zXdinv[8];
      du1[3] += N[i] * dirdq[1] * zXdinv[6];
      du1[4] += N[i] * dirdq[1] * zXdinv[7];
      du1[5] += N[i] * dirdq[1] * zXdinv[8];
      du1[6] += N[i] * dirdq[2] * zXdinv[6];
      du1[7] += N[i] * dirdq[2] * zXdinv[7];
      du1[8] += N[i] * dirdq[2] * zXdinv[8];
      matMatMult(du1, T, tmp);
      matTransMatMult(T, tmp, du1);

      du0 += 9;
      du1 += 9;
      dirdq += 3;
    }
  }

  // Compute the derivatives of the strain w.r.t. the displacement
  // variables. This code takes advantage of the sparsity of the
  // derivatives to simplify the computations
  for (int i = 0; i < 7 * NUM_NODES; i++) {
    const TacsScalar *dU0i = &dU0[9 * i];
    const TacsScalar *dU1i = &dU1[9 * i];
    for (int j = i; j < 7 * NUM_NODES; j++) {
      const TacsScalar *dU0j = &dU0[9 * j];
      const TacsScalar *dU1j = &dU1[9 * j];

      // Compute the real indices
      int ii = 8 * (i / 7) + (i % 7);
      int jj = 8 * (j / 7) + (j % 7);
      int idx = 8 * NUM_NODES * ii + jj;
      int sym = 8 * NUM_NODES * jj + ii;

      // Compute the derivative of the in-plane strain
      TacsScalar b[6];
      b[0] = dU0i[0] * dU0j[0] + dU0i[3] * dU0j[3] + dU0i[6] * dU0j[6];
      b[1] = dU0i[1] * dU0j[1] + dU0i[4] * dU0j[4] + dU0i[7] * dU0j[7];
      b[2] = (dU0i[0] * dU0j[1] + dU0i[3] * dU0j[4] + dU0i[6] * dU0j[7] +
              dU0j[0] * dU0i[1] + dU0j[3] * dU0i[4] + dU0j[6] * dU0i[7]);

      b[3] = (dU0i[0] * dU1j[0] + dU0i[3] * dU1j[3] + dU0i[6] * dU1j[6] +
              dU0j[0] * dU1i[0] + dU0j[3] * dU1i[3] + dU0j[6] * dU1i[6]);
      b[4] = (dU0i[1] * dU1j[1] + dU0i[4] * dU1j[4] + dU0i[7] * dU1j[7] +
              dU0j[1] * dU1i[1] + dU0j[4] * dU1i[4] + dU0j[7] * dU1i[7]);
      b[5] = (dU0i[0] * dU1j[1] + dU0i[3] * dU1j[4] + dU0i[6] * dU1j[7] +
              dU1i[0] * dU0j[1] + dU1i[3] * dU0j[4] + dU1i[6] * dU0j[7] +
              dU0j[0] * dU1i[1] + dU0j[3] * dU1i[4] + dU0j[6] * dU1i[7] +
              dU1j[0] * dU0i[1] + dU1j[3] * dU0i[4] + dU1j[6] * dU0i[7]);

      TacsScalar Jadd = scale * (b[0] * s[0] + b[1] * s[1] + b[2] * s[2] +
                                 b[3] * s[3] + b[4] * s[4] + b[5] * s[5]);

      // Add the values symmetrically
      J[idx] += Jadd;
      if (ii != jj) {
        J[sym] += Jadd;
      }
    }
  }

  // Add the contributions from the second derivatives of the
  // quaternions
  for (int i = 0; i < NUM_NODES; i++) {
    // d = N[i]*C(q_{i})^{T}*n - these terms only get added along
    // each diagonal in the matrix

    // Extract the normal from the frame
    TacsScalar normal[3];
    normal[0] = Xr[2];
    normal[1] = Xr[5];
    normal[2] = Xr[8];
    Xr += 9;

    // Compute the second derivative w.r.t. the quaternion
    TacsScalar dCtndq[3 * 9];
    computeQtr2ndDeriv(normal, dCtndq);
    const TacsScalar *dC = dCtndq;

    // Compute the partials derivatives w.r.t. eta,eps
    for (int ii = 0; ii < 9; ii++) {
      // Compute dU0 = T^{T}*dUr*Xdinv*T
      dU0[0] = N[i] * dC[0] * Xdinv[6];
      dU0[1] = N[i] * dC[0] * Xdinv[7];
      dU0[2] = N[i] * dC[0] * Xdinv[8];
      dU0[3] = N[i] * dC[1] * Xdinv[6];
      dU0[4] = N[i] * dC[1] * Xdinv[7];
      dU0[5] = N[i] * dC[1] * Xdinv[8];
      dU0[6] = N[i] * dC[2] * Xdinv[6];
      dU0[7] = N[i] * dC[2] * Xdinv[7];
      dU0[8] = N[i] * dC[2] * Xdinv[8];
      matMatMult(dU0, T, tmp);
      matTransMatMult(T, tmp, dU0);

      // Add the contributions from the other
      // derivatives
      TacsScalar drdq[9];
      drdq[0] = Na[i] * dC[0];
      drdq[3] = Na[i] * dC[1];
      drdq[6] = Na[i] * dC[2];
      drdq[1] = Nb[i] * dC[0];
      drdq[4] = Nb[i] * dC[1];
      drdq[7] = Nb[i] * dC[2];
      drdq[2] = drdq[5] = drdq[8] = 0.0;
      matMatMult(drdq, Xdinv, dU1);

      dU1[0] += N[i] * dC[0] * zXdinv[6];
      dU1[1] += N[i] * dC[0] * zXdinv[7];
      dU1[2] += N[i] * dC[0] * zXdinv[8];
      dU1[3] += N[i] * dC[1] * zXdinv[6];
      dU1[4] += N[i] * dC[1] * zXdinv[7];
      dU1[5] += N[i] * dC[1] * zXdinv[8];
      dU1[6] += N[i] * dC[2] * zXdinv[6];
      dU1[7] += N[i] * dC[2] * zXdinv[7];
      dU1[8] += N[i] * dC[2] * zXdinv[8];
      matMatMult(dU1, T, tmp);
      matTransMatMult(T, tmp, dU1);

      // Compute the derivative of the in-plane strain
      TacsScalar b[6];
      // Compute the derivative of the in-plane strain
      b[0] = dU0[0] + (U0[0] * dU0[0] + U0[3] * dU0[3] + U0[6] * dU0[6]);
      b[1] = dU0[4] + (U0[1] * dU0[1] + U0[4] * dU0[4] + U0[7] * dU0[7]);
      b[2] = dU0[1] + dU0[3] +
             (U0[0] * dU0[1] + U0[3] * dU0[4] + U0[6] * dU0[7] +
              dU0[0] * U0[1] + dU0[3] * U0[4] + dU0[6] * U0[7]);

      // Compute the derivative of the bending strain
      b[3] = dU1[0] + (U0[0] * dU1[0] + U0[3] * dU1[3] + U0[6] * dU1[6] +
                       dU0[0] * U1[0] + dU0[3] * U1[3] + dU0[6] * U1[6]);
      b[4] = dU1[4] + (U0[1] * dU1[1] + U0[4] * dU1[4] + U0[7] * dU1[7] +
                       dU0[1] * U1[1] + dU0[4] * U1[4] + dU0[7] * U1[7]);
      b[5] =
          dU1[1] + dU1[3] +
          (U0[0] * dU1[1] + U0[3] * dU1[4] + U0[6] * dU1[7] + U1[0] * dU0[1] +
           U1[3] * dU0[4] + U1[6] * dU0[7] + dU0[0] * U1[1] + dU0[3] * U1[4] +
           dU0[6] * U1[7] + dU1[0] * U0[1] + dU1[3] * U0[4] + dU1[6] * U0[7]);

      TacsScalar Jadd = scale * (b[0] * s[0] + b[1] * s[1] + b[2] * s[2] +
                                 b[3] * s[3] + b[4] * s[4] + b[5] * s[5]);

      if (ii < 3) {
        int iv = 8 * i + 3;
        int jv = 8 * i + 4 + ii;
        J[(8 * NUM_NODES) * iv + jv] += Jadd;
        J[(8 * NUM_NODES) * jv + iv] += Jadd;
      } else {
        if (ii == 3) {
          int iv = 8 * i + 4;
          J[(8 * NUM_NODES + 1) * iv] += Jadd;
        } else if (ii == 4) {
          int iv = 8 * i + 4, jv = 8 * i + 5;
          J[(8 * NUM_NODES) * iv + jv] += Jadd;
          J[(8 * NUM_NODES) * jv + iv] += Jadd;
        } else if (ii == 5) {
          int iv = 8 * i + 4, jv = 8 * i + 6;
          J[(8 * NUM_NODES) * iv + jv] += Jadd;
          J[(8 * NUM_NODES) * jv + iv] += Jadd;
        } else if (ii == 6) {
          int iv = 8 * i + 5;
          J[(8 * NUM_NODES + 1) * iv] += Jadd;
        } else if (ii == 7) {
          int iv = 8 * i + 5, jv = 8 * i + 6;
          J[(8 * NUM_NODES) * iv + jv] += Jadd;
          J[(8 * NUM_NODES) * jv + iv] += Jadd;
        } else if (ii == 8) {
          int iv = 8 * i + 6;
          J[(8 * NUM_NODES + 1) * iv] += Jadd;
        }
      }

      dC += 3;
    }
  }
}

/*
  Compute the value of the tensorial strain at the tying points
  in the element.

  This code evaluates the tensorial shear strain values at the
  tying points which consist of the 2-point Gauss quadrature
  points in one direction, and the nodal locations along the other
  direction.

  input:
  X:     the initial values of the nodal coordinates
  vars:  the values of the variables
  dr:    the director values at every node in the element

  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::computeTyingStrain(TacsScalar g13[], TacsScalar g23[],
                               const TacsScalar X[], const TacsScalar Xr[],
                               const TacsScalar vars[],
                               const TacsScalar dir[]) {
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t, t, -t, t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0, s, s};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    g13[k] = 0.5 * (vecDot(Xa, d) + vecDot(Ua, fn) + vecDot(Ua, d));
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0, s, -s, 0.0, s};
  const double g23_vpts[] = {-t, -t, -t, t, t, t};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    g23[k] = 0.5 * (vecDot(Xb, d) + vecDot(Ub, fn) + vecDot(Ub, d));
  }
}

/*
  Add the derivative of the tying strain computation to the output
  vectors Xd, Xrd and dird
*/
void MITC9::addComputeTyingStrainSens(
    TacsScalar Xd[], TacsScalar Xrd[], TacsScalar dird[],
    const TacsScalar g13d[], const TacsScalar g23d[], const TacsScalar X[],
    const TacsScalar Xr[], const TacsScalar vars[], const TacsScalar dir[]) {
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t, t, -t, t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0, s, s};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    // Extract the partials of g13
    TacsScalar Xad[3], fnd[3];
    TacsScalar scale = 0.5 * g13d[k];
    Xad[0] = scale * d[0];
    Xad[1] = scale * d[1];
    Xad[2] = scale * d[2];

    fnd[0] = scale * Ua[0];
    fnd[1] = scale * Ua[1];
    fnd[2] = scale * Ua[2];

    TacsScalar dd[3];
    dd[0] = scale * (Xa[0] + Ua[0]);
    dd[1] = scale * (Xa[1] + Ua[1]);
    dd[2] = scale * (Xa[2] + Ua[2]);

    // Add the normal frame component
    addFrameNormalSens(fnd, N, Xrd);

    // Add the contributions to dird and Xa
    for (int j = 0; j < NUM_NODES; j++) {
      Xd[3 * j] += Na[j] * Xad[0];
      Xd[1 + 3 * j] += Na[j] * Xad[1];
      Xd[2 + 3 * j] += Na[j] * Xad[2];

      dird[3 * j] += N[j] * dd[0];
      dird[1 + 3 * j] += N[j] * dd[1];
      dird[2 + 3 * j] += N[j] * dd[2];
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0, s, -s, 0.0, s};
  const double g23_vpts[] = {-t, -t, -t, t, t, t};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3], Xb[3];
    innerProduct(Na, X, Xa);
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    // Extract the partials of g13
    TacsScalar Xbd[3], fnd[3];
    TacsScalar scale = 0.5 * g23d[k];
    Xbd[0] = scale * d[0];
    Xbd[1] = scale * d[1];
    Xbd[2] = scale * d[2];

    fnd[0] = scale * Ub[0];
    fnd[1] = scale * Ub[1];
    fnd[2] = scale * Ub[2];

    TacsScalar dd[3];
    dd[0] = scale * (Xb[0] + Ub[0]);
    dd[1] = scale * (Xb[1] + Ub[1]);
    dd[2] = scale * (Xb[2] + Ub[2]);

    // Add the normal frame component
    addFrameNormalSens(fnd, N, Xrd);

    // Add the contributions to dird and Xa
    for (int j = 0; j < NUM_NODES; j++) {
      Xd[3 * j] += Nb[j] * Xbd[0];
      Xd[1 + 3 * j] += Nb[j] * Xbd[1];
      Xd[2 + 3 * j] += Nb[j] * Xbd[2];

      dird[3 * j] += N[j] * dd[0];
      dird[1 + 3 * j] += N[j] * dd[1];
      dird[2 + 3 * j] += N[j] * dd[2];
    }
  }
}

/*
  Compute the value of the tensorial strain and the derivative of the
  tensorial strain w.r.t. the state variables at the tying points in
  the element.

  This code evaluates the tensorial shear strain values at the tying
  points which consist of the 2-point Gauss quadrature points in one
  direction, and the nodal locations along the other direction.

  input:
  X:     the initial values of the nodal coordinates
  vars:  the values of the variables
  dr:    the director values at every node in the element

  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::computeTyingBmat(TacsScalar g13[], TacsScalar g23[],
                             TacsScalar B13[], TacsScalar B23[],
                             const TacsScalar X[], const TacsScalar Xr[],
                             const TacsScalar vars[], const TacsScalar dir[],
                             const TacsScalar dirdq[]) {
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t, t, -t, t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0, s, s};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    g13[k] = 0.5 * (vecDot(Xa, d) + vecDot(Ua, fn) + vecDot(Ua, d));

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xa[0] + Ua[0];
    t[1] = Xa[1] + Ua[1];
    t[2] = Xa[2] + Ua[2];

    TacsScalar *b13 = &B13[8 * NUM_NODES * k];
    const TacsScalar *dq = dirdq;
    for (int i = 0; i < NUM_NODES; i++) {
      // Compute the derivatives w.r.t. U
      b13[0] = 0.5 * Na[i] * (fn[0] + d[0]);
      b13[1] = 0.5 * Na[i] * (fn[1] + d[1]);
      b13[2] = 0.5 * Na[i] * (fn[2] + d[2]);

      // Compute the derivatives w.r.t. q
      b13[3] = 0.5 * N[i] * (dq[0] * t[0] + dq[1] * t[1] + dq[2] * t[2]);
      b13[4] = 0.5 * N[i] * (dq[3] * t[0] + dq[4] * t[1] + dq[5] * t[2]);
      b13[5] = 0.5 * N[i] * (dq[6] * t[0] + dq[7] * t[1] + dq[8] * t[2]);
      b13[6] = 0.5 * N[i] * (dq[9] * t[0] + dq[10] * t[1] + dq[11] * t[2]);
      b13[7] = 0.0;

      b13 += 8;
      dq += 12;
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0, s, -s, 0.0, s};
  const double g23_vpts[] = {-t, -t, -t, t, t, t};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    g23[k] = 0.5 * (vecDot(Xb, d) + vecDot(Ub, fn) + vecDot(Ub, d));

    TacsScalar t[3];
    t[0] = Xb[0] + Ub[0];
    t[1] = Xb[1] + Ub[1];
    t[2] = Xb[2] + Ub[2];

    TacsScalar *b23 = &B23[8 * NUM_NODES * k];
    const TacsScalar *dq = dirdq;
    for (int i = 0; i < NUM_NODES; i++) {
      // Compute the derivatives w.r.t. U
      b23[0] = 0.5 * Nb[i] * (fn[0] + d[0]);
      b23[1] = 0.5 * Nb[i] * (fn[1] + d[1]);
      b23[2] = 0.5 * Nb[i] * (fn[2] + d[2]);

      // Compute the derivatives w.r.t. q
      b23[3] = 0.5 * N[i] * (dq[0] * t[0] + dq[1] * t[1] + dq[2] * t[2]);
      b23[4] = 0.5 * N[i] * (dq[3] * t[0] + dq[4] * t[1] + dq[5] * t[2]);
      b23[5] = 0.5 * N[i] * (dq[6] * t[0] + dq[7] * t[1] + dq[8] * t[2]);
      b23[6] = 0.5 * N[i] * (dq[9] * t[0] + dq[10] * t[1] + dq[11] * t[2]);
      b23[7] = 0.0;

      b23 += 8;
      dq += 12;
    }
  }
}

/*
  Add the contributions from the derivative of the tying strain
  derivative matrix to the output vectors.
*/
void MITC9::addComputeTyingBmatSens(
    TacsScalar Xd[], TacsScalar Xrd[], TacsScalar dird[], TacsScalar dirdqd[],
    const TacsScalar g13d[], const TacsScalar g23d[], const TacsScalar B13d[],
    const TacsScalar B23d[], const TacsScalar X[], const TacsScalar Xr[],
    const TacsScalar vars[], const TacsScalar dir[], const TacsScalar dirdq[]) {
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t, t, -t, t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0, s, s};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    // Extract the partials of g13
    TacsScalar Xad[3], fnd[3];
    TacsScalar scale = 0.5 * g13d[k];
    Xad[0] = scale * d[0];
    Xad[1] = scale * d[1];
    Xad[2] = scale * d[2];

    fnd[0] = scale * Ua[0];
    fnd[1] = scale * Ua[1];
    fnd[2] = scale * Ua[2];

    TacsScalar dd[3];
    dd[0] = scale * (Xa[0] + Ua[0]);
    dd[1] = scale * (Xa[1] + Ua[1]);
    dd[2] = scale * (Xa[2] + Ua[2]);

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xa[0] + Ua[0];
    t[1] = Xa[1] + Ua[1];
    t[2] = Xa[2] + Ua[2];

    const TacsScalar *b13d = &B13d[8 * NUM_NODES * k];
    const TacsScalar *dq = dirdq;
    TacsScalar *dqd = dirdqd;
    for (int i = 0; i < NUM_NODES; i++) {
      // Add contributions from the derivatives w.r.t. U
      fnd[0] += 0.5 * Na[i] * b13d[0];
      fnd[1] += 0.5 * Na[i] * b13d[1];
      fnd[2] += 0.5 * Na[i] * b13d[2];

      dd[0] += 0.5 * Na[i] * b13d[0];
      dd[1] += 0.5 * Na[i] * b13d[1];
      dd[2] += 0.5 * Na[i] * b13d[2];

      // Add the contributions from the derivatives w.r.t. q
      Xad[0] += 0.5 * N[i] *
                (dq[0] * b13d[3] + dq[3] * b13d[4] + dq[6] * b13d[5] +
                 dq[9] * b13d[6]);
      Xad[1] += 0.5 * N[i] *
                (dq[1] * b13d[3] + dq[4] * b13d[4] + dq[7] * b13d[5] +
                 dq[10] * b13d[6]);
      Xad[2] += 0.5 * N[i] *
                (dq[2] * b13d[3] + dq[5] * b13d[4] + dq[8] * b13d[5] +
                 dq[11] * b13d[6]);

      dqd[0] += 0.5 * N[i] * t[0] * b13d[3];
      dqd[1] += 0.5 * N[i] * t[1] * b13d[3];
      dqd[2] += 0.5 * N[i] * t[2] * b13d[3];

      dqd[3] += 0.5 * N[i] * t[0] * b13d[4];
      dqd[4] += 0.5 * N[i] * t[1] * b13d[4];
      dqd[5] += 0.5 * N[i] * t[2] * b13d[4];

      dqd[6] += 0.5 * N[i] * t[0] * b13d[5];
      dqd[7] += 0.5 * N[i] * t[1] * b13d[5];
      dqd[8] += 0.5 * N[i] * t[2] * b13d[5];

      dqd[9] += 0.5 * N[i] * t[0] * b13d[6];
      dqd[10] += 0.5 * N[i] * t[1] * b13d[6];
      dqd[11] += 0.5 * N[i] * t[2] * b13d[6];

      b13d += 8;
      dq += 12;
      dqd += 12;
    }

    // Add the normal frame component
    addFrameNormalSens(fnd, N, Xrd);

    // Add the contributions to dird and Xa
    for (int j = 0; j < NUM_NODES; j++) {
      Xd[3 * j] += Na[j] * Xad[0];
      Xd[1 + 3 * j] += Na[j] * Xad[1];
      Xd[2 + 3 * j] += Na[j] * Xad[2];

      dird[3 * j] += N[j] * dd[0];
      dird[1 + 3 * j] += N[j] * dd[1];
      dird[2 + 3 * j] += N[j] * dd[2];
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0, s, -s, 0.0, s};
  const double g23_vpts[] = {-t, -t, -t, t, t, t};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[9], Na[9], Nb[9];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // Compute the frame normal at the current point
    TacsScalar fn[3];
    computeFrameNormal(N, Xr, fn);

    // Interpolate to compute the director at the current point
    TacsScalar d[3];
    innerProduct(N, dir, d);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    // Extract the partials of g13
    TacsScalar Xbd[3], fnd[3];
    TacsScalar scale = 0.5 * g23d[k];
    Xbd[0] = scale * d[0];
    Xbd[1] = scale * d[1];
    Xbd[2] = scale * d[2];

    fnd[0] = scale * Ub[0];
    fnd[1] = scale * Ub[1];
    fnd[2] = scale * Ub[2];

    TacsScalar dd[3];
    dd[0] = scale * (Xb[0] + Ub[0]);
    dd[1] = scale * (Xb[1] + Ub[1]);
    dd[2] = scale * (Xb[2] + Ub[2]);

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xb[0] + Ub[0];
    t[1] = Xb[1] + Ub[1];
    t[2] = Xb[2] + Ub[2];

    const TacsScalar *b23d = &B23d[8 * NUM_NODES * k];
    const TacsScalar *dq = dirdq;
    TacsScalar *dqd = dirdqd;
    for (int i = 0; i < NUM_NODES; i++) {
      // Add contributions from the derivatives w.r.t. U
      fnd[0] += 0.5 * Nb[i] * b23d[0];
      fnd[1] += 0.5 * Nb[i] * b23d[1];
      fnd[2] += 0.5 * Nb[i] * b23d[2];

      dd[0] += 0.5 * Nb[i] * b23d[0];
      dd[1] += 0.5 * Nb[i] * b23d[1];
      dd[2] += 0.5 * Nb[i] * b23d[2];

      // Add the contributions from the derivatives w.r.t. q
      Xbd[0] += 0.5 * N[i] *
                (dq[0] * b23d[3] + dq[3] * b23d[4] + dq[6] * b23d[5] +
                 dq[9] * b23d[6]);
      Xbd[1] += 0.5 * N[i] *
                (dq[1] * b23d[3] + dq[4] * b23d[4] + dq[7] * b23d[5] +
                 dq[10] * b23d[6]);
      Xbd[2] += 0.5 * N[i] *
                (dq[2] * b23d[3] + dq[5] * b23d[4] + dq[8] * b23d[5] +
                 dq[11] * b23d[6]);

      dqd[0] += 0.5 * N[i] * t[0] * b23d[3];
      dqd[1] += 0.5 * N[i] * t[1] * b23d[3];
      dqd[2] += 0.5 * N[i] * t[2] * b23d[3];

      dqd[3] += 0.5 * N[i] * t[0] * b23d[4];
      dqd[4] += 0.5 * N[i] * t[1] * b23d[4];
      dqd[5] += 0.5 * N[i] * t[2] * b23d[4];

      dqd[6] += 0.5 * N[i] * t[0] * b23d[5];
      dqd[7] += 0.5 * N[i] * t[1] * b23d[5];
      dqd[8] += 0.5 * N[i] * t[2] * b23d[5];

      dqd[9] += 0.5 * N[i] * t[0] * b23d[6];
      dqd[10] += 0.5 * N[i] * t[1] * b23d[6];
      dqd[11] += 0.5 * N[i] * t[2] * b23d[6];

      b23d += 8;
      dq += 12;
      dqd += 12;
    }

    // Add the normal frame component
    addFrameNormalSens(fnd, N, Xrd);

    // Add the contributions to dird and Xa
    for (int j = 0; j < NUM_NODES; j++) {
      Xd[3 * j] += Nb[j] * Xbd[0];
      Xd[1 + 3 * j] += Nb[j] * Xbd[1];
      Xd[2 + 3 * j] += Nb[j] * Xbd[2];

      dird[3 * j] += N[j] * dd[0];
      dird[1 + 3 * j] += N[j] * dd[1];
      dird[2 + 3 * j] += N[j] * dd[2];
    }
  }
}

/*
  Compute the value of the tensorial strain and the derivative of the
  tensorial strain w.r.t. the state variables at the tying points in
  the element.

  This code evaluates the tensorial shear strain values at the tying
  points which consist of the 2-point Gauss quadrature points in one
  direction, and the nodal locations along the other direction.

  input:
  X:     the initial values of the nodal coordinates
  vars:  the values of the variables
  dr:    the director values at every node in the element

  output:
  g13:   the values of the tensorial strain at the tying points
  g23:   the values of the tensorial strain at the tying points
*/
void MITC9::addTyingGmat(TacsScalar J[], const TacsScalar w13[],
                         const TacsScalar w23[], const TacsScalar X[],
                         const TacsScalar Xr[], const TacsScalar vars[],
                         const TacsScalar dir[], const TacsScalar dirdq[]) {
  const double s = 0.774596669241483;
  const double t = 0.577350269189626;

  const int iv[] = {3, 3, 3, 4, 4, 4, 5, 5, 6};
  const int jv[] = {4, 5, 6, 4, 5, 6, 5, 6, 6};

  // The tying points where the strain will be evaluated
  const double g13_upts[] = {-t, t, -t, t, -t, t};
  const double g13_vpts[] = {-s, -s, 0.0, 0.0, s, s};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    computeShapeFunc(g13_upts[k], g13_vpts[k], N);
    computeShapeFunc(g13_upts[k], g13_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xa[3];
    innerProduct(Na, X, Xa);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ua[3];
    innerProduct8(Na, vars, Ua);

    // Temp vector for computing the derivative
    TacsScalar t[3];
    t[0] = Xa[0] + Ua[0];
    t[1] = Xa[1] + Ua[1];
    t[2] = Xa[2] + Ua[2];

    const TacsScalar *xr = Xr;
    for (int i = 0; i < NUM_NODES; i++) {
      // Extract the normal from the frame
      TacsScalar normal[3];
      normal[0] = xr[2];
      normal[1] = xr[5];
      normal[2] = xr[8];
      xr += 9;

      // Compute the second derivative of the quaternion
      TacsScalar dCtndq[3 * 9];
      computeQtr2ndDeriv(normal, dCtndq);

      // Compute the second derivatives from the quaternions alone
      const TacsScalar *dC = dCtndq;
      for (int ii = 0; ii < 9; ii++) {
        TacsScalar Jadd = 0.5 * w13[k] * N[i] * vecDot(t, dC);
        J[(8 * NUM_NODES) * (8 * i + iv[ii]) + (8 * i + jv[ii])] += Jadd;
        if (iv[ii] != jv[ii]) {
          J[(8 * NUM_NODES) * (8 * i + jv[ii]) + (8 * i + iv[ii])] += Jadd;
        }
        dC += 3;
      }

      // Set the derivative of the director w.r.t. q
      const TacsScalar *dq = dirdq;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int ii = 0; ii < 4; ii++) {    // Loop over quaternions
          for (int jj = 0; jj < 3; jj++) {  // Loop over displacements
            TacsScalar Jadd = 0.5 * w13[k] * N[j] * Na[i] * dq[3 * ii + jj];
            J[8 * NUM_NODES * (8 * i + jj) + 8 * j + 3 + ii] += Jadd;
            J[8 * NUM_NODES * (8 * j + 3 + ii) + 8 * i + jj] += Jadd;
          }
        }
        dq += 12;
      }
    }
  }

  // The tying points where the strain will be evaluated
  const double g23_upts[] = {-s, 0.0, s, -s, 0.0, s};
  const double g23_vpts[] = {-t, -t, -t, t, t, t};

  for (int k = 0; k < 6; k++) {
    // Evaluate the element shape functions
    double N[NUM_NODES], Na[NUM_NODES], Nb[NUM_NODES];
    computeShapeFunc(g23_upts[k], g23_vpts[k], N);
    computeShapeFunc(g23_upts[k], g23_vpts[k], Na, Nb);

    // The derivative of the shell surface coordinates w.r.t. u and v
    // directions
    TacsScalar Xb[3];
    innerProduct(Nb, X, Xb);

    // The derivatives of the displacements w.r.t. the u and v
    // directions
    TacsScalar Ub[3];
    innerProduct8(Nb, vars, Ub);

    TacsScalar t[3];
    t[0] = Xb[0] + Ub[0];
    t[1] = Xb[1] + Ub[1];
    t[2] = Xb[2] + Ub[2];

    const TacsScalar *xr = Xr;
    for (int i = 0; i < NUM_NODES; i++) {
      // Extract the normal from the frame
      TacsScalar normal[3];
      normal[0] = xr[2];
      normal[1] = xr[5];
      normal[2] = xr[8];
      xr += 9;

      // Compute the second derivative of the quaternion
      TacsScalar dCtndq[3 * 9];
      computeQtr2ndDeriv(normal, dCtndq);

      // Compute the second derivatives from the quaternions alone
      const TacsScalar *dC = dCtndq;
      for (int ii = 0; ii < 9; ii++) {
        TacsScalar Jadd = 0.5 * w23[k] * N[i] * vecDot(t, dC);
        J[(8 * NUM_NODES) * (8 * i + iv[ii]) + (8 * i + jv[ii])] += Jadd;
        if (iv[ii] != jv[ii]) {
          J[(8 * NUM_NODES) * (8 * i + jv[ii]) + (8 * i + iv[ii])] += Jadd;
        }
        dC += 3;
      }

      // Set the derivative of the director w.r.t. q
      const TacsScalar *dq = dirdq;
      for (int j = 0; j < NUM_NODES; j++) {
        for (int ii = 0; ii < 4; ii++) {    // Loop over quaternions
          for (int jj = 0; jj < 3; jj++) {  // Loop over displacements
            TacsScalar Jadd = 0.5 * w23[k] * N[j] * Nb[i] * dq[3 * ii + jj];
            J[8 * NUM_NODES * (8 * i + jj) + 8 * j + 3 + ii] += Jadd;
            J[8 * NUM_NODES * (8 * j + 3 + ii) + 8 * i + jj] += Jadd;
          }
        }
        dq += 12;
      }
    }
  }
}

/*
  Add the assumed strain distribution interpolated from the tying
  points to the strain at the present point.

  This code interpolates the strain from the tying points, which
  interpolates the tensorial components of the strain, and adds these
  values to the strain. This involes an interpolation and then a
  transformation. Note that not all components of the transformation
  are required. In general, the transformation takes the form:

  e = A^{T}*g*A

  where g are the tensorial components of the strain and A is the
  transformation matrix (which is not necessarily orthonormal!)

  input:
  N13:  the shape functions for the g13 tensorial strain
  N23:  the shape functions for the g23 tensorial strain
  g13:  the g13 tensorial strain at the tying points
  g23:  the g23 tensorial strain at the tying points
  T:    the transformation matrix

  output:
  e:    the shear strain components are over-written
*/
void MITC9::addTyingStrain(TacsScalar e[], const double N13[],
                           const double N23[], const TacsScalar g13[],
                           const TacsScalar g23[], const TacsScalar Xdinv[],
                           const TacsScalar T[]) {
  // Compute the strain using the assumed strain distribution
  // and the strain evaluated at the tying points
  TacsScalar G13, G23;
  G13 = (g13[0] * N13[0] + g13[1] * N13[1] + g13[2] * N13[2] + g13[3] * N13[3] +
         g13[4] * N13[4] + g13[5] * N13[5]);
  G23 = (g23[0] * N23[0] + g23[1] * N23[1] + g23[2] * N23[2] + g23[3] * N23[3] +
         g23[4] * N23[4] + g23[5] * N23[5]);

  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0] * T[0] + Xdinv[1] * T[3] + Xdinv[2] * T[6];
  A12 = Xdinv[0] * T[1] + Xdinv[1] * T[4] + Xdinv[2] * T[7];

  A21 = Xdinv[3] * T[0] + Xdinv[4] * T[3] + Xdinv[5] * T[6];
  A22 = Xdinv[3] * T[1] + Xdinv[4] * T[4] + Xdinv[5] * T[7];

  // Compute and set the final strain values
  // e = 2*A^{T}*G
  e[6] = 2.0 * (A12 * G13 + A22 * G23);
  e[7] = 2.0 * (A11 * G13 + A21 * G23);
}

/*
  Add the derivative of the tying strain with respect to the inputs
  to the output data
*/
void MITC9::addTyingStrainSens(TacsScalar g13d[], TacsScalar g23d[],
                               TacsScalar Xdinvd[], TacsScalar Td[],
                               TacsScalar scale, const TacsScalar eSens[],
                               const double N13[], const double N23[],
                               const TacsScalar g13[], const TacsScalar g23[],
                               const TacsScalar Xdinv[], const TacsScalar T[]) {
  // Compute the strain using the assumed strain distribution
  // and the strain evaluated at the tying points
  TacsScalar G13, G23;
  G13 = (g13[0] * N13[0] + g13[1] * N13[1] + g13[2] * N13[2] + g13[3] * N13[3] +
         g13[4] * N13[4] + g13[5] * N13[5]);
  G23 = (g23[0] * N23[0] + g23[1] * N23[1] + g23[2] * N23[2] + g23[3] * N23[3] +
         g23[4] * N23[4] + g23[5] * N23[5]);

  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0] * T[0] + Xdinv[1] * T[3] + Xdinv[2] * T[6];
  A12 = Xdinv[0] * T[1] + Xdinv[1] * T[4] + Xdinv[2] * T[7];

  A21 = Xdinv[3] * T[0] + Xdinv[4] * T[3] + Xdinv[5] * T[6];
  A22 = Xdinv[3] * T[1] + Xdinv[4] * T[4] + Xdinv[5] * T[7];

  // Compute and set the final strain values
  // e = 2*A^{T}*G
  // e[6] = 2.0*(A12*G13 + A22*G23);
  // e[7] = 2.0*(A11*G13 + A21*G23);
  TacsScalar G13d = 2.0 * scale * (A12 * eSens[6] + A11 * eSens[7]);
  TacsScalar G23d = 2.0 * scale * (A22 * eSens[6] + A21 * eSens[7]);

  TacsScalar A11d = 2.0 * scale * G13 * eSens[7];
  TacsScalar A12d = 2.0 * scale * G13 * eSens[6];
  TacsScalar A21d = 2.0 * scale * G23 * eSens[7];
  TacsScalar A22d = 2.0 * scale * G23 * eSens[6];

  // Add the derivative to the transformation terms
  Td[0] += Xdinv[0] * A11d + Xdinv[3] * A21d;
  Td[1] += Xdinv[0] * A12d + Xdinv[3] * A22d;
  Td[3] += Xdinv[1] * A11d + Xdinv[4] * A21d;
  Td[4] += Xdinv[1] * A12d + Xdinv[4] * A22d;
  Td[6] += Xdinv[2] * A11d + Xdinv[5] * A21d;
  Td[7] += Xdinv[2] * A12d + Xdinv[5] * A22d;

  // Add the terms to the Xdinv
  Xdinvd[0] += T[0] * A11d + T[1] * A12d;
  Xdinvd[1] += T[3] * A11d + T[4] * A12d;
  Xdinvd[2] += T[6] * A11d + T[7] * A12d;

  Xdinvd[3] += T[0] * A21d + T[1] * A22d;
  Xdinvd[4] += T[3] * A21d + T[4] * A22d;
  Xdinvd[5] += T[6] * A21d + T[7] * A22d;

  // Add the terms from the tying strain interpolation
  for (int k = 0; k < 6; k++) {
    g13d[k] += N13[k] * G13d;
    g23d[k] += N23[k] * G23d;
  }
}

/*
  Add the tying coefficients required to compute the geometric
  stiffness matrix term
*/
void MITC9::addTyingGmatWeights(TacsScalar w13[], TacsScalar w23[],
                                const TacsScalar scalar, const TacsScalar s[],
                                const double N13[], const double N23[],
                                const TacsScalar Xdinv[],
                                const TacsScalar T[]) {
  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0] * T[0] + Xdinv[1] * T[3] + Xdinv[2] * T[6];
  A12 = Xdinv[0] * T[1] + Xdinv[1] * T[4] + Xdinv[2] * T[7];

  A21 = Xdinv[3] * T[0] + Xdinv[4] * T[3] + Xdinv[5] * T[6];
  A22 = Xdinv[3] * T[1] + Xdinv[4] * T[4] + Xdinv[5] * T[7];

  // Add the contributions to the weights
  for (int k = 0; k < 6; k++) {
    w13[k] += 2.0 * scalar * (s[6] * A12 + s[7] * A11) * N13[k];
    w23[k] += 2.0 * scalar * (s[6] * A22 + s[7] * A21) * N23[k];
  }
}

/*
  Add the contribution from the assumed strain derivatives
  interpolated from the tying points to the derivative of the strain
  at the present point.

  Note that this code is analagous to the addTyingStrain function
  except for the B matrix.

  input:
  N13:  the shape functions for the g13 tensorial strain
  N23:  the shape functions for the g23 tensorial strain
  B13:  the g13 tensorial strain at the tying points
  B23:  the g23 tensorial strain at the tying points
  T:    the transformation matrix

  output:
  B:    the derivative of the strain
*/
void MITC9::addTyingBmat(TacsScalar B[], const double N13[], const double N23[],
                         const TacsScalar B13[], const TacsScalar B23[],
                         const TacsScalar Xdinv[], const TacsScalar T[]) {
  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0] * T[0] + Xdinv[1] * T[3] + Xdinv[2] * T[6];
  A12 = Xdinv[0] * T[1] + Xdinv[1] * T[4] + Xdinv[2] * T[7];

  A21 = Xdinv[3] * T[0] + Xdinv[4] * T[3] + Xdinv[5] * T[6];
  A22 = Xdinv[3] * T[1] + Xdinv[4] * T[4] + Xdinv[5] * T[7];

  const int offset = 8 * NUM_NODES;
  for (int k = 0; k < offset; k++) {
    // Compute the strain using the assumed strain distribution
    // and the strain evaluated at the tying points
    TacsScalar G13, G23;
    G13 = (B13[0] * N13[0] + B13[offset] * N13[1] + B13[2 * offset] * N13[2] +
           B13[3 * offset] * N13[3] + B13[4 * offset] * N13[4] +
           B13[5 * offset] * N13[5]);
    G23 = (B23[0] * N23[0] + B23[offset] * N23[1] + B23[2 * offset] * N23[2] +
           B23[3 * offset] * N23[3] + B23[4 * offset] * N23[4] +
           B23[5 * offset] * N23[5]);

    // Compute and set the final strain values
    // e = 2*A^{T}*G
    B[6] = 2.0 * (A12 * G13 + A22 * G23);
    B[7] = 2.0 * (A11 * G13 + A21 * G23);

    B += 8;
    B13 += 1;
    B23 += 1;
  }
}

/*
  Add the sensitivities from the tying strain components
*/
void MITC9::addTyingBmatSens(TacsScalar B13d[], TacsScalar B23d[],
                             TacsScalar Xdinvd[], TacsScalar Td[],
                             const TacsScalar eSens[], const TacsScalar psi[],
                             const double N13[], const double N23[],
                             const TacsScalar B13[], const TacsScalar B23[],
                             const TacsScalar Xdinv[], const TacsScalar T[]) {
  // Compute the coefficients for the strain transformation. Note
  // that the remaining values in the matrix A = Xdinv*T are either
  // zero or unity.
  TacsScalar A11, A12, A21, A22;

  // A = Xdinv*T
  A11 = Xdinv[0] * T[0] + Xdinv[1] * T[3] + Xdinv[2] * T[6];
  A12 = Xdinv[0] * T[1] + Xdinv[1] * T[4] + Xdinv[2] * T[7];

  A21 = Xdinv[3] * T[0] + Xdinv[4] * T[3] + Xdinv[5] * T[6];
  A22 = Xdinv[3] * T[1] + Xdinv[4] * T[4] + Xdinv[5] * T[7];

  TacsScalar A11d = 0.0, A12d = 0.0;
  TacsScalar A21d = 0.0, A22d = 0.0;

  const int offset = 8 * NUM_NODES;
  for (int k = 0; k < offset; k++) {
    // Compute the strain using the assumed strain distribution
    // and the strain evaluated at the tying points
    TacsScalar G13, G23;
    G13 = (B13[0] * N13[0] + B13[offset] * N13[1] + B13[2 * offset] * N13[2] +
           B13[3 * offset] * N13[3] + B13[4 * offset] * N13[4] +
           B13[5 * offset] * N13[5]);
    G23 = (B23[0] * N23[0] + B23[offset] * N23[1] + B23[2 * offset] * N23[2] +
           B23[3 * offset] * N23[3] + B23[4 * offset] * N23[4] +
           B23[5 * offset] * N23[5]);

    // Compute and set the final strain values
    // e = 2*A^{T}*G
    TacsScalar G13d = 2.0 * psi[0] * (A12 * eSens[6] + A11 * eSens[7]);
    TacsScalar G23d = 2.0 * psi[0] * (A22 * eSens[6] + A21 * eSens[7]);

    A11d += 2.0 * psi[0] * G13 * eSens[7];
    A12d += 2.0 * psi[0] * G13 * eSens[6];
    A21d += 2.0 * psi[0] * G23 * eSens[7];
    A22d += 2.0 * psi[0] * G23 * eSens[6];

    for (int i = 0; i < 6; i++) {
      B13d[offset * i] += G13d * N13[i];
      B23d[offset * i] += G23d * N23[i];
    }

    B13++;
    B23++;
    B13d++;
    B23d++;
    psi++;
  }

  // Add the derivative to the transformation terms
  Td[0] += Xdinv[0] * A11d + Xdinv[3] * A21d;
  Td[1] += Xdinv[0] * A12d + Xdinv[3] * A22d;
  Td[3] += Xdinv[1] * A11d + Xdinv[4] * A21d;
  Td[4] += Xdinv[1] * A12d + Xdinv[4] * A22d;
  Td[6] += Xdinv[2] * A11d + Xdinv[5] * A21d;
  Td[7] += Xdinv[2] * A12d + Xdinv[5] * A22d;

  // Add the terms to the Xdinv
  Xdinvd[0] += T[0] * A11d + T[1] * A12d;
  Xdinvd[1] += T[3] * A11d + T[4] * A12d;
  Xdinvd[2] += T[6] * A11d + T[7] * A12d;

  Xdinvd[3] += T[0] * A21d + T[1] * A22d;
  Xdinvd[4] += T[3] * A21d + T[4] * A22d;
  Xdinvd[5] += T[6] * A21d + T[7] * A22d;
}

/*
  Evaluate the in-plane rotation penalty term

  This term is used to add an artificial drilling stiffness to the
  plate to enable the use of full rotations at the element nodes.

  The rotation penalty term is calculated as follows:

  rot = Xa^{T}*C*(Xb + Ub) - Xb^{T}*C*(Xa + Ua)

  In this code, the approximate rotation matrix is interpolated from
  the nodes of the mesh and is computed directly from the directors
*/
TacsScalar MITC9::computeRotPenalty(const double N[], const TacsScalar Xa[],
                                    const TacsScalar Xb[],
                                    const TacsScalar Ua[],
                                    const TacsScalar Ub[],
                                    const TacsScalar vars[]) {
  TacsScalar Ci[9];
  Ci[0] = Ci[1] = Ci[2] = 0.0;
  Ci[3] = Ci[4] = Ci[5] = 0.0;
  Ci[6] = Ci[7] = Ci[8] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the quaternions
    const TacsScalar *q = &vars[8 * i + 3];

    // Compute the full rotation matrix
    TacsScalar C[9];
    C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
    C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
    C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

    C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
    C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
    C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

    C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
    C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
    C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

    // Interpolate the rotation matrix
    Ci[0] += N[i] * C[0];
    Ci[1] += N[i] * C[1];
    Ci[2] += N[i] * C[2];
    Ci[3] += N[i] * C[3];
    Ci[4] += N[i] * C[4];
    Ci[5] += N[i] * C[5];
    Ci[6] += N[i] * C[6];
    Ci[7] += N[i] * C[7];
    Ci[8] += N[i] * C[8];
  }

  // Compute va = C*(Xa + Ua)
  TacsScalar va[3];
  matMult(Ci, Xa, va);
  matMultAdd(Ci, Ua, va);

  // Compute vb = C*(Xa + Ua)
  TacsScalar vb[3];
  matMult(Ci, Xb, vb);
  matMultAdd(Ci, Ub, vb);

  // Compute the penalty term
  return vecDot(Xa, vb) - vecDot(Xb, va);
}

/*
  Evaluate the derivative of the in-plane rotation penalty term

  This term is used to add an artificial drilling stiffness to the
  plate to enable the use of full rotations at the element nodes.

  The rotation penalty term is calculated as follows:

  rot = Xa^{T}*C*(Xb + Ub) - Xb^{T}*C*(Xa + Ua)

  In this code, the approximate rotation matrix is interpolated from
  the nodes of the mesh and is computed directly from the directors
*/
TacsScalar MITC9::computeBRotPenalty(
    TacsScalar brot[], const double N[], const double Na[], const double Nb[],
    const TacsScalar Xa[], const TacsScalar Xb[], const TacsScalar Ua[],
    const TacsScalar Ub[], const TacsScalar vars[]) {
  TacsScalar Ci[9];
  Ci[0] = Ci[1] = Ci[2] = 0.0;
  Ci[3] = Ci[4] = Ci[5] = 0.0;
  Ci[6] = Ci[7] = Ci[8] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the quaternions
    const TacsScalar *q = &vars[8 * i + 3];

    // Compute the rotation matrix
    TacsScalar C[9];
    C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
    C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
    C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

    C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
    C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
    C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

    C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
    C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
    C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

    // Interpolate the rotation matrix to the nodes
    Ci[0] += N[i] * C[0];
    Ci[1] += N[i] * C[1];
    Ci[2] += N[i] * C[2];
    Ci[3] += N[i] * C[3];
    Ci[4] += N[i] * C[4];
    Ci[5] += N[i] * C[5];
    Ci[6] += N[i] * C[6];
    Ci[7] += N[i] * C[7];
    Ci[8] += N[i] * C[8];
  }

  // Compute va = Ci*(Xa + Ua)
  TacsScalar va[3], ta[3];
  ta[0] = Xa[0] + Ua[0];
  ta[1] = Xa[1] + Ua[1];
  ta[2] = Xa[2] + Ua[2];
  matMult(Ci, ta, va);

  // Compute vb = Ci*(Xb + Ub)
  TacsScalar vb[3], tb[3];
  tb[0] = Xb[0] + Ub[0];
  tb[1] = Xb[1] + Ub[1];
  tb[2] = Xb[2] + Ub[2];
  matMult(Ci, tb, vb);

  // Compute wa = Ci^{T}*Xa and wb = Ci^{T}*Xb
  TacsScalar wa[3], wb[3];
  matMultTrans(Ci, Xa, wa);
  matMultTrans(Ci, Xb, wb);

  // Add the values to the Brot array
  TacsScalar *b = brot;
  for (int i = 0; i < NUM_NODES; i++) {
    // wa*Nb - wb*Na
    b[0] = Nb[i] * wa[0] - Na[i] * wb[0];
    b[1] = Nb[i] * wa[1] - Na[i] * wb[1];
    b[2] = Nb[i] * wa[2] - Na[i] * wb[2];

    // Add the terms from the derivatives w.r.t. Ci. Note that
    // each term adds a contribution which takes the following form:
    //  N[i]*(Xa^{T}*D*tb - Xb^{T}*D*ta)
    const TacsScalar *q = &vars[8 * i + 3];
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0 * q[3];
    Q[2] = -2.0 * q[2];

    Q[3] = -2.0 * q[3];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[1];

    Q[6] = 2.0 * q[2];
    Q[7] = -2.0 * q[1];
    Q[8] = 0.0;
    b[3] = N[i] * (mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0 * q[2];
    Q[2] = 2.0 * q[3];

    Q[3] = 2.0 * q[2];
    Q[4] = -4.0 * q[1];
    Q[5] = 2.0 * q[0];

    Q[6] = 2.0 * q[3];
    Q[7] = -2.0 * q[0];
    Q[8] = -4.0 * q[1];
    b[4] = N[i] * (mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[2]
    Q[0] = -4.0 * q[2];
    Q[1] = 2.0 * q[1];
    Q[2] = -2.0 * q[0];

    Q[3] = 2.0 * q[1];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[3];

    Q[6] = 2.0 * q[0];
    Q[7] = 2.0 * q[3];
    Q[8] = -4.0 * q[2];
    b[5] = N[i] * (mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    // Derivative w.r.t. q[3]
    Q[0] = -4.0 * q[3];
    Q[1] = 2.0 * q[0];
    Q[2] = 2.0 * q[1];

    Q[3] = -2.0 * q[0];
    Q[4] = -4.0 * q[3];
    Q[5] = 2.0 * q[2];

    Q[6] = 2.0 * q[1];
    Q[7] = 2.0 * q[2];
    Q[8] = 0.0;
    b[6] = N[i] * (mat3x3Inner(Q, Xa, tb) - mat3x3Inner(Q, Xb, ta));

    b[7] = 0.0;

    b += 8;
  }

  // Compute the penalty term
  return vecDot(Xa, vb) - vecDot(Xb, va);
}

/*
  Add the derivative of the rotational penalty term with respect to
  the inputs Xa/Xb
*/
void MITC9::addBRotPenaltySens(TacsScalar _Xad[], TacsScalar _Xbd[],
                               const TacsScalar rotd, const TacsScalar scale,
                               const TacsScalar psi[], const double N[],
                               const double Na[], const double Nb[],
                               const TacsScalar Xa[], const TacsScalar Xb[],
                               const TacsScalar Ua[], const TacsScalar Ub[],
                               const TacsScalar vars[]) {
  TacsScalar Ci[9];
  Ci[0] = Ci[1] = Ci[2] = 0.0;
  Ci[3] = Ci[4] = Ci[5] = 0.0;
  Ci[6] = Ci[7] = Ci[8] = 0.0;

  for (int i = 0; i < NUM_NODES; i++) {
    // Set the pointer to the quaternions
    const TacsScalar *q = &vars[8 * i + 3];

    // Compute the rotation matrix
    TacsScalar C[9];
    C[0] = 1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);
    C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
    C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

    C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
    C[4] = 1.0 - 2.0 * (q[1] * q[1] + q[3] * q[3]);
    C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

    C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
    C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
    C[8] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

    // Interpolate the rotation matrix to the nodes
    Ci[0] += N[i] * C[0];
    Ci[1] += N[i] * C[1];
    Ci[2] += N[i] * C[2];
    Ci[3] += N[i] * C[3];
    Ci[4] += N[i] * C[4];
    Ci[5] += N[i] * C[5];
    Ci[6] += N[i] * C[6];
    Ci[7] += N[i] * C[7];
    Ci[8] += N[i] * C[8];
  }

  // Compute va = Ci*(Xa + Ua)
  TacsScalar va[3], ta[3];
  ta[0] = Xa[0] + Ua[0];
  ta[1] = Xa[1] + Ua[1];
  ta[2] = Xa[2] + Ua[2];
  matMult(Ci, ta, va);

  // Compute vb = Ci*(Xb + Ub)
  TacsScalar vb[3], tb[3];
  tb[0] = Xb[0] + Ub[0];
  tb[1] = Xb[1] + Ub[1];
  tb[2] = Xb[2] + Ub[2];
  matMult(Ci, tb, vb);

  TacsScalar tad[3], tbd[3];
  tad[0] = tad[1] = tad[2] = 0.0;
  tbd[0] = tbd[1] = tbd[2] = 0.0;

  // Compute wa = Ci^{T}*Xa and wb = Ci^{T}*Xb
  TacsScalar wa[3], wb[3];
  matMultTrans(Ci, Xa, wa);
  matMultTrans(Ci, Xb, wb);

  // Accumulate the derivative contributions
  TacsScalar wad[3], wbd[3];
  wad[0] = wad[1] = wad[2] = 0.0;
  wbd[0] = wbd[1] = wbd[2] = 0.0;

  TacsScalar Xad[3], Xbd[3];
  Xad[0] = Xad[1] = Xad[2] = 0.0;
  Xbd[0] = Xbd[1] = Xbd[2] = 0.0;

  // Add the values from the
  const TacsScalar *p = psi;
  for (int i = 0; i < NUM_NODES; i++) {
    // wa*Nb - wb*Na
    wad[0] += Nb[i] * p[0];
    wad[1] += Nb[i] * p[1];
    wad[2] += Nb[i] * p[2];

    wbd[0] -= Na[i] * p[0];
    wbd[1] -= Na[i] * p[1];
    wbd[2] -= Na[i] * p[2];

    // Add the terms from the derivatives w.r.t. Ci. Note that
    // each term adds a contribution which takes the following form:
    //  N[i]*(Xa^{T}*D*tb - Xb^{T}*D*ta)
    const TacsScalar *q = &vars[8 * i + 3];
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0 * q[3];
    Q[2] = -2.0 * q[2];

    Q[3] = -2.0 * q[3];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[1];

    Q[6] = 2.0 * q[2];
    Q[7] = -2.0 * q[1];
    Q[8] = 0.0;

    // Perform the required matrix-matrix multiplications
    TacsScalar t[3];
    TacsScalar s = N[i] * p[3];
    matMultTrans(Q, Xa, t);
    vecAxpy(s, t, tbd);
    matMult(Q, tb, t);
    vecAxpy(s, t, Xad);

    matMultTrans(Q, Xb, t);
    vecAxpy(-s, t, tad);
    matMult(Q, ta, t);
    vecAxpy(-s, t, Xbd);

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0 * q[2];
    Q[2] = 2.0 * q[3];

    Q[3] = 2.0 * q[2];
    Q[4] = -4.0 * q[1];
    Q[5] = 2.0 * q[0];

    Q[6] = 2.0 * q[3];
    Q[7] = -2.0 * q[0];
    Q[8] = -4.0 * q[1];

    s = N[i] * p[4];
    matMultTrans(Q, Xa, t);
    vecAxpy(s, t, tbd);
    matMult(Q, tb, t);
    vecAxpy(s, t, Xad);

    matMultTrans(Q, Xb, t);
    vecAxpy(-s, t, tad);
    matMult(Q, ta, t);
    vecAxpy(-s, t, Xbd);

    // Derivative w.r.t. q[2]
    Q[0] = -4.0 * q[2];
    Q[1] = 2.0 * q[1];
    Q[2] = -2.0 * q[0];

    Q[3] = 2.0 * q[1];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[3];

    Q[6] = 2.0 * q[0];
    Q[7] = 2.0 * q[3];
    Q[8] = -4.0 * q[2];

    s = N[i] * p[5];
    matMultTrans(Q, Xa, t);
    vecAxpy(s, t, tbd);
    matMult(Q, tb, t);
    vecAxpy(s, t, Xad);

    matMultTrans(Q, Xb, t);
    vecAxpy(-s, t, tad);
    matMult(Q, ta, t);
    vecAxpy(-s, t, Xbd);

    // Derivative w.r.t. q[3]
    Q[0] = -4.0 * q[3];
    Q[1] = 2.0 * q[0];
    Q[2] = 2.0 * q[1];

    Q[3] = -2.0 * q[0];
    Q[4] = -4.0 * q[3];
    Q[5] = 2.0 * q[2];

    Q[6] = 2.0 * q[1];
    Q[7] = 2.0 * q[2];
    Q[8] = 0.0;

    s = N[i] * p[6];
    matMultTrans(Q, Xa, t);
    vecAxpy(s, t, tbd);
    matMult(Q, tb, t);
    vecAxpy(s, t, Xad);

    matMultTrans(Q, Xb, t);
    vecAxpy(-s, t, tad);
    matMult(Q, ta, t);
    vecAxpy(-s, t, Xbd);

    p += 8;
  }

  // Compute the penalty term
  matMultAdd(Ci, wad, Xad);
  matMultAdd(Ci, wbd, Xbd);

  // Add the result scaled by the value to the output
  TacsScalar t[3];
  matMultTrans(Ci, Xb, t);
  _Xad[0] += scale * (Xad[0] + tad[0]) + rotd * (vb[0] - t[0]);
  _Xad[1] += scale * (Xad[1] + tad[1]) + rotd * (vb[1] - t[1]);
  _Xad[2] += scale * (Xad[2] + tad[2]) + rotd * (vb[2] - t[2]);

  matMultTrans(Ci, Xa, t);
  _Xbd[0] += scale * (Xbd[0] + tbd[0]) + rotd * (t[0] - va[0]);
  _Xbd[1] += scale * (Xbd[1] + tbd[1]) + rotd * (t[1] - va[1]);
  _Xbd[2] += scale * (Xbd[2] + tbd[2]) + rotd * (t[2] - va[2]);
}

/*
  Add the second derivative of the in-plane penalty term to the
  stiffness matrix.
*/
void MITC9::addGRotMat(TacsScalar J[], const TacsScalar scale, const double N[],
                       const double Na[], const double Nb[],
                       const TacsScalar Xa[], const TacsScalar Xb[],
                       const TacsScalar Ua[], const TacsScalar Ub[],
                       const TacsScalar vars[]) {
  // Compute va = Ci*(Xa + Ua)
  TacsScalar ta[3];
  ta[0] = Xa[0] + Ua[0];
  ta[1] = Xa[1] + Ua[1];
  ta[2] = Xa[2] + Ua[2];

  // Compute vb = Ci*(Xb + Ub)
  TacsScalar tb[3];
  tb[0] = Xb[0] + Ub[0];
  tb[1] = Xb[1] + Ub[1];
  tb[2] = Xb[2] + Ub[2];

  // Compute the second derivative w.r.t. the quaternion
  TacsScalar dCXa[3 * 9], dCXb[3 * 9];
  computeQtr2ndDeriv(Xa, dCXa);
  computeQtr2ndDeriv(Xb, dCXb);

  // Pre-compute terms for the second derivatives
  TacsScalar Wa[12 * NUM_NODES], Wb[12 * NUM_NODES];
  TacsScalar *wa = Wa, *wb = Wb;

  for (int i = 0; i < NUM_NODES; i++) {
    const TacsScalar *q = &vars[8 * i + 3];
    TacsScalar Q[9];
    Q[0] = 0.0;
    Q[1] = 2.0 * q[3];
    Q[2] = -2.0 * q[2];

    Q[3] = -2.0 * q[3];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[1];

    Q[6] = 2.0 * q[2];
    Q[7] = -2.0 * q[1];
    Q[8] = 0.0;
    matMultTrans(Q, Xa, &wa[0]);
    matMultTrans(Q, Xb, &wb[0]);

    // Derivative w.r.t. q[1]
    Q[0] = 0.0;
    Q[1] = 2.0 * q[2];
    Q[2] = 2.0 * q[3];

    Q[3] = 2.0 * q[2];
    Q[4] = -4.0 * q[1];
    Q[5] = 2.0 * q[0];

    Q[6] = 2.0 * q[3];
    Q[7] = -2.0 * q[0];
    Q[8] = -4.0 * q[1];
    matMultTrans(Q, Xa, &wa[3]);
    matMultTrans(Q, Xb, &wb[3]);

    // Derivative w.r.t. q[2]
    Q[0] = -4.0 * q[2];
    Q[1] = 2.0 * q[1];
    Q[2] = -2.0 * q[0];

    Q[3] = 2.0 * q[1];
    Q[4] = 0.0;
    Q[5] = 2.0 * q[3];

    Q[6] = 2.0 * q[0];
    Q[7] = 2.0 * q[3];
    Q[8] = -4.0 * q[2];
    matMultTrans(Q, Xa, &wa[6]);
    matMultTrans(Q, Xb, &wb[6]);

    // Derivative w.r.t. q[3]
    Q[0] = -4.0 * q[3];
    Q[1] = 2.0 * q[0];
    Q[2] = 2.0 * q[1];

    Q[3] = -2.0 * q[0];
    Q[4] = -4.0 * q[3];
    Q[5] = 2.0 * q[2];

    Q[6] = 2.0 * q[1];
    Q[7] = 2.0 * q[2];
    Q[8] = 0.0;
    matMultTrans(Q, Xa, &wa[9]);
    matMultTrans(Q, Xb, &wb[9]);

    // Increment the pointers
    wa += 12;
    wb += 12;
  }

  const int iv[] = {3, 3, 3, 4, 4, 4, 5, 5, 6};
  const int jv[] = {4, 5, 6, 4, 5, 6, 5, 6, 6};

  // Add the contributions from the second derivatives of the quaternions
  for (int i = 0; i < NUM_NODES; i++) {
    // Add the second derivative coupling between the drilling
    // rotation and the in-plane displacement
    wa = Wa;
    wb = Wb;
    for (int j = 0; j < NUM_NODES; j++) {
      for (int ii = 0; ii < 4; ii++) {    // Loop over the quaternions
        for (int jj = 0; jj < 3; jj++) {  // Loop over the displacements
          J[8 * NUM_NODES * (8 * i + jj) + 8 * j + 3 + ii] +=
              scale * N[j] *
              (Nb[i] * wa[3 * ii + jj] - Na[i] * wb[3 * ii + jj]);

          J[8 * NUM_NODES * (8 * j + 3 + ii) + 8 * i + jj] +=
              scale * N[j] *
              (Nb[i] * wa[3 * ii + jj] - Na[i] * wb[3 * ii + jj]);
        }
      }
      wa += 12;
      wb += 12;
    }

    // Set pointers to the second derivatives
    const TacsScalar *dCa = dCXa, *dCb = dCXb;

    // Compute the second derivatives from the quaternions alone
    for (int ii = 0; ii < 9; ii++) {
      TacsScalar Jadd = scale * N[i] * (vecDot(tb, dCa) - vecDot(ta, dCb));
      J[(8 * NUM_NODES) * (8 * i + iv[ii]) + (8 * i + jv[ii])] += Jadd;
      if (iv[ii] != jv[ii]) {
        J[(8 * NUM_NODES) * (8 * i + jv[ii]) + (8 * i + iv[ii])] += Jadd;
      }
      dCa += 3;
      dCb += 3;
    }
  }
}

/*
  Get the constitutive object
*/
TACSConstitutive *MITC9::getConstitutive() { return stiff; }

/*
  Return the number of quadrature points
*/
int MITC9::getNumGaussPts() { return ORDER * ORDER; }

/*
  Return the quadrature points and weights
*/
double MITC9::getGaussWtsPts(const int num, double pt[]) {
  int m = (int)(num / ORDER);
  int n = num % ORDER;
  pt[0] = gaussPts[n];
  pt[1] = gaussPts[m];

  return gaussWts[n] * gaussWts[m];
}

/*
  Get the values of the shape functions
*/
void MITC9::getShapeFunctions(const double pt[], double N[]) {
  computeShapeFunc(pt[0], pt[1], N);
}

/*
  Retrieve the determinant of the Jacobian transformation matrix
*/
TacsScalar MITC9::getDetJacobian(const double pt[], const TacsScalar X[]) {
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the determinant of the Jacobian transformation
  return det3x3(Xd);
}

/*
  Evaluate the determinant of the Jacobian transformation w.r.t.
  the node locations
*/
TacsScalar MITC9::getDetJacobianXptSens(TacsScalar *hXptSens, const double pt[],
                                        const TacsScalar X[]) {
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the determinant of the Jacobian transformation
  TacsScalar h = det3x3(Xd);

  // Compute the derivative of the determinant
  TacsScalar Xdd[9];
  det3x3Sens(Xd, Xdd);

  // Extract/add the sensitivities from the frame
  TacsScalar fnd[3];
  fnd[0] = Xdd[2];
  fnd[1] = Xdd[5];
  fnd[2] = Xdd[8];

  // Add the contributions to Xad and Xbd
  TacsScalar Xad[3], Xbd[3];
  Xad[0] = Xdd[0];
  Xad[1] = Xdd[3];
  Xad[2] = Xdd[6];
  Xbd[0] = Xdd[1];
  Xbd[1] = Xdd[4];
  Xbd[2] = Xdd[7];

  // // Compute the frame normal
  TacsScalar Xrd[9 * NUM_NODES];
  memset(Xrd, 0, 9 * NUM_NODES * sizeof(TacsScalar));
  addFrameNormalSens(fnd, N, Xrd);

  // Add the derivatives the shape function directions
  for (int k = 0; k < NUM_NODES; k++) {
    hXptSens[3 * k] = Na[k] * Xad[0] + Nb[k] * Xbd[0];
    hXptSens[3 * k + 1] = Na[k] * Xad[1] + Nb[k] * Xbd[1];
    hXptSens[3 * k + 2] = Na[k] * Xad[2] + Nb[k] * Xbd[2];
  }

  addFramesSens(hXptSens, Xrd, X);

  return h;
}

/*
  Evaluate the strain at a parametric point within the element
*/
void MITC9::getStrain(TacsScalar e[], const double pt[], const TacsScalar X[],
                      const TacsScalar vars[]) {
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3 * NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9];
  computeTransform(T, Xa, Xb);

  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);

  // Compute the displacement-based strain
  evalStrain(e, Ur, dr, Xdinv, zXdinv, T);

  // Add the contribution from the tying strain
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
}

/*
  Add the derivative of the product of the array esens with the strain
  with respect to the state variables
*/
void MITC9::addStrainSVSens(TacsScalar sens[], const double pt[],
                            const TacsScalar scale, const TacsScalar esens[],
                            const TacsScalar X[], const TacsScalar vars[]) {
  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the derivatives of the directors
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];
  computeDirectors(dir, vars, Xr);
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the derivative of the tying strain
  TacsScalar g13[6], g23[6];
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  // Set the parametric locations
  const double u = pt[0];
  const double v = pt[1];

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the derivatives of the shape functions
  TacsScalar Xdinv[9];
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9];
  computeTransform(T, Xa, Xb);

  // Compute the displacement-based strain
  TacsScalar e[8], B[64 * NUM_NODES];
  evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

  // Add the contribution from the tying straint
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
  addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

  const TacsScalar *b = B;
  for (int ii = 0; ii < 8 * NUM_NODES; ii++) {
    sens[ii] += scale * (b[0] * esens[0] + b[1] * esens[1] + b[2] * esens[2] +
                         b[3] * esens[3] + b[4] * esens[4] + b[5] * esens[5] +
                         b[6] * esens[6] + b[7] * esens[7]);
    b += 8;
  }
}

/*
  Add the derivative of the strain w.r.t. the node locations
*/
void MITC9::addStrainXptSens(TacsScalar eXpt[], const double pt[],
                             const TacsScalar scale, const TacsScalar eSens[],
                             const TacsScalar X[], const TacsScalar vars[]) {
  // Set the u/v locations
  const double u = pt[0];
  const double v = pt[1];

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3 * NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the inverse of the transformation matrix
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9];
  computeTransform(T, Xa, Xb);

  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);

  // Compute the derivative w.r.t. Ur, dr, Xdinv, zXdinv and T
  TacsScalar Urd[9], drd[9], Xdinvd[9], zXdinvd[9], Td[9];
  evalStrainSens(Urd, drd, Xdinvd, zXdinvd, Td, scale, eSens, Ur, dr, Xdinv,
                 zXdinv, T);

  // Add the contribution from the tying strain
  TacsScalar g13d[6], g23d[6];
  memset(g13d, 0, 6 * sizeof(TacsScalar));
  memset(g23d, 0, 6 * sizeof(TacsScalar));
  addTyingStrainSens(g13d, g23d, Xdinvd, Td, scale, eSens, N13, N23, g13, g23,
                     Xdinv, T);

  // Add the contributions to dird
  TacsScalar dird[3 * NUM_NODES];
  for (int k = 0; k < NUM_NODES; k++) {
    dird[3 * k] = Na[k] * drd[0] + Nb[k] * drd[1] + N[k] * Urd[2];
    dird[3 * k + 1] = Na[k] * drd[3] + Nb[k] * drd[4] + N[k] * Urd[5];
    dird[3 * k + 2] = Na[k] * drd[6] + Nb[k] * drd[7] + N[k] * Urd[8];
  }

  // Compute the derivatives Xad and Xbd
  TacsScalar Xad[3], Xbd[3];
  computeTransformSens(Xad, Xbd, Td, Xa, Xb);

  // Compute the derivatives Xrd
  TacsScalar Xrd[9 * NUM_NODES];
  memset(Xrd, 0, 9 * NUM_NODES * sizeof(TacsScalar));
  addNormalRateMatSens(Xrd, Xdinvd, zXdinvd, Na, Nb, Xr, Xdinv);

  // Compute the inverse of the transformation matrix
  TacsScalar Xdd[9];
  inv3x3Sens(Xdd, Xdinvd, Xdinv);

  // Extract/add the sensitivities from the frame
  TacsScalar fnd[3];
  fnd[0] = Xdd[2];
  fnd[1] = Xdd[5];
  fnd[2] = Xdd[8];

  // Add the contributions to Xad and Xbd
  Xad[0] += Xdd[0];
  Xad[1] += Xdd[3];
  Xad[2] += Xdd[6];
  Xbd[0] += Xdd[1];
  Xbd[1] += Xdd[4];
  Xbd[2] += Xdd[7];

  // // Compute the frame normal
  addFrameNormalSens(fnd, N, Xrd);

  // Add the derivatives the shape function directions
  for (int k = 0; k < NUM_NODES; k++) {
    eXpt[3 * k] += Na[k] * Xad[0] + Nb[k] * Xbd[0];
    eXpt[3 * k + 1] += Na[k] * Xad[1] + Nb[k] * Xbd[1];
    eXpt[3 * k + 2] += Na[k] * Xad[2] + Nb[k] * Xbd[2];
  }

  // input: g13d, g23d,
  addComputeTyingStrainSens(eXpt, Xrd, dird, g13d, g23d, X, Xr, vars, dir);
  addDirectorsSens(Xrd, dird, vars);
  addFramesSens(eXpt, Xrd, X);
}

/*
  Determine the number of nodes and elements for visualization
  generated by the data in this element. Note that higher-order
  elements are broken down into bi-linear elements for visualization.

  output:
  nelems:  the number of visualization elements
  nnodes:  the number of nodes used for visualization
  ncsr:    the number of entries in a CSR-type data structure used
  to store the connectivity
*/
void MITC9::addOutputCount(int *nelems, int *nnodes, int *ncsr) {
  *nelems += 4;
  *nnodes += 9;
  *ncsr += 16;
}

/*
  Get the output data from this element and place it in a real
  array for visualization later. The values generated for visualization
  are determined by a bit-wise selection variable 'out_type' which is
  can be used to simultaneously write out different data. Note that this
  is why the bitwise operation & is used below.

  The output may consist of the following:
  - the nodal locations
  - the displacements and rotations
  - the strains or strains within the element
  - extra variables that are used for optimization

  output:
  data:     the data to write to the file (eventually)

  input:
  out_type: the bit-wise variable used to specify what data to generate
  vars:     the element variables
  Xpts:     the element nodal locations
*/
void MITC9::getOutputData(unsigned int out_type, double *data, int ld_data,
                          const TacsScalar Xpts[], const TacsScalar vars[]) {
  for (int m = 0, p = 0; m < 3; m++) {
    for (int n = 0; n < 3; n++, p++) {
      double pt[2];
      pt[0] = -1.0 + 1.0 * n;
      pt[1] = -1.0 + 1.0 * m;

      TacsScalar strain[8], stress[8];
      getStrain(strain, pt, Xpts, vars);

      int index = 0;
      if (out_type & TACSElement::OUTPUT_NODES) {
        for (int k = 0; k < 3; k++) {
          data[index + k] = TacsRealPart(Xpts[3 * p + k]);
        }
        index += 3;
      }
      if (out_type & TACSElement::OUTPUT_DISPLACEMENTS) {
        for (int k = 0; k < NUM_DISPS; k++) {
          data[index + k] = TacsRealPart(vars[NUM_DISPS * p + k]);
        }
        index += NUM_DISPS;
      }
      if (out_type & TACSElement::OUTPUT_STRAINS) {
        // Add the term due to the potential energy
        for (int k = 0; k < NUM_STRESSES; k++) {
          data[index + k] = TacsRealPart(strain[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_STRESSES) {
        // Evaluate the stiffness at the current point
        // and then calculate the stress
        stiff->calculateStress(pt, strain, stress);

        for (int k = 0; k < NUM_STRESSES; k++) {
          data[index + k] = TacsRealPart(stress[k]);
        }
        index += NUM_STRESSES;
      }
      if (out_type & TACSElement::OUTPUT_EXTRAS) {
        // Compute the failure value
        TacsScalar lambda;
        stiff->failure(pt, strain, &lambda);
        data[index] = TacsRealPart(lambda);

        // Compute the buckling constraint value
        TacsScalar bval;
        stiff->buckling(strain, &bval);
        data[index + 1] = TacsRealPart(bval);

        data[index + 2] = TacsRealPart(stiff->getDVOutputValue(0, pt));
        data[index + 3] = TacsRealPart(stiff->getDVOutputValue(1, pt));

        index += NUM_EXTRAS;
      }
      if (out_type & TACSElement::OUTPUT_COORDINATES) {
        index += 9;
      }

      data += ld_data;
    }
  }
}

/*
  Get the element connectivity for visualization purposes. Since each
  element consists of a series of sub-elements used for visualization,
  we also need the connectivity of these visualization elements.

  output:
  con:  the connectivity of the local visualization elements contributed
  by this finite-element

  input:
  node:  the node offset number - so that this connectivity is more or
  less global
*/
void MITC9::getOutputConnectivity(int *con, int node) {
  int p = 0;
  for (int m = 0; m < 2; m++) {
    for (int n = 0; n < 2; n++) {
      con[4 * p] = node + n + 3 * m;
      con[4 * p + 1] = node + n + 1 + 3 * m;
      con[4 * p + 2] = node + n + 1 + 3 * (m + 1);
      con[4 * p + 3] = node + n + 3 * (m + 1);
      p++;
    }
  }
}

/*
  Test the implementation of the strain by comparing against a
  rigid-body displacement and rotation.
*/
void MITC9::testStrain(const TacsScalar X[]) {
  TacsScalar vars[8 * NUM_NODES];
  memset(vars, 0, sizeof(vars));

  const double u = 0.125;
  const double v = -0.63;

  // Pick a random set of values for the quaternions
  TacsScalar q[4];
  q[1] = 0.25;
  q[2] = 0.5;
  q[3] = 0.125;
  q[0] = sqrt(1.0 - vecDot(&q[1], &q[1]));

  // Compute the rotation matrix C = rot - I
  TacsScalar C[9];
  C[0] = -2.0 * (q[2] * q[2] + q[3] * q[3]);
  C[1] = 2.0 * (q[1] * q[2] + q[3] * q[0]);
  C[2] = 2.0 * (q[1] * q[3] - q[2] * q[0]);

  C[3] = 2.0 * (q[2] * q[1] - q[3] * q[0]);
  C[4] = -2.0 * (q[1] * q[1] + q[3] * q[3]);
  C[5] = 2.0 * (q[2] * q[3] + q[1] * q[0]);

  C[6] = 2.0 * (q[3] * q[1] + q[2] * q[0]);
  C[7] = 2.0 * (q[3] * q[2] - q[1] * q[0]);
  C[8] = -2.0 * (q[1] * q[1] + q[2] * q[2]);

  // Set the rigid displacement
  TacsScalar u0[3] = {1.25, -2.5, -4.0};

  // Compute the variables and set the quaternion values
  for (int k = 0; k < NUM_NODES; k++) {
    // Compute the displacements
    matMultTrans(C, &X[3 * k], &vars[8 * k]);
    for (int i = 0; i < 3; i++) {
      vars[8 * k + i] += u0[i];
    }

    // Copy the values of the quaternions
    memcpy(&vars[8 * k + 3], q, 4 * sizeof(TacsScalar));
  }

  // Compute the strain given the rigid rotation/translation
  // -------------------------------------------------------

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3 * NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(u, v, N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(u, v, Na, Nb);

  // Compute the derivative along the shape function
  // directions
  TacsScalar Xa[3], Xb[3];
  innerProduct(Na, X, Xa);
  innerProduct(Nb, X, Xb);

  // Compute the frame normal
  TacsScalar fn[3];
  computeFrameNormal(N, Xr, fn);

  // Evaluate the derivatives in the locally-aligned frame
  TacsScalar Xd[9], Xdinv[9];
  assembleFrame(Xa, Xb, fn, Xd);

  // Compute the Jacobian transformation
  inv3x3(Xd, Xdinv);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(u, v, N13, N23);

  // Compute the through-thickness derivative of [X,r]^{-1}
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  // Compute the derivatives of Ua/Ub along the given directions
  TacsScalar Ur[9], Ua[3], Ub[3], d[3];
  innerProduct8(Na, vars, Ua);
  innerProduct8(Nb, vars, Ub);
  innerProduct(N, dir, d);
  assembleFrame(Ua, Ub, d, Ur);

  // Now compute the derivatives of the director along each
  // coordinate direction
  TacsScalar dr[9], da[3], db[3], zero[3];
  innerProduct(Na, dir, da);
  innerProduct(Nb, dir, db);
  zero[0] = zero[1] = zero[2] = 0.0;
  assembleFrame(da, db, zero, dr);

  // Compute the rotation penalty
  computeRotPenalty(N, Xa, Xb, Ua, Ub, vars);

  // Compute the transformation to the locally-aligned frame
  TacsScalar T[9];
  computeTransform(T, Xa, Xb);

  // Compute the displacement-based strain
  TacsScalar e[8];
  evalStrain(e, Ur, dr, Xdinv, zXdinv, T);

  // Add the contribution from the tying strain
  addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);

  // Write out the error components
  TacsScalar fd[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  writeErrorComponents(stdout, "strain after rigid rotation", e, fd, 8);

  // Compute the derivatives of the directors
  TacsScalar dirdq[12 * NUM_NODES];
  computeDirectorDeriv(dirdq, vars, Xr);

  // Compute the B matrix
  TacsScalar B[64 * NUM_NODES];
  evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);

  // Add the tying strain
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
  addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

  // Compute the derivative of the strain w.r.t.
  double dh = 1e-6;

  for (int k = 0; k < 8 * NUM_NODES; k++) {
    TacsScalar vtmp = vars[k];
    vars[k] = vtmp + dh;

    // Compute the directors at the nodes
    computeDirectors(dir, vars, Xr);

    // Compute the tensorial shear strain at the tying points
    computeTyingStrain(g13, g23, X, Xr, vars, dir);

    // Compute the derivatives of Ua/Ub along the given directions
    innerProduct8(Na, vars, Ua);
    innerProduct8(Nb, vars, Ub);
    innerProduct(N, dir, d);
    assembleFrame(Ua, Ub, d, Ur);

    // Now compute the derivatives of the director along each
    // coordinate direction
    innerProduct(Na, dir, da);
    innerProduct(Nb, dir, db);
    zero[0] = zero[1] = zero[2] = 0.0;
    assembleFrame(da, db, zero, dr);

    evalStrain(fd, Ur, dr, Xdinv, zXdinv, T);

    // Add the contribution from the tying strain
    addTyingStrain(fd, N13, N23, g13, g23, Xdinv, T);

    for (int i = 0; i < 8; i++) {
      fd[i] = (fd[i] - e[i]) / dh;
    }

    vars[k] = vtmp;

    // Write out the error components
    char descript[64];
    sprintf(descript, "B%d", k);
    writeErrorComponents(stdout, descript, &B[8 * k], fd, 8);
  }
}

/*
  Test helper member functions which are needed to compute the
  geometric derivatives
*/
void MITC9::testXptSens(double dh) {
  testInv3x3Sens(dh);
  testStrainSens(dh);
  testTransformSens(dh);
  testNormalRateSens(dh);
  testTyingStrainSens(dh);
  testBmatSens(dh);
  testTyingBmatSens(dh);
  testBrotSens(dh);
  testFrameSens(dh);
}

/*
  Test the derivative of the 3x3 inverse computation
*/
void MITC9::testInv3x3Sens(double dh) {
  TacsScalar A[9], Ainvd[9];
  for (int i = 0; i < 9; i++) {
    A[i] = 1.0 + 2.0 * rand() / RAND_MAX;
    Ainvd[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute the inverse of the 3x3 matrix
  TacsScalar Ainv[9];
  inv3x3(A, Ainv);

  TacsScalar fval = 0.0;
  for (int j = 0; j < 9; j++) {
    fval += Ainv[j] * Ainvd[j];
  }

  // Compute the derivative
  TacsScalar Ad[9];
  inv3x3Sens(Ad, Ainvd, Ainv);

  // Form the finite-difference approximation for each component
  TacsScalar fd[9];
  for (int i = 0; i < 9; i++) {
    TacsScalar At = A[i];

    // Compute the finite-difference value
#ifdef TACS_USE_COMPLEX
    A[i] = A[i] + TacsScalar(0.0, dh);
#else
    A[i] = At + dh;
#endif
    inv3x3(A, Ainv);

    // Compute the function and the finite-difference result
    fd[i] = 0.0;
    for (int j = 0; j < 9; j++) {
      fd[i] += Ainv[j] * Ainvd[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    // Reset the value of A
    A[i] = At;
  }

  // Write out the error components
  writeErrorComponents(stdout, "inv3x3Sens", Ad, fd, 9);
}

/*
  Test the derivative of the strain
*/
void MITC9::testStrainSens(double dh) {
  // Select a random point and compute the shape functions
  double pt[2] = {-0.135, 0.223};
  double N13[6], N23[6];
  computeTyingFunc(pt[0], pt[1], N13, N23);

  // Compute the derivative of the strain
  TacsScalar Ur[9], dr[9];
  TacsScalar Xdinv[9], zXdinv[9];
  TacsScalar T[9];

  // Assign random inputs
  for (int i = 0; i < 9; i++) {
    Ur[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    dr[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Xdinv[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    zXdinv[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    T[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Assign random tying strain
  TacsScalar g13[6], g23[6];
  for (int i = 0; i < 6; i++) {
    g13[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    g23[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Evaluate the strain
  TacsScalar strain[8];
  evalStrain(strain, Ur, dr, Xdinv, zXdinv, T);
  addTyingStrain(strain, N13, N23, g13, g23, Xdinv, T);

  // Set a random perturbation vector
  TacsScalar scale = 1.25;
  TacsScalar eSens[8];
  for (int i = 0; i < 8; i++) {
    eSens[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  TacsScalar fval = 0.0;
  for (int j = 0; j < 8; j++) {
    fval += scale * strain[j] * eSens[j];
  }

  // Compute the derivative of the strain w.r.t. each
  // input component of the strain expression
  TacsScalar Urd[9], drd[9], Xdinvd[9], zXdinvd[9], Td[9];
  TacsScalar g13d[6], g23d[6];
  memset(g13d, 0, 6 * sizeof(TacsScalar));
  memset(g23d, 0, 6 * sizeof(TacsScalar));
  evalStrainSens(Urd, drd, Xdinvd, zXdinvd, Td, scale, eSens, Ur, dr, Xdinv,
                 zXdinv, T);
  addTyingStrainSens(g13d, g23d, Xdinvd, Td, scale, eSens, N13, N23, g13, g23,
                     Xdinv, T);

  // Compute the derivative using finite-difference
  TacsScalar e[8];

  // Compute the derivative w.r.t. Ur
  TacsScalar fd[9];
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = Ur[i];
#ifdef TACS_USE_COMPLEX
    Ur[i] = Ur[i] + TacsScalar(0.0, dh);
#else
    Ur[i] = Ur[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Ur[i] = tmp;
  }
  writeErrorComponents(stdout, "Urd", Urd, fd, 9);

  // Compute the derivative w.r.t. drd
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = dr[i];
#ifdef TACS_USE_COMPLEX
    dr[i] = dr[i] + TacsScalar(0.0, dh);
#else
    dr[i] = dr[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dr[i] = tmp;
  }
  writeErrorComponents(stdout, "drd", drd, fd, 9);

  // Compute the derivative w.r.t. Xdinv
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = Xdinv[i];
#ifdef TACS_USE_COMPLEX
    Xdinv[i] = Xdinv[i] + TacsScalar(0.0, dh);
#else
    Xdinv[i] = Xdinv[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xdinv[i] = tmp;
  }
  writeErrorComponents(stdout, "Xdinvd", Xdinvd, fd, 9);

  // Compute the derivative w.r.t. zXdinvd
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = zXdinv[i];
#ifdef TACS_USE_COMPLEX
    zXdinv[i] = zXdinv[i] + TacsScalar(0.0, dh);
#else
    zXdinv[i] = zXdinv[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    zXdinv[i] = tmp;
  }
  writeErrorComponents(stdout, "zXdinvd", zXdinvd, fd, 9);

  // Compute the derivative w.r.t. Td
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = T[i];
#ifdef TACS_USE_COMPLEX
    T[i] = T[i] + TacsScalar(0.0, dh);
#else
    T[i] = T[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    T[i] = tmp;
  }
  writeErrorComponents(stdout, "Td", Td, fd, 9);

  // Compute the derivative w.r.t. g13/g23
  for (int i = 0; i < 6; i++) {
    TacsScalar tmp = g13[i];
#ifdef TACS_USE_COMPLEX
    g13[i] = g13[i] + TacsScalar(0.0, dh);
#else
    g13[i] = g13[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    g13[i] = tmp;
  }
  writeErrorComponents(stdout, "g13d", g13d, fd, 6);

  for (int i = 0; i < 6; i++) {
    TacsScalar tmp = g23[i];
#ifdef TACS_USE_COMPLEX
    g23[i] = g23[i] + TacsScalar(0.0, dh);
#else
    g23[i] = g23[i] + dh;
#endif
    evalStrain(e, Ur, dr, Xdinv, zXdinv, T);
    addTyingStrain(e, N13, N23, g13, g23, Xdinv, T);
    fd[i] = 0.0;
    for (int j = 0; j < 8; j++) {
      fd[i] += scale * e[j] * eSens[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    g23[i] = tmp;
  }
  writeErrorComponents(stdout, "g23d", g23d, fd, 6);
}

/*
  Test the derivative of the transformation to the locally-aligned
  strain axis
*/
void MITC9::testTransformSens(double dh) {
  // Compute the transform as a function of Xa/Xb
  TacsScalar Xa[3], Xb[3];
  for (int i = 0; i < 3; i++) {
    Xa[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Xb[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  TacsScalar T[9], Td[9];
  computeTransform(T, Xa, Xb);
  for (int i = 0; i < 3; i++) {
    Td[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute a linear function of the transformation matrix components
  TacsScalar fval = 0.0;
  for (int i = 0; i < 9; i++) {
    fval += T[i] * Td[i];
  }

  // Compute the derivative of the transformation matrix
  TacsScalar Xad[3], Xbd[3];
  computeTransformSens(Xad, Xbd, Td, Xa, Xb);

  TacsScalar fd[3];
  for (int i = 0; i < 3; i++) {
    TacsScalar tmp = Xa[i];
#ifdef TACS_USE_COMPLEX
    Xa[i] = Xa[i] + TacsScalar(0.0, dh);
#else
    Xa[i] = Xa[i] + dh;
#endif
    TacsScalar T2[9];
    computeTransform(T2, Xa, Xb);
    fd[i] = 0.0;
    for (int j = 0; j < 9; j++) {
      fd[i] += T2[j] * Td[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xa[i] = tmp;
  }
  writeErrorComponents(stdout, "Xad", Xad, fd, 3);

  for (int i = 0; i < 3; i++) {
    TacsScalar tmp = Xb[i];
#ifdef TACS_USE_COMPLEX
    Xb[i] = Xb[i] + TacsScalar(0.0, dh);
#else
    Xb[i] = Xb[i] + dh;
#endif
    TacsScalar T2[9];
    computeTransform(T2, Xa, Xb);
    fd[i] = 0.0;
    for (int j = 0; j < 9; j++) {
      fd[i] += T2[j] * Td[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xb[i] = tmp;
  }
  writeErrorComponents(stdout, "Xbd", Xbd, fd, 3);
}

/*
  Test the sensitivity of the tying strain
*/
void MITC9::testTyingStrainSens(double dh) {
  // Set random points/variables for testing
  TacsScalar X[3 * NUM_NODES];
  TacsScalar vars[8 * NUM_NODES];
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    X[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 8 * NUM_NODES; i++) {
    vars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Set the derivatives
  TacsScalar g13d[6], g23d[6];
  for (int i = 0; i < 6; i++) {
    g13d[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    g23d[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute the reference frames at the nodes
  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  // Compute the directors at the nodes
  TacsScalar dir[3 * NUM_NODES];
  computeDirectors(dir, vars, Xr);

  // Compute the tensorial shear strain at the tying points
  TacsScalar g13[6], g23[6];
  computeTyingStrain(g13, g23, X, Xr, vars, dir);

  TacsScalar fval = 0.0;
  for (int i = 0; i < 6; i++) {
    fval += g13[i] * g13d[i] + g23[i] * g23d[i];
  }

  // The output from the computations
  TacsScalar dird[3 * NUM_NODES];
  TacsScalar Xrd[9 * NUM_NODES];
  TacsScalar Xd[3 * NUM_NODES];
  memset(dird, 0, 3 * NUM_NODES * sizeof(TacsScalar));
  memset(Xrd, 0, 9 * NUM_NODES * sizeof(TacsScalar));
  memset(Xd, 0, 3 * NUM_NODES * sizeof(TacsScalar));

  // input: g13d, g23d,
  addComputeTyingStrainSens(Xd, Xrd, dird, g13d, g23d, X, Xr, vars, dir);

  // Compute the derivative w.r.t. Ur
  TacsScalar fd[9 * NUM_NODES];
  for (int i = 0; i < 9 * NUM_NODES; i++) {
    TacsScalar tmp = Xr[i];
#ifdef TACS_USE_COMPLEX
    Xr[i] = Xr[i] + TacsScalar(0.0, dh);
#else
    Xr[i] = Xr[i] + dh;
#endif
    computeTyingStrain(g13, g23, X, Xr, vars, dir);
    fd[i] = 0.0;
    for (int j = 0; j < 6; j++) {
      fd[i] += g13[j] * g13d[j] + g23[j] * g23d[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xr[i] = tmp;
  }
  writeErrorComponents(stdout, "Xrd", Xrd, fd, 9 * NUM_NODES);

  // Compute the derivative w.r.t. drd
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    TacsScalar tmp = dir[i];
#ifdef TACS_USE_COMPLEX
    dir[i] = dir[i] + TacsScalar(0.0, dh);
#else
    dir[i] = dir[i] + dh;
#endif
    computeTyingStrain(g13, g23, X, Xr, vars, dir);
    fd[i] = 0.0;
    for (int j = 0; j < 6; j++) {
      fd[i] += g13[j] * g13d[j] + g23[j] * g23d[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dir[i] = tmp;
  }
  writeErrorComponents(stdout, "dird", dird, fd, 3 * NUM_NODES);

  // Compute the derivative w.r.t. Xdinv
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    TacsScalar tmp = X[i];
#ifdef TACS_USE_COMPLEX
    X[i] = X[i] + TacsScalar(0.0, dh);
#else
    X[i] = X[i] + dh;
#endif
    computeTyingStrain(g13, g23, X, Xr, vars, dir);
    fd[i] = 0.0;
    for (int j = 0; j < 6; j++) {
      fd[i] += g13[j] * g13d[j] + g23[j] * g23d[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    X[i] = tmp;
  }
  writeErrorComponents(stdout, "Xd", Xd, fd, 3 * NUM_NODES);
}

/*
  Test the implementation of the derivative of the through-thickness
  derivative d(Xdinv)/dz
*/
void MITC9::testNormalRateSens(double dh) {
  double Na[9], Nb[9];
  computeShapeFunc(-0.231, -0.734, Na, Nb);

  TacsScalar Xdinv[9];
  for (int i = 0; i < 9; i++) {
    Xdinv[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  TacsScalar Xr[9 * NUM_NODES];
  for (int i = 0; i < 9 * NUM_NODES; i++) {
    Xr[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute the normal rate matrix
  TacsScalar zXdinv[9];
  computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);

  TacsScalar Xdinvd[9], zXdinvd[9];
  for (int i = 0; i < 9; i++) {
    zXdinvd[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Xdinvd[i] = 0.0;
  }

  TacsScalar Xrd[9 * NUM_NODES];
  memset(Xrd, 0, 9 * NUM_NODES * sizeof(TacsScalar));
  addNormalRateMatSens(Xrd, Xdinvd, zXdinvd, Na, Nb, Xr, Xdinv);

  TacsScalar fval = 0.0;
  for (int j = 0; j < 9; j++) {
    fval += zXdinvd[j] * zXdinv[j];
  }

  // Compute the derivative of a function of zXdinvd
  // w.r.t. the components of Xr
  TacsScalar fd[9 * NUM_NODES];

  for (int i = 0; i < 9 * NUM_NODES; i++) {
    TacsScalar tmp = Xr[i];
#ifdef TACS_USE_COMPLEX
    Xr[i] = tmp + TacsScalar(0.0, dh);
#else
    Xr[i] = tmp + dh;
#endif
    computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);
    fd[i] = 0.0;
    for (int j = 0; j < 9; j++) {
      fd[i] += zXdinvd[j] * zXdinv[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xr[i] = tmp;
  }
  writeErrorComponents(stdout, "Xrd", Xrd, fd, 9 * NUM_NODES);

  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = Xdinv[i];
#ifdef TACS_USE_COMPLEX
    Xdinv[i] = tmp + TacsScalar(0.0, dh);
#else
    Xdinv[i] = tmp + dh;
#endif
    computeNormalRateMat(Na, Nb, Xr, Xdinv, zXdinv);
    fd[i] = 0.0;
    for (int j = 0; j < 9; j++) {
      fd[i] += zXdinvd[j] * zXdinv[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xdinv[i] = tmp;
  }
  writeErrorComponents(stdout, "Xdinv", Xdinvd, fd, 9);
}

/*
  Test the derivative of the Bmat computation with respect to the
  input parameters
*/
void MITC9::testBmatSens(double dh) {
  // Select a random point
  double pt[2] = {0.234, -0.783};

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(pt[0], pt[1], N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(pt[0], pt[1], Na, Nb);

  // Evaluate the tying strain interpolation
  double N13[6], N23[6];
  computeTyingFunc(pt[0], pt[1], N13, N23);

  // Inputs required for the bmatrix
  TacsScalar Ur[9], dr[9];
  TacsScalar Xdinv[9], zXdinv[9];
  TacsScalar T[9];

  // Assign random inputs
  for (int i = 0; i < 9; i++) {
    Ur[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    dr[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Xdinv[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    zXdinv[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    T[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute random values for the derivatives of the director
  TacsScalar dirdq[12 * NUM_NODES];
  for (int i = 0; i < 12 * NUM_NODES; i++) {
    dirdq[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Set random tying strain values
  TacsScalar B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  for (int i = 0; i < 48 * NUM_NODES; i++) {
    B13[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    B23[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Allocate random input vectors
  TacsScalar eSens[8], psi[8 * NUM_NODES];
  for (int i = 0; i < 8; i++) {
    eSens[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 8 * NUM_NODES; i++) {
    psi[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Evaluate the bmat/strain
  TacsScalar e[8], B[64 * NUM_NODES];
  evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
  addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);

  TacsScalar fval = 0.0;
  for (int k = 0; k < 8 * NUM_NODES; k++) {
    for (int j = 0; j < 8; j++) {
      fval += eSens[j] * B[j + 8 * k] * psi[k];
    }
  }

  // Compute the derivative of the strain w.r.t. each
  // input component of the strain expression
  TacsScalar Urd[9], drd[9], Xdinvd[9], zXdinvd[9], Td[9];
  TacsScalar dirdqd[12 * NUM_NODES];
  memset(Urd, 0, 9 * sizeof(TacsScalar));
  memset(drd, 0, 9 * sizeof(TacsScalar));
  memset(Xdinvd, 0, 9 * sizeof(TacsScalar));
  memset(zXdinvd, 0, 9 * sizeof(TacsScalar));
  memset(Td, 0, 9 * sizeof(TacsScalar));
  memset(dirdqd, 0, 12 * NUM_NODES * sizeof(TacsScalar));
  addBmatSens(Urd, drd, Xdinvd, zXdinvd, Td, dirdqd, eSens, psi, N, Na, Nb, Ur,
              dr, Xdinv, zXdinv, T, dirdq);

  // Add the derivative of the tying components of strain
  TacsScalar B13d[6 * 8 * NUM_NODES], B23d[6 * 8 * NUM_NODES];
  memset(B13d, 0, 6 * 8 * NUM_NODES * sizeof(TacsScalar));
  memset(B23d, 0, 6 * 8 * NUM_NODES * sizeof(TacsScalar));
  addTyingBmatSens(B13d, B23d, Xdinvd, Td, eSens, psi, N13, N23, B13, B23,
                   Xdinv, T);

  // Compute the derivative w.r.t. Ur
  TacsScalar fd[6 * 8 * NUM_NODES];
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = Ur[i];
#ifdef TACS_USE_COMPLEX
    Ur[i] = Ur[i] + TacsScalar(0.0, dh);
#else
    Ur[i] = Ur[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Ur[i] = tmp;
  }
  writeErrorComponents(stdout, "Urd", Urd, fd, 9);

  // Compute the derivative w.r.t. drd
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = dr[i];
#ifdef TACS_USE_COMPLEX
    dr[i] = dr[i] + TacsScalar(0.0, dh);
#else
    dr[i] = dr[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dr[i] = tmp;
  }
  writeErrorComponents(stdout, "drd", drd, fd, 9);

  // Compute the derivative w.r.t. Xdinv
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = Xdinv[i];
#ifdef TACS_USE_COMPLEX
    Xdinv[i] = Xdinv[i] + TacsScalar(0.0, dh);
#else
    Xdinv[i] = Xdinv[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xdinv[i] = tmp;
  }
  writeErrorComponents(stdout, "Xdinvd", Xdinvd, fd, 9);

  // Compute the derivative w.r.t. zXdinvd
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = zXdinv[i];
#ifdef TACS_USE_COMPLEX
    zXdinv[i] = zXdinv[i] + TacsScalar(0.0, dh);
#else
    zXdinv[i] = zXdinv[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    zXdinv[i] = tmp;
  }
  writeErrorComponents(stdout, "zXdinvd", zXdinvd, fd, 9);

  // Compute the derivative w.r.t. Td
  for (int i = 0; i < 9; i++) {
    TacsScalar tmp = T[i];
#ifdef TACS_USE_COMPLEX
    T[i] = T[i] + TacsScalar(0.0, dh);
#else
    T[i] = T[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    T[i] = tmp;
  }
  writeErrorComponents(stdout, "Td", Td, fd, 9);

  // Compute the derivative w.r.t. dirdq
  for (int i = 0; i < 12 * NUM_NODES; i++) {
    TacsScalar tmp = dirdq[i];
#ifdef TACS_USE_COMPLEX
    dirdq[i] = dirdq[i] + TacsScalar(0.0, dh);
#else
    dirdq[i] = dirdq[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dirdq[i] = tmp;
  }
  writeErrorComponents(stdout, "dirdq", dirdqd, fd, 12 * NUM_NODES);

  for (int i = 0; i < 6 * 8 * NUM_NODES; i++) {
    TacsScalar tmp = B13[i];
#ifdef TACS_USE_COMPLEX
    B13[i] = B13[i] + TacsScalar(0.0, dh);
#else
    B13[i] = B13[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    B13[i] = tmp;
  }
  writeErrorComponents(stdout, "B13", B13d, fd, 6 * 8 * NUM_NODES);

  for (int i = 0; i < 6 * 8 * NUM_NODES; i++) {
    TacsScalar tmp = B23[i];
#ifdef TACS_USE_COMPLEX
    B23[i] = B23[i] + TacsScalar(0.0, dh);
#else
    B23[i] = B23[i] + dh;
#endif
    evalBmat(e, B, N, Na, Nb, Ur, dr, Xdinv, zXdinv, T, dirdq);
    addTyingBmat(B, N13, N23, B13, B23, Xdinv, T);
    fd[i] = 0.0;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      for (int j = 0; j < 8; j++) {
        fd[i] += eSens[j] * B[j + 8 * k] * psi[k];
      }
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    B23[i] = tmp;
  }
  writeErrorComponents(stdout, "B23", B23d, fd, 6 * 8 * NUM_NODES);
}

/*
  Test the derivative of the tying strain with respect to the nodes
*/
void MITC9::testTyingBmatSens(double dh) {
  // Assign random inputs
  TacsScalar vars[8 * NUM_NODES];
  TacsScalar X[3 * NUM_NODES], Xr[9 * NUM_NODES];
  TacsScalar dir[3 * NUM_NODES], dirdq[12 * NUM_NODES];

  // Assign random inputs
  for (int i = 0; i < 8 * NUM_NODES; i++) {
    vars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    X[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    dir[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 9 * NUM_NODES; i++) {
    Xr[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 12 * NUM_NODES; i++) {
    dirdq[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Compute random input perturbations
  TacsScalar g13d[6], g23d[6];
  for (int i = 0; i < 6; i++) {
    g13d[i] = 1.0 * rand() / RAND_MAX;
    g23d[i] = 1.0 * rand() / RAND_MAX;
  }

  TacsScalar B13d[6 * 8 * NUM_NODES], B23d[6 * 8 * NUM_NODES];
  for (int i = 0; i < 6 * 8 * NUM_NODES; i++) {
    B13d[i] = 1.0 * rand() / RAND_MAX;
    B23d[i] = 1.0 * rand() / RAND_MAX;
  }

  // Compute the derivatives
  TacsScalar Xd[3 * NUM_NODES], Xrd[9 * NUM_NODES];
  TacsScalar dird[3 * NUM_NODES], dirdqd[12 * NUM_NODES];
  memset(Xd, 0, sizeof(Xd));
  memset(Xrd, 0, sizeof(Xrd));
  memset(dird, 0, sizeof(dird));
  memset(dirdqd, 0, sizeof(dirdqd));

  addComputeTyingBmatSens(Xd, Xrd, dird, dirdqd, g13d, g23d, B13d, B23d, X, Xr,
                          vars, dir, dirdq);

  // Evaluate the function
  TacsScalar g13[6], g23[6], B13[6 * 8 * NUM_NODES], B23[6 * 8 * NUM_NODES];
  computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);

  TacsScalar fval = 0.0;
  for (int k = 0; k < 6; k++) {
    fval += g13d[k] * g13[k] + g23d[k] * g23[k];
  }
  for (int k = 0; k < 6 * 8 * NUM_NODES; k++) {
    fval += B13d[k] * B13[k] + B23d[k] * B23[k];
  }

  TacsScalar fd[12 * NUM_NODES];
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    TacsScalar tmp = X[i];
#ifdef TACS_USE_COMPLEX
    X[i] = X[i] + TacsScalar(0.0, dh);
#else
    X[i] = X[i] + dh;
#endif
    computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
    fd[i] = 0.0;
    for (int k = 0; k < 6; k++) {
      fd[i] += g13d[k] * g13[k] + g23d[k] * g23[k];
    }
    for (int k = 0; k < 6 * 8 * NUM_NODES; k++) {
      fd[i] += B13d[k] * B13[k] + B23d[k] * B23[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    X[i] = tmp;
  }
  writeErrorComponents(stdout, "Xd", Xd, fd, 3 * NUM_NODES);

  for (int i = 0; i < 9 * NUM_NODES; i++) {
    TacsScalar tmp = Xr[i];
#ifdef TACS_USE_COMPLEX
    Xr[i] = Xr[i] + TacsScalar(0.0, dh);
#else
    Xr[i] = Xr[i] + dh;
#endif
    computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
    fd[i] = 0.0;
    for (int k = 0; k < 6; k++) {
      fd[i] += g13d[k] * g13[k] + g23d[k] * g23[k];
    }
    for (int k = 0; k < 6 * 8 * NUM_NODES; k++) {
      fd[i] += B13d[k] * B13[k] + B23d[k] * B23[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xr[i] = tmp;
  }
  writeErrorComponents(stdout, "Xrd", Xrd, fd, 9 * NUM_NODES);

  for (int i = 0; i < 3 * NUM_NODES; i++) {
    TacsScalar tmp = dir[i];
#ifdef TACS_USE_COMPLEX
    dir[i] = dir[i] + TacsScalar(0.0, dh);
#else
    dir[i] = dir[i] + dh;
#endif
    computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
    fd[i] = 0.0;
    for (int k = 0; k < 6; k++) {
      fd[i] += g13d[k] * g13[k] + g23d[k] * g23[k];
    }
    for (int k = 0; k < 6 * 8 * NUM_NODES; k++) {
      fd[i] += B13d[k] * B13[k] + B23d[k] * B23[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dir[i] = tmp;
  }
  writeErrorComponents(stdout, "dird", dird, fd, 3 * NUM_NODES);

  for (int i = 0; i < 12 * NUM_NODES; i++) {
    TacsScalar tmp = dirdq[i];
#ifdef TACS_USE_COMPLEX
    dirdq[i] = dirdq[i] + TacsScalar(0.0, dh);
#else
    dirdq[i] = dirdq[i] + dh;
#endif
    computeTyingBmat(g13, g23, B13, B23, X, Xr, vars, dir, dirdq);
    fd[i] = 0.0;
    for (int k = 0; k < 6; k++) {
      fd[i] += g13d[k] * g13[k] + g23d[k] * g23[k];
    }
    for (int k = 0; k < 6 * 8 * NUM_NODES; k++) {
      fd[i] += B13d[k] * B13[k] + B23d[k] * B23[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    dirdq[i] = tmp;
  }
  writeErrorComponents(stdout, "dirdqd", dirdqd, fd, 12 * NUM_NODES);
}

/*
  Test the rotational penalty derivative
*/
void MITC9::testBrotSens(double dh) {
  // Assign random inputs
  TacsScalar vars[8 * NUM_NODES];
  TacsScalar Xa[3], Xb[3], Ua[3], Ub[3];

  // Assign random inputs
  for (int i = 0; i < 8 * NUM_NODES; i++) {
    vars[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 3; i++) {
    Xa[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Xb[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Ua[i] = -1.0 + 2.0 * rand() / RAND_MAX;
    Ub[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Select a random point
  double pt[2] = {0.763, -0.314};

  // Evaluate the shape functions
  double N[NUM_NODES];
  computeShapeFunc(pt[0], pt[1], N);

  // Evaluate the derivatives of the shape functions
  double Na[NUM_NODES], Nb[NUM_NODES];
  computeShapeFunc(pt[0], pt[1], Na, Nb);

  // Set up the input parameters
  TacsScalar scale = 1.0 * rand() / RAND_MAX;
  TacsScalar rotd = rand() / RAND_MAX;
  TacsScalar psi[8 * NUM_NODES];
  for (int i = 0; i < 8 * NUM_NODES; i++) {
    psi[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  // Evaluate the derivative
  TacsScalar Xad[3], Xbd[3];
  memset(Xad, 0, 3 * sizeof(TacsScalar));
  memset(Xbd, 0, 3 * sizeof(TacsScalar));
  addBRotPenaltySens(Xad, Xbd, rotd, scale, psi, N, Na, Nb, Xa, Xb, Ua, Ub,
                     vars);

  // Compute the shape function and its derivative
  TacsScalar rot;
  TacsScalar Brot[8 * NUM_NODES];
  rot = computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);

  TacsScalar fval = scale * rot * rotd;
  for (int k = 0; k < 8 * NUM_NODES; k++) {
    fval += scale * Brot[k] * psi[k];
  }

  // Compute the finite-difference approximations
  TacsScalar fd[3];
  for (int i = 0; i < 3; i++) {
    TacsScalar tmp = Xa[i];
#ifdef TACS_USE_COMPLEX
    Xa[i] = Xa[i] + TacsScalar(0.0, dh);
#else
    Xa[i] = Xa[i] + dh;
#endif
    rot = computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);
    fd[i] = scale * rot * rotd;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      fd[i] += scale * Brot[k] * psi[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xa[i] = tmp;
  }
  writeErrorComponents(stdout, "Xad", Xad, fd, 3);

  for (int i = 0; i < 3; i++) {
    TacsScalar tmp = Xb[i];
#ifdef TACS_USE_COMPLEX
    Xb[i] = Xb[i] + TacsScalar(0.0, dh);
#else
    Xb[i] = Xb[i] + dh;
#endif
    rot = computeBRotPenalty(Brot, N, Na, Nb, Xa, Xb, Ua, Ub, vars);
    fd[i] = scale * rot * rotd;
    for (int k = 0; k < 8 * NUM_NODES; k++) {
      fd[i] += scale * Brot[k] * psi[k];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    Xb[i] = tmp;
  }
  writeErrorComponents(stdout, "Xbd", Xbd, fd, 3);
}

/*
  Test the sensitivity of the frame assembly code
*/
void MITC9::testFrameSens(double dh) {
  TacsScalar X[3 * NUM_NODES], Xrd[9 * NUM_NODES];
  for (int i = 0; i < 3 * NUM_NODES; i++) {
    X[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }
  for (int i = 0; i < 9 * NUM_NODES; i++) {
    Xrd[i] = -1.0 + 2.0 * rand() / RAND_MAX;
  }

  TacsScalar Xr[9 * NUM_NODES];
  computeFrames(Xr, X);

  TacsScalar fval = 0.0;
  for (int j = 0; j < 9 * NUM_NODES; j++) {
    fval += Xr[j] * Xrd[j];
  }

  TacsScalar Xd[3 * NUM_NODES];
  memset(Xd, 0, 3 * NUM_NODES * sizeof(TacsScalar));
  addFramesSens(Xd, Xrd, X);

  // Compute the derivative of a function of zXdinvd
  // w.r.t. the components of Xr
  TacsScalar fd[3 * NUM_NODES];

  for (int i = 0; i < 3 * NUM_NODES; i++) {
    TacsScalar tmp = X[i];
#ifdef TACS_USE_COMPLEX
    X[i] = tmp + TacsScalar(0.0, dh);
#else
    X[i] = tmp + dh;
#endif
    computeFrames(Xr, X);
    fd[i] = 0.0;
    for (int j = 0; j < 9 * NUM_NODES; j++) {
      fd[i] += Xr[j] * Xrd[j];
    }
#ifdef TACS_USE_COMPLEX
    fd[i] = TacsImagPart(fd[i]) / dh;
#else
    fd[i] = (fd[i] - fval) / dh;
#endif
    X[i] = tmp;
  }
  writeErrorComponents(stdout, "Xd", Xd, fd, 3 * NUM_NODES);
}
