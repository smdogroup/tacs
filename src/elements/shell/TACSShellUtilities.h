#ifndef TACS_SHELL_UTILITIES_H
#define TACS_SHELL_UTILITIES_H

#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"
#include "TACSShellElementTransform.h"

inline void TacsShellAssembleFrame(const TacsScalar Xxi[], const TacsScalar n[],
                                   TacsScalar Xd[]) {
  Xd[0] = Xxi[0];
  Xd[1] = Xxi[1];
  Xd[2] = n[0];

  Xd[3] = Xxi[2];
  Xd[4] = Xxi[3];
  Xd[5] = n[1];

  Xd[6] = Xxi[4];
  Xd[7] = Xxi[5];
  Xd[8] = n[2];
}

inline void TacsShellAssembleFrame(const TacsScalar nxi[], TacsScalar Xdz[]) {
  Xdz[0] = nxi[0];
  Xdz[1] = nxi[1];
  Xdz[2] = 0.0;

  Xdz[3] = nxi[2];
  Xdz[4] = nxi[3];
  Xdz[5] = 0.0;

  Xdz[6] = nxi[4];
  Xdz[7] = nxi[5];
  Xdz[8] = 0.0;
}

inline void TacsShellAssembleFrame(const TacsScalar a[], const TacsScalar b[],
                                   const TacsScalar c[], TacsScalar Xd[]) {
  if (a) {
    Xd[0] = a[0];
    Xd[3] = a[1];
    Xd[6] = a[2];
  } else {
    Xd[0] = 0.0;
    Xd[3] = 0.0;
    Xd[6] = 0.0;
  }

  if (b) {
    Xd[1] = b[0];
    Xd[4] = b[1];
    Xd[7] = b[2];
  } else {
    Xd[1] = 0.0;
    Xd[4] = 0.0;
    Xd[7] = 0.0;
  }

  if (c) {
    Xd[2] = c[0];
    Xd[5] = c[1];
    Xd[8] = c[2];
  } else {
    Xd[2] = 0.0;
    Xd[5] = 0.0;
    Xd[8] = 0.0;
  }
}

inline void TacsShellExtractFrame(const TacsScalar Xd[], TacsScalar Xxi[],
                                  TacsScalar n[]) {
  Xxi[0] = Xd[0];
  Xxi[1] = Xd[1];
  n[0] = Xd[2];

  Xxi[2] = Xd[3];
  Xxi[3] = Xd[4];
  n[1] = Xd[5];

  Xxi[4] = Xd[6];
  Xxi[5] = Xd[7];
  n[2] = Xd[8];
}

inline void TacsShellExtractFrame(const TacsScalar Xd[], TacsScalar Xxi[]) {
  Xxi[0] = Xd[0];
  Xxi[1] = Xd[1];

  Xxi[2] = Xd[3];
  Xxi[3] = Xd[4];

  Xxi[4] = Xd[6];
  Xxi[5] = Xd[7];
}

inline void TacsShellExtractFrameSens(const TacsScalar d2u0d[],
                                      TacsScalar d2u0xi[]) {
  // Extract the second derivatives
  for (int j = 0; j < 6; j++) {
    int jj = 3 * (j / 2) + (j % 2);
    for (int i = 0; i < 6; i++) {
      int ii = 3 * (i / 2) + (i % 2);

      d2u0xi[6 * j + i] = d2u0d[9 * jj + ii];
    }
  }
}

inline void TacsShellExtractFrameSens(const TacsScalar d2u0d[],
                                      TacsScalar d2u0xi[], TacsScalar d2d0[],
                                      TacsScalar d2d0u0xi[]) {
  // Extract the second derivatives
  for (int j = 0; j < 6; j++) {
    int jj = 3 * (j / 2) + (j % 2);
    for (int i = 0; i < 6; i++) {
      int ii = 3 * (i / 2) + (i % 2);

      d2u0xi[6 * j + i] = d2u0d[9 * jj + ii];
    }
  }

  for (int j = 0; j < 3; j++) {
    int jj = 3 * j + 2;
    for (int i = 0; i < 3; i++) {
      int ii = 3 * i + 2;
      d2d0[3 * j + i] = d2u0d[9 * jj + ii];
    }
  }

  for (int j = 0; j < 3; j++) {
    int jj = 3 * j + 2;
    for (int i = 0; i < 6; i++) {
      int ii = 3 * (i / 2) + (i % 2);
      d2d0u0xi[6 * j + i] = d2u0d[9 * jj + ii];
    }
  }
}

inline void TacsShellExtractFrameMixedSens(const TacsScalar d2u0du1d[],
                                           TacsScalar d2d0xiu0xi[],
                                           TacsScalar d2d0d0xi[]) {
  // Extract the second derivatives
  for (int j = 0; j < 6; j++) {
    int jj = 3 * (j / 2) + (j % 2);
    for (int i = 0; i < 6; i++) {
      int ii = 3 * (i / 2) + (i % 2);

      d2d0xiu0xi[6 * i + j] = d2u0du1d[9 * jj + ii];
    }
  }

  for (int j = 0; j < 3; j++) {
    int jj = 3 * j + 2;
    for (int i = 0; i < 6; i++) {
      int ii = 3 * (i / 2) + (i % 2);
      d2d0d0xi[6 * j + i] = d2u0du1d[9 * jj + ii];
    }
  }
}

/*
  Compute the second derivative of the transformation

  C = T^{T}*A*B
  dA = T*dC*B^{T}
*/
inline void mat3x3TransMatMatHessian(const TacsScalar T[], const TacsScalar B[],
                                     const TacsScalar d2C[], TacsScalar d2A[]) {
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for (int i = 0; i < 9; i++) {
    mat3x3MatTransMult(&d2C[9 * i], B, t1);
    mat3x3MatMult(T, t1, &tmp2[9 * i]);
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      t1[j] = tmp2[i + 9 * j];
    }

    mat3x3MatTransMult(t1, B, t2);
    mat3x3MatMult(T, t2, t1);

    for (int j = 0; j < 9; j++) {
      d2A[i + 9 * j] = t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAdd(const TacsScalar T[],
                                        const TacsScalar B[],
                                        const TacsScalar d2C[],
                                        TacsScalar d2A[]) {
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for (int i = 0; i < 9; i++) {
    mat3x3MatTransMult(&d2C[9 * i], B, t1);
    mat3x3MatMult(T, t1, &tmp2[9 * i]);
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      t1[j] = tmp2[i + 9 * j];
    }

    mat3x3MatTransMult(t1, B, t2);
    mat3x3MatMult(T, t2, t1);

    for (int j = 0; j < 9; j++) {
      d2A[i + 9 * j] += t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessian(const TacsScalar T1[],
                                     const TacsScalar T2[],
                                     const TacsScalar B1[],
                                     const TacsScalar B2[],
                                     const TacsScalar d2C[], TacsScalar d2A[]) {
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for (int i = 0; i < 9; i++) {
    mat3x3MatTransMult(&d2C[9 * i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9 * i]);
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      t1[j] = tmp2[i + 9 * j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for (int j = 0; j < 9; j++) {
      d2A[i + 9 * j] = t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAdd(
    const TacsScalar T1[], const TacsScalar T2[], const TacsScalar B1[],
    const TacsScalar B2[], const TacsScalar d2C[], TacsScalar d2A[]) {
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for (int i = 0; i < 9; i++) {
    mat3x3MatTransMult(&d2C[9 * i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9 * i]);
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      t1[j] = tmp2[i + 9 * j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for (int j = 0; j < 9; j++) {
      d2A[i + 9 * j] += t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAddSymm(
    const TacsScalar T1[], const TacsScalar T2[], const TacsScalar B1[],
    const TacsScalar B2[], const TacsScalar d2C[], TacsScalar d2A[]) {
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for (int i = 0; i < 9; i++) {
    mat3x3MatTransMult(&d2C[9 * i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9 * i]);
  }

  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      t1[j] = tmp2[i + 9 * j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for (int j = 0; j < 9; j++) {
      d2A[i + 9 * j] += t1[j];
      d2A[j + 9 * i] += t1[j];
    }
  }
}

/**
  Compute the frame normals at each of the nodes of the shell element

  This code computes the in-plane coordinates derivatives to evaluate
  the shell normal.

  @param Xpts The node locations for the shell element
  @param fn The frame normals computed at each node
  @param Xdn The derivatives at the node
  @param fnorm Optional: the norm of the cross-product
*/
template <class basis>
void TacsShellComputeNodeNormals(const TacsScalar Xpts[], TacsScalar fn[],
                                 TacsScalar Xdn[] = NULL,
                                 TacsScalar fnorm[] = NULL) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    TacsScalar Xxi[6];
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

    TacsScalar a[3], b[3];
    a[0] = Xxi[0];
    a[1] = Xxi[2];
    a[2] = Xxi[4];

    b[0] = Xxi[1];
    b[1] = Xxi[3];
    b[2] = Xxi[5];

    // Compute the normal direction at the point
    crossProduct(a, b, &fn[3 * i]);

    // Compute the 2-norm of the vector in the normal direction
    TacsScalar norm = sqrt(vec3Dot(&fn[3 * i], &fn[3 * i]));

    // Save the 2-norm value if the fnorm argument is not NULL
    if (fnorm) {
      fnorm[i] = norm;
    }

    // Scale the normal direction
    if (norm != 0.0) {
      vec3Scale(1.0 / norm, &fn[3 * i]);
    }

    if (Xdn) {
      TacsShellAssembleFrame(Xxi, &fn[3 * i], &Xdn[9 * i]);
    }
  }
}

/**
  Sensitivity of TacsShellComputeNodeNormals with respect to Xpts.

  Given a seed dfn on the per-node normalized frame normals fn computed by
  TacsShellComputeNodeNormals, accumulate dXpts += d(dfn^T*fn)/d(Xpts).

  @param Xpts The node locations for the shell element
  @param dfn The seed on the per-node normal directions, size 3*NUM_NODES
  @param dXpts The accumulated sensitivity with respect to Xpts
*/
template <class basis>
void TacsShellAddNodeNormalsSens(const TacsScalar Xpts[],
                                 const TacsScalar dfn[], TacsScalar dXpts[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Recompute the forward quantities needed for the reverse sweep
    TacsScalar Xxi[6];
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

    TacsScalar a[3], b[3];
    a[0] = Xxi[0];
    a[1] = Xxi[2];
    a[2] = Xxi[4];

    b[0] = Xxi[1];
    b[1] = Xxi[3];
    b[2] = Xxi[5];

    TacsScalar raw[3];
    crossProduct(a, b, raw);

    // Reverse sweep: fn = normalize(raw)
    TacsScalar draw[3] = {dfn[3 * i], dfn[3 * i + 1], dfn[3 * i + 2]};
    TacsScalar sAnrm;
    vec3NormalizeSens(raw, &sAnrm, draw);

    // raw = cross(a, b): da += cross(b, draw), db += cross(draw, a)
    TacsScalar da[3], db[3];
    crossProduct(b, draw, da);
    crossProduct(draw, a, db);

    TacsScalar dXxi[6];
    dXxi[0] = da[0];
    dXxi[2] = da[1];
    dXxi[4] = da[2];

    dXxi[1] = db[0];
    dXxi[3] = db[1];
    dXxi[5] = db[2];

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dXpts);
  }
}

/**
  Compute the displacement gradient of the constant and through-thickness
  rate of change of the displacements.

  @param pt The parametric point
  @param Xpts The node locations for the element
  @param vars The element variables
  @param fn The frame normal directions at each node
  @param d The director field at each node
  @param Xxi The in-plane coordinate derivatives
  @param n0 The interpolated frame normal direction
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param u0x Derivative of the displacement in the local x coordinates
  @param u1x Derivative of the through-thickness disp. in local x coordinates
*/
template <int vars_per_node, class basis>
TacsScalar TacsShellComputeDispGrad(const double pt[], const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar fn[], const TacsScalar d[],
                                    const TacsScalar Xxi[],
                                    const TacsScalar n0[], const TacsScalar T[],
                                    TacsScalar XdinvT[], TacsScalar XdinvzT[],
                                    TacsScalar u0x[], TacsScalar u1x[]) {
  // Compute n,xi = [dn/dxi1; dn/dxi2]
  TacsScalar nxi[6];
  basis::template interpFieldsGrad<3, 3>(pt, fn, nxi);

  // Assemble the terms Xd = [Xxi; n] and Xdz
  TacsScalar Xd[9], Xdz[9];
  TacsShellAssembleFrame(Xxi, n0, Xd);
  TacsShellAssembleFrame(nxi, Xdz);

  // Compute the inverse of the 3x3 Jacobian transformation
  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  // Compute negXdinvXdz = -Xdinv*Xdz
  TacsScalar negXdinvXdz[9];
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for (int i = 0; i < 9; i++) {
    negXdinvXdz[i] *= -1.0;
  }

  // Compute XdinvT = Xdinv*T
  mat3x3MatMult(Xdinv, T, XdinvT);

  // Compute Xdinvz = -Xdinv*Xdz*Xdinv*T
  mat3x3MatMult(negXdinvXdz, XdinvT, XdinvzT);

  // Compute the director field and the gradient of the director
  // field at the specified point
  TacsScalar d0[3], d0xi[6];
  basis::template interpFields<3, 3>(pt, d, d0);
  basis::template interpFieldsGrad<3, 3>(pt, d, d0xi);

  // Compute the gradient of the displacement solution at the quadrature points
  TacsScalar u0xi[6];
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);

  // Compute the derivative u0,x
  TacsShellAssembleFrame(u0xi, d0, u0x);  // Use u0x to store [u0,xi; d0]

  // u1x = T^{T}*u1d*XdinvT + T^{T}*u0d*XdinvzT
  TacsScalar tmp[9];
  TacsShellAssembleFrame(d0xi, u1x);  // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1x, XdinvT, tmp);
  mat3x3MatMultAdd(u0x, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1x);

  // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
  // u0x = T^{T}*u0d*Xdinv*T
  mat3x3MatMult(u0x, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0x);

  return detXd;
}

/**
  Compute the displacement gradient and the derivative of the displacement
  gradient with respect to input variables

  @param pt The parametric point
  @param Xpts The node locations for the element
  @param vars The element variables
  @param fn The frame normal directions at each node
  @param d The director field at each node
  @param Xxi The in-plane coordinate derivatives
  @param n0 The interpolated frame normal direction
  @param T The transformation to local coordinates
  @param varsd The derivative of the vars
  @param dd The derivative of the director field
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param u0x Derivative of the displacement in the local x coordinates
  @param u1x Derivative of the through-thickness disp. in local x coordinates
  @param u0xd Derivative of u0x
  @param u1xd Derivative of u1x
*/
template <int vars_per_node, class basis>
TacsScalar TacsShellComputeDispGradDeriv(
    const double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar fn[], const TacsScalar d[], const TacsScalar Xxi[],
    const TacsScalar n0[], const TacsScalar T[], const TacsScalar varsd[],
    const TacsScalar dd[], TacsScalar XdinvT[], TacsScalar XdinvzT[],
    TacsScalar u0x[], TacsScalar u1x[], TacsScalar u0xd[], TacsScalar u1xd[]) {
  // Compute n,xi = [dn/dxi1; dn/dxi2]
  TacsScalar nxi[6];
  basis::template interpFieldsGrad<3, 3>(pt, fn, nxi);

  // Assemble the terms Xd = [Xxi; n] and Xdz
  TacsScalar Xd[9], Xdz[9];
  TacsShellAssembleFrame(Xxi, n0, Xd);
  TacsShellAssembleFrame(nxi, Xdz);

  // Compute the inverse of the 3x3 Jacobian transformation
  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  // Compute negXdinvXdz = -Xdinv*Xdz
  TacsScalar negXdinvXdz[9];
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for (int i = 0; i < 9; i++) {
    negXdinvXdz[i] *= -1.0;
  }

  // Compute XdinvT = Xdinv*T
  mat3x3MatMult(Xdinv, T, XdinvT);

  // Compute Xdinvz = -Xdinv*Xdz*Xdinv*T
  mat3x3MatMult(negXdinvXdz, XdinvT, XdinvzT);

  // Compute the director field and the gradient of the director
  // field at the specified point
  TacsScalar d0[3], d0xi[6], d0d[3], d0xid[6];
  basis::template interpFields<3, 3>(pt, d, d0);
  basis::template interpFieldsGrad<3, 3>(pt, d, d0xi);
  basis::template interpFields<3, 3>(pt, dd, d0d);
  basis::template interpFieldsGrad<3, 3>(pt, dd, d0xid);

  // Compute the gradient of the displacement solution at the quadrature points
  TacsScalar u0xi[6], u0xid[6];
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, u0xid);

  // Compute the derivative u0,x
  TacsShellAssembleFrame(u0xi, d0, u0x);  // Use u0x to store [u0,xi; d0]
  TacsShellAssembleFrame(u0xid, d0d, u0xd);

  // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
  TacsScalar tmp[9];
  TacsShellAssembleFrame(d0xi, u1x);  // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1x, XdinvT, tmp);
  mat3x3MatMultAdd(u0x, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1x);

  TacsShellAssembleFrame(d0xid, u1xd);  // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1xd, XdinvT, tmp);
  mat3x3MatMultAdd(u0xd, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1xd);

  // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
  // u0x = T^{T}*u0d*Xdinv*T
  mat3x3MatMult(u0x, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0x);

  mat3x3MatMult(u0xd, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0xd);

  return detXd;
}

/**
  Sensitivity of TacsShellComputeDispGrad with respect to Xxi, n0, fn, T and
  d (the state-variable/vars direction is handled separately by
  TacsShellAddDispGradSens).

  Given adjoint seeds du0x/du1x on the two matrix outputs of
  TacsShellComputeDispGrad and ddetXd on its scalar detXd return value,
  accumulate the corresponding sensitivities into dXxi, dn0, dfn, dT, dd.

  @param pt The parametric point
  @param Xpts The node locations for the element (unused by the forward
              function's own body; kept for signature symmetry)
  @param vars The element variables (fixed; only used to rebuild u0xi)
  @param fn The frame normal directions at each node
  @param d The director field at each node
  @param Xxi The in-plane coordinate derivatives
  @param n0 The interpolated frame normal direction
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param du0x Seed on u0x
  @param du1x Seed on u1x
  @param ddetXd Seed on detXd
  @param dXxi Xxi sensitivity (accumulated)
  @param dn0 n0 sensitivity (accumulated)
  @param dfn Nodal fn sensitivity (accumulated)
  @param dT T sensitivity (accumulated)
  @param dd Nodal director-field sensitivity (accumulated)
*/
template <int vars_per_node, class basis>
void TacsShellComputeDispGradXptSens(
    const double pt[], const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar fn[], const TacsScalar d[], const TacsScalar Xxi[],
    const TacsScalar n0[], const TacsScalar T[], const TacsScalar XdinvT[],
    const TacsScalar XdinvzT[], const TacsScalar du0x[],
    const TacsScalar du1x[], const TacsScalar ddetXd, TacsScalar dXxi[],
    TacsScalar dn0[], TacsScalar dfn[], TacsScalar dT[], TacsScalar dd[]) {
  // Recompute the forward quantities needed for the reverse sweep (mirrors
  // TacsShellComputeDispGrad's own body)
  TacsScalar nxi[6];
  basis::template interpFieldsGrad<3, 3>(pt, fn, nxi);

  TacsScalar Xd[9], Xdz[9];
  TacsShellAssembleFrame(Xxi, n0, Xd);
  TacsShellAssembleFrame(nxi, Xdz);

  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  TacsScalar negXdinvXdz[9];
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for (int i = 0; i < 9; i++) {
    negXdinvXdz[i] *= -1.0;
  }

  TacsScalar d0[3], d0xi[6];
  basis::template interpFields<3, 3>(pt, d, d0);
  basis::template interpFieldsGrad<3, 3>(pt, d, d0xi);

  TacsScalar u0xi[6];
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);

  TacsScalar U[9], V[9];
  TacsShellAssembleFrame(u0xi, d0, U);
  TacsShellAssembleFrame(d0xi, V);

  TacsScalar tmp1[9], tmp2[9];
  mat3x3MatMult(V, XdinvT, tmp1);
  mat3x3MatMultAdd(U, XdinvzT, tmp1);
  mat3x3MatMult(U, XdinvT, tmp2);

  // ---- Reverse sweep ----
  TacsScalar dtmp1[9], dtmp2[9], tmp[9];

  // u0x = T^{T}*tmp2
  mat3x3MatMult(T, du0x, dtmp2);
  mat3x3MatTransMult(tmp2, du0x, tmp);
  for (int i = 0; i < 9; i++) dT[i] += tmp[i];

  // u1x = T^{T}*tmp1
  mat3x3MatMult(T, du1x, dtmp1);
  mat3x3MatTransMult(tmp1, du1x, tmp);
  for (int i = 0; i < 9; i++) dT[i] += tmp[i];

  // tmp2 = U*XdinvT
  TacsScalar dU[9], dV[9], dXdinvT[9], dXdinvzT[9];
  mat3x3MatTransMult(dtmp2, XdinvT, dU);
  mat3x3TransMatMult(U, dtmp2, dXdinvT);

  // tmp1 = V*XdinvT + U*XdinvzT
  mat3x3MatTransMult(dtmp1, XdinvT, dV);
  mat3x3TransMatMult(V, dtmp1, tmp);
  for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

  mat3x3MatTransMult(dtmp1, XdinvzT, tmp);
  for (int i = 0; i < 9; i++) dU[i] += tmp[i];
  mat3x3TransMatMult(U, dtmp1, dXdinvzT);

  // XdinvzT = negXdinvXdz*XdinvT
  TacsScalar dnegXdinvXdz[9];
  mat3x3MatTransMult(dXdinvzT, XdinvT, dnegXdinvXdz);
  mat3x3TransMatMult(negXdinvXdz, dXdinvzT, tmp);
  for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

  // XdinvT = Xdinv*T
  TacsScalar dXdinv[9];
  mat3x3MatTransMult(dXdinvT, T, dXdinv);
  mat3x3TransMatMult(Xdinv, dXdinvT, tmp);
  for (int i = 0; i < 9; i++) dT[i] += tmp[i];

  // negXdinvXdz = -Xdinv*Xdz => the seed on the unscaled product Xdinv*Xdz
  // is -dnegXdinvXdz
  TacsScalar dCprime[9], dXdz[9];
  for (int i = 0; i < 9; i++) dCprime[i] = -dnegXdinvXdz[i];
  mat3x3MatTransMult(dCprime, Xdz, tmp);
  for (int i = 0; i < 9; i++) dXdinv[i] += tmp[i];
  mat3x3TransMatMult(Xdinv, dCprime, dXdz);

  // detXd = det(Xd): dXd += ddetXd*detXd*Xdinv^{T}
  TacsScalar dXd[9];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dXd[3 * i + j] = ddetXd * detXd * Xdinv[3 * j + i];
    }
  }

  // Xdinv = inv3x3(Xd): dXd += -Xdinv^{T}*dXdinv*Xdinv^{T}
  TacsScalar dXd_inv[9];
  inv3x3Sens(Xdinv, dXdinv, dXd_inv);
  for (int i = 0; i < 9; i++) dXd[i] += dXd_inv[i];

  // Xd = assembleFrame(Xxi, n0)
  dXxi[0] += dXd[0];
  dXxi[1] += dXd[1];
  dXxi[2] += dXd[3];
  dXxi[3] += dXd[4];
  dXxi[4] += dXd[6];
  dXxi[5] += dXd[7];

  dn0[0] += dXd[2];
  dn0[1] += dXd[5];
  dn0[2] += dXd[8];

  // Xdz = assembleFrame(nxi) (no third column, no fn-direction term there)
  TacsScalar dnxi[6];
  dnxi[0] = dXdz[0];
  dnxi[1] = dXdz[1];
  dnxi[2] = dXdz[3];
  dnxi[3] = dXdz[4];
  dnxi[4] = dXdz[6];
  dnxi[5] = dXdz[7];

  basis::template addInterpFieldsGradTranspose<3, 3>(pt, dnxi, dfn);

  // U = assembleFrame(u0xi, d0): only the d0 (third column) direction is
  // differentiated here (u0xi is the vars direction, handled elsewhere)
  TacsScalar dd0[3];
  dd0[0] = dU[2];
  dd0[1] = dU[5];
  dd0[2] = dU[8];

  // V = assembleFrame(d0xi)
  TacsScalar dd0xi[6];
  dd0xi[0] = dV[0];
  dd0xi[1] = dV[1];
  dd0xi[2] = dV[3];
  dd0xi[3] = dV[4];
  dd0xi[4] = dV[6];
  dd0xi[5] = dV[7];

  basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, dd0xi, dd);
}

/**
  Add/accumulate the contributions to the residual from the coefficients
  of u0x, u1x and Ct

  @param pt The parametric point
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param du0x Coefficients for u0x
  @param du1x Coefficients for u1x
  @param res The residual
  @param dd Residual intermediate for the director field
*/
template <int vars_per_node, class basis>
void TacsShellAddDispGradSens(const double pt[], const TacsScalar T[],
                              const TacsScalar XdinvT[],
                              const TacsScalar XdinvzT[],
                              const TacsScalar du0x[], const TacsScalar du1x[],
                              TacsScalar res[], TacsScalar dd[]) {
  // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
  TacsScalar du0d[9], tmp[9];
  mat3x3MatTransMult(du1x, XdinvzT, tmp);
  mat3x3MatTransMultAdd(du0x, XdinvT, tmp);
  mat3x3MatMult(T, tmp, du0d);

  // Compute du1d = T*du1x*XdinvT^{T}
  TacsScalar du1d[9];
  mat3x3MatTransMult(du1x, XdinvT, tmp);
  mat3x3MatMult(T, tmp, du1d);

  // du0d = [du0xi; dd0]
  TacsScalar du0xi[6], dd0[3];
  TacsShellExtractFrame(du0d, du0xi, dd0);

  TacsScalar dd0xi[6];
  TacsShellExtractFrame(du1d, dd0xi);

  // Compute the director field and the gradient of the director
  // field at the specified point
  basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, dd0xi, dd);

  // Compute the gradient of the displacement solution at the quadrature points
  if (res) {
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, du0xi,
                                                                   res);
  }
}

/**
  Add/accumulate the contribution to the Jacobian matrix from the coefficients
  of u0x, u1x

  @param alpha Scale factor applied to the Jacobian entries
  @param pt The parametric point
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param d2u0x Coefficients for the second derivative of u0x
  @param d2u1x Coefficients for the second derivative of u1x
  @param d2u0xu1x Coefficients for mixed partials of u0x/u1x
  @param mat The Jacobian matrix
  @param d2d The Jacobian matrix intermediate for the director field
  @param d2du The mixed Jacobian matrix intermediate for displacements/d
*/
template <int vars_per_node, class basis>
void TacsShellAddDispGradHessian(const double pt[], const TacsScalar T[],
                                 const TacsScalar XdinvT[],
                                 const TacsScalar XdinvzT[],
                                 const TacsScalar d2u0x[],
                                 const TacsScalar d2u1x[],
                                 const TacsScalar d2u0xu1x[], TacsScalar mat[],
                                 TacsScalar d2d[], TacsScalar d2du[]) {
  // d2u0d = d2u0x*[d(u0x)/d(u0d)]^2 +
  //         d2u1x*[d(u1x)/d(u0d)]^2 +
  //         d2u0xu1x*[d(u0x)/d(u0d)*d(u1x)/d(u0x)]
  TacsScalar d2u0d[81];
  mat3x3TransMatMatHessian(T, XdinvT, d2u0x, d2u0d);
  mat3x3TransMatMatHessianAdd(T, XdinvzT, d2u1x, d2u0d);
  mat3x3TransMatMatHessianAddSymm(T, T, XdinvzT, XdinvT, d2u0xu1x, d2u0d);

  // d2u1d = d2u1x*[d(u1x)/d(u1d)]^2
  TacsScalar d2u1d[81];
  mat3x3TransMatMatHessian(T, XdinvT, d2u1x, d2u1d);

  // d2u0du1d = d2u1x*[d(u1x)/d(u0d)*d(u1x)/d(u1d)] +
  //            d2u0xu1x*[d(u0x)/d(u0d)*d(u1x)/d(u1d)]
  TacsScalar d2u0du1d[81];
  mat3x3TransMatMatHessian(T, XdinvT, d2u0xu1x, d2u0du1d);
  mat3x3TransMatMatHessianAdd(T, T, XdinvT, XdinvzT, d2u1x, d2u0du1d);

  // Extract the frame derivatives
  TacsScalar d2u0xi[36], d2d0[9], d2d0u0xi[18];
  TacsShellExtractFrameSens(d2u0d, d2u0xi, d2d0, d2d0u0xi);

  TacsScalar d2d0xi[36];
  TacsShellExtractFrameSens(d2u1d, d2d0xi);

  TacsScalar d2d0xiu0xi[36], d2d0d0xi[18];
  TacsShellExtractFrameMixedSens(d2u0du1d, d2d0xiu0xi, d2d0d0xi);

  // df/(d(u0d)d(u1d))
  basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2d0, d2d);
  basis::template addInterpGradOuterProduct<3, 3, 3, 3>(pt, d2d0xi, d2d);
  basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt, d2d0d0xi,
                                                             d2d0d0xi, d2d);

  // df/(d(d0)d(u0xi))
  basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt, d2d0u0xi, NULL,
                                                             d2du);
  basis::template addInterpGradOuterProduct<3, 3, 3, 3>(pt, d2d0xiu0xi, d2du);

  // Add the contribution to the Jacobian matrix
  if (mat) {
    basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3,
                                              3>(pt, d2u0xi, mat);
  }
}

/**
  Compute the drilling strain penalty at each node

  @param transform Transformation object
  @param Xdn The frame derivatives at each node
  @param fn The frame normals at each node
  @param vars The state variable values
  @param XdinvTn Computed inverse frame times transformation at each node
  @param Tn The transformation at each node
  @param u0xn The derivative of the displacements at each node
  @param Ctn The rotation matrix at each node
  @param etn The drill strain penalty value at each node
*/
template <int vars_per_node, int offset, class basis, class director,
          class model>
void TacsShellComputeDrillStrain(TACSShellTransform *transform,
                                 const TacsScalar Xdn[], const TacsScalar fn[],
                                 const TacsScalar vars[], TacsScalar XdinvTn[],
                                 TacsScalar Tn[], TacsScalar u0xn[],
                                 TacsScalar Ctn[], TacsScalar etn[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the transformation at the node
    TacsScalar Xxi[6];
    TacsShellExtractFrame(&Xdn[9 * i], Xxi);
    transform->computeTransform(Xxi, &fn[3 * i], &Tn[9 * i]);

    // Compute the field gradient at the node
    TacsScalar u0xi[6];
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);

    // Compute the inverse transformation
    TacsScalar Xdinv[9];
    inv3x3(&Xdn[9 * i], Xdinv);

    // Compute XdinvT = Xdinv*T
    mat3x3MatMult(Xdinv, &Tn[9 * i], &XdinvTn[9 * i]);
    TacsShellAssembleFrame(u0xi, &u0xn[9 * i]);  // Use u0x to store [u0,xi; 0]

    // Compute the rotation matrix at this node
    TacsScalar C[9], tmp[9];
    director::template computeRotationMat<vars_per_node, offset, 1>(
        &vars[vars_per_node * i], C);

    // Compute Ct = T^{T}*C*T
    mat3x3TransMatMult(&Tn[9 * i], C, tmp);
    mat3x3MatMult(tmp, &Tn[9 * i], &Ctn[9 * i]);

    // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
    // u0x = T^{T}*u0d*Xdinv*T
    mat3x3MatMult(&u0xn[9 * i], &XdinvTn[9 * i], tmp);
    mat3x3TransMatMult(&Tn[9 * i], tmp, &u0xn[9 * i]);

    etn[i] = director::evalDrillStrain(&u0xn[9 * i], &Ctn[9 * i]);
  }
}

/**
  Compute the drilling strain penalty at each node

  @param transform Transformation object
  @param Xdn The frame derivatives at each node
  @param fn The frame normals at each node
  @param vars The state variable values
  @param varsd The derivative of the state variable values
  @param XdinvTn Computed inverse frame times transformation at each node
  @param Tn The transformation at each node
  @param u0xn The derivative of the displacements at each node
  @param Ctn The rotation matrix at each node
  @param etn The drill strain penalty value at each node
*/
template <int vars_per_node, int offset, class basis, class director,
          class model>
void TacsShellComputeDrillStrainDeriv(
    TACSShellTransform *transform, const TacsScalar Xdn[],
    const TacsScalar fn[], const TacsScalar vars[], const TacsScalar varsd[],
    TacsScalar XdinvTn[], TacsScalar Tn[], TacsScalar u0xn[], TacsScalar Ctn[],
    TacsScalar etn[], TacsScalar etnd[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the transformation at the node
    TacsScalar Xxi[6];
    TacsShellExtractFrame(&Xdn[9 * i], Xxi);
    transform->computeTransform(Xxi, &fn[3 * i], &Tn[9 * i]);

    // Compute the field gradient at the node
    TacsScalar u0xi[6], u0xid[6];
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, u0xid);

    // Compute the inverse transformation
    TacsScalar Xdinv[9];
    inv3x3(&Xdn[9 * i], Xdinv);

    // Compute XdinvT = Xdinv*T
    TacsScalar u0xnd[9];
    mat3x3MatMult(Xdinv, &Tn[9 * i], &XdinvTn[9 * i]);
    TacsShellAssembleFrame(u0xi, &u0xn[9 * i]);  // Use u0x to store [u0,xi; 0]
    TacsShellAssembleFrame(u0xid, u0xnd);

    // Compute the rotation matrix at this node
    TacsScalar C[9], Cd[9], tmp[9];
    director::template computeRotationMatDeriv<vars_per_node, offset, 1>(
        &vars[vars_per_node * i], &varsd[vars_per_node * i], C, Cd);

    // Compute Ct = T^{T}*C*T
    mat3x3TransMatMult(&Tn[9 * i], C, tmp);
    mat3x3MatMult(tmp, &Tn[9 * i], &Ctn[9 * i]);

    TacsScalar Ctnd[9];
    mat3x3TransMatMult(&Tn[9 * i], Cd, tmp);
    mat3x3MatMult(tmp, &Tn[9 * i], Ctnd);

    // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
    // u0x = T^{T}*u0d*Xdinv*T
    mat3x3MatMult(&u0xn[9 * i], &XdinvTn[9 * i], tmp);
    mat3x3TransMatMult(&Tn[9 * i], tmp, &u0xn[9 * i]);

    mat3x3MatMult(u0xnd, &XdinvTn[9 * i], tmp);
    mat3x3TransMatMult(&Tn[9 * i], tmp, u0xnd);

    etn[i] = director::evalDrillStrainDeriv(&u0xn[9 * i], &Ctn[9 * i], u0xnd,
                                            Ctnd, &etnd[i]);
  }
}

/**
  Add the derivative of the drilling strain penalty to the residual

  @param transform Transformation object
  @param Xdn The frame derivatives at each node
  @param fn The frame normals at each node
  @param vars The state variable values
  @param XdinvTn Computed inverse frame times transformation at each node
  @param Tn The transformation at each node
  @param u0xn The derivative of the displacements at each node
  @param Ctn The rotation matrix at each node
  @param detn The derivative of the drill strain penalty value at each node
  @param res The element residual
*/
template <int vars_per_node, int offset, class basis, class director,
          class model>
void TacsShellAddDrillStrainSens(const TacsScalar Xdn[], const TacsScalar fn[],
                                 const TacsScalar vars[],
                                 const TacsScalar XdinvTn[],
                                 const TacsScalar Tn[], const TacsScalar u0xn[],
                                 const TacsScalar Ctn[],
                                 const TacsScalar detn[], TacsScalar res[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    TacsScalar du0x[9], dCt[9];
    director::evalDrillStrainSens(detn[i], &u0xn[9 * i], &Ctn[9 * i], du0x,
                                  dCt);

    // Compute dCpt = T*dCt*T^{T}
    TacsScalar dCpt[9], tmp[9];
    mat3x3MatMult(&Tn[9 * i], dCt, tmp);
    mat3x3MatTransMult(tmp, &Tn[9 * i], dCpt);

    director::template addRotationMatResidual<vars_per_node, offset, 1>(
        &vars[i * vars_per_node], dCpt, &res[i * vars_per_node]);

    // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
    TacsScalar du0d[9];
    mat3x3MatTransMult(du0x, &XdinvTn[9 * i], tmp);
    mat3x3MatMult(&Tn[9 * i], tmp, du0d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6];
    TacsShellExtractFrame(du0d, du0xi);

    // Compute the gradient of the displacement solution at the quadrature
    // points
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, du0xi,
                                                                   res);
  }
}

/**
  Add the Xpts/fn sensitivity of the drilling strain penalty

  Given a per-node adjoint seed detn on the drill strain penalty etn computed
  by TacsShellComputeDrillStrain(Deriv), accumulate the corresponding
  sensitivity into dfdXpts (through each node's own Xxi) and dfn (through
  each node's own fn, closed later by TacsShellAddNodeNormalsSens). No
  vars-direction work is added here (that is TacsShellAddDrillStrainSens's
  job, unchanged).

  @param transform Transformation object
  @param Xdn The frame derivatives at each node
  @param fn The frame normals at each node
  @param vars The state variable values
  @param XdinvTn Computed inverse frame times transformation at each node
  @param Tn The transformation at each node
  @param u0xn The derivative of the displacements at each node
  @param Ctn The rotation matrix at each node
  @param detn The adjoint seed on the drill strain penalty at each node
  @param dfdXpts The accumulated sensitivity with respect to Xpts
  @param dfn The accumulated sensitivity with respect to the nodal fn field
*/
template <int vars_per_node, int offset, class basis, class director,
          class model>
void TacsShellAddDrillStrainXptSens(
    TACSShellTransform *transform, const TacsScalar Xdn[],
    const TacsScalar fn[], const TacsScalar vars[],
    const TacsScalar XdinvTn[], const TacsScalar Tn[], const TacsScalar u0xn[],
    const TacsScalar Ctn[], const TacsScalar detn[], TacsScalar dfdXpts[],
    TacsScalar dfn[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Recompute the forward quantities needed for the reverse sweep (mirrors
    // TacsShellComputeDrillStrain's per-node body)
    TacsScalar Xxi[6];
    TacsShellExtractFrame(&Xdn[9 * i], Xxi);

    TacsScalar Xdinv[9];
    inv3x3(&Xdn[9 * i], Xdinv);

    TacsScalar u0xi[6];
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);
    TacsScalar u0xnLocal[9];
    TacsShellAssembleFrame(u0xi, u0xnLocal);  // [u0,xi; 0]

    TacsScalar C[9];
    director::template computeRotationMat<vars_per_node, offset, 1>(
        &vars[vars_per_node * i], C);

    const TacsScalar *T = &Tn[9 * i];
    const TacsScalar *XdinvT = &XdinvTn[9 * i];

    // ---- Reverse sweep ----
    TacsScalar du0x[9], dCt[9];
    director::evalDrillStrainSens(detn[i], &u0xn[9 * i], &Ctn[9 * i], du0x,
                                  dCt);

    TacsScalar dT[9];
    memset(dT, 0, sizeof(dT));

    // u0x = T^{T}*tmp, tmp = u0xnLocal*XdinvT
    TacsScalar tmp[9], dtmp[9];
    mat3x3MatMult(u0xnLocal, XdinvT, tmp);
    mat3x3MatMult(T, du0x, dtmp);
    TacsScalar dTloc[9];
    mat3x3MatTransMult(tmp, du0x, dTloc);
    for (int k = 0; k < 9; k++) dT[k] += dTloc[k];

    // tmp = u0xnLocal*XdinvT (only need dXdinvT; u0xnLocal is vars-direction,
    // not differentiated here)
    TacsScalar dXdinvT[9];
    mat3x3TransMatMult(u0xnLocal, dtmp, dXdinvT);

    // Ct = T^{T}*M, M = C*T
    TacsScalar M[9];
    mat3x3MatMult(C, T, M);
    TacsScalar dM[9];
    mat3x3MatMult(T, dCt, dM);
    mat3x3MatTransMult(M, dCt, dTloc);
    for (int k = 0; k < 9; k++) dT[k] += dTloc[k];

    // M = C*T
    mat3x3TransMatMult(C, dM, dTloc);
    for (int k = 0; k < 9; k++) dT[k] += dTloc[k];

    // XdinvT = Xdinv*T
    TacsScalar dXdinv[9];
    mat3x3MatTransMult(dXdinvT, T, dXdinv);
    mat3x3TransMatMult(Xdinv, dXdinvT, dTloc);
    for (int k = 0; k < 9; k++) dT[k] += dTloc[k];

    // Xdinv = inv3x3(Xdn_i): dXdn += -Xdinv^{T}*dXdinv*Xdinv^{T}
    TacsScalar dXdn[9];
    inv3x3Sens(Xdinv, dXdinv, dXdn);

    // T = transform->computeTransform(Xxi, fn_i)
    TacsScalar dXxi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    TacsScalar dfnLocal[3] = {0.0, 0.0, 0.0};
    transform->addTransformSens(Xxi, &fn[3 * i], dT, dXxi, dfnLocal);

    // Xdn_i = assembleFrame(Xxi, fn_i): fold in the Xdinv-adjoint contribution
    dXxi[0] += dXdn[0];
    dXxi[1] += dXdn[1];
    dXxi[2] += dXdn[3];
    dXxi[3] += dXdn[4];
    dXxi[4] += dXdn[6];
    dXxi[5] += dXdn[7];

    dfnLocal[0] += dXdn[2];
    dfnLocal[1] += dXdn[5];
    dfnLocal[2] += dXdn[8];

    // Xxi_i = basis::interpFieldsGrad(pt_i, Xpts) directly at node i's own
    // parametric point
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);

    dfn[3 * i] += dfnLocal[0];
    dfn[3 * i + 1] += dfnLocal[1];
    dfn[3 * i + 2] += dfnLocal[2];
  }
}

/**
  Add the first and second derivatives of the drilling strain penalty
  to the residual and Jacobian matrix

  @param Xdn The frame derivatives at each node
  @param fn The frame normals at each node
  @param vars The state variable values
  @param XdinvTn Computed inverse frame times transformation at each node
  @param Tn The transformation at each node
  @param u0xn The derivative of the displacements at each node
  @param Ctn The rotation matrix at each node
  @param detn The derivative of the drill strain penalty value at each node
  @param d2etn The second derivative of the drill strain penalty value at each
  node
  @param res The element residual
  @param mat The element Jacobian matrix
*/
template <int vars_per_node, int offset, class basis, class director,
          class model>
void TacsShellAddDrillStrainHessian(
    const TacsScalar Xdn[], const TacsScalar fn[], const TacsScalar vars[],
    const TacsScalar XdinvTn[], const TacsScalar Tn[], const TacsScalar u0xn[],
    const TacsScalar Ctn[], const TacsScalar detn[], const TacsScalar d2etn[],
    TacsScalar res[], TacsScalar mat[]) {
  // Zero the derivative of the rotation
  const int size = vars_per_node * basis::NUM_NODES;
  const int num_nodes = basis::NUM_NODES;

  TacsScalar drot[num_nodes * size];
  memset(drot, 0, num_nodes * size * sizeof(TacsScalar));

  for (int i = 0; i < num_nodes; i++) {
    // Set the
    TacsScalar *t = &drot[i * size];

    // Get the node location
    double pt[2];
    basis::getNodePoint(i, pt);

    TacsScalar du0x[9], dCt[9];
    director::evalDrillStrainSens(1.0, &u0xn[9 * i], &Ctn[9 * i], du0x, dCt);

    // Compute dCpt = T*dCt*T^{T}
    TacsScalar dCpt[9], tmp[9];
    mat3x3MatMult(&Tn[9 * i], dCt, tmp);
    mat3x3MatTransMult(tmp, &Tn[9 * i], dCpt);

    director::template addRotationMatResidual<vars_per_node, offset, 1>(
        &vars[i * vars_per_node], dCpt, &t[i * vars_per_node]);

    // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
    TacsScalar du0d[9];
    mat3x3MatTransMult(du0x, &XdinvTn[9 * i], tmp);
    mat3x3MatMult(&Tn[9 * i], tmp, du0d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6];
    TacsShellExtractFrame(du0d, du0xi);

    // Compute the gradient of the displacement solution at the quadrature
    // points
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, du0xi,
                                                                   t);
  }

  // Add the contribution to the residual from the drilling rotation
  if (res) {
    for (int i = 0; i < num_nodes; i++) {
      const TacsScalar *t = &drot[i * size];
      for (int j = 0; j < size; j++) {
        res[j] += detn[i] * t[j];
      }
    }
  }

  // Add the contributions to the Jacobian
  for (int i = 0; i < num_nodes; i++) {
    TacsScalar t[size];
    memset(
        t, 0,
        size * sizeof(TacsScalar));  // *const TacsScalar *t1 = &drot[i*size];

    for (int j = 0; j < num_nodes; j++) {
      const TacsScalar *t1 = &drot[j * size];
      for (int ii = 0; ii < size; ii++) {
        t[ii] += d2etn[i * num_nodes + j] * t1[ii];
      }
    }

    const TacsScalar *t2 = &drot[i * size];
    for (int ii = 0; ii < size; ii++) {
      for (int jj = 0; jj < size; jj++) {
        mat[ii * size + jj] += t[ii] * t2[jj];
      }
    }
  }

  for (int i = 0; i < num_nodes; i++) {
    // Get the node location
    double pt[2];
    basis::getNodePoint(i, pt);

    TacsScalar d2u0x[81], d2Ct[81], d2Ctu0x[81];
    director::evalDrillStrainHessian(detn[i], &u0xn[9 * i], &Ctn[9 * i], d2u0x,
                                     d2Ct, d2Ctu0x);

    // Compute the second derivative w.r.t. u0d
    TacsScalar d2u0d[81];
    mat3x3TransMatMatHessian(&Tn[9 * i], &XdinvTn[9 * i], d2u0x, d2u0d);

    // // Compute the second derivative w.r.t. Cpt
    TacsScalar d2Cpt[81];
    mat3x3TransMatMatHessian(&Tn[9 * i], &Tn[9 * i], d2Ct, d2Cpt);

    // d2C0u0d = Ctu0x*[d(u0x)/d(u0d)]*[d(Ct)/d(Cpt)]
    TacsScalar d2Cptu0d[81];
    mat3x3TransMatMatHessianAdd(&Tn[9 * i], &Tn[9 * i], &XdinvTn[9 * i],
                                &Tn[9 * i], d2Ctu0x, d2Cptu0d);

    // Extract the shell frame
    TacsScalar d2u0xi[36];
    TacsShellExtractFrameSens(d2u0d, d2u0xi);

    basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3,
                                              3>(pt, d2u0xi, mat);
  }
}

/**
  Add/accumulate the contribution to the Jacobians from the coupling between the
  tying strain and the displacement gradient

  @param pt The parametric point
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param d2e0tyu0x The second derivative of the tying strain and u0x
  @param d2e0tyu1x The second derivative of the tying strain and u1x
  @param d2etyu The second derivative of the tying strain and displacements
  @param d2etyd The second derivative of the tying strain and the director
*/
template <class basis>
void TacsShellAddTyingDispCoupling(const double pt[], const TacsScalar T[],
                                   const TacsScalar XdinvT[],
                                   const TacsScalar XdinvzT[],
                                   const TacsScalar d2e0tyu0x[],
                                   const TacsScalar d2e0tyu1x[],
                                   TacsScalar d2etyu[], TacsScalar d2etyd[]) {
  const int usize = 3 * basis::NUM_NODES;
  const int dsize = 3 * basis::NUM_NODES;

  // Compute d2gtyu0d, d2gtyu1d
  TacsScalar d2gtyu0d[54], d2gtyu1d[54];
  TacsScalar tmp0d[54], tmp1d[54];
  for (int k = 0; k < 6; k++) {
    // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
    TacsScalar tmp[9];
    mat3x3MatTransMult(&d2e0tyu1x[9 * k], XdinvzT, tmp);
    mat3x3MatTransMultAdd(&d2e0tyu0x[9 * k], XdinvT, tmp);
    mat3x3MatMult(T, tmp, &tmp0d[9 * k]);

    // Compute du1d = T*du1x*XdinvT^{T}
    mat3x3MatTransMult(&d2e0tyu1x[9 * k], XdinvT, tmp);
    mat3x3MatMult(T, tmp, &tmp1d[9 * k]);
  }

  // Perform the sensitivity transformation
  for (int k = 0; k < 9; k++) {
    TacsScalar t0[6], out[6];
    for (int kk = 0; kk < 6; kk++) {
      t0[kk] = tmp0d[9 * kk + k];
    }
    mat3x3SymmTransformTransSens(XdinvT, t0, out);
    for (int kk = 0; kk < 6; kk++) {
      d2gtyu0d[9 * kk + k] = out[kk];
    }
  }

  for (int k = 0; k < 9; k++) {
    TacsScalar t0[6], out[6];
    for (int kk = 0; kk < 6; kk++) {
      t0[kk] = tmp1d[9 * kk + k];
    }
    mat3x3SymmTransformTransSens(XdinvT, t0, out);
    for (int kk = 0; kk < 6; kk++) {
      d2gtyu1d[9 * kk + k] = out[kk];
    }
  }

  TacsScalar d2gtyu[6 * usize], d2gtyd[6 * dsize];
  memset(d2gtyu, 0, 6 * usize * sizeof(TacsScalar));
  memset(d2gtyd, 0, 6 * dsize * sizeof(TacsScalar));

  for (int k = 0; k < 6; k++) {
    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6], dd0[3];
    TacsShellExtractFrame(&d2gtyu0d[9 * k], du0xi, dd0);

    TacsScalar dd0xi[6];
    TacsShellExtractFrame(&d2gtyu1d[9 * k], dd0xi);

    // Compute the director field and the gradient of the director
    // field at the specified point
    basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, &d2gtyd[dsize * k]);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, dd0xi,
                                                       &d2gtyd[dsize * k]);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, du0xi,
                                                       &d2gtyu[usize * k]);
  }

  // Add the values into d2etyu and d2etyd
  for (int k = 0; k < usize; k++) {
    TacsScalar t1[6], t2[basis::NUM_TYING_POINTS];
    memset(t2, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    for (int kk = 0; kk < 6; kk++) {
      t1[kk] = d2gtyu[usize * kk + k];
    }

    basis::addInterpTyingStrainTranspose(pt, t1, t2);

    for (int kk = 0; kk < basis::NUM_TYING_POINTS; kk++) {
      d2etyu[kk * usize + k] += t2[kk];
    }
  }

  for (int k = 0; k < dsize; k++) {
    TacsScalar t1[6], t2[basis::NUM_TYING_POINTS];
    memset(t2, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    for (int kk = 0; kk < 6; kk++) {
      t1[kk] = d2gtyd[dsize * kk + k];
    }

    basis::addInterpTyingStrainTranspose(pt, t1, t2);

    for (int kk = 0; kk < basis::NUM_TYING_POINTS; kk++) {
      d2etyd[kk * dsize + k] += t2[kk];
    }
  }
}

/*
  Test the implementation of the shell terms for a given basis
*/
template <int vars_per_node, class basis>
int TacsTestShellUtilities(double dh = 1e-7, int test_print_level = 2,
                           double test_fail_atol = 1e-5,
                           double test_fail_rtol = 1e-5) {
  const int size = vars_per_node * basis::NUM_NODES;
  const int usize = 3 * basis::NUM_NODES;
  const int dsize = 3 * basis::NUM_NODES;
  const int xsize = 3 * basis::NUM_NODES;

  double pt[2];
  TacsGenerateRandomArray(pt, 2);

  TacsScalar T[9];
  TacsGenerateRandomArray(T, 9);

  TacsScalar Xpts[xsize], fn[xsize];
  TacsGenerateRandomArray(Xpts, xsize);
  TacsGenerateRandomArray(fn, xsize);

  TacsScalar n0[3], Xxi[6];
  basis::template interpFields<3, 3>(pt, fn, n0);
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

  TacsScalar vars[size], d[dsize];
  TacsGenerateRandomArray(vars, size);
  TacsGenerateRandomArray(d, dsize);

  // Generate random perturbations for the linear terms
  TacsScalar du0x[9], du1x[9];
  TacsGenerateRandomArray(du0x, 9);
  TacsGenerateRandomArray(du1x, 9);

  // Generate the derivative terms
  TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
  TacsGenerateRandomArray(d2u0x, 81);
  TacsGenerateRandomArray(d2u1x, 81);
  TacsGenerateRandomArray(d2u0xu1x, 81);

  // Symmetrize the random matrices
  for (int i = 0; i < 9; i++) {
    for (int j = i + 1; j < 9; j++) {
      d2u0x[9 * i + j] = d2u0x[9 * j + i];
      d2u1x[9 * i + j] = d2u1x[9 * j + i];
    }
  }

  TacsScalar res[size], mat[size * size];
  memset(res, 0, size * sizeof(TacsScalar));
  memset(mat, 0, size * size * sizeof(TacsScalar));

  TacsScalar dd[dsize], d2d[dsize * dsize], d2du[usize * dsize];
  memset(dd, 0, dsize * sizeof(TacsScalar));
  memset(d2d, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2du, 0, usize * dsize * sizeof(TacsScalar));

  TacsScalar XdinvT[9], XdinvzT[9];
  TacsScalar u0x[9], u1x[9];
  TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, d, Xxi, n0,
                                                 T, XdinvT, XdinvzT, u0x, u1x);

  TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT, du0x,
                                                 du1x, res, dd);
  TacsShellAddDispGradHessian<vars_per_node, basis>(
      pt, T, XdinvT, XdinvzT, d2u0x, d2u1x, d2u0xu1x, mat, d2d, d2du);

  // Now, check the result
  TacsScalar fdmat[size * size], fddu[dsize * usize];
  for (int k = 0; k < size; k++) {
    TacsScalar varst[size];
    memcpy(varst, vars, size * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh;
#endif  // TACS_USE_COMPLEX

    // Compute the pertubation
    TacsScalar u0xt[9], u1xt[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, varst, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0xt, u1xt);

    // d2u0xu1x[9*i + j] = p2f/(p(u0x[i]) p(u1x[j]))
    // d2Ctu0x[9*i + j] = p2f/(p(Ct[i]) p(u0x[j]))

    // Compute the perturbed values based on the outputs
    TacsScalar du0xt[9], du1xt[9];
    for (int i = 0; i < 9; i++) {
      du0xt[i] = du0x[i];
      du1xt[i] = du1x[i];
      for (int j = 0; j < 9; j++) {
        du0xt[i] += d2u0x[9 * i + j] * (u0xt[j] - u0x[j]) +
                    d2u0xu1x[9 * i + j] * (u1xt[j] - u1x[j]);

        du1xt[i] += d2u1x[9 * i + j] * (u1xt[j] - u1x[j]) +
                    d2u0xu1x[9 * j + i] * (u0xt[j] - u0x[j]);
      }
    }

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size * sizeof(TacsScalar));
    memset(ddt, 0, dsize * sizeof(TacsScalar));
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                   du0xt, du1xt, rest, ddt);

    // Place the result in the arrays...
    for (int j = 0; j < size; j++) {
#ifdef TACS_USE_COMPLEX
      fdmat[k + size * j] = TacsImagPart(rest[j]) / dh;
#else
      fdmat[k + size * j] = (rest[j] - res[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }

    if (k % vars_per_node < 3) {
      // Compute the u-index
      int index = 3 * (k / vars_per_node) + k % vars_per_node;

      for (int j = 0; j < dsize; j++) {
#ifdef TACS_USE_COMPLEX
        fddu[index + usize * j] = TacsImagPart(ddt[j]) / dh;
#else
        fddu[index + usize * j] = (ddt[j] - dd[j]) / dh;
#endif  // TACS_USE_COMPLEX
      }
    }
  }

  // Variables to store the max error and indices
  int max_err_index, max_rel_index;
  double max_err, max_rel;

  // Keep track of the failure flag
  int fail = 0;

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size * size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size * size, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size * size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the error
  max_err = TacsGetMaxError(d2du, fddu, usize * dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2du, fddu, usize * dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. d and vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2du", d2du, fddu, usize * dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  TacsScalar fddd[dsize * dsize];
  for (int k = 0; k < dsize; k++) {
    TacsScalar dt[dsize];
    memcpy(dt, d, dsize * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    dt[k] = d[k] + TacsScalar(0.0, dh);
#else
    dt[k] = d[k] + dh;
#endif  // TACS_USE_COMPLEX

    // Compute the pertubation
    TacsScalar u0xt[9], u1xt[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, dt, Xxi, n0, T, XdinvT, XdinvzT, u0xt, u1xt);

    // d2u0xu1x[9*i + j] = p2f/(p(u0x[i]) p(u1x[j]))
    // d2Ctu0x[9*i + j] = p2f/(p(Ct[i]) p(u0x[j]))

    // Compute the perturbed values based on the outputs
    TacsScalar du0xt[9], du1xt[9];
    for (int i = 0; i < 9; i++) {
      du0xt[i] = du0x[i];
      du1xt[i] = du1x[i];
      for (int j = 0; j < 9; j++) {
        du0xt[i] += d2u0x[9 * i + j] * (u0xt[j] - u0x[j]) +
                    d2u0xu1x[9 * i + j] * (u1xt[j] - u1x[j]);

        du1xt[i] += d2u1x[9 * i + j] * (u1xt[j] - u1x[j]) +
                    d2u0xu1x[9 * j + i] * (u0xt[j] - u0x[j]);
      }
    }

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size * sizeof(TacsScalar));
    memset(ddt, 0, dsize * sizeof(TacsScalar));
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                   du0xt, du1xt, rest, ddt);

    for (int j = 0; j < dsize; j++) {
#ifdef TACS_USE_COMPLEX
      fddd[k + dsize * j] = TacsImagPart(ddt[j]) / dh;
#else
      fddd[k + dsize * j] = (ddt[j] - dd[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2d, fddd, dsize * dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d, fddd, dsize * dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2d", d2d, fddd, dsize * dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif  // TACS_SHELL_UTILITIES_H
