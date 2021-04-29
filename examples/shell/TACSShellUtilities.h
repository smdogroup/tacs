#ifndef TACS_SHELL_UTILITIES_H
#define TACS_SHELL_UTILITIES_H

#include "TACSElementAlgebra.h"

inline void assembleFrame( const TacsScalar Xxi[],
                           const TacsScalar n[],
                           TacsScalar Xd[] ){
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

inline void assembleFrame( const TacsScalar nxi[],
                           TacsScalar Xdz[] ){
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

inline void extractFrame( const TacsScalar Xd[],
                          TacsScalar Xxi[],
                          TacsScalar n[] ){
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

inline void extractFrame( const TacsScalar Xd[],
                          TacsScalar Xxi[] ){
  Xxi[0] = Xd[0];
  Xxi[1] = Xd[1];

  Xxi[2] = Xd[3];
  Xxi[3] = Xd[4];

  Xxi[4] = Xd[6];
  Xxi[5] = Xd[7];
}

inline void extractFrameSens( const TacsScalar d2u0d[],
                              TacsScalar d2u0xi[] ){
  // Extract the second derivatives
  for ( int j = 0; j < 6; j++ ){
    int jj = 3*(j / 2) + (j % 2);
    for ( int i = 0; i < 6; i++ ){
      int ii = 3*(i / 2) + (i % 2);

      d2u0xi[6*j + i] = d2u0d[9*jj + ii];
    }
  }
}

inline void extractFrameSens( const TacsScalar d2u0d[],
                              TacsScalar d2u0xi[],
                              TacsScalar d2d0[],
                              TacsScalar d2d0u0xi[] ){
  // Extract the second derivatives
  for ( int j = 0; j < 6; j++ ){
    int jj = 3*(j / 2) + (j % 2);
    for ( int i = 0; i < 6; i++ ){
      int ii = 3*(i / 2) + (i % 2);

      d2u0xi[6*j + i] = d2u0d[9*jj + ii];
    }
  }

  for ( int j = 0; j < 3; j++ ){
    int jj = 3*j + 2;
    for ( int i = 0; i < 3; i++ ){
      int ii = 3*i + 2;
      d2d0[3*j + i] = d2u0d[9*jj + ii];
    }
  }

  for ( int j = 0; j < 3; j++ ){
    int jj = 3*j + 2;
    for ( int i = 0; i < 6; i++ ){
      int ii = 3*(i / 2) + (i % 2);
      d2d0u0xi[6*j + i] = d2u0d[9*jj + ii];
    }
  }
}

inline void extractFrameMixedSens( const TacsScalar d2u0du1d[],
                                   TacsScalar d2d0xiu0xi[],
                                   TacsScalar d2d0d0xi[] ){
  // Extract the second derivatives
  for ( int j = 0; j < 6; j++ ){
    int jj = 3*(j / 2) + (j % 2);
    for ( int i = 0; i < 6; i++ ){
      int ii = 3*(i / 2) + (i % 2);

      d2d0xiu0xi[6*i + j] = d2u0du1d[9*jj + ii];
    }
  }

  for ( int j = 0; j < 3; j++ ){
    int jj = 3*j + 2;
    for ( int i = 0; i < 6; i++ ){
      int ii = 3*(i / 2) + (i % 2);
      d2d0d0xi[6*j + i] = d2u0du1d[9*jj + ii];
    }
  }
}

/*
  Compute the second derivative of the transformation

  C = T^{T}*A*B
  dA = T*dC*B^{T}

  Given d
*/
inline void mat3x3TransMatMatHessian( const TacsScalar T[],
                                      const TacsScalar B[],
                                      const TacsScalar d2C[],
                                      TacsScalar d2A[] ){
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for ( int i = 0; i < 9; i++ ){
    mat3x3MatTransMult(&d2C[9*i], B, t1);
    mat3x3MatMult(T, t1, &tmp2[9*i]);
  }

  for ( int i = 0; i < 9; i++ ){
    for ( int j = 0; j < 9; j++ ){
      t1[j] = tmp2[i + 9*j];
    }

    mat3x3MatTransMult(t1, B, t2);
    mat3x3MatMult(T, t2, t1);

    for ( int j = 0; j < 9; j++ ){
      d2A[i + 9*j] = t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAdd( const TacsScalar T[],
                                         const TacsScalar B[],
                                         const TacsScalar d2C[],
                                         TacsScalar d2A[] ){
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for ( int i = 0; i < 9; i++ ){
    mat3x3MatTransMult(&d2C[9*i], B, t1);
    mat3x3MatMult(T, t1, &tmp2[9*i]);
  }

  for ( int i = 0; i < 9; i++ ){
    for ( int j = 0; j < 9; j++ ){
      t1[j] = tmp2[i + 9*j];
    }

    mat3x3MatTransMult(t1, B, t2);
    mat3x3MatMult(T, t2, t1);

    for ( int j = 0; j < 9; j++ ){
      d2A[i + 9*j] += t1[j];
    }
  }
}


inline void mat3x3TransMatMatHessian( const TacsScalar T1[],
                                      const TacsScalar T2[],
                                      const TacsScalar B1[],
                                      const TacsScalar B2[],
                                      const TacsScalar d2C[],
                                      TacsScalar d2A[] ){
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for ( int i = 0; i < 9; i++ ){
    mat3x3MatTransMult(&d2C[9*i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9*i]);
  }

  for ( int i = 0; i < 9; i++ ){
    for ( int j = 0; j < 9; j++ ){
      t1[j] = tmp2[i + 9*j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for ( int j = 0; j < 9; j++ ){
      d2A[i + 9*j] = t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAdd( const TacsScalar T1[],
                                         const TacsScalar T2[],
                                         const TacsScalar B1[],
                                         const TacsScalar B2[],
                                         const TacsScalar d2C[],
                                         TacsScalar d2A[] ){
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for ( int i = 0; i < 9; i++ ){
    mat3x3MatTransMult(&d2C[9*i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9*i]);
  }

  for ( int i = 0; i < 9; i++ ){
    for ( int j = 0; j < 9; j++ ){
      t1[j] = tmp2[i + 9*j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for ( int j = 0; j < 9; j++ ){
      d2A[i + 9*j] += t1[j];
    }
  }
}

inline void mat3x3TransMatMatHessianAddSymm( const TacsScalar T1[],
                                             const TacsScalar T2[],
                                             const TacsScalar B1[],
                                             const TacsScalar B2[],
                                             const TacsScalar d2C[],
                                             TacsScalar d2A[] ){
  // Compute the second derivatives
  TacsScalar tmp2[81], t1[9], t2[9];
  for ( int i = 0; i < 9; i++ ){
    mat3x3MatTransMult(&d2C[9*i], B1, t1);
    mat3x3MatMult(T1, t1, &tmp2[9*i]);
  }

  for ( int i = 0; i < 9; i++ ){
    for ( int j = 0; j < 9; j++ ){
      t1[j] = tmp2[i + 9*j];
    }

    mat3x3MatTransMult(t1, B2, t2);
    mat3x3MatMult(T2, t2, t1);

    for ( int j = 0; j < 9; j++ ){
      d2A[i + 9*j] += t1[j];
      d2A[j + 9*i] += t1[j];
    }
  }
}

/**
  Compute the frame normals at each of the nodes of the shell element

  This code computes the in-plane coordinates derivatives to evaluate
  the shell normal.

  @param Xpts The node locations for the shell element
  @param fn The frame normals computed at each node
  @param fnorm Optional: the norm of the cross-product
*/
template <class basis>
void getNodeNormals( const TacsScalar Xpts[],
                     TacsScalar fn[],
                     TacsScalar fnorm[]=NULL ){
  for ( int i = 0; i < basis::NUM_NODES; i++ ){
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
    crossProduct(a, b, &fn[3*i]);

    // Compute the 2-norm of the vector in the normal direction
    TacsScalar norm = sqrt(vec3Dot(&fn[3*i], &fn[3*i]));

    // Save the 2-norm value if the fnorm argument is not NULL
    if (fnorm){
      fnorm[i] = norm;
    }

    // Scale the normal direction
    if (norm != 0.0){
      vec3Scale(1.0/norm, &fn[3*i]);
    }
  }
}

/**
  Compute the displacement gradient of the constant and through-thickness
  rate of change of the displacements.

  @param pt The parametric point
  @param Xpts The node locations for the element
  @param vars The element variables
  @param fn The frame normal directions at each node
  @param C The rotation matrix stored at each node
  @param d The director field at each node
  @param Xxi The in-plane coordinate derivatives
  @param n0 The interpolated frame normal direction
  @param T The transformation to local coordinates
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param u0x Derivative of the displacement in the local x coordinates
  @param u1x Derivative of the through-thickness disp. in local x coordinates
  @param Ct The interpolated rotation matrix pre/post multiplied by T
*/
template <int vars_per_node, class basis>
TacsScalar computeDispGrad( const double pt[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar fn[],
                            const TacsScalar C[],
                            const TacsScalar d[],
                            const TacsScalar Xxi[],
                            const TacsScalar n0[],
                            const TacsScalar T[],
                            TacsScalar XdinvT[],
                            TacsScalar XdinvzT[],
                            TacsScalar u0x[],
                            TacsScalar u1x[],
                            TacsScalar Ct[] ){
  // Compute n,xi = [dn/dxi1; dn/dxi2]
  TacsScalar nxi[6];
  basis::template interpFieldsGrad<3, 3>(pt, fn, nxi);

  // Assemble the terms Xd = [Xxi; n] and Xdz
  TacsScalar Xd[9], Xdz[9];
  assembleFrame(Xxi, n0, Xd);
  assembleFrame(nxi, Xdz);

  // Compute the inverse of the 3x3 Jacobian transformation
  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  // Compute negXdinvXdz = -Xdinv*Xdz
  TacsScalar negXdinvXdz[9];
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for ( int i = 0; i < 9; i++ ){
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
  assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

  // u1x = T^{T}*u1d*XdinvT + T^{T}*u0d*XdinvzT
  TacsScalar tmp[9];
  assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1x, XdinvT, tmp);
  mat3x3MatMultAdd(u0x, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1x);

  // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
  // u0x = T^{T}*u0d*Xdinv*T
  mat3x3MatMult(u0x, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0x);

  // Compute the interpolation of the entries of the C matrix
  TacsScalar Cpt[9];
  basis::template interpFields<9, 9>(pt, C, Cpt);

  // Compute Ct = T^{T}*Cpt*T
  mat3x3TransMatMult(T, Cpt, tmp);
  mat3x3MatMult(tmp, T, Ct);

  return detXd;
}

/**
  Compute the displacement gradient and the derivative of the displacement
  gradient with respect to input variables

  @param pt The parametric point
  @param Xpts The node locations for the element
  @param vars The element variables
  @param fn The frame normal directions at each node
  @param C The rotation matrix stored at each node
  @param d The director field at each node
  @param Xxi The in-plane coordinate derivatives
  @param n0 The interpolated frame normal direction
  @param T The transformation to local coordinates
  @param varsd The derivative of the vars
  @param dd The derivative of the director field
  @param Cd The derivative of the rotation matrix
  @param XdinvT Product of inverse of the Jacobian trans. and T
  @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
  @param u0x Derivative of the displacement in the local x coordinates
  @param u1x Derivative of the through-thickness disp. in local x coordinates
  @param Ct The interpolated rotation matrix pre/post multiplied by T
  @param u0xd Derivative of u0x
  @param u1xd Derivative of u1x
  @param Ctd Derivative of Ct
*/
template <int vars_per_node, class basis>
TacsScalar computeDispGradDeriv( const double pt[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar fn[],
                                 const TacsScalar C[],
                                 const TacsScalar d[],
                                 const TacsScalar Xxi[],
                                 const TacsScalar n0[],
                                 const TacsScalar T[],
                                 const TacsScalar varsd[],
                                 const TacsScalar dd[],
                                 const TacsScalar Cd[],
                                 TacsScalar XdinvT[],
                                 TacsScalar XdinvzT[],
                                 TacsScalar u0x[],
                                 TacsScalar u1x[],
                                 TacsScalar Ct[],
                                 TacsScalar u0xd[],
                                 TacsScalar u1xd[],
                                 TacsScalar Ctd[] ){
  // Compute n,xi = [dn/dxi1; dn/dxi2]
  TacsScalar nxi[6];
  basis::template interpFieldsGrad<3, 3>(pt, fn, nxi);

  // Assemble the terms Xd = [Xxi; n] and Xdz
  TacsScalar Xd[9], Xdz[9];
  assembleFrame(Xxi, n0, Xd);
  assembleFrame(nxi, Xdz);

  // Compute the inverse of the 3x3 Jacobian transformation
  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  // Compute negXdinvXdz = -Xdinv*Xdz
  TacsScalar negXdinvXdz[9];
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for ( int i = 0; i < 9; i++ ){
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
  assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]
  assembleFrame(u0xid, d0d, u0xd);

  // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
  TacsScalar tmp[9];
  assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1x, XdinvT, tmp);
  mat3x3MatMultAdd(u0x, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1x);

  assembleFrame(d0xid, u1xd); // Use u1x to store [d0,xi; 0]
  mat3x3MatMult(u1xd, XdinvT, tmp);
  mat3x3MatMultAdd(u0xd, XdinvzT, tmp);
  mat3x3TransMatMult(T, tmp, u1xd);

  // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
  // u0x = T^{T}*u0d*Xdinv*T
  mat3x3MatMult(u0x, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0x);

  mat3x3MatMult(u0xd, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0xd);

  // Compute the interpolation of the entries of the C matrix
  TacsScalar Cpt[9], Cptd[9];
  basis::template interpFields<9, 9>(pt, C, Cpt);
  basis::template interpFields<9, 9>(pt, Cd, Cptd);

  // Compute Ct = T^{T}*Cpt*T
  mat3x3TransMatMult(T, Cpt, tmp);
  mat3x3MatMult(tmp, T, Ct);

  mat3x3TransMatMult(T, Cptd, tmp);
  mat3x3MatMult(tmp, T, Ctd);

  return detXd;
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
  @param dCt Coefficients for Ct
  @param res The residual
  @param dd Residual intermediate for the director field
  @param dC Residual intermediate for the interpolated rotations
*/
template <int vars_per_node, class basis>
void addDispGradSens( const double pt[],
                      const TacsScalar T[],
                      const TacsScalar XdinvT[],
                      const TacsScalar XdinvzT[],
                      const TacsScalar du0x[],
                      const TacsScalar du1x[],
                      const TacsScalar dCt[],
                      TacsScalar res[],
                      TacsScalar dd[],
                      TacsScalar dC[] ){
  // Compute dCpt = T*dCt*T^{T}
  TacsScalar dCpt[9], tmp[9];
  mat3x3MatMult(T, dCt, tmp);
  mat3x3MatTransMult(tmp, T, dCpt);
  basis::template addInterpFieldsTranspose<9, 9>(pt, dCpt, dC);

  // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
  TacsScalar du0d[9];
  mat3x3MatTransMult(du1x, XdinvzT, tmp);
  mat3x3MatTransMultAdd(du0x, XdinvT, tmp);
  mat3x3MatMult(T, tmp, du0d);

  // Compute du1d = T*du1x*XdinvT^{T}
  TacsScalar du1d[9];
  mat3x3MatTransMult(du1x, XdinvT, tmp);
  mat3x3MatMult(T, tmp, du1d);

  // du0d = [du0xi; dd0]
  TacsScalar du0xi[6], dd0[3];
  extractFrame(du0d, du0xi, dd0);

  TacsScalar dd0xi[6];
  extractFrame(du1d, dd0xi);

  // Compute the director field and the gradient of the director
  // field at the specified point
  basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, dd0xi, dd);

  // Compute the gradient of the displacement solution at the quadrature points
  basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, du0xi, res);
}

/**
  Add/accumulate the contribution to the Jacobian matrix from the coefficients
  of u0x, u1x

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
void addDispGradHessian( const double pt[],
                         const TacsScalar T[],
                         const TacsScalar XdinvT[],
                         const TacsScalar XdinvzT[],
                         const TacsScalar d2u0x[],
                         const TacsScalar d2u1x[],
                         const TacsScalar d2u0xu1x[],
                         TacsScalar mat[],
                         TacsScalar d2d[],
                         TacsScalar d2du[] ){
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
  extractFrameSens(d2u0d, d2u0xi, d2d0, d2d0u0xi);

  TacsScalar d2d0xi[36];
  extractFrameSens(d2u1d, d2d0xi);

  TacsScalar d2d0xiu0xi[36], d2d0d0xi[18];
  extractFrameMixedSens(d2u0du1d, d2d0xiu0xi, d2d0d0xi);

  // df/(d(u0d)d(u1d))
  basis::template addInterpFieldsOuterProduct<3, 3>(pt, d2d0, d2d);
  basis::template addInterpGradOuterProduct<3, 3>(pt, d2d0xi, d2d);
  basis::template addInterpGradMixedOuterProduct<3, 3>(pt, d2d0d0xi, d2d0d0xi, d2d);

  // df/(d(d0)d(u0xi))
  basis::template addInterpGradMixedOuterProduct<3, 3>(pt, d2d0u0xi, NULL, d2du);
  basis::template addInterpGradOuterProduct<3, 3>(pt, d2d0xiu0xi, d2du);

  // Add the contribution to the Jacobian matrix
  basis::template addInterpGradOuterProduct<vars_per_node, 3>(pt, d2u0xi, mat);
}

/**
  Evaluate the tensorial components of the strain tensor at the specific
  quadrature point

  gty = [g11  g12  g13]
        [sym  g22  g23]
        [sym  sym  g33]

  As a result: gty[0] = g11, gty[1] = g12, gty[2] = g13, gty[3] = g22
  and gty[4] = g23, with gty[5] = 0.0

  @param pt The quadrature point
  @param ety The strain computed at the tying points
  @param gty The interpolated tying strain
*/
template <class basis>
void interpTyingStrain( const double pt[],
                        const TacsScalar ety[],
                        TacsScalar gty[] ){
  // Set the values into the strain tensor
  const int index[] = {0, 3, 1, 4, 2};
  const int num_tying_fields = 5;
  for ( int field = 0; field < num_tying_fields; field++ ){
    gty[index[field]] = basis::interpTying(field, pt, ety);
    ety += basis::getNumTyingPoints(field);
  }
  gty[5] = 0.0;
}

/**
  Add the derivative of the tying strain to the residual

  @param pt The quadrature point
  @param dgty The derivative of the interpolated strain
  @param dety The output derivative of the strain at the tying points
*/
template <class basis>
void addInterpTyingStrainTranspose( const double pt[],
                                    const TacsScalar dgty[],
                                    TacsScalar dety[] ){
  // Set the values into the strain tensor
  const int index[] = {0, 3, 1, 4, 2};
  const int num_tying_fields = 5;
  for ( int field = 0; field < num_tying_fields; field++ ){
    basis::addInterpTyingTranspose(field, pt, dgty[index[field]], dety);
    dety += basis::getNumTyingPoints(field);
  }
}

/**
  Add the second derivative of the tying strain at the tying points

  @param pt The quadrature point
  @param d2gty The second derivative of the interpolated strain
  @param d2ety The second derivatives of the strain at the tying points
*/
template <class basis>
void addInterpTyingStrainHessian( const double pt[],
                                  const TacsScalar d2gty[],
                                  TacsScalar d2ety[] ){
  // Set the values into the strain tensor
  const int index[] = {0, 3, 1, 4, 2};
  const int num_strains = 6;
  const int num_tying_fields = 5;
  for ( int field1 = 0; field1 < num_tying_fields; field1++ ){
    for ( int field2 = 0; field2 < num_tying_fields; field2++ ){
      TacsScalar value = d2gty[num_strains*index[field1] + index[field2]];
      basis::addInterpTyingOuterProduct(field1, field2, pt, value, d2ety);
      d2ety += basis::getNumTyingPoints(field1)*basis::getNumTyingPoints(field2);
    }
  }
}

/*
  Test the implementation of the shell terms for a given basis
*/
template <int vars_per_node, class basis>
int TacsTestShellUtilities( double dh=1e-7,
                            int test_print_level=2,
                            double test_fail_atol=1e-5,
                            double test_fail_rtol=1e-5 ){
  const int size = vars_per_node*basis::NUM_NODES;
  const int usize = 3*basis::NUM_NODES;
  const int dsize = 3*basis::NUM_NODES;
  const int csize = 9*basis::NUM_NODES;
  const int xsize = 3*basis::NUM_NODES;

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

  TacsScalar vars[size], C[csize], d[dsize];
  TacsGenerateRandomArray(vars, size);
  TacsGenerateRandomArray(C, csize);
  TacsGenerateRandomArray(d, dsize);

  // Generate random perturbations for the linear terms
  TacsScalar du0x[9], du1x[9], dCt[9];
  TacsGenerateRandomArray(du0x, 9);
  TacsGenerateRandomArray(du1x, 9);
  TacsGenerateRandomArray(dCt, 9);

  // Generate the derivative terms
  TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
  TacsScalar d2Ct[81], d2Ctu0x[81];
  TacsGenerateRandomArray(d2u0x, 81);
  TacsGenerateRandomArray(d2u1x, 81);
  TacsGenerateRandomArray(d2u0xu1x, 81);
  TacsGenerateRandomArray(d2Ct, 81);
  TacsGenerateRandomArray(d2Ctu0x, 81);

  // Symmetrize the random matrices
  for ( int i = 0; i < 9; i++ ){
    for ( int j = i+1; j < 9; j++ ){
      d2u0x[9*i + j] = d2u0x[9*j + i];
      d2u1x[9*i + j] = d2u1x[9*j + i];
      d2Ct[9*i + j] = d2Ct[9*j + i];
    }
  }

  TacsScalar res[size], mat[size*size];
  memset(res, 0, size*sizeof(TacsScalar));
  memset(mat, 0, size*size*sizeof(TacsScalar));

  TacsScalar dd[dsize], d2d[dsize*dsize], d2du[usize*dsize];
  memset(dd, 0, dsize*sizeof(TacsScalar));
  memset(d2d, 0, dsize*dsize*sizeof(TacsScalar));
  memset(d2du, 0, usize*dsize*sizeof(TacsScalar));

  TacsScalar dC[csize];
  memset(dC, 0, csize*sizeof(TacsScalar));

  TacsScalar XdinvT[9], XdinvzT[9];
  TacsScalar u0x[9], u1x[9], Ct[9];
  computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                        XdinvT, XdinvzT, u0x, u1x, Ct);

  addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                        du0x, du1x, dCt, res, dd, dC);
  addDispGradHessian<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                           d2u0x, d2u1x, d2u0xu1x,
                                           mat, d2d, d2du);

  // Now, check the result
  TacsScalar fdmat[size*size], fddu[dsize*usize];
  for ( int k = 0; k < size; k++ ){
    TacsScalar varst[size];
    memcpy(varst, vars, size*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[k] = vars[k] + TacsScalar(0.0, dh);
#else
    varst[k] = vars[k] + dh;
#endif // TACS_USE_COMPLEX

    // Compute the pertubation
    TacsScalar u0xt[9], u1xt[9], Ctt[9];
    computeDispGrad<vars_per_node, basis>(pt, Xpts, varst, fn, C, d, Xxi, n0, T,
                                          XdinvT, XdinvzT, u0xt, u1xt, Ctt);

    // d2u0xu1x[9*i + j] = p2f/(p(u0x[i]) p(u1x[j]))
    // d2Ctu0x[9*i + j] = p2f/(p(Ct[i]) p(u0x[j]))

    // Compute the perturbed values based on the outputs
    TacsScalar du0xt[9], du1xt[9];
    for ( int i = 0; i < 9; i++ ){
      du0xt[i] = du0x[i];
      du1xt[i] = du1x[i];
      for ( int j = 0; j < 9; j++ ){
        du0xt[i] += d2u0x[9*i + j]*(u0xt[j] - u0x[j]) +
                    d2u0xu1x[9*i + j]*(u1xt[j] - u1x[j]);

        du1xt[i] += d2u1x[9*i + j]*(u1xt[j] - u1x[j]) +
                    d2u0xu1x[9*j + i]*(u0xt[j] - u0x[j]);
      }
    }

    TacsScalar rest[size], ddt[dsize], dCvt[csize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));
    memset(dCvt, 0, csize*sizeof(TacsScalar));
    addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                          du0xt, du1xt, dCt, rest, ddt, dCvt);

    // Place the result in the arrays...
    for ( int j = 0; j < size; j++ ){
#ifdef TACS_USE_COMPLEX
      fdmat[k + size*j] = TacsImagPart(rest[j])/dh;
#else
      fdmat[k + size*j] = (rest[j] - res[j])/dh;
#endif // TACS_USE_COMPLEX
    }

    if (k % vars_per_node < 3){
      // Compute the u-index
      int index = 3*(k / vars_per_node) + k % vars_per_node;

      for ( int j = 0; j < dsize; j++ ){
#ifdef TACS_USE_COMPLEX
        fddu[index + usize*j] = TacsImagPart(ddt[j])/dh;
#else
        fddu[index + usize*j] = (ddt[j] - dd[j])/dh;
#endif // TACS_USE_COMPLEX
      }
    }
  }

  // Variables to store the max error and indices
  int max_err_index, max_rel_index;
  double max_err, max_rel;

  // Keep track of the failure flag
  int fail = 0;

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size*size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size*size, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size*size);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the error
  max_err = TacsGetMaxError(d2du, fddu, usize*dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2du, fddu, usize*dsize, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. d and vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2du", d2du, fddu, usize*dsize);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Now, check the result
  TacsScalar fddd[dsize*dsize];
  for ( int k = 0; k < dsize; k++ ){
    TacsScalar dt[dsize];
    memcpy(dt, d, dsize*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    dt[k] = d[k] + TacsScalar(0.0, dh);
#else
    dt[k] = d[k] + dh;
#endif // TACS_USE_COMPLEX

    // Compute the pertubation
    TacsScalar u0xt[9], u1xt[9], Ctt[9];
    computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, dt, Xxi, n0, T,
                                          XdinvT, XdinvzT, u0xt, u1xt, Ctt);

    // d2u0xu1x[9*i + j] = p2f/(p(u0x[i]) p(u1x[j]))
    // d2Ctu0x[9*i + j] = p2f/(p(Ct[i]) p(u0x[j]))

    // Compute the perturbed values based on the outputs
    TacsScalar du0xt[9], du1xt[9];
    for ( int i = 0; i < 9; i++ ){
      du0xt[i] = du0x[i];
      du1xt[i] = du1x[i];
      for ( int j = 0; j < 9; j++ ){
        du0xt[i] += d2u0x[9*i + j]*(u0xt[j] - u0x[j]) +
                    d2u0xu1x[9*i + j]*(u1xt[j] - u1x[j]);

        du1xt[i] += d2u1x[9*i + j]*(u1xt[j] - u1x[j]) +
                    d2u0xu1x[9*j + i]*(u0xt[j] - u0x[j]);
      }
    }

    TacsScalar rest[size], ddt[dsize], dCvt[csize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));
    memset(dCvt, 0, csize*sizeof(TacsScalar));
    addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                          du0xt, du1xt, dCt, rest, ddt, dCvt);

    for ( int j = 0; j < dsize; j++ ){
#ifdef TACS_USE_COMPLEX
      fddd[k + dsize*j] = TacsImagPart(ddt[j])/dh;
#else
      fddd[k + dsize*j] = (ddt[j] - dd[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2d, fddd, dsize*dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d, fddd, dsize*dsize, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2d", d2d, fddd, dsize*dsize);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif // TACS_SHELL_UTILITIES_H