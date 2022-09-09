#ifndef TACS_BEAM_UTILITIES_H
#define TACS_BEAM_UTILITIES_H

#include "TACSElementAlgebra.h"
#include "a2d.h"

/*
  Compute the frame normals at each node location

  @param Xpts The node locations for the elements
  @param axis The coordinates of the reference axis
  @param fn1 The first normal direction
  @param fn2 The second normal direction
*/
template <class basis>
void TacsBeamComputeNodeNormals(const TacsScalar Xpts[], const A2D::Vec3& axis,
                                TacsScalar fn1[], TacsScalar fn2[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    A2D::Vec3 X0xi;
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);

    // Normalize the first direction.
    A2D::Vec3 t1;
    A2D::Vec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::Vec3 t2_dir;
    A2D::Scalar dot;
    A2D::Vec3Dot dott1(axis, t1, dot);
    A2D::Vec3Axpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Compute the t2 direction
    A2D::Vec3 t2;
    A2D::Vec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::Vec3 t3;
    A2D::Vec3CrossProduct cross(t1, t2, t3);

    fn1[0] = t2.x[0];
    fn1[1] = t2.x[1];
    fn1[2] = t2.x[2];

    fn2[0] = t3.x[0];
    fn2[1] = t3.x[1];
    fn2[2] = t3.x[2];

    fn1 += 3;
    fn2 += 3;
  }
}

/*
  Compute the frame normals at each node location

  @param Xpts The node locations for the elements
  @param axis The coordinates of the reference axis
  @param fn1 The first normal direction
  @param fn2 The second normal direction
*/
template <class basis>
void TacsBeamAddNodeNormalsSens(const TacsScalar Xpts[], const A2D::Vec3& axis,
                                const TacsScalar dfn1[],
                                const TacsScalar dfn2[], TacsScalar dXpts[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    A2D::ADVec3 X0xi;
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);

    // Normalize the first direction.
    A2D::ADVec3 t1;
    A2D::ADVec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::ADVec3 t2_dir;
    A2D::ADScalar dot;
    A2D::Vec3ADVecDot dott1(axis, t1, dot);
    A2D::ADVec3VecADScalarAxpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Compute the t2 direction
    A2D::ADVec3 t2(NULL, dfn1);
    A2D::ADVec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::ADVec3 t3(NULL, dfn2);
    A2D::ADVec3CrossProduct cross(t1, t2, t3);

    cross.reverse();
    normalizet2.reverse();
    axpy.reverse();
    dott1.reverse();
    normalizet1.reverse();

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, X0xi.xd, dXpts);

    dfn1 += 3;
    dfn2 += 3;
  }
}

// /**
//   Compute the displacement gradient of the constant and through-thickness
//   rate of change of the displacements.

//   @param pt The parametric point
//   @param Xpts The node locations for the element
//   @param vars The element variables
//   @param fn The frame normal directions at each node
//   @param d The director field at each node
//   @param Xxi The in-plane coordinate derivatives
//   @param n0 The interpolated frame normal direction
//   @param T The transformation to local coordinates
//   @param XdinvT Product of inverse of the Jacobian trans. and T
//   @param XdinvzT Product of z-derivative of Jac. trans. inv. and T
//   @param u0x Derivative of the displacement in the local x coordinates
//   @param u1x Derivative of the through-thickness disp. in local x coordinates
// */
// template <int vars_per_node, class basis>
// TacsScalar TacsBeamComputeDispGrad( const double pt[],
//                                     const TacsScalar Xpts[],
//                                     const TacsScalar vars[],
//                                     const TacsScalar fn1[],
//                                     const TacsScalar fn2[],
//                                     const TacsScalar d1[],
//                                     const TacsScalar d2[],
//                                     const TacsScalar Xxi[],
//                                     const TacsScalar n1[],
//                                     const TacsScalar n2[],
//                                     const TacsScalar T[],
//                                     TacsScalar XdinvT[],
//                                     TacsScalar Xdinvz1T[],
//                                     TacsScalar Xdinvz2T[],
//                                     TacsScalar u0x[],
//                                     TacsScalar d1x[],
//                                     TacsScalar d2x[] ){
//   // Assemble the reference frame
//   TacsScalar Xd[9];
//   TacsShellAssembleFrame(Xxi, n1, n2, Xd);

//   // Compute the inverse of the 3x3 Jacobian transformation
//   TacsScalar Xdinv[9];
//   TacsScalar detXd = inv3x3(Xd, Xdinv);

//   // Assemble the derivative of the reference frame
//   TacsScalar Ud[9];
//   TacsShellAssembleFrame(u0xi, d1, d2, Ud);

//   // Compute n,xi = [dn/dxi1; dn/dxi2]
//   TacsScalar n1xi[3], n2xi[3];
//   basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi);
//   basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi);

//   // Assemble the terms Xd = [Xxi; n1; n2] and Xdz
//   TacsScalar Xdz1[9], Xdz2[9];
//   TacsShellAssembleFrame(n1xi, Xdz1);
//   TacsShellAssembleFrame(n2xi, Xdz2);

//   // Compute negXdinvXdz = -Xdinv*Xdz
//   TacsScalar negXdinvXdz1[9], negXdinvXdz2[9];
//   mat3x3MatMult(Xdinv, Xdz1, negXdinvXdz1);
//   for ( int i = 0; i < 9; i++ ){
//     negXdinvXdz1[i] *= -1.0;
//   }
//   mat3x3MatMult(Xdinv, Xdz2, negXdinvXdz2);
//   for ( int i = 0; i < 9; i++ ){
//     negXdinvXdz2[i] *= -1.0;
//   }

//   // Compute XdinvT = Xdinv*T
//   mat3x3MatMult(Xdinv, T, XdinvT);

//   // Compute Xdinvz = -Xdinv*Xdz*Xdinv*T
//   mat3x3MatMult(negXdinvXdz1, XdinvT, Xdinvz1T);
//   mat3x3MatMult(negXdinvXdz2, XdinvT, Xdinvz2T);

//   // Compute the director field and the gradient of the director
//   // field at the specified point
//   TacsScalar d01[3], d02[3], d1xi[3], d2xi[3];
//   basis::template interpFields<3, 3>(pt, d1, d01);
//   basis::template interpFields<3, 3>(pt, d2, d02);
//   basis::template interpFieldsGrad<3, 3>(pt, d1, d1xi);
//   basis::template interpFieldsGrad<3, 3>(pt, d2, d2xi);

//   // d1x = T^{T}*(d1xi*e1^{T} + Ur*z1Xdinv)*T*e1
//   TacsScalar scale = Xdinv[0]*T[0] + Xdinv[1]*T[3] + Xdinv[2]*T[6];
//   tmp[0] = T[0];
//   tmp[1] = T[3];
//   tmp[2] = T[6];
//   mat3x3Mult(z1Xdinv, tmp, d1x); // tmp = z1Xdinv*T*e1
//   mat3x3Mult(u0xi, d1x, tmp); // tmp = Ur*z1Xdinv*T*e1
//   vec3Axpy(scale, d1a, tmp);
//   mat3x3MultTrans(T, tmp, Td1a);

//   // Td2a = T^{T}*d2a*e1^{T}*Xdinv*T*e1
//   TacsScalar Td2a[3], z2Te1[3];
//   tmp[0] = T[0];
//   tmp[1] = T[3];
//   tmp[2] = T[6];
//   mat3x3Mult(z2Xdinv, tmp, z2Te1); // tmp = z1Xdinv*T*e1
//   mat3x3Mult(Ur, z2Te1, tmp); // tmp = Ur*z1Xdinv*T*e1
//   vec3Axpy(S[0], d2a, tmp);
//   mat3x3MultTrans(T, tmp, Td2a);

//   // Compute the gradient of the displacement solution at the quadrature
//   points TacsScalar u0[3], u0xi[3]; basis::template
//   interpFields<vars_per_node, 3>(pt, vars, u0); basis::template
//   interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi);

//   // Compute the derivative u0,x
//   TacsShellAssembleFrame(u0xi, d01, d01, u0x); // Use u0x to store [u0,xi;
//   d1, d2]

//   // d1x = T^{T}*d1*XdinvT + T^{T}*u0*XdinvzT
//   // TacsScalar tmp[9];
//   // TacsShellAssembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
//   // mat3x3MatMult(u1x, XdinvT, tmp);
//   // mat3x3MatMultAdd(u0x, XdinvzT, tmp);
//   // mat3x3TransMatMult(T, tmp, u1x);

//   // // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
//   // // u0x = T^{T}*u0d*Xdinv*T
//   // mat3x3MatMult(u0x, XdinvT, tmp);
//   // mat3x3TransMatMult(T, tmp, u0x);

//   return detXd;
// }

/*
  Test the implementation of the shell terms for a given basis
*/
template <int vars_per_node, class basis>
int TacsTestBeamUtilities(double dh = 1e-7, int test_print_level = 2,
                          double test_fail_atol = 1e-5,
                          double test_fail_rtol = 1e-5) {
  int fail = 0;
  /*
  const int size = vars_per_node*basis::NUM_NODES;
  const int usize = 3*basis::NUM_NODES;
  const int dsize = 3*basis::NUM_NODES;
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
  for ( int i = 0; i < 9; i++ ){
    for ( int j = i+1; j < 9; j++ ){
      d2u0x[9*i + j] = d2u0x[9*j + i];
      d2u1x[9*i + j] = d2u1x[9*j + i];
    }
  }

  TacsScalar res[size], mat[size*size];
  memset(res, 0, size*sizeof(TacsScalar));
  memset(mat, 0, size*size*sizeof(TacsScalar));

  TacsScalar dd[dsize], d2d[dsize*dsize], d2du[usize*dsize];
  memset(dd, 0, dsize*sizeof(TacsScalar));
  memset(d2d, 0, dsize*dsize*sizeof(TacsScalar));
  memset(d2du, 0, usize*dsize*sizeof(TacsScalar));

  TacsScalar XdinvT[9], XdinvzT[9];
  TacsScalar u0x[9], u1x[9];
  TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, d, Xxi, n0,
T, XdinvT, XdinvzT, u0x, u1x);

  TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                 du0x, du1x, res, dd);
  TacsShellAddDispGradHessian<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
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
    TacsScalar u0xt[9], u1xt[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, varst, fn, d, Xxi,
n0, T, XdinvT, XdinvzT, u0xt, u1xt);

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

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                   du0xt, du1xt, rest, ddt);

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
    TacsScalar u0xt[9], u1xt[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, dt, Xxi,
n0, T, XdinvT, XdinvzT, u0xt, u1xt);

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

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                   du0xt, du1xt, rest, ddt);

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
  */

  return fail;
}

#endif  // TACS_BEAM_UTILITIES_H
