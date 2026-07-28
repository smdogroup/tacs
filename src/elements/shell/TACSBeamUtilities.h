#ifndef TACS_BEAM_UTILITIES_H
#define TACS_BEAM_UTILITIES_H

#include <cstring>

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
void TacsBeamComputeNodeNormals(const TacsScalar Xpts[], const A2D::Vec3 &axis,
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
void TacsBeamAddNodeNormalsSens(const TacsScalar Xpts[], const A2D::Vec3 &axis,
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

/*
  Add the cross-Hessian contribution between the two director fields (d1, d2)
  to the element Jacobian. Shell elements have a single director, so there is
  no precedent: director::addDirectorJacobian handles each director's own
  self-Hessian and its coupling to the translational DOFs, but has no slot for
  a cross term between the two directors.

  Beam's two directors share the same per-node rotation vector q but use
  different reference normals fn1/fn2. For TACSLinearizedRotation the director
  map is d = crossProduct(q, t) = -skew(t)*q, so d(d_node)/d(q_node) =
  -skew(t_node), a fixed per-node congruence factor. Propagating the
  director-space cross-Hessian d2d1d2 through this congruence on both sides
  gives, per node pair (i, j):

    mat_block(q_i, q_j) -= skew(fn1_i) * d2d1d2[i,j] * skew(fn2_j)

  computed via mat3x3SkewMatSkewTransform (TACSElementAlgebra.h). Both
  directors share q, so this accumulates additively into the same
  (rotational, rotational) sub-block that the two director::addDirectorJacobian
  calls also write into.

  Exact only for TACSLinearizedRotation, where the director-to-rotation map is
  a fixed linear congruence. TACSQuadraticRotation/TACSQuaternionRotation make
  this map state-dependent and would need an additional cross-director
  curvature correction.

  @param vars The full variable vector (unused for TACSLinearizedRotation's
  fixed congruence factor; retained for a possible nonlinear-director extension)
  @param fn1 The first reference normal direction at each node
  @param fn2 The second reference normal direction at each node
  @param d2d1d2 The (d1-index, d2-index) cross-Hessian, dsize x dsize
  row-major, already alpha/gamma-scaled by the caller
  @param mat The element Jacobian matrix
*/
template <int vars_per_node, int offset, int num_nodes>
void TacsBeamAddCrossDirectorJacobian(const TacsScalar vars[],
                                      const TacsScalar fn1[],
                                      const TacsScalar fn2[],
                                      const TacsScalar d2d1d2[],
                                      TacsScalar mat[]) {
  const int dsize = 3 * num_nodes;
  const int nvars = vars_per_node * num_nodes;

  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      TacsScalar block[9];
      for (int p = 0; p < 3; p++) {
        for (int q = 0; q < 3; q++) {
          block[3 * p + q] = d2d1d2[dsize * (3 * i + p) + (3 * j + q)];
        }
      }

      TacsScalar D[9];
      mat3x3SkewMatSkewTransform(&fn1[3 * i], block, &fn2[3 * j], D);

      for (int p = 0; p < 3; p++) {
        for (int q = 0; q < 3; q++) {
          int row = offset + vars_per_node * i + p;
          int col = offset + vars_per_node * j + q;
          mat[row * nvars + col] -= D[3 * p + q];

          // Reusing the SAME D[3*p+q] value (rather than a separately
          // computed transpose) for the (j,i) slot is valid only because
          // this model's d2d1d2 3x3 node-pair blocks are always symmetric
          // (block[p,q] == block[q,p]): the static contribution is a
          // rank-1 t1 (tangent) outer product with itself, and the
          // dynamics contribution is rho[5]*I3 -- both symmetric 3x3
          // blocks. A kinematics change that made d2d1d2's node-pair blocks
          // non-symmetric would need a separate transposed computation
          // here, not this shortcut.
          int rowT = offset + vars_per_node * j + q;
          int colT = offset + vars_per_node * i + p;
          mat[rowT * nvars + colT] -= D[3 * p + q];
        }
      }
    }
  }
}

/*
  Zero the second-order (.xp/.xh, .Ap/.Ah) fields on every second-order A2D
  node before a per-DOF hforward/hreverse sweep iteration. These fields are
  "+=" accumulators that are not auto-reset, and most are not overwritten by
  each iteration's own seed step (e.g. the translational sweep seeds only
  u0xi.xp), so stale values from a previous iteration must be cleared here.

  @param u0xi The beam-axis displacement gradient (du0/dxi) node
  @param d01 The first director node
  @param d02 The second director node
  @param d01xi The parametric gradient of the first director node
  @param d02xi The parametric gradient of the second director node
  @param u0d The displacement gradient matrix node
  @param u0dXdinvT The displacement gradient matrix times XdinvT node
  @param u0x The displacement gradient in the local frame node
  @param d1t The first director rate before the frame transform node
  @param d1x The first director gradient in the local frame node
  @param d2t The second director rate before the frame transform node
  @param d2x The second director gradient in the local frame node
*/
static inline void TacsBeamZeroSecondOrderNodes(
    A2D::ADVec3 &u0xi, A2D::ADVec3 &d01, A2D::ADVec3 &d02, A2D::ADVec3 &d01xi,
    A2D::ADVec3 &d02xi, A2D::ADMat3x3 &u0d, A2D::ADMat3x3 &u0dXdinvT,
    A2D::ADMat3x3 &u0x, A2D::ADVec3 &d1t, A2D::ADVec3 &d1x, A2D::ADVec3 &d2t,
    A2D::ADVec3 &d2x) {
  memset(u0xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(u0xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(d01.xp, 0, 3 * sizeof(TacsScalar));
  memset(d01.xh, 0, 3 * sizeof(TacsScalar));
  memset(d02.xp, 0, 3 * sizeof(TacsScalar));
  memset(d02.xh, 0, 3 * sizeof(TacsScalar));
  memset(d01xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(d01xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(d02xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(d02xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(u0d.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0d.Ah, 0, 9 * sizeof(TacsScalar));
  memset(u0dXdinvT.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0dXdinvT.Ah, 0, 9 * sizeof(TacsScalar));
  memset(u0x.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0x.Ah, 0, 9 * sizeof(TacsScalar));
  memset(d1t.xp, 0, 3 * sizeof(TacsScalar));
  memset(d1t.xh, 0, 3 * sizeof(TacsScalar));
  memset(d1x.xp, 0, 3 * sizeof(TacsScalar));
  memset(d1x.xh, 0, 3 * sizeof(TacsScalar));
  memset(d2t.xp, 0, 3 * sizeof(TacsScalar));
  memset(d2t.xh, 0, 3 * sizeof(TacsScalar));
  memset(d2x.xp, 0, 3 * sizeof(TacsScalar));
  memset(d2x.xh, 0, 3 * sizeof(TacsScalar));
}

/*
  Contract the fixed, per-quadrature-point material Hessian blocks
  against the forward-propagated seed direction (u0xp/d1xp/d2xp)
  to produce the Hessian-vector product at the u0x/d1x/d2x
  nodes (u0xh/d1xh/d2xh), which the following hreverse() sweep propagates back
  to the vars-space DOFs. This is the "H*p" step of the per-DOF sweep: a pure
  linear-algebra contraction, not a differentiation. All array shapes and
  indices are row-major and match model::evalStrainHessian's layout, e.g.
  d2u0xd1x[3*i+j] = d^2/d(u0x_i)d(d1x_j).

  @param d2u0x The u0x self-Hessian block (9x9, row-major)
  @param d2d1x The d1x self-Hessian block (3x3, row-major)
  @param d2d2x The d2x self-Hessian block (3x3, row-major)
  @param d2u0xd1x The mixed u0x-d1x Hessian block (9x3, row-major)
  @param d2u0xd2x The mixed u0x-d2x Hessian block (9x3, row-major)
  @param d2d1xd2x The mixed d1x-d2x Hessian block (3x3, row-major)
  @param u0xp The forward-propagated u0x seed direction (input)
  @param d1xp The forward-propagated d1x seed direction (input)
  @param d2xp The forward-propagated d2x seed direction (input)
  @param u0xh The u0x Hessian-vector product, accumulated (output)
  @param d1xh The d1x Hessian-vector product, accumulated (output)
  @param d2xh The d2x Hessian-vector product, accumulated (output)
*/
static inline void TacsBeamContractStrainHessian(
    const TacsScalar d2u0x[81], const TacsScalar d2d1x[9],
    const TacsScalar d2d2x[9], const TacsScalar d2u0xd1x[27],
    const TacsScalar d2u0xd2x[27], const TacsScalar d2d1xd2x[9],
    const TacsScalar u0xp[9], const TacsScalar d1xp[3],
    const TacsScalar d2xp[3], TacsScalar u0xh[9], TacsScalar d1xh[3],
    TacsScalar d2xh[3]) {
  for (int i = 0; i < 9; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 9; j++) {
      val += d2u0x[9 * i + j] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2u0xd1x[3 * i + j] * d1xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2u0xd2x[3 * i + j] * d2xp[j];
    }
    u0xh[i] += val;
  }

  for (int i = 0; i < 3; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 3; j++) {
      val += d2d1x[3 * i + j] * d1xp[j];
    }
    for (int j = 0; j < 9; j++) {
      val += d2u0xd1x[3 * j + i] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2d1xd2x[3 * i + j] * d2xp[j];
    }
    d1xh[i] += val;
  }

  for (int i = 0; i < 3; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 3; j++) {
      val += d2d2x[3 * i + j] * d2xp[j];
    }
    for (int j = 0; j < 9; j++) {
      val += d2u0xd2x[3 * j + i] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2d1xd2x[3 * j + i] * d1xp[j];
    }
    d2xh[i] += val;
  }
}

#endif  // TACS_BEAM_UTILITIES_H
