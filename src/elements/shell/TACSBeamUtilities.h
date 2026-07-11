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
  Add the cross-Hessian contribution between the two director fields (d1,
  d2) to the element Jacobian matrix (SPEC.md sec 1.3.2, "Gap 2" -- the
  feature's first genuinely novel derivation; no shell precedent, since
  shell has only one director).

  director::addDirectorJacobian handles each director's own self-Hessian
  (and its coupling to the translational DOFs) but is called once per
  director and has no parameter slot for a genuine CROSS term between the
  two directors -- this closure fills that gap.

  Beam's two directors are both parametrized by the SAME per-node
  rotational sub-vector q (vars[offset..offset+NUM_PARAMETERS-1]), via
  DIFFERENT reference normals fn1/fn2 (TACSBeamElementModel::
  computeTyingStrain's own d1/d2 arrays are each produced by a separate
  director::computeDirectorRates call sharing the same q). For
  TACSLinearizedRotation (the only director class test_beam_element.py
  exercises), TACSDirector.h's computeDirectorRates computes
  d = crossProduct(q, t) = -skew(t)*q, so d(d_node)/d(q_node) = -skew(t_node)
  exactly -- a fixed, state-independent per-node congruence factor (no
  chain-rule re-derivation needed elsewhere; do not duplicate this
  derivative inline, it belongs to the director class).

  Propagating d2d1d2 (the (d1-index, d2-index) cross-Hessian in raw
  "director space", dsize x dsize row-major -- row = d1's compact index,
  col = d2's compact index -- fed by the per-DOF hforward/hreverse sweep
  in addJacobian, capturing the leakage that occurs because u0d's column
  assembly couples u0xi/d01/d02 together) through this congruence on BOTH
  sides gives, per node pair (i, j):

    mat_block(q_i, q_j) += (-skew(fn1_i))^T * d2d1d2[i,j] * (-skew(fn2_j))
                          = -skew(fn1_i) * d2d1d2[i,j] * skew(fn2_j)

  (the two minus signs from the congruence factors cancel one of the two
  transposes, leaving a single overall minus sign -- this is the identical
  sign structure TACSDirector.h's own addDirectorJacobian uses for the
  analogous SAME-director self term, "jac1 -= mat3x3SkewMatSkewTransform(
  ti, d, tj, ...)", TACSDirector.h:406-418, confirming this derivation
  against existing, tested code rather than a fresh, unverified sign
  convention). mat3x3SkewMatSkewTransform(a, B, c, D) computes
  D = skew(a)*B*skew(c) directly (TACSElementAlgebra.h) -- negate its
  output to get the formula above.

  Because both directors share the SAME per-node q, this block
  accumulates ADDITIVELY into the SAME (rotational, rotational) sub-block
  of mat[] that the two separate director::addDirectorJacobian calls (for
  d1/fn1 and d2/fn2) also write into -- this is expected, not a collision:
  the total Hessian w.r.t. (q_i, q_j) receives contributions from every
  path (via d1 and via d2), all correctly combined via +=.

  Location note: placed here (TACSBeamUtilities.h) rather than
  TACSDirector.h because this closure needs BOTH director classes' chain
  rules simultaneously (fn1's and fn2's), not one class's own method --
  TACSDirector.h's per-director-class structure has no natural home for a
  cross-class operation.

  Scope boundary (SPEC.md sec 1.3.2, explicit and documented, not a silent
  omission): this closure is EXACT only for TACSLinearizedRotation, where
  d(director)/d(rotation param) is a fixed linear congruence. For
  TACSQuadraticRotation/TACSQuaternionRotation this map is itself nonlinear
  in the state (TACSDirector.h:621-/1481-), so an additional
  "cross-director curvature" correction would be needed for exactness on
  those two director classes -- test_beam_element.py does not exercise
  them today. Do not silently extend this closure to other director
  classes without adding that correction.

  @param vars The full variable vector (unused for TACSLinearizedRotation's
  fixed congruence factor; retained in the signature for parity with a
  possible future nonlinear-director extension, per SPEC.md sec 1.3.2)
  @param fn1 The first reference normal direction at each node
  @param fn2 The second reference normal direction at each node
  @param d2d1d2 The (d1-index, d2-index) cross-Hessian, dsize x dsize
  row-major, already alpha/gamma-scaled by the caller (mirroring
  director::addDirectorJacobian's own "caller pre-scales d2d/d2du by
  alpha" contract)
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

          int rowT = offset + vars_per_node * j + q;
          int colT = offset + vars_per_node * i + p;
          mat[rowT * nvars + colT] -= D[3 * p + q];
        }
      }
    }
  }
}

#endif  // TACS_BEAM_UTILITIES_H
