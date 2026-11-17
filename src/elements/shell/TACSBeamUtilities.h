#ifndef TACS_BEAM_UTILITIES_H
#define TACS_BEAM_UTILITIES_H

#include "TACSA2DUtilities.h"
#include "TACSElementAlgebra.h"

/*
  Compute the frame normals at each node location

  @param Xpts The node locations for the elements
  @param axis The coordinates of the reference axis
  @param fn1 The first normal direction
  @param fn2 The second normal direction
*/
template <class basis>
void TacsBeamComputeNodeNormals(const TacsScalar Xpts[],
                                const A2D::Vec<TacsScalar, 3> &axis,
                                TacsScalar fn1[], TacsScalar fn2[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    A2D::Vec<TacsScalar, 3> X0xi;
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));

    // Normalize the first direction.
    A2D::Vec<TacsScalar, 3> t1;
    A2D::VecNormalize(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    TacsScalar dot;
    A2D::VecDot(axis, t1, dot);
    A2D::Vec<TacsScalar, 3> t2_dir;
    A2D::VecSum(TacsScalar(1.0), axis, -dot, t1, t2_dir);

    // Compute the t2 direction
    A2D::Vec<TacsScalar, 3> t2;
    A2D::VecNormalize(t2_dir, t2);

    // Compute the n2 direction
    A2D::Vec<TacsScalar, 3> t3;
    A2D::VecCross(t1, t2, t3);

    fn1[0] = t2(0);
    fn1[1] = t2(1);
    fn1[2] = t2(2);

    fn2[0] = t3(0);
    fn2[1] = t3(1);
    fn2[2] = t3(2);

    fn1 += 3;
    fn2 += 3;
  }
}

/*
  Add the sensitivity of the frame normals at each node location

  @param Xpts The node locations for the elements
  @param axis The coordinates of the reference axis
  @param dfn1 The seed for the first normal direction
  @param dfn2 The seed for the second normal direction
  @param dXpts The derivative w.r.t. the node locations
*/
template <class basis>
void TacsBeamAddNodeNormalsSens(const TacsScalar Xpts[],
                                const A2D::Vec<TacsScalar, 3> &axis,
                                const TacsScalar dfn1[],
                                const TacsScalar dfn2[], TacsScalar dXpts[]) {
  for (int i = 0; i < basis::NUM_NODES; i++) {
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0xi;
    basis::template interpFieldsGrad<3, 3>(pt, Xpts,
                                           A2D::get_data(X0xi.value()));

    A2D::ADObj<TacsScalar> dot;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> t1, w, t2_dir, t2, t3;

    auto stack = A2D::MakeStack(
        // Normalize the first direction
        A2D::VecNormalize(X0xi, t1),
        // t2_dir = axis - dot(t1, axis) * t1
        A2D::VecDot(axis, t1, dot), A2D::VecScale(dot, t1, w),
        A2D::VecSum(TacsScalar(1.0), axis, TacsScalar(-1.0), w, t2_dir),
        // Compute the t2 direction
        A2D::VecNormalize(t2_dir, t2),
        // Compute the n2 direction
        A2D::VecCross(t1, t2, t3));

    // Seed the outputs and propagate the derivatives back to X0xi
    for (int j = 0; j < 3; j++) {
      t2.bvalue()(j) = dfn1[j];
      t3.bvalue()(j) = dfn2[j];
    }
    stack.reverse();

    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(X0xi.bvalue()), dXpts);

    dfn1 += 3;
    dfn2 += 3;
  }
}

#endif  // TACS_BEAM_UTILITIES_H
