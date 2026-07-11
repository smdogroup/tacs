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

#endif  // TACS_BEAM_UTILITIES_H
