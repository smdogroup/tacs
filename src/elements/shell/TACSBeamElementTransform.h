#ifndef TACS_BEAM_ELEMENT_TRANSFORM_H
#define TACS_BEAM_ELEMENT_TRANSFORM_H

#include "TACSElement.h"
#include "a2d.h"

/*
  Compute the transformation from the local coordinates
*/
class TACSBeamTransform : public TACSObject {
 public:
  /*
    Given the local beam element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform(const TacsScalar Xxi[], TacsScalar T[]) = 0;
  virtual void addTransformSens(const TacsScalar X0xi_vals[],
                                const TacsScalar dTvals[],
                                TacsScalar dX0xi[]) = 0;
  virtual A2D::Vec3 &getRefAxis() = 0;
};

/*
  Compute the transformation
*/
class TACSBeamRefAxisTransform : public TACSBeamTransform {
 public:
  TACSBeamRefAxisTransform(const TacsScalar axis_dir[]) {
    A2D::Vec3 axdir(axis_dir);
    A2D::Vec3Normalize normalize(axdir, axis);
  }

  void computeTransform(const TacsScalar X0xi_vals[], TacsScalar Tvals[]) {
    // Normalize the first direction.
    A2D::Vec3 X0xi(X0xi_vals);
    A2D::Vec3 t1;
    A2D::Vec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::Vec3 t2_dir;
    A2D::Scalar dot;
    A2D::Vec3Dot dott1(axis, t1, dot);
    A2D::Vec3Axpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Check if ref axis is parallel to beam
    if (abs(TacsRealPart(dot.value)) > 1.0 - SMALL_NUM) {
      fprintf(stderr,
              "TACSBeamRefAxisTransform: Error, user-provided reference axis "
              "is parallel to beam axis. "
              "Element behavior may be ill-conditioned.\n");
    }

    // Compute the t2 direction
    A2D::Vec3 t2;
    A2D::Vec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::Vec3 t3;
    A2D::Vec3CrossProduct cross(t1, t2, t3);

    // Assemble the reference frame
    A2D::Mat3x3 T;
    A2D::Mat3x3FromThreeVec3 assembleT(t1, t2, t3, T);

    for (int i = 0; i < 9; i++) {
      Tvals[i] = T.A[i];
    }
  }

  void addTransformSens(const TacsScalar X0xi_vals[], const TacsScalar dTvals[],
                        TacsScalar dX0xi[]) {
    // Normalize the first direction.
    A2D::ADVec3 X0xi(X0xi_vals);
    A2D::ADVec3 t1;
    A2D::ADVec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::ADVec3 t2_dir;
    A2D::ADScalar dot;
    A2D::Vec3ADVecDot dott1(axis, t1, dot);
    A2D::ADVec3VecADScalarAxpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Compute the t2 direction
    A2D::ADVec3 t2;
    A2D::ADVec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::ADVec3 t3;
    A2D::ADVec3CrossProduct cross(t1, t2, t3);

    // Assemble the referece frame
    A2D::ADMat3x3 T(NULL, dTvals);  // Set the seeds for T
    A2D::ADMat3x3FromThreeADVec3 assembleT(t1, t2, t3, T);

    // Reverse the operations to get the derivative w.r.t. X0
    assembleT.reverse();
    cross.reverse();
    normalizet2.reverse();
    axpy.reverse();
    dott1.reverse();
    normalizet1.reverse();

    for (int i = 0; i < 3; i++) {
      dX0xi[i] += X0xi.xd[i];
    }
  }
  A2D::Vec3 &getRefAxis() { return axis; }

  void getRefAxis(TacsScalar _axis[]) {
    _axis[0] = axis.x[0];
    _axis[1] = axis.x[1];
    _axis[2] = axis.x[2];
  }

 private:
  A2D::Vec3 axis;
  /* Tolerance for colinearity test in between beam axis and ref axis */
  const double SMALL_NUM = 1e-8;
};

#endif  // TACS_BEAM_ELEMENT_TRANSFORM_H