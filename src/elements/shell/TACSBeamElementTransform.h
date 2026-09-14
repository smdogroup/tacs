#ifndef TACS_BEAM_ELEMENT_TRANSFORM_H
#define TACS_BEAM_ELEMENT_TRANSFORM_H

#include "TACSA2DUtilities.h"
#include "TACSElement.h"

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
  virtual A2D::Vec<TacsScalar, 3> &getRefAxis() = 0;
};

/*
  Compute the transformation
*/
class TACSBeamRefAxisTransform : public TACSBeamTransform {
 public:
  TACSBeamRefAxisTransform(const TacsScalar axis_dir[]) {
    A2D::Vec<TacsScalar, 3> axdir(axis_dir);
    A2D::VecNormalize(axdir, axis);
  }

  void computeTransform(const TacsScalar X0xi_vals[], TacsScalar Tvals[]) {
    // Normalize the first direction.
    A2D::Vec<TacsScalar, 3> X0xi(X0xi_vals);
    A2D::Vec<TacsScalar, 3> t1;
    A2D::VecNormalize(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    TacsScalar dot;
    A2D::VecDot(axis, t1, dot);
    A2D::Vec<TacsScalar, 3> t2_dir;
    A2D::VecSum(TacsScalar(1.0), axis, -dot, t1, t2_dir);

    // Check if ref axis is parallel to beam
    if (fabs(TacsRealPart(dot)) > 1.0 - SMALL_NUM) {
      fprintf(stderr,
              "TACSBeamRefAxisTransform: Error, user-provided reference axis "
              "is parallel to beam axis. "
              "Element behavior may be ill-conditioned.\n");
    }

    // Compute the t2 direction
    A2D::Vec<TacsScalar, 3> t2;
    A2D::VecNormalize(t2_dir, t2);

    // Compute the n2 direction
    A2D::Vec<TacsScalar, 3> t3;
    A2D::VecCross(t1, t2, t3);

    // Assemble the reference frame
    A2D::Mat<TacsScalar, 3, 3> T;
    A2D::MatFromThreeVec(t1, t2, t3, T);

    for (int i = 0; i < 9; i++) {
      Tvals[i] = A2D::get_data(T)[i];
    }
  }

  void addTransformSens(const TacsScalar X0xi_vals[], const TacsScalar dTvals[],
                        TacsScalar dX0xi[]) {
    // Normalize the first direction.
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0xi(X0xi_vals);

    A2D::ADObj<TacsScalar> dot;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> t1, w, t2_dir, t2, t3;
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> T;

    auto stack = A2D::MakeStack(
        // Normalize the first direction
        A2D::VecNormalize(X0xi, t1),
        // t2_dir = axis - dot(t1, axis) * t1
        A2D::VecDot(axis, t1, dot), A2D::VecScale(dot, t1, w),
        A2D::VecSum(TacsScalar(1.0), axis, TacsScalar(-1.0), w, t2_dir),
        // Compute the t2 direction
        A2D::VecNormalize(t2_dir, t2),
        // Compute the n2 direction
        A2D::VecCross(t1, t2, t3),
        // Assemble the reference frame
        A2D::MatFromThreeVec(t1, t2, t3, T));

    // Set the seeds for T, then reverse the operations to get the
    // derivative w.r.t. X0xi
    for (int i = 0; i < 9; i++) {
      A2D::get_data(T.bvalue())[i] = dTvals[i];
    }
    stack.reverse();

    for (int i = 0; i < 3; i++) {
      dX0xi[i] += X0xi.bvalue()(i);
    }
  }
  A2D::Vec<TacsScalar, 3> &getRefAxis() { return axis; }

  void getRefAxis(TacsScalar _axis[]) {
    _axis[0] = axis(0);
    _axis[1] = axis(1);
    _axis[2] = axis(2);
  }

 private:
  A2D::Vec<TacsScalar, 3> axis;
  /* Tolerance for colinearity test in between beam axis and ref axis */
  const double SMALL_NUM = 1e-8;
};

#endif  // TACS_BEAM_ELEMENT_TRANSFORM_H
