#ifndef TACS_SHELL_ELEMENT_TRANSFORM_H
#define TACS_SHELL_ELEMENT_TRANSFORM_H

#include "TACSElementAlgebra.h"
#include "TACSObject.h"

/*
  Compute the transformation from the local coordinates
  to
*/
class TACSShellTransform : public TACSObject {
 public:
  /*
    Given the local shell element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform(const TacsScalar Xxi[], const TacsScalar n0[],
                                TacsScalar T[]) = 0;
};

class TACSShellNaturalTransform : public TACSShellTransform {
 public:
  TACSShellNaturalTransform() {}

  void computeTransform(const TacsScalar Xxi[], const TacsScalar n0[],
                        TacsScalar T[]) {
    TacsScalar n[3];
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];

    // Scale by the normal
    TacsScalar inv = 1.0 / sqrt(vec3Dot(n, n));
    vec3Scale(inv, n);

    TacsScalar t1[3];
    t1[0] = Xxi[0];
    t1[1] = Xxi[2];
    t1[2] = Xxi[4];

    TacsScalar d = vec3Dot(n, t1);
    t1[0] = t1[0] - d * n[0];
    t1[0] = t1[0] - d * n[0];
    t1[0] = t1[0] - d * n[0];

    inv = 1.0 / sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    TacsScalar t2[3];
    crossProduct(n, t1, t2);

    /*

    // Compute the transformation
    TacsScalar t1[3], t2[3];
    t1[0] = Xxi[0];
    t1[1] = Xxi[2];
    t1[2] = Xxi[4];

    t2[0] = Xxi[1];
    t2[1] = Xxi[3];
    t2[2] = Xxi[5];

    // Compute the normal direction
    TacsScalar n[3];
    crossProduct(t1, t2, n);

    // Normalize the normal direction
    TacsScalar invNorm = 1.0/sqrt(vec3Dot(n, n));
    vec3Scale(invNorm, n);

    // Normalize the 1-direction of the element
    TacsScalar inv = 1.0/sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    // Take the cross product to determine the 2-direction
    crossProduct(n, t1, t2);
    */

    // Set the components of the transformation
    T[0] = t1[0];
    T[3] = t1[1];
    T[6] = t1[2];

    T[1] = t2[0];
    T[4] = t2[1];
    T[7] = t2[2];

    T[2] = n[0];
    T[5] = n[1];
    T[8] = n[2];
  }
};

class TACSShellRefAxisTransform : public TACSShellTransform {
 public:
  TACSShellRefAxisTransform(const TacsScalar _axis[]) {
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];

    TacsScalar norm = sqrt(vec3Dot(axis, axis));
    TacsScalar invNorm = 0.0;
    if (norm != 0.0) {
      invNorm = 1.0 / norm;
    }
    vec3Scale(invNorm, axis);
  }

  void computeTransform(const TacsScalar Xxi[], const TacsScalar n0[],
                        TacsScalar T[]) {
    TacsScalar n[3];
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];

    // Scale by the normal
    TacsScalar inv = 1.0 / sqrt(vec3Dot(n, n));
    vec3Scale(inv, n);

    // Compute the dot product with
    TacsScalar an = vec3Dot(axis, n);

    // Take the component of the reference axis perpendicular
    // to the surface
    TacsScalar t1[3];
    t1[0] = axis[0] - an * n[0];
    t1[1] = axis[1] - an * n[1];
    t1[2] = axis[2] - an * n[2];

    // Normalize the new direction
    inv = 1.0 / sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    // Take the cross product to determine the 2-direction
    TacsScalar t2[3];
    crossProduct(n, t1, t2);

    /*
        // Compute the transformation
        TacsScalar t1[3], t2[3];
        t1[0] = Xxi[0];
        t1[1] = Xxi[2];
        t1[2] = Xxi[4];

        t2[0] = Xxi[1];
        t2[1] = Xxi[3];
        t2[2] = Xxi[5];
    */

    /*
    // Compute the transformation
    TacsScalar t1[3], t2[3];
    t1[0] = Xxi[0];
    t1[1] = Xxi[2];
    t1[2] = Xxi[4];

    t2[0] = Xxi[1];
    t2[1] = Xxi[3];
    t2[2] = Xxi[5];

    // Compute the normal direction
    TacsScalar n[3];
    crossProduct(t1, t2, n);

    // Normalize the normal direction
    TacsScalar invNorm = 1.0/sqrt(vec3Dot(n, n));
    vec3Scale(invNorm, n);

    // Compute the dot product with
    TacsScalar an = vec3Dot(axis, n);

    // Take the component of the reference axis perpendicular
    // to the surface
    t1[0] = axis[0] - an*n[0];
    t1[1] = axis[1] - an*n[1];
    t1[2] = axis[2] - an*n[2];

    // Normalize the new direction
    TacsScalar inv = 1.0/sqrt(vec3Dot(t1, t1));
    vec3Scale(inv, t1);

    // Take the cross product to determine the 2-direction
    crossProduct(n, t1, t2);
    */

    // Set the components of the transformation
    T[0] = t1[0];
    T[3] = t1[1];
    T[6] = t1[2];

    T[1] = t2[0];
    T[4] = t2[1];
    T[7] = t2[2];

    T[2] = n[0];
    T[5] = n[1];
    T[8] = n[2];
  }

 private:
  TacsScalar axis[3];
};

#endif  // TACS_SHELL_ELEMENT_TRANSFORM_H