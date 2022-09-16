#ifndef TACS_SPRING_ELEMENT_TRANSFORM_H
#define TACS_SPRING_ELEMENT_TRANSFORM_H

#include "TACSElementAlgebra.h"
#include "TACSObject.h"

/*
  Compute the transformation from the local coordinates
  to
*/
class TACSSpringTransform : public TACSObject {
 public:
  /*
    Given the local spring element reference frame Xf, compute the
    transformation from the global coordinates to the spring-aligned local axis.
  */
  virtual void computeTransform(const TacsScalar Xpts[], TacsScalar T[]) = 0;

  /*
    Compute the derivative of the transform with respect to the given
    nodal coordinate.
  */
  virtual void computeTransformXptSens(const TacsScalar Xpts[], int component,
                                       TacsScalar t[], TacsScalar tSens[]) {
    memset(tSens, 0, 9 * sizeof(TacsScalar));
    computeTransform(Xpts, t);
  }
};

class TACSSpringIdentityTransform : public TACSSpringTransform {
 public:
  TACSSpringIdentityTransform() {}

  void computeTransform(const TacsScalar Xpts[], TacsScalar T[]) {
    memset(T, 0, 9 * sizeof(TacsScalar));
    T[0] = T[4] = T[8] = 1.0;
  }
};

class TACSSpringRefAxisTransform : public TACSSpringTransform {
 public:
  TACSSpringRefAxisTransform(const TacsScalar _axis[]) {
    ref_dir[0] = _axis[0];
    ref_dir[1] = _axis[1];
    ref_dir[2] = _axis[2];

    vec3Normalize(ref_dir);
  }

  /*
    Compute the local transformation from the global axis to the local
    reference frame for the spring element.

    The first direction is along the axis of the spring:

    t[0] = (Xpts[3:6] - Xpts[0:3])/||Xpts[3:6] - Xpts[0:3]||

    The third direction is perpendicular to the reference axis and the
    direction along the length of the spring

    t[6] = (t[0] x ref_dir)/||t[0] x ref_dir||

    The second direction is perpendicular to both the t[6] and t[0]

    t[3] = t[6] x t[0]
  */
  void computeTransform(const TacsScalar Xpts[], TacsScalar t[]) {
    t[0] = Xpts[3] - Xpts[0];
    t[1] = Xpts[4] - Xpts[1];
    t[2] = Xpts[5] - Xpts[2];
    vec3Normalize(t);

    if (vec3Dot(ref_dir, &t[0]) == 1.0) {
      fprintf(stderr,
              "TACSSpringElement: Error, reference axis and spring axis are "
              "parallel\n");
    } else {
      crossProduct(&t[0], ref_dir, &t[6]);
      vec3Normalize(&t[6]);
      crossProduct(&t[6], &t[0], &t[3]);
    }
  }

  void computeTransformXptSens(const TacsScalar Xpts[], int component,
                               TacsScalar t[], TacsScalar tSens[]) {
    TacsScalar LSens;

    t[0] = Xpts[3] - Xpts[0];
    t[1] = Xpts[4] - Xpts[1];
    t[2] = Xpts[5] - Xpts[2];

    tSens[0] = tSens[1] = tSens[2] = 0.0;
    if (component < 3) {
      tSens[component] = -1.0;
    } else {
      tSens[component - 3] = 1.0;
    }

    vec3NormalizeSens(t, &LSens, tSens);

    crossProduct(&t[0], ref_dir, &t[6]);
    crossProduct(&tSens[0], ref_dir, &tSens[6]);
    TacsScalar temp;
    vec3NormalizeSens(&t[6], &temp, &tSens[6]);
    crossProductSens(&t[6], &t[0], &tSens[6], &tSens[0], &t[3], &tSens[3]);
  }

 private:
  TacsScalar ref_dir[3];
};

class TACSSpringRefFrameTransform : public TACSSpringTransform {
 public:
  TACSSpringRefFrameTransform(TacsScalar _ref_axis_i[],
                              TacsScalar _ref_axis_j[]) {
    transform[0] = _ref_axis_i[0];
    transform[1] = _ref_axis_i[1];
    transform[2] = _ref_axis_i[2];
    vec3Normalize(transform);
    crossProduct(&transform[0], _ref_axis_j, &transform[6]);
    vec3Normalize(&transform[6]);
    crossProduct(&transform[6], &transform[0], &transform[3]);
  }

  void computeTransform(const TacsScalar Xpts[], TacsScalar T[]) {
    memcpy(T, transform, 9 * sizeof(TacsScalar));
  }

 private:
  TacsScalar transform[9];
};

#endif  // TACS_SPRING_ELEMENT_TRANSFORM_H