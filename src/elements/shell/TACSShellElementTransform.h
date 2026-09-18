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

  /*
    Given the forward inputs Xxi/n0 (as passed to computeTransform) and an
    upstream seed dT on the transform T, accumulate the sensitivity of T
    with respect to Xxi and n0 into dXxi/dn0 (accumulate, not overwrite).
  */
  virtual void addTransformSens(const TacsScalar Xxi[], const TacsScalar n0[],
                                const TacsScalar dT[], TacsScalar dXxi[],
                                TacsScalar dn0[]) = 0;
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

  /*
    Adjoint of computeTransform() above: given the forward inputs Xxi/n0 and
    an upstream seed dT on T, accumulate dXxi/dn0.
  */
  void addTransformSens(const TacsScalar Xxi[], const TacsScalar n0[],
                        const TacsScalar dT[], TacsScalar dXxi[],
                        TacsScalar dn0[]) {
    // Recompute the forward quantities needed for the reverse sweep
    TacsScalar n[3];
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];

    TacsScalar nnrm = sqrt(vec3Dot(n, n));
    TacsScalar ninv = 1.0 / nnrm;
    vec3Scale(ninv, n);

    TacsScalar t1_orig[3];
    t1_orig[0] = Xxi[0];
    t1_orig[1] = Xxi[2];
    t1_orig[2] = Xxi[4];

    TacsScalar d = vec3Dot(n, t1_orig);

    TacsScalar t1[3];
    t1[0] = t1_orig[0] - 3.0 * d * n[0];
    t1[1] = t1_orig[1];
    t1[2] = t1_orig[2];

    // Guard against an exactly-degenerate in-plane projection, matching the
    // same zero-guarded-normalize idiom applied in TACSShellRefAxisTransform
    // below (see that class's addTransformSens for the full NaN-propagation
    // rationale): without this guard, a t1 that lands exactly on the zero
    // vector produces t1n = 0.0*Inf = NaN, which then survives multiplying
    // an exactly-zero incoming seed (NaN*0.0 = NaN, not 0.0).
    TacsScalar t1nrm = sqrt(vec3Dot(t1, t1));
    TacsScalar t1inv = 0.0;
    if (t1nrm != 0.0) {
      t1inv = 1.0 / t1nrm;
    }
    TacsScalar t1n[3];
    t1n[0] = t1[0] * t1inv;
    t1n[1] = t1[1] * t1inv;
    t1n[2] = t1[2] * t1inv;

    // Reverse sweep: T = [t1n | t2 | n] (columns)
    TacsScalar dt1n[3] = {dT[0], dT[3], dT[6]};
    TacsScalar dt2[3] = {dT[1], dT[4], dT[7]};
    TacsScalar dn[3] = {dT[2], dT[5], dT[8]};

    // t2 = cross(n, t1n): dn += cross(t1n, dt2), dt1n += cross(dt2, n)
    TacsScalar tmp[3];
    crossProduct(t1n, dt2, tmp);
    dn[0] += tmp[0];
    dn[1] += tmp[1];
    dn[2] += tmp[2];
    crossProduct(dt2, n, tmp);
    dt1n[0] += tmp[0];
    dt1n[1] += tmp[1];
    dt1n[2] += tmp[2];

    // t1n = normalize(t1)
    TacsScalar t1copy[3] = {t1[0], t1[1], t1[2]};
    TacsScalar dt1[3] = {dt1n[0], dt1n[1], dt1n[2]};
    TacsScalar sAnrm;
    vec3NormalizeSens(t1copy, &sAnrm, dt1);

    // t1 = [t1_orig[0] - 3*d*n[0], t1_orig[1], t1_orig[2]]
    TacsScalar dt1_orig[3] = {0.0, 0.0, 0.0};
    TacsScalar dd = 0.0;
    dt1_orig[0] += dt1[0];
    dd += -3.0 * n[0] * dt1[0];
    dn[0] += -3.0 * d * dt1[0];
    dt1_orig[1] += dt1[1];
    dt1_orig[2] += dt1[2];

    // d = dot(n, t1_orig)
    vec3Axpy(dd, t1_orig, dn);
    vec3Axpy(dd, n, dt1_orig);

    // t1_orig = [Xxi[0], Xxi[2], Xxi[4]]
    dXxi[0] += dt1_orig[0];
    dXxi[2] += dt1_orig[1];
    dXxi[4] += dt1_orig[2];

    // n = normalize(n0)
    TacsScalar n0copy[3] = {n0[0], n0[1], n0[2]};
    vec3NormalizeSens(n0copy, &sAnrm, dn);

    dn0[0] += dn[0];
    dn0[1] += dn[1];
    dn0[2] += dn[2];
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

  void getRefAxis(TacsScalar _axis[]) {
    _axis[0] = axis[0];
    _axis[1] = axis[1];
    _axis[2] = axis[2];
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

    // Check if ref axis is parallel with normal
    if (abs(TacsRealPart(an)) > 1.0 - SMALL_NUM) {
      fprintf(stderr,
              "TACSShellRefAxisTransform: Error, user-provided reference axis "
              "is perpendicular to shell. "
              "Element behavior may be ill-conditioned.\n");
    }

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

  /*
    Adjoint of computeTransform() above: given the forward inputs Xxi/n0 and
    an upstream seed dT on T, accumulate dXxi/dn0. Note computeTransform()
    above never uses Xxi (the reference axis, not Xxi, supplies the in-plane
    direction), so this override contributes nothing to dXxi.
  */
  void addTransformSens(const TacsScalar Xxi[], const TacsScalar n0[],
                        const TacsScalar dT[], TacsScalar dXxi[],
                        TacsScalar dn0[]) {
    // Recompute the forward quantities needed for the reverse sweep
    TacsScalar n[3];
    n[0] = n0[0];
    n[1] = n0[1];
    n[2] = n0[2];

    TacsScalar nnrm = sqrt(vec3Dot(n, n));
    TacsScalar ninv = 1.0 / nnrm;
    vec3Scale(ninv, n);

    TacsScalar an = vec3Dot(axis, n);

    TacsScalar t1[3];
    t1[0] = axis[0] - an * n[0];
    t1[1] = axis[1] - an * n[1];
    t1[2] = axis[2] - an * n[2];

    // Guard against an exactly-degenerate in-plane projection (reference
    // axis exactly parallel to the shell normal -- the same ill-conditioned
    // case computeTransform() warns about above). Without this guard,
    // t1inv = 1.0/0.0 = +Inf and t1n = t1*t1inv degrades to 0.0*Inf = NaN
    // even though t1 is exactly the zero vector; that NaN then survives a
    // multiply by an exactly-zero incoming seed (NaN*0.0 = NaN, not 0.0 in
    // IEEE 754) a few lines below, corrupting dn/dt1n even when the caller
    // passed dT = 0 (as every geometry-only addPointQuantityXptSens branch
    // does). Mirrors the same zero-guarded-normalize idiom already used by
    // vec3NormalizeSens (TACSElementAlgebra.h) and by this class's own
    // constructor/TacsShellComputeNodeNormals.
    TacsScalar t1nrm = sqrt(vec3Dot(t1, t1));
    TacsScalar t1inv = 0.0;
    if (t1nrm != 0.0) {
      t1inv = 1.0 / t1nrm;
    }
    TacsScalar t1n[3];
    t1n[0] = t1[0] * t1inv;
    t1n[1] = t1[1] * t1inv;
    t1n[2] = t1[2] * t1inv;

    // Reverse sweep: T = [t1n | t2 | n] (columns)
    TacsScalar dt1n[3] = {dT[0], dT[3], dT[6]};
    TacsScalar dt2[3] = {dT[1], dT[4], dT[7]};
    TacsScalar dn[3] = {dT[2], dT[5], dT[8]};

    // t2 = cross(n, t1n): dn += cross(t1n, dt2), dt1n += cross(dt2, n)
    TacsScalar tmp[3];
    crossProduct(t1n, dt2, tmp);
    dn[0] += tmp[0];
    dn[1] += tmp[1];
    dn[2] += tmp[2];
    crossProduct(dt2, n, tmp);
    dt1n[0] += tmp[0];
    dt1n[1] += tmp[1];
    dt1n[2] += tmp[2];

    // t1n = normalize(t1)
    TacsScalar t1copy[3] = {t1[0], t1[1], t1[2]};
    TacsScalar dt1[3] = {dt1n[0], dt1n[1], dt1n[2]};
    TacsScalar sAnrm;
    vec3NormalizeSens(t1copy, &sAnrm, dt1);

    // t1[i] = axis[i] - an*n[i]  (axis is a constant, no sensitivity path)
    TacsScalar dan = -vec3Dot(dt1, n);
    vec3Axpy(-an, dt1, dn);

    // an = dot(axis, n)
    vec3Axpy(dan, axis, dn);

    // n = normalize(n0)
    TacsScalar n0copy[3] = {n0[0], n0[1], n0[2]};
    vec3NormalizeSens(n0copy, &sAnrm, dn);

    dn0[0] += dn[0];
    dn0[1] += dn[1];
    dn0[2] += dn[2];
  }

 private:
  TacsScalar axis[3];
  /* Tolerance for colinearity test in between shell normal and ref axis */
  const double SMALL_NUM = 1e-8;
};

#endif  // TACS_SHELL_ELEMENT_TRANSFORM_H