#ifndef A2D_OBJS_H
#define A2D_OBJS_H

#include "TACSObject.h"

namespace A2D {

/*
  Scalar type
*/
class Scalar {
 public:
  Scalar() {}
  Scalar(const Scalar &a) { value = a.value; }
  Scalar(const TacsScalar a) { value = a; }
  TacsScalar value;
};

/*
  Active scalar type

  Second-order extension (feature-beam-element-methods, SPEC.md sec 1.2.2):
  valuep/valueh are additive fields for the forward-over-reverse
  hforward()/hreverse() scheme -- valuep is a caller-chosen forward "seed"
  direction, valueh is the resulting (projected) second-order reverse
  accumulator. Extended in place rather than via a new parallel type: every
  existing ADScalar construction call site in this codebase (grepped across
  TACSBeamElement.h/TACSBeamUtilities.h) uses only the (), (TacsScalar), or
  (TacsScalar, TacsScalar) constructors, none of which change shape here --
  the new fields are simply zero-initialized alongside value/valued, so no
  first-order call site is affected.
*/
class ADScalar {
 public:
  ADScalar() {
    value = 0.0;
    valued = 0.0;
    valuep = 0.0;
    valueh = 0.0;
  }
  ADScalar(const TacsScalar &a) {
    value = a;
    valued = 0.0;
    valuep = 0.0;
    valueh = 0.0;
  }
  ADScalar(const TacsScalar &a, const TacsScalar &ad) {
    value = a;
    valued = ad;
    valuep = 0.0;
    valueh = 0.0;
  }
  ADScalar(const ADScalar &a) {
    value = a.value;
    valued = a.valued;
    valuep = a.valuep;
    valueh = a.valueh;
  }

  TacsScalar value;
  TacsScalar valued;
  TacsScalar valuep;
  TacsScalar valueh;
};

/*
  Passive vector type
*/
class Vec3 {
 public:
  Vec3() {
    for (int i = 0; i < 3; i++) {
      x[i] = 0.0;
    }
  }
  Vec3(const TacsScalar vx, const TacsScalar vy, const TacsScalar vz) {
    x[0] = vx;
    x[1] = vy;
    x[2] = vz;
  }
  Vec3(const TacsScalar a[]) {
    for (int i = 0; i < 3; i++) {
      x[i] = a[i];
    }
  }
  Vec3(const Vec3 &a) {
    for (int i = 0; i < 3; i++) {
      x[i] = a.x[i];
    }
  }

  TacsScalar x[3];
};

/*
  Active vector type

  Second-order extension (feature-beam-element-methods, SPEC.md sec 1.2.2):
  xp[3]/xh[3] are additive fields for the forward-over-reverse
  hforward()/hreverse() scheme -- xp is a caller-chosen forward "seed"
  direction, xh is the resulting (projected) second-order reverse
  accumulator, mirroring the existing x/xd (primal/first-order-adjoint)
  convention. Extended in place rather than via a new parallel type: every
  existing ADVec3 construction call site in this codebase (grepped across
  TACSBeamElement.h/TACSBeamUtilities.h) uses only the default constructor,
  so the new fields (zero-initialized everywhere, including inside the
  existing 1-2 array-argument constructors) cannot affect any first-order
  call site.
*/
class ADVec3 {
 public:
  ADVec3() {
    for (int i = 0; i < 3; i++) {
      x[i] = 0.0;
      xd[i] = 0.0;
      xp[i] = 0.0;
      xh[i] = 0.0;
    }
  }
  ADVec3(const TacsScalar vx, const TacsScalar vy, const TacsScalar vz) {
    x[0] = vx;
    x[1] = vy;
    x[2] = vz;
    xd[0] = xd[1] = xd[2] = 0.0;
    xp[0] = xp[1] = xp[2] = 0.0;
    xh[0] = xh[1] = xh[2] = 0.0;
  }
  ADVec3(const TacsScalar a[]) {
    for (int i = 0; i < 3; i++) {
      x[i] = a[i];
      xd[i] = 0.0;
      xp[i] = 0.0;
      xh[i] = 0.0;
    }
  }
  ADVec3(const TacsScalar a[], const TacsScalar ad[]) {
    if (a) {
      for (int i = 0; i < 3; i++) {
        x[i] = a[i];
      }
    } else {
      for (int i = 0; i < 3; i++) {
        x[i] = 0.0;
      }
    }
    if (ad) {
      for (int i = 0; i < 3; i++) {
        xd[i] = ad[i];
      }
    } else {
      for (int i = 0; i < 3; i++) {
        xd[i] = 0.0;
      }
    }
    for (int i = 0; i < 3; i++) {
      xp[i] = 0.0;
      xh[i] = 0.0;
    }
  }
  ADVec3(const ADVec3 &a) {
    for (int i = 0; i < 3; i++) {
      x[i] = a.x[i];
      xd[i] = a.xd[i];
      xp[i] = a.xp[i];
      xh[i] = a.xh[i];
    }
  }

  TacsScalar x[3], xd[3], xp[3], xh[3];
};

/*
  Passive symmetric 3x3 matrix
*/
class Symm3x3 {
 public:
  Symm3x3() {
    for (int i = 0; i < 6; i++) {
      A[i] = 0.0;
    }
  }
  Symm3x3(const TacsScalar a[]) {
    for (int i = 0; i < 6; i++) {
      A[i] = a[i];
    }
  }
  Symm3x3(const Symm3x3 &a) {
    for (int i = 0; i < 6; i++) {
      A[i] = a.A[i];
    }
  }

  TacsScalar A[6];
};

/*
  Active symmetric 3x3 matrix class
*/
class ADSymm3x3 {
 public:
  ADSymm3x3() {
    for (int i = 0; i < 6; i++) {
      A[i] = 0.0;
      Ad[i] = 0.0;
    }
  }
  ADSymm3x3(const TacsScalar a[]) {
    for (int i = 0; i < 6; i++) {
      A[i] = a[i];
      Ad[i] = 0.0;
    }
  }
  ADSymm3x3(const TacsScalar a[], const TacsScalar ad[]) {
    if (a) {
      for (int i = 0; i < 6; i++) {
        A[i] = a[i];
      }
    } else {
      for (int i = 0; i < 6; i++) {
        A[i] = 0.0;
      }
    }
    if (ad) {
      for (int i = 0; i < 6; i++) {
        Ad[i] = ad[i];
      }
    } else {
      for (int i = 0; i < 6; i++) {
        Ad[i] = 0.0;
      }
    }
  }
  ADSymm3x3(const ADSymm3x3 &a) {
    for (int i = 0; i < 6; i++) {
      A[i] = a.A[i];
      Ad[i] = a.Ad[i];
    }
  }

  TacsScalar A[6], Ad[6];
};

/*
  Passive 3x2 matrix class
*/
class Mat3x2 {
 public:
  Mat3x2() {
    for (int i = 0; i < 6; i++) {
      A[i] = 0.0;
    }
  }
  Mat3x2(const TacsScalar a[]) {
    for (int i = 0; i < 6; i++) {
      A[i] = a[i];
    }
  }
  Mat3x2(const Mat3x2 &a) {
    for (int i = 0; i < 6; i++) {
      A[i] = a.A[i];
    }
  }

  TacsScalar A[6];
};

/*
  Active 3x2 matrix class
*/
class ADMat3x2 {
 public:
  ADMat3x2() {
    for (int i = 0; i < 6; i++) {
      A[i] = 0.0;
      Ad[i] = 0.0;
    }
  }
  ADMat3x2(const TacsScalar a[]) {
    for (int i = 0; i < 6; i++) {
      A[i] = a[i];
      Ad[i] = 0.0;
    }
  }
  ADMat3x2(const TacsScalar a[], const TacsScalar ad[]) {
    if (a) {
      for (int i = 0; i < 6; i++) {
        A[i] = a[i];
      }
    } else {
      for (int i = 0; i < 6; i++) {
        A[i] = 0.0;
      }
    }
    if (ad) {
      for (int i = 0; i < 6; i++) {
        Ad[i] = ad[i];
      }
    } else {
      for (int i = 0; i < 6; i++) {
        Ad[i] = 0.0;
      }
    }
  }
  ADMat3x2(const ADMat3x2 &a) {
    for (int i = 0; i < 6; i++) {
      A[i] = a.A[i];
      Ad[i] = a.Ad[i];
    }
  }

  TacsScalar A[6];
  TacsScalar Ad[6];
};

/*
  Passive 3x3 matrix class
*/
class Mat3x3 {
 public:
  Mat3x3() {
    for (int i = 0; i < 9; i++) {
      A[i] = 0.0;
    }
  }
  Mat3x3(const TacsScalar a[]) {
    for (int i = 0; i < 9; i++) {
      A[i] = a[i];
    }
  }
  Mat3x3(const Mat3x3 &a) {
    for (int i = 0; i < 9; i++) {
      A[i] = a.A[i];
    }
  }

  TacsScalar A[9];
};

/*
  Active 3x3 matrix class

  Second-order extension (feature-beam-element-methods, SPEC.md sec 1.2.2):
  Ap[9]/Ah[9] are additive fields for the forward-over-reverse
  hforward()/hreverse() scheme -- Ap is a caller-chosen forward "seed"
  direction, Ah is the resulting (projected) second-order reverse
  accumulator, mirroring the existing A/Ad (primal/first-order-adjoint)
  convention. Extended in place rather than via a new parallel type: every
  existing ADMat3x3 construction call site in this codebase (grepped across
  TACSBeamElement.h/TACSBeamUtilities.h) uses only the default constructor,
  so the new fields (zero-initialized everywhere, including inside the
  existing 1-2 array-argument constructors) cannot affect any first-order
  call site.
*/
class ADMat3x3 {
 public:
  ADMat3x3() {
    for (int i = 0; i < 9; i++) {
      A[i] = 0.0;
      Ad[i] = 0.0;
      Ap[i] = 0.0;
      Ah[i] = 0.0;
    }
  }
  ADMat3x3(const TacsScalar a[]) {
    for (int i = 0; i < 9; i++) {
      A[i] = a[i];
      Ad[i] = 0.0;
      Ap[i] = 0.0;
      Ah[i] = 0.0;
    }
  }
  ADMat3x3(const TacsScalar a[], const TacsScalar ad[]) {
    if (a) {
      for (int i = 0; i < 9; i++) {
        A[i] = a[i];
      }
    } else {
      for (int i = 0; i < 9; i++) {
        A[i] = 0.0;
      }
    }
    if (ad) {
      for (int i = 0; i < 9; i++) {
        Ad[i] = ad[i];
      }
    } else {
      for (int i = 0; i < 9; i++) {
        Ad[i] = 0.0;
      }
    }
    for (int i = 0; i < 9; i++) {
      Ap[i] = 0.0;
      Ah[i] = 0.0;
    }
  }
  ADMat3x3(const ADMat3x3 &a) {
    for (int i = 0; i < 9; i++) {
      A[i] = a.A[i];
      Ad[i] = a.Ad[i];
      Ap[i] = a.Ap[i];
      Ah[i] = a.Ah[i];
    }
  }

  TacsScalar A[9], Ad[9], Ap[9], Ah[9];
};

}  // namespace A2D

#endif  // A2D_OBJS_H
