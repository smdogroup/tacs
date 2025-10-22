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
  Scalar(const Scalar& a) { value = a.value; }
  Scalar(const TacsScalar a) { value = a; }
  TacsScalar value;
};

/*
  Active scalar type
*/
class ADScalar {
 public:
  ADScalar() {
    value = 0.0;
    valued = 0.0;
  }
  ADScalar(const TacsScalar& a) {
    value = a;
    valued = 0.0;
  }
  ADScalar(const TacsScalar& a, const TacsScalar& ad) {
    value = a;
    valued = ad;
  }
  ADScalar(const ADScalar& a) {
    value = a.value;
    valued = a.valued;
  }

  TacsScalar value;
  TacsScalar valued;
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
  Vec3(const Vec3& a) {
    for (int i = 0; i < 3; i++) {
      x[i] = a.x[i];
    }
  }

  TacsScalar x[3];
};

/*
  Active vector type
*/
class ADVec3 {
 public:
  ADVec3() {
    for (int i = 0; i < 3; i++) {
      x[i] = 0.0;
      xd[i] = 0.0;
    }
  }
  ADVec3(const TacsScalar vx, const TacsScalar vy, const TacsScalar vz) {
    x[0] = vx;
    x[1] = vy;
    x[2] = vz;
    xd[0] = xd[1] = xd[2] = 0.0;
  }
  ADVec3(const TacsScalar a[]) {
    for (int i = 0; i < 3; i++) {
      x[i] = a[i];
      xd[i] = 0.0;
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
  }
  ADVec3(const ADVec3& a) {
    for (int i = 0; i < 3; i++) {
      x[i] = a.x[i];
      xd[i] = a.xd[i];
    }
  }

  TacsScalar x[3], xd[3];
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
  Symm3x3(const Symm3x3& a) {
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
  ADSymm3x3(const ADSymm3x3& a) {
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
  Mat3x2(const Mat3x2& a) {
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
  ADMat3x2(const ADMat3x2& a) {
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
  Mat3x3(const Mat3x3& a) {
    for (int i = 0; i < 9; i++) {
      A[i] = a.A[i];
    }
  }

  TacsScalar A[9];
};

/*
  Active 3x3 matrix class
*/
class ADMat3x3 {
 public:
  ADMat3x3() {
    for (int i = 0; i < 9; i++) {
      A[i] = 0.0;
      Ad[i] = 0.0;
    }
  }
  ADMat3x3(const TacsScalar a[]) {
    for (int i = 0; i < 9; i++) {
      A[i] = a[i];
      Ad[i] = 0.0;
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
  }
  ADMat3x3(const ADMat3x3& a) {
    for (int i = 0; i < 9; i++) {
      A[i] = a.A[i];
      Ad[i] = a.Ad[i];
    }
  }

  TacsScalar A[9], Ad[9];
};

}  // namespace A2D

#endif  // A2D_OBJS_H
