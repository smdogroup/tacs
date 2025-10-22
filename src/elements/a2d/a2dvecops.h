#ifndef A2D_VEC_OPS_H
#define A2D_VEC_OPS_H

#include "a2dobjs.h"
#include "a2dveccore.h"

namespace A2D {

/*
  Vec3Norm
*/
class Vec3Norm {
 public:
  Vec3Norm(const Vec3& x, Scalar& alpha) { alpha.value = Vec3NormCore(x.x); }
};

class ADVec3Norm {
 public:
  ADVec3Norm(ADVec3& x, ADScalar& alpha) : x(x), alpha(alpha) {
    alpha.value = Vec3NormCore(x.x);
  }
  void forward() {
    if (alpha.value != 0.0) {
      alpha.valued = Vec3DotCore(x.x, x.xd) / alpha.value;
    }
  }
  void reverse() {
    if (alpha.value != 0.0) {
      TacsScalar ainv = alpha.valued / alpha.value;
      x.xd[0] += ainv * x.x[0];
      x.xd[1] += ainv * x.x[1];
      x.xd[2] += ainv * x.x[2];
    }
  }

  ADVec3& x;
  ADScalar& alpha;
};

/*
  Vec3Scale
*/
class Vec3Scale {
 public:
  Vec3Scale(const Scalar& alpha, Vec3& x, Vec3& v) {
    Vec3ScaleCore(alpha.value, x.x, v.x);
  }
};

class ADVec3Scale {
 public:
  ADVec3Scale(ADScalar& alpha, ADVec3& x, ADVec3& v)
      : alpha(alpha), x(x), v(v) {
    Vec3ScaleCore(alpha.value, x.x, v.x);
  }
  void forward() {
    v.xd[0] = alpha.valued * x.x[0] + alpha.value * x.xd[0];
    v.xd[1] = alpha.valued * x.x[1] + alpha.value * x.xd[1];
    v.xd[2] = alpha.valued * x.x[2] + alpha.value * x.xd[2];
  }
  void reverse() {
    alpha.valued += Vec3DotCore(v.xd, x.x);
    x.xd[0] += alpha.value * v.xd[0];
    x.xd[1] += alpha.value * v.xd[1];
    x.xd[2] += alpha.value * v.xd[2];
  }

  ADScalar& alpha;
  ADVec3& x;
  ADVec3& v;
};

/*
  Vec3Axpy
*/
class Vec3Axpy {
 public:
  Vec3Axpy(const Scalar& alpha, const Vec3& x, const Vec3& y, Vec3& v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  Vec3Axpy(const TacsScalar scale, const Scalar& alpha, const Vec3& x,
           const Vec3& y, Vec3& v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
};

class Vec3VecADScalarAxpy {
 public:
  Vec3VecADScalarAxpy(ADScalar& alpha, const Vec3& x, const Vec3& y, ADVec3& v)
      : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  Vec3VecADScalarAxpy(const TacsScalar scale, ADScalar& alpha, const Vec3& x,
                      const Vec3& y, ADVec3& v)
      : scale(scale), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
  void forward() {
    v.xd[0] = scale * alpha.valued * x.x[0];
    v.xd[1] = scale * alpha.valued * x.x[1];
    v.xd[2] = scale * alpha.valued * x.x[2];
  }
  void reverse() { alpha.valued += scale * Vec3DotCore(x.x, v.xd); }

  const TacsScalar scale;
  ADScalar& alpha;
  const Vec3& x;
  const Vec3& y;
  ADVec3& v;
};

class ADVec3VecADScalarAxpy {
 public:
  ADVec3VecADScalarAxpy(ADScalar& alpha, ADVec3& x, const Vec3& y, ADVec3& v)
      : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  ADVec3VecADScalarAxpy(const TacsScalar scale, ADScalar& alpha, ADVec3& x,
                        const Vec3& y, ADVec3& v)
      : scale(scale), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
  void forward() {
    v.xd[0] = scale * (alpha.valued * x.x[0] + alpha.value * x.xd[0]);
    v.xd[1] = scale * (alpha.valued * x.x[1] + alpha.value * x.xd[1]);
    v.xd[2] = scale * (alpha.valued * x.x[2] + alpha.value * x.xd[2]);
  }
  void reverse() {
    alpha.valued += scale * Vec3DotCore(x.x, v.xd);
    x.xd[0] += scale * alpha.value * v.xd[0];
    x.xd[1] += scale * alpha.value * v.xd[1];
    x.xd[2] += scale * alpha.value * v.xd[2];
  }

  const TacsScalar scale;
  ADScalar& alpha;
  ADVec3& x;
  const Vec3& y;
  ADVec3& v;
};

class ADVec3ADVecScalarAxpy {
 public:
  ADVec3ADVecScalarAxpy(const Scalar& alpha, ADVec3& x, ADVec3& y, ADVec3& v)
      : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  ADVec3ADVecScalarAxpy(const TacsScalar scale, const Scalar& alpha, ADVec3& x,
                        ADVec3& y, ADVec3& v)
      : scale(scale), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
  void forward() {
    v.xd[0] = scale * (alpha.value * x.xd[0]) + y.xd[0];
    v.xd[1] = scale * (alpha.value * x.xd[1]) + y.xd[1];
    v.xd[2] = scale * (alpha.value * x.xd[2]) + y.xd[2];
  }
  void reverse() {
    x.xd[0] += scale * alpha.value * v.xd[0];
    x.xd[1] += scale * alpha.value * v.xd[1];
    x.xd[2] += scale * alpha.value * v.xd[2];
    y.xd[0] += v.xd[0];
    y.xd[1] += v.xd[1];
    y.xd[2] += v.xd[2];
  }

  const TacsScalar scale;
  const Scalar& alpha;
  ADVec3& x;
  ADVec3& y;
  ADVec3& v;
};

class ADVec3Axpy {
 public:
  ADVec3Axpy(ADScalar& alpha, ADVec3& x, ADVec3& y, ADVec3& v)
      : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  ADVec3Axpy(const TacsScalar scale, ADScalar& alpha, ADVec3& x, ADVec3& y,
             ADVec3& v)
      : scale(scale), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
  void forward() {
    v.xd[0] = scale * (alpha.valued * x.x[0] + alpha.value * x.xd[0]) + y.xd[0];
    v.xd[1] = scale * (alpha.valued * x.x[1] + alpha.value * x.xd[1]) + y.xd[1];
    v.xd[2] = scale * (alpha.valued * x.x[2] + alpha.value * x.xd[2]) + y.xd[2];
  }
  void reverse() {
    alpha.valued += scale * Vec3DotCore(x.x, v.xd);
    x.xd[0] += scale * alpha.value * v.xd[0];
    x.xd[1] += scale * alpha.value * v.xd[1];
    x.xd[2] += scale * alpha.value * v.xd[2];
    y.xd[0] += v.xd[0];
    y.xd[1] += v.xd[1];
    y.xd[2] += v.xd[2];
  }

  const TacsScalar scale;
  ADScalar& alpha;
  ADVec3& x;
  ADVec3& y;
  ADVec3& v;
};

/*
  Vec3Dot
*/
class Vec3Dot {
 public:
  Vec3Dot(const Vec3& x, const Vec3& y, Scalar& alpha) {
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  Vec3Dot(const TacsScalar scale, const Vec3& x, const Vec3& y, Scalar& alpha) {
    alpha.value = scale * Vec3DotCore(x.x, y.x);
  }
};

class Vec3ADVecDot {
 public:
  Vec3ADVecDot(const Vec3& x, ADVec3& y, ADScalar& alpha)
      : scale(1.0), x(x), y(y), alpha(alpha) {
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  Vec3ADVecDot(const TacsScalar scale, const Vec3& x, ADVec3& y,
               ADScalar& alpha)
      : scale(scale), x(x), y(y), alpha(alpha) {
    alpha.value = scale * Vec3DotCore(x.x, y.x);
  }
  void forward() { alpha.valued = scale * Vec3DotCore(x.x, y.xd); }
  void reverse() {
    TacsScalar s = scale * alpha.valued;
    y.xd[0] += s * x.x[0];
    y.xd[1] += s * x.x[1];
    y.xd[2] += s * x.x[2];
  }

  const TacsScalar scale;
  const Vec3& x;
  ADVec3& y;
  ADScalar& alpha;
};

class ADVec3Dot {
 public:
  ADVec3Dot(ADVec3& x, ADVec3& y, ADScalar& alpha)
      : scale(1.0), x(x), y(y), alpha(alpha) {
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  ADVec3Dot(const TacsScalar scale, ADVec3& x, ADVec3& y, ADScalar& alpha)
      : scale(scale), x(x), y(y), alpha(alpha) {
    alpha.value = scale * Vec3DotCore(x.x, y.x);
  }
  void forward() {
    alpha.valued = scale * (Vec3DotCore(x.x, y.xd) + Vec3DotCore(x.xd, y.x));
  }
  void reverse() {
    TacsScalar s = scale * alpha.valued;
    x.xd[0] += s * y.x[0];
    x.xd[1] += s * y.x[1];
    x.xd[2] += s * y.x[2];

    y.xd[0] += s * x.x[0];
    y.xd[1] += s * x.x[1];
    y.xd[2] += s * x.x[2];
  }

  const TacsScalar scale;
  ADVec3& x;
  ADVec3& y;
  ADScalar& alpha;
};

/*
  Vec3CrossProduct
*/
class Vec3CrossProduct {
 public:
  Vec3CrossProduct(const Vec3& x, const Vec3& y, Vec3& v) {
    Vec3CrossProductCore(x.x, y.x, v.x);
  }
};

class ADVec3CrossProduct {
 public:
  ADVec3CrossProduct(ADVec3& x, ADVec3& y, ADVec3& v) : x(x), y(y), v(v) {
    Vec3CrossProductCore(x.x, y.x, v.x);
  }
  void forward() {
    Vec3CrossProductCore(x.xd, y.x, v.xd);
    Vec3CrossProductAddCore(x.x, y.xd, v.xd);
  }
  void reverse() {
    Vec3CrossProductAddCore(y.x, v.xd, x.xd);
    Vec3CrossProductAddCore(v.xd, x.x, y.xd);
  }

  ADVec3& x;
  ADVec3& y;
  ADVec3& v;
};

/*
  Vec3 normalize
*/
class Vec3Normalize {
 public:
  Vec3Normalize(const Vec3& x, Vec3& y) {
    TacsScalar alpha = Vec3DotCore(x.x, x.x);
    if (alpha != 0.0) {
      TacsScalar inv = 1.0 / sqrt(alpha);
      y.x[0] = inv * x.x[0];
      y.x[1] = inv * x.x[1];
      y.x[2] = inv * x.x[2];
    }
  }
};

class ADVec3Normalize {
 public:
  ADVec3Normalize(ADVec3& x, ADVec3& y) : x(x), y(y) {
    alpha = Vec3DotCore(x.x, x.x);
    if (alpha != 0.0) {
      inv = 1.0 / sqrt(alpha);
      y.x[0] = inv * x.x[0];
      y.x[1] = inv * x.x[1];
      y.x[2] = inv * x.x[2];
    } else {
      inv = 0.0;
    }
  }
  void forward() {
    TacsScalar beta = Vec3DotCore(x.x, x.xd);
    TacsScalar scale = inv * inv * inv;
    y.xd[0] = (alpha * x.xd[0] - beta * x.x[0]) * scale;
    y.xd[1] = (alpha * x.xd[1] - beta * x.x[1]) * scale;
    y.xd[2] = (alpha * x.xd[2] - beta * x.x[2]) * scale;
  }
  void reverse() {
    TacsScalar beta = Vec3DotCore(x.x, y.xd);
    TacsScalar scale = inv * inv * inv;
    x.xd[0] += (alpha * y.xd[0] - beta * x.x[0]) * scale;
    x.xd[1] += (alpha * y.xd[1] - beta * x.x[1]) * scale;
    x.xd[2] += (alpha * y.xd[2] - beta * x.x[2]) * scale;
  }

  TacsScalar alpha, inv;
  ADVec3& x;
  ADVec3& y;
};

class Mat3x2ToVec3 {
 public:
  Mat3x2ToVec3(const Mat3x2& A, Vec3& x, Vec3& y) {
    x.x[0] = A.A[0];
    x.x[1] = A.A[2];
    x.x[2] = A.A[4];

    y.x[0] = A.A[1];
    y.x[1] = A.A[3];
    y.x[2] = A.A[5];
  }
};

class ADMat3x2ToADVec3 {
 public:
  ADMat3x2ToADVec3(ADMat3x2& A, ADVec3& x, ADVec3& y) : A(A), x(x), y(y) {
    x.x[0] = A.A[0];
    x.x[1] = A.A[2];
    x.x[2] = A.A[4];

    y.x[0] = A.A[1];
    y.x[1] = A.A[3];
    y.x[2] = A.A[5];
  }
  void forward() {
    x.xd[0] = A.Ad[0];
    x.xd[1] = A.Ad[2];
    x.xd[2] = A.Ad[4];

    y.xd[0] = A.Ad[1];
    y.xd[1] = A.Ad[3];
    y.xd[2] = A.Ad[5];
  }
  void reverse() {
    A.Ad[0] += x.xd[0];
    A.Ad[2] += x.xd[1];
    A.Ad[4] += x.xd[2];

    A.Ad[1] += y.xd[0];
    A.Ad[3] += y.xd[1];
    A.Ad[5] += y.xd[2];
  }

  ADMat3x2& A;
  ADVec3& x;
  ADVec3& y;
};

/*
  Matrix-vector product
*/
class Mat3x3VecMult {
 public:
  Mat3x3VecMult(const Mat3x3& A, const Vec3& x, Vec3& y) {
    Mat3x3VecMultCore(A.A, x.x, y.x);
  }
};

class ADMat3x3VecMult {
 public:
  ADMat3x3VecMult(ADMat3x3& A, const Vec3& x, ADVec3& y) : A(A), x(x), y(y) {
    Mat3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() { Mat3x3VecMultCore(A.Ad, x.x, y.xd); }
  void reverse() { Vec3OuterProductAddCore(y.xd, x.x, A.Ad); }

  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class Mat3x3ADVecMult {
 public:
  Mat3x3ADVecMult(const Mat3x3& A, ADVec3& x, ADVec3& y) : A(A), x(x), y(y) {
    Mat3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() { Mat3x3VecMultAddCore(A.A, x.xd, y.xd); }
  void reverse() { MatTrans3x3VecMultAddCore(A.A, y.xd, x.xd); }

  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMat3x3ADVecMult {
 public:
  ADMat3x3ADVecMult(ADMat3x3& A, ADVec3& x, ADVec3& y) : A(A), x(x), y(y) {
    Mat3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() {
    Mat3x3VecMultCore(A.Ad, x.x, y.xd);
    Mat3x3VecMultAddCore(A.A, x.xd, y.xd);
  }
  void reverse() {
    MatTrans3x3VecMultAddCore(A.A, y.xd, x.xd);
    Vec3OuterProductAddCore(y.xd, x.x, A.Ad);
  }

  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

/*
  Transpose matrix-vector product
*/
class MatTrans3x3VecMult {
 public:
  MatTrans3x3VecMult(const Mat3x3& A, const Vec3& x, Vec3& y) {
    MatTrans3x3VecMultCore(A.A, x.x, y.x);
  }
};

class ADMatTrans3x3VecMult {
 public:
  ADMatTrans3x3VecMult(ADMat3x3& A, const Vec3& x, ADVec3& y)
      : A(A), x(x), y(y) {
    MatTrans3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() { MatTrans3x3VecMultCore(A.Ad, x.x, y.xd); }
  void reverse() { Vec3OuterProductAddCore(x.x, y.xd, A.Ad); }

  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class MatTrans3x3ADVecMult {
 public:
  MatTrans3x3ADVecMult(const Mat3x3& A, ADVec3& x, ADVec3& y)
      : A(A), x(x), y(y) {
    MatTrans3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() { MatTrans3x3VecMultAddCore(A.A, x.xd, y.xd); }
  void reverse() { Mat3x3VecMultAddCore(A.A, y.xd, x.xd); }

  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMatTrans3x3ADVecMult {
 public:
  ADMatTrans3x3ADVecMult(ADMat3x3& A, ADVec3& x, ADVec3& y) : A(A), x(x), y(y) {
    MatTrans3x3VecMultCore(A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultCore(A.Ad, x.x, y.xd);
    MatTrans3x3VecMultAddCore(A.A, x.xd, y.xd);
  }
  void reverse() {
    Mat3x3VecMultAddCore(A.A, y.xd, x.xd);
    Vec3OuterProductAddCore(x.x, y.xd, A.Ad);
  }

  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

/*
  Matrix-vector product
*/
class Mat3x3VecMultScale {
 public:
  Mat3x3VecMultScale(const Scalar& scale, const Mat3x3& A, const Vec3& x,
                     Vec3& y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
};

class ADMat3x3VecMultScale {
 public:
  ADMat3x3VecMultScale(const Scalar& scale, ADMat3x3& A, const Vec3& x,
                       ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() { Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd); }
  void reverse() { Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad); }

  const Scalar& scale;
  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class Mat3x3ADVecMultScale {
 public:
  Mat3x3ADVecMultScale(const Scalar& scale, const Mat3x3& A, ADVec3& x,
                       ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() { Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd); }
  void reverse() {
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
  }

  const Scalar& scale;
  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMat3x3ADVecMultScale {
 public:
  ADMat3x3ADVecMultScale(const Scalar& scale, ADMat3x3& A, ADVec3& x, ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
  }
  void reverse() {
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
  }

  const Scalar& scale;
  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class Mat3x3VecMultADScale {
 public:
  Mat3x3VecMultADScale(ADScalar& scale, const Mat3x3& A, const Vec3& x,
                       ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() { Mat3x3VecMultScaleCore(scale.valued, A.A, x.x, y.xd); }
  void reverse() { scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x); }

  ADScalar& scale;
  const Mat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class ADMat3x3VecMultADScale {
 public:
  ADMat3x3VecMultADScale(ADScalar& scale, ADMat3x3& A, const Vec3& x, ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
    scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
  }

  ADScalar& scale;
  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class Mat3x3ADVecMultADScale {
 public:
  Mat3x3ADVecMultADScale(ADScalar& scale, const Mat3x3& A, ADVec3& x, ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.xd, y.xd);
    Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
  }

  ADScalar& scale;
  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMat3x3ADVecMultADScale {
 public:
  ADMat3x3ADVecMultADScale(ADScalar& scale, ADMat3x3& A, ADVec3& x, ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    Mat3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    Mat3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    Mat3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    Mat3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    Vec3OuterProductAddScaleCore(scale.value, y.xd, x.x, A.Ad);
    scale.valued += Mat3x3InnerProductCore(A.A, y.xd, x.x);
  }

  ADScalar& scale;
  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

/*
  Transpose matrix-vector product
*/
class MatTrans3x3VecMultScale {
 public:
  MatTrans3x3VecMultScale(const Scalar& scale, const Mat3x3& A, const Vec3& x,
                          Vec3& y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
};

class ADMatTrans3x3VecMultScale {
 public:
  ADMatTrans3x3VecMultScale(const Scalar& scale, ADMat3x3& A, const Vec3& x,
                            ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() { MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd); }
  void reverse() { Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad); }

  const Scalar& scale;
  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class MatTrans3x3ADVecMultScale {
 public:
  MatTrans3x3ADVecMultScale(const Scalar& scale, const Mat3x3& A, ADVec3& x,
                            ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
  }
  void reverse() { Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd); }

  const Scalar& scale;
  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMatTrans3x3ADVecMultScale {
 public:
  ADMatTrans3x3ADVecMultScale(const Scalar& scale, ADMat3x3& A, ADVec3& x,
                              ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
  }
  void reverse() {
    Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
  }

  const Scalar& scale;
  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class MatTrans3x3VecMultADScale {
 public:
  MatTrans3x3VecMultADScale(ADScalar& scale, const Mat3x3& A, const Vec3& x,
                            ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() { MatTrans3x3VecMultScaleCore(scale.valued, A.A, x.x, y.xd); }
  void reverse() { scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd); }

  ADScalar& scale;
  const Mat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class ADMatTrans3x3VecMultADScale {
 public:
  ADMatTrans3x3VecMultADScale(ADScalar& scale, ADMat3x3& A, const Vec3& x,
                              ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
    scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
  }

  ADScalar& scale;
  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
};

class MatTrans3x3ADVecMultADScale {
 public:
  MatTrans3x3ADVecMultADScale(ADScalar& scale, const Mat3x3& A, ADVec3& x,
                              ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.xd, y.xd);
    MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
  }

  ADScalar& scale;
  const Mat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

class ADMatTrans3x3ADVecMultADScale {
 public:
  ADMatTrans3x3ADVecMultADScale(ADScalar& scale, ADMat3x3& A, ADVec3& x,
                                ADVec3& y)
      : scale(scale), A(A), x(x), y(y) {
    MatTrans3x3VecMultScaleCore(scale.value, A.A, x.x, y.x);
  }
  void forward() {
    MatTrans3x3VecMultScaleCore(scale.value, A.Ad, x.x, y.xd);
    MatTrans3x3VecMultAddScaleCore(scale.value, A.A, x.xd, y.xd);
    MatTrans3x3VecMultAddScaleCore(scale.valued, A.A, x.x, y.xd);
  }
  void reverse() {
    Mat3x3VecMultAddScaleCore(scale.value, A.A, y.xd, x.xd);
    Vec3OuterProductAddScaleCore(scale.value, x.x, y.xd, A.Ad);
    scale.valued += Mat3x3InnerProductCore(A.A, x.x, y.xd);
  }

  ADScalar& scale;
  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
};

/*
  Inner products alpha = x^{T} * A * y
*/
class Mat3x3VecVecInnerProduct {
 public:
  Mat3x3VecVecInnerProduct(const Mat3x3& A, const Vec3& x, const Vec3& y,
                           Scalar& alpha) {
    alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
  }
};

class ADMat3x3VecVecInnerProduct {
 public:
  ADMat3x3VecVecInnerProduct(ADMat3x3& A, const Vec3& x, const Vec3& y,
                             ADScalar& alpha)
      : A(A), x(x), y(y), alpha(alpha) {
    alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
  }
  void forward() { alpha.valued = Mat3x3InnerProductCore(A.Ad, x.x, y.x); }
  void reverse() { Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad); }

  ADMat3x3& A;
  const Vec3& x;
  const Vec3& y;
  ADScalar& alpha;
};

class ADMat3x3VecADVecInnerProduct {
 public:
  ADMat3x3VecADVecInnerProduct(ADMat3x3& A, const Vec3& x, ADVec3& y,
                               ADScalar& alpha)
      : A(A), x(x), y(y), alpha(alpha) {
    alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
  }
  void forward() {
    alpha.valued = Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
                   Mat3x3InnerProductCore(A.A, x.x, y.xd);
  }
  void reverse() {
    MatTrans3x3VecMultAddScaleCore(alpha.valued, A.A, x.x, y.xd);
    Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
  }

  ADMat3x3& A;
  const Vec3& x;
  ADVec3& y;
  ADScalar& alpha;
};

class ADMat3x3ADVecVecInnerProduct {
 public:
  ADMat3x3ADVecVecInnerProduct(ADMat3x3& A, ADVec3& x, const Vec3& y,
                               ADScalar& alpha)
      : A(A), x(x), y(y), alpha(alpha) {
    alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
  }
  void forward() {
    alpha.valued = Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
                   Mat3x3InnerProductCore(A.A, x.xd, y.x);
  }
  void reverse() {
    Mat3x3VecMultAddScaleCore(alpha.valued, A.A, y.x, x.xd);
    Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
  }

  ADMat3x3& A;
  ADVec3& x;
  const Vec3& y;
  ADScalar& alpha;
};

class ADMat3x3ADVecADVecInnerProduct {
 public:
  ADMat3x3ADVecADVecInnerProduct(ADMat3x3& A, ADVec3& x, ADVec3& y,
                                 ADScalar& alpha)
      : A(A), x(x), y(y), alpha(alpha) {
    alpha.value = Mat3x3InnerProductCore(A.A, x.x, y.x);
  }
  void forward() {
    alpha.valued = Mat3x3InnerProductCore(A.Ad, x.x, y.x) +
                   Mat3x3InnerProductCore(A.A, x.xd, y.x) +
                   Mat3x3InnerProductCore(A.A, x.x, y.xd);
  }
  void reverse() {
    Mat3x3VecMultAddScaleCore(alpha.valued, A.A, y.x, x.xd);
    MatTrans3x3VecMultAddScaleCore(alpha.valued, A.A, x.x, y.xd);
    Vec3OuterProductAddScaleCore(alpha.valued, x.x, y.x, A.Ad);
  }

  ADMat3x3& A;
  ADVec3& x;
  ADVec3& y;
  ADScalar& alpha;
};

}  // namespace A2D

#endif  // A2D_VEC_OPS_H
