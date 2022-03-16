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
  Vec3Norm( const Vec3& x, Scalar& alpha ){
    alpha.value = Vec3NormCore(x.x);
  }
};

class ADVec3Norm {
public:
  ADVec3Norm( ADVec3& x, ADScalar& alpha ) : x(x), alpha(alpha) {
    alpha.value = Vec3NormCore(x.x);
  }
  void forward(){
    if (alpha.value != 0.0){
      alpha.valued = Vec3DotCore(x.x, x.xd)/alpha.value;
    }
  }
  void reverse(){
    if (alpha.value != 0.0){
      TacsScalar ainv = alpha.valued/alpha.value;
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
  Vec3Scale( const Scalar& alpha, Vec3& x, Vec3& v ){
    Vec3ScaleCore(alpha.value, x.x, v.x);
  }
};

class ADVec3Scale {
public:
  ADVec3Scale( ADScalar& alpha, ADVec3& x, ADVec3& v ) : alpha(alpha), x(x), v(v) {
    Vec3ScaleCore(alpha.value, x.x, v.x);
  }
  void forward(){
    v.xd[0] = alpha.valued * x.x[0] + alpha.value * x.xd[0];
    v.xd[1] = alpha.valued * x.x[1] + alpha.value * x.xd[1];
    v.xd[2] = alpha.valued * x.x[2] + alpha.value * x.xd[2];
  }
  void reverse(){
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
  Vec3Axpy( const Scalar& alpha, const Vec3& x, Vconst ec3& y, Vec3& v ){
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  Vec3Axpy( const TacsScalar scale, const Scalar& alpha, const Vec3& x, const Vec3& y, Vec3& v ){
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
};

class ADVec3Axpy {
public:
  ADVec3Axpy( ADScalar& alpha, ADVec3& x, ADVec3& y, ADVec3& v ) : scale(1.0), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(alpha.value, x.x, y.x, v.x);
  }
  ADVec3Axpy( const TacsScalar scale, ADScalar& alpha, ADVec3& x, ADVec3& y, ADVec3& v ) : scale(scale), alpha(alpha), x(x), y(y), v(v) {
    Vec3AXPYCore(scale * alpha.value, x.x, y.x, v.x);
  }
  void forward(){
    v.xd[0] = scale * (alpha.valued * x.x[0] + alpha.value * x.xd[0]) + y.xd[0];
    v.xd[1] = scale * (alpha.valued * x.x[1] + alpha.value * x.xd[1]) + y.xd[1];
    v.xd[2] = scale * (alpha.valued * x.x[2] + alpha.value * x.xd[2]) + y.xd[2];
  }
  void reverse(){
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
  Vec3Dot( const Vec3& x, const Vec3& y, Scalar& alpha ){
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  Vec3Dot( const TacsScalar scale, const Vec3& x, const Vec3& y, Scalar& alpha ){
    alpha.value = scale * Vec3DotCore(x.x, y.x);
  }
};

class ADVec3Dot {
public:
  ADVec3Dot( ADVec3& x, ADVec3& y, ADScalar& alpha ) : scale(1.0), x(x), y(y), alpha(alpha) {
    alpha.value = Vec3DotCore(x.x, y.x);
  }
  ADVec3Dot( const TacsScalar scale, ADVec3& x, ADVec3& y, ADScalar& alpha ) : scale(scale), x(x), y(y), alpha(alpha) {
    alpha.value = scale * Vec3DotCore(x.x, y.x);
  }
  void forward(){
    alpha.valued = scale * (Vec3DotCore(x.x, y.xd) + Vec3DotCore(x.xd, y.x));
  }
  void reverse(){
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
  Vec3CrossProduct( const Vec3& x, const Vec3& y, Vec3& v ){
    Vec3CrossProductCore(x.x, y.x, v.x);
  }
};

class ADVec3CrossProduct {
public:
  ADVec3CrossProduct( ADVec3& x, ADVec3& y, ADVec3& v ) : x(x), y(y), v(v) {
    Vec3CrossProductCore(x.x, y.x, v.x);
  }
  void forward(){
    Vec3CrossProductCore(x.xd, y.x, v.xd);
    Vec3CrossProductAddCore(x.x, y.xd, v.xd);
  }
  void reverse(){
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
  Vec3Normalize( const Vec3& x, Vec3& y ){
    TacsScalar alpha = Vec3DotCore(x.x, x.x);
    if (alpha != 0.0){
      TacsScalar inv = 1.0/sqrt(alpha);
      y.x[0] = inv * x.x[0];
      y.x[1] = inv * x.x[1];
      y.x[2] = inv * x.x[2];
    }
  }
};

class ADVec3Normalize {
public:
  ADVec3Normalize( ADVec3& x, ADVec3& y ) : x(x), y(y) {
    alpha = Vec3DotCore(x.x, x.x);
    if (alpha != 0.0){
      inv = 1.0/sqrt(alpha);
      y.x[0] = inv * x.x[0];
      y.x[1] = inv * x.x[1];
      y.x[2] = inv * x.x[2];
    }
    else {
      inv = 0.0;
    }
  }
  void forward(){
    TacsScalar beta = Vec3DotCore(x.x, x.xd);
    TacsScalar scale = inv * inv * inv;
    y.xd[0] = (alpha * x.xd[0] - beta * x.x[0]) * scale;
    y.xd[1] = (alpha * x.xd[1] - beta * x.x[1]) * scale;
    y.xd[2] = (alpha * x.xd[2] - beta * x.x[2]) * scale;
  }
  void reverse(){
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
  Mat3x2ToVec3( const Mat3x2& A, Vec3& x, Vec3& y ){
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
  ADMat3x2ToADVec3( ADMat3x2& A, ADVec3& x, ADVec3& y ) : A(A), x(x), y(y) {
    x.x[0] = A.A[0];
    x.x[1] = A.A[2];
    x.x[2] = A.A[4];

    y.x[0] = A.A[1];
    y.x[1] = A.A[3];
    y.x[2] = A.A[5];
  }
  void forward(){
    x.xd[0] = A.Ad[0];
    x.xd[1] = A.Ad[2];
    x.xd[2] = A.Ad[4];

    y.xd[0] = A.Ad[1];
    y.xd[1] = A.Ad[3];
    y.xd[2] = A.Ad[5];
  }
  void reverse(){
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

// /*
//   Vec3SymmOuterProduct
// */
// class Vec3SymmOuterProduct {
// public:
//   Vec3SymmOuterProduct( Scalar& alpha, Vec3& x, Symm3x3& S ){
//     Vec3SymmOuterProductCore(alpha.value, x.x, S.A);
//   }
// };

// class ADVec3SymmOuterProduct {
// public:
//   ADVec3SymmOuterProduct( ADScalar& alpha, ADVec3& x, ADSymm3x3& S ){
//     Vec3SymmOuterProductCore(alpha.value, x.x, S.A);
//   }
//   void computeDeriv(){

//   }

//   ADScalar& alpha;
//   ADVec3& x;
//   ADSymm3x3& S;
// };

// /*
//   Vec3OuterProduct
// */
// class Vec3OuterProduct {
// public:
//   Vec3OuterProduct( ADScalar& alpha, Vec3& x, Vec3& y, Mat3x3& A ){
//     Vec3OuterPproductCore(alpha.value, x.x, y.x, A.A);
//   }
// };

// class ADVec3OuterProduct {
// public:
//   ADVec3OuterProduct( ADScalar& alpha, Vec3& x, Vec3& y, Mat3x3& A ) : alpha(alpha), x(x), y(y), A(A) {
//     Vec3OuterPproductCore(alpha.value, x.x, y.x, A.A);
//   }
//   void computeDeriv(){

//   }
//   ADScalar& alpha;
//   Vec3& x;
//   Vec3& y;
//   Mat3x3& A;
// };

} // namespace AD

#endif // A2D_VEC_OPS_H
