#ifndef TACS_A2D_UTILITIES_H
#define TACS_A2D_UTILITIES_H

#include "a2dcore.h"

/*
  TACS-side extensions to the A2D library (extern/a2d).

  These operations follow the upstream A2D expression conventions so they
  compose with A2D::MakeStack. They live in the A2D namespace so stacks read
  uniformly at the point of use.
*/
namespace A2D {

/*
  Assemble a matrix from three column vectors: C = [x | y | z]
*/
template <typename T>
A2D_FUNCTION void MatFromThreeVec(const Vec<T, 3> &x, const Vec<T, 3> &y,
                                  const Vec<T, 3> &z, Mat<T, 3, 3> &C) {
  for (int i = 0; i < 3; i++) {
    C(i, 0) = x(i);
    C(i, 1) = y(i);
    C(i, 2) = z(i);
  }
}

template <class xtype, class ytype, class ztype, class Ctype>
class MatFromThreeVecExpr {
 public:
  MatFromThreeVecExpr(xtype &x, ytype &y, ztype &z, Ctype &C)
      : x(x), y(y), z(z), C(C) {}

  A2D_FUNCTION void eval() {
    for (int i = 0; i < 3; i++) {
      C.value()(i, 0) = x.value()(i);
      C.value()(i, 1) = y.value()(i);
      C.value()(i, 2) = z.value()(i);
    }
  }

  A2D_FUNCTION void bzero() { C.bzero(); }

  template <ADorder forder>
  A2D_FUNCTION void forward() {
    constexpr ADseed seed = conditional_value < ADseed,
                     forder == ADorder::FIRST, ADseed::b, ADseed::p > ::value;
    auto &xp = GetSeed<seed>::get_obj(x);
    auto &yp = GetSeed<seed>::get_obj(y);
    auto &zp = GetSeed<seed>::get_obj(z);
    auto &Cp = GetSeed<seed>::get_obj(C);
    for (int i = 0; i < 3; i++) {
      Cp(i, 0) = xp(i);
      Cp(i, 1) = yp(i);
      Cp(i, 2) = zp(i);
    }
  }

  A2D_FUNCTION void reverse() {
    auto &xb = GetSeed<ADseed::b>::get_obj(x);
    auto &yb = GetSeed<ADseed::b>::get_obj(y);
    auto &zb = GetSeed<ADseed::b>::get_obj(z);
    auto &Cb = GetSeed<ADseed::b>::get_obj(C);
    for (int i = 0; i < 3; i++) {
      xb(i) += Cb(i, 0);
      yb(i) += Cb(i, 1);
      zb(i) += Cb(i, 2);
    }
  }

  A2D_FUNCTION void hzero() { C.hzero(); }

  A2D_FUNCTION void hreverse() {
    auto &xh = GetSeed<ADseed::h>::get_obj(x);
    auto &yh = GetSeed<ADseed::h>::get_obj(y);
    auto &zh = GetSeed<ADseed::h>::get_obj(z);
    auto &Ch = GetSeed<ADseed::h>::get_obj(C);
    for (int i = 0; i < 3; i++) {
      xh(i) += Ch(i, 0);
      yh(i) += Ch(i, 1);
      zh(i) += Ch(i, 2);
    }
  }

  xtype &x;
  ytype &y;
  ztype &z;
  Ctype &C;
};

template <class xtype, class ytype, class ztype, class Ctype>
A2D_FUNCTION auto MatFromThreeVec(ADObj<xtype> &x, ADObj<ytype> &y,
                                  ADObj<ztype> &z, ADObj<Ctype> &C) {
  return MatFromThreeVecExpr<ADObj<xtype>, ADObj<ytype>, ADObj<ztype>,
                             ADObj<Ctype>>(x, y, z, C);
}

template <class xtype, class ytype, class ztype, class Ctype>
A2D_FUNCTION auto MatFromThreeVec(A2DObj<xtype> &x, A2DObj<ytype> &y,
                                  A2DObj<ztype> &z, A2DObj<Ctype> &C) {
  return MatFromThreeVecExpr<A2DObj<xtype>, A2DObj<ytype>, A2DObj<ztype>,
                             A2DObj<Ctype>>(x, y, z, C);
}

}  // namespace A2D

#endif  // TACS_A2D_UTILITIES_H
