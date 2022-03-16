#ifndef A2D_VEC_CORE_H
#define A2D_VEC_CORE_H

namespace A2D {

inline TacsScalar Vec3DotCore( const TacsScalar x[],
                               const TacsScalar y[] ){
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

inline TacsScalar Vec3NormCore( const TacsScalar x[] ){
  return sqrt(Vec3DotCore(x, x));
}

inline void Vec3ScaleCore( const TacsScalar alpha, const TacsScalar x[], TacsScalar v[] ){
  v[0] = alpha * x[0];
  v[1] = alpha * x[1];
  v[2] = alpha * x[2];
}

inline void Vec3AXPYCore( const TacsScalar alpha, const TacsScalar x[],
                          const TacsScalar y[], TacsScalar v[] ){
  v[0] = alpha * x[0] + y[0];
  v[1] = alpha * x[1] + y[1];
  v[2] = alpha * x[2] + y[2];
}

inline void Vec3CrossProductCore( const TacsScalar x[], const TacsScalar y[],
                                  TacsScalar v[] ){

  v[0] = x[1] * y[2] - x[2] * y[1];
  v[1] = x[2] * y[0] - x[0] * y[2];
  v[2] = x[0] * y[1] - x[1] * y[0];
}

inline void Vec3CrossProductAddCore( const TacsScalar x[], const TacsScalar y[],
                                     TacsScalar v[] ){

  v[0] += x[1] * y[2] - x[2] * y[1];
  v[1] += x[2] * y[0] - x[0] * y[2];
  v[2] += x[0] * y[1] - x[1] * y[0];
}

} // namespace AD

#endif // A2D_VEC_CORE_H
