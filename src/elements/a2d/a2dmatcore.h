#ifndef A2D_MAT_CORE_H
#define A2D_MAT_CORE_H

namespace A2D {

inline TacsScalar Mat3x3DetCore(const TacsScalar A[]) {
  return (A[8] * (A[0] * A[4] - A[3] * A[1]) -
          A[7] * (A[0] * A[5] - A[3] * A[2]) +
          A[6] * (A[1] * A[5] - A[2] * A[4]));
}

inline TacsScalar Mat3x3DetDerivForwardCore(const TacsScalar A[],
                                            TacsScalar Ad[]) {
  return (Ad[0] * (A[8] * A[4] - A[7] * A[5]) +
          Ad[1] * (A[6] * A[5] - A[8] * A[3]) +
          Ad[2] * (A[7] * A[3] - A[6] * A[4]) +
          Ad[3] * (A[7] * A[2] - A[8] * A[1]) +
          Ad[4] * (A[8] * A[0] - A[6] * A[2]) +
          Ad[5] * (A[6] * A[1] - A[7] * A[0]) +
          Ad[6] * (A[1] * A[5] - A[2] * A[4]) +
          Ad[7] * (A[3] * A[2] - A[0] * A[5]) +
          Ad[8] * (A[0] * A[4] - A[3] * A[1]));
}

inline void Mat3x3DetDerivReverseCore(const TacsScalar detd,
                                      const TacsScalar A[], TacsScalar Ad[]) {
  Ad[0] += (A[8] * A[4] - A[7] * A[5]) * detd;
  Ad[1] += (A[6] * A[5] - A[8] * A[3]) * detd;
  Ad[2] += (A[7] * A[3] - A[6] * A[4]) * detd;
  Ad[3] += (A[7] * A[2] - A[8] * A[1]) * detd;
  Ad[4] += (A[8] * A[0] - A[6] * A[2]) * detd;
  Ad[5] += (A[6] * A[1] - A[7] * A[0]) * detd;
  Ad[6] += (A[1] * A[5] - A[2] * A[4]) * detd;
  Ad[7] += (A[3] * A[2] - A[0] * A[5]) * detd;
  Ad[8] += (A[0] * A[4] - A[3] * A[1]) * detd;
}

inline TacsScalar Symm3x3DetCore(const TacsScalar S[]) {
  return (S[5] * (S[0] * S[3] - S[1] * S[1]) -
          S[4] * (S[0] * S[4] - S[1] * S[2]) +
          S[2] * (S[1] * S[4] - S[2] * S[3]));
}

inline TacsScalar Symm3x3DetDerivForwardCore(const TacsScalar S[],
                                             const TacsScalar Sd[]) {
  return (Sd[0] * (S[5] * S[3] - S[4] * S[4]) +
          2.0 * Sd[1] * (S[2] * S[4] - S[5] * S[1]) +
          2.0 * Sd[2] * (S[1] * S[4] - S[3] * S[2]) +
          Sd[3] * (S[5] * S[0] - S[2] * S[2]) +
          2.0 * Sd[4] * (S[1] * S[2] - S[0] * S[4]) +
          Sd[5] * (S[0] * S[3] - S[1] * S[1]));
}

inline void Symm3x3DetDerivReverseCore(const TacsScalar detd,
                                       const TacsScalar S[], TacsScalar Sd[]) {
  Sd[0] += (S[5] * S[3] - S[4] * S[4]) * detd;
  Sd[1] += 2.0 * (S[2] * S[4] - S[5] * S[1]) * detd;
  Sd[2] += 2.0 * (S[1] * S[4] - S[3] * S[2]) * detd;
  Sd[3] += (S[5] * S[0] - S[2] * S[2]) * detd;
  Sd[4] += 2.0 * (S[1] * S[2] - S[0] * S[4]) * detd;
  Sd[5] += (S[0] * S[3] - S[1] * S[1]) * detd;
}

inline void Symm3x3SymmMultCore(const TacsScalar A[], const TacsScalar B[],
                                TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] = (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] = (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] = (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                     const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = scale * (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] = scale * (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = scale * (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] = scale * (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultAddCore(const TacsScalar A[], const TacsScalar B[],
                                   TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] += (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] += (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] += (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultSubCore(const TacsScalar A[], const TacsScalar B[],
                                   TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] -= (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] -= (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] -= (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] -= (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] -= (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] -= (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] -= (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3SymmMultAddScaleCore(TacsScalar scale, const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += scale * (A[1] * B[1] + A[3] * B[3] + A[4] * B[4]);
  C[5] += scale * (A[1] * B[2] + A[3] * B[4] + A[4] * B[5]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += scale * (A[2] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[8] += scale * (A[2] * B[2] + A[4] * B[4] + A[5] * B[5]);
}

inline void Symm3x3MatMultCore(const TacsScalar A[], const TacsScalar B[],
                               TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] = (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] = (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] = (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] = (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] = (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                    const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] = scale * (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] = scale * (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] = scale * (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] = scale * (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultAddCore(const TacsScalar A[], const TacsScalar B[],
                                  TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] += (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] += (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] += (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] += (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] += (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultSubCore(const TacsScalar A[], const TacsScalar B[],
                                  TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] -= (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] -= (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] -= (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] -= (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] -= (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] -= (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] -= (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] -= (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatMultAddScaleCore(TacsScalar scale, const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[3] + A[4] * B[6]);
  C[4] += scale * (A[1] * B[1] + A[3] * B[4] + A[4] * B[7]);
  C[5] += scale * (A[1] * B[2] + A[3] * B[5] + A[4] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[7] += scale * (A[2] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[8] += scale * (A[2] * B[2] + A[4] * B[5] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultCore(const TacsScalar A[], const TacsScalar B[],
                                    TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] = (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] = (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] = (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                         const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] = scale * (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] = scale * (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] = scale * (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] = scale * (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultAddCore(const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] += (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] += (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] += (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultSubCore(const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] -= (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] -= (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] -= (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] -= (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] -= (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] -= (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] -= (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Symm3x3MatTransMultAddScaleCore(TacsScalar scale,
                                            const TacsScalar A[],
                                            const TacsScalar B[],
                                            TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[3] * B[1] + A[4] * B[2]);
  C[4] += scale * (A[1] * B[3] + A[3] * B[4] + A[4] * B[5]);
  C[5] += scale * (A[1] * B[6] + A[3] * B[7] + A[4] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[7] += scale * (A[2] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[8] += scale * (A[2] * B[6] + A[4] * B[7] + A[5] * B[8]);
}

inline void Mat3x3SymmMultCore(const TacsScalar A[], const TacsScalar B[],
                               TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] = (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] = (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] = (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                    const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = scale * (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] = scale * (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] = scale * (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultAddCore(const TacsScalar A[], const TacsScalar B[],
                                  TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] += (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] += (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] += (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultSubCore(const TacsScalar A[], const TacsScalar B[],
                                  TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] -= (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] -= (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] -= (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] -= (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] -= (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] -= (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] -= (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void Mat3x3SymmMultAddScaleCore(TacsScalar scale, const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[3] + A[2] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[4] + A[2] * B[5]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += scale * (A[3] * B[1] + A[4] * B[3] + A[5] * B[4]);
  C[5] += scale * (A[3] * B[2] + A[4] * B[4] + A[5] * B[5]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[6] * B[1] + A[7] * B[3] + A[8] * B[4]);
  C[8] += scale * (A[6] * B[2] + A[7] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultCore(const TacsScalar A[], const TacsScalar B[],
                                    TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] = (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] = (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] = (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] = (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] = (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                         const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = scale * (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] = scale * (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = scale * (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] = scale * (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] = scale * (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultAddCore(const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] += (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] += (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] += (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] += (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] += (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultSubCore(const TacsScalar A[],
                                       const TacsScalar B[], TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] -= (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] -= (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] -= (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] -= (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] -= (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] -= (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] -= (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] -= (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void MatTrans3x3SymmMultAddScaleCore(TacsScalar scale,
                                            const TacsScalar A[],
                                            const TacsScalar B[],
                                            TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += scale * (A[0] * B[1] + A[3] * B[3] + A[6] * B[4]);
  C[2] += scale * (A[0] * B[2] + A[3] * B[4] + A[6] * B[5]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += scale * (A[1] * B[1] + A[4] * B[3] + A[7] * B[4]);
  C[5] += scale * (A[1] * B[2] + A[4] * B[4] + A[7] * B[5]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[2] * B[1] + A[5] * B[3] + A[8] * B[4]);
  C[8] += scale * (A[2] * B[2] + A[5] * B[4] + A[8] * B[5]);
}

inline void Mat3x3MatMultCore(const TacsScalar A[], const TacsScalar B[],
                              TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] = (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] = (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] = (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] = (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] = (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                   const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] = scale * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] = scale * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] = scale * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] = scale * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultAddCore(const TacsScalar A[], const TacsScalar B[],
                                 TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] += (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] += (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] += (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] += (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] += (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultSubCore(const TacsScalar A[], const TacsScalar B[],
                                 TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] -= (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] -= (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] -= (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] -= (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] -= (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] -= (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] -= (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] -= (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatMultAddScaleCore(TacsScalar scale, const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[3] + A[2] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[1] * B[4] + A[2] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[1] * B[5] + A[2] * B[8]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[3] + A[5] * B[6]);
  C[4] += scale * (A[3] * B[1] + A[4] * B[4] + A[5] * B[7]);
  C[5] += scale * (A[3] * B[2] + A[4] * B[5] + A[5] * B[8]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[3] + A[8] * B[6]);
  C[7] += scale * (A[6] * B[1] + A[7] * B[4] + A[8] * B[7]);
  C[8] += scale * (A[6] * B[2] + A[7] * B[5] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultCore(const TacsScalar A[], const TacsScalar B[],
                                   TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] = (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] = (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] = (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] = scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] = scale * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] = scale * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] = scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] = scale * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultAddCore(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] += (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] += (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] += (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultSubCore(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] -= (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] -= (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] -= (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] -= (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] -= (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] -= (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] -= (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] -= (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void Mat3x3MatTransMultAddScaleCore(TacsScalar scale,
                                           const TacsScalar A[],
                                           const TacsScalar B[],
                                           TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[1] * B[4] + A[2] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[1] * B[7] + A[2] * B[8]);
  C[3] += scale * (A[3] * B[0] + A[4] * B[1] + A[5] * B[2]);
  C[4] += scale * (A[3] * B[3] + A[4] * B[4] + A[5] * B[5]);
  C[5] += scale * (A[3] * B[6] + A[4] * B[7] + A[5] * B[8]);
  C[6] += scale * (A[6] * B[0] + A[7] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[6] * B[3] + A[7] * B[4] + A[8] * B[5]);
  C[8] += scale * (A[6] * B[6] + A[7] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultCore(const TacsScalar A[], const TacsScalar B[],
                                   TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] = (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] = (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] = (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] = (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] = (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] = (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] = (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] = (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultScaleCore(TacsScalar scale, const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] = scale * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] = scale * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] = scale * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] = scale * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] = scale * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] = scale * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultAddCore(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] += (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] += (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] += (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] += (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] += (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] += (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] += (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] += (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultSubCore(const TacsScalar A[],
                                      const TacsScalar B[], TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] -= (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] -= (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] -= (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] -= (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] -= (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] -= (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] -= (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] -= (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatMultAddScaleCore(TacsScalar scale,
                                           const TacsScalar A[],
                                           const TacsScalar B[],
                                           TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[3] * B[3] + A[6] * B[6]);
  C[1] += scale * (A[0] * B[1] + A[3] * B[4] + A[6] * B[7]);
  C[2] += scale * (A[0] * B[2] + A[3] * B[5] + A[6] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[3] + A[7] * B[6]);
  C[4] += scale * (A[1] * B[1] + A[4] * B[4] + A[7] * B[7]);
  C[5] += scale * (A[1] * B[2] + A[4] * B[5] + A[7] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[3] + A[8] * B[6]);
  C[7] += scale * (A[2] * B[1] + A[5] * B[4] + A[8] * B[7]);
  C[8] += scale * (A[2] * B[2] + A[5] * B[5] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultCore(const TacsScalar A[],
                                        const TacsScalar B[], TacsScalar C[]) {
  C[0] = (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] = (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] = (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] = (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] = (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] = (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultScaleCore(TacsScalar scale,
                                             const TacsScalar A[],
                                             const TacsScalar B[],
                                             TacsScalar C[]) {
  C[0] = scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] = scale * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] = scale * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] = scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] = scale * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] = scale * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] = scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] = scale * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] = scale * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultAddCore(const TacsScalar A[],
                                           const TacsScalar B[],
                                           TacsScalar C[]) {
  C[0] += (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] += (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] += (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] += (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] += (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] += (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultSubCore(const TacsScalar A[],
                                           const TacsScalar B[],
                                           TacsScalar C[]) {
  C[0] -= (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] -= (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] -= (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] -= (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] -= (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] -= (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] -= (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] -= (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] -= (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline void MatTrans3x3MatTransMultAddScaleCore(TacsScalar scale,
                                                const TacsScalar A[],
                                                const TacsScalar B[],
                                                TacsScalar C[]) {
  C[0] += scale * (A[0] * B[0] + A[3] * B[1] + A[6] * B[2]);
  C[1] += scale * (A[0] * B[3] + A[3] * B[4] + A[6] * B[5]);
  C[2] += scale * (A[0] * B[6] + A[3] * B[7] + A[6] * B[8]);
  C[3] += scale * (A[1] * B[0] + A[4] * B[1] + A[7] * B[2]);
  C[4] += scale * (A[1] * B[3] + A[4] * B[4] + A[7] * B[5]);
  C[5] += scale * (A[1] * B[6] + A[4] * B[7] + A[7] * B[8]);
  C[6] += scale * (A[2] * B[0] + A[5] * B[1] + A[8] * B[2]);
  C[7] += scale * (A[2] * B[3] + A[5] * B[4] + A[8] * B[5]);
  C[8] += scale * (A[2] * B[6] + A[5] * B[7] + A[8] * B[8]);
}

inline TacsScalar Mat3x3InverseCore(const TacsScalar A[], TacsScalar Ainv[]) {
  TacsScalar det =
      (A[8] * (A[0] * A[4] - A[3] * A[1]) - A[7] * (A[0] * A[5] - A[3] * A[2]) +
       A[6] * (A[1] * A[5] - A[2] * A[4]));
  TacsScalar detinv = 1.0 / det;

  Ainv[0] = (A[4] * A[8] - A[5] * A[7]) * detinv;
  Ainv[1] = -(A[1] * A[8] - A[2] * A[7]) * detinv;
  Ainv[2] = (A[1] * A[5] - A[2] * A[4]) * detinv;

  Ainv[3] = -(A[3] * A[8] - A[5] * A[6]) * detinv;
  Ainv[4] = (A[0] * A[8] - A[2] * A[6]) * detinv;
  Ainv[5] = -(A[0] * A[5] - A[2] * A[3]) * detinv;

  Ainv[6] = (A[3] * A[7] - A[4] * A[6]) * detinv;
  Ainv[7] = -(A[0] * A[7] - A[1] * A[6]) * detinv;
  Ainv[8] = (A[0] * A[4] - A[1] * A[3]) * detinv;

  return det;
}

inline void Mat3x3InverseDerivForwardCore(const TacsScalar Ainv[],
                                          const TacsScalar Ad[],
                                          TacsScalar Bd[]) {
  TacsScalar t[9];
  Mat3x3MatMultCore(Ainv, Ad, t);
  Mat3x3MatMultScaleCore(-1.0, t, Ainv, Bd);
}

inline void Mat3x3InverseDerivReverseCore(const TacsScalar Ainv[],
                                          const TacsScalar Bd[],
                                          TacsScalar Ad[]) {
  TacsScalar t[9];
  MatTrans3x3MatMultCore(Ainv, Bd, t);
  Mat3x3MatTransMultAddScaleCore(-1.0, t, Ainv, Ad);
}

inline TacsScalar Symm3x3MatMultTraceCore(const TacsScalar S[],
                                          const TacsScalar T[]) {
  return ((S[0] * T[0] + S[3] * T[3] + S[5] * T[5]) +
          2.0 * (S[1] * T[1] + S[2] * T[2] + S[4] * T[4]));
}

inline void Symm3x3MatMultTraceReverseCore(const TacsScalar scale,
                                           const TacsScalar S[],
                                           TacsScalar Td[]) {
  Td[0] += scale * S[0];
  Td[3] += scale * S[3];
  Td[5] += scale * S[5];
  Td[1] += 2.0 * scale * S[1];
  Td[2] += 2.0 * scale * S[2];
  Td[4] += 2.0 * scale * S[4];
}

inline void Mat3x3LinearGreenStrainCore(const TacsScalar Ux[], TacsScalar E[]) {
  // E = 0.5*(Ux + Ux^{T});
  E[0] = Ux[0];
  E[3] = Ux[4];
  E[5] = Ux[8];
  E[1] = 0.5 * (Ux[1] + Ux[3]);
  E[2] = 0.5 * (Ux[2] + Ux[6]);
  E[4] = 0.5 * (Ux[5] + Ux[7]);
}

inline void Mat3x3LinearGreenStrainReverseCore(const TacsScalar Ed[],
                                               TacsScalar Uxd[]) {
  Uxd[0] += Ed[0];
  Uxd[1] += Ed[1];
  Uxd[2] += Ed[2];
  Uxd[3] += Ed[1];
  Uxd[4] += Ed[3];
  Uxd[5] += Ed[4];
  Uxd[6] += Ed[2];
  Uxd[7] += Ed[4];
  Uxd[8] += Ed[5];
}

inline void Mat3x3GreenStrainCore(const TacsScalar Ux[], TacsScalar E[]) {
  // E = 0.5*(Ux + Ux^{T} + Ux^{T} * Ux)
  E[0] = Ux[0] + 0.5 * (Ux[0] * Ux[0] + Ux[3] * Ux[3] + Ux[6] * Ux[6]);
  E[3] = Ux[4] + 0.5 * (Ux[1] * Ux[1] + Ux[4] * Ux[4] + Ux[7] * Ux[7]);
  E[5] = Ux[8] + 0.5 * (Ux[2] * Ux[2] + Ux[5] * Ux[5] + Ux[8] * Ux[8]);
  E[1] = 0.5 * (Ux[1] + Ux[3] + Ux[0] * Ux[1] + Ux[3] * Ux[4] + Ux[6] * Ux[7]);
  E[2] = 0.5 * (Ux[2] + Ux[6] + Ux[0] * Ux[2] + Ux[3] * Ux[5] + Ux[6] * Ux[8]);
  E[4] = 0.5 * (Ux[5] + Ux[7] + Ux[1] * Ux[2] + Ux[4] * Ux[5] + Ux[7] * Ux[8]);
}

inline void Mat3x3GreenStrainForwardCore(const TacsScalar Ux[],
                                         const TacsScalar Uxd[],
                                         TacsScalar Ed[]) {
  Ed[0] = Uxd[0] + Ux[0] * Uxd[0] + Ux[3] * Uxd[3] + Ux[6] * Uxd[6];
  Ed[3] = Uxd[4] + Ux[1] * Uxd[1] + Ux[4] * Uxd[4] + Ux[7] * Uxd[7];
  Ed[5] = Uxd[8] + Ux[2] * Uxd[2] + Ux[5] * Uxd[5] + Ux[8] * Uxd[8];

  Ed[1] =
      0.5 * (Uxd[1] + Uxd[3] + Ux[0] * Uxd[1] + Ux[3] * Uxd[4] +
             Ux[6] * Uxd[7] + Uxd[0] * Ux[1] + Uxd[3] * Ux[4] + Uxd[6] * Ux[7]);
  Ed[2] =
      0.5 * (Uxd[2] + Uxd[6] + Ux[0] * Uxd[2] + Ux[3] * Uxd[5] +
             Ux[6] * Uxd[8] + Uxd[0] * Ux[2] + Uxd[3] * Ux[5] + Uxd[6] * Ux[8]);
  Ed[4] =
      0.5 * (Uxd[5] + Uxd[7] + Ux[1] * Uxd[2] + Ux[4] * Uxd[5] +
             Ux[7] * Uxd[8] + Uxd[1] * Ux[2] + Uxd[4] * Ux[5] + Uxd[7] * Ux[8]);
}

inline void Mat3x3GreenStrainReverseCore(const TacsScalar Ux[],
                                         const TacsScalar Ed[],
                                         TacsScalar Uxd[]) {
  TacsScalar ux0 = Ux[0] + 1.0;
  TacsScalar ux4 = Ux[4] + 1.0;
  TacsScalar ux8 = Ux[8] + 1.0;
  TacsScalar e1 = 0.5 * Ed[1];
  TacsScalar e2 = 0.5 * Ed[2];
  TacsScalar e4 = 0.5 * Ed[4];

  // Uxd = (I + Ux) * E
  Uxd[0] += ux0 * Ed[0] + Ux[1] * e1 + Ux[2] * e2;
  Uxd[1] += ux0 * e1 + Ux[1] * Ed[3] + Ux[2] * e4;
  Uxd[2] += ux0 * e2 + Ux[1] * e4 + Ux[2] * Ed[5];
  Uxd[3] += Ux[3] * Ed[0] + ux4 * e1 + Ux[5] * e2;
  Uxd[4] += Ux[3] * e1 + ux4 * Ed[3] + Ux[5] * e4;
  Uxd[5] += Ux[3] * e2 + ux4 * e4 + Ux[5] * Ed[5];
  Uxd[6] += Ux[6] * Ed[0] + Ux[7] * e1 + ux8 * e2;
  Uxd[7] += Ux[6] * e1 + Ux[7] * Ed[3] + ux8 * e4;
  Uxd[8] += Ux[6] * e2 + Ux[7] * e4 + ux8 * Ed[5];
}

inline void Symm3x3IsotropicConstitutiveCore(const TacsScalar mu,
                                             const TacsScalar lambda,
                                             const TacsScalar E[],
                                             TacsScalar S[]) {
  TacsScalar tr = lambda * (E[0] + E[3] + E[5]);
  TacsScalar mu2 = 2.0 * mu;
  S[0] = mu2 * E[0] + tr;
  S[1] = mu2 * E[1];
  S[2] = mu2 * E[2];
  S[3] = mu2 * E[3] + tr;
  S[4] = mu2 * E[4];
  S[5] = mu2 * E[5] + tr;
}

inline void Symm3x3IsotropicConstitutiveReverseCore(const TacsScalar mu,
                                                    const TacsScalar lambda,
                                                    const TacsScalar Sd[],
                                                    TacsScalar Ed[]) {
  TacsScalar tr = lambda * (Sd[0] + Sd[3] + Sd[5]);
  TacsScalar mu2 = 2.0 * mu;
  Ed[0] += mu2 * Sd[0] + tr;
  Ed[1] += mu2 * Sd[1];
  Ed[2] += mu2 * Sd[2];
  Ed[3] += mu2 * Sd[3] + tr;
  Ed[4] += mu2 * Sd[4];
  Ed[5] += mu2 * Sd[5] + tr;
}

}  // namespace A2D

#endif  // A2D_MAT_CORE_H
