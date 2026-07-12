// Op-level second-order verification for A2D::ADMat3x3ADMatMult
// (feature-beam-element-methods, SPEC-phase-7.md sec 2.2/sec 6.1). Standalone,
// header-only. Formulas ported from the E7 experiment's verified prototype
// (docs/plans/feature-beam-element-methods/scripts/004_e7_zero_seed_role.cpp,
// TestADMat3x3ADMatMult, scale=1 case), generalized here to also cover the
// scaled constructor (mirrors ADMat3x3MatMult's own TacsRealPart(scale)==1.0
// branch convention).
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3admatmult_second_order.cpp -o t20_real && ./t20_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3admatmult_second_order.cpp -o t20_complex && ./t20_complex
// clang-format on
//
// C = scale*A*B is bilinear in its two active inputs (A, B), so both
// hforward and hreverse carry cross terms between A and B (same shape as
// ADVec3Dot's own bilinear test).

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(20233);

  const TacsScalar scale = 1.37;
  ADMat3x3 A, B, C;
  for (int i = 0; i < 9; i++) {
    A.A[i] = randReal();
    B.A[i] = randReal();
  }
  TacsScalar Ap[9], Bp[9];
  for (int i = 0; i < 9; i++) {
    Ap[i] = randReal();
    Bp[i] = randReal();
  }

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADMat3x3 At = A, Bt = B;
    for (int i = 0; i < 9; i++) {
      At.Ap[i] = Ap[i];
      Bt.Ap[i] = Bp[i];
    }
    ADMat3x3 Ct;
    ADMat3x3ADMatMult op(scale, At, Bt, Ct);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADMat3x3 Ah, Bh;
    for (int i = 0; i < 9; i++) {
      Ah.A[i] = A.A[i] + TacsScalar(0.0, dh) * Ap[i];
      Bh.A[i] = B.A[i] + TacsScalar(0.0, dh) * Bp[i];
    }
    ADMat3x3 Ch;
    ADMat3x3ADMatMult(scale, Ah, Bh, Ch);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsImagPart(Ch.A[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADMat3x3 Aph, Bph, Amh, Bmh;
    for (int i = 0; i < 9; i++) {
      Aph.A[i] = A.A[i] + dh * Ap[i];
      Bph.A[i] = B.A[i] + dh * Bp[i];
      Amh.A[i] = A.A[i] - dh * Ap[i];
      Bmh.A[i] = B.A[i] - dh * Bp[i];
    }
    ADMat3x3 Cph, Cmh;
    ADMat3x3ADMatMult(scale, Aph, Bph, Cph);
    ADMat3x3ADMatMult(scale, Amh, Bmh, Cmh);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsRealPart((Cph.A[i] - Cmh.A[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#endif
    printf(
        "  [ADMat3x3ADMatMult] hforward vs directional-derivative-of-primal: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) hreverse (.Ah) vs hand-computed formula ---------------------------
  {
    TacsScalar Cd[9], Ch[9];
    for (int i = 0; i < 9; i++) {
      Cd[i] = randReal();
      Ch[i] = randReal();
    }

    ADMat3x3 At = A, Bt = B;
    for (int i = 0; i < 9; i++) {
      At.Ap[i] = Ap[i];
      Bt.Ap[i] = Bp[i];
    }
    ADMat3x3 Ct;
    ADMat3x3ADMatMult op(scale, At, Bt, Ct);
    op.hforward();
    for (int i = 0; i < 9; i++) {
      Ct.Ad[i] = Cd[i];
      Ct.Ah[i] = Ch[i];
    }
    op.reverse();
    op.hreverse();

    // Expected: A.Ah = scale*(Ch*B.A^T + Cd*Bp^T)
    //           B.Ah = scale*(A.A^T*Ch + Ap^T*Cd)
    double expectedAh[9] = {0}, expectedBh[9] = {0};
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        double sA = 0.0, sB = 0.0;
        for (int j = 0; j < 3; j++) {
          sA += TacsRealPart(Ch[3 * i + j]) * TacsRealPart(B.A[3 * k + j]);
          sA += TacsRealPart(Cd[3 * i + j]) * TacsRealPart(Bp[3 * k + j]);
          sB += TacsRealPart(A.A[3 * j + i]) * TacsRealPart(Ch[3 * j + k]);
          sB += TacsRealPart(Ap[3 * j + i]) * TacsRealPart(Cd[3 * j + k]);
        }
        expectedAh[3 * i + k] = TacsRealPart(scale) * sA;
        expectedBh[3 * i + k] = TacsRealPart(scale) * sB;
      }
    }

    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(At.Ah[i]) - expectedAh[i]));
      maxerr = std::max(maxerr, fabs(TacsRealPart(Bt.Ah[i]) - expectedBh[i]));
    }
    printf(
        "  [ADMat3x3ADMatMult] hreverse (.Ah) vs hand-computed formula: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (3) First-order forward()/reverse() call sites are unaffected --------
  {
    ADMat3x3 A0, B0, C0;
    for (int i = 0; i < 9; i++) {
      A0.A[i] = randReal();
      A0.Ad[i] = randReal();
      B0.A[i] = randReal();
      B0.Ad[i] = randReal();
    }
    ADMat3x3ADMatMult op0(A0, B0, C0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
