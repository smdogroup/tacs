// Op-level second-order verification for A2D::MatTrans3x3ADMatMult
// (feature-beam-element-methods, SPEC.md sec 1.2/sec 6.6). Standalone,
// header-only.
//
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_mattrans3x3admatmult_second_order.cpp -o t18_real && ./t18_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_mattrans3x3admatmult_second_order.cpp -o t18_complex && ./t18_complex
//
// C = scale*A^T*B is linear in its single active input B (A is a passive
// Mat3x3) -- same verification shape as ADMat3x3MatMult (Task 1.5).

#include "TACSObject.h"
#include "a2d.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(9753);

  const TacsScalar scale = 2.1;
  Mat3x3 A;
  ADMat3x3 B, C;
  for (int i = 0; i < 9; i++) {
    A.A[i] = randReal();
    B.A[i] = randReal();
  }
  TacsScalar Bp[9];
  for (int i = 0; i < 9; i++) Bp[i] = randReal();

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADMat3x3 Bt = B;
    for (int i = 0; i < 9; i++) Bt.Ap[i] = Bp[i];
    ADMat3x3 Ct;
    MatTrans3x3ADMatMult op(scale, A, Bt, Ct);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADMat3x3 Bh;
    for (int i = 0; i < 9; i++) Bh.A[i] = B.A[i] + TacsScalar(0.0, dh) * Bp[i];
    ADMat3x3 Ch;
    MatTrans3x3ADMatMult(scale, A, Bh, Ch);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsImagPart(Ch.A[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADMat3x3 Bph, Bmh;
    for (int i = 0; i < 9; i++) {
      Bph.A[i] = B.A[i] + dh * Bp[i];
      Bmh.A[i] = B.A[i] - dh * Bp[i];
    }
    ADMat3x3 Cph, Cmh;
    MatTrans3x3ADMatMult(scale, A, Bph, Cph);
    MatTrans3x3ADMatMult(scale, A, Bmh, Cmh);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsRealPart((Cph.A[i] - Cmh.A[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#endif
    printf(
        "  [MatTrans3x3ADMatMult] hforward vs "
        "directional-derivative-of-primal: %.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) zero-self-Hessian: C.Ah = 0, B.Ap arbitrary -> B.Ah bit-exact 0 --
  {
    ADMat3x3 Bt = B;
    for (int i = 0; i < 9; i++) Bt.Ap[i] = Bp[i];
    ADMat3x3 Ct;
    MatTrans3x3ADMatMult op(scale, A, Bt, Ct);
    op.hforward();
    // Ct.Ah left at its zero-initialized default.
    op.hreverse();

    double maxabs = 0.0;
    for (int i = 0; i < 9; i++) maxabs = std::max(maxabs, fabs(TacsRealPart(Bt.Ah[i])));
    printf(
        "  [MatTrans3x3ADMatMult] zero-self-Hessian check (expect exactly "
        "0.0): %.3e\n",
        maxabs);
    if (maxabs != 0.0) fail = 1;
  }

  // --- (3) hreverse vs hand-computed scale*A*Ch ------------------------------
  {
    TacsScalar Ch[9];
    for (int i = 0; i < 9; i++) Ch[i] = randReal();

    ADMat3x3 Bt = B;
    for (int i = 0; i < 9; i++) Bt.Ap[i] = Bp[i];
    ADMat3x3 Ct;
    MatTrans3x3ADMatMult op(scale, A, Bt, Ct);
    op.hforward();
    for (int i = 0; i < 9; i++) Ct.Ah[i] = Ch[i];
    op.hreverse();

    double expected[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        double s = 0.0;
        for (int j = 0; j < 3; j++) {
          s += TacsRealPart(A.A[3 * i + j]) * TacsRealPart(Ch[3 * j + k]);
        }
        expected[3 * i + k] += TacsRealPart(scale) * s;
      }
    }

    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(Bt.Ah[i]) - expected[i]));
    }
    printf(
        "  [MatTrans3x3ADMatMult] hreverse vs hand-computed scale*A*Ch: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (4) First-order forward()/reverse() call sites are unaffected --------
  {
    ADMat3x3 B0, C0;
    for (int i = 0; i < 9; i++) {
      B0.A[i] = randReal();
      B0.Ad[i] = randReal();
    }
    MatTrans3x3ADMatMult op0(A, B0, C0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
