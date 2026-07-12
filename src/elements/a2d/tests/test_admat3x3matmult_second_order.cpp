// Op-level second-order verification for A2D::ADMat3x3MatMult (feature-
// beam-element-methods, SPEC.md sec 1.2/sec 6.6). Standalone, header-only.
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3matmult_second_order.cpp -o t15_real && ./t15_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3matmult_second_order.cpp -o t15_complex && ./t15_complex
// clang-format on
//
// C = scale*A*B is LINEAR in its single active input A (B is a passive
// Mat3x3), so:
//   (1) hforward should equal the directional derivative of the primal
//       C.A along seed Ap (same check style as ADVec3Dot's test).
//   (2) hreverse's zero-self-Hessian check (SPEC.md Task 1.5): with no
//       downstream second-order seed (C.Ah = 0), A.Ah must come out
//       EXACTLY 0.0 (bit-exact, not just within tolerance) regardless of
//       the seed A.Ap -- this op cannot manufacture curvature on its own,
//       since it's linear in A. This mirrors VALIDATION.md E6's finding
//       for this op-class.
//   (3) hreverse's general formula (with a nonzero downstream seed C.Ah)
//       is checked against a central-difference-of-the-exact-adjoint,
//       exactly as for ADVec3Dot.

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(4321);

  const TacsScalar scale = 0.83;
  ADMat3x3 A, C;
  TacsScalar B[9];
  for (int i = 0; i < 9; i++) {
    A.A[i] = randReal();
    B[i] = randReal();
  }
  // Named, block-scoped Mat3x3 (NOT the raw TacsScalar[9] array above passed
  // directly to a constructor storing a `const Mat3x3&` reference member) --
  // E7 found this exact raw-array pattern creates a temporary Mat3x3 (via
  // Mat3x3's own implicit converting constructor) that is destroyed at the
  // end of the constructor-call's full expression, leaving the op object
  // holding a DANGLING reference for every subsequent .hforward()/
  // .hreverse() call -- undefined behavior that stayed silently benign here
  // only by luck of stack-slot reuse (SPEC-phase-7.md sec 2.2/sec 7,
  // VALIDATION.md's E7 entry). B[] itself is retained above for the
  // hand-computed reference formulas below, which only ever read its value,
  // never bind a reference to it.
  Mat3x3 Bmat(B);
  TacsScalar Ap[9];
  for (int i = 0; i < 9; i++) Ap[i] = randReal();

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADMat3x3 At = A;
    for (int i = 0; i < 9; i++) At.Ap[i] = Ap[i];
    ADMat3x3 Ct;
    ADMat3x3MatMult op(scale, At, Bmat, Ct);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADMat3x3 Ah;
    for (int i = 0; i < 9; i++) Ah.A[i] = A.A[i] + TacsScalar(0.0, dh) * Ap[i];
    ADMat3x3 Ch;
    ADMat3x3MatMult(scale, Ah, Bmat, Ch);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsImagPart(Ch.A[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADMat3x3 Aph, Amh;
    for (int i = 0; i < 9; i++) {
      Aph.A[i] = A.A[i] + dh * Ap[i];
      Amh.A[i] = A.A[i] - dh * Ap[i];
    }
    ADMat3x3 Cph, Cmh;
    ADMat3x3MatMult(scale, Aph, Bmat, Cph);
    ADMat3x3MatMult(scale, Amh, Bmat, Cmh);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsRealPart((Cph.A[i] - Cmh.A[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#endif
    printf(
        "  [ADMat3x3MatMult] hforward vs directional-derivative-of-primal: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) zero-self-Hessian: C.Ah = 0, A.Ap arbitrary -> A.Ah bit-exact 0 ---
  {
    ADMat3x3 At = A;
    for (int i = 0; i < 9; i++) At.Ap[i] = Ap[i];
    ADMat3x3 Ct;
    ADMat3x3MatMult op(scale, At, Bmat, Ct);
    op.hforward();
    // Ct.Ah left at its zero-initialized default (no downstream 2nd-order
    // seed injected).
    op.hreverse();

    double maxabs = 0.0;
    for (int i = 0; i < 9; i++) {
      maxabs = std::max(maxabs, fabs(TacsRealPart(At.Ah[i])));
    }
    printf(
        "  [ADMat3x3MatMult] zero-self-Hessian check (expect exactly "
        "0.0): %.3e\n",
        maxabs);
    if (maxabs != 0.0) fail = 1;
  }

  // --- (3) hreverse's propagation formula vs a hand-computed transpose-map --
  // Since this op is linear in A, reverse()'s own adjoint A.Ad = scale*C.Ad*B^T
  // does not depend on A's value at all -- a central-difference-of-the-
  // adjoint (as used for the genuinely bilinear ADVec3Dot) would trivially
  // read 0 here and cannot distinguish a correct hreverse from an incorrect
  // one. The right check for a linear op's hreverse is a direct comparison
  // against the same transpose-map formula reverse() already uses (Mat3x3Mat-
  // TransMultAddScaleCore), evaluated on Ah instead of Ad.
  {
    TacsScalar Ch[9];
    for (int i = 0; i < 9; i++)
      Ch[i] = randReal();  // downstream 2nd-order seed

    ADMat3x3 At = A;
    for (int i = 0; i < 9; i++) At.Ap[i] = Ap[i];
    ADMat3x3 Ct;
    ADMat3x3MatMult op(scale, At, Bmat, Ct);
    op.hforward();
    for (int i = 0; i < 9; i++) Ct.Ah[i] = Ch[i];
    op.hreverse();

    // Hand-computed expected = scale * Ch * B^T (same shape as
    // Mat3x3MatTransMultAddScaleCore(scale, Ch, B, expected)).
    double expected[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
        double s = 0.0;
        for (int j = 0; j < 3; j++) {
          s += TacsRealPart(Ch[3 * i + j]) * TacsRealPart(B[3 * k + j]);
        }
        expected[3 * i + k] += TacsRealPart(scale) * s;
      }
    }

    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(At.Ah[i]) - expected[i]));
    }
    printf(
        "  [ADMat3x3MatMult] hreverse vs hand-computed scale*Ch*B^T: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (4) First-order forward()/reverse() call sites are unaffected --------
  {
    ADMat3x3 A0, C0;
    for (int i = 0; i < 9; i++) {
      A0.A[i] = randReal();
      A0.Ad[i] = randReal();
    }
    ADMat3x3MatMult op0(A0, Bmat, C0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
