// Op-level second-order verification for A2D::ADMatTrans3x3ADVecMultADScale
// (feature-beam-element-methods, SPEC-phase-7.md sec 2.2/sec 6.1).
// Standalone, header-only. Formulas derived by generalizing the bilinear
// "reverse-shape + one-cross-term" pattern E7 verified for
// ADMat3x3ADMatMult/ADMatTrans3x3ADMatMult to this op's genuinely TRILINEAR
// form (y = scale.value*A^T*x, all three of scale/A/x are AD-typed) -- NOT
// independently re-run through E7's own harness (SPEC-phase-7.md sec 2.2/sec
// 5) -- this test is that independent verification, with ALL THREE operands
// genuinely seeded (not the restricted zero-seed sub-case G4 actually uses,
// where scale/A stay Xpts-role-only -- this test verifies the op's general
// correctness, matching the convention every other op-level test in this
// directory already follows).
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admattrans3x3advecmultadscale_second_order.cpp -o t23_real && ./t23_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admattrans3x3advecmultadscale_second_order.cpp -o t23_complex && ./t23_complex
// clang-format on
//
// y = scale.value*A^T*x is TRILINEAR in (scale, A, x), so hforward is
// forward()'s own formula with the p-seed substituted, and hreverse is
// reverse()'s own formula (using .xh in place of .xd) PLUS TWO cross terms
// per input (one per OTHER factor, since there are two other factors in a
// trilinear form, unlike the bilinear ops' single cross term).

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(23561);

  ADScalar scale(0.83);
  ADMat3x3 A;
  ADVec3 x, y;
  for (int i = 0; i < 9; i++) A.A[i] = randReal();
  for (int i = 0; i < 3; i++) x.x[i] = randReal();

  TacsScalar scalep = randReal();
  TacsScalar Ap[9];
  for (int i = 0; i < 9; i++) Ap[i] = randReal();
  TacsScalar xp[3];
  for (int i = 0; i < 3; i++) xp[i] = randReal();

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADScalar scalet = scale;
    scalet.valuep = scalep;
    ADMat3x3 At = A;
    for (int i = 0; i < 9; i++) At.Ap[i] = Ap[i];
    ADVec3 xt = x;
    for (int i = 0; i < 3; i++) xt.xp[i] = xp[i];
    ADVec3 yt;
    ADMatTrans3x3ADVecMultADScale op(scalet, At, xt, yt);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADScalar scaleh(scale.value + TacsScalar(0.0, dh) * scalep);
    ADMat3x3 Ah;
    for (int i = 0; i < 9; i++) Ah.A[i] = A.A[i] + TacsScalar(0.0, dh) * Ap[i];
    ADVec3 xh;
    for (int i = 0; i < 3; i++) xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
    ADVec3 yh;
    ADMatTrans3x3ADVecMultADScale(scaleh, Ah, xh, yh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsImagPart(yh.x[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xp[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADScalar scaleph(scale.value + dh * scalep), scalemh(scale.value -
                                                          dh * scalep);
    ADMat3x3 Aph, Amh;
    for (int i = 0; i < 9; i++) {
      Aph.A[i] = A.A[i] + dh * Ap[i];
      Amh.A[i] = A.A[i] - dh * Ap[i];
    }
    ADVec3 xph, xmh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
    }
    ADVec3 yph, ymh;
    ADMatTrans3x3ADVecMultADScale(scaleph, Aph, xph, yph);
    ADMatTrans3x3ADVecMultADScale(scalemh, Amh, xmh, ymh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsRealPart((yph.x[i] - ymh.x[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xp[i]) - ref));
    }
#endif
    printf(
        "  [ADMatTrans3x3ADVecMultADScale] hforward vs "
        "directional-derivative-of-primal: %.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) hreverse (.valueh/.Ah/.xh) vs hand-computed formula ---------------
  {
    TacsScalar yd[3], yh_seed[3];
    for (int i = 0; i < 3; i++) {
      yd[i] = randReal();
      yh_seed[i] = randReal();
    }

    ADScalar scalet = scale;
    scalet.valuep = scalep;
    ADMat3x3 At = A;
    for (int i = 0; i < 9; i++) At.Ap[i] = Ap[i];
    ADVec3 xt = x;
    for (int i = 0; i < 3; i++) xt.xp[i] = xp[i];
    ADVec3 yt;
    ADMatTrans3x3ADVecMultADScale op(scalet, At, xt, yt);
    op.hforward();
    for (int i = 0; i < 3; i++) {
      yt.xd[i] = yd[i];
      yt.xh[i] = yh_seed[i];
    }
    op.reverse();
    op.hreverse();

    // Expected (trilinear y = scale.value*A^T*x):
    //   x.xh    += scale.value*A.A*y.xh   + scale.valuep*A.A*y.xd
    //            +  scale.value*A.Ap*y.xd
    //   A.Ah    += scale.value*outer(x.x, y.xh) + scale.valuep*outer(x.x, y.xd)
    //            + scale.value*outer(xp, y.xd)
    //   scale.valueh += (A.A^T*x.x).y.xh + (A.Ap^T*x.x).y.xd
    //                  + (A.A^T*xp).y.xd
    double sVal = TacsRealPart(scale.value);
    double sP = TacsRealPart(scalep);
    double expectedXh[3] = {0, 0, 0};
    double expectedAh[9] = {0};
    double expectedScaleH = 0.0;

    // x.xh: scale.value*A.A*y.xh + scale.valuep*A.A*y.xd + scale.value*A.Ap*y.xd
    // (A applied UN-transposed, matching Mat3x3VecMultAddScaleCore's own
    // convention -- the pullback of y=scale*A^T*x's linear map back to x is
    // scale*A, not scale*A^T; confirmed against the shipped reverse()).
    for (int i = 0; i < 3; i++) {
      double s1 = 0.0, s2 = 0.0, s3 = 0.0;
      for (int j = 0; j < 3; j++) {
        s1 += TacsRealPart(A.A[3 * i + j]) * TacsRealPart(yh_seed[j]);
        s2 += TacsRealPart(A.A[3 * i + j]) * TacsRealPart(yd[j]);
        s3 += TacsRealPart(Ap[3 * i + j]) * TacsRealPart(yd[j]);
      }
      expectedXh[i] = sVal * s1 + sP * s2 + sVal * s3;
    }
    // A.Ah: outer(x.x, y.xh)*scale.value + outer(x.x, y.xd)*scale.valuep +
    // outer(xp, y.xd)*scale.value  -- A.Ad = scale*outer(x, y) shape,
    // consistent with ADMat3x3VecVecInnerProduct/AXPY-style outer-product
    // convention used elsewhere: A[i,j] += scale * x[i] * y[j]
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double v1 = TacsRealPart(x.x[i]) * TacsRealPart(yh_seed[j]);
        double v2 = TacsRealPart(x.x[i]) * TacsRealPart(yd[j]);
        double v3 = TacsRealPart(xp[i]) * TacsRealPart(yd[j]);
        expectedAh[3 * i + j] = sVal * v1 + sP * v2 + sVal * v3;
      }
    }
    // scale.valueh
    for (int i = 0; i < 3; i++) {
      double Ax_i = 0.0, Apx_i = 0.0, Axp_i = 0.0;
      for (int j = 0; j < 3; j++) {
        Ax_i += TacsRealPart(A.A[3 * j + i]) * TacsRealPart(x.x[j]);
        Apx_i += TacsRealPart(Ap[3 * j + i]) * TacsRealPart(x.x[j]);
        Axp_i += TacsRealPart(A.A[3 * j + i]) * TacsRealPart(xp[j]);
      }
      expectedScaleH += Ax_i * TacsRealPart(yh_seed[i]);
      expectedScaleH += Apx_i * TacsRealPart(yd[i]);
      expectedScaleH += Axp_i * TacsRealPart(yd[i]);
    }

    double maxerr = fabs(TacsRealPart(scalet.valueh) - expectedScaleH);
    for (int i = 0; i < 3; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - expectedXh[i]));
    }
    for (int i = 0; i < 9; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(At.Ah[i]) - expectedAh[i]));
    }
    printf(
        "  [ADMatTrans3x3ADVecMultADScale] hreverse vs hand-computed "
        "formula: %.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (3) First-order forward()/reverse() call sites are unaffected --------
  {
    ADScalar scale0(0.4, 0.1);
    ADMat3x3 A0;
    ADVec3 x0, y0;
    for (int i = 0; i < 9; i++) {
      A0.A[i] = randReal();
      A0.Ad[i] = randReal();
    }
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      x0.xd[i] = randReal();
    }
    ADMatTrans3x3ADVecMultADScale op0(scale0, A0, x0, y0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
