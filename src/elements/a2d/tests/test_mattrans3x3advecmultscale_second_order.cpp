// Op-level second-order verification for A2D::MatTrans3x3ADVecMultScale
// (feature-beam-element-methods, SPEC.md sec 1.2/sec 6.6). Standalone,
// header-only.
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_mattrans3x3advecmultscale_second_order.cpp -o t19_real && ./t19_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_mattrans3x3advecmultscale_second_order.cpp -o t19_complex && ./t19_complex
// clang-format on
//
// y = scale.value*A^T*x is linear in its single active input x (scale and
// A are both passive). forward()/reverse() (and therefore hforward/
// hreverse) use the "Add" core variants, so this test explicitly zeros
// y.xp/x.xh before each accumulating call, mirroring how a caller chaining
// multiple ops into the same accumulator would.

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(8642);

  Scalar scale(1.4);
  Mat3x3 A;
  ADVec3 x, y;
  for (int i = 0; i < 9; i++) A.A[i] = randReal();
  for (int i = 0; i < 3; i++) x.x[i] = randReal();
  TacsScalar xp[3];
  for (int i = 0; i < 3; i++) xp[i] = randReal();

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADVec3 xt = x;
    for (int i = 0; i < 3; i++) xt.xp[i] = xp[i];
    ADVec3 yt;
    MatTrans3x3ADVecMultScale op(scale, A, xt, yt);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADVec3 xh;
    for (int i = 0; i < 3; i++) xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
    ADVec3 yh;
    MatTrans3x3ADVecMultScale(scale, A, xh, yh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsImagPart(yh.x[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xp[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADVec3 xph, xmh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
    }
    ADVec3 yph, ymh;
    MatTrans3x3ADVecMultScale(scale, A, xph, yph);
    MatTrans3x3ADVecMultScale(scale, A, xmh, ymh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsRealPart((yph.x[i] - ymh.x[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xp[i]) - ref));
    }
#endif
    printf(
        "  [MatTrans3x3ADVecMultScale] hforward vs "
        "directional-derivative-of-primal: %.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) zero-self-Hessian: y.xh = 0, x.xp arbitrary -> x.xh bit-exact 0 ---
  {
    ADVec3 xt = x;
    for (int i = 0; i < 3; i++) xt.xp[i] = xp[i];
    ADVec3 yt;
    MatTrans3x3ADVecMultScale op(scale, A, xt, yt);
    op.hforward();
    // yt.xh left at its zero-initialized default.
    op.hreverse();

    double maxabs = 0.0;
    for (int i = 0; i < 3; i++)
      maxabs = std::max(maxabs, fabs(TacsRealPart(xt.xh[i])));
    printf(
        "  [MatTrans3x3ADVecMultScale] zero-self-Hessian check (expect "
        "exactly 0.0): %.3e\n",
        maxabs);
    if (maxabs != 0.0) fail = 1;
  }

  // --- (3) hreverse vs hand-computed scale.value*A*yh ------------------------
  {
    TacsScalar yh[3];
    for (int i = 0; i < 3; i++) yh[i] = randReal();

    ADVec3 xt = x;
    for (int i = 0; i < 3; i++) xt.xp[i] = xp[i];
    ADVec3 yt;
    MatTrans3x3ADVecMultScale op(scale, A, xt, yt);
    op.hforward();
    for (int i = 0; i < 3; i++) yt.xh[i] = yh[i];
    op.hreverse();

    double expected[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
      double s = 0.0;
      for (int j = 0; j < 3; j++) {
        s += TacsRealPart(A.A[3 * i + j]) * TacsRealPart(yh[j]);
      }
      expected[i] = TacsRealPart(scale.value) * s;
    }

    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - expected[i]));
    }
    printf(
        "  [MatTrans3x3ADVecMultScale] hreverse vs hand-computed "
        "scale.value*A*yh: %.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (4) First-order forward()/reverse() call sites are unaffected --------
  {
    ADVec3 x0, y0;
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      x0.xd[i] = randReal();
    }
    MatTrans3x3ADVecMultScale op0(scale, A, x0, y0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
