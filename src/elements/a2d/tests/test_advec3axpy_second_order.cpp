// Op-level second-order verification for A2D::ADVec3Axpy (feature-beam-
// element-methods, SPEC-phase-7.md sec 2.2/sec 6.1). Standalone,
// header-only. Formulas derived by the same pattern already E7-verified for
// ADMat3x3ADMatMult/ADMatTrans3x3ADMatMult, NOT independently re-run through
// E7's own harness (SPEC-phase-7.md sec 2.2/sec 5) -- this test is that
// independent verification.
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3axpy_second_order.cpp -o t22_real && ./t22_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3axpy_second_order.cpp -o t22_complex && ./t22_complex
// clang-format on
//
// v = scale*alpha.value*x + y is bilinear in (alpha, x) and purely additive
// (linear, no cross term) in y -- so hforward is forward()'s own formula
// with the p-seed substituted, and hreverse is reverse()'s own formula
// (using .xh in place of .xd) plus ONE cross term for the (alpha, x)
// bilinear pair (y contributes no cross term since d^2v/d(alpha)d(y) =
// d^2v/dx*dy = 0 -- y enters v additively, not multiplicatively).

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(22459);

  const TacsScalar scale = -1.0;
  ADScalar alpha(0.71);
  ADVec3 x, y, v;
  for (int i = 0; i < 3; i++) {
    x.x[i] = randReal();
    y.x[i] = randReal();
  }
  TacsScalar alphap = randReal();
  TacsScalar xp[3], yp[3];
  for (int i = 0; i < 3; i++) {
    xp[i] = randReal();
    yp[i] = randReal();
  }

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADScalar alphat = alpha;
    alphat.valuep = alphap;
    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADVec3 vt;
    ADVec3Axpy op(scale, alphat, xt, yt, vt);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADScalar alphah(alpha.value + TacsScalar(0.0, dh) * alphap);
    ADVec3 xh, yh;
    for (int i = 0; i < 3; i++) {
      xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
      yh.x[i] = y.x[i] + TacsScalar(0.0, dh) * yp[i];
    }
    ADVec3 vh;
    ADVec3Axpy(scale, alphah, xh, yh, vh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsImagPart(vh.x[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(vt.xp[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADScalar alphaph(alpha.value + dh * alphap),
        alphamh(alpha.value - dh * alphap);
    ADVec3 xph, yph, xmh, ymh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      yph.x[i] = y.x[i] + dh * yp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
      ymh.x[i] = y.x[i] - dh * yp[i];
    }
    ADVec3 vph, vmh;
    ADVec3Axpy(scale, alphaph, xph, yph, vph);
    ADVec3Axpy(scale, alphamh, xmh, ymh, vmh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsRealPart((vph.x[i] - vmh.x[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(vt.xp[i]) - ref));
    }
#endif
    printf(
        "  [ADVec3Axpy] hforward vs directional-derivative-of-primal: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) hreverse (.xh/.valueh) vs hand-computed formula
  // --------------------
  {
    TacsScalar vd[3], vh_seed[3];
    for (int i = 0; i < 3; i++) {
      vd[i] = randReal();
      vh_seed[i] = randReal();
    }

    ADScalar alphat = alpha;
    alphat.valuep = alphap;
    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADVec3 vt;
    ADVec3Axpy op(scale, alphat, xt, yt, vt);
    op.hforward();
    for (int i = 0; i < 3; i++) {
      vt.xd[i] = vd[i];
      vt.xh[i] = vh_seed[i];
    }
    op.reverse();
    op.hreverse();

    // Expected:
    //   alpha.valueh = scale*(x.x . v.xh) + scale*(xp . v.xd)
    //   x.xh = scale*(alpha.value*v.xh + alpha.valuep*v.xd)
    //   y.xh = v.xh
    double s = TacsRealPart(scale);
    double expectedAlphaH = 0.0, expectedXh[3] = {0, 0, 0},
           expectedYh[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
      expectedAlphaH += s * TacsRealPart(x.x[i]) * TacsRealPart(vh_seed[i]);
      expectedAlphaH += s * TacsRealPart(xp[i]) * TacsRealPart(vd[i]);
      expectedXh[i] =
          s * (TacsRealPart(alpha.value) * TacsRealPart(vh_seed[i]) +
               TacsRealPart(alphap) * TacsRealPart(vd[i]));
      expectedYh[i] = TacsRealPart(vh_seed[i]);
    }

    double maxerr = fabs(TacsRealPart(alphat.valueh) - expectedAlphaH);
    for (int i = 0; i < 3; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - expectedXh[i]));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xh[i]) - expectedYh[i]));
    }
    printf(
        "  [ADVec3Axpy] hreverse (.xh/.valueh) vs hand-computed formula: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (3) First-order forward()/reverse() call sites are unaffected --------
  {
    ADScalar alpha0(0.5, 0.2);
    ADVec3 x0, y0, v0;
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      x0.xd[i] = randReal();
      y0.x[i] = randReal();
      y0.xd[i] = randReal();
    }
    ADVec3Axpy op0(scale, alpha0, x0, y0, v0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
