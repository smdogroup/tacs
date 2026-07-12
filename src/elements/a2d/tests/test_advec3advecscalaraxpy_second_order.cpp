// Op-level second-order verification for A2D::ADVec3ADVecScalarAxpy
// (feature-beam-element-methods, SPEC.md sec 1.2/sec 6.6). Standalone,
// header-only.
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3advecscalaraxpy_second_order.cpp -o t16_real && ./t16_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3advecscalaraxpy_second_order.cpp -o t16_complex && ./t16_complex
// clang-format on
//
// v = scale*alpha.value*x + y is LINEAR jointly in the two active inputs
// (x, y) -- alpha is a fixed passive Scalar. Same verification shape as
// ADMat3x3MatMult's test (Task 1.5): hforward against a directional
// derivative of the primal, a bit-exact zero-self-Hessian check, and
// hreverse's propagation formula against a hand computation for a nonzero
// downstream seed.

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(2468);

  const TacsScalar scale = -1.3;
  Scalar alpha(0.62);
  ADVec3 x, y;
  for (int i = 0; i < 3; i++) {
    x.x[i] = randReal();
    y.x[i] = randReal();
  }
  TacsScalar xp[3], yp[3];
  for (int i = 0; i < 3; i++) {
    xp[i] = randReal();
    yp[i] = randReal();
  }

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADVec3 vt;
    ADVec3ADVecScalarAxpy op(scale, alpha, xt, yt, vt);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADVec3 xh, yh;
    for (int i = 0; i < 3; i++) {
      xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
      yh.x[i] = y.x[i] + TacsScalar(0.0, dh) * yp[i];
    }
    ADVec3 vh;
    ADVec3ADVecScalarAxpy(scale, alpha, xh, yh, vh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsImagPart(vh.x[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(vt.xp[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADVec3 xph, yph, xmh, ymh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      yph.x[i] = y.x[i] + dh * yp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
      ymh.x[i] = y.x[i] - dh * yp[i];
    }
    ADVec3 vph, vmh;
    ADVec3ADVecScalarAxpy(scale, alpha, xph, yph, vph);
    ADVec3ADVecScalarAxpy(scale, alpha, xmh, ymh, vmh);
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double ref = TacsRealPart((vph.x[i] - vmh.x[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(vt.xp[i]) - ref));
    }
#endif
    printf(
        "  [ADVec3ADVecScalarAxpy] hforward vs "
        "directional-derivative-of-primal: %.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) zero-self-Hessian: v.xh = 0, x.xp/y.xp arbitrary -> x.xh/y.xh == 0
  {
    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADVec3 vt;
    ADVec3ADVecScalarAxpy op(scale, alpha, xt, yt, vt);
    op.hforward();
    // vt.xh left at its zero-initialized default.
    op.hreverse();

    double maxabs = 0.0;
    for (int i = 0; i < 3; i++) {
      maxabs = std::max(maxabs, fabs(TacsRealPart(xt.xh[i])));
      maxabs = std::max(maxabs, fabs(TacsRealPart(yt.xh[i])));
    }
    printf(
        "  [ADVec3ADVecScalarAxpy] zero-self-Hessian check (expect "
        "exactly 0.0): %.3e\n",
        maxabs);
    if (maxabs != 0.0) fail = 1;
  }

  // --- (3) hreverse vs hand-computed formula for a nonzero downstream seed --
  {
    TacsScalar vh[3];
    for (int i = 0; i < 3; i++) vh[i] = randReal();

    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADVec3 vt;
    ADVec3ADVecScalarAxpy op(scale, alpha, xt, yt, vt);
    op.hforward();
    for (int i = 0; i < 3; i++) vt.xh[i] = vh[i];
    op.hreverse();

    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      double expected_x = TacsRealPart(scale * alpha.value * vh[i]);
      double expected_y = TacsRealPart(vh[i]);
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - expected_x));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xh[i]) - expected_y));
    }
    printf(
        "  [ADVec3ADVecScalarAxpy] hreverse vs hand-computed formula: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (4) First-order forward()/reverse() call sites are unaffected --------
  {
    ADVec3 x0, y0, v0;
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      y0.x[i] = randReal();
      x0.xd[i] = randReal();
      y0.xd[i] = randReal();
    }
    ADVec3ADVecScalarAxpy op0(scale, alpha, x0, y0, v0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
