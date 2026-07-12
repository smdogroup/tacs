// Op-level second-order verification for A2D::ADVec3Dot (feature-beam-
// element-methods, SPEC.md sec 1.2/sec 6.6). Standalone, header-only:
// compiles directly against src/elements/a2d/*.h + TACSObject.h, no TACS
// build/libtacs linkage required (per VALIDATION.md's E6 experiment).
//
// NOTE: this file includes "TACSObject.h" before "a2d.h" to work around a
// header-ordering fragility in a2dmatops.h/a2dmatcore.h (they use
// TacsScalar before a2dobjs.h -- the file that #includes TACSObject.h -- is
// reached); this is a standalone-compilation-only issue, not a defect this
// feature is responsible for fixing (see PLAN.md's Phase 1 preamble).
//
// clang-format off
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3dot_second_order.cpp -o t14_real && ./t14_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_advec3dot_second_order.cpp -o t14_complex && ./t14_complex
// clang-format on
//
// Verification strategy (SPEC.md sec 6.6, ported from VALIDATION.md's E6):
//   1. hforward is checked against a directional derivative of the PRIMAL
//      relation alpha.value = scale*(x.x . y.x) along a random seed
//      direction (xp, yp): complex-step in the complex build (exact to
//      round-off), central difference in the real build (exact here too,
//      since alpha.value is an exact bilinear form in (x, y) -- the central-
//      difference formula has zero truncation error for a quadratic).
//   2. hreverse is checked against a central difference of the EXACT
//      reverse() adjoint (fixed alpha.valued = 1) as the primal inputs are
//      perturbed along the seed direction -- this is VALIDATION.md E6's
//      "central-difference-of-exact-adjoint" technique, which avoids
//      needing a nested-complex-step type for a genuine second derivative.

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "TACSObject.h"
#include "a2d.h"

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(1234);

  const TacsScalar scale = 1.7;
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

  // --- (1) hforward vs directional derivative of the primal -----------------
  {
    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADScalar alpha;
    ADVec3Dot op(scale, xt, yt, alpha);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADVec3 xh, yh;
    for (int i = 0; i < 3; i++) {
      xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
      yh.x[i] = y.x[i] + TacsScalar(0.0, dh) * yp[i];
    }
    ADScalar alphah;
    ADVec3Dot(scale, xh, yh, alphah);
    double ref = TacsImagPart(alphah.value) / dh;
#else
    double dh = 1e-6;
    ADVec3 xph, yph, xmh, ymh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      yph.x[i] = y.x[i] + dh * yp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
      ymh.x[i] = y.x[i] - dh * yp[i];
    }
    ADScalar alphap, alpham;
    ADVec3Dot(scale, xph, yph, alphap);
    ADVec3Dot(scale, xmh, ymh, alpham);
    double ref = TacsRealPart((alphap.value - alpham.value) / (2.0 * dh));
#endif
    double err = fabs(TacsRealPart(alpha.valuep) - ref);
    printf("  [ADVec3Dot] hforward vs directional-derivative-of-primal: %.3e\n",
           err);
    if (err > 1e-6) fail = 1;
  }

  // --- (2) hreverse vs central-diff-of-exact-adjoint -------------------------
  {
    auto adjointAt = [&](const TacsScalar xv[3], const TacsScalar yv[3],
                         TacsScalar gradx[3], TacsScalar grady[3]) {
      ADVec3 xt(xv), yt(yv);
      ADScalar at;
      ADVec3Dot op(scale, xt, yt, at);
      at.valued = 1.0;
      op.reverse();
      for (int i = 0; i < 3; i++) {
        gradx[i] = xt.xd[i];
        grady[i] = yt.xd[i];
      }
    };

    double dh = 1e-5;
    double Href[6];  // d/dt [adjoint(x+t*xp, y+t*yp)] at t=0, ordered (x,y)
    {
      TacsScalar xph[3], yph[3], xmh[3], ymh[3];
      for (int i = 0; i < 3; i++) {
        xph[i] = x.x[i] + dh * xp[i];
        yph[i] = y.x[i] + dh * yp[i];
        xmh[i] = x.x[i] - dh * xp[i];
        ymh[i] = y.x[i] - dh * yp[i];
      }
      TacsScalar gxp[3], gyp[3], gxm[3], gym[3];
      adjointAt(xph, yph, gxp, gyp);
      adjointAt(xmh, ymh, gxm, gym);
      for (int i = 0; i < 3; i++) {
        Href[i] = TacsRealPart((gxp[i] - gxm[i]) / (2.0 * dh));
        Href[3 + i] = TacsRealPart((gyp[i] - gym[i]) / (2.0 * dh));
      }
    }

    ADVec3 xt = x, yt = y;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
    }
    ADScalar at;
    ADVec3Dot op(scale, xt, yt, at);
    at.valued = 1.0;
    at.valueh = 0.0;
    op.hreverse();

    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - Href[i]));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xh[i]) - Href[3 + i]));
    }
    printf("  [ADVec3Dot] hreverse vs central-diff-of-exact-adjoint: %.3e\n",
           maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (3) First-order forward()/reverse() call sites are unaffected --------
  {
    ADVec3 x0, y0;
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      y0.x[i] = randReal();
      x0.xd[i] = randReal();
    }
    ADScalar a0;
    ADVec3Dot op0(scale, x0, y0, a0);
    op0.forward();
    op0.reverse();
    // No assertion beyond "still compiles and runs" -- regression guard for
    // the existing first-order call sites in TACSBeamElement.h.
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
