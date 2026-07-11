// Op-level second-order verification for A2D::ADMat3x3FromThreeADVec3
// (feature-beam-element-methods, SPEC.md sec 1.2/sec 6.6). Standalone,
// header-only.
//
// Build (real mode):
//   mpicxx -std=c++11 -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3fromthreeadvec3_second_order.cpp -o t17_real && ./t17_real
// Build (complex mode):
//   mpicxx -std=c++11 -DTACS_USE_COMPLEX -I<repo>/src -I<repo>/src/elements/a2d \
//     test_admat3x3fromthreeadvec3_second_order.cpp -o t17_complex && ./t17_complex
//
// C's columns are literally x/y/z's components -- a trivial linear
// assignment (not even a product of two active quantities), so this op's
// hforward/hreverse carry no curvature at all: same verification shape as
// Tasks 1.5/1.6 (directional-derivative check, bit-exact zero-self-
// Hessian, hand-computed hreverse formula).

#include "TACSObject.h"
#include "a2d.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace A2D;

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

int main() {
  srand(3579);

  ADVec3 x, y, z;
  for (int i = 0; i < 3; i++) {
    x.x[i] = randReal();
    y.x[i] = randReal();
    z.x[i] = randReal();
  }
  TacsScalar xp[3], yp[3], zp[3];
  for (int i = 0; i < 3; i++) {
    xp[i] = randReal();
    yp[i] = randReal();
    zp[i] = randReal();
  }

  int fail = 0;

  // --- (1) hforward vs directional derivative of the primal ------------------
  {
    ADVec3 xt = x, yt = y, zt = z;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
      zt.xp[i] = zp[i];
    }
    ADMat3x3 Ct;
    ADMat3x3FromThreeADVec3 op(xt, yt, zt, Ct);
    op.hforward();

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    ADVec3 xh, yh, zh;
    for (int i = 0; i < 3; i++) {
      xh.x[i] = x.x[i] + TacsScalar(0.0, dh) * xp[i];
      yh.x[i] = y.x[i] + TacsScalar(0.0, dh) * yp[i];
      zh.x[i] = z.x[i] + TacsScalar(0.0, dh) * zp[i];
    }
    ADMat3x3 Ch;
    ADMat3x3FromThreeADVec3(xh, yh, zh, Ch);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsImagPart(Ch.A[i]) / dh;
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#else
    double dh = 1e-6;
    ADVec3 xph, yph, zph, xmh, ymh, zmh;
    for (int i = 0; i < 3; i++) {
      xph.x[i] = x.x[i] + dh * xp[i];
      yph.x[i] = y.x[i] + dh * yp[i];
      zph.x[i] = z.x[i] + dh * zp[i];
      xmh.x[i] = x.x[i] - dh * xp[i];
      ymh.x[i] = y.x[i] - dh * yp[i];
      zmh.x[i] = z.x[i] - dh * zp[i];
    }
    ADMat3x3 Cph, Cmh;
    ADMat3x3FromThreeADVec3(xph, yph, zph, Cph);
    ADMat3x3FromThreeADVec3(xmh, ymh, zmh, Cmh);
    double maxerr = 0.0;
    for (int i = 0; i < 9; i++) {
      double ref = TacsRealPart((Cph.A[i] - Cmh.A[i]) / (2.0 * dh));
      maxerr = std::max(maxerr, fabs(TacsRealPart(Ct.Ap[i]) - ref));
    }
#endif
    printf(
        "  [ADMat3x3FromThreeADVec3] hforward vs "
        "directional-derivative-of-primal: %.3e\n",
        maxerr);
    if (maxerr > 1e-6) fail = 1;
  }

  // --- (2) zero-self-Hessian: C.Ah = 0, seeds arbitrary -> x/y/z.xh == 0 -----
  {
    ADVec3 xt = x, yt = y, zt = z;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
      zt.xp[i] = zp[i];
    }
    ADMat3x3 Ct;
    ADMat3x3FromThreeADVec3 op(xt, yt, zt, Ct);
    op.hforward();
    // Ct.Ah left at its zero-initialized default.
    op.hreverse();

    double maxabs = 0.0;
    for (int i = 0; i < 3; i++) {
      maxabs = std::max(maxabs, fabs(TacsRealPart(xt.xh[i])));
      maxabs = std::max(maxabs, fabs(TacsRealPart(yt.xh[i])));
      maxabs = std::max(maxabs, fabs(TacsRealPart(zt.xh[i])));
    }
    printf(
        "  [ADMat3x3FromThreeADVec3] zero-self-Hessian check (expect "
        "exactly 0.0): %.3e\n",
        maxabs);
    if (maxabs != 0.0) fail = 1;
  }

  // --- (3) hreverse vs hand-computed formula for a nonzero downstream seed --
  {
    TacsScalar Ch[9];
    for (int i = 0; i < 9; i++) Ch[i] = randReal();

    ADVec3 xt = x, yt = y, zt = z;
    for (int i = 0; i < 3; i++) {
      xt.xp[i] = xp[i];
      yt.xp[i] = yp[i];
      zt.xp[i] = zp[i];
    }
    ADMat3x3 Ct;
    ADMat3x3FromThreeADVec3 op(xt, yt, zt, Ct);
    op.hforward();
    for (int i = 0; i < 9; i++) Ct.Ah[i] = Ch[i];
    op.hreverse();

    double expected_x[3] = {TacsRealPart(Ch[0]), TacsRealPart(Ch[3]),
                            TacsRealPart(Ch[6])};
    double expected_y[3] = {TacsRealPart(Ch[1]), TacsRealPart(Ch[4]),
                            TacsRealPart(Ch[7])};
    double expected_z[3] = {TacsRealPart(Ch[2]), TacsRealPart(Ch[5]),
                            TacsRealPart(Ch[8])};
    double maxerr = 0.0;
    for (int i = 0; i < 3; i++) {
      maxerr = std::max(maxerr, fabs(TacsRealPart(xt.xh[i]) - expected_x[i]));
      maxerr = std::max(maxerr, fabs(TacsRealPart(yt.xh[i]) - expected_y[i]));
      maxerr = std::max(maxerr, fabs(TacsRealPart(zt.xh[i]) - expected_z[i]));
    }
    printf(
        "  [ADMat3x3FromThreeADVec3] hreverse vs hand-computed formula: "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-10) fail = 1;
  }

  // --- (4) First-order forward()/reverse() call sites are unaffected --------
  {
    ADVec3 x0, y0, z0;
    for (int i = 0; i < 3; i++) {
      x0.x[i] = randReal();
      y0.x[i] = randReal();
      z0.x[i] = randReal();
      x0.xd[i] = randReal();
      y0.xd[i] = randReal();
      z0.xd[i] = randReal();
    }
    ADMat3x3 C0;
    ADMat3x3FromThreeADVec3 op0(x0, y0, z0, C0);
    op0.forward();
    op0.reverse();
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
