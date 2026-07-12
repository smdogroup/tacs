// Director-hook-level verification for TACSDirector.h's two new
// second-derivative hooks (feature-beam-element-methods, SPEC-phase-7.md
// sec 3.1/3.2/6.2): addDirectorHessianProduct and
// addDirectorHessianRefNormalSens, on TACSQuadraticRotation and
// TACSQuaternionRotation. Unlike src/elements/a2d/tests/*, this is NOT
// header-only (TACSDirector.h pulls in TACSElementVerification.h ->
// TACSElement.h -> the full TACS include graph) -- it links against the
// already-built libtacs.so, per SPEC-phase-7.md sec 6.2's own allowance
// ("not necessarily header-only/no-build, since TACSDirector.h is already
// included by the existing TACS build"). NOT part of
// src/elements/a2d/tests/run_tests.sh's glob (a different test class);
// build/run manually per the recipe below.
//
// clang-format off
// Build (real mode, from the repo root, after `python3 setup.py build_ext
// --inplace` has produced lib/libtacs.so):
//   mpicxx -std=c++11 -Isrc -Isrc/bpmat -Isrc/elements -Isrc/elements/a2d \
//     -Isrc/elements/dynamics -Isrc/elements/basis -Isrc/elements/shell \
//     -Isrc/constitutive -Isrc/functions -Isrc/io -Iextern/AMD/Include \
//     -Iextern/UFconfig -Iextern/metis/include \
//     src/elements/shell/tests/test_director_hessian_second_order.cpp \
//     -Llib -ltacs -Wl,-rpath,$(pwd)/lib -o director_hessian_test
//   ./director_hessian_test
// Build (complex mode): add -DTACS_USE_COMPLEX and link against a
// complex-mode libtacs.so instead (this repo builds real/complex into
// separate trees; see HANDOFF's own real/complex rebuild-loop notes).
// clang-format on
//
// Two checks per director class:
// (1) Regression-safe-refactor check (SPEC sec 3.1's own required FIRST
//     step): addDirectorHessianProduct, called with qa=qb=e_k/e_l over all
//     (k,l) basis-vector pairs and scattered into a 3x3 (resp. 4x4) block,
//     must reproduce addRotationMatJacobian's own existing i==j diagonal
//     block bit-for-bit (both use the SAME dC; the new hook reconstructs
//     dC internally as outer(dd, t)).
// (2) Complex-step check of addDirectorHessianRefNormalSens: perturbing t
//     by i*dh must match addDirectorHessianProduct's own value's imaginary
//     part / dh (since the Xpts-adjoint returned is exactly d(val)/dt).

#include <cmath>
#include <cstdio>

#include "TACSDirector.h"

static double randReal() { return 2.0 * (double(rand()) / RAND_MAX) - 1.0; }

// Reconstruct addRotationMatJacobian's own i==j diagonal block value at
// (row,col) via dC = outer(dd,t), so the SAME formula is exercised twice by
// two different code paths (existing addRotationMatJacobian and the new
// standalone hook), matching bit-for-bit is the actual regression claim.
static void buildDC(const TacsScalar dd[3], const TacsScalar t[3],
                     TacsScalar dC[9]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dC[3 * i + j] = dd[i] * t[j];
    }
  }
}

int main() {
  srand(74001);
  int fail = 0;
  const int vars_per_node = 6;
  const int offset = 3;
  const int num_nodes = 1;

  TacsScalar t[3], dd[3];
  for (int i = 0; i < 3; i++) {
    t[i] = randReal();
    dd[i] = randReal();
  }

  // =========================================================================
  // TACSQuadraticRotation
  // =========================================================================
  {
    TacsScalar dC[9];
    buildDC(dd, t, dC);
    // Existing diagonal-block formula (TACSDirector.h:806-818), reproduced
    // directly here as the independent reference (not calling
    // addRotationMatJacobian itself, which requires a full vars[]/mat[]
    // system context) -- this IS the exact formula the new hook must match.
    TacsScalar jacRef[9] = {0};
    jacRef[0] -= dC[4] + dC[8];
    jacRef[1] += 0.5 * (dC[1] + dC[3]);
    jacRef[2] += 0.5 * (dC[2] + dC[6]);
    jacRef[3] += 0.5 * (dC[1] + dC[3]);
    jacRef[4] -= dC[0] + dC[8];
    jacRef[5] += 0.5 * (dC[5] + dC[7]);
    jacRef[6] += 0.5 * (dC[2] + dC[6]);
    jacRef[7] += 0.5 * (dC[5] + dC[7]);
    jacRef[8] -= dC[0] + dC[4];

    TacsScalar vars[offset + 3] = {0};
    double maxerr = 0.0;
    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++) {
        TacsScalar qa[offset + 3] = {0}, qb[offset + 3] = {0};
        qa[offset + k] = 1.0;
        qb[offset + l] = 1.0;
        TacsScalar mat[1] = {0.0};
        TACSQuadraticRotation::addDirectorHessianProduct<vars_per_node, offset,
                                                          num_nodes>(t, dd, qa,
                                                                     qb, mat);
        maxerr =
            std::max(maxerr, fabs(TacsRealPart(mat[0]) -
                                   TacsRealPart(jacRef[3 * k + l])));
      }
    }
    printf(
        "  [TACSQuadraticRotation] addDirectorHessianProduct vs "
        "addRotationMatJacobian's own diagonal block (regression check): "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-13) fail = 1;

    // Complex-step check of addDirectorHessianRefNormalSens against
    // addDirectorHessianProduct's own value, perturbing t.
    TacsScalar qa[offset + 3] = {0}, qb[offset + 3] = {0};
    for (int i = 0; i < 3; i++) {
      qa[offset + i] = randReal();
      qb[offset + i] = randReal();
    }
    TacsScalar dt[3] = {0, 0, 0};
    TACSQuadraticRotation::addDirectorHessianRefNormalSens<vars_per_node,
                                                            offset, num_nodes>(
        t, dd, qa, qb, dt);

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    double maxerr2 = 0.0;
    for (int p = 0; p < 3; p++) {
      TacsScalar tph[3];
      for (int i = 0; i < 3; i++) tph[i] = t[i];
      tph[p] += TacsScalar(0.0, dh);
      TacsScalar matph[1] = {0.0};
      TACSQuadraticRotation::addDirectorHessianProduct<vars_per_node, offset,
                                                        num_nodes>(tph, dd, qa,
                                                                   qb, matph);
      double ref = TacsImagPart(matph[0]) / dh;
      maxerr2 = std::max(maxerr2, fabs(TacsRealPart(dt[p]) - ref));
    }
#else
    double dh = 1e-6;
    double maxerr2 = 0.0;
    for (int p = 0; p < 3; p++) {
      TacsScalar tph[3], tmh[3];
      for (int i = 0; i < 3; i++) {
        tph[i] = t[i];
        tmh[i] = t[i];
      }
      tph[p] += dh;
      tmh[p] -= dh;
      TacsScalar matph[1] = {0.0}, matmh[1] = {0.0};
      TACSQuadraticRotation::addDirectorHessianProduct<vars_per_node, offset,
                                                        num_nodes>(tph, dd, qa,
                                                                   qb, matph);
      TACSQuadraticRotation::addDirectorHessianProduct<vars_per_node, offset,
                                                        num_nodes>(tmh, dd, qa,
                                                                   qb, matmh);
      double ref =
          TacsRealPart((matph[0] - matmh[0]) / (2.0 * dh));
      maxerr2 = std::max(maxerr2, fabs(TacsRealPart(dt[p]) - ref));
    }
#endif
    printf(
        "  [TACSQuadraticRotation] addDirectorHessianRefNormalSens vs "
        "directional-derivative-of-addDirectorHessianProduct w.r.t. t: "
        "%.3e\n",
        maxerr2);
    if (maxerr2 > 1e-6) fail = 1;
  }

  // =========================================================================
  // TACSQuaternionRotation (5 vars: q[0..3] quaternion + lambda; offset/mat
  // size is 4 here, matching NUM_PARAMETERS)
  // =========================================================================
  {
    TacsScalar dC[9];
    buildDC(dd, t, dC);
    TacsScalar jacRef[16] = {0};
    jacRef[1] += 2.0 * (dC[5] - dC[7]);
    jacRef[2] += 2.0 * (dC[6] - dC[2]);
    jacRef[3] += 2.0 * (dC[1] - dC[3]);
    jacRef[4] += 2.0 * (dC[5] - dC[7]);
    jacRef[5] -= 4.0 * (dC[4] + dC[8]);
    jacRef[6] += 2.0 * (dC[1] + dC[3]);
    jacRef[7] += 2.0 * (dC[2] + dC[6]);
    jacRef[8] += 2.0 * (dC[6] - dC[2]);
    jacRef[9] += 2.0 * (dC[1] + dC[3]);
    jacRef[10] -= 4.0 * (dC[0] + dC[8]);
    jacRef[11] += 2.0 * (dC[5] + dC[7]);
    jacRef[12] += 2.0 * (dC[1] - dC[3]);
    jacRef[13] += 2.0 * (dC[2] + dC[6]);
    jacRef[14] += 2.0 * (dC[5] + dC[7]);
    jacRef[15] -= 4.0 * (dC[0] + dC[4]);

    const int qoffset = 3;  // quaternion offset used for this probe
    double maxerr = 0.0;
    for (int k = 0; k < 4; k++) {
      for (int l = 0; l < 4; l++) {
        TacsScalar qa[qoffset + 4] = {0}, qb[qoffset + 4] = {0};
        qa[qoffset + k] = 1.0;
        qb[qoffset + l] = 1.0;
        TacsScalar mat[1] = {0.0};
        TACSQuaternionRotation::addDirectorHessianProduct<vars_per_node,
                                                           qoffset, num_nodes>(
            t, dd, qa, qb, mat);
        maxerr = std::max(
            maxerr,
            fabs(TacsRealPart(mat[0]) - TacsRealPart(jacRef[4 * k + l])));
      }
    }
    printf(
        "  [TACSQuaternionRotation] addDirectorHessianProduct vs "
        "addRotationMatJacobian's own diagonal block (regression check): "
        "%.3e\n",
        maxerr);
    if (maxerr > 1e-13) fail = 1;

    TacsScalar qa[qoffset + 4] = {0}, qb[qoffset + 4] = {0};
    for (int i = 0; i < 4; i++) {
      qa[qoffset + i] = randReal();
      qb[qoffset + i] = randReal();
    }
    TacsScalar dt[3] = {0, 0, 0};
    TACSQuaternionRotation::addDirectorHessianRefNormalSens<vars_per_node,
                                                             qoffset,
                                                             num_nodes>(
        t, dd, qa, qb, dt);

#ifdef TACS_USE_COMPLEX
    double dh = 1e-30;
    double maxerr2 = 0.0;
    for (int p = 0; p < 3; p++) {
      TacsScalar tph[3];
      for (int i = 0; i < 3; i++) tph[i] = t[i];
      tph[p] += TacsScalar(0.0, dh);
      TacsScalar matph[1] = {0.0};
      TACSQuaternionRotation::addDirectorHessianProduct<vars_per_node,
                                                         qoffset, num_nodes>(
          tph, dd, qa, qb, matph);
      double ref = TacsImagPart(matph[0]) / dh;
      maxerr2 = std::max(maxerr2, fabs(TacsRealPart(dt[p]) - ref));
    }
#else
    double dh = 1e-6;
    double maxerr2 = 0.0;
    for (int p = 0; p < 3; p++) {
      TacsScalar tph[3], tmh[3];
      for (int i = 0; i < 3; i++) {
        tph[i] = t[i];
        tmh[i] = t[i];
      }
      tph[p] += dh;
      tmh[p] -= dh;
      TacsScalar matph[1] = {0.0}, matmh[1] = {0.0};
      TACSQuaternionRotation::addDirectorHessianProduct<vars_per_node,
                                                         qoffset, num_nodes>(
          tph, dd, qa, qb, matph);
      TACSQuaternionRotation::addDirectorHessianProduct<vars_per_node,
                                                         qoffset, num_nodes>(
          tmh, dd, qa, qb, matmh);
      double ref = TacsRealPart((matph[0] - matmh[0]) / (2.0 * dh));
      maxerr2 = std::max(maxerr2, fabs(TacsRealPart(dt[p]) - ref));
    }
#endif
    printf(
        "  [TACSQuaternionRotation] addDirectorHessianRefNormalSens vs "
        "directional-derivative-of-addDirectorHessianProduct w.r.t. t: "
        "%.3e\n",
        maxerr2);
    if (maxerr2 > 1e-6) fail = 1;
  }

  printf(fail ? "FAIL\n" : "PASS\n");
  return fail;
}
