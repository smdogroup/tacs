/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSTensorProductBasisImpl.h"

void TACSInterpAllTensor3DInterp5(const int m, const double N[],
                                  const double Nx[], const TacsScalar values[],
                                  TacsScalar out[]) {
  memset(out, 0, 4 * m * 125 * sizeof(TacsScalar));

  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    const TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar *g = &out[m];
        for (int p = 0; p < m; p++) {
          g[0] += n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                         n1x[3] * v[3 * m] + n1x[4] * v[4 * m]);
          TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
                           n1[3] * v[3 * m] + n1[4] * v[4 * m]);
          out[p] += n23 * t1;
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }

    out += 4 * m;
  }
}

void TACSInterpAllTensor3DInterp5VarsPerNode1(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 125 * sizeof(TacsScalar));

  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    const TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[1] += n23 * (n1x[0] * v[0] + n1x[1] * v[1] + n1x[2] * v[2] +
                         n1x[3] * v[3] + n1x[4] * v[4]);
        TacsScalar t1 = (n1[0] * v[0] + n1[1] * v[1] + n1[2] * v[2] +
                         n1[3] * v[3] + n1[4] * v[4]);
        out[0] += n23 * t1;
        out[2] += n2x3 * t1;
        out[3] += n23x * t1;
        v += 5;
      }
    }
    out += 4;
  }
}

void TACSInterpAllTensor3DInterp5VarsPerNode3(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 3 * 125 * sizeof(TacsScalar));

  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    const TacsScalar *v = values;
    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar t1;
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[3] += n23 * (n1x[0] * v[0] + n1x[1] * v[3] + n1x[2] * v[6] +
                         n1x[3] * v[9] + n1x[4] * v[12]);
        t1 = (n1[0] * v[0] + n1[1] * v[3] + n1[2] * v[6] + n1[3] * v[9] +
              n1[4] * v[12]);
        out[0] += n23 * t1;
        out[4] += n2x3 * t1;
        out[5] += n23x * t1;

        out[6] += n23 * (n1x[0] * v[1] + n1x[1] * v[4] + n1x[2] * v[7] +
                         n1x[3] * v[10] + n1x[4] * v[13]);
        t1 = (n1[0] * v[1] + n1[1] * v[4] + n1[2] * v[7] + n1[3] * v[10] +
              n1[4] * v[13]);
        out[1] += n23 * t1;
        out[7] += n2x3 * t1;
        out[8] += n23x * t1;

        out[9] += n23 * (n1x[0] * v[2] + n1x[1] * v[5] + n1x[2] * v[8] +
                         n1x[3] * v[11] + n1x[4] * v[14]);
        t1 = (n1[0] * v[2] + n1[1] * v[5] + n1[2] * v[8] + n1[3] * v[11] +
              n1[4] * v[14]);
        out[2] += n23 * t1;
        out[10] += n2x3 * t1;
        out[11] += n23x * t1;

        v += 15;
      }
    }
    out += 12;
  }
}

void TACSInterpAllTensor3DInterp5VarsPerNode4(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 4 * 125 * sizeof(TacsScalar));

  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    const TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        TacsScalar t1;
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[4] += n23 * (n1x[0] * v[0] + n1x[1] * v[4] + n1x[2] * v[8] +
                         n1x[3] * v[12] + n1x[4] * v[16]);
        t1 = (n1[0] * v[0] + n1[1] * v[4] + n1[2] * v[8] + n1[3] * v[12] +
              n1[4] * v[16]);
        out[0] += n23 * t1;
        out[5] += n2x3 * t1;
        out[6] += n23x * t1;

        out[7] += n23 * (n1x[0] * v[1] + n1x[1] * v[5] + n1x[2] * v[9] +
                         n1x[3] * v[13] + n1x[4] * v[17]);
        t1 = (n1[0] * v[1] + n1[1] * v[5] + n1[2] * v[9] + n1[3] * v[13] +
              n1[4] * v[17]);
        out[1] += n23 * t1;
        out[8] += n2x3 * t1;
        out[9] += n23x * t1;

        out[10] += n23 * (n1x[0] * v[2] + n1x[1] * v[6] + n1x[2] * v[10] +
                          n1x[3] * v[14] + n1x[4] * v[18]);
        t1 = (n1[0] * v[2] + n1[1] * v[6] + n1[2] * v[10] + n1[3] * v[14] +
              n1[4] * v[18]);
        out[2] += n23 * t1;
        out[11] += n2x3 * t1;
        out[12] += n23x * t1;

        out[13] += n23 * (n1x[0] * v[3] + n1x[1] * v[7] + n1x[2] * v[11] +
                          n1x[3] * v[15] + n1x[4] * v[19]);
        t1 = (n1[0] * v[3] + n1[1] * v[7] + n1[2] * v[11] + n1[3] * v[15] +
              n1[4] * v[19]);
        out[3] += n23 * t1;
        out[14] += n2x3 * t1;
        out[15] += n23x * t1;

        v += 20;
      }
    }

    out += 16;
  }
}

void TacsAddAllTransTensor3DInterp5(const int m, const double N[],
                                    const double Nx[], const TacsScalar in[],
                                    TacsScalar values[]) {
  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        const TacsScalar *g = &in[m];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * in[p] + n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          g += 3;
          v++;
        }
        v += 4 * m;
      }
    }
    in += 4 * m;
  }
}

void TacsAddAllTransTensor3DInterp5VarsPerNode1(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a = n23 * in[0] + n2x3 * in[2] + n23x * in[3];
        TacsScalar b = n23 * in[1];
        v[0] += a * n1[0] + b * n1x[0];
        v[1] += a * n1[1] + b * n1x[1];
        v[2] += a * n1[2] + b * n1x[2];
        v[3] += a * n1[3] + b * n1x[3];
        v[4] += a * n1[4] + b * n1x[4];
        v += 5;
      }
    }
    in += 4;
  }
}

void TacsAddAllTransTensor3DInterp5VarsPerNode3(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    TacsScalar *v = values;
    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a1 = n23 * in[0] + n2x3 * in[4] + n23x * in[5];
        TacsScalar b1 = n23 * in[3];

        TacsScalar a2 = n23 * in[1] + n2x3 * in[7] + n23x * in[8];
        TacsScalar b2 = n23 * in[6];

        TacsScalar a3 = n23 * in[2] + n2x3 * in[10] + n23x * in[11];
        TacsScalar b3 = n23 * in[9];

        v[0] += a1 * n1[0] + b1 * n1x[0];
        v[1] += a2 * n1[0] + b2 * n1x[0];
        v[2] += a3 * n1[0] + b3 * n1x[0];

        v[3] += a1 * n1[1] + b1 * n1x[1];
        v[4] += a2 * n1[1] + b2 * n1x[1];
        v[5] += a3 * n1[1] + b3 * n1x[1];

        v[6] += a1 * n1[2] + b1 * n1x[2];
        v[7] += a2 * n1[2] + b2 * n1x[2];
        v[8] += a3 * n1[2] + b3 * n1x[2];

        v[9] += a1 * n1[3] + b1 * n1x[3];
        v[10] += a2 * n1[3] + b2 * n1x[3];
        v[11] += a3 * n1[3] + b3 * n1x[3];

        v[12] += a1 * n1[4] + b1 * n1x[4];
        v[13] += a2 * n1[4] + b2 * n1x[4];
        v[14] += a3 * n1[4] + b3 * n1x[4];

        v += 15;
      }
    }
    in += 12;
  }
}

void TacsAddAllTransTensor3DInterp5VarsPerNode4(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 125; n++) {
    const double *n1 = &N[5 * (n % 5)];
    const double *n2 = &N[5 * ((n % 25) / 5)];
    const double *n3 = &N[5 * (n / 25)];
    const double *n1x = &Nx[5 * (n % 5)];
    const double *n2x = &Nx[5 * ((n % 25) / 5)];
    const double *n3x = &Nx[5 * (n / 25)];

    TacsScalar *v = values;

    for (int k = 0; k < 5; k++) {
      for (int j = 0; j < 5; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a1 = n23 * in[0] + n2x3 * in[5] + n23x * in[6];
        TacsScalar b1 = n23 * in[4];

        TacsScalar a2 = n23 * in[1] + n2x3 * in[8] + n23x * in[9];
        TacsScalar b2 = n23 * in[7];

        TacsScalar a3 = n23 * in[2] + n2x3 * in[11] + n23x * in[12];
        TacsScalar b3 = n23 * in[10];

        TacsScalar a4 = n23 * in[3] + n2x3 * in[14] + n23x * in[15];
        TacsScalar b4 = n23 * in[13];

        v[0] += a1 * n1[0] + b1 * n1x[0];
        v[1] += a2 * n1[0] + b2 * n1x[0];
        v[2] += a3 * n1[0] + b3 * n1x[0];
        v[3] += a4 * n1[0] + b4 * n1x[0];

        v[4] += a1 * n1[1] + b1 * n1x[1];
        v[5] += a2 * n1[1] + b2 * n1x[1];
        v[6] += a3 * n1[1] + b3 * n1x[1];
        v[7] += a4 * n1[1] + b4 * n1x[1];

        v[8] += a1 * n1[2] + b1 * n1x[2];
        v[9] += a2 * n1[2] + b2 * n1x[2];
        v[10] += a3 * n1[2] + b3 * n1x[2];
        v[11] += a4 * n1[2] + b4 * n1x[2];

        v[12] += a1 * n1[3] + b1 * n1x[3];
        v[13] += a2 * n1[3] + b2 * n1x[3];
        v[14] += a3 * n1[3] + b3 * n1x[3];
        v[15] += a4 * n1[3] + b4 * n1x[3];

        v[16] += a1 * n1[4] + b1 * n1x[4];
        v[17] += a2 * n1[4] + b2 * n1x[4];
        v[18] += a3 * n1[4] + b3 * n1x[4];
        v[19] += a4 * n1[4] + b4 * n1x[4];

        v += 20;
      }
    }
    in += 16;
  }
}

void TACSInterpAllTensor3DInterp6(const int m, const double N[],
                                  const double Nx[], const TacsScalar values[],
                                  TacsScalar out[]) {
  memset(out, 0, 4 * m * 216 * sizeof(TacsScalar));
  for (int n = 0; n < 216; n++) {
    const double *n1 = &N[6 * (n % 6)];
    const double *n2 = &N[6 * ((n % 36) / 6)];
    const double *n3 = &N[6 * (n / 36)];
    const double *n1x = &Nx[6 * (n % 6)];
    const double *n2x = &Nx[6 * ((n % 36) / 6)];
    const double *n3x = &Nx[6 * (n / 36)];

    const TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar *g = &out[m];
        for (int p = 0; p < m; p++) {
          g[0] +=
              n23 * (n1x[0] * v[0] + n1x[1] * v[m] + n1x[2] * v[2 * m] +
                     n1x[3] * v[3 * m] + n1x[4] * v[4 * m] + n1x[5] * v[5 * m]);
          TacsScalar t1 =
              (n1[0] * v[0] + n1[1] * v[m] + n1[2] * v[2 * m] +
               n1[3] * v[3 * m] + n1[4] * v[4 * m] + n1[5] * v[5 * m]);
          out[p] += n23 * t1;
          g[1] += n2x3 * t1;
          g[2] += n23x * t1;
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }

    out += 4 * m;
  }
}

void TACSInterpAllTensor3DInterp6VarsPerNode1(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 216 * sizeof(TacsScalar));

  for (int n = 0; n < 36; n++) {
    const double *n11 = &N[0];
    const double *n12 = &N[6];
    const double *n13 = &N[12];
    const double *n14 = &N[18];
    const double *n15 = &N[24];
    const double *n16 = &N[30];

    const double *n1x1 = &Nx[0];
    const double *n1x2 = &Nx[6];
    const double *n1x3 = &Nx[12];
    const double *n1x4 = &Nx[18];
    const double *n1x5 = &Nx[24];
    const double *n1x6 = &Nx[30];

    const double *n2 = &N[6 * (n % 6)];
    const double *n3 = &N[6 * (n / 6)];
    const double *n2x = &Nx[6 * (n % 6)];
    const double *n3x = &Nx[6 * (n / 6)];

    const TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[1] += n23 * (n1x1[0] * v[0] + n1x1[1] * v[1] + n1x1[2] * v[2] +
                         n1x1[3] * v[3] + n1x1[4] * v[4] + n1x1[5] * v[5]);
        TacsScalar t1 = (n11[0] * v[0] + n11[1] * v[1] + n11[2] * v[2] +
                         n11[3] * v[3] + n11[4] * v[4] + n11[5] * v[5]);
        out[0] += n23 * t1;
        out[2] += n2x3 * t1;
        out[3] += n23x * t1;

        out[5] += n23 * (n1x2[0] * v[0] + n1x2[1] * v[1] + n1x2[2] * v[2] +
                         n1x2[3] * v[3] + n1x2[4] * v[4] + n1x2[5] * v[5]);
        TacsScalar t2 = (n12[0] * v[0] + n12[1] * v[1] + n12[2] * v[2] +
                         n12[3] * v[3] + n12[4] * v[4] + n12[5] * v[5]);
        out[4] += n23 * t2;
        out[6] += n2x3 * t2;
        out[7] += n23x * t2;

        out[9] += n23 * (n1x3[0] * v[0] + n1x3[1] * v[1] + n1x3[2] * v[2] +
                         n1x3[3] * v[3] + n1x3[4] * v[4] + n1x3[5] * v[5]);
        TacsScalar t3 = (n13[0] * v[0] + n13[1] * v[1] + n13[2] * v[2] +
                         n13[3] * v[3] + n13[4] * v[4] + n13[5] * v[5]);
        out[8] += n23 * t3;
        out[10] += n2x3 * t3;
        out[11] += n23x * t3;

        out[13] += n23 * (n1x4[0] * v[0] + n1x4[1] * v[1] + n1x4[2] * v[2] +
                          n1x4[3] * v[3] + n1x4[4] * v[4] + n1x4[5] * v[5]);
        TacsScalar t4 = (n14[0] * v[0] + n14[1] * v[1] + n14[2] * v[2] +
                         n14[3] * v[3] + n14[4] * v[4] + n14[5] * v[5]);
        out[12] += n23 * t4;
        out[14] += n2x3 * t4;
        out[15] += n23x * t4;

        out[17] += n23 * (n1x5[0] * v[0] + n1x5[1] * v[1] + n1x5[2] * v[2] +
                          n1x5[3] * v[3] + n1x5[4] * v[4] + n1x5[5] * v[5]);
        TacsScalar t5 = (n15[0] * v[0] + n15[1] * v[1] + n15[2] * v[2] +
                         n15[3] * v[3] + n15[4] * v[4] + n15[5] * v[5]);
        out[16] += n23 * t5;
        out[18] += n2x3 * t5;
        out[19] += n23x * t5;

        out[21] += n23 * (n1x6[0] * v[0] + n1x6[1] * v[1] + n1x6[2] * v[2] +
                          n1x6[3] * v[3] + n1x6[4] * v[4] + n1x6[5] * v[5]);
        TacsScalar t6 = (n16[0] * v[0] + n16[1] * v[1] + n16[2] * v[2] +
                         n16[3] * v[3] + n16[4] * v[4] + n16[5] * v[5]);
        out[20] += n23 * t6;
        out[22] += n2x3 * t6;
        out[23] += n23x * t6;

        v += 6;
      }
    }
    out += 24;
  }
}

void TACSInterpAllTensor3DInterp6VarsPerNode3(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 3 * 216 * sizeof(TacsScalar));

  for (int n = 0; n < 216; n++) {
    const double *n1 = &N[6 * (n % 6)];
    const double *n2 = &N[6 * ((n % 36) / 6)];
    const double *n3 = &N[6 * (n / 36)];
    const double *n1x = &Nx[6 * (n % 6)];
    const double *n2x = &Nx[6 * ((n % 36) / 6)];
    const double *n3x = &Nx[6 * (n / 36)];

    const TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar t1;
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[3] += n23 * (n1x[0] * v[0] + n1x[1] * v[3] + n1x[2] * v[6] +
                         n1x[3] * v[9] + n1x[4] * v[12] + n1x[5] * v[15]);
        t1 = (n1[0] * v[0] + n1[1] * v[3] + n1[2] * v[6] + n1[3] * v[9] +
              n1[4] * v[12] + n1[5] * v[15]);
        out[0] += n23 * t1;
        out[4] += n2x3 * t1;
        out[5] += n23x * t1;

        out[6] += n23 * (n1x[0] * v[1] + n1x[1] * v[4] + n1x[2] * v[7] +
                         n1x[3] * v[10] + n1x[4] * v[13] + n1x[5] * v[16]);
        t1 = (n1[0] * v[1] + n1[1] * v[4] + n1[2] * v[7] + n1[3] * v[10] +
              n1[4] * v[13] + n1[5] * v[16]);
        out[1] += n23 * t1;
        out[7] += n2x3 * t1;
        out[8] += n23x * t1;

        out[9] += n23 * (n1x[0] * v[2] + n1x[1] * v[5] + n1x[2] * v[8] +
                         n1x[3] * v[11] + n1x[4] * v[14] + n1x[5] * v[17]);
        t1 = (n1[0] * v[2] + n1[1] * v[5] + n1[2] * v[8] + n1[3] * v[11] +
              n1[4] * v[14] + n1[5] * v[17]);
        out[2] += n23 * t1;
        out[10] += n2x3 * t1;
        out[11] += n23x * t1;

        v += 18;
      }
    }
    out += 12;
  }
}

void TACSInterpAllTensor3DInterp6VarsPerNode4(const double N[],
                                              const double Nx[],
                                              const TacsScalar values[],
                                              TacsScalar out[]) {
  memset(out, 0, 4 * 4 * 216 * sizeof(TacsScalar));

  for (int n = 0; n < 216; n++) {
    const double *n1 = &N[6 * (n % 6)];
    const double *n2 = &N[6 * ((n % 36) / 6)];
    const double *n3 = &N[6 * (n / 36)];
    const double *n1x = &Nx[6 * (n % 6)];
    const double *n2x = &Nx[6 * ((n % 36) / 6)];
    const double *n3x = &Nx[6 * (n / 36)];

    const TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        TacsScalar t1;
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        out[4] += n23 * (n1x[0] * v[0] + n1x[1] * v[4] + n1x[2] * v[8] +
                         n1x[3] * v[12] + n1x[4] * v[16] + n1x[5] * v[20]);
        t1 = (n1[0] * v[0] + n1[1] * v[4] + n1[2] * v[8] + n1[3] * v[12] +
              n1[4] * v[16] + n1[5] * v[20]);
        out[0] += n23 * t1;
        out[5] += n2x3 * t1;
        out[6] += n23x * t1;

        out[7] += n23 * (n1x[0] * v[1] + n1x[1] * v[5] + n1x[2] * v[9] +
                         n1x[3] * v[13] + n1x[4] * v[17] + n1x[5] * v[21]);
        t1 = (n1[0] * v[1] + n1[1] * v[5] + n1[2] * v[9] + n1[3] * v[13] +
              n1[4] * v[17] + n1[5] * v[21]);
        out[1] += n23 * t1;
        out[8] += n2x3 * t1;
        out[9] += n23x * t1;

        out[10] += n23 * (n1x[0] * v[2] + n1x[1] * v[6] + n1x[2] * v[10] +
                          n1x[3] * v[14] + n1x[4] * v[18] + n1x[5] * v[22]);
        t1 = (n1[0] * v[2] + n1[1] * v[6] + n1[2] * v[10] + n1[3] * v[14] +
              n1[4] * v[18] + n1[5] * v[22]);
        out[2] += n23 * t1;
        out[11] += n2x3 * t1;
        out[12] += n23x * t1;

        out[13] += n23 * (n1x[0] * v[3] + n1x[1] * v[7] + n1x[2] * v[11] +
                          n1x[3] * v[15] + n1x[4] * v[19] + n1x[5] * v[23]);
        t1 = (n1[0] * v[3] + n1[1] * v[7] + n1[2] * v[11] + n1[3] * v[15] +
              n1[4] * v[19] + n1[5] * v[23]);
        out[3] += n23 * t1;
        out[14] += n2x3 * t1;
        out[15] += n23x * t1;

        v += 24;
      }
    }
    out += 16;
  }
}

void TacsAddAllTransTensor3DInterp6(const int m, const double N[],
                                    const double Nx[], const TacsScalar in[],
                                    TacsScalar values[]) {
  for (int n = 0; n < 216; n++) {
    const double *n1 = &N[6 * (n % 6)];
    const double *n2 = &N[6 * ((n % 36) / 6)];
    const double *n3 = &N[6 * (n / 36)];
    const double *n1x = &Nx[6 * (n % 6)];
    const double *n2x = &Nx[6 * ((n % 36) / 6)];
    const double *n3x = &Nx[6 * (n / 36)];

    TacsScalar *v = values;
    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        const TacsScalar *g = &in[m];
        for (int p = 0; p < m; p++) {
          TacsScalar a = n23 * in[p] + n2x3 * g[1] + n23x * g[2];
          TacsScalar b = n23 * g[0];
          v[0] += a * n1[0] + b * n1x[0];
          v[m] += a * n1[1] + b * n1x[1];
          v[2 * m] += a * n1[2] + b * n1x[2];
          v[3 * m] += a * n1[3] + b * n1x[3];
          v[4 * m] += a * n1[4] + b * n1x[4];
          v[5 * m] += a * n1[5] + b * n1x[5];
          g += 3;
          v++;
        }
        v += 5 * m;
      }
    }

    in += 4 * m;
  }
}

void TacsAddAllTransTensor3DInterp6VarsPerNode1(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 36; n++) {
    const double *n11 = &N[0];
    const double *n12 = &N[6];
    const double *n13 = &N[12];
    const double *n14 = &N[18];
    const double *n15 = &N[24];
    const double *n16 = &N[30];

    const double *n1x1 = &Nx[0];
    const double *n1x2 = &Nx[6];
    const double *n1x3 = &Nx[12];
    const double *n1x4 = &Nx[18];
    const double *n1x5 = &Nx[24];
    const double *n1x6 = &Nx[30];

    const double *n2 = &N[6 * (n % 6)];
    const double *n3 = &N[6 * (n / 6)];
    const double *n2x = &Nx[6 * (n % 6)];
    const double *n3x = &Nx[6 * (n / 6)];

    TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a[6], b[6];
        a[0] = n23 * in[0] + n2x3 * in[2] + n23x * in[3];
        b[0] = n23 * in[1];

        a[1] = n23 * in[4] + n2x3 * in[6] + n23x * in[7];
        b[1] = n23 * in[5];

        a[2] = n23 * in[8] + n2x3 * in[10] + n23x * in[11];
        b[2] = n23 * in[9];

        a[3] = n23 * in[12] + n2x3 * in[14] + n23x * in[15];
        b[3] = n23 * in[13];

        a[4] = n23 * in[16] + n2x3 * in[18] + n23x * in[19];
        b[4] = n23 * in[17];

        a[5] = n23 * in[20] + n2x3 * in[22] + n23x * in[23];
        b[5] = n23 * in[21];

        v[0] +=
            (a[0] * n11[0] + b[0] * n1x1[0] + a[1] * n12[0] + b[1] * n1x2[0] +
             a[2] * n13[0] + b[2] * n1x3[0] + a[3] * n14[0] + b[3] * n1x4[0] +
             a[4] * n15[0] + b[4] * n1x5[0] + a[5] * n16[0] + b[5] * n1x6[0]);
        v[1] +=
            (a[0] * n11[1] + b[0] * n1x1[1] + a[1] * n12[1] + b[1] * n1x2[1] +
             a[2] * n13[1] + b[2] * n1x3[1] + a[3] * n14[1] + b[3] * n1x4[1] +
             a[4] * n15[1] + b[4] * n1x5[1] + a[5] * n16[1] + b[5] * n1x6[1]);
        v[2] +=
            (a[0] * n11[2] + b[0] * n1x1[2] + a[1] * n12[2] + b[1] * n1x2[2] +
             a[2] * n13[2] + b[2] * n1x3[2] + a[3] * n14[2] + b[3] * n1x4[2] +
             a[4] * n15[2] + b[4] * n1x5[2] + a[5] * n16[2] + b[5] * n1x6[2]);
        v[3] +=
            (a[0] * n11[3] + b[0] * n1x1[3] + a[1] * n12[3] + b[1] * n1x2[3] +
             a[2] * n13[3] + b[2] * n1x3[3] + a[3] * n14[3] + b[3] * n1x4[3] +
             a[4] * n15[3] + b[4] * n1x5[3] + a[5] * n16[3] + b[5] * n1x6[3]);
        v[4] +=
            (a[0] * n11[4] + b[0] * n1x1[4] + a[1] * n12[4] + b[1] * n1x2[4] +
             a[2] * n13[4] + b[2] * n1x3[4] + a[3] * n14[4] + b[3] * n1x4[4] +
             a[4] * n15[4] + b[4] * n1x5[4] + a[5] * n16[4] + b[5] * n1x6[4]);
        v[5] +=
            (a[0] * n11[5] + b[0] * n1x1[5] + a[1] * n12[5] + b[1] * n1x2[5] +
             a[2] * n13[5] + b[2] * n1x3[5] + a[3] * n14[5] + b[3] * n1x4[5] +
             a[4] * n15[5] + b[4] * n1x5[5] + a[5] * n16[5] + b[5] * n1x6[5]);
        v += 6;
      }
    }
    in += 24;
  }
}

void TacsAddAllTransTensor3DInterp6VarsPerNode3(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 36; n++) {
    const double *n11 = &N[0];
    const double *n12 = &N[6];
    const double *n13 = &N[12];
    const double *n14 = &N[18];
    const double *n15 = &N[24];
    const double *n16 = &N[30];

    const double *n1x1 = &Nx[0];
    const double *n1x2 = &Nx[6];
    const double *n1x3 = &Nx[12];
    const double *n1x4 = &Nx[18];
    const double *n1x5 = &Nx[24];
    const double *n1x6 = &Nx[30];

    const double *n2 = &N[6 * (n % 6)];
    const double *n3 = &N[6 * (n / 6)];
    const double *n2x = &Nx[6 * (n % 6)];
    const double *n3x = &Nx[6 * (n / 6)];

    TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a[18], b[18];
        a[0] = n23 * in[0] + n2x3 * in[4] + n23x * in[5];
        b[0] = n23 * in[3];
        a[1] = n23 * in[1] + n2x3 * in[7] + n23x * in[8];
        b[1] = n23 * in[6];
        a[2] = n23 * in[2] + n2x3 * in[10] + n23x * in[11];
        b[2] = n23 * in[9];
        a[3] = n23 * in[12] + n2x3 * in[16] + n23x * in[17];
        b[3] = n23 * in[15];
        a[4] = n23 * in[13] + n2x3 * in[19] + n23x * in[20];
        b[4] = n23 * in[18];
        a[5] = n23 * in[14] + n2x3 * in[22] + n23x * in[23];
        b[5] = n23 * in[21];
        a[6] = n23 * in[24] + n2x3 * in[28] + n23x * in[29];
        b[6] = n23 * in[27];
        a[7] = n23 * in[25] + n2x3 * in[31] + n23x * in[32];
        b[7] = n23 * in[30];
        a[8] = n23 * in[26] + n2x3 * in[34] + n23x * in[35];
        b[8] = n23 * in[33];
        a[9] = n23 * in[36] + n2x3 * in[40] + n23x * in[41];
        b[9] = n23 * in[39];
        a[10] = n23 * in[37] + n2x3 * in[43] + n23x * in[44];
        b[10] = n23 * in[42];
        a[11] = n23 * in[38] + n2x3 * in[46] + n23x * in[47];
        b[11] = n23 * in[45];
        a[12] = n23 * in[48] + n2x3 * in[52] + n23x * in[53];
        b[12] = n23 * in[51];
        a[13] = n23 * in[49] + n2x3 * in[55] + n23x * in[56];
        b[13] = n23 * in[54];
        a[14] = n23 * in[50] + n2x3 * in[58] + n23x * in[59];
        b[14] = n23 * in[57];
        a[15] = n23 * in[60] + n2x3 * in[64] + n23x * in[65];
        b[15] = n23 * in[63];
        a[16] = n23 * in[61] + n2x3 * in[67] + n23x * in[68];
        b[16] = n23 * in[66];
        a[17] = n23 * in[62] + n2x3 * in[70] + n23x * in[71];
        b[17] = n23 * in[69];

        v[0] += (a[0] * n11[0] + b[0] * n1x1[0] + a[3] * n12[0] +
                 b[3] * n1x2[0] + a[6] * n13[0] + b[6] * n1x3[0] +
                 a[9] * n14[0] + b[9] * n1x4[0] + a[12] * n15[0] +
                 b[12] * n1x5[0] + a[15] * n16[0] + b[15] * n1x6[0]);
        v[1] += (a[1] * n11[0] + b[1] * n1x1[0] + a[4] * n12[0] +
                 b[4] * n1x2[0] + a[7] * n13[0] + b[7] * n1x3[0] +
                 a[10] * n14[0] + b[10] * n1x4[0] + a[13] * n15[0] +
                 b[13] * n1x5[0] + a[16] * n16[0] + b[16] * n1x6[0]);
        v[2] += (a[2] * n11[0] + b[2] * n1x1[0] + a[5] * n12[0] +
                 b[5] * n1x2[0] + a[8] * n13[0] + b[8] * n1x3[0] +
                 a[11] * n14[0] + b[11] * n1x4[0] + a[14] * n15[0] +
                 b[14] * n1x5[0] + a[17] * n16[0] + b[17] * n1x6[0]);
        v[3] += (a[0] * n11[1] + b[0] * n1x1[1] + a[3] * n12[1] +
                 b[3] * n1x2[1] + a[6] * n13[1] + b[6] * n1x3[1] +
                 a[9] * n14[1] + b[9] * n1x4[1] + a[12] * n15[1] +
                 b[12] * n1x5[1] + a[15] * n16[1] + b[15] * n1x6[1]);
        v[4] += (a[1] * n11[1] + b[1] * n1x1[1] + a[4] * n12[1] +
                 b[4] * n1x2[1] + a[7] * n13[1] + b[7] * n1x3[1] +
                 a[10] * n14[1] + b[10] * n1x4[1] + a[13] * n15[1] +
                 b[13] * n1x5[1] + a[16] * n16[1] + b[16] * n1x6[1]);
        v[5] += (a[2] * n11[1] + b[2] * n1x1[1] + a[5] * n12[1] +
                 b[5] * n1x2[1] + a[8] * n13[1] + b[8] * n1x3[1] +
                 a[11] * n14[1] + b[11] * n1x4[1] + a[14] * n15[1] +
                 b[14] * n1x5[1] + a[17] * n16[1] + b[17] * n1x6[1]);
        v[6] += (a[0] * n11[2] + b[0] * n1x1[2] + a[3] * n12[2] +
                 b[3] * n1x2[2] + a[6] * n13[2] + b[6] * n1x3[2] +
                 a[9] * n14[2] + b[9] * n1x4[2] + a[12] * n15[2] +
                 b[12] * n1x5[2] + a[15] * n16[2] + b[15] * n1x6[2]);
        v[7] += (a[1] * n11[2] + b[1] * n1x1[2] + a[4] * n12[2] +
                 b[4] * n1x2[2] + a[7] * n13[2] + b[7] * n1x3[2] +
                 a[10] * n14[2] + b[10] * n1x4[2] + a[13] * n15[2] +
                 b[13] * n1x5[2] + a[16] * n16[2] + b[16] * n1x6[2]);
        v[8] += (a[2] * n11[2] + b[2] * n1x1[2] + a[5] * n12[2] +
                 b[5] * n1x2[2] + a[8] * n13[2] + b[8] * n1x3[2] +
                 a[11] * n14[2] + b[11] * n1x4[2] + a[14] * n15[2] +
                 b[14] * n1x5[2] + a[17] * n16[2] + b[17] * n1x6[2]);
        v[9] += (a[0] * n11[3] + b[0] * n1x1[3] + a[3] * n12[3] +
                 b[3] * n1x2[3] + a[6] * n13[3] + b[6] * n1x3[3] +
                 a[9] * n14[3] + b[9] * n1x4[3] + a[12] * n15[3] +
                 b[12] * n1x5[3] + a[15] * n16[3] + b[15] * n1x6[3]);
        v[10] += (a[1] * n11[3] + b[1] * n1x1[3] + a[4] * n12[3] +
                  b[4] * n1x2[3] + a[7] * n13[3] + b[7] * n1x3[3] +
                  a[10] * n14[3] + b[10] * n1x4[3] + a[13] * n15[3] +
                  b[13] * n1x5[3] + a[16] * n16[3] + b[16] * n1x6[3]);
        v[11] += (a[2] * n11[3] + b[2] * n1x1[3] + a[5] * n12[3] +
                  b[5] * n1x2[3] + a[8] * n13[3] + b[8] * n1x3[3] +
                  a[11] * n14[3] + b[11] * n1x4[3] + a[14] * n15[3] +
                  b[14] * n1x5[3] + a[17] * n16[3] + b[17] * n1x6[3]);
        v[12] += (a[0] * n11[4] + b[0] * n1x1[4] + a[3] * n12[4] +
                  b[3] * n1x2[4] + a[6] * n13[4] + b[6] * n1x3[4] +
                  a[9] * n14[4] + b[9] * n1x4[4] + a[12] * n15[4] +
                  b[12] * n1x5[4] + a[15] * n16[4] + b[15] * n1x6[4]);
        v[13] += (a[1] * n11[4] + b[1] * n1x1[4] + a[4] * n12[4] +
                  b[4] * n1x2[4] + a[7] * n13[4] + b[7] * n1x3[4] +
                  a[10] * n14[4] + b[10] * n1x4[4] + a[13] * n15[4] +
                  b[13] * n1x5[4] + a[16] * n16[4] + b[16] * n1x6[4]);
        v[14] += (a[2] * n11[4] + b[2] * n1x1[4] + a[5] * n12[4] +
                  b[5] * n1x2[4] + a[8] * n13[4] + b[8] * n1x3[4] +
                  a[11] * n14[4] + b[11] * n1x4[4] + a[14] * n15[4] +
                  b[14] * n1x5[4] + a[17] * n16[4] + b[17] * n1x6[4]);
        v[15] += (a[0] * n11[5] + b[0] * n1x1[5] + a[3] * n12[5] +
                  b[3] * n1x2[5] + a[6] * n13[5] + b[6] * n1x3[5] +
                  a[9] * n14[5] + b[9] * n1x4[5] + a[12] * n15[5] +
                  b[12] * n1x5[5] + a[15] * n16[5] + b[15] * n1x6[5]);
        v[16] += (a[1] * n11[5] + b[1] * n1x1[5] + a[4] * n12[5] +
                  b[4] * n1x2[5] + a[7] * n13[5] + b[7] * n1x3[5] +
                  a[10] * n14[5] + b[10] * n1x4[5] + a[13] * n15[5] +
                  b[13] * n1x5[5] + a[16] * n16[5] + b[16] * n1x6[5]);
        v[17] += (a[2] * n11[5] + b[2] * n1x1[5] + a[5] * n12[5] +
                  b[5] * n1x2[5] + a[8] * n13[5] + b[8] * n1x3[5] +
                  a[11] * n14[5] + b[11] * n1x4[5] + a[14] * n15[5] +
                  b[14] * n1x5[5] + a[17] * n16[5] + b[17] * n1x6[5]);

        v += 18;
      }
    }
    in += 72;
  }
}

void TacsAddAllTransTensor3DInterp6VarsPerNode4(const double N[],
                                                const double Nx[],
                                                const TacsScalar in[],
                                                TacsScalar values[]) {
  for (int n = 0; n < 36; n++) {
    const double *n11 = &N[0];
    const double *n12 = &N[6];
    const double *n13 = &N[12];
    const double *n14 = &N[18];
    const double *n15 = &N[24];
    const double *n16 = &N[30];

    const double *n1x1 = &Nx[0];
    const double *n1x2 = &Nx[6];
    const double *n1x3 = &Nx[12];
    const double *n1x4 = &Nx[18];
    const double *n1x5 = &Nx[24];
    const double *n1x6 = &Nx[30];

    const double *n2 = &N[6 * (n % 6)];
    const double *n3 = &N[6 * (n / 6)];
    const double *n2x = &Nx[6 * (n % 6)];
    const double *n3x = &Nx[6 * (n / 6)];

    TacsScalar *v = values;

    for (int k = 0; k < 6; k++) {
      for (int j = 0; j < 6; j++) {
        double n23 = n2[j] * n3[k];
        double n2x3 = n2x[j] * n3[k];
        double n23x = n2[j] * n3x[k];

        TacsScalar a[24], b[24];
        a[0] = n23 * in[0] + n2x3 * in[5] + n23x * in[6];
        b[0] = n23 * in[4];
        a[1] = n23 * in[1] + n2x3 * in[8] + n23x * in[9];
        b[1] = n23 * in[7];
        a[2] = n23 * in[2] + n2x3 * in[11] + n23x * in[12];
        b[2] = n23 * in[10];
        a[3] = n23 * in[3] + n2x3 * in[14] + n23x * in[15];
        b[3] = n23 * in[13];
        a[4] = n23 * in[16] + n2x3 * in[21] + n23x * in[22];
        b[4] = n23 * in[20];
        a[5] = n23 * in[17] + n2x3 * in[24] + n23x * in[25];
        b[5] = n23 * in[23];
        a[6] = n23 * in[18] + n2x3 * in[27] + n23x * in[28];
        b[6] = n23 * in[26];
        a[7] = n23 * in[19] + n2x3 * in[30] + n23x * in[31];
        b[7] = n23 * in[29];
        a[8] = n23 * in[32] + n2x3 * in[37] + n23x * in[38];
        b[8] = n23 * in[36];
        a[9] = n23 * in[33] + n2x3 * in[40] + n23x * in[41];
        b[9] = n23 * in[39];
        a[10] = n23 * in[34] + n2x3 * in[43] + n23x * in[44];
        b[10] = n23 * in[42];
        a[11] = n23 * in[35] + n2x3 * in[46] + n23x * in[47];
        b[11] = n23 * in[45];
        a[12] = n23 * in[48] + n2x3 * in[53] + n23x * in[54];
        b[12] = n23 * in[52];
        a[13] = n23 * in[49] + n2x3 * in[56] + n23x * in[57];
        b[13] = n23 * in[55];
        a[14] = n23 * in[50] + n2x3 * in[59] + n23x * in[60];
        b[14] = n23 * in[58];
        a[15] = n23 * in[51] + n2x3 * in[62] + n23x * in[63];
        b[15] = n23 * in[61];
        a[16] = n23 * in[64] + n2x3 * in[69] + n23x * in[70];
        b[16] = n23 * in[68];
        a[17] = n23 * in[65] + n2x3 * in[72] + n23x * in[73];
        b[17] = n23 * in[71];
        a[18] = n23 * in[66] + n2x3 * in[75] + n23x * in[76];
        b[18] = n23 * in[74];
        a[19] = n23 * in[67] + n2x3 * in[78] + n23x * in[79];
        b[19] = n23 * in[77];
        a[20] = n23 * in[80] + n2x3 * in[85] + n23x * in[86];
        b[20] = n23 * in[84];
        a[21] = n23 * in[81] + n2x3 * in[88] + n23x * in[89];
        b[21] = n23 * in[87];
        a[22] = n23 * in[82] + n2x3 * in[91] + n23x * in[92];
        b[22] = n23 * in[90];
        a[23] = n23 * in[83] + n2x3 * in[94] + n23x * in[95];
        b[23] = n23 * in[93];

        v[0] += (a[0] * n11[0] + b[0] * n1x1[0] + a[4] * n12[0] +
                 b[4] * n1x2[0] + a[8] * n13[0] + b[8] * n1x3[0] +
                 a[12] * n14[0] + b[12] * n1x4[0] + a[16] * n15[0] +
                 b[16] * n1x5[0] + a[20] * n16[0] + b[20] * n1x6[0]);
        v[1] += (a[1] * n11[0] + b[1] * n1x1[0] + a[5] * n12[0] +
                 b[5] * n1x2[0] + a[9] * n13[0] + b[9] * n1x3[0] +
                 a[13] * n14[0] + b[13] * n1x4[0] + a[17] * n15[0] +
                 b[17] * n1x5[0] + a[21] * n16[0] + b[21] * n1x6[0]);
        v[2] += (a[2] * n11[0] + b[2] * n1x1[0] + a[6] * n12[0] +
                 b[6] * n1x2[0] + a[10] * n13[0] + b[10] * n1x3[0] +
                 a[14] * n14[0] + b[14] * n1x4[0] + a[18] * n15[0] +
                 b[18] * n1x5[0] + a[22] * n16[0] + b[22] * n1x6[0]);
        v[3] += (a[3] * n11[0] + b[3] * n1x1[0] + a[7] * n12[0] +
                 b[7] * n1x2[0] + a[11] * n13[0] + b[11] * n1x3[0] +
                 a[15] * n14[0] + b[15] * n1x4[0] + a[19] * n15[0] +
                 b[19] * n1x5[0] + a[23] * n16[0] + b[23] * n1x6[0]);
        v[4] += (a[0] * n11[1] + b[0] * n1x1[1] + a[4] * n12[1] +
                 b[4] * n1x2[1] + a[8] * n13[1] + b[8] * n1x3[1] +
                 a[12] * n14[1] + b[12] * n1x4[1] + a[16] * n15[1] +
                 b[16] * n1x5[1] + a[20] * n16[1] + b[20] * n1x6[1]);
        v[5] += (a[1] * n11[1] + b[1] * n1x1[1] + a[5] * n12[1] +
                 b[5] * n1x2[1] + a[9] * n13[1] + b[9] * n1x3[1] +
                 a[13] * n14[1] + b[13] * n1x4[1] + a[17] * n15[1] +
                 b[17] * n1x5[1] + a[21] * n16[1] + b[21] * n1x6[1]);
        v[6] += (a[2] * n11[1] + b[2] * n1x1[1] + a[6] * n12[1] +
                 b[6] * n1x2[1] + a[10] * n13[1] + b[10] * n1x3[1] +
                 a[14] * n14[1] + b[14] * n1x4[1] + a[18] * n15[1] +
                 b[18] * n1x5[1] + a[22] * n16[1] + b[22] * n1x6[1]);
        v[7] += (a[3] * n11[1] + b[3] * n1x1[1] + a[7] * n12[1] +
                 b[7] * n1x2[1] + a[11] * n13[1] + b[11] * n1x3[1] +
                 a[15] * n14[1] + b[15] * n1x4[1] + a[19] * n15[1] +
                 b[19] * n1x5[1] + a[23] * n16[1] + b[23] * n1x6[1]);
        v[8] += (a[0] * n11[2] + b[0] * n1x1[2] + a[4] * n12[2] +
                 b[4] * n1x2[2] + a[8] * n13[2] + b[8] * n1x3[2] +
                 a[12] * n14[2] + b[12] * n1x4[2] + a[16] * n15[2] +
                 b[16] * n1x5[2] + a[20] * n16[2] + b[20] * n1x6[2]);
        v[9] += (a[1] * n11[2] + b[1] * n1x1[2] + a[5] * n12[2] +
                 b[5] * n1x2[2] + a[9] * n13[2] + b[9] * n1x3[2] +
                 a[13] * n14[2] + b[13] * n1x4[2] + a[17] * n15[2] +
                 b[17] * n1x5[2] + a[21] * n16[2] + b[21] * n1x6[2]);
        v[10] += (a[2] * n11[2] + b[2] * n1x1[2] + a[6] * n12[2] +
                  b[6] * n1x2[2] + a[10] * n13[2] + b[10] * n1x3[2] +
                  a[14] * n14[2] + b[14] * n1x4[2] + a[18] * n15[2] +
                  b[18] * n1x5[2] + a[22] * n16[2] + b[22] * n1x6[2]);
        v[11] += (a[3] * n11[2] + b[3] * n1x1[2] + a[7] * n12[2] +
                  b[7] * n1x2[2] + a[11] * n13[2] + b[11] * n1x3[2] +
                  a[15] * n14[2] + b[15] * n1x4[2] + a[19] * n15[2] +
                  b[19] * n1x5[2] + a[23] * n16[2] + b[23] * n1x6[2]);
        v[12] += (a[0] * n11[3] + b[0] * n1x1[3] + a[4] * n12[3] +
                  b[4] * n1x2[3] + a[8] * n13[3] + b[8] * n1x3[3] +
                  a[12] * n14[3] + b[12] * n1x4[3] + a[16] * n15[3] +
                  b[16] * n1x5[3] + a[20] * n16[3] + b[20] * n1x6[3]);
        v[13] += (a[1] * n11[3] + b[1] * n1x1[3] + a[5] * n12[3] +
                  b[5] * n1x2[3] + a[9] * n13[3] + b[9] * n1x3[3] +
                  a[13] * n14[3] + b[13] * n1x4[3] + a[17] * n15[3] +
                  b[17] * n1x5[3] + a[21] * n16[3] + b[21] * n1x6[3]);
        v[14] += (a[2] * n11[3] + b[2] * n1x1[3] + a[6] * n12[3] +
                  b[6] * n1x2[3] + a[10] * n13[3] + b[10] * n1x3[3] +
                  a[14] * n14[3] + b[14] * n1x4[3] + a[18] * n15[3] +
                  b[18] * n1x5[3] + a[22] * n16[3] + b[22] * n1x6[3]);
        v[15] += (a[3] * n11[3] + b[3] * n1x1[3] + a[7] * n12[3] +
                  b[7] * n1x2[3] + a[11] * n13[3] + b[11] * n1x3[3] +
                  a[15] * n14[3] + b[15] * n1x4[3] + a[19] * n15[3] +
                  b[19] * n1x5[3] + a[23] * n16[3] + b[23] * n1x6[3]);
        v[16] += (a[0] * n11[4] + b[0] * n1x1[4] + a[4] * n12[4] +
                  b[4] * n1x2[4] + a[8] * n13[4] + b[8] * n1x3[4] +
                  a[12] * n14[4] + b[12] * n1x4[4] + a[16] * n15[4] +
                  b[16] * n1x5[4] + a[20] * n16[4] + b[20] * n1x6[4]);
        v[17] += (a[1] * n11[4] + b[1] * n1x1[4] + a[5] * n12[4] +
                  b[5] * n1x2[4] + a[9] * n13[4] + b[9] * n1x3[4] +
                  a[13] * n14[4] + b[13] * n1x4[4] + a[17] * n15[4] +
                  b[17] * n1x5[4] + a[21] * n16[4] + b[21] * n1x6[4]);
        v[18] += (a[2] * n11[4] + b[2] * n1x1[4] + a[6] * n12[4] +
                  b[6] * n1x2[4] + a[10] * n13[4] + b[10] * n1x3[4] +
                  a[14] * n14[4] + b[14] * n1x4[4] + a[18] * n15[4] +
                  b[18] * n1x5[4] + a[22] * n16[4] + b[22] * n1x6[4]);
        v[19] += (a[3] * n11[4] + b[3] * n1x1[4] + a[7] * n12[4] +
                  b[7] * n1x2[4] + a[11] * n13[4] + b[11] * n1x3[4] +
                  a[15] * n14[4] + b[15] * n1x4[4] + a[19] * n15[4] +
                  b[19] * n1x5[4] + a[23] * n16[4] + b[23] * n1x6[4]);
        v[20] += (a[0] * n11[5] + b[0] * n1x1[5] + a[4] * n12[5] +
                  b[4] * n1x2[5] + a[8] * n13[5] + b[8] * n1x3[5] +
                  a[12] * n14[5] + b[12] * n1x4[5] + a[16] * n15[5] +
                  b[16] * n1x5[5] + a[20] * n16[5] + b[20] * n1x6[5]);
        v[21] += (a[1] * n11[5] + b[1] * n1x1[5] + a[5] * n12[5] +
                  b[5] * n1x2[5] + a[9] * n13[5] + b[9] * n1x3[5] +
                  a[13] * n14[5] + b[13] * n1x4[5] + a[17] * n15[5] +
                  b[17] * n1x5[5] + a[21] * n16[5] + b[21] * n1x6[5]);
        v[22] += (a[2] * n11[5] + b[2] * n1x1[5] + a[6] * n12[5] +
                  b[6] * n1x2[5] + a[10] * n13[5] + b[10] * n1x3[5] +
                  a[14] * n14[5] + b[14] * n1x4[5] + a[18] * n15[5] +
                  b[18] * n1x5[5] + a[22] * n16[5] + b[22] * n1x6[5]);
        v[23] += (a[3] * n11[5] + b[3] * n1x1[5] + a[7] * n12[5] +
                  b[7] * n1x2[5] + a[11] * n13[5] + b[11] * n1x3[5] +
                  a[15] * n14[5] + b[15] * n1x4[5] + a[19] * n15[5] +
                  b[19] * n1x5[5] + a[23] * n16[5] + b[23] * n1x6[5]);
        v += 24;
      }
    }
    in += 96;
  }
}
