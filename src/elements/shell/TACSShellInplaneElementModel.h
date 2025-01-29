#ifndef TACS_SHELL_INPLANE_ELEMENT_MODEL_H
#define TACS_SHELL_INPLANE_ELEMENT_MODEL_H

#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"
#include "TACSShellElementQuadBasis.h"

class TACSShellInplaneLinearModel {
 public:
  /**
    Compute the tensorial components of the tying strain

    G = 0.5*(X,eta^{T}*U,eta + U,eta^{T}*X,eta)

    The derivative with respect to the frame gives

    X,eta = [X,xi ; n]

    The derivative with respect to the displacements gives

    u,eta = [u,xi ; d]

    @param Xxi Derivatives of the node locations with respect to xi
    @param n The interpolated frame normal
    @param Uxi Derivatives of the displacements with respect to xi
    @param d The interpolated director field
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrain(const TacsScalar Xpts[], const TacsScalar fn[],
                                 const TacsScalar vars[], const TacsScalar d[],
                                 TacsScalar ety[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6], d0[3], n0[3];

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, fn, n0);

      ety[index] = 0.0;
      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        ety[index] = 0.5 * (Xxi[1] * d0[0] + Xxi[3] * d0[1] + Xxi[5] * d0[2] +
                            n0[0] * Uxi[1] + n0[1] * Uxi[3] + n0[2] * Uxi[5]);
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[2] * d0[1] + Xxi[4] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[2] + n0[2] * Uxi[4]);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose(
      const TacsScalar Xpts[], const TacsScalar fn[], const TacsScalar vars[],
      const TacsScalar d[], const TacsScalar dety[], TacsScalar res[],
      TacsScalar dd[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6], dUxi[6];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

      TacsScalar d0[3], dd0[3], n0[3];
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, fn, n0);

      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi[0] = 0.0;
        dUxi[1] = 0.5 * dety[index] * n0[0];
        dUxi[2] = 0.0;
        dUxi[3] = 0.5 * dety[index] * n0[1];
        dUxi[4] = 0.0;
        dUxi[5] = 0.5 * dety[index] * n0[2];

        dd0[0] = 0.5 * dety[index] * Xxi[1];
        dd0[1] = 0.5 * dety[index] * Xxi[3];
        dd0[2] = 0.5 * dety[index] * Xxi[5];
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi[0] = 0.5 * dety[index] * n0[0];
        dUxi[1] = 0.0;
        dUxi[2] = 0.5 * dety[index] * n0[1];
        dUxi[3] = 0.0;
        dUxi[4] = 0.5 * dety[index] * n0[2];
        dUxi[5] = 0.0;

        dd0[0] = 0.5 * dety[index] * Xxi[0];
        dd0[1] = 0.5 * dety[index] * Xxi[2];
        dd0[2] = 0.5 * dety[index] * Xxi[4];
      }

      basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
      if (res) {
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi,
                                                                       res);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian(
      const TacsScalar alpha, const TacsScalar Xpts[], const TacsScalar fn[],
      const TacsScalar vars[], const TacsScalar d[], const TacsScalar dety[],
      const TacsScalar d2ety[], const TacsScalar d2etyu[],
      const TacsScalar d2etyd[], TacsScalar mat[], TacsScalar d2d[],
      TacsScalar d2du[]) {
    // Initialize the data
    TacsScalar n0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity;

    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      n0 += 3;
      Xxi += 6;
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity;
    for (int i1 = 0; i1 < basis::NUM_TYING_POINTS; i1++, n01 += 3, Xxi1 += 6) {
      // Get the field index
      const TacsShellTyingStrainComponent f1 = basis::getTyingField(i1);

      // Get the tying point parametric location
      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar du2[3 * basis::NUM_NODES], dd2[3 * basis::NUM_NODES];
      memset(du2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      TacsScalar *n02 = n0ty, *Xxi2 = Xxity;
      for (int i2 = 0; i2 < basis::NUM_TYING_POINTS;
           i2++, n02 += 3, Xxi2 += 6) {
        // Get the field index
        const TacsShellTyingStrainComponent f2 = basis::getTyingField(i2);

        // Get the tying point parametric location
        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        TacsScalar value = d2ety[basis::NUM_TYING_POINTS * i1 + i2];

        TacsScalar dUxi2[6], dd02[3];
        if (f2 == TACS_SHELL_G23_COMPONENT) {
          // Compute g23 = e2^{T}*G*e3
          dUxi2[0] = 0.0;
          dUxi2[1] = 0.5 * value * n02[0];
          dUxi2[2] = 0.0;
          dUxi2[3] = 0.5 * value * n02[1];
          dUxi2[4] = 0.0;
          dUxi2[5] = 0.5 * value * n02[2];

          dd02[0] = 0.5 * value * Xxi2[1];
          dd02[1] = 0.5 * value * Xxi2[3];
          dd02[2] = 0.5 * value * Xxi2[5];
        } else if (f2 == TACS_SHELL_G13_COMPONENT) {
          // Compute g13 = e1^{T}*G*e3
          dUxi2[0] = 0.5 * value * n02[0];
          dUxi2[1] = 0.0;
          dUxi2[2] = 0.5 * value * n02[1];
          dUxi2[3] = 0.0;
          dUxi2[4] = 0.5 * value * n02[2];
          dUxi2[5] = 0.0;

          dd02[0] = 0.5 * value * Xxi2[0];
          dd02[1] = 0.5 * value * Xxi2[2];
          dd02[2] = 0.5 * value * Xxi2[4];
        }

        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
      }

      TacsScalar du1[3 * basis::NUM_NODES], dd1[3 * basis::NUM_NODES];
      memset(du1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      // Store the the derivative information for the first point
      TacsScalar dUxi1[6], dd01[3];
      if (f1 == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi1[0] = 0.0;
        dUxi1[1] = 0.5 * n01[0];
        dUxi1[2] = 0.0;
        dUxi1[3] = 0.5 * n01[1];
        dUxi1[4] = 0.0;
        dUxi1[5] = 0.5 * n01[2];

        dd01[0] = 0.5 * Xxi1[1];
        dd01[1] = 0.5 * Xxi1[3];
        dd01[2] = 0.5 * Xxi1[5];
      } else if (f1 == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi1[0] = 0.5 * n01[0];
        dUxi1[1] = 0.0;
        dUxi1[2] = 0.5 * n01[1];
        dUxi1[3] = 0.0;
        dUxi1[4] = 0.5 * n01[2];
        dUxi1[5] = 0.0;

        dd01[0] = 0.5 * Xxi1[0];
        dd01[1] = 0.5 * Xxi1[2];
        dd01[2] = 0.5 * Xxi1[4];
      }

      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

      const TacsScalar *etd = &d2etyd[3 * basis::NUM_NODES * i1];
      const TacsScalar *etu = &d2etyu[3 * basis::NUM_NODES * i1];
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2d[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * dd2[j] + dd1[i] * etd[j] + etd[i] * dd1[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * du2[j] + dd1[i] * etu[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] += etd[i] * du1[j];
        }
      }

      const int nvars = vars_per_node * basis::NUM_NODES;
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        int ii = vars_per_node * (i / 3) + (i % 3);
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          int jj = vars_per_node * (j / 3) + (j % 3);
          mat[nvars * ii + jj] +=
              du1[i] * du2[j] + du1[i] * etu[j] + etu[i] * du1[j];
        }
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessianDeriv(
      const TacsScalar alpha, const TacsScalar Xpts[], const TacsScalar fn[],
      const TacsScalar vars[], const TacsScalar d[], const TacsScalar dety[],
      const TacsScalar d2ety[], const TacsScalar d2etyu[],
      const TacsScalar d2etyd[], const TacsScalar vars_d[],
      const TacsScalar d_d[], const TacsScalar dety_d[],
      const TacsScalar d2ety_d[], const TacsScalar d2etyu_d[],
      const TacsScalar d2etyd_d[], TacsScalar mat[], TacsScalar d2d[],
      TacsScalar d2du[], TacsScalar mat_d[], TacsScalar d2d_d[],
      TacsScalar d2du_d[]) {
    // Initialize the data
    TacsScalar n0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity;

    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      n0 += 3;
      Xxi += 6;
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity;
    for (int i1 = 0; i1 < basis::NUM_TYING_POINTS; i1++, n01 += 3, Xxi1 += 6) {
      // Get the field index
      const TacsShellTyingStrainComponent f1 = basis::getTyingField(i1);

      // Get the tying point parametric location
      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar du2[3 * basis::NUM_NODES], dd2[3 * basis::NUM_NODES];
      memset(du2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      TacsScalar du2d[3 * basis::NUM_NODES], dd2d[3 * basis::NUM_NODES];
      memset(du2d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      TacsScalar *n02 = n0ty, *Xxi2 = Xxity;
      for (int i2 = 0; i2 < basis::NUM_TYING_POINTS;
           i2++, n02 += 3, Xxi2 += 6) {
        // Get the field index
        const TacsShellTyingStrainComponent f2 = basis::getTyingField(i2);

        // Get the tying point parametric location
        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        TacsScalar value = d2ety[basis::NUM_TYING_POINTS * i1 + i2];
        TacsScalar valued = d2ety_d[basis::NUM_TYING_POINTS * i1 + i2];

        TacsScalar dUxi2[6], dd02[3];
        TacsScalar dUxi2d[6], dd02d[3];
        if (f2 == TACS_SHELL_G23_COMPONENT) {
          // Compute g23 = e2^{T}*G*e3
          dUxi2[0] = 0.0;
          dUxi2[1] = 0.5 * value * n02[0];
          dUxi2[2] = 0.0;
          dUxi2[3] = 0.5 * value * n02[1];
          dUxi2[4] = 0.0;
          dUxi2[5] = 0.5 * value * n02[2];

          dUxi2d[0] = 0.0;
          dUxi2d[1] = 0.5 * valued * n02[0];
          dUxi2d[2] = 0.0;
          dUxi2d[3] = 0.5 * valued * n02[1];
          dUxi2d[4] = 0.0;
          dUxi2d[5] = 0.5 * valued * n02[2];

          dd02[0] = 0.5 * value * Xxi2[1];
          dd02[1] = 0.5 * value * Xxi2[3];
          dd02[2] = 0.5 * value * Xxi2[5];

          dd02d[0] = 0.5 * valued * Xxi2[1];
          dd02d[1] = 0.5 * valued * Xxi2[3];
          dd02d[2] = 0.5 * valued * Xxi2[5];
        } else if (f2 == TACS_SHELL_G13_COMPONENT) {
          // Compute g13 = e1^{T}*G*e3
          dUxi2[0] = 0.5 * value * n02[0];
          dUxi2[1] = 0.0;
          dUxi2[2] = 0.5 * value * n02[1];
          dUxi2[3] = 0.0;
          dUxi2[4] = 0.5 * value * n02[2];
          dUxi2[5] = 0.0;

          dUxi2d[0] = 0.5 * valued * n02[0];
          dUxi2d[1] = 0.0;
          dUxi2d[2] = 0.5 * valued * n02[1];
          dUxi2d[3] = 0.0;
          dUxi2d[4] = 0.5 * valued * n02[2];
          dUxi2d[5] = 0.0;

          dd02[0] = 0.5 * value * Xxi2[0];
          dd02[1] = 0.5 * value * Xxi2[2];
          dd02[2] = 0.5 * value * Xxi2[4];

          dd02d[0] = 0.5 * valued * Xxi2[0];
          dd02d[1] = 0.5 * valued * Xxi2[2];
          dd02d[2] = 0.5 * valued * Xxi2[4];
        }

        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02d, dd2d);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2d, du2d);
      }

      TacsScalar du1[3 * basis::NUM_NODES], dd1[3 * basis::NUM_NODES];
      memset(du1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      // Store the the derivative information for the first point
      TacsScalar dUxi1[6], dd01[3];
      if (f1 == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi1[0] = 0.0;
        dUxi1[1] = 0.5 * n01[0];
        dUxi1[2] = 0.0;
        dUxi1[3] = 0.5 * n01[1];
        dUxi1[4] = 0.0;
        dUxi1[5] = 0.5 * n01[2];

        dd01[0] = 0.5 * Xxi1[1];
        dd01[1] = 0.5 * Xxi1[3];
        dd01[2] = 0.5 * Xxi1[5];
      } else if (f1 == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi1[0] = 0.5 * n01[0];
        dUxi1[1] = 0.0;
        dUxi1[2] = 0.5 * n01[1];
        dUxi1[3] = 0.0;
        dUxi1[4] = 0.5 * n01[2];
        dUxi1[5] = 0.0;

        dd01[0] = 0.5 * Xxi1[0];
        dd01[1] = 0.5 * Xxi1[2];
        dd01[2] = 0.5 * Xxi1[4];
      }

      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

      const TacsScalar *etd = &d2etyd[3 * basis::NUM_NODES * i1];
      const TacsScalar *etu = &d2etyu[3 * basis::NUM_NODES * i1];
      const TacsScalar *etdd = &d2etyd_d[3 * basis::NUM_NODES * i1];
      const TacsScalar *etud = &d2etyu_d[3 * basis::NUM_NODES * i1];
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2d[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * dd2[j] + dd1[i] * etd[j] + etd[i] * dd1[j];
          d2d_d[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * etdd[j] + etdd[i] * dd1[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * du2[j] + dd1[i] * etu[j];
          d2du_d[3 * basis::NUM_NODES * i + j] += dd1[i] * etud[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] += etd[i] * du1[j];
          d2du_d[3 * basis::NUM_NODES * i + j] += etdd[i] * du1[j];
        }
      }

      const int nvars = vars_per_node * basis::NUM_NODES;
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        int ii = vars_per_node * (i / 3) + (i % 3);
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          int jj = vars_per_node * (j / 3) + (j % 3);
          if (mat) {
            mat[nvars * ii + jj] +=
                du1[i] * du2[j] + du1[i] * etu[j] + etu[i] * du1[j];
          }
          mat_d[nvars * ii + jj] += du1[i] * etud[j] + etud[i] * du1[j];
        }
      }
    }
  }

  /*
    Compute the directional derivative
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv(
      const TacsScalar Xpts[], const TacsScalar fn[], const TacsScalar vars[],
      const TacsScalar d[], const TacsScalar varsd[], const TacsScalar dd[],
      TacsScalar ety[], TacsScalar etyd[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6], Uxid[6];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, Uxid);

      TacsScalar d0[3], d0d[3], n0[3];
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, dd, d0d);
      basis::template interpFields<3, 3>(pt, fn, n0);

      ety[index] = 0.0;
      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        ety[index] = 0.5 * (Xxi[1] * d0[0] + Xxi[3] * d0[1] + Xxi[5] * d0[2] +
                            n0[0] * Uxi[1] + n0[1] * Uxi[3] + n0[2] * Uxi[5]);
        etyd[index] =
            0.5 * (Xxi[1] * d0d[0] + Xxi[3] * d0d[1] + Xxi[5] * d0d[2] +
                   n0[0] * Uxid[1] + n0[1] * Uxid[3] + n0[2] * Uxid[5]);
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[2] * d0[1] + Xxi[4] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[2] + n0[2] * Uxi[4]);
        etyd[index] =
            0.5 * (Xxi[0] * d0d[0] + Xxi[2] * d0d[1] + Xxi[4] * d0d[2] +
                   n0[0] * Uxid[0] + n0[1] * Uxid[2] + n0[2] * Uxid[4]);
      }
    }
  }

  /*
    Evaluate the strain as a function of the displacement derivatives
    and interpolated strain from the tensorial components
  */
  static inline void evalStrain(const TacsScalar u0x[], const TacsScalar u1x[],
                                const TacsScalar e0ty[], TacsScalar e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0];
    e[1] = u0x[4];
    e[2] = u0x[1] + u0x[3];

    // Compute the bending strain
    e[3] = u1x[0];
    e[4] = u1x[4];
    e[5] = u1x[1] + u1x[3];

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4];
    e[7] = 2.0 * e0ty[2];
  }

  /**
    Evaluate the derivative of the strain
  */
  static inline void evalStrainSens(const TacsScalar scale,
                                    const TacsScalar dfde[],
                                    const TacsScalar u0x[],
                                    const TacsScalar u1x[], TacsScalar du0x[],
                                    TacsScalar du1x[], TacsScalar de0ty[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = 0.0;
    de0ty[1] = 0.0;
    de0ty[2] = 2.0 * scale * dfde[7];
    de0ty[3] = 0.0;
    de0ty[4] = 2.0 * scale * dfde[6];
    de0ty[5] = 0.0;

    // Compute the derivative with respect to u0x
    du0x[0] = scale * dfde[0];
    du0x[1] = scale * dfde[2];
    du0x[2] = 0.0;
    du0x[3] = scale * dfde[2];
    du0x[4] = scale * dfde[1];
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;

    // Compute the derivative with respect to u1x
    du1x[0] = scale * dfde[3];
    du1x[1] = scale * dfde[5];
    du1x[2] = 0.0;
    du1x[3] = scale * dfde[5];
    du1x[4] = scale * dfde[4];
    du1x[5] = 0.0;
    du1x[6] = 0.0;
    du1x[7] = 0.0;
    du1x[8] = 0.0;
  }

  static inline void evalStrainSensDeriv(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar u0x[],
      const TacsScalar u1x[], const TacsScalar dfded[], const TacsScalar u0xd[],
      const TacsScalar u1xd[], TacsScalar du0x[], TacsScalar du1x[],
      TacsScalar de0ty[], TacsScalar du0xd[], TacsScalar du1xd[],
      TacsScalar de0tyd[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = 0.0;
    de0ty[1] = 0.0;
    de0ty[2] = 2.0 * scale * dfde[7];
    de0ty[3] = 0.0;
    de0ty[4] = 2.0 * scale * dfde[6];
    de0ty[5] = 0.0;

    de0tyd[0] = 0.0;
    de0tyd[1] = 0.0;
    de0tyd[2] = 2.0 * scale * dfded[7];
    de0tyd[3] = 0.0;
    de0tyd[4] = 2.0 * scale * dfded[6];
    de0tyd[5] = 0.0;

    // Compute the derivative with respect to u0x
    du0x[0] = scale * dfde[0];
    du0x[1] = scale * dfde[2];
    du0x[2] = 0.0;
    du0x[3] = scale * dfde[2];
    du0x[4] = scale * dfde[1];
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;

    du0xd[0] = scale * dfded[0];
    du0xd[1] = scale * dfded[2];
    du0xd[2] = 0.0;
    du0xd[3] = scale * dfded[2];
    du0xd[4] = scale * dfded[1];
    du0xd[5] = 0.0;
    du0xd[6] = 0.0;
    du0xd[7] = 0.0;
    du0xd[8] = 0.0;

    // Compute the derivative with respect to u1x
    du1x[0] = scale * dfde[3];
    du1x[1] = scale * dfde[5];
    du1x[2] = 0.0;
    du1x[3] = scale * dfde[5];
    du1x[4] = scale * dfde[4];
    du1x[5] = 0.0;
    du1x[6] = 0.0;
    du1x[7] = 0.0;
    du1x[8] = 0.0;

    du1xd[0] = scale * dfded[3];
    du1xd[1] = scale * dfded[5];
    du1xd[2] = 0.0;
    du1xd[3] = scale * dfded[5];
    du1xd[4] = scale * dfded[4];
    du1xd[5] = 0.0;
    du1xd[6] = 0.0;
    du1xd[7] = 0.0;
    du1xd[8] = 0.0;
  }

  static inline void evalStrainDeriv(
      const TacsScalar u0x[], const TacsScalar u1x[], const TacsScalar e0ty[],
      const TacsScalar u0xd[], const TacsScalar u1xd[],
      const TacsScalar e0tyd[], TacsScalar e[], TacsScalar ed[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0];
    e[1] = u0x[4];
    e[2] = u0x[1] + u0x[3];

    // Compute the bending strain
    e[3] = u1x[0];
    e[4] = u1x[4];
    e[5] = u1x[1] + u1x[3];

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4];
    e[7] = 2.0 * e0ty[2];

    // Evaluate the in-plane strains from the tying strain expressions
    ed[0] = u0xd[0];
    ed[1] = u0xd[4];
    ed[2] = u0xd[1] + u0xd[3];

    // Compute the bending strain
    ed[3] = u1xd[0];
    ed[4] = u1xd[4];
    ed[5] = u1xd[1] + u1xd[3];

    // Add the components of the shear strain
    ed[6] = 2.0 * e0tyd[4];
    ed[7] = 2.0 * e0tyd[2];
  }

  static inline void evalStrainHessian(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar Cs[],
      const TacsScalar u0x[], const TacsScalar u1x[], const TacsScalar e0ty[],
      TacsScalar d2u0x[], TacsScalar d2u1x[], TacsScalar d2u0xu1x[],
      TacsScalar d2e0ty[], TacsScalar d2e0tyu0x[], TacsScalar d2e0tyu1x[]) {
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2u1x, 0, 81 * sizeof(TacsScalar));
    memset(d2u0xu1x, 0, 81 * sizeof(TacsScalar));
    memset(d2e0ty, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54 * sizeof(TacsScalar));

    // Compute the second derivatives
    TacsScalar *d2;
    d2 = d2u0x;
    d2[0] = scale * A[0];
    d2[4] = scale * A[1];
    d2[1] = scale * A[2];
    d2[3] = scale * A[2];

    d2 = &d2u0x[4 * 9];
    d2[0] = scale * A[1];
    d2[4] = scale * A[3];
    d2[1] = scale * A[4];
    d2[3] = scale * A[4];

    d2 = &d2u0x[9];
    d2[0] = scale * A[2];
    d2[4] = scale * A[4];
    d2[1] = scale * A[5];
    d2[3] = scale * A[5];

    d2 = &d2u0x[3 * 9];
    d2[0] = scale * A[2];
    d2[4] = scale * A[4];
    d2[1] = scale * A[5];
    d2[3] = scale * A[5];

    // e[3] = u1x[0];
    // e[4] = u1x[4];
    // e[5] = u1x[1] + u1x[3];
    d2 = d2u1x;
    d2[0] = scale * D[0];
    d2[4] = scale * D[1];
    d2[1] = scale * D[2];
    d2[3] = scale * D[2];

    d2 = &d2u1x[4 * 9];
    d2[0] = scale * D[1];
    d2[4] = scale * D[3];
    d2[1] = scale * D[4];
    d2[3] = scale * D[4];

    d2 = &d2u1x[9];
    d2[0] = scale * D[2];
    d2[4] = scale * D[4];
    d2[1] = scale * D[5];
    d2[3] = scale * D[5];

    d2 = &d2u1x[3 * 9];
    d2[0] = scale * D[2];
    d2[4] = scale * D[4];
    d2[1] = scale * D[5];
    d2[3] = scale * D[5];

    d2 = d2u0xu1x;
    d2[0] = scale * B[0];
    d2[4] = scale * B[1];
    d2[1] = scale * B[2];
    d2[3] = scale * B[2];

    d2 = &d2u0xu1x[4 * 9];
    d2[0] = scale * B[1];
    d2[4] = scale * B[3];
    d2[1] = scale * B[4];
    d2[3] = scale * B[4];

    d2 = &d2u0xu1x[9];
    d2[0] = scale * B[2];
    d2[4] = scale * B[4];
    d2[1] = scale * B[5];
    d2[3] = scale * B[5];

    d2 = &d2u0xu1x[3 * 9];
    d2[0] = scale * B[2];
    d2[4] = scale * B[4];
    d2[1] = scale * B[5];
    d2[3] = scale * B[5];

    // e[6] = 2.0*e0ty[4];
    // e[7] = 2.0*e0ty[2];
    d2 = &d2e0ty[4 * 6];
    d2[4] = 4.0 * scale * As[0];
    d2[2] = 4.0 * scale * As[1];

    d2 = &d2e0ty[2 * 6];
    d2[4] = 4.0 * scale * As[1];
    d2[2] = 4.0 * scale * As[2];
  }

  static inline void evalStrainHessianDeriv(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar Cs[],
      const TacsScalar u0x[], const TacsScalar u1x[], const TacsScalar e0ty[],
      const TacsScalar dfded[], const TacsScalar u0xd[],
      const TacsScalar u1xd[], const TacsScalar e0tyd[], TacsScalar d2u0x[],
      TacsScalar d2u1x[], TacsScalar d2u0xu1x[], TacsScalar d2e0ty[],
      TacsScalar d2e0tyu0x[], TacsScalar d2e0tyu1x[], TacsScalar d2u0xd[],
      TacsScalar d2u1xd[], TacsScalar d2u0xu1xd[], TacsScalar d2e0tyd[],
      TacsScalar d2e0tyu0xd[], TacsScalar d2e0tyu1xd[]) {
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2u1x, 0, 81 * sizeof(TacsScalar));
    memset(d2u0xu1x, 0, 81 * sizeof(TacsScalar));
    memset(d2e0ty, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54 * sizeof(TacsScalar));
    memset(d2u0xd, 0, 81 * sizeof(TacsScalar));
    memset(d2u1xd, 0, 81 * sizeof(TacsScalar));
    memset(d2u0xu1xd, 0, 81 * sizeof(TacsScalar));
    memset(d2e0tyd, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0xd, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1xd, 0, 54 * sizeof(TacsScalar));

    // Compute the second derivatives
    TacsScalar *d2;
    d2 = d2u0x;
    d2[0] = scale * A[0];
    d2[4] = scale * A[1];
    d2[1] = scale * A[2];
    d2[3] = scale * A[2];

    d2 = &d2u0x[4 * 9];
    d2[0] = scale * A[1];
    d2[4] = scale * A[3];
    d2[1] = scale * A[4];
    d2[3] = scale * A[4];

    d2 = &d2u0x[9];
    d2[0] = scale * A[2];
    d2[4] = scale * A[4];
    d2[1] = scale * A[5];
    d2[3] = scale * A[5];

    d2 = &d2u0x[3 * 9];
    d2[0] = scale * A[2];
    d2[4] = scale * A[4];
    d2[1] = scale * A[5];
    d2[3] = scale * A[5];

    // e[3] = u1x[0];
    // e[4] = u1x[4];
    // e[5] = u1x[1] + u1x[3];
    d2 = d2u1x;
    d2[0] = scale * D[0];
    d2[4] = scale * D[1];
    d2[1] = scale * D[2];
    d2[3] = scale * D[2];

    d2 = &d2u1x[4 * 9];
    d2[0] = scale * D[1];
    d2[4] = scale * D[3];
    d2[1] = scale * D[4];
    d2[3] = scale * D[4];

    d2 = &d2u1x[9];
    d2[0] = scale * D[2];
    d2[4] = scale * D[4];
    d2[1] = scale * D[5];
    d2[3] = scale * D[5];

    d2 = &d2u1x[3 * 9];
    d2[0] = scale * D[2];
    d2[4] = scale * D[4];
    d2[1] = scale * D[5];
    d2[3] = scale * D[5];

    d2 = d2u0xu1x;
    d2[0] = scale * B[0];
    d2[4] = scale * B[1];
    d2[1] = scale * B[2];
    d2[3] = scale * B[2];

    d2 = &d2u0xu1x[4 * 9];
    d2[0] = scale * B[1];
    d2[4] = scale * B[3];
    d2[1] = scale * B[4];
    d2[3] = scale * B[4];

    d2 = &d2u0xu1x[9];
    d2[0] = scale * B[2];
    d2[4] = scale * B[4];
    d2[1] = scale * B[5];
    d2[3] = scale * B[5];

    d2 = &d2u0xu1x[3 * 9];
    d2[0] = scale * B[2];
    d2[4] = scale * B[4];
    d2[1] = scale * B[5];
    d2[3] = scale * B[5];

    // e[6] = 2.0*e0ty[4];
    // e[7] = 2.0*e0ty[2];
    d2 = &d2e0ty[4 * 6];
    d2[4] = 4.0 * scale * As[0];
    d2[2] = 4.0 * scale * As[1];

    d2 = &d2e0ty[2 * 6];
    d2[4] = 4.0 * scale * As[1];
    d2[2] = 4.0 * scale * As[2];
  }
};

class TACSShellInplaneNonlinearModel {
 public:
  /**
    Compute the tensorial components of the tying strain

    G = 0.5*(X,eta^{T}*U,eta + U,eta^{T}*X,eta)

    The derivative with respect to the frame gives

    X,eta = [X,xi ; n]

    The derivative with respect to the displacements gives

    u,eta = [u,xi ; d]

    @param Xxi Derivatives of the node locations with respect to xi
    @param n The interpolated frame normal
    @param Uxi Derivatives of the displacements with respect to xi
    @param d The interpolated director field
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrain(const TacsScalar Xpts[], const TacsScalar fn[],
                                 const TacsScalar vars[], const TacsScalar d[],
                                 TacsScalar ety[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6];
      TacsScalar d0[3], n0[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, fn, n0);

      ety[index] = 0.0;
      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        ety[index] =
            0.5 * (Xxi[1] * d0[0] + Xxi[3] * d0[1] + Xxi[5] * d0[2] +
                   (n0[0] + d0[0]) * Uxi[1] + (n0[1] + d0[1]) * Uxi[3] +
                   (n0[2] + d0[2]) * Uxi[5]);
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        ety[index] =
            0.5 * (Xxi[0] * d0[0] + Xxi[2] * d0[1] + Xxi[4] * d0[2] +
                   (n0[0] + d0[0]) * Uxi[0] + (n0[1] + d0[1]) * Uxi[2] +
                   (n0[2] + d0[2]) * Uxi[4]);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose(
      const TacsScalar Xpts[], const TacsScalar fn[], const TacsScalar vars[],
      const TacsScalar d[], const TacsScalar dety[], TacsScalar res[],
      TacsScalar dd[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6], dUxi[6];
      TacsScalar n0[3], d0[3], dd0[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, fn, n0);

      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi[0] = 0.0;
        dUxi[1] = 0.5 * dety[index] * (n0[0] + d0[0]);
        dUxi[2] = 0.0;
        dUxi[3] = 0.5 * dety[index] * (n0[1] + d0[1]);
        dUxi[4] = 0.0;
        dUxi[5] = 0.5 * dety[index] * (n0[2] + d0[2]);

        dd0[0] = 0.5 * dety[index] * (Xxi[1] + Uxi[1]);
        dd0[1] = 0.5 * dety[index] * (Xxi[3] + Uxi[3]);
        dd0[2] = 0.5 * dety[index] * (Xxi[5] + Uxi[5]);
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi[0] = 0.5 * dety[index] * (n0[0] + d0[0]);
        dUxi[1] = 0.0;
        dUxi[2] = 0.5 * dety[index] * (n0[1] + d0[1]);
        dUxi[3] = 0.0;
        dUxi[4] = 0.5 * dety[index] * (n0[2] + d0[2]);
        dUxi[5] = 0.0;

        dd0[0] = 0.5 * dety[index] * (Xxi[0] + Uxi[0]);
        dd0[1] = 0.5 * dety[index] * (Xxi[2] + Uxi[2]);
        dd0[2] = 0.5 * dety[index] * (Xxi[4] + Uxi[4]);
      }
      basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);

      if (res) {
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi,
                                                                       res);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian(
      const TacsScalar alpha, const TacsScalar Xpts[], const TacsScalar fn[],
      const TacsScalar vars[], const TacsScalar d[], const TacsScalar dety[],
      const TacsScalar d2ety[], const TacsScalar d2etyu[],
      const TacsScalar d2etyd[], TacsScalar mat[], TacsScalar d2d[],
      TacsScalar d2du[]) {
    // Initialize the data
    TacsScalar n0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar d0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Uxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity, *d0 = d0ty, *Uxi = Uxity;

    // Pre-compute terms needed at each tying point
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFields<3, 3>(pt, d, d0);

      n0 += 3;
      Xxi += 6;
      d0 += 3;
      Uxi += 6;
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity, *d01 = d0ty, *Uxi1 = Uxity;
    for (int i1 = 0; i1 < basis::NUM_TYING_POINTS;
         i1++, n01 += 3, Xxi1 += 6, d01 += 3, Uxi1 += 6) {
      // Get the field index
      const TacsShellTyingStrainComponent f1 = basis::getTyingField(i1);

      // Get the tying point parametric location
      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar du2[3 * basis::NUM_NODES], dd2[3 * basis::NUM_NODES];
      memset(du2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      TacsScalar *n02 = n0ty, *Xxi2 = Xxity, *d02 = d0ty, *Uxi2 = Uxity;
      for (int i2 = 0; i2 < basis::NUM_TYING_POINTS;
           i2++, n02 += 3, Xxi2 += 6, d02 += 3, Uxi2 += 6) {
        // Get the field index
        const TacsShellTyingStrainComponent f2 = basis::getTyingField(i2);

        // Get the tying point parametric location
        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        const TacsScalar value = d2ety[basis::NUM_TYING_POINTS * i1 + i2];

        TacsScalar dUxi2[6], dd02[3];
        if (f2 == TACS_SHELL_G23_COMPONENT) {
          // Compute g23 = e2^{T}*G*e3
          dUxi2[0] = 0.0;
          dUxi2[1] = 0.5 * value * (n02[0] + d02[0]);
          dUxi2[2] = 0.0;
          dUxi2[3] = 0.5 * value * (n02[1] + d02[1]);
          dUxi2[4] = 0.0;
          dUxi2[5] = 0.5 * value * (n02[2] + d02[2]);

          dd02[0] = 0.5 * value * (Xxi2[1] + Uxi2[1]);
          dd02[1] = 0.5 * value * (Xxi2[3] + Uxi2[3]);
          dd02[2] = 0.5 * value * (Xxi2[5] + Uxi2[5]);
        } else if (f2 == TACS_SHELL_G13_COMPONENT) {
          // Compute g13 = e1^{T}*G*e3
          dUxi2[0] = 0.5 * value * (n02[0] + d02[0]);
          dUxi2[1] = 0.0;
          dUxi2[2] = 0.5 * value * (n02[1] + d02[1]);
          dUxi2[3] = 0.0;
          dUxi2[4] = 0.5 * value * (n02[2] + d02[2]);
          dUxi2[5] = 0.0;

          dd02[0] = 0.5 * value * (Xxi2[0] + Uxi2[0]);
          dd02[1] = 0.5 * value * (Xxi2[2] + Uxi2[2]);
          dd02[2] = 0.5 * value * (Xxi2[4] + Uxi2[4]);
        }

        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
      }

      TacsScalar du1[3 * basis::NUM_NODES];
      TacsScalar dd1[3 * basis::NUM_NODES];
      memset(du1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      // Store the the derivative information for the first point
      TacsScalar dUxi1[6], d2Uxi[36];
      memset(d2Uxi, 0, 36 * sizeof(TacsScalar));

      TacsScalar dd01[3], d2dUxi[18];
      memset(d2dUxi, 0, 18 * sizeof(TacsScalar));

      if (f1 == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi1[0] = 0.0;
        dUxi1[1] = 0.5 * (n01[0] + d01[0]);
        dUxi1[2] = 0.0;
        dUxi1[3] = 0.5 * (n01[1] + d01[1]);
        dUxi1[4] = 0.0;
        dUxi1[5] = 0.5 * (n01[2] + d01[2]);

        dd01[0] = 0.5 * (Xxi1[1] + Uxi1[1]);
        dd01[1] = 0.5 * (Xxi1[3] + Uxi1[3]);
        dd01[2] = 0.5 * (Xxi1[5] + Uxi1[5]);

        d2dUxi[1] = 0.5 * alpha * dety[i1];
        d2dUxi[9] = 0.5 * alpha * dety[i1];
        d2dUxi[17] = 0.5 * alpha * dety[i1];
      } else if (f1 == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi1[0] = 0.5 * (n01[0] + d01[0]);
        dUxi1[1] = 0.0;
        dUxi1[2] = 0.5 * (n01[1] + d01[1]);
        dUxi1[3] = 0.0;
        dUxi1[4] = 0.5 * (n01[2] + d01[2]);
        dUxi1[5] = 0.0;

        dd01[0] = 0.5 * (Xxi1[0] + Uxi1[0]);
        dd01[1] = 0.5 * (Xxi1[2] + Uxi1[2]);
        dd01[2] = 0.5 * (Xxi1[4] + Uxi1[4]);

        d2dUxi[0] = 0.5 * alpha * dety[i1];
        d2dUxi[8] = 0.5 * alpha * dety[i1];
        d2dUxi[16] = 0.5 * alpha * dety[i1];
      }

      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

      basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt1, d2dUxi,
                                                                 NULL, d2du);
      basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3,
                                                3>(pt1, d2Uxi, mat);

      const TacsScalar *etd = &d2etyd[3 * basis::NUM_NODES * i1];
      const TacsScalar *etu = &d2etyu[3 * basis::NUM_NODES * i1];
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2d[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * dd2[j] + dd1[i] * etd[j] + etd[i] * dd1[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * du2[j] + dd1[i] * etu[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] += etd[i] * du1[j];
        }
      }

      const int nvars = vars_per_node * basis::NUM_NODES;
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        int ii = vars_per_node * (i / 3) + i % 3;
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          int jj = vars_per_node * (j / 3) + j % 3;
          mat[nvars * ii + jj] +=
              du1[i] * du2[j] + du1[i] * etu[j] + etu[i] * du1[j];
        }
      }
    }
  }

  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv(
      const TacsScalar Xpts[], const TacsScalar fn[], const TacsScalar vars[],
      const TacsScalar d[], const TacsScalar varsd[], const TacsScalar dd[],
      TacsScalar ety[], TacsScalar etyd[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsShellTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[6], Xxi[6], Uxid[6];
      TacsScalar n0[3], d0[3], d0d[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, Uxid);
      basis::template interpFields<3, 3>(pt, d, d0);
      basis::template interpFields<3, 3>(pt, dd, d0d);
      basis::template interpFields<3, 3>(pt, fn, n0);

      ety[index] = 0.0;
      etyd[index] = 0.0;

      if (field == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        ety[index] =
            0.5 * (Xxi[1] * d0[0] + Xxi[3] * d0[1] + Xxi[5] * d0[2] +
                   (n0[0] + d0[0]) * Uxi[1] + (n0[1] + d0[1]) * Uxi[3] +
                   (n0[2] + d0[2]) * Uxi[5]);
        etyd[index] =
            0.5 * (Xxi[1] * d0d[0] + Xxi[3] * d0d[1] + Xxi[5] * d0d[2] +
                   (n0[0] + d0[0]) * Uxid[1] + d0d[0] * Uxi[1] +
                   (n0[1] + d0[1]) * Uxid[3] + d0d[1] * Uxi[3] +
                   (n0[2] + d0[2]) * Uxid[5] + d0d[2] * Uxi[5]);
      } else if (field == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        ety[index] =
            0.5 * (Xxi[0] * d0[0] + Xxi[2] * d0[1] + Xxi[4] * d0[2] +
                   (n0[0] + d0[0]) * Uxi[0] + (n0[1] + d0[1]) * Uxi[2] +
                   (n0[2] + d0[2]) * Uxi[4]);
        etyd[index] =
            0.5 * (Xxi[0] * d0d[0] + Xxi[2] * d0d[1] + Xxi[4] * d0d[2] +
                   (n0[0] + d0[0]) * Uxid[0] + d0d[0] * Uxi[0] +
                   (n0[1] + d0[1]) * Uxid[2] + d0d[1] * Uxi[2] +
                   (n0[2] + d0[2]) * Uxid[4] + d0d[2] * Uxi[4]);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessianDeriv(
      const TacsScalar alpha, const TacsScalar Xpts[], const TacsScalar fn[],
      const TacsScalar vars[], const TacsScalar d[], const TacsScalar dety[],
      const TacsScalar d2ety[], const TacsScalar d2etyu[],
      const TacsScalar d2etyd[], const TacsScalar vars_d[],
      const TacsScalar d_d[], const TacsScalar dety_d[],
      const TacsScalar d2ety_d[], const TacsScalar d2etyu_d[],
      const TacsScalar d2etyd_d[], TacsScalar mat[], TacsScalar d2d[],
      TacsScalar d2du[], TacsScalar mat_d[], TacsScalar d2d_d[],
      TacsScalar d2du_d[]) {
    // Initialize the data
    TacsScalar n0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar d0ty[3 * basis::NUM_TYING_POINTS];
    TacsScalar Uxity[6 * basis::NUM_TYING_POINTS];
    TacsScalar d0tyd[3 * basis::NUM_TYING_POINTS];
    TacsScalar Uxityd[6 * basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity, *d0 = d0ty, *Uxi = Uxity;
    TacsScalar *d0d = d0tyd, *Uxid = Uxityd;

    // Pre-compute terms needed at each tying point
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFields<3, 3>(pt, d, d0);

      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars_d, Uxid);
      basis::template interpFields<3, 3>(pt, d_d, d0d);

      n0 += 3;
      Xxi += 6;
      d0 += 3;
      Uxi += 6;

      d0d += 3;
      Uxid += 6;
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity, *d01 = d0ty, *Uxi1 = Uxity,
               *d01d = d0tyd, *Uxi1d = Uxityd;
    for (int i1 = 0; i1 < basis::NUM_TYING_POINTS; i1++, n01 += 3, Xxi1 += 6,
             d01 += 3, Uxi1 += 6, d01d += 3, Uxi1d += 6) {
      // Get the field index
      const TacsShellTyingStrainComponent f1 = basis::getTyingField(i1);

      // Get the tying point parametric location
      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar du2[3 * basis::NUM_NODES], dd2[3 * basis::NUM_NODES];
      memset(du2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      TacsScalar du2d[3 * basis::NUM_NODES], dd2d[3 * basis::NUM_NODES];
      memset(du2d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd2d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      TacsScalar *n02 = n0ty, *Xxi2 = Xxity, *d02 = d0ty, *Uxi2 = Uxity,
                 *d02d = d0tyd, *Uxi2d = Uxityd;
      for (int i2 = 0; i2 < basis::NUM_TYING_POINTS; i2++, n02 += 3, Xxi2 += 6,
               d02 += 3, Uxi2 += 6, d02d += 3, Uxi2d += 6) {
        // Get the field index
        const TacsShellTyingStrainComponent f2 = basis::getTyingField(i2);

        // Get the tying point parametric location
        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        const TacsScalar value = d2ety[basis::NUM_TYING_POINTS * i1 + i2];
        const TacsScalar valued = d2ety_d[basis::NUM_TYING_POINTS * i1 + i2];

        TacsScalar dUxi2[6], dd02[3];
        TacsScalar dUxi2d[6], dd02d[3];
        if (f2 == TACS_SHELL_G23_COMPONENT) {
          // Compute g23 = e2^{T}*G*e3
          dUxi2[0] = 0.0;
          dUxi2[1] = 0.5 * value * (n02[0] + d02[0]);
          dUxi2[2] = 0.0;
          dUxi2[3] = 0.5 * value * (n02[1] + d02[1]);
          dUxi2[4] = 0.0;
          dUxi2[5] = 0.5 * value * (n02[2] + d02[2]);

          dUxi2d[0] = 0.0;
          dUxi2d[1] =
              0.5 * value * (d02d[0]) + 0.5 * valued * (n02[0] + d02[0]);
          dUxi2d[2] = 0.0;
          dUxi2d[3] =
              0.5 * value * (d02d[1]) + 0.5 * valued * (n02[1] + d02[1]);
          dUxi2d[4] = 0.0;
          dUxi2d[5] =
              0.5 * value * (d02d[2]) + 0.5 * valued * (n02[2] + d02[2]);

          dd02[0] = 0.5 * value * (Xxi2[1] + Uxi2[1]);
          dd02[1] = 0.5 * value * (Xxi2[3] + Uxi2[3]);
          dd02[2] = 0.5 * value * (Xxi2[5] + Uxi2[5]);

          dd02d[0] =
              0.5 * value * (Uxi2d[1]) + 0.5 * valued * (Xxi2[1] + Uxi2[1]);
          dd02d[1] =
              0.5 * value * (Uxi2d[3]) + 0.5 * valued * (Xxi2[3] + Uxi2[3]);
          dd02d[2] =
              0.5 * value * (Uxi2d[5]) + 0.5 * valued * (Xxi2[5] + Uxi2[5]);
        } else if (f2 == TACS_SHELL_G13_COMPONENT) {
          // Compute g13 = e1^{T}*G*e3
          dUxi2[0] = 0.5 * value * (n02[0] + d02[0]);
          dUxi2[1] = 0.0;
          dUxi2[2] = 0.5 * value * (n02[1] + d02[1]);
          dUxi2[3] = 0.0;
          dUxi2[4] = 0.5 * value * (n02[2] + d02[2]);
          dUxi2[5] = 0.0;

          dUxi2d[0] =
              0.5 * value * (d02d[0]) + 0.5 * valued * (n02[0] + d02[0]);
          dUxi2d[1] = 0.0;
          dUxi2d[2] =
              0.5 * value * (d02d[1]) + 0.5 * valued * (n02[1] + d02[1]);
          dUxi2d[3] = 0.0;
          dUxi2d[4] =
              0.5 * value * (d02d[2]) + 0.5 * valued * (n02[2] + d02d[2]);
          dUxi2d[5] = 0.0;

          dd02[0] = 0.5 * value * (Xxi2[0] + Uxi2[0]);
          dd02[1] = 0.5 * value * (Xxi2[2] + Uxi2[2]);
          dd02[2] = 0.5 * value * (Xxi2[4] + Uxi2[4]);

          dd02d[0] =
              0.5 * value * (Uxi2d[0]) + 0.5 * valued * (Xxi2[0] + Uxi2[0]);
          dd02d[1] =
              0.5 * value * (Uxi2d[2]) + 0.5 * valued * (Xxi2[2] + Uxi2[2]);
          dd02d[2] =
              0.5 * value * (Uxi2d[4]) + 0.5 * valued * (Xxi2[4] + Uxi2[4]);
        }

        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02d, dd2d);
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2d, du2d);
      }

      TacsScalar du1[3 * basis::NUM_NODES];
      TacsScalar dd1[3 * basis::NUM_NODES];
      memset(du1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      TacsScalar du1d[3 * basis::NUM_NODES];
      TacsScalar dd1d[3 * basis::NUM_NODES];
      memset(du1d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
      memset(dd1d, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

      // Store the the derivative information for the first point
      TacsScalar dUxi1[6], d2Uxi[36];
      memset(d2Uxi, 0, 36 * sizeof(TacsScalar));
      TacsScalar dUxi1d[6], d2Uxid[36];
      memset(d2Uxid, 0, 36 * sizeof(TacsScalar));

      TacsScalar dd01[3], d2dUxi[18];
      memset(d2dUxi, 0, 18 * sizeof(TacsScalar));
      TacsScalar dd01d[3], d2dUxid[18];
      memset(d2dUxid, 0, 18 * sizeof(TacsScalar));

      if (f1 == TACS_SHELL_G23_COMPONENT) {
        // Compute g23 = e2^{T}*G*e3
        dUxi1[0] = 0.0;
        dUxi1[1] = 0.5 * (n01[0] + d01[0]);
        dUxi1[2] = 0.0;
        dUxi1[3] = 0.5 * (n01[1] + d01[1]);
        dUxi1[4] = 0.0;
        dUxi1[5] = 0.5 * (n01[2] + d01[2]);

        dUxi1d[0] = 0.0;
        dUxi1d[1] = 0.5 * (d01d[0]);
        dUxi1d[2] = 0.0;
        dUxi1d[3] = 0.5 * (d01d[1]);
        dUxi1d[4] = 0.0;
        dUxi1d[5] = 0.5 * (d01d[2]);

        dd01[0] = 0.5 * (Xxi1[1] + Uxi1[1]);
        dd01[1] = 0.5 * (Xxi1[3] + Uxi1[3]);
        dd01[2] = 0.5 * (Xxi1[5] + Uxi1[5]);

        dd01d[0] = 0.5 * (Uxi1d[1]);
        dd01d[1] = 0.5 * (Uxi1d[3]);
        dd01d[2] = 0.5 * (Uxi1d[5]);

        d2dUxi[1] = 0.5 * alpha * dety[i1];
        d2dUxi[9] = 0.5 * alpha * dety[i1];
        d2dUxi[17] = 0.5 * alpha * dety[i1];

        d2dUxid[1] = 0.5 * alpha * dety_d[i1];
        d2dUxid[9] = 0.5 * alpha * dety_d[i1];
        d2dUxid[17] = 0.5 * alpha * dety_d[i1];
      } else if (f1 == TACS_SHELL_G13_COMPONENT) {
        // Compute g13 = e1^{T}*G*e3
        dUxi1[0] = 0.5 * (n01[0] + d01[0]);
        dUxi1[1] = 0.0;
        dUxi1[2] = 0.5 * (n01[1] + d01[1]);
        dUxi1[3] = 0.0;
        dUxi1[4] = 0.5 * (n01[2] + d01[2]);
        dUxi1[5] = 0.0;

        dUxi1d[0] = 0.5 * (d01d[0]);
        dUxi1d[1] = 0.0;
        dUxi1d[2] = 0.5 * (d01d[1]);
        dUxi1d[3] = 0.0;
        dUxi1d[4] = 0.5 * (d01d[2]);
        dUxi1d[5] = 0.0;

        dd01[0] = 0.5 * (Xxi1[0] + Uxi1[0]);
        dd01[1] = 0.5 * (Xxi1[2] + Uxi1[2]);
        dd01[2] = 0.5 * (Xxi1[4] + Uxi1[4]);

        dd01d[0] = 0.5 * (Uxi1d[0]);
        dd01d[1] = 0.5 * (Uxi1d[2]);
        dd01d[2] = 0.5 * (Uxi1d[4]);

        d2dUxi[0] = 0.5 * alpha * dety[i1];
        d2dUxi[8] = 0.5 * alpha * dety[i1];
        d2dUxi[16] = 0.5 * alpha * dety[i1];

        d2dUxid[0] = 0.5 * alpha * dety_d[i1];
        d2dUxid[8] = 0.5 * alpha * dety_d[i1];
        d2dUxid[16] = 0.5 * alpha * dety_d[i1];
      }

      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);
      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01d, dd1d);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1d, du1d);

      basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt1, d2dUxi,
                                                                 NULL, d2du);
      if (mat) {
        basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node,
                                                  3, 3>(pt1, d2Uxi, mat);
      }
      basis::template addInterpGradMixedOuterProduct<3, 3, 3, 3>(pt1, d2dUxid,
                                                                 NULL, d2du_d);
      basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 3,
                                                3>(pt1, d2Uxid, mat_d);

      const TacsScalar *etd = &d2etyd[3 * basis::NUM_NODES * i1];
      const TacsScalar *etu = &d2etyu[3 * basis::NUM_NODES * i1];
      const TacsScalar *etdd = &d2etyd_d[3 * basis::NUM_NODES * i1];
      const TacsScalar *etud = &d2etyu_d[3 * basis::NUM_NODES * i1];
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2d[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * dd2[j] + dd1[i] * etd[j] + etd[i] * dd1[j];
          d2d_d[3 * basis::NUM_NODES * i + j] +=
              dd1d[i] * dd2[j] + dd1d[i] * etd[j] + etdd[i] * dd1[j] +
              dd1[i] * dd2d[j] + dd1[i] * etdd[j] + etd[i] * dd1d[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] +=
              dd1[i] * du2[j] + dd1[i] * etu[j];
          d2du_d[3 * basis::NUM_NODES * i + j] +=
              dd1d[i] * du2[j] + dd1d[i] * etu[j] + dd1[i] * du2d[j] +
              dd1[i] * etud[j];
        }
      }

      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          d2du[3 * basis::NUM_NODES * i + j] += etd[i] * du1[j];
          d2du_d[3 * basis::NUM_NODES * i + j] +=
              etdd[i] * du1[j] + etd[i] * du1d[j];
        }
      }

      const int nvars = vars_per_node * basis::NUM_NODES;
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        int ii = vars_per_node * (i / 3) + i % 3;
        for (int j = 0; j < 3 * basis::NUM_NODES; j++) {
          int jj = vars_per_node * (j / 3) + j % 3;
          if (mat) {
            mat[nvars * ii + jj] +=
                du1[i] * du2[j] + du1[i] * etu[j] + etu[i] * du1[j];
          }
          mat_d[nvars * ii + jj] += du1d[i] * du2[j] + du1d[i] * etu[j] +
                                    etud[i] * du1[j] + du1[i] * du2d[j] +
                                    du1[i] * etud[j] + etu[i] * du1d[j];
        }
      }
    }
  }

  /*
    Evaluate the strain as a function of the displacement derivatives
    and interpolated strain from the tensorial components
  */
  static void evalStrain(const TacsScalar u0x[], const TacsScalar u1x[],
                         const TacsScalar e0ty[], TacsScalar e[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0] + 0.5 * (u0x[0] * u0x[0] + u0x[3] * u0x[3] + u0x[6] * u0x[6]);
    e[1] = u0x[4] + 0.5 * (u0x[1] * u0x[1] + u0x[4] * u0x[4] + u0x[7] * u0x[7]);
    e[2] =
        u0x[1] + u0x[3] + (u0x[0] * u0x[1] + u0x[3] * u0x[4] + u0x[6] * u0x[7]);

    // Compute the bending strain
    e[3] = u1x[0] + (u0x[0] * u1x[0] + u0x[3] * u1x[3] + u0x[6] * u1x[6]);
    e[4] = u1x[4] + (u0x[1] * u1x[1] + u0x[4] * u1x[4] + u0x[7] * u1x[7]);
    e[5] = u1x[1] + u1x[3] +
           (u0x[0] * u1x[1] + u0x[3] * u1x[4] + u0x[6] * u1x[7] +
            u1x[0] * u0x[1] + u1x[3] * u0x[4] + u1x[6] * u0x[7]);

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4];
    e[7] = 2.0 * e0ty[2];
  }

  /**
    Evaluate the derivative of the strain
  */
  static void evalStrainSens(const TacsScalar scale, const TacsScalar dfde[],
                             const TacsScalar u0x[], const TacsScalar u1x[],
                             TacsScalar du0x[], TacsScalar du1x[],
                             TacsScalar de0ty[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = 0.0;
    de0ty[1] = 0.0;
    de0ty[2] = 2.0 * scale * dfde[7];
    de0ty[3] = 0.0;
    de0ty[4] = 2.0 * scale * dfde[6];
    de0ty[5] = 0.0;

    // Derivative with respect to u0x
    du0x[0] = scale * (dfde[0] * (u0x[0] + 1.0) + dfde[2] * u0x[1] +
                       dfde[3] * u1x[0] + dfde[5] * u1x[1]);
    du0x[1] = scale * (dfde[1] * u0x[1] + dfde[2] * (u0x[0] + 1.0) +
                       dfde[4] * u1x[1] + dfde[5] * u1x[0]);
    du0x[2] = 0.0;
    du0x[3] = scale * (dfde[0] * u0x[3] + dfde[2] * (u0x[4] + 1.0) +
                       dfde[3] * u1x[3] + dfde[5] * u1x[4]);
    du0x[4] = scale * (dfde[1] * (u0x[4] + 1.0) + dfde[2] * u0x[3] +
                       dfde[4] * u1x[4] + dfde[5] * u1x[3]);
    du0x[5] = 0.0;
    du0x[6] = scale * (dfde[0] * u0x[6] + dfde[2] * u0x[7] + dfde[3] * u1x[6] +
                       dfde[5] * u1x[7]);
    du0x[7] = scale * (dfde[1] * u0x[7] + dfde[2] * u0x[6] + dfde[4] * u1x[7] +
                       dfde[5] * u1x[6]);
    du0x[8] = 0.0;

    du1x[0] = scale * (dfde[3] * (u0x[0] + 1.0) + dfde[5] * u0x[1]);
    du1x[1] = scale * (dfde[4] * u0x[1] + dfde[5] * (u0x[0] + 1.0));
    du1x[2] = 0.0;
    du1x[3] = scale * (dfde[3] * u0x[3] + dfde[5] * (u0x[4] + 1.0));
    du1x[4] = scale * (dfde[4] * (u0x[4] + 1.0) + dfde[5] * u0x[3]);
    du1x[5] = 0.0;
    du1x[6] = scale * (dfde[3] * u0x[6] + dfde[5] * u0x[7]);
    du1x[7] = scale * (dfde[4] * u0x[7] + dfde[5] * u0x[6]);
    du1x[8] = 0.0;
  }

  static inline void evalStrainSensDeriv(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar u0x[],
      const TacsScalar u1x[], const TacsScalar u0xd[], const TacsScalar u1xd[],
      const TacsScalar dfded[], TacsScalar du0x[], TacsScalar du1x[],
      TacsScalar de0ty[], TacsScalar du0xd[], TacsScalar du1xd[],
      TacsScalar de0tyd[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = 0.0;
    de0ty[1] = 0.0;
    de0ty[2] = 2.0 * scale * dfde[7];
    de0ty[3] = 0.0;
    de0ty[4] = 2.0 * scale * dfde[6];
    de0ty[5] = 0.0;

    de0tyd[0] = 0.0;
    de0tyd[1] = 0.0;
    de0tyd[2] = 2.0 * scale * dfded[7];
    de0tyd[3] = 0.0;
    de0tyd[4] = 2.0 * scale * dfded[6];
    de0tyd[5] = 0.0;

    // Derivative with respect to u0x
    du0x[0] = scale * (dfde[0] * (u0x[0] + 1.0) + dfde[2] * u0x[1] +
                       dfde[3] * u1x[0] + dfde[5] * u1x[1]);
    du0x[1] = scale * (dfde[1] * u0x[1] + dfde[2] * (u0x[0] + 1.0) +
                       dfde[4] * u1x[1] + dfde[5] * u1x[0]);
    du0x[2] = 0.0;
    du0x[3] = scale * (dfde[0] * u0x[3] + dfde[2] * (u0x[4] + 1.0) +
                       dfde[3] * u1x[3] + dfde[5] * u1x[4]);
    du0x[4] = scale * (dfde[1] * (u0x[4] + 1.0) + dfde[2] * u0x[3] +
                       dfde[4] * u1x[4] + dfde[5] * u1x[3]);
    du0x[5] = 0.0;
    du0x[6] = scale * (dfde[0] * u0x[6] + dfde[2] * u0x[7] + dfde[3] * u1x[6] +
                       dfde[5] * u1x[7]);
    du0x[7] = scale * (dfde[1] * u0x[7] + dfde[2] * u0x[6] + dfde[4] * u1x[7] +
                       dfde[5] * u1x[6]);
    du0x[8] = 0.0;

    du0xd[0] =
        scale * (dfde[0] * u0xd[0] + dfde[2] * u0xd[1] + dfde[3] * u1xd[0] +
                 dfde[5] * u1xd[1] + dfded[0] * (u0x[0] + 1.0) +
                 dfded[2] * u0x[1] + dfded[3] * u1x[0] + dfded[5] * u1x[1]);
    du0xd[1] = scale * (dfde[1] * u0xd[1] + dfde[2] * u0xd[0] +
                        dfde[4] * u1xd[1] + dfde[5] * u1xd[0] +
                        dfded[1] * u0x[1] + dfded[2] * (u0x[0] + 1.0) +
                        dfded[4] * u1x[1] + dfded[5] * u1x[0]);
    du0xd[2] = 0.0;
    du0xd[3] = scale * (dfde[0] * u0xd[3] + dfde[2] * u0xd[4] +
                        dfde[3] * u1xd[3] + dfde[5] * u1xd[4] +
                        dfded[0] * u0x[3] + dfded[2] * (u0x[4] + 1.0) +
                        dfded[3] * u1x[3] + dfded[5] * u1x[4]);
    du0xd[4] =
        scale * (dfde[1] * u0xd[4] + dfde[2] * u0xd[3] + dfde[4] * u1xd[4] +
                 dfde[5] * u1xd[3] + dfded[1] * (u0x[4] + 1.0) +
                 dfded[2] * u0x[3] + dfded[4] * u1x[4] + dfded[5] * u1x[3]);
    du0xd[5] = 0.0;
    du0xd[6] =
        scale * (dfde[0] * u0xd[6] + dfde[2] * u0xd[7] + dfde[3] * u1xd[6] +
                 dfde[5] * u1xd[7] + dfded[0] * u0x[6] + dfded[2] * u0x[7] +
                 dfded[3] * u1x[6] + dfded[5] * u1x[7]);
    du0xd[7] =
        scale * (dfde[1] * u0xd[7] + dfde[2] * u0xd[6] + dfde[4] * u1xd[7] +
                 dfde[5] * u1xd[6] + dfded[1] * u0x[7] + dfded[2] * u0x[6] +
                 dfded[4] * u1x[7] + dfded[5] * u1x[6]);
    du0xd[8] = 0.0;

    du1x[0] = scale * (dfde[3] * (u0x[0] + 1.0) + dfde[5] * u0x[1]);
    du1x[1] = scale * (dfde[4] * u0x[1] + dfde[5] * (u0x[0] + 1.0));
    du1x[2] = 0.0;
    du1x[3] = scale * (dfde[3] * u0x[3] + dfde[5] * (u0x[4] + 1.0));
    du1x[4] = scale * (dfde[4] * (u0x[4] + 1.0) + dfde[5] * u0x[3]);
    du1x[5] = 0.0;
    du1x[6] = scale * (dfde[3] * u0x[6] + dfde[5] * u0x[7]);
    du1x[7] = scale * (dfde[4] * u0x[7] + dfde[5] * u0x[6]);
    du1x[8] = 0.0;

    du1xd[0] = scale * (dfde[3] * u0xd[0] + dfde[5] * u0xd[1] +
                        dfded[3] * (u0x[0] + 1.0) + dfded[5] * u0x[1]);
    du1xd[1] = scale * (dfde[4] * u0xd[1] + dfde[5] * u0xd[0] +
                        dfded[4] * u0x[1] + dfded[5] * (u0x[0] + 1.0));
    du1xd[2] = 0.0;
    du1xd[3] = scale * (dfde[3] * u0xd[3] + dfde[5] * u0xd[4] +
                        dfded[3] * u0x[3] + dfded[5] * (u0x[4] + 1.0));
    du1xd[4] = scale * (dfde[4] * u0xd[4] + dfde[5] * u0xd[3] +
                        dfded[4] * (u0x[4] + 1.0) + dfded[5] * u0x[3]);
    du1xd[5] = 0.0;
    du1xd[6] = scale * (dfde[3] * u0xd[6] + dfde[5] * u0xd[7] +
                        dfded[3] * u0x[6] + dfded[5] * u0x[7]);
    du1xd[7] = scale * (dfde[4] * u0xd[7] + dfde[5] * u0xd[6] +
                        dfded[4] * u0x[7] + dfded[5] * u0x[6]);
    du1xd[8] = 0.0;
  }

  static void evalStrainDeriv(const TacsScalar u0x[], const TacsScalar u1x[],
                              const TacsScalar e0ty[], const TacsScalar u0xd[],
                              const TacsScalar u1xd[], const TacsScalar e0tyd[],
                              TacsScalar e[], TacsScalar ed[]) {
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0] + 0.5 * (u0x[0] * u0x[0] + u0x[3] * u0x[3] + u0x[6] * u0x[6]);
    e[1] = u0x[4] + 0.5 * (u0x[1] * u0x[1] + u0x[4] * u0x[4] + u0x[7] * u0x[7]);
    e[2] =
        u0x[1] + u0x[3] + (u0x[0] * u0x[1] + u0x[3] * u0x[4] + u0x[6] * u0x[7]);

    // Compute the bending strain
    e[3] = u1x[0] + (u0x[0] * u1x[0] + u0x[3] * u1x[3] + u0x[6] * u1x[6]);
    e[4] = u1x[4] + (u0x[1] * u1x[1] + u0x[4] * u1x[4] + u0x[7] * u1x[7]);
    e[5] = u1x[1] + u1x[3] +
           (u0x[0] * u1x[1] + u0x[3] * u1x[4] + u0x[6] * u1x[7] +
            u1x[0] * u0x[1] + u1x[3] * u0x[4] + u1x[6] * u0x[7]);

    // Add the components of the shear strain
    e[6] = 2.0 * e0ty[4];
    e[7] = 2.0 * e0ty[2];

    // Evaluate the in-plane strains from the tying strain expressions
    ed[0] = u0xd[0] + (u0x[0] * u0xd[0] + u0x[3] * u0xd[3] + u0x[6] * u0xd[6]);
    ed[1] = u0xd[4] + (u0x[1] * u0xd[1] + u0x[4] * u0xd[4] + u0x[7] * u0xd[7]);
    ed[2] = u0xd[1] + u0xd[3] +
            (u0xd[0] * u0x[1] + u0xd[3] * u0x[4] + u0xd[6] * u0x[7] +
             u0x[0] * u0xd[1] + u0x[3] * u0xd[4] + u0x[6] * u0xd[7]);

    // Compute the bending strain
    ed[3] = u1xd[0] + (u0xd[0] * u1x[0] + u0xd[3] * u1x[3] + u0xd[6] * u1x[6] +
                       u0x[0] * u1xd[0] + u0x[3] * u1xd[3] + u0x[6] * u1xd[6]);
    ed[4] = u1xd[4] + (u0xd[1] * u1x[1] + u0xd[4] * u1x[4] + u0xd[7] * u1x[7] +
                       u0x[1] * u1xd[1] + u0x[4] * u1xd[4] + u0x[7] * u1xd[7]);
    ed[5] = u1xd[1] + u1xd[3] +
            (u0xd[0] * u1x[1] + u0xd[3] * u1x[4] + u0xd[6] * u1x[7] +
             u1xd[0] * u0x[1] + u1xd[3] * u0x[4] + u1xd[6] * u0x[7] +
             u0x[0] * u1xd[1] + u0x[3] * u1xd[4] + u0x[6] * u1xd[7] +
             u1x[0] * u0xd[1] + u1x[3] * u0xd[4] + u1x[6] * u0xd[7]);

    // Add the components of the shear strain
    ed[6] = 2.0 * e0tyd[4];
    ed[7] = 2.0 * e0tyd[2];
  }

  static void evalStrainHessian(const TacsScalar scale, const TacsScalar s[],
                                const TacsScalar Cs[], const TacsScalar u0x[],
                                const TacsScalar u1x[], const TacsScalar e0ty[],
                                TacsScalar d2u0x[], TacsScalar d2u1x[],
                                TacsScalar d2u0xu1x[], TacsScalar d2e0ty[],
                                TacsScalar d2e0tyu0x[],
                                TacsScalar d2e0tyu1x[]) {
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    // Compute the second derivatives
    memset(d2e0ty, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54 * sizeof(TacsScalar));

    TacsScalar *d2 = &d2e0ty[0];
    d2 = &d2e0ty[4 * 6];
    d2[4] = 4.0 * scale * As[0];
    d2[2] = 4.0 * scale * As[1];

    d2 = &d2e0ty[2 * 6];
    d2[4] = 4.0 * scale * As[1];
    d2[2] = 4.0 * scale * As[2];

    TacsScalar Du0x[36];
    Du0x[0] =
        A[0] * (u0x[0] + 1.0) + A[2] * u0x[1] + B[0] * u1x[0] + B[2] * u1x[1];
    Du0x[1] =
        A[1] * (u0x[0] + 1.0) + A[4] * u0x[1] + B[1] * u1x[0] + B[4] * u1x[1];
    Du0x[2] =
        A[2] * (u0x[0] + 1.0) + A[5] * u0x[1] + B[2] * u1x[0] + B[5] * u1x[1];
    Du0x[3] =
        B[0] * (u0x[0] + 1.0) + B[2] * u0x[1] + D[0] * u1x[0] + D[2] * u1x[1];
    Du0x[4] =
        B[1] * (u0x[0] + 1.0) + B[4] * u0x[1] + D[1] * u1x[0] + D[4] * u1x[1];
    Du0x[5] =
        B[2] * (u0x[0] + 1.0) + B[5] * u0x[1] + D[2] * u1x[0] + D[5] * u1x[1];
    Du0x[6] =
        A[1] * u0x[1] + A[2] * (u0x[0] + 1.0) + B[1] * u1x[1] + B[2] * u1x[0];
    Du0x[7] =
        A[3] * u0x[1] + A[4] * (u0x[0] + 1.0) + B[3] * u1x[1] + B[4] * u1x[0];
    Du0x[8] =
        A[4] * u0x[1] + A[5] * (u0x[0] + 1.0) + B[4] * u1x[1] + B[5] * u1x[0];
    Du0x[9] =
        B[1] * u0x[1] + B[2] * (u0x[0] + 1.0) + D[1] * u1x[1] + D[2] * u1x[0];
    Du0x[10] =
        B[3] * u0x[1] + B[4] * (u0x[0] + 1.0) + D[3] * u1x[1] + D[4] * u1x[0];
    Du0x[11] =
        B[4] * u0x[1] + B[5] * (u0x[0] + 1.0) + D[4] * u1x[1] + D[5] * u1x[0];
    Du0x[12] =
        A[0] * u0x[3] + A[2] * (u0x[4] + 1.0) + B[0] * u1x[3] + B[2] * u1x[4];
    Du0x[13] =
        A[1] * u0x[3] + A[4] * (u0x[4] + 1.0) + B[1] * u1x[3] + B[4] * u1x[4];
    Du0x[14] =
        A[2] * u0x[3] + A[5] * (u0x[4] + 1.0) + B[2] * u1x[3] + B[5] * u1x[4];
    Du0x[15] =
        B[0] * u0x[3] + B[2] * (u0x[4] + 1.0) + D[0] * u1x[3] + D[2] * u1x[4];
    Du0x[16] =
        B[1] * u0x[3] + B[4] * (u0x[4] + 1.0) + D[1] * u1x[3] + D[4] * u1x[4];
    Du0x[17] =
        B[2] * u0x[3] + B[5] * (u0x[4] + 1.0) + D[2] * u1x[3] + D[5] * u1x[4];
    Du0x[18] =
        A[1] * (u0x[4] + 1.0) + A[2] * u0x[3] + B[1] * u1x[4] + B[2] * u1x[3];
    Du0x[19] =
        A[3] * (u0x[4] + 1.0) + A[4] * u0x[3] + B[3] * u1x[4] + B[4] * u1x[3];
    Du0x[20] =
        A[4] * (u0x[4] + 1.0) + A[5] * u0x[3] + B[4] * u1x[4] + B[5] * u1x[3];
    Du0x[21] =
        B[1] * (u0x[4] + 1.0) + B[2] * u0x[3] + D[1] * u1x[4] + D[2] * u1x[3];
    Du0x[22] =
        B[3] * (u0x[4] + 1.0) + B[4] * u0x[3] + D[3] * u1x[4] + D[4] * u1x[3];
    Du0x[23] =
        B[4] * (u0x[4] + 1.0) + B[5] * u0x[3] + D[4] * u1x[4] + D[5] * u1x[3];
    Du0x[24] = A[0] * u0x[6] + A[2] * u0x[7] + B[0] * u1x[6] + B[2] * u1x[7];
    Du0x[25] = A[1] * u0x[6] + A[4] * u0x[7] + B[1] * u1x[6] + B[4] * u1x[7];
    Du0x[26] = A[2] * u0x[6] + A[5] * u0x[7] + B[2] * u1x[6] + B[5] * u1x[7];
    Du0x[27] = B[0] * u0x[6] + B[2] * u0x[7] + D[0] * u1x[6] + D[2] * u1x[7];
    Du0x[28] = B[1] * u0x[6] + B[4] * u0x[7] + D[1] * u1x[6] + D[4] * u1x[7];
    Du0x[29] = B[2] * u0x[6] + B[5] * u0x[7] + D[2] * u1x[6] + D[5] * u1x[7];
    Du0x[30] = A[1] * u0x[7] + A[2] * u0x[6] + B[1] * u1x[7] + B[2] * u1x[6];
    Du0x[31] = A[3] * u0x[7] + A[4] * u0x[6] + B[3] * u1x[7] + B[4] * u1x[6];
    Du0x[32] = A[4] * u0x[7] + A[5] * u0x[6] + B[4] * u1x[7] + B[5] * u1x[6];
    Du0x[33] = B[1] * u0x[7] + B[2] * u0x[6] + D[1] * u1x[7] + D[2] * u1x[6];
    Du0x[34] = B[3] * u0x[7] + B[4] * u0x[6] + D[3] * u1x[7] + D[4] * u1x[6];
    Du0x[35] = B[4] * u0x[7] + B[5] * u0x[6] + D[4] * u1x[7] + D[5] * u1x[6];

    TacsScalar Du1x[36];
    Du1x[0] = B[0] * (u0x[0] + 1.0) + B[2] * u0x[1];
    Du1x[1] = B[1] * (u0x[0] + 1.0) + B[4] * u0x[1];
    Du1x[2] = B[2] * (u0x[0] + 1.0) + B[5] * u0x[1];
    Du1x[3] = D[0] * (u0x[0] + 1.0) + D[2] * u0x[1];
    Du1x[4] = D[1] * (u0x[0] + 1.0) + D[4] * u0x[1];
    Du1x[5] = D[2] * (u0x[0] + 1.0) + D[5] * u0x[1];
    Du1x[6] = B[1] * u0x[1] + B[2] * (u0x[0] + 1.0);
    Du1x[7] = B[3] * u0x[1] + B[4] * (u0x[0] + 1.0);
    Du1x[8] = B[4] * u0x[1] + B[5] * (u0x[0] + 1.0);
    Du1x[9] = D[1] * u0x[1] + D[2] * (u0x[0] + 1.0);
    Du1x[10] = D[3] * u0x[1] + D[4] * (u0x[0] + 1.0);
    Du1x[11] = D[4] * u0x[1] + D[5] * (u0x[0] + 1.0);
    Du1x[12] = B[0] * u0x[3] + B[2] * (u0x[4] + 1.0);
    Du1x[13] = B[1] * u0x[3] + B[4] * (u0x[4] + 1.0);
    Du1x[14] = B[2] * u0x[3] + B[5] * (u0x[4] + 1.0);
    Du1x[15] = D[0] * u0x[3] + D[2] * (u0x[4] + 1.0);
    Du1x[16] = D[1] * u0x[3] + D[4] * (u0x[4] + 1.0);
    Du1x[17] = D[2] * u0x[3] + D[5] * (u0x[4] + 1.0);
    Du1x[18] = B[1] * (u0x[4] + 1.0) + B[2] * u0x[3];
    Du1x[19] = B[3] * (u0x[4] + 1.0) + B[4] * u0x[3];
    Du1x[20] = B[4] * (u0x[4] + 1.0) + B[5] * u0x[3];
    Du1x[21] = D[1] * (u0x[4] + 1.0) + D[2] * u0x[3];
    Du1x[22] = D[3] * (u0x[4] + 1.0) + D[4] * u0x[3];
    Du1x[23] = D[4] * (u0x[4] + 1.0) + D[5] * u0x[3];
    Du1x[24] = B[0] * u0x[6] + B[2] * u0x[7];
    Du1x[25] = B[1] * u0x[6] + B[4] * u0x[7];
    Du1x[26] = B[2] * u0x[6] + B[5] * u0x[7];
    Du1x[27] = D[0] * u0x[6] + D[2] * u0x[7];
    Du1x[28] = D[1] * u0x[6] + D[4] * u0x[7];
    Du1x[29] = D[2] * u0x[6] + D[5] * u0x[7];
    Du1x[30] = B[1] * u0x[7] + B[2] * u0x[6];
    Du1x[31] = B[3] * u0x[7] + B[4] * u0x[6];
    Du1x[32] = B[4] * u0x[7] + B[5] * u0x[6];
    Du1x[33] = D[1] * u0x[7] + D[2] * u0x[6];
    Du1x[34] = D[3] * u0x[7] + D[4] * u0x[6];
    Du1x[35] = D[4] * u0x[7] + D[5] * u0x[6];

    d2u0x[0] = scale * (Du0x[0] * (u0x[0] + 1.0) + Du0x[2] * u0x[1] +
                        Du0x[3] * u1x[0] + Du0x[5] * u1x[1] + s[0]);
    d2u0x[1] = scale * (Du0x[1] * u0x[1] + Du0x[2] * (u0x[0] + 1.0) +
                        Du0x[4] * u1x[1] + Du0x[5] * u1x[0] + s[2]);
    d2u0x[2] = 0.0;
    d2u0x[3] = scale * (Du0x[0] * u0x[3] + Du0x[2] * (u0x[4] + 1.0) +
                        Du0x[3] * u1x[3] + Du0x[5] * u1x[4]);
    d2u0x[4] = scale * (Du0x[1] * (u0x[4] + 1.0) + Du0x[2] * u0x[3] +
                        Du0x[4] * u1x[4] + Du0x[5] * u1x[3]);
    d2u0x[5] = 0.0;
    d2u0x[6] = scale * (Du0x[0] * u0x[6] + Du0x[2] * u0x[7] + Du0x[3] * u1x[6] +
                        Du0x[5] * u1x[7]);
    d2u0x[7] = scale * (Du0x[1] * u0x[7] + Du0x[2] * u0x[6] + Du0x[4] * u1x[7] +
                        Du0x[5] * u1x[6]);
    d2u0x[8] = 0.0;

    d2u0x[9] = scale * (Du0x[11] * u1x[1] + Du0x[6] * (u0x[0] + 1.0) +
                        Du0x[8] * u0x[1] + Du0x[9] * u1x[0] + s[2]);
    d2u0x[10] = scale * (Du0x[10] * u1x[1] + Du0x[11] * u1x[0] +
                         Du0x[7] * u0x[1] + Du0x[8] * (u0x[0] + 1.0) + s[1]);
    d2u0x[11] = 0.0;
    d2u0x[12] = scale * (Du0x[11] * u1x[4] + Du0x[6] * u0x[3] +
                         Du0x[8] * (u0x[4] + 1.0) + Du0x[9] * u1x[3]);
    d2u0x[13] = scale * (Du0x[10] * u1x[4] + Du0x[11] * u1x[3] +
                         Du0x[7] * (u0x[4] + 1.0) + Du0x[8] * u0x[3]);
    d2u0x[14] = 0.0;
    d2u0x[15] = scale * (Du0x[11] * u1x[7] + Du0x[6] * u0x[6] +
                         Du0x[8] * u0x[7] + Du0x[9] * u1x[6]);
    d2u0x[16] = scale * (Du0x[10] * u1x[7] + Du0x[11] * u1x[6] +
                         Du0x[7] * u0x[7] + Du0x[8] * u0x[6]);
    d2u0x[17] = 0.0;

    d2u0x[18] = 0.0;
    d2u0x[19] = 0.0;
    d2u0x[20] = 0.0;
    d2u0x[21] = 0.0;
    d2u0x[22] = 0.0;
    d2u0x[23] = 0.0;
    d2u0x[24] = 0.0;
    d2u0x[25] = 0.0;
    d2u0x[26] = 0.0;

    d2u0x[27] = scale * (Du0x[12] * (u0x[0] + 1.0) + Du0x[14] * u0x[1] +
                         Du0x[15] * u1x[0] + Du0x[17] * u1x[1]);
    d2u0x[28] = scale * (Du0x[13] * u0x[1] + Du0x[14] * (u0x[0] + 1.0) +
                         Du0x[16] * u1x[1] + Du0x[17] * u1x[0]);
    d2u0x[29] = 0.0;
    d2u0x[30] = scale * (Du0x[12] * u0x[3] + Du0x[14] * (u0x[4] + 1.0) +
                         Du0x[15] * u1x[3] + Du0x[17] * u1x[4] + s[0]);
    d2u0x[31] = scale * (Du0x[13] * (u0x[4] + 1.0) + Du0x[14] * u0x[3] +
                         Du0x[16] * u1x[4] + Du0x[17] * u1x[3] + s[2]);
    d2u0x[32] = 0.0;
    d2u0x[33] = scale * (Du0x[12] * u0x[6] + Du0x[14] * u0x[7] +
                         Du0x[15] * u1x[6] + Du0x[17] * u1x[7]);
    d2u0x[34] = scale * (Du0x[13] * u0x[7] + Du0x[14] * u0x[6] +
                         Du0x[16] * u1x[7] + Du0x[17] * u1x[6]);
    d2u0x[35] = 0.0;

    d2u0x[36] = scale * (Du0x[18] * (u0x[0] + 1.0) + Du0x[20] * u0x[1] +
                         Du0x[21] * u1x[0] + Du0x[23] * u1x[1]);
    d2u0x[37] = scale * (Du0x[19] * u0x[1] + Du0x[20] * (u0x[0] + 1.0) +
                         Du0x[22] * u1x[1] + Du0x[23] * u1x[0]);
    d2u0x[38] = 0.0;
    d2u0x[39] = scale * (Du0x[18] * u0x[3] + Du0x[20] * (u0x[4] + 1.0) +
                         Du0x[21] * u1x[3] + Du0x[23] * u1x[4] + s[2]);
    d2u0x[40] = scale * (Du0x[19] * (u0x[4] + 1.0) + Du0x[20] * u0x[3] +
                         Du0x[22] * u1x[4] + Du0x[23] * u1x[3] + s[1]);
    d2u0x[41] = 0.0;
    d2u0x[42] = scale * (Du0x[18] * u0x[6] + Du0x[20] * u0x[7] +
                         Du0x[21] * u1x[6] + Du0x[23] * u1x[7]);
    d2u0x[43] = scale * (Du0x[19] * u0x[7] + Du0x[20] * u0x[6] +
                         Du0x[22] * u1x[7] + Du0x[23] * u1x[6]);
    d2u0x[44] = 0.0;

    d2u0x[45] = 0.0;
    d2u0x[46] = 0.0;
    d2u0x[47] = 0.0;
    d2u0x[48] = 0.0;
    d2u0x[49] = 0.0;
    d2u0x[50] = 0.0;
    d2u0x[51] = 0.0;
    d2u0x[52] = 0.0;
    d2u0x[53] = 0.0;

    d2u0x[54] = scale * (Du0x[24] * (u0x[0] + 1.0) + Du0x[26] * u0x[1] +
                         Du0x[27] * u1x[0] + Du0x[29] * u1x[1]);
    d2u0x[55] = scale * (Du0x[25] * u0x[1] + Du0x[26] * (u0x[0] + 1.0) +
                         Du0x[28] * u1x[1] + Du0x[29] * u1x[0]);
    d2u0x[56] = 0.0;
    d2u0x[57] = scale * (Du0x[24] * u0x[3] + Du0x[26] * (u0x[4] + 1.0) +
                         Du0x[27] * u1x[3] + Du0x[29] * u1x[4]);
    d2u0x[58] = scale * (Du0x[25] * (u0x[4] + 1.0) + Du0x[26] * u0x[3] +
                         Du0x[28] * u1x[4] + Du0x[29] * u1x[3]);
    d2u0x[59] = 0.0;
    d2u0x[60] = scale * (Du0x[24] * u0x[6] + Du0x[26] * u0x[7] +
                         Du0x[27] * u1x[6] + Du0x[29] * u1x[7] + s[0]);
    d2u0x[61] = scale * (Du0x[25] * u0x[7] + Du0x[26] * u0x[6] +
                         Du0x[28] * u1x[7] + Du0x[29] * u1x[6] + s[2]);
    d2u0x[62] = 0.0;

    d2u0x[63] = scale * (Du0x[30] * (u0x[0] + 1.0) + Du0x[32] * u0x[1] +
                         Du0x[33] * u1x[0] + Du0x[35] * u1x[1]);
    d2u0x[64] = scale * (Du0x[31] * u0x[1] + Du0x[32] * (u0x[0] + 1.0) +
                         Du0x[34] * u1x[1] + Du0x[35] * u1x[0]);
    d2u0x[65] = 0.0;
    d2u0x[66] = scale * (Du0x[30] * u0x[3] + Du0x[32] * (u0x[4] + 1.0) +
                         Du0x[33] * u1x[3] + Du0x[35] * u1x[4]);
    d2u0x[67] = scale * (Du0x[31] * (u0x[4] + 1.0) + Du0x[32] * u0x[3] +
                         Du0x[34] * u1x[4] + Du0x[35] * u1x[3]);
    d2u0x[68] = 0.0;
    d2u0x[69] = scale * (Du0x[30] * u0x[6] + Du0x[32] * u0x[7] +
                         Du0x[33] * u1x[6] + Du0x[35] * u1x[7] + s[2]);
    d2u0x[70] = scale * (Du0x[31] * u0x[7] + Du0x[32] * u0x[6] +
                         Du0x[34] * u1x[7] + Du0x[35] * u1x[6] + s[1]);
    d2u0x[71] = 0.0;

    d2u0x[72] = 0.0;
    d2u0x[73] = 0.0;
    d2u0x[74] = 0.0;
    d2u0x[75] = 0.0;
    d2u0x[76] = 0.0;
    d2u0x[77] = 0.0;
    d2u0x[78] = 0.0;
    d2u0x[79] = 0.0;
    d2u0x[80] = 0.0;

    d2u0xu1x[0] = scale * (Du0x[3] * (u0x[0] + 1.0) + Du0x[5] * u0x[1] + s[3]);
    d2u0xu1x[1] = scale * (Du0x[4] * u0x[1] + Du0x[5] * (u0x[0] + 1.0) + s[5]);
    d2u0xu1x[2] = 0.0;
    d2u0xu1x[3] = scale * (Du0x[3] * u0x[3] + Du0x[5] * (u0x[4] + 1.0));
    d2u0xu1x[4] = scale * (Du0x[4] * (u0x[4] + 1.0) + Du0x[5] * u0x[3]);
    d2u0xu1x[5] = 0.0;
    d2u0xu1x[6] = scale * (Du0x[3] * u0x[6] + Du0x[5] * u0x[7]);
    d2u0xu1x[7] = scale * (Du0x[4] * u0x[7] + Du0x[5] * u0x[6]);
    d2u0xu1x[8] = 0.0;

    d2u0xu1x[9] = scale * (Du0x[11] * u0x[1] + Du0x[9] * (u0x[0] + 1.0) + s[5]);
    d2u0xu1x[10] =
        scale * (Du0x[10] * u0x[1] + Du0x[11] * (u0x[0] + 1.0) + s[4]);
    d2u0xu1x[11] = 0.0;
    d2u0xu1x[12] = scale * (Du0x[11] * (u0x[4] + 1.0) + Du0x[9] * u0x[3]);
    d2u0xu1x[13] = scale * (Du0x[10] * (u0x[4] + 1.0) + Du0x[11] * u0x[3]);
    d2u0xu1x[14] = 0.0;
    d2u0xu1x[15] = scale * (Du0x[11] * u0x[7] + Du0x[9] * u0x[6]);
    d2u0xu1x[16] = scale * (Du0x[10] * u0x[7] + Du0x[11] * u0x[6]);
    d2u0xu1x[17] = 0.0;

    d2u0xu1x[18] = 0.0;
    d2u0xu1x[19] = 0.0;
    d2u0xu1x[20] = 0.0;
    d2u0xu1x[21] = 0.0;
    d2u0xu1x[22] = 0.0;
    d2u0xu1x[23] = 0.0;
    d2u0xu1x[24] = 0.0;
    d2u0xu1x[25] = 0.0;
    d2u0xu1x[26] = 0.0;

    d2u0xu1x[27] = scale * (Du0x[15] * (u0x[0] + 1.0) + Du0x[17] * u0x[1]);
    d2u0xu1x[28] = scale * (Du0x[16] * u0x[1] + Du0x[17] * (u0x[0] + 1.0));
    d2u0xu1x[29] = 0.0;
    d2u0xu1x[30] =
        scale * (Du0x[15] * u0x[3] + Du0x[17] * (u0x[4] + 1.0) + s[3]);
    d2u0xu1x[31] =
        scale * (Du0x[16] * (u0x[4] + 1.0) + Du0x[17] * u0x[3] + s[5]);
    d2u0xu1x[32] = 0.0;
    d2u0xu1x[33] = scale * (Du0x[15] * u0x[6] + Du0x[17] * u0x[7]);
    d2u0xu1x[34] = scale * (Du0x[16] * u0x[7] + Du0x[17] * u0x[6]);
    d2u0xu1x[35] = 0.0;

    d2u0xu1x[36] = scale * (Du0x[21] * (u0x[0] + 1.0) + Du0x[23] * u0x[1]);
    d2u0xu1x[37] = scale * (Du0x[22] * u0x[1] + Du0x[23] * (u0x[0] + 1.0));
    d2u0xu1x[38] = 0.0;
    d2u0xu1x[39] =
        scale * (Du0x[21] * u0x[3] + Du0x[23] * (u0x[4] + 1.0) + s[5]);
    d2u0xu1x[40] =
        scale * (Du0x[22] * (u0x[4] + 1.0) + Du0x[23] * u0x[3] + s[4]);
    d2u0xu1x[41] = 0.0;
    d2u0xu1x[42] = scale * (Du0x[21] * u0x[6] + Du0x[23] * u0x[7]);
    d2u0xu1x[43] = scale * (Du0x[22] * u0x[7] + Du0x[23] * u0x[6]);
    d2u0xu1x[44] = 0.0;

    d2u0xu1x[45] = 0.0;
    d2u0xu1x[46] = 0.0;
    d2u0xu1x[47] = 0.0;
    d2u0xu1x[48] = 0.0;
    d2u0xu1x[49] = 0.0;
    d2u0xu1x[50] = 0.0;
    d2u0xu1x[51] = 0.0;
    d2u0xu1x[52] = 0.0;
    d2u0xu1x[53] = 0.0;

    d2u0xu1x[54] = scale * (Du0x[27] * (u0x[0] + 1.0) + Du0x[29] * u0x[1]);
    d2u0xu1x[55] = scale * (Du0x[28] * u0x[1] + Du0x[29] * (u0x[0] + 1.0));
    d2u0xu1x[56] = 0.0;
    d2u0xu1x[57] = scale * (Du0x[27] * u0x[3] + Du0x[29] * (u0x[4] + 1.0));
    d2u0xu1x[58] = scale * (Du0x[28] * (u0x[4] + 1.0) + Du0x[29] * u0x[3]);
    d2u0xu1x[59] = 0.0;
    d2u0xu1x[60] = scale * (Du0x[27] * u0x[6] + Du0x[29] * u0x[7] + s[3]);
    d2u0xu1x[61] = scale * (Du0x[28] * u0x[7] + Du0x[29] * u0x[6] + s[5]);
    d2u0xu1x[62] = 0.0;

    d2u0xu1x[63] = scale * (Du0x[33] * (u0x[0] + 1.0) + Du0x[35] * u0x[1]);
    d2u0xu1x[64] = scale * (Du0x[34] * u0x[1] + Du0x[35] * (u0x[0] + 1.0));
    d2u0xu1x[65] = 0.0;
    d2u0xu1x[66] = scale * (Du0x[33] * u0x[3] + Du0x[35] * (u0x[4] + 1.0));
    d2u0xu1x[67] = scale * (Du0x[34] * (u0x[4] + 1.0) + Du0x[35] * u0x[3]);
    d2u0xu1x[68] = 0.0;
    d2u0xu1x[69] = scale * (Du0x[33] * u0x[6] + Du0x[35] * u0x[7] + s[5]);
    d2u0xu1x[70] = scale * (Du0x[34] * u0x[7] + Du0x[35] * u0x[6] + s[4]);
    d2u0xu1x[71] = 0.0;

    d2u0xu1x[72] = 0.0;
    d2u0xu1x[73] = 0.0;
    d2u0xu1x[74] = 0.0;
    d2u0xu1x[75] = 0.0;
    d2u0xu1x[76] = 0.0;
    d2u0xu1x[77] = 0.0;
    d2u0xu1x[78] = 0.0;
    d2u0xu1x[79] = 0.0;
    d2u0xu1x[80] = 0.0;

    d2u1x[0] = scale * (Du1x[3] * (u0x[0] + 1.0) + Du1x[5] * u0x[1]);
    d2u1x[1] = scale * (Du1x[4] * u0x[1] + Du1x[5] * (u0x[0] + 1.0));
    d2u1x[2] = 0.0;
    d2u1x[3] = scale * (Du1x[3] * u0x[3] + Du1x[5] * (u0x[4] + 1.0));
    d2u1x[4] = scale * (Du1x[4] * (u0x[4] + 1.0) + Du1x[5] * u0x[3]);
    d2u1x[5] = 0.0;
    d2u1x[6] = scale * (Du1x[3] * u0x[6] + Du1x[5] * u0x[7]);
    d2u1x[7] = scale * (Du1x[4] * u0x[7] + Du1x[5] * u0x[6]);
    d2u1x[8] = 0.0;

    d2u1x[9] = scale * (Du1x[11] * u0x[1] + Du1x[9] * (u0x[0] + 1.0));
    d2u1x[10] = scale * (Du1x[10] * u0x[1] + Du1x[11] * (u0x[0] + 1.0));
    d2u1x[11] = 0.0;
    d2u1x[12] = scale * (Du1x[11] * (u0x[4] + 1.0) + Du1x[9] * u0x[3]);
    d2u1x[13] = scale * (Du1x[10] * (u0x[4] + 1.0) + Du1x[11] * u0x[3]);
    d2u1x[14] = 0.0;
    d2u1x[15] = scale * (Du1x[11] * u0x[7] + Du1x[9] * u0x[6]);
    d2u1x[16] = scale * (Du1x[10] * u0x[7] + Du1x[11] * u0x[6]);
    d2u1x[17] = 0.0;

    d2u1x[18] = 0.0;
    d2u1x[19] = 0.0;
    d2u1x[20] = 0.0;
    d2u1x[21] = 0.0;
    d2u1x[22] = 0.0;
    d2u1x[23] = 0.0;
    d2u1x[24] = 0.0;
    d2u1x[25] = 0.0;
    d2u1x[26] = 0.0;

    d2u1x[27] = scale * (Du1x[15] * (u0x[0] + 1.0) + Du1x[17] * u0x[1]);
    d2u1x[28] = scale * (Du1x[16] * u0x[1] + Du1x[17] * (u0x[0] + 1.0));
    d2u1x[29] = 0.0;
    d2u1x[30] = scale * (Du1x[15] * u0x[3] + Du1x[17] * (u0x[4] + 1.0));
    d2u1x[31] = scale * (Du1x[16] * (u0x[4] + 1.0) + Du1x[17] * u0x[3]);
    d2u1x[32] = 0.0;
    d2u1x[33] = scale * (Du1x[15] * u0x[6] + Du1x[17] * u0x[7]);
    d2u1x[34] = scale * (Du1x[16] * u0x[7] + Du1x[17] * u0x[6]);
    d2u1x[35] = 0.0;

    d2u1x[36] = scale * (Du1x[21] * (u0x[0] + 1.0) + Du1x[23] * u0x[1]);
    d2u1x[37] = scale * (Du1x[22] * u0x[1] + Du1x[23] * (u0x[0] + 1.0));
    d2u1x[38] = 0.0;
    d2u1x[39] = scale * (Du1x[21] * u0x[3] + Du1x[23] * (u0x[4] + 1.0));
    d2u1x[40] = scale * (Du1x[22] * (u0x[4] + 1.0) + Du1x[23] * u0x[3]);
    d2u1x[41] = 0.0;
    d2u1x[42] = scale * (Du1x[21] * u0x[6] + Du1x[23] * u0x[7]);
    d2u1x[43] = scale * (Du1x[22] * u0x[7] + Du1x[23] * u0x[6]);
    d2u1x[44] = 0.0;

    d2u1x[45] = 0.0;
    d2u1x[46] = 0.0;
    d2u1x[47] = 0.0;
    d2u1x[48] = 0.0;
    d2u1x[49] = 0.0;
    d2u1x[50] = 0.0;
    d2u1x[51] = 0.0;
    d2u1x[52] = 0.0;
    d2u1x[53] = 0.0;

    d2u1x[54] = scale * (Du1x[27] * (u0x[0] + 1.0) + Du1x[29] * u0x[1]);
    d2u1x[55] = scale * (Du1x[28] * u0x[1] + Du1x[29] * (u0x[0] + 1.0));
    d2u1x[56] = 0.0;
    d2u1x[57] = scale * (Du1x[27] * u0x[3] + Du1x[29] * (u0x[4] + 1.0));
    d2u1x[58] = scale * (Du1x[28] * (u0x[4] + 1.0) + Du1x[29] * u0x[3]);
    d2u1x[59] = 0.0;
    d2u1x[60] = scale * (Du1x[27] * u0x[6] + Du1x[29] * u0x[7]);
    d2u1x[61] = scale * (Du1x[28] * u0x[7] + Du1x[29] * u0x[6]);
    d2u1x[62] = 0.0;

    d2u1x[63] = scale * (Du1x[33] * (u0x[0] + 1.0) + Du1x[35] * u0x[1]);
    d2u1x[64] = scale * (Du1x[34] * u0x[1] + Du1x[35] * (u0x[0] + 1.0));
    d2u1x[65] = 0.0;
    d2u1x[66] = scale * (Du1x[33] * u0x[3] + Du1x[35] * (u0x[4] + 1.0));
    d2u1x[67] = scale * (Du1x[34] * (u0x[4] + 1.0) + Du1x[35] * u0x[3]);
    d2u1x[68] = 0.0;
    d2u1x[69] = scale * (Du1x[33] * u0x[6] + Du1x[35] * u0x[7]);
    d2u1x[70] = scale * (Du1x[34] * u0x[7] + Du1x[35] * u0x[6]);
    d2u1x[71] = 0.0;

    d2u1x[72] = 0.0;
    d2u1x[73] = 0.0;
    d2u1x[74] = 0.0;
    d2u1x[75] = 0.0;
    d2u1x[76] = 0.0;
    d2u1x[77] = 0.0;
    d2u1x[78] = 0.0;
    d2u1x[79] = 0.0;
    d2u1x[80] = 0.0;
  }

  static inline void evalStrainHessianDeriv(
      const TacsScalar scale, const TacsScalar s[], const TacsScalar Cs[],
      const TacsScalar u0x[], const TacsScalar u1x[], const TacsScalar e0ty[],
      const TacsScalar sd[], const TacsScalar u0xd[], const TacsScalar u1xd[],
      const TacsScalar e0tyd[], TacsScalar d2u0x[], TacsScalar d2u1x[],
      TacsScalar d2u0xu1x[], TacsScalar d2e0ty[], TacsScalar d2e0tyu0x[],
      TacsScalar d2e0tyu1x[], TacsScalar d2u0xd[], TacsScalar d2u1xd[],
      TacsScalar d2u0xu1xd[], TacsScalar d2e0tyd[], TacsScalar d2e0tyu0xd[],
      TacsScalar d2e0tyu1xd[]) {
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    // Compute the second derivatives
    memset(d2e0ty, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyd, 0, 36 * sizeof(TacsScalar));
    memset(d2e0tyu0xd, 0, 54 * sizeof(TacsScalar));
    memset(d2e0tyu1xd, 0, 54 * sizeof(TacsScalar));

    TacsScalar *d2 = &d2e0ty[0];
    d2 = &d2e0ty[4 * 6];
    d2[4] = 4.0 * scale * As[0];
    d2[2] = 4.0 * scale * As[1];

    d2 = &d2e0ty[2 * 6];
    d2[4] = 4.0 * scale * As[1];
    d2[2] = 4.0 * scale * As[2];

    TacsScalar Du0x[36];
    Du0x[0] =
        A[0] * (u0x[0] + 1.0) + A[2] * u0x[1] + B[0] * u1x[0] + B[2] * u1x[1];
    Du0x[1] =
        A[1] * (u0x[0] + 1.0) + A[4] * u0x[1] + B[1] * u1x[0] + B[4] * u1x[1];
    Du0x[2] =
        A[2] * (u0x[0] + 1.0) + A[5] * u0x[1] + B[2] * u1x[0] + B[5] * u1x[1];
    Du0x[3] =
        B[0] * (u0x[0] + 1.0) + B[2] * u0x[1] + D[0] * u1x[0] + D[2] * u1x[1];
    Du0x[4] =
        B[1] * (u0x[0] + 1.0) + B[4] * u0x[1] + D[1] * u1x[0] + D[4] * u1x[1];
    Du0x[5] =
        B[2] * (u0x[0] + 1.0) + B[5] * u0x[1] + D[2] * u1x[0] + D[5] * u1x[1];
    Du0x[6] =
        A[1] * u0x[1] + A[2] * (u0x[0] + 1.0) + B[1] * u1x[1] + B[2] * u1x[0];
    Du0x[7] =
        A[3] * u0x[1] + A[4] * (u0x[0] + 1.0) + B[3] * u1x[1] + B[4] * u1x[0];
    Du0x[8] =
        A[4] * u0x[1] + A[5] * (u0x[0] + 1.0) + B[4] * u1x[1] + B[5] * u1x[0];
    Du0x[9] =
        B[1] * u0x[1] + B[2] * (u0x[0] + 1.0) + D[1] * u1x[1] + D[2] * u1x[0];
    Du0x[10] =
        B[3] * u0x[1] + B[4] * (u0x[0] + 1.0) + D[3] * u1x[1] + D[4] * u1x[0];
    Du0x[11] =
        B[4] * u0x[1] + B[5] * (u0x[0] + 1.0) + D[4] * u1x[1] + D[5] * u1x[0];
    Du0x[12] =
        A[0] * u0x[3] + A[2] * (u0x[4] + 1.0) + B[0] * u1x[3] + B[2] * u1x[4];
    Du0x[13] =
        A[1] * u0x[3] + A[4] * (u0x[4] + 1.0) + B[1] * u1x[3] + B[4] * u1x[4];
    Du0x[14] =
        A[2] * u0x[3] + A[5] * (u0x[4] + 1.0) + B[2] * u1x[3] + B[5] * u1x[4];
    Du0x[15] =
        B[0] * u0x[3] + B[2] * (u0x[4] + 1.0) + D[0] * u1x[3] + D[2] * u1x[4];
    Du0x[16] =
        B[1] * u0x[3] + B[4] * (u0x[4] + 1.0) + D[1] * u1x[3] + D[4] * u1x[4];
    Du0x[17] =
        B[2] * u0x[3] + B[5] * (u0x[4] + 1.0) + D[2] * u1x[3] + D[5] * u1x[4];
    Du0x[18] =
        A[1] * (u0x[4] + 1.0) + A[2] * u0x[3] + B[1] * u1x[4] + B[2] * u1x[3];
    Du0x[19] =
        A[3] * (u0x[4] + 1.0) + A[4] * u0x[3] + B[3] * u1x[4] + B[4] * u1x[3];
    Du0x[20] =
        A[4] * (u0x[4] + 1.0) + A[5] * u0x[3] + B[4] * u1x[4] + B[5] * u1x[3];
    Du0x[21] =
        B[1] * (u0x[4] + 1.0) + B[2] * u0x[3] + D[1] * u1x[4] + D[2] * u1x[3];
    Du0x[22] =
        B[3] * (u0x[4] + 1.0) + B[4] * u0x[3] + D[3] * u1x[4] + D[4] * u1x[3];
    Du0x[23] =
        B[4] * (u0x[4] + 1.0) + B[5] * u0x[3] + D[4] * u1x[4] + D[5] * u1x[3];
    Du0x[24] = A[0] * u0x[6] + A[2] * u0x[7] + B[0] * u1x[6] + B[2] * u1x[7];
    Du0x[25] = A[1] * u0x[6] + A[4] * u0x[7] + B[1] * u1x[6] + B[4] * u1x[7];
    Du0x[26] = A[2] * u0x[6] + A[5] * u0x[7] + B[2] * u1x[6] + B[5] * u1x[7];
    Du0x[27] = B[0] * u0x[6] + B[2] * u0x[7] + D[0] * u1x[6] + D[2] * u1x[7];
    Du0x[28] = B[1] * u0x[6] + B[4] * u0x[7] + D[1] * u1x[6] + D[4] * u1x[7];
    Du0x[29] = B[2] * u0x[6] + B[5] * u0x[7] + D[2] * u1x[6] + D[5] * u1x[7];
    Du0x[30] = A[1] * u0x[7] + A[2] * u0x[6] + B[1] * u1x[7] + B[2] * u1x[6];
    Du0x[31] = A[3] * u0x[7] + A[4] * u0x[6] + B[3] * u1x[7] + B[4] * u1x[6];
    Du0x[32] = A[4] * u0x[7] + A[5] * u0x[6] + B[4] * u1x[7] + B[5] * u1x[6];
    Du0x[33] = B[1] * u0x[7] + B[2] * u0x[6] + D[1] * u1x[7] + D[2] * u1x[6];
    Du0x[34] = B[3] * u0x[7] + B[4] * u0x[6] + D[3] * u1x[7] + D[4] * u1x[6];
    Du0x[35] = B[4] * u0x[7] + B[5] * u0x[6] + D[4] * u1x[7] + D[5] * u1x[6];

    TacsScalar Du0xd[36];
    Du0xd[0] =
        A[0] * u0xd[0] + A[2] * u0xd[1] + B[0] * u1xd[0] + B[2] * u1xd[1];
    Du0xd[1] =
        A[1] * u0xd[0] + A[4] * u0xd[1] + B[1] * u1xd[0] + B[4] * u1xd[1];
    Du0xd[2] =
        A[2] * u0xd[0] + A[5] * u0xd[1] + B[2] * u1xd[0] + B[5] * u1xd[1];
    Du0xd[3] =
        B[0] * u0xd[0] + B[2] * u0xd[1] + D[0] * u1xd[0] + D[2] * u1xd[1];
    Du0xd[4] =
        B[1] * u0xd[0] + B[4] * u0xd[1] + D[1] * u1xd[0] + D[4] * u1xd[1];
    Du0xd[5] =
        B[2] * u0xd[0] + B[5] * u0xd[1] + D[2] * u1xd[0] + D[5] * u1xd[1];
    Du0xd[6] =
        A[1] * u0xd[1] + A[2] * u0xd[0] + B[1] * u1xd[1] + B[2] * u1xd[0];
    Du0xd[7] =
        A[3] * u0xd[1] + A[4] * u0xd[0] + B[3] * u1xd[1] + B[4] * u1xd[0];
    Du0xd[8] =
        A[4] * u0xd[1] + A[5] * u0xd[0] + B[4] * u1xd[1] + B[5] * u1xd[0];
    Du0xd[9] =
        B[1] * u0xd[1] + B[2] * u0xd[0] + D[1] * u1xd[1] + D[2] * u1xd[0];
    Du0xd[10] =
        B[3] * u0xd[1] + B[4] * u0xd[0] + D[3] * u1xd[1] + D[4] * u1xd[0];
    Du0xd[11] =
        B[4] * u0xd[1] + B[5] * u0xd[0] + D[4] * u1xd[1] + D[5] * u1xd[0];
    Du0xd[12] =
        A[0] * u0xd[3] + A[2] * u0xd[4] + B[0] * u1xd[3] + B[2] * u1xd[4];
    Du0xd[13] =
        A[1] * u0xd[3] + A[4] * u0xd[4] + B[1] * u1xd[3] + B[4] * u1xd[4];
    Du0xd[14] =
        A[2] * u0xd[3] + A[5] * u0xd[4] + B[2] * u1xd[3] + B[5] * u1xd[4];
    Du0xd[15] =
        B[0] * u0xd[3] + B[2] * u0xd[4] + D[0] * u1xd[3] + D[2] * u1xd[4];
    Du0xd[16] =
        B[1] * u0xd[3] + B[4] * u0xd[4] + D[1] * u1xd[3] + D[4] * u1xd[4];
    Du0xd[17] =
        B[2] * u0xd[3] + B[5] * u0xd[4] + D[2] * u1xd[3] + D[5] * u1xd[4];
    Du0xd[18] =
        A[1] * u0xd[4] + A[2] * u0xd[3] + B[1] * u1xd[4] + B[2] * u1xd[3];
    Du0xd[19] =
        A[3] * u0xd[4] + A[4] * u0xd[3] + B[3] * u1xd[4] + B[4] * u1xd[3];
    Du0xd[20] =
        A[4] * u0xd[4] + A[5] * u0xd[3] + B[4] * u1xd[4] + B[5] * u1xd[3];
    Du0xd[21] =
        B[1] * u0xd[4] + B[2] * u0xd[3] + D[1] * u1xd[4] + D[2] * u1xd[3];
    Du0xd[22] =
        B[3] * u0xd[4] + B[4] * u0xd[3] + D[3] * u1xd[4] + D[4] * u1xd[3];
    Du0xd[23] =
        B[4] * u0xd[4] + B[5] * u0xd[3] + D[4] * u1xd[4] + D[5] * u1xd[3];
    Du0xd[24] =
        A[0] * u0xd[6] + A[2] * u0xd[7] + B[0] * u1xd[6] + B[2] * u1xd[7];
    Du0xd[25] =
        A[1] * u0xd[6] + A[4] * u0xd[7] + B[1] * u1xd[6] + B[4] * u1xd[7];
    Du0xd[26] =
        A[2] * u0xd[6] + A[5] * u0xd[7] + B[2] * u1xd[6] + B[5] * u1xd[7];
    Du0xd[27] =
        B[0] * u0xd[6] + B[2] * u0xd[7] + D[0] * u1xd[6] + D[2] * u1xd[7];
    Du0xd[28] =
        B[1] * u0xd[6] + B[4] * u0xd[7] + D[1] * u1xd[6] + D[4] * u1xd[7];
    Du0xd[29] =
        B[2] * u0xd[6] + B[5] * u0xd[7] + D[2] * u1xd[6] + D[5] * u1xd[7];
    Du0xd[30] =
        A[1] * u0xd[7] + A[2] * u0xd[6] + B[1] * u1xd[7] + B[2] * u1xd[6];
    Du0xd[31] =
        A[3] * u0xd[7] + A[4] * u0xd[6] + B[3] * u1xd[7] + B[4] * u1xd[6];
    Du0xd[32] =
        A[4] * u0xd[7] + A[5] * u0xd[6] + B[4] * u1xd[7] + B[5] * u1xd[6];
    Du0xd[33] =
        B[1] * u0xd[7] + B[2] * u0xd[6] + D[1] * u1xd[7] + D[2] * u1xd[6];
    Du0xd[34] =
        B[3] * u0xd[7] + B[4] * u0xd[6] + D[3] * u1xd[7] + D[4] * u1xd[6];
    Du0xd[35] =
        B[4] * u0xd[7] + B[5] * u0xd[6] + D[4] * u1xd[7] + D[5] * u1xd[6];

    TacsScalar Du1x[36];
    Du1x[0] = B[0] * (u0x[0] + 1.0) + B[2] * u0x[1];
    Du1x[1] = B[1] * (u0x[0] + 1.0) + B[4] * u0x[1];
    Du1x[2] = B[2] * (u0x[0] + 1.0) + B[5] * u0x[1];
    Du1x[3] = D[0] * (u0x[0] + 1.0) + D[2] * u0x[1];
    Du1x[4] = D[1] * (u0x[0] + 1.0) + D[4] * u0x[1];
    Du1x[5] = D[2] * (u0x[0] + 1.0) + D[5] * u0x[1];
    Du1x[6] = B[1] * u0x[1] + B[2] * (u0x[0] + 1.0);
    Du1x[7] = B[3] * u0x[1] + B[4] * (u0x[0] + 1.0);
    Du1x[8] = B[4] * u0x[1] + B[5] * (u0x[0] + 1.0);
    Du1x[9] = D[1] * u0x[1] + D[2] * (u0x[0] + 1.0);
    Du1x[10] = D[3] * u0x[1] + D[4] * (u0x[0] + 1.0);
    Du1x[11] = D[4] * u0x[1] + D[5] * (u0x[0] + 1.0);
    Du1x[12] = B[0] * u0x[3] + B[2] * (u0x[4] + 1.0);
    Du1x[13] = B[1] * u0x[3] + B[4] * (u0x[4] + 1.0);
    Du1x[14] = B[2] * u0x[3] + B[5] * (u0x[4] + 1.0);
    Du1x[15] = D[0] * u0x[3] + D[2] * (u0x[4] + 1.0);
    Du1x[16] = D[1] * u0x[3] + D[4] * (u0x[4] + 1.0);
    Du1x[17] = D[2] * u0x[3] + D[5] * (u0x[4] + 1.0);
    Du1x[18] = B[1] * (u0x[4] + 1.0) + B[2] * u0x[3];
    Du1x[19] = B[3] * (u0x[4] + 1.0) + B[4] * u0x[3];
    Du1x[20] = B[4] * (u0x[4] + 1.0) + B[5] * u0x[3];
    Du1x[21] = D[1] * (u0x[4] + 1.0) + D[2] * u0x[3];
    Du1x[22] = D[3] * (u0x[4] + 1.0) + D[4] * u0x[3];
    Du1x[23] = D[4] * (u0x[4] + 1.0) + D[5] * u0x[3];
    Du1x[24] = B[0] * u0x[6] + B[2] * u0x[7];
    Du1x[25] = B[1] * u0x[6] + B[4] * u0x[7];
    Du1x[26] = B[2] * u0x[6] + B[5] * u0x[7];
    Du1x[27] = D[0] * u0x[6] + D[2] * u0x[7];
    Du1x[28] = D[1] * u0x[6] + D[4] * u0x[7];
    Du1x[29] = D[2] * u0x[6] + D[5] * u0x[7];
    Du1x[30] = B[1] * u0x[7] + B[2] * u0x[6];
    Du1x[31] = B[3] * u0x[7] + B[4] * u0x[6];
    Du1x[32] = B[4] * u0x[7] + B[5] * u0x[6];
    Du1x[33] = D[1] * u0x[7] + D[2] * u0x[6];
    Du1x[34] = D[3] * u0x[7] + D[4] * u0x[6];
    Du1x[35] = D[4] * u0x[7] + D[5] * u0x[6];

    TacsScalar Du1xd[36];
    Du1xd[0] = B[0] * u0xd[0] + B[2] * u0xd[1];
    Du1xd[1] = B[1] * u0xd[0] + B[4] * u0xd[1];
    Du1xd[2] = B[2] * u0xd[0] + B[5] * u0xd[1];
    Du1xd[3] = D[0] * u0xd[0] + D[2] * u0xd[1];
    Du1xd[4] = D[1] * u0xd[0] + D[4] * u0xd[1];
    Du1xd[5] = D[2] * u0xd[0] + D[5] * u0xd[1];
    Du1xd[6] = B[1] * u0xd[1] + B[2] * u0xd[0];
    Du1xd[7] = B[3] * u0xd[1] + B[4] * u0xd[0];
    Du1xd[8] = B[4] * u0xd[1] + B[5] * u0xd[0];
    Du1xd[9] = D[1] * u0xd[1] + D[2] * u0xd[0];
    Du1xd[10] = D[3] * u0xd[1] + D[4] * u0xd[0];
    Du1xd[11] = D[4] * u0xd[1] + D[5] * u0xd[0];
    Du1xd[12] = B[0] * u0xd[3] + B[2] * u0xd[4];
    Du1xd[13] = B[1] * u0xd[3] + B[4] * u0xd[4];
    Du1xd[14] = B[2] * u0xd[3] + B[5] * u0xd[4];
    Du1xd[15] = D[0] * u0xd[3] + D[2] * u0xd[4];
    Du1xd[16] = D[1] * u0xd[3] + D[4] * u0xd[4];
    Du1xd[17] = D[2] * u0xd[3] + D[5] * u0xd[4];
    Du1xd[18] = B[1] * u0xd[4] + B[2] * u0xd[3];
    Du1xd[19] = B[3] * u0xd[4] + B[4] * u0xd[3];
    Du1xd[20] = B[4] * u0xd[4] + B[5] * u0xd[3];
    Du1xd[21] = D[1] * u0xd[4] + D[2] * u0xd[3];
    Du1xd[22] = D[3] * u0xd[4] + D[4] * u0xd[3];
    Du1xd[23] = D[4] * u0xd[4] + D[5] * u0xd[3];
    Du1xd[24] = B[0] * u0xd[6] + B[2] * u0xd[7];
    Du1xd[25] = B[1] * u0xd[6] + B[4] * u0xd[7];
    Du1xd[26] = B[2] * u0xd[6] + B[5] * u0xd[7];
    Du1xd[27] = D[0] * u0xd[6] + D[2] * u0xd[7];
    Du1xd[28] = D[1] * u0xd[6] + D[4] * u0xd[7];
    Du1xd[29] = D[2] * u0xd[6] + D[5] * u0xd[7];
    Du1xd[30] = B[1] * u0xd[7] + B[2] * u0xd[6];
    Du1xd[31] = B[3] * u0xd[7] + B[4] * u0xd[6];
    Du1xd[32] = B[4] * u0xd[7] + B[5] * u0xd[6];
    Du1xd[33] = D[1] * u0xd[7] + D[2] * u0xd[6];
    Du1xd[34] = D[3] * u0xd[7] + D[4] * u0xd[6];
    Du1xd[35] = D[4] * u0xd[7] + D[5] * u0xd[6];

    d2u0x[0] = scale * (Du0x[0] * (u0x[0] + 1.0) + Du0x[2] * u0x[1] +
                        Du0x[3] * u1x[0] + Du0x[5] * u1x[1] + s[0]);
    d2u0x[1] = scale * (Du0x[1] * u0x[1] + Du0x[2] * (u0x[0] + 1.0) +
                        Du0x[4] * u1x[1] + Du0x[5] * u1x[0] + s[2]);
    d2u0x[2] = 0.0;
    d2u0x[3] = scale * (Du0x[0] * u0x[3] + Du0x[2] * (u0x[4] + 1.0) +
                        Du0x[3] * u1x[3] + Du0x[5] * u1x[4]);
    d2u0x[4] = scale * (Du0x[1] * (u0x[4] + 1.0) + Du0x[2] * u0x[3] +
                        Du0x[4] * u1x[4] + Du0x[5] * u1x[3]);
    d2u0x[5] = 0.0;
    d2u0x[6] = scale * (Du0x[0] * u0x[6] + Du0x[2] * u0x[7] + Du0x[3] * u1x[6] +
                        Du0x[5] * u1x[7]);
    d2u0x[7] = scale * (Du0x[1] * u0x[7] + Du0x[2] * u0x[6] + Du0x[4] * u1x[7] +
                        Du0x[5] * u1x[6]);
    d2u0x[8] = 0.0;

    d2u0x[9] = scale * (Du0x[11] * u1x[1] + Du0x[6] * (u0x[0] + 1.0) +
                        Du0x[8] * u0x[1] + Du0x[9] * u1x[0] + s[2]);
    d2u0x[10] = scale * (Du0x[10] * u1x[1] + Du0x[11] * u1x[0] +
                         Du0x[7] * u0x[1] + Du0x[8] * (u0x[0] + 1.0) + s[1]);
    d2u0x[11] = 0.0;
    d2u0x[12] = scale * (Du0x[11] * u1x[4] + Du0x[6] * u0x[3] +
                         Du0x[8] * (u0x[4] + 1.0) + Du0x[9] * u1x[3]);
    d2u0x[13] = scale * (Du0x[10] * u1x[4] + Du0x[11] * u1x[3] +
                         Du0x[7] * (u0x[4] + 1.0) + Du0x[8] * u0x[3]);
    d2u0x[14] = 0.0;
    d2u0x[15] = scale * (Du0x[11] * u1x[7] + Du0x[6] * u0x[6] +
                         Du0x[8] * u0x[7] + Du0x[9] * u1x[6]);
    d2u0x[16] = scale * (Du0x[10] * u1x[7] + Du0x[11] * u1x[6] +
                         Du0x[7] * u0x[7] + Du0x[8] * u0x[6]);
    d2u0x[17] = 0.0;

    d2u0x[18] = 0.0;
    d2u0x[19] = 0.0;
    d2u0x[20] = 0.0;
    d2u0x[21] = 0.0;
    d2u0x[22] = 0.0;
    d2u0x[23] = 0.0;
    d2u0x[24] = 0.0;
    d2u0x[25] = 0.0;
    d2u0x[26] = 0.0;

    d2u0x[27] = scale * (Du0x[12] * (u0x[0] + 1.0) + Du0x[14] * u0x[1] +
                         Du0x[15] * u1x[0] + Du0x[17] * u1x[1]);
    d2u0x[28] = scale * (Du0x[13] * u0x[1] + Du0x[14] * (u0x[0] + 1.0) +
                         Du0x[16] * u1x[1] + Du0x[17] * u1x[0]);
    d2u0x[29] = 0.0;
    d2u0x[30] = scale * (Du0x[12] * u0x[3] + Du0x[14] * (u0x[4] + 1.0) +
                         Du0x[15] * u1x[3] + Du0x[17] * u1x[4] + s[0]);
    d2u0x[31] = scale * (Du0x[13] * (u0x[4] + 1.0) + Du0x[14] * u0x[3] +
                         Du0x[16] * u1x[4] + Du0x[17] * u1x[3] + s[2]);
    d2u0x[32] = 0.0;
    d2u0x[33] = scale * (Du0x[12] * u0x[6] + Du0x[14] * u0x[7] +
                         Du0x[15] * u1x[6] + Du0x[17] * u1x[7]);
    d2u0x[34] = scale * (Du0x[13] * u0x[7] + Du0x[14] * u0x[6] +
                         Du0x[16] * u1x[7] + Du0x[17] * u1x[6]);
    d2u0x[35] = 0.0;

    d2u0x[36] = scale * (Du0x[18] * (u0x[0] + 1.0) + Du0x[20] * u0x[1] +
                         Du0x[21] * u1x[0] + Du0x[23] * u1x[1]);
    d2u0x[37] = scale * (Du0x[19] * u0x[1] + Du0x[20] * (u0x[0] + 1.0) +
                         Du0x[22] * u1x[1] + Du0x[23] * u1x[0]);
    d2u0x[38] = 0.0;
    d2u0x[39] = scale * (Du0x[18] * u0x[3] + Du0x[20] * (u0x[4] + 1.0) +
                         Du0x[21] * u1x[3] + Du0x[23] * u1x[4] + s[2]);
    d2u0x[40] = scale * (Du0x[19] * (u0x[4] + 1.0) + Du0x[20] * u0x[3] +
                         Du0x[22] * u1x[4] + Du0x[23] * u1x[3] + s[1]);
    d2u0x[41] = 0.0;
    d2u0x[42] = scale * (Du0x[18] * u0x[6] + Du0x[20] * u0x[7] +
                         Du0x[21] * u1x[6] + Du0x[23] * u1x[7]);
    d2u0x[43] = scale * (Du0x[19] * u0x[7] + Du0x[20] * u0x[6] +
                         Du0x[22] * u1x[7] + Du0x[23] * u1x[6]);
    d2u0x[44] = 0.0;

    d2u0x[45] = 0.0;
    d2u0x[46] = 0.0;
    d2u0x[47] = 0.0;
    d2u0x[48] = 0.0;
    d2u0x[49] = 0.0;
    d2u0x[50] = 0.0;
    d2u0x[51] = 0.0;
    d2u0x[52] = 0.0;
    d2u0x[53] = 0.0;

    d2u0x[54] = scale * (Du0x[24] * (u0x[0] + 1.0) + Du0x[26] * u0x[1] +
                         Du0x[27] * u1x[0] + Du0x[29] * u1x[1]);
    d2u0x[55] = scale * (Du0x[25] * u0x[1] + Du0x[26] * (u0x[0] + 1.0) +
                         Du0x[28] * u1x[1] + Du0x[29] * u1x[0]);
    d2u0x[56] = 0.0;
    d2u0x[57] = scale * (Du0x[24] * u0x[3] + Du0x[26] * (u0x[4] + 1.0) +
                         Du0x[27] * u1x[3] + Du0x[29] * u1x[4]);
    d2u0x[58] = scale * (Du0x[25] * (u0x[4] + 1.0) + Du0x[26] * u0x[3] +
                         Du0x[28] * u1x[4] + Du0x[29] * u1x[3]);
    d2u0x[59] = 0.0;
    d2u0x[60] = scale * (Du0x[24] * u0x[6] + Du0x[26] * u0x[7] +
                         Du0x[27] * u1x[6] + Du0x[29] * u1x[7] + s[0]);
    d2u0x[61] = scale * (Du0x[25] * u0x[7] + Du0x[26] * u0x[6] +
                         Du0x[28] * u1x[7] + Du0x[29] * u1x[6] + s[2]);
    d2u0x[62] = 0.0;

    d2u0x[63] = scale * (Du0x[30] * (u0x[0] + 1.0) + Du0x[32] * u0x[1] +
                         Du0x[33] * u1x[0] + Du0x[35] * u1x[1]);
    d2u0x[64] = scale * (Du0x[31] * u0x[1] + Du0x[32] * (u0x[0] + 1.0) +
                         Du0x[34] * u1x[1] + Du0x[35] * u1x[0]);
    d2u0x[65] = 0.0;
    d2u0x[66] = scale * (Du0x[30] * u0x[3] + Du0x[32] * (u0x[4] + 1.0) +
                         Du0x[33] * u1x[3] + Du0x[35] * u1x[4]);
    d2u0x[67] = scale * (Du0x[31] * (u0x[4] + 1.0) + Du0x[32] * u0x[3] +
                         Du0x[34] * u1x[4] + Du0x[35] * u1x[3]);
    d2u0x[68] = 0.0;
    d2u0x[69] = scale * (Du0x[30] * u0x[6] + Du0x[32] * u0x[7] +
                         Du0x[33] * u1x[6] + Du0x[35] * u1x[7] + s[2]);
    d2u0x[70] = scale * (Du0x[31] * u0x[7] + Du0x[32] * u0x[6] +
                         Du0x[34] * u1x[7] + Du0x[35] * u1x[6] + s[1]);
    d2u0x[71] = 0.0;

    d2u0x[72] = 0.0;
    d2u0x[73] = 0.0;
    d2u0x[74] = 0.0;
    d2u0x[75] = 0.0;
    d2u0x[76] = 0.0;
    d2u0x[77] = 0.0;
    d2u0x[78] = 0.0;
    d2u0x[79] = 0.0;
    d2u0x[80] = 0.0;

    d2u0xd[0] = scale * (Du0x[0] * u0xd[0] + Du0x[2] * u0xd[1] +
                         Du0x[3] * u1xd[0] + Du0x[5] * u1xd[1] + sd[0]);
    d2u0xd[1] = scale * (Du0x[1] * u0xd[1] + Du0x[2] * u0xd[0] +
                         Du0x[4] * u1xd[1] + Du0x[5] * u1xd[0] + sd[2]);
    d2u0xd[2] = 0.0;
    d2u0xd[3] = scale * (Du0x[0] * u0xd[3] + Du0x[2] * u0xd[4] +
                         Du0x[3] * u1xd[3] + Du0x[5] * u1xd[4]);
    d2u0xd[4] = scale * (Du0x[1] * u0xd[4] + Du0x[2] * u0xd[3] +
                         Du0x[4] * u1xd[4] + Du0x[5] * u1xd[3]);
    d2u0xd[5] = 0.0;
    d2u0xd[6] = scale * (Du0x[0] * u0xd[6] + Du0x[2] * u0xd[7] +
                         Du0x[3] * u1xd[6] + Du0x[5] * u1xd[7]);
    d2u0xd[7] = scale * (Du0x[1] * u0xd[7] + Du0x[2] * u0xd[6] +
                         Du0x[4] * u1xd[7] + Du0x[5] * u1xd[6]);
    d2u0xd[8] = 0.0;

    d2u0xd[9] = scale * (Du0x[11] * u1xd[1] + Du0x[6] * u0xd[0] +
                         Du0x[8] * u0xd[1] + Du0x[9] * u1xd[0] + sd[2]);
    d2u0xd[10] = scale * (Du0x[10] * u1xd[1] + Du0x[11] * u1xd[0] +
                          Du0x[7] * u0xd[1] + Du0x[8] * u0xd[0] + sd[1]);
    d2u0xd[11] = 0.0;
    d2u0xd[12] = scale * (Du0x[11] * u1xd[4] + Du0x[6] * u0xd[3] +
                          Du0x[8] * u0xd[4] + Du0x[9] * u1xd[3]);
    d2u0xd[13] = scale * (Du0x[10] * u1xd[4] + Du0x[11] * u1xd[3] +
                          Du0x[7] * u0xd[4] + Du0x[8] * u0xd[3]);
    d2u0xd[14] = 0.0;
    d2u0xd[15] = scale * (Du0x[11] * u1xd[7] + Du0x[6] * u0xd[6] +
                          Du0x[8] * u0xd[7] + Du0x[9] * u1xd[6]);
    d2u0xd[16] = scale * (Du0x[10] * u1xd[7] + Du0x[11] * u1xd[6] +
                          Du0x[7] * u0xd[7] + Du0x[8] * u0xd[6]);
    d2u0xd[17] = 0.0;

    d2u0xd[18] = 0.0;
    d2u0xd[19] = 0.0;
    d2u0xd[20] = 0.0;
    d2u0xd[21] = 0.0;
    d2u0xd[22] = 0.0;
    d2u0xd[23] = 0.0;
    d2u0xd[24] = 0.0;
    d2u0xd[25] = 0.0;
    d2u0xd[26] = 0.0;

    d2u0xd[27] = scale * (Du0x[12] * u0xd[0] + Du0x[14] * u0xd[1] +
                          Du0x[15] * u1xd[0] + Du0x[17] * u1xd[1]);
    d2u0xd[28] = scale * (Du0x[13] * u0xd[1] + Du0x[14] * u0xd[0] +
                          Du0x[16] * u1xd[1] + Du0x[17] * u1xd[0]);
    d2u0xd[29] = 0.0;
    d2u0xd[30] = scale * (Du0x[12] * u0xd[3] + Du0x[14] * u0xd[4] +
                          Du0x[15] * u1xd[3] + Du0x[17] * u1xd[4] + sd[0]);
    d2u0xd[31] = scale * (Du0x[13] * u0xd[4] + Du0x[14] * u0xd[3] +
                          Du0x[16] * u1xd[4] + Du0x[17] * u1xd[3] + sd[2]);
    d2u0xd[32] = 0.0;
    d2u0xd[33] = scale * (Du0x[12] * u0xd[6] + Du0x[14] * u0xd[7] +
                          Du0x[15] * u1xd[6] + Du0x[17] * u1xd[7]);
    d2u0xd[34] = scale * (Du0x[13] * u0xd[7] + Du0x[14] * u0xd[6] +
                          Du0x[16] * u1xd[7] + Du0x[17] * u1xd[6]);
    d2u0xd[35] = 0.0;

    d2u0xd[36] = scale * (Du0x[18] * u0xd[0] + Du0x[20] * u0xd[1] +
                          Du0x[21] * u1xd[0] + Du0x[23] * u1xd[1]);
    d2u0xd[37] = scale * (Du0x[19] * u0xd[1] + Du0x[20] * u0xd[0] +
                          Du0x[22] * u1xd[1] + Du0x[23] * u1xd[0]);
    d2u0xd[38] = 0.0;
    d2u0xd[39] = scale * (Du0x[18] * u0xd[3] + Du0x[20] * u0xd[4] +
                          Du0x[21] * u1xd[3] + Du0x[23] * u1xd[4] + sd[2]);
    d2u0xd[40] = scale * (Du0x[19] * u0xd[4] + Du0x[20] * u0xd[3] +
                          Du0x[22] * u1xd[4] + Du0x[23] * u1xd[3] + sd[1]);
    d2u0xd[41] = 0.0;
    d2u0xd[42] = scale * (Du0x[18] * u0xd[6] + Du0x[20] * u0xd[7] +
                          Du0x[21] * u1xd[6] + Du0x[23] * u1xd[7]);
    d2u0xd[43] = scale * (Du0x[19] * u0xd[7] + Du0x[20] * u0xd[6] +
                          Du0x[22] * u1xd[7] + Du0x[23] * u1xd[6]);
    d2u0xd[44] = 0.0;

    d2u0xd[45] = 0.0;
    d2u0xd[46] = 0.0;
    d2u0xd[47] = 0.0;
    d2u0xd[48] = 0.0;
    d2u0xd[49] = 0.0;
    d2u0xd[50] = 0.0;
    d2u0xd[51] = 0.0;
    d2u0xd[52] = 0.0;
    d2u0xd[53] = 0.0;

    d2u0xd[54] = scale * (Du0x[24] * u0xd[0] + Du0x[26] * u0xd[1] +
                          Du0x[27] * u1xd[0] + Du0x[29] * u1xd[1]);
    d2u0xd[55] = scale * (Du0x[25] * u0xd[1] + Du0x[26] * u0xd[0] +
                          Du0x[28] * u1xd[1] + Du0x[29] * u1xd[0]);
    d2u0xd[56] = 0.0;
    d2u0xd[57] = scale * (Du0x[24] * u0xd[3] + Du0x[26] * u0xd[4] +
                          Du0x[27] * u1xd[3] + Du0x[29] * u1xd[4]);
    d2u0xd[58] = scale * (Du0x[25] * u0xd[4] + Du0x[26] * u0xd[3] +
                          Du0x[28] * u1xd[4] + Du0x[29] * u1xd[3]);
    d2u0xd[59] = 0.0;
    d2u0xd[60] = scale * (Du0x[24] * u0xd[6] + Du0x[26] * u0xd[7] +
                          Du0x[27] * u1xd[6] + Du0x[29] * u1xd[7] + sd[0]);
    d2u0xd[61] = scale * (Du0x[25] * u0xd[7] + Du0x[26] * u0xd[6] +
                          Du0x[28] * u1xd[7] + Du0x[29] * u1xd[6] + sd[2]);
    d2u0xd[62] = 0.0;

    d2u0xd[63] = scale * (Du0x[30] * u0xd[0] + Du0x[32] * u0xd[1] +
                          Du0x[33] * u1xd[0] + Du0x[35] * u1xd[1]);
    d2u0xd[64] = scale * (Du0x[31] * u0xd[1] + Du0x[32] * u0xd[0] +
                          Du0x[34] * u1xd[1] + Du0x[35] * u1xd[0]);
    d2u0xd[65] = 0.0;
    d2u0xd[66] = scale * (Du0x[30] * u0xd[3] + Du0x[32] * u0xd[4] +
                          Du0x[33] * u1xd[3] + Du0x[35] * u1xd[4]);
    d2u0xd[67] = scale * (Du0x[31] * u0xd[4] + Du0x[32] * u0xd[3] +
                          Du0x[34] * u1xd[4] + Du0x[35] * u1xd[3]);
    d2u0xd[68] = 0.0;
    d2u0xd[69] = scale * (Du0x[30] * u0xd[6] + Du0x[32] * u0xd[7] +
                          Du0x[33] * u1xd[6] + Du0x[35] * u1xd[7] + sd[2]);
    d2u0xd[70] = scale * (Du0x[31] * u0xd[7] + Du0x[32] * u0xd[6] +
                          Du0x[34] * u1xd[7] + Du0x[35] * u1xd[6] + sd[1]);
    d2u0xd[71] = 0.0;

    d2u0xd[72] = 0.0;
    d2u0xd[73] = 0.0;
    d2u0xd[74] = 0.0;
    d2u0xd[75] = 0.0;
    d2u0xd[76] = 0.0;
    d2u0xd[77] = 0.0;
    d2u0xd[78] = 0.0;
    d2u0xd[79] = 0.0;
    d2u0xd[80] = 0.0;

    d2u0xd[0] += scale * (Du0xd[0] * (u0x[0] + 1.0) + Du0xd[2] * u0x[1] +
                          Du0xd[3] * u1x[0] + Du0xd[5] * u1x[1]);
    d2u0xd[1] += scale * (Du0xd[1] * u0x[1] + Du0xd[2] * (u0x[0] + 1.0) +
                          Du0xd[4] * u1x[1] + Du0xd[5] * u1x[0]);
    d2u0xd[3] += scale * (Du0xd[0] * u0x[3] + Du0xd[2] * (u0x[4] + 1.0) +
                          Du0xd[3] * u1x[3] + Du0xd[5] * u1x[4]);
    d2u0xd[4] += scale * (Du0xd[1] * (u0x[4] + 1.0) + Du0xd[2] * u0x[3] +
                          Du0xd[4] * u1x[4] + Du0xd[5] * u1x[3]);
    d2u0xd[6] += scale * (Du0xd[0] * u0x[6] + Du0xd[2] * u0x[7] +
                          Du0xd[3] * u1x[6] + Du0xd[5] * u1x[7]);
    d2u0xd[7] += scale * (Du0xd[1] * u0x[7] + Du0xd[2] * u0x[6] +
                          Du0xd[4] * u1x[7] + Du0xd[5] * u1x[6]);

    d2u0xd[9] += scale * (Du0xd[11] * u1x[1] + Du0xd[6] * (u0x[0] + 1.0) +
                          Du0xd[8] * u0x[1] + Du0xd[9] * u1x[0]);
    d2u0xd[10] += scale * (Du0xd[10] * u1x[1] + Du0xd[11] * u1x[0] +
                           Du0xd[7] * u0x[1] + Du0xd[8] * (u0x[0] + 1.0));
    d2u0xd[12] += scale * (Du0xd[11] * u1x[4] + Du0xd[6] * u0x[3] +
                           Du0xd[8] * (u0x[4] + 1.0) + Du0xd[9] * u1x[3]);
    d2u0xd[13] += scale * (Du0xd[10] * u1x[4] + Du0xd[11] * u1x[3] +
                           Du0xd[7] * (u0x[4] + 1.0) + Du0xd[8] * u0x[3]);
    d2u0xd[15] += scale * (Du0xd[11] * u1x[7] + Du0xd[6] * u0x[6] +
                           Du0xd[8] * u0x[7] + Du0xd[9] * u1x[6]);
    d2u0xd[16] += scale * (Du0xd[10] * u1x[7] + Du0xd[11] * u1x[6] +
                           Du0xd[7] * u0x[7] + Du0xd[8] * u0x[6]);

    d2u0xd[27] += scale * (Du0xd[12] * (u0x[0] + 1.0) + Du0xd[14] * u0x[1] +
                           Du0xd[15] * u1x[0] + Du0xd[17] * u1x[1]);
    d2u0xd[28] += scale * (Du0xd[13] * u0x[1] + Du0xd[14] * (u0x[0] + 1.0) +
                           Du0xd[16] * u1x[1] + Du0xd[17] * u1x[0]);
    d2u0xd[30] += scale * (Du0xd[12] * u0x[3] + Du0xd[14] * (u0x[4] + 1.0) +
                           Du0xd[15] * u1x[3] + Du0xd[17] * u1x[4]);
    d2u0xd[31] += scale * (Du0xd[13] * (u0x[4] + 1.0) + Du0xd[14] * u0x[3] +
                           Du0xd[16] * u1x[4] + Du0xd[17] * u1x[3]);
    d2u0xd[33] += scale * (Du0xd[12] * u0x[6] + Du0xd[14] * u0x[7] +
                           Du0xd[15] * u1x[6] + Du0xd[17] * u1x[7]);
    d2u0xd[34] += scale * (Du0xd[13] * u0x[7] + Du0xd[14] * u0x[6] +
                           Du0xd[16] * u1x[7] + Du0xd[17] * u1x[6]);

    d2u0xd[36] += scale * (Du0xd[18] * (u0x[0] + 1.0) + Du0xd[20] * u0x[1] +
                           Du0xd[21] * u1x[0] + Du0xd[23] * u1x[1]);
    d2u0xd[37] += scale * (Du0xd[19] * u0x[1] + Du0xd[20] * (u0x[0] + 1.0) +
                           Du0xd[22] * u1x[1] + Du0xd[23] * u1x[0]);
    d2u0xd[39] += scale * (Du0xd[18] * u0x[3] + Du0xd[20] * (u0x[4] + 1.0) +
                           Du0xd[21] * u1x[3] + Du0xd[23] * u1x[4]);
    d2u0xd[40] += scale * (Du0xd[19] * (u0x[4] + 1.0) + Du0xd[20] * u0x[3] +
                           Du0xd[22] * u1x[4] + Du0xd[23] * u1x[3]);
    d2u0xd[42] += scale * (Du0xd[18] * u0x[6] + Du0xd[20] * u0x[7] +
                           Du0xd[21] * u1x[6] + Du0xd[23] * u1x[7]);
    d2u0xd[43] += scale * (Du0xd[19] * u0x[7] + Du0xd[20] * u0x[6] +
                           Du0xd[22] * u1x[7] + Du0xd[23] * u1x[6]);

    d2u0xd[54] += scale * (Du0xd[24] * (u0x[0] + 1.0) + Du0xd[26] * u0x[1] +
                           Du0xd[27] * u1x[0] + Du0xd[29] * u1x[1]);
    d2u0xd[55] += scale * (Du0xd[25] * u0x[1] + Du0xd[26] * (u0x[0] + 1.0) +
                           Du0xd[28] * u1x[1] + Du0xd[29] * u1x[0]);
    d2u0xd[57] += scale * (Du0xd[24] * u0x[3] + Du0xd[26] * (u0x[4] + 1.0) +
                           Du0xd[27] * u1x[3] + Du0xd[29] * u1x[4]);
    d2u0xd[58] += scale * (Du0xd[25] * (u0x[4] + 1.0) + Du0xd[26] * u0x[3] +
                           Du0xd[28] * u1x[4] + Du0xd[29] * u1x[3]);
    d2u0xd[60] += scale * (Du0xd[24] * u0x[6] + Du0xd[26] * u0x[7] +
                           Du0xd[27] * u1x[6] + Du0xd[29] * u1x[7]);
    d2u0xd[61] += scale * (Du0xd[25] * u0x[7] + Du0xd[26] * u0x[6] +
                           Du0xd[28] * u1x[7] + Du0xd[29] * u1x[6]);
    d2u0xd[63] += scale * (Du0xd[30] * (u0x[0] + 1.0) + Du0xd[32] * u0x[1] +
                           Du0xd[33] * u1x[0] + Du0xd[35] * u1x[1]);
    d2u0xd[64] += scale * (Du0xd[31] * u0x[1] + Du0xd[32] * (u0x[0] + 1.0) +
                           Du0xd[34] * u1x[1] + Du0xd[35] * u1x[0]);
    d2u0xd[66] += scale * (Du0xd[30] * u0x[3] + Du0xd[32] * (u0x[4] + 1.0) +
                           Du0xd[33] * u1x[3] + Du0xd[35] * u1x[4]);
    d2u0xd[67] += scale * (Du0xd[31] * (u0x[4] + 1.0) + Du0xd[32] * u0x[3] +
                           Du0xd[34] * u1x[4] + Du0xd[35] * u1x[3]);
    d2u0xd[69] += scale * (Du0xd[30] * u0x[6] + Du0xd[32] * u0x[7] +
                           Du0xd[33] * u1x[6] + Du0xd[35] * u1x[7]);
    d2u0xd[70] += scale * (Du0xd[31] * u0x[7] + Du0xd[32] * u0x[6] +
                           Du0xd[34] * u1x[7] + Du0xd[35] * u1x[6]);

    d2u0xu1x[0] = scale * (Du0x[3] * (u0x[0] + 1.0) + Du0x[5] * u0x[1] + s[3]);
    d2u0xu1x[1] = scale * (Du0x[4] * u0x[1] + Du0x[5] * (u0x[0] + 1.0) + s[5]);
    d2u0xu1x[2] = 0.0;
    d2u0xu1x[3] = scale * (Du0x[3] * u0x[3] + Du0x[5] * (u0x[4] + 1.0));
    d2u0xu1x[4] = scale * (Du0x[4] * (u0x[4] + 1.0) + Du0x[5] * u0x[3]);
    d2u0xu1x[5] = 0.0;
    d2u0xu1x[6] = scale * (Du0x[3] * u0x[6] + Du0x[5] * u0x[7]);
    d2u0xu1x[7] = scale * (Du0x[4] * u0x[7] + Du0x[5] * u0x[6]);
    d2u0xu1x[8] = 0.0;

    d2u0xu1x[9] = scale * (Du0x[11] * u0x[1] + Du0x[9] * (u0x[0] + 1.0) + s[5]);
    d2u0xu1x[10] =
        scale * (Du0x[10] * u0x[1] + Du0x[11] * (u0x[0] + 1.0) + s[4]);
    d2u0xu1x[11] = 0.0;
    d2u0xu1x[12] = scale * (Du0x[11] * (u0x[4] + 1.0) + Du0x[9] * u0x[3]);
    d2u0xu1x[13] = scale * (Du0x[10] * (u0x[4] + 1.0) + Du0x[11] * u0x[3]);
    d2u0xu1x[14] = 0.0;
    d2u0xu1x[15] = scale * (Du0x[11] * u0x[7] + Du0x[9] * u0x[6]);
    d2u0xu1x[16] = scale * (Du0x[10] * u0x[7] + Du0x[11] * u0x[6]);
    d2u0xu1x[17] = 0.0;

    d2u0xu1x[18] = 0.0;
    d2u0xu1x[19] = 0.0;
    d2u0xu1x[20] = 0.0;
    d2u0xu1x[21] = 0.0;
    d2u0xu1x[22] = 0.0;
    d2u0xu1x[23] = 0.0;
    d2u0xu1x[24] = 0.0;
    d2u0xu1x[25] = 0.0;
    d2u0xu1x[26] = 0.0;

    d2u0xu1x[27] = scale * (Du0x[15] * (u0x[0] + 1.0) + Du0x[17] * u0x[1]);
    d2u0xu1x[28] = scale * (Du0x[16] * u0x[1] + Du0x[17] * (u0x[0] + 1.0));
    d2u0xu1x[29] = 0.0;
    d2u0xu1x[30] =
        scale * (Du0x[15] * u0x[3] + Du0x[17] * (u0x[4] + 1.0) + s[3]);
    d2u0xu1x[31] =
        scale * (Du0x[16] * (u0x[4] + 1.0) + Du0x[17] * u0x[3] + s[5]);
    d2u0xu1x[32] = 0.0;
    d2u0xu1x[33] = scale * (Du0x[15] * u0x[6] + Du0x[17] * u0x[7]);
    d2u0xu1x[34] = scale * (Du0x[16] * u0x[7] + Du0x[17] * u0x[6]);
    d2u0xu1x[35] = 0.0;

    d2u0xu1x[36] = scale * (Du0x[21] * (u0x[0] + 1.0) + Du0x[23] * u0x[1]);
    d2u0xu1x[37] = scale * (Du0x[22] * u0x[1] + Du0x[23] * (u0x[0] + 1.0));
    d2u0xu1x[38] = 0.0;
    d2u0xu1x[39] =
        scale * (Du0x[21] * u0x[3] + Du0x[23] * (u0x[4] + 1.0) + s[5]);
    d2u0xu1x[40] =
        scale * (Du0x[22] * (u0x[4] + 1.0) + Du0x[23] * u0x[3] + s[4]);
    d2u0xu1x[41] = 0.0;
    d2u0xu1x[42] = scale * (Du0x[21] * u0x[6] + Du0x[23] * u0x[7]);
    d2u0xu1x[43] = scale * (Du0x[22] * u0x[7] + Du0x[23] * u0x[6]);
    d2u0xu1x[44] = 0.0;

    d2u0xu1x[45] = 0.0;
    d2u0xu1x[46] = 0.0;
    d2u0xu1x[47] = 0.0;
    d2u0xu1x[48] = 0.0;
    d2u0xu1x[49] = 0.0;
    d2u0xu1x[50] = 0.0;
    d2u0xu1x[51] = 0.0;
    d2u0xu1x[52] = 0.0;
    d2u0xu1x[53] = 0.0;

    d2u0xu1x[54] = scale * (Du0x[27] * (u0x[0] + 1.0) + Du0x[29] * u0x[1]);
    d2u0xu1x[55] = scale * (Du0x[28] * u0x[1] + Du0x[29] * (u0x[0] + 1.0));
    d2u0xu1x[56] = 0.0;
    d2u0xu1x[57] = scale * (Du0x[27] * u0x[3] + Du0x[29] * (u0x[4] + 1.0));
    d2u0xu1x[58] = scale * (Du0x[28] * (u0x[4] + 1.0) + Du0x[29] * u0x[3]);
    d2u0xu1x[59] = 0.0;
    d2u0xu1x[60] = scale * (Du0x[27] * u0x[6] + Du0x[29] * u0x[7] + s[3]);
    d2u0xu1x[61] = scale * (Du0x[28] * u0x[7] + Du0x[29] * u0x[6] + s[5]);
    d2u0xu1x[62] = 0.0;

    d2u0xu1x[63] = scale * (Du0x[33] * (u0x[0] + 1.0) + Du0x[35] * u0x[1]);
    d2u0xu1x[64] = scale * (Du0x[34] * u0x[1] + Du0x[35] * (u0x[0] + 1.0));
    d2u0xu1x[65] = 0.0;
    d2u0xu1x[66] = scale * (Du0x[33] * u0x[3] + Du0x[35] * (u0x[4] + 1.0));
    d2u0xu1x[67] = scale * (Du0x[34] * (u0x[4] + 1.0) + Du0x[35] * u0x[3]);
    d2u0xu1x[68] = 0.0;
    d2u0xu1x[69] = scale * (Du0x[33] * u0x[6] + Du0x[35] * u0x[7] + s[5]);
    d2u0xu1x[70] = scale * (Du0x[34] * u0x[7] + Du0x[35] * u0x[6] + s[4]);
    d2u0xu1x[71] = 0.0;

    d2u0xu1x[72] = 0.0;
    d2u0xu1x[73] = 0.0;
    d2u0xu1x[74] = 0.0;
    d2u0xu1x[75] = 0.0;
    d2u0xu1x[76] = 0.0;
    d2u0xu1x[77] = 0.0;
    d2u0xu1x[78] = 0.0;
    d2u0xu1x[79] = 0.0;
    d2u0xu1x[80] = 0.0;

    d2u0xu1xd[0] = scale * (Du0x[3] * u0xd[0] + Du0x[5] * u0xd[1] + sd[3]);
    d2u0xu1xd[1] = scale * (Du0x[4] * u0xd[1] + Du0x[5] * u0xd[0] + sd[5]);
    d2u0xu1xd[2] = 0.0;
    d2u0xu1xd[3] = scale * (Du0x[3] * u0xd[3] + Du0x[5] * u0xd[4]);
    d2u0xu1xd[4] = scale * (Du0x[4] * u0xd[4] + Du0x[5] * u0xd[3]);
    d2u0xu1xd[5] = 0.0;
    d2u0xu1xd[6] = scale * (Du0x[3] * u0xd[6] + Du0x[5] * u0xd[7]);
    d2u0xu1xd[7] = scale * (Du0x[4] * u0xd[7] + Du0x[5] * u0xd[6]);
    d2u0xu1xd[8] = 0.0;

    d2u0xu1xd[9] = scale * (Du0x[11] * u0xd[1] + Du0x[9] * u0xd[0] + sd[5]);
    d2u0xu1xd[10] = scale * (Du0x[10] * u0xd[1] + Du0x[11] * u0xd[0] + sd[4]);
    d2u0xu1xd[11] = 0.0;
    d2u0xu1xd[12] = scale * (Du0x[11] * u0xd[4] + Du0x[9] * u0xd[3]);
    d2u0xu1xd[13] = scale * (Du0x[10] * u0xd[4] + Du0x[11] * u0xd[3]);
    d2u0xu1xd[14] = 0.0;
    d2u0xu1xd[15] = scale * (Du0x[11] * u0xd[7] + Du0x[9] * u0xd[6]);
    d2u0xu1xd[16] = scale * (Du0x[10] * u0xd[7] + Du0x[11] * u0xd[6]);
    d2u0xu1xd[17] = 0.0;

    d2u0xu1xd[18] = 0.0;
    d2u0xu1xd[19] = 0.0;
    d2u0xu1xd[20] = 0.0;
    d2u0xu1xd[21] = 0.0;
    d2u0xu1xd[22] = 0.0;
    d2u0xu1xd[23] = 0.0;
    d2u0xu1xd[24] = 0.0;
    d2u0xu1xd[25] = 0.0;
    d2u0xu1xd[26] = 0.0;

    d2u0xu1xd[27] = scale * (Du0x[15] * u0xd[0] + Du0x[17] * u0xd[1]);
    d2u0xu1xd[28] = scale * (Du0x[16] * u0xd[1] + Du0x[17] * u0xd[0]);
    d2u0xu1xd[29] = 0.0;
    d2u0xu1xd[30] = scale * (Du0x[15] * u0xd[3] + Du0x[17] * u0xd[4] + sd[3]);
    d2u0xu1xd[31] = scale * (Du0x[16] * u0xd[4] + Du0x[17] * u0xd[3] + sd[5]);
    d2u0xu1xd[32] = 0.0;
    d2u0xu1xd[33] = scale * (Du0x[15] * u0xd[6] + Du0x[17] * u0xd[7]);
    d2u0xu1xd[34] = scale * (Du0x[16] * u0xd[7] + Du0x[17] * u0xd[6]);
    d2u0xu1xd[35] = 0.0;

    d2u0xu1xd[36] = scale * (Du0x[21] * u0xd[0] + Du0x[23] * u0xd[1]);
    d2u0xu1xd[37] = scale * (Du0x[22] * u0xd[1] + Du0x[23] * u0xd[0]);
    d2u0xu1xd[38] = 0.0;
    d2u0xu1xd[39] = scale * (Du0x[21] * u0xd[3] + Du0x[23] * u0xd[4] + sd[5]);
    d2u0xu1xd[40] = scale * (Du0x[22] * u0xd[4] + Du0x[23] * u0xd[3] + sd[4]);
    d2u0xu1xd[41] = 0.0;
    d2u0xu1xd[42] = scale * (Du0x[21] * u0xd[6] + Du0x[23] * u0xd[7]);
    d2u0xu1xd[43] = scale * (Du0x[22] * u0xd[7] + Du0x[23] * u0xd[6]);
    d2u0xu1xd[44] = 0.0;

    d2u0xu1xd[45] = 0.0;
    d2u0xu1xd[46] = 0.0;
    d2u0xu1xd[47] = 0.0;
    d2u0xu1xd[48] = 0.0;
    d2u0xu1xd[49] = 0.0;
    d2u0xu1xd[50] = 0.0;
    d2u0xu1xd[51] = 0.0;
    d2u0xu1xd[52] = 0.0;
    d2u0xu1xd[53] = 0.0;

    d2u0xu1xd[54] = scale * (Du0x[27] * u0xd[0] + Du0x[29] * u0xd[1]);
    d2u0xu1xd[55] = scale * (Du0x[28] * u0xd[1] + Du0x[29] * u0xd[0]);
    d2u0xu1xd[56] = 0.0;
    d2u0xu1xd[57] = scale * (Du0x[27] * u0xd[3] + Du0x[29] * u0xd[4]);
    d2u0xu1xd[58] = scale * (Du0x[28] * u0xd[4] + Du0x[29] * u0xd[3]);
    d2u0xu1xd[59] = 0.0;
    d2u0xu1xd[60] = scale * (Du0x[27] * u0xd[6] + Du0x[29] * u0xd[7] + sd[3]);
    d2u0xu1xd[61] = scale * (Du0x[28] * u0xd[7] + Du0x[29] * u0xd[6] + sd[5]);
    d2u0xu1xd[62] = 0.0;

    d2u0xu1xd[63] = scale * (Du0x[33] * u0xd[0] + Du0x[35] * u0xd[1]);
    d2u0xu1xd[64] = scale * (Du0x[34] * u0xd[1] + Du0x[35] * u0xd[0]);
    d2u0xu1xd[65] = 0.0;
    d2u0xu1xd[66] = scale * (Du0x[33] * u0xd[3] + Du0x[35] * u0xd[4]);
    d2u0xu1xd[67] = scale * (Du0x[34] * u0xd[4] + Du0x[35] * u0xd[3]);
    d2u0xu1xd[68] = 0.0;
    d2u0xu1xd[69] = scale * (Du0x[33] * u0xd[6] + Du0x[35] * u0xd[7] + sd[5]);
    d2u0xu1xd[70] = scale * (Du0x[34] * u0xd[7] + Du0x[35] * u0xd[6] + sd[4]);
    d2u0xu1xd[71] = 0.0;

    d2u0xu1xd[72] = 0.0;
    d2u0xu1xd[73] = 0.0;
    d2u0xu1xd[74] = 0.0;
    d2u0xu1xd[75] = 0.0;
    d2u0xu1xd[76] = 0.0;
    d2u0xu1xd[77] = 0.0;
    d2u0xu1xd[78] = 0.0;
    d2u0xu1xd[79] = 0.0;
    d2u0xu1xd[80] = 0.0;

    d2u0xu1xd[0] += scale * (Du0xd[3] * (u0x[0] + 1.0) + Du0xd[5] * u0x[1]);
    d2u0xu1xd[1] += scale * (Du0xd[4] * u0x[1] + Du0xd[5] * (u0x[0] + 1.0));
    d2u0xu1xd[3] += scale * (Du0xd[3] * u0x[3] + Du0xd[5] * (u0x[4] + 1.0));
    d2u0xu1xd[4] += scale * (Du0xd[4] * (u0x[4] + 1.0) + Du0xd[5] * u0x[3]);
    d2u0xu1xd[6] += scale * (Du0xd[3] * u0x[6] + Du0xd[5] * u0x[7]);
    d2u0xu1xd[7] += scale * (Du0xd[4] * u0x[7] + Du0xd[5] * u0x[6]);

    d2u0xu1xd[9] += scale * (Du0xd[11] * u0x[1] + Du0xd[9] * (u0x[0] + 1.0));
    d2u0xu1xd[10] += scale * (Du0xd[10] * u0x[1] + Du0xd[11] * (u0x[0] + 1.0));
    d2u0xu1xd[12] += scale * (Du0xd[11] * (u0x[4] + 1.0) + Du0xd[9] * u0x[3]);
    d2u0xu1xd[13] += scale * (Du0xd[10] * (u0x[4] + 1.0) + Du0xd[11] * u0x[3]);
    d2u0xu1xd[15] += scale * (Du0xd[11] * u0x[7] + Du0xd[9] * u0x[6]);
    d2u0xu1xd[16] += scale * (Du0xd[10] * u0x[7] + Du0xd[11] * u0x[6]);

    d2u0xu1xd[27] += scale * (Du0xd[15] * (u0x[0] + 1.0) + Du0xd[17] * u0x[1]);
    d2u0xu1xd[28] += scale * (Du0xd[16] * u0x[1] + Du0xd[17] * (u0x[0] + 1.0));
    d2u0xu1xd[30] += scale * (Du0xd[15] * u0x[3] + Du0xd[17] * (u0x[4] + 1.0));
    d2u0xu1xd[31] += scale * (Du0xd[16] * (u0x[4] + 1.0) + Du0xd[17] * u0x[3]);
    d2u0xu1xd[33] += scale * (Du0xd[15] * u0x[6] + Du0xd[17] * u0x[7]);
    d2u0xu1xd[34] += scale * (Du0xd[16] * u0x[7] + Du0xd[17] * u0x[6]);

    d2u0xu1xd[36] += scale * (Du0xd[21] * (u0x[0] + 1.0) + Du0xd[23] * u0x[1]);
    d2u0xu1xd[37] += scale * (Du0xd[22] * u0x[1] + Du0xd[23] * (u0x[0] + 1.0));
    d2u0xu1xd[39] += scale * (Du0xd[21] * u0x[3] + Du0xd[23] * (u0x[4] + 1.0));
    d2u0xu1xd[40] += scale * (Du0xd[22] * (u0x[4] + 1.0) + Du0xd[23] * u0x[3]);
    d2u0xu1xd[42] += scale * (Du0xd[21] * u0x[6] + Du0xd[23] * u0x[7]);
    d2u0xu1xd[43] += scale * (Du0xd[22] * u0x[7] + Du0xd[23] * u0x[6]);

    d2u0xu1xd[54] += scale * (Du0xd[27] * (u0x[0] + 1.0) + Du0xd[29] * u0x[1]);
    d2u0xu1xd[55] += scale * (Du0xd[28] * u0x[1] + Du0xd[29] * (u0x[0] + 1.0));
    d2u0xu1xd[57] += scale * (Du0xd[27] * u0x[3] + Du0xd[29] * (u0x[4] + 1.0));
    d2u0xu1xd[58] += scale * (Du0xd[28] * (u0x[4] + 1.0) + Du0xd[29] * u0x[3]);
    d2u0xu1xd[60] += scale * (Du0xd[27] * u0x[6] + Du0xd[29] * u0x[7]);
    d2u0xu1xd[61] += scale * (Du0xd[28] * u0x[7] + Du0xd[29] * u0x[6]);

    d2u0xu1xd[63] += scale * (Du0xd[33] * (u0x[0] + 1.0) + Du0xd[35] * u0x[1]);
    d2u0xu1xd[64] += scale * (Du0xd[34] * u0x[1] + Du0xd[35] * (u0x[0] + 1.0));
    d2u0xu1xd[65] += 0.0;
    d2u0xu1xd[66] += scale * (Du0xd[33] * u0x[3] + Du0xd[35] * (u0x[4] + 1.0));
    d2u0xu1xd[67] += scale * (Du0xd[34] * (u0x[4] + 1.0) + Du0xd[35] * u0x[3]);
    d2u0xu1xd[68] += 0.0;
    d2u0xu1xd[69] += scale * (Du0xd[33] * u0x[6] + Du0xd[35] * u0x[7]);
    d2u0xu1xd[70] += scale * (Du0xd[34] * u0x[7] + Du0xd[35] * u0x[6]);

    d2u1x[0] = scale * (Du1x[3] * (u0x[0] + 1.0) + Du1x[5] * u0x[1]);
    d2u1x[1] = scale * (Du1x[4] * u0x[1] + Du1x[5] * (u0x[0] + 1.0));
    d2u1x[2] = 0.0;
    d2u1x[3] = scale * (Du1x[3] * u0x[3] + Du1x[5] * (u0x[4] + 1.0));
    d2u1x[4] = scale * (Du1x[4] * (u0x[4] + 1.0) + Du1x[5] * u0x[3]);
    d2u1x[5] = 0.0;
    d2u1x[6] = scale * (Du1x[3] * u0x[6] + Du1x[5] * u0x[7]);
    d2u1x[7] = scale * (Du1x[4] * u0x[7] + Du1x[5] * u0x[6]);
    d2u1x[8] = 0.0;

    d2u1x[9] = scale * (Du1x[11] * u0x[1] + Du1x[9] * (u0x[0] + 1.0));
    d2u1x[10] = scale * (Du1x[10] * u0x[1] + Du1x[11] * (u0x[0] + 1.0));
    d2u1x[11] = 0.0;
    d2u1x[12] = scale * (Du1x[11] * (u0x[4] + 1.0) + Du1x[9] * u0x[3]);
    d2u1x[13] = scale * (Du1x[10] * (u0x[4] + 1.0) + Du1x[11] * u0x[3]);
    d2u1x[14] = 0.0;
    d2u1x[15] = scale * (Du1x[11] * u0x[7] + Du1x[9] * u0x[6]);
    d2u1x[16] = scale * (Du1x[10] * u0x[7] + Du1x[11] * u0x[6]);
    d2u1x[17] = 0.0;

    d2u1x[18] = 0.0;
    d2u1x[19] = 0.0;
    d2u1x[20] = 0.0;
    d2u1x[21] = 0.0;
    d2u1x[22] = 0.0;
    d2u1x[23] = 0.0;
    d2u1x[24] = 0.0;
    d2u1x[25] = 0.0;
    d2u1x[26] = 0.0;

    d2u1x[27] = scale * (Du1x[15] * (u0x[0] + 1.0) + Du1x[17] * u0x[1]);
    d2u1x[28] = scale * (Du1x[16] * u0x[1] + Du1x[17] * (u0x[0] + 1.0));
    d2u1x[29] = 0.0;
    d2u1x[30] = scale * (Du1x[15] * u0x[3] + Du1x[17] * (u0x[4] + 1.0));
    d2u1x[31] = scale * (Du1x[16] * (u0x[4] + 1.0) + Du1x[17] * u0x[3]);
    d2u1x[32] = 0.0;
    d2u1x[33] = scale * (Du1x[15] * u0x[6] + Du1x[17] * u0x[7]);
    d2u1x[34] = scale * (Du1x[16] * u0x[7] + Du1x[17] * u0x[6]);
    d2u1x[35] = 0.0;

    d2u1x[36] = scale * (Du1x[21] * (u0x[0] + 1.0) + Du1x[23] * u0x[1]);
    d2u1x[37] = scale * (Du1x[22] * u0x[1] + Du1x[23] * (u0x[0] + 1.0));
    d2u1x[38] = 0.0;
    d2u1x[39] = scale * (Du1x[21] * u0x[3] + Du1x[23] * (u0x[4] + 1.0));
    d2u1x[40] = scale * (Du1x[22] * (u0x[4] + 1.0) + Du1x[23] * u0x[3]);
    d2u1x[41] = 0.0;
    d2u1x[42] = scale * (Du1x[21] * u0x[6] + Du1x[23] * u0x[7]);
    d2u1x[43] = scale * (Du1x[22] * u0x[7] + Du1x[23] * u0x[6]);
    d2u1x[44] = 0.0;

    d2u1x[45] = 0.0;
    d2u1x[46] = 0.0;
    d2u1x[47] = 0.0;
    d2u1x[48] = 0.0;
    d2u1x[49] = 0.0;
    d2u1x[50] = 0.0;
    d2u1x[51] = 0.0;
    d2u1x[52] = 0.0;
    d2u1x[53] = 0.0;

    d2u1x[54] = scale * (Du1x[27] * (u0x[0] + 1.0) + Du1x[29] * u0x[1]);
    d2u1x[55] = scale * (Du1x[28] * u0x[1] + Du1x[29] * (u0x[0] + 1.0));
    d2u1x[56] = 0.0;
    d2u1x[57] = scale * (Du1x[27] * u0x[3] + Du1x[29] * (u0x[4] + 1.0));
    d2u1x[58] = scale * (Du1x[28] * (u0x[4] + 1.0) + Du1x[29] * u0x[3]);
    d2u1x[59] = 0.0;
    d2u1x[60] = scale * (Du1x[27] * u0x[6] + Du1x[29] * u0x[7]);
    d2u1x[61] = scale * (Du1x[28] * u0x[7] + Du1x[29] * u0x[6]);
    d2u1x[62] = 0.0;

    d2u1x[63] = scale * (Du1x[33] * (u0x[0] + 1.0) + Du1x[35] * u0x[1]);
    d2u1x[64] = scale * (Du1x[34] * u0x[1] + Du1x[35] * (u0x[0] + 1.0));
    d2u1x[65] = 0.0;
    d2u1x[66] = scale * (Du1x[33] * u0x[3] + Du1x[35] * (u0x[4] + 1.0));
    d2u1x[67] = scale * (Du1x[34] * (u0x[4] + 1.0) + Du1x[35] * u0x[3]);
    d2u1x[68] = 0.0;
    d2u1x[69] = scale * (Du1x[33] * u0x[6] + Du1x[35] * u0x[7]);
    d2u1x[70] = scale * (Du1x[34] * u0x[7] + Du1x[35] * u0x[6]);
    d2u1x[71] = 0.0;

    d2u1x[72] = 0.0;
    d2u1x[73] = 0.0;
    d2u1x[74] = 0.0;
    d2u1x[75] = 0.0;
    d2u1x[76] = 0.0;
    d2u1x[77] = 0.0;
    d2u1x[78] = 0.0;
    d2u1x[79] = 0.0;
    d2u1x[80] = 0.0;

    d2u1xd[0] = scale * (Du1x[3] * u0xd[0] + Du1x[5] * u0xd[1]);
    d2u1xd[1] = scale * (Du1x[4] * u0xd[1] + Du1x[5] * u0xd[0]);
    d2u1xd[2] = 0.0;
    d2u1xd[3] = scale * (Du1x[3] * u0xd[3] + Du1x[5] * u0xd[4]);
    d2u1xd[4] = scale * (Du1x[4] * u0xd[4] + Du1x[5] * u0xd[3]);
    d2u1xd[5] = 0.0;
    d2u1xd[6] = scale * (Du1x[3] * u0xd[6] + Du1x[5] * u0xd[7]);
    d2u1xd[7] = scale * (Du1x[4] * u0xd[7] + Du1x[5] * u0xd[6]);
    d2u1xd[8] = 0.0;

    d2u1xd[9] = scale * (Du1x[11] * u0xd[1] + Du1x[9] * u0xd[0]);
    d2u1xd[10] = scale * (Du1x[10] * u0xd[1] + Du1x[11] * u0xd[0]);
    d2u1xd[11] = 0.0;
    d2u1xd[12] = scale * (Du1x[11] * u0xd[4] + Du1x[9] * u0xd[3]);
    d2u1xd[13] = scale * (Du1x[10] * u0xd[4] + Du1x[11] * u0xd[3]);
    d2u1xd[14] = 0.0;
    d2u1xd[15] = scale * (Du1x[11] * u0xd[7] + Du1x[9] * u0xd[6]);
    d2u1xd[16] = scale * (Du1x[10] * u0xd[7] + Du1x[11] * u0xd[6]);
    d2u1xd[17] = 0.0;

    d2u1xd[18] = 0.0;
    d2u1xd[19] = 0.0;
    d2u1xd[20] = 0.0;
    d2u1xd[21] = 0.0;
    d2u1xd[22] = 0.0;
    d2u1xd[23] = 0.0;
    d2u1xd[24] = 0.0;
    d2u1xd[25] = 0.0;
    d2u1xd[26] = 0.0;

    d2u1xd[27] = scale * (Du1x[15] * u0xd[0] + Du1x[17] * u0xd[1]);
    d2u1xd[28] = scale * (Du1x[16] * u0xd[1] + Du1x[17] * u0xd[0]);
    d2u1xd[29] = 0.0;
    d2u1xd[30] = scale * (Du1x[15] * u0xd[3] + Du1x[17] * u0xd[4]);
    d2u1xd[31] = scale * (Du1x[16] * u0xd[4] + Du1x[17] * u0xd[3]);
    d2u1xd[32] = 0.0;
    d2u1xd[33] = scale * (Du1x[15] * u0xd[6] + Du1x[17] * u0xd[7]);
    d2u1xd[34] = scale * (Du1x[16] * u0xd[7] + Du1x[17] * u0xd[6]);
    d2u1xd[35] = 0.0;

    d2u1xd[36] = scale * (Du1x[21] * u0xd[0] + Du1x[23] * u0xd[1]);
    d2u1xd[37] = scale * (Du1x[22] * u0xd[1] + Du1x[23] * u0xd[0]);
    d2u1xd[38] = 0.0;
    d2u1xd[39] = scale * (Du1x[21] * u0xd[3] + Du1x[23] * u0xd[4]);
    d2u1xd[40] = scale * (Du1x[22] * u0xd[4] + Du1x[23] * u0xd[3]);
    d2u1xd[41] = 0.0;
    d2u1xd[42] = scale * (Du1x[21] * u0xd[6] + Du1x[23] * u0xd[7]);
    d2u1xd[43] = scale * (Du1x[22] * u0xd[7] + Du1x[23] * u0xd[6]);
    d2u1xd[44] = 0.0;

    d2u1xd[45] = 0.0;
    d2u1xd[46] = 0.0;
    d2u1xd[47] = 0.0;
    d2u1xd[48] = 0.0;
    d2u1xd[49] = 0.0;
    d2u1xd[50] = 0.0;
    d2u1xd[51] = 0.0;
    d2u1xd[52] = 0.0;
    d2u1xd[53] = 0.0;

    d2u1xd[54] = scale * (Du1x[27] * u0xd[0] + Du1x[29] * u0xd[1]);
    d2u1xd[55] = scale * (Du1x[28] * u0xd[1] + Du1x[29] * u0xd[0]);
    d2u1xd[56] = 0.0;
    d2u1xd[57] = scale * (Du1x[27] * u0xd[3] + Du1x[29] * u0xd[4]);
    d2u1xd[58] = scale * (Du1x[28] * u0xd[4] + Du1x[29] * u0xd[3]);
    d2u1xd[59] = 0.0;
    d2u1xd[60] = scale * (Du1x[27] * u0xd[6] + Du1x[29] * u0xd[7]);
    d2u1xd[61] = scale * (Du1x[28] * u0xd[7] + Du1x[29] * u0xd[6]);
    d2u1xd[62] = 0.0;

    d2u1xd[63] = scale * (Du1x[33] * u0xd[0] + Du1x[35] * u0xd[1]);
    d2u1xd[64] = scale * (Du1x[34] * u0xd[1] + Du1x[35] * u0xd[0]);
    d2u1xd[65] = 0.0;
    d2u1xd[66] = scale * (Du1x[33] * u0xd[3] + Du1x[35] * u0xd[4]);
    d2u1xd[67] = scale * (Du1x[34] * u0xd[4] + Du1x[35] * u0xd[3]);
    d2u1xd[68] = 0.0;
    d2u1xd[69] = scale * (Du1x[33] * u0xd[6] + Du1x[35] * u0xd[7]);
    d2u1xd[70] = scale * (Du1x[34] * u0xd[7] + Du1x[35] * u0xd[6]);
    d2u1xd[71] = 0.0;

    d2u1xd[72] = 0.0;
    d2u1xd[73] = 0.0;
    d2u1xd[74] = 0.0;
    d2u1xd[75] = 0.0;
    d2u1xd[76] = 0.0;
    d2u1xd[77] = 0.0;
    d2u1xd[78] = 0.0;
    d2u1xd[79] = 0.0;
    d2u1xd[80] = 0.0;

    d2u1xd[0] += scale * (Du1xd[3] * (u0x[0] + 1.0) + Du1xd[5] * u0x[1]);
    d2u1xd[1] += scale * (Du1xd[4] * u0x[1] + Du1xd[5] * (u0x[0] + 1.0));
    d2u1xd[3] += scale * (Du1xd[3] * u0x[3] + Du1xd[5] * (u0x[4] + 1.0));
    d2u1xd[4] += scale * (Du1xd[4] * (u0x[4] + 1.0) + Du1xd[5] * u0x[3]);
    d2u1xd[6] += scale * (Du1xd[3] * u0x[6] + Du1xd[5] * u0x[7]);
    d2u1xd[7] += scale * (Du1xd[4] * u0x[7] + Du1xd[5] * u0x[6]);

    d2u1xd[9] += scale * (Du1xd[11] * u0x[1] + Du1xd[9] * (u0x[0] + 1.0));
    d2u1xd[10] += scale * (Du1xd[10] * u0x[1] + Du1xd[11] * (u0x[0] + 1.0));
    d2u1xd[12] += scale * (Du1xd[11] * (u0x[4] + 1.0) + Du1xd[9] * u0x[3]);
    d2u1xd[13] += scale * (Du1xd[10] * (u0x[4] + 1.0) + Du1xd[11] * u0x[3]);
    d2u1xd[15] += scale * (Du1xd[11] * u0x[7] + Du1xd[9] * u0x[6]);
    d2u1xd[16] += scale * (Du1xd[10] * u0x[7] + Du1xd[11] * u0x[6]);

    d2u1xd[27] += scale * (Du1xd[15] * (u0x[0] + 1.0) + Du1xd[17] * u0x[1]);
    d2u1xd[28] += scale * (Du1xd[16] * u0x[1] + Du1xd[17] * (u0x[0] + 1.0));
    d2u1xd[30] += scale * (Du1xd[15] * u0x[3] + Du1xd[17] * (u0x[4] + 1.0));
    d2u1xd[31] += scale * (Du1xd[16] * (u0x[4] + 1.0) + Du1xd[17] * u0x[3]);
    d2u1xd[33] += scale * (Du1xd[15] * u0x[6] + Du1xd[17] * u0x[7]);
    d2u1xd[34] += scale * (Du1xd[16] * u0x[7] + Du1xd[17] * u0x[6]);

    d2u1xd[36] += scale * (Du1xd[21] * (u0x[0] + 1.0) + Du1xd[23] * u0x[1]);
    d2u1xd[37] += scale * (Du1xd[22] * u0x[1] + Du1xd[23] * (u0x[0] + 1.0));
    d2u1xd[39] += scale * (Du1xd[21] * u0x[3] + Du1xd[23] * (u0x[4] + 1.0));
    d2u1xd[40] += scale * (Du1xd[22] * (u0x[4] + 1.0) + Du1xd[23] * u0x[3]);
    d2u1xd[42] += scale * (Du1xd[21] * u0x[6] + Du1xd[23] * u0x[7]);
    d2u1xd[43] += scale * (Du1xd[22] * u0x[7] + Du1xd[23] * u0x[6]);

    d2u1xd[54] += scale * (Du1xd[27] * (u0x[0] + 1.0) + Du1xd[29] * u0x[1]);
    d2u1xd[55] += scale * (Du1xd[28] * u0x[1] + Du1xd[29] * (u0x[0] + 1.0));
    d2u1xd[57] += scale * (Du1xd[27] * u0x[3] + Du1xd[29] * (u0x[4] + 1.0));
    d2u1xd[58] += scale * (Du1xd[28] * (u0x[4] + 1.0) + Du1xd[29] * u0x[3]);
    d2u1xd[60] += scale * (Du1xd[27] * u0x[6] + Du1xd[29] * u0x[7]);
    d2u1xd[61] += scale * (Du1xd[28] * u0x[7] + Du1xd[29] * u0x[6]);

    d2u1xd[63] += scale * (Du1xd[33] * (u0x[0] + 1.0) + Du1xd[35] * u0x[1]);
    d2u1xd[64] += scale * (Du1xd[34] * u0x[1] + Du1xd[35] * (u0x[0] + 1.0));
    d2u1xd[66] += scale * (Du1xd[33] * u0x[3] + Du1xd[35] * (u0x[4] + 1.0));
    d2u1xd[67] += scale * (Du1xd[34] * (u0x[4] + 1.0) + Du1xd[35] * u0x[3]);
    d2u1xd[69] += scale * (Du1xd[33] * u0x[6] + Du1xd[35] * u0x[7]);
    d2u1xd[70] += scale * (Du1xd[34] * u0x[7] + Du1xd[35] * u0x[6]);
  }
};

#endif  // TACS_SHELL_INPLANE_ELEMENT_MODEL_H
