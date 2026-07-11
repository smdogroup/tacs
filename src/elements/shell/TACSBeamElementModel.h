#ifndef TACS_BEAM_ELEMENT_MODEL_H
#define TACS_BEAM_ELEMENT_MODEL_H

#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"

class TACSBeamLinearModel {
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
  static void computeTyingStrain(const TacsScalar Xpts[],
                                 const TacsScalar fn1[], const TacsScalar fn2[],
                                 const TacsScalar vars[], const TacsScalar d1[],
                                 const TacsScalar d2[], TacsScalar ety[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[3], Xxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

      ety[index] = 0.0;
      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        basis::template interpFields<3, 3>(pt, fn1, n0);

        // Compute g12 = e1^{T}*G*e2
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        basis::template interpFields<3, 3>(pt, fn2, n0);

        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose(
      const TacsScalar Xpts[], const TacsScalar fn1[], const TacsScalar fn2[],
      const TacsScalar vars[], const TacsScalar d1[], const TacsScalar d2[],
      const TacsScalar dety[], TacsScalar res[], TacsScalar dd1[],
      TacsScalar dd2[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Xxi[3], dUxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar dd0[3], n0[3];
        basis::template interpFields<3, 3>(pt, fn1, n0);

        // Compute g12 = e1^{T}*G*e2
        dUxi[0] = 0.5 * dety[index] * n0[0];
        dUxi[1] = 0.5 * dety[index] * n0[1];
        dUxi[2] = 0.5 * dety[index] * n0[2];

        dd0[0] = 0.5 * dety[index] * Xxi[0];
        dd0[1] = 0.5 * dety[index] * Xxi[1];
        dd0[2] = 0.5 * dety[index] * Xxi[2];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd1);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar dd0[3], n0[3];
        basis::template interpFields<3, 3>(pt, fn2, n0);

        // Compute g13 = e1^{T}*G*e3
        dUxi[0] = 0.5 * dety[index] * n0[0];
        dUxi[1] = 0.5 * dety[index] * n0[1];
        dUxi[2] = 0.5 * dety[index] * n0[2];

        dd0[0] = 0.5 * dety[index] * Xxi[0];
        dd0[1] = 0.5 * dety[index] * Xxi[1];
        dd0[2] = 0.5 * dety[index] * Xxi[2];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd2);
      }

      if (res) {
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi,
                                                                       res);
      }
    }
  }

  /*
    Compute the directional derivative
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv(
      const TacsScalar Xpts[], const TacsScalar fn1[], const TacsScalar fn2[],
      const TacsScalar vars[], const TacsScalar d1[], const TacsScalar d2[],
      const TacsScalar psi[], const TacsScalar d1psi[],
      const TacsScalar d2psi[], TacsScalar ety[], TacsScalar etypsi[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[3], Xxi[3], Uxipsi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, Uxipsi);

      ety[index] = 0.0;
      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar d0[3], d0psi[3], n0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        basis::template interpFields<3, 3>(pt, d1psi, d0psi);
        basis::template interpFields<3, 3>(pt, fn1, n0);

        // Compute g12 = e1^{T}*G*e2
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
        etypsi[index] =
            0.5 * (Xxi[0] * d0psi[0] + Xxi[1] * d0psi[1] + Xxi[2] * d0psi[2] +
                   n0[0] * Uxipsi[0] + n0[1] * Uxipsi[1] + n0[2] * Uxipsi[2]);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], d0psi[3], n0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        basis::template interpFields<3, 3>(pt, d2psi, d0psi);
        basis::template interpFields<3, 3>(pt, fn2, n0);

        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
        etypsi[index] =
            0.5 * (Xxi[0] * d0psi[0] + Xxi[1] * d0psi[1] + Xxi[2] * d0psi[2] +
                   n0[0] * Uxipsi[0] + n0[1] * Uxipsi[1] + n0[2] * Uxipsi[2]);
      }
    }
  }

  /*
    Compute the derivative of the tying strain w.r.t. quantities that
    are computed from the node locations
  */
  template <int vars_per_node, class basis>
  static void addTyingStrainXptSens(
      const TacsScalar Xpts[], const TacsScalar fn1[], const TacsScalar fn2[],
      const TacsScalar vars[], const TacsScalar d1[], const TacsScalar d2[],
      const TacsScalar dety[], TacsScalar dfdXpts[], TacsScalar dfn1[],
      TacsScalar dfn2[], TacsScalar dd1[], TacsScalar dd2[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Xxi[3], Uxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

      TacsScalar dXxi[3];
      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar d0[3], dd0[3], dn0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        dXxi[0] = 0.5 * d0[0] * dety[index];
        dXxi[1] = 0.5 * d0[1] * dety[index];
        dXxi[2] = 0.5 * d0[2] * dety[index];

        dn0[0] = 0.5 * Uxi[0] * dety[index];
        dn0[1] = 0.5 * Uxi[1] * dety[index];
        dn0[2] = 0.5 * Uxi[2] * dety[index];

        dd0[0] = 0.5 * Xxi[0] * dety[index];
        dd0[1] = 0.5 * Xxi[1] * dety[index];
        dd0[2] = 0.5 * Xxi[2] * dety[index];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn1);
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd1);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], dd0[3], dn0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        dXxi[0] = 0.5 * d0[0] * dety[index];
        dXxi[1] = 0.5 * d0[1] * dety[index];
        dXxi[2] = 0.5 * d0[2] * dety[index];

        dn0[0] = 0.5 * Uxi[0] * dety[index];
        dn0[1] = 0.5 * Uxi[1] * dety[index];
        dn0[2] = 0.5 * Uxi[2] * dety[index];

        dd0[0] = 0.5 * Xxi[0] * dety[index];
        dd0[1] = 0.5 * Xxi[1] * dety[index];
        dd0[2] = 0.5 * Xxi[2] * dety[index];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn2);
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd2);
      }

      basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
    }
  }

  /*
    Compute the derivative of the tying strain w.r.t. quantities that
    are computed from the node locations
  */
  template <int vars_per_node, class basis>
  static void addTyingStrainDerivXptSens(
      const TacsScalar Xpts[], const TacsScalar fn1[], const TacsScalar fn2[],
      const TacsScalar vars[], const TacsScalar d1[], const TacsScalar d2[],
      const TacsScalar psi[], const TacsScalar d1psi[],
      const TacsScalar d2psi[], const TacsScalar dety[],
      const TacsScalar detypsi[], TacsScalar dfdXpts[], TacsScalar dfn1[],
      TacsScalar dfn2[], TacsScalar dd1[], TacsScalar dd2[],
      TacsScalar dd1psi[], TacsScalar dd2psi[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Xxi[3], Uxi[3], Uxipsi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, Uxipsi);

      TacsScalar dXxi[3];
      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar d0[3], dd0[3], dn0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        dXxi[0] = 0.5 * d0[0] * dety[index];
        dXxi[1] = 0.5 * d0[1] * dety[index];
        dXxi[2] = 0.5 * d0[2] * dety[index];

        basis::template interpFields<3, 3>(pt, d1psi, d0);
        dXxi[0] += 0.5 * d0[0] * detypsi[index];
        dXxi[1] += 0.5 * d0[1] * detypsi[index];
        dXxi[2] += 0.5 * d0[2] * detypsi[index];

        // Add contributions from the director fields
        dn0[0] = 0.5 * (Uxi[0] * dety[index] + Uxipsi[0] * detypsi[index]);
        dn0[1] = 0.5 * (Uxi[1] * dety[index] + Uxipsi[1] * detypsi[index]);
        dn0[2] = 0.5 * (Uxi[2] * dety[index] + Uxipsi[2] * detypsi[index]);
        basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn1);

        dd0[0] = 0.5 * Xxi[0] * dety[index];
        dd0[1] = 0.5 * Xxi[1] * dety[index];
        dd0[2] = 0.5 * Xxi[2] * dety[index];
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd1);

        dd0[0] = 0.5 * Xxi[0] * detypsi[index];
        dd0[1] = 0.5 * Xxi[1] * detypsi[index];
        dd0[2] = 0.5 * Xxi[2] * detypsi[index];
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd1psi);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], dd0[3], dn0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        dXxi[0] = 0.5 * d0[0] * dety[index];
        dXxi[1] = 0.5 * d0[1] * dety[index];
        dXxi[2] = 0.5 * d0[2] * dety[index];

        basis::template interpFields<3, 3>(pt, d2psi, d0);
        dXxi[0] += 0.5 * d0[0] * detypsi[index];
        dXxi[1] += 0.5 * d0[1] * detypsi[index];
        dXxi[2] += 0.5 * d0[2] * detypsi[index];

        // Add contributions from the director fields
        dn0[0] = 0.5 * (Uxi[0] * dety[index] + Uxipsi[0] * detypsi[index]);
        dn0[1] = 0.5 * (Uxi[1] * dety[index] + Uxipsi[1] * detypsi[index]);
        dn0[2] = 0.5 * (Uxi[2] * dety[index] + Uxipsi[2] * detypsi[index]);
        basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn2);

        dd0[0] = 0.5 * Xxi[0] * dety[index];
        dd0[1] = 0.5 * Xxi[1] * dety[index];
        dd0[2] = 0.5 * Xxi[2] * dety[index];
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd2);

        dd0[0] = 0.5 * Xxi[0] * detypsi[index];
        dd0[1] = 0.5 * Xxi[1] * detypsi[index];
        dd0[2] = 0.5 * Xxi[2] * detypsi[index];
        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd2psi);
      }

      basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
    }
  }

  static inline void evalStrain(const TacsScalar u0x[], const TacsScalar d1x[],
                                const TacsScalar d2x[], const TacsScalar e0ty[],
                                TacsScalar e[]) {
    // Axial strain
    e[0] = u0x[0];

    // Torsional component of the strain
    e[1] = 0.5 * (d1x[2] - d2x[1]);

    // Bending components of the strain
    e[2] = d1x[0];
    e[3] = d2x[0];

    // Add the tying shear strain components
    e[4] = e0ty[0];
    e[5] = e0ty[1];
  }

  static inline void evalStrainDeriv(
      const TacsScalar u0x[], const TacsScalar d1x[], const TacsScalar d2x[],
      const TacsScalar e0ty[], const TacsScalar u0xd[], const TacsScalar d1xd[],
      const TacsScalar d2xd[], const TacsScalar e0tyd[], TacsScalar e[],
      TacsScalar ed[]) {
    // Axial strain
    e[0] = u0x[0];

    // Torsional component of the strain
    e[1] = 0.5 * (d1x[2] - d2x[1]);

    // Bending components of the strain
    e[2] = d1x[0];
    e[3] = d2x[0];

    // Add the tying shear strain components
    e[4] = e0ty[0];
    e[5] = e0ty[1];

    // Axial strain
    ed[0] = u0xd[0];

    // Torsional component of the strain
    ed[1] = 0.5 * (d1xd[2] - d2xd[1]);

    // Bending components of the strain
    ed[2] = d1xd[0];
    ed[3] = d2xd[0];

    // Add the tying shear strain components
    ed[4] = e0tyd[0];
    ed[5] = e0tyd[1];
  }

  static inline void evalStrainSens(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar u0x[],
      const TacsScalar d1x[], const TacsScalar d2x[], const TacsScalar e0ty[],
      TacsScalar du0x[], TacsScalar dd1x[], TacsScalar dd2x[],
      TacsScalar de0ty[]) {
    du0x[0] = scale * dfde[0];
    du0x[1] = du0x[2] = 0.0;
    du0x[3] = du0x[4] = du0x[5] = 0.0;
    du0x[6] = du0x[7] = du0x[8] = 0.0;

    dd1x[0] = scale * dfde[2];
    dd1x[1] = 0.0;
    dd1x[2] = 0.5 * scale * dfde[1];

    dd2x[0] = scale * dfde[3];
    dd2x[1] = -0.5 * scale * dfde[1];
    dd2x[2] = 0.0;

    de0ty[0] = scale * dfde[4];
    de0ty[1] = scale * dfde[5];
  }

  /**
    Evaluate the Hessian of the strain energy-like quadratic form
    f = 0.5*s.e (s = Cs*e via TACSBeamConstitutive::computeStress) with
    respect to the four strain-map inputs (u0x, d1x, d2x, e0ty).

    Because e = L(u0x, d1x, d2x, e0ty) (see evalStrain, unchanged) is a
    fixed *linear* map with constant coefficients, this Hessian is a pure
    algebraic pullback of the packed-symmetric Cs[21] through L -- no
    differentiation is performed at runtime, and the result does not
    depend on the values of s/u0x/d1x/d2x/e0ty, only on Cs and scale (they
    remain in the signature to match the call-site contract and shell's
    analogous evalStrainHessian). Every Cs row/col index below was derived
    directly from TACSBeamConstitutive::computeStress's packed-symmetric
    index arithmetic (TACSBeamConstitutive.h:73-87), not transcribed from
    SPEC.md a second time; each of the seven blocks is verified
    independently against a central-difference-of-f ground truth (exact
    for this quadratic f) in the Task 1.2 op-level scratch check.
  */
  static void evalStrainHessian(const TacsScalar scale, const TacsScalar s[],
                                const TacsScalar Cs[], const TacsScalar u0x[],
                                const TacsScalar d1x[], const TacsScalar d2x[],
                                const TacsScalar e0ty[], TacsScalar d2u0x[],
                                TacsScalar d2d1x[], TacsScalar d2d2x[],
                                TacsScalar d2e0ty[], TacsScalar d2u0xd1x[],
                                TacsScalar d2u0xd2x[], TacsScalar d2d1xd2x[]) {
    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2d1x, 0, 9 * sizeof(TacsScalar));
    memset(d2d2x, 0, 9 * sizeof(TacsScalar));
    memset(d2e0ty, 0, 4 * sizeof(TacsScalar));
    memset(d2u0xd1x, 0, 27 * sizeof(TacsScalar));
    memset(d2u0xd2x, 0, 27 * sizeof(TacsScalar));
    memset(d2d1xd2x, 0, 9 * sizeof(TacsScalar));

    // e[0] = u0x[0] (axial) -- only u0x[0] enters e, so d2u0x is rank-1.
    d2u0x[0] = scale * Cs[0];

    // e[1] = 0.5*(d1x[2] - d2x[1]) (torsion), e[2] = d1x[0] (bend1).
    d2d1x[0] = scale * Cs[11];                  // bend1-bend1
    d2d1x[2] = d2d1x[6] = 0.5 * scale * Cs[7];  // bend1-torsion
    d2d1x[8] = 0.25 * scale * Cs[6];            // torsion-torsion

    // e[3] = d2x[0] (bend2); torsion cross term picks up the -0.5 sign from
    // evalStrainSens's dd2x[1] = -0.5*scale*dfde[1].
    d2d2x[0] = scale * Cs[15];                   // bend2-bend2
    d2d2x[1] = d2d2x[3] = -0.5 * scale * Cs[8];  // bend2-torsion
    d2d2x[4] = 0.25 * scale * Cs[6];             // torsion-torsion

    // d1x-d2x cross block (bend1-bend2 and the two torsion cross terms).
    d2d1xd2x[0] = scale * Cs[12];         // bend1 - bend2
    d2d1xd2x[1] = -0.5 * scale * Cs[7];   // bend1 - torsion(d2x)
    d2d1xd2x[6] = 0.5 * scale * Cs[8];    // torsion(d1x) - bend2
    d2d1xd2x[7] = -0.25 * scale * Cs[6];  // torsion-torsion cross

    // e[4] = e0ty[0], e[5] = e0ty[1]: direct 2x2 sub-block of Cs at {4,5}.
    d2e0ty[0] = scale * Cs[18];
    d2e0ty[1] = d2e0ty[2] = scale * Cs[19];
    d2e0ty[3] = scale * Cs[20];

    // Axial (u0x[0]) vs bend1/torsion(d1x).
    d2u0xd1x[0] = scale * Cs[2];
    d2u0xd1x[2] = 0.5 * scale * Cs[1];

    // Axial (u0x[0]) vs bend2/torsion(d2x); torsion cross term picks up the
    // same sign flip as d2d2x above.
    d2u0xd2x[0] = scale * Cs[3];
    d2u0xd2x[1] = -0.5 * scale * Cs[1];
  }

  /**
    Convert the tying-point-space tying-strain Hessian (d2ety, computed by
    basis::addInterpTyingStrainHessian from the quadrature-point-level
    strain Hessian) into mat[]'s vars-space entries, plus the director-
    space Hessian accumulators (d2d1/d2d2) and their cross terms with the
    translational field (d2d1u/d2d2u) that director::addDirectorJacobian
    consumes (SPEC.md sec 1.3.1, "Gap 1" -- found missing during Phase 2
    implementation; mirrors TACSShellElementModel::addComputeTyingStrainHessian,
    TACSShellElementModel.h:159-390, adapted to beam's two-director shape).

    Beam's version is simpler than shell's template: beam has exactly two
    tying fields (G12, G13), and -- confirmed directly from
    computeTyingStrain's formula above -- G12 depends only on d1/fn1 and
    G13 depends only on d2/fn2; the two fields never cross-couple in their
    FIRST derivative structure. Shell's arbitrary-tying-field-pair double
    loop therefore reduces, for beam, to a per-field loop (i1, i2 both G12,
    or both G13); cross-field (G12, G13) pairs are skipped. This relies on
    d2ety being genuinely zero for cross-field pairs, which in turn
    requires Cs's e0ty-e0ty off-diagonal entry (Cs[19], the G12-G13 shear
    coupling term feeding evalStrainHessian's d2e0ty[1]/d2e0ty[2]) to be
    zero -- confirmed directly from TACSIsoTubeBeamConstitutive::
    evalTangentStiffness (TACSIsoTubeBeamConstitutive.cpp:190-212, the only
    constitutive class test_beam_element.py exercises), which memsets Cs
    to zero and never sets index 19 (a circular/symmetric tube has no
    shear-direction coupling). A future constitutive model with nonzero
    Cs[19] would need this loop extended back to shell's full arbitrary-
    pair structure -- not needed today.

    Shell's more general version also takes external tying-strain-vs-
    other-strain cross-Hessian inputs (d2etyu/d2etyd, populated from a real
    material coupling via TacsShellAddTyingDispCoupling). Beam's
    evalStrainHessian above never produces any cross term between e0ty and
    (u0x, d1x, d2x) -- so beam's analogous d2etyu/d2etyd1/d2etyd2 inputs
    are always zero-filled placeholders at every call site in this
    feature, kept in the signature only for structural parity with
    shell's template (and possible future extensibility).

    @param alpha Unused (kept for signature parity with shell's template;
    every scale factor entering this closure is already baked into d2ety
    by the caller, via evalStrainHessian's own alpha*detXd scale argument)
    @param Xpts The element node locations
    @param fn1 The first reference normal direction at each node
    @param fn2 The second reference normal direction at each node
    @param vars The full variable vector (unused; kept for signature parity)
    @param d1 The first director field at each node (unused directly here;
    kept for signature parity, mirrors addComputeTyingStrainTranspose)
    @param d2 The second director field at each node (unused directly here)
    @param dety The first derivative of the tying strain (unused directly
    here; kept for signature parity)
    @param d2ety The NUM_TYING_POINTS x NUM_TYING_POINTS tying-point-space
    Hessian (already alpha-scaled)
    @param d2etyu Zero-filled placeholder, tying-strain vs u0xi cross
    Hessian (NUM_TYING_POINTS x dsize)
    @param d2etyd1 Zero-filled placeholder, tying-strain vs d1 cross
    Hessian (NUM_TYING_POINTS x dsize)
    @param d2etyd2 Zero-filled placeholder, tying-strain vs d2 cross
    Hessian (NUM_TYING_POINTS x dsize)
    @param mat The element Jacobian matrix (receives the (u0xi, u0xi) block
    directly)
    @param d2d1 The d1-director-space Hessian accumulator (dsize x dsize)
    @param d2d2 The d2-director-space Hessian accumulator (dsize x dsize)
    @param d2d1u The (d1, u0xi) cross-Hessian accumulator (dsize x dsize)
    @param d2d2u The (d2, u0xi) cross-Hessian accumulator (dsize x dsize)
  */
  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian(
      const TacsScalar alpha, const TacsScalar Xpts[], const TacsScalar fn1[],
      const TacsScalar fn2[], const TacsScalar vars[], const TacsScalar d1[],
      const TacsScalar d2[], const TacsScalar dety[], const TacsScalar d2ety[],
      const TacsScalar d2etyu[], const TacsScalar d2etyd1[],
      const TacsScalar d2etyd2[], TacsScalar mat[], TacsScalar d2d1[],
      TacsScalar d2d2[], TacsScalar d2d1u[], TacsScalar d2d2u[]) {
    const int dsize = 3 * basis::NUM_NODES;
    const int nvars = vars_per_node * basis::NUM_NODES;

    for (int i1 = 0; i1 < basis::NUM_TYING_POINTS; i1++) {
      const TacsBeamTyingStrainComponent f1 = basis::getTyingField(i1);
      const TacsScalar *n1field = (f1 == TACS_BEAM_G12_COMPONENT) ? fn1 : fn2;
      TacsScalar *d2d = (f1 == TACS_BEAM_G12_COMPONENT) ? d2d1 : d2d2;
      TacsScalar *d2du = (f1 == TACS_BEAM_G12_COMPONENT) ? d2d1u : d2d2u;

      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar Xxi1[3];
      basis::template interpFieldsGrad<3, 3>(pt1, Xpts, Xxi1);

      // Accumulate the i2-summed (d2ety[i1, :]-weighted) gradients for the
      // matching field only -- cross-field (G12, G13) pairs are skipped
      // (see the class-level comment above).
      TacsScalar du2[dsize], dd2[dsize];
      memset(du2, 0, dsize * sizeof(TacsScalar));
      memset(dd2, 0, dsize * sizeof(TacsScalar));

      for (int i2 = 0; i2 < basis::NUM_TYING_POINTS; i2++) {
        const TacsBeamTyingStrainComponent f2 = basis::getTyingField(i2);
        if (f2 != f1) {
          continue;
        }

        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        TacsScalar Xxi2[3], n02[3];
        basis::template interpFieldsGrad<3, 3>(pt2, Xpts, Xxi2);
        basis::template interpFields<3, 3>(pt2, n1field, n02);

        TacsScalar value = d2ety[basis::NUM_TYING_POINTS * i1 + i2];

        TacsScalar dUxi2[3], dd02[3];
        dUxi2[0] = 0.5 * value * n02[0];
        dUxi2[1] = 0.5 * value * n02[1];
        dUxi2[2] = 0.5 * value * n02[2];

        dd02[0] = 0.5 * value * Xxi2[0];
        dd02[1] = 0.5 * value * Xxi2[1];
        dd02[2] = 0.5 * value * Xxi2[2];

        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
        basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
      }

      // The (unweighted) gradient of ety[i1] itself, for the outer product
      // against the i2-summed vectors above.
      TacsScalar n01[3];
      basis::template interpFields<3, 3>(pt1, n1field, n01);

      TacsScalar dUxi1[3], dd01[3];
      dUxi1[0] = 0.5 * n01[0];
      dUxi1[1] = 0.5 * n01[1];
      dUxi1[2] = 0.5 * n01[2];

      dd01[0] = 0.5 * Xxi1[0];
      dd01[1] = 0.5 * Xxi1[1];
      dd01[2] = 0.5 * Xxi1[2];

      TacsScalar du1[dsize], dd1[dsize];
      memset(du1, 0, dsize * sizeof(TacsScalar));
      memset(dd1, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);
      basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);

      for (int i = 0; i < dsize; i++) {
        for (int j = 0; j < dsize; j++) {
          d2d[dsize * i + j] += dd1[i] * dd2[j];
          d2du[dsize * i + j] += dd1[i] * du2[j];
        }
      }

      for (int i = 0; i < dsize; i++) {
        int ii = vars_per_node * (i / 3) + (i % 3);
        for (int j = 0; j < dsize; j++) {
          int jj = vars_per_node * (j / 3) + (j % 3);
          mat[nvars * ii + jj] += du1[i] * du2[j];
        }
      }
    }
  }
};

/*

class TACSBeamNonlinearModel {
 public:

  template <int vars_per_node, class basis>
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn1[],
                                  const TacsScalar fn2[],
                                  const TacsScalar vars[],
                                  const TacsScalar d1[],
                                  const TacsScalar d2[],
                                  TacsScalar ety[] ){
    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Uxi[3], Xxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

      ety[index] = 0.0;
      if (field == TACS_BEAM_G12_COMPONENT){
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        basis::template interpFields<3, 3>(pt, fn1, n0);

        // Compute g12 = e1^{T}*G*e2
        ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[1]*d0[1] + Xxi[2]*d0[2] +
                          n0[0]*Uxi[0] + n0[1]*Uxi[1] + n0[2]*Uxi[2]);
      }
      else { // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        basis::template interpFields<3, 3>(pt, fn2, n0);

        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[1]*d0[1] + Xxi[2]*d0[2] +
                          n0[0]*Uxi[0] + n0[1]*Uxi[1] + n0[2]*Uxi[2]);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn1[],
                                              const TacsScalar fn2[],
                                              const TacsScalar vars[],
                                              const TacsScalar d1[],
                                              const TacsScalar d2[],
                                              const TacsScalar dety[],
                                              TacsScalar res[],
                                              TacsScalar dd1[],
                                              TacsScalar dd2[] ){
    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
      // Get the field index
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      // Interpolate the field value
      TacsScalar Xxi[3], dUxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

      if (field == TACS_BEAM_G12_COMPONENT){
        TacsScalar dd0[3], n0[3];
        basis::template interpFields<3, 3>(pt, fn1, n0);

        // Compute g12 = e1^{T}*G*e2
        dUxi[0] = 0.5*dety[index]*n0[0];
        dUxi[1] = 0.5*dety[index]*n0[1];
        dUxi[2] = 0.5*dety[index]*n0[2];

        dd0[0] = 0.5*dety[index]*Xxi[0];
        dd0[1] = 0.5*dety[index]*Xxi[1];
        dd0[2] = 0.5*dety[index]*Xxi[2];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd1);
      }
      else { // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar dd0[3], n0[3];
        basis::template interpFields<3, 3>(pt, fn2, n0);

        // Compute g13 = e1^{T}*G*e3
        dUxi[0] = 0.5*dety[index]*n0[0];
        dUxi[1] = 0.5*dety[index]*n0[1];
        dUxi[2] = 0.5*dety[index]*n0[2];

        dd0[0] = 0.5*dety[index]*Xxi[0];
        dd0[1] = 0.5*dety[index]*Xxi[1];
        dd0[2] = 0.5*dety[index]*Xxi[2];

        basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd2);
      }

      if (res){
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi,
res);
      }
    }
  }

  static inline void evalStrain( const TacsScalar u0x[],
                                 const TacsScalar d1x[],
                                 const TacsScalar d2x[],
                                 const TacsScalar e0ty[],
                                 TacsScalar e[] ){
    // Axial strain
    e[0] = u0x[0] + 0.5*(u0x[0]*u0x[0] + u0x[1]*u0x[1] + u0x[2]*u0x[2]);

    // Torsional strain
    e[1] = 0.5*(d1x[2] - d2x[1]) +
      (d1x[0]*u0x[2] + d1x[1]*u0x[5] + d1x[2]*u0x[8]) -
      (d2x[0]*u0x[1] + d2x[1]*u0x[4] + d2x[2]*u0x[7]));

    // Compute the bending components of the strain
    e[2] = d1x[0] + (u0x[0]*d1x[0] + u0x[1]*d1x[1] + u0x[2]*d1x[2]);
    e[3] = d2x[0] + (u0x[0]*d2x[0] + u0x[1]*d2x[1] + u0x[2]*d2x[2]);

    // Add the tying strain
    e[4] = e0ty[0];
    e[5] = e0ty[1];
  }

  static inline void evalStrainSens( const TacsScalar scale,
                                     const TacsScalar dfde[],
                                     const TacsScalar u0x[],
                                     const TacsScalar d1x[],
                                     const TacsScalar d2x[],
                                     TacsScalar du0x[],
                                     TacsScalar dd1x[],
                                     TacsScalar dd2x[] ){


  }
};
*/

template <int vars_per_node, class basis, class model>
int TacsTestBeamModelDerivatives(double dh = 1e-7, int test_print_level = 2,
                                 double test_fail_atol = 1e-5,
                                 double test_fail_rtol = 1e-5) {
  // Set the failure flag
  int fail = 0;

  // Set random values for the constitutive data and inputs
  TacsScalar Cs[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  TacsScalar u0x[9], d1x[3], d2x[3], e0ty[2];
  TacsScalar detXd;

  // Set random data
  TacsGenerateRandomArray(Cs,
                          TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES);
  TacsGenerateRandomArray(u0x, 9);
  TacsGenerateRandomArray(d1x, 3);
  TacsGenerateRandomArray(d2x, 3);
  TacsGenerateRandomArray(e0ty, 2);
  TacsGenerateRandomArray(&detXd, 1);

  // Compute the strain
  TacsScalar e[6];
  model::evalStrain(u0x, d1x, d2x, e0ty, e);

  // Compute the stress
  TacsScalar s[6];
  TACSBeamConstitutive::computeStress(Cs, e, s);

  // Compute the derivative of the product of the stress and strain
  // with respect to u0x, d1x, d2x and e0ty
  TacsScalar du0x[9], dd1x[3], dd2x[3], de0ty[2];
  model::evalStrainSens(detXd, s, u0x, d1x, d2x, e0ty, du0x, dd1x, dd2x, de0ty);

  TacsScalar f0 = 0.0;
  for (int j = 0; j < 6; j++) {
    f0 += 0.5 * detXd * e[j] * s[j];
  }

  // Compute against the derivatives for the strain
  TacsScalar fdu0x[9];
  for (int i = 0; i < 9; i++) {
    TacsScalar u0xt[9], et[6], st[6];
    memcpy(u0xt, u0x, 9 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u0xt[i] = u0x[i] + TacsScalar(0.0, dh);
#else
    u0xt[i] = u0x[i] + dh;
#endif  // TACS_USE_COMPLEX
    model::evalStrain(u0xt, d1x, d2x, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for (int j = 0; j < 6; j++) {
      f1 += 0.5 * detXd * et[j] * st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdu0x[i] = TacsImagPart(f1) / dh;
#else
    fdu0x[i] = (f1 - f0) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(du0x, fdu0x, 9, &max_err_index);
  double max_rel = TacsGetMaxRelError(du0x, fdu0x, 9, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "du0x", du0x, fdu0x, 9);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdd1x[3];
  for (int i = 0; i < 3; i++) {
    TacsScalar d1xt[3], et[6], st[6];
    memcpy(d1xt, d1x, 3 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    d1xt[i] = d1x[i] + TacsScalar(0.0, dh);
#else
    d1xt[i] = d1x[i] + dh;
#endif  // TACS_USE_COMPLEX
    model::evalStrain(u0x, d1xt, d2x, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for (int j = 0; j < 6; j++) {
      f1 += 0.5 * detXd * et[j] * st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdd1x[i] = TacsImagPart(f1) / dh;
#else
    fdd1x[i] = (f1 - f0) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dd1x, fdd1x, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd1x, fdd1x, 3, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. d1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dd1x", dd1x, fdd1x, 3);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdd2x[3];
  for (int i = 0; i < 3; i++) {
    TacsScalar d2xt[3], et[6], st[6];
    memcpy(d2xt, d2x, 3 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    d2xt[i] = d2x[i] + TacsScalar(0.0, dh);
#else
    d2xt[i] = d2x[i] + dh;
#endif  // TACS_USE_COMPLEX
    model::evalStrain(u0x, d1x, d2xt, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for (int j = 0; j < 6; j++) {
      f1 += 0.5 * detXd * et[j] * st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdd2x[i] = TacsImagPart(f1) / dh;
#else
    fdd2x[i] = (f1 - f0) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dd2x, fdd2x, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd2x, fdd2x, 3, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. d2x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "dd2x", dd2x, fdd2x, 3);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fde0ty[2];
  for (int i = 0; i < 2; i++) {
    TacsScalar e0tyt[2], et[6], st[6];
    memcpy(e0tyt, e0ty, 2 * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    e0tyt[i] = e0ty[i] + TacsScalar(0.0, dh);
#else
    e0tyt[i] = e0ty[i] + dh;
#endif  // TACS_USE_COMPLEX
    model::evalStrain(u0x, d1x, d2x, e0tyt, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for (int j = 0; j < 6; j++) {
      f1 += 0.5 * detXd * et[j] * st[j];
    }

#ifdef TACS_USE_COMPLEX
    fde0ty[i] = TacsImagPart(f1) / dh;
#else
    fde0ty[i] = (f1 - f0) / dh;
#endif  // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(de0ty, fde0ty, 2, &max_err_index);
  max_rel = TacsGetMaxRelError(de0ty, fde0ty, 2, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the derivative w.r.t. e0ty\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "de0ty", de0ty, fde0ty, 2);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  /*
  TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
  TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
  model::evalStrainHessian(detXd, s, Cs, u0x, u1x, e0ty,
                           d2u0x, d2u1x, d2u0xu1x,
                           d2e0ty, d2e0tyu0x, d2e0tyu1x);

  // Compute against the derivatives for the strain
  TacsScalar fd2u0x[81], fd2u0xu1x[81];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u0xt[9], et[9], st[9];
    memcpy(u0xt, u0x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u0xt[i] = u0x[i] + TacsScalar(0.0, dh);
#else
    u0xt[i] = u0x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0xt, u1x, e0ty, et);
    et[8] = 0.0;
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6];
    model::evalStrainSens(detXd, st, u0xt, u1x, du0xt, du1xt, de0tyt);

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2u0x[9*i + j] = TacsImagPart(du0xt[j])/dh;
      fd2u0xu1x[9*i + j] = TacsImagPart(du1xt[j])/dh;
#else
      fd2u0x[9*i + j] = (du0xt[j] - du0x[j])/dh;
      fd2u0xu1x[9*i + j] = (du1xt[j] - du1x[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2u0x, fd2u0x, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0x, fd2u0x, 81, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2u0x", d2u0x, fd2u0x, 81);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Compute the error
  max_err = TacsGetMaxError(d2u0xu1x, fd2u0xu1x, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0xu1x, fd2u0xu1x, 81, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. u0x and u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2u0xu1x", d2u0xu1x, fd2u0xu1x, 81);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Compute against the derivatives for the strain
  TacsScalar fd2u1x[81];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u1xt[9], et[9], st[9];
    memcpy(u1xt, u1x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u1xt[i] = u1x[i] + TacsScalar(0.0, dh);
#else
    u1xt[i] = u1x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1xt, e0ty,  et);
    et[8] = 0.0;
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6];
    model::evalStrainSens(detXd, st, u0x, u1xt, du0xt, du1xt, de0tyt);

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2u1x[9*i + j] = TacsImagPart(du1xt[j])/dh;
#else
      fd2u1x[9*i + j] = (du1xt[j] - du1x[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2u1x, fd2u1x, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u1x, fd2u1x, 81, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2u1x", d2u0x, fd2u0x, 81);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  TacsScalar fd2e0ty[36], fd2e0tyu0x[54], fd2e0tyu1x[54];
  for ( int i = 0; i < 6; i++ ){
    TacsScalar e0tyt[6], et[9], st[9];
    memcpy(e0tyt, e0ty, 6*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    e0tyt[i] = e0ty[i] + TacsScalar(0.0, dh);
#else
    e0tyt[i] = e0ty[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, e0tyt, et);
    et[8] = 0.0;
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6];
    model::evalStrainSens(detXd, st, u0x, u1x, du0xt, du1xt, de0tyt);

    for ( int j = 0; j < 6; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2e0ty[6*i + j] = TacsImagPart(de0tyt[j])/dh;
#else
      fd2e0ty[6*i + j] = (de0tyt[j] - de0ty[j])/dh;
#endif // TACS_USE_COMPLEX
    }

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2e0tyu0x[9*i + j] = TacsImagPart(du0xt[j])/dh;
      fd2e0tyu1x[9*i + j] = TacsImagPart(du1xt[j])/dh;
#else
      fd2e0tyu0x[9*i + j] = (du0xt[j] - du0x[j])/dh;
      fd2e0tyu1x[9*i + j] = (du1xt[j] - du1x[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2e0ty, fd2e0ty, 36, &max_err_index);
  max_rel = TacsGetMaxRelError(d2e0ty, fd2e0ty, 36, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. e0ty\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2e0ty", d2e0ty, fd2e0ty, 36);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Compute the error
  max_err = TacsGetMaxError(d2e0tyu0x, fd2e0tyu0x, 54, &max_err_index);
  max_rel = TacsGetMaxRelError(d2e0tyu0x, fd2e0tyu0x, 54, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. e0ty and u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2e0tyu0x", d2e0tyu0x, fd2e0tyu0x, 54);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

    // Compute the error
  max_err = TacsGetMaxError(d2e0tyu1x, fd2e0tyu1x, 54, &max_err_index);
  max_rel = TacsGetMaxRelError(d2e0tyu1x, fd2e0tyu1x, 54, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. e0ty and u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2e0tyu1x", d2e0tyu1x, fd2e0tyu1x, 54);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }
  */

  return fail;
}

#endif  // TACS_BEAM_ELEMENT_MODEL_H
