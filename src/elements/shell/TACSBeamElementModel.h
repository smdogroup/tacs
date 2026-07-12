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

class TACSBeamNonlinearModel {
 public:
  /*
    Tying strain (shear g12/g13) is left in its LINEAR form for this
    nonlinear model, byte-identical to TACSBeamLinearModel's own
    computeTyingStrain/addComputeTyingStrainTranspose/computeTyingStrainDeriv
    above -- confirmed directly from this (formerly commented-out) block's
    own computeTyingStrain body, which was already an exact copy of the
    linear model's before this task touched anything. Only evalStrain's
    axial/torsion/bending components pick up the geometric (u0x-quadratic
    and u0x-director-bilinear) nonlinearity; e0ty stays a linear function
    of (Xxi, Uxi, d, fn) throughout. computeTyingStrainDeriv is duplicated
    here (rather than shared) only so `model::computeTyingStrainDeriv` is
    callable identically for either model class from generic templated
    call sites (Task 5.2), per PLAN.md Task 5.1's explicit instruction.
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrain(const TacsScalar Xpts[],
                                 const TacsScalar fn1[], const TacsScalar fn2[],
                                 const TacsScalar vars[], const TacsScalar d1[],
                                 const TacsScalar d2[], TacsScalar ety[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      double pt[2];
      basis::getTyingPoint(index, pt);

      TacsScalar Uxi[3], Xxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

      ety[index] = 0.0;
      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d1, d0);
        basis::template interpFields<3, 3>(pt, fn1, n0);

        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
      } else {  // if (field == TACS_BEAM_G13_COMPONENT){
        TacsScalar d0[3], n0[3];
        basis::template interpFields<3, 3>(pt, d2, d0);
        basis::template interpFields<3, 3>(pt, fn2, n0);

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
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      double pt[2];
      basis::getTyingPoint(index, pt);

      TacsScalar Xxi[3], dUxi[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

      if (field == TACS_BEAM_G12_COMPONENT) {
        TacsScalar dd0[3], n0[3];
        basis::template interpFields<3, 3>(pt, fn1, n0);

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

  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv(
      const TacsScalar Xpts[], const TacsScalar fn1[], const TacsScalar fn2[],
      const TacsScalar vars[], const TacsScalar d1[], const TacsScalar d2[],
      const TacsScalar psi[], const TacsScalar d1psi[],
      const TacsScalar d2psi[], TacsScalar ety[], TacsScalar etypsi[]) {
    for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
      const TacsBeamTyingStrainComponent field = basis::getTyingField(index);

      double pt[2];
      basis::getTyingPoint(index, pt);

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

        ety[index] = 0.5 * (Xxi[0] * d0[0] + Xxi[1] * d0[1] + Xxi[2] * d0[2] +
                            n0[0] * Uxi[0] + n0[1] * Uxi[1] + n0[2] * Uxi[2]);
        etypsi[index] =
            0.5 * (Xxi[0] * d0psi[0] + Xxi[1] * d0psi[1] + Xxi[2] * d0psi[2] +
                   n0[0] * Uxipsi[0] + n0[1] * Uxipsi[1] + n0[2] * Uxipsi[2]);
      }
    }
  }

  /**
    Nonlinear strain measure (SPEC.md sec 2.4.4.a). u0x/d1x/d2x are indexed
    exactly as in TACSBeamLinearModel above: u0x[3*b+a] = d(u_a)/d(xi_b)
    with xi_0 the beam-axis parametric direction and xi_1/xi_2 the two
    transverse directions, so u0x[0..2] = d(u)/dxi_0 (the axial-direction
    gradient of the full displacement vector), u0x[3*1+j] = d(u)/dxi_1,
    u0x[3*2+j] = d(u)/dxi_2 -- confirmed by matching e[0]'s Green-Lagrange
    axial-strain form (du/dxi0 + 0.5*|du/dxi0|^2) against the linear
    model's e[0] = u0x[0] (its small-strain special case).

    e[0]: Green-Lagrange axial strain along the beam axis.
    e[1]: torsion strain, linear part unchanged from the linear model,
    plus the two Wagner-type bilinear cross terms d1x.(du/dxi2) and
    d2x.(du/dxi1) that couple twist to transverse-direction stretching.
    e[2]/e[3]: bending strains, linear part unchanged, plus the bilinear
    cross term (du/dxi0).d1x (resp. d2x) that couples curvature to axial
    stretching.
    e[4]/e[5]: tying shear strain, unchanged/linear (see class comment).
  */
  static inline void evalStrain(const TacsScalar u0x[], const TacsScalar d1x[],
                                const TacsScalar d2x[], const TacsScalar e0ty[],
                                TacsScalar e[]) {
    // Axial strain
    e[0] = u0x[0] + 0.5 * (u0x[0] * u0x[0] + u0x[1] * u0x[1] + u0x[2] * u0x[2]);

    // Torsional component of the strain
    e[1] = 0.5 * (d1x[2] - d2x[1]) +
           (d1x[0] * u0x[2] + d1x[1] * u0x[5] + d1x[2] * u0x[8]) -
           (d2x[0] * u0x[1] + d2x[1] * u0x[4] + d2x[2] * u0x[7]);

    // Bending components of the strain
    e[2] = d1x[0] + (u0x[0] * d1x[0] + u0x[1] * d1x[1] + u0x[2] * d1x[2]);
    e[3] = d2x[0] + (u0x[0] * d2x[0] + u0x[1] * d2x[1] + u0x[2] * d2x[2]);

    // Add the tying shear strain components
    e[4] = e0ty[0];
    e[5] = e0ty[1];
  }

  static inline void evalStrainDeriv(
      const TacsScalar u0x[], const TacsScalar d1x[], const TacsScalar d2x[],
      const TacsScalar e0ty[], const TacsScalar u0xd[], const TacsScalar d1xd[],
      const TacsScalar d2xd[], const TacsScalar e0tyd[], TacsScalar e[],
      TacsScalar ed[]) {
    evalStrain(u0x, d1x, d2x, e0ty, e);

    // Directional derivative of the axial strain
    ed[0] = u0xd[0] + (u0x[0] * u0xd[0] + u0x[1] * u0xd[1] + u0x[2] * u0xd[2]);

    // Directional derivative of the torsional strain
    ed[1] = 0.5 * (d1xd[2] - d2xd[1]) +
            (d1xd[0] * u0x[2] + d1x[0] * u0xd[2] + d1xd[1] * u0x[5] +
             d1x[1] * u0xd[5] + d1xd[2] * u0x[8] + d1x[2] * u0xd[8]) -
            (d2xd[0] * u0x[1] + d2x[0] * u0xd[1] + d2xd[1] * u0x[4] +
             d2x[1] * u0xd[4] + d2xd[2] * u0x[7] + d2x[2] * u0xd[7]);

    // Directional derivative of the bending strains
    ed[2] = d1xd[0] + (u0xd[0] * d1x[0] + u0x[0] * d1xd[0] + u0xd[1] * d1x[1] +
                       u0x[1] * d1xd[1] + u0xd[2] * d1x[2] + u0x[2] * d1xd[2]);
    ed[3] = d2xd[0] + (u0xd[0] * d2x[0] + u0x[0] * d2xd[0] + u0xd[1] * d2x[1] +
                       u0x[1] * d2xd[1] + u0xd[2] * d2x[2] + u0x[2] * d2xd[2]);

    // Add the tying shear strain components
    ed[4] = e0tyd[0];
    ed[5] = e0tyd[1];
  }

  /**
    Sensitivity of f = scale*dfde.e w.r.t. (u0x, d1x, d2x, e0ty). Widened
    to take/produce e0ty/de0ty (unlike this class's pre-Task-5.1 empty
    stub, which had neither in its signature) to match
    TACSBeamLinearModel::evalStrainSens's own shape -- confirmed necessary
    by this signature mismatch surfacing as a real compile error against
    TacsTestBeamModelDerivatives's existing (unchanged) call site, this
    task's RED step.
  */
  static inline void evalStrainSens(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar u0x[],
      const TacsScalar d1x[], const TacsScalar d2x[], const TacsScalar e0ty[],
      TacsScalar du0x[], TacsScalar dd1x[], TacsScalar dd2x[],
      TacsScalar de0ty[]) {
    du0x[0] = scale *
              (dfde[0] * (1.0 + u0x[0]) + dfde[2] * d1x[0] + dfde[3] * d2x[0]);
    du0x[1] = scale * (dfde[0] * u0x[1] - dfde[1] * d2x[0] + dfde[2] * d1x[1] +
                       dfde[3] * d2x[1]);
    du0x[2] = scale * (dfde[0] * u0x[2] + dfde[1] * d1x[0] + dfde[2] * d1x[2] +
                       dfde[3] * d2x[2]);
    du0x[3] = 0.0;
    du0x[4] = -scale * dfde[1] * d2x[1];
    du0x[5] = scale * dfde[1] * d1x[1];
    du0x[6] = 0.0;
    du0x[7] = -scale * dfde[1] * d2x[2];
    du0x[8] = scale * dfde[1] * d1x[2];

    dd1x[0] = scale * (dfde[1] * u0x[2] + dfde[2] * (1.0 + u0x[0]));
    dd1x[1] = scale * (dfde[1] * u0x[5] + dfde[2] * u0x[1]);
    dd1x[2] = scale * (dfde[1] * (0.5 + u0x[8]) + dfde[2] * u0x[2]);

    dd2x[0] = scale * (-dfde[1] * u0x[1] + dfde[3] * (1.0 + u0x[0]));
    dd2x[1] = scale * (-dfde[1] * (0.5 + u0x[4]) + dfde[3] * u0x[1]);
    dd2x[2] = scale * (-dfde[1] * u0x[7] + dfde[3] * u0x[2]);

    de0ty[0] = scale * dfde[4];
    de0ty[1] = scale * dfde[5];
  }

  /**
    Directional derivative (along u0xd/d1xd/d2xd/e0tyd, with dfde also
    perturbed by dfded) of evalStrainSens's output -- mechanical product-
    rule differentiation of every term in evalStrainSens above, verified
    against complex-step in the Task 5.1 scratch check (see PLAN.md's
    progress note for this task).
  */
  static inline void evalStrainSensDeriv(
      const TacsScalar scale, const TacsScalar dfde[], const TacsScalar u0x[],
      const TacsScalar d1x[], const TacsScalar d2x[], const TacsScalar e0ty[],
      const TacsScalar dfded[], const TacsScalar u0xd[],
      const TacsScalar d1xd[], const TacsScalar d2xd[],
      const TacsScalar e0tyd[], TacsScalar du0x[], TacsScalar dd1x[],
      TacsScalar dd2x[], TacsScalar de0ty[], TacsScalar du0xd[],
      TacsScalar dd1xd[], TacsScalar dd2xd[], TacsScalar de0tyd[]) {
    evalStrainSens(scale, dfde, u0x, d1x, d2x, e0ty, du0x, dd1x, dd2x, de0ty);

    du0xd[0] = scale * (dfded[0] * (1.0 + u0x[0]) + dfde[0] * u0xd[0] +
                        dfded[2] * d1x[0] + dfde[2] * d1xd[0] +
                        dfded[3] * d2x[0] + dfde[3] * d2xd[0]);
    du0xd[1] =
        scale * (dfded[0] * u0x[1] + dfde[0] * u0xd[1] - dfded[1] * d2x[0] -
                 dfde[1] * d2xd[0] + dfded[2] * d1x[1] + dfde[2] * d1xd[1] +
                 dfded[3] * d2x[1] + dfde[3] * d2xd[1]);
    du0xd[2] =
        scale * (dfded[0] * u0x[2] + dfde[0] * u0xd[2] + dfded[1] * d1x[0] +
                 dfde[1] * d1xd[0] + dfded[2] * d1x[2] + dfde[2] * d1xd[2] +
                 dfded[3] * d2x[2] + dfde[3] * d2xd[2]);
    du0xd[3] = 0.0;
    du0xd[4] = -scale * (dfded[1] * d2x[1] + dfde[1] * d2xd[1]);
    du0xd[5] = scale * (dfded[1] * d1x[1] + dfde[1] * d1xd[1]);
    du0xd[6] = 0.0;
    du0xd[7] = -scale * (dfded[1] * d2x[2] + dfde[1] * d2xd[2]);
    du0xd[8] = scale * (dfded[1] * d1x[2] + dfde[1] * d1xd[2]);

    dd1xd[0] = scale * (dfded[1] * u0x[2] + dfde[1] * u0xd[2] +
                        dfded[2] * (1.0 + u0x[0]) + dfde[2] * u0xd[0]);
    dd1xd[1] = scale * (dfded[1] * u0x[5] + dfde[1] * u0xd[5] +
                        dfded[2] * u0x[1] + dfde[2] * u0xd[1]);
    dd1xd[2] = scale * (dfded[1] * (0.5 + u0x[8]) + dfde[1] * u0xd[8] +
                        dfded[2] * u0x[2] + dfde[2] * u0xd[2]);

    dd2xd[0] = scale * (-dfded[1] * u0x[1] - dfde[1] * u0xd[1] +
                        dfded[3] * (1.0 + u0x[0]) + dfde[3] * u0xd[0]);
    dd2xd[1] = scale * (-dfded[1] * (0.5 + u0x[4]) - dfde[1] * u0xd[4] +
                        dfded[3] * u0x[1] + dfde[3] * u0xd[1]);
    dd2xd[2] = scale * (-dfded[1] * u0x[7] - dfde[1] * u0xd[7] +
                        dfded[3] * u0x[2] + dfde[3] * u0xd[2]);

    de0tyd[0] = scale * dfded[4];
    de0tyd[1] = scale * dfded[5];
  }

  /**
    Hessian of f = 0.5*s.e (s = Cs*e) w.r.t. (u0x, d1x, d2x, e0ty), for the
    genuinely nonlinear e[0..3] above. Unlike TACSBeamLinearModel's
    evalStrainHessian (a pure algebraic pullback of Cs, independent of the
    strain-state values, since e is linear there), this Hessian has two
    parts everywhere e is not linear:
      (1) the "material" part J^T*Cs*J, where J = de/dx is now
          state-dependent (e.g. de[0]/du0x[0] = 1+u0x[0]); and
      (2) a "geometric" part sum_i s[i]*(d2e[i]/dx2), nonzero only where
          e[i] is itself at least bilinear in x -- which, since e[0..3]
          above are each at most bilinear (e[0] is quadratic in u0x alone;
          e[1..3] are bilinear cross terms between u0x and d1x/d2x, with
          no self-quadratic d1x/d1x, d2x/d2x, or d1x/d2x term at all),
          means every one of these second derivatives is a CONSTANT
          matrix, independent of the current state (confirmed directly:
          differentiating evalStrain's e[0..3] twice never reproduces any
          of u0x/d1x/d2x). This is what makes evalStrainHessianDeriv below
          tractable without a third round of state-dependent algebra: only
          the *weights* (s vs. sd) change, never the constant geometric
          shape.

    Implemented via small (<=9-wide) local Jacobian-row loops rather than
    fully unrolled scalar assignments (a deviation from
    TACSBeamLinearModel::evalStrainHessian's fully-unrolled style,
    documented here per this feature's convention of flagging style
    deviations at their point of introduction): the nonlinear u0x block is
    9x9 (81 entries after the cross-blocks), an order of magnitude larger
    than the linear model's, and the loop form is verified directly
    against complex-step rather than transcribed by hand.
  */
  static void evalStrainHessian(const TacsScalar scale, const TacsScalar s[],
                                const TacsScalar Cs[], const TacsScalar u0x[],
                                const TacsScalar d1x[], const TacsScalar d2x[],
                                const TacsScalar e0ty[], TacsScalar d2u0x[],
                                TacsScalar d2d1x[], TacsScalar d2d2x[],
                                TacsScalar d2e0ty[], TacsScalar d2u0xd1x[],
                                TacsScalar d2u0xd2x[], TacsScalar d2d1xd2x[]) {
    // Jacobian rows (e0..e3 only -- e4/e5 are linear in e0ty alone and
    // handled separately below) of de/du0x, de/dd1x, de/dd2x.
    TacsScalar Ju0x[4][9], Jd1x[4][3], Jd2x[4][3];
    memset(Ju0x, 0, sizeof(Ju0x));
    memset(Jd1x, 0, sizeof(Jd1x));
    memset(Jd2x, 0, sizeof(Jd2x));

    Ju0x[0][0] = 1.0 + u0x[0];
    Ju0x[0][1] = u0x[1];
    Ju0x[0][2] = u0x[2];

    Ju0x[1][1] = -d2x[0];
    Ju0x[1][2] = d1x[0];
    Ju0x[1][4] = -d2x[1];
    Ju0x[1][5] = d1x[1];
    Ju0x[1][7] = -d2x[2];
    Ju0x[1][8] = d1x[2];

    Ju0x[2][0] = d1x[0];
    Ju0x[2][1] = d1x[1];
    Ju0x[2][2] = d1x[2];

    Ju0x[3][0] = d2x[0];
    Ju0x[3][1] = d2x[1];
    Ju0x[3][2] = d2x[2];

    Jd1x[1][0] = u0x[2];
    Jd1x[1][1] = u0x[5];
    Jd1x[1][2] = 0.5 + u0x[8];
    Jd1x[2][0] = 1.0 + u0x[0];
    Jd1x[2][1] = u0x[1];
    Jd1x[2][2] = u0x[2];

    Jd2x[1][0] = -u0x[1];
    Jd2x[1][1] = -0.5 - u0x[4];
    Jd2x[1][2] = -u0x[7];
    Jd2x[3][0] = 1.0 + u0x[0];
    Jd2x[3][1] = u0x[1];
    Jd2x[3][2] = u0x[2];

    // 4x4 sub-block of Cs over the strain indices {0 (axial), 1 (torsion),
    // 2 (bend1), 3 (bend2)} -- packed-symmetric indices per
    // TACSBeamConstitutive::computeStress.
    TacsScalar C4[4][4];
    C4[0][0] = Cs[0];
    C4[0][1] = C4[1][0] = Cs[1];
    C4[0][2] = C4[2][0] = Cs[2];
    C4[0][3] = C4[3][0] = Cs[3];
    C4[1][1] = Cs[6];
    C4[1][2] = C4[2][1] = Cs[7];
    C4[1][3] = C4[3][1] = Cs[8];
    C4[2][2] = Cs[11];
    C4[2][3] = C4[3][2] = Cs[12];
    C4[3][3] = Cs[15];

    memset(d2u0x, 0, 81 * sizeof(TacsScalar));
    memset(d2d1x, 0, 9 * sizeof(TacsScalar));
    memset(d2d2x, 0, 9 * sizeof(TacsScalar));
    memset(d2e0ty, 0, 4 * sizeof(TacsScalar));
    memset(d2u0xd1x, 0, 27 * sizeof(TacsScalar));
    memset(d2u0xd2x, 0, 27 * sizeof(TacsScalar));
    memset(d2d1xd2x, 0, 9 * sizeof(TacsScalar));

    for (int a = 0; a < 9; a++) {
      for (int b = 0; b < 9; b++) {
        TacsScalar val = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            val += Ju0x[i][a] * C4[i][j] * Ju0x[j][b];
          }
        }
        d2u0x[9 * a + b] = scale * val;
      }
    }
    // Geometric term: only e[0]'s pure-quadratic u0x self-term contributes.
    d2u0x[9 * 0 + 0] += scale * s[0];
    d2u0x[9 * 1 + 1] += scale * s[0];
    d2u0x[9 * 2 + 2] += scale * s[0];

    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        TacsScalar vald1 = 0.0, vald2 = 0.0, valcross = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            vald1 += Jd1x[i][a] * C4[i][j] * Jd1x[j][b];
            vald2 += Jd2x[i][a] * C4[i][j] * Jd2x[j][b];
            valcross += Jd1x[i][a] * C4[i][j] * Jd2x[j][b];
          }
        }
        d2d1x[3 * a + b] = scale * vald1;
        d2d2x[3 * a + b] = scale * vald2;
        d2d1xd2x[3 * a + b] = scale * valcross;
      }
    }
    // No geometric contribution to d2d1x/d2d2x/d2d1xd2x: none of e[0..3]
    // is quadratic in d1x/d2x alone or bilinear between them.

    for (int a = 0; a < 9; a++) {
      for (int b = 0; b < 3; b++) {
        TacsScalar val1 = 0.0, val2 = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            val1 += Ju0x[i][a] * C4[i][j] * Jd1x[j][b];
            val2 += Ju0x[i][a] * C4[i][j] * Jd2x[j][b];
          }
        }
        d2u0xd1x[3 * a + b] = scale * val1;
        d2u0xd2x[3 * a + b] = scale * val2;
      }
    }
    // Geometric term for (u0x, d1x): e[1]'s d1x0-u0x2/d1x1-u0x5/d1x2-u0x8
    // bilinear terms, plus e[2]'s u0x.d1x bilinear term.
    d2u0xd1x[3 * 0 + 0] += scale * s[2];
    d2u0xd1x[3 * 1 + 1] += scale * s[2];
    d2u0xd1x[3 * 2 + 2] += scale * s[2];
    d2u0xd1x[3 * 2 + 0] += scale * s[1];
    d2u0xd1x[3 * 5 + 1] += scale * s[1];
    d2u0xd1x[3 * 8 + 2] += scale * s[1];

    // Geometric term for (u0x, d2x): e[1]'s -d2x0-u0x1/-d2x1-u0x4/-d2x2-u0x7
    // bilinear terms, plus e[3]'s u0x.d2x bilinear term.
    d2u0xd2x[3 * 0 + 0] += scale * s[3];
    d2u0xd2x[3 * 1 + 1] += scale * s[3];
    d2u0xd2x[3 * 2 + 2] += scale * s[3];
    d2u0xd2x[3 * 1 + 0] -= scale * s[1];
    d2u0xd2x[3 * 4 + 1] -= scale * s[1];
    d2u0xd2x[3 * 7 + 2] -= scale * s[1];

    // e[4] = e0ty[0], e[5] = e0ty[1]: direct 2x2 sub-block of Cs at {4,5},
    // exactly as in the linear model (e0ty stays linear here too).
    d2e0ty[0] = scale * Cs[18];
    d2e0ty[1] = d2e0ty[2] = scale * Cs[19];
    d2e0ty[3] = scale * Cs[20];
  }

  /**
    Directional derivative (along the "path" direction u0xd/d1xd/d2xd,
    with the stress itself also perturbed to sd = Cs*ed) of
    evalStrainHessian's output. Per the class comment on evalStrainHessian,
    every geometric-term shape above is a CONSTANT matrix -- so its
    directional derivative needs only s -> sd, no new shape derivation.
    The material term J^T*Cs*J does need differentiating, via
    d(J^T*Cs*J)/dh = Jd^T*Cs*J + J^T*Cs*Jd, where Jd (the Jacobian rows'
    own directional derivatives) is exactly the same row structure as J
    with each state variable replaced by its "d"-direction counterpart
    (mechanical, since every J entry above is at most linear in the
    state) -- verified against complex-step in the Task 5.1/5.2 scratch
    checks, not hand-verified term-by-term a second time.
  */
  static void evalStrainHessianDeriv(
      const TacsScalar scale, const TacsScalar s[], const TacsScalar Cs[],
      const TacsScalar u0x[], const TacsScalar d1x[], const TacsScalar d2x[],
      const TacsScalar e0ty[], const TacsScalar sd[], const TacsScalar u0xd[],
      const TacsScalar d1xd[], const TacsScalar d2xd[],
      const TacsScalar e0tyd[], TacsScalar d2u0x[], TacsScalar d2d1x[],
      TacsScalar d2d2x[], TacsScalar d2e0ty[], TacsScalar d2u0xd1x[],
      TacsScalar d2u0xd2x[], TacsScalar d2d1xd2x[], TacsScalar d2u0xd[],
      TacsScalar d2d1xd[], TacsScalar d2d2xd[], TacsScalar d2e0tyd[],
      TacsScalar d2u0xd1xd[], TacsScalar d2u0xd2xd[], TacsScalar d2d1xd2xd[]) {
    // Base (unperturbed) Hessian blocks -- identical to evalStrainHessian.
    evalStrainHessian(scale, s, Cs, u0x, d1x, d2x, e0ty, d2u0x, d2d1x, d2d2x,
                      d2e0ty, d2u0xd1x, d2u0xd2x, d2d1xd2x);

    TacsScalar Ju0x[4][9], Jd1x[4][3], Jd2x[4][3];
    TacsScalar Ju0xd[4][9], Jd1xd[4][3], Jd2xd[4][3];
    memset(Ju0x, 0, sizeof(Ju0x));
    memset(Jd1x, 0, sizeof(Jd1x));
    memset(Jd2x, 0, sizeof(Jd2x));
    memset(Ju0xd, 0, sizeof(Ju0xd));
    memset(Jd1xd, 0, sizeof(Jd1xd));
    memset(Jd2xd, 0, sizeof(Jd2xd));

    Ju0x[0][0] = 1.0 + u0x[0];
    Ju0x[0][1] = u0x[1];
    Ju0x[0][2] = u0x[2];
    Ju0x[1][1] = -d2x[0];
    Ju0x[1][2] = d1x[0];
    Ju0x[1][4] = -d2x[1];
    Ju0x[1][5] = d1x[1];
    Ju0x[1][7] = -d2x[2];
    Ju0x[1][8] = d1x[2];
    Ju0x[2][0] = d1x[0];
    Ju0x[2][1] = d1x[1];
    Ju0x[2][2] = d1x[2];
    Ju0x[3][0] = d2x[0];
    Ju0x[3][1] = d2x[1];
    Ju0x[3][2] = d2x[2];

    Jd1x[1][0] = u0x[2];
    Jd1x[1][1] = u0x[5];
    Jd1x[1][2] = 0.5 + u0x[8];
    Jd1x[2][0] = 1.0 + u0x[0];
    Jd1x[2][1] = u0x[1];
    Jd1x[2][2] = u0x[2];

    Jd2x[1][0] = -u0x[1];
    Jd2x[1][1] = -0.5 - u0x[4];
    Jd2x[1][2] = -u0x[7];
    Jd2x[3][0] = 1.0 + u0x[0];
    Jd2x[3][1] = u0x[1];
    Jd2x[3][2] = u0x[2];

    // Directional derivatives of the same Jacobian rows: each entry above
    // is at most linear in the state, so its own directional derivative is
    // the same formula with the state replaced by its "d" counterpart.
    Ju0xd[0][0] = u0xd[0];
    Ju0xd[0][1] = u0xd[1];
    Ju0xd[0][2] = u0xd[2];
    Ju0xd[1][1] = -d2xd[0];
    Ju0xd[1][2] = d1xd[0];
    Ju0xd[1][4] = -d2xd[1];
    Ju0xd[1][5] = d1xd[1];
    Ju0xd[1][7] = -d2xd[2];
    Ju0xd[1][8] = d1xd[2];
    Ju0xd[2][0] = d1xd[0];
    Ju0xd[2][1] = d1xd[1];
    Ju0xd[2][2] = d1xd[2];
    Ju0xd[3][0] = d2xd[0];
    Ju0xd[3][1] = d2xd[1];
    Ju0xd[3][2] = d2xd[2];

    Jd1xd[1][0] = u0xd[2];
    Jd1xd[1][1] = u0xd[5];
    Jd1xd[1][2] = u0xd[8];
    Jd1xd[2][0] = u0xd[0];
    Jd1xd[2][1] = u0xd[1];
    Jd1xd[2][2] = u0xd[2];

    Jd2xd[1][0] = -u0xd[1];
    Jd2xd[1][1] = -u0xd[4];
    Jd2xd[1][2] = -u0xd[7];
    Jd2xd[3][0] = u0xd[0];
    Jd2xd[3][1] = u0xd[1];
    Jd2xd[3][2] = u0xd[2];

    TacsScalar C4[4][4];
    C4[0][0] = Cs[0];
    C4[0][1] = C4[1][0] = Cs[1];
    C4[0][2] = C4[2][0] = Cs[2];
    C4[0][3] = C4[3][0] = Cs[3];
    C4[1][1] = Cs[6];
    C4[1][2] = C4[2][1] = Cs[7];
    C4[1][3] = C4[3][1] = Cs[8];
    C4[2][2] = Cs[11];
    C4[2][3] = C4[3][2] = Cs[12];
    C4[3][3] = Cs[15];

    memset(d2u0xd, 0, 81 * sizeof(TacsScalar));
    memset(d2d1xd, 0, 9 * sizeof(TacsScalar));
    memset(d2d2xd, 0, 9 * sizeof(TacsScalar));
    memset(d2e0tyd, 0, 4 * sizeof(TacsScalar));
    memset(d2u0xd1xd, 0, 27 * sizeof(TacsScalar));
    memset(d2u0xd2xd, 0, 27 * sizeof(TacsScalar));
    memset(d2d1xd2xd, 0, 9 * sizeof(TacsScalar));

    for (int a = 0; a < 9; a++) {
      for (int b = 0; b < 9; b++) {
        TacsScalar val = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            val += Ju0xd[i][a] * C4[i][j] * Ju0x[j][b] +
                   Ju0x[i][a] * C4[i][j] * Ju0xd[j][b];
          }
        }
        d2u0xd[9 * a + b] = scale * val;
      }
    }
    d2u0xd[9 * 0 + 0] += scale * sd[0];
    d2u0xd[9 * 1 + 1] += scale * sd[0];
    d2u0xd[9 * 2 + 2] += scale * sd[0];

    for (int a = 0; a < 3; a++) {
      for (int b = 0; b < 3; b++) {
        TacsScalar vald1 = 0.0, vald2 = 0.0, valcross = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            vald1 += Jd1xd[i][a] * C4[i][j] * Jd1x[j][b] +
                     Jd1x[i][a] * C4[i][j] * Jd1xd[j][b];
            vald2 += Jd2xd[i][a] * C4[i][j] * Jd2x[j][b] +
                     Jd2x[i][a] * C4[i][j] * Jd2xd[j][b];
            valcross += Jd1xd[i][a] * C4[i][j] * Jd2x[j][b] +
                        Jd1x[i][a] * C4[i][j] * Jd2xd[j][b];
          }
        }
        d2d1xd[3 * a + b] = scale * vald1;
        d2d2xd[3 * a + b] = scale * vald2;
        d2d1xd2xd[3 * a + b] = scale * valcross;
      }
    }

    for (int a = 0; a < 9; a++) {
      for (int b = 0; b < 3; b++) {
        TacsScalar val1 = 0.0, val2 = 0.0;
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
            val1 += Ju0xd[i][a] * C4[i][j] * Jd1x[j][b] +
                    Ju0x[i][a] * C4[i][j] * Jd1xd[j][b];
            val2 += Ju0xd[i][a] * C4[i][j] * Jd2x[j][b] +
                    Ju0x[i][a] * C4[i][j] * Jd2xd[j][b];
          }
        }
        d2u0xd1xd[3 * a + b] = scale * val1;
        d2u0xd2xd[3 * a + b] = scale * val2;
      }
    }
    d2u0xd1xd[3 * 0 + 0] += scale * sd[2];
    d2u0xd1xd[3 * 1 + 1] += scale * sd[2];
    d2u0xd1xd[3 * 2 + 2] += scale * sd[2];
    d2u0xd1xd[3 * 2 + 0] += scale * sd[1];
    d2u0xd1xd[3 * 5 + 1] += scale * sd[1];
    d2u0xd1xd[3 * 8 + 2] += scale * sd[1];

    d2u0xd2xd[3 * 0 + 0] += scale * sd[3];
    d2u0xd2xd[3 * 1 + 1] += scale * sd[3];
    d2u0xd2xd[3 * 2 + 2] += scale * sd[3];
    d2u0xd2xd[3 * 1 + 0] -= scale * sd[1];
    d2u0xd2xd[3 * 4 + 1] -= scale * sd[1];
    d2u0xd2xd[3 * 7 + 2] -= scale * sd[1];

    // d2e0ty is a Cs-only constant (e0ty stays linear) -- its directional
    // derivative is exactly zero.
    d2e0tyd[0] = d2e0tyd[1] = d2e0tyd[2] = d2e0tyd[3] = 0.0;
  }
};

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

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

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

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

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

  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the analytic Hessian blocks (ported from
  // TacsTestShellModelDerivatives's structure, ~TACSShellElementModel.h:4110
  // onward, adapted to beam's smaller 4-nonlinear-strain kinematics -- see
  // PLAN.md Task 5.1).
  TacsScalar d2u0x[81], d2d1x[9], d2d2x[9], d2e0ty[4], d2u0xd1x[27],
      d2u0xd2x[27], d2d1xd2x[9];
  model::evalStrainHessian(detXd, s, Cs, u0x, d1x, d2x, e0ty, d2u0x, d2d1x,
                           d2d2x, d2e0ty, d2u0xd1x, d2u0xd2x, d2d1xd2x);

  // FD of du0x/dd1x/dd2x/de0ty w.r.t. u0x gives the u0x row of every block
  TacsScalar fd2u0x[81], fd2u0xd1x[27], fd2u0xd2x[27];
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

    TacsScalar du0xt[9], dd1xt[3], dd2xt[3], de0tyt[2];
    model::evalStrainSens(detXd, st, u0xt, d1x, d2x, e0ty, du0xt, dd1xt, dd2xt,
                          de0tyt);

    for (int j = 0; j < 9; j++) {
#ifdef TACS_USE_COMPLEX
      fd2u0x[9 * i + j] = TacsImagPart(du0xt[j]) / dh;
#else
      fd2u0x[9 * i + j] = (du0xt[j] - du0x[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
    for (int j = 0; j < 3; j++) {
#ifdef TACS_USE_COMPLEX
      fd2u0xd1x[3 * i + j] = TacsImagPart(dd1xt[j]) / dh;
      fd2u0xd2x[3 * i + j] = TacsImagPart(dd2xt[j]) / dh;
#else
      fd2u0xd1x[3 * i + j] = (dd1xt[j] - dd1x[j]) / dh;
      fd2u0xd2x[3 * i + j] = (dd2xt[j] - dd2x[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  max_err = TacsGetMaxError(d2u0x, fd2u0x, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0x, fd2u0x, 81, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2u0x", d2u0x, fd2u0x, 81);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  max_err = TacsGetMaxError(d2u0xd1x, fd2u0xd1x, 27, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0xd1x, fd2u0xd1x, 27, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. u0x and d1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2u0xd1x", d2u0xd1x, fd2u0xd1x, 27);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  max_err = TacsGetMaxError(d2u0xd2x, fd2u0xd2x, 27, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0xd2x, fd2u0xd2x, 27, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. u0x and d2x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2u0xd2x", d2u0xd2x, fd2u0xd2x, 27);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // FD of dd1x/dd2x w.r.t. d1x gives d2d1x and (transposed) d2u0xd1x/d2d1xd2x
  TacsScalar fd2d1x[9], fd2d1xd2x[9];
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

    TacsScalar du0xt[9], dd1xt[3], dd2xt[3], de0tyt[2];
    model::evalStrainSens(detXd, st, u0x, d1xt, d2x, e0ty, du0xt, dd1xt, dd2xt,
                          de0tyt);

    for (int j = 0; j < 3; j++) {
#ifdef TACS_USE_COMPLEX
      fd2d1x[3 * i + j] = TacsImagPart(dd1xt[j]) / dh;
      fd2d1xd2x[3 * i + j] = TacsImagPart(dd2xt[j]) / dh;
#else
      fd2d1x[3 * i + j] = (dd1xt[j] - dd1x[j]) / dh;
      fd2d1xd2x[3 * i + j] = (dd2xt[j] - dd2x[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  max_err = TacsGetMaxError(d2d1x, fd2d1x, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d1x, fd2d1x, 9, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. d1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2d1x", d2d1x, fd2d1x, 9);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  max_err = TacsGetMaxError(d2d1xd2x, fd2d1xd2x, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d1xd2x, fd2d1xd2x, 9, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. d1x and d2x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2d1xd2x", d2d1xd2x, fd2d1xd2x, 9);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // FD of dd2x w.r.t. d2x gives d2d2x
  TacsScalar fd2d2x[9];
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

    TacsScalar du0xt[9], dd1xt[3], dd2xt[3], de0tyt[2];
    model::evalStrainSens(detXd, st, u0x, d1x, d2xt, e0ty, du0xt, dd1xt, dd2xt,
                          de0tyt);

    for (int j = 0; j < 3; j++) {
#ifdef TACS_USE_COMPLEX
      fd2d2x[3 * i + j] = TacsImagPart(dd2xt[j]) / dh;
#else
      fd2d2x[3 * i + j] = (dd2xt[j] - dd2x[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  max_err = TacsGetMaxError(d2d2x, fd2d2x, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d2x, fd2d2x, 9, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. d2x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2d2x", d2d2x, fd2d2x, 9);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // FD of de0ty w.r.t. e0ty gives d2e0ty
  TacsScalar fd2e0ty[4];
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

    TacsScalar du0xt[9], dd1xt[3], dd2xt[3], de0tyt[2];
    model::evalStrainSens(detXd, st, u0x, d1x, d2x, e0tyt, du0xt, dd1xt, dd2xt,
                          de0tyt);

    for (int j = 0; j < 2; j++) {
#ifdef TACS_USE_COMPLEX
      fd2e0ty[2 * i + j] = TacsImagPart(de0tyt[j]) / dh;
#else
      fd2e0ty[2 * i + j] = (de0tyt[j] - de0ty[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  max_err = TacsGetMaxError(d2e0ty, fd2e0ty, 4, &max_err_index);
  max_rel = TacsGetMaxRelError(d2e0ty, fd2e0ty, 4, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. e0ty\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2e0ty", d2e0ty, fd2e0ty, 4);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}
/*
  FD/CS checks of TACSBeamNonlinearModel's "...Deriv"-suffixed family
  (evalStrainDeriv, evalStrainSensDeriv, evalStrainHessianDeriv), added per
  Phase 5's review feedback: these are exactly the methods
  TACSBeamElement::getMatType's TACS_GEOMETRIC_STIFFNESS_MATRIX branch
  depends on, and (unlike evalStrain/evalStrainSens/evalStrainHessian,
  covered by TacsTestBeamModelDerivatives above) had no committed,
  reproducible regression coverage before this function -- only an
  uncommitted scratch complex-step driver, per this task's own PLAN.md
  note.

  This is a SEPARATE function from TacsTestBeamModelDerivatives (rather
  than folded into it) because that function is templated on `model` and
  also instantiated with TACSBeamLinearModel, which does not have (and, by
  this feature's own design -- mirroring TACSShellLinearModel's identical
  omission -- is not meant to have) any of these three methods; hardcoding
  this function to TACSBeamNonlinearModel avoids breaking that other
  instantiation.
*/
int TacsTestBeamNonlinearModelDerivFamily(double dh = 1e-7,
                                          int test_print_level = 2,
                                          double test_fail_atol = 1e-5,
                                          double test_fail_rtol = 1e-5) {
  int fail = 0;

  TacsScalar Cs[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  TacsScalar u0x[9], d1x[3], d2x[3], e0ty[2];
  TacsScalar detXd;
  TacsGenerateRandomArray(Cs,
                          TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES);
  TacsGenerateRandomArray(u0x, 9);
  TacsGenerateRandomArray(d1x, 3);
  TacsGenerateRandomArray(d2x, 3);
  TacsGenerateRandomArray(e0ty, 2);
  TacsGenerateRandomArray(&detXd, 1);

  TacsScalar e[6];
  TACSBeamNonlinearModel::evalStrain(u0x, d1x, d2x, e0ty, e);

  TacsScalar s[6];
  TACSBeamConstitutive::computeStress(Cs, e, s);

  int max_err_index, max_rel_index;
  double max_err, max_rel;

  // A random direction stands in for the "path" direction getMatType
  // actually uses; the check itself only needs evalStrain/evalStrainSens/
  // evalStrainHessian's own consistency with their "Deriv" counterparts,
  // not any particular physical direction.
  TacsScalar u0xd[9], d1xd[3], d2xd[3], e0tyd[2];
  TacsGenerateRandomArray(u0xd, 9);
  TacsGenerateRandomArray(d1xd, 3);
  TacsGenerateRandomArray(d2xd, 3);
  TacsGenerateRandomArray(e0tyd, 2);

  // --- evalStrainDeriv: ed should match the directional derivative of
  // evalStrain along (u0xd, d1xd, d2xd, e0tyd). ---
  TacsScalar ed[6];
  TACSBeamNonlinearModel::evalStrainDeriv(u0x, d1x, d2x, e0ty, u0xd, d1xd, d2xd,
                                          e0tyd, e, ed);

  TacsScalar u0xt2[9], d1xt2[3], d2xt2[3], e0tyt2[2], et2[6];
  for (int i = 0; i < 9; i++) {
#ifdef TACS_USE_COMPLEX
    u0xt2[i] = u0x[i] + TacsScalar(0.0, dh) * u0xd[i];
#else
    u0xt2[i] = u0x[i] + dh * u0xd[i];
#endif  // TACS_USE_COMPLEX
  }
  for (int i = 0; i < 3; i++) {
#ifdef TACS_USE_COMPLEX
    d1xt2[i] = d1x[i] + TacsScalar(0.0, dh) * d1xd[i];
    d2xt2[i] = d2x[i] + TacsScalar(0.0, dh) * d2xd[i];
#else
    d1xt2[i] = d1x[i] + dh * d1xd[i];
    d2xt2[i] = d2x[i] + dh * d2xd[i];
#endif  // TACS_USE_COMPLEX
  }
  for (int i = 0; i < 2; i++) {
#ifdef TACS_USE_COMPLEX
    e0tyt2[i] = e0ty[i] + TacsScalar(0.0, dh) * e0tyd[i];
#else
    e0tyt2[i] = e0ty[i] + dh * e0tyd[i];
#endif  // TACS_USE_COMPLEX
  }
  TACSBeamNonlinearModel::evalStrain(u0xt2, d1xt2, d2xt2, e0tyt2, et2);

  TacsScalar fded[6];
  for (int i = 0; i < 6; i++) {
#ifdef TACS_USE_COMPLEX
    fded[i] = TacsImagPart(et2[i]) / dh;
#else
    fded[i] = (et2[i] - e[i]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(ed, fded, 6, &max_err_index);
  max_rel = TacsGetMaxRelError(ed, fded, 6, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainDeriv\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "ed", ed, fded, 6);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // --- evalStrainSensDeriv: directional derivative of evalStrainSens's
  // output along (dfded, u0xd, d1xd, d2xd, e0tyd). ---
  TacsScalar dfded[6];
  TacsGenerateRandomArray(dfded, 6);

  TacsScalar du0x2[9], dd1x2[3], dd2x2[3], de0ty2[2];
  TacsScalar du0xd[9], dd1xd[3], dd2xd[3], de0tyd[2];
  TACSBeamNonlinearModel::evalStrainSensDeriv(
      detXd, s, u0x, d1x, d2x, e0ty, dfded, u0xd, d1xd, d2xd, e0tyd, du0x2,
      dd1x2, dd2x2, de0ty2, du0xd, dd1xd, dd2xd, de0tyd);

  TacsScalar st2[6];
  for (int i = 0; i < 6; i++) {
#ifdef TACS_USE_COMPLEX
    st2[i] = s[i] + TacsScalar(0.0, dh) * dfded[i];
#else
    st2[i] = s[i] + dh * dfded[i];
#endif  // TACS_USE_COMPLEX
  }

  TacsScalar du0xt2[9], dd1xt2[3], dd2xt2[3], de0tyt2[2];
  TACSBeamNonlinearModel::evalStrainSens(
      detXd, st2, u0xt2, d1xt2, d2xt2, e0tyt2, du0xt2, dd1xt2, dd2xt2, de0tyt2);

  TacsScalar fdu0xd[9], fdd1xd[3], fdd2xd[3], fde0tyd[2];
  for (int i = 0; i < 9; i++) {
#ifdef TACS_USE_COMPLEX
    fdu0xd[i] = TacsImagPart(du0xt2[i]) / dh;
#else
    fdu0xd[i] = (du0xt2[i] - du0x2[i]) / dh;
#endif  // TACS_USE_COMPLEX
  }
  for (int i = 0; i < 3; i++) {
#ifdef TACS_USE_COMPLEX
    fdd1xd[i] = TacsImagPart(dd1xt2[i]) / dh;
    fdd2xd[i] = TacsImagPart(dd2xt2[i]) / dh;
#else
    fdd1xd[i] = (dd1xt2[i] - dd1x2[i]) / dh;
    fdd2xd[i] = (dd2xt2[i] - dd2x2[i]) / dh;
#endif  // TACS_USE_COMPLEX
  }
  for (int i = 0; i < 2; i++) {
#ifdef TACS_USE_COMPLEX
    fde0tyd[i] = TacsImagPart(de0tyt2[i]) / dh;
#else
    fde0tyd[i] = (de0tyt2[i] - de0ty2[i]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(du0xd, fdu0xd, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(du0xd, fdu0xd, 9, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainSensDeriv (du0xd)\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  max_err = TacsGetMaxError(dd1xd, fdd1xd, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd1xd, fdd1xd, 3, &max_rel_index);
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainSensDeriv (dd1xd)\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
    fprintf(stderr, "\n");
  }

  max_err = TacsGetMaxError(dd2xd, fdd2xd, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd2xd, fdd2xd, 3, &max_rel_index);
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainSensDeriv (dd2xd)\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
    fprintf(stderr, "\n");
  }

  max_err = TacsGetMaxError(de0tyd, fde0tyd, 2, &max_err_index);
  max_rel = TacsGetMaxRelError(de0tyd, fde0tyd, 2, &max_rel_index);
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainSensDeriv (de0tyd)\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
    fprintf(stderr, "\n");
  }

  // --- evalStrainHessianDeriv: directional derivative of
  // evalStrainHessian's 7 output blocks along (dfded (reused as the
  // stress-direction stand-in), u0xd, d1xd, d2xd, e0tyd). Only the
  // "d2u0xd" block (the largest, 81-entry, and the one carrying the
  // geometric constant-shape terms getMatType relies on most directly) is
  // checked here, to keep this addition bounded -- the same pattern
  // extends identically to the other 6 blocks. ---
  TacsScalar d2u0xH[81], d2d1xH[9], d2d2xH[9], d2e0tyH[4];
  TacsScalar d2u0xd1xH[27], d2u0xd2xH[27], d2d1xd2xH[9];
  TacsScalar d2u0xHd[81], d2d1xHd[9], d2d2xHd[9], d2e0tyHd[4];
  TacsScalar d2u0xd1xHd[27], d2u0xd2xHd[27], d2d1xd2xHd[9];
  TACSBeamNonlinearModel::evalStrainHessianDeriv(
      detXd, s, Cs, u0x, d1x, d2x, e0ty, dfded, u0xd, d1xd, d2xd, e0tyd, d2u0xH,
      d2d1xH, d2d2xH, d2e0tyH, d2u0xd1xH, d2u0xd2xH, d2d1xd2xH, d2u0xHd,
      d2d1xHd, d2d2xHd, d2e0tyHd, d2u0xd1xHd, d2u0xd2xHd, d2d1xd2xHd);

  TacsScalar d2u0xH_t[81], d2d1xH_t[9], d2d2xH_t[9], d2e0tyH_t[4];
  TacsScalar d2u0xd1xH_t[27], d2u0xd2xH_t[27], d2d1xd2xH_t[9];
  TACSBeamNonlinearModel::evalStrainHessian(
      detXd, st2, Cs, u0xt2, d1xt2, d2xt2, e0tyt2, d2u0xH_t, d2d1xH_t, d2d2xH_t,
      d2e0tyH_t, d2u0xd1xH_t, d2u0xd2xH_t, d2d1xd2xH_t);

  TacsScalar fd2u0xHd[81];
  for (int i = 0; i < 81; i++) {
#ifdef TACS_USE_COMPLEX
    fd2u0xHd[i] = TacsImagPart(d2u0xH_t[i]) / dh;
#else
    fd2u0xHd[i] = (d2u0xH_t[i] - d2u0xH[i]) / dh;
#endif  // TACS_USE_COMPLEX
  }

  max_err = TacsGetMaxError(d2u0xHd, fd2u0xHd, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2u0xHd, fd2u0xHd, 81, &max_rel_index);
  if (test_print_level > 0) {
    fprintf(stderr, "Testing evalStrainHessianDeriv (d2u0xd)\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2u0xHd", d2u0xHd, fd2u0xHd, 81);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }
  fail = fail || (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif  // TACS_BEAM_ELEMENT_MODEL_H
