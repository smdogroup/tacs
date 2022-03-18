#ifndef TACS_BEAM_ELEMENT_MODEL_H
#define TACS_BEAM_ELEMENT_MODEL_H

#include "TACSElementAlgebra.h"
#include "TACSBeamConstitutive.h"
#include "TACSElementVerification.h"
#include "TACSBeamElementBasis.h"

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
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi, res);
      }
    }
  }

  static inline void evalStrain( const TacsScalar u0x[],
                                 const TacsScalar d1x[],
                                 const TacsScalar d2x[],
                                 const TacsScalar e0ty[],
                                 TacsScalar e[] ){
    // Axial strain
    e[0] = u0x[0];

    // Torsional component of the strain
    e[1] = 0.5*(d1x[2] - d2x[1]);

    // Bending components of the strain
    e[2] = d1x[0];
    e[3] = d2x[0];

    // Add the tying shear strain components
    e[4] = e0ty[0];
    e[5] = e0ty[1];
  }

  static inline void evalStrainSens( const TacsScalar scale,
                                     const TacsScalar dfde[],
                                     const TacsScalar u0x[],
                                     const TacsScalar d1x[],
                                     const TacsScalar d2x[],
                                     const TacsScalar e0ty[],
                                     TacsScalar du0x[],
                                     TacsScalar dd1x[],
                                     TacsScalar dd2x[],
                                     TacsScalar de0ty[] ){
    du0x[0] = scale*dfde[0];
    du0x[1] = du0x[2] = 0.0;
    du0x[3] = du0x[4] = du0x[5] = 0.0;
    du0x[6] = du0x[7] = du0x[8] = 0.0;

    dd1x[0] = scale*dfde[2];
    dd1x[1] = 0.0;
    dd1x[2] = 0.5*scale*dfde[1];
    dd1x[3] = dd1x[4] = dd1x[5] = 0.0;
    dd1x[6] = dd1x[7] = dd1x[8] = 0.0;

    dd2x[0] = scale*dfde[3];
    dd2x[1] = -0.5*scale*dfde[1];
    dd2x[2] = 0.0;
    dd2x[3] = dd2x[4] = dd2x[5] = 0.0;
    dd2x[6] = dd2x[7] = dd2x[8] = 0.0;

    de0ty[0] = scale*dfde[4];
    de0ty[1] = scale*dfde[5];
  }

  static void evalStrainHessian( const TacsScalar scale,
                                 const TacsScalar s[],
                                 const TacsScalar Cs[],
                                 const TacsScalar u0x[],
                                 const TacsScalar d1x[],
                                 const TacsScalar d2x[],
                                 const TacsScalar e0ty[],
                                 TacsScalar d2u0x[],
                                 TacsScalar d2d1x[],
                                 TacsScalar d2d2x[],
                                 TacsScalar d2e0ty[],
                                 TacsScalar d2u0xd1x[],
                                 TacsScalar d2u0xd2x[],
                                 TacsScalar d2d1xd2x[] ){


  }
};

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
        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi, res);
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
int TacsTestBeamModelDerivatives( double dh=1e-7,
                                  int test_print_level=2,
                                  double test_fail_atol=1e-5,
                                  double test_fail_rtol=1e-5 ){
  // Set the failure flag
  int fail = 0;

  // Set random values for the constitutive data and inputs
  TacsScalar Cs[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  TacsScalar u0x[9], u1x[9], u2x[9], e0ty[2];
  TacsScalar detXd;

  // Set random data
  TacsGenerateRandomArray(Cs, TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES);
  TacsGenerateRandomArray(u0x, 9);
  TacsGenerateRandomArray(u1x, 9);
  TacsGenerateRandomArray(u2x, 9);
  TacsGenerateRandomArray(e0ty, 2);
  TacsGenerateRandomArray(&detXd, 1);

  // Compute the strain
  TacsScalar e[6];
  model::evalStrain(u0x, u1x, u2x, e0ty, e);

  // Compute the stress
  TacsScalar s[6];
  TACSBeamConstitutive::computeStress(Cs, e, s);

  // Compute the derivative of the product of the stress and strain
  // with respect to u0x, u1x and e0ty
  TacsScalar du0x[9], du1x[9], du2x[9], de0ty[2];
  model::evalStrainSens(detXd, s, u0x, u1x, u2x, e0ty, du0x, du1x, du2x, de0ty);

  TacsScalar f0 = 0.0;
  for ( int j = 0; j < 6; j++ ){
    f0 += 0.5*detXd*e[j]*s[j];
  }

  // Compute against the derivatives for the strain
  TacsScalar fdu0x[9];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u0xt[9], et[6], st[6];
    memcpy(u0xt, u0x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u0xt[i] = u0x[i] + TacsScalar(0.0, dh);
#else
    u0xt[i] = u0x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0xt, u1x, u2x, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdu0x[i] = TacsImagPart(f1)/dh;
#else
    fdu0x[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  int max_err_index, max_rel_index;
  double max_err = TacsGetMaxError(du0x, fdu0x, 9, &max_err_index);
  double max_rel = TacsGetMaxRelError(du0x, fdu0x, 9, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "du0x", du0x, fdu0x, 9);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdu1x[9];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u1xt[9], et[6], st[6];
    memcpy(u1xt, u1x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u1xt[i] = u1x[i] + TacsScalar(0.0, dh);
#else
    u1xt[i] = u1x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1xt, u2x, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdu1x[i] = TacsImagPart(f1)/dh;
#else
    fdu1x[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(du1x, fdu1x, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(du1x, fdu1x, 9, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "du1x", du1x, fdu1x, 9);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdu2x[9];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u2xt[9], et[6], st[6];
    memcpy(u2xt, u2x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u2xt[i] = u2x[i] + TacsScalar(0.0, dh);
#else
    u2xt[i] = u2x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, u2xt, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdu2x[i] = TacsImagPart(f1)/dh;
#else
    fdu2x[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(du2x, fdu2x, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(du2x, fdu2x, 9, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "du2x", du2x, fdu2x, 9);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fde0ty[2];
  for ( int i = 0; i < 2; i++ ){
    TacsScalar e0tyt[2], et[6], st[6];
    memcpy(e0tyt, e0ty, 2*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    e0tyt[i] = e0ty[i] + TacsScalar(0.0, dh);
#else
    e0tyt[i] = e0ty[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, u2x, e0tyt, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fde0ty[i] = TacsImagPart(f1)/dh;
#else
    fde0ty[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(de0ty, fde0ty, 2, &max_err_index);
  max_rel = TacsGetMaxRelError(de0ty, fde0ty, 2, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. e0ty\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "de0ty", de0ty, fde0ty, 2);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

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

#endif // TACS_BEAM_ELEMENT_MODEL_H
