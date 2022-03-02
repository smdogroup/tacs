#ifndef TACS_BEAM_ELEMENT_MODEL_H
#define TACS_BEAM_ELEMENT_MODEL_H

#include "TACSElementAlgebra.h"
#include "TACSBeamConstitutive.h"
#include "TACSElementVerification.h"

class TACSBeamLinearModel {
 public:

  /*
    Compute the tensorial components of the tying strain

  */

  /*
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




    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn[],
                                              const TacsScalar vars[],
                                              const TacsScalar d[],
                                              const TacsScalar dety[],
                                              TacsScalar res[],
                                              TacsScalar dd[] ){

  }
  */

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
    du0x[1] = 0.0;
    du0x[2] = 0.0;

    dd1x[0] = scale*dfde[2];
    dd1x[1] = 0.0;
    dd1x[2] = 0.5*scale*dfde[1];
    
    dd2x[0] = scale*dfde[3];
    dd2x[1] = -0.5*scale*dfde[1];
    dd2x[2] = 0.0;

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


/*
class TACSBeamNonlinearModel {
 public:

  template <int vars_per_node, class basis>
  static void computeTyingStrain( ){


  }


  static inline void evalStrain( const TacsScalar u0x[],
                                 const TacsScalar d1x[],
                                 const TacsScalar d2x[],
                                 const TacsScalar e0ty[],
                                 TacsScalar e[] ){
    // Axial strain
    e[0] = u0x[0] + 0.5*(u0x[0]*u0x[0] + u0x[1]*u0x[1] + u0x[2]*u0x[2]);

    // Torsional strain
    e[1] = 0.5*(d1x[2] - d2x[1]);

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
  TacsScalar u0x[3], d1x[3], d2x[3], e0ty[2];
  TacsScalar detXd;

  // Set random data
  TacsGenerateRandomArray(Cs, TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES);
  TacsGenerateRandomArray(u0x, 3);
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
  // with respect to u0x, u1x and e0ty
  TacsScalar du0x[3], dd1x[3], dd2x[3], de0ty[2];
  model::evalStrainSens(detXd, s, u0x, d1x, d2x, e0ty, du0x, dd1x, dd2x, de0ty);

  TacsScalar f0 = 0.0;
  for ( int j = 0; j < 6; j++ ){
    f0 += 0.5*detXd*e[j]*s[j];
  }

  // Compute against the derivatives for the strain
  TacsScalar fdu0x[3];
  for ( int i = 0; i < 3; i++ ){
    TacsScalar u0xt[3], et[6], st[6];
    memcpy(u0xt, u0x, 3*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u0xt[i] = u0x[i] + TacsScalar(0.0, dh);
#else
    u0xt[i] = u0x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0xt, d1x, d2x, e0ty, et);
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
  double max_err = TacsGetMaxError(du0x, fdu0x, 3, &max_err_index);
  double max_rel = TacsGetMaxRelError(du0x, fdu0x, 3, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "du0x", du0x, fdu0x, 3);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdd1x[3];
  for ( int i = 0; i < 3; i++ ){
    TacsScalar d1xt[3], et[6], st[6];
    memcpy(d1xt, d1x, 3*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    d1xt[i] = d1x[i] + TacsScalar(0.0, dh);
#else
    d1xt[i] = d1x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, d1xt, d2x, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdd1x[i] = TacsImagPart(f1)/dh;
#else
    fdd1x[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dd1x, fdd1x, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd1x, fdd1x, 3, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "dd1x", dd1x, fdd1x, 3);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdd2x[3];
  for ( int i = 0; i < 3; i++ ){
    TacsScalar d2xt[3], et[6], st[6];
    memcpy(d2xt, d2x, 3*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    d2xt[i] = d2x[i] + TacsScalar(0.0, dh);
#else
    d2xt[i] = d2x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, d1x, d2xt, e0ty, et);
    TACSBeamConstitutive::computeStress(Cs, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 6; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdd2x[i] = TacsImagPart(f1)/dh;
#else
    fdd2x[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dd2x, fdd2x, 3, &max_err_index);
  max_rel = TacsGetMaxRelError(dd2x, fdd2x, 3, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. u1x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "dd2x", dd2x, fdd2x, 3);
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
    model::evalStrain(u0x, d1x, d2x, e0tyt, et);
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

