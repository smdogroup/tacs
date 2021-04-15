#ifndef TACS_SHELL_ELEMENT_MODEL_H
#define TACS_SHELL_ELEMENT_MODEL_H

#include "TACSShellConstitutive.h"
#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"

class TACSShellLinearModel {
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
  template <class basis>
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn[],
                                  const int vars_per_node,
                                  const TacsScalar vars[],
                                  const TacsScalar d[],
                                  TacsScalar ety[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);

        ety[index] = 0.0;
        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          ety[index] = (Uxi[0]*Xxi[0] + Uxi[2]*Xxi[2] + Uxi[4]*Xxi[4]);
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          ety[index] = (Uxi[1]*Xxi[1] + Uxi[3]*Xxi[3] + Uxi[5]*Xxi[5]);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          ety[index] = 0.5*(Uxi[0]*Xxi[1] + Uxi[2]*Xxi[3] + Uxi[4]*Xxi[5] +
                            Uxi[1]*Xxi[0] + Uxi[3]*Xxi[2] + Uxi[5]*Xxi[4]);
        }
        else {
          TacsScalar d0[3], n0[3];
          basis::interpFields(pt, 3, d, 3, d0);
          basis::interpFields(pt, 3, fn, 3, n0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                              n0[0]*Uxi[1] + n0[1]*Uxi[3] + n0[2]*Uxi[5]);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                              n0[0]*Uxi[0] + n0[1]*Uxi[2] + n0[2]*Uxi[4]);
          }
        }
      }
    }
  }

  template <class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn[],
                                              const int vars_per_node,
                                              const TacsScalar vars[],
                                              const TacsScalar d[],
                                              const TacsScalar dety[],
                                              TacsScalar res[],
                                              TacsScalar dd[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6], dUxi[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);

        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          dUxi[0] = dety[index]*Xxi[0];
          dUxi[1] = 0.0;
          dUxi[2] = dety[index]*Xxi[2];
          dUxi[3] = 0.0;
          dUxi[4] = dety[index]*Xxi[4];
          dUxi[5] = 0.0;

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          dUxi[0] = 0.0;
          dUxi[1] = dety[index]*Xxi[1];
          dUxi[2] = 0.0;
          dUxi[3] = dety[index]*Xxi[3];
          dUxi[4] = 0.0;
          dUxi[5] = dety[index]*Xxi[5];

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          dUxi[0] = 0.5*dety[index]*Xxi[1];
          dUxi[1] = 0.5*dety[index]*Xxi[0];
          dUxi[2] = 0.5*dety[index]*Xxi[3];
          dUxi[3] = 0.5*dety[index]*Xxi[2];
          dUxi[4] = 0.5*dety[index]*Xxi[5];
          dUxi[5] = 0.5*dety[index]*Xxi[4];

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else {
          TacsScalar d0[3], dd0[3], n0[3];
          basis::interpFields(pt, 3, d, 3, d0);
          basis::interpFields(pt, 3, fn, 3, n0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi[0] = 0.0;
            dUxi[1] = 0.5*dety[index]*n0[0];
            dUxi[2] = 0.0;
            dUxi[3] = 0.5*dety[index]*n0[1];
            dUxi[4] = 0.0;
            dUxi[5] = 0.5*dety[index]*n0[2];
            basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);

            dd0[0] = 0.5*dety[index]*Xxi[1];
            dd0[1] = 0.5*dety[index]*Xxi[3];
            dd0[2] = 0.5*dety[index]*Xxi[5];
            basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi[0] = 0.5*dety[index]*n0[0];
            dUxi[1] = 0.0;
            dUxi[2] = 0.5*dety[index]*n0[1];
            dUxi[3] = 0.0;
            dUxi[4] = 0.5*dety[index]*n0[2];
            dUxi[5] = 0.0;
            basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);

            dd0[0] = 0.5*dety[index]*Xxi[0];
            dd0[1] = 0.5*dety[index]*Xxi[2];
            dd0[2] = 0.5*dety[index]*Xxi[4];
            basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
          }
        }
      }
    }
  }

  template <class basis>
  static void addComputeTyingStrainHessian( const TacsScalar Xpts[],
                                            const TacsScalar fn[],
                                            const int vars_per_node,
                                            const TacsScalar vars[],
                                            const TacsScalar d[],
                                            const TacsScalar d2ety[],
                                            TacsScalar mat[],
                                            TacsScalar d2d[],
                                            TacsScalar d2du[] ){
    const int num_tying_fields = 5;
    for ( int f1 = 0; f1 < num_tying_fields; f1++ ){
      const int nty1 = basis::getNumTyingPoints(f1);

      for ( int f2 = 0; f2 < num_tying_fields; f2++ ){
        const int nty2 = basis::getNumTyingPoints(f2);

        for ( int ty1 = 0; ty1 < nty1; ty1++ ){
          double pt1[2];
          basis::getTyingPoint(f1, ty1, pt1);

          TacsScalar Xxi1[6];
          basis::interpFieldsGrad(pt1, 3, Xpts, 3, Xxi1);

          // Store the the derivative information for the first point
          TacsScalar dUxi1[6], dd01[3];

          if (f1 == 0){
            // Compute g11 = e1^{T}*G*e1
            dUxi1[0] = Xxi1[0];
            dUxi1[1] = 0.0;
            dUxi1[2] = Xxi1[2];
            dUxi1[3] = 0.0;
            dUxi1[4] = Xxi1[4];
            dUxi1[5] = 0.0;
          }
          else if (f1 == 1){
            // Compute g22 = e2^{T}*G*e2
            dUxi1[0] = 0.0;
            dUxi1[1] = Xxi1[1];
            dUxi1[2] = 0.0;
            dUxi1[3] = Xxi1[3];
            dUxi1[4] = 0.0;
            dUxi1[5] = Xxi1[5];
          }
          else if (f1 == 2){
            // Compute g12 = e2^{T}*G*e1
            dUxi1[0] = 0.5*Xxi1[1];
            dUxi1[1] = 0.5*Xxi1[0];
            dUxi1[2] = 0.5*Xxi1[3];
            dUxi1[3] = 0.5*Xxi1[2];
            dUxi1[4] = 0.5*Xxi1[5];
            dUxi1[5] = 0.5*Xxi1[4];
          }
          else {
            TacsScalar n0[3];
            basis::interpFields(pt1, 3, fn, 3, n0);

            if (f1 == 3){
              // Compute g23 = e2^{T}*G*e3
              dUxi1[0] = 0.0;
              dUxi1[1] = 0.5*n0[0];
              dUxi1[2] = 0.0;
              dUxi1[3] = 0.5*n0[1];
              dUxi1[4] = 0.0;
              dUxi1[5] = 0.5*n0[2];

              dd01[0] = 0.5*Xxi1[1];
              dd01[1] = 0.5*Xxi1[3];
              dd01[2] = 0.5*Xxi1[5];
            }
            else if (f1 == 4){
              // Compute g13 = e1^{T}*G*e3
              dUxi1[0] = 0.5*n0[0];
              dUxi1[1] = 0.0;
              dUxi1[2] = 0.5*n0[1];
              dUxi1[3] = 0.0;
              dUxi1[4] = 0.5*n0[2];
              dUxi1[5] = 0.0;

              dd01[0] = 0.5*Xxi1[0];
              dd01[1] = 0.5*Xxi1[2];
              dd01[2] = 0.5*Xxi1[4];
            }
          }

          for ( int ty2 = 0; ty2 < nty2; ty2++, d2ety++ ){
            double pt2[2];
            basis::getTyingPoint(f2, ty2, pt2);

            TacsScalar Xxi2[6];
            basis::interpFieldsGrad(pt2, 3, Xpts, 3, Xxi2);

            // Store the the derivative information for the second point
            TacsScalar dUxi2[6], dd02[3];

            if (f2 == 0){
              // Compute g11 = e1^{T}*G*e1
              dUxi2[0] = Xxi2[0];
              dUxi2[1] = 0.0;
              dUxi2[2] = Xxi2[2];
              dUxi2[3] = 0.0;
              dUxi2[4] = Xxi2[4];
              dUxi2[5] = 0.0;
            }
            else if (f2 == 1){
              // Compute g22 = e2^{T}*G*e2
              dUxi2[0] = 0.0;
              dUxi2[1] = Xxi2[1];
              dUxi2[2] = 0.0;
              dUxi2[3] = Xxi2[3];
              dUxi2[4] = 0.0;
              dUxi2[5] = Xxi2[5];
            }
            else if (f2 == 2){
              // Compute g12 = e2^{T}*G*e1
              dUxi2[0] = 0.5*Xxi2[1];
              dUxi2[1] = 0.5*Xxi2[0];
              dUxi2[2] = 0.5*Xxi2[3];
              dUxi2[3] = 0.5*Xxi2[2];
              dUxi2[4] = 0.5*Xxi2[5];
              dUxi2[5] = 0.5*Xxi2[4];
            }
            else {
              TacsScalar n0[3];
              basis::interpFields(pt2, 3, fn, 3, n0);

              if (f2 == 3){
                // Compute g23 = e2^{T}*G*e3
                dUxi2[0] = 0.0;
                dUxi2[1] = 0.5*n0[0];
                dUxi2[2] = 0.0;
                dUxi2[3] = 0.5*n0[1];
                dUxi2[4] = 0.0;
                dUxi2[5] = 0.5*n0[2];

                dd02[0] = 0.5*Xxi2[1];
                dd02[1] = 0.5*Xxi2[3];
                dd02[2] = 0.5*Xxi2[5];
              }
              else if (f2 == 4){
                // Compute g13 = e1^{T}*G*e3
                dUxi2[0] = 0.5*n0[0];
                dUxi2[1] = 0.0;
                dUxi2[2] = 0.5*n0[1];
                dUxi2[3] = 0.0;
                dUxi2[4] = 0.5*n0[2];
                dUxi2[5] = 0.0;

                dd02[0] = 0.5*Xxi2[0];
                dd02[1] = 0.5*Xxi2[2];
                dd02[2] = 0.5*Xxi2[4];
              }
            }

            TacsScalar d2Uxi[36];
            for ( int i = 0; i < 6; i++ ){
              for ( int j = 0; j < 6; j++ ){
                d2Uxi[6*i + j] = d2ety[0]*dUxi1[i]*dUxi2[j];
              }
            }
            basis::addInterpGradOuterProduct(pt1, pt2, 3, d2Uxi, vars_per_node, mat);

            if (f1 >= 3){
              TacsScalar d2d0Uxi[18];
              for ( int i = 0; i < 3; i++ ){
                for ( int j = 0; j < 6; j++ ){
                  d2d0Uxi[6*i + j] = d2ety[0]*dd01[i]*dUxi2[j];
                }
              }
              basis::addInterpGradMixedOuterProduct(pt1, pt2, 3, d2d0Uxi, NULL, 3, d2du);
            }
            if (f1 >= 3 && f2 >= 3){
              TacsScalar d2d0[9];
              for ( int i = 0; i < 3; i++ ){
                for ( int j = 0; j < 3; j++ ){
                  d2d0[3*i + j] = d2ety[0]*dd01[i]*dd02[j];
                }
              }
              basis::addInterpFieldsOuterProduct(pt1, pt2, 3, d2d0, 3, d2d);
            }
          }
        }
      }
    }
  }

  /*
    Compute the directional derivative
  */
  template <class basis>
  static void computeTyingStrainDeriv( const TacsScalar Xpts[],
                                       const TacsScalar fn[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       const TacsScalar d[],
                                       const TacsScalar varsd[],
                                       const TacsScalar dd[],
                                       TacsScalar ety[],
                                       TacsScalar etyd[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6], Uxid[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);
        basis::interpFieldsGrad(pt, vars_per_node, varsd, 3, Uxid);

        ety[index] = 0.0;
        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          ety[index] = (Uxi[0]*Xxi[0] + Uxi[2]*Xxi[2] + Uxi[4]*Xxi[4]);
          etyd[index] = (Uxid[0]*Xxi[0] + Uxid[2]*Xxi[2] + Uxid[4]*Xxi[4]);
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          ety[index] = (Uxi[1]*Xxi[1] + Uxi[3]*Xxi[3] + Uxi[5]*Xxi[5]);
          etyd[index] = (Uxid[1]*Xxi[1] + Uxid[3]*Xxi[3] + Uxid[5]*Xxi[5]);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          ety[index] = 0.5*(Uxi[0]*Xxi[1] + Uxi[2]*Xxi[3] + Uxi[4]*Xxi[5] +
                            Uxi[1]*Xxi[0] + Uxi[3]*Xxi[2] + Uxi[5]*Xxi[4]);
          etyd[index] = 0.5*(Uxid[0]*Xxi[1] + Uxid[2]*Xxi[3] + Uxid[4]*Xxi[5] +
                             Uxid[1]*Xxi[0] + Uxid[3]*Xxi[2] + Uxid[5]*Xxi[4]);
        }
        else {
          TacsScalar d0[3], d0d[3], n0[3];
          basis::interpFields(pt, 3, d, 3, d0);
          basis::interpFields(pt, 3, dd, 3, d0d);
          basis::interpFields(pt, 3, fn, 3, n0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                              n0[0]*Uxi[1] + n0[1]*Uxi[3] + n0[2]*Uxi[5]);
            etyd[index] = 0.5*(Xxi[1]*d0d[0] + Xxi[3]*d0d[1] + Xxi[5]*d0d[2] +
                               n0[0]*Uxid[1] + n0[1]*Uxid[3] + n0[2]*Uxid[5]);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                              n0[0]*Uxi[0] + n0[1]*Uxi[2] + n0[2]*Uxi[4]);
            etyd[index] = 0.5*(Xxi[0]*d0d[0] + Xxi[2]*d0d[1] + Xxi[4]*d0d[2] +
                               n0[0]*Uxid[0] + n0[1]*Uxid[2] + n0[2]*Uxid[4]);
          }
        }
      }
    }
  }

  /*
    Evaluate the strain as a function of the displacement derivatives
    and interpolated strain from the tensorial components
  */
  static void evalStrain( const TacsScalar u0x[],
                          const TacsScalar u1x[],
                          const TacsScalar e0ty[],
                          const TacsScalar Ct[],
                          TacsScalar e[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0];
    e[1] = e0ty[3];
    e[2] = 2.0*e0ty[1];

    // Compute the bending strain
    e[3] = u1x[0];
    e[4] = u1x[4];
    e[5] = u1x[1] + u1x[3];

    // Add the components of the shear strain
    e[6] = 2.0*e0ty[4];
    e[7] = 2.0*e0ty[2];

    // Compute the rotational penalty
    e[8] = 0.5*(Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  /**
    Evaluate the derivative of the strain
  */
  static void evalStrainSens( const TacsScalar scale,
                              const TacsScalar dfde[],
                              const TacsScalar u0x[],
                              const TacsScalar u1x[],
                              const TacsScalar Ct[],
                              TacsScalar du0x[],
                              TacsScalar du1x[],
                              TacsScalar de0ty[],
                              TacsScalar dCt[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = scale*dfde[0];
    de0ty[1] = 2.0*scale*dfde[2];
    de0ty[2] = 2.0*scale*dfde[7];
    de0ty[3] = scale*dfde[1];
    de0ty[4] = 2.0*scale*dfde[6];
    de0ty[5] = 0.0;

    // Linear strain relationships
    // Derivative with respect to u0x
    du0x[0] = 0.0;
    du0x[1] = -0.5*scale*dfde[8];
    du0x[2] = 0.0;
    du0x[3] = 0.5*scale*dfde[8];
    du0x[4] = 0.0;
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;

    // Compute the derivative with respect to U1
    du1x[0] = scale*dfde[3];
    du1x[1] = scale*dfde[5];
    du1x[2] = 0.0;
    du1x[3] = scale*dfde[5];
    du1x[4] = scale*dfde[4];
    du1x[5] = 0.0;
    du1x[6] = 0.0;
    du1x[7] = 0.0;
    du1x[8] = 0.0;

    dCt[0] = 0.0;
    dCt[1] = -0.5*scale*dfde[8];
    dCt[2] = 0.0;
    dCt[3] = 0.5*scale*dfde[8];
    dCt[4] = 0.0;
    dCt[5] = 0.0;
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;
  }

  static void evalStrainDeriv( const TacsScalar u0x[],
                               const TacsScalar u1x[],
                               const TacsScalar e0ty[],
                               const TacsScalar Ct[],
                               const TacsScalar u0xd[],
                               const TacsScalar u1xd[],
                               const TacsScalar e0tyd[],
                               const TacsScalar Ctd[],
                               TacsScalar e[],
                               TacsScalar ed[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0];
    e[1] = e0ty[3];
    e[2] = 2.0*e0ty[1];

    // Compute the bending strain
    e[3] = u1x[0];
    e[4] = u1x[4];
    e[5] = u1x[1] + u1x[3];

    // Add the components of the shear strain
    e[6] = 2.0*e0ty[4];
    e[7] = 2.0*e0ty[2];

    // Compute the rotational penalty
    e[8] = 0.5*(Ct[3] + u0x[3] - Ct[1] - u0x[1]);

    // Evaluate the in-plane strains from the tying strain expressions
    ed[0] = e0tyd[0];
    ed[1] = e0tyd[3];
    ed[2] = 2.0*e0tyd[1];

    // Compute the bending strain
    ed[3] = u1xd[0];
    ed[4] = u1xd[4];
    ed[5] = u1xd[1] + u1xd[3];

    // Add the components of the shear strain
    ed[6] = 2.0*e0tyd[4];
    ed[7] = 2.0*e0tyd[2];

    // Compute the rotational penalty
    ed[8] = 0.5*(Ctd[3] + u0xd[3] - Ctd[1] - u0xd[1]);
  }

  static void evalStrainHessian( const TacsScalar scale,
                                 const TacsScalar dfde[],
                                 const TacsScalar Cs[],
                                 const TacsScalar u0x[],
                                 const TacsScalar u1x[],
                                 const TacsScalar Ct[],
                                 TacsScalar d2u0x[],
                                 TacsScalar d2u1x[],
                                 TacsScalar d2u0xu1x[],
                                 TacsScalar d2e0ty[],
                                 TacsScalar d2e0tyu0x[],
                                 TacsScalar d2e0tyu1x[],
                                 TacsScalar d2Ct[],
                                 TacsScalar d2Ctu0x[] ){
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    memset(d2u0x, 0, 81*sizeof(TacsScalar));
    memset(d2u1x, 0, 81*sizeof(TacsScalar));
    memset(d2u0xu1x, 0, 81*sizeof(TacsScalar));
    memset(d2e0ty, 0, 36*sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54*sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54*sizeof(TacsScalar));
    memset(d2Ct, 0, 81*sizeof(TacsScalar));
    memset(d2Ctu0x, 0, 81*sizeof(TacsScalar));

    // Compute the second derivatives
    // e[3] = u1x[0];
    // e[4] = u1x[4];
    // e[5] = u1x[1] + u1x[3];
    TacsScalar *d2;
    d2 = d2u1x;
    d2[0] = scale*D[0];
    d2[4] = scale*D[1];
    d2[1] = scale*D[2];
    d2[3] = scale*D[2];

    d2 = &d2u1x[4*9];
    d2[0] = scale*D[1];
    d2[4] = scale*D[3];
    d2[1] = scale*D[4];
    d2[3] = scale*D[4];

    d2 = &d2u1x[9];
    d2[0] = scale*D[2];
    d2[4] = scale*D[4];
    d2[1] = scale*D[5];
    d2[3] = scale*D[5];

    d2 = &d2u1x[3*9];
    d2[0] = scale*D[2];
    d2[4] = scale*D[4];
    d2[1] = scale*D[5];
    d2[3] = scale*D[5];

    // Evaluate the in-plane strains from the tying strain expressions
    // e[0] = e0ty[0];
    // e[1] = e0ty[3];
    // e[2] = 2.0*e0ty[1];
    d2 = &d2e0ty[0];
    d2[0] = scale*A[0];
    d2[3] = scale*A[1];
    d2[1] = 2.0*scale*A[2];

    d2 = &d2e0ty[3*6];
    d2[0] = scale*A[1];
    d2[3] = scale*A[3];
    d2[1] = 2.0*scale*A[4];

    d2 = &d2e0ty[6];
    d2[0] = 2.0*scale*A[2];
    d2[3] = 2.0*scale*A[4];
    d2[1] = 4.0*scale*A[5];

    // e[6] = 2.0*e0ty[4];
    // e[7] = 2.0*e0ty[2];
    d2 = &d2e0ty[4*6];
    d2[4] = 4.0*scale*As[0];
    d2[2] = 4.0*scale*As[1];

    d2 = &d2e0ty[2*6];
    d2[4] = 4.0*scale*As[1];
    d2[2] = 4.0*scale*As[2];

    // Evaluate the cross-coupling derivatives
    d2 = &d2e0tyu1x[0];
    d2[0] = scale*B[0];
    d2[4] = scale*B[1];
    d2[1] = scale*B[2];
    d2[3] = scale*B[2];

    d2 = &d2e0tyu1x[3*9];
    d2[0] = scale*B[1];
    d2[4] = scale*B[3];
    d2[1] = scale*B[4];
    d2[3] = scale*B[4];

    d2 = &d2e0tyu1x[9];
    d2[0] = 2.0*scale*B[2];
    d2[4] = 2.0*scale*B[4];
    d2[1] = 2.0*scale*B[5];
    d2[3] = 2.0*scale*B[5];

    // Compute the contribution from the drilling strain
    // e[8] = 0.5*(Ct[3] + u0x[3] - Ct[1] - u0x[1]);
    d2 = &d2Ct[3*9];
    d2[3] = 0.25*scale*drill;
    d2[1] = -0.25*scale*drill;

    d2 = &d2Ct[9];
    d2[3] = -0.25*scale*drill;
    d2[1] = 0.25*scale*drill;

    d2 = &d2u0x[3*9];
    d2[3] = 0.25*scale*drill;
    d2[1] = -0.25*scale*drill;

    d2 = &d2u0x[9];
    d2[3] = -0.25*scale*drill;
    d2[1] = 0.25*scale*drill;

    d2 = &d2Ctu0x[3*9];
    d2[3] = 0.25*scale*drill;
    d2[1] = -0.25*scale*drill;

    d2 = &d2Ctu0x[9];
    d2[3] = -0.25*scale*drill;
    d2[1] = 0.25*scale*drill;
  }
};


class TACSShellNonlinearModel {
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
  template <class basis>
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn[],
                                  const int vars_per_node,
                                  const TacsScalar vars[],
                                  const TacsScalar d[],
                                  TacsScalar ety[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);

        ety[index] = 0.0;
        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          ety[index] = (Uxi[0]*Xxi[0] + Uxi[2]*Xxi[2] + Uxi[4]*Xxi[4] +
                        0.5*(Uxi[0]*Uxi[0] + Uxi[2]*Uxi[2] + Uxi[4]*Uxi[4]));
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          ety[index] = (Uxi[1]*Xxi[1] + Uxi[3]*Xxi[3] + Uxi[5]*Xxi[5] +
                        0.5*(Uxi[1]*Uxi[1] + Uxi[3]*Uxi[3] + Uxi[5]*Uxi[5]));
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          ety[index] = 0.5*(Uxi[0]*Xxi[1] + Uxi[2]*Xxi[3] + Uxi[4]*Xxi[5] +
                            Uxi[1]*Xxi[0] + Uxi[3]*Xxi[2] + Uxi[5]*Xxi[4] +
                            Uxi[0]*Uxi[1] + Uxi[2]*Uxi[3] + Uxi[4]*Uxi[5]);
        }
        else {
          TacsScalar n0[3], d0[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                              (n0[0] + d0[0])*Uxi[1] +
                              (n0[1] + d0[1])*Uxi[3] +
                              (n0[2] + d0[2])*Uxi[5]);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                              (n0[0] + d0[0])*Uxi[0] +
                              (n0[1] + d0[1])*Uxi[2] +
                              (n0[2] + d0[2])*Uxi[4]);
          }
        }
      }
    }
  }

  template <class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn[],
                                              const int vars_per_node,
                                              const TacsScalar vars[],
                                              const TacsScalar d[],
                                              const TacsScalar dety[],
                                              TacsScalar res[],
                                              TacsScalar dd[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6], dUxi[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);

        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          dUxi[0] = dety[index]*(Xxi[0] + Uxi[0]);
          dUxi[1] = 0.0;
          dUxi[2] = dety[index]*(Xxi[2] + Uxi[2]);
          dUxi[3] = 0.0;
          dUxi[4] = dety[index]*(Xxi[4] + Uxi[4]);
          dUxi[5] = 0.0;

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          dUxi[0] = 0.0;
          dUxi[1] = dety[index]*(Xxi[1] + Uxi[1]);
          dUxi[2] = 0.0;
          dUxi[3] = dety[index]*(Xxi[3] + Uxi[3]);
          dUxi[4] = 0.0;
          dUxi[5] = dety[index]*(Xxi[5] + Uxi[5]);

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          dUxi[0] = 0.5*dety[index]*(Xxi[1] + Uxi[1]);
          dUxi[1] = 0.5*dety[index]*(Xxi[0] + Uxi[0]);
          dUxi[2] = 0.5*dety[index]*(Xxi[3] + Uxi[3]);
          dUxi[3] = 0.5*dety[index]*(Xxi[2] + Uxi[2]);
          dUxi[4] = 0.5*dety[index]*(Xxi[5] + Uxi[5]);
          dUxi[5] = 0.5*dety[index]*(Xxi[4] + Uxi[4]);

          basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);
        }
        else {
          TacsScalar n0[3], d0[3], dd0[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi[0] = 0.0;
            dUxi[1] = 0.5*dety[index]*(n0[0] + d0[0]);
            dUxi[2] = 0.0;
            dUxi[3] = 0.5*dety[index]*(n0[1] + d0[1]);
            dUxi[4] = 0.0;
            dUxi[5] = 0.5*dety[index]*(n0[2] + d0[2]);
            basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);

            dd0[0] = 0.5*dety[index]*(Xxi[1] + Uxi[1]);
            dd0[1] = 0.5*dety[index]*(Xxi[3] + Uxi[3]);
            dd0[2] = 0.5*dety[index]*(Xxi[5] + Uxi[5]);
            basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi[0] = 0.5*dety[index]*(n0[0] + d0[0]);
            dUxi[1] = 0.0;
            dUxi[2] = 0.5*dety[index]*(n0[1] + d0[1]);
            dUxi[3] = 0.0;
            dUxi[4] = 0.5*dety[index]*(n0[2] + d0[2]);
            dUxi[5] = 0.0;
            basis::addInterpFieldsGradTranspose(pt, 3, dUxi, vars_per_node, res);

            dd0[0] = 0.5*dety[index]*(Xxi[0] + Uxi[0]);
            dd0[1] = 0.5*dety[index]*(Xxi[2] + Uxi[2]);
            dd0[2] = 0.5*dety[index]*(Xxi[4] + Uxi[4]);
            basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
          }
        }
      }
    }
  }

  template <class basis>
  static void computeTyingStrainDeriv( const TacsScalar Xpts[],
                                       const TacsScalar fn[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       const TacsScalar d[],
                                       const TacsScalar varsd[],
                                       const TacsScalar dd[],
                                       TacsScalar ety[],
                                       TacsScalar etyd[] ){
    const int num_tying_fields = 5;
    for ( int field = 0, index = 0; field < num_tying_fields; field++ ){
      const int num_tying_points = basis::getNumTyingPoints(field);
      for ( int ty = 0; ty < num_tying_points; ty++, index++ ){
        double pt[2];
        basis::getTyingPoint(field, ty, pt);

        // Interpolate the field value
        TacsScalar Uxi[6], Xxi[6], Uxid[6];
        basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
        basis::interpFieldsGrad(pt, vars_per_node, vars, 3, Uxi);
        basis::interpFieldsGrad(pt, vars_per_node, varsd, 3, Uxid);

        ety[index] = 0.0;
        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          ety[index] = (Uxi[0]*Xxi[0] + Uxi[2]*Xxi[2] + Uxi[4]*Xxi[4] +
                        0.5*(Uxi[0]*Uxi[0] + Uxi[2]*Uxi[2] + Uxi[4]*Uxi[4]));
          etyd[index] = (Uxid[0]*Xxi[0] + Uxid[2]*Xxi[2] + Uxid[4]*Xxi[4] +
                         Uxi[0]*Uxid[0] + Uxi[2]*Uxid[2] + Uxi[4]*Uxid[4]);
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          ety[index] = (Uxi[1]*Xxi[1] + Uxi[3]*Xxi[3] + Uxi[5]*Xxi[5] +
                        0.5*(Uxi[1]*Uxi[1] + Uxi[3]*Uxi[3] + Uxi[5]*Uxi[5]));
          etyd[index] = (Uxid[1]*Xxi[1] + Uxid[3]*Xxi[3] + Uxid[5]*Xxi[5] +
                         Uxi[1]*Uxid[1] + Uxi[3]*Uxid[3] + Uxi[5]*Uxid[5]);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          ety[index] = 0.5*(Uxi[0]*Xxi[1] + Uxi[2]*Xxi[3] + Uxi[4]*Xxi[5] +
                            Uxi[1]*Xxi[0] + Uxi[3]*Xxi[2] + Uxi[5]*Xxi[4] +
                            Uxi[0]*Uxi[1] + Uxi[2]*Uxi[3] + Uxi[4]*Uxi[5]);
          etyd[index] = 0.5*(Uxid[0]*Xxi[1] + Uxid[2]*Xxi[3] + Uxid[4]*Xxi[5] +
                             Uxid[1]*Xxi[0] + Uxid[3]*Xxi[2] + Uxid[5]*Xxi[4] +
                             Uxid[0]*Uxi[1] + Uxid[2]*Uxi[3] + Uxid[4]*Uxi[5] +
                             Uxi[0]*Uxid[1] + Uxi[2]*Uxid[3] + Uxi[4]*Uxid[5]);
        }
        else {
          TacsScalar n0[3], d0[3], d0d[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);
          basis::interpFields(pt, 3, dd, 3, d0d);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                              (n0[0] + d0[0])*Uxi[1] +
                              (n0[1] + d0[1])*Uxi[3] +
                              (n0[2] + d0[2])*Uxi[5]);
            etyd[index] = 0.5*(Xxi[1]*d0d[0] + Xxi[3]*d0d[1] + Xxi[5]*d0d[2] +
                               (n0[0] + d0[0])*Uxid[1] + d0d[0]*Uxi[1] +
                               (n0[1] + d0[1])*Uxid[3] + d0d[1]*Uxi[3] +
                               (n0[2] + d0[2])*Uxid[5] + d0d[2]*Uxi[5]);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                              (n0[0] + d0[0])*Uxi[0] +
                              (n0[1] + d0[1])*Uxi[2] +
                              (n0[2] + d0[2])*Uxi[4]);
            etyd[index] = 0.5*(Xxi[0]*d0d[0] + Xxi[2]*d0d[1] + Xxi[4]*d0d[2] +
                              (n0[0] + d0[0])*Uxid[0] + d0d[0]*Uxi[0] +
                              (n0[1] + d0[1])*Uxid[2] + d0d[1]*Uxi[2] +
                              (n0[2] + d0[2])*Uxid[4] + d0d[2]*Uxi[4]);
          }
        }
      }
    }
  }

  /*
    Evaluate the strain as a function of the displacement derivatives
    and interpolated strain from the tensorial components
  */
  static void evalStrain( const TacsScalar u0x[],
                          const TacsScalar u1x[],
                          const TacsScalar e0ty[],
                          const TacsScalar Ct[],
                          TacsScalar e[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0];
    e[1] = e0ty[3];
    e[2] = 2.0*e0ty[1];

    // Compute the bending strain
    e[3] = u1x[0] + (u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6]);
    e[4] = u1x[4] + (u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7]);
    e[5] = u1x[1] + u1x[3] + (u0x[0]*u1x[1] + u0x[3]*u1x[4] + u0x[6]*u1x[7] +
                              u1x[0]*u0x[1] + u1x[3]*u0x[4] + u1x[6]*u0x[7]);

    // Add the components of the shear strain
    e[6] = 2.0*e0ty[4];
    e[7] = 2.0*e0ty[2];

    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    e[8] =
      (Ct[3]*(1.0 + u0x[0]) + Ct[4]*u0x[3] + Ct[5]*u0x[6]) -
      (Ct[0]*u0x[1] + Ct[1]*(1.0 + u0x[4]) + Ct[2]*u0x[7]);
  }

  /**
    Evaluate the derivative of the strain
  */
  static void evalStrainSens( const TacsScalar scale,
                              const TacsScalar dfde[],
                              const TacsScalar u0x[],
                              const TacsScalar u1x[],
                              const TacsScalar Ct[],
                              TacsScalar du0x[],
                              TacsScalar du1x[],
                              TacsScalar de0ty[],
                              TacsScalar dCt[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = scale*dfde[0];
    de0ty[1] = 2.0*scale*dfde[2];
    de0ty[2] = 2.0*scale*dfde[7];
    de0ty[3] = scale*dfde[1];
    de0ty[4] = 2.0*scale*dfde[6];
    de0ty[5] = 0.0;

    // Derivative with respect to u0x
    du0x[0] = scale*((dfde[3]*u1x[0] + dfde[5]*u1x[1]) + Ct[3]*dfde[8]);
    du0x[1] = scale*((dfde[4]*u1x[1] + dfde[5]*u1x[0]) - Ct[0]*dfde[8]);
    du0x[2] = 0.0;
    du0x[3] = scale*((dfde[3]*u1x[3] + dfde[5]*u1x[4]) + Ct[4]*dfde[8]);
    du0x[4] = scale*((dfde[4]*u1x[4] + dfde[5]*u1x[3]) - Ct[1]*dfde[8]);
    du0x[5] = 0.0;
    du0x[6] = scale*((dfde[3]*u1x[6] + dfde[5]*u1x[7]) + Ct[5]*dfde[8]);
    du0x[7] = scale*((dfde[4]*u1x[7] + dfde[5]*u1x[6]) - Ct[2]*dfde[8]);
    du0x[8] = 0.0;

    // Compute the derivative with respect to U1
    du1x[0] = scale*(dfde[3]*(1.0 + u0x[0]) + dfde[5]*u0x[1]);
    du1x[1] = scale*(dfde[5]*(1.0 + u0x[0]) + dfde[4]*u0x[1]);
    du1x[2] = 0.0;
    du1x[3] = scale*(dfde[5]*(1.0 + u0x[4]) + dfde[3]*u0x[3]);
    du1x[4] = scale*(dfde[4]*(1.0 + u0x[4]) + dfde[5]*u0x[3]);
    du1x[5] = 0.0;
    du1x[6] = scale*(dfde[3]*u0x[6] + dfde[5]*u0x[7]);
    du1x[7] = scale*(dfde[4]*u0x[7] + dfde[5]*u0x[6]);
    du1x[8] = 0.0;

    dCt[0] = -scale*u0x[0]*dfde[8];
    dCt[1] = -scale*(1.0 + u0x[4])*dfde[8];
    dCt[2] = -scale*u0x[7]*dfde[8];
    dCt[3] = scale*(1.0 + u0x[0])*dfde[8];
    dCt[4] = scale*u0x[3]*dfde[8];
    dCt[5] = scale*u0x[6]*dfde[8];
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;
  }

  static void evalStrainDeriv( const TacsScalar u0x[],
                               const TacsScalar u1x[],
                               const TacsScalar e0ty[],
                               const TacsScalar Ct[],
                               const TacsScalar u0xd[],
                               const TacsScalar u1xd[],
                               const TacsScalar e0tyd[],
                               const TacsScalar Ctd[],
                               TacsScalar e[],
                               TacsScalar ed[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = e0ty[0];
    e[1] = e0ty[3];
    e[2] = 2.0*e0ty[1];

    // Compute the bending strain
    e[3] = u1x[0] + (u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6]);
    e[4] = u1x[4] + (u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7]);
    e[5] = u1x[1] + u1x[3] + (u0x[0]*u1x[1] + u0x[3]*u1x[4] + u0x[6]*u1x[7] +
                              u1x[0]*u0x[1] + u1x[3]*u0x[4] + u1x[6]*u0x[7]);

    // Add the components of the shear strain
    e[6] = 2.0*e0ty[4];
    e[7] = 2.0*e0ty[2];

    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    e[8] =
      (Ct[3]*(1.0 + u0x[0]) + Ct[4]*u0x[3] + Ct[5]*u0x[6]) -
      (Ct[0]*u0x[1] + Ct[1]*(1.0 + u0x[4]) + Ct[2]*u0x[7]);

    // Evaluate the in-plane strains from the tying strain expressions
    ed[0] = e0tyd[0];
    ed[1] = e0tyd[3];
    ed[2] = 2.0*e0tyd[1];

    // Compute the bending strain
    ed[3] = u1xd[0] + (u0xd[0]*u1x[0] + u0xd[3]*u1x[3] + u0xd[6]*u1x[6] +
                       u0x[0]*u1xd[0] + u0x[3]*u1xd[3] + u0x[6]*u1xd[6]);
    ed[4] = u1xd[4] + (u0xd[1]*u1x[1] + u0xd[4]*u1x[4] + u0xd[7]*u1x[7] +
                       u0x[1]*u1xd[1] + u0x[4]*u1xd[4] + u0x[7]*u1xd[7]);
    ed[5] = u1xd[1] + u1xd[3] + (u0xd[0]*u1x[1] + u0xd[3]*u1x[4] + u0xd[6]*u1x[7] +
                                 u1xd[0]*u0x[1] + u1xd[3]*u0x[4] + u1xd[6]*u0x[7] +
                                 u0x[0]*u1xd[1] + u0x[3]*u1xd[4] + u0x[6]*u1xd[7] +
                                 u1x[0]*u0xd[1] + u1x[3]*u0xd[4] + u1x[6]*u0xd[7]);

    // Add the components of the shear strain
    ed[6] = 2.0*e0tyd[4];
    ed[7] = 2.0*e0tyd[2];

    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    ed[8] =
      ((Ctd[3]*(1.0 + u0x[0]) + Ctd[4]*u0x[3] + Ctd[5]*u0x[6]) +
       (Ct[3]*u0xd[0] + Ct[4]*u0xd[3] + Ct[5]*u0xd[6])) -
      ((Ctd[0]*u0x[1] + Ctd[1]*(1.0 + u0x[4]) + Ctd[2]*u0x[7]) +
       (Ct[0]*u0xd[1] + Ct[1]*u0xd[4] + Ct[2]*u0xd[7]));
  }
};

template <int vars_per_node, class basis, class model>
int TacsTestShellModelDerivatives( double dh=1e-7,
                                   int test_print_level=2,
                                   double test_fail_atol=1e-5,
                                   double test_fail_rtol=1e-5 ){
  // Set the failure flag
  int fail = 0;

  // Set random values for the constitutive data and inputs
  TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  TacsScalar u0x[9], u1x[9], e0ty[6], Ct[9];
  TacsScalar detXd;

  // Set random data
  TacsGenerateRandomArray(Cs, TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES);
  TacsGenerateRandomArray(u0x, 9);
  TacsGenerateRandomArray(u1x, 9);
  TacsGenerateRandomArray(e0ty, 6);
  TacsGenerateRandomArray(Ct, 9);
  TacsGenerateRandomArray(&detXd, 1);

  // Compute the strain
  TacsScalar e[9];
  model::evalStrain(u0x, u1x, e0ty, Ct, e);

  // Compute the stress
  TacsScalar s[9];
  TacsScalar drill;
  const TacsScalar *A, *B, *D, *As;
  TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);
  TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);

  // Compute the derivative of the product of the stress and strain
  // with respect to u0x, u1x and e0ty
  TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
  model::evalStrainSens(detXd, s, u0x, u1x, Ct,
                        du0x, du1x, de0ty, dCt);

  TacsScalar f0 = 0.0;
  for ( int j = 0; j < 9; j++ ){
    f0 += 0.5*detXd*e[j]*s[j];
  }

  // Compute against the derivatives for the strain
  TacsScalar fdu0x[9];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar u0xt[9], et[9], st[9];
    memcpy(u0xt, u0x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u0xt[i] = u0x[i] + TacsScalar(0.0, dh);
#else
    u0xt[i] = u0x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0xt, u1x, e0ty, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 9; j++ ){
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
    TacsScalar u1xt[9], et[9], st[9];
    memcpy(u1xt, u1x, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    u1xt[i] = u1x[i] + TacsScalar(0.0, dh);
#else
    u1xt[i] = u1x[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1xt, e0ty, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 9; j++ ){
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
  TacsScalar fde0ty[6];
  for ( int i = 0; i < 6; i++ ){
    TacsScalar e0tyt[6], et[9], st[9];
    memcpy(e0tyt, e0ty, 6*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    e0tyt[i] = e0ty[i] + TacsScalar(0.0, dh);
#else
    e0tyt[i] = e0ty[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, e0tyt, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 9; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fde0ty[i] = TacsImagPart(f1)/dh;
#else
    fde0ty[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(de0ty, fde0ty, 6, &max_err_index);
  max_rel = TacsGetMaxRelError(de0ty, fde0ty, 6, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. e0ty\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "de0ty", de0ty, fde0ty, 6);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute against the derivatives for the strain
  TacsScalar fdCt[9];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar Ctt[9], et[9], st[9];
    memcpy(Ctt, Ct, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    Ctt[i] = Ct[i] + TacsScalar(0.0, dh);
#else
    Ctt[i] = Ct[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, e0ty, Ctt, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar f1 = 0.0;
    for ( int j = 0; j < 9; j++ ){
      f1 += 0.5*detXd*et[j]*st[j];
    }

#ifdef TACS_USE_COMPLEX
    fdCt[i] = TacsImagPart(f1)/dh;
#else
    fdCt[i] = (f1 - f0)/dh;
#endif // TACS_USE_COMPLEX
  }

  // Compute the error
  max_err = TacsGetMaxError(dCt, fdCt, 9, &max_err_index);
  max_rel = TacsGetMaxRelError(dCt, fdCt, 9, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the derivative w.r.t. Ct\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "dCt", du1x, fdu1x, 9);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
  TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
  TacsScalar d2Ct[81], d2Ctu0x[81];
  model::evalStrainHessian(detXd, s, Cs, u0x, u1x, Ct,
                           d2u0x, d2u1x, d2u0xu1x,
                           d2e0ty, d2e0tyu0x, d2e0tyu1x,
                           d2Ct, d2Ctu0x);

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
    model::evalStrain(u0xt, u1x, e0ty, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6], dCtt[9];
    model::evalStrainSens(detXd, st, u0xt, u1x, Ct, du0xt, du1xt, de0tyt, dCtt);

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
    model::evalStrain(u0x, u1xt, e0ty, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6], dCtt[9];
    model::evalStrainSens(detXd, st, u0x, u1xt, Ct, du0xt, du1xt, de0tyt, dCtt);

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
    model::evalStrain(u0x, u1x, e0tyt, Ct, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6], dCtt[9];
    model::evalStrainSens(detXd, st, u0x, u1x, Ct, du0xt, du1xt, de0tyt, dCtt);

    for ( int j = 0; j < 6; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2e0ty[6*i + j] = TacsImagPart(de0tyt[j])/dh;
#else
      fd2e0ty[6*i + j] = (de0tyt[j] - de0ty[j])/dh;
#endif // TACS_USE_COMPLEX
    }

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2e0tyu0x[9*i + j] = TacsImagPart(du0x[j])/dh;
      fd2e0tyu1x[9*i + j] = TacsImagPart(du1x[j])/dh;
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

  TacsScalar fd2Ct[81], fd2Ctu0x[81];
  for ( int i = 0; i < 9; i++ ){
    TacsScalar Ctt[9], et[9], st[9];
    memcpy(Ctt, Ct, 9*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    Ctt[i] = Ct[i] + TacsScalar(0.0, dh);
#else
    Ctt[i] = Ct[i] + dh;
#endif // TACS_USE_COMPLEX
    model::evalStrain(u0x, u1x, e0ty, Ctt, et);
    TACSShellConstitutive::computeStress(A, B, D, As, drill, et, st);

    TacsScalar du0xt[9], du1xt[9], de0tyt[6], dCtt[9];
    model::evalStrainSens(detXd, st, u0x, u1x, Ctt, du0xt, du1xt, de0tyt, dCtt);

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2Ct[9*i + j] = TacsImagPart(dCtt[j])/dh;
#else
      fd2Ct[9*i + j] = (dCtt[j] - dCt[j])/dh;
#endif // TACS_USE_COMPLEX
    }

    for ( int j = 0; j < 9; j++ ){
#ifdef TACS_USE_COMPLEX
      fd2Ctu0x[9*i + j] = TacsImagPart(du0x[j])/dh;
#else
      fd2Ctu0x[9*i + j] = (du0xt[j] - du0x[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2Ct, fd2Ct, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2Ct, fd2Ct, 81, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. Ct\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2Ct", d2Ct, fd2Ct, 81);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  // Compute the error
  max_err = TacsGetMaxError(d2Ctu0x, fd2Ctu0x, 81, &max_err_index);
  max_rel = TacsGetMaxRelError(d2Ctu0x, fd2Ctu0x, 81, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. Ct and u0x\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2Ctu0x", d2Ctu0x, fd2Ctu0x, 81);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  return fail;
}

#endif // TACS_SHELL_ELEMENT_MODEL_H
