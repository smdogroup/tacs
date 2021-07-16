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
  template <int vars_per_node, class basis>
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

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
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, fn, n0);

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

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          dUxi[0] = dety[index]*Xxi[0];
          dUxi[1] = 0.0;
          dUxi[2] = dety[index]*Xxi[2];
          dUxi[3] = 0.0;
          dUxi[4] = dety[index]*Xxi[4];
          dUxi[5] = 0.0;
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          dUxi[0] = 0.0;
          dUxi[1] = dety[index]*Xxi[1];
          dUxi[2] = 0.0;
          dUxi[3] = dety[index]*Xxi[3];
          dUxi[4] = 0.0;
          dUxi[5] = dety[index]*Xxi[5];
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          dUxi[0] = 0.5*dety[index]*Xxi[1];
          dUxi[1] = 0.5*dety[index]*Xxi[0];
          dUxi[2] = 0.5*dety[index]*Xxi[3];
          dUxi[3] = 0.5*dety[index]*Xxi[2];
          dUxi[4] = 0.5*dety[index]*Xxi[5];
          dUxi[5] = 0.5*dety[index]*Xxi[4];
        }
        else {
          TacsScalar d0[3], dd0[3], n0[3];
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, fn, n0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi[0] = 0.0;
            dUxi[1] = 0.5*dety[index]*n0[0];
            dUxi[2] = 0.0;
            dUxi[3] = 0.5*dety[index]*n0[1];
            dUxi[4] = 0.0;
            dUxi[5] = 0.5*dety[index]*n0[2];

            dd0[0] = 0.5*dety[index]*Xxi[1];
            dd0[1] = 0.5*dety[index]*Xxi[3];
            dd0[2] = 0.5*dety[index]*Xxi[5];
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi[0] = 0.5*dety[index]*n0[0];
            dUxi[1] = 0.0;
            dUxi[2] = 0.5*dety[index]*n0[1];
            dUxi[3] = 0.0;
            dUxi[4] = 0.5*dety[index]*n0[2];
            dUxi[5] = 0.0;

            dd0[0] = 0.5*dety[index]*Xxi[0];
            dd0[1] = 0.5*dety[index]*Xxi[2];
            dd0[2] = 0.5*dety[index]*Xxi[4];
          }

          basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
        }

        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi, res);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian( const TacsScalar Xpts[],
                                            const TacsScalar fn[],
                                            const TacsScalar vars[],
                                            const TacsScalar d[],
                                            const TacsScalar d2ety[],
                                            TacsScalar mat[],
                                            TacsScalar d2d[],
                                            TacsScalar d2du[] ){
    // Set the number of tying fields
    const int num_tying_fields = 5;

    // Initialize the data
    TacsScalar n0ty[3*basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6*basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity;

    // Pre-compute terms needed at each tying point
    for ( int f1 = 0; f1 < num_tying_fields; f1++ ){
      const int nty1 = basis::getNumTyingPoints(f1);

      for ( int ty1 = 0; ty1 < nty1; ty1++, n0 += 3, Xxi += 6 ){
        double pt[2];
        basis::getTyingPoint(f1, ty1, pt);

        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFields<3, 3>(pt, fn, n0);
      }
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity;
    for ( int f1 = 0, base = 0; f1 < num_tying_fields;
      base += basis::NUM_TYING_POINTS*basis::getNumTyingPoints(f1), f1++ ){
      const int nty1 = basis::getNumTyingPoints(f1);

      for ( int ty1 = 0; ty1 < nty1; ty1++, n01 += 3, Xxi1 += 6 ){

        TacsScalar du2[3*basis::NUM_NODES];
        TacsScalar dd2[3*basis::NUM_NODES];
        memset(du2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
        memset(dd2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

        TacsScalar *n02 = n0ty, *Xxi2 = Xxity;
        for ( int f2 = 0, offset = base; f2 < num_tying_fields;
          offset += nty1*basis::getNumTyingPoints(f2), f2++ ){
          const int nty2 = basis::getNumTyingPoints(f2);

          for ( int ty2 = 0; ty2 < nty2; ty2++, n02 += 3, Xxi2 += 6 ){
            double pt2[2];
            basis::getTyingPoint(f2, ty2, pt2);

            TacsScalar value = d2ety[offset + ty2 + ty1*nty2];

            TacsScalar dUxi2[6];
            if (f2 == 0){
              // Compute g11 = e1^{T}*G*e1
              dUxi2[0] = value*Xxi2[0];
              dUxi2[1] = 0.0;
              dUxi2[2] = value*Xxi2[2];
              dUxi2[3] = 0.0;
              dUxi2[4] = value*Xxi2[4];
              dUxi2[5] = 0.0;
            }
            else if (f2 == 1){
              // Compute g22 = e2^{T}*G*e2
              dUxi2[0] = 0.0;
              dUxi2[1] = value*Xxi2[1];
              dUxi2[2] = 0.0;
              dUxi2[3] = value*Xxi2[3];
              dUxi2[4] = 0.0;
              dUxi2[5] = value*Xxi2[5];
            }
            else if (f2 == 2){
              // Compute g12 = e2^{T}*G*e1
              dUxi2[0] = 0.5*value*Xxi2[1];
              dUxi2[1] = 0.5*value*Xxi2[0];
              dUxi2[2] = 0.5*value*Xxi2[3];
              dUxi2[3] = 0.5*value*Xxi2[2];
              dUxi2[4] = 0.5*value*Xxi2[5];
              dUxi2[5] = 0.5*value*Xxi2[4];
            }
            else {
              TacsScalar dd02[3];
              if (f2 == 3){
                // Compute g23 = e2^{T}*G*e3
                dUxi2[0] = 0.0;
                dUxi2[1] = 0.5*value*n02[0];
                dUxi2[2] = 0.0;
                dUxi2[3] = 0.5*value*n02[1];
                dUxi2[4] = 0.0;
                dUxi2[5] = 0.5*value*n02[2];

                dd02[0] = 0.5*value*Xxi2[1];
                dd02[1] = 0.5*value*Xxi2[3];
                dd02[2] = 0.5*value*Xxi2[5];
              }
              else if (f2 == 4){
                // Compute g13 = e1^{T}*G*e3
                dUxi2[0] = 0.5*value*n02[0];
                dUxi2[1] = 0.0;
                dUxi2[2] = 0.5*value*n02[1];
                dUxi2[3] = 0.0;
                dUxi2[4] = 0.5*value*n02[2];
                dUxi2[5] = 0.0;

                dd02[0] = 0.5*value*Xxi2[0];
                dd02[1] = 0.5*value*Xxi2[2];
                dd02[2] = 0.5*value*Xxi2[4];
              }

              basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
            }

            basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
          }
        }

        TacsScalar du1[3*basis::NUM_NODES];
        memset(du1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

        double pt1[2];
        basis::getTyingPoint(f1, ty1, pt1);

        // Store the the derivative information for the first point
        TacsScalar dUxi1[6];

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
          TacsScalar dd1[3*basis::NUM_NODES];
          memset(dd1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

          TacsScalar dd01[3];
          if (f1 == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi1[0] = 0.0;
            dUxi1[1] = 0.5*n01[0];
            dUxi1[2] = 0.0;
            dUxi1[3] = 0.5*n01[1];
            dUxi1[4] = 0.0;
            dUxi1[5] = 0.5*n01[2];

            dd01[0] = 0.5*Xxi1[1];
            dd01[1] = 0.5*Xxi1[3];
            dd01[2] = 0.5*Xxi1[5];
          }
          else if (f1 == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi1[0] = 0.5*n01[0];
            dUxi1[1] = 0.0;
            dUxi1[2] = 0.5*n01[1];
            dUxi1[3] = 0.0;
            dUxi1[4] = 0.5*n01[2];
            dUxi1[5] = 0.0;

            dd01[0] = 0.5*Xxi1[0];
            dd01[1] = 0.5*Xxi1[2];
            dd01[2] = 0.5*Xxi1[4];
          }

          basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);

          for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
            for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
              d2d[3*basis::NUM_NODES*i + j] += dd1[i]*dd2[j];
            }
          }

          for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
            for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
              d2du[3*basis::NUM_NODES*i + j] += dd1[i]*du2[j];
            }
          }
        }

        basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

        const int nvars = vars_per_node*basis::NUM_NODES;
        for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
          int ii = vars_per_node*(i / 3) + i % 3;
          for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
            int jj = vars_per_node*(j / 3) + j % 3;
            mat[nvars*ii + jj] += du1[i]*du2[j];
          }
        }
      }
    }
  }

  /*
    Compute the directional derivative
  */
  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv( const TacsScalar Xpts[],
                                       const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, Uxid);

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
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, dd, d0d);
          basis::template interpFields<3, 3>(pt, fn, n0);

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
  }

  /**
    Evaluate the derivative of the strain
  */
  static void evalStrainSens( const TacsScalar scale,
                              const TacsScalar dfde[],
                              const TacsScalar u0x[],
                              const TacsScalar u1x[],
                              TacsScalar du0x[],
                              TacsScalar du1x[],
                              TacsScalar de0ty[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = scale*dfde[0];
    de0ty[1] = 2.0*scale*dfde[2];
    de0ty[2] = 2.0*scale*dfde[7];
    de0ty[3] = scale*dfde[1];
    de0ty[4] = 2.0*scale*dfde[6];
    de0ty[5] = 0.0;

    du0x[0] = 0.0;
    du0x[1] = 0.0;
    du0x[2] = 0.0;
    du0x[3] = 0.0;
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
                                 const TacsScalar e0ty[],
                                 TacsScalar d2u0x[],
                                 TacsScalar d2u1x[],
                                 TacsScalar d2u0xu1x[],
                                 TacsScalar d2e0ty[],
                                 TacsScalar d2e0tyu0x[],
                                 TacsScalar d2e0tyu1x[] ){
    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    memset(d2u0x, 0, 81*sizeof(TacsScalar));
    memset(d2u1x, 0, 81*sizeof(TacsScalar));
    memset(d2u0xu1x, 0, 81*sizeof(TacsScalar));
    memset(d2e0ty, 0, 36*sizeof(TacsScalar));
    memset(d2e0tyu0x, 0, 54*sizeof(TacsScalar));
    memset(d2e0tyu1x, 0, 54*sizeof(TacsScalar));

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
  }

  static TacsScalar evalDrillStrain( const TacsScalar u0x[],
                                     const TacsScalar Ct[] ){
    // Compute the rotational penalty
    return 0.5*(Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  }

  static void evalDrillStrainSens( TacsScalar scale,
                                   const TacsScalar u0x[],
                                   const TacsScalar Ct[],
                                   TacsScalar du0x[],
                                   TacsScalar dCt[] ){
    dCt[0] = 0.0;
    dCt[1] = -0.5*scale;
    dCt[2] = 0.0;
    dCt[3] = 0.5*scale;
    dCt[4] = 0.0;
    dCt[5] = 0.0;
    dCt[6] = 0.0;
    dCt[7] = 0.0;
    dCt[8] = 0.0;

    du0x[0] = 0.0;
    du0x[1] = -0.5*scale;
    du0x[2] = 0.0;
    du0x[3] = 0.5*scale;
    du0x[4] = 0.0;
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;
  }

  // static void evalDrillStrainHessian( TacsScalar det,
  //                                     TacsScalar d2et,
  //                                     const TacsScalar u0x[],
  //                                     const TacsScalar Ct[],
  //                                     TacsScalar d2u0x[],
  //                                     TacsScalar d2Ct[],
  //                                     TacsScalar d2Ctu0x[] ){
  //   memset(d2u0x, 0, 81*sizeof(TacsScalar));
  //   memset(d2Ct, 0, 81*sizeof(TacsScalar));
  //   memset(d2Ctu0x, 0, 81*sizeof(TacsScalar));

  //   TacsScalar *d2;

  //   // Compute the contribution from the drilling strain
  //   // e[8] = 0.5*(Ct[3] + u0x[3] - Ct[1] - u0x[1]);
  //   d2 = &d2Ct[3*9];
  //   d2[3] = 0.25*scale*drill;
  //   d2[1] = -0.25*scale*drill;

  //   d2 = &d2Ct[9];
  //   d2[3] = -0.25*scale*drill;
  //   d2[1] = 0.25*scale*drill;

  //   d2 = &d2u0x[3*9];
  //   d2[3] = 0.25*scale*drill;
  //   d2[1] = -0.25*scale*drill;

  //   d2 = &d2u0x[9];
  //   d2[3] = -0.25*scale*drill;
  //   d2[1] = 0.25*scale*drill;

  //   d2 = &d2Ctu0x[3*9];
  //   d2[3] = 0.25*scale*drill;
  //   d2[1] = -0.25*scale*drill;

  //   d2 = &d2Ctu0x[9];
  //   d2[3] = -0.25*scale*drill;
  //   d2[1] = 0.25*scale*drill;
  // }





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
  template <int vars_per_node, class basis>
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

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
          TacsScalar d0[3], n0[3];
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, fn, n0);

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

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainTranspose( const TacsScalar Xpts[],
                                              const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);

        if (field == 0){
          // Compute g11 = e1^{T}*G*e1
          dUxi[0] = dety[index]*(Xxi[0] + Uxi[0]);
          dUxi[1] = 0.0;
          dUxi[2] = dety[index]*(Xxi[2] + Uxi[2]);
          dUxi[3] = 0.0;
          dUxi[4] = dety[index]*(Xxi[4] + Uxi[4]);
          dUxi[5] = 0.0;
        }
        else if (field == 1){
          // Compute g22 = e2^{T}*G*e2
          dUxi[0] = 0.0;
          dUxi[1] = dety[index]*(Xxi[1] + Uxi[1]);
          dUxi[2] = 0.0;
          dUxi[3] = dety[index]*(Xxi[3] + Uxi[3]);
          dUxi[4] = 0.0;
          dUxi[5] = dety[index]*(Xxi[5] + Uxi[5]);
        }
        else if (field == 2){
          // Compute g12 = e2^{T}*G*e1
          dUxi[0] = 0.5*dety[index]*(Xxi[1] + Uxi[1]);
          dUxi[1] = 0.5*dety[index]*(Xxi[0] + Uxi[0]);
          dUxi[2] = 0.5*dety[index]*(Xxi[3] + Uxi[3]);
          dUxi[3] = 0.5*dety[index]*(Xxi[2] + Uxi[2]);
          dUxi[4] = 0.5*dety[index]*(Xxi[5] + Uxi[5]);
          dUxi[5] = 0.5*dety[index]*(Xxi[4] + Uxi[4]);
        }
        else {
          TacsScalar n0[3], d0[3], dd0[3];
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, fn, n0);

          if (field == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi[0] = 0.0;
            dUxi[1] = 0.5*dety[index]*(n0[0] + d0[0]);
            dUxi[2] = 0.0;
            dUxi[3] = 0.5*dety[index]*(n0[1] + d0[1]);
            dUxi[4] = 0.0;
            dUxi[5] = 0.5*dety[index]*(n0[2] + d0[2]);

            dd0[0] = 0.5*dety[index]*(Xxi[1] + Uxi[1]);
            dd0[1] = 0.5*dety[index]*(Xxi[3] + Uxi[3]);
            dd0[2] = 0.5*dety[index]*(Xxi[5] + Uxi[5]);
          }
          else if (field == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi[0] = 0.5*dety[index]*(n0[0] + d0[0]);
            dUxi[1] = 0.0;
            dUxi[2] = 0.5*dety[index]*(n0[1] + d0[1]);
            dUxi[3] = 0.0;
            dUxi[4] = 0.5*dety[index]*(n0[2] + d0[2]);
            dUxi[5] = 0.0;

            dd0[0] = 0.5*dety[index]*(Xxi[0] + Uxi[0]);
            dd0[1] = 0.5*dety[index]*(Xxi[2] + Uxi[2]);
            dd0[2] = 0.5*dety[index]*(Xxi[4] + Uxi[4]);
          }
          basis::template addInterpFieldsTranspose<3, 3>(pt, dd0, dd);
        }

        basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi, res);
      }
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian( const TacsScalar Xpts[],
                                            const TacsScalar fn[],
                                            const TacsScalar vars[],
                                            const TacsScalar d[],
                                            const TacsScalar d2ety[],
                                            TacsScalar mat[],
                                            TacsScalar d2d[],
                                            TacsScalar d2du[] ){
    // Set the number of tying fields
    const int num_tying_fields = 5;

    // Initialize the data
    TacsScalar n0ty[3*basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6*basis::NUM_TYING_POINTS];
    TacsScalar d0ty[3*basis::NUM_TYING_POINTS];
    TacsScalar Uxity[6*basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity, *d0 = d0ty, *Uxi = Uxity;

    // Pre-compute terms needed at each tying point
    for ( int f1 = 0; f1 < num_tying_fields; f1++ ){
      const int nty1 = basis::getNumTyingPoints(f1);

      for ( int ty1 = 0; ty1 < nty1; ty1++, n0 += 3, Xxi += 6, d0 += 3, Uxi += 6 ){
        double pt[2];
        basis::getTyingPoint(f1, ty1, pt);

        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFields<3, 3>(pt, fn, n0);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
        basis::template interpFields<3, 3>(pt, d, d0);
      }
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity, *d01 = d0ty, *Uxi1 = Uxity;
    for ( int f1 = 0, base = 0; f1 < num_tying_fields;
      base += basis::NUM_TYING_POINTS*basis::getNumTyingPoints(f1), f1++ ){
      const int nty1 = basis::getNumTyingPoints(f1);

      for ( int ty1 = 0; ty1 < nty1; ty1++, n01 += 3, Xxi1 += 6, d01 += 3, Uxi1 += 6 ){

        TacsScalar du2[3*basis::NUM_NODES];
        TacsScalar dd2[3*basis::NUM_NODES];
        memset(du2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
        memset(dd2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

        TacsScalar *n02 = n0ty, *Xxi2 = Xxity, *d02 = d0ty, *Uxi2 = Uxity;
        for ( int f2 = 0, offset = base; f2 < num_tying_fields;
          offset += nty1*basis::getNumTyingPoints(f2), f2++ ){
          const int nty2 = basis::getNumTyingPoints(f2);

          for ( int ty2 = 0; ty2 < nty2; ty2++, n02 += 3, Xxi2 += 6, d02 += 3, Uxi2 += 6 ){
            double pt2[2];
            basis::getTyingPoint(f2, ty2, pt2);

            TacsScalar value = d2ety[offset + ty2 + ty1*nty2];

            TacsScalar dUxi2[6];
            if (f2 == 0){
              // Compute g11 = e1^{T}*G*e1
              dUxi2[0] = value*(Xxi2[0] + Uxi2[0]);
              dUxi2[1] = 0.0;
              dUxi2[2] = value*(Xxi2[2] + Uxi2[2]);
              dUxi2[3] = 0.0;
              dUxi2[4] = value*(Xxi2[4] + Uxi2[4]);
              dUxi2[5] = 0.0;
            }
            else if (f2 == 1){
              // Compute g22 = e2^{T}*G*e2
              dUxi2[0] = 0.0;
              dUxi2[1] = value*(Xxi2[1] + Uxi2[1]);
              dUxi2[2] = 0.0;
              dUxi2[3] = value*(Xxi2[3] + Uxi2[3]);
              dUxi2[4] = 0.0;
              dUxi2[5] = value*(Xxi2[5] + Uxi2[5]);
            }
            else if (f2 == 2){
              // Compute g12 = e2^{T}*G*e1
              dUxi2[0] = 0.5*value*(Xxi2[1] + Uxi2[1]);
              dUxi2[1] = 0.5*value*(Xxi2[0] + Uxi2[0]);
              dUxi2[2] = 0.5*value*(Xxi2[3] + Uxi2[3]);
              dUxi2[3] = 0.5*value*(Xxi2[2] + Uxi2[2]);
              dUxi2[4] = 0.5*value*(Xxi2[5] + Uxi2[5]);
              dUxi2[5] = 0.5*value*(Xxi2[4] + Uxi2[4]);
            }
            else {
              TacsScalar dd02[3];
              if (f2 == 3){
                // Compute g23 = e2^{T}*G*e3
                dUxi2[0] = 0.0;
                dUxi2[1] = 0.5*value*(n02[0] + d02[0]);
                dUxi2[2] = 0.0;
                dUxi2[3] = 0.5*value*(n02[1] + d02[1]);
                dUxi2[4] = 0.0;
                dUxi2[5] = 0.5*value*(n02[2] + d02[2]);

                dd02[0] = 0.5*value*(Xxi2[1] + Uxi2[1]);
                dd02[1] = 0.5*value*(Xxi2[3] + Uxi2[3]);
                dd02[2] = 0.5*value*(Xxi2[5] + Uxi2[5]);
              }
              else if (f2 == 4){
                // Compute g13 = e1^{T}*G*e3
                dUxi2[0] = 0.5*value*(n02[0] + d02[0]);
                dUxi2[1] = 0.0;
                dUxi2[2] = 0.5*value*(n02[1] + d02[1]);
                dUxi2[3] = 0.0;
                dUxi2[4] = 0.5*value*(n02[2] + d02[2]);
                dUxi2[5] = 0.0;

                dd02[0] = 0.5*value*(Xxi2[0] + Uxi2[0]);
                dd02[1] = 0.5*value*(Xxi2[2] + Uxi2[2]);
                dd02[2] = 0.5*value*(Xxi2[4] + Uxi2[4]);
              }

              basis::template addInterpFieldsTranspose<3, 3>(pt2, dd02, dd2);
            }

            basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
          }
        }

        TacsScalar du1[3*basis::NUM_NODES];
        memset(du1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

        double pt1[2];
        basis::getTyingPoint(f1, ty1, pt1);

        // Store the the derivative information for the first point
        TacsScalar dUxi1[6];

        if (f1 == 0){
          // Compute g11 = e1^{T}*G*e1
          dUxi1[0] = (Xxi1[0] + Uxi1[0]);
          dUxi1[1] = 0.0;
          dUxi1[2] = (Xxi1[2] + Uxi1[2]);
          dUxi1[3] = 0.0;
          dUxi1[4] = (Xxi1[4] + Uxi1[4]);
          dUxi1[5] = 0.0;
        }
        else if (f1 == 1){
          // Compute g22 = e2^{T}*G*e2
          dUxi1[0] = 0.0;
          dUxi1[1] = (Xxi1[1] + Uxi1[1]);
          dUxi1[2] = 0.0;
          dUxi1[3] = (Xxi1[3] + Uxi1[3]);
          dUxi1[4] = 0.0;
          dUxi1[5] = (Xxi1[5] + Uxi1[5]);
        }
        else if (f1 == 2){
          // Compute g12 = e2^{T}*G*e1
          dUxi1[0] = 0.5*(Xxi1[1] + Uxi1[1]);
          dUxi1[1] = 0.5*(Xxi1[0] + Uxi1[0]);
          dUxi1[2] = 0.5*(Xxi1[3] + Uxi1[3]);
          dUxi1[3] = 0.5*(Xxi1[2] + Uxi1[2]);
          dUxi1[4] = 0.5*(Xxi1[5] + Uxi1[5]);
          dUxi1[5] = 0.5*(Xxi1[4] + Uxi1[4]);
        }
        else {
          TacsScalar dd1[3*basis::NUM_NODES];
          memset(dd1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

          TacsScalar dd01[3];
          if (f1 == 3){
            // Compute g23 = e2^{T}*G*e3
            dUxi1[0] = 0.0;
            dUxi1[1] = 0.5*(n01[0] + d01[0]);
            dUxi1[2] = 0.0;
            dUxi1[3] = 0.5*(n01[1] + d01[1]);
            dUxi1[4] = 0.0;
            dUxi1[5] = 0.5*(n01[2] + d01[2]);

            dd01[0] = 0.5*(Xxi1[1] + Uxi1[1]);
            dd01[1] = 0.5*(Xxi1[3] + Uxi1[3]);
            dd01[2] = 0.5*(Xxi1[5] + Uxi1[5]);
          }
          else if (f1 == 4){
            // Compute g13 = e1^{T}*G*e3
            dUxi1[0] = 0.5*(n01[0] + d01[0]);
            dUxi1[1] = 0.0;
            dUxi1[2] = 0.5*(n01[1] + d01[1]);
            dUxi1[3] = 0.0;
            dUxi1[4] = 0.5*(n01[2] + d01[2]);
            dUxi1[5] = 0.0;

            dd01[0] = 0.5*(Xxi1[0] + Uxi1[0]);
            dd01[1] = 0.5*(Xxi1[2] + Uxi1[2]);
            dd01[2] = 0.5*(Xxi1[4] + Uxi1[4]);
          }

          basis::template addInterpFieldsTranspose<3, 3>(pt1, dd01, dd1);

          for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
            for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
              d2d[3*basis::NUM_NODES*i + j] += dd1[i]*dd2[j];
            }
          }

          for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
            for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
              d2du[3*basis::NUM_NODES*i + j] += dd1[i]*du2[j];
            }
          }
        }

        basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

        const int nvars = vars_per_node*basis::NUM_NODES;
        for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
          int ii = vars_per_node*(i / 3) + i % 3;
          for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
            int jj = vars_per_node*(j / 3) + j % 3;
            mat[nvars*ii + jj] += du1[i]*du2[j];
          }
        }
      }
    }
  }

  template <int vars_per_node, class basis>
  static void computeTyingStrainDeriv( const TacsScalar Xpts[],
                                       const TacsScalar fn[],
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
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, Uxi);
        basis::template interpFieldsGrad<vars_per_node, 3>(pt, varsd, Uxid);

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
          basis::template interpFields<3, 3>(pt, d, d0);
          basis::template interpFields<3, 3>(pt, dd, d0d);
          basis::template interpFields<3, 3>(pt, fn, n0);

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
  }

  /**
    Evaluate the derivative of the strain
  */
  static void evalStrainSens( const TacsScalar scale,
                              const TacsScalar dfde[],
                              const TacsScalar u0x[],
                              const TacsScalar u1x[],
                              TacsScalar du0x[],
                              TacsScalar du1x[],
                              TacsScalar de0ty[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = scale*dfde[0];
    de0ty[1] = 2.0*scale*dfde[2];
    de0ty[2] = 2.0*scale*dfde[7];
    de0ty[3] = scale*dfde[1];
    de0ty[4] = 2.0*scale*dfde[6];
    de0ty[5] = 0.0;

    // Derivative with respect to u0x
    du0x[0] = scale*(dfde[3]*u1x[0] + dfde[5]*u1x[1]);
    du0x[1] = scale*(dfde[4]*u1x[1] + dfde[5]*u1x[0]);
    du0x[2] = 0.0;
    du0x[3] = scale*(dfde[3]*u1x[3] + dfde[5]*u1x[4]);
    du0x[4] = scale*(dfde[4]*u1x[4] + dfde[5]*u1x[3]);
    du0x[5] = 0.0;
    du0x[6] = scale*(dfde[3]*u1x[6] + dfde[5]*u1x[7]);
    du0x[7] = scale*(dfde[4]*u1x[7] + dfde[5]*u1x[6]);
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
  }

  static TacsScalar evalDrillStrain( const TacsScalar u0x[],
                                     const TacsScalar Ct[] ){
    // e2^{T}*Ct*(e1 + u_{,x}*e1) - e1^{T}*Ct*(e2 + u_{,x}*e2)
    return ((Ct[3]*(1.0 + u0x[0]) + Ct[4]*u0x[3] + Ct[5]*u0x[6]) -
            (Ct[0]*u0x[1] + Ct[1]*(1.0 + u0x[4]) + Ct[2]*u0x[7]));
  }


  static void evalDrillStrainSens( TacsScalar scale,
                                   const TacsScalar u0x[],
                                   const TacsScalar Ct[],
                                   TacsScalar du0x[],
                                   TacsScalar dCt[] ){
    // Derivative with respect to u0x
    du0x[0] = scale*Ct[3];
    du0x[1] = -scale*Ct[0];
    du0x[2] = 0.0;
    du0x[3] = scale*Ct[4];
    du0x[4] = -scale*Ct[1];
    du0x[5] = 0.0;
    du0x[6] = scale*Ct[5];
    du0x[7] = -scale*Ct[2];
    du0x[8] = 0.0;

    dCt[0] = -scale*u0x[1];
    dCt[1] = -scale*(1.0 + u0x[4]);
    dCt[2] = -scale*u0x[7];
    dCt[3] = scale*(1.0 + u0x[0]);
    dCt[4] = scale*u0x[3];
    dCt[5] = scale*u0x[6];
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

  static void evalStrainHessian( const TacsScalar scale,
                                 const TacsScalar dfde[],
                                 const TacsScalar Cs[],
                                 const TacsScalar u0x[],
                                 const TacsScalar u1x[],
                                 const TacsScalar e0ty[],
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


    // Compute the second derivatives
    memset(d2e0ty, 0, 36*sizeof(TacsScalar));

    TacsScalar *d2 = &d2e0ty[0];
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

    d2 = &d2e0ty[4*6];
    d2[4] = 4.0*scale*As[0];
    d2[2] = 4.0*scale*As[1];

    d2 = &d2e0ty[2*6];
    d2[4] = 4.0*scale*As[1];
    d2[2] = 4.0*scale*As[2];

    d2u0x[0] = scale*(Ct[3]*Ct[3]*drill + u1x[0]*(D[0]*u1x[0] + D[2]*u1x[1]) + u1x[1]*(D[2]*u1x[0] + D[5]*u1x[1]));
    d2u0x[1] = scale*(-Ct[0]*Ct[3]*drill + D[1]*u1x[0]*u1x[1] + D[2]*u1x[0]*u1x[0] + D[4]*u1x[1]*u1x[1] + D[5]*u1x[0]*u1x[1]);
    d2u0x[2] = 0.0;
    d2u0x[3] = scale*(Ct[3]*Ct[4]*drill + D[0]*u1x[0]*u1x[3] + D[2]*u1x[0]*u1x[4] + D[2]*u1x[1]*u1x[3] + D[5]*u1x[1]*u1x[4]);
    d2u0x[4] = scale*(-Ct[1]*Ct[3]*drill + D[1]*u1x[0]*u1x[4] + D[2]*u1x[0]*u1x[3] + D[4]*u1x[1]*u1x[4] + D[5]*u1x[1]*u1x[3]);
    d2u0x[5] = 0.0;
    d2u0x[6] = scale*(Ct[3]*Ct[5]*drill + D[0]*u1x[0]*u1x[6] + D[2]*u1x[0]*u1x[7] + D[2]*u1x[1]*u1x[6] + D[5]*u1x[1]*u1x[7]);
    d2u0x[7] = scale*(-Ct[2]*Ct[3]*drill + D[1]*u1x[0]*u1x[7] + D[2]*u1x[0]*u1x[6] + D[4]*u1x[1]*u1x[7] + D[5]*u1x[1]*u1x[6]);
    d2u0x[8] = 0.0;

    d2u0x[9] = scale*(-Ct[0]*Ct[3]*drill + D[1]*u1x[0]*u1x[1] + D[2]*u1x[0]*u1x[0] + D[4]*u1x[1]*u1x[1] + D[5]*u1x[0]*u1x[1]);
    d2u0x[10] = scale*(Ct[0]*Ct[0]*drill + u1x[0]*(D[4]*u1x[1] + D[5]*u1x[0]) + u1x[1]*(D[3]*u1x[1] + D[4]*u1x[0]));
    d2u0x[11] = 0.0;
    d2u0x[12] = scale*(-Ct[0]*Ct[4]*drill + D[1]*u1x[1]*u1x[3] + D[2]*u1x[0]*u1x[3] + D[4]*u1x[1]*u1x[4] + D[5]*u1x[0]*u1x[4]);
    d2u0x[13] = scale*(Ct[0]*Ct[1]*drill + D[3]*u1x[1]*u1x[4] + D[4]*u1x[0]*u1x[4] + D[4]*u1x[1]*u1x[3] + D[5]*u1x[0]*u1x[3]);
    d2u0x[14] = 0.0;
    d2u0x[15] = scale*(-Ct[0]*Ct[5]*drill + D[1]*u1x[1]*u1x[6] + D[2]*u1x[0]*u1x[6] + D[4]*u1x[1]*u1x[7] + D[5]*u1x[0]*u1x[7]);
    d2u0x[16] = scale*(Ct[0]*Ct[2]*drill + D[3]*u1x[1]*u1x[7] + D[4]*u1x[0]*u1x[7] + D[4]*u1x[1]*u1x[6] + D[5]*u1x[0]*u1x[6]);
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

    d2u0x[27] = scale*(Ct[3]*Ct[4]*drill + D[0]*u1x[0]*u1x[3] + D[2]*u1x[0]*u1x[4] + D[2]*u1x[1]*u1x[3] + D[5]*u1x[1]*u1x[4]);
    d2u0x[28] = scale*(-Ct[0]*Ct[4]*drill + D[1]*u1x[1]*u1x[3] + D[2]*u1x[0]*u1x[3] + D[4]*u1x[1]*u1x[4] + D[5]*u1x[0]*u1x[4]);
    d2u0x[29] = 0.0;
    d2u0x[30] = scale*(Ct[4]*Ct[4]*drill + u1x[3]*(D[0]*u1x[3] + D[2]*u1x[4]) + u1x[4]*(D[2]*u1x[3] + D[5]*u1x[4]));
    d2u0x[31] = scale*(-Ct[1]*Ct[4]*drill + D[1]*u1x[3]*u1x[4] + D[2]*u1x[3]*u1x[3] + D[4]*u1x[4]*u1x[4] + D[5]*u1x[3]*u1x[4]);
    d2u0x[32] = 0.0;
    d2u0x[33] = scale*(Ct[4]*Ct[5]*drill + D[0]*u1x[3]*u1x[6] + D[2]*u1x[3]*u1x[7] + D[2]*u1x[4]*u1x[6] + D[5]*u1x[4]*u1x[7]);
    d2u0x[34] = scale*(-Ct[2]*Ct[4]*drill + D[1]*u1x[3]*u1x[7] + D[2]*u1x[3]*u1x[6] + D[4]*u1x[4]*u1x[7] + D[5]*u1x[4]*u1x[6]);
    d2u0x[35] = 0.0;

    d2u0x[36] = scale*(-Ct[1]*Ct[3]*drill + D[1]*u1x[0]*u1x[4] + D[2]*u1x[0]*u1x[3] + D[4]*u1x[1]*u1x[4] + D[5]*u1x[1]*u1x[3]);
    d2u0x[37] = scale*(Ct[0]*Ct[1]*drill + D[3]*u1x[1]*u1x[4] + D[4]*u1x[0]*u1x[4] + D[4]*u1x[1]*u1x[3] + D[5]*u1x[0]*u1x[3]);
    d2u0x[38] = 0.0;
    d2u0x[39] = scale*(-Ct[1]*Ct[4]*drill + D[1]*u1x[3]*u1x[4] + D[2]*u1x[3]*u1x[3] + D[4]*u1x[4]*u1x[4] + D[5]*u1x[3]*u1x[4]);
    d2u0x[40] = scale*(Ct[1]*Ct[1]*drill + u1x[3]*(D[4]*u1x[4] + D[5]*u1x[3]) + u1x[4]*(D[3]*u1x[4] + D[4]*u1x[3]));
    d2u0x[41] = 0.0;
    d2u0x[42] = scale*(-Ct[1]*Ct[5]*drill + D[1]*u1x[4]*u1x[6] + D[2]*u1x[3]*u1x[6] + D[4]*u1x[4]*u1x[7] + D[5]*u1x[3]*u1x[7]);
    d2u0x[43] = scale*(Ct[1]*Ct[2]*drill + D[3]*u1x[4]*u1x[7] + D[4]*u1x[3]*u1x[7] + D[4]*u1x[4]*u1x[6] + D[5]*u1x[3]*u1x[6]);
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

    d2u0x[54] = scale*(Ct[3]*Ct[5]*drill + D[0]*u1x[0]*u1x[6] + D[2]*u1x[0]*u1x[7] + D[2]*u1x[1]*u1x[6] + D[5]*u1x[1]*u1x[7]);
    d2u0x[55] = scale*(-Ct[0]*Ct[5]*drill + D[1]*u1x[1]*u1x[6] + D[2]*u1x[0]*u1x[6] + D[4]*u1x[1]*u1x[7] + D[5]*u1x[0]*u1x[7]);
    d2u0x[56] = 0.0;
    d2u0x[57] = scale*(Ct[4]*Ct[5]*drill + D[0]*u1x[3]*u1x[6] + D[2]*u1x[3]*u1x[7] + D[2]*u1x[4]*u1x[6] + D[5]*u1x[4]*u1x[7]);
    d2u0x[58] = scale*(-Ct[1]*Ct[5]*drill + D[1]*u1x[4]*u1x[6] + D[2]*u1x[3]*u1x[6] + D[4]*u1x[4]*u1x[7] + D[5]*u1x[3]*u1x[7]);
    d2u0x[59] = 0.0;
    d2u0x[60] = scale*(Ct[5]*Ct[5]*drill + u1x[6]*(D[0]*u1x[6] + D[2]*u1x[7]) + u1x[7]*(D[2]*u1x[6] + D[5]*u1x[7]));
    d2u0x[61] = scale*(-Ct[2]*Ct[5]*drill + D[1]*u1x[6]*u1x[7] + D[2]*u1x[6]*u1x[6] + D[4]*u1x[7]*u1x[7] + D[5]*u1x[6]*u1x[7]);
    d2u0x[62] = 0.0;

    d2u0x[63] = scale*(-Ct[2]*Ct[3]*drill + D[1]*u1x[0]*u1x[7] + D[2]*u1x[0]*u1x[6] + D[4]*u1x[1]*u1x[7] + D[5]*u1x[1]*u1x[6]);
    d2u0x[64] = scale*(Ct[0]*Ct[2]*drill + D[3]*u1x[1]*u1x[7] + D[4]*u1x[0]*u1x[7] + D[4]*u1x[1]*u1x[6] + D[5]*u1x[0]*u1x[6]);
    d2u0x[65] = 0.0;
    d2u0x[66] = scale*(-Ct[2]*Ct[4]*drill + D[1]*u1x[3]*u1x[7] + D[2]*u1x[3]*u1x[6] + D[4]*u1x[4]*u1x[7] + D[5]*u1x[4]*u1x[6]);
    d2u0x[67] = scale*(Ct[1]*Ct[2]*drill + D[3]*u1x[4]*u1x[7] + D[4]*u1x[3]*u1x[7] + D[4]*u1x[4]*u1x[6] + D[5]*u1x[3]*u1x[6]);
    d2u0x[68] = 0.0;
    d2u0x[69] = scale*(-Ct[2]*Ct[5]*drill + D[1]*u1x[6]*u1x[7] + D[2]*u1x[6]*u1x[6] + D[4]*u1x[7]*u1x[7] + D[5]*u1x[6]*u1x[7]);
    d2u0x[70] = scale*(Ct[2]*Ct[2]*drill + u1x[6]*(D[4]*u1x[7] + D[5]*u1x[6]) + u1x[7]*(D[3]*u1x[7] + D[4]*u1x[6]));
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

    d2u0xu1x[0] = scale*(B[0]*e0ty[0] + B[1]*e0ty[3] + 2.0*B[2]*e0ty[1] + D[0]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[1]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[2]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[1]*(D[2]*u1x[0] + D[5]*u1x[1]) + 0.5*u1x[0]*(D[0]*(u0x[0] + 1.0) + D[2]*u0x[1]) + 0.5*u1x[1]*(D[2]*(u0x[0] + 1.0) + D[5]*u0x[1]) + 0.5*(u0x[0] + 1.0)*(D[0]*u1x[0] + D[2]*u1x[1]));
    d2u0xu1x[1] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[1]*(D[1]*u1x[0] + D[4]*u1x[1]) + 0.5*u1x[0]*(D[1]*u0x[1] + D[2]*(u0x[0] + 1.0)) + 0.5*u1x[1]*(D[4]*u0x[1] + D[5]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[2]*u1x[0] + D[5]*u1x[1]));
    d2u0xu1x[2] = 0.0;
    d2u0xu1x[3] = scale*(D[0]*u0x[3]*u1x[0] + D[2]*u0x[3]*u1x[1] + D[2]*u0x[4]*u1x[0] + D[2]*u1x[0] + D[5]*u0x[4]*u1x[1] + D[5]*u1x[1]);
    d2u0xu1x[4] = scale*(D[1]*u0x[4]*u1x[0] + D[1]*u1x[0] + D[2]*u0x[3]*u1x[0] + D[4]*u0x[4]*u1x[1] + D[4]*u1x[1] + D[5]*u0x[3]*u1x[1]);
    d2u0xu1x[5] = 0.0;
    d2u0xu1x[6] = scale*(D[0]*u0x[6]*u1x[0] + D[2]*u0x[6]*u1x[1] + D[2]*u0x[7]*u1x[0] + D[5]*u0x[7]*u1x[1]);
    d2u0xu1x[7] = scale*(D[1]*u0x[7]*u1x[0] + D[2]*u0x[6]*u1x[0] + D[4]*u0x[7]*u1x[1] + D[5]*u0x[6]*u1x[1]);
    d2u0xu1x[8] = 0.0;

    d2u0xu1x[9] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[1]*(D[4]*u1x[1] + D[5]*u1x[0]) + 0.5*u1x[0]*(D[2]*(u0x[0] + 1.0) + D[5]*u0x[1]) + 0.5*u1x[1]*(D[1]*(u0x[0] + 1.0) + D[4]*u0x[1]) + 0.5*(u0x[0] + 1.0)*(D[1]*u1x[1] + D[2]*u1x[0]));
    d2u0xu1x[10] = scale*(B[1]*e0ty[0] + B[3]*e0ty[3] + 2.0*B[4]*e0ty[1] + D[1]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[3]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[4]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[1]*(D[3]*u1x[1] + D[4]*u1x[0]) + 0.5*u1x[0]*(D[4]*u0x[1] + D[5]*(u0x[0] + 1.0)) + 0.5*u1x[1]*(D[3]*u0x[1] + D[4]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[4]*u1x[1] + D[5]*u1x[0]));
    d2u0xu1x[11] = 0.0;
    d2u0xu1x[12] = scale*(D[1]*u0x[3]*u1x[1] + D[2]*u0x[3]*u1x[0] + D[4]*u0x[4]*u1x[1] + D[4]*u1x[1] + D[5]*u0x[4]*u1x[0] + D[5]*u1x[0]);
    d2u0xu1x[13] = scale*(D[3]*u0x[4]*u1x[1] + D[3]*u1x[1] + D[4]*u0x[3]*u1x[1] + D[4]*u0x[4]*u1x[0] + D[4]*u1x[0] + D[5]*u0x[3]*u1x[0]);
    d2u0xu1x[14] = 0.0;
    d2u0xu1x[15] = scale*(D[1]*u0x[6]*u1x[1] + D[2]*u0x[6]*u1x[0] + D[4]*u0x[7]*u1x[1] + D[5]*u0x[7]*u1x[0]);
    d2u0xu1x[16] = scale*(D[3]*u0x[7]*u1x[1] + D[4]*u0x[6]*u1x[1] + D[4]*u0x[7]*u1x[0] + D[5]*u0x[6]*u1x[0]);
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

    d2u0xu1x[27] = scale*(D[0]*u0x[0]*u1x[3] + D[0]*u1x[3] + D[2]*u0x[0]*u1x[4] + D[2]*u0x[1]*u1x[3] + D[2]*u1x[4] + D[5]*u0x[1]*u1x[4]);
    d2u0xu1x[28] = scale*(D[1]*u0x[1]*u1x[3] + D[2]*u0x[0]*u1x[3] + D[2]*u1x[3] + D[4]*u0x[1]*u1x[4] + D[5]*u0x[0]*u1x[4] + D[5]*u1x[4]);
    d2u0xu1x[29] = 0.0;
    d2u0xu1x[30] = scale*(B[0]*e0ty[0] + B[1]*e0ty[3] + 2.0*B[2]*e0ty[1] + D[0]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[1]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[2]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[3]*(D[0]*u1x[3] + D[2]*u1x[4]) + 0.5*u1x[3]*(D[0]*u0x[3] + D[2]*(u0x[4] + 1.0)) + 0.5*u1x[4]*(D[2]*u0x[3] + D[5]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[2]*u1x[3] + D[5]*u1x[4]));
    d2u0xu1x[31] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[3]*(D[2]*u1x[3] + D[5]*u1x[4]) + 0.5*u1x[3]*(D[1]*(u0x[4] + 1.0) + D[2]*u0x[3]) + 0.5*u1x[4]*(D[4]*(u0x[4] + 1.0) + D[5]*u0x[3]) + 0.5*(u0x[4] + 1.0)*(D[1]*u1x[3] + D[4]*u1x[4]));
    d2u0xu1x[32] = 0.0;
    d2u0xu1x[33] = scale*(D[0]*u0x[6]*u1x[3] + D[2]*u0x[6]*u1x[4] + D[2]*u0x[7]*u1x[3] + D[5]*u0x[7]*u1x[4]);
    d2u0xu1x[34] = scale*(D[1]*u0x[7]*u1x[3] + D[2]*u0x[6]*u1x[3] + D[4]*u0x[7]*u1x[4] + D[5]*u0x[6]*u1x[4]);
    d2u0xu1x[35] = 0.0;

    d2u0xu1x[36] = scale*(D[1]*u0x[0]*u1x[4] + D[1]*u1x[4] + D[2]*u0x[0]*u1x[3] + D[2]*u1x[3] + D[4]*u0x[1]*u1x[4] + D[5]*u0x[1]*u1x[3]);
    d2u0xu1x[37] = scale*(D[3]*u0x[1]*u1x[4] + D[4]*u0x[0]*u1x[4] + D[4]*u0x[1]*u1x[3] + D[4]*u1x[4] + D[5]*u0x[0]*u1x[3] + D[5]*u1x[3]);
    d2u0xu1x[38] = 0.0;
    d2u0xu1x[39] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[3]*(D[1]*u1x[4] + D[2]*u1x[3]) + 0.5*u1x[3]*(D[2]*u0x[3] + D[5]*(u0x[4] + 1.0)) + 0.5*u1x[4]*(D[1]*u0x[3] + D[4]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[4]*u1x[4] + D[5]*u1x[3]));
    d2u0xu1x[40] = scale*(B[1]*e0ty[0] + B[3]*e0ty[3] + 2.0*B[4]*e0ty[1] + D[1]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[3]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[4]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[3]*(D[4]*u1x[4] + D[5]*u1x[3]) + 0.5*u1x[3]*(D[4]*(u0x[4] + 1.0) + D[5]*u0x[3]) + 0.5*u1x[4]*(D[3]*(u0x[4] + 1.0) + D[4]*u0x[3]) + 0.5*(u0x[4] + 1.0)*(D[3]*u1x[4] + D[4]*u1x[3]));
    d2u0xu1x[41] = 0.0;
    d2u0xu1x[42] = scale*(D[1]*u0x[6]*u1x[4] + D[2]*u0x[6]*u1x[3] + D[4]*u0x[7]*u1x[4] + D[5]*u0x[7]*u1x[3]);
    d2u0xu1x[43] = scale*(D[3]*u0x[7]*u1x[4] + D[4]*u0x[6]*u1x[4] + D[4]*u0x[7]*u1x[3] + D[5]*u0x[6]*u1x[3]);
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

    d2u0xu1x[54] = scale*(D[0]*u0x[0]*u1x[6] + D[0]*u1x[6] + D[2]*u0x[0]*u1x[7] + D[2]*u0x[1]*u1x[6] + D[2]*u1x[7] + D[5]*u0x[1]*u1x[7]);
    d2u0xu1x[55] = scale*(D[1]*u0x[1]*u1x[6] + D[2]*u0x[0]*u1x[6] + D[2]*u1x[6] + D[4]*u0x[1]*u1x[7] + D[5]*u0x[0]*u1x[7] + D[5]*u1x[7]);
    d2u0xu1x[56] = 0.0;
    d2u0xu1x[57] = scale*(D[0]*u0x[3]*u1x[6] + D[2]*u0x[3]*u1x[7] + D[2]*u0x[4]*u1x[6] + D[2]*u1x[6] + D[5]*u0x[4]*u1x[7] + D[5]*u1x[7]);
    d2u0xu1x[58] = scale*(D[1]*u0x[4]*u1x[6] + D[1]*u1x[6] + D[2]*u0x[3]*u1x[6] + D[4]*u0x[4]*u1x[7] + D[4]*u1x[7] + D[5]*u0x[3]*u1x[7]);
    d2u0xu1x[59] = 0.0;
    d2u0xu1x[60] = scale*(B[0]*e0ty[0] + B[1]*e0ty[3] + 2.0*B[2]*e0ty[1] + D[0]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[1]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[2]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[6]*(D[0]*u1x[6] + D[2]*u1x[7]) + 0.5*u0x[7]*(D[2]*u1x[6] + D[5]*u1x[7]) + 0.5*u1x[6]*(D[0]*u0x[6] + D[2]*u0x[7]) + 0.5*u1x[7]*(D[2]*u0x[6] + D[5]*u0x[7]));
    d2u0xu1x[61] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[6]*(D[2]*u1x[6] + D[5]*u1x[7]) + 0.5*u0x[7]*(D[1]*u1x[6] + D[4]*u1x[7]) + 0.5*u1x[6]*(D[1]*u0x[7] + D[2]*u0x[6]) + 0.5*u1x[7]*(D[4]*u0x[7] + D[5]*u0x[6]));
    d2u0xu1x[62] = 0.0;

    d2u0xu1x[63] = scale*(D[1]*u0x[0]*u1x[7] + D[1]*u1x[7] + D[2]*u0x[0]*u1x[6] + D[2]*u1x[6] + D[4]*u0x[1]*u1x[7] + D[5]*u0x[1]*u1x[6]);
    d2u0xu1x[64] = scale*(D[3]*u0x[1]*u1x[7] + D[4]*u0x[0]*u1x[7] + D[4]*u0x[1]*u1x[6] + D[4]*u1x[7] + D[5]*u0x[0]*u1x[6] + D[5]*u1x[6]);
    d2u0xu1x[65] = 0.0;
    d2u0xu1x[66] = scale*(D[1]*u0x[3]*u1x[7] + D[2]*u0x[3]*u1x[6] + D[4]*u0x[4]*u1x[7] + D[4]*u1x[7] + D[5]*u0x[4]*u1x[6] + D[5]*u1x[6]);
    d2u0xu1x[67] = scale*(D[3]*u0x[4]*u1x[7] + D[3]*u1x[7] + D[4]*u0x[3]*u1x[7] + D[4]*u0x[4]*u1x[6] + D[4]*u1x[6] + D[5]*u0x[3]*u1x[6]);
    d2u0xu1x[68] = 0.0;
    d2u0xu1x[69] = scale*(B[2]*e0ty[0] + B[4]*e0ty[3] + 2.0*B[5]*e0ty[1] + D[2]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[4]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[5]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[6]*(D[1]*u1x[7] + D[2]*u1x[6]) + 0.5*u0x[7]*(D[4]*u1x[7] + D[5]*u1x[6]) + 0.5*u1x[6]*(D[2]*u0x[6] + D[5]*u0x[7]) + 0.5*u1x[7]*(D[1]*u0x[6] + D[4]*u0x[7]));
    d2u0xu1x[70] = scale*(B[1]*e0ty[0] + B[3]*e0ty[3] + 2.0*B[4]*e0ty[1] + D[1]*(u0x[0]*u1x[0] + u0x[3]*u1x[3] + u0x[6]*u1x[6] + u1x[0]) + D[3]*(u0x[1]*u1x[1] + u0x[4]*u1x[4] + u0x[7]*u1x[7] + u1x[4]) + D[4]*(u0x[0]*u1x[1] + u0x[1]*u1x[0] + u0x[3]*u1x[4] + u0x[4]*u1x[3] + u0x[6]*u1x[7] + u0x[7]*u1x[6] + u1x[1] + u1x[3]) + 0.5*u0x[6]*(D[4]*u1x[7] + D[5]*u1x[6]) + 0.5*u0x[7]*(D[3]*u1x[7] + D[4]*u1x[6]) + 0.5*u1x[6]*(D[4]*u0x[7] + D[5]*u0x[6]) + 0.5*u1x[7]*(D[3]*u0x[7] + D[4]*u0x[6]));
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

    d2u1x[0] = scale*(u0x[1]*(D[2]*(u0x[0] + 1.0) + D[5]*u0x[1]) + (u0x[0] + 1.0)*(D[0]*(u0x[0] + 1.0) + D[2]*u0x[1]));
    d2u1x[1] = scale*(0.5*u0x[1]*(D[1]*(u0x[0] + 1.0) + D[4]*u0x[1]) + 0.5*u0x[1]*(D[4]*u0x[1] + D[5]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[1]*u0x[1] + D[2]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[2]*(u0x[0] + 1.0) + D[5]*u0x[1]));
    d2u1x[2] = 0.0;
    d2u1x[3] = scale*(D[0]*u0x[0]*u0x[3] + D[0]*u0x[3] + D[2]*u0x[0]*u0x[4] + D[2]*u0x[0] + D[2]*u0x[1]*u0x[3] + D[2]*u0x[4] + D[2] + D[5]*u0x[1]*u0x[4] + D[5]*u0x[1]);
    d2u1x[4] = scale*(D[1]*u0x[0]*u0x[4] + D[1]*u0x[0] + D[1]*u0x[4] + D[1] + D[2]*u0x[0]*u0x[3] + D[2]*u0x[3] + D[4]*u0x[1]*u0x[4] + D[4]*u0x[1] + D[5]*u0x[1]*u0x[3]);
    d2u1x[5] = 0.0;
    d2u1x[6] = scale*(D[0]*u0x[0]*u0x[6] + D[0]*u0x[6] + D[2]*u0x[0]*u0x[7] + D[2]*u0x[1]*u0x[6] + D[2]*u0x[7] + D[5]*u0x[1]*u0x[7]);
    d2u1x[7] = scale*(D[1]*u0x[0]*u0x[7] + D[1]*u0x[7] + D[2]*u0x[0]*u0x[6] + D[2]*u0x[6] + D[4]*u0x[1]*u0x[7] + D[5]*u0x[1]*u0x[6]);
    d2u1x[8] = 0.0;

    d2u1x[9] = scale*(0.5*u0x[1]*(D[1]*(u0x[0] + 1.0) + D[4]*u0x[1]) + 0.5*u0x[1]*(D[4]*u0x[1] + D[5]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[1]*u0x[1] + D[2]*(u0x[0] + 1.0)) + 0.5*(u0x[0] + 1.0)*(D[2]*(u0x[0] + 1.0) + D[5]*u0x[1]));
    d2u1x[10] = scale*(u0x[1]*(D[3]*u0x[1] + D[4]*(u0x[0] + 1.0)) + (u0x[0] + 1.0)*(D[4]*u0x[1] + D[5]*(u0x[0] + 1.0)));
    d2u1x[11] = 0.0;
    d2u1x[12] = scale*(D[1]*u0x[1]*u0x[3] + D[2]*u0x[0]*u0x[3] + D[2]*u0x[3] + D[4]*u0x[1]*u0x[4] + D[4]*u0x[1] + D[5]*u0x[0]*u0x[4] + D[5]*u0x[0] + D[5]*u0x[4] + D[5]);
    d2u1x[13] = scale*(D[3]*u0x[1]*u0x[4] + D[3]*u0x[1] + D[4]*u0x[0]*u0x[4] + D[4]*u0x[0] + D[4]*u0x[1]*u0x[3] + D[4]*u0x[4] + D[4] + D[5]*u0x[0]*u0x[3] + D[5]*u0x[3]);
    d2u1x[14] = 0.0;
    d2u1x[15] = scale*(D[1]*u0x[1]*u0x[6] + D[2]*u0x[0]*u0x[6] + D[2]*u0x[6] + D[4]*u0x[1]*u0x[7] + D[5]*u0x[0]*u0x[7] + D[5]*u0x[7]);
    d2u1x[16] = scale*(D[3]*u0x[1]*u0x[7] + D[4]*u0x[0]*u0x[7] + D[4]*u0x[1]*u0x[6] + D[4]*u0x[7] + D[5]*u0x[0]*u0x[6] + D[5]*u0x[6]);
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

    d2u1x[27] = scale*(D[0]*u0x[0]*u0x[3] + D[0]*u0x[3] + D[2]*u0x[0]*u0x[4] + D[2]*u0x[0] + D[2]*u0x[1]*u0x[3] + D[2]*u0x[4] + D[2] + D[5]*u0x[1]*u0x[4] + D[5]*u0x[1]);
    d2u1x[28] = scale*(D[1]*u0x[1]*u0x[3] + D[2]*u0x[0]*u0x[3] + D[2]*u0x[3] + D[4]*u0x[1]*u0x[4] + D[4]*u0x[1] + D[5]*u0x[0]*u0x[4] + D[5]*u0x[0] + D[5]*u0x[4] + D[5]);
    d2u1x[29] = 0.0;
    d2u1x[30] = scale*(u0x[3]*(D[0]*u0x[3] + D[2]*(u0x[4] + 1.0)) + (u0x[4] + 1.0)*(D[2]*u0x[3] + D[5]*(u0x[4] + 1.0)));
    d2u1x[31] = scale*(0.5*u0x[3]*(D[1]*(u0x[4] + 1.0) + D[2]*u0x[3]) + 0.5*u0x[3]*(D[2]*u0x[3] + D[5]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[1]*u0x[3] + D[4]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[4]*(u0x[4] + 1.0) + D[5]*u0x[3]));
    d2u1x[32] = 0.0;
    d2u1x[33] = scale*(D[0]*u0x[3]*u0x[6] + D[2]*u0x[3]*u0x[7] + D[2]*u0x[4]*u0x[6] + D[2]*u0x[6] + D[5]*u0x[4]*u0x[7] + D[5]*u0x[7]);
    d2u1x[34] = scale*(D[1]*u0x[3]*u0x[7] + D[2]*u0x[3]*u0x[6] + D[4]*u0x[4]*u0x[7] + D[4]*u0x[7] + D[5]*u0x[4]*u0x[6] + D[5]*u0x[6]);
    d2u1x[35] = 0.0;

    d2u1x[36] = scale*(D[1]*u0x[0]*u0x[4] + D[1]*u0x[0] + D[1]*u0x[4] + D[1] + D[2]*u0x[0]*u0x[3] + D[2]*u0x[3] + D[4]*u0x[1]*u0x[4] + D[4]*u0x[1] + D[5]*u0x[1]*u0x[3]);
    d2u1x[37] = scale*(D[3]*u0x[1]*u0x[4] + D[3]*u0x[1] + D[4]*u0x[0]*u0x[4] + D[4]*u0x[0] + D[4]*u0x[1]*u0x[3] + D[4]*u0x[4] + D[4] + D[5]*u0x[0]*u0x[3] + D[5]*u0x[3]);
    d2u1x[38] = 0.0;
    d2u1x[39] = scale*(0.5*u0x[3]*(D[1]*(u0x[4] + 1.0) + D[2]*u0x[3]) + 0.5*u0x[3]*(D[2]*u0x[3] + D[5]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[1]*u0x[3] + D[4]*(u0x[4] + 1.0)) + 0.5*(u0x[4] + 1.0)*(D[4]*(u0x[4] + 1.0) + D[5]*u0x[3]));
    d2u1x[40] = scale*(u0x[3]*(D[4]*(u0x[4] + 1.0) + D[5]*u0x[3]) + (u0x[4] + 1.0)*(D[3]*(u0x[4] + 1.0) + D[4]*u0x[3]));
    d2u1x[41] = 0.0;
    d2u1x[42] = scale*(D[1]*u0x[4]*u0x[6] + D[1]*u0x[6] + D[2]*u0x[3]*u0x[6] + D[4]*u0x[4]*u0x[7] + D[4]*u0x[7] + D[5]*u0x[3]*u0x[7]);
    d2u1x[43] = scale*(D[3]*u0x[4]*u0x[7] + D[3]*u0x[7] + D[4]*u0x[3]*u0x[7] + D[4]*u0x[4]*u0x[6] + D[4]*u0x[6] + D[5]*u0x[3]*u0x[6]);
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

    d2u1x[54] = scale*(D[0]*u0x[0]*u0x[6] + D[0]*u0x[6] + D[2]*u0x[0]*u0x[7] + D[2]*u0x[1]*u0x[6] + D[2]*u0x[7] + D[5]*u0x[1]*u0x[7]);
    d2u1x[55] = scale*(D[1]*u0x[1]*u0x[6] + D[2]*u0x[0]*u0x[6] + D[2]*u0x[6] + D[4]*u0x[1]*u0x[7] + D[5]*u0x[0]*u0x[7] + D[5]*u0x[7]);
    d2u1x[56] = 0.0;
    d2u1x[57] = scale*(D[0]*u0x[3]*u0x[6] + D[2]*u0x[3]*u0x[7] + D[2]*u0x[4]*u0x[6] + D[2]*u0x[6] + D[5]*u0x[4]*u0x[7] + D[5]*u0x[7]);
    d2u1x[58] = scale*(D[1]*u0x[4]*u0x[6] + D[1]*u0x[6] + D[2]*u0x[3]*u0x[6] + D[4]*u0x[4]*u0x[7] + D[4]*u0x[7] + D[5]*u0x[3]*u0x[7]);
    d2u1x[59] = 0.0;
    d2u1x[60] = scale*(D[0]*u0x[6]*u0x[6] + 2.0*D[2]*u0x[6]*u0x[7] + D[5]*u0x[7]*u0x[7]);
    d2u1x[61] = scale*(D[1]*u0x[6]*u0x[7] + D[2]*u0x[6]*u0x[6] + D[4]*u0x[7]*u0x[7] + D[5]*u0x[6]*u0x[7]);
    d2u1x[62] = 0.0;

    d2u1x[63] = scale*(D[1]*u0x[0]*u0x[7] + D[1]*u0x[7] + D[2]*u0x[0]*u0x[6] + D[2]*u0x[6] + D[4]*u0x[1]*u0x[7] + D[5]*u0x[1]*u0x[6]);
    d2u1x[64] = scale*(D[3]*u0x[1]*u0x[7] + D[4]*u0x[0]*u0x[7] + D[4]*u0x[1]*u0x[6] + D[4]*u0x[7] + D[5]*u0x[0]*u0x[6] + D[5]*u0x[6]);
    d2u1x[65] = 0.0;
    d2u1x[66] = scale*(D[1]*u0x[3]*u0x[7] + D[2]*u0x[3]*u0x[6] + D[4]*u0x[4]*u0x[7] + D[4]*u0x[7] + D[5]*u0x[4]*u0x[6] + D[5]*u0x[6]);
    d2u1x[67] = scale*(D[3]*u0x[4]*u0x[7] + D[3]*u0x[7] + D[4]*u0x[3]*u0x[7] + D[4]*u0x[4]*u0x[6] + D[4]*u0x[6] + D[5]*u0x[3]*u0x[6]);
    d2u1x[68] = 0.0;
    d2u1x[69] = scale*(D[1]*u0x[6]*u0x[7] + D[2]*u0x[6]*u0x[6] + D[4]*u0x[7]*u0x[7] + D[5]*u0x[6]*u0x[7]);
    d2u1x[70] = scale*(D[3]*u0x[7]*u0x[7] + 2.0*D[4]*u0x[6]*u0x[7] + D[5]*u0x[6]*u0x[6]);
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

    d2e0tyu0x[0] = scale*(B[0]*u1x[0] + B[2]*u1x[1]);
    d2e0tyu0x[1] = scale*(B[1]*u1x[1] + B[2]*u1x[0]);
    d2e0tyu0x[2] = 0.0;
    d2e0tyu0x[3] = scale*(B[0]*u1x[3] + B[2]*u1x[4]);
    d2e0tyu0x[4] = scale*(B[1]*u1x[4] + B[2]*u1x[3]);
    d2e0tyu0x[5] = 0.0;
    d2e0tyu0x[6] = scale*(B[0]*u1x[6] + B[2]*u1x[7]);
    d2e0tyu0x[7] = scale*(B[1]*u1x[7] + B[2]*u1x[6]);
    d2e0tyu0x[8] = 0.0;

    d2e0tyu0x[9] = 2.0*scale*(B[2]*u1x[0] + B[5]*u1x[1]);
    d2e0tyu0x[10] = 2.0*scale*(B[4]*u1x[1] + B[5]*u1x[0]);
    d2e0tyu0x[11] = 0.0;
    d2e0tyu0x[12] = 2.0*scale*(B[2]*u1x[3] + B[5]*u1x[4]);
    d2e0tyu0x[13] = 2.0*scale*(B[4]*u1x[4] + B[5]*u1x[3]);
    d2e0tyu0x[14] = 0.0;
    d2e0tyu0x[15] = 2.0*scale*(B[2]*u1x[6] + B[5]*u1x[7]);
    d2e0tyu0x[16] = 2.0*scale*(B[4]*u1x[7] + B[5]*u1x[6]);
    d2e0tyu0x[17] = 0.0;

    d2e0tyu0x[18] = 0.0;
    d2e0tyu0x[19] = 0.0;
    d2e0tyu0x[20] = 0.0;
    d2e0tyu0x[21] = 0.0;
    d2e0tyu0x[22] = 0.0;
    d2e0tyu0x[23] = 0.0;
    d2e0tyu0x[24] = 0.0;
    d2e0tyu0x[25] = 0.0;
    d2e0tyu0x[26] = 0.0;

    d2e0tyu0x[27] = scale*(B[1]*u1x[0] + B[4]*u1x[1]);
    d2e0tyu0x[28] = scale*(B[3]*u1x[1] + B[4]*u1x[0]);
    d2e0tyu0x[29] = 0.0;
    d2e0tyu0x[30] = scale*(B[1]*u1x[3] + B[4]*u1x[4]);
    d2e0tyu0x[31] = scale*(B[3]*u1x[4] + B[4]*u1x[3]);
    d2e0tyu0x[32] = 0.0;
    d2e0tyu0x[33] = scale*(B[1]*u1x[6] + B[4]*u1x[7]);
    d2e0tyu0x[34] = scale*(B[3]*u1x[7] + B[4]*u1x[6]);
    d2e0tyu0x[35] = 0.0;

    d2e0tyu0x[36] = 0.0;
    d2e0tyu0x[37] = 0.0;
    d2e0tyu0x[38] = 0.0;
    d2e0tyu0x[39] = 0.0;
    d2e0tyu0x[40] = 0.0;
    d2e0tyu0x[41] = 0.0;
    d2e0tyu0x[42] = 0.0;
    d2e0tyu0x[43] = 0.0;
    d2e0tyu0x[44] = 0.0;

    d2e0tyu0x[45] = 0.0;
    d2e0tyu0x[46] = 0.0;
    d2e0tyu0x[47] = 0.0;
    d2e0tyu0x[48] = 0.0;
    d2e0tyu0x[49] = 0.0;
    d2e0tyu0x[50] = 0.0;
    d2e0tyu0x[51] = 0.0;
    d2e0tyu0x[52] = 0.0;
    d2e0tyu0x[53] = 0.0;

    d2e0tyu1x[0] = scale*(B[0]*(u0x[0] + 1.0) + B[2]*u0x[1]);
    d2e0tyu1x[1] = scale*(B[1]*u0x[1] + B[2]*(u0x[0] + 1.0));
    d2e0tyu1x[2] = 0.0;
    d2e0tyu1x[3] = scale*(B[0]*u0x[3] + B[2]*(u0x[4] + 1.0));
    d2e0tyu1x[4] = scale*(B[1]*(u0x[4] + 1.0) + B[2]*u0x[3]);
    d2e0tyu1x[5] = 0.0;
    d2e0tyu1x[6] = scale*(B[0]*u0x[6] + B[2]*u0x[7]);
    d2e0tyu1x[7] = scale*(B[1]*u0x[7] + B[2]*u0x[6]);
    d2e0tyu1x[8] = 0.0;

    d2e0tyu1x[9] = 2.0*scale*(B[2]*(u0x[0] + 1.0) + B[5]*u0x[1]);
    d2e0tyu1x[10] = 2.0*scale*(B[4]*u0x[1] + B[5]*(u0x[0] + 1.0));
    d2e0tyu1x[11] = 0.0;
    d2e0tyu1x[12] = 2.0*scale*(B[2]*u0x[3] + B[5]*(u0x[4] + 1.0));
    d2e0tyu1x[13] = 2.0*scale*(B[4]*(u0x[4] + 1.0) + B[5]*u0x[3]);
    d2e0tyu1x[14] = 0.0;
    d2e0tyu1x[15] = 2.0*scale*(B[2]*u0x[6] + B[5]*u0x[7]);
    d2e0tyu1x[16] = 2.0*scale*(B[4]*u0x[7] + B[5]*u0x[6]);
    d2e0tyu1x[17] = 0.0;

    d2e0tyu1x[18] = 0.0;
    d2e0tyu1x[19] = 0.0;
    d2e0tyu1x[20] = 0.0;
    d2e0tyu1x[21] = 0.0;
    d2e0tyu1x[22] = 0.0;
    d2e0tyu1x[23] = 0.0;
    d2e0tyu1x[24] = 0.0;
    d2e0tyu1x[25] = 0.0;
    d2e0tyu1x[26] = 0.0;

    d2e0tyu1x[27] = scale*(B[1]*(u0x[0] + 1.0) + B[4]*u0x[1]);
    d2e0tyu1x[28] = scale*(B[3]*u0x[1] + B[4]*(u0x[0] + 1.0));
    d2e0tyu1x[29] = 0.0;
    d2e0tyu1x[30] = scale*(B[1]*u0x[3] + B[4]*(u0x[4] + 1.0));
    d2e0tyu1x[31] = scale*(B[3]*(u0x[4] + 1.0) + B[4]*u0x[3]);
    d2e0tyu1x[32] = 0.0;
    d2e0tyu1x[33] = scale*(B[1]*u0x[6] + B[4]*u0x[7]);
    d2e0tyu1x[34] = scale*(B[3]*u0x[7] + B[4]*u0x[6]);
    d2e0tyu1x[35] = 0.0;

    d2e0tyu1x[36] = 0.0;
    d2e0tyu1x[37] = 0.0;
    d2e0tyu1x[38] = 0.0;
    d2e0tyu1x[39] = 0.0;
    d2e0tyu1x[40] = 0.0;
    d2e0tyu1x[41] = 0.0;
    d2e0tyu1x[42] = 0.0;
    d2e0tyu1x[43] = 0.0;
    d2e0tyu1x[44] = 0.0;

    d2e0tyu1x[45] = 0.0;
    d2e0tyu1x[46] = 0.0;
    d2e0tyu1x[47] = 0.0;
    d2e0tyu1x[48] = 0.0;
    d2e0tyu1x[49] = 0.0;
    d2e0tyu1x[50] = 0.0;
    d2e0tyu1x[51] = 0.0;
    d2e0tyu1x[52] = 0.0;
    d2e0tyu1x[53] = 0.0;

    d2Ct[0] = scale*(drill*u0x[1]*u0x[1]);
    d2Ct[1] = scale*(drill*u0x[1]*(u0x[4] + 1.0));
    d2Ct[2] = scale*(drill*u0x[1]*u0x[7]);
    d2Ct[3] = scale*(-drill*u0x[1]*(u0x[0] + 1.0));
    d2Ct[4] = scale*(-drill*u0x[1]*u0x[3]);
    d2Ct[5] = scale*(-drill*u0x[1]*u0x[6]);
    d2Ct[6] = 0.0;
    d2Ct[7] = 0.0;
    d2Ct[8] = 0.0;

    d2Ct[9] = scale*(drill*u0x[1]*(u0x[4] + 1.0));
    d2Ct[10] = scale*(0.5*drill*(u0x[4] + 1.0)*(2.0*u0x[4] + 2.0));
    d2Ct[11] = scale*(drill*u0x[7]*(u0x[4] + 1.0));
    d2Ct[12] = scale*(-0.5*drill*(u0x[0] + 1.0)*(2.0*u0x[4] + 2.0));
    d2Ct[13] = scale*(-drill*u0x[3]*(u0x[4] + 1.0));
    d2Ct[14] = scale*(-drill*u0x[6]*(u0x[4] + 1.0));
    d2Ct[15] = 0.0;
    d2Ct[16] = 0.0;
    d2Ct[17] = 0.0;

    d2Ct[18] = scale*(drill*u0x[1]*u0x[7]);
    d2Ct[19] = scale*(drill*u0x[7]*(u0x[4] + 1.0));
    d2Ct[20] = scale*(drill*u0x[7]*u0x[7]);
    d2Ct[21] = scale*(-drill*u0x[7]*(u0x[0] + 1.0));
    d2Ct[22] = scale*(-drill*u0x[3]*u0x[7]);
    d2Ct[23] = scale*(-drill*u0x[6]*u0x[7]);
    d2Ct[24] = 0.0;
    d2Ct[25] = 0.0;
    d2Ct[26] = 0.0;

    d2Ct[27] = scale*(-drill*u0x[1]*(u0x[0] + 1.0));
    d2Ct[28] = scale*(-0.5*drill*(2.0*u0x[0] + 2.0)*(u0x[4] + 1.0));
    d2Ct[29] = scale*(-drill*u0x[7]*(u0x[0] + 1.0));
    d2Ct[30] = scale*(0.5*drill*(u0x[0] + 1.0)*(2.0*u0x[0] + 2.0));
    d2Ct[31] = scale*(drill*u0x[3]*(u0x[0] + 1.0));
    d2Ct[32] = scale*(drill*u0x[6]*(u0x[0] + 1.0));
    d2Ct[33] = 0.0;
    d2Ct[34] = 0.0;
    d2Ct[35] = 0.0;

    d2Ct[36] = scale*(-drill*u0x[1]*u0x[3]);
    d2Ct[37] = scale*(-drill*u0x[3]*(u0x[4] + 1.0));
    d2Ct[38] = scale*(-drill*u0x[3]*u0x[7]);
    d2Ct[39] = scale*(drill*u0x[3]*(u0x[0] + 1.0));
    d2Ct[40] = scale*(drill*u0x[3]*u0x[3]);
    d2Ct[41] = scale*(drill*u0x[3]*u0x[6]);
    d2Ct[42] = 0.0;
    d2Ct[43] = 0.0;
    d2Ct[44] = 0.0;

    d2Ct[45] = scale*(-drill*u0x[1]*u0x[6]);
    d2Ct[46] = scale*(-drill*u0x[6]*(u0x[4] + 1.0));
    d2Ct[47] = scale*(-drill*u0x[6]*u0x[7]);
    d2Ct[48] = scale*(drill*u0x[6]*(u0x[0] + 1.0));
    d2Ct[49] = scale*(drill*u0x[3]*u0x[6]);
    d2Ct[50] = scale*(drill*u0x[6]*u0x[6]);
    d2Ct[51] = 0.0;
    d2Ct[52] = 0.0;
    d2Ct[53] = 0.0;

    d2Ct[54] = 0.0;
    d2Ct[55] = 0.0;
    d2Ct[56] = 0.0;
    d2Ct[57] = 0.0;
    d2Ct[58] = 0.0;
    d2Ct[59] = 0.0;
    d2Ct[60] = 0.0;
    d2Ct[61] = 0.0;
    d2Ct[62] = 0.0;

    d2Ct[63] = 0.0;
    d2Ct[64] = 0.0;
    d2Ct[65] = 0.0;
    d2Ct[66] = 0.0;
    d2Ct[67] = 0.0;
    d2Ct[68] = 0.0;
    d2Ct[69] = 0.0;
    d2Ct[70] = 0.0;
    d2Ct[71] = 0.0;

    d2Ct[72] = 0.0;
    d2Ct[73] = 0.0;
    d2Ct[74] = 0.0;
    d2Ct[75] = 0.0;
    d2Ct[76] = 0.0;
    d2Ct[77] = 0.0;
    d2Ct[78] = 0.0;
    d2Ct[79] = 0.0;
    d2Ct[80] = 0.0;

    d2Ctu0x[0] = scale*(-Ct[3]*drill*u0x[1]);
    d2Ctu0x[1] = scale*(drill*(2.0*Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0) + Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) - Ct[4]*u0x[3] - Ct[5]*u0x[6]));
    d2Ctu0x[2] = 0.0;
    d2Ctu0x[3] = scale*(-Ct[4]*drill*u0x[1]);
    d2Ctu0x[4] = scale*(Ct[1]*drill*u0x[1]);
    d2Ctu0x[5] = 0.0;
    d2Ctu0x[6] = scale*(-Ct[5]*drill*u0x[1]);
    d2Ctu0x[7] = scale*(Ct[2]*drill*u0x[1]);
    d2Ctu0x[8] = 0.0;

    d2Ctu0x[9] = scale*(-Ct[3]*drill*(u0x[4] + 1.0));
    d2Ctu0x[10] = scale*(Ct[0]*drill*(u0x[4] + 1.0));
    d2Ctu0x[11] = 0.0;
    d2Ctu0x[12] = scale*(-Ct[4]*drill*(u0x[4] + 1.0));
    d2Ctu0x[13] = scale*(drill*(Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0) + 0.5*Ct[1]*(2.0*u0x[4] + 2.0) + Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) - Ct[4]*u0x[3] - Ct[5]*u0x[6]));
    d2Ctu0x[14] = 0.0;
    d2Ctu0x[15] = scale*(-Ct[5]*drill*(u0x[4] + 1.0));
    d2Ctu0x[16] = scale*(Ct[2]*drill*(u0x[4] + 1.0));
    d2Ctu0x[17] = 0.0;

    d2Ctu0x[18] = scale*(-Ct[3]*drill*u0x[7]);
    d2Ctu0x[19] = scale*(Ct[0]*drill*u0x[7]);
    d2Ctu0x[20] = 0.0;
    d2Ctu0x[21] = scale*(-Ct[4]*drill*u0x[7]);
    d2Ctu0x[22] = scale*(Ct[1]*drill*u0x[7]);
    d2Ctu0x[23] = 0.0;
    d2Ctu0x[24] = scale*(-Ct[5]*drill*u0x[7]);
    d2Ctu0x[25] = scale*(drill*(Ct[0]*u0x[1] + Ct[1]*(u0x[4] + 1.0) + 2.0*Ct[2]*u0x[7] - Ct[3]*(u0x[0] + 1.0) - Ct[4]*u0x[3] - Ct[5]*u0x[6]));
    d2Ctu0x[26] = 0.0;

    d2Ctu0x[27] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) - Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + 0.5*Ct[3]*(2.0*u0x[0] + 2.0) + Ct[4]*u0x[3] + Ct[5]*u0x[6]));
    d2Ctu0x[28] = scale*(-Ct[0]*drill*(u0x[0] + 1.0));
    d2Ctu0x[29] = 0.0;
    d2Ctu0x[30] = scale*(Ct[4]*drill*(u0x[0] + 1.0));
    d2Ctu0x[31] = scale*(-Ct[1]*drill*(u0x[0] + 1.0));
    d2Ctu0x[32] = 0.0;
    d2Ctu0x[33] = scale*(Ct[5]*drill*(u0x[0] + 1.0));
    d2Ctu0x[34] = scale*(-Ct[2]*drill*(u0x[0] + 1.0));
    d2Ctu0x[35] = 0.0;

    d2Ctu0x[36] = scale*(Ct[3]*drill*u0x[3]);
    d2Ctu0x[37] = scale*(-Ct[0]*drill*u0x[3]);
    d2Ctu0x[38] = 0.0;
    d2Ctu0x[39] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) - Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + 2.0*Ct[4]*u0x[3] + Ct[5]*u0x[6]));
    d2Ctu0x[40] = scale*(-Ct[1]*drill*u0x[3]);
    d2Ctu0x[41] = 0.0;
    d2Ctu0x[42] = scale*(Ct[5]*drill*u0x[3]);
    d2Ctu0x[43] = scale*(-Ct[2]*drill*u0x[3]);
    d2Ctu0x[44] = 0.0;

    d2Ctu0x[45] = scale*(Ct[3]*drill*u0x[6]);
    d2Ctu0x[46] = scale*(-Ct[0]*drill*u0x[6]);
    d2Ctu0x[47] = 0.0;
    d2Ctu0x[48] = scale*(Ct[4]*drill*u0x[6]);
    d2Ctu0x[49] = scale*(-Ct[1]*drill*u0x[6]);
    d2Ctu0x[50] = 0.0;
    d2Ctu0x[51] = scale*(drill*(-Ct[0]*u0x[1] - Ct[1]*(u0x[4] + 1.0) - Ct[2]*u0x[7] + Ct[3]*(u0x[0] + 1.0) + Ct[4]*u0x[3] + 2.0*Ct[5]*u0x[6]));
    d2Ctu0x[52] = scale*(-Ct[2]*drill*u0x[6]);
    d2Ctu0x[53] = 0.0;

    d2Ctu0x[54] = 0.0;
    d2Ctu0x[55] = 0.0;
    d2Ctu0x[56] = 0.0;
    d2Ctu0x[57] = 0.0;
    d2Ctu0x[58] = 0.0;
    d2Ctu0x[59] = 0.0;
    d2Ctu0x[60] = 0.0;
    d2Ctu0x[61] = 0.0;

    d2Ctu0x[62] = 0.0;
    d2Ctu0x[63] = 0.0;
    d2Ctu0x[64] = 0.0;
    d2Ctu0x[65] = 0.0;
    d2Ctu0x[66] = 0.0;
    d2Ctu0x[67] = 0.0;
    d2Ctu0x[68] = 0.0;
    d2Ctu0x[69] = 0.0;
    d2Ctu0x[70] = 0.0;
    d2Ctu0x[71] = 0.0;

    d2Ctu0x[72] = 0.0;
    d2Ctu0x[73] = 0.0;
    d2Ctu0x[74] = 0.0;
    d2Ctu0x[75] = 0.0;
    d2Ctu0x[76] = 0.0;
    d2Ctu0x[77] = 0.0;
    d2Ctu0x[78] = 0.0;
    d2Ctu0x[79] = 0.0;
    d2Ctu0x[80] = 0.0;
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
    TacsPrintErrorComponents(stderr, "dCt", dCt, fdCt, 9);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
  TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
  TacsScalar d2Ct[81], d2Ctu0x[81];
  model::evalStrainHessian(detXd, s, Cs, u0x, u1x, e0ty, Ct,
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
