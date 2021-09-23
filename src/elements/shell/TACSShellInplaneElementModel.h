#ifndef TACS_SHELL_INPLANE_ELEMENT_MODEL_H
#define TACS_SHELL_INPLANE_ELEMENT_MODEL_H

#include "TACSShellElementQuadBasis.h"
#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"

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
  static void computeTyingStrain( const TacsScalar Xpts[],
                                  const TacsScalar fn[],
                                  const TacsScalar vars[],
                                  const TacsScalar d[],
                                  TacsScalar ety[] ){
    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
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

      if (field == TACS_SHELL_G23_COMPONENT){
        // Compute g23 = e2^{T}*G*e3
        ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                          n0[0]*Uxi[1] + n0[1]*Uxi[3] + n0[2]*Uxi[5]);
      }
      else if (field == TACS_SHELL_G13_COMPONENT){
        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                          n0[0]*Uxi[0] + n0[1]*Uxi[2] + n0[2]*Uxi[4]);
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
    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
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

      if (field == TACS_SHELL_G23_COMPONENT){
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
      else if (field == TACS_SHELL_G13_COMPONENT){
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
      basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, dUxi, res);
    }
  }

  template <int vars_per_node, class basis>
  static void addComputeTyingStrainHessian( const TacsScalar alpha,
                                            const TacsScalar Xpts[],
                                            const TacsScalar fn[],
                                            const TacsScalar vars[],
                                            const TacsScalar d[],
                                            const TacsScalar dety[],
                                            const TacsScalar d2ety[],
                                            const TacsScalar d2etyu[],
                                            const TacsScalar d2etyd[],
                                            TacsScalar mat[],
                                            TacsScalar d2d[],
                                            TacsScalar d2du[] ){
    // Initialize the data
    TacsScalar n0ty[3*basis::NUM_TYING_POINTS];
    TacsScalar Xxity[6*basis::NUM_TYING_POINTS];
    TacsScalar *n0 = n0ty, *Xxi = Xxity;

    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
      // Get the tying point parametric location
      double pt[2];
      basis::getTyingPoint(index, pt);

      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      n0 += 3;
      Xxi += 6;
    }

    TacsScalar *n01 = n0ty, *Xxi1 = Xxity;
    for ( int i1 = 0; i1 < basis::NUM_TYING_POINTS; i1++, n01 += 3, Xxi1 += 6 ){
      // Get the field index
      const TacsShellTyingStrainComponent f1 = basis::getTyingField(i1);

      // Get the tying point parametric location
      double pt1[2];
      basis::getTyingPoint(i1, pt1);

      TacsScalar du2[3*basis::NUM_NODES], dd2[3*basis::NUM_NODES];
      memset(du2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
      memset(dd2, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

      TacsScalar *n02 = n0ty, *Xxi2 = Xxity;
      for ( int i2 = 0; i2 < basis::NUM_TYING_POINTS; i2++, n02 += 3, Xxi2 += 6 ){
        // Get the field index
        const TacsShellTyingStrainComponent f2 = basis::getTyingField(i2);

        // Get the tying point parametric location
        double pt2[2];
        basis::getTyingPoint(i2, pt2);

        TacsScalar value = d2ety[basis::NUM_TYING_POINTS*i1 + i2];

        TacsScalar dUxi2[6], dd02[3];
        if (f2 == TACS_SHELL_G23_COMPONENT){
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
        else if (f2 == TACS_SHELL_G13_COMPONENT){
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
        basis::template addInterpFieldsGradTranspose<3, 3>(pt2, dUxi2, du2);
      }

      TacsScalar du1[3*basis::NUM_NODES], dd1[3*basis::NUM_NODES];
      memset(du1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
      memset(dd1, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

      // Store the the derivative information for the first point
      TacsScalar dUxi1[6], dd01[3];
      if (f1 == TACS_SHELL_G23_COMPONENT){
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
      else if (f1 == TACS_SHELL_G13_COMPONENT){
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
      basis::template addInterpFieldsGradTranspose<3, 3>(pt1, dUxi1, du1);

      const TacsScalar *etd = &d2etyd[3*basis::NUM_NODES*i1];
      const TacsScalar *etu = &d2etyu[3*basis::NUM_NODES*i1];
      for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
        for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
          d2d[3*basis::NUM_NODES*i + j] += dd1[i]*dd2[j] + dd1[i]*etd[j] + etd[i]*dd1[j];
        }
      }

      for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
        for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
          d2du[3*basis::NUM_NODES*i + j] += dd1[i]*du2[j] + dd1[i]*etu[j];
        }
      }

      for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
        for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
          d2du[3*basis::NUM_NODES*i + j] += etd[i]*du1[j];
        }
      }

      const int nvars = vars_per_node*basis::NUM_NODES;
      for ( int i = 0; i < 3*basis::NUM_NODES; i++ ){
        int ii = vars_per_node*(i / 3) + (i % 3);
        for ( int j = 0; j < 3*basis::NUM_NODES; j++ ){
          int jj = vars_per_node*(j / 3) + (j % 3);
          mat[nvars*ii + jj] += du1[i]*du2[j] + du1[i]*etu[j] + etu[i]*du1[j];
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
    for ( int index = 0; index < basis::NUM_TYING_POINTS; index++ ){
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
      if (field == TACS_SHELL_G23_COMPONENT){
        // Compute g23 = e2^{T}*G*e3
        ety[index] = 0.5*(Xxi[1]*d0[0] + Xxi[3]*d0[1] + Xxi[5]*d0[2] +
                          n0[0]*Uxi[1] + n0[1]*Uxi[3] + n0[2]*Uxi[5]);
        etyd[index] = 0.5*(Xxi[1]*d0d[0] + Xxi[3]*d0d[1] + Xxi[5]*d0d[2] +
                            n0[0]*Uxid[1] + n0[1]*Uxid[3] + n0[2]*Uxid[5]);
      }
      else if (field == TACS_SHELL_G13_COMPONENT){
        // Compute g13 = e1^{T}*G*e3
        ety[index] = 0.5*(Xxi[0]*d0[0] + Xxi[2]*d0[1] + Xxi[4]*d0[2] +
                          n0[0]*Uxi[0] + n0[1]*Uxi[2] + n0[2]*Uxi[4]);
        etyd[index] = 0.5*(Xxi[0]*d0d[0] + Xxi[2]*d0d[1] + Xxi[4]*d0d[2] +
                            n0[0]*Uxid[0] + n0[1]*Uxid[2] + n0[2]*Uxid[4]);
      }
    }
  }

  /*
    Evaluate the strain as a function of the displacement derivatives
    and interpolated strain from the tensorial components
  */
  static inline void evalStrain( const TacsScalar u0x[],
                                 const TacsScalar u1x[],
                                 const TacsScalar e0ty[],
                                 TacsScalar e[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0];
    e[1] = u0x[4];
    e[2] = u0x[1] + u0x[3];

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
  static inline void evalStrainSens( const TacsScalar scale,
                                     const TacsScalar dfde[],
                                     const TacsScalar u0x[],
                                     const TacsScalar u1x[],
                                     TacsScalar du0x[],
                                     TacsScalar du1x[],
                                     TacsScalar de0ty[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    de0ty[0] = 0.0;
    de0ty[1] = 0.0;
    de0ty[2] = 2.0*scale*dfde[7];
    de0ty[3] = 0.0;
    de0ty[4] = 2.0*scale*dfde[6];
    de0ty[5] = 0.0;

    // Compute the derivative with respect to u0x
    du0x[0] = scale*dfde[0];
    du0x[1] = scale*dfde[2];
    du0x[2] = 0.0;
    du0x[3] = scale*dfde[2];
    du0x[4] = scale*dfde[1];
    du0x[5] = 0.0;
    du0x[6] = 0.0;
    du0x[7] = 0.0;
    du0x[8] = 0.0;

    // Compute the derivative with respect to u1x
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

  static inline void evalStrainDeriv( const TacsScalar u0x[],
                                      const TacsScalar u1x[],
                                      const TacsScalar e0ty[],
                                      const TacsScalar u0xd[],
                                      const TacsScalar u1xd[],
                                      const TacsScalar e0tyd[],
                                      TacsScalar e[],
                                      TacsScalar ed[] ){
    // Evaluate the in-plane strains from the tying strain expressions
    e[0] = u0x[0];
    e[1] = u0x[4];
    e[2] = u0x[1] + u0x[3];

    // Compute the bending strain
    e[3] = u1x[0];
    e[4] = u1x[4];
    e[5] = u1x[1] + u1x[3];

    // Add the components of the shear strain
    e[6] = 2.0*e0ty[4];
    e[7] = 2.0*e0ty[2];

    // Evaluate the in-plane strains from the tying strain expressions
    ed[0] = u0xd[0];
    ed[1] = u0xd[4];
    ed[2] = u0xd[1] + u0xd[3];

    // Compute the bending strain
    ed[3] = u1xd[0];
    ed[4] = u1xd[4];
    ed[5] = u1xd[1] + u1xd[3];

    // Add the components of the shear strain
    ed[6] = 2.0*e0tyd[4];
    ed[7] = 2.0*e0tyd[2];
  }

  static inline void evalStrainHessian( const TacsScalar scale,
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
    TacsScalar *d2;
    d2 = d2u0x;
    d2[0] = scale*A[0];
    d2[4] = scale*A[1];
    d2[1] = scale*A[2];
    d2[3] = scale*A[2];

    d2 = &d2u0x[4*9];
    d2[0] = scale*A[1];
    d2[4] = scale*A[3];
    d2[1] = scale*A[4];
    d2[3] = scale*A[4];

    d2 = &d2u0x[9];
    d2[0] = scale*A[2];
    d2[4] = scale*A[4];
    d2[1] = scale*A[5];
    d2[3] = scale*A[5];

    d2 = &d2u0x[3*9];
    d2[0] = scale*A[2];
    d2[4] = scale*A[4];
    d2[1] = scale*A[5];
    d2[3] = scale*A[5];

    // e[3] = u1x[0];
    // e[4] = u1x[4];
    // e[5] = u1x[1] + u1x[3];
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

    d2 = d2u0xu1x;
    d2[0] = scale*B[0];
    d2[4] = scale*B[1];
    d2[1] = scale*B[2];
    d2[3] = scale*B[2];

    d2 = &d2u0xu1x[4*9];
    d2[0] = scale*B[1];
    d2[4] = scale*B[3];
    d2[1] = scale*B[4];
    d2[3] = scale*B[4];

    d2 = &d2u0xu1x[9];
    d2[0] = scale*B[2];
    d2[4] = scale*B[4];
    d2[1] = scale*B[5];
    d2[3] = scale*B[5];

    d2 = &d2u0xu1x[3*9];
    d2[0] = scale*B[2];
    d2[4] = scale*B[4];
    d2[1] = scale*B[5];
    d2[3] = scale*B[5];

    // e[6] = 2.0*e0ty[4];
    // e[7] = 2.0*e0ty[2];
    d2 = &d2e0ty[4*6];
    d2[4] = 4.0*scale*As[0];
    d2[2] = 4.0*scale*As[1];

    d2 = &d2e0ty[2*6];
    d2[4] = 4.0*scale*As[1];
    d2[2] = 4.0*scale*As[2];
  }
};

#endif // TACS_SHELL_INPLANE_ELEMENT_MODEL_H
