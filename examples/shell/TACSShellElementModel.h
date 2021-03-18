#ifndef TACS_SHELL_ELEMENT_MODEL_H
#define TACS_SHELL_ELEMENT_MODEL_H

#include "TACSElementAlgebra.h"

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
          TacsScalar n0[3], d0[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);

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
          TacsScalar n0[3], d0[3], dd0[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);

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
          TacsScalar n0[3], d0[3], d0d[3];
          basis::interpFields(pt, 3, fn, 3, n0);
          basis::interpFields(pt, 3, d, 3, d0);
          basis::interpFields(pt, 3, dd, 3, d0d);

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

  /**
    Evaluate the tensorial components of the strain tensor at the specific
    quadrature point

    gty = [g11  g12  g13]
          [sym  g22  g23]
          [sym  sym  g33]

    As a result: gty[0] = g11, gty[1] = g12, gty[2] = g13, gty[3] = g22
    and gty[4] = g23, with gty[5] = 0.0

    @param pt The quadrature point
    @param ety The strain computed at the tying points
    @param gty The interpolated tying strain
  */
  template <class basis>
  static void interpTyingStrain( const double pt[],
                                const TacsScalar ety[],
                                TacsScalar gty[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      gty[index[field]] = basis::interpTying(field, pt, ety);
      ety += basis::getNumTyingPoints(field);
    }
    gty[5] = 0.0; // g33 = 0.0
  }

  /*
    Add the derivative of the tying strain to the given component field

  */
  template <class basis>
  static void addInterpTyingStrainTranspose( const double pt[],
                                             const TacsScalar dgty[],
                                             TacsScalar dety[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      basis::addInterpTyingTranspose(field, pt, dgty[index[field]], dety);
      dety += basis::getNumTyingPoints(field);
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

  /**
    Evaluate the tensorial components of the strain tensor at the specific
    quadrature point

    gty = [g11  g12  g13]
          [sym  g22  g23]
          [sym  sym  g33]

    As a result: gty[0] = g11, gty[1] = g12, gty[2] = g13, gty[3] = g22
    and gty[4] = g23, with gty[5] = 0.0

    @param pt The quadrature point
    @param ety The strain computed at the tying points
    @param gty The interpolated tying strain
  */
  template <class basis>
  static void interpTyingStrain( const double pt[],
                                const TacsScalar ety[],
                                TacsScalar gty[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      gty[index[field]] = basis::interpTying(field, pt, ety);
      ety += basis::getNumTyingPoints(field);
    }
    gty[5] = 0.0;
  }

  /*
    Add the derivative of the tying strain to the given component field

  */
  template <class basis>
  static void addInterpTyingStrainTranspose( const double pt[],
                                             const TacsScalar dgty[],
                                             TacsScalar dety[] ){
    // Set the values into the strain tensor
    const int index[] = {0, 3, 1, 4, 2};
    const int num_tying_fields = 5;
    for ( int field = 0; field < num_tying_fields; field++ ){
      basis::addInterpTyingTranspose(field, pt, dgty[index[field]], dety);
      dety += basis::getNumTyingPoints(field);
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

#endif // TACS_SHELL_ELEMENT_MODEL_H
