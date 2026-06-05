/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_MATERIAL_PROPERTIES_H
#define TACS_MATERIAL_PROPERTIES_H

/*
  The following classes contain the data necessary to compute the
  constitutive relationships or failure data.
*/

#include "TACSObject.h"
#include "TacsUtilities.h"

enum MaterialType { TACS_ISOTROPIC_MATERIAL, TACS_ANISOTROPIC_MATERIAL };

/**
   This class stores the mechanical and thermal material properties for
   isotropic and anisotropic materials.

   The goal of this class is to store a set of material properties
   that can be queried by constitutive classes for beams, shells,
   plane stress and solid elements. The minimum set of properties
   consists of isotropic mechanical properties, with zero thermal
   expansion.
*/
class TACSMaterialProperties : public TACSObject {
 public:
  TACSMaterialProperties(TacsScalar _rho, TacsScalar _specific_heat,
                         TacsScalar _E, TacsScalar _nu, TacsScalar _ys,
                         TacsScalar _alpha, TacsScalar _kappa);
  TACSMaterialProperties(
      TacsScalar _rho, TacsScalar _specific_heat, TacsScalar _E1,
      TacsScalar _E2, TacsScalar _E3, TacsScalar _nu12, TacsScalar _nu13,
      TacsScalar _nu23, TacsScalar _G12, TacsScalar _G13, TacsScalar _G23,
      TacsScalar _T1 = 0.0, TacsScalar _C1 = 0.0, TacsScalar _T2 = 0.0,
      TacsScalar _C2 = 0.0, TacsScalar _T3 = 0.0, TacsScalar _C3 = 0.0,
      TacsScalar _S12 = 0.0, TacsScalar _S13 = 0.0, TacsScalar _S23 = 0.0,
      TacsScalar _alpha1 = 0.0, TacsScalar _alpha2 = 0.0,
      TacsScalar _alpha3 = 0.0, TacsScalar _kappa1 = 0.0,
      TacsScalar _kappa2 = 0.0, TacsScalar _kappa3 = 0.0,
      TacsScalar _b_tt = 0.0, TacsScalar _b_tl = 0.0, TacsScalar _muWF = 0.0,
      TacsScalar _mu3W = 0.0, TacsScalar _mu3F = 0.0, TacsScalar _m = 0.0);
  ~TACSMaterialProperties() {}

  // Get the material type
  MaterialType getMaterialType();

  // Extract material property values
  TacsScalar getDensity();
  TacsScalar getSpecificHeat();

  // Set material property values
  void setDensity(TacsScalar _rho);
  void setSpecificHeat(TacsScalar _specific_heat);

  // Extract the coefficients
  void getIsotropicProperties(TacsScalar *_E, TacsScalar *_nu);
  void getOrthotropicProperties(TacsScalar *_E1, TacsScalar *_E2,
                                TacsScalar *_E3, TacsScalar *_nu12,
                                TacsScalar *_nu13, TacsScalar *_nu23,
                                TacsScalar *_G12, TacsScalar *_G13,
                                TacsScalar *_G23);
  void getStrengthProperties(TacsScalar *_T1, TacsScalar *_C1, TacsScalar *_T2,
                             TacsScalar *_C2, TacsScalar *_T3, TacsScalar *_C3,
                             TacsScalar *_S12, TacsScalar *_S13,
                             TacsScalar *_S23, TacsScalar *_b_tt,
                             TacsScalar *_b_tl, TacsScalar *_muWF,
                             TacsScalar *_mu3W, TacsScalar *_mu3F,
                             TacsScalar *_m);
  void getCoefThermalExpansion(TacsScalar *_a1, TacsScalar *_a2,
                               TacsScalar *_a3);
  void getThermalConductivity(TacsScalar *_k1, TacsScalar *_k2,
                              TacsScalar *_k3);

  // Evaluate the constitutive relationships
  void evalTangentStiffness3D(TacsScalar C[]);
  void evalTangentStiffness2D(TacsScalar C[]);

  // Evaluate the thermal conductivity matrices
  void evalTangentHeatFlux3D(TacsScalar Kc[]);
  void evalTangentHeatFlux2D(TacsScalar Kc[]);

  // Evaluate the thermal strain
  void evalThermalStrain3D(TacsScalar e[]);
  void evalThermalStrain2D(TacsScalar e[]);

  // Given the strain state, compute the stresses
  void evalStress2D(const TacsScalar e[], TacsScalar s[]);
  void evalStress3D(const TacsScalar e[], TacsScalar s[]);

  // Compute the failure index using the von Mises stress
  TacsScalar vonMisesFailure3D(const TacsScalar stress[]);
  TacsScalar vonMisesFailure3DStressSens(const TacsScalar stress[],
                                         TacsScalar sens[]);

  // The von Mises failure criteria for plane stress problems
  TacsScalar vonMisesFailure2D(const TacsScalar stress[]);
  TacsScalar vonMisesFailure2DStressSens(const TacsScalar stress[],
                                         TacsScalar sens[]);

 private:
  MaterialType mat_type;

  // Density
  TacsScalar rho;

  // Specific heat
  TacsScalar specific_heat;

  // Isotropic properties
  TacsScalar E, nu, G;  // Modulus and Poisson ratio

  // Coefficient of thermal expansion
  TacsScalar alpha;

  // Yield stress value
  TacsScalar ys;

  // Thermal conductivity coefficient
  TacsScalar kappa;

  // The orthotropic material properties
  TacsScalar E1, E2, E3, nu12, nu13, nu23, G12, G13, G23;

  // The strength coefficients
  TacsScalar T1, C1, T2, C2, T3, C3, S12, S13, S23, b_tt, b_tl, muWF, mu3W,
      mu3F, m;

  // The thermal coefficients of expansion
  TacsScalar alpha1, alpha2, alpha3;

  // Thermal conductivity coefficients along each axis
  TacsScalar kappa1, kappa2, kappa3;
};

/*
  The following class holds the material stiffness and strength
  properties for an orthotropic ply. This class is used by several
  constitutive classes within TACS.

  The stiffness of the ply in the global axis is supplied by the
  functions calculateAbar/calculateQbar, where the argument 'angle' is
  provided in radians. The failure properties are calculated either
  based on a Tsai-Wu failure criterion, the Cuntze failure criterion,
  or a smoothed maximum strain failure criterion, where the smoothing
  is performed using a KS function.

  The interaction coefficient for the Tsai-Wu failure criterion is set
  to zero by default. If a value of C, the failure stress under
  combined in-plane loading, is supplied, the interaction coefficient
  is determined. Be careful - the value can easily fall outside
  acceptable bounds - these are tested during initialization.
*/
class TACSOrthotropicPly : public TACSObject {
 public:
  // Create OrthoPly with a full set of orthotropic material relationships
  // ---------------------------------------------------------------------
  TACSOrthotropicPly(TacsScalar _plyThickness,
                     TACSMaterialProperties *_properties);

  enum CompositeFailureCriterion {
    MAX_STRAIN,        // Smoothed max-strain criterion (KS)
    TSAI_WU,           // Standard Tsai-Wu failure index
    TSAI_WU_MODIFIED,  // Modified Tsai-Wu returning a strength ratio (default)
    CUNTZE_UD,         // Cuntze criterion for unidirectional plies
    CUNTZE_WOVEN       // Cuntze criterion for woven plies
  };

  void setKSWeight(TacsScalar _ksWeight);
  void setFailureCriterion(CompositeFailureCriterion fc) {
    failureCriterion = fc;
  }

  // Retrieve the material properties
  // --------------------------------
  TacsScalar getDensity();
  TacsScalar getPlyThickness();
  void getStiffness(TacsScalar *_E1, TacsScalar *_E2, TacsScalar *_nu12,
                    TacsScalar *_G12, TacsScalar *_G23, TacsScalar *_G13);
  void getLaminateStiffness(TacsScalar *_Q11, TacsScalar *_Q12,
                            TacsScalar *_Q22, TacsScalar *_Q44,
                            TacsScalar *_Q55, TacsScalar *_Q66);
  void getStrength(TacsScalar *_Xt, TacsScalar *_Xc, TacsScalar *_Yt,
                   TacsScalar *_Yc, TacsScalar *_S12, TacsScalar *_C);
  void getStrainStrength(TacsScalar *_eXt, TacsScalar *_eXc, TacsScalar *_eYt,
                         TacsScalar *_eYc, TacsScalar *_eS12);
  void getTsaiWu(TacsScalar *_F1, TacsScalar *_F2, TacsScalar *_F11,
                 TacsScalar *_F12, TacsScalar *_F22, TacsScalar *_F66);
  void getCuntzeConstants(TacsScalar *_b_tl, TacsScalar *_muWF, TacsScalar *_m);
  void getLaminateInvariants(TacsScalar *U1, TacsScalar *U2, TacsScalar *U3,
                             TacsScalar *U4, TacsScalar *U5, TacsScalar *U6);

  // Calculate the Abar and Qbar matrices and their derivatives
  // Taken from Jones, Mechanics of composite materials pg. 51
  // ----------------------------------------------------------
  void calculateAbar(TacsScalar angle, TacsScalar Abar[]);
  void calculateAbarAngleSens(TacsScalar angle, TacsScalar Abar[]);
  void calculateQbar(TacsScalar angle, TacsScalar Qbar[]);
  void calculateQbarAngleSens(TacsScalar angle, TacsScalar Qbar[]);

  // Given the strain in the laminate coordinates, determine the
  // stress in the laminate coordinates
  // -----------------------------------------------------------
  void calculateStress(TacsScalar angle, const TacsScalar strain[],
                       TacsScalar stress[]);

  // Given the strain in the laminate coordinates, determine the failure load
  // ------------------------------------------------------------------------
  TacsScalar failure(TacsScalar angle, const TacsScalar strain[]);
  TacsScalar failureStrainSens(TacsScalar angle, const TacsScalar strain[],
                               TacsScalar sens[]);
  TacsScalar failureAngleSens(TacsScalar angle, const TacsScalar strain[],
                              TacsScalar *failSens);

  // Given the strain and stress in the local coordinates,
  // determine the failure modes of the Cuntze failure criterion
  // as well as the global failure value.
  // -----------------------------------------------------------
  TacsScalar CuntzeUD_FailureModes(const TacsScalar e[], const TacsScalar s[],
                                   TacsScalar *eff_ff1, TacsScalar *eff_ff2,
                                   TacsScalar *eff_iff1, TacsScalar *eff_iff2,
                                   TacsScalar *eff_iff3);

  TacsScalar CuntzeWoven_FailureModes(
      const TacsScalar e[], const TacsScalar s[], TacsScalar *eff_ff1,
      TacsScalar *eff_ff2, TacsScalar *eff_ff3, TacsScalar *eff_ff4,
      TacsScalar *eff_iff1, TacsScalar *eff_iff2, TacsScalar *eff_iff3,
      TacsScalar *eff_iff4, TacsScalar *eff_iff5);

  // Calculate the failure load fraction for given
  // constant and linear strain components
  // ---------------------------------------------
  TacsScalar calculateFailLoad(TacsScalar angle, const TacsScalar cstrain[],
                               const TacsScalar lstrain[]);
  TacsScalar calculateFailLoadStrainSens(TacsScalar angle,
                                         const TacsScalar cstrain[],
                                         const TacsScalar lstrain[],
                                         TacsScalar cSens[],
                                         TacsScalar lSens[]);
  TacsScalar calculateFailLoadAngleSens(TacsScalar angle,
                                        const TacsScalar cstrain[],
                                        const TacsScalar lstrain[],
                                        TacsScalar *posSens);
  void getPlyStress(const TacsScalar strain[], TacsScalar stress[]);

  // Transform the stress and strain between global/local frames
  // -----------------------------------------------------------
  /*
    For a ply at an angle t, the stress in the local ply frame is computed from
    the stress in the global frame as:

    s_local = T(t) * s_global

    where T is the transformation matrix:

    T(t) =
        [c^2  s^2        2*s*c ]
        [s^2  c^2        -2*s*c]
        [s*c (c^2 - s^2) s*c   ]

    where c = cos(t) and s = sin(t)

    The strain in the local ply frame is computed from the transposed
    transformation matrix (because we use engineering shear strain rather than
    tensorial shear strain):

    e_local = T(t)^T * e_global

    The transformation from the local ply frame to the global frame can be
    computed by simply replacing t with -t in the transformation matrix:

    s_global = T(-t) * s_local

    e_global = T(-t)^T * e_local

    Sensitivities w.r.t local or global stresses and strains can be transformed
    using the transpose of the relevant transformation matrix, e.g:

    e_local = T(t)^T * e_global
    e_global = T(-t)^T * e_local
    df/de_local = T(t) * df/de_global
    df/de_global = T(-t) * df/de_local
  */

  // Generic transformations
  // -----------------------
  /**
   * @brief Compute out = T(angle) * in
   *
   * @param angle Transformation angle
   * @param in 3-component input vector
   * @param out 3-component output vector
   */
  static void applyTransform(const TacsScalar angle, const TacsScalar in[],
                             TacsScalar out[]);

  /**
   * @brief Compute out = T(angle)^T * in
   *
   * @param angle Transformation angle
   * @param in 3-component input vector
   * @param out 3-component output vector
   */
  static void applyTransformTranspose(const TacsScalar angle,
                                      const TacsScalar in[], TacsScalar out[]);

  // Stress transformations
  // ----------------------
  /**
   * @brief Given the stress in the global frame, determine the stress in the
   * local ply frame
   *
   * @param angle Angle between ply and global frame in radians
   * @param global Stress in the global frame (s11, s22, s12)
   * @param plyStress Stress in the local ply frame (s1, s2, s12)
   */
  static void transformStressGlobal2Ply(const TacsScalar angle,
                                        const TacsScalar global[],
                                        TacsScalar plyStress[]) {
    applyTransform(angle, global, plyStress);
  }

  /**
   * @brief Given the stress in the local ply frame, determine the stress in the
   * global frame
   *
   * @param angle Angle between ply and global frame in radians
   * @param plyStress Stress in the local ply frame (s1, s2, s12)
   * @param global Stress in the global frame (s11, s22, s12)
   */
  static void transformStressPly2Global(const TacsScalar angle,
                                        const TacsScalar plyStress[],
                                        TacsScalar global[]) {
    applyTransform(-angle, plyStress, global);
  }

  /**
   * @brief Compute the derivative of the ply-frame stress with respect to the
   * ply angle
   *
   * @param angle Angle between ply and global frame in radians
   * @param global Stress in the global frame (s11, s22, s12)
   * @param plyStress Derivative of the ply-frame stress with respect to the ply
   * angle (ds1/dt, ds2/dt, ds12/dt)
   */
  static void transformStressGlobal2PlyAngleSens(const TacsScalar angle,
                                                 const TacsScalar global[],
                                                 TacsScalar plyStress[]);

  /**
   * @brief Compute the derivative of the global-frame stress with respect to
   * the ply angle
   *
   * @param angle Angle between ply and global frame in radians
   * @param plyStress Stress in the local ply frame (s1, s2, s12)
   * @param global Derivative of the global-frame stress with respect to the ply
   * angle (ds11/dt, ds22/dt, ds12/dt)
   */
  static void transformStressPly2GlobalAngleSens(const TacsScalar angle,
                                                 const TacsScalar plyStress[],
                                                 TacsScalar global[]);

  // Strain transformations
  // ----------------------
  /**
   * @brief Given the strain in the global frame, determine the strain in the
   * local ply frame
   *
   * @param angle Angle between ply and global frame in radians
   * @param global Strain in the global frame (e11, e22, gamma12)
   * @param plyStrain Strain in the local ply frame (e1, e2, gamma12)
   */
  static void transformStrainGlobal2Ply(const TacsScalar angle,
                                        const TacsScalar global[],
                                        TacsScalar plyStrain[]) {
    applyTransformTranspose(angle, global, plyStrain);
  }

  /**
   * @brief Given the strain in the local ply frame, determine the strain in the
   * global frame
   *
   * @param angle Angle between ply and global frame in radians
   * @param plyStrain Strain in the local ply frame (e1, e2, gamma12)
   * @param global Strain in the global frame (e11, e22, gamma12)
   */
  static void transformStrainPly2Global(const TacsScalar angle,
                                        const TacsScalar plyStrain[],
                                        TacsScalar global[]) {
    applyTransformTranspose(-angle, plyStrain, global);
  }

  /**
   * @brief Compute the derivative of the ply-frame strain with respect to the
   * ply angle
   *
   * @param angle Angle between ply and global frame in radians
   * @param global Strain in the global frame (e11, e22, gamma12)
   * @param plyStrain Derivative of the ply-frame strain with respect to the ply
   * angle (de1/dt, de2/dt, dgamma12/dt)
   */
  static void transformStrainGlobal2PlyAngleSens(const TacsScalar angle,
                                                 const TacsScalar global[],
                                                 TacsScalar plyStrain[]);

  /**
   * @brief Compute the derivative of the global-frame strain with respect to
   * the ply angle
   *
   * @param angle Angle between ply and global frame in radians
   * @param plyStrain Strain in the local ply frame (e1, e2, gamma12)
   * @param global Derivative of the global-frame strain with respect to the ply
   * angle (de11/dt, de22/dt, dgamma12/dt)
   */
  static void transformStrainPly2GlobalAngleSens(const TacsScalar angle,
                                                 const TacsScalar plyStrain[],
                                                 TacsScalar global[]);

  // Strain sensitivity transformations
  // ----------------------------------
  /**
   * @brief Transform a sensitivity w.r.t the global strain to a sensitivity
   * w.r.t the ply strain
   *
   * @param angle Angle between ply and global frame in radians
   * @param globalStrainSens Sensitivity w.r.t. strain in the global frame
   * (df/de11, df/de22, df/dgamma12)
   * @param plyStrainSens Sensitivity w.r.t. strain in the local ply frame
   * (df/de1, df/de2, df/dgamma12)
   */
  static void transformStrainSensGlobal2Ply(const TacsScalar angle,
                                            const TacsScalar globalStrainSens[],
                                            TacsScalar plyStrainSens[]) {
    applyTransform(angle, globalStrainSens, plyStrainSens);
  }

  /**
   * @brief Transform a sensitivity w.r.t the ply strain to a sensitivity w.r.t
   * the global strain
   *
   * @param angle Angle between ply and global frame in radians
   * @param plyStrainSens Sensitivity w.r.t. strain in the local ply frame
   * (df/de1, df/de2, df/dgamma12)
   * @param globalStrainSens Sensitivity w.r.t. strain in the global frame
   * (df/de11, df/de22, df/dgamma12)
   */
  static void transformStrainSensPly2Global(const TacsScalar angle,
                                            const TacsScalar plyStrainSens[],
                                            TacsScalar globalStrainSens[]) {
    applyTransform(-angle, plyStrainSens, globalStrainSens);
  }

  // Get or print other information
  // ------------------------------
  void testFailSens(double dh, TacsScalar angle);
  void printProperties();

  const char *getObjectName();

 private:
  static const char *name;

  TACSMaterialProperties *properties;

  // The stiffness properties
  TacsScalar Q11, Q12, Q22, Q44, Q55, Q66;
  TacsScalar C12, C16, C26, C66;

  TacsScalar plyThickness;
  TacsScalar rho;
  TacsScalar E1, E2;
  TacsScalar nu12, nu21;
  TacsScalar G12, G23, G13;

  // Keep track of which failure criterion to use
  CompositeFailureCriterion failureCriterion;

  // The stress-based strength properties
  TacsScalar Xt, Xc;
  TacsScalar Yt, Yc;
  TacsScalar S12;
  TacsScalar C;

  // The coefficients for the Tsai-Wu failure criterion
  TacsScalar F1, F2;
  TacsScalar F11, F12, F22, F66;

  // The material properties for the Cuntze failure criterion
  TacsScalar b_tt, b_tl;
  TacsScalar muWF, mu3W, mu3F;
  TacsScalar m;

  double Cuntze_IFF3_ksWeight;

  // The strain-based strength properties
  TacsScalar eXt, eXc;
  TacsScalar eYt, eYc;
  TacsScalar eS12;

  // The KS weight for the maximum-strain failure criterion
  TacsScalar ksWeight;
};

#endif  // TACS_MATERIAL_PROPERTIES_H
