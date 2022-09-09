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
  TACSMaterialProperties(TacsScalar _rho, TacsScalar _specific_heat,
                         TacsScalar _E1, TacsScalar _E2, TacsScalar _E3,
                         TacsScalar _nu12, TacsScalar _nu13, TacsScalar _nu23,
                         TacsScalar _G12, TacsScalar _G13, TacsScalar _G23,
                         TacsScalar _T1 = 0.0, TacsScalar _C1 = 0.0,
                         TacsScalar _T2 = 0.0, TacsScalar _C2 = 0.0,
                         TacsScalar _T3 = 0.0, TacsScalar _C3 = 0.0,
                         TacsScalar _S12 = 0.0, TacsScalar _S13 = 0.0,
                         TacsScalar _S23 = 0.0, TacsScalar _alpha1 = 0.0,
                         TacsScalar _alpha2 = 0.0, TacsScalar _alpha3 = 0.0,
                         TacsScalar _kappa1 = 0.0, TacsScalar _kappa2 = 0.0,
                         TacsScalar _kappa3 = 0.0);
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
                             TacsScalar *_S23);
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
  TacsScalar T1, C1, T2, C2, T3, C3, S12, S13, S23;

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
  based on a Tsai-Wu failure criterion or a smoothed maximum strain
  failure criterion, where the smoothing is performed using a KS
  function.

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

  void setKSWeight(TacsScalar _ksWeight);
  void setUseMaxStrainCriterion();
  void setUseTsaiWuCriterion();

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

  // Transform the stress and strain between global/local frames
  // -----------------------------------------------------------
  void getPlyStress(const TacsScalar strain[], TacsScalar stress[]);
  void transformStressGlobal2Ply(TacsScalar angle, const TacsScalar global[],
                                 TacsScalar plyStress[]);
  void transformStressPly2Global(TacsScalar angle, const TacsScalar plyStress[],
                                 TacsScalar global[]);
  void transformStressGlobal2PlyAngleSens(TacsScalar angle,
                                          const TacsScalar global[],
                                          TacsScalar plyStress[]);
  void transformStressPly2GlobalAngleSens(TacsScalar angle,
                                          const TacsScalar plyStress[],
                                          TacsScalar global[]);
  void transformStrainGlobal2Ply(TacsScalar angle, const TacsScalar global[],
                                 TacsScalar plyStrain[]);
  void transformStrainPly2Global(TacsScalar angle, const TacsScalar plyStrain[],
                                 TacsScalar global[]);
  void transformStrainGlobal2PlyAngleSens(TacsScalar angle,
                                          const TacsScalar global[],
                                          TacsScalar plyStrain[]);
  void transformStrainPly2GlobalAngleSens(TacsScalar angle,
                                          const TacsScalar plyStrain[],
                                          TacsScalar global[]);

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
  int useTsaiWuCriterion;

  // The stress-based strength properties
  TacsScalar Xt, Xc;
  TacsScalar Yt, Yc;
  TacsScalar S12;
  TacsScalar C;

  // The coefficients for the Tsai-Wu failure criterion
  TacsScalar F1, F2;
  TacsScalar F11, F12, F22, F66;

  // The strain-based strength properties
  TacsScalar eXt, eXc;
  TacsScalar eYt, eYc;
  TacsScalar eS12;

  // The KS weight for the maximum-strain failure criterion
  TacsScalar ksWeight;
};

#endif  // TACS_MATERIAL_PROPERTIES_H
