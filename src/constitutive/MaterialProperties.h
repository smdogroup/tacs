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
  acceptible bounds - these are tested during initialization.
*/
class OrthoPly : public TACSObject {
 public:
  // Create OrthoPly with a full set of orthotropic material relationships
  // ---------------------------------------------------------------------
  OrthoPly( TacsScalar _plyThickness, TacsScalar _rho, 
	    TacsScalar _E1, TacsScalar _E2, TacsScalar _nu12, 
	    TacsScalar _G12, TacsScalar _G23, TacsScalar _G13,
	    TacsScalar _Xt, TacsScalar _Xc, TacsScalar _Yt, TacsScalar _Yc, 
	    TacsScalar _S12, TacsScalar C = 0.0 );

  // Create OrthoPly with an equivalent von Mises relationships
  // ----------------------------------------------------------
  OrthoPly( TacsScalar _plyThickness, TacsScalar _rho, 
            TacsScalar E, TacsScalar nu, TacsScalar ys );

  void setKSWeight( TacsScalar _ksWeight );
  void setUseMaxStrainCriterion();
  void setUseTsaiWuCriterion();

  // Retrieve the material properties
  // --------------------------------
  TacsScalar getRho();
  TacsScalar getPlyThickness();
  void getStiffness( TacsScalar *_E1, TacsScalar *_E2, TacsScalar *_nu12,
                     TacsScalar *_G12, TacsScalar *_G23, TacsScalar *_G13 );
  void getLaminateStiffness( TacsScalar *_Q11, TacsScalar *_Q12, 
                             TacsScalar *_Q22, TacsScalar *_Q44,
                             TacsScalar *_Q55, TacsScalar *_Q66 );
  void getStrength( TacsScalar *_Xt, TacsScalar *_Xc, TacsScalar *_Yt, 
                    TacsScalar *_Yc, TacsScalar *_S12, TacsScalar *_C );
  void getStrainStrength( TacsScalar *_eXt, TacsScalar *_eXc, TacsScalar *_eYt, 
                          TacsScalar *_eYc, TacsScalar *_eS12 );
  void getTsaiWu( TacsScalar *_F1, TacsScalar *_F2, 
                  TacsScalar *_F11, TacsScalar *_F12, TacsScalar *_F22,
                  TacsScalar *_F66 );
  void getLaminateInvariants( TacsScalar *U1, TacsScalar *U2, 
			      TacsScalar *U3, TacsScalar *U4,
			      TacsScalar *U5, TacsScalar *U6 );

  // Calculate the Abar and Qbar matrices and their derivatives
  // Taken from Jones, Mechanics of composite materials pg. 51
  // ----------------------------------------------------------
  void calculateAbar( TacsScalar Abar[], TacsScalar angle );
  void calculateAbarAngleSens( TacsScalar Abar[], TacsScalar angle );
  void calculateQbar( TacsScalar Qbar[], TacsScalar angle );
  void calculateQbarAngleSens( TacsScalar Qbar[], TacsScalar angle );

  // Given the strain in the laminate coordinates, determine the 
  // stress in the laminate coordinates
  // -----------------------------------------------------------
  void calculateStress( TacsScalar stress[], const TacsScalar strain[], 
			TacsScalar angle );

  // Given the strain in the laminate coordinates, determine the failure load
  // ------------------------------------------------------------------------
  TacsScalar failure( TacsScalar angle, 
                      const TacsScalar strain[] );
  TacsScalar failureStrainSens( TacsScalar sens[], 
                                TacsScalar angle, 
                                const TacsScalar strain[] );
  TacsScalar failureAngleSens( TacsScalar * failSens, 
                               TacsScalar angle, 
                               const TacsScalar strain[] );

  // Calculate the failure load fraction for given 
  // constant and linear strain components
  // ---------------------------------------------
  TacsScalar calculateFailLoad( TacsScalar angle, 
                                const TacsScalar cstrain[], 
                                const TacsScalar lstrain[] );
  TacsScalar calculateFailLoadStrainSens( TacsScalar cSens[], 
                                          TacsScalar lSens[], 
                                          TacsScalar angle, 
                                          const TacsScalar cstrain[], 
                                          const TacsScalar lstrain[] );
  TacsScalar calculateFailLoadAngleSens( TacsScalar * posSens, 
                                         TacsScalar angle, 
                                         const TacsScalar cstrain[], 
                                         const TacsScalar lstrain[] );

  // Transform the stress and strain between global/local frames
  // -----------------------------------------------------------
  void getPlyStress( TacsScalar stress[], const TacsScalar strain[] );
  void transformStressGlobal2Ply( TacsScalar plyStress[], 
				  const TacsScalar global[], 
				  TacsScalar angle );
  void transformStressPly2Global( TacsScalar global[], 
				  const TacsScalar plyStress[], 
				  TacsScalar angle );
  void transformStressGlobal2PlyAngleSens( TacsScalar plyStress[], 
					   const TacsScalar global[], 
					   TacsScalar angle );
  void transformStressPly2GlobalAngleSens( TacsScalar global[], 
					   const TacsScalar plyStress[], 
					   TacsScalar angle );
  void transformStrainGlobal2Ply( TacsScalar plyStrain[], 
				  const TacsScalar global[], 
				  TacsScalar angle );
  void transformStrainPly2Global( TacsScalar global[], 
				  const TacsScalar plyStrain[], 
				  TacsScalar angle );
  void transformStrainGlobal2PlyAngleSens( TacsScalar plyStrain[], 
					   const TacsScalar global[], 
					   TacsScalar angle );
  void transformStrainPly2GlobalAngleSens( TacsScalar global[], 
					   const TacsScalar plyStrain[], 
					   TacsScalar angle );
  // Get or print other information
  // ------------------------------
  void testFailSens( double dh, TacsScalar angle );
  void printProperties();

  const char * TACSObjectName();

 private:
  static const char * name;

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

#endif
