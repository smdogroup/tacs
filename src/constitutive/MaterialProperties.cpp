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

#include "MaterialProperties.h"
#include "YSlibrary.h"

const char * OrthoPly::name = "OrthoPly";

const static double LINEAR_STRESS_CUTOFF = 1e-15;
const static double HUGE_FAILURE_LOAD = 1e20;

/*
  The OrthoPly material class. 

  This uses a Tsai-Wu tensor failure criteria. 

  Inputs are:
  E1, E2, nu12, G12 = In-plane stiffness properties
  G23, G13 = Out of plane shear stiffness

  Xt, Xc = Tensile and compressive fiber-aligned failure loads
  Yt, Yc = Tensile and compressive off-axis failure loads
  S12 = In plane shear failure load
  C = Interaction strength such that sigma_1 = sigma_2 = C
*/
OrthoPly::OrthoPly( TacsScalar _plyThickness, TacsScalar _rho, 
                    TacsScalar _E1, TacsScalar _E2, TacsScalar _nu12, 
                    TacsScalar _G12, TacsScalar _G23, TacsScalar _G13,
                    TacsScalar _Xt, TacsScalar _Xc, 
                    TacsScalar _Yt, TacsScalar _Yc, 
                    TacsScalar _S12, TacsScalar _C ){
  plyThickness = _plyThickness;
  rho = _rho;
  E1 = _E1;
  E2 = _E2;
  nu12 = _nu12;
  G12 = _G12;
  G23 = _G23;
  G13 =_G13;
  
  nu21 = nu12*E2/E1;
  
  Q11 = E1/(1.0 - nu12*nu21);
  Q22 = E2/(1.0 - nu12*nu21);
  Q12 = nu12*E2/(1.0 - nu12*nu21);
  
  Q44 = G23;
  Q55 = G13;
  Q66 = G12;

  // Definitions from Jones, Mechanics of Composite Materials 
  C12 = (Q11 + Q22 - 4.0*Q66);
  C16 = (Q11 - Q12 - 2.0*Q66);
  C26 = (Q12 - Q22 + 2.0*Q66);
  C66 = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66);

  // Record the failure data. If the coefficients
  // are negative, make them positive
  Xt = _Xt;
  Xc = _Xc;
  if (TacsRealPart(Xc) < 0.0){ 
    Xc *= -1.0; 
  }

  Yt = _Yt;
  Yc = _Yc;
  if (TacsRealPart(Yc) < 0.0){ 
    Yc *= -1.0; 
  }

  S12 = _S12;
  C = _C;

  // Determine the strain strength values based on the 
  // supplied stress strength values
  eXt = Xt/E1;
  eXc = Xc/E1;

  eYt = Yt/E2;
  eYc = Yc/E2;

  eS12 = S12/G12;

  // By default, use Tsai-Wu
  useTsaiWuCriterion = 1;

  // Compute the coefficients for the Tsai-Wu failure criteria
  F1 = (Xc - Xt)/(Xt*Xc);
  F2 = (Yc - Yt)/(Yt*Yc);
  
  F11 = 1.0/(Xt*Xc);
  F22 = 1.0/(Yt*Yc);

  F66 = 1.0/(S12*S12);
  if (TacsRealPart(C) != 0.0){
    F12 = 0.5*(1.0 - (F1 + F2)*C - (F11 + F22)*C*C)/(C*C);
  }
  else {
    F12 = 0.0; // Assume no interaction
  }

  // Check the stability criterion
  if (TacsRealPart(F12*F12) >= TacsRealPart(F11*F22)){
    fprintf(stderr, "OrthoPly: Value of C = %e results in non-physical \
F12 = %e. Setting F12 = 0.\n", TacsRealPart(C), TacsRealPart(F12));
    fprintf(stderr, "OrthoPly: Tsai-Wu coefficients: F11: %e, F22: %e, \
F66: %e, F1: %e, F2: %e\n", TacsRealPart(F11), TacsRealPart(F22), 
            TacsRealPart(F66), TacsRealPart(F1), TacsRealPart(F2));
    F12 = 0.0;      
  }

  // Set the default value of the KS penalty
  ksWeight = 100.0;
}

/*
  Set the orthotropic material properties with an isotropic 
  material and
*/
OrthoPly::OrthoPly( TacsScalar _plyThickness, TacsScalar _rho, 
                    TacsScalar E, TacsScalar nu, TacsScalar ys ){
  plyThickness = _plyThickness;
  rho = _rho;
  E1 = E;
  E2 = E;
  nu12 = nu;
  G12 = 0.5*E/(1.0 + nu);
  G23 = G12;
  G13 = G12;
  
  nu21 = nu12*E2/E1;
  
  Q11 = E1/(1.0 - nu12*nu21);
  Q22 = E2/(1.0 - nu12*nu21);
  Q12 = nu12*E2/(1.0 - nu12*nu21);
  
  Q44 = G23;
  Q55 = G13;
  Q66 = G12;

  // Definitions from Jones, Mechanics of Composite Materials 
  C12 = (Q11 + Q22 - 4.0*Q66);
  C16 = (Q11 - Q12 - 2.0*Q66);
  C26 = (Q12 - Q22 + 2.0*Q66);
  C66 = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66);

  // Record the failure data
  Xt = Xc = ys;
  Yt = Yc = ys;
  S12 = ys/sqrt(3);
  C = ys;

  // Determine the strain strength values based on the 
  // supplied stress strength values
  eXt = Xt/E1;
  eXc = Xc/E1;

  eYt = Yt/E2;
  eYc = Yc/E2;

  eS12 = S12/G12;

  // By default, use Tsai-Wu
  useTsaiWuCriterion = 1;

  // Set the coefficients for the Tsai-Wu failure criteria so that they
  // match the plane-stress von Mises failure criteria
  F1 = 0.0;
  F2 = 0.0;
  
  F11 = 1.0/(ys*ys);
  F22 = 1.0/(ys*ys);

  F66 = 3.0/(ys*ys);
  F12 = -0.5/(ys*ys);

  // Set the default value of the KS penalty
  ksWeight = 100.0;
}

/*
  Set the KS penalty factor
*/
void OrthoPly::setKSWeight( TacsScalar _ksWeight ){
  ksWeight = _ksWeight;
}

/*
  Set the failure criteria to use
*/
void OrthoPly::setUseMaxStrainCriterion(){
  useTsaiWuCriterion = 0;
}

void OrthoPly::setUseTsaiWuCriterion(){
  useTsaiWuCriterion = 1;
}

/*
  Get the density of the material
*/
TacsScalar OrthoPly::getRho(){ return rho; }

/*
  Get the thickness of a single laminae 
*/
TacsScalar OrthoPly::getPlyThickness(){ return plyThickness; }

/*
  Get the stiffness constants
*/
void OrthoPly::getStiffness( TacsScalar *_E1, TacsScalar *_E2,
                             TacsScalar *_nu12, TacsScalar *_G12, 
                             TacsScalar *_G23, TacsScalar *_G13 ){
  *_E1 = E1;
  *_E2 = E2;
  *_nu12 = nu12;
  *_G12 = G12;
  *_G23 = G23;
  *_G13 = G13;  
}

/*
  Get the lamination stiffness objects
*/
void OrthoPly::getLaminateStiffness( TacsScalar *_Q11, TacsScalar *_Q12, 
                                     TacsScalar *_Q22, TacsScalar *_Q44,
                                     TacsScalar *_Q55, TacsScalar *_Q66 ){
  *_Q11 = Q11;
  *_Q12 = Q12;
  *_Q22 = Q22;
  *_Q44 = Q44;
  *_Q55 = Q55;
  *_Q66 = Q66;
}

/*
  Get the strength parameters
*/
void OrthoPly::getStrength( TacsScalar *_Xt, TacsScalar *_Xc, TacsScalar *_Yt, 
                            TacsScalar *_Yc, TacsScalar *_S12, TacsScalar *_C ){
  *_Xt = Xt;
  *_Xc = Xc;
  *_Yt = Yt;
  *_Yc = Yc;
  *_S12 = S12;
  *_C = C;
}

/*
  Get the strength parameters for strain
*/
void OrthoPly::getStrainStrength( TacsScalar *_eXt, TacsScalar *_eXc, 
                                  TacsScalar *_eYt, TacsScalar *_eYc, 
                                  TacsScalar *_eS12 ){
  *_eXt = eXt;
  *_eXc = eXc;
  *_eYt = eYt;
  *_eYc = eYc;
  *_eS12 = eS12;
}

/*
  Get the Tsai-Wu failure constants
*/
void OrthoPly::getTsaiWu( TacsScalar *_F1, TacsScalar *_F2, 
                          TacsScalar *_F11, TacsScalar *_F12, TacsScalar *_F22,
                          TacsScalar *_F66 ){
  *_F1 = F1;
  *_F2 = F2;
  *_F11 = F11;
  *_F12 = F12;
  *_F22 = F22;
  *_F66 = F66;
}

/*
  Retrieve the stiffness invariants for the laminate
*/
void OrthoPly::getLaminateInvariants( TacsScalar *U1, TacsScalar *U2, 
                                      TacsScalar *U3, TacsScalar *U4,
                                      TacsScalar *U5, TacsScalar *U6 ){
  *U1 = 0.125*(3.0*Q11 + 3.0*Q22 + 2.0*Q12 + 4.0*Q66);
  *U2 = 0.5*(Q11 - Q22);
  *U3 = 0.125*(Q11 + Q22 - 2.0*Q12 - 4.0*Q66);
  *U4 = 0.125*(Q11 + Q22 + 6.0*Q12 - 4.0*Q66);
  *U5 = 0.5*(Q44 + Q55);
  *U6 = 0.5*(Q44 - Q55);
}

/*!
  Calculate the in-plane stresses in the ply-coordinate axis
*/
void OrthoPly::getPlyStress( TacsScalar stress[], const TacsScalar strain[] ){
  stress[0] = Q11*strain[0] + Q12*strain[1];
  stress[1] = Q12*strain[0] + Q22*strain[1];
  stress[2] = Q66*strain[2];
}

/*
  Compute the rotation
*/
void OrthoPly::calculateAbar( TacsScalar Abar[], TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);
  
  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;

  Abar[0] = cos2*Q44 + sin2*Q55;
  Abar[1] = cos1*sin1*(Q55 - Q44);
  Abar[2] = sin2*Q44 + cos2*Q55;
}

/*
  Calculate the derivative of the Abar matrix w.r.t. the angle 
*/
void OrthoPly::calculateAbarAngleSens( TacsScalar Abar[], TacsScalar angle ){
  TacsScalar cos_2 = cos(2.0*angle);
  TacsScalar sin_2 = sin(2.0*angle);

  Abar[0] = - sin_2*Q44 + sin_2*Q55;
  Abar[1] = cos_2*(Q55 - Q44);
  Abar[2] = sin_2*Q44 - sin_2*Q55;
}

/*
  Taken from Jones, Mechanics of composite materials pg. 51
*/
void OrthoPly::calculateQbar( TacsScalar Qbar[], TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);
  
  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;
  
  TacsScalar cos4 = cos2*cos2;
  TacsScalar sin4 = sin2*sin2;
    
  // First row of the matrix
  Qbar[0] = Q11*cos4 + 2.0*(Q12 + 2.0*Q66 )*sin2*cos2 + Q22*sin4;
  Qbar[1] = C12*sin2*cos2 + Q12*(sin4 + cos4 );
  Qbar[2] = C16*sin1*cos2*cos1 + C26*sin2*sin1*cos1;
  
  // Second row of the matrix    
  Qbar[3] = Q11*sin4 + 2.0*(Q12 + 2.0*Q66 )*sin2*cos2 + Q22*cos4;
  Qbar[4] = C16*sin2*sin1*cos1 + C26*sin1*cos2*cos1;
  
  // Last row of the matrix
  Qbar[5] = C66*sin2*cos2 + Q66*(sin4 + cos4 );
} 

/*
  Calculate the sensitivity of the Qbar matrix w.r.t. the angle
*/
void OrthoPly::calculateQbarAngleSens( TacsScalar Qbar[], TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;

  TacsScalar s_cos2 = - 2.0*cos1*sin1;
  TacsScalar s_sin2 =   2.0*sin1*cos1;
  
  TacsScalar s_cos4 = - 4.0*cos2*cos1*sin1;
  TacsScalar s_sin4 =   4.0*sin2*sin1*cos1;
  
  TacsScalar s_sin2cos2 = (sin2*s_cos2 + s_sin2*cos2);
  
  TacsScalar s_sin3cos1 =   3.0*sin2*cos2 - sin2*sin2;
  TacsScalar s_cos3sin1 = - 3.0*cos2*sin2 + cos2*cos2;
  
  // First row of the matrix
  Qbar[0] = Q11*s_cos4 + 2.0*(Q12 + 2.0*Q66)*s_sin2cos2 +  Q22*s_sin4;
  Qbar[1] = C12*s_sin2cos2 + Q12*(s_sin4 + s_cos4);  
  Qbar[2] = C16*s_cos3sin1 + C26*s_sin3cos1;
  
  // Second row of the matrix    
  Qbar[3] = Q11*s_sin4 + 2.0*(Q12 + 2.0*Q66)*s_sin2cos2 + Q22*s_cos4;
  Qbar[4] = C16*s_sin3cos1 + C26*s_cos3sin1;
  
  // Last row of the matrix
  Qbar[5] = C66*s_sin2cos2 + Q66*(s_sin4 + s_cos4 );
} 

void OrthoPly::calculateStress( TacsScalar stress[], const TacsScalar strain[], 
                                TacsScalar angle ){
  TacsScalar strn[3], strs[3];
  transformStrainGlobal2Ply(strn, strain, angle);
  getPlyStress(strs, strn);
  transformStressPly2Global(stress, strs, angle);
}

/*
  Calculate the failure criteria

  fail (the return value) is such that
  
  fail <  1.0 ==> Okay
  fail >= 1.0 ==> Material failure

  One of two failure criteria are used:

  1. The Tsai-Wu failure criterion:

  F(s) = 
  F1*s[0] + F2*s[1]
  F11*s[0]**2 + F22*s[1]**2 + 
  2.0*F12*s[0]*s[1] + F66*s[2]**2 <= 1.0

  2. The maximum strain failure criteria:

  KS(e_{i}/e_max^{+/-}, ksWeight) <= 1.0

  where e_max^{+/-} are the postive and negative strains at failure
  for component i. This ensures that the ratio of the maximum strain
  over the strain at failure is below 1.
*/
TacsScalar OrthoPly::failure( TacsScalar angle, 
                              const TacsScalar strain[] ){
  TacsScalar e[3]; // Ply strain
  transformStrainGlobal2Ply(e, strain, angle);
  
  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion){
    TacsScalar s[3]; // Ply stress
    getPlyStress(s, e);
    
    fail = (F11*s[0]*s[0] + F22*s[1]*s[1] + 2.0*F12*s[0]*s[1] + 
            F66*s[2]*s[2] + F1*s[0] + F2*s[1]);
  }
  else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6];
    f[0] = e[0]/eXt;   f[1] = -e[0]/eXc;
    f[2] = e[1]/eYt;   f[3] = -e[1]/eYc;
    f[4] = e[2]/eS12;  f[5] = -e[2]/eS12;

    TacsScalar max = f[0];
    for ( int k = 1; k < 6; k++ ){
      if (TacsRealPart(f[k]) > TacsRealPart(max)){ 
        max = f[k]; 
      }
    }

    TacsScalar ksSum = 0.0;
    for ( int k = 0; k < 6; k++ ){
      ksSum += exp(ksWeight*(f[k] - max));
    }
    fail = max + log(ksSum)/ksWeight;
  }

  return fail;
}

TacsScalar OrthoPly::failureStrainSens( TacsScalar sens[], 
                                        TacsScalar angle, 
                                        const TacsScalar strain[] ){
  TacsScalar e[3]; // Ply strain
  transformStrainGlobal2Ply(e, strain, angle);

  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion){
    TacsScalar s[3]; // Ply stress
    getPlyStress(s, e);
    
     fail = (F11*s[0]*s[0] + F22*s[1]*s[1] + 2.0*F12*s[0]*s[1] + 
             F66*s[2]*s[2] + F1*s[0] + F2*s[1]);
    
    sens[0] = F1 + 2.0*F11*s[0] + 2.0*F12*s[1];
    sens[1] = F2 + 2.0*F22*s[1] + 2.0*F12*s[0];
    sens[2] = 2.0*F66*s[2];

    TacsScalar sSens[3];
    getPlyStress(sSens, sens);
    
    transformStressPly2Global(sens, sSens, angle);
  }
  else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6], fexp[6];
    f[0] = e[0]/eXt;   f[1] = -e[0]/eXc;
    f[2] = e[1]/eYt;   f[3] = -e[1]/eYc;
    f[4] = e[2]/eS12;  f[5] = -e[2]/eS12;

    TacsScalar max = f[0];
    for ( int k = 1; k < 6; k++ ){
      if (TacsRealPart(f[k]) > TacsRealPart(max)){
        max = f[k]; 
      }
    }

    TacsScalar ksSum = 0.0;
    for ( int k = 0; k < 6; k++ ){
      fexp[k] = exp(ksWeight*(f[k] - max));
      ksSum += fexp[k];
    }

    fail = max + log(ksSum)/ksWeight;

    TacsScalar sSens[3];
    sSens[0] = (fexp[0]*eXc - fexp[1]*eXt)/(ksSum*eXt*eXc);
    sSens[1] = (fexp[2]*eYc - fexp[3]*eYt)/(ksSum*eYt*eYc);
    sSens[2] = (fexp[4] - fexp[5])/(ksSum*eS12);

    transformStressPly2Global(sens, sSens, angle);
  }

  return fail;
}

TacsScalar OrthoPly::failureAngleSens( TacsScalar * failSens, 
                                       TacsScalar angle, 
                                       const TacsScalar strain[] ){
 
  TacsScalar e[3], se[3]; // The ply strain
  transformStrainGlobal2Ply(e, strain, angle);

  TacsScalar fail = 0.0;
  if (useTsaiWuCriterion){
    TacsScalar s[3]; // The ply stress
    getPlyStress(s, e);
    
    fail = (F11*s[0]*s[0] + F22*s[1]*s[1] + 2.0*F12*s[0]*s[1] + 
            F66*s[2]*s[2] + F1*s[0] + F2*s[1]);
    
    // Compute the sensitivity of the transformation
    transformStrainGlobal2PlyAngleSens(se, strain, angle);
    
    // Compute the sensitivity of the stress
    TacsScalar ss[3];
    getPlyStress(ss, se);
    
    *failSens = (2.0*(F11*s[0]*ss[0] + F22*s[1]*ss[1] + 
                      F12*(ss[0]*s[1] + s[0]*ss[1]) +
                      F66*s[2]*ss[2]) + 
                 F1*ss[0] + F2*ss[1]);
  }
  else {
    // Calculate the values of each of the failure criteria
    TacsScalar f[6], fs[6];
    f[0] =  e[0]/eXt;   f[1] = -e[0]/eXc;
    f[2] =  e[1]/eYt;   f[3] = -e[1]/eYc;
    f[4] =  e[2]/eS12;  f[5] = -e[2]/eS12;

    // Compute the sensitivity of the transformation
    transformStrainGlobal2PlyAngleSens(se, strain, angle);
    fs[0] =  se[0]/eXt;   fs[1] = -se[0]/eXc;
    fs[2] =  se[1]/eYt;   fs[3] = -se[1]/eYc;
    fs[4] =  se[2]/eS12;  fs[5] = -se[2]/eS12;

    TacsScalar max = f[0];
    for ( int k = 1; k < 6; k++ ){
      if (TacsRealPart(f[k]) > TacsRealPart(max)){
        max = f[k]; 
      }
    }

    TacsScalar ksSum = 0.0;
    TacsScalar fSens = 0.0;
    for ( int k = 0; k < 6; k++ ){
      TacsScalar fexp = exp(ksWeight*(f[k] - max));
      fSens += fexp*fs[k];
      ksSum += fexp;
    }

    fail = max + log(ksSum)/ksWeight;

    *failSens = fSens/ksSum;
  }  

  return fail;
}

/*
  Calculate the failure load based on constant and linear components
  of the strain. The components of the strain are in the laminate
  coordinates (global coordinates) but just the in-plane
  contributions. The angle specified in radians is angle that the ply
  is oriented at - to get to the lamina coordinate system. See Jones -
  Mechanics of Composite Materials for details.

  returns:
  the multiple of lstrain at which failure occurs
  
  input:
  angle = [radians] the orientation angle of the lamina
  cstrain = the constant in-plane strain components
  lstrain = the linear in-plane strain components

  This works by computing F(cstrain, P*lstrain) = 0 for some P. This
  turns out to be a quadratic in P.
*/
TacsScalar OrthoPly::calculateFailLoad( TacsScalar angle, 
                                        const TacsScalar cstrain[], 
                                        const TacsScalar lstrain[] ){

  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3]; 
  transformStrainGlobal2Ply(cstn, cstrain, angle);
  transformStrainGlobal2Ply(lstn, lstrain, angle);
  
  TacsScalar cstr[3], lstr[3]; 
  getPlyStress(cstr, cstn);
  getPlyStress(lstr, lstn);

  TacsScalar c = (F1*cstr[0] + F2*cstr[1] + 
                  F11*cstr[0]*cstr[0] + F22*cstr[1]*cstr[1] + 
                  2.0*F12*cstr[0]*cstr[1] + 
                  F66*cstr[2]*cstr[2]) - 1.0;
  TacsScalar b = (F1*lstr[0] + F2*lstr[1] + 
                  2.0*F11*cstr[0]*lstr[0] + 2.0*F22*cstr[1]*lstr[1] + 
                  2.0*F12*(cstr[0]*lstr[1] + cstr[1]*lstr[0]) +
                  2.0*F66*cstr[2]*lstr[2]);
  TacsScalar a = (F11*lstr[0]*lstr[0] + F22*lstr[1]*lstr[1] + 
                  2.0*F12*lstr[0]*lstr[1] + 
                  F66*lstr[2]*lstr[2]);
  TacsScalar pos = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF*TacsRealPart(F11 + F22)){
    pos = HUGE_FAILURE_LOAD;
  }
  else if (TacsRealPart(c) >= 0.0){
    pos = 0.0;
  }
  else {
    TacsScalar discrim = b*b - 4.0*a*c;
    if (TacsRealPart(discrim) < 0.0){
      pos = 0.0;
    }
    else {
      discrim = sqrt(discrim);
      
      if (TacsRealPart(b) >= 0.0){
        pos = 2.0*c/(b + discrim);
      }
      else { // b < 0.0
        pos = 0.5*(-b + discrim)/a;
      }
    }
  }
  
  return pos;
}

/*
  Determine the sensitivity of the failure load calculation above to
  the strain components cstrain and lstrain. Here the arguments angle
  [radians], and the in-plane strain components cstrain[] and
  lstrain[] are the same as above. The output cSens[], lSens[] are
  the sensitivity of the failure load w.r.t. to the constant and linear
  strain components.
*/
TacsScalar OrthoPly::calculateFailLoadStrainSens( TacsScalar cSens[], 
                                                  TacsScalar lSens[], 
                                                  TacsScalar angle, 
                                                  const TacsScalar cstrain[], 
                                                  const TacsScalar lstrain[] ){
  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3]; 
  transformStrainGlobal2Ply(cstn, cstrain, angle);
  transformStrainGlobal2Ply(lstn, lstrain, angle);

  // The constant and linearly varying components of stress - in the ply axis
  TacsScalar cstr[3], lstr[3]; 
  getPlyStress(cstr, cstn);
  getPlyStress(lstr, lstn);
 
  TacsScalar c = (F1*cstr[0] + F2*cstr[1] + 
                  F11*cstr[0]*cstr[0] + F22*cstr[1]*cstr[1] + 
                  2.0*F12*cstr[0]*cstr[1] + 
                  F66*cstr[2]*cstr[2]) - 1.0;
  TacsScalar b = (F1*lstr[0] + F2*lstr[1] + 
                  2.0*F11*cstr[0]*lstr[0] + 2.0*F22*cstr[1]*lstr[1] + 
                  2.0*F12*(cstr[0]*lstr[1] + cstr[1]*lstr[0]) +
                  2.0*F66*cstr[2]*lstr[2]);  
  TacsScalar a = (F11*lstr[0]*lstr[0] + F22*lstr[1]*lstr[1] + 
                  2.0*F12*lstr[0]*lstr[1] + 
                  F66*lstr[2]*lstr[2]);

  TacsScalar pos = 0.0;
  TacsScalar pa = 0.0, pb = 0.0, pc = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF*TacsRealPart(F11 + F22)){
    pos = HUGE_FAILURE_LOAD;
  }
  else if (TacsRealPart(c) >= 0.0){
    pos = 0.0;
  }
  else {
    TacsScalar discrim = b*b - 4.0*a*c;
    if (TacsRealPart(discrim) < 0.0){
      pos = 0.0;
    }
    else {
      discrim = sqrt(discrim);      

      if (TacsRealPart(b) >= 0.0){
        pos = 2.0*c/(b + discrim);
      }
      else {
        pos = 0.5*(-b + discrim)/a;
      }

      pa = -(c + pos*discrim)/(discrim*a);
      pb = - pos/discrim;
      pc = - 1.0/discrim;
    }
  }
  
  // Now, multiply the derivative of dp/da, dp/db, dp/dc by 
  // da/dcstr, db/dcstr, 

  cSens[0] = (pc*(F1 + 2.0*F11*cstr[0] + 2.0*F12*cstr[1]) +
              pb*(2.0*F11*lstr[0] + 2.0*F12*lstr[1]));
  cSens[1] = (pc*(F2 + 2.0*F12*cstr[0] + 2.0*F22*cstr[1]) +
              pb*(2.0*F12*lstr[0] + 2.0*F22*lstr[1]));
  cSens[2] = (pc*(2.0*F66*cstr[2]) +
              pb*(2.0*F66*lstr[2]));

  lSens[0] = (pb*(F1 + 2.0*F11*cstr[0] + 2.0*F12*cstr[1]) +
              pa*(2.0*F11*lstr[0] + 2.0*F12*lstr[1]));
  lSens[1] = (pb*(F2 + 2.0*F12*cstr[1] + 2.0*F22*cstr[1]) +
              pa*(2.0*F12*lstr[0] + 2.0*F22*lstr[1]));
  lSens[2] = (pb*(2.0*F66*cstr[2]) + 
              pa*(2.0*F66*lstr[2]));
  
  TacsScalar cstrSens[3], lstrSens[3];
  getPlyStress(cstrSens, cSens);
  getPlyStress(lstrSens, lSens);

  transformStressPly2Global(cSens, cstrSens, angle);
  transformStressPly2Global(lSens, lstrSens, angle);

  return pos;
}

/*
  Determine the sensitivity of the failure load calculation to to ply
  angle. Here, the arguments are the same. The return value posSens is
  the sensitivity of the failure load w.r.t. the ply angle.
*/
TacsScalar OrthoPly::calculateFailLoadAngleSens( TacsScalar * posSens, 
                                                 TacsScalar angle, 
                                                 const TacsScalar cstrain[], 
                                                 const TacsScalar lstrain[] ){
 
  // The constant and linearly varying components of the strain
  TacsScalar cstn[3], lstn[3]; 
  transformStrainGlobal2Ply(cstn, cstrain, angle);
  transformStrainGlobal2Ply(lstn, lstrain, angle);
  
  // The constant and linearly varying components of stress - in the ply axis
  TacsScalar cstr[3], lstr[3];
  getPlyStress(cstr, cstn);
  getPlyStress(lstr, lstn);

  // Now, determine the sensitivity of the transformation
  transformStrainGlobal2PlyAngleSens(cstn, cstrain, angle);
  transformStrainGlobal2PlyAngleSens(lstn, lstrain, angle);

  // The constant and linearly varying sensitivites w.r.t. the strain
  TacsScalar scstr[3], slstr[3];
  getPlyStress(scstr, cstn);
  getPlyStress(slstr, lstn);

  TacsScalar c = (F1*cstr[0] + F2*cstr[1] + 
                  F11*cstr[0]*cstr[0] + F22*cstr[1]*cstr[1] + 
                  2.0*F12*cstr[0]*cstr[1] + 
                  F66*cstr[2]*cstr[2]) - 1.0;
  TacsScalar b = (F1*lstr[0] + F2*lstr[1] + 
                  2.0*F11*cstr[0]*lstr[0] + 2.0*F22*cstr[1]*lstr[1] + 
                  2.0*F12*(cstr[0]*lstr[1] + cstr[1]*lstr[0]) +
                  2.0*F66*cstr[2]*lstr[2]);
  TacsScalar a = (F11*lstr[0]*lstr[0] + F22*lstr[1]*lstr[1] + 
                  2.0*F12*lstr[0]*lstr[1] + 
                  F66*lstr[2]*lstr[2]);

  *posSens = 0.0;
  TacsScalar pos = 0.0;
  TacsScalar pa = 0.0, pb = 0.0, pc = 0.0;

  if (fabs(TacsRealPart(a)) < LINEAR_STRESS_CUTOFF*TacsRealPart(F11 + F22)){
    pos = HUGE_FAILURE_LOAD;
  }
  else if (TacsRealPart(c) >= 0.0){
    pos = 0.0;
  }
  else {
    TacsScalar discrim = b*b - 4.0*a*c;
    if (TacsRealPart(discrim) < 0.0){
      pos = 0.0;
    }
    else {
      discrim = sqrt(discrim);

      if (TacsRealPart(b) >= 0.0){
        pos = 2.0*c/(b + discrim);
      }
      else {
        pos = 0.5*(-b + discrim)/a;
      }

      pa = -(c + pos*discrim)/(discrim*a);
      pb = - pos/discrim;
      pc = - 1.0/discrim;
    }
  }

  TacsScalar sc = (F1*scstr[0] + F2*scstr[1] + 
                   2.0*F11*cstr[0]*scstr[0] + 2.0*F22*cstr[1]*scstr[1] + 
                   2.0*F12*(cstr[0]*scstr[1] + scstr[0]*cstr[1] ) +
                   2.0*F66*cstr[2]*scstr[2]);
  TacsScalar sb = (F1*slstr[0] + F2*slstr[1] + 
                   2.0*F11*(cstr[0]*slstr[0] + scstr[0]*lstr[0]) + 
                   2.0*F22*(cstr[1]*slstr[1] + scstr[1]*lstr[1]) +
                   2.0*F12*(cstr[0]*slstr[1] + scstr[1]*lstr[0] +
                            scstr[0]*lstr[1] + cstr[1]*slstr[0]) +
                   2.0*F66*(cstr[2]*slstr[2] + scstr[2]*lstr[2]));
  TacsScalar sa = (2.0*F11*lstr[0]*slstr[0] + 2.0*F22*lstr[1]*slstr[1] + 
                   2.0*F12*(lstr[0]*slstr[1] + slstr[0]*lstr[1]) +
                   2.0*F66*lstr[2]*slstr[2]);

  *posSens = pa*sa + pb*sb + pc*sc;
  return pos;
}

/*! 
  The following function tests the accuracy of the implementation
  of the sensitivities for the failure calculation and the sensitivity
  of the transformation equations 
*/
void OrthoPly::testFailSens( double dh, TacsScalar angle ){
  // Select various strain components ...
  TacsScalar strain[3];

  printf("\nTesting failure sensitivity for angle = %5.2f \n", 
         TacsRealPart(angle));

  for ( int k = 0; k < 3; k++ ){
    strain[k] = -1.0;

    // Calculate the failure load
    TacsScalar p = failure(angle,strain);
    printf("Failure criteria = %15.8e \n", TacsRealPart(p));

    // Calculate the sensitivity of the failure load
    TacsScalar sens[3];
    failureStrainSens(sens, angle, strain);

    // Compare the result to a finite-difference calculation
    for ( int j = 0; j < 3; j++ ){
      // Test the constant component
      TacsScalar val = strain[j];
      strain[j] = val + dh;
      TacsScalar p1 = failure(angle, strain);

      strain[j] = val - dh;
      TacsScalar p2 = failure(angle, strain);
      strain[j] = val;

      TacsScalar fd = 0.5*(p1 - p2)/dh;
      printf("sens[%d] FD: %15.8e An: %15.8e Error: %10.3e \n", 
             j, TacsRealPart(fd), TacsRealPart(sens[j]), 
             TacsRealPart((fd - sens[j])/sens[j]));
    }
    
    // Calculate the sensitivity w.r.t. the angle
    TacsScalar paSens;
    failureAngleSens(&paSens, angle, strain);
    TacsScalar p1 = failure(angle+dh, strain);
    TacsScalar p2 = failure(angle-dh, strain);
    TacsScalar fd = 0.5*(p1 - p2)/dh;

    printf("Angle sensitivity FD: %15.8e An: %15.8e Error: %10.3e \n", 
           TacsRealPart(fd), TacsRealPart(paSens), TacsRealPart((fd - paSens)/paSens));
  }
}

/*
  Print the stiffness and failure properties
*/
void OrthoPly::printProperties(){
  printf("\nStiffness properties \n");
  printf("E1   = %15.5e \n", TacsRealPart(E1));
  printf("E2   = %15.5e \n", TacsRealPart(E2));
  printf("nu12 = %15.5e \n", TacsRealPart(nu12));
  printf("nu21 = %15.5e \n", TacsRealPart(nu21));
  printf("G12  = %15.5e \n", TacsRealPart(G12));
  printf("G23  = %15.5e \n", TacsRealPart(G23));
  printf("G13  = %15.5e \n", TacsRealPart(G13));

  printf("\nFailure Properties \n");
  printf("Xt   = %15.5e \n", TacsRealPart(Xt));
  printf("Xc   = %15.5e \n", TacsRealPart(Xc));
  printf("Yt   = %15.5e \n", TacsRealPart(Yt));
  printf("Yc   = %15.5e \n", TacsRealPart(Yc));
  printf("S12  = %15.5e \n", TacsRealPart(S12));
  printf("C    = %15.5e \n", TacsRealPart(C));

  printf("\nStrain Failure Properties \n");
  printf("eXt  = %15.5e \n", TacsRealPart(eXt));
  printf("eXc  = %15.5e \n", TacsRealPart(eXc));
  printf("eYt  = %15.5e \n", TacsRealPart(eYt));
  printf("eYc  = %15.5e \n", TacsRealPart(eYc));
  printf("eS12 = %15.5e \n", TacsRealPart(eS12));

  printf("\nTsai-Wu tensor coefficients\n");
  printf("F1   = %15.5e \n", TacsRealPart(F1));
  printf("F2   = %15.5e \n", TacsRealPart(F2));
  printf("F11  = %15.5e \n", TacsRealPart(F11));
  printf("F12  = %15.5e \n", TacsRealPart(F12));
  printf("F22  = %15.5e \n", TacsRealPart(F22));
  printf("F66  = %15.5e \n", TacsRealPart(F66));
}

// First, the stress transformations
// Transform stress from the global to the local frame  
void OrthoPly::transformStressGlobal2Ply( TacsScalar plyStress[], 
                                          const TacsScalar global[], 
                                          TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;

  plyStress[0] = cos2*global[0] + sin2*global[1] + 2.0*sin1*cos1*global[2];
  plyStress[1] = sin2*global[0] + cos2*global[1] - 2.0*sin1*cos1*global[2];
  plyStress[2] = sin1*cos1*(- global[0] + global[1]) + 
    (cos2 - sin2)*global[2];
}

// Transform stress from the ply frame to the global frame
void OrthoPly::transformStressPly2Global( TacsScalar global[], 
                                          const TacsScalar plyStress[], 
                                          TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;

  global[0] = cos2*plyStress[0] + sin2*plyStress[1] - 2.0*sin1*cos1*plyStress[2];
  global[1] = sin2*plyStress[0] + cos2*plyStress[1] + 2.0*sin1*cos1*plyStress[2];
  global[2] = sin1*cos1*(plyStress[0] - plyStress[1]) + (cos2 - sin2)*plyStress[2];
}

// The sensitivity of the transformation of 
void OrthoPly::transformStressGlobal2PlyAngleSens( TacsScalar plyStress[], 
                                                   const TacsScalar global[], 
                                                   TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = - 2.0*cos1*sin1; // sensitivity of cos**2
  TacsScalar s_sin2 =   2.0*sin1*cos1; // sensitivity of sin**2  
  TacsScalar s_sincos = cos1*cos1 - sin1*sin1; // sensitivity of cos*sin

  plyStress[0] = s_cos2*global[0] + s_sin2*global[1] + 2.0*s_sincos*global[2];
  plyStress[1] = s_sin2*global[0] + s_cos2*global[1] - 2.0*s_sincos*global[2];
  plyStress[2] = s_sincos*(global[1] - global[0]) + (s_cos2 - s_sin2)*global[2];
}

  // Transform stress from the ply frame to the global frame
void OrthoPly::transformStressPly2GlobalAngleSens( TacsScalar global[], 
                                                   const TacsScalar plyStress[], 
                                                   TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = - 2.0*cos1*sin1; // sensitivity of cos**2
  TacsScalar s_sin2 =   2.0*sin1*cos1; // sensitivity of sin**2  
  TacsScalar s_sincos = cos1*cos1 - sin1*sin1; // sensitivity of cos*sin

  global[0] = s_cos2*plyStress[0] + s_sin2*plyStress[1] - 2.0*s_sincos*plyStress[2];
  global[1] = s_sin2*plyStress[0] + s_cos2*plyStress[1] + 2.0*s_sincos*plyStress[2];
  global[2] = s_sincos*(plyStress[0] - plyStress[1]) + 
    (s_cos2 - s_sin2)*plyStress[2];
}


// Next, the strain transformations
// Transform strain from the global to the local frame
void OrthoPly::transformStrainGlobal2Ply( TacsScalar plyStrain[], 
                                          const TacsScalar global[], 
                                          TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);
  
  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;

  plyStrain[0] = cos2*global[0] + sin2*global[1] + sin1*cos1*global[2];
  plyStrain[1] = sin2*global[0] + cos2*global[1] - sin1*cos1*global[2];
  plyStrain[2] = 2.0*sin1*cos1*(global[1] - global[0]) + (cos2 - sin2)*global[2];
}

// Transform stress from the ply frame to the global frame
void OrthoPly::transformStrainPly2Global( TacsScalar global[], 
                                          const TacsScalar plyStrain[], 
                                          TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar cos2 = cos1*cos1;
  TacsScalar sin2 = sin1*sin1;
    
  global[0] = cos2*plyStrain[0] + sin2*plyStrain[1] - sin1*cos1*plyStrain[2];
  global[1] = sin2*plyStrain[0] + cos2*plyStrain[1] + sin1*cos1*plyStrain[2];
  global[2] = 2.0*sin1*cos1*(plyStrain[0] - plyStrain[1]) + 
    (cos2 - sin2)*plyStrain[2];
}

// The sensitivity of the transformation to the 
void OrthoPly::transformStrainGlobal2PlyAngleSens( TacsScalar plyStrain[], 
                                                   const TacsScalar global[], 
                                                   TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = - 2.0*cos1*sin1; // sensitivity of cos**2
  TacsScalar s_sin2 =   2.0*sin1*cos1; // sensitivity of sin**2  
  TacsScalar s_sincos = cos1*cos1 - sin1*sin1; // sensitivity of cos*sin

  plyStrain[0] = s_cos2*global[0] + s_sin2*global[1] + s_sincos*global[2];
  plyStrain[1] = s_sin2*global[0] + s_cos2*global[1] - s_sincos*global[2];
  plyStrain[2] = 2.0*s_sincos*(- global[0] + global[1]) + 
    (s_cos2 - s_sin2)*global[2];
}

// Transform stress from the ply frame to the global frame
void OrthoPly::transformStrainPly2GlobalAngleSens( TacsScalar global[], 
                                                   const TacsScalar plyStrain[],
                                                   TacsScalar angle ){
  TacsScalar cos1 = cos(angle);
  TacsScalar sin1 = sin(angle);

  TacsScalar s_cos2 = - 2.0*cos1*sin1; // sensitivity of cos**2
  TacsScalar s_sin2 =   2.0*sin1*cos1; // sensitivity of sin**2  
  TacsScalar s_sincos = cos1*cos1 - sin1*sin1; // sensitivity of cos*sin

  global[0] = s_cos2*plyStrain[0] + s_sin2*plyStrain[1] - s_sincos*plyStrain[2];
  global[1] = s_sin2*plyStrain[0] + s_cos2*plyStrain[1] + s_sincos*plyStrain[2];
  global[2] = 2.0*s_sincos*(plyStrain[0] - plyStrain[1]) + 
    (s_cos2 - s_sin2)*plyStrain[2];
}

const char*OrthoPly::TACSObjectName(){ return name; }  
