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

#include "SolidStiffness.h"
#include "YSlibrary.h"

const char * SolidStiffness::constName = "SolidStiffness";

const char * SolidStiffness::constitutiveName(){
  return constName;
}
  
/*
  SolidStiffness member function definitions: Definitions for a simple
  beam element
*/
SolidStiffness::SolidStiffness( TacsScalar _rho, TacsScalar _E,
                                TacsScalar _nu, TacsScalar _ys, int _eNum ){
  rho = _rho;
  E = _E;
  nu = _nu;
  ys = _ys;
  
  TacsScalar D = E/((1.0 + nu)*(1.0 - 2.0*nu));
  C[0] = C[3] = C[5] = (1.0 - nu)*D;
  C[1] = C[2] = C[4] = nu*D;

  G23 = G13 = G12 = 0.5*E/(1.0 + nu);

  eNum = _eNum;
}

/*
  Default constructor that just assigns zero properties
*/
SolidStiffness::SolidStiffness(){
  rho = 0.0;
  E = 0.0;
  nu = 0.0;
  C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
  G23 = G13 = G12 = 0.0;
  eNum = -1;
}

int SolidStiffness::getNumStresses(){ return NUM_STRESSES; }

/*
  Set the design variable values from the vector
*/
void SolidStiffness::setDesignVars( const TacsScalar dvs[], int dvLen ){
  if (eNum >= 0 && eNum < dvLen){
    E = dvs[eNum];

    TacsScalar D = E/((1.0 + nu)*(1.0 - 2.0*nu));
    C[0] = C[3] = C[5] = (1.0 - nu)*D;
    C[1] = C[2] = C[4] = nu*D;

    G23 = G13 = G12 = 0.5*E/(1.0 + nu);
  }
}

/*
  Get the design variable values from the object
*/
void SolidStiffness::getDesignVars( TacsScalar dvs[], int dvLen ){
  if (eNum >= 0 && eNum < dvLen){
    dvs[eNum] = E;
  }
}

/*
  Given the strain, compute the stress at the given parametric point
*/
void SolidStiffness::calculateStress( const double gpt[], 
                                      const TacsScalar strain[],
                                      TacsScalar stress[] ){
  calcStress(strain, stress);
}

/*
  Add the derivative of the product of the stress with an arbitrary
  input vector psi, given the strain
*/
void SolidStiffness::addStressDVSens( const double pt[],
                                      const TacsScalar e[],
                                      TacsScalar alpha,
                                      const TacsScalar psi[],
                                      TacsScalar dvSens[], int dvLen ){

  if (eNum >= 0 && eNum < dvLen){
    TacsScalar C0 = alpha*(1.0 - nu)/((1.0 + nu)*(1.0 - 2.0*nu));
    TacsScalar C1 = alpha*nu/((1.0 + nu)*(1.0 - 2.0*nu));
    TacsScalar G  = alpha*0.5/(1.0 + nu);

    dvSens[eNum] += psi[0]*(C0*e[0] + C1*e[1] + C1*e[2]);
    dvSens[eNum] += psi[1]*(C1*e[0] + C0*e[1] + C1*e[2]);
    dvSens[eNum] += psi[2]*(C1*e[0] + C1*e[1] + C0*e[2]);

    dvSens[eNum] += psi[3]*G*e[3];
    dvSens[eNum] += psi[4]*G*e[4];
    dvSens[eNum] += psi[5]*G*e[5];
  }
}

/*
  Compute the failure function
*/
void SolidStiffness::failure( const double pt[], 
                              const TacsScalar e[], 
                              TacsScalar *fail ){
  TacsScalar s[6];
  calcStress(e, s);
  *fail = VonMisesFailure3D(s, ys);
}

/*
  Evaluate the derivative of the failure function with respect to the
  strain
*/
void SolidStiffness::failureStrainSens( const double pt[], 
                                        const TacsScalar e[],
                                        TacsScalar sens[] ){
  TacsScalar s[6], ssens[6];
  calcStress(e, s);
  VonMisesFailure3DStressSens(ssens, s, ys);
  calcStress(ssens, sens);  
}

/*
  Compute the derivative with respect to the design variables
*/
void SolidStiffness::addFailureDVSens( const double pt[],
                                       const TacsScalar e[],
                                       TacsScalar alpha,
                                       TacsScalar dvSens[], int dvLen ){
  if (eNum >= 0 && eNum < dvLen ){
    TacsScalar s[6], ssens[6];
    calcStress(e, s);
    VonMisesFailure3DStressSens(ssens, s, ys);

    TacsScalar C0 = alpha*(1.0 - nu)/((1.0 + nu)*(1.0 - 2.0*nu));
    TacsScalar C1 = alpha*nu/((1.0 + nu)*(1.0 - 2.0*nu));
    TacsScalar G  = alpha*0.5/(1.0 + nu);

    dvSens[eNum] += ssens[0]*(C0*e[0] + C1*e[1] + C1*e[2]);
    dvSens[eNum] += ssens[1]*(C1*e[0] + C0*e[1] + C1*e[2]);
    dvSens[eNum] += ssens[2]*(C1*e[0] + C1*e[1] + C0*e[2]);

    dvSens[eNum] += ssens[3]*G*e[3];
    dvSens[eNum] += ssens[4]*G*e[4];
    dvSens[eNum] += ssens[5]*G*e[5];
  }
}
  
