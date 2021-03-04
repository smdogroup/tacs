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

#include "TACSShellConstitutive.h"
#include "tacslapack.h"

const char* TACSShellConstitutive::constName = "TACSShellConstitutive";

/*
  Create the shell constitutive
*/
TACSShellConstitutive::TACSShellConstitutive( TACSMaterialProperties *props,
                                              TacsScalar _t,
                                              int _tNum,
                                              TacsScalar _tlb,
                                              TacsScalar _tub ){
  properties = props;
  if (properties){
    properties->incref();
  }

  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
}

TACSShellConstitutive::~TACSShellConstitutive(){
  if (properties){
    properties->decref();
  }
}

int TACSShellConstitutive::getNumStresses(){
  return NUM_STRESSES;
}

// Retrieve the global design variable numbers
int TACSShellConstitutive::getDesignVarNums( int elemIndex,
                                             int dvLen, int dvNums[] ){
  if (tNum >= 0){
    if (dvNums && dvLen >= 1){
      dvNums[0] = tNum;
    }
    return 1;
  }
  return 0;
}

// Set the element design variable from the design vector
int TACSShellConstitutive::setDesignVars( int elemIndex,
                                          int dvLen,
                                          const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSShellConstitutive::getDesignVars( int elemIndex,
                                          int dvLen,
                                          TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSShellConstitutive::getDesignVarRange( int elemIndex,
                                              int dvLen,
                                              TacsScalar lb[],
                                              TacsScalar ub[] ){
  if (tNum >= 0 && dvLen >= 1){
    if (lb){ lb[0] = tlb; }
    if (ub){ ub[0] = tub; }
    return 1;
  }
  return 0;
}

// Evaluate the material density
TacsScalar TACSShellConstitutive::evalDensity( int elemIndex,
                                               const double pt[],
                                               const TacsScalar X[] ){
  if (properties){
    return t*properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSShellConstitutive::addDensityDVSens( int elemIndex,
                                              TacsScalar scale,
                                              const double pt[],
                                              const TacsScalar X[],
                                              int dvLen,
                                              TacsScalar dfdx[] ){
  if (properties && tNum >= 0){
    dfdx[0] += scale*properties->getDensity();
  }
}

// Evaluate the specific heat
TacsScalar TACSShellConstitutive::evalSpecificHeat( int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[] ){
  if (properties){
    return properties->getSpecificHeat();
  }
  return 0.0;
}

// Evaluate the stresss
void TACSShellConstitutive::evalStress( int elemIndex,
                                        const double pt[],
                                        const TacsScalar X[],
                                        const TacsScalar e[],
                                        TacsScalar s[] ){
  if (properties){
    TacsScalar A[6], B[6], D[6], As[3], drill;

    // Compute the tangent stiffness matrix
    properties->evalTangentStiffness2D(A);

    // The bending-stretch coupling matrix is zero in this case
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

    // Scale the in-plane matrix and bending stiffness
    // matrix by the appropriate quantities
    TacsScalar I = t*t*t/12.0;
    for ( int i = 0; i < 6; i++ ){
      D[i] = I*A[i];
      A[i] *= t;
    }

    // Set the through-thickness shear stiffness
    As[0] = As[2] = (5.0/6.0)*A[5];
    As[1] = 0.0;

    drill = 0.5*DRILLING_REGULARIZATION*(As[0] + As[2]);

    // Evaluate the stress
    evalStress(A, B, D, As, drill, e, s);
  }
  else {
    s[0] = s[1] = s[2] = 0.0;
    s[3] = s[4] = s[5] = 0.0;
    s[6] = s[7] = s[8] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSShellConstitutive::evalTangentStiffness( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  TacsScalar C[] ){
  if (properties){
    TacsScalar *A = &C[0];
    TacsScalar *B = &C[6];
    TacsScalar *D = &C[12];
    TacsScalar *As = &C[18];

    // Compute the tangent stiffness matrix
    properties->evalTangentStiffness2D(A);

    // The bending-stretch coupling matrix is zero in this case
    B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

    // Scale the in-plane matrix and bending stiffness
    // matrix by the appropriate quantities
    TacsScalar I = t*t*t/12.0;
    for ( int i = 0; i < 6; i++ ){
      D[i] = I*A[i];
      A[i] *= t;
    }

    // Set the through-thickness shear stiffness
    As[0] = As[2] = (5.0/6.0)*A[5];
    As[1] = 0.0;

    C[21] = 0.5*DRILLING_REGULARIZATION*(As[0] + As[2]);
  }
  else {
    memset(C, 0, 22*sizeof(TacsScalar));
  }
}

// Extract the tangent stiffness components from the matrix
void TACSShellConstitutive::extractTangenttStiffness( const TacsScalar *C,
                                                      const TacsScalar **A,
                                                      const TacsScalar **B,
                                                      const TacsScalar **D,
                                                      const TacsScalar **As,
                                                      TacsScalar *drill ){
  if (A){ *A = &C[0]; }
  if (B){ *B = &C[6]; }
  if (D){ *D = &C[12]; }
  if (As){ *As = &C[18]; }
  if (drill){ *drill = C[21]; }
}

// Add the contribution
void TACSShellConstitutive::addStressDVSens( int elemIndex,
                                             TacsScalar scale,
                                             const double pt[],
                                             const TacsScalar X[],
                                             const TacsScalar strain[],
                                             const TacsScalar psi[],
                                             int dvLen, TacsScalar dfdx[] ){}

// Evaluate the thermal strain
void TACSShellConstitutive::evalThermalStrain( int elemIndex,
                                               const double pt[],
                                               const TacsScalar X[],
                                               TacsScalar theta,
                                               TacsScalar e[] ){
  if (properties){
    properties->evalThermalStrain2D(e);
    e[0] *= theta;
    e[1] *= theta;
    e[2] *= theta;

    e[3] = e[4] = e[5] = 0.0;
    e[6] = e[7] = e[8] = 0.0;
  }
  else {
    e[0] = e[1] = e[2] = 0.0;
    e[3] = e[4] = e[5] = 0.0;
    e[6] = e[7] = e[8] = 0.0;
  }
}

// Evaluate the heat flux, given the thermal gradient
void TACSShellConstitutive::evalHeatFlux( int elemIndex,
                                          const double pt[],
                                          const TacsScalar X[],
                                          const TacsScalar grad[],
                                          TacsScalar flux[] ){
  if (properties){
    TacsScalar Kc[3];
    properties->evalTangentHeatFlux2D(Kc);
    flux[0] = t*(Kc[0]*grad[0] + Kc[1]*grad[1]);
    flux[1] = t*(Kc[1]*grad[0] + Kc[2]*grad[1]);
  }
}

// Evaluate the tangent of the heat flux
void TACSShellConstitutive::evalTangentHeatFlux( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 TacsScalar Kc[] ){
  if (properties){
    properties->evalTangentHeatFlux2D(Kc);
    Kc[0] *= t;  Kc[1] *= t;  Kc[2] *= t;
  }
}

/*
  Return the constitutive name
*/
const char* TACSShellConstitutive::getObjectName(){
  return constName;
}

/*
  Set the default drilling regularization value
*/
double TACSShellConstitutive::DRILLING_REGULARIZATION = 10.0;

/*
  Set the drilling stiffness regularization parameter
*/
void TACSShellConstitutive::setDrillingRegularization( double kval ){
  DRILLING_REGULARIZATION = kval;
}
