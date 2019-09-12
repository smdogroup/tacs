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

#include "TACSPlaneStressConstitutive.h"

const char* TACSPlaneStressConstitutive::psName = "TACSPlaneStressConstitutive";

const char* TACSPlaneStressConstitutive::getObjectName(){
  return psName;
}

/*
  PlaneStressStiffness member function definitions
*/
TACSPlaneStressConstitutive::TACSPlaneStressConstitutive( TACSMaterialProperties *props,
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

TACSPlaneStressConstitutive::~TACSPlaneStressConstitutive(){
  if (properties){
    properties->decref();
  }
}

int TACSPlaneStressConstitutive::getNumStresses(){
  return NUM_STRESSES;
}

// Retrieve the global design variable numbers
int TACSPlaneStressConstitutive::getDesignVarNums( int elemIndex,
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
void TACSPlaneStressConstitutive::setDesignVars( int elemIndex,
                                                 int dvLen,
                                                 const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
  }
}

// Get the element design variables values
void TACSPlaneStressConstitutive::getDesignVars( int elemIndex,
                                                 int dvLen,
                                                 TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
  }
}

// Get the lower and upper bounds for the design variable values
void TACSPlaneStressConstitutive::getDesignVarRange( int elemIndex,
                                                     int dvLen,
                                                     TacsScalar lb[],
                                                     TacsScalar ub[] ){
  if (tNum >= 0 && dvLen >= 1){
    if (lb){ lb[0] = tlb; }
    if (ub){ ub[0] = tub; }
  }
}

/**
  Evaluate the stresss
*/
void TACSPlaneStressConstitutive::evalStress( int elemIndex,
                                              const double pt[],
                                              const TacsScalar X[],
                                              const TacsScalar e[],
                                              TacsScalar s[] ){
  TacsScalar C[6];
  if (properties){
    properties->evalTangentStiffness2D(C);

    s[0] = t*(C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = t*(C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = t*(C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);
  }
  else {
    s[0] = s[1] = s[2] = 0.0;
  }
}

/**
  Evaluate the tangent stiffness
*/
void TACSPlaneStressConstitutive::evalTangentStiffness( int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        TacsScalar C[] ){
  if (properties){
    properties->evalTangentStiffness2D(C);
    C[0] *= t;  C[1] *= t;  C[2] *= t;
    C[3] *= t;  C[4] *= t;  C[5] *= t;
  }
  else {
    C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
  }
}

void TACSPlaneStressConstitutive::addStressDVSens( int elemIndex,
                                                   const double pt[],
                                                   const TacsScalar X[],
                                                   const TacsScalar e[],
                                                   TacsScalar scale,
                                                   const TacsScalar psi[],
                                                   int dvLen,
                                                   TacsScalar dvSens[] ){
  if (properties && tNum >= 0){
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    // Compute the derivative w.r.t. the design vector
    dvSens[0] += scale*(s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2]);
  }
}

/**
  Evaluate the thermal strain
*/
void TACSPlaneStressConstitutive::evalThermalStrain( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     TacsScalar e[] ){
  if (properties){
    properties->evalThermalStrain2D(e);
    e[0] *= t;  e[1] *= t;  e[2] *= t;
  }
  else {
    e[0] = e[1] = e[2] = 0.0;
  }
}

// Evaluate the material density
TacsScalar TACSPlaneStressConstitutive::evalDensity( int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[] ){
  if (properties){
    return t*properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSPlaneStressConstitutive::addDensityDVSens( int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar scale,
                                                    int dvLen,
                                                    TacsScalar dvSens[] ){
  if (properties && tNum >= 0){
    dvSens[0] += scale*properties->getDensity();
  }
}

// Evaluate the material failure index
TacsScalar TACSPlaneStressConstitutive::failure( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar e[] ){
  if (properties){
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = t*(C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = t*(C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = t*(C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    return properties->vonMisesFailure2D(s);
  }
  return 0.0;
}


// Evaluate the derivative of the failure criteria w.r.t. strain
TacsScalar TACSPlaneStressConstitutive::failureStrainSens( int elemIndex,
                                                           const double pt[],
                                                           const TacsScalar X[],
                                                           const TacsScalar e[],
                                                           TacsScalar dfde[] ){
  if (properties){
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3];
    s[0] = t*(C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s[1] = t*(C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s[2] = t*(C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    TacsScalar sens[3];
    TacsScalar fail = properties->vonMisesFailure2DStressSens(s, sens);

    dfde[0] = t*(C[0]*sens[0] + C[1]*sens[1] + C[2]*sens[2]);
    dfde[1] = t*(C[1]*sens[0] + C[3]*sens[1] + C[4]*sens[2]);
    dfde[2] = t*(C[2]*sens[0] + C[4]*sens[1] + C[5]*sens[2]);

    return fail;
  }
  return 0.0;
}

// Add the derivative of the failure w.r.t. design variables
void TACSPlaneStressConstitutive::addFailureDVSens( int elemIndex,
                                                    const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar e[],
                                                    TacsScalar scale,
                                                    int dvLen, TacsScalar dvSens[] ){

  if (properties && tNum >= 0){
    TacsScalar C[6];
    properties->evalTangentStiffness2D(C);

    TacsScalar s[3], s0[3];
    s0[0] = (C[0]*e[0] + C[1]*e[1] + C[2]*e[2]);
    s0[1] = (C[1]*e[0] + C[3]*e[1] + C[4]*e[2]);
    s0[2] = (C[2]*e[0] + C[4]*e[1] + C[5]*e[2]);

    s[0] = t*s0[0];
    s[1] = t*s0[1];
    s[2] = t*s0[2];

    TacsScalar sens[3];
    properties->vonMisesFailure2DStressSens(s, sens);

    // Compute the derivative w.r.t. the design vector
    dvSens[0] += scale*(sens[0]*s0[0] + sens[1]*s0[1] + sens[2]*s0[2]);
  }
}
