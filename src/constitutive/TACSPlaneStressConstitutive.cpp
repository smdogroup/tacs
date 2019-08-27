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
TACSPlaneStressConstitutive::TACSPlaneStressConstitutive( TACSMaterialProperties *props ){
  properties = props;
  if (properties){
    properties->incref();
  }
}

TACSPlaneStressConstitutive::~TACSPlaneStressConstitutive(){
  if (properties){
    properties->decref();
  }
}

int TACSPlaneStressConstitutive::getNumStresses(){
  return NUM_STRESSES;
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

    s[0] = C[0]*e[0] + C[1]*e[1] + C[2]*e[2];
    s[1] = C[1]*e[0] + C[3]*e[1] + C[4]*e[2];
    s[2] = C[2]*e[0] + C[4]*e[1] + C[5]*e[2];
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
  }
  else {
    C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
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
    return properties->getDensity();
  }

  return 0.0;
}

// Evaluate the material failure index
TacsScalar TACSPlaneStressConstitutive::failure( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar e[] ){
  if (properties){
    TacsScalar s[3];
    evalStress(elemIndex, pt, X, e, s);
    return properties->vonMisesFailure2D(s);
  }
  return 0.0;
}
