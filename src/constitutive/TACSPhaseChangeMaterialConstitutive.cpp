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

#include "TACSPhaseChangeMaterialConstitutive.h"

const char* TACSPhaseChangeMaterialConstitutive::psName = "TACSPhaseChangeMaterialConstitutive";

const char* TACSPhaseChangeMaterialConstitutive::getObjectName(){
  return psName;
}

/*
  PhaseChangeMaterialStiffness member function definitions
*/
TACSPhaseChangeMaterialConstitutive::TACSPhaseChangeMaterialConstitutive( TACSMaterialProperties *solid_props,
                                                                          TACSMaterialProperties *liquid_props,
                                                                          TacsScalar _lh,
                                                                          TacsScalar _mt,
                                                                          TacsScalar _t,
                                                                          int _tNum,
                                                                          TacsScalar _tlb,
                                                                          TacsScalar _tub ){
  solid_properties = solid_props;
  if (solid_properties){
    solid_properties->incref();
  }
  liquid_properties = liquid_props;
  if (liquid_properties){
    liquid_properties->incref();
  }
  lh = _lh;
  mt = _mt;
  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
}

TACSPhaseChangeMaterialConstitutive::~TACSPhaseChangeMaterialConstitutive(){
  if (solid_properties){
    solid_properties->decref();
  }
  if (liquid_properties){
    liquid_properties->decref();
  }
}

int TACSPhaseChangeMaterialConstitutive::getNumStresses(){
  return 0;
}

// Retrieve the global design variable numbers
int TACSPhaseChangeMaterialConstitutive::getDesignVarNums( int elemIndex,
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
int TACSPhaseChangeMaterialConstitutive::setDesignVars( int elemIndex,
                                                        int dvLen,
                                                        const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSPhaseChangeMaterialConstitutive::getDesignVars( int elemIndex,
                                                        int dvLen,
                                                        TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSPhaseChangeMaterialConstitutive::getDesignVarRange( int elemIndex,
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

// Evaluate the temperature at a given element
TacsScalar TACSPhaseChangeMaterialConstitutive::evalTemperature( int elemIndex,
                                                                 const double pt[],
                                                                 const TacsScalar X[],
                                                                 const TacsScalar U ){
  if (solid_properties && liquid_properties){
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar rhol = liquid_properties->getDensity();
    TacsScalar Ut = cs*mt;
    TacsScalar Um = (cs*mt + lh);

    if (U < Ut){
      return U/(cs);
    }
    else if ((U >= Ut) && (U < Um)){
      return mt;
    }
      return (U-Um)/(cl) + mt;
  }
  return 0.0;
}

// Check if the element is undergoing phase change
int TACSPhaseChangeMaterialConstitutive::checkPhaseChange( const TacsScalar U ){
  if (solid_properties){
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar Ut = cs*mt;
    TacsScalar Um = (cs*mt + lh);
    if ((U >= Ut) && (U < Um)){
      return 1;
    }
    return 0;
  }
  return 0;
}

// Evaluate the material's phase (0=solid, 1=liquid)
int TACSPhaseChangeMaterialConstitutive::evalPhase( const TacsScalar U ){
  if (solid_properties){
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar Um = (cs*mt + lh);
    if (U <  Um){
      return 0;
    }
    return 1;
  }
  return 0;
}

// Evaluate the material density
TacsScalar TACSPhaseChangeMaterialConstitutive::evalDensity( int elemIndex,
                                                             const double pt[],
                                                             const TacsScalar X[] ){
  if (solid_properties && liquid_properties){
    return t*solid_properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSPhaseChangeMaterialConstitutive::addDensityDVSens( int elemIndex,
                                                            TacsScalar scale,
                                                            const double pt[],
                                                            const TacsScalar X[],
                                                            int dvLen,
                                                            TacsScalar dfdx[] ){
  if (solid_properties && liquid_properties && tNum >= 0){
      dfdx[0] += scale*solid_properties->getDensity();
  }
}

TacsScalar TACSPhaseChangeMaterialConstitutive::evalSpecificHeat( int elemIndex,
                                                                  const double pt[],
                                                                  const TacsScalar X[] ){
  return 0.0;
}

void TACSPhaseChangeMaterialConstitutive::evalStress( int elemIndex,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      const TacsScalar strain[],
                                                      TacsScalar stress[] ){}

void TACSPhaseChangeMaterialConstitutive::evalTangentStiffness( int elemIndex,
                                                                const double pt[],
                                                                const TacsScalar X[],
                                                                TacsScalar C[] ){}

// Evaluate the heat flux, given the thermal gradient
void TACSPhaseChangeMaterialConstitutive::evalHeatFlux( int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        const TacsScalar grad[],
                                                        TacsScalar flux[],
                                                        const TacsScalar U ){
  if (solid_properties && liquid_properties){
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar rhol = liquid_properties->getDensity();
    TacsScalar Ut = cs*mt;
    TacsScalar Um = (cs*mt + lh);
    TacsScalar t_sc;

    TacsScalar Kc[3];
    if (evalPhase(U)){
      liquid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = rhos/rhol;
    }
    else{
      solid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = 1.0;
    }

    if (U < Ut){
      Kc[0] *= 1.0/(cs);
      Kc[1] *= 1.0/(cs);
      Kc[2] *= 1.0/(cs);
    }
    else if ((U >= Ut) && (U <= Um)){
      Kc[0] *= 0.0;
      Kc[1] *= 0.0;
      Kc[2] *= 0.0;
    }
    else{
      Kc[0] *= 1.0/(cl);
      Kc[1] *= 1.0/(cl);
      Kc[2] *= 1.0/(cl);
    }

    flux[0] = t*t_sc*(Kc[0]*grad[0]+Kc[1]*grad[1]);
    flux[1] = t*t_sc*(Kc[1]*grad[0]+Kc[2]*grad[1]);
  }
}

// Evaluate the tangent of the heat flux
void TACSPhaseChangeMaterialConstitutive::evalTangentHeatFlux( int elemIndex,
                                                               const double pt[],
                                                               const TacsScalar X[],
                                                               TacsScalar Kc[],
                                                               const TacsScalar U ){
  if (solid_properties && liquid_properties){
    TacsScalar t_sc;
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar rhol = liquid_properties->getDensity();
    TacsScalar Ut = rhos*cs*mt;
    TacsScalar Um = rhos*(cs*mt + lh);
    if (evalPhase(U)){
      liquid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = rhos/rhol;
    }
    else{
      solid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = 1.0;
    }

    if (U < Ut){
      Kc[0] *= 1.0/(cs);
      Kc[1] *= 1.0/(cs);
      Kc[2] *= 1.0/(cs);
    }
    else if ((U >= Ut) && (U <= Um)){
      Kc[0] *= 0.0;
      Kc[1] *= 0.0;
      Kc[2] *= 0.0;
    }
    else{
      Kc[0] *= 1.0/(cl);
      Kc[1] *= 1.0/(cl);
      Kc[2] *= 1.0/(cl);
    }

    Kc[0] *= t*t_sc;  Kc[1] *= t*t_sc;  Kc[2] *= t*t_sc;
  }
}

// Add the derivative of the heat flux
void TACSPhaseChangeMaterialConstitutive::addHeatFluxDVSens( int elemIndex,
                                                             TacsScalar scale,
                                                             const double pt[],
                                                             const TacsScalar X[],
                                                             const TacsScalar grad[],
                                                             const TacsScalar psi[],
                                                             int dvLen,
                                                             TacsScalar dfdx[],
                                                             const TacsScalar U ){
  if (solid_properties && liquid_properties && tNum >= 0){
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar rhol = liquid_properties->getDensity();
    TacsScalar Ut = cs*mt;
    TacsScalar Um = (cs*mt + lh);
    TacsScalar t_sc;

    TacsScalar Kc[3];
    if (evalPhase(U)){
      liquid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = rhos/rhol;
    }
    else{
      solid_properties->evalTangentHeatFlux2D(Kc);
      t_sc = 1.0;
    }

    if (U < Ut){
      Kc[0] *= 1.0/(cs);
    }
    else if ((U >= Ut) && (U <= Um)){
      Kc[0] *= 0.0;
    }
    else{
      Kc[0] *= 1.0/(cl);
    }

    dfdx[0] += scale*Kc[0]*(psi[0]*grad[0] + psi[1]*grad[1]);
  }
}