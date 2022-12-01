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

/*
The phase change material (PCM) is defined to be solid below the melting
temperature (Tm) and liquid above it. Phase transition (from solid to liquid or
liquid to solid) requires energy, defined by the material's specific latent heat
(lh). This implementation  approximates this by defining a transition region
between [Tm-dT, Tm+dT], where dT is a small temperature with a default value
of 1.0. In this transition region, the latent heat is added to the material's
specific heat. In the transition region, density and conductivity are linearly
interpolated between the solid and liquid values. The piecewise-linear
transition boundaries are smoothly-approximated using an arctangent function to
make these values differentiable.
*/

#include "TACSPhaseChangeMaterialConstitutive.h"

const char* TACSPhaseChangeMaterialConstitutive::psName =
    "TACSPhaseChangeMaterialConstitutive";

const char* TACSPhaseChangeMaterialConstitutive::getObjectName() {
  return psName;
}

/*
  PhaseChangeMaterialStiffness member function definitions
*/
TACSPhaseChangeMaterialConstitutive::TACSPhaseChangeMaterialConstitutive(
    TACSMaterialProperties* solid_props, TACSMaterialProperties* liquid_props,
    TacsScalar _lh, TacsScalar _Tm, TacsScalar _dT, TacsScalar _t, int _tNum,
    TacsScalar _tlb, TacsScalar _tub) {
  solid_properties = solid_props;
  if (solid_properties) {
    solid_properties->incref();
  }
  liquid_properties = liquid_props;
  if (liquid_properties) {
    liquid_properties->incref();
  }
  lh = _lh;  // material's specific latent heat
  Tm = _Tm;  // melting temperature
  t = _t;
  tNum = _tNum;
  tlb = _tlb;
  tub = _tub;
  dT = _dT;  // small temperature value defining the transition region: [Tm-dT,
             // Tm+dT]
  b = 2.0 * tan(0.99 * M_PI /
                2.0);  // constant used to evaluate transition boundaries
}

TACSPhaseChangeMaterialConstitutive::~TACSPhaseChangeMaterialConstitutive() {
  if (solid_properties) {
    solid_properties->decref();
  }
  if (liquid_properties) {
    liquid_properties->decref();
  }
}

int TACSPhaseChangeMaterialConstitutive::getNumStresses() { return 0; }

// Retrieve the global design variable numbers
int TACSPhaseChangeMaterialConstitutive::getDesignVarNums(int elemIndex,
                                                          int dvLen,
                                                          int dvNums[]) {
  if (tNum >= 0) {
    if (dvNums && dvLen >= 1) {
      dvNums[0] = tNum;
    }
    return 1;
  }
  return 0;
}

// Set the element design variable from the design vector
int TACSPhaseChangeMaterialConstitutive::setDesignVars(int elemIndex, int dvLen,
                                                       const TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSPhaseChangeMaterialConstitutive::getDesignVars(int elemIndex, int dvLen,
                                                       TacsScalar dvs[]) {
  if (tNum >= 0 && dvLen >= 1) {
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSPhaseChangeMaterialConstitutive::getDesignVarRange(int elemIndex,
                                                           int dvLen,
                                                           TacsScalar lb[],
                                                           TacsScalar ub[]) {
  if (tNum >= 0 && dvLen >= 1) {
    if (lb) {
      lb[0] = tlb;
    }
    if (ub) {
      ub[0] = tub;
    }
    return 1;
  }
  return 0;
}

// Compute the phase change coefficient
TacsScalar TACSPhaseChangeMaterialConstitutive::evalTransitionCoef(
    const TacsScalar T) {
  TacsScalar phi = atan(b / dT * (T - Tm)) / M_PI + 0.5;
  return phi;
}

// Compute the phase change coefficient
TacsScalar TACSPhaseChangeMaterialConstitutive::evalTransitionCoefSVSens(
    const TacsScalar T) {
  return (b / dT / (1.0 + pow(b / dT * (T - Tm), 2))) / M_PI;
}

// Evaluate the material's phase (0=solid, 1=liquid)
int TACSPhaseChangeMaterialConstitutive::evalPhase(const TacsScalar T) {
  if (TacsRealPart(T) < TacsRealPart(Tm)) {
    return 0;
  }
  return 1;
}

// Evaluate the material density
TacsScalar TACSPhaseChangeMaterialConstitutive::evalDensity(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar u[]) {
  if (solid_properties && liquid_properties) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar rhol = liquid_properties->getDensity();
    return t * (rhos + (rhol - rhos) * B);
  }
  return 0.0;
}

// Add the derivative of the density
void TACSPhaseChangeMaterialConstitutive::addDensityDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    int dvLen, TacsScalar dfdx[], const TacsScalar u[]) {
  if (solid_properties && liquid_properties && tNum >= 0) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar rhol = liquid_properties->getDensity();
    dfdx[0] += scale * (rhos + (rhol - rhos) * B);
  }
}

void TACSPhaseChangeMaterialConstitutive::addDensitySVSens(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar dfdu[],
    const TacsScalar u[]) {
  if (solid_properties && liquid_properties && tNum >= 0) {
    TacsScalar T = u[0];
    TacsScalar dBdT = evalTransitionCoefSVSens(T);
    TacsScalar rhos = solid_properties->getDensity();
    TacsScalar rhol = liquid_properties->getDensity();
    dfdu[0] += t * (rhol - rhos) * dBdT;
  }
}

TacsScalar TACSPhaseChangeMaterialConstitutive::evalSpecificHeat(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar u[]) {
  if (solid_properties && liquid_properties) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar f1 = -(T - Tm) * (T - Tm) / (dT * dT);
    TacsScalar f2 = pow(dT * dT, 0.5) / 3.14159265359;
    TacsScalar D = exp(f1 / f2);
    return cs + (cl - cs) * B + lh * D;
  }
  return 0.0;
}

void TACSPhaseChangeMaterialConstitutive::addSpecificHeatSVSens(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar dfdu[],
    const TacsScalar u[]) {
  if (solid_properties && liquid_properties && tNum >= 0) {
    TacsScalar T = u[0];
    TacsScalar dBdT = evalTransitionCoefSVSens(T);
    TacsScalar cs = solid_properties->getSpecificHeat();
    TacsScalar cl = liquid_properties->getSpecificHeat();
    TacsScalar f1 = -(T - Tm) * (T - Tm) / (dT * dT);
    TacsScalar f2 = pow(dT * dT, 0.5) / 3.14159265359;
    TacsScalar dDdT = (-2.0 / f2) * exp(f1 / f2) * (T - Tm) / (dT * dT);
    dfdu[0] += (cl - cs) * dBdT + lh * dDdT;
  }
}

void TACSPhaseChangeMaterialConstitutive::evalStress(int elemIndex,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar strain[],
                                                     TacsScalar stress[]) {}

void TACSPhaseChangeMaterialConstitutive::evalTangentStiffness(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar C[]) {}

// Evaluate the heat flux, given the thermal gradient
void TACSPhaseChangeMaterialConstitutive::evalHeatFlux(
    int elemIndex, const double pt[], const TacsScalar X[],
    const TacsScalar grad[], TacsScalar flux[], const TacsScalar u[]) {
  if (solid_properties && liquid_properties) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);

    TacsScalar Kcs[3], Kcl[3], Kc[3];
    solid_properties->evalTangentHeatFlux2D(Kcs);
    liquid_properties->evalTangentHeatFlux2D(Kcl);
    Kc[0] = Kcs[0] + (Kcl[0] - Kcs[0]) * B;
    Kc[1] = Kcs[1] + (Kcl[1] - Kcs[1]) * B;
    Kc[2] = Kcs[2] + (Kcl[2] - Kcs[2]) * B;

    flux[0] = t * (Kc[0] * grad[0] + Kc[1] * grad[1]);
    flux[1] = t * (Kc[1] * grad[0] + Kc[2] * grad[1]);
  }
}

// Evaluate the tangent of the heat flux
void TACSPhaseChangeMaterialConstitutive::evalTangentHeatFlux(
    int elemIndex, const double pt[], const TacsScalar X[], TacsScalar Kc[],
    const TacsScalar u[]) {
  if (solid_properties && liquid_properties) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);

    TacsScalar Kcs[3], Kcl[3];
    solid_properties->evalTangentHeatFlux2D(Kcs);
    liquid_properties->evalTangentHeatFlux2D(Kcl);
    Kc[0] = Kcs[0] + (Kcl[0] - Kcs[0]) * B;
    Kc[1] = Kcs[1] + (Kcl[1] - Kcs[1]) * B;
    Kc[2] = Kcs[2] + (Kcl[2] - Kcs[2]) * B;

    Kc[0] *= t;
    Kc[1] *= t;
    Kc[2] *= t;
  }
}

// Add the derivative of the heat flux
void TACSPhaseChangeMaterialConstitutive::addHeatFluxDVSens(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar grad[], const TacsScalar psi[], int dvLen,
    TacsScalar dfdx[], const TacsScalar u[]) {
  if (solid_properties && liquid_properties && tNum >= 0) {
    TacsScalar T = u[0];
    TacsScalar B = evalTransitionCoef(T);

    TacsScalar Kcs[3], Kcl[3], Kc[3];
    solid_properties->evalTangentHeatFlux2D(Kcs);
    liquid_properties->evalTangentHeatFlux2D(Kcl);
    Kc[0] = Kcs[0] + (Kcl[0] - Kcs[0]) * B;
    Kc[1] = Kcs[1] + (Kcl[1] - Kcs[1]) * B;
    Kc[2] = Kcs[2] + (Kcl[2] - Kcs[2]) * B;

    dfdx[0] += scale * (psi[0] * (Kc[0] * grad[0] + Kc[1] * grad[1]) +
                        psi[1] * (Kc[1] * grad[0] + Kc[2] * grad[1]));
  }
}

void TACSPhaseChangeMaterialConstitutive::addKappaSVSens(int elemIndex,
                                                         const double pt[],
                                                         const TacsScalar X[],
                                                         TacsScalar dfdu[],
                                                         const TacsScalar u[]) {
  if (solid_properties && liquid_properties && tNum >= 0) {
    TacsScalar T = u[0];
    TacsScalar dBdT = evalTransitionCoefSVSens(T);

    TacsScalar Kcs[3], Kcl[3];
    solid_properties->evalTangentHeatFlux2D(Kcs);
    liquid_properties->evalTangentHeatFlux2D(Kcl);
    dfdu[0] = t * (Kcl[0] - Kcs[0]) * dBdT;
    dfdu[1] = t * (Kcl[1] - Kcs[1]) * dBdT;
    dfdu[2] = t * (Kcl[2] - Kcs[2]) * dBdT;
  }
}