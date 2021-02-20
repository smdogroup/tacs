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

#include "TACSCompositeShellConstitutive.h"
#include "tacslapack.h"

const char* TACSCompositeShellConstitutive::constName = "TACSCompositeShellConstitutive";

/*
  Create the shell constitutive
*/
TACSCompositeShellConstitutive::TACSCompositeShellConstitutive( TACSOrthotropicPly **_ply_props,
                                                                const TacsScalar *_ply_thickness,
                                                                const TacsScalar *_ply_angles,
                                                                int _num_plies ){
  num_plies = _num_plies;
  ply_thickness = new TacsScalar[ num_plies ];
  ply_angles = new TacsScalar[ num_plies ];
  ply_props = new TACSOrthotropicPly[ num_plies ];

  for ( int i = 0; i < num_plies; i++ ){
    ply_props[i] = _ply_props[i];
    ply_props[i]->incref();

    ply_thickness[i] = _ply_thickness[i];
    ply_angles[i] = _ply_angles[i];
  }
}

TACSCompositeShellConstitutive::~TACSCompositeShellConstitutive(){
  for ( int i = 0; i < num_plies; i++ ){
    ply_props[i]->decref();
  }
  delete [] ply_props;
  delete [] ply_thickness;
  delete [] ply_angles;
}


// Retrieve the global design variable numbers
int TACSCompositeShellConstitutive::getDesignVarNums( int elemIndex,
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
int TACSCompositeShellConstitutive::setDesignVars( int elemIndex,
                                                   int dvLen,
                                                   const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSCompositeShellConstitutive::getDesignVars( int elemIndex,
                                                   int dvLen,
                                                   TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSCompositeShellConstitutive::getDesignVarRange( int elemIndex,
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
TacsScalar TACSCompositeShellConstitutive::evalDensity( int elemIndex,
                                                        const double pt[],
                                                        const TacsScalar X[] ){
  if (properties){
    return t*properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSCompositeShellConstitutive::addDensityDVSens( int elemIndex,
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
TacsScalar TACSCompositeShellConstitutive::evalSpecificHeat( int elemIndex,
                                                             const double pt[],
                                                             const TacsScalar X[] ){
  if (properties){
    return properties->getSpecificHeat();
  }
  return 0.0;
}

// Evaluate the stresss
void TACSCompositeShellConstitutive::evalStress( int elemIndex,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar e[],
                                                 TacsScalar s[] ){
  TacsScalar A[6], B[6], D[6], As[3], drill;

  // Zero the stiffness matrices
  for ( int k = 0; k < 6; k++ ){
    A[k] = B[k] = D[k] = 0.0;
  }

  for ( int k = 0; k < 3; k++ ){
    As[k] = 0.0;
  }

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for ( int i = 0; i < num_plies; i++ ){
    t += ply_thickness[i];
  }

  // Compute the contribution to the stiffness from each layer
  TacsScalar t0 = -0.5*t;
  for ( int k = 0; k < num_plies; k++ ){
    TacsScalar Qbar[6], Abar[3];
    ply_props[k]->calculateQbar(Qbar, ply_angles[k]);
    ply_props[k]->calculateAbar(Abar, ply_angles[k]);

    TacsScalar t1 = t0 + ply_thickness[k];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5*(t1*t1 - t0*t0);
    TacsScalar d = 1.0/3.0*(t1*t1*t1 - t0*t0*t0);

    for ( int i = 0; i < 6; i++ ){
      A[i] += a*Qbar[i];
      B[i] += b*Qbar[i];
      D[i] += d*Qbar[i];
    }

    for ( int i = 0; i < 3; i++ ){
      As[i] += kcorr*a*Abar[i];
    }

    // Update the position of the bottom interface
    t0 = t1;
  }

  drill = 0.5*DRILLING_REGULARIZATION*(As[0] + As[2]);

  // Evaluate the stress
  evalStress(A, B, D, As, drill, e, s);
}

// Evaluate the tangent stiffness
void TACSCompositeShellConstitutive::evalTangentStiffness( int elemIndex,
                                                           const double pt[],
                                                           const TacsScalar X[],
                                                           TacsScalar C[] ){
  TacsScalar *A = &C[0];
  TacsScalar *B = &C[6];
  TacsScalar *D = &C[12];
  TacsScalar *As = &C[18];

  // Zero the stiffness matrices
  for ( int k = 0; k < 6; k++ ){
    A[k] = B[k] = D[k] = 0.0;
  }

  for ( int k = 0; k < 3; k++ ){
    As[k] = 0.0;
  }

  // Compute the total thickness of the laminate
  TacsScalar t = 0.0;
  for ( int i = 0; i < num_plies; i++ ){
    t += ply_thickness[i];
  }

  // Compute the contribution to the stiffness from each layer
  TacsScalar t0 = -0.5*t;
  for ( int k = 0; k < num_plies; k++ ){
    TacsScalar Qbar[6], Abar[3];
    ply_props[k]->calculateQbar(Qbar, ply_angles[k]);
    ply_props[k]->calculateAbar(Abar, ply_angles[k]);

    TacsScalar t1 = t0 + ply_thickness[k];

    TacsScalar a = (t1 - t0);
    TacsScalar b = 0.5*(t1*t1 - t0*t0);
    TacsScalar d = 1.0/3.0*(t1*t1*t1 - t0*t0*t0);

    for ( int i = 0; i < 6; i++ ){
      A[i] += a*Qbar[i];
      B[i] += b*Qbar[i];
      D[i] += d*Qbar[i];
    }

    for ( int i = 0; i < 3; i++ ){
      As[i] += kcorr*a*Abar[i];
    }

    // Update the position of the bottom interface
    t0 = t1;
  }

  C[21] = 0.5*DRILLING_REGULARIZATION*(As[0] + As[2]);
}

// Add the contribution
void TACSCompositeShellConstitutive::addStressDVSens( int elemIndex,
                                                      TacsScalar scale,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      const TacsScalar strain[],
                                                      const TacsScalar psi[],
                                                      int dvLen, TacsScalar dfdx[] ){}

// Evaluate the thermal strain
void TACSCompositeShellConstitutive::evalThermalStrain( int elemIndex,
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
void TACSCompositeShellConstitutive::evalHeatFlux( int elemIndex,
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
void TACSCompositeShellConstitutive::evalTangentHeatFlux( int elemIndex,
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
const char* TACSCompositeShellConstitutive::getObjectName(){
  return constName;
}
