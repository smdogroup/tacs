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

#include "TACSIsoShellConstitutive.h"
#include "TACSElementAlgebra.h"

const char* TACSIsoShellConstitutive::constName = "TACSIsoShellConstitutive";

/*
  Create the shell constitutive object
*/
TACSIsoShellConstitutive::TACSIsoShellConstitutive( TACSMaterialProperties *props,
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
  kcorr = 5.0/6.0;
}

TACSIsoShellConstitutive::~TACSIsoShellConstitutive(){
  if (properties){
    properties->decref();
  }
}

// Retrieve the global design variable numbers
int TACSIsoShellConstitutive::getDesignVarNums( int elemIndex,
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
int TACSIsoShellConstitutive::setDesignVars( int elemIndex,
                                             int dvLen,
                                             const TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    t = dvs[0];
    return 1;
  }
  return 0;
}

// Get the element design variables values
int TACSIsoShellConstitutive::getDesignVars( int elemIndex,
                                             int dvLen,
                                             TacsScalar dvs[] ){
  if (tNum >= 0 && dvLen >= 1){
    dvs[0] = t;
    return 1;
  }
  return 0;
}

// Get the lower and upper bounds for the design variable values
int TACSIsoShellConstitutive::getDesignVarRange( int elemIndex,
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
TacsScalar TACSIsoShellConstitutive::evalDensity( int elemIndex,
                                                  const double pt[],
                                                  const TacsScalar X[] ){
  if (properties){
    return t*properties->getDensity();
  }
  return 0.0;
}

// Add the derivative of the density
void TACSIsoShellConstitutive::addDensityDVSens( int elemIndex,
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
TacsScalar TACSIsoShellConstitutive::evalSpecificHeat( int elemIndex,
                                                       const double pt[],
                                                       const TacsScalar X[] ){
  if (properties){
    return properties->getSpecificHeat();
  }
  return 0.0;
}

// Evaluate the stresss
void TACSIsoShellConstitutive::evalStress( int elemIndex,
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
    computeStress(A, B, D, As, drill, e, s);
  }
  else {
    s[0] = s[1] = s[2] = 0.0;
    s[3] = s[4] = s[5] = 0.0;
    s[6] = s[7] = s[8] = 0.0;
  }
}

// Evaluate the tangent stiffness
void TACSIsoShellConstitutive::evalTangentStiffness( int elemIndex,
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

// Add the contribution
void TACSIsoShellConstitutive::addStressDVSens( int elemIndex,
                                                TacsScalar scale,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar e[],
                                                const TacsScalar psi[],
                                                int dvLen, TacsScalar dfdx[] ){
  if (properties && tNum >= 0){
    // Compute the tangent stiffness matrix
    TacsScalar A[6];
    properties->evalTangentStiffness2D(A);

    TacsScalar dI = t*t/4.0;

    dfdx[0] += scale*(
      mat3x3SymmInner(A, &psi[0], &e[0]) +
      dI*mat3x3SymmInner(A, &psi[3], &e[3]) +
      (5.0/6.0)*A[5]*(psi[6]*e[6] + psi[7]*e[7] +
        DRILLING_REGULARIZATION*psi[8]*e[8]));
  }
}

// Evaluate the thermal strain
void TACSIsoShellConstitutive::evalThermalStrain( int elemIndex,
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
void TACSIsoShellConstitutive::evalHeatFlux( int elemIndex,
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
void TACSIsoShellConstitutive::evalTangentHeatFlux( int elemIndex,
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
const char* TACSIsoShellConstitutive::getObjectName(){
  return constName;
}

/*
  Add the derivative of the product of the stress to the design
  derivative array
*/
/*
  void isoFSDTStiffness::addStiffnessDVSens( const double pt[],
  const TacsScalar e[],
  const TacsScalar psi[],
  TacsScalar rotPsi,
  TacsScalar fdvSens[], int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
  // Compute the derivative of the stiffness coefficients
  TacsScalar A = E/(1.0 - nu*nu);
  TacsScalar D = t*t*A/4.0;

  // Store the derivative of the stress values
  TacsScalar s[8];

  // Compute the in-plane resultants
  s[0] = A*(e[0] + nu*e[1]);
  s[1] = A*(e[1] + nu*e[0]);
  s[2] = G*e[2];

  // Compute the bending moments
  s[3] = D*(e[3] + nu*e[4]);
  s[4] = D*(e[4] + nu*e[3]);
  s[5] = 0.5*D*(1.0 - nu)*e[5];

  // Compute the shear resultants
  s[6] = kcorr*G*e[6];
  s[7] = kcorr*G*e[7];

  TacsScalar ksens = DRILLING_REGULARIZATION*G;

  // Add the result to the design variable vector
  fdvSens[tNum] +=
  (s[0]*psi[0] + s[1]*psi[1] + s[2]*psi[2] +
  s[3]*psi[3] + s[4]*psi[4] + s[5]*psi[5] +
  s[6]*psi[6] + s[7]*psi[7] + rotPsi*ksens);
  }
  }
*/

/*
  Compute the von Mises failure criterion on the upper and lower
  surfaces of the plate model
*/
/*
  void isoFSDTStiffness::failure( const double gpt[],
  const TacsScalar strain[],
  TacsScalar * fail ){
  TacsScalar stress[3];
  TacsScalar ht = 0.5*t;

  // Determine whether the failure will occur on the top or the bottom
  // Test the top of the plate
  calculatePlaneStress(stress, ht, strain);
  TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

  // Test the bottom of the plate
  calculatePlaneStress(stress, -ht, strain);
  TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

  *fail = (TacsRealPart(failTop) > TacsRealPart(failBot) ?
  failTop : failBot);
  }
*/
/*
  Compute the derivative of the von Mises failure criterion on the
  upper/lower surfaces with respect to the strain values
*/
/*
  void isoFSDTStiffness::failureStrainSens( const double gpt[],
  const TacsScalar strain[],
  TacsScalar sens[] ){
  TacsScalar stress[3];
  TacsScalar ht = 0.5*t;

  // Determine whether the failure will occur on the top or the bottom
  // Test the top of the plate
  calculatePlaneStress(stress, ht, strain);
  TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

  // Test the bottom of the plate
  calculatePlaneStress(stress, -ht, strain);
  TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

  if (TacsRealPart(failTop) > TacsRealPart(failBot)){
  // Test the top of the plate
  TacsScalar stressSens[3];
  calculatePlaneStress(stress, ht, strain);
  VonMisesFailurePlaneStressSens(stressSens, stress,
  yieldStress);

  calculatePlaneStressTranspose(sens, ht, stressSens);
  }
  else {
  // Test the bottom of the plate
  TacsScalar stressSens[3];
  calculatePlaneStress(stress, -ht, strain);
  VonMisesFailurePlaneStressSens(stressSens, stress,
  yieldStress);

  calculatePlaneStressTranspose(sens, -ht, stressSens);
  }
  }
*/
/*
  Add the derivative of the failure sensitivity on the upper and lower
  surfaces with respect to the design variable
*/
/*
  void isoFSDTStiffness::addFailureDVSens( const double pt[],
  const TacsScalar strain[],
  TacsScalar alpha,
  TacsScalar dvSens[], int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
  TacsScalar stress[3];
  TacsScalar ht = 0.5*t;

  // Determine whether the failure will occur on the top or the bottom
  // Test the top of the plate
  calculatePlaneStress(stress, ht, strain);
  TacsScalar failTop = VonMisesFailurePlaneStress(stress, yieldStress);

  // Test the bottom of the plate
  calculatePlaneStress(stress, -ht, strain);
  TacsScalar failBot = VonMisesFailurePlaneStress(stress, yieldStress);

  if (TacsRealPart(failTop) > TacsRealPart(failBot)){
  // Test the top of the plate
  TacsScalar stressSens[3];
  calculatePlaneStress(stress, ht, strain);
  VonMisesFailurePlaneStressSens(stressSens, stress,
  yieldStress);
  dvSens[tNum] +=
  0.5*alpha*(calculatePlaneStressTSensProduct(stressSens, strain));
  }
  else {
  TacsScalar stressSens[3];
  calculatePlaneStress(stress, -ht, strain);
  VonMisesFailurePlaneStressSens(stressSens, stress,
  yieldStress);
  dvSens[tNum] +=
  -0.5*alpha*(calculatePlaneStressTSensProduct(stressSens, strain));
  }
  }
  }
*/

/*
  Return the thickness as the design variable for this constitutive object
*/
/*
  TacsScalar isoFSDTStiffness::getDVOutputValue( int dv_index,
  const double pt[] ){
  if (dv_index == 0){
  return t;
  }
  if (dv_index == 1){
  return tNum;
  }
  return 0.0;
  }
*/
