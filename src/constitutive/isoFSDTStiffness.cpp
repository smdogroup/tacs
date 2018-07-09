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

#include <stdio.h>
#include "isoFSDTStiffness.h"
#include "YSlibrary.h"

const char * isoFSDTStiffness::constName = "isoFSDTStiffness";

isoFSDTStiffness::isoFSDTStiffness( TacsScalar _rho, TacsScalar _E,
                                    TacsScalar _nu, TacsScalar _kcorr,
                                    TacsScalar _yieldStress,
                                    TacsScalar _thickness, int _tNum,
                                    TacsScalar _minThickness,
                                    TacsScalar _maxThickness ){
  // Set the material properties
  rho = _rho;
  E = _E;
  nu = _nu;
  G = 0.5*E/(1.0 + nu);
  kcorr = _kcorr;
  yieldStress = _yieldStress;

  // Copy over the thickness data
  t = _thickness;
  tNum = _tNum;
  minThickness = _minThickness;
  maxThickness = _maxThickness;
}

/*
  Set the design variable values from the vector
*/
void isoFSDTStiffness::setDesignVars( const TacsScalar dvs[], int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
    t = dvs[tNum];
  }
}

/*
  Get the design variable value from the object
*/
void isoFSDTStiffness::getDesignVars( TacsScalar dvs[], int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
    dvs[tNum] = t;
  }
}

/*
  Get the lower and upper bounds for the design variable values
*/
void isoFSDTStiffness::getDesignVarRange( TacsScalar lowerBound[],
                                          TacsScalar upperBound[],
                                          int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
    lowerBound[tNum] = minThickness;
    upperBound[tNum] = maxThickness;
  }
}

/*
  Get the areal mass and the moment of mass per unit area
*/
void isoFSDTStiffness::getPointwiseMass( const double pt[],
                                         TacsScalar mass[] ){
  mass[0] = rho*t;
  mass[1] = (rho*t*t*t)/12.0;
}

/*
  Add the derivative of mass to the given design variable vector
*/
void isoFSDTStiffness::addPointwiseMassDVSens( const double pt[],
                                               const TacsScalar alpha[],
                                               TacsScalar fdvSens[], int dvLen ){
  if (tNum >= 0 && tNum < dvLen){
    fdvSens[tNum] += rho*(alpha[0] + 0.25*alpha[1]*t*t);
  }
}

/*
  Compute the stiffness at the given parametric location within the
  element

  Note that the ABD matrices have the following layout:

  A = [ A[0] A[1] A[2] ]
  .   [ A[1] A[3] A[4] ]
  .   [ A[2] A[4] A[5] ]
*/
TacsScalar isoFSDTStiffness::getStiffness( const double pt[],
                                           TacsScalar A[], TacsScalar B[],
                                           TacsScalar D[], TacsScalar As[] ){
  TacsScalar a = t*E/(1.0 - nu*nu);
  TacsScalar d = t*t*a/12.0;

  // Compute the in-plane stiffness
  A[2] = A[4] = 0.0;
  A[0] = A[3] = a;
  A[1] = nu*a;
  A[5] = G*t;

  // Set the bending-stretching coupling to zero
  B[0] = B[1] = B[2] = B[3] = B[4] = B[5] = 0.0;

  // Compute the bending moments
  D[2] = D[4] = 0.0;
  D[0] = D[3] = d;
  D[1] = nu*d;
  D[5] = 0.5*d*(1.0 - nu);

  // Compute the shear resultants
  As[0] = kcorr*G*t;
  As[1] = 0.0;
  As[2] = kcorr*G*t;

  return DRILLING_REGULARIZATION*G*t;
}

/*
  Add the derivative of the product of the stress to the design
  derivative array
*/
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

/*
  Compute the von Mises failure criterion on the upper and lower
  surfaces of the plate model
*/
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

/*
  Compute the derivative of the von Mises failure criterion on the
  upper/lower surfaces with respect to the strain values
*/
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

/*
  Add the derivative of the failure sensitivity on the upper and lower
  surfaces with respect to the design variable
*/
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

/*
  Return the thickness as the design variable for this constitutive object
*/
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
