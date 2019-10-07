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

#ifndef TACS_ISO_FSDT_STIFFNESS_H
#define TACS_ISO_FSDT_STIFFNESS_H

#include "FSDTStiffness.h"

/*
  This object defines the stiffness for a plate or shell. It is
  derived from the FSDTStiffness class.

  This particular variable thickness constitutive relationship defines
  the necessary methods for calculating a variable thickness plate.
*/
class isoFSDTStiffness : public FSDTStiffness {
 public:
  static const int NUM_STRESSES = FSDTStiffness::NUM_STRESSES;

  isoFSDTStiffness( TacsScalar _rho, TacsScalar _E, TacsScalar _nu,
                    TacsScalar _kcorr, TacsScalar _yieldStress,
                    TacsScalar _thickness, int _tNum=-1,
                    TacsScalar _minThickness=1e-3,
                    TacsScalar _maxThickness=50.0 );
  ~isoFSDTStiffness(){}

  // Functions for design variable control
  // -------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lowerBound[],
                          TacsScalar upperBound[], int numDVs );

  // Functions required by TACSConstitutive
  // --------------------------------------
  void getPointwiseMass( const double pt[], TacsScalar mass[] );
  void addPointwiseMassDVSens( const double pt[],
                               const TacsScalar alpha[],
                               TacsScalar dvSens[], int dvLen );

  // Functions required by FSDTStiffness
  // -----------------------------------
  TacsScalar getStiffness( const double pt[],
                           TacsScalar A[], TacsScalar B[],
                           TacsScalar D[], TacsScalar As[] );
  void addStiffnessDVSens( const double pt[],
                           const TacsScalar e[], const TacsScalar psi[],
                           TacsScalar rotPsi,
                           TacsScalar fdvSens[], int dvLen );

  // Functions to compute the failure properties
  // -------------------------------------------
  void failure( const double pt[],
                const TacsScalar strain[],
                TacsScalar *fail );
  void failureStrainSens( const double pt[],
                          const TacsScalar strain[],
                          TacsScalar sens[] );
  void addFailureDVSens( const double pt[],
                         const TacsScalar strain[],
                         TacsScalar alpha,
                         TacsScalar dvSens[], int dvLen );

  // Retrieve the design variable for plotting purposes
  // --------------------------------------------------
  TacsScalar getDVOutputValue( int dv_index, const double pt[] );

  // Return the name of the constitutive object
  // ------------------------------------------
  const char *constitutiveName(){ return constName; }

 private:
  // Calculate the state of plane stress within the element
  // ------------------------------------------------------
  inline void calculatePlaneStress( TacsScalar pstr[], const TacsScalar ht,
                                    const TacsScalar strain[] ){
    const TacsScalar D = E/(1.0 - nu*nu);
    pstr[0] = D*((strain[0] + ht*strain[3]) + nu*(strain[1] + ht*strain[4]));
    pstr[1] = D*(nu*(strain[0] + ht*strain[3]) + (strain[1] + ht*strain[4]));
    pstr[2] = G*(strain[2] + ht*strain[5]);
  }

  inline void calculatePlaneStressTranspose( TacsScalar strain[],
                                             const TacsScalar ht,
                                             const TacsScalar pstr[] ){
    const TacsScalar D = E/(1.0 - nu*nu);
    strain[0] = D*(pstr[0] + nu*pstr[1]);
    strain[1] = D*(nu*pstr[0] + pstr[1]);
    strain[2] = G*pstr[2];

    strain[3] = ht*D*(pstr[0] + nu*pstr[1]);
    strain[4] = ht*D*(nu*pstr[0] + pstr[1]);
    strain[5] = G*ht*pstr[2];

    strain[6] = 0.0;
    strain[7] = 0.0;
  }

  TacsScalar calculatePlaneStressTSensProduct( const TacsScalar pstr[],
                                               const TacsScalar strain[] ){
    return E/(1.0 - nu*nu)*(pstr[0]*(strain[3] + nu*strain[4]) +
                            pstr[1]*(nu*strain[3] + strain[4])) +
      G*pstr[2]*strain[5];
  }

  // The material properties required for this object
  TacsScalar kcorr; // The shear correction factor
  TacsScalar rho; // The material density
  TacsScalar E, nu, G; // The stiffness properties
  TacsScalar yieldStress; // The material yield stress

  // The thickness information required by this object
  int tNum; // The design variable number
  TacsScalar t; // The thickness
  TacsScalar minThickness, maxThickness; // The min/max thickness values

  // The name of this object
  static const char * constName;
};

#endif // TACS_ISO_FSDT_STIFFNESS_H
