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

#ifndef TACS_SOLID_STIFFNESS_H
#define TACS_SOLID_STIFFNESS_H

#include "TACSConstitutive.h"

/*
  Compute the stiffness matrix associated with a linear solid
  isotropic material.
*/
class SolidStiffness : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 6;
  SolidStiffness( TacsScalar _rho, TacsScalar _E,
                  TacsScalar _nu, TacsScalar _ys,
                  int _eNum=-1 );
  virtual ~SolidStiffness(){}

  // Functions for design variable control
  // -------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  
  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 1; }
  void getPointwiseMass( const double gpt[], TacsScalar mass[] ){
    mass[0] = rho;
  }

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress( const double pt[], 
                        const TacsScalar strain[],
                        TacsScalar stress[] );
  void addStressDVSens( const double pt[], const TacsScalar strain[],
                        TacsScalar alpha, const TacsScalar psi[],
                        TacsScalar dvSens[], int dvLen );

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
  
  // Extra info about the constitutive class
  // ---------------------------------------
  const char *constitutiveName();
  
 protected:
  SolidStiffness();

  inline void calcStress( const TacsScalar e[], TacsScalar s[] );

  // Mass, stiffness and Poisson ratio
  TacsScalar rho;
  TacsScalar E, nu;
  TacsScalar ys;

  // Design variables
  int eNum;

  // The stiffness parameters
  TacsScalar C[6];
  TacsScalar G23, G13, G12;

 private:
  static const char *constName;
};

inline void SolidStiffness::calcStress( const TacsScalar e[],
					TacsScalar s[] ){
  s[0] = C[0]*e[0] + C[1]*e[1] + C[2]*e[2];
  s[1] = C[1]*e[0] + C[3]*e[1] + C[4]*e[2];
  s[2] = C[2]*e[0] + C[4]*e[1] + C[5]*e[2];

  s[3] = G23*e[3];
  s[4] = G13*e[4];
  s[5] = G12*e[5];
}

#endif // TACS_SOLID_STIFFNESS_H
