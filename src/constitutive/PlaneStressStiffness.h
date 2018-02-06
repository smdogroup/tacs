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

#ifndef TACS_PLANE_STRESS_STIFFNESS_H
#define TACS_PLANE_STRESS_STIFFNESS_H

#include "TACSConstitutive.h"

/*
  This is the base class for the plane stress constitutive objects. 
  
  All objects performing plane stress analysis should utilize this class. 
*/
class PlaneStressStiffness : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 3;
  PlaneStressStiffness();
  PlaneStressStiffness( TacsScalar _rho, TacsScalar E, TacsScalar nu );
  PlaneStressStiffness( TacsScalar _rho, TacsScalar E1, 
			TacsScalar E2, TacsScalar G12, TacsScalar nu12 );
  virtual ~PlaneStressStiffness(){}

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress( const double pt[], 
                        const TacsScalar strain[],
			TacsScalar stress[] );

  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 1; }
  void getPointwiseMass( const double pt[], TacsScalar mass[] ){
    mass[0] = rho;
  }

  // Extra info about the constitutive class
  // ---------------------------------------
  const char *constitutiveName();

 protected:
  // The stiffness matrix
  TacsScalar Cmat[6]; 
  TacsScalar rho;
 private:
  static const char *constName;
};

#endif
