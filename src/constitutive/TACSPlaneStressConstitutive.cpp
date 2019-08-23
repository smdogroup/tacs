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

#include "PlaneStressStiffness.h"

const char * PlaneStressStiffness::constName = "PlaneStressStiffness";

const char * PlaneStressStiffness::constitutiveName(){ 
  return constName; 
}

/*
  PlaneStressStiffness member function definitions
*/
PlaneStressStiffness::PlaneStressStiffness( TacsScalar _rho, TacsScalar E,
					    TacsScalar nu ){
  rho = _rho;

  Cmat[0] = Cmat[3] = E/(1.0 - nu*nu);
  Cmat[5] = 0.5*E/(1.0 + nu);

  Cmat[1] = nu*Cmat[0];
  Cmat[2] = 0.0;
  Cmat[4] = 0.0;
}

PlaneStressStiffness::PlaneStressStiffness( TacsScalar _rho, TacsScalar E1,
					    TacsScalar E2, TacsScalar G12, 
					    TacsScalar nu12 ){
  rho = _rho;

  TacsScalar nu21 = E2/E1*nu12;
  TacsScalar S = (1.0 - nu12*nu21);

  Cmat[0] = E1/S;
  Cmat[3] = E2/S;
  Cmat[5] = G12;

  Cmat[1] = nu21 * E1/S;
  Cmat[2] = 0.0;
  Cmat[4] = 0.0;
}

PlaneStressStiffness::PlaneStressStiffness(){
  Cmat[0] = Cmat[1] = Cmat[2] = 0.0;
  Cmat[3] = Cmat[4] = Cmat[5] = 0.0;

  rho = 0.0;
}

int PlaneStressStiffness::getNumStresses(){ return NUM_STRESSES; }

void PlaneStressStiffness::calculateStress( const double pt[], 
					    const TacsScalar strain[],
					    TacsScalar stress[] ){
  stress[0] = Cmat[0]*strain[0] + Cmat[1]*strain[1] + Cmat[2]*strain[2];
  stress[1] = Cmat[1]*strain[0] + Cmat[3]*strain[1] + Cmat[4]*strain[2];
  stress[2] = Cmat[2]*strain[0] + Cmat[4]*strain[1] + Cmat[5]*strain[2];
}
