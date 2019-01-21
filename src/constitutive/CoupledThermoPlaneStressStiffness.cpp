#include "CoupledThermoPlaneStressStiffness.h"
/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * CoupledThermoPlaneStressStiffness::constName = "CoupledThermoPlaneStressStiffness";

const char * CoupledThermoPlaneStressStiffness::constitutiveName(){ 
  return constName; 
}
/*
  CoupledThermoPlaneStressStiffness member function definitions
*/
CoupledThermoPlaneStressStiffness::CoupledThermoPlaneStressStiffness( TacsScalar _rho, 
                                                                      TacsScalar E,
                                                                      TacsScalar nu,
                                                                      TacsScalar _alpha,
                                                                      TacsScalar _Tref,
                                                                      TacsScalar _kcond
                                                                      ){
  // Density of the element
  rho = _rho;
  // Constitutive matrix for the structural model
  Cmat[0] = Cmat[3] = E/(1.0 - nu*nu);
  Cmat[5] = 0.5*E/(1.0 + nu);

  Cmat[1] = nu*Cmat[0];
  Cmat[2] = 0.0;
  Cmat[4] = 0.0;
  
  // Constitutive matrix for the thermal analysis
  Tmat[0] = Tmat[2] = _kcond;
  Tmat[1] = 0.0;
  
  // Initialize the thermal parameter
  alpha = _alpha;
  Tref = _Tref;
  xw = 0.0;
}

CoupledThermoPlaneStressStiffness::CoupledThermoPlaneStressStiffness(){
  Cmat[0] = Cmat[1] = Cmat[2] = 0.0;
  Cmat[3] = Cmat[4] = Cmat[5] = 0.0;

  Tmat[0] = Tmat[1] = Tmat[2] = 0.0;

  rho = 0.0;
  alpha = 0.0;
  Tref = 0.0;
}

CoupledThermoPlaneStressStiffness::~CoupledThermoPlaneStressStiffness(){
  Cmat[0] = Cmat[1] = Cmat[2] = Cmat[3] = Cmat[4] = Cmat[5] = 0.0;
}
