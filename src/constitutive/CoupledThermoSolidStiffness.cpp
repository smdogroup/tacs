#include "CoupledThermoSolidStiffness.h"
/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * CoupledThermoSolidStiffness::constName = "CoupledThermoSolidStiffness";

const char * CoupledThermoSolidStiffness::constitutiveName(){ 
  return constName; 
}
/*
  CoupledThermoSolidStiffness member function definitions
*/
CoupledThermoSolidStiffness::CoupledThermoSolidStiffness( TacsScalar _rho, 
                                                          TacsScalar E,
                                                          TacsScalar nu,
                                                          TacsScalar _alpha,
                                                          TacsScalar _Tref,
                                                          TacsScalar _kcond ){
  // Density of the element
  rho = _rho;
  // Constitutive matrix for the structural model
  TacsScalar D = E/((1.0 + nu)*(1.0 - 2.0*nu));
  C[0] = C[3] = C[5] = (1.0 - nu)*D;
  C[1] = C[2] = C[4] = nu*D;

  G23 = G13 = G12 = 0.5*E/(1.0 + nu);
  
  // Constitutive matrix for the thermal analysis
  Tmat[0] = Tmat[2] = _kcond;
  Tmat[1] = _kcond;
  
  // Initialize the thermal parameter
  alpha = _alpha;
  Tref = _Tref;
  xw = 0.0;
}

CoupledThermoSolidStiffness::CoupledThermoSolidStiffness(){
  C[0] = C[1] = C[2] = 0.0;
  C[3] = C[4] = C[5] = 0.0;

  Tmat[0] = Tmat[1] = Tmat[2] = 0.0;

  rho = 0.0;
  alpha = 0.0;
  Tref = 0.0;
}

CoupledThermoSolidStiffness::~CoupledThermoSolidStiffness(){
  C[0] = C[1] = C[2] = C[3] = C[4] = C[5] = 0.0;
}
