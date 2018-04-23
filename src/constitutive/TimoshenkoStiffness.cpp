/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#include "TimoshenkoStiffness.h"

const char *TimoshenkoStiffness::constName = "TimoshenkoStiffness";

const char *TimoshenkoStiffness::constitutiveName(){
  return constName;
}

/*
  Constructor for Timoshenko beam theory based constitutive object.
  EA               : axial stiffness
  EI22, EI33, EI23 : bending stiffness
  GJ               : torsional stiffness
  kG22, kG33, kG23 : shearing stiffness
  m00              : mass per unit span 
  m11, m22, m33    : moments of inertia such that m11 = m22 + m33
  xm2, xm3         : cross sectional center of mass location
  xc2, xc3         : cross sectional centroid
  xk2, xk3         : cross sectional shear center
  muS              : viscous damping coefficient
*/
TimoshenkoStiffness::TimoshenkoStiffness( const TacsScalar _axis[],
                                          TacsScalar EA, 
                                          TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
                                          TacsScalar GJ,
                                          TacsScalar kG22, TacsScalar kG33, TacsScalar kG23,
                                          TacsScalar m00,
                                          TacsScalar m11, TacsScalar m22, TacsScalar m33,
                                          TacsScalar xm2, TacsScalar xm3,
                                          TacsScalar xc2, TacsScalar xc3,
                                          TacsScalar xk2, TacsScalar xk3,
                                          TacsScalar muS ){
  // Set the reference axis and normalize it
  axis[0] = _axis[0];
  axis[1] = _axis[1];
  axis[2] = _axis[2];
  TacsScalar tmp = 1.0/sqrt(axis[0]*axis[0] + 
                            axis[1]*axis[1] + 
                            axis[2]*axis[2]);
  axis[0] *= tmp;
  axis[1] *= tmp;
  axis[2] *= tmp;

  // Set the entries of the stiffness matrix
  memset(C, 0, 36*sizeof(TacsScalar));

  // row 1 for axial force
  C[0] = EA;
  C[2] = xc3*EA;
  C[3] = -xc2*EA;

  // row 2 for twisting moment
  C[7] = GJ + xk2*xk2*kG33 + xk3*xk3*kG22 + 2.0*xk2*xk3*kG23; 
  C[10] = -xk2*kG23 - xk3*kG22; 
  C[11] = xk2*kG33 + xk3*kG23; 

  // row 3 for bending moment about axis 2
  C[12] = C[2];
  C[14] = EI22 + xc3*xc3*EA; 
  C[15] = -(EI23 + xc2*xc3*EA);
  
  // row 4 for bending moment about axis 3
  C[18] = C[3];
  C[20] = C[15];
  C[21] = EI33 + xc2*xc2*EA;   

  // row 5 for shear 2
  C[25] = C[10];
  C[28] = kG22;
  C[29] = -kG23;

  // row 6 for shear 3
  C[31] = C[11];
  C[34] = C[29];
  C[35] = kG33;
  
  // Set the entries of the density matrix
  rho[0] = m00;
  rho[1] = m11;
  rho[2] = m22;
  rho[3] = m00*xm2*xm3;
}

/*
  Set the diagonal components of the stiffness matrix and the mass
  moments of the cross-section.
*/
TimoshenkoStiffness::TimoshenkoStiffness( TacsScalar rhoA, 
                                          TacsScalar rhoIy,
                                          TacsScalar rhoIz, 
                                          TacsScalar rhoIyz,
                                          TacsScalar EA, 
                                          TacsScalar GJ,
                                          TacsScalar EIy, 
                                          TacsScalar EIz,
                                          TacsScalar kGAy, 
                                          TacsScalar kGAz,
                                          const TacsScalar _axis[] ){
  // Set the reference axis and normalize it
  axis[0] = _axis[0];
  axis[1] = _axis[1];
  axis[2] = _axis[2];
  TacsScalar tmp = 1.0/sqrt(axis[0]*axis[0] + 
                            axis[1]*axis[1] + 
                            axis[2]*axis[2]);
  axis[0] *= tmp;
  axis[1] *= tmp;
  axis[2] *= tmp;

  // Set the entries of the stiffness matrix
  memset(C, 0, 36*sizeof(TacsScalar));
  C[0] = EA;
  C[7] = GJ;
  C[14] = EIy;
  C[21] = EIz;
  C[28] = kGAy;
  C[35] = kGAz;

  // Set the entries of the density matrix
  rho[0] = rhoA;
  rho[1] = rhoIy;
  rho[2] = rhoIz;
  rho[3] = rhoIyz;
}

/*
  Set the full stiffness matrix
*/
TimoshenkoStiffness::TimoshenkoStiffness( const TacsScalar _rho[], 
                                          const TacsScalar _C[],
                                          const TacsScalar _axis[] ){
  setData(_rho, _C, _axis);
}

/*
  Set the stiffness data
*/
void TimoshenkoStiffness::setData( const TacsScalar _rho[],
                                   const TacsScalar _C[],
                                   const TacsScalar _axis[] ){
  // Set the reference axis and normalize it
  axis[0] = _axis[0];
  axis[1] = _axis[1];
  axis[2] = _axis[2];
  TacsScalar tmp = 1.0/sqrt(axis[0]*axis[0] + 
                            axis[1]*axis[1] + 
                            axis[2]*axis[2]);
  axis[0] *= tmp;
  axis[1] *= tmp;
  axis[2] *= tmp;

  // Copy the density/stiffness matrix
  memcpy(rho, _rho, 4*sizeof(TacsScalar));
  memcpy(C, _C, 36*sizeof(TacsScalar));
}

TimoshenkoStiffness::~TimoshenkoStiffness(){}

/*
  Get the number of stress components
*/
int TimoshenkoStiffness::getNumStresses(){
  return 6;
}

/*
  Compute the stress, given the strain
*/
void TimoshenkoStiffness::calculateStress( const double pt[], 
                                           const TacsScalar strain[],
                                           TacsScalar stress[] ){
  calcStress(strain, stress);
}
