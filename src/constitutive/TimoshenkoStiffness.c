#include "TimoshenkoStiffness.h"

const char *TimoshenkoStiffness::constName = "TimoshenkoStiffness";

const char *TimoshenkoStiffness::constitutiveName(){
  return constName;
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
