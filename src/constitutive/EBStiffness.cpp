#include "EBStiffness.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.
*/

const char * EBStiffness::constName = "EBStiffness";
const char * EBStiffness::constitutiveName() {
  return constName;
}

int EBStiffness::ownsDesignVar( const int dvNum ) const {
  return 0;
}

int EBStiffness::getNumDesignVars() const {
  return 0;
}

int EBStiffness::getDesignVarNums( int * dvNums, int * dvIndex, int dvLen ) const {
  return 0;
}

EBStiffness::EBStiffness( TacsScalar _rho, TacsScalar _E, TacsScalar _G,
			  TacsScalar _A, TacsScalar _Ix, TacsScalar _Iy,
			  TacsScalar _J, TacsScalar _ref_dir[3],
			  enum EBBeamReferenceDirection _ref_type){

  // Save the input parameters
  rho = _rho;
  E = _E;
  G = _G;
  A = _A;
  Ix = _Ix;
  Iy = _Iy;
  J = _J;
  ref_dir[0] = _ref_dir[0];
  ref_dir[1] = _ref_dir[1];
  ref_dir[2] = _ref_dir[2];
  ref_type = _ref_type;

  for ( int k = 0; k < 10; k++ ){
    C[k] = 0.0;
  }

  C[0] = E*A;
  C[4] = E*Ix;
  C[7] = E*Iy;
  C[9] = G*J;

  for ( int k = 0; k < 6; k++ ){
    mass[k] = 0.0;
  }

  mass[0] = rho*A;
  mass[3] = rho*Ix;
  mass[5] = rho*Iy;
}

EBStiffness::EBStiffness(){
  for ( int k = 0; k < 10; k++ ){
    C[k] = 0.0;
  }

  for ( int k = 0; k < 6; k++ ){
    mass[k] = 0.0;
  }
}

int EBStiffness::getNumStresses() { return NUM_STRESSES; }

void EBStiffness::calculateStress( const double gpt[],
				   const TacsScalar strain[],
				   TacsScalar stress[] ){
  TacsScalar Ct[10];
  getStiffness(gpt, Ct);
  calcStress(Ct, strain, stress);
}

void EBStiffness::calculateStressDVSens( int dvNum,
                                         const double gpt[],
					 const TacsScalar strain[],
					 TacsScalar stress[] ){
  TacsScalar sCt[10];
  getStiffnessDVSens(dvNum, gpt, sCt);
  calcStress(sCt, strain, stress);
}


void EBStiffness::getStiffness( const double pt[], TacsScalar Ct[] ){
  memcpy(Ct, C, 10*sizeof(TacsScalar));
}

void EBStiffness::getStiffnessDVSens( int dvNum, const double pt[],
                                      TacsScalar Ct[] ){
  memset(Ct, 0, 10*sizeof(TacsScalar));
}

void EBStiffness::getPointwiseMass( const double pt[], TacsScalar _mass[] ){
  memcpy(_mass, mass, 6*sizeof(TacsScalar));
}

void EBStiffness::pointwiseMassDVSens( int dvNum, const double gpt[],
                                       TacsScalar massDVSens[] ){
  memset(massDVSens, 0, 6*sizeof(TacsScalar));
}
