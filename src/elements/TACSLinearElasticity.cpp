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

#include "TACSLinearElasticity.h"


TACSLinearEleasticity2D::TACSLinearEleasticity2D( PlaneStressStiffness *_con,
                                                  ElementStrainType _strain_type ){
  con = _con;
  con->incref();
  strain_type = _strain_type;
}

TACSLinearEleasticity2D::~TACSLinearEleasticity2D(){
  con->decref();
}
  
int TACSLinearEleasticity2D::getSpatialDim(){
  return 3;
}

int TACSLinearEleasticity2D::getVarsPerNode(){
  return 3;
}

void TACSLinearEleasticity2D::evalWeakIntegrand( const double time,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar U[],
                                                 const TacsScalar Udot[],
                                                 const TacsScalar Uddot[],
                                                 const TacsScalar Ux[],
                                                 TacsScalar DUt[],
                                                 TacsScalar DUx[] ){
  // Evaluate the density
  TacsScalar rho = con->evalDensity(pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Uddot[0];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Uddot[1];

  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[2]*Ux[2]);
    e[1] = Ux[3] + 0.5*(Ux[1]*Ux[1] + Ux[3]*Ux[3]);
    e[2] = Ux[1] + Ux[2] + (Ux[0]*Ux[1] + Ux[2]*Ux[3]);
  }

  // Evaluate the stress
  TacsScalar s[3];
  con->evalStress(pt, X, e, s);

  DUx[0] = 0.0;
  DUx[1] = s[0];
  DUx[2] = s[2];

  DUx[4] = 0.0;
  DUx[5] = s[2];
  DUx[6] = s[1];
}
  
void TACSLinearEleasticity2D::evalIntegrandDeriv( const double time,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar U[],
                                                  const TacsScalar Udot[],
                                                  const TacsScalar Uddot[],
                                                  const TacsScalar Ux[],
                                                  TacsScalar DUt[],
                                                  TacsScalar DUx[],
                                                  int *_DDt_nnz,
                                                  const int *_DDt_pairs[],
                                                  TacsScalar DDUt[],
                                                  int *_DDUx_nnz,
                                                  const int *_DDUx_pairs[],
                                                  TacsScalar DDUx[] ){

  // Evaluate the density
  TacsScalar rho = con->evalDensity(pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Uddot[0];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Uddot[1];

  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0];
    e[1] = Ux[3];
    e[2] = Ux[1] + Ux[2];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[2]*Ux[2]);
    e[1] = Ux[3] + 0.5*(Ux[1]*Ux[1] + Ux[3]*Ux[3]);
    e[2] = Ux[1] + Ux[2] + (Ux[0]*Ux[1] + Ux[2]*Ux[3]);
  }

  // Evaluate the stress
  TacsScalar s[3];
  con->evalStress(pt, X, e, s);

  DUx[0] = 0.0;
  DUx[1] = s[0];
  DUx[2] = s[2];

  DUx[4] = 0.0;
  DUx[5] = s[2];
  DUx[6] = s[1];

  TacsScalar C[36];
  con->evalTangentStiffness(pt, X, C);

  // Use a dense matrix
  _DDUx_num_non_zeros = -1;
  _DDUx_non_zero_pairs = NULL;

  if (strain_type == TACS_LINEAR_STRAIN){

  }
  else {

  }
}
  
TACSLinearEleasticity3D::TACSLinearEleasticity3D( TACSConstitutive *_con ){
  con = _con;
  con->incref();
}

TACSLinearEleasticity3D::~TACSLinearEleasticity3D(){
  con->decref();
}
  
int TACSLinearEleasticity3D::getSpatialDim(){
  return 3;
}

int TACSLinearEleasticity3D::getVarsPerNode(){
  return 3;
}

void TACSLinearEleasticity3D::evalWeakIntegrand( const double time,
                                                 const double pt[],
                                                 const TacsScalar X[],
                                                 const TacsScalar U[],
                                                 const TacsScalar Udot[],
                                                 const TacsScalar Uddot[],
                                                 const TacsScalar Ux[],
                                                 TacsScalar DUt[],
                                                 TacsScalar DUx[] ){
  // Evaluate the density
  TacsScalar rho = con->evalDensity(pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Uddot[0];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Uddot[1];

  DUt[6] = 0.0;
  DUt[7] = 0.0;
  DUt[8] = rho*Uddot[2];

  TacsScalar e[6];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[3]*Ux[3] + Ux[6]*Ux[6]);
    e[1] = Ux[4] + 0.5*(Ux[1]*Ux[1] + Ux[4]*Ux[4] + Ux[7]*Ux[7]);
    e[2] = Ux[8] + 0.5*(Ux[2]*Ux[2] + Ux[5]*Ux[5] + Ux[8]*Ux[8]);
    
    e[3] = Ux[5] + Ux[7] + (Ux[1]*Ux[2] + Ux[4]*Ux[5] + Ux[7]*Ux[8]);
    e[4] = Ux[2] + Ux[6] + (Ux[0]*Ux[2] + Ux[3]*Ux[5] + Ux[6]*Ux[8]);
    e[5] = Ux[1] + Ux[3] + (Ux[0]*Ux[1] + Ux[3]*Ux[4] + Ux[6]*Ux[7]);
  }

  // Evaluate the stress
  TacsScalar s[6];
  con->evalStress(pt, X, e, s);

  DUx[0] = 0.0;
  DUx[1] = s[0];
  DUx[2] = s[5];
  DUx[3] = s[4];

  DUx[4] = 0.0;
  DUx[5] = s[5];
  DUx[6] = s[1];
  DUx[7] = s[3];

  DUx[8] = 0.0;
  DUx[9] = s[4];
  DUx[10] = s[3];
  DUx[11] = s[2];
}

void TACSLinearEleasticity3D::evalIntegrandDeriv(const double time,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar U[],
                                                  const TacsScalar Udot[],
                                                  const TacsScalar Uddot[],
                                                  const TacsScalar Ux[],
                                                  TacsScalar DUt[],
                                                  TacsScalar DUx[],
                                                  int *_DDt_nnz,
                                                  const int *_DDt_pairs[],
                                                  TacsScalar DDUt[],
                                                  int *_DDUx_nnz,
                                                  const int *_DDUx_pairs[],
                                                  TacsScalar DDUx[] ){
  // Evaluate the density
  TacsScalar rho = con->evalDensity(pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Uddot[0];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Uddot[1];

  DUt[6] = 0.0;
  DUt[7] = 0.0;
  DUt[8] = rho*Uddot[2];

  *_DDt_num_non_zeros = 3;
  *_DDt_non_zero_pairs = DDt_non_zero_pairs;

  DDUt[0] = rho;
  DDUt[1] = rho;
  DDUt[2] = rho;

  TacsScalar e[6];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0];
    e[1] = Ux[4];
    e[2] = Ux[8];

    e[3] = Ux[5] + Ux[7];
    e[4] = Ux[2] + Ux[6];
    e[5] = Ux[1] + Ux[3];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[3]*Ux[3] + Ux[6]*Ux[6]);
    e[1] = Ux[4] + 0.5*(Ux[1]*Ux[1] + Ux[4]*Ux[4] + Ux[7]*Ux[7]);
    e[2] = Ux[8] + 0.5*(Ux[2]*Ux[2] + Ux[5]*Ux[5] + Ux[8]*Ux[8]);
    
    e[3] = Ux[5] + Ux[7] + (Ux[1]*Ux[2] + Ux[4]*Ux[5] + Ux[7]*Ux[8]);
    e[4] = Ux[2] + Ux[6] + (Ux[0]*Ux[2] + Ux[3]*Ux[5] + Ux[6]*Ux[8]);
    e[5] = Ux[1] + Ux[3] + (Ux[0]*Ux[1] + Ux[3]*Ux[4] + Ux[6]*Ux[7]);
  }

  // Evaluate the stress
  TacsScalar s[6];
  con->evalStress(pt, X, e, s);

  DUx[0] = 0.0;
  DUx[1] = s[0];
  DUx[2] = s[5];
  DUx[3] = s[4];

  DUx[4] = 0.0;
  DUx[5] = s[5];
  DUx[6] = s[1];
  DUx[7] = s[3];

  DUx[8] = 0.0;
  DUx[9] = s[4];
  DUx[10] = s[3];
  DUx[11] = s[2];

  TacsScalar C[36];
  con->evalTangentStiffness(pt, X, C);

  // Use a dense matrix
  _DDUx_num_non_zeros = -1;
  _DDUx_non_zero_pairs = NULL;

  if (strain_type == TACS_LINEAR_STRAIN){

  }
  else {

  }
}
