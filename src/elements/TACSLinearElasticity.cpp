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

TACSLinearElasticity2D::TACSLinearElasticity2D( TACSPlaneStressConstitutive *_stiff,
                                                ElementStrainType _strain_type ){
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
}

TACSLinearElasticity2D::~TACSLinearElasticity2D(){
  stiff->decref();
}

const int TACSLinearElasticity2D::DDUt_pairs[] =
  {2, 2, 5, 5};

const int TACSLinearElasticity2D::DDUx_pairs[] =
  {1, 1, 1, 2, 1, 4, 1, 5,
   2, 1, 2, 2, 2, 4, 2, 5,
   4, 1, 4, 2, 4, 4, 4, 5,
   5, 1, 5, 2, 5, 4, 5, 5};

int TACSLinearElasticity2D::getSpatialDim(){
  return 3;
}

int TACSLinearElasticity2D::getVarsPerNode(){
  return 2;
}

void TACSLinearElasticity2D::evalWeakIntegrand( int elemIndex,
                                                const double time,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar U[],
                                                const TacsScalar Udot[],
                                                const TacsScalar Uddot[],
                                                const TacsScalar Ux[],
                                                TacsScalar DUt[],
                                                TacsScalar DUx[] ){
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

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
  stiff->evalStress(elemIndex, pt, X, e, s);

  DUx[0] = 0.0;
  DUx[1] = s[0];
  DUx[2] = s[2];

  DUx[4] = 0.0;
  DUx[5] = s[2];
  DUx[6] = s[1];
}

void TACSLinearElasticity2D::evalWeakJacobian( int elemIndex,
                                               const double time,
                                               const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar U[],
                                               const TacsScalar Udot[],
                                               const TacsScalar Uddot[],
                                               const TacsScalar Ux[],
                                               TacsScalar DUt[],
                                               TacsScalar DUx[],
                                               int *DDUt_nnz,
                                               const int *_DDUt_pairs[],
                                               TacsScalar DDUt[],
                                               int *DDUx_nnz,
                                               const int *_DDUx_pairs[],
                                               TacsScalar DDUx[] ){
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

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
  stiff->evalStress(elemIndex, pt, X, e, s);

  DUx[0] = 0.0;  // u
  DUx[1] = s[0]; // u,x
  DUx[2] = s[2]; // u,y

  DUx[4] = 0.0;  // v
  DUx[5] = s[2]; // v,x
  DUx[6] = s[1]; // v,y

  TacsScalar C[21];
  stiff->evalTangentStiffness(elemIndex, pt, X, C);

  // Set the non-zero terms in the Jacobian
  *DDUt_nnz = 2;
  *_DDUt_pairs = DDUt_pairs;
  *DDUx_nnz = 16;
  *_DDUx_pairs = DDUx_pairs;

  // Set the acceleration terms
  DDUt[0] = rho;
  DDUt[1] = rho;

  // s = C*e
  if (strain_type == TACS_LINEAR_STRAIN){
    // Index:       1            5            2     4
    // s[0] = C[0]*(u,x) + C[1]*(v,y) + C[2]*(u,y + v,x)
    // s[1] = C[1]*(u,x) + C[3]*(v,y) + C[4]*(u,y + v,x)
    // s[2] = C[2]*(u,x) + C[4]*(v,y) + C[5]*(u,y + v,x)

    // i == 1 (s[0])
    DDUx[0] = C[0]; // j == 1
    DDUx[1] = C[2]; // j == 2
    DDUx[2] = C[2]; // j == 4
    DDUx[3] = C[1]; // j == 5

    // i == 2 (s[2])
    DDUx[4] = C[2]; // j == 1
    DDUx[5] = C[5]; // j == 2
    DDUx[6] = C[5]; // j == 4
    DDUx[7] = C[4]; // j == 5

    // i == 4 (s[2])
    DDUx[8] = C[2]; // j == 1
    DDUx[9] = C[5]; // j == 2
    DDUx[10] = C[5]; // j == 4
    DDUx[11] = C[4]; // j == 5

    // i == 5 (s[1])
    DDUx[12] = C[1]; // j == 1
    DDUx[13] = C[4]; // j == 2
    DDUx[14] = C[4]; // j == 4
    DDUx[15] = C[3]; // j == 5
  }
  else {
    // i == 1 (s[0])
    DDUx[0] = C[0]; // j == 1
    DDUx[1] = C[2]; // j == 2
    DDUx[2] = C[2]; // j == 4
    DDUx[3] = C[1]; // j == 5

    // i == 2 (s[2])
    DDUx[4] = C[2]; // j == 1
    DDUx[5] = C[5]; // j == 2
    DDUx[6] = C[5]; // j == 4
    DDUx[7] = C[4]; // j == 5

    // i == 4 (s[2])
    DDUx[8] = C[2]; // j == 1
    DDUx[9] = C[5]; // j == 2
    DDUx[10] = C[5]; // j == 4
    DDUx[11] = C[4]; // j == 5

    // i == 5 (s[1])
    DDUx[12] = C[1]; // j == 1
    DDUx[13] = C[4]; // j == 2
    DDUx[14] = C[4]; // j == 4
    DDUx[15] = C[3]; // j == 5
  }
}

/*
TACSLinearElasticity3D::TACSLinearElasticity3D( TACSConstitutive *_con ){
  con = _con;
  con->incref();
}

TACSLinearElasticity3D::~TACSLinearElasticity3D(){
  con->decref();
}

int TACSLinearElasticity3D::getSpatialDim(){
  return 3;
}

int TACSLinearElasticity3D::getVarsPerNode(){
  return 3;
}

void TACSLinearElasticity3D::evalWeakIntegrand( const double time,
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

void TACSLinearElasticity3D::evalIntegrandDeriv(const double time,
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
*/
