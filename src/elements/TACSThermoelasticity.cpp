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

#include "TACSThermoelasticity.h"

TACSLinearThermoelasticity2D::TACSLinearThermoelasticity2D( TACSPlaneStressConstitutive *_stiff,
                                                            ElementStrainType _strain_type ){
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
}

TACSLinearThermoelasticity2D::~TACSLinearThermoelasticity2D(){
  stiff->decref();
}

// 0;   1;    2;   3;   4; 5;   6;    7;   8;   9;10;  11;   12;  13;  14
// u; u,t; u,tt; u,x; u,y; v; v,t; v,tt; v,x; v,y; T; T,t; T,tt; T,x; T,y

const int TACSLinearThermoelasticity2D::linear_Jac_pairs[] =
  {2, 2, 7, 7,
   3, 3, 3, 4, 3, 8, 3, 9,
   4, 3, 4, 4, 4, 8, 4, 9,
   8, 3, 8, 4, 8, 8, 8, 9,
   9, 3, 9, 4, 9, 8, 9, 9};

int TACSLinearThermoelasticity2D::getSpatialDim(){
  return 2;
}

int TACSLinearThermoelasticity2D::getVarsPerNode(){
  return 3;
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearThermoelasticity2D::getDesignVarNums( int elemIndex, int dvLen,
                                                    int dvNums[] ){
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
void TACSLinearThermoelasticity2D::setDesignVars( int elemIndex, int dvLen,
                                                  const TacsScalar dvs[] ){
  stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
void TACSLinearThermoelasticity2D::getDesignVars( int elemIndex, int dvLen,
                                                  TacsScalar dvs[] ){
  stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
void TACSLinearThermoelasticity2D::getDesignVarRange( int elemIndex, int dvLen,
                                                      TacsScalar lb[],
                                                      TacsScalar ub[] ){
  stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearThermoelasticity2D::evalWeakIntegrand( int elemIndex,
                                                      int n,
                                                      const double time,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      const TacsScalar Ut[],
                                                      const TacsScalar Ux[],
                                                      TacsScalar DUt[],
                                                      TacsScalar DUx[] ){
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Ut[2];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Ut[5];

  DUt[6] = 0.0;
  DUt[7] = 0.0;
  DUt[8] = rho*Cp*Ut[8];

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

  DUx[0] = s[0];
  DUx[1] = s[2];

  DUx[2] = s[2];
  DUx[3] = s[1];
}

void TACSLinearThermoelasticity2D::evalWeakJacobian( int elemIndex,
                                                     int n,
                                                     const double time,
                                                     const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar Ut[],
                                                     const TacsScalar Ux[],
                                                     TacsScalar DUt[],
                                                     TacsScalar DUx[],
                                                     int *Jac_nnz,
                                                     const int *Jac_pairs[],
                                                     TacsScalar Jac[] ){
  // Evaluate the density
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Ut[2];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Ut[5];

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

  DUx[0] = s[0]; // u,x
  DUx[1] = s[2]; // u,y

  DUx[2] = s[2]; // v,x
  DUx[3] = s[1]; // v,y

  TacsScalar C[6];
  stiff->evalTangentStiffness(elemIndex, pt, X, C);

  // Set the non-zero terms in the Jacobian
  *Jac_nnz = 18;
  *Jac_pairs = linear_Jac_pairs;

  // Set the acceleration terms
  Jac[0] = rho;
  Jac[1] = rho;

  // s = C*e
  if (strain_type == TACS_LINEAR_STRAIN){
    // Index:         3            9            4     8
    // s[0] = C[0]*(u,x) + C[1]*(v,y) + C[2]*(u,y + v,x)
    // s[1] = C[1]*(u,x) + C[3]*(v,y) + C[4]*(u,y + v,x)
    // s[2] = C[2]*(u,x) + C[4]*(v,y) + C[5]*(u,y + v,x)

    // i == 1 (s[0])
    Jac[2] = C[0]; // j == 3
    Jac[3] = C[2]; // j == 4
    Jac[4] = C[2]; // j == 8
    Jac[5] = C[1]; // j == 9

    // i == 2 (s[2])
    Jac[6] = C[2]; // j == 3
    Jac[7] = C[5]; // j == 4
    Jac[8] = C[5]; // j == 8
    Jac[9] = C[4]; // j == 9

    // i == 4 (s[2])
    Jac[10] = C[2]; // j == 3
    Jac[11] = C[5]; // j == 4
    Jac[12] = C[5]; // j == 8
    Jac[13] = C[4]; // j == 9

    // i == 5 (s[1])
    Jac[14] = C[1]; // j == 3
    Jac[15] = C[4]; // j == 4
    Jac[16] = C[4]; // j == 8
    Jac[17] = C[3]; // j == 9
  }
}

/*
  Get the data for visualization at a given point
*/
void TACSLinearThermoelasticity2D::getOutputData( int elemIndex,
                                                  const double time,
                                                  ElementType etype,
                                                  int write_flag,
                                                  const double pt[],
                                                  const TacsScalar X[],
                                                  const TacsScalar Ut[],
                                                  const TacsScalar Ux[],
                                                  int ld_data,
                                                  TacsScalar *data ){
  if (etype == TACS_PLANE_STRESS_ELEMENT){
    if (write_flag & TACS_OUTPUT_NODES){
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
      data[0] = Ut[0];
      data[1] = Ut[3];
      data += 2;
    }

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

    if (write_flag & TACS_OUTPUT_STRAINS){
      data[0] = e[0];
      data[1] = e[1];
      data[2] = e[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_STRESSES){
      TacsScalar s[3];
      stiff->evalStress(elemIndex, pt, X, e, s);
      data[0] = s[0];
      data[1] = s[1];
      data[2] = s[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS){
      data[0] = stiff->failure(elemIndex, pt, X, e);
      data[1] = 0.0;
      data[2] = 0.0;
      data += 3;
    }
  }
}
