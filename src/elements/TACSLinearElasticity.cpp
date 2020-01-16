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

// 0;   1;    2;   3;   4; 5;   6;    7;   8;   9;
// u; u,t; u,tt; u,x; u,y; v; v,t; v,tt; v,x; v,y;

const int TACSLinearElasticity2D::linear_Jac_pairs[] =
  {2, 2, 7, 7,
   3, 3, 3, 4, 3, 8, 3, 9,
   4, 3, 4, 4, 4, 8, 4, 9,
   8, 3, 8, 4, 8, 8, 8, 9,
   9, 3, 9, 4, 9, 8, 9, 9};

int TACSLinearElasticity2D::getSpatialDim(){
  return 2;
}

int TACSLinearElasticity2D::getVarsPerNode(){
  return 2;
}

int TACSLinearElasticity2D::getDesignVarsPerNode(){
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearElasticity2D::getDesignVarNums( int elemIndex, int dvLen,
                                              int dvNums[] ){
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearElasticity2D::setDesignVars( int elemIndex, int dvLen,
                                           const TacsScalar dvs[] ){
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearElasticity2D::getDesignVars( int elemIndex, int dvLen,
                                           TacsScalar dvs[] ){
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearElasticity2D::getDesignVarRange( int elemIndex, int dvLen,
                                               TacsScalar lb[],
                                               TacsScalar ub[] ){
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearElasticity2D::evalWeakIntegrand( int elemIndex,
                                                const double time,
                                                int n,
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

void TACSLinearElasticity2D::evalWeakJacobian( int elemIndex,
                                               const double time,
                                               int n,
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
  DUt[2] = rho*Ut[2]; // u,tt

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Ut[5]; // v,tt

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
  Add the design variable derivative of the product of the adjoint
  vector with the weak form of the residual
*/
void TACSLinearElasticity2D::addWeakAdjProduct( int elemIndex,
                                                const double time,
                                                int n,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar Ut[],
                                                const TacsScalar Ux[],
                                                const TacsScalar Psi[],
                                                const TacsScalar Psix[],
                                                TacsScalar scale,
                                                int dvLen,
                                                TacsScalar *fdvSens ){
  // Evaluate the density
  TacsScalar rho_coef = scale*(Ut[2]*Psi[0] + Ut[5]*Psi[1]);
  stiff->addDensityDVSens(elemIndex, pt, X, rho_coef, dvLen, fdvSens);

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

  TacsScalar phi[3];
  phi[0] = Psix[0];
  phi[1] = Psix[3];
  phi[2] = Psix[1] + Psix[2];
  stiff->addStressDVSens(elemIndex, pt, X, e, scale, phi, dvLen, fdvSens);
}

void TACSLinearElasticity2D::evalWeakAdjXptSensProduct( int elemIndex,
                                                        const double time,
                                                        int n,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        const TacsScalar Ut[],
                                                        const TacsScalar Ux[],
                                                        const TacsScalar Psi[],
                                                        const TacsScalar Psix[],
                                                        TacsScalar *product,
                                                        TacsScalar dfdX[],
                                                        TacsScalar dfdUx[],
                                                        TacsScalar dfdPsix[] ){
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;
  dfdUx[0] = dfdUx[1] = dfdUx[2] = dfdUx[3] = 0.0;
  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = dfdPsix[3] = 0.0;

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

  TacsScalar phi[3];
  phi[0] = Psix[0];
  phi[1] = Psix[3];
  phi[2] = Psix[1] + Psix[2];

  TacsScalar t1[3], t2[3];
  stiff->evalStress(elemIndex, pt, X, e, t1);
  stiff->evalStress(elemIndex, pt, X, phi, t2);

  *product = t2[0]*e[0] + t2[1]*e[1] + t2[2]*e[2];

  if (strain_type == TACS_LINEAR_STRAIN){
    dfdPsix[0] = t1[0];
    dfdPsix[3] = t1[1];
    dfdPsix[1] = t1[2];
    dfdPsix[2] = t1[2];

    dfdUx[0] = t2[0];
    dfdUx[3] = t2[1];
    dfdUx[1] = t2[2];
    dfdUx[2] = t2[2];
  }
}

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearElasticity2D::evalPointQuantity( int elemIndex,
                                               const int quantityType,
                                               const double time,
                                               int n, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar Xd[],
                                               const TacsScalar Ut[],
                                               const TacsScalar Ux[],
                                               TacsScalar *quantity ){
  if (quantityType == TACS_FAILURE_INDEX){
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

    *quantity = stiff->evalFailure(elemIndex, pt, X, e);

    return 1;
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    *quantity = stiff->evalDensity(elemIndex, pt, X);

    return 1;
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    // Evaluate the strain energy density
    *quantity = (e[0]*s[0] + e[1]*s[1] + e[2]*s[2]);

    return 1;
  }

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSLinearElasticity2D::addPointQuantityDVSens( int elemIndex,
                                                     const int quantityType,
                                                     const double time,
                                                     TacsScalar scale,
                                                     int n, const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar Xd[],
                                                     const TacsScalar Ut[],
                                                     const TacsScalar Ux[],
                                                     const TacsScalar dfdq[],
                                                     int dvLen,
                                                     TacsScalar dfdx[] ){
  if (quantityType == TACS_FAILURE_INDEX){
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

    stiff->addFailureDVSens(elemIndex, pt, X, e, scale*dfdq[0],
                            dvLen, dfdx);
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    stiff->addDensityDVSens(elemIndex, pt, X, scale*dfdq[0], dvLen, dfdx);
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);
    stiff->addStressDVSens(elemIndex, pt, X, e, scale*dfdq[0],
                           e, dvLen, dfdx);
  }
}

/*
  Evaluate the derivatives of the point-wise quantity of interest
  with respect to X, Ut and Ux.
*/
void TACSLinearElasticity2D::evalPointQuantitySens( int elemIndex,
                                                    const int quantityType,
                                                    const double time,
                                                    int n, const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar Xd[],
                                                    const TacsScalar Ut[],
                                                    const TacsScalar Ux[],
                                                    const TacsScalar dfdq[],
                                                    TacsScalar dfdX[],
                                                    TacsScalar dfdXd[],
                                                    TacsScalar dfdUt[],
                                                    TacsScalar dfdUx[] ){
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = 0.0;
  dfdXd[2] = dfdXd[3] = 0.0;

  dfdUt[0] = dfdUt[1] = dfdUt[2] = 0.0;
  dfdUt[3] = dfdUt[4] = dfdUt[5] = 0.0;

  dfdUx[0] = dfdUx[1] = 0.0;
  dfdUx[2] = dfdUx[3] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX){
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

    TacsScalar sens[3];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    if (strain_type == TACS_LINEAR_STRAIN){
      dfdUx[0] = dfdq[0]*sens[0];
      dfdUx[3] = dfdq[0]*sens[1];

      dfdUx[1] = dfdq[0]*sens[2];
      dfdUx[2] = dfdq[0]*sens[2];
    }
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[3];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN){
      dfdUx[0] = 2.0*dfdq[0]*s[0];
      dfdUx[3] = 2.0*dfdq[0]*s[1];

      dfdUx[1] = 2.0*dfdq[0]*s[2];
      dfdUx[2] = 2.0*dfdq[0]*s[2];
    }
  }
}

/*
  Get the data for visualization at a given point
*/
void TACSLinearElasticity2D::getOutputData( int elemIndex,
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
      data[0] = stiff->evalFailure(elemIndex, pt, X, e);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 4;
    }
  }
}

TACSLinearElasticity3D::TACSLinearElasticity3D( TACSSolidConstitutive *_stiff,
                                                ElementStrainType _strain_type ){
  stiff = _stiff;
  stiff->incref();
  strain_type = _strain_type;
}

TACSLinearElasticity3D::~TACSLinearElasticity3D(){
  stiff->decref();
}

// 0;   1;    2;   3;   4;   5;
// u; u,t; u,tt; u,x; u,y; u,z;

// 6;   7;    8;   9;  10;  11;
// v; v,t; v,tt; v,x; v,y; v,z;

//12;  13;   14;  15;  16;  17;
// w; w,t; w,tt; w,x; w,y; w,z;

const int TACSLinearElasticity3D::linear_Jac_pairs[] =
  {2, 2, 8, 8, 14, 14,
   3, 3, 3, 4, 3, 5, 3, 9, 3, 10, 3, 11, 3, 15, 3, 16, 3, 17,
   4, 3, 4, 4, 4, 5, 4, 9, 4, 10, 4, 11, 4, 15, 4, 16, 4, 17,
   5, 3, 5, 4, 5, 5, 5, 9, 5, 10, 5, 11, 5, 15, 5, 16, 5, 17,
   9, 3, 9, 4, 9, 5, 9, 9, 9, 10, 9, 11, 9, 15, 9, 16, 9, 17,
   10, 3, 10, 4, 10, 5, 10, 9, 10, 10, 10, 11, 10, 15, 10, 16, 10, 17,
   11, 3, 11, 4, 11, 5, 11, 9, 11, 10, 11, 11, 11, 15, 11, 16, 11, 17,
   15, 3, 15, 4, 15, 5, 15, 9, 15, 10, 15, 11, 15, 15, 15, 16, 15, 17,
   16, 3, 16, 4, 16, 5, 16, 9, 16, 10, 16, 11, 16, 15, 16, 16, 16, 17,
   17, 3, 17, 4, 17, 5, 17, 9, 17, 10, 17, 11, 17, 15, 17, 16, 17, 17};

int TACSLinearElasticity3D::getSpatialDim(){
  return 3;
}

int TACSLinearElasticity3D::getVarsPerNode(){
  return 3;
}

int TACSLinearElasticity3D::getDesignVarsPerNode(){
  return stiff->getDesignVarsPerNode();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSLinearElasticity3D::getDesignVarNums( int elemIndex, int dvLen,
                                              int dvNums[] ){
  return stiff->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSLinearElasticity3D::setDesignVars( int elemIndex, int dvLen,
                                           const TacsScalar dvs[] ){
  return stiff->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSLinearElasticity3D::getDesignVars( int elemIndex, int dvLen,
                                           TacsScalar dvs[] ){
  return stiff->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSLinearElasticity3D::getDesignVarRange( int elemIndex, int dvLen,
                                               TacsScalar lb[],
                                               TacsScalar ub[] ){
  return stiff->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSLinearElasticity3D::evalWeakIntegrand( int elemIndex,
                                                const double time,
                                                int n,
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
  DUt[8] = rho*Ut[8];

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
  stiff->evalStress(elemIndex, pt, X, e, s);

  DUx[0] = s[0];
  DUx[1] = s[5];
  DUx[2] = s[4];

  DUx[3] = s[5];
  DUx[4] = s[1];
  DUx[5] = s[3];

  DUx[6] = s[4];
  DUx[7] = s[3];
  DUx[8] = s[2];
}

void TACSLinearElasticity3D::evalWeakJacobian( int elemIndex,
                                               const double time,
                                               int n,
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

  DUt[6] = 0.0;
  DUt[7] = 0.0;
  DUt[8] = rho*Ut[8];

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
  stiff->evalStress(elemIndex, pt, X, e, s);

  DUx[0] = s[0];
  DUx[1] = s[5];
  DUx[2] = s[4];
    
  DUx[3] = s[5];
  DUx[4] = s[1];
  DUx[5] = s[3];
    
  DUx[6] = s[4];
  DUx[7] = s[3];
  DUx[8] = s[2];

  TacsScalar C[21];
  stiff->evalTangentStiffness(elemIndex, pt, X, C);

  // Set nonzero Jacobian terms
  *Jac_nnz = 84;
  *Jac_pairs = linear_Jac_pairs;

  // Acceleration terms
  Jac[0] = rho;
  Jac[1] = rho;
  Jac[2] = rho;

  if (strain_type == TACS_LINEAR_STRAIN){
    // u,x u,y u,z v,x v,y v,z w,x w,y w,z
    // e[0] = Ux[0]; e[1] = Ux[4]; e[2] = Ux[8];
    // e[3] = Ux[5] + Ux[7]; e[4] = Ux[2] + Ux[6]; e[5] = Ux[1] + Ux[3];

    // s =  [s0 s5 s4]
    //      [s5 s1 s3]
    //      [s4 s3 s2]

    // Index:
    // s[0] = C[0]*(u,x) + C[1]*(v,y) + C[2]*(w,z) + C[3]*(v,z + w,y)
    //                   + C[4]*(u,z + w,x) + C[5]*(u,y + v,x)
    // s[1] = C[1]*(u,x) + C[6]*(v,y) + C[7]*(w,z) + C[8]*(v,z + w,y)
    //                   + C[9]*(u,z + w,x) + C[10]*(u,y + v,x)
    // s[2] = C[2]*(u,x) + C[7]*(v,y) + C[11]*(w,z) + C[12]*(v,z + w,y)
    //                   + C[13]*(u,z + w,x) + C[14]*(u,y + v,x)
    // s[3] = C[3]*(u,x) + C[8]*(v,y) + C[12]*(w,z) + C[15]*(v,z + w,y)
    //                   + C[16]*(u,z + w,x) + C[17]*(u,y + v,x)
    // s[4] = C[4]*(u,x) + C[9]*(v,y) + C[13]*(w,z) + C[16]*(v,z + w,y)
    //                   + C[18]*(u,z + w,x) + C[19]*(u,y + v,x)
    // s[5] = C[5]*(u,x) + C[10]*(v,y) + C[14]*(w,z) + C[17]*(v,z + w,y)
    //                   + C[19]*(u,z + w,x) + C[20]*(u,y + v,x)

    // 0;   1;    2;   3;   4;   5;
    // u; u,t; u,tt; u,x; u,y; u,z;

    // 6;   7;    8;   9;  10;  11;
    // v; v,t; v,tt; v,x; v,y; v,z;

    //12;  13;   14;  15;  16;  17;
    // w; w,t; w,tt; w,x; w,y; w,z;
      
    // s[0]
    Jac[3] = C[0]; // u,x 3
    Jac[4] = C[5]; // u,y 4
    Jac[5] = C[4]; // u,z 5
    Jac[6] = C[5]; // v,x 9
    Jac[7] = C[1]; // v,y 10
    Jac[8] = C[3]; // v,z 11
    Jac[9] = C[4]; // w,x 15
    Jac[10] = C[3]; // w,y 16
    Jac[11] = C[2]; // w,z 17

    // s[5]
    Jac[12] = C[5];
    Jac[13] = C[20];
    Jac[14] = C[19];
    Jac[15] = C[20];
    Jac[16] = C[10];
    Jac[17] = C[17];
    Jac[18] = C[19];
    Jac[19] = C[17];
    Jac[20] = C[14];

    // s[4]
    Jac[21] = C[4];
    Jac[22] = C[19];
    Jac[23] = C[18];
    Jac[24] = C[19];
    Jac[25] = C[9];
    Jac[26] = C[16];
    Jac[27] = C[18];
    Jac[28] = C[16];
    Jac[29] = C[13];

    // s[5]
    Jac[30] = C[5];
    Jac[31] = C[20];
    Jac[32] = C[19];
    Jac[33] = C[20];
    Jac[34] = C[10];
    Jac[35] = C[17];
    Jac[36] = C[19];
    Jac[37] = C[17];
    Jac[38] = C[14];

    // s[1]
    Jac[39] = C[1];
    Jac[40] = C[10];
    Jac[41] = C[9];
    Jac[42] = C[10];
    Jac[43] = C[6];
    Jac[44] = C[8];
    Jac[45] = C[9];
    Jac[46] = C[8];
    Jac[47] = C[7];

    // s[3]
    Jac[48] = C[3];
    Jac[49] = C[17];
    Jac[50] = C[16];
    Jac[51] = C[17];
    Jac[52] = C[8];
    Jac[53] = C[15];
    Jac[54] = C[16];
    Jac[55] = C[15];
    Jac[56] = C[12];

    // s[4]
    Jac[57] = C[4];
    Jac[58] = C[19];
    Jac[59] = C[18];
    Jac[60] = C[19];
    Jac[61] = C[9];
    Jac[62] = C[16];
    Jac[63] = C[18];
    Jac[64] = C[16];
    Jac[65] = C[13];

    // s[3]
    Jac[66] = C[3];
    Jac[67] = C[17];
    Jac[68] = C[16];
    Jac[69] = C[17];
    Jac[70] = C[8];
    Jac[71] = C[15];
    Jac[72] = C[16];
    Jac[73] = C[15];
    Jac[74] = C[12];
    
    // s[2]
    Jac[75] = C[2];
    Jac[76] = C[14];
    Jac[77] = C[13];
    Jac[78] = C[14];
    Jac[79] = C[7];
    Jac[80] = C[12];
    Jac[81] = C[13];
    Jac[82] = C[12];
    Jac[83] = C[11];
  }
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSLinearElasticity3D::addWeakAdjProduct( int elemIndex,
                                                const double time,
                                                int n,
                                                const double pt[],
                                                const TacsScalar X[],
                                                const TacsScalar Ut[],
                                                const TacsScalar Ux[],
                                                const TacsScalar Psi[],
                                                const TacsScalar Psix[],
                                                TacsScalar scale,
                                                int dvLen,
                                                TacsScalar *fdvSens ){
  // Evaluate the density
  TacsScalar rho_coef = scale*(Ut[2]*Psi[0] + Ut[5]*Psi[1] + Ut[8]*Psi[2]);
  stiff->addDensityDVSens(elemIndex, pt, X, rho_coef, dvLen, fdvSens);

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

  TacsScalar phi[6];
  phi[0] = Psix[0];
  phi[1] = Psix[4];
  phi[2] = Psix[8];
  phi[3] = Psix[5] + Psix[7];
  phi[4] = Psix[2] + Psix[6];
  phi[5] = Psix[1] + Psix[3];

  stiff->addStressDVSens(elemIndex, pt, X, e, scale, phi, dvLen, fdvSens);
}

void TACSLinearElasticity3D::evalWeakAdjXptSensProduct( int elemIndex,
                                                        const double time,
                                                        int n,
                                                        const double pt[],
                                                        const TacsScalar X[],
                                                        const TacsScalar Ut[],
                                                        const TacsScalar Ux[],
                                                        const TacsScalar Psi[],
                                                        const TacsScalar Psix[],
                                                        TacsScalar *product,
                                                        TacsScalar dfdX[],
                                                        TacsScalar dfdUx[],
                                                        TacsScalar dfdPsix[] ){
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;
  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;
  dfdPsix[0] = dfdPsix[1] = dfdPsix[2] = 0.0;
  dfdPsix[3] = dfdPsix[4] = dfdPsix[5] = 0.0;
  dfdPsix[6] = dfdPsix[7] = dfdPsix[8] = 0.0;

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

  TacsScalar phi[6];
  phi[0] = Psix[0];
  phi[1] = Psix[4];
  phi[2] = Psix[8];
  phi[3] = Psix[5] + Psix[7];
  phi[4] = Psix[2] + Psix[6];
  phi[5] = Psix[1] + Psix[3];

  TacsScalar t1[6], t2[6];
  stiff->evalStress(elemIndex, pt, X, e, t1);
  stiff->evalStress(elemIndex, pt, X, phi, t2);

  *product = (t2[0]*e[0] + t2[1]*e[1] + t2[2]*e[2] +
              t2[3]*e[3] + t2[4]*e[4] + t2[5]*e[5]);

  if (strain_type == TACS_LINEAR_STRAIN){
    dfdUx[0] = t2[0];
    dfdUx[4] = t2[1];
    dfdUx[8] = t2[2];
    dfdUx[5] = t2[3];
    dfdUx[7] = t2[3];
    dfdUx[2] = t2[4];
    dfdUx[6] = t2[4];
    dfdUx[1] = t2[5];
    dfdUx[3] = t2[5];

    dfdPsix[0] = t1[0];
    dfdPsix[4] = t1[1];
    dfdPsix[8] = t1[2];
    dfdPsix[5] = t1[3];
    dfdPsix[7] = t1[3];
    dfdPsix[2] = t1[4];
    dfdPsix[6] = t1[4];
    dfdPsix[1] = t1[5];
    dfdPsix[3] = t1[5];
  }
}
/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearElasticity3D::evalPointQuantity( int elemIndex,
                                               const int quantityType,
                                               const double time,
                                               int n, const double pt[],
                                               const TacsScalar X[],
                                               const TacsScalar Xd[],
                                               const TacsScalar Ut[],
                                               const TacsScalar Ux[],
                                               TacsScalar *quantity ){
  if (quantityType == TACS_FAILURE_INDEX){
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

    *quantity = stiff->evalFailure(elemIndex, pt, X, e);

    return 1;
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    *quantity = stiff->evalDensity(elemIndex, pt, X);

    return 1;
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    // Evaluate the strain energy density
    *quantity = (e[0]*s[0] + e[1]*s[1] + e[2]*s[2] +
                 e[3]*s[3] + e[4]*s[4] + e[5]*s[5]);

    return 1;
  }

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSLinearElasticity3D::addPointQuantityDVSens( int elemIndex,
                                                     const int quantityType,
                                                     const double time,
                                                     TacsScalar scale,
                                                     int n, const double pt[],
                                                     const TacsScalar X[],
                                                     const TacsScalar Xd[],
                                                     const TacsScalar Ut[],
                                                     const TacsScalar Ux[],
                                                     const TacsScalar dfdq[],
                                                     int dvLen,
                                                     TacsScalar dfdx[] ){
  if (quantityType == TACS_FAILURE_INDEX){
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

    stiff->addFailureDVSens(elemIndex, pt, X, e, scale*dfdq[0],
                            dvLen, dfdx);
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    stiff->addDensityDVSens(elemIndex, pt, X, scale*dfdq[0], dvLen, dfdx);
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);
    stiff->addStressDVSens(elemIndex, pt, X, e, scale*dfdq[0],
                           e, dvLen, dfdx);
  }
}

/*
  Evaluate the derivatives of the point-wise quantity of interest
  with respect to X, Ut and Ux.
*/
void TACSLinearElasticity3D::evalPointQuantitySens( int elemIndex,
                                                    const int quantityType,
                                                    const double time,
                                                    int n, const double pt[],
                                                    const TacsScalar X[],
                                                    const TacsScalar Xd[],
                                                    const TacsScalar Ut[],
                                                    const TacsScalar Ux[],
                                                    const TacsScalar dfdq[],
                                                    TacsScalar dfdX[],
                                                    TacsScalar dfdXd[],
                                                    TacsScalar dfdUt[],
                                                    TacsScalar dfdUx[] ){
  dfdX[0] = dfdX[1] = dfdX[2] = 0.0;

  dfdXd[0] = dfdXd[1] = dfdXd[2] = 0.0;
  dfdXd[3] = dfdXd[4] = dfdXd[5] = 0.0;
  dfdXd[6] = dfdXd[7] = dfdXd[8] = 0.0;

  dfdUt[0] = dfdUt[1] = dfdUt[2] = 0.0;
  dfdUt[3] = dfdUt[4] = dfdUt[5] = 0.0;
  dfdUt[6] = dfdUt[7] = dfdUt[8] = 0.0;

  dfdUx[0] = dfdUx[1] = dfdUx[2] = 0.0;
  dfdUx[3] = dfdUx[4] = dfdUx[5] = 0.0;
  dfdUx[6] = dfdUx[7] = dfdUx[8] = 0.0;

  if (quantityType == TACS_FAILURE_INDEX){
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

    TacsScalar sens[6];
    stiff->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    if (strain_type == TACS_LINEAR_STRAIN){
      dfdUx[0] = dfdq[0]*sens[0];
      dfdUx[4] = dfdq[0]*sens[1];
      dfdUx[8] = dfdq[0]*sens[2];

      dfdUx[5] = dfdq[0]*sens[3];
      dfdUx[7] = dfdq[0]*sens[3];

      dfdUx[2] = dfdq[0]*sens[4];
      dfdUx[6] = dfdq[0]*sens[4];

      dfdUx[1] = dfdq[0]*sens[5];
      dfdUx[3] = dfdq[0]*sens[5];
    }
  }
  else if (quantityType == TACS_STRAIN_ENERGY_DENSITY){
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

    TacsScalar s[6];
    stiff->evalStress(elemIndex, pt, X, e, s);

    if (strain_type == TACS_LINEAR_STRAIN){
      dfdUx[0] = 2.0*dfdq[0]*s[0];
      dfdUx[4] = 2.0*dfdq[0]*s[1];
      dfdUx[8] = 2.0*dfdq[0]*s[2];

      dfdUx[5] = 2.0*dfdq[0]*s[3];
      dfdUx[7] = 2.0*dfdq[0]*s[3];

      dfdUx[2] = 2.0*dfdq[0]*s[4];
      dfdUx[6] = 2.0*dfdq[0]*s[4];

      dfdUx[1] = 2.0*dfdq[0]*s[5];
      dfdUx[3] = 2.0*dfdq[0]*s[5];
    }
  }
}

/*
  Get the data for visualization at a given point
*/
void TACSLinearElasticity3D::getOutputData( int elemIndex,
                                            const double time,
                                            ElementType etype,
                                            int write_flag,
                                            const double pt[],
                                            const TacsScalar X[],
                                            const TacsScalar Ut[],
                                            const TacsScalar Ux[],
                                            int ld_data,
                                            TacsScalar *data ){
  if (etype == TACS_SOLID_ELEMENT){
    if (write_flag & TACS_OUTPUT_NODES){
      // doesn't this depend whether it's linear/quadratic/etc?
      data[0] = X[0];
      data[1] = X[1];
      data[2] = X[2];
      data += 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
      data[0] = Ut[0];
      data[1] = Ut[3];
      data[3] = Ut[6];
      data += 3;
    }

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

    if (write_flag & TACS_OUTPUT_STRAINS){
      data[0] = e[0];
      data[1] = e[1];
      data[2] = e[2];
      data[3] = e[3];
      data[4] = e[4];
      data[5] = e[5];
      data += 6;
    }
    if (write_flag & TACS_OUTPUT_STRESSES){
      TacsScalar s[6];
      stiff->evalStress(elemIndex, pt, X, e, s);
      data[0] = s[0];
      data[1] = s[1];
      data[2] = s[2];
      data[3] = s[3];
      data[4] = s[4];
      data[5] = s[5];
      data += 6;
    }
    if (write_flag & TACS_OUTPUT_EXTRAS){
      data[0] = stiff->evalFailure(elemIndex, pt, X, e);
      data[1] = stiff->evalDesignFieldValue(elemIndex, pt, X, 0);
      data[2] = stiff->evalDesignFieldValue(elemIndex, pt, X, 1);
      data[3] = stiff->evalDesignFieldValue(elemIndex, pt, X, 2);
      data += 3;
    }
  }
}
