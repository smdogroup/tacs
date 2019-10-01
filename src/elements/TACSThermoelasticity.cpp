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
  {2, 2, 7, 7, 11, 11,
   3, 3, 3, 4, 3, 8, 3, 9, 3, 10,
   4, 3, 4, 4, 4, 8, 4, 9, 4, 10,
   8, 3, 8, 4, 8, 8, 8, 9, 8, 10,
   9, 3, 9, 4, 9, 8, 9, 9, 9, 10,
   13, 13, 13, 14,
   14, 13, 14, 14};

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
                                                      const double time,
                                                      int n,
                                                      const double pt[],
                                                      const TacsScalar X[],
                                                      const TacsScalar Ut[],
                                                      const TacsScalar Ux[],
                                                      TacsScalar DUt[],
                                                      TacsScalar DUx[] ){
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Ut[2];

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Ut[5];

  DUt[6] = 0.0;
  DUt[7] = c*rho*Ut[7];
  DUt[8] = 0.0;

  // Compute the thermal strain components
  TacsScalar theta = Ut[6]; // The temperature value
  TacsScalar et[3];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0] - et[0];
    e[1] = Ux[3] - et[1];
    e[2] = Ux[1] + Ux[2] - et[2];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[2]*Ux[2]) - et[0];
    e[1] = Ux[3] + 0.5*(Ux[1]*Ux[1] + Ux[3]*Ux[3]) - et[1];
    e[2] = Ux[1] + Ux[2] + (Ux[0]*Ux[1] + Ux[2]*Ux[3]) - et[2];
  }

  // Evaluate the components of the stress
  TacsScalar s[3];
  stiff->evalStress(elemIndex, pt, X, e, s);
  DUx[0] = s[0];
  DUx[1] = s[2];

  DUx[2] = s[2];
  DUx[3] = s[1];

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2], flux[2];
  grad[0] = Ux[4];
  grad[1] = Ux[5];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[4] = flux[0];
  DUx[5] = flux[1];
}

void TACSLinearThermoelasticity2D::evalWeakJacobian( int elemIndex,
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
  // Evaluate the density and specific heat
  TacsScalar rho = stiff->evalDensity(elemIndex, pt, X);
  TacsScalar c = stiff->evalSpecificHeat(elemIndex, pt, X);

  DUt[0] = 0.0;
  DUt[1] = 0.0;
  DUt[2] = rho*Ut[2]; // u,tt

  DUt[3] = 0.0;
  DUt[4] = 0.0;
  DUt[5] = rho*Ut[5]; // v,tt

  DUt[6] = 0.0;
  DUt[7] = c*rho*Ut[7];
  DUt[8] = 0.0;

  // Compute the thermal strain components
  TacsScalar theta = Ut[6]; // The temperature value
  TacsScalar et[3];
  stiff->evalThermalStrain(elemIndex, pt, X, theta, et);

  // Compute the mechanical strain e = 0.5*(u,x + u,x^{T}) - et
  TacsScalar e[3];
  if (strain_type == TACS_LINEAR_STRAIN){
    e[0] = Ux[0] - et[0];
    e[1] = Ux[3] - et[1];
    e[2] = Ux[1] + Ux[2] - et[2];
  }
  else {
    e[0] = Ux[0] + 0.5*(Ux[0]*Ux[0] + Ux[2]*Ux[2]) - et[0];
    e[1] = Ux[3] + 0.5*(Ux[1]*Ux[1] + Ux[3]*Ux[3]) - et[1];
    e[2] = Ux[1] + Ux[2] + (Ux[0]*Ux[1] + Ux[2]*Ux[3]) - et[2];
  }

  // Evaluate the components of the stress
  TacsScalar s[3];
  stiff->evalStress(elemIndex, pt, X, e, s);
  DUx[0] = s[0];
  DUx[1] = s[2];

  DUx[2] = s[2];
  DUx[3] = s[1];

  // Compute the thermal flux from the thermal gradient
  TacsScalar grad[2], flux[2];
  grad[0] = Ux[4];
  grad[1] = Ux[5];

  // Add the flux components to the heat transfer portion
  // of the governing equations
  stiff->evalHeatFlux(elemIndex, pt, X, grad, flux);
  DUx[4] = flux[0];
  DUx[5] = flux[1];

  // Set the non-zero terms in the Jacobian
  *Jac_nnz = 27;
  *Jac_pairs = linear_Jac_pairs;

  // Set the time-dependent terms
  Jac[0] = rho;
  Jac[1] = rho;
  Jac[2] = c*rho;

  // Compute the unit strain
  TacsScalar C[6], Kc[3];
  stiff->evalThermalStrain(elemIndex, pt, X, 1.0, et);
  stiff->evalTangentStiffness(elemIndex, pt, X, C);
  stiff->evalTangentHeatFlux(elemIndex, pt, X, Kc);

  // Add the terms for linear thermoelasticity
  if (strain_type == TACS_LINEAR_STRAIN){
    // s[0] = C[0]*(u,x - theta*et[0]) + C[1]*(v,y - theta*et[1])
    //      + C[2]*(u,y + v,x - theta*et[2])
    // s[1] = C[1]*(u,x) + C[3]*(v,y) + C[4]*(u,y + v,x)
    // s[2] = C[2]*(u,x - theta*et[0]) + C[4]*(v,y - theta*et[1])
    //      + C[5]*(u,y + v,x - theta*et[2])

    // i == 3 (s[0])
    Jac[3] = C[0]; // j == 3
    Jac[4] = C[2]; // j == 4
    Jac[5] = C[2]; // j == 8
    Jac[6] = C[1]; // j == 9
    Jac[7] = -(C[0]*et[0] + C[1]*et[1] + C[2]*et[2]); // j == 10

    // i == 4 (s[2])
    Jac[8] = C[2]; // j == 3
    Jac[9] = C[5]; // j == 4
    Jac[10] = C[5]; // j == 8
    Jac[11] = C[4]; // j == 9
    Jac[12] = -(C[2]*et[0] + C[4]*et[1] + C[5]*et[2]); // j == 10

    // i == 8 (s[2])
    Jac[13] = C[2]; // j == 3
    Jac[14] = C[5]; // j == 4
    Jac[15] = C[5]; // j == 8
    Jac[16] = C[4]; // j == 9
    Jac[17] = -(C[2]*et[0] + C[4]*et[1] + C[5]*et[2]); // j == 10

    // i == 9 (s[1])
    Jac[18] = C[1]; // j == 3
    Jac[19] = C[4]; // j == 4
    Jac[20] = C[4]; // j == 8
    Jac[21] = C[3]; // j == 9
    Jac[22] = -(C[1]*et[0] + C[3]*et[1] + C[4]*et[2]); // j == 10
  }

  // i == 13
  Jac[23] = Kc[0]; // j == 13
  Jac[24] = Kc[1]; // j == 14

  // i == 14
  Jac[25] = Kc[1]; // j == 13
  Jac[26] = Kc[2]; // j == 14
}

/*
  Add the product of the adjoint vector times the weak form of the adjoint
  equations to the design variable components
*/
void TACSLinearThermoelasticity2D::addWeakAdjProduct( int elemIndex,
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

/*
  Evaluate a specified pointwise quantity of interest
*/
int TACSLinearThermoelasticity2D::evalPointQuantity( int elemIndex,
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

  return 0;
}

/*
  Add the derivative of the point-wise quantity of interest w.r.t.
  design variables to the design vector
*/
void TACSLinearThermoelasticity2D::addPointQuantityDVSens( int elemIndex,
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
}

/*
   Evaluate the derivatives of the point-wise quantity of interest
   with respect to X, Ut and Ux.
*/
void TACSLinearThermoelasticity2D::evalPointQuantitySens( int elemIndex,
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
      data[0] = stiff->evalFailure(elemIndex, pt, X, e);
      data[1] = 0.0;
      data[2] = 0.0;
      data += 3;
    }
  }
}
