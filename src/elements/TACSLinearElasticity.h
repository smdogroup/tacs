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

#ifndef TACS_LINEAR_ELASTICITY_H
#define TACS_LINEAR_ELASTICITY_H

#include "TACSElementModel.h"
#include "TACSPlaneStressConstitutive.h"
// #include "TACSSolidConstitutive.h"

enum ElementStrainType { TACS_LINEAR_STRAIN,
                         TACS_NONLINEAR_STRAIN };

class TACSLinearElasticity2D : public TACSElementModel {
 public:
  TACSLinearElasticity2D( TACSPlaneStressConstitutive *_con,
                          ElementStrainType strain_type );
  ~TACSLinearElasticity2D();

  int getSpatialDim();
  int getVarsPerNode();

  /**
    Retrieve the global design variable numbers associated with this element
  */
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  /**
    Set the element design variables from the design vector
  */
  void setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  /**
    Get the element design variables values
  */
  void getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  /**
    Get the lower and upper bounds for the design variable values
  */
  void getDesignVarRange( int elemIndex, int dvLen,
                          TacsScalar lb[], TacsScalar ub[] );

  /**
    Evaluate the coefficients of the weak form integrand
  */
  void evalWeakIntegrand( int elemIndex, int n, const double time,
                          const double pt[], const TacsScalar X[],
                          const TacsScalar Ut[], const TacsScalar Ux[],
                          TacsScalar DUt[], TacsScalar DUx[] );

  /**
    Evaluate the derivatives of the weak form coefficients
  */
  void evalWeakJacobian( int elemIndex, int n, const double time,
                          const double pt[], const TacsScalar X[],
                          const TacsScalar Ut[], const TacsScalar Ux[],
                          TacsScalar DUt[], TacsScalar DUx[],
                          int *Jac_nnz, const int *_Jac_pairs[],
                          TacsScalar Jac[] );

  void addWeakAdjProduct( int elemIndex, int n, const double time,
                          const double pt[], const TacsScalar X[],
                          const TacsScalar Ut[], const TacsScalar Ux[],
                          const TacsScalar Psi[], const TacsScalar Psix[],
                          TacsScalar scale, int dvLen, TacsScalar *fdvSens );

  /**
    Get the output for a single node in the mesh
  */
  void getOutputData( int elemIndex, const double time,
                      ElementType etype, int write_flag,
                      const double pt[], const TacsScalar X[],
                      const TacsScalar Ut[], const TacsScalar Ux[],
                      int ld_data, TacsScalar *data );

 private:
  ElementStrainType strain_type;
  TACSPlaneStressConstitutive *stiff;

  // Constant member data
  static const int linear_Jac_pairs[36];
};

/*

class TACSLinearElasticity3D : public TACSElementModel {
 public:
  TACSLinearElasticity3D( SolidStiffness *_con );
  ~TACSLinearElasticity3D();

  int getSpatialDim();
  int getVarsPerNode();
  void evalWeakIntegrand( int elemIndex,
                          const double time,
                          const double pt[],
                          const TacsScalar X[],
                          const TacsScalar U[],
                          const TacsScalar Udot[],
                          const TacsScalar Uddot[],
                          const TacsScalar Ux[],
                          TacsScalar DUt[],
                          TacsScalar DUx[] );
  void evalIntegrandDeriv( int elemIndex,
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
                           TacsScalar DDUx[] );

 private:
  ElementStrainType strain_type;
  SolidStiffness *stiff;
  static const int DDt_non_zero_pairs[6] = {2, 2, 5, 5, 8, 8};
};
*/

#endif // TACS_LINEAR_ELASTICITY_H
