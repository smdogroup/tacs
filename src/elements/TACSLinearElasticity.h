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
#include "TACSSolidConstitutive.h"

enum ElementStrainType { TACS_LINEAR_STRAIN,
                         TACS_NONLINEAR_STRAIN };

class TACSLinearEleasticity2D : public TACSElementModel {
 public:
  TACSLinearEleasticity2D( PlaneStressStiffness *_con,
                           ElementStrainType strain_type );
  ~TACSLinearEleasticity2D();
  
  int getSpatialDim();
  int getVarsPerNode();
  void evalWeakIntegrand( const double time,
                          const double pt[],
                          const TacsScalar X[],
                          const TacsScalar U[],
                          const TacsScalar Udot[],
                          const TacsScalar Uddot[],
                          const TacsScalar Ux[],
                          TacsScalar DUt[],
                          TacsScalar DUx[] );
  void evalIntegrandDeriv( const double time,
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
  PlaneStressStiffness *stiff;

  // Constant member data
  static const int DDUt_pairs[4];
  static const int DDUx_pairs[32];
};

/*

class TACSLinearEleasticity3D : public TACSElementModel {
 public:
  TACSLinearEleasticity3D( SolidStiffness *_con );
  ~TACSLinearEleasticity3D();
  
  int getSpatialDim();
  int getVarsPerNode();
  void evalWeakIntegrand( const double time,
                          const double pt[],
                          const TacsScalar X[],
                          const TacsScalar U[],
                          const TacsScalar Udot[],
                          const TacsScalar Uddot[],
                          const TacsScalar Ux[],
                          TacsScalar DUt[],
                          TacsScalar DUx[] );
  void evalIntegrandDeriv( const double time,
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

