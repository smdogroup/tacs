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

#ifndef TACS_THERMOELASTICITY_H
#define TACS_THERMOELASTICITY_H

#include "TACSElementModel.h"
#include "TACSLinearElasticity.h"
#include "TACSPlaneStressConstitutive.h"

class TACSLinearThermoelasticity2D : public TACSElementModel {
 public:
  TACSLinearThermoelasticity2D( TACSPlaneStressConstitutive *_con,
                                ElementStrainType strain_type );
  ~TACSLinearThermoelasticity2D();

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

#endif // TACS_THERMOELASTICITY_H
