/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_POINT_MASS_CONSTITUTIVE_H
#define TACS_POINT_MASS_CONSTITUTIVE_H

#include "TACSGeneralMassConstitutive.h"

/**
  This is the base class for the traditional point mass constitutive objects
  with no translation-rotation coupling. Assumes 6 dofs.

*/
class TACSPointMassConstitutive : public TACSGeneralMassConstitutive {
 public:
  TACSPointMassConstitutive(
      TacsScalar _m, TacsScalar _I11 = 0.0, TacsScalar _I22 = 0.0,
      TacsScalar _I33 = 0.0, TacsScalar _I12 = 0.0, TacsScalar _I13 = 0.0,
      TacsScalar _I23 = 0.0, int _mNum = -1, TacsScalar _mlb = 0.0,
      TacsScalar _mub = 1e20, int _I11Num = -1, TacsScalar _I11lb = 0.0,
      TacsScalar _I11ub = 1e20, int _I22Num = -1, TacsScalar _I22lb = 0.0,
      TacsScalar _I22ub = 1e20, int _I33Num = -1, TacsScalar _I33lb = 0.0,
      TacsScalar _I33ub = 1e20, int _I12Num = -1, TacsScalar _I12lb = -1e20,
      TacsScalar _I12ub = 1e20, int _I13Num = -1, TacsScalar _I13lb = -1e20,
      TacsScalar _I13ub = 1e20, int _I23Num = -1, TacsScalar _I23lb = -1e20,
      TacsScalar _I23ub = 1e20);

  // Extra info about the constitutive class
  const char *getObjectName();

  // Retrieve the global design variable numbers
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  // Set the element design variable from the design vector
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  // Get the element design variables values
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  // Evaluate the mass matrix
  void evalMassMatrix(int elemIndex, const double pt[], const TacsScalar X[],
                      TacsScalar M[]);

  // Evaluate the mass matrix inner product sens
  void addMassMatrixDVSensInnerProduct(int elemIndex, TacsScalar scale,
                                       const double pt[], const TacsScalar X[],
                                       const TacsScalar psi[],
                                       const TacsScalar phi[], int dvLen,
                                       TacsScalar dfdx[]);

 private:
  static const char *name;

  TacsScalar m, I11, I22, I33, I12, I13, I23;
  int mNum, I11Num, I22Num, I33Num, I12Num, I13Num, I23Num;
  TacsScalar mlb, mub, I11lb, I11ub, I22lb, I22ub, I33lb, I33ub, I12lb, I12ub,
      I13lb, I13ub, I23lb, I23ub;
};

#endif  // TACS_POINT_MASS_CONSTITUTIVE_H
