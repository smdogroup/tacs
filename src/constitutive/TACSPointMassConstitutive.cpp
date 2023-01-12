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

#include "TACSPointMassConstitutive.h"

const char* TACSPointMassConstitutive::name = "TACSPointMassConstitutive";

const char* TACSPointMassConstitutive::getObjectName() { return name; }

/*
  PointMassConstitutive member function definitions
*/
TACSPointMassConstitutive::TACSPointMassConstitutive(
    TacsScalar _m, TacsScalar _I11, TacsScalar _I22, TacsScalar _I33,
    TacsScalar _I12, TacsScalar _I13, TacsScalar _I23, int _mNum,
    TacsScalar _mlb, TacsScalar _mub, int _I11Num, TacsScalar _I11lb,
    TacsScalar _I11ub, int _I22Num, TacsScalar _I22lb, TacsScalar _I22ub,
    int _I33Num, TacsScalar _I33lb, TacsScalar _I33ub, int _I12Num,
    TacsScalar _I12lb, TacsScalar _I12ub, int _I13Num, TacsScalar _I13lb,
    TacsScalar _I13ub, int _I23Num, TacsScalar _I23lb, TacsScalar _I23ub) {
  m = _m;
  I11 = _I11;
  I22 = _I22;
  I33 = _I33;
  I12 = _I12;
  I13 = _I13;
  I23 = _I23;

  mNum = _mNum;
  I11Num = _I11Num;
  I22Num = _I22Num;
  I33Num = _I33Num;
  I12Num = _I12Num;
  I13Num = _I13Num;
  I23Num = _I23Num;

  mlb = _mlb;
  mub = _mub;
  I11lb = _I11lb;
  I11ub = _I11ub;
  I22lb = _I22lb;
  I22ub = _I22ub;
  I33lb = _I33lb;
  I33ub = _I33ub;
  I12lb = _I12lb;
  I12ub = _I12ub;
  I13lb = _I13lb;
  I13ub = _I13ub;
  I23lb = _I23lb;
  I23ub = _I23ub;
}

// Retrieve the global design variable numbers
int TACSPointMassConstitutive::getDesignVarNums(int elemIndex, int dvLen,
                                                int dvNums[]) {
  int index = 0;
  if (mNum >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = mNum;
    }
    index++;
  }
  if (I11Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I11Num;
    }
    index++;
  }
  if (I22Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I22Num;
    }
    index++;
  }
  if (I33Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I33Num;
    }
    index++;
  }
  if (I12Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I12Num;
    }
    index++;
  }
  if (I13Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I13Num;
    }
    index++;
  }
  if (I23Num >= 0) {
    if (dvNums && dvLen > index) {
      dvNums[index] = I23Num;
    }
    index++;
  }
  return index;
}

// Set the element design variable from the design vector
int TACSPointMassConstitutive::setDesignVars(int elemIndex, int dvLen,
                                             const TacsScalar dvs[]) {
  int index = 0;
  if (mNum >= 0) {
    m = dvs[index];
    index++;
  }
  if (I11Num >= 0) {
    I11 = dvs[index];
    index++;
  }
  if (I22Num >= 0) {
    I22 = dvs[index];
    index++;
  }
  if (I33Num >= 0) {
    I33 = dvs[index];
    index++;
  }
  if (I12Num >= 0) {
    I12 = dvs[index];
    index++;
  }
  if (I13Num >= 0) {
    I13 = dvs[index];
    index++;
  }
  if (I23Num >= 0) {
    I23 = dvs[index];
    index++;
  }
  return index;
}

// Get the element design variables values
int TACSPointMassConstitutive::getDesignVars(int elemIndex, int dvLen,
                                             TacsScalar dvs[]) {
  int index = 0;
  if (mNum >= 0) {
    dvs[index] = m;
    index++;
  }
  if (I11Num >= 0) {
    dvs[index] = I11;
    index++;
  }
  if (I22Num >= 0) {
    dvs[index] = I22;
    index++;
  }
  if (I33Num >= 0) {
    dvs[index] = I33;
    index++;
  }
  if (I12Num >= 0) {
    dvs[index] = I12;
    index++;
  }
  if (I13Num >= 0) {
    dvs[index] = I13;
    index++;
  }
  if (I23Num >= 0) {
    dvs[index] = I23;
    index++;
  }
  return index;
}

// Get the lower and upper bounds for the design variable values
int TACSPointMassConstitutive::getDesignVarRange(int elemIndex, int dvLen,
                                                 TacsScalar lb[],
                                                 TacsScalar ub[]) {
  int index = 0;
  if (mNum >= 0) {
    lb[index] = mlb;
    ub[index] = mub;
    index++;
  }
  if (I11Num >= 0) {
    lb[index] = I11lb;
    ub[index] = I11ub;
    index++;
  }
  if (I22Num >= 0) {
    lb[index] = I22lb;
    ub[index] = I22ub;
    index++;
  }
  if (I33Num >= 0) {
    lb[index] = I33lb;
    ub[index] = I33ub;
    index++;
  }
  if (I12Num >= 0) {
    lb[index] = I12lb;
    ub[index] = I12ub;
    index++;
  }
  if (I13Num >= 0) {
    lb[index] = I13lb;
    ub[index] = I13ub;
    index++;
  }
  if (I23Num >= 0) {
    lb[index] = I23lb;
    ub[index] = I23ub;
    index++;
  }
  return index;
}

// Evaluate the mass matrix
void TACSPointMassConstitutive::evalMassMatrix(int elemIndex, const double pt[],
                                               const TacsScalar X[],
                                               TacsScalar M[]) {
  memset(M, 0, 21 * sizeof(TacsScalar));
  M[0] = M[6] = M[11] = m;
  M[15] = I11;
  M[18] = I22;
  M[20] = I33;
  M[16] = I12;
  M[17] = I13;
  M[19] = I23;
}

// Evaluate the mass matrix inner product sens
void TACSPointMassConstitutive::addMassMatrixDVSensInnerProduct(
    int elemIndex, TacsScalar scale, const double pt[], const TacsScalar X[],
    const TacsScalar psi[], const TacsScalar phi[], int dvLen,
    TacsScalar dfdx[]) {
  int index = 0;
  if (mNum >= 0) {
    dfdx[index] +=
        scale * (psi[0] * phi[0] + psi[1] * phi[1] + psi[2] * phi[2]);
    index++;
  }

  if (I11Num >= 0) {
    dfdx[index] += scale * (psi[3] * phi[3]);
    index++;
  }

  if (I22Num >= 0) {
    dfdx[index] += scale * (psi[4] * phi[4]);
    index++;
  }

  if (I33Num >= 0) {
    dfdx[index] += scale * (psi[5] * phi[5]);
    index++;
  }

  if (I12Num >= 0) {
    dfdx[index] += scale * (psi[3] * phi[4] + psi[4] * phi[3]);
    index++;
  }

  if (I13Num >= 0) {
    dfdx[index] += scale * (psi[3] * phi[5] + psi[5] * phi[3]);
    index++;
  }

  if (I23Num >= 0) {
    dfdx[index] += scale * (psi[4] * phi[5] + psi[5] * phi[4]);
    index++;
  }
}
