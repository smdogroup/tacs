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
    TacsScalar m, TacsScalar I11, TacsScalar I22, TacsScalar I33,
    TacsScalar I12, TacsScalar I13, TacsScalar I23) {
  M[0] = M[6] = M[11] = m;
  M[15] = I11;
  M[18] = I22;
  M[20] = I33;
  M[16] = I12;
  M[17] = I13;
  M[19] = I23;
}
