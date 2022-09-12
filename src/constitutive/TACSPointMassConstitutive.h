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
  TACSPointMassConstitutive(TacsScalar m, TacsScalar I11 = 0.0,
                            TacsScalar I22 = 0.0, TacsScalar I33 = 0.0,
                            TacsScalar I12 = 0.0, TacsScalar I13 = 0.0,
                            TacsScalar I23 = 0.0);

  // Extra info about the constitutive class
  const char *getObjectName();

 private:
  static const char *name;
};

#endif  // TACS_POINT_MASS_CONSTITUTIVE_H
