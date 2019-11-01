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

#ifndef TACS_CONSTITUTIVE_VERIFICATION_H
#define TACS_CONSTITUTIVE_VERIFICATION_H

#include "TACSConstitutive.h"
#include "TACSElementVerification.h"

/**
  Test the residual implementation against the Lagrangian equations
  of motion. This relies on the element implementation of the kinetic
  and potential energies.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutive( TACSConstitutive *con,
                          int elemIndex,
                          double dh=1e-7,
                          int test_print_level=2,
                          double test_fail_atol=1e-5,
                          double test_fail_rtol=1e-5 );

#endif // TACS_CONSTITUTIVE_VERIFICATION_H
