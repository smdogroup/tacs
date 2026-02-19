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
  Test the density design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveDensity(TACSConstitutive *con, int elemIndex,
                                const double pt[], const TacsScalar X[],
                                int ndvs, const TacsScalar *dvs,
                                double dh = 1e-7, int test_print_level = 2,
                                double test_fail_atol = 1e-5,
                                double test_fail_rtol = 1e-5);

/**
  Test the specific heat design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveSpecificHeat(TACSConstitutive *con, int elemIndex,
                                     const double pt[], const TacsScalar X[],
                                     int ndvs, const TacsScalar *dvs, double dh,
                                     int test_print_level,
                                     double test_fail_atol,
                                     double test_fail_rtol);

/**
  Test the heat flux design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveHeatFlux(TACSConstitutive *con, int elemIndex,
                                 const double pt[], const TacsScalar X[],
                                 int ndvs, const TacsScalar *dvs, double dh,
                                 int test_print_level, double test_fail_atol,
                                 double test_fail_rtol);

/**
  Test the stress design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveStress(TACSConstitutive *con, int elemIndex,
                               const double pt[], const TacsScalar X[],
                               int ndvs, const TacsScalar *dvs, double dh,
                               int test_print_level, double test_fail_atol,
                               double test_fail_rtol);

/**
  Test the stress design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveThermalStrain(TACSConstitutive *con, int elemIndex,
                                      const double pt[], const TacsScalar X[],
                                      int ndvs, const TacsScalar *dvs,
                                      double dh, int test_print_level,
                                      double test_fail_atol,
                                      double test_fail_rtol);

/**
  Test the failure design variable sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param ndvs The length of the design variable array
  @param dvs The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveFailure(TACSConstitutive *con, int elemIndex,
                                const double pt[], const TacsScalar X[],
                                int ndvs, const TacsScalar *dvs, double dh,
                                int test_print_level, double test_fail_atol,
                                double test_fail_rtol);

/**
  Test the failure strain sensitivity implementation.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param pt The Gauss point location
  @param X The spatial location
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutiveFailureStrainSens(TACSConstitutive *con, int elemIndex,
                                          const double pt[],
                                          const TacsScalar X[], double dh,
                                          int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol);

/**
  Test the implementation of the constitutive object.

  @param con The constitutive object to test
  @param elemIndex The element index for the constitutive test
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestConstitutive(TACSConstitutive *con, int elemIndex, double dh = 1e-7,
                         int test_print_level = 2, double test_fail_atol = 1e-5,
                         double test_fail_rtol = 1e-5);

#endif  // TACS_CONSTITUTIVE_VERIFICATION_H
