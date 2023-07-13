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

#ifndef TACS_ELEMENT_VERIFICATION_H
#define TACS_ELEMENT_VERIFICATION_H

#include "TACSElement.h"
#include "TACSElementBasis.h"
#include "TACSElementModel.h"

/**
  Assign variables randomly to an array. This is useful for
  testing various things.
 */
void TacsSeedRandomGenerator(int seed);
void TacsGenerateRandomArray(TacsReal *array, int size, TacsReal lower = -1.0,
                             TacsReal upper = 1.0);
void TacsGenerateRandomArray(TacsComplex *array, int size,
                             TacsComplex lower = -1.0, TacsComplex upper = 1.0);

/*
  Find the largest absolute value of the difference between an array of test
  values and an array of reference values
*/
double TacsGetMaxError(TacsScalar *testVals, TacsScalar *refVals, int size,
                       int *max_index);

/*
  Find the maximum relative error between an array of test
  values and an array of reference values
*/
double TacsGetMaxRelError(TacsScalar *testVals, TacsScalar *refVals, int size,
                          int *max_index);

/*
  Print out the values and the relative errors
*/
void TacsPrintErrorComponents(FILE *fp, const char *descript,
                              TacsScalar *testVals, TacsScalar *refVals,
                              int size);

/**
 * @brief Assert that the values in testVals are close to the values in refVals.
 * The closeness is computed using a hybrid of absolute and relative errors that
 * should avoid false positives when the values are close to zero. Basically, it
 * computes a relative error for "large" values and an absolute error for
 * "small" values. The criterion is identical to numpy's `assert_allclose`
 * function:
 *
 * Failure is raised if::
 *   abs(`test` - `ref`) > (`atol` + `rtol` * abs(`ref`))
 *
 * @param testVals Array of test values
 * @param refVals Array of reference values
 * @param size Number of values in each array
 * @param atol Absolute tolerance
 * @param rtol Relative tolerance
 * @return bool True if all values are close, false otherwise
 */
bool TacsAssertAllClose(TacsScalar *testVals, TacsScalar *refVals, int size,
                        double atol, double rtol);

/**
 * @brief Assert that the values in testVals are close to the values in refVals.
 * The closeness is computed using a hybrid of absolute and relative errors that
 * should avoid false positives when the values are close to zero. Basically, it
 * computes a relative error for "large" values and an absolute error for
 * "small" values. The criterion is identical to numpy's `assert_allclose`
 * function:
 *
 * Failure is raised if::
 *   abs(`test` - `ref`) > (`atol` + `rtol` * abs(`ref`))
 *
 * @param testVals Array of test values
 * @param refVals Array of reference values
 * @param size Number of values in each array
 * @param atol Absolute tolerance
 * @param rtol Relative tolerance
 * @return bool True if all values are close, false otherwise
 */
bool TacsAssertAllClose(TacsScalar *testVals, TacsScalar *refVals, int size,
                        double atol, double rtol);

/*
  Perturb the input variables in the forward sense
*/
void TacsForwardDiffPerturb(TacsScalar *out, int size, const TacsScalar *orig,
                            const TacsScalar *pert, double dh);

/*
  Perturb the variables in the backward sense
*/
void TacsBackwardDiffPerturb(TacsScalar *out, int size, const TacsScalar *orig,
                             const TacsScalar *pert, double dh);

/*
  Form the forward approximation
*/
void TacsFormDiffApproximate(TacsScalar *forward, const TacsScalar *backward,
                             int size, TacsScalar dh);

/**
  Test the residual implementation against the Lagrangian equations
  of motion. This relies on the element implementation of the kinetic
  and potential energies.

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementResidual(TACSElement *element, int elemIndex, double time,
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            double dh = 1e-7, int test_print_level = 2,
                            double test_fail_atol = 1e-5,
                            double test_fail_rtol = 1e-5);

/**
  Test the Jacobian matrix implementation against the residual.

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementJacobian(TACSElement *element, int elemIndex, double time,
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            int col = -1, double dh = 1e-7,
                            int test_print_level = 2,
                            double test_fail_atol = 1e-5,
                            double test_fail_rtol = 1e-5);

/**
  Test the Jacobian matrix implementation against the residual.

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementMatFreeJacobian(TACSElement *element, int elemIndex,
                                   double time, const TacsScalar Xpts[],
                                   const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[], int col = -1,
                                   double dh = 1e-7, int test_print_level = 2,
                                   double test_fail_atol = 1e-5,
                                   double test_fail_rtol = 1e-5);

/**
  Test the adjoint-residual product implementation

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dvLen The length of the design variable array
  @param x The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestAdjResProduct(TACSElement *element, int elemIndex, double time,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          int dvLen, const TacsScalar *x, double dh = 1e-7,
                          int test_print_level = 2,
                          double test_fail_atol = 1e-5,
                          double test_fail_rtol = 1e-5);

/**
  Test the adjoint-residual product implementation

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestAdjResXptProduct(TACSElement *element, int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             const TacsScalar dvars[],
                             const TacsScalar ddvars[], double dh = 1e-7,
                             int test_print_level = 2,
                             double test_fail_atol = 1e-5,
                             double test_fail_rtol = 1e-5);

/**
  Test the matrix design variable sensitivity implementation

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dvLen The length of the design variable array
  @param x The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementMatDVSens(TACSElement *element, ElementMatrixType elemType,
                             int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             int dvLen, const TacsScalar *x, double dh = 1e-7,
                             int test_print_level = 2,
                             double test_fail_atol = 1e-5,
                             double test_fail_rtol = 1e-5);

/**
 Test the matrix nodal coordinate sensitivity implementation

 @param element The element object
 @param time Simulation time
 @param Xpts The element nodal variables
 @param vars The element state variables
 @param dvars The time derivatives of the state variables
 @param ddvars The second time derivatives of the state variables
 @param dh The finite-difference step size
 @param test_print_level The output level
 @param test_fail_atol The test absolute tolerance
 @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementMatXptSens(TACSElement *element, ElementMatrixType elemType,
                              int elemIndex, double time,
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              double dh = 1e-7, int test_print_level = 2,
                              double test_fail_atol = 1e-5,
                              double test_fail_rtol = 1e-5);

/**
  Test the matrix state variable sensitivity implementation

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementMatSVSens(TACSElement *element, ElementMatrixType elemType,
                             int elemIndex, double time,
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             double dh = 1e-7, int test_print_level = 2,
                             double test_fail_atol = 1e-5,
                             double test_fail_rtol = 1e-5);

/**
  Test the quantity output design variable sensitivities

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dvLen The length of the design variable array
  @param x The design variable array
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementQuantityDVSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh = 1e-7, int test_print_level = 2,
    double test_fail_atol = 1e-5, double test_fail_rtol = 1e-5);

/**
  Test the quantity output state variable sensitivities

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementQuantitySVSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh = 1e-7, int test_print_level = 2,
    double test_fail_atol = 1e-5, double test_fail_rtol = 1e-5);

/**
  Test the quantity output element node sensitivities

  @param element The element object
  @param time Simulation time
  @param Xpts The element nodal variables
  @param vars The element state variables
  @param dvars The time derivatives of the state variables
  @param ddvars The second time derivatives of the state variables
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementQuantityXptSens(
    TACSElement *element, int elemIndex, int quantityType, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], double dh = 1e-7, int test_print_level = 2,
    double test_fail_atol = 1e-5, double test_fail_rtol = 1e-5);

/**
  Test if the basis function derivatives are implemented correct

  @param basis The TACSElementBasis to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementBasisFunctions(TACSElementBasis *basis, double dh,
                                  int test_print_level, double test_fail_atol,
                                  double test_fail_rtol);

/**
  Test if the basis face normals are implemented correct

  @param basis The TACSElementBasis to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementBasisFaceNormals(TACSElementBasis *basis, double dh,
                                    int test_print_level, double test_fail_atol,
                                    double test_fail_rtol);

/**
  Test if the basis jacobian transform is implemented correct

  @param basis The TACSElementBasis to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementBasisJacobianTransform(TACSElementBasis *basis, double dh,
                                          int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol);

/**
  Test if the basis function derivatives are implemented correct

  @param basis The TACSElementBasis to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementBasis(TACSElementBasis *basis, double dh = 1e-7,
                         int test_print_level = 2, double test_fail_atol = 1e-5,
                         double test_fail_rtol = 1e-5);

/**
  Test the element model Jacobian

  @param model The TACSElementModel instance to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementModelJacobian(TACSElementModel *model, int elemIndex,
                                 const double time, double dh,
                                 int test_print_level, double test_fail_atol,
                                 double test_fail_rtol);

/**
  Test the element model adjoint xpt sens product

  @param model The TACSElementModel instance to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementModelAdjXptSensProduct(TACSElementModel *model,
                                          int elemIndex, const double time,
                                          double dh, int test_print_level,
                                          double test_fail_atol,
                                          double test_fail_rtol);

/**
  Test the element model functions

  @param model The TACSElementModel instance to check
  @param dh The finite-difference step size
  @param test_print_level The output level
  @param test_fail_atol The test absolute tolerance
  @param test_fail_rtol The test relative tolerance
*/
int TacsTestElementModel(TACSElementModel *model, int elemIndex,
                         const double time, double dh = 1e-7,
                         int test_print_level = 2, double test_fail_atol = 1e-5,
                         double test_fail_rtol = 1e-5);

#endif  // TACS_ELEMENT_VERIFICATION_H
