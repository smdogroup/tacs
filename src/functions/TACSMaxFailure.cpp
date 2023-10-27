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

#include "TACSMaxFailure.h"

#include "TACSAssembler.h"

/*
  Initialize the TACSMaxFailure class properties
*/
TACSMaxFailure::TACSMaxFailure(TACSAssembler *_assembler, double _alpha,
                               double _safetyFactor)
    : TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::SINGLE_STAGE, 0) {
  alpha = _alpha;
  safetyFactor = _safetyFactor;

  // Initialize the maximum failure value and KS sum to default values
  // that will be overwritten later.
  maxFail = -1e20;
}

TACSMaxFailure::~TACSMaxFailure() {}

/*
  TACSMaxFailure function name
*/
const char *TACSMaxFailure::funcName = "TACSMaxFailure";

/*
  Return the function name
*/
const char *TACSMaxFailure::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSMaxFailure::getFunctionValue() { return maxFail; }

/*
  Initialize the internal values stored within the KS function
*/
void TACSMaxFailure::initEvaluation(EvaluationType ftype) { maxFail = -1e20; }

/*
  Reduce the function values across all MPI processes
*/
void TACSMaxFailure::finalEvaluation(EvaluationType ftype) {
  // Distribute the values of the KS function computed on this domain
  TacsScalar temp = maxFail;
  MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                assembler->getMPIComm());
}

/*
  Perform the element-wise evaluation of the TACSMaxFailure function.
*/
void TACSMaxFailure::elementWiseEval(EvaluationType ftype, int elemIndex,
                                     TACSElement *element, double time,
                                     TacsScalar scale, const TacsScalar Xpts[],
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[],
                                     const TacsScalar ddvars[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate the failure index, and check whether it is an
    // undefined quantity of interest on this element
    TacsScalar fail = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &fail);

    // Scale failure value by safety factor
    fail *= safetyFactor;

    // Check whether the quantity requested is defined or not
    if (count >= 1) {
      // Set the maximum failure load
      if (TacsRealPart(fail) > TacsRealPart(maxFail)) {
        maxFail = fail;
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSMaxFailure::getElementSVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar alpha,
    TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdu[]) {
  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->getNumVariables();
  memset(dfdu, 0, numVars * sizeof(TacsScalar));

  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar fail = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &fail);

    // Scale failure value by safety factor
    fail *= safetyFactor;

    if (count >= 1 && fail == maxFail) {
      TacsScalar dfdq = safetyFactor;
      element->addPointQuantitySVSens(elemIndex, TACS_FAILURE_INDEX, time,
                                      alpha, beta, gamma, i, pt, Xpts, vars,
                                      dvars, ddvars, &dfdq, dfdu);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSMaxFailure::getElementXptSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdXpts[]) {
  // Zero the derivative of the function w.r.t. the element node
  // locations
  int numNodes = element->getNumNodes();
  memset(dfdXpts, 0, 3 * numNodes * sizeof(TacsScalar));

  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar fail = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &fail);

    // Scale failure value by safety factor
    fail *= safetyFactor;

    if (count >= 1 && fail == maxFail) {
      TacsScalar dfdq = safetyFactor;
      TacsScalar dfddetXd = 0.0;
      element->addPointQuantityXptSens(elemIndex, TACS_FAILURE_INDEX, time,
                                       scale, i, pt, Xpts, vars, dvars, ddvars,
                                       dfddetXd, &dfdq, dfdXpts);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSMaxFailure::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar fail = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, TACS_FAILURE_INDEX, time, i, pt,
                                   Xpts, vars, dvars, ddvars, &detXd, &fail);

    // Scale failure value by safety factor
    fail *= safetyFactor;

    if (count >= 1 && fail == maxFail) {
      TacsScalar dfdq = safetyFactor;
      element->addPointQuantityDVSens(elemIndex, TACS_FAILURE_INDEX, time,
                                      scale, i, pt, Xpts, vars, dvars, ddvars,
                                      &dfdq, dvLen, dfdx);
    }
  }
}
