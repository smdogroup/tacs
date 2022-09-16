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

#include "TACSKSFailure.h"

#include "TACSAssembler.h"

/*
  Initialize the TACSKSFailure class properties
*/
TACSKSFailure::TACSKSFailure(TACSAssembler *_assembler, double _ksWeight,
                             double _alpha, double _safetyFactor)
    : TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::TWO_STAGE, 0) {
  ksWeight = _ksWeight;
  alpha = _alpha;
  safetyFactor = _safetyFactor;
  ksType = CONTINUOUS;

  // Initialize the maximum failure value and KS sum to default values
  // that will be overwritten later.
  maxFail = -1e20;
  ksFailSum = 0.0;
  invPnorm = 0.0;
}

TACSKSFailure::~TACSKSFailure() {}

/*
  TACSKSFailure function name
*/
const char *TACSKSFailure::funcName = "TACSKSFailure";

/*
  Set the KS aggregation type
*/
void TACSKSFailure::setKSFailureType(enum KSFailureType type) { ksType = type; }

/*
  Retrieve the KS aggregation weight
*/
double TACSKSFailure::getParameter() { return ksWeight; }

/*
  Set the KS aggregation parameter
*/
void TACSKSFailure::setParameter(double _ksWeight) { ksWeight = _ksWeight; }

/*
  Return the function name
*/
const char *TACSKSFailure::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSKSFailure::getFunctionValue() {
  // Compute the final value of the KS function on all processors
  TacsScalar ksFail = maxFail + log(ksFailSum / alpha) / ksWeight;

  return ksFail;
}

/*
  Retrieve the maximum value
*/
TacsScalar TACSKSFailure::getMaximumFailure() { return maxFail; }

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSFailure::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    maxFail = -1e20;
  } else if (ftype == TACSFunction::INTEGRATE) {
    ksFailSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSFailure::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxFail;
    MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                  assembler->getMPIComm());
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksFailSum;
    MPI_Allreduce(&temp, &ksFailSum, 1, TACS_MPI_TYPE, MPI_SUM,
                  assembler->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS) {
      if (ksFailSum != 0.0) {
        invPnorm = pow(ksFailSum, (1.0 - ksWeight) / ksWeight);
      }
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSFailure function.
*/
void TACSKSFailure::elementWiseEval(EvaluationType ftype, int elemIndex,
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
      if (ftype == TACSFunction::INITIALIZE) {
        // Set the maximum failure load
        if (TacsRealPart(fail) > TacsRealPart(maxFail)) {
          maxFail = fail;
        }
      } else {
        // Add the failure load to the sum
        if (ksType == DISCRETE) {
          TacsScalar fexp = exp(ksWeight * (fail - maxFail));
          ksFailSum += scale * fexp;
        } else if (ksType == CONTINUOUS) {
          TacsScalar fexp = exp(ksWeight * (fail - maxFail));
          ksFailSum += scale * weight * detXd * fexp;
        } else if (ksType == PNORM_DISCRETE) {
          TacsScalar fpow = pow(fabs(TacsRealPart(fail / maxFail)), ksWeight);
          ksFailSum += scale * fpow;
        } else if (ksType == PNORM_CONTINUOUS) {
          TacsScalar fpow = pow(fabs(TacsRealPart(fail / maxFail)), ksWeight);
          ksFailSum += scale * weight * detXd * fpow;
        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSFailure::getElementSVSens(
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

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx
        dfdq = exp(ksWeight * (fail - maxFail)) / ksFailSum;
      } else if (ksType == CONTINUOUS) {
        dfdq = exp(ksWeight * (fail - maxFail)) / ksFailSum;
        dfdq *= weight * detXd;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm;
        dfdq *= weight * detXd;
      }

      // Scale failure sens by safety factor
      dfdq *= safetyFactor;

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
void TACSKSFailure::getElementXptSens(
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

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      TacsScalar dfddetXd = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx
        dfdq = exp(ksWeight * (fail - maxFail)) / ksFailSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact = exp(ksWeight * (fail - maxFail)) / ksFailSum;
        dfddetXd = weight * expfact / ksWeight;
        dfdq = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm * weight * detXd;
        dfddetXd = fail * fpow * invPnorm * weight;
      }

      // Scale failure sens by safety factor
      dfdq *= safetyFactor;

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
void TACSKSFailure::addElementDVSens(
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

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx
        dfdq = exp(ksWeight * (fail - maxFail)) / ksFailSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact = exp(ksWeight * (fail - maxFail)) / ksFailSum;
        dfdq = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(fail / maxFail)), ksWeight - 2.0);
        dfdq = fail * fpow * invPnorm * weight * detXd;
      }

      // Scale failure sens by safety factor
      dfdq *= safetyFactor;

      element->addPointQuantityDVSens(elemIndex, TACS_FAILURE_INDEX, time,
                                      scale, i, pt, Xpts, vars, dvars, ddvars,
                                      &dfdq, dvLen, dfdx);
    }
  }
}
