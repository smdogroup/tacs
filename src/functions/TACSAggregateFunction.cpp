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

#include "TACSAggregateFunction.h"

#include "TACSAssembler.h"

/*
  Initialize the TACSAggregateFunction class properties
*/
TACSAggregateFunction::TACSAggregateFunction(TACSAssembler *_assembler,
                                             int outputType, double _aggWeight,
                                             double _alpha,
                                             double _funcScaleFactor)
    : TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::TWO_STAGE, 0) {
  this->outputType = outputType;
  this->aggWeight = _aggWeight;
  this->alpha = _alpha;
  this->funcScaleFactor = _funcScaleFactor;
  this->aggType = CONTINUOUS;

  // Initialize the maximum failure value and KS sum to default values
  // that will be overwritten later.
  this->maxVal = -1e20;
  this->aggSum = 0.0;
  this->invPnorm = 0.0;
}

TACSAggregateFunction::~TACSAggregateFunction() {}

/*
  TACSAggregateFunction function name
*/
const char *TACSAggregateFunction::funcName = "TACSAggregateFunction";

/*
  Set the KS aggregation type
*/
void TACSAggregateFunction::setAggregationType(enum AggregationType type) {
  this->aggType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSAggregateFunction::getParameter() { return this->aggWeight; }

/*
  Set the KS aggregation parameter
*/
void TACSAggregateFunction::setParameter(double _aggWeight) {
  this->aggWeight = _aggWeight;
}

/*
  Return the function name
*/
const char *TACSAggregateFunction::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSAggregateFunction::getFunctionValue() {
  // Compute the final value of the aggregate function on all processors
  return this->maxVal + log(this->aggSum / this->alpha) / this->aggWeight;
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSAggregateFunction::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    this->maxVal = -1e20;
  } else if (ftype == TACSFunction::INTEGRATE) {
    this->aggSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSAggregateFunction::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = this->maxVal;
    MPI_Allreduce(&temp, &this->maxVal, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                  assembler->getMPIComm());
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = this->aggSum;
    MPI_Allreduce(&temp, &this->aggSum, 1, TACS_MPI_TYPE, MPI_SUM,
                  assembler->getMPIComm());

    // Compute the P-norm quantity if needed
    this->invPnorm = 0.0;
    if (this->aggType == PNORM_DISCRETE || this->aggType == PNORM_CONTINUOUS) {
      if (this->aggSum != 0.0) {
        this->invPnorm =
            pow(this->aggSum, (1.0 - this->aggWeight) / this->aggWeight);
      }
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSAggregateFunction function.
*/
void TACSAggregateFunction::elementWiseEval(
    EvaluationType ftype, int elemIndex, TACSElement *element, double time,
    TacsScalar scale, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  int numQuadPoints = element->getNumQuadraturePoints();
  for (int i = 0; i < numQuadPoints; i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate the failure index, and check whether it is an
    // undefined quantity of interest on this element
    TacsScalar funcVal = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, outputType, time, i, pt, Xpts,
                                   vars, dvars, ddvars, &detXd, &funcVal);

    // Scale failure value by safety factor
    funcVal *= this->funcScaleFactor;

    if (this->aggType == CONTINUOUS || this->aggType == PNORM_CONTINUOUS) {
      // Scale the scaling factor by the integration factor
      scale *= weight * detXd;
    }

    // Check whether the quantity requested is defined or not
    if (count >= 1) {
      this->addFuncVal(ftype, funcVal, scale);
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSAggregateFunction::getElementSVSens(
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

    TacsScalar funcVal = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, outputType, time, i, pt, Xpts,
                                   vars, dvars, ddvars, &detXd, &funcVal);

    // Scale failure value by safety factor
    funcVal *= this->funcScaleFactor;

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = this->computeFuncSens(funcVal);

      if (this->aggType == CONTINUOUS || this->aggType == PNORM_CONTINUOUS) {
        // Scale the scaling factor by the integration factor
        dfdq *= weight * detXd;
      }

      element->addPointQuantitySVSens(elemIndex, outputType, time, alpha, beta,
                                      gamma, i, pt, Xpts, vars, dvars, ddvars,
                                      &dfdq, dfdu);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSAggregateFunction::getElementXptSens(
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

    TacsScalar funcVal = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, outputType, time, i, pt, Xpts,
                                   vars, dvars, ddvars, &detXd, &funcVal);

    // Scale failure value by safety factor
    funcVal *= this->funcScaleFactor;

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      TacsScalar dfddetXd = 0.0;
      if (this->aggType == DISCRETE) {
        // d(log(this->aggSum))/dx = 1/(this->aggSum)*d(funcVal)/dx
        dfdq = exp(this->aggWeight * (funcVal - this->maxVal)) / this->aggSum;
      } else if (this->aggType == CONTINUOUS) {
        TacsScalar expfact =
            exp(this->aggWeight * (funcVal - this->maxVal)) / this->aggSum;
        dfddetXd = weight * expfact / this->aggWeight;
        dfdq = weight * detXd * expfact;
      } else if (this->aggType == PNORM_DISCRETE) {
        TacsScalar fpow = pow(fabs(TacsRealPart(funcVal / this->maxVal)),
                              this->aggWeight - 2.0);
        dfdq = funcVal * fpow * this->invPnorm;
      } else if (this->aggType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow = pow(fabs(TacsRealPart(funcVal / this->maxVal)),
                              this->aggWeight - 2.0);
        dfdq = funcVal * fpow * this->invPnorm * weight * detXd;
        dfddetXd = funcVal * fpow * this->invPnorm * weight;
      }

      // Scale failure sens by safety factor
      dfdq *= this->funcScaleFactor;

      element->addPointQuantityXptSens(elemIndex, outputType, time, scale, i,
                                       pt, Xpts, vars, dvars, ddvars, dfddetXd,
                                       &dfdq, dfdXpts);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSAggregateFunction::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar funcVal = 0.0, detXd = 0.0;
    int count =
        element->evalPointQuantity(elemIndex, outputType, time, i, pt, Xpts,
                                   vars, dvars, ddvars, &detXd, &funcVal);

    // Scale failure value by safety factor
    funcVal *= this->funcScaleFactor;

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = this->computeFuncSens(funcVal);

      if (this->aggType == CONTINUOUS || this->aggType == PNORM_CONTINUOUS) {
        // Scale the scaling factor by the integration factor
        dfdq *= weight * detXd;
      }

      element->addPointQuantityDVSens(elemIndex, outputType, time, scale, i, pt,
                                      Xpts, vars, dvars, ddvars, &dfdq, dvLen,
                                      dfdx);
    }
  }
}

void TACSAggregateFunction::addFuncVal(EvaluationType ftype, TacsScalar funcVal,
                                       TacsScalar scale) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Set the maximum value
    if (TacsRealPart(funcVal) > TacsRealPart(this->maxVal)) {
      this->maxVal = funcVal;
    }
  } else {
    // Add the function value to the sum
    TacsScalar f;
    if (this->aggType == DISCRETE) {
      f = exp(this->aggWeight * (funcVal - this->maxVal));
    } else if (this->aggType == CONTINUOUS) {
      f = exp(this->aggWeight * (funcVal - this->maxVal));
    } else if (this->aggType == PNORM_DISCRETE) {
      f = pow(fabs(TacsRealPart(funcVal / this->maxVal)), this->aggWeight);
    } else if (this->aggType == PNORM_CONTINUOUS) {
      f = pow(fabs(TacsRealPart(funcVal / this->maxVal)), this->aggWeight);
    }
    this->aggSum += scale * f;
  }
}

TacsScalar TACSAggregateFunction::computeFuncSens(TacsScalar funcVal) {
  // Compute the sensitivity contribution
  TacsScalar dfdq;
  if (this->aggType == DISCRETE || this->aggType == CONTINUOUS) {
    // d(log(this->aggSum))/dx = 1/(this->aggSum)*d(funcVal)/dx
    dfdq = exp(this->aggWeight * (funcVal - this->maxVal)) / this->aggSum;
  } else if (this->aggType == PNORM_DISCRETE ||
             this->aggType == PNORM_CONTINUOUS) {
    TacsScalar fpow =
        pow(fabs(TacsRealPart(funcVal / this->maxVal)), this->aggWeight - 2.0);
    dfdq = funcVal * fpow * this->invPnorm;
  }

  // Scale failure sens by safety factor
  dfdq *= this->funcScaleFactor;

  return dfdq;
}
