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

#include "TACSKSTemperature.h"

#include "TACSAssembler.h"

/*
  Initialize the TACSKSTemperature class properties
*/
TACSKSTemperature::TACSKSTemperature(TACSAssembler *_assembler,
                                     double _ksWeight, double _alpha)
    : TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::TWO_STAGE, 0) {
  ksWeight = _ksWeight;
  alpha = _alpha;
  ksType = CONTINUOUS;

  // Initialize the maximum temperature value and KS sum to default values
  // that will be overwritten later.
  maxTemp = -1e20;
  ksTempSum = 0.0;
  invPnorm = 0.0;
}

TACSKSTemperature::~TACSKSTemperature() {}

/*
  TACSKSTemperature function name
*/
const char *TACSKSTemperature::funcName = "TACSKSTemperature";

/*
  Set the KS aggregation type
*/
void TACSKSTemperature::setKSTemperatureType(enum KSTemperatureType type) {
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSKSTemperature::getParameter() { return ksWeight; }

/*
  Set the KS aggregation parameter
*/
void TACSKSTemperature::setParameter(double _ksWeight) { ksWeight = _ksWeight; }

/*
  Return the function name
*/
const char *TACSKSTemperature::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSKSTemperature::getFunctionValue() {
  // Compute the final value of the KS function on all processors
  TacsScalar ksTemp = maxTemp + log(ksTempSum / alpha) / ksWeight;

  return ksTemp;
}

/*
  Retrieve the maximum value
*/
TacsScalar TACSKSTemperature::getMaximumTemperature() { return maxTemp; }

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSTemperature::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    maxTemp = -1e20;
  } else if (ftype == TACSFunction::INTEGRATE) {
    ksTempSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSTemperature::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxTemp;
    MPI_Allreduce(&temp, &maxTemp, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                  assembler->getMPIComm());
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksTempSum;
    MPI_Allreduce(&temp, &ksTempSum, 1, TACS_MPI_TYPE, MPI_SUM,
                  assembler->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS) {
      if (ksTempSum != 0.0) {
        invPnorm = pow(ksTempSum, (1.0 - ksWeight) / ksWeight);
      }
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSTemperature function.
*/
void TACSKSTemperature::elementWiseEval(
    EvaluationType ftype, int elemIndex, TACSElement *element, double time,
    TacsScalar scale, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate the temperature, and check whether it is an
    // undefined quantity of interest on this element
    TacsScalar temperature = 0.0, detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_TEMPERATURE, time, i,
                                           pt, Xpts, vars, dvars, ddvars,
                                           &detXd, &temperature);

    // Check whether the quantity requested is defined or not
    if (count >= 1) {
      if (ftype == TACSFunction::INITIALIZE) {
        // Set the maximum temperature
        if (TacsRealPart(temperature) > TacsRealPart(maxTemp)) {
          maxTemp = temperature;
        }
      } else {
        // Add the temperature to the sum
        if (ksType == DISCRETE) {
          TacsScalar fexp = exp(ksWeight * (temperature - maxTemp));
          ksTempSum += scale * fexp;
        } else if (ksType == CONTINUOUS) {
          TacsScalar fexp = exp(ksWeight * (temperature - maxTemp));
          ksTempSum += scale * weight * detXd * fexp;
        } else if (ksType == PNORM_DISCRETE) {
          TacsScalar fpow =
              pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight);
          ksTempSum += scale * fpow;
        } else if (ksType == PNORM_CONTINUOUS) {
          TacsScalar fpow =
              pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight);
          ksTempSum += scale * weight * detXd * fpow;
        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSTemperature::getElementSVSens(
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

    TacsScalar temperature = 0.0, detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_TEMPERATURE, time, i,
                                           pt, Xpts, vars, dvars, ddvars,
                                           &detXd, &temperature);

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar ksPtWeight = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksTempSum))/dx = 1/(ksTempSum)*d(temperature)/dx
        ksPtWeight = exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
      } else if (ksType == CONTINUOUS) {
        ksPtWeight = exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
        ksPtWeight *= weight * detXd;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        ksPtWeight = temperature * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        ksPtWeight = temperature * fpow * invPnorm;
        ksPtWeight *= weight * detXd;
      }

      TacsScalar dfdq = ksPtWeight;
      element->addPointQuantitySVSens(elemIndex, TACS_TEMPERATURE, time, alpha,
                                      beta, gamma, i, pt, Xpts, vars, dvars,
                                      ddvars, &dfdq, dfdu);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSKSTemperature::getElementXptSens(
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

    TacsScalar temperature = 0.0, detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_TEMPERATURE, time, i,
                                           pt, Xpts, vars, dvars, ddvars,
                                           &detXd, &temperature);

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      TacsScalar dfddetXd = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksTempSum))/dx = 1/(ksTempSum)*d(temperature)/dx
        dfdq = exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact =
            exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
        dfddetXd = weight * expfact / ksWeight;
        dfdq = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        dfdq = temperature * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        dfdq = temperature * fpow * invPnorm * weight * detXd;
        dfddetXd = temperature * fpow * invPnorm * weight;
      }

      element->addPointQuantityXptSens(elemIndex, TACS_TEMPERATURE, time, scale,
                                       i, pt, Xpts, vars, dvars, ddvars,
                                       dfddetXd, &dfdq, dfdXpts);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSTemperature::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar temperature = 0.0, detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_TEMPERATURE, time, i,
                                           pt, Xpts, vars, dvars, ddvars,
                                           &detXd, &temperature);

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar dfdq = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksTempSum))/dx = 1/(ksTempSum)*d(temperature)/dx
        dfdq = exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact =
            exp(ksWeight * (temperature - maxTemp)) / ksTempSum;
        dfdq = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        dfdq = temperature * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(temperature / maxTemp)), ksWeight - 2.0);
        dfdq = temperature * fpow * invPnorm * weight * detXd;
      }

      element->addPointQuantityDVSens(elemIndex, TACS_TEMPERATURE, time, scale,
                                      i, pt, Xpts, vars, dvars, ddvars, &dfdq,
                                      dvLen, dfdx);
    }
  }
}
