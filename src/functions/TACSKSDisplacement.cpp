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

#include "TACSKSDisplacement.h"

#include "TACSAssembler.h"

/*
  Initialize the TACSKSDisplacement class properties
*/
TACSKSDisplacement::TACSKSDisplacement(TACSAssembler *_assembler,
                                       double _ksWeight, const double _dir[],
                                       double _alpha)
    : TACSFunction(_assembler, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::TWO_STAGE, 0) {
  ksWeight = _ksWeight;
  dir[0] = _dir[0];
  dir[1] = _dir[1];
  dir[2] = _dir[2];
  alpha = _alpha;
  ksType = CONTINUOUS;

  // Initialize the maximum displacement value and KS sum to default values
  // that will be overwritten later.
  maxDisp = -1e20;
  ksDispSum = 0.0;
  invPnorm = 0.0;
}

TACSKSDisplacement::~TACSKSDisplacement() {}

/*
  TACSKSDisplacement function name
*/
const char *TACSKSDisplacement::funcName = "TACSKSDisplacement";

/*
  Set the KS aggregation type
*/
void TACSKSDisplacement::setKSDisplacementType(enum KSDisplacementType type) {
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSKSDisplacement::getParameter() { return ksWeight; }

/*
  Set the KS aggregation parameter
*/
void TACSKSDisplacement::setParameter(double _ksWeight) {
  ksWeight = _ksWeight;
}

/*
  Return the function name
*/
const char *TACSKSDisplacement::getObjectName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSKSDisplacement::getFunctionValue() {
  // Compute the final value of the KS function on all processors
  TacsScalar ksDisp = maxDisp + log(ksDispSum / alpha) / ksWeight;

  return ksDisp;
}

/*
  Retrieve the maximum value
*/
TacsScalar TACSKSDisplacement::getMaximumDisplacement() { return maxDisp; }

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSDisplacement::initEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    maxDisp = -1e20;
  } else if (ftype == TACSFunction::INTEGRATE) {
    ksDispSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSDisplacement::finalEvaluation(EvaluationType ftype) {
  if (ftype == TACSFunction::INITIALIZE) {
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxDisp;
    MPI_Allreduce(&temp, &maxDisp, 1, TACS_MPI_TYPE, TACS_MPI_MAX,
                  assembler->getMPIComm());
  } else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksDispSum;
    MPI_Allreduce(&temp, &ksDispSum, 1, TACS_MPI_TYPE, MPI_SUM,
                  assembler->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS) {
      if (ksDispSum != 0.0) {
        invPnorm = pow(ksDispSum, (1.0 - ksWeight) / ksWeight);
      }
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSDisplacement function.
*/
void TACSKSDisplacement::elementWiseEval(
    EvaluationType ftype, int elemIndex, TACSElement *element, double time,
    TacsScalar scale, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    // Evaluate the element displacements,
    // and project it against the user-defined vector
    TacsScalar dispVec[3], detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                           time, i, pt, Xpts, vars, dvars,
                                           ddvars, &detXd, dispVec);
    // project displacement along user-defined direction
    TacsScalar dispProj = 0.0;
    for (int i = 0; i < count; i++) {
      dispProj += dispVec[i] * dir[i];
    }

    // Check whether the quantity requested is defined or not
    if (count >= 1) {
      if (ftype == TACSFunction::INITIALIZE) {
        // Set the maximum displacement
        if (TacsRealPart(dispProj) > TacsRealPart(maxDisp)) {
          maxDisp = dispProj;
        }
      } else {
        // Add the displacement to the sum
        if (ksType == DISCRETE) {
          TacsScalar fexp = exp(ksWeight * (dispProj - maxDisp));
          ksDispSum += scale * fexp;
        } else if (ksType == CONTINUOUS) {
          TacsScalar fexp = exp(ksWeight * (dispProj - maxDisp));
          ksDispSum += scale * weight * detXd * fexp;
        } else if (ksType == PNORM_DISCRETE) {
          TacsScalar fpow =
              pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight);
          ksDispSum += scale * fpow;
        } else if (ksType == PNORM_CONTINUOUS) {
          TacsScalar fpow =
              pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight);
          ksDispSum += scale * weight * detXd * fpow;
        }
      }
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSDisplacement::getElementSVSens(
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

    TacsScalar dispVec[3], detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                           time, i, pt, Xpts, vars, dvars,
                                           ddvars, &detXd, dispVec);
    // project displacement along user-defined direction
    TacsScalar dispProj = 0.0;
    for (int i = 0; i < count; i++) {
      dispProj += dispVec[i] * dir[i];
    }

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar factor = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksDispSum))/dx = 1/(ksDispSum)*d(dispProj)/dx
        factor = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
      } else if (ksType == CONTINUOUS) {
        factor = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
        factor *= weight * detXd;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm;
        factor *= weight * detXd;
      }
      TacsScalar dfdq[3];
      dfdq[0] = factor * dir[0];
      dfdq[1] = factor * dir[1];
      dfdq[2] = factor * dir[2];

      element->addPointQuantitySVSens(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                      time, alpha, beta, gamma, i, pt, Xpts,
                                      vars, dvars, ddvars, dfdq, dfdu);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSKSDisplacement::getElementXptSens(
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

    TacsScalar dispVec[3], detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                           time, i, pt, Xpts, vars, dvars,
                                           ddvars, &detXd, dispVec);
    // project displacement along user-defined direction
    TacsScalar dispProj = 0.0;
    for (int i = 0; i < count; i++) {
      dispProj += dispVec[i] * dir[i];
    }

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar factor = 0.0;
      TacsScalar dfddetXd = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksDispSum))/dx = 1/(ksDispSum)*d(dispProj)/dx
        factor = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
        dfddetXd = weight * expfact / ksWeight;
        factor = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm * weight * detXd;
        dfddetXd = dispProj * fpow * invPnorm * weight;
      }

      TacsScalar dfdq[3];
      dfdq[0] = factor * dir[0];
      dfdq[1] = factor * dir[1];
      dfdq[2] = factor * dir[2];
      element->addPointQuantityXptSens(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                       time, scale, i, pt, Xpts, vars, dvars,
                                       ddvars, dfddetXd, dfdq, dfdXpts);
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSDisplacement::addElementDVSens(
    int elemIndex, TACSElement *element, double time, TacsScalar scale,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  for (int i = 0; i < element->getNumQuadraturePoints(); i++) {
    double pt[3];
    double weight = element->getQuadraturePoint(i, pt);

    TacsScalar dispVec[3], detXd = 0.0;
    int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                           time, i, pt, Xpts, vars, dvars,
                                           ddvars, &detXd, dispVec);
    // project displacement along user-defined direction
    TacsScalar dispProj = 0.0;
    for (int i = 0; i < count; i++) {
      dispProj += dispVec[i] * dir[i];
    }

    if (count >= 1) {
      // Compute the sensitivity contribution
      TacsScalar factor = 0.0;
      if (ksType == DISCRETE) {
        // d(log(ksDispSum))/dx = 1/(ksDispSum)*d(dispProj)/dx
        factor = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
      } else if (ksType == CONTINUOUS) {
        TacsScalar expfact = exp(ksWeight * (dispProj - maxDisp)) / ksDispSum;
        factor = weight * detXd * expfact;
      } else if (ksType == PNORM_DISCRETE) {
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm;
      } else if (ksType == PNORM_CONTINUOUS) {
        // Get the determinant of the Jacobian
        TacsScalar fpow =
            pow(fabs(TacsRealPart(dispProj / maxDisp)), ksWeight - 2.0);
        factor = dispProj * fpow * invPnorm * weight * detXd;
      }

      TacsScalar dfdq[3];
      dfdq[0] = factor * dir[0];
      dfdq[1] = factor * dir[1];
      dfdq[2] = factor * dir[2];
      element->addPointQuantityDVSens(elemIndex, TACS_ELEMENT_DISPLACEMENT,
                                      time, scale, i, pt, Xpts, vars, dvars,
                                      ddvars, dfdq, dvLen, dfdx);
    }
  }
}
