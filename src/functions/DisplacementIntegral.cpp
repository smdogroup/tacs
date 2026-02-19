/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2018 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "DisplacementIntegral.h"

#include "TACSAssembler.h"

/*
  The context for the TACSDisplacementIntegral function
*/
class DisplacementIntCtx : public TACSFunctionCtx {
 public:
  DisplacementIntCtx(TACSFunction *func, int maxNodes) {
    // Allocate the working array
    N = new double[maxNodes];
    value = 0.0;
  }
  ~DisplacementIntCtx() { delete[] N; }

  // Data to be used for the function computation
  TacsScalar value;
  double *N;
};

/*
  Initialize the TACSDisplacementIntegral class properties
*/
TACSDisplacementIntegral::TACSDisplacementIntegral(TACSAssembler *_tacs,
                                                   const TacsScalar _dir[])
    : TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
                   TACSFunction::SINGLE_STAGE, 0) {
  maxNumNodes = _tacs->getMaxElementNodes();
  dir[0] = _dir[0];
  dir[1] = _dir[1];
  dir[2] = _dir[2];
  value = 0.0;
}

TACSDisplacementIntegral::~TACSDisplacementIntegral() {}

/*
  TACSDisplacementIntegral function name
*/
const char *TACSDisplacementIntegral::funcName = "TACSDisplacementIntegral";

/*
  Return the function name
*/
const char *TACSDisplacementIntegral::functionName() { return funcName; }

/*
  Retrieve the function value
*/
TacsScalar TACSDisplacementIntegral::getFunctionValue() { return value; }

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSDisplacementIntegral::createFunctionCtx() {
  return new DisplacementIntCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSDisplacementIntegral::initEvaluation(EvaluationType ftype) {
  value = 0.0;
}

/*
  Reduce the function values across all MPI processes
*/
void TACSDisplacementIntegral::finalEvaluation(EvaluationType ftype) {
  TacsScalar temp = value;
  MPI_Allreduce(&temp, &value, 1, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());
}

/*
  Initialize the context for either integration or initialization
*/
void TACSDisplacementIntegral::initThread(const double tcoef,
                                          EvaluationType ftype,
                                          TACSFunctionCtx *fctx) {
  DisplacementIntCtx *ctx = dynamic_cast<DisplacementIntCtx *>(fctx);
  if (ctx) {
    ctx->value = 0.0;
  }
}

/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void TACSDisplacementIntegral::elementWiseEval(
    EvaluationType ftype, TACSElement *element, int elemNum,
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  DisplacementIntCtx *ctx = dynamic_cast<DisplacementIntCtx *>(fctx);
  if (ctx) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      element->getShapeFunctions(pt, ctx->N);

      // Evaluate the dot-product with the displacements
      const double *N = ctx->N;
      const TacsScalar *d = vars;

      TacsScalar value = 0.0;
      for (int j = 0; j < numNodes; j++) {
        if (numDisps == 1) {
          value += dir[0] * N[0] * d[0];
        } else if (numDisps == 2) {
          value += N[0] * (dir[0] * d[0] + dir[1] * d[1]);
        } else {
          value += N[0] * (dir[0] * d[0] + dir[1] * d[1] + dir[2] * d[2]);
        }
        d += numDisps;
        N++;
      }

      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      h *= weight;

      ctx->value += h * value;
    }
  }
}

/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TACSDisplacementIntegral::finalThread(const double tcoef,
                                           EvaluationType ftype,
                                           TACSFunctionCtx *fctx) {
  DisplacementIntCtx *ctx = dynamic_cast<DisplacementIntCtx *>(fctx);
  if (ctx) {
    value += ctx->value;
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSDisplacementIntegral::getElementSVSens(
    double alpha, double beta, double gamma, TacsScalar *elemSVSens,
    TACSElement *element, int elemNum, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TACSFunctionCtx *fctx) {
  DisplacementIntCtx *ctx = dynamic_cast<DisplacementIntCtx *>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars * sizeof(TacsScalar));

  if (ctx) {
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for (int i = 0; i < numGauss; i++) {
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      element->getShapeFunctions(pt, ctx->N);

      // Evaluate the dot-product with the displacements
      const double *N = ctx->N;
      const TacsScalar *d = vars;

      TacsScalar value = 0.0;
      for (int j = 0; j < numNodes; j++) {
        if (numDisps == 1) {
          value += dir[0] * N[0] * d[0];
        } else if (numDisps == 2) {
          value += N[0] * (dir[0] * d[0] + dir[1] * d[1]);
        } else {
          value += N[0] * (dir[0] * d[0] + dir[1] * d[1] + dir[2] * d[2]);
        }
        d += numDisps;
        N++;
      }

      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      h *= weight;

      // Reset the shape function pointer and run through the
      // element nodes again to set the derivative
      N = ctx->N;
      TacsScalar *s = elemSVSens;
      for (int j = 0; j < numNodes; j++) {
        if (numDisps == 1) {
          s[0] += h * dir[0] * N[0];
        } else if (numDisps == 2) {
          s[0] += h * dir[0] * N[0];
          s[1] += h * dir[1] * N[0];
        } else {
          s[0] += h * dir[0] * N[0];
          s[1] += h * dir[1] * N[0];
          s[2] += h * dir[2] * N[0];
        }
        s += numDisps;
        N++;
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSDisplacementIntegral::getElementXptSens(
    const double tcoef, TacsScalar fXptSens[], TACSElement *element,
    int elemNum, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    TACSFunctionCtx *fctx) {}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSDisplacementIntegral::addElementDVSens(
    const double tcoef, TacsScalar *fdvSens, int numDVs, TACSElement *element,
    int elemNum, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[],
    TACSFunctionCtx *fctx) {}
