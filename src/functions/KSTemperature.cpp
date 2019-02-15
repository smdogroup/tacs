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

#include "KSTemperature.h"
#include "TACSAssembler.h"

/*
  The context for the TACSKSTemperature function
*/
class KSTemperatureCtx : public TACSFunctionCtx {
 public:
  KSTemperatureCtx( TACSFunction *func,
                    int maxNodes ){
    // Allocate the working array
    N = new double[ maxNodes ];
    ksSum = 0.0;
    maxValue = -1e20;
  }
  ~KSTemperatureCtx(){
    delete [] N;
  }

  // Data to be used for the function computation
  TacsScalar ksSum;
  TacsScalar maxValue;
  double *N;
};

/*
  Initialize the TACSKSTemperature class properties
*/
TACSKSTemperature::TACSKSTemperature( TACSAssembler *_tacs,
                                      double _ksWeight,
                                      KSTemperatureType _ksType ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  ksWeight = _ksWeight;
  ksType = _ksType;
  maxValue = -1e20;
  ksSum = 0.0;
  invPnorm = 0.0;
}

TACSKSTemperature::~TACSKSTemperature(){}

/*
  TACSKSTemperature function name
*/
const char * TACSKSTemperature::funcName = "TACSKSTemperature";

/*
  Return the function name
*/
const char *TACSKSTemperature::functionName(){
  return funcName;
}

/*
  Set the displacement aggregate type
*/
void TACSKSTemperature::setKSDispType( KSTemperatureType _ksType ){
  ksType = _ksType;
}

/*
  Retrieve the function value
*/
TacsScalar TACSKSTemperature::getFunctionValue(){
  // Compute the final value of the KS function on all processors
  if (ksType == CONTINUOUS || ksType == DISCRETE){
    return maxValue + log(ksSum)/ksWeight;
  }
  else {
    return maxValue*pow(ksSum, 1.0/ksWeight);
  }
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSKSTemperature::createFunctionCtx(){
  return new KSTemperatureCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSTemperature::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    maxValue = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    ksSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSTemperature::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain
    TacsScalar temp = maxValue;
    MPI_Allreduce(&temp, &maxValue, 1, TACS_MPI_TYPE,
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksSum;
    MPI_Allreduce(&temp, &ksSum, 1, TACS_MPI_TYPE,
                  MPI_SUM, tacs->getMPIComm());

    // Compute the P-norm quantity if needed
    invPnorm = 0.0;
    if (ksType == PNORM_DISCRETE || ksType == PNORM_CONTINUOUS){
      if (ksSum != 0.0){
        invPnorm = pow(ksSum, (1.0 - ksWeight)/ksWeight);
      }
    }
  }
}

/*
  Initialize the context for either integration or initialization
*/
void TACSKSTemperature::initThread( const double tcoef,
                                     EvaluationType ftype,
                                     TACSFunctionCtx *fctx ){
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx*>(fctx);
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      ctx->maxValue = -1e20;
      ctx->ksSum = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSTemperature function.
*/
void TACSKSTemperature::elementWiseEval( EvaluationType ftype,
                                          TACSElement *element,
                                          int elemNum,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[],
                                          TACSFunctionCtx *fctx ){
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx*>(fctx);
  if (ctx){
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    if (ftype == TACSFunction::INITIALIZE){
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, ctx->N);

        // Evaluate the dot-product with the displacements
        const double *N = ctx->N;
        const TacsScalar *d = vars;

        TacsScalar value = 0.0;
        for ( int j = 0; j < numNodes; j++ ){
          value += N[0]*d[numDisps-1];
          d += numDisps;
          N++;
        }

        if (TacsRealPart(value) > TacsRealPart(ctx->maxValue)){
          ctx->maxValue = value;
        }
      }
    }
    else {
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
        element->getShapeFunctions(pt, ctx->N);

        // Evaluate the dot-product with the displacements
        const double *N = ctx->N;
        const TacsScalar *d = vars;

        TacsScalar value = 0.0;
        for ( int j = 0; j < numNodes; j++ ){
          value += N[0]*d[numDisps-1];
          d += numDisps;
          N++;
        }

        // Add up the contribution from the quadrature
        TacsScalar h = element->getDetJacobian(pt, Xpts);
        if (ksType == CONTINUOUS){
          ctx->ksSum += h*weight*exp(ksWeight*(value - maxValue));
        }
        else if (ksType == DISCRETE){
          ctx->ksSum += exp(ksWeight*(value - maxValue));
        }
        else if (ksType == PNORM_CONTINUOUS){
          ctx->ksSum += 
            h*weight*pow(fabs(TacsRealPart(value/maxValue)), ksWeight);
        }
        else if (ksType == PNORM_DISCRETE){
          ctx->ksSum += pow(fabs(TacsRealPart(value/maxValue)), ksWeight);
        }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TACSKSTemperature::finalThread( const double tcoef,
                                      EvaluationType ftype,
                                      TACSFunctionCtx *fctx ){
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx*>(fctx);
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      if (TacsRealPart(ctx->maxValue) > TacsRealPart(maxValue)){
        maxValue = ctx->maxValue;
      }
    }
    else {
      ksSum += tcoef*ctx->ksSum;
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSKSTemperature::getElementSVSens( double alpha,
                                           double beta,
                                           double gamma,
                                           TacsScalar *elemSVSens,
                                           TACSElement *element,
                                           int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  KSTemperatureCtx *ctx = dynamic_cast<KSTemperatureCtx*>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  if (ctx){
    // Get the number of quadrature points for this element
    const int numGauss = element->getNumGaussPts();
    const int numDisps = element->numDisplacements();
    const int numNodes = element->numNodes();

    // With the first iteration, find the maximum over the domain
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
      element->getShapeFunctions(pt, ctx->N);

      // Evaluate the dot-product with the displacements
      const double *N = ctx->N;
      const TacsScalar *d = vars;

      TacsScalar value = 0.0;
      for ( int j = 0; j < numNodes; j++ ){
        value += N[0]*d[numDisps-1];
        d += numDisps;
        N++;
      }

      // Add up the contribution from the quadrature
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      TacsScalar ptWeight = 0.0;

      if (ksType == CONTINUOUS){
        ptWeight = alpha*h*weight*exp(ksWeight*(value - maxValue))/ksSum;
      }
      else if (ksType == DISCRETE){
        ptWeight = alpha*exp(ksWeight*(value - maxValue))/ksSum;
      }
      else if (ksType == PNORM_CONTINUOUS){
        ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)), ksWeight-2.0);
        ptWeight *= alpha*h*weight*invPnorm;
      }
      else if (ksType == PNORM_DISCRETE){
        ptWeight = value*pow(fabs(TacsRealPart(value/maxValue)), ksWeight-2.0);
        ptWeight *= alpha*ksWeight*invPnorm;
      }

      // Reset the shape function pointer and run through the
      // element nodes again to set the derivative
      N = ctx->N;
      TacsScalar *s = elemSVSens;
      for ( int j = 0; j < numNodes; j++ ){
        s[numDisps-1] += ptWeight*N[0];
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
void TACSKSTemperature::getElementXptSens( const double tcoef,
                                            TacsScalar fXptSens[],
                                            TACSElement *element,
                                            int elemNum,
                                            const TacsScalar Xpts[],
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[],
                                            const TacsScalar ddvars[],
                                            TACSFunctionCtx *fctx ){}

/*
  Determine the derivative of the function with respect to
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSTemperature::addElementDVSens( const double tcoef,
                                           TacsScalar *fdvSens, int numDVs,
                                           TACSElement *element, int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){}
