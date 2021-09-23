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

#include "KSDisplacement.h"
#include "TACSAssembler.h"

/*
  The context for the TACSKSDisplacement function
*/
class KSDisplacementCtx : public TACSFunctionCtx {
 public:
  KSDisplacementCtx( TACSFunction *func,
                     int maxNodes ){
    // Allocate the working array
    N = new double[ maxNodes ];
    ksSum = 0.0;
    maxValue = -1e20;
  }
  ~KSDisplacementCtx(){
    delete [] N;
  }

  // Data to be used for the function computation
  TacsScalar ksSum;
  TacsScalar maxValue;
  double *N;
};

/*
  Initialize the TACSKSDisplacement class properties
*/
TACSKSDisplacement::TACSKSDisplacement( TACSAssembler *_tacs,
                                        double _ksWeight,
                                        const TacsScalar _dir[],
                                        KSDisplacementType _ksType ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  dir[0] = _dir[0];
  dir[1] = _dir[1];
  dir[2] = _dir[2];
  ksWeight = _ksWeight;
  ksType = _ksType;
  maxValue = -1e20;
  ksSum = 0.0;
  invPnorm = 0.0;
}

TACSKSDisplacement::~TACSKSDisplacement(){}

/*
  TACSKSDisplacement function name
*/
const char * TACSKSDisplacement::funcName = "TACSKSDisplacement";

/*
  Return the function name
*/
const char *TACSKSDisplacement::functionName(){
  return funcName;
}

/*
  Set the displacement aggregate type
*/
void TACSKSDisplacement::setKSDispType( KSDisplacementType _ksType ){
  ksType = _ksType;
}

/*
  Retrieve the function value
*/
TacsScalar TACSKSDisplacement::getFunctionValue(){
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
TACSFunctionCtx *TACSKSDisplacement::createFunctionCtx(){
  return new KSDisplacementCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSDisplacement::initEvaluation( EvaluationType ftype ){
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
void TACSKSDisplacement::finalEvaluation( EvaluationType ftype ){
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
void TACSKSDisplacement::initThread( const double tcoef,
                                     EvaluationType ftype,
                                     TACSFunctionCtx *fctx ){
  KSDisplacementCtx *ctx = dynamic_cast<KSDisplacementCtx*>(fctx);
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      ctx->maxValue = -1e20;
      ctx->ksSum = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSDisplacement function.
*/
void TACSKSDisplacement::elementWiseEval( EvaluationType ftype,
                                          TACSElement *element,
                                          int elemNum,
                                          const TacsScalar Xpts[],
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[],
                                          TACSFunctionCtx *fctx ){
  KSDisplacementCtx *ctx = dynamic_cast<KSDisplacementCtx*>(fctx);
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
          if (numDisps == 1){
            value += dir[0]*N[0]*d[0];
          }
          else if (numDisps == 2){
            value += N[0]*(dir[0]*d[0] + dir[1]*d[1]);
          }
          else {
            value += N[0]*(dir[0]*d[0] + dir[1]*d[1] + dir[2]*d[2]);
          }
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
          if (numDisps == 1){
            value += dir[0]*N[0]*d[0];
          }
          else if (numDisps == 2){
            value += N[0]*(dir[0]*d[0] + dir[1]*d[1]);
          }
          else {
            value += N[0]*(dir[0]*d[0] + dir[1]*d[1] + dir[2]*d[2]);
          }
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
void TACSKSDisplacement::finalThread( const double tcoef,
                                      EvaluationType ftype,
                                      TACSFunctionCtx *fctx ){
  KSDisplacementCtx *ctx = dynamic_cast<KSDisplacementCtx*>(fctx);
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
void TACSKSDisplacement::getElementSVSens( double alpha,
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
  KSDisplacementCtx *ctx = dynamic_cast<KSDisplacementCtx*>(fctx);

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
        if (numDisps == 1){
          value += dir[0]*N[0]*d[0];
        }
        else if (numDisps == 2){
          value += N[0]*(dir[0]*d[0] + dir[1]*d[1]);
        }
        else {
          value += N[0]*(dir[0]*d[0] + dir[1]*d[1] + dir[2]*d[2]);
        }
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
        if (numDisps == 1){
          s[0] += ptWeight*dir[0]*N[0];
        }
        else if (numDisps == 2){
          s[0] += ptWeight*dir[0]*N[0];
          s[1] += ptWeight*dir[1]*N[0];
        }
        else {
          s[0] += ptWeight*dir[0]*N[0];
          s[1] += ptWeight*dir[1]*N[0];
          s[2] += ptWeight*dir[2]*N[0];
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
void TACSKSDisplacement::getElementXptSens( const double tcoef,
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
void TACSKSDisplacement::addElementDVSens( const double tcoef,
                                           TacsScalar *fdvSens, int numDVs,
                                           TACSElement *element, int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){}
