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
class DispIntegralCtx : public TACSFunctionCtx {
 public:
  DispIntegralCtx( TACSFunction *func,
                 int maxNodes ){
    // Allocate the working array
    work = new double[ maxNodes ];
  }
  ~DispIntegralCtx(){
    delete [] work;
  }

  // Data to be used for the function computation
  TacsScalar integral
  TacsScalar *work;
};

/*
  Initialize the TACSDisplacementIntegral class properties
*/
TACSDisplacementIntegral::TACSDisplacementIntegral( TACSAssembler *_tacs ):
  TACSFunction(_tacs){
  maxNumNodes = _tacs->getMaxElementNodes();
}

TACSDisplacementIntegral::~TACSDisplacementIntegral(){}

/*
  TACSDisplacementIntegral function name
*/
const char * TACSDisplacementIntegral::funcName = "TACSDisplacementIntegral";

/*
  Return the function name
*/
const char *TACSDisplacementIntegral::functionName(){
  return funcName;
}

/*
  Retrieve the function value
*/
TacsScalar TACSDisplacementIntegral::getFunctionValue(){
  // Compute the final value of the KS function on all processors
  return integral;
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSDisplacementIntegral::createFunctionCtx(){
  return new DispIntegralCtx(this, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSDisplacementIntegral::initEvaluation( EvaluationType ftype ){
  integral = 0.0;
}

/*
  Reduce the function values across all MPI processes
*/
void TACSDisplacementIntegral::finalEvaluation( EvaluationType ftype ){
  // Distribute the values of the KS function computed on this domain
  TacsScalar temp = integral;
  MPI_Allreduce(&temp, &integral, 1, TACS_MPI_TYPE,
                TACS_MPI_MAX, tacs->getMPIComm());
}

/*
  Initialize the context for either integration or initialization
*/
void TACSDisplacementIntegral::initThread( const double tcoef,
                                           EvaluationType ftype,
                                           TACSFunctionCtx *fctx ){
  DispIntegralCtx *ctx = dynamic_cast<DispIntegralCtx*>(fctx);
  if (ctx){
    ctx->integral = 0.0;
  }
}

/*
  Perform the element-wise evaluation of the TACSDisplacementIntegral function.
*/
void TACSDisplacementIntegral::elementWiseEval( EvaluationType ftype,
                                                TACSElement *element,
                                                int elemNum,
                                                const TacsScalar Xpts[],
                                                const TacsScalar vars[],
                                                const TacsScalar dvars[],
                                                const TacsScalar ddvars[],
                                                TACSFunctionCtx *fctx ){
  DispIntegralCtx *ctx = dynamic_cast<DispIntegralCtx*>(fctx);

  if (ctx){
    // Temporary values to hold the
    double *N = ctx->work;

    // Get the number of quadrature points for this element
    int numGauss = element->getNumGaussPts();

    // With the first iteration, find the maximum over the domain
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      element->getGaussWtsPts(i, pt);

      // Get the strain
      element->getShapeFunctions(pt, N);

      // // Scale the strain by the load factor
      // for ( int k = 0; k < numDisps; k++ ){
      //   [k] *= loadFactor;
      // }

      ctx->integral += value;
    }
  }
}

/*
  For each thread used to evaluate the function, call the
  post-evaluation code once.
*/
void TACSDisplacementIntegral::finalThread( const double tcoef,
                                            EvaluationType ftype,
                                            TACSFunctionCtx *fctx ){
  DispIntegralCtx *ctx = dynamic_cast<DispIntegralCtx*>(fctx);
  if (ctx){
    integral += tcoef*ctx->integral;
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.
*/
void TACSDisplacementIntegral::getElementSVSens( double alpha,
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
  DispIntegralCtx *ctx = dynamic_cast<DispIntegralCtx*>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  if (ctx){

  }
}

/*
  Determine the derivative of the function with respect to
  the element nodal locations
*/
void TACSDisplacementIntegral::getElementXptSens( const double tcoef,
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
void TACSDisplacementIntegral::addElementDVSens( const double tcoef,
                                                 TacsScalar *fdvSens, int numDVs,
                                                 TACSElement *element, int elemNum,
                                                 const TacsScalar Xpts[],
                                                 const TacsScalar vars[],
                                                 const TacsScalar dvars[],
                                                 const TacsScalar ddvars[],
                                                 TACSFunctionCtx *fctx ){}
