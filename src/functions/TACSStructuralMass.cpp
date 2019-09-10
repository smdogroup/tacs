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

#include "StructuralMass.h"
#include "TACSAssembler.h"

/*
  The context for the TACSKSFailure function
*/
class StructuralMassCtx : public TACSFunctionCtx {
 public:
  StructuralMassCtx( TACSFunction *mfunc,
                     int maxNodes ){
    mass = 0.0;
    hXptSens = new TacsScalar[ 3*maxNodes ];
  }
  ~StructuralMassCtx(){
    delete [] hXptSens;
  }

  // Data to be used for the function computation
  TacsScalar mass;
  TacsScalar *hXptSens;
};

/*
  Allocate the structural mass function
*/
TACSStructuralMass::TACSStructuralMass( TACSAssembler *_tacs ):
TACSFunction(_tacs){
  totalMass = 0.0;
  maxNumNodes = _tacs->getMaxElementNodes();
}

/*
  Destructor for the structural mass
*/
TACSStructuralMass::~TACSStructuralMass(){}

const char *TACSStructuralMass::funcName = "StructuralMass";

/*
  The structural mass function name
*/
const char *TACSStructuralMass::functionName(){
  return funcName;
}

/*
  Get the function name
*/
TacsScalar TACSStructuralMass::getFunctionValue(){
  return totalMass;
}

/*
  Create the function context
*/
TACSFunctionCtx *TACSStructuralMass::createFunctionCtx(){
  return new StructuralMassCtx(this, maxNumNodes);
}

/*
  Initialize the mass to zero
*/
void TACSStructuralMass::initEvaluation( EvaluationType ftype ){
  totalMass = 0.0;
}

/*
  Sum the mass across all MPI processes
*/
void TACSStructuralMass::finalEvaluation( EvaluationType ftype ){
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE,
                MPI_SUM, tacs->getMPIComm());
}

/*
  Initialize the context for threaded execution
*/
void TACSStructuralMass::initThread( double tcoef,
                                     EvaluationType ftype,
                                     TACSFunctionCtx *fctx ){
  StructuralMassCtx *ctx = dynamic_cast<StructuralMassCtx*>(fctx);
  if (ctx){
    ctx->mass = 0.0;
  }
}

/*
  Evaluate the mass for each element in the domain
*/
void TACSStructuralMass::elementWiseEval(  EvaluationType ftype,
                                           TACSElement *element,
                                           int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  StructuralMassCtx *ctx = dynamic_cast<StructuralMassCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // If the element does not define a constitutive class,
  // return without adding any contribution to the function
  if (ctx && constitutive){
    int numGauss = element->getNumGaussPts();
    for ( int i = 0; i < numGauss; i++ ){
      TacsScalar ptmass[6];
      double pt[3];
      double gauss_weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      constitutive->getPointwiseMass(pt, ptmass);
      ctx->mass += gauss_weight*h*ptmass[0];
    }
  }
}

/*
  Add the contribution from the mass from all threads
*/
void TACSStructuralMass::finalThread( double tcoef,
                                      EvaluationType ftype,
                                      TACSFunctionCtx *fctx ){
  StructuralMassCtx *ctx = dynamic_cast<StructuralMassCtx*>(fctx);
  if (ctx){
    totalMass += tcoef*ctx->mass;
  }
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void TACSStructuralMass::getElementXptSens( double tcoef,
                                            TacsScalar fXptSens[],
                                            TACSElement *element,
                                            int elemNum,
                                            const TacsScalar Xpts[],
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[],
                                            const TacsScalar ddvars[],
                                            TACSFunctionCtx *fctx ){
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // Get the objects
  StructuralMassCtx *ctx = dynamic_cast<StructuralMassCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // If the element does not define a constitutive class,
  // return without adding any contribution to the function
  if (ctx && constitutive){
    TacsScalar *hXptSens = ctx->hXptSens;

    int numGauss = element->getNumGaussPts();
    // Add the sensitivity due to det of the Jacobian
    for ( int i = 0; i < numGauss; i++ ){
      double pt[3]; // The gauss point
      TacsScalar gauss_weight = element->getGaussWtsPts(i, pt);
      element->getDetJacobianXptSens(hXptSens, pt, Xpts);

      TacsScalar ptmass[6];
      constitutive->getPointwiseMass(pt, ptmass);

      for ( int k = 0; k < 3*numNodes; k++ ){
        fXptSens[k] += tcoef*gauss_weight*hXptSens[k]*ptmass[0];
      }
    }
  }
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void TACSStructuralMass::addElementDVSens( double tcoef,
                                           TacsScalar *fdvSens, int numDVs,
                                           TACSElement *element, int elemNum,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[],
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  // Get the objects
  StructuralMassCtx *ctx = dynamic_cast<StructuralMassCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // If the element does not define a constitutive class,
  // return without adding any contribution to the function
  if (ctx && constitutive){
    // The coefficients on the mass moments
    TacsScalar alpha[6] = {1.0, 0.0, 0.0,
                           0.0, 0.0, 0.0};
    int numGauss = element->getNumGaussPts();

    // Add the sensitivity from the first mass moment
    for ( int i = 0; i < numGauss; i++ ){
      double pt[3];
      TacsScalar gauss_weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = element->getDetJacobian(pt, Xpts);

      alpha[0] = tcoef*gauss_weight*h;
      constitutive->addPointwiseMassDVSens(pt, alpha, fdvSens, numDVs);
    }
  }
}
