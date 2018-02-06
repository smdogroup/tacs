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

#include "Compliance.h"
#include "TACSAssembler.h"

/*
  The context for the TACSKSFailure function
*/
class ComplianceCtx : public TACSFunctionCtx {
 public:
  ComplianceCtx( TACSFunction *ks,
                 int maxStrains,
                 int maxNodes ){
    // Set the compliance
    compliance = 0.0;

    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];
    
    // Set the pointers into the work array
    strain = &work[0];
    stress = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
  }
  ~ComplianceCtx(){
    delete [] work;
  }

  // Data to be used for the function computation
  TacsScalar compliance;
  TacsScalar *strain, *stress;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

/*
  Initialize the Compliance class properties
*/
TACSCompliance::TACSCompliance( TACSAssembler *_tacs ): 
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::SINGLE_STAGE, 0){
  maxNumNodes = _tacs->getMaxElementNodes();
  maxNumStresses = _tacs->getMaxElementStrains();
}

TACSCompliance::~TACSCompliance(){}

const char *TACSCompliance::funcName = "TACSCompliance";

const char *TACSCompliance::functionName(){ 
  return funcName; 
}

/*
  Create the context for this function
*/
TACSFunctionCtx *TACSCompliance::createFunctionCtx(){
  return new ComplianceCtx(this, maxNumStresses, maxNumNodes);
}

/*
  Retrieve the function value
*/
TacsScalar TACSCompliance::getFunctionValue(){
  return compliance;
}

/*
  Set the compliance to zero on all MPI processes
*/
void TACSCompliance::initEvaluation( EvaluationType ftype ){
  compliance = 0.0;
}

/*
  Sum the compliance across all MPI processors
*/
void TACSCompliance::finalEvaluation( EvaluationType ftype ){
  // Distribute the values of the KS function computed on this domain    
  TacsScalar temp = compliance;
  MPI_Allreduce(&temp, &compliance, 1, TACS_MPI_TYPE,
                MPI_SUM, tacs->getMPIComm());
}

/*
  Set up the threaded execution of the evaluation function
*/
void TACSCompliance::initThread( const double tcoef,
                                 EvaluationType ftype,
                                 TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  if (ctx){
    ctx->compliance = 0.0;
  }
}

/*
  Evaluate the compliance contributed by this element
*/
void TACSCompliance::elementWiseEval( EvaluationType ftype,
                                      TACSElement *element, int elemNum,
                                      const TacsScalar Xpts[], 
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[], 
                                      const TacsScalar ddvars[],
                                      TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  if (ctx && constitutive){
    int numGauss = element->getNumGaussPts();
    int numStresses = element->numStresses();

    // Pointer to the strain/stress
    TacsScalar *strain = ctx->strain;
    TacsScalar *stress = ctx->stress;

    // With the first iteration, find the minimum over the domain
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = element->getDetJacobian(pt, Xpts);
      
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
      constitutive->calculateStress(pt, strain, stress);
      
      // Calculate the compliance
      TacsScalar SEdensity = 0.0;
      for ( int k = 0; k < numStresses; k++ ){
        SEdensity += stress[k]*strain[k];
      }
      ctx->compliance += weight*h*SEdensity;
    }
  }
}

/*
  Add the contribution from the thread to the compliance
*/
void TACSCompliance::finalThread( const double tcoef,
                                  EvaluationType ftype,
                                  TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  if (ctx){
    compliance += tcoef*ctx->compliance;
  }
}

/*
  These functions are used to determine the sensitivity of the
  function to the state variables.
*/
void TACSCompliance::getElementSVSens( double alpha, double beta, double gamma, 
                                       TacsScalar *elemSVSens, 
                                       TACSElement *element, int elemNum,
                                       const TacsScalar Xpts[], 
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[], 
                                       const TacsScalar ddvars[],
                                       TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // Zero the contribution from this element
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (ctx && constitutive){
    int numGauss = element->getNumGaussPts();
    
    // Set the stress/strain arrays
    TacsScalar *strain = ctx->strain;
    TacsScalar *stress = ctx->stress;
    
    for ( int i = 0; i < numGauss; i++ ){
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);
    
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
      constitutive->calculateStress(pt, strain, stress);
       
      // Add the sensitivity of the compliance to the strain 
      // c = e^{T} * D * e    
      // dc/du = 2.0 * e^{T} * D * de/du
      element->addStrainSVSens(elemSVSens, pt, 2.0*h*alpha, stress, 
                               Xpts, vars);
    }
  }
}

/*
  Retrieve the element contribution to the derivative of the function
  w.r.t. the element nodes
*/
void TACSCompliance::getElementXptSens( const double tcoef, 
                                        TacsScalar fXptSens[],
                                        TACSElement *element, int elemNum,
                                        const TacsScalar Xpts[], 
                                        const TacsScalar vars[],
                                        const TacsScalar dvars[], 
                                        const TacsScalar ddvars[],
                                        TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // Zero the sensitivity
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (ctx && constitutive){
    int numGauss = element->getNumGaussPts();  
    int numStresses = element->numStresses();

    // Set the stress/strain arrays
    TacsScalar *strain = ctx->strain;
    TacsScalar *stress = ctx->stress;
    TacsScalar *hXptSens = ctx->hXptSens;
    
    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

      // Scale the weight by the linear coefficient
      weight *= tcoef;

      // Add contribution to the sensitivity from the strain calculation
      element->getStrain(strain, pt, Xpts, vars);
   
      // Get the stress
      constitutive->calculateStress(pt, strain, stress);
        
      // Compute the strain energy density
      TacsScalar SEdensity = 0.0;
      for ( int k = 0; k < numStresses; k++ ){
        SEdensity += strain[k]*stress[k];
      }

      // Add the contribution from the change in geometry
      for ( int k = 0; k < 3*numNodes; k++ ){
        fXptSens[k] += weight*hXptSens[k]*SEdensity;
      }

      // Add the terms from the derivative of the strain w.r.t. nodes
      element->addStrainXptSens(fXptSens, pt, 2.0*h*weight,
                                stress, Xpts, vars);
    }
  }
}

/*
  Evaluate the derivative of the compliance w.r.t. the material
  design variables
*/
void TACSCompliance::addElementDVSens( const double tcoef, 
                                       TacsScalar *fdvSens, int numDVs,
                                       TACSElement *element, int elemNum,
                                       const TacsScalar Xpts[], 
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[], 
                                       const TacsScalar ddvars[],
                                       TACSFunctionCtx *fctx ){
  ComplianceCtx *ctx = dynamic_cast<ComplianceCtx*>(fctx);
  TACSConstitutive *constitutive = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (ctx && constitutive){
    int numGauss = element->getNumGaussPts();

    // Set the stress/strain arrays
    TacsScalar *strain = ctx->strain;
    
    for ( int i = 0; i < numGauss; i++ ){
      // Get the quadrature point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);	

      // Get the strain at the current point
      element->getStrain(strain, pt, Xpts, vars);
      constitutive->addStressDVSens(pt, strain, tcoef*h, strain, 
                                    fdvSens, numDVs);
    }
  }
}
