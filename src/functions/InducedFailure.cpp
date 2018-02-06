/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#include "InducedFailure.h"
#include "TACSAssembler.h"

/*
  The context for the TACSInducedFailure function
*/
class InducedFailureCtx : public TACSFunctionCtx {
 public:
  InducedFailureCtx( TACSFunction *ks,
                     int maxStrains,
                     int maxNodes ){
    max_fail = -1e20;
    fail_numer = 0.0;
    fail_denom = 0.0;

    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];
    
    // Set the pointers into the work array
    strain = &work[0];
    failSens = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
  }
  ~InducedFailureCtx(){
    delete [] work;
  }

  // Data to be used for the function computation
  TacsScalar max_fail;
  TacsScalar fail_numer;
  TacsScalar fail_denom;
  TacsScalar *strain;
  TacsScalar *failSens;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

/*
  Initialize the TACSInducedFailure class properties

  Evaluate the Induced only on the elements specified
*/
TACSInducedFailure::TACSInducedFailure( TACSAssembler *_tacs, 
                                        double _P,
                                        InducedConstitutiveFunction func ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){
  // Set the penalization information
  P = _P;
  norm_type = EXPONENTIAL;
  load_factor = 1.0;

  // Set the constitutive evaluation type
  con_type = func;

  // Get the maximum size of the elements from TACS
  max_nodes = _tacs->getMaxElementNodes();
  max_stresses = _tacs->getMaxElementStrains();

  // Initialize values
  max_fail = 0.0;
  fail_numer = 0.0;
  fail_denom = 0.0;
}

/*
  Delete all the allocated data
*/
TACSInducedFailure::~TACSInducedFailure(){}

/*
  The name of the function class
*/
const char * TACSInducedFailure::funcName = "TACSInducedFailure";

/*
  Set the value of P
*/
void TACSInducedFailure::setParameter( double _P ){
  P = _P;
}

/*
  Retrieve the value of P
*/
double TACSInducedFailure::getParameter(){
  return P;
}

/*
  Set the type of p-norm to use
*/
void TACSInducedFailure::setInducedType( enum InducedNormType type ){
  norm_type = type;
}

/*
  Set the load factor - load factor of applied before testing the
  failure criterion 
*/
void TACSInducedFailure::setLoadFactor( TacsScalar _load_factor ){
  if (TacsRealPart(_load_factor) >= 1.0){ 
    load_factor = _load_factor;
  }
}

/*
  Retrieve the function value
*/
TacsScalar TACSInducedFailure::getFunctionValue(){
  // Compute the final value of the function on all processors
  return max_fail*fail_numer/fail_denom;
}

/*
  Retrieve the function name
*/
const char *TACSInducedFailure::functionName(){ return funcName; }

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSInducedFailure::createFunctionCtx(){
  return new InducedFailureCtx(this, max_stresses, max_nodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSInducedFailure::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    max_fail = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    fail_numer = 0.0;
    fail_denom = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSInducedFailure::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = max_fail;
    MPI_Allreduce(&temp, &max_fail, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (ftype == TACSFunction::INTEGRATE){
    // Find the sum of the ks contributions from all processes
    TacsScalar in[2], out[2];
    in[0] = fail_numer;
    in[1] = fail_denom;

    MPI_Allreduce(in, out, 2, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

    fail_numer = out[0];
    fail_denom = out[1];
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void TACSInducedFailure::initThread( const double tcoef,
                                     EvaluationType ftype,
                                     TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      ctx->max_fail = -1e20;
    }
    else if (ftype == TACSFunction::INTEGRATE){
      ctx->fail_numer = 0.0;
      ctx->fail_denom = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSInducedFailure
  function.
*/
void TACSInducedFailure::elementWiseEval( EvaluationType ftype,
                                          TACSElement *element, 
                                          int elemNum,
                                          const TacsScalar Xpts[], 
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[], 
                                          const TacsScalar ddvars[],
                                          TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  if (ctx){
    // Retrieve the number of stress components for this element
    int numStresses = element->numStresses();
    
    // Get the number of quadrature points for this element
    int numGauss = element->getNumGaussPts();
    
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    
    // If the element does not define a constitutive class, 
    // return without adding any contribution to the function
    if (constitutive){
      // Set pointers into the buffer
      TacsScalar *strain = ctx->strain;
      
      if (ftype == TACSFunction::INITIALIZE){
        // With the first iteration, find the maximum over the domain
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          element->getGaussWtsPts(i, pt);
          
          // Get the strain
          element->getStrain(strain, pt, Xpts, vars);

          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= load_factor;
          }      
          
          // Determine the failure criteria
          TacsScalar fail;
          if (con_type == FAILURE){
            constitutive->failure(pt, strain, &fail);
          }
          else {
            constitutive->buckling(strain, &fail);
          }

          // Set the maximum failure load
          if (TacsRealPart(fail) > TacsRealPart(ctx->max_fail)){
            ctx->max_fail = fail;
          }
        }
      }
      else {
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          double weight = element->getGaussWtsPts(i, pt);
          
          // Get the determinant of the Jacobian
          TacsScalar h = element->getDetJacobian(pt, Xpts);
          
          // Get the strain
          element->getStrain(strain, pt, Xpts, vars);
          
          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= load_factor;
          }

          // Determine the failure criteria again
          TacsScalar fail;
          if (con_type == FAILURE){
            constitutive->failure(pt, strain, &fail);
          }
          else {
            constitutive->buckling(strain, &fail);
          }

          if (norm_type == POWER){
            TacsScalar fp = pow(fabs(fail/max_fail), P);
            ctx->fail_numer += weight*h*(fail/max_fail)*fp;
            ctx->fail_denom += weight*h*fp;
          }
          else if (norm_type == DISCRETE_POWER){
            TacsScalar fp = pow(fabs(fail/max_fail), P);
            ctx->fail_numer += (fail/max_fail)*fp;
            ctx->fail_denom += fp;
          }
          else if (norm_type == POWER_SQUARED){
            TacsScalar fp = pow(fabs(fail/max_fail), P);
            ctx->fail_numer += weight*h*(fail*fail/max_fail)*fp;
            ctx->fail_denom += weight*h*fp;
          }
          else if (norm_type == DISCRETE_POWER_SQUARED){
            TacsScalar fp = pow(fabs(fail/max_fail), P);
            ctx->fail_numer += (fail*fail/max_fail)*fp;
            ctx->fail_denom += fp;
          }
          else if (norm_type == EXPONENTIAL){
            TacsScalar efp = exp(P*(fail - max_fail));
            ctx->fail_numer += weight*h*(fail/max_fail)*efp;
            ctx->fail_denom += weight*h*efp;
          }
          else if (norm_type == DISCRETE_EXPONENTIAL){
            TacsScalar efp = exp(P*(fail - max_fail));
            ctx->fail_numer += (fail/max_fail)*efp;
            ctx->fail_denom += efp;
          }
          else if (norm_type == EXPONENTIAL_SQUARED){
            TacsScalar efp = exp(P*(fail - max_fail));
            ctx->fail_numer += weight*h*(fail*fail/max_fail)*efp;
            ctx->fail_denom += weight*h*efp;
          }
          else if (norm_type == DISCRETE_EXPONENTIAL_SQUARED){
            TacsScalar efp = exp(P*(fail - max_fail));
            ctx->fail_numer += (fail*fail/max_fail)*efp;
            ctx->fail_denom += efp;
          }
        }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void TACSInducedFailure::finalThread( const double tcoef,
                                      EvaluationType ftype,
                                      TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      if (TacsRealPart(ctx->max_fail) > TacsRealPart(max_fail)){
        max_fail = ctx->max_fail;
      }
    }
    else if (ftype == TACSFunction::INTEGRATE){
      fail_numer += tcoef*ctx->fail_numer;
      fail_denom += tcoef*ctx->fail_denom;
    }
  }
}

/*
  Determine the derivative of the P-norm function w.r.t. the state
  variables over this element.
*/
void TACSInducedFailure::getElementSVSens( double alpha, double beta, 
                                           double gamma, 
                                           TacsScalar *elemSVSens, 
                                           TACSElement *element, 
                                           int elemNum,
                                           const TacsScalar Xpts[], 
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[], 
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  // Get the constitutive object
  TACSConstitutive *constitutive = element->getConstitutive();

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  if (ctx && constitutive){
    // Get the number of stress components and total number of
    // variables for this element.
    int numStresses = element->numStresses();

    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();

    // Set pointers into the buffer
    TacsScalar *strain = ctx->strain;
    TacsScalar *failSens = ctx->failSens;

    for ( int i = 0; i < numGauss; i++ ){
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
        
      // Get the determinant of the Jacobian
      TacsScalar h = element->getDetJacobian(pt, Xpts);
        
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
        
      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;      
      }
        
      TacsScalar fail;
      if (con_type == FAILURE){
        // Evaluate the failure criter and its derivative
        constitutive->failure(pt, strain, &fail);
        constitutive->failureStrainSens(pt, strain, failSens);
      }
      else {
        constitutive->buckling(strain, &fail);
        constitutive->bucklingStrainSens(strain, failSens);
      }

      // Compute the derivative of the induced aggregation with
      // respect to the failure function
      TacsScalar s = 0.0;
      if (norm_type == POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
          
        s = weight*h*((1.0 + P)*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);    
      }
      else if (norm_type == DISCRETE_POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
          
        s = ((1.0 + P)*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);    
      }
      else if (norm_type == POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
          
        s = weight*h*((2.0 + P)*fail*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);    
      }
      else if (norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);

        s = ((2.0 + P)*fail*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);    
      }
      else if (norm_type == EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = weight*h*(((1.0 + P*fail)*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = (((1.0 + P*fail)*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = weight*h*(((2.0 + P*fail)*fail*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = (((2.0 + P*fail)*fail*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }

      // Determine the sensitivity of the state variables to SV
      element->addStrainSVSens(elemSVSens, pt, alpha*load_factor*s,
                               failSens, Xpts, vars);
    }
  }
}

/*
  Determine the derivative of the function with respect to the element
  nodal locations 
*/
void TACSInducedFailure::getElementXptSens( const double tcoef, 
                                            TacsScalar fXptSens[],
                                            TACSElement *element, 
                                            int elemNum,
                                            const TacsScalar Xpts[], 
                                            const TacsScalar vars[],
                                            const TacsScalar dvars[], 
                                            const TacsScalar ddvars[],
                                            TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();

  // Zero the sensitivity w.r.t. the nodes
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  if (ctx && constitutive){
    // Get the number of stress components, the total number of
    // variables, and the total number of nodes
    int numStresses = element->numStresses();

    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();
    
    // Set pointers into the buffer
    TacsScalar *strain = ctx->strain;
    TacsScalar *failSens = ctx->failSens;
    TacsScalar *hXptSens = ctx->hXptSens;

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);

      // Get the derivative of the determinant of the Jacobian
      // w.r.t. the nodes
      TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;
      }

      // Determine the strain failure criteria
      TacsScalar fail; 
      if (con_type == FAILURE){
        constitutive->failure(pt, strain, &fail);
        constitutive->failureStrainSens(pt, strain, failSens);
      }
      else {
        constitutive->buckling(strain, &fail);
        constitutive->bucklingStrainSens(strain, failSens);
      }

      // Compute the derivative of the induced aggregation with respect
      // to the failure function
      TacsScalar s = 0.0;
      TacsScalar sx = 0.0;
      if (norm_type == POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);

        s = weight*h*((1.0 + P)*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);
        sx = ((fail*fail_denom - 
               max_fail*fail_numer)*weight*fp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);

        s = ((1.0 + P)*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);    
        sx = 0.0;
      }
      else if (norm_type == POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);

        s = weight*h*((2.0 + P)*fail*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);
        sx = ((fail*fail*fail_denom - 
               max_fail*fail_numer)*weight*fp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);

        s = ((2.0 + P)*fail*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);    
        sx = 0.0;
      }
      else if (norm_type == EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = weight*h*(((1.0 + P*fail)*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
        sx = ((fail*fail_denom - 
               max_fail*fail_numer)*weight*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = (((1.0 + P*fail)*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
        sx = 0.0;
      }
      else if (norm_type == EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = weight*h*(((2.0 + P*fail)*fail*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
        sx = ((fail*fail*fail_denom - 
               max_fail*fail_numer)*weight*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));

        s = (((2.0 + P*fail)*fail*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
        sx = 0.0;
      }

      // Scale the scaling factors by the coefficient
      sx *= tcoef;
      s *= tcoef;

      // Add the contribution from the derivative of the determinant
      // of the Jacobian w.r.t. the nodes
      for ( int j = 0; j < 3*numNodes; j++ ){
        fXptSens[j] += sx*hXptSens[j];
      }

      // Add the contribution from the derivative of the strain w.r.t.
      // the nodal locations
      element->addStrainXptSens(fXptSens, pt, s*load_factor,
                                failSens, Xpts, vars);
    }
  }
}

/*
  Determine the derivative of the function with respect to the design
  variables defined by the element - usually just the
  constitutive/material design variables.  
*/
void TACSInducedFailure::addElementDVSens( const double tcoef, 
                                           TacsScalar *fdvSens, int numDVs,
                                           TACSElement *element, int elemNum,
                                           const TacsScalar Xpts[], 
                                           const TacsScalar vars[],
                                           const TacsScalar dvars[], 
                                           const TacsScalar ddvars[],
                                           TACSFunctionCtx *fctx ){
  InducedFailureCtx *ctx = dynamic_cast<InducedFailureCtx*>(fctx);

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();
  
  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (ctx && constitutive){
    // Get the number of stress components, the total number of
    // variables, and the total number of nodes
    int numStresses = element->numStresses();
    
    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();
    
    // Set pointers into the buffer
    TacsScalar *strain = ctx->strain;

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
    
      // Get the determinant of the Jacobian
      TacsScalar h = element->getDetJacobian(pt, Xpts);
    
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
    
      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;        
      }

      // Determine the strain failure criteria
      TacsScalar fail;
      if (con_type == FAILURE){
        constitutive->failure(pt, strain, &fail);
      }
      else {
        constitutive->buckling(strain, &fail);
      }

      // Compute the derivative of the induced aggregation with
      // respect to the failure function
      TacsScalar s = 0.0;
      if (norm_type == POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
      
        s = weight*h*((1.0 + P)*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_POWER){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
      
        s = ((1.0 + P)*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);
      }
      else if (norm_type == POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
      
        s = weight*h*((2.0 + P)*fail*g*fail_denom - 
                      P*fail_numer)*fp/(g*fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = fail/max_fail;
        TacsScalar fp = pow(fabs(g), P);
	  
        s = ((2.0 + P)*fail*g*fail_denom - 
             P*fail_numer)*fp/(g*fail_denom*fail_denom);
      }
      else if (norm_type == EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));
	  
        s = weight*h*(((1.0 + P*fail)*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(fail - max_fail));
	  
        s = (((1.0 + P*fail)*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));
	  
        s = weight*h*(((2.0 + P*fail)*fail*fail_denom -
                       P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      else if (norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(fail - max_fail));
	  
        s = (((2.0 + P*fail)*fail*fail_denom -
              P*max_fail*fail_numer)*efp)/(fail_denom*fail_denom);
      }
      
      // Scale the coefficient by the input
      s *= tcoef;

      if (con_type == FAILURE){
        // Add the contribution to the sensitivity of the failure load
        constitutive->addFailureDVSens(pt, strain, s,
                                       fdvSens, numDVs);
      }
      else {
        constitutive->addBucklingDVSens(strain, s, fdvSens, numDVs);
      }
    }
  }
}
