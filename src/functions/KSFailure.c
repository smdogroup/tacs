#include "KSFailure.h"
#include "TACSAssembler.h"

/*
   KS function implementation

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  The context for the TACSKSFailure function
*/
class KSFunctionCtx : public TACSFunctionCtx {
 public:
  KSFunctionCtx( TACSFunction *ks,
                 int maxStrains,
                 int maxNodes ){
    maxFail = -1e20;
    ksFailSum = 0.0;

    // Allocate the working array
    work = new TacsScalar[2*maxStrains + 3*maxNodes];
    
    // Set the pointers into the work array
    strain = &work[0];
    failSens = &work[maxStrains];
    hXptSens = &work[2*maxStrains];
  }
  ~KSFunctionCtx(){
    delete [] work;
  }

  // Data to be used for the function computation
  TacsScalar maxFail;
  TacsScalar ksFailSum;
  TacsScalar *strain;
  TacsScalar *failSens;
  TacsScalar *hXptSens;
  TacsScalar *work;
};

/*
  Initialize the TACSKSFailure class properties
*/
TACSKSFailure::TACSKSFailure( TACSAssembler *_tacs,
                              double _ksWeight, 
                              KSConstitutiveFunction func,
                              double _alpha ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN,
             TACSFunction::TWO_STAGE, 0){
  ksWeight = _ksWeight;
  alpha = _alpha;
  conType = func;
  ksType = DISCRETE;
  loadFactor = 1.0;

  // Get the max number of nodes/stresses 
  maxNumNodes = _tacs->getMaxElementNodes();
  maxNumStrains = _tacs->getMaxElementStrains();
}

TACSKSFailure::~TACSKSFailure(){}

/*
  TACSKSFailure function name
*/
const char * TACSKSFailure::funcName = "TACSKSFailure";

/*
  Set the KS aggregation type
*/
void TACSKSFailure::setKSFailureType( enum KSFailureType type ){
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double TACSKSFailure::getParameter(){
  return ksWeight; 
}

/*
  Set the KS aggregation parameter
*/
void TACSKSFailure::setParameter( double _ksWeight ){ 
  ksWeight = _ksWeight; 
}

/*
  Set the load factor to some value greater than or equal to 1.0
*/
void TACSKSFailure::setLoadFactor( TacsScalar _loadFactor ){
  if (RealPart(_loadFactor) >= 1.0){ 
    loadFactor = _loadFactor;
  }
}

/*
  Return the function name
*/
const char *TACSKSFailure::functionName(){ 
  return funcName; 
}

/*
  Retrieve the function value
*/
TacsScalar TACSKSFailure::getFunctionValue(){
  // Compute the final value of the KS function on all processors
  TacsScalar ksFail = maxFail + log(ksFailSum/alpha)/ksWeight;
  return ksFail;
}

/*
  Allocate and return the function-specific context
*/
TACSFunctionCtx *TACSKSFailure::createFunctionCtx(){
  return new KSFunctionCtx(this, maxNumStrains, maxNumNodes);
}

/*
  Initialize the internal values stored within the KS function
*/
void TACSKSFailure::initEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    maxFail = -1e20;
  }
  else if (ftype == TACSFunction::INTEGRATE){
    ksFailSum = 0.0;
  }
}

/*
  Reduce the function values across all MPI processes
*/
void TACSKSFailure::finalEvaluation( EvaluationType ftype ){
  if (ftype == TACSFunction::INITIALIZE){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = maxFail;
    MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else {
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksFailSum;
    MPI_Allreduce(&temp, &ksFailSum, 1, TACS_MPI_TYPE, 
		  MPI_SUM, tacs->getMPIComm());
  }
}

/*
  Initialize the context for either integration or initialization
*/
void TACSKSFailure::initThread( const double tcoef,
                                EvaluationType ftype,
                                TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      ctx->maxFail = -1e20;
      ctx->ksFailSum = 0.0;
    }
  }
}

/*
  Perform the element-wise evaluation of the TACSKSFailure function.
*/
void TACSKSFailure::elementWiseEval( EvaluationType ftype,
                                     TACSElement *element, int elemNum,
                                     const TacsScalar Xpts[], 
                                     const TacsScalar vars[],
                                     const TacsScalar dvars[], 
                                     const TacsScalar ddvars[],
                                     TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  if (ctx){
    // Retrieve the number of stress components for this element
    int numStresses = element->numStresses();
    
    // Get the number of quadrature points for this element
    int numGauss = element->getNumGaussPts();
    
    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();
    
    if (constitutive){
      // Set the strain buffer
      TacsScalar *strain = ctx->strain;
      
      if (ftype == TACSFunction::INITIALIZE){
        // With the first iteration, find the maximum over the domain
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          element->getGaussWtsPts(i, pt);
          
          // Get the strain
          element->getStrain(strain, pt, Xpts, vars);
       
          // Scale the strain by the load factor
          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= loadFactor;
          }      
          
          // Determine the failure criteria
          TacsScalar fail;
          if (conType == FAILURE){
            constitutive->failure(pt, strain, &fail);
          }
          else {
            constitutive->buckling(strain, &fail);
          }

          // Set the maximum failure load
          if (RealPart(fail) > RealPart(ctx->maxFail)){
            ctx->maxFail = fail;
          }
        }
      }
      else {
        for ( int i = 0; i < numGauss; i++ ){
          // Get the Gauss points one at a time
          double pt[3];
          double weight = element->getGaussWtsPts(i, pt);
          
          // Get the strain
          element->getStrain(strain, pt, Xpts, vars);
          
          // Scale the strain by the load factor
          for ( int k = 0; k < numStresses; k++ ){
            strain[k] *= loadFactor;
          }
        
          // Determine the failure criteria again
          TacsScalar fail;
          if (conType == FAILURE){
            constitutive->failure(pt, strain, &fail);
          }
          else {
            constitutive->buckling(strain, &fail);
          }

          // Add the failure load to the sum
          TacsScalar fexp = exp(ksWeight*(fail - maxFail));
          
          if (ksType == DISCRETE){
            ctx->ksFailSum += fexp;
          }
          else {
            // Get the determinant of the Jacobian
            TacsScalar h = element->getDetJacobian(pt, Xpts);
            ctx->ksFailSum += h*weight*fexp;
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
void TACSKSFailure::finalThread( const double tcoef,
                                 EvaluationType ftype,
                                 TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);
  
  if (ctx){
    if (ftype == TACSFunction::INITIALIZE){
      if (RealPart(ctx->maxFail) > RealPart(maxFail)){
        maxFail = ctx->maxFail;
      }
    }
    else {
      ksFailSum += tcoef*ctx->ksFailSum;
    }
  }
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.  
*/
void TACSKSFailure::getElementSVSens( double alpha, double beta, double gamma, 
                                      TacsScalar *elemSVSens, 
                                      TACSElement *element, int elemNum,
                                      const TacsScalar Xpts[], 
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[], 
                                      const TacsScalar ddvars[],
                                      TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Zero the derivative of the function w.r.t. the element state
  // variables
  int numVars = element->numVariables();
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  if (ctx){
    // Get the number of stress components and total number of variables
    // for this element.
    int numStresses = element->numStresses();
    
    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();
    
    // Get the constitutive object
    TACSConstitutive *constitutive = element->getConstitutive();
        
    if (constitutive){
      // Set pointers into the buffer
      TacsScalar *strain = ctx->strain;
      TacsScalar *failSens = ctx->failSens;
      
      for ( int i = 0; i < numGauss; i++ ){
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);
        
        // Get the strain
        element->getStrain(strain, pt, Xpts, vars);
        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;      
        }
        
        TacsScalar fail;
        if (conType == FAILURE){
          // Evaluate the failure criteria and its derivative 
          constitutive->failure(pt, strain, &fail);
          constitutive->failureStrainSens(pt, strain, failSens);
        }
        else {
          // Compute the buckling criteria and its derivative
          constitutive->buckling(strain, &fail);
          constitutive->bucklingStrainSens(strain, failSens);
        }

        // Compute the sensitivity contribution
        TacsScalar ksPtWeight = 0.0;
        if (ksType == DISCRETE){
          // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx 
          ksPtWeight = loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
        }
        else {
          // Get the determinant of the Jacobian
          TacsScalar h = element->getDetJacobian(pt, Xpts);
          
          ksPtWeight = h*weight*loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
        }

        // Determine the sensitivity of the state variables to SV
        element->addStrainSVSens(elemSVSens, pt, alpha*ksPtWeight, 
                                 failSens, Xpts, vars);
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void TACSKSFailure::getElementXptSens( const double tcoef, 
                                       TacsScalar fXptSens[],
                                       TACSElement *element, int elemNum,
                                       const TacsScalar Xpts[], 
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[], 
                                       const TacsScalar ddvars[],
                                       TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Zero the sensitivity w.r.t. the nodes
  int numNodes = element->numNodes();
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  if (ctx){
    // Get the number of stress components, the total number of
    // variables, and the total number of nodes
    int numStresses = element->numStresses();
    int numVars = element->numVariables();

    // Get the quadrature scheme information
    int numGauss = element->getNumGaussPts();

    // Get the constitutive object for this element
    TACSConstitutive *constitutive = element->getConstitutive();

    if (constitutive){
      // Set pointers into the buffer
      TacsScalar *strain = ctx->strain;
      TacsScalar *failSens = ctx->failSens;
      TacsScalar *hXptSens = ctx->hXptSens;
  
      for ( int i = 0; i < numGauss; i++ ){
        // Get the gauss point
        double pt[3];
        double weight = element->getGaussWtsPts(i, pt);

        // Get the strain at the current point within the element
        element->getStrain(strain, pt, Xpts, vars);

        // Multiply by the load factor
        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;
        }

        // Determine the strain failure criteria
        TacsScalar fail; 
        if (conType == FAILURE){
          // Determine the sensitivity of the failure criteria to 
          // the design variables and stresses
          constitutive->failure(pt, strain, &fail);
          constitutive->failureStrainSens(pt, strain, failSens);
        }
        else {
          constitutive->buckling(strain, &fail);
          constitutive->bucklingStrainSens(strain, failSens);
        }

        // Compute the sensitivity contribution
        if (ksType == DISCRETE){
          // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx 
          TacsScalar ksPtWeight = 
            loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
      
          element->addStrainXptSens(fXptSens, pt, tcoef*ksPtWeight, failSens,
                                    Xpts, vars);
        }
        else {
          // Get the derivative of the determinant of the Jacobian
          // w.r.t. the nodes
          TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

          // Compute the derivative of the KS functional
          TacsScalar ksExp = exp(ksWeight*(fail - maxFail))/ksFailSum;
          TacsScalar ksHptWeight = tcoef*weight*ksExp/ksWeight;
          TacsScalar ksPtWeight = h*weight*loadFactor*ksExp;

          for ( int j = 0; j < 3*numNodes; j++ ){
            fXptSens[j] += ksHptWeight*hXptSens[j];
          }

          element->addStrainXptSens(fXptSens, pt, tcoef*ksPtWeight, failSens,
                                    Xpts, vars);
        }
      }
    }
  }
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void TACSKSFailure::addElementDVSens( const double tcoef, 
                                      TacsScalar *fdvSens, int numDVs,
                                      TACSElement *element, int elemNum,
                                      const TacsScalar Xpts[], 
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[], 
                                      const TacsScalar ddvars[],
                                      TACSFunctionCtx *fctx ){
  KSFunctionCtx *ctx = dynamic_cast<KSFunctionCtx*>(fctx);

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();

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

      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= loadFactor;        
      }
    
      // Determine the strain failure criteria
      TacsScalar fail; 
      if (conType == FAILURE){
        constitutive->failure(pt, strain, &fail);
      }
      else {
        constitutive->buckling(strain, &fail);
      }

      // Add contribution from the design variable sensitivity 
      // of the failure calculation     
      // Compute the sensitivity contribution
      TacsScalar ksPtWeight = 0.0;
      if (ksType == DISCRETE){
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx 
        ksPtWeight = exp(ksWeight*(fail - maxFail))/ksFailSum;
      }
      else {
        // Get the determinant of the Jacobian
        TacsScalar h = element->getDetJacobian(pt, Xpts);
      
        ksPtWeight = h*weight*exp(ksWeight*(fail - maxFail))/ksFailSum;
      }

      // Add the derivative of the criteria w.r.t. design variables
      if (conType == FAILURE){
        constitutive->addFailureDVSens(pt, strain, tcoef*ksPtWeight,
                                       fdvSens, numDVs);
      }
      else {
        constitutive->addBucklingDVSens(strain, tcoef*ksPtWeight,
                                        fdvSens, numDVs);
      }
    }
  }
}
