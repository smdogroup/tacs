#include "KSFailure.h"
/*
   KS function implementation

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * KSFailure::funcName = "KSFailure";

/*
  Initialize the KSFailure class properties
*/
KSFailure::KSFailure( TACSAssembler * _tacs, 
		      int _elementNums[], int _numElements, 
		      double _ksWeight, double _alpha ):
TACSFunction(_tacs, _elementNums, _numElements, _numElements, 2){ 
  ksWeight = _ksWeight;
  alpha = _alpha;
  ksType = DISCRETE;
  loadFactor = 1.0;
  maxNumNodes = 0;
  maxNumStresses = 0;
}

KSFailure::KSFailure( TACSAssembler * _tacs, 
		      double _ksWeight, double _alpha ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, 0, 2){ 
  ksWeight = _ksWeight;
  alpha = _alpha;
  ksType = DISCRETE;
  loadFactor = 1.0;
  maxNumNodes = 0;
  maxNumStresses = 0;
}

KSFailure::~KSFailure(){}

/*
  Set the KS aggregation type
*/
void KSFailure::setKSFailureType( enum KSFailureType type ){
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double KSFailure::getParameter(){
  return ksWeight; 
}

/*
  Set the KS aggregation parameter
*/
void KSFailure::setParameter( double _ksWeight ){ 
  ksWeight = _ksWeight; 
}

/*
  Set the load factor to some value greater than or equal to 1.0
*/
void KSFailure::setLoadFactor( TacsScalar _loadFactor ){
  if (_loadFactor >= 1.0){ 
    loadFactor = _loadFactor;
  }
}

const char * KSFailure::functionName(){ return funcName; }

/*
  Determine the maximum number of stresses/variables and local design
  variables
*/
void KSFailure::preInitialize(){
  this->initialized(); // No further initialization necessary
}

void KSFailure::elementWiseInitialize( TACSElement * element, int elemNum ){
  int numStresses = element->numStresses();
  if ( numStresses > maxNumStresses ){
    maxNumStresses = numStresses;
  }

  int numNodes = element->numNodes();
  if ( numNodes > maxNumNodes ){
    maxNumNodes = numNodes;
  }
}  

void KSFailure::postInitialize(){}

/*
  Initialize the internal values stored within the KS function
*/
void KSFailure::preEval( const int iter ){
  if (iter == 0){ 
    maxFail = -1e20;
    ksFailSum = 0.0;
    ksFailWeightSum = 0.0; 
    ksFail = 0.0;
  }
}

/*
  Determine the size of the work + iwork arrays
*/
void KSFailure::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 2 + maxNumStresses;
}

/*
  Initialize the components of the work array

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
  work[2]: the sum of f[i]*exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSFailure::preEvalThread( const int iter, 
			       int * iwork, TacsScalar * work ){
  work[0] = -1e20; // Initialize the maximum failure function value
  work[1] = 0.0;   // Initialize the sum of exp(ks_weight*(f[i] - max{f[i]}))
}

/*
  Perform the element-wise evaluation of the KSFailure function.

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
  work[2]: the sum of f[i]*exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSFailure::elementWiseEval( const int iter, 
				 TACSElement * element, int elemNum,
				 const TacsScalar Xpts[],
				 const TacsScalar vars[],
				 int * iwork, TacsScalar * work ){
  // Retrieve the number of stress components for this element
  int numStresses = element->numStresses();

  // Get the number of quadrature points for this element
  int numGauss = element->getNumGaussPts();

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();

  if (constitutive){
    // Set the strain buffer
    TacsScalar * strain = &work[3];
    
    if (iter == 0){    
      // With the first iteration, find the maximum over the domain
      for ( int i = 0; i < numGauss; i++ ){
        // Get the Gauss points one at a time
        double pt[3];
        element->getGaussWtsPts(i, pt);
        
        // Get the strain
        element->getStrain(strain, pt, Xpts, vars);
      
        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;
        }      
        
        // Determine the failure criteria
        TacsScalar fail;
        constitutive->failure(pt, strain, &fail);
        
        // Set the maximum failure load
        if (fail > work[0]){
          work[0] = fail;
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
        
        for ( int k = 0; k < numStresses; k++ ){
          strain[k] *= loadFactor;        
        }
        
        // Determine the failure criteria again
        TacsScalar fail;
        constitutive->failure(pt, strain, &fail);
        
        // Add the failure load to the sum
        TacsScalar fexp = exp(ksWeight*(fail - maxFail));
        
        if (ksType == DISCRETE){
          work[1] += fexp;
        }
        else {
          // Get the determinant of the Jacobian
          TacsScalar h = element->getDetJacobian(pt, Xpts);
          work[1] += h*weight*fexp;
        }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void KSFailure::postEvalThread( const int iter,
				int * iwork, TacsScalar * work ){
  if (iter == 0){
    if (work[0] > maxFail){
      maxFail = work[0];
    }
  }
  else if (iter == 1){
    ksFailSum += work[1];
  }
}

/*
  Reduce the function values across all MPI processes
*/
void KSFailure::postEval( const int iter ){
  if (iter == 0){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = maxFail;
    MPI_Allreduce(&temp, &maxFail, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (iter == 1){
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksFailSum;
    MPI_Allreduce(&temp, &ksFailSum, 1, TACS_MPI_TYPE, 
		  MPI_SUM, tacs->getMPIComm());

    // Compute the final value of the KS function on all processors
    ksFail = maxFail + log(ksFailSum/alpha)/ksWeight;
  }
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int KSFailure::getSVSensWorkSize(){
  return 2*maxNumStresses;
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.  
*/
void KSFailure::elementWiseSVSens( TacsScalar * elemSVSens, 
				   TACSElement * element, int elemNum,
				   const TacsScalar Xpts[],
				   const TacsScalar vars[],
				   TacsScalar * work ){
  // Get the number of stress components and total number of variables
  // for this element.
  int numStresses = element->numStresses();
  int numVars = element->numVariables();

  // Get the quadrature scheme information
  int numGauss = element->getNumGaussPts();

  // Get the constitutive object
  TACSConstitutive *constitutive = element->getConstitutive();
  
  // Zero the derivative of the function w.r.t. the element state variables
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  if (constitutive){
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];
    TacsScalar * failSens = &work[maxNumStresses];
  
    for ( int i = 0; i < numGauss; i++ ){
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);
        
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
            
      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= loadFactor;      
      }

      // Determine the strain failure criteria
      TacsScalar fail;
      constitutive->failure(pt, strain, &fail);
    
      // Determine the sensitivity of the failure criteria to the 
      // design variables and stresses
      constitutive->failureStrainSens(pt, strain, failSens);
    
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
      element->addStrainSVSens(elemSVSens, pt, ksPtWeight, 
                               failSens, Xpts, vars);
    }
  }
}

/*
  Return the size of the work array for XptSens function
*/
int KSFailure::getXptSensWorkSize(){
  return 2*maxNumStresses + 3*maxNumNodes;
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void KSFailure::elementWiseXptSens( TacsScalar fXptSens[],
				    TACSElement * element, int elemNum,
				    const TacsScalar Xpts[],
				    const TacsScalar vars[],
				    TacsScalar * work ){
  // Get the number of stress components, the total number of
  // variables, and the total number of nodes
  int numStresses = element->numStresses();
  int numVars = element->numVariables();
  int numNodes = element->numNodes();

  // Get the quadrature scheme information
  int numGauss = element->getNumGaussPts();

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();

  // Zero the sensitivity w.r.t. the nodes
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  if (constitutive){
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];
    TacsScalar * failSens = &work[maxNumStresses];
    TacsScalar * hXptSens = &work[2*maxNumStresses];
  
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
      constitutive->failure(pt, strain, &fail);

      // Determine the sensitivity of the failure criteria to 
      // the design variables and stresses
      constitutive->failureStrainSens(pt, strain, failSens);

      // Compute the sensitivity contribution
      if (ksType == DISCRETE){
        // d(log(ksFailSum))/dx = 1/(ksFailSum)*d(fail)/dx 
        TacsScalar ksPtWeight = 
          loadFactor*exp(ksWeight*(fail - maxFail))/ksFailSum;
      
        element->addStrainXptSens(fXptSens, pt, ksPtWeight, failSens,
                                  Xpts, vars);
      }
      else {
        // Get the derivative of the determinant of the Jacobian
        // w.r.t. the nodes
        TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

        // Compute the derivative of the KS functional
        TacsScalar ksExp = exp(ksWeight*(fail - maxFail))/ksFailSum;
        TacsScalar ksHptWeight = weight*ksExp/ksWeight;
        TacsScalar ksPtWeight = h*weight*loadFactor*ksExp;

        for ( int j = 0; j < 3*numNodes; j++ ){
          fXptSens[j] += ksHptWeight*hXptSens[j];
        }

        element->addStrainXptSens(fXptSens, pt, ksPtWeight, failSens,
                                  Xpts, vars);
      }
    }
  }
}

/*
  Return the size of the maximum work array required for DVSens
*/
int KSFailure::getDVSensWorkSize(){
  return maxNumStresses;
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void KSFailure::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
				   TACSElement * element, int elemNum,
				   const TacsScalar Xpts[],
				   const TacsScalar vars[],
				   TacsScalar * work ){ 
  // Get the number of stress components, the total number of
  // variables, and the total number of nodes
  int numStresses = element->numStresses();

  // Get the quadrature scheme information
  int numGauss = element->getNumGaussPts();

  // Get the constitutive object for this element
  TACSConstitutive *constitutive = element->getConstitutive();

  if (constitutive){
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];

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
      constitutive->failure(pt, strain, &fail);

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

      constitutive->addFailureDVSens(pt, strain, ksPtWeight,
                                     fdvSens, numDVs);
    }
  }
}
  
