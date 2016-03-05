#include "KSBuckling.h"
/*
   KS function implementation

  Copyright (c) 2013 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * KSBuckling::funcName = "KSBuckling";

/*
  Initialize the KSBuckling class properties
*/
KSBuckling::KSBuckling( TACSAssembler * _tacs, 
			int _elementNums[], int _numElements, 
			TacsScalar _ksWeight ):
TACSFunction(_tacs, _elementNums, _numElements, _numElements, 2){ 
  ksWeight = _ksWeight;
  loadFactor = 1.0;
  maxNumNodes = 0;
  maxNumStresses = 0;
}

KSBuckling::KSBuckling( TACSAssembler * _tacs, TacsScalar _ksWeight ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, 0, 2){ 
  ksWeight = _ksWeight;
  loadFactor = 1.0;
  maxNumNodes = 0;
  maxNumStresses = 0;
}

KSBuckling::~KSBuckling(){}

TacsScalar KSBuckling::getParameter(){
  return ksWeight; 
}

/*
  Set the KS aggregation parameter
*/
void KSBuckling::setParameter( TacsScalar _ksWeight ){ 
  ksWeight = _ksWeight; 
}

/*
  Set the load factor to some value greater than or equal to 1.0
*/
void KSBuckling::setLoadFactor( TacsScalar _loadFactor ){
  if (_loadFactor >= 1.0){ 
    loadFactor = _loadFactor;
  }
}

const char * KSBuckling::functionName(){ return funcName; }

/*
  Determine the maximum number of stresses/variables and local design
  variables
*/
void KSBuckling::preInitialize(){
  this->initialized(); // No further initialization necessary
}

void KSBuckling::elementWiseInitialize( TACSElement * element, int elemNum ){
  int numStresses = element->numStresses();
  if ( numStresses > maxNumStresses ){
    maxNumStresses = numStresses;
  }

  int numNodes = element->numNodes();
  if ( numNodes > maxNumNodes ){
    maxNumNodes = numNodes;
  }
}  

void KSBuckling::postInitialize(){}

/*
  Initialize the internal values stored within the KS function
*/
void KSBuckling::preEval( const int iter ){
  if (iter == 0){ 
    maxBuckling = -1e20;
    ksBucklingSum = 0.0;
    ksBuckling = 0.0;
  }
}

/*
  Determine the size of the work + iwork arrays
*/
void KSBuckling::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 2 + maxNumStresses;
}

/*
  Initialize the components of the work array

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise buckling value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSBuckling::preEvalThread( const int iter, 
				int * iwork, TacsScalar * work ){
  work[0] = -1e20; // Initialize the maximum buckling function value
  work[1] = 0.0;   // Initialize the sum of exp(ks_weight*(f[i] - max{f[i]}))
}

/*
  Perform the element-wise evaluation of the KSBuckling function.

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise buckling value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSBuckling::elementWiseEval( const int iter, 
				  TACSElement * element, int elemNum,
				  const TacsScalar Xpts[],
				  const TacsScalar vars[],
				  int * iwork, TacsScalar * work ){
  int numGauss = element->getNumGaussPts();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set the strain buffer
  TacsScalar * strain = &work[2];
 
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
      
      // Determine the buckling criteria
      TacsScalar buckling;
      constitutive->buckling(strain, &buckling);

      // Set the maximum buckling load
      if (buckling > work[0]){
	work[0] = buckling;
      }
    }
  }
  else {
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      element->getGaussWtsPts(i, pt);

      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= loadFactor;        
      }

      // Determine the buckling criteria again
      TacsScalar buckling;
      constitutive->buckling(strain, &buckling);

      // Add the buckling load to the sum
      work[1] += exp(ksWeight*(buckling - maxBuckling));
    }
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void KSBuckling::postEvalThread( const int iter,
				 int * iwork, TacsScalar * work ){
  if (iter == 0){
    if (work[0] > maxBuckling){
      maxBuckling = work[0];
    }
  }
  else if (iter == 1){
    ksBucklingSum += work[1];
  }
}

/*
  Reduce the function values across all MPI processes
*/
void KSBuckling::postEval( const int iter ){
  if (iter == 0){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = maxBuckling;
    MPI_Allreduce(&temp, &maxBuckling, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (iter == 1){
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksBucklingSum;
    MPI_Allreduce(&temp, &ksBucklingSum, 1, TACS_MPI_TYPE, 
                  MPI_SUM, tacs->getMPIComm());

    // Compute the final value of the KS function on all processors
    ksBuckling = maxBuckling + log(ksBucklingSum)/ksWeight;
  }
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int KSBuckling::getSVSensWorkSize(){
  return 2*maxNumStresses;
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.  
*/
void KSBuckling::elementWiseSVSens( TacsScalar * elemSVSens, 
				    TACSElement * element, int elemNum,
				    const TacsScalar Xpts[],
				    const TacsScalar vars[],
				    TacsScalar * work ){
  int numGauss = element->getNumGaussPts();
  int numVars = element->numVariables();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();
  
  // Zero the derivative of the function w.r.t. the element state variables
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set pointers into the buffer
  TacsScalar * strain = &work[0];
  TacsScalar * bucklingSens = &work[maxNumStresses];
  
  for ( int i = 0; i < numGauss; i++ ){
    double pt[3];
    element->getGaussWtsPts(i, pt);
        
    // Get the strain
    element->getStrain(strain, pt, Xpts, vars);
            
    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= loadFactor;      
    }

    // Determine the strain buckling criteria
    TacsScalar buckling;
    constitutive->buckling(strain, &buckling);
    
    // Determine the sensitivity of the buckling criteria to the 
    // design variables and stresses
    constitutive->bucklingStrainSens(strain, bucklingSens);
    
    // ksBucklingSum += exp(ksWeight*(buckling - maxBuckling));
    TacsScalar ksLocal = exp(ksWeight*(buckling - maxBuckling));
               
    // Determine the sensitivity of the state variables to SV
    element->addStrainSVSens(elemSVSens, pt, 
			     loadFactor*(ksLocal/ksBucklingSum),
			     bucklingSens, Xpts, vars);
  }
}

/*
  Return the size of the work array for XptSens function
*/
int KSBuckling::getXptSensWorkSize(){
  return 2*maxNumStresses;
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void KSBuckling::elementWiseXptSens( TacsScalar fXptSens[],
				     TACSElement * element, int elemNum,
				     const TacsScalar Xpts[],
				     const TacsScalar vars[],
				     TacsScalar * work ){
  int numGauss = element->getNumGaussPts();
  int numNodes = element->numNodes();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  // Zero the sensitivity w.r.t. the nodes
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set pointers into the buffer
  TacsScalar * strain = &work[0];
  TacsScalar * bucklingSens = &work[maxNumStresses];
  
  for ( int i = 0; i < numGauss; i++ ){
    // Get the gauss point
    double pt[3];
    element->getGaussWtsPts(i, pt);

    // Get the strain at the current point within the element
    element->getStrain(strain, pt, Xpts, vars);

    // Multiply by the load factor
    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= loadFactor;
    }

    // Determine the strain buckling criteria
    TacsScalar buckling; 
    constitutive->buckling(strain, &buckling);

    // Determine the sensitivity of the buckling criteria to 
    // the design variables and stresses
    constitutive->bucklingStrainSens(strain, bucklingSens);

    // ksBucklingSum += exp(ksWeight*(buckling - maxBuckling));
    TacsScalar ksLocal = loadFactor*exp(ksWeight*(buckling - maxBuckling));
	
    ksLocal = (ksLocal/ksBucklingSum);
    element->addStrainXptSens(fXptSens, pt, ksLocal, bucklingSens,
			      Xpts, vars);
  }
}

/*
  Return the size of the maximum work array required for DVSens
*/
int KSBuckling::getDVSensWorkSize(){
  return maxNumStresses;
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void KSBuckling::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
				    TACSElement * element, int elemNum,
				    const TacsScalar Xpts[],
				    const TacsScalar vars[],
				    TacsScalar * work ){ 
  int numGauss = element->getNumGaussPts();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set pointers into the buffer
  TacsScalar * strain = &work[0];

  for ( int i = 0; i < numGauss; i++ ){
    // Get the gauss point
    double pt[3];
    element->getGaussWtsPts(i, pt);
    
    // Get the strain
    element->getStrain(strain, pt, Xpts, vars);
    
    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= loadFactor;        
    }
    
    // Determine the strain buckling criteria
    TacsScalar buckling; 
    constitutive->buckling(strain, &buckling);
    
    // Add contribution from the design variable sensitivity 
    // of the buckling calculation     
    // ksBucklingSum += exp(ksWeight*(buckling - maxBuckling));
    TacsScalar ksLocal = exp(ksWeight*(buckling - maxBuckling));
	
    // Add the contribution to the sensitivity of the buckling load
    constitutive->addBucklingDVSens(strain, (ksLocal/ksBucklingSum),
				    fdvSens, numDVs);
  }
}
