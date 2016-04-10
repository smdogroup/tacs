#include "Compliance.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * Compliance::funcName = "Compliance";

/*
  Initialize the Compliance class properties
*/
Compliance::Compliance( TACSAssembler * _tacs, 
                        int _elementNums[], int _numElements ):
TACSFunction(_tacs, _elementNums, _numElements){
  maxNumNodes = 0;
  maxNumStresses = 0;
}

Compliance::Compliance( TACSAssembler * _tacs ): 
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN){ 
  maxNumNodes = 0;
  maxNumStresses = 0;
}

Compliance::~Compliance(){}

const char * Compliance::functionName(){ return funcName; }

void Compliance::preInitialize(){
  this->initialized(); // No further initialization necessary
}

void Compliance::elementWiseInitialize( TACSElement * element, int elemNum ){
  int numStresses = element->numStresses();
  if (numStresses > maxNumStresses){
    maxNumStresses = numStresses;
  }

  int numNodes = element->numNodes();
  if (numNodes > maxNumNodes){
    maxNumNodes = numNodes;
  }
}
  
void Compliance::postInitialize(){}

/*
  Get the sizes of the work arrays required for evaluation
*/
void Compliance::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 2*maxNumStresses + 1;
}

/*
  Set the compliance to zero on all MPI processes
*/
void Compliance::preEval( const int iter ){
  compliance = 0.0;
}

/*
  Set up the threaded execution of the evaluation function

  The work array contains the following information:
  work[0]: The compliance evaluated by this thread
*/
void Compliance::preEvalThread( const int iter, 
                                int * iwork, TacsScalar * work ){
  work[0] = 0.0;
}

/*
  Evaluate the compliance contributed by this element and add the
  result to work[0].
*/
void Compliance::elementWiseEval( const int iter, 
                                  TACSElement * element, int elemNum,
				  const TacsScalar Xpts[],
				  const TacsScalar vars[], 
                                  int * iwork, TacsScalar * work ){
  double pt[3];
  int numGauss = element->getNumGaussPts();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Pointer to the strain/stress
  TacsScalar * strain = &work[1];
  TacsScalar * stress = &work[1 + maxNumStresses];

  // With the first iteration, find the minimum over the domain
  for ( int i = 0; i < numGauss; i++ ){
    // Get the Gauss points one at a time
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
    work[0] += weight*h*SEdensity;
  }
}

/*
  Add the contribution from the thread to the compliance
*/
void Compliance::postEvalThread( const int iter,
                                 int * iwork, TacsScalar * work ){
  compliance += work[0];
}

/*
  Sum the compliance across all MPI processors
*/
void Compliance::postEval( const int iter ){
  // Distribute the values of the KS function computed on this domain    
  TacsScalar temp = compliance;
  MPI_Allreduce(&temp, &compliance, 1, TACS_MPI_TYPE,
                MPI_SUM, tacs->getMPIComm());
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int Compliance::getSVSensWorkSize(){
  return 2*maxNumStresses;
}

/*
  These functions are used to determine the sensitivity of the function to the
  state variables.
*/
void Compliance::elementWiseSVSens( TacsScalar * elemSVSens, 
				    TACSElement * element, int elemNum,
				    const TacsScalar Xpts[],
				    const TacsScalar vars[], 
                                    TacsScalar * work ){
  double pt[3];
  int numGauss = element->getNumGaussPts();
  int numVars = element->numVariables();
  TACSConstitutive * constitutive = element->getConstitutive();
  
  // Zero the contribution from this element
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  
  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set the stress/strain arrays
    TacsScalar * strain = &work[0];
    TacsScalar * stress = &work[maxNumStresses];
    
    for ( int i = 0; i < numGauss; i++ ){
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);
    
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
      constitutive->calculateStress(pt, strain, stress);
       
      // Add the sensitivity of the compliance to the strain c = e^{T} * D * e    
      // dc/du = 2.0 * e^{T} * D * de/du
      element->addStrainSVSens(elemSVSens, pt, 2.0*h, stress, 
                               Xpts, vars);
    }
  }
}

/*
  Return the size of the work array for the XptSens calculations
*/
int Compliance::getXptSensWorkSize(){
  return 2*maxNumStresses + 3*maxNumNodes;
}

/*
  Retrieve the element contribution to the derivative of the function
  w.r.t. the element nodes
*/
void Compliance::elementWiseXptSens( TacsScalar fXptSens[],
				     TACSElement * element, int elemNum,
				     const TacsScalar Xpts[],
				     const TacsScalar vars[],
                                     TacsScalar * work  ){
  int numGauss = element->getNumGaussPts();
  int numNodes = element->numNodes();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set the stress/strain arrays
    TacsScalar * strain = &work[0];
    TacsScalar * stress = &work[maxNumStresses];
    TacsScalar * hXptSens = &work[2*maxNumStresses];
  
    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = element->getDetJacobianXptSens(hXptSens, pt, Xpts);

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
      element->addStrainXptSens(fXptSens, pt, 2.0*h,
                                stress, Xpts, vars);
    }
  }
}

/*
  Return the size of the work array for the XptSens calculations
*/
int Compliance::getDVSensWorkSize(){
  return maxNumStresses;
}

/*
  Evaluate the derivative of the compliance w.r.t. the material
  design variables
*/
void Compliance::elementWiseDVSens( TacsScalar fdvSens[], int numDVs, 
				    TACSElement * element, int elemNum,
				    const TacsScalar Xpts[],
				    const TacsScalar vars[],
                                    TacsScalar * work ){
  int numGauss = element->getNumGaussPts();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set the stress/strain arrays
    TacsScalar * strain = &work[0];
    
    for ( int i = 0; i < numGauss; i++ ){
      // Get the quadrature point
      double pt[3];
      TacsScalar weight = element->getGaussWtsPts(i, pt);
      TacsScalar h = weight*element->getDetJacobian(pt, Xpts);	

      // Get the strain at the current point
      element->getStrain(strain, pt, Xpts, vars);
      constitutive->addStressDVSens(pt, strain, h, strain, 
                                    fdvSens, numDVs);
    }
  }
}
