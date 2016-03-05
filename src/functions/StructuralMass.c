#include "StructuralMass.h"

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

const char * StructuralMass::funcName = "StructuralMass";

StructuralMass::StructuralMass( TACSAssembler * _tacs ):
TACSFunction(_tacs){
  totalMass = 0.0;
  maxNumNodes = 0;
}

StructuralMass::~StructuralMass(){}

const char * StructuralMass::functionName(){ return funcName; }
  
void StructuralMass::preInitialize(){
  this->initialized(); // No further initialization necessary
}

void StructuralMass::elementWiseInitialize( TACSElement * element, int elemNum ){
  int numNodes = element->numNodes();
  if ( numNodes > maxNumNodes ){
    maxNumNodes = numNodes;
  }
}  

void StructuralMass::postInitialize(){}

/*
  Get the sizes of the work arrays required for evaluation
*/
void StructuralMass::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 1;
}

/*
  Set the total mass to zero on all MPI processes
*/
void StructuralMass::preEval( const int _iter ){
  totalMass = TacsScalar(0.0);
}

/*
  Set the first element of the work array to zero - the total
  mass for all threads
*/
void StructuralMass::preEvalThread( const int iter, 
                                    int * iwork, TacsScalar * work ){
  work[0] = 0.0;
}
  
/*
  Evaluate the mass for each element in the domain
*/
void StructuralMass::elementWiseEval( const int iter, 
				      TACSElement * element, int elemNum,
				      const TacsScalar Xpts[],
				      const TacsScalar vars[], 
                                      int * iwork, TacsScalar * work ){
  int numGauss = element->getNumGaussPts();
  TACSConstitutive * material = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!material){
    return;
  }

  for ( int i = 0; i < numGauss; i++ ){
    TacsScalar ptmass[6];
    double pt[3];
    double gauss_weight = element->getGaussWtsPts(i, pt);
    TacsScalar h = element->getJacobian(pt, Xpts);
    material->pointwiseMass(pt, ptmass);

    work[0] += gauss_weight*h*ptmass[0];
  }
}

/*
  Add the contribution from the mass from all threads
*/
void StructuralMass::postEvalThread( const int iter,
                                     int * iwork, TacsScalar * work ){
  totalMass += work[0];
}

/*
  Sum the mass across all MPI processes
*/
void StructuralMass::postEval( int _iter ){
  TacsScalar temp = totalMass;
  MPI_Allreduce(&temp, &totalMass, 1, TACS_MPI_TYPE, 
                MPI_SUM, tacs->getMPIComm());
}

/*
  Determine the size of the array required for computing the 
  derivative of the Jacobian w.r.t. the nodal locations.
*/
int StructuralMass::getXptSensWorkSize(){
  return 3*maxNumNodes;
}

/*
  Determine the derivative of the mass w.r.t. the element nodal
  locations.
*/
void StructuralMass::elementWiseXptSens( TacsScalar fXptSens[],
					 TACSElement * element, int elemNum,
					 const TacsScalar Xpts[],
					 const TacsScalar vars[], 
                                         TacsScalar * work ){
  double pt[3]; // The gauss point
  int numGauss = element->getNumGaussPts();
  int numNodes = element->numNodes();
  TACSConstitutive * material = element->getConstitutive();

  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!material){
    return;
  }

  TacsScalar * hXptSens = &work[0];

  // Add the sensitivity due to the material
  TacsScalar ptmass[6];
  for ( int i = 0; i < numGauss; i++ ){
    TacsScalar gauss_weight = element->getGaussWtsPts(i, pt);
    element->getJacobianXptSens(hXptSens, pt, Xpts);

    material->pointwiseMass(pt, ptmass);

    for ( int k = 0; k < 3*numNodes; k++ ){
      fXptSens[k] += gauss_weight*hXptSens[k]*ptmass[0];
    }
  }
}

/*
  The size of the work array required for the DVSens calcs = 0
*/
int StructuralMass::getDVSensWorkSize(){
  return 0;
}

/*
  Determine the derivative of the mass w.r.t. the material
  design variables
*/
void StructuralMass::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
					TACSElement * element, int elemNum,
					const TacsScalar Xpts[],
                                        const TacsScalar vars[],  
					TacsScalar * work ){
  double pt[3];
  int numGauss = element->getNumGaussPts();
  TACSConstitutive * material = element->getConstitutive();

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!material){
    return;
  }

  // The coefficients on the mass moments
  TacsScalar alpha[6] = {1.0, 0.0, 0.0,
			 0.0, 0.0, 0.0};
  
  // Add the sensitivity from the first mass moment
  for ( int i = 0; i < numGauss; i++ ){
    TacsScalar gauss_weight = element->getGaussWtsPts(i, pt);
    TacsScalar h = element->getJacobian(pt, Xpts);
    
    alpha[0] = gauss_weight*h;
    material->addPointwiseMassDVSens(pt, alpha, fdvSens, numDVs);
  }
}
