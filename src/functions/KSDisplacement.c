#include "KSDisplacement.h"

/*
  KS function implementation

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * KSDisplacement::funcName = "KSDisplacement";

/*
  Initialize the KSDisplacement class properties
*/
KSDisplacement::KSDisplacement( TACSAssembler * _tacs,
                                const TacsScalar d[],
                                int _elementNums[], int _numElements, 
                                double _ksWeight, double _alpha ):
TACSFunction(_tacs, _elementNums, _numElements, _numElements, 2){ 
  ksWeight = _ksWeight;

  // Set the displacement vector
  int ndisp = tacs->getVarsPerNode();
  if (ndisp > MAX_DISPLACEMENTS){
    ndisp = MAX_DISPLACEMENTS;
  }
  for ( int k = 0; k < ndisp; k++ ){
    dir[k] = d[k];
  }

  // Set the domain scaling and the KS function type
  alpha = _alpha;
  ksType = DISCRETE;
  maxNumNodes = 0;
}

KSDisplacement::KSDisplacement( TACSAssembler * _tacs,
                                const TacsScalar d[],
                                double _ksWeight, double _alpha ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, 0, 2){ 
  ksWeight = _ksWeight;
  
  // Set the displacement vector
  int ndisp = tacs->getVarsPerNode();
  if (ndisp > MAX_DISPLACEMENTS){
    ndisp = MAX_DISPLACEMENTS;
  }
  for ( int k = 0; k < ndisp; k++ ){
    dir[k] = d[k];
  }

  // Set the domain scaling and the KS function type
  alpha = _alpha;
  ksType = DISCRETE;
  maxNumNodes = 0;
}

KSDisplacement::~KSDisplacement(){}

/*
  Set the KS aggregation type
*/
void KSDisplacement::setKSDisplacementType( enum KSDisplacementType type ){
  ksType = type;
}

/*
  Retrieve the KS aggregation weight
*/
double KSDisplacement::getParameter(){
  return ksWeight; 
}

/*
  Set the KS aggregation parameter
*/
void KSDisplacement::setParameter( double _ksWeight ){ 
  ksWeight = _ksWeight; 
}

const char * KSDisplacement::functionName(){ return funcName; }

/*
  Determine the maximum number of stresses/variables and local design
  variables
*/
void KSDisplacement::preInitialize(){
  this->initialized(); // No further initialization necessary
}

void KSDisplacement::elementWiseInitialize( TACSElement * element, 
                                            int elemNum ){
  int numNodes = element->numNodes();
  if ( numNodes > maxNumNodes ){
    maxNumNodes = numNodes;
  }
}  

void KSDisplacement::postInitialize(){}

/*
  Initialize the internal values stored within the KS function
*/
void KSDisplacement::preEval( const int iter ){
  if (iter == 0){ 
    maxDisp = -1e20;
    ksDispSum = 0.0;
    ksDispWeightSum = 0.0; 
    ksDisp = 0.0;
  }
}

/*
  Determine the size of the work + iwork arrays
*/
void KSDisplacement::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 2 + maxNumNodes;
}

/*
  Initialize the components of the work array

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
  work[2]: the sum of f[i]*exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSDisplacement::preEvalThread( const int iter, 
                                    int * iwork, 
                                    TacsScalar * work ){
  work[0] = -1e20; // Initialize the maximum failure function value
  work[1] = 0.0;   // Initialize the sum of exp(ks_weight*(f[i] - max{f[i]}))
}

/*
  Perform the element-wise evaluation of the KSDisplacement function.

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the sum of exp(ks_weight*(f[i] - max{f[i]}))
  work[2]: the sum of f[i]*exp(ks_weight*(f[i] - max{f[i]}))
*/
void KSDisplacement::elementWiseEval( const int iter, 
                                      TACSElement * element, int elemNum,
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      int * iwork, 
                                      TacsScalar * work ){
  // Retrieve the number of stress components for this element
  const int nnodes = element->numNodes();
  const int varsPerNode = element->numDisplacements();
  int ndisp = varsPerNode;
  if (ndisp > MAX_DISPLACEMENTS){
    ndisp = MAX_DISPLACEMENTS;
  }

  // Get the number of quadrature points for this element
  int numGauss = element->getNumGaussPts();
    
  // Set the shape function memroy
  double *N = reinterpret_cast<double*>(&work[2]);

  if (iter == 0){    
    // With the first iteration, find the maximum over the domain
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      element->getGaussWtsPts(i, pt);

      // Retrieve the shape functions
      element->getShapeFunctions(pt, N);

      // Compute the displacement
      TacsScalar disp = 0.0;
      const TacsScalar *v = vars;
      for ( int j = 0; j < nnodes; j++ ){
        for ( int k = 0; k < ndisp; k++ ){
          disp += N[j]*v[k]*dir[k];
        }
        v += varsPerNode;
      }
        
      // Set the maximum failure load
      if (disp > work[0]){
        work[0] = disp;
      }
    }
  }
  else {
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      double weight = element->getGaussWtsPts(i, pt);

      // Retrieve the shape functions
      element->getShapeFunctions(pt, N);

      // Compute the displacement
      TacsScalar disp = 0.0;
      const TacsScalar *v = vars;
      for ( int j = 0; j < nnodes; j++ ){
        for ( int k = 0; k < ndisp; k++ ){
          disp += N[j]*v[k]*dir[k];
        }
        v += varsPerNode;
      }
        
      // Add the failure load to the sum
      TacsScalar fexp = exp(ksWeight*(disp - maxDisp));
        
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

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void KSDisplacement::postEvalThread( const int iter,
                                     int * iwork, TacsScalar * work ){
  if (iter == 0){
    if (work[0] > maxDisp){
      maxDisp = work[0];
    }
  }
  else if (iter == 1){
    ksDispSum += work[1];
  }
}

/*
  Reduce the function values across all MPI processes
*/
void KSDisplacement::postEval( const int iter ){
  if (iter == 0){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = maxDisp;
    MPI_Allreduce(&temp, &maxDisp, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (iter == 1){
    // Find the sum of the ks contributions from all processes
    TacsScalar temp = ksDispSum;
    MPI_Allreduce(&temp, &ksDispSum, 1, TACS_MPI_TYPE, 
		  MPI_SUM, tacs->getMPIComm());

    printf("maxDisp = %15.5f\n", maxDisp);

    // Compute the final value of the KS function on all processors
    ksDisp = maxDisp + log(ksDispSum/alpha)/ksWeight;
  }
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int KSDisplacement::getSVSensWorkSize(){
  return maxNumNodes;
}

/*
  These functions are used to determine the sensitivity of the
  function with respect to the state variables.  
*/
void KSDisplacement::elementWiseSVSens( TacsScalar * elemSVSens, 
                                        TACSElement * element, int elemNum,
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        TacsScalar * work ){
  // Retrieve the number of stress components for this element
  const int nnodes = element->numNodes();
  const int varsPerNode = element->numDisplacements();
  int ndisp = varsPerNode;
  if (ndisp > MAX_DISPLACEMENTS){
    ndisp = MAX_DISPLACEMENTS;
  }

  // Get the number of quadrature points for this element
  int numGauss = element->getNumGaussPts();
    
  // Set the shape function memroy
  double *N = reinterpret_cast<double*>(work);

  // Zero the derivative of the function w.r.t. the element state variables
  memset(elemSVSens, 0, element->numVariables()*sizeof(TacsScalar));

  // With the first iteration, find the maximum over the domain
  for ( int i = 0; i < numGauss; i++ ){
    // Get the Gauss points one at a time
    double pt[3];
    double weight = element->getGaussWtsPts(i, pt);

    // Retrieve the shape functions
    element->getShapeFunctions(pt, N);
    
    // Compute the displacement
    TacsScalar disp = 0.0;
    const TacsScalar *v = vars;
    for ( int j = 0; j < nnodes; j++ ){
      for ( int k = 0; k < ndisp; k++ ){
        disp += N[j]*v[k]*dir[k];
      }
      v += varsPerNode;
    }
    
    // Compute the sensitivity contribution
    TacsScalar ksPtWeight = 0.0;
    if (ksType == DISCRETE){
      // d(log(ksDispSum))/dx = 1/(ksDispSum)*d(fail)/dx 
      ksPtWeight = exp(ksWeight*(disp - maxDisp))/ksDispSum;
    }
    else {
      // Get the determinant of the Jacobian
      TacsScalar h = element->getDetJacobian(pt, Xpts);      
      ksPtWeight = h*weight*exp(ksWeight*(disp - maxDisp))/ksDispSum;
    }

    // Add the sensitivity
    TacsScalar *sens = elemSVSens;
    for ( int j = 0; j < nnodes; j++ ){
      for ( int k = 0; k < ndisp; k++ ){
        sens[k] += ksPtWeight*N[j]*dir[k];
      }
      sens += varsPerNode;
    }
  }
}

/*
  Return the size of the work array for XptSens function
*/
int KSDisplacement::getXptSensWorkSize(){ 
  return 0;
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void KSDisplacement::elementWiseXptSens( TacsScalar fXptSens[],
                                         TACSElement * element, int elemNum,
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         TacsScalar * work ){
  // Zero the sensitivity w.r.t. the nodes
  int nnodes = element->numNodes();
  memset(fXptSens, 0, 3*nnodes*sizeof(TacsScalar));
}

/*
  Return the size of the maximum work array required for DVSens
*/
int KSDisplacement::getDVSensWorkSize(){
  return 0;
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void KSDisplacement::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
                                        TACSElement * element, int elemNum,
                                        const TacsScalar Xpts[],
                                        const TacsScalar vars[],
                                        TacsScalar * work ){}
