#include "InducedBuckling.h"

/*
  Copyright (c) 2012 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * InducedBuckling::funcName = "InducedBuckling";

/*
  Initialize the InducedBuckling class properties

  Evaluate the Induced only on the elements specified
*/
InducedBuckling::InducedBuckling( TACSAssembler * _tacs, 
				int _elementNums[], 
				int _numElements, double _P ):
TACSFunction(_tacs, _elementNums, _numElements, _numElements, 2){ 
  P = _P;
  norm_type = EXPONENTIAL;
  load_factor = 1.0;

  max_nodes = 0;
  max_stresses = 0;

  maxBuckling = 0.0;
  bucklingNumer = 0.0;
  bucklingDenom = 0.0;
  func_val = 0.0;
}

/*
  Evaluate the P-norm over the entire domain
*/
InducedBuckling::InducedBuckling( TACSAssembler * _tacs, double _P ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, 0, 2){ 
  P = _P;
  norm_type = EXPONENTIAL;
  load_factor = 1.0;

  max_nodes = 0;
  max_stresses = 0;

  maxBuckling = 0.0;
  bucklingNumer = 0.0;
  bucklingDenom = 0.0;
  func_val = 0.0;
}

/*
  Delete all the allocated data
*/
InducedBuckling::~InducedBuckling(){}

/*
  Set the value of P
*/
void InducedBuckling::setParameter( double _P ){
  P = _P;
}

/*
  Retrieve the value of P
*/
double InducedBuckling::getParameter(){
  return P;
}

/*
  Set the type of p-norm to use
*/
void InducedBuckling::setInducedType( enum InducedNormType type ){
  norm_type = type;
}

/*
  Set the load factor - load factor of applied before testing the
  buckling criterion 
*/
void InducedBuckling::setLoadFactor( TacsScalar _load_factor ){
  if (_load_factor >= 1.0){ 
    load_factor = _load_factor;
  }
}

/*
  Retrieve the function name
*/
const char * InducedBuckling::functionName(){ return funcName; }

/*
  Perform the pre-initialization before each function in the domain is
  called with elementWiseInitialize
*/
void InducedBuckling::preInitialize(){
  this->initialized(); // No further initialization necessary
}

/*
  Determine the maximum number of stresses and maximum number of nodes
  in any element 
*/
void InducedBuckling::elementWiseInitialize( TACSElement * element, int elemNum ){
  int numStresses = element->numStresses();
  if ( numStresses > max_stresses ){
    max_stresses = numStresses;
  }

  int numNodes = element->numNodes();
  if ( numNodes > max_nodes ){
    max_nodes = numNodes;
  }
}  

/*
  Allocate data required to evaluate the function and its derivatives
*/
void InducedBuckling::postInitialize(){}

/*
  Perform pre-evaluation operations
  
  Iter denotes whether this is the first time through evaluating the
  function or the second time through.
*/
void InducedBuckling::preEval( const int iter ){
  if (iter == 0){
    maxBuckling = -1e20;
    func_val = 0.0;
    bucklingNumer = 0.0;
    bucklingDenom = 0.0;
  }
}

/*
  Determine the size of the work + iwork arrays
*/
void InducedBuckling::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 3 + max_stresses;
}

/*
  Initialize the components of the work array

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise buckling value
  work[1]: the integral of (f, g(f, P))
  work[2]: the integral of g(f, P)

  g(f, P) = f^{P}, or exp(P*f)
*/
void InducedBuckling::preEvalThread( const int iter, 
				    int * iwork, TacsScalar * work ){
  work[0] = -1e20; // Initialize the maximum buckling function value
  work[1] = 0.0; 
  work[2] = 0.0;
}

/*
  Perform the element-wise evaluation of the InducedBuckling function.

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise buckling value
  work[1]: the integral of (f, g(f, P))
  work[2]: the integral of g(f, P)
*/
void InducedBuckling::elementWiseEval( const int iter,
				      TACSElement * element, int elemNum,
				      const TacsScalar Xpts[],
				      const TacsScalar vars[],
				      int * iwork, TacsScalar * work ){
  // Get the number of quadrature points
  int numGauss = element->getNumGaussPts();

  // Get the number of stresses
  int numStresses = element->numStresses();

  // Retrieve the constitutive object
  TACSConstitutive * constitutive = element->getConstitutive();
 
  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set pointers into the buffer
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
          strain[k] *= load_factor;
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
          strain[k] *= load_factor;
        }

        // Determine the buckling criteria again
        TacsScalar buckling;
        constitutive->buckling(strain, &buckling);

        if (norm_type == POWER || 
            norm_type == DISCRETE_POWER){
          TacsScalar fp = pow(fabs(buckling/maxBuckling), P);
          work[1] += (buckling/maxBuckling)*fp;
          work[2] += fp;
        }
        else if (norm_type == POWER_SQUARED || 
                 norm_type == DISCRETE_POWER_SQUARED){
          TacsScalar fp = pow(fabs(buckling/maxBuckling), P);
          work[1] += (buckling*buckling/maxBuckling)*fp;
          work[2] += fp;
        }
        else if (norm_type == EXPONENTIAL || 
                 norm_type == DISCRETE_EXPONENTIAL){
          TacsScalar efp = exp(P*(buckling - maxBuckling));
          work[1] += (buckling/maxBuckling)*efp;
          work[2] += efp;
        }
        else if (norm_type == EXPONENTIAL_SQUARED || 
                 norm_type == DISCRETE_EXPONENTIAL_SQUARED){
          TacsScalar efp = exp(P*(buckling - maxBuckling));
          work[1] += (buckling*buckling/maxBuckling)*efp;
          work[2] += efp;
        }
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void InducedBuckling::postEvalThread( const int iter,
                                      int * iwork, 
                                      TacsScalar * work ){
  if (iter == 0){
    if (work[0] > maxBuckling){
      maxBuckling = work[0];
    }
  }
  else if (iter == 1){
    bucklingNumer += work[1];
    bucklingDenom += work[2];
  }
}

/*
  Reduce the function values across all MPI processes
*/
void InducedBuckling::postEval( const int iter ){
  if (iter == 0){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = maxBuckling;
    MPI_Allreduce(&temp, &maxBuckling, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (iter == 1){
    // Find the sum of the ks contributions from all processes
    TacsScalar in[2], out[2];
    in[0] = bucklingNumer;
    in[1] = bucklingDenom;

    MPI_Allreduce(in, out, 2, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

    bucklingNumer = out[0];
    bucklingDenom = out[1];

    // Compute the function value
    func_val = maxBuckling*bucklingNumer/bucklingDenom;
  }
}

/*
  Get the function value
*/
TacsScalar InducedBuckling::getValue(){
  return func_val;
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int InducedBuckling::getSVSensWorkSize(){
  return 2*max_stresses;
}

/*
  Determine the derivative of the P-norm function w.r.t. the state
  variables over this element.
*/
void InducedBuckling::elementWiseSVSens( TacsScalar * elemSVSens, 
					 TACSElement * element, int elemNum,
					 const TacsScalar Xpts[],
					 const TacsScalar vars[],
					 TacsScalar * work ){
  // Retrieve the number of quadrature points
  int numGauss = element->getNumGaussPts();

  // Get the number of variables for this element
  int numVars = element->numVariables();

  // Get the number of stresses for this element
  int numStresses = element->numStresses();

  // Retrieve the constitutive object
  TACSConstitutive * constitutive = element->getConstitutive();
  
  // Zero the derivative w.r.t. the state variables
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];
    TacsScalar * bucklingSens = &work[max_stresses];

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      element->getGaussWtsPts(i, pt);
        
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
            
      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;      
      }

      // Determine the strain buckling criteria
      TacsScalar buckling;
      constitutive->buckling(strain, &buckling);
    
      // Determine the sensitivity of the buckling criteria to the 
      // design variables and stresses
      constitutive->bucklingStrainSens(strain, bucklingSens);
    
      // Compute the derivative of the induced aggregation with respect
      // to the buckling function
      TacsScalar s = 0.0;
      if (norm_type == POWER || norm_type == DISCRETE_POWER){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);

        s = ((1.0 + P)*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);    
      }
      else if (norm_type == POWER_SQUARED || 
               norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);

        s = ((2.0 + P)*buckling*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);    
      }
      else if (norm_type == EXPONENTIAL || 
               norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(buckling - maxBuckling));

        s = (((1.0 + P*buckling)*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }
      else if (norm_type == EXPONENTIAL_SQUARED || 
               norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(buckling - maxBuckling));

        s = (((2.0 + P*buckling)*buckling*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }

      // Determine the sensitivity of the state variables to SV
      element->addStrainSVSens(elemSVSens, pt, load_factor*s,
                               bucklingSens, Xpts, vars);
    }
  }
}

/*
  Return the size of the work array for XptSens function
*/
int InducedBuckling::getXptSensWorkSize(){
  return 2*max_stresses + 3*max_nodes;
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void InducedBuckling::elementWiseXptSens( TacsScalar fXptSens[],
					 TACSElement * element, int elemNum,
					 const TacsScalar Xpts[],
					 const TacsScalar vars[],
					 TacsScalar * work ){
  // Get the number of quadrature points
  int numGauss = element->getNumGaussPts();

  // Get the number of nodes
  int numNodes = element->numNodes();

  // Get the number of stresses
  int numStresses = element->numStresses();

  // Retrieve the constitutive object
  TACSConstitutive * constitutive = element->getConstitutive();

  // Zero the sensitivity w.r.t. the nodes
  memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];
    TacsScalar * bucklingSens = &work[max_stresses];

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      double pt[3];
      element->getGaussWtsPts(i, pt);

      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;
      }

      // Determine the strain buckling criteria
      TacsScalar buckling; 
      constitutive->buckling(strain, &buckling);

      // Determine the sensitivity of the buckling criteria to 
      // the design variables and stresses
      constitutive->bucklingStrainSens(strain, bucklingSens);

      // Compute the derivative of the induced aggregation with respect
      // to the buckling function
      TacsScalar s = 0.0;
      if (norm_type == POWER || norm_type == DISCRETE_POWER){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);

        s = ((1.0 + P)*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);
      }
      else if (norm_type == POWER_SQUARED || 
               norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);

        s = ((2.0 + P)*buckling*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);
      }
      else if (norm_type == EXPONENTIAL || 
               norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(buckling - maxBuckling));

        s = (((1.0 + P*buckling)*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }
      else if (norm_type == EXPONENTIAL_SQUARED || 
               norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(buckling - maxBuckling));

        s = (((2.0 + P*buckling)*buckling*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }

      element->addStrainXptSens(fXptSens, pt, s*load_factor,
                                bucklingSens, Xpts, vars);
    }
  }
}

/*
  Return the size of the maximum work array required for DVSens
*/
int InducedBuckling::getDVSensWorkSize(){
  return max_stresses;
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void InducedBuckling::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
					TACSElement * element, int elemNum,
					const TacsScalar Xpts[],
					const TacsScalar vars[],
					TacsScalar * work ){ 
  double pt[3];
  int numGauss = element->getNumGaussPts();
  int numStresses = element->numStresses();
  TACSConstitutive * constitutive = element->getConstitutive();
  
  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (constitutive){  
    // Set pointers into the buffer
    TacsScalar * strain = &work[0];

    for ( int i = 0; i < numGauss; i++ ){
      // Get the gauss point
      element->getGaussWtsPts(i, pt);
    
      // Get the strain
      element->getStrain(strain, pt, Xpts, vars);
    
      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;        
      }

      // Determine the strain buckling criteria
      TacsScalar buckling; 
      constitutive->buckling(strain, &buckling);

      // Compute the derivative of the induced aggregation with
      // respect to the buckling function
      TacsScalar s = 0.0;
      if (norm_type == POWER || norm_type == DISCRETE_POWER){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);
      
        s = ((1.0 + P)*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);
      }
      else if (norm_type == POWER_SQUARED || 
               norm_type == DISCRETE_POWER_SQUARED){
        TacsScalar g = buckling/maxBuckling;
        TacsScalar fp = pow(fabs(g), P);
      
        s = ((2.0 + P)*buckling*g*bucklingDenom - 
             P*bucklingNumer)*fp/(g*bucklingDenom*bucklingDenom);
      }
      else if (norm_type == EXPONENTIAL || 
               norm_type == DISCRETE_EXPONENTIAL){
        TacsScalar efp = exp(P*(buckling - maxBuckling));
	  
        s = (((1.0 + P*buckling)*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }
      else if (norm_type == EXPONENTIAL_SQUARED || 
               norm_type == DISCRETE_EXPONENTIAL_SQUARED){
        TacsScalar efp = exp(P*(buckling - maxBuckling));
	  
        s = (((2.0 + P*buckling)*buckling*bucklingDenom -
              P*maxBuckling*bucklingNumer)*efp)/(bucklingDenom*bucklingDenom);
      }

      // Add the contribution to the sensitivity of the buckling load
      constitutive->addBucklingDVSens(strain, s,
                                      fdvSens, numDVs);
    }
  }
}
