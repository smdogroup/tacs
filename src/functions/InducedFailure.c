#include "InducedFailure.h"

/*
  Copyright (c) 2012 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
const char * InducedFailure::funcName = "InducedFailure";

/*
  Initialize the InducedFailure class properties

  Evaluate the Induced only on the elements specified
*/
InducedFailure::InducedFailure( TACSAssembler * _tacs, 
				int _elementNums[], 
				int _numElements, double _P ):
TACSFunction(_tacs, _elementNums, _numElements, _numElements, 2){ 
  P = _P;
  norm_type = EXPONENTIAL;
  load_factor = 1.0;

  max_nodes = 0;
  max_stresses = 0;
  scheme_elevation = 0;
  quad_type = GAUSS_QUADRATURE;

  max_fail = 0.0;
  fail_numer = 0.0;
  fail_denom = 0.0;
  func_val = 0.0;
}

/*
  Evaluate the P-norm over the entire domain
*/
InducedFailure::InducedFailure( TACSAssembler * _tacs, double _P ):
TACSFunction(_tacs, TACSFunction::ENTIRE_DOMAIN, 0, 2){ 
  P = _P;
  norm_type = EXPONENTIAL;
  load_factor = 1.0;

  max_nodes = 0;
  max_stresses = 0;
  scheme_elevation = 0;
  quad_type = GAUSS_QUADRATURE;

  max_fail = 0.0;
  fail_numer = 0.0;
  fail_denom = 0.0;
  func_val = 0.0;
}

/*
  Delete all the allocated data
*/
InducedFailure::~InducedFailure(){}

/*
  Set the value of P
*/
void InducedFailure::setParameter( double _P ){
  P = _P;
}

/*
  Retrieve the value of P
*/
double InducedFailure::getParameter(){
  return P;
}

/*
  Set the type of p-norm to use
*/
void InducedFailure::setInducedType( enum InducedNormType type ){
  norm_type = type;
}

/*
  Set the load factor - load factor of applied before testing the
  failure criterion 
*/
void InducedFailure::setLoadFactor( TacsScalar _load_factor ){
  if (_load_factor >= 1.0){ 
    load_factor = _load_factor;
  }
}

/*
  Set the type of quadrature scheme to use. This is optionally
  ignored by individual elements so be careful.
*/
void InducedFailure::setQuadratureType( enum QuadratureType _quad_type ){
  quad_type = _quad_type;
}

/*
  Set the Gauss quadrature elevation 
*/
void InducedFailure::setQuadratureElevation( int elev ){
  if (elev >= 0){
    scheme_elevation = elev;
  }
}

/*
  Retrieve the function name
*/
const char * InducedFailure::functionName(){ return funcName; }

/*
  Perform the pre-initialization before each function in the domain is
  called with elementWiseInitialize
*/
void InducedFailure::preInitialize(){
  this->initialized(); // No further initialization necessary
}

/*
  Determine the maximum number of stresses and maximum number of nodes
  in any element 
*/
void InducedFailure::elementWiseInitialize( TACSElement * element, int elemNum ){
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
void InducedFailure::postInitialize(){}

/*
  Perform pre-evaluation operations
  
  Iter denotes whether this is the first time through evaluating the
  function or the second time through.
*/
void InducedFailure::preEval( const int iter ){
  if (iter == 0){
    max_fail = -1e20;
    func_val = 0.0;
    fail_numer = 0.0;
    fail_denom = 0.0;
  }
}

/*
  Determine the size of the work + iwork arrays
*/
void InducedFailure::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0;
  *work = 3 + max_stresses;
}

/*
  Initialize the components of the work array

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the integral of (f, g(f, P))
  work[2]: the integral of g(f, P)

  g(f, P) = f^{P}, or exp(P*f)
*/
void InducedFailure::preEvalThread( const int iter, 
				    int * iwork, TacsScalar * work ){
  work[0] = -1e20; // Initialize the maximum failure function value
  work[1] = 0.0; 
  work[2] = 0.0;
}

/*
  Perform the element-wise evaluation of the InducedFailure function.

  The content of these buffers consist of the following:
  work[0]: the maximum pointwise failure value
  work[1]: the integral of exp(ks_weight*(f - max{f}))
*/
void InducedFailure::elementWiseEval( const int iter,
				      TACSElement * element, int elemNum,
				      const TacsScalar vars[],
				      const TacsScalar Xpts[],
				      int * iwork, TacsScalar * work ){
  // Retrieve the number of stress components for this element
  int numStresses = element->numStresses();

  // Retrieve the quadrature scheme
  int scheme = element->getGaussPtScheme();
  scheme += scheme_elevation;

  // Get the number of quadrature points for this element
  int numGauss = element->getNumGaussPts(scheme);

  // Get the constitutive object for this element
  TACSConstitutive * constitutive = element->getConstitutive();
 
  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set pointers into the buffer
  TacsScalar * strain = &work[3];

  if (iter == 0){    
    // With the first iteration, find the maximum over the domain
    for ( int i = 0; i < numGauss; i++ ){
      // Get the Gauss points one at a time
      double pt[3];
      element->getGaussWtsPts(quad_type, scheme, i, pt);

      // Get the strain
      element->getPtwiseStrain(strain, pt, vars, Xpts);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;
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
      double weight = element->getGaussWtsPts(quad_type, scheme, i, pt);

      // Get the determinant of the Jacobian
      TacsScalar h = element->getJacobian(pt, Xpts);

      // Get the strain
      element->getPtwiseStrain(strain, pt, vars, Xpts);

      for ( int k = 0; k < numStresses; k++ ){
        strain[k] *= load_factor;
      }

      // Determine the failure criteria again
      TacsScalar fail;
      constitutive->failure(pt, strain, &fail);

      if (norm_type == POWER){
	TacsScalar fp = pow(fabs(fail/max_fail), P);
        work[1] += weight*h*(fail/max_fail)*fp;
	work[2] += weight*h*fp;
      }
      else if (norm_type == DISCRETE_POWER){
	TacsScalar fp = pow(fabs(fail/max_fail), P);
        work[1] += (fail/max_fail)*fp;
	work[2] += fp;
      }
      else if (norm_type == POWER_SQUARED){
	TacsScalar fp = pow(fabs(fail/max_fail), P);
        work[1] += weight*h*(fail*fail/max_fail)*fp;
	work[2] += weight*h*fp;
      }
      else if (norm_type == DISCRETE_POWER_SQUARED){
	TacsScalar fp = pow(fabs(fail/max_fail), P);
        work[1] += (fail*fail/max_fail)*fp;
	work[2] += fp;
      }
      else if (norm_type == EXPONENTIAL){
	TacsScalar efp = exp(P*(fail - max_fail));
        work[1] += weight*h*(fail/max_fail)*efp;
	work[2] += weight*h*efp;
      }
      else if (norm_type == DISCRETE_EXPONENTIAL){
	TacsScalar efp = exp(P*(fail - max_fail));
        work[1] += (fail/max_fail)*efp;
	work[2] += efp;
      }
      else if (norm_type == EXPONENTIAL_SQUARED){
	TacsScalar efp = exp(P*(fail - max_fail));
        work[1] += weight*h*(fail*fail/max_fail)*efp;
	work[2] += weight*h*efp;
      }
      else if (norm_type == DISCRETE_EXPONENTIAL_SQUARED){
	TacsScalar efp = exp(P*(fail - max_fail));
        work[1] += (fail*fail/max_fail)*efp;
	work[2] += efp;
      }
    }
  }
}

/*
  For each thread used to evaluate the function, call the 
  post-evaluation code once.
*/
void InducedFailure::postEvalThread( const int iter,
				     int * iwork, 
				     TacsScalar * work ){
  if (iter == 0){
    if (work[0] > max_fail){
      max_fail = work[0];
    }
  }
  else if (iter == 1){
    fail_numer += work[1];
    fail_denom += work[2];
  }
}

/*
  Reduce the function values across all MPI processes
*/
void InducedFailure::postEval( const int iter ){
  if (iter == 0){
    // Distribute the values of the KS function computed on this domain    
    TacsScalar temp = max_fail;
    MPI_Allreduce(&temp, &max_fail, 1, TACS_MPI_TYPE, 
                  TACS_MPI_MAX, tacs->getMPIComm());
  }
  else if (iter == 1){
    // Find the sum of the ks contributions from all processes
    TacsScalar in[2], out[2];
    in[0] = fail_numer;
    in[1] = fail_denom;

    MPI_Allreduce(in, out, 2, TACS_MPI_TYPE, MPI_SUM, tacs->getMPIComm());

    fail_numer = out[0];
    fail_denom = out[1];

    // Compute the function value
    func_val = max_fail*fail_numer/fail_denom;
  }
}

/*
  Get the function value
*/
TacsScalar InducedFailure::getValue(){
  return func_val;
}

/*
  Return the size of the buffer for the state variable sensitivities
*/
int InducedFailure::getSVSensWorkSize(){
  return 2*max_stresses;
}

/*
  Determine the derivative of the P-norm function w.r.t. the state
  variables over this element.
*/
void InducedFailure::elementWiseSVSens( TacsScalar * elemSVSens, 
					TACSElement * element, int elemNum,
					const TacsScalar vars[],
					const TacsScalar Xpts[],
					TacsScalar * work ){
  // Get the number of stress components and total number of variables
  // for this element.
  int numStresses = element->numStresses();
  int numVars = element->numVariables();

  // Get the quadrature scheme information
  int scheme = element->getGaussPtScheme();
  scheme += scheme_elevation;
  int numGauss = element->getNumGaussPts(scheme);

  // Get the constitutive object
  TACSConstitutive * constitutive = element->getConstitutive();
  
  // Zero the derivative w.r.t. the state variables
  memset(elemSVSens, 0, numVars*sizeof(TacsScalar));

  // If the element does not define a constitutive class, 
  // return without adding any contribution to the function
  if (!constitutive){
    return;
  }

  // Set pointers into the buffer
  TacsScalar * strain = &work[0];
  TacsScalar * failSens = &work[max_stresses];

  for ( int i = 0; i < numGauss; i++ ){
    double pt[3];
    double weight = element->getGaussWtsPts(quad_type, scheme, i, pt);
        
    // Get the determinant of the Jacobian
    TacsScalar h = element->getJacobian(pt, Xpts);

    // Get the strain
    element->getPtwiseStrain(strain, pt, vars, Xpts);
            
    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= load_factor;      
    }

    // Determine the strain failure criteria
    TacsScalar fail;
    constitutive->failure(pt, strain, &fail);
    
    // Determine the sensitivity of the failure criteria to the 
    // design variables and stresses
    constitutive->failureStrainSens(pt, strain, failSens);
    
    // Compute the derivative of the induced aggregation with respect
    // to the failure function
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
    element->addPtwiseStrainSVSens(elemSVSens, pt, load_factor*s,
                                   failSens, vars, Xpts);
  }
}

/*
  Return the size of the work array for XptSens function
*/
int InducedFailure::getXptSensWorkSize(){
  return (2 + 3*max_nodes)*max_stresses + 3*max_nodes;
}

/*
  Determine the derivative of the function with respect to 
  the element nodal locations
*/
void InducedFailure::elementWiseXptSens( TacsScalar fXptSens[],
					 TACSElement * element, int elemNum,
					 const TacsScalar vars[],
					 const TacsScalar Xpts[],
					 TacsScalar * work ){
  // Get the number of stress components, the total number of
  // variables, and the total number of nodes
  int numStresses = element->numStresses();
  int numVars = element->numVariables();
  int numNodes = element->numNodes();

  // Get the quadrature scheme information
  int scheme = element->getGaussPtScheme();
  scheme += scheme_elevation;
  int numGauss = element->getNumGaussPts(scheme);

  // Get the constitutive object for this element
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
  TacsScalar * failSens = &work[max_stresses];
  TacsScalar * hXptSens = &work[2*max_stresses];
  TacsScalar * strainXptSens = &work[2*max_stresses + 3*max_nodes];

  for ( int i = 0; i < numGauss; i++ ){
    // Get the gauss point
    double pt[3];
    double weight = element->getGaussWtsPts(quad_type, scheme, i, pt);

    // Get the derivative of the determinant of the Jacobian w.r.t. the nodes
    TacsScalar h = element->getJacobianXptSens(hXptSens, pt, Xpts);

    // Get the strain
    element->getPtwiseStrainXptSens(strain, strainXptSens, pt,
                                    vars, Xpts);

    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= load_factor;
    }

    // Determine the strain failure criteria
    TacsScalar fail; 
    constitutive->failure(pt, strain, &fail);

    // Determine the sensitivity of the failure criteria to 
    // the design variables and stresses
    constitutive->failureStrainSens(pt, strain, failSens);

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

    for ( int j = 0; j < 3*numNodes; j++ ){
      fXptSens[j] += sx*hXptSens[j];

      for ( int k = 0; k < numStresses; k++ ){
	fXptSens[j] += 
	  s*load_factor*(failSens[k]*strainXptSens[k + j*numStresses]);
      }
    }
  }
}

/*
  Return the size of the maximum work array required for DVSens
*/
int InducedFailure::getDVSensWorkSize(){
  return max_stresses;
}

/*
  Determine the derivative of the function with respect to 
  the design variables defined by the element - usually just
  the constitutive/material design variables.
*/
void InducedFailure::elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
					TACSElement * element, int elemNum,
					const TacsScalar vars[],
					const TacsScalar Xpts[],
					TacsScalar * work ){ 
  // Get the number of stress components, the total number of
  // variables, and the total number of nodes
  int numStresses = element->numStresses();

  // Get the quadrature scheme information
  int scheme = element->getGaussPtScheme();
  scheme += scheme_elevation;
  int numGauss = element->getNumGaussPts(scheme);

  // Get the constitutive object for this element
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
    double weight = element->getGaussWtsPts(quad_type, scheme, i, pt);
    
    // Get the determinant of the Jacobian
    TacsScalar h = element->getJacobian(pt, Xpts);
    
    // Get the strain
    element->getPtwiseStrain(strain, pt, vars, Xpts);
    
    for ( int k = 0; k < numStresses; k++ ){
      strain[k] *= load_factor;        
    }

    // Determine the strain failure criteria
    TacsScalar fail; 
    constitutive->failure(pt, strain, &fail);

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

    // Add the contribution to the sensitivity of the failure load
    constitutive->addFailureDVSens(pt, strain, s,
				   fdvSens, numDVs);
  }
}
