#ifndef TACS_KS_FAILURE_H
#define TACS_KS_FAILURE_H

/*
  Compute the KS function in TACS  
  
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSFunction.h"
#include "TACSConstitutive.h"
#include "TACSElement.h"

/*
  The following class implements the methods from TACSFunction.h
  necessary to calculate the KS function of either a stress or strain
  failure criteria over the domain of some finite element model.

  Each class should only ever be passed to a single instance of
  TACS. If the KS function needs to be calculated for separate
  instances, this should be handled by separate instances of
  KSFailure.

  The failure load is calculated using the strain-based failure
  criteria from the base Constitutive class which requires linear and
  constant components of the strain to determine the failure load.
 
  The arguments to the KSFailure class are:

  ksWeight:  the ks weight used in the calculation

  optional arguments: 

  elementNums, numElements -- these specify a subdomain of the TACS
  model over which the KS function should be calculated by passing in
  the element numbers and number of elements in the subdomain.

  note: if no subdomain is specified, the calculation takes place over
  all the elements in the model 
*/
class KSFailure : public TACSFunction {
 public:
  enum KSFailureType { DISCRETE, CONTINUOUS };

  KSFailure( TACSAssembler * _tacs, 
	     int _elementNums[], int _numElements, 
	     double ksWeight, double alpha = 1.0 );
  KSFailure( TACSAssembler * _tacs, 
	     double ksWeight, double alpha = 1.0 );
  ~KSFailure();

  // Retrieve the name of the function
  // ---------------------------------
  const char * functionName();

  // Set parameters for the KS function
  // ----------------------------------
  void setKSFailureType( enum KSFailureType type );
  double getParameter();
  void setParameter( double _ksWeight );
  void setLoadFactor( TacsScalar _loadFactor );

  // Functions for initialization
  // ----------------------------
  void preInitialize();
  void elementWiseInitialize( TACSElement * element, int elemNum );
  void postInitialize();

  // Functions for evaluation
  // ------------------------
  void getEvalWorkSizes( int * iwork, int * work );
  void preEval( const int iter );
  void preEvalThread( const int iter, int * iwork, TacsScalar * work );
  void elementWiseEval( const int iter, TACSElement * element, int elemNum,
                        const TacsScalar Xpts[], 
			const TacsScalar vars[],
			int * iwork, TacsScalar * work );
  void postEvalThread( const int iter, int * iwork, TacsScalar * work );
  void postEval( const int iter );

  // Return the value of the function
  // --------------------------------
  TacsScalar getValue(){ return ksFail; }

  // State variable sensitivities
  // ----------------------------
  int getSVSensWorkSize();
  void elementWiseSVSens( TacsScalar * elemSVSens, 
                          TACSElement * element, int elemNum,
                          const TacsScalar Xpts[],
			  const TacsScalar vars[], 
			  TacsScalar * work );

  // Design variable sensitivity evaluation
  // --------------------------------------
  int getDVSensWorkSize();
  void elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
                          TACSElement * element, int elemNum,
                          const TacsScalar Xpts[],
			  const TacsScalar vars[], 
			  TacsScalar * work );

  // Nodal sensitivities
  // -------------------
  int getXptSensWorkSize();
  void elementWiseXptSens( TacsScalar fXptSens[],
			   TACSElement * element, int elemNum, 
			   const TacsScalar Xpts[], 
			   const TacsScalar vars[],
			   TacsScalar * work );
  
 private:
  // The type of aggregation to use
  enum KSFailureType ksType;

  // The weight on the ks function value
  double ksWeight;

  // The integral scaling value
  double alpha;

  // Load factor applied to the strain
  TacsScalar loadFactor; 

  // The maximum number of nodes/stresses in any given element
  int maxNumNodes, maxNumStresses;

  // The name of the function
  static const char * funcName;

  // The maximum failure value, the sum of exp(ksWeight*(f[i] - maxFail)
  // and the value of the KS function
  TacsScalar ksFailSum, ksFailWeightSum, maxFail;
  TacsScalar ksFail;
};

#endif 
