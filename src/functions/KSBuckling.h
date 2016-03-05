#ifndef TACS_KS_BUCKLING_H
#define TACS_KS_BUCKLING_H

/*
  Compute the KS function for buckling in TACS  
  
  Copyright (c) 2013 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSFunction.h"

/*
  The following class implements the methods from TACSFunction.h
  necessary to calculate the KS function of either a stress or strain
  buckling criteria over the domain of some finite element model.

  Each class should only ever be passed to a single instance of
  TACS. If the KS function needs to be calculated for separate
  instances, this should be handled by separate instances of
  KSBuckling.

  The buckling load is calculated using the strain-based buckling
  criteria from the base Constitutive class to determine the buckling
  load.
 
  The arguments to the KSBuckling class are:

  ksWeight -- the ks weight used in the calculation

  optional: 

  elementNums, numElements -- these specify a subdomain of the TACS
  model over which the KS function should be calculated by passing in
  the element numbers and number of elements in the subdomain.

  note: if no subdomain is specified, the calculation takes place over
  all the elements in the model 
*/
class KSBuckling : public TACSFunction {
 public:
  KSBuckling( TACSAssembler * _tacs, 
	      int _elementNums[], int _numElements, 
	      TacsScalar ksWeight);
  KSBuckling( TACSAssembler * _tacs, TacsScalar ksWeight );
  ~KSBuckling();

  TacsScalar getParameter();
  void setParameter( TacsScalar _ksWeight );
  void setLoadFactor( TacsScalar _loadFactor );
  
  const char * functionName();

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
  TacsScalar getValue(){ return ksBuckling; }

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
  // The weight on the ks function value
  TacsScalar ksWeight;

  // Load factor applied to the strain
  TacsScalar loadFactor; 

  // The maximum number of nodes/stresses in any given element
  int maxNumNodes, maxNumStresses;

  // The name of the function
  static const char * funcName;

  // The maximum buckling value, the sum of exp(ksWeight*(f[i] - maxBuckling)
  // and the value of the KS function
  TacsScalar ksBucklingSum, maxBuckling;
  TacsScalar ksBuckling;
};

#endif 
