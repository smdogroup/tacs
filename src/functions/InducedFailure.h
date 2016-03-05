#ifndef TACS_INDUCED_FAILURE_H
#define TACS_INDUCED_FAILURE_H

/*
  Compute an aggregated function using an induced norm approach
  
  Copyright (c) 2014 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSFunction.h"

/*
  Compute the induced aggregation norm of the failure function. This
  is more accurate than the KS functional/p-norm functional for
  sufficiently large P. The practical performance of this method
  within the context of an optimization problem is still an open
  question. The induced aggregation can be computed as follows:

  f_i = int_{A} f * g(f, P) dA/ int_{A} g(f, P) dA

  This is always less than or equal to the supremum or infinity norm
  of the function f. This can be shown by Holder's inequality. 

  The induced aggregation technique relies on selecting a function
  g(f, P) that approximates the action of a delta function as P -->
  inifty. In addition, g(f, P) should be easily computable. Two 
  possible selectionsn are:
  
  g(f, P) = e^{P*f}

  and
  
  g(f, P) = |f|^{P}

  Both of these options are implemented below. Note that the first
  method produces an approximate of max f, while the second produces
  an approximate of max |f|.

  The first method is related to the KS functional. In particular, the
  derivative of the KS functional w.r.t. P produces the following:

  d(KS)/d(P) = 1/P*(f_i - KS)

  The sign of this factor can be used to determine whether the KS
  functional is conservative or not - it is usually not conservative,
  expect in very specific circumstances. If the sign is positive --
  for sufficiently large P -- the KS functional is not conservative,
  while if it is negative, the KS functional is conservative.
*/

class InducedFailure : public TACSFunction {
 public:
  enum InducedNormType { EXPONENTIAL, POWER, 
			 EXPONENTIAL_SQUARED, POWER_SQUARED,
			 DISCRETE_EXPONENTIAL, 
			 DISCRETE_POWER, 
			 DISCRETE_EXPONENTIAL_SQUARED,
			 DISCRETE_POWER_SQUARED };

  InducedFailure( TACSAssembler * _tacs, int _elementNums[], 
		  int _numElements, double _P );
  InducedFailure( TACSAssembler * _tacs, double _P );
  ~InducedFailure();

  // Retrieve the name of the function
  // ---------------------------------
  const char * functionName();

  // Set parameters to control how the induced failure functions are evaluated
  // -------------------------------------------------------------------------
  void setParameter( double _P );
  double getParameter();
  void setInducedType( enum InducedNormType _norm_type );
  void setLoadFactor( TacsScalar _loadFactor );
  
  // Set the special quadrature information
  // --------------------------------------
  void setQuadratureType( enum QuadratureType _quad_type );
  void setQuadratureElevation( int elev );

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
                        const TacsScalar vars[], const TacsScalar Xpts[], 
			int * iwork, TacsScalar * work );
  void postEvalThread( const int iter, int * iwork, TacsScalar * work );
  void postEval( const int iter );

  // Return the value of the function
  // --------------------------------
  TacsScalar getValue();

  // State variable sensitivity determination
  // ----------------------------------------
  int getSVSensWorkSize();
  void elementWiseSVSens( TacsScalar * elemSVSens, 
                          TACSElement * element, int elemNum,
                          const TacsScalar vars[], const TacsScalar Xpts[], 
			  TacsScalar * work );

  // Design variable sensitivity evaluation
  // --------------------------------------
  int getDVSensWorkSize();
  void elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
                          TACSElement * element, int elemNum,
                          const TacsScalar vars[], const TacsScalar Xpts[], 
			  TacsScalar * work );

  int getXptSensWorkSize();
  void elementWiseXptSens( TacsScalar fXptSens[],
			   TACSElement * element, int elemNum, 
                           const TacsScalar vars[], const TacsScalar Xpts[], 
			   TacsScalar * work );

 private:
  // The type of norm to evaluate
  InducedNormType norm_type;

  TacsScalar load_factor; // Load factor applied to the strain
  int max_nodes, max_stresses; // The max number of nodes/stresses

  // Record the quadrature scheme to use 
  enum QuadratureType quad_type;

  // Set the quadrature scheme elevation
  int scheme_elevation;

  // The name of the function
  static const char * funcName;

  TacsScalar max_fail; // The maximum failure function at a Gauss point
  TacsScalar fail_numer, fail_denom; // The numerator and denominator
  TacsScalar func_val;  // The value of the P-norm
  
  // The P in the P-norm
  double P;
};

#endif 
