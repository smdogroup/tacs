#ifndef COMPLIANCE_H
#define COMPLIANCE_H

/*
  Calculate the compliance in the model

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSFunction.h"
#include "TACSConstitutive.h"
#include "TACSElement.h"

/*
  Evaluate the compliance of the structure. 

  This evaluates the compliance within the structure based on an
  integration of the strain energy in each element, not by the product
  of the load vector with the displacement.
*/
class Compliance : public TACSFunction {
 public:
  Compliance( TACSAssembler * _tacs, 
              int _elementNums[], int _numElements );
  Compliance( TACSAssembler * _tacs );
  ~Compliance();

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
  TacsScalar getValue(){ return compliance; }

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
  // The maximum number of nodes/stresses within any given element
  int maxNumNodes, maxNumStresses;

  // The name of the function
  static const char * funcName;

  // The compliance value
  TacsScalar compliance;  
};

#endif
