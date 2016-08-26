#ifndef TACS_STRUCTURAL_MASS_H
#define TACS_STRUCTURAL_MASS_H

/*
  Compute the structural mass
  
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
#include "TACSFunction.h"

/*
  Evaluate the structural mass of the structure
*/
class StructuralMass : public TACSFunction {
 public:
  StructuralMass( TACSAssembler * _tacs );
  ~StructuralMass();
  
  const char *functionName();

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
  TacsScalar getValue(){ return totalMass; }

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
  // Max. number of nodes
  int maxNumNodes;

  // The total mass of all elements in the specified domain
  TacsScalar totalMass; 

  static const char * funcName;
};

#endif 
