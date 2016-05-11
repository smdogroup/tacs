#ifndef TACS_KS_DISPLACEMENT_H
#define TACS_KS_DISPLACEMENT_H

/*
  Compute the KS function of the displacement
  
  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSFunction.h"

/*
  The following function computes the approximate maximum displacement
  along a specified direction. 
*/
class KSDisplacement : public TACSFunction {
 public:
  static const int MAX_DISPLACEMENTS = 8;
  enum KSDisplacementType { DISCRETE, CONTINUOUS };

  KSDisplacement( TACSAssembler * _tacs, const TacsScalar d[],
                  int _elementNums[], int _numElements, 
                  double ksWeight, double alpha = 1.0 );
  KSDisplacement( TACSAssembler * _tacs, const TacsScalar d[],
                  double ksWeight, double alpha = 1.0 );
  ~KSDisplacement();

  // Retrieve the name of the function
  // ---------------------------------
  const char * functionName();

  // Set parameters for the KS function
  // ----------------------------------
  void setKSDisplacementType( enum KSDisplacementType type );
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
  TacsScalar getValue(){ return ksDisp; }

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
  enum KSDisplacementType ksType;

  // The direction for the displacement
  TacsScalar dir[MAX_DISPLACEMENTS];

  // The weight on the ks function value
  double ksWeight;

  // The integral scaling value
  double alpha;

  // The maximum number of nodes/stresses in any given element
  int maxNumNodes;

  // The name of the function
  static const char * funcName;

  // The maximum failure value, the sum of exp(ksWeight*(f[i] - maxDisp)
  // and the value of the KS function
  TacsScalar ksDispSum, ksDispWeightSum, maxDisp;
  TacsScalar ksDisp;
};

#endif 
