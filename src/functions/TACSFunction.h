#ifndef TACS_FUNCTION_H
#define TACS_FUNCTION_H

class TACSFunction;

#include "TACSObject.h"
#include "TACSElement.h"
#include "TACSAssembler.h"

/*!
  Copyright (c) 2013 Graeme Kennedy. All rights reserved. 
*/


/*
  Base class for the TACSFunctionCtx. Each context is function-specific
  and is designed to store information required to 

  It's implementation is designed to be opaque to the user, but its
  data is required when evaluating the function. It is used to store
  information for each thread in the function/gradient evaluation.
*/
class TACSFunctionCtx {
 public:
  TACSFunctionCtx( TACSFunction *_me ){ me = _me; }
  TACSFunction *me;
};

/*
  TACSFunction is the base class used to calculate the values of 
  functions of interest within TACS. This class also defines the 
  methods required for gradient evaluation. This class should be
  used for objectives and constraints for optimization problems
  posed using TACS.

  Thread-safe implementation
  --------------------------
  This base class has been re-defined to enable thread-safe execution
  of certain performance-critical methods. The main issue is that the
  function class must store data internally in order to accumulate
  information that will be used to evaluate the final function value
  within a local member of TACSFunction. This is not thread-safe since
  multiple threads could update the same member.

  Instead, the code has been modified to use a buffer that is allocated
  once for each function for each active thread. This buffer is allocated
  externally so that the buffers for multiple functions can be allocated
  simultaneously within a larger chunk of memory or allocated once and used
  for multiple functions.

  The maximum buffer size for a function can be queried using the
  getMinBufferSize() call. This can only be called if the function has
  been initialized. Note that the initialization calls are NOT
  thread-safe!

  The functions that are thread-safe are:
  elementWiseEval(), elementWiseSVSens(), elementWiseDVSens() and
  elementWiseXptSens()

  The way this object works is with the following sequence of calls:

  1. If the function is not initialized (!isInitialized()), then call
  preInitialize(), for each element in the function domain call
  elementWiseInitialize(), then call postInitialize(). Note that these
  functions are not thread-safe.

  2. Obtain the size of the buffer sizes required for the function by
  calling getEvalWorkSizes( int * iwork, int * work ). Call preEval()
  once on each MPI process, call preEvalThread() once for each thread
  in a mutex-protected mode, call the thread-safe code
  elementWiseEval(), once for each element in the domain, then call
  postEvalThread() for each thread in a mutex-protected mode and
  postEval on all MPI process. Subsequent calls to getValue() must
  return the same value on all processes.

  3. To evaluate the gradient of the function w.r.t. either the state
  variables, design variables or nodes call elementWiseSVSens(),
  elementWiseDVSens() and elementWiseXptSens(). Note that these are all
  thread-safe. Note that the results must be summed across all threads/
  MPI processes. Also note that the function may use internal values stored
  from a previous function call. As a reult, it may be necessary to evaluate
  the function before evaluating the derivatives.

  Note: You cannot mix calling sequences. That is you cannot call 
  elementWiseDVSens() before finishing the ENTIRE evaluation sequence in
  2. Otherwise the work arrays will not contain the correct data.
*/
class TACSFunction : public TACSObject {
 public:
  enum FunctionDomainType { ENTIRE_DOMAIN, SUB_DOMAIN, NO_DOMAIN };
  enum FunctionEvaluationType { SINGLE_STAGE, TWO_STAGE };
  enum FunctionStageType { INITIALIZATION, INTEGRATION };

  TACSFunction( TACSAssembler *_tacs, 
		FunctionDomain _funcDomain=ENTIRE_DOMAIN,
                FunctionEvaluationType _funcEval=SINGLE_STAGE,
		int _maxElems=0 );
  virtual ~TACSFunction();

  virtual const char * functionName() = 0;
  const char * TACSObjectName();

  // Functions for setting/adjusting the domain
  // ------------------------------------------
  FunctionDomainType getDomainType();
  FunctionEvaluationType getFunctionEvalType();

  // Set the function domain by adding or setting element numbers
  // ------------------------------------------------------------
  void setDomain( int _elemNums[], int _numElems );
  void addDomain( int elemNums[], int numElems );

  // Retrieve information about the domain
  // -------------------------------------
  int getElementNums( const int **_elemNums );

  // Return associated TACSAssembler object
  // ---------------------------------------
  TACSAssembler *getTACS();

  // Get the integer/scalar work sizes required for this function
  // ------------------------------------------------------------
  virtual TACSFunctionCtx *createFuncCtx();
  virtual void destroyFuncCtx( TACSFunctionCtx **ctx );
 
  // Evaluate the function - initialize/call once per element/post-processing
  // ------------------------------------------------------------------------
  virtual void preEval( FunctionEvaluationType ftype ){}
  virtual void elementWiseEval( FunctionEvaluationType ftype,
				TACSElement *element, int elemNum,
				const TacsScalar Xpts[],
                                const TacsScalar vars[], 
				TACSFunctionCtx *ctx );
  virtual void postEval( FunctionEvaluationType ftype ){}

  // Return the value of the function
  // --------------------------------
  virtual TacsScalar getValue() = 0;

  // Calculate the sensitivity w.r.t. the state variables
  // ----------------------------------------------------
  virtual void elementWiseSVSens( TacsScalar *elemSVSens, 
				  TACSElement *element, 
				  int elemNum, 
				  const TacsScalar Xpts[], 
				  const TacsScalar vars[], 
				  TACSFunctionCtx *ctx ){
    int numVars = element->numVariables();
    memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  }

  // Calculate the sensitivity w.r.t. the design variables
  // -----------------------------------------------------
  virtual void elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
				  TACSElement *element, int elemNum,
				  const TacsScalar Xpts[],
                                  const TacsScalar vars[], 
				  TACSFunctionCtx *ctx ){}

  // Calculate the sensitivity w.r.t. the node locations
  // ---------------------------------------------------
  virtual void elementWiseXptSens( TacsScalar fXptSens[],
				   TACSElement *element, int elemNum,
				   const TacsScalar Xpts[],
                                   const TacsScalar vars[], 
				   TACSFunctionCtx *ctx ){
    int numNodes = element->numNodes();
    memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));
  }

 protected:
  TACSAssembler *tacs;
  
 private:
  // Store the function domain type
  FunctionDomainType funcDomain; 
  FunctionEvaluationType funcEvalType;

  // Store the element domain information
  int maxElems; // maximum size of currently allocated elemNums array
  int numElems; // number of elements actually stored in elemNums
  int *elemNums; // sorted array of element numbers
};

#endif
