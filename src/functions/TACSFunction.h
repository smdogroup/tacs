#ifndef TACS_FUNCTION_H
#define TACS_FUNCTION_H

class TACSFunction;

#include "TACSObject.h"
#include "TACSElement.h"
#include "TACSAssembler.h"
#include "FElibrary.h"

/*!
  Copyright (c) 2013 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.

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
  enum FunctionDomain { ENTIRE_DOMAIN, SUB_DOMAIN, NO_DOMAIN };
  
  TACSFunction( TACSAssembler * _tacs, 
		FunctionDomain _funcDomain = ENTIRE_DOMAIN, 
		int _maxElems = 0, int _numIterations = 1 );
  TACSFunction( TACSAssembler * _tacs, int _elemNums[], int _numElems, 
		int maxElems = -1, int _numIterations = 1 );

  virtual ~TACSFunction();

  virtual const char * functionName() = 0;
  const char * TACSObjectName();

  // Functions for setting/adjusting the domain
  // ------------------------------------------
  FunctionDomain getDomain(); 
  void setDomainSize( int _maxElems );
  void setMaxDomainSize();
  void setDomain( int _elemNums[], int _numElems );
  void addDomain( int elemNums[], int numElems );
  int getNumElements();  
  int getElements( const int ** _elemNums );

  //! Return associated TACSAssembler object
  // ---------------------------------------
  TACSAssembler * getTACS();

  // Get the number of iterations required through pre/element-wise/post calls
  // -------------------------------------------------------------------------
  int getNumIterations(){ return numIterations; }

  // Test if the function has been initialized
  // -----------------------------------------
  int isInitialized() const { return initFlag; }

  // Perform initialization of the functions
  // ---------------------------------------
  virtual void preInitialize(){}
  virtual void elementWiseInitialize( TACSElement * element, int elemNum ){}
  virtual void postInitialize(){}

  // Get the integer/scalar work sizes required for this function
  // ------------------------------------------------------------
  
  // Evaluate the function - initialize/call once per element/post-processing
  // ------------------------------------------------------------------------
  virtual void getEvalWorkSizes( int * iwork, int * work );
  virtual void preEval( const int iter ){} 
  virtual void preEvalThread( const int iter, 
                              int * iwork, TacsScalar * work ){}
  virtual void elementWiseEval( const int iter, TACSElement * element, int elemNum,
				const TacsScalar elemVars[], 
				const TacsScalar Xpts[],
                                int * iwork, TacsScalar * work ){}
  virtual void postEvalThread( const int iter, 
                               int * iwork, TacsScalar * work ){}
  virtual void postEval( const int iter ){}

  // Return the value of the function
  // --------------------------------
  virtual TacsScalar getValue() = 0;

  // Calculate the sensitivity w.r.t. the state variables
  // ----------------------------------------------------
  virtual int getSVSensWorkSize();
  virtual void elementWiseSVSens( TacsScalar * elemSVSens, TACSElement * element, 
				  int elemNum, const TacsScalar elemVars[], 
				  const TacsScalar Xpts[], TacsScalar * work ){
    int numVars = element->numVariables();
    memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  }

  // Calculate the sensitivity w.r.t. the design variables
  // -----------------------------------------------------
  virtual int getDVSensWorkSize();
  virtual void elementWiseDVSens( TacsScalar fdvSens[], int numDVs,
				  TACSElement * element, int elemNum,
				  const TacsScalar elemVars[], 
				  const TacsScalar Xpts[],
                                  TacsScalar * work ){}

  virtual int getXptSensWorkSize();
  virtual void elementWiseXptSens( TacsScalar fXptSens[],
				   TACSElement * element, int elemNum,
				   const TacsScalar elemVars[], 
				   const TacsScalar Xpts[],
                                   TacsScalar * work ){
    int numNodes = element->numNodes();
    memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));
  }

 protected:
  //! The function has been initialized - no further initialization necessary 
  void initialized(){ initFlag = 1; } 

  // Number of iterations through pre/element-wise/post eval required
  // ---------------------------------------------------------------- 
  int numIterations;
  TACSAssembler * tacs;
  
 private:
  FunctionDomain funcDomain; 
  int maxElems;   // The maximum size of currently allocated elemNums array
  int numElems;   // The number of elements actually stored in elemNums
  int * elemNums; // A list of element numbers
  int initFlag;   // if (initFlag) the function has been initialized
};

#endif
