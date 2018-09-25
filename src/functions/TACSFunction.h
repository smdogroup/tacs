/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#ifndef TACS_FUNCTION_H
#define TACS_FUNCTION_H

class TACSAssembler;

#include "TACSObject.h"
#include "TACSElement.h"

/*
  Base class for the TACSFunctionCtx. Each context is function-specific
  and is designed to store information required to 

  It's implementation is designed to be opaque to the user, but its
  data is required when evaluating the function. It is used to store
  information for each thread in the function/gradient evaluation.
*/
class TACSFunctionCtx {
 public:
  TACSFunctionCtx(){}
  virtual ~TACSFunctionCtx(){}
};

/*
  TACSFunction is the base class used to calculate the values of
  functions of interest within TACS. This class also defines the
  methods required for gradient evaluation. This class should be used
  for objectives and constraints for optimization problems posed using
  TACS.

  Thread-safe implementation
  --------------------------

  This base class is designed to support thread-safe execution of
  certain performance-critical methods. The main issue is that the
  function class must accumulate/integrate data that will be used to
  evaluate the final function value within a local member of
  TACSFunction. This is not thread-safe since multiple threads could
  update the same memory address.

  To circumvent this issue, the code uses a function context for each
  active thread. Each context is thread-specific and data is
  accumulated first within each context then a mutex is used to ensure
  thread-safe writing to the data within the function class itself.

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
  enum DomainType { ENTIRE_DOMAIN, SUB_DOMAIN, NO_DOMAIN };
  enum StageType { SINGLE_STAGE, TWO_STAGE };
  enum EvaluationType { INITIALIZE, INTEGRATE };

  TACSFunction( TACSAssembler *_tacs, 
                DomainType _funcDomain=ENTIRE_DOMAIN,
                StageType _funcStages=SINGLE_STAGE,
                int _maxElems=0 );
  virtual ~TACSFunction();

  virtual const char *functionName() = 0;
  const char *TACSObjectName();

  // Functions for setting/adjusting the domain
  // ------------------------------------------
  DomainType getDomainType();
  StageType getStageType();

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

  // Create the function context for evaluation
  // ------------------------------------------
  virtual TACSFunctionCtx *createFunctionCtx() = 0;

  // Collective calls on the TACS MPI Comm
  // -------------------------------------
  virtual void initEvaluation( EvaluationType ftype ){}
  virtual void finalEvaluation( EvaluationType ftype ){}

  // Functions for integration over the structural domain on each thread
  // -------------------------------------------------------------------
  virtual void initThread( double tcoef,
                           EvaluationType ftype,
                           TACSFunctionCtx *ctx ){}
  virtual void elementWiseEval( EvaluationType ftype,
                                TACSElement *element, int elemNum,
                                const TacsScalar Xpts[], const TacsScalar vars[],
                                const TacsScalar dvars[], const TacsScalar ddvars[],
                                TACSFunctionCtx *ctx ){}
  virtual void finalThread( double tcoef, 
                            EvaluationType ftype,
                            TACSFunctionCtx *ctx ){}

  // Return the value of the function
  // --------------------------------
  virtual TacsScalar getFunctionValue() = 0;

  // State variable sensitivities
  // ----------------------------
  virtual void getElementSVSens( double alpha, double beta, double gamma, 
                                 TacsScalar *elemSVSens, 
                                 TACSElement *element, int elemNum,
                                 const TacsScalar Xpts[], const TacsScalar vars[],
                                 const TacsScalar dvars[], const TacsScalar ddvars[],
                                 TACSFunctionCtx *ctx ){
    int numVars = element->numVariables();
    memset(elemSVSens, 0, numVars*sizeof(TacsScalar));
  }

  // Design variable sensitivity evaluation
  // --------------------------------------
  virtual void addElementDVSens( double tcoef, TacsScalar *fdvSens, int numDVs,
                                 TACSElement *element, int elemNum,
                                 const TacsScalar Xpts[], const TacsScalar vars[],
                                 const TacsScalar dvars[], const TacsScalar ddvars[],
                                 TACSFunctionCtx *ctx ){}

  // Nodal sensitivities
  // -------------------
  virtual void getElementXptSens( double tcoef, TacsScalar fXptSens[],
                                  TACSElement *element, int elemNum,
                                  const TacsScalar Xpts[], const TacsScalar vars[],
                                  const TacsScalar dvars[], const TacsScalar ddvars[],
                                  TACSFunctionCtx *ctx ){
    int numNodes = element->numNodes();
    memset(fXptSens, 0, 3*numNodes*sizeof(TacsScalar));
  }

 protected:
  TACSAssembler *tacs;
  
 private:
  // Store the function domain type
  DomainType funcDomain; 
  StageType funcStageType;

  // Store the element domain information
  int maxElems; // maximum size of currently allocated elemNums array
  int numElems; // number of elements actually stored in elemNums
  int *elemNums; // sorted array of element numbers
};

#endif
