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

#include "TACSElement.h"
#include "TACSObject.h"

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

  TACSFunction(TACSAssembler *_assembler,
               DomainType _funcDomain = ENTIRE_DOMAIN,
               StageType _funcStages = SINGLE_STAGE, int _maxElems = 0);
  virtual ~TACSFunction();

  /**
     Get the object name
  */
  const char *getObjectName();

  /**
     Get the type of integration domain

     @return The enum type of domain
  */
  DomainType getDomainType();

  /**
     Get the stage type of this function: Either one or two stage

     Some functions (such as aggregation functionals) require a
     two-stage integration strategy for numerical stability.

     @return The enum type indicating whether this is a one or two stage func.
  */
  StageType getStageType();

  /**
     Set the element numbers within the domain to integrate

     @param numElems The number of elements in the array
     @param elemNums The local element numbers for the domain
  */
  void setDomain(int numElems, const int _elemNums[]);

  /**
     Add the element numbers to the domain

     @param numElems The number of elements in the array
     @param elemNums The local element numbers for the domain
  */
  void addDomain(int numElems, const int elemNums[]);

  /**
     Retrieve the element domain from the function

     @param elemNums The element numbers defining the domain
     @return The numer of elements in the domain
  */
  int getElementNums(const int **_elemNums);

  /**
     Return the TACSAssembler object associated with this function
  */
  TACSAssembler *getAssembler();

  /**
     Initialize the function for the given type of evaluation

     This call is collective on all processors in the assembler.
  */
  virtual void initEvaluation(EvaluationType ftype) {}

  /**
     Perform an element-wise integration over this element.

     Note that this is not a collective call and should be called once
     for each element within the integration domain.

     @param ftype The type of evaluation
     @param elemIndex The local element index
     @param element The TACSElement object
     @param time The simulation time
     @param scale The scalar integration factor to apply
     @param Xpts The element node locations
     @param vars The element DOF
     @param dvars The first time derivatives of the element DOF
     @param ddvars The second time derivatives of the element DOF
  */
  virtual void elementWiseEval(EvaluationType ftype, int elemIndex,
                               TACSElement *element, double time,
                               TacsScalar scale, const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[]) {}

  /**
     Finalize the function evaluation for the specified eval type.

     This call is collective on all processors in the assembler.
  */
  virtual void finalEvaluation(EvaluationType ftype) {}

  /**
     Get the value of the function
  */
  virtual TacsScalar getFunctionValue() = 0;

  /**
     Evaluate the derivative of the function w.r.t. state variables

     @param elemIndex The local element index
     @param element The TACSElement object
     @param time The simulation time
     @param alpha Coefficient for the DOF derivative
     @param beta Coefficient for the first time DOF derivative
     @param gamma Coefficient for the second time DOF derivative
     @param Xpts The element node locations
     @param vars The element DOF
     @param dvars The first time derivatives of the element DOF
     @param ddvars The second time derivatives of the element DOF
  */
  virtual void getElementSVSens(int elemIndex, TACSElement *element,
                                double time, TacsScalar alpha, TacsScalar beta,
                                TacsScalar gamma, const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[],
                                const TacsScalar ddvars[], TacsScalar dfdu[]) {
    int numVars = element->getNumVariables();
    memset(dfdu, 0, numVars * sizeof(TacsScalar));
  }

  /**
     Add the derivative of the function w.r.t. the design variables

     The design variables *must* be the same set of variables defined
     in the element. The TACSFunction class cannot define new design
     variables!

     @param elemIndex The local element index
     @param element The TACSElement object
     @param time The simulation time
     @param Xpts The element node locations
     @param vars The element DOF
     @param dvars The first time derivatives of the element DOF
     @param ddvars The second time derivatives of the element DOF
  */
  virtual void addElementDVSens(int elemIndex, TACSElement *element,
                                double time, TacsScalar scale,
                                const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[],
                                const TacsScalar ddvars[], int dvLen,
                                TacsScalar dfdx[]) {}

  /**
     Evaluate the derivative of the function w.r.t. the node locations

     @param elemIndex The local element index
     @param element The TACSElement object
     @param time The simulation time
     @param scale The scalar integration factor to apply
     @param Xpts The element node locations
     @param vars The element DOF
     @param dvars The first time derivatives of the element DOF
     @param ddvars The second time derivatives of the element DOF
  */
  virtual void getElementXptSens(int elemIndex, TACSElement *element,
                                 double time, TacsScalar scale,
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[],
                                 TacsScalar dfdXpts[]) {
    int numNodes = element->getNumNodes();
    memset(dfdXpts, 0, 3 * numNodes * sizeof(TacsScalar));
  }

 protected:
  TACSAssembler *assembler;

 private:
  // Store the function domain type
  DomainType funcDomain;
  StageType funcStageType;

  // Store the element domain information
  int maxElems;   // maximum size of currently allocated elemNums array
  int numElems;   // number of elements actually stored in elemNums
  int *elemNums;  // sorted array of element numbers
};

#endif  // TACS_FUNCTION_H
