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

#ifndef TACS_KS_DISPLACEMENT_H
#define TACS_KS_DISPLACEMENT_H

/*
  Compute the KS function in TACS
*/

#include "TACSFunction.h"

/*
  The following class implements the methods from TACSFunction.h
  necessary to calculate the KS function of either a stress or strain
  failure criteria over the domain of some finite element model.

  Each class should only ever be passed to a single instance of
  TACS. If the KS function needs to be calculated for separate
  instances, this should be handled by separate instances of
  KSDispure.

  The failure load is calculated using the strain-based failure
  criteria from the base Constitutive class which requires linear and
  constant components of the strain to determine the failure load.

  The arguments to the KSDispure class are:

  ksWeight:  the ks weight used in the calculation

  optional arguments:

  elementNums, numElements -- these specify a subdomain of the TACS
  model over which the KS function should be calculated by passing in
  the element numbers and number of elements in the subdomain.

  note: if no subdomain is specified, the calculation takes place over
  all the elements in the model
*/
class TACSKSDisplacement : public TACSFunction {
 public:
  TACSKSDisplacement( TACSAssembler *_assembler, double ksWeight,
                      double alpha=1.0 );
  ~TACSKSDisplacement();

  // Retrieve the name of the function
  // ---------------------------------
  const char* getObjectName();

  // Set parameters for the KS function
  // ----------------------------------
  double getParameter();
  void setParameter( double _ksWeight );

  /**
     Initialize the function for the given type of evaluation
  */
  void initEvaluation( EvaluationType ftype );

  /**
     Perform an element-wise integration over this element.
  */
  void elementWiseEval( EvaluationType ftype,
                        int elemIndex, TACSElement *element,
                        double time, TacsScalar scale,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[] );

  /**
     Finalize the function evaluation for the specified eval type.
  */
  void finalEvaluation( EvaluationType ftype );

  /**
     Get the value of the function
  */
  TacsScalar getFunctionValue();

  /**
     Evaluate the derivative of the function w.r.t. state variables
  */
  void getElementSVSens( int elemIndex, TACSElement *element, double time,
                         TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TacsScalar *elemSVSens );

  /**
     Add the derivative of the function w.r.t. the design variables
  */
  void addElementDVSens( int elemIndex, TACSElement *element,
                         double time, TacsScalar scale,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         int dvLen, TacsScalar dfdx[] );

  /**
     Evaluate the derivative of the function w.r.t. the node locations
  */
  void getElementXptSens( int elemIndex, TACSElement *element,
                          double time, TacsScalar scale,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TacsScalar fXptSens[] );

 private:
  // The weight on the ks function value
  double ksWeight;

  // The integral scaling value
  double alpha;

  // The name of the function
  static const char *funcName;

  // The maximum failure value, the sum of exp(ksWeight*(f[i] - maxDisp)
  // and the value of the KS function
  TacsScalar ksDispSum, maxDisp;

  // Used for the case when this is used to evaluate the p-norm
  TacsScalar invPnorm;
};

#endif // TACS_KS_FAILURE_H
