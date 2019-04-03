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

#ifndef TACS_KS_FAILURE_H
#define TACS_KS_FAILURE_H

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
  KSFailure.

  The failure load is calculated using the strain-based failure
  criteria from the base Constitutive class which requires linear and
  constant components of the strain to determine the failure load.

  The arguments to the KSFailure class are:

  ksWeight:  the ks weight used in the calculation

  optional arguments:

  elementNums, numElements -- these specify a subdomain of the TACS
  model over which the KS function should be calculated by passing in
  the element numbers and number of elements in the subdomain.

  note: if no subdomain is specified, the calculation takes place over
  all the elements in the model
*/
class TACSKSFailure : public TACSFunction {
 public:
  enum KSFailureType { DISCRETE, CONTINUOUS,
                       PNORM_DISCRETE, PNORM_CONTINUOUS };
  enum KSConstitutiveFunction { FAILURE, BUCKLING };

  TACSKSFailure( TACSAssembler * _tacs, double ksWeight,
                 KSConstitutiveFunction func=FAILURE,
                 double alpha=1.0 );
  ~TACSKSFailure();

  // Retrieve the name of the function
  // ---------------------------------
  const char *functionName();

  // Set parameters for the KS function
  // ----------------------------------
  void setKSFailureType( enum KSFailureType type );
  double getParameter();
  void setParameter( double _ksWeight );
  void setLoadFactor( TacsScalar _loadFactor );

  // Set the value of the failure offset for numerical stability
  // -----------------------------------------------------------
  void setMaxFailOffset( TacsScalar _maxFail ){
    maxFail = _maxFail;
  }

  // Create the function context for evaluation
  // ------------------------------------------
  TACSFunctionCtx *createFunctionCtx();

  // Collective calls on the TACS MPI Comm
  // -------------------------------------
  void initEvaluation( EvaluationType ftype );
  void finalEvaluation( EvaluationType ftype );

  // Functions for integration over the structural domain on each thread
  // -------------------------------------------------------------------
  void initThread( double tcoef,
                   EvaluationType ftype,
                   TACSFunctionCtx *ctx );
  void elementWiseEval( EvaluationType ftype,
                        TACSElement *element, int elemNum,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TACSFunctionCtx *ctx );
  void finalThread( double tcoef,
                    EvaluationType ftype,
                    TACSFunctionCtx *ctx );

  // Return the value of the function
  // --------------------------------
  TacsScalar getFunctionValue();
  TacsScalar getMaximumFailure();

  // State variable sensitivities
  // ----------------------------
  void getElementSVSens( double alpha, double beta, double gamma,
                         TacsScalar *elemSVSens,
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx );

  // Design variable sensitivity evaluation
  // --------------------------------------
  void addElementDVSens( double tcoef, TacsScalar *fdvSens, int numDVs,
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx );

  // Nodal sensitivities
  // -------------------
  void getElementXptSens( double tcoef, TacsScalar fXptSens[],
                          TACSElement *element, int elemNum,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TACSFunctionCtx *ctx );
 private:
  // The type of aggregation to use
  KSFailureType ksType;

  // The constitutive function to use
  KSConstitutiveFunction conType;

  // The weight on the ks function value
  double ksWeight;

  // The integral scaling value
  double alpha;

  // Load factor applied to the strain
  TacsScalar loadFactor;

  // The maximum number of nodes/stresses in any given element
  int maxNumNodes, maxNumStrains;

  // The name of the function
  static const char *funcName;

  // The maximum failure value, the sum of exp(ksWeight*(f[i] - maxFail)
  // and the value of the KS function
  TacsScalar ksFailSum, maxFail;

  // Used for the case when this is used to evaluate the p-norm
  TacsScalar invPnorm;
};

#endif // TACS_KS_FAILURE_H
