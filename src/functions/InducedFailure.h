/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#ifndef TACS_INDUCED_FAILURE_H
#define TACS_INDUCED_FAILURE_H

/*
  Compute an aggregated function using an induced norm approach
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

class TACSInducedFailure : public TACSFunction {
 public:
  enum InducedNormType { EXPONENTIAL, POWER, 
                         EXPONENTIAL_SQUARED, POWER_SQUARED,
                         DISCRETE_EXPONENTIAL, 
                         DISCRETE_POWER, 
                         DISCRETE_EXPONENTIAL_SQUARED,
                         DISCRETE_POWER_SQUARED };
  enum InducedConstitutiveFunction { FAILURE, BUCKLING };

  TACSInducedFailure( TACSAssembler *_tacs, double _P,
                      InducedConstitutiveFunction func=FAILURE );
  ~TACSInducedFailure();

  // Retrieve the name of the function
  // ---------------------------------
  const char * functionName();

  // Set parameters to control how the induced functions are evaluated
  // -----------------------------------------------------------------
  void setParameter( double _P );
  double getParameter();
  void setInducedType( InducedNormType _norm_type );
  void setLoadFactor( TacsScalar _loadFactor );

  // Set the value of the failure offset for numerical stability
  // -----------------------------------------------------------
  void setMaxFailOffset( TacsScalar _max_fail ){
    max_fail = _max_fail;
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
  // The type of norm to evaluate
  InducedNormType norm_type;
  InducedConstitutiveFunction con_type;

  TacsScalar load_factor; // Load factor applied to the strain
  int max_nodes, max_stresses; // The max number of nodes/stresses

  TacsScalar max_fail; // The maximum failure function at a Gauss point
  TacsScalar fail_numer, fail_denom; // The numerator and denominator
  
  // The P in the P-norm
  double P;

  // The name of the function
  static const char *funcName;
};

#endif // TACS_INDUCED_FAILURE_H
