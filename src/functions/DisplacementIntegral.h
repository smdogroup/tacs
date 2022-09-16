/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2018 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_DISPLACEMENT_INTEGRAL_H
#define TACS_DISPLACEMENT_INTEGRAL_H

#include "TACSFunction.h"

/*
  Compute the KS functional of the displacement along a given direction
*/
class TACSDisplacementIntegral : public TACSFunction {
 public:
  TACSDisplacementIntegral(TACSAssembler *_tacs, const TacsScalar _dir[]);
  ~TACSDisplacementIntegral();

  // Retrieve the name of the function
  // ---------------------------------
  const char *functionName();

  // Create the function context for evaluation
  // ------------------------------------------
  TACSFunctionCtx *createFunctionCtx();

  // Collective calls on the TACS MPI Comm
  // -------------------------------------
  void initEvaluation(EvaluationType ftype);
  void finalEvaluation(EvaluationType ftype);

  // Functions for integration over the structural domain on each thread
  // -------------------------------------------------------------------
  void initThread(double tcoef, EvaluationType ftype, TACSFunctionCtx *ctx);
  void elementWiseEval(EvaluationType ftype, TACSElement *element, int elemNum,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[], const TacsScalar ddvars[],
                       TACSFunctionCtx *ctx);
  void finalThread(double tcoef, EvaluationType ftype, TACSFunctionCtx *ctx);

  // Return the value of the function
  // --------------------------------
  TacsScalar getFunctionValue();

  // State variable sensitivities
  // ----------------------------
  void getElementSVSens(double alpha, double beta, double gamma,
                        TacsScalar *elemSVSens, TACSElement *element,
                        int elemNum, const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TACSFunctionCtx *ctx);

  // Design variable sensitivity evaluation
  // --------------------------------------
  void addElementDVSens(double tcoef, TacsScalar *fdvSens, int numDVs,
                        TACSElement *element, int elemNum,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TACSFunctionCtx *ctx);

  // Nodal sensitivities
  // -------------------
  void getElementXptSens(double tcoef, TacsScalar fXptSens[],
                         TACSElement *element, int elemNum,
                         const TacsScalar Xpts[], const TacsScalar vars[],
                         const TacsScalar dvars[], const TacsScalar ddvars[],
                         TACSFunctionCtx *ctx);

 private:
  // The name of the function
  static const char *funcName;

  // The direction
  TacsScalar dir[3];

  // The value of the KS weight
  TacsScalar value;

  // The max number of nodes
  int maxNumNodes;
};

#endif  // TACS_DISPLACEMENT_INTEGRAL_H
