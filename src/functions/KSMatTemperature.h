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

#ifndef TACS_KS_MAT_TEMPERATURE_H
#define TACS_KS_MAT_TEMPERATURE_H

#include "KSTemperature.h"
#include "TACSFunction.h"
/*
  Compute the KS functional of the displacement along a given direction
*/
class TACSKSMatTemperature : public TACSFunction {
 public:
  /* enum KSTemperatureType { DISCRETE, CONTINUOUS, */
  /*                          PNORM_DISCRETE, PNORM_CONTINUOUS }; */

  TACSKSMatTemperature(TACSAssembler *_tacs, double _ksWeight,
                       TACSKSTemperature::KSTemperatureType _ksType =
                           TACSKSTemperature::CONTINUOUS,
                       int _nmats = 1);
  ~TACSKSMatTemperature();

  // Retrieve the name of the function
  // ---------------------------------
  const char *functionName();

  // Create the function context for evaluation
  // ------------------------------------------
  TACSFunctionCtx *createFunctionCtx();

  // Set the type of displacement aggregate
  // --------------------------------------
  void setKSDispType(TACSKSTemperature::KSTemperatureType _ksType);
  void setNumMats(int _nmats) { nmats = _nmats; }

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

  // The value of the KS weight
  double ksWeight;

  // Intermediate values in the functional evaluation
  TacsScalar ksSum;
  TacsScalar *maxValue;
  TacsScalar invPnorm;

  // Set the type of constraint aggregate
  TACSKSTemperature::KSTemperatureType ksType;

  // The max number of nodes
  int maxNumNodes;

  // Whether the domain is a plane stress or 3d
  int is_2d, is_3d;
  int nmats;
};

#endif  // TACS_KS_DISPLACEMENT_H
