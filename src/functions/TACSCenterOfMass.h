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

#ifndef TACS_CENTER_OF_MASS_H
#define TACS_CENTER_OF_MASS_H

/*
  Compute the structural mass
*/

#include "TACSFunction.h"

/*
  Evaluate the structural mass of the structure
*/
class TACSCenterOfMass : public TACSFunction {
 public:
  TACSCenterOfMass(TACSAssembler *_assembler, const double _dir[]);
  ~TACSCenterOfMass();

  const char *getObjectName();

  /**
    Member functions to integrate the function value
  */
  void initEvaluation(EvaluationType ftype);
  void elementWiseEval(EvaluationType ftype, int elemIndex,
                       TACSElement *element, double time, TacsScalar scale,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[], const TacsScalar ddvars[]);
  void finalEvaluation(EvaluationType ftype);

  /**
    Return the value of the function
  */
  TacsScalar getFunctionValue();

  /**
     Add the derivative of the function w.r.t. the design variables
  */
  void addElementDVSens(int elemIndex, TACSElement *element, double time,
                        TacsScalar scale, const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dfdx[]);

  /**
     Evaluate the derivative of the function w.r.t. the node locations
  */
  void getElementXptSens(int elemIndex, TACSElement *element, double time,
                         TacsScalar scale, const TacsScalar Xpts[],
                         const TacsScalar vars[], const TacsScalar dvars[],
                         const TacsScalar ddvars[], TacsScalar fXptSens[]);

 private:
  // The total mass of all elements in the specified domain
  TacsScalar totalMass, massMoment;
  // cg projection direction
  double dir[3];

  static const char *funcName;
};

#endif  // TACS_CENTER_OF_MASS_H
