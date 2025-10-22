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

#ifndef TACS_AVERAGE_TEMPERATURE_H
#define TACS_AVERAGE_TEMPERATURE_H

/*
  Calculate the average temperature
*/
#include "TACSFunction.h"

/*
  Evaluate the spatial average of the temperature of the model.
*/
class TACSAverageTemperature : public TACSFunction {
 public:
  TACSAverageTemperature(TACSAssembler *_assembler, TacsScalar _volume = 1.0);
  ~TACSAverageTemperature();

  /**
    Get the object/function name
  */
  const char *getObjectName();

  /**
     Initialize the function for the given type of evaluation
  */
  void initEvaluation(EvaluationType ftype);

  /**
     Perform an element-wise integration over this element.
  */
  void elementWiseEval(EvaluationType ftype, int elemIndex,
                       TACSElement *element, double time, TacsScalar scale,
                       const TacsScalar Xpts[], const TacsScalar vars[],
                       const TacsScalar dvars[], const TacsScalar ddvars[]);

  /**
     Finalize the function evaluation for the specified eval type.
  */
  void finalEvaluation(EvaluationType ftype);

  /**
     Get the value of the function
  */
  TacsScalar getFunctionValue();

  /**
     Evaluate the derivative of the function w.r.t. state variables
  */
  void getElementSVSens(int elemIndex, TACSElement *element, double time,
                        TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[], const TacsScalar ddvars[],
                        TacsScalar dfdu[]);

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
  // The name of the function
  static const char *funcName;

  // The volume (area) and the integral of temperature
  TacsScalar inv_volume;
  TacsScalar integral_temp;
};

#endif  // TACS_AVERAGE_TEMPERATURE_H
