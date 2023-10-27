/*
=============================================================================
Generic aggregate function class for TACS
=============================================================================
@File    :   TACSAggregateFunction.h
@Date    :   2023/10/26
@Author  :   Alasdair Christison Gray
@Description :
*/

#pragma once

// =============================================================================
// Standard Library Includes
// =============================================================================

// =============================================================================
// Extension Includes
// =============================================================================
#include "TACSFunction.h"

/**
 * @brief
 *
 */
class TACSAggregateFunction : public TACSFunction {
 public:
  enum AggregationType {
    DISCRETE,
    CONTINUOUS,
    PNORM_DISCRETE,
    PNORM_CONTINUOUS
  };

  TACSAggregateFunction(TACSAssembler *_assembler, int outputType,
                        double _aggWeight, double _alpha = 1.0,
                        double _funcScaleFactor = 1.0);
  ~TACSAggregateFunction();

  /**
    Get the object/function name
  */
  const char *getObjectName();

  // Set parameters for the aggregation function
  // ----------------------------------
  void setAggregationType(enum AggregationType type);
  double getParameter();
  void setParameter(double _aggWeight);

  // Set the value of the failure offset for numerical stability
  // -----------------------------------------------------------
  void setMaxValOffset(TacsScalar _maxVal) { maxVal = _maxVal; }

  /*
    Retrieve the maximum value
  */
  TacsScalar getMaximumValue() { return maxVal; }

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
                        TacsScalar *elemSVSens);

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
  // The type of output being aggregated
  int outputType;

  // The type of aggregation to use
  AggregationType aggType;

  // The aggregation weight parameter
  double aggWeight;

  // The integral scaling value
  double alpha;

  // The safety factor
  double funcScaleFactor;

  // The name of the function
  static const char *funcName;

  // The maximum failure value, the sum of exp(aggWeight*(f[i] - maxVal)
  // and the value of the KS function
  TacsScalar aggSum, maxVal;

  // Used for the case when this is used to evaluate the p-norm
  TacsScalar invPnorm;

  void addFuncVal(EvaluationType ftype, TacsScalar funcVal, TacsScalar scale);

  TacsScalar computeFuncSens(TacsScalar funcVal);
};
