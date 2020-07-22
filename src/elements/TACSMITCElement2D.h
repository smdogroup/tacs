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

#ifndef TACS_ELEMENT_2D_H
#define TACS_ELEMENT_2D_H

#include "TACSElement.h"
#include "TACSMITCModel.h"
#include "TACSMITCBasis.h"

class TACSMITCElement2D : public TACSElement {
 public:
  static const int MAX_VARS_PER_NODE = 8;
  static const int MAX_NUM_TYING_FIELDS = 5;
  static const int MAX_TOTAL_TYING_POINTS = 5*16;

  TACSMITCElement2D( TACSMITCModel *_model, TACSMITCBasis *_basis );
  ~TACSMITCElement2D();

  // Get the layout properties of the element
  int getVarsPerNode();
  int getNumNodes();
  int getDesignVarsPerNode();
  ElementLayout getLayoutType();
  TACSElementBasis* getElementBasis();

  /**
    Retrieve the global design variable numbers associated with this element
  */
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  /**
    Set the element design variables from the design vector
  */
  int setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  /**
    Get the element design variables values
  */
  int getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  /**
    Get the lower and upper bounds for the design variable values
  */
  int getDesignVarRange( int elemIndex, int dvLen,
                         TacsScalar lb[], TacsScalar ub[] );

  /**
    Add the residual to the provided vector
  */
  void addResidual( int elemIndex, double time, const TacsScalar *Xpts,
                    const TacsScalar *vars, const TacsScalar *dvars,
                    const TacsScalar *ddvars, TacsScalar *res );

  /**
    Evaluate a point-wise quantity of interest.
  */
  int evalPointQuantity( int elemIndex, int quantityType, double time,
                         int n, double pt[], const TacsScalar Xpts[],
                         const TacsScalar vars[], const TacsScalar dvars[],
                         const TacsScalar ddvars[], TacsScalar *quantity );

  /**
    Add the derivative of the point quantity w.r.t. the design variables
  */
  void addPointQuantityDVSens( int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[], const TacsScalar ddvars[],
                               const TacsScalar dfdq[], int dvLen, TacsScalar dfdx[] );

  /**
    Add the derivative of the point quantity w.r.t. the state variables
  */
  void addPointQuantitySVSens( int elemIndex, int quantityType, double time,
                               TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                               int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[], const TacsScalar ddvars[],
                               const TacsScalar dfdq[], TacsScalar dfdu[] );

  /**
    Add the derivative of the point quantity w.r.t. the node locations
  */
  void addPointQuantityXptSens( int elemIndex, int quantityType, double time,
                                TacsScalar scale, int n, double pt[],
                                const TacsScalar Xpts[], const TacsScalar vars[],
                                const TacsScalar dvars[], const TacsScalar ddvars[],
                                const TacsScalar dfdq[], TacsScalar dfdXpts[] );

  /**
    Compute the output data for visualization
  */
  void getOutputData( int elemIndex, ElementType etype, int write_flag,
                      const TacsScalar Xpts[], const TacsScalar vars[],
                      const TacsScalar dvars[], const TacsScalar ddvars[],
                      int ld_data, TacsScalar *data );

 private:
  TACSMITCModel *model;
  TACSMITCBasis *basis;
};

#endif // TACS_ELEMENT_2D_H
