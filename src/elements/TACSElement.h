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

#ifndef TACS_ELEMENT_H
#define TACS_ELEMENT_H

/*
  Basic TACSElement definition

  The purpose of this file is to provide an interface for creating and
  storing different instances of the finite elements that will be used
  by TACS.  This is what should be extended when including more
  elements and not the underlying TACS implementation itself.
*/

#include "TACSObject.h"
#include "TACSElementTypes.h"

// The TACSElement base class
class TACSElement : public TACSOptObject {
 public:
  TACSElement( int _componentNum=0 ){
    componentNum = _componentNum;
  }
  virtual ~TACSElement(){}

  // Get the number of displacements, stresses, nodes, etc.
  // ------------------------------------------------------
  virtual int getNumDisplacements() = 0; // Degrees of freedom per node
  virtual int getNumNodes() = 0; // Number of nodes for this element

  // Identifies whether the nodes are associated with the multipliers
  virtual void getMultiplierIndex( int *multiplier ){
    *multiplier = -1;
  }

  // Get the number of extras and element type information
  // -----------------------------------------------------
  virtual ElementLayout getElementLayout(){ return TACS_LAYOUT_NONE; }

  // Retrieve the initial conditions and add the derivative
  // ------------------------------------------------------
  virtual void getInitConditions( TacsScalar vars[],
                                  TacsScalar dvars[],
                                  TacsScalar ddvars[],
                                  const TacsScalar Xpts[] ){
    memset(vars, 0, numVariables()*sizeof(TacsScalar));
    memset(dvars, 0, numVariables()*sizeof(TacsScalar));
    memset(ddvars, 0, numVariables()*sizeof(TacsScalar));
  }

  // Add the product of the initial condition with the adjoint variables
  // -------------------------------------------------------------------
  virtual void addInitConditionAdjResProduct( TacsScalar fdvSens[], int dvLen,
                                              const TacsScalar adjVars[],
                                              const TacsScalar adjDVars[],
                                              const TacsScalar Xpts[] ){}
  virtual void getInitConditionAdjResXptProduct( TacsScalar fXptSens[],
                                                 const TacsScalar adjVars[],
                                                 const TacsScalar adjDVars[],
                                                 const TacsScalar Xpts[] ){
    memset(fXptSens, 0, 3*numNodes()*sizeof(TacsScalar));
  }

  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  virtual void computeEnergies( double time,
                                TacsScalar *_Te,
                                TacsScalar *_Pe,
                                const TacsScalar Xpts[],
                                const TacsScalar vars[],
                                const TacsScalar dvars[] ){
    *_Te = 0.0;
    *_Pe = 0.0;
  }

  // Compute the residual of the governing equations
  // -----------------------------------------------
  virtual void addResidual( double time, TacsScalar res[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] ) = 0;

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  virtual void addJacobian( double time, TacsScalar res[], TacsScalar J[],
                            double alpha, double beta, double gamma,
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] );

  // Add the product of the adjoint variables with the derivative of the residual
  // ----------------------------------------------------------------------------
  virtual void addAdjResProduct( double time, double scale,
                                 TacsScalar dvSens[], int dvLen,
                                 const TacsScalar psi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 const TacsScalar dvars[],
                                 const TacsScalar ddvars[] ){}
  virtual void addAdjResXptProduct( double time, double scale,
                                    TacsScalar fXptSens[],
                                    const TacsScalar psi[],
                                    const TacsScalar Xpts[],
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[] ){}

  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  virtual void getMatType( ElementMatrixType matType,
                           TacsScalar mat[],
                           const TacsScalar Xpts[],
                           const TacsScalar vars[] ){
    int size = numVariables()*numVariables();
    memset(mat, 0, size*sizeof(TacsScalar));
  }

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  virtual void addMatDVSensInnerProduct( ElementMatrixType matType,
                                         double scale,
                                         TacsScalar dvSens[], int dvLen,
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[] ){}

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  virtual void getMatSVSensInnerProduct( ElementMatrixType matType,
                                         TacsScalar res[],
                                         const TacsScalar psi[],
                                         const TacsScalar phi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[] ){
    memset(res, 0, numVariables()*sizeof(TacsScalar));
  }

  // Functions for retrieving data from the element for visualization
  // ----------------------------------------------------------------
  void setComponentNum( int comp_num ){ componentNum = comp_num; }
  int getComponentNum(){ return componentNum; }
  virtual void getOutputData( ElementType etype, int write_flag,
                              double *data, int ld_data,
                              const TacsScalar Xpts[],
                              const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[] ){}

 private:
  int componentNum;
};

#endif // TACS_ELEMENT_H
