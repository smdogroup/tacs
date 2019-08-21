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

#ifndef TACS_ELEMENT_3D_H
#define TACS_ELEMENT_3D_H

#include "TACSElementModel.h"
#include "TACSElementBasis.h"

template <int VARS_PER_NODE>
class TACSElement3D : public TACSElement {
 public:
  TACSElement3D( TACSElementModel *_model, TACSElementBasis *_basis ){
    model = _model;
    basis = _basis;
  }
  
  int getNumDisplacements(){ return VARS_PER_NODE; }
  int getNumNodes(){ return basis->getNumNodes(); }

  // Get the variable information
  // ----------------------------
  ElementLayout getLayoutType(){ basis->getLayoutType(); }
  
  // The design variable query functions
  // -----------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lowerBound[], 
                          TacsScalar upperBound[], int numDVs );
  
  // Compute the kinetic and potential energy within the element
  // -----------------------------------------------------------
  void computeEnergies( double time, TacsScalar *_Te, TacsScalar *_Pe,
                        const TacsScalar Xpts[], const TacsScalar vars[],
                        const TacsScalar dvars[] );

  // Compute the residual of the governing equations
  // -----------------------------------------------
  void addResidual( double time, TacsScalar res[], const TacsScalar Xpts[],
                    const TacsScalar vars[], const TacsScalar dvars[],
                    const TacsScalar ddvars[] );

  // Compute the Jacobian of the governing equations
  // -----------------------------------------------
  void addJacobian( double time, TacsScalar res[], TacsScalar J[],
                    double alpha, double beta, double gamma,
                    const TacsScalar Xpts[], const TacsScalar vars[],
                    const TacsScalar dvars[], const TacsScalar ddvars[] );

  // Add the product of the adjoint with the derivative of the design variables
  // --------------------------------------------------------------------------
  void addAdjResProduct( double time, double scale,
                         TacsScalar dvSens[], int dvLen,
                         const TacsScalar psi[], const TacsScalar Xpts[],
                         const TacsScalar vars[], const TacsScalar dvars[],
                         const TacsScalar ddvars[] );

  // Add the product of the adjoint with the derivative of the design variables
  // -------------------------------------------------------------------------- 
  void addAdjResXptProduct( double time, double scale, 
                            TacsScalar fXptSens[],
                            const TacsScalar psi[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[] );
    
  // Retrieve a specific time-independent matrix from the element
  // ------------------------------------------------------------
  void getMatType( ElementMatrixType matType, 
                   TacsScalar mat[], 
                   const TacsScalar Xpts[],
                   const TacsScalar vars[] );

  // Compute the derivative of the inner product w.r.t. design variables
  // -------------------------------------------------------------------
  void addMatDVSensInnerProduct( ElementMatrixType matType, 
                                 double scale,
                                 TacsScalar dvSens[], int dvLen,
                                 const TacsScalar psi[], 
                                 const TacsScalar phi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[] );

  // Compute the derivative of the inner product w.r.t. vars[]
  // ---------------------------------------------------------
  void getMatSVSensInnerProduct( ElementMatrixType matType, 
                                 TacsScalar res[],
                                 const TacsScalar psi[], 
                                 const TacsScalar phi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[] );

};

/*
  Add the residual to the provided vector

  @param time The simulation time 
  @param Xpts The element node locations
  @param vars The element state variable values
  @param dvars The time derivative of the element state variables
  @param ddvars The second time derivative of the element state variables
  @param res The element residual output 
*/
template <int VARS_PER_NODE>
void TACSElement3D<VARS_PER_NODE>::addResidual( double time,
                                                const TacsScalar *Xpts,
                                                const TacsScalar *vars,
                                                const TacsScalar *dvars,
                                                const TacsScalar *ddvars, 
                                                TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double weight = basis->getQuadratureWeight(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], J[9];
    TacsScalar U[VARS_PER_NODE], Udot[VARS_PER_NODE], Uddot[VARS_PER_NODE];
    TacsScalar Ux[3*VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(n, Xpts, VARS_PER_NODE,
                                              vars, dvars, ddvars,
                                              X, J, U, Udot, Uddot, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[VARS_PER_NODE*3], DUx[6*VARS_PER_NODE];
    model->evalWeakIntegrand(time, pt, X, U, Udot, Uddot, DUt, DUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormResidual(n, pt, detJ, J, VARS_PER_NODE, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays

  @param time The simulation time 
  @param Xpts The element node locations
  @param vars The element state variable values
  @param dvars The time derivative of the element state variables
  @param ddvars The second time derivative of the element state variables
  @param res The element residual output
  @param mat The element Jacobian output 
*/
template <int VARS_PER_NODE>
void TACSElement3D<VARS_PER_NODE>::addJacobian( double time,
                                                double alpha,
                                                double beta,
                                                double gamma,
                                                const TacsScalar *Xpts,
                                                const TacsScalar *vars,
                                                const TacsScalar *dvars,
                                                const TacsScalar *ddvars, 
                                                TacsScalar *res,
                                                TacsScalar *mat ){
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double weight = basis->getQuadratureWeight(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], J[9];
    TacsScalar U[VARS_PER_NODE], Udot[VARS_PER_NODE], Uddot[VARS_PER_NODE];
    TacsScalar Ux[3*VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(n, Xpts, VARS_PER_NODE,
                                              vars, dvars, ddvars,
                                              X, J, U, Udot, Uddot, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    int DDUt_nnz, DDUx_nnz;
    const int *DDUt_pairs, *DDUx_pairs; 
    TacsScalar DUt[VARS_PER_NODE*3], DUx[6*VARS_PER_NODE];
    TacsScalar DDUt[9*VARS_PER_NODE*VARS_PER_NODE];
    TacsScalar DDUx[32*VARS_PER_NODE*VARS_PER_NODE];
    model->evalWeakIntegrandJacobian(time, pt, X, U, Udot, Uddot, DUt, DUx,
                                     &DDUt_nnz, &DDUt_pairs, DDUt,
                                     &DDUx_nnz, &DDUx_pairs, DDUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormJacobian(n, pt, detJ, J, VARS_PER_NODE, DUt, DUx,
                               alpha, beta, gamma,
                               DDUt_nnz, DDUt_pairs, DDUt,
                               DDUx_nnz, DDUx_paris, DDUx, res, mat);
  }
}

#endif // TACS_3D_ELEMENT_H
