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

#include "TACSElement2D.h"

TACSElement2D::TACSElement2D( TACSElementModel *_model,
                              TACSElementBasis *_basis ){
  model = _model;
  basis = _basis;
}

TACSElement2D::~TACSElement2D(){}

// Get the layout properties of the element
int TACSElement2D::getVarsPerNode(){
  return model->getVarsPerNode();
}

int TACSElement2D::getNumNodes(){
  return basis->getNumNodes();
}

ElementLayout TACSElement2D::getLayoutType(){
  return basis->getLayoutType();
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSElement2D::getDesignVarNums( int elemIndex, int dvLen, int dvNums[] ){
  return model->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
void TACSElement2D::setDesignVars( int elemIndex, int dvLen,
                                   const TacsScalar dvs[] ){
  model->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
void TACSElement2D::getDesignVars( int elemIndex, int dvLen,
                                   TacsScalar dvs[] ){
  model->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
void TACSElement2D::getDesignVarRange( int elemIndex, int dvLen,
                                       TacsScalar lb[], TacsScalar ub[] ){
  model->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

/*
  Add the residual to the provided vector
*/
void TACSElement2D::addResidual( int elemIndex,
                                 double time,
                                 const TacsScalar *Xpts,
                                 const TacsScalar *vars,
                                 const TacsScalar *dvars,
                                 const TacsScalar *ddvars,
                                 TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();
  const int vars_per_node = model->getVarsPerNode();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar U[MAX_VARS_PER_NODE], Udot[MAX_VARS_PER_NODE];
    TacsScalar Uddot[MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ud, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3*MAX_VARS_PER_NODE], DUx[3*MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, pt, X, U, Udot, Uddot, Ux, DUt, DUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormResidual(n, pt, detJ, J, vars_per_node, DUt, DUx, res);
  }
}

/*
  Add the residual and Jacobians to the arrays
*/
void TACSElement2D::addJacobian( int elemIndex,
                                 double time,
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
  const int vars_per_node = model->getVarsPerNode();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar U[MAX_VARS_PER_NODE], Udot[MAX_VARS_PER_NODE];
    TacsScalar Uddot[MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ud, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    int DDUt_nnz, DDUx_nnz;
    const int *DDUt_pairs, *DDUx_pairs;
    TacsScalar DUt[3*MAX_VARS_PER_NODE], DUx[3*MAX_VARS_PER_NODE];
    TacsScalar DDUt[9*MAX_VARS_PER_NODE*MAX_VARS_PER_NODE];
    TacsScalar DDUx[9*MAX_VARS_PER_NODE*MAX_VARS_PER_NODE];
    model->evalWeakJacobian(elemIndex, time, pt,
                            X, U, Udot, Uddot, Ux, DUt, DUx,
                            &DDUt_nnz, &DDUt_pairs, DDUt,
                            &DDUx_nnz, &DDUx_pairs, DDUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormJacobian(n, pt, detJ, J, vars_per_node, DUt, DUx,
                               alpha, beta, gamma,
                               DDUt_nnz, DDUt_pairs, DDUt,
                               DDUx_nnz, DDUx_pairs, DDUx, res, mat);
  }
}

// Functions for the adjoint
void TACSElement2D::addAdjResProduct( int elemIndex,
                                      double time,
                                      double scale,
                                      const TacsScalar psi[],
                                      const TacsScalar Xpts[],
                                      const TacsScalar vars[],
                                      const TacsScalar dvars[],
                                      const TacsScalar ddvars[],
                                      int dvLen,
                                      TacsScalar dvSens[] ){
  /*
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();
  const int vars_per_node = model->getVarsPerNode();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar U[MAX_VARS_PER_NODE], Udot[MAX_VARS_PER_NODE];
    TacsScalar Uddot[MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ud, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3*MAX_VARS_PER_NODE], DUx[3*MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, pt, X, U, Udot, Uddot, Ux, DUt, DUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormResidual(n, pt, detJ, J, vars_per_node, DUt, DUx, res);
  }
  */
}

void TACSElement2D::addAdjResXptProduct( int elemIndex,
                                         double time,
                                         double scale,
                                         const TacsScalar psi[],
                                         const TacsScalar Xpts[],
                                         const TacsScalar vars[],
                                         const TacsScalar dvars[],
                                         const TacsScalar ddvars[],
                                         TacsScalar fXptSens[] ){
  /*
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();
  const int vars_per_node = model->getVarsPerNode();

  // Loop over each quadrature point and add the residual contribution
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar U[MAX_VARS_PER_NODE], Udot[MAX_VARS_PER_NODE];
    TacsScalar Uddot[MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ud, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[3*MAX_VARS_PER_NODE], DUx[3*MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, time, pt, X, U, Udot, Uddot, Ux, DUt, DUx);

    // Add the weak form of the residual at this point
    basis->addWeakFormResidual(n, pt, detJ, J, vars_per_node, DUt, DUx, res);
  }
  */
}

/*
  Get the element data for the basis
*/
void TACSElement2D::getOutputData( int elemIndex,
                                   ElementType etype,
                                   int write_flag,
                                   const TacsScalar Xpts[],
                                   const TacsScalar vars[],
                                   const TacsScalar dvars[],
                                   const TacsScalar ddvars[],
                                   int ld_data,
                                   TacsScalar *data ){
  int num_vis_nodes = TacsGetNumVisNodes(basis->getLayoutType());
  const int vars_per_node = model->getVarsPerNode();

  // Write out the output data
  for ( int i = 0; i < num_vis_nodes; i++ ){
    double pt[3];
    basis->getVisPoint(i, pt);

    // Get the field gradient information
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar U[MAX_VARS_PER_NODE], Udot[MAX_VARS_PER_NODE];
    TacsScalar Uddot[MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    basis->getFieldGradient(pt, Xpts, vars_per_node,
                            vars, dvars, ddvars,
                            X, Xd, J, U, Udot, Uddot, Ud, Ux);

    // Evaluate the output from the data
    model->getOutputData(elemIndex, etype, write_flag, pt,
                         X, U, Udot, Uddot, Ux, ld_data, data);

    data += ld_data;
  }
}
