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
    TacsScalar Ux[3*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    TacsScalar DUt[MAX_VARS_PER_NODE*3], DUx[6*MAX_VARS_PER_NODE];
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
    TacsScalar Ux[3*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, U, Udot, Uddot, Ux);

    // Multiply the weight by the quadrature point
    detJ *= weight;

    // Evaluate the weak form of the model
    int DDUt_nnz, DDUx_nnz;
    const int *DDUt_pairs, *DDUx_pairs;
    TacsScalar DUt[MAX_VARS_PER_NODE*3], DUx[6*MAX_VARS_PER_NODE];
    TacsScalar DDUt[9*MAX_VARS_PER_NODE*MAX_VARS_PER_NODE];
    TacsScalar DDUx[32*MAX_VARS_PER_NODE*MAX_VARS_PER_NODE];
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
