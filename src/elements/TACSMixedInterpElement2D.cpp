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

#include "TACSMixedInterpElement2D.h"

TACSMixedInterpElement2D::TACSMixedInterpElement2D( TASCMixedInterpElementModel *_model,
                                                    TACSMixedInterpElementBasis *_basis ){
  model = _model;  model->incref();
  basis = _basis;  basis->incref();
}

TACSMixedInterpElement2D::~TACSMixedInterpElement2D(){
  model->decref();
  basis->decref();
}

// Get the layout properties of the element
int TACSMixedInterpElement2D::getVarsPerNode(){
  return model->getVarsPerNode();
}

int TACSMixedInterpElement2D::getNumNodes(){
  return basis->getNumNodes();
}

int TACSMixedInterpElement2D::getDesignVarsPerNode(){
  return model->getDesignVarsPerNode();
}

ElementLayout TACSMixedInterpElement2D::getLayoutType(){
  return basis->getLayoutType();
}

TACSElementBasis* TACSMixedInterpElement2D::getElementBasis(){
  return basis;
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSMixedInterpElement2D::getDesignVarNums( int elemIndex,
                                                int dvLen,
                                                int dvNums[] ){
  return model->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSMixedInterpElement2D::setDesignVars( int elemIndex,
                                             int dvLen,
                                             const TacsScalar dvs[] ){
  return model->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSMixedInterpElement2D::getDesignVars( int elemIndex,
                                             int dvLen,
                                             TacsScalar dvs[] ){
  return model->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSMixedInterpElement2D::getDesignVarRange( int elemIndex,
                                                 int dvLen,
                                                 TacsScalar lb[],
                                                 TacsScalar ub[] ){
  return model->getDesignVarRange(elemIndex, dvLen, lb, ub);
}


void TACSMixedInterpElement2D::addResidual( int elemIndex,
                                            double time,
                                            const TacsScalar *Xpts,
                                            const TacsScalar *vars,
                                            const TacsScalar *dvars,
                                            const TacsScalar *ddvars,
                                            TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();
  const int vars_per_node = model->getVarsPerNode();

  // Get the number of tying points
  int ntying = basis->getNumTyingPoints();
  for ( int ty = 0; ty < ntying; ty++ ){
    basis->getTyingPoint(ty, pt);

    // Compute the field gradient at the tying point. Note that this computes the
    // gradient w.r.t. parametric coordiantes only. No coordinate
    // transformation is applied.
    TacsScalar X[3], Xd[4];
    TacsScalar U[MAX_VARS_PER_NODE], Ud[2*MAX_VARS_PER_NODE];
    basis->getFieldGradient(ty, pt, Xpts, vars_per_node, vars, X, Xd, U, Ud);

    // Evaluate the tying quantity based on the tying values
    model->evalTyingQuantity(ty, pt, U, Ud, &qty[ty]);
  }



  // Loop over each quadrature point
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar Ut[3*MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detJ = basis->getFieldGradient(pt, Xpts, vars_per_node,
                                              vars, dvars, ddvars,
                                              X, Xd, J, Ut, Ud, Ux);

    // Evaluate the tying quantities
    basis->evalTyingField(n, pt, qty, Uty);

    // Evaluate the weak integrand for the quantity of interest..
    model->evalWeakIntegrand(elemIndex, n, time, pt, X,
                             Ut, Ux, Uty, DUt, DUx, &DUty[ntying*n]);

    // Add the weak form of the residual at this point
    basis->addWeakFormResidual(n, pt, detJ, J, vars_per_node, DUt, DUx, res);
  }
}


void TACSMixedInterpElement2D::addJacobian( int elemIndex,
                                            double time,
                                            TacsScalar alpha,
                                            TacsScalar beta,
                                            TacsScalar gamma,
                                            const TacsScalar *Xpts,
                                            const TacsScalar *vars,
                                            const TacsScalar *dvars,
                                            const TacsScalar *ddvars,
                                            TacsScalar *res,
                                            TacsScalar *mat ){

}


// Functions for the adjoint
void TACSMixedInterpElement2D::addAdjResProduct( int elemIndex,
                                                 double time,
                                                 TacsScalar scale,
                                                 const TacsScalar psi[],
                                                 const TacsScalar Xpts[],
                                                 const TacsScalar vars[],
                                                 const TacsScalar dvars[],
                                                 const TacsScalar ddvars[],
                                                 int dvLen,
                                                 TacsScalar dfdx[] ){}

void TACSMixedInterpElement2D::addAdjResXptProduct( int elemIndex,
                                                    double time,
                                                    TacsScalar scale,
                                                    const TacsScalar psi[],
                                                    const TacsScalar Xpts[],
                                                    const TacsScalar vars[],
                                                    const TacsScalar dvars[],
                                                    const TacsScalar ddvars[],
                                                    TacsScalar dfdXpts[] ){}

/**
   Compute a specific type of element matrix (mass, stiffness, geometric
   stiffness, etc.)
*/
void TACSMixedInterpElement2D::getMatType( ElementMatrixType matType,
                                           int elemIndex,
                                           double time,
                                           const TacsScalar Xpts[],
                                           const TacsScalar vars[],
                                           TacsScalar mat[] ){}

/**
   Add the derivative of the product of a specific matrix w.r.t.
   the design variables
*/
void TACSMixedInterpElement2D::addMatDVSensInnerProduct( ElementMatrixType matType,
                                                         int elemIndex,
                                                         double time,
                                                         TacsScalar scale,
                                                         const TacsScalar psi[],
                                                         const TacsScalar phi[],
                                                         const TacsScalar Xpts[],
                                                         const TacsScalar vars[],
                                                         int dvLen,
                                                         TacsScalar dfdx[] ){}

/**
   Compute the derivative of the product of a specific matrix w.r.t.
   the input variables (vars).
*/
void TACSMixedInterpElement2D::getMatSVSensInnerProduct( ElementMatrixType matType,
                                                         int elemIndex,
                                                         double time,
                                                         const TacsScalar psi[],
                                                         const TacsScalar phi[],
                                                         const TacsScalar Xpts[],
                                                         const TacsScalar vars[],
                                                         TacsScalar dfdu[] ){}

/**
   Evaluate a point-wise quantity of interest.
*/
int TACSMixedInterpElement2D::evalPointQuantity( int elemIndex,
                                                 int quantityType,
                                                 double time,
                                                 int n, double pt[],
                                                 const TacsScalar Xpts[],
                                                 const TacsScalar vars[],
                                                 const TacsScalar dvars[],
                                                 const TacsScalar ddvars[],
                                                 TacsScalar *quantity ){
  const int vars_per_node = model->getVarsPerNode();
  TacsScalar X[3], Xd[4], J[4];
  TacsScalar Ut[3*MAX_VARS_PER_NODE];
  TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
  basis->getFieldGradient(n, pt, Xpts, vars_per_node, vars, dvars, ddvars,
                          X, Xd, J, Ut, Ud, Ux);

  return model->evalPointQuantity(elemIndex, quantityType, time, n, pt,
                                  X, Xd, Ut, Ux, quantity);
}

/**
   Add the derivative of the point quantity w.r.t. the design variables
*/
void TACSMixedInterpElement2D::addPointQuantityDVSens( int elemIndex,
                                                       int quantityType,
                                                       double time,
                                                       TacsScalar scale,
                                                       int n, double pt[],
                                                       const TacsScalar Xpts[],
                                                       const TacsScalar vars[],
                                                       const TacsScalar dvars[],
                                                       const TacsScalar ddvars[],
                                                       const TacsScalar dfdq[],
                                                       int dvLen, TacsScalar dfdx[] ){
  const int vars_per_node = model->getVarsPerNode();
  TacsScalar X[3], Xd[4], J[4];
  TacsScalar Ut[3*MAX_VARS_PER_NODE];
  TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
  basis->getFieldGradient(n, pt, Xpts, vars_per_node, vars, dvars, ddvars,
                          X, Xd, J, Ut, Ud, Ux);

  model->addPointQuantityDVSens(elemIndex, quantityType, time, scale, n, pt,
                                X, Xd, Ut, Ux, dfdq, dvLen, dfdx);
}

/**
   Add the derivative of the point quantity w.r.t. the state variables
*/
void TACSMixedInterpElement2D::addPointQuantitySVSens( int elemIndex,
                                                       int quantityType,
                                                       double time,
                                                       TacsScalar alpha,
                                                       TacsScalar beta,
                                                       TacsScalar gamma,
                                                       int n, double pt[],
                                                       const TacsScalar Xpts[],
                                                       const TacsScalar vars[],
                                                       const TacsScalar dvars[],
                                                       const TacsScalar ddvars[],
                                                       const TacsScalar dfdq[],
                                                       TacsScalar dfdu[] ){
  const int vars_per_node = model->getVarsPerNode();
  TacsScalar X[3], Xd[4], J[4];
  TacsScalar Ut[3*MAX_VARS_PER_NODE];
  TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
  basis->getFieldGradient(n, pt, Xpts, vars_per_node, vars, dvars, ddvars,
                          X, Xd, J, Ut, Ud, Ux);

  // Evaluate the derivative of the function with respect to X, Ut, Ux
  TacsScalar dfdX[3], dfdXd[4];
  TacsScalar dfdUt[3*MAX_VARS_PER_NODE], dfdUx[2*MAX_VARS_PER_NODE];
  model->evalPointQuantitySens(elemIndex, quantityType, time, n, pt,
                               X, Xd, Ut, Ux, dfdq, dfdX, dfdXd, dfdUt, dfdUx);

  basis->addFieldGradientSVSens(n, pt, Xpts, vars_per_node, Xd, J, Ud,
                                dfdUt, dfdUx, dfdu);
}

/**
   Add the derivative of the point quantity w.r.t. the node locations
*/
void TACSMixedInterpElement2D::addPointQuantityXptSens( int elemIndex,
                                                        int quantityType,
                                                        double time,
                                                        TacsScalar scale,
                                                        int n, double pt[],
                                                        const TacsScalar Xpts[],
                                                        const TacsScalar vars[],
                                                        const TacsScalar dvars[],
                                                        const TacsScalar ddvars[],
                                                        const TacsScalar dfdq[],
                                                        TacsScalar dfdXpts[] ){
  const int vars_per_node = model->getVarsPerNode();
  TacsScalar X[3], Xd[4], J[4];
  TacsScalar Ut[3*MAX_VARS_PER_NODE];
  TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
  basis->getFieldGradient(n, pt, Xpts, vars_per_node, vars, dvars, ddvars,
                          X, Xd, J, Ut, Ud, Ux);

  // Evaluate the derivative of the function with respect to X, Ut, Ux
  TacsScalar dfdX[3], dfdXd[4];
  TacsScalar dfdUt[3*MAX_VARS_PER_NODE], dfdUx[2*MAX_VARS_PER_NODE];
  model->evalPointQuantitySens(elemIndex, quantityType, time, n, pt,
                               X, Xd, Ut, Ux, dfdq, dfdX, dfdXd, dfdUt, dfdUx);

  basis->addFieldGradientXptSens(n, pt, Xpts, vars_per_node, Xd, J, Ud,
                                 0.0, dfdX, dfdXd, NULL, dfdUx, dfdXpts);
}

/*
  Get the element data for the basis
*/
void TACSMixedInterpElement2D::getOutputData( int elemIndex,
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
  for ( int n = 0; n < num_vis_nodes; n++ ){
    double pt[3];
    basis->getVisPoint(n, pt);

    // Get the field gradient information
    TacsScalar X[3], Xd[4], J[4];
    TacsScalar Ut[3*MAX_VARS_PER_NODE];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    basis->getFieldGradient(-1, pt, Xpts, vars_per_node,
                            vars, dvars, ddvars,
                            X, Xd, J, Ut, Ud, Ux);

    // Evaluate the output from the data
    double time = 0.0;
    model->getOutputData(elemIndex, time, etype, write_flag, pt,
                         X, Ut, Ux, ld_data, data);

    data += ld_data;
  }
}
