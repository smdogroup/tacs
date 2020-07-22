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

#include "TACSMITCElement2D.h"

TACSMITCElement2D::TACSMITCElement2D( TACSMITCModel *_model,
                                      TACSMITCBasis *_basis ){
  model = _model;  model->incref();
  basis = _basis;  basis->incref();
}

TACSMITCElement2D::~TACSMITCElement2D(){
  model->decref();
  basis->decref();
}

// Get the layout properties of the element
int TACSMITCElement2D::getVarsPerNode(){
  return model->getVarsPerNode();
}

int TACSMITCElement2D::getNumNodes(){
  return basis->getNumNodes();
}

int TACSMITCElement2D::getDesignVarsPerNode(){
  return model->getDesignVarsPerNode();
}

ElementLayout TACSMITCElement2D::getLayoutType(){
  return basis->getLayoutType();
}

TACSElementBasis* TACSMITCElement2D::getElementBasis(){
  return basis;
}

/*
  Retrieve the global design variable numbers associated with this element
*/
int TACSMITCElement2D::getDesignVarNums( int elemIndex,
                                         int dvLen,
                                         int dvNums[] ){
  return model->getDesignVarNums(elemIndex, dvLen, dvNums);
}

/*
  Set the element design variables from the design vector
*/
int TACSMITCElement2D::setDesignVars( int elemIndex,
                                      int dvLen,
                                      const TacsScalar dvs[] ){
  return model->setDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the element design variables values
*/
int TACSMITCElement2D::getDesignVars( int elemIndex,
                                      int dvLen,
                                      TacsScalar dvs[] ){
  return model->getDesignVars(elemIndex, dvLen, dvs);
}

/*
  Get the lower and upper bounds for the design variable values
*/
int TACSMITCElement2D::getDesignVarRange( int elemIndex,
                                          int dvLen,
                                          TacsScalar lb[],
                                          TacsScalar ub[] ){
  return model->getDesignVarRange(elemIndex, dvLen, lb, ub);
}

void TACSMITCElement2D::addResidual( int elemIndex,
                                     double time,
                                     const TacsScalar *Xpts,
                                     const TacsScalar *vars,
                                     const TacsScalar *dvars,
                                     const TacsScalar *ddvars,
                                     TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = basis->getNumQuadraturePoints();
  const int vars_per_node = model->getVarsPerNode();

  // Compute the values at the tying points
  TacsScalar qty[MAX_TOTAL_TYING_POINTS];

  // Compute the residual contribution associated with each tying point
  TacsScalar rty[MAX_TOTAL_TYING_POINTS];
  memset(rty, 0, MAX_TOTAL_TYING_POINTS*sizeof(TacsScalar));

  // Get the number of tying field values
  const int nfields = basis->getNumTyingFields();
  for ( int field = 0, tyindex = 0; field < nfields; field++ ){
    // Loop over the tying points field values
    const int ntying = basis->getNumTyingPoints(field);
    for ( int ty = 0; ty < ntying; ty++, tyindex++ ){
      double pt[2];
      basis->getTyingPoint(field, ty, pt);

      // Compute the field gradient at the tying point. Note that this computes the
      // gradient w.r.t. parametric coordiantes only. No coordinate
      // transformation is applied.
      TacsScalar Xd[6];
      TacsScalar U[MAX_VARS_PER_NODE], Ud[2*MAX_VARS_PER_NODE];
      basis->getTyingPointGradient(field, ty, pt, Xpts, vars_per_node, vars, Xd, U, Ud);

      // Evaluate the tying quantity based on the tying values
      qty[tyindex] = model->evalTyingQuantity(field, ty, pt, Xd, U, Ud);
    }
  }

  // Loop over each quadrature points
  for ( int n = 0; n < nquad; n++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = basis->getQuadraturePoint(n, pt);

    // Get the solution field and the solution field gradient and the
    // Jacobian transformation
    TacsScalar X[3], Xd[6], J[4];
    TacsScalar Ut[3*MAX_VARS_PER_NODE + MAX_NUM_TYING_FIELDS];
    TacsScalar Ud[2*MAX_VARS_PER_NODE], Ux[2*MAX_VARS_PER_NODE];
    TacsScalar detXd = basis->getFieldGradient(n, pt, Xpts, vars_per_node,
                                               vars, dvars, ddvars,
                                               X, Xd, J, Ut, Ud, Ux);

    // Evaluate the tying quantities
    TacsScalar *Uty = &Ut[3*vars_per_node];
    basis->interpTyingField(n, pt, qty, Uty);

    // Evaluate the weak integrand for the quantity of interest..
    TacsScalar DUt[3*MAX_VARS_PER_NODE + MAX_NUM_TYING_FIELDS];
    TacsScalar DUx[2*MAX_VARS_PER_NODE];
    model->evalWeakIntegrand(elemIndex, n, time, pt, X, Xd,
                             Ut, Ux, DUt, DUx);

    // Multiply the weight by the quadrature point
    detXd *= weight;

    // Add the weak form of the residual
    basis->addWeakResidual(n, pt, detXd, J, vars_per_node, DUt, DUx, res, rty);
  }

  // Get the number of tying field values
  for ( int field = 0, tyindex = 0; field < nfields; field++ ){
    // Loop over the tying points field value
    const int ntying = basis->getNumTyingPoints(field);
    for ( int ty = 0; ty < ntying; ty++, tyindex++ ){
      double pt[2];
      basis->getTyingPoint(field, ty, pt);

      // Evaluate the point and field quantity
      TacsScalar Xd[6];
      TacsScalar U[MAX_VARS_PER_NODE], Ud[2*MAX_VARS_PER_NODE];
      basis->getTyingPointGradient(field, ty, pt, Xpts, vars_per_node, vars, Xd, U, Ud);

      // Evaluate the tying quantity based on the tying values
      TacsScalar DU[MAX_VARS_PER_NODE], DUd[2*MAX_VARS_PER_NODE];
      model->evalTyingQuantitySVSens(field, ty, pt, rty[tyindex], Xd, U, Ud, DU, DUd);

      // Add the weak form derivative
      basis->addTyingResidual(field, ty, pt, vars_per_node, DU, DUd, res);
    }
  }
}

/**
   Evaluate a point-wise quantity of interest.
*/
int TACSMITCElement2D::evalPointQuantity( int elemIndex,
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
void TACSMITCElement2D::addPointQuantityDVSens( int elemIndex,
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
void TACSMITCElement2D::addPointQuantitySVSens( int elemIndex,
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
void TACSMITCElement2D::addPointQuantityXptSens( int elemIndex,
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
void TACSMITCElement2D::getOutputData( int elemIndex,
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
