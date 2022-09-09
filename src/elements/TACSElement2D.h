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
#include "TACSElementBasis.h"
#include "TACSElementModel.h"

class TACSElement2D : public TACSElement {
 public:
  static const int MAX_VARS_PER_NODE = 8;

  TACSElement2D(TACSElementModel *_model, TACSElementBasis *_basis);
  ~TACSElement2D();

  // Get the layout properties of the element
  int getVarsPerNode();
  int getNumNodes();
  int getDesignVarsPerNode();
  ElementLayout getLayoutType();
  ElementType getElementType();
  TACSElementBasis *getElementBasis();
  TACSElementModel *getElementModel();
  TACSElement *createElementTraction(int faceIndex, const TacsScalar t[]);
  TACSElement *createElementPressure(int faceIndex, TacsScalar p);
  TACSElement *createElementInertialForce(const TacsScalar inertiaVec[]);
  int getNumQuadraturePoints();
  double getQuadratureWeight(int n);
  double getQuadraturePoint(int n, double pt[]);
  int getNumElementFaces();
  int getNumFaceQuadraturePoints(int face);
  double getFaceQuadraturePoint(int face, int n, double pt[], double tangent[]);

  /**
    Retrieve the global design variable numbers associated with this element
  */
  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]);

  /**
    Set the element design variables from the design vector
  */
  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]);

  /**
    Get the element design variables values
  */
  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]);

  /**
    Get the lower and upper bounds for the design variable values
  */
  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]);

  /**
    Add the residual to the provided vector
  */
  void addResidual(int elemIndex, double time, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res);

  /**
    Add the residual and Jacobians to the arrays
  */
  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res, TacsScalar *mat);

  /**
    Add the derivative of the product of the adjoint variables w.r.t.
    the material design variables
  */
  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dvSens[]);

  /**
    Add the derivative of the product of the adjoint variables and the
    residuals with respect to the node locations
  */
  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

  /**
    Get the size of the data for the matrix-vector product
  */
  void getMatVecDataSizes(ElementMatrixType matType, int elemIndex,
                          int *_data_size, int *_temp_size);

  /**
    Get the data for a matrix vector product. When data is NULL, the function
    returns the size of the required array
  */
  void getMatVecProductData(ElementMatrixType matType, int elemIndex,
                            double time, TacsScalar alpha, TacsScalar beta,
                            TacsScalar gamma, const TacsScalar Xpts[],
                            const TacsScalar vars[], const TacsScalar dvars[],
                            const TacsScalar ddvars[], TacsScalar data[]);

  /**
    Compute the matrix-vector product with the given element data
  */
  void addMatVecProduct(ElementMatrixType matType, int elemIndex,
                        const TacsScalar data[], TacsScalar temp[],
                        const TacsScalar px[], TacsScalar py[]);

  /**
    Compute a specific type of element matrix (mass, stiffness, geometric
    stiffness, etc.)
  */
  void getMatType(ElementMatrixType matType, int elemIndex, double time,
                  const TacsScalar Xpts[], const TacsScalar vars[],
                  TacsScalar mat[]);

  /**
    Add the derivative of the product of a specific matrix w.r.t.
    the design variables
  */
  void addMatDVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, TacsScalar scale,
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[], int dvLen,
                                TacsScalar dfdx[]);

  /**
    Compute the derivative of the product of a specific matrix w.r.t.
    the input variables (vars).
  */
  void getMatSVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, const TacsScalar psi[],
                                const TacsScalar phi[], const TacsScalar Xpts[],
                                const TacsScalar vars[], TacsScalar dfdu[]);

  /**
    Evaluate a point-wise quantity of interest.
  */
  int evalPointQuantity(int elemIndex, int quantityType, double time, int n,
                        double pt[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar *detXd,
                        TacsScalar *quantity);

  /**
    Add the derivative of the point quantity w.r.t. the design variables
  */
  void addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                              TacsScalar scale, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], int dvLen,
                              TacsScalar dfdx[]);

  /**
    Add the derivative of the point quantity w.r.t. the state variables
  */
  void addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                              TacsScalar alpha, TacsScalar beta,
                              TacsScalar gamma, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], TacsScalar dfdu[]);

  /**
    Add the derivative of the point quantity w.r.t. the node locations
  */
  void addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfddetXd,
                               const TacsScalar dfdq[], TacsScalar dfdXpts[]);

  /**
    Compute the output data for visualization
  */
  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

 private:
  TACSElementModel *model;
  TACSElementBasis *basis;
};

#endif  // TACS_ELEMENT_2D_H
