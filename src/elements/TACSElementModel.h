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

#ifndef TACS_ELEMENT_MODEL_H
#define TACS_ELEMENT_MODEL_H

#include "TACSConstitutive.h"
#include "TACSElementTypes.h"
#include "TACSObject.h"

/**
  TACSElementModel defines a physical model class independent of a
  finite element basis
*/
class TACSElementModel : public TACSObject {
 public:
  /**
    Returns the spatial dimension of the element: 1, 2 or 3

    @return The number of computational parameters associated with the model
  */
  virtual int getNumParameters() = 0;

  /**
    Returns the number of degrees of freedom per node

    @return Degrees of freedom per node
  */
  virtual int getVarsPerNode() = 0;

  /**
    Return the number of design variables per "design node"
  */
  virtual int getDesignVarsPerNode() { return 1; }

  /**
    Retrieve the global design variable numbers associated with this element

    Note when the dvNums argument is NULL, then the result is a query
    on the number of design variables and the array is not set.

    @param dvLen The length of the array dvNums
    @param dvNums An array of the design variable numbers for this element
    @return The number of design variable numbers defined by the model
  */
  virtual int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
    return 0;
  }

  /**
    Set the element design variables from the design vector

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
    @return The number of design variable numbers defined by the model
  */
  virtual int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    return 0;
  }

  /**
    Get the element design variables values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
    @return The number of design variable numbers defined by the model
  */
  virtual int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    return 0;
  }

  /**
    Get the lower and upper bounds for the design variable values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param lowerBound The design variable lower bounds
    @param lowerBound The design variable upper bounds
    @return The number of design variable numbers defined by the model
  */
  virtual int getDesignVarRange(int elemIndex, int dvLen,
                                TacsScalar lowerBound[],
                                TacsScalar upperBound[]) {
    return 0;
  }

  /**
    Get underlying tacs constitutive object for this model

    @return The TACSConstitutive object used by this model
  */
  virtual TACSConstitutive *getConstitutive() { return NULL; }

  /**
    Evaluate the point-wise integrand for the weak form of the governing
    equation.

    The weak form consists of two groups of components, the coefficients
    of time-dependent terms (up to second-order), and the coefficients of
    the spatial derivative terms (only first-order).

    Note that we assume separability between the spatial derivatives and the
    temporal derivatives, so that DUt[] does not depend on Ux, and DUx does
    not depend on Udot or Uddot.

    The parameter *DUt* contains the time coefficients in the weak form in
    a matrix of size (vars_per_node x 3). Each column in the matrix represents
    the zero-th, first and second time derivatives with the rows representing
    each variable. Therefore, the weak form for a problem with the variable
    components U = (u, v) will have the following form:

    int_{Area} (DUt[0]*du + DUt[1]*d(dot{u}) + DUt[2]*d(ddot{u}) +
                DUt[3]*dv + DUt[4]*d(dot{v}) + DUt[5]*d(ddot{v}) +
                spatial terms) dA = 0

    The parameter *DUx* contains the spatial derivative components of the
    weak form in a matrix of size (vars_per_node x (spatial_dim + 1)).
    The first component represents the coefficient of the variable, while
    the second, third and possibly fourth component represent the remaining
    spatial derivative coefficients. A problem with the variable
    components U = (u, v) with a spatial dimension of two will have the
    following weak form:

    int_{Area} (time dependent terms +
                DUx[0]*du + DUx[1]*d(du/dx) + DUx[2]*d(du/dy)) +
                DUx[3]*dv + DUx[4]*d(dv/dx) + DUx[5]*d(dv/dy)) dA = 0

    Note that the coefficients DUt[0] and DUx[0], and DUt[3] and DUx[3],
    are both for the displacement u and v, respectively. This means that
    the implementation is not unique.

    @param elemIndex The local element index
    @param time The simulation time
    @param n The quadrature point index
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Xd The derivative physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param DUt Coefficients of the time-dependent weak form
    @param DUx Coefficients of the spatial-derivative weak form
  */
  virtual void evalWeakIntegrand(int elemIndex, const double time, int n,
                                 const double pt[], const TacsScalar X[],
                                 const TacsScalar Xd[], const TacsScalar Ut[],
                                 const TacsScalar Ux[], TacsScalar DUt[],
                                 TacsScalar DUx[]) = 0;

  /**
     Add the product of the adjoint with the derivative w.r.t. design
     variables to the design variable vector.

     @param elemIndex The local element index
     @param time The simulation time
     @param scale The scaling factor applied to the derivative
     @param n The quadrature point index
     @param pt The parametric position of the quadrature point
     @param X The physical position of the quadrature point
     @param Xd The derivative physical position of the quadrature point
     @param Ut Values of the state variables and their 1st/2nd time derivs
     @param Ux The spatial derivatives of the state variables
     @param Psi The adjoint variable values
     @param Psix The spatial derivatives of the adjoint variable values
     @param dvLen The length of the design variable vector
     @param fdvSens The derivative vector
  */
  virtual void addWeakAdjProduct(int elemIndex, const double time,
                                 TacsScalar scale, int n, const double pt[],
                                 const TacsScalar X[], const TacsScalar Xd[],
                                 const TacsScalar Ut[], const TacsScalar Ux[],
                                 const TacsScalar Psi[],
                                 const TacsScalar Psix[], int dvLen,
                                 TacsScalar *fdvSens) {}

  /**
    Compute the spatial derivatives of the point-wise product of the adjoint

    @param elemIndex The local element index
    @param time The simulation time
    @param n The quadrature point index
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Xd The derivative physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param Psi The adjoint variable values
    @param Psix The spatial derivatives of the adjoint variable values
    @param product The product of the adjoint and element residual vector
    @param dfdX The derivative of the product w.r.t. X
    @param dfdXd The derivative of the product w.r.t. Xd
    @param dfdUx The derivative of the product with respect to Ux
    @param dfdPsix The derivative of the product with respect to Psix
  */
  virtual void evalWeakAdjXptSensProduct(
      int elemIndex, const double time, int n, const double pt[],
      const TacsScalar X[], const TacsScalar Xd[], const TacsScalar Ut[],
      const TacsScalar Ux[], const TacsScalar Psi[], const TacsScalar Psix[],
      TacsScalar *product, TacsScalar dfdX[], TacsScalar dfdXd[],
      TacsScalar dfdUx[], TacsScalar dfdPsix[]) {
    *product = 0.0;
    dfdX[0] = dfdX[1] = dfdX[2] = 0.0;
    const int vars_per_node = getVarsPerNode();
    const int num_params = getNumParameters();

    for (int i = 0; i < num_params * vars_per_node; i++) {
      dfdUx[i] = 0.0;
      dfdPsix[i] = 0.0;
    }

    for (int i = 0; i < 3 * num_params; i++) {
      dfdXd[i] = 0.0;
    }
  }

  /**
    Get the non-zero pattern of the element type

    @param matType The type of matrix to use
    @param elemIndex The element index
    @param Jac_nnz Number of non-zeros
    @param Jac_pairs Indicex pairs of the non-zero entries in the Jacobian
    matrix
  */
  virtual void getWeakMatrixNonzeros(ElementMatrixType matType, int elemIndex,
                                     int *Jac_nnz, const int *Jac_pairs[]) = 0;

  /**
    Evaluate the Jacobian of the point-wise integrand for the weak form of
    the governing equations.

    The following code computes the weak form coefficients and their
    derivatives with respect to each of the input components. The
    descriptions of the terms DUt and DUx are the same as the
    evalWeakIntegrand() function described above.

    The parameter Jac contains a sparse matrix representation of the
    the derivatives of the coefficients in DUt and DUx. The matrix
    contains (3 + spatial_dim)*vars_per_node rows and columns.

    For instance, for the 2D problem (spatial_dim = 2) with the variables
    U = (u, v), the Jac matrix would contain 10 x 10 entries. The rows of the
    matrix (corresponding to DUt and DUx) are ordered first by variable, then
    by derivative. The columns of the matrix are ordered in a similar manner
    so that for this case:

    Index:     0;       1;      2;      3;      4;
    rows:  DUt[0]; DUt[1]; DUt[2]; DUx[0]; DUx[1];
    cols:      u;     u,t;   u,tt;    u,x;    u,y;

    Index:      5;      6;      7;      8;      9;
    rows:  DUt[3]; DUt[4]; DUt[5]; DUx[2]; DUx[3];
    cols:       v;    v,t;   v,tt;    v,x;    v,y;

    However, the Jacobian matrix of the terms DUt/DUx w.r.t. Ut and Ux is
    often sparse. For this reason, the sparsity pattern is returned in a
    pair-wise format with in Jac_pairs which stores the (row, column) entries
    that are non-zero.

    @param matType The element matrix type
    @param elemIndex The local element index
    @param time The simulation time
    @param n The quadrature point index
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Xd The derivative physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param DUt Coefficients of the time-dependent weak form
    @param DUx Coefficients of the spatial-derivative weak form
    @param Jac Jacobian entries of the weak form
  */
  virtual void evalWeakMatrix(ElementMatrixType matType, int elemIndex,
                              const double time, int n, const double pt[],
                              const TacsScalar X[], const TacsScalar Xd[],
                              const TacsScalar Ut[], const TacsScalar Ux[],
                              TacsScalar DUt[], TacsScalar DUx[],
                              TacsScalar Jac[]) = 0;

  /**
    Add the derivative of the matrix inner product with respect to the
    design variables to the element design vector.

    @param elemIndex The local element index
    @param time The simulation time
    @param scale Scaling parameter
    @param n The quadrature point index
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Xd The derivative physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param Psi Values of the left vector
    @param Psix Spatial derivatives of the left vector
    @param Phi Value of the right vector
    @param Phix Spatial derivative of the right vector
    @param dvLen The length of the design sensitivity array
    @param dfdx The derivative array being computed
  */
  virtual void addWeakMatDVSens(ElementMatrixType matType, int elemIndex,
                                const double time, TacsScalar scale, int n,
                                const double pt[], const TacsScalar X[],
                                const TacsScalar Xd[], const TacsScalar Ut[],
                                const TacsScalar Ux[], const TacsScalar Psi[],
                                const TacsScalar Psix[], const TacsScalar Phi[],
                                const TacsScalar Phix[], int dvLen,
                                TacsScalar dfdx[]) {}

  /**
    Evaluate the derivative of the matrix inner product with respect to the
    state variables U and Ux.

    @param elemIndex The local element index
    @param time The simulation time
    @param scale Scaling parameter
    @param n The quadrature point index
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Xd The derivative physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param Psi Values of the left vector
    @param Psix Spatial derivatives of the left vector
    @param Phi Value of the right vector
    @param Phix Spatial derivative of the right vector
    @param dfdU Derivative of the inner product w.r.t. U
    @param dfdUx Derivative of the inner product w.r.t. Ux
  */
  virtual void evalWeakMatSVSens(
      ElementMatrixType matType, int elemIndex, const double time,
      TacsScalar scale, int n, const double pt[], const TacsScalar X[],
      const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
      const TacsScalar Psi[], const TacsScalar Psix[], const TacsScalar Phi[],
      const TacsScalar Phix[], TacsScalar dfdU[], TacsScalar dfdUx[]) {
    const int vars_per_node = getVarsPerNode();
    const int num_params = getNumParameters();

    for (int i = 0; i < 3 * vars_per_node; i++) {
      dfdU[i] = 0.0;
    }
    for (int i = 0; i < num_params * vars_per_node; i++) {
      dfdUx[i] = 0.0;
    }
  }

  /**
     Evaluate a point-wise quantity of interest at a quadrature point

     This function computes a local pointwise quantity of interest
     within the element. The quantity type is determined based on the
     argument quatityType value. Note that specific integers are not
     assigned, so you can make up new values and add them as
     needed.

     @param elemIndex The local element index
     @param quantityType Integer indicating the type of pointwise quantity
     @param time The simulation time
     @param n The quadrature point index
     @param pt The parametric position of the quadrature point
     @param X The physical position of the quadrature point
     @param Xd The derivative of the position w.r.t. the parameters
     @param Ut Values of the state variables and their 1st/2nd time derivs
     @param Ux The spatial derivatives of the state variables
     @param quantity The quantity of interest
     @return Length of the quantity computed (a scalar or physical vector)
  */
  virtual int evalPointQuantity(int elemIndex, const int quantityType,
                                const double time, int n, const double pt[],
                                const TacsScalar X[], const TacsScalar Xd[],
                                const TacsScalar Ut[], const TacsScalar Ux[],
                                TacsScalar *quantity) {
    return 0;
  }

  /**
     Add the derivative of the quantity w.r.t. the design variables

     @param elemIndex The local element index
     @param quantityType Integer indicating the type of pointwise quantity
     @param time The simulation time
     @param scale Scale factor applied to the derivative
     @param n The quadrature point index
     @param pt The parametric position of the quadrature point
     @param X The physical position of the quadrature point
     @param Xd The derivative of the position w.r.t. the parameters
     @param Ut Values of the state variables and their 1st/2nd time derivs
     @param Ux The spatial derivatives of the state variables
     @param dfdq The derivative of the function of interest w.r.t. the quantity
     @param dvLen The length of the design sensitivity array
     @param dfdx The derivative array being computed
  */
  virtual void addPointQuantityDVSens(
      int elemIndex, const int quantityType, const double time,
      TacsScalar scale, int n, const double pt[], const TacsScalar X[],
      const TacsScalar Xd[], const TacsScalar Ut[], const TacsScalar Ux[],
      const TacsScalar dfdq[], int dvLen, TacsScalar dfdx[]) {}

  /**
     Evaluate the derivatives of the point-wise quantity of interest
     with respect to X, Ut and Ux.

     @param elemIndex The local element index
     @param quantityType Integer indicating the type of pointwise quantity
     @param time The simulation time
     @param n The quadrature point index
     @param pt The parametric position of the quadrature point
     @param X The physical position of the quadrature point
     @param Xd The derivative of the position w.r.t. the parameters
     @param Ut Values of the state variables and their 1st/2nd time derivs
     @param Ux The spatial derivatives of the state variables
     @param dfdq The derivative of the function of interest w.r.t. the quantity
     @param dfdX The derivative of the function w.r.t. X
     @param dfdX The derivative of the function w.r.t. Ut
     @param dfdX The derivative of the function w.r.t. Ux
  */
  virtual void evalPointQuantitySens(
      int elemIndex, const int quantityType, const double time, int n,
      const double pt[], const TacsScalar X[], const TacsScalar Xd[],
      const TacsScalar Ut[], const TacsScalar Ux[], const TacsScalar dfdq[],
      TacsScalar dfdX[], TacsScalar dfdXd[], TacsScalar dfdUt[],
      TacsScalar dfdUx[]) {}

  /**
    Generate a line of output for a single visualization point

    @param elemIndex The local element index
    @param etype The class of element output to generate
    @param write_flag The flag to indicate which components to write
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param U The values of the state variables
    @param Udot The time derivatives of the state variables
    @param Uddot The second time derivatives of the state variables
    @param Ux The spatial derivatives of the state variables
  */
  virtual void getOutputData(int elemIndex, const double time,
                             ElementType etype, int write_flag,
                             const double pt[], const TacsScalar X[],
                             const TacsScalar Ut[], const TacsScalar Ux[],
                             int ld_data, TacsScalar *data) {}
};

#endif  // TACS_ELEMENT_MODEL_H
