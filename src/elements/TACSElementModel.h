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

#include "TACSObject.h"
#include "TACSElementTypes.h"

/**
  TACSElementModel defines a physical model class independent of a
  finite element basis
*/
class TACSElementModel {
 public:
  /**
    Returns the spatial dimension of the element: 1, 2 or 3

    @return Degrees of freedom per node
  */
  virtual int getSpatialDim() = 0;

  /**
    Returns the number of degrees of freedom per node

    @return Degrees of freedom per node
  */
  virtual int getVarsPerNode() = 0;

  /**
    Retrieve the global design variable numbers associated with this element

    Note when the dvNums argument is NULL, then the result is a query
    on the number of design variables and the array is not set.

    @param dvLen The length of the array dvNums
    @param dvNums An array of the design variable numbers for this element
    @return The number of design variable numbers defined by the element
  */
  virtual int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] ){
    return 0;
  }

  /**
    Set the element design variables from the design vector

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
  */
  virtual void setDesignVars( int elemIndex,
                              int dvLen, const TacsScalar dvs[] ){}

  /**
    Get the element design variables values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param dvs The design variable values
  */
  virtual void getDesignVars( int elemIndex,
                              int dvLen, TacsScalar dvs[] ){}

  /**
    Get the lower and upper bounds for the design variable values

    @param elemIndex The local element index
    @param dvLen The length of the design array
    @param lowerBound The design variable lower bounds
    @param lowerBound The design variable upper bounds
  */
  virtual void getDesignVarRange( int elemIndex, int dvLen,
                                  TacsScalar lowerBound[],
                                  TacsScalar upperBound[] ){}

  /**
    Evaluate the point-wise integrand for the weak form of the governing
    equations of motion.

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
    @param n The quadrature point index
    @param time The simulation time
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param DUt Coefficients of the time-dependent weak form
    @param DUx Coefficients of the spatial-derivative weak form
  */
  virtual void evalWeakIntegrand( int elemIndex,
                                  int n,
                                  const double time,
                                  const double pt[],
                                  const TacsScalar X[],
                                  const TacsScalar Ut[],
                                  const TacsScalar Ux[],
                                  TacsScalar DUt[],
                                  TacsScalar DUx[] ) = 0;

  /**
    Evaluate the point-wise integrand for the weak form of the governing
    equations of motion.

    The following code computes the weak form coefficients and their
    derivatives with respect to each of the input components.

    The descriptions of the terms DUt and DUx are the same as the
    evalWeakIntegrand() function described above.

    The parameter *DDUt* contains a sparse matrix representation of the
    the derivatives of the coefficients in DUt with respect to the
    displacement and their first and second time derivatives. A complete
    matrix DDUt would be a matrix of dimension
    (3*vars_per_node X 3*vars_per_node).

    For instance, for the problem U = (u, v), the full DDUt would be a
    6 X 6 matrix whose first row would contain

    DDUt[0] = d(DUt[0])/du
    DDUt[1] = d(DUt[0])/d(dot{u})
    DDUt[2] = d(DUt[0])/d(ddot{u})
    DDUt[3] = d(DUt[0])/dv
    DDUt[4] = d(DUt[0])/d(dot{v})
    DDUt[5] = d(DUt[0])/d(ddot{v})

    However, this matrix is usually very sparse. For this reason, the
    a simple non-zero pattern is returned in DDT_non_zero_paris[], which
    contains the entries of the matrix which are non-zero. For instance
    if: DUt[2] = rho*Uddot[0] and DUt[5] = rho*Uddot[1], then the output
    would be

    DDt_num_non_zeros = 2
    DDt_non_zero_pairs = [2, 2, 5, 5]
    DDt[0] = rho
    DDt[1] = rho

    @param elemIndex The local element index
    @param n The quadrature point index
    @param time The simulation time
    @param pt The parametric position of the quadrature point
    @param X The physical position of the quadrature point
    @param Ut Values of the state variables and their 1st/2nd time derivs
    @param Ux The spatial derivatives of the state variables
    @param DUt Coefficients of the time-dependent weak form
    @param DUx Coefficients of the spatial-derivative weak form
    @param DDUt_num_non_zeros Number of non-zeros (negative for dense matrix)
    @param DDUt_non_zero_pairs Non-zero Jacobian matrix pairs for DDt
    @param DDUt Jacobian of the time-dependent weak form
  */
  virtual void evalWeakJacobian( int elemIndex,
                                 int n,
                                 const double time,
                                 const double pt[],
                                 const TacsScalar X[],
                                 const TacsScalar Ut[],
                                 const TacsScalar Ux[],
                                 TacsScalar DUt[],
                                 TacsScalar DUx[],
                                 int *Jac_nnz,
                                 const int *Jac_Pairs[],
                                 TacsScalar Jac[] ) = 0;

  /*

  */
  virtual void addWeakAdjProduct( int elemIndex,
                                  int n,
                                  const double time,
                                  const double pt[],
                                  const TacsScalar X[],
                                  const TacsScalar Ut[],
                                  const TacsScalar Ux[],
                                  const TacsScalar Psi[],
                                  const TacsScalar Psix[],
                                  TacsScalar scale,
                                  int dvLen,
                                  TacsScalar *fdvSens ){}

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
  virtual void getOutputData( int elemIndex,
                              const double time,
                              ElementType etype,
                              int write_flag,
                              const double pt[],
                              const TacsScalar X[],
                              const TacsScalar Ut[],
                              const TacsScalar Ux[],
                              int ld_data,
                              TacsScalar *data ) = 0;
};

#endif // TACS_ELEMENT_MODEL_H
