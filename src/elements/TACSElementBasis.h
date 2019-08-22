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

#ifndef TACS_ELEMENT_BASIS_H
#define TACS_ELEMENT_BASIS_H

#include "TACSObject.h"

/*
  This virtual base class defines the interface for the basis functions
  and quadrature schemes used in most TACS elements.

  These are designed to provide common quadrature, interpolation, and
  transformation computations needed for finite element computations.
  This is also designed to capture
*/

class TACSElementBasis {
 public:
  static const int MAX_BASIS_SIZE = 216;

  /**
    Get the layout type
  */
  virtual getLayoutType(){}

  /**
    Get the number of basis functions
  */
  virtual int getNumNodes() = 0;
  
  /**
    Get the spatial dimension of parameter input space = 1, 2, or 3
  */
  virtual int getNumParameters() = 0;
  
  /**
    Get the number of quadrature points for the volume/area of the element
  */
  virtual int getNumQuadraturePoints() = 0;

  /**
    Get the quadrature weight for the n-th quadrature point

    @param n The quadrature point index
    @return The quadrature weight value
  */
  virtual double getQuadratureWeight( int n ) = 0;

  /**
    Get the parametric location of the n-th quadrature point

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @return The quadrature weight value
  */
  virtual double getQuadraturePoint( int n, double pt[] ) = 0;

  /**
    Get the number of faces or edges for the element

    @return The number of faces/edges for the basis
  */
  virtual int getNumElementFaces() = 0;

  /**
    Get the number of quadrature points for the given face

    @param face The face/edge index
    @return The number of quadrature points for the face
  */
  virtual int getNumFaceQuadraturePoints( int face ) = 0;

  /**
    Get the quadrature point for the given face/edge

    The quadrature point and weight are in the original parameter space
    (not parametrized along an edge or face). The tangent parameter
    direction(s) correspond to the directions in parameter space along
    the specified surface. In the case when the parameter space is
    of dimention 1, 2, or 3, there are respectively 0, 1 and 2 tagents
    stored in row major order so that for the 3D case:

    tangent = [d1[0], d1[1], d1[2], d2[0], d2[1], d2[2]]

    Note that the tangents obey the right-hand rule so that
    crossProduct(Xd*d1, Xd*d2) gives an outward-facing normal direction.

    @param face The face/edge index
    @param n The quadrautre point index
    @param pt The quadrature point
    @param tangent Parametric direction(s) parallel to the face
    @return The quadrature weight for the face
  */
  virtual double getFaceQuadraturePoint( int face, int n, double pt[],
                                         double tangent[] ) = 0;

  /**
    Get the field values at the specified quadrature point

    @param n The quadrature point
    @param Xpts The element node locations
    @param vars The element variable values
    @param X The computed coordinate location at quadrature point n
    @param U The computed field values at quadrature point n
  */
  virtual void getFieldValues( const double pt[],
                               const TacsScalar Xpts[],
                               const int vars_per_node,
                               const TacsScalar vars[],
                               TacsScalar X[],
                               TacsScalar U[] );

  /**
    Get the field values at the specified quadrature point

    @param n The quadrature point
    @param Xpts The element node locations
    @param vars The element variable values
    @param X The computed coordinate location at quadrature point n
    @param U The computed field values at quadrature point n
  */
  virtual void getFieldValues( int n,
                               const TacsScalar Xpts[],
                               const int vars_per_node,
                               const TacsScalar vars[],
                               TacsScalar X[],
                               TacsScalar U[] );

  /**
    Get the gradient of the field at the quadrature point

    Note that all matrices are row-major order.

    @param n The quadrature point index
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param vars The element state variables
    @param dvars The first time derivative of the element state vars
    @param ddvars The second time derivative of the element state vars
    @param X The physical quadrature point
    @param U The variables at the quadrature point
    @param Udot The time derivative of the variables
    @param Uddot The second time derivative of the variables
    @param Ux The spatial derivatives of the variables at the parametric location 
  */
  virtual TacsScalar getFieldGradient( const double pt[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[],
                                       const TacsScalar ddvars[],
                                       TacsScalar X[],
                                       TacsScalar Xd[],
                                       TacsScalar J[],
                                       TacsScalar U[],
                                       TacsScalar Udot[],
                                       TacsScalar Uddot[],
                                       TacsScalar Ux[] );

  /**
    Add the weak form of the governing equations to the residual

    @param n The quadrautre point index
    @param pt The quadrature point value
    @param weight The quadrature weight
    @param J The Jacobian coordinate transformation
    @param vars_per_node The number of variables per node
    @param DUt The coefficients of the temporal part of the weak form
    @param DUx The coefficients of the spatial part of the weak form
  */
  virtual void addWeakFormResidual( int n, const double pt[],
                                    TacsScalar weight,
                                    const TacsScalar J[],
                                    const int vars_per_node,
                                    const TacsScalar DUt[],
                                    const TacsScalar DUx[],
                                    TacsScalar res[] );

  /**
    Add the weak form of the governing equations to the residual

    @param n The quadrautre point index
    @param pt The quadrature point value
    @param weight The quadrature weight
    @param J The Jacobian coordinate transformation
    @param vars_per_node The number of variables per node
    @param DUt The coefficients of the temporal part of the weak form
    @param DUx The coefficients of the spatial part of the weak form
  */
  virtual void addWeakFormJacobian( int n, const double pt[],
                                    TacsScalar weight,
                                    const TacsScalar J[],
                                    const int vars_per_node,
                                    const TacsScalar DUt[],
                                    const TacsScalar DUx[],
                                    double alpha, double beta, double gamma,
                                    const int DDUt_nnz,
                                    const int *DDUt_pairs,
                                    const TacsScalar *DDUt,
                                    const int DDUx_nnz,
                                    const int *DDUx_paris,
                                    const TacsScalar *DDUx,
                                    TacsScalar res[],
                                    TacsScalar mat[] );

  /**
    Evaluate basis functions at a parametric point
    
    @param pt The parametric point
    @parma N The shape function values
  */
  virtual void computeBasis( const double pt[], double N[] ) = 0;

  /**
    Compute the derivative of the basis functions with respect to the
    parametric coordinates. This is stored by basis function in
    coordinate order, i.e. Nx = [N[0],pt[0], N[0],pt[1], N[0],pt[2],
    N[1],pt[0] ...
    
    @param pt The parametric point
    @parma N The shape function values
    @param Nxi The derivative of the shape functions w.r.t. the parameters   
  */
  virtual void computeBasisGradient( const double pt[], double N[],
                                     double Nxi[] ) = 0;

  /**
    Compute the derivative of the basis functions with respect to 

    @param n Index of the parametric point
    @parma N The shape function values
  */
  virtual void computeBasis( int n, double N[] ){
    double pt[3];
    getQuadraturePoint(n, pt);
    computeBasis(pt, N);
  }

  /**
    Compute the derivative of the basis functions with respect to 

    @param n Index of the parametric point
    @parma N The shape function values
  */
  virtual void computeBasisGradient( int n, double N[], double Nxi[] ){
    double pt[3];
    getQuadraturePoint(n, pt);
    computeBasis(pt, N, Nxi);
  }

 protected:
  static void computeFieldValues( const int num_nodes,
                                  const double N[],
                                  const TacsScalar Xpts[],
                                  const int vars_per_node,
                                  const TacsScalar vars[],
                                  TacsScalar X[],
                                  TacsScalar U[] );
  static void computeFieldGradient( const int num_params,
                                    const int num_nodes,
                                    const double N[],
                                    const double Nxi[],
                                    const TacsScalar Xpts[],
                                    const int vars_per_node,
                                    const TacsScalar vars[],
                                    const TacsScalar dvars[],
                                    const TacsScalar ddvars[],
                                    TacsScalar X[],
                                    TacsScalar Xd[],
                                    TacsScalar J[],
                                    TacsScalar U[],
                                    TacsScalar Udot[],
                                    TacsScalar Uddot[],
                                    TacsScalar Ux[] );
  static void addWeakFormResidual( const int num_params,
                                   const int num_nodes,
                                   const double N[],
                                   const double Nxi[],
                                   const TacsScalar J[],
                                   TacsScalar weight,
                                   const int vars_per_node,
                                   const TacsScalar DUt[],
                                   const TacsScalar DUx[],
                                   TacsScalar res[] );
  static void addWeakFormJacobian( const int num_params,
                                   const int num_nodes,
                                   const double N[],
                                   const TacsScalar Nx[],
                                   const TacsScalar J[],
                                   TacsScalar weight,
                                   const int vars_per_node,
                                   const TacsScalar DUt[],
                                   const TacsScalar DUx[],
                                   double alpha, double beta, double gamma,
                                   const int DDUt_nnz,
                                   const int *DDUt_pairs,
                                   const TacsScalar *DDUt,
                                   const int DDUx_nnz,
                                   const int *DDUx_paris,
                                   const TacsScalar *DDUx,
                                   TacsScalar *res,
                                   TacsScalar *mat );
};

#endif // TACS_ELEMENT_BASIS_H
