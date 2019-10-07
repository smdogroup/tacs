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
#include "TACSElementTypes.h"

/*
  This virtual base class defines the interface for the basis functions
  and quadrature schemes used in most TACS elements.

  These are designed to provide common quadrature, interpolation, and
  transformation computations needed for finite element computations.
  This is also designed to capture
*/

class TACSElementBasis : public TACSObject {
 public:
  static const int MAX_BASIS_SIZE = 216;

  /**
    Get the layout type

    @return The element layout type
  */
  virtual ElementLayout getLayoutType();

  /**
    Get the parametric point visualization point

    Note that the number of visualization points must be consistent with
    the number of points defined by the call to TacsGetNumVisNodes()
    with the corresponding element layout type.

    @param n Index for the parametric point for visualization
    @param pt Parametric point location within the element
  */
  virtual void getVisPoint( int n, double pt[] );

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
    Get the face normal at a specified face quadrature point.

    This function returns a 2 or 3-vector depending on the dimension
    of the problem. Note that this function can only be used to evaluate
    the face normal at locations defined by the basis function class.

    @param face The face/edge index
    @param n The quadrautre point index
    @param Xpts The node locations
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param n The face (or edge) normal
  */
  virtual TacsScalar getFaceNormal( int face, int n,
                                    const TacsScalar Xpts[],
                                    TacsScalar Xd[],
                                    TacsScalar normal[] );

  /**
    Add the derivative of the face normal into the nodal sensitivities

  */
  virtual void addFaceNormalXptSens( int face, int n,
                                     const TacsScalar Xpts[],
                                     const TacsScalar A,
                                     const TacsScalar normal[],
                                     const TacsScalar dfdA,
                                     const TacsScalar dfdXd[],
                                     const TacsScalar dfdn[],
                                     TacsScalar dfdXpts[] ){}

  /**
    Get the Jacobian transformation from computational to physical
    coordinates.

    This code returns the determinant of the transformation and
    returns both the Xd, the derivative of physical coordinates
    w.r.t. the parameters, and J, the inverse of Xd. This code can be
    used for computing the volume/area quadratures.

    @param pt The quadrature point
    @param Xpts The node locations
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @return The determinant of Xd
  */
  virtual TacsScalar getJacobianTransform( const double pt[],
                                           const TacsScalar Xpts[],
                                           TacsScalar Xd[],
                                           TacsScalar J[] );

  /**
    Compute the derivative of the Jacobian transformation

    This code adds the contribution to the derivative of the Jacobian
    transformation to the output dfdXpts. Note that the inputs here
    should be generated by a call to getJacobianTransform first.  Note
    that dfdXd and dfdJ may be passed as NULL.

    @param pt The quadrature point
    @param Xpts The node locations
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param dfddetJ The derivative of the function w.r.t. detJ
    @param dfdXd The derivative of the function w.r.t. Xd
    @param dfdJ The derivative of the function w.r.t. J
    @param dfdXpts The output derivative of the function w.r.t. Xpts
  */
  virtual void addJacobianTransformSens( const double pt[],
                                         const TacsScalar Xd[],
                                         const TacsScalar J[],
                                         TacsScalar dfddetJ,
                                         const TacsScalar dfXd[],
                                         const TacsScalar dfdJ[],
                                         TacsScalar dfdXpts[] );
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
    @param Ut The variables and time derivatives at the quadrature point
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param Ux The derivative of the variables w.r.t. the spatial coords
    @return The determinant of the matrix Xd
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
                                       TacsScalar Ut[],
                                       TacsScalar Ud[],
                                       TacsScalar Ux[] );
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
    @param Ut The variables and time derivatives at the quadrature point
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param Ux The derivative of the variables w.r.t. the spatial coords
    @param Psi The adjoint variable values (no time derivatives!)
    @param Psid The derivatives of the adjoint variables w.r.t. parameters
    @param Psix The spatial derivatives of the adjoint variables
    @return The determinant of the matrix Xd
  */
  virtual TacsScalar getFieldGradient( const double pt[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar vars[],
                                       const TacsScalar dvars[],
                                       const TacsScalar ddvars[],
                                       const TacsScalar psi[],
                                       TacsScalar X[],
                                       TacsScalar Xd[],
                                       TacsScalar J[],
                                       TacsScalar Ut[],
                                       TacsScalar Ud[],
                                       TacsScalar Ux[],
                                       TacsScalar Psi[],
                                       TacsScalar Psid[],
                                       TacsScalar Psix[] );

  /**
    Add contributions to the gradient
  */
  virtual void addFieldGradientXptSens( const double pt[],
                                        const TacsScalar Xpts[],
                                        const int vars_per_node,
                                        const TacsScalar Xd[],
                                        const TacsScalar J[],
                                        const TacsScalar Ud[],
                                        const TacsScalar dfddetJ,
                                        const TacsScalar dfdX[],
                                        const TacsScalar dfdXd[],
                                        const TacsScalar dfdJ[],
                                        const TacsScalar dfdUx[],
                                        TacsScalar dfdXpts[] );

  /**
    Add contributions to the gradient
  */
  virtual void addFieldGradientSVSens( const double pt[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar Xd[],
                                       const TacsScalar J[],
                                       const TacsScalar Ud[],
                                       const TacsScalar dfdUt[],
                                       const TacsScalar dfdUx[],
                                       TacsScalar dfdu[] );

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
    @param alpha Coefficient for the Jacobian w.r.t. states
    @param beta Coefficient for the Jacobian w.r.t. first time deriv states
    @param gamma Coefficient for the Jacobian w.r.t. second time deriv states
    @param Jac_nnz Number of non-zero Jacobian entries
    @param Jac_paris The (i,j) locations of the Jacobian entries
    @param Jac The Jacobian values
    @param res The resulting residual
    @param mat The Jacobian matrix
  */
  virtual void addWeakFormJacobian( int n,
                                    const double pt[],
                                    TacsScalar weight,
                                    const TacsScalar J[],
                                    const int vars_per_node,
                                    const TacsScalar DUt[],
                                    const TacsScalar DUx[],
                                    TacsScalar alpha,
                                    TacsScalar beta,
                                    TacsScalar gamma,
                                    const int Jac_nnz,
                                    const int *Jac_pairs,
                                    const TacsScalar *Jac,
                                    TacsScalar *res,
                                    TacsScalar *mat );

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
    computeBasisGradient(pt, N, Nxi);
  }

 protected:
  static TacsScalar computeFaceNormal( const int num_params,
                                       const int num_nodes,
                                       const double Nxi[],
                                       const TacsScalar Xpts[],
                                       const double tangents[],
                                       TacsScalar Xd[],
                                       TacsScalar n[] );
  static TacsScalar computeJacobianTransform( const int num_params,
                                              const int num_nodes,
                                              const double Nxi[],
                                              const TacsScalar Xpts[],
                                              TacsScalar Xd[],
                                              TacsScalar J[] );
  static void addJacobianTransformSens( const int num_params,
                                        const int num_nodes,
                                        const double Nxi[],
                                        const TacsScalar Xd[],
                                        const TacsScalar J[],
                                        TacsScalar dfddetJ,
                                        const TacsScalar dfXd[],
                                        const TacsScalar dfdJ[],
                                        TacsScalar dfdXpts[] );
  static void computeFieldValues( const int num_nodes, const double N[],
                                  const TacsScalar Xpts[],
                                  const int vars_per_node,
                                  const TacsScalar vars[],
                                  TacsScalar X[], TacsScalar U[] );
  static TacsScalar computeFieldGradient( const int num_params,
                                          const int num_nodes,
                                          const double N[], const double Nxi[],
                                          const TacsScalar Xpts[],
                                          const int vars_per_node,
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[],
                                          TacsScalar X[], TacsScalar Xd[],
                                          TacsScalar J[], TacsScalar Ut[],
                                          TacsScalar Ud[], TacsScalar Ux[] );
  static TacsScalar computeFieldGradient( const int num_params,
                                          const int num_nodes,
                                          const double N[], const double Nxi[],
                                          const TacsScalar Xpts[],
                                          const int vars_per_node,
                                          const TacsScalar vars[],
                                          const TacsScalar dvars[],
                                          const TacsScalar ddvars[],
                                          const TacsScalar psi[],
                                          TacsScalar X[], TacsScalar Xd[],
                                          TacsScalar J[], TacsScalar Ut[],
                                          TacsScalar Ud[], TacsScalar Ux[],
                                          TacsScalar Psi[], TacsScalar Psid[],
                                          TacsScalar Psix[] );
  static void addFieldGradientSVSens( const int num_params,
                                      const int num_nodes,
                                      const double N[],
                                      const double Nxi[],
                                      const int vars_per_node,
                                      const TacsScalar J[],
                                      const TacsScalar dfdUt[],
                                      const TacsScalar dfdUx[],
                                      TacsScalar dfdu[] );
  static void addFieldGradientXptSens( const int num_params,
                                       const int num_nodes,
                                       const double N[],
                                       const double Nxi[],
                                       const TacsScalar Xpts[],
                                       const int vars_per_node,
                                       const TacsScalar Xd[],
                                       const TacsScalar J[],
                                       const TacsScalar Ud[],
                                       const TacsScalar dfddetJ,
                                       const TacsScalar dfdX[],
                                       const TacsScalar dfdXd[],
                                       const TacsScalar dfdJ[],
                                       const TacsScalar dfdUx[],
                                       TacsScalar dfdXpts[] );
  static void addWeakFormResidual( const int num_params, const int num_nodes,
                                   const double N[], const double Nxi[],
                                   TacsScalar weight, const TacsScalar J[],
                                   const int vars_per_node,
                                   const TacsScalar DUt[], const TacsScalar DUx[],
                                   TacsScalar res[] );
  static void addWeakFormJacobian( const int num_params, const int num_nodes,
                                   const TacsScalar N[], const TacsScalar Nx[],
                                   TacsScalar weight, const TacsScalar J[],
                                   const int vars_per_node,
                                   const TacsScalar DUt[], const TacsScalar DUx[],
                                   TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                                   const int Jac_nnz, const int *Jac_pairs,
                                   const TacsScalar *Jac,
                                   TacsScalar *res, TacsScalar *mat );
};

#endif // TACS_ELEMENT_BASIS_H
