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
    @param normal The face (or edge) normal
    @return The area contribution
  */
  TacsScalar getFaceNormal( int face, int n,
                            const TacsScalar Xpts[],
                            TacsScalar X[],
                            TacsScalar Xd[],
                            TacsScalar normal[] );

  /**
    Add the derivative of the face normal into the nodal sensitivities

    @param face The face/edge index
    @param n The quadrautre point index
    @param Xpts The node locations
    @param A The area contribution (computed from forward code)
    @param normal The face normal
    @param dfdA The input derivative of the function w.r.t. area
    @param dfdXd The derivative of the function w.r.t. Xd
    @param dfdn The derivative of the function w.r.t. surface normal
    @param dfdXpts The output derivative w.r.t. the node locations
  */
  void addFaceNormalXptSens( int face, int n,
                             const TacsScalar A,
                             const TacsScalar Xd[],
                             const TacsScalar normal[],
                             const TacsScalar dfdA,
                             const TacsScalar dfdX[],
                             const TacsScalar dfdXd[],
                             const TacsScalar dfdn[],
                             TacsScalar dfdXpts[] );

  /**
    Get the Jacobian transformation from computational to physical
    coordinates.

    This code returns the determinant of the transformation and
    returns both the Xd, the derivative of physical coordinates
    w.r.t. the parameters, and J, the inverse of Xd. This code can be
    used for computing the volume/area quadratures.

    @param n The index of the quadrature point
    @param pt The quadrature point
    @param Xpts The node locations
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @return The determinant of Xd
  */
  TacsScalar getJacobianTransform( int n,
                                   const double pt[],
                                   const TacsScalar Xpts[],
                                   TacsScalar Xd[],
                                   TacsScalar J[] );

  /**
    Compute the derivative of the Jacobian transformation

    This code adds the contribution to the derivative of the Jacobian
    transformation to the output dfdXpts. Note that the inputs here
    should be generated by a call to getJacobianTransform first.  Note
    that dfdXd and dfdJ may be passed as NULL.

    @param n The index of the quadrature point
    @param pt The quadrature point
    @param Xpts The node locations
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param dfddetJ The derivative of the function w.r.t. detJ
    @param dfdXd The derivative of the function w.r.t. Xd
    @param dfdJ The derivative of the function w.r.t. J
    @param dfdXpts The output derivative of the function w.r.t. Xpts
  */
  void addJacobianTransformXptSens( int n,
                                    const double pt[],
                                    const TacsScalar Xd[],
                                    const TacsScalar J[],
                                    TacsScalar dfddetJ,
                                    const TacsScalar dfXd[],
                                    const TacsScalar dfdJ[],
                                    TacsScalar dfdXpts[] );

  /**
    Get the field values at the specified quadrature point

    @param n The index of the quadrature point
    @param pt The quadrature point
    @param Xpts The element node locations
    @param vars The element variable values
    @param X The computed coordinate location at quadrature point n
    @param U The computed field values at quadrature point n
  */
  void getFieldValues( int n, const double pt[],
                       const TacsScalar Xpts[],
                       const int vars_per_node,
                       const TacsScalar vars[],
                       TacsScalar X[],
                       TacsScalar U[] );

  /**
    Get the gradient of the field at the quadrature point.

    Note that all matrices are row-major order.

    The arguments Xd and Ud are often only intermediate values, but are returned
    here so that they can be passed in to the differentiated code.

    @param n The index of the quadrature point
    @param pt The parametric location
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param vars The element state variables
    @param dvars The first time derivative of the element state vars
    @param ddvars The second time derivative of the element state vars
    @param X The physical quadrature point
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param Ut The variables and time derivatives at the quadrature point
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param Ux The derivative of the variables w.r.t. the spatial coords
    @return The determinant of the matrix Xd
  */
  TacsScalar getFieldGradient( int n,
                               const double pt[],
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
    Add the derivative of the field gradient terms with respect to
    the state variable values.

    @param n The index of the quadrature point
    @param pt The parametric location
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param dfdUt The derivative of the function w.r.t. Ut
    @param dfdUx The derivative of the function w.r.t. Ux
    @param dfdu The output derivative of the function w.r.t. state variables
  */
  void addFieldGradientSVSens( int n,
                               const double pt[],
                               const TacsScalar Xpts[],
                               const int vars_per_node,
                               const TacsScalar Xd[],
                               const TacsScalar J[],
                               const TacsScalar Ud[],
                               const TacsScalar dfdUt[],
                               TacsScalar dfdUx[],
                               TacsScalar dfdu[] );

  /**
    Add the contributions to the derivative of the node locations from
    the computation of the field gradient.

    This member function adds the contribution to the gradient from derivatives
    of a function with respect to the determinant, the node location, the
    derivatives of the nodes with respect to the parametric coordinates, and
    the terms in the Jacobian transformation, and the derivatives of the
    field quantities.

    The terms dfdX, dfdXd, dfdJ and dfdUx may be passed in as NULL if they
    are zero.

    @param n The index of the quadrature point
    @param pt The parametric location
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param dfddetJ The derivative of the determinant
    @param dfdX The derivative w.r.t. node location
    @param dfdXd The derivative w.r.t. the components of Xd
    @param dfdJ The derivative w.r.t. the components of J
    @param dfdUx The derivative w.r.t. the components of Ux
    @param dfdXpts The output derivative of the function w.r.t. node locations
  */
  void addFieldGradientXptSens( int n,
                                const double pt[],
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
    Get the gradient of the field at the quadrature point.

    This function returns the values and derivatives of the state variables,
    and a corresponding adjoint vector. This is used in computing the
    derivative of the adjoint-residual product.

    Note that the adjoint variable vector Psi is of length vars_per_node,
    (not 3*var_per_node), since it only contains the adjoint values, and
    not their first and second time derivatives.

    @param n The index of the quadrature point
    @param pt The parametric location
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param vars The element state variables
    @param dvars The first time derivative of the element state vars
    @param ddvars The second time derivative of the element state vars
    @param X The physical quadrature point
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param Ut The variables and time derivatives at the quadrature point
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param Ux The derivative of the variables w.r.t. the spatial coords
    @param Psi The adjoint variable values (no time derivatives!)
    @param Psid The derivatives of the adjoint variables w.r.t. parameters
    @param Psix The spatial derivatives of the adjoint variables
    @return The determinant of the matrix Xd
  */
  TacsScalar getFieldGradient( int n,
                               const double pt[],
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
    Add the contributions to the derivative of the node locations from
    the computation of the field and adjoint gradients.

    The terms dfdX, dfdXd, dfdJ, dfdPsix or dfdUx may be passed in as NULL if they
    are zero. Note that if either dfdPsix or dfdUx is NULL, then the other must
    be as well. (If not you should use the other version of this code.)

    @param n The index of the quadrature point
    @param pt The parametric location
    @param Xpts The element node locations
    @param vars_per_node The number of degrees of freedom per node
    @param Xd The derivative of the physical node location w.r.t. parameters
    @param J The Jacobian transformation (inverse of Xd)
    @param Ud The derivative of the variables w.r.t. the parametric coords
    @param Psid The derivatives of the adjoint variables w.r.t. parameters
    @param dfddetJ The derivative of the determinant
    @param dfdX The derivative w.r.t. node location
    @param dfdXd The derivative w.r.t. the components of Xd
    @param dfdJ The derivative w.r.t. the components of J
    @param dfdUx The derivative w.r.t. the components of Ux
    @param dfdPsix The derivative w.r.t. the components of Psix
    @param dfdXpts The output derivative of the function w.r.t. node locations
  */
  void addFieldGradientXptSens( int n,
                                const double pt[],
                                const TacsScalar Xpts[],
                                const int vars_per_node,
                                const TacsScalar Xd[],
                                const TacsScalar J[],
                                const TacsScalar Ud[],
                                const TacsScalar Psid[],
                                const TacsScalar dfddetJ,
                                const TacsScalar dfdX[],
                                const TacsScalar dfdXd[],
                                const TacsScalar dfdJ[],
                                const TacsScalar dfdUx[],
                                const TacsScalar dfdPsix[],
                                TacsScalar dfdXpts[] );

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
  void addWeakResidual( int n, const double pt[],
                        TacsScalar weight,
                        const TacsScalar J[],
                        const int vars_per_node,
                        TacsScalar DUt[],
                        TacsScalar DUx[],
                        TacsScalar res[] );

  /**
    Scale the terms in the Jacobian matrix

    @param weight The quadrature weight
    @param alpha Coefficient for the Jacobian w.r.t. states
    @param beta Coefficient for the Jacobian w.r.t. first time deriv states
    @param gamma Coefficient for the Jacobian w.r.t. second time deriv states
    @param Jac_nnz Number of non-zero Jacobian entries
    @param Jac_paris The (i,j) locations of the Jacobian entries
    @param Jac The Jacobian values
  */
  void scaleWeakMatrix( const TacsScalar weight,
                        const TacsScalar alpha,
                        const TacsScalar beta,
                        const TacsScalar gamma,
                        const int Jac_nnz,
                        const int *Jac_pairs,
                        TacsScalar *Jac );

  /**
    Add the entries from the matrix of the weak form residual to a matrix

    @param n The quadrautre point index
    @param pt The quadrature point value
    @param J The Jacobian coordinate transformation
    @param vars_per_node The number of variables per node
    @param Jac_nnz Number of non-zero Jacobian entries
    @param Jac_paris The (i,j) locations of the Jacobian entries
    @param Jac The Jacobian values
    @param mat The Jacobian matrix
  */
  void addWeakMatrix( int n, const double pt[],
                      const TacsScalar J[],
                      const int vars_per_node,
                      const int Jac_nnz,
                      const int *Jac_pairs,
                      const TacsScalar *Jac,
                      TacsScalar *mat );

  /**
    Add the entries from the Jacobian of the weak form residual to a matrix

    @param n The quadrautre point index
    @param pt The quadrature point value
    @param J The Jacobian coordinate transformation
    @param vars_per_node The number of variables per node
    @param Jac_nnz Number of non-zero Jacobian entries
    @param Jac_paris The (i,j) locations of the Jacobian entries
    @param Jac The Jacobian values
    @param temp A temporary array of at least size (2*getNumParameters()+1)*vars_per_node
    @param px The input vector
    @param py The output vector
  */
  void addMatVecProduct( int n, const double pt[],
                         const TacsScalar J[],
                         const int vars_per_node,
                         const int Jac_nnz,
                         const int *Jac_pairs,
                         const TacsScalar *Jac,
                         TacsScalar temp[],
                         const TacsScalar *px,
                         TacsScalar *py );

  // Temporary before wrecking everything...
  virtual void computeBasis( const double pt[], double N[] ) = 0;
  virtual void computeBasisGradient( const double pt[], double N[], double Nxi[] ) = 0;

  /**
    Interpolate the specified number of fields

    This function computes the following for i = 0, num_fields-1

    field[incr*i] = sum_{j} N[j]*values[num_fields*j + i]

    @param n The quadrature point index
    @param pt The parametric point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param incr The increment between locations in the field array
    @param field The field values
  */
  virtual void interpFields( const int n,
                             const double pt[],
                             const int num_fields,
                             const TacsScalar values[],
                             const int incr,
                             TacsScalar field[] ){
    const int num_nodes = getNumNodes();
    double N[256];
    computeBasis(pt, N);
    for ( int i = 0; i < num_fields; i++ ){
      field[incr*i] = 0.0;
      for ( int j = 0; j < num_nodes; j++ ){
        field[incr*i] += values[num_fields*j + i]*N[j];
      }
    }
  }

  /**
    Compute the interpolation field for three different interpolants
    simultaneously. This is common when assemblying the temporal derivatives

    This function computes the following for i = 0, num_fields-1

    field[3*i] = sum_{j} N[j]*vals1[num_fields*j + i]
    field[3*i+1] = sum_{j} N[j]*vals2[num_fields*j + i]
    field[3*i+2] = sum_{j} N[j]*vals3[num_fields*j + i]

    @param n The quadrature point index
    @param pt The parametric point
    @param num_fields The number of fields to interpolate
    @param vals1 The first array of values at the nodes
    @param vals2 The second array of values at the nodes
    @param vals3 The third array of values at the nodes
    @param field The field values
  */
  virtual void interpFields( const int n,
                             const double pt[],
                             const int num_fields,
                             const TacsScalar val1[],
                             const TacsScalar val2[],
                             const TacsScalar val3[],
                             TacsScalar field[] ){
    interpFields(n, pt, num_fields, val1, 3, &field[0]);
    interpFields(n, pt, num_fields, val2, 3, &field[1]);
    interpFields(n, pt, num_fields, val3, 3, &field[2]);
  }

  /**
    Compute the interpolate to a quadrature point on the face

    By default, this evaluates the interpFields function, but if you have
    a specific shape function implementation, indexed on the quadrature
    point, this will not work.

    @param face The face index
    @param n The quadrature point index on this face
    @param num_fields The number of fields to interpolate
    @param values The values to interpolate
    @param incr The increment between locations in the field array
    @param field The field values
  */
  virtual void interpFaceFields( const int face,
                                 const int n,
                                 const double pt[],
                                 const int num_fields,
                                 const TacsScalar values[],
                                 const int incr,
                                 TacsScalar field[] ){
    interpFields(-1, pt, num_fields, values, incr, field);
  }

  /**
    Add the transpose of the interpolation operation to the vector

    This function computes the following for i = 0, num_fields-1,
    and j = 0, num_nodes-1

    values[num_fields*j + i] += N[j]*field[incr*i]

    @param n The quadrature point index
    @param pt The parametric point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param incr The increment between locations in the field array
    @param field The field values
  */
  virtual void addInterpFieldsTranspose( const int n,
                                         const double pt[],
                                         const int incr,
                                         const TacsScalar field[],
                                         const int num_fields,
                                         TacsScalar values[] ){
    const int num_nodes = getNumNodes();
    double N[256];
    computeBasis(pt, N);
    for ( int i = 0; i < num_fields; i++ ){
      for ( int j = 0; j < num_nodes; j++ ){
        values[num_fields*j + i] += field[incr*i]*N[j];
      }
    }
  }

  /**
    Add the transpose of the interpolation operation to the vector
    on the face quadrature points.

    @param n The quadrature point index
    @param pt The parametric point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param incr The increment between locations in the field array
    @param field The field values
  */
  virtual void addInterpFaceFieldsTranspose( const int face,
                                             const int n,
                                             const double pt[],
                                             const int incr,
                                             const TacsScalar field[],
                                             const int num_fields,
                                             TacsScalar values[] ){
    addInterpFieldsTranspose(-1, pt, incr, field, num_fields, values);
  }

  /**
    Compute the gradient of the fields in the computational space

    This function must compute

    grad[num_params*i + j] = sum_{k} N_{k,j}*values[num_fields*k + i]

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param grad The gradient of the field in the computational space
  */
  virtual void interpFieldsGrad( const int n,
                                 const double pt[],
                                 const int num_fields,
                                 const TacsScalar values[],
                                 TacsScalar grad[] ){
    const int num_nodes = getNumNodes();
    const int num_params = getNumParameters();
    double N[256], Nxi[3*256];
    computeBasisGradient(pt, N, Nxi);
    for ( int i = 0; i < num_fields; i++ ){
      for ( int j = 0; j < num_params; j++ ){
        grad[num_params*i + j] = 0.0;
      }
      for ( int k = 0; k < num_nodes; k++ ){
        for ( int j = 0; j < num_params; j++ ){
          grad[num_params*i + j] += values[num_fields*k + i]*Nxi[num_params*k + j];
        }
      }
    }
  }

  /**
    Compute the interpolate to a quadrature point on the face

    By default, this evaluates the interpFields function, but if you have
    a specific shape function implementation, indexed on the quadrature
    point, this will not work.

    @param face The face index
    @param n The quadrature point index on this face
    @param num_fields The number of fields to interpolate
    @param values The values to interpolate
    @param grad The gradient of the field in the computational space
  */
  virtual void interpFaceFieldsGrad( const int face,
                                     const int n,
                                     const double pt[],
                                     const int num_fields,
                                     const TacsScalar values[],
                                     TacsScalar grad[] ){
    interpFieldsGrad(-1, pt, num_fields, values, grad);
  }

  /**
    Add the transpose of the gradient interpolation to the vector

    This function computes the following for i = 0, num_fields-1,
    j = 1, num_params-1, and j = 0, num_nodes-1

    values[num_fields*k + i] += N_{k,j}*grad[incr*i + j]

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param grad The gradient of the field in the computational space
  */
  virtual void addInterpFieldsGradTranspose( int n,
                                             const double pt[],
                                             const int num_fields,
                                             const TacsScalar grad[],
                                             TacsScalar values[] ){
    const int num_nodes = getNumNodes();
    const int num_params = getNumParameters();
    double N[256], Nxi[3*256];
    computeBasisGradient(pt, N, Nxi);
    for ( int i = 0; i < num_fields; i++ ){
      for ( int k = 0; k < num_nodes; k++ ){
        for ( int j = 0; j < num_params; j++ ){
          values[num_fields*k + i] += grad[num_params*i + j]*Nxi[num_params*k + j];
        }
      }
    }
  }

  /**
    Add the transpose of the gradient interpolation to the vector
    at quadrature points on the specified face

    @param face The face index
    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param num_fields The number of fields to interpolate
    @param values The values of the interpolant at the nodes
    @param grad The gradient of the field in the computational space
  */
  virtual void addInterpFaceFieldsGradTranspose( const int face,
                                                 int n,
                                                 const double pt[],
                                                 const int num_fields,
                                                 const TacsScalar grad[],
                                                 TacsScalar values[] ){
    addInterpFieldsGradTranspose(-1, pt, num_fields, grad, values);
  }

  /*
    Add the outer-product of the shape functions to the matrix

    mat[row_incr*i + col_incr*j] += scale*N[i]*N[j]

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param weight The weight factor added to the matrix
    @param row_incr The row increment applied to the matrix
    @param col_incr The column increment applied to the matrix
    @param mat The element matrix
  */
  virtual void addInterpOuterProduct( const int n,
                                      const double pt[],
                                      const TacsScalar weight,
                                      const int row_incr,
                                      const int col_incr,
                                      TacsScalar *mat ){
    const int num_nodes = getNumNodes();
    double N[256];
    computeBasis(pt, N);
    for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
      for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
        mat[0] += weight*N[i]*N[j];
      }
    }
  }

  /*
    Add the outer-product of the shape functions to the matrix

    mat[row_incr*i + col_incr*j] += scale*N[i]*N[j]

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param weight The weight factor added to the matrix
    @param row_incr The row increment applied to the matrix
    @param col_incr The column increment applied to the matrix
    @param mat The element matrix
  */
  virtual void addInterpGradOuterProduct( const int n,
                                          const double pt[],
                                          const int transpose,
                                          const TacsScalar weight,
                                          const TacsScalar scale[],
                                          const int row_incr,
                                          const int col_incr,
                                          TacsScalar *mat ){
    const int num_nodes = getNumNodes();
    const int num_params = getNumParameters();
    double N[256], Nxi[3*256];
    computeBasisGradient(pt, N, Nxi);

    if (transpose){
      if (num_params == 1){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] +=
              weight*N[j]*(Nxi[i]*scale[0]);
          }
        }
      }
      else if (num_params == 2){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] += weight*N[j]*(Nxi[2*i]*scale[0] + Nxi[2*i+1]*scale[1]);
          }
        }
      }
      else if (num_params == 3){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] += weight*N[j]*(Nxi[3*i]*scale[0] +
                                   Nxi[3*i+1]*scale[1] +
                                   Nxi[3*i+2]*scale[2]);
          }
        }
      }
    }
    else {
      if (num_params == 1){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] +=
              weight*N[i]*(Nxi[j]*scale[0]);
          }
        }
      }
      else if (num_params == 2){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] += weight*N[i]*(Nxi[2*j]*scale[0] + Nxi[2*j+1]*scale[1]);
          }
        }
      }
      else if (num_params == 3){
        for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
          for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
            mat[0] += weight*N[i]*(Nxi[3*j]*scale[0] +
                                   Nxi[3*j+1]*scale[1] +
                                   Nxi[3*j+2]*scale[2]);
          }
        }
      }
    }
  }

  /*
    Add the outer-product of the shape functions to the matrix

    mat[row_incr*i + col_incr*j] += scale*N[i]*N[j]

    @param n The quadrature point index
    @param pt The parametric location of the quadrature point
    @param weight The weight factor added to the matrix
    @param row_incr The row increment applied to the matrix
    @param col_incr The column increment applied to the matrix
    @param mat The element matrix
  */
  virtual void addInterpGradGradOuterProduct( const int n,
                                              const double pt[],
                                              const TacsScalar weight,
                                              const TacsScalar iscale[],
                                              const TacsScalar jscale[],
                                              const int row_incr,
                                              const int col_incr,
                                              TacsScalar *mat ){
    const int num_nodes = getNumNodes();
    const int num_params = getNumParameters();
    double N[256], Nxi[3*256];
    computeBasisGradient(pt, N, Nxi);

    if (num_params == 1){
      for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
        for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
          mat[0] +=
            weight*(Nxi[i]*iscale[0])*(Nxi[j]*jscale[0]);
        }
      }
    }
    else if (num_params == 2){
      for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
        for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
          mat[0] +=
            weight*(Nxi[2*i]*iscale[0] + Nxi[2*i+1]*iscale[1])*
                  (Nxi[2*j]*jscale[0] + Nxi[2*j+1]*jscale[1]);
        }
      }
    }
    else if (num_params == 3){
      for ( int i = 0; i < num_nodes; i++, mat += row_incr ){
        for ( int j = 0; j < num_nodes; j++, mat += col_incr ){
          mat[0] +=
            weight*(Nxi[3*i]*iscale[0] + Nxi[3*i+1]*iscale[1] + Nxi[3*i+2]*iscale[2])*
                   (Nxi[3*j]*jscale[0] + Nxi[3*j+1]*jscale[1] + Nxi[3*j+2]*jscale[2]);
        }
      }
    }
  }
};

#endif // TACS_ELEMENT_BASIS_H
