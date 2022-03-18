#ifndef TACS_BEAM_ELEMENT_H
#define TACS_BEAM_ELEMENT_H

#include "TACSBeamElementModel.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamUtilities.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSBeamConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"
#include "a2d.h"

/*
  Compute the transformation from the local coordinates
*/
class TACSBeamTransform : public TACSObject {
 public:
  /*
    Given the local beam element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ) = 0;
  virtual void computeTransformSens( const TacsScalar X0xi_vals[],
                                     const TacsScalar dTvals[],
                                     TacsScalar dX0xi[] ) = 0;
  virtual A2D::Vec3& getRefAxis() = 0;
};

/*
  Compute the transformation
*/
class TACSBeamRefAxisTransform : public TACSBeamTransform {
 public:
  TACSBeamRefAxisTransform( const TacsScalar axis_dir[] ){
    A2D::Vec3 axdir(axis_dir);
    A2D::Vec3Normalize normalize(axdir, axis);
  }

  void computeTransform( const TacsScalar X0xi_vals[],
                         TacsScalar Tvals[] ){
    // Normalize the first direction.
    A2D::Vec3 X0xi(X0xi_vals);
    A2D::Vec3 t1;
    A2D::Vec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::Vec3 t2_dir;
    A2D::Scalar dot;
    A2D::Vec3Dot dott1(axis, t1, dot);
    A2D::Vec3Axpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Compute the t2 direction
    A2D::Vec3 t2;
    A2D::Vec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::Vec3 t3;
    A2D::Vec3CrossProduct cross(t1, t2, t3);

    // Assemble the referece frame
    A2D::Mat3x3 T;
    A2D::Mat3x3FromThreeVec3 assembleT(t1, t2, t3, T);

    for ( int i = 0; i < 9; i++ ){
      Tvals[i] = T.A[i];
    }
  }

  void computeTransformSens( const TacsScalar X0xi_vals[],
                             const TacsScalar dTvals[],
                             TacsScalar dX0xi[] ){
    // Normalize the first direction.
    A2D::ADVec3 X0xi(X0xi_vals);
    A2D::ADVec3 t1;
    A2D::ADVec3Normalize normalizet1(X0xi, t1);

    // t2_dir = axis - dot(t1, axis) * t1
    A2D::ADVec3 t2_dir;
    A2D::ADScalar dot;
    A2D::Vec3ADVecDot dott1(axis, t1, dot);
    A2D::ADVec3VecAxpy axpy(-1.0, dot, t1, axis, t2_dir);

    // Compute the t2 direction
    A2D::ADVec3 t2;
    A2D::ADVec3Normalize normalizet2(t2_dir, t2);

    // Compute the n2 direction
    A2D::ADVec3 t3;
    A2D::ADVec3CrossProduct cross(t1, t2, t3);

    // Assemble the referece frame
    A2D::ADMat3x3 T(NULL, dTvals); // Set the seeds for T
    A2D::ADMat3x3FromThreeADVec3 assembleT(t1, t2, t3, T);

    // Reverse the operations to get the derivative w.r.t. X0
    assembleT.reverse();
    cross.reverse();
    normalizet2.reverse();
    axpy.reverse();
    dott1.reverse();
    normalizet1.reverse();

    for ( int i = 0; i < 3; i++ ){
      dX0xi[i] = X0xi.xd[i];
    }
  }
  A2D::Vec3& getRefAxis(){
    return axis;
  }

 private:
  A2D::Vec3 axis;
};

/*
  The position in the beam is parametrized using the coordinates d =
  (xi, z1, z2) and the n1(xi) and n2(xi) coordinates that are
  orthogonal to the tangent of the beam at the nodes. The normalize
  directions n1(xi) are interpolated from the nodes using the same
  interpolation scheme for x0(xi) and are therefore not normal at
  quadrature points. The position is given as

  x = x0(xi) + z1 * n1(xi) + z2 * n2(xi)

  The derivative of the position w.r.t. the parameters is

  x,d = [x0,xi + z1 * n1,xi + z2 * n2,xi | n1 | n2 ]

  Imporantatly, the curvature of the beam leads to the second
  derivatives

  x,d,z1 = n1,xi * e1^{T}
  x,d,z2 = n2,xi * e1^{T}

  To compute the derivatives in the global reference frame, it is
  necessary to apply a Jacobian transformation using the inverse
  matrix

  x,d0^{-1} = [x0,xi | n1 | n2 ]^{-1}

  Additionally, we require the derivative of this transformation through
  the thickness that takes the form

  x,d0^{-1},z1 = - x,d0^{-1} * x,d,z1 * x,d0^{-1}
  x,d0^{-1},z2 = - x,d0^{-1} * x,d,z2 * x,d0^{-1}

  The displacement of the beam is parametrized based on the reference
  line displacement u0(xi) and the directors d1(xi) and d2(xi),
  parametrized in a consistent manner with the shell element.  The
  displacement field is given as

  u = u0(xi) + z1 * d1(xi) + z2 * d2(xi)

  The derivative of the displacement field w.r.t. the natural parameters is

  u,d = [u0,xi + z1 * n1,xi + z2 * n2,xi | d1 | d2 ]

  Again, higher-order derivatives are required

  u,d,z1 = [d1,xi | 0 | 0 ] = d1,xi * e1^{T}
  u,d,z2 = [d2,xi | 0 | 0 ] = d2,xi * e1^{T}

  Now, the displacement gradient in a beam-attached local coordinate
  frame can be found, where T is a transformation from the local to
  global reference frame

  u,x = T^{T} * u,d * [ x,d ]^{-1} * T

  However, both u,d and x,d depend on the through-beam coordinates z1
  and z2. As a result, we use the approximation

  u,x ~= u0x + z1 * d1x * e1^{T} + z2 * d2x * e1^{T} + O(z**2)

  The term u0x is

  u0x = T^{T} * [u0,xi | d1 | d2] * [x0,xi | n1 | n2]^{-1} * T

  The term d1x is

  d1x = T^{T} * (d1,xi * e1^{T} * x,d0^{-1} +
  .              u0,xi * e1^{T} * x,d0^{-1},z1 ) * T^{T} * e1


  Note that the second coefficient for d1x can be simplifed to

  e1^{T} * x,d0^{-1},z1 * T^{T} * e1
  .   = - e1^{T} * x,d0^{-1} * x,d,z1 * x,d0^{-1} * T^{T} * e1
  .   = - (e1^{T} * x,d0^{-1} * n1,xi) * (e1^{T} * x,d0^{-1} * T^{T} * e1)

  The term d2x is

  d2x = T^{T} * (d2,xi * e1^{T} * x,d0^{-1} +
  .              u0,xi * e1^{T} * x,d0^{-1},z2 ) * T^{T} * e1

  These expressions can be simplifed to the following:

  d1x = s0 * T^{T} * (d1,xi - sz1 * u0,xi)
  d2x = s0 * T^{T} * (d2,xi - sz2 * u0,xi)

  where:

  s0 = e1^{T} * x,d0^{-1} * T^{T} * e1
  sz1 = e1^{T} * x,d0^{-1} * n1,xi
  sz2 = e1^{T} * x,d0^{-1} * n2,xi
*/
template <class quadrature, class basis, class director, class model>
class TACSBeamElement : public TACSElement {
 public:
  // Offset within the solution vector to the roational
  // parametrization defined via the director class. Here the offset
  // is 3 corresponding to the (u, v, w) displacements of the beam.
  static const int offset = 3;

  // The number of variables defined at each node of the shell
  // element.  There are 3 mid-surface displacements plus however many
  // parameters are defined by the director class for the
  // parametrization.
  static const int vars_per_node = offset + director::NUM_PARAMETERS;

  // The number of nodes for this element. This is derived from the
  // basis function class. This is just a handy re-definition since
  // this constant is used in many locations within the element.
  static const int num_nodes = basis::NUM_NODES;

  TACSBeamElement( TACSBeamTransform *_transform,
                   TACSBeamConstitutive *_con ){
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  int getVarsPerNode(){
    return vars_per_node;
  }
  int getNumNodes(){
    return basis::NUM_NODES;
  }

  ElementLayout getLayoutType(){
    return basis::getLayoutType();
  }

  int getNumQuadraturePoints(){
    return quadrature::getNumQuadraturePoints();
  }

  double getQuadratureWeight( int n ){
    return quadrature::getQuadratureWeight(n);
  }

  double getQuadraturePoint( int n, double pt[] ){
    return quadrature::getQuadraturePoint(n, pt);
  }

  int getNumElementFaces(){
    return quadrature::getNumFaces();
  }

  int getNumFaceQuadraturePoints( int face ){
    return quadrature::getNumFaceQuadraturePoints(face);
  }

  double getFaceQuadraturePoint( int face, int n, double pt[],
                                 double tangent[] ){
    return quadrature::getFaceQuadraturePoint(face, n, pt, tangent);
  }

  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] ){
    return con->getDesignVarNums(elemIndex, dvLen, dvNums);
  }

  int setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] ){
    return con->setDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] ){
    return con->getDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVarRange( int elemIndex, int dvLen, TacsScalar lb[], TacsScalar ub[] ){
    return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
  }

  void computeEnergies( int elemIndex,
                        double time,
                        const TacsScalar Xpts[],
                        const TacsScalar vars[],
                        const TacsScalar dvars[],
                        TacsScalar *Te,
                        TacsScalar *Pe );

  void addResidual( int elemIndex,
                    double time,
                    const TacsScalar *Xpts,
                    const TacsScalar *vars,
                    const TacsScalar *dvars,
                    const TacsScalar *ddvars,
                    TacsScalar *res );

  void addAdjResXptProduct( int elemIndex,
                            double time,
                            TacsScalar scale,
                            const TacsScalar psi[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[],
                            TacsScalar fXptSens[] );

  // void addJacobian( int elemIndex, double time,
  //                   TacsScalar alpha,
  //                   TacsScalar beta,
  //                   TacsScalar gamma,
  //                   const TacsScalar Xpts[],
  //                   const TacsScalar vars[],
  //                   const TacsScalar dvars[],
  //                   const TacsScalar ddvars[],
  //                   TacsScalar res[],
  //                   TacsScalar mat[] );

  // void addAdjResProduct( int elemIndex, double time,
  //                        TacsScalar scale,
  //                        const TacsScalar psi[],
  //                        const TacsScalar Xpts[],
  //                        const TacsScalar vars[],
  //                        const TacsScalar dvars[],
  //                        const TacsScalar ddvars[],
  //                        int dvLen,
  //                        TacsScalar dfdx[] );

  // int evalPointQuantity( int elemIndex, int quantityType,
  //                        double time,
  //                        int n, double pt[],
  //                        const TacsScalar Xpts[],
  //                        const TacsScalar vars[],
  //                        const TacsScalar dvars[],
  //                        const TacsScalar ddvars[],
  //                        TacsScalar *detXd,
  //                        TacsScalar *quantity );

  // void addPointQuantityDVSens( int elemIndex, int quantityType,
  //                              double time,
  //                              TacsScalar scale,
  //                              int n, double pt[],
  //                              const TacsScalar Xpts[],
  //                              const TacsScalar vars[],
  //                              const TacsScalar dvars[],
  //                              const TacsScalar ddvars[],
  //                              const TacsScalar dfdq[],
  //                              int dvLen,
  //                              TacsScalar dfdx[] );

  // void addPointQuantitySVSens( int elemIndex, int quantityType,
  //                              double time,
  //                              TacsScalar alpha,
  //                              TacsScalar beta,
  //                              TacsScalar gamma,
  //                              int n, double pt[],
  //                              const TacsScalar Xpts[],
  //                              const TacsScalar vars[],
  //                              const TacsScalar dvars[],
  //                              const TacsScalar ddvars[],
  //                              const TacsScalar dfdq[],
  //                              TacsScalar dfdu[] );

  void getOutputData( int elemIndex,
                      ElementType etype,
                      int write_flag,
                      const TacsScalar Xpts[],
                      const TacsScalar vars[],
                      const TacsScalar dvars[],
                      const TacsScalar ddvars[],
                      int ld_data,
                      TacsScalar *data );

 private:
  // Set sizes for the different components
  static const int usize = 3*basis::NUM_NODES;
  static const int dsize = 3*basis::NUM_NODES;
  static const int csize = 9*basis::NUM_NODES;

  TACSBeamTransform *transform;
  TACSBeamConstitutive *con;
};


/*
  Compute the kinetic and potential energies of the shell
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
  computeEnergies( int elemIndex,
                   double time,
                   const TacsScalar *Xpts,
                   const TacsScalar *vars,
                   const TacsScalar *dvars,
                   TacsScalar *Telem, TacsScalar *Uelem ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec3& axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3*basis::NUM_NODES], fn2[3*basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template
    computeDirectorRates<vars_per_node, offset,
                         basis::NUM_NODES>(vars, dvars, fn1, d1, d1dot);
  director::template
    computeDirectorRates<vars_per_node, offset,
                         basis::NUM_NODES>(vars, dvars, fn2, d2, d2dot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template
    computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars, d1, d2, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::Mat3x3 T;

    // Parametric location
    A2D::Vec3 X0;

    // Tangent to the beam
    A2D::Vec3 X0xi;

    // Interpolated normal directions
    A2D::Vec3 n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec3 n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::Vec3 u0xi, d01, d02, d01xi, d02xi;

    // Compute X, X,xi and the interpolated normal
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Compute the transformation at the quadrature point
    transform->computeTransform(X0xi.x, T.A);

    // Compute the inverse
    A2D::Mat3x3 Xd, Xdinv;
    A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::Mat3x3Inverse invXd(Xd, Xdinv);

    // Compute the determinant of the transform
    A2D::Scalar detXd;
    A2D::Mat3x3Det computedetXd(weight, Xd, detXd);

    // Compute XdinvT = Xdinv * T
    A2D::Mat3x3 XdinvT;
    A2D::Mat3x3MatMult multXinvT(Xdinv, T, XdinvT);

    // Assemble u0d
    A2D::Mat3x3 u0d;
    A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::Mat3x3 u0dXdinvT, u0x;
    A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

    // Compute s0, sz1 and sz2
    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1({1.0, 0.0, 0.0});
    A2D::Mat3x3InnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3InnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3InnerProduct innersz1(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::Vec3 d1t, d1x;
    A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d1xi, d1t);
    A2D::MatTrans3x3MultScale matmult2dx(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::Vec3 d2t, d2x;
    A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d2xi, d2t);
    A2D::MatTrans3x3MultScale matmultd2x(s0, T, d2t, d2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2]; // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2];
    e0ty[0] = gty[0];
    e0ty[1] = gty[1];

    // Compute the set of strain components
    TacsScalar e[6]; // The components of the strain
    model::evalStrain(u0x.A, u1x.A, u2x.A, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    Ue += 0.5 * detXd.value * (s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                               s[3]*e[3] + s[4]*e[4] + s[5]*e[5]);
  }

  *Telem = Te;
  *Uelem = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
  addResidual( int elemIndex,
               double time,
               const TacsScalar *Xpts,
               const TacsScalar *vars,
               const TacsScalar *dvars,
               const TacsScalar *ddvars,
               TacsScalar *res ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec3& axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3*basis::NUM_NODES], fn2[3*basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1d[dsize], d2d[dsize];
  TacsScalar d1dot[dsize], d2dot[dsize];
  TacsScalar d1ddot[dsize], d2ddot[dsize];
  memset(d1d, 0, dsize*sizeof(TacsScalar));
  memset(d2d, 0, dsize*sizeof(TacsScalar));
  director::template
    computeDirectorRates<vars_per_node, offset,
                         basis::NUM_NODES>(vars, dvars, ddvars, fn1, d1, d1dot, d1ddot);
  director::template
    computeDirectorRates<vars_per_node, offset,
                         basis::NUM_NODES>(vars, dvars, ddvars, fn2, d2, d2dot, d2ddot);

  // Compute the tying strain values
  TacsScalar ety[basis::NUM_TYING_POINTS], dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));
  model::template
    computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars, d1, d2, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::Mat3x3 T;

    // Parametric location
    A2D::Vec3 X0;

    // Tangent to the beam
    A2D::Vec3 X0xi;

    // Interpolated normal directions
    A2D::Vec3 n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec3 n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::ADVec3 u0xi, d01, d02, d01xi, d02xi;

    // Compute X, X,xi and the interpolated normal
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Compute the transformation at the quadrature point
    transform->computeTransform(X0xi.x, T.A);

    // Compute the inverse
    A2D::Mat3x3 Xd, Xdinv;
    A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::Mat3x3Inverse invXd(Xd, Xdinv);

    // Compute the determinant of the transform
    A2D::Scalar detXd;
    A2D::Mat3x3Det computedetXd(weight, Xd, detXd);

    // Compute XdinvT = Xdinv * T
    A2D::Mat3x3 XdinvT;
    A2D::Mat3x3MatMult multXinvT(Xdinv, T, XdinvT);

    // Assemble the matrix Xdz1 = [n1,xi | 0 | 0] and Xdz2 = [n2,xi | 0 | 0 ]
    A2D::Mat3x3 Xdz1, Xdz2;
    A2D::Mat3x3FromVec3 assembleXdz1(n1xi, Xdz1);
    A2D::Mat3x3FromVec3 assembleXdz2(n2xi, Xdz2);

    // Compute Xdinvz1T = - Xdinv * Xdz1 * Xdinv * T
    A2D::Mat3x3 Xdinvz1T, Xdz1XdinvT;
    A2D::Mat3x3MatMult multXdz1XdinvT(Xdz1, XdinvT, Xdz1XdinvT);
    A2D::Mat3x3MatMult multXdinvz1T(-1.0, Xdinv, Xdz1XdinvT, Xdinvz1T);

    // Compute Xdinvz2T = - Xdinv * Xdz2 * Xdinv * T
    A2D::Mat3x3 Xdinvz2T, Xdz2XdinvT;
    A2D::Mat3x3MatMult multXdz2XdinvT(Xdz2, XdinvT, Xdz2XdinvT);
    A2D::Mat3x3MatMult multXdinvz2T(-1.0, Xdinv, Xdz2XdinvT, Xdinvz2T);

    // Assemble u0d, u1d and u2d
    A2D::ADMat3x3 u0d, u1d, u2d;
    A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);
    A2D::ADMat3x3FromADVec3 assembleu1d(d01xi, u1d);
    A2D::ADMat3x3FromADVec3 assembleu2d(d02xi, u2d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::ADMat3x3 u0dXdinvT, u0x;
    A2D::ADMat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

    // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
    A2D::ADMat3x3 u1dXdinvT, u1x;
    A2D::ADMat3x3MatMult multu1d(u1d, XdinvT, u1dXdinvT);
    A2D::ADMat3x3MatMultAdd multu1dadd(u0d, Xdinvz1T, u1dXdinvT);
    A2D::MatTrans3x3ADMatMult multu1x(T, u1dXdinvT, u1x);

    // Compute u2x = T^{T} * (u2d * XdinvT + u0d * XdinvzT)
    A2D::ADMat3x3 u2dXdinvT, u2x;
    A2D::ADMat3x3MatMult multu2d(u2d, XdinvT, u2dXdinvT);
    A2D::ADMat3x3MatMultAdd multu2dadd(u0d, Xdinvz2T, u2dXdinvT);
    A2D::MatTrans3x3ADMatMult multu2x(T, u2dXdinvT, u2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2]; // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], de0ty[2];
    e0ty[0] = gty[0];
    e0ty[1] = gty[1];

    // Evaluate the strain
    TacsScalar e[6];
    model::evalStrain(u0x.A, u1x.A, u2x.A, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    // Evaluate the strain and strain derivatives from the
    model::evalStrainSens(detXd.value, s, u0x.A, u1x.A, u2x.A, e0ty,
                          u0x.Ad, u1x.Ad, u2x.Ad, de0ty);

    // Reverse the operations for the derivative w.r.t. state variables
    multu2x.reverse();
    multu2dadd.reverse();
    multu2d.reverse();
    multu1x.reverse();
    multu1dadd.reverse();
    multu1d.reverse();
    multu0x.reverse();
    multu0d.reverse();
    assembleu2d.reverse();
    assembleu1d.reverse();
    assembleu0d.reverse();

    // Add the residual contributions back to the element
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.xd, res);

    // Add the constributions back to the derivative
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, d2d);

    // Add the contributions to the tying strain
    TacsScalar dgty[2];
    dgty[0] = de0ty[0];
    dgty[1] = de0ty[1];

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);
  }

  // Add the contributions from the tying strain
  model::template
    addComputeTyingStrainTranspose<vars_per_node, basis>(Xpts, fn1, fn2, vars, d1, d2,
                                                         dety, res, d1d, d2d);

  // Add the contributions to the director field
  director::template
    addDirectorResidual<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn1, d1d, res);
  director::template
    addDirectorResidual<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn2, d2d, res);

  // Add the contribution from the rotation constraint (defined by the
  // rotational parametrization) - if any
  director::template
    addRotationConstraint<vars_per_node, offset, num_nodes>(vars, res);
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
  addAdjResXptProduct( int elemIndex,
                       double time,
                       TacsScalar scale,
                       const TacsScalar psi[],
                       const TacsScalar Xpts[],
                       const TacsScalar vars[],
                       const TacsScalar dvars[],
                       const TacsScalar ddvars[],
                       TacsScalar fXptSens[] ){
  // // Compute the number of quadrature points
  // const int nquad = quadrature::getNumQuadraturePoints();

  // // Get the reference axis
  // const A2D::Vec3& axis = transform->getRefAxis();

  // // Compute the normal directions
  // TacsScalar fn1[3*basis::NUM_NODES], fn2[3*basis::NUM_NODES];
  // TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // // Compute the frame normal and directors at each node
  // TacsScalar d1[dsize], d2[dsize], d1d[dsize], d2d[dsize];
  // TacsScalar d1dot[dsize], d2dot[dsize];
  // TacsScalar d1ddot[dsize], d2ddot[dsize];
  // memset(d1d, 0, dsize*sizeof(TacsScalar));
  // memset(d2d, 0, dsize*sizeof(TacsScalar));
  // director::template
  //   computeDirectorRates<vars_per_node, offset,
  //                        basis::NUM_NODES>(vars, dvars, ddvars, fn1, d1, d1dot, d1ddot);
  // director::template
  //   computeDirectorRates<vars_per_node, offset,
  //                        basis::NUM_NODES>(vars, dvars, ddvars, fn2, d2, d2dot, d2ddot);

  // // Compute the tying strain values
  // TacsScalar ety[basis::NUM_TYING_POINTS], dety[basis::NUM_TYING_POINTS];
  // memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));
  // model::template
  //   computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars, d1, d2, ety);

  // // Loop over each quadrature point and add the residual contribution
  // for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
  //   // Get the quadrature weight
  //   double pt[3];
  //   double weight = quadrature::getQuadraturePoint(quad_index, pt);

  //   // The transformation to the local beam coordinates
  //   A2D::Mat3x3 T;

  //   // Parametric location
  //   A2D::Vec3 X0;

  //   // Tangent to the beam
  //   A2D::Vec3 X0xi;

  //   // Interpolated normal directions
  //   A2D::Vec3 n1, n2;

  //   // Derivatives of the interpolated normal directions
  //   A2D::Vec3 n1xi, n2xi;

  //   // The values of the director fields and their derivatives
  //   A2D::ADVec3 u0xi, d01, d02, d01xi, d02xi;

  //   // Compute X, X,xi and the interpolated normal
  //   basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
  //   basis::template interpFields<3, 3>(pt, Xpts, X0.x);
  //   basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
  //   basis::template interpFields<3, 3>(pt, fn1, n1.x);
  //   basis::template interpFields<3, 3>(pt, fn2, n2.x);
  //   basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
  //   basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);
  //   basis::template interpFields<3, 3>(pt, d1, d01.x);
  //   basis::template interpFields<3, 3>(pt, d2, d02.x);
  //   basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
  //   basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

  //   // Compute the transformation at the quadrature point
  //   transform->computeTransform(X0xi.x, T.A);

  //   // Compute the inverse
  //   A2D::Mat3x3 Xd, Xdinv;
  //   A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
  //   A2D::Mat3x3Inverse invXd(Xd, Xdinv);

  //   // Compute the determinant of the transform
  //   A2D::Scalar detXd;
  //   A2D::Mat3x3Det computedetXd(weight, Xd, detXd);

  //   // Compute XdinvT = Xdinv * T
  //   A2D::Mat3x3 XdinvT;
  //   A2D::Mat3x3MatMult multXinvT(Xdinv, T, XdinvT);

  //   // Assemble the matrix Xdz1 = [n1,xi | 0 | 0] and Xdz2 = [n2,xi | 0 | 0 ]
  //   A2D::Mat3x3 Xdz1, Xdz2;
  //   A2D::Mat3x3FromVec3 assembleXdz1(n1xi, Xdz1);
  //   A2D::Mat3x3FromVec3 assembleXdz2(n2xi, Xdz2);

  //   // Compute Xdinvz1T = - Xdinv * Xdz1 * Xdinv * T
  //   A2D::Mat3x3 Xdinvz1T, Xdz1XdinvT;
  //   A2D::Mat3x3MatMult multXdz1XdinvT(Xdz1, XdinvT, Xdz1XdinvT);
  //   A2D::Mat3x3MatMult multXdinvz1T(-1.0, Xdinv, Xdz1XdinvT, Xdinvz1T);

  //   // Compute Xdinvz2T = - Xdinv * Xdz2 * Xdinv * T
  //   A2D::Mat3x3 Xdinvz2T, Xdz2XdinvT;
  //   A2D::Mat3x3MatMult multXdz2XdinvT(Xdz2, XdinvT, Xdz2XdinvT);
  //   A2D::Mat3x3MatMult multXdinvz2T(-1.0, Xdinv, Xdz2XdinvT, Xdinvz2T);

  //   // Assemble u0d, u1d and u2d
  //   A2D::ADMat3x3 u0d, u1d, u2d;
  //   A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);
  //   A2D::ADMat3x3FromADVec3 assembleu1d(d01xi, u1d);
  //   A2D::ADMat3x3FromADVec3 assembleu2d(d02xi, u2d);

  //   // Compute u0x = T^{T} * u0d * XdinvT
  //   A2D::ADMat3x3 u0dXdinvT, u0x;
  //   A2D::ADMat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
  //   A2D::MatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

  //   // Compute u1x = T^{T} * (u1d * XdinvT + u0d * XdinvzT)
  //   A2D::ADMat3x3 u1dXdinvT, u1x;
  //   A2D::ADMat3x3MatMult multu1d(u1d, XdinvT, u1dXdinvT);
  //   A2D::ADMat3x3MatMultAdd multu1dadd(u0d, Xdinvz1T, u1dXdinvT);
  //   A2D::MatTrans3x3ADMatMult multu1x(T, u1dXdinvT, u1x);

  //   // Compute u2x = T^{T} * (u2d * XdinvT + u0d * XdinvzT)
  //   A2D::ADMat3x3 u2dXdinvT, u2x;
  //   A2D::ADMat3x3MatMult multu2d(u2d, XdinvT, u2dXdinvT);
  //   A2D::ADMat3x3MatMultAdd multu2dadd(u0d, Xdinvz2T, u2dXdinvT);
  //   A2D::MatTrans3x3ADMatMult multu2x(T, u2dXdinvT, u2x);

  //   // Evaluate the tying components of the strain
  //   TacsScalar gty[2]; // The components of the tying strain
  //   basis::interpTyingStrain(pt, ety, gty);

  //   // Transform the tying strain to the local coordinates
  //   TacsScalar de0ty[2], e0ty[2] = {0.0, 0.0};
  //   // mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

  //   // Evaluate the strain
  //   TacsScalar e[6];
  //   model::evalStrain(u0x.A, u1x.A, u2x.A, e0ty, e);

  //   // Compute the corresponding stresses
  //   TacsScalar s[6];
  //   con->evalStress(elemIndex, pt, X0.x, e, s);

  //   // Evaluate the strain and strain derivatives from the
  //   model::evalStrainSens(detXd.value, s, u0x.A, u1x.A, u2x.A, e0ty,
  //                         u0x.Ad, u1x.Ad, u2x.Ad, de0ty);

  //   // Reverse the operations for the derivative w.r.t. state variables
  //   multu2x.reverse();
  //   multu2dadd.reverse();
  //   multu2d.reverse();
  //   multu1x.reverse();
  //   multu1dadd.reverse();
  //   multu1d.reverse();
  //   multu0x.reverse();
  //   multu0d.reverse();
  //   assembleu2d.reverse();
  //   assembleu1d.reverse();
  //   assembleu0d.reverse();

  //   // Add the residual contributions back to the element
  //   basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.xd, res);

  //   // Add the constributions back to the derivative
  //   basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, d1d);
  //   basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, d2d);
  //   basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, d1d);
  //   basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, d2d);

  //   // // Add the contributions to the tying strain
  //   // TacsScalar dgty[6];
  //   // mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

  //   // // Evaluate the tying strain
  //   // basis::addInterpTyingStrainTranspose(pt, dgty, dety);
  // }

  // // // Add the contributions from the tying strain
  // // model::template
  // //   addComputeTyingStrainTranspose<vars_per_node, basis>(Xpts, fn, vars,
  // //                                                        d, dety, res, dd);

  // // Add the contributions to the director field
  // director::template
  //   addDirectorResidual<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn1, d1d, res);
  // director::template
  //   addDirectorResidual<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn2, d2d, res);

  // // Add the contribution from the rotation constraint (defined by the
  // // rotational parametrization) - if any
  // director::template
  //   addRotationConstraint<vars_per_node, offset, num_nodes>(vars, res);
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
  getOutputData( int elemIndex,
                 ElementType etype,
                 int write_flag,
                 const TacsScalar Xpts[],
                 const TacsScalar vars[],
                 const TacsScalar dvars[],
                 const TacsScalar ddvars[],
                 int ld_data,
                 TacsScalar *data ){

  // // Get the number of nodes associated with the visualization
  // int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

  // // Compute the number of quadrature points
  // const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // // Compute the node normal directions
  // TacsScalar fn[3*basis::NUM_NODES];
  // getNodeNormals(Xpts, fn);

  // // Compute the frame normal and directors at each node
  // TacsScalar C[9*basis::NUM_NODES];
  // TacsScalar d[3*basis::NUM_NODES];
  // TacsScalar ddot[3*basis::NUM_NODES];
  // for ( int i = 0, offset = 4; i < basis::NUM_NODES; i++, offset += vars_per_node ){
  //   director::computeDirectorRates(&vars[offset], &dvars[offset], &fn[3*i],
  //                                  &C[9*i], &d[3*i], &ddot[3*i]);
  // }

  // // Set the total number of tying points needed for this element
  // TacsScalar ety[basis::NUM_TYING_POINTS];
  // model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // // Loop over each quadrature point and add the residual contribution
  // for ( int index = 0; index < num_vis_nodes; index++ ){
  //   // Get the quadrature weight
  //   double pt[3];
  //   basis::getNodePoint(index, pt);

  //   // Evaluate the displacement gradient at the point
  //   TacsScalar X[3], T[9];
  //   TacsScalar XdinvT[9], negXdinvXdz[9];
  //   TacsScalar u0x[9], u1x[9], Ct[9];
  //   computeDispGrad(pt, Xpts, vars, fn, C, d,
  //                   X, T, XdinvT, negXdinvXdz,
  //                   u0x, u1x, Ct);

  //   // Evaluate the tying components of the strain
  //   TacsScalar gty[6]; // The symmetric components of the tying strain
  //   model::template interpTyingStrain<basis>(pt, ety, gty);

  //   // Compute the symmetric parts of the tying strain
  //   TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
  //   mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

  //   // Compute the set of strain components
  //   TacsScalar e[9]; // The components of the strain
  //   model::evalStrain(u0x, u1x, e0ty, Ct, e);

  //   // Evaluate the temperature and temperature gradient
  //   TacsScalar t;
  //   basis::interpFields(pt, vars_per_node, &vars[3], 1, &t);

  //   // Compute the thermal strain
  //   TacsScalar et[9];
  //   con->evalThermalStrain(elemIndex, pt, X, t, et);

  //   // Compute the mechanical strain (and stress)
  //   TacsScalar em[9];
  //   for ( int i = 0; i < 9; i++ ){
  //     em[i] = e[i] - et[i];
  //   }

  //   // Compute the corresponding stresses
  //   TacsScalar s[9];
  //   con->evalStress(elemIndex, pt, X, em, s);

  //   if (etype == TACS_BEAM_OR_SHELL_ELEMENT){
  //     if (write_flag & TACS_OUTPUT_NODES){
  //       data[0] = X[0];
  //       data[1] = X[1];
  //       data[2] = X[2];
  //       data += 3;
  //     }
  //     if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
  //       int len = vars_per_node;
  //       if (len > 6){
  //         len = 6;
  //       }
  //       for ( int i = 0; i < len; i++ ){
  //         data[i] = vars[i + vars_per_node*index];
  //       }
  //       for ( int i = len; i < 6; i++ ){
  //         data[i] = 0.0;
  //       }
  //       data += 6;
  //     }
  //     if (write_flag & TACS_OUTPUT_STRAINS){
  //       for ( int i = 0; i < 9; i++ ){
  //         data[i] = e[i];
  //       }
  //       data += 9;
  //     }
  //     if (write_flag & TACS_OUTPUT_STRESSES){
  //       for ( int i = 0; i < 9; i++ ){
  //         data[i] = s[i];
  //       }
  //       data += 9;
  //     }
  //     if (write_flag & TACS_OUTPUT_EXTRAS){
  //       data[0] = con->evalFailure(elemIndex, pt, X, e);
  //       data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
  //       data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
  //       data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
  //       data += 4;
  //     }
  //   }
  // }
}

#endif // TACS_SHELL_ELEMENT_H
