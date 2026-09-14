#ifndef TACS_BEAM_ELEMENT_H
#define TACS_BEAM_ELEMENT_H

#include "TACSA2DUtilities.h"
#include "TACSBeamCentrifugalForce.h"
#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementModel.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamElementTransform.h"
#include "TACSBeamInertialForce.h"
#include "TACSBeamTraction.h"
#include "TACSBeamUtilities.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"

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


  Note that the second coefficient for d1x can be simplified to

  e1^{T} * x,d0^{-1},z1 * T^{T} * e1
  .   = - e1^{T} * x,d0^{-1} * x,d,z1 * x,d0^{-1} * T^{T} * e1
  .   = - (e1^{T} * x,d0^{-1} * n1,xi) * (e1^{T} * x,d0^{-1} * T^{T} * e1)

  The term d2x is

  d2x = T^{T} * (d2,xi * e1^{T} * x,d0^{-1} +
  .              u0,xi * e1^{T} * x,d0^{-1},z2 ) * T^{T} * e1

  These expressions can be simplified to the following:

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
  // Offset within the solution vector to the rotational
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

  TACSBeamElement(TACSBeamTransform *_transform, TACSBeamConstitutive *_con) {
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  ~TACSBeamElement() {
    if (transform) {
      transform->decref();
    }

    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSBeamElement"; }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return basis::NUM_NODES; }

  ElementLayout getLayoutType() { return basis::getLayoutType(); }

  ElementType getElementType() { return TACS_BEAM_OR_SHELL_ELEMENT; }

  int getNumQuadraturePoints() { return quadrature::getNumQuadraturePoints(); }

  double getQuadratureWeight(int n) {
    return quadrature::getQuadratureWeight(n);
  }

  double getQuadraturePoint(int n, double pt[]) {
    return quadrature::getQuadraturePoint(n, pt);
  }

  int getNumElementFaces() { return quadrature::getNumElementFaces(); }

  int getNumFaceQuadraturePoints(int face) {
    return quadrature::getNumFaceQuadraturePoints(face);
  }

  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return quadrature::getFaceQuadraturePoint(face, n, pt, tangent);
  }

  int getDesignVarNums(int elemIndex, int dvLen, int dvNums[]) {
    return con->getDesignVarNums(elemIndex, dvLen, dvNums);
  }

  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    return con->setDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    return con->getDesignVars(elemIndex, dvLen, dvs);
  }

  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]) {
    return con->getDesignVarRange(elemIndex, dvLen, lb, ub);
  }

  TACSElement *createElementTraction(int faceIndex, const TacsScalar t[]) {
    return new TACSBeamTraction<vars_per_node, quadrature, basis>(t);
  }

  TACSElement *createElementInertialForce(const TacsScalar inertiaVec[]) {
    return new TACSBeamInertialForce<vars_per_node, quadrature, basis>(
        transform, con, inertiaVec);
  }

  TACSElement *createElementCentrifugalForce(const TacsScalar omegaVec[],
                                             const TacsScalar rotCenter[],
                                             const bool first_order = false) {
    return new TACSBeamCentrifugalForce<vars_per_node, quadrature, basis>(
        transform, con, omegaVec, rotCenter);
  }

  void computeEnergies(int elemIndex, double time, const TacsScalar Xpts[],
                       const TacsScalar vars[], const TacsScalar dvars[],
                       TacsScalar *Te, TacsScalar *Pe);

  void addResidual(int elemIndex, double time, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res);

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

  void getMatType(ElementMatrixType matType, int elemIndex, double time,
                  const TacsScalar Xpts[], const TacsScalar vars[],
                  TacsScalar mat[]);

  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dfdx[]);

  void addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                           const TacsScalar psi[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], TacsScalar fXptSens[]);

  int evalPointQuantity(int elemIndex, int quantityType, double time, int n,
                        double pt[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar *detXd,
                        TacsScalar *quantity);

  void addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                              TacsScalar scale, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], int dvLen,
                              TacsScalar dfdx[]);

  void addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                              TacsScalar alpha, TacsScalar beta,
                              TacsScalar gamma, int n, double pt[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              const TacsScalar dvars[],
                              const TacsScalar ddvars[],
                              const TacsScalar dfdq[], TacsScalar dfdu[]);

  void addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfddetXd,
                               const TacsScalar dfdq[], TacsScalar dfdXpts[]);

  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

 private:
  // Set sizes for the different components
  static const int usize = 3 * basis::NUM_NODES;
  static const int dsize = 3 * basis::NUM_NODES;
  static const int csize = 9 * basis::NUM_NODES;

  TACSBeamTransform *transform;
  TACSBeamConstitutive *con;
};

/*
  Compute the kinetic and potential energies of the shell
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::computeEnergies(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, TacsScalar *Telem, TacsScalar *Uelem) {
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn1,
                                                            d1, d1dot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn2,
                                                            d2, d2dot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::Mat<TacsScalar, 3, 3> T;

    // Parametric location
    A2D::Vec<TacsScalar, 3> X0;

    // Tangent to the beam
    A2D::Vec<TacsScalar, 3> X0xi;

    // Interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::Vec<TacsScalar, 3> u0xi, d01, d02, d01xi, d02xi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                       A2D::get_data(u0xi));
    basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
    basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
    basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
    basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
    basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
    basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
    basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
    basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

    // Compute the transformation at the quadrature point
    transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

    // Compute the inverse
    A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
    A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
    A2D::MatInv(Xd, Xdinv);

    // Compute the determinant of the transform
    TacsScalar detXd;
    A2D::MatDet(Xd, detXd);
    detXd *= weight;

    // Compute XdinvT = Xdinv * T
    A2D::Mat<TacsScalar, 3, 3> XdinvT;
    A2D::MatMatMult(Xdinv, T, XdinvT);

    // Assemble u0d
    A2D::Mat<TacsScalar, 3, 3> u0d;
    A2D::MatFromThreeVec(u0xi, d01, d02, u0d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::Mat<TacsScalar, 3, 3> u0dXdinvT, u0x;
    A2D::MatMatMult(u0d, XdinvT, u0dXdinvT);
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                               u0x);

    // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
    // and sz2 = e1^{T} * Xdinv * n2xi
    const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
    A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
    TacsScalar s0, sz1, sz2;
    A2D::MatVecMult(XdinvT, e1, ts0);
    A2D::VecDot(e1, ts0, s0);
    A2D::MatVecMult(Xdinv, n1xi, tsz1);
    A2D::VecDot(e1, tsz1, sz1);
    A2D::MatVecMult(Xdinv, n2xi, tsz2);
    A2D::VecDot(e1, tsz2, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::Vec<TacsScalar, 3> d1t, d1s, d1x;
    A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s);
    A2D::VecScale(s0, d1s, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::Vec<TacsScalar, 3> d2t, d2s, d2x;
    A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s);
    A2D::VecScale(s0, d2s, d2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2];
    e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
    e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

    // Compute the set of strain components
    TacsScalar e[6];  // The components of the strain
    model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x),
                      A2D::get_data(d2x), e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);

    Ue += 0.5 * detXd *
          (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] + s[4] * e[4] +
           s[5] * e[5]);

    // Evaluate the velocities
    A2D::Vec<TacsScalar, 3> u0dot, d01dot, d02dot;
    basis::template interpFields<vars_per_node, 3>(pt, dvars,
                                                   A2D::get_data(u0dot));
    basis::template interpFields<3, 3>(pt, d1dot, A2D::get_data(d01dot));
    basis::template interpFields<3, 3>(pt, d2dot, A2D::get_data(d02dot));

    // dot{u} = \dot{u0} + z1 * dot{d1} + z2 * dot{d2}
    TacsScalar u0d0, u0d10, u0d20, d1d10, d2d20, d1d20;
    A2D::VecDot(u0dot, u0dot, u0d0);
    A2D::VecDot(u0dot, d01dot, u0d10);
    A2D::VecDot(u0dot, d02dot, u0d20);
    A2D::VecDot(d01dot, d01dot, d1d10);
    A2D::VecDot(d02dot, d02dot, d2d20);
    A2D::VecDot(d01dot, d02dot, d1d20);

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), rho);

    Te += 0.5 * detXd *
          (rho[0] * u0d0 + 2.0 * rho[1] * u0d10 + 2.0 * rho[2] * u0d20 +
           rho[3] * d1d10 + rho[4] * d2d20 + 2.0 * rho[5] * d1d20);
  }

  *Telem = Te;
  *Uelem = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(
      vars, dvars, ddvars, fn1, d1, d1dot, d1ddot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(
      vars, dvars, ddvars, fn2, d2, d2dot, d2ddot);

  // Add the contributions to the derivative
  TacsScalar d1d[dsize], d2d[dsize];
  memset(d1d, 0, dsize * sizeof(TacsScalar));
  memset(d2d, 0, dsize * sizeof(TacsScalar));

  // Compute the tying strain values
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::Mat<TacsScalar, 3, 3> T;

    // Parametric location
    A2D::Vec<TacsScalar, 3> X0;

    // Tangent to the beam
    A2D::Vec<TacsScalar, 3> X0xi;

    // Interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0xi, d01, d02, d01xi, d02xi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                       A2D::get_data(u0xi));
    basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
    basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
    basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
    basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
    basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
    basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
    basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
    basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

    // Compute the transformation at the quadrature point
    transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

    // Compute the inverse
    A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
    A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
    A2D::MatInv(Xd, Xdinv);

    // Compute the determinant of the transform
    TacsScalar detXd;
    A2D::MatDet(Xd, detXd);
    detXd *= weight;

    // Compute XdinvT = Xdinv * T
    A2D::Mat<TacsScalar, 3, 3> XdinvT;
    A2D::MatMatMult(Xdinv, T, XdinvT);

    // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
    // and sz2 = e1^{T} * Xdinv * n2xi
    const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
    A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
    TacsScalar s0, sz1, sz2;
    A2D::MatVecMult(XdinvT, e1, ts0);
    A2D::VecDot(e1, ts0, s0);
    A2D::MatVecMult(Xdinv, n1xi, tsz1);
    A2D::VecDot(e1, tsz1, sz1);
    A2D::MatVecMult(Xdinv, n2xi, tsz2);
    A2D::VecDot(e1, tsz2, sz2);

    // Assemble u0d and compute u0x = T^{T} * u0d * XdinvT,
    // d1x = s0 * T^{T} * (d1xi - sz1 * u0xi) and
    // d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> u0d, u0dXdinvT, u0x;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> d1t, d1s, d1x, d2t, d2s, d2x;
    auto strain_stack = A2D::MakeStack(
        A2D::MatFromThreeVec(u0xi, d01, d02, u0d),
        A2D::MatMatMult(u0d, XdinvT, u0dXdinvT),
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                                   u0x),
        A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s),
        A2D::VecScale(s0, d1s, d1x),
        A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s),
        A2D::VecScale(s0, d2s, d2x));

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], de0ty[2];
    e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
    e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

    // Evaluate the strain
    TacsScalar e[6];
    model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x),
                      A2D::get_data(d2x), e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);

    // Evaluate the strain and strain derivatives from the
    model::evalStrainSens(detXd, s, A2D::get_data(u0x), A2D::get_data(d1x),
                          A2D::get_data(d2x), e0ty, A2D::get_data(u0x.bvalue()),
                          A2D::get_data(d1x.bvalue()),
                          A2D::get_data(d2x.bvalue()), de0ty);

    // Convert the contributions to the tying strain
    TacsScalar dgty[2];
    dgty[0] = 2.0 * A2D::get_data(XdinvT)[0] * de0ty[0];
    dgty[1] = 2.0 * A2D::get_data(XdinvT)[0] * de0ty[1];

    strain_stack.reverse();

    // Add the residual contributions back to the element
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(
        pt, A2D::get_data(u0xi.bvalue()), res);

    // Add the constributions back to the derivative
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01.bvalue()), d1d);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02.bvalue()), d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d01xi.bvalue()), d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d02xi.bvalue()), d2d);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    // Evaluate the accelerations
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0ddot, d01ddot, d02ddot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars,
                                                   A2D::get_data(u0ddot));
    basis::template interpFields<3, 3>(pt, d1ddot, A2D::get_data(d01ddot));
    basis::template interpFields<3, 3>(pt, d2ddot, A2D::get_data(d02ddot));

    // dot{u}(xi, z1, z2) = dot{u0} + z1 * dot{d1} + z2 * dot{d2}
    A2D::ADObj<TacsScalar> u0d0, u0d10, u0d20, d1d10, d2d20, d1d20;
    auto accel_stack = A2D::MakeStack(A2D::VecDot(u0ddot, u0ddot, u0d0),
                                      A2D::VecDot(u0ddot, d01ddot, u0d10),
                                      A2D::VecDot(u0ddot, d02ddot, u0d20),
                                      A2D::VecDot(d01ddot, d01ddot, d1d10),
                                      A2D::VecDot(d02ddot, d02ddot, d2d20),
                                      A2D::VecDot(d01ddot, d02ddot, d1d20));

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), rho);

    u0d0.bvalue() = 0.5 * rho[0] * detXd;
    u0d10.bvalue() = rho[1] * detXd;
    u0d20.bvalue() = rho[2] * detXd;
    d1d10.bvalue() = 0.5 * rho[3] * detXd;
    d2d20.bvalue() = 0.5 * rho[4] * detXd;
    d1d20.bvalue() = rho[5] * detXd;

    accel_stack.reverse();

    basis::template addInterpFieldsTranspose<vars_per_node, 3>(
        pt, A2D::get_data(u0ddot.bvalue()), res);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01ddot.bvalue()), d1d);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02ddot.bvalue()), d2d);
  }

  // Add the contributions from the tying strain
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn1, fn2, vars, d1, d2, dety, res, d1d, d2d);

  // Add the contributions to the director field
  director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn1, d1d, res);
  director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn2, d2d, res);

  // Add the contribution from the rotation constraint (defined by the
  // rotational parametrization) - if any
  director::template addRotationConstraint<vars_per_node, offset, num_nodes>(
      vars, res);
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::getMatType(
    ElementMatrixType matType, int elemIndex, double time,
    const TacsScalar Xpts[], const TacsScalar vars[], TacsScalar mat[]) {
  memset(mat, 0,
         vars_per_node * num_nodes * vars_per_node * num_nodes *
             sizeof(TacsScalar));
  TacsScalar alpha, beta, gamma;
  alpha = beta = gamma = 0.0;
  // Set alpha or gamma based on if this is a stiffness or mass matrix
  if (matType == TACS_STIFFNESS_MATRIX) {
    alpha = 1.0;
  } else if (matType == TACS_MASS_MATRIX) {
    gamma = 1.0;
  } else {  // TACS_GEOMETRIC_STIFFNESS_MATRIX
    // Not implemented
    return;
  }
  // Create dummy residual vector
  TacsScalar res[vars_per_node * num_nodes];
  memset(res, 0, vars_per_node * num_nodes * sizeof(TacsScalar));
  // Add appropriate Jacobian to matrix
  addJacobian(elemIndex, time, alpha, beta, gamma, Xpts, vars, vars, vars, res,
              mat);
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize], d1psi[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize], d2psi[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               num_nodes>(
      vars, dvars, ddvars, psi, fn1, d1, d1dot, d1ddot, d1psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               num_nodes>(
      vars, dvars, ddvars, psi, fn2, d2, d2dot, d2ddot, d2psi);

  // Compute the tying strain values
  TacsScalar ety[basis::NUM_TYING_POINTS], etypsi[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn1, fn2, vars, d1, d2, psi, d1psi, d2psi, ety, etypsi);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::Mat<TacsScalar, 3, 3> T;

    // Parametric location
    A2D::Vec<TacsScalar, 3> X0;

    // Tangent to the beam
    A2D::Vec<TacsScalar, 3> X0xi;

    // Interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::Vec<TacsScalar, 3> u0xi, d01, d02, d01xi, d02xi;
    A2D::Vec<TacsScalar, 3> u0xipsi, d01psi, d02psi, d01xipsi, d02xipsi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                       A2D::get_data(u0xi));
    basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
    basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
    basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
    basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

    // Interpolate the adjoint solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi,
                                                       A2D::get_data(u0xipsi));
    basis::template interpFields<3, 3>(pt, d1psi, A2D::get_data(d01psi));
    basis::template interpFields<3, 3>(pt, d2psi, A2D::get_data(d02psi));
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, A2D::get_data(d01xipsi));
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, A2D::get_data(d02xipsi));

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
    basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
    basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
    basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
    basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

    // Compute the transformation at the quadrature point
    transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

    // Compute the inverse
    A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
    A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
    A2D::MatInv(Xd, Xdinv);

    // Compute the determinant of the transform
    TacsScalar detXd;
    A2D::MatDet(Xd, detXd);
    detXd *= weight;

    // Compute XdinvT = Xdinv * T
    A2D::Mat<TacsScalar, 3, 3> XdinvT;
    A2D::MatMatMult(Xdinv, T, XdinvT);

    // Assemble u0d and u0psi
    A2D::Mat<TacsScalar, 3, 3> u0d, u0dpsi;
    A2D::MatFromThreeVec(u0xi, d01, d02, u0d);
    A2D::MatFromThreeVec(u0xipsi, d01psi, d02psi, u0dpsi);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::Mat<TacsScalar, 3, 3> u0dXdinvT, u0x;
    A2D::MatMatMult(u0d, XdinvT, u0dXdinvT);
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                               u0x);

    // Compute u0xpsi = ^{T} * u0dpsi * XdinvT
    A2D::Mat<TacsScalar, 3, 3> u0dXdinvTpsi, u0xpsi;
    A2D::MatMatMult(u0dpsi, XdinvT, u0dXdinvTpsi);
    A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvTpsi,
                                                               u0xpsi);

    // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
    // and sz2 = e1^{T} * Xdinv * n2xi
    const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
    A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
    TacsScalar s0, sz1, sz2;
    A2D::MatVecMult(XdinvT, e1, ts0);
    A2D::VecDot(e1, ts0, s0);
    A2D::MatVecMult(Xdinv, n1xi, tsz1);
    A2D::VecDot(e1, tsz1, sz1);
    A2D::MatVecMult(Xdinv, n2xi, tsz2);
    A2D::VecDot(e1, tsz2, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::Vec<TacsScalar, 3> d1t, d1s, d1x;
    A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s);
    A2D::VecScale(s0, d1s, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::Vec<TacsScalar, 3> d2t, d2s, d2x;
    A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s);
    A2D::VecScale(s0, d2s, d2x);

    // Compute d1xpsi = s0 * T^{T} * (d1xipsi - sz1 * u0xipsi)
    A2D::Vec<TacsScalar, 3> d1tpsi, d1spsi, d1xpsi;
    A2D::VecSum(-sz1, u0xipsi, TacsScalar(1.0), d01xipsi, d1tpsi);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1tpsi, d1spsi);
    A2D::VecScale(s0, d1spsi, d1xpsi);

    // Compute d2xpsi = s0 * T^{T} * (d2xipsi - sz2 * u0xipsi)
    A2D::Vec<TacsScalar, 3> d2tpsi, d2spsi, d2xpsi;
    A2D::VecSum(-sz2, u0xipsi, TacsScalar(1.0), d02xipsi, d2tpsi);
    A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2tpsi, d2spsi);
    A2D::VecScale(s0, d2spsi, d2xpsi);

    // Evaluate the tying components of the strain
    TacsScalar gty[2], gtypsi[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);
    basis::interpTyingStrain(pt, etypsi, gtypsi);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], e0typsi[2];
    e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
    e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];
    e0typsi[0] = 2.0 * A2D::get_data(XdinvT)[0] * gtypsi[0];
    e0typsi[1] = 2.0 * A2D::get_data(XdinvT)[0] * gtypsi[1];

    // // Evaluate the strain and the adjoint strain
    TacsScalar e[6], epsi[6];
    model::evalStrainDeriv(A2D::get_data(u0x), A2D::get_data(d1x),
                           A2D::get_data(d2x), e0ty, A2D::get_data(u0xpsi),
                           A2D::get_data(d1xpsi), A2D::get_data(d2xpsi),
                           e0typsi, e, epsi);

    // Add the product of the derivative of the stress
    con->addStressDVSens(elemIndex, scale * detXd, pt, A2D::get_data(X0), e,
                         epsi, dvLen, dfdx);

    // Evaluate the accelerations
    A2D::Vec<TacsScalar, 3> u0dot, d01dot, d02dot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars,
                                                   A2D::get_data(u0dot));
    basis::template interpFields<3, 3>(pt, d1ddot, A2D::get_data(d01dot));
    basis::template interpFields<3, 3>(pt, d2ddot, A2D::get_data(d02dot));

    A2D::Vec<TacsScalar, 3> u0dotpsi, d01dotpsi, d02dotpsi;
    basis::template interpFields<vars_per_node, 3>(pt, psi,
                                                   A2D::get_data(u0dotpsi));
    basis::template interpFields<3, 3>(pt, d1psi, A2D::get_data(d01dotpsi));
    basis::template interpFields<3, 3>(pt, d2psi, A2D::get_data(d02dotpsi));

    // Compute the dot-products
    TacsScalar u0psi, u0psid1, u0psid2, u0d1psi, u0d2psi;
    A2D::VecDot(u0dot, u0dotpsi, u0psi);
    A2D::VecDot(u0dotpsi, d01dot, u0psid1);
    A2D::VecDot(u0dot, d01dotpsi, u0d1psi);
    A2D::VecDot(u0dotpsi, d02dot, u0psid2);
    A2D::VecDot(u0dot, d02dotpsi, u0d2psi);

    TacsScalar d1d1psi, d2d2psi, d1psid2, d1d2psi;
    A2D::VecDot(d01dot, d01dotpsi, d1d1psi);
    A2D::VecDot(d02dot, d02dotpsi, d2d2psi);
    A2D::VecDot(d01dotpsi, d02dot, d1psid2);
    A2D::VecDot(d01dot, d02dotpsi, d1d2psi);

    // Add derivatives from the mass moments
    TacsScalar rho[6];
    TacsScalar alpha = scale * detXd;
    rho[0] = alpha * u0psi;
    rho[1] = alpha * (u0psid1 + u0d1psi);
    rho[2] = alpha * (u0psid2 + u0d2psi);
    rho[3] = alpha * d1d1psi;
    rho[4] = alpha * d2d2psi;
    rho[5] = alpha * (d1psid2 + d1d2psi);

    con->addMassMomentsDVSens(elemIndex, pt, A2D::get_data(X0), rho, dvLen,
                              dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::addAdjResXptProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar dfdXpts[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Derivatives w.r.t. the frame normals
  TacsScalar dfn1[3 * basis::NUM_NODES], dfn2[3 * basis::NUM_NODES];
  memset(dfn1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
  memset(dfn2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize], d1psi[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize], d2psi[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               num_nodes>(
      vars, dvars, ddvars, psi, fn1, d1, d1dot, d1ddot, d1psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               num_nodes>(
      vars, dvars, ddvars, psi, fn2, d2, d2dot, d2ddot, d2psi);

  // Derivatives of the adjoint-residual product w.r.t. d1/d2 and d1psi/d2psi
  TacsScalar dd1[dsize], dd1psi[dsize];
  TacsScalar dd2[dsize], dd2psi[dsize];
  memset(dd1, 0, dsize * sizeof(TacsScalar));
  memset(dd2, 0, dsize * sizeof(TacsScalar));
  memset(dd1psi, 0, dsize * sizeof(TacsScalar));
  memset(dd2psi, 0, dsize * sizeof(TacsScalar));

  // Compute the tying strain values
  TacsScalar ety[basis::NUM_TYING_POINTS], etypsi[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn1, fn2, vars, d1, d2, psi, d1psi, d2psi, ety, etypsi);

  TacsScalar etyd[basis::NUM_TYING_POINTS], etypsid[basis::NUM_TYING_POINTS];
  memset(etyd, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(etypsid, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // The transformation to the local beam coordinates
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> T;

    // Parametric location
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0;

    // Tangent to the beam
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0xi;

    // Interpolated normal directions
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0xi, d01, d02, d01xi, d02xi;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0xipsi, d01psi, d02psi, d01xipsi,
        d02xipsi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                       A2D::get_data(u0xi));
    basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
    basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
    basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
    basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

    // Interpolate the adjoint solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi,
                                                       A2D::get_data(u0xipsi));
    basis::template interpFields<3, 3>(pt, d1psi, A2D::get_data(d01psi));
    basis::template interpFields<3, 3>(pt, d2psi, A2D::get_data(d02psi));
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, A2D::get_data(d01xipsi));
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, A2D::get_data(d02xipsi));

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
    basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
    basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
    basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
    basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

    // Compute the transformation at the quadrature point
    transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

    // Assemble Xd, compute its inverse and (weighted) determinant, compute
    // XdinvT = Xdinv * T, u0x = T^{T} * u0d * XdinvT (and the psi
    // counterpart), the scale factors s0 = e1^{T} * XdinvT * e1,
    // sz1 = e1^{T} * Xdinv * n1xi and sz2 = e1^{T} * Xdinv * n2xi, and
    // dkx = s0 * T^{T} * (dkxi - szk * u0xi) for k = 1, 2 (and the psi
    // counterparts)
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> Xd, Xdinv, XdinvT;
    A2D::ADObj<TacsScalar> detXd0, detXd;
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> u0d, u0dpsi, u0dXdinvT, u0x;
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> u0dXdinvTpsi, u0xpsi;
    const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
    A2D::Vec<TacsScalar, 3> e1(e1_data);
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> ts0, tsz1, tsz2;
    A2D::ADObj<TacsScalar> s0, sz1, sz2;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> w1, d1t, d1s, d1x;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> w2, d2t, d2s, d2x;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> w1psi, d1tpsi, d1spsi, d1xpsi;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> w2psi, d2tpsi, d2spsi, d2xpsi;

    auto strain_stack = A2D::MakeStack(
        A2D::MatFromThreeVec(X0xi, n1, n2, Xd), A2D::MatInv(Xd, Xdinv),
        A2D::MatDet(Xd, detXd0), A2D::Eval(detXd0 * weight, detXd),
        A2D::MatMatMult(Xdinv, T, XdinvT),
        A2D::MatFromThreeVec(u0xi, d01, d02, u0d),
        A2D::MatFromThreeVec(u0xipsi, d01psi, d02psi, u0dpsi),
        A2D::MatMatMult(u0d, XdinvT, u0dXdinvT),
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                                   u0x),
        A2D::MatMatMult(u0dpsi, XdinvT, u0dXdinvTpsi),
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(
            T, u0dXdinvTpsi, u0xpsi),
        A2D::MatVecMult(XdinvT, e1, ts0), A2D::VecDot(e1, ts0, s0),
        A2D::MatVecMult(Xdinv, n1xi, tsz1), A2D::VecDot(e1, tsz1, sz1),
        A2D::MatVecMult(Xdinv, n2xi, tsz2), A2D::VecDot(e1, tsz2, sz2),
        A2D::VecScale(sz1, u0xi, w1),
        A2D::VecSum(TacsScalar(1.0), d01xi, TacsScalar(-1.0), w1, d1t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s),
        A2D::VecScale(s0, d1s, d1x), A2D::VecScale(sz2, u0xi, w2),
        A2D::VecSum(TacsScalar(1.0), d02xi, TacsScalar(-1.0), w2, d2t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s),
        A2D::VecScale(s0, d2s, d2x), A2D::VecScale(sz1, u0xipsi, w1psi),
        A2D::VecSum(TacsScalar(1.0), d01xipsi, TacsScalar(-1.0), w1psi, d1tpsi),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1tpsi, d1spsi),
        A2D::VecScale(s0, d1spsi, d1xpsi), A2D::VecScale(sz2, u0xipsi, w2psi),
        A2D::VecSum(TacsScalar(1.0), d02xipsi, TacsScalar(-1.0), w2psi, d2tpsi),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2tpsi, d2spsi),
        A2D::VecScale(s0, d2spsi, d2xpsi));

    // Evaluate the tying components of the strain
    TacsScalar gty[2], gtypsi[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);
    basis::interpTyingStrain(pt, etypsi, gtypsi);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], e0typsi[2];
    e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
    e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];
    e0typsi[0] = 2.0 * A2D::get_data(XdinvT)[0] * gtypsi[0];
    e0typsi[1] = 2.0 * A2D::get_data(XdinvT)[0] * gtypsi[1];

    // // Evaluate the strain and the adjoint strain
    TacsScalar e[6], epsi[6];
    model::evalStrainDeriv(A2D::get_data(u0x), A2D::get_data(d1x),
                           A2D::get_data(d2x), e0ty, A2D::get_data(u0xpsi),
                           A2D::get_data(d1xpsi), A2D::get_data(d2xpsi),
                           e0typsi, e, epsi);

    // Compute the stress due to the strain
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);

    // Compute the psioint stress - assuming a linear relationship for
    // stress/strain
    TacsScalar spsi[6];
    con->evalStress(elemIndex, pt, A2D::get_data(X0), epsi, spsi);

    // Evaluate the sensitivities
    TacsScalar e0tyd[2], e0typsid[2];
    model::evalStrainSens(
        scale * detXd.value(), spsi, A2D::get_data(u0x), A2D::get_data(d1x),
        A2D::get_data(d2x), e0ty, A2D::get_data(u0x.bvalue()),
        A2D::get_data(d1x.bvalue()), A2D::get_data(d2x.bvalue()), e0tyd);
    model::evalStrainSens(scale * detXd.value(), s, A2D::get_data(u0x),
                          A2D::get_data(d1x), A2D::get_data(d2x), e0ty,
                          A2D::get_data(u0xpsi.bvalue()),
                          A2D::get_data(d1xpsi.bvalue()),
                          A2D::get_data(d2xpsi.bvalue()), e0typsid);
    detXd.bvalue() = scale * (e[0] * spsi[0] + e[1] * spsi[1] + e[2] * spsi[2] +
                              e[3] * spsi[3] + e[4] * spsi[4] + e[5] * spsi[5]);

    // Apply the tying strain transformation
    TacsScalar gtyd[2], gtypsid[2];
    gtyd[0] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[0];
    gtyd[1] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[1];
    gtypsid[0] = 2.0 * A2D::get_data(XdinvT)[0] * e0typsid[0];
    gtypsid[1] = 2.0 * A2D::get_data(XdinvT)[0] * e0typsid[1];

    A2D::get_data(XdinvT.bvalue())[0] +=
        2.0 * (gty[0] * e0tyd[0] + gty[1] * e0tyd[1] + e0typsid[0] * gtypsi[0] +
               e0typsid[1] * gtypsi[1]);

    // Evaluate the accelerations
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0ddot, d01ddot, d02ddot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars,
                                                   A2D::get_data(u0ddot));
    basis::template interpFields<3, 3>(pt, d1ddot, A2D::get_data(d01ddot));
    basis::template interpFields<3, 3>(pt, d2ddot, A2D::get_data(d02ddot));

    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0ddotpsi, d01ddotpsi, d02ddotpsi;
    basis::template interpFields<vars_per_node, 3>(pt, psi,
                                                   A2D::get_data(u0ddotpsi));
    basis::template interpFields<3, 3>(pt, d1psi, A2D::get_data(d01ddotpsi));
    basis::template interpFields<3, 3>(pt, d2psi, A2D::get_data(d02ddotpsi));

    // Compute the dot-products
    A2D::ADObj<TacsScalar> u0psi, u0psid1, u0psid2, u0d1psi, u0d2psi;
    A2D::ADObj<TacsScalar> d1d1psi, d2d2psi, d1psid2, d1d2psi;
    auto accel_stack =
        A2D::MakeStack(A2D::VecDot(u0ddot, u0ddotpsi, u0psi),
                       A2D::VecDot(u0ddotpsi, d01ddot, u0psid1),
                       A2D::VecDot(u0ddot, d01ddotpsi, u0d1psi),
                       A2D::VecDot(u0ddotpsi, d02ddot, u0psid2),
                       A2D::VecDot(u0ddot, d02ddotpsi, u0d2psi),
                       A2D::VecDot(d01ddot, d01ddotpsi, d1d1psi),
                       A2D::VecDot(d02ddot, d02ddotpsi, d2d2psi),
                       A2D::VecDot(d01ddotpsi, d02ddot, d1psid2),
                       A2D::VecDot(d01ddot, d02ddotpsi, d1d2psi));

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), rho);

    // Add the contribution from the adjoint-residual product from the
    // dynamics
    detXd.bvalue() +=
        scale *
        (rho[0] * u0psi.value() + rho[1] * (u0psid1.value() + u0d1psi.value()) +
         rho[2] * (u0psid2.value() + u0d2psi.value()) +
         rho[3] * d1d1psi.value() + rho[4] * d2d2psi.value() +
         rho[5] * (d1psid2.value() + d1d2psi.value()));

    // Set the seeds for the dot-products
    TacsScalar alpha = scale * detXd.value();
    u0psi.bvalue() = alpha * rho[0];
    u0psid1.bvalue() = alpha * rho[1];
    u0d1psi.bvalue() = alpha * rho[1];

    u0psid2.bvalue() = alpha * rho[2];
    u0d2psi.bvalue() = alpha * rho[2];
    d1d1psi.bvalue() = alpha * rho[3];
    d2d2psi.bvalue() = alpha * rho[4];
    d1psid2.bvalue() = alpha * rho[5];
    d1d2psi.bvalue() = alpha * rho[5];

    // Reverse the dot-products, then the strain expressions
    accel_stack.reverse();
    strain_stack.reverse();

    // Reverse the transformation sensitivities
    transform->addTransformSens(A2D::get_data(X0xi), A2D::get_data(T.bvalue()),
                                A2D::get_data(X0xi.bvalue()));

    // Add the sensitivities to the input fields...
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(X0.bvalue()), dfdXpts);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(X0xi.bvalue()), dfdXpts);

    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(n1.bvalue()), dfn1);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(n2.bvalue()), dfn2);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(n1xi.bvalue()), dfn1);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(n2xi.bvalue()), dfn2);

    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01.bvalue()), dd1);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02.bvalue()), dd2);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d01xi.bvalue()), dd1);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d02xi.bvalue()), dd2);

    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01psi.bvalue()), dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02psi.bvalue()), dd2psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d01xipsi.bvalue()), dd1psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d02xipsi.bvalue()), dd2psi);

    // Add the contributions from the dynamics
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01ddot.bvalue()), dd1);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02ddot.bvalue()), dd2);

    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01ddotpsi.bvalue()), dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02ddotpsi.bvalue()), dd2psi);

    // Add the contributions to the tying strain
    basis::addInterpTyingStrainTranspose(pt, gtyd, etyd);
    basis::addInterpTyingStrainTranspose(pt, gtypsid, etypsid);
  }

  // Add the sensitivity contributions from the tying strain
  model::template addTyingStrainDerivXptSens<vars_per_node, basis>(
      Xpts, fn1, fn2, vars, d1, d2, psi, d1psi, d2psi, etyd, etypsid, dfdXpts,
      dfn1, dfn2, dd1, dd2, dd1psi, dd2psi);

  // Add the contributions from the derivative of the director
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn1, dd1, dd1psi, dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn2, dd2, dd2psi, dfn2);

  // Add the contributions from the node normals
  TacsBeamAddNodeNormalsSens<basis>(Xpts, axis, dfn1, dfn2, dfdXpts);
}

template <class quadrature, class basis, class director, class model>
int TACSBeamElement<quadrature, basis, director, model>::evalPointQuantity(
    int elemIndex, int quantityType, double time, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *detXdval, TacsScalar *quantity) {
  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn1,
                                                            d1, d1dot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn2,
                                                            d2, d2dot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  // The transformation to the local beam coordinates
  A2D::Mat<TacsScalar, 3, 3> T;

  // Parametric location
  A2D::Vec<TacsScalar, 3> X0;

  // Tangent to the beam
  A2D::Vec<TacsScalar, 3> X0xi;

  // Interpolated normal directions
  A2D::Vec<TacsScalar, 3> n1, n2;

  // Derivatives of the interpolated normal directions
  A2D::Vec<TacsScalar, 3> n1xi, n2xi;

  // The values of the director fields and their derivatives
  A2D::Vec<TacsScalar, 3> u0xi, d01, d02, d01xi, d02xi;

  // Interpolate the solution fields
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                     A2D::get_data(u0xi));
  basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
  basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
  basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
  basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

  // Compute X, X,xi and the interpolated normal
  basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
  basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
  basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
  basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
  basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

  // Compute the transformation at the quadrature point
  transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

  // Compute the inverse
  A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
  A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
  A2D::MatInv(Xd, Xdinv);

  // Compute the determinant of the transform
  TacsScalar detXd;
  A2D::MatDet(Xd, detXd);

  if (detXdval) {
    *detXdval = detXd;
  }

  if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = con->evalDensity(elemIndex, pt, A2D::get_data(X0));
    }
    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      // Compute the interpolated displacements
      basis::template interpFields<vars_per_node, 3>(pt, vars, quantity);
    }
    return 3;
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    if (quantity) {
      TacsScalar mass_moment[6];
      con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), mass_moment);
      TacsScalar density = mass_moment[0];

      for (int i = 0; i < 3; i++) {
        quantity[i] = density * A2D::get_data(X0)[i] -
                      mass_moment[1] * A2D::get_data(n1)[i] -
                      mass_moment[2] * A2D::get_data(n2)[i];
      }
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      TacsScalar I0[6] = {0.0};

      // Evaluate the self MOI
      TacsScalar moments[6];
      con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), moments);
      TacsScalar density = moments[0];
      I0[3] = moments[4] - moments[2] * moments[2] / density;
      I0[4] = -moments[5] + moments[1] * moments[2] / density;
      I0[5] = moments[3] - moments[1] * moments[1] / density;
      // Compute T*I0*T^{T}
      mat3x3SymmTransform(A2D::get_data(T), I0, quantity);
      TacsScalar dXcg[3];
      for (int i = 0; i < 3; i++) {
        dXcg[i] = A2D::get_data(X0)[i] - (moments[1] * A2D::get_data(n1)[i] +
                                          moments[2] * A2D::get_data(n2)[i]) /
                                             density;
      }

      // Use parallel axis theorem to move MOI to origin
      quantity[0] += density * (dXcg[1] * dXcg[1] + dXcg[2] * dXcg[2]);
      quantity[1] += -density * dXcg[0] * dXcg[1];
      quantity[2] += -density * dXcg[0] * dXcg[2];
      quantity[3] += density * (dXcg[0] * dXcg[0] + dXcg[2] * dXcg[2]);
      quantity[4] += -density * dXcg[2] * dXcg[1];
      quantity[5] += density * (dXcg[0] * dXcg[0] + dXcg[1] * dXcg[1]);
    }

    return 6;
  }

  // Compute XdinvT = Xdinv * T
  A2D::Mat<TacsScalar, 3, 3> XdinvT;
  A2D::MatMatMult(Xdinv, T, XdinvT);

  // Assemble u0d
  A2D::Mat<TacsScalar, 3, 3> u0d;
  A2D::MatFromThreeVec(u0xi, d01, d02, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  A2D::Mat<TacsScalar, 3, 3> u0dXdinvT, u0x;
  A2D::MatMatMult(u0d, XdinvT, u0dXdinvT);
  A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT, u0x);

  // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
  // and sz2 = e1^{T} * Xdinv * n2xi
  const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
  A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
  TacsScalar s0, sz1, sz2;
  A2D::MatVecMult(XdinvT, e1, ts0);
  A2D::VecDot(e1, ts0, s0);
  A2D::MatVecMult(Xdinv, n1xi, tsz1);
  A2D::VecDot(e1, tsz1, sz1);
  A2D::MatVecMult(Xdinv, n2xi, tsz2);
  A2D::VecDot(e1, tsz2, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  A2D::Vec<TacsScalar, 3> d1t, d1s, d1x;
  A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t);
  A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s);
  A2D::VecScale(s0, d1s, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  A2D::Vec<TacsScalar, 3> d2t, d2s, d2x;
  A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t);
  A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s);
  A2D::VecScale(s0, d2s, d2x);

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
  e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x), A2D::get_data(d2x),
                    e0ty, e);

  if (quantityType == TACS_FAILURE_INDEX) {
    if (quantity) {
      *quantity = con->evalFailure(elemIndex, pt, A2D::get_data(X0), e);
    }
    return 1;
  }

  if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      TacsScalar s[6];
      con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);
      *quantity = 0.0;
      for (int i = 0; i < 6; i++) {
        *quantity += e[i] * s[i];
      }
    }
    return 1;
  }

  return 0;
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
    addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                           TacsScalar scale, int n, double pt[],
                           const TacsScalar Xpts[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[],
                           const TacsScalar dfdq[], int dvLen,
                           TacsScalar dfdx[]) {
  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn1,
                                                            d1, d1dot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn2,
                                                            d2, d2dot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  // The transformation to the local beam coordinates
  A2D::Mat<TacsScalar, 3, 3> T;

  // Parametric location
  A2D::Vec<TacsScalar, 3> X0;

  // Tangent to the beam
  A2D::Vec<TacsScalar, 3> X0xi;

  // Interpolated normal directions
  A2D::Vec<TacsScalar, 3> n1, n2;

  // Derivatives of the interpolated normal directions
  A2D::Vec<TacsScalar, 3> n1xi, n2xi;

  // The values of the director fields and their derivatives
  A2D::Vec<TacsScalar, 3> u0xi, d01, d02, d01xi, d02xi;

  // Interpolate the solution fields
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                     A2D::get_data(u0xi));
  basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
  basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
  basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
  basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

  // Compute X, X,xi and the interpolated normal
  basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
  basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
  basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
  basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
  basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

  // Compute the transformation at the quadrature point
  transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

  // Compute the inverse
  A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
  A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
  A2D::MatInv(Xd, Xdinv);

  // Compute the determinant of the transform
  TacsScalar detXd;
  A2D::MatDet(Xd, detXd);

  if (quantityType == TACS_ELEMENT_DENSITY) {
    con->addDensityDVSens(elemIndex, scale * dfdq[0], pt, A2D::get_data(X0),
                          dvLen, dfdx);
    return;
  }

  // Compute XdinvT = Xdinv * T
  A2D::Mat<TacsScalar, 3, 3> XdinvT;
  A2D::MatMatMult(Xdinv, T, XdinvT);

  // Assemble u0d
  A2D::Mat<TacsScalar, 3, 3> u0d;
  A2D::MatFromThreeVec(u0xi, d01, d02, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  A2D::Mat<TacsScalar, 3, 3> u0dXdinvT, u0x;
  A2D::MatMatMult(u0d, XdinvT, u0dXdinvT);
  A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT, u0x);

  // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
  // and sz2 = e1^{T} * Xdinv * n2xi
  const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
  A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
  TacsScalar s0, sz1, sz2;
  A2D::MatVecMult(XdinvT, e1, ts0);
  A2D::VecDot(e1, ts0, s0);
  A2D::MatVecMult(Xdinv, n1xi, tsz1);
  A2D::VecDot(e1, tsz1, sz1);
  A2D::MatVecMult(Xdinv, n2xi, tsz2);
  A2D::VecDot(e1, tsz2, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  A2D::Vec<TacsScalar, 3> d1t, d1s, d1x;
  A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t);
  A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s);
  A2D::VecScale(s0, d1s, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  A2D::Vec<TacsScalar, 3> d2t, d2s, d2x;
  A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t);
  A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s);
  A2D::VecScale(s0, d2s, d2x);

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
  e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x), A2D::get_data(d2x),
                    e0ty, e);

  if (quantityType == TACS_FAILURE_INDEX) {
    // Add the sensitivity contribution from the design variables
    con->addFailureDVSens(elemIndex, dfdq[0] * scale, pt, A2D::get_data(X0), e,
                          dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);
    con->addStressDVSens(elemIndex, scale * dfdq[0], pt, A2D::get_data(X0), e,
                         e, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar dfdmom[6] = {0.0};

    for (int i = 0; i < 3; i++) {
      dfdmom[0] += scale * dfdq[i] * A2D::get_data(X0)[i];
      dfdmom[1] += -scale * dfdq[i] * A2D::get_data(n1)[i];
      dfdmom[2] += -scale * dfdq[i] * A2D::get_data(n2)[i];
    }

    con->addMassMomentsDVSens(elemIndex, pt, A2D::get_data(X0), dfdmom, dvLen,
                              dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar dfdI0[6] = {0.0};

    // Evaluate the self MOI
    TacsScalar moments[6];
    con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), moments);
    TacsScalar density = moments[0];

    TacsScalar dfdmoments[6] = {0.0};
    mat3x3SymmTransformSens(A2D::get_data(T), dfdq, dfdI0);
    dfdmoments[0] = scale *
                    (moments[2] * moments[2] * dfdI0[3] +
                     -moments[1] * moments[2] * dfdI0[4] +
                     moments[1] * moments[1] * dfdI0[5]) /
                    density / density;
    dfdmoments[1] = scale *
                    (-2.0 * moments[1] * dfdI0[5] + moments[2] * dfdI0[4]) /
                    density;
    dfdmoments[2] = scale *
                    (-2.0 * moments[2] * dfdI0[3] + moments[1] * dfdI0[4]) /
                    density;
    dfdmoments[3] = scale * dfdI0[5];
    dfdmoments[4] = scale * dfdI0[3];
    dfdmoments[5] = -scale * dfdI0[4];

    TacsScalar dXcg[3], dXcgdrho[3], dXcgdmom1[3], dXcgdmom2[3];
    for (int i = 0; i < 3; i++) {
      dXcg[i] = A2D::get_data(X0)[i] - (moments[1] * A2D::get_data(n1)[i] +
                                        moments[2] * A2D::get_data(n2)[i]) /
                                           density;
      dXcgdrho[i] = (moments[1] * A2D::get_data(n1)[i] +
                     moments[2] * A2D::get_data(n2)[i]) /
                    density / density;
      dXcgdmom2[i] = -A2D::get_data(n2)[i] / density;
      dXcgdmom1[i] = -A2D::get_data(n1)[i] / density;
    }

    // Use parallel axis theorem to move MOI to origin
    dfdmoments[0] +=
        scale * dfdq[0] *
        (dXcg[1] * dXcg[1] + dXcg[2] * dXcg[2] +
         2.0 * density * (dXcg[1] * dXcgdrho[1] + dXcg[2] * dXcgdrho[2]));
    dfdmoments[0] -= scale * dfdq[1] *
                     (dXcg[0] * dXcg[1] + density * (dXcgdrho[0] * dXcg[1] +
                                                     dXcg[0] * dXcgdrho[1]));
    dfdmoments[0] -= scale * dfdq[2] *
                     (dXcg[0] * dXcg[2] + density * (dXcgdrho[0] * dXcg[2] +
                                                     dXcg[0] * dXcgdrho[2]));
    dfdmoments[0] +=
        scale * dfdq[3] *
        (dXcg[0] * dXcg[0] + dXcg[2] * dXcg[2] +
         2.0 * density * (dXcg[0] * dXcgdrho[0] + dXcg[2] * dXcgdrho[2]));
    dfdmoments[0] -= scale * dfdq[4] *
                     (dXcg[2] * dXcg[1] + density * (dXcgdrho[1] * dXcg[2] +
                                                     dXcg[1] * dXcgdrho[2]));
    dfdmoments[0] +=
        scale * dfdq[5] *
        (dXcg[0] * dXcg[0] + dXcg[1] * dXcg[1] +
         2.0 * density * (dXcg[0] * dXcgdrho[0] + dXcg[1] * dXcgdrho[1]));

    dfdmoments[1] += scale * dfdq[0] * density * 2.0 *
                     (dXcg[1] * dXcgdmom1[1] + dXcg[2] * dXcgdmom1[2]);
    dfdmoments[1] += -scale * dfdq[1] * density *
                     (dXcgdmom1[0] * dXcg[1] + dXcg[0] * dXcgdmom1[1]);
    dfdmoments[1] += -scale * dfdq[2] * density *
                     (dXcgdmom1[0] * dXcg[2] + dXcg[0] * dXcgdmom1[2]);
    dfdmoments[1] += scale * dfdq[3] * density * 2.0 *
                     (dXcg[0] * dXcgdmom1[0] + dXcg[2] * dXcgdmom1[2]);
    dfdmoments[1] += -scale * dfdq[4] * density *
                     (dXcgdmom1[2] * dXcg[1] + dXcg[2] * dXcgdmom1[1]);
    dfdmoments[1] += scale * dfdq[5] * density * 2.0 *
                     (dXcg[0] * dXcgdmom1[0] + dXcg[1] * dXcgdmom1[1]);

    dfdmoments[2] += scale * dfdq[0] * density * 2.0 *
                     (dXcg[1] * dXcgdmom2[1] + dXcg[2] * dXcgdmom2[2]);
    dfdmoments[2] += -scale * dfdq[1] * density *
                     (dXcgdmom2[0] * dXcg[1] + dXcg[0] * dXcgdmom2[1]);
    dfdmoments[2] += -scale * dfdq[2] * density *
                     (dXcgdmom2[0] * dXcg[2] + dXcg[0] * dXcgdmom2[2]);
    dfdmoments[2] += scale * dfdq[3] * density * 2.0 *
                     (dXcg[0] * dXcgdmom2[0] + dXcg[2] * dXcgdmom2[2]);
    dfdmoments[2] += -scale * dfdq[4] * density *
                     (dXcgdmom2[2] * dXcg[1] + dXcg[2] * dXcgdmom2[1]);
    dfdmoments[2] += scale * dfdq[5] * density * 2.0 *
                     (dXcg[0] * dXcgdmom2[0] + dXcg[1] * dXcgdmom2[1]);

    con->addMassMomentsDVSens(elemIndex, pt, A2D::get_data(X0), dfdmoments,
                              dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
    addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                           TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                           int n, double pt[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], const TacsScalar dfdq[],
                           TacsScalar dfdu[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Get the reference axis
    const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

    // Compute the normal directions
    TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
    TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

    // Compute the frame normal and directors at each node
    TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn1,
                                                              d1, d1dot);
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn2,
                                                              d2, d2dot);

    // Add the contributions to the derivative
    TacsScalar d1d[dsize], d2d[dsize];
    memset(d1d, 0, dsize * sizeof(TacsScalar));
    memset(d2d, 0, dsize * sizeof(TacsScalar));

    // Compute the tying strain values
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2,
                                                             vars, d1, d2, ety);

    TacsScalar dety[basis::NUM_TYING_POINTS];
    memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    // The transformation to the local beam coordinates
    A2D::Mat<TacsScalar, 3, 3> T;

    // Parametric location
    A2D::Vec<TacsScalar, 3> X0;

    // Tangent to the beam
    A2D::Vec<TacsScalar, 3> X0xi;

    // Interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::Vec<TacsScalar, 3> n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0xi, d01, d02, d01xi, d02xi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                       A2D::get_data(u0xi));
    basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
    basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
    basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
    basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
    basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
    basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
    basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
    basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

    // Compute the transformation at the quadrature point
    transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

    // Compute the inverse
    A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
    A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
    A2D::MatInv(Xd, Xdinv);

    // Compute the determinant of the transform
    TacsScalar detXd;
    A2D::MatDet(Xd, detXd);

    // Compute XdinvT = Xdinv * T
    A2D::Mat<TacsScalar, 3, 3> XdinvT;
    A2D::MatMatMult(Xdinv, T, XdinvT);

    // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
    // and sz2 = e1^{T} * Xdinv * n2xi
    const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
    A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
    TacsScalar s0, sz1, sz2;
    A2D::MatVecMult(XdinvT, e1, ts0);
    A2D::VecDot(e1, ts0, s0);
    A2D::MatVecMult(Xdinv, n1xi, tsz1);
    A2D::VecDot(e1, tsz1, sz1);
    A2D::MatVecMult(Xdinv, n2xi, tsz2);
    A2D::VecDot(e1, tsz2, sz2);

    // Assemble u0d and compute u0x = T^{T} * u0d * XdinvT,
    // d1x = s0 * T^{T} * (d1xi - sz1 * u0xi) and
    // d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> u0d, u0dXdinvT, u0x;
    A2D::ADObj<A2D::Vec<TacsScalar, 3>> d1t, d1s, d1x, d2t, d2s, d2x;
    auto strain_stack = A2D::MakeStack(
        A2D::MatFromThreeVec(u0xi, d01, d02, u0d),
        A2D::MatMatMult(u0d, XdinvT, u0dXdinvT),
        A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                                   u0x),
        A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s),
        A2D::VecScale(s0, d1s, d1x),
        A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t),
        A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s),
        A2D::VecScale(s0, d2s, d2x));

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2];
    e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
    e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

    // Evaluate the strain
    TacsScalar e[6];
    model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x),
                      A2D::get_data(d2x), e0ty, e);

    TacsScalar esens[6];
    if (quantityType == TACS_FAILURE_INDEX) {
      // Compute the sensitivity of the failure index w.r.t. the strain
      con->evalFailureStrainSens(elemIndex, pt, A2D::get_data(X0), e, esens);
    } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
      // Compute the sensitivity of the strain energy density w.r.t. the strain
      con->evalStress(elemIndex, pt, A2D::get_data(X0), e, esens);
      for (int i = 0; i < 6; i++) {
        esens[i] *= 2.0;
      }
    }

    // Evaluate the strain and strain derivatives from the
    TacsScalar e0tyd[2];
    model::evalStrainSens(
        alpha * dfdq[0], esens, A2D::get_data(u0x), A2D::get_data(d1x),
        A2D::get_data(d2x), e0ty, A2D::get_data(u0x.bvalue()),
        A2D::get_data(d1x.bvalue()), A2D::get_data(d2x.bvalue()), e0tyd);

    // Convert the contributions to the tying strain
    TacsScalar gtyd[2];
    gtyd[0] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[0];
    gtyd[1] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[1];

    strain_stack.reverse();

    // Add the residual contributions back to the element
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(
        pt, A2D::get_data(u0xi.bvalue()), dfdu);

    // Add the constributions back to the derivative
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d01.bvalue()), d1d);
    basis::template addInterpFieldsTranspose<3, 3>(
        pt, A2D::get_data(d02.bvalue()), d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d01xi.bvalue()), d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(
        pt, A2D::get_data(d02xi.bvalue()), d2d);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, gtyd, dety);

    // Add the contributions from the tying strain
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn1, fn2, vars, d1, d2, dety, dfdu, d1d, d2d);

    // Add the contributions to the director field
    director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
        vars, dvars, ddvars, fn1, d1d, dfdu);
    director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
        vars, dvars, ddvars, fn2, d2d, dfdu);
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    // Compute the interpolated displacements
    TacsScalar scale[3];
    scale[0] = alpha * dfdq[0];
    scale[1] = alpha * dfdq[1];
    scale[2] = alpha * dfdq[2];
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, scale, dfdu);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
    addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                            TacsScalar scale, int n, double pt[],
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[], const TacsScalar ddvars[],
                            const TacsScalar dfddetXd, const TacsScalar dfdq[],
                            TacsScalar dfdXpts[]) {
  // Get the reference axis
  const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

  // Compute the normal directions
  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Derivatives w.r.t. the frame normals
  TacsScalar dfn1[3 * basis::NUM_NODES], dfn2[3 * basis::NUM_NODES];
  memset(dfn1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
  memset(dfn2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn1,
                                                            d1, d1dot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(vars, dvars, fn2,
                                                            d2, d2dot);

  // Derivatives w.r.t. the d1 and d2 fields
  TacsScalar dd1[dsize], dd2[dsize];
  memset(dd1, 0, dsize * sizeof(TacsScalar));
  memset(dd2, 0, dsize * sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  TacsScalar etyd[basis::NUM_TYING_POINTS];
  memset(etyd, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // The transformation to the local beam coordinates
  A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> T;

  // Parametric location
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0;

  // Tangent to the beam
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> X0xi;

  // Interpolated normal directions
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> n1, n2;

  // Derivatives of the interpolated normal directions
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> n1xi, n2xi;

  // The values of the director fields and their derivatives
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> u0xi, d01, d02, d01xi, d02xi;

  // Compute X, X,xi and the interpolated normal
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                     A2D::get_data(u0xi));
  basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
  basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
  basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
  basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
  basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));
  basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
  basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
  basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
  basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

  // Compute the transformation at the quadrature point
  transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

  // Assemble Xd, compute its inverse and determinant, compute
  // XdinvT = Xdinv * T, u0x = T^{T} * u0d * XdinvT, the scale factors
  // s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi and
  // sz2 = e1^{T} * Xdinv * n2xi, and dkx = s0 * T^{T} * (dkxi - szk * u0xi)
  // for k = 1, 2
  A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> Xd, Xdinv, XdinvT;
  A2D::ADObj<TacsScalar> detXd;
  A2D::ADObj<A2D::Mat<TacsScalar, 3, 3>> u0d, u0dXdinvT, u0x;
  const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
  A2D::Vec<TacsScalar, 3> e1(e1_data);
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> ts0, tsz1, tsz2;
  A2D::ADObj<TacsScalar> s0, sz1, sz2;
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> w1, d1t, d1s, d1x;
  A2D::ADObj<A2D::Vec<TacsScalar, 3>> w2, d2t, d2s, d2x;

  auto strain_stack = A2D::MakeStack(
      A2D::MatFromThreeVec(X0xi, n1, n2, Xd), A2D::MatInv(Xd, Xdinv),
      A2D::MatDet(Xd, detXd), A2D::MatMatMult(Xdinv, T, XdinvT),
      A2D::MatFromThreeVec(u0xi, d01, d02, u0d),
      A2D::MatMatMult(u0d, XdinvT, u0dXdinvT),
      A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                                 u0x),
      A2D::MatVecMult(XdinvT, e1, ts0), A2D::VecDot(e1, ts0, s0),
      A2D::MatVecMult(Xdinv, n1xi, tsz1), A2D::VecDot(e1, tsz1, sz1),
      A2D::MatVecMult(Xdinv, n2xi, tsz2), A2D::VecDot(e1, tsz2, sz2),
      A2D::VecScale(sz1, u0xi, w1),
      A2D::VecSum(TacsScalar(1.0), d01xi, TacsScalar(-1.0), w1, d1t),
      A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s),
      A2D::VecScale(s0, d1s, d1x), A2D::VecScale(sz2, u0xi, w2),
      A2D::VecSum(TacsScalar(1.0), d02xi, TacsScalar(-1.0), w2, d2t),
      A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s),
      A2D::VecScale(s0, d2s, d2x));

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * A2D::get_data(XdinvT)[0] * gty[0];
  e0ty[1] = 2.0 * A2D::get_data(XdinvT)[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x), A2D::get_data(d2x),
                    e0ty, e);

  // Evaluate the failure sensitivity contribution
  TacsScalar esens[6] = {0.0};
  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the sensitivity of the failure index w.r.t. the strain
    con->evalFailureStrainSens(elemIndex, pt, A2D::get_data(X0), e, esens);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the sensitivity of the strain energy density w.r.t. the strain
    con->evalStress(elemIndex, pt, A2D::get_data(X0), e, esens);
    for (int i = 0; i < 6; i++) {
      esens[i] *= 2.0;
    }
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    // Compute the sensitivity of the strain energy density w.r.t. the strain
    TacsScalar mass_moment[6];
    con->evalMassMoments(elemIndex, pt, A2D::get_data(X0), mass_moment);
    TacsScalar density = mass_moment[0];

    for (int i = 0; i < 3; i++) {
      A2D::get_data(X0.bvalue())[i] = density * dfdq[i];
      A2D::get_data(n1.bvalue())[i] = mass_moment[2] * dfdq[i];
      A2D::get_data(n2.bvalue())[i] = -mass_moment[1] * dfdq[i];
    }
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TACSElement::addPointQuantityXptSens(elemIndex, quantityType, time, scale,
                                         n, pt, Xpts, vars, dvars, ddvars,
                                         dfddetXd, dfdq, dfdXpts);
    return;
  }

  // Evaluate the strain and strain derivatives from the
  TacsScalar e0tyd[2];
  model::evalStrainSens(
      scale * dfdq[0], esens, A2D::get_data(u0x), A2D::get_data(d1x),
      A2D::get_data(d2x), e0ty, A2D::get_data(u0x.bvalue()),
      A2D::get_data(d1x.bvalue()), A2D::get_data(d2x.bvalue()), e0tyd);
  detXd.bvalue() = scale * dfddetXd;

  // Apply the tying strain transformation
  TacsScalar gtyd[2];
  gtyd[0] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[0];
  gtyd[1] = 2.0 * A2D::get_data(XdinvT)[0] * e0tyd[1];

  A2D::get_data(XdinvT.bvalue())[0] +=
      2.0 * (gty[0] * e0tyd[0] + gty[1] * e0tyd[1]);

  strain_stack.reverse();

  // Reverse the transformation sensitivities
  transform->addTransformSens(A2D::get_data(X0xi), A2D::get_data(T.bvalue()),
                              A2D::get_data(X0xi.bvalue()));

  // Add the sensitivities to the input fields...
  basis::template addInterpFieldsTranspose<3, 3>(pt, A2D::get_data(X0.bvalue()),
                                                 dfdXpts);
  basis::template addInterpFieldsGradTranspose<3, 3>(
      pt, A2D::get_data(X0xi.bvalue()), dfdXpts);
  basis::template addInterpFieldsTranspose<3, 3>(pt, A2D::get_data(n1.bvalue()),
                                                 dfn1);
  basis::template addInterpFieldsTranspose<3, 3>(pt, A2D::get_data(n2.bvalue()),
                                                 dfn2);
  basis::template addInterpFieldsGradTranspose<3, 3>(
      pt, A2D::get_data(n1xi.bvalue()), dfn1);
  basis::template addInterpFieldsGradTranspose<3, 3>(
      pt, A2D::get_data(n2xi.bvalue()), dfn2);
  basis::template addInterpFieldsTranspose<3, 3>(
      pt, A2D::get_data(d01.bvalue()), dd1);
  basis::template addInterpFieldsTranspose<3, 3>(
      pt, A2D::get_data(d02.bvalue()), dd2);
  basis::template addInterpFieldsGradTranspose<3, 3>(
      pt, A2D::get_data(d01xi.bvalue()), dd1);
  basis::template addInterpFieldsGradTranspose<3, 3>(
      pt, A2D::get_data(d02xi.bvalue()), dd2);

  // Add the derivative contributions to the tying strain
  basis::addInterpTyingStrainTranspose(pt, gtyd, etyd);

  // Add the sensitivity contributions from the tying strain
  model::template addTyingStrainXptSens<vars_per_node, basis>(
      Xpts, fn1, fn2, vars, d1, d2, etyd, dfdXpts, dfn1, dfn2, dd1, dd2);

  // vars, dvars, dd1, dd1dot -> varsd, dvarsd and dfn1
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(vars, fn1, dd1,
                                                                dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(vars, fn2, dd2,
                                                                dfn2);

  // Add the contributions from the node normals
  TacsBeamAddNodeNormalsSens<basis>(Xpts, axis, dfn1, dfn2, dfdXpts);
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::getOutputData(
    int elemIndex, ElementType etype, int write_flag, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int ld_data, TacsScalar *data) {
  if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    // Get the number of nodes associated with the visualization
    int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

    // Get the reference axis
    const A2D::Vec<TacsScalar, 3> &axis = transform->getRefAxis();

    // Compute the normal directions
    TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
    TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

    // Compute the frame normal and directors at each node
    TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn1,
                                                              d1, d1dot);
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn2,
                                                              d2, d2dot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2,
                                                             vars, d1, d2, ety);

    // Loop over each quadrature point and add the residual contribution
    for (int index = 0; index < num_vis_nodes; index++) {
      // Get the quadrature weight
      double pt[3];
      basis::getNodePoint(index, pt);

      // The transformation to the local beam coordinates
      A2D::Mat<TacsScalar, 3, 3> T;

      // Parametric location
      A2D::Vec<TacsScalar, 3> X0;

      // Tangent to the beam
      A2D::Vec<TacsScalar, 3> X0xi;

      // Interpolated normal directions
      A2D::Vec<TacsScalar, 3> n1, n2;

      // Derivatives of the interpolated normal directions
      A2D::Vec<TacsScalar, 3> n1xi, n2xi;

      // The values of the director fields and their derivatives
      A2D::Vec<TacsScalar, 3> u0xi, d01, d02, d01xi, d02xi;

      // Interpolate the solution fields
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars,
                                                         A2D::get_data(u0xi));
      basis::template interpFields<3, 3>(pt, d1, A2D::get_data(d01));
      basis::template interpFields<3, 3>(pt, d2, A2D::get_data(d02));
      basis::template interpFieldsGrad<3, 3>(pt, d1, A2D::get_data(d01xi));
      basis::template interpFieldsGrad<3, 3>(pt, d2, A2D::get_data(d02xi));

      // Interpolate the geometry fields
      basis::template interpFields<3, 3>(pt, Xpts, A2D::get_data(X0));
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, A2D::get_data(X0xi));
      basis::template interpFields<3, 3>(pt, fn1, A2D::get_data(n1));
      basis::template interpFields<3, 3>(pt, fn2, A2D::get_data(n2));
      basis::template interpFieldsGrad<3, 3>(pt, fn1, A2D::get_data(n1xi));
      basis::template interpFieldsGrad<3, 3>(pt, fn2, A2D::get_data(n2xi));

      // Compute the transformation at the quadrature point
      transform->computeTransform(A2D::get_data(X0xi), A2D::get_data(T));

      // Compute the inverse
      A2D::Mat<TacsScalar, 3, 3> Xd, Xdinv;
      A2D::MatFromThreeVec(X0xi, n1, n2, Xd);
      A2D::MatInv(Xd, Xdinv);

      // Compute the determinant of the transform
      TacsScalar detXd;
      A2D::MatDet(Xd, detXd);

      // Compute XdinvT = Xdinv * T
      A2D::Mat<TacsScalar, 3, 3> XdinvT;
      A2D::MatMatMult(Xdinv, T, XdinvT);

      // Assemble u0d
      A2D::Mat<TacsScalar, 3, 3> u0d;
      A2D::MatFromThreeVec(u0xi, d01, d02, u0d);

      // Compute u0x = T^{T} * u0d * XdinvT
      A2D::Mat<TacsScalar, 3, 3> u0dXdinvT, u0x;
      A2D::MatMatMult(u0d, XdinvT, u0dXdinvT);
      A2D::MatMatMult<A2D::MatOp::TRANSPOSE, A2D::MatOp::NORMAL>(T, u0dXdinvT,
                                                                 u0x);

      // Compute s0 = e1^{T} * XdinvT * e1, sz1 = e1^{T} * Xdinv * n1xi
      // and sz2 = e1^{T} * Xdinv * n2xi
      const TacsScalar e1_data[3] = {1.0, 0.0, 0.0};
      A2D::Vec<TacsScalar, 3> e1(e1_data), ts0, tsz1, tsz2;
      TacsScalar s0, sz1, sz2;
      A2D::MatVecMult(XdinvT, e1, ts0);
      A2D::VecDot(e1, ts0, s0);
      A2D::MatVecMult(Xdinv, n1xi, tsz1);
      A2D::VecDot(e1, tsz1, sz1);
      A2D::MatVecMult(Xdinv, n2xi, tsz2);
      A2D::VecDot(e1, tsz2, sz2);

      // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
      A2D::Vec<TacsScalar, 3> d1t, d1s, d1x;
      A2D::VecSum(-sz1, u0xi, TacsScalar(1.0), d01xi, d1t);
      A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d1t, d1s);
      A2D::VecScale(s0, d1s, d1x);

      // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
      A2D::Vec<TacsScalar, 3> d2t, d2s, d2x;
      A2D::VecSum(-sz2, u0xi, TacsScalar(1.0), d02xi, d2t);
      A2D::MatVecMult<A2D::MatOp::TRANSPOSE>(T, d2t, d2s);
      A2D::VecScale(s0, d2s, d2x);

      // Evaluate the tying components of the strain
      TacsScalar gty[2];  // The components of the tying strain
      basis::interpTyingStrain(pt, ety, gty);

      // Transform the tying strain to the local coordinates
      TacsScalar e0ty[2];
      e0ty[0] = gty[0];
      e0ty[1] = gty[1];

      // Compute the set of strain components
      TacsScalar e[6];  // The components of the strain
      model::evalStrain(A2D::get_data(u0x), A2D::get_data(d1x),
                        A2D::get_data(d2x), e0ty, e);

      // Compute the corresponding stresses
      TacsScalar s[6];
      con->evalStress(elemIndex, pt, A2D::get_data(X0), e, s);

      if (write_flag & TACS_OUTPUT_NODES) {
        data[0] = A2D::get_data(X0)[0];
        data[1] = A2D::get_data(X0)[1];
        data[2] = A2D::get_data(X0)[2];
        data += 3;
      }
      if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
        int len = vars_per_node;
        if (len > 6) {
          len = 6;
        }
        for (int i = 0; i < len; i++) {
          data[i] = vars[i + vars_per_node * index];
        }
        for (int i = len; i < 6; i++) {
          data[i] = 0.0;
        }
        data += 6;
      }
      if (write_flag & TACS_OUTPUT_STRAINS) {
        for (int i = 0; i < 6; i++) {
          data[i] = e[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_STRESSES) {
        for (int i = 0; i < 6; i++) {
          data[i] = s[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_EXTRAS) {
        for (int failInd = 0; failInd < 7; failInd++) {
          data[failInd] = con->evalFailureFieldValue(
              elemIndex, pt, A2D::get_data(X0), e, failInd);
        }
        for (int dvInd = 0; dvInd < 7; dvInd++) {
          data[dvInd + 7] = con->evalDesignFieldValue(elemIndex, pt,
                                                      A2D::get_data(X0), dvInd);
        }
        data += 14;
      }
      if (write_flag & TACS_OUTPUT_COORDINATE_FRAME) {
        data[0] = A2D::get_data(T)[0];
        data[1] = A2D::get_data(T)[3];
        data[2] = A2D::get_data(T)[6];

        data[3] = A2D::get_data(T)[1];
        data[4] = A2D::get_data(T)[4];
        data[5] = A2D::get_data(T)[7];

        data[6] = A2D::get_data(T)[2];
        data[7] = A2D::get_data(T)[5];
        data[8] = A2D::get_data(T)[8];
        data += 9;
      }
    }
  }
}

#endif  // TACS_BEAM_ELEMENT_H
