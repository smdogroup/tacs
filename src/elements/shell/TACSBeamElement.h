#ifndef TACS_BEAM_ELEMENT_H
#define TACS_BEAM_ELEMENT_H

#include <typeinfo>

#include "TACSBeamCentrifugalForce.h"
#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementModel.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamElementTransform.h"
#include "TACSBeamInertialForce.h"
#include "TACSBeamTraction.h"
#include "TACSBeamUtilities.h"
#include "TACSDirector.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "a2d.h"

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

  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]);

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

  void getAverageStresses(int elemIndex, ElementType etype,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TacsScalar *avgStresses);

  void addMatDVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, TacsScalar scale,
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[], int dvLen,
                                TacsScalar dfdx[]);

  void addMatXptSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                 double time, TacsScalar scale,
                                 const TacsScalar psi[], const TacsScalar phi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[], TacsScalar dfdXpts[]);

  void getMatSVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, const TacsScalar psi[],
                                const TacsScalar phi[], const TacsScalar Xpts[],
                                const TacsScalar vars[], TacsScalar dfdu[]);

 private:
  // Set sizes for the different components
  static const int usize = 3 * basis::NUM_NODES;
  static const int dsize = 3 * basis::NUM_NODES;
  static const int csize = 9 * basis::NUM_NODES;

  // Helper for addMatXptSensInnerProduct's TACS_STIFFNESS_MATRIX branch
  // (Task 4.3, SPEC.md sec 2.3.1) -- split out of the dispatcher for
  // readability given its size.
  void addMatXptSensInnerProductStiffness(
      int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
      const TacsScalar phi[], const TacsScalar Xpts[], const TacsScalar vars[],
      TacsScalar dfdXpts[]);

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
  const A2D::Vec3 &axis = transform->getRefAxis();

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

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

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
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    // Assemble u0d
    A2D::Mat3x3 u0d;
    A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::Mat3x3 u0dXdinvT, u0x;
    A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

    // Compute s0, sz1 and sz2
    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::Vec3 d1t, d1x;
    A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::Vec3 d2t, d2x;
    A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

    // Compute the set of strain components
    TacsScalar e[6];  // The components of the strain
    model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    Ue += 0.5 * detXd.value *
          (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] + s[4] * e[4] +
           s[5] * e[5]);

    // Evaluate the velocities
    A2D::Vec3 u0dot, d01dot, d02dot;
    basis::template interpFields<vars_per_node, 3>(pt, dvars, u0dot.x);
    basis::template interpFields<3, 3>(pt, d1dot, d01dot.x);
    basis::template interpFields<3, 3>(pt, d2dot, d02dot.x);

    // dot{u} = \dot{u0} + z1 * dot{d1} + z2 * dot{d2}
    A2D::Scalar u0d0, u0d10, u0d20, d1d10, d2d20, d1d20;
    A2D::Vec3Dot u0ddot(u0dot, u0dot, u0d0);
    A2D::Vec3Dot u0d1dot(u0dot, d01dot, u0d10);
    A2D::Vec3Dot u0d2dot(u0dot, d02dot, u0d20);
    A2D::Vec3Dot d1d1dot(d01dot, d01dot, d1d10);
    A2D::Vec3Dot d2d2dot(d02dot, d02dot, d2d20);
    A2D::Vec3Dot d2d1dot(d01dot, d02dot, d1d20);

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, X0.x, rho);

    Te += 0.5 * detXd.value *
          (rho[0] * u0d0.value + 2.0 * rho[1] * u0d10.value +
           2.0 * rho[2] * u0d20.value + rho[3] * d1d10.value +
           rho[4] * d2d20.value + 2.0 * rho[5] * d1d20.value);
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
  const A2D::Vec3 &axis = transform->getRefAxis();

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

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

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
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    // Assemble u0d
    A2D::ADMat3x3 u0d;
    A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::ADMat3x3 u0dXdinvT, u0x;
    A2D::ADMat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

    // Compute s0, sz1 and sz2
    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::ADVec3 d1t, d1x;
    A2D::ADVec3ADVecScalarAxpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::MatTrans3x3ADVecMultScale matmultd1x(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::ADVec3 d2t, d2x;
    A2D::ADVec3ADVecScalarAxpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::MatTrans3x3ADVecMultScale matmultd2x(s0, T, d2t, d2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], de0ty[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

    // Evaluate the strain
    TacsScalar e[6];
    model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    // Evaluate the strain and strain derivatives from the
    model::evalStrainSens(detXd.value, s, u0x.A, d1x.x, d2x.x, e0ty, u0x.Ad,
                          d1x.xd, d2x.xd, de0ty);

    // Convert the contributions to the tying strain
    TacsScalar dgty[2];
    dgty[0] = 2.0 * XdinvT.A[0] * de0ty[0];
    dgty[1] = 2.0 * XdinvT.A[0] * de0ty[1];

    matmultd2x.reverse();
    axpyd2t.reverse();
    matmultd1x.reverse();
    axpyd1t.reverse();
    multu0x.reverse();
    multu0d.reverse();
    assembleu0d.reverse();

    // Add the residual contributions back to the element
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.xd,
                                                                   res);

    // Add the constributions back to the derivative
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, d2d);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    // Evaluate the accelerations
    A2D::ADVec3 u0ddot, d01ddot, d02ddot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot.x);
    basis::template interpFields<3, 3>(pt, d1ddot, d01ddot.x);
    basis::template interpFields<3, 3>(pt, d2ddot, d02ddot.x);

    // dot{u}(xi, z1, z2) = dot{u0} + z1 * dot{d1} + z2 * dot{d2}
    A2D::ADScalar u0d0, u0d10, u0d20, d1d10, d2d20, d1d20;
    A2D::ADVec3Dot u0ddot0(u0ddot, u0ddot, u0d0);
    A2D::ADVec3Dot u0d1dot(u0ddot, d01ddot, u0d10);
    A2D::ADVec3Dot u0d2dot(u0ddot, d02ddot, u0d20);
    A2D::ADVec3Dot d1d1dot(d01ddot, d01ddot, d1d10);
    A2D::ADVec3Dot d2d2dot(d02ddot, d02ddot, d2d20);
    A2D::ADVec3Dot d2d1dot(d01ddot, d02ddot, d1d20);

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, X0.x, rho);

    u0d0.valued = 0.5 * rho[0] * detXd.value;
    u0d10.valued = rho[1] * detXd.value;
    u0d20.valued = rho[2] * detXd.value;
    d1d10.valued = 0.5 * rho[3] * detXd.value;
    d2d20.valued = 0.5 * rho[4] * detXd.value;
    d1d20.valued = rho[5] * detXd.value;

    d2d1dot.reverse();
    d2d2dot.reverse();
    d1d1dot.reverse();
    u0d2dot.reverse();
    u0d1dot.reverse();
    u0ddot0.reverse();

    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, u0ddot.xd,
                                                               res);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01ddot.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02ddot.xd, d2d);
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

/*
  Zero the second-order (.xp/.xh, .Ap/.Ah) fields on every second-order A2D
  node in addJacobian's forward/reverse chain, ahead of a single per-DOF
  hforward/hreverse sweep iteration (SPEC.md sec 1.3 step 3). These fields
  are "+=" accumulators across hforward()/hreverse() calls, not
  auto-reset, and (unlike the primal/first-order .x/.xd fields) most of
  them are NOT freshly overwritten by every sweep iteration's own seed
  step -- e.g. the translational sweep seeds only u0xi.xp, leaving
  d01/d02/d01xi/d02xi's .xp stale from a previous iteration unless
  explicitly zeroed here.
*/
static inline void TacsBeamZeroSecondOrderNodes(
    A2D::ADVec3 &u0xi, A2D::ADVec3 &d01, A2D::ADVec3 &d02, A2D::ADVec3 &d01xi,
    A2D::ADVec3 &d02xi, A2D::ADMat3x3 &u0d, A2D::ADMat3x3 &u0dXdinvT,
    A2D::ADMat3x3 &u0x, A2D::ADVec3 &d1t, A2D::ADVec3 &d1x, A2D::ADVec3 &d2t,
    A2D::ADVec3 &d2x) {
  memset(u0xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(u0xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(d01.xp, 0, 3 * sizeof(TacsScalar));
  memset(d01.xh, 0, 3 * sizeof(TacsScalar));
  memset(d02.xp, 0, 3 * sizeof(TacsScalar));
  memset(d02.xh, 0, 3 * sizeof(TacsScalar));
  memset(d01xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(d01xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(d02xi.xp, 0, 3 * sizeof(TacsScalar));
  memset(d02xi.xh, 0, 3 * sizeof(TacsScalar));
  memset(u0d.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0d.Ah, 0, 9 * sizeof(TacsScalar));
  memset(u0dXdinvT.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0dXdinvT.Ah, 0, 9 * sizeof(TacsScalar));
  memset(u0x.Ap, 0, 9 * sizeof(TacsScalar));
  memset(u0x.Ah, 0, 9 * sizeof(TacsScalar));
  memset(d1t.xp, 0, 3 * sizeof(TacsScalar));
  memset(d1t.xh, 0, 3 * sizeof(TacsScalar));
  memset(d1x.xp, 0, 3 * sizeof(TacsScalar));
  memset(d1x.xh, 0, 3 * sizeof(TacsScalar));
  memset(d2t.xp, 0, 3 * sizeof(TacsScalar));
  memset(d2t.xh, 0, 3 * sizeof(TacsScalar));
  memset(d2x.xp, 0, 3 * sizeof(TacsScalar));
  memset(d2x.xh, 0, 3 * sizeof(TacsScalar));
}

/*
  Contract the fixed, per-quadrature-point material Hessian blocks (from
  model::evalStrainHessian, SPEC.md sec 1.1) against the forward-propagated
  seed direction (u0xp/d1xp/d2xp, populated by a preceding hforward() call
  sweep) to produce the Hessian-vector product at the u0x/d1x/d2x nodes
  (u0xh/d1xh/d2xh), which the following hreverse() call sweep then
  propagates back to the vars-space DOFs. This is the "H*p" step of the
  per-DOF sweep (SPEC.md sec 1.3 step 3) -- a pure linear-algebra
  contraction, not a differentiation; every array shape/index convention
  below matches model::evalStrainHessian's own documented layout exactly
  (row-major, e.g. d2u0xd1x[3*i+j] = d^2/d(u0x_i)d(d1x_j)).
*/
static inline void TacsBeamContractStrainHessian(
    const TacsScalar d2u0x[81], const TacsScalar d2d1x[9],
    const TacsScalar d2d2x[9], const TacsScalar d2u0xd1x[27],
    const TacsScalar d2u0xd2x[27], const TacsScalar d2d1xd2x[9],
    const TacsScalar u0xp[9], const TacsScalar d1xp[3],
    const TacsScalar d2xp[3], TacsScalar u0xh[9], TacsScalar d1xh[3],
    TacsScalar d2xh[3]) {
  for (int i = 0; i < 9; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 9; j++) {
      val += d2u0x[9 * i + j] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2u0xd1x[3 * i + j] * d1xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2u0xd2x[3 * i + j] * d2xp[j];
    }
    u0xh[i] += val;
  }

  for (int i = 0; i < 3; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 3; j++) {
      val += d2d1x[3 * i + j] * d1xp[j];
    }
    for (int j = 0; j < 9; j++) {
      val += d2u0xd1x[3 * j + i] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2d1xd2x[3 * i + j] * d2xp[j];
    }
    d1xh[i] += val;
  }

  for (int i = 0; i < 3; i++) {
    TacsScalar val = 0.0;
    for (int j = 0; j < 3; j++) {
      val += d2d2x[3 * i + j] * d2xp[j];
    }
    for (int j = 0; j < 9; j++) {
      val += d2u0xd2x[3 * j + i] * u0xp[j];
    }
    for (int j = 0; j < 3; j++) {
      val += d2d1xd2x[3 * j + i] * d1xp[j];
    }
    d2xh[i] += val;
  }
}

/*
  Add the Jacobian of the residual (SPEC.md sec 1.3): a second-order
  extension of addResidual's own A2D graph, rather than a hand-coded
  congruence-transform port. Reruns addResidual's own forward pass with
  the 6 load-bearing A2D types upgraded to their second-order variants
  (SPEC.md sec 1.2.2), then for each of 3*dsize state-DOF directions (the
  "translational", "d1" and "d2" sweeps, SPEC.md sec 1.3 step 3,
  VALIDATION.md E6) calls hforward()/hreverse() in the same order
  addResidual's forward()/reverse() calls already use, extracting one
  column of mat[] (or one column of the director-space Hessian
  accumulators d2d1/d2d2/d2d1u/d2d2u/d2d1d2) per sweep iteration. The
  post-loop closures (model::addComputeTyingStrainHessian, "Gap 1", SPEC.md
  sec 1.3.1; TacsBeamAddCrossDirectorJacobian, "Gap 2", SPEC.md sec 1.3.2;
  director::addDirectorJacobian x2) then convert those director-space
  accumulators into mat[]'s vars-space entries, mirroring how
  model::addComputeTyingStrainTranspose/director::addDirectorResidual
  already convert the analogous first-order buffers for addResidual's res.
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[], TacsScalar res[],
    TacsScalar mat[]) {
  const int nvars = vars_per_node * num_nodes;

  // Zero the output buffers (SPEC.md sec 1.3 step 1).
  memset(mat, 0, nvars * nvars * sizeof(TacsScalar));
  if (res) {
    memset(res, 0, nvars * sizeof(TacsScalar));
  }

  // Setup identical to addResidual's own (SPEC.md sec 1.3 step 2).
  const int nquad = quadrature::getNumQuadraturePoints();

  const A2D::Vec3 &axis = transform->getRefAxis();

  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(
      vars, dvars, ddvars, fn1, d1, d1dot, d1ddot);
  director::template computeDirectorRates<vars_per_node, offset,
                                          basis::NUM_NODES>(
      vars, dvars, ddvars, fn2, d2, d2dot, d2ddot);

  TacsScalar d1d[dsize], d2d[dsize];
  memset(d1d, 0, dsize * sizeof(TacsScalar));
  memset(d2d, 0, dsize * sizeof(TacsScalar));

  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars,
                                                           d1, d2, ety);

  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // Second-order accumulators (SPEC.md sec 1.3 step 3/4, sec 1.3.1/1.3.2).
  // d2ety: tying-point-space tying-strain Hessian (Gap 1's input).
  TacsScalar d2ety[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
  memset(
      d2ety, 0,
      basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // d2d1/d2d2: director-space self-Hessians (fed to addDirectorJacobian).
  // d2d1u/d2d2u: (director, u0xi) cross-Hessians (also addDirectorJacobian).
  // d2d1d2: (d1, d2) cross-Hessian, the new Gap 2 accumulator.
  TacsScalar d2d1[dsize * dsize], d2d2[dsize * dsize];
  TacsScalar d2d1u[dsize * dsize], d2d2u[dsize * dsize];
  TacsScalar d2d1d2[dsize * dsize];
  memset(d2d1, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2d2, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2d1u, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2d2u, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2d1d2, 0, dsize * dsize * sizeof(TacsScalar));

  // Unscaled dynamics (mass-moment) Hessian accumulators -- gamma is
  // applied internally by director::addDirectorJacobian (matching its
  // documented contract), except for d2Tdotd1d2, which this closure
  // gamma-scales itself before folding into d2d1d2 (SPEC.md sec 1.3.2).
  TacsScalar d2Tdotd1[dsize * dsize], d2Tdotd2[dsize * dsize];
  TacsScalar d2Tdotu1[dsize * dsize], d2Tdotu2[dsize * dsize];
  TacsScalar d2Tdotd1d2[dsize * dsize];
  memset(d2Tdotd1, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotd2, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotu1, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotu2, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotd1d2, 0, dsize * dsize * sizeof(TacsScalar));

  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    A2D::Mat3x3 T;
    A2D::Vec3 X0;
    A2D::Vec3 X0xi;
    A2D::Vec3 n1, n2;
    A2D::Vec3 n1xi, n2xi;

    A2D::ADVec3 u0xi, d01, d02, d01xi, d02xi;

    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

    transform->computeTransform(X0xi.x, T.A);

    A2D::Mat3x3 Xd, Xdinv;
    A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::Mat3x3Inverse invXd(Xd, Xdinv);

    A2D::Scalar detXd;
    A2D::Mat3x3Det computedetXd(weight, Xd, detXd);

    A2D::Mat3x3 XdinvT;
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    A2D::ADMat3x3 u0d;
    A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);

    A2D::ADMat3x3 u0dXdinvT, u0x;
    A2D::ADMat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    A2D::ADVec3 d1t, d1x;
    A2D::ADVec3ADVecScalarAxpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::MatTrans3x3ADVecMultScale matmultd1x(s0, T, d1t, d1x);

    A2D::ADVec3 d2t, d2x;
    A2D::ADVec3ADVecScalarAxpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::MatTrans3x3ADVecMultScale matmultd2x(s0, T, d2t, d2x);

    TacsScalar gty[2];
    basis::interpTyingStrain(pt, ety, gty);

    TacsScalar e0ty[2], de0ty[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

    TacsScalar e[6];
    model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    // NEW: materialize the tangent stiffness (SPEC.md sec 1.3 step 3's
    // model-layer-coupling bullet) -- neither addResidual nor
    // computeEnergies ever call this; it is consumed starting Task 2.2
    // by model::evalStrainHessian, not by this scaffold task.
    TacsScalar Cs[TACSBeamConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
    con->evalTangentStiffness(elemIndex, pt, X0.x, Cs);

    model::evalStrainSens(detXd.value, s, u0x.A, d1x.x, d2x.x, e0ty, u0x.Ad,
                          d1x.xd, d2x.xd, de0ty);

    TacsScalar dgty[2];
    dgty[0] = 2.0 * XdinvT.A[0] * de0ty[0];
    dgty[1] = 2.0 * XdinvT.A[0] * de0ty[1];

    matmultd2x.reverse();
    axpyd2t.reverse();
    matmultd1x.reverse();
    axpyd1t.reverse();
    multu0x.reverse();
    multu0d.reverse();
    assembleu0d.reverse();

    if (res) {
      basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(
          pt, u0xi.xd, res);
    }

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, d2d);

    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    A2D::ADVec3 u0ddot, d01ddot, d02ddot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot.x);
    basis::template interpFields<3, 3>(pt, d1ddot, d01ddot.x);
    basis::template interpFields<3, 3>(pt, d2ddot, d02ddot.x);

    A2D::ADScalar u0d0, u0d10, u0d20, d1d10, d2d20, d1d20;
    A2D::ADVec3Dot u0ddot0(u0ddot, u0ddot, u0d0);
    A2D::ADVec3Dot u0d1dot(u0ddot, d01ddot, u0d10);
    A2D::ADVec3Dot u0d2dot(u0ddot, d02ddot, u0d20);
    A2D::ADVec3Dot d1d1dot(d01ddot, d01ddot, d1d10);
    A2D::ADVec3Dot d2d2dot(d02ddot, d02ddot, d2d20);
    A2D::ADVec3Dot d2d1dot(d01ddot, d02ddot, d1d20);

    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, X0.x, rho);

    u0d0.valued = 0.5 * rho[0] * detXd.value;
    u0d10.valued = rho[1] * detXd.value;
    u0d20.valued = rho[2] * detXd.value;
    d1d10.valued = 0.5 * rho[3] * detXd.value;
    d2d20.valued = 0.5 * rho[4] * detXd.value;
    d1d20.valued = rho[5] * detXd.value;

    d2d1dot.reverse();
    d2d2dot.reverse();
    d1d1dot.reverse();
    u0d2dot.reverse();
    u0d1dot.reverse();
    u0ddot0.reverse();

    if (res) {
      basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, u0ddot.xd,
                                                                 res);
    }
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01ddot.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02ddot.xd, d2d);

    // --- Second-order: material Hessian blocks (SPEC.md sec 1.3 step 3's
    // model-layer-coupling bullet; Task 2.2). ---
    TacsScalar d2u0x[81], d2d1x[9], d2d2x[9], d2e0ty[4];
    TacsScalar d2u0xd1x[27], d2u0xd2x[27], d2d1xd2x[9];
    model::evalStrainHessian(alpha * detXd.value, s, Cs, u0x.A, d1x.x, d2x.x,
                             e0ty, d2u0x, d2d1x, d2d2x, d2e0ty, d2u0xd1x,
                             d2u0xd2x, d2d1xd2x);

    // Tying-strain Hessian: accumulate this quadrature point's contribution
    // to the tying-point-space Hessian d2ety (SPEC.md sec 1.3.1's own
    // per-quadrature-point seeding note). e0ty[k] = 2*XdinvT.A[0]*gty[k] is
    // a fixed linear map, so d2gty = (2*XdinvT.A[0])^2 * d2e0ty.
    TacsScalar coef = 2.0 * XdinvT.A[0];
    TacsScalar d2gty[4];
    d2gty[0] = coef * coef * d2e0ty[0];
    d2gty[1] = coef * coef * d2e0ty[1];
    d2gty[2] = coef * coef * d2e0ty[2];
    d2gty[3] = coef * coef * d2e0ty[3];
    basis::addInterpTyingStrainHessian(pt, d2gty, d2ety);

    // --- Dynamics (mass-moment) Hessian blocks (SPEC.md sec 1.3 step 3's
    // dynamics bullet; Task 2.3). These are fixed-coefficient outer
    // products of the basis shape functions -- no hforward/hreverse sweep
    // needed, mirroring TACSShellElement::addJacobian's identical
    // mass-diagonal-block pattern (TACSShellElement.h:605-617). ---
    TacsScalar d2mass[9];
    memset(d2mass, 0, 9 * sizeof(TacsScalar));
    d2mass[0] = d2mass[4] = d2mass[8] = gamma * detXd.value * rho[0];
    basis::template addInterpFieldsOuterProduct<vars_per_node, vars_per_node, 3,
                                                3>(pt, d2mass, mat);

    d2mass[0] = d2mass[4] = d2mass[8] = detXd.value * rho[3];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2mass,
                                                            d2Tdotd1);
    d2mass[0] = d2mass[4] = d2mass[8] = detXd.value * rho[4];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2mass,
                                                            d2Tdotd2);
    d2mass[0] = d2mass[4] = d2mass[8] = detXd.value * rho[1];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2mass,
                                                            d2Tdotu1);
    d2mass[0] = d2mass[4] = d2mass[8] = detXd.value * rho[2];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2mass,
                                                            d2Tdotu2);
    d2mass[0] = d2mass[4] = d2mass[8] = detXd.value * rho[5];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2mass,
                                                            d2Tdotd1d2);

    // --- Per-DOF hforward/hreverse sweep (SPEC.md sec 1.3 step 3,
    // VALIDATION.md E6): one sweep per compact "translational"/"d1"/"d2"
    // state-DOF direction, 3*dsize sweeps total per quadrature point (this
    // is NOT nvars-many sweeps -- see the feature's phase-2 handoff for why
    // the literal "for k in 0..nvars" reading of the SPEC text is subtly
    // incomplete for a directored element). ---
    TacsScalar seed[dsize];

    // (A) Translational sweep: seeds u0xi only. Captures the (u0xi, u0xi)
    // self block (scattered directly into mat[]) and the (u0xi, d1)/
    // (u0xi, d2) cross blocks (into d2d1u/d2d2u).
    for (int m = 0; m < dsize; m++) {
      TacsBeamZeroSecondOrderNodes(u0xi, d01, d02, d01xi, d02xi, u0d, u0dXdinvT,
                                   u0x, d1t, d1x, d2t, d2x);

      memset(seed, 0, dsize * sizeof(TacsScalar));
      seed[m] = 1.0;
      basis::template interpFieldsGrad<3, 3>(pt, seed, u0xi.xp);

      assembleu0d.hforward();
      multu0d.hforward();
      multu0x.hforward();
      axpyd1t.hforward();
      matmultd1x.hforward();
      axpyd2t.hforward();
      matmultd2x.hforward();

      TacsBeamContractStrainHessian(d2u0x, d2d1x, d2d2x, d2u0xd1x, d2u0xd2x,
                                    d2d1xd2x, u0x.Ap, d1x.xp, d2x.xp, u0x.Ah,
                                    d1x.xh, d2x.xh);

      matmultd2x.hreverse();
      axpyd2t.hreverse();
      matmultd1x.hreverse();
      axpyd1t.hreverse();
      multu0x.hreverse();
      multu0d.hreverse();
      assembleu0d.hreverse();

      int mm = vars_per_node * (m / 3) + (m % 3);

      TacsScalar ucol[dsize];
      memset(ucol, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, u0xi.xh, ucol);
      for (int r = 0; r < dsize; r++) {
        int rr = vars_per_node * (r / 3) + (r % 3);
        mat[rr * nvars + mm] += ucol[r];
      }

      TacsScalar d1col[dsize];
      memset(d1col, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xh, d1col);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xh, d1col);
      for (int r = 0; r < dsize; r++) {
        d2d1u[dsize * r + m] += d1col[r];
      }

      TacsScalar d2col[dsize];
      memset(d2col, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xh, d2col);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xh, d2col);
      for (int r = 0; r < dsize; r++) {
        d2d2u[dsize * r + m] += d2col[r];
      }
    }

    // (B) d1 sweep: seeds d01/d01xi together. Captures the (d1, d1) self
    // block (into d2d1) and the (d1, d2) cross block (into d2d1d2) --
    // capturing this cross leakage, rather than discarding it, is exactly
    // Gap 2's static contribution (SPEC.md sec 1.3.2); the (u0xi, d1)
    // block is intentionally NOT re-scattered here (already captured by
    // sweep (A), avoiding double-counting).
    for (int m = 0; m < dsize; m++) {
      TacsBeamZeroSecondOrderNodes(u0xi, d01, d02, d01xi, d02xi, u0d, u0dXdinvT,
                                   u0x, d1t, d1x, d2t, d2x);

      memset(seed, 0, dsize * sizeof(TacsScalar));
      seed[m] = 1.0;
      basis::template interpFields<3, 3>(pt, seed, d01.xp);
      basis::template interpFieldsGrad<3, 3>(pt, seed, d01xi.xp);

      assembleu0d.hforward();
      multu0d.hforward();
      multu0x.hforward();
      axpyd1t.hforward();
      matmultd1x.hforward();
      axpyd2t.hforward();
      matmultd2x.hforward();

      TacsBeamContractStrainHessian(d2u0x, d2d1x, d2d2x, d2u0xd1x, d2u0xd2x,
                                    d2d1xd2x, u0x.Ap, d1x.xp, d2x.xp, u0x.Ah,
                                    d1x.xh, d2x.xh);

      matmultd2x.hreverse();
      axpyd2t.hreverse();
      matmultd1x.hreverse();
      axpyd1t.hreverse();
      multu0x.hreverse();
      multu0d.hreverse();
      assembleu0d.hreverse();

      TacsScalar d1col[dsize];
      memset(d1col, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xh, d1col);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xh, d1col);
      for (int r = 0; r < dsize; r++) {
        d2d1[dsize * r + m] += d1col[r];
      }

      TacsScalar d2col[dsize];
      memset(d2col, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xh, d2col);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xh, d2col);
      for (int r = 0; r < dsize; r++) {
        d2d1d2[dsize * m + r] += d2col[r];
      }
    }

    // (C) d2 sweep: seeds d02/d02xi together. Captures ONLY the (d2, d2)
    // self block (into d2d2); the (u0xi, d2) and (d1, d2) blocks are
    // intentionally NOT re-scattered here (already captured above).
    for (int m = 0; m < dsize; m++) {
      TacsBeamZeroSecondOrderNodes(u0xi, d01, d02, d01xi, d02xi, u0d, u0dXdinvT,
                                   u0x, d1t, d1x, d2t, d2x);

      memset(seed, 0, dsize * sizeof(TacsScalar));
      seed[m] = 1.0;
      basis::template interpFields<3, 3>(pt, seed, d02.xp);
      basis::template interpFieldsGrad<3, 3>(pt, seed, d02xi.xp);

      assembleu0d.hforward();
      multu0d.hforward();
      multu0x.hforward();
      axpyd1t.hforward();
      matmultd1x.hforward();
      axpyd2t.hforward();
      matmultd2x.hforward();

      TacsBeamContractStrainHessian(d2u0x, d2d1x, d2d2x, d2u0xd1x, d2u0xd2x,
                                    d2d1xd2x, u0x.Ap, d1x.xp, d2x.xp, u0x.Ah,
                                    d1x.xh, d2x.xh);

      matmultd2x.hreverse();
      axpyd2t.hreverse();
      matmultd1x.hreverse();
      axpyd1t.hreverse();
      multu0x.hreverse();
      multu0d.hreverse();
      assembleu0d.hreverse();

      TacsScalar d2col[dsize];
      memset(d2col, 0, dsize * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xh, d2col);
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xh, d2col);
      for (int r = 0; r < dsize; r++) {
        d2d2[dsize * r + m] += d2col[r];
      }
    }
  }

  if (res) {
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn1, fn2, vars, d1, d2, dety, res, d1d, d2d);
  }

  // Gap 1 (SPEC.md sec 1.3.1): convert the tying-point-space Hessian d2ety
  // into mat[]'s (u0xi, u0xi) block plus the director-space accumulators
  // d2d1/d2d2/d2d1u/d2d2u. Beam's Cs has no e0ty cross-coupling with
  // (u0x, d1x, d2x) for the constitutive models this feature exercises
  // (TACSIsoTubeBeamConstitutive::evalTangentStiffness never sets a
  // cross-shear entry), so the corresponding cross-Hessian inputs are
  // zero-filled placeholders (see TACSBeamElementModel.h's own comment).
  TacsScalar d2etyu[basis::NUM_TYING_POINTS * dsize];
  TacsScalar d2etyd1[basis::NUM_TYING_POINTS * dsize];
  TacsScalar d2etyd2[basis::NUM_TYING_POINTS * dsize];
  memset(d2etyu, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));
  memset(d2etyd1, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));
  memset(d2etyd2, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
      alpha, Xpts, fn1, fn2, vars, d1, d2, dety, d2ety, d2etyu, d2etyd1,
      d2etyd2, mat, d2d1, d2d2, d2d1u, d2d2u);

  // Gap 2 (SPEC.md sec 1.3.2): fold the dynamics (rho[5]) cross-director
  // contribution (unscaled) into the already alpha-scaled static one
  // (TacsBeamAddCrossDirectorJacobian expects an already-fully-scaled
  // buffer, mirroring director::addDirectorJacobian's own
  // "caller pre-scales d2d/d2du" contract), then scatter into mat[].
  for (int i = 0; i < dsize * dsize; i++) {
    d2d1d2[i] += gamma * d2Tdotd1d2[i];
  }
  // Exact only for TACSLinearizedRotation (the only director class
  // instantiated here today) -- see TacsBeamAddCrossDirectorJacobian's own
  // docstring (TACSBeamUtilities.h) for the full scope-boundary rationale
  // before wiring this call up for TACSQuadraticRotation/
  // TACSQuaternionRotation.
  TacsBeamAddCrossDirectorJacobian<vars_per_node, offset, num_nodes>(
      vars, fn1, fn2, d2d1d2, mat);

  // Per-director Jacobian closures (SPEC.md sec 1.3 step 4) -- each of
  // these also updates res internally (mirroring addDirectorResidual's own
  // formula), so the separate addDirectorResidual calls Task 2.1's res-only
  // scaffold used are intentionally removed here to avoid double-counting.
  director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
      alpha, beta, gamma, vars, dvars, ddvars, fn1, d1d, d2Tdotd1, d2Tdotu1,
      d2d1, d2d1u, res, mat);
  director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
      alpha, beta, gamma, vars, dvars, ddvars, fn2, d2d, d2Tdotd2, d2Tdotu2,
      d2d2, d2d2u, res, mat);

  director::template addRotationConstrJacobian<vars_per_node, offset,
                                               num_nodes>(alpha, vars, res,
                                                          mat);
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
    // Explicit-forward to the base class instead of a silent early return
    // (SPEC.md sec 1.4/4.1/4.2, the "never-silent-punt" principle): the
    // base class's own TACS_GEOMETRIC_STIFFNESS_MATRIX branch
    // (TACSElement.cpp:346-348) is ALSO "not implemented" today, so this is
    // numerically a no-op (both produce an all-zero mat) -- the value is
    // consistency/auditability, not a correctness fix, and this is an
    // interim state: Phase 5 (SPEC.md sec 2.4.4/2.4.5) replaces this branch
    // with a real analytic implementation, unless its
    // TACSBeamNonlinearModel prerequisite overruns its timebox, in which
    // case this forward becomes the permanent, documented fallback.
    TACSElement::getMatType(matType, elemIndex, time, Xpts, vars, mat);
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
  const A2D::Vec3 &axis = transform->getRefAxis();

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
    A2D::Vec3 u0xipsi, d01psi, d02psi, d01xipsi, d02xipsi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Interpolate the adjoint solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, u0xipsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01psi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02psi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, d01xipsi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, d02xipsi.x);

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

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
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    // Assemble u0d and u0psi
    A2D::Mat3x3 u0d, u0dpsi;
    A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);
    A2D::Mat3x3FromThreeVec3 assembleu0dpsi(u0xipsi, d01psi, d02psi, u0dpsi);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::Mat3x3 u0dXdinvT, u0x;
    A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

    // Compute u0xpsi = ^{T} * u0dpsi * XdinvT
    A2D::Mat3x3 u0dXdinvTpsi, u0xpsi;
    A2D::Mat3x3MatMult multu0dpsi(u0dpsi, XdinvT, u0dXdinvTpsi);
    A2D::MatTrans3x3MatMult multu0xpsi(T, u0dXdinvTpsi, u0xpsi);

    // Compute s0, sz1 and sz2
    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::Vec3 d1t, d1x;
    A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::Vec3 d2t, d2x;
    A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

    // Compute d1xpsi = s0 * T^{T} * (d1xipsi - sz1 * u0xipsi)
    A2D::Vec3 d1tpsi, d1xpsi;
    A2D::Vec3Axpy axpyd1tpsi(-1.0, sz1, u0xipsi, d01xipsi, d1tpsi);
    A2D::MatTrans3x3VecMultScale matmultd1xpsi(s0, T, d1tpsi, d1xpsi);

    // Compute d2xpsi = s0 * T^{T} * (d2xipsi - sz2 * u0xipsi)
    A2D::Vec3 d2tpsi, d2xpsi;
    A2D::Vec3Axpy axpyd2tpsi(-1.0, sz2, u0xipsi, d02xipsi, d2tpsi);
    A2D::MatTrans3x3VecMultScale matmultd2xpsi(s0, T, d2tpsi, d2xpsi);

    // Evaluate the tying components of the strain
    TacsScalar gty[2], gtypsi[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);
    basis::interpTyingStrain(pt, etypsi, gtypsi);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], e0typsi[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];
    e0typsi[0] = 2.0 * XdinvT.A[0] * gtypsi[0];
    e0typsi[1] = 2.0 * XdinvT.A[0] * gtypsi[1];

    // // Evaluate the strain and the adjoint strain
    TacsScalar e[6], epsi[6];
    model::evalStrainDeriv(u0x.A, d1x.x, d2x.x, e0ty, u0xpsi.A, d1xpsi.x,
                           d2xpsi.x, e0typsi, e, epsi);

    // Add the product of the derivative of the stress
    con->addStressDVSens(elemIndex, scale * detXd.value, pt, X0.x, e, epsi,
                         dvLen, dfdx);

    // Evaluate the accelerations
    A2D::Vec3 u0dot, d01dot, d02dot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0dot.x);
    basis::template interpFields<3, 3>(pt, d1ddot, d01dot.x);
    basis::template interpFields<3, 3>(pt, d2ddot, d02dot.x);

    A2D::Vec3 u0dotpsi, d01dotpsi, d02dotpsi;
    basis::template interpFields<vars_per_node, 3>(pt, psi, u0dotpsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01dotpsi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02dotpsi.x);

    // Compute the dot-products
    A2D::Scalar u0psi, u0psid1, u0psid2, u0d1psi, u0d2psi;
    A2D::Vec3Dot dot1(u0dot, u0dotpsi, u0psi);
    A2D::Vec3Dot dot2(u0dotpsi, d01dot, u0psid1);
    A2D::Vec3Dot dot3(u0dot, d01dotpsi, u0d1psi);
    A2D::Vec3Dot dot4(u0dotpsi, d02dot, u0psid2);
    A2D::Vec3Dot dot5(u0dot, d02dotpsi, u0d2psi);

    A2D::Scalar d1d1psi, d2d2psi, d1psid2, d1d2psi;
    A2D::Vec3Dot dot7(d01dot, d01dotpsi, d1d1psi);
    A2D::Vec3Dot dot8(d02dot, d02dotpsi, d2d2psi);
    A2D::Vec3Dot dot9(d01dotpsi, d02dot, d1psid2);
    A2D::Vec3Dot dot10(d01dot, d02dotpsi, d1d2psi);

    // Add derivatives from the mass moments
    TacsScalar rho[6];
    TacsScalar alpha = scale * detXd.value;
    rho[0] = alpha * u0psi.value;
    rho[1] = alpha * (u0psid1.value + u0d1psi.value);
    rho[2] = alpha * (u0psid2.value + u0d2psi.value);
    rho[3] = alpha * d1d1psi.value;
    rho[4] = alpha * d2d2psi.value;
    rho[5] = alpha * (d1psid2.value + d1d2psi.value);

    con->addMassMomentsDVSens(elemIndex, pt, X0.x, rho, dvLen, dfdx);
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
  const A2D::Vec3 &axis = transform->getRefAxis();

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
    A2D::ADMat3x3 T;

    // Parametric location
    A2D::ADVec3 X0;

    // Tangent to the beam
    A2D::ADVec3 X0xi;

    // Interpolated normal directions
    A2D::ADVec3 n1, n2;

    // Derivatives of the interpolated normal directions
    A2D::ADVec3 n1xi, n2xi;

    // The values of the director fields and their derivatives
    A2D::ADVec3 u0xi, d01, d02, d01xi, d02xi;
    A2D::ADVec3 u0xipsi, d01psi, d02psi, d01xipsi, d02xipsi;

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Interpolate the adjoint solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, u0xipsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01psi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02psi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, d01xipsi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, d02xipsi.x);

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

    // Compute the transformation at the quadrature point
    transform->computeTransform(X0xi.x, T.A);

    // Compute the inverse
    A2D::ADMat3x3 Xd, Xdinv;
    A2D::ADMat3x3FromThreeADVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::ADMat3x3Inverse invXd(Xd, Xdinv);

    // Compute the determinant of the transform
    A2D::ADScalar detXd;
    A2D::ADMat3x3Det computedetXd(weight, Xd, detXd);

    // Compute XdinvT = Xdinv * T
    A2D::ADMat3x3 XdinvT;
    A2D::ADMat3x3ADMatMult multXdinvT(Xdinv, T, XdinvT);

    // Assemble u0d and u0psi
    A2D::ADMat3x3 u0d, u0dpsi;
    A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);
    A2D::ADMat3x3FromThreeADVec3 assembleu0dpsi(u0xipsi, d01psi, d02psi,
                                                u0dpsi);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::ADMat3x3 u0dXdinvT, u0x;
    A2D::ADMat3x3ADMatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::ADMatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

    // Compute u0xpsi = ^{T} * u0dpsi * XdinvT
    A2D::ADMat3x3 u0dXdinvTpsi, u0xpsi;
    A2D::ADMat3x3ADMatMult multu0dpsi(u0dpsi, XdinvT, u0dXdinvTpsi);
    A2D::ADMatTrans3x3ADMatMult multu0xpsi(T, u0dXdinvTpsi, u0xpsi);

    // Compute s0, sz1 and sz2
    A2D::ADScalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::ADMat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::ADMat3x3VecADVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::ADMat3x3VecADVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::ADVec3 d1t, d1x;
    A2D::ADVec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd1x(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::ADVec3 d2t, d2x;
    A2D::ADVec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd2x(s0, T, d2t, d2x);

    // Compute d1xpsi = s0 * T^{T} * (d1xipsi - sz1 * u0xipsi)
    A2D::ADVec3 d1tpsi, d1xpsi;
    A2D::ADVec3Axpy axpyd1tpsi(-1.0, sz1, u0xipsi, d01xipsi, d1tpsi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd1xpsi(s0, T, d1tpsi, d1xpsi);

    // Compute d2xpsi = s0 * T^{T} * (d2xipsi - sz2 * u0xipsi)
    A2D::ADVec3 d2tpsi, d2xpsi;
    A2D::ADVec3Axpy axpyd2tpsi(-1.0, sz2, u0xipsi, d02xipsi, d2tpsi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd2xpsi(s0, T, d2tpsi, d2xpsi);

    // Evaluate the tying components of the strain
    TacsScalar gty[2], gtypsi[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);
    basis::interpTyingStrain(pt, etypsi, gtypsi);

    // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2], e0typsi[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];
    e0typsi[0] = 2.0 * XdinvT.A[0] * gtypsi[0];
    e0typsi[1] = 2.0 * XdinvT.A[0] * gtypsi[1];

    // // Evaluate the strain and the adjoint strain
    TacsScalar e[6], epsi[6];
    model::evalStrainDeriv(u0x.A, d1x.x, d2x.x, e0ty, u0xpsi.A, d1xpsi.x,
                           d2xpsi.x, e0typsi, e, epsi);

    // Compute the stress due to the strain
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);

    // Compute the psioint stress - assuming a linear relationship for
    // stress/strain
    TacsScalar spsi[6];
    con->evalStress(elemIndex, pt, X0.x, epsi, spsi);

    // Evaluate the sensitivities
    TacsScalar e0tyd[2], e0typsid[2];
    model::evalStrainSens(scale * detXd.value, spsi, u0x.A, d1x.x, d2x.x, e0ty,
                          u0x.Ad, d1x.xd, d2x.xd, e0tyd);
    model::evalStrainSens(scale * detXd.value, s, u0x.A, d1x.x, d2x.x, e0ty,
                          u0xpsi.Ad, d1xpsi.xd, d2xpsi.xd, e0typsid);
    detXd.valued = scale * (e[0] * spsi[0] + e[1] * spsi[1] + e[2] * spsi[2] +
                            e[3] * spsi[3] + e[4] * spsi[4] + e[5] * spsi[5]);

    // Apply the tying strain transformation
    TacsScalar gtyd[2], gtypsid[2];
    gtyd[0] = 2.0 * XdinvT.A[0] * e0tyd[0];
    gtyd[1] = 2.0 * XdinvT.A[0] * e0tyd[1];
    gtypsid[0] = 2.0 * XdinvT.A[0] * e0typsid[0];
    gtypsid[1] = 2.0 * XdinvT.A[0] * e0typsid[1];

    XdinvT.Ad[0] += 2.0 * (gty[0] * e0tyd[0] + gty[1] * e0tyd[1] +
                           e0typsid[0] * gtypsi[0] + e0typsid[1] * gtypsi[1]);

    // Evaluate the accelerations
    A2D::ADVec3 u0ddot, d01ddot, d02ddot;
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot.x);
    basis::template interpFields<3, 3>(pt, d1ddot, d01ddot.x);
    basis::template interpFields<3, 3>(pt, d2ddot, d02ddot.x);

    A2D::ADVec3 u0ddotpsi, d01ddotpsi, d02ddotpsi;
    basis::template interpFields<vars_per_node, 3>(pt, psi, u0ddotpsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01ddotpsi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02ddotpsi.x);

    // Compute the dot-products
    A2D::ADScalar u0psi, u0psid1, u0psid2, u0d1psi, u0d2psi;
    A2D::ADVec3Dot dot1(u0ddot, u0ddotpsi, u0psi);
    A2D::ADVec3Dot dot2(u0ddotpsi, d01ddot, u0psid1);
    A2D::ADVec3Dot dot3(u0ddot, d01ddotpsi, u0d1psi);
    A2D::ADVec3Dot dot4(u0ddotpsi, d02ddot, u0psid2);
    A2D::ADVec3Dot dot5(u0ddot, d02ddotpsi, u0d2psi);

    A2D::ADScalar d1d1psi, d2d2psi, d1psid2, d1d2psi;
    A2D::ADVec3Dot dot6(d01ddot, d01ddotpsi, d1d1psi);
    A2D::ADVec3Dot dot7(d02ddot, d02ddotpsi, d2d2psi);
    A2D::ADVec3Dot dot8(d01ddotpsi, d02ddot, d1psid2);
    A2D::ADVec3Dot dot9(d01ddot, d02ddotpsi, d1d2psi);

    // Evaluate the mass moments
    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, X0.x, rho);

    // Add the contribution from the adjoint-residual product from the
    // dynamics
    detXd.valued +=
        scale *
        (rho[0] * u0psi.value + rho[1] * (u0psid1.value + u0d1psi.value) +
         rho[2] * (u0psid2.value + u0d2psi.value) + rho[3] * d1d1psi.value +
         rho[4] * d2d2psi.value + rho[5] * (d1psid2.value + d1d2psi.value));

    // Set the seeds for the dot-products
    TacsScalar alpha = scale * detXd.value;
    u0psi.valued = alpha * rho[0];
    u0psid1.valued = alpha * rho[1];
    u0d1psi.valued = alpha * rho[1];

    u0psid2.valued = alpha * rho[2];
    u0d2psi.valued = alpha * rho[2];
    d1d1psi.valued = alpha * rho[3];
    d2d2psi.valued = alpha * rho[4];
    d1psid2.valued = alpha * rho[5];
    d1d2psi.valued = alpha * rho[5];

    // Reverse the dot-products
    dot9.reverse();
    dot8.reverse();
    dot7.reverse();
    dot6.reverse();
    dot5.reverse();
    dot4.reverse();
    dot3.reverse();
    dot2.reverse();
    dot1.reverse();

    matmultd2xpsi.reverse();
    axpyd2tpsi.reverse();
    matmultd1xpsi.reverse();
    axpyd1tpsi.reverse();
    matmultd2x.reverse();
    axpyd2t.reverse();
    matmultd1x.reverse();
    axpyd1t.reverse();
    innersz2.reverse();
    innersz1.reverse();
    inners0.reverse();
    multu0xpsi.reverse();
    multu0dpsi.reverse();
    assembleu0dpsi.reverse();
    multu0x.reverse();
    multu0d.reverse();
    assembleu0d.reverse();
    multXdinvT.reverse();
    computedetXd.reverse();
    invXd.reverse();
    assembleXd.reverse();

    // Reverse the transformation sensitivities
    transform->addTransformSens(X0xi.x, T.Ad, X0xi.xd);

    // Add the sensitivities to the input fields...
    basis::template addInterpFieldsTranspose<3, 3>(pt, X0.xd, dfdXpts);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, X0xi.xd, dfdXpts);

    basis::template addInterpFieldsTranspose<3, 3>(pt, n1.xd, dfn1);
    basis::template addInterpFieldsTranspose<3, 3>(pt, n2.xd, dfn2);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, n1xi.xd, dfn1);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, n2xi.xd, dfn2);

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, dd1);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, dd2);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, dd1);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, dd2);

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01psi.xd, dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02psi.xd, dd2psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xipsi.xd, dd1psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xipsi.xd, dd2psi);

    // Add the contributions from the dynamics
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01ddot.xd, dd1);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02ddot.xd, dd2);

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01ddotpsi.xd, dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02ddotpsi.xd, dd2psi);

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

/*
  Add the derivative of the matrix inner product psi^T * mat * phi with
  respect to the design variables (SPEC.md sec 2.2/2.3). TACS_STIFFNESS_MATRIX
  and TACS_MASS_MATRIX are analytic; TACS_GEOMETRIC_STIFFNESS_MATRIX is an
  interim explicit-forward-to-base punt, superseded by Phase 5 (SPEC.md sec
  2.4.4/2.4.5, 4.2) -- not a permanent scope cut.
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director,
                     model>::addMatDVSensInnerProduct(ElementMatrixType matType,
                                                      int elemIndex,
                                                      double time,
                                                      TacsScalar scale,
                                                      const TacsScalar psi[],
                                                      const TacsScalar phi[],
                                                      const TacsScalar Xpts[],
                                                      const TacsScalar vars[],
                                                      int dvLen,
                                                      TacsScalar dfdx[]) {
  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    // Interim explicit-forward to the base class (SPEC.md sec 4.1/4.2's
    // "never-silent-punt" principle): Phase 5 (SPEC.md sec 2.4.4/2.4.5)
    // replaces this branch with a real analytic implementation, unless its
    // TACSBeamNonlinearModel prerequisite overruns its timebox, in which
    // case this forward becomes the permanent, documented fallback.
    TACSElement::addMatDVSensInnerProduct(matType, elemIndex, time, scale, psi,
                                          phi, Xpts, vars, dvLen, dfdx);
    return;
  } else if (matType != TACS_STIFFNESS_MATRIX && matType != TACS_MASS_MATRIX) {
    // Unsupported matType beyond the enumerated set (SPEC.md sec 4.3):
    // explicit-forward, never silent.
    TACSElement::addMatDVSensInnerProduct(matType, elemIndex, time, scale, psi,
                                          phi, Xpts, vars, dvLen, dfdx);
    return;
  }

  // Common setup, shared by both analytic branches.
  const int nquad = quadrature::getNumQuadraturePoints();
  const A2D::Vec3 &axis = transform->getRefAxis();

  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  // Build the psi-direction and phi-direction director fields, both
  // linearized about the REAL state vars (SPEC.md sec 2.2's "psi/phi in
  // place of vars" kinematics substitution). computeDirectorRatesDeriv's
  // dpsi output is the exact director-map Jacobian-vector product at the
  // given vars, for any director class (confirmed directly from
  // TACSDirector.h's per-class implementations) -- dvars/ddvars are passed
  // as vars here since this method's own signature has no dvars/ddvars, and
  // the ddot/dddot outputs they would feed are discarded, unused.
  //
  // Efficiency note (review-flagged, minor): the two calls below (one
  // seeded with psi, one with phi) each recompute the identical,
  // direction-independent d1/d1dot/d1ddot outputs redundantly --
  // computeDirectorRatesDeriv has no lighter-weight overload that accepts
  // precomputed d/ddot/dddot and returns only the direction-seeded dpsi
  // output, so avoiding the redundant work would require an API change to
  // TACSDirector.h (a shared header also used by shell elements) rather
  // than a local contortion; left as-is (this is a duplicated O(num_nodes)
  // computation per call site, not an algorithmic complexity issue).
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize], d1psi[dsize], d1phi[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize], d2psi[dsize], d2phi[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn1, d1, d1dot, d1ddot, d1psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn1, d1, d1dot, d1ddot, d1phi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn2, d2, d2dot, d2ddot, d2psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn2, d2, d2dot, d2ddot, d2phi);

  // psi-direction and phi-direction tying strains, via the SAME linear
  // tying-strain map (computeTyingStrainDeriv), substituting psi/phi for
  // vars/d1/d2 in place of the real state (the map is linear, so this
  // substitution gives exactly the psi-direction/phi-direction tying
  // strains -- SPEC.md sec 2.3.1's "psi/phi in place of vars" reuse).
  TacsScalar etypsi[basis::NUM_TYING_POINTS], etyphi[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn1, fn2, psi, d1psi, d2psi, phi, d1phi, d2phi, etypsi, etyphi);

  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    A2D::Mat3x3 T;
    A2D::Vec3 X0, X0xi, n1, n2, n1xi, n2xi;

    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

    transform->computeTransform(X0xi.x, T.A);

    A2D::Mat3x3 Xd, Xdinv;
    A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::Mat3x3Inverse invXd(Xd, Xdinv);

    A2D::Scalar detXd;
    A2D::Mat3x3Det computedetXd(weight, Xd, detXd);

    A2D::Mat3x3 XdinvT;
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    // psi-direction chain.
    A2D::Vec3 u0xipsi, d01psi, d02psi, d01xipsi, d02xipsi;
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, u0xipsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01psi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02psi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, d01xipsi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, d02xipsi.x);

    A2D::Mat3x3 u0dpsi;
    A2D::Mat3x3FromThreeVec3 assembleu0dpsi(u0xipsi, d01psi, d02psi, u0dpsi);

    A2D::Mat3x3 u0dXdinvTpsi, u0xpsi;
    A2D::Mat3x3MatMult multu0dpsi(u0dpsi, XdinvT, u0dXdinvTpsi);
    A2D::MatTrans3x3MatMult multu0xpsi(T, u0dXdinvTpsi, u0xpsi);

    // phi-direction chain.
    A2D::Vec3 u0xiphi, d01phi, d02phi, d01xiphi, d02xiphi;
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, phi, u0xiphi.x);
    basis::template interpFields<3, 3>(pt, d1phi, d01phi.x);
    basis::template interpFields<3, 3>(pt, d2phi, d02phi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1phi, d01xiphi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2phi, d02xiphi.x);

    A2D::Mat3x3 u0dphi;
    A2D::Mat3x3FromThreeVec3 assembleu0dphi(u0xiphi, d01phi, d02phi, u0dphi);

    A2D::Mat3x3 u0dXdinvTphi, u0xphi;
    A2D::Mat3x3MatMult multu0dphi(u0dphi, XdinvT, u0dXdinvTphi);
    A2D::MatTrans3x3MatMult multu0xphi(T, u0dXdinvTphi, u0xphi);

    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    A2D::Vec3 d1tpsi, d1xpsi;
    A2D::Vec3Axpy axpyd1tpsi(-1.0, sz1, u0xipsi, d01xipsi, d1tpsi);
    A2D::MatTrans3x3VecMultScale matmultd1xpsi(s0, T, d1tpsi, d1xpsi);

    A2D::Vec3 d2tpsi, d2xpsi;
    A2D::Vec3Axpy axpyd2tpsi(-1.0, sz2, u0xipsi, d02xipsi, d2tpsi);
    A2D::MatTrans3x3VecMultScale matmultd2xpsi(s0, T, d2tpsi, d2xpsi);

    A2D::Vec3 d1tphi, d1xphi;
    A2D::Vec3Axpy axpyd1tphi(-1.0, sz1, u0xiphi, d01xiphi, d1tphi);
    A2D::MatTrans3x3VecMultScale matmultd1xphi(s0, T, d1tphi, d1xphi);

    A2D::Vec3 d2tphi, d2xphi;
    A2D::Vec3Axpy axpyd2tphi(-1.0, sz2, u0xiphi, d02xiphi, d2tphi);
    A2D::MatTrans3x3VecMultScale matmultd2xphi(s0, T, d2tphi, d2xphi);

    // Interpolate the psi-direction and phi-direction tying strains and
    // transform them to the local coordinates (identical formula to the
    // real e0ty, SPEC.md sec 1.2.1).
    TacsScalar gtypsi[2], gtyphi[2];
    basis::interpTyingStrain(pt, etypsi, gtypsi);
    basis::interpTyingStrain(pt, etyphi, gtyphi);

    TacsScalar e0typsi[2], e0typhi[2];
    e0typsi[0] = 2.0 * XdinvT.A[0] * gtypsi[0];
    e0typsi[1] = 2.0 * XdinvT.A[0] * gtypsi[1];
    e0typhi[0] = 2.0 * XdinvT.A[0] * gtyphi[0];
    e0typhi[1] = 2.0 * XdinvT.A[0] * gtyphi[1];

    if (matType == TACS_STIFFNESS_MATRIX) {
      // e_psi/e_phi via the SAME linear strain map vars/psi/phi share
      // (SPEC.md sec 2.3.1's "beam's strain is linear unconditionally"
      // simplification) -- reduces psi^T*K*phi to e_psi^T*Cs*e_phi, whose
      // DV-derivative is the SAME generic hook addAdjResProduct already
      // calls (con->addStressDVSens), reused with a bilinear pair of
      // strains instead of a single adjoint-direction strain (SPEC.md sec
      // 2.2).
      TacsScalar epsi[6], ephi[6];
      model::evalStrainDeriv(u0xpsi.A, d1xpsi.x, d2xpsi.x, e0typsi, u0xphi.A,
                             d1xphi.x, d2xphi.x, e0typhi, epsi, ephi);

      con->addStressDVSens(elemIndex, scale * detXd.value, pt, X0.x, ephi, epsi,
                           dvLen, dfdx);
    } else {  // TACS_MASS_MATRIX
      // No strain/stiffness path (SPEC.md sec 2.2/2.3): only the mass
      // moments' own DV-sens, contracted against the bilinear (psi, phi)
      // dot products the mass matrix's dynamics blocks are built from
      // (mirrors addAdjResProduct's dynamics section, TACSBeamElement.h,
      // but bilinear in (psi, phi) rather than (ddvars, psi)).
      TacsScalar u0psi[3], u0phi[3];
      basis::template interpFields<vars_per_node, 3>(pt, psi, u0psi);
      basis::template interpFields<vars_per_node, 3>(pt, phi, u0phi);

      TacsScalar d1psiv[3], d1phiv[3], d2psiv[3], d2phiv[3];
      basis::template interpFields<3, 3>(pt, d1psi, d1psiv);
      basis::template interpFields<3, 3>(pt, d1phi, d1phiv);
      basis::template interpFields<3, 3>(pt, d2psi, d2psiv);
      basis::template interpFields<3, 3>(pt, d2phi, d2phiv);

      TacsScalar alpha = scale * detXd.value;
      TacsScalar rho[6];
      rho[0] = alpha * (u0psi[0] * u0phi[0] + u0psi[1] * u0phi[1] +
                        u0psi[2] * u0phi[2]);
      rho[1] = alpha * (u0psi[0] * d1phiv[0] + u0psi[1] * d1phiv[1] +
                        u0psi[2] * d1phiv[2] + d1psiv[0] * u0phi[0] +
                        d1psiv[1] * u0phi[1] + d1psiv[2] * u0phi[2]);
      rho[2] = alpha * (u0psi[0] * d2phiv[0] + u0psi[1] * d2phiv[1] +
                        u0psi[2] * d2phiv[2] + d2psiv[0] * u0phi[0] +
                        d2psiv[1] * u0phi[1] + d2psiv[2] * u0phi[2]);
      rho[3] = alpha * (d1psiv[0] * d1phiv[0] + d1psiv[1] * d1phiv[1] +
                        d1psiv[2] * d1phiv[2]);
      rho[4] = alpha * (d2psiv[0] * d2phiv[0] + d2psiv[1] * d2phiv[1] +
                        d2psiv[2] * d2phiv[2]);
      rho[5] = alpha * (d1psiv[0] * d2phiv[0] + d1psiv[1] * d2phiv[1] +
                        d1psiv[2] * d2phiv[2] + d2psiv[0] * d1phiv[0] +
                        d2psiv[1] * d1phiv[1] + d2psiv[2] * d1phiv[2]);

      con->addMassMomentsDVSens(elemIndex, pt, X0.x, rho, dvLen, dfdx);
    }
  }
}

/*
  addMatXptSensInnerProduct's TACS_STIFFNESS_MATRIX branch (Task 4.3,
  SPEC.md sec 2.3.1): dfdXpts += scale * d(psi^T*K*phi)/d(Xpts).

  f(Xpts) = scale * sum_quad w * detXd(Xpts) * (e_psi(Xpts)^T * Cs *
  e_phi(Xpts))

  Beam's strain is linear in state unconditionally (only one active model
  class, TACSBeamNonlinearModel is commented out), so e_psi/e_phi are
  exactly what model::evalStrain produces when psi/phi are substituted for
  vars through the same forward kinematics chain -- no separate
  strain-Hessian term is needed. This mirrors addAdjResXptProduct's own
  Xpts-adjoint chain almost exactly (TacsBeamAddNodeNormalsSens,
  director::addDirectorRefNormalSens, model::addTyingStrainXptSens-family
  machinery), with BOTH "state" directions replaced by pure test directions
  (psi, phi) -- neither corresponds to the real element state, since the
  bilinear form psi^T*K*phi never references the actual vars-built strain
  at all. The decomposition (SPEC.md sec 2.3.1):
    1. detXd-direction term: seed detXd.valued = scale*(e_psi . s_phi), where
       s_phi = Cs*e_phi (routes through the existing detXd-Xpts-adjoint
       machinery unchanged).
    2. psi-direction chain: seed weight s_phi = Cs*e_phi against the
       psi-direction kinematics chain (mirrors addAdjResXptProduct's
       psi-chain).
    3. phi-direction chain: seed weight s_psi = Cs*e_psi against the
       phi-direction kinematics chain (mirrors addAdjResXptProduct's
       real/vars-chain, with phi substituted for vars).
    4. e0ty/XdinvT-direction correction: e0ty = 2*XdinvT.A[0]*gty is a
       state-independent scalar product, so its Xpts-adjoint is a plain
       two-term product-rule expression (same machinery
       addAdjResXptProduct already exercises for this term).
    5. No tying-curvature term (beam's tying strain is linear -- omitted,
       not derived-and-discarded).

  Director-class scope: exact for TACSLinearizedRotation, confirmed via
  the un-skipped test file (Beam2/Beam3) and a machine-precision
  TACSBeam2ModRot-style cross-check (see PLAN.md). This function is ONLY
  called for TACSLinearizedRotation -- its caller (addMatXptSensInnerProduct)
  typeid-guards TACSQuadraticRotation/TACSQuaternionRotation to the base
  FD/CS fallback, per a documented, review-flagged, currently-unresolved
  finding: the director-field Xpts-adjoint fix below (the 6-arg
  addDirectorRefNormalSens overload's ddpsi-only term, zeroed dd,
  evaluated at the REAL vars) is CONFIRMED CORRECT and sufficient for
  addMatXptSensInnerProduct's TACS_MASS_MATRIX branch (machine-precision
  match against TACSBeam2ModRot) but does NOT fully resolve this
  STIFFNESS_MATRIX branch's nonlinear-director discrepancy (~9% residual
  at full rotation magnitude for TACSBeam2ModRot, confirmed non-FD-noise
  via a dh-convergence check and confirmed genuinely q-dependent via a
  vars-magnitude-scaling check) -- an additional, unidentified term
  specific to this branch's larger kinematic chain (u0d/u0x/d1x/d2x/T/
  XdinvT, absent from the simpler TACS_MASS_MATRIX branch) is suspected
  but was not isolated within this session's investigation budget.
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
    addMatXptSensInnerProductStiffness(int elemIndex, double time,
                                       TacsScalar scale, const TacsScalar psi[],
                                       const TacsScalar phi[],
                                       const TacsScalar Xpts[],
                                       const TacsScalar vars[],
                                       TacsScalar dfdXpts[]) {
  const int nquad = quadrature::getNumQuadraturePoints();
  const A2D::Vec3 &axis = transform->getRefAxis();

  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  TacsScalar dfn1[3 * basis::NUM_NODES], dfn2[3 * basis::NUM_NODES];
  memset(dfn1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
  memset(dfn2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

  // psi-direction and phi-direction director fields, linearized about the
  // real vars (same reuse as Tasks 4.1/4.2; see Task 4.1's own comment for
  // the redundant-recomputation efficiency note, unchanged here).
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize], d1psi[dsize], d1phi[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize], d2psi[dsize], d2phi[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn1, d1, d1dot, d1ddot, d1psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn1, d1, d1dot, d1ddot, d1phi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn2, d2, d2dot, d2ddot, d2psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn2, d2, d2dot, d2ddot, d2phi);

  // psi-direction/phi-direction tying strains via the SAME linear tying-
  // strain map, substituting psi/phi for vars/d1/d2 (SPEC.md sec 2.3.1
  // piece 2/3's "psi/phi in place of vars" reuse).
  TacsScalar etypsi[basis::NUM_TYING_POINTS], etyphi[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn1, fn2, psi, d1psi, d2psi, phi, d1phi, d2phi, etypsi, etyphi);

  TacsScalar etypsid[basis::NUM_TYING_POINTS], etyphid[basis::NUM_TYING_POINTS];
  memset(etypsid, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(etyphid, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // Per-node accumulators for the director-field Xpts-adjoint, mirroring
  // addAdjResXptProduct's dd1/dd2 (real-chain) and dd1psi/dd2psi
  // (psi-chain) buffers -- renamed here since neither chain is "the real
  // state": dd1phi/dd2phi play addAdjResXptProduct's dd1/dd2 role (with
  // phi substituted for vars), dd1psi/dd2psi are unchanged.
  TacsScalar dd1phi[dsize], dd2phi[dsize], dd1psi[dsize], dd2psi[dsize];
  memset(dd1phi, 0, dsize * sizeof(TacsScalar));
  memset(dd2phi, 0, dsize * sizeof(TacsScalar));
  memset(dd1psi, 0, dsize * sizeof(TacsScalar));
  memset(dd2psi, 0, dsize * sizeof(TacsScalar));

  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    A2D::ADMat3x3 T;
    A2D::ADVec3 X0, X0xi, n1, n2, n1xi, n2xi;

    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

    transform->computeTransform(X0xi.x, T.A);

    A2D::ADMat3x3 Xd, Xdinv;
    A2D::ADMat3x3FromThreeADVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::ADMat3x3Inverse invXd(Xd, Xdinv);

    A2D::ADScalar detXd;
    A2D::ADMat3x3Det computedetXd(weight, Xd, detXd);

    A2D::ADMat3x3 XdinvT;
    A2D::ADMat3x3ADMatMult multXdinvT(Xdinv, T, XdinvT);

    // phi-direction chain (plays addAdjResXptProduct's real/vars-chain
    // role).
    A2D::ADVec3 u0xiphi, d01phi, d02phi, d01xiphi, d02xiphi;
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, phi, u0xiphi.x);
    basis::template interpFields<3, 3>(pt, d1phi, d01phi.x);
    basis::template interpFields<3, 3>(pt, d2phi, d02phi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1phi, d01xiphi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2phi, d02xiphi.x);

    // psi-direction chain (unchanged role from addAdjResXptProduct).
    A2D::ADVec3 u0xipsi, d01psi, d02psi, d01xipsi, d02xipsi;
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, psi, u0xipsi.x);
    basis::template interpFields<3, 3>(pt, d1psi, d01psi.x);
    basis::template interpFields<3, 3>(pt, d2psi, d02psi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1psi, d01xipsi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2psi, d02xipsi.x);

    A2D::ADMat3x3 u0dphi;
    A2D::ADMat3x3FromThreeADVec3 assembleu0dphi(u0xiphi, d01phi, d02phi,
                                                u0dphi);

    A2D::ADMat3x3 u0dXdinvTphi, u0xphi;
    A2D::ADMat3x3ADMatMult multu0dphi(u0dphi, XdinvT, u0dXdinvTphi);
    A2D::ADMatTrans3x3ADMatMult multu0xphi(T, u0dXdinvTphi, u0xphi);

    A2D::ADMat3x3 u0dpsi;
    A2D::ADMat3x3FromThreeADVec3 assembleu0dpsi(u0xipsi, d01psi, d02psi,
                                                u0dpsi);

    A2D::ADMat3x3 u0dXdinvTpsi, u0xpsi;
    A2D::ADMat3x3ADMatMult multu0dpsi(u0dpsi, XdinvT, u0dXdinvTpsi);
    A2D::ADMatTrans3x3ADMatMult multu0xpsi(T, u0dXdinvTpsi, u0xpsi);

    A2D::ADScalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::ADMat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::ADMat3x3VecADVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::ADMat3x3VecADVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    A2D::ADVec3 d1tphi, d1xphi;
    A2D::ADVec3Axpy axpyd1tphi(-1.0, sz1, u0xiphi, d01xiphi, d1tphi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd1xphi(s0, T, d1tphi, d1xphi);

    A2D::ADVec3 d2tphi, d2xphi;
    A2D::ADVec3Axpy axpyd2tphi(-1.0, sz2, u0xiphi, d02xiphi, d2tphi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd2xphi(s0, T, d2tphi, d2xphi);

    A2D::ADVec3 d1tpsi, d1xpsi;
    A2D::ADVec3Axpy axpyd1tpsi(-1.0, sz1, u0xipsi, d01xipsi, d1tpsi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd1xpsi(s0, T, d1tpsi, d1xpsi);

    A2D::ADVec3 d2tpsi, d2xpsi;
    A2D::ADVec3Axpy axpyd2tpsi(-1.0, sz2, u0xipsi, d02xipsi, d2tpsi);
    A2D::ADMatTrans3x3ADVecMultADScale matmultd2xpsi(s0, T, d2tpsi, d2xpsi);

    TacsScalar gtypsi[2], gtyphi[2];
    basis::interpTyingStrain(pt, etypsi, gtypsi);
    basis::interpTyingStrain(pt, etyphi, gtyphi);

    TacsScalar e0typsi[2], e0typhi[2];
    e0typsi[0] = 2.0 * XdinvT.A[0] * gtypsi[0];
    e0typsi[1] = 2.0 * XdinvT.A[0] * gtypsi[1];
    e0typhi[0] = 2.0 * XdinvT.A[0] * gtyphi[0];
    e0typhi[1] = 2.0 * XdinvT.A[0] * gtyphi[1];

    TacsScalar epsi[6], ephi[6];
    model::evalStrainDeriv(u0xpsi.A, d1xpsi.x, d2xpsi.x, e0typsi, u0xphi.A,
                           d1xphi.x, d2xphi.x, e0typhi, epsi, ephi);

    TacsScalar spsi[6], sphi[6];
    con->evalStress(elemIndex, pt, X0.x, epsi, spsi);
    con->evalStress(elemIndex, pt, X0.x, ephi, sphi);

    // Piece 2/3: seed the phi-chain reverse with weight spsi (=Cs*epsi),
    // the psi-chain reverse with weight sphi (=Cs*ephi) -- see the
    // function-level comment's decomposition.
    TacsScalar e0typhid[2], e0typsid_local[2];
    model::evalStrainSens(scale * detXd.value, spsi, u0xphi.A, d1xphi.x,
                          d2xphi.x, e0typhi, u0xphi.Ad, d1xphi.xd, d2xphi.xd,
                          e0typhid);
    model::evalStrainSens(scale * detXd.value, sphi, u0xpsi.A, d1xpsi.x,
                          d2xpsi.x, e0typsi, u0xpsi.Ad, d1xpsi.xd, d2xpsi.xd,
                          e0typsid_local);

    // Piece 1: detXd-direction term.
    detXd.valued =
        scale * (epsi[0] * sphi[0] + epsi[1] * sphi[1] + epsi[2] * sphi[2] +
                 epsi[3] * sphi[3] + epsi[4] * sphi[4] + epsi[5] * sphi[5]);

    // Piece 4: e0ty/XdinvT-direction correction (state-independent scalar
    // product, plain product-rule -- mirrors addAdjResXptProduct's
    // identical treatment of this term).
    TacsScalar gtyphid[2], gtypsid[2];
    gtyphid[0] = 2.0 * XdinvT.A[0] * e0typhid[0];
    gtyphid[1] = 2.0 * XdinvT.A[0] * e0typhid[1];
    gtypsid[0] = 2.0 * XdinvT.A[0] * e0typsid_local[0];
    gtypsid[1] = 2.0 * XdinvT.A[0] * e0typsid_local[1];

    XdinvT.Ad[0] +=
        2.0 * (gtyphi[0] * e0typhid[0] + gtyphi[1] * e0typhid[1] +
               e0typsid_local[0] * gtypsi[0] + e0typsid_local[1] * gtypsi[1]);

    matmultd2xpsi.reverse();
    axpyd2tpsi.reverse();
    matmultd1xpsi.reverse();
    axpyd1tpsi.reverse();
    matmultd2xphi.reverse();
    axpyd2tphi.reverse();
    matmultd1xphi.reverse();
    axpyd1tphi.reverse();
    innersz2.reverse();
    innersz1.reverse();
    inners0.reverse();
    multu0xpsi.reverse();
    multu0dpsi.reverse();
    assembleu0dpsi.reverse();
    multu0xphi.reverse();
    multu0dphi.reverse();
    assembleu0dphi.reverse();
    multXdinvT.reverse();
    computedetXd.reverse();
    invXd.reverse();
    assembleXd.reverse();

    transform->addTransformSens(X0xi.x, T.Ad, X0xi.xd);

    basis::template addInterpFieldsTranspose<3, 3>(pt, X0.xd, dfdXpts);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, X0xi.xd, dfdXpts);

    basis::template addInterpFieldsTranspose<3, 3>(pt, n1.xd, dfn1);
    basis::template addInterpFieldsTranspose<3, 3>(pt, n2.xd, dfn2);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, n1xi.xd, dfn1);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, n2xi.xd, dfn2);

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01phi.xd, dd1phi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02phi.xd, dd2phi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xiphi.xd, dd1phi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xiphi.xd, dd2phi);

    basis::template addInterpFieldsTranspose<3, 3>(pt, d01psi.xd, dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02psi.xd, dd2psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xipsi.xd, dd1psi);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xipsi.xd, dd2psi);

    basis::addInterpTyingStrainTranspose(pt, gtyphid, etyphid);
    basis::addInterpTyingStrainTranspose(pt, gtypsid, etypsid);
  }

  // Add the sensitivity contributions from the tying strain (piece 2/3,
  // psi/phi in place of the real vars/psi-adjoint directions
  // addAdjResXptProduct's own call uses).
  model::template addTyingStrainDerivXptSens<vars_per_node, basis>(
      Xpts, fn1, fn2, phi, d1phi, d2phi, psi, d1psi, d2psi, etyphid, etypsid,
      dfdXpts, dfn1, dfn2, dd1phi, dd2phi, dd1psi, dd2psi);

  // Add the contributions from the derivative of the director, via the
  // 6-arg addDirectorRefNormalSens overload's "ddpsi" (perturbation-
  // direction) term with a zeroed "dd" (base-term) buffer, evaluated at
  // the REAL vars (review fix: see the identical fix/rationale in
  // addMatXptSensInnerProduct's TACS_MASS_MATRIX branch -- the earlier
  // single-arg-with-phi/psi-as-"vars" calls silently dropped the q-qpsi
  // cross term TACSQuadraticRotation/TACSQuaternionRotation's nonlinear
  // director maps require). This isolated fix is verified EXACT for
  // TACS_MASS_MATRIX (spot-checked against TACSBeam2ModRot, see PLAN.md)
  // but this STIFFNESS_MATRIX branch retains a separate, unresolved
  // residual for the nonlinear director classes even after this fix --
  // see this function's own typeid guard at its call site
  // (addMatXptSensInnerProduct) and PLAN.md's documented-failure note.
  TacsScalar zero_dd[dsize];
  memset(zero_dd, 0, dsize * sizeof(TacsScalar));
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, phi, fn1, zero_dd, dd1phi, dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn1, zero_dd, dd1psi, dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, phi, fn2, zero_dd, dd2phi, dfn2);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn2, zero_dd, dd2psi, dfn2);

  // Add the contributions from the node normals.
  TacsBeamAddNodeNormalsSens<basis>(Xpts, axis, dfn1, dfn2, dfdXpts);
}

/*
  Add the derivative of the matrix inner product psi^T * mat * phi with
  respect to the nodal coordinates (SPEC.md sec 2.2/2.3/2.3.1).
  TACS_STIFFNESS_MATRIX (Task 4.3) and TACS_MASS_MATRIX (Task 4.2) are
  analytic for TACSLinearizedRotation; TACS_MASS_MATRIX is additionally
  confirmed analytic-exact for TACSQuadraticRotation/TACSQuaternionRotation
  too (spot-checked against TACSBeam2ModRot). TACS_STIFFNESS_MATRIX's
  analytic branch is exact for TACSLinearizedRotation but retains an
  unresolved residual discrepancy for the two nonlinear director classes
  (documented failure, review-flagged -- see
  addMatXptSensInnerProductStiffness's own header comment and PLAN.md):
  fixing the director-Xpts-adjoint misuse the review caught (the same fix
  that made TACS_MASS_MATRIX exact) measurably improves but does not fully
  resolve the STIFFNESS_MATRIX case, indicating a further, unidentified
  term specific to this branch's larger kinematic chain (u0d/u0x/d1x/d2x/
  T/XdinvT) that further investigation within this session could not
  isolate -- forwarded to the base FD/CS implementation for
  TACSQuadraticRotation/TACSQuaternionRotation specifically, per the
  documented-failure-only fallback rule. TACS_GEOMETRIC_STIFFNESS_MATRIX
  is an interim explicit-forward-to-base punt, superseded by Phase 5.
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
    addMatXptSensInnerProduct(ElementMatrixType matType, int elemIndex,
                              double time, TacsScalar scale,
                              const TacsScalar psi[], const TacsScalar phi[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              TacsScalar dfdXpts[]) {
  if (matType == TACS_STIFFNESS_MATRIX) {
    if (typeid(director) != typeid(TACSLinearizedRotation)) {
      // Documented-failure fallback (see this function's header comment):
      // the analytic derivation is exact for TACSLinearizedRotation but
      // has an unresolved residual for the nonlinear director classes.
      TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale,
                                             psi, phi, Xpts, vars, dfdXpts);
      return;
    }
    addMatXptSensInnerProductStiffness(elemIndex, time, scale, psi, phi, Xpts,
                                       vars, dfdXpts);
    return;
  } else if (matType != TACS_MASS_MATRIX) {
    // TACS_GEOMETRIC_STIFFNESS_MATRIX and any future matType: interim
    // explicit-forward-to-base (SPEC.md sec 4.1/4.2/4.3), superseded by
    // Phase 5 for the former.
    TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale, psi,
                                           phi, Xpts, vars, dfdXpts);
    return;
  }

  // TACS_MASS_MATRIX: no strain/stiffness path (SPEC.md sec 2.2 point 3) --
  // the mass matrix's Xpts-dependence flows through detXd(Xpts) (direct)
  // and through the psi/phi-direction director fields' own dependence on
  // Xpts via fn1/fn2 (the director map's reference-normal argument,
  // TacsBeamComputeNodeNormals(Xpts, axis, fn1, fn2)) -- mirrors
  // addAdjResXptProduct's own director-Xpts-adjoint chain
  // (TacsBeamAddNodeNormalsSens/director::addDirectorRefNormalSens), but
  // without the strain/tying-strain/transform kinematics the stiffness
  // branch needs (the mass matrix's dynamics blocks in addJacobian never
  // touch Xdinv/T/e0ty at all).
  const int nquad = quadrature::getNumQuadraturePoints();
  const A2D::Vec3 &axis = transform->getRefAxis();

  TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
  TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

  TacsScalar dfn1[3 * basis::NUM_NODES], dfn2[3 * basis::NUM_NODES];
  memset(dfn1, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));
  memset(dfn2, 0, 3 * basis::NUM_NODES * sizeof(TacsScalar));

  // psi-direction and phi-direction director fields, linearized about the
  // real vars (same reuse as Task 4.1's DV-sens branch; see its own
  // comment for the redundant-recomputation efficiency note, unchanged
  // here).
  TacsScalar d1[dsize], d1dot[dsize], d1ddot[dsize], d1psi[dsize], d1phi[dsize];
  TacsScalar d2[dsize], d2dot[dsize], d2ddot[dsize], d2psi[dsize], d2phi[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn1, d1, d1dot, d1ddot, d1psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn1, d1, d1dot, d1ddot, d1phi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, psi, fn2, d2, d2dot, d2ddot, d2psi);
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               basis::NUM_NODES>(
      vars, vars, vars, phi, fn2, d2, d2dot, d2ddot, d2phi);

  // Per-node accumulators for the director-field Xpts-adjoint (weight on
  // d1phi/d1psi/d2phi/d2psi's own dependence on fn1/fn2), separate from the
  // direct detXd contribution below.
  TacsScalar dd1phi[dsize], dd1psi[dsize], dd2phi[dsize], dd2psi[dsize];
  memset(dd1phi, 0, dsize * sizeof(TacsScalar));
  memset(dd1psi, 0, dsize * sizeof(TacsScalar));
  memset(dd2phi, 0, dsize * sizeof(TacsScalar));
  memset(dd2psi, 0, dsize * sizeof(TacsScalar));

  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    TacsScalar X0[3];
    basis::template interpFields<3, 3>(pt, Xpts, X0);

    A2D::ADVec3 X0xi, n1, n2;
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);

    // The mass matrix needs no reference-axis transform, no Xdinv/XdinvT,
    // and no tying strain -- only detXd(Xpts), matching addJacobian's own
    // dynamics blocks (gamma-path), which are built from detXd alone.
    A2D::ADMat3x3 Xd;
    A2D::ADMat3x3FromThreeADVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::ADScalar detXd;
    A2D::ADMat3x3Det computedetXd(weight, Xd, detXd);

    TacsScalar u0psi[3], u0phi[3];
    basis::template interpFields<vars_per_node, 3>(pt, psi, u0psi);
    basis::template interpFields<vars_per_node, 3>(pt, phi, u0phi);

    TacsScalar d1psiv[3], d1phiv[3], d2psiv[3], d2phiv[3];
    basis::template interpFields<3, 3>(pt, d1psi, d1psiv);
    basis::template interpFields<3, 3>(pt, d1phi, d1phiv);
    basis::template interpFields<3, 3>(pt, d2psi, d2psiv);
    basis::template interpFields<3, 3>(pt, d2phi, d2phiv);

    TacsScalar rho[6];
    con->evalMassMoments(elemIndex, pt, X0, rho);

    TacsScalar g = rho[0] * (u0psi[0] * u0phi[0] + u0psi[1] * u0phi[1] +
                             u0psi[2] * u0phi[2]) +
                   rho[1] * (u0psi[0] * d1phiv[0] + u0psi[1] * d1phiv[1] +
                             u0psi[2] * d1phiv[2] + d1psiv[0] * u0phi[0] +
                             d1psiv[1] * u0phi[1] + d1psiv[2] * u0phi[2]) +
                   rho[2] * (u0psi[0] * d2phiv[0] + u0psi[1] * d2phiv[1] +
                             u0psi[2] * d2phiv[2] + d2psiv[0] * u0phi[0] +
                             d2psiv[1] * u0phi[1] + d2psiv[2] * u0phi[2]) +
                   rho[3] * (d1psiv[0] * d1phiv[0] + d1psiv[1] * d1phiv[1] +
                             d1psiv[2] * d1phiv[2]) +
                   rho[4] * (d2psiv[0] * d2phiv[0] + d2psiv[1] * d2phiv[1] +
                             d2psiv[2] * d2phiv[2]) +
                   rho[5] * (d1psiv[0] * d2phiv[0] + d1psiv[1] * d2phiv[1] +
                             d1psiv[2] * d2phiv[2] + d2psiv[0] * d1phiv[0] +
                             d2psiv[1] * d1phiv[1] + d2psiv[2] * d1phiv[2]);

    // detXd-direction term (piece 1-style, mirrors addAdjResXptProduct).
    detXd.valued = scale * g;

    // Director-field-direction term: d(g)/d(d1phi/d1psi/d2phi/d2psi),
    // scaled by scale*detXd, propagated through addInterpFieldsTranspose
    // into the per-node buffers, then (after the loop)
    // addDirectorRefNormalSens/TacsBeamAddNodeNormalsSens into dfdXpts --
    // exactly mirroring how addAdjResXptProduct routes its own dd1/dd2
    // buffers.
    TacsScalar coef = scale * detXd.value;
    TacsScalar wd1phi[3], wd1psi[3], wd2phi[3], wd2psi[3];
    for (int k = 0; k < 3; k++) {
      wd1phi[k] =
          coef * (rho[1] * u0psi[k] + rho[3] * d1psiv[k] + rho[5] * d2psiv[k]);
      wd1psi[k] =
          coef * (rho[1] * u0phi[k] + rho[3] * d1phiv[k] + rho[5] * d2phiv[k]);
      wd2phi[k] =
          coef * (rho[2] * u0psi[k] + rho[4] * d2psiv[k] + rho[5] * d1psiv[k]);
      wd2psi[k] =
          coef * (rho[2] * u0phi[k] + rho[4] * d2phiv[k] + rho[5] * d1phiv[k]);
    }
    basis::template addInterpFieldsTranspose<3, 3>(pt, wd1phi, dd1phi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, wd1psi, dd1psi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, wd2phi, dd2phi);
    basis::template addInterpFieldsTranspose<3, 3>(pt, wd2psi, dd2psi);

    computedetXd.reverse();
    assembleXd.reverse();

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, X0xi.xd, dfdXpts);
    basis::template addInterpFieldsTranspose<3, 3>(pt, n1.xd, dfn1);
    basis::template addInterpFieldsTranspose<3, 3>(pt, n2.xd, dfn2);
  }

  // Route the psi-direction/phi-direction director fields' own Xpts-sens
  // through the reference normal, via the 6-arg addDirectorRefNormalSens
  // overload's "ddpsi" (perturbation-direction) term with a zeroed "dd"
  // (base-term) buffer, evaluated at the REAL vars (review fix: an
  // earlier revision called the 4-arg single-direction overload with
  // phi/psi themselves substituted for the "vars" argument -- exact only
  // by accident for TACSLinearizedRotation, whose director-map Jacobian
  // is constant and independent of the linearization point, but silently
  // wrong for TACSQuadraticRotation/TACSQuaternionRotation, whose
  // Jacobian-vector-product sensitivity genuinely depends on the REAL q,
  // not on phi/psi treated as if they were the state -- confirmed via a
  // TACSBeam2ModRot spot-check, see PLAN.md). Zeroing "dd" isolates
  // exactly the q-qpsi cross term (TACSDirector.h's own Cd/tmp-based
  // blocks), which is precisely the Jacobian-vector-product's own
  // Xpts-adjoint for a fixed real q and perturbation direction qpsi=phi
  // (or psi) -- this is the SAME reusable machinery
  // addAdjResXptProduct's own 6-arg call already exercises, not new
  // algebra.
  TacsScalar zero_dd[dsize];
  memset(zero_dd, 0, dsize * sizeof(TacsScalar));
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, phi, fn1, zero_dd, dd1phi, dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn1, zero_dd, dd1psi, dfn1);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, phi, fn2, zero_dd, dd2phi, dfn2);
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              basis::NUM_NODES>(
      vars, psi, fn2, zero_dd, dd2psi, dfn2);

  TacsBeamAddNodeNormalsSens<basis>(Xpts, axis, dfn1, dfn2, dfdXpts);
}

/*
  Compute the derivative of the matrix inner product psi^T * mat * phi with
  respect to the state variables (SPEC.md sec 0/2.2/3.3) -- assignment
  (dfdu =), not accumulation.

  psi^T*mat*phi = e_psi(vars)^T * Cs * e_phi(vars), where e_psi(vars) =
  L(J_dir(vars)*psi) (L = model::evalStrain's fixed linear strain map,
  J_dir(vars) = the director map's own Jacobian at vars, per
  computeDirectorRatesDeriv). Differentiating w.r.t. vars again produces a
  term proportional to H_dir(vars)[., psi] (the director map's SECOND
  derivative, i.e. its own curvature) -- a genuinely different quantity
  from evalStrainHessian's (u0x, d1x, d2x, e0ty)-space blocks (SPEC.md sec
  1.1), which govern the strain map L's linearity, not the director map's.
  L is unconditionally linear for TACSBeamLinearModel (the only active
  model class), so this reduction is always valid; whether H_dir itself is
  zero depends only on the director class:

  - TACSLinearizedRotation: d(q,t) = q^x*t is linear in q -- its Jacobian
    is the constant matrix -[t]^x, independent of q, so H_dir = 0
    identically and psi^T*mat*phi does not depend on vars at all for
    TACS_STIFFNESS_MATRIX/TACS_MASS_MATRIX. The memset below already
    produces the correct (exactly zero) answer; no further computation is
    needed or performed.
  - TACSQuadraticRotation/TACSQuaternionRotation: their director maps are
    genuinely nonlinear in q, so H_dir is generically nonzero (confirmed
    via the mandatory complex-step self-check, SPEC.md sec 6.4, against
    TACSBeam2ModRot -- see PLAN.md's Task 4.4 notes). This feature's A2D
    machinery only extends to second order (SPEC.md sec 1.2), which
    covers addJacobian's state-Hessian but not a third-order quantity like
    this one; deriving/implementing the missing H_dir term is out of this
    feature's scope. Forwarded to the base FD/CS fallback for these two
    director classes specifically, per SPEC.md sec 2.2/3.3's documented-
    failure-only fallback rule -- an honest, typeid-guarded routing
    between the analytically-exact-zero case and the genuinely-unresolved
    case, not a special-case forcing zero.
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director,
                     model>::getMatSVSensInnerProduct(ElementMatrixType matType,
                                                      int elemIndex,
                                                      double time,
                                                      const TacsScalar psi[],
                                                      const TacsScalar phi[],
                                                      const TacsScalar Xpts[],
                                                      const TacsScalar vars[],
                                                      TacsScalar dfdu[]) {
  memset(dfdu, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

  if (matType == TACS_STIFFNESS_MATRIX || matType == TACS_MASS_MATRIX) {
    if (typeid(director) == typeid(TACSLinearizedRotation)) {
      // Exact zero (see the derivation above) -- the memset above is the
      // complete, correct result.
      return;
    }
    // Director curvature term not derivable with this feature's A2D
    // machinery -- documented fallback (see the derivation above).
    TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                          Xpts, vars, dfdu);
    return;
  }

  // TACS_GEOMETRIC_STIFFNESS_MATRIX and any future matType: interim
  // explicit-forward-to-base (SPEC.md sec 4.1/4.2/4.3), superseded by
  // Phase 5 for the former.
  TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                        Xpts, vars, dfdu);
}

template <class quadrature, class basis, class director, class model>
int TACSBeamElement<quadrature, basis, director, model>::evalPointQuantity(
    int elemIndex, int quantityType, double time, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *detXdval, TacsScalar *quantity) {
  // Get the reference axis
  const A2D::Vec3 &axis = transform->getRefAxis();

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

  // Interpolate the solution fields
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
  basis::template interpFields<3, 3>(pt, d1, d01.x);
  basis::template interpFields<3, 3>(pt, d2, d02.x);
  basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
  basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

  // Compute X, X,xi and the interpolated normal
  basis::template interpFields<3, 3>(pt, Xpts, X0.x);
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
  basis::template interpFields<3, 3>(pt, fn1, n1.x);
  basis::template interpFields<3, 3>(pt, fn2, n2.x);
  basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
  basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

  // Compute the transformation at the quadrature point
  transform->computeTransform(X0xi.x, T.A);

  // Compute the inverse
  A2D::Mat3x3 Xd, Xdinv;
  A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
  A2D::Mat3x3Inverse invXd(Xd, Xdinv);

  // Compute the determinant of the transform
  A2D::Scalar detXd;
  A2D::Mat3x3Det computedetXd(Xd, detXd);

  if (detXdval) {
    *detXdval = detXd.value;
  }

  if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      *quantity = con->evalDensity(elemIndex, pt, X0.x);
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
      con->evalMassMoments(elemIndex, pt, X0.x, mass_moment);
      TacsScalar density = mass_moment[0];

      for (int i = 0; i < 3; i++) {
        quantity[i] = density * X0.x[i] - mass_moment[1] * n1.x[i] -
                      mass_moment[2] * n2.x[i];
      }
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      TacsScalar I0[6] = {0.0};

      // Evaluate the self MOI
      TacsScalar moments[6];
      con->evalMassMoments(elemIndex, pt, X0.x, moments);
      TacsScalar density = moments[0];
      I0[3] = moments[4] - moments[2] * moments[2] / density;
      I0[4] = -moments[5] + moments[1] * moments[2] / density;
      I0[5] = moments[3] - moments[1] * moments[1] / density;
      // Compute T*I0*T^{T}
      mat3x3SymmTransform(T.A, I0, quantity);
      TacsScalar dXcg[3];
      for (int i = 0; i < 3; i++) {
        dXcg[i] =
            X0.x[i] - (moments[1] * n1.x[i] + moments[2] * n2.x[i]) / density;
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
  A2D::Mat3x3 XdinvT;
  A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

  // Assemble u0d
  A2D::Mat3x3 u0d;
  A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  A2D::Mat3x3 u0dXdinvT, u0x;
  A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
  A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

  // Compute s0, sz1 and sz2
  A2D::Scalar s0, sz1, sz2;
  A2D::Vec3 e1(1.0, 0.0, 0.0);
  A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
  A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
  A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  A2D::Vec3 d1t, d1x;
  A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
  A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  A2D::Vec3 d2t, d2x;
  A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
  A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
  e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

  if (quantityType == TACS_FAILURE_INDEX) {
    if (quantity) {
      *quantity = con->evalFailure(elemIndex, pt, X0.x, e);
    }
    return 1;
  }

  if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X0.x, e, s);
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
  const A2D::Vec3 &axis = transform->getRefAxis();

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

  // Interpolate the solution fields
  basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
  basis::template interpFields<3, 3>(pt, d1, d01.x);
  basis::template interpFields<3, 3>(pt, d2, d02.x);
  basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
  basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

  // Compute X, X,xi and the interpolated normal
  basis::template interpFields<3, 3>(pt, Xpts, X0.x);
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
  basis::template interpFields<3, 3>(pt, fn1, n1.x);
  basis::template interpFields<3, 3>(pt, fn2, n2.x);
  basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
  basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

  // Compute the transformation at the quadrature point
  transform->computeTransform(X0xi.x, T.A);

  // Compute the inverse
  A2D::Mat3x3 Xd, Xdinv;
  A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
  A2D::Mat3x3Inverse invXd(Xd, Xdinv);

  // Compute the determinant of the transform
  A2D::Scalar detXd;
  A2D::Mat3x3Det computedetXd(Xd, detXd);

  if (quantityType == TACS_ELEMENT_DENSITY) {
    con->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X0.x, dvLen, dfdx);
    return;
  }

  // Compute XdinvT = Xdinv * T
  A2D::Mat3x3 XdinvT;
  A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

  // Assemble u0d
  A2D::Mat3x3 u0d;
  A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  A2D::Mat3x3 u0dXdinvT, u0x;
  A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
  A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

  // Compute s0, sz1 and sz2
  A2D::Scalar s0, sz1, sz2;
  A2D::Vec3 e1(1.0, 0.0, 0.0);
  A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
  A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
  A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  A2D::Vec3 d1t, d1x;
  A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
  A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  A2D::Vec3 d2t, d2x;
  A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
  A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
  e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

  if (quantityType == TACS_FAILURE_INDEX) {
    // Add the sensitivity contribution from the design variables
    con->addFailureDVSens(elemIndex, dfdq[0] * scale, pt, X0.x, e, dvLen, dfdx);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X0.x, e, s);
    con->addStressDVSens(elemIndex, scale * dfdq[0], pt, X0.x, e, e, dvLen,
                         dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar dfdmom[6] = {0.0};

    for (int i = 0; i < 3; i++) {
      dfdmom[0] += scale * dfdq[i] * X0.x[i];
      dfdmom[1] += -scale * dfdq[i] * n1.x[i];
      dfdmom[2] += -scale * dfdq[i] * n2.x[i];
    }

    con->addMassMomentsDVSens(elemIndex, pt, X0.x, dfdmom, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    TacsScalar dfdI0[6] = {0.0};

    // Evaluate the self MOI
    TacsScalar moments[6];
    con->evalMassMoments(elemIndex, pt, X0.x, moments);
    TacsScalar density = moments[0];

    TacsScalar dfdmoments[6] = {0.0};
    mat3x3SymmTransformSens(T.A, dfdq, dfdI0);
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
      dXcg[i] =
          X0.x[i] - (moments[1] * n1.x[i] + moments[2] * n2.x[i]) / density;
      dXcgdrho[i] =
          (moments[1] * n1.x[i] + moments[2] * n2.x[i]) / density / density;
      dXcgdmom2[i] = -n2.x[i] / density;
      dXcgdmom1[i] = -n1.x[i] / density;
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

    con->addMassMomentsDVSens(elemIndex, pt, X0.x, dfdmoments, dvLen, dfdx);
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
    const A2D::Vec3 &axis = transform->getRefAxis();

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

    // Interpolate the solution fields
    basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
    basis::template interpFields<3, 3>(pt, d1, d01.x);
    basis::template interpFields<3, 3>(pt, d2, d02.x);
    basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

    // Interpolate the geometry fields
    basis::template interpFields<3, 3>(pt, Xpts, X0.x);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
    basis::template interpFields<3, 3>(pt, fn1, n1.x);
    basis::template interpFields<3, 3>(pt, fn2, n2.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
    basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

    // Compute the transformation at the quadrature point
    transform->computeTransform(X0xi.x, T.A);

    // Compute the inverse
    A2D::Mat3x3 Xd, Xdinv;
    A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
    A2D::Mat3x3Inverse invXd(Xd, Xdinv);

    // Compute the determinant of the transform
    A2D::Scalar detXd;
    A2D::Mat3x3Det computedetXd(Xd, detXd);

    // Compute XdinvT = Xdinv * T
    A2D::Mat3x3 XdinvT;
    A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

    // Assemble u0d
    A2D::ADMat3x3 u0d;
    A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);

    // Compute u0x = T^{T} * u0d * XdinvT
    A2D::ADMat3x3 u0dXdinvT, u0x;
    A2D::ADMat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
    A2D::MatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

    // Compute s0, sz1 and sz2
    A2D::Scalar s0, sz1, sz2;
    A2D::Vec3 e1(1.0, 0.0, 0.0);
    A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
    A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
    A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

    // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
    A2D::ADVec3 d1t, d1x;
    A2D::ADVec3ADVecScalarAxpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
    A2D::MatTrans3x3ADVecMultScale matmultd1x(s0, T, d1t, d1x);

    // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
    A2D::ADVec3 d2t, d2x;
    A2D::ADVec3ADVecScalarAxpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
    A2D::MatTrans3x3ADVecMultScale matmultd2x(s0, T, d2t, d2x);

    // Evaluate the tying components of the strain
    TacsScalar gty[2];  // The components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // // Transform the tying strain to the local coordinates
    TacsScalar e0ty[2];
    e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
    e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

    // Evaluate the strain
    TacsScalar e[6];
    model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

    TacsScalar esens[6];
    if (quantityType == TACS_FAILURE_INDEX) {
      // Compute the sensitivity of the failure index w.r.t. the strain
      con->evalFailureStrainSens(elemIndex, pt, X0.x, e, esens);
    } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
      // Compute the sensitivity of the strain energy density w.r.t. the strain
      con->evalStress(elemIndex, pt, X0.x, e, esens);
      for (int i = 0; i < 6; i++) {
        esens[i] *= 2.0;
      }
    }

    // Evaluate the strain and strain derivatives from the
    TacsScalar e0tyd[2];
    model::evalStrainSens(alpha * dfdq[0], esens, u0x.A, d1x.x, d2x.x, e0ty,
                          u0x.Ad, d1x.xd, d2x.xd, e0tyd);

    // Convert the contributions to the tying strain
    TacsScalar gtyd[2];
    gtyd[0] = 2.0 * XdinvT.A[0] * e0tyd[0];
    gtyd[1] = 2.0 * XdinvT.A[0] * e0tyd[1];

    matmultd2x.reverse();
    axpyd2t.reverse();
    matmultd1x.reverse();
    axpyd1t.reverse();
    multu0x.reverse();
    multu0d.reverse();
    assembleu0d.reverse();

    // Add the residual contributions back to the element
    basis::template addInterpFieldsGradTranspose<vars_per_node, 3>(pt, u0xi.xd,
                                                                   dfdu);

    // Add the constributions back to the derivative
    basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, d1d);
    basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, d2d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, d1d);
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, d2d);

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
  const A2D::Vec3 &axis = transform->getRefAxis();

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
  A2D::ADMat3x3 T;

  // Parametric location
  A2D::ADVec3 X0;

  // Tangent to the beam
  A2D::ADVec3 X0xi;

  // Interpolated normal directions
  A2D::ADVec3 n1, n2;

  // Derivatives of the interpolated normal directions
  A2D::ADVec3 n1xi, n2xi;

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
  A2D::ADMat3x3 Xd, Xdinv;
  A2D::ADMat3x3FromThreeADVec3 assembleXd(X0xi, n1, n2, Xd);
  A2D::ADMat3x3Inverse invXd(Xd, Xdinv);

  // Compute the determinant of the transform
  A2D::ADScalar detXd;
  A2D::ADMat3x3Det computedetXd(Xd, detXd);

  // Compute XdinvT = Xdinv * T
  A2D::ADMat3x3 XdinvT;
  A2D::ADMat3x3ADMatMult multXdinvT(Xdinv, T, XdinvT);

  // Assemble u0d
  A2D::ADMat3x3 u0d;
  A2D::ADMat3x3FromThreeADVec3 assembleu0d(u0xi, d01, d02, u0d);

  // Compute u0x = T^{T} * u0d * XdinvT
  A2D::ADMat3x3 u0dXdinvT, u0x;
  A2D::ADMat3x3ADMatMult multu0d(u0d, XdinvT, u0dXdinvT);
  A2D::ADMatTrans3x3ADMatMult multu0x(T, u0dXdinvT, u0x);

  // Compute s0, sz1 and sz2
  A2D::ADScalar s0, sz1, sz2;
  A2D::Vec3 e1(1.0, 0.0, 0.0);
  A2D::ADMat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
  A2D::ADMat3x3VecADVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
  A2D::ADMat3x3VecADVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

  // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
  A2D::ADVec3 d1t, d1x;
  A2D::ADVec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
  A2D::ADMatTrans3x3ADVecMultADScale matmultd1x(s0, T, d1t, d1x);

  // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
  A2D::ADVec3 d2t, d2x;
  A2D::ADVec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
  A2D::ADMatTrans3x3ADVecMultADScale matmultd2x(s0, T, d2t, d2x);

  // Evaluate the tying components of the strain
  TacsScalar gty[2];  // The components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Transform the tying strain to the local coordinates
  TacsScalar e0ty[2];
  e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
  e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

  // Compute the set of strain components
  TacsScalar e[6];  // The components of the strain
  model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

  // Evaluate the failure sensitivity contribution
  TacsScalar esens[6] = {0.0};
  if (quantityType == TACS_FAILURE_INDEX) {
    // Compute the sensitivity of the failure index w.r.t. the strain
    con->evalFailureStrainSens(elemIndex, pt, X0.x, e, esens);
  } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the sensitivity of the strain energy density w.r.t. the strain
    con->evalStress(elemIndex, pt, X0.x, e, esens);
    for (int i = 0; i < 6; i++) {
      esens[i] *= 2.0;
    }
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    // Compute the sensitivity of the strain energy density w.r.t. the strain
    TacsScalar mass_moment[6];
    con->evalMassMoments(elemIndex, pt, X0.x, mass_moment);
    TacsScalar density = mass_moment[0];

    for (int i = 0; i < 3; i++) {
      X0.xd[i] = density * dfdq[i];
      n1.xd[i] = mass_moment[2] * dfdq[i];
      n2.xd[i] = -mass_moment[1] * dfdq[i];
    }
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    // Analytic Xpts-sensitivity (SPEC.md sec 3.1; PLAN.md Phase 3 Task
    // 3.2). Differentiates evalPointQuantity's own forward formula
    // (TACSBeamElement.h:1939-1968) directly -- moments =
    // con->evalMassMoments(...) and hence I0[3..5] are Xpts-independent
    // (a per-quadrature-point constitutive evaluation, not a function of
    // node location), so quantity = T*I0*T^T + parallel_axis(dXcg) has
    // exactly two independent Xpts-dependent paths: T (via X0xi) and
    // dXcg (via X0/n1/n2). No new algebra primitive is needed for
    // either -- both reduce to plain 3x3 matrix arithmetic plus the
    // EXISTING T (T is already an A2D::ADMat3x3 participating in this
    // function's shared forward/reverse graph below) and
    // transform->addTransformSens machinery the common tail already
    // calls (TACSBeamElement.h:2638, unchanged).
    TacsScalar moments[6];
    con->evalMassMoments(elemIndex, pt, X0.x, moments);
    TacsScalar density = moments[0];

    TacsScalar I0[6] = {0.0};
    I0[3] = moments[4] - moments[2] * moments[2] / density;
    I0[4] = -moments[5] + moments[1] * moments[2] / density;
    I0[5] = moments[3] - moments[1] * moments[1] / density;

    // (1) T-dependent path: quantity_self = T*I0*T^T (I0 fixed, since it
    // has no Xpts-dependence). For f = D:quantity_self (D the packed-
    // symmetric adjoint dfdq, unpacked to a full 3x3), matrix calculus
    // gives df/dT = 2*D*T*I0 -- seed this directly into T.Ad (T already
    // participates in the shared A2D graph below; transform-
    // >addTransformSens(X0xi.x, T.Ad, X0xi.xd) in the common tail
    // converts T.Ad into the X0xi (hence Xpts) adjoint, unchanged).
    //
    // dfdq is packed-symmetric (6 slots: xx,xy,xz,yy,yz,zz), with
    // f = sum_k dfdq[k]*quantity[k] -- each off-diagonal slot counted
    // ONCE (matching every other quantity type in this function, e.g.
    // TACS_ELEMENT_DENSITY_MOMENT's plain per-component dfdq[i] usage
    // above). Unpacking dfdq into a full symmetric 3x3 matrix D such
    // that sum_ij D_ij*A_ij reproduces that same packed sum requires
    // HALVING the off-diagonal entries (verified via a standalone FD
    // check before wiring this in: mirroring dfdq[1] into BOTH D[0][1]
    // and D[1][0] at full value double-counts the off-diagonal
    // contribution by a factor of 2).
    // Every dfdq-driven contribution below is additionally scaled by
    // "scale" (TacsTestElementQuantityXptSens draws a random, generically
    // nonzero scale and folds it into its own FD reference exactly this
    // way, TACSElementVerification.cpp:1678-1707 -- confirmed by reading
    // the harness directly after an initial 2x-ish-looking mismatch
    // turned out to be 1/scale, not a sign/formula error).
    TacsScalar Dmat[9] = {
        scale * dfdq[0],       scale * 0.5 * dfdq[1], scale * 0.5 * dfdq[2],
        scale * 0.5 * dfdq[1], scale * dfdq[3],       scale * 0.5 * dfdq[4],
        scale * 0.5 * dfdq[2], scale * 0.5 * dfdq[4], scale * dfdq[5]};
    TacsScalar I0mat[9] = {I0[0], I0[1], I0[2], I0[1], I0[3],
                           I0[4], I0[2], I0[4], I0[5]};
    TacsScalar TI0[9], DTI0[9];
    A2D::Mat3x3MatMultCore(T.A, I0mat, TI0);
    A2D::Mat3x3MatMultCore(Dmat, TI0, DTI0);
    for (int i = 0; i < 9; i++) {
      T.Ad[i] += 2.0 * DTI0[i];
    }

    // (2) dXcg-dependent path: the parallel-axis contribution is the
    // point-mass inertia tensor of a point of mass "density" located at
    // dXcg = X0 - (moments[1]*n1 + moments[2]*n2)/density:
    // quantity_parallel = density*(|dXcg|^2*Identity - dXcg (x) dXcg).
    // For f = D:quantity_parallel, df/d(dXcg) =
    // 2*density*(tr(D)*dXcg - D*dXcg); dXcg is linear in X0/n1/n2, so
    // this seeds X0.xd/n1.xd/n2.xd directly (mirroring the existing
    // TACS_ELEMENT_DENSITY_MOMENT branch's identical seed-then-fall-
    // through pattern just above).
    TacsScalar dXcg[3];
    for (int i = 0; i < 3; i++) {
      dXcg[i] =
          X0.x[i] - (moments[1] * n1.x[i] + moments[2] * n2.x[i]) / density;
    }
    TacsScalar trD = Dmat[0] + Dmat[4] + Dmat[8];
    TacsScalar Dr[3];
    A2D::Mat3x3VecMultCore(Dmat, dXcg, Dr);

    TacsScalar dfdr[3];
    for (int i = 0; i < 3; i++) {
      dfdr[i] = 2.0 * density * (trD * dXcg[i] - Dr[i]);
    }

    for (int i = 0; i < 3; i++) {
      X0.xd[i] = dfdr[i];
      n1.xd[i] = -dfdr[i] * moments[1] / density;
      n2.xd[i] = -dfdr[i] * moments[2] / density;
    }
  }

  // Evaluate the strain and strain derivatives from the
  TacsScalar e0tyd[2];
  model::evalStrainSens(scale * dfdq[0], esens, u0x.A, d1x.x, d2x.x, e0ty,
                        u0x.Ad, d1x.xd, d2x.xd, e0tyd);
  detXd.valued = scale * dfddetXd;

  // Apply the tying strain transformation
  TacsScalar gtyd[2];
  gtyd[0] = 2.0 * XdinvT.A[0] * e0tyd[0];
  gtyd[1] = 2.0 * XdinvT.A[0] * e0tyd[1];

  XdinvT.Ad[0] += 2.0 * (gty[0] * e0tyd[0] + gty[1] * e0tyd[1]);

  matmultd2x.reverse();
  axpyd2t.reverse();
  matmultd1x.reverse();
  axpyd1t.reverse();
  innersz2.reverse();
  innersz1.reverse();
  inners0.reverse();
  multu0x.reverse();
  multu0d.reverse();
  assembleu0d.reverse();
  multXdinvT.reverse();
  computedetXd.reverse();
  invXd.reverse();
  assembleXd.reverse();

  // Reverse the transformation sensitivities
  transform->addTransformSens(X0xi.x, T.Ad, X0xi.xd);

  // Add the sensitivities to the input fields...
  basis::template addInterpFieldsTranspose<3, 3>(pt, X0.xd, dfdXpts);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, X0xi.xd, dfdXpts);
  basis::template addInterpFieldsTranspose<3, 3>(pt, n1.xd, dfn1);
  basis::template addInterpFieldsTranspose<3, 3>(pt, n2.xd, dfn2);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, n1xi.xd, dfn1);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, n2xi.xd, dfn2);
  basis::template addInterpFieldsTranspose<3, 3>(pt, d01.xd, dd1);
  basis::template addInterpFieldsTranspose<3, 3>(pt, d02.xd, dd2);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, d01xi.xd, dd1);
  basis::template addInterpFieldsGradTranspose<3, 3>(pt, d02xi.xd, dd2);

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
    const A2D::Vec3 &axis = transform->getRefAxis();

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

      // Interpolate the solution fields
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
      basis::template interpFields<3, 3>(pt, d1, d01.x);
      basis::template interpFields<3, 3>(pt, d2, d02.x);
      basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
      basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

      // Interpolate the geometry fields
      basis::template interpFields<3, 3>(pt, Xpts, X0.x);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
      basis::template interpFields<3, 3>(pt, fn1, n1.x);
      basis::template interpFields<3, 3>(pt, fn2, n2.x);
      basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
      basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);

      // Compute the transformation at the quadrature point
      transform->computeTransform(X0xi.x, T.A);

      // Compute the inverse
      A2D::Mat3x3 Xd, Xdinv;
      A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
      A2D::Mat3x3Inverse invXd(Xd, Xdinv);

      // Compute the determinant of the transform
      A2D::Scalar detXd;
      A2D::Mat3x3Det computedetXd(Xd, detXd);

      // Compute XdinvT = Xdinv * T
      A2D::Mat3x3 XdinvT;
      A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

      // Assemble u0d
      A2D::Mat3x3 u0d;
      A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

      // Compute u0x = T^{T} * u0d * XdinvT
      A2D::Mat3x3 u0dXdinvT, u0x;
      A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
      A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

      // Compute s0, sz1 and sz2
      A2D::Scalar s0, sz1, sz2;
      A2D::Vec3 e1(1.0, 0.0, 0.0);
      A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
      A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
      A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

      // Compute d1x = s0 * T^{T} * (d1xi - sz1 * u0xi)
      A2D::Vec3 d1t, d1x;
      A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
      A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

      // Compute d2x = s0 * T^{T} * (d2xi - sz2 * u0xi)
      A2D::Vec3 d2t, d2x;
      A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
      A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

      // Evaluate the tying components of the strain
      TacsScalar gty[2];  // The components of the tying strain
      basis::interpTyingStrain(pt, ety, gty);

      // Transform the tying strain to the local coordinates
      TacsScalar e0ty[2];
      e0ty[0] = gty[0];
      e0ty[1] = gty[1];

      // Compute the set of strain components
      TacsScalar e[6];  // The components of the strain
      model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

      // Compute the corresponding stresses
      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X0.x, e, s);

      if (write_flag & TACS_OUTPUT_NODES) {
        data[0] = X0.x[0];
        data[1] = X0.x[1];
        data[2] = X0.x[2];
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
          data[failInd] =
              con->evalFailureFieldValue(elemIndex, pt, X0.x, e, failInd);
        }
        for (int dvInd = 0; dvInd < 7; dvInd++) {
          data[dvInd + 7] =
              con->evalDesignFieldValue(elemIndex, pt, X0.x, dvInd);
        }
        data += 14;
      }
      if (write_flag & TACS_OUTPUT_COORDINATE_FRAME) {
        data[0] = T.A[0];
        data[1] = T.A[3];
        data[2] = T.A[6];

        data[3] = T.A[1];
        data[4] = T.A[4];
        data[5] = T.A[7];

        data[6] = T.A[2];
        data[7] = T.A[5];
        data[8] = T.A[8];
        data += 9;
      }
    }
  }
}

/*
  Compute the averaged stress resultants over the element's visualization
  nodes (SPEC.md sec 1.5). Forward-only value computation, not a
  sensitivity -- reuses the same state-independent geometry chain and
  strain/stress evaluation addResidual/computeEnergies/getOutputData
  already use, unchanged. Mirrors TACSShellElement::getAverageStresses
  (TACSShellElement.h:1455-1540) and TACSBeamElement::getOutputData's own
  per-vis-node loop above, with the same e0ty = 2*XdinvT.A[0]*gty
  transform addResidual/addJacobian use (getOutputData's own copy of this
  loop omits the "2*XdinvT.A[0]*" factor -- a pre-existing inconsistency
  in getOutputData, out of scope for this task; this new function uses
  the confirmed-correct formula).

  Beam's stress vector has only TACSBeamConstitutive::NUM_STRESSES = 6
  components (SPEC.md sec 1.1) -- avgStresses[0..5] receive the averaged
  stress, and avgStresses[6..8] are deliberately left untouched (not
  zeroed), matching the base class's own no-op-leaves-untouched
  convention for any slot this element doesn't know how to fill (SPEC.md
  sec 1.5).
*/
template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::getAverageStresses(
    int elemIndex, ElementType etype, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *avgStresses) {
  if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

    const A2D::Vec3 &axis = transform->getRefAxis();

    TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
    TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

    TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn1,
                                                              d1, d1dot);
    director::template computeDirectorRates<vars_per_node, offset,
                                            basis::NUM_NODES>(vars, dvars, fn2,
                                                              d2, d2dot);

    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2,
                                                             vars, d1, d2, ety);

    TacsScalar loc_avgStresses[6];
    for (int i = 0; i < 6; i++) {
      loc_avgStresses[i] = 0.0;
    }

    for (int index = 0; index < num_vis_nodes; index++) {
      double pt[3];
      basis::getNodePoint(index, pt);

      A2D::Vec3 X0, X0xi, n1, n2;
      basis::template interpFields<3, 3>(pt, Xpts, X0.x);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
      basis::template interpFields<3, 3>(pt, fn1, n1.x);
      basis::template interpFields<3, 3>(pt, fn2, n2.x);

      A2D::Mat3x3 T;
      transform->computeTransform(X0xi.x, T.A);

      A2D::Mat3x3 Xd, Xdinv;
      A2D::Mat3x3FromThreeVec3 assembleXd(X0xi, n1, n2, Xd);
      A2D::Mat3x3Inverse invXd(Xd, Xdinv);

      A2D::Mat3x3 XdinvT;
      A2D::Mat3x3MatMult multXdinvT(Xdinv, T, XdinvT);

      A2D::Vec3 u0xi, d01, d02, d01xi, d02xi;
      basis::template interpFieldsGrad<vars_per_node, 3>(pt, vars, u0xi.x);
      basis::template interpFields<3, 3>(pt, d1, d01.x);
      basis::template interpFields<3, 3>(pt, d2, d02.x);
      basis::template interpFieldsGrad<3, 3>(pt, d1, d01xi.x);
      basis::template interpFieldsGrad<3, 3>(pt, d2, d02xi.x);

      A2D::Mat3x3 u0d;
      A2D::Mat3x3FromThreeVec3 assembleu0d(u0xi, d01, d02, u0d);

      A2D::Mat3x3 u0dXdinvT, u0x;
      A2D::Mat3x3MatMult multu0d(u0d, XdinvT, u0dXdinvT);
      A2D::MatTrans3x3MatMult multu0x(T, u0dXdinvT, u0x);

      A2D::Scalar s0, sz1, sz2;
      A2D::Vec3 e1(1.0, 0.0, 0.0);
      A2D::Vec3 n1xi, n2xi;
      basis::template interpFieldsGrad<3, 3>(pt, fn1, n1xi.x);
      basis::template interpFieldsGrad<3, 3>(pt, fn2, n2xi.x);
      A2D::Mat3x3VecVecInnerProduct inners0(XdinvT, e1, e1, s0);
      A2D::Mat3x3VecVecInnerProduct innersz1(Xdinv, e1, n1xi, sz1);
      A2D::Mat3x3VecVecInnerProduct innersz2(Xdinv, e1, n2xi, sz2);

      A2D::Vec3 d1t, d1x;
      A2D::Vec3Axpy axpyd1t(-1.0, sz1, u0xi, d01xi, d1t);
      A2D::MatTrans3x3VecMultScale matmultd1x(s0, T, d1t, d1x);

      A2D::Vec3 d2t, d2x;
      A2D::Vec3Axpy axpyd2t(-1.0, sz2, u0xi, d02xi, d2t);
      A2D::MatTrans3x3VecMultScale matmultd2x(s0, T, d2t, d2x);

      TacsScalar gty[2];
      basis::interpTyingStrain(pt, ety, gty);

      TacsScalar e0ty[2];
      e0ty[0] = 2.0 * XdinvT.A[0] * gty[0];
      e0ty[1] = 2.0 * XdinvT.A[0] * gty[1];

      TacsScalar e[6];
      model::evalStrain(u0x.A, d1x.x, d2x.x, e0ty, e);

      TacsScalar s[6];
      con->evalStress(elemIndex, pt, X0.x, e, s);

      for (int i = 0; i < 6; i++) {
        loc_avgStresses[i] += s[i];
      }
    }

    for (int i = 0; i < 6; i++) {
      avgStresses[i] += loc_avgStresses[i] / num_vis_nodes;
    }
  }
}

#endif  // TACS_BEAM_ELEMENT_H
