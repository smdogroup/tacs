#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

#include "TACSDirector.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSElementVerification.h"
#include "TACSShellCentrifugalForce.h"
#include "TACSShellConstitutive.h"
#include "TACSShellElementModel.h"
#include "TACSShellElementTransform.h"
#include "TACSShellInertialForce.h"
#include "TACSShellPressure.h"
#include "TACSShellTraction.h"
#include "TACSShellUtilities.h"

template <class quadrature, class basis, class director, class model>
class TACSShellElement : public TACSElement {
 public:
  // Offset within the solution vector to the rotational
  // parametrization defined via the director class. Here the offset
  // is 3 corresponding to the (u, v, w) displacements of the
  // mid-surface of the shell.
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

  TACSShellElement(TACSShellTransform *_transform,
                   TACSShellConstitutive *_con) {
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  ~TACSShellElement() {
    if (transform) {
      transform->decref();
    }

    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSShellElement"; }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return num_nodes; }

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
    return new TACSShellTraction<vars_per_node, quadrature, basis>(t);
  }

  TACSElement *createElementPressure(int faceIndex, TacsScalar p) {
    return new TACSShellPressure<vars_per_node, quadrature, basis>(p);
  }

  TACSElement *createElementInertialForce(const TacsScalar inertiaVec[]) {
    return new TACSShellInertialForce<vars_per_node, quadrature, basis>(
        con, inertiaVec);
  }

  TACSElement *createElementCentrifugalForce(const TacsScalar omega[],
                                             const TacsScalar rotCenter[]) {
    return new TACSShellCentrifugalForce<vars_per_node, quadrature, basis>(
        con, omega, rotCenter);
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

  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

 private:
  // Set sizes for the different components
  static const int usize = 3 * num_nodes;
  static const int dsize = 3 * num_nodes;
  static const int csize = 9 * num_nodes;

  TACSShellTransform *transform;
  TACSShellConstitutive *con;
};

/*
  Compute the kinetic and potential energies of the shell
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::computeEnergies(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, TacsScalar *Te, TacsScalar *Ue) {
  // Zero the kinetic and potential energies
  TacsScalar Telem = 0.0;
  TacsScalar Uelem = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar etn[num_nodes];
  TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
  TacsScalar u0xn[9 * num_nodes], Ctn[csize];
  TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
      transform, Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  // Compute the director rates
  TacsScalar d[dsize], ddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, fn, d, ddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                           ety);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9], et;
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);
    basis::template interpFields<1, 1>(pt, etn, &et);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9];
    TacsScalar detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9];  // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = et;

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    Uelem +=
        0.5 * detXd *
        (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] + s[4] * e[4] +
         s[5] * e[5] + s[6] * e[6] + s[7] * e[7] + s[8] * e[8]);

    // Evaluate the mass moments
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    // Compute the velocities and the director velocities
    TacsScalar u0dot[3], d0dot[3];
    basis::template interpFields<vars_per_node, 3>(pt, dvars, u0dot);
    basis::template interpFields<3, 3>(pt, ddot, d0dot);

    Telem += 0.5 * detXd *
             (moments[0] * vec3Dot(u0dot, u0dot) +
              2.0 * moments[1] * vec3Dot(u0dot, d0dot) +
              moments[2] * vec3Dot(d0dot, d0dot));
  }

  *Te = Telem;
  *Ue = Uelem;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::addResidual(
    int elemIndex, double time, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar res[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Derivative of the director field and matrix at each point
  TacsScalar dd[dsize];
  memset(dd, 0, 3 * num_nodes * sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes], detn[num_nodes];
  memset(detn, 0, num_nodes * sizeof(TacsScalar));

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
  TacsScalar u0xn[9 * num_nodes], Ctn[csize];
  TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
      transform, Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn, d, ddot, dddot);

  // Compute the tying strain values
  TacsScalar ety[basis::NUM_TYING_POINTS], dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                           ety);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9], et;
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);
    basis::template interpFields<1, 1>(pt, etn, &et);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9];
    TacsScalar detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9];  // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = et;

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    // Add the contribution to the drilling strain
    TacsScalar det = detXd * s[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det, detn);

    // Add the contributions to the residual from du0x, du1x and dCt
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT, du0x,
                                                   du1x, res, dd);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    // Evaluate the mass moments
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    // Evaluate the second time derivatives
    TacsScalar u0ddot[3], d0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot);
    basis::template interpFields<3, 3>(pt, dddot, d0ddot);

    // Add the contributions to the derivative
    TacsScalar du0dot[3];
    du0dot[0] = detXd * (moments[0] * u0ddot[0] + moments[1] * d0ddot[0]);
    du0dot[1] = detXd * (moments[0] * u0ddot[1] + moments[1] * d0ddot[1]);
    du0dot[2] = detXd * (moments[0] * u0ddot[2] + moments[1] * d0ddot[2]);
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, du0dot, res);

    TacsScalar dd0dot[3];
    dd0dot[0] = detXd * (moments[1] * u0ddot[0] + moments[2] * d0ddot[0]);
    dd0dot[1] = detXd * (moments[1] * u0ddot[1] + moments[2] * d0ddot[1]);
    dd0dot[2] = detXd * (moments[1] * u0ddot[2] + moments[2] * d0ddot[2]);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dd0dot, dd);
  }

  // Add the contribution to the residual from the drill strain
  TacsShellAddDrillStrainSens<vars_per_node, offset, basis, director, model>(
      Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, detn, res);

  // Add the contributions from the tying strain
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, res, dd);

  // Add the contributions to the director field
  director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn, dd, res);

  // Add the contribution from the rotation constraint (defined by the
  // rotational parametrization) - if any
  director::template addRotationConstraint<vars_per_node, offset, num_nodes>(
      vars, res);
}

/*
  Add the contributions to the residual and Jacobian matrix
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::addJacobian(
    int elemIndex, double time, TacsScalar alpha, TacsScalar beta,
    TacsScalar gamma, const TacsScalar Xpts[], const TacsScalar vars[],
    const TacsScalar dvars[], const TacsScalar ddvars[], TacsScalar res[],
    TacsScalar mat[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Derivative of the director field
  TacsScalar dd[dsize];
  memset(dd, 0, dsize * sizeof(TacsScalar));

  // Second derivatives required for the director
  TacsScalar d2d[dsize * dsize], d2du[usize * dsize];
  TacsScalar d2Tdotd[dsize * dsize], d2Tdotu[usize * dsize];
  memset(d2d, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2du, 0, usize * dsize * sizeof(TacsScalar));
  memset(d2Tdotd, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2Tdotu, 0, usize * dsize * sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  TacsScalar d2ety[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
  TacsScalar d2etyu[basis::NUM_TYING_POINTS * usize];
  TacsScalar d2etyd[basis::NUM_TYING_POINTS * dsize];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(
      d2ety, 0,
      basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(d2etyu, 0, basis::NUM_TYING_POINTS * usize * sizeof(TacsScalar));
  memset(d2etyd, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes], detn[num_nodes];
  TacsScalar d2etn[num_nodes * num_nodes];
  memset(detn, 0, num_nodes * sizeof(TacsScalar));
  memset(d2etn, 0, num_nodes * num_nodes * sizeof(TacsScalar));

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
  TacsScalar u0xn[9 * num_nodes], Ctn[csize];
  TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
      transform, Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
      vars, dvars, ddvars, fn, d, ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                           ety);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9], et;
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);
    basis::template interpFields<1, 1>(pt, etn, &et);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9];
    TacsScalar detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9];  // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = et;

    // Compute the tangent stiffness matrix
    TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
    con->evalTangentStiffness(elemIndex, pt, X, Cs);

    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    // Compute the stress based on the tangent stiffness
    TacsScalar s[9];
    TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
    TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
    model::evalStrainHessian(alpha * detXd, s, Cs, u0x, u1x, e0ty, d2u0x, d2u1x,
                             d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x);

    // Add the contributions to the residual from du0x and du1x
    TacsScalar det = detXd * s[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det, detn);

    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT, du0x,
                                                   du1x, res, dd);

    // Add the contribution from the drilling stiffness
    TacsScalar d2et = detXd * alpha * Cs[21];
    basis::template addInterpFieldsOuterProduct<1, 1, 1, 1>(pt, &d2et, d2etn);

    // Add the contributions to the residual from du0x, du1x and dCt
    TacsShellAddDispGradHessian<vars_per_node, basis>(
        pt, T, XdinvT, XdinvzT, d2u0x, d2u1x, d2u0xu1x, mat, d2d, d2du);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6], d2gty[36];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
    mat3x3SymmTransformTransHessian(XdinvT, d2e0ty, d2gty);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);
    basis::addInterpTyingStrainHessian(pt, d2gty, d2ety);

    // Add the coupling between the displacement and tying strain
    TacsShellAddTyingDispCoupling<basis>(pt, T, XdinvT, XdinvzT, d2e0tyu0x,
                                         d2e0tyu1x, d2etyu, d2etyd);

    // Evaluate the mass moments
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    // Evaluate the second time derivatives
    TacsScalar u0ddot[3], d0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot);
    basis::template interpFields<3, 3>(pt, dddot, d0ddot);

    // Add the contributions to the derivative
    TacsScalar du0dot[3];
    du0dot[0] = detXd * (moments[0] * u0ddot[0] + moments[1] * d0ddot[0]);
    du0dot[1] = detXd * (moments[0] * u0ddot[1] + moments[1] * d0ddot[1]);
    du0dot[2] = detXd * (moments[0] * u0ddot[2] + moments[1] * d0ddot[2]);
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, du0dot, res);

    TacsScalar dd0dot[3];
    dd0dot[0] = detXd * (moments[1] * u0ddot[0] + moments[2] * d0ddot[0]);
    dd0dot[1] = detXd * (moments[1] * u0ddot[1] + moments[2] * d0ddot[1]);
    dd0dot[2] = detXd * (moments[1] * u0ddot[2] + moments[2] * d0ddot[2]);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dd0dot, dd);

    TacsScalar d2u0dot[9];
    memset(d2u0dot, 0, 9 * sizeof(TacsScalar));
    d2u0dot[0] = d2u0dot[4] = d2u0dot[8] = gamma * detXd * moments[0];
    basis::template addInterpFieldsOuterProduct<vars_per_node, vars_per_node, 3,
                                                3>(pt, d2u0dot, mat);

    TacsScalar d2Td[9];
    memset(d2Td, 0, 9 * sizeof(TacsScalar));
    d2Td[0] = d2Td[4] = d2Td[8] = detXd * moments[2];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2Td, d2Tdotd);

    d2Td[0] = d2Td[4] = d2Td[8] = detXd * moments[1];
    basis::template addInterpFieldsOuterProduct<3, 3, 3, 3>(pt, d2Td, d2Tdotu);
  }

  // Add the contribution to the residual from the drill strain
  TacsShellAddDrillStrainHessian<vars_per_node, offset, basis, director, model>(
      Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, detn, d2etn, res, mat);

  // Add the residual from the tying strain
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, res, dd);

  // Add the second order terms from the tying strain
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
      alpha, Xpts, fn, vars, d, dety, d2ety, d2etyu, d2etyd, mat, d2d, d2du);

  // Add the contributions to the stiffness matrix
  director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
      alpha, beta, gamma, vars, dvars, ddvars, fn, dd, d2Tdotd, d2Tdotu, d2d,
      d2du, res, mat);

  // Add the constraint associated with the rotational parametrization (if any)
  director::template addRotationConstrJacobian<vars_per_node, offset,
                                               num_nodes>(alpha, vars, res,
                                                          mat);
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::getMatType(
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
    // Not implimented
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
void TACSShellElement<quadrature, basis, director, model>::addAdjResProduct(
    int elemIndex, double time, TacsScalar scale, const TacsScalar psi[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int dvLen, TacsScalar dfdx[]) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar etn[num_nodes], etnd[num_nodes];
  TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
  TacsScalar u0xn[9 * num_nodes], Ctn[csize];
  TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                   model>(transform, Xdn, fn, vars, psi,
                                          XdinvTn, Tn, u0xn, Ctn, etn, etnd);

  // Compute the director rates and their derivatives
  TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                               num_nodes>(
      vars, dvars, ddvars, psi, fn, d, ddot, dddot, dd);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS], etyd[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn, vars, d, psi, dd, ety, etyd);

  // Loop over each quadrature point and add the residual contribution
  for (int quad_index = 0; quad_index < nquad; quad_index++) {
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9], et, etd;
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);
    basis::template interpFields<1, 1>(pt, etn, &et);
    basis::template interpFields<1, 1>(pt, etnd, &etd);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], u0xd[9], u1xd[9];
    TacsScalar detXd = TacsShellComputeDispGradDeriv<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, psi, dd, XdinvT, XdinvzT, u0x, u1x,
        u0xd, u1xd);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6], gtyd[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);
    basis::interpTyingStrain(pt, etyd, gtyd);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6], e0tyd[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
    mat3x3SymmTransformTranspose(XdinvT, gtyd, e0tyd);

    // Compute the set of strain components
    TacsScalar e[9];   // The components of the strain
    TacsScalar ed[9];  // The directional derivative components of the strain
    model::evalStrainDeriv(u0x, u1x, e0ty, u0xd, u1xd, e0tyd, e, ed);
    e[8] = et;
    ed[8] = etd;

    // The directional derivative of the strain along the adjoint direction
    con->addStressDVSens(elemIndex, scale * detXd, pt, X, e, ed, dvLen, dfdx);

    // Evaluate the second time derivatives
    TacsScalar u0ddot[3], d0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot);
    basis::template interpFields<3, 3>(pt, dddot, d0ddot);

    TacsScalar du0ddot[3], dd0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, psi, du0ddot);
    basis::template interpFields<3, 3>(pt, dd, dd0ddot);

    TacsScalar coef[3];
    coef[0] = scale * detXd * vec3Dot(u0ddot, du0ddot);
    coef[1] =
        scale * detXd * (vec3Dot(u0ddot, dd0ddot) + vec3Dot(du0ddot, d0ddot));
    coef[2] = scale * detXd * vec3Dot(d0ddot, dd0ddot);

    // Add the contribution from the dynamics
    con->addMassMomentsDVSens(elemIndex, pt, X, coef, dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
int TACSShellElement<quadrature, basis, director, model>::evalPointQuantity(
    int elemIndex, int quantityType, double time, int n, double pt[],
    const TacsScalar Xpts[], const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *detXd, TacsScalar *quantity) {
  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn);

  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    if (quantity) {
      // Compute the director rates
      TacsScalar d[dsize], ddot[dsize];
      director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
          vars, dvars, fn, d, ddot);

      // Set the total number of tying points needed for this element
      TacsScalar ety[basis::NUM_TYING_POINTS];
      model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars,
                                                               d, ety);

      // Compute X, X,xi and the interpolated normal n0
      TacsScalar X[3], Xxi[6], n0[3], T[9];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      // Compute the transformation at the quadrature point
      transform->computeTransform(Xxi, n0, T);

      // Evaluate the displacement gradient at the point
      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9];
      *detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);

      // Evaluate the tying components of the strain
      TacsScalar gty[6];  // The symmetric components of the tying strain
      basis::interpTyingStrain(pt, ety, gty);

      // Compute the symmetric parts of the tying strain
      TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

      // Compute the set of strain components
      TacsScalar e[9];  // The components of the strain
      model::evalStrain(u0x, u1x, e0ty, e);
      e[8] = 0.0;

      if (quantityType == TACS_FAILURE_INDEX) {
        *quantity = con->evalFailure(elemIndex, pt, X, e);
      } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
        TacsScalar s[9];
        con->evalStress(elemIndex, pt, X, e, s);
        *quantity = 0.0;
        for (int i = 0; i < 9; i++) {
          *quantity += e[i] * s[i];
        }
      }
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    if (quantity) {
      TacsScalar Xxi[6], n0[3], X[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);
      *quantity = con->evalDensity(elemIndex, pt, X);
    }

    return 1;
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    if (quantity) {
      TacsScalar Xxi[6], n0[3], X[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);

      // Compute the interpolated displacements
      basis::template interpFields<vars_per_node, 3>(pt, vars, quantity);
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    if (quantity) {
      TacsScalar Xxi[6], n0[3], X[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);
      TacsScalar density = con->evalDensity(elemIndex, pt, X);

      quantity[0] = density * X[0];
      quantity[1] = density * X[1];
      quantity[2] = density * X[2];
    }

    return 3;
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    if (quantity) {
      // Compute X, X,xi and the interpolated normal n0
      TacsScalar X[3], Xxi[6], n0[3], T[9];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      // Compute the transformation at the quadrature point
      transform->computeTransform(Xxi, n0, T);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);
      TacsScalar density = con->evalDensity(elemIndex, pt, X);

      TacsScalar I0[6] = {0.0};

      // Evaluate the self MOI
      TacsScalar moments[3];
      con->evalMassMoments(elemIndex, pt, X, moments);
      I0[0] = I0[3] = moments[2];
      // Compute T*I0*T^{T}
      mat3x3SymmTransform(T, I0, quantity);

      // Use parallel axis theorem to move MOI to origin
      quantity[0] += density * (X[1] * X[1] + X[2] * X[2]);
      quantity[1] += -density * X[0] * X[1];
      quantity[2] += -density * X[0] * X[2];
      quantity[3] += density * (X[0] * X[0] + X[2] * X[2]);
      quantity[4] += -density * X[2] * X[1];
      quantity[5] += density * (X[0] * X[0] + X[1] * X[1]);
    }

    return 6;
  }

  return 0;
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                           TacsScalar scale, int n, double pt[],
                           const TacsScalar Xpts[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[],
                           const TacsScalar dfdq[], int dvLen,
                           TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn);

    // Compute the director rates
    TacsScalar d[dsize], ddot[dsize];
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        vars, dvars, fn, d, ddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                             ety);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9];  // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = 0.0;

    if (quantityType == TACS_FAILURE_INDEX) {
      con->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, e, dvLen, dfdx);
    } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
      TacsScalar s[9];
      con->evalStress(elemIndex, pt, X, e, s);
      con->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);
    }
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    TacsScalar X[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);

    con->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    TacsScalar X[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);

    TacsScalar dfdm = 0.0;

    for (int i = 0; i < 3; i++) {
      dfdm += scale * dfdq[i] * X[i];
    }

    con->addDensityDVSens(elemIndex, dfdm, pt, X, dvLen, dfdx);
  } else if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    TacsScalar Xd[9];
    TacsShellAssembleFrame(Xxi, n0, Xd);

    TacsScalar dfdI0[6] = {0.0};

    // Evaluate the self MOI
    TacsScalar dfdmoments[3] = {0.0};
    mat3x3SymmTransformSens(T, dfdq, dfdI0);
    dfdmoments[2] = scale * (dfdI0[0] + dfdI0[3]);

    con->addMassMomentsDVSens(elemIndex, pt, X, dfdmoments, dvLen, dfdx);

    TacsScalar dfdm = 0.0;

    // Use parallel axis theorem to move MOI to origin
    dfdm += scale * dfdq[0] * (X[1] * X[1] + X[2] * X[2]);
    dfdm -= scale * dfdq[1] * X[0] * X[1];
    dfdm -= scale * dfdq[2] * X[0] * X[2];
    dfdm += scale * dfdq[3] * (X[0] * X[0] + X[2] * X[2]);
    dfdm -= scale * dfdq[4] * X[2] * X[1];
    dfdm += scale * dfdq[5] * (X[0] * X[0] + X[1] * X[1]);

    con->addDensityDVSens(elemIndex, dfdm, pt, X, dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                           TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                           int n, double pt[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], const TacsScalar dfdq[],
                           TacsScalar dfdu[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Derivative of the director field
    TacsScalar dd[dsize];
    memset(dd, 0, 3 * num_nodes * sizeof(TacsScalar));

    // Zero the contributions to the tying strain derivatives
    TacsScalar dety[basis::NUM_TYING_POINTS];
    memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar d[dsize], ddot[dsize], dddot[dsize];
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        vars, dvars, ddvars, fn, d, ddot, dddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                             ety);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9];
    TacsShellComputeDispGrad<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);

    // Evaluate the tying components of the strain
    TacsScalar gty[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9];  // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = 0.0;

    TacsScalar sens[9];
    if (quantityType == TACS_FAILURE_INDEX) {
      // Compute the sensitivity of the failure index w.r.t. the strain
      con->evalFailureStrainSens(elemIndex, pt, X, e, sens);
    } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
      // Compute the sensitivity of the strain energy density w.r.t. the strain
      con->evalStress(elemIndex, pt, X, e, sens);
      for (int i = 0; i < 9; i++) {
        sens[i] *= 2.0;
      }
    }

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(alpha * dfdq[0], sens, u0x, u1x, du0x, du1x, de0ty);

    // Add the contributions to the residual from du0x, du1x and dCt
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT, du0x,
                                                   du1x, dfdu, dd);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn, vars, d, dety, dfdu, dd);

    // Add the contributions to the director field
    director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
        vars, dvars, ddvars, fn, dd, dfdu);
  } else if (quantityType == TACS_ELEMENT_DISPLACEMENT) {
    // Compute the interpolated displacements
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, dfdq, dfdu);
  }
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::getOutputData(
    int elemIndex, ElementType etype, int write_flag, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int ld_data, TacsScalar *data) {
  if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
    // Get the number of nodes associated with the visualization
    int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    // Store information about the transformation and derivatives at each node
    // for the drilling degrees of freedom
    TacsScalar etn[num_nodes];
    TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
    TacsScalar u0xn[9 * num_nodes], Ctn[csize];
    TacsShellComputeDrillStrain<vars_per_node, offset, basis, director, model>(
        transform, Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

    TacsScalar d[dsize], ddot[dsize], dddot[dsize];
    director::template computeDirectorRates<vars_per_node, offset, num_nodes>(
        vars, dvars, ddvars, fn, d, ddot, dddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                             ety);

    // Loop over each quadrature point and add the residual contribution
    for (int index = 0; index < num_vis_nodes; index++) {
      // Get the quadrature weight
      double pt[3];
      basis::getNodePoint(index, pt);

      // Compute X, X,xi and the interpolated normal n0
      TacsScalar X[3], Xxi[6], n0[3], T[9], et;
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFields<1, 1>(pt, etn, &et);

      // Compute the transformation at the quadrature point
      transform->computeTransform(Xxi, n0, T);

      // Evaluate the displacement gradient at the point
      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9];
      TacsShellComputeDispGrad<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);

      // Evaluate the tying components of the strain
      TacsScalar gty[6];  // The symmetric components of the tying strain
      basis::interpTyingStrain(pt, ety, gty);

      // Compute the symmetric parts of the tying strain
      TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

      // Compute the set of strain components
      TacsScalar e[9];  // The components of the strain
      model::evalStrain(u0x, u1x, e0ty, e);
      e[8] = et;

      // Compute the corresponding stresses
      TacsScalar s[9];
      con->evalStress(elemIndex, pt, X, e, s);

      if (write_flag & TACS_OUTPUT_NODES) {
        data[0] = X[0];
        data[1] = X[1];
        data[2] = X[2];
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
        for (int i = 0; i < 9; i++) {
          data[i] = e[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_STRESSES) {
        for (int i = 0; i < 9; i++) {
          data[i] = s[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_EXTRAS) {
        data[0] = con->evalFailure(elemIndex, pt, X, e);
        data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
        data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
        data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
        data += 4;
      }
    }
  }
}

template <int vars_per_node, class basis, class model>
int TacsTestShellTyingStrain(double dh = 1e-7, int test_print_level = 2,
                             double test_fail_atol = 1e-5,
                             double test_fail_rtol = 1e-5) {
  const int size = vars_per_node * basis::NUM_NODES;
  const int usize = 3 * basis::NUM_NODES;
  const int dsize = 3 * basis::NUM_NODES;

  TacsScalar Xpts[3 * basis::NUM_NODES], fn[3 * basis::NUM_NODES];
  TacsGenerateRandomArray(Xpts, 3 * basis::NUM_NODES);
  TacsGenerateRandomArray(fn, 3 * basis::NUM_NODES);

  TacsScalar d[dsize], vars[size];
  TacsGenerateRandomArray(d, dsize);
  TacsGenerateRandomArray(vars, size);

  TacsScalar XdinvT[9];
  TacsGenerateRandomArray(XdinvT, 9);

  TacsScalar de0ty[6], d2e0ty[36];
  TacsGenerateRandomArray(de0ty, 6);
  TacsGenerateRandomArray(d2e0ty, 36);

  double pt[2];
  TacsGenerateRandomArray(pt, 2);

  TacsScalar dety[basis::NUM_TYING_POINTS];
  TacsScalar d2ety[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
  TacsScalar d2etyu[basis::NUM_TYING_POINTS * usize];
  TacsScalar d2etyd[basis::NUM_TYING_POINTS * dsize];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(
      d2ety, 0,
      basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(d2etyu, 0, basis::NUM_TYING_POINTS * usize * sizeof(TacsScalar));
  memset(d2etyd, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));

  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d,
                                                           ety);

  // Evaluate the tying components of the strain
  TacsScalar gty[6];  // The symmetric components of the tying strain
  basis::interpTyingStrain(pt, ety, gty);

  // Compute the symmetric parts of the tying strain
  TacsScalar e0ty[6];  // e0ty = XdinvT^{T}*gty*XdinvT
  mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

  // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
  TacsScalar dgty[6], d2gty[36];
  mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
  mat3x3SymmTransformTransHessian(XdinvT, d2e0ty, d2gty);

  // Evaluate the tying strain
  basis::addInterpTyingStrainTranspose(pt, dgty, dety);
  basis::addInterpTyingStrainHessian(pt, d2gty, d2ety);

  TacsScalar res[size], dd[dsize];
  memset(res, 0, size * sizeof(TacsScalar));
  memset(dd, 0, dsize * sizeof(TacsScalar));

  TacsScalar mat[size * size], d2d[dsize * dsize], d2du[dsize * usize];
  memset(mat, 0, size * size * sizeof(TacsScalar));
  memset(d2d, 0, dsize * dsize * sizeof(TacsScalar));
  memset(d2du, 0, dsize * usize * sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, res, dd);
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
      1.0, Xpts, fn, vars, d, dety, d2ety, d2etyu, d2etyd, mat, d2d, d2du);

  TacsScalar fdmat[size * size], fdd2du[dsize * usize];
  for (int i = 0; i < size; i++) {
    TacsScalar varst[size];
    memcpy(varst, vars, size * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[i] = vars[i] + TacsScalar(0.0, dh);
#else
    varst[i] = vars[i] + dh;
#endif  // TACS_USE_COMPLEX

    // Perturb the variables
    TacsScalar etyt[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, varst, d,
                                                             etyt);

    // Evaluate the tying components of the strain
    TacsScalar gtyt[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, etyt, gtyt);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0tyt[6];
    mat3x3SymmTransformTranspose(XdinvT, gtyt, e0tyt);

    TacsScalar de0tyt[6];
    for (int j = 0; j < 6; j++) {
      de0tyt[j] = de0ty[j];
      for (int k = 0; k < 6; k++) {
        de0tyt[j] += d2e0ty[6 * j + k] * (e0tyt[k] - e0ty[k]);
      }
    }

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgtyt[6];
    mat3x3SymmTransformTransSens(XdinvT, de0tyt, dgtyt);

    TacsScalar detyt[basis::NUM_TYING_POINTS];
    memset(detyt, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgtyt, detyt);

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size * sizeof(TacsScalar));
    memset(ddt, 0, dsize * sizeof(TacsScalar));

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn, varst, d, detyt, rest, ddt);

    for (int j = 0; j < size; j++) {
#ifdef TACS_USE_COMPLEX
      fdmat[size * j + i] = TacsImagPart(rest[j]) / dh;
#else
      fdmat[size * j + i] = (rest[j] - res[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }

    if (i % vars_per_node < 3) {
      int index = 3 * (i / vars_per_node) + i % vars_per_node;
      for (int j = 0; j < dsize; j++) {
#ifdef TACS_USE_COMPLEX
        fdd2du[usize * j + index] = TacsImagPart(ddt[j]) / dh;
#else
        fdd2du[usize * j + index] = (ddt[j] - dd[j]) / dh;
#endif  // TACS_USE_COMPLEX
      }
    }
  }

  int fail = 0;
  double max_err, max_rel;
  int max_err_index, max_rel_index;

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size * size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size * size, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size * size);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the error
  max_err = TacsGetMaxError(d2du, fdd2du, dsize * usize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2du, fdd2du, dsize * usize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. vars and d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2du", d2du, fdd2du, dsize * usize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  TacsScalar fdd2d[dsize * dsize];
  for (int i = 0; i < dsize; i++) {
    TacsScalar dt[size];
    memcpy(dt, d, dsize * sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    dt[i] = d[i] + TacsScalar(0.0, dh);
#else
    dt[i] = d[i] + dh;
#endif  // TACS_USE_COMPLEX

    // Perturb the variables
    TacsScalar etyt[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, dt,
                                                             etyt);

    // Evaluate the tying components of the strain
    TacsScalar gtyt[6];  // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, etyt, gtyt);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0tyt[6];
    mat3x3SymmTransformTranspose(XdinvT, gtyt, e0tyt);

    TacsScalar de0tyt[6];
    for (int j = 0; j < 6; j++) {
      de0tyt[j] = de0ty[j];
      for (int k = 0; k < 6; k++) {
        de0tyt[j] += d2e0ty[6 * j + k] * (e0tyt[k] - e0ty[k]);
      }
    }

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgtyt[6];
    mat3x3SymmTransformTransSens(XdinvT, de0tyt, dgtyt);

    TacsScalar detyt[basis::NUM_TYING_POINTS];
    memset(detyt, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgtyt, detyt);

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size * sizeof(TacsScalar));
    memset(ddt, 0, dsize * sizeof(TacsScalar));

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn, vars, dt, detyt, rest, ddt);

    for (int j = 0; j < dsize; j++) {
#ifdef TACS_USE_COMPLEX
      fdd2d[dsize * j + i] = TacsImagPart(ddt[j]) / dh;
#else
      fdd2d[dsize * j + i] = (ddt[j] - dd[j]) / dh;
#endif  // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2d, fdd2d, dsize * dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d, fdd2d, dsize * dsize, &max_rel_index);

  if (test_print_level > 0) {
    fprintf(stderr, "Testing the second derivative w.r.t. d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n", max_err,
            max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n", max_rel,
            max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1) {
    TacsPrintErrorComponents(stderr, "d2d", d2d, fdd2d, dsize * dsize);
  }
  if (test_print_level) {
    fprintf(stderr, "\n");
  }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif  // TACS_SHELL_ELEMENT_H
