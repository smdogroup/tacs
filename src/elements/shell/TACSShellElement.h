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
#include "TACSShellInplaneElementModel.h"
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
                                             const TacsScalar rotCenter[],
                                             const bool first_order = false) {
    return new TACSShellCentrifugalForce<vars_per_node, quadrature, basis>(
        con, omega, rotCenter, first_order);
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

  void addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                               TacsScalar scale, int n, double pt[],
                               const TacsScalar Xpts[], const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfddetXd,
                               const TacsScalar dfdq[], TacsScalar dfdXpts[]);

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

  void addMatDVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, TacsScalar scale,
                                const TacsScalar psi[], const TacsScalar phi[],
                                const TacsScalar Xpts[],
                                const TacsScalar vars[], int dvLen,
                                TacsScalar dfdx[]);

  void addMatXptSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                 double time, TacsScalar scale,
                                 const TacsScalar psi[],
                                 const TacsScalar phi[],
                                 const TacsScalar Xpts[],
                                 const TacsScalar vars[],
                                 TacsScalar dfdXpts[]);

  void getMatSVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                                double time, const TacsScalar psi[],
                                const TacsScalar phi[], const TacsScalar Xpts[],
                                const TacsScalar vars[], TacsScalar dfdu[]);

  void getOutputData(int elemIndex, ElementType etype, int write_flag,
                     const TacsScalar Xpts[], const TacsScalar vars[],
                     const TacsScalar dvars[], const TacsScalar ddvars[],
                     int ld_data, TacsScalar *data);

  void getAverageStresses(int elemIndex, ElementType etype,
                          const TacsScalar Xpts[], const TacsScalar vars[],
                          const TacsScalar dvars[], const TacsScalar ddvars[],
                          TacsScalar *avgStresses);

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
    if (res) {
      basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, du0dot,
                                                                 res);
    }

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
  // Create dummy zero vector for dvars and ddvars
  TacsScalar zeros[vars_per_node * num_nodes];
  memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

  // Set alpha or gamma based on if this is a stiffness or mass matrix
  if (matType == TACS_STIFFNESS_MATRIX) {
    alpha = 1.0;
  } else if (matType == TACS_MASS_MATRIX) {
    gamma = 1.0;
  } else {  // TACS_GEOMETRIC_STIFFNESS_MATRIX
    // Approximate geometric stiffness using directional derivative of
    // tangential stiffness projected along path of current state vars
    alpha = 1.0;

    // deriv direction
    const TacsScalar *path = vars;

    // Compute the number of quadrature points
    const int nquad = quadrature::getNumQuadraturePoints();

    // Derivative of the director field
    TacsScalar dd[dsize];
    memset(dd, 0, dsize * sizeof(TacsScalar));

    // Second derivatives required for the director
    TacsScalar d2d[dsize * dsize], d2du[usize * dsize];
    TacsScalar d2dd[dsize * dsize], d2dud[usize * dsize];
    TacsScalar d2Tdotd[dsize * dsize], d2Tdotu[usize * dsize];
    memset(d2d, 0, dsize * dsize * sizeof(TacsScalar));
    memset(d2du, 0, usize * dsize * sizeof(TacsScalar));
    memset(d2dd, 0, dsize * dsize * sizeof(TacsScalar));
    memset(d2dud, 0, usize * dsize * sizeof(TacsScalar));
    memset(d2Tdotd, 0, dsize * dsize * sizeof(TacsScalar));
    memset(d2Tdotu, 0, usize * dsize * sizeof(TacsScalar));

    // Zero the contributions to the tying strain derivatives
    TacsScalar dety[basis::NUM_TYING_POINTS];
    TacsScalar d2ety[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
    TacsScalar d2etyu[basis::NUM_TYING_POINTS * usize];
    TacsScalar d2etyd[basis::NUM_TYING_POINTS * dsize];
    TacsScalar dety_d[basis::NUM_TYING_POINTS];
    TacsScalar d2ety_d[basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS];
    TacsScalar d2etyud[basis::NUM_TYING_POINTS * usize];
    TacsScalar d2etydd[basis::NUM_TYING_POINTS * dsize];
    memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(
        d2ety, 0,
        basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(d2etyu, 0, basis::NUM_TYING_POINTS * usize * sizeof(TacsScalar));
    memset(d2etyd, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));
    memset(dety_d, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(
        d2ety_d, 0,
        basis::NUM_TYING_POINTS * basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(d2etyud, 0, basis::NUM_TYING_POINTS * usize * sizeof(TacsScalar));
    memset(d2etydd, 0, basis::NUM_TYING_POINTS * dsize * sizeof(TacsScalar));

    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    // Compute the drill strain penalty at each node
    TacsScalar etn[num_nodes], etnd[num_nodes], detn[num_nodes],
        detnd[num_nodes];
    TacsScalar d2etn[num_nodes * num_nodes];
    memset(detn, 0, num_nodes * sizeof(TacsScalar));
    memset(d2etn, 0, num_nodes * num_nodes * sizeof(TacsScalar));
    memset(detnd, 0, num_nodes * sizeof(TacsScalar));

    // Store information about the transformation and derivatives at each node
    // for the drilling degrees of freedom
    TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
    TacsScalar u0xn[9 * num_nodes], Ctn[csize];
    TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                     model>(transform, Xdn, fn, zeros, path,
                                            XdinvTn, Tn, u0xn, Ctn, etn, etnd);

    TacsScalar d[dsize], ddot[dsize], dddot[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        zeros, zeros, zeros, path, fn, d, ddot, dddot, dd);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS], etyd[basis::NUM_TYING_POINTS];
    if (typeid(model) == typeid(TACSShellLinearModel)) {
      TACSShellNonlinearModel::template computeTyingStrainDeriv<vars_per_node,
                                                                basis>(
          Xpts, fn, zeros, d, path, dd, ety, etyd);
    } else if (typeid(model) == typeid(TACSShellInplaneLinearModel)) {
      TACSShellInplaneNonlinearModel::template computeTyingStrainDeriv<
          vars_per_node, basis>(Xpts, fn, zeros, d, path, dd, ety, etyd);
    } else {
      model::template computeTyingStrainDeriv<vars_per_node, basis>(
          Xpts, fn, zeros, d, path, dd, ety, etyd);
    }

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
          pt, Xpts, zeros, fn, d, Xxi, n0, T, path, dd, XdinvT, XdinvzT, u0x,
          u1x, u0xd, u1xd);
      detXd *= weight;

      // Evaluate the tying components of the strain
      TacsScalar gty[6],
          gtyd[6];  // The symmetric components of the tying strain
      basis::interpTyingStrain(pt, ety, gty);
      basis::interpTyingStrain(pt, etyd, gtyd);

      // Compute the symmetric parts of the tying strain
      TacsScalar e0ty[6], e0tyd[6];  // e0ty = XdinvT^{T}*gty*XdinvT
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
      mat3x3SymmTransformTranspose(XdinvT, gtyd, e0tyd);

      // Compute the set of strain components
      TacsScalar e[9];   // The components of the strain
      TacsScalar ed[9];  // The directional derivative components of the strain
      if (typeid(model) == typeid(TACSShellLinearModel)) {
        TACSShellNonlinearModel::evalStrainDeriv(u0x, u1x, e0ty, u0xd, u1xd,
                                                 e0tyd, e, ed);
      } else if (typeid(model) == typeid(TACSShellInplaneLinearModel)) {
        TACSShellInplaneNonlinearModel::evalStrainDeriv(u0x, u1x, e0ty, u0xd,
                                                        u1xd, e0tyd, e, ed);
      } else {
        model::evalStrainDeriv(u0x, u1x, e0ty, u0xd, u1xd, e0tyd, e, ed);
      }
      e[8] = et;
      ed[8] = etd;

      // Compute the tangent stiffness matrix
      TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
      con->evalTangentStiffness(elemIndex, pt, X, Cs);

      TacsScalar drill;
      const TacsScalar *A, *B, *D, *As;
      TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As,
                                                     &drill);

      // Compute the stress based on the tangent stiffness
      TacsScalar s[9], sd[9];
      TACSShellConstitutive::computeStress(A, B, D, As, drill, e, s);
      TACSShellConstitutive::computeStress(A, B, D, As, drill, ed, sd);

      // Compute the derivative of the product of the stress and strain
      // with respect to u0x, u1x and e0ty
      TacsScalar du0x[9], du1x[9], de0ty[6], du0xd[9], du1xd[9], de0tyd[6];
      if (typeid(model) == typeid(TACSShellLinearModel)) {
        TACSShellNonlinearModel::evalStrainSensDeriv(
            detXd, s, u0x, u1x, sd, u0xd, u1xd, du0x, du1x, de0ty, du0xd, du1xd,
            de0tyd);
      } else if (typeid(model) == typeid(TACSShellInplaneLinearModel)) {
        TACSShellInplaneNonlinearModel::evalStrainSensDeriv(
            detXd, s, u0x, u1x, sd, u0xd, u1xd, du0x, du1x, de0ty, du0xd, du1xd,
            de0tyd);
      } else {
        model::evalStrainSensDeriv(detXd, s, u0x, u1x, sd, u0xd, u1xd, du0x,
                                   du1x, de0ty, du0xd, du1xd, de0tyd);
      }

      TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
      TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
      TacsScalar d2u0xd[81], d2u1xd[81], d2u0xu1xd[81];
      TacsScalar d2e0tyd[36], d2e0tyu0xd[54], d2e0tyu1xd[54];
      if (typeid(model) == typeid(TACSShellLinearModel)) {
        TACSShellNonlinearModel::evalStrainHessianDeriv(
            alpha * detXd, s, Cs, u0x, u1x, e0ty, sd, u0xd, u1xd, e0tyd, d2u0x,
            d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x, d2u0xd, d2u1xd,
            d2u0xu1xd, d2e0tyd, d2e0tyu0xd, d2e0tyu1xd);
      } else if (typeid(model) == typeid(TACSShellInplaneLinearModel)) {
        TACSShellInplaneNonlinearModel::evalStrainHessianDeriv(
            alpha * detXd, s, Cs, u0x, u1x, e0ty, sd, u0xd, u1xd, e0tyd, d2u0x,
            d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x, d2u0xd, d2u1xd,
            d2u0xu1xd, d2e0tyd, d2e0tyu0xd, d2e0tyu1xd);
      } else {
        model::evalStrainHessianDeriv(
            alpha * detXd, s, Cs, u0x, u1x, e0ty, sd, u0xd, u1xd, e0tyd, d2u0x,
            d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x, d2u0xd, d2u1xd,
            d2u0xu1xd, d2e0tyd, d2e0tyu0xd, d2e0tyu1xd);
      }

      // Add the contributions to the residual from du0x and du1x
      TacsScalar det = detXd * s[8];
      basis::template addInterpFieldsTranspose<1, 1>(pt, &det, detn);
      TacsScalar detd = detXd * sd[8];
      basis::template addInterpFieldsTranspose<1, 1>(pt, &detd, detnd);

      // Add the contribution from the drilling stiffness
      TacsScalar d2et = detXd * alpha * Cs[21];
      basis::template addInterpFieldsOuterProduct<1, 1, 1, 1>(pt, &d2et, d2etn);

      // Add the contributions to the residual from du0x, du1x and dCt
      TacsShellAddDispGradHessian<vars_per_node, basis>(
          pt, T, XdinvT, XdinvzT, d2u0x, d2u1x, d2u0xu1x, NULL, d2d, d2du);
      TacsShellAddDispGradHessian<vars_per_node, basis>(
          pt, T, XdinvT, XdinvzT, d2u0xd, d2u1xd, d2u0xu1xd, mat, d2dd, d2dud);

      // Compute the of the tying strain w.r.t. derivative w.r.t. the
      // coefficients
      TacsScalar dgty[6], d2gty[36], dgtyd[6], d2gtyd[36];
      mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
      mat3x3SymmTransformTransHessian(XdinvT, d2e0ty, d2gty);
      mat3x3SymmTransformTransSens(XdinvT, de0tyd, dgtyd);
      mat3x3SymmTransformTransHessian(XdinvT, d2e0tyd, d2gtyd);

      // Evaluate the tying strain
      basis::addInterpTyingStrainTranspose(pt, dgty, dety);
      basis::addInterpTyingStrainHessian(pt, d2gty, d2ety);
      basis::addInterpTyingStrainTranspose(pt, dgtyd, dety_d);
      basis::addInterpTyingStrainHessian(pt, d2gtyd, d2ety_d);

      // Add the coupling between the displacement and tying strain
      TacsShellAddTyingDispCoupling<basis>(pt, T, XdinvT, XdinvzT, d2e0tyu0x,
                                           d2e0tyu1x, d2etyu, d2etyd);
      TacsShellAddTyingDispCoupling<basis>(pt, T, XdinvT, XdinvzT, d2e0tyu0xd,
                                           d2e0tyu1xd, d2etyud, d2etydd);
    }

    // Add the second order terms from the tying strain
    if (typeid(model) == typeid(TACSShellLinearModel)) {
      TACSShellNonlinearModel::template addComputeTyingStrainHessianDeriv<
          vars_per_node, basis>(alpha, Xpts, fn, zeros, d, dety, d2ety, d2etyu,
                                d2etyd, path, dd, dety_d, d2ety_d, d2etyud,
                                d2etydd, NULL, d2d, d2du, mat, d2dd, d2dud);
    } else if (typeid(model) == typeid(TACSShellInplaneLinearModel)) {
      TACSShellInplaneNonlinearModel::
          template addComputeTyingStrainHessianDeriv<vars_per_node, basis>(
              alpha, Xpts, fn, zeros, d, dety, d2ety, d2etyu, d2etyd, path, dd,
              dety_d, d2ety_d, d2etyud, d2etydd, NULL, d2d, d2du, mat, d2dd,
              d2dud);
    } else {
      model::template addComputeTyingStrainHessianDeriv<vars_per_node, basis>(
          alpha, Xpts, fn, zeros, d, dety, d2ety, d2etyu, d2etyd, path, dd,
          dety_d, d2ety_d, d2etyud, d2etydd, NULL, d2d, d2du, mat, d2dd, d2dud);
    }

    // Add the contributions to the stiffness matrix
    director::template addDirectorJacobian<vars_per_node, offset, num_nodes>(
        alpha, beta, gamma, path, zeros, zeros, fn, dd, d2Tdotd, d2Tdotu, d2dd,
        d2dud, NULL, mat);

    return;
  }
  // Add appropriate Jacobian to matrix
  addJacobian(elemIndex, time, alpha, beta, gamma, Xpts, vars, zeros, zeros,
              NULL, mat);
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
void TACSShellElement<quadrature, basis, director, model>::
    addAdjResXptProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], TacsScalar fXptSens[]) {
  // The dynamics/rotary-inertia closure below (dd_xpt_ddot's substitution of
  // ddvars for vars into the single-seed addDirectorRefNormalSens overload)
  // is exact only for TACSLinearizedRotation, whose dddot = crossProduct(
  // ddvars-rotation, fn) has the same bilinear structure as d = crossProduct(
  // vars-rotation, fn). For TACSQuadraticRotation/TACSQuaternionRotation,
  // dddot genuinely couples q/qdot/qddot together (see the closure's own
  // comment below) and this substitution is only an approximation - fine
  // under the real-mode rtol=1e-2 gate, but not exact under complex-step
  // (verified: forwarding here is what keeps Complex-Ubuntu/Complex-MacOS CI
  // green for *NonlinearShellModRot elements). Forward to the base-class
  // FD/CS implementation for those director classes, mirroring
  // addMatDVSensInnerProduct's/getMatSVSensInnerProduct's identical guard; a
  // fully exact per-director closure remains a fast-follow candidate.
  if (typeid(director) != typeid(TACSLinearizedRotation)) {
    TACSElement::addAdjResXptProduct(elemIndex, time, scale, psi, Xpts, vars,
                                     dvars, ddvars, fXptSens);
    return;
  }

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Accumulator for the nodal frame-normal-direction sensitivity, closed at
  // the end via TacsShellAddNodeNormalsSens
  TacsScalar dfn[3 * num_nodes];
  memset(dfn, 0, 3 * num_nodes * sizeof(TacsScalar));

  // Store information about the transformation and derivatives at each node
  // for the drilling degrees of freedom
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

  // Xpts-direction accumulators for the director field: dd_xpt is the
  // primal-direction seed, dd_xpt_psi is the psi-direction seed. Both are
  // needed since the residual is bilinear in the director field (mirrors
  // TACSBeamElement.h's dd1/dd1psi pairing, TACSBeamElement.h:930-935).
  TacsScalar dd_xpt[dsize], dd_xpt_psi[dsize];
  memset(dd_xpt, 0, dsize * sizeof(TacsScalar));
  memset(dd_xpt_psi, 0, dsize * sizeof(TacsScalar));

  // Xpts-direction accumulator for the *acceleration*-direction director
  // field (dddot = crossProduct(ddvars-rotation, fn)), needed by the
  // dynamics/rotary-inertia term below. This is a THIRD, distinct seed
  // buffer (not dd_xpt/dd_xpt_psi) because dddot is parametrized by ddvars,
  // not vars or psi, so it needs its own addDirectorRefNormalSens closure
  // (the single-seed 4-argument overload, called with ddvars substituted for
  // vars -- exact, not an approximation, since d = crossProduct(q,t) is
  // linear in q for any q).
  TacsScalar dd_xpt_ddot[dsize];
  memset(dd_xpt_ddot, 0, dsize * sizeof(TacsScalar));

  // Reverse-mode adjoint seed accumulators for the tying strain (primal and
  // psi-direction), closed after the quadrature loop via
  // model::addTyingStrainXptSens
  TacsScalar dety[basis::NUM_TYING_POINTS], dety_psi[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
  memset(dety_psi, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

  // Reverse-mode adjoint seed accumulators for the drill strain (primal and
  // psi-direction), closed after the quadrature loop via
  // TacsShellAddDrillStrainXptSens / TacsShellAddDrillStrainXptSensDeriv
  TacsScalar detn_c2[num_nodes], detn_c1[num_nodes];
  memset(detn_c2, 0, num_nodes * sizeof(TacsScalar));
  memset(detn_c1, 0, num_nodes * sizeof(TacsScalar));

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
    TacsScalar gty[6], gtyd[6];  // The symmetric components of the tying
                                 // strain (primal and psi-direction)
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

    // Compute the stresses corresponding to the primal and psi-direction
    // strains (stress is linear in strain for a fixed tangent stiffness)
    TacsScalar s[9], sd[9];
    con->evalStress(elemIndex, pt, X, e, s);
    con->evalStress(elemIndex, pt, X, ed, sd);

    // Reverse-mode seeds on the drill strain at this quadrature point
    // (both chains, mirroring the getMatType Deriv path's
    // "det = detXd*s[8]; addInterpFieldsTranspose(pt,&det,detn);" pattern,
    // TACSShellElement.h:856-857), accumulated into per-node buffers and
    // closed post-loop via TacsShellAddDrillStrainXptSens (chain 2) /
    // TacsShellAddDrillStrainXptSensDeriv (chain 1)
    TacsScalar det_c2 = scale * detXd * sd[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det_c2, detn_c2);
    TacsScalar det_c1 = scale * detXd * s[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det_c1, detn_c1);

    // The seed on detXd from the strain-energy-like term, mirroring beam's
    // detXd.valued = scale*(e[0]*spsi[0] + ...) at TACSBeamElement.h:1084-1085
    TacsScalar ddetXd_total =
        scale * (e[0] * sd[0] + e[1] * sd[1] + e[2] * sd[2] + e[3] * sd[3] +
                 e[4] * sd[4] + e[5] * sd[5] + e[6] * sd[6] + e[7] * sd[7] +
                 e[8] * sd[8]);

    // Evaluate the second time derivatives and add the dynamics contribution
    // to the seed on detXd (detXd's own Xpts-dependence; the director-field
    // fields d0ddot/dd0ddot used below have a SEPARATE Xpts-dependence via
    // fn, handled just after)
    TacsScalar u0ddot[3], d0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, ddvars, u0ddot);
    basis::template interpFields<3, 3>(pt, dddot, d0ddot);

    TacsScalar du0ddot[3], dd0ddot[3];
    basis::template interpFields<vars_per_node, 3>(pt, psi, du0ddot);
    basis::template interpFields<3, 3>(pt, dd, dd0ddot);

    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    ddetXd_total +=
        scale * (moments[0] * vec3Dot(u0ddot, du0ddot) +
                 moments[1] * (vec3Dot(u0ddot, dd0ddot) +
                               vec3Dot(du0ddot, d0ddot)) +
                 moments[2] * vec3Dot(d0ddot, dd0ddot));

    // The rotary-inertia moments (moments[1], moments[2]) couple through the
    // director-derivative fields d0ddot/dd0ddot, which themselves depend on
    // fn (hence Xpts) via crossProduct(rotation-dof, fn) -- this is a SEPARATE
    // Xpts-dependence from detXd's (folded into ddetXd_total above), and was
    // previously missing entirely. moments[0] (translational mass) has no
    // such dependence (u0ddot/du0ddot are plain translational DOFs, not
    // director-derived), so no seed is needed for that term.
    TacsScalar seed_dd0ddot[3], seed_d0ddot[3];
    for (int i = 0; i < 3; i++) {
      seed_dd0ddot[i] = scale * detXd *
                       (moments[1] * u0ddot[i] + moments[2] * d0ddot[i]);
      seed_d0ddot[i] = scale * detXd *
                      (moments[1] * du0ddot[i] + moments[2] * dd0ddot[i]);
    }
    // dd0ddot = interpFields(dd) is the psi-direction rate field -- the same
    // field already closed via dd_xpt_psi/addDirectorRefNormalSens's ddpsi
    // slot, so its seed accumulates into that same buffer.
    basis::template addInterpFieldsTranspose<3, 3>(pt, seed_dd0ddot,
                                                   dd_xpt_psi);
    // d0ddot = interpFields(dddot) is the ddvars-direction acceleration
    // field -- distinct from both the vars- and psi-direction fields, so it
    // needs its own buffer, closed post-loop below.
    basis::template addInterpFieldsTranspose<3, 3>(pt, seed_d0ddot,
                                                   dd_xpt_ddot);

    TacsScalar dXxi[6], dn0[3], dT[9];
    memset(dXxi, 0, sizeof(dXxi));
    memset(dn0, 0, sizeof(dn0));
    memset(dT, 0, sizeof(dT));

    // ---- Chain 2 (primal-direction adjoint, seed = sd) ----
    // NOTE: TacsShellAddDispGradSens is deliberately NOT called here -- its
    // "dd" output is mathematically identical to TacsShellComputeDispGradXptSens's
    // own "dd" output (both reduce to the same T*(du0x*XdinvT^T+du1x*XdinvzT^T)
    // chain applied to the same seed), so calling both and accumulating into
    // the same buffer would double-count this contribution (verified
    // numerically in isolation). TacsShellComputeDispGradXptSens alone
    // supplies the complete dd contribution from this seed.
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(scale * detXd, sd, u0x, u1x, du0x, du1x, de0ty);

    TacsShellComputeDispGradXptSens<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, du0x, du1x,
        ddetXd_total * weight, dXxi, dn0, dfn, dT, dd_xpt);

    // Tying-strain (primal-direction) contribution: accumulate the seed on
    // the tying-point strain, closed after the loop
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);

    // ---- Chain 1 (psi-direction adjoint, seed = s) ----
    // u0xd/u1xd (psi-direction) were built from psi/dd (not vars/d) by
    // TacsShellComputeDispGradDeriv, so the reverse sweep's internal
    // recompute of u0xi/d0/d0xi must also use psi/dd here to match.
    // PLUS the u0x/u1x Hessian-coupling term: for a strain measure that is
    // itself nonlinear in (u0x,u1x) (TACSShellNonlinearModel/
    // TACSShellInplaneNonlinearModel), the Jacobian J(u0x,u1x) used above
    // depends on the base point, so ed = J(u0x,u1x)*(u0xd,u1xd,e0tyd)'s own
    // Xpts-dependence has a second contribution beyond chain 1 above: the
    // base point's Xpts-dependence feeding through dJ/d(u0x,u1x). This is
    // exactly what evalStrainSensDeriv's second (u0xd,u1xd,e0tyd)-slot
    // output gives when called with a zero seed-direction (dfded=0): a
    // Hessian . (u0xd,u1xd) contraction against the seed s. It is
    // identically zero for the two linear model classes (evalStrainSens
    // there ignores its own u0x/u1x arguments entirely -- verified
    // numerically), so no model-class branch is needed; it naturally
    // vanishes.
    TacsScalar du0xd[9], du1xd[9], de0tyd[6];
    TacsScalar dfded_zero[9];
    memset(dfded_zero, 0, sizeof(dfded_zero));
    TacsScalar du0x_hess[9], du1x_hess[9], de0ty_hess[6];
    model::evalStrainSensDeriv(scale * detXd, s, u0x, u1x, dfded_zero, u0xd,
                               u1xd, du0xd, du1xd, de0tyd, du0x_hess,
                               du1x_hess, de0ty_hess);

    TacsShellComputeDispGradXptSens<vars_per_node, basis>(
        pt, Xpts, psi, fn, dd, Xxi, n0, T, XdinvT, XdinvzT, du0xd, du1xd, 0.0,
        dXxi, dn0, dfn, dT, dd_xpt_psi);

    TacsScalar dgtyd[6];
    mat3x3SymmTransformTransSens(XdinvT, de0tyd, dgtyd);
    basis::addInterpTyingStrainTranspose(pt, dgtyd, dety_psi);

    // The Hessian-coupling term routes through the PRIMAL (vars,d) chain
    // (it is the base point's own Xpts-dependence), accumulating into
    // dd_xpt alongside chain 2 -- no separate detXd seed here (chain 2
    // above already carries the complete one). de0ty_hess is identically
    // zero (e0ty enters the strain linearly, no Hessian coupling there), so
    // no separate tying-transform-Sens/dgty step is needed for this term.
    TacsShellComputeDispGradXptSens<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, du0x_hess,
        du1x_hess, 0.0, dXxi, dn0, dfn, dT, dd_xpt);

    // ---- XdinvT-direction contribution from the tying-strain transform
    // (e0ty = XdinvT^{T}*gty*XdinvT, e0tyd = XdinvT^{T}*gtyd*XdinvT) ----
    // Not captured by TacsShellComputeDispGradXptSens, whose contract is
    // limited to the u0x/u1x/detXd chain; hand-derived here mirroring the
    // mat3x3MatMult adjoint rules already used above. NOTE: the seed's
    // off-diagonal (shear) entries must be halved when expanded from the
    // compact 6-component symmetric storage to a full 3x3 matrix -- the
    // compact form counts each off-diagonal term once, but expanding it
    // symmetrically into both (i,j) and (j,i) slots of the full matrix and
    // then using the plain (non-symmetric-aware) mat3x3MatMult adjoint rule
    // would otherwise double-count it (verified numerically in isolation).
    {
      TacsScalar gty_full[9], gtyd_full[9], de0ty_full[9], de0tyd_full[9];
      gty_full[0] = gty[0];
      gty_full[1] = gty[1];
      gty_full[2] = gty[2];
      gty_full[3] = gty[1];
      gty_full[4] = gty[3];
      gty_full[5] = gty[4];
      gty_full[6] = gty[2];
      gty_full[7] = gty[4];
      gty_full[8] = gty[5];

      gtyd_full[0] = gtyd[0];
      gtyd_full[1] = gtyd[1];
      gtyd_full[2] = gtyd[2];
      gtyd_full[3] = gtyd[1];
      gtyd_full[4] = gtyd[3];
      gtyd_full[5] = gtyd[4];
      gtyd_full[6] = gtyd[2];
      gtyd_full[7] = gtyd[4];
      gtyd_full[8] = gtyd[5];

      de0ty_full[0] = de0ty[0];
      de0ty_full[1] = 0.5 * de0ty[1];
      de0ty_full[2] = 0.5 * de0ty[2];
      de0ty_full[3] = 0.5 * de0ty[1];
      de0ty_full[4] = de0ty[3];
      de0ty_full[5] = 0.5 * de0ty[4];
      de0ty_full[6] = 0.5 * de0ty[2];
      de0ty_full[7] = 0.5 * de0ty[4];
      de0ty_full[8] = de0ty[5];

      de0tyd_full[0] = de0tyd[0];
      de0tyd_full[1] = 0.5 * de0tyd[1];
      de0tyd_full[2] = 0.5 * de0tyd[2];
      de0tyd_full[3] = 0.5 * de0tyd[1];
      de0tyd_full[4] = de0tyd[3];
      de0tyd_full[5] = 0.5 * de0tyd[4];
      de0tyd_full[6] = 0.5 * de0tyd[2];
      de0tyd_full[7] = 0.5 * de0tyd[4];
      de0tyd_full[8] = de0tyd[5];

      TacsScalar dXdinvT[9];
      memset(dXdinvT, 0, sizeof(dXdinvT));
      TacsScalar W[9], dW[9], tmp[9];

      // chain 2: e0ty = XdinvT^{T}*W, W = gty_full*XdinvT
      mat3x3MatMult(gty_full, XdinvT, W);
      mat3x3MatMult(W, de0ty_full, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];
      mat3x3MatMult(XdinvT, de0ty_full, dW);
      mat3x3MatMult(gty_full, dW, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

      // chain 1: e0tyd = XdinvT^{T}*Wd, Wd = gtyd_full*XdinvT
      mat3x3MatMult(gtyd_full, XdinvT, W);
      mat3x3MatMult(W, de0tyd_full, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];
      mat3x3MatMult(XdinvT, de0tyd_full, dW);
      mat3x3MatMult(gtyd_full, dW, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

      // Propagate dXdinvT back through XdinvT = Xdinv*T and
      // Xdinv = inv3x3(Xd), Xd = assembleFrame(Xxi, n0)
      TacsScalar Xd[9], Xdinv[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      inv3x3(Xd, Xdinv);

      TacsScalar dXdinv[9];
      mat3x3MatTransMult(dXdinvT, T, dXdinv);
      mat3x3TransMatMult(Xdinv, dXdinvT, tmp);
      for (int i = 0; i < 9; i++) dT[i] += tmp[i];

      TacsScalar dXd[9];
      inv3x3Sens(Xdinv, dXdinv, dXd);

      dXxi[0] += dXd[0];
      dXxi[1] += dXd[1];
      dXxi[2] += dXd[3];
      dXxi[3] += dXd[4];
      dXxi[4] += dXd[6];
      dXxi[5] += dXd[7];

      dn0[0] += dXd[2];
      dn0[1] += dXd[5];
      dn0[2] += dXd[8];
    }

    // Fold the T-direction seed back onto Xxi/n0
    transform->addTransformSens(Xxi, n0, dT, dXxi, dn0);

    // Close the loop from Xxi/n0 onto the element-level Xpts output and the
    // nodal frame-normal field
    basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, fXptSens);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn);
  }

  // Add the sensitivity contributions from the drilling strain (primal and
  // psi-direction). The psi-direction chain needs the joint Cd =
  // computeRotationMatDeriv(vars,psi) [C(q) is nonlinear in q, unlike the
  // linear u0x/u1x/tying-strain interpolations, so a bare vars->psi
  // substitution does not work here] -- handled by the
  // TacsShellAddDrillStrainXptSensDeriv sibling helper.
  TacsShellAddDrillStrainXptSens<vars_per_node, offset, basis, director,
                                 model>(transform, Xdn, fn, vars, XdinvTn, Tn,
                                        u0xn, Ctn, detn_c2, fXptSens, dfn);
  TacsShellAddDrillStrainXptSensDeriv<vars_per_node, offset, basis, director,
                                      model>(transform, Xdn, fn, vars, psi,
                                             XdinvTn, Tn, u0xn, Ctn, detn_c1,
                                             fXptSens, dfn);

  // Add the sensitivity contributions from the tying strain (primal and
  // psi-direction). The psi-direction call substitutes psi/dd for vars/d:
  // dety_psi is a seed on etyd (= computeTyingStrainDeriv's psi-direction
  // output, itself built by replaying computeTyingStrain's own formula with
  // psi/dd in place of vars/d), so the reverse sweep must differentiate
  // THAT same substituted map to match (verified numerically in isolation;
  // mirrors the identical substitution TacsShellComputeDispGradXptSens
  // needs above).
  model::template addTyingStrainXptSens<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, fXptSens, dfn, dd_xpt);
  model::template addTyingStrainXptSens<vars_per_node, basis>(
      Xpts, fn, psi, dd, dety_psi, fXptSens, dfn, dd_xpt_psi);

  // Correct the psi-direction call above for model classes where
  // computeTyingStrain has a genuine (director-field, Uxi) cross term (the
  // two nonlinear model classes): the substituted call's single "dd" output
  // silently drops the primal-direction contribution and mis-attributes part
  // of the psi-direction one. This is a no-op for the two linear model
  // classes, where the substitution above is already exact (see
  // TACSShellLinearModel::addTyingStrainXptSensDeriv).
  model::template addTyingStrainXptSensDeriv<vars_per_node, basis>(
      vars, psi, dety_psi, dd_xpt, dd_xpt_psi);

  // Add the contributions from the derivative of the director
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              num_nodes>(
      vars, psi, fn, dd_xpt, dd_xpt_psi, dfn);

  // Close the dynamics/rotary-inertia term's acceleration-direction seed
  // (dd_xpt_ddot, accumulated above). This point is only ever reached for
  // TACSLinearizedRotation (the typeid guard at the top of this function
  // forwards every other director class to the base-class FD/CS
  // implementation, precisely because this closure is NOT exact for them --
  // see that guard's comment for why). For TACSLinearizedRotation,
  // dddot = crossProduct(ddvars-rotation, fn) has exactly the same bilinear
  // structure as d = crossProduct(vars-rotation, fn), so substituting
  // ddvars for vars into the single-seed overload is EXACT (verified: this
  // combination reaches ~1e-8-1e-10 relative error under both the real-mode
  // pytest gate and complex step).
  director::template addDirectorRefNormalSens<vars_per_node, offset,
                                              num_nodes>(ddvars, fn,
                                                         dd_xpt_ddot, dfn);

  // Add the contributions from the node normals
  TacsShellAddNodeNormalsSens<basis>(Xpts, dfn, fXptSens);
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    addPointQuantityXptSens(int elemIndex, int quantityType, double time,
                            TacsScalar scale, int n, double pt[],
                            const TacsScalar Xpts[], const TacsScalar vars[],
                            const TacsScalar dvars[],
                            const TacsScalar ddvars[],
                            const TacsScalar dfddetXd,
                            const TacsScalar dfdq[], TacsScalar dfdXpts[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY) {
    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn);

    // Accumulator for the nodal frame-normal-direction sensitivity, closed
    // at the end via TacsShellAddNodeNormalsSens
    TacsScalar dfn[3 * num_nodes];
    memset(dfn, 0, 3 * num_nodes * sizeof(TacsScalar));

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

    // Compute the sensitivity of the quantity w.r.t. the strain (esens
    // mirrors beam's esens construction, TACSBeamElement.h:1971-1981,
    // generalized to shell's 9-component strain, with drill e[8] fixed at
    // 0 exactly as evalPointQuantity's own strain-based branch does)
    TacsScalar esens[9];
    if (quantityType == TACS_FAILURE_INDEX) {
      con->evalFailureStrainSens(elemIndex, pt, X, e, esens);
    } else {  // quantityType == TACS_STRAIN_ENERGY_DENSITY
      con->evalStress(elemIndex, pt, X, e, esens);
      for (int i = 0; i < 9; i++) {
        esens[i] *= 2.0;
      }
    }

    // The detXd-direction seed is always included, regardless of branch
    TacsScalar ddetXd_total = scale * dfddetXd;

    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(scale * dfdq[0], esens, u0x, u1x, du0x, du1x,
                          de0ty);

    TacsScalar dXxi[6], dn0[3], dT[9], dd[dsize];
    memset(dXxi, 0, sizeof(dXxi));
    memset(dn0, 0, sizeof(dn0));
    memset(dT, 0, sizeof(dT));
    memset(dd, 0, dsize * sizeof(TacsScalar));

    TacsShellComputeDispGradXptSens<vars_per_node, basis>(
        pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, du0x, du1x,
        ddetXd_total, dXxi, dn0, dfn, dT, dd);

    // Tying-strain contribution: accumulate the seed on the tying-point
    // strain, then close it into the Xpts/fn/director-field sensitivities
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
    TacsScalar dety[basis::NUM_TYING_POINTS];
    memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);
    model::template addTyingStrainXptSens<vars_per_node, basis>(
        Xpts, fn, vars, d, dety, dfdXpts, dfn, dd);

    // XdinvT-direction contribution from the tying-strain transform itself
    // (e0ty = XdinvT^{T}*gty*XdinvT) -- NOT captured by
    // TacsShellComputeDispGradXptSens, whose contract is limited to the
    // u0x/u1x/detXd chain; hand-derived here mirroring the mat3x3MatMult
    // adjoint rules addAdjResXptProduct already uses for the identical
    // e0ty=XdinvT^{T}*gty*XdinvT relationship ("XdinvT-direction
    // contribution from the tying-strain transform" block above, single-seed
    // specialization of that two-chain derivation). NOTE: the seed's
    // off-diagonal (shear) entries must be halved when expanded from the
    // compact 6-component symmetric storage to a full 3x3 matrix -- the
    // compact form counts each off-diagonal term once, but expanding it
    // symmetrically into both (i,j) and (j,i) slots of the full matrix and
    // then using the plain (non-symmetric-aware) mat3x3MatMult adjoint rule
    // would otherwise double-count it.
    {
      TacsScalar gty_full[9], de0ty_full[9];
      gty_full[0] = gty[0];
      gty_full[1] = gty[1];
      gty_full[2] = gty[2];
      gty_full[3] = gty[1];
      gty_full[4] = gty[3];
      gty_full[5] = gty[4];
      gty_full[6] = gty[2];
      gty_full[7] = gty[4];
      gty_full[8] = gty[5];

      de0ty_full[0] = de0ty[0];
      de0ty_full[1] = 0.5 * de0ty[1];
      de0ty_full[2] = 0.5 * de0ty[2];
      de0ty_full[3] = 0.5 * de0ty[1];
      de0ty_full[4] = de0ty[3];
      de0ty_full[5] = 0.5 * de0ty[4];
      de0ty_full[6] = 0.5 * de0ty[2];
      de0ty_full[7] = 0.5 * de0ty[4];
      de0ty_full[8] = de0ty[5];

      TacsScalar dXdinvT[9];
      memset(dXdinvT, 0, sizeof(dXdinvT));
      TacsScalar W[9], dW[9], tmp[9];

      // e0ty = XdinvT^{T}*W, W = gty_full*XdinvT
      mat3x3MatMult(gty_full, XdinvT, W);
      mat3x3MatMult(W, de0ty_full, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];
      mat3x3MatMult(XdinvT, de0ty_full, dW);
      mat3x3MatMult(gty_full, dW, tmp);
      for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

      // Propagate dXdinvT back through XdinvT = Xdinv*T and
      // Xdinv = inv3x3(Xd), Xd = assembleFrame(Xxi, n0)
      TacsScalar Xd[9], Xdinv[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      inv3x3(Xd, Xdinv);

      TacsScalar dXdinv[9];
      mat3x3MatTransMult(dXdinvT, T, dXdinv);
      mat3x3TransMatMult(Xdinv, dXdinvT, tmp);
      for (int i = 0; i < 9; i++) dT[i] += tmp[i];

      TacsScalar dXd[9];
      inv3x3Sens(Xdinv, dXdinv, dXd);

      dXxi[0] += dXd[0];
      dXxi[1] += dXd[1];
      dXxi[2] += dXd[3];
      dXxi[3] += dXd[4];
      dXxi[4] += dXd[6];
      dXxi[5] += dXd[7];

      dn0[0] += dXd[2];
      dn0[1] += dXd[5];
      dn0[2] += dXd[8];
    }

    // Add the contributions from the derivative of the director (the
    // single-seed 4-argument overload -- exact, since only one adjoint
    // direction, "dd", is closed here, unlike addAdjResXptProduct's
    // bilinear psi-direction pairing)
    director::template addDirectorRefNormalSens<vars_per_node, offset,
                                                num_nodes>(vars, fn, dd, dfn);

    // transform->addTransformSens is needed here since dT is generally
    // nonzero for this branch (T enters u0x/u1x's own construction)
    transform->addTransformSens(Xxi, n0, dT, dXxi, dn0);

    basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn);

    // Add the contributions from the node normals (final step, every
    // branch touches fn/n0 in some form, even if only through detXd)
    TacsShellAddNodeNormalsSens<basis>(Xpts, dfn, dfdXpts);
    return;
  }

  if (quantityType == TACS_ELEMENT_MOMENT_OF_INERTIA) {
    // Explicit punt, mirroring beam's own punt exactly (TACSBeamElement.h:
    // 1993-1997): this quantity's forward formula additionally depends on T
    // (via mat3x3SymmTransform) on top of everything DENSITY_MOMENT depends
    // on, plus the parallel-axis dXcg term -- deriving this analytically is
    // disproportionate effort for one quantity type when a correct,
    // precedented fallback exists.
    TACSElement::addPointQuantityXptSens(elemIndex, quantityType, time, scale,
                                         n, pt, Xpts, vars, dvars, ddvars,
                                         dfddetXd, dfdq, dfdXpts);
    return;
  }

  if (quantityType != TACS_ELEMENT_DENSITY &&
      quantityType != TACS_ELEMENT_DISPLACEMENT &&
      quantityType != TACS_ELEMENT_DENSITY_MOMENT &&
      quantityType != TACS_ELEMENT_ENCLOSED_VOLUME) {
    // Unrecognized quantity type: evalPointQuantity returns 0 (no value
    // defined) for anything not in the 7-row coverage table, so there is
    // nothing to differentiate here -- both this method and the base
    // class's FD loop would compute nothing, so a plain no-op (not a
    // forward to base) is the correct, non-wasteful match.
    return;
  }

  // Geometry-only branches: TACS_ELEMENT_DENSITY/TACS_ELEMENT_DISPLACEMENT
  // depend on Xpts only through the shared detXd-direction term (step 4
  // below); TACS_ELEMENT_DENSITY_MOMENT/TACS_ELEMENT_ENCLOSED_VOLUME also
  // have a direct X/n0 seed (mirrors beam's DENSITY_MOMENT branch setting
  // X0.xd/n1.xd/n2.xd directly, TACSBeamElement.h:1988-1991).
  TacsScalar fn[3 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn);

  TacsScalar dfn[3 * num_nodes];
  memset(dfn, 0, 3 * num_nodes * sizeof(TacsScalar));

  TacsScalar X[3], Xxi[6], n0[3];
  basis::template interpFields<3, 3>(pt, Xpts, X);
  basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
  basis::template interpFields<3, 3>(pt, fn, n0);

  if (quantityType == TACS_ELEMENT_DENSITY_MOMENT) {
    // quantity[0..2] = density*X + moments[1]*n0
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);
    TacsScalar density = moments[0];

    TacsScalar dX[3], dn0_direct[3];
    for (int i = 0; i < 3; i++) {
      dX[i] = scale * dfdq[i] * density;
      dn0_direct[i] = scale * dfdq[i] * moments[1];
    }
    basis::template addInterpFieldsTranspose<3, 3>(pt, dX, dfdXpts);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dn0_direct, dfn);
  } else if (quantityType == TACS_ELEMENT_ENCLOSED_VOLUME) {
    // quantity[0] = (X.n0)/3
    TacsScalar dX[3], dn0_direct[3];
    for (int i = 0; i < 3; i++) {
      dX[i] = scale * dfdq[0] * n0[i] / 3.0;
      dn0_direct[i] = scale * dfdq[0] * X[i] / 3.0;
    }
    basis::template addInterpFieldsTranspose<3, 3>(pt, dX, dfdXpts);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dn0_direct, dfn);
  }
  // TACS_ELEMENT_DENSITY/TACS_ELEMENT_DISPLACEMENT contribute no direct
  // seed on X/n0/quantity beyond the shared detXd-direction term below.

  // Shared detXd-direction term: the ONLY Xpts-dependence for
  // TACS_ELEMENT_DENSITY/TACS_ELEMENT_DISPLACEMENT, and an additional term
  // for the other two branches above. No du0x/du1x seed applies here, so
  // this takes the lighter-weight direct path through
  // Xd = assembleFrame(Xxi, n0); detXd = det3x3(Xd) -- the same chain
  // TacsShellComputeDispGradXptSens's own detXd-direction step reduces to
  // when its du0x/du1x seeds are zero (SPEC's documented "implementer's
  // choice": avoids the wasted work of computing the director field/tying
  // strain forward recompute that TacsShellComputeDispGradXptSens's
  // contract otherwise requires as inputs).
  TacsScalar Xd[9], Xdinv[9];
  TacsShellAssembleFrame(Xxi, n0, Xd);
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  TacsScalar ddetXd_total = scale * dfddetXd;
  TacsScalar dXd[9];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dXd[3 * i + j] = ddetXd_total * detXd * Xdinv[3 * j + i];
    }
  }

  TacsScalar dXxi[6], dn0[3];
  dXxi[0] = dXd[0];
  dXxi[1] = dXd[1];
  dXxi[2] = dXd[3];
  dXxi[3] = dXd[4];
  dXxi[4] = dXd[6];
  dXxi[5] = dXd[7];

  dn0[0] = dXd[2];
  dn0[1] = dXd[5];
  dn0[2] = dXd[8];

  // transform->addTransformSens is only needed for branches whose dT is
  // nonzero (the strain-based branches above); safe to call unconditionally
  // with a zeroed dT here since none of these geometry-only branches touch T
  TacsScalar dT[9];
  memset(dT, 0, sizeof(dT));
  transform->addTransformSens(Xxi, n0, dT, dXxi, dn0);

  basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
  basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn);

  // Add the contributions from the node normals (final step, every branch
  // touches fn/n0 in some form, even if only through detXd)
  TacsShellAddNodeNormalsSens<basis>(Xpts, dfn, dfdXpts);
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
      TacsScalar moments[3];
      con->evalMassMoments(elemIndex, pt, X, moments);
      TacsScalar density = moments[0];

      quantity[0] = density * X[0] + moments[1] * n0[0];
      quantity[1] = density * X[1] + moments[1] * n0[1];
      quantity[2] = density * X[2] + moments[1] * n0[2];
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

      TacsScalar I0[6] = {0.0};

      // Evaluate the self MOI
      TacsScalar moments[3];
      con->evalMassMoments(elemIndex, pt, X, moments);
      TacsScalar density = moments[0];
      I0[0] = I0[3] = moments[2] - moments[1] * moments[1] / density;
      // Compute T*I0*T^{T}
      mat3x3SymmTransform(T, I0, quantity);
      TacsScalar dXcg[3];
      for (int i = 0; i < 3; i++) {
        dXcg[i] = X[i] + moments[1] / density * n0[i];
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

  else if (quantityType == TACS_ELEMENT_ENCLOSED_VOLUME) {
    if (quantity) {
      // Compute X, X,xi and the interpolated normal n0
      TacsScalar Xxi[6], n0[3], X[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);

      // Compute 1/3*int[x * n]dA
      // This can be shown to equivalent to the volume through Gauss' Theorem
      quantity[0] = (X[0] * n0[0] + X[1] * n0[1] + X[2] * n0[2]) / 3.0;
    }

    return 1;
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
    // Compute the node normal directions
    TacsScalar fn[3 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn);

    // Compute X and the interpolated normal n0
    TacsScalar X[3], n0[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFields<3, 3>(pt, fn, n0);

    TacsScalar dfdmoments[3] = {0.0};

    for (int i = 0; i < 3; i++) {
      dfdmoments[0] += scale * dfdq[i] * X[i];
      dfdmoments[1] += scale * dfdq[i] * n0[i];
    }

    con->addMassMomentsDVSens(elemIndex, pt, X, dfdmoments, dvLen, dfdx);
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
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);
    TacsScalar density = moments[0];

    // Evaluate the self MOI
    TacsScalar dfdmoments[3];
    mat3x3SymmTransformSens(T, dfdq, dfdI0);
    dfdmoments[2] = scale * (dfdI0[0] + dfdI0[3]);
    dfdmoments[1] = -scale * 2.0 * moments[1] / density * (dfdI0[0] + dfdI0[3]);
    dfdmoments[0] = scale * moments[1] * moments[1] / density / density *
                    (dfdI0[0] + dfdI0[3]);

    TacsScalar dXcg[3];
    for (int i = 0; i < 3; i++) {
      dXcg[i] = X[i] + moments[1] / density * n0[i];
    }

    // Use parallel axis theorem to move MOI to origin
    dfdmoments[0] +=
        scale * dfdq[0] *
        (dXcg[1] * dXcg[1] + dXcg[2] * dXcg[2] -
         2.0 * moments[1] / density * (dXcg[1] * n0[1] + dXcg[2] * n0[2]));
    dfdmoments[0] -=
        scale * dfdq[1] *
        (dXcg[0] * dXcg[1] -
         moments[1] / density * (dXcg[0] * n0[1] + dXcg[1] * n0[0]));
    dfdmoments[0] -=
        scale * dfdq[2] *
        (dXcg[0] * dXcg[2] -
         moments[1] / density * (dXcg[0] * n0[2] + dXcg[2] * n0[0]));
    dfdmoments[0] +=
        scale * dfdq[3] *
        (dXcg[0] * dXcg[0] + dXcg[2] * dXcg[2] -
         2.0 * moments[1] / density * (dXcg[0] * n0[0] + dXcg[2] * n0[2]));
    dfdmoments[0] -=
        scale * dfdq[4] *
        (dXcg[2] * dXcg[1] -
         moments[1] / density * (dXcg[1] * n0[2] + dXcg[2] * n0[1]));
    dfdmoments[0] +=
        scale * dfdq[5] *
        (dXcg[0] * dXcg[0] + dXcg[1] * dXcg[1] -
         2.0 * moments[1] / density * (dXcg[0] * n0[0] + dXcg[1] * n0[1]));

    dfdmoments[1] +=
        scale * dfdq[0] * 2.0 * (dXcg[1] * n0[1] + dXcg[2] * n0[2]);
    dfdmoments[1] -= scale * dfdq[1] * (dXcg[0] * n0[1] + dXcg[1] * n0[0]);
    dfdmoments[1] -= scale * dfdq[2] * (dXcg[0] * n0[2] + dXcg[2] * n0[0]);
    dfdmoments[1] +=
        scale * dfdq[3] * 2.0 * (dXcg[0] * n0[0] + dXcg[2] * n0[2]);
    dfdmoments[1] -= scale * dfdq[4] * (dXcg[1] * n0[2] + dXcg[2] * n0[1]);
    dfdmoments[1] +=
        scale * dfdq[5] * 2.0 * (dXcg[0] * n0[0] + dXcg[1] * n0[1]);

    con->addMassMomentsDVSens(elemIndex, pt, X, dfdmoments, dvLen, dfdx);
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
  Bilinear contraction of the 6 Hessian blocks returned by
  model::evalStrainHessian against two directional strain-gradient tuples
  (psi_u0x,psi_u1x,psi_e0ty) and (phi_u0x,phi_u1x,phi_e0ty).

  Row/column conventions (confirmed by inspection of TacsShellAddDispGradHessian
  and TacsShellAddTyingDispCoupling, the routines that scatter these same
  blocks into addJacobian's mat[]):
    d2u0x, d2u1x: symmetric 9x9 self-blocks - contract once.
    d2u0xu1x: (row=u0x, col=u1x) cross-block - contract BOTH orderings,
      since the full assembled Hessian has this block AND its transpose
      (psi_u0x . d2u0xu1x . phi_u1x) + (phi_u0x . d2u0xu1x . psi_u1x).
    d2e0ty: symmetric 6x6 self-block - contract once.
    d2e0tyu0x, d2e0tyu1x: (row=e0ty(6), col=disp(9)) cross-blocks - contract
      both orderings, same reasoning as d2u0xu1x.
*/
static inline TacsScalar TacsShellContractStrainHessian(
    const TacsScalar d2u0x[], const TacsScalar d2u1x[],
    const TacsScalar d2u0xu1x[], const TacsScalar d2e0ty[],
    const TacsScalar d2e0tyu0x[], const TacsScalar d2e0tyu1x[],
    const TacsScalar psi_u0x[], const TacsScalar psi_u1x[],
    const TacsScalar psi_e0ty[], const TacsScalar phi_u0x[],
    const TacsScalar phi_u1x[], const TacsScalar phi_e0ty[]) {
  TacsScalar val = 0.0;
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      val += psi_u0x[i] * d2u0x[9 * i + j] * phi_u0x[j];
      val += psi_u1x[i] * d2u1x[9 * i + j] * phi_u1x[j];
      val += psi_u0x[i] * d2u0xu1x[9 * i + j] * phi_u1x[j];
      val += phi_u0x[i] * d2u0xu1x[9 * i + j] * psi_u1x[j];
    }
  }
  for (int k = 0; k < 6; k++) {
    for (int l = 0; l < 6; l++) {
      val += psi_e0ty[k] * d2e0ty[6 * k + l] * phi_e0ty[l];
    }
    for (int j = 0; j < 9; j++) {
      val += psi_e0ty[k] * d2e0tyu0x[9 * k + j] * phi_u0x[j];
      val += phi_e0ty[k] * d2e0tyu0x[9 * k + j] * psi_u0x[j];
      val += psi_e0ty[k] * d2e0tyu1x[9 * k + j] * phi_u1x[j];
      val += phi_e0ty[k] * d2e0tyu1x[9 * k + j] * psi_u1x[j];
    }
  }
  return val;
}

/*
  Compute E''(psi,phi): the strain's own second directional derivative along
  (psi,phi), i.e. the 9-component vector whose k-th entry is
  psiU^T * (d^2 e_k / dU^2) * phiU, where U = (u0x,u1x,e0ty).

  model::evalStrainHessian's Hessian blocks are additive in (Cs-quadratic
  term) + (s-linear term); calling it with Cs=0 kills the Cs-quadratic part,
  leaving exactly the s-linear part, which is itself linear in the passed-in
  "stress" argument. Looping over the 8 unit-stress directions therefore
  extracts E''(psi,phi) component-by-component without any model-class
  special-casing (linear strain models return exactly zero, since their
  evalStrainHessian has no s-dependence at all).

  This does NOT capture the tying strain's own vars-curvature (e0ty(vars) is
  quadratic in vars for the nonlinear model classes, independent of
  model::evalStrain's u0x/u1x nonlinearity) - see
  TacsShellAddTyingStrainCurvature below for that additive piece.
*/
template <class model>
void TacsShellAddStrainHessianBilinear(
    const TacsScalar u0x[], const TacsScalar u1x[], const TacsScalar e0ty[],
    const TacsScalar psi_u0x[], const TacsScalar psi_u1x[],
    const TacsScalar psi_e0ty[], const TacsScalar phi_u0x[],
    const TacsScalar phi_u1x[], const TacsScalar phi_e0ty[],
    TacsScalar Epp[]) {
  TacsScalar Cs_zero[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
  memset(Cs_zero, 0, sizeof(Cs_zero));

  for (int k = 0; k < 8; k++) {
    TacsScalar sK[9];
    memset(sK, 0, sizeof(sK));
    sK[k] = 1.0;

    TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
    TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
    model::evalStrainHessian(1.0, sK, Cs_zero, u0x, u1x, e0ty, d2u0x, d2u1x,
                             d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x);

    Epp[k] += TacsShellContractStrainHessian(
        d2u0x, d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x, psi_u0x, psi_u1x,
        psi_e0ty, phi_u0x, phi_u1x, phi_e0ty);
  }
}

/*
  One-sided Hessian-vector product: out += H*dir, where H is the same
  block-matrix Hessian TacsShellContractStrainHessian contracts bilinearly
  (psi^T*H*phi as a scalar) - here instead producing the vector H*dir (one
  index contracted against "dir", the other left free as the output "out",
  in (u0x,u1x,e0ty)-space). Row/column/transpose conventions match
  TacsShellContractStrainHessian exactly (verified by construction: summing
  psi[.] . out[.] over all three sub-spaces, for out=H*phi, reproduces
  TacsShellContractStrainHessian(psi,phi) term-for-term).

  ADDS into out_u0x/out_u1x/out_e0ty (does not zero them first), so that two
  calls (e.g. weighted by s_phi against psi's direction, and by s_psi against
  phi's direction) can be accumulated into one seed before scattering.
*/
static inline void TacsShellAddStrainHessianVector(
    const TacsScalar d2u0x[], const TacsScalar d2u1x[],
    const TacsScalar d2u0xu1x[], const TacsScalar d2e0ty[],
    const TacsScalar d2e0tyu0x[], const TacsScalar d2e0tyu1x[],
    const TacsScalar dir_u0x[], const TacsScalar dir_u1x[],
    const TacsScalar dir_e0ty[], TacsScalar out_u0x[], TacsScalar out_u1x[],
    TacsScalar out_e0ty[]) {
  for (int i = 0; i < 9; i++) {
    TacsScalar au0x = 0.0, au1x = 0.0;
    for (int j = 0; j < 9; j++) {
      au0x += d2u0x[9 * i + j] * dir_u0x[j] + d2u0xu1x[9 * i + j] * dir_u1x[j];
      au1x += d2u1x[9 * i + j] * dir_u1x[j] + d2u0xu1x[9 * j + i] * dir_u0x[j];
    }
    for (int k = 0; k < 6; k++) {
      au0x += d2e0tyu0x[9 * k + i] * dir_e0ty[k];
      au1x += d2e0tyu1x[9 * k + i] * dir_e0ty[k];
    }
    out_u0x[i] += au0x;
    out_u1x[i] += au1x;
  }
  for (int k = 0; k < 6; k++) {
    TacsScalar ae0ty = 0.0;
    for (int l = 0; l < 6; l++) {
      ae0ty += d2e0ty[6 * k + l] * dir_e0ty[l];
    }
    for (int j = 0; j < 9; j++) {
      ae0ty += d2e0tyu0x[9 * k + j] * dir_u0x[j] +
               d2e0tyu1x[9 * k + j] * dir_u1x[j];
    }
    out_e0ty[k] += ae0ty;
  }
}

/*
  Compute the tying strain's own vars-curvature: the per-tying-point
  quantity c[index] = psi^T * (d^2 ety[index] / dvars^2) * phi.

  model::computeTyingStrainDeriv(Xpts,fn,vars,d,varsd,dd,ety,etyd) is, for
  fixed (varsd,dd), an AFFINE function of (vars,d): a constant term from the
  varsd-only piece of each tying-strain field's formula, plus a term
  bilinear in (vars-role,varsd-role) (e.g. g11's 0.5*|Uxi(vars)|^2). The
  bilinear cross term is recovered exactly - via the polarization identity,
  not a finite difference, so there is no truncation error or step-size
  parameter - as:
     etyd(vars=phi, d=dd_phi; varsd=psi, dd=dd_psi)
   - etyd(vars=0,   d=0;      varsd=psi, dd=dd_psi)
  which cancels the constant (varsd-only) term exactly and leaves precisely
  the bilinear cross-contribution. This generalizes to any model class's
  tying-strain formula without special-casing: linear models' etyd has no
  bilinear part at all, so this difference is exactly zero for them.
*/
template <int vars_per_node, class basis, class model>
void TacsShellAddTyingStrainCurvature(const TacsScalar Xpts[],
                                      const TacsScalar fn[],
                                      const TacsScalar psi[],
                                      const TacsScalar dd_psi[],
                                      const TacsScalar phi[],
                                      const TacsScalar dd_phi[],
                                      TacsScalar c_ety[]) {
  const int dsize = 3 * basis::NUM_NODES;
  TacsScalar zeros_v[vars_per_node * basis::NUM_NODES];
  TacsScalar zeros_d[3 * basis::NUM_NODES];
  memset(zeros_v, 0, vars_per_node * basis::NUM_NODES * sizeof(TacsScalar));
  memset(zeros_d, 0, dsize * sizeof(TacsScalar));

  TacsScalar ety_fwd[basis::NUM_TYING_POINTS], etyd_fwd[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn, phi, dd_phi, psi, dd_psi, ety_fwd, etyd_fwd);

  TacsScalar ety_base[basis::NUM_TYING_POINTS],
      etyd_base[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
      Xpts, fn, zeros_v, zeros_d, psi, dd_psi, ety_base, etyd_base);

  for (int index = 0; index < basis::NUM_TYING_POINTS; index++) {
    c_ety[index] = etyd_fwd[index] - etyd_base[index];
  }
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    addMatDVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                             double time, TacsScalar scale,
                             const TacsScalar psi[], const TacsScalar phi[],
                             const TacsScalar Xpts[], const TacsScalar vars[],
                             int dvLen, TacsScalar dfdx[]) {
  // The geometric stiffness matrix requires a hand-unrolled third
  // derivative through getMatType's directional-derivative code path -
  // out of scope for this feature (documented scope cut, SPEC.md sec 4.2).
  // Forward to the base-class FD/CS implementation.
  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TACSElement::addMatDVSensInnerProduct(matType, elemIndex, time, scale,
                                          psi, phi, Xpts, vars, dvLen, dfdx);
    return;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    // psi^T*K*phi is the energy Hessian's bilinear form:
    //   psi^T*K*phi = e_psi^T*Cs*e_phi + s:E''(psi,phi)
    // where e_psi/e_phi are the strain's directional derivatives along
    // psi/phi (the same TacsShellComputeDispGradDeriv + evalStrainDeriv
    // pattern addAdjResProduct already uses for a single adjoint direction),
    // s = Cs*e is the primal stress, and E''(psi,phi) is the strain's own
    // second directional derivative (nonzero only where the strain is
    // nonlinear in the state). d/dx of each term is a generic
    // con->addStressDVSens call, exactly mirroring how addAdjResProduct
    // calls it for the single-direction case.
    //
    // E''(psi,phi) has two additive sources, both required and both handled
    // without model-class special-casing (HANDOFF-task-4.md's follow-up):
    //  (a) model::evalStrain's own nonlinearity in (u0x,u1x) - extracted via
    //      TacsShellAddStrainHessianBilinear (evalStrainHessian w/ Cs=0).
    //  (b) the tying strain's own nonlinearity in vars (e.g.
    //      TACSShellNonlinearModel::computeTyingStrain's 0.5*|Uxi|^2 term) -
    //      a genuinely separate source from (a), extracted via
    //      TacsShellAddTyingStrainCurvature's polarization identity and
    //      routed through model::evalStrain's (linear) e0ty pass-through.
    //
    // This is exact for TACSLinearizedRotation (director field linear in
    // vars, verified against getMatType to ~1e-15 relative). For
    // TACSQuadraticRotation/TACSQuaternionRotation, the director's own
    // curvature (D(vars,t) nonlinear in vars, HANDOFF-task-4.md's "second,
    // independent risk") is a genuinely separate, unimplemented additive
    // term - empirically confirmed (via TacsTestElementMatDVSens) to push
    // TACSQuadraticRotation's residual over this method's real-mode
    // rtol=1e-2 gate for some element/seed combinations (~1e-3 relative
    // error for others - inconsistent, not a fixed small ratio, i.e.
    // genuinely missing physics rather than roundoff). TACSQuaternionRotation
    // has the same structural gap (its director field is likewise nonlinear
    // in vars) but is not exercised by this repo's test suite at all, so it
    // cannot be empirically confirmed either way; both nonlinear-director
    // classes are forwarded to the base-class FD/CS implementation for
    // correctness. This is a narrower, new pattern in this file (an
    // analytic path retained for the one verified-exact director class,
    // fallback only for the others) - NOT the same shape as
    // getMatSVSensInnerProduct's TACS_MASS_MATRIX fix, which forwards to
    // base unconditionally for every director class (no typeid split at
    // all) rather than keeping any director on an analytic path.
    // TACSLinearizedRotation keeps the exact analytic path below.
    if (typeid(director) != typeid(TACSLinearizedRotation)) {
      TACSElement::addMatDVSensInnerProduct(matType, elemIndex, time, scale,
                                            psi, phi, Xpts, vars, dvLen, dfdx);
      return;
    }

    const int nquad = quadrature::getNumQuadraturePoints();

    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar zeros[vars_per_node * num_nodes];
    memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

    TacsScalar etn[num_nodes], etnd_psi[num_nodes], etnd_phi[num_nodes];
    TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
    TacsScalar u0xn[9 * num_nodes], Ctn[csize];
    TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                     model>(transform, Xdn, fn, vars, psi,
                                            XdinvTn, Tn, u0xn, Ctn, etn,
                                            etnd_psi);
    {
      TacsScalar etn_tmp[num_nodes];
      TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                       model>(transform, Xdn, fn, vars, phi,
                                              XdinvTn, Tn, u0xn, Ctn, etn_tmp,
                                              etnd_phi);
    }

    TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd_psi[dsize], dd_phi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, psi, fn, d, ddot, dddot, dd_psi);
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, phi, fn, d, ddot, dddot, dd_phi);

    TacsScalar ety[basis::NUM_TYING_POINTS], etyd_psi[basis::NUM_TYING_POINTS],
        etyd_phi[basis::NUM_TYING_POINTS];
    model::template computeTyingStrainDeriv<vars_per_node, basis>(
        Xpts, fn, vars, d, psi, dd_psi, ety, etyd_psi);
    {
      TacsScalar ety_tmp[basis::NUM_TYING_POINTS];
      model::template computeTyingStrainDeriv<vars_per_node, basis>(
          Xpts, fn, vars, d, phi, dd_phi, ety_tmp, etyd_phi);
    }

    TacsScalar c_ety[basis::NUM_TYING_POINTS];
    TacsShellAddTyingStrainCurvature<vars_per_node, basis, model>(
        Xpts, fn, psi, dd_psi, phi, dd_phi, c_ety);

    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar X[3], Xxi[6], n0[3], T[9], et, etd_psi, etd_phi;
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFields<1, 1>(pt, etn, &et);
      basis::template interpFields<1, 1>(pt, etnd_psi, &etd_psi);
      basis::template interpFields<1, 1>(pt, etnd_phi, &etd_phi);

      transform->computeTransform(Xxi, n0, T);

      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9], u0xd_psi[9], u1xd_psi[9];
      TacsScalar detXd = TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, psi, dd_psi, XdinvT, XdinvzT, u0x,
          u1x, u0xd_psi, u1xd_psi);
      detXd *= weight;

      TacsScalar u0x_tmp[9], u1x_tmp[9], u0xd_phi[9], u1xd_phi[9];
      TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, phi, dd_phi, XdinvT, XdinvzT,
          u0x_tmp, u1x_tmp, u0xd_phi, u1xd_phi);

      TacsScalar gty[6], gtyd_psi[6], gtyd_phi[6];
      basis::interpTyingStrain(pt, ety, gty);
      basis::interpTyingStrain(pt, etyd_psi, gtyd_psi);
      basis::interpTyingStrain(pt, etyd_phi, gtyd_phi);

      TacsScalar e0ty[6], e0tyd_psi[6], e0tyd_phi[6];
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_psi, e0tyd_psi);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_phi, e0tyd_phi);

      TacsScalar e[9], e_psi[9], e_phi[9];
      model::evalStrain(u0x, u1x, e0ty, e);
      e[8] = et;
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_psi, u1xd_psi, e0tyd_psi, e,
                             e_psi);
      e_psi[8] = etd_psi;
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_phi, u1xd_phi, e0tyd_phi, e,
                             e_phi);
      e_phi[8] = etd_phi;

      // Cs-quadratic term: d/dx[e_psi^T*Cs*e_phi]
      con->addStressDVSens(elemIndex, scale * detXd, pt, X, e_phi, e_psi,
                           dvLen, dfdx);

      // Kinematic term: d/dx[e^T*Cs*E''(psi,phi)]
      TacsScalar Epp[9];
      memset(Epp, 0, sizeof(Epp));
      TacsShellAddStrainHessianBilinear<model>(
          u0x, u1x, e0ty, u0xd_psi, u1xd_psi, e0tyd_psi, u0xd_phi, u1xd_phi,
          e0tyd_phi, Epp);

      TacsScalar gty_curv[6];
      basis::interpTyingStrain(pt, c_ety, gty_curv);
      TacsScalar e0ty_curv[6];
      mat3x3SymmTransformTranspose(XdinvT, gty_curv, e0ty_curv);
      TacsScalar zeros9[9];
      memset(zeros9, 0, sizeof(zeros9));
      TacsScalar Epp_ety[9];
      model::evalStrain(zeros9, zeros9, e0ty_curv, Epp_ety);
      for (int k = 0; k < 8; k++) {
        Epp[k] += Epp_ety[k];
      }
      Epp[8] = 0.0;

      con->addStressDVSens(elemIndex, scale * detXd, pt, X, e, Epp, dvLen,
                           dfdx);
    }
    return;
  } else if (matType == TACS_MASS_MATRIX) {
    // The mass matrix has no strain/stiffness path at all - its only DV
    // dependence is through con->evalMassMoments. psi^T*M*phi is a bilinear
    // form in the (translational, director-directional-derivative) fields
    // interpolated from psi and phi respectively, exactly mirroring the
    // dynamics term addAdjResProduct already accumulates from (dvars,ddvars)
    // and the psi-direction director derivative - here both "velocity" and
    // "adjoint direction" roles are played by psi and phi respectively.
    const int nquad = quadrature::getNumQuadraturePoints();

    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar zeros[vars_per_node * num_nodes];
    memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

    // Directional derivative of the director field along psi and along phi
    TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd_psi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, psi, fn, d, ddot, dddot, dd_psi);

    TacsScalar d2[dsize], ddot2[dsize], dddot2[dsize], dd_phi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, phi, fn, d2, ddot2, dddot2, dd_phi);

    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar X[3], Xxi[6], n0[3], T[9];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      transform->computeTransform(Xxi, n0, T);

      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9];
      TacsScalar detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);
      detXd *= weight;

      TacsScalar psi_u0[3], phi_u0[3], psi_d[3], phi_d[3];
      basis::template interpFields<vars_per_node, 3>(pt, psi, psi_u0);
      basis::template interpFields<vars_per_node, 3>(pt, phi, phi_u0);
      basis::template interpFields<3, 3>(pt, dd_psi, psi_d);
      basis::template interpFields<3, 3>(pt, dd_phi, phi_d);

      TacsScalar coef[3];
      coef[0] = scale * detXd * vec3Dot(psi_u0, phi_u0);
      coef[1] = scale * detXd *
                (vec3Dot(psi_u0, phi_d) + vec3Dot(psi_d, phi_u0));
      coef[2] = scale * detXd * vec3Dot(psi_d, phi_d);

      con->addMassMomentsDVSens(elemIndex, pt, X, coef, dvLen, dfdx);
    }
    return;
  } else {
    // Unsupported/unknown matType - forward to the base class rather than
    // leaving dfdx untouched.
    TACSElement::addMatDVSensInnerProduct(matType, elemIndex, time, scale,
                                          psi, phi, Xpts, vars, dvLen, dfdx);
    return;
  }
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    addMatXptSensInnerProduct(ElementMatrixType matType, int elemIndex,
                              double time, TacsScalar scale,
                              const TacsScalar psi[], const TacsScalar phi[],
                              const TacsScalar Xpts[], const TacsScalar vars[],
                              TacsScalar dfdXpts[]) {
  // The geometric stiffness matrix requires a hand-unrolled third
  // derivative through getMatType's directional-derivative code path -
  // out of scope for this feature (documented scope cut, SPEC.md sec 4.2).
  // Forward to the base-class FD/CS implementation.
  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale,
                                           psi, phi, Xpts, vars, dfdXpts);
    return;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    // dfdXpts = d/dXpts[psi^T*K*phi], the Xpts-adjoint of Task 4.1/4.3's
    //   psi^T*K*phi = e_psi^T*Cs*e_phi + s:E''(psi,phi)
    // decomposition (see HANDOFF-task-4.md's "Task 4.2" sections for the
    // full derivation). Neither e_psi nor e_phi plays the "primal e" role
    // addAdjResXptProduct's own bilinear form has (both are directional
    // derivatives, model::evalStrainDeriv-shaped) - Task 4.1's own e^T*Cs*Epp
    // kinematic term is what supplies a genuine "primal e" factor. Since
    // model::evalStrain's Hessian d^2e/dU^2 is a vars-independent CONSTANT
    // (Task 4.3's finding - every model's strain formula is at most
    // quadratic in U), Epp(psi,phi)'s own Xpts-dependence flows only through
    // psiU/phiU's geometric chains and the tying-curvature term's fn-
    // dependence - no third-derivative object exists anywhere here.
    //
    // Per quadrature point, f_q = scale*w*detXd*(e_psi.Cs.e_phi + s:Epp)
    // decomposes into seeds on three disp-grad-Xpts chain invocations
    // (psi-direction, phi-direction, primal) plus the standard closures:
    //  1. psi-direction chain: seed = evalStrainSens(s_phi=Cs*e_phi)
    //     [piece 1a] + H(s,Cs=0).phiU [piece 2b, nonlinear only].
    //  2. phi-direction chain: seed = evalStrainSens(s_psi=Cs*e_psi)
    //     [piece 1a] + H(s,Cs=0).psiU [piece 2b, nonlinear only].
    //  3. primal (vars,d) chain: seed = Hessian-coupling
    //     H(s_phi,Cs=0).psiU + H(s_psi,Cs=0).phiU [piece 1b, nonlinear only]
    //     + evalStrainSens(s_Epp=Cs*Epp_total) [piece 2a, nonlinear only].
    //  4. detXd term: ddetXd += scale*w*(e_psi.Cs.e_phi + s:Epp), folded into
    //     the chain calls' ddetXd slot exactly as addAdjResXptProduct does.
    //  5. Tying-curvature Xpts closure [nonlinear tying models only]: c_ety's
    //     own Xpts-dependence, via the polarization identity on the EXISTING
    //     addTyingStrainXptSens/XptSensDeriv pair.
    //  6. XdinvT-direction correction (e0ty = XdinvT^{T}*gty*XdinvT
    //     transform): the "chain-1-shaped sub-block" addAdjResXptProduct
    //     uses, applied TWICE - once for (gtyd_psi, de0tyd-from-s_phi), once
    //     for (gtyd_phi, de0tyd-from-s_psi) - no "chain 2" analog exists
    //     since neither e_psi nor e_phi is a primal-e-shaped factor.
    //
    // This is exact for TACSLinearizedRotation (mirrors Task 4.1/4.3's
    // guard); TACSQuadraticRotation/TACSQuaternionRotation have the same
    // unimplemented director-curvature gap flagged throughout this file's
    // other TACS_STIFFNESS_MATRIX branches, so they forward to base.
    if (typeid(director) != typeid(TACSLinearizedRotation)) {
      TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale,
                                             psi, phi, Xpts, vars, dfdXpts);
      return;
    }

    // STAGE-A-ONLY GUARD (HANDOFF-task-4.md staging): pieces 1b/2a/2b/5
    // above (all nonlinear-model-only) are not yet implemented - restrict
    // the analytic path to the two linear model classes for now, where
    // those pieces vanish identically (E''=0, tying strain linear in vars).
    // Nonlinear-model elements temporarily forward to base; this guard is
    // removed once pieces 1b/2a/2b/5 land.
    if (typeid(model) != typeid(TACSShellLinearModel) &&
        typeid(model) != typeid(TACSShellInplaneLinearModel)) {
      TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale,
                                             psi, phi, Xpts, vars, dfdXpts);
      return;
    }

    const int nquad = quadrature::getNumQuadraturePoints();

    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar dfn[3 * num_nodes];
    memset(dfn, 0, 3 * num_nodes * sizeof(TacsScalar));

    TacsScalar zeros[vars_per_node * num_nodes];
    memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

    TacsScalar etn[num_nodes], etnd_psi[num_nodes], etnd_phi[num_nodes];
    TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
    TacsScalar u0xn[9 * num_nodes], Ctn[csize];
    TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                     model>(transform, Xdn, fn, vars, psi,
                                            XdinvTn, Tn, u0xn, Ctn, etn,
                                            etnd_psi);
    {
      TacsScalar etn_tmp[num_nodes];
      TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                       model>(transform, Xdn, fn, vars, phi,
                                              XdinvTn, Tn, u0xn, Ctn, etn_tmp,
                                              etnd_phi);
    }

    TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd_psi[dsize],
        dd_phi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, psi, fn, d, ddot, dddot, dd_psi);
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, phi, fn, d, ddot, dddot, dd_phi);

    TacsScalar ety[basis::NUM_TYING_POINTS], etyd_psi[basis::NUM_TYING_POINTS],
        etyd_phi[basis::NUM_TYING_POINTS];
    model::template computeTyingStrainDeriv<vars_per_node, basis>(
        Xpts, fn, vars, d, psi, dd_psi, ety, etyd_psi);
    {
      TacsScalar ety_tmp[basis::NUM_TYING_POINTS];
      model::template computeTyingStrainDeriv<vars_per_node, basis>(
          Xpts, fn, vars, d, phi, dd_phi, ety_tmp, etyd_phi);
    }

    // Xpts-direction accumulators for the two direction chains (psi, phi) -
    // no primal-direction accumulator exists yet (Stage A has no primal
    // chain: pieces 1b/2a, the only sources that would populate one, are
    // nonlinear-model-only and guarded off above).
    TacsScalar dd_xpt_psi[dsize], dd_xpt_phi[dsize];
    memset(dd_xpt_psi, 0, dsize * sizeof(TacsScalar));
    memset(dd_xpt_phi, 0, dsize * sizeof(TacsScalar));

    TacsScalar dety_psi[basis::NUM_TYING_POINTS],
        dety_phi[basis::NUM_TYING_POINTS];
    memset(dety_psi, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(dety_phi, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    TacsScalar detn_c1_psi[num_nodes], detn_c1_phi[num_nodes];
    memset(detn_c1_psi, 0, num_nodes * sizeof(TacsScalar));
    memset(detn_c1_phi, 0, num_nodes * sizeof(TacsScalar));

    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar X[3], Xxi[6], n0[3], T[9], et, etd_psi, etd_phi;
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFields<1, 1>(pt, etn, &et);
      basis::template interpFields<1, 1>(pt, etnd_psi, &etd_psi);
      basis::template interpFields<1, 1>(pt, etnd_phi, &etd_phi);

      transform->computeTransform(Xxi, n0, T);

      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9], u0xd_psi[9], u1xd_psi[9];
      TacsScalar detXd = TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, psi, dd_psi, XdinvT, XdinvzT, u0x,
          u1x, u0xd_psi, u1xd_psi);
      detXd *= weight;

      TacsScalar u0x_tmp[9], u1x_tmp[9], u0xd_phi[9], u1xd_phi[9];
      TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, phi, dd_phi, XdinvT, XdinvzT,
          u0x_tmp, u1x_tmp, u0xd_phi, u1xd_phi);

      TacsScalar gty[6], gtyd_psi[6], gtyd_phi[6];
      basis::interpTyingStrain(pt, ety, gty);
      basis::interpTyingStrain(pt, etyd_psi, gtyd_psi);
      basis::interpTyingStrain(pt, etyd_phi, gtyd_phi);

      TacsScalar e0ty[6], e0tyd_psi[6], e0tyd_phi[6];
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_psi, e0tyd_psi);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_phi, e0tyd_phi);

      TacsScalar e_dummy[9], e_psi[9], e_phi[9];
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_psi, u1xd_psi, e0tyd_psi,
                             e_dummy, e_psi);
      e_psi[8] = etd_psi;
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_phi, u1xd_phi, e0tyd_phi,
                             e_dummy, e_phi);
      e_phi[8] = etd_phi;

      TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
      con->evalTangentStiffness(elemIndex, pt, X, Cs);
      TacsScalar drill;
      const TacsScalar *A, *B, *D, *As;
      TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As,
                                                      &drill);

      TacsScalar s_phi[9], s_psi[9];
      TACSShellConstitutive::computeStress(A, B, D, As, drill, e_phi, s_phi);
      TACSShellConstitutive::computeStress(A, B, D, As, drill, e_psi, s_psi);

      // Piece 4: the detXd-direction seed (e_psi.Cs.e_phi = dot(e_psi,s_phi))
      TacsScalar ddetXd_total =
          scale * (e_psi[0] * s_phi[0] + e_psi[1] * s_phi[1] +
                   e_psi[2] * s_phi[2] + e_psi[3] * s_phi[3] +
                   e_psi[4] * s_phi[4] + e_psi[5] * s_phi[5] +
                   e_psi[6] * s_phi[6] + e_psi[7] * s_phi[7] +
                   e_psi[8] * s_phi[8]);

      // Drill contribution: both e_psi[8]/e_phi[8] are directional (etd)
      // outputs (no primal et ever appears in this bilinear form, since
      // Epp[8] is always 0 - the drill strain has no Hessian coupling for
      // TACSLinearizedRotation), so only the two direction chains need a
      // drill seed - no primal-chain drill term.
      TacsScalar det_c1_psi = scale * detXd * s_phi[8];
      basis::template addInterpFieldsTranspose<1, 1>(pt, &det_c1_psi,
                                                     detn_c1_psi);
      TacsScalar det_c1_phi = scale * detXd * s_psi[8];
      basis::template addInterpFieldsTranspose<1, 1>(pt, &det_c1_phi,
                                                     detn_c1_phi);

      TacsScalar dXxi[6], dn0[3], dT[9];
      memset(dXxi, 0, sizeof(dXxi));
      memset(dn0, 0, sizeof(dn0));
      memset(dT, 0, sizeof(dT));

      // ---- psi-direction chain (piece 1a, seed weighted by s_phi) ----
      TacsScalar du0x_psi[9], du1x_psi[9], de0tyd_psi[6];
      model::evalStrainSens(scale * detXd, s_phi, u0x, u1x, du0x_psi,
                            du1x_psi, de0tyd_psi);

      TacsShellComputeDispGradXptSens<vars_per_node, basis>(
          pt, Xpts, psi, fn, dd_psi, Xxi, n0, T, XdinvT, XdinvzT, du0x_psi,
          du1x_psi, ddetXd_total * weight, dXxi, dn0, dfn, dT, dd_xpt_psi);

      TacsScalar dgtyd_psi[6];
      mat3x3SymmTransformTransSens(XdinvT, de0tyd_psi, dgtyd_psi);
      basis::addInterpTyingStrainTranspose(pt, dgtyd_psi, dety_psi);

      // ---- phi-direction chain (piece 1a, seed weighted by s_psi) ----
      TacsScalar du0x_phi[9], du1x_phi[9], de0tyd_phi[6];
      model::evalStrainSens(scale * detXd, s_psi, u0x, u1x, du0x_phi,
                            du1x_phi, de0tyd_phi);

      TacsShellComputeDispGradXptSens<vars_per_node, basis>(
          pt, Xpts, phi, fn, dd_phi, Xxi, n0, T, XdinvT, XdinvzT, du0x_phi,
          du1x_phi, 0.0, dXxi, dn0, dfn, dT, dd_xpt_phi);

      TacsScalar dgtyd_phi[6];
      mat3x3SymmTransformTransSens(XdinvT, de0tyd_phi, dgtyd_phi);
      basis::addInterpTyingStrainTranspose(pt, dgtyd_phi, dety_phi);

      // ---- Piece 6: XdinvT-direction correction, applied twice (no
      // "chain 2" sub-block - neither e_psi nor e_phi is primal-e-shaped) ----
      {
        TacsScalar gtyd_psi_full[9], gtyd_phi_full[9];
        TacsScalar de0tyd_psi_full[9], de0tyd_phi_full[9];

        gtyd_psi_full[0] = gtyd_psi[0];
        gtyd_psi_full[1] = gtyd_psi[1];
        gtyd_psi_full[2] = gtyd_psi[2];
        gtyd_psi_full[3] = gtyd_psi[1];
        gtyd_psi_full[4] = gtyd_psi[3];
        gtyd_psi_full[5] = gtyd_psi[4];
        gtyd_psi_full[6] = gtyd_psi[2];
        gtyd_psi_full[7] = gtyd_psi[4];
        gtyd_psi_full[8] = gtyd_psi[5];

        gtyd_phi_full[0] = gtyd_phi[0];
        gtyd_phi_full[1] = gtyd_phi[1];
        gtyd_phi_full[2] = gtyd_phi[2];
        gtyd_phi_full[3] = gtyd_phi[1];
        gtyd_phi_full[4] = gtyd_phi[3];
        gtyd_phi_full[5] = gtyd_phi[4];
        gtyd_phi_full[6] = gtyd_phi[2];
        gtyd_phi_full[7] = gtyd_phi[4];
        gtyd_phi_full[8] = gtyd_phi[5];

        de0tyd_psi_full[0] = de0tyd_psi[0];
        de0tyd_psi_full[1] = 0.5 * de0tyd_psi[1];
        de0tyd_psi_full[2] = 0.5 * de0tyd_psi[2];
        de0tyd_psi_full[3] = 0.5 * de0tyd_psi[1];
        de0tyd_psi_full[4] = de0tyd_psi[3];
        de0tyd_psi_full[5] = 0.5 * de0tyd_psi[4];
        de0tyd_psi_full[6] = 0.5 * de0tyd_psi[2];
        de0tyd_psi_full[7] = 0.5 * de0tyd_psi[4];
        de0tyd_psi_full[8] = de0tyd_psi[5];

        de0tyd_phi_full[0] = de0tyd_phi[0];
        de0tyd_phi_full[1] = 0.5 * de0tyd_phi[1];
        de0tyd_phi_full[2] = 0.5 * de0tyd_phi[2];
        de0tyd_phi_full[3] = 0.5 * de0tyd_phi[1];
        de0tyd_phi_full[4] = de0tyd_phi[3];
        de0tyd_phi_full[5] = 0.5 * de0tyd_phi[4];
        de0tyd_phi_full[6] = 0.5 * de0tyd_phi[2];
        de0tyd_phi_full[7] = 0.5 * de0tyd_phi[4];
        de0tyd_phi_full[8] = de0tyd_phi[5];

        TacsScalar dXdinvT[9];
        memset(dXdinvT, 0, sizeof(dXdinvT));
        TacsScalar W[9], dW[9], tmp[9];

        // sub-block 1: (gtyd_psi, de0tyd_psi)
        mat3x3MatMult(gtyd_psi_full, XdinvT, W);
        mat3x3MatMult(W, de0tyd_psi_full, tmp);
        for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];
        mat3x3MatMult(XdinvT, de0tyd_psi_full, dW);
        mat3x3MatMult(gtyd_psi_full, dW, tmp);
        for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

        // sub-block 2: (gtyd_phi, de0tyd_phi)
        mat3x3MatMult(gtyd_phi_full, XdinvT, W);
        mat3x3MatMult(W, de0tyd_phi_full, tmp);
        for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];
        mat3x3MatMult(XdinvT, de0tyd_phi_full, dW);
        mat3x3MatMult(gtyd_phi_full, dW, tmp);
        for (int i = 0; i < 9; i++) dXdinvT[i] += tmp[i];

        // Propagate dXdinvT back through XdinvT = Xdinv*T and
        // Xdinv = inv3x3(Xd), Xd = assembleFrame(Xxi, n0)
        TacsScalar Xd[9], Xdinv[9];
        TacsShellAssembleFrame(Xxi, n0, Xd);
        inv3x3(Xd, Xdinv);

        TacsScalar dXdinv[9];
        mat3x3MatTransMult(dXdinvT, T, dXdinv);
        mat3x3TransMatMult(Xdinv, dXdinvT, tmp);
        for (int i = 0; i < 9; i++) dT[i] += tmp[i];

        TacsScalar dXd[9];
        inv3x3Sens(Xdinv, dXdinv, dXd);

        dXxi[0] += dXd[0];
        dXxi[1] += dXd[1];
        dXxi[2] += dXd[3];
        dXxi[3] += dXd[4];
        dXxi[4] += dXd[6];
        dXxi[5] += dXd[7];

        dn0[0] += dXd[2];
        dn0[1] += dXd[5];
        dn0[2] += dXd[8];
      }

      // Fold the T-direction seed back onto Xxi/n0
      transform->addTransformSens(Xxi, n0, dT, dXxi, dn0);

      // Close the loop from Xxi/n0 onto the element-level Xpts output and
      // the nodal frame-normal field
      basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
      basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn);
    }

    // Close the drilling-strain contributions (both direction chains - no
    // primal-chain drill term, see the comment inside the loop above)
    TacsShellAddDrillStrainXptSensDeriv<vars_per_node, offset, basis, director,
                                       model>(transform, Xdn, fn, vars, psi,
                                              XdinvTn, Tn, u0xn, Ctn,
                                              detn_c1_psi, dfdXpts, dfn);
    TacsShellAddDrillStrainXptSensDeriv<vars_per_node, offset, basis, director,
                                       model>(transform, Xdn, fn, vars, phi,
                                              XdinvTn, Tn, u0xn, Ctn,
                                              detn_c1_phi, dfdXpts, dfn);

    // Close the tying-strain contributions (both direction chains)
    model::template addTyingStrainXptSens<vars_per_node, basis>(
        Xpts, fn, psi, dd_psi, dety_psi, dfdXpts, dfn, dd_xpt_psi);
    model::template addTyingStrainXptSens<vars_per_node, basis>(
        Xpts, fn, phi, dd_phi, dety_phi, dfdXpts, dfn, dd_xpt_phi);

    // Close the director-field contributions (both direction chains, no
    // primal-role seed - mirrors the TACS_MASS_MATRIX branch's own
    // zero-dd/two-call pattern)
    TacsScalar zero_dd[dsize];
    memset(zero_dd, 0, dsize * sizeof(TacsScalar));
    director::template addDirectorRefNormalSens<vars_per_node, offset,
                                                num_nodes>(vars, psi, fn,
                                                           zero_dd, dd_xpt_psi,
                                                           dfn);
    director::template addDirectorRefNormalSens<vars_per_node, offset,
                                                num_nodes>(vars, phi, fn,
                                                           zero_dd, dd_xpt_phi,
                                                           dfn);

    TacsShellAddNodeNormalsSens<basis>(Xpts, dfn, dfdXpts);
    return;
  } else if (matType == TACS_MASS_MATRIX) {
    // Mass has no strain path. psi^T*M*phi's Xpts-dependence comes from two
    // sources: (a) the shared detXd geometric factor (same detXd-only chain
    // used for the geometry-only addPointQuantityXptSens branches above),
    // and (b) the reference direction t=fn used inside the director
    // Jacobian - psi_d = D_i(vars, t)*psi and phi_d = D_i(vars, t)*phi are
    // each linear in t for fixed vars, so their t-sensitivity closes through
    // director::addDirectorRefNormalSens's psi-paired overload.
    const int nquad = quadrature::getNumQuadraturePoints();

    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar zeros[vars_per_node * num_nodes];
    memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

    TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd_psi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, psi, fn, d, ddot, dddot, dd_psi);

    TacsScalar d2[dsize], ddot2[dsize], dddot2[dsize], dd_phi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, phi, fn, d2, ddot2, dddot2, dd_phi);

    TacsScalar dfn[3 * num_nodes];
    memset(dfn, 0, 3 * num_nodes * sizeof(TacsScalar));

    // Nodal seeds accumulated across all quadrature points, paired with
    // psi's/phi's own directional-derivative closure respectively
    TacsScalar dd_seed_psi[dsize], dd_seed_phi[dsize];
    memset(dd_seed_psi, 0, dsize * sizeof(TacsScalar));
    memset(dd_seed_phi, 0, dsize * sizeof(TacsScalar));

    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar X[3], Xxi[6], n0[3], T[9];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      transform->computeTransform(Xxi, n0, T);

      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9];
      TacsScalar detXd = TacsShellComputeDispGrad<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, XdinvT, XdinvzT, u0x, u1x);
      detXd *= weight;

      TacsScalar psi_u0[3], phi_u0[3], psi_d[3], phi_d[3];
      basis::template interpFields<vars_per_node, 3>(pt, psi, psi_u0);
      basis::template interpFields<vars_per_node, 3>(pt, phi, phi_u0);
      basis::template interpFields<3, 3>(pt, dd_psi, psi_d);
      basis::template interpFields<3, 3>(pt, dd_phi, phi_d);

      TacsScalar moments[3];
      con->evalMassMoments(elemIndex, pt, X, moments);

      // (a) detXd-direction term
      TacsScalar coefTotal =
          scale * (moments[0] * vec3Dot(psi_u0, phi_u0) +
                   moments[1] * (vec3Dot(psi_u0, phi_d) +
                                vec3Dot(psi_d, phi_u0)) +
                   moments[2] * vec3Dot(psi_d, phi_d));

      TacsScalar Xd[9], Xdinv[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      TacsScalar detXd_raw = inv3x3(Xd, Xdinv);

      TacsScalar ddetXd_total = coefTotal * weight;
      TacsScalar dXd[9];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          dXd[3 * i + j] = ddetXd_total * detXd_raw * Xdinv[3 * j + i];
        }
      }

      TacsScalar dXxi[6], dn0[3];
      dXxi[0] = dXd[0];
      dXxi[1] = dXd[1];
      dXxi[2] = dXd[3];
      dXxi[3] = dXd[4];
      dXxi[4] = dXd[6];
      dXxi[5] = dXd[7];
      dn0[0] = dXd[2];
      dn0[1] = dXd[5];
      dn0[2] = dXd[8];

      TacsScalar dT[9];
      memset(dT, 0, sizeof(dT));
      transform->addTransformSens(Xxi, n0, dT, dXxi, dn0);

      basis::template addInterpFieldsGradTranspose<3, 3>(pt, dXxi, dfdXpts);
      basis::template addInterpFieldsTranspose<3, 3>(pt, dn0, dfn);

      // (b) t=fn-direction term through the director Jacobian
      TacsScalar seed_paired_with_psi[3], seed_paired_with_phi[3];
      for (int i = 0; i < 3; i++) {
        seed_paired_with_psi[i] =
            scale * detXd * (moments[1] * phi_u0[i] + moments[2] * phi_d[i]);
        seed_paired_with_phi[i] =
            scale * detXd * (moments[1] * psi_u0[i] + moments[2] * psi_d[i]);
      }
      basis::template addInterpFieldsTranspose<3, 3>(pt, seed_paired_with_psi,
                                                      dd_seed_psi);
      basis::template addInterpFieldsTranspose<3, 3>(pt, seed_paired_with_phi,
                                                      dd_seed_phi);
    }

    TacsScalar zero_dd[dsize];
    memset(zero_dd, 0, dsize * sizeof(TacsScalar));
    director::template addDirectorRefNormalSens<vars_per_node, offset,
                                                num_nodes>(
        vars, psi, fn, zero_dd, dd_seed_psi, dfn);
    director::template addDirectorRefNormalSens<vars_per_node, offset,
                                                num_nodes>(
        vars, phi, fn, zero_dd, dd_seed_phi, dfn);

    TacsShellAddNodeNormalsSens<basis>(Xpts, dfn, dfdXpts);
    return;
  } else {
    // Unsupported/unknown matType - forward to the base class rather than
    // leaving dfdXpts untouched.
    TACSElement::addMatXptSensInnerProduct(matType, elemIndex, time, scale,
                                           psi, phi, Xpts, vars, dfdXpts);
    return;
  }
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
    getMatSVSensInnerProduct(ElementMatrixType matType, int elemIndex,
                             double time, const TacsScalar psi[],
                             const TacsScalar phi[], const TacsScalar Xpts[],
                             const TacsScalar vars[], TacsScalar dfdu[]) {
  // getMatSVSensInnerProduct is an assignment (dfdu =), not an accumulation
  // (SPEC.md sec 3.5) - zero the output unconditionally as the first
  // statement, before the geometric-stiffness punt check.
  memset(dfdu, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

  // The geometric stiffness matrix requires a hand-unrolled third
  // derivative through getMatType's directional-derivative code path -
  // out of scope for this feature (documented scope cut, SPEC.md sec 4.2).
  // Forward to the base-class FD/CS implementation.
  if (matType == TACS_GEOMETRIC_STIFFNESS_MATRIX) {
    TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                          Xpts, vars, dfdu);
    return;
  } else if (matType == TACS_STIFFNESS_MATRIX) {
    // dfdu = d/dvars[psi^T*K*phi], differentiating Task 4.1's
    //   psi^T*K*phi = e_psi^T*Cs*e_phi + s:E''(psi,phi)
    // decomposition w.r.t. vars instead of design variables (see
    // HANDOFF-task-4.md's "Task 4.3" section for the full derivation this
    // mirrors). Unlike the Xpts-adjoint (addMatXptSensInnerProduct), no
    // geometric adjoint chain is needed: model::evalStrain's Hessian d^2e/dU^2
    // is a vars-independent constant for every model class here (every
    // model's strain formula is at most quadratic in U=(u0x,u1x,e0ty)), so
    // E''(psi,phi) itself has zero vars-dependence - only the state-dependent
    // stress-like weights, and psiU/phiU's own e0ty-component (nonlinear
    // tying-strain models only), carry vars-dependence. Three additive
    // pieces, applying the product rule to d/dvars[e_psi^T*Cs*e_phi] (pieces
    // 1 and 3) plus d/dvars[e^T*Cs*Epp] (piece 2):
    //  (1) "Cs-quadratic term's Hessian-coupling": e_psi = J(U)*psiU, and
    //      J(U)=dE/dU itself depends on the primal U(vars) for nonlinear
    //      models. d/dvars[J(U)*psiU] has a term (dJ/dU*dU/dvars)*psiU -
    //      a one-sided Hessian-vector product (TacsShellAddStrainHessianVector,
    //      new helper above) using the same Cs=0 evalStrainHessian extraction
    //      Task 4.1's TacsShellAddStrainHessianBilinear already validated,
    //      weighted by s_phi=Cs*e_phi (for psiU's own direction) and by
    //      s_psi=Cs*e_phi (for phiU's own direction, the symmetric term) -
    //      scattered via the ordinary vars-direction TacsShellAddDispGradSens/
    //      addComputeTyingStrainTranspose (no Xpts machinery needed, since
    //      T/XdinvT/XdinvzT are Xpts-only and unaffected by which "direction"
    //      the seed came from).
    //  (2) Kinematic term's weight: s_Epp=Cs*Epp_total (Epp_total from Task
    //      4.1, vars-independent as established above) routed through the
    //      ordinary single-direction model::evalStrainSens(detXd, s_Epp,
    //      u0x, u1x, ...) pattern addResidual/addAdjResProduct already use,
    //      scattered the same way.
    //  (3) psiU/phiU's own e0ty-component vars-curvature: e0tyd_psi (part of
    //      psiU) is affine, not constant, in vars for nonlinear tying-strain
    //      models - the same fact TacsShellAddTyingStrainCurvature (Task 4.1)
    //      already exploits, but there contracted against a second FIXED
    //      direction (phi); here it must be differentiated w.r.t. a generic
    //      vars direction and scattered. Closed via the same polarization
    //      identity a second time: accumulate a per-tying-point weight (the
    //      e0ty-component of J^T*s_phi / J^T*s_psi, i.e. exactly
    //      model::evalStrainSens's own de0ty output) across the quadrature
    //      loop, then close via model::addComputeTyingStrainTranspose(Xpts,
    //      fn, psi, dd_psi, weight, ...) MINUS addComputeTyingStrainTranspose
    //      (Xpts, fn, 0, 0, weight, ...) - the zero-baseline subtraction is
    //      required (unlike the Xpts-adjoint's analogous substitution trick)
    //      because addComputeTyingStrainTranspose's own gradient formula has
    //      a vars-independent "+Xxi"-type term (depends only on geometry,
    //      not on the substituted vars/d argument) that does not otherwise
    //      cancel; substituting vars=psi (or phi) into the linearization-
    //      point argument recovers exactly the bilinear cross-term's own
    //      vars-gradient once that Xxi-only piece is subtracted off.
    //
    // Exact for TACSLinearizedRotation (director field linear in vars, so
    // u0x/u1x/psiU/phiU's non-e0ty components are vars-independent and the
    // drill strain - also linear in vars for this director - contributes no
    // Hessian-coupling of its own); TACSQuadraticRotation/
    // TACSQuaternionRotation have the same unimplemented director-curvature
    // gap flagged throughout this file's other TACS_STIFFNESS_MATRIX
    // branches, so they forward to base, mirroring
    // addMatDVSensInnerProduct's guard exactly.
    if (typeid(director) != typeid(TACSLinearizedRotation)) {
      TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                            Xpts, vars, dfdu);
      return;
    }

    const int nquad = quadrature::getNumQuadraturePoints();

    TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
    TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

    TacsScalar zeros[vars_per_node * num_nodes];
    memset(zeros, 0, vars_per_node * num_nodes * sizeof(TacsScalar));

    TacsScalar etn[num_nodes], etnd_psi[num_nodes], etnd_phi[num_nodes];
    TacsScalar XdinvTn[9 * num_nodes], Tn[9 * num_nodes];
    TacsScalar u0xn[9 * num_nodes], Ctn[csize];
    TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                     model>(transform, Xdn, fn, vars, psi,
                                            XdinvTn, Tn, u0xn, Ctn, etn,
                                            etnd_psi);
    {
      TacsScalar etn_tmp[num_nodes];
      TacsShellComputeDrillStrainDeriv<vars_per_node, offset, basis, director,
                                       model>(transform, Xdn, fn, vars, phi,
                                              XdinvTn, Tn, u0xn, Ctn, etn_tmp,
                                              etnd_phi);
    }

    TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd_psi[dsize], dd_phi[dsize];
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, psi, fn, d, ddot, dddot, dd_psi);
    director::template computeDirectorRatesDeriv<vars_per_node, offset,
                                                 num_nodes>(
        vars, zeros, zeros, phi, fn, d, ddot, dddot, dd_phi);

    TacsScalar ety[basis::NUM_TYING_POINTS], etyd_psi[basis::NUM_TYING_POINTS],
        etyd_phi[basis::NUM_TYING_POINTS];
    model::template computeTyingStrainDeriv<vars_per_node, basis>(
        Xpts, fn, vars, d, psi, dd_psi, ety, etyd_psi);
    {
      TacsScalar ety_tmp[basis::NUM_TYING_POINTS];
      model::template computeTyingStrainDeriv<vars_per_node, basis>(
          Xpts, fn, vars, d, phi, dd_phi, ety_tmp, etyd_phi);
    }

    // Task 4.1's tying-strain-curvature term (vars-independent, reused as-is
    // for piece 2's Epp_total)
    TacsScalar c_ety[basis::NUM_TYING_POINTS];
    TacsShellAddTyingStrainCurvature<vars_per_node, basis, model>(
        Xpts, fn, psi, dd_psi, phi, dd_phi, c_ety);

    // Accumulators closed after the quadrature loop, exactly mirroring
    // addResidual's own dety/dd accumulate-then-close pattern.
    TacsScalar dd_total[dsize];
    memset(dd_total, 0, dsize * sizeof(TacsScalar));
    TacsScalar dety_total[basis::NUM_TYING_POINTS];
    memset(dety_total, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    // Piece 3's per-tying-point closure weights: weight_psi pairs with the
    // vars=psi substitution (sourced from s_phi's de0ty component),
    // weight_phi pairs with the vars=phi substitution (sourced from s_psi's).
    TacsScalar weight_psi[basis::NUM_TYING_POINTS],
        weight_phi[basis::NUM_TYING_POINTS];
    memset(weight_psi, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));
    memset(weight_phi, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

    TacsScalar Cs_zero[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
    memset(Cs_zero, 0, sizeof(Cs_zero));

    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar X[3], Xxi[6], n0[3], T[9], et, etd_psi, etd_phi;
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);
      basis::template interpFields<1, 1>(pt, etn, &et);
      basis::template interpFields<1, 1>(pt, etnd_psi, &etd_psi);
      basis::template interpFields<1, 1>(pt, etnd_phi, &etd_phi);

      transform->computeTransform(Xxi, n0, T);

      TacsScalar XdinvT[9], XdinvzT[9];
      TacsScalar u0x[9], u1x[9], u0xd_psi[9], u1xd_psi[9];
      TacsScalar detXd = TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, psi, dd_psi, XdinvT, XdinvzT, u0x,
          u1x, u0xd_psi, u1xd_psi);
      detXd *= weight;

      TacsScalar u0x_tmp[9], u1x_tmp[9], u0xd_phi[9], u1xd_phi[9];
      TacsShellComputeDispGradDeriv<vars_per_node, basis>(
          pt, Xpts, vars, fn, d, Xxi, n0, T, phi, dd_phi, XdinvT, XdinvzT,
          u0x_tmp, u1x_tmp, u0xd_phi, u1xd_phi);

      TacsScalar gty[6], gtyd_psi[6], gtyd_phi[6];
      basis::interpTyingStrain(pt, ety, gty);
      basis::interpTyingStrain(pt, etyd_psi, gtyd_psi);
      basis::interpTyingStrain(pt, etyd_phi, gtyd_phi);

      TacsScalar e0ty[6], e0tyd_psi[6], e0tyd_phi[6];
      mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_psi, e0tyd_psi);
      mat3x3SymmTransformTranspose(XdinvT, gtyd_phi, e0tyd_phi);

      TacsScalar e[9], e_psi[9], e_phi[9];
      model::evalStrain(u0x, u1x, e0ty, e);
      e[8] = et;
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_psi, u1xd_psi, e0tyd_psi, e,
                             e_psi);
      e_psi[8] = etd_psi;
      model::evalStrainDeriv(u0x, u1x, e0ty, u0xd_phi, u1xd_phi, e0tyd_phi, e,
                             e_phi);
      e_phi[8] = etd_phi;

      TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
      con->evalTangentStiffness(elemIndex, pt, X, Cs);
      TacsScalar drill;
      const TacsScalar *A, *B, *D, *As;
      TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As,
                                                      &drill);

      TacsScalar s_phi[9], s_psi[9];
      TACSShellConstitutive::computeStress(A, B, D, As, drill, e_phi, s_phi);
      TACSShellConstitutive::computeStress(A, B, D, As, drill, e_psi, s_psi);

      // Piece 1: one-sided Hessian-vector products, summed into one seed.
      TacsScalar out_u0x[9], out_u1x[9], out_e0ty[6];
      memset(out_u0x, 0, sizeof(out_u0x));
      memset(out_u1x, 0, sizeof(out_u1x));
      memset(out_e0ty, 0, sizeof(out_e0ty));

      TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
      TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
      model::evalStrainHessian(detXd, s_phi, Cs_zero, u0x, u1x, e0ty, d2u0x,
                               d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x);
      TacsShellAddStrainHessianVector(d2u0x, d2u1x, d2u0xu1x, d2e0ty,
                                      d2e0tyu0x, d2e0tyu1x, u0xd_psi,
                                      u1xd_psi, e0tyd_psi, out_u0x, out_u1x,
                                      out_e0ty);

      model::evalStrainHessian(detXd, s_psi, Cs_zero, u0x, u1x, e0ty, d2u0x,
                               d2u1x, d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x);
      TacsShellAddStrainHessianVector(d2u0x, d2u1x, d2u0xu1x, d2e0ty,
                                      d2e0tyu0x, d2e0tyu1x, u0xd_phi,
                                      u1xd_phi, e0tyd_phi, out_u0x, out_u1x,
                                      out_e0ty);

      // Piece 2: kinematic term's weight (Epp_total is vars-independent).
      TacsScalar Epp[9];
      memset(Epp, 0, sizeof(Epp));
      TacsShellAddStrainHessianBilinear<model>(
          u0x, u1x, e0ty, u0xd_psi, u1xd_psi, e0tyd_psi, u0xd_phi, u1xd_phi,
          e0tyd_phi, Epp);

      TacsScalar gty_curv[6];
      basis::interpTyingStrain(pt, c_ety, gty_curv);
      TacsScalar e0ty_curv[6];
      mat3x3SymmTransformTranspose(XdinvT, gty_curv, e0ty_curv);
      TacsScalar zeros9[9];
      memset(zeros9, 0, sizeof(zeros9));
      TacsScalar Epp_ety[9];
      model::evalStrain(zeros9, zeros9, e0ty_curv, Epp_ety);
      for (int k = 0; k < 8; k++) {
        Epp[k] += Epp_ety[k];
      }
      Epp[8] = 0.0;

      TacsScalar s_Epp[9];
      TACSShellConstitutive::computeStress(A, B, D, As, drill, Epp, s_Epp);

      TacsScalar du0x2[9], du1x2[9], de0ty2[6];
      model::evalStrainSens(detXd, s_Epp, u0x, u1x, du0x2, du1x2, de0ty2);
      for (int k = 0; k < 9; k++) {
        out_u0x[k] += du0x2[k];
        out_u1x[k] += du1x2[k];
      }
      for (int k = 0; k < 6; k++) {
        out_e0ty[k] += de0ty2[k];
      }

      // Scatter pieces 1+2's u0x/u1x part directly (ordinary vars-direction
      // dispgrad-sens, no Xpts machinery needed).
      TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                     out_u0x, out_u1x, dfdu,
                                                     dd_total);

      // ...and the e0ty part, accumulated into dety_total, closed once after
      // the loop with the REAL primal (vars,d) as the linearization point -
      // exactly like addResidual's own dety accumulator.
      TacsScalar dgty[6];
      mat3x3SymmTransformTransSens(XdinvT, out_e0ty, dgty);
      basis::addInterpTyingStrainTranspose(pt, dgty, dety_total);

      // Piece 3: psiU/phiU's own e0ty-component vars-curvature closure
      // weights, accumulated per tying point.
      TacsScalar junk_u0x[9], junk_u1x[9];
      TacsScalar weight_from_s_phi[6], weight_from_s_psi[6];
      model::evalStrainSens(detXd, s_phi, u0x, u1x, junk_u0x, junk_u1x,
                            weight_from_s_phi);
      model::evalStrainSens(detXd, s_psi, u0x, u1x, junk_u0x, junk_u1x,
                            weight_from_s_psi);

      TacsScalar dgty_psi[6], dgty_phi[6];
      mat3x3SymmTransformTransSens(XdinvT, weight_from_s_phi, dgty_psi);
      mat3x3SymmTransformTransSens(XdinvT, weight_from_s_psi, dgty_phi);
      basis::addInterpTyingStrainTranspose(pt, dgty_psi, weight_psi);
      basis::addInterpTyingStrainTranspose(pt, dgty_phi, weight_phi);
    }

    // Close pieces 1+2's e0ty accumulator at the real primal (vars,d) point.
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
        Xpts, fn, vars, d, dety_total, dfdu, dd_total);

    // Close piece 3: substituted-linearization-point calls, subtracting the
    // zero baseline to isolate the bilinear cross-term's own vars-gradient
    // (see the derivation note above the branch for why the subtraction is
    // required here, unlike the Xpts-adjoint's analogous substitution).
    {
      TacsScalar zeros_d[dsize];
      memset(zeros_d, 0, dsize * sizeof(TacsScalar));

      TacsScalar dfdu_tmp[vars_per_node * num_nodes], dd_tmp[dsize];

      memset(dfdu_tmp, 0, vars_per_node * num_nodes * sizeof(TacsScalar));
      memset(dd_tmp, 0, dsize * sizeof(TacsScalar));
      model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
          Xpts, fn, psi, dd_psi, weight_psi, dfdu_tmp, dd_tmp);
      for (int k = 0; k < vars_per_node * num_nodes; k++) {
        dfdu[k] += dfdu_tmp[k];
      }
      for (int k = 0; k < dsize; k++) {
        dd_total[k] += dd_tmp[k];
      }

      memset(dfdu_tmp, 0, vars_per_node * num_nodes * sizeof(TacsScalar));
      memset(dd_tmp, 0, dsize * sizeof(TacsScalar));
      model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
          Xpts, fn, zeros, zeros_d, weight_psi, dfdu_tmp, dd_tmp);
      for (int k = 0; k < vars_per_node * num_nodes; k++) {
        dfdu[k] -= dfdu_tmp[k];
      }
      for (int k = 0; k < dsize; k++) {
        dd_total[k] -= dd_tmp[k];
      }

      memset(dfdu_tmp, 0, vars_per_node * num_nodes * sizeof(TacsScalar));
      memset(dd_tmp, 0, dsize * sizeof(TacsScalar));
      model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
          Xpts, fn, phi, dd_phi, weight_phi, dfdu_tmp, dd_tmp);
      for (int k = 0; k < vars_per_node * num_nodes; k++) {
        dfdu[k] += dfdu_tmp[k];
      }
      for (int k = 0; k < dsize; k++) {
        dd_total[k] += dd_tmp[k];
      }

      memset(dfdu_tmp, 0, vars_per_node * num_nodes * sizeof(TacsScalar));
      memset(dd_tmp, 0, dsize * sizeof(TacsScalar));
      model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
          Xpts, fn, zeros, zeros_d, weight_phi, dfdu_tmp, dd_tmp);
      for (int k = 0; k < vars_per_node * num_nodes; k++) {
        dfdu[k] -= dfdu_tmp[k];
      }
      for (int k = 0; k < dsize; k++) {
        dd_total[k] -= dd_tmp[k];
      }
    }

    // Close the director-space accumulator into dfdu. TACSLinearizedRotation's
    // director Jacobian is vars-independent, so this is exact regardless of
    // what vars/dvars/ddvars are passed - mirrors the MASS_MATRIX branches'
    // own zeros-for-dvars/ddvars convention (this method has no dvars/ddvars
    // parameters of its own).
    director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
        vars, zeros, zeros, fn, dd_total, dfdu);
    return;
  } else if (matType == TACS_MASS_MATRIX) {
    // NOTE: for TACSLinearizedRotation the mass matrix is genuinely
    // state-independent (its rotational-DOF Jacobian is a fixed function of
    // the reference direction t only) and psi^T*M*phi's vars-derivative is
    // exactly zero, matching MITCShell's comment cited in VALIDATION.md
    // Claim 2. However, this was verified NOT to generalize:
    // TACSQuadraticRotation and TACSQuaternionRotation build the
    // rotational-rotational mass block from a director Jacobian D_i(vars,t)
    // that is itself a nonlinear function of vars, so mat's entries (and
    // hence psi^T*mat*phi) genuinely depend on vars for those two director
    // classes even though the mass moments/detXd do not (empirically
    // confirmed against the FD/CS harness: Quad4ShellModRot/Quad4Quaternion
    // both fail an exact-zero implementation while Quad4Shell - the
    // TACSLinearizedRotation default - passes). Deriving the analytic
    // third-derivative-like correction for the nonlinear-rotation directors
    // is out of scope for this pass (see HANDOFF-task-4.md); forward to the
    // base-class FD/CS implementation for correctness rather than assert an
    // incorrect zero.
    TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                          Xpts, vars, dfdu);
    return;
  } else {
    // Unsupported/unknown matType - forward to the base class rather than
    // leaving dfdu untouched.
    TACSElement::getMatSVSensInnerProduct(matType, elemIndex, time, psi, phi,
                                          Xpts, vars, dfdu);
    return;
  }
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::getAverageStresses(
    int elemIndex, ElementType etype, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], TacsScalar *avgStresses) {
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

    TacsScalar loc_avgStresses[9];
    // memset(loc_avgStresses,0,9);
    for (int i = 0; i < 9; i++) {
      loc_avgStresses[i] = 0.0;
    }

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

      for (int i = 0; i < 9; i++) {
        loc_avgStresses[i] += s[i];
      }
    }

    // average the average stresses among the quadrature points
    for (int i = 0; i < 9; i++) {
      loc_avgStresses[i] /= num_vis_nodes;
      avgStresses[i] += loc_avgStresses[i];
    }
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
        for (int failInd = 0; failInd < 7; failInd++) {
          data[failInd] =
              con->evalFailureFieldValue(elemIndex, pt, X, e, failInd);
        }
        for (int dvInd = 0; dvInd < 7; dvInd++) {
          data[dvInd + 7] = con->evalDesignFieldValue(elemIndex, pt, X, dvInd);
        }
        data += 14;
      }
      if (write_flag & TACS_OUTPUT_COORDINATE_FRAME) {
        data[0] = T[0];
        data[1] = T[3];
        data[2] = T[6];

        data[3] = T[1];
        data[4] = T[4];
        data[5] = T[7];

        data[6] = T[2];
        data[7] = T[5];
        data[8] = T[8];
        data += 9;
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