#ifndef TACS_THERMAL_SHELL_ELEMENT_H
#define TACS_THERMAL_SHELL_ELEMENT_H

#include "TACSDirector.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSElementVerification.h"
#include "TACSShellConstitutive.h"
#include "TACSShellElementModel.h"
#include "TACSShellElementTransform.h"
#include "TACSShellUtilities.h"

template <class quadrature, class basis, class director, class model>
class TACSThermalShellElement : public TACSElement {
 public:
  // Offset within the solution vector to the roational
  // parametrization defined via the director class. Here the offset
  // is 3 corresponding to the (u, v, w) displacements of the
  // mid-surface of the shell.
  static const int offset = 3;

  // The number of variables defined at each node of the shell
  // element.  There are 3 mid-surface displacements, the temperature
  // plus however many parameters are defined by the director class
  // for the parametrization.
  static const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // The number of nodes for this element. This is derived from the
  // basis function class. This is just a handy re-definition since
  // this constant is used in many locations within the element.
  static const int num_nodes = basis::NUM_NODES;

  // The thermal degree of freedom. This comes last
  static const int thermal_dof = vars_per_node - 1;

  TACSThermalShellElement(TACSShellTransform *_transform,
                          TACSShellConstitutive *_con) {
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  ~TACSThermalShellElement() {
    if (transform) {
      transform->decref();
    }

    if (con) {
      con->decref();
    }
  }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return num_nodes; }

  const char *getObjectName() { return "TACSThermalShellElement"; }

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
  static const int size = vars_per_node * num_nodes;
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
void TACSThermalShellElement<quadrature, basis, director,
                             model>::computeEnergies(int elemIndex, double time,
                                                     const TacsScalar *Xpts,
                                                     const TacsScalar *vars,
                                                     const TacsScalar *dvars,
                                                     TacsScalar *_Te,
                                                     TacsScalar *_Ue) {
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes];

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
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

    Ue += 0.5 * detXd *
          (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] + s[4] * e[4] +
           s[5] * e[5] + s[6] * e[6] + s[7] * e[7] + s[8] * e[8]);

    // Evaluate the mass moments
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    // Compute the velocities and the director velocities
    TacsScalar u0dot[3], d0dot[3];
    basis::template interpFields<vars_per_node, 3>(pt, dvars, u0dot);
    basis::template interpFields<3, 3>(pt, ddot, d0dot);

    Te += 0.5 * detXd *
          (moments[0] * vec3Dot(u0dot, u0dot) +
           2.0 * moments[1] * vec3Dot(u0dot, d0dot) +
           moments[2] * vec3Dot(d0dot, d0dot));
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSThermalShellElement<quadrature, basis, director, model>::addResidual(
    int elemIndex, double time, const TacsScalar *Xpts, const TacsScalar *vars,
    const TacsScalar *dvars, const TacsScalar *ddvars, TacsScalar *res) {
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Derivative of the director field and matrix at each point
  TacsScalar dd[dsize];
  memset(dd, 0, 3 * num_nodes * sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

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

    // Evaluate the temperature and temperature gradient
    TacsScalar t, txi[2];
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);
    basis::template interpFieldsGrad<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                       txi);

    // Transform to the local component of the heat flux
    TacsScalar tx[2];  // tx = txi*Xdinv*T
    tx[0] = XdinvT[0] * txi[0] + XdinvT[1] * txi[1];
    tx[1] = XdinvT[3] * txi[0] + XdinvT[4] * txi[1];

    // Compute the heat flux
    TacsScalar q[2];
    con->evalHeatFlux(elemIndex, pt, X, tx, q);

    TacsScalar qxi[2];
    qxi[0] = detXd * (XdinvT[0] * q[0] + XdinvT[3] * q[1]);
    qxi[1] = detXd * (XdinvT[1] * q[0] + XdinvT[4] * q[1]);
    basis::template addInterpFieldsGradTranspose<vars_per_node, 1>(
        pt, qxi, &res[thermal_dof]);

    // Compute the thermal strain
    TacsScalar eth[9];
    con->evalThermalStrain(elemIndex, pt, X, t, eth);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for (int i = 0; i < 9; i++) {
      em[i] = e[i] - eth[i];
    }

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, em, s);

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

  // Set the total number of tying points needed for this element
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
void TACSThermalShellElement<quadrature, basis, director, model>::addJacobian(
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

    // Evaluate the temperature and temperature gradient
    TacsScalar t, txi[2];
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);
    basis::template interpFieldsGrad<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                       txi);

    // Transform to the local component of the heat flux
    TacsScalar tx[2];  // tx = txi*Xdinv*T
    tx[0] = XdinvT[0] * txi[0] + XdinvT[1] * txi[1];
    tx[1] = XdinvT[3] * txi[0] + XdinvT[4] * txi[1];

    // Compute the heat flux
    TacsScalar q[2];
    con->evalHeatFlux(elemIndex, pt, X, tx, q);

    TacsScalar qxi[2];
    qxi[0] = detXd * (XdinvT[0] * q[0] + XdinvT[3] * q[1]);
    qxi[1] = detXd * (XdinvT[1] * q[0] + XdinvT[4] * q[1]);
    basis::template addInterpFieldsGradTranspose<vars_per_node, 1>(
        pt, qxi, &res[thermal_dof]);

    // Set the terms in the Jacobian matrix
    TacsScalar Kt[3];
    con->evalTangentHeatFlux(elemIndex, pt, X, Kt);

    // Compute the terms for the thermal stiffness matrix
    // [ XdinvT[0], XdinvT[3] ][ Kt[0], Kt[1] ][ XdinvT[0], XdinvT[1] ]
    // [ XdinvT[1], XdinvT[4] ][ Kt[1], Kt[2] ][ XdinvT[3], XdinvT[4] ]
    TacsScalar Ktmp[4];
    Ktmp[0] = Kt[0] * XdinvT[0] + Kt[1] * XdinvT[3];
    Ktmp[1] = Kt[0] * XdinvT[1] + Kt[1] * XdinvT[4];
    Ktmp[2] = Kt[1] * XdinvT[0] + Kt[2] * XdinvT[3];
    Ktmp[3] = Kt[1] * XdinvT[1] + Kt[2] * XdinvT[4];

    TacsScalar q2xi[4];
    q2xi[0] = alpha * detXd * (XdinvT[0] * Ktmp[0] + XdinvT[3] * Ktmp[2]);
    q2xi[1] = alpha * detXd * (XdinvT[0] * Ktmp[1] + XdinvT[3] * Ktmp[3]);
    q2xi[2] = alpha * detXd * (XdinvT[1] * Ktmp[0] + XdinvT[4] * Ktmp[2]);
    q2xi[3] = alpha * detXd * (XdinvT[1] * Ktmp[1] + XdinvT[4] * Ktmp[3]);

    basis::template addInterpGradOuterProduct<vars_per_node, vars_per_node, 1,
                                              1>(
        pt, q2xi, &mat[thermal_dof * (size + 1)]);

    // Compute the thermal strain
    TacsScalar eth[9];
    con->evalThermalStrain(elemIndex, pt, X, 1.0, eth);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for (int i = 0; i < 9; i++) {
      em[i] = e[i] - t * eth[i];
    }

    // Compute the tangent stiffness matrix
    TacsScalar Cs[TACSShellConstitutive::NUM_TANGENT_STIFFNESS_ENTRIES];
    con->evalTangentStiffness(elemIndex, pt, X, Cs);

    TacsScalar drill;
    const TacsScalar *A, *B, *D, *As;
    TACSShellConstitutive::extractTangentStiffness(Cs, &A, &B, &D, &As, &drill);

    // Compute the stress based on the tangent stiffness
    TacsScalar s[9];
    TACSShellConstitutive::computeStress(A, B, D, As, drill, em, s);

    {
      TacsScalar tmp[size], ddth[dsize];
      memset(tmp, 0, size * sizeof(TacsScalar));
      memset(ddth, 0, dsize * sizeof(TacsScalar));

      TacsScalar dthety[basis::NUM_TYING_POINTS];
      memset(dthety, 0, basis::NUM_TYING_POINTS * sizeof(TacsScalar));

      // Compute the thermal strain
      TacsScalar sth[9];
      TACSShellConstitutive::computeStress(A, B, D, As, drill, eth, sth);

      // Compute the derivative of the product of the stress and strain
      // with respect to u0x, u1x and e0ty
      TacsScalar dthu0x[9], dthu1x[9], dthe0ty[6];
      model::evalStrainSens(-alpha * detXd, sth, u0x, u1x, dthu0x, dthu1x,
                            dthe0ty);

      // Add the contributions to the residual from du0x, du1x and dCt
      TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                     dthu0x, dthu1x, tmp, ddth);

      // Compute the of the tying strain w.r.t. derivative w.r.t. the
      // coefficients
      TacsScalar dthgty[6];
      mat3x3SymmTransformTransSens(XdinvT, dthe0ty, dthgty);

      // Evaluate the tying strain
      basis::addInterpTyingStrainTranspose(pt, dthgty, dthety);

      // Add the contributions from the tying strain
      model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
          Xpts, fn, vars, d, dthety, tmp, ddth);

      // Add the contributions to the director field
      director::template addDirectorResidual<vars_per_node, offset, num_nodes>(
          vars, dvars, ddvars, fn, ddth, tmp);

      TacsScalar N[num_nodes], one = 1.0;
      memset(N, 0, num_nodes * sizeof(TacsScalar));
      basis::template addInterpFieldsTranspose<1, 1>(pt, &one, N);

      for (int j = 0; j < num_nodes; j++) {
        for (int i = 0; i < size; i++) {
          mat[i * size + (vars_per_node * j + thermal_dof)] += tmp[i] * N[j];
        }
      }
    }

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
    TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
    model::evalStrainHessian(alpha * detXd, s, Cs, u0x, u1x, e0ty, d2u0x, d2u1x,
                             d2u0xu1x, d2e0ty, d2e0tyu0x, d2e0tyu1x);

    // Add the components from the drilling penalty
    TacsScalar det = detXd * s[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det, detn);

    // Add the contributions to the residual from du0x and du1x
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

    // Add the coupling between the displacement and tying strain
    TacsShellAddTyingDispCoupling<basis>(pt, T, XdinvT, XdinvzT, d2e0tyu0x,
                                         d2e0tyu1x, d2etyu, d2etyd);

    // Evaluate the tying strain
    basis::addInterpTyingStrainTranspose(pt, dgty, dety);
    basis::addInterpTyingStrainHessian(pt, d2gty, d2ety);

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
void TACSThermalShellElement<quadrature, basis, director, model>::getMatType(
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
void TACSThermalShellElement<quadrature, basis, director, model>::
    addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                     const TacsScalar psi[], const TacsScalar Xpts[],
                     const TacsScalar vars[], const TacsScalar dvars[],
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

    // Evaluate the temperature and temperature gradient
    TacsScalar t, txi[2], txid[2];
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);
    basis::template interpFieldsGrad<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                       txi);
    basis::template interpFieldsGrad<vars_per_node, 1>(pt, &psi[thermal_dof],
                                                       txid);

    // Transform to the local component of the heat flux
    TacsScalar tx[2];  // tx = txi*Xdinv*T
    tx[0] = XdinvT[0] * txi[0] + XdinvT[1] * txi[1];
    tx[1] = XdinvT[3] * txi[0] + XdinvT[4] * txi[1];

    TacsScalar txd[2];
    txd[0] = XdinvT[0] * txid[0] + XdinvT[1] * txid[1];
    txd[1] = XdinvT[3] * txid[0] + XdinvT[4] * txid[1];

    // Add the contributions from the heat flux
    con->addHeatFluxDVSens(elemIndex, scale * detXd, pt, X, tx, txd, dvLen,
                           dfdx);

    // Add the contribution from the thermal strain
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, ed, s);
    con->addThermalStrainDVSens(elemIndex, pt, X, -t * scale * detXd, s, dvLen,
                                dfdx);

    // Compute the thermal strain
    TacsScalar eth[9];
    con->evalThermalStrain(elemIndex, pt, X, t, eth);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for (int i = 0; i < 9; i++) {
      em[i] = e[i] - eth[i];
    }

    // The directional derivative of the strain along the adjoint direction
    con->addStressDVSens(elemIndex, scale * detXd, pt, X, em, ed, dvLen, dfdx);

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
int TACSThermalShellElement<quadrature, basis, director,
                            model>::evalPointQuantity(int elemIndex,
                                                      int quantityType,
                                                      double time, int n,
                                                      double pt[],
                                                      const TacsScalar Xpts[],
                                                      const TacsScalar vars[],
                                                      const TacsScalar dvars[],
                                                      const TacsScalar ddvars[],
                                                      TacsScalar *detXd,
                                                      TacsScalar *quantity) {
  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn);

  if (quantityType == TACS_ELEMENT_DENSITY) {
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
  } else if (quantityType == TACS_HEAT_FLUX) {
    if (quantity) {
      // Compute X, X,xi
      TacsScalar X[3], Xxi[6], n0[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      // Compute the transformation at the quadrature point
      TacsScalar T[9];
      transform->computeTransform(Xxi, n0, T);

      // Evaluate the displacement gradient at the point
      TacsScalar Xd[9], Xdinv[9], XdinvT[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = inv3x3(Xd, Xdinv);
      mat3x3MatMult(Xdinv, T, XdinvT);

      // Evaluate the temperature and temperature gradient
      TacsScalar txi[2];
      basis::template interpFieldsGrad<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                         txi);

      // Transform to the local component of the heat flux
      TacsScalar tx[2];  // tx = txi*Xdinv*T
      tx[0] = XdinvT[0] * txi[0] + XdinvT[1] * txi[1];
      tx[1] = XdinvT[3] * txi[0] + XdinvT[4] * txi[1];

      con->evalHeatFlux(elemIndex, pt, X, tx, quantity);
    }

    return 2;
  } else if (quantityType == TACS_TEMPERATURE) {
    if (quantity) {
      TacsScalar Xxi[6], n0[3];
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
      basis::template interpFields<3, 3>(pt, fn, n0);

      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n0, Xd);
      *detXd = det3x3(Xd);

      TacsScalar t;
      basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                     &t);
      *quantity = t;
    }

    return 1;
  } else if (quantityType == TACS_FAILURE_INDEX ||
             quantityType == TACS_STRAIN_ENERGY_DENSITY ||
             quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
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

      // Evaluate the temperature
      TacsScalar t;
      basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof],
                                                     &t);

      // Compute the thermal strain
      TacsScalar eth[9];
      con->evalThermalStrain(elemIndex, pt, X, t, eth);

      // Compute the mechanical strain (and stress)
      TacsScalar em[9];
      for (int i = 0; i < 9; i++) {
        em[i] = e[i] - eth[i];
      }

      if (quantityType == TACS_FAILURE_INDEX) {
        *quantity = con->evalFailure(elemIndex, pt, X, em);
      } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
        TacsScalar s[9];
        con->evalStress(elemIndex, pt, X, em, s);

        *quantity = (s[0] * em[0] + s[1] * em[1] + s[2] * em[2] + s[3] * em[3] +
                     s[4] * em[4] + s[5] * em[5] + s[6] * em[6] + s[7] * em[7] +
                     s[8] * em[8]);
      } else if (quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
        TacsScalar s[9];
        con->evalStress(elemIndex, pt, X, e, s);

        *quantity = (s[0] * e[0] + s[1] * e[1] + s[2] * e[2] + s[3] * e[3] +
                     s[4] * e[4] + s[5] * e[5] + s[6] * e[6] + s[7] * e[7] +
                     s[8] * e[8]);
      }
    }

    return 1;
  }

  return 0;
}

template <class quadrature, class basis, class director, class model>
void TACSThermalShellElement<quadrature, basis, director, model>::
    addPointQuantityDVSens(int elemIndex, int quantityType, double time,
                           TacsScalar scale, int n, double pt[],
                           const TacsScalar Xpts[], const TacsScalar vars[],
                           const TacsScalar dvars[], const TacsScalar ddvars[],
                           const TacsScalar dfdq[], int dvLen,
                           TacsScalar dfdx[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY ||
      quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
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

    // Evaluate the temperature
    TacsScalar t;
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);

    // Compute the thermal strain
    TacsScalar eth[9];
    con->evalThermalStrain(elemIndex, pt, X, t, eth);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for (int i = 0; i < 9; i++) {
      em[i] = e[i] - eth[i];
    }

    if (quantityType == TACS_FAILURE_INDEX) {
      con->addFailureDVSens(elemIndex, scale * dfdq[0], pt, X, em, dvLen, dfdx);
    } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
      con->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, em, em, dvLen,
                           dfdx);
    } else {  // quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY
      con->addStressDVSens(elemIndex, scale * dfdq[0], pt, X, e, e, dvLen,
                           dfdx);
    }
  } else if (quantityType == TACS_ELEMENT_DENSITY) {
    TacsScalar X[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);

    con->addDensityDVSens(elemIndex, scale * dfdq[0], pt, X, dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSThermalShellElement<quadrature, basis, director, model>::
    addPointQuantitySVSens(int elemIndex, int quantityType, double time,
                           TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                           int n, double pt[], const TacsScalar Xpts[],
                           const TacsScalar vars[], const TacsScalar dvars[],
                           const TacsScalar ddvars[], const TacsScalar dfdq[],
                           TacsScalar dfdu[]) {
  if (quantityType == TACS_FAILURE_INDEX ||
      quantityType == TACS_STRAIN_ENERGY_DENSITY ||
      quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY) {
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

    // Evaluate the temperature
    TacsScalar t;
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);

    // Compute the derivative of the output with respect to the strain
    // components in the local frame
    TacsScalar sens[9];
    TacsScalar escale = 1.0;  // Scale factor for the sensitivity

    if (quantityType == TACS_FAILURE_INDEX) {
      // Compute the thermal strain
      TacsScalar eth[9];
      con->evalThermalStrain(elemIndex, pt, X, 1.0, eth);

      // Compute the mechanical strain (and stress)
      TacsScalar em[9];
      for (int i = 0; i < 9; i++) {
        em[i] = e[i] - t * eth[i];
      }

      // Compute the sensitivity of the failure index w.r.t. the strain
      con->evalFailureStrainSens(elemIndex, pt, X, em, sens);
      escale = alpha * dfdq[0];

      // Add contribution from the thermal part
      TacsScalar scale =
          -alpha * dfdq[0] *
          (sens[0] * eth[0] + sens[1] * eth[1] + sens[2] * eth[2] +
           sens[3] * eth[3] + sens[4] * eth[4] + sens[5] * eth[5] +
           sens[6] * eth[6] + sens[7] * eth[7] + sens[8] * eth[8]);
      basis::template addInterpFieldsTranspose<vars_per_node, 1>(
          pt, &scale, &dfdu[thermal_dof]);
    } else if (quantityType == TACS_STRAIN_ENERGY_DENSITY) {
      // Compute the thermal strain
      TacsScalar eth[9];
      con->evalThermalStrain(elemIndex, pt, X, 1.0, eth);

      // Compute the mechanical strain (and stress)
      TacsScalar em[9];
      for (int i = 0; i < 9; i++) {
        em[i] = e[i] - t * eth[i];
      }

      // Compute the sensitivity
      con->evalStress(elemIndex, pt, X, em, sens);
      escale = 2.0 * alpha * dfdq[0];

      // Add contribution from the thermal part
      TacsScalar scale =
          -2.0 * alpha * dfdq[0] *
          (sens[0] * eth[0] + sens[1] * eth[1] + sens[2] * eth[2] +
           sens[3] * eth[3] + sens[4] * eth[4] + sens[5] * eth[5] +
           sens[6] * eth[6] + sens[7] * eth[7] + sens[8] * eth[8]);
      basis::template addInterpFieldsTranspose<vars_per_node, 1>(
          pt, &scale, &dfdu[thermal_dof]);
    } else {  // quantityType == TACS_TOTAL_STRAIN_ENERGY_DENSITY
      con->evalStress(elemIndex, pt, X, e, sens);
      escale = 2.0 * alpha * dfdq[0];
    }

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(escale, sens, u0x, u1x, du0x, du1x, de0ty);

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
  } else if (quantityType == TACS_TEMPERATURE) {
    TacsScalar scale = alpha * dfdq[0];
    basis::template addInterpFieldsTranspose<vars_per_node, 1>(
        pt, &scale, &dfdu[thermal_dof]);
  }
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSThermalShellElement<quadrature, basis, director, model>::getOutputData(
    int elemIndex, ElementType etype, int write_flag, const TacsScalar Xpts[],
    const TacsScalar vars[], const TacsScalar dvars[],
    const TacsScalar ddvars[], int ld_data, TacsScalar *data) {
  // Get the number of nodes associated with the visualization
  int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

  // Compute the node normal directions
  TacsScalar fn[3 * num_nodes], Xdn[9 * num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes];

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
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

    // Evaluate the temperature and temperature gradient
    TacsScalar t;
    basis::template interpFields<vars_per_node, 1>(pt, &vars[thermal_dof], &t);

    // Compute the thermal strain
    TacsScalar eth[9];
    con->evalThermalStrain(elemIndex, pt, X, t, eth);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for (int i = 0; i < 9; i++) {
      em[i] = e[i] - eth[i];
    }

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, em, s);

    if (etype == TACS_BEAM_OR_SHELL_ELEMENT) {
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

#endif  // TACS_THERMAL_SHELL_ELEMENT_H
