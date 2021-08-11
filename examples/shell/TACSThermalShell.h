#ifndef TACS_THERMAL_SHELL_ELEMENT_H
#define TACS_THERMAL_SHELL_ELEMENT_H

#include "TACSShellElementModel.h"
#include "TACSShellElementBasis.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSShellUtilities.h"
#include "TACSShellConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"
#include "TACSDirector.h"
#include "TACSShellElementTransform.h"
#include "TACSElementVerification.h"

template <class quadrature, class basis, class director, class model>
class TACSThermalShellElement : public TACSElement {
 public:
  // Offset within the solution vector to the roational
  // parametrization defined via the director class. Here the offset
  // is 4 corresponding to the (u, v, w) displacements of the
  // mid-surface of the shell and the shell temperature.
  static const int offset = 4;

  // The number of variables defined at each node of the shell
  // element.  There are 3 mid-surface displacements, the temperature
  // plus however many parameters are defined by the director class
  // for the parametrization.
  static const int vars_per_node = offset + director::NUM_PARAMETERS;

  // The number of nodes for this element. This is derived from the
  // basis function class. This is just a handy re-definition since
  // this constant is used in many locations within the element.
  static const int num_nodes = basis::NUM_NODES;

  TACSThermalShellElement( TACSShellTransform *_transform,
                           TACSShellConstitutive *_con ){
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  int getVarsPerNode(){
    return vars_per_node;
  }
  int getNumNodes(){
    return num_nodes;
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
    return quadrature::getNumElementFaces();
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
  static const int usize = 3*num_nodes;
  static const int dsize = 3*num_nodes;
  static const int csize = 9*num_nodes;

  void computeDrillStrain( const TacsScalar Xdn[],
                           const TacsScalar fn[],
                           const TacsScalar vars[],
                           TacsScalar XdinvTn[],
                           TacsScalar Tn[],
                           TacsScalar u0xn[],
                           TacsScalar Ctn[],
                           TacsScalar etn[] );

  TACSShellTransform *transform;
  TACSShellConstitutive *con;
};

/*
  Compute the kinetic and potential energies of the shell
*/
template <class quadrature, class basis, class director, class model>
void TACSThermalShell<quadrature, basis, director, model>::
  computeEnergies( int elemIndex,
                   double time,
                   const TacsScalar *Xpts,
                   const TacsScalar *vars,
                   const TacsScalar *dvars,
                   TacsScalar *_Te, TacsScalar *_Ue ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3*num_nodes], Xdn[9*num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes];

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar XdinvTn[9*num_nodes], Tn[9*num_nodes];
  TacsScalar u0xn[9*num_nodes], Ctn[csize];
  computeDrillStrain(Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  // Compute the director rates
  TacsScalar d[dsize], ddot[dsize];
  director::template
    computeDirectorRates<vars_per_node, offset, num_nodes>(vars, dvars, fn, d, ddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template
    computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
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
    TacsScalar detXd =
      TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, d, Xxi, n0, T,
                                                     XdinvT, XdinvzT, u0x, u1x);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = et;

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    Ue += 0.5*detXd*(s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                     s[3]*e[3] + s[4]*e[4] + s[5]*e[5] +
                     s[6]*e[6] + s[7]*e[7] + s[8]*e[8]);

    // Evaluate the mass moments
    TacsScalar moments[3];
    con->evalMassMoments(elemIndex, pt, X, moments);

    // Compute the velocities and the director velocities
    TacsScalar u0dot[3], d0dot[3];
    basis::template interpFields<vars_per_node, 3>(pt, dvars, u0dot);
    basis::template interpFields<3, 3>(pt, ddot, d0dot);

    Te += 0.5*detXd*(moments[0]*vec3Dot(u0dot, u0dot) +
                     2.0*moments[1]*vec3Dot(u0dot, d0dot) +
                     moments[2]*vec3Dot(d0dot, d0dot));
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSThermalShell<quadrature, basis, director, model>::
  addResidual( int elemIndex,
               double time,
               const TacsScalar *Xpts,
               const TacsScalar *vars,
               const TacsScalar *dvars,
               const TacsScalar *ddvars,
               TacsScalar *res ){

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Derivative of the director field and matrix at each point
  TacsScalar dd[dsize], dTdot[dsize];
  memset(dd, 0, 3*num_nodes*sizeof(TacsScalar));
  memset(dTdot, 0, 3*num_nodes*sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3*num_nodes], Xdn[9*num_nodes];
  TacsShellComputeNodeNormals<basis>(Xpts, fn, Xdn);

  // Compute the drill strain penalty at each node
  TacsScalar etn[num_nodes], detn[num_nodes];
  memset(detn, 0, num_nodes*sizeof(TacsScalar));

  // Store information about the transformation and derivatives at each node for
  // the drilling degrees of freedom
  TacsScalar XdinvTn[9*num_nodes], Tn[9*num_nodes];
  TacsScalar u0xn[9*num_nodes], Ctn[csize];
  computeDrillStrain(Xdn, fn, vars, XdinvTn, Tn, u0xn, Ctn, etn);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template
    computeDirectorRates<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn,
                                                           d, ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
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
    TacsScalar detXd =
      TacsShellComputeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, d, Xxi, n0, T,
                                                     XdinvT, XdinvzT, u0x, u1x);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    basis::interpTyingStrain(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, e);
    e[8] = et;

    // Evaluate the temperature and temperature gradient
    TacsScalar t, txi[2];
    basis::template interpFields<vars_per_node, 1>(pt, &vars[3], &t);
    basis::template interpFieldsGrad<vars_per_node, 1>(pt, &vars[3], txi);

    // Transform to the local component of the heat flux
    TacsScalar tx[2]; // tx = txi*Xdinv*T
    tx[0] = XdinvT[0]*txi[0] + XdinvT[1]*txi[1];
    tx[1] = XdinvT[3]*txi[0] + XdinvT[4]*txi[1];

    // Compute the heat flux
    TacsScalar q[2];
    con->evalHeatFlux(elemIndex, pt, X, tx, q);

    TacsScalar qxi[2];
    qxi[0] = detXd*(XdinvT[0]*q[0] + XdinvT[3]*q[1]);
    qxi[1] = detXd*(XdinvT[1]*q[0] + XdinvT[4]*q[1]);
    basis::addInterpFieldsGradTranspose(pt, 1, qxi, vars_per_node, &res[3]);

    // Compute the thermal strain
    TacsScalar et[9];
    con->evalThermalStrain(elemIndex, pt, X, t, et);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for ( int i = 0; i < 9; i++ ){
      em[i] = e[i] - et[i];
    }

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, em, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6];
    model::evalStrainSens(detXd, s, u0x, u1x, du0x, du1x, de0ty);

    // Add the contribution to the drilling strain
    TacsScalar det = detXd*s[8];
    basis::template addInterpFieldsTranspose<1, 1>(pt, &det, detn);

    // Add the contributions to the residual from du0x, du1x and dCt
    TacsShellAddDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                                   du0x, du1x, res, dd);

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
    du0dot[0] = detXd*(moments[0]*u0ddot[0] + moments[1]*d0ddot[0]);
    du0dot[1] = detXd*(moments[0]*u0ddot[1] + moments[1]*d0ddot[1]);
    du0dot[2] = detXd*(moments[0]*u0ddot[2] + moments[1]*d0ddot[2]);
    basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, du0dot, res);

    TacsScalar dd0dot[3];
    dd0dot[0] = detXd*(moments[1]*u0ddot[0] + moments[2]*d0ddot[0]);
    dd0dot[1] = detXd*(moments[1]*u0ddot[1] + moments[2]*d0ddot[1]);
    dd0dot[2] = detXd*(moments[1]*u0ddot[2] + moments[2]*d0ddot[2]);
    basis::template addInterpFieldsTranspose<3, 3>(pt, dd0dot, dTdot);
  }

  for ( int i = 0; i < num_nodes; i++ ){
    double pt[2];
    basis::getNodePoint(i, pt);

    TacsScalar du0x[9], dCt[9];
    model::evalDrillStrainSens(detn[i], &u0xn[9*i], &Ctn[9*i], du0x, dCt);

    // Compute dCpt = T*dCt*T^{T}
    TacsScalar dCpt[9], tmp[9];
    mat3x3MatMult(&Tn[9*i], dCt, tmp);
    mat3x3MatTransMult(tmp, &Tn[9*i], dCpt);

    director::template
      addRotationMatResidual<vars_per_node, offset, 1>(&vars[i*vars_per_node],
                                                       dCpt, &res[i*vars_per_node]);

    // Compute du0d = T*du0x*XdinvT^{T} + T*du1x*XdinvzT^{T}
    TacsScalar du0d[9];
    mat3x3MatTransMult(du0x, &XdinvTn[9*i], tmp);
    mat3x3MatMult(&Tn[9*i], tmp, du0d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6];
    TacsShellExtractFrame(du0d, du0xi);

    // Compute the gradient of the displacement solution at the quadrature points
    basis::template
      addInterpFieldsGradTranspose<vars_per_node, 3>(pt, du0xi, res);
  }

  // Set the total number of tying points needed for this element
  model::template
    addComputeTyingStrainTranspose<vars_per_node, basis>(Xpts, fn, vars,
                                                         d, dety, res, dd);

  // Add the contributions to the director field
  director::template
    addDirectorResidual<vars_per_node, offset, num_nodes>(vars, dvars, ddvars, fn,
                                                          dTdot, dd, res);

  // Add the contribution from the rotation constraint (defined by the
  // rotational parametrization) - if any
  director::template
    addRotationConstraint<vars_per_node, offset, num_nodes>(vars, res);
}

/*
  Add the residual to the provided vector
*/
/*
template <class quadrature, class basis, class director, class model>
void TACSThermalShell<quadrature, basis, director, model>::
  addJacobian( int elemIndex, double time,
               TacsScalar alpha,
               TacsScalar beta,
               TacsScalar gamma,
               const TacsScalar Xpts[],
               const TacsScalar vars[],
               const TacsScalar dvars[],
               const TacsScalar ddvars[],
               TacsScalar res[],
               TacsScalar mat[] ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 3 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar C[9*basis::NUM_NODES];
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  TacsScalar dddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &ddvars[offset], &fn[3*i], &C[9*i],
                                   &d[3*i], &ddot[3*i], &dddot[3*i]);
  }

  // Derivative of the director field
  TacsScalar dd[3*basis::NUM_NODES];
  memset(dd, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));

  TacsScalar dC[9*basis::NUM_NODES];
  memset(dC, 0, 9*basis::NUM_NODES*sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Evaluate the displacement gradient at the point
    TacsScalar X[3], T[9];
    TacsScalar XdinvT[9], negXdinvXdz[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    TacsScalar detXd = computeDispGrad(pt, Xpts, vars, fn, C, d,
                                       X, T, XdinvT, negXdinvXdz,
                                       u0x, u1x, Ct);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
    TacsScalar d2u0x[45], d2u1x[45], d2e0ty[21], d2Ct[45];
    TacsScalar d2u0x_u1x[81], d2u0x_d2Ct[81], d2u0x_e0ty[81], d2u1x_e0ty[54];
    model::evalStrainSens(detXd, s, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

    d2e0ty[21];

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    model::template addInterpTyingStrainTranspose<basis>(pt, dgty, dety);

    // Compute dCpt = T*dC*T^{T}
    TacsScalar dCpt[9], tmp[9];
    mat3x3MatMult(T, dCt, tmp);
    mat3x3MatTransMult(tmp, T, dCpt);
    basis::addInterpFieldsTranspose(pt, 9, dCpt, 9, dC);

    // Compute du0d = T*(du0x*XdinvT^{T} + du1x*XdinvT^{T}*negXdinvXdz^{T})
    TacsScalar du0d[9];
    mat3x3MatTransMult(du1x, XdinvT, du0d);
    mat3x3MatTransMult(du0d, negXdinvXdz, tmp);
    mat3x3MatTransMultAdd(du0x, XdinvT, tmp);
    mat3x3MatMult(T, tmp, du0d);

    // Compute du1d = T*du1x*XdinvT^{T}
    TacsScalar du1d[9];
    mat3x3MatTransMult(du1x, XdinvT, tmp);
    mat3x3MatMult(T, tmp, du1d);

    // du0d = [du0xi; dd0]
    TacsScalar du0xi[6], dd0[3];
    extractFrame(du0d, du0xi, dd0);

    TacsScalar dd0xi[6];
    extractFrame(du1d, dd0xi);

    // Compute the director field and the gradient of the director
    // field at the specified point
    basis::addInterpFieldsTranspose(pt, 3, dd0, 3, dd);
    basis::addInterpFieldsGradTranspose(pt, 3, dd0xi, 3, dd);

    // Compute the gradient of the displacement solution at the quadrature points
    basis::addInterpFieldsGradTranspose(pt, 3, du0xi, vars_per_node, res);

    TacsScalar d2u0xi[21]; // 21 coefficients
    for ( int j = 0, index = 0; j < 6; j++ ){
      for ( int i = 0; i <= j; i++, index++ ){
        basis::addInterpGradGradOuterProduct(pt, jac, row_incr, col_incr, mat);
      }
    }


    basis::addInterpGradGradOuterProduct(pt, 3, d2u0xi, vars_per_node, )








  }

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<basis>(Xpts, fn, vars_per_node, vars,
                                                        d, dety, res, dd);

  // Add the contributions to the director field
  for ( int i = 0, offset = 3; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::addDirectorJacobian(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &dC[9*i], &dd[3*i], &res[offset]);
  }
}
*/

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSThermalShell<quadrature, basis, director, model>::
  getOutputData( int elemIndex,
                 ElementType etype,
                 int write_flag,
                 const TacsScalar Xpts[],
                 const TacsScalar vars[],
                 const TacsScalar dvars[],
                 const TacsScalar ddvars[],
                 int ld_data,
                 TacsScalar *data ){
  // Get the number of nodes associated with the visualization
  int num_vis_nodes = TacsGetNumVisNodes(basis::getLayoutType());

  // Compute the number of quadrature points
  const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar C[9*basis::NUM_NODES];
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 4; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset], &fn[3*i],
                                   &C[9*i], &d[3*i], &ddot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int index = 0; index < num_vis_nodes; index++ ){
    // Get the quadrature weight
    double pt[3];
    basis::getNodePoint(index, pt);

    // Evaluate the displacement gradient at the point
    TacsScalar X[3], T[9];
    TacsScalar XdinvT[9], negXdinvXdz[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    computeDispGrad(pt, Xpts, vars, fn, C, d,
                    X, T, XdinvT, negXdinvXdz,
                    u0x, u1x, Ct);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    model::template interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Evaluate the temperature and temperature gradient
    TacsScalar t;
    basis::interpFields(pt, vars_per_node, &vars[3], 1, &t);

    // Compute the thermal strain
    TacsScalar et[9];
    con->evalThermalStrain(elemIndex, pt, X, t, et);

    // Compute the mechanical strain (and stress)
    TacsScalar em[9];
    for ( int i = 0; i < 9; i++ ){
      em[i] = e[i] - et[i];
    }

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, em, s);

    if (etype == TACS_BEAM_OR_SHELL_ELEMENT){
      if (write_flag & TACS_OUTPUT_NODES){
        data[0] = X[0];
        data[1] = X[1];
        data[2] = X[2];
        data += 3;
      }
      if (write_flag & TACS_OUTPUT_DISPLACEMENTS){
        int len = vars_per_node;
        if (len > 6){
          len = 6;
        }
        for ( int i = 0; i < len; i++ ){
          data[i] = vars[i + vars_per_node*index];
        }
        for ( int i = len; i < 6; i++ ){
          data[i] = 0.0;
        }
        data += 6;
      }
      if (write_flag & TACS_OUTPUT_STRAINS){
        for ( int i = 0; i < 9; i++ ){
          data[i] = e[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_STRESSES){
        for ( int i = 0; i < 9; i++ ){
          data[i] = s[i];
        }
        data += 9;
      }
      if (write_flag & TACS_OUTPUT_EXTRAS){
        data[0] = con->evalFailure(elemIndex, pt, X, e);
        data[1] = con->evalDesignFieldValue(elemIndex, pt, X, 0);
        data[2] = con->evalDesignFieldValue(elemIndex, pt, X, 1);
        data[3] = con->evalDesignFieldValue(elemIndex, pt, X, 2);
        data += 4;
      }
    }
  }
}

#endif // TACS_THERMAL_SHELL_ELEMENT_H
