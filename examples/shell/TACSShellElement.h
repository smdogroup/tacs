#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

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
class TACSShellElement : public TACSElement {
 public:
  static const int disp_offset = 3;
  static const int vars_per_node = disp_offset + director::NUM_PARAMETERS;

  TACSShellElement( TACSShellTransform *_transform,
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

  void addJacobian( int elemIndex, double time,
                    TacsScalar alpha,
                    TacsScalar beta,
                    TacsScalar gamma,
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[],
                    TacsScalar res[],
                    TacsScalar mat[] );

  void addAdjResProduct( int elemIndex, double time,
                         TacsScalar scale,
                         const TacsScalar psi[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[],
                         int dvLen,
                         TacsScalar dfdx[] );

  int evalPointQuantity( int elemIndex, int quantityType,
                         double time,
                         int n, double pt[],
                         const TacsScalar Xpts[],
                         const TacsScalar vars[],
                         const TacsScalar dvars[],
                         const TacsScalar ddvars[],
                         TacsScalar *detXd,
                         TacsScalar *quantity );

  void addPointQuantityDVSens( int elemIndex, int quantityType,
                               double time,
                               TacsScalar scale,
                               int n, double pt[],
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfdq[],
                               int dvLen,
                               TacsScalar dfdx[] );

  void addPointQuantitySVSens( int elemIndex, int quantityType,
                               double time,
                               TacsScalar alpha,
                               TacsScalar beta,
                               TacsScalar gamma,
                               int n, double pt[],
                               const TacsScalar Xpts[],
                               const TacsScalar vars[],
                               const TacsScalar dvars[],
                               const TacsScalar ddvars[],
                               const TacsScalar dfdq[],
                               TacsScalar dfdu[] );

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

  TACSShellTransform *transform;
  TACSShellConstitutive *con;
};

/*
  Compute the kinetic and potential energies of the shell
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
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
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar C[csize];
  director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

  TacsScalar d[dsize], ddot[dsize];
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, fn, d, ddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    TacsScalar detXd =
      computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                            XdinvT, XdinvzT, u0x, u1x, Ct);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

    Ue += 0.5*detXd*(s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                     s[3]*e[3] + s[4]*e[4] + s[5]*e[5] +
                     s[6]*e[6] + s[7]*e[7] + s[8]*e[8]);
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addResidual( int elemIndex,
               double time,
               const TacsScalar Xpts[],
               const TacsScalar vars[],
               const TacsScalar dvars[],
               const TacsScalar ddvars[],
               TacsScalar res[] ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Derivative of the director field and matrix at each point
  TacsScalar dd[dsize], dC[csize];
  memset(dd, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
  memset(dC, 0, 9*basis::NUM_NODES*sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  // Compute the rotation matrix and directors at each node
  TacsScalar C[csize];
  director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, fn, d, ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    TacsScalar detXd =
      computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                            XdinvT, XdinvzT, u0x, u1x, Ct);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

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
    model::evalStrainSens(detXd, s, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

    // Add the contributions to the residual from du0x, du1x and dCt
    addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                          du0x, du1x, dCt, res, dd, dC);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    addInterpTyingStrainTranspose<basis>(pt, dgty, dety);
  }

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
    Xpts, fn, vars, d, dety, res, dd);

  // Add the contributions to the director field
  director::template
    addRotationMatResidual<vars_per_node, disp_offset, basis::NUM_NODES>(vars, dC, res);
  director::template addDirectorResidual<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, fn, dd, res);
}

/*
  Add the contributions to the residual and Jacobian matrix
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
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

  // Derivative of the director field
  TacsScalar dd[dsize], d2d[dsize*dsize];
  memset(dd, 0, dsize*sizeof(TacsScalar));
  memset(d2d, 0, dsize*dsize*sizeof(TacsScalar));

  // Derivative of the rotation matrix at each point
  TacsScalar dC[csize], d2C[csize*csize];
  memset(dC, 0, csize*sizeof(TacsScalar));
  memset(d2C, 0, csize*csize*sizeof(TacsScalar));

  // Coupling derivatives
  TacsScalar d2du[usize*dsize], d2Cu[usize*csize], d2Cd[csize*dsize];;
  memset(d2du, 0, usize*dsize*sizeof(TacsScalar));
  memset(d2Cu, 0, usize*csize*sizeof(TacsScalar));
  memset(d2Cd, 0, csize*dsize*sizeof(TacsScalar));

  // Zero the contributions to the tying strain derivatives
  TacsScalar dety[basis::NUM_TYING_POINTS];
  TacsScalar d2ety[basis::NUM_TYING_POINTS*basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));
  memset(d2ety, 0, basis::NUM_TYING_POINTS*basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  // Compute the rotation matrix and directors at each node
  TacsScalar C[csize];
  director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, fn, d, ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    TacsScalar detXd =
      computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                            XdinvT, XdinvzT, u0x, u1x, Ct);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

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
    TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
    model::evalStrainSens(detXd, s, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

    TacsScalar d2u0x[81], d2u1x[81], d2u0xu1x[81];
    TacsScalar d2e0ty[36], d2e0tyu0x[54], d2e0tyu1x[54];
    TacsScalar d2Ct[81], d2Ctu0x[81];
    model::evalStrainHessian(alpha*detXd, s, Cs, u0x, u1x, e0ty, Ct,
                             d2u0x, d2u1x, d2u0xu1x,
                             d2e0ty, d2e0tyu0x, d2e0tyu1x,
                             d2Ct, d2Ctu0x);

    // Add the contributions to the residual from du0x, du1x and dCt
    if (res){
      addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                            du0x, du1x, dCt, res, dd, dC);
    }

    // Add the contributions to the residual from du0x, du1x and dCt
    addDispGradHessian<vars_per_node, basis>(
      pt, T, XdinvT, XdinvzT, d2u0x, d2u1x, d2u0xu1x,
      d2Ct, d2Ctu0x, mat, d2d, d2du, d2C, d2Cd, d2Cu);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6], d2gty[36];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
    mat3x3SymmTransformTransHessian(XdinvT, d2e0ty, d2gty);

    // Evaluate the tying strain
    addInterpTyingStrainTranspose<basis>(pt, dgty, dety);
    addInterpTyingStrainHessian<basis>(pt, d2gty, d2ety);
  }

  // Add the residual from the tying strain
  if (res){
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, res, dd);
  }

  // Add the second order terms from the tying strain
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
    Xpts, fn, vars, d, d2ety, mat, d2d, d2du);

  if (res){
    // Add the contributions to the director field
    director::template
      addRotationMatResidual<vars_per_node, disp_offset, basis::NUM_NODES>(vars, dC, res);
    director::template addDirectorResidual<vars_per_node, disp_offset, basis::NUM_NODES>(
      vars, dvars, ddvars, fn, dd, res);
  }

  // Add the contributions to the stiffness matrix
  director::template addDirectorJacobian<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, fn, d2d, d2du, mat);
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addAdjResProduct( int elemIndex, double time,
                    TacsScalar scale,
                    const TacsScalar psi[],
                    const TacsScalar Xpts[],
                    const TacsScalar vars[],
                    const TacsScalar dvars[],
                    const TacsScalar ddvars[],
                    int dvLen,
                    TacsScalar dfdx[] ){
  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  // Compute the frame normal and directors at each node
  // Compute the rotation matrix and directors at each node
  TacsScalar C[csize], Cd[csize];
  director::template computeRotationMatDeriv<vars_per_node, disp_offset, basis::NUM_NODES>(vars, psi, C, Cd);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize], dd[dsize];
  director::template computeDirectorRatesDeriv<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, psi, fn, d, ddot, dddot, dd);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS], etyd[basis::NUM_TYING_POINTS];
  model::template computeTyingStrainDeriv<vars_per_node, basis>(
    Xpts, fn, vars, d, psi, dd, ety, etyd);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    TacsScalar u0xd[9], u1xd[9], Ctd[9];
    TacsScalar detXd =
      computeDispGradDeriv<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                                 psi, dd, Cd, XdinvT, XdinvzT,
                                                 u0x, u1x, Ct, u0xd, u1xd, Ctd);
    detXd *= weight;

    // Evaluate the tying components of the strain
    TacsScalar gty[6], gtyd[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);
    interpTyingStrain<basis>(pt, etyd, gtyd);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6], e0tyd[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);
    mat3x3SymmTransformTranspose(XdinvT, gtyd, e0tyd);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    TacsScalar ed[9]; // The directional derivative components of the strain
    model::evalStrainDeriv(u0x, u1x, e0ty, Ct, u0xd, u1xd, e0tyd, Ctd, e, ed);

    // The directional derivative of the strain along the adjoint direction
    con->addStressDVSens(elemIndex, scale*detXd, pt, X, e, ed, dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
int TACSShellElement<quadrature, basis, director, model>::
  evalPointQuantity( int elemIndex, int quantityType,
                     double time,
                     int n, double pt[],
                     const TacsScalar Xpts[],
                     const TacsScalar vars[],
                     const TacsScalar dvars[],
                     const TacsScalar ddvars[],
                     TacsScalar *detXd,
                     TacsScalar *quantity ){
  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  if (quantityType == TACS_FAILURE_INDEX){
    // Compute the frame normal and directors at each node
    TacsScalar C[csize], d[dsize], ddot[dsize], dddot[dsize];
    director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);
    director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
      vars, dvars, ddvars, fn, d, ddot, dddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    *detXd =
      computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                            XdinvT, XdinvzT, u0x, u1x, Ct);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    *quantity = con->evalFailure(elemIndex, pt, X, e);

    return 1;
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    TacsScalar Xxi[6], n0[3], X[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    TacsScalar Xd[9];
    assembleFrame(Xxi, n0, Xd);
    *detXd = det3x3(Xd);

    *quantity = con->evalDensity(elemIndex, pt, X);

    return 1;
  }

  return 0;
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addPointQuantityDVSens( int elemIndex, int quantityType,
                          double time,
                          TacsScalar scale,
                          int n, double pt[],
                          const TacsScalar Xpts[],
                          const TacsScalar vars[],
                          const TacsScalar dvars[],
                          const TacsScalar ddvars[],
                          const TacsScalar dfdq[],
                          int dvLen,
                          TacsScalar dfdx[] ){
  if (quantityType == TACS_FAILURE_INDEX){
    // Compute the node normal directions
    TacsScalar fn[3*basis::NUM_NODES];
    getNodeNormals<basis>(Xpts, fn);

    // Compute the frame normal and directors at each node
    TacsScalar C[csize];
    director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

    TacsScalar d[dsize], ddot[dsize], dddot[dsize];
    director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
      vars, dvars, ddvars, fn, d, ddot, dddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                          XdinvT, XdinvzT, u0x, u1x, Ct);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    con->addFailureDVSens(elemIndex, scale*dfdq[0], pt, X, e, dvLen, dfdx);
  }
  else if (quantityType == TACS_ELEMENT_DENSITY){
    TacsScalar X[3];
    basis::template interpFields<3, 3>(pt, Xpts, X);

    con->addDensityDVSens(elemIndex, scale*dfdq[0], pt, X, dvLen, dfdx);
  }
}

template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
  addPointQuantitySVSens( int elemIndex, int quantityType,
                          double time,
                          TacsScalar alpha,
                          TacsScalar beta,
                          TacsScalar gamma,
                          int n, double pt[],
                          const TacsScalar Xpts[],
                          const TacsScalar vars[],
                          const TacsScalar dvars[],
                          const TacsScalar ddvars[],
                          const TacsScalar dfdq[],
                          TacsScalar dfdu[] ){
  if (quantityType == TACS_FAILURE_INDEX){
    // Derivative of the director field
    TacsScalar dd[dsize], dC[csize];
    memset(dd, 0, dsize*sizeof(TacsScalar));
    memset(dC, 0, csize*sizeof(TacsScalar));

    // Compute the node normal directions
    TacsScalar fn[3*basis::NUM_NODES];
    getNodeNormals<basis>(Xpts, fn);

    // Compute the frame normal and directors at each node
    TacsScalar C[csize];
    director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

    TacsScalar d[dsize], ddot[dsize], dddot[dsize];
    director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
      vars, dvars, ddvars, fn, d, ddot, dddot);

    // Set the total number of tying points needed for this element
    TacsScalar ety[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

    // Zero the contributions to the
    TacsScalar dety[basis::NUM_TYING_POINTS];
    memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                          XdinvT, XdinvzT, u0x, u1x, Ct);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Compute the sensitivity of the failure index w.r.t. the strain
    TacsScalar sens[9];
    con->evalFailureStrainSens(elemIndex, pt, X, e, sens);

    // Compute the derivative of the product of the stress and strain
    // with respect to u0x, u1x and e0ty
    TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
    model::evalStrainSens(alpha*dfdq[0], sens, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

    addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
                                          du0x, du1x, dCt, dfdu, dd, dC);

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgty[6];
    mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

    // Evaluate the tying strain
    addInterpTyingStrainTranspose<basis>(pt, dgty, dety);

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, d, dety, dfdu, dd);

    // Add the contributions to the director field
    director::template addRotationMatResidual<vars_per_node, disp_offset, basis::NUM_NODES>(vars, dC, dfdu);
    director::template addDirectorResidual<vars_per_node, disp_offset, basis::NUM_NODES>(
      vars, dvars, ddvars, fn, dd, dfdu);
  }
}

/*
  Get the element data for the basis
*/
template <class quadrature, class basis, class director, class model>
void TACSShellElement<quadrature, basis, director, model>::
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

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn);

  // Compute the rotation matrix and directors at each node
  TacsScalar C[csize];
  director::template computeRotationMat<vars_per_node, disp_offset, basis::NUM_NODES>(vars, C);

  TacsScalar d[dsize], ddot[dsize], dddot[dsize];
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
    vars, dvars, ddvars, fn, d, ddot, dddot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int index = 0; index < num_vis_nodes; index++ ){
    // Get the quadrature weight
    double pt[3];
    basis::getNodePoint(index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[6], n0[3], T[9];
    basis::template interpFields<3, 3>(pt, Xpts, X);
    basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
    basis::template interpFields<3, 3>(pt, fn, n0);

    // Compute the transformation at the quadrature point
    transform->computeTransform(Xxi, n0, T);

    // Evaluate the displacement gradient at the point
    TacsScalar XdinvT[9], XdinvzT[9];
    TacsScalar u0x[9], u1x[9], Ct[9];
    computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
                                          XdinvT, XdinvzT, u0x, u1x, Ct);

    // Evaluate the tying components of the strain
    TacsScalar gty[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, ety, gty);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
    mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

    // Compute the set of strain components
    TacsScalar e[9]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Compute the corresponding stresses
    TacsScalar s[9];
    con->evalStress(elemIndex, pt, X, e, s);

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

template <int vars_per_node, class basis, class model>
int TacsTestShellTyingStrain( double dh=1e-7,
                              int test_print_level=2,
                              double test_fail_atol=1e-5,
                              double test_fail_rtol=1e-5 ){
  const int size = vars_per_node*basis::NUM_NODES;
  const int usize = 3*basis::NUM_NODES;
  const int dsize = 3*basis::NUM_NODES;

  TacsScalar Xpts[3*basis::NUM_NODES], fn[3*basis::NUM_NODES];
  TacsGenerateRandomArray(Xpts, 3*basis::NUM_NODES);
  TacsGenerateRandomArray(fn, 3*basis::NUM_NODES);

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
  TacsScalar d2ety[basis::NUM_TYING_POINTS*basis::NUM_TYING_POINTS];
  memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));
  memset(d2ety, 0, basis::NUM_TYING_POINTS*basis::NUM_TYING_POINTS*sizeof(TacsScalar));

  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

  // Evaluate the tying components of the strain
  TacsScalar gty[6]; // The symmetric components of the tying strain
  interpTyingStrain<basis>(pt, ety, gty);

  // Compute the symmetric parts of the tying strain
  TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
  mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

  // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
  TacsScalar dgty[6], d2gty[36];
  mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);
  mat3x3SymmTransformTransHessian(XdinvT, d2e0ty, d2gty);

  // Evaluate the tying strain
  addInterpTyingStrainTranspose<basis>(pt, dgty, dety);
  addInterpTyingStrainHessian<basis>(pt, d2gty, d2ety);

  TacsScalar res[size], dd[dsize];
  memset(res, 0, size*sizeof(TacsScalar));
  memset(dd, 0, dsize*sizeof(TacsScalar));

  TacsScalar mat[size*size], d2d[dsize*dsize], d2du[dsize*usize];
  memset(mat, 0, size*size*sizeof(TacsScalar));
  memset(d2d, 0, dsize*dsize*sizeof(TacsScalar));
  memset(d2du, 0, dsize*usize*sizeof(TacsScalar));

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
    Xpts, fn, vars, d, dety, res, dd);
  model::template addComputeTyingStrainHessian<vars_per_node, basis>(
    Xpts, fn, vars, d, d2ety, mat, d2d, d2du);

  TacsScalar fdmat[size*size], fdd2du[dsize*usize];
  for ( int i = 0; i < size; i++ ){
    TacsScalar varst[size];
    memcpy(varst, vars, size*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    varst[i] = vars[i] + TacsScalar(0.0, dh);
#else
    varst[i] = vars[i] + dh;
#endif // TACS_USE_COMPLEX

    // Perturb the variables
    TacsScalar etyt[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, varst, d, etyt);

    // Evaluate the tying components of the strain
    TacsScalar gtyt[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, etyt, gtyt);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0tyt[6];
    mat3x3SymmTransformTranspose(XdinvT, gtyt, e0tyt);

    TacsScalar de0tyt[6];
    for ( int j = 0; j < 6; j++ ){
      de0tyt[j] = de0ty[j];
      for ( int k = 0; k < 6; k++ ){
        de0tyt[j] += d2e0ty[6*j + k]*(e0tyt[k] - e0ty[k]);
      }
    }

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgtyt[6];
    mat3x3SymmTransformTransSens(XdinvT, de0tyt, dgtyt);

    TacsScalar detyt[basis::NUM_TYING_POINTS];
    memset(detyt, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

    // Evaluate the tying strain
    addInterpTyingStrainTranspose<basis>(pt, dgtyt, detyt);

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, varst, d, detyt, rest, ddt);

    for ( int j = 0; j < size; j++ ){
#ifdef TACS_USE_COMPLEX
      fdmat[size*j + i] = TacsImagPart(rest[j])/dh;
#else
      fdmat[size*j + i] = (rest[j] - res[j])/dh;
#endif // TACS_USE_COMPLEX
    }

    if (i % vars_per_node < 3){
      int index = 3*(i / vars_per_node) + i % vars_per_node;
      for ( int j = 0; j < dsize; j++ ){
#ifdef TACS_USE_COMPLEX
        fdd2du[usize*j + index] = TacsImagPart(ddt[j])/dh;
#else
        fdd2du[usize*j + index] = (ddt[j] - dd[j])/dh;
#endif // TACS_USE_COMPLEX
      }
    }
  }

  int fail = 0;
  double max_err, max_rel;
  int max_err_index, max_rel_index;

  // Compute the error
  max_err = TacsGetMaxError(mat, fdmat, size*size, &max_err_index);
  max_rel = TacsGetMaxRelError(mat, fdmat, size*size, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. vars\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "mat", mat, fdmat, size*size);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  // Compute the error
  max_err = TacsGetMaxError(d2du, fdd2du, dsize*usize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2du, fdd2du, dsize*usize, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. vars and d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2du", d2du, fdd2du, dsize*usize);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);


  TacsScalar fdd2d[dsize*dsize];
  for ( int i = 0; i < dsize; i++ ){
    TacsScalar dt[size];
    memcpy(dt, d, dsize*sizeof(TacsScalar));

#ifdef TACS_USE_COMPLEX
    dt[i] = d[i] + TacsScalar(0.0, dh);
#else
    dt[i] = d[i] + dh;
#endif // TACS_USE_COMPLEX

    // Perturb the variables
    TacsScalar etyt[basis::NUM_TYING_POINTS];
    model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, dt, etyt);

    // Evaluate the tying components of the strain
    TacsScalar gtyt[6]; // The symmetric components of the tying strain
    interpTyingStrain<basis>(pt, etyt, gtyt);

    // Compute the symmetric parts of the tying strain
    TacsScalar e0tyt[6];
    mat3x3SymmTransformTranspose(XdinvT, gtyt, e0tyt);

    TacsScalar de0tyt[6];
    for ( int j = 0; j < 6; j++ ){
      de0tyt[j] = de0ty[j];
      for ( int k = 0; k < 6; k++ ){
        de0tyt[j] += d2e0ty[6*j + k]*(e0tyt[k] - e0ty[k]);
      }
    }

    // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
    TacsScalar dgtyt[6];
    mat3x3SymmTransformTransSens(XdinvT, de0tyt, dgtyt);

    TacsScalar detyt[basis::NUM_TYING_POINTS];
    memset(detyt, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

    // Evaluate the tying strain
    addInterpTyingStrainTranspose<basis>(pt, dgtyt, detyt);

    TacsScalar rest[size], ddt[dsize];
    memset(rest, 0, size*sizeof(TacsScalar));
    memset(ddt, 0, dsize*sizeof(TacsScalar));

    // Set the total number of tying points needed for this element
    model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
      Xpts, fn, vars, dt, detyt, rest, ddt);

    for ( int j = 0; j < dsize; j++ ){
#ifdef TACS_USE_COMPLEX
      fdd2d[dsize*j + i] = TacsImagPart(ddt[j])/dh;
#else
      fdd2d[dsize*j + i] = (ddt[j] - dd[j])/dh;
#endif // TACS_USE_COMPLEX
    }
  }

  // Compute the error
  max_err = TacsGetMaxError(d2d, fdd2d, dsize*dsize, &max_err_index);
  max_rel = TacsGetMaxRelError(d2d, fdd2d, dsize*dsize, &max_rel_index);

  if (test_print_level > 0){
    fprintf(stderr, "Testing the second derivative w.r.t. d\n");
    fprintf(stderr, "Max Err: %10.4e in component %d.\n",
            max_err, max_err_index);
    fprintf(stderr, "Max REr: %10.4e in component %d.\n",
            max_rel, max_rel_index);
  }
  // Print the error if required
  if (test_print_level > 1){
    TacsPrintErrorComponents(stderr, "d2d", d2d, fdd2d, dsize*dsize);
  }
  if (test_print_level){ fprintf(stderr, "\n"); }

  fail = (max_err > test_fail_atol || max_rel > test_fail_rtol);

  return fail;
}

#endif // TACS_SHELL_ELEMENT_H
