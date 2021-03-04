#ifndef TACS_BEAM_ELEMENT_H
#define TACS_BEAM_ELEMENT_H

#include "TACSBeamElementModel.h"
#include "TACSBeamElementBasis.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSBeamConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"

/*
  Compute the transformation from the local coordinates
  to
*/
class TACSBeamTransform : public TACSObject {
 public:
  /*
    Given the local beam element reference frame Xf, compute the
    transformation from the global coordinates to the shell-aligned local axis.
  */
  virtual void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ) = 0;
};

class TACSBeamRefAxisTransform : public TACSShellTransform {
 public:
  TACSBeamRefAxisTransform( const TacsScalar _axis[] ){
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];

    TacsScalar norm = sqrt(vec3Dot(axis, axis));
    TacsScalar invNorm = 0.0;
    if (norm != 0.0){
      invNorm = 1.0/norm;
    }
    vec3Scale(invNorm, axis);
  }

  void computeTransform( const TacsScalar Xxi[], TacsScalar T[] ){
    // Compute the reference frame
    TacsScalar t[3];
    TacsScalar tnorm = sqrt(vec3Dot(Xxi, Xxi));
    TacsScalar tinv = 1.0/tnorm;
    t[0] = tinv*Xa[0];
    t[1] = tinv*Xa[1];
    t[2] = tinv*Xa[2];

    // Compute the first direction in the plane
    TacsScalar n1[3];
    TacsScalar tdot = vec3Dot(t, axis);
    n1[0] = axis[0] - tdot*t[0];
    n1[1] = axis[1] - tdot*t[1];
    n1[2] = axis[2] - tdot*t[2];

    // Compute the norm
    TacsScalar n1inv = 1.0/sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
    n1[0] *= n1inv;
    n1[1] *= n1inv;
    n1[2] *= n1inv;

    // Compute the cross product
    TacsScalar n2[3];
    crossProduct(1.0, t, n1, n2);

    // Set the components of the transformation matrix    // Set the components of the transformation
    T[0] = t[0];
    T[3] = t[1];
    T[6] = t[2];

    T[1] = n1[0];
    T[4] = n1[1];
    T[7] = n1[2];

    T[2] = n2[0];
    T[5] = n2[1];
    T[8] = n2[2];
  }

 private:
  TacsScalar axis[3];
};


template <class quadrature, class basis, class director, class model>
class TACSBeamElement : public TACSElement {
 public:
  TACSBeamElement( TACSBeamTransform *_transform,
                   TACSBeamConstitutive *_con ){
    transform = _transform;
    transform->incref();

    con = _con;
    con->incref();
  }

  int getVarsPerNode(){
    return 4 + director::NUM_PARAMETERS;
  }
  int getNumNodes(){
    return basis::NUM_NODES;
  }

  ElementLayout getLayoutType(){
    return basis::getLayoutType();
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

/*
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
*/

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
  TACSBeamTransform *transform;
  TACSBeamConstitutive *con;

  void getNodeNormals( const TacsScalar Xpts[],
                       TacsScalar fn[],
                       TacsScalar fnorm[]=NULL );

  TacsScalar computeDispGrad( const double pt[],
                              const TacsScalar Xpts[],
                              const TacsScalar vars[],
                              const TacsScalar fn[],
                              const TacsScalar d[],
                              const TacsScalar C[],
                              TacsScalar X[],
                              TacsScalar T[],
                              TacsScalar XdinvT[],
                              TacsScalar negXdinvXdz[],
                              TacsScalar u0x[],
                              TacsScalar u1x[],
                              TacsScalar Ct[] );

  static void assembleFrame( const TacsScalar Xxi[],
                             const TacsScalar n[],
                             TacsScalar Xd[] ){
    Xd[0] = Xxi[0];
    Xd[1] = Xxi[1];
    Xd[2] = n[0];

    Xd[3] = Xxi[2];
    Xd[4] = Xxi[3];
    Xd[5] = n[1];

    Xd[6] = Xxi[4];
    Xd[7] = Xxi[5];
    Xd[8] = n[2];
  }

  static void assembleFrame( const TacsScalar nxi[],
                             TacsScalar Xdz[] ){
    Xdz[0] = nxi[0];
    Xdz[1] = nxi[1];
    Xdz[2] = 0.0;

    Xdz[3] = nxi[2];
    Xdz[4] = nxi[3];
    Xdz[5] = 0.0;

    Xdz[6] = nxi[4];
    Xdz[7] = nxi[5];
    Xdz[8] = 0.0;
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[],
                            TacsScalar n[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];
    n[0] = Xd[2];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];
    n[1] = Xd[5];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
    n[2] = Xd[8];
  }

  static void extractFrame( const TacsScalar Xd[],
                            TacsScalar Xxi[] ){
    Xxi[0] = Xd[0];
    Xxi[1] = Xd[1];

    Xxi[2] = Xd[3];
    Xxi[3] = Xd[4];

    Xxi[4] = Xd[6];
    Xxi[5] = Xd[7];
  }
};

template <class quadrature, class basis, class director, class model>
void TACSBeamElement<quadrature, basis, director, model>::
  getNodeNormals( const TacsScalar Xpts[],
                  TacsScalar fn[],
                  TacsScalar fnorm[] ){
  for ( int i = 0; i < basis::NUM_NODES; i++ ){
    double pt[2];
    basis::getNodePoint(i, pt);

    // Compute the derivative X,xi at each node
    TacsScalar Xxi[6];
    basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);

    TacsScalar a[3], b[3];
    a[0] = Xxi[0];
    a[1] = Xxi[2];
    a[2] = Xxi[4];

    b[0] = Xxi[1];
    b[1] = Xxi[3];
    b[2] = Xxi[5];

    // Compute the normal direction at the point
    crossProduct(a, b, &fn[3*i]);

    // Compute the 2-norm of the vector in the normal direction
    TacsScalar norm = sqrt(vec3Dot(&fn[3*i], &fn[3*i]));

    // Save the 2-norm value if the fnorm argument is not NULL
    if (fnorm){
      fnorm[i] = norm;
    }

    // Scale the normal direction
    if (norm != 0.0){
      vec3Scale(1.0/norm, &fn[3*i]);
    }
  }
}

/*
  Compute the displacement gradient of the constant and through-thickness
  rate of change of the displacements.
*/
template <class quadrature, class basis, class director, class model>
TacsScalar TACSBeamElement<quadrature, basis, director, model>::
  computeDispGrad( const double pt[],
                   const TacsScalar Xpts[],
                   const TacsScalar vars[],
                   const TacsScalar fn[],
                   const TacsScalar C[],
                   const TacsScalar d[],
                   TacsScalar X[],
                   TacsScalar T[],
                   TacsScalar XdinvT[],
                   TacsScalar negXdinvXdz[],
                   TacsScalar u0x[],
                   TacsScalar u1x[],
                   TacsScalar Ct[] ){
  const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // Compute X,xi = [dX/dxi1 ; dX/dxi2] and n,xi = [dn/dxi1; dn/dxi2]
  TacsScalar Xxi[6], n[3], nxi[6];
  basis::interpFields(pt, 3, Xpts, 3, X);
  basis::interpFields(pt, 3, fn, 3, n);
  basis::interpFieldsGrad(pt, 3, Xpts, 3, Xxi);
  basis::interpFieldsGrad(pt, 3, fn, 3, nxi);

  // Compute the transformation at the quadrature point
  transform->computeTransform(Xxi, T);

  // Assemble the terms Xd = [Xxi; n] and Xdz
  TacsScalar Xd[9], Xdz[9];
  assembleFrame(Xxi, n, Xd);
  assembleFrame(nxi, Xdz);

  // Compute the inverse of the 3x3 Jacobian transformation
  TacsScalar Xdinv[9];
  TacsScalar detXd = inv3x3(Xd, Xdinv);

  // Compute negXdinvXdz = -Xdinv*Xdz
  mat3x3MatMult(Xdinv, Xdz, negXdinvXdz);
  for ( int i = 0; i < 9; i++ ){
    negXdinvXdz[i] *= -1.0;
  }

  // Compute XdinvT = Xdinv*T
  mat3x3MatMult(Xdinv, T, XdinvT);

  // Compute the director field and the gradient of the director
  // field at the specified point
  TacsScalar d0[3], d0xi[6];
  basis::interpFields(pt, 3, d, 3, d0);
  basis::interpFieldsGrad(pt, 3, d, 3, d0xi);

  // Compute the gradient of the displacement solution at the quadrature points
  TacsScalar u0xi[6];
  basis::interpFieldsGrad(pt, vars_per_node, vars, 3, u0xi);

  // Input: u0xi, d0, d0xi, T, negXdinvXdz, XdinvT
  // Output: u0x, u1x

  // Set u0x = [u0,xi ; d]
  assembleFrame(u0xi, d0, u0x); // Use u0x to store [u0,xi; d0]

  // u1x = T^{T}*(u0d*(-Xdinv*Xdz) + u1d)*Xdinv*T
  TacsScalar tmp[9];
  assembleFrame(d0xi, u1x); // Use u1x to store [d0,xi; 0]
  mat3x3MatMultAdd(u0x, negXdinvXdz, u1x);
  mat3x3TransMatMult(T, u1x, tmp);
  mat3x3MatMult(tmp, XdinvT, u1x);

  // Compute the transformation u0x = T^{T}*ueta*Xdinv*T
  // u0x = T^{T}*u0d*Xdinv*T
  mat3x3MatMult(u0x, XdinvT, tmp);
  mat3x3TransMatMult(T, tmp, u0x);

  // Compute the interpolation of the entries of the C matrix
  TacsScalar Cpt[9];
  basis::interpFields(pt, 9, C, 9, Cpt);

  // Compute Ct = T^{T}*Cpt*T
  mat3x3TransMatMult(T, Cpt, tmp);
  mat3x3MatMult(tmp, T, Ct);

  return detXd;
}

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
                   TacsScalar *_Te, TacsScalar *_Ue ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();
  const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn1[3*basis::NUM_NODES], fn2[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[3*basis::NUM_NODES];
  TacsScalar d1dot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 4; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::computeDirectorRates(&vars[offset], &dvars[offset],
                                   &fn1[3*i], &d1[3*i], &d1dot[3*i]);
  }

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<basis>(Xpts, fn, vars_per_node, vars, d, ety);

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
  const int vars_per_node = 4 + director::NUM_PARAMETERS;

  // Compute the node normal directions
  TacsScalar fn[3*basis::NUM_NODES];
  getNodeNormals(Xpts, fn);

  // Compute the frame normal and directors at each node
  TacsScalar C[9*basis::NUM_NODES];
  TacsScalar d[3*basis::NUM_NODES];
  TacsScalar ddot[3*basis::NUM_NODES];
  TacsScalar dddot[3*basis::NUM_NODES];
  for ( int i = 0, offset = 4; i < basis::NUM_NODES; i++, offset += vars_per_node ){
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

  // Zero the contributions to the
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

    // Evaluate the temperature and temperature gradient
    TacsScalar t, txi[2];
    basis::interpFields(pt, vars_per_node, &vars[3], 1, &t);
    basis::interpFieldsGrad(pt, vars_per_node, &vars[3], 1, txi);

    // Transform to the local component of the heat flux
    TacsScalar tx[2]; // tx = txi*Xdinv*T
    tx[0] = XdinvT[0]*txi[0] + XdinvT[1]*txi[1];
    tx[1] = XdinvT[3]*txi[0] + XdinvT[4]*txi[1];

    TacsScalar q[2];
    con->evalHeatFlux(elemIndex, pt, X, tx, q);

    TacsScalar qxi[2];
    qxi[0] = detXd*(XdinvT[0]*q[0] + XdinvT[3]*q[1]);
    qxi[1] = detXd*(XdinvT[1]*q[0] + XdinvT[4]*q[1]);
    basis::addInterpFieldsGradTranspose(pt, 1, qxi, vars_per_node, &res[3]);

    // This is adding a constant heat source to the element. This
    // is only for testing purposes and should be removed eventually
    TacsScalar rhs = -10.0*detXd;
    basis::addInterpFieldsTranspose(pt, 1, &rhs, vars_per_node, &res[3]);

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
    TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
    model::evalStrainSens(detXd, s, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

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
  }

  // Set the total number of tying points needed for this element
  model::template addComputeTyingStrainTranspose<basis>(Xpts, fn, vars_per_node, vars,
                                                        d, dety, res, dd);

  // Add the contributions to the director field
  for ( int i = 0, offset = 4; i < basis::NUM_NODES; i++, offset += vars_per_node ){
    director::addDirectorResidual(&vars[offset], &dvars[offset], &ddvars[offset],
                                   &fn[3*i], &dC[9*i], &dd[3*i], &res[offset]);
  }
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

#endif // TACS_SHELL_ELEMENT_H
