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

    // Set the components of the transformation matrix
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
  static const int disp_offset = 3;
  static const int vars_per_node = disp_offset + director::NUM_PARAMETERS;

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

  // void addResidual( int elemIndex,
  //                   double time,
  //                   const TacsScalar *Xpts,
  //                   const TacsScalar *vars,
  //                   const TacsScalar *dvars,
  //                   const TacsScalar *ddvars,
  //                   TacsScalar *res );

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
                   TacsScalar *_Te, TacsScalar *_Ue ){
  // Zero the kinetic and potential energies
  TacsScalar Te = 0.0;
  TacsScalar Ue = 0.0;

  // Compute the number of quadrature points
  const int nquad = quadrature::getNumQuadraturePoints();

  // Compute the node normal directions
  TacsScalar fn1[3*basis::NUM_NODES], fn2[3*basis::NUM_NODES];
  getNodeNormals<basis>(Xpts, fn1, fn2);

  // Compute the frame normal and directors at each node
  TacsScalar d1[dsize], d2[dsize], d1dot[dsize], d2dot[dsize];
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(vars, dvars, fn1, d1, d1dot);
  director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(vars, dvars, fn2, d2, d2dot);

  // Set the total number of tying points needed for this element
  TacsScalar ety[basis::NUM_TYING_POINTS];
  model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn1, fn2, vars, d1, d2, ety);

  // Loop over each quadrature point and add the residual contribution
  for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
    // Get the quadrature weight
    double pt[3];
    double weight = quadrature::getQuadraturePoint(quad_index, pt);

    // Compute X, X,xi and the interpolated normal n0
    TacsScalar X[3], Xxi[3], n0[3], T[9];
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
    TacsScalar e[6]; // The components of the strain
    model::evalStrain(u0x, u1x, e0ty, Ct, e);

    // Compute the corresponding stresses
    TacsScalar s[6];
    con->evalStress(elemIndex, pt, X, e, s);

    Ue += 0.5*detXd*(s[0]*e[0] + s[1]*e[1] + s[2]*e[2] +
                     s[3]*e[3] + s[4]*e[4] + s[5]*e[5]);
  }

  *_Te = Te;
  *_Ue = Ue;
}

/*
  Add the residual to the provided vector
*/
// template <class quadrature, class basis, class director, class model>
// void TACSBeamElement<quadrature, basis, director, model>::
//   addResidual( int elemIndex,
//                double time,
//                const TacsScalar *Xpts,
//                const TacsScalar *vars,
//                const TacsScalar *dvars,
//                const TacsScalar *ddvars,
//                TacsScalar *res ){
//   // Compute the number of quadrature points
//   const int nquad = quadrature::getNumQuadraturePoints();

//   // Derivative of the director field and matrix at each point
//   TacsScalar dd[dsize], dC[csize];
//   memset(dd, 0, 3*basis::NUM_NODES*sizeof(TacsScalar));
//   memset(dC, 0, 9*basis::NUM_NODES*sizeof(TacsScalar));

//   // Zero the contributions to the tying strain derivatives
//   TacsScalar dety[basis::NUM_TYING_POINTS];
//   memset(dety, 0, basis::NUM_TYING_POINTS*sizeof(TacsScalar));

//   // Compute the node normal directions
//   TacsScalar fn[3*basis::NUM_NODES];
//   getNodeNormals<basis>(Xpts, fn);

//   // Compute the frame normal and directors at each node
//   TacsScalar C[csize], d[dsize], ddot[dsize], dddot[dsize];
//   director::template computeDirectorRates<vars_per_node, disp_offset, basis::NUM_NODES>(
//     vars, dvars, ddvars, fn, C, d, ddot, dddot);

//   // Set the total number of tying points needed for this element
//   TacsScalar ety[basis::NUM_TYING_POINTS];
//   model::template computeTyingStrain<vars_per_node, basis>(Xpts, fn, vars, d, ety);

//   // Loop over each quadrature point and add the residual contribution
//   for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
//     // Get the quadrature weight
//     double pt[3];
//     double weight = quadrature::getQuadraturePoint(quad_index, pt);

//     // Compute X, X,xi and the interpolated normal n0
//     TacsScalar X[3], Xxi[6], n0[3], T[9];
//     basis::template interpFields<3, 3>(pt, Xpts, X);
//     basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);
//     basis::template interpFields<3, 3>(pt, fn, n0);

//     // Compute the transformation at the quadrature point
//     transform->computeTransform(Xxi, n0, T);

//     // Evaluate the displacement gradient at the point
//     TacsScalar XdinvT[9], XdinvzT[9];
//     TacsScalar u0x[9], u1x[9], Ct[9];
//     TacsScalar detXd =
//       computeDispGrad<vars_per_node, basis>(pt, Xpts, vars, fn, C, d, Xxi, n0, T,
//                                             XdinvT, XdinvzT, u0x, u1x, Ct);
//     detXd *= weight;

//     // Evaluate the tying components of the strain
//     TacsScalar gty[6]; // The symmetric components of the tying strain
//     interpTyingStrain<basis>(pt, ety, gty);

//     // Compute the symmetric parts of the tying strain
//     TacsScalar e0ty[6]; // e0ty = XdinvT^{T}*gty*XdinvT
//     mat3x3SymmTransformTranspose(XdinvT, gty, e0ty);

//     // Compute the set of strain components
//     TacsScalar e[9]; // The components of the strain
//     model::evalStrain(u0x, u1x, e0ty, Ct, e);

//     // Compute the corresponding stresses
//     TacsScalar s[9];
//     con->evalStress(elemIndex, pt, X, e, s);

//     // Compute the derivative of the product of the stress and strain
//     // with respect to u0x, u1x and e0ty
//     TacsScalar du0x[9], du1x[9], de0ty[6], dCt[9];
//     model::evalStrainSens(detXd, s, u0x, u1x, Ct, du0x, du1x, de0ty, dCt);

//     // Add the contributions to the residual from du0x, du1x and dCt
//     addDispGradSens<vars_per_node, basis>(pt, T, XdinvT, XdinvzT,
//                                           du0x, du1x, dCt, res, dd, dC);

//     // Compute the of the tying strain w.r.t. derivative w.r.t. the coefficients
//     TacsScalar dgty[6];
//     mat3x3SymmTransformTransSens(XdinvT, de0ty, dgty);

//     // Evaluate the tying strain
//     addInterpTyingStrainTranspose<basis>(pt, dgty, dety);
//   }

//   // Set the total number of tying points needed for this element
//   model::template addComputeTyingStrainTranspose<vars_per_node, basis>(
//     Xpts, fn, vars, d, dety, res, dd);

//   // Add the contributions to the director field
//   director::template addDirectorResidual<vars_per_node, disp_offset, basis::NUM_NODES>(
//     vars, dvars, ddvars, fn, dC, dd, res);
// }

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
