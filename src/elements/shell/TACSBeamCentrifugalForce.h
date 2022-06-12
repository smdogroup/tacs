#ifndef TACS_BEAM_CENTRIFUGAL_FORCE_H
#define TACS_BEAM_CENTRIFUGAL_FORCE_H

#include "TACSBeamElementBasis.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamUtilities.h"
#include "TACSGaussQuadrature.h"
#include "TACSElementAlgebra.h"
#include "TACSBeamConstitutive.h"
#include "TACSElement.h"
#include "TACSElementTypes.h"
#include "a2d.h"

template <int vars_per_node, class quadrature, class basis>
class TACSBeamCentrifugalForce : public TACSElement {
 public:
  TACSBeamCentrifugalForce( TACSBeamConstitutive *_con,
                            const TacsScalar _omegaVec[],
                            const TacsScalar _rotCenter[] ){
    con = _con;
    con->incref();
    memcpy(omegaVec, _omegaVec, 3*sizeof(TacsScalar));
    memcpy(rotCenter, _rotCenter, 3*sizeof(TacsScalar));
  }

  ~TACSBeamCentrifugalForce(){
    if (con){
      con->decref();
    }
  }

  const char* getObjectName(){
    return "TACSBeamCentrifugalForce";
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

  void addResidual( int elemIndex,
                    double time,
                    const TacsScalar *Xpts,
                    const TacsScalar *vars,
                    const TacsScalar *dvars,
                    const TacsScalar *ddvars,
                    TacsScalar *res ){

    // Compute the number of quadrature points
    const int nquad = quadrature::getNumQuadraturePoints();

    // Loop over each quadrature point and add the residual contribution
    for ( int quad_index = 0; quad_index < nquad; quad_index++ ){
      // Get the quadrature weight
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      // Tangent to the beam
      A2D::Vec3 X0, X0xi;

      // Compute X, X,xi and the interpolated normal
      basis::template interpFields<3, 3>(pt, Xpts, X0.x);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);

      // Compute the determinant of the transform
      A2D::Scalar detXd;
      A2D::Vec3Norm(X0xi, detXd);

      TacsScalar mass = con->evalDensity(elemIndex, pt, X0.x);

      TacsScalar r[3], wxr[3], ac[3];

      // Create vector pointing from rotation center to element gpt
      r[0] = X0.x[0] - rotCenter[0];
      r[1] = X0.x[1] - rotCenter[1];
      r[2] = X0.x[2] - rotCenter[2];

      // Compute omega x r
      crossProduct(omegaVec, r, wxr);

      // Compute centrifugal acceleration
      crossProduct(omegaVec, wxr, ac);

      // Compute the traction
      TacsScalar tr[3];
      tr[0] = detXd.value * weight * mass * ac[0];
      tr[1] = detXd.value * weight * mass * ac[1];
      tr[2] = detXd.value * weight * mass * ac[2];

      basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, tr, res);
    }
  }

 private:
  TacsScalar omegaVec[3], rotCenter[3];
  TACSBeamConstitutive* con;
};

#endif // TACS_BEAM_CENTRIFUGAL_FORCE_H