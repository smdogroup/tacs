#ifndef TACS_BEAM_TRACTION_H
#define TACS_BEAM_TRACTION_H

#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamUtilities.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "a2d.h"

template <int vars_per_node, class quadrature, class basis>
class TACSBeamTraction : public TACSElement {
 public:
  TACSBeamTraction(const TacsScalar _t[], int useConstTrac = 1) {
    // Traction is constant across element
    if (useConstTrac) {
      for (int i = 0; i < basis::NUM_NODES; i++) {
        for (int j = 0; j < 3; j++) {
          t[3 * i + j] = _t[j];
        }
      }
    }
    // Traction varies across element
    else {
      for (int i = 0; i < 3 * basis::NUM_NODES; i++) {
        t[i] = _t[i];
      }
    }
  }

  const char *getObjectName() { return "TACSBeamTraction"; }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return basis::NUM_NODES; }

  ElementLayout getLayoutType() { return basis::getLayoutType(); }

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

  void addResidual(int elemIndex, double time, const TacsScalar *Xpts,
                   const TacsScalar *vars, const TacsScalar *dvars,
                   const TacsScalar *ddvars, TacsScalar *res) {
    // Compute the number of quadrature points
    const int nquad = quadrature::getNumQuadraturePoints();

    // Loop over each quadrature point and add the residual contribution
    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      // Get the quadrature weight
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      // Tangent to the beam
      A2D::Vec3 X0xi;

      // Compute X, X,xi and the interpolated normal
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);

      // Compute the determinant of the transform
      A2D::Scalar detXd;
      A2D::Vec3Norm(X0xi, detXd);

      // Compute the traction
      TacsScalar tr[3];
      basis::template interpFields<3, 3>(pt, t, tr);
      tr[0] *= -detXd.value * weight;
      tr[1] *= -detXd.value * weight;
      tr[2] *= -detXd.value * weight;

      basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, tr, res);
    }
  }

 private:
  TacsScalar t[3 * basis::NUM_NODES];
};

#endif  // TACS_BEAM_TRACTION_H