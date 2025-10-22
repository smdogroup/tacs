
#ifndef TACS_SHELL_TRACTION_H
#define TACS_SHELL_TRACTION_H

#include "TACSElementAlgebra.h"
#include "TACSShellUtilities.h"

template <int vars_per_node, class quadrature, class basis>
class TACSShellTraction : public TACSElement {
 public:
  TACSShellTraction(const TacsScalar _t[], int useConstTrac = 1) {
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

  const char *getObjectName() { return "TACSShellTraction"; }

  int getVarsPerNode() { return vars_per_node; }
  int getNumNodes() { return basis::NUM_NODES; }

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

    // Compute the node normal directions
    TacsScalar fn[3 * basis::NUM_NODES];
    TacsShellComputeNodeNormals<basis>(Xpts, fn);

    // Loop over each quadrature point and add the residual contribution
    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      // Get the quadrature weight
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      TacsScalar Xxi[6], n[3];
      basis::template interpFields<3, 3>(pt, fn, n);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

      // Assemble the terms Xd = [Xxi; n] and Xdz
      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n, Xd);

      // Compute the inverse of the 3x3 Jacobian transformation
      TacsScalar detXd = det3x3(Xd);
      detXd *= weight;

      // Interpolate the traction
      TacsScalar tr[3];
      basis::template interpFields<3, 3>(pt, t, tr);

      // Scale the traction
      tr[0] *= -detXd;
      tr[1] *= -detXd;
      tr[2] *= -detXd;

      basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, tr, res);
    }
  }

 private:
  TacsScalar t[3 * basis::NUM_NODES];
};

#endif  // TACS_SHELL_TRACTION_H