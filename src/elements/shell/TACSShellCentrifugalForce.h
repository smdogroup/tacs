
#ifndef TACS_SHELL_CENTRIFUGAL_FORCE_H
#define TACS_SHELL_CENTRIFUGAL_FORCE_H

#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"
#include "TACSShellUtilities.h"

template <int vars_per_node, class quadrature, class basis>
class TACSShellCentrifugalForce : public TACSElement {
 public:
  TACSShellCentrifugalForce(TACSShellConstitutive *_con,
                            const TacsScalar _omegaVec[],
                            const TacsScalar _rotCenter[],
                            bool _first_order = false) {
    con = _con;
    con->incref();
    memcpy(omegaVec, _omegaVec, 3 * sizeof(TacsScalar));
    memcpy(rotCenter, _rotCenter, 3 * sizeof(TacsScalar));
    first_order = _first_order;
  }

  ~TACSShellCentrifugalForce() {
    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSShellCentrifugalForce"; }

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

      TacsScalar Xxi[6], n[3], X[3], U[3];
      memset(U, 0, 3 * sizeof(TacsScalar));
      basis::template interpFields<3, 3>(pt, Xpts, X);
      if (first_order) {
        basis::template interpFields<vars_per_node, 3>(pt, vars, U);
      }
      basis::template interpFields<3, 3>(pt, fn, n);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

      // Assemble the terms Xd = [Xxi; n] and Xdz
      TacsScalar Xd[9];
      TacsShellAssembleFrame(Xxi, n, Xd);

      // Compute the inverse of the 3x3 Jacobian transformation
      TacsScalar detXd = det3x3(Xd);
      detXd *= weight;

      TacsScalar moments[3];
      con->evalMassMoments(elemIndex, pt, X, moments);
      TacsScalar mass = moments[0];

      TacsScalar r[3], wxr[3], ac[3];

      // Create vector pointing from rotation center to element gpt
      r[0] = X[0] - rotCenter[0] + U[0];
      r[1] = X[1] - rotCenter[1] + U[1];
      r[2] = X[2] - rotCenter[2] + U[2];

      // Account for shell mass offset
      if (!first_order) {
        for (int i = 0; i < 3; i++) {
          r[i] -= moments[1] / mass * n[i];
        }
      }

      // Compute omega x r
      crossProduct(omegaVec, r, wxr);

      // Compute centrifugal acceleration
      crossProduct(omegaVec, wxr, ac);

      // Compute the traction
      TacsScalar tr[vars_per_node] = {0.0};
      tr[0] = detXd * mass * ac[0];
      tr[1] = detXd * mass * ac[1];
      tr[2] = detXd * mass * ac[2];
      // Add moment terms if theres a shell offset
      if (!first_order) {
        crossProductAdd(-detXd * moments[1], n, ac, &tr[3]);
      }

      basis::template addInterpFieldsTranspose<vars_per_node, vars_per_node>(
          pt, tr, res);
    }
  }

  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]) {
    if (res || first_order) {
      // Compute the number of quadrature points
      const int nquad = quadrature::getNumQuadraturePoints();

      // Compute the node normal directions
      TacsScalar fn[3 * basis::NUM_NODES];
      TacsShellComputeNodeNormals<basis>(Xpts, fn);

      // Loop over each quadrature point and add the residual/jacobian
      // contributions
      for (int quad_index = 0; quad_index < nquad; quad_index++) {
        // Get the quadrature weight
        double pt[3];
        double weight = quadrature::getQuadraturePoint(quad_index, pt);

        TacsScalar Xxi[6], n[3], X[3], U[3];
        memset(U, 0, 3 * sizeof(TacsScalar));
        basis::template interpFields<3, 3>(pt, Xpts, X);
        if (first_order) {
          basis::template interpFields<vars_per_node, 3>(pt, vars, U);
        }
        basis::template interpFields<3, 3>(pt, fn, n);
        basis::template interpFieldsGrad<3, 3>(pt, Xpts, Xxi);

        // Assemble the terms Xd = [Xxi; n] and Xdz
        TacsScalar Xd[9];
        TacsShellAssembleFrame(Xxi, n, Xd);

        // Compute the determinant of the 3x3 Jacobian transformation
        TacsScalar scale = det3x3(Xd);
        scale *= weight;

        TacsScalar mass = con->evalDensity(elemIndex, pt, X);
        scale *= mass;

        // Add the contribution to the residual
        if (res) {
          TacsScalar r[3], wxr[3], ac[3];

          // Create vector pointing from rotation center to element gpt
          r[0] = X[0] - rotCenter[0] + U[0];
          r[1] = X[1] - rotCenter[1] + U[1];
          r[2] = X[2] - rotCenter[2] + U[2];

          // Compute omega x r
          crossProduct(omegaVec, r, wxr);

          // Compute centrifugal acceleration
          crossProduct(omegaVec, wxr, ac);

          // Compute the traction
          TacsScalar tr[3];
          tr[0] = scale * ac[0];
          tr[1] = scale * ac[1];
          tr[2] = scale * ac[2];
          basis::template addInterpFieldsTranspose<vars_per_node, 3>(pt, tr,
                                                                     res);
        }

        if (first_order) {
          // Compute the jacobian contribution
          // The Jacobian of the centrifugal force w.r.t the
          // location/displacement is: dtrdU = detXd * mass *
          // [[-w_y^2 - w_z^2,  w_x*w_y, w_x*w_z             ],
          //  [ w_x*w_y,       -w_x^2 - w_z^2,  w_y*w_z      ],
          //  [ w_x*w_z,        w_y*w_z,       -w_x^2 - w_y^2]]

          TacsScalar dtrdU[9];
          scale *= alpha;
          TacsScalar ww0 = omegaVec[0] * omegaVec[0];
          TacsScalar ww1 = omegaVec[1] * omegaVec[1];
          TacsScalar ww2 = omegaVec[2] * omegaVec[2];

          dtrdU[0] = (-ww1 - ww2) * scale;
          dtrdU[1] = (omegaVec[0] * omegaVec[1]) * scale;
          dtrdU[2] = (omegaVec[0] * omegaVec[2]) * scale;

          dtrdU[3] = dtrdU[1];
          dtrdU[4] = (-ww0 - ww2) * scale;
          dtrdU[5] = (omegaVec[1] * omegaVec[2]) * scale;

          dtrdU[6] = dtrdU[2];
          dtrdU[7] = dtrdU[5];
          dtrdU[8] = (-ww0 - ww1) * scale;

          // Add the contribution to the Jacobian, N^T * dtrdU * N
          basis::template addInterpFieldsOuterProduct<vars_per_node,
                                                      vars_per_node, 3, 3>(
              pt, dtrdU, mat);
        }
      }
    }
  }

 private:
  TacsScalar omegaVec[3], rotCenter[3];
  TACSShellConstitutive *con;
  bool first_order;
};

#endif  // TACS_SHELL_CENTRIFUGAL_FORCE_H
