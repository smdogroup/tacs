
#ifndef TACS_SHELL_INERTIAL_FORCE_H
#define TACS_SHELL_INERTIAL_FORCE_H

#include "TACSElementAlgebra.h"
#include "TACSShellConstitutive.h"
#include "TACSShellUtilities.h"

template <int vars_per_node, class quadrature, class basis>
class TACSShellInertialForce : public TACSElement {
 public:
  TACSShellInertialForce(TACSShellConstitutive *_con,
                         const TacsScalar _inertiaVec[],
                         const int *_inertiaVecDVNums = NULL) {
    con = _con;
    con->incref();
    memcpy(inertiaVec, _inertiaVec, 3 * sizeof(TacsScalar));
    if (_inertiaVecDVNums) {
      memcpy(inertiaVecDVNums, _inertiaVecDVNums, 3 * sizeof(int));
    } else {
      inertiaVecDVNums[0] = inertiaVecDVNums[1] = inertiaVecDVNums[2] = -1;
    }
  }

  ~TACSShellInertialForce() {
    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSShellInertialForce"; }

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
    int num = con->getDesignVarNums(elemIndex, dvLen, dvNums);
    for (int i = 0; i < 3; i++) {
      if (inertiaVecDVNums[i] >= 0) {
        if (dvNums && num < dvLen) {
          dvNums[num] = inertiaVecDVNums[i];
        }
        num++;
      }
    }
    return num;
  }

  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    int num = con->setDesignVars(elemIndex, dvLen, dvs);
    for (int i = 0; i < 3; i++) {
      if (inertiaVecDVNums[i] >= 0) {
        if (num < dvLen) {
          inertiaVec[i] = dvs[num];
        }
        num++;
      }
    }
    return num;
  }

  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    int num = con->getDesignVars(elemIndex, dvLen, dvs);
    for (int i = 0; i < 3; i++) {
      if (inertiaVecDVNums[i] >= 0) {
        if (dvs && num < dvLen) {
          dvs[num] = inertiaVec[i];
        }
        num++;
      }
    }
    return num;
  }

  int getDesignVarRange(int elemIndex, int dvLen, TacsScalar lb[],
                        TacsScalar ub[]) {
    int num = con->getDesignVarRange(elemIndex, dvLen, lb, ub);
    for (int i = 0; i < 3; i++) {
      if (inertiaVecDVNums[i] >= 0) {
        if (num < dvLen) {
          lb[num] = -1e20;
          ub[num] = 1e20;
        }
        num++;
      }
    }
    return num;
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

      TacsScalar Xxi[6], n[3], X[3];
      basis::template interpFields<3, 3>(pt, Xpts, X);
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

      // Compute the traction
      TacsScalar tr[vars_per_node] = {0.0};
      tr[0] = -detXd * mass * inertiaVec[0];
      tr[1] = -detXd * mass * inertiaVec[1];
      tr[2] = -detXd * mass * inertiaVec[2];
      // Add moment terms if theres a shell offset
      crossProductAdd(detXd * moments[1], n, inertiaVec, &tr[3]);

      basis::template addInterpFieldsTranspose<vars_per_node, vars_per_node>(
          pt, tr, res);
    }
  }

 private:
  TacsScalar inertiaVec[3];
  int inertiaVecDVNums[3];
  TACSShellConstitutive *con;
};

#endif  // TACS_SHELL_INERTIAL_FORCE_H