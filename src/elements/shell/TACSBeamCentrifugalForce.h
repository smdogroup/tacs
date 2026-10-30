#ifndef TACS_BEAM_CENTRIFUGAL_FORCE_H
#define TACS_BEAM_CENTRIFUGAL_FORCE_H

#include <algorithm>

#include "TACSBeamConstitutive.h"
#include "TACSBeamElementBasis.h"
#include "TACSBeamElementQuadrature.h"
#include "TACSBeamElementTransform.h"
#include "TACSBeamUtilities.h"
#include "TACSElement.h"
#include "TACSElementAlgebra.h"
#include "TACSElementTypes.h"
#include "TACSGaussQuadrature.h"
#include "a2d.h"

template <int vars_per_node, class quadrature, class basis>
class TACSBeamCentrifugalForce : public TACSElement {
 public:
  TACSBeamCentrifugalForce(TACSBeamTransform *_transform,
                           TACSBeamConstitutive *_con,
                           const TacsScalar _omegaVec[],
                           const TacsScalar _rotCenter[],
                           const int *_omegaDVNums = NULL,
                           const int *_rotCenterDVNums = NULL) {
    transform = _transform;
    transform->incref();
    con = _con;
    con->incref();
    std::copy_n(_omegaVec, 3, omegaVec);
    std::copy_n(_rotCenter, 3, rotCenter);
    if (_omegaDVNums) {
      std::copy_n(_omegaDVNums, 3, omegaDVNums);
    } else {
      omegaDVNums[0] = omegaDVNums[1] = omegaDVNums[2] = -1;
    }
    if (_rotCenterDVNums) {
      std::copy_n(_rotCenterDVNums, 3, rotCenterDVNums);
    } else {
      rotCenterDVNums[0] = rotCenterDVNums[1] = rotCenterDVNums[2] = -1;
    }
  }

  ~TACSBeamCentrifugalForce() {
    if (transform) {
      transform->decref();
    }
    if (con) {
      con->decref();
    }
  }

  const char *getObjectName() { return "TACSBeamCentrifugalForce"; }

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
      if (omegaDVNums[i] >= 0) {
        if (dvNums && num < dvLen) {
          dvNums[num] = omegaDVNums[i];
        }
        num++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (dvNums && num < dvLen) {
          dvNums[num] = rotCenterDVNums[i];
        }
        num++;
      }
    }
    return num;
  }

  int setDesignVars(int elemIndex, int dvLen, const TacsScalar dvs[]) {
    int num = con->setDesignVars(elemIndex, dvLen, dvs);
    for (int i = 0; i < 3; i++) {
      if (omegaDVNums[i] >= 0) {
        if (num < dvLen) {
          omegaVec[i] = dvs[num];
        }
        num++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (num < dvLen) {
          rotCenter[i] = dvs[num];
        }
        num++;
      }
    }
    return num;
  }

  int getDesignVars(int elemIndex, int dvLen, TacsScalar dvs[]) {
    int num = con->getDesignVars(elemIndex, dvLen, dvs);
    for (int i = 0; i < 3; i++) {
      if (omegaDVNums[i] >= 0) {
        if (dvs && num < dvLen) {
          dvs[num] = omegaVec[i];
        }
        num++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (dvs && num < dvLen) {
          dvs[num] = rotCenter[i];
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
      if (omegaDVNums[i] >= 0) {
        if (num < dvLen) {
          lb[num] = -TACS_LARGE_DV_BOUND;
          ub[num] = TACS_LARGE_DV_BOUND;
        }
        num++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (num < dvLen) {
          lb[num] = -TACS_LARGE_DV_BOUND;
          ub[num] = TACS_LARGE_DV_BOUND;
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

    // Get the reference axis
    const A2D::Vec3 &axis = transform->getRefAxis();

    // Compute the normal directions
    TacsScalar fn1[3 * basis::NUM_NODES], fn2[3 * basis::NUM_NODES];
    TacsBeamComputeNodeNormals<basis>(Xpts, axis, fn1, fn2);

    // Loop over each quadrature point and add the residual contribution
    for (int quad_index = 0; quad_index < nquad; quad_index++) {
      // Get the quadrature weight
      double pt[3];
      double weight = quadrature::getQuadraturePoint(quad_index, pt);

      // Tangent to the beam
      A2D::Vec3 X0, X0xi;

      // Interpolated normal directions
      A2D::ADVec3 n1, n2;

      // Compute X, X,xi and the interpolated normal
      basis::template interpFields<3, 3>(pt, Xpts, X0.x);
      basis::template interpFieldsGrad<3, 3>(pt, Xpts, X0xi.x);
      basis::template interpFields<3, 3>(pt, fn1, n1.x);
      basis::template interpFields<3, 3>(pt, fn2, n2.x);

      // Compute the determinant of the transform
      A2D::Scalar detXd;
      A2D::Vec3Norm(X0xi, detXd);

      TacsScalar mass_moment[6];
      con->evalMassMoments(elemIndex, pt, X0.x, mass_moment);
      TacsScalar mass = mass_moment[0];

      TacsScalar r[3], wxr[3], ac[3];

      // Create vector pointing from rotation center to element gpt
      r[0] = X0.x[0] - rotCenter[0];
      r[1] = X0.x[1] - rotCenter[1];
      r[2] = X0.x[2] - rotCenter[2];

      // Account for beam offset
      for (int i = 0; i < 3; i++) {
        r[i] -= (mass_moment[1] * n1.x[i] + mass_moment[2] * n2.x[i]) / mass;
      }

      // Compute omega x r
      crossProduct(omegaVec, r, wxr);

      // Compute centrifugal acceleration
      crossProduct(omegaVec, wxr, ac);

      // Compute the traction
      TacsScalar tr[vars_per_node] = {0.0};
      tr[0] = detXd.value * weight * mass * ac[0];
      tr[1] = detXd.value * weight * mass * ac[1];
      tr[2] = detXd.value * weight * mass * ac[2];
      // Add moment terms if theres a beam offset
      crossProductAdd(-detXd.value * weight * mass_moment[1], n1.x, ac, &tr[3]);
      crossProductAdd(-detXd.value * weight * mass_moment[2], n2.x, ac, &tr[3]);

      basis::template addInterpFieldsTranspose<vars_per_node, vars_per_node>(
          pt, tr, res);
    }
  }

 private:
  TacsScalar omegaVec[3], rotCenter[3];
  int omegaDVNums[3], rotCenterDVNums[3];
  TACSBeamTransform *transform;
  TACSBeamConstitutive *con;
};

#endif  // TACS_BEAM_CENTRIFUGAL_FORCE_H