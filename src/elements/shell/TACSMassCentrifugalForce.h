#ifndef TACS_MASS_CENTRIFUGAL_FORCE_H
#define TACS_MASS_CENTRIFUGAL_FORCE_H

#include "TACSElement.h"
#include "TACSGeneralMassConstitutive.h"

class TACSMassCentrifugalForce : public TACSElement {
 public:
  TACSMassCentrifugalForce(TACSGeneralMassConstitutive *_con,
                           const TacsScalar *_omegaVec,
                           const TacsScalar *_rotCenter,
                           const int *_omegaDVNums = NULL,
                           const int *_rotCenterDVNums = NULL);
  ~TACSMassCentrifugalForce();

  // Get the element properties and names
  // ------------------------------------
  const char *getObjectName();
  int getVarsPerNode();
  int getNumNodes();
  int getNumQuadraturePoints() { return 1; }
  double getQuadratureWeight(int n) { return 1.0; }
  double getQuadraturePoint(int n, double pt[]) { return 1.0; }
  int getNumElementFaces() { return 0; }
  int getNumFaceQuadraturePoints(int face) { return 0; }
  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return 0.0;
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
          lb[num] = -1e20;
          ub[num] = 1e20;
        }
        num++;
      }
    }
    for (int i = 0; i < 3; i++) {
      if (rotCenterDVNums[i] >= 0) {
        if (num < dvLen) {
          lb[num] = -1e20;
          ub[num] = 1e20;
        }
        num++;
      }
    }
    return num;
  }

  // Functions for analysis
  // ----------------------

  void addResidual(int elemIndex, double time, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[]);

  void addJacobian(int elemIndex, double time, TacsScalar alpha,
                   TacsScalar beta, TacsScalar gamma, const TacsScalar Xpts[],
                   const TacsScalar vars[], const TacsScalar dvars[],
                   const TacsScalar ddvars[], TacsScalar res[],
                   TacsScalar mat[]);

  void addAdjResProduct(int elemIndex, double time, TacsScalar scale,
                        const TacsScalar psi[], const TacsScalar Xpts[],
                        const TacsScalar vars[], const TacsScalar dvars[],
                        const TacsScalar ddvars[], int dvLen,
                        TacsScalar dfdx[]);

 private:
  static const int NUM_DISPS = 6;
  static const int NUM_NODES = 1;
  static const int NUM_VARIABLES = NUM_NODES * NUM_DISPS;

  TACSGeneralMassConstitutive *con;
  TacsScalar omegaVec[3], rotCenter[3];
  int omegaDVNums[3], rotCenterDVNums[3];
};

#endif  // TACS_MASS_CENTRIFUGAL_FORCE_H
