#ifndef TACS_MASS_CENTRIFUGAL_FORCE_H
#define TACS_MASS_CENTRIFUGAL_FORCE_H

#include "TACSElement.h"
#include "TACSGeneralMassConstitutive.h"

class TACSMassCentrifugalForce : public TACSElement {
 public:
  TACSMassCentrifugalForce(TACSGeneralMassConstitutive* _con,
                           const TacsScalar* _omegaVec,
                           const TacsScalar* _rotCenter);
  ~TACSMassCentrifugalForce();

  // Get the element properties and names
  // ------------------------------------
  const char* getObjectName();
  int getVarsPerNode();
  int getNumNodes();
  int getDesignVarsPerNode() { return 0; }
  int getNumQuadraturePoints() { return 1; }
  double getQuadratureWeight(int n) { return 1.0; }
  double getQuadraturePoint(int n, double pt[]) { return 1.0; }
  int getNumElementFaces() { return 0; }
  int getNumFaceQuadraturePoints(int face) { return 0; }
  double getFaceQuadraturePoint(int face, int n, double pt[],
                                double tangent[]) {
    return 0.0;
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

 private:
  static const int NUM_DISPS = 6;
  static const int NUM_NODES = 1;
  static const int NUM_VARIABLES = NUM_NODES * NUM_DISPS;

  TACSGeneralMassConstitutive* con;
  TacsScalar omegaVec[3], rotCenter[3];
};

#endif  // TACS_MASS_CENTRIFUGAL_FORCE_H
