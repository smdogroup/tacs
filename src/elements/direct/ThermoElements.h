#include "CoupledThermoPlaneStressStiffness.h"
#include "CoupledThermoSolidStiffness.h"
#include "TACSElement.h"

class ThermoQuad : public TACSElement {
 public:
  ThermoQuad(int component) : TACSElement(component){};
  virtual void getShapeFunctions(const double pt[], double N[]) = 0;
  virtual void getShapeFunctions(const double pt[], double N[], double Na[],
                                 double Nb[]) = 0;
  virtual void getBT(TacsScalar strain[], const double pt[],
                     const TacsScalar Xpts[], const TacsScalar vars[]) = 0;
  virtual void addBTSVSens(TacsScalar strainSVSens[], const double pt[],
                           const TacsScalar scale,
                           const TacsScalar strainSens[],
                           const TacsScalar Xpts[],
                           const TacsScalar vars[]) = 0;
  virtual void addEffStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                                  const TacsScalar scale,
                                  const TacsScalar strainSens[],
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[], int vars_j = 0) = 0;
  virtual void getTemperature(TacsScalar T[], const double N[],
                              const TacsScalar vars[]) = 0;
};

class ThermoSolid : public TACSElement {
 public:
  ThermoSolid(int component) : TACSElement(component){};
  virtual void getShapeFunctions(const double pt[], double N[]) = 0;
  virtual void getShapeFunctions(const double pt[], double N[], double Na[],
                                 double Nb[], double Nc[]) = 0;
  virtual void getBT(TacsScalar strain[], const double pt[],
                     const TacsScalar Xpts[], const TacsScalar vars[]) = 0;
  virtual void addBTSVSens(TacsScalar strainSVSens[], const double pt[],
                           const TacsScalar scale,
                           const TacsScalar strainSens[],
                           const TacsScalar Xpts[],
                           const TacsScalar vars[]) = 0;
  virtual void addEffStrainSVSens(TacsScalar strainSVSens[], const double pt[],
                                  const TacsScalar scale,
                                  const TacsScalar strainSens[],
                                  const TacsScalar Xpts[],
                                  const TacsScalar vars[], int vars_j = 0) = 0;
  virtual void getTemperature(TacsScalar T[], const double N[],
                              const TacsScalar vars[]) = 0;
};
