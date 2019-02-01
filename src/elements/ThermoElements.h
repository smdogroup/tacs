#include "TACS3DCoupledThermoElement.h"
#include "TACS2DCoupledThermoElement.h"
#include "TACSElement.h"
#include "CoupledThermoPlaneStressStiffness.h"
#include "CoupledThermoSolidStiffness.h"

class ThermoQuad : public TACSElement {
 public:
  ThermoQuad( CoupledThermoPlaneStressStiffness * _stiff, 
              ElementBehaviorType type=LINEAR, 
              int _componentNum=0 ){}
  virtual void getShapeFunctions( const double pt[], double N[],
                                  double Na[], double Nb[] );
  virtual void getBT( TacsScalar strain[], const double pt[], 
                      const TacsScalar Xpts[], const TacsScalar vars[] );
  virtual void addBTSVSens( TacsScalar strainSVSens[],
                            const double pt[],
                            const TacsScalar scale,
                            const TacsScalar strainSens[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[] );
};

class ThermoSolid : public TACSElement {
 public:
  ThermoSolid( CoupledThermoSolidStiffness * _stiff, 
               ElementBehaviorType type=LINEAR, 
               int _componentNum=0 ){}
  virtual void getShapeFunctions( const double pt[], double N[],
                                  double Na[], double Nb[] );
  virtual void getBT( TacsScalar strain[], const double pt[], 
                      const TacsScalar Xpts[], const TacsScalar vars[] );
  virtual void addBTSVSens( TacsScalar strainSVSens[],
                            const double pt[],
                            const TacsScalar scale,
                            const TacsScalar strainSens[],
                            const TacsScalar Xpts[],
                            const TacsScalar vars[] );
};
