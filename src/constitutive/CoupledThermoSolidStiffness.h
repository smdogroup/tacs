#ifndef TACS_COUPLED_THERMO_SOLID_STIFFNESS_H
#define TACS_COUPLED_THERMO_SOLID_STIFFNESS_H

/*
  Copyright (c) 2017 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/
#include "SolidStiffness.h"
/*
  This is the plane stress constitutive objects with thermal loading. 
*/
class CoupledThermoSolidStiffness : public SolidStiffness {
 public:
  static const int NUM_STRESSES = 6;
  CoupledThermoSolidStiffness();
  CoupledThermoSolidStiffness( TacsScalar _rho, TacsScalar E,
                               TacsScalar nu, TacsScalar _alpha, 
                               TacsScalar Tref, TacsScalar kcond );

  ~CoupledThermoSolidStiffness();
  // Calculate the product B*u
  // ---------------------------------------------------------------
  virtual void calculateStress( const double pt[],
                                const TacsScalar strain[],
                                TacsScalar stress[] ) = 0;
  virtual void addStressDVSens( const double pt[],
                                const TacsScalar e[],
                                TacsScalar alpha,
                                const TacsScalar psi[],
                                TacsScalar fdvSens[],
                                int dvLen ) = 0;

  // Compute the design dependent thermal loading stress
  //---------------------------------------------------
  virtual TacsScalar getEffThermalAlpha(int var_j) = 0;
  virtual void calculateThermal( const double pt[],
                                 const TacsScalar strain[],
                                 TacsScalar stress[] ) = 0;
  virtual void addThermalDVSens( const double pt[],
                                 const TacsScalar e[],
                                 TacsScalar alpha,
                                 const TacsScalar psi[],
                                 TacsScalar fdvSens[],
                                 int dvLen ) = 0;
  virtual void calculateConduction( const double pt[],
                                    const TacsScalar strain[],
                                    TacsScalar stress[] ) = 0;
  virtual void addConductionDVSens(const double pt[],
                                   const TacsScalar e[],
                                   TacsScalar alpha,
                                   const TacsScalar psi[],
                                   TacsScalar fdvSens[],
                                   int dvLen) = 0;
  // Functions that evaluate the failure function
  virtual void failure( const double pt[],
                        const TacsScalar T[],
                        const TacsScalar strain[],
                        TacsScalar * fail ) = 0;
  virtual void addFailureDVSens( const double pt[],
                                 const TacsScalar T[],
                                 const TacsScalar strain[],
                                 TacsScalar alpha,
                                 TacsScalar dvSens[], int dvLen ) = 0;
  virtual void failureStrainSens(const double pt[],
                                 const TacsScalar T[],
                                 const TacsScalar strain[],
                                 TacsScalar sens[], 
                                 int vars_j=0 ) = 0;
  virtual int getVarsPerNode() = 0;
  // Functions that evaluate the heat flux at a location
  virtual void heatflux( const double pt[],
                         const TacsScalar normal[],
                         const TacsScalar strain[],
                         TacsScalar * qn ) = 0;
  virtual void addHeatFluxDVSens( const double pt[],
                                  const TacsScalar normal[],
                                  const TacsScalar strain[],
                                  TacsScalar alpha,
                                  TacsScalar dvSens[], int dvLen ) = 0;
  virtual void heatfluxStrainSens( const double pt[],
                                   const TacsScalar normal[],
                                   TacsScalar sens[] ) = 0;
  // Functions that evaluate the maximum temperature of an element
  virtual void maxtemp( const double pt[],
			const TacsScalar max_temp,
			TacsScalar *fail,
			int vars_j=0 ) = 0;
  virtual void addMaxTempDVSens ( const double pt[], 
				  const TacsScalar max_temp,
				  TacsScalar alpha,
                                  TacsScalar dvSens[], int dvLen,
				  int vars_j=0 ) = 0;
  virtual void maxtempStrainSens( const double pt[],
				  const TacsScalar max_temp,
				  TacsScalar sens[],
				  int vars_j=0 ) = 0;
  
  // Extra info about the constitutive class
  // ---------------------------------------
  const char *constitutiveName();
 protected:
  // The stiffness matrix
  TacsScalar Tmat[3]; 
  // The stiffness parameters
  TacsScalar C[6];
  TacsScalar G23, G13, G12;
  TacsScalar rho, alpha,xw, Tref;
 private:
  static const char *constName;
};
#endif // TACS_THERMO_PLANE_STRESS_STIFFNESS_H
