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

  //Compute the design dependent thermal loading stress
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
