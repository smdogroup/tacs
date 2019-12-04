/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef TACS_TIMOSHENKO_STIFFNESS_H
#define TACS_TIMOSHENKO_STIFFNESS_H

/*
  Base class for the Timoshenko beam constitutive object
*/

#include "TACSConstitutive.h"

class TACSTimoshenkoConstitutive : public TACSConstitutive {
 public:
  TACSTimoshenkoConstitutive( const TacsScalar _axis[],
                              TacsScalar EA,
                              TacsScalar EI22, TacsScalar EI33, TacsScalar EI23,
                              TacsScalar GJ,
                              TacsScalar kG22, TacsScalar kG33, TacsScalar kG23,
                              TacsScalar m00,
                              TacsScalar m11, TacsScalar m22, TacsScalar m33,
                              TacsScalar xm2, TacsScalar xm3,
                              TacsScalar xc2, TacsScalar xc3,
                              TacsScalar xk2, TacsScalar xk3,
                              TacsScalar muS );
  TACSTimoshenkoConstitutive( TacsScalar rhoA, TacsScalar rhoIy,
                              TacsScalar rhoIz, TacsScalar rhoIyz,
                              TacsScalar EA, TacsScalar GJ,
                              TacsScalar EIy, TacsScalar EIz,
                              TacsScalar kGAy, TacsScalar kGAz,
                              const TacsScalar axis[] );
  TACSTimoshenkoConstitutive( const TacsScalar rho[],
                              const TacsScalar C[],
                              const TacsScalar axis[] );
  virtual ~TACSTimoshenkoConstitutive();

  /**
     Set the mass moments, stiffness matrix and reference axis into
     the constitutive object

     Note that a null argument is ignored on input.

     @param _rho The mass moments about the axis
     @param _C The stiffness matrix
     @param _axis The reference axis used to compute the stiffness
  */
  void setProperties( const TacsScalar _rho[], const TacsScalar _C[],
                      const TacsScalar _axis[] );


  /**
     Get the mass moments, stiffness matrix and reference axis into
     the constitutive object

     Note that a null argument is ignored.

     @param _rho The mass moments about the axis
     @param _C The stiffness matrix
     @param _axis The reference axis used to compute the stiffness
  */
  void getProperties( TacsScalar _rho[], TacsScalar _C[],
                      TacsScalar _axis[] );

  /**
    Get the reference axis for the beam.

    This reference axis defines the direction of the stiffest bending
    component of the beam.

    @return Constant pointer to the three axis components
  */
  const TacsScalar *getRefAxis(){ return axis; }

  /**
    Get the mass moments of the beam

    The mass moments consist of the mass per unit area and the
  */
  virtual void getMassMoments( int elemIndex, const double pt[], TacsScalar moments[] ){
    moments[0] = rho[0];
    moments[1] = rho[1];
    moments[2] = rho[2];
    moments[3] = rho[3];
  }

  virtual void addMassMomentsDVSens( int elemIndex, const double pt[],
                                     TacsScalar scale, const TacsScalar psi[],
                                     int dvLen, TacsScalar dfdx[] ){}
  
  // Get the number of stress components
  int getNumStresses();

  // Evaluate material properties
  TacsScalar evalDensity( int elemIndex, const double pt[],
                          const TacsScalar X[] );
  TacsScalar evalSpecificHeat( int elemIndex, const double pt[],
                               const TacsScalar X[] );

  // Evaluate the stress and the tangent stiffness matrix
  void evalStress( int elemIndex, const double pt[],
                   const TacsScalar X[], const TacsScalar e[],
                   TacsScalar s[] );
  void evalTangentStiffness( int elemIndex, const double pt[],
                             const TacsScalar X[], TacsScalar C[] );

  // Class name
  const char *getObjectName();

  /**
    Compute a matrix-vector product using the tangent stiffness matrix

    @param Ct The tangent stiffness matrix
    @param e The component of the strain
    @param s The components of the stress
  */
  static inline void computeStress( const TacsScalar Ct[],
                                    const TacsScalar e[],
                                    TacsScalar s[] ){
    s[0] = Ct[0]*e[0]  + Ct[1]*e[1]  + Ct[2]*e[2]  + Ct[3]*e[3]  + Ct[4]*e[4]  + Ct[5]*e[5];
    s[1] = Ct[6]*e[0]  + Ct[7]*e[1]  + Ct[8]*e[2]  + Ct[9]*e[3]  + Ct[10]*e[4] + Ct[11]*e[5];
    s[2] = Ct[12]*e[0] + Ct[13]*e[1] + Ct[14]*e[2] + Ct[15]*e[3] + Ct[16]*e[4] + Ct[17]*e[5];
    s[3] = Ct[18]*e[0] + Ct[19]*e[1] + Ct[20]*e[2] + Ct[21]*e[3] + Ct[22]*e[4] + Ct[23]*e[5];
    s[4] = Ct[24]*e[0] + Ct[25]*e[1] + Ct[26]*e[2] + Ct[27]*e[3] + Ct[28]*e[4] + Ct[29]*e[5];
    s[5] = Ct[30]*e[0] + Ct[31]*e[1] + Ct[32]*e[2] + Ct[33]*e[3] + Ct[34]*e[4] + Ct[35]*e[5];
  }

 protected:
  // The constitutive matrix
  TacsScalar C[36];

  // The moments of the density
  TacsScalar rho[4];

  // The reference axis - parallel with the local z-direction of the
  // beam. This direction cannot be parallel with the beam.
  TacsScalar axis[3];

  // Set the constitutive name
  static const char *constName;
};

#endif // TACS_TIMOSHENKO_STIFFNESS_H
