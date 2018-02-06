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

class TimoshenkoStiffness : public TACSConstitutive {
 public:
  TimoshenkoStiffness( const TacsScalar _axis[],
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
  TimoshenkoStiffness( TacsScalar rhoA, TacsScalar rhoIy,
                       TacsScalar rhoIz, TacsScalar rhoIyz,
                       TacsScalar EA, TacsScalar GJ,
                       TacsScalar EIy, TacsScalar EIz,
                       TacsScalar kGAy, TacsScalar kGAz,
                       const TacsScalar axis[] );
  TimoshenkoStiffness( const TacsScalar rho[], 
                       const TacsScalar C[],
                       const TacsScalar axis[] );
  virtual ~TimoshenkoStiffness();

  // Retrieve the reference axis
  // ---------------------------
  const TacsScalar *getRefAxis(){ return axis; }

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress( const double pt[], 
                        const TacsScalar strain[],
			TacsScalar stress[] );
  
  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 4; }
  void getPointwiseMass( const double pt[], TacsScalar mass[] ){
    mass[0] = rho[0];
    mass[1] = rho[1];
    mass[2] = rho[2];
    mass[3] = rho[3];
  }

  // Extra info about the constitutive class
  // ---------------------------------------
  const char *constitutiveName();
  
 protected:
  // Compute the stress using an inline function
  inline void calcStress( const TacsScalar e[], TacsScalar s[] );

  void setData( const TacsScalar _rho[],
                const TacsScalar _C[],
                const TacsScalar _axis[] );

 private:
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

/*
  Given the strain, compute the stress
*/
inline void TimoshenkoStiffness::calcStress( const TacsScalar e[],
                                             TacsScalar s[] ){
  s[0] = C[0]*e[0]  + C[1]*e[1]  + C[2]*e[2]  + C[3]*e[3]  + C[4]*e[4]  + C[5]*e[5];
  s[1] = C[6]*e[0]  + C[7]*e[1]  + C[8]*e[2]  + C[9]*e[3]  + C[10]*e[4] + C[11]*e[5];
  s[2] = C[12]*e[0] + C[13]*e[1] + C[14]*e[2] + C[15]*e[3] + C[16]*e[4] + C[17]*e[5];
  s[3] = C[18]*e[0] + C[19]*e[1] + C[20]*e[2] + C[21]*e[3] + C[22]*e[4] + C[23]*e[5];
  s[4] = C[24]*e[0] + C[25]*e[1] + C[26]*e[2] + C[27]*e[3] + C[28]*e[4] + C[29]*e[5];
  s[5] = C[30]*e[0] + C[31]*e[1] + C[32]*e[2] + C[33]*e[3] + C[34]*e[4] + C[35]*e[5];
}

#endif // TACS_TIMOSHENKO_STIFFNESS_H
