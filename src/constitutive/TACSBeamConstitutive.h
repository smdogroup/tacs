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

#ifndef TACS_BEAM_CONSTITUTIVE_H
#define TACS_BEAM_CONSTITUTIVE_H

/*
  Base class for the Timoshenko beam constitutive object
*/

#include "TACSConstitutive.h"

class TACSBeamConstitutive : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 6;
  static const int NUM_TANGENT_STIFFNESS_ENTRIES = 21;

  TACSBeamConstitutive(){}

  /**
    Get the cross-sectional mass per unit area and the second moments
    of mass for the cross section

    moments = [rho * A, Iz1z1, Iz2z2, Iz1z2 ]

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The moments of the mass
  */
  virtual void evalMassMoments( int elemIndex,
                                const double pt[],
                                const TacsScalar X[],
                                TacsScalar moments[] ) = 0;

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param scale Scale factor for the moments
    @param dvLen the length of the sensitivity array
    @param dfdx The sensitivity array
  */
  virtual void addMassMomentsDVSens( int elemIndex,
                                     const double pt[],
                                     const TacsScalar X[],
                                     const TacsScalar scale[],
                                     int dvLen, TacsScalar dfdx[] ){}

  // Get the number of stress components
  int getNumStresses(){
    return NUM_STRESSES;
  }

  // Return the class name
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
    s[0] = Ct[0]*e[0] + Ct[1]*e[1]  + Ct[2]*e[2]  + Ct[3]*e[3]  + Ct[4]*e[4]  + Ct[5]*e[5];
    s[1] = Ct[1]*e[0] + Ct[6]*e[1]  + Ct[7]*e[2]  + Ct[8]*e[3]  + Ct[9]*e[4]  + Ct[10]*e[5];
    s[2] = Ct[2]*e[0] + Ct[7]*e[1]  + Ct[11]*e[2] + Ct[12]*e[3] + Ct[13]*e[4] + Ct[14]*e[5];
    s[3] = Ct[3]*e[0] + Ct[8]*e[1]  + Ct[12]*e[2] + Ct[15]*e[3] + Ct[16]*e[4] + Ct[17]*e[5];
    s[4] = Ct[4]*e[0] + Ct[9]*e[1]  + Ct[13]*e[2] + Ct[16]*e[3] + Ct[18]*e[4] + Ct[19]*e[5];
    s[5] = Ct[5]*e[0] + Ct[10]*e[1] + Ct[14]*e[2] + Ct[17]*e[3] + Ct[19]*e[4] + Ct[20]*e[5];
  }

 private:
  // Set the constitutive name
  static const char *constName;
};

#endif // TACS_BEAM_CONSTITUTIVE_H
