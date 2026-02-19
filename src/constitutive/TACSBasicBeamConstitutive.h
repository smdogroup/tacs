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

#ifndef TACS_BASIC_BEAM_CONSTITUTIVE_H
#define TACS_BASIC_BEAM_CONSTITUTIVE_H

/*
  Base class for the Timoshenko beam constitutive object
*/

#include "TACSBeamConstitutive.h"
#include "TACSMaterialProperties.h"

class TACSBasicBeamConstitutive : public TACSBeamConstitutive {
 public:
  TACSBasicBeamConstitutive(TacsScalar EA, TacsScalar EI22, TacsScalar EI33,
                            TacsScalar EI23, TacsScalar GJ, TacsScalar kG22,
                            TacsScalar kG33, TacsScalar kG23, TacsScalar m00,
                            TacsScalar m11, TacsScalar m22, TacsScalar m33,
                            TacsScalar xm2, TacsScalar xm3, TacsScalar xc2,
                            TacsScalar xc3, TacsScalar xk2, TacsScalar xk3,
                            TacsScalar muS);
  TACSBasicBeamConstitutive(TacsScalar rhoA, TacsScalar rhoIy, TacsScalar rhoIz,
                            TacsScalar rhoIyz, TacsScalar EA, TacsScalar GJ,
                            TacsScalar EIy, TacsScalar EIz, TacsScalar kGAy,
                            TacsScalar kGAz);
  TACSBasicBeamConstitutive(TACSMaterialProperties *properties, TacsScalar A,
                            TacsScalar J, TacsScalar Iy, TacsScalar Iz,
                            TacsScalar Iyz, TacsScalar ky = 5.0 / 6.0,
                            TacsScalar kz = 5.0 / 6.0);

  ~TACSBasicBeamConstitutive();

  /**
    Get the cross-sectional mass per unit area and the second moments
    of mass for the cross section

    moments = [rho * A, Iz1z1, Iz2z2, Iz1z2 ]

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The moments of the mass
  */
  void evalMassMoments(int elemIndex, const double pt[], const TacsScalar X[],
                       TacsScalar moments[]) {
    moments[0] = rho[0];
    moments[1] = rho[1];
    moments[2] = rho[2];
    moments[3] = rho[3];
    moments[4] = rho[4];
    moments[5] = rho[5];
  }

  /**
    Add the derivative of the pointwise mass times the given scalar

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @param scale Scale factor for the moments
    @param dvLen the length of the sensitivity array
    @param dfdx The sensitivity array
  */
  virtual void addMassMomentsDVSens(int elemIndex, const double pt[],
                                    const TacsScalar X[],
                                    const TacsScalar scale[], int dvLen,
                                    TacsScalar dfdx[]) {}

  /**
    Evaluate the mass per unit length of the beam
  */
  TacsScalar evalDensity(int elemIndex, const double pt[],
                         const TacsScalar X[]) {
    return rho[0];
  }

  /**
     Compute the stress at the material point, given the strain value
  */
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar e[], TacsScalar s[]) {
    computeStress(C, e, s);
  }

  /**
     Compute (or in this case copy) the values of the tangent stiffness matrix
  */
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C0[]) {
    memcpy(C0, C, NUM_TANGENT_STIFFNESS_ENTRIES * sizeof(TacsScalar));
  }

  /**
     Evaluate the specific heat
  */
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]) {
    return 0.0;
  }

  // The name of the constitutive object
  const char *getObjectName();

 private:
  // The constitutive matrix
  TacsScalar C[36];

  // The moments of the density
  TacsScalar rho[6];

  // Material properties class
  TACSMaterialProperties *props;

  // The object name
  static const char *constName;
};

#endif  // TACS_BASIC_BEAM_CONSTITUTIVE_H
