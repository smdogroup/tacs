/*
  The stiffness object for the variety of springs

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Copyright (C) 2020 Aerion Technologies Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.
*/
#ifndef TACS_GENERAL_SPRING_CONSTITUTIVE_H
#define TACS_GENERAL_SPRING_CONSTITUTIVE_H

#include "TACSConstitutive.h"

/**
  This is the base class for the fully general spring constitutive objects.
  Assumes 6 dofs (3 translations + 3 rotations).
  The spring properties of this object are specified using a symmetric 6 x 6
  stiffness matrix
*/
class TACSGeneralSpringConstitutive : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 6;

  TACSGeneralSpringConstitutive(TacsScalar _C[]);

  // Evaluate the material density
  virtual TacsScalar evalDensity(int elemIndex, const double pt[],
                                 const TacsScalar X[]) {
    return 0.0;
  }

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat(int elemIndex, const double pt[],
                              const TacsScalar X[]) {
    return 0.0;
  }

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void evalStress(int elemIndex, const double pt[], const TacsScalar X[],
                  const TacsScalar strain[], TacsScalar stress[]);

  // Evaluate the stiffness matrix
  // -----------------------------
  void evalTangentStiffness(int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar Ct[]);

  const char* constitutiveName();

 protected:
  TACSGeneralSpringConstitutive();

  TacsScalar C[21];

  // Calculate the stress resultants
  inline void calcStress(const TacsScalar Ct[], const TacsScalar e[],
                         TacsScalar s[]);

 private:
  static const char* constName;
};

/*
  The constitutive matrix is,

  [ Ct[ 0] Ct[ 1] Ct[ 2] Ct[ 3] Ct[ 4] Ct[ 5] ]
  [ Ct[ 1] Ct[ 6] Ct[ 7] Ct[ 8] Ct[ 9] Ct[10] ]
  [ Ct[ 2] Ct[ 7] Ct[11] Ct[12] Ct[13] Ct[14] ]
  [ Ct[ 3] Ct[ 8] Ct[12] Ct[15] Ct[16] Ct[17] ]
  [ Ct[ 4] Ct[ 9] Ct[13] Ct[16] Ct[18] Ct[19] ]
  [ Ct[ 5] Ct[10] Ct[14] Ct[17] Ct[19] Ct[20] ]
*/
inline void TACSGeneralSpringConstitutive::calcStress(const TacsScalar Ct[],
                                                      const TacsScalar e[],
                                                      TacsScalar s[]) {
  s[0] = Ct[0] * e[0] + Ct[1] * e[1] + Ct[2] * e[2] + Ct[3] * e[3] +
         Ct[4] * e[4] + Ct[5] * e[5];
  s[1] = Ct[1] * e[0] + Ct[6] * e[1] + Ct[7] * e[2] + Ct[8] * e[3] +
         Ct[9] * e[4] + Ct[10] * e[5];
  s[2] = Ct[2] * e[0] + Ct[7] * e[1] + Ct[11] * e[2] + Ct[12] * e[3] +
         Ct[13] * e[4] + Ct[14] * e[5];
  s[3] = Ct[3] * e[0] + Ct[8] * e[1] + Ct[12] * e[2] + Ct[15] * e[3] +
         Ct[16] * e[4] + Ct[17] * e[5];
  s[4] = Ct[4] * e[0] + Ct[9] * e[1] + Ct[13] * e[2] + Ct[16] * e[3] +
         Ct[18] * e[4] + Ct[19] * e[5];
  s[5] = Ct[5] * e[0] + Ct[10] * e[1] + Ct[14] * e[2] + Ct[17] * e[3] +
         Ct[19] * e[4] + Ct[20] * e[5];
}

#endif  // TACS_GENERAL_SPRING_CONSTITUTIVE_H
