/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#ifndef EB_BEAM_STIFFNESS
#define EB_BEAM_STIFFNESS

/*
  The stiffness object for the variety of different beams
*/

#include "TACSConstitutive.h"

/*
  The following class defines the stiffness for Euler--Bernoulli
  beam theory

  The strong axis/weak axis flag indicate how the reference direction
  is used within the Euler--Bernoulli beam elements. The strong axis
  indicates that the component of the axis orthogonal to the beam
  direction is the first x-direction. The weak axis indicates that
  the orthogonal component is the y-axis.
*/
class EBStiffness : public TACSConstitutive {
 public:
  enum EBReferenceDirection { STRONG_AXIS, WEAK_AXIS };
  static const int NUM_STRESSES = 4;

  EBStiffness(TacsScalar rho, TacsScalar E, TacsScalar G, TacsScalar A,
              TacsScalar Ix, TacsScalar Iy, TacsScalar J,
              TacsScalar _ref_dir[3],
              EBReferenceDirection _ref_type = WEAK_AXIS);
  virtual ~EBStiffness() {}

  // The reference direction of the beam
  // -----------------------------------
  const EBReferenceDirection ref_type;
  const TacsScalar *getRefDir() { return ref_dir; }

  // Retrieve the stiffness
  // ----------------------
  virtual void getStiffness(const double pt[], TacsScalar Ct[]);

  // Calculate the stress
  // --------------------
  int getNumStresses();
  void calculateStress(const double pt[], const TacsScalar strain[],
                       TacsScalar stress[]);
  void addStressDVSens(const double pt[], const TacsScalar strain[],
                       TacsScalar alpha, const TacsScalar psi[],
                       TacsScalar fdvSens, int numDVs);

  // Return the mass moments
  // -----------------------
  int getNumMassMoments();
  void getPointwiseMass(const double pt[], TacsScalar mass[]);
  void addPointwiseMassDVSens(const double pt[], const TacsScalar alpha[],
                              TacsScalar dvSens[], int dvLen);

  const char *constitutiveName();

 protected:
  EBStiffness();
  TacsScalar ref_dir[3];
  TacsScalar mass[6];
  TacsScalar C[10];

  // Calculate the stress resultants
  inline void calcStress(const TacsScalar Ct[], const TacsScalar e[],
                         TacsScalar s[]);

 private:
  static const char *constName;
};

/*
  The constitutive matrix is,

  [ Ct[0] Ct[1] Ct[2] Ct[3] ]
  [ Ct[1] Ct[4] Ct[5] Ct[6] ]
  [ Ct[2] Ct[5] Ct[7] Ct[8] ]
  [ Ct[3] Ct[6] Ct[8] Ct[9] ]
*/
inline void EBStiffness::calcStress(const TacsScalar Ct[], const TacsScalar e[],
                                    TacsScalar s[]) {
  s[0] = Ct[0] * e[0] + Ct[1] * e[1] + Ct[2] * e[2] + Ct[3] * e[3];
  s[1] = Ct[1] * e[0] + Ct[4] * e[1] + Ct[5] * e[2] + Ct[6] * e[3];
  s[2] = Ct[2] * e[0] + Ct[5] * e[1] + Ct[7] * e[2] + Ct[8] * e[3];
  s[3] = Ct[3] * e[0] + Ct[6] * e[1] + Ct[8] * e[2] + Ct[9] * e[3];
}

#endif
