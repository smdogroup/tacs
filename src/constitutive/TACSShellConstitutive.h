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

#ifndef TACS_SHELL_CONSTITUTIVE_H
#define TACS_SHELL_CONSTITUTIVE_H

#include "TACSConstitutive.h"
#include "TACSMaterialProperties.h"

/**
  This constitutive class defines the stiffness properties for a
  first-order shear deformation theory type element. This class
  is derived from the TACSConstitutive object, but is still
  a pure virtual base class.
*/
class TACSShellConstitutive : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 9;
  static const int NUM_TANGENT_STIFFNESS_ENTRIES = 22;

  TACSShellConstitutive() {}
  virtual ~TACSShellConstitutive() {}

  // Get the number of stresses
  int getNumStresses();

  /*
    Get the integrals of the density through the thickness

    moments = int_{-t/2}^{t/2} [1, z, z^2] rho(z) dz

    @param elemIndex The local element index
    @param pt The parametric location
    @param X The point location
    @return The moments of the mass
  */
  virtual void evalMassMoments(int elemIndex, const double pt[],
                               const TacsScalar X[], TacsScalar moments[]) = 0;

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

  // Set the drilling regularization value
  static void setDrillingRegularization(double kval);

  // Extract the tangent
  static void extractTangentStiffness(const TacsScalar *C, const TacsScalar **A,
                                      const TacsScalar **B,
                                      const TacsScalar **D,
                                      const TacsScalar **As, TacsScalar *drill);

  // Once the stiffness matrices have been evaluated, use this
  // function to compute the stress given the strain components
  static inline void computeStress(const TacsScalar A[], const TacsScalar B[],
                                   const TacsScalar D[], const TacsScalar As[],
                                   const TacsScalar drill, const TacsScalar e[],
                                   TacsScalar s[]);

  // The name of the constitutive object
  const char *getObjectName();

 protected:
  // The drilling regularization constant
  static double DRILLING_REGULARIZATION;

 private:
  // The object name
  static const char *constName;
};

/*
  Given the stiffness matrices, compute the stress based on the values
  of the strain. This computation takes the form:

  [ s[0:3] ] = [ A B 0  ][ e[0:3] ]
  [ s[3:6] ] = [ B D 0  ][ e[3:6] ]
  [ s[6:8] ] = [ 0 0 As ][ e[6:8] ]

  Each matrix in the ABD matrix is symmetric and is stored as follows:

  [A] = [ A[0] A[1] A[2] ]
  .     [ A[1] A[3] A[4] ]
  .     [ A[2] A[4] A[5] ]

  The shear stiffness matrix takes the form:

  [As] = [ As[0] As[1] ]
  .      [ As[1] As[2] ]
*/
inline void TACSShellConstitutive::computeStress(
    const TacsScalar A[], const TacsScalar B[], const TacsScalar D[],
    const TacsScalar As[], const TacsScalar drill, const TacsScalar e[],
    TacsScalar s[]) {
  s[0] = A[0] * e[0] + A[1] * e[1] + A[2] * e[2] + B[0] * e[3] + B[1] * e[4] +
         B[2] * e[5];
  s[1] = A[1] * e[0] + A[3] * e[1] + A[4] * e[2] + B[1] * e[3] + B[3] * e[4] +
         B[4] * e[5];
  s[2] = A[2] * e[0] + A[4] * e[1] + A[5] * e[2] + B[2] * e[3] + B[4] * e[4] +
         B[5] * e[5];

  s[3] = B[0] * e[0] + B[1] * e[1] + B[2] * e[2] + D[0] * e[3] + D[1] * e[4] +
         D[2] * e[5];
  s[4] = B[1] * e[0] + B[3] * e[1] + B[4] * e[2] + D[1] * e[3] + D[3] * e[4] +
         D[4] * e[5];
  s[5] = B[2] * e[0] + B[4] * e[1] + B[5] * e[2] + D[2] * e[3] + D[4] * e[4] +
         D[5] * e[5];

  s[6] = As[0] * e[6] + As[1] * e[7];
  s[7] = As[1] * e[6] + As[2] * e[7];

  s[8] = drill * e[8];
}

#endif
