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

enum TACSShellCoordinateTransform {
  TACS_NATURAL_SHELL_COORDINATES,
  TACS_REFERENCE_AXIS_COORDINATES };

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

  TACSShellConstitutive( TACSMaterialProperties *props,
                         TacsScalar _t=1.0, int _tNum=-1,
                         TacsScalar _tlb=0.0, TacsScalar _tub=1.0 );
  ~TACSShellConstitutive();

  // Set the reference axis used for the transformation
  void setRefAxis( const TacsScalar _axis[] );

  // Get the type of transformation
  enum TACSShellCoordinateTransform getTransformType(){
    return transform_type;
  }

  // Return the reference axis itself
  const TacsScalar *getRefAxis(){
    return axis;
  }

  // Get the number of stresses
  int getNumStresses();

  // Retrieve the global design variable numbers
  int getDesignVarNums( int elemIndex, int dvLen, int dvNums[] );

  // Set the element design variable from the design vector
  int setDesignVars( int elemIndex, int dvLen, const TacsScalar dvs[] );

  // Get the element design variables values
  int getDesignVars( int elemIndex, int dvLen, TacsScalar dvs[] );

  // Get the lower and upper bounds for the design variable values
  int getDesignVarRange( int elemIndex, int dvLen,
                         TacsScalar lb[], TacsScalar ub[] );

  // Evaluate the material density
  TacsScalar evalDensity( int elemIndex, const double pt[],
                          const TacsScalar X[] );

  // Add the derivative of the density
  void addDensityDVSens( int elemIndex, TacsScalar scale,
                         const double pt[], const TacsScalar X[],
                         int dvLen, TacsScalar dfdx[] );

  // Evaluate the specific heat
  TacsScalar evalSpecificHeat( int elemIndex, const double pt[],
                               const TacsScalar X[] );

  // Evaluate the stresss
  void evalStress( int elemIndex, const double pt[], const TacsScalar X[],
                   const TacsScalar strain[], TacsScalar stress[] );

  // Evaluate the tangent stiffness
  void evalTangentStiffness( int elemIndex, const double pt[],
                             const TacsScalar X[], TacsScalar C[] );

  // Add the contribution
  void addStressDVSens( int elemIndex, TacsScalar scale,
                        const double pt[], const TacsScalar X[],
                        const TacsScalar strain[], const TacsScalar psi[],
                        int dvLen, TacsScalar dfdx[] );

  // Evaluate the thermal strain
  void evalThermalStrain( int elemIndex, const double pt[],
                          const TacsScalar X[], TacsScalar theta,
                          TacsScalar strain[] );

  // Evaluate the heat flux, given the thermal gradient
  void evalHeatFlux( int elemIndex, const double pt[],
                     const TacsScalar X[], const TacsScalar grad[],
                     TacsScalar flux[] );

  // Evaluate the tangent of the heat flux
  void evalTangentHeatFlux( int elemIndex, const double pt[],
                            const TacsScalar X[], TacsScalar C[] );

  // The name of the constitutive object
  const char *getObjectName();

  // Set the drilling regularization value
  void setDrillingRegularization( double kval );

  // Extract the tangent
  static void extractTangenttStiffness( const TacsScalar *C,
                                        const TacsScalar **A, const TacsScalar **B,
                                        const TacsScalar **D, const TacsScalar **As,
                                        TacsScalar *drill );


  // Once the stiffness matrices have been evaluated, use this
  // function to compute the stress given the strain components
  inline void evalStress( const TacsScalar A[], const TacsScalar B[],
                          const TacsScalar D[], const TacsScalar As[],
                          const TacsScalar drill, const TacsScalar e[],
                          TacsScalar s[] );

 protected:
  // Material properties class
  TACSMaterialProperties *properties;

  // The drilling regularization constant
  static double DRILLING_REGULARIZATION;

  // Reference axis information
  TACSShellCoordinateTransform transform_type;
  TacsScalar axis[3]; // The reference axis

 private:
  // Store information about the design variable
  TacsScalar t, tlb, tub;
  int tNum;

  // The object name
  static const char * constName;
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
inline void TACSShellConstitutive::evalStress( const TacsScalar A[],
                                               const TacsScalar B[],
                                               const TacsScalar D[],
                                               const TacsScalar As[],
                                               const TacsScalar drill,
                                               const TacsScalar e[],
                                               TacsScalar s[] ){
  s[0] = A[0]*e[0]+A[1]*e[1]+A[2]*e[2] + B[0]*e[3]+B[1]*e[4]+B[2]*e[5];
  s[1] = A[1]*e[0]+A[3]*e[1]+A[4]*e[2] + B[1]*e[3]+B[3]*e[4]+B[4]*e[5];
  s[2] = A[2]*e[0]+A[4]*e[1]+A[5]*e[2] + B[2]*e[3]+B[4]*e[4]+B[5]*e[5];

  s[3] = B[0]*e[0]+B[1]*e[1]+B[2]*e[2] + D[0]*e[3]+D[1]*e[4]+D[2]*e[5];
  s[4] = B[1]*e[0]+B[3]*e[1]+B[4]*e[2] + D[1]*e[3]+D[3]*e[4]+D[4]*e[5];
  s[5] = B[2]*e[0]+B[4]*e[1]+B[5]*e[2] + D[2]*e[3]+D[4]*e[4]+D[5]*e[5];

  s[6] = As[0]*e[6]+As[1]*e[7];
  s[7] = As[1]*e[6]+As[2]*e[7];

  s[8] = drill*e[8];
}

#endif
