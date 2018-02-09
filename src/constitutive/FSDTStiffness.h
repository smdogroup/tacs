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

#ifndef FSDT_STIFFNESS_H
#define FSDT_STIFFNESS_H

#include "TACSConstitutive.h"

/*!  
  This constitutive class defines the stiffness properties for a
  first-order shear deformation theory type element. This class
  is derived from the TACSConstitutive object, but is still
  a pure virtual base class. 
*/
class FSDTStiffness : public TACSConstitutive {
 public:
  // Set the number of stresses and number of mass moments
  static const int NUM_STRESSES = 8;
  static const int NUM_MASS_MOMENTS = 2;

  // The type of local transformation to apply at each point 
  // -------------------------------------------------------
  enum FSDTTransformType { NATURAL, REFERENCE_AXIS };
  
  // The default constructor for the base FSDTStiffness class
  // --------------------------------------------------------
  FSDTStiffness();
  virtual ~FSDTStiffness(){}

  // Return the number of mass stress and number of mass moments
  // -----------------------------------------------------------
  int getNumStresses();
  int getNumMassMoments();

  // Set the drilling stiffness regularization parameter
  // This modifies the parameter for all FSDT stiffness objects!
  // -----------------------------------------------------------
  void setDrillingRegularization( double kval );

  // Set the reference axis used for the transformation
  // --------------------------------------------------
  void setRefAxis( const TacsScalar _axis[] ){
    transform_type = REFERENCE_AXIS;
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];
  }

  // Return the type of transformation
  // ---------------------------------
  enum FSDTTransformType getTransformType(){ 
    return transform_type;
  }

  // Return the reference axis itself
  // --------------------------------
  const TacsScalar *getRefAxis(){ 
    return axis;
  }

  // The name of the constitutive object
  // -----------------------------------
  const char *constitutiveName();

  // Evaluate the stiffness matrix at the given parametric point
  // -----------------------------------------------------------
  virtual TacsScalar getStiffness( const double pt[], 
                                   TacsScalar A[], TacsScalar B[],
                                   TacsScalar D[], TacsScalar As[] ) = 0;
  
  // Add the derivative of the product of the stiffness matrix with the vectors
  // --------------------------------------------------------------------------
  virtual void addStiffnessDVSens( const double pt[],
                                   const TacsScalar e[], const TacsScalar psi[],
                                   TacsScalar rotPsi, 
                                   TacsScalar fdvSens[], int dvLen ){}

  // Once the stiffness matrices have been evaluated, use this
  // function to compute the stress given the strain components
  // ----------------------------------------------------------
  inline void calculateStress( const TacsScalar A[], const TacsScalar B[],
                               const TacsScalar D[], const TacsScalar As[],
                               const TacsScalar e[], TacsScalar s[] );

  // The following required TACSConstitutive functions are implemented
  // -----------------------------------------------------------------
  void calculateStress( const double pt[], 
                        const TacsScalar strain[], 
                        TacsScalar stress[] );
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
                        TacsScalar alpha, const TacsScalar psi[], 
                        TacsScalar dvSens[], int dvLen );

  // The following function prints the stiffness 
  // -------------------------------------------
  void printStiffness( const double pt[]=NULL );
  
 protected:
  // The drilling regularization constant
  static double DRILLING_REGULARIZATION;

  // Reference axis information
  enum FSDTTransformType transform_type;
  TacsScalar axis[3]; // The reference axis

 private:
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
inline void FSDTStiffness::calculateStress( const TacsScalar A[], 
                                            const TacsScalar B[],
                                            const TacsScalar D[], 
                                            const TacsScalar As[],
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
}

#endif
