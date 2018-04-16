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

#ifndef TACS_CONSTITUTIVE_H
#define TACS_CONSTITUTIVE_H

#include <stdio.h>
#include "TACSObject.h"

/*!  
  This class defines the basic behaviour required for all
  constitutive objects in TACS

  The object defines mass, stiffness, failure and buckling behavior.
  These properties do not all need to be defined to implement the
  object.  The failure and buckling properties are not needed unless
  the object will be used by TACSFunction objects that evaluate
  failure or buckling criterion.

  The TACSConstitutive object also implements methods to visualize the
  failure and buckling envelopes. In these methods, two stress
  coordinate axes are defined and these are used to compute a 2D
  representation of the failure or buckling envelope.
*/
class TACSConstitutive : public TACSOptObject {
 public:
  TACSConstitutive(){}
  virtual ~TACSConstitutive(){}

  // Define the object name 
  // ----------------------
  virtual const char * constitutiveName() = 0;
  const char * TACSObjectName();
  
  // Return the number of stress and strain components
  // -------------------------------------------------
  virtual int getNumStresses() = 0;
  
  // Return the stress as a function of the strain at the Gauss point
  // ----------------------------------------------------------------
  virtual void calculateStress( const double pt[],
                                const TacsScalar strain[],
                                TacsScalar stress[] ) = 0;

  // Add the derivative of the stress times an input vector to dvSens
  // ----------------------------------------------------------------
  virtual void addStressDVSens( const double pt[], const TacsScalar strain[],
                                TacsScalar alpha, const TacsScalar psi[],
                                TacsScalar dvSens[], int dvLen ){}

  // Return the number of mass moments
  // ---------------------------------
  virtual int getNumMassMoments() = 0;

  // Return the mass moments
  // -----------------------
  virtual void getPointwiseMass( const double pt[],
                                 TacsScalar mass[] ) = 0;

  // Add the derivative of the pointwise mass times the given scalar
  // ---------------------------------------------------------------
  virtual void addPointwiseMassDVSens( const double pt[],
                                       const TacsScalar alpha[],
                                       TacsScalar dvSens[], int dvLen ){}

  // Evaluate the failure function at a quadrature point
  // ---------------------------------------------------
  virtual void failure( const double pt[], 
                        const TacsScalar strain[],
                        TacsScalar * fail ){ *fail = 0.0; }

  // Evaluate the derivative of the failure point w.r.t. the strain
  // --------------------------------------------------------------
  virtual void failureStrainSens( const double pt[], 
                                  const TacsScalar strain[],
                                  TacsScalar sens[] ){
    memset(sens, 0, getNumStresses()*sizeof(TacsScalar));
  }

  // Add the derivative of the failure w.r.t. design variables
  // ---------------------------------------------------------
  virtual void addFailureDVSens( const double pt[],
                                 const TacsScalar strain[],
                                 TacsScalar alpha,
                                 TacsScalar dvSens[], int dvLen ){}
  
  // Apply a buckling criterion
  // --------------------------
  virtual void buckling( const TacsScalar strain[],
                         TacsScalar *bval ){ *bval = 0.0; }

  // Compute the derivative of the buckling value w.r.t. the strain
  // --------------------------------------------------------------
  virtual void bucklingStrainSens( const TacsScalar strain[],
                                   TacsScalar *sens ){
    memset(sens, 0, getNumStresses()*sizeof(TacsScalar));
  }

  // Add the derivative of the failure w.r.t. design variables
  // ---------------------------------------------------------
  virtual void addBucklingDVSens( const TacsScalar strain[],
                                  TacsScalar alpha,
                                  TacsScalar dvSens[], int dvLen ){}

  // Return a design variable value for visualization
  // ------------------------------------------------
  virtual TacsScalar getDVOutputValue( int dvIndex, const double pt[] ){ 
    return 0.0; 
  }

  // Write out a two-dimensional representation of the failure envelope
  // ------------------------------------------------------------------
  void writeFailureEnvelope( const char *file_name, int npts,
                             const double pt[],
                             const TacsScalar x_stress[],
                             const TacsScalar y_stress[] );

  // Write out a two-dimensional buckling envelope to a file
  // -------------------------------------------------------
  void writeBucklingEnvelope( const char *file_name, int npts,
                              const double pt[],
                              const TacsScalar x_stress[],
                              const TacsScalar y_stress[],
                              double theta_min, double theta_max );
};

#endif // TACS_CONSTITUTIVE_H
