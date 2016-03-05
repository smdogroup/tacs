#ifndef TACS_CONSTITUTIVE_H
#define TACS_CONSTITUTIVE_H

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include <stdio.h>
#include "TACSObject.h"

/*!  
  This class defines the basic behaviour required for all
  constitutive relations in TACS

  The following functionality is determined by this class:

  Design variable information:
  ----------------------------  
  - Determine the design variables that are owned by this object
  - Set/get the design variables from a global list of design
  variables. The ordering of the design variables must be handled at
  the time of initialization.
  - Get the design variable ranges 

  TACSConstitutive relationship
  -------------------------
  - Calculate the stress based on the strain
  - Calculate the sensitivity of the stress to the strain

  Mass information:
  -----------------
  - Determine the pointwise mass of this object. 
  - Determine the derivative of the pointwise mass w.r.t. the design variables
  
  Failure loads:
  -------------- 
  - Determine the failure load at a point in the element due to a
  combination of constant/linear loads 

  Additional failure/buckling envelope visualization
  --------------------------------------------------
  - Make a plot in stress space of the failure envelope
  - Likewise, create a plot for the buckling envelope
*/
class TACSConstitutive : public TACSSparseConObject {
 public:
  TACSConstitutive(){}
  virtual ~TACSConstitutive(){}

  virtual const char * constitutiveName() const = 0;

  const char * TACSObjectName();
  
  // Calculate the stress and strain information
  // -------------------------------------------

  // Return the number of stress and strain components
  virtual int getNumStresses() const = 0;
  
  // Return the stress as a function of the strain at the Gauss point 
  virtual void calculateStress( const double pt[], const TacsScalar strain[], 
				TacsScalar stress[] ) = 0;

  // Return the sensitivity of the stress w.r.t. the design variable  
  virtual void calculateStressDVSens( int dvNum, const double pt[], 
				      const TacsScalar strain[], 
				      TacsScalar stress[] );

  // Add the derivative of the stress times an input vector to dvSens
  virtual void addStressDVSens( const double pt[], const TacsScalar strain[], 
				TacsScalar alpha, const TacsScalar psi[], 
				TacsScalar dvSens[], int dvLen );

  // Return information about the mass
  // ---------------------------------

  // Return the number of mass moments
  virtual int getNumMassMoments();

  // Return the mass moments
  virtual void pointwiseMass( const double pt[], TacsScalar mass[] );

  // Return the sensitivity of the mass moments
  virtual void pointwiseMassDVSens( int dvNum, const double pt[],
                                    TacsScalar massDVSens[] );

  // Add the derivative of the pointwise mass times the given scalar 
  virtual void addPointwiseMassDVSens( const double pt[], 
				       const TacsScalar alpha[],
				       TacsScalar dvSens[], int dvLen );

  /*
    Calculate a failure constraint at a Gauss point within the
    structure.  The failure constraint is applied to constant and
    linearly varying components of the strain.

    The failure criteria is that 

    fail <  1 ==> Okay
    fail >= 1 ==> Failure occurs

    pt:     The Gauss point to use - in parameter space
    fail:   The failure function
    strain: The constant componetns of the strain
  */
  virtual void failure( const double pt[], 
                        const TacsScalar strain[],
                        TacsScalar * fail );

  /*
    Determine the first derivative of the failure criteria w.r.t. the
    constant and linearly varying strain components.

    gpt:     The Gauss point 
    sens:    The derivative w.r.t. the strain components
    strain:  The components of the strain
  */
  virtual void failureStrainSens( const double pt[], 
				  const TacsScalar strain[],
                                  TacsScalar sens[] );

  /*  
    Determine the sensitivity of the failure criteria to the design
    variables.

    dvNum:    The design variable number
    pt:       The Gauss point to evaluate the failure criteria 
    strain:   The components of the strain
    failSens: The derivative of the failure criteria
  */
  virtual void failureDVSens( int dvNum, const double pt[], 
			      const TacsScalar strain[],
                              TacsScalar * failSens );

  /*
    Add the derivative of the failure function with respect to the
    design variables to the input vector

    pt:       The Gauss point to evaluate the failure criteria 
    strain:   The components of the strain
    alpha:    Scalar multiple for the derivative
    dvSens:   The vector in which to add the derivative
    dvLen:    The length of the design variable vector
  */
  virtual void addFailureDVSens( const double pt[], 
				 const TacsScalar strain[],
				 TacsScalar alpha,
				 TacsScalar dvSens[], int dvLen );
  
  /*
    Apply an appropriate buckling constraint based on the provided
    strain values. Not all constitutive classes need to provide a 
    buckling constraint.
  
    buckling: the normalized buckling criteria such that 
    bval < 1.0: no buckling
    bval >= 1.0: buckling will occur

    strain: the value of the strain to be used within the buckling
    criteria
  */
  virtual void buckling( const TacsScalar strain[],
                         TacsScalar * bval );

  /*
    The derivative of the buckling calculation w.r.t. the strain
    components.

    bvalSens: the sensitivity of the buckling calculation to the
    strain at the current point

    strain: the strain used for the buckling calculation
  */
  virtual void bucklingStrainSens( const TacsScalar strain[],
                                   TacsScalar * bvalSens );

  /*
    The derivative of the buckling calculation w.r.t. the design
    variables owned by the constitutive property.

    bvalSens: the sensitivity of the buckling calculation to the
    set design variable
    
    strain: the strain used for the buckling calculation
  */
  virtual void bucklingDVSens( int dvNum, 
                               const TacsScalar strain[],
                               TacsScalar * bvalSens );

  /*
    Add the derivative of the buckling criteria with respect to the
    design variable

    pt:       The Gauss point to evaluate the failure criteria 
    strain:   The components of the strain
    alpha:    Scalar multiple for the derivative
    dvSens:   The vector in which to add the derivative
    dvLen:    The length of the design variable vector
  */
  virtual void addBucklingDVSens( const TacsScalar strain[],
				  TacsScalar alpha,
				  TacsScalar dvSens[], int dvLen );

  /*
    Return a design variable value for visualization purposes.
    This only provides a scalar at the given Gauss point.
  */
  virtual TacsScalar getDVOutputValue( int dvIndex, 
                                       const double pt[] ){ return 0.0; }

  /*
    Write out a two-dimensional representation of the failure envelope
    to a file
  */
  void writeFailureEnvelope( const char * file_name, int npts,
			     const double pt[], 
			     const TacsScalar x_stress[],
			     const TacsScalar y_stress[] );

  /*
    Write out a two-dimensional buckling envelope to a file
  */
  void writeBucklingEnvelope( const char * file_name, int npts,
			      const double pt[], 
			      const TacsScalar x_stress[],
			      const TacsScalar y_stress[],
			      double theta_min, double theta_max );
};

# endif
