#ifndef FSDT_STIFFNESS_H
#define FSDT_STIFFNESS_H

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSConstitutive.h"

/*!
  This class defines the stiffness properties for an FSDT element
  - it is derived from the Constitutive pure virtual base class
*/
class FSDTStiffness : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 8;

  enum FSDTTransformType { NATURAL, REFERENCE_AXIS };
  
  FSDTStiffness();
  virtual ~FSDTStiffness(){}

  // Set the type of transformation 
  // ------------------------------
  void setRefAxis( const TacsScalar _axis[] ){
    transform_type = REFERENCE_AXIS;
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];
  }

  // Return the type of transformation
  // ---------------------------------
  enum FSDTTransformType getTransformType() const { 
    return transform_type; 
  }

  // Return the reference axis itself
  // --------------------------------
  const TacsScalar * getRefAxis(){ 
    return axis; 
  }

  // Does the stiffness vary parametrically in the element?
  // Note that this defaults to true.
  // ------------------------------------------------------
  int isVariableStiffness(){ return is_param_variable; }

  // The standard methods required from the constitutive class
  // ---------------------------------------------------------
  
  // Calculate the stress
  // --------------------
  int getNumStresses() const;
  void calculateStress( const double gpt[], const TacsScalar strain[],
                        TacsScalar stress[] );
  void calculateStressDVSens( int dvNum, const double gpt[],
                              const TacsScalar strain[], TacsScalar stress[] );
  void addStressDVSens( const double pt[], const TacsScalar strain[], 
			TacsScalar alpha, const TacsScalar psi[], 
			TacsScalar dvSens[], int dvLen );

  // Retrieve the stiffness information for this object
  // --------------------------------------------------
  virtual TacsScalar getStiffness( const double gpt[],
                                   TacsScalar A[], TacsScalar B[],
                                   TacsScalar D[], TacsScalar As[] ) = 0;
  virtual TacsScalar getStiffnessDVSens( int dvNum, const double gpt[],
                                         TacsScalar sA[], TacsScalar sB[],
                                         TacsScalar sD[], TacsScalar sAs[] );
  virtual void addStiffnessDVSens( const double pt[], TacsScalar alpha, 
				   const TacsScalar psi[], const TacsScalar phi[], 
				   const TacsScalar psi_rot, const TacsScalar phi_rot,
				   TacsScalar dvSens[], int dvLen );

  // Based on the provided values of A, B, D, As, compute the the
  // stress s[] based on the values of the strain e[]
  // ------------------------------------------------------------
  inline void calcStress( const TacsScalar A[], const TacsScalar B[],
                          const TacsScalar D[], const TacsScalar As[],
                          const TacsScalar e[], TacsScalar s[] );

  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 3; }
  void pointwiseMassDVSens( int dvNum, const double gpt[],
                            TacsScalar massDVSens[] );
  void addPointwiseMassDVSens( const double pt[], 
			       const TacsScalar alpha[],
			       TacsScalar dvSens[], int dvLen );
  virtual TacsScalar getDensity(){ return 0.0; }

  //  These are the default implementations
  // --------------------------------------
  void addFailureDVSens( const double pt[], const TacsScalar strain[],
			 TacsScalar alpha, TacsScalar dvSens[], int dvLen );
  void addBucklingDVSens( const TacsScalar strain[], TacsScalar alpha,
			  TacsScalar dvSens[], int dvLen );

  // Set the drilling stiffness regularization parameter
  // This modifies the parameter for all FSDT stiffness objects!
  // -----------------------------------------------------------
  void setDrillingRegularization( double kval );

  // Extra info about the constitutive class
  // ---------------------------------------
  void printInfo();
  const char * constitutiveName() const;

  // These are functions for gradient-based design optimization. They
  // are designed to be used directly from the object, not through the
  // regular TACS interface. 
  // -----------------------------------------------------------------
  TacsScalar computeFailure( const double pt[], const TacsScalar s[] );
  TacsScalar computeBuckling( const double pt[], const TacsScalar s[] );
  void computeFailureDVSens( const double pt[], const TacsScalar s[], 
                             TacsScalar fdvSens[], int numDVs );
  void computeBucklingDVSens( const double pt[], const TacsScalar s[],
                              TacsScalar fdvSens[], int numDVs );

 protected:
  // Does this stiffness object vary parametrically within an element?
  int is_param_variable;
  
  // The drilling regularization constant
  static double DRILLING_REGULARIZATION;

  // Reference axis information
  enum FSDTTransformType transform_type;
  TacsScalar axis[3]; // The reference axis

 private:
  void computeStrainFromStress( const double pt[], const TacsScalar s[],
                                TacsScalar e[] );
  void computeStrainFromStressDVSens( int dvNum, 
                                      const double pt[], 
                                      const TacsScalar e[],
                                      TacsScalar es[] );

  static const char * constName;
};

/*
  Calculate the stress based on the strain
  
  s = [ A B 0 ]
  .   [ B D 0 ]
  .   [ 0 0 As]

  [A] = [ A[0] A[1] A[2] ]
  .     [ A[1] A[3] A[4] ]
  .     [ A[2] A[4] A[5] ]

  Same with B, D
*/
inline void FSDTStiffness::calcStress( const TacsScalar A[], 
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
