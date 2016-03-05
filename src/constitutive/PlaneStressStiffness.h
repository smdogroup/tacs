#ifndef PLANE_STRESS_STIFFNESS_H
#define PLANE_STRESS_STIFFNESS_H

/*
  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSConstitutive.h"

/*
  This is the base class for the plane stress constitutive objects. 
  
  All objects performing plane stress analysis should utilize this class. 
*/
class PlaneStressStiffness : public TACSConstitutive {
 public:
  static const int NUM_STRESSES = 3;
  PlaneStressStiffness( TacsScalar _rho, TacsScalar E, TacsScalar nu );
  PlaneStressStiffness( TacsScalar _rho, TacsScalar E1, 
			TacsScalar E2, TacsScalar G12, TacsScalar nu12 );
  virtual ~PlaneStressStiffness(){}

  // Calculate the stress
  // --------------------
  int getNumStresses() const;
  void calculateStress( const double pt[], const TacsScalar strain[],
			TacsScalar stress[] );

  // Return the mass moments
  // -----------------------
  int getNumMassMoments(){ return 1; }
  void pointwiseMass( const double pt[], TacsScalar mass[] ){
    mass[0] = rho;
  }

  // Extra info about the constitutive class
  // ---------------------------------------
  const char * constitutiveName() const;

 protected:
  PlaneStressStiffness();

  // The stiffness matrix
  TacsScalar Cmat[6]; 

  TacsScalar rho;
 private:
  static const char * constName;
};

#endif
