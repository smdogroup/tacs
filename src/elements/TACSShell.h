#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

/*
  Copyright (c) 2015 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSElement.h"
#include "FSDTStiffness.h"

/*
  The following class defines the base class for all Shell elements
  used in TACS. This class defines some of the operations required 
  by the generic TACSElement base class.
*/
class TACSShell : public TACSElement {
 public:
  // Define some constants for this element type
  static const int NUM_DISPS = 6;
  static const int NUM_STRESSES = 8;
  static const int NUM_EXTRAS = 4;

  TACSShell( FSDTStiffness * _stiff,
	     int componentNum );
  ~TACSShell();

  // Functions to access the design variable information
  // ---------------------------------------------------
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lowerBound[], 
			  TacsScalar upperBound[], int numDVs );

  // Functions to determine the variable names and quantities
  // --------------------------------------------------------
  const char * elementName();
  const char * displacementName( int i );
  const char * stressName( int i );
  const char * strainName( int i );
  const char * extraName( int i );
  int numDisplacements();
  int numStresses();
  int numExtras();
  ElementType getElementType();

 protected:
  // The constitutive object pointer
  FSDTStiffness *stiff;

 private:
  // The names of the displacements, stresses etc.
  static const char * elemName;
  static const char * dispNames[6];
  static const char * stressNames[8];
  static const char * strainNames[8];
  static const char * extraNames[4];
};

#endif
