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

#ifndef TACS_SHELL_ELEMENT_H
#define TACS_SHELL_ELEMENT_H

#include "FSDTStiffness.h"
#include "TACSElement.h"

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

  TACSShell(FSDTStiffness *_stiff, int componentNum);
  ~TACSShell();

  // Functions to access the design variable information
  // ---------------------------------------------------
  void setDesignVars(const TacsScalar dvs[], int numDVs);
  void getDesignVars(TacsScalar dvs[], int numDVs);
  void getDesignVarRange(TacsScalar lowerBound[], TacsScalar upperBound[],
                         int numDVs);

  // Functions to determine the variable names and quantities
  // --------------------------------------------------------
  const char *elementName();
  const char *displacementName(int i);
  const char *stressName(int i);
  const char *strainName(int i);
  const char *extraName(int i);
  int numDisplacements();
  int numStresses();
  int numExtras();
  ElementType getElementType();

 protected:
  // The constitutive object pointer
  FSDTStiffness *stiff;

 private:
  // The names of the displacements, stresses etc.
  static const char *elemName;
  static const char *dispNames[6];
  static const char *stressNames[8];
  static const char *strainNames[8];
  static const char *extraNames[4];
};

#endif
