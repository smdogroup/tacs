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

#include "TACSShell.h"

TACSShell::TACSShell(FSDTStiffness* _stiff, int componentNum)
    : TACSElement(componentNum) {
  stiff = _stiff;
  stiff->incref();
}

TACSShell::~TACSShell() { stiff->decref(); }

/*
  Set up the internal static data for the names of the element,
  displacements, stresses, strains and extra variables, respectively.
*/
const char* TACSShell::elemName = "TACSShell";

const char* TACSShell::dispNames[] = {"u0", "v0", "w0", "rotx", "roty", "rotz"};

const char* TACSShell::stressNames[] = {"sx0", "sy0",  "sxy0", "sx1",
                                        "sy1", "sxy1", "syz0", "sxz0"};

const char* TACSShell::strainNames[] = {"ex0", "ey0",  "exy0", "ex1",
                                        "ey1", "exy1", "eyz0", "exz0"};

const char* TACSShell::extraNames[] = {"lambda", "buckling", "dv1", "dv2"};

/*
  Define the functions that return the element names, displacements,
  stresses, strains and extra variables.
*/
const char* TACSShell::elementName() { return elemName; }

const char* TACSShell::displacementName(int i) {
  if (i >= 0 && i < NUM_DISPS) {
    return dispNames[i];
  }
  return NULL;
}

const char* TACSShell::stressName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return stressNames[i];
  }
  return NULL;
}

const char* TACSShell::strainName(int i) {
  if (i >= 0 && i < NUM_STRESSES) {
    return strainNames[i];
  }
  return NULL;
}

const char* TACSShell::extraName(int i) {
  if (i >= 0 && i < NUM_EXTRAS) {
    return extraNames[i];
  }
  return NULL;
}

/*
  Define the functions that return the number of displacements, stresses,
  nodes, variables and extra variables
*/
int TACSShell::numDisplacements() { return NUM_DISPS; }

int TACSShell::numStresses() { return NUM_STRESSES; }

int TACSShell::numExtras() { return NUM_EXTRAS; }

ElementType TACSShell::getElementType() { return TACS_SHELL; }

/*
  Set the values of the design variables
*/
void TACSShell::setDesignVars(const TacsScalar dvs[], int numDVs) {
  stiff->setDesignVars(dvs, numDVs);
}

/*
  Get the design variable values -- populate the array dvs[]
*/
void TACSShell::getDesignVars(TacsScalar dvs[], int numDVs) {
  stiff->getDesignVars(dvs, numDVs);
}

/*
  Populate the arrays lowerBound[] and upperBound[]
*/
void TACSShell::getDesignVarRange(TacsScalar lowerBound[],
                                  TacsScalar upperBound[], int numDVs) {
  stiff->getDesignVarRange(lowerBound, upperBound, numDVs);
}
