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

#ifndef TACS_AUX_ELEMENTS_H
#define TACS_AUX_ELEMENTS_H

/*
  TACSAuxElements implementation
*/

#include "TACSElement.h"
#include "TACSObject.h"

/*
  The following class defines a single auxiliary element and
  its associated element number.

  This class does not need to be instantiated by the user - this is a
  class used within TACSAuxElements/TACSAssembler to facilitate the
  sorting/searching of auxiliary elements.
*/
class TACSAuxElem {
 public:
  TACSElement *elem;
  int num;
};

/*
  The TACSAuxiliaryElements class

  This class provides a way to add extra elements - often defining
  body-loads, tractions etc. for elements within the TACSAssembler
  class. These elements are restricted to have the exact same non-zero
  pattern as the existing elements that are set in the TACSAssembler
  object.
*/
class TACSAuxElements : public TACSObject {
 public:
  TACSAuxElements(int _num_elems = 100);
  ~TACSAuxElements();

  // Sort the list of elements
  // -------------------------
  void sort();

  // Add a single element
  // --------------------
  void addElement(int num, TACSElement *elem);

  // Add an entire list of elements
  // ------------------------------
  void addElements(int nums[], TACSElement **elem, int num_elems);

  // Get the elements and sort them (if they are not already)
  // --------------------------------------------------------
  int getAuxElements(TACSAuxElem **elems);

  // Functions to control the design variables
  // -----------------------------------------
  void getDesignVars(int numDVs, TacsScalar dvs[]);
  void setDesignVars(int numDVs, const TacsScalar dvs[]);
  void getDesignVarRange(int numDVs, TacsScalar lb[], TacsScalar ub[]);

  // Print the name of the TACSObject
  // --------------------------------
  const char *TACSObjectName() { return auxName; }

 private:
  // Keep track of whether the element list has been sorted
  int is_sorted;

  // The maxinum size of the internal array - this will be expanded
  // automatically as more elements are added
  int max_elements;

  // Keep track of the number of added elements
  int num_elements;

  // The auxiliary elements
  TACSAuxElem *aux;

  // The auxiliary object name
  static const char *auxName;
};

#endif  // TACS_AUXILIARY_ELEMENTS_H
