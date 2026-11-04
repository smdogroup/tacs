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

#include "TACSAuxElements.h"

#include <vector>

static int compare_elems(const void *a, const void *b) {
  const TACSAuxElem *ao = static_cast<const TACSAuxElem *>(a);
  const TACSAuxElem *bo = static_cast<const TACSAuxElem *>(b);

  return ao->num - bo->num;
}

const char *TACSAuxElements::auxName = "TACSAuxElements";

/*
  Initialize the auxiliary element object

  input:
  num_elems:  an estimate of the number of elements that will be added
*/
TACSAuxElements::TACSAuxElements(int _num_elems) {
  max_elements = (_num_elems < 100 ? 100 : _num_elems);
  aux = new TACSAuxElem[max_elements];
  num_elements = 0;
  is_sorted = 0;
}

/*
  Free all the auxiliary element data
*/
TACSAuxElements::~TACSAuxElements() {
  for (int i = 0; i < num_elements; i++) {
    aux[i].elem->decref();
  }
  delete[] aux;
}

/*
  Sort the list of auxiliary elements so that their numbers are listed
  in ascending order. This is required within the TACSAssembler object
  for efficient residual assembly operations.
*/
void TACSAuxElements::sort() {
  if (!is_sorted) {
    qsort(aux, num_elements, sizeof(TACSAuxElem), compare_elems);
    is_sorted = 1;
  }
}

/*
  Add an element to the auxiliary element list

  This function can be called at anytime after the TACSAssembler
  function has been initilized. All elements added here must occupy
  the same non-zero pattern as those elements they are matched with in
  TACSAssembler. In essence, the number of element variables for each
  element must be identical.

  input:
  num:   the TACSAssembler element number
  elem:  the TACSElement pointer
*/
void TACSAuxElements::addElement(int num, TACSElement *elem) {
  // The array is not large enough to handle a new element
  // Create a new array and copy over the values
  if (num_elements >= max_elements) {
    max_elements *= 2;
    TACSAuxElem *tmp = new TACSAuxElem[max_elements];
    memcpy(tmp, aux, num_elements * sizeof(TACSAuxElem));
    delete[] aux;
    aux = tmp;
  }

  // Insert the new element into the list
  elem->incref();
  aux[num_elements].elem = elem;
  aux[num_elements].num = num;
  num_elements++;

  // The list is not sorted anymore
  is_sorted = 0;
}

/*
  Add an entire list of elements

  This function adds a number of TACSElement objects to the auxiliary
  list of elements all at the same time. The code automatically
  extends the internally stored array if it is not long enough.

  input:
  nums:        an array of the TACSAssembler element numbers
  elem:        an array of axuiliary elements
  num_elems:   the number of elements in the arrays
*/
void TACSAuxElements::addElements(int nums[], TACSElement **elem,
                                  int num_elems) {
  // The array is not large enough to handle a new element
  // Create a new array and copy over the values
  if (num_elements + num_elems >= max_elements) {
    max_elements = 2 * max_elements + num_elems;
    TACSAuxElem *tmp = new TACSAuxElem[max_elements];
    memcpy(tmp, aux, num_elements * sizeof(TACSAuxElem));
    delete[] aux;
    aux = tmp;
  }

  // Insert the new element into the list
  for (int k = 0; k < num_elems; k++) {
    elem[k]->incref();
    aux[num_elements].elem = elem[k];
    aux[num_elements].num = nums[k];
    num_elements++;
  }

  // The list is not sorted anymore
  is_sorted = 0;
}

/*
  Get the elements and sort them (if they are not already)
*/
int TACSAuxElements::getAuxElements(TACSAuxElem **_elems) {
  *_elems = aux;
  return num_elements;
}

/*
  Set the design variables in all auxiliary elements from a design vector

  This is the only implementation of the auxiliary element design variable
  update: TACSAssembler::setDesignVars uses it for the auxiliary elements that
  are attached to the assembler. Problems that own several TACSAuxElements
  objects (for instance one per time step of a transient problem) call it on
  each object directly.

  This call is collective, because it distributes the non-local design
  variable values. Every rank must make it, even a rank with no auxiliary
  elements.

  Each design variable number in the design vector has a block of
  dvs->getBlockSize() values (the designVarsPerNode of the assembler), and the
  values of an element's design variables are passed to the element block by
  block. The load elements (pressure, traction, inertial and centrifugal
  forces) use one value per design variable. With more than one value per
  design variable, a load element with a single design variable therefore
  gets the first value of its block, and its sensitivity goes into the first
  value. A load element with several design variables reads its values from
  consecutive entries that do not line up with the blocks, so it is not
  supported with more than one value per design variable.

  input:
  dvs:  the design variable vector
*/
void TACSAuxElements::setDesignVars(TACSBVec *dvs) {
  dvs->beginDistributeValues();
  dvs->endDistributeValues();
  setDistributedDesignVars(dvs);
}

/*
  Set the design variables in all auxiliary elements from a design vector
  whose non-local values have already been distributed

  This call is not collective.

  input:
  dvs:  the distributed design variable vector
*/
void TACSAuxElements::setDistributedDesignVars(TACSBVec *dvs) {
  const int valuesPerDV = dvs->getBlockSize();
  std::vector<int> dvNums;
  std::vector<TacsScalar> dvVals;

  for (int i = 0; i < num_elements; i++) {
    TACSElement *elem = aux[i].elem;
    const int elemNum = aux[i].num;
    const int numDVs = elem->getDesignVarNums(elemNum, 0, NULL);
    if (numDVs > 0) {
      dvNums.resize(numDVs);
      dvVals.resize(valuesPerDV * numDVs);
      elem->getDesignVarNums(elemNum, numDVs, dvNums.data());
      const int fail = dvs->getValues(numDVs, dvNums.data(), dvVals.data());
      if (fail != 0) {
        fprintf(stderr,
                "TACSAuxElements: aux element %d uses a design variable that "
                "is not in the external design variable list of this rank\n",
                elemNum);
        continue;
      }
      elem->setDesignVars(elemNum, numDVs, dvVals.data());
    }
  }
}
