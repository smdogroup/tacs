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
  dv_buffer_size = 0;
  dv_nums = NULL;
  dv_vals = NULL;
}

/*
  Free all the auxiliary element data
*/
TACSAuxElements::~TACSAuxElements() {
  for (int i = 0; i < num_elements; i++) {
    aux[i].elem->decref();
  }
  delete[] aux;
  delete[] dv_nums;
  delete[] dv_vals;
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

  This performs the same work as the auxiliary element loop in
  TACSAssembler::setDesignVars, but on this object alone. It is needed for
  problems that own several TACSAuxElements objects (for instance one per time
  step of a transient problem) and therefore cannot rely on the single
  TACSAuxElements instance currently attached to the assembler.

  input:
  dvs:  the design variable vector
*/
void TACSAuxElements::setDesignVars(TACSBVec *dvs) {
  if (num_elements == 0) {
    return;
  }

  // Distribute the non-local design variable values
  dvs->beginDistributeValues();
  dvs->endDistributeValues();

  // Make sure the scratch arrays are large enough to hold the design variables
  // of any one of the auxiliary elements
  int maxDVs = 0;
  for (int i = 0; i < num_elements; i++) {
    int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, 0, NULL);
    if (numDVs > maxDVs) {
      maxDVs = numDVs;
    }
  }
  if (maxDVs > dv_buffer_size) {
    delete[] dv_nums;
    delete[] dv_vals;
    dv_buffer_size = maxDVs;
    dv_nums = new int[dv_buffer_size];
    dv_vals = new TacsScalar[dv_buffer_size];
  }

  for (int i = 0; i < num_elements; i++) {
    int numDVs =
        aux[i].elem->getDesignVarNums(aux[i].num, dv_buffer_size, dv_nums);
    if (numDVs > 0) {
      dvs->getValues(numDVs, dv_nums, dv_vals);
      aux[i].elem->setDesignVars(aux[i].num, numDVs, dv_vals);
    }
  }
}
