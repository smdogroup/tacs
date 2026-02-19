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

#include "TACSFunction.h"

#include "TACSAssembler.h"
#include "TacsUtilities.h"

/*
  Base TACSFunction implementation
*/

/*
  Create the function base class with the default domain - typically
  the entire finite-element mesh.

  input:
  assembler:     the TACSAssembler object
  funcDomain:    the domain type
  funcEval:      the type of evaluation to use
  maxElems:      the maximum number of elements expected
*/
TACSFunction::TACSFunction(TACSAssembler *_assembler, DomainType _funcDomain,
                           StageType _funcStages, int _maxElems) {
  assembler = _assembler;
  assembler->incref();

  // Set function domain and function evaluation type
  funcDomain = _funcDomain;
  funcStageType = _funcStages;

  // Set the domain elements
  maxElems = (_maxElems > 0 ? _maxElems : 0);
  numElems = 0;
  elemNums = NULL;
}

/*
  Destroy the TACSFunction base class
*/
TACSFunction::~TACSFunction() {
  if (elemNums) {
    delete[] elemNums;
  }
  assembler->decref();
}

/*
  Retrieve the type of domain specified by this object
*/
enum TACSFunction::DomainType TACSFunction::getDomainType() {
  return funcDomain;
}

/*
  Retrieve the type of function
*/
enum TACSFunction::StageType TACSFunction::getStageType() {
  return funcStageType;
}

/*
  Overwrite the domain in the function with a new set of elements.
  This reallocates the existing array if it is not long enough.

  input:
  elemNums: the element numbers used to set the domain
  numElems: the number of elements to add
*/
void TACSFunction::setDomain(int _numElems, const int _elemNums[]) {
  if (funcDomain == NO_DOMAIN) {
    fprintf(stderr, "Cannot set function domain for %s\n", getObjectName());
    return;
  } else {
    funcDomain = SUB_DOMAIN;

    if (_numElems > maxElems) {
      if (elemNums) {
        delete[] elemNums;
      }
      numElems = _numElems;
      maxElems = _numElems;
      elemNums = new int[maxElems];
    } else {
      numElems = _numElems;
    }

    memcpy(elemNums, _elemNums, numElems * sizeof(int));
    numElems = TacsUniqueSort(numElems, elemNums);
  }
}

/*
  Add the elements in the list to the domain. If the existing array is
  not long enough, create a new one that is large enough for
  everything. This creates an array that is exactly large enough - if
  many elements are going to be added, you should calculate in advance
  how many and allocate enough to start with.

  input:
  elemNums: the element numbers to add to the domain
  numElems: the number of elements to add
*/
void TACSFunction::addDomain(int _numElems, const int _elemNums[]) {
  if (funcDomain == NO_DOMAIN) {
    fprintf(stderr, "Cannot add function domain for %s\n", getObjectName());
    return;
  } else {
    funcDomain = SUB_DOMAIN;

    if (_numElems + numElems > maxElems) {
      maxElems = _numElems + numElems;
      int *temp = new int[maxElems];

      int i = 0;
      for (; i < numElems; i++) {
        temp[i] = elemNums[i];
      }
      for (; i < maxElems; i++) {
        temp[i] = _elemNums[i - numElems];
      }
      numElems = i;

      if (elemNums) {
        delete[] elemNums;
      }

      elemNums = temp;
    } else {
      for (int i = 0; i < _numElems; i++, numElems++) {
        elemNums[numElems] = _elemNums[i];
      }
    }

    numElems = TacsUniqueSort(numElems, elemNums);
  }
}

/*
  Get the elements in the domain of this object
*/
int TACSFunction::getElementNums(const int **_elemNums) {
  if (_elemNums) {
    *_elemNums = elemNums;
  }
  return numElems;
}

/*
  Get the TACSAssembler object associated with this function
*/
TACSAssembler *TACSFunction::getAssembler() { return assembler; }

/*
  Retrieve the object name
*/
const char *TACSFunction::getObjectName() { return "TACSFunction"; }
