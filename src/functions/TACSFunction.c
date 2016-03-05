#include "TACSFunction.h"

/*
  Base TACSFunction implementation

  Copyright (c) 2010 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  Create the function base class with the default domain - typically
  the entire finite-element mesh.

  input:
  tacs:          the TACSAssembler object
  funcDomain:    the domain type
  maxElems:      the maximum number of elements expected
  numIterations: iterations required to evaluate the function
*/
TACSFunction::TACSFunction( TACSAssembler * _tacs, FunctionDomain _funcDomain, 
			    int _maxElems, int _numIterations ): 
numIterations(_numIterations){
  tacs = _tacs;
  tacs->incref();

  initFlag = 0;
  funcDomain = _funcDomain;
  
  maxElems = (_maxElems > 0 ? _maxElems : 0);
  numElems = 0;
  elemNums = NULL;
}

/*
  Create the function base class with the specified element domain.
  The element domain is indicated by the element numbers.

  input:
  tacs:          the TACSAssembler object
  elemNums:      the element numbers
  numElems:      the number of elements in elemNums
  maxElems:      the maximum number of elements expected
  numIterations: iterations required to evaluate the function
*/ 
TACSFunction::TACSFunction( TACSAssembler * _tacs, int _elemNums[], 
			    int _numElems, int _maxElems, 
			    int _numIterations ): 
numIterations(_numIterations){
  tacs = _tacs;
  tacs->incref();

  initFlag = 0;
  funcDomain = SUB_DOMAIN;

  maxElems = (_maxElems > _numElems ? _maxElems : _numElems);
  numElems = _numElems;
  elemNums = new int[ maxElems ];
  
  for ( int i = 0; i < numElems; i++ ){
    elemNums[i] = _elemNums[i];
  }

  numElems = FElibrary::uniqueSort(elemNums, numElems);
}

/*
  Destroy the TACSFunction base class
*/
TACSFunction::~TACSFunction(){
  if (elemNums){
    delete [] elemNums;
  }
  tacs->decref();
}

/*
  Get the sizes of the work arrays used for the function evaluations
  The arrays are allocated externally in order to be thread-safe.  
*/
void TACSFunction::getEvalWorkSizes( int * iwork, int * work ){
  *iwork = 0, *work = 0;
}

int TACSFunction::getSVSensWorkSize(){
  return 0;
}

int TACSFunction::getDVSensWorkSize(){
  return 0;
}

int TACSFunction::getXptSensWorkSize(){
  return 0;
}

/*
  Extend the size of the domain to maxElems.

  input:
  maxElems:  the new maximum size of the domain
*/
void TACSFunction::setDomainSize( int _maxElems ){
  if (funcDomain == NO_DOMAIN){
    fprintf(stderr, "Cannot set function domain for %s\n",
	    this->functionName());
    return;
  }
  else if (maxElems < _maxElems){
    funcDomain = SUB_DOMAIN;
    maxElems = _maxElems;

    int * temp = new int[ maxElems ];
      
    int i = 0;
    for ( ; i < numElems; i++ ){
      temp[i] = elemNums[i];
    }

    if ( elemNums ){
      delete [] elemNums;
    }

    elemNums = temp;
  }
}

/*
  Extend the size of the domain to tacs->getNumElements().
  All that should ever be required!
*/
void TACSFunction::setMaxDomainSize(){
  setDomainSize(tacs->getNumElements());
}

/*
  Overwrite the domain in the function with a new set of elements.
  This reallocates the existing array if it is not long enough.  
 
  input:
  elemNums: the element numbers used to set the domain
  numElems: the number of elements to add
*/
void TACSFunction::setDomain( int _elemNums[], int _numElems ){
  if (funcDomain == NO_DOMAIN){
    fprintf(stderr, "Cannot set function domain for %s\n",
	    this->functionName());
    return;
  }
  else {
    funcDomain = SUB_DOMAIN;

    if (_numElems > maxElems){
      if (elemNums){
	delete [] elemNums;
      }
      numElems = _numElems;
      maxElems = _numElems;
      elemNums = new int[ maxElems ];
    }
    else {
      numElems = _numElems;
    }

    memcpy(elemNums, _elemNums, numElems*sizeof(int));
    numElems = FElibrary::uniqueSort(elemNums, numElems);
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
void TACSFunction::addDomain( int _elemNums[], int _numElems ){
  if (funcDomain == NO_DOMAIN){
    fprintf(stderr, "Cannot add function domain for %s\n",
	    this->functionName());
    return;
  }
  else {
    funcDomain = SUB_DOMAIN;

    if (_numElems + numElems > maxElems){
      maxElems = _numElems + numElems;
      int * temp = new int[ maxElems ];
      
      int i = 0;
      for ( ; i < numElems; i++ ){
	temp[i] = elemNums[i];
      }

      for ( ; i < maxElems; i++ ){
	temp[i] = _elemNums[i-numElems];
      }

      if ( elemNums ){
	delete [] elemNums;
      }

      elemNums = temp;
    }
    else {
      for ( int i = 0; i < _numElems; i++, numElems++ ){
	elemNums[numElems] = _elemNums[i];
      }      
    }

    numElems = FElibrary::uniqueSort(elemNums, numElems);
  }
}

/*
  Retrieve the name of the function object
*/
const char * TACSFunction::TACSObjectName(){ return this->functionName(); }

/*
  Retrieve the type of domain specified by this object
*/
enum TACSFunction::FunctionDomain TACSFunction::getDomain(){
  return funcDomain; 
}

/*
  Retrieve the number of elements in the domain of this object
*/
int TACSFunction::getNumElements(){
  return numElems;
}

/*
  Get the elements in the domain of this object
*/
int TACSFunction::getElements( const int ** _elemNums ){
  *_elemNums = elemNums;
  return numElems;
}

TACSAssembler * TACSFunction::getTACS(){ return tacs; }
