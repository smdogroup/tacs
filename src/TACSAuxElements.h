#ifndef TACS_AUX_ELEMENTS_H
#define TACS_AUX_ELEMENTS_H

/*
  TACSAuxElements implementation

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

#include "TACSObject.h"
#include "TACSElement.h"

/*
  The following is th
*/
class TACSAuxElem {
 public:
  TACSElement *elem;
  int num;
};

/*
  The TACSAuxiliaryElements class 

  This class is a simple way to add tractions/body forces into 
  TACS 

  a class for adding extra
  forces/elements/contributions to the Jacobian.

  The following class defines a group of extra elements that can be
  added to the TACSAssembler object. These additional elements can
  only contribute existing the non-zero
  
*/
class TACSAuxElements : public TACSOptObject {
 public:
  TACSAuxElements( int _num_elems );
  ~TACSAuxElements();

  // Sort the list of elements
  // -------------------------
  void sort();

  // Add a single element
  // --------------------
  void addElement( int num, TACSElement *elem );

  // Add an entire list of elements
  // ------------------------------
  void addElements( int nums[], TACSElement **elem, int num_elems );

  // Get the elements and sort them (if they are not already)
  // --------------------------------------------------------
  int getAuxElements( TACSAuxElem **elems );

  // Functions to control the design variables
  // -----------------------------------------
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[], int numDVs );
  
  // Print the name of the TACSObject
  const char* TACSObjectName(){ return auxName; }

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

#endif // TACS_AUXILIARY_ELEMENTS_H
