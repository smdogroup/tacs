#ifndef TACS_FH5_H
#define TACS_FH5_H

/*
  Create an FH5 file from the TACSAssembler object

  Copyright (c) 2012 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/  

#include "FH5.h"
#include "TACSAssembler.h"

class TACSToFH5 : public TACSObject {
 public:
  TACSToFH5( TACSAssembler * _tacs,
             enum ElementType _elem_type,
             unsigned int _out_type );
  ~TACSToFH5();

  // Set the group name for each zone
  // --------------------------------
  void setComponentName( int comp_num, const char * group_name );
  
  // Write the data to a file
  // ------------------------
  void writeToFile( const char * filename );

 private:
  // Get a character string of the variable names
  char * getElementVarNames();

  TACSAssembler  * tacs;
  enum ElementType elem_type;
  unsigned int write_flag;
  int nvals;           // Number of total values per point
  int ndisplacements;  // Number of displacements
  int nstresses;       // Number of stresses and strains
  int nextras;         // Number of extra values
  int ncoordinates;    // Number of coordinates per point (always 9)

  int num_components; // The number of components in the model
  char ** component_names; // The names of each of the components
  char * variable_names; // The names of all the variables
};

#endif // TACS_FH5
