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

#ifndef TACS_TO_FH5_H
#define TACS_TO_FH5_H

/*
  Create an FH5 file from the TACSAssembler object
*/  

#include "FH5.h"
#include "TACSAssembler.h"

class TACSToFH5 : public TACSObject {
 public:
  TACSToFH5( TACSAssembler *_tacs,
             enum ElementType _elem_type,
             unsigned int _out_type );
  ~TACSToFH5();

  // Set the group name for each zone
  // --------------------------------
  void setComponentName( int comp_num, const char *group_name );
  
  // Write the data to a file
  // ------------------------
  void writeToFile( const char *filename );

 private:
  // Get a character string of the variable names
  char *getElementVarNames();

  TACSAssembler *tacs;
  enum ElementType elem_type;
  unsigned int write_flag;
  int nvals;           // Number of total values per point
  int ndisplacements;  // Number of displacements
  int nstresses;       // Number of stresses and strains
  int nextras;         // Number of extra values
  int ncoordinates;    // Number of coordinates per point (always 9)

  int num_components; // The number of components in the model
  char **component_names; // The names of each of the components
  char *variable_names; // The names of all the variables
};

#endif // TACS_TO_FH5
