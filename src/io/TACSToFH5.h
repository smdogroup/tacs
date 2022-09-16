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

#include "TACSAssembler.h"
#include "TACSFH5.h"

/**
  Write out the solution data to a binary file in a TACS-specific
  format.

  This .f5 file format is specific to TACS, but is written in
  parallel.  Data recorded in the file can later be accessed and
  converted to formats for visualization.
*/
class TACSToFH5 : public TACSObject {
 public:
  TACSToFH5(TACSAssembler *assembler, ElementType elem_type, int write_flag);
  ~TACSToFH5();

  // Set the group name for each zone
  void setComponentName(int comp_num, const char *group_name);

  // Write the data to a file
  int writeToFile(const char *filename);

 private:
  // Get a character string of the variable names
  char *getElementVarNames(int flag);

  // Write the connectivity information to a file
  int writeConnectivity(TACSFH5File *file);

  // The Assembler object
  TACSAssembler *assembler;

  // Parameters to control how data is written to the file
  ElementType elem_type;   // Write flag type
  int write_flag;          // Keep track of which data to write
  int element_write_flag;  // Element-wise data write flag
  int nvals;               // The number of element-wise values

  int num_components;      // The number of components in the model
  char **component_names;  // The names of each of the components
  char *variable_names;    // The names of all the variables
};

#endif  // TACS_TO_FH5
