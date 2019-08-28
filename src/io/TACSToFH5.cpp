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

#include "TACSToFH5.h"

/**
   Create the TACSToFH5 object.

   This class creates an .f5 file in parallel using a binary output
   format.  This format can then be converted into output for viewing
   or used for post-processing in some cases.

   The object creates two types of output: (1) element-wise
   independent output or (2) continuous output with nodal averaging of
   quantities of interest. (1) may be useful if you want to look at
   the stress or strain independently in each element, (2) can produce
   a smaller output file and may be best for large-scale computations.

   @param assembler The TACSAssembler object
   @param elem_type The type of element output to generate
   @param write_flag XOR flag indicating classes of output to write
   @param average_node_data Average the node data or use element-independent
*/
TACSToFH5::TACSToFH5( TACSAssembler *_assembler,
                      ElementType _elem_type,
                      int _write_flag,
                      int _average_node_data ){
  assembler = _assembler;
  assembler->incref();

  // Record the options
  elem_type = _elem_type;
  write_flag = _write_flag;
  average_node_data = _average_node_data;

  // Count up the number of values that will be output for each point
  // in the mesh
  nvals = TacsGetTotalOutputCount(elem_type, write_flag);

  // Get a comma separated list of the variable names
  variable_names = getElementVarNames();

  // Retrieve the number of components
  num_components = assembler->getNumComponents();

  // Allocate space for the component names
  component_names = new char*[ num_components ];
  memset(component_names, 0, num_components*sizeof(char*));

  for ( int k = 0; k < num_components; k++ ){
    char comp_name[128];
    sprintf(comp_name, "Component %d", k);
    setComponentName(k, comp_name);
  }
}

/**
   Free the FH5 object
*/
TACSToFH5::~TACSToFH5(){
  assembler->decref();

  // Deallocate the comma separated list of variable names
  delete [] variable_names;

  // Deallocate the component names
  for ( int k = 0; k < num_components; k++ ){
    if (component_names[k]){
      delete [] component_names[k];
    }
  }
  delete [] component_names;
}

/**
   Set the specified component name for a group of elements

   @param comp_num The component number to set
   @param comp_name The component name to apply
*/
void TACSToFH5::setComponentName( int comp_num, const char *comp_name ){
  if (comp_num >= 0 && comp_num < num_components){
    // If the name already exists, over-write it
    if (component_names[comp_num]){
      delete [] component_names[comp_num];
    }

    // Allocate space for the name
    size_t len = strlen(comp_name)+1;
    component_names[comp_num] = new char[ len ];

    // Copy the component name
    strcpy(component_names[comp_num], comp_name);
  }
}

/**
   Write the data stored in the TACSAssembler object to a file

   Write the data in parallel to a binary file format with the
   specified name. This code writes out the component names, component
   numbers, layout types, element data and global connectivity (when
   node-averaging data is specified). Computing and writting all this
   data can be expensive and take considerable diskspace, so care
   should be exercised when creating these files.

   @param filename The name of the file to write
*/
void TACSToFH5::writeToFile( const char *filename ){
  int rank, size;
  MPI_Comm_rank(assembler->getMPIComm(), &rank);
  MPI_Comm_size(assembler->getMPIComm(), &size);

  // Create the FH5 file object for writting
  FH5File *file = new FH5File(assembler->getMPIComm());
  file->incref();

  // Open the file - if possible for writing
  int write_err = file->createFile(filename, component_names,
                                   num_components);

  if (write_err){
    file->decref();
    if (rank == 0){
      fprintf(stderr, "[%d] TACSToFH5 error: Could not create file\n",
              rank);
    }
    return;
  }

  // Record the layout types and component numbers
  int num_elements = assembler->getNumElements();
  int *comp_nums = new int[ num_elements ];
  int *layout_types = new int[ num_elements ];

  // Get the array of elements
  TACSElement **elements = assembler->getElements();

  // Set the layout types and the component numbers
  for ( int i = 0; i < num_elements; i++ ){
    layout_types[i] = elements[i]->getLayoutType();
    comp_nums[i] = elements[i]->getComponentNum();
  }

  // Write the component numbers to a zone
  int dim1 = num_elements;
  int dim2 = 1;
  char comp_name[] = "components";
  file->writeZoneData(comp_name, comp_name, FH5File::FH5_INT,
                      comp_nums, dim1, dim2);
  delete [] comp_nums;

  // Write the layout types to a new zone
  char layout_name[] = "ltypes";
  file->writeZoneData(layout_name, layout_name, FH5File::FH5_INT,
                      layout_types, dim1, dim2);
  delete [] layout_types;

  if (average_node_data){
    // Write out the data to a file
    TacsScalar *data;
    assembler->getElementOutputData(elem_type, write_flag,
                                    &dim1, &dim2, &data);

    // Convert the data to float
    float *float_data = new float[ dim1*dim2 ];
    for ( int i = 0; i < dim1*dim2; i++ ){
      float_data[i] = data[i];
    }
    delete [] data;

    // Write the data with a time stamp from the simulation in TACS
    char data_name[128];
    double t = assembler->getSimulationTime();
    sprintf(data_name, "data t=%.10e", t);
    file->writeZoneData(data_name, variable_names,
                        FH5File::FH5_FLOAT, float_data, dim1, dim2);
    delete [] float_data;
  }
  else {
    // Create the continuous output data
    TACSBVec *data_vec =
      assembler->getNodeAverageOutputData(elem_type, write_flag);
    data_vec->incref();

    /*
    float *float_data = NULL;

    // Convert the data to float
    float_data = new float[ dim1*dim2 ];
    for ( int i = 0; i < dim1*dim2; i++ ){
      float_data[i] = data[i];
    }
    delete [] data;

    // Write the data with a time stamp from the simulation in TACS
    char data_name[128];
    double t = assembler->getSimulationTime();
    sprintf(data_name, "data t=%.10e", t);
    file->writeZoneData(data_name, variable_names,
                        FH5File::FH5_FLOAT, float_data, dim1, dim2);
    delete [] float_data;
    */

    data_vec->decref();
  }

  file->close();
  file->decref();
}

/**
  Create a comma-separated list of the element variable names
*/
char *TACSToFH5::getElementVarNames(){
  // Find the first variable name
  char *elem_vars = NULL;
  char *output_names[5] = { NULL, NULL, NULL, NULL, NULL };

  int out_types[5] =
    { TACS_OUTPUT_NODES,
      TACS_OUTPUT_DISPLACEMENTS,
      TACS_OUTPUT_STRAINS,
      TACS_OUTPUT_STRESSES,
      TACS_OUTPUT_EXTRAS };

  for ( int k = 0; k < 5; k++ ){
    if (write_flag & out_types[k]){
      const char *stemp = NULL;
      int nd = TacsGetOutputComponentCount(elem_type, out_types[k]);
      size_t str_len = 1;
      for ( int i = 0; i < nd; i++ ){
        stemp = TacsGetOutputComponentName(elem_type, out_types[k], i);
        str_len += strlen(stemp) + 1;
      }
      char *temp = new char[ str_len ];
      if (nd > 0){
        stemp = TacsGetOutputComponentName(elem_type, out_types[k], 0);
        strcpy(temp, stemp);
        for ( int i = 1; i < nd; i++ ){
          stemp = TacsGetOutputComponentName(elem_type, out_types[k], i);
          size_t len = strlen(temp);
          sprintf(&(temp[len]), ",%s", stemp);
        }
      }
      output_names[k] = temp;
    }
  }

  // Count up the size of the elem_vars string
  int elem_size = 1; // Extra space for either a comma or \0
  for ( int k = 0; k < 5; k++ ){
    if (output_names[k]){
      elem_size += strlen(output_names[k])+1;
    }
  }

  elem_vars = new char[ elem_size ];

  // Copy the first zone into the list directly
  int k = 0;
  for ( ; k < 5; k++ ){
    if (output_names[k]){
      strcpy(elem_vars, output_names[k]);
      k++;
      break;
    }
  }

  // For subsequent non-zero zones - add a comma before adding
  // the remainder of the list
  for ( ; k < 5; k++ ){
    if (output_names[k]){
      int len = strlen(elem_vars);
      sprintf(&elem_vars[len], ",%s", output_names[k]);
    }
  }

  for ( int k = 0; k < 5; k++ ){
    if (output_names[k]){ delete [] output_names[k]; }
  }

  return elem_vars;
}
