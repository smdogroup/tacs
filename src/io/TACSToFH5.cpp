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
*/
TACSToFH5::TACSToFH5(TACSAssembler *_assembler, ElementType _elem_type,
                     int _write_flag) {
  assembler = _assembler;
  assembler->incref();

  // Record the options
  elem_type = _elem_type;
  write_flag = _write_flag;

  // Form a flag that masks off the connectivity, nodes and displacements
  // which are handled separately from the remaining variables
  element_write_flag =
      write_flag & (~(TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                      TACS_OUTPUT_DISPLACEMENTS));

  // Count up the number of values that will be output for each point
  // in the mesh
  nvals = TacsGetTotalOutputCount(elem_type, element_write_flag);

  // Get a comma separated list of the variable names
  variable_names = getElementVarNames(element_write_flag);

  // Retrieve the number of components
  num_components = assembler->getNumComponents();

  // Allocate space for the component names
  component_names = new char *[num_components];
  memset(component_names, 0, num_components * sizeof(char *));

  for (int k = 0; k < num_components; k++) {
    char comp_name[128];
    sprintf(comp_name, "Component %d", k);
    setComponentName(k, comp_name);
  }
}

/**
   Free the FH5 object
*/
TACSToFH5::~TACSToFH5() {
  assembler->decref();

  // Deallocate the comma separated list of variable names
  delete[] variable_names;

  // Deallocate the component names
  for (int k = 0; k < num_components; k++) {
    if (component_names[k]) {
      delete[] component_names[k];
    }
  }
  delete[] component_names;
}

/**
   Set the specified component name for a group of elements

   @param comp_num The component number to set
   @param comp_name The component name to apply
*/
void TACSToFH5::setComponentName(int comp_num, const char *comp_name) {
  if (comp_num >= 0 && comp_num < num_components) {
    // If the name already exists, over-write it
    if (component_names[comp_num]) {
      delete[] component_names[comp_num];
    }

    // Allocate space for the name
    size_t len = strlen(comp_name) + 1;
    component_names[comp_num] = new char[len];

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
int TACSToFH5::writeToFile(const char *filename) {
  int rank, size;
  MPI_Comm_rank(assembler->getMPIComm(), &rank);
  MPI_Comm_size(assembler->getMPIComm(), &size);

  // Create the FH5 file object for writting
  TACSFH5File *file = new TACSFH5File(assembler->getMPIComm());
  file->incref();

  // Open the file - if possible for writing
  int fail = file->createFile(filename, num_components, component_names);

  if (fail) {
    file->decref();
    if (rank == 0) {
      fprintf(stderr, "[%d] TACSToFH5 error: Could not create file\n", rank);
    }
    return 1;
  }

  if (write_flag & TACS_OUTPUT_CONNECTIVITY) {
    writeConnectivity(file);
  }

  // Write out the nodes and solution vector to a file (continuous)
  if (write_flag & TACS_OUTPUT_NODES ||
      write_flag & TACS_OUTPUT_DISPLACEMENTS) {
    int vars_per_node = assembler->getVarsPerNode();

    // Find the maximum string length
    int str_len = strlen("X,Y,Z") + 1;
    int nd = TacsGetOutputComponentCount(elem_type, TACS_OUTPUT_DISPLACEMENTS);
    int k = 0;
    for (; (k < nd && k < vars_per_node); k++) {
      const char *stemp =
          TacsGetOutputComponentName(elem_type, TACS_OUTPUT_DISPLACEMENTS, k);
      str_len += strlen(stemp) + 1;
    }
    for (; k < vars_per_node; k++) {
      char stemp[64];
      sprintf(stemp, "v%d", k);
      str_len += strlen(stemp) + 1;
    }

    char *var_names = new char[str_len];
    var_names[0] = '\0';
    if (write_flag & TACS_OUTPUT_NODES) {
      sprintf(var_names, "X,Y,Z");
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      str_len = strlen(var_names);
      nd = TacsGetOutputComponentCount(elem_type, TACS_OUTPUT_DISPLACEMENTS);
      k = 0;
      for (; (k < nd && k < vars_per_node); k++) {
        const char *stemp =
            TacsGetOutputComponentName(elem_type, TACS_OUTPUT_DISPLACEMENTS, k);
        size_t len = strlen(var_names);
        if (k == 0 && !(write_flag & TACS_OUTPUT_NODES)) {
          sprintf(&(var_names[len]), "%s", stemp);
        } else {
          sprintf(&(var_names[len]), ",%s", stemp);
        }
      }
      for (; k < vars_per_node; k++) {
        size_t len = strlen(var_names);
        sprintf(&(var_names[len]), ",v%d", k);
      }
    }

    // Get the data from tacs
    TACSBVec *q, *X;
    assembler->getVariables(&q);
    assembler->getNodes(&X);

    // Get the arrays
    TacsScalar *ans_array, *Xarray;
    q->getArray(&ans_array);
    X->getArray(&Xarray);

    // Compute the first length of the array
    int nnodes = assembler->getNumOwnedNodes();
    int ndep = assembler->getNumDependentNodes();
    int dim1 = nnodes + ndep;

    // Check the dimension of the data
    int dim2 = 0, offset = 0;
    if (write_flag & TACS_OUTPUT_NODES) {
      dim2 = 3;
      offset = 3;
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      dim2 += assembler->getVarsPerNode();
    }

    // Allocate the float data
    float *float_data = new float[dim1 * dim2];

    // Write out the file
    if (write_flag & TACS_OUTPUT_NODES) {
      for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < 3; j++) {
          float_data[dim2 * i + j] = TacsRealPart(Xarray[3 * i + j]);
        }
      }
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < vars_per_node; j++) {
          float_data[dim2 * i + offset + j] =
              TacsRealPart(ans_array[vars_per_node * i + j]);
        }
      }
    }

    // Get the dependent parts of the arrays
    q->getDepArray(&ans_array);
    X->getDepArray(&Xarray);

    // Write the dependent nodes
    if (write_flag & TACS_OUTPUT_NODES) {
      for (int i = 0; i < ndep; i++) {
        for (int j = 0; j < 3; j++) {
          float_data[dim2 * (i + nnodes) + j] = TacsRealPart(Xarray[3 * i + j]);
        }
      }
    }
    if (write_flag & TACS_OUTPUT_DISPLACEMENTS) {
      for (int i = 0; i < ndep; i++) {
        for (int j = 0; j < vars_per_node; j++) {
          float_data[dim2 * (i + nnodes) + offset + j] =
              TacsRealPart(ans_array[vars_per_node * i + j]);
        }
      }
    }

    // Write the data with a time stamp from the simulation in TACS
    char data_name[128];
    double t = assembler->getSimulationTime();
    sprintf(data_name, "continuous data t=%.10e", t);
    file->writeZoneData(data_name, var_names, TACSFH5File::FH5_FLOAT, dim1,
                        dim2, float_data);
    delete[] float_data;
    delete[] var_names;
  }

  if (nvals > 0) {
    // Write out the data to a file
    TacsScalar *data;
    int dim1, dim2;
    assembler->getElementOutputData(elem_type, element_write_flag, &dim1, &dim2,
                                    &data);

    // Convert the data to float
    float *float_data = new float[dim1 * dim2];
    for (int i = 0; i < dim1 * dim2; i++) {
      float_data[i] = TacsRealPart(data[i]);
    }
    delete[] data;

    // Write the data with a time stamp from the simulation in TACS
    char data_name[128];
    double t = assembler->getSimulationTime();
    sprintf(data_name, "element data t=%.10e", t);
    file->writeZoneData(data_name, variable_names, TACSFH5File::FH5_FLOAT, dim1,
                        dim2, float_data);
    delete[] float_data;
  }

  file->close();
  file->decref();

  return 0;
}

/**
  Write out the connectivity information to the file
*/
int TACSToFH5::writeConnectivity(TACSFH5File *file) {
  int mpi_rank, mpi_size;
  MPI_Comm comm = assembler->getMPIComm();
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  // Record the layout types and component numbers
  int num_elements = assembler->getNumElements();
  int *comp_nums = new int[num_elements];
  int *layout_types = new int[num_elements];

  // Get the array of elements
  TACSElement **elements = assembler->getElements();

  // Set the layout types and the component numbers
  for (int i = 0; i < num_elements; i++) {
    layout_types[i] = elements[i]->getLayoutType();
    comp_nums[i] = elements[i]->getComponentNum();
  }

  // Write the component numbers to a zone
  int dim1 = num_elements;
  int dim2 = 1;
  char comp_name[] = "components";
  file->writeZoneData(comp_name, comp_name, TACSFH5File::FH5_INT, dim1, dim2,
                      comp_nums);
  delete[] comp_nums;

  // Write the layout types to a new zone
  char layout_name[] = "ltypes";
  file->writeZoneData(layout_name, layout_name, TACSFH5File::FH5_INT, dim1,
                      dim2, layout_types);
  delete[] layout_types;

  // Copy over the connectivity
  const int *ptr, *conn;
  assembler->getElementConnectivity(&ptr, &conn);

  // Modify the pointer so that it is consistent across processors
  int *ptr_copy = new int[num_elements + 1];
  memcpy(ptr_copy, ptr, (num_elements + 1) * sizeof(int));

  int offset = 0;
  if (mpi_rank > 0) {
    MPI_Status status;
    MPI_Recv(&offset, 1, MPI_INT, mpi_rank - 1, 1, comm, &status);
  }

  if (offset > 0) {
    for (int i = 0; i <= num_elements; i++) {
      ptr_copy[i] += offset;
    }
  }
  offset = ptr_copy[num_elements];

  if (mpi_rank < mpi_size - 1) {
    MPI_Send(&offset, 1, MPI_INT, mpi_rank + 1, 1, comm);
  }

  dim1 = num_elements;
  if (mpi_rank == mpi_size - 1) {
    dim1++;
  }

  char ptr_name[] = "ptr";
  file->writeZoneData(ptr_name, ptr_name, TACSFH5File::FH5_INT, dim1, dim2,
                      ptr_copy);
  delete[] ptr_copy;

  // Get the ownership range for each group of nodes
  const int *ownerRange;
  TACSNodeMap *nodeMap = assembler->getNodeMap();
  nodeMap->getOwnerRange(&ownerRange);

  // Get the number of nodes and dependent nodes
  int nnodes = assembler->getNumOwnedNodes();
  int ndep = assembler->getNumDependentNodes();

  // Count up the number of new nodes (including dependent nodes) for each
  // processor. Reset the connectivity so that both dependent and indepdent
  // nodes are all positive.
  int node_count = nnodes + ndep;
  int *new_owner_range = new int[mpi_size + 1];
  MPI_Allgather(&node_count, 1, MPI_INT, &new_owner_range[1], 1, MPI_INT, comm);
  new_owner_range[0] = 0;
  for (int k = 0; k < mpi_size; k++) {
    new_owner_range[k + 1] += new_owner_range[k];
  }

  // Create the copy of the connectivity
  int conn_size = ptr[num_elements];
  int *conn_copy = new int[conn_size];

  for (int i = 0; i < conn_size; i++) {
    // Get the global node number
    if (conn[i] < 0) {
      int dep = -conn[i] - 1;
      conn_copy[i] = new_owner_range[mpi_rank] + nnodes + dep;
    } else {
      if (conn[i] >= ownerRange[mpi_rank] &&
          conn[i] < ownerRange[mpi_rank + 1]) {
        conn_copy[i] =
            (conn[i] - ownerRange[mpi_rank] + new_owner_range[mpi_rank]);
      } else {
        for (int j = 0; j < mpi_size; j++) {
          if (conn[i] >= ownerRange[j] && conn[i] < ownerRange[j + 1]) {
            conn_copy[i] = (conn[i] - ownerRange[j] + new_owner_range[j]);
            break;
          }
        }
      }
    }
  }

  dim1 = conn_size;
  dim2 = 1;
  char conn_name[] = "connectivity";
  file->writeZoneData(conn_name, conn_name, TACSFH5File::FH5_INT, dim1, dim2,
                      conn_copy);
  delete[] conn_copy;
  delete[] new_owner_range;

  return 0;
}

/**
  Create a comma-separated list of the element variable names
*/
char *TACSToFH5::getElementVarNames(int flag) {
  // Find the first variable name
  char *elem_vars = NULL;
  char *output_names[3] = {NULL, NULL, NULL};

  int out_types[3] = {TACS_OUTPUT_STRAINS, TACS_OUTPUT_STRESSES,
                      TACS_OUTPUT_EXTRAS};

  for (int k = 0; k < 3; k++) {
    if (flag & out_types[k]) {
      const char *stemp = NULL;
      int nd = TacsGetOutputComponentCount(elem_type, out_types[k]);
      size_t str_len = 1;
      for (int i = 0; i < nd; i++) {
        stemp = TacsGetOutputComponentName(elem_type, out_types[k], i);
        str_len += strlen(stemp) + 1;
      }
      char *temp = new char[str_len];
      if (nd > 0) {
        stemp = TacsGetOutputComponentName(elem_type, out_types[k], 0);
        strcpy(temp, stemp);
        for (int i = 1; i < nd; i++) {
          stemp = TacsGetOutputComponentName(elem_type, out_types[k], i);
          size_t len = strlen(temp);
          sprintf(&(temp[len]), ",%s", stemp);
        }
      }
      output_names[k] = temp;
    }
  }

  // Count up the size of the elem_vars string
  int elem_size = 1;  // Extra space for either a comma or \0
  for (int k = 0; k < 3; k++) {
    if (output_names[k]) {
      elem_size += strlen(output_names[k]) + 1;
    }
  }

  elem_vars = new char[elem_size];

  // Copy the first zone into the list directly
  int k = 0;
  for (; k < 3; k++) {
    if (output_names[k]) {
      strcpy(elem_vars, output_names[k]);
      k++;
      break;
    }
  }

  // For subsequent non-zero zones - add a comma before adding
  // the remainder of the list
  for (; k < 3; k++) {
    if (output_names[k]) {
      int len = strlen(elem_vars);
      sprintf(&elem_vars[len], ",%s", output_names[k]);
    }
  }

  for (int k = 0; k < 3; k++) {
    if (output_names[k]) {
      delete[] output_names[k];
    }
  }

  return elem_vars;
}
