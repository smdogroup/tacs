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

#include "TACSMeshLoader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*!
  This is an interface for reading NASTRAN-style files.

  Description of the Nastran input file:

  It is a fixed-width format. There are five sections in a Nastran
  file - often called an input deck.

  1. Nastran statement
  2. File management section
  3. Executive control section
  CEND
  4. Case control section - load cases defined
  BEGIN BULK
  5. Bulk data section
  ENDDATA

  The first two sections are optional. The delimiters are not
  optional.

  Format of bulk data:
  - Each line must have 80 columns.
  - A bulk data entry may span multiple lines
  - Three types of data formats: integer, real, character string
  - Each field in the file has a specific input format - the quick
  reference guide specifies the format for each type of element/data
  type.

  - There are three types of field formats: small field, large field,
  and free field format

  Small field format:
  - Each line is divided into 10 fields, 8 columns wide

  Large field format:
  - Each line is divided into 6 fields, the first is 8 columns wide,
  the next 4 are 16 wide and the last is 8 columns wide. The large field
  format is signified by an asterisk after the keyword.

  Free field format:
  - Fields separated by commas, only precise to eight digits - higher
  precision truncated.

  Continuation entries:

  Method 1: +M101, +M101 on the preceeding and continuation entries.
  Method 2: Last entry of preceeding and first entry of the continuation
  line left blank.

  - Input data in fields 1, 10 must be left justified
  - All real numbers must have a decimal place - including zero
  - Default values may be left blank
  - Comments start with a dollar sign
*/

/*
  Functions for sorting a list such that:

  arg_sort_list[list[i]] is in ascending order
*/
static const int *arg_sort_list = NULL;

static int compare_arg_sort(const void *a, const void *b) {
  return arg_sort_list[*(int *)a] - arg_sort_list[*(int *)b];
}

/*
  Read a line from the buffer.

  Return the number of read characters. Do not exceed the buffer
  length.

  Given the line buffer 'line', and the size of the line buffer line_len
*/
static int read_buffer_line(char *line, size_t line_len, size_t *loc,
                            char *buffer, size_t buffer_len) {
  size_t i = 0;
  for (; (i < line_len) && (*loc < buffer_len); i++, (*loc)++) {
    if (buffer[*loc] == '\n') {
      break;
    }
    line[i] = buffer[*loc];
  }
  if (i < line_len) {
    line[i] = '\0';
  }

  // Read until the end of the line
  while ((*loc < buffer_len) && (buffer[*loc] != '\n')) {
    (*loc)++;
  }

  (*loc)++;  // Increment so that loc is at the next location

  return i;
}

/*
  Reverse look up.

  Given the sorted array list[arg[k]], find k, such that

  list[arg[k]] = var
*/
static int find_index_arg_sorted(int var, int size, const int *list,
                                 const int *args) {
  // Binary search an array to find k such that list[k] = var,
  // where the array list[args[k]] is sorted in ascending
  // order
  int high = size - 1;
  int low = 0;
  int high_val = list[args[high]];
  int low_val = list[args[low]];

  // Check if the index is at the end points
  if (var == low_val) {
    return low;
  } else if (var < low_val) {
    return -1;
  }

  if (var == high_val) {
    return high;
  } else if (var > high_val) {
    return -1;
  }

  int mid = low + (int)((high - low) / 2);

  // While there are values left in the list
  while (low != mid) {
    int mid_val = list[args[mid]];
    if (mid_val == var) {
      return mid;
    }

    if (var < mid_val) {
      high = mid;
      high_val = mid_val;
    } else {
      low = mid;
      low_val = mid_val;
    }

    mid = low + (int)((high - low) / 2);
  }

  return -1;
}

/*
  Convert a Nastran-style number with an exponent to a double.
*/
static double bdf_atof(char *str) {
  // First, check if the string contains an E/e or D/d - if so, convert it
  int slen = strlen(str);
  for (int i = 0; i < slen; i++) {
    if (str[i] == 'e' || str[i] == 'E') {
      return atof(str);
    } else if (str[i] == 'd' || str[i] == 'D') {
      str[i] = 'e';
      return atof(str);
    }
  }

  // Convert the special Nastran number format without e/E or d/D
  // 4.-4 or 5.34+2 etc.
  char temp[24];
  int i = 0, j = 0;
  while (i < slen && str[i] == ' ') {
    i++;
  }
  if (i == slen) {
    return 0.0;
  }

  // Take care of a leading minus sign
  if (str[i] == '-') {
    temp[j] = str[i];
    j++, i++;
  }

  for (; i < slen; i++, j++) {
    // Add the character 'e' before the exponent
    if (str[i] == '-' || str[i] == '+') {
      temp[j] = 'e';
      j++;
    }
    temp[j] = str[i];
  }
  temp[j] = '\0';

  return atof(temp);
}

/*
  Parse the long-field format

  The long field format is split into the fixed widths:

  0 --- 8 ---24--- 40--- 56--- 72 ---80

  GRID* num  coord x     y
        z

  This code ignores the coordinate info
*/
static void parse_node_long_field(char *line, char *line2, int *node, double *x,
                                  double *y, double *z) {
  char Node[32], X[32], Y[32], Z[32];

  strncpy(Node, &line[8], 16);
  Node[16] = '\0';
  strncpy(X, &line[40], 16);
  X[16] = '\0';
  strncpy(Y, &line[56], 16);
  Y[16] = '\0';
  strncpy(Z, &line2[8], 16);
  Z[16] = '\0';

  *node = atoi(Node);
  *x = bdf_atof(X);
  *y = bdf_atof(Y);
  *z = bdf_atof(Z);
}

/*
  Parse the short-field or free-field comma separated format

  The short field format is fixed-width as follows:

  0 --- 8 --- 16 --- 24 --- 32 --- 40
  GRID  num   coord  x      y      z
*/
static void parse_node_short_free_field(char *line, int *node, double *x,
                                        double *y, double *z) {
  char field[5][32];

  // Look for a comma
  int len = strlen(line);
  int comma_format = 0;
  for (int i = 0; (i < len) && (i < 80); i++) {
    if (line[i] == ',') {
      comma_format = 1;
      break;
    }
  }

  if (comma_format) {
    int start = 8;
    int end = 8;
    for (int i = 0; i < 5; i++) {
      start = end;
      while (end < len && line[end] != ',') end++;
      int flen = end - start;
      strncpy(field[i], &line[start], flen);
      field[i][flen] = '\0';
      end++;
    }
  } else {  // Short-format, fixed width
    strncpy(field[0], &line[8], 8);
    field[0][8] = '\0';
    strncpy(field[2], &line[24], 8);
    field[2][8] = '\0';
    strncpy(field[3], &line[32], 8);
    field[3][8] = '\0';
    strncpy(field[4], &line[40], 8);
    field[4][8] = '\0';
  }

  *node = atoi(field[0]);
  *x = bdf_atof(field[2]);
  *y = bdf_atof(field[3]);
  *z = bdf_atof(field[4]);
}

/*
  Parse a single element entry in the file.

  The element entries are fixed-width. The first entry consists
  of an element type
*/
static int parse_element_field(size_t *loc, char *buffer,
                               const size_t buffer_len, const int entry_width,
                               const int max_num_nodes, int *elem_num,
                               int *component_num, int *node_nums,
                               int *num_nodes) {
  int fail = 0;
  *num_nodes = -1;
  char line[81];  // Space for the line read from the buffer
  char node[17];  // Space for the node number

  // Offset to find the first entry - always at position 8
  int entry = 8;

  size_t temp_loc = *loc;
  int line_len =
      read_buffer_line(line, sizeof(line), &temp_loc, buffer, buffer_len);

  // This is a hard failure - no element/component defined
  if (line_len <= 0) {
    fail = 1;
    return fail;
  }

  // Set the buffer location
  *loc = temp_loc;

  // Find the element number
  strncpy(node, &line[entry], entry_width);
  node[entry_width] = '\0';
  *elem_num = atoi(node);
  entry += entry_width;

  // The element indices must be positive
  if (*elem_num <= 0) {
    fail = 1;
    return fail;
  }

  // Find the component number
  strncpy(node, &line[entry], entry_width);
  node[entry_width] = '\0';
  *component_num = atoi(node);
  entry += entry_width;

  // The component indices must be positive
  if (*component_num <= 0) {
    fail = 1;
    return fail;
  }

  // Keep track of the number of nodes found
  int n = 0;

  while (n < max_num_nodes) {
    for (; (n < max_num_nodes) && (entry < line_len); entry += entry_width) {
      // Parse the fixed-width entry containing the node number
      strncpy(node, &line[entry], entry_width);
      node[entry_width] = '\0';

      // Check if the entry is entirely blank. Skip the update if it
      // is completely blank
      int count = 0;
      while (count < entry_width && isspace(node[count])) count++;
      if (count < entry_width) {
        int temp = atoi(node);
        if (temp > 0) {
          if (node_nums) {
            node_nums[n] = temp;
          }
          n++;
        } else if (temp <= 0) {
          *num_nodes = n;
          fail = 0;
          return fail;
        }
      }
    }

    if ((n < max_num_nodes) && (entry >= line_len)) {
      // Try to read the next line of the file
      temp_loc = *loc;
      line_len =
          read_buffer_line(line, sizeof(line), &temp_loc, buffer, buffer_len);

      // This is not a failure - element index/component defined, but
      // if either of these are true, then the next line does not
      // contain a continuation of this element
      if (line_len <= 0) {
        *num_nodes = n;
        fail = 0;
        return fail;
      } else if (!(line[0] == '*' || line[0] == ' ')) {
        // The first character of a continuation line must be "*" or " ".
        // Otherwise, this next line is a new element that must be
        // parsed separately
        *num_nodes = n;
        fail = 0;
        return fail;
      }

      // The next line is a continuation of this element.
      // Set the end of the new buffer accordingly
      *loc = temp_loc;

      // Set the new entry location for the next line. Note that
      // this is offset by 8 regardless of the entry width
      entry = 8;
    }
  }

  // This read did not fail
  fail = 0;
  *num_nodes = n;
  return fail;
}

/*
  The TACSMeshLoader class

  To load a mesh, you first pass in the communicator on which TACS
  will be allocated. You can then scan a file, set the element objects
  and create the TACSAssembler object.

  This constructor simply sets all data to NULL and stores the
  communicator for later use. Note that the file is only scanned on
  the root processor.
*/
TACSMeshLoader::TACSMeshLoader(MPI_Comm _comm) {
  comm = _comm;

  // Initialize everything to zero
  num_nodes = num_elements = 0;
  num_bcs = 0;
  file_node_nums = NULL;
  file_elem_nums = NULL;
  node_arg_sort_list = NULL;
  elem_arg_sort_list = NULL;

  elem_node_conn = elem_node_ptr = NULL;
  elem_component = NULL;
  Xpts = NULL;
  bc_vals = NULL;
  bc_nodes = bc_vars = bc_ptr = NULL;

  num_components = 0;
  elements = NULL;
  component_elems = NULL;
  component_descript = NULL;

  // Set the creator object to NULL
  creator = NULL;
}

/*
  Destroy existing data that may have been allocated

  Note that most of this data will only be allocated on the root
  processor
*/
TACSMeshLoader::~TACSMeshLoader() {
  if (elem_node_conn) {
    delete[] elem_node_conn;
  }
  if (elem_node_ptr) {
    delete[] elem_node_ptr;
  }
  if (elem_component) {
    delete[] elem_component;
  }
  if (Xpts) {
    delete[] Xpts;
  }
  if (bc_nodes) {
    delete[] bc_nodes;
  }
  if (bc_vars) {
    delete[] bc_vars;
  }
  if (bc_vals) {
    delete[] bc_vals;
  }
  if (bc_ptr) {
    delete[] bc_ptr;
  }

  if (elements) {
    for (int k = 0; k < num_components; k++) {
      if (elements[k]) {
        elements[k]->decref();
      }
    }
    delete[] elements;
  }
  if (component_elems) {
    delete[] component_elems;
  }
  if (component_descript) {
    delete[] component_descript;
  }
  if (file_node_nums) {
    delete[] file_node_nums;
  }
  if (file_elem_nums) {
    delete[] file_elem_nums;
  }
  if (node_arg_sort_list) {
    delete[] node_arg_sort_list;
  }
  if (elem_arg_sort_list) {
    delete[] elem_arg_sort_list;
  }

  // Free the creator object
  if (creator) {
    creator->decref();
  }
}

/*
  Get the number of components defined by the data
*/
int TACSMeshLoader::getNumComponents() { return num_components; }

/*
  Set the element associated with a given component number
*/
void TACSMeshLoader::setElement(int component_num, TACSElement *_element) {
  if (_element && (component_num >= 0) && (component_num < num_components)) {
    _element->incref();
    _element->setComponentNum(component_num);
    elements[component_num] = _element;
  }
}

/*
  Get the component description from the file
*/
const char *TACSMeshLoader::getComponentDescript(int comp_num) {
  if (component_descript && (comp_num >= 0) && (comp_num < num_components)) {
    return &component_descript[33 * comp_num];
  }
  return NULL;
}

/*
  Retrieve the element description corresponding to the component number
*/
const char *TACSMeshLoader::getElementDescript(int comp_num) {
  if (component_elems && (comp_num >= 0) && (comp_num < num_components)) {
    return &component_elems[9 * comp_num];
  }
  return NULL;  // No associated element
}

/*
  This functions scans a Nastran BDF file - only scanning in
  information from the bulk data section.

  Only the element types, boundary conditions, connectivitiy and GRID
  entries are scanned.  Any entries associated with constitutive
  properties are ignored.
*/
int TACSMeshLoader::scanBDFFile(const char *file_name) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  int fail = 0;

  const int root = 0;
  if (rank == root) {
    FILE *fp = fopen(file_name, "r");
    if (!fp) {
      fprintf(stderr, "TACSMeshLoader: Unable to open file %s\n", file_name);
      fail = 1;
      MPI_Abort(comm, fail);
      return fail;
    }

    // Count up the number of nodes, elements and size of connectivity data
    num_nodes = 0;
    num_elements = 0;
    num_components = 0;
    num_bcs = 0;

    // The size of the connectivity arrays
    int elem_conn_size = 0;
    int bc_vars_size = 0;

    // Each line can only be 80 characters long
    char line[81];

    // Determine the size of the file
    fseek(fp, 0, SEEK_END);
    size_t buffer_len = ftell(fp);
    rewind(fp);

    // Allocate enough space to store the entire file
    char *buffer = new char[buffer_len];
    if (fread(buffer, 1, buffer_len, fp) != buffer_len) {
      fprintf(stderr, "[%d] TACSMeshLoader: Problem reading file %s\n", rank,
              file_name);
      MPI_Abort(comm, 1);
      return 1;
    }
    fclose(fp);

    // Keep track of where the current point in the buffer is
    size_t buffer_loc = 0;
    read_buffer_line(line, sizeof(line), &buffer_loc, buffer, buffer_len);

    // Flags which indicate where the bulk data begins
    int in_bulk = 0;
    int bulk_start = 0;

    // Scan the file for the begin bulk location. If none exists, then
    // the whole file is treated as bulk data.
    while (buffer_loc < buffer_len) {
      if (strncmp(line, "BEGIN BULK", 10) == 0) {
        in_bulk = 1;
        bulk_start = buffer_loc;
      }
      read_buffer_line(line, sizeof(line), &buffer_loc, buffer, buffer_len);
    }

    // Treat the entire file as bulk data
    if (!in_bulk) {
      in_bulk = 1;
    }

    // Set the starting location
    buffer_loc = bulk_start;

    // Keep track of the maximum number of nodes defined by
    // any of the elements
    int max_element_conn = -1;

    while (buffer_loc < buffer_len) {
      // Read the first line of the buffer
      size_t buffer_temp_loc = buffer_loc;
      if (!read_buffer_line(line, sizeof(line), &buffer_temp_loc, buffer,
                            buffer_len)) {
        fail = 1;
        break;
      }

      if (line[0] != '$') {
        if (in_bulk && (strncmp(line, "END BULK", 8) == 0 ||
                        strncmp(line, "ENDDATA", 7) == 0)) {
          buffer_temp_loc = buffer_len;
        } else if (strncmp(line, "GRID*", 5) == 0) {
          char line2[81];
          if (!read_buffer_line(line2, sizeof(line2), &buffer_temp_loc, buffer,
                                buffer_len)) {
            fail = 1;
            break;
          }
          int node;
          double x, y, z;
          parse_node_long_field(line, line2, &node, &x, &y, &z);
          num_nodes++;
        } else if (strncmp(line, "GRID", 4) == 0) {
          int node;
          double x, y, z;
          parse_node_short_free_field(line, &node, &x, &y, &z);
          num_nodes++;
        } else if (strncmp(line, "SPC", 3) == 0) {
          bc_vars_size += 8;
          num_bcs++;
        } else {
          // Check the library of elements
          int max_num_conn = -1;
          int entry_width = 8;

          // Loop over the number of types and determine the number of
          // nodes
          int index = -1;
          for (int k = 0; k < TacsMeshLoaderNumElementTypes; k++) {
            int len = strlen(TacsMeshLoaderElementTypes[k]);
            if (strncmp(line, TacsMeshLoaderElementTypes[k], len) == 0) {
              max_num_conn = TacsMeshLoaderElementLimits[k][1];
              index = k;

              // Check if we should use the extended width or not
              if (line[len] == '*') {
                entry_width = 16;
              }
              break;
            }
          }

          if (index >= 0) {
            // Find the number of entries in the element
            int elem_num, component_num, num_conn;
            buffer_temp_loc = buffer_loc;
            int fail = parse_element_field(&buffer_temp_loc, buffer, buffer_len,
                                           entry_width, max_num_conn, &elem_num,
                                           &component_num, NULL, &num_conn);
            if (fail) {
              break;
            }

            // Check if the number of nodes is within the prescribed limits
            if (num_conn < TacsMeshLoaderElementLimits[index][0]) {
              fprintf(stderr,
                      "TACSMeshLoader: Number of nodes for element %s "
                      "not within limits\n",
                      TacsMeshLoaderElementTypes[index]);
              fail = 1;
              break;
            }

            elem_conn_size += num_conn;
            num_elements++;
            if (num_conn > max_element_conn) {
              max_element_conn = num_conn;
            }
            if (component_num > num_components) {
              num_components = component_num;
            }
          } else {
            fprintf(stderr,
                    "TACSMeshLoader: Element not recognized. Line\n %s\n",
                    line);
          }
        }
      }

      buffer_loc = buffer_temp_loc;
    }

    // Allocate space to store the node numbers
    int *temp_nodes = new int[max_element_conn];

    // Allocate space for everything
    file_node_nums = new int[num_nodes];
    double *file_Xpts = new double[3 * num_nodes];

    // Element type information
    file_elem_nums = new int[num_elements];
    int *file_comp = new int[num_elements];

    // The connectivity information
    int *file_conn = new int[elem_conn_size];
    int *file_conn_ptr = new int[num_elements + 1];
    file_conn_ptr[0] = 0;

    // Boundary condition information
    bc_nodes = new int[num_bcs];
    bc_ptr = new int[num_bcs + 1];
    bc_vals = new TacsScalar[bc_vars_size];
    bc_vars = new int[bc_vars_size];
    bc_ptr[0] = 0;

    // Allocate space for storing the component names
    component_elems = new char[9 * num_components];
    component_descript = new char[33 * num_components];
    memset(component_elems, '\0', 9 * num_components * sizeof(char));
    memset(component_descript, '\0', 33 * num_components * sizeof(char));

    // Reset the sizes of the things to be read in
    num_nodes = 0;
    num_elements = 0;
    num_bcs = 0;
    elem_conn_size = 0;
    bc_vars_size = 0;

    // Rewind to the beginning of the bulk section and allocate everything
    buffer_loc = bulk_start;

    // Keep track of the component numbers loaded from an
    // ICEM-generated bdf file
    int component_counter = 0;

    while (buffer_loc < buffer_len) {
      // Read the first line of the buffer
      size_t buffer_temp_loc = buffer_loc;
      if (!read_buffer_line(line, sizeof(line), &buffer_temp_loc, buffer,
                            buffer_len)) {
        fail = 1;
        break;
      }

      if (strncmp(line, "$       Shell", 13) == 0) {
        // A standard icem output - description of each
        // component. This is very useful for describing what the
        // components actually are with a string.
        // Again use a fixed width format
        char comp[33];
        int comp_num = component_counter;
        component_counter++;

        strncpy(comp, &line[41], 32);
        comp[32] = '\0';
        // Remove white space
        if (comp_num >= 0 && comp_num < num_components) {
          sscanf(comp, "%s", &component_descript[33 * comp_num]);
        }
      }
      if (line[0] != '$') {  // A comment line
        if (in_bulk && (strncmp(line, "END BULK", 8) == 0 ||
                        strncmp(line, "ENDDATA", 7) == 0)) {
          buffer_temp_loc = buffer_len;
        } else if (strncmp(line, "GRID*", 5) == 0) {
          char line2[81];
          if (!read_buffer_line(line2, sizeof(line2), &buffer_temp_loc, buffer,
                                buffer_len)) {
            fail = 1;
            break;
          }
          int node;
          double x, y, z;
          parse_node_long_field(line, line2, &node, &x, &y, &z);
          file_node_nums[num_nodes] = node - 1;  // Get the C ordering
          file_Xpts[3 * num_nodes] = x;
          file_Xpts[3 * num_nodes + 1] = y;
          file_Xpts[3 * num_nodes + 2] = z;
          num_nodes++;
        } else if (strncmp(line, "GRID", 4) == 0) {
          int node;
          double x, y, z;
          parse_node_short_free_field(line, &node, &x, &y, &z);
          file_node_nums[num_nodes] = node - 1;  // Get the C ordering
          file_Xpts[3 * num_nodes] = x;
          file_Xpts[3 * num_nodes + 1] = y;
          file_Xpts[3 * num_nodes + 2] = z;
          num_nodes++;
        } else if (strncmp(line, "SPC", 3) == 0) {
          // This is a variable-length format. Read in grid points until
          // zero is reached. This is a fixed-width format
          // SPC SID  G1  C  D

          // Read in the nodal value
          char node[9];
          strncpy(node, &line[16], 8);
          node[8] = '\0';
          bc_nodes[num_bcs] = atoi(node) - 1;

          strncpy(node, &line[32], 8);
          node[8] = '\0';
          double val = bdf_atof(node);

          // Read in the dof that will be constrained
          for (int k = 24; k < 32; k++) {
            char dofs[9] = "12345678";

            for (int j = 0; j < 8; j++) {
              if (dofs[j] == line[k]) {
                bc_vars[bc_vars_size] = j;
                bc_vals[bc_vars_size] = val;
                bc_vars_size++;
                break;
              }
            }
          }

          bc_ptr[num_bcs + 1] = bc_vars_size;
          num_bcs++;
        } else {
          // Check the library of elements
          int max_num_conn = -1;
          int entry_width = 8;

          // Loop over the number of types and determine the number of
          // nodes
          int index = -1;
          for (int k = 0; k < TacsMeshLoaderNumElementTypes; k++) {
            int len = strlen(TacsMeshLoaderElementTypes[k]);
            if (strncmp(line, TacsMeshLoaderElementTypes[k], len) == 0) {
              max_num_conn = TacsMeshLoaderElementLimits[k][1];
              index = k;

              // Check if we should use the extended width or not
              if (line[len] == '*') {
                entry_width = 16;
              }
              break;
            }
          }

          if (index >= 0) {
            // Find the number of entries in the element
            int elem_num, component_num, num_conn;
            buffer_temp_loc = buffer_loc;
            int fail = parse_element_field(
                &buffer_temp_loc, buffer, buffer_len, entry_width, max_num_conn,
                &elem_num, &component_num, temp_nodes, &num_conn);
            if (fail) {
              break;
            }

            if (strncmp(line, "CQUAD4", 6) == 0 ||
                strncmp(line, "CQUADR", 6) == 0) {
              file_conn[elem_conn_size] = temp_nodes[0] - 1;
              file_conn[elem_conn_size + 1] = temp_nodes[1] - 1;
              file_conn[elem_conn_size + 2] = temp_nodes[3] - 1;
              file_conn[elem_conn_size + 3] = temp_nodes[2] - 1;
            } else if (strncmp(line, "CQUAD9", 6) == 0 ||
                       strncmp(line, "CQUAD", 5) == 0) {
              file_conn[elem_conn_size] = temp_nodes[0] - 1;
              file_conn[elem_conn_size + 1] = temp_nodes[4] - 1;
              file_conn[elem_conn_size + 2] = temp_nodes[1] - 1;
              file_conn[elem_conn_size + 3] = temp_nodes[7] - 1;
              file_conn[elem_conn_size + 4] = temp_nodes[8] - 1;
              file_conn[elem_conn_size + 5] = temp_nodes[5] - 1;
              file_conn[elem_conn_size + 6] = temp_nodes[3] - 1;
              file_conn[elem_conn_size + 7] = temp_nodes[6] - 1;
              file_conn[elem_conn_size + 8] = temp_nodes[2] - 1;
            } else if (strncmp(line, "CHEXA", 5) == 0) {
              file_conn[elem_conn_size] = temp_nodes[0] - 1;
              file_conn[elem_conn_size + 1] = temp_nodes[1] - 1;
              file_conn[elem_conn_size + 2] = temp_nodes[3] - 1;
              file_conn[elem_conn_size + 3] = temp_nodes[2] - 1;
              file_conn[elem_conn_size + 4] = temp_nodes[4] - 1;
              file_conn[elem_conn_size + 5] = temp_nodes[5] - 1;
              file_conn[elem_conn_size + 6] = temp_nodes[7] - 1;
              file_conn[elem_conn_size + 7] = temp_nodes[6] - 1;
            } else {
              for (int k = 0; k < num_conn; k++) {
                file_conn[elem_conn_size + k] = temp_nodes[k] - 1;
              }
            }

            // Set the node numbers
            file_elem_nums[num_elements] = elem_num - 1;
            file_comp[num_elements] = component_num - 1;

            elem_conn_size += num_conn;
            file_conn_ptr[num_elements + 1] = elem_conn_size;
            num_elements++;

            if (component_elems[9 * (component_num - 1)] == '\0') {
              if (strncmp(line, "CTETRA", 6) == 0 && num_conn == 10) {
                strcpy(&component_elems[9 * (component_num - 1)], "CTETRA10");
              } else {
                strcpy(&component_elems[9 * (component_num - 1)],
                       TacsMeshLoaderElementTypes[index]);
              }
            }
          } else {
            fprintf(stderr,
                    "TACSMeshLoader: Element not recognized. Line\n %s\n",
                    line);
          }
        }
      }

      buffer_loc = buffer_temp_loc;
    }

    delete[] buffer;
    delete[] temp_nodes;

    if (fail) {
      MPI_Abort(comm, fail);
      return fail;
    }

    // Arg sort the list of nodes
    node_arg_sort_list = new int[num_nodes];
    for (int k = 0; k < num_nodes; k++) {
      node_arg_sort_list[k] = k;
    }

    arg_sort_list = file_node_nums;
    qsort(node_arg_sort_list, num_nodes, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Arg sort the list of elements
    elem_arg_sort_list = new int[num_elements];
    for (int k = 0; k < num_elements; k++) {
      elem_arg_sort_list[k] = k;
    }

    arg_sort_list = file_elem_nums;
    qsort(elem_arg_sort_list, num_elements, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Now file_node_nums[node_arg_sort_list[k]] and
    // file_elem_nums[elem_arg_sort_list[k] are sorted.

    // Create the output for the nodes
    Xpts = new TacsScalar[3 * num_nodes];

    for (int k = 0; k < num_nodes; k++) {
      int n = node_arg_sort_list[k];
      for (int j = 0; j < 3; j++) {
        Xpts[3 * k + j] = file_Xpts[3 * n + j];
      }
    }
    delete[] file_Xpts;

    // Read in the connectivity array and store the information
    elem_node_conn = new int[elem_conn_size];
    elem_node_ptr = new int[num_elements + 1];
    elem_component = new int[num_elements];

    // Now, loop over all of the elements and all of the connectivity
    // list and find the
    elem_node_ptr[0] = 0;
    for (int k = 0, n = 0; k < num_elements; k++) {
      int e = elem_arg_sort_list[k];

      for (int j = file_conn_ptr[e]; j < file_conn_ptr[e + 1]; j++, n++) {
        int node_num = file_conn[j];

        // Find node_num in the list
        int node = find_index_arg_sorted(node_num, num_nodes, file_node_nums,
                                         node_arg_sort_list);
        if (node < 0) {
          elem_node_conn[n] = -1;
          fail = 1;
        } else {
          elem_node_conn[n] = node;
        }
      }

      elem_component[k] = file_comp[e];
      elem_node_ptr[k + 1] = n;
    }

    // Find the boundary condition nodes
    for (int k = 0; k < num_bcs; k++) {
      int node = find_index_arg_sorted(bc_nodes[k], num_nodes, file_node_nums,
                                       node_arg_sort_list);
      if (node < 0) {
        fail = 1;
        bc_nodes[k] = -1;
      } else {
        bc_nodes[k] = node;
      }
    }

    // Free data that has been allocated locally
    delete[] file_comp;
    delete[] file_conn;
    delete[] file_conn_ptr;
  }

  // Distribute the component numbers and descritpions
  MPI_Bcast(&num_components, 1, MPI_INT, root, comm);
  if (rank != root) {
    component_elems = new char[9 * num_components];
    component_descript = new char[33 * num_components];
  }
  MPI_Bcast(component_elems, 9 * num_components, MPI_CHAR, root, comm);
  MPI_Bcast(component_descript, 33 * num_components, MPI_CHAR, root, comm);

  elements = new TACSElement *[num_components];
  for (int k = 0; k < num_components; k++) {
    elements[k] = NULL;
  }

  return fail;
}

/*
  Retrieve the number of nodes in the model
*/
int TACSMeshLoader::getNumNodes() { return num_nodes; }

/*
  Create a TACSToFH5 file creation object
*/
TACSToFH5 *TACSMeshLoader::createTACSToFH5(TACSAssembler *tacs,
                                           ElementType elem_type,
                                           int write_flag) {
  // Set the component numbers in the elements
  for (int k = 0; k < num_components; k++) {
    elements[k]->setComponentNum(k);
  }

  TACSToFH5 *f5 = new TACSToFH5(tacs, elem_type, write_flag);
  for (int k = 0; k < num_components; k++) {
    if (strlen(&component_descript[33 * k]) == 0) {
      char name[64];
      sprintf(name, "Component %d", k);
      f5->setComponentName(k, name);
    } else {
      f5->setComponentName(k, &component_descript[33 * k]);
    }
  }

  return f5;
}

/*
  Create a distributed version of TACS
*/
TACSAssembler *TACSMeshLoader::createTACS(
    int vars_per_node, TACSAssembler::OrderingType order_type,
    TACSAssembler::MatrixOrderingType mat_type) {
  // Set the root processor
  const int root = 0;

  // Get the rank of the current processor
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Allocate the TACS creator
  creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  // Set the ordering type and matrix type
  creator->setReorderingType(order_type, mat_type);

  if (rank == root) {
    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elements, elem_node_ptr,
                                   elem_node_conn, elem_component);

    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes, bc_ptr, bc_vars, bc_vals);

    // Set the nodal locations
    creator->setNodes(Xpts);

    // Free things that are no longer required
    delete[] elem_node_ptr;
    elem_node_ptr = NULL;
    delete[] elem_node_conn;
    elem_node_conn = NULL;
    delete[] elem_component;
    elem_component = NULL;

    // Free the boundary conditions
    delete[] bc_nodes;
    bc_nodes = NULL;
    delete[] bc_ptr;
    bc_ptr = NULL;
    delete[] bc_vars;
    bc_vars = NULL;
    delete[] bc_vals;
    bc_vals = NULL;
  }

  // This call must occur on all processor
  creator->setElements(num_components, elements);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();

  return tacs;
}

/*
  Retrieve the number of elements owned by this processes
*/
int TACSMeshLoader::getNumElements() { return num_elements; }

/*
  Set the function domain

  Given the function, and the set of component numbers that define the
  domain of interest, set the element numbers in the function that
*/
void TACSMeshLoader::addFunctionDomain(TACSFunction *function, int num_comps,
                                       int comp_nums[]) {
  if (creator) {
    int *elems;
    int num_elems = creator->getElementIdNums(num_comps, comp_nums, &elems);
    function->addDomain(num_elems, elems);
    delete[] elems;
  }
}

/*
  Add the auxiliary element to the given domain specified by the
  component number
*/
void TACSMeshLoader::addAuxElement(TACSAuxElements *aux, int component_num,
                                   TACSElement *element) {
  if (creator) {
    int *elems;
    int num_elems = creator->getElementIdNums(1, &component_num, &elems);
    for (int i = 0; i < num_elems; i++) {
      aux->addElement(elems[i], element);
    }
    delete[] elems;
  }
}

/**
  Given node numbers from the original file on the root processor,
  find the corresponding global node numbers in the given assembler object.

  Note that the node numbers are assumed to be 1-based as is the case in the
  original file format. In addition, the node array is over-written by a
  temporary ordering. The number of nodes and their numbers are returned in
  a newly allocated array.
*/
void TACSMeshLoader::getAssemblerNodeNums(TACSAssembler *assembler,
                                          int num_nodes, int *node_nums,
                                          int *num_new_nodes, int **new_nodes) {
  *num_new_nodes = 0;
  *new_nodes = NULL;

  if (creator) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    int index = 0;
    if (rank == 0) {
      // Convert from the BDF order, to the local TACSMeshLoader order
      for (int k = 0; k < num_nodes; k++) {
        int node_num = node_nums[k] - 1;
        int node = find_index_arg_sorted(node_num, num_nodes, file_node_nums,
                                         node_arg_sort_list);
        if (node >= 0) {
          node_nums[index] = node;
          index++;
        }
      }
    }

    creator->getAssemblerNodeNums(assembler, index, node_nums, num_new_nodes,
                                  new_nodes);
  }
}

/*
  Get the element connectivity and node locations
*/
void TACSMeshLoader::getConnectivity(int *_num_nodes, int *_num_elements,
                                     const int **_elem_node_ptr,
                                     const int **_elem_node_conn,
                                     const int **_elem_component,
                                     const TacsScalar **_Xpts) {
  if (_num_nodes) {
    *_num_nodes = num_nodes;
  }
  if (_num_elements) {
    *_num_elements = num_elements;
  }
  if (_elem_node_ptr) {
    *_elem_node_ptr = elem_node_ptr;
  }
  if (_elem_node_conn) {
    *_elem_node_conn = elem_node_conn;
  }
  if (_elem_component) {
    *_elem_component = elem_component;
  }
  if (_Xpts) {
    *_Xpts = Xpts;
  }
}

/*
  Get the boundary conditions and data
*/
void TACSMeshLoader::getBCs(int *_num_bcs, const int **_bc_nodes,
                            const int **_bc_vars, const int **_bc_ptr,
                            const TacsScalar **_bc_vals) {
  if (_num_bcs) {
    *_num_bcs = num_bcs;
  }
  if (_bc_nodes) {
    *_bc_nodes = bc_nodes;
  }
  if (_bc_vars) {
    *_bc_vars = bc_vars;
  }
  if (_bc_ptr) {
    *_bc_ptr = bc_ptr;
  }
  if (_bc_vals) {
    *_bc_vals = bc_vals;
  }
}
