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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "TACSMeshLoader.h"
#include "FElibrary.h"

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
static const int * arg_sort_list = NULL;

static int compare_arg_sort( const void * a, const void * b ){
  return arg_sort_list[*(int*)a] - arg_sort_list[*(int*)b];
}

/*
  Read a line from the buffer.

  Return the number of read characters. Do not exceed the buffer
  length.

  Given the line buffer 'line', and the size of the line buffer line_len
*/
static int read_buffer_line( char *line, size_t line_len, 
                             size_t *loc, char *buffer, 
                             size_t buffer_len ){
  size_t i = 0;
  for ( ; (i < line_len) && (*loc < buffer_len); i++, (*loc)++ ){
    if (buffer[*loc] == '\n'){
      break;
    }
    line[i] = buffer[*loc];    
  }
  if (i < line_len){
    line[i] = '\0';
  }

  // Read until the end of the line
  while ((*loc < buffer_len) && (buffer[*loc] != '\n')){
    (*loc)++;
  }

  (*loc)++; // Increment so that loc is at the next location
  
  return i;
}

/*
  Reverse look up.

  Given the sorted array list[arg[k]], find k, such that 

  list[arg[k]] = var
*/
static int find_index_arg_sorted( int var, int size, 
                                  const int * list, 
                                  const int * args ){
  // Binary search an array to find k such that list[k] = var,
  // where the array list[args[k]] is sorted in ascending
  // order
  int high = size-1;
  int low = 0;
  int high_val = list[args[high]];
  int low_val = list[args[low]];

  // Check if the index is at the end points
  if (var == low_val){
    return low;
  }
  else if (var < low_val){
    return -1;
  }

  if (var == high_val){
    return high;
  }
  else if (var > high_val){
    return -1;
  }

  int mid = low + (int)((high - low)/2);
      
  // While there are values left in the list
  while (low != mid){
    int mid_val = list[args[mid]];
    if (mid_val == var){
      return mid;
    }   
    
    if (var < mid_val){
      high = mid;
      high_val = mid_val;
    }
    else {
      low = mid;
      low_val = mid_val;
    }
    
    mid = low + (int)((high - low)/2);
  }      
  
  return -1;
}

/*
  Convert a Nastran-style number with an exponent to a double.
*/
static double bdf_atof( char * str ){
  // First, check if the string contains an E/e or D/d - if so, convert it
  int slen = strlen(str);
  for ( int i = 0; i < slen; i++ ){
    if (str[i] == 'e' || str[i] == 'E'){
      return atof(str);
    }
    else if (str[i] == 'd' || str[i] == 'D'){
      str[i] = 'e';
      return atof(str);
    }
  }

  // Convert the special Nastran number format without e/E or d/D 
  // 4.-4 or 5.34+2 etc.
  char temp[24];
  int i = 0, j = 0;
  while (i < slen && str[i] == ' '){ i++; }
  if (i == slen){ 
    return 0.0;
  }

  // Take care of a leading minus sign
  if (str[i] == '-' ){ 
    temp[j] = str[i];
    j++, i++; 
  }

  for ( ; i < slen; i++, j++ ){
    // Add the character 'e' before the exponent
    if (str[i] == '-' || str[i] == '+'){
      temp[j] = 'e'; j++; 
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
static void parse_node_long_field( char *line, char *line2, int *node, 
                                   double *x, double *y, double *z ){
  char Node[32], X[32], Y[32], Z[32];

  strncpy(Node, &line[8], 16); Node[16] = '\0';
  strncpy(X, &line[40], 16);   X[16] = '\0';
  strncpy(Y, &line[56], 16);   Y[16] = '\0';
  strncpy(Z, &line2[8], 16);   Z[16] = '\0';

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
static void parse_node_short_free_field( char *line, int *node, 
                                         double *x, double *y, double *z ){
  char field[5][32];

  // Look for a comma 
  int len = strlen(line);
  int comma_format = 0;
  for ( int i = 0; (i < len) && (i < 80); i++ ){
    if (line[i] == ','){ 
      comma_format = 1;
      break;
    }
  }

  if (comma_format){
    int start = 8;
    int end = 8;
    for ( int i = 0; i < 5; i++ ){
      start = end;
      while (end < len && line[end] != ',')
        end++;
      int flen = end-start;
      strncpy(field[i], &line[start], flen);
      field[i][flen] = '\0';
      end++;
    }
  }
  else { // Short-format, fixed width
    strncpy(field[0], &line[8], 8);   field[0][8] = '\0';
    strncpy(field[2], &line[24], 8);  field[2][8] = '\0';
    strncpy(field[3], &line[32], 8);  field[3][8] = '\0';
    strncpy(field[4], &line[40], 8);  field[4][8] = '\0';
  }

  *node = atoi(field[0]);
  *x = bdf_atof(field[2]);
  *y = bdf_atof(field[3]);
  *z = bdf_atof(field[4]);
}

static void parse_element_field( char line[], 
                                 int * elem_num, int * component_num,
                                 int * node_nums, int num_nodes ){
  char node[9];
  int entry = 8;

  strncpy(node, &line[entry], 8);
  node[8] = '\0';
  *elem_num = atoi(node);
  entry += 8;
  
  strncpy(node, &line[entry], 8);
  node[8] = '\0';
  *component_num = atoi(node);
  entry += 8;
  
  if (*component_num <= 0){
    fprintf(stderr, 
            "Error: The component numbers must be strictly positive\n");
  }
  
  for ( int n = 0; n < num_nodes && entry < 80; entry += 8, n++ ){
    // Parse the line containing the entry
    strncpy(node, &line[entry], 8);
    node[8] = '\0';
    node_nums[n] = atoi(node);
  }
}

static void parse_element_field2( char line1[], char line2[],
                                  int * elem_num, int * component_num,
                                  int * node_nums, int num_nodes, int width=8 ){

  int n = 0; // The number of parsed nodes
  char node[17];
  for ( int m = 0; m < 2; m++ ){
    int entry = width;
    const char * line = line1;
    if (m == 1){
      line = line2;
    }

    if (n == 0){ 
      strncpy(node, &line[entry], width);
      node[width] = '\0';
      *elem_num = atoi(node);
      entry += width;

      strncpy(node, &line[entry], width);
      node[width] = '\0';
      *component_num = atoi(node);
      entry += width;
    }
    
    for ( ; n < num_nodes && entry < 72; entry += width, n++ ){
      // Parse the line containing the entry
      strncpy(node, &line[entry], width);
      node[width] = '\0';
      node_nums[n] = atoi(node);
    }
  }
}


static void parse_element_field3( char line1[], char line2[], char line3[],
                                  int * elem_num, int * component_num,
                                  int * node_nums, int num_nodes, int width=8 ){  
  int n = 0; // The number of parsed nodes
  char node[17];

  for ( int m = 0; m < 3; m++ ){
    int entry = width;
    const char * line = line1;
    if (m == 1){
      line = line2;
    }
    else if (m == 2){
      line = line3;
    }

    if (n == 0){ 
      strncpy(node, &line[entry], width);
      node[width] = '\0';
      *elem_num = atoi(node);
      entry += width;

      strncpy(node, &line[entry], width);
      node[width] = '\0';
      *component_num = atoi(node);
      entry += width;
    }
    
    for ( ; n < num_nodes && entry < 72; entry += width, n++ ){
      // Parse the line containing the entry
      strncpy(node, &line[entry], width);
      node[width] = '\0';
      node_nums[n] = atoi(node);
    }
  }
}

static void parse_element_field4( char line1[], char line2[], char line3[],
                                  char line4[],
                                  int * elem_num, int * component_num,
                                  int * node_nums, int num_nodes ){  
  int n = 0; // The number of parsed nodes
  char node[9];

  for ( int m = 0; m < 4; m++ ){
    int entry = 8;
    const char * line = line1;
    if (m == 1){
      line = line2;
    }
    else if (m == 2){
      line = line3;
    }
    else if (m == 3){
      line = line4;
    }

    if (n == 0){ 
      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      *elem_num = atoi(node);
      entry += 8;

      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      *component_num = atoi(node);
      entry += 8;
    }

    for ( ; n < num_nodes && entry < 72; entry += 8, n++ ){
      // Parse the line containing the entry
      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      node_nums[n] = atoi(node);
    }
  }
}

static void parse_element_field9( char line1[], char line2[], char line3[],
                                  char line4[], char line5[], char line6[],
                                  char line7[], char line8[], char line9[],
                                  int * elem_num, int * component_num,
                                  int * node_nums, int num_nodes ){  
  int n = 0; // The number of parsed nodes
  char node[9];

  for ( int m = 0; m < 4; m++ ){
    int entry = 8;
    const char * line = line1;
    if (m == 1){ line = line2; }
    else if (m == 2){ line = line3; }
    else if (m == 3){ line = line4; }
    else if (m == 4){ line = line5; }
    else if (m == 5){ line = line6; }
    else if (m == 6){ line = line7; }
    else if (m == 7){ line = line8; }
    else if (m == 8){ line = line9; }
    
    if (n == 0){ 
      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      *elem_num = atoi(node);
      entry += 8;

      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      *component_num = atoi(node);
      entry += 8;
    }

    for ( ; n < num_nodes && entry < 72; entry += 8, n++ ){
      // Parse the line containing the entry
      strncpy(node, &line[entry], 8);
      node[8] = '\0';
      node_nums[n] = atoi(node);
    }
  }
}

/*
  Converts the connectivity information loaded from BDF file to
  coordinate ordering used in TACS.
*/
static void convert_to_coordinate( int * coord, int * orig ){
  coord[0] = orig[0];
  coord[1] = orig[4];
  coord[2] = orig[1];

  coord[3] = orig[7];
  coord[4] = orig[8];
  coord[5] = orig[5];

  coord[6] = orig[3];
  coord[7] = orig[6];
  coord[8] = orig[2];  
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
TACSMeshLoader::TACSMeshLoader( MPI_Comm _comm ){
  comm = _comm;

  // Initialize everything to zero
  num_nodes = num_elements = 0;
  num_bcs = 0;
  node_nums = NULL;
  Xpts_unsorted = NULL;
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
  
  // Default is not to convert to coordinate, the supplied BDF is
  // assumed in order
  convertToCoordinate = 0;
}

/*
  Destroy existing data that may have been allocated

  Note that most of this data will only be allocated on the root
  processor 
*/
TACSMeshLoader::~TACSMeshLoader(){
  if (elem_node_conn){ delete [] elem_node_conn; }
  if (elem_node_ptr){ delete [] elem_node_ptr; }
  if (elem_component){ delete [] elem_component; }
  if (Xpts){ delete [] Xpts; }
  if (bc_nodes){ delete [] bc_nodes; }
  if (bc_vars){ delete [] bc_vars; }
  if (bc_vals){ delete [] bc_vals; }
  if (bc_ptr){ delete [] bc_ptr; }

  if (elements){
    for ( int k = 0; k < num_components; k++ ){
      if (elements[k]){ elements[k]->decref(); }
    }
    delete [] elements;
  }
  if (component_elems){ delete [] component_elems; }
  if (component_descript){ delete [] component_descript; }
  if (Xpts_unsorted){ delete [] Xpts_unsorted;}
  if (node_nums) {delete [] node_nums;}

  // Free the creator object
  if (creator){ creator->decref(); }
}

/*
  Get the number of components defined by the data
*/
int TACSMeshLoader::getNumComponents(){
  return num_components;
}

/*
  Set the element associated with a given component number
*/
void TACSMeshLoader::setElement( int component_num, 
                                 TACSElement *_element ){
  if (_element && (component_num >= 0) && (component_num < num_components)){
    _element->incref();
    _element->setComponentNum(component_num);
    elements[component_num] = _element;
  }
}

/*
  Set whether to convert to coordinate ordering before creating
  elements
*/
void TACSMeshLoader::setConvertToCoordinate( int flag ){
  convertToCoordinate = flag;
}

/*
  Get the component description from the file
*/
const char *TACSMeshLoader::getComponentDescript( int comp_num ){
  if (component_descript && (comp_num >= 0) && 
      (comp_num < num_components)){
    return &component_descript[33*comp_num];
  }
  return NULL;
}

/*
  Retrieve the element description corresponding to the component number
*/
const char *TACSMeshLoader::getElementDescript( int comp_num ){
  if (component_elems && (comp_num >= 0) && 
      (comp_num < num_components)){
    return &component_elems[9*comp_num];
  }
  return NULL; // No associated element
}

/*
  This functions scans a Nastran BDF file - only scanning in
  information from the bulk data section.

  Only the element types, boundary conditions, connectivitiy and GRID
  entries are scanned.  Any entries associated with constitutive
  properties are ignored.
*/
int TACSMeshLoader::scanBDFFile( const char * file_name ){
  int rank;
  MPI_Comm_rank(comm, &rank);
  int fail = 0;

  const int root = 0;
  if (rank == root){
    FILE * fp = fopen(file_name, "r"); 
    if (!fp){ 
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
    int elem_con_size = 0;
    int bc_vars_size = 0;

    // Each line can only be 80 characters long
    char line[9][80];

    // Determine the size of the file
    fseek(fp, 0, SEEK_END);
    size_t buffer_len = ftell(fp);
    rewind(fp);

    // Allocate enough space to store the entire file
    char * buffer = new char[buffer_len];
    if (fread(buffer, 1, buffer_len, fp) != buffer_len){
      fprintf(stderr, "[%d] TACSMeshLoader: Problem reading file %s\n",
              rank, file_name);
      MPI_Abort(comm, 1);
      return 1;
    }
    fclose(fp);

    // Keep track of where the current point in the buffer is
    size_t buffer_loc = 0;
    read_buffer_line(line[0], sizeof(line[0]), 
                     &buffer_loc, buffer, buffer_len);

    // Flags which indicate where the bulk data begins
    int in_bulk = 0;
    int bulk_start = 0;
   
    // Scan the file for the begin bulk location. If none exists, then
    // the whole file is treated as bulk data.
    while (buffer_loc < buffer_len){
      if (strncmp(line[0], "BEGIN BULK", 10) == 0){
        in_bulk = 1;
        bulk_start = buffer_loc;
      }
      read_buffer_line(line[0], sizeof(line[0]), 
                       &buffer_loc, buffer, buffer_len);
    }

    // Treat the entire file as bulk data 
    if (!in_bulk){
      in_bulk = 1;
    }

    // Set the starting location
    buffer_loc = bulk_start;

    while (buffer_loc < buffer_len){
      if (line[0][0] != '$'){ // This is not a comment line
        int node;
        double x, y, z;

        // Check for GRID or GRID*
        if (strncmp(line[0], "GRID*", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1;
            break;
          }
          parse_node_long_field(line[0], line[1], &node, &x, &y, &z);
          num_nodes++;
        }
        else if (strncmp(line[0], "GRID", 4) == 0){
          parse_node_short_free_field(line[0], &node, &x, &y, &z);
          num_nodes++;
        }
        else if (strncmp(line[0], "CBAR", 4) == 0){
          // Read in the component number and nodes associated with
          // this element
          int elem_num, component_num;
          int nodes[2]; // Should have at most four nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 2);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 2;
          num_elements++;
        }
        
	else if (strncmp(line[0], "CHEXA*", 6) == 0){
	  for ( int i = 1; i < 3; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }
          // Read in the component number and nodes associated with
          // this element
          int elem_num, component_num, nodes[8];
	  parse_element_field3(line[0], line[1], line[2],
                               &elem_num, &component_num, 
                               nodes, 8, 16);
	  if (component_num > num_components){
            num_components = component_num;
          }
	  elem_con_size += 8;
          num_elements++;
	}
        else if (strncmp(line[0], "CHEXA", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated with
          // this element
          int elem_num, component_num, nodes[8]; 
          parse_element_field2(line[0], line[1], 
                               &elem_num, &component_num, 
                               nodes, 8);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 8;
          num_elements++;
        }
        else if (strncmp(line[0], "CHEXA27", 7) == 0){
          for ( int i = 1; i < 4; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }

          // Read in the component number and nodes associated with
          // this element
          int elem_num, component_num, nodes[27];
          parse_element_field4(line[0], line[1], line[2], line[3],
                               &elem_num, &component_num, 
                               nodes, 27);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 27;
          num_elements++;
        }
        else if (strncmp(line[0], "CHEXA64", 7) == 0){
          for ( int i = 1; i < 9; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num, nodes[27];
          parse_element_field9(line[0], line[1], line[2], line[3], line[4],
                               line[5], line[6], line[7], line[8],
                               &elem_num, &component_num, 
                               nodes, 64);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 64;
          num_elements++;
        }
        else if (strncmp(line[0], "CQUAD16", 7) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }
          if (!read_buffer_line(line[2], sizeof(line[2]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[16]; // Should have at most four nodes
          parse_element_field3(line[0], line[1], line[2],
                               &elem_num, &component_num,
                               nodes, 16);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 16;
          num_elements++;
        }
        else if (strncmp(line[0], "CQUAD9", 6) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[9]; // Should have at most 9 nodes          
          if (!convertToCoordinate) {
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 nodes, 9);
          } 
          else {
            int tmp[9];
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 tmp, 9);
            // convert to coordinate ordering for gmsh
            convert_to_coordinate(&nodes[0], &tmp[0]);
          }

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 9;
          num_elements++;
        }
        else if (strncmp(line[0], "CQUAD4*", 7) == 0 ){
          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[4]; // Should have at most four nodes
          parse_element_field2(line[0], line[1],
                              &elem_num, &component_num,
                               nodes, 4, 16);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 4;
          num_elements++;
        }
        else if (strncmp(line[0], "CQUAD4", 6) == 0 || 
                 strncmp(line[0], "CQUADR", 6) == 0){
          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[4]; // Should have at most four nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 4);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 4;
          num_elements++;
        }  
        else if (strncmp(line[0], "CQUAD", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[9]; // Should have at most four nodes
          if (!convertToCoordinate) {
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 nodes, 9);
          } 
          else {
            int tmp[9];
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 tmp, 9);
            // convert to coordinate ordering for gmsh
            convert_to_coordinate(&nodes[0], &tmp[0]);
          }

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 9;
          num_elements++;
        }      
        else if (strncmp(line[0], "CTRIA3", 6) == 0){
          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[3]; // Should have at most three nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 3);

          if (component_num > num_components){
            num_components = component_num;
          }

          elem_con_size += 3;
          num_elements++;
        }  
        else if (strncmp(line[0], "SPC", 3) == 0){
          bc_vars_size += 8;
          num_bcs++;
        }
        else if (strncmp(line[0], "FFORCE", 6) == 0){     
          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[1]; // Should have at most four nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 1);

          if (component_num > num_components){
            num_components = component_num;
          }
          elem_con_size += 1;
          num_elements++;
        }
      }

      read_buffer_line(line[0], sizeof(line[0]), 
                       &buffer_loc, buffer, buffer_len);
    }
  
    // Allocate space for everything    
    node_nums = new int[ num_nodes ];
    Xpts_unsorted = new double[ 3*num_nodes ];
    
    // Element type information
    int * elem_nums = new int[ num_elements ];
    int * elem_comp = new int[ num_elements ];
    
    // The connectivity information
    int * elem_con = new int[ elem_con_size ];
    int * elem_con_ptr = new int[ num_elements+1 ];
    elem_con_ptr[0] = 0;
    
    // Boundary condition information
    bc_nodes = new int[ num_bcs ];
    bc_ptr = new int[ num_bcs+1 ];
    bc_vals = new TacsScalar[ bc_vars_size ];
    bc_vars = new int[ bc_vars_size ];
    bc_ptr[0] = 0;

    // Allocate space for storing the component names
    component_elems = new char[ 9*num_components ];
    component_descript = new char[ 33*num_components ];
    memset(component_elems, '\0', 9*num_components*sizeof(char));
    memset(component_descript, '\0', 33*num_components*sizeof(char));

    // Reset the sizes of the things to be read in
    num_nodes = 0;
    num_elements = 0;
    num_bcs = 0;
    elem_con_size = 0;
    bc_vars_size = 0;

    // Rewind to the beginning of the bulk section and allocate everything
    buffer_loc = bulk_start;
    read_buffer_line(line[0], sizeof(line[0]), 
                     &buffer_loc, buffer, buffer_len);

    // Keep track of the component numbers loaded from an
    // ICEM-generated bdf file
    int component_counter = 0;

    while (buffer_loc < buffer_len){        
      int node;
      double x, y, z;

      if (strncmp(line[0], "$       Shell", 13) == 0){
        // A standard icem output - description of each
        // component. This is very useful for describing what the
        // components actually are with a string.
        // Again use a fixed width format
        char comp[33];
        int comp_num = component_counter; 
        component_counter++;

        strncpy(comp, &line[0][41], 32);
        comp[32] = '\0';
        // Remove white space
        if (comp_num >= 0 && comp_num < num_components){
          sscanf(comp, "%s", &component_descript[33*comp_num]);
        }
      }
      if (line[0][0] != '$'){ // A comment line
        // Check for GRID or GRID*
        if (strncmp(line[0], "GRID*", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1;
            break;
          }
          parse_node_long_field(line[0], line[1], &node, &x, &y, &z);
          node_nums[num_nodes] = node-1; // Get the C ordering
          Xpts_unsorted[3*num_nodes]   = x;
          Xpts_unsorted[3*num_nodes+1] = y;
          Xpts_unsorted[3*num_nodes+2] = z;
          num_nodes++;
        }
        else if (strncmp(line[0], "GRID", 4) == 0){
          parse_node_short_free_field(line[0], &node, &x, &y, &z);
          node_nums[num_nodes] = node-1; // Get the C ordering
          Xpts_unsorted[3*num_nodes]   = x;
          Xpts_unsorted[3*num_nodes+1] = y;
          Xpts_unsorted[3*num_nodes+2] = z;     
          num_nodes++;
        }        
        else if (strncmp(line[0], "CBAR", 4) == 0){
          // Read in the component number and nodes associated with
          // this element
          int elem_num, component_num;
          int nodes[2]; // Should have at most two nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 2);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;

          elem_con[elem_con_size] = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;

          elem_con_size += 2;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CBAR");
          }
        }
        
	else if (strncmp(line[0], "CHEXA*", 6) == 0){
	  for ( int i = 1; i < 3; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }
          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num, nodes[8]; 
          parse_element_field3(line[0], line[1], line[2], 
                               &elem_num, &component_num, 
                               nodes, 8, 16);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;

          elem_con[elem_con_size]   = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;
          elem_con[elem_con_size+2] = nodes[3]-1;
          elem_con[elem_con_size+3] = nodes[2]-1;
          elem_con[elem_con_size+4] = nodes[4]-1;
          elem_con[elem_con_size+5] = nodes[5]-1;
          elem_con[elem_con_size+6] = nodes[7]-1;
          elem_con[elem_con_size+7] = nodes[6]-1;
	  
          elem_con_size += 8;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;
	  if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CHEXA*");
          }
        }
        else if (strncmp(line[0], "CHEXA", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num, nodes[8]; 
          parse_element_field2(line[0], line[1], 
                               &elem_num, &component_num, 
                               nodes, 8);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;

          elem_con[elem_con_size]   = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;
          elem_con[elem_con_size+2] = nodes[3]-1;
          elem_con[elem_con_size+3] = nodes[2]-1;
          elem_con[elem_con_size+4] = nodes[4]-1;
          elem_con[elem_con_size+5] = nodes[5]-1;
          elem_con[elem_con_size+6] = nodes[7]-1;
          elem_con[elem_con_size+7] = nodes[6]-1;
         
          elem_con_size += 8;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CHEXA");
          }
        }
        else if (strncmp(line[0], "CHEXA27", 7) == 0){
          for ( int i = 1; i < 4; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }

          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num, nodes[27];
          parse_element_field4(line[0], line[1], line[2], line[3],
                               &elem_num, &component_num, 
                               nodes, 27);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          
          for ( int k = 0; k < 27; k++ ){
            elem_con[elem_con_size+k] = nodes[k]-1;
          }

          elem_con_size += 27;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CHEXA64");
          }
        }
        else if (strncmp(line[0], "CHEXA64", 7) == 0){
          for ( int i = 1; i < 9; i++ ){
            if (!read_buffer_line(line[i], sizeof(line[i]), 
                                  &buffer_loc, buffer, buffer_len)){
              fail = 1; break;
            }
          }

          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num, nodes[64];
          parse_element_field9(line[0], line[1], line[2], line[3], line[4],
                               line[5], line[6], line[7], line[8],
                               &elem_num, &component_num, 
                               nodes, 64);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          
          for ( int k = 0; k < 64; k++ ){
            elem_con[elem_con_size+k] = nodes[k]-1;
          }

          elem_con_size += 64;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CHEXA64");
          }
        }
        else if (strncmp(line[0], "CQUAD16", 7) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }
          if (!read_buffer_line(line[2], sizeof(line[2]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[16]; // Should have at most four nodes
          parse_element_field3(line[0], line[1], line[2],
                               &elem_num, &component_num,
                               nodes, 16);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          
          for ( int k = 0; k < 16; k++ ){
            elem_con[elem_con_size+k] = nodes[k]-1;
          }

          elem_con_size += 16;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CQUAD16");
          }
        }
        else if (strncmp(line[0], "CQUAD9", 6) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num;
          int nodes[9]; // Should have at most 9 nodes          
          if (!convertToCoordinate) {
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 nodes, 9);
          } 
          else {
            int tmp[9];
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 tmp, 9);
            // convert to coordinate ordering for gmsh
            convert_to_coordinate(&nodes[0], &tmp[0]);
          }

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          
          for ( int k = 0; k < 9; k++ ){
            elem_con[elem_con_size+k] = nodes[k]-1;
          }

          elem_con_size += 9;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CQUAD9");
          }
        }
        else if (strncmp(line[0], "CQUAD4*", 7) == 0 ){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1;
            break;
          }
          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num;
          int nodes[4]; // Should have at most four nodes
          parse_element_field2(line[0], line[1],
                               &elem_num, &component_num,
                               nodes, 4, 16);
          // Add the element to the connectivity list
          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          elem_con[elem_con_size]   = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;
          elem_con[elem_con_size+2] = nodes[3]-1;
          elem_con[elem_con_size+3] = nodes[2]-1;
          elem_con_size += 4;

          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CQUAD4");
          }
        }
        else if (strncmp(line[0], "CQUAD4", 6) == 0 ||
                 strncmp(line[0], "CQUADR", 6) == 0){
          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num;
          int nodes[4]; // Should have at most four nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 4);
          // Add the element to the connectivity list
          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;

          elem_con[elem_con_size]   = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;
          elem_con[elem_con_size+2] = nodes[3]-1;
          elem_con[elem_con_size+3] = nodes[2]-1;
          elem_con_size += 4;

          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CQUAD4");
          }
        }        
        else if (strncmp(line[0], "CQUAD", 5) == 0){
          if (!read_buffer_line(line[1], sizeof(line[1]), 
                                &buffer_loc, buffer, buffer_len)){
            fail = 1; break;
          }

          // Read in the component number and nodes associated 
          // with this element
          int elem_num, component_num;
          int nodes[9]; // Should have at most four nodes
          if (!convertToCoordinate) {
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 nodes, 9);
          } 
          else {
            int tmp[9];
            parse_element_field2(line[0], line[1],
                                 &elem_num, &component_num,
                                 tmp, 9);
            // convert to coordinate ordering for gmsh
            convert_to_coordinate(&nodes[0], &tmp[0]);
          }

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;
          
          elem_con[elem_con_size] = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[4]-1;
          elem_con[elem_con_size+2] = nodes[1]-1;

          elem_con[elem_con_size+3] = nodes[7]-1;
          elem_con[elem_con_size+4] = nodes[8]-1;
          elem_con[elem_con_size+5] = nodes[5]-1;

          elem_con[elem_con_size+6] = nodes[3]-1;
          elem_con[elem_con_size+7] = nodes[6]-1;
          elem_con[elem_con_size+8] = nodes[2]-1;          

          elem_con_size += 9;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CQUAD");
          }
        }
        else if (strncmp(line[0], "CTRIA3", 6) == 0){
          // Read in the component number and nodes associated
          // with this element
          int elem_num, component_num;
          int nodes[3]; // Should have at most four nodes
          parse_element_field(line[0],
                              &elem_num, &component_num,
                              nodes, 3);
          
          // Add the element to the connectivity list
          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;

          elem_con[elem_con_size]   = nodes[0]-1;
          elem_con[elem_con_size+1] = nodes[1]-1;
          elem_con[elem_con_size+2] = nodes[2]-1;
          elem_con_size += 3;

          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "CTRIA3");
          }
        }
        else if (strncmp(line[0], "FFORCE", 6) == 0){
          // Read in the component number and nodes associated 
          // with the following force
          int elem_num, component_num;
          int nodes[1]; // Should have at most four nodes
          parse_element_field(line[0], &elem_num, &component_num,
                              nodes, 1);

          elem_nums[num_elements] = elem_num-1;
          elem_comp[num_elements] = component_num-1;          
          elem_con[elem_con_size] = nodes[0]-1;

          elem_con_size++;
          elem_con_ptr[num_elements+1] = elem_con_size;
          num_elements++;

          if (component_elems[9*(component_num-1)] == '\0'){
            strcpy(&component_elems[9*(component_num-1)], "FFORCE");
          }
        }
        else if (strncmp(line[0], "SPC", 3) == 0){
          // This is a variable-length format. Read in grid points until 
          // zero is reached. This is a fixed-width format
          // SPC SID  G1  C  D
          
          // Read in the nodal value
          char node[9];
          strncpy(node, &line[0][16], 8);
          node[8] = '\0';
          bc_nodes[num_bcs] = atoi(node)-1;
          
          strncpy(node, &line[0][32], 8);
          node[8] = '\0';
          double val = bdf_atof(node);
          
          // Read in the dof that will be constrained
          for ( int k = 24; k < 32; k++ ){
            char dofs[9] = "12345678";
            
            for ( int j = 0; j < 8; j++ ){
              if (dofs[j] == line[0][k]){
                bc_vars[bc_vars_size] = j;
                bc_vals[bc_vars_size] = val;
                bc_vars_size++;
                break;
              }
            }
          }
          
          bc_ptr[num_bcs+1] = bc_vars_size;
          num_bcs++;
        }
      }

      read_buffer_line(line[0], sizeof(line[0]), 
                       &buffer_loc, buffer, buffer_len);
    }

    delete [] buffer;

    if (fail){
      delete [] elem_nums;
      delete [] elem_comp;
      delete [] elem_con;
      delete [] elem_con_ptr;
      MPI_Abort(comm, fail);
      return fail;
    }

    // Arg sort the list of nodes
    int *node_args = new int[ num_nodes ]; 
    for ( int k = 0; k < num_nodes; k++ ){ 
      node_args[k] = k; 
    }

    arg_sort_list = node_nums;
    qsort(node_args, num_nodes, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Arg sort the list of elements
    int *elem_args = new int[ num_elements ];
    for ( int k = 0; k < num_elements; k++ ){
      elem_args[k] = k;
    }

    arg_sort_list = elem_nums;
    qsort(elem_args, num_elements, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Now node_nums[node_args[k]] and elem_nums[elem_args[k] are sorted.

    // Create the output for the nodes
    Xpts = new TacsScalar[3*num_nodes];

    for ( int k = 0; k < num_nodes; k++ ){
      int n = node_args[k];
      for ( int j = 0; j < 3; j++ ){
        Xpts[3*k+j] = Xpts_unsorted[3*n+j];
      }
    }

    // Read in the connectivity array and store the information
    elem_node_conn = new int[ elem_con_size ];
    elem_node_ptr = new int[ num_elements+1 ];
    elem_component = new int[ num_elements ];

    elem_node_ptr[0] = 0;
    for ( int k = 0, n = 0; k < num_elements; k++ ){
      int e = elem_args[k];

      for ( int j = elem_con_ptr[e]; j < elem_con_ptr[e+1]; j++, n++ ){
        int node_num = elem_con[j];

        // Find node_num in the list
        int node = find_index_arg_sorted(node_num, num_nodes,
                                         node_nums, node_args);
        if (node < 0){
          elem_node_conn[n] = -1;
          fail = 1;
        }
        else {
          elem_node_conn[n] = node;
        }
      }

      elem_component[k] = elem_comp[e];
      elem_node_ptr[k+1] = n;
    }

    // Find the boundary condition nodes
    for ( int k = 0; k < num_bcs; k++ ){
      int node = find_index_arg_sorted(bc_nodes[k], num_nodes,
                                       node_nums, node_args);    
      if (node < 0){
        fail = 1;
        bc_nodes[k] = -1;
      }
      else {
        bc_nodes[k] = node;
      }
    }

    // Free data that has been allocated locally
    delete [] elem_nums;
    delete [] elem_comp;
    delete [] elem_con;
    delete [] elem_con_ptr;
    delete [] elem_args;
    delete [] node_args;
  }

  // Distribute the component numbers and descritpions 
  MPI_Bcast(&num_components, 1, MPI_INT, root, comm);
  if (rank != root){
    component_elems = new char[9*num_components];
    component_descript = new char[33*num_components];
  }
  MPI_Bcast(component_elems, 9*num_components, MPI_CHAR, root, comm);
  MPI_Bcast(component_descript, 33*num_components, MPI_CHAR, root, comm);

  elements = new TACSElement*[ num_components ];
  for ( int k = 0; k < num_components; k++ ){
    elements[k] = NULL;
  }

  return fail;
}

/*
  Retrieve the number of nodes in the model
*/
int TACSMeshLoader::getNumNodes(){
  return num_nodes;
}

/*
  Create a TACSToFH5 file creation object
*/
TACSToFH5 *TACSMeshLoader::createTACSToFH5( TACSAssembler *tacs,
                                            ElementType elem_type,
                                            unsigned int write_flag ){
  // Set the component numbers in the elements
  for ( int k = 0; k < num_components; k++ ){
    elements[k]->setComponentNum(k);
  }

  TACSToFH5 * f5 = new TACSToFH5(tacs, elem_type, write_flag);
  for ( int k = 0; k < num_components; k++ ){
    if (strlen(&component_descript[33*k]) == 0){
      char name[64];
      sprintf(name, "Component %d", k);
      f5->setComponentName(k, name);
    }
    else {
      f5->setComponentName(k, &component_descript[33*k]);
    }
  }

  return f5;
}

/*
  Create a distributed version of TACS
*/
TACSAssembler *TACSMeshLoader::createTACS( int vars_per_node,
                                           TACSAssembler::OrderingType 
                                           order_type, 
                                           TACSAssembler::MatrixOrderingType 
                                           mat_type ){
  // Set the root processor
  const int root = 0;

  // Get the rank of the current processor
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Allocate the TACS creator
  creator = new TACSCreator(comm, vars_per_node);
  creator->incref();

  if (rank == root){
    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elements,
                                   elem_node_ptr, elem_node_conn,
                                   elem_component);
    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes, 
                                   bc_ptr, bc_vars, bc_vals);
    
    // Set the nodal locations
    creator->setNodes(Xpts);

    // Free things that are no longer required
    delete [] elem_node_ptr;   elem_node_ptr = NULL;
    delete [] elem_node_conn;  elem_node_conn = NULL;
    delete [] elem_component;  elem_component = NULL;

    // Free the boundary conditions
    delete [] bc_nodes;   bc_nodes = NULL;
    delete [] bc_ptr;     bc_ptr = NULL;
    delete [] bc_vars;    bc_vars = NULL;
    delete [] bc_vals;    bc_vals = NULL;
  }

  // This call must occur on all processor
  creator->setElements(elements, num_components);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();

  return tacs;
}

/*
  Retrieve the number of elements owned by this processes
*/
int TACSMeshLoader::getNumElements(){
  return num_elements;
}

/*
  Set the function domain

  Given the function, and the set of component numbers that define the
  domain of interest, set the element numbers in the function that
*/
void TACSMeshLoader::addFunctionDomain( TACSFunction * function,
                                        int comp_nums[], int num_comps ){
  if (creator){
    int *elems;
    int num_elems = creator->getElementIdNums(comp_nums, num_comps,
                                              &elems);
    function->addDomain(elems, num_elems);
    delete [] elems;
  }  
}

/*
  Add the auxiliary element to the given domain specified by the
  component number
*/
void TACSMeshLoader::addAuxElement( TACSAuxElements *aux, int component_num,
                                    TACSElement *element ){
  if (creator){
    int *elems;
    int num_elems = creator->getElementIdNums(&component_num, 1,
                                              &elems);
    for ( int i = 0; i < num_elems; i++ ){
      aux->addElement(elems[i], element);
    }
    delete [] elems;
  }
}
