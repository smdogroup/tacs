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

#ifndef FH5_INCLUDE_H
#define FH5_INCLUDE_H

#include "TACSObject.h"

/*
  Create a file that contains information about a finite-element
  problem.
*/

class TACSFH5File : public TACSObject {
 public:
  // Data types accepted by FH5: Note that float comes last
  // for backwards compatibility
  enum FH5DataType { FH5_INT = 0, FH5_DOUBLE = 1, FH5_FLOAT = 2 };

  // Create the FH5 object
  TACSFH5File(MPI_Comm _comm);
  ~TACSFH5File();

  // Create an output file
  int createFile(const char *file_name, int num_components,
                 char **component_names);
  int writeZoneData(char *zone_name, char *var_names, FH5DataType data_name,
                    int dim1, int dim2, void *data, int *dim1_range = NULL);
  void close();

  // Open a file for reading input
  int openFile(const char *file_name);

  // Retrieve the component names
  int getNumComponents();
  char *getComponentName(int comp);

  // Retrieve zone data
  void firstZone();
  int nextZone();
  int getZoneInfo(const char **zone_name, const char **var_names,
                  FH5DataType *_dtype, int *dim1, int *dim2);
  int getZoneData(const char **zone_name, const char **var_names,
                  FH5DataType *_dtype, int *dim1, int *dim2, void **data);

 private:
  // Store information about the location of the data within the file
  class FH5FileInfo {
   public:
    FH5FileInfo() {
      zone_name = NULL;
      next = NULL;
      var_names = NULL;
      dtype = -1;
      dim1 = dim2 = 0;
      data_offset = 0;
    }
    ~FH5FileInfo() {
      if (zone_name) {
        delete[] zone_name;
      }
      if (var_names) {
        delete[] var_names;
      }
    }
    int dtype;
    char *zone_name;
    char *var_names;
    int dim1, dim2;
    size_t data_offset;
    FH5FileInfo *next;
  } * root, *tip, *current;

  // Scan the file and record the header information
  int scanFH5File();
  void deleteFH5FileInfo();

  int num_comp;       // The number of components
  char **comp_names;  // The component names

  int file_for_writing;    // Is this file for writing?
  MPI_Comm comm;           // The communicator over which the
  MPI_File fp;             // The MPI file pointer
  MPI_Offset file_offset;  // The offset into the file
  MPI_Offset file_end;     // The offset at the end of the file

  // Serial file containing the FE solution
  FILE *rfp;
};

#endif  // FH5_INCLUDE_H
