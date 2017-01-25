#ifndef FH5_INCLUDE_H
#define FH5_INCLUDE_H

/*
  Copyright (c) 2012 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.

  Create a file that contains information about a finite-element
  problem. 
*/

#include "TACSObject.h"

class FH5File : public TACSObject {
 public:
  enum FH5DataNames { FH5_INT = 0,
                      FH5_DOUBLE };

  FH5File( MPI_Comm _comm );
  ~FH5File();

  // Create a file for output
  // ------------------------
  int createFile( const char *file_name,
                  char **component_names, 
                  int num_components );
  int writeZoneData( char *zone_name, 
                     enum FH5DataNames data_name, 
                     char *var_names,
                     void *data, int dim1, int dim2 );
  void close();

  // Open a file for reading input
  // -----------------------------
  int openFile( const char *file_name );

  // Retrieve the component names
  // ----------------------------
  int getNumComponents();
  char *getComponentName( int comp );

  // Retrieve zone data
  // ------------------
  void firstZone();
  int nextZone(); 
  int getZoneInfo( const char **zone_name,
                   const char **var_names,
                   int *dim1, int *dim2 );
  int getZoneData( const char **zone_name,
                   const char **var_names,
                   void **data, int *dim1, int *dim2 );
 private:
  // Store information about the location of the data within the file
  class FH5FileInfo {
  public:
    FH5FileInfo(){
      zone_name = NULL;
      next = NULL;
      var_names = NULL;
      dtype = -1;
      dim1 = dim2 = 0;
      data_offset = 0;
    }
    ~FH5FileInfo(){
      if (zone_name){ delete [] zone_name; }
      if (var_names){ delete [] var_names; }
    }
    int dtype;
    char *zone_name;
    char *var_names;
    int dim1, dim2;
    int data_offset;
    FH5FileInfo *next;
  } *root, *tip, *current;

  // Scan the file and record the header information
  void scanFH5File();
  void deleteFH5FileInfo();
  
  int num_comp; // The number of components
  char **comp_names; // The component names

  int file_for_writing; // Is this file for writing?
  MPI_Comm comm; // The communicator over which the 
  MPI_File fp; // The MPI file pointer
  MPI_Offset file_offset; // The offset into the file 
  MPI_Offset file_end; // The offset at the end of the file

  // Serial file containing the FE solution
  FILE *rfp;
};

#endif // FH5_INCLUDE_H
