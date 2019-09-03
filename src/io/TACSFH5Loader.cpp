/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSFH5Loader.h"

TACSFH5Loader::TACSFH5Loader(){
  comp_nums = NULL;
  ltypes = NULL;
  ptr = NULL;
  conn = NULL;

  num_elements = -1;
  conn_size = -1;

  num_nodes_continuous = -1;
  num_vals_continuous = -1;
  num_nodes_element = -1;
  num_vals_element = -1;

  data_file = NULL;
}

TACSFH5Loader::~TACSFH5Loader(){
  if (comp_nums){ delete [] comp_nums; }
  if (ltypes){ delete [] ltypes; }
  if (ptr){ delete [] ptr; }
  if (conn){ delete [] conn; }
  if (data_file){ data_file->decref(); }
}

int TACSFH5Loader::loadData( const char *conn_fname,
                             const char *data_fname ){
  // Load in the data for the connectivity
  TACSFH5File *conn_file = new TACSFH5File(MPI_COMM_SELF);
  conn_file->incref();

  // Open the connectivity file
  int fail = conn_file->openFile(conn_fname);
  if (fail){
    conn_file->decref();
    return fail;
  }

  conn_file->firstZone();
  int iterate = 1;
  while (iterate){
    const char *zone_name, *var_names;
    TACSFH5File::FH5DataType dtype;
    int dim1, dim2;

    // Get the name of the zone and its dimensions
    conn_file->getZoneInfo(&zone_name, &var_names, &dtype, &dim1, &dim2);

    if (strcmp("components", zone_name) == 0){
      void *idata;
      conn_file->getZoneData(NULL, NULL, NULL, NULL, NULL, &idata);
      comp_nums = (int*)idata;
      if (num_elements < 0){
        num_elements = dim1;
      }
      else if (num_elements != dim1){
        fprintf(stderr, "TACSFH5Loader: Number of elements is inconsistent\n");
        return 1;
      }
    }
    else if (strcmp("ltypes", zone_name) == 0){
      void *idata;
      conn_file->getZoneData(NULL, NULL, NULL, NULL, NULL, &idata);
      ltypes = (int*)idata;
      if (num_elements < 0){
        num_elements = dim1;
      }
      else if (num_elements != dim1){
        fprintf(stderr, "TACSFH5Loader: Number of elements is inconsistent\n");
        return 1;
      }
    }
    else if (strcmp("ptr", zone_name) == 0){
      void *idata;
      conn_file->getZoneData(NULL, NULL, NULL, NULL, NULL, &idata);
      ptr = (int*)idata;
      if (num_elements < 0){
        num_elements = dim1-1;
      }
      else if (num_elements != dim1-1){
        fprintf(stderr, "TACSFH5Loader: Number of elements is inconsistent\n");
        return 1;
      }
    }
    else if (strcmp("connectivity", zone_name) == 0){
      void *idata;
      conn_file->getZoneData(NULL, NULL, NULL, NULL, NULL, &idata);
      conn = (int*)idata;
      conn_size = dim1;
    }

    if (!conn_file->nextZone()){
      iterate = 0;
    }
  }

  conn_file->close();
  conn_file->decref();

  if (num_elements >= 0 && conn_size >= 0){
    if (!data_fname){
      data_fname = conn_fname;
    }

    // Load in the data for the connectivity
    data_file = new TACSFH5File(MPI_COMM_SELF);
    data_file->incref();

    // Open the connectivity file
    fail = data_file->openFile(data_fname);
    if (fail){
      data_file->decref();
    }

    iterate = 1;
    while (iterate){
      const char *zone_name, *var_names;
      TACSFH5File::FH5DataType dtype;
      int dim1, dim2;

      // Get the name of the zone and its dimensions
      data_file->getZoneInfo(&zone_name, &var_names, &dtype, &dim1, &dim2);

      if (strncmp("continuous data", zone_name, 15) == 0){
        num_nodes_continuous = dim1;
        num_vals_continuous = dim2;
      }
      else if (strncmp("element data", zone_name, 12) == 0){
        num_nodes_element = dim1;
        num_vals_element = dim2;
      }

      if (!data_file->nextZone()){
        iterate = 0;
      }
    }
  }

  return 0;
}

/**
   Get the number of components defined in this file
*/
int TACSFH5Loader::getNumComponents(){
  if (data_file){
    return data_file->getNumComponents();
  }
  return 0;
}

/**
   Return the component name

   @param comp The component number
   @return The component name
*/
char *TACSFH5Loader::getComponentName( int comp ){
  if (data_file){
    return data_file->getComponentName(comp);
  }
  return NULL;
}
