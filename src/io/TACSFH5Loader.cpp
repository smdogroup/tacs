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
  data_file = NULL;

  comp_nums = NULL;
  ltypes = NULL;
  ptr = NULL;
  conn = NULL;

  num_elements = -1;
  conn_size = -1;

  continuous_data = NULL;
  continuous_zone = NULL;
  continuous_vars = NULL;
  num_nodes_continuous = -1;
  num_vals_continuous = -1;

  element_data = NULL;
  element_zone = NULL;
  element_vars = NULL;
  num_nodes_element = -1;
  num_vals_element = -1;
}

TACSFH5Loader::~TACSFH5Loader(){
  if (comp_nums){ delete [] comp_nums; }
  if (ltypes){ delete [] ltypes; }
  if (ptr){ delete [] ptr; }
  if (conn){ delete [] conn; }
  if (continuous_data){ delete [] continuous_data; }
  if (element_data){ delete [] element_data; }
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
    conn_file = NULL;
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
      data_file = NULL;
      return fail;
    }

    iterate = 1;
    while (iterate){
      const char *zone_name, *var_names;
      TACSFH5File::FH5DataType dtype;
      int dim1, dim2;

      // Get the name of the zone and its dimensions
      data_file->getZoneInfo(&zone_name, &var_names, &dtype, &dim1, &dim2);

      if (strncmp("continuous data", zone_name, 15) == 0){
        void *fdata;
        data_file->getZoneData(&continuous_zone, &continuous_vars,
                               NULL, NULL, NULL, &fdata);
        num_nodes_continuous = dim1;
        num_vals_continuous = dim2;
        continuous_data = (float*)fdata;
      }
      else if (strncmp("element data", zone_name, 12) == 0){
        void *fdata;
        data_file->getZoneData(&element_zone, &element_vars,
                               NULL, NULL, NULL, &fdata);
        num_nodes_element = dim1;
        num_vals_element = dim2;
        element_data = (float*)fdata;
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

void TACSFH5Loader::getConnectivity( int *_num_elements,
                                     int **_comp_nums, int **_ltypes,
                                     int **_ptr, int **_conn ){
  if (_num_elements){ *_num_elements = num_elements; }
  if (_comp_nums){ *_comp_nums = comp_nums; }
  if (_ltypes){ *_ltypes = ltypes; }
  if (_ptr){ *_ptr = ptr; }
  if (_conn){ *_conn = conn; }
}

void TACSFH5Loader::getContinuousData( const char **zone_name,
                                       const char **var_names,
                                       int *dim1, int *dim2, float **data ){
  if (zone_name){ *zone_name = continuous_zone; }
  if (var_names){ *var_names = continuous_vars; }
  if (dim1){ *dim1 = num_nodes_continuous; }
  if (dim2){ *dim2 = num_vals_continuous; }
  if (data){ *data = continuous_data; }
}

void TACSFH5Loader::getElementData( const char **zone_name,
                                    const char **var_names,
                                    int *dim1, int *dim2, float **data ){
  if (zone_name){ *zone_name = element_zone; }
  if (var_names){ *var_names = element_vars; }
  if (dim1){ *dim1 = num_nodes_element; }
  if (dim2){ *dim2 = num_vals_element; }
  if (data){ *data = element_data; }
}
