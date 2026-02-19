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

#ifndef TACS_FH5_LOADER_H
#define TACS_FH5_LOADER_H

#include "TACSElementTypes.h"
#include "TACSFH5.h"

/**
   Load data from a file created by the TACSToFH5.

   This code reads in the connectivity and element-wise data created
   via TACSToFH5 into memory. The data can then be accessed via member
   functions. You can copy out the data if you so desire, but only one
   copy of the data is ever stored by TACSFH5Loader.
*/
class TACSFH5Loader : public TACSObject {
 public:
  TACSFH5Loader();
  ~TACSFH5Loader();

  // Load data from the file
  int loadData(const char *conn_file, const char *data_file = NULL);

  // Get the component names/data from the file
  int getNumComponents();
  char *getComponentName(int comp);

  // Get the component data
  void getConnectivity(int *num_elements, int **comp_nums, int **ltypes,
                       int **ptr, int **conn);
  void getContinuousData(const char **zone_name, const char **var_names,
                         int *dim1, int *dim2, float **data);
  void getElementData(const char **zone_name, const char **var_names, int *dim1,
                      int *dim2, float **data);

  // Methods for post-processing data
  void getElementDataAsContinuous(int index, float *data);
  void computeValueMask(ElementLayout layout, int use_continuous_data,
                        int index, float lower, float upper, int *mask);
  void computePlanarMask(ElementLayout layout, const float base[],
                         const float normal[], int *mask);
  void getUnmatchedEdgesAndFaces(ElementLayout layout, const int *mask,
                                 int *_num_edges, int **_edges, int *_num_faces,
                                 int **_faces);
  void getIsoSurfaces(ElementLayout layout, const int *mask, float isoval,
                      int index, float *_data, int *_ntris, float **_verts);
  void getUnmatchedEdgesAndFaces(ElementLayout layout, const int *mask,
                                 int index, const float *data, int *_num_points,
                                 float **_points, float **_values,
                                 int *_num_edges, int **_edges, int *_ntris,
                                 int **_verts);

 private:
  void computeNodeToElement(ElementLayout layout, const int *mask,
                            int num_sub_indices, const int *sub_indices,
                            int **_node_to_element_ptr, int **_node_to_element);

  // Things associated with the types of elements
  int num_elements;
  int *comp_nums, *ltypes;

  // Data entries associated with the connectivity
  int conn_size;
  int *ptr, *conn;

  // Continuous data
  const char *continuous_zone, *continuous_vars;
  int num_nodes_continuous, num_vals_continuous;
  float *continuous_data;

  // Element data
  const char *element_zone, *element_vars;
  int num_nodes_element, num_vals_element;
  float *element_data;

  // Open file that contains the
  TACSFH5File *data_file;
};

#endif  // TACS_FH5_LOADER_H
