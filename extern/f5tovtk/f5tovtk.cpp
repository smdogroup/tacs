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

/*
  An FH5 to VTK converter. This only works for specially written FH5
  files. (This is designed primarily to work with TACS).
*/

// Include FH5 header files
#include "TACSFH5Loader.h"
#include "TACSElementTypes.h"

const int VTK_VERTEX = 1;
const int VTK_LINE = 3;
const int VTK_TRIANGLE = 5;
const int VTK_QUAD = 9;
const int VTK_TETRA = 10;
const int VTK_HEXAHEDRON = 12;

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Convert hdf5 file argv[1] to
  if (argc == 1){
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }

  for ( int iter = 1; iter < argc; iter++ ){
    char *infile = new char[ strlen(argv[iter])+1 ];
    strcpy(infile, argv[iter]);

    // Set the output file
    char *outfile = new char[ strlen(infile)+5 ];
    int len = strlen(infile);
    int i = len-1;
    for ( ; i >= 0; i-- ){
      if (infile[i] == '.'){ break; }     
    }
    strcpy(outfile, infile);
    strcpy(&outfile[i], ".vtk");

    printf("Trying to convert FH5 file %s to vtk file %s\n", 
           infile, outfile);

    // Create the loader object
    TACSFH5Loader *loader = new TACSFH5Loader();
    loader->incref();

    int fail = loader->loadData(infile);
    if (fail){
      fprintf(stderr, "Failed to open the file %s\n", infile);
      return (1);
    }

    int num_elements;
    int *comp_nums, *ltypes, *ptr, *conn;
    loader->getConnectivity(&num_elements, &comp_nums, &ltypes, &ptr, &conn);

    const char *cname, *cvars;
    int cdim1, cdim2;
    float *cdata;
    loader->getContinuousData(&cname, &cvars, &cdim1, &cdim2, &cdata);

    const char *ename, *evars;
    int edim1, edim2;
    float *edata;
    loader->getElementData(&ename, &evars, &edim1, &edim2, &edata);

    // Open the output file
    FILE *fp = fopen(outfile, "w");
    if (!fp){
      fprintf(stderr, "Failed to open the output file %s\n", outfile);
      return (1);
    }

    // Write out the vtk file
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d float\n", cdim1);

    const float *d = cdata;
    for ( int k = 0; k < cdim1; k++ ){
      fprintf(fp, "%e %e %e\n", d[0], d[1], d[2]);
      d += cdim2;
    }

    int num_basic_elements = 0;
    int basic_conn_size = 0;
    for ( int k = 0; k < num_elements; k++ ){
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
      num_basic_elements += ntypes;
      basic_conn_size += nconn;
    }

    int *basic_ltypes = new int[ num_basic_elements ];
    int *basic_conn = new int[ basic_conn_size ];

    int *btypes = basic_ltypes;
    int *bconn = basic_conn;
    for ( int k = 0; k < num_elements; k++ ){
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
      TacsConvertVisLayoutToBasic(ltype, &conn[ptr[k]],
                                  btypes, bconn);
      btypes += ntypes;
      bconn += nconn;
    }

    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", num_basic_elements,
            num_basic_elements + basic_conn_size);

    int basic_conn_offset = 0;
    for ( int k = 0; k < num_basic_elements; k++ ){
      ElementLayout ltype = (ElementLayout)basic_ltypes[k];
      int conn_size = TacsGetNumVisNodes(ltype);
      fprintf(fp, "%d ", conn_size);
      if (basic_ltypes[k] == TACS_QUAD_ELEMENT){
        const int convert[] = {0, 1, 3, 2};
        for ( int j = 0; j < conn_size; j++ ){
          fprintf(fp, "%d ", basic_conn[basic_conn_offset + convert[j]]);
        }
        basic_conn_offset += 4;
      }
      else if (basic_ltypes[k] == TACS_HEXA_ELEMENT){
        const int convert[] = {0, 1, 3, 2, 4, 5, 8, 7};
        for ( int j = 0; j < conn_size; j++ ){
          fprintf(fp, "%d ", basic_conn[basic_conn_offset + convert[j]]);
        }
        basic_conn_offset += 8;
      }
      else {
        for ( int j = 0; j < conn_size; j++, basic_conn_offset++ ){
          fprintf(fp, "%d ", basic_conn[basic_conn_offset]);
        }
      }
      fprintf(fp, "\n");
    }

    // All tetrahedrals...
    fprintf(fp, "\nCELL_TYPES %d\n", num_basic_elements);
    for ( int k = 0; k < num_basic_elements; k++ ){
      if (basic_ltypes[k] == TACS_POINT_ELEMENT){
        fprintf(fp, "%d\n", VTK_VERTEX);
      }
      else if (basic_ltypes[k] == TACS_LINE_ELEMENT){
        fprintf(fp, "%d\n", VTK_LINE);
      }
      else if (basic_ltypes[k] == TACS_TRI_ELEMENT){
        fprintf(fp, "%d\n", VTK_TRIANGLE);
      }
      else if (basic_ltypes[k] == TACS_QUAD_ELEMENT){
        fprintf(fp, "%d\n", VTK_QUAD);
      }
      else if (basic_ltypes[k] == TACS_TETRA_ELEMENT){
        fprintf(fp, "%d\n", VTK_TETRA);
      }
      else if (basic_ltypes[k] == TACS_HEXA_ELEMENT){
        fprintf(fp, "%d\n", VTK_HEXAHEDRON);
      }
    }
    delete [] basic_conn;
    delete [] basic_ltypes;

    // Print out the rest as fields one-by-one
    fprintf(fp, "POINT_DATA %d\n", cdim1);

    for ( int j = 0; j < cdim2; j++ ){
      char name[256];
      int index = 0;
      while (strlen(cvars) > 0 && cvars[0] != ','){
        name[index] = cvars[0];
        index++; cvars++;
      }
      name[index] = '\0';
      cvars++;

      // Write out the zone names
      if (j >= 3){
        fprintf(fp, "SCALARS %s float 1\n", name);
        fprintf(fp, "LOOKUP_TABLE default\n");

        for ( int k = 0; k < cdim1; k++ ){
          double d = cdata[cdim2*k + j];
          fprintf(fp, "%.3e\n", d);
        }
      }
    }

    // Count up the number of times each node is referred to
    // in the discontinuous element-wise data
    float *counts = new float[ cdim1 ];
    memset(counts, 0, cdim1*sizeof(float));
    for ( int j = 0; j < ptr[num_elements]; j++ ){
      counts[conn[j]] += 1.0;
    }
    for ( int i = 0; i < cdim1; i++ ){
      if (counts[i] != 0.0){
        counts[i] = 1.0/counts[i];
      }
    }

    // For each component, average the nodal data
    float *data = new float[ cdim1 ];
    for ( int j = 0; j < edim2; j++ ){
      char name[256];
      int index = 0;
      while (strlen(evars) > 0 && evars[0] != ','){
        name[index] = evars[0];
        index++; evars++;
      }
      name[index] = '\0';
      evars++;

      // Nodally average the data
      memset(data, 0, cdim1*sizeof(float));
      for ( int k = 0; k < ptr[num_elements]; k++ ){
        data[conn[k]] += counts[conn[k]]*edata[edim2*k + j];
      }

      // Write out the zone names
      fprintf(fp, "SCALARS %s float 1\n", name);
      fprintf(fp, "LOOKUP_TABLE default\n");

      for ( int k = 0; k < cdim1; k++ ){
        fprintf(fp, "%.3e\n", data[k]);
      }
    }

    delete [] counts;
    delete [] data;

    fclose(fp);

    /*

    // Retrieve all the data from the file including the variables,
    // connectivity and component numbers
    double solution_time = 0.0;
    int *element_comp_num = NULL;
    int *conn = NULL;
    double *data = NULL;
    float *float_data = NULL;
    int conn_dim = 0, num_elements = 0;
    int num_points = 0, num_variables = 0;

    // Keep the variable names
    char *vars = NULL;

    // Load in the first zone
    file->firstZone();
    do {
      // Find the zone corresponding to all the data
      const char *zone_name, *var_names;
      FH5File::FH5DataType dtype;
      int dim1, dim2;

      if (!file->getZoneInfo(&zone_name, &var_names, &dtype, &dim1, &dim2)){
        fprintf(stderr, "Error, zone not defined\n");
        break;
      }
      
      if (strcmp(zone_name, "components") == 0){
        void *vdata;
        if (file->getZoneData(&zone_name, &var_names, &dtype,
                              &vdata, &dim1, &dim2)){
          element_comp_num = (int*)vdata;
        }
      }
      else if (strcmp(zone_name, "connectivity") == 0){
        num_elements = dim1;
        conn_dim = dim2;
        void *vdata;
        if (file->getZoneData(&zone_name, &var_names, &dtype,
                              &vdata, &dim1, &dim2)){
          conn = (int*)vdata;
        }
      }
      else if (strncmp(zone_name, "data", 4) == 0){
        // Try to retrieve the solution time - this may fail if an older
        // version of the F5 file is used
        if (!(sscanf(zone_name, "data t=%lf", &solution_time) == 1)){
          solution_time = 0.0;
        }

        // Initialize the tecplot file with the variables
        vars = new char[ strlen(var_names)+1 ];
        strcpy(vars, var_names);
   
        // Retrieve the data
        void *vdata;
        if (file->getZoneData(&zone_name, &var_names, &dtype,
                              &vdata, &dim1, &dim2)){
          num_points = dim1;
          num_variables = dim2;
          if (dtype == FH5File::FH5_DOUBLE){
            data = (double*)vdata;
          }
          else if (dtype == FH5File::FH5_FLOAT){
            float_data = (float*)vdata;
          }
        }
      }
    } while (file->nextZone());
      
    if (!(element_comp_num && conn && (data || float_data))){
      fprintf(stderr, 
              "Error, data, connectivity or component numbers not defined in file\n");
    }

    // Write out the vtk file
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
      
    // Write out the points
    fprintf(fp, "POINTS %d float\n", num_points);

    if (data){
      const double *d = data;
      if (plot_displaced_shape){
        for ( int k = 0; k < num_points; k++ ){
          fprintf(fp, "%e %e %e\n", d[0]+d[3], d[1]+d[4], d[2]+d[5]);
          d += num_variables;
        }
      }
      else {
        for ( int k = 0; k < num_points; k++ ){
          fprintf(fp, "%e %e %e\n", d[0], d[1], d[2]);
          d += num_variables;
        }
      }
    }
    else if (float_data){
      const float *d = float_data;
      if (plot_displaced_shape){
        for ( int k = 0; k < num_points; k++ ){
          fprintf(fp, "%e %e %e\n", d[0]+d[3], d[1]+d[4], d[2]+d[5]);
          d += num_variables;
        }
      }
      else {
        for ( int k = 0; k < num_points; k++ ){
          fprintf(fp, "%e %e %e\n", d[0], d[1], d[2]);
          d += num_variables;
        }
      }
    }

    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", 
            num_elements, (conn_dim+1)*num_elements);
    const int *c = conn;
    for ( int k = 0; k < num_elements; k++ ){
      fprintf(fp, "%d ", conn_dim);
      for ( int j = 0; j < conn_dim; j++ ){
        fprintf(fp, "%d ", c[0]);
        c++;
      }
      fprintf(fp, "\n");
    }

    int vtk_elem_id = 0;
    if (conn_dim == 2){ 
      vtk_elem_id = 3; // VTK_LINE;
    }
    else if (conn_dim == 8){ 
      vtk_elem_id = 12; // VTK_HEXAHEDRON
    }
    else {
      vtk_elem_id = 9; // VTK_QUAD
    }

    // All quadrilaterals
    fprintf(fp, "\nCELL_TYPES %d\n", num_elements);
    for ( int k = 0; k < num_elements; k++ ){
      fprintf(fp, "%d\n", vtk_elem_id);
    }

    // Print out the rest as fields one-by-one
    fprintf(fp, "POINT_DATA %d\n", num_points);

    const char *ptr = vars;
    for ( int j = 0; j < num_variables; j++ ){
      char name[256];
      int index = 0;
      while (strlen(ptr) > 0 && ptr[0] != ','){ 
        name[index] = ptr[0];
        index++; ptr++;
      }
      name[index] = '\0';
      ptr++;

      // Write out the zone names
      if (j >= 3){
        fprintf(fp, "SCALARS %s float 1\n", name);
        fprintf(fp, "LOOKUP_TABLE default\n");
      
        if (data){
          for ( int k = 0; k < num_points; k++ ){
            fprintf(fp, "%.3e\n", data[num_variables*k + j]);
          }
        }
        else if (float_data){
          for ( int k = 0; k < num_points; k++ ){
            double d = float_data[num_variables*k + j];
            fprintf(fp, "%.3e\n", d);
          }          
        }
      }
    }
    fclose(fp);

    // Free the variable names
    delete [] vars;

    file->close();
    file->decref();

    // Clean up memory
    delete [] data;
    delete [] conn;
    delete [] element_comp_num;
    */

    loader->decref();

    delete [] infile;
    delete [] outfile;
  }

  MPI_Finalize();

  return (0);
}

