/*
  An FH5 to VTK converter. This only works for specially written FH5
  files. (This is designed primarily to work with TACS).

  Copyright (c) 2017 Graeme Kennedy. All rights reserved. 
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Include FH5 header files
#include "FH5.h"

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Test if the files exist/you have permission to modify them
  char *infile;
  char *outfile;

  // Convert hdf5 file argv[1] to tecplot file argv[2]
  if (argc == 1){
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }
  
  infile = new char[ strlen(argv[1])+1 ];
  strcpy(infile, argv[1]);

  // Set the output file
  outfile = new char[ strlen(infile)+5 ];
  int len = strlen(infile);
  int i = len-1;
  for ( ; i >= 0; i-- ){
    if (infile[i] == '.'){ break; }     
  }
  strcpy(outfile, infile);
  strcpy(&outfile[i], ".vtk");

  printf("Trying to convert FH5 file %s to vtk file %s\n", 
	 infile, outfile);

  // Open the output file
  FILE *fp = fopen(outfile, "w");
  if (!fp){
    fprintf(stderr, "Failed to open the output file %s\n", outfile);
    return (1);
  }

  // Open the FH5 file for reading
  FH5File * file = new FH5File(MPI_COMM_SELF);
  file->incref();

  if (!file->openFile(infile)){
    fprintf(stderr, "Failed to open the file %s\n", infile);
    return (1);
  }

  // Retrieve all the data from the file including the variables,
  // connectivity and component numbers
  double solution_time = 0.0;
  int *element_comp_num = NULL;
  int *conn = NULL;
  double *data = NULL;
  int conn_dim = 0, num_elements = 0;
  int num_points = 0, num_variables = 0;

  // Keep the variable names
  char *vars = NULL;

  // Load in the first zone
  file->firstZone();
  do {
    // Find the zone corresponding to all the data
    const char *zone_name, *var_names;
    int dim1, dim2;

    if (!file->getZoneInfo(&zone_name, &var_names, &dim1, &dim2)){
      fprintf(stderr, "Error, zone not defined\n");
      break;
    }
    
    if (strcmp(zone_name, "components") == 0){
      void *vdata;
      if (file->getZoneData(&zone_name, &var_names, &vdata, 
                            &dim1, &dim2)){
        element_comp_num = (int*)vdata;
      }
    }
    else if (strcmp(zone_name, "connectivity") == 0){
      num_elements = dim1;
      conn_dim = dim2;
      void *vdata;
      if (file->getZoneData(&zone_name, &var_names, &vdata, 
                            &dim1, &dim2)){
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
      if (file->getZoneData(&zone_name, &var_names, &vdata, 
                            &dim1, &dim2)){
        num_points = dim1;
        num_variables = dim2;
        data = (double*)vdata;
      }
    }
  } while (file->nextZone());
    
  if (!(element_comp_num && conn && data)){
    fprintf(stderr, 
            "Error, data, connectivity or component numbers not defined in file\n");
  }

  // Write out the vtk file
  fprintf(fp, "# vtk DataFile Version 3.0\n");
  fprintf(fp, "vtk output\nASCII\n");
  fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    
  // Write out the points
  fprintf(fp, "POINTS %d float\n", num_points);
  const double *d = data;
  for ( int k = 0; k < num_points; k++ ){
    fprintf(fp, "%e %e %e\n", d[0]+d[3], d[1]+d[4], d[2]+d[5]);
    d += num_variables;
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
    
      for ( int k = 0; k < num_points; k++ ){
        fprintf(fp, "%e\n", data[num_variables*k + j]);
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

  delete [] infile;
  delete [] outfile;

  MPI_Finalize();

  return (0);
}

