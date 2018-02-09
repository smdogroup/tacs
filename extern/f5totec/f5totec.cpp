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
  An FH5 to tecplot converter. This only works for specially written
  FH5 files. (This is designed primarily to work with TACS).
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Include FH5 header files
#include "FH5.h"

// Include Tecplot header files
#include "TECIO.h"

enum FileType { FULL=0, GRID=1, SOLUTION=2 };

enum ZoneType { ORDERED=0, FELINESEG, FETRIANGLE, 
                FEQUADRILATERAL, FETETRAHEDRON, FEBRICK, 
                FEPOLYGON, FEPOLYHEDRA };

/*
  Initialize a data file.

  data_info == A description of the data
  var_names == Comma separated variable names
  file_name == Name of the file
  dir_name  == Name of the directory
  file_type == Type of data
*/
int create_tec_file( char *data_info, char *var_names,
                     char *file_name, char *dir_name,
                     enum FileType _file_type ){
  INTEGER4 file_type = _file_type;
  INTEGER4 debug = 0;
  INTEGER4 variables_are_double = 0;
  return TECINI112(data_info, var_names, file_name, dir_name,
                   &file_type, &debug, &variables_are_double);
}

/*
  A method to create a zone without the extra stuff
  
  zone_name    == The name of the zone to use
  zone_type    == One of LINESEG, TRIANGLE, QUAD, BRICK etc
  num_points   == The number of points
  num_elements == The number of elements
*/
int create_fe_tec_zone( char *zone_name, ZoneType _zone_type,
                        int _num_points, int _num_elements,
                        int use_strands=0,
                        double solution_time=0.0 ){
  if (_zone_type == ORDERED ||
      _zone_type == FEPOLYGON ||
      _zone_type == FEPOLYHEDRA ){
    fprintf(stderr, "Cannot create finite element zone with given \
zone type\n");
    return -1;
  }

  INTEGER4 zone_type = _zone_type;
  INTEGER4 num_points = _num_points;
  INTEGER4 num_elements = _num_elements;
  INTEGER4 num_faces = 0; // For all zones allowed here
  INTEGER4 icmax = 0, jcmax = 0, kcmax = 0; // Ignored
  INTEGER4 strand_id = 0;
  INTEGER4 parent_zone = 0;
  INTEGER4 is_block = 1; // Apparently this always needs to be 1
  // These are only for cell-based finite element data - we use node-based
  INTEGER4 num_face_connections = 0;
  INTEGER4 face_neighbour_mode = 0;
  INTEGER4 total_num_face_nodes = 0;
  INTEGER4 num_connected_boundary_faces = 0;
  INTEGER4 total_num_boundary_connections = 0;
  INTEGER4 *passive_var_list = NULL; // No passive variables
  INTEGER4 *value_location = NULL; // All values are nodal values
  INTEGER4 *share_var_from_zone = NULL;
  INTEGER4 share_con_from_zone = 0;

  // If we're using strands, set the strand ID
  if (use_strands){
    strand_id = 1;
  }

  return TECZNE112(zone_name, &zone_type, &num_points, &num_elements, 
                   &num_faces, &icmax, &jcmax, &kcmax,
                   &solution_time, &strand_id, &parent_zone, &is_block,
                   &num_face_connections, &face_neighbour_mode,
                   &total_num_face_nodes, &num_connected_boundary_faces,
                   &total_num_boundary_connections, 
                   passive_var_list, value_location,
                   share_var_from_zone, &share_con_from_zone);
}

/*
  Write data to a tecplot file

  len  == Length of the data
  data == The array of data
*/
int write_tec_double_data( int _len, double *data ){
  INTEGER4 len = _len;
  INTEGER4 is_double = 1;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write float data to a tecplot file
*/
int write_tec_float_data( int _len, float *data ){
  INTEGER4 len = _len;
  INTEGER4 is_double = 0;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write the connectivity data
*/
int write_con_data( int * con_data ){
  return TECNOD112(con_data);
}

/*
  End the file output
*/
int close_tec_file(){
  return TECEND112();
}

int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);

  // Check if we're going to use strands or not
  int use_strands = 0;
  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "--use_strands") == 0){
      use_strands = 1;
    }
  }

  // Convert hdf5 file argv[1] to tecplot file argv[2]
  if (argc == 1){
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }
  
  for ( int k = 1; k < argc; k++ ){
    char *infile = new char[ strlen(argv[k])+1 ];
    strcpy(infile, argv[k]);

    FILE *fp = fopen(infile, "r");
    if (fp){
      fclose(fp);
    }
    else {
      delete [] infile;
      continue;
    }

    // Set the output file
    char *outfile = new char[ strlen(infile)+5 ];
    int len = strlen(infile);
    int i = len-1;
    for ( ; i >= 0; i-- ){
      if (infile[i] == '.'){ break; }     
    }
    strcpy(outfile, infile);
    strcpy(&outfile[i], ".plt");
    
    printf("Trying to convert FH5 file %s to tecplot file %s\n", 
           infile, outfile);

    char data_info[] = "Created by f5totec";
    char dir_name[] = "."; // We'll be working with the current directory
    int tec_init = 0;      // Tecplot initialized flag

    // Open the FH5 file for reading
    FH5File * file = new FH5File(MPI_COMM_SELF);
    file->incref();
    
    if (!file->openFile(infile)){
      fprintf(stderr, "Failed to open the file %s\n", infile);
      return (1);
    }

    // Retrieve all the data from the file including the 
    // variables, connectivity and component numbers
    double solution_time = 0.0;
    int *element_comp_num = NULL;
    int *conn = NULL;
    double *data = NULL;
    float *float_data = NULL;
    int conn_dim = 0, num_elements = 0, num_points = 0, num_variables = 0;
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
        char *vars = new char[ strlen(var_names)+1 ];
        strcpy(vars, var_names);
        create_tec_file(data_info, vars,
                        outfile, dir_name, FULL);
        tec_init = 1;
        delete [] vars;
 
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
      fprintf(stderr, "Error, data, connectivity or \
component numbers not defined in file\n");
    }

    // Set the element type to use
    ZoneType zone_type;
    if (conn_dim == 2){      
      zone_type = FELINESEG; 
    }
    else if (conn_dim == 8){ 
      zone_type = FEBRICK; 
    }
    else {
      zone_type = FEQUADRILATERAL; 
    }

    int num_comp = file->getNumComponents();

    int *reduced_points = new int[ num_points ];
    int *reduced_conn = new int[ conn_dim*num_elements ];
    double *reduced_data = NULL;
    float *reduced_float_data = NULL;
    if (data){ 
      reduced_data = new double[ num_points ];
    }
    else if (float_data){
      reduced_float_data = new float[ num_points ];      
    }

    for ( int k = 0; k < num_comp; k++ ){
      // Count up the number of elements that use the connectivity
      char *comp_name = file->getComponentName(k);
      //printf("Converting zone %d: %s at time %g\n", 
      // k, comp_name, solution_time);

      memset(reduced_points, 0, num_points*sizeof(int));
      memset(reduced_conn, 0, conn_dim*num_elements*sizeof(int));

      int npts = 1, nelems = 0;
      // Count up the number of points/elements in this sub-domain
      for ( int i = 0; i < num_elements; i++ ){
        if (element_comp_num[i] == k){
          // Add this element to the reduced connectivity
          for ( int j = 0; j < conn_dim; j++ ){
            int pt = conn[conn_dim*i + j];
          
            // If a reduced numbering has not been applied to this point,
            // create a new number for it
            if (reduced_points[pt] == 0){
              reduced_points[pt] = npts;
              npts++;
            }
          
            // Set the reduced connectivity
            reduced_conn[conn_dim*nelems + j] = reduced_points[pt];
          }
          nelems++;
        }
      }

      // Since we started at npts = 1, we have one more point
      // than the actual number of points.
      npts--;

      if (nelems > 0 && npts > 0){
        // Create the zone with the solution time
        create_fe_tec_zone(comp_name, zone_type, npts, nelems,
                           use_strands, solution_time);

        // Retrieve the data
        for ( int j = 0; j < num_variables; j++ ){
          if (reduced_data){
            for ( int i = 0; i < num_points; i++ ){
              if (reduced_points[i] > 0){
                reduced_data[reduced_points[i]-1] = data[i*num_variables + j];
              }
            }
            write_tec_double_data(npts, reduced_data);
          }
          else if (reduced_float_data){
            for ( int i = 0; i < num_points; i++ ){
              if (reduced_points[i] > 0){
                reduced_float_data[reduced_points[i]-1] = 
                  float_data[i*num_variables + j];
              }
            }
            write_tec_float_data(npts, reduced_float_data);
          }
        }
      
        // Now, write the connectivity
        write_con_data(reduced_conn);
      }
    }

    if (tec_init){
      close_tec_file();
    }
    file->close();
    file->decref();

    // Clean up memory
    delete [] reduced_points;
    delete [] reduced_conn;
    delete [] reduced_data;

    delete [] data;
    delete [] conn;
    delete [] element_comp_num;

    delete [] infile;
    delete [] outfile;
  }

  MPI_Finalize();

  return (0);
}

