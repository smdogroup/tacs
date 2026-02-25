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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Include FH5 header files
#include "TACSFH5Loader.h"
// Include TACS element types
#include "TACSElementTypes.h"

// Include Tecplot header files
#include "TECIO.h"

enum FileType { FULL = 0, GRID = 1, SOLUTION = 2 };

enum ZoneType {
  ORDERED = 0,
  FELINESEG,
  FETRIANGLE,
  FEQUADRILATERAL,
  FETETRAHEDRON,
  FEBRICK,
  FEPOLYGON,
  FEPOLYHEDRA
};

/*
  Initialize a data file.

  data_info == A description of the data
  var_names == Comma separated variable names
  file_name == Name of the file
  dir_name  == Name of the directory
  file_type == Type of data
*/
int create_tec_file(char *data_info, char *var_names, char *file_name,
                    char *dir_name, enum FileType _file_type) {
  INTEGER4 file_type = _file_type;
  INTEGER4 debug = 0;
  INTEGER4 variables_are_double = 0;
  return TECINI112(data_info, var_names, file_name, dir_name, &file_type,
                   &debug, &variables_are_double);
}

/*
  A method to create a zone without the extra stuff

  zone_name    == The name of the zone to use
  zone_type    == One of LINESEG, TRIANGLE, QUAD, BRICK etc
  num_points       == The number of points
  num_elements     == The number of elements
  value_location   == The location of the values (0 for element-wise, 1 for
  nodal)
*/
int create_fe_tec_zone(char *zone_name, ZoneType _zone_type, int _num_points,
                       int _num_elements, int *_value_location = NULL,
                       int use_strands = 0, double solution_time = 0.0) {
  if (_zone_type == ORDERED || _zone_type == FEPOLYGON ||
      _zone_type == FEPOLYHEDRA) {
    fprintf(stderr,
            "Cannot create finite element zone with given \
zone type\n");
    return -1;
  }

  INTEGER4 zone_type = _zone_type;
  INTEGER4 num_points = _num_points;
  INTEGER4 num_elements = _num_elements;
  INTEGER4 num_faces = 0;                    // For all zones allowed here
  INTEGER4 icmax = 0, jcmax = 0, kcmax = 0;  // Ignored
  INTEGER4 strand_id = 0;
  INTEGER4 parent_zone = 0;
  INTEGER4 is_block = 1;  // Apparently this always needs to be 1
  // These are only for cell-based finite element data - we use node-based
  INTEGER4 num_face_connections = 0;
  INTEGER4 face_neighbour_mode = 0;
  INTEGER4 total_num_face_nodes = 0;
  INTEGER4 num_connected_boundary_faces = 0;
  INTEGER4 total_num_boundary_connections = 0;
  INTEGER4 *passive_var_list = NULL;  // No passive variables
  INTEGER4 *value_location =
      _value_location;  // Array defining if value is nodal or element-wise
  INTEGER4 *share_var_from_zone = NULL;
  INTEGER4 share_con_from_zone = 0;

  // If we're using strands, set the strand ID
  if (use_strands) {
    strand_id = 1;
  }

  return TECZNE112(zone_name, &zone_type, &num_points, &num_elements,
                   &num_faces, &icmax, &jcmax, &kcmax, &solution_time,
                   &strand_id, &parent_zone, &is_block, &num_face_connections,
                   &face_neighbour_mode, &total_num_face_nodes,
                   &num_connected_boundary_faces,
                   &total_num_boundary_connections, passive_var_list,
                   value_location, share_var_from_zone, &share_con_from_zone);
}

/*
  Write data to a tecplot file

  len  == Length of the data
  data == The array of data
*/
int write_tec_double_data(int _len, double *data) {
  INTEGER4 len = _len;
  INTEGER4 is_double = 1;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write float data to a tecplot file
*/
int write_tec_float_data(int _len, float *data) {
  INTEGER4 len = _len;
  INTEGER4 is_double = 0;
  return TECDAT112(&len, data, &is_double);
}

/*
  Write the connectivity data
*/
int write_con_data(int *con_data) { return TECNOD112(con_data); }

/*
  End the file output
*/
int close_tec_file() { return TECEND112(); }

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Check if we're going to use strands or not
  int use_strands = 0;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "--use_strands") == 0) {
      use_strands = 1;
    }
  }

  // Convert hdf5 file argv[1] to tecplot file argv[2]
  if (argc == 1) {
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }

  for (int k = 1; k < argc; k++) {
    char *infile = new char[strlen(argv[k]) + 1];
    strcpy(infile, argv[k]);

    FILE *fp = fopen(infile, "r");
    if (fp) {
      fclose(fp);
    } else {
      delete[] infile;
      continue;
    }

    // Set the output file
    char *outfile = new char[strlen(infile) + 5];
    int len = strlen(infile);
    int i = len - 1;
    for (; i >= 0; i--) {
      if (infile[i] == '.') {
        break;
      }
    }
    strcpy(outfile, infile);
    strcpy(&outfile[i], ".plt");

    printf("Trying to convert FH5 file %s to tecplot file %s\n", infile,
           outfile);

    char data_info[] = "Created by f5totec";
    char dir_name[] = ".";  // We'll be working with the current directory
    int tec_init = 0;       // Tecplot initialized flag

    // Create the loader object
    TACSFH5Loader *loader = new TACSFH5Loader();
    loader->incref();

    int fail = loader->loadData(infile);
    if (fail) {
      fprintf(stderr, "Failed to open the file %s\n", infile);
      return (1);
    }

    int num_elements;
    int *element_comp_num, *ltypes, *ptr, *conn;
    loader->getConnectivity(&num_elements, &element_comp_num, &ltypes, &ptr,
                            &conn);

    const char *cname, *var_names;
    int num_points, num_variables;
    float *cdata;
    loader->getContinuousData(&cname, &var_names, &num_points, &num_variables,
                              &cdata);

    const char *ename, *evar_names;
    int edim1, num_evariables;
    float *edata;
    loader->getElementData(&ename, &evar_names, &edim1, &num_evariables,
                           &edata);

    double solution_time = 0.0;

    // Initialize the tecplot file with the variables
    // Concatenate continuous and element variable names
    char *all_vars = new char[strlen(var_names) + strlen(evar_names) + 2];
    strcpy(all_vars, var_names);
    all_vars[strlen(var_names)] = ',';
    strcpy(&all_vars[strlen(var_names) + 1], evar_names);
    create_tec_file(data_info, all_vars, outfile, dir_name, FULL);
    tec_init = 1;
    delete[] all_vars;

    // For each element, average the values at the nodes to get a single element
    // value
    float *avg_edata = new float[num_elements * num_evariables];
    memset(avg_edata, 0, num_elements * num_evariables * sizeof(float));
    for (int i = 0; i < num_elements; i++) {
      int nnodes = ptr[i + 1] - ptr[i];
      for (int j = 0; j < num_evariables; j++) {
        for (int k = ptr[i]; k < ptr[i + 1]; k++) {
          avg_edata[num_evariables * i + j] += edata[num_evariables * k + j];
        }
        avg_edata[num_evariables * i + j] /= nnodes;
      }
    }

    if (!(element_comp_num && conn && cdata)) {
      fprintf(stderr,
              "Error, data, connectivity or \
              component numbers not defined in file\n");
    }

    // Setup visualization elements
    int num_basic_elements = 0;
    int basic_conn_size = 0;
    for (int k = 0; k < num_elements; k++) {
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
      // For triangular elements we'll add a 4th dummy node
      // to the connectivity that's just the third node repeated
      // This way triangles can be treated as degenerate quads from Tecplot's
      // perspective
      if (ltype == TACS_TRI_ELEMENT || ltype == TACS_TRI_QUADRATIC_ELEMENT ||
          ltype == TACS_TRI_CUBIC_ELEMENT) {
        nconn = nconn + ntypes;
      }
      // Plot higher order tetrahedral elements as a single linear element
      else if (ltype == TACS_TETRA_QUADRATIC_ELEMENT ||
               ltype == TACS_TETRA_CUBIC_ELEMENT) {
        nconn = 4;
        ntypes = 1;
      }
      num_basic_elements += ntypes;
      basic_conn_size += nconn;
    }

    int *basic_ltypes = new int[num_basic_elements];
    int *basic_element_comp_num = new int[num_basic_elements];
    int *basic_conn = new int[basic_conn_size];
    int *basic_element_global_ptr = new int[num_basic_elements];

    int *btypes = basic_ltypes;
    int *belem_comp_num = basic_element_comp_num;
    int *bconn = basic_conn;
    int *belem_global_ptr = basic_element_global_ptr;
    for (int k = 0; k < num_elements; k++) {
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      // Add our dummy nodes for triangular elements
      if (ltype == TACS_TRI_ELEMENT || ltype == TACS_TRI_QUADRATIC_ELEMENT ||
          ltype == TACS_TRI_CUBIC_ELEMENT) {
        TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
        int *tri_conn = new int[nconn];
        TacsConvertVisLayoutToBasic(ltype, &conn[ptr[k]], btypes, tri_conn);

        const int convert[] = {0, 1, 2, 2};
        for (int jj = 0; jj < ntypes; jj++) {
          for (int ii = 0; ii < 4; ii++) {
            bconn[4 * jj + ii] = tri_conn[3 * jj + convert[ii]];
          }
          btypes[jj] = TACS_QUAD_ELEMENT;
        }
        nconn = nconn + ntypes;
        delete[] tri_conn;
      }
      // Plot only the first four nodes (conrners) of higher order tets
      else if (ltype == TACS_TETRA_QUADRATIC_ELEMENT ||
               ltype == TACS_TETRA_CUBIC_ELEMENT) {
        memcpy(bconn, &conn[ptr[k]], 4 * sizeof(int));
        nconn = 4;
        ntypes = 1;
        btypes[0] = TACS_TETRA_ELEMENT;
      } else {
        TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
        TacsConvertVisLayoutToBasic(ltype, &conn[ptr[k]], btypes, bconn);
      }
      // Set the basic element component to match the parent
      for (int ii = 0; ii < ntypes; ii++) {
        belem_comp_num[ii] = element_comp_num[k];
      }
      btypes += ntypes;
      belem_comp_num += ntypes;
      bconn += nconn;
      for (int ii = 0; ii < ntypes; ii++) {
        belem_global_ptr[ii] = k;
      }
      belem_global_ptr += ntypes;
    }

    int num_comp = loader->getNumComponents();

    int *reduced_points = new int[num_points];
    int *reduced_conn = new int[basic_conn_size];
    int *value_location = new int[num_variables + num_evariables];
    for (int j = 0; j < num_variables; j++) {
      value_location[j] = 1;  // 1 = nodal
    }
    for (int j = 0; j < num_evariables; j++) {
      value_location[num_variables + j] = 0;  // 0 = cell-centered
    }

    for (int k = 0; k < num_comp; k++) {
      // Count up the number of elements that use the connectivity
      char *comp_name = loader->getComponentName(k);
      // printf("Converting zone %d: %s at time %g\n",
      //  k, comp_name, solution_time);

      memset(reduced_points, 0, num_points * sizeof(int));
      memset(reduced_conn, 0, basic_conn_size * sizeof(int));

      int npts = 1, nelems = 0;
      int zone_btype = -1;
      int basic_conn_offset = 0;
      // Count up the number of points/elements in this sub-domain
      for (int i = 0; i < num_basic_elements; i++) {
        ElementLayout ltype = (ElementLayout)basic_ltypes[i];
        int conn_size = TacsGetNumVisNodes(ltype);

        if (basic_element_comp_num[i] == k) {
          // Make sure all elements in this zone are the same type
          if (zone_btype == -1) {
            zone_btype = basic_ltypes[i];
          } else if (zone_btype != basic_ltypes[i]) {
            fprintf(stderr, "Component %d has conflicting element types\n", k);
            return (1);
          }

          int pt;
          for (int j = 0; j < conn_size; j++) {
            // Add this element to the reduced connectivity
            if (basic_ltypes[i] == TACS_QUAD_ELEMENT) {
              const int convert[] = {0, 1, 3, 2};
              pt = basic_conn[basic_conn_offset + convert[j]];
            } else if (basic_ltypes[i] == TACS_HEXA_ELEMENT) {
              const int convert[] = {0, 1, 3, 2, 4, 5, 7, 6};
              pt = basic_conn[basic_conn_offset + convert[j]];
            } else {
              pt = basic_conn[basic_conn_offset + j];
            }

            // If a reduced numbering has not been applied to this point,
            // create a new number for it
            if (reduced_points[pt] == 0) {
              reduced_points[pt] = npts;
              npts++;
            }

            // Set the reduced connectivity
            reduced_conn[conn_size * nelems + j] = reduced_points[pt];
          }

          nelems++;
        }
        basic_conn_offset += conn_size;
      }

      // Since we started at npts = 1, we have one more point
      // than the actual number of points.
      npts--;

      // Set the element type to use
      ZoneType zone_type;
      if (zone_btype == TACS_LINE_ELEMENT) {
        zone_type = FELINESEG;
      } else if (zone_btype == TACS_QUAD_ELEMENT) {
        zone_type = FEQUADRILATERAL;
      } else if (zone_btype == TACS_TETRA_ELEMENT) {
        zone_type = FETETRAHEDRON;
      } else if (zone_btype == TACS_HEXA_ELEMENT) {
        zone_type = FEBRICK;
      } else {
        fprintf(stderr,
                "Component %d has unsupported element types for f5totec\n", k);
        return (1);
      }

      float *reduced_float_data = NULL;
      reduced_float_data = new float[npts];

      float *element_float_data = NULL;
      element_float_data = new float[nelems];

      if (nelems > 0 && npts > 0) {
        // Create the zone with the solution time
        create_fe_tec_zone(comp_name, zone_type, npts, nelems, value_location,
                           use_strands, solution_time);

        // Retrieve the continuous data
        for (int j = 0; j < num_variables; j++) {
          for (int i = 0; i < num_points; i++) {
            if (reduced_points[i] > 0) {
              reduced_float_data[reduced_points[i] - 1] =
                  cdata[i * num_variables + j];
            }
          }
          write_tec_float_data(npts, reduced_float_data);
        }

        // Retrieve the element data as cell-centered values
        int count = 0;
        for (int j = 0; j < num_evariables; j++) {
          for (int i = 0; i < num_basic_elements; i++) {
            if (basic_element_comp_num[i] == k) {
              int elem_idx = basic_element_global_ptr[i];
              element_float_data[count] =
                  avg_edata[num_evariables * elem_idx + j];
              count++;
            }
          }
          write_tec_float_data(nelems, element_float_data);
          count = 0;
        }

        // Now, write the connectivity
        write_con_data(reduced_conn);
      }
      // Clean up memory
      delete[] reduced_float_data;
      delete[] element_float_data;
    }

    if (tec_init) {
      close_tec_file();
    }

    loader->decref();

    // Clean up memory
    delete[] reduced_points;
    delete[] reduced_conn;
    delete[] avg_edata;
    delete[] value_location;

    delete[] basic_ltypes;
    delete[] basic_conn;
    delete[] basic_element_comp_num;
    delete[] basic_element_global_ptr;

    delete[] infile;
    delete[] outfile;
  }

  MPI_Finalize();

  return (0);
}
