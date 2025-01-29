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

/**
 * @brief Convert from TACS element node ordering to tecplot node ordering
 *
 * @param eltype Element type
 * @param ind TACS node number
 * @return int Tecplot node number
 */
int convertLocalNodeInd(const ElementLayout eltype, const int ind) {
  if (eltype == TACS_QUAD_ELEMENT) {
    const int convert[] = {0, 1, 3, 2};
    return convert[ind];
  } else if (eltype == TACS_HEXA_ELEMENT) {
    const int convert[] = {0, 1, 3, 2, 4, 5, 7, 6};
    return convert[ind];
  } else {
    return ind;
  }
}

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
  num_points   == The number of points
  num_elements == The number of elements
*/
int create_fe_tec_zone(char *zone_name, ZoneType _zone_type, int _num_points,
                       int _num_elements, int use_strands = 0,
                       double solution_time = 0.0) {
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
  INTEGER4 *value_location = NULL;    // All values are nodal values
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
    int *element_comp_num, *eltypes, *ptr, *conn;
    loader->getConnectivity(&num_elements, &element_comp_num, &eltypes, &ptr,
                            &conn);

    const char *cname, *var_names;
    int num_points, num_cvars_per_node;
    float *cdata;
    loader->getContinuousData(&cname, &var_names, &num_points,
                              &num_cvars_per_node, &cdata);

    const char *ename, *evar_names;
    int num_enodes, num_evars_per_node;
    float *edata;
    loader->getElementData(&ename, &evar_names, &num_enodes,
                           &num_evars_per_node, &edata);

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

    // Count up the number of times each node is referred to
    // in the discontinuous element-wise data
    float *counts = new float[num_points];
    memset(counts, 0, num_points * sizeof(float));
    for (int j = 0; j < ptr[num_elements]; j++) {
      counts[conn[j]] += 1.0;
    }
    for (int i = 0; i < num_points; i++) {
      if (counts[i] != 0.0) {
        counts[i] = 1.0 / counts[i];
      }
    }

    // Nodally average the data
    float *avg_edata = new float[num_points * num_evars_per_node];
    memset(avg_edata, 0, num_points * num_evars_per_node * sizeof(float));
    for (int j = 0; j < num_evars_per_node; j++) {
      for (int k = 0; k < ptr[num_elements]; k++) {
        avg_edata[num_evars_per_node * conn[k] + j] +=
            counts[conn[k]] * edata[num_evars_per_node * k + j];
      }
    }

    delete[] counts;

    if (!(element_comp_num && conn && cdata)) {
      fprintf(stderr,
              "Error, data, connectivity or \
              component numbers not defined in file\n");
    }

    // Setup visualization elements
    int num_basic_elements = 0;
    int basic_conn_size = 0;
    for (int k = 0; k < num_elements; k++) {
      int numNewElems = 0, numNewNodes = 0;
      const ElementLayout eltype = (ElementLayout)eltypes[k];
      TacsConvertVisLayoutToBasicCount(eltype, &numNewElems, &numNewNodes);
      // For triangular elements we'll add a 4th dummy node
      // to the connectivity that's just the third node repeated
      // This way triangles can be treated as degenerate quads from Tecplot's
      // perspective
      if (eltype == TACS_TRI_ELEMENT || eltype == TACS_TRI_QUADRATIC_ELEMENT ||
          eltype == TACS_TRI_CUBIC_ELEMENT) {
        numNewNodes = numNewNodes + numNewElems;
      }
      // Plot higher order tetrahedral elements as a single linear element
      else if (eltype == TACS_TETRA_QUADRATIC_ELEMENT ||
               eltype == TACS_TETRA_CUBIC_ELEMENT) {
        numNewNodes = 4;
        numNewElems = 1;
      }
      num_basic_elements += numNewElems;
      basic_conn_size += numNewNodes;
    }

    int *const basic_element_data_map =
        new int[num_basic_elements];  // For each basic element, where does the
                                      // data for the corresponding original
                                      // element start in the element data array
    int *elem_data_map_ptr = basic_element_data_map;

    int *basic_eltypes = new int[num_basic_elements];
    int *basic_element_comp_num = new int[num_basic_elements];
    int *basic_conn = new int[basic_conn_size];
    int *const basic_node_map =
        new int[basic_conn_size];  // Which node in the original element is this
                                   // node in the basic element
    int *basic_node_map_ptr = basic_node_map;

    int *btypes = basic_eltypes;
    int *belem_comp_num = basic_element_comp_num;
    int *bconn = basic_conn;

    int origNodeCount = 0;

    for (int k = 0; k < num_elements; k++) {
      int numNewElems = 0, numNewNodes = 0;
      ElementLayout eltype = (ElementLayout)eltypes[k];
      const int nodes_per_elem = TacsGetNumVisNodes(eltype);
      // Add our dummy nodes for triangular elements
      if (eltype == TACS_TRI_ELEMENT || eltype == TACS_TRI_QUADRATIC_ELEMENT ||
          eltype == TACS_TRI_CUBIC_ELEMENT) {
        TacsConvertVisLayoutToBasicCount(eltype, &numNewElems, &numNewNodes);
        int *tri_conn = new int[numNewNodes];
        TacsConvertVisLayoutToBasic(eltype, &conn[ptr[k]], btypes, tri_conn);

        const int convert[] = {0, 1, 2, 2};
        for (int jj = 0; jj < numNewElems; jj++) {
          for (int ii = 0; ii < 4; ii++) {
            bconn[4 * jj + ii] = tri_conn[3 * jj + convert[ii]];
          }
          btypes[jj] = TACS_QUAD_ELEMENT;
        }
        numNewNodes = numNewNodes + numNewElems;
        delete[] tri_conn;
      }
      // Plot only the first four nodes (conrners) of higher order tets
      else if (eltype == TACS_TETRA_QUADRATIC_ELEMENT ||
               eltype == TACS_TETRA_CUBIC_ELEMENT) {
        memcpy(bconn, &conn[ptr[k]], 4 * sizeof(int));
        numNewNodes = 4;
        numNewElems = 1;
        btypes[0] = TACS_TETRA_ELEMENT;
      } else {
        TacsConvertVisLayoutToBasicCount(eltype, &numNewElems, &numNewNodes);
        TacsConvertVisLayoutToBasic(eltype, &conn[ptr[k]], btypes, bconn);
      }
      // Set the basic element component to match the parent
      for (int ii = 0; ii < numNewElems; ii++) {
        belem_comp_num[ii] = element_comp_num[k];
      }
      // Copy data to the basic element data
      for (int ii = 0; ii < numNewNodes; ii++) {
        // Find which node this is in the original element
        int elemLocalInd = -1;
        for (int jj = 0; jj < nodes_per_elem; jj++) {
          if (conn[ptr[k] + jj] == bconn[ii]) {
            elemLocalInd = jj;
            break;
          }
        }
        if (elemLocalInd == -1) {
          fprintf(stderr, "Error, node not found in original element\n");
          return 1;
        }
        basic_node_map_ptr[ii] = elemLocalInd;
      }

      for (int jj = 0; jj < numNewElems; jj++) {
        elem_data_map_ptr[jj] = origNodeCount * num_evars_per_node;
      }
      origNodeCount += nodes_per_elem;
      elem_data_map_ptr += numNewElems;

      btypes += numNewElems;
      belem_comp_num += numNewElems;
      bconn += numNewNodes;
      basic_node_map_ptr += numNewNodes;
    }

    int num_comp = loader->getNumComponents();

    int *reduced_points = new int[num_points];
    int *reduced_conn = new int[basic_conn_size];
    float *reduced_float_data = NULL;
    reduced_float_data = new float[num_points];
    int *const compNodeGlobalInd =
        new int[basic_conn_size];  // Map from the duplicated node index within
                                   // a component to the index of that node in
                                   // the original full mesh
    int *const compBasicElemInd =
        new int[num_basic_elements];  // Map from the index of a basic element
                                      // within a component to the index of that
                                      // basic element in the full mesh

    for (int k = 0; k < num_comp; k++) {
      // Count up the number of elements that use the connectivity
      char *comp_name = loader->getComponentName(k);
      // printf("Converting zone %d: %s at time %g\n",
      //  k, comp_name, solution_time);

      memset(reduced_points, 0, num_points * sizeof(int));
      memset(compNodeGlobalInd, 0, basic_conn_size * sizeof(int));
      memset(reduced_conn, 0, basic_conn_size * sizeof(int));
      memset(compBasicElemInd, 0, num_basic_elements * sizeof(int));

      int npts = 1, num_elements = 0;
      int zone_btype = -1;
      int basic_conn_offset = 0;
      // Count up the number of points/elements in this sub-domain
      for (int i = 0; i < num_basic_elements; i++) {
        ElementLayout eltype = (ElementLayout)basic_eltypes[i];
        int num_node_per_elem = TacsGetNumVisNodes(eltype);

        if (basic_element_comp_num[i] == k) {
          // Make sure all elements in this zone are the same type
          if (zone_btype == -1) {
            zone_btype = basic_eltypes[i];
          } else if (zone_btype != basic_eltypes[i]) {
            fprintf(stderr, "Component %d has conflicting element types\n", k);
            return (1);
          }
          compBasicElemInd[num_elements] = i;

          int pt;
          for (int j = 0; j < num_node_per_elem; j++) {
            // Add this element to the reduced connectivity
            pt = basic_conn[basic_conn_offset + convertLocalNodeInd(eltype, j)];
            compNodeGlobalInd[num_elements * num_node_per_elem + j] = pt;

            // If a reduced numbering has not been applied to this point,
            // create a new number for it
            if (reduced_points[pt] == 0) {
              reduced_points[pt] = npts;
              npts++;
            }

            // Set the reduced connectivity
            reduced_conn[num_node_per_elem * num_elements + j] =
                reduced_points[pt];
          }

          num_elements++;
        }
        basic_conn_offset += num_node_per_elem;
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

      if (num_elements > 0 && npts > 0) {
        const bool duplicateNodes = true;
        if (duplicateNodes) {
          // In this case, we duplicate nodes so that every element has its own
          // nodes (e.g first element has nodes 0->n, second elements has nodes
          // n+1->2n etc) and we can plot discontinuous data

          // Create the zone with the solution time
          const int nodes_per_elem =
              TacsGetNumVisNodes((ElementLayout)zone_btype);
          const int num_dup_nodes = num_elements * nodes_per_elem;
          create_fe_tec_zone(comp_name, zone_type, num_dup_nodes, num_elements,
                             use_strands, solution_time);

          // Allocate a new array for the nodal data
          float *const nodalData = new float[num_dup_nodes];

          // Continuous data
          for (int jj = 0; jj < num_cvars_per_node; jj++) {
            for (int elemInd = 0; elemInd < num_elements; elemInd++) {
              for (int nodeInd = 0; nodeInd < nodes_per_elem; nodeInd++) {
                const int dupNodeInd = elemInd * nodes_per_elem + nodeInd;
                const int globalNodeInd = compNodeGlobalInd[dupNodeInd];
                nodalData[dupNodeInd] =
                    cdata[globalNodeInd * num_cvars_per_node + jj];
              }
            }
            write_tec_float_data(num_dup_nodes, nodalData);
          }
          // Element data
          // The TACS data is stored such that all variable values for a given
          // node are contiguous, tecplot wants the data to be stored such that
          // all values for a given variable are contiguous
          for (int jj = 0; jj < num_evars_per_node; jj++) {
            for (int elemInd = 0; elemInd < num_elements; elemInd++) {
              const int globalElemInd =
                  compBasicElemInd[elemInd];  // Global basic element index
              const int elemDataStartInd =
                  basic_element_data_map[globalElemInd];  // Start of data for
                                                          // this element in the
                                                          // element data array
              const float *const origElementData =
                  &edata[elemDataStartInd];  // Start of data for this element

              for (int nodeInd = 0; nodeInd < nodes_per_elem; nodeInd++) {
                // At this point I know:
                // - Which basic element I'm handling
                // - Which original element that basic element came from
                // - Where in the element data array does the data for this
                // original element start
                // - Which node (in TACS ordering) within the basic element I'm
                // handling
                // - Which node of the basic element this in the tecplot
                // ordering
                // - Which node in the original element is this

                const int dupNodeInd = elemInd * nodes_per_elem + nodeInd;

                const int basicElemLocalNodeInd =
                    convertLocalNodeInd((ElementLayout)zone_btype, nodeInd);

                const int originalElemLocalNodeInd =
                    basic_node_map[globalElemInd * nodes_per_elem +
                                   basicElemLocalNodeInd];

                nodalData[dupNodeInd] =
                    origElementData[originalElemLocalNodeInd *
                                        num_evars_per_node +
                                    jj];
              }
            }
            write_tec_float_data(num_dup_nodes, nodalData);
          }

          // Now create the connectivity, because of the way we duplicate the
          // nodal values, the connectivity is just 1,2,3,4...
          int *const dupConn = new int[num_dup_nodes];
          for (int i = 0; i < num_dup_nodes; i++) {
            dupConn[i] = i + 1;
          }
          write_con_data(dupConn);
          delete[] nodalData;
          delete[] dupConn;
        } else {
          // Create the zone with the solution time
          create_fe_tec_zone(comp_name, zone_type, npts, num_elements,
                             use_strands, solution_time);

          // Retrieve the continuous data
          for (int j = 0; j < num_cvars_per_node; j++) {
            for (int i = 0; i < num_points; i++) {
              if (reduced_points[i] > 0) {
                reduced_float_data[reduced_points[i] - 1] =
                    cdata[i * num_cvars_per_node + j];
              }
            }
            write_tec_float_data(npts, reduced_float_data);
          }

          // Retrieve the element data
          for (int j = 0; j < num_evars_per_node; j++) {
            for (int i = 0; i < num_points; i++) {
              if (reduced_points[i] > 0) {
                reduced_float_data[reduced_points[i] - 1] =
                    avg_edata[num_evars_per_node * i + j];
              }
            }
            write_tec_float_data(npts, reduced_float_data);
          }

          // Now, write the connectivity
          write_con_data(reduced_conn);
        }
      }
    }

    if (tec_init) {
      close_tec_file();
    }

    loader->decref();

    // Clean up memory
    delete[] reduced_points;
    delete[] reduced_conn;
    delete[] reduced_float_data;
    delete[] avg_edata;

    delete[] basic_eltypes;
    delete[] basic_conn;
    delete[] basic_element_comp_num;
    delete[] basic_node_map;
    delete[] basic_element_data_map;
    delete[] compNodeGlobalInd;
    delete[] compBasicElemInd;

    delete[] infile;
    delete[] outfile;
  }

  MPI_Finalize();

  return (0);
}
