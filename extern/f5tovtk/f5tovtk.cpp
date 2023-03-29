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
#include "TACSElementTypes.h"
#include "TACSFH5Loader.h"

const int VTK_VERTEX = 1;
const int VTK_LINE = 3;
const int VTK_TRIANGLE = 5;
const int VTK_QUAD = 9;
const int VTK_TETRA = 10;
const int VTK_HEXAHEDRON = 12;
const int VTK_QUADRATIC_TRIANGLE = 22;
const int VTK_QUADRATIC_TETRA = 24;

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Convert hdf5 file argv[1] to
  if (argc == 1) {
    fprintf(stderr, "Error, no input files\n");
    return (1);
  }

  for (int iter = 1; iter < argc; iter++) {
    char *infile = new char[strlen(argv[iter]) + 1];
    strcpy(infile, argv[iter]);

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
    strcpy(&outfile[i], ".vtk");

    printf("Trying to convert FH5 file %s to vtk file %s\n", infile, outfile);

    // Create the loader object
    TACSFH5Loader *loader = new TACSFH5Loader();
    loader->incref();

    int fail = loader->loadData(infile);
    if (fail) {
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
    if (!fp) {
      fprintf(stderr, "Failed to open the output file %s\n", outfile);
      return (1);
    }

    // Write out the vtk file
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\nASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // Write out the points
    fprintf(fp, "POINTS %d double\n", cdim1);

    const float *d = cdata;
    for (int k = 0; k < cdim1; k++) {
      fprintf(fp, "%e %e %e\n", d[0], d[1], d[2]);
      d += cdim2;
    }

    int num_basic_elements = 0;
    int basic_conn_size = 0;
    for (int k = 0; k < num_elements; k++) {
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      if (ltype == TACS_TRI_QUADRATIC_ELEMENT) {
        ntypes = 1;
        nconn = 6;
      } else if (ltype == TACS_TETRA_QUADRATIC_ELEMENT) {
        ntypes = 1;
        nconn = 10;
      } else {
        TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
      }
      num_basic_elements += ntypes;
      basic_conn_size += nconn;
    }

    int *basic_ltypes = new int[num_basic_elements];
    int *basic_conn = new int[basic_conn_size];

    int *btypes = basic_ltypes;
    int *bconn = basic_conn;
    for (int k = 0; k < num_elements; k++) {
      int ntypes = 0, nconn = 0;
      ElementLayout ltype = (ElementLayout)ltypes[k];
      if (ltype == TACS_TRI_QUADRATIC_ELEMENT) {
        btypes[0] = ltype;
        ntypes = 1;
        nconn = 6;
        memcpy(bconn, &conn[ptr[k]], 6 * sizeof(int));
      } else if (ltype == TACS_TETRA_QUADRATIC_ELEMENT) {
        btypes[0] = ltype;
        ntypes = 1;
        nconn = 10;
        memcpy(bconn, &conn[ptr[k]], 10 * sizeof(int));
      } else {
        TacsConvertVisLayoutToBasicCount(ltype, &ntypes, &nconn);
        TacsConvertVisLayoutToBasic(ltype, &conn[ptr[k]], btypes, bconn);
      }
      btypes += ntypes;
      bconn += nconn;
    }

    // Write out the cell values
    fprintf(fp, "\nCELLS %d %d\n", num_basic_elements,
            num_basic_elements + basic_conn_size);

    int basic_conn_offset = 0;
    for (int k = 0; k < num_basic_elements; k++) {
      ElementLayout ltype = (ElementLayout)basic_ltypes[k];
      int conn_size = TacsGetNumVisNodes(ltype);
      fprintf(fp, "%d ", conn_size);
      if (basic_ltypes[k] == TACS_QUAD_ELEMENT) {
        const int convert[] = {0, 1, 3, 2};
        for (int j = 0; j < conn_size; j++) {
          fprintf(fp, "%d ", basic_conn[basic_conn_offset + convert[j]]);
        }
        basic_conn_offset += 4;
      } else if (basic_ltypes[k] == TACS_HEXA_ELEMENT) {
        const int convert[] = {0, 1, 3, 2, 4, 5, 7, 6};
        for (int j = 0; j < conn_size; j++) {
          fprintf(fp, "%d ", basic_conn[basic_conn_offset + convert[j]]);
        }
        basic_conn_offset += 8;
      } else {
        for (int j = 0; j < conn_size; j++, basic_conn_offset++) {
          fprintf(fp, "%d ", basic_conn[basic_conn_offset]);
        }
      }
      fprintf(fp, "\n");
    }

    // All tetrahedrals...
    fprintf(fp, "\nCELL_TYPES %d\n", num_basic_elements);
    for (int k = 0; k < num_basic_elements; k++) {
      if (basic_ltypes[k] == TACS_POINT_ELEMENT) {
        fprintf(fp, "%d\n", VTK_VERTEX);
      } else if (basic_ltypes[k] == TACS_LINE_ELEMENT) {
        fprintf(fp, "%d\n", VTK_LINE);
      } else if (basic_ltypes[k] == TACS_TRI_ELEMENT) {
        fprintf(fp, "%d\n", VTK_TRIANGLE);
      } else if (basic_ltypes[k] == TACS_TRI_QUADRATIC_ELEMENT) {
        fprintf(fp, "%d\n", VTK_QUADRATIC_TRIANGLE);
      } else if (basic_ltypes[k] == TACS_QUAD_ELEMENT) {
        fprintf(fp, "%d\n", VTK_QUAD);
      } else if (basic_ltypes[k] == TACS_TETRA_ELEMENT) {
        fprintf(fp, "%d\n", VTK_TETRA);
      } else if (basic_ltypes[k] == TACS_TETRA_QUADRATIC_ELEMENT) {
        fprintf(fp, "%d\n", VTK_QUADRATIC_TETRA);
      } else if (basic_ltypes[k] == TACS_HEXA_ELEMENT) {
        fprintf(fp, "%d\n", VTK_HEXAHEDRON);
      }
    }
    delete[] basic_conn;
    delete[] basic_ltypes;

    // Print out the rest as fields one-by-one
    fprintf(fp, "POINT_DATA %d\n", cdim1);

    for (int j = 0; j < cdim2; j++) {
      char name[256];
      int index = 0;
      while (strlen(cvars) > 0 && cvars[0] != ',') {
        name[index] = cvars[0];
        index++;
        cvars++;
      }
      name[index] = '\0';
      cvars++;

      // Write out the zone names
      if (j >= 3) {
        fprintf(fp, "SCALARS %s double 1\n", name);
        fprintf(fp, "LOOKUP_TABLE default\n");

        for (int k = 0; k < cdim1; k++) {
          double d = cdata[cdim2 * k + j];
          // If the value is smaller than 10^-15, set it to 0
          // so Paraview won't throw an error
          // if (abs(d) < 1e-15){
          //   d = 0.0;
          // }
          fprintf(fp, "%.3e\n", d);
        }
      }
    }

    // Count up the number of times each node is referred to
    // in the discontinuous element-wise data
    float *counts = new float[cdim1];
    memset(counts, 0, cdim1 * sizeof(float));
    for (int j = 0; j < ptr[num_elements]; j++) {
      counts[conn[j]] += 1.0;
    }
    for (int i = 0; i < cdim1; i++) {
      if (counts[i] != 0.0) {
        counts[i] = 1.0 / counts[i];
      }
    }

    // For each component, average the nodal data
    float *data = new float[cdim1];
    for (int j = 0; j < edim2; j++) {
      char name[256];
      int index = 0;
      while (strlen(evars) > 0 && evars[0] != ',') {
        name[index] = evars[0];
        index++;
        evars++;
      }
      name[index] = '\0';
      evars++;

      // Nodally average the data
      memset(data, 0, cdim1 * sizeof(float));
      for (int k = 0; k < ptr[num_elements]; k++) {
        data[conn[k]] += counts[conn[k]] * edata[edim2 * k + j];
      }

      // Write out the zone names
      fprintf(fp, "SCALARS %s double 1\n", name);
      fprintf(fp, "LOOKUP_TABLE default\n");

      for (int k = 0; k < cdim1; k++) {
        fprintf(fp, "%.3e\n", data[k]);
      }
    }

    delete[] counts;
    delete[] data;

    fclose(fp);

    loader->decref();

    delete[] infile;
    delete[] outfile;
  }

  MPI_Finalize();

  return (0);
}
