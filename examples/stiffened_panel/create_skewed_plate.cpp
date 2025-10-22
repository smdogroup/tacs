#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
  Get the coordinates of the deformed plate
*/
void get_coordinates(double* x, double* y, const double a, const double b,
                     const double beta, double u, double v) {
  double uvals[] = {0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0};

  double vvals[] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

  double defect = 0.08;

  uvals[1] = 0.5 - defect;
  uvals[4] = 0.5 + defect;
  uvals[7] = 0.5 - defect;

  vvals[3] = 0.5 - defect;
  vvals[4] = 0.5 + defect;
  vvals[5] = 0.5 - defect;

  int ii = 0;
  if (u < 0.5) {
    u = 4.0 * u - 1.0;
    ii = 0;
  } else {
    u = 4.0 * (u - 0.5) - 1.0;
    ii = 1;
  }

  int jj = 0;
  if (v < 0.5) {
    v = 4.0 * v - 1.0;
    jj = 0;
  } else {
    v = 4.0 * (v - 0.5) - 1.0;
    jj = 1;
  }

  double N[4] = {0.25 * (1.0 + u) * (1.0 + v), 0.25 * (1.0 - u) * (1.0 + v),
                 0.25 * (1.0 + u) * (1.0 - v), 0.25 * (1.0 - u) * (1.0 - v)};

  u = (N[0] * uvals[ii + 3 * jj] + N[1] * uvals[ii + 1 + 3 * jj] +
       N[2] * uvals[ii + 3 * (jj + 1)] + N[3] * uvals[ii + 1 + 3 * (jj + 1)]);

  v = (N[0] * vvals[ii + 3 * jj] + N[1] * vvals[ii + 1 + 3 * jj] +
       N[2] * vvals[ii + 3 * (jj + 1)] + N[3] * vvals[ii + 1 + 3 * (jj + 1)]);

  // Compute the coordinates
  *y = b * v;
  *x = a * u + beta * b * v;
}

/*
  Create a .bdf file for a skewed plate for buckling verification.

  input:
  file_name: the file name
  Nx, Ny: the number of elements along the x and y direction
  a, b: the length and width of the plate
  theta: the skew angle of the plate
*/

void write_skewed_plate_file(const char* file_name, int Nx, int Ny, double a,
                             double b, double theta, double gamma,
                             const char* bc_type) {
  double beta = tan(theta);

  const int elem_order = 4;
  int nx = (elem_order - 1) * Nx + 1;
  int ny = (elem_order - 1) * Ny + 1;

  int coord_id = 0;
  int coord_disp = 0;
  int seid = 0;
  const char* spc = " ";

  FILE* fp = fopen(file_name, "w");
  fprintf(fp, "SOL 103\n");
  fprintf(fp, "CEND\n");
  fprintf(fp, "BEGIN BULK\n");

  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      int node_num = i + j * nx + 1;

      double x = 0.0, y = 0.0;
      double u = (1.0 * i) / (nx - 1);
      double v = (1.0 * j) / (ny - 1);
      get_coordinates(&x, &y, a, b, beta, u, v);

      fprintf(fp, "%-8s%16d%16d%16.9e%16.9e*       \n", "GRID*", node_num,
              coord_id, x, y);
      fprintf(fp, "*       %16.9e%16d%16s%16d        \n", 0.0, coord_disp, spc,
              seid);
    }
  }

  // Write the elements out to the file
  int nodes[16], elem = 1;
  for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
      for (int n = 0; n < elem_order; n++) {
        for (int m = 0; m < elem_order; m++) {
          nodes[elem_order * n + m] =
              (elem_order - 1) * i + m + ((elem_order - 1) * j + n) * nx + 1;
        }
      }

      int part_id = 1;
      fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n", "CQUAD16", elem, part_id,
              nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5]);
      fprintf(fp, "        %8d%8d%8d%8d%8d%8d%8d%8d\n", nodes[6], nodes[7],
              nodes[8], nodes[9], nodes[10], nodes[11], nodes[12], nodes[13]);
      fprintf(fp, "        %8d%8d\n", nodes[14], nodes[15]);
      elem += 1;
    }
  }

  if (strcmp(bc_type, "shear") == 0) {
    // The shear boundary conditions
    // Set the BCs at the x = 0 edge
    for (int j = 0; j < ny; j++) {
      int node_num = j * nx + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "3", 0.0);

      // Calculate the displacements u, v
      double x = 0.0, y = 0.0;
      double u = 0.0;
      double v = (1.0 * j) / (ny - 1);
      get_coordinates(&x, &y, a, b, beta, u, v);

      u = 0.5 * gamma * y;
      v = 0.5 * gamma * x;
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "2", v);
    }

    // Set the BCs at the x = a edge
    for (int j = 0; j < ny; j++) {
      int node_num = (j + 1) * nx;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "3", 0.0);

      // Calculate the displacements u, v
      double x = 0.0, y = 0.0;
      double u = 1.0;
      double v = (1.0 * j) / (ny - 1);
      get_coordinates(&x, &y, a, b, beta, u, v);

      u = 0.5 * gamma * y;
      v = 0.5 * gamma * x;
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "2", v);
    }

    // Set the BCs at the y = 0 edge
    for (int i = 1; i < nx - 1; i++) {
      int node_num = i + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "35", 0.0);

      // Calculate the displacements u, v
      double x = 0.0, y = 0.0;
      double u = (1.0 * i) / (nx - 1);
      double v = 0.0;
      get_coordinates(&x, &y, a, b, beta, u, v);

      u = 0.5 * gamma * y;
      v = 0.5 * gamma * x;
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "2", v);
    }

    // Set the BCs at the y = b edge
    for (int i = 1; i < nx - 1; i++) {
      int node_num = i + nx * (ny - 1) + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "35", 0.0);

      // Calculate the displacements u, v
      double x = 0.0, y = 0.0;
      double u = (1.0 * i) / (nx - 1);
      double v = 1.0;
      get_coordinates(&x, &y, a, b, beta, u, v);

      u = 0.5 * gamma * y;
      v = 0.5 * gamma * x;
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "2", v);
    }
  } else {
    // The end-load boundary conditions
    // Set the BCs at the x = 0 edge
    // u, v, w, rotx, roty, rotz
    fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, 1, "12345", 0.0);

    for (int j = 1; j < ny; j++) {
      int node_num = j * nx + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "345", 0.0);

      double y = (b * j) / (ny - 1);
      double x = beta * y;
      double u = gamma * x;

      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
    }

    // Set the BCs at the x = a edge
    for (int j = 0; j < ny; j++) {
      int node_num = (j + 1) * nx;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "345", 0.0);

      double y = (b * j) / (ny - 1);
      double x = beta * y + a;
      double u = gamma * x;
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "1", u);
    }

    // Set the BCs at the y = 0 edge
    for (int i = 1; i < nx - 1; i++) {
      int node_num = i + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "345", 0.0);
    }

    // Set the BCs at the y = b edge
    for (int i = 1; i < nx - 1; i++) {
      int node_num = i + nx * (ny - 1) + 1;

      // w, rotx, roty, rotz
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node_num, "345", 0.0);
    }
  }

  fprintf(fp, "END BULK\n");
  fclose(fp);
}

int main(int argc, char* argv[]) {
  int Nx = 16;
  int Ny = 16;

  int id = 0;
  double theta = 0.0;
  if (argc >= 2) {
    theta = atof(argv[1]);
    theta *= (M_PI / 180.0);
  }
  if (argc >= 3) {
    id = atoi(argv[2]);
  }

  double a = 5.0;
  double b = a * cos(theta);
  double gamma = 0.001;

  // if (id == 0){ Nx = -1.0, Nxy = 0.0; }
  // else if (id == 1){ Nx = 0.0, Nxy = 1.0; }
  // else if (id == 2){ Nx = 0.0, Nxy = -1.0; }

  if (id == 0) {
    write_skewed_plate_file("skewed_plate.bdf", Nx, Ny, a, b, theta, gamma,
                            "normal");
  } else if (id == 1) {
    write_skewed_plate_file("skewed_plate.bdf", Nx, Ny, a, b, theta, gamma,
                            "shear");
  } else if (id == 2) {
    write_skewed_plate_file("skewed_plate.bdf", Nx, Ny, a, b, theta, -gamma,
                            "shear");
  }
  return (0);
}
