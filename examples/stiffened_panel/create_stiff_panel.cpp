#include <math.h>
#include <stdio.h>
#include <string.h>

/*
  Create a .bdf file for a stiffened plate

  input:
  file_name: the .bdf file name containing the stiffened panel model
  Nx: the number of elements in the x direction
  nrepeat: the number of times the geometry repeats itself

  geometric input:
  Lx: the length of the panel in the x direction
  b: the width of the repeating segment
  wb: the width of the base of the stiffener
  hs: the height of the stiffener


  The stiffened plate consists of a series of stiffeners that have
  the following geometry:

  |-->  b/2 <---|-->  b/2 <---|
  .         |-- wb--|

  |----------=======----------|  -
  .             |                |
  .             |                hs
  .             |                |
  .                              -


  This repeats 'nrepeat' times.

  The aspect ratios of the elements are chosen to be approximately 1
*/
void write_stiffened_panel_axial(const char* file_name, int Nx, int nrepeat,
                                 double Lx, double b, double wb, double hs,
                                 double angle, int write_shear_file) {
  // Check to make sure everything makes sense
  if (nrepeat < 1) {
    nrepeat = 1;
  }
  if (wb >= b) {
    printf(
        "The width of the base must be less "
        "than the width of the segment\n");
  }

  const int elem_order = 3;
  FILE* fp = fopen(file_name, "w");

  fprintf(fp, "SOL 103\n");
  fprintf(fp, "CEND\n");
  fprintf(fp, "BEGIN BULK\n");

  // Compute the number of elements on the surface

  // Nyb = the number of elements over the segment (b - wb)/2
  int Nyb = int((Nx * (0.5 * (b - wb)) / Lx));
  if (Nyb < 1) {
    Nyb = 1;
  }

  // Nywb = the number of elements over the segment wb/2
  int Nywb = int((Nx * (0.5 * wb)) / Lx);
  if (Nywb < 1) {
    Nywb = 1;
  }

  // Nhs = the number of elements along the stiffener
  int Nhs = int((Nx * hs) / Lx);
  if (Nhs < 1) {
    Nhs = 1;
  }

  // The number of elements along the length
  // int nelems = (Nhs + 2*(Nyb + Nywb))*Nx*nrepeat;

  int Ny = 2 * (Nyb + Nywb) * nrepeat;

  // Determine the ordering of the elements on the surface
  int nx = (elem_order - 1) * Nx + 1;
  int ny = (elem_order - 1) * Ny + 1;
  int nhs = (elem_order - 1) * Nhs + 1;

  // Create lists of node numbers
  int* surf_node_nums = new int[nx * ny];
  int* stiff_node_nums = new int[nrepeat * nx * nhs];

  double beta = tan(angle);
  double* Xpanel = new double[nx];
  double* Ysurf = new double[ny];
  double* Zstiff = new double[nhs];

  int nnodes = 1;
  // Order the nodes on the surface
  for (int j = 0; j < ny; j++) {
    for (int i = 0; i < nx; i++) {
      surf_node_nums[ny * i + j] = nnodes;
      nnodes++;
    }
  }

  // Order the nodes on the stiffeners
  for (int jj = 0; jj < nrepeat; jj++) {
    int surf_node_offset = (elem_order - 1) * (Nyb + Nywb) * (2 * jj + 1);

    for (int j = 0; j < nhs - 1; j++) {
      for (int i = 0; i < nx; i++) {
        stiff_node_nums[jj * nx * nhs + nhs * i + j] = nnodes;
        nnodes++;
      }
    }

    for (int i = 0; i < nx; i++) {
      stiff_node_nums[jj * nx * nhs + nhs * (i + 1) - 1] =
          surf_node_nums[ny * i + surf_node_offset];
    }
  }

  // Set the nodal locations
  for (int i = 0; i < nx; i++) {
    Xpanel[i] = (Lx * i) / (nx - 1);
  }

  for (int k = 0; k < nhs; k++) {
    Zstiff[k] = -hs + (hs * k) / (nhs - 1);
  }

  for (int j = 0, jj = 0; jj < nrepeat; jj++) {
    double offset = b * jj;
    for (int k = 0; k < (elem_order - 1) * Nyb; k++, j++) {
      Ysurf[j] = offset + (0.5 * (b - wb) * k) / ((elem_order - 1) * Nyb);
    }

    offset += 0.5 * (b - wb);
    for (int k = 0; k < 2 * (elem_order - 1) * Nywb; k++, j++) {
      Ysurf[j] = offset + (wb * k) / (2 * (elem_order - 1) * Nywb);
    }

    offset += wb;
    for (int k = 0; k < (elem_order - 1) * Nyb; k++, j++) {
      Ysurf[j] = offset + (0.5 * (b - wb) * k) / ((elem_order - 1) * Nyb);
    }
  }

  Ysurf[ny - 1] = b * nrepeat;

  /*
    Write out the nodal locations
  */

  int coord_id = 0;
  int coord_disp = 0;
  int seid = 0;
  const char* spc = " ";

  for (int jj = 0, j = 0; jj < nrepeat; jj++) {
    // Write out the surface locations
    for (int k = 0; k < 2 * (elem_order - 1) * (Nyb + Nywb); k++, j++) {
      for (int i = 0; i < nx; i++) {
        int node_num = surf_node_nums[ny * i + j];

        fprintf(fp, "%-8s%16d%16d%16.9e%16.9e*       \n", "GRID*", node_num,
                coord_id, Xpanel[i] + beta * Ysurf[j], Ysurf[j]);
        fprintf(fp, "*       %16.9e%16d%16s%16d        \n", 0.0, coord_disp,
                spc, seid);
      }
    }

    // Write out the stiffener locations
    double y = Ysurf[(elem_order - 1) * (Nyb + Nywb) * (2 * jj + 1)];
    for (int i = 0; i < nx; i++) {
      for (int k = 0; k < (elem_order - 1) * Nhs; k++) {
        int node_num = stiff_node_nums[jj * nx * nhs + nhs * i + k];

        fprintf(fp, "%-8s%16d%16d%16.9e%16.9e*       \n", "GRID*", node_num,
                coord_id, Xpanel[i] + beta * y, y);
        fprintf(fp, "*       %16.9e%16d%16s%16d        \n", Zstiff[k],
                coord_disp, spc, seid);
      }
    }
  }

  // Write the last set of nodes
  for (int i = 0; i < nx; i++) {
    int node_num = surf_node_nums[ny * i + ny - 1];

    fprintf(fp, "%-8s%16d%16d%16.9e%16.9e*       \n", "GRID*", node_num,
            coord_id, Xpanel[i] + beta * Ysurf[ny - 1], Ysurf[ny - 1]);
    fprintf(fp, "*       %16.9e%16d%16s%16d        \n", 0.0, coord_disp, spc,
            seid);
  }

  delete[] Xpanel;
  delete[] Ysurf;
  delete[] Zstiff;

  /*
    Write the element connectivity
  */

  int elem = 1;

  for (int jj = 0, j = 0; jj < nrepeat; jj++) {
    int nodes[16];

    // Write the elements on the skin
    for (int k = 0; k < 2 * (Nyb + Nywb); k++, j++) {
      for (int i = 0; i < Nx; i++) {
        for (int n = 0; n < elem_order; n++) {
          for (int m = 0; m < elem_order; m++) {
            nodes[elem_order * n + m] =
                surf_node_nums[ny * ((elem_order - 1) * i + m) +
                               (elem_order - 1) * j + n];
          }
        }

        int part_id = 1;
        if (k >= Nyb && k < Nyb + 2 * Nywb) {
          part_id = 2;
        }

        fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n", "CQUAD9", elem, part_id,
                nodes[0], nodes[2], nodes[8], nodes[6], nodes[1], nodes[5]);
        fprintf(fp, "        %8d%8d%8d\n", nodes[7], nodes[3], nodes[4]);
        elem += 1;
      }
    }

    // Write the nodes on the stiffener
    for (int k = 0; k < Nhs; k++) {
      for (int i = 0; i < Nx; i++) {
        for (int n = 0; n < elem_order; n++) {
          for (int m = 0; m < elem_order; m++) {
            nodes[elem_order * n + m] =
                stiff_node_nums[nhs * nx * jj +
                                nhs * ((elem_order - 1) * i + m) +
                                (elem_order - 1) * k + n];
          }
        }

        int part_id = 3;
        fprintf(fp, "%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n", "CQUAD9", elem, part_id,
                nodes[0], nodes[2], nodes[8], nodes[6], nodes[1], nodes[5]);
        fprintf(fp, "        %8d%8d%8d\n", nodes[7], nodes[3], nodes[4]);
        elem += 1;
      }
    }
  }

  if (write_shear_file) {
    double scale = -0.001;
    // Clamp w, rotx at one side
    for (int j = 0; j < ny; j++) {
      int node = surf_node_nums[j];
      double v = scale * (Xpanel[0] + beta * Ysurf[j]);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "345",
              0.0);  // w, rotx, roty
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "2", v);  // v
    }

    for (int jj = 0; jj < nrepeat; jj++) {
      for (int j = 0; j < nhs; j++) {
        int node = stiff_node_nums[jj * nx * nhs + j];
        double u =
            scale * Ysurf[(elem_order - 1) * (Nyb + Nywb) * (2 * jj + 1)];
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", u);  // u
      }
    }

    // Clamp/set bcs at the other end
    for (int j = 0; j < ny; j++) {
      int node = surf_node_nums[ny * (nx - 1) + j];
      double v = scale * (Xpanel[nx - 1] + beta * Ysurf[j]);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "345", 0.0);  // w
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "2", v);      // v
    }

    for (int jj = 0; jj < nrepeat; jj++) {
      for (int j = 0; j < nhs; j++) {
        int node = stiff_node_nums[jj * nx * nhs + nhs * (nx - 1) + j];
        double u =
            scale * Ysurf[(elem_order - 1) * (Nyb + Nywb) * (2 * jj + 1)];
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", u);  // u
      }
    }

    // Set simply supported boundary conditions along the edge
    for (int i = 0; i < nx; i++) {
      int node = surf_node_nums[ny * i];
      double u = scale * Ysurf[0];
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "345", 0.0);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", u);

      node = surf_node_nums[ny * (i + 1) - 1];
      u = scale * Ysurf[ny - 1];
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "345", 0.0);
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", u);
    }
  } else {
    // Clamp w, rotx at one side
    for (int j = 0; j < ny; j++) {
      int node = surf_node_nums[j];
      if (j == 0) {
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "12345",
                0.0);  // u, w
      } else {
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1345",
                0.0);  // u, w
      }
    }

    for (int jj = 0; jj < nrepeat; jj++) {
      for (int j = 0; j < nhs; j++) {
        int node = stiff_node_nums[jj * nx * nhs + j];
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", 0.0);  // u
      }
    }

    // Clamp/set bcs at the other end
    for (int j = 0; j < ny; j++) {
      int node = surf_node_nums[ny * (nx - 1) + j];
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "345", 0.0);  // w
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", -1.0);   // u
    }

    for (int jj = 0; jj < nrepeat; jj++) {
      for (int j = 0; j < nhs; j++) {
        int node = stiff_node_nums[jj * nx * nhs + nhs * (nx - 1) + j];
        fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, "1", -1.0);  // u
      }
    }

    // Set simply supported boundary conditions along the edge
    const char* left_spc = "345";
    const char* right_spc = "345";

    for (int i = 0; i < nx; i++) {
      int node = surf_node_nums[ny * i];
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, left_spc, 0.0);
      node = surf_node_nums[ny * (i + 1) - 1];
      fprintf(fp, "%-8s%8d%8d%8s%8.6f\n", "SPC", 1, node, right_spc, 0.0);
    }
  }

  delete[] surf_node_nums;
  delete[] stiff_node_nums;

  fprintf(fp, "END BULK\n");
  fclose(fp);
}

int main(int argc, char* argv[]) {
  int Nx = 45;
  int nrepeat = 4;
  double Lx = 450.0;
  double b = 110.0;
  double wb = 35.0;
  double hs = 20.0;
  double angle = -15.0 / 180.0 * M_PI;

  int write_shear_file = 0;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "shear") == 0) {
      write_shear_file = 1;
    }
    if (sscanf(argv[k], "Nx=%d", &Nx) == 1) {
      if (Nx < 4) {
        Nx = 4;
      }
      if (Nx > 200) {
        Nx = 200;
      }
    }
  }

  if (write_shear_file) {
    write_stiffened_panel_axial("shear_stiffened_panel.bdf", Nx, nrepeat, Lx, b,
                                wb, hs, angle, write_shear_file);
  } else {
    write_stiffened_panel_axial("axial_stiffened_panel.bdf", Nx, nrepeat, Lx, b,
                                wb, hs, angle, write_shear_file);
  }

  return (0);
}
