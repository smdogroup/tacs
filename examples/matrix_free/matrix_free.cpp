/*
  This code tests the implementation of higher-order elements with
  the matrix-free matrix matrix-vector product implementation.
*/

#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElement2D.h"
#include "TACSElement3D.h"
#include "TACSHeatConduction.h"
#include "TACSHexaBasis.h"
#include "TACSLagrangeInterpolation.h"
#include "TACSLinearElasticity.h"
#include "TACSMatrixFreeMat.h"
#include "TACSMg.h"
#include "TACSQuadBasis.h"
#include "TACSThermoelasticity.h"
#include "TACSToFH5.h"

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int varsPerNode, int order, int nx, int ny,
                     int nz, TACSAssembler **_assembler,
                     TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (!(varsPerNode == 1 || varsPerNode == 3 || varsPerNode == 4)) {
    varsPerNode = 4;
  }

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);
  creator->incref();

  if (rank == 0) {
    // Set the number of elements
    int numNodes = ((order - 1) * nx + 1) * ((order - 1) * ny + 1) *
                   ((order - 1) * nz + 1);
    int numElements = nx * ny * nz;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[order * order * order * numElements];

    // Set the ids
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    for (int n = 0; n < numElements; n++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = (n % (nx * ny)) % nx;
      int j = (n % (nx * ny)) / nx;
      int k = n / (nx * ny);

      // Set the node connectivity
      for (int kk = 0; kk < order; kk++) {
        for (int jj = 0; jj < order; jj++) {
          for (int ii = 0; ii < order; ii++) {
            conn[order * order * order * n + ii + order * jj +
                 order * order * kk] =
                ((order - 1) * i + ii) +
                ((order - 1) * j + jj) * ((order - 1) * nx + 1) +
                ((order - 1) * k + kk) * ((order - 1) * nx + 1) *
                    ((order - 1) * ny + 1);
          }
        }
      }
      ptr[n + 1] = order * order * order * (n + 1);
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    // We're over-counting one of the nodes on each edge
    int numBcs = ((order - 1) * nx + 1) * ((order - 1) * ny + 1);
    int *bcNodes = new int[numBcs];

    for (int i = 0; i < numBcs; i++) {
      bcNodes[i] = i;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete[] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];
    double hx = 1.0 / nx;
    double hy = 1.0 / ny;
    double hz = 1.0 / nz;

    const double *pts = NULL;
    if (order == 2) {
      pts = TacsGaussLobattoPoints2;
    }
    if (order == 3) {
      pts = TacsGaussLobattoPoints3;
    } else if (order == 4) {
      pts = TacsGaussLobattoPoints4;
    } else if (order == 5) {
      pts = TacsGaussLobattoPoints5;
    } else if (order == 6) {
      pts = TacsGaussLobattoPoints6;
    }

    for (int k = 0; k < (order - 1) * nz + 1; k++) {
      for (int j = 0; j < (order - 1) * ny + 1; j++) {
        for (int i = 0; i < (order - 1) * nx + 1; i++) {
          int node = i + j * ((order - 1) * nx + 1) +
                     k * ((order - 1) * nx + 1) * ((order - 1) * ny + 1);
          Xpts[3 * node] =
              hx * (i / (order - 1) + 0.5 * (1.0 + pts[i % (order - 1)]));
          Xpts[3 * node + 1] =
              hy * (j / (order - 1) + 0.5 * (1.0 + pts[j % (order - 1)]));
          Xpts[3 * node + 2] =
              hz * (k / (order - 1) + 0.5 * (1.0 + pts[k % (order - 1)]));
        }
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete[] Xpts;
  }

  // Create the isotropic material class
  TacsScalar rho = 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E1 = 70e3;
  TacsScalar E2 = 70e-3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props1 =
      new TACSMaterialProperties(rho, specific_heat, E1, nu, ys, cte, kappa);
  TACSMaterialProperties *props2 =
      new TACSMaterialProperties(rho, specific_heat, E2, nu, ys, cte, kappa);

  // Create the stiffness object
  TACSSolidConstitutive *stiff1 = new TACSSolidConstitutive(props1);
  TACSSolidConstitutive *stiff2 = new TACSSolidConstitutive(props2);

  // Create the model class
  TACSElementModel *model1 = NULL, *model2 = NULL;
  if (varsPerNode == 1) {
    model1 = new TACSHeatConduction3D(stiff1);
    model2 = new TACSHeatConduction3D(stiff2);
  } else if (varsPerNode == 3) {
    model1 = new TACSLinearElasticity3D(stiff1, TACS_LINEAR_STRAIN);
    model2 = new TACSLinearElasticity3D(stiff2, TACS_LINEAR_STRAIN);
  } else if (varsPerNode == 4) {
    model1 = new TACSLinearThermoelasticity3D(stiff1, TACS_LINEAR_STRAIN);
    model2 = new TACSLinearThermoelasticity3D(stiff2, TACS_LINEAR_STRAIN);
  }

  // Create the element class
  TACSElementBasis *basis = NULL;
  if (order == 2) {
    basis = new TACSLinearHexaBasis();
  } else if (order == 3) {
    basis = new TACSQuadraticHexaBasis();
  } else if (order == 4) {
    basis = new TACSCubicHexaBasis();
  } else if (order == 5) {
    basis = new TACSQuarticHexaBasis();
  } else if (order == 6) {
    basis = new TACSQuinticHexaBasis();
  }
  TACSElement3D *linear_element1 = new TACSElement3D(model1, basis);
  TACSElement3D *linear_element2 = new TACSElement3D(model2, basis);

  // Set the one element
  TACSElement *elems[2];
  elems[0] = linear_element1;
  elems[1] = linear_element2;
  creator->setElements(2, elems);

  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();
  assembler->incref();

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

/*
  Do a direct comparison in time and result against the in-memory
  matrix-vector products
*/
void testMatrixVectorProducts(MPI_Comm comm, int varsPerNode, int order, int nx,
                              int ny, int nz) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  TACSAssembler *assembler;
  TACSCreator *creator;

  createAssembler(comm, varsPerNode, order, nx, ny, nz, &assembler, &creator);

  TACSParallelMat *mat = assembler->createMat();
  mat->incref();

  TACSBVec *x = assembler->createVec();
  TACSBVec *y_mat = assembler->createVec();
  TACSBVec *y_free = assembler->createVec();
  x->setRand(-1.0, 1.0);
  assembler->applyBCs(x);

  assembler->setVariables(x);
  if (rank == 0) {
    assembler->testElement(0, 2);
  }

  double alpha = 1.0, beta = 0.0, gamma = 0.0;

  double tassemble = MPI_Wtime();
  assembler->assembleJacobian(alpha, beta, gamma, NULL, mat);
  tassemble = MPI_Wtime() - tassemble;
  printf("Assembly time: %e\n", tassemble);

  double tprod = MPI_Wtime();
  for (int i = 0; i < 20; i++) {
    mat->mult(x, y_mat);
  }
  tprod = MPI_Wtime() - tprod;
  printf("Matrix-vector product time: %e\n", tprod);

  TACSMatrixFreeMat *mat_free = new TACSMatrixFreeMat(assembler);
  double tassemble_free = MPI_Wtime();
  mat_free->assembleMatrixFreeData(TACS_JACOBIAN_MATRIX, alpha, beta, gamma);
  tassemble_free = MPI_Wtime() - tassemble_free;
  printf("Matrix-free assembly time: %e\n", tassemble_free);

  double tprod_free = MPI_Wtime();
  for (int i = 0; i < 20; i++) {
    mat_free->mult(x, y_free);
  }
  tprod_free = MPI_Wtime() - tprod_free;
  printf("Matrix-free matrix-vector product time: %e\n", tprod_free);

  assembler->applyBCs(y_mat);
  y_mat->axpy(-1.0, y_free);
  TacsScalar norm = y_mat->norm();
  if (rank == 0) {
    printf("Residual norm of the difference: %e\n", TacsRealPart(norm));
  }

  mat->decref();

  assembler->decref();
  creator->decref();
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Set the dimension of the mesh
  int nx = 5;
  int ny = 5;
  int nz = 5;

  // Set the largest mesh order and number of variables at each node
  int order = 6;
  int varsPerNode = 4;

  // This flag indicates whether to perform a test against the
  // in-memory matrix-vector products
  int test_products = 0;

  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "varsPerNode=1") == 0) {
      varsPerNode = 1;
    } else if (strcmp(argv[i], "varsPerNode=3") == 0) {
      varsPerNode = 3;
    } else if (strcmp(argv[i], "varsPerNode=4") == 0) {
      varsPerNode = 4;
    } else if (strcmp(argv[i], "order=2") == 0) {
      order = 2;
    } else if (strcmp(argv[i], "order=3") == 0) {
      order = 3;
    } else if (strcmp(argv[i], "order=4") == 0) {
      order = 4;
    } else if (strcmp(argv[i], "order=5") == 0) {
      order = 5;
    } else if (strcmp(argv[i], "order=6") == 0) {
      order = 6;
    } else if (strcmp(argv[i], "test") == 0) {
      test_products = 1;
    }
    if (sscanf(argv[i], "nx=%d", &nx) == 0) {
      if (nx < 1) {
        nx = 1;
      }
      if (nx > 100) {
        nx = 100;
      }
    }
    if (sscanf(argv[i], "ny=%d", &ny) == 0) {
      if (ny < 1) {
        ny = 1;
      }
      if (ny > 100) {
        ny = 100;
      }
    }
    if (sscanf(argv[i], "nz=%d", &nz) == 0) {
      if (nz < 1) {
        nz = 1;
      }
      if (nz > 100) {
        nz = 100;
      }
    }
  }

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  if (test_products) {
    testMatrixVectorProducts(comm, varsPerNode, order, nx, ny, nz);
  }

  TACSAssembler *assembler[5];
  TACSCreator *creator[5];

  // Get the multigrid levels
  for (int i = 0; i < order - 1; i++) {
    int mesh_order = order - i;
    createAssembler(comm, varsPerNode, mesh_order, nx, ny, nz, &assembler[i],
                    &creator[i]);
  }

  // Create the multigrid object
  int nlevels = order - 1;
  TACSMg *mg = new TACSMg(comm, nlevels);
  mg->incref();

  // Create the interpolation operators
  for (int level = 0; level < nlevels - 1; level++) {
    // Allocate the interpolation object
    TACSBVecInterp *interp = NULL;
    interp = new TACSBVecInterp(assembler[level + 1], assembler[level]);

    int fine_order = order - level;
    int coarse_order = order - level - 1;

    if (rank == 0) {
      // Retrieve the node numbers
      const int *nodes, *coarse_nodes;
      creator[level]->getNodeNums(&nodes);
      creator[level + 1]->getNodeNums(&coarse_nodes);

      const double *fine_knots = NULL;
      if (fine_order == 3) {
        fine_knots = TacsGaussLobattoPoints3;
      } else if (fine_order == 4) {
        fine_knots = TacsGaussLobattoPoints4;
      } else if (fine_order == 5) {
        fine_knots = TacsGaussLobattoPoints5;
      } else if (fine_order == 6) {
        fine_knots = TacsGaussLobattoPoints6;
      }

      const double *coarse_knots = NULL;
      if (coarse_order == 2) {
        coarse_knots = TacsGaussLobattoPoints2;
      } else if (coarse_order == 3) {
        coarse_knots = TacsGaussLobattoPoints3;
      } else if (coarse_order == 4) {
        coarse_knots = TacsGaussLobattoPoints4;
      } else if (coarse_order == 5) {
        coarse_knots = TacsGaussLobattoPoints5;
      }

      // In this case, we cheat and set the entire interpolation
      // operator just from the root processor. The operator will be
      // communicated to the appropriate procs during initialization
      // on all procs.
      for (int nindex = 0, k = 0; k < (fine_order - 1) * nz + 1; k++) {
        for (int j = 0; j < (fine_order - 1) * ny + 1; j++) {
          for (int i = 0; i < (fine_order - 1) * nx + 1; i++, nindex++) {
            int istart, iend;
            int jstart, jend;
            int kstart, kend;
            double ni[6], nj[6], nk[6];

            // Find the indices of the coarse and fine element
            int ix = i / (fine_order - 1);
            int iy = j / (fine_order - 1);
            int iz = k / (fine_order - 1);
            if (ix >= nx) {
              ix = nx - 1;
            }
            if (iy >= ny) {
              iy = ny - 1;
            }
            if (iz >= nz) {
              iz = nz - 1;
            }

            if (i == (fine_order - 1) * nx) {
              istart = coarse_order - 1;
              iend = coarse_order;
              ni[0] = 1.0;
            } else if (i % (fine_order - 1) == 0) {
              istart = 0;
              iend = 1;
              ni[0] = 1.0;
            } else {
              istart = 0;
              iend = coarse_order;
              double u = fine_knots[i % (fine_order - 1)];
              TacsLagrangeShapeFunctions(coarse_order, u, coarse_knots, ni);
            }
            if (j == (fine_order - 1) * ny) {
              jstart = coarse_order - 1;
              jend = coarse_order;
              nj[0] = 1.0;
            } else if (j % (fine_order - 1) == 0) {
              jstart = 0;
              jend = 1;
              nj[0] = 1.0;
            } else {
              jstart = 0;
              jend = coarse_order;
              double u = fine_knots[j % (fine_order - 1)];
              TacsLagrangeShapeFunctions(coarse_order, u, coarse_knots, nj);
            }
            if (k == (fine_order - 1) * nz) {
              kstart = coarse_order - 1;
              kend = coarse_order;
              nk[0] = 1.0;
            } else if (k % (fine_order - 1) == 0) {
              kstart = 0;
              kend = 1;
              nk[0] = 1.0;
            } else {
              kstart = 0;
              kend = coarse_order;
              double u = fine_knots[k % (fine_order - 1)];
              TacsLagrangeShapeFunctions(coarse_order, u, coarse_knots, nk);
            }

            // Construct the interpolation
            int count = 0;
            TacsScalar N[216];
            int vars[216];
            for (int kk = kstart; kk < kend; kk++) {
              for (int jj = jstart; jj < jend; jj++) {
                for (int ii = istart; ii < iend; ii++) {
                  N[count] = ni[ii] * nj[jj] * nk[kk];

                  int index = (ix * (coarse_order - 1) + ii) +
                              (iy * (coarse_order - 1) + jj) *
                                  ((coarse_order - 1) * nx + 1) +
                              (iz * (coarse_order - 1) + kk) *
                                  ((coarse_order - 1) * nx + 1) *
                                  ((coarse_order - 1) * ny + 1);
                  vars[count] = coarse_nodes[index];
                  count++;
                }
              }
            }

            int node = nodes[nindex];
            interp->addInterp(node, N, vars, count);
          }
        }
      }
    }

    // Initialize the interpolation object. This is a collective
    // call that distributes the interpolation operator.
    interp->initialize();

    // Set the matrices to use
    int use_galerkin = 0;
    TACSMat *mat = NULL;
    TACSPc *smoother = NULL;

    if (fine_order > 3) {
      // Create the matrix-free matrix object
      mat = new TACSMatrixFreeMat(assembler[level]);

      // Create the smoother
      int cheb_degree = 3;
      double lower = 1.0 / 10.0, upper = 1.1;
      int smooth_iters = 1;
      smoother = new TACSChebyshevSmoother(mat, cheb_degree, lower, upper,
                                           smooth_iters);
    }

    mg->setLevel(level, assembler[level], interp, 1, use_galerkin, mat,
                 smoother);
  }

  // Set the model at the lowest grid level
  int use_galerkin = 1;
  mg->setLevel(nlevels - 1, assembler[nlevels - 1], NULL, 1, use_galerkin);

  // We no longer require any of the creator objects
  for (int i = 0; i < nlevels; i++) {
    creator[i]->decref();
  }

  // Create the residual and solution vectors on the finest TACS mesh
  TACSBVec *force = assembler[0]->createVec();
  force->incref();
  TACSBVec *res = assembler[0]->createVec();
  res->incref();
  TACSBVec *ans = assembler[0]->createVec();
  ans->incref();

  // Allocate the GMRES solution method
  int gmres_iters = 25;
  int nrestart = 8;
  int is_flexible = 0;
  GMRES *ksm = new GMRES(mg->getMat(0), mg, gmres_iters, nrestart, is_flexible);
  ksm->incref();

  // Set a monitor to check on solution progress
  int freq = 1;
  ksm->setMonitor(new KSMPrintStdout("GMRES", rank, freq));

  // The initial time
  double t0 = MPI_Wtime();

  // Assemble the Jacobian matrix for each level
  mg->assembleJacobian(1.0, 0.0, 0.0, res);

  force->zeroEntries();
  TacsScalar *force_array;
  int size = force->getArray(&force_array);
  for (int i = 1; i < size; i += varsPerNode) {
    force_array[i] = 1.0;
  }
  assembler[0]->applyBCs(force);

  // Factor the preconditioner
  mg->factor();

  // Compute the solution using GMRES
  ksm->solve(force, ans);

  t0 = MPI_Wtime() - t0;

  // Set the variables into TACS
  assembler[0]->setVariables(ans);

  // Compute the residual
  TACSMat *matrix = mg->getMat(0);
  matrix->mult(ans, res);
  res->axpy(-1.0, force);
  TacsScalar res_norm = res->norm();
  if (rank == 0) {
    printf("||R||: %15.5e\n", TacsRealPart(res_norm));
    printf("Solution time: %e\n", t0);
  }

  // Output for visualization
  ElementType etype = TACS_SOLID_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler[0], etype, write_flag);
  f5->incref();
  f5->writeToFile("volume.f5");

  // Free the memory
  f5->decref();
  ans->decref();
  res->decref();
  force->decref();
  ksm->decref();
  for (int i = 0; i < nlevels; i++) {
    assembler[i]->decref();
  }

  MPI_Finalize();
  return (0);
}
