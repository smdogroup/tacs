#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElement2D.h"
#include "TACSHeatConduction.h"
#include "TACSKSTemperature.h"
#include "TACSMg.h"
#include "TACSQuadBasis.h"
#include "TACSSpectralIntegrator.h"
#include "TACSToFH5.h"

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int nx, int ny, TACSAssembler **_assembler,
                     TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the number of nodes/elements on this proc
  int varsPerNode = 1;

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0) {
    // Set the number of elements
    int numNodes = (nx + 1) * (ny + 1);
    int numElements = nx * ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[4 * numElements];

    // Set the ids as alternating
    memset(ids, 0, numElements * sizeof(int));
    for (int k = 0; k < numElements; k++) {
      if (k % nx < nx / 2 && (k / nx) < nx / 2) {
        ids[k] = 1;
      } else if (k % nx >= nx / 2 && (k / nx) >= nx / 2) {
        ids[k] = 1;
      }
    }

    ptr[0] = 0;
    for (int k = 0; k < numElements; k++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k / nx;

      // Set the node connectivity
      conn[4 * k] = i + j * (nx + 1);
      conn[4 * k + 1] = i + 1 + j * (nx + 1);
      conn[4 * k + 2] = i + (j + 1) * (nx + 1);
      conn[4 * k + 3] = i + 1 + (j + 1) * (nx + 1);
      ptr[k + 1] = 4 * (k + 1);
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    // We're over-counting one of the nodes on each edge
    int numBcs = (nx + 1);
    int *bcNodes = new int[numBcs];

    for (int i = 0; i < (nx + 1); i++) {
      bcNodes[i] = i;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete[] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];
    for (int j = 0; j < ny + 1; j++) {
      for (int i = 0; i < nx + 1; i++) {
        int node = i + (nx + 1) * j;
        Xpts[3 * node] = (1.0 * i) / nx;
        Xpts[3 * node + 1] = (1.0 * j) / ny;
        Xpts[3 * node + 2] = 0.0;
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
  TacsScalar t1 = 1.0, t2 = 1.0;
  int t1Num = 0, t2Num = 1;
  TACSPlaneStressConstitutive *stiff1 =
      new TACSPlaneStressConstitutive(props1, t1, t1Num);
  TACSPlaneStressConstitutive *stiff2 =
      new TACSPlaneStressConstitutive(props2, t2, t2Num);

  // Create the model class
  TACSHeatConduction2D *model1 = new TACSHeatConduction2D(stiff1);
  TACSHeatConduction2D *model2 = new TACSHeatConduction2D(stiff2);

  // Create the element class
  TACSElementBasis *linear_basis = new TACSLinearQuadBasis();
  TACSElement2D *linear_element1 = new TACSElement2D(model1, linear_basis);
  TACSElement2D *linear_element2 = new TACSElement2D(model2, linear_basis);

  // Set the one element
  TACSElement *elems[2];
  elems[0] = linear_element1;
  elems[1] = linear_element2;
  creator->setElements(2, elems);

  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

/*
  The following code illustrates the use of the geometric multigrid
  capabilities in TACS combined with a spectral time integration technique.
*/
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Number of different levels
  int nlevels = 3;
  const int max_nlevels = 8;
  TACSAssembler *assembler[max_nlevels];
  TACSCreator *creator[max_nlevels];

  // Set the dimension of the largest meshes
  int nx = 128;
  int ny = 128;

  // Use GMRES (if not, use GCROT)
  int use_gmres = 1;

  // Get the global size of the mesh from the input
  for (int k = 0; k < argc; k++) {
    int xpow;
    if (sscanf(argv[k], "xpow=%d", &xpow) == 1) {
      if (xpow < 5) {
        xpow = 5;
      }
      if (xpow > 10) {
        xpow = 10;
      }
      nx = 1 << xpow;
      ny = 1 << xpow;
    }
    if (sscanf(argv[k], "nlevels=%d", &nlevels) == 1) {
      if (nlevels < 2) {
        nlevels = 2;
      }
      if (nlevels > max_nlevels) {
        nlevels = max_nlevels;
      }
    }
    if (strcmp("use_gcrot", argv[k]) == 0) {
      use_gmres = 0;
    }
  }

  // Create the TACS/Creator objects for all levels
  for (int i = 0; i < nlevels; i++) {
    double t0 = MPI_Wtime();
    int Nx = nx / (1 << i), Ny = ny / (1 << i);
    createAssembler(comm, Nx, Ny, &assembler[i], &creator[i]);
    t0 = MPI_Wtime() - t0;
    assembler[i]->incref();
    creator[i]->incref();
    if (rank == 0) {
      printf("Assembler creation time for level %d: %e\n", i, t0);
    }
  }

  // Allocate the interpolation objects for all remaining levels
  TACSBVecInterp *interp[max_nlevels - 1];

  // Create the interpolation operators
  for (int level = 0; level < nlevels - 1; level++) {
    // Allocate the interpolation object
    interp[level] = new TACSBVecInterp(assembler[level + 1], assembler[level]);

    if (rank == 0) {
      // Retrieve the node numbers
      const int *nodes, *coarse_nodes;
      creator[level]->getNodeNums(&nodes);
      creator[level + 1]->getNodeNums(&coarse_nodes);

      // Set the weights and variables for
      // each node in the fine mesh
      int vars[4];
      TacsScalar w[4];

      // Compute the number of nodes on the fine and coarse
      // meshes
      int nx_fine = nx / (1 << level);
      int ny_fine = ny / (1 << level);
      int nx_coarse = nx / (1 << (level + 1));

      // In this case, we cheat and set the entire interpolation
      // operator just from the root processor. The operator will be
      // communicated to the appropriate procs during initialization
      // on all procs.
      for (int j = 0; j < ny_fine + 1; j++) {
        for (int i = 0; i < nx_fine + 1; i++) {
          int nvars = 0;
          if ((i % 2 == 0) && (j % 2 == 0)) {
            w[0] = 1.0;
            vars[0] = coarse_nodes[(i / 2) + (j / 2) * (nx_coarse + 1)];
            nvars = 1;
          } else if (i % 2 == 0) {
            w[0] = w[1] = 0.5;
            vars[0] = coarse_nodes[(i / 2) + (j / 2) * (nx_coarse + 1)];
            vars[1] = coarse_nodes[(i / 2) + (j / 2 + 1) * (nx_coarse + 1)];
            nvars = 2;
          } else if (j % 2 == 0) {
            w[0] = w[1] = 0.5;
            vars[0] = coarse_nodes[(i / 2) + (j / 2) * (nx_coarse + 1)];
            vars[1] = coarse_nodes[(i / 2 + 1) + (j / 2) * (nx_coarse + 1)];
            nvars = 2;
          } else {
            w[0] = w[1] = w[2] = w[3] = 0.25;
            vars[0] = coarse_nodes[(i / 2) + (j / 2) * (nx_coarse + 1)];
            vars[1] = coarse_nodes[(i / 2 + 1) + (j / 2) * (nx_coarse + 1)];
            vars[2] = coarse_nodes[(i / 2) + (j / 2 + 1) * (nx_coarse + 1)];
            vars[3] = coarse_nodes[(i / 2 + 1) + (j / 2 + 1) * (nx_coarse + 1)];
            nvars = 4;
          }

          int node = nodes[i + (nx_fine + 1) * j];
          interp[level]->addInterp(node, w, vars, nvars);
        }
      }
    }

    // Initialize the interpolation object. This is a collective
    // call that distributes the interpolation operator.
    interp[level]->initialize();
  }

  // Now create the spectral integrator class
  double tfinal = 300.0;  // Final time for the simulation
  int N = 32;
  TACSSpectralIntegrator *spectral =
      new TACSSpectralIntegrator(assembler[0], tfinal, N);
  spectral->incref();

  // Create the matrix
  TACSLinearSpectralMat *mat = spectral->createLinearMat();
  mat->incref();

  // Create the multigrid preconditioner for the spectral problem
  int coarsen_time[max_nlevels - 1];

  for (int i = 0; i < nlevels - 1; i++) {
    coarsen_time[i] = 0;
  }

  // Only coarsen the time twice and not at the first level
  for (int i = 3; i <= 4 && i < nlevels - 1; i++) {
    coarsen_time[i] = 1;
  }

  TACSLinearSpectralMg *mg =
      new TACSLinearSpectralMg(mat, nlevels, assembler, interp, coarsen_time);
  mg->incref();

  // Form the spectral problem and factor it
  spectral->assembleMat(mat);
  mg->factor();

  // Solve the problem
  TACSSpectralVec *res = spectral->createVec();
  res->incref();
  TACSSpectralVec *rhs = spectral->createVec();
  rhs->incref();
  TACSSpectralVec *ans = spectral->createVec();
  ans->incref();

  // Set the values of the right-hand-side
  for (int i = N / 2; i < N; i++) {
    rhs->getVec(i)->set(100.0);
    assembler[0]->applyBCs(rhs->getVec(i));
  }

  TACSKsm *ksm = NULL;

  if (use_gmres) {
    // Allocate the GMRES solution method
    int gmres_subspace = 50;
    int nrestart = 8;
    int is_flexible = 0;
    ksm = new GMRES(mat, mg, gmres_subspace, nrestart, is_flexible);
  } else {
    int gmres_subspace = 25;
    int outer = 10;
    int max_outer = 5 * outer;
    int is_flexible = 1;
    ksm = new GCROT(mat, mg, outer, max_outer, gmres_subspace, is_flexible);
  }
  ksm->incref();

  // Set a monitor to check on solution progress
  int freq = 5;
  ksm->setMonitor(new KSMPrintStdout("GMRES", rank, freq));
  double rtol = 1e-12, atol = 1e-30;
  ksm->setTolerances(rtol, atol);

  // Compute the solution using GMRES
  ksm->solve(rhs, ans);

  // Set the variables into TACS
  spectral->setVariables(ans);
  spectral->assembleRes(res);
  res->axpy(-1.0, rhs);

  // Evaluate a function of interest
  double ksweight = 100;
  TACSFunction *func = new TACSKSTemperature(assembler[0], ksweight);
  func->incref();

  TacsScalar fval;
  spectral->evalFunctions(1, &func, &fval);

  TacsScalar res_norm = res->norm();
  if (rank == 0) {
    printf("Maximum temperature: %25.15e\n", TacsRealPart(fval));
    printf("||R||: %25.15e\n", TacsRealPart(res_norm));
  }

  TACSSpectralVec *dfdu = spectral->createVec();
  dfdu->incref();

  // Assemble the system of equations
  spectral->evalSVSens(func, dfdu);

  // Assemble the transpose matrix and factor it
  spectral->assembleMat(mat, TACS_MAT_TRANSPOSE);
  mg->factor();

  // Solve for the adjoint variables
  ksm->solve(dfdu, ans);

  TACSBVec *dfdx = assembler[0]->createDesignVec();
  dfdx->incref();

  spectral->addAdjointResProduct(-1.0, ans, dfdx);

  dfdx->beginSetValues(TACS_ADD_VALUES);
  dfdx->endSetValues(TACS_ADD_VALUES);

  TacsScalar *dfdx_array;
  int len = dfdx->getArray(&dfdx_array);
  for (int i = 0; i < len; i++) {
    printf("dfdx[%d] = %25.15e\n", i, TacsRealPart(dfdx_array[i]));
  }

  TACSBVec *x = assembler[0]->createDesignVec();
  assembler[0]->getDesignVars(x);
  TACSBVec *px = assembler[0]->createDesignVec();
  px->setRand(-1.0, 1.0);
  double dh = 1e-6;
  x->axpy(dh, px);

  // Reset the design variable values
  assembler[0]->setDesignVars(x);

  // Assemble the matrix and factor it
  spectral->assembleMat(mat);
  mg->factor();

  // Compute the solution using GMRES
  ksm->solve(rhs, ans);

  // Set the variables into TACS
  spectral->setVariables(ans);

  TacsScalar fval1;
  spectral->evalFunctions(1, &func, &fval1);

  TacsScalar fd = (fval1 - fval) / dh;
  TacsScalar dot = dfdx->dot(px);

  if (rank == 0) {
    printf("FD:   %25.15e\n", TacsRealPart(fd));
    printf("An:   %25.15e\n", TacsRealPart(dot));
    printf("Err:  %25.15e\n", TacsRealPart((fd - dot) / fd));
  }

  dfdx->decref();
  dfdu->decref();

  // Output for visualization
  ElementType etype = TACS_PLANE_STRESS_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler[0], etype, write_flag);
  f5->incref();

  TACSBVec *u = assembler[0]->createVec();
  u->incref();

  int nviz = 100;
  for (int i = 0; i < nviz; i++) {
    char filename[256];
    snprintf(filename, sizeof(filename), "plate%d.f5", i);

    double time = i * tfinal / (nviz - 1);
    spectral->computeSolutionAndDeriv(time, ans, u);
    assembler[0]->setVariables(u);
    f5->writeToFile(filename);
  }

  u->decref();

  // Free the memory
  f5->decref();
  ans->decref();
  res->decref();
  ksm->decref();
  for (int i = 0; i < nlevels; i++) {
    assembler[i]->decref();
  }
  spectral->decref();

  MPI_Finalize();
  return (0);
}
