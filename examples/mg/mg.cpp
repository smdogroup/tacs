#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElement2D.h"
#include "TACSLinearElasticity.h"
#include "TACSMg.h"
#include "TACSQuadBasis.h"
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
  int varsPerNode = 2;

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
    int numBcs = 4 * (nx + 1);
    int *bcNodes = new int[numBcs];

    for (int i = 0; i < (nx + 1); i++) {
      bcNodes[4 * i] = i;
      bcNodes[4 * i + 1] = i + (nx + 1) * ny;
      bcNodes[4 * i + 2] = i * (nx + 1);
      bcNodes[4 * i + 3] = (i + 1) * (nx + 1) - 1;
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
  TACSPlaneStressConstitutive *stiff1 = new TACSPlaneStressConstitutive(props1);
  TACSPlaneStressConstitutive *stiff2 = new TACSPlaneStressConstitutive(props2);

  // Create the model class
  TACSLinearElasticity2D *model1 =
      new TACSLinearElasticity2D(stiff1, TACS_LINEAR_STRAIN);
  TACSLinearElasticity2D *model2 =
      new TACSLinearElasticity2D(stiff2, TACS_LINEAR_STRAIN);

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
  capabilities in TACS using a flat plate example.

  The main work required is the assembly of the inter-grid operators
  themselves. The TACSBVecInterp class is used to simplify the set up
  of these operators which are stored explicitly.

  The multigrid operations are performed using the TACSMg object which
  inherits from the TACSPc interface. It can be used directly as a
  preconditioner and needs to be "factored" in the same manner as
  regular preconditioners.

  TACSMg has a number of member functions which facilitate setting
  state and design variables on all multigrid levels and matrix
  assembly.  Matrix assembly on all levels can be performed by calling
  the assembleJacobian() and the assembleMatType() functions directly
  in the TACSMg interface.
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
  }

  // Create the multigrid object
  TACSMg *mg = new TACSMg(comm, nlevels);
  mg->incref();

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

  double tmg = MPI_Wtime();

  // Create the matrix for the finest grid level
  TACSParallelMat *mat = assembler[0]->createMat();
  mat->incref();

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

    // Set the multigrid information at this level
    double tlev = MPI_Wtime();
    int use_galerkin = 1;
    mg->setLevel(level, assembler[level], interp[level], 1, use_galerkin);
    tlev = MPI_Wtime() - tlev;
    if (rank == 0) {
      printf("Initialization time for level %d: %e\n", level, tlev);
    }
  }

  // Set the model at the lowest grid level
  int use_galerkin = 1;
  mg->setLevel(nlevels - 1, assembler[nlevels - 1], NULL, 1, use_galerkin);

  // We no longer require any of the creator objects
  for (int i = 0; i < nlevels; i++) {
    creator[i]->decref();
  }

  tmg = MPI_Wtime() - tmg;
  if (rank == 0) {
    printf("TACSMg creation time: %e\n", tmg);
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
  for (int i = 1; i < size; i += 2) {
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
  ElementType etype = TACS_PLANE_STRESS_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler[0], etype, write_flag);
  f5->incref();
  f5->writeToFile("plate.f5");

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
