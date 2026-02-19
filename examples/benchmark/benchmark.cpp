#include "TACSAssembler.h"
#include "TACSCompliance.h"
#include "TACSElement2D.h"
#include "TACSLinearElasticity.h"
#include "TACSQuadBasis.h"
#include "TACSTraction2D.h"

/*
  Create the material properties object
*/
void createMaterialProperties(TACSMaterialProperties **props) {
  TacsScalar rho = 1.0;
  TacsScalar specific_heat = 0.0;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 1.0;
  TacsScalar alpha = 0.0;
  TacsScalar kappa = 0.0;

  *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, alpha, kappa);
}

/*
  Create a plane stress stiffness object with default properties
*/
void createPlaneStressElement(TACSMaterialProperties *props, int order, int num,
                              TACSElement **elem) {
  TACSPlaneStressConstitutive *stiff =
      new TACSPlaneStressConstitutive(props, 1.0, num);

  TACSLinearElasticity2D *model =
      new TACSLinearElasticity2D(stiff, TACS_LINEAR_STRAIN);

  if (order == 2) {
    TACSElementBasis *basis = new TACSLinearQuadBasis();
    *elem = new TACSElement2D(model, basis);
  } else if (order == 3) {
    TACSElementBasis *basis = new TACSQuadraticQuadBasis();
    *elem = new TACSElement2D(model, basis);
  } else {  // if (order == 4){
    TACSElementBasis *basis = new TACSCubicQuadBasis();
    *elem = new TACSElement2D(model, basis);
  }
}

/*
  Create a 2D TACS model that is either a plane stress or shell model

  input:
  comm:         the MPI communicator
  varsPerNode:  the number of unknowns per node
  nx:           number of elements in the x-direction
  ny:           number of elements in the y-direction
  firstNode:    the first owned node number on this partition
  lastNode:     the last node owned by this partition
  firstElem:    the first element index
  lastElem:     the last element index
  noptions:     the number of options
  opts:         the option values
*/
TACSAssembler *create2DModel(MPI_Comm comm, int varsPerNode, int nx, int ny,
                             int order, int firstNode, int lastNode,
                             int firstElem, int lastElem, int noptions,
                             const char *opts[]) {
  // Set the number of nodes
  int numOwnedNodes = lastNode - firstNode;
  int numElements = lastElem - firstElem;

  // There are no dependent nodes in this problem
  TACSAssembler *assembler =
      new TACSAssembler(comm, varsPerNode, numOwnedNodes, numElements);

  // The elements are ordered as (i + j*nx)
  int *ptr = new int[numElements + 1];
  int *conn = new int[order * order * numElements];

  int *c = conn;
  ptr[0] = 0;
  for (int k = 0, elem = firstElem; elem < lastElem; k++, elem++) {
    // Back out the i, j coordinates from the corresponding
    // element number
    int i = elem % nx;
    int j = elem / nx;

    // Set the node connectivity
    for (int jj = 0; jj < order; jj++) {
      for (int ii = 0; ii < order; ii++) {
        c[0] = (order - 1) * i + ii +
               ((order - 1) * nx + 1) * ((order - 1) * j + jj);
        c++;
      }
    }
    ptr[k + 1] = c - conn;
  }

  // Set the connectivity
  assembler->setElementConnectivity(ptr, conn);
  delete[] conn;
  delete[] ptr;

  // Create and set the elements
  TACSElement **elements = new TACSElement *[numElements];

  TACSMaterialProperties *props;
  createMaterialProperties(&props);
  for (int k = 0, elem = firstElem; elem < lastElem; k++, elem++) {
    // Create a plane stress model
    createPlaneStressElement(props, order, elem, &elements[k]);
  }

  // Set the elements into the mesh
  assembler->setElements(elements);
  delete[] elements;

  // Set the boundary conditions - this will only record the
  // boundary conditions on its own nodes
  for (int i = 0; i < (order - 1) * nx + 1; i++) {
    int nodes[4];
    nodes[0] = i;
    nodes[1] = i + ((order - 1) * nx + 1) * ((order - 1) * ny);
    nodes[2] = i * ((order - 1) * nx + 1);
    nodes[3] = (i + 1) * ((order - 1) * nx + 1) - 1;
    assembler->addBCs(4, nodes);
  }

  // Reorder the nodal variables
  int reorder = 0;
  enum TACSAssembler::OrderingType order_type = TACSAssembler::ND_ORDER;
  enum TACSAssembler::MatrixOrderingType mat_type =
      TACSAssembler::APPROXIMATE_SCHUR;

  for (int k = 0; k < noptions; k++) {
    if (strcmp(opts[k], "AMD") == 0) {
      order_type = TACSAssembler::AMD_ORDER;
      reorder = 1;
    } else if (strcmp(opts[k], "RCM") == 0) {
      order_type = TACSAssembler::RCM_ORDER;
      reorder = 1;
    } else if (strcmp(opts[k], "ND") == 0) {
      order_type = TACSAssembler::ND_ORDER;
      reorder = 1;
    } else if (strcmp(opts[k], "TACS_AMD") == 0) {
      order_type = TACSAssembler::TACS_AMD_ORDER;
      reorder = 1;
    } else if (strcmp(opts[k], "DirectSchur") == 0) {
      mat_type = TACSAssembler::DIRECT_SCHUR;
      reorder = 1;
    } else if (strcmp(opts[k], "ApproximateSchur") == 0) {
      mat_type = TACSAssembler::APPROXIMATE_SCHUR;
      reorder = 1;
    } else if (strcmp(opts[k], "AdditiveSchwarz") == 0) {
      mat_type = TACSAssembler::ADDITIVE_SCHWARZ;
      reorder = 1;
    }
  }

  // If we're going to use the FEMat class, then we don't need to
  // perform a reordering
  if (reorder) {
    assembler->computeReordering(order_type, mat_type);
  }

  // Perform initialization - cannot add any more elements/vars etc
  assembler->initialize();

  // Create the node vector
  TACSBVec *X = assembler->createNodeVec();
  X->incref();

  // Get the local node locations
  TacsScalar *Xpts = NULL;
  X->getArray(&Xpts);
  for (int k = 0, node = firstNode; node < lastNode; k += 3, node++) {
    int i = node % ((order - 1) * nx + 1);
    int j = node / ((order - 1) * nx + 1);
    Xpts[k] = 250.0 * i / ((order - 1) * nx);
    Xpts[k + 1] = 250.0 * j / ((order - 1) * ny);
  }

  // Reorder the vector if required
  assembler->reorderVec(X);

  // Set the node locations
  assembler->setNodes(X);
  X->decref();

  return assembler;
}

/*
  Solve the problem with the specified options
*/
void testSolve(TACSAssembler *assembler, int noptions, const char *opts[]) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Create the preconditioner
  int isflexible = 1;
  TACSMat *mat = NULL;
  TACSPc *pc = NULL;

  // Set the default options
  double fill = 10.0;
  int levFill = 5;
  int inner_gmres = 10;
  double sor_omega = 1.0;
  int sor_iters = 2;
  int sor_symmetric = 1;
  int gmres_iters = 50;

  // Check if any of the options are set
  for (int k = 0; k < noptions; k++) {
    if (sscanf(opts[k], "levFill=%d", &levFill) == 1) {
      if (levFill < 0) {
        levFill = 0;
      }
    }
    if (sscanf(opts[k], "inner_gmres=%d", &inner_gmres) == 1) {
      if (inner_gmres < 1) {
        inner_gmres = 1;
      }
    }
    if (sscanf(opts[k], "gmres_iters=%d", &gmres_iters) == 1) {
      if (gmres_iters < 1) {
        gmres_iters = 1;
      }
    }
  }

  // Create the matrix and associated preconditioner
  for (int k = 0; k < noptions; k++) {
    if (strcmp(opts[k], "ApproximateSchur") == 0) {
      TACSParallelMat *_mat = assembler->createMat();
      double inner_rtol = -1.0;
      double inner_atol = 1e-10;
      pc = new TACSApproximateSchur(_mat, levFill, fill, inner_gmres,
                                    inner_rtol, inner_atol);
      mat = _mat;
      break;
    } else if (strcmp(opts[k], "DirectSchur") == 0) {
      TACSSchurMat *_mat = assembler->createSchurMat();
      int reorder_schur = 1;
      pc = new TACSSchurPc(_mat, levFill, fill, reorder_schur);
      mat = _mat;
      break;
    } else if (strcmp(opts[k], "GaussSeidel") == 0) {
      int zero_guess = 0;  // Zero the initial guess for psor
      TACSParallelMat *_mat = assembler->createMat();
      pc = new TACSGaussSeidel(_mat, zero_guess, sor_omega, sor_iters,
                               sor_symmetric);
      mat = _mat;
    }
  }

  // Create the additive Schwarz preconditioner
  if (!pc || !mat) {
    isflexible = 0;
    TACSParallelMat *_mat = assembler->createMat();
    pc = new TACSAdditiveSchwarz(_mat, levFill, fill);
    mat = _mat;
  }
  mat->incref();
  pc->incref();

  // Create the vectors to be used in the solution
  TACSBVec *ans = assembler->createVec();
  TACSBVec *rhs = assembler->createVec();
  ans->incref();
  rhs->incref();
  assembler->zeroVariables();

  double start = MPI_Wtime();
  assembler->assembleJacobian(1.0, 0.0, 0.0, rhs, mat);
  double stop = MPI_Wtime();
  if (rank == 0) {
    printf("Matrix assembly time: %10.4f\n", stop - start);
  }

  // Set a right-hand-side
  rhs->set(1.0);
  assembler->applyBCs(rhs);

  TacsZeroNumFlops();
  start = MPI_Wtime();
  pc->factor();
  stop = MPI_Wtime();
  double flops = TacsGetNumFlops();
  if (rank == 0) {
    printf("Factor time: %10.4f\n", stop - start);
  }

  printf("[%d] FLOPS: %15.5e FLOP rate: %15.5e\n", rank, flops,
         flops / (stop - start));

  // Set up the problem
  GMRES *solver = new GMRES(mat, pc, gmres_iters, 2, isflexible);
  solver->incref();
  solver->setTolerances(1e-10, 1e-30);
  solver->setMonitor(new KSMPrintStdout(" Iteration", rank, 5));

  start = MPI_Wtime();
  solver->solve(rhs, ans);
  stop = MPI_Wtime();
  if (rank == 0) {
    printf("Solve time: %10.4f\n", stop - start);
  }

  solver->decref();
  ans->scale(-1.0);
  assembler->setVariables(ans);

  assembler->assembleRes(rhs);
  TacsScalar rhs_norm = rhs->norm();
  if (rank == 0) {
    printf("Residual norm: %10.4e\n", TacsRealPart(rhs_norm));
  }

  // Evaluate the compliance
  TACSFunction *comp = new TACSCompliance(assembler);
  comp->incref();
  TacsScalar compVal = 0.0;
  assembler->evalFunctions(1, &comp, &compVal);
  if (rank == 0) {
    printf("Compliance: %25.12f\n", TacsRealPart(compVal));
  }
  comp->decref();

  pc->decref();
  mat->decref();
  ans->decref();
  rhs->decref();
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  // Get the options
  int size, rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Set up the TACSAssembler model
  int varsPerNode = 2;
  int order = 3;

  // Number of elements in the x/y directions
  int nx = 75;
  int ny = 75;

  // Retrieve the options
  int noptions = 6;
  const char *opts[] = {"AMD",   "DirectSchur", "nx=50",
                        "ny=50", "order=3",     "levFill=1000"};

  for (int k = 0; k < noptions; k++) {
    if (sscanf(opts[k], "nx=%d", &nx) == 1) {
    }
    if (sscanf(opts[k], "ny=%d", &ny) == 1) {
    }
    if (sscanf(opts[k], "order=%d", &order) == 1) {
      if (order < 2) {
        order = 2;
      }
      if (order > 4) {
        order = 4;
      }
    }
  }

  // Determine the total number of nodes/elements
  int nelems = nx * ny;
  int nnodes = (((order - 1) * nx + 1) * ((order - 1) * ny + 1));

  // Determin the partition of the nodes (this won't be very good...)
  int elemsPerProc = nelems / size;
  int nodesPerProc = nnodes / size;
  int firstNode = rank * nodesPerProc;
  int lastNode = (rank + 1) * nodesPerProc;
  int firstElem = rank * elemsPerProc;
  int lastElem = (rank + 1) * elemsPerProc;
  if (rank == size - 1) {
    lastNode = nnodes;
    lastElem = nelems;
  }

  // Create the TACSAssembler object
  TACSAssembler *assembler =
      create2DModel(comm, varsPerNode, nx, ny, order, firstNode, lastNode,
                    firstElem, lastElem, noptions, opts);
  assembler->incref();

  // Test solve the first level for now
  testSolve(assembler, noptions, opts);

  assembler->decref();
  MPI_Finalize();

  return (0);
}
