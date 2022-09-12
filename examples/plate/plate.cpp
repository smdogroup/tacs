/*
  Pressure loaded plate
*/
#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSToFH5.h"

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the number of nodes/elements on this proc
  int varsPerNode = 6;

  // Set the number of elements in the x/y direction
  int order = 3;
  int nx = 25, ny = 25;

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);
  creator->incref();

  if (rank == 0) {
    int numNodes = ((order - 1) * nx + 1) * ((order - 1) * ny + 1);
    int numElements = nx * ny;

    // The elements are ordered as (i + j*nx)
    int *ptr = new int[numElements + 1];
    int *conn = new int[order * order * numElements];
    int *ids = new int[numElements];
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    int *c = conn;
    for (int j = 0, k = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++, k++) {
        // Set the node connectivity
        for (int jj = 0, kk = 0; jj < order; jj++) {
          for (int ii = 0; ii < order; ii++, kk++) {
            c[kk] = ((order - 1) * i + ii) +
                    ((order - 1) * j + jj) * ((order - 1) * nx + 1);
          }
        }
        c += order * order;
        ptr[k + 1] = order * order * (k + 1);
      }
    }

    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);

    // We're over-counting one of the nodes on each edge
    int numBcs = 4 * ((order - 1) * nx + 1);
    int *bcNodes = new int[numBcs];

    for (int i = 0; i < ((order - 1) * nx + 1); i++) {
      bcNodes[4 * i] = i;
      bcNodes[4 * i + 1] = i + ((order - 1) * nx + 1) * (order - 1) * ny;
      bcNodes[4 * i + 2] = i * ((order - 1) * nx + 1);
      bcNodes[4 * i + 3] = (i + 1) * ((order - 1) * nx + 1) - 1;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete[] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];
    for (int j = 0; j < (order - 1) * ny + 1; j++) {
      for (int i = 0; i < (order - 1) * nx + 1; i++) {
        int node = i + ((order - 1) * nx + 1) * j;
        Xpts[3 * node] = (1.0 * i) / ((order - 1) * nx);
        Xpts[3 * node + 1] = (1.0 * j) / ((order - 1) * ny);
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
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {1.0, 0.0, 0.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TacsScalar t = 0.01;
  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t);

  TACSElement *element = NULL;
  if (order == 2) {
    element = new TACSQuad4Shell(transform, con);

  } else {  // order == 3
    element = new TACSQuad9Shell(transform, con);
  }

  // Set the elements
  creator->setElements(1, &element);

  // Create the TACSAssembler object
  TACSAssembler *assembler = creator->createTACS();
  creator->decref();
  assembler->incref();

  // Create the preconditioner
  TACSBVec *res = assembler->createVec();
  TACSBVec *ans = assembler->createVec();

  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();

  // Set the parameters for the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;

  // Allocate, assemble and factor the TACSSchurMat and preconditioner
  TACSSchurMat *schur_mat = assembler->createSchurMat();
  schur_mat->incref();

  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *schur_pc = new TACSSchurPc(schur_mat, lev, fill, reorder_schur);
  schur_pc->incref();

  double schur_time = MPI_Wtime();
  assembler->assembleJacobian(alpha, beta, gamma, res, schur_mat);
  schur_pc->factor();
  schur_time = MPI_Wtime() - schur_time;

  if (rank == 0) {
    printf("TACSSchurMat/TACSSchurPc assembly and factorization time: %15.5e\n",
           schur_time);
  }

  // Allocate the GMRES object with the TACSSchurMat matrix
  int gmres_iters = 10;  // Number of GMRES iterations
  int nrestart = 2;      // Number of allowed restarts
  int is_flexible = 1;   // Is a flexible preconditioner?
  GMRES *schur_gmres =
      new GMRES(schur_mat, schur_pc, gmres_iters, nrestart, is_flexible);
  TacsScalar *res_array;
  int size = res->getArray(&res_array);
  for (int i = 2; i < size; i += varsPerNode) {
    res_array[i] = 1.0;
  }

  assembler->applyBCs(res);
  schur_gmres->solve(res, ans);
  assembler->setVariables(ans);

  // Create an TACSToFH5 object for writing output to files
  ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("plate.f5");
  f5->decref();

  assembler->decref();

  MPI_Finalize();
  return 0;
}
