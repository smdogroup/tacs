#include "TACSAssembler.h"
#include "TACSCreator.h"
#include "TACSElementAlgebra.h"
#include "TACSElementVerification.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSShellElementDefs.h"
#include "TACSShellTraction.h"
#include "TACSToFH5.h"
#include "tacslapack.h"

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler(MPI_Comm comm, int nx, int ny, TACSElement *element,
                     TACSAssembler **_assembler, TACSCreator **_creator) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the number of nodes/elements on this proc
  int varsPerNode = element->getVarsPerNode();

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0) {
    // Set the number of elements
    int numNodes = (2 * nx + 1) * (2 * ny + 1);
    int numElements = nx * ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[numElements];
    int *ptr = new int[numElements + 1];
    int *conn = new int[9 * numElements];

    // Set the element identifiers to all zero
    memset(ids, 0, numElements * sizeof(int));

    ptr[0] = 0;
    for (int k = 0; k < numElements; k++) {
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k / nx;

      // Set the node connectivity
      for (int jj = 0; jj < 3; jj++) {
        for (int ii = 0; ii < 3; ii++) {
          conn[9 * k + ii + 3 * jj] = 2 * i + ii + (2 * j + jj) * (2 * nx + 1);
          ptr[k + 1] = 9 * (k + 1);
        }
      }
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements, ptr, conn, ids);
    delete[] conn;
    delete[] ptr;
    delete[] ids;

    // We're over-counting one of the nodes on each edge
    int numBcs = 4 * (2 * nx + 1);
    int *bcNodes = new int[numBcs];

    for (int i = 0; i < (2 * nx + 1); i++) {
      bcNodes[4 * i] = i;
      bcNodes[4 * i + 1] = i + (2 * nx + 1) * (2 * ny);
      bcNodes[4 * i + 2] = i * (2 * nx + 1);
      bcNodes[4 * i + 3] = (i + 1) * (2 * nx + 1) - 1;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete[] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[3 * numNodes];
    for (int j = 0; j < 2 * ny + 1; j++) {
      for (int i = 0; i < 2 * nx + 1; i++) {
        int node = i + (2 * nx + 1) * j;
        TacsScalar x = (0.5 * i) / nx;
        TacsScalar y = (0.5 * j) / ny;
        TacsScalar z = x * (1.0 - x) * y * (1.0 - y);
        Xpts[3 * node] = x;
        Xpts[3 * node + 1] = y;
        Xpts[3 * node + 2] = z;
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete[] Xpts;
  }

  // Set the one element
  creator->setElements(1, &element);

  // Set the reordering type
  creator->setReorderingType(TACSAssembler::MULTICOLOR_ORDER,
                             TACSAssembler::GAUSS_SEIDEL);

  // Create TACS
  TACSAssembler *assembler = creator->createTACS();

  // Set the pointers
  *_assembler = assembler;
  *_creator = creator;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TacsScalar rho = 2700.0;
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e3;
  TacsScalar nu = 0.3;
  TacsScalar ys = 270.0;
  TacsScalar cte = 24.0e-6;
  TacsScalar kappa = 230.0;
  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TacsScalar axis[] = {0.0, 1.0, 0.0};
  TACSShellTransform *transform = new TACSShellRefAxisTransform(axis);

  TacsScalar t = 0.01;
  int t_num = 0;
  TACSShellConstitutive *con = new TACSIsoShellConstitutive(props, t, t_num);

  TACSElement *shell9 = new TACSQuad9Shell(transform, con);
  shell9->incref();

  int nx = 20, ny = 20;
  TACSAssembler *assembler;
  TACSCreator *creator;
  createAssembler(comm, nx, ny, shell9, &assembler, &creator);
  assembler->incref();
  creator->incref();

  // Free the creator object
  creator->decref();

  // Create matrix and vectors
  TACSBVec *ans = assembler->createVec();  // displacements and rotations
  TACSBVec *f = assembler->createVec();    // loads
  TACSBVec *res = assembler->createVec();  // The residual
  TACSSchurMat *mat = assembler->createSchurMat();  // stiffness matrix

  // Increment reference count to the matrix/vectors
  ans->incref();
  f->incref();
  res->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 10000;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
  pc->incref();

  // Set all the entries in load vector to specified value
  TacsScalar *force_vals;
  int size = f->getArray(&force_vals);
  for (int k = 2; k < size; k += assembler->getVarsPerNode()) {
    force_vals[k] += 100.0;
  }
  assembler->applyBCs(f);

  // Assemble and factor the stiffness/Jacobian matrix. Factor the
  // Jacobian and solve the linear system for the displacements
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  assembler->assembleJacobian(alpha, beta, gamma, res, mat);
  pc->factor();  // LU factorization of stiffness matrix

  res->axpy(-1.0, f);  // Compute res - f
  pc->applyFactor(res, ans);

  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Output for visualization
  ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;
  int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_CONNECTIVITY |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("plate_solution.f5");
  assembler->decref();

  shell9->decref();

  MPI_Finalize();
}