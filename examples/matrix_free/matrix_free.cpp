/*
  This code tests the implementation of higher-order elements with
  the matrix-free matrix matrix-vector product implementation.
*/

#include "TACSLinearElasticity.h"
#include "TACSQuadBasis.h"
#include "TACSElement2D.h"

#include "TACSAssembler.h"
#include "TACSMatrixFreeMat.h"
#include "TACSMg.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"

/*
  Create the TACSAssembler object and return the associated TACS
  creator object
*/
void createAssembler( MPI_Comm comm, int order, int nx, int ny,
                      TACSAssembler **_assembler, TACSCreator **_creator ){
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set the number of nodes/elements on this proc
  int varsPerNode = 2;

  // Set up the creator object
  TACSCreator *creator = new TACSCreator(comm, varsPerNode);

  if (rank == 0){
    // Set the number of elements
    int numNodes = ((order - 1)*nx + 1)*((order - 1)*ny + 1);
    int numElements = nx*ny;

    // Allocate the input arrays into the creator object
    int *ids = new int[ numElements ];
    int *ptr = new int[ numElements+1 ];
    int *conn = new int[ order*order*numElements ];

    // Set the ids
    memset(ids, 0, numElements*sizeof(int));
    for ( int k = 0; k < numElements; k++ ){
      if (k % nx < nx/2 && (k/nx) < nx/2){
        ids[k] = 1;
      }
      else if (k % nx >= nx/2 && (k/nx) >= nx/2){
        ids[k] = 1;
      }
    }

    ptr[0] = 0;
    for ( int k = 0; k < numElements; k++ ){
      // Back out the i, j coordinates from the corresponding
      // element number
      int i = k % nx;
      int j = k/nx;

      // Set the node connectivity
      for ( int jj = 0; jj < order; jj++ ){
        for ( int ii = 0; ii < order; ii++ ){
          conn[order*order*k + ii + order*jj] =
            ((order-1)*i + ii) +
            ((order-1)*j + jj)*((order-1)*nx + 1);
        }
      }
      ptr[k+1] = order*order*(k+1);
    }

    // Set the connectivity
    creator->setGlobalConnectivity(numNodes, numElements,
                                   ptr, conn, ids);
    delete [] conn;
    delete [] ptr;
    delete [] ids;

    // We're over-counting one of the nodes on each edge
    int numBcs = 4*((order - 1)*nx + 1);
    int *bcNodes = new int[ numBcs ];

    for ( int i = 0; i < ((order-1)*nx + 1); i++ ){
      bcNodes[4*i] = i;
      bcNodes[4*i+1] = i + ((order-1)*nx + 1)*ny;
      bcNodes[4*i+2] = i*((order-1)*nx + 1);
      bcNodes[4*i+3] = (i + 1)*((order-1)*nx + 1) - 1;
    }

    // Set the boundary conditions
    creator->setBoundaryConditions(numBcs, bcNodes);
    delete [] bcNodes;

    // Set the node locations
    TacsScalar *Xpts = new TacsScalar[ 3*numNodes ];
    for ( int j = 0; j < (order - 1)*ny + 1; j++ ){
      for ( int i = 0; i < (order - 1)*nx + 1; i++ ){
        int node = i + ((order - 1)*nx + 1)*j;
        Xpts[3*node] = (1.0*i)/nx;
        Xpts[3*node+1] = (1.0*j)/ny;
        Xpts[3*node+2] = 0.0;
      }
    }

    // Set the nodal locations
    creator->setNodes(Xpts);
    delete [] Xpts;
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
  TACSPlaneStressConstitutive *stiff1 =
    new TACSPlaneStressConstitutive(props1);
  TACSPlaneStressConstitutive *stiff2 =
    new TACSPlaneStressConstitutive(props2);

  // Create the model class
  TACSLinearElasticity2D *model1 =
    new TACSLinearElasticity2D(stiff1, TACS_LINEAR_STRAIN);
  TACSLinearElasticity2D *model2 =
    new TACSLinearElasticity2D(stiff2, TACS_LINEAR_STRAIN);

  // Create the element class
  TACSElementBasis *basis = NULL;
  if (order == 2){
    basis = new TACSLinearQuadBasis();
  }
  else if (order == 3){
    basis = new TACSQuadraticQuadBasis();
  }
  else if (order == 4){
    basis = new TACSCubicQuadBasis();
  }
  else if (order == 5){
    basis = new TACSQuarticQuadBasis();
  }
  else if (order == 6){
    basis = new TACSQuinticQuadBasis();
  }
  TACSElement2D *linear_element1 = new TACSElement2D(model1, basis);
  TACSElement2D *linear_element2 = new TACSElement2D(model2, basis);

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

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  TACSAssembler *assembler;
  TACSCreator *creator;

  // Set the dimension of the mesh
  int nx = 64;
  int ny = 64;
  int order = 6;

  createAssembler(comm, order, nx, ny, &assembler, &creator);
  assembler->incref();
  creator->incref();

  TACSParallelMat *mat = assembler->createMat();
  mat->incref();

  TACSBVec *x = assembler->createVec();
  TACSBVec *y_mat = assembler->createVec();
  TACSBVec *y_free = assembler->createVec();
  x->setRand(-1.0, 1.0);
  assembler->applyBCs(x);

  assembler->setVariables(x);
  assembler->testElement(0, 2);

  double tassemble = MPI_Wtime();
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  assembler->assembleJacobian(alpha, beta, gamma, NULL, mat);
  tassemble = MPI_Wtime() - tassemble;
  printf("Assembly time: %e\n", tassemble);

  double tprod = MPI_Wtime();
  for ( int i = 0; i < 20; i++ ){
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
  for ( int i = 0; i < 20; i++ ){
    mat_free->mult(x, y_free);
  }
  tprod_free = MPI_Wtime() - tprod_free;
  printf("Matrix-free matrix-vector product time: %e\n", tprod_free);

  assembler->applyBCs(y_mat);
  y_mat->axpy(-1.0, y_free);
  TacsScalar norm = y_mat->norm();
  if (rank == 0){
    printf("Residual norm of the difference: %e\n", norm);
  }

  mat->decref();

  assembler->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}