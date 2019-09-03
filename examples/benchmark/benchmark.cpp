#include <stdio.h>

#include "TACSAssembler.h"
#include "GSEP.h"
#include "isoFSDTStiffness.h"
#include "PlaneStressQuad.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "Compliance.h"
#include "KSFailure.h"

/*
  Create a plane stress stiffness object with default properties
*/
void createPlaneStressElement( int order, int num, TACSElement **elem ){
  TacsScalar rho = 1.0, E = 70000, mu = 0.3;
  PlaneStressStiffness *stiff =
    new PlaneStressStiffness(rho, E, E, 0.5*E/(1.0 + mu), mu);

  if (order == 2){
    *elem = new PlaneStressQuad<2>(stiff);
  }
  else if (order == 3){
    *elem = new PlaneStressQuad<3>(stiff);
  }
  else if (order == 4){
    *elem = new PlaneStressQuad<4>(stiff);
  }
}

/*
  Create a shell stiffness object
*/
void createShellElement( int order, int num, TACSElement **elem ){
  TacsScalar thickness = 2.25;
  TacsScalar rho = 2750.0, E = 70000.0, mu = 0.3, kcorr = 5.0/6.0;
  TacsScalar ys = 1.0;

  FSDTStiffness *stiff = new isoFSDTStiffness(rho, E, mu, kcorr, ys,
                                              thickness, num);
  *elem = NULL;
  if (order == 2){
    *elem = new MITCShell<2>(stiff);
  }
  else if (order == 3){
    *elem = new MITCShell<3>(stiff);
  }
  else if (order == 4){
    *elem = new MITCShell<4>(stiff);
  }
}

/*
  Create the shell traction
*/
void createShellTraction( int order, int num, TACSElement **elem ){
  TacsScalar tx = 0.0, ty = 0.0, tz = 10.0;
  *elem = NULL;
  if (order == 2){
    *elem = new TACSShellTraction<2>(tx, ty, tz);
  }
  else if (order == 3){
    *elem = new TACSShellTraction<3>(tx, ty, tz);
  }
  else if (order == 4){
    *elem = new TACSShellTraction<4>(tx, ty, tz);
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
TACSAssembler *create2DModel( MPI_Comm comm, int varsPerNode,
                              int nx, int ny, int order,
                              int firstNode, int lastNode,
                              int firstElem, int lastElem,
                              int noptions, const char *opts[] ){
  // Set the number of nodes
  int numOwnedNodes = lastNode - firstNode;
  int numElements = lastElem - firstElem;

  // There are no dependent nodes in this problem
  TACSAssembler *tacs = new TACSAssembler(comm, varsPerNode,
                                          numOwnedNodes, numElements);

  // The elements are ordered as (i + j*nx)
  int *ptr = new int[ numElements+1 ];
  int *conn = new int[ order*order*numElements ];

  int *c = conn;
  ptr[0] = 0;
  for ( int k = 0, elem = firstElem; elem < lastElem; k++, elem++ ){
    // Back out the i, j coordinates from the corresponding
    // element number
    int i = elem % nx;
    int j = elem/nx;

    // Set the node connectivity
    for ( int jj = 0; jj < order; jj++ ){
      for ( int ii = 0; ii < order; ii++ ){
        c[0] = (order-1)*i+ii + ((order-1)*nx+1)*((order-1)*j+jj);
        c++;
      }
    }
    ptr[k+1] = c - conn;
  }

  // Set the connectivity
  tacs->setElementConnectivity(ptr, conn);
  delete [] conn;
  delete [] ptr;

  // Create and set the elements
  TACSElement **elements = new TACSElement*[ numElements ];
  TACSAuxElements *aux = new TACSAuxElements(numElements);

  for ( int k = 0, elem = firstElem; elem < lastElem; k++, elem++ ){
    if (varsPerNode == 2){
      // Create a plane stress model
      createPlaneStressElement(order, elem, &elements[k]);
    }
    else { // This must be a shell
      createShellElement(order, elem, &elements[k]);
      TACSElement *trac;
      createShellTraction(order, elem, &trac);
      aux->addElement(k, trac);
    }
  }

  // Set the elements into the mesh
  tacs->setElements(elements);
  delete [] elements;

  // Set the boundary conditions - this will only record the
  // boundary conditions on its own nodes
  for ( int i = 0; i < (order-1)*nx+1; i++ ){
    int nodes[4];
    nodes[0] = i;
    nodes[1] = i + ((order-1)*nx+1)*((order-1)*ny);
    nodes[2] = i*((order-1)*nx+1);
    nodes[3] = (i+1)*((order-1)*nx+1)-1;
    tacs->addBCs(4, nodes);
  }

  // Reorder the nodal variables
  int reorder = 0;
  enum TACSAssembler::OrderingType order_type = TACSAssembler::ND_ORDER;
  enum TACSAssembler::MatrixOrderingType mat_type =
    TACSAssembler::APPROXIMATE_SCHUR;

  for ( int k = 0; k < noptions; k++ ){
    if (strcmp(opts[k], "AMD") == 0){
      order_type = TACSAssembler::AMD_ORDER; reorder = 1;
    }
    else if (strcmp(opts[k], "RCM") == 0){
      order_type = TACSAssembler::RCM_ORDER; reorder = 1;
    }
    else if (strcmp(opts[k], "ND") == 0){
      order_type = TACSAssembler::ND_ORDER; reorder = 1;
    }
    else if (strcmp(opts[k], "TACS_AMD") == 0){
      order_type = TACSAssembler::TACS_AMD_ORDER; reorder = 1;
    }
    else if (strcmp(opts[k], "DirectSchur") == 0){
      mat_type = TACSAssembler::DIRECT_SCHUR; reorder = 1;
    }
    else if (strcmp(opts[k], "ApproximateSchur") == 0){
      mat_type = TACSAssembler::APPROXIMATE_SCHUR; reorder = 1;
    }
    else if (strcmp(opts[k], "AdditiveSchwarz") == 0){
      mat_type = TACSAssembler::ADDITIVE_SCHWARZ; reorder = 1;
    }
  }

  // If we're going to use the FEMat class, then we don't need to
  // perform a reordering
  if (reorder){
    tacs->computeReordering(order_type, mat_type);
  }

  // Perform initialization - cannot add any more elements/vars etc
  tacs->initialize();

  // Create the node vector
  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Get the local node locations
  TacsScalar *Xpts = NULL;
  X->getArray(&Xpts);
  for ( int k = 0, node = firstNode; node < lastNode; k += 3, node++ ){
    int i = node % ((order-1)*nx+1);
    int j = node/((order-1)*nx+1);
    Xpts[k] = 250.0*i/((order-1)*nx);
    Xpts[k+1] = 250.0*j/((order-1)*ny);
  }

  // Reorder the vector if required
  tacs->reorderVec(X);

  // Set the node locations
  tacs->setNodes(X);
  X->decref();

  // Set the shell order
  tacs->setAuxElements(aux);

  return tacs;
}

/*
  Test the eigenvalue solver using some of the given matrices
*/
void testEigenSolver( TACSAssembler *tacs,
                      int noptions, const char *opts[] ){
  int gmres_iters = 50;
  int levFill = 5;
  double fill = 10.0;
  int inner_iters = 10;

  // Check if any of the options are set
  for ( int k = 0; k < noptions; k++ ){
    if (sscanf(opts[k], "levFill=%d", &levFill) == 1){
      if (levFill < 0){ levFill = 0; }
    }
    if (sscanf(opts[k], "gmres_iters=%d", &gmres_iters) == 1){
      if (gmres_iters < 1){ gmres_iters = 1; }
    }
    if (sscanf(opts[k], "inner_iters=%d", &inner_iters) == 1){
      if (inner_iters < 1){ inner_iters = 1; }
    }
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  TACSBVec *ans = tacs->createVec();  ans->incref();
  TACSBVec *rhs = tacs->createVec();  rhs->incref();
  TACSBVec *res = tacs->createVec();  res->incref();

  TACSPMat *mat = tacs->createMat();       mat->incref();
  TACSPMat *mass_mat = tacs->createMat();  mass_mat->incref();
  TACSPMat *k_mat = tacs->createMat();     k_mat->incref();

  double start = MPI_Wtime();
  tacs->assembleMatType(STIFFNESS_MATRIX, k_mat);
  tacs->assembleMatType(MASS_MATRIX, mass_mat);
  double stop = MPI_Wtime();
  printf("Time to assemble mass and stiffness matrices: %10.4f\n", stop-start);

  // The initial problem is to find (x,lambda),
  // K x = lambda M x
  // (K - sigma M)^{-1} B x = 1.0/(lambda - sigma) x

  TacsScalar sigma = 0.01;
  mat->copyValues(k_mat);
  mat->axpy(-sigma, mass_mat);

  // Set up the spectral transformation
  double rtol = -1.0; // Max-out on iterations
  double atol = 1e-10;
  TACSApproximateSchur *pc = new TACSApproximateSchur(mat, levFill, fill,
                                                      inner_iters, rtol, atol);
  pc->incref();
  pc->factor();

  int isflexible = 1;
  int outer_iters = 30;
  int nrestart = outer_iters;
  GCROT *solver = new GCROT(mat, pc, outer_iters,
                            nrestart, gmres_iters, isflexible);
  solver->incref();
  solver->setTolerances(1e-10, 1e-30);

  int max_eigs = 40;
  EPOperator *ep_op = new EPGeneralizedShiftInvert(sigma, solver, mass_mat);
  SEP *sep = new SEP(ep_op, max_eigs, SEP::FULL);
  sep->incref();
  sep->setTolerances(1e-12, SEP::SMALLEST, 8);
  sep->solve();

  printf("||I - Q^{T}*Q||_F :   %25.10e \n",
         TacsRealPart(sep->checkOrthogonality()));

  for ( int k = 0; k < 10; k++ ){
    TacsScalar error;
    TacsScalar lambda = sep->extractEigenvector(k, ans, &error);

    k_mat->mult(ans, rhs);
    mass_mat->mult(ans, res);
    rhs->axpy(- lambda, res);

    TacsScalar ans_norm = ans->norm();
    TacsScalar rhs_norm = rhs->norm();

    if (rank == 0){
      printf("lambda:                 %25.10e \n", TacsRealPart(lambda));
      printf("||x||:                  %25.10e \n", TacsRealPart(ans_norm));
      printf("||A*x - lambda*B*x||:   %25.10e \n", TacsRealPart(rhs_norm));
      printf("error:                  %25.10e \n", TacsRealPart(error));
      printf("\n");
    }
  }

  res->decref();
  rhs->decref();
  ans->decref();

  mat->decref();
  mass_mat->decref();
  k_mat->decref();

  pc->decref();
  solver->decref();
  sep->decref();
}

/*
  Solve the problem with the specified options
*/
void testSolve( TACSAssembler *tacs,
                int noptions, const char *opts[] ){
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
  for ( int k = 0; k < noptions; k++ ){
    if (sscanf(opts[k], "levFill=%d", &levFill) == 1){
      if (levFill < 0){ levFill = 0; }
    }
    if (sscanf(opts[k], "inner_gmres=%d", &inner_gmres) == 1){
      if (inner_gmres < 1){ inner_gmres = 1; }
    }
    if (sscanf(opts[k], "gmres_iters=%d", &gmres_iters) == 1){
      if (gmres_iters < 1){ gmres_iters = 1; }
    }
  }

  // Create the matrix and associated preconditioner
  for ( int k = 0; k < noptions; k++ ){
    if (strcmp(opts[k], "ApproximateSchur") == 0){
      TACSPMat *_mat = tacs->createMat();
      double inner_rtol = -1.0;
      double inner_atol = 1e-10;
      pc = new TACSApproximateSchur(_mat, levFill, fill,
                                    inner_gmres, inner_rtol, inner_atol);
      mat = _mat;
      break;
    }
    else if (strcmp(opts[k], "DirectSchur") == 0){
      FEMat *_mat = tacs->createFEMat();
      int reorder_schur = 1;
      pc = new PcScMat(_mat, levFill, fill, reorder_schur);
      mat = _mat;
      break;
    }
    else if (strcmp(opts[k], "GaussSeidel") == 0){
      int zero_guess = 0; // Zero the initial guess for psor
      TACSPMat *_mat = tacs->createMat();
      pc = new TACSGaussSeidel(_mat, zero_guess,
                               sor_omega, sor_iters, sor_symmetric);
      mat = _mat;
    }
  }

  // Create the additive Schwarz preconditioner
  if (!pc || !mat){
    isflexible = 0;
    TACSPMat *_mat = tacs->createMat();
    pc = new TACSAdditiveSchwarz(_mat, levFill, fill);
    mat = _mat;
  }
  mat->incref();
  pc->incref();

  // Create the vectors to be used in the solution
  TACSBVec *ans = tacs->createVec();
  TACSBVec *rhs = tacs->createVec();
  ans->incref();
  rhs->incref();
  tacs->zeroVariables();

  double start = MPI_Wtime();
  tacs->assembleJacobian(1.0, 0.0, 0.0, rhs, mat);
  double stop = MPI_Wtime();
  if (rank == 0){
    printf("Matrix assembly time: %10.4f\n", stop-start);
  }

  TacsZeroNumFlops();
  start = MPI_Wtime();
  pc->factor();
  stop = MPI_Wtime();
  double flops = TacsGetNumFlops();
  if (rank == 0){
    printf("Factor time: %10.4f\n", stop-start);
  }

  printf("[%d] FLOPS: %15.5e FLOP rate: %15.5e\n", rank, flops, flops/(stop-start));

  // Set up the problem
  GMRES *solver = new GMRES(mat, pc, gmres_iters, 2, isflexible);
  solver->incref();
  solver->setTolerances(1e-10, 1e-30);
  solver->setMonitor(new KSMPrintStdout(" Iteration", rank, 5));

  start = MPI_Wtime();
  solver->solve(rhs, ans);
  stop = MPI_Wtime();
  if (rank == 0){ printf("Solve time: %10.4f\n", stop-start); }

  solver->decref();
  ans->scale(-1.0);
  tacs->setVariables(ans);

  tacs->assembleRes(rhs);
  TacsScalar rhs_norm = rhs->norm();
  if (rank == 0){
    printf("Residual norm: %10.4e\n", TacsRealPart(rhs_norm));
  }

  // Evaluate the compliance
  TACSFunction *comp = new TACSCompliance(tacs);
  comp->incref();
  TacsScalar compVal = 0.0;
  tacs->evalFunctions(&comp, 1, &compVal);
  if (rank == 0){
    printf("Compliance: %25.12f\n", TacsRealPart(compVal));
  }
  comp->decref();

  pc->decref();
  mat->decref();
  ans->decref();
  rhs->decref();
}

/*
  Test the underlying BCSRMat code
*/
void testBCSRMat( TACSAssembler *tacs ){
  int rank;
  MPI_Comm_rank(tacs->getMPIComm(), &rank);

  FEMat *mat = tacs->createFEMat();
  mat->incref();

  // Assemble the Jacobian
  tacs->assembleJacobian(1.0, 0.0, 0.0, NULL, mat);

  if (rank == 0){
    printf("Running BCSRMat timing tests on the root processor\n");

    // Obtain the BCSRMat classes from the FEMat class
    BCSRMat *B, *E, *F, *C;
    mat->getBCSRMat(&B, &E, &F, &C);
    int size = B->getBlockSize()*B->getRowDim();

    TacsScalar *x = new TacsScalar[ size ];
    TacsScalar *y = new TacsScalar[ size ];

    // Set the x values to 1.0
    for ( int k = 0; k < size; k++ ){
      x[k] = 1.0;
    }

    int max_num_threads = 8;

    // Test the time required for matrix-vector products with various
    // numbers of pthreads
    for ( int p = 1; p <= max_num_threads; p++ ){
      tacs->setNumThreads(p);
      if (p == 1){ TacsZeroNumFlops(); }
      double t0 = MPI_Wtime();
      for ( int k = 0; k < 100; k++ ){
        B->mult(x, y);
      }
      t0 = MPI_Wtime() - t0;

      printf("Time for 100 mat-vec products, num_threads %2d: %15.6f\n",
             p, t0);
      if (p == 1){
        double flops = TacsGetNumFlops();
        printf("FLOPS: %15.5e FLOP rate: %15.5e\n", flops, flops/t0);
      }
    }

    // Compute the time required for matrix factorizing with various numbers
    // of pthreads
    double fill = 10.0;
    int levFill = 1000;

    BCSRMat *Bpc = new BCSRMat(MPI_COMM_SELF, B, levFill, fill);
    Bpc->incref();

    for ( int p = 1; p <= max_num_threads; p++ ){
      tacs->setNumThreads(p);
      double t0 = MPI_Wtime();
      if (p == 1){ TacsZeroNumFlops(); }
      // Copy values and factor the matrix
      Bpc->copyValues(B);
      Bpc->factor();

      t0 = MPI_Wtime() - t0;

      printf("Time to factor the matrix, num_threads %2d: %15.6f\n",
             p, t0);
      if (p == 1){
        double flops = TacsGetNumFlops();
        printf("FLOPS: %15.5e FLOP rate: %15.5e\n", flops, flops/t0);
      }
    }

    // Test the time required for matrix-vector products with various
    // numbers of pthreads
    TacsScalar *xt = new TacsScalar[ size ];

    for ( int p = 1; p <= max_num_threads; p++ ){
      tacs->setNumThreads(p);
      double t0 = MPI_Wtime();
      if (p == 1){ TacsZeroNumFlops(); }
      for ( int k = 0; k < 25; k++ ){
        Bpc->applyFactor(y, x);
      }
      t0 = MPI_Wtime() - t0;

      printf("Time for 25 U^{-1}L^{-1}x, num_threads %2d: %15.6f\n",
             p, t0);
      if (p == 1){
        double flops = TacsGetNumFlops();
        printf("FLOPS: %15.5e FLOP rate: %15.5e\n", flops, flops/t0);
      }

      TacsScalar err = 0.0;
      for ( int i = 0; i < size; i++ ){
        err += (x[i] - 1.0)*(x[i] - 1.0);
      }

      for ( int i = 0; i < size; i++ ){
        if (fabs(TacsRealPart(x[i]) - 1.0) > 1e-8){
          printf("x[%d] = %15.8e\n", i, TacsRealPart(x[i]));
        }
      }

      printf("||x - e||_2: %15.4e\n", TacsRealPart(sqrt(err)));
    }

    delete [] x;
    delete [] y;
    delete [] xt;

    Bpc->decref();

    // Compute the number of matrix-vector products

    BCSRMat *A = new BCSRMat(MPI_COMM_SELF,
                             B, NULL, B, NULL, B, levFill, fill);
    A->incref();

    for ( int p = 1; p <= max_num_threads; p++ ){
      tacs->setNumThreads(p);
      double t0 = MPI_Wtime();
      if (p == 1){ TacsZeroNumFlops(); }

      for ( int k = 0; k < 5; k++ ){
        double alpha = 1.0;
        A->matMultAdd(alpha, B, B);
      }
      t0 = MPI_Wtime() - t0;

      printf("Time for 5 mat-mat products, num_threads %2d: %15.6f\n",
             p, t0);
      if (p == 1){
        double flops = TacsGetNumFlops();
        printf("FLOPS: %15.5e FLOP rate: %15.5e\n", flops, flops/t0);
      }
    }

    A->decref();
  }

  mat->decref();
}

/*
  Test the threaded implementation of the adjoint-residual products
*/
void testDVSensThreads( TACSAssembler *tacs, int numDVs ){
  int rank;
  MPI_Comm_rank(tacs->getMPIComm(), &rank);

  // Maximum number of threads to test
  int max_num_threads = 8;

  // Create a vector of functions
  const int numFuncs = 10;
  TacsScalar funcVals[numFuncs];
  TACSFunction *funcs[numFuncs];
  TACSBVec *adjoints[numFuncs];
  for ( int k = 0; k < numFuncs; k++ ){
    double ks_weight = 50.0;
    funcs[k] = new TACSKSFailure(tacs, ks_weight);
    funcs[k]->incref();

    adjoints[k] = tacs->createVec();
    adjoints[k]->incref();
    adjoints[k]->setRand(-1.0, 1.0);
    tacs->applyBCs(adjoints[k]);
  }

  // Time the evaluation of the partial derivatives
  TacsScalar *fdvSens = new TacsScalar[ numFuncs*numDVs ];
  memset(fdvSens, 0, numFuncs*numDVs*sizeof(TacsScalar));

  // Evaluate the functions
  tacs->evalFunctions(funcs, numFuncs, funcVals);

  // Test the time required for matrix-vector products with various
  // numbers of pthreads
  for ( int p = 1; p <= max_num_threads; p++ ){
    tacs->setNumThreads(p);
    double t0 = MPI_Wtime();
    for ( int k = 0; k < 5; k++ ){
      tacs->addDVSens(1.0, funcs, numFuncs, fdvSens, numDVs);
    }
    t0 = MPI_Wtime() - t0;

    if (rank == 0){
      printf("Time for 5 calls to TACSAssembler::addDVSens(), \
num_threads %2d: %15.6f\n", p, t0);
    }
  }

  // Test the time required for matrix-vector products with various
  // numbers of pthreads
  for ( int p = 1; p <= max_num_threads; p++ ){
    tacs->setNumThreads(p);
    double t0 = MPI_Wtime();
    for ( int k = 0; k < 5; k++ ){
      tacs->addAdjointResProducts(1.0, adjoints, numFuncs,
                                  fdvSens, numDVs);
    }
    t0 = MPI_Wtime() - t0;

    if (rank == 0){
      printf("Time for 5 calls to TACSAssembler::addAdjointResProducts(), \
num_threads %2d: %15.6f\n", p, t0);
    }
  }

  for ( int k = 0; k < numFuncs; k++ ){
    funcs[k]->decref();
    adjoints[k]->decref();
  }

  delete [] fdvSens;
}

int main( int argc, char *argv[] ){
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
  int noptions = 7;
  const char *opts[] =
    {"AMD", "DirectSchur", "nx=50", "ny=50",
     "varsPerNode=6", "order=3", "levFill=1000"};

  for ( int k = 0; k < noptions; k++ ){
    if (sscanf(opts[k], "nx=%d", &nx) == 1){}
    if (sscanf(opts[k], "ny=%d", &ny) == 1){}
    if (sscanf(opts[k], "order=%d", &order) == 1){
      if (order < 2){ order = 2; }
      if (order > 4){ order = 4; }
    }
    if (sscanf(opts[k], "varsPerNode=%d", &varsPerNode) == 1){
      if (!(varsPerNode == 2 || varsPerNode == 6)){
        varsPerNode = 2;
      }
    }
  }

  // Determine the total number of nodes/elements
  int nelems = nx*ny;
  int nnodes = (((order-1)*nx+1)*((order-1)*ny+1));

  // Determin the partition of the nodes (this won't be very good...)
  int elemsPerProc = nelems/size;
  int nodesPerProc = nnodes/size;
  int firstNode = rank*nodesPerProc;
  int lastNode = (rank+1)*nodesPerProc;
  int firstElem = rank*elemsPerProc;
  int lastElem = (rank+1)*elemsPerProc;
  if (rank == size-1){
    lastNode = nnodes;
    lastElem = nelems;
  }

  // Create the TACSAssembler object
  TACSAssembler *tacs =
    create2DModel(comm, varsPerNode, nx, ny, order,
                  firstNode, lastNode, firstElem, lastElem,
                  noptions, opts);
  tacs->incref();

  // Test the BCSRMat implementation with threads
  testBCSRMat(tacs);

  int max_num_threads = 8;
  for ( int k = 1; k <= max_num_threads; k++ ){
    if (rank == 0){
      printf("Solving with numPthreads = %d\n", k);
    }
    tacs->setNumThreads(k);

    // Test solve the first level for now
    testSolve(tacs, noptions, opts);
  }

  tacs->decref();
  MPI_Finalize();

  return (0);
}
