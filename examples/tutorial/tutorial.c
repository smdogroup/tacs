/*
  A pressure-loaded plate example for TACS.

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
*/

#include "TACSAssembler.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"
#include "TACSToFH5.h"
#include "KSFailure.h"

/*
  The following example demonstrates the use of TACS on a pressure
  loaded plate. 

  This code uses the TACSAssembler interface directly. Other creation
  code (TACSCreator/TACSMeshLoader) can also be used to generate
  TACSAssembler instances. Once a TACSAssembler instance has been
  created and initialized, it should be able to be used interchangeably
  withe 
  
  Note: This code does not intelligently partition the mesh. You could
  use TACSCreator to perform the partitioning for you to achieve
  better results, but this code is designed to demonstrate the
  TACSAssembler interface itself.

  The command line inputs (nx, ny) provide the number of elements
  along the x and y directions, respectively.
*/
int main( int argc, char * argv[] ){
  MPI_Init(&argc, &argv);
  
  // Find the MPI rank and size
  MPI_Comm tacs_comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(tacs_comm, &rank);
  MPI_Comm_size(tacs_comm, &size);

  // Set the dimensions of the plate
  double Lx = 1.0;
  double Ly = 1.0;

  // Get the global size of the mesh from the input
  int nx = 30, ny = 30;
  for ( int k = 0; k < argc; k++ ){
    int _nx, _ny;
    if (sscanf(argv[k], "nx=%d", &_nx) == 1){
      nx = (_nx > 2 ? _nx : 2);
    }
    if (sscanf(argv[k], "ny=%d", &_ny) == 1){
      ny = (_ny > 2 ? _ny : 2);
    }
  }

  // Look for a finite-difference step interval
  double dh = 1e-6;
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "dh=%lf", &dh) == 1){
      // Check that the step size makes sense
      if (dh < 0){ dh = 1e-6; }
    }
  }

  /*
    To create TACS we need the following information:

    1. The communicator
    2. The number of variables per node (the same across the entire
    mesh)
    3. The number of nodes that are owned by this processor
    4. The number of elements that are owned by this processor
    5. The number of dependent nodes (nodes that depend linearly
    on other nodes)

    In this example, nx and ny are the number of elements in the
    global element mesh. Note that TACS deals exclusively with global
    node numbers to try to make things easier for the user.

    The ownership range of each processor (the range of node numbers
    that belong to each processor) is calculated using
  */

  // We know in advance that the number of unknowns per node is
  // going to be equal to 6 (You can find this value by checking
  // with element->numDisplacements() which returns the number
  // of displacements (or unknowns) per node)
  int varsPerNode = 6;  

  int nodesPerProc = ((nx+1)*(ny+1))/size;
  int elemsPerProc = (nx*ny)/size;

  int numOwnedNodes = nodesPerProc;
  int numElements = elemsPerProc;

  // On the last rank, adjust the ownership so we get the
  // total that we need
  if (rank == size-1){
    numOwnedNodes = (nx+1)*(ny+1) - nodesPerProc*(size-1);
    numElements = nx*ny - elemsPerProc*(size-1);
  }

  // There are no dependent nodes in this problem
  int numDependentNodes = 0;
  TACSAssembler *tacs = new TACSAssembler(tacs_comm, varsPerNode,
                                          numOwnedNodes, numElements,
                                          numDependentNodes);
  tacs->incref(); // Increase the reference count to TACSAssembler

  // Set the global element index for the first and last element 
  // in the partition
  int firstElem = rank*elemsPerProc;
  int firstNode = rank*nodesPerProc;

  int lastElem = (rank+1)*elemsPerProc;
  int lastNode = (rank+1)*nodesPerProc;
  if (rank == size-1){
    lastElem = nx*ny;
    lastNode = (nx+1)*(ny+1);
  }

  /*
    The element connectivity defines the mapping between the element
    and its corresponding nodes. The node numbers are global. Since
    the number of nodes per element may vary, we also provide a
    pointer into the element connectivity array denoting the begining
    location of each element node list. This data is passed in to
    TACSAssembler directly.

    In this case we know that we only ever have 4 nodes per element.
  */

  // The elements are ordered as (i + j*nx)
  int *ptr = new int[ numElements+1 ];
  int *conn = new int[ 4*numElements ];

  ptr[0] = 0;
  for ( int k = 0, elem = firstElem; elem < lastElem; k++, elem++ ){
    // Back out the i, j coordinates from the corresponding
    // element number
    int i = elem % nx;
    int j = elem/nx;

    // Set the node connectivity
    conn[4*k] = i + j*(nx+1);
    conn[4*k+1] = i+1 + j*(nx+1);
    conn[4*k+2] = i + (j+1)*(nx+1);
    conn[4*k+3] = i+1 + (j+1)*(nx+1);
    ptr[k+1] = 4*(k+1);
  }
  
  // Set the connectivity
  tacs->setElementConnectivity(conn, ptr);
  delete [] conn;
  delete [] ptr;

  // Create and set the elements
  TACSElement **elements = new TACSElement*[ numElements ];

  // Create the auxiliary element class - we'll use this to apply
  // surface tractions
  TACSAuxElements *aux = new TACSAuxElements(numElements);
  
  for ( int k = 0, elem = firstElem; elem < lastElem; k++, elem++ ){
    // Create the constitutive objects
    TacsScalar rho = 2500.0; // Not used
    TacsScalar E = 70e9;
    TacsScalar nu = 0.3;
    TacsScalar kcorr = 5.0/6.0; // The shear correction factor
    TacsScalar yield_stress = 464.0e6;
    TacsScalar thickness = 0.005;

    // Set the thickness design variable = the element number
    int tNum = elem;
    
    // Create the stiffness object
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, 
                                                   yield_stress, thickness, 
                                                   tNum);

    // Create the shell element    
    elements[k] = new MITCShell<2>(stiff);

    // Create a surface traction associated with this element and add
    // it to the auxilary elements. Note that the element number must
    // correspond to the local element number used for this processor.
    TacsScalar tx = 0.0, ty = 0.0, tz = -1e5; 
    TACSShellTraction<2> *trac = new TACSShellTraction<2>(tx, ty, tz);
    aux->addElement(k, trac);
  }

  // Set the elements into the mesh
  tacs->setElements(elements);
  delete [] elements;

  // Set the boundary conditions - this will only record the
  // boundary conditions on its own nodes
  for ( int i = 0; i < nx+1; i++ ){
    int nodes[] = {i, i + (nx+1)*ny, i*(nx+1), (i+1)*(nx+1)-1};
    tacs->addBCs(4, nodes);
  }

  // Reorder the nodal variables
  int use_fe_mat = 1;
  TACSAssembler::OrderingType order_type = TACSAssembler::ND_ORDER;
  TACSAssembler::MatrixOrderingType mat_type = 
    TACSAssembler::APPROXIMATE_SCHUR;

  for ( int k = 0; k < argc; k++ ){
    if (strcmp(argv[k], "AMD") == 0){ 
      order_type = TACSAssembler::AMD_ORDER;
    }
    else if (strcmp(argv[k], "RCM") == 0){ 
      order_type = TACSAssembler::RCM_ORDER;
    }
    else if (strcmp(argv[k], "ND") == 0){ 
      order_type = TACSAssembler::ND_ORDER;
    }
    else if (strcmp(argv[k], "TACS_AMD") == 0){
      order_type = TACSAssembler::TACS_AMD_ORDER;
    }
    else if (strcmp(argv[k], "DirectSchur") == 0){ 
      mat_type = TACSAssembler::DIRECT_SCHUR;
    }
    else if (strcmp(argv[k], "ApproximateSchur") == 0){ 
      mat_type = TACSAssembler::APPROXIMATE_SCHUR; 
      use_fe_mat = 0;
    }
    else if (strcmp(argv[k], "AdditiveSchwarz") == 0){ 
      mat_type = TACSAssembler::ADDITIVE_SCHWARZ; 
      use_fe_mat = 0;
    }
  }

  // If we're going to use the FEMat class, then we don't need to
  // perform a reordering
  if (!use_fe_mat){
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
    int i = node % (nx+1);
    int j = node/(nx+1);
    Xpts[k] = i*Lx/nx;
    Xpts[k+1] = j*Ly/ny;
  }

  // Reorder the vector if required
  tacs->reorderVec(X);  

  // Set the node locations
  tacs->setNodes(X);

  // Set the auxiliary elements
  tacs->setAuxElements(aux);

  // Solve the problem and set the variables into TACS
  TACSMat *kmat = NULL;
  TACSMat *mmat = NULL;
  TACSPc *pc = NULL;

  // Depending on the input options, solve the 
  int lev_fill = 5; // ILU(k) fill in
  int fill = 8.0; // Expected number of non-zero entries

  // Options for the ApproximateSchur preconditioner class
  int inner_gmres_iters = 10; 
  double inner_rtol = 1e-4, inner_atol = 1e-30;

  // These calls compute the symbolic factorization and allocate
  // the space required for the preconditioners
  if (use_fe_mat){
    // Set the level of fill to be large
    lev_fill = 1000;

    // Create the FE matrix
    FEMat *_kmat = tacs->createFEMat(order_type);
    FEMat *_mmat = tacs->createFEMat();
    int reorder_schur = 1;
    pc = new PcScMat(_kmat, lev_fill, fill, reorder_schur);
    kmat = _kmat;
    mmat = _mmat;
  }
  else {
    // Adjust the level of fill based on the input argument
    for ( int k = 0; k < argc; k++ ){
      int _lev_fill;
      if (sscanf(argv[k], "lev_fill=%d", &_lev_fill) == 1){
        lev_fill = _lev_fill;
      }
    }

    // Create the distributed matrix class
    TACSDistMat *_kmat = tacs->createMat();
    TACSDistMat *_mmat = tacs->createMat();
    pc = new TACSApproximateSchur(_kmat, lev_fill, fill, 
                                  inner_gmres_iters, inner_rtol, inner_atol);
    kmat = _kmat;
    mmat = _mmat;
  }
  mmat->incref();
  kmat->incref();
  pc->incref();

  // Assemble the stiffness matrix and residual
  TACSBVec *res = tacs->createVec();  res->incref();
  TACSBVec *ans = tacs->createVec();  ans->incref();
  TACSBVec *tmp = tacs->createVec();  tmp->incref();

  /*
    Assemble the Jacobian of governing equations. Note that the alpha,
    beta and gamma coefficients represent the coefficients of the
    structural governing equations for the variables, their first
    derivatives and their second derivatives, respectively such that:

    kmat = J = alpha*dR/du + beta*dR/d(dot{u}) + gamma*dR/d(ddot{u})

    In general, when the governing equations are non-symmetric, the
    assembleJacobian code takes a matrix orientation argument which
    permits the assembly of the transpose of the governing equations
    such that

    kmat = J^{T}

    This capability is required for the adjoint equations.
  */
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, res, kmat);

  // This call copies then factors the matrix
  double t0 = MPI_Wtime();
  pc->factor(); 
  t0 = MPI_Wtime() - t0;
  printf("[%d] Factor time %f \n", rank, t0);

  // Now, set up the solver
  int use_gmres = 1;
  int gmres_iters = 80; 
  int nrestart = 2; // Number of allowed restarts
  int is_flexible = 1; // Is a flexible preconditioner?

  // Arguments for the ApproximateSchur preconditioner
  int outer_iters = 15; // Outer subspace size
  int max_outer_iters = 45; // Maximum number of outer iterations

  // Create the Krylov Subspace Method (KSM) object
  TACSKsm *ksm = NULL;
  int freq = 1;
  if (use_gmres){
    ksm = new GMRES(kmat, pc, gmres_iters, nrestart, is_flexible);
    ksm->setMonitor(new KSMPrintStdout("GMRES", rank, freq));
  }
  else {
    ksm = new GCROT(kmat, pc, outer_iters, max_outer_iters,
                    gmres_iters, is_flexible);
    ksm->setMonitor(new KSMPrintStdout("GCROT", rank, freq));
  }
  ksm->incref();  

  // Test the actual residual
  ksm->solve(res, ans);
  kmat->mult(ans, tmp);
  tmp->axpy(-1.0, res);
  TacsScalar norm = tmp->norm();
  if (rank == 0){
    printf("|Ax - b|: %15.5e\n", TacsRealPart(norm));
  }

  // Assemble the residual and print the result
  ans->scale(-1.0);
  tacs->setVariables(ans);
  tacs->assembleRes(res);
  norm = res->norm();  
  if (rank == 0){
    printf("|R|: %15.5e\n", TacsRealPart(norm));
  }

  // Output for visualization
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();
  f5->writeToFile("tutorial_output.f5");
  f5->decref();

  /*
    Now we will use TACS to compute the total derivative of a function
    of interest with respect to design variables. The design variable
    numbers are set in the constitutive object. 

    The user is in charge of the design variable numbering. It can be
    either local to each processor or global across all processors.
    When it is a global ordering, as it is in this tutorial, the user
    is responsible for suming the contributions to the gradient across
    all processors in the communicator. In more advanced applications,
    the design vector can be distributed, in which case the the design
    variable numbers are local. In this case, the user can reduce the
    design variable values across all processors in the mesh in a
    consistent manner.

    In this tutorial we will not use geometric design variables which
    modify the nodal locations within the mesh. (The geometric
    derivative evaluation routines are denoted "XptSens" and operate
    on TACSBVec objects of the same dimension/parallel layout as the
    vectors created during the createNodeVec() call.)
  */

  // The number of design variable values
  int numDesignVars = nx*ny;

  // The function that we will use: The KS failure function evaluated
  // over all the elements in the mesh
  double ksRho = 100.0;
  TACSFunction *func = new TACSKSFailure(tacs, ksRho);
  func->incref();

  // Allocate an array for the design variable values
  TacsScalar *x = new TacsScalar[ numDesignVars ];
  memset(x, 0, numDesignVars*sizeof(TacsScalar));
  tacs->getDesignVars(x, numDesignVars);

  /*
    Since the design variable vector is a global vector for this case,
    we find the maximum design variable value. Note that here we use
    TACS_MPI_MAX which is defined for real and complex numbers
    (MPI_MAX is not defined for complex numbers).
  */
  MPI_Allreduce(MPI_IN_PLACE, x, numDesignVars, TACS_MPI_TYPE,
                TACS_MPI_MAX, tacs_comm);
  
  // Now, set the design variables values to ensure that they are
  // locally consistent
  tacs->setDesignVars(x, numDesignVars); 

  // Evaluate the function
  TacsScalar ksFuncVal = 0.0;
  tacs->evalFunctions(&func, 1, &ksFuncVal);

  // Now, compute the total derivative of the function of interest
  TacsScalar *dfdx = new TacsScalar[ numDesignVars ];
  memset(dfdx, 0, numDesignVars*sizeof(TacsScalar));

  // Evaluate the partial derivative w.r.t. the design variables
  tacs->addDVSens(1.0, &func, 1, dfdx, numDesignVars);

  // Evaluate the partial derivative
  TACSBVec *dfdu = tacs->createVec();
  dfdu->incref();

  // Add the partial derivative of the function w.r.t. the state
  // variables. In this case, we do not need to zero dfdu since it is
  // zeroed on initialization, however, in general it is good practice
  // to zero them unless you're absolutely sure...
  dfdu->zeroEntries();
  tacs->addSVSens(alpha, beta, gamma, &func, 1, &dfdu);

  // Solve for the adjoint variables
  ksm->solve(dfdu, ans);
  
  // Compute the total derivative
  tacs->addAdjointResProducts(-1.0, &ans, 1, dfdx, numDesignVars);

  // Add up the contributions across all processors
  MPI_Allreduce(MPI_IN_PLACE, dfdx, numDesignVars, TACS_MPI_TYPE, 
                MPI_SUM, tacs_comm);

  // Now check with a finite-difference projected derivative
  TacsScalar proj_deriv = 0.0;
  for ( int k = 0; k < numDesignVars; k++ ){
    proj_deriv += fabs(dfdx[k]);
#ifdef TACS_USE_COMPLEX
    if (TacsRealPart(dfdx[k]) > 0){
      x[k] = x[k] + TacsScalar(0.0, dh);
    }
    else {
      x[k] = x[k] - TacsScalar(0.0, dh);
    }
#else
    if (dfdx[k] > 0){
      x[k] += dh;
    }
    else {
      x[k] -= dh;
    }
#endif
  }

  // Set the new design variable values
  tacs->setDesignVars(x, numDesignVars);

  // Evaluate the problem again and compare the results
  tacs->zeroVariables();
  tacs->assembleJacobian(alpha, beta, gamma, res, kmat);

  // Factor the preconditioner
  pc->factor();

  // Solve the problem
  ksm->solve(res, ans);
  ans->scale(-1.0);
  tacs->setVariables(ans);
  
  // Evaluate the function
  TacsScalar ksFuncVal1 = 0.0;
  tacs->evalFunctions(&func, 1, &ksFuncVal1);

  // Compare with finite-difference and print the result
  if (rank == 0){
    TacsScalar fd = 0.0;
#ifdef TACS_USE_COMPLEX
    fd = TacsImagPart(ksFuncVal1)/dh;
#else
    fd = (ksFuncVal1 - ksFuncVal)/dh;
#endif
    printf("The %s is %15.8f \n", func->functionName(), 
           TacsRealPart(ksFuncVal));
    printf("The projected derivative is             %20.8e \n", 
           TacsRealPart(proj_deriv));
    printf("The finite-difference approximation is  %20.8e \n", 
           TacsRealPart(fd));
    printf("The relative error is                   %20.5e \n", 
	   fabs(TacsRealPart((fd - proj_deriv)/fd)));
  } 

  // Clean up data required by the adjoint
  delete [] x;
  delete [] dfdx;
  func->decref();
  dfdu->decref();

  // Clean up data required for the analysis (and adjoint too)
  ksm->decref();
  pc->decref();
  kmat->decref();
  mmat->decref();
  ans->decref();
  res->decref();
  tmp->decref();
  tacs->decref();

  MPI_Finalize();

  return (0);
}
