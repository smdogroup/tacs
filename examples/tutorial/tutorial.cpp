/*
  A pressure-loaded plate example for TACS.

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved.
*/

#include "TACSAssembler.h"
#include "TACSElement2D.h"
#include "TACSHeatFlux.h"
#include "TACSInducedFailure.h"
#include "TACSKSFailure.h"
#include "TACSLinearElasticity.h"
#include "TACSQuadBasis.h"
#include "TACSStructuralMass.h"
#include "TACSThermoelasticity.h"
#include "TACSToFH5.h"

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
int main(int argc, char *argv[]) {
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
  for (int k = 0; k < argc; k++) {
    int _nx, _ny;
    if (sscanf(argv[k], "nx=%d", &_nx) == 1) {
      nx = (_nx > 2 ? _nx : 2);
    }
    if (sscanf(argv[k], "ny=%d", &_ny) == 1) {
      ny = (_ny > 2 ? _ny : 2);
    }
  }

  // Look for a finite-difference step interval
  double dh = 1e-6;
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "dh=%lf", &dh) == 1) {
      // Check that the step size makes sense
      if (dh < 0) {
        dh = 1e-6;
      }
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
  // going to be equal to 2 (You can find this value by checking
  // with element->getVarsPerNode() which returns the number
  // of unknowns per node)
  int varsPerNode = 3;

  int nodesPerProc = ((nx + 1) * (ny + 1)) / size;
  int elemsPerProc = (nx * ny) / size;

  int numOwnedNodes = nodesPerProc;
  int numElements = elemsPerProc;

  // On the last rank, adjust the ownership so we get the
  // total that we need
  if (rank == size - 1) {
    numOwnedNodes = (nx + 1) * (ny + 1) - nodesPerProc * (size - 1);
    numElements = nx * ny - elemsPerProc * (size - 1);
  }

  // There are no dependent nodes in this problem
  int numDependentNodes = 0;
  TACSAssembler *assembler = new TACSAssembler(
      tacs_comm, varsPerNode, numOwnedNodes, numElements, numDependentNodes);
  assembler->incref();  // Increase the reference count to TACSAssembler

  // Set the global element index for the first and last element
  // in the partition
  int firstElem = rank * elemsPerProc;
  int firstNode = rank * nodesPerProc;

  int lastElem = (rank + 1) * elemsPerProc;
  int lastNode = (rank + 1) * nodesPerProc;
  if (rank == size - 1) {
    lastElem = nx * ny;
    lastNode = (nx + 1) * (ny + 1);
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
  int *ptr = new int[numElements + 1];
  int *conn = new int[4 * numElements];

  ptr[0] = 0;
  for (int k = 0, elem = firstElem; elem < lastElem; k++, elem++) {
    // Back out the i, j coordinates from the corresponding
    // element number
    int i = elem % nx;
    int j = elem / nx;

    // Set the node connectivity
    conn[4 * k] = i + j * (nx + 1);
    conn[4 * k + 1] = i + 1 + j * (nx + 1);
    conn[4 * k + 2] = i + (j + 1) * (nx + 1);
    conn[4 * k + 3] = i + 1 + (j + 1) * (nx + 1);
    ptr[k + 1] = 4 * (k + 1);
  }

  // Set the connectivity
  assembler->setElementConnectivity(ptr, conn);
  delete[] conn;
  delete[] ptr;

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

  // Create the element. This element class consists of a constitutive
  // object, which stores information about the material, a model class
  // which computes the variational form of the governing equations based
  // on input from the basis class, which contains the basis functions and
  // quadrature scheme for the element.
  TACSElementBasis *linear_basis = new TACSLinearQuadBasis();

  // Create the plane stress constitutive object
  TACSPlaneStressConstitutive *stiff = new TACSPlaneStressConstitutive(props);
  stiff->incref();

  // Create and set the elements
  TACSElement **elements = new TACSElement *[numElements];

  // Create the auxiliary element class - we'll use this to apply
  // surface tractions
  // TACSAuxElements *aux = new TACSAuxElements(numElements);

  for (int k = 0, elem = firstElem; elem < lastElem; k++, elem++) {
    // Set the thickness design variable = the element number
    TacsScalar t = 1.0;
    int tNum = elem;

    // Create the plane stress constitutive object
    TACSPlaneStressConstitutive *stiff =
        new TACSPlaneStressConstitutive(props, t, tNum);

    // Create the element class
    TACSLinearThermoelasticity2D *model =
        new TACSLinearThermoelasticity2D(stiff, TACS_LINEAR_STRAIN);
    elements[k] = new TACSElement2D(model, linear_basis);

    // Create a surface traction associated with this element and add
    // it to the auxilary elements. Note that the element number must
    // correspond to the local element number used for this processor.
    // TacsScalar tx = 0.0, ty = 0.0, tz = -1e5;
    // TACSShellTraction<2> *trac = new TACSShellTraction<2>(tx, ty, tz);
    // aux->addElement(k, trac);
  }

  // Set the elements into the mesh
  assembler->setElements(elements);
  delete[] elements;

  // Set the boundary conditions - this will only record the
  // boundary conditions on its own nodes
  for (int i = 0; i < nx + 1; i++) {
    int nodes[] = {i, i + (nx + 1) * ny, i * (nx + 1), (i + 1) * (nx + 1) - 1};
    TacsScalar values[] = {0.0, 0.0, 1.0 * i};
    int vars[] = {0, 1, 2};
    assembler->addBCs(4, nodes, 3, vars, values);
  }

  // Reorder the nodal variables
  int use_schur_mat = 1;
  TACSAssembler::OrderingType order_type = TACSAssembler::ND_ORDER;
  TACSAssembler::MatrixOrderingType mat_type = TACSAssembler::APPROXIMATE_SCHUR;

  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "AMD") == 0) {
      order_type = TACSAssembler::AMD_ORDER;
    } else if (strcmp(argv[k], "RCM") == 0) {
      order_type = TACSAssembler::RCM_ORDER;
    } else if (strcmp(argv[k], "ND") == 0) {
      order_type = TACSAssembler::ND_ORDER;
    } else if (strcmp(argv[k], "TACS_AMD") == 0) {
      order_type = TACSAssembler::TACS_AMD_ORDER;
    } else if (strcmp(argv[k], "DirectSchur") == 0) {
      mat_type = TACSAssembler::DIRECT_SCHUR;
    } else if (strcmp(argv[k], "ApproximateSchur") == 0) {
      mat_type = TACSAssembler::APPROXIMATE_SCHUR;
      use_schur_mat = 0;
    } else if (strcmp(argv[k], "AdditiveSchwarz") == 0) {
      mat_type = TACSAssembler::ADDITIVE_SCHWARZ;
      use_schur_mat = 0;
    }
  }

  // If we're going to use the FEMat class, then we don't need to
  // perform a reordering
  if (!use_schur_mat) {
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
    int i = node % (nx + 1);
    int j = node / (nx + 1);
    Xpts[k] = i * Lx / nx;
    Xpts[k + 1] = j * Ly / ny;
  }

  // Reorder the vector if required
  assembler->reorderVec(X);

  // Set the node locations
  assembler->setNodes(X);

  // Set the auxiliary elements
  // assembler->setAuxElements(aux);

  // Solve the problem and set the variables into TACS
  TACSMat *kmat = NULL;
  TACSMat *mmat = NULL;
  TACSPc *pc = NULL;

  // Depending on the input options, solve the
  int lev_fill = 5;  // ILU(k) fill in
  int fill = 8.0;    // Expected number of non-zero entries

  // Options for the ApproximateSchur preconditioner class
  int inner_gmres_iters = 10;
  double inner_rtol = 1e-4, inner_atol = 1e-30;

  // These calls compute the symbolic factorization and allocate
  // the space required for the preconditioners
  if (use_schur_mat) {
    // Set the level of fill to be large
    lev_fill = 1000;

    // Create the FE matrix
    TACSSchurMat *_kmat = assembler->createSchurMat(order_type);
    TACSSchurMat *_mmat = assembler->createSchurMat();
    int reorder_schur = 1;
    pc = new TACSSchurPc(_kmat, lev_fill, fill, reorder_schur);
    kmat = _kmat;
    mmat = _mmat;
  } else {
    // Adjust the level of fill based on the input argument
    for (int k = 0; k < argc; k++) {
      int _lev_fill;
      if (sscanf(argv[k], "lev_fill=%d", &_lev_fill) == 1) {
        lev_fill = _lev_fill;
      }
    }

    // Create the distributed matrix class
    TACSParallelMat *_kmat = assembler->createMat();
    TACSParallelMat *_mmat = assembler->createMat();
    pc = new TACSApproximateSchur(_kmat, lev_fill, fill, inner_gmres_iters,
                                  inner_rtol, inner_atol);
    kmat = _kmat;
    mmat = _mmat;
  }
  mmat->incref();
  kmat->incref();
  pc->incref();

  // Allocate space for the vectors
  TACSBVec *force = assembler->createVec();
  force->incref();
  TACSBVec *res = assembler->createVec();
  res->incref();
  TACSBVec *ans = assembler->createVec();
  ans->incref();
  TACSBVec *tmp = assembler->createVec();
  tmp->incref();

  // Set all components of the vector to 1.0 and apply boundary
  // conditions
  force->set(1.0);
  assembler->setBCs(force);

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
  assembler->assembleJacobian(alpha, beta, gamma, res, kmat);
  res->axpy(-1.0, force);  // res = Ku - f

  // This call copies then factors the matrix
  double t0 = MPI_Wtime();
  pc->factor();
  t0 = MPI_Wtime() - t0;
  printf("[%d] Factor time %f \n", rank, t0);

  // Now, set up the solver
  int use_gmres = 1;
  int gmres_iters = 80;
  int nrestart = 2;     // Number of allowed restarts
  int is_flexible = 1;  // Is a flexible preconditioner?

  // Arguments for the ApproximateSchur preconditioner
  int outer_iters = 15;      // Outer subspace size
  int max_outer_iters = 45;  // Maximum number of outer iterations

  // Create the Krylov Subspace Method (KSM) object
  TACSKsm *ksm = NULL;
  int freq = 1;
  if (use_gmres) {
    ksm = new GMRES(kmat, pc, gmres_iters, nrestart, is_flexible);
    ksm->setMonitor(new KSMPrintStdout("GMRES", rank, freq));
  } else {
    ksm = new GCROT(kmat, pc, outer_iters, max_outer_iters, gmres_iters,
                    is_flexible);
    ksm->setMonitor(new KSMPrintStdout("GCROT", rank, freq));
  }
  ksm->incref();

  // Test the actual residual
  ksm->solve(res, ans);
  kmat->mult(ans, tmp);
  tmp->axpy(-1.0, res);
  TacsScalar norm = tmp->norm();
  if (rank == 0) {
    printf("|Ax - b|: %15.5e\n", TacsRealPart(norm));
  }

  // Assemble the residual and print the result
  ans->scale(-1.0);
  assembler->setVariables(ans);
  assembler->assembleRes(res);
  res->axpy(-1.0, force);  // res = Ku - f
  norm = res->norm();
  if (rank == 0) {
    printf("|R|:      %15.5e\n", TacsRealPart(norm));
  }

  // Output for visualization
  ElementType etype = TACS_PLANE_STRESS_ELEMENT;
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(assembler, etype, write_flag);
  f5->incref();
  f5->writeToFile("tutorial.f5");
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

  // The function that we will use: The KS failure function evaluated
  // over all the elements in the mesh
  double ksRho = 10.0;
  TACSKSFailure *ksfunc = new TACSKSFailure(assembler, ksRho);
  ksfunc->setKSFailureType(TACSKSFailure::CONTINUOUS);
  TACSFunction *func = ksfunc;
  func->incref();

  // Allocate an array for the design variable values
  TACSBVec *x = assembler->createDesignVec();
  x->incref();
  assembler->getDesignVars(x);
  assembler->getNodes(X);

  // Evaluate the function
  TacsScalar ksFuncVal = 0.0;
  assembler->evalFunctions(1, &func, &ksFuncVal);

  // Now, compute the total derivative of the function of interest
  TACSBVec *dfdx = assembler->createDesignVec();
  TACSBVec *dfdX = assembler->createNodeVec();
  dfdx->incref();
  dfdX->incref();

  // Evaluate the partial derivative w.r.t. the design variables
  assembler->addDVSens(1.0, 1, &func, &dfdx);
  assembler->addXptSens(1.0, 1, &func, &dfdX);

  // Evaluate the partial derivative
  TACSBVec *dfdu = assembler->createVec();
  dfdu->incref();

  // Add the partial derivative of the function w.r.t. the state
  // variables. In this case, we do not need to zero dfdu since it is
  // zeroed on initialization, however, in general it is good practice
  // to zero them unless you're absolutely sure...
  dfdu->zeroEntries();
  assembler->addSVSens(alpha, beta, gamma, 1, &func, &dfdu);

  // Solve for the adjoint variables
  assembler->assembleJacobian(alpha, beta, gamma, res, kmat,
                              TACS_MAT_TRANSPOSE);
  pc->factor();
  ksm->solve(dfdu, ans);

  // Compute the total derivative
  assembler->addAdjointResProducts(-1.0, 1, &ans, &dfdx);
  assembler->addAdjointResXptSensProducts(-1.0, 1, &ans, &dfdX);

  dfdx->beginSetValues(TACS_ADD_VALUES);
  dfdX->beginSetValues(TACS_ADD_VALUES);
  dfdx->endSetValues(TACS_ADD_VALUES);
  dfdX->endSetValues(TACS_ADD_VALUES);

  // Now check with a finite-difference projected derivative
  TACSBVec *px = assembler->createDesignVec();
  TACSBVec *pX = assembler->createNodeVec();
  px->incref();
  pX->incref();

  // Set the entries in the px vector based on the value
  // of the dfdx vector
  TacsScalar *px_array, *dfdx_array;
  int num_design_vars = px->getArray(&px_array);
  dfdx->getArray(&dfdx_array);

  for (int i = 0; i < num_design_vars; i++) {
    if (TacsRealPart(dfdx_array[i]) >= 0.0) {
      px_array[i] = 1.0;
    } else {
      px_array[i] = -1.0;
    }
  }

  TacsScalar *pX_array, *dfdX_array;
  int local_node_size = pX->getArray(&pX_array);
  dfdX->getArray(&dfdX_array);

  for (int i = 0; i < local_node_size; i++) {
    if (TacsRealPart(dfdX_array[i]) >= 0.0) {
      pX_array[i] = 1.0;
    } else {
      pX_array[i] = -1.0;
    }
  }

  // Compute the projected derivative along the directin px
  TacsScalar proj_deriv = px->dot(dfdx);
  TacsScalar proj_node_deriv = pX->dot(dfdX);

  TACSBVec *xtemp = assembler->createDesignVec();
  xtemp->incref();

  xtemp->copyValues(x);
#ifdef TACS_USE_COMPLEX
  xtemp->axpy(TacsScalar(0.0, dh), px);
#else
  xtemp->axpy(dh, px);
#endif

  // Set the new design variable values
  assembler->setDesignVars(xtemp);
  xtemp->decref();

  // Evaluate the problem again and compare the results
  assembler->zeroVariables();
  assembler->assembleJacobian(alpha, beta, gamma, res, kmat);
  res->axpy(-1.0, force);  // res = Ku - f

  // Factor the preconditioner
  pc->factor();

  // Solve the problem
  ksm->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Evaluate the function
  TacsScalar ksFuncVal1 = 0.0;
  assembler->evalFunctions(1, &func, &ksFuncVal1);

  assembler->setDesignVars(x);

  // Compare with finite-difference and print the result
  if (rank == 0) {
    TacsScalar fd = 0.0;
#ifdef TACS_USE_COMPLEX
    fd = TacsImagPart(ksFuncVal1) / dh;
#else
    fd = (ksFuncVal1 - ksFuncVal) / dh;
#endif
    printf("The %s function value is %15.8f \n", func->getObjectName(),
           TacsRealPart(ksFuncVal));
    printf("The projected derivative is             %20.8e \n",
           TacsRealPart(proj_deriv));
    printf("The finite-difference approximation is  %20.8e \n",
           TacsRealPart(fd));
    printf("The relative error is                   %20.5e \n",
           fabs(TacsRealPart((fd - proj_deriv) / fd)));
  }

#ifdef TACS_USE_COMPLEX
  X->axpy(TacsScalar(0.0, dh), pX);
#else
  X->axpy(dh, pX);
#endif
  // Set the new design variable values
  assembler->setNodes(X);

  // Evaluate the problem again and compare the results
  assembler->zeroVariables();
  assembler->assembleJacobian(alpha, beta, gamma, res, kmat);
  res->axpy(-1.0, force);  // res = Ku - f

  // Factor the preconditioner
  pc->factor();

  // Solve the problem
  ksm->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Evaluate the function
  ksFuncVal1 = 0.0;
  assembler->evalFunctions(1, &func, &ksFuncVal1);

  // Compare with finite-difference and print the result
  if (rank == 0) {
    TacsScalar fd = 0.0;
#ifdef TACS_USE_COMPLEX
    fd = TacsImagPart(ksFuncVal1) / dh;
#else
    fd = (ksFuncVal1 - ksFuncVal) / dh;
#endif
    printf("The %s function value is %15.8f \n", func->getObjectName(),
           TacsRealPart(ksFuncVal));
    printf("The projected derivative is             %20.8e \n",
           TacsRealPart(proj_node_deriv));
    printf("The finite-difference approximation is  %20.8e \n",
           TacsRealPart(fd));
    printf("The relative error is                   %20.5e \n",
           fabs(TacsRealPart((fd - proj_node_deriv) / fd)));
  }

  // Clean up data required by the adjoint
  func->decref();
  dfdu->decref();
  x->decref();
  X->decref();
  px->decref();
  pX->decref();
  dfdx->decref();
  dfdX->decref();

  // Clean up data required for the analysis (and adjoint too)
  ksm->decref();
  pc->decref();
  kmat->decref();
  mmat->decref();
  ans->decref();
  res->decref();
  tmp->decref();
  force->decref();
  assembler->decref();

  MPI_Finalize();

  return (0);
}
