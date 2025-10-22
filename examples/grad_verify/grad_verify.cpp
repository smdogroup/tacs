#include "TACSAssembler.h"
#include "TACSCompliance.h"
#include "TACSInducedFailure.h"
#include "TACSIsoShellConstitutive.h"
#include "TACSKSFailure.h"
#include "TACSMeshLoader.h"
#include "TACSShellElementDefs.h"
#include "TACSStructuralMass.h"

/*
  The following code takes as input a bdf file (which contains only
  shell elements) and computes the gradient of a number of functions
  of interest using either finite-difference, complex-step (when TACS
  is compiled in complex mode) or the adjoint method. This code can be
  used for gradient verification purposes.

  Useage:

  ./shell_grad_verify input_file.bdf [options]

  where [options] may consist of any of the following:

  test             perform additional tests on the functions/elements
  fd               indicate whether to use finite-differece or not
  fd=%f            specify the finite-difference (or CS) interval
  num_threads=%d   specify the number of pthreads
*/
int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Set the communicator to use for this example
  MPI_Comm comm = MPI_COMM_WORLD;

  // The number of pthreads to use
  int num_threads = 1;

  for (int k = 1; k < argc; k++) {
    if (sscanf(argv[k], "num_threads=%d", &num_threads) == 1) {
      break;
    }
  }

  // Check if the input file exists, if not quit
  const char *bdf_file = argv[1];
  FILE *file_test = NULL;
  if ((file_test = fopen(bdf_file, "r"))) {
    fclose(file_test);
  } else {
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
      fprintf(stderr, "Input file %s not found\n", bdf_file);
    }

    MPI_Finalize();
    return 1;
  }

  // Create the TACS mesh loader class and load in the data file
  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();

  // Scan the BDF file to get the input data
  mesh->scanBDFFile(bdf_file);

  // Retrieve the MPI rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set properties needed to create stiffness object
  TacsScalar rho = 2500.0;  // density, kg/m^3
  TacsScalar specific_heat = 921.096;
  TacsScalar E = 70e9;       // elastic modulus, Pa
  TacsScalar nu = 0.3;       // poisson's ratio
  TacsScalar ys = 350e6;     // yield stress, Pa
  TacsScalar cte = 24.0e-6;  // Coefficient of thermal expansion
  TacsScalar kappa = 230.0;  // Thermal conductivity

  TACSMaterialProperties *props =
      new TACSMaterialProperties(rho, specific_heat, E, nu, ys, cte, kappa);

  TACSShellTransform *transform = new TACSShellNaturalTransform();

  // Set up the elements for each component using a new design
  // variable number for each stiffness object. This code assumes
  // that the mesh is entirely composed of shell elements. This
  // code applies fixed material properties that may not be
  // appropriate for some cases (for instance if the mesh dimensions
  // are not in meters). However, this code can still be be
  // use for gradient verification purposes.
  int num_comp = mesh->getNumComponents();
  for (int k = 0; k < num_comp; k++) {
    const char *descriptor = mesh->getElementDescript(k);
    TacsScalar thickness = 0.01;
    int thickness_index = k;
    TacsScalar min_thickness = 0.01;
    TacsScalar max_thickness = 0.20;
    TACSShellConstitutive *con = new TACSIsoShellConstitutive(
        props, thickness, thickness_index, min_thickness, max_thickness);

    // Initialize element object
    TACSElement *elem = TacsCreateShellByName(descriptor, transform, con);

    // Set the element object into the mesh
    mesh->setElement(k, elem);
  }

  // The number of variables per node in the mesh
  int vars_per_node = 6;

  // Create the TACSAssembler object
  TACSAssembler *assembler = mesh->createTACS(vars_per_node);
  assembler->incref();

  // Set the number of threads
  assembler->setNumThreads(num_threads);

  // Create an TACSToFH5 object for writing output to files
  int write_flag = (TACS_OUTPUT_CONNECTIVITY | TACS_OUTPUT_NODES |
                    TACS_OUTPUT_DISPLACEMENTS | TACS_OUTPUT_STRAINS |
                    TACS_OUTPUT_STRESSES | TACS_OUTPUT_EXTRAS);
  TACSToFH5 *f5 =
      new TACSToFH5(assembler, TACS_BEAM_OR_SHELL_ELEMENT, write_flag);
  f5->incref();

  // Set a gravity load
  TACSAuxElements *aux = new TACSAuxElements(assembler->getNumElements());

  // // Add a traction in the negative z-direction
  // TacsScalar tx = 0.0, ty = 0.0, tz = -1e3;
  // for ( int k = 0; k < mesh->getNumComponents(); k++ ){
  //   // Get the element description
  //   const char *descript = mesh->getElementDescript(k);

  //   // Allocate the associate gravity traction
  //   TACSElement *elem = NULL;
  //   if (strcmp(elem_descript, "CQUAD4") == 0){
  //     elem = new TACSShellTraction<2>(tx, ty, tz);
  //   }
  //   else if (strcmp(elem_descript, "CQUAD9") == 0 ||
  //            strcmp(elem_descript, "CQUAD") == 0){
  //     elem = new TACSShellTraction<3>(tx, ty, tz);
  //   }
  //   else if (strcmp(elem_descript, "CQUAD16") == 0){
  //     elem = new TACSShellTraction<4>(tx, ty, tz);
  //   }

  //   // Add the element traction to the mesh
  //   mesh->addAuxElement(aux, k, elem);
  // }

  // Add the tractions/loads to the TACSAssembler object
  assembler->setAuxElements(aux);

  //  Set up the ks functions
  static const int NUM_FUNCS = 12;
  TACSFunction *funcs[NUM_FUNCS];

  // Place functions into the func list
  funcs[0] = new TACSStructuralMass(assembler);
  funcs[1] = new TACSCompliance(assembler);

  // Set the discrete and continuous KS functions
  TACSKSFailure *func = new TACSKSFailure(assembler, 20.0);
  func->setKSFailureType(TACSKSFailure::DISCRETE);
  funcs[2] = func;

  func = new TACSKSFailure(assembler, 20.0);
  func->setKSFailureType(TACSKSFailure::CONTINUOUS);
  funcs[3] = func;

  // Set the induced norm failure types
  TACSInducedFailure *ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
  funcs[4] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL);
  funcs[5] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL_SQUARED);
  funcs[6] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL_SQUARED);
  funcs[7] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::POWER);
  funcs[8] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER);
  funcs[9] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::POWER_SQUARED);
  funcs[10] = ifunc;

  ifunc = new TACSInducedFailure(assembler, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER_SQUARED);
  funcs[11] = ifunc;

  for (int k = 0; k < NUM_FUNCS; k++) {
    funcs[k]->incref();
  }

  // Decrease the reference count to the mesh loader class
  mesh->decref();

  // Create the structural matrix/preconditioner/KSM to use in conjunction
  // with the aerostructural solution
  TACSSchurMat *mat = assembler->createSchurMat();
  TACSBVec *ans = assembler->createVec();
  TACSBVec *res = assembler->createVec();
  ans->incref();
  res->incref();

  // Create the preconditioner
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  TACSSchurPc *pc = new TACSSchurPc(mat, lev, fill, reorder_schur);
  pc->setMonitorFactorFlag(0);
  pc->setAlltoallAssemblyFlag(1);

  // Create GMRES object
  int gmres_iters = 4, nrestart = 0, isflexible = 0;
  GMRES *gmres = new GMRES(mat, pc, gmres_iters, nrestart, isflexible);
  gmres->incref();

  // Set the GMRES tolerances
  double rtol = 1e-12, atol = 1e-30;
  gmres->setTolerances(rtol, atol);

  // Alllocate the design variable vector/arrays
  TACSBVec *xvals = assembler->createDesignVec();
  xvals->incref();

  TACSBVec *xtemp = assembler->createDesignVec();
  xtemp->incref();

  // Set a random perturbation vector
  TACSBVec *pert = assembler->createDesignVec();
  pert->incref();
  pert->setRand(1.0, 2.0);

  TACSBVec *dfdx[NUM_FUNCS];
  for (int k = 0; k < NUM_FUNCS; k++) {
    dfdx[k] = assembler->createDesignVec();
    dfdx[k]->incref();
  }

  // Set the design variables
  assembler->getDesignVars(xvals);

  // Scan any remaining arguments that may be required
  int test_flag = 0;
  double dh = 1e-6;
  for (int k = 0; k < argc; k++) {
    if (sscanf(argv[k], "dh=%lf", &dh) == 1) {
    }
    if (strcmp(argv[k], "test") == 0) {
      test_flag = 1;
    }
  }

  // Print out the finite-difference interval
  if (rank == 0) {
#ifdef TACS_USE_COMPLEX
    printf("Complex-step interval: %le\n", dh);
#else
    printf("Finite-difference interval: %le\n", dh);
#endif
  }

  // First solve the governing equations
  assembler->zeroVariables();
  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();
  gmres->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Write the solution to an .f5 file
  if (f5) {
    // Find the location of the first '.' in the string
    char *outfile = new char[strlen(bdf_file) + 5];
    int len = strlen(bdf_file);
    int i = len - 1;
    for (; i >= 0; i--) {
      if (bdf_file[i] == '.') {
        break;
      }
    }

    // Copy the entire input string to the filename
    strcpy(outfile, bdf_file);

    // Over-write the last part of the string
    strcpy(&outfile[i], ".f5");

    // Write out the file
    f5->writeToFile(outfile);
    delete[] outfile;
  }

  // Evaluate each of the functions
  TacsScalar fvals[NUM_FUNCS];
  assembler->evalFunctions(NUM_FUNCS, funcs, fvals);

  // Evaluate the partial derivative of each function w.r.t. x
  assembler->addDVSens(1.0, NUM_FUNCS, funcs, dfdx);

  // Create the adjoint variables for each function of interest
  TACSBVec *adjoints[NUM_FUNCS];
  TACSBVec *dfdu[NUM_FUNCS];
  for (int j = 0; j < NUM_FUNCS; j++) {
    adjoints[j] = assembler->createVec();
    adjoints[j]->incref();
    dfdu[j] = assembler->createVec();
    dfdu[j]->incref();
  }

  // Evaluate the partial derivative of each function of interest
  // w.r.t. the state variables and compute the adjoint
  assembler->addSVSens(1.0, 0.0, 0.0, NUM_FUNCS, funcs, dfdu);

  // Solve all of the adjoint equations
  for (int j = 0; j < NUM_FUNCS; j++) {
    gmres->solve(dfdu[j], adjoints[j]);
  }

  // Evaluate the adjoint-residual product
  assembler->addAdjointResProducts(-1.0, NUM_FUNCS, adjoints, dfdx);

  // Delete all of the adjoints
  for (int j = 0; j < NUM_FUNCS; j++) {
    adjoints[j]->decref();
    dfdu[j]->decref();
  }

  // Add up the contributions across all processors
  for (int k = 0; k < NUM_FUNCS; k++) {
    dfdx[k]->beginSetValues(TACS_ADD_VALUES);
  }
  for (int k = 0; k < NUM_FUNCS; k++) {
    dfdx[k]->endSetValues(TACS_ADD_VALUES);
  }

  // Test the function that will be used
  if (test_flag) {
    // Test the derivatives of the functions of interest w.r.t. the
    // state and design variables
    for (int k = 0; k < NUM_FUNCS; k++) {
      assembler->testFunction(funcs[k], dh);
    }
  }

  TacsScalar fd[NUM_FUNCS], dfdp[NUM_FUNCS];
  for (int k = 0; k < NUM_FUNCS; k++) {
    dfdp[k] = dfdx[k]->dot(pert);
  }

#if TACS_USE_COMPLEX
  // (xvals[k] + dh)
  xtemp->copyValues(xvals);
  xtemp->axpy(TacsScalar(0.0, dh), pert);

  // Set the design variables
  assembler->setDesignVars(xtemp);

  // Solve the problem
  assembler->zeroVariables();
  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();
  gmres->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Evaluate all of the functions of interest
  assembler->evalFunctions(NUM_FUNCS, funcs, fd);

  // Compute their derivative based on the complex component
  // of the function value
  for (int j = 0; j < NUM_FUNCS; j++) {
    fd[j] = TacsImagPart(fvals[j]) / dh;
  }

#else   // !TACS_USE_COMPLEX
  // Test the structural sensitivities
  // (xvals[ok] + dh)
  xtemp->copyValues(xvals);
  xtemp->axpy(dh, pert);

  // Set the perturbed value of the design variables
  assembler->setDesignVars(xtemp);

  // Solve the problem
  assembler->zeroVariables();
  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();
  gmres->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Evaluate the function at the perturbed solution
  assembler->evalFunctions(NUM_FUNCS, funcs, fvals);

  // (xvals[k] - dh)
  xtemp->copyValues(xvals);
  xtemp->axpy(-dh, pert);

  // Set the perturbed values of the design variables
  assembler->setDesignVars(xtemp);

  // Solve the finite-element problem at the perturbed values
  // of the design variables
  assembler->zeroVariables();
  assembler->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();
  gmres->solve(res, ans);
  ans->scale(-1.0);
  assembler->setVariables(ans);

  // Complete the finite-difference computation for each function
  assembler->evalFunctions(NUM_FUNCS, funcs, fd);
  for (int k = 0; k < NUM_FUNCS; k++) {
    fd[k] = 0.5 * (fvals[k] - fd[k]) / dh;
  }
#endif  // TACS_USE_COMPLEX

  if (rank == 0) {
    printf("Structural sensitivities\n");
    for (int k = 0; k < NUM_FUNCS; k++) {
      printf("Sensitivities for funcion %s\n", funcs[k]->getObjectName());
      printf("%25s %25s %25s\n", "Adjoint", "FD/CS", "Abs. error");
      printf("%25.15e %25.15e %25.15e\n", TacsRealPart(dfdp[k]),
             TacsRealPart(fd[k]), TacsRealPart((dfdp[k] - fd[k]) / fd[k]));
    }
  }

  ans->decref();
  res->decref();
  assembler->decref();
  gmres->decref();

  for (int k = 0; k < NUM_FUNCS; k++) {
    funcs[k]->incref();
  }

  MPI_Finalize();
  return (0);
}
