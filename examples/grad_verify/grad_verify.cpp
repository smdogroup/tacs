#include "TACSAssembler.h"
#include "TACSMeshLoader.h"
#include "MITCShell.h"
#include "TACSShellTraction.h"
#include "isoFSDTStiffness.h"
#include "Compliance.h"
#include "KSFailure.h"
#include "InducedFailure.h"
#include "StructuralMass.h"

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
int main( int argc, char * argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Set the communicator to use for this example
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // The number of pthreads to use
  int num_threads = 1;

  for ( int k = 1; k < argc; k++ ){
    if (sscanf(argv[k], "num_threads=%d", &num_threads) == 1){
      break;
    }
  }

  // Check if the input file exists, if not quit
  const char * bdf_file = argv[1];
  FILE *file_test = NULL;
  if ((file_test = fopen(bdf_file, "r"))){
    fclose(file_test);
  }
  else {
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0){
      fprintf(stderr, "Input file %s not found\n", bdf_file);
    }

    MPI_Finalize();
    return 1;
  }

  // Create the TACS mesh loader class and load in the data file
  TACSMeshLoader * mesh = new TACSMeshLoader(comm);
  mesh->incref();

  // Scan the BDF file to get the input data
  mesh->scanBDFFile(bdf_file);

  // Retrieve the MPI rank
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Set up the elements for each component using a new design
  // variable number for each stiffness object. This code assumes
  // that the mesh is entirely composed of shell elements. This
  // code applies fixed material properties that may not be
  // appropriate for some cases (for instance if the mesh dimensions
  // are not in meters). However, this code can still be be
  // use for gradient verification purposes.
  int num_comp = mesh->getNumComponents();
  for ( int k = 0; k < mesh->getNumComponents(); k++ ){
    // Set up the constitutive relationship
    TacsScalar rho = 2700.0, E = 35e8, nu = 0.3, kcorr = 0.8333;
    TacsScalar ys = 434.0e6, t = 0.01, t_min = 0.001, t_max = 0.01;
    int t_num = k;

    FSDTStiffness *con =
      new isoFSDTStiffness(rho, E, nu, kcorr,
                           ys, t, t_num, t_min, t_max);

    // Split up the element id to show the domain decomposition.
    int elem_id = rank;

    // Get the description of the element type to determine what
    // order of shell we're dealing with
    const char *elem_descript = mesh->getElementDescript(k); 

    // Allocate the element
    TACSElement *elem = NULL;
    if (strcmp(elem_descript, "CQUAD4") == 0){
      elem = new MITCShell<2>(con, LINEAR, elem_id);
    }
    else if (strcmp(elem_descript, "CQUAD9") == 0 ||
             strcmp(elem_descript, "CQUAD") == 0){
      elem = new MITCShell<3>(con, LINEAR, elem_id);
    }
    else if (strcmp(elem_descript, "CQUAD16") == 0){
      elem = new MITCShell<4>(con, LINEAR, elem_id);
    }
    
    // Set the element object into the mesh
    mesh->setElement(k, elem);
  }
  
  // The number of variables per node in the mesh
  int vars_per_node = 6;

  // Create the TACSAssembler object
  TACSAssembler * tacs = mesh->createTACS(vars_per_node);
  tacs->incref();

  // Set the number of threads
  tacs->setNumThreads(num_threads);

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_SHELL, write_flag);
  f5->incref();

  // Set a gravity load
  TACSAuxElements *aux = new TACSAuxElements(tacs->getNumElements());
  
  // Add a traction in the negative z-direction
  TacsScalar tx = 0.0, ty = 0.0, tz = -1e3;
  for ( int k = 0; k < mesh->getNumComponents(); k++ ){
    // Get the element description
    const char *elem_descript = mesh->getElementDescript(k); 

    // Allocate the associate gravity traction
    TACSElement *elem = NULL;
    if (strcmp(elem_descript, "CQUAD4") == 0){
      elem = new TACSShellTraction<2>(tx, ty, tz);
    }
    else if (strcmp(elem_descript, "CQUAD9") == 0 ||
             strcmp(elem_descript, "CQUAD") == 0){
      elem = new TACSShellTraction<3>(tx, ty, tz);
    }
    else if (strcmp(elem_descript, "CQUAD16") == 0){
      elem = new TACSShellTraction<4>(tx, ty, tz);
    }
    
    // Add the element traction to the mesh
    mesh->addAuxElement(aux, k, elem);
  }

  // Add the tractions/loads to the TACSAssembler object
  tacs->setAuxElements(aux);
  
  //  Set up the ks functions
  static const int NUM_FUNCS = 12;
  TACSFunction *funcs[NUM_FUNCS];

  // Place functions into the func list
  funcs[0] = new TACSStructuralMass(tacs);
  funcs[1] = new TACSCompliance(tacs);

  // Set the discrete and continuous KS functions
  TACSKSFailure *func = new TACSKSFailure(tacs, 20.0);
  func->setKSFailureType(TACSKSFailure::DISCRETE);
  funcs[2] = func;

  func = new TACSKSFailure(tacs, 20.0);
  func->setKSFailureType(TACSKSFailure::CONTINUOUS);
  funcs[3] = func;

  // Set the induced norm failure types
  TACSInducedFailure * ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL);
  funcs[4] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL);
  funcs[5] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_EXPONENTIAL_SQUARED);
  funcs[6] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::EXPONENTIAL_SQUARED);
  funcs[7] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::POWER);
  funcs[8] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER);
  funcs[9] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::POWER_SQUARED);
  funcs[10] = ifunc;

  ifunc = new TACSInducedFailure(tacs, 20.0);
  ifunc->setInducedType(TACSInducedFailure::DISCRETE_POWER_SQUARED);
  funcs[11] = ifunc;
  
  for ( int k = 0; k < NUM_FUNCS; k++ ){
    funcs[k]->incref();
  }

  // Decrease the reference count to the mesh loader class
  mesh->decref();

  // Create the structural matrix/preconditioner/KSM to use in conjunction
  // with the aerostructural solution
  FEMat *mat = tacs->createFEMat();
  TACSBVec *ans = tacs->createVec();
  TACSBVec *res = tacs->createVec();
  ans->incref();
  res->incref();

  // Create the preconditioner
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->setMonitorFactorFlag(0);
  pc->setAlltoallAssemblyFlag(1);
  
  // Create GMRES object
  int gmres_iters = 4, nrestart = 0, isflexible = 0;
  GMRES *gmres = new GMRES(mat, pc, gmres_iters, nrestart, isflexible);
  gmres->incref();

  // Set the GMRES tolerances
  double rtol = 1e-12, atol = 1e-30;
  gmres->setTolerances(rtol, atol);

  // Now set up the information for the design variable sensitivities
  int num_design_vars = num_comp;

  // Alllocate the design variable vector/arrays
  TacsScalar *xvals = new TacsScalar[ num_design_vars ];
  TacsScalar *xtemp = new TacsScalar[ num_design_vars ];

  // Adjoint-based gradient evaluation
  TacsScalar *dfdx = new TacsScalar[ NUM_FUNCS*num_design_vars ];
  memset(dfdx, 0, NUM_FUNCS*num_design_vars*sizeof(TacsScalar));

  // Finite-difference gradient
  double *fd = new double[ NUM_FUNCS*num_design_vars ];
  memset(fd, 0, NUM_FUNCS*num_design_vars*sizeof(double));

  // Set the design variable values
  for ( int i = 0; i < num_design_vars; i++ ){
    xvals[i] = -1e20;
  }

  // Ensure consistency of the design variable values
  tacs->getDesignVars(xvals, num_design_vars);
  MPI_Allreduce(MPI_IN_PLACE, xvals, num_design_vars, TACS_MPI_TYPE,
                TACS_MPI_MAX, comm);

  // Set the design variables
  tacs->setDesignVars(xvals, num_design_vars);

  // The maximum number of gradient components to test
  // using finite-difference/complex-step
  int num_grad_comp = 10;
  num_grad_comp = (num_design_vars > num_grad_comp ? 
		   num_grad_comp : num_design_vars);

  // Scan any remaining arguments that may be required
  int test_flag = 0;
  double dh = 1e-6;
  for ( int k = 0; k < argc; k++ ){
    if (sscanf(argv[k], "dh=%lf", &dh) == 1){}
    else if (sscanf(argv[k], "num_grad_comp=%d", 
		    &num_grad_comp) == 1){
      num_grad_comp = (num_design_vars > num_grad_comp ? 
		       num_grad_comp : num_design_vars);
    }
    if (strcmp(argv[k], "test") == 0){
      test_flag = 1;
    }
  }

  // Print out the finite-difference interval
  if (rank == 0){
#ifdef TACS_USE_COMPLEX
    printf("Complex-step interval: %le\n", dh);
#else
    printf("Finite-difference interval: %le\n", dh);
#endif
  }

  // First solve the governing equations
  tacs->zeroVariables();
  tacs->assembleJacobian(1.0, 0.0, 0.0, res, mat);
  pc->factor();
  gmres->solve(res, ans);
  ans->scale(-1.0);
  tacs->setVariables(ans);
  
  // Write the solution to an .f5 file
  if (f5){
    // Find the location of the first '.' in the string
    char * outfile = new char[ strlen(bdf_file)+5 ];
    int len = strlen(bdf_file);
    int i = len-1;
    for ( ; i >= 0; i-- ){
      if (bdf_file[i] == '.'){ 
        break; 
      }     
    }
    
    // Copy the entire input string to the filename
    strcpy(outfile, bdf_file);
    
    // Over-write the last part of the string
    strcpy(&outfile[i], ".f5");
    
    // Write out the file
    f5->writeToFile(outfile);
    delete [] outfile;
  }

  // Evaluate each of the functions
  TacsScalar fvals[NUM_FUNCS];
  tacs->evalFunctions(funcs, NUM_FUNCS, fvals);
  
  // Evaluate the partial derivative of each function w.r.t. x
  tacs->addDVSens(1.0, funcs, NUM_FUNCS, 
                  dfdx, num_design_vars);
  
  // Create the adjoint variables for each function of interest
  TACSBVec *adjoints[NUM_FUNCS];
  TACSBVec *dfdu[NUM_FUNCS];
  for ( int j = 0; j < NUM_FUNCS; j++ ){
    adjoints[j] = tacs->createVec();
    adjoints[j]->incref();
    dfdu[j] = tacs->createVec();
    dfdu[j]->incref();
  }
  
  // Evaluate the partial derivative of each function of interest
  // w.r.t. the state variables and compute the adjoint
  tacs->addSVSens(1.0, 0.0, 0.0, funcs, NUM_FUNCS, dfdu);
  
  // Solve all of the adjoint equations
  for ( int j = 0; j < NUM_FUNCS; j++ ){
    gmres->solve(dfdu[j], adjoints[j]);
  }
  
  // Evaluate the adjoint-residual product
  tacs->addAdjointResProducts(-1.0, adjoints, NUM_FUNCS,
                              dfdx, num_design_vars);
  
  // Delete all of the adjoints
  for ( int j = 0; j < NUM_FUNCS; j++ ){
    adjoints[j]->decref();
    dfdu[j]->decref();
  }
  
  // Add up the contributions across all processors
  MPI_Allreduce(MPI_IN_PLACE, dfdx, NUM_FUNCS*num_design_vars,
                TACS_MPI_TYPE, MPI_SUM, comm);
  
  // Test the function that will be used
  if (test_flag){
    // Test the derivatives of the functions of interest w.r.t. the
    // state and design variables
    for ( int k = 0; k < NUM_FUNCS; k++ ){
      tacs->testFunction(funcs[k], num_design_vars, dh);
    }
  }

#if TACS_USE_COMPLEX
  // Compute the total derivative using the complex-step method
  for ( int k = 0; k < num_grad_comp; k++ ){
    // (xvals[k] + dh)
    memcpy(xtemp, xvals, num_design_vars*sizeof(TacsScalar));
    xtemp[k] += TacsScalar(0.0, dh);

    // Set the design variables
    tacs->setDesignVars(xtemp, num_design_vars);

    // Solve the problem
    tacs->zeroVariables();
    tacs->assembleJacobian(1.0, 0.0, 0.0, res, mat);
    pc->factor();
    gmres->solve(res, ans);
    ans->scale(-1.0);
    tacs->setVariables(ans);

    // Evaluate all of the functions of interest
    TacsScalar fvals[NUM_FUNCS];
    tacs->evalFunctions(funcs, NUM_FUNCS, fvals);

    // Compute their derivative based on the complex component
    // of the function value
    for ( int j = 0; j < NUM_FUNCS; j++ ){
      fd[k + j*num_design_vars] = TacsImagPart(fvals[j])/dh;
    }
  }

#else // !TACS_USE_COMPLEX
  // Test the structural sensitivities
  for ( int k = 0; k < num_grad_comp; k++ ){
    // (xvals[ok] + dh)
    memcpy(xtemp, xvals, num_design_vars*sizeof(TacsScalar));
    xtemp[k] += dh;
    
    // Set the perturbed value of the design variables
    tacs->setDesignVars(xtemp, num_design_vars);
      
    // Solve the problem
    tacs->zeroVariables();
    tacs->assembleJacobian(1.0, 0.0, 0.0, res, mat);
    pc->factor();
    gmres->solve(res, ans);
    ans->scale(-1.0);
    tacs->setVariables(ans);

    // Evaluate the function at the perturbed solution
    TacsScalar fvals[NUM_FUNCS];
    tacs->evalFunctions(funcs, NUM_FUNCS, fvals);
    for ( int j = 0; j < NUM_FUNCS; j++ ){
      fd[k + j*num_design_vars] = fvals[j];
    }

    // (xvals[k] - dh)
    memcpy(xtemp, xvals, num_design_vars*sizeof(TacsScalar));
    xtemp[k] -= dh;

    // Set the perturbed values of the design variables
    tacs->setDesignVars(xtemp, num_design_vars);

    // Solve the finite-element problem at the perturbed values
    // of the design variables
    tacs->zeroVariables();
    tacs->assembleJacobian(1.0, 0.0, 0.0, res, mat);
    pc->factor();
    gmres->solve(res, ans);
    ans->scale(-1.0);
    tacs->setVariables(ans);

    // Complete the finite-difference computation for each function
    tacs->evalFunctions(funcs, NUM_FUNCS, fvals);
    for ( int j = 0; j < NUM_FUNCS; j++ ){
      dfdx[k + j*num_design_vars] = 0.5*(dfdx[k + j*num_design_vars] -
                                         fvals[j])/dh;
    }
  }
#endif // TACS_USE_COMPLEX

  if (rank == 0){
    printf("Structural sensitivities\n");
    for ( int j = 0; j < NUM_FUNCS; j++ ){
      printf("Sensitivities for funcion %s\n",
             funcs[j]->functionName());
      printf("%25s %25s %25s\n",
             "Adjoint", "FD/CS", "Abs. error");
      for ( int k = 0; k < num_grad_comp; k++ ){
        printf("%25.15e %25.15e %25.15e\n", 
               TacsRealPart(dfdx[k + j*num_design_vars]),
               fd[k + j*num_design_vars],
               TacsRealPart(dfdx[k + j*num_design_vars]) -
               fd[k + j*num_design_vars]);
      }
    }
  }

  if (test_flag){
    int elem_num = 0;
    int print_level = 2;
    tacs->testElement(elem_num, print_level);
    tacs->testConstitutive(elem_num, print_level);
  }

  delete [] xvals;
  delete [] xtemp;
  delete [] dfdx;

  ans->decref();
  res->decref();
  tacs->decref();
  gmres->decref();

  for ( int k = 0; k < NUM_FUNCS; k++ ){
    funcs[k]->incref();
  }

  MPI_Finalize();
  return (0);
}
