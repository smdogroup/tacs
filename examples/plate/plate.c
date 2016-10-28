#include "TACSIntegrator.h"

#include "TACSMeshLoader.h"
#include "MITC9.h"
#include "isoFSDTStiffness.h"
#include "TACSShellTraction.h"

#include "KSFailure.h"
#include "StructuralMass.h"
#include "Compliance.h"

/*
  Code for testing adjoints with plate example. Use command line
  arguments as necessary.
 
  BDF1 BDF2 BDF3    : for BDF integrators
  DIRK2 DIRK3 DIRK4 : for DIRK integrators
  ABM1-6            : for ABM integrators
  NBG               : for Newmark integrator

  test_gradient : to do a complex step verification of the adjoint gradient
  test_element  : to test the element implementation
  num_funcs     : 1 to 3 for the adjoint 
  num_threads   : number of threads
  write_solution: write solution to f5 frequency
  print_level: 0, 1, 2
*/
int main( int argc, char **argv ){

  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank; 
  MPI_Comm_rank(comm, &rank); 
  
  // Parse command line arguments
  int num_funcs = 1;
  int num_threads = 1;
  int test_gradient = 0;
  int test_element = 0;
  int write_solution = 0;
  int print_level = 1;
  enum IntegratorType type = BDF2;
  for ( int i = 0; i < argc; i++ ){

    // Determine whether or not to test gradients with complex step
    if (strcmp("help", argv[i]) == 0){
      if (rank ==0){ 
        printf("TACS time-dependent analysis of a plate located in the folder as plate.bdf\n\n");
        printf("BDF1-3, DIRK2-4, ABM1-6, NBG : Selects the integrator to use\n");
        printf("test_gradient                : Complex-step verification of adjoint gradient\n");
        printf("num_funcs=1,2,3              : Number of functions for adjoint problem\n");
        printf("num_threads=1,2,3...         : Number of threads to use\n");
        printf("print_level=0,1,2            : Controls the amount of information to print\n");
        printf("write_solution=0,1,2...      : Controls the frequency of f5 file output\n\n");
      }
      MPI_Finalize();
      return 0;
    }
    
    // Backward Difference Formulae
    if (strcmp("BDF1", argv[i]) == 0){
      type = BDF1;
    } else if (strcmp("BDF2", argv[i]) == 0){
      type = BDF2;
    } else if (strcmp("BDF3", argv[i]) == 0){
      type = BDF3;

      // Adams-Bashforth-Moulton
    } else if (strcmp("ABM1", argv[i]) == 0){
      type = ABM1;
    } else if (strcmp("ABM2", argv[i]) == 0){
      type = ABM2;
    } else if (strcmp("ABM3", argv[i]) == 0){
      type = ABM3;
    } else if (strcmp("ABM4", argv[i]) == 0){
      type = ABM4;
    } else if (strcmp("ABM5", argv[i]) == 0){
      type = ABM5;
    } else if (strcmp("ABM6", argv[i]) == 0){
      type = ABM6;

      // Diagonally Implicit Runge Kutta
    } else if (strcmp("DIRK2", argv[i]) == 0){
      type = DIRK2;
    } else if (strcmp("DIRK3", argv[i]) == 0){
      type = DIRK3;
    } else if (strcmp("DIRK4", argv[i]) == 0){
      type = DIRK4;

      // Newmark-Beta-Gamma method
    } else if (strcmp("NBG", argv[i]) == 0){
      type = NBG;
    }

    // Determine the number of functions for adjoint
    if (sscanf(argv[i], "num_funcs=%d", &num_funcs) == 1){
      if (num_funcs < 0){ num_funcs = 1; }
      if (num_funcs > 3){ num_funcs = 3; }
      if (rank == 0){ printf("Number of functions : %d\n", num_funcs); }
    }

    // How frequent to write the f5 files
    if (sscanf(argv[i], "write_solution=%d", &write_solution) == 1){
      if (write_solution < 0){ write_solution = 0; }
      if (rank == 0){ printf("Write solution freq : %d\n", write_solution); }
    }

    // Set the print level
    if (sscanf(argv[i], "print_level=%d", &print_level) == 1){
      if (print_level < 0){ print_level = 1; }
      if (print_level > 3){ print_level = 3; }
      if (rank == 0){ printf("Print level : %d\n", print_level); }
    }

    // Determine the number of threads
    if (sscanf(argv[i], "num_threads=%d", &num_threads) == 1){
      if (num_threads < 0){ num_threads = 1; }
      if (num_threads > 24){ num_threads = 24; }
      if (rank ==0){ printf("Number of threads : %d\n", num_threads); }
    }

    // Determine whether or not to test gradients with complex step
    if (strcmp("test_gradient", argv[i]) == 0){
      test_gradient = 1;
      if (rank ==0){ printf("Enabled complex-step verification of gradients...\n"); }
    }
  }

  /*-----------------------------------------------------------------*/
  /*----------------Load Mesh and Setup TACS ------------------------*/
  /*-----------------------------------------------------------------*/

  // Write name of BDF file to be load to char array
  const char *filename = "plate.bdf";

  // Create the mesh loader object and load file
  TACSMeshLoader *mesh = new TACSMeshLoader(comm);
  mesh->incref();

  mesh->scanBDFFile(filename);

  // Get number of components prescribed in BDF file
  int num_components = mesh->getNumComponents();

  // Set properties for structural elements
  double rho   = 2500.0;  // density, kg/m^3
  double E     = 70e9;    // elastic modulus, Pa
  double nu    = 0.3;     // poisson's ratio
  double kcorr = 5.0/6.0; // shear correction factor
  double ys    = 350e6;   // yield stress, Pa

  // Set properties for dynamics
  TacsScalar g[]          = {0.0, 0.0, -9.81};
  TacsScalar v_init[]     = {0.0, 0.0, 1.e-2}; // not being used currently
  TacsScalar omega_init[] = {0.0, 0.0, 0.0};   // not being used currently

  /* TacsScalar v_init[] = {0.1, 0.1, 0.1};  */
  /* TacsScalar omega_init[] = {0.3, 0.1, 0.2}; */

  TACSGibbsVector *gravity = new TACSGibbsVector(g);  gravity->incref();
  TACSGibbsVector *v0      = new TACSGibbsVector(v_init); v0->incref();
  TACSGibbsVector *omega0  = new TACSGibbsVector(omega_init); omega0->incref();

  int vars_per_node;
  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ) {
    const char       *descriptor    = mesh->getElementDescript(i);
    double            min_thickness = 0.01;
    double            max_thickness = 0.1;
    double            thickness     = 0.05;
    isoFSDTStiffness *stiff         = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
						    thickness, i, 
						    min_thickness, max_thickness); 
    stiff->incref();
    
    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type defined in
    // the descriptor
    if (strcmp(descriptor, "CQUAD") == 0){
      element = new MITC9(stiff, gravity);
      element->incref();
    }
    else {
      printf("[%d] TACS Warning: Unsupported element %s in BDF file\n", rank, descriptor);
    }

    // Set the number of displacements
    vars_per_node = element->numDisplacements();
    mesh->setElement(i, element);

    stiff->decref();
    element->decref();
  }

  // Create tacs assembler from mesh loader object
  TACSAssembler *tacs = mesh->createTACS(vars_per_node);
  tacs->incref();

  /*-----------------------------------------------------------------*/
  /*-------------------------Setup Forces----------------------------*/
  /*-----------------------------------------------------------------*/
  
  /*
  // Create the traction
  TACSElement *trac = new TACSShellTraction<3>(0.0, 0.0, 0.01);
  trac->incref();
  
  // Create the auxiliary element class
  int nelems = tacs->getNumElements();
  TACSAuxElements *aux = new TACSAuxElements(nelems);
  aux->incref();

  // Associate the traction with each auxiliary element
  for ( int i = 0; i < nelems; i++ ){
  aux->addElement(i, trac);
  }
  
  // Set the auxiliary element in to TACS
  tacs->setAuxElements(aux);
  */

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);

  TACSToFH5 * f5 = new TACSToFH5(tacs, SHELL, write_flag);
  f5->incref();

  /*-----------------------------------------------------------------*/
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  int num_dvs = num_components;

  // Create functions of interest  
  TACSFunction * func[num_funcs];
  if (num_funcs == 1){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
  }
  else if (num_funcs == 2){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
    func[1] = new TACSCompliance(tacs); func[1]->incref();
  } 
  else if (num_funcs == 3){
    func[0] = new TACSKSFailure(tacs, 100.0); func[0]->incref();
    func[1] = new TACSCompliance(tacs); func[1]->incref();
    func[2] = new TACSStructuralMass(tacs); func[2]->incref();
  }

  TacsScalar *funcVals     = new TacsScalar[num_funcs]; // adjoint
  TacsScalar *funcValsTmp  = new TacsScalar[num_funcs]; // CSD
  TacsScalar *funcVals1    = new TacsScalar[num_funcs]; // forward/reverse

  TacsScalar *dfdx    = new TacsScalar[num_funcs*num_dvs]; // adjoint
  TacsScalar *dfdx1   = new TacsScalar[num_funcs*num_dvs]; // CSD
  TacsScalar *dfdxTmp = new TacsScalar[num_funcs*num_dvs]; // forward/reverse

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[ num_dvs ]; 
  x[0] = 0.03; 

  // Set paramters for time marching
  double tinit             = 0.0;
  double tfinal            = 100.0e-3; 
  int    num_steps_per_sec = 1000;

  TACSIntegrator *obj = TACSIntegrator::getInstance(tacs, tinit, tfinal, 
                                                    num_steps_per_sec, 
                                                    type);
  obj->incref();
  
  // Set options
  obj->setJacAssemblyFreq(1);
  obj->setAbsTol(1.0e-10);
  obj->setRelTol(1.0e-8);
  obj->setPrintLevel(print_level);
  obj->configureOutput(f5, write_solution, "output/plate_%04d.f5");
  
  // Set functions of interest for adjoint solve
  obj->setFunction(func, num_funcs);

  // Get the adjoint gradient
  double t0 = MPI_Wtime();
  obj->getFuncGrad(num_dvs, x, funcVals, dfdx);
  t0 = MPI_Wtime() - t0;

  // Print a summary of time taken
  if (rank == 0){
    obj->printWallTime(t0, 2);

    // Print the adjoint derivative values
    for( int j = 0; j < num_funcs; j++) {
      printf("[%d] Adj NEW  func: %d fval: %15.8e dfdx:", rank, j, RealPart(funcVals[j]));
      for ( int n = 0; n < num_dvs; n++) {
        printf(" %15.8e ",  RealPart(dfdx[n+j*num_dvs]));
      }
      printf("\n");
    }
    printf("\n");
  }

  // Test the adjoint derivatives if sought
  if (test_gradient){

    if ( rank == 0) { 
      printf("Finding complex-step gradient...\n");
    }

    // Complex step verification
    obj->getFDFuncGrad(num_dvs, x, funcValsTmp, dfdxTmp, 1.0e-16);

    if ( rank == 0) { 
      // Print complex step derivative values
      for( int j = 0; j < num_funcs; j++) {
        printf("[%d] CSD      func: %d fval: %15.8e dfdx:", rank, j, RealPart(funcValsTmp[j]));
        for ( int n = 0; n < num_dvs; n++) {
          printf(" %15.8e ",  RealPart(dfdxTmp[n+j*num_dvs]));
        }
        printf("\n");
      }
      printf("\n");

      // Print the differences between complex step and adjoint derivtives
      for ( int j = 0; j < num_funcs; j++ ) {
        printf("[%d] Error Adj NEW  func: %d ferror: %15.8e dfdx error:", rank, j, RealPart(funcValsTmp[j])-RealPart(funcVals[j]) );
        for ( int n = 0; n < num_dvs; n++ ) {
          printf(" %15.8e ",  RealPart(dfdxTmp[j*num_dvs+n]) -  RealPart(dfdx[j*num_dvs+n]) );
        }
        printf("\n");
      }
      printf("\n");
    }
  }

  obj->decref();

  mesh->decref();
  gravity->decref();
  v0->decref();
  omega0->decref();

  if (num_funcs == 1){
    func[0]->decref();
  }
  else if (num_funcs == 2){
    func[0]->decref();
    func[1]->decref();
  } 
  else if (num_funcs == 3){
    func[0]->decref();
    func[1]->decref();
    func[2]->decref();
  }

  tacs->decref();
  
  delete [] dfdx;
  delete [] dfdx1;
  delete [] dfdxTmp;

  delete [] x;

  delete [] funcVals;
  delete [] funcVals1;
  delete [] funcValsTmp;

  MPI_Finalize();
  return 0;
}

