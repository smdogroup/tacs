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

  test_element : for testing the element
*/
int main( int argc, char **argv ){

  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;

  // Parse command line arguments
  int test_element = 0;
  enum IntegratorType type = BDF1;
  for ( int i = 0; i < argc; i++ ){
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

    if (strcmp("test_element", argv[i]) == 0){
      test_element = 1;
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
  double rho = 2500.0;      // density, kg/m^3
  double E = 70e9;          // elastic modulus, Pa
  double nu = 0.3;          // poisson's ratio
  double kcorr = 5.0/6.0;   // shear correction factor
  double ys = 350e6;        // yield stress, Pa

  // Set properties for dynamics
  TacsScalar g[] = {0.0, 0.0, -9.81};
  TacsScalar v_init[] = {0.0, 0.0, 0.0};
  TacsScalar omega_init[] = {0.0, 0.0, 0.0};

  /* TacsScalar v_init[] = {0.1, 0.1, 0.1};  */
  /* TacsScalar omega_init[] = {0.3, 0.1, 0.2}; */

  TACSGibbsVector *gravity = new TACSGibbsVector(g); 
  gravity->incref();

  TACSGibbsVector *v0 = new TACSGibbsVector(v_init); 
  v0->incref();

  TACSGibbsVector *omega0 = new TACSGibbsVector(omega_init); 
  omega0->incref();

  int vars_per_node;
  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ) {
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.1;
    double thickness = 0.05;
    isoFSDTStiffness *stiff =  new isoFSDTStiffness(rho, E, nu, kcorr, ys,
						    thickness, i, 
						    min_thickness, max_thickness); 
    stiff->incref();
    
    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type defined in
    // the descriptor
    if (strcmp(descriptor, "CQUAD") == 0){
      element = new MITC9(stiff, gravity, v0, omega0);
      element->incref();
    }
    else {
      printf("TACS Warning: Unsupported element %s in BDF file\n", descriptor);
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

  if (test_element) {
    // Extract the element
    TacsScalar Xpts[3*9];
    TACSElement *elem = tacs->getElement(0, Xpts, NULL, NULL, NULL);
    elem->incref();

    // Test the element;
    TestElement *test = new TestElement(elem, Xpts);
    test->incref();
    test->setPrintLevel(2);
    test->testResidual();
    for ( int k = 0; k < elem->numVariables(); k++ ){
      test->testJacobian(k);
    }
    test->decref();
    elem->decref();
  }
  
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

  /*-----------------------------------------------------------------*/
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  int num_dvs = num_components;

  // Create functions of interest  
  static const int NUM_FUNCS = 3;
  TACSFunction * func[NUM_FUNCS];
  
  func[0] = new Compliance(tacs);
  func[0]->incref();
  
  func[1] = new StructuralMass(tacs);
  func[1]->incref();

  func[2] = new KSFailure(tacs, 100.0, 1.0);
  func[2]->incref();

  TacsScalar *funcVals     = new TacsScalar[NUM_FUNCS]; // adjoint
  TacsScalar *funcValsTmp  = new TacsScalar[NUM_FUNCS]; // CSD
  TacsScalar *funcVals1    = new TacsScalar[NUM_FUNCS]; // forward/reverse

  TacsScalar *dfdx = new TacsScalar[NUM_FUNCS*num_dvs];  // adjoint
  TacsScalar *dfdx1 = new TacsScalar[NUM_FUNCS*num_dvs]; // CSD
  TacsScalar *dfdxTmp = new TacsScalar[NUM_FUNCS*num_dvs]; // forward/reverse

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[ num_dvs ]; 
  x[0] = 0.03; 

  // Set paramters for time marching
  double tinit = 0.0, tfinal = 10.e-3; int num_steps_per_sec = 1000;

  TACSIntegrator *obj =  TACSIntegrator::getInstance(tacs, tinit, tfinal, 
                                                     num_steps_per_sec, 
                                                     type);
  obj->incref();

  // Set options
  obj->setJacAssemblyFreq(1);
  obj->setPrintLevel(1);
  obj->configureOutput(NULL, 1, "plate_%04d.f5");
  
  // Set functions of interest
  obj->setFunction(func, NUM_FUNCS);

  // COMPLEX STEP
  obj->getFDFuncGrad(num_dvs, x, funcValsTmp, dfdxTmp, 1.0e-8);

  // ADJOINT NEW
  obj->getFuncGrad(num_dvs, x, funcVals, dfdx);

  for( int j = 0; j < NUM_FUNCS; j++) {
    printf("CSD      func: %d fval: %15.8e dfdx:", j, RealPart(funcValsTmp[j]));
    for ( int n = 0; n < num_dvs; n++) {
      printf(" %15.8e ",  RealPart(dfdxTmp[n+j*num_dvs]));
    }
    printf("\n");
  }
  printf("\n");
 
  for( int j = 0; j < NUM_FUNCS; j++) {
    printf("Adj NEW  func: %d fval: %15.8e dfdx:", j, RealPart(funcVals[j]));
    for ( int n = 0; n < num_dvs; n++) {
      printf(" %15.8e ",  RealPart(dfdx[n+j*num_dvs]));
    }
    printf("\n");
  }
  printf("\n");

  // Error Summary
  for ( int j = 0; j < NUM_FUNCS; j++ ) {
    printf("Error Adj NEW  func: %d ferror: %15.8e dfdx error:", j, RealPart(funcValsTmp[j])-RealPart(funcVals[j]) );
    for ( int n = 0; n < num_dvs; n++ ) {
      printf(" %15.8e ",  RealPart(dfdxTmp[j*num_dvs+n]) -  RealPart(dfdx[j*num_dvs+n]) );
    }
    printf("\n");
  }
  printf("\n");

  /* // ADJOINT FORWARD REVERSE */
  /* for( int j = 0; j < NUM_FUNCS; j++) { */
  /*   funcVals1[j] = obj->forward(x, num_dvs, func[j]); */
  /*   obj->reverse(&dfdx1[j], num_dvs, func[j]); */
  /*   printf("BDF Adj FR   func: %d fval: %15.8e dfdx:", j, RealPart(funcVals1[j])); */
  /*   for ( int n = 0; n < num_dvs; n++) { */
  /*     printf(" %15.8e ",  RealPart(dfdx1[j*num_dvs+n])); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

  /* // Error Summary */
  /* for ( int j = 0; j < NUM_FUNCS; j++ ) { */
  /*   printf("Error Adj FR   func: %d ferror: %15.8e dfdx error:", j, RealPart(funcValsTmp[j])-RealPart(funcVals1[j]) ); */
  /*   for ( int n = 0; n < num_dvs; n++ ) { */
  /*     printf(" %15.8e ",  RealPart(dfdxTmp[j*num_dvs+n]) -  RealPart(dfdx1[j*num_dvs+n]) ); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("\n"); */

  obj->decref();

  mesh->decref();
  gravity->decref();
  v0->decref();
  omega0->decref();

  func[0]->decref();
  func[1]->decref();
  func[2]->decref();
 
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

