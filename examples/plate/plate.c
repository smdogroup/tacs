#include "TACSIntegrator.h"
#include "TACSMeshLoader.h"
#include "MITC9.h"
#include "isoFSDTStiffness.h"
#include "KSFailure.h"
#include "TACSShellTraction.h"
#include "Compliance.h"
/*
  Function that sets up the rotor blade dynamics in TACS and
  integrates over time, solve the adjoint.

  ./plate BDF/DIRK
*/
int main( int argc, char **argv ){

  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  // Choice of integrator from command line
  int use_bdf = 0, test_element = 0;;
  for ( int i = 0; i < argc; i++ ){
    if (strcmp("BDF", argv[i]) == 0){
      use_bdf = 1;
    } else if (strcmp("test_element", argv[i]) == 0){
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

  // Set properties needed to create stiffness object 
  double rho = 2500.0;      // density, kg/m^3
  double E = 70e9;          // elastic modulus, Pa
  double nu = 0.3;          // poisson's ratio
  double kcorr = 5.0/6.0;   // shear correction factor
  double ys = 350e6;        // yield stress, Pa

  int vars_per_node = 8;

  // Set up the global Gibbs vectors
  TacsScalar g[] = {0.0, 0.0, -9.81};
  TacsScalar v_init[] = {0.25, 0.25, 0.25};
  TacsScalar omega_init[] = {0.0, 0.0, 0.0};

  TACSGibbsVector *gravity = new TACSGibbsVector(g); gravity->incref();
  TACSGibbsVector *v0 = new TACSGibbsVector(v_init); v0->incref();
  TACSGibbsVector *omega0 = new TACSGibbsVector(omega_init); omega0->incref();

  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ) {
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.1;
    double thickness = 0.05;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
						   thickness, i, min_thickness, max_thickness); 

    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type defined in
    // the descriptor
    if (strcmp(descriptor, "CQUAD") == 0){
      element = new MITC9(stiff, gravity, v0, omega0);
    }
    else {
      printf("TACS Warning: Unsupported element %s in BDF file\n", descriptor);
    }

    // Set the number of displacements
    vars_per_node = element->numDisplacements();
    mesh->setElement(i, element);
  }

  // Create tacs assembler from mesh loader object
  TACSAssembler *tacs = mesh->createTACS(vars_per_node);
  tacs->incref();

  if (test_element) {
    // Extract the element
    TacsScalar Xpts[3*9];
    TACSElement *elem = tacs->getElement(0, Xpts, NULL, NULL, NULL);
    
    // Test the element;
    TestElement *test = new TestElement(elem, Xpts);
    test->incref();
    test->setPrintLevel(2);
    test->testResidual();
    for ( int k = 0; k < elem->numVariables(); k++ ){
      test->testJacobian(k);
    }
    test->decref();
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

  /*-----------------------------------------------------------------*/
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  int num_dvs = num_components;
  
  static const int NUM_FUNCS = 1;
  TACSFunction * func[NUM_FUNCS];
  
  //  func[0] = new KSFailure(tacs, 100.0, 1.0);
  func[0] = new Compliance(tacs);
  func[0]->incref();
   
  // Evaluate the function
  TacsScalar *funcVals = new TacsScalar[NUM_FUNCS];

  TacsScalar *dfdx = new TacsScalar[NUM_FUNCS*num_dvs];
  memset(dfdx, 0, NUM_FUNCS*num_dvs*sizeof(TacsScalar));
  
  TacsScalar *dfdxTmp = new TacsScalar[NUM_FUNCS*num_dvs];
  memset(dfdxTmp, 0, NUM_FUNCS*num_dvs*sizeof(TacsScalar));

  // Create an array of design variables
  TacsScalar *x = new TacsScalar[ num_dvs ];
  x[0] = 0.03; 

  TacsIntegrator *integrator = NULL;

  double tinit = 0.0, tfinal = 0.02;
  int num_steps_per_sec = 1000;

  int num_stages = 3, max_bdf_order = 3;
  
   if (!use_bdf) { 

    integrator = new TacsDIRKIntegrator( tacs, tinit, tfinal, num_steps_per_sec, num_stages ); 
    integrator->incref(); 

    integrator->integrate();
   
  } else { 

    integrator = new TacsBDFIntegrator( tacs, tinit, tfinal, num_steps_per_sec, max_bdf_order );
    integrator->incref();
    
    integrator->setJacAssemblyFreq(1);

    integrator->getAdjointGradient(func, NUM_FUNCS, num_dvs, x, funcVals, dfdx);
    printf("Compliance = %e\n", RealPart(funcVals[0])); 

    integrator->getApproxGradient(func, NUM_FUNCS, num_dvs, x, funcVals, dfdxTmp, 1.0e-8);
    printf("Compliance = %e\n", RealPart(funcVals[0])); 

    for ( int i = 0; i < num_dvs; i++) {
      printf("dfdx[%d] = %e %e %e \n", i, RealPart(dfdx[i]), RealPart(dfdxTmp[i]), (RealPart(dfdx[i]) - RealPart(dfdxTmp[i])));
    } 

  }

  integrator->decref();
  func[0]->decref();

  delete [] dfdx;
  delete [] dfdxTmp;
  delete [] x;
  delete [] funcVals;
 
  gravity->decref();
  v0->decref();
  omega0->decref();
  mesh->decref();
  tacs->decref();

  // aux->decref();
  // trac->decref();
  
  MPI_Finalize();

  return 0;
}
