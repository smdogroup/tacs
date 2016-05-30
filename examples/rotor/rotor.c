#include "TACSIntegrator.h"
#include "TACSMeshLoader.h"
#include "MITC9.h"
#include "MITCShell.h"
#include "isoFSDTStiffness.h"
#include "KSFailure.h"
#include "TACSShellTraction.h"

/*
  Function that sets up the rotor blade dynamics in TACS and
  integrates over time, solve the adjoint.
*/
int main( int argc, char **argv ){

  // Intialize MPI and declare communicator
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  /*-----------------------------------------------------------------*/
  /*----------------Load Mesh and Setup TACS ------------------------*/
  /*-----------------------------------------------------------------*/

  // Write name of BDF file to be load to char array
  const char *filename = "rotor.bdf";

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

  int vars_per_node = 6;

  // Set up the global Gibbs vectors
  TacsScalar g[] = {0.0, 0.0, -9.81};
  TacsScalar v_init[] = {0.0, 0.0, 0.0};
  TacsScalar omega_init[] = {0.0, 0.0, 0.0}; 
  TACSGibbsVector *gravity = new TACSGibbsVector(g); gravity->incref();
  TACSGibbsVector *v0 = new TACSGibbsVector(v_init); v0->incref();
  TACSGibbsVector *omega0 = new TACSGibbsVector(omega_init); omega0->incref();

  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ){
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.20;
    double thickness = 0.01;
    isoFSDTStiffness *stiff = new isoFSDTStiffness(rho, E, nu, kcorr, ys,
						   thickness, i, min_thickness, max_thickness); 
   
    // Initialize element object
    TACSElement *element = NULL;

    // Create element object using constituitive information and type defined in
    // the descriptor
    if ( strcmp(descriptor, "CQUAD") == 0 || 
         strcmp(descriptor, "CQUADR") == 0 ||
         strcmp(descriptor, "CQUAD4") == 0) {
      
      element = new MITCShell<2>(stiff, LINEAR, i);
      
      // if (strcmp(descriptor, "CQUAD") == 0){
      // element = new MITC9(stiff, gravity, v0, omega0);
      
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
  mesh->decref();

  /*
  // Extract the element
  TacsScalar Xpts[3*9];
  TACSElement *elem = tacs->getElement(0, Xpts, NULL, NULL, NULL);

  // Test the element;
  TestElement *test = new TestElement(elem, Xpts);
  test->setPrintLevel(2);
  test->testResidual();
  for ( int k = 0; k < elem->numVariables(); k++ ){
  test->testJacobian(k);
  }
  */

  /*-----------------------------------------------------------------*/
  /*-------------------------Setup Forces----------------------------*/
  /*-----------------------------------------------------------------*/

  // Create the traction
  TACSElement *trac = new TACSShellTraction<2>(1.0, 1.0, 1.0);
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

  TacsIntegrator *integrator = NULL;

  double tinit = 0.0, tfinal = 0.3;
  int num_steps_per_sec = 10, num_stages = 1, max_bdf_order = 2;

  // Create functions for adjoint solve
  TACSFunction *func = new KSFailure(tacs, 100.0, 1.0);

  // Evaluate the function
  TacsScalar *funcVals = new TacsScalar[1];
  tacs->evalFunctions( &func, 1, funcVals);

  printf("KSFailure:%f\n=", funcVals[0]);

  //*****************************************************************//
  // Integrate using DIRK
  //*****************************************************************//

  integrator = new TacsDIRKIntegrator(tacs, tinit, tfinal,
				      num_steps_per_sec, num_stages);
  integrator->incref();

  // Set optional parameters
  integrator->setRelTol(1.0e-7);
  integrator->setAbsTol(1.0e-12);
  integrator->setMaxNewtonIters(24);
  integrator->setPrintLevel(2);
  integrator->setJacAssemblyFreq(1);

  // Solve the equations over time
  printf(">> Integrating using DIRK\n");

  integrator->integrate();

  // Set the function and solve for the adjoint variables
  printf(">> Adjoint solve using DIRK\n");

  integrator->setFunction(&func, 1);
  integrator->adjointSolve();

  integrator->decref();

  //*****************************************************************//
  // Integrate using BDF
  //*****************************************************************//

  integrator = new TacsBDFIntegrator(tacs, tinit, tfinal, 
				     num_steps_per_sec, max_bdf_order);
  integrator->incref();
  
  func = new KSFailure(tacs, 100.0, 1.0);

  // Set optional parameters
  integrator->setRelTol(1.0e-7);
  integrator->setAbsTol(1.0e-12);
  integrator->setMaxNewtonIters(24);
  integrator->setPrintLevel(2);
  integrator->setJacAssemblyFreq(1);

  // Solve the equations over time
  printf(">> Integrating using BDF\n");

  integrator->integrate();

  // Set the function and solve for the adjoint variables
  printf(">> Adjoint solve using BDF\n");

  integrator->setFunction(&func, 1);
  integrator->adjointSolve();
  
  integrator->decref();

  // Delete objects
  trac->decref();
  aux->decref();

  gravity->decref();
  v0->decref();
  omega0->decref();
 
  tacs->decref();
  MPI_Finalize();

  return 0;
}
