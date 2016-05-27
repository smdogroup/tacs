#include "TACSIntegrator.h"
#include "TACSMeshLoader.h"
#include "MITC9.h"
#include "isoFSDTStiffness.h"
#include "KSFailure.h"
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

  int vars_per_node = 6;

  // Set up the global Gibbs vectors
  TacsScalar g[] = {0.0, 0.0, -9.81};
  TacsScalar v_init[] = {0.0, 0.0, 0.0};
  TacsScalar omega_init[] = {0.0, 0.0, 0.0};

  TACSGibbsVector *gravity = new TACSGibbsVector(g);
  TACSGibbsVector *v0 = new TACSGibbsVector(v_init);
  TACSGibbsVector *omega0 = new TACSGibbsVector(omega_init);

  // Loop over components, creating constituitive object for each
  for ( int i = 0; i < num_components; i++ ){
    const char *descriptor = mesh->getElementDescript(i);
    double min_thickness = 0.01;
    double max_thickness = 0.20;
    double thickness = 0.07;
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
  /*------------------ Time Integration and Adjoint Solve -----------*/
  /*-----------------------------------------------------------------*/

  TacsIntegrator *integrator = NULL;

  double tinit = 0.0, tfinal = 1.0;
  int num_steps_per_sec = 10, num_stages = 3, max_bdf_order = 3;

  // Create functions for adjoint solve
  TACSFunction *func = new KSFailure(tacs, 100.0, 1.0);

  /*
  // Integrate using DIRK
  integrator = new TacsDIRKIntegrator(tacs, tinit, tfinal,
  num_steps_per_sec, num_stages);
  integrator->incref();

  integrator->setPrintLevel(1);
  integrator->integrate();
  integrator->writeSolution("dirk.dat");
  integrator->setFunction(&func, 1);
  integrator->adjointSolve();

  integrator->decref();
  */
  
  // Integrate using BDF
  integrator = new TacsBDFIntegrator(tacs, tinit, tfinal, 
				     num_steps_per_sec, max_bdf_order);
  integrator->incref();

  integrator->setPrintLevel(1);
  integrator->integrate();
  integrator->writeSolution("bdf.dat");
  integrator->setFunction(&func, 1);
  integrator->adjointSolve();
  
  integrator->decref();
  

  tacs->decref();
  MPI_Finalize();

  return 0;
}
