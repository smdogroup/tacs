#include "TACSIntegrator.h"
#include "TACSAssembler.h"
#include "RigidBody.h"

/*
  Function to test the rigid body dynamics implementation
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  TacsScalar mass = 4.0;
  TacsScalar c[] = {0.5, 0.2, -0.1};
  TacsScalar J[] = {1.0, -0.1, 0.25,
                    2.0, 0.1,
                    0.75};

  // Zero the Jacobian coefficients
  double alpha = 1.26, beta = 0.35, gamma = 4.34;

  TACSRigidBody *bodyA = new TACSRigidBody(mass, c, J);
  bodyA->incref();
  
  // bodyA->testResidual(1e-4);

  // Modify the inertial properties
  mass *= 2.0;
  J[0] += 1.0;
  TACSRigidBody *bodyB = new TACSRigidBody(mass, c, J);
  bodyB->incref();

  // Create the constraint
  TACSSphericalConstraint *con = new TACSSphericalConstraint();
  con->incref();

  TACSRevoluteConstraint *rev = new TACSRevoluteConstraint();

  // TestElement *test = new TestElement(rev);
  /*
    TestElement *test = new TestElement(rev);
    test->setPrintLevel(2);
    for ( int i = 0; i < 24; i++ ){
    test->testJacobian(i);
    }
  */

  // Set up the TACSAssembler object
  int num_nodes = 3;
  int vars_per_node = 8;
  int num_elems = 3;

  TACSAssembler *tacs = 
    new TACSAssembler(MPI_COMM_WORLD, vars_per_node,
                      num_nodes, num_elems);
  tacs->incref();

  // Set the elements
  TACSElement *elements[] = {bodyA, bodyB, con};
  tacs->setElements(elements);

  // Set the connectivity
  int conn[] = {0, 1, 0, 1, 2};
  int ptr[] = {0, 1, 2, 5};
  tacs->setElementConnectivity(conn, ptr);
  tacs->initialize();

  // Create the TACSIntegrator object
  double t_init = 0.0, t_final = 1.0;
  int steps_per_second = 100; int num_stages = 2;
  TACSDIRKIntegrator *dirk = 
    new TACSDIRKIntegrator(tacs, t_init, t_final,
                           steps_per_second, num_stages);
  dirk->incref();

  // Set optional parameters
  dirk->setRelTol(1.0e-10);
  dirk->setAbsTol(1.0e-14);
  dirk->setMaxNewtonIters(24);
  dirk->setPrintLevel(1);
  dirk->setJacAssemblyFreq(1);
  dirk->setUseLapack(0);
  
  // Integrate and write solution to file
  dirk->integrate();
  dirk->writeSolution("solutionDIRK.dat");
  dirk->decref();
  
  TACSBDFIntegrator *bdf = 
    new TACSBDFIntegrator(tacs, t_init, t_final,
                          steps_per_second, 2);
  bdf->incref();

  // Set optional parameters
  bdf->setRelTol(1.0e-10);
  bdf->setAbsTol(1.0e-14);
  bdf->setMaxNewtonIters(24);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(0);

  // Integrate and write solution to file
  bdf->integrate();
  bdf->writeSolution("solutionBDF.dat");
  bdf->decref();

  // Delete objects
  tacs->decref();
  bodyA->decref();
  bodyB->decref();
  con->decref();

  MPI_Finalize();
  return (0);
}
