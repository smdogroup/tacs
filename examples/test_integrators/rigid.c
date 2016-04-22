#include "TACSIntegrator.h"
#include "TACSAssembler.h"
#include "RigidBody.h"

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
  // bodyA->testJacobian(1e-6, alpha, beta, gamma);

  // Modify the inertial properties
  mass *= 2.0;
  J[0] += 1.0;
  TACSRigidBody *bodyB = new TACSRigidBody(mass, c, J);
  bodyB->incref();

  // Create the constraint
  TACSSphericalConstraint *con = new TACSSphericalConstraint();
  con->incref();

  // con->testJacobian(1e-6, alpha, beta, gamma);

  // Set up the TACSAssembler object
  int num_nodes = 3;
  int vars_per_node = 8;
  int num_elems = 3;
  int max_csr = 2;

  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD,
                                          num_nodes, vars_per_node,
                                          num_elems, num_nodes,
                                          max_csr);
  tacs->incref();

  // Add the nodes and the connectivity
  int conn[] = {0, 1, 2};
  tacs->addNode(0, 0);
  tacs->addNode(1, 1);
  tacs->addNode(2, 2);
  
  // Add the connectivity
  tacs->addElement(bodyA, &conn[0], 1);
  tacs->addElement(bodyB, &conn[1], 1);
  tacs->addElement(con, conn, 3);

  tacs->finalize();

  // Create the TACSIntegrator object
  double t_init = 0.0, t_final = 1.0;
  int steps_per_second = 100;
  int num_stages = 2;
  TacsDIRKIntegrator *dirk = new TacsDIRKIntegrator(tacs, t_init, t_final,
                                                    steps_per_second, num_stages);
  dirk->incref();

  TacsBDFIntegrator *bdf = new TacsBDFIntegrator(tacs, t_init, t_final,
                                                 steps_per_second, 2);
  bdf->incref();

  // bdf->integrate();

  // Integrate the equations of motion forward in time
  dirk->integrate();
  dirk->writeSolution("solution.dat");

  bdf->decref();
  dirk->decref();
  tacs->decref();
  bodyA->decref();
  bodyB->decref();
  con->decref();

  MPI_Finalize();
  return (0);
}
