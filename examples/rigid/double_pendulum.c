#include "TACSIntegrator.h"
#include "TACSAssembler.h"
#include "RigidBody.h"

/*
  Function to test the rigid body dynamics implementation
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // The acceleration due to gravity in global frame of reference
  const TacsScalar grav[3] = {0.0, 0.0, -10.0}; 
  TACSGibbsVector *gravVec = new TACSGibbsVector(grav); gravVec->incref();

  // Construct the frame of reference
  const TacsScalar rA0[3] = {0.0, 0.0, 0.0}; // The base point
  const TacsScalar rA1[3] = {0.5, 0.0, 0.0}; // The first coordinate direction
  const TacsScalar rA2[3] = {0.0, 0.5, 0.0}; // The second coordinate direction
  TACSGibbsVector *rA0Vec = new TACSGibbsVector(rA0); rA0Vec->incref();
  TACSGibbsVector *rA1Vec = new TACSGibbsVector(rA1); rA1Vec->incref();
  TACSGibbsVector *rA2Vec = new TACSGibbsVector(rA2); rA2Vec->incref();
  TACSRefFrame *refFrameA = new TACSRefFrame(rA0Vec, rA1Vec, rA2Vec);

  // Define the inertial properties
  const TacsScalar mA    = 1.0;
  const TacsScalar cA[3] = {0.5*1.0, 0.0, 0.0};
  const TacsScalar JA[6] = {0.025, 0.0, 0.0,
                            1.0/3.0, 0.0,
                            1.0/3.0};
  
  // Define dynamics properties
  const TacsScalar rAInit[3]     = {0.0, 0.0, 0.0}; // The initial position
  const TacsScalar vAInit[3]     = {0.0, 1.0, 0.0}; // The initial velocity
  const TacsScalar omegaAInit[3] = {0.0, 1.0, 0.0}; // The initial angular velocity
  TACSGibbsVector *rAInitVec     = new TACSGibbsVector(rAInit); rAInitVec->incref();
  TACSGibbsVector *vAInitVec     = new TACSGibbsVector(vAInit); vAInitVec->incref();
  TACSGibbsVector *omegaAInitVec = new TACSGibbsVector(omegaAInit); omegaAInitVec->incref();

  // Construct a rigid body
  TACSRigidBody *bodyA = new  TACSRigidBody(refFrameA,
                                            mA, cA, JA,
                                            rAInitVec, vAInitVec, omegaAInitVec, gravVec);
  
  //------------------------------------------------------------------//
  //                    Setup the second body                         //
  //------------------------------------------------------------------//

  // Construct the frame of reference
  const TacsScalar rB0[3] = {1.0, -1.0, 0.0}; // The base point
  const TacsScalar rB1[3] = {2.0, -1.0, 0.0}; // The first coordinate direction
  const TacsScalar rB2[3] = {1.0, 1.0, 0.0};  // The second coordinate direction
  TACSGibbsVector *rB0Vec = new TACSGibbsVector(rB0); rB0Vec->incref();
  TACSGibbsVector *rB1Vec = new TACSGibbsVector(rB1); rB1Vec->incref();
  TACSGibbsVector *rB2Vec = new TACSGibbsVector(rB2); rB2Vec->incref();
  TACSRefFrame *refFrameB = new TACSRefFrame(rB0Vec, rB1Vec, rB2Vec);

  // Define the inertial properties
  const TacsScalar mB    = 2.0;
  const TacsScalar cB[3] = {1.0*2.0, 0.0, 0.0};
  const TacsScalar JB[6] = {0.05, 0.0, 0.0,
                            8.0/3.0, 0.0,
                            8.0/3.0};
  
  // Define dynamics properties
  const TacsScalar rBInit[3]     = {0.0, 0.0, 0.0};   // The initial position
  const TacsScalar vBInit[3]     = {0.0, 1.0, 0.0};   // The initial velocity
  const TacsScalar omegaBInit[3] = {0.0, 1.0, 0.0};   // The initial angular velocity
  TACSGibbsVector *rBInitVec     = new TACSGibbsVector(rBInit); rBInitVec->incref();
  TACSGibbsVector *vBInitVec     = new TACSGibbsVector(vBInit); vBInitVec->incref();
  TACSGibbsVector *omegaBInitVec = new TACSGibbsVector(omegaBInit); omegaBInitVec->incref();
  
  TACSRigidBody *bodyB = new  TACSRigidBody(refFrameB, 
                                            mB, cB, JB,
                                            rBInitVec, vBInitVec, omegaBInitVec, gravVec);
  bodyB->incref();

  //-----------------------------------------------------------------------//
  //                   Setup Revolute Joint                                //
  //-----------------------------------------------------------------------//
 
  // Create the constraint
  TacsScalar revPoint[3] = {1.0, 0.0, 0.0};
  TacsScalar revAxis[3]  = {0.0, 1.0, 0.0};
  TACSGibbsVector *revPointVec = new TACSGibbsVector(revPoint); revPointVec->incref();
  TACSGibbsVector *revAxisVec  = new TACSGibbsVector(revAxis); revAxisVec->incref();

  // Construct the spherical constraint
  // TACSSphericalConstraint *con = new TACSSphericalConstraint(bodyA, bodyB, revPointVec);
  TACSRevoluteConstraint *con = new TACSRevoluteConstraint(bodyA, bodyB, revPointVec, revAxisVec);
  con->incref();

  //---------------------------------------------------------------------//
  //                 Set up the TACSAssembler object                     //
  //---------------------------------------------------------------------//

  int num_nodes     = 3;
  int vars_per_node = 8;
  int num_elems     = 3;

  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD, vars_per_node,
                                          num_nodes, num_elems);
  tacs->incref();
  
  // Set the elements
  TACSElement *elements[] = {bodyA, bodyB, con};
  tacs->setElements(elements);
  
  // Set the connectivity
  int conn[] = {0, 1, 0, 1, 2};
  int ptr[]  = {0, 1, 2, 5};
  tacs->setElementConnectivity(conn, ptr);

  tacs->initialize();

  //--------------------------------------------------------------------//
  //                   Create the TACSIntegrator object                 //
  //--------------------------------------------------------------------//
 // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS);
  TACSToFH5 * f5 = new TACSToFH5(tacs, RIGID, write_flag);
  f5->incref();

  double tinit = 0.0, tfinal = 0.1;
  int steps_per_second = 1000; int num_stages = 2;
  int max_bdf_order = 2;
  TACSBDFIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal,
                                                 steps_per_second, max_bdf_order);
  bdf->incref();
  
  // Set optional parameters
  bdf->setRelTol(1.0e-10);
  bdf->setAbsTol(1.0e-12);
  bdf->setMaxNewtonIters(24);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(0);
  bdf->configureOutput(f5, 1, "output/pendulum_%04d.f5");

  // Integrate and write solution to file
  bdf->integrate();
  bdf->writeSolution("solutionBDF.dat");
  bdf->decref();

  /*
  TACSDIRKIntegrator *dirk = 
  new TACSDIRKIntegrator(tacs, tinit, tfinal,
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
  */
  /*
  // Delete objects
  bodyA->decref();
  bodyB->decref();
  con->decref();

  r0Vec->decref();  
  r1Vec->decref();
  r2Vec->decref();

  gravVec->decref();
  rInitVec->decref();
  vInitVec->decref();

  tacs->decref();
  */
  MPI_Finalize();
  return (0);
}


/*

// bodyA->testResidual(1e-4);
  
// Modify the inertial properties
mass *= 2.0;
J[0] += 1.0;
TACSRigidBody *bodyB;// = new TACSRigidBody(mass, c, J);
bodyB->incref();

// Create the constraint
TACSSphericalConstraint *con;// = new TACSSphericalConstraint();
*/
//  con->incref();

//  TACSRevoluteConstraint *rev = new TACSRevoluteConstraint();

// TestElement *test = new TestElement(rev);
/*
  TestElement *test = new TestElement(rev);
  test->setPrintLevel(2);
  for ( int i = 0; i < 24; i++ ){
  test->testJacobian(i);
  }
*/
