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
  const TacsScalar grav[3] = {0.0, 0.0, -9.8};
  TACSGibbsVector *gravVec = new TACSGibbsVector(grav); gravVec->incref();

  // Construct the frame of reference
  const TacsScalar rA0[3] = {0.0, 0.0, 0.0}; // The base point
  const TacsScalar rA1[3] = {1.0, 0.0, 0.0}; // The first coordinate direction
  const TacsScalar rA2[3] = {0.0, 1.0, 0.0}; // The second coordinate direction
  TACSGibbsVector *rA0Vec = new TACSGibbsVector(rA0); rA0Vec->incref();
  TACSGibbsVector *rA1Vec = new TACSGibbsVector(rA1); rA1Vec->incref();
  TACSGibbsVector *rA2Vec = new TACSGibbsVector(rA2); rA2Vec->incref();
  TACSRefFrame *refFrameA = new TACSRefFrame(rA0Vec, rA1Vec, rA2Vec);
  refFrameA->incref();

  // Define the inertial properties
  const TacsScalar mA    = 1.0;
  const TacsScalar cA[3] = {0.0, 0.0, 0.0};
  const TacsScalar JA[6] = {1.0/3.0, 0.0, 0.0,
                            1.0/3.0, 0.0,
                            1.0/3.0};
  
  // Define dynamics properties
  const TacsScalar rAInit[3]     = {0.0, 2.0, 0.0}; // The initial position
  const TacsScalar vAInit[3]     = {0.0, 0.0, 0.0}; // The initial velocity
  const TacsScalar omegaAInit[3] = {0.0, 0.0, 0.0}; // The initial angular velocity
  TACSGibbsVector *rAInitVec     = new TACSGibbsVector(rAInit); rAInitVec->incref();
  TACSGibbsVector *vAInitVec     = new TACSGibbsVector(vAInit); vAInitVec->incref();
  TACSGibbsVector *omegaAInitVec = new TACSGibbsVector(omegaAInit); omegaAInitVec->incref();

  // Create visualization
  TACSRigidBodyViz *vizA = new TACSRigidBodyViz(0.5, 5.0, 0.5);
  vizA->incref();

  // Construct a rigid body
  TACSRigidBody *bodyA = new  TACSRigidBody(refFrameA,
                                            mA, cA, JA,
                                            rAInitVec, vAInitVec, omegaAInitVec, gravVec);
  bodyA->setVisualization(vizA);

  bodyA->incref();

  //------------------------------------------------------------------//
  //                    Setup the second body                         //
  //------------------------------------------------------------------//

  // Construct the frame of reference
  /* const TacsScalar rB0[3] = {0.0, 1.0, 0.0};  // The base point */
  /* const TacsScalar rB1[3] = {2.0, -1.0, 0.0}; // The first coordinate direction */
  /* const TacsScalar rB2[3] = {1.0, 1.0, 0.0};  // The second coordinate direction */
  /* TACSGibbsVector *rB0Vec = new TACSGibbsVector(rB0); rB0Vec->incref(); */
  /* TACSGibbsVector *rB1Vec = new TACSGibbsVector(rB1); rB1Vec->incref(); */
  /* TACSGibbsVector *rB2Vec = new TACSGibbsVector(rB2); rB2Vec->incref(); */
  /* TACSRefFrame *refFrameB = new TACSRefFrame(rB0Vec, rB1Vec, rB2Vec); */
  /* refFrameB->incref(); */

  // Define the inertial properties
  const TacsScalar mB    = 2.0;
  const TacsScalar cB[3] = {0.0, 0.0, 0.0};
  const TacsScalar JB[6] = {8.0/3.0, 0.0, 0.0,
                            8.0/3.0, 0.0,
                            8.0/3.0};

  // Define dynamics properties
  const TacsScalar rBInit[3]     = {0.0, 5.0, 1.0};   // The initial position
  const TacsScalar vBInit[3]     = {0.0, 0.0, 0.0};   // The initial velocity
  const TacsScalar omegaBInit[3] = {0.0, 0.0, 0.0};   // The initial angular velocity
  TACSGibbsVector *rBInitVec     = new TACSGibbsVector(rBInit); rBInitVec->incref();
  TACSGibbsVector *vBInitVec     = new TACSGibbsVector(vBInit); vBInitVec->incref();
  TACSGibbsVector *omegaBInitVec = new TACSGibbsVector(omegaBInit); omegaBInitVec->incref();

  // Create visualization
  TACSRigidBodyViz *vizB = new TACSRigidBodyViz(1.0);
  vizB->incref();
  
  TACSRigidBody *bodyB = new  TACSRigidBody(refFrameA, 
                                            mB, cB, JB,
                                            rBInitVec, vBInitVec, omegaBInitVec, gravVec);
  bodyB->incref();

  bodyB->setVisualization(vizB);

  //-----------------------------------------------------------------------//
  //                   Setup Revolute Joint                                //
  //-----------------------------------------------------------------------//
 
  // Create the constraint
  TacsScalar revPoint[3] = {0.0, 1.0, 0.0};
  TacsScalar revAxis[3]  = {0.0, 1.0, 0.0};
  TACSGibbsVector *revPointVec = new TACSGibbsVector(revPoint); revPointVec->incref();
  TACSGibbsVector *revAxisVec  = new TACSGibbsVector(revAxis); revAxisVec->incref();

  // Construct the spherical constraint
  TACSSphericalConstraint *con = new TACSSphericalConstraint(bodyA, bodyB, revPointVec);
  // TACSRevoluteConstraint *con = new TACSRevoluteConstraint(bodyA, bodyB, revPointVec, revAxisVec);
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
  TACSElement *elements[] = {con, bodyA, bodyB};
  tacs->setElements(elements);
  
  // Set the connectivity
  int conn[] = {0, 1, 2, 0, 1};
  int ptr[]  = {0, 3, 4, 5};
  tacs->setElementConnectivity(conn, ptr);

  tacs->initialize();

  //--------------------------------------------------------------------//
  //                   Create the TACSIntegrator object                 //
  //--------------------------------------------------------------------//

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, RIGID, write_flag);
  f5->incref();

  double tinit            = 0.0;
  double tfinal           = 1.0;
  int    steps_per_second = 200; 
  int    num_stages       = 2;
  int    max_bdf_order    = 2;
  TACSBDFIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal,
                                                 steps_per_second, 
                                                 max_bdf_order);
  bdf->incref();
  
  // Set optional parameters
  bdf->setRelTol(1.0e-10);
  bdf->setAbsTol(1.0e-12);
  bdf->setMaxNewtonIters(24);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(0);
  bdf->configureOutput(f5, 1, "double-pendulum-output/pendulum_%04d.f5");

  // Integrate and write solution to file
  bdf->integrate();
  bdf->writeSolution("solutionBDF.dat");

  // Delete objects
  rA0Vec->decref();  
  rA1Vec->decref();
  rA2Vec->decref();
  refFrameA->decref();

  gravVec->decref();
  rAInitVec->decref();
  vAInitVec->decref();
  omegaAInitVec->decref();

  /* rB0Vec->decref();   */
  /* rB1Vec->decref(); */
  /* rB2Vec->decref(); */
  /* refFrameB->decref(); */

  rBInitVec->decref();
  vBInitVec->decref();
  omegaBInitVec->decref();

  bodyA->decref();
  vizA->decref();

  bodyB->decref();
  vizB->decref();

  con->decref();

  tacs->decref();
  f5->decref();

  bdf->decref();

  MPI_Finalize();
  return (0);
}
