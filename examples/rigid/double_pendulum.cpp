#include "TACSAssembler.h"
#include "TACSIntegrator.h"
#include "TACSKinematicConstraints.h"
#include "TACSRigidBody.h"

/*
  Function to test the rigid body dynamics implementation
*/
int main(int argc, char *argv[]) {
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int use_revolute = 0;
  for (int k = 0; k < argc; k++) {
    if (strcmp(argv[k], "revolute") == 0) {
      use_revolute = 1;
    }
  }

  // The acceleration due to gravity in global frame of reference
  TACSGibbsVector *gravVec = new TACSGibbsVector(0.0, 0.0, -9.8);

  // Define the zero vector
  TACSGibbsVector *zero = new TACSGibbsVector(0.0, 0.0, 0.0);

  // Construct the frame of reference
  TACSGibbsVector *rA0Vec =
      new TACSGibbsVector(0.0, 0.0, 0.0);  // The base point
  TACSGibbsVector *rA1Vec =
      new TACSGibbsVector(1.0, 0.0, 0.0);  // The first coordinate
  TACSGibbsVector *rA2Vec =
      new TACSGibbsVector(0.0, 1.0, 0.0);  // The second coordinate
  TACSRefFrame *refFrameA = new TACSRefFrame(rA0Vec, rA1Vec, rA2Vec);

  // Define the inertial properties
  const TacsScalar mA = 1.0;
  const TacsScalar cA[3] = {0.0, 0.0, 0.0};
  const TacsScalar JA[6] = {1.0 / 3.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 1.0 / 3.0};

  // Define dynamics properties
  TACSGibbsVector *rAInitVec = new TACSGibbsVector(0.0, 2.5, 0.0);

  // Construct a rigid body
  TACSRigidBody *bodyA =
      new TACSRigidBody(refFrameA, mA, cA, JA, rAInitVec, zero, zero, gravVec);
  bodyA->incref();
  bodyA->setComponentNum(0);

  // Define the inertial properties
  const TacsScalar mB = 2.0;
  const TacsScalar cB[3] = {0.0, 0.0, 0.0};
  const TacsScalar JB[6] = {8.0 / 3.0, 0.0, 0.0, 8.0 / 3.0, 0.0, 8.0 / 3.0};

  // Define dynamics properties
  TACSGibbsVector *rBInitVec = new TACSGibbsVector(0.0, 5.5, 0.0);

  // Construct the second rigid body
  TACSRigidBody *bodyB =
      new TACSRigidBody(refFrameA, mB, cB, JB, rBInitVec, zero, zero, gravVec);
  bodyA->setComponentNum(1);
  bodyB->incref();

  // Create the constraint points
  TACSGibbsVector *basePt = new TACSGibbsVector(0.0, 0.0, 0.0);
  TACSGibbsVector *touchAB = new TACSGibbsVector(0.0, 5.0, 0.0);

  // Create a revolute axis
  TACSGibbsVector *rev = new TACSGibbsVector(1.0, 1.0, 1.0);

  // Set the constraints
  TACSElement *conA, *conB;
  if (use_revolute) {
    // Construct the revolute constraints
    conA = new TACSRevoluteConstraint(bodyA, basePt, rev);
    conB = new TACSRevoluteConstraint(bodyA, bodyB, touchAB, rev);
  } else {
    // Construct the spherical constraint
    conA = new TACSSphericalConstraint(bodyA, basePt);
    conB = new TACSSphericalConstraint(bodyA, bodyB, touchAB);
  }
  conA->incref();
  conB->incref();

  //---------------------------------------------------------------------//
  //                 Set up the TACSAssembler object                     //
  //---------------------------------------------------------------------//

  int num_nodes = 4;
  int vars_per_node = 8;
  int num_elems = 4;
  TACSAssembler *tacs =
      new TACSAssembler(MPI_COMM_WORLD, vars_per_node, num_nodes, num_elems);
  tacs->incref();

  // Set the elements
  TACSElement *elements[] = {bodyA, bodyB, conA, conB};
  tacs->setElements(elements);

  // Set the connectivity
  int conn[] = {0, 1, 0, 2, 0, 1, 3};
  int ptr[] = {0, 1, 2, 4, 7};
  tacs->setElementConnectivity(ptr, conn);
  tacs->initialize();

  tacs->testElement(2, 2);

  //--------------------------------------------------------------------//
  //                   Create the TACSIntegrator object                 //
  //--------------------------------------------------------------------//

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACS_OUTPUT_NODES | TACS_OUTPUT_DISPLACEMENTS);
  ElementType etype = TACS_BEAM_OR_SHELL_ELEMENT;  // How to set rigid type?
  TACSToFH5 *f5 = new TACSToFH5(tacs, etype, write_flag);
  f5->incref();

  double tinit = 0.0;
  double tfinal = 8.0;
  double steps_per_second = 50.0;
  int num_stages = 2;
  int max_bdf_order = 2;
  TACSBDFIntegrator *bdf = new TACSBDFIntegrator(
      tacs, tinit, tfinal, steps_per_second, max_bdf_order);
  bdf->incref();

  // Set optional parameters
  bdf->setRelTol(1.0e-8);
  bdf->setAbsTol(1.0e-12);
  bdf->setMaxNewtonIters(24);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(0);
  bdf->setOutputFrequency(1);

  // Integrate and write solution to file
  bdf->integrate();

  // Delete objects
  bodyA->decref();
  bodyB->decref();
  conA->decref();
  conB->decref();

  tacs->decref();
  f5->decref();
  bdf->decref();

  MPI_Finalize();
  return (0);
}
