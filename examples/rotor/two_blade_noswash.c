#include "TACSIntegrator.h"
#include "TACSAssembler.h"
#include "RigidBody.h"

/*
  Creates a rotor assembly example using rigid bodies and constraints
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // The acceleration due to gravity in global frame of reference
  TACSGibbsVector *gravVec = new TACSGibbsVector(0.0, 0.0, -10.0);
  gravVec->incref();
  
  TACSGibbsVector *omega = new TACSGibbsVector(0.0, 0.0, 1.0);
  omega->incref();

  // Define the zero vector
  TACSGibbsVector *zero = new TACSGibbsVector(0.0, 0.0, 0.0);
  zero->incref();

  // Define the inertial properties of the shaft
  const TacsScalar mA    = 1.0;
  const TacsScalar cA[3] = {0.0, 0.0, 0.0};
  const TacsScalar JA[6] = {1.0/3.0, 0.0, 0.0,
                            1.0/3.0, 0.0,
                            1.0/3.0};

  const double pitch_offset = 0.3/4.0;

  //------------------------------------------------------------------//
  //                          Shaft                                   //
  //------------------------------------------------------------------//
  
  // Shaft is to be 2m long

  TACSGibbsVector  *rShaft   = new TACSGibbsVector(0.0, 0.0, 1.0); 
  rShaft->incref();
  TACSGibbsVector  *rShaftx  = new TACSGibbsVector(1.0, 0.0, 1.0); 
  rShaftx->incref();
  TACSGibbsVector  *rShafty  = new TACSGibbsVector(0.0, 1.0, 1.0); 
  rShafty->incref();
  TACSRefFrame     *fShaft   = new TACSRefFrame(rShaft, rShaftx, rShafty); 
  fShaft->incref();
  TACSRigidBodyViz *vizShaft = new TACSRigidBodyViz(0.1, 0.1, 2.0);
  vizShaft->incref();
  TACSRigidBody    *shaft    = new TACSRigidBody(fShaft,
                                                 mA, cA, JA,
                                                 rShaft, zero, omega, gravVec);
  shaft->incref();
  shaft->setVisualization(vizShaft);

  //------------------------------------------------------------------//
  //                 Upper swash plate (USP)
  //------------------------------------------------------------------//
  
  // USP is located 1.5m from the bottom of the shaft

  TACSGibbsVector  *rUSP   = new TACSGibbsVector(0.0, 0.0, 1.5); 
  rUSP->incref();
  TACSGibbsVector  *rUSPx  = new TACSGibbsVector(1.0, 0.0, 1.5);
  rUSPx->incref();
  TACSGibbsVector  *rUSPy  = new TACSGibbsVector(0.0, 1.0, 1.5);
  rUSPy->incref();
  TACSRefFrame     *fUSP   = new TACSRefFrame(rUSP, rUSPx, rUSPy);
  fUSP->incref();
  TACSRigidBodyViz *vizUSP = new TACSRigidBodyViz(0.3, 0.3, 0.1);
  vizUSP->incref();
  TACSRigidBody    *usp    = new TACSRigidBody(fUSP,
                                               mA, cA, JA,
                                               rUSP, zero, omega, gravVec);
  usp->incref();
  usp->setVisualization(vizUSP);

  //------------------------------------------------------------------//
  //                 Blade 3 (into the plane)
  //------------------------------------------------------------------//

  // blade is 5m long
  
  TACSGibbsVector  *rBlade3   = new TACSGibbsVector(0.0, 2.5, 2.0);
  rBlade3->incref();
  TACSGibbsVector  *rBlade3x  = new TACSGibbsVector(1.0, 2.5, 2.0);
  rBlade3x->incref();
  TACSGibbsVector  *rBlade3y  = new TACSGibbsVector(0.0, 2.5 + 1.0, 2.0);
  rBlade3y->incref();
  TACSRefFrame     *fBlade3   = new TACSRefFrame(rBlade3, rBlade3x, rBlade3y);
  fBlade3->incref();
  TACSRigidBodyViz *vizBlade3 = new TACSRigidBodyViz(0.2, 5.0, 0.05);
  vizBlade3->incref();
  TACSRigidBody    *blade3    = new TACSRigidBody(fBlade3,
                                                  mA, cA, JA,
                                                  rBlade3, zero, omega, gravVec);
  blade3->incref();
  blade3->setVisualization(vizBlade3);

  //------------------------------------------------------------------//
  //                 Pitch Link 3 (into the plane)
  //------------------------------------------------------------------//
  
  // pitch link is half a meter long

  TACSGibbsVector  *rPitch3   = new TACSGibbsVector(-pitch_offset + 0.0, 0.15, 1.75);
  rPitch3->incref();
  TACSGibbsVector  *rPitch3x  = new TACSGibbsVector(-pitch_offset + 1.0, 0.15, 1.75);
  rPitch3x->incref();
  TACSGibbsVector  *rPitch3y  = new TACSGibbsVector(-pitch_offset + 0.0, 1.15, 1.75);
  rPitch3y->incref();
  TACSRefFrame     *fPitch3   = new TACSRefFrame(rPitch3, rPitch3x, rPitch3y);
  fPitch3->incref();
  TACSRigidBodyViz *vizPitch3 = new TACSRigidBodyViz(0.02, 0.02, 0.5);
  vizPitch3->incref();
  TACSRigidBody    *pitch3    = new TACSRigidBody(fPitch3,
                                                  mA, cA, JA,
                                                  rPitch3, zero, omega, gravVec);
  pitch3->incref();
  pitch3->setVisualization(vizPitch3);

  //------------------------------------------------------------------//
  //                 Blade 4 (out of paper)
  //------------------------------------------------------------------//

  // blade is 5m long
  
  TACSGibbsVector  *rBlade4   = new TACSGibbsVector(0.0, -2.5, 2.0);  
  rBlade4->incref();
  TACSGibbsVector  *rBlade4x  = new TACSGibbsVector(1.0, -2.5, 2.0);
  rBlade4x->incref();
  TACSGibbsVector  *rBlade4y  = new TACSGibbsVector(0.0, -2.5 + 1.0, 2.0);
  rBlade4y->incref();
  TACSRefFrame     *fBlade4   = new TACSRefFrame(rBlade4, rBlade4x, rBlade4y);
  fBlade4->incref();
  TACSRigidBodyViz *vizBlade4 = new TACSRigidBodyViz(0.2, 5.0, 0.05);
  vizBlade4->incref();
  TACSRigidBody    *blade4    = new TACSRigidBody(fBlade4,
                                                  mA, cA, JA,
                                                  rBlade4, zero, omega, gravVec);
  blade4->incref();
  blade4->setVisualization(vizBlade4);

  //------------------------------------------------------------------//
  //                 Pitch Link 4 (out of paper)
  //------------------------------------------------------------------//
  
  // pitch link is half a meter long

  TACSGibbsVector  *rPitch4   = new TACSGibbsVector(pitch_offset + 0.0, -0.15, 1.75);
  rPitch4->incref();
  TACSGibbsVector  *rPitch4x  = new TACSGibbsVector(pitch_offset + 1.0, -0.15, 1.75);
  rPitch4x->incref();
  TACSGibbsVector  *rPitch4y  = new TACSGibbsVector(pitch_offset + 0.0, -0.15 + 1.0, 1.75);
  rPitch4y->incref();
  TACSRefFrame     *fPitch4   = new TACSRefFrame(rPitch4, rPitch4x, rPitch4y);
  fPitch4->incref();
  TACSRigidBodyViz *vizPitch4 = new TACSRigidBodyViz(0.02, 0.02, 0.5);
  vizPitch4->incref();
  TACSRigidBody    *pitch4    = new TACSRigidBody(fPitch4,
                                                  mA, cA, JA,
                                                  rPitch4, zero, omega, gravVec);
  pitch4->incref();
  pitch4->setVisualization(vizPitch4);

  //------------------------------------------------------------------//
  //              Constraints
  //------------------------------------------------------------------//

  TACSGibbsVector         *xrev           = new TACSGibbsVector(1.0, 0.0, 0.0);
  xrev->incref();

  TACSGibbsVector         *yrev           = new TACSGibbsVector(0.0, 1.0, 0.0);
  yrev->incref();

  TACSGibbsVector         *zrev           = new TACSGibbsVector(0.0, 0.0, 1.0);
  zrev->incref();
    
  // shaft with base
  TACSGibbsVector         *pnt_shaft_base = new TACSGibbsVector(0.0, 0.0, 0.0);
  pnt_shaft_base->incref();
  TACSRevoluteConstraint *rev_shaft_base = new TACSRevoluteConstraint(shaft, pnt_shaft_base, zrev);
  rev_shaft_base->incref();

  // shaft with upper swash plate
  TACSGibbsVector         *pnt_shaft_usp = new TACSGibbsVector(0.0, 0.0, 1.5);
  pnt_shaft_usp->incref();
  TACSSphericalConstraint *sph_shaft_usp = new TACSSphericalConstraint(shaft, usp, pnt_shaft_usp);
  sph_shaft_usp->incref();

  // inner blade 3 (along +y direction)
  TACSGibbsVector         *pnt_shaft_blade3  = new TACSGibbsVector(0.0, 0.0, 2.0);
  pnt_shaft_blade3->incref();
  TACSRevoluteConstraint *rev_shaft_blade3  = new TACSRevoluteConstraint(shaft, blade3, pnt_shaft_blade3, yrev);
  rev_shaft_blade3->incref();
  TACSGibbsVector         *pnt_blade3_pitch3 = new TACSGibbsVector(-pitch_offset, 0.15, 2.0);
  pnt_blade3_pitch3->incref();
  TACSSphericalConstraint *rev_blade3_pitch3 = new TACSSphericalConstraint(blade3, pitch3, pnt_blade3_pitch3);
  rev_blade3_pitch3->incref();
  TACSGibbsVector         *pnt_pitch3_usp    = new TACSGibbsVector(-pitch_offset, 0.15, 1.5);
  pnt_pitch3_usp->incref();
  TACSSphericalConstraint *rev_pitch3_usp    = new TACSSphericalConstraint(pitch3, usp, pnt_pitch3_usp);
  rev_pitch3_usp->incref();

  // outer blade 4  (along -y direction)
  TACSGibbsVector         *pnt_shaft_blade4  = new TACSGibbsVector(0.0, 0.0, 2.0);
  pnt_shaft_blade4->incref();
  TACSRevoluteConstraint *rev_shaft_blade4  = new TACSRevoluteConstraint(shaft, blade4, pnt_shaft_blade4, yrev);
  rev_shaft_blade4->incref();
  TACSGibbsVector         *pnt_blade4_pitch4 = new TACSGibbsVector(pitch_offset, -0.15, 2.0);
  pnt_blade4_pitch4->incref();
  TACSSphericalConstraint *rev_blade4_pitch4 = new TACSSphericalConstraint(blade4, pitch4, pnt_blade4_pitch4);
  rev_blade4_pitch4->incref();
  TACSGibbsVector         *pnt_pitch4_usp    = new TACSGibbsVector(pitch_offset, -0.15, 1.5);
  pnt_pitch4_usp->incref();
  TACSSphericalConstraint *rev_pitch4_usp    = new TACSSphericalConstraint(pitch4, usp, pnt_pitch4_usp);
  rev_pitch4_usp->incref();

  //------------------------------------------------------------------//
  //              Set up the TACSAssembler object                     //
  //------------------------------------------------------------------//

  int num_nodes     = 12;  // Number of finite element nodes
  int vars_per_node = 8;   // 
  int num_elems     = 12;

  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD, vars_per_node,
                                          num_nodes, num_elems);
  tacs->incref();

  /*
  // Set all the component numbers for visualization
  shaft->setComponentNum(0);
  usp->setComponentNum(1);
  blade1->setComponentNum(2);
  blade2->setComponentNum(3);
  blade3->setComponentNum(4);
  blade4->setComponentNum(5);
  pitch1->setComponentNum(6);
  pitch2->setComponentNum(7);
  pitch3->setComponentNum(8);
  pitch4->setComponentNum(9);
  */

  // Set the elements
  TACSElement *elements[] = { // bodies
                             shaft, 
                             usp, 
                             blade3, pitch3, 
                             blade4, pitch4,
                             // constraints
                             rev_shaft_base, 
                             sph_shaft_usp, 
                             rev_shaft_blade3, rev_blade3_pitch3, 
                             rev_shaft_blade4, rev_blade4_pitch4
  };
  tacs->setElements(elements); 
  
  // Set the connectivity
  int conn[] = {0, 1, 2, 3, 4, 5,
                0,    6,  // shaft_fix
                0, 1, 7,  // shaft and usp
                0, 2, 8,  // blade1 and shaft
                2, 3, 9,  // blade1 and pitch1
                0, 4, 10, // blade2 and shaft
                4, 5, 11  // blade2 and pitch2
               };

  int ptr[]  = {
    // bodies
    0, 1, 2, 3, 4, 5,
    // constraints
    6, 
    8, 11, 
    14, 17, 
    20, 23,
    26};
    
  tacs->setElementConnectivity(conn, ptr);
  tacs->initialize();

  //------------------------------------------------------------------//
  //                 Create the TACSIntegrator object                 //
  //------------------------------------------------------------------//

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, RIGID, write_flag);
  f5->incref();

  double tinit            = 0.0;
  double tfinal           = 1.0;
  int    steps_per_second = 300; 
  int    num_stages       = 2;
  int    max_bdf_order    = 1;
  TACSIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal,
                                               steps_per_second, max_bdf_order);
  bdf->incref();
  
  // Set optional parameters
  bdf->setOrderingType(TACSAssembler::NATURAL_ORDER);
  bdf->setRelTol(1.0e-13);
  bdf->setAbsTol(1.0e-14);
  bdf->setMaxNewtonIters(5);
  bdf->setPrintLevel(1);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(1);
  bdf->configureOutput(f5, 1, "output/rotor_%04d.f5");

  // Integrate and write solution to file
  bdf->integrate();
  bdf->writeSolution("solutionDIRK.dat");

  // Decref objects
  gravVec->decref();
  omega->decref();
  zero->decref();

  // shaft 
  rShaft->decref();
  rShaftx->decref();
  rShaftx->decref();
  fShaft->decref();
  vizShaft->decref();
  shaft->decref();

  // upper swash plate
  rUSP->decref();  
  rUSPx->decref();
  rUSPy->decref();
  fUSP->decref();
  vizUSP->decref();
  usp->decref();

  // inner blade 3
  rBlade3->decref();
  rBlade3x->decref();
  rBlade3y->decref();
  fBlade3->decref();
  vizBlade3->decref();
  blade3->decref();

  // pitch link blade 3
  rPitch3->decref();
  rPitch3x->decref();
  rPitch3y->decref();
  fPitch3->decref();
  vizPitch3->decref();
  pitch3->decref();

  // inner blade 4
  rBlade4->decref();
  rBlade4x->decref();
  rBlade4y->decref();
  fBlade4->decref();
  vizBlade4->decref();
  blade4->decref();

  // pitch link blade 4
  rPitch4->decref();
  rPitch4x->decref();
  rPitch4y->decref();
  fPitch4->decref();
  vizPitch4->decref();
  pitch4->decref();

  xrev->decref();
  yrev->decref();
  zrev->decref();
  
  // shaft base constraint
  pnt_shaft_base->decref();
  rev_shaft_base->decref();

  // shaft usp constraint
  pnt_shaft_usp->decref();
  pnt_shaft_usp->decref();

  // shaft blade 3 constraint
  pnt_shaft_blade3->decref();
  rev_shaft_blade3->decref();
  // blade 3 and pitch link 3
  pnt_blade3_pitch3->decref();
  rev_blade3_pitch3->decref();
  // pitch3 usp constraint
  pnt_pitch3_usp->decref();
  rev_pitch3_usp->decref();

  // shaft blade 4 constraint
  pnt_shaft_blade4->decref();
  rev_shaft_blade4->decref();
  // blade 4 pitch 4 constraint
  pnt_blade4_pitch4->decref();
  rev_blade4_pitch4->decref();
  // pitch 4 usp constraint
  pnt_pitch4_usp->decref();
  rev_pitch4_usp->decref();

  // tacs
  tacs->decref();

  // f5
  f5->decref();

  // integrator
  bdf->decref();

  MPI_Finalize();
  return (0);
}
