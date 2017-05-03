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
  TACSGibbsVector  *rUSPy  = new TACSGibbsVector(0.5, 1.0, 1.5);
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
  //                 Blade 1
  //------------------------------------------------------------------//
  
  // blade is 5m long

  TACSGibbsVector  *rBlade1   = new TACSGibbsVector(2.5 + 0.0, 0.0, 2.0);
  rBlade1->incref();
  TACSGibbsVector  *rBlade1x  = new TACSGibbsVector(2.5 + 1.0, 0.0, 2.0);
  rBlade1x->incref();
  TACSGibbsVector  *rBlade1y  = new TACSGibbsVector(2.5 + 0.0, 1.0, 2.0);
  rBlade1y->incref();
  TACSRefFrame     *fBlade1   = new TACSRefFrame(rBlade1, rBlade1x, rBlade1y);
  fBlade1->incref();
  TACSRigidBodyViz *vizBlade1 = new TACSRigidBodyViz(5.0, 0.2, 0.05);
  vizBlade1->incref();
  TACSRigidBody    *blade1    = new TACSRigidBody(fBlade1,
                                                  mA, cA, JA,
                                                  rBlade1, zero, omega, gravVec);
  blade1->incref();
  blade1->setVisualization(vizBlade1);

  //------------------------------------------------------------------//
  //                 Pitch Link 1
  //------------------------------------------------------------------//

  // pitch link is half a meter long

  TACSGibbsVector  *rPitch1   = new TACSGibbsVector(0.15 + 0.0, pitch_offset + 0.0, 1.75);
  rPitch1->incref();
  TACSGibbsVector  *rPitch1x  = new TACSGibbsVector(0.15 + 1.0, pitch_offset + 0.0, 1.75);
  rPitch1x->incref();
  TACSGibbsVector  *rPitch1y  = new TACSGibbsVector(0.15 + 0.0, pitch_offset + 1.0, 1.75);
  rPitch1y->incref();
  TACSRefFrame     *fPitch1   = new TACSRefFrame(rPitch1, rPitch1x, rPitch1y);
  fPitch1->incref();
  TACSRigidBodyViz *vizPitch1 = new TACSRigidBodyViz(0.02, 0.02, 0.5);
  vizPitch1->incref();
  TACSRigidBody    *pitch1    = new  TACSRigidBody(fPitch1,
                                                   mA, cA, JA,
                                                   rPitch1, zero, omega, gravVec);
  pitch1->incref();
  pitch1->setVisualization(vizPitch1);

  //------------------------------------------------------------------//
  //                 Blade 2 (left hand side)
  //------------------------------------------------------------------//
  
  // blade is 5m long

  TACSGibbsVector  *rBlade2   = new TACSGibbsVector(-2.5, 0.0, 2.0);
  rBlade2->incref();
  TACSGibbsVector  *rBlade2x  = new TACSGibbsVector(-2.5 + 1.0, 0.0, 2.0);
  rBlade2x->incref();
  TACSGibbsVector  *rBlade2y  = new TACSGibbsVector(-2.5, 1.0, 2.0);
  rBlade2y->incref();
  TACSRefFrame     *fBlade2   = new TACSRefFrame(rBlade2, rBlade2x, rBlade2y);
  fBlade2->incref();
  TACSRigidBodyViz *vizBlade2 = new TACSRigidBodyViz(5.0, 0.2, 0.05);
  vizBlade2->incref();
  TACSRigidBody    *blade2    = new TACSRigidBody(fBlade2,
                                                  mA, cA, JA,
                                                  rBlade2, zero, omega, gravVec);
  blade2->incref();
  blade2->setVisualization(vizBlade2);

  //------------------------------------------------------------------//
  //                 Pitch Link 2 (left hand side)
  //------------------------------------------------------------------//

  // pitch link is half a meter long
  
  TACSGibbsVector  *rPitch2   = new TACSGibbsVector(-0.15 + 0.0, -pitch_offset + 0.0, 1.75);
  rPitch2->incref();
  TACSGibbsVector  *rPitch2x  = new TACSGibbsVector(-0.15 + 1.0, -pitch_offset + 0.0, 1.75);
  rPitch2x->incref();
  TACSGibbsVector  *rPitch2y  = new TACSGibbsVector(-0.15 + 0.0, -pitch_offset + 1.0, 1.75);
  rPitch2y->incref();
  TACSRefFrame     *fPitch2   = new TACSRefFrame(rPitch2, rPitch2x, rPitch2y);
  fPitch2->incref();
  TACSRigidBodyViz *vizPitch2 = new TACSRigidBodyViz(0.02, 0.02, 0.5);
  vizPitch2->incref();
  TACSRigidBody    *pitch2    = new  TACSRigidBody(fPitch2,
                                                   mA, cA, JA,
                                                   rPitch2, zero, omega, gravVec);
  pitch2->incref();
  pitch2->setVisualization(vizPitch2);

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

  // right blade 1 (along +x direction)
  TACSGibbsVector         *pnt_shaft_blade1  = new TACSGibbsVector(0.0, 0.0, 2.0);
  pnt_shaft_blade1->incref();
  TACSRevoluteConstraint *rev_shaft_blade1  = new TACSRevoluteConstraint(shaft, blade1, pnt_shaft_blade1, xrev);
  rev_shaft_blade1->incref();
  TACSGibbsVector         *pnt_blade1_pitch1 = new TACSGibbsVector(0.15, pitch_offset, 2.0);
  pnt_blade1_pitch1->incref();
  TACSSphericalConstraint *rev_blade1_pitch1 = new TACSSphericalConstraint(blade1, pitch1, pnt_blade1_pitch1);
  rev_blade1_pitch1->incref();
  TACSGibbsVector         *pnt_pitch1_usp    = new TACSGibbsVector(0.15, pitch_offset, 1.5);
  pnt_pitch1_usp->incref();
  TACSSphericalConstraint *rev_pitch1_usp    = new TACSSphericalConstraint(pitch1, usp, pnt_pitch1_usp);
  rev_pitch1_usp->incref();

  // left blade 2 (along -x direction)
  TACSGibbsVector         *pnt_shaft_blade2  = new TACSGibbsVector(0.0, 0.0, 2.0);
  pnt_shaft_blade2->incref();
  TACSRevoluteConstraint *rev_shaft_blade2  = new TACSRevoluteConstraint(shaft, blade2, pnt_shaft_blade2, xrev);
  rev_shaft_blade2->incref();
  TACSGibbsVector         *pnt_blade2_pitch2 = new TACSGibbsVector(-0.15, -pitch_offset, 2.0);
  pnt_blade2_pitch2->incref();
  TACSSphericalConstraint *rev_blade2_pitch2 = new TACSSphericalConstraint(blade2, pitch2, pnt_blade2_pitch2);
  rev_blade2_pitch2->incref();
  TACSGibbsVector         *pnt_pitch2_usp    = new TACSGibbsVector(-0.15, -pitch_offset, 1.5);
  pnt_pitch2_usp->incref();
  TACSSphericalConstraint *rev_pitch2_usp    = new TACSSphericalConstraint(pitch2, usp, pnt_pitch2_usp);
  rev_pitch2_usp->incref();

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

  int num_nodes     = 20;  // Number of finite element nodes
  int vars_per_node = 8;   // 
  int num_elems     = 20;

  TACSAssembler *tacs = new TACSAssembler(MPI_COMM_WORLD, vars_per_node,
                                          num_nodes, num_elems);
  tacs->incref();

  // Set all the component numbers for visualization
  /*
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
    blade1, pitch1, 
    blade2, pitch2, 
    blade3, pitch3, 
    blade4, pitch4,
    // constraints
    rev_shaft_base, 
    sph_shaft_usp, 
    rev_shaft_blade1, rev_blade1_pitch1, 
    rev_shaft_blade2, rev_blade2_pitch2, 
    rev_shaft_blade3, rev_blade3_pitch3, 
    rev_shaft_blade4, rev_blade4_pitch4 
  };
  tacs->setElements(elements); 
  
  // Set the connectivity
  int conn[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
                0,    10, // shaft_fix
                0, 1, 11, // shaft and usp

                0, 2, 12, // blade1 and shaft
                2, 3, 13, // blade1 and pitch1

                0, 4, 14, // blade2 and shaft
                4, 5, 15, // blade2 and pitch2

                0, 6, 16, // blade3 and shaft
                6, 7, 17, // blade3 and pitch3

                0, 8, 18, // blade4 and shaft
                8, 9, 19  // blade4 and pitch4
  };

  int ptr[]  = {
    // bodies
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 
    // constraints
    10,                         // base
    12,                         // shaft with swash plate
    15, 18, 
    21, 24, 
    27, 30, 
    33, 36, 
    39};
    
  tacs->setElementConnectivity(conn, ptr);
  tacs->initialize();

  //------------------------------------------------------------------//
  //                 Create the TACSIntegrator object                 //
  //------------------------------------------------------------------//

  // Create an TACSToFH5 object for writing output to files
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, TACS_RIGID, write_flag);
  f5->incref();

  double tinit            = 0.0;
  double tfinal           = 1.0;
  double steps_per_second = 300.0;
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
  bdf->setOutputFrequency(1);

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

  // right blade
  rBlade1->decref();
  rBlade1x->decref();
  rBlade1y->decref();
  fBlade1->decref();
  vizBlade1->decref();  
  blade1->decref();

  // pitch link blade 1
  rPitch1->decref();
  rPitch1x->decref();
  rPitch1y->decref();
  fPitch1->decref();
  vizPitch1->decref();
  pitch1->decref();

  // left blade 2
  rBlade2->decref();
  rBlade2x->decref();
  rBlade2y->decref();
  fBlade2->decref();
  vizBlade2->decref();
  blade2->decref();

  // pitch link blade 2
  rPitch2->decref();
  rPitch2x->decref();
  rPitch2y->decref();
  fPitch2->decref();
  vizPitch2->decref();
  pitch2->decref();

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
  
  // shaft blade1 constraint
  pnt_shaft_blade1->decref();
  rev_shaft_blade1->decref();
  // blade 1 pitch link constraint
  pnt_blade1_pitch1->decref();
  rev_blade1_pitch1->decref();
  // pitch 1 and usp constraint
  pnt_pitch1_usp->decref();
  rev_pitch1_usp->decref();

  // shaft blade 2 constraint
  pnt_shaft_blade2->decref();
  rev_shaft_blade2->decref();
  // blade 2 pitch 2 constraint
  pnt_blade2_pitch2->decref();
  rev_blade2_pitch2->decref();  
  // pitch2 usp constraint
  pnt_pitch2_usp->decref();
  rev_pitch2_usp->decref();

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
