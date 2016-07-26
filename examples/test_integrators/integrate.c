#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "TACSDummyElement.h"
#include "TACSIntegrator.h"

/*
  Function that tests the BDF and DIRK integration schemes within TACS
  on a test function that is implemented in TACSDummyElement class
*/
int main( int argc, char *argv[] ){
  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*-----------------------------------------------------------------*/
  /*                    Create TACS Element                          */
  /*-----------------------------------------------------------------*/

  TACSElement *elem = new TACSDummyElement();
  elem->incref();
  
  /*-----------------------------------------------------------------*/
  /*                    Create TACS Creator                          */
  /*-----------------------------------------------------------------*/

  int vars_per_node = 2;
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, vars_per_node);
  creator->incref();

  creator->setElements(&elem, 1);
  creator->setBoundaryConditions(0, NULL, NULL, NULL);

  // Points into the starting index of the connectivity array and is of
  // length num_elements+1
  const int elem_node_ptr[] = {0,1}; 
  const int elem_node_conn[] = {0};

  // Associate elements with element id numbers
  int elem_id_nums[] = {0}; 
  int num_nodes = 1, num_elements = 1;
  creator->setGlobalConnectivity(num_nodes, num_elements,
				 elem_node_ptr, elem_node_conn, elem_id_nums);
  
  const TacsScalar xpts[] = {0,0,0};
  creator->setNodes(xpts);

  /*-----------------------------------------------------------------*/
  /*                    Create TACS Assembler                        */
  /*-----------------------------------------------------------------*/
  
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();
 
  double tinit = 0.0; double tfinal = 10.0;
  int num_steps_per_sec = 10;

  /*-----------------------------------------------------------------*/
  /*                    Test DIRK Scheme                             */
  /*-----------------------------------------------------------------*/

  int num_stages = 2;
  TACSDIRKIntegrator *dirk = new TACSDIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec,
						num_stages);
  dirk->incref();

  // Set optional parameters
  dirk->setRelTol(1.0e-12);
  dirk->setAbsTol(1.0e-14);
  dirk->setMaxNewtonIters(24);
  dirk->setPrintLevel(0);
  dirk->setJacAssemblyFreq(1);
  dirk->setUseLapack(0);

  // Integrate and write solution to file
  dirk->integrate();
  //  dirk->forward(NULL, NULL, NULL);
  dirk->writeSolution("dirk.dat");

  dirk->decref();
  
  //-----------------------------------------------------------------//
  //                    Test BDF Scheme                             //
  //-----------------------------------------------------------------//
  
  int max_bdf_order = 3;
  TACSBDFIntegrator *bdf = new TACSBDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec,
					      max_bdf_order);
  bdf->incref();

  // Set optional parameters
  bdf->setRelTol(1.0e-12);
  bdf->setAbsTol(1.0e-14);
  bdf->setMaxNewtonIters(24);
  bdf->setPrintLevel(0);
  bdf->setJacAssemblyFreq(1);
  bdf->setUseLapack(0);

  // Integrate and write solution to file
  bdf->integrate();
  // bdf->forward(NULL, NULL, NULL);
  bdf->writeSolution("bdf.dat");

  bdf->decref();

  //-----------------------------------------------------------------//
  //                    Test ABM Scheme                             //
  //-----------------------------------------------------------------//
  
  int max_abm_order = 6;
  TACSABMIntegrator *abm = new TACSABMIntegrator(tacs, tinit, tfinal, num_steps_per_sec,
					      max_abm_order);
  abm->incref();

  // Set optional parameters
  abm->setRelTol(1.0e-12);
  abm->setAbsTol(1.0e-14);
  abm->setMaxNewtonIters(24);
  abm->setPrintLevel(1);
  abm->setJacAssemblyFreq(1);
  abm->setUseLapack(0);

  // Integrate and write solution to file
  abm->integrate();
  abm->writeSolution("abm.dat");

  abm->decref();

  //-----------------------------------------------------------------//
  //                    Test NBG Scheme                             //
  //-----------------------------------------------------------------//
  
  TACSNBGIntegrator *nbg = new TACSNBGIntegrator(tacs, tinit, tfinal, num_steps_per_sec);
  nbg->incref();

  // Set optional parameters
  nbg->setRelTol(1.0e-12);
  nbg->setAbsTol(1.0e-14);
  nbg->setMaxNewtonIters(24);
  nbg->setPrintLevel(1);
  nbg->setJacAssemblyFreq(1);
  nbg->setUseLapack(0);

  // Integrate and write solution to file
  nbg->integrate();
  nbg->writeSolution("nbg.dat");

  nbg->decref();

  // Deallocate objects
  tacs->decref();
  elem->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
