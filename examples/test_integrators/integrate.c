#include "TACSCreator.h"
#include "TACSAssembler.h"
#include "TACSDummyElement.h"
#include "TACSIntegrator.h"
#include "KSFailure.h"
/*
  Function that tests the BDF and DIRK integration schemes within TACS
  on a test function that is implemented in TACSDummyElement class
*/
int main( int argc, char *argv[] ){
  
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
 
  /*-----------------------------------------------------------------*/
  /*                    Test DIRK Scheme                             */
  /*-----------------------------------------------------------------*/

  printf(">> Testing DIRK\n");  
 
  double tinit = 0.0;
  double tfinal = 25.0;

  int num_steps_per_sec = 10;
  int num_stages = 1;
  
  TacsIntegrator *dirk = new TacsDIRKIntegrator(tacs, tinit, tfinal, num_steps_per_sec,
						num_stages);
  dirk->incref();

  dirk->integrate();
  dirk->writeSolution("dirk.dat");

  dirk->decref();
  
  printf(">> Testing DIRK complete\n");  

  //-----------------------------------------------------------------//
  //                    Test BDF Scheme                             //
  //-----------------------------------------------------------------//
  
  printf(">> Testing BDF\n");  

  int max_bdf_order = 3;
  TacsIntegrator *bdf = new TacsBDFIntegrator(tacs, tinit, tfinal, num_steps_per_sec,
					      max_bdf_order);
  bdf->incref();

  bdf->integrate();
  bdf->writeSolution("bdf.dat");

  // Create an objective function
  TACSFunction *func =  new KSFailure(tacs, elem_id_nums, num_elements, 1.0, 100.0);
  func->incref();

  bdf->setFunction(&func, 1);
  bdf->adjointSolve();

  bdf->decref();

  printf(">> BDF complete \n");

  // Deallocate objects
  func->decref();

  tacs->decref();
  elem->decref();
  creator->decref();

  MPI_Finalize();

  return (0);
}
