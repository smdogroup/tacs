#include "PlaneStressQuad.h"
#include "PlaneStressBspline.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"

int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);

  // Get the rank of the processor
  int rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Allocate the TACS creator
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, 2);
  creator->incref();

  // Create the stiffness object
  PlaneStressStiffness *stiff = new PlaneStressStiffness(2700.0,
                                                         70.0e9, 0.3);
  stiff->incref();
  
  // Create the plane stress element
  TACSElement *elem = new PlaneStressQuad<2>(stiff);
  elem->incref();
  /* double Tu[] = {0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0}; */
  /* double Tv[] = {0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0}; */
  /* // Number of knot intervals without the repeating knots */
  /* int Lu = 3, Lv = 3; */
  /* TACSElement *elem[9]; */
  /* for (int i = 0; i < 9; i++){ */
  /*   elem[i] = new PlaneStressBspline(stiff,Tu,Tv, */
  /*                                    Lu, Lv, NULL, LINEAR, */
  /*                                    i/3, i); */
  /*   elem[i]->incref(); */
  /* }   */
  
  // Creating mesh on root processor
  if (rank == 0){
    int num_nodes = 16;
    int num_elems = 9;
    int ndep_nodes = 0;

    int conn[] = {0,1,4,5,
                  1,2,5,6,
                  2,3,6,7,
                  4,5,8,9,
                  5,6,9,10,
                  6,7,10,11,
                  8,9,12,13,
                  9,10,13,14,
                  10,11,14,15};
    int ptr[] = {0,4,8,12,16,20,24,28,32,36};
    int elem_ids[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    //int elem_ids[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    int dep_conn[] = {4, 5,
                      4, 7,
                      5, 8,
                      7, 8,
                      4, 5, 7, 8};
    int dep_ptr[] = {0, 2, 4, 6, 8, 12};
    double dep_weights[] = {0.5, 0.5, 0.5, 0.5,
                            0.5, 0.5, 0.5, 0.5,
                            0.25, 0.25, 0.25, 0.25};    
    int num_bcs = 4;
    int bc_nodes[] = {0, 4, 8, 12};
    TacsScalar Xpts[] = {0.0, 0.0, 0.0,
                         1.0, 0.0, 0.0,
                         2.0, 0.0, 0.0,
                         3.0, 0.0, 0.0,
                         0.0, 1.0, 0.0,
                         1.0, 1.0, 0.0,
                         2.0, 1.0, 0.0,
                         3.0, 1.0, 0.0,
                         0.0, 2.0, 0.0,
                         1.0, 2.0, 0.0,
                         2.0, 2.0, 0.0,
                         3.0, 2.0, 0.0,
                         0.0, 3.0, 0.0,
                         1.0, 3.0, 0.0,
                         2.0, 3.0, 0.0,
                         3.0, 3.0, 0.0};
    
    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elems,
                                   ptr, conn,
                                   elem_ids);
    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes);

    // Set the nodal locations
    creator->setNodes(Xpts);

    // Set the dependent nodes
    /* creator->setDependentNodes(ndep_nodes, dep_ptr, */
    /*                            dep_conn, dep_weights); */

    /* // Free all the allocated data */
    /* delete [] Xpts; */
    /* delete [] elem_node_ptr; */
    /* delete [] elem_node_conn; */
    /* delete [] elem_id_nums; */
    /* delete [] bc_nodes; */
    /* delete [] anodes; */
    /* delete [] dep_ptr; */
    /* delete [] dep_conn; */
    /* delete [] dep_weights; */
  }
  
  // This call must occur on all processor
  creator->setElements(&elem, 1);
  
  // Set the reordering typr
  creator->setReorderingType(TACSAssembler::TACS_AMD_ORDER,
                             TACSAssembler::APPROXIMATE_SCHUR);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();
  /* tacs->testElement(0,2); */
  /* tacs->decref(); */
  /* exit(0); */
  // Create the preconditioner
  TACSBVec *res = tacs->createVec();
  TACSBVec *ans = tacs->createVec();
  FEMat *mat = tacs->createFEMat();
  
  // Increment the reference count to the matrix/vectors
  res->incref();
  ans->incref();
  mat->incref();

  // Allocate the factorization
  int lev = 4500;
  double fill = 10.0;
  int reorder_schur = 1;
  PcScMat *pc = new PcScMat(mat, lev, fill, reorder_schur);
  pc->incref();

  // Assemble and factor the stiffness/Jacobian matrix
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  tacs->assembleJacobian(alpha, beta, gamma, res, mat);
  mat->applyBCs();
  pc->factor();
  
  // Number of GMRES iterations
  int gmres_iters = 10;
  // Number of allowed restartions
  int nrestart = 2;
  // Is it a flexible preconditioner
  int is_flexible = 1;
  GMRES *gmres = new GMRES(mat, pc, gmres_iters,
                           nrestart, is_flexible);
  gmres->incref();
  
  res->set(1.0);
  res->applyBCs();
  gmres->solve(res, ans);
  
  tacs->setVariables(ans);
  TacsScalar *ans_array;
  int ans_size = ans->getArray(&ans_array);
 
  // Create TACSToFH5 object for writing to tecplot
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
  f5->writeToFile("plate_simple.f5");
  
  // Free everything
  f5->decref();
  
  // Decrease the reference count to the linear algebra objects
  gmres->decref();
  pc->decref();
  mat->decref();
  ans->decref();
  res->decref();

  // Decrease the reference count to everything else
  stiff->decref();
  /* if (elem){ */
  /*   for (int i = 0; i < 9; i++){ */
  /*     elem[i]->decref(); */
  /*   } */
  /* } */
  elem->decref();
  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
