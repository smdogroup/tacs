#include "PlaneStressBspline.h"
#include "PlaneStressBsplineAll.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"
// Convert patch mesh to control point mesh
void controlPointMesh(int *num_nodes_x, int *num_nodes_y,
                      int *num_nodes,
                      int *num_patches_x, int *num_patches_y,
                      int **ptr, int **conn, int order){
  *num_nodes_x = *num_patches_x+(order-1);
  *num_nodes_y = *num_patches_y+(order-1);
  *num_nodes = (*num_nodes_x)*(*num_nodes_y);
 
  int num_patches = (*num_patches_x)*(*num_patches_y);
  int *_ptr = new int[num_patches+1];
  int *_conn = new int[num_patches*order*order];
  for (int i = 0; i < num_patches+1; i++){
    _ptr[i] = i*order*order;
  }
  int node = 0;
  int k = 0;
  if (order == 3){
    for (int j = 0; j < *num_patches_y; j++){
      for (int i = 0; i < *num_patches_x; i++){
        if (i == 0){
          node = j*(*num_nodes_x);
        }
        _conn[9*k] = node;
        _conn[9*k+1] = node+1;
        _conn[9*k+2] = node+2;
      
        _conn[9*k+3] = node+*num_nodes_x;
        _conn[9*k+4] = node+*num_nodes_x+1;
        _conn[9*k+5] = node+*num_nodes_x+2;
      
        _conn[9*k+6] = node+*num_nodes_x*2;
        _conn[9*k+7] = node+*num_nodes_x*2+1;
        _conn[9*k+8] = node+*num_nodes_x*2+2;
        node++;
        k++;
      }
    }
  }
  else if (order == 4){
    for (int j = 0; j < *num_patches_y; j++){
      for (int i = 0; i < *num_patches_x; i++){
        if (i == 0){
          node = j*(*num_nodes_x);
        }
        _conn[16*k] = node;
        _conn[16*k+1] = node+1;
        _conn[16*k+2] = node+2;
        _conn[16*k+3] = node+3;
        
        _conn[16*k+4] = node+*num_nodes_x;
        _conn[16*k+5] = node+*num_nodes_x+1;
        _conn[16*k+6] = node+*num_nodes_x+2;
        _conn[16*k+7] = node+*num_nodes_x+3;
        
        _conn[16*k+8] = node+*num_nodes_x*2;
        _conn[16*k+9] = node+*num_nodes_x*2+1;
        _conn[16*k+10] = node+*num_nodes_x*2+2;
        _conn[16*k+11] = node+*num_nodes_x*2+3;
        
        _conn[16*k+12] = node+*num_nodes_x*3;
        _conn[16*k+13] = node+*num_nodes_x*3+1;
        _conn[16*k+14] = node+*num_nodes_x*3+2;
        _conn[16*k+15] = node+*num_nodes_x*3+3;
        node++;
        k++;
      }
    }
  }
  *ptr = _ptr;
  *conn = _conn;
  //delete [] _ptr;
  //delete [] _conn;
}
int main( int argc, char *argv[] ){
  MPI_Init(&argc, &argv);
  
  // Get the rank of the processor
  int rank,size; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Allocate the TACS creator
  TACSCreator *creator = new TACSCreator(MPI_COMM_WORLD, 2);
  creator->incref();

  // Create the stiffness object
  PlaneStressStiffness *stiff = new PlaneStressStiffness(2700.0,
                                                         70.0e9, 0.3);
  stiff->incref();
  
  // Number of knot intervals without the repeating knots
  int Lu = 4, Lv = 4;
  TACSElement *elem[Lu*Lv];
  int ind = 0;
  int order = 3;
  // Order of the Bspline
  for (int k = 0; k < argc; k++){
    int _order;
    if (sscanf(argv[k], "order=%d", &_order) == 1){
      order = _order;
    }
  }
  printf("Order is %d \n", order);
  if (order == 4){
    for (int j = 0; j < Lv; j++){
      for (int i = 0; i < Lu; i++){
        double Tu[] = {0.0,0.0,0.0,0.0,1.0,2.0,3.0,4.0,4.0,4.0,4.0};
        double Tv[] = {0.0,0.0,0.0,0.0,1.0,2.0,3.0,4.0,4.0,4.0,4.0};
        
        elem[ind] = new PlaneStressBsplineAll<4>(stiff,Tu,Tv,
                                                 Lu, Lv,
                                                 LINEAR,
                                                 0, ind);
        
        elem[ind]->incref();
        ind++;
      }
    }
  }
  else if (order == 3){
    for (int j = 0; j < Lv; j++){
      for (int i = 0; i < Lu; i++){
        double Tu[] = {0.0,0.0,0.0,1.0,2.0,3.0,4.0,4.0,4.0};
        double Tv[] = {0.0,0.0,0.0,1.0,2.0,3.0,4.0,4.0,4.0};
        elem[ind] = new PlaneStressBsplineAll<3>(stiff,Tu,Tv,
                                                 Lu, Lv,
                                                 LINEAR,
                                                 0, ind); 
       
        /* elem[ind] = new PlaneStressBspline(stiff,Tu,Tv, */
        /*                                    Lu, Lv, */
        /*                                    LINEAR, */
        /*                                    0, ind); */
        elem[ind]->incref();
        ind++;
      }
    }
  }
  // Create mesh on root processor
  if (rank == 0){
    int num_nodes, num_nodes_x, num_nodes_y;
    int num_elems_x = 1*Lu, num_elems_y = 1*Lv;
    int num_elems = num_elems_x*num_elems_y;
    int *ptr_c, *conn_c;
    int elem_ids[Lu*Lv];
    for (int i = 0; i < Lu*Lv; i++){
      elem_ids[i] = i*1;
    }
    controlPointMesh(&num_nodes_x, &num_nodes_y,
                     &num_nodes, &num_elems_x,
                     &num_elems_y,
                     &ptr_c, &conn_c, order);
    /* for (int i = 0; i < Lu*Lv; i++){ */
    /*   for (int k = 0; k < 16; k++){ */
    /*     printf("%d ", conn_c[16*i+k]); */
    /*   } */
    /*   printf("\n"); */
    /* } */

    creator->setGlobalConnectivity(num_nodes, num_elems, 
                                   ptr_c, conn_c,elem_ids);
    
    int num_bcs_c = num_nodes_y*1;
    int bc_nodes_c[num_bcs_c];
    for (int i = 0; i < num_bcs_c; i++){
      bc_nodes_c[i] = i*num_nodes_y;
    }    
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs_c, bc_nodes_c);
    if (order == 4){
      TacsScalar Xpts_c[] = {0.0, 0.0, 0.0,
                             0.667, 0.0, 0.0,
                             1.333, 0.0, 0.0,
                             2.0, 0.0, 0.0,
                             2.667, 0.0, 0.0,
                             3.333, 0.0, 0.0,
                             4.0, 0.0, 0.0,
                             0.0, 0.667, 0.0,
                             0.667, 0.667, 0.0,
                             1.333, 0.667, 0.0,
                             2.0, 0.667, 0.0,
                             2.667, 0.667, 0.0,
                             3.333, 0.667, 0.0,
                             4.0, 0.667, 0.0,
                             0.0, 1.333, 0.0,
                             0.667, 1.333, 0.0,
                             1.333, 1.333, 0.0,
                             2.0, 1.333, 0.0,
                             2.667, 1.333, 0.0,
                             3.333, 1.333, 0.0,
                             4.0, 1.333, 0.0,
                             0.0, 2.0, 0.0,
                             0.667, 2.0, 0.0,
                             1.333, 2.0, 0.0,
                             2.0, 2.0, 0.0,
                             2.667, 2.0, 0.0,
                             3.333, 2.0, 0.0,
                             4.0, 2.0, 0.0,
                             0.0, 2.667, 0.0,
                             0.667, 2.667, 0.0,
                             1.333, 2.667, 0.0,
                             2.0, 2.667, 0.0,
                             2.667, 2.667, 0.0,
                             3.333, 2.667, 0.0,
                             4.0, 2.667, 0.0,
                             0.0, 3.333, 0.0,
                             0.667, 3.333, 0.0,
                             1.333, 3.333, 0.0,
                             2.0, 3.333, 0.0,
                             2.667, 3.333, 0.0,
                             3.333, 3.333, 0.0,
                             4.0, 3.333, 0.0,
                             0.0, 4.0, 0.0,
                             0.667, 4.0, 0.0,
                             1.333, 4.0, 0.0,
                             2.0, 4.0, 0.0,
                             2.667, 4.0, 0.0,
                             3.333, 4.0, 0.0,
                             4.0, 4.0, 0.0};
      // Set the nodal locations
      creator->setNodes(Xpts_c);
    }
    if (order == 3){
      TacsScalar Xpts_c[] = {0.0, 0.0, 0.0,
                             0.8, 0.0, 0.0,
                             1.6, 0.0, 0.0,
                             2.4, 0.0, 0.0,
                             3.2, 0.0, 0.0,
                             4.0, 0.0, 0.0,
                             0.0, 0.8, 0.0,
                             0.8, 0.8, 0.0,
                             1.6, 0.8, 0.0,
                             2.4, 0.8, 0.0,
                             3.2, 0.8, 0.0,
                             4.0, 0.8, 0.0,
                             0.0, 1.6, 0.0,
                             0.8, 1.6, 0.0,
                             1.6, 1.6, 0.0,
                             2.4, 1.6, 0.0,
                             3.2, 1.6, 0.0,
                             4.0, 1.6, 0.0,
                             0.0, 2.4, 0.0,
                             0.8, 2.4, 0.0,
                             1.6, 2.4, 0.0,
                             2.4, 2.4, 0.0,
                             3.2, 2.4, 0.0,
                             4.0, 2.4, 0.0,
                             0.0, 3.2, 0.0,
                             0.8, 3.2, 0.0,
                             1.6, 3.2, 0.0,
                             2.4, 3.2, 0.0,
                             3.2, 3.2, 0.0,
                             4.0, 3.2, 0.0,
                             0.0, 4.0, 0.0,
                             0.8, 4.0, 0.0,
                             1.6, 4.0, 0.0,
                             2.4, 4.0, 0.0,
                             3.2, 4.0, 0.0,
                             4.0, 4.0, 0.0};
      // Set the nodal locations
      creator->setNodes(Xpts_c);
    }
    /* TacsScalar Xpts_c[] = {0.0, 0.0, 0.0, */
    /*                        0.75, 0.0, 0.0, */
    /*                        1.5, 0.0, 0.0, */
    /*                        2.25, 0.0, 0.0, */
    /*                        3.0, 0.0, 0.0, */
    /*                        0.0, 0.75, 0.0, */
    /*                        0.75, 0.75, 0.0, */
    /*                        1.5, 0.75, 0.0, */
    /*                        2.25, 0.75, 0.0, */
    /*                        3.0, 0.75, 0.0, */
    /*                        0.0, 1.5, 0.0, */
    /*                        0.75, 1.5, 0.0, */
    /*                        1.5, 1.5, 0.0, */
    /*                        2.25, 1.5, 0.0, */
    /*                        3.0, 1.5, 0.0, */
    /*                        0.0, 2.25, 0.0, */
    /*                        0.75, 2.25, 0.0, */
    /*                        1.5, 2.25, 0.0, */
    /*                        2.25, 2.25, 0.0, */
    /*                        3.0, 2.25, 0.0, */
    /*                        0.0, 3.0, 0.0, */
    /*                        0.75, 3.0, 0.0, */
    /*                        1.5, 3.0, 0.0, */
    /*                        2.25, 3.0, 0.0, */
    /*                        3.0, 3.0, 0.0};     */

  }
  // This call must occur on all processor
  creator->setElements(&elem[0], Lu*Lv);
  // Set the reordering typr
  creator->setReorderingType(TACSAssembler::NATURAL_ORDER,
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
 
  /* TacsScalar *ans_array; */
  /* int ans_size = ans->getArray(&ans_array); */
  /* for (int i = 0; i < ans_size; i++){ */
  /*   printf("%e \n", ans_array[i]); */
  /* } */
  /* printf("Ax = %e \n", ans->norm()); */
  /* exit(0); */

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
  TacsScalar *res_array;
  int res_size = res->getArray(&res_array);
  res_array[res_size-1] = 1.0;
  //res->set(1.0);
  res->applyBCs();
  gmres->solve(res, ans);
  tacs->setVariables(ans);
  /* TacsScalar *ans_array;  */
  /* int ans_size = ans->getArray(&ans_array); */
  /* for (int i = 0; i < ans_size; i++){ */
  /*   printf("%e \n", ans_array[i]); */
  /* } */
  // Create TACSToFH5 object for writing to tecplot
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
  char filename[256];
  sprintf(filename, "plate_test%d.f5",order);
  f5->writeToFile(filename);
  
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
  if (elem){
    for (int i = 0; i < 1; i++){
      elem[i]->decref();
    }
  }
  creator->decref();

  MPI_Finalize();
  return 0;
}
