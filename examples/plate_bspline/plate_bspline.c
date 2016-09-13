#include "PlaneStressQuad.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"

int main ( int argc, char *argv[] ){
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

  // Only set the mesh/boundary conditions etc. on the root processor
  if (rank == 0){
    int num_nodes = 36;
    int num_elems = 52;
    int dep_nodes = 33;
    int nx = 5, ny = 5;
    
    int conn[] = {0, 1, 6, 7,
                  1, 2, 7, 8,
                  2, 3, 8, 9,
                  3, 4, 9, 10,
                  4, 5, 10, 11,
                  6, 7, 12, 13,
                  7, 8, 13, 14, 
                  8, 9, 14, 15,
                  9, 10, 15, 16,
                  10, 11, 16, 17,
                  12, 13, 18, 19,
                  13, 14, 19, 20,
                  18, 19, 24, 25,
                  19, 20, 25, 26,
                  24, 25, 30, 31,
                  25, 26, 31, 32,
                  14, -1, -4, -5,
                  -1, 15, -5, -6,
                  15, -2, -6, -7,
                  -2, 16, -7, -8,
                  16, -3, -8, -9,
                  -3, 17, -9, -10,
                  -4, -5, 20, -11,
                  -5, -6, -11, 21,
                  -6, -7, 21, -12,
                  -7, -8, -12, 22,
                  -8, -9, 22, -13, 
                  -9, -10, -13, 23,
                  20, -11, -14, -15,
                  -11, 21, -15, -16,
                  21, -12, -16, -17,
                  -12, 22, -17, -18,
                  22, -13, -18, -19,
                  -13, 23, -19, -20,
                  -14, -15, 26, -21,
                  -15, -16, -21, 27, 
                  -16, -17, 27, -22,
                  -17, -18, -22, 28,
                  -18, -19, 28, -23,
                  -19, -20, -23, 29,
                  26, -21, -24, -25,
                  -21, 27, -25, -26,
                  27, -22, -26, -27,
                  -22, 28, -27, -28,
                  28, -23, -28, -29,
                  -23, 29, -29, -30,
                  -24, -25, 32, -31,
                  -25, -26, -31, 33,
                  -26, -27, 33, -32,
                  -27, -28, -32, 34,
                  -28, -29, 34, -33,
                  -29, -30, -33, 35};

    int *elem_node_ptr = new int[num_elems+1];
    elem_node_ptr[0] = 0;
    for (int i = 0; i < num_elems; i++){
      elem_node_ptr[i+1] = elem_node_ptr[i]+4;
    }
    int *dep_ptr = new int[dep_nodes+1];
    int *dep_conn = new int[4*dep_nodes];
    double *dep_weights = new double[4*dep_nodes];
    int d_node = 14;
    for (int jdex = 0; jdex < dep_nodes; jdex+=10){
      // 1st dependent node
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;

      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+6;
      dep_conn[4*jdex+3] = d_node+7;
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node, d_node+1, d_node+6, d_node+7); */
      jdex++;
      d_node++;
      // 2nd dependent node
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;

      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+6;
      dep_conn[4*jdex+3] = d_node+7;
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node, d_node+1, d_node+6, d_node+7); */
      jdex++;
      d_node++;
      // 3rd dependent node
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;

      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+6;
      dep_conn[4*jdex+3] = d_node+7;
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node, d_node+1, d_node+6, d_node+7); */
      jdex-=2;
      d_node+=4;
    }
    d_node = 26;
    for (int jdex = 30; jdex < dep_nodes; jdex++){
     
      dep_weights[4*jdex] = 0.0;
      dep_weights[4*jdex+1] = 0.0;
      dep_weights[4*jdex+2] = 0.5;
      dep_weights[4*jdex+3] = 0.5;

      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+6;
      dep_conn[4*jdex+3] = d_node+7;
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node, d_node+1, d_node+6, d_node+7); */
      d_node++;
    }
    // Reset to initial node count
    d_node = 14;
    // Assign weights to dependent nodes on 2 independent nodes
    for (int jdex = 3; jdex < dep_nodes-3; jdex+=2){
      if ((jdex-1) % 10 == 0){
        jdex += 2;
        d_node+=2;
      }      
      // Set the dependent weights
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+6;
      dep_conn[4*jdex+2] = d_node+1;
      dep_conn[4*jdex+3] = d_node+7;
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node, d_node+6, d_node+1, d_node+7); */
      d_node++;
    }
    dep_weights[4*29] = 0.0;
    dep_weights[4*29+1] = 0.0;
    dep_weights[4*29+2] = 0.5;
    dep_weights[4*29+3] = 0.5;
    
    dep_conn[4*29] = 28;
    dep_conn[4*29+1] = 34;
    dep_conn[4*29+2] = 29;
    dep_conn[4*29+3] = 35;    
    
    // Reset to initial node count
    d_node = 14;
    // Assign weights to dependents nodes on 4 independent nodes
    for (int jdex = 4; jdex < dep_nodes-4; jdex+=2){
      if (jdex % 10 == 0){
        jdex += 4;
        d_node+=3;
      }
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*         d_node,d_node+1, d_node+6, d_node+7); */
      // Set the dependent weights
      dep_weights[4*jdex] = 0.25;
      dep_weights[4*jdex+1] = 0.25;
      dep_weights[4*jdex+2] = 0.25;
      dep_weights[4*jdex+3] = 0.25;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+6;
      dep_conn[4*jdex+3] = d_node+7;
      d_node++;
    }
    dep_ptr[0] = 0;
    for (int i = 0; i < dep_nodes; i++){
      dep_ptr[i+1] = dep_ptr[i] + 4;
    }
    
    // Set the identity numbers
    int *elem_id_nums = new int[num_elems];
    memset(elem_id_nums, 0, num_elems*sizeof(int));
    
    // Set the boundary conditions
    int num_bcs = 6;
    int bc_nodes[] = {0, 6, 12, 18, 24, 30};

    // Set the length along x and y direction
    double Lx = 5.0, Ly = 5.0;

    // Create the nodes
    TacsScalar *Xpts = new TacsScalar[3*num_nodes];
        
    int node = 0;
    // Loop over the nodes
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){        
        Xpts[3*node] = Lx*i/nx;
        Xpts[3*node+1] = Ly*j/ny;
        Xpts[3*node+2] = 0.0;
        node++;
      }
    }
   
    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elems,
                                   elem_node_ptr, conn,
                                   elem_id_nums);
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes);

    // Set the nodal locations
    creator->setNodes(Xpts);
        
    // Set the dependent nodes
    creator->setDependentNodes(dep_nodes, dep_ptr,
                               dep_conn, dep_weights);

    // Free all the allocated data
    delete [] elem_id_nums;
    delete [] elem_node_ptr;
    delete [] dep_ptr;
    delete [] dep_weights;
    delete [] dep_conn;
    delete [] Xpts;
  }
  // This call must occur on all processor
  creator->setElements(&elem, 1);
  
  // Set the reordering typr
  creator->setReorderingType(TACSAssembler::TACS_AMD_ORDER,
                             TACSAssembler::APPROXIMATE_SCHUR);

  // Create the TACSAssembler object
  TACSAssembler *tacs = creator->createTACS();
  tacs->incref();

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
  
  // Create TACSToFH5 object for writing to tecplot
  unsigned int write_flag = (TACSElement::OUTPUT_NODES |
                             TACSElement::OUTPUT_DISPLACEMENTS |
                             TACSElement::OUTPUT_STRAINS |
                             TACSElement::OUTPUT_STRESSES |
                             TACSElement::OUTPUT_EXTRAS);
  TACSToFH5 *f5 = new TACSToFH5(tacs, PLANE_STRESS, write_flag);
  f5->incref();
  f5->writeToFile("plate.f5");
  
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
  elem->decref();
  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
