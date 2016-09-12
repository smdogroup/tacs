#include "PlaneStressQuad.h"
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

  // Creating mesh on root processor
  if (rank == 0){
    // Create a mesh of quad elements
    int nx = 2, ny = 2;
    int num_elems = nx*ny;
    
    // Allocate and number the nodes
    int lena = (nx+1)*(ny+1)+5;
    int *anodes = new int[lena];
    memset(anodes, 0, lena*sizeof(int));

    int dep_nodes = 0;
    // Set dependent nodes and weights
    for (int i = (nx+1)*(ny+1); i < lena; i++){
      anodes[i] = -dep_nodes-1;
      dep_nodes++;
    }
    // Number the nodes
    int num_nodes = 0;
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){
        int index = i+(nx+1)*j;
        anodes[i+(nx+1)*j] = num_nodes;
        num_nodes++;
      }
    }
    // Allocate dependent node data structure
    int *dep_ptr = new int[dep_nodes+1];
    int *dep_conn = new int[4*dep_nodes];
    double *dep_weights = new double[4*dep_nodes];

    dep_ptr[0] = 0;
    for (int i = 0; i < dep_nodes; i++){
      dep_ptr[i+1] = dep_ptr[i]+4;
    }
    int d_node = 4;
    // For 1st dependent node
    dep_weights[0] = 0.5;
    dep_weights[1] = 0.5;
    dep_weights[2] = 0.0;
    dep_weights[3] = 0.0;
    
    dep_conn[0] = 4;
    dep_conn[1] = 5;
    dep_conn[2] = 7; 
    dep_conn[3] = 8;
    
    // For 2nd dependent node
    dep_weights[1*4] = 0.5;
    dep_weights[1*4+1] = 0.0;
    dep_weights[1*4+2] = 0.5;
    dep_weights[1*4+3] = 0.0;
    
    dep_conn[1*4] = 4;
    dep_conn[1*4+1] = 5;
    dep_conn[1*4+2] = 7; 
    dep_conn[1*4+3] = 8;

    // For 3rd dependent node
    dep_weights[2*4] = 0.25;
    dep_weights[2*4+1] = 0.25;
    dep_weights[2*4+2] = 0.25;
    dep_weights[2*4+3] = 0.25;
    
    dep_conn[2*4] = 4;
    dep_conn[2*4+1] = 5;
    dep_conn[2*4+2] = 7; 
    dep_conn[2*4+3] = 8;

    // For 4th dependent node
    dep_weights[3*4] = 0.0;
    dep_weights[3*4+1] = 0.5;
    dep_weights[3*4+2] = 0.0;
    dep_weights[3*4+3] = 0.5;
    
    dep_conn[3*4] = 4;
    dep_conn[3*4+1] = 5;
    dep_conn[3*4+2] = 7; 
    dep_conn[3*4+3] = 8;
    
    // For 5th dependent node
    dep_weights[4*4] = 0.0;
    dep_weights[4*4+1] = 0.0;
    dep_weights[4*4+2] = 0.5;
    dep_weights[4*4+3] = 0.5;
    
    dep_conn[4*4] = 4;
    dep_conn[4*4+1] = 5;
    dep_conn[4*4+2] = 7; 
    dep_conn[4*4+3] = 8;
      
    // Set the length along the x and y direction
    double Lx = 2.0, Ly = 2.0;
    
    // Create the nodes
    TacsScalar *Xpts = new TacsScalar[3*num_nodes];
    
    // Loop over the nodes
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){
        int node = anodes[i+(nx+1)*j];
        if (node >= 0){
          Xpts[3*node] = Lx*i/nx;
          Xpts[3*node+1] = Ly*j/ny;
          Xpts[3*node+2] = 0.0;
        }
      }
    }
    
    // Set up the element connectivity array
    int *elem_node_conn = new int[4*num_elems];
    int *elem_node_ptr = new int[num_elems+1];
    // Element counter
    int n = 0;
    int *conn = elem_node_conn;
    elem_node_ptr[0] = 0;
    // Add the element for the mesh
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        conn[0] = anodes[i+(nx+1)*j];
        conn[1] = anodes[i+1+(nx+1)*j];
        conn[2] = anodes[i+(nx+1)*(j+1)];
        conn[3] = anodes[i+1+(nx+1)*(j+1)];
        conn += 4;
        elem_node_ptr[n+1] = elem_node_ptr[n]+4;
        n++;
      }
    }
    // Set the identity numbers
    int *elem_id_nums = new int[num_elems];
    memset(elem_id_nums, 0, num_elems*sizeof(int));
    
    // Set the boundary conditions
    int num_bcs = ny+1;
    int *bc_nodes = new int[num_bcs];
    for (int j = 0; j < ny+1; j++){
      bc_nodes[j] = anodes[j*(nx+1)];
    }
    // Set the connectivity
    creator->setGlobalConnectivity(num_nodes, num_elems,
                                   elem_node_ptr, elem_node_conn,
                                   elem_id_nums);
    // Set the boundary conditions
    creator->setBoundaryConditions(num_bcs, bc_nodes);

    // Set the nodal locations
    creator->setNodes(Xpts);

    // Set the dependent nodes
    creator->setDependentNodes(dep_nodes, dep_ptr,
                               dep_conn, dep_weights);

    // Free all the allocated data
    delete [] Xpts;
    delete [] elem_node_ptr;
    delete [] elem_node_conn;
    delete [] elem_id_nums;
    delete [] bc_nodes;
    delete [] anodes;
    delete [] dep_ptr;
    delete [] dep_conn;
    delete [] dep_weights;
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
  /* // Set node numbering as displacement */
  /* TacsScalar *res_array; */
  /* int res_size = ans->getArray(&res_array); */
  /* for (int i = 0; i < res_size; i++){ */
  /*   res_array[i] = 1.0*i; */
  /* } */
  tacs->setVariables(ans);
    
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
  elem->decref();
  tacs->decref();
  creator->decref();

  MPI_Finalize();
  return (0);
}
