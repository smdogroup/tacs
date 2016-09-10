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
    // Create a mesh of quad elements
    int nx = 5, ny = 5;
    int num_elems = 2*nx*ny;
    int nx_ref = 6, ny_ref = 6;
    
    // Allocate and number the nodes
    int lena = (nx+1)*(ny+1);
    int lenb = (nx+1)*(ny+1)+(nx_ref+1)*(ny_ref+1)-(nx_ref/2+1)*(ny_ref/2+1);    
    int *anodes = new int[lena];
    int *bnodes = new int[lenb];
    
    memset(anodes, 0, lena*sizeof(int));
    memset(bnodes, 0, lenb*sizeof(int));
    
    int dep_nodes = 0;
    // Set dependent nodes and weights
    for (int i = (nx+1)*(ny+1); i < lenb; i++){
      bnodes[i] = -dep_nodes-1;
      dep_nodes++;
    }
    // Number the nodes for the a-block
    int num_nodes = 0;
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){
        int index = i+(nx+1)*j;
        anodes[i+(nx+1)*j] = num_nodes;
        num_nodes++;
      }
    }
    // Copy over the overlapping node indices
    for (int j = 0; j < ny+1; j++){
      bnodes[(nx+1)*j] = anodes[(nx+1)*(j+1)-1];
    }
    // Number all the remaining indices
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){
        int index = i+(nx+1)*j;
        if (bnodes[index] == 0){
          bnodes[index] = num_nodes;
          num_nodes++;
        }
      }
    }
    // Allocate dependent node data structure
    int *dep_ptr = new int[dep_nodes+1];
    int *dep_conn = new int[4*dep_nodes];
    double *dep_weights = new double[4*dep_nodes];    
    
    dep_ptr[0] = 0;
    int d_node = 41;
    for (int jdex = 0; jdex < dep_nodes; jdex+= 10){
      // Set the dependent weights (1st dependent node)
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+nx;
      dep_conn[4*jdex+3] = d_node+1+nx;
      jdex++;
      d_node++;
      
      // Set the dependent weights (2nd dependent node)
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+nx;
      dep_conn[4*jdex+3] = d_node+1+nx;
      jdex++;
      d_node++;
      
      // Set the dependent weights (3rd dependent node)
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+nx;
      dep_conn[4*jdex+3] = d_node+1+nx;
      jdex-=2;
      d_node+=3;
    }
    // Reset to initial node count
    d_node = 41;
    // Assign weights to dependent nodes on 2 independent nodes
    for (int jdex = 3; jdex < dep_nodes-3; jdex+=2){
      if ((jdex-1) % 10 == 0){
        jdex += 2;
        d_node++;
      }
      /* printf("dep_node: %d d_node: %d %d\n",-jdex-1, d_node,d_node+nx); */
      // Set the dependent weights
      dep_weights[4*jdex] = 0.5;
      dep_weights[4*jdex+1] = 0.5;
      dep_weights[4*jdex+2] = 0.0;
      dep_weights[4*jdex+3] = 0.0;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+nx;
      dep_conn[4*jdex+2] = d_node+1;
      dep_conn[4*jdex+3] = d_node+1+nx;
      d_node++;
    }
    
    // Reset to initial node count
    d_node = 41;
    // Assign weights to dependents nodes on 4 independent nodes
    for (int jdex = 4; jdex < dep_nodes-3; jdex+=2){
      if (jdex % 10 == 0){
        jdex += 4;
        d_node+=2;
      }
      /* printf("dep_node: %d d_node: %d %d %d %d\n",-jdex-1, */
      /*        d_node,d_node+1, d_node+nx, d_node+1+nx); */
      // Set the dependent weights
      dep_weights[4*jdex] = 0.25;
      dep_weights[4*jdex+1] = 0.25;
      dep_weights[4*jdex+2] = 0.25;
      dep_weights[4*jdex+3] = 0.25;
      // Set the dependent node connectivity
      dep_conn[4*jdex] = d_node;
      dep_conn[4*jdex+1] = d_node+1;
      dep_conn[4*jdex+2] = d_node+nx;
      dep_conn[4*jdex+3] = d_node+1+nx;
      d_node++;
    }
   
    for (int i = 0; i < dep_nodes; i++){
      dep_ptr[i+1] = dep_ptr[i] + 4;
    }
    // Set the length along x and y direction
    double Lx = 5.0, Ly = 5.0;

    // Create the nodes
    TacsScalar *Xpts = new TacsScalar[3*num_nodes];
    
    // Loop over the a-nodes
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
    // Loop over the b-nodes
    for (int j = 0; j < ny+1; j++){
      for (int i = 0; i < nx+1; i++){
        int node = bnodes[i+(nx+1)*j];
        if (node >= 0){
          Xpts[3*node] = Lx+Lx*i/nx;
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
    // Add the element from the a-block
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        /* printf("Anodes: %d %d %d %d \n",  */
        /*        anodes[i+(nx+1)*j], anodes[i+1+(nx+1)*j], */
        /*        anodes[i+1+(nx+1)*(j+1)], anodes[i+(nx+1)*(j+1)]); */
        conn[0] = anodes[i+(nx+1)*j];
        conn[1] = anodes[i+1+(nx+1)*j];
        conn[2] = anodes[i+(nx+1)*(j+1)];
        conn[3] = anodes[i+1+(nx+1)*(j+1)];
        conn += 4;
        elem_node_ptr[n+1] = elem_node_ptr[n]+4;
        n++;
      }
    }
    // Add the element from the b-block
    for (int j = 0; j < ny; j++){
      for (int i = 0; i < nx; i++){
        /* printf("Bnodes: %d %d %d %d \n",  */
        /*        bnodes[i+(nx+1)*j], bnodes[i+1+(nx+1)*j], */
        /*        bnodes[i+1+(nx+1)*(j+1)], bnodes[i+(nx+1)*(j+1)]); */
        conn[0] = bnodes[i+(nx+1)*j];
        conn[1] = bnodes[i+1+(nx+1)*j];
        conn[2] = bnodes[i+(nx+1)*(j+1)];
        conn[3] = bnodes[i+1+(nx+1)*(j+1)];
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
    
    /* for (int i = 0; i <= dep_nodes; i++){ */
    /* /\*   printf("conn[%d]: %d %d %d %d \n",-i-1,  *\/ */
    /* /\*          dep_conn[4*i], dep_conn[4*i+1],dep_conn[4*i+2],dep_conn[4*i+3]); *\/ */
    /* /\*   printf("weights[%d]: %e %e %e %e \n",-i-1, *\/ */
    /* /\*          dep_weights[4*i],dep_weights[4*i+1],dep_weights[4*i+2],dep_weights[4*i+3]); *\/ */
      
    /*   printf("ptr[%d]: %d \n", i, dep_ptr[i]); */
    /* } */
    /* exit(0); */
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
    delete [] bnodes;
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
