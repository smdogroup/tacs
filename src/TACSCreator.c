#include "TACSCreator.h"
#include "FElibrary.h"
#include "tacsmetis.h"

/*
  Extend an integer array 
*/
static void extend_int_array( int ** array, int old_len, 
			      int new_len ){
  int * temp = new int[ new_len ];
  memcpy(temp, *array, sizeof(int)*old_len);
  delete [] *array;
  *array = temp;
}

/*
  Functions for sorting a list such that:

  arg_sort_list[list[i]] is in ascending order
*/
static const int * arg_sort_list = NULL;

static int compare_arg_sort( const void * a, const void * b ){
  return arg_sort_list[*(int*)a] - arg_sort_list[*(int*)b];
}

/*
  Allocate the TACSCreator object
*/
TACSCreator::TACSCreator( MPI_Comm _comm, 
			  int _vars_per_node ){
  // Set the communicator
  comm = _comm;
  root_rank = 0;

  // Set the number of variables per node in the mesh
  vars_per_node = _vars_per_node;

  // Element connectivity information
  num_nodes = 0;
  num_elements = 0;
  elem_id_nums = NULL;
  elem_node_ptr = NULL;
  elem_node_conn = NULL;

  // Set the initial point locations
  Xpts = NULL;

  // Boundary condition information
  bc_nodes = NULL;
  bc_vars = NULL;
  bc_ptr = NULL;
}

/*
  Deallocate the TACSCreator object
*/
TACSCreator::~TACSCreator(){
  if (elem_id_nums){ delete [] elem_id_nums; }
  if (elem_node_ptr){ delete [] elem_node_ptr; }
  if (elem_node_conn){ delete [] elem_node_conn; }
  if (Xpts){ delete [] Xpts; }
  if (bc_nodes){ delete [] bc_nodes; }
  if (bc_vars){ delete [] bc_vars; }
  if (bc_ptr){ delete [] bc_ptr; }
}

/*
  Set the global element connectivity arrays into the class
*/
void TACSCreator::setGlobalConnectivity( int _num_nodes, int _num_elements,
					 int *_elem_node_ptr, 
					 int *_elem_node_conn,
					 int *_elem_id_nums ){
  num_elements = _num_elements;
  num_nodes = _num_nodes;

  // Copy the pointer array into the elem_node_conn array
  elem_node_ptr = new int[ num_elements+1 ];
  memcpy(elem_node_ptr, _elem_node_ptr, (num_elements+1)*sizeof(int));

  // Copy over the element node connectivity
  elem_node_conn = new int[ elem_node_ptr[num_elements] ];
  memcpy(elem_node_conn, _elem_node_conn, 
	 elem_node_ptr[num_elements]*sizeof(int));

  // Copy over the element ids
  elem_id_nums = new int[ num_elements ];
  memcpy(elem_id_nums, _elem_id_nums, num_elements*sizeof(int));
}

/*
  Set the boundary condition data
*/
void TACSCreator::setBoundaryConditions( int _num_bcs, 
					 int *_bc_nodes, int *_bc_vars,
					 int *_bc_ptr ){
  // Set the number of boundary conditions and the node numbers
  // to which they correspond
  num_bcs = _num_bcs;
  bc_nodes = new int[ num_bcs ];
  memcpy(bc_nodes, _bc_nodes, num_bcs*sizeof(int));

  // Allocate space for the boundary condition pointer
  bc_ptr = new int[ num_bcs+1 ];

  if (_bc_vars){
    // Copy over the boundary condition variable pointer
    memcpy(bc_ptr, _bc_ptr, (num_bcs+1)*sizeof(int));
  
    // Allocate the number of variables
    bc_vars = new int[ bc_ptr[num_bcs] ];
    memcpy(bc_vars, _bc_vars, bc_ptr[num_bcs]*sizeof(int));
  }
  else {
    // Since the bc_vars array is input as NULL, assume that
    // all the variables at this node are fully restrained.
    bc_vars = new int[ num_bcs*vars_per_node ];
    bc_ptr[0] = 0;
    for ( int i = 0; i < num_bcs; i++ ){
      bc_ptr[i+1] = bc_ptr[i] + vars_per_node;

      // Restrain all the variables at the given node
      for ( int j = 0; j < vars_per_node; j++ ){
	bc_vars[bc_ptr[i]+j] = j;
      }
    }
  }
}

/*
  Create the instance of the TACSAssembler object and return it.

  This code partitions the mesh, calls for the elements to be
  allocated and returns a valid instance of the TACSAssembler object.

*/
TACSAssembler *TACSCreator::createTACS( enum TACSAssembler::OrderingType order_type,
					enum TACSAssembler::MatrixOrderingType mat_type ){
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // These arrays are only significant on the root processor
  int *partition = NULL, *new_nodes = NULL;
  int *owned_elements = NULL, *owned_nodes = NULL;

  if (rank == root_rank){   
    // Partition the mesh using the serial code on the root
    // processor.
    partition = new int[ num_elements ];
    new_nodes = new int[ num_nodes ];
    owned_elements = new int[size];
    owned_nodes = new int[size];
    
    splitSerialMesh(size, partition, new_nodes, 
		    owned_elements, owned_nodes);    
  }
  
  // The number of local and owned nodes and the number of
  // elements for each processor in the mesh
  int num_local_nodes = 0;
  int num_owned_nodes = 0;
  int num_owned_elements = 0; 
  MPI_Scatter(owned_nodes, 1, MPI_INT, 
              &num_owned_nodes, 1, MPI_INT, root_rank, comm);
  MPI_Scatter(owned_elements, 1, MPI_INT, 
              &num_owned_elements, 1, MPI_INT, root_rank, comm);

  // Local element connectivity information
  int *local_elem_id_nums = new int[ num_owned_elements ];

  // Allocate space for the portion of the element connectivity on 
  // this processor
  int *local_elem_node_ptr = new int[ num_owned_elements+1 ];
  int *local_elem_node_conn = NULL;

  // Loacal nodal information
  int *local_tacs_nodes = NULL;
  double *Xpts_local = NULL;

  // For each processor, send the information to the owner
  if (rank == root_rank){
    // Reset the nodes for the boundary conditions so that they
    // correspond to the new ordering
    for ( int j = 0; j < num_bcs; j++ ){
      bc_nodes[j] = new_nodes[bc_nodes[j]];
    }

    // Find the inverse mapping between the new and old node
    // numbers so that it's faster to access
    int *inv_new_nodes = new int[ num_nodes ];
    for ( int i = 0; i < num_nodes; i++ ){
      inv_new_nodes[new_nodes[i]] = i;
    }

    // Set up the element partition so that it is sorted so that
    // we can find 
    int *elem_part = new int[ num_elements ];
    for ( int k = 0; k < num_elements; k++ ){
      elem_part[k] = k;
    }

    // Now partition[elem_part[k]] is sorted
    arg_sort_list = partition;
    qsort(elem_part, num_elements, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Compute the local CSR data structure on this process
    // Use an upper bound for the memory requirements
    int *tacs_nodes = new int[ elem_node_ptr[num_elements] ];
    int *elem_ptr = new int[ num_elements+1 ];
    int *elem_conn = new int[ elem_node_ptr[num_elements] ];


    int *elem_ids = new int[ num_elements ];
    TacsScalar * xpts = new TacsScalar[ 3*num_nodes ];

    int start = 0;
    for ( int k = 0; k < size; k++ ){
      int end = start + owned_elements[k];

      int local_node_size = 0;
      elem_ptr[0] = 0;

      // Cycle through partition k
      int n = 0;
      for ( int j = start; j < end; j++, n++ ){
        int elem = elem_part[j];

	// Add the element
        for ( int i = elem_node_ptr[elem]; 
	      i < elem_node_ptr[elem+1]; i++, local_node_size++ ){
          int node = elem_node_conn[i];
          tacs_nodes[local_node_size] = new_nodes[node];
	  elem_conn[local_node_size] = new_nodes[node];          
        }

	elem_comp[n] = elem_component[elem];
	elem_ptr[n+1] = local_node_size;
      }

      // Uniquify the list
      local_node_size = FElibrary::uniqueSort(tacs_nodes, local_node_size);

      // tacs_nodes is sorted and defines the local ordering 
      // tacs_nodes[i] -> local node i      
      for ( int j = 0; j < elem_ptr[n]; j++ ){
	int * item = (int*)bsearch(&elem_conn[j], tacs_nodes, local_node_size,
				   sizeof(int), FElibrary::comparator);
	if (item){
	  elem_conn[j] = item - tacs_nodes;
	}
	else {
	  fprintf(stderr, 
		  "[%d] TACSCreator could not find %d node in local list\n",
		  rank, elem_conn[j]);
	  MPI_Abort(comm, 1);
	}
      }

      // Copy over the data from the local node numbers
      for ( int j = 0; j < local_node_size; j++ ){
	int node = inv_new_nodes[tacs_nodes[j]];
	xpts[3*j] = Xpts[3*node];
	xpts[3*j+1] = Xpts[3*node+1];
	xpts[3*j+2] = Xpts[3*node+2];
      }
    
      // Create the CSR data structure required
      if (k == root_rank){
	// Copy the values over from the nodes
	num_local_nodes = local_node_size;
	local_tacs_nodes = new int[ num_local_nodes ];
	Xpts_local = new double[ 3*num_local_nodes ];
	
	memcpy(local_tacs_nodes, tacs_nodes, num_local_nodes*sizeof(int));
	memcpy(Xpts_local, xpts, 3*num_local_nodes*sizeof(double));

	// Copy over values for the elements
	memcpy(local_component_nums, elem_comp, 
	       num_owned_elements*sizeof(int));
	memcpy(local_elem_node_ptr, elem_ptr, 
	       (num_owned_elements+1)*sizeof(int));

	local_elem_node_con = new int[local_elem_node_ptr[num_owned_elements]];
	memcpy(local_elem_node_con, elem_conn,
	       local_elem_node_ptr[num_owned_elements]*sizeof(int));
      }
      else {
        // Send the data to the other process
	MPI_Send(&local_node_size, 1, MPI_INT, k, 1, comm);

	MPI_Send(tacs_nodes, local_node_size, MPI_INT, k, 2, comm);
	MPI_Send(xpts, 3*local_node_size, MPI_DOUBLE, k, 3, comm);

	// Send the element data
        MPI_Send(elem_comp, owned_elements[k], MPI_INT, k, 4, comm);
        MPI_Send(elem_ptr, owned_elements[k]+1, MPI_INT, k, 5, comm);

        MPI_Send(elem_conn, elem_ptr[owned_elements[k]], MPI_INT, k, 6, comm);
      }

      start += owned_elements[k];
    }

    delete [] elem_part;
    delete [] inv_new_nodes;
    delete [] tacs_nodes;
    delete [] elem_ptr;
    delete [] elem_conn;
    delete [] elem_comp;
    delete [] xpts;
  }
  else {
    // Recv the data from the root process
    MPI_Status status;
    MPI_Recv(&num_local_nodes, 1, MPI_INT, 
	     root_rank, 1, comm, &status);

    // Allocate space for the incoming data
    local_tacs_nodes = new int[ num_local_nodes ];
    Xpts_local = new double[ 3*num_local_nodes ];

    MPI_Recv(local_tacs_nodes, num_local_nodes, MPI_INT, 
	     root_rank, 2, comm, &status);
    MPI_Recv(Xpts_local, 3*num_local_nodes, MPI_DOUBLE, 
	     root_rank, 3, comm, &status);

    // Receive the element data
    MPI_Recv(local_component_nums, num_owned_elements, MPI_INT, 
	     root_rank, 4, comm, &status);
    MPI_Recv(local_elem_node_ptr, num_owned_elements+1, MPI_INT, 
	     root_rank, 5, comm, &status);

    int con_size = local_elem_node_ptr[num_owned_elements];
    local_elem_node_con = new int[con_size];
    MPI_Recv(local_elem_node_con, con_size, MPI_INT, 
	     root_rank, 6, comm, &status);
  }
  
  int node_max_csr_size = local_elem_node_ptr[num_owned_elements];  
  TACSAssembler * tacs = new TACSAssembler(comm, num_owned_nodes, vars_per_node,
                                           num_owned_elements, num_local_nodes,
                                           node_max_csr_size);

  // Sort out the boundary conditions
  // Broadcast the boundary condition information
  MPI_Bcast(bc_nodes, num_bcs, MPI_INT, root_rank, comm);

  // Get the local node numbers for the boundary conditions
  for ( int k = 0; k < num_bcs; k++ ){
    int * item = (int*)bsearch(&bc_nodes[k], local_tacs_nodes, num_local_nodes, 
			       sizeof(int), FElibrary::comparator);
    if (item){
      bc_nodes[k] = item - local_tacs_nodes;
    }
    else {
      bc_nodes[k] = -1;
    }
  }

  // Add the node numbers - this steals the reference
  tacs->addNodes(&local_tacs_nodes);
  
  // Add the elements
  for ( int k = 0; k < num_owned_elements; k++ ){
    TACSElement * element = elements[local_component_nums[k]];
    if (!element){
      fprintf(stderr, 
              "[%d] TACSMeshLoader: Element undefined for component %d\n",
              rank, local_component_nums[k]);
      MPI_Abort(comm, 1);
      return NULL;
    }

    // Add the element node numbers
    int start = local_elem_node_ptr[k];
    int end = local_elem_node_ptr[k+1];
    tacs->addElement(element, &local_elem_node_con[start], end-start);
  }

  tacs->computeReordering(order_type, mat_type);

  // Finalize the ordering
  tacs->finalize();

  // Set the nodes
  TacsScalar * x;
  tacs->getNodeArray(&x);
  
  for ( int k = 0; k < 3*num_local_nodes; k++ ){
    x[k] = Xpts_local[k];
  }

  // Set the boundar conditions
  int bvars[6];
  TacsScalar bvals[6];
  for ( int k = 0; k < num_bcs; k++ ){
    if (bc_nodes[k] >= 0){
      int nbcs = bc_ptr[k+1] - bc_ptr[k];
      int n = 0;
      for ( int j = 0; j < nbcs; j++ ){
        if (bc_con[bc_ptr[k] + j] < vars_per_node){
          bvars[n] = bc_con[bc_ptr[k] + j];
          bvals[n] = bc_vals[bc_ptr[k] + j];
          n++;
        }
      }
      if (n > 0){
        tacs->addBC(bc_nodes[k], bvars, bvals, n);
      }
    }
  }
  
  if (rank == root_rank){
    delete [] new_nodes;
    delete [] owned_elements;
    delete [] owned_nodes;
  }

  delete [] partition;
  delete [] local_elem_node_ptr;
  delete [] local_elem_node_con;
  delete [] Xpts_local;

  return tacs;
}

/*
  Partition the mesh stored on the root processor for parallel
  computations.

  This function is only called on the root processor. The function
  first forms the dual mesh with an element->element data structure.
  The function then calls METIS to partition the mesh.

  input:
  split_size:      the number of segments in the partition

  output:
  elem_partition:  the element->processor assignment
  new_nodes:       node index i has a new index new_nodes[i]
  owned_elements:  number of elements owned by processor i
  owned_nodes:     number of nodes owned by processor i
*/
void TACSCreator::splitSerialMesh( int split_size, 
				   int *elem_partition, 
				   int *new_nodes,
				   int *owned_elements, 
				   int *owned_nodes ){
  // Compute the node to element CSR data structure
  // node_elem_conn[node_elem_ptr[node]:node_elem_ptr[node+1]] are the
  // elements that contain node, 'node'
  int *node_elem_ptr = new int[ num_nodes+1 ];
  memset(node_elem_ptr, 0, (num_nodes+1)*sizeof(int));
    
  for ( int i = 0; i < num_elements; i++ ){
    int end = elem_node_ptr[i+1]; 
    for ( int j = elem_node_ptr[i]; j < end; j++ ){
      int node = elem_node_conn[j];
      node_elem_ptr[node+1]++;
    }
  }

  // Determine the size of the node to element array
  for ( int i = 0; i < num_nodes; i++ ){
    node_elem_ptr[i+1] += node_elem_ptr[i];
  }
  int *node_elem_conn = new int[ node_elem_ptr[num_nodes] ];

  // Fill in the entries of the node->element data structure
  for ( int i = 0; i < num_elements; i++ ){
    int end = elem_node_ptr[i+1];
    for ( int j = elem_node_ptr[i]; j < end; j++ ){
      int node = elem_node_conn[j];
      node_elem_conn[node_elem_ptr[node]] = i;
      node_elem_ptr[node]++;
    }
  }
  
  // Reset the node_elem_ptr array to the correct range
  for ( int i = num_nodes; i > 0; i-- ){
    node_elem_ptr[i] = node_elem_ptr[i-1];
  }
  node_elem_ptr[0] = 0;

  // Now we set up the element to element connectivity.  For this to
  // work within METIS, you have to remove the diagonal contribution
  // so that there is no self-reference in the graph.
  int *elem_ptr = new int[ num_elements+1 ];
  elem_ptr[0] = 0;

  // Information to keep track of how big the data structure is
  int elem_conn_size = 0;

  // Estimate the maximum size of the connectivity data
  const int ROW_SIZE_EST = 27;
  int max_elem_conn_size = ROW_SIZE_EST*num_elements;
  int *elem_conn = new int[ max_elem_conn_size ];

  // Assemble things one row at a time
  int *row = new int[ num_elements ];

  // Go through the elements and merge the arrays
  for ( int i = 0; i < num_elements; i++ ){
    int row_size = 0;
    
    // Add the element -> element connectivity
    for ( int j = elem_node_ptr[i]; j < elem_node_ptr[i+1]; j++ ){
      int node = elem_node_conn[j];

      // Extend the row array while keeping it sorted 
      int start = node_elem_ptr[node];
      int size = node_elem_ptr[node+1] - start;
      row_size = FElibrary::mergeArrays(row, row_size, 
                                        &node_elem_conn[start], size);
    }

    // Check if adding this new row to the data structure will excee
    if (elem_conn_size + row_size > max_elem_conn_size){
      max_elem_conn_size += 0.5*max_elem_conn_size;
      if (max_elem_conn_size < elem_conn_size + row_size){
        max_elem_conn_size += elem_conn_size + row_size;
      }
      extend_int_array(&elem_conn, elem_conn_size, 
                       max_elem_conn_size);
    }

    // Add the elements - minus the diagonal entry
    for ( int j = 0; j < row_size; j++ ){
      if (row[j] != i){ // Not the diagonal 
        elem_con[elem_con_size] = row[j];
        elem_con_size++;
      }
    }
    elem_ptr[i+1] = elem_con_size;
  }

  // Free the memory for data that is not needed 
  delete [] row;
  delete [] node_elem_ptr;
  delete [] node_elem_con;

  // Partition the mesh using METIS.
  if (split_size > 1){
    int options[5];
    options[0] = 0; // use the default options
    int wgtflag = 0; // weights are on the verticies
    int numflag = 0; // C style numbering 
    int edgecut = -1;        
      
    int *vwgts = NULL; // Weights on the vertices 
    int *adjwgts = NULL;  // Weights on the edges or adjacency
      
    if (split_size < 8){
      METIS_PartGraphRecursive(&num_elements, elem_ptr, elem_conn, 
                               vwgts, adjwgts, 
                               &wgtflag, &numflag, &split_size, 
                               options, &edgecut, elem_partition);
    }
    else {
      METIS_PartGraphKway(&num_elements, elem_ptr, elem_conn, 
                          vwgts, adjwgts, 
                          &wgtflag, &numflag, &split_size, 
                          options, &edgecut, elem_partition);
    }
  }
  else {
    // If there is no split, just assign all elements to the
    // root processor
    for ( int k = 0; k < num_elements; k++ ){
      elem_partition[k] = 0;
    }
  }

  delete [] elem_conn;
  delete [] elem_ptr;

  // Now, re-order the variables so that they are almost contiguous
  // over each processor 
  memset(owned_nodes, 0, split_size*sizeof(int));
  memset(owned_elements, 0, split_size*sizeof(int));

  // Set up an array that tracks the mapping from the old
  // node ordering to the new node ordering. The array stores
  // the information such that the old node with index i is
  // stored at the index new_nodes[i].

  // First, treat new_nodes as an array of flags that indicate
  // whether this node has been counted yet.
  memset(new_nodes, 0, num_nodes*sizeof(int));

  // Keep a count that is the new global ordering of the nodes
  // within the finite-element mesh.
  int count = 0;
  for ( int j = 0; j < num_elements; j++ ){
    // The owner of element j
    int owner = elem_partition[j];
    owned_elements[owner]++;

    // Assign the un-assigned nodes associated with element j
    // to the owner of element j.
    for ( int i = elem_node_ptr[j]; i < elem_node_ptr[j+1]; i++ ){
      int node = elem_node_conn[i];
      if (!new_nodes[node]){
        new_nodes[node] = 1;
        owned_nodes[owner]++;
        count++;
      }
    }
  }

  // Reset the new_nodes array and prepare to assign the new
  // node numbers in a contiguous fashion.
  for ( int k = 0; k < num_nodes; k++ ){
    new_nodes[k] = -1;
  }

  // Find the offset to the ordering of the nodes for each partition
  // such that split_offset[i] is the first node number for
  // the i-th processor
  int *split_offset = new int[ split_size ];
  split_offset[0] = 0;
  for ( int k = 1; k < split_size; k++ ){
    split_offset[k] = split_offset[k-1] + owned_nodes[k-1];
  }

  // Now assign the new_nodes with the correct offsets. As a result,
  // the nodes will be ordered properly.
  for ( int j = 0; j < num_elements; j++ ){
    int owner = elem_partition[j];
    for ( int i = elem_node_ptr[j]; i < elem_node_ptr[j+1]; i++ ){
      int node = elem_node_con[i];
      if (new_nodes[node] < 0){
        new_nodes[node] = split_offset[owner];
        split_offset[owner]++;
      }
    }
  }

  delete [] split_offset;
}

