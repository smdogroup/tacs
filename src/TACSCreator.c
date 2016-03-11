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
  num_dependent_nodes = 0;
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

  if (elements){
    for ( int i = 0; i < num_elem_ids; i++ ){
      if (elements[i]){ elements[i]->decref(); }
    }
    delete [] elements;
  }
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
  Set the elements
*/
void TACSCreator::setElements( TACSElement **_elements, int _num_elem_ids ){
  num_elem_ids = _num_elem_ids;
  elements = new TACSElement*[ num_elem_ids ];
  memcpy(elements, _elements, num_elem_ids*sizeof(TACSElement*));
  for ( int i = 0; i < num_elem_ids; i++ ){
    if (elements[i]){ elements[i]->incref(); }
  }
}

/*
  Set the nodal locations
*/
void TACSCreator::setNodes( TacsScalar *_Xpts ){
  Xpts = new TacsScalar[ 3*num_nodes ];
  memcpy(Xpts, _Xpts, 3*num_nodes*sizeof(TacsScalar));
}

/*
  Create the instance of the TACSAssembler object and return it.

  This code partitions the mesh, calls for the elements to be
  allocated and returns a valid instance of the TACSAssembler object.

*/
TACSAssembler* TACSCreator::createTACS( enum TACSAssembler::OrderingType order_type,
					enum TACSAssembler::MatrixOrderingType mat_type ){
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // These arrays are only significant on the root processor
  int *partition = NULL, *new_nodes = NULL;
  int *new_dep_nodes = NULL, *owned_elements = NULL;
  int *owned_nodes = NULL, *owned_dep_nodes = NULL;

  if (rank == root_rank){   
    // Partition the mesh using the serial code on the root
    // processor.
    partition = new int[ num_elements ];
    new_nodes = new int[ num_nodes ];
    new_dep_nodes = new int[ num_dependent_nodes ];
    owned_elements = new int[ size ];
    owned_nodes = new int[ size ];
    owned_dep_nodes = new int[ size ];
    
    splitSerialMesh(size, partition, new_nodes, new_dep_nodes,
		    owned_elements, owned_nodes, owned_dep_nodes);    
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
    // we can quickly find the elements associated with a processor
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
    TacsScalar *xpts = new TacsScalar[ 3*num_nodes ];

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

	elem_ids[n] = elem_id_nums[elem];
	elem_ptr[n+1] = local_node_size;
      }

      // Uniquify the list
      local_node_size = FElibrary::uniqueSort(tacs_nodes, local_node_size);

      // tacs_nodes is sorted and defines the local ordering 
      // tacs_nodes[i] -> local node i      
      for ( int j = 0; j < elem_ptr[n]; j++ ){
	int * item = (int*)bsearch(&elem_conn[j], tacs_nodes, 
				   local_node_size,
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
	memcpy(local_elem_id_nums, elem_ids, 
	       num_owned_elements*sizeof(int));
	memcpy(local_elem_node_ptr, elem_ptr, 
	       (num_owned_elements+1)*sizeof(int));

	local_elem_node_conn = new int[local_elem_node_ptr[num_owned_elements]];
	memcpy(local_elem_node_conn, elem_conn,
	       local_elem_node_ptr[num_owned_elements]*sizeof(int));
      }
      else {
        // Send the data to the other process
	MPI_Send(&local_node_size, 1, MPI_INT, k, 1, comm);

	MPI_Send(tacs_nodes, local_node_size, MPI_INT, k, 2, comm);
	MPI_Send(xpts, 3*local_node_size, MPI_DOUBLE, k, 3, comm);

	// Send the element data
        MPI_Send(elem_ids, owned_elements[k], MPI_INT, k, 4, comm);
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
    delete [] elem_ids;
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
    MPI_Recv(local_elem_id_nums, num_owned_elements, MPI_INT, 
	     root_rank, 4, comm, &status);
    MPI_Recv(local_elem_node_ptr, num_owned_elements+1, MPI_INT, 
	     root_rank, 5, comm, &status);

    int con_size = local_elem_node_ptr[num_owned_elements];
    local_elem_node_conn = new int[con_size];
    MPI_Recv(local_elem_node_conn, con_size, MPI_INT, 
	     root_rank, 6, comm, &status);
  }
  
  int node_max_csr_size = local_elem_node_ptr[num_owned_elements];  
  TACSAssembler * tacs = new TACSAssembler(comm, num_owned_nodes, vars_per_node,
                                           num_owned_elements, num_local_nodes,
                                           node_max_csr_size);

  // Add the node numbers - this steals the reference
  tacs->addNodes(&local_tacs_nodes);
  
  // Add the elements
  for ( int k = 0; k < num_owned_elements; k++ ){
    TACSElement * element = elements[local_elem_id_nums[k]];
    if (!element){
      fprintf(stderr, 
              "[%d] TACSMeshLoader: Element undefined for component %d\n",
              rank, local_elem_id_nums[k]);
      MPI_Abort(comm, 1);
      return NULL;
    }

    // Add the element node numbers
    int start = local_elem_node_ptr[k];
    int end = local_elem_node_ptr[k+1];
    tacs->addElement(element, &local_elem_node_conn[start], end-start);
  }

  tacs->computeReordering(order_type, mat_type);

  // Finalize the ordering
  tacs->finalize();

  // Set the node locations in TACS
  TacsScalar * x;
  tacs->getNodeArray(&x);
  for ( int k = 0; k < 3*num_local_nodes; k++ ){
    x[k] = Xpts_local[k];
  }

  // Broadcast the boundary condition information
  MPI_Bcast(&num_bcs, 1, MPI_INT, root_rank, comm);

  // Broadcast all the boundary conditions to all processors
  if (rank != root_rank){
    if (bc_nodes){ delete [] bc_nodes; }
    bc_nodes = new int[ num_bcs ];
  }
  MPI_Bcast(bc_nodes, num_bcs, MPI_INT, root_rank, comm);

  // Broacast the boundary condition pointer array 
  if (rank != root_rank){
    if (bc_ptr){ delete [] bc_ptr; }
    bc_ptr = new int[ num_bcs+1 ];
  }
  MPI_Bcast(bc_ptr, num_bcs+1, MPI_INT, root_rank, comm);

  if (rank != root_rank){
    if (bc_vars){ delete [] bc_vars; }
    bc_vars = new int[ bc_ptr[num_bcs] ]; 
  }
  MPI_Bcast(bc_vars, bc_ptr[num_bcs], MPI_INT, root_rank, comm);
  
  // Get the local node numbers for the boundary conditions
  const int *tacs_nodes = NULL;
  tacs->getTacsNodeNums(&tacs_nodes);

  for ( int k = 0; k < num_bcs; k++ ){
    int * item = (int*)bsearch(&bc_nodes[k], tacs_nodes, num_local_nodes, 
			       sizeof(int), FElibrary::comparator);
    if (item){
      bc_nodes[k] = item - local_tacs_nodes;
    }
    else {
      bc_nodes[k] = -1;
    }
  }

  // Allocate the arrays to store the variable values
  int *bvars = new int[ vars_per_node ];
  TacsScalar *bvals = new TacsScalar[ vars_per_node ];

  // Set the boundary conditions
  for ( int k = 0; k < num_bcs; k++ ){
    if (bc_nodes[k] >= 0){
      int nbcs = bc_ptr[k+1] - bc_ptr[k];
      int n = 0;
      for ( int j = 0; j < nbcs; j++ ){
        if (bc_vars[bc_ptr[k] + j] < vars_per_node){
          bvars[n] = bc_vars[bc_ptr[k] + j];
          bvals[n] = 0.0;
          n++;
        }
      }
      if (n > 0){
        tacs->addBC(bc_nodes[k], bvars, bvals, n);
      }
    }
  }

  // Free the bvars/bvals arrays
  delete [] bvars;
  delete [] bvals;
  
  // Free memory only allocated on the root processor
  if (rank == root_rank){
    delete [] new_nodes;
    delete [] new_dep_nodes;
    delete [] owned_elements;
    delete [] owned_nodes;
    delete [] owned_dep_nodes;
  }

  // Free all the remaining memory
  delete [] partition;
  delete [] local_elem_node_ptr;
  delete [] local_elem_node_conn;
  delete [] local_elem_id_nums;
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
  partition:       the element->processor assignment
  new_nodes:       node index i has a new index new_nodes[i]
  new_dep_nodes:   dependent node i has new index new_dep_nodes[i]
  owned_elements:  number of elements owned by processor i
  owned_nodes:     number of nodes owned by processor i
  owned_dep_nodes: number of dependent nodes owned by processor i
*/
void TACSCreator::splitSerialMesh( int split_size, 
				   int *partition, 
				   int *new_nodes,
				   int *new_dep_nodes,
				   int *owned_elements, 
				   int *owned_nodes,
				   int *owned_dep_nodes ){
  // Compute the node to element CSR data structure for both
  // the indepedent and dependent nodes
  int *node_elem_ptr = new int[ num_nodes+1 ];
  int *dep_node_elem_ptr = new int[ num_dependent_nodes+1 ];
  memset(node_elem_ptr, 0, (num_nodes+1)*sizeof(int));
  memset(dep_node_elem_ptr, 0, (num_dependent_nodes+1)*sizeof(int));
    
  for ( int i = 0; i < num_elements; i++ ){
    int end = elem_node_ptr[i+1]; 
    for ( int j = elem_node_ptr[i]; j < end; j++ ){
      int node = elem_node_conn[j];
      if (node >= 0){
	// This is an independent node
	node_elem_ptr[node+1]++;
      }
      else {
	// This is a dependent node. Compute the dependent node
	// number and add it to the dependent node connectivity
	// counter
	node = -node-1;
	dep_node_elem_ptr[node+1]++;
      }
    }
  }

  // Determine the size of the node to element array
  for ( int i = 0; i < num_nodes; i++ ){
    node_elem_ptr[i+1] += node_elem_ptr[i];
  }
  int *node_elem_conn = new int[ node_elem_ptr[num_nodes] ];

  // Determine the size of the dependent node to element array
  for ( int i = 0; i < num_dependent_nodes; i++ ){
    dep_node_elem_ptr[i+1] += dep_node_elem_ptr[i];
  }
  int *dep_node_elem_conn = new int[ dep_node_elem_ptr[num_dependent_nodes] ];

  // Fill in the entries of the node->element data structure
  for ( int i = 0; i < num_elements; i++ ){
    int end = elem_node_ptr[i+1];
    for ( int j = elem_node_ptr[i]; j < end; j++ ){
      int node = elem_node_conn[j];
      if (node >= 0){
	// This is an independent node
	node_elem_conn[node_elem_ptr[node]] = i;
	node_elem_ptr[node]++;
      }
      else {
	// This is a dependent node. Compute the dependent node
	// number and add the element to the connectivity list
	node = -node-1;
	dep_node_elem_conn[dep_node_elem_ptr[node]] = i;
	dep_node_elem_ptr[node]++;
      }
    }
  }
  
  // Reset the node_elem_ptr array to the correct values
  for ( int i = num_nodes; i > 0; i-- ){
    node_elem_ptr[i] = node_elem_ptr[i-1];
  }
  node_elem_ptr[0] = 0;

  // Reset the dep_elem_ptr array to the correct values
  for ( int i = num_dependent_nodes; i > 0; i-- ){
    dep_node_elem_ptr[i] = dep_node_elem_ptr[i-1];
  }
  dep_node_elem_ptr[0] = 0;

  // Now set up an element mask so that elements linked through a
  // common dependent node must lie on the same processor
  int *elem_mask = new int[ num_elements ];
  int *inv_elem_mask = new int[ num_elements ];
  int *elem_queue = new int[ num_elements ];
  for ( int i = 0; i < num_elements; i++ ){
    elem_mask[i] = -1;
  }
  
  // Search through the dependent nodes and label all elements that
  // share a common dependent node
  int num_masked_elems = 0;
  int inv_mask_index = 0;
  for ( int i = 0; i < num_dependent_nodes; i++ ){
    int queue_count = 0;
    for ( int jp = dep_node_elem_ptr[i]; 
	  jp < dep_node_elem_ptr[i+1]; jp++ ){

      // Get the dependent element variable
      int elem = dep_node_elem_conn[jp];
      if (elem_mask[elem] < 0){
	elem_queue[queue_count] = elem;
	elem_mask[elem] = num_masked_elems;
	queue_count++;

	inv_elem_mask[inv_mask_index] = elem;
	inv_mask_index++;
      }
    }

    // Increment the mask number only if we find a new
    // dependent node that is not yet masked
    int new_flag = (queue_count > 0);

    while (queue_count > 0){
      // Search for adjacent elements that share another
      // common element
      queue_count--;
      int elem = elem_queue[queue_count];
      for ( int jp = elem_node_ptr[elem]; 
	    jp < elem_node_ptr[elem+1]; jp++ ){
	// Get the node and check if it is a dependent node
	int node = elem_node_conn[jp];

	if (node < 0){
	  // Compute the dependent node index
	  node = -node-1;
	
	  // Now find the elements associated with this dependent node
	  // and add them to the queue
	  for ( int kp = dep_node_elem_ptr[node]; 
		kp < dep_node_elem_ptr[node+1]; kp++ ){

	    // Get the dependent element variable
	    int new_elem = dep_node_elem_conn[kp];
	    if (elem_mask[new_elem] < 0){
	      elem_queue[queue_count] = new_elem;
	      elem_mask[new_elem] = num_masked_elems;
	      queue_count++;

	      inv_elem_mask[inv_mask_index] = new_elem;
	      inv_mask_index++;
	    }
	  }
	}
      }
    }

    // If a new unmasked element/node was found, then
    // increment the number of masked elements
    if (new_flag){
      num_masked_elems++;
    }
  }

  // Set the remaining nodes
  for ( int i = 0; i < num_elements; i++ ){
    if (elem_mask[i] < 0){
      elem_mask[i] = num_masked_elems;
      num_masked_elems++;

      inv_elem_mask[inv_mask_index] = i;
      inv_mask_index++;
    }
  }

  // Delete the element queue
  delete [] elem_queue;
  
  // Now that the element mask has been created, we don't need
  // the dependent element connectivity information anymore
  delete [] dep_node_elem_ptr;
  delete [] dep_node_elem_conn;

  // Now we set up the element to element connectivity using the 
  // set of maksed elements. For this to work within METIS, we have
  // to remove the diagonal contribution so that there is no
  // self-reference in the graph.
  int *elem_ptr = new int[ num_masked_elems+1 ];
  elem_ptr[0] = 0;

  // Information to keep track of how big the data structure is
  int elem_conn_size = 0;

  // Estimate the maximum size of the connectivity data
  const int ROW_SIZE_EST = 27;
  int max_elem_conn_size = ROW_SIZE_EST*num_elements;
  int *elem_conn = new int[ max_elem_conn_size ];

  // Assemble things one row at a time
  int *row = new int[ num_elements ];

  for ( int i = 0, index = 0; i < num_masked_elems; i++ ){    
    // Set the size of the new row in the data structure to zero
    int row_size = 0;

    // While the mask remains unchanged - still equal to the element
    // mask number i - continue to add the connectivity to the current
    // row without changing the mask number.
    while (index < num_elements &&
	   elem_mask[inv_elem_mask[index]] == i){
      int elem = inv_elem_mask[index];

      // From the element number, find the independent nodes 
      // that are associated with this element number
      int jpend = elem_node_ptr[elem+1];
      for ( int jp = elem_node_ptr[elem]; jp < jpend; jp++ ){
	int node = elem_node_conn[jp];

	if (node >= 0){
	  // Add all the masked components to the array
	  int kpend = node_elem_ptr[node+1];
	  for ( int kp = node_elem_ptr[node]; kp < kpend; kp++ ){

	    // Get the index of the new masked element
	    int new_elem = elem_mask[node_elem_conn[kp]];

	    // Check if adding another element will exceed the size
	    // of the array. If it will, sort the array and remove
	    // duplicates. This is guaranteed to make additional 
	    // room.
	    if (row_size >= num_elements){
	      row_size = FElibrary::uniqueSort(row, row_size);
	    }

	    // Append the masked element index to the list
	    row[row_size] = new_elem;
	    row_size++;
	  }
	}
      }
      
      // Increment the index counter
      index++;
    }

    // Sort the element indices and remove duplicates from the
    // array
    row_size = FElibrary::uniqueSort(row, row_size);

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
        elem_conn[elem_conn_size] = row[j];
        elem_conn_size++;
      }
    }
    elem_ptr[i+1] = elem_conn_size;
  }

  // Free the inverse element mask
  delete [] inv_elem_mask;

  // Free the memory for data that is not needed 
  delete [] row;
  
  // Free the node to element pointer
  delete [] node_elem_ptr;
  delete [] node_elem_conn;

  // Partition the mesh using METIS.
  if (split_size > 1){
    // Allocate the element partition
    int *elem_part = new int[ num_masked_elems ];
    
    int options[5];
    options[0] = 0; // use the default options
    int wgtflag = 0; // weights are on the verticies
    int numflag = 0; // C style numbering 
    int edgecut = -1; 
    int *vwgts = NULL; // Weights on the vertices 
    int *adjwgts = NULL;  // Weights on the edges or adjacency
      
    if (split_size < 8){
      METIS_PartGraphRecursive(&num_masked_elems, elem_ptr, elem_conn, 
                               vwgts, adjwgts, 
                               &wgtflag, &numflag, &split_size, 
                               options, &edgecut, elem_part);
    }
    else {
      METIS_PartGraphKway(&num_masked_elems, elem_ptr, elem_conn, 
                          vwgts, adjwgts, 
                          &wgtflag, &numflag, &split_size, 
                          options, &edgecut, elem_part);
    }

    // Un-mask the partition
    for ( int i = 0; i < num_elements; i++ ){
      partition[i] = elem_part[elem_mask[i]];
    }

    delete [] elem_part;
  }
  else {
    // If there is no split, just assign all elements to the
    // root processor
    for ( int k = 0; k < num_elements; k++ ){
      partition[k] = 0;
    }
  }

  // Free the masked element connectivity and element pointer
  delete [] elem_conn;
  delete [] elem_ptr;

  // Now, re-order the variables so that they are almost contiguous
  // over each processor 
  memset(owned_nodes, 0, split_size*sizeof(int));
  memset(owned_dep_nodes, 0, split_size*sizeof(int));
  memset(owned_elements, 0, split_size*sizeof(int));

  // Set up an array that tracks the mapping from the old
  // node ordering to the new node ordering. The array stores
  // the information such that the old node with index i is
  // stored at the index new_nodes[i].

  // First, treat new_nodes as an array of flags that indicate
  // whether this node has been counted yet.
  memset(new_nodes, 0, num_nodes*sizeof(int));
  memset(new_dep_nodes, 0, num_dependent_nodes*sizeof(int));

  // Keep a count that is the new global ordering of the nodes
  // within the finite-element mesh.
  int count = 0;
  for ( int j = 0; j < num_elements; j++ ){
    // The owner of element j
    int owner = partition[elem_mask[j]];
    owned_elements[owner]++;

    // Assign the un-assigned nodes associated with element j
    // to the owner of element j.
    for ( int i = elem_node_ptr[j]; i < elem_node_ptr[j+1]; i++ ){
      int node = elem_node_conn[i];
      if (node >= 0){
	if (!new_nodes[node]){
	  new_nodes[node] = 1;
	  owned_nodes[owner]++;
	  count++;
	}
      }
      else {
	node = -node-1;
	if (!new_dep_nodes[node]){
	  new_dep_nodes[node] = 1;
	  owned_dep_nodes[owner]++;
	  count++;
	}
      }
    }
  }

  // Reset the new_nodes array and prepare to assign the new
  // node numbers in a contiguous fashion.
  for ( int k = 0; k < num_nodes; k++ ){
    new_nodes[k] = -1;
  }
  for ( int k = 0; k < num_dependent_nodes; k++ ){
    new_dep_nodes[k] = -1;
  }

  // Find the offset to the ordering of the nodes for each partition
  // such that split_offset[i] is the first node number for
  // the i-th processor
  int *split_offset = new int[ split_size ];
  split_offset[0] = 0;
  for ( int k = 1; k < split_size; k++ ){
    split_offset[k] = split_offset[k-1] + owned_nodes[k-1];
  }

  int *split_dep_offset = new int[ split_size ];
  split_dep_offset[0] = 0;
  for ( int k = 1; k < split_size; k++ ){
    split_dep_offset[k] = split_dep_offset[k-1] + owned_dep_nodes[k-1];
  }

  // Now assign the new_nodes with the correct offsets. As a result,
  // the nodes will be ordered properly.
  for ( int j = 0; j < num_elements; j++ ){
    int owner = partition[j];

    for ( int i = elem_node_ptr[j]; i < elem_node_ptr[j+1]; i++ ){
      int node = elem_node_conn[i];
      if (node >= 0){
	if (new_nodes[node] < 0){
	  new_nodes[node] = split_offset[owner];
	  split_offset[owner]++;
	}
      }
      else {
	node = -node-1;
	if (new_dep_nodes[node] < 0){
	  new_dep_nodes[node] = split_dep_offset[owner];
	  split_dep_offset[owner]++;
	}
      }
    }
  }

  // Free the element mask
  delete [] elem_mask;

  // Free the split offset and split dependent node offset arrays
  delete [] split_offset;
  delete [] split_dep_offset;
}

