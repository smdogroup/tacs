/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2014 Georgia Tech Research Corporation

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/

#include "TACSCreator.h"

#include "TacsUtilities.h"
#include "tacsmetis.h"

/*
  Extend an integer array
*/
static void extend_int_array(int **array, int old_len, int new_len) {
  int *temp = new int[new_len];
  memcpy(temp, *array, sizeof(int) * old_len);
  delete[] * array;
  *array = temp;
}

/*
  Functions for sorting a list such that once the list is sorted:

  arg_sort_list[list[i]] is in ascending order
*/
static const int *arg_sort_list = NULL;

static int compare_arg_sort(const void *a, const void *b) {
  const int aval = arg_sort_list[*(int *)a];
  const int bval = arg_sort_list[*(int *)b];

  // Break ties using the index of the list
  if (aval == bval) {
    return *(const int *)(a) - *(const int *)(b);
  }

  return aval - bval;
}

/**
  Allocate the TACSCreator object

  @param _comm The MPI_Comm object
  @param _vars_per_node The number of variables at each node
*/
TACSCreator::TACSCreator(MPI_Comm _comm, int _vars_per_node) {
  // Set the communicator
  comm = _comm;
  root_rank = 0;

  // Set the type of reordering to use. By default
  // this reordering is not applied to the mesh.
  use_reordering = 0;
  order_type = TACSAssembler::NATURAL_ORDER;
  mat_type = TACSAssembler::DIRECT_SCHUR;

  // Set the number of variables per node in the mesh
  vars_per_node = _vars_per_node;

  // Element connectivity information
  num_nodes = 0;
  num_elements = 0;
  num_dependent_nodes = 0;

  // Set the element connectivity and nodes
  elem_id_nums = NULL;
  elem_node_ptr = NULL;
  elem_node_conn = NULL;

  // Set the dependent node data
  dep_node_ptr = NULL;
  dep_node_conn = NULL;
  dep_node_weights = NULL;

  // Set the initial point locations
  Xpts = NULL;

  // Boundary condition information
  num_bcs = 0;
  bc_nodes = NULL;
  bc_ptr = NULL;
  bc_vars = NULL;
  bc_vals = NULL;

  // Set the dependent node pointers to zero
  new_nodes = NULL;

  // The element partition
  partition = NULL;
  owned_elements = NULL;
  owned_nodes = NULL;

  // Set the elements array to NULL
  elements = NULL;
  element_creator = NULL;

  // Information about the partitioned mesh
  num_owned_elements = 0;
  num_owned_nodes = 0;
  local_elem_id_nums = NULL;
}

/**
  Deallocate the TACSCreator object
*/
TACSCreator::~TACSCreator() {
  if (elem_id_nums) {
    delete[] elem_id_nums;
  }
  if (elem_node_ptr) {
    delete[] elem_node_ptr;
  }
  if (elem_node_conn) {
    delete[] elem_node_conn;
  }
  if (Xpts) {
    delete[] Xpts;
  }
  if (bc_nodes) {
    delete[] bc_nodes;
  }
  if (bc_vars) {
    delete[] bc_vars;
  }
  if (bc_vals) {
    delete[] bc_vals;
  }
  if (bc_ptr) {
    delete[] bc_ptr;
  }
  if (dep_node_ptr) {
    delete[] dep_node_ptr;
  }
  if (dep_node_conn) {
    delete[] dep_node_conn;
  }
  if (dep_node_weights) {
    delete[] dep_node_weights;
  }
  if (new_nodes) {
    delete[] new_nodes;
  }
  if (partition) {
    delete[] partition;
  }
  if (owned_elements) {
    delete[] owned_elements;
  }
  if (owned_nodes) {
    delete[] owned_nodes;
  }
  if (local_elem_id_nums) {
    delete[] local_elem_id_nums;
  }

  if (elements) {
    for (int i = 0; i < num_elem_ids; i++) {
      if (elements[i]) {
        elements[i]->decref();
      }
    }
    delete[] elements;
  }
}

/*
  Set the global element connectivity arrays into the class.

  @param _num_nodes The number of nodes in the mesh
  @param _num_elements The number of
*/
void TACSCreator::setGlobalConnectivity(int _num_nodes, int _num_elements,
                                        const int *_elem_node_ptr,
                                        const int *_elem_node_conn,
                                        const int *_elem_id_nums) {
  num_elements = _num_elements;
  num_nodes = _num_nodes;

  // Copy the pointer array into the elem_node_conn array
  elem_node_ptr = new int[num_elements + 1];
  memcpy(elem_node_ptr, _elem_node_ptr, (num_elements + 1) * sizeof(int));

  // Copy over the element node connectivity
  elem_node_conn = new int[elem_node_ptr[num_elements]];
  memcpy(elem_node_conn, _elem_node_conn,
         elem_node_ptr[num_elements] * sizeof(int));

  // Copy over the element ids
  elem_id_nums = new int[num_elements];
  memcpy(elem_id_nums, _elem_id_nums, num_elements * sizeof(int));
}

/*
  Set the dependent node information
*/
void TACSCreator::setDependentNodes(int num_dep_nodes, const int *_dep_node_ptr,
                                    const int *_dep_node_conn,
                                    const double *_dep_node_weights) {
  // Set the number of dependent nodes
  num_dependent_nodes = num_dep_nodes;

  // Allocate space for the dependent node pointer data
  dep_node_ptr = new int[num_dependent_nodes + 1];
  memcpy(dep_node_ptr, _dep_node_ptr, (num_dependent_nodes + 1) * sizeof(int));

  // Copy over the connectivity
  int conn_len = dep_node_ptr[num_dependent_nodes];
  dep_node_conn = new int[conn_len];
  dep_node_weights = new double[conn_len];
  memcpy(dep_node_conn, _dep_node_conn, conn_len * sizeof(int));
  memcpy(dep_node_weights, _dep_node_weights, conn_len * sizeof(double));
}

/*
  Set the boundary condition data
*/
void TACSCreator::setBoundaryConditions(int _num_bcs, const int *_bc_nodes,
                                        const int *_bc_ptr, const int *_bc_vars,
                                        const TacsScalar *_bc_vals) {
  // Set the number of boundary conditions and the node numbers
  // to which they correspond
  num_bcs = _num_bcs;
  bc_nodes = new int[num_bcs];
  memcpy(bc_nodes, _bc_nodes, num_bcs * sizeof(int));

  // Allocate space for the boundary condition pointer
  bc_ptr = new int[num_bcs + 1];

  if (_bc_ptr && _bc_vars) {
    // Copy over the boundary condition variable pointer
    memcpy(bc_ptr, _bc_ptr, (num_bcs + 1) * sizeof(int));

    // Allocate the number of variables
    bc_vars = new int[bc_ptr[num_bcs]];
    memcpy(bc_vars, _bc_vars, bc_ptr[num_bcs] * sizeof(int));

    // Allocate and set the variable values
    bc_vals = new TacsScalar[bc_ptr[num_bcs]];
    if (_bc_vals) {
      memcpy(bc_vals, _bc_vals, bc_ptr[num_bcs] * sizeof(TacsScalar));
    } else {
      memset(bc_vals, 0, bc_ptr[num_bcs] * sizeof(TacsScalar));
    }
  } else {
    // Since the bc_vars array is input as NULL, assume that
    // all the variables at this node are fully restrained.
    bc_vars = new int[num_bcs * vars_per_node];
    bc_ptr[0] = 0;
    for (int i = 0; i < num_bcs; i++) {
      bc_ptr[i + 1] = bc_ptr[i] + vars_per_node;

      // Restrain all the variables at the given node
      for (int j = 0; j < vars_per_node; j++) {
        bc_vars[bc_ptr[i] + j] = j;
      }
    }

    // Allocate the values array and set it equal to zero
    bc_vals = new TacsScalar[bc_ptr[num_bcs]];
    memset(bc_vals, 0, bc_ptr[num_bcs] * sizeof(TacsScalar));
  }
}

/*
  Set the elements with an array indexed by element id
*/
void TACSCreator::setElements(int _num_elem_ids, TACSElement **_elements) {
  num_elem_ids = _num_elem_ids;
  elements = new TACSElement *[num_elem_ids];
  memcpy(elements, _elements, num_elem_ids * sizeof(TACSElement *));
  for (int i = 0; i < num_elem_ids; i++) {
    if (elements[i]) {
      elements[i]->incref();
    }
  }
}

/*
  Set the element creator callback function
*/
void TACSCreator::setElementCreator(TACSElement *(*func)(int, int)) {
  element_creator = func;
}

/*
  Set the nodal locations
*/
void TACSCreator::setNodes(const TacsScalar *_Xpts) {
  Xpts = new TacsScalar[3 * num_nodes];
  memcpy(Xpts, _Xpts, 3 * num_nodes * sizeof(TacsScalar));
}

/*
  Set the type of ordering to use
*/
void TACSCreator::setReorderingType(
    TACSAssembler::OrderingType _order_type,
    TACSAssembler::MatrixOrderingType _mat_type) {
  use_reordering = 1;
  order_type = _order_type;
  mat_type = _mat_type;
}

/*
  Get the new node numbers
*/
int TACSCreator::getNodeNums(const int **_new_nodes) {
  *_new_nodes = new_nodes;
  return num_nodes;
}

/*
  Get the element partition
*/
int TACSCreator::getElementPartition(const int **_partition) {
  if (_partition) {
    *_partition = partition;
  }
  return num_elements;
}

/*
  Get the number of nodes owned by each processor
*/
void TACSCreator::getNumOwnedNodes(int **_owned_nodes) {
  *_owned_nodes = owned_nodes;
}

/*
  Get the number of elements owned by each processor
*/
void TACSCreator::getNumOwnedElements(int **_owned_elements) {
  *_owned_elements = owned_elements;
}

/*
  Input the node numbers on the root processor in the original order,
  and get out the distributed node numbers on the final mesh
*/
void TACSCreator::getAssemblerNodeNums(TACSAssembler *assembler,
                                       int num_orig_nodes,
                                       const int *_orig_nodes,
                                       int *num_dist_nodes, int **_tacs_nodes) {
  // Figure out how to distribute the nodes
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // The array of original nodes - only relevant on the root proc
  int *orig_nodes = NULL;
  int *ext_ptr = NULL;

  // The array of new node numbers
  int *tacs_nodes = NULL;

  // Get the variable map from TACSAssembler
  TACSNodeMap *nodeMap = assembler->getNodeMap();
  const int *owner_range = NULL;
  nodeMap->getOwnerRange(&owner_range);

  if (rank == root_rank) {
    // First allocate an array of nodes that we will overwrite
    orig_nodes = new int[num_orig_nodes];
    memcpy(orig_nodes, _orig_nodes, num_orig_nodes * sizeof(int));

    // Convert the original node ordering to the new node ordering
    for (int i = 0; i < num_orig_nodes; i++) {
      if (orig_nodes[i] >= 0 && orig_nodes[i] < num_nodes) {
        orig_nodes[i] = new_nodes[orig_nodes[i]];
      } else {
        orig_nodes[i] = 0;
      }
    }

    // Sort the nodes uniquely
    num_orig_nodes = TacsUniqueSort(num_orig_nodes, orig_nodes);

    // Match the intervals for the external node numbers
    ext_ptr = new int[size + 1];
    TacsMatchIntervals(size, owner_range, num_orig_nodes, orig_nodes, ext_ptr);

    // Send the processors the information, one at a time
    for (int i = 0; i < size; i++) {
      if (i != root_rank) {
        int count = ext_ptr[i + 1] - ext_ptr[i];
        int tag = 1;
        MPI_Send(&orig_nodes[ext_ptr[i]], count, MPI_INT, i, tag, comm);
      }
    }

    // Get the number of distributed nodes for the root proc
    int count = ext_ptr[root_rank + 1] - ext_ptr[root_rank];
    *num_dist_nodes = count;

    // Copy the nodes on the root proc
    tacs_nodes = new int[count];
    memcpy(tacs_nodes, &orig_nodes[ext_ptr[root_rank]], count * sizeof(int));
  } else {
    MPI_Status status;
    int tag = 1;
    MPI_Probe(root_rank, tag, comm, &status);

    // Get the count that we are about to recv
    MPI_Get_count(&status, MPI_INT, num_dist_nodes);
    int count = *num_dist_nodes;
    if (count > 0) {
      tacs_nodes = new int[count];
    }

    // Recv the nodes from the root
    MPI_Recv(tacs_nodes, count, MPI_INT, root_rank, tag, comm,
             MPI_STATUS_IGNORE);
  }

  // Apply the reordering in TACS
  assembler->reorderNodes(*num_dist_nodes, tacs_nodes);

  *_tacs_nodes = tacs_nodes;
}

/*
  Create the instance of the TACSAssembler object and return it.

  This code partitions the mesh, calls for the elements to be
  allocated and returns a valid instance of the TACSAssembler object.
*/
TACSAssembler *TACSCreator::createTACS() {
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  if (rank == root_rank && !partition) {
    // Partition the mesh using the serial code on the root
    // processor.
    partitionMesh(size);
  }

  // The number of local and owned nodes and the number of
  // elements for each processor in the mesh
  num_owned_nodes = 0;
  num_owned_elements = 0;
  MPI_Scatter(owned_nodes, 1, MPI_INT, &num_owned_nodes, 1, MPI_INT, root_rank,
              comm);
  MPI_Scatter(owned_elements, 1, MPI_INT, &num_owned_elements, 1, MPI_INT,
              root_rank, comm);

  // The number of local nodes
  int num_local_dep_nodes = 0;

  // Allocate space for the portion of the element connectivity on
  // this processor
  int *local_elem_node_ptr = new int[num_owned_elements + 1];
  int *local_elem_node_conn = NULL;

  // This will be used later to determine which elements belong to
  // which domain within the finite-element mesh
  local_elem_id_nums = new int[num_owned_elements];

  // Loacal nodal information
  TacsScalar *Xpts_local = NULL;

  // Local dependent node information
  int *local_dep_node_ptr = NULL;
  int *local_dep_node_conn = NULL;
  double *local_dep_node_weights = NULL;

  // For each processor, send the information to the owner
  if (rank == root_rank) {
    // Reset the nodes for the boundary conditions so that they
    // correspond to the new ordering
    for (int j = 0; j < num_bcs; j++) {
      bc_nodes[j] = new_nodes[bc_nodes[j]];
    }

    // Find the inverse mapping between the new and old node
    // numbers so that it's faster to access
    int *inv_new_nodes = new int[num_nodes];
    for (int i = 0; i < num_nodes; i++) {
      inv_new_nodes[new_nodes[i]] = i;
    }

    // Set up the element partition so that it is sorted so that
    // we can quickly find the elements associated with a processor
    int *elem_part = new int[num_elements];
    for (int k = 0; k < num_elements; k++) {
      elem_part[k] = k;
    }

    // Now partition[elem_part[k]] is sorted
    arg_sort_list = partition;
    qsort(elem_part, num_elements, sizeof(int), compare_arg_sort);
    arg_sort_list = NULL;

    // Determine a sorted ordering of the new dependent nodes
    // that can be used to map depedent nodes to process owners
    int *dep_nodes_part = NULL;
    if (num_dependent_nodes > 0) {
      dep_nodes_part = new int[num_dependent_nodes];
      for (int k = 0; k < num_dependent_nodes; k++) {
        dep_nodes_part[k] = k;
      }
    }

    // Compute the local CSR data structure on this process
    // Use an upper bound for the memory requirements
    int max_conn_size = elem_node_ptr[num_elements];
    int *dep_nodes = new int[num_dependent_nodes];
    int *dep_node_flags = new int[num_dependent_nodes];
    int *inv_dep_nodes = new int[num_dependent_nodes];

    // The local element connectivity
    int *elem_ptr = new int[num_elements + 1];
    int *elem_conn = new int[max_conn_size];
    int *elem_ids = new int[num_elements];
    TacsScalar *xpts = new TacsScalar[3 * num_nodes];

    // Allocate space for the locally-reduced dependent node data
    int *dep_ptr = NULL;
    int *dep_conn = NULL;
    double *dep_weights = NULL;
    if (num_dependent_nodes > 0) {
      dep_ptr = new int[num_dependent_nodes + 1];
      dep_conn = new int[dep_node_ptr[num_dependent_nodes]];
      dep_weights = new double[dep_node_ptr[num_dependent_nodes]];
    }

    // Set the dependent node flags to -1
    for (int i = 0; i < num_dependent_nodes; i++) {
      dep_node_flags[i] = -1;
    }

    // Loop over all the processors on the root processor and
    // create the local connectivity, element info, and dependent node
    // data for the k-th processor. Then send it to processor k, except
    // in the case when k == root_rank in which case make a local
    // copy of the data without sending anything.
    int start = 0;
    int start_node = 0;
    for (int k = 0; k < size; k++) {
      // Keep track of the starting and ending element
      // index within the sorted list
      int end = start + owned_elements[k];
      int end_node = start_node + owned_nodes[k];

      // Keep track of the number of local nodes and local depdendent
      // nodes that are indexed by processor k. This will be larger
      // than owned_nodes[k].
      int local_dep_node_size = 0;

      // Set the first index within the element pointer to zero
      elem_ptr[0] = 0;

      // Cycle through the element partition k, and assign
      // the connectivity for the new
      int n = 0;
      int local_conn_size = 0;
      for (int j = start; j < end; j++, n++) {
        int elem = elem_part[j];

        // Add the connectivity from the element 'elem'
        // that belongs to partition k
        int iend = elem_node_ptr[elem + 1];
        for (int i = elem_node_ptr[elem]; i < iend; i++) {
          // Get add the node indices associated with the element
          int node = elem_node_conn[i];
          if (node >= 0) {
            // Add the global connectivity index - this will
            // be converted to a local index later.
            elem_conn[local_conn_size] = new_nodes[node];
            local_conn_size++;
          } else {
            // Convert to the dependent node index
            node = -node - 1;

            // If this dependent node has not been assigned,
            // create a new number for it
            if (dep_node_flags[node] != k) {
              dep_node_flags[node] = k;
              dep_nodes[node] = local_dep_node_size;
              inv_dep_nodes[local_dep_node_size] = node;
              local_dep_node_size++;
            }

            // Assign the dependent node number
            elem_conn[local_conn_size] = -dep_nodes[node] - 1;
            local_conn_size++;
          }
        }

        // Set the element
        elem_ids[n] = elem_id_nums[elem];
        elem_ptr[n + 1] = local_conn_size;
      }

      // Loop over all of the local dependent nodes and find their
      // corresponding independent nodes and copy over the weight
      // data.
      if (local_dep_node_size > 0) {
        // Keep track of the offset into the pointer array
        int dep_conn_len = 0;
        dep_ptr[0] = 0;

        for (int j = 0; j < local_dep_node_size; j++) {
          int dep_node = inv_dep_nodes[j];

          // Loop over the dependent nodes within this list
          for (int i = dep_node_ptr[dep_node]; i < dep_node_ptr[dep_node + 1];
               i++) {
            dep_conn[dep_conn_len] = new_nodes[dep_node_conn[i]];
            dep_weights[dep_conn_len] = dep_node_weights[i];
            dep_conn_len++;
          }

          // Set the depdendent node pointer data
          dep_ptr[j + 1] = dep_conn_len;
        }
      }

      // Copy over the node data from the old indices
      for (int i = 0, j = start_node; j < end_node; i++, j++) {
        int node = inv_new_nodes[j];
        xpts[3 * i] = Xpts[3 * node];
        xpts[3 * i + 1] = Xpts[3 * node + 1];
        xpts[3 * i + 2] = Xpts[3 * node + 2];
      }

      // Create the CSR data structure required
      if (k == root_rank) {
        // Copy the values over from the nodes
        num_local_dep_nodes = local_dep_node_size;
        Xpts_local = new TacsScalar[3 * num_owned_nodes];
        memcpy(Xpts_local, xpts, 3 * num_owned_nodes * sizeof(TacsScalar));

        if (num_local_dep_nodes > 0) {
          // Copy over the data for the dependent nodes
          local_dep_node_ptr = new int[num_local_dep_nodes + 1];
          memcpy(local_dep_node_ptr, dep_ptr,
                 (num_local_dep_nodes + 1) * sizeof(int));

          // Allocate the local dependent node connectivity and weights
          int conn_len = local_dep_node_ptr[num_local_dep_nodes];
          local_dep_node_conn = new int[conn_len];
          local_dep_node_weights = new double[conn_len];
          memcpy(local_dep_node_conn, dep_conn, conn_len * sizeof(int));
          memcpy(local_dep_node_weights, dep_weights,
                 conn_len * sizeof(double));
        }

        // Copy over values for the elements
        memcpy(local_elem_id_nums, elem_ids, num_owned_elements * sizeof(int));
        memcpy(local_elem_node_ptr, elem_ptr,
               (num_owned_elements + 1) * sizeof(int));
        local_elem_node_conn = new int[local_elem_node_ptr[num_owned_elements]];
        memcpy(local_elem_node_conn, elem_conn,
               local_elem_node_ptr[num_owned_elements] * sizeof(int));
      } else {
        // Send the data to the other process
        MPI_Send(&local_dep_node_size, 1, MPI_INT, k, 2, comm);
        MPI_Send(xpts, 3 * owned_nodes[k], TACS_MPI_TYPE, k, 4, comm);

        // Send the dependent node data to processor k
        if (local_dep_node_size > 0) {
          MPI_Send(dep_ptr, local_dep_node_size + 1, MPI_INT, k, 5, comm);
          MPI_Send(dep_conn, dep_ptr[local_dep_node_size], MPI_INT, k, 6, comm);
          MPI_Send(dep_weights, dep_ptr[local_dep_node_size], MPI_DOUBLE, k, 7,
                   comm);
        }

        // Send the element data
        MPI_Send(elem_ids, owned_elements[k], MPI_INT, k, 8, comm);
        MPI_Send(elem_ptr, owned_elements[k] + 1, MPI_INT, k, 9, comm);
        MPI_Send(elem_conn, elem_ptr[owned_elements[k]], MPI_INT, k, 10, comm);
      }

      // Increment the pointers to the start of the new nodes
      start += owned_elements[k];
      start_node += owned_nodes[k];
    }

    delete[] inv_new_nodes;
    delete[] elem_part;
    delete[] dep_nodes;
    delete[] dep_node_flags;
    delete[] inv_dep_nodes;
    delete[] elem_ptr;
    delete[] elem_conn;
    delete[] elem_ids;
    delete[] xpts;
    if (dep_nodes_part) {
      delete[] dep_nodes_part;
    }
    if (dep_ptr) {
      delete[] dep_ptr;
    }
    if (dep_conn) {
      delete[] dep_conn;
    }
    if (dep_weights) {
      delete[] dep_weights;
    }
  } else {
    // Recv the data from the root process
    MPI_Status status;
    MPI_Recv(&num_local_dep_nodes, 1, MPI_INT, root_rank, 2, comm, &status);

    // Allocate space for the incoming data
    Xpts_local = new TacsScalar[3 * num_owned_nodes];
    MPI_Recv(Xpts_local, 3 * num_owned_nodes, TACS_MPI_TYPE, root_rank, 4, comm,
             &status);

    // Receive the dependent node data if there are any dependent
    // nodes on this processor
    if (num_local_dep_nodes > 0) {
      local_dep_node_ptr = new int[num_local_dep_nodes + 1];
      MPI_Recv(local_dep_node_ptr, num_local_dep_nodes + 1, MPI_INT, root_rank,
               5, comm, &status);

      // Allocate the local connectivity and weight arrays
      int conn_len = local_dep_node_ptr[num_local_dep_nodes];
      local_dep_node_conn = new int[conn_len];
      local_dep_node_weights = new double[conn_len];
      MPI_Recv(local_dep_node_conn, conn_len, MPI_INT, root_rank, 6, comm,
               &status);
      MPI_Recv(local_dep_node_weights, conn_len, MPI_DOUBLE, root_rank, 7, comm,
               &status);
    }

    // Receive the element data
    MPI_Recv(local_elem_id_nums, num_owned_elements, MPI_INT, root_rank, 8,
             comm, &status);
    MPI_Recv(local_elem_node_ptr, num_owned_elements + 1, MPI_INT, root_rank, 9,
             comm, &status);

    int con_size = local_elem_node_ptr[num_owned_elements];
    local_elem_node_conn = new int[con_size];
    MPI_Recv(local_elem_node_conn, con_size, MPI_INT, root_rank, 10, comm,
             &status);
  }

  // Broadcast the boundary condition information
  MPI_Bcast(&num_bcs, 1, MPI_INT, root_rank, comm);

  if (num_bcs > 0) {
    // Broadcast all the boundary conditions to all processors
    if (rank != root_rank) {
      if (bc_nodes) {
        delete[] bc_nodes;
      }
      bc_nodes = new int[num_bcs];
    }
    MPI_Bcast(bc_nodes, num_bcs, MPI_INT, root_rank, comm);

    // Broacast the boundary condition pointer array
    if (rank != root_rank) {
      if (bc_ptr) {
        delete[] bc_ptr;
      }
      bc_ptr = new int[num_bcs + 1];
    }
    MPI_Bcast(bc_ptr, num_bcs + 1, MPI_INT, root_rank, comm);

    if (rank != root_rank) {
      if (bc_vars) {
        delete[] bc_vars;
      }
      bc_vars = new int[bc_ptr[num_bcs]];
    }
    MPI_Bcast(bc_vars, bc_ptr[num_bcs], MPI_INT, root_rank, comm);

    if (rank != root_rank) {
      if (bc_vals) {
        delete[] bc_vals;
      }
      bc_vals = new TacsScalar[bc_ptr[num_bcs]];
    }
    MPI_Bcast(bc_vals, bc_ptr[num_bcs], TACS_MPI_TYPE, root_rank, comm);
  }

  TACSAssembler *tacs =
      new TACSAssembler(comm, vars_per_node, num_owned_nodes,
                        num_owned_elements, num_local_dep_nodes);

  // Set the dependent node data
  if (num_local_dep_nodes > 0) {
    tacs->setDependentNodes(local_dep_node_ptr, local_dep_node_conn,
                            local_dep_node_weights);
  }

  // Set the connectivity
  tacs->setElementConnectivity(local_elem_node_ptr, local_elem_node_conn);

  // Add the elements
  if (elements) {
    TACSElement **elems = new TACSElement *[num_owned_elements];
    for (int k = 0; k < num_owned_elements; k++) {
      elems[k] = elements[local_elem_id_nums[k]];
      if (!elems[k]) {
        fprintf(stderr,
                "[%d] TACSCreator: Element undefined for element ID %d\n", rank,
                local_elem_id_nums[k]);
        MPI_Abort(comm, 1);
        return NULL;
      }
    }

    // Set the elements and free the array
    tacs->setElements(elems);
    delete[] elems;
  } else if (element_creator) {
    TACSElement **elems = new TACSElement *[num_owned_elements];
    for (int k = 0; k < num_owned_elements; k++) {
      elems[k] = element_creator(k, local_elem_id_nums[k]);
      if (!elems[k]) {
        fprintf(stderr, "[%d] TACSCreator: Callback failed for element ID %d\n",
                rank, local_elem_id_nums[k]);
        MPI_Abort(comm, 1);
        return NULL;
      }
    }

    // Set the elements and free the array
    tacs->setElements(elems);
    delete[] elems;
  } else {
    fprintf(stderr, "[%d] TACSCreator: Elements and callback not defined\n",
            rank);
    MPI_Abort(comm, 1);
    return NULL;
  }

  // Allocate the arrays to store the variable values
  int *bvars = new int[vars_per_node];
  TacsScalar *bvals = new TacsScalar[vars_per_node];

  // Set the boundary conditions
  for (int k = 0; k < num_bcs; k++) {
    if (bc_nodes[k] >= 0) {
      int nbcs = bc_ptr[k + 1] - bc_ptr[k];
      int n = 0;
      for (int j = 0; j < nbcs; j++) {
        if (bc_vars[bc_ptr[k] + j] < vars_per_node) {
          bvars[n] = bc_vars[bc_ptr[k] + j];
          bvals[n] = bc_vals[bc_ptr[k] + j];
          n++;
        }
      }
      if (n > 0) {
        tacs->addBCs(1, &bc_nodes[k], n, bvars, bvals);
      }
    }
  }

  // Free the bvars/bvals arrays
  delete[] bvars;
  delete[] bvals;

  // We no longer need the boundary condition information
  delete[] bc_nodes;
  bc_nodes = NULL;
  delete[] bc_ptr;
  bc_ptr = NULL;
  delete[] bc_vars;
  bc_vars = NULL;
  delete[] bc_vals;
  bc_vals = NULL;

  // Use the reordering if the flag has been set in the
  // TACSCreator object
  if (use_reordering) {
    tacs->computeReordering(order_type, mat_type);
  }

  // Finish the initialization of TACS
  tacs->initialize();

  // Create the new node vector
  TACSBVec *X = tacs->createNodeVec();
  X->incref();

  // Copy the node locations to the vector
  TacsScalar *Xpt_vals;
  X->getArray(&Xpt_vals);
  memcpy(Xpt_vals, Xpts_local, 3 * num_owned_nodes * sizeof(TacsScalar));

  // Reorder the node vector
  if (use_reordering) {
    tacs->reorderVec(X);
  }

  // Set the node locations
  tacs->setNodes(X);
  X->decref();

  // Free all the remaining memory
  delete[] local_elem_node_ptr;
  delete[] local_elem_node_conn;
  delete[] Xpts_local;

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
  part:            (optional) the specified partition
*/
void TACSCreator::partitionMesh(int split_size, const int *part) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank != root_rank) {
    return;
  }

  // Deallocate data that may have been allocated on a
  // previous call
  if (new_nodes) {
    delete[] new_nodes;
  }
  if (owned_elements) {
    delete[] owned_elements;
  }
  if (owned_nodes) {
    delete[] owned_nodes;
  }

  // Allocate the arrays to store the number of owned elements and
  // nodes
  int mpi_size;
  MPI_Comm_size(comm, &mpi_size);

  // Set the split size to ensure that it is less than the
  // comm size
  if (split_size <= 0 || split_size > mpi_size) {
    split_size = mpi_size;
  }

  if (part) {
    // If the partition is specified, then delete the old
    // one if it is defined
    if (partition) {
      delete[] partition;
    }
    partition = NULL;

    // Check whether the suggested partition is legitimate
    int legit = 1;
    for (int i = 0; i < num_elements; i++) {
      if (part[i] < 0 || part[i] >= split_size) {
        legit = 0;
        break;
      }
    }

    // If the partition is legit, copy it over
    if (legit) {
      partition = new int[num_elements];
      memcpy(partition, part, num_elements * sizeof(int));
    }
  }

  if (!partition) {
    // Allocate space for the new node numbers and the new dependent
    // node numbers
    partition = new int[num_elements];

    // Compute the node to element CSR data structure for both
    // the indepedent and dependent nodes
    int *node_elem_ptr = new int[num_nodes + 1];
    memset(node_elem_ptr, 0, (num_nodes + 1) * sizeof(int));

    for (int i = 0; i < num_elements; i++) {
      int end = elem_node_ptr[i + 1];
      for (int j = elem_node_ptr[i]; j < end; j++) {
        int node = elem_node_conn[j];
        if (node >= 0) {
          // This is an independent node
          node_elem_ptr[node + 1]++;
        }
      }
    }

    // Determine the size of the node to element array
    for (int i = 0; i < num_nodes; i++) {
      node_elem_ptr[i + 1] += node_elem_ptr[i];
    }
    int *node_elem_conn = new int[node_elem_ptr[num_nodes]];

    // Fill in the entries of the node->element data structure
    for (int i = 0; i < num_elements; i++) {
      int end = elem_node_ptr[i + 1];
      for (int j = elem_node_ptr[i]; j < end; j++) {
        int node = elem_node_conn[j];
        if (node >= 0) {
          // This is an independent node
          node_elem_conn[node_elem_ptr[node]] = i;
          node_elem_ptr[node]++;
        }
      }
    }

    // Reset the node_elem_ptr array to the correct values
    for (int i = num_nodes; i > 0; i--) {
      node_elem_ptr[i] = node_elem_ptr[i - 1];
    }
    node_elem_ptr[0] = 0;

    // Now we set up the element to element connectivity using the
    // set of maksed elements. For this to work within METIS, we have
    // to remove the diagonal contribution so that there is no
    // self-reference in the graph.
    int *elem_ptr = new int[num_elements + 1];
    elem_ptr[0] = 0;

    // Information to keep track of how big the data structure is
    int elem_conn_size = 0;

    // Estimate the maximum size of the connectivity data
    const int ROW_SIZE_EST = 27;
    int max_elem_conn_size = ROW_SIZE_EST * num_elements;
    int *elem_conn = new int[max_elem_conn_size];

    // Assemble things one row at a time
    int *row = new int[num_elements + 1];

    for (int i = 0; i < num_elements; i++) {
      // Set the size of the new row in the data structure to zero
      int row_size = 0;

      // From the element number, find the independent nodes
      // that are associated with this element number
      int jpend = elem_node_ptr[i + 1];
      for (int jp = elem_node_ptr[i]; jp < jpend; jp++) {
        int node = elem_node_conn[jp];

        if (node >= 0) {
          // Add the element->element data
          int kpend = node_elem_ptr[node + 1];
          for (int kp = node_elem_ptr[node]; kp < kpend; kp++) {
            // Get the index of the new element
            int new_elem = node_elem_conn[kp];

            // Check if adding another element will exceed the size
            // of the array. If it will, sort the array and remove
            // duplicates. This is guaranteed to make additional
            // room.
            if (row_size >= num_elements) {
              row_size = TacsUniqueSort(row_size, row);
            }

            // Append the element index to the list
            row[row_size] = new_elem;
            row_size++;
          }
        }
      }

      // Sort the element indices and remove duplicates from the
      // array
      row_size = TacsUniqueSort(row_size, row);

      // Check if adding this new row to the data structure will excee
      if (elem_conn_size + row_size > max_elem_conn_size) {
        max_elem_conn_size += 0.5 * max_elem_conn_size;
        if (max_elem_conn_size < elem_conn_size + row_size) {
          max_elem_conn_size += elem_conn_size + row_size;
        }
        extend_int_array(&elem_conn, elem_conn_size, max_elem_conn_size);
      }

      // Add the elements - minus the diagonal entry
      for (int j = 0; j < row_size; j++) {
        if (row[j] != i) {  // Not the diagonal
          elem_conn[elem_conn_size] = row[j];
          elem_conn_size++;
        }
      }
      elem_ptr[i + 1] = elem_conn_size;
    }

    // Free the memory for data that is not needed
    delete[] row;

    // Free the node to element pointer
    delete[] node_elem_ptr;
    delete[] node_elem_conn;

    // Partition the mesh using METIS.
    if (split_size > 1) {
      int ncon = 1;  // "It should be at least 1"??

      // Set the default options
      int options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(options);

      // Use 0-based numbering
      options[METIS_OPTION_NUMBERING] = 0;

      // The objective value in METIS
      int objval = 0;

      if (split_size < 8) {
        METIS_PartGraphRecursive(&num_elements, &ncon, elem_ptr, elem_conn,
                                 NULL, NULL, NULL, &split_size, NULL, NULL,
                                 options, &objval, partition);
      } else {
        METIS_PartGraphKway(&num_elements, &ncon, elem_ptr, elem_conn, NULL,
                            NULL, NULL, &split_size, NULL, NULL, options,
                            &objval, partition);
      }
    } else {
      // If there is no split, just assign all elements to the
      // root processor
      for (int k = 0; k < num_elements; k++) {
        partition[k] = 0;
      }
    }

    // Free the element connectivity and element pointer
    delete[] elem_conn;
    delete[] elem_ptr;
  }

  // Now, re-order the node numbers so that they are almost contiguous
  // over each processor
  new_nodes = new int[num_nodes];
  owned_elements = new int[mpi_size];
  owned_nodes = new int[mpi_size];

  // Set the number of nodes/elements to zero
  memset(new_nodes, 0, num_nodes * sizeof(int));
  memset(owned_elements, 0, mpi_size * sizeof(int));
  memset(owned_nodes, 0, mpi_size * sizeof(int));

  // Set up an array that tracks the mapping from the old
  // node ordering to the new node ordering. The array stores
  // the information such that the old node with index i is
  // stored at the index new_nodes[i].

  // First, treat new_nodes as an array of flags that indicate
  // whether this node has been counted yet. Keep track of the new
  // global ordering of the nodes within the finite-element mesh.
  for (int j = 0; j < num_elements; j++) {
    // The owner of element j
    int owner = partition[j];
    owned_elements[owner]++;

    // Assign the un-assigned nodes associated with element j
    // to the owner of element j.
    for (int i = elem_node_ptr[j]; i < elem_node_ptr[j + 1]; i++) {
      int node = elem_node_conn[i];
      if (node >= 0) {
        if (!new_nodes[node]) {
          new_nodes[node] = 1;
          owned_nodes[owner]++;
        }
      }
    }
  }

  // Reset the new_nodes array and prepare to assign the new
  // node numbers in a contiguous fashion.
  for (int k = 0; k < num_nodes; k++) {
    new_nodes[k] = -1;
  }

  // Find the offset to the ordering of the nodes for each partition
  // such that split_offset[i] is the first node number for
  // the i-th processor
  int *split_offset = new int[split_size];
  split_offset[0] = 0;
  for (int k = 1; k < split_size; k++) {
    split_offset[k] = split_offset[k - 1] + owned_nodes[k - 1];
  }

  // Now assign the new_nodes with the correct offsets. As a result,
  // the nodes will be ordered properly.
  for (int j = 0; j < num_elements; j++) {
    int owner = partition[j];

    for (int i = elem_node_ptr[j]; i < elem_node_ptr[j + 1]; i++) {
      int node = elem_node_conn[i];
      if (node >= 0) {
        if (new_nodes[node] < 0) {
          new_nodes[node] = split_offset[owner];
          split_offset[owner]++;
        }
      }
    }
  }

  // Free the split offset and split dependent node offset arrays
  delete[] split_offset;
}

/*
  Retrieve the element numbers on each processor corresponding to the
  given component numbers.
*/
int TACSCreator::getElementIdNums(int num_ids, int ids[], int **elem_nums) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (!local_elem_id_nums) {
    fprintf(stderr,
            "[%d] TACSCreator: Cannot get elements until "
            "mesh is partitioned\n",
            rank);
    *elem_nums = NULL;
    return 0;
  }

  // Sort the component numbers on input
  num_ids = TacsUniqueSort(num_ids, ids);

  int *all_elems = new int[num_owned_elements];
  int elem_size = 0;
  for (int k = 0; k < num_owned_elements; k++) {
    if (TacsSearchArray(local_elem_id_nums[k], num_ids, ids)) {
      all_elems[elem_size] = k;
      elem_size++;
    }
  }

  // Truncate the array to the correct size
  *elem_nums = new int[elem_size];
  memcpy(*elem_nums, all_elems, elem_size * sizeof(int));
  delete[] all_elems;

  return elem_size;
}
