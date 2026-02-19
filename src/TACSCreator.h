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

#ifndef TACS_CREATOR_H
#define TACS_CREATOR_H

#include "TACSAssembler.h"

/**
  The following file contains the TACSCreator object which can be used
  to create TACSAssembler object for finite-element analysis.

  This code automatically redistributes a serial mesh in parallel. It cannot
  handle extremely large meshes, but is useful for many moderate-scale
  applications.

  The user must specify the connectivity, boundary conditions,
  elements, node locations and optionally any dependent nodes that are
  defined.  The elements may be provided in a list, or a callback
  function can be used to create them. In the case of the callback,
  the function is called repeatedly for every element within the
  finite-element mesh (not just once per element-id number).

  The user may wish to modify the ordering of TACS. This can be done
  by specifying the reordering type prior to calling createTACS().

  The new node numbers and new element partition can be retrieved from
  the creator object using the getNodeNums()/getElementPartion().
  Note that it is guaranteed that on each partiton, the elements will
  be numbered in ascending global order. This can be used to remap the
  distributed element order back to the original element order.
*/
class TACSCreator : public TACSObject {
 public:
  TACSCreator(MPI_Comm comm, int _vars_per_node);
  ~TACSCreator();

  // Set the connectivity for the global mesh
  // ----------------------------------------
  void setGlobalConnectivity(int _num_nodes, int _num_elements,
                             const int *_elem_node_ptr,
                             const int *_elem_node_conn,
                             const int *_elem_id_nums);

  // Set the boundary conditions
  // ---------------------------
  void setBoundaryConditions(int _num_bcs, const int *_bc_nodes,
                             const int *_bc_ptr = NULL,
                             const int *_bc_vars = NULL,
                             const TacsScalar *_bc_vals = NULL);

  // Set the dependent node connectivity and weights
  // -----------------------------------------------
  void setDependentNodes(int num_dep_nodes, const int *_dep_node_ptr,
                         const int *_dep_node_conn,
                         const double *_dep_node_weights);

  // Set the nodal locations
  // -----------------------
  void setNodes(const TacsScalar *_Xpts);

  // Set the type of ordering to use
  // -------------------------------
  void setReorderingType(TACSAssembler::OrderingType _order_type,
                         TACSAssembler::MatrixOrderingType _mat_type);

  // Partition the mesh
  // ------------------
  void partitionMesh(int split_size = 0, const int *part = NULL);

  // Set the elements into TACS creator
  // ----------------------------------
  void setElements(int _num_elems, TACSElement **_elements);
  void setElementCreator(TACSElement *(*func)(int, int));

  // Create the TACSAssembler object
  // -------------------------------
  TACSAssembler *createTACS();

  // Get local element numbers with the given set of element-id numbers
  // ------------------------------------------------------------------
  int getElementIdNums(int num_ids, int *ids, int **elem_nums);

  // Convert from the list of nodes from the original serial ordering
  // ----------------------------------------------------------------
  void getAssemblerNodeNums(TACSAssembler *tacs, int num_orig_nodes,
                            const int *orig_nodes, int *num_dist_nodes,
                            int **new_nodes);

  // Get the new node numbers and element partition on the root proc
  // ---------------------------------------------------------------
  int getNodeNums(const int **_new_nodes);
  int getElementPartition(const int **_partition);
  void getNumOwnedNodes(int **_owned_nodes);
  void getNumOwnedElements(int **_owned_elements);

 private:
  // The magic element-generator function pointer
  TACSElement *(*element_creator)(int local, int elem_id);

  // Set the type of reordering to use
  int use_reordering;
  TACSAssembler::OrderingType order_type;
  TACSAssembler::MatrixOrderingType mat_type;

  // The number of variables per node in the mesh
  int vars_per_node;

  // The MPI communicator
  MPI_Comm comm;
  int root_rank;

  // The global connectivity information
  int num_nodes, num_elements;

  // The dependent node data, connectivity and weights
  int num_dependent_nodes;
  int *dep_node_ptr, *dep_node_conn;
  double *dep_node_weights;

  // The element connectivity
  int *elem_node_ptr, *elem_node_conn;

  // Unique global element identifier
  int *elem_id_nums;

  // Boundary conditions
  int num_bcs;
  int *bc_nodes, *bc_vars, *bc_ptr;
  TacsScalar *bc_vals;

  // The node locations
  TacsScalar *Xpts;

  // Elements and the corresponding element id numbers
  int num_elem_ids;
  TACSElement **elements;

  // Keep the number of owned nodes/elements
  int *owned_nodes, *owned_elements;

  // The new node numbers for the independent nodes
  int *new_nodes;

  // The element partition
  int *partition;

  // Local information about the partitioned mesh
  int num_owned_elements, num_owned_nodes;
  int *local_elem_id_nums;
};

#endif  // TACS_CREATOR_H
