#ifndef TACS_CREATOR_H
#define TACS_CREATOR_H

/*
  The following file contains the TACSCreator object which can be used
  to create TACSAssembler object for finite-element analysis.

*/

#include "TACSAssembler.h"

class TACSCreator : public TACSObject {
 public:
  TACSCreator( MPI_Comm comm, int _vars_per_node );
  ~TACSCreator();

  // Set the connectivity for the global mesh
  // ----------------------------------------
  void setGlobalConnectivity( int _num_nodes, int _num_elements,
			      int *_elem_node_ptr, int *_elem_node_conn,
			      int *_elem_id_nums );

  // Set the boundary conditions
  // ---------------------------
  void setBoundaryConditions( int _num_bcs, int *_bc_nodes, int *_bc_vars,
			      int *_bc_ptr );

  // Set the elements into TACS creator
  // ----------------------------------
  void setElements( TACSElement **_elements, int _num_elem_ids );

  // Set the nodal locations
  // -----------------------
  void setNodes( TacsScalar *_Xpts );

  // Create the TACSAssembler object
  // -------------------------------
  TACSAssembler *createTACS( enum TACSAssembler::OrderingType order_type =
			     TACSAssembler::ND_ORDER,
			     enum TACSAssembler::MatrixOrderingType mat_type =
			     TACSAssembler::DIRECT_SCHUR );

 private:
  // Partition the mesh stored internally
  void splitSerialMesh( int split_size, int *elem_partition, 
			int *new_nodes, int *owned_elements, 
			int *owned_nodes );

  // The number of variables per node in the mesh
  int vars_per_node;

  // The MPI communicator
  MPI_Comm comm;
  int root_rank;

  // The global connectivity information
  int num_nodes, num_elements;

  // The element connectivity
  int *elem_node_ptr, *elem_node_conn;
  
  // Unique global element identifier
  int *elem_id_nums; 

  // Boundary conditions
  int num_bcs;
  int *bc_nodes, *bc_vars, *bc_ptr;

  // The node locations
  TacsScalar *Xpts;

  // Elements corresponding to the 
  int num_elem_ids;
  TACSElement **elements;
};

#endif // TACS_CREATOR_H
