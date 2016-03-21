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
			      const int *_elem_node_ptr, 
                              const int *_elem_node_conn,
			      const int *_elem_id_nums );

  // Set the boundary conditions
  // ---------------------------
  void setBoundaryConditions( int _num_bcs, 
                              const int *_bc_nodes, 
			      const int *_bc_vars, 
                              const int *_bc_ptr );

  // Set the dependent node connectivity and weights
  // -----------------------------------------------
  void setDependentNodes( int num_dep_nodes, 
                          const int *_dep_node_ptr,
			  const int *_dep_node_conn, 
			  const TacsScalar *_dep_node_weights );

  // Set the elements into TACS creator
  // ----------------------------------
  void setElements( TACSElement **_elements, int _num_elem_ids );

  // Set the nodal locations
  // -----------------------
  void setNodes( const TacsScalar *_Xpts );

  // Create the TACSAssembler object
  // -------------------------------
  TACSAssembler *createTACS( enum TACSAssembler::OrderingType order_type =
			     TACSAssembler::ND_ORDER,
			     enum TACSAssembler::MatrixOrderingType mat_type =
			     TACSAssembler::DIRECT_SCHUR );

  // Get the new node numbers
  // ------------------------
  int getNodeNums( const int **_new_nodes ){ 
    *_new_nodes = new_nodes;
    return num_nodes; 
  }

 private:
  // Partition the mesh stored internally
  void splitSerialMesh( int split_size, int *partition, 
			int *owned_elements, int *owned_nodes );

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
  TacsScalar *dep_node_weights;

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

  // The new node numbers for the independent nodes
  int *new_nodes;
};

#endif // TACS_CREATOR_H
