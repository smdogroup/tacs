#ifndef TACS_MESH_LOADER_H
#define TACS_MESH_LOADER_H

/*
  Copyright (c) 2011 Graeme Kennedy. All rights reserved.
  Not for commercial purposes.

  This class provides a limited capability to read in nodal and
  connectivity information from a .bdf file - and could be extended
  to read in other formats. The class can also be used to distribute
  the mesh over a set of processors. 

  This class could easily be extended to take in the mesh and
  connectivity from memory on the root processor and then distribute
  it to other processors.
  
  The limited capabilities of reading in data from a Nastran file are:
  1. Reading GRID and GRID* entries for the physical locations of the 
  nodes.
  2. Reading in connectivity information from CQUAD4 elements.
  3. Reading SPC entries (single point constraint) for the application
  of boundary conditions.

  The loader does not understand the different load-case capabilities
  that can be placed within a Nastran file. The elements must be passed 
  in to the object based on the component number.
*/

#include "mpi.h"
#include "TACSObject.h"
#include "TACSElement.h"
#include "TACSAssembler.h"
#include "TACSToFH5.h"

class TACSMeshLoader : public TACSObject {
 public:
  TACSMeshLoader( MPI_Comm _comm );
  ~TACSMeshLoader();

  // Read a BDF file for input
  // -------------------------
  int scanBDFFile( const char * file_name );

  // Get information about the mesh after scanning
  // ---------------------------------------------
  int getNumComponents();
  const char* getComponentDescript( int comp_num );
  const char* getElementDescript( int comp_num );

  // Set the elements corresponding to each of the component numbers 
  // ---------------------------------------------------------------
  void setElement( int component_num, TACSElement * _element );

  // Retrieve the element numbers corresponding to the given 
  // component numbers
  // -------------------------------------------------------
  int getNumNodes();
  int getNumElements();

  // Create a TACSToFH5 file writer
  // ------------------------------
  TACSToFH5* createTACSToFH5( TACSAssembler * tacs,
			      enum ElementType elem_type,
			      unsigned int write_flag );

  // Distribute the mesh and create TACS
  // -----------------------------------
  TACSAssembler* createTACS( int vars_per_node,
			     enum TACSAssembler::OrderingType order_type = 
			     TACSAssembler::NATURAL_ORDER, 
			     enum TACSAssembler::MatrixOrderingType mat_type = 
			     TACSAssembler::DIRECT_SCHUR);

  // Set the domain of a structural function with component numbers
  // --------------------------------------------------------------
  // void setFunctionDomain( TACSFunction * function, 
  //                         int comp_nums[], int num_comps );

  void getConnectivity( int *_num_nodes, int *_num_elements,
                        const int **_elem_node_ptr, 
                        const int **_elem_node_conn,
                        const double **_Xpts ){
    if (_num_nodes){ *_num_nodes = num_nodes; }
    if (_num_elements){ *_num_elements = num_elements; }
    if (_elem_node_ptr){ *_elem_node_ptr = elem_node_ptr; }
    if (_elem_node_conn){ *_elem_node_conn = elem_node_conn; }
    if (_Xpts){ *_Xpts = Xpts; }
  }

 private:
  // Communicator for all processors
  MPI_Comm comm;

  // The element corresponding to each of the component numbers
  TACSElement **elements;

  // Original BDF mesh information: Note that the original
  // ordering may not be contiguous. The node numbers associated
  // with the original ordering are stored in node_nums.
  int *node_nums;
  double *Xpts_unsorted; 

  // Reduced set of contiguous nodes 
  double *Xpts;

  // The mesh and element connectivity
  int num_nodes, num_elements;
  int *elem_node_conn, *elem_node_ptr;
  int *elem_component;

  // Store information about the components
  int num_components;
  char *component_elems;
  char *component_descript;

  // The boundary conditions
  int num_bcs;
  int *bc_nodes, *bc_vars, *bc_ptr;
  double *bc_vals;
};

#endif
