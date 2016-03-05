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
  int scanBdfFile( const char * file_name );

  // Get information about the mesh after scanning
  // ---------------------------------------------
  int getNumComponents();
  const char * getComponentDescript( int comp_num );
  const char * getElementDescript( int comp_num );

  // Set the elements corresponding to each of the component numbers 
  // ---------------------------------------------------------------
  void setElement( int component_num, TACSElement * _element );

  // Retrieve the element numbers corresponding to the given 
  // component numbers
  // -------------------------------------------------------
  int getNumNodes();
  int getNumElements();
  void getComponentNums( int comp_nums[], int num_elements );
  int getComponentElementNums( int ** elem_nums,
                               int comp_nums[], int num_comps );
  void setFunctionDomain( TACSFunction * function, 
			  int comp_nums[], int num_comps );

  // Create a TACSToFH5 file writer
  // ------------------------------
  TACSToFH5 * createTACSToFH5( TACSAssembler * tacs,
                               enum ElementType elem_type,
                               unsigned int write_flag );

  // Distribute the mesh and create TACS
  // -----------------------------------
  TACSAssembler * createTACS( int vars_per_node,
			      int num_load_cases,
			      enum TACSAssembler::OrderingType order_type = 
			      TACSAssembler::NO_ORDER, 
			      enum TACSAssembler::MatrixOrderingType mat_type = 
			      TACSAssembler::DIRECT_SCHUR);
  
  // Just create TACS - do not distribute the mesh
  // ---------------------------------------------
  TACSAssembler * createSerialTACS( int split_size,
                                    int vars_per_node,
                                    int num_load_cases );

  // Add an element traction class to all elements in the given component
  // --------------------------------------------------------------------
  void addTraction( TACSSurfaceTraction * st,
                    int component_num, TACSElementTraction * traction );

  // Functions required for writing a BDF file
  // -----------------------------------------
  void getNumElementsForComps( int * numElem, int sizeNumComp );
  void getConnectivity( int * conn, int sizeConn );
  void getElementComponents( int * compIDs, int sizeCompIds );
  int getTotalNumElements();
  void getConnectivityForComp( int compID, int * conn, int sizeConn );
  void getNodes( int * nodeList, int nNodes, double * pts, int nPts );
  void getOrigNodes( double * xOrig, int n );
  void getOrigNodeNums( int * nodeNumsOrig, int n );
  void writeBDF( const char * filename, TacsScalar *bdfNodes, 
		 int * nodeNums, int nBDFNodes );

  // Write the data from the constitutive objects
  // --------------------------------------------
  void writeBDFConstitutive( const char *filename );

 private:
  // Split the mesh and reorder (but do not distribute)
  void splitMesh( int split_size, 
                  int * elem_partition, int * new_nodes,
                  int * num_owned_elements, int * num_owned_nodes );

  // Write the constitutive data to the given file pointer
  void writeBDFConstitutive( FILE *fp );

  MPI_Comm comm;

  // The element corresponding to each of the component numbers
  TACSElement ** elements;

  // Original BDF mesh information:
  int * node_nums;
  double * Xpts_unsorted;

  // The mesh and element connectivity
  int num_nodes, num_elements;
  int *elem_node_con, *elem_node_ptr;
  int *elem_component;
  double *Xpts;

  // Store information about the components
  int num_components;
  char * component_elems;
  char * component_descript;

  int num_owned_elements;
  int * local_component_nums;

  // The boundary conditions
  int num_bcs;
  int *bc_nodes, *bc_con, *bc_ptr;
  int *orig_bc_nodes, *orig_bc_con, *orig_bc_ptr;
  double *bc_vals, *orig_bc_vals;
};

#endif
