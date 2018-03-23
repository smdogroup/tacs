/*
  This file is part of TACS: The Toolkit for the Analysis of Composite
  Structures, a parallel finite-element code for structural and
  multidisciplinary design optimization.

  Copyright (C) 2010 University of Toronto
  Copyright (C) 2012 University of Michigan
  Copyright (C) 2014 Georgia Tech Research Corporation
  Additional copyright (C) 2010 Graeme J. Kennedy and Joaquim
  R.R.A. Martins All rights reserved.

  TACS is licensed under the Apache License, Version 2.0 (the
  "License"); you may not use this software except in compliance with
  the License.  You may obtain a copy of the License at
  
  http://www.apache.org/licenses/LICENSE-2.0 
*/

#ifndef TACS_MESH_LOADER_H
#define TACS_MESH_LOADER_H

/*
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

#include "TACSObject.h"
#include "TACSElement.h"
#include "TACSAssembler.h"
#include "TACSToFH5.h"
#include "TACSAuxElements.h"
#include "TACSCreator.h"

class TACSMeshLoader : public TACSObject {
 public:
  TACSMeshLoader( MPI_Comm _comm );
  ~TACSMeshLoader();

  // Read a BDF file for input
  // -------------------------
  int scanBDFFile( const char *file_name );

  // Get information about the mesh after scanning
  // ---------------------------------------------
  int getNumComponents();
  const char *getComponentDescript( int comp_num );
  const char *getElementDescript( int comp_num );

  // Set the elements corresponding to each of the component numbers 
  // ---------------------------------------------------------------
  void setElement( int component_num, TACSElement *_element );
  void setConvertToCoordinate( int flag );

  // Retrieve the element numbers corresponding to the given 
  // component numbers
  // -------------------------------------------------------
  int getNumNodes();
  int getNumElements();

  // Create a TACSToFH5 file writer
  // ------------------------------
  TACSToFH5 *createTACSToFH5( TACSAssembler *tacs,
                              enum ElementType elem_type,
                              unsigned int write_flag );

  // Distribute the mesh and create TACS
  // -----------------------------------
  TACSAssembler *createTACS( int vars_per_node,
                             TACSAssembler::OrderingType order_type = 
                             TACSAssembler::NATURAL_ORDER,
                             TACSAssembler::MatrixOrderingType mat_type = 
                             TACSAssembler::DIRECT_SCHUR);

  // Set the domain of a structural function with component numbers
  // --------------------------------------------------------------
  void addFunctionDomain( TACSFunction *function, 
                          int comp_nums[], int num_comps );

  // Add the auxiliary element to the given component
  // ------------------------------------------------
  void addAuxElement( TACSAuxElements *aux, int component_num,
                      TACSElement *_element );

  void getConnectivity( int *_num_nodes, int *_num_elements,
                        const int **_elem_node_ptr, 
                        const int **_elem_node_conn,
                        const int **_elem_component,
                        const TacsScalar **_Xpts ){
    if (_num_nodes){ *_num_nodes = num_nodes; }
    if (_num_elements){ *_num_elements = num_elements; }
    if (_elem_node_ptr){ *_elem_node_ptr = elem_node_ptr; }
    if (_elem_node_conn){ *_elem_node_conn = elem_node_conn; }
    if (_elem_component){ *_elem_component = elem_component; }
    if (_Xpts){ *_Xpts = Xpts; }
  }
  void getBCs( int *_num_bcs, const int **_bc_nodes, const int **_bc_vars, 
               const int **_bc_ptr, const TacsScalar **_bc_vals ){
    if (_num_bcs){ *_num_bcs = num_bcs; }
    if (_bc_nodes){ *_bc_nodes = bc_nodes; }
    if (_bc_vars){ *_bc_vars = bc_vars; }
    if (_bc_ptr){ *_bc_ptr = bc_ptr; }
    if (_bc_vals){ *_bc_vals = bc_vals; }
  }

 private:
  // Communicator for all processors
  MPI_Comm comm;

  // The underlying creator object
  TACSCreator *creator;

  // The element corresponding to each of the component numbers
  TACSElement **elements;

  // Original BDF mesh information: Note that the original
  // ordering may not be contiguous. The node numbers associated
  // with the original ordering are stored in node_nums.
  int *node_nums;
  double *Xpts_unsorted; 

  // Reduced set of contiguous nodes 
  TacsScalar *Xpts;

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
  TacsScalar *bc_vals;

  // Convert to TACS's coordinate ordering if necessary (need to
  // convert if loading BDF from gmsh output)
  int convertToCoordinate;
};

#endif // TACS_MESH_LOADER_H
