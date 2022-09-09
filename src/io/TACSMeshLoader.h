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

const int TacsMeshLoaderNumElementTypes = 10;

const char *TacsMeshLoaderElementTypes[] = {
    "CBAR",  "CQUADR",  "CQUAD4", "CQUAD8", "CQUAD9",
    "CQUAD", "CHEXA27", "CHEXA",  "CTRIA3", "CTETRA"};

// Lower and upper limits for the number of nodes
const int TacsMeshLoaderElementLimits[][2] = {{2, 2},    // CBAR
                                              {4, 4},    // CQUADR
                                              {4, 4},    // CQUAD4
                                              {8, 8},    // CQUAD8
                                              {9, 9},    // CQUAD9
                                              {9, 9},    // CQUAD
                                              {27, 27},  // CHEXA27
                                              {8, 8},    // CHEXA
                                              {3, 3},    // CTRIA3
                                              {4, 10}};  // CTETRA

/*
  This class provides a limited capability to read in nodal and
  connectivity information from a .bdf file - and could be extended
  to read in other formats. The class can also be used to distribute
  the mesh over a set of processors.

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

#include "TACSAuxElements.h"
#include "TACSCreator.h"
#include "TACSToFH5.h"

class TACSMeshLoader : public TACSObject {
 public:
  TACSMeshLoader(MPI_Comm _comm);
  ~TACSMeshLoader();

  // Read a BDF file for input
  // -------------------------
  int scanBDFFile(const char *file_name);

  // Get information about the mesh after scanning
  // ---------------------------------------------
  int getNumComponents();
  const char *getComponentDescript(int comp_num);
  const char *getElementDescript(int comp_num);

  // Set the elements corresponding to each of the component numbers
  // ---------------------------------------------------------------
  void setElement(int component_num, TACSElement *_element);

  // Retrieve the element numbers corresponding to the given
  // component numbers
  // -------------------------------------------------------
  int getNumNodes();
  int getNumElements();

  // Create a TACSToFH5 file writer
  // ------------------------------
  TACSToFH5 *createTACSToFH5(TACSAssembler *tacs, ElementType elem_type,
                             int write_flag);

  // Distribute the mesh and create TACS
  // -----------------------------------
  TACSAssembler *createTACS(
      int vars_per_node,
      TACSAssembler::OrderingType order_type = TACSAssembler::NATURAL_ORDER,
      TACSAssembler::MatrixOrderingType mat_type = TACSAssembler::DIRECT_SCHUR);

  // Set the domain of a structural function with component numbers
  // --------------------------------------------------------------
  void addFunctionDomain(TACSFunction *function, int num_comps,
                         int comp_nums[]);

  // Add the auxiliary element to the given component
  // ------------------------------------------------
  void addAuxElement(TACSAuxElements *aux, int component_num,
                     TACSElement *_element);

  // Get the node numbers in the Assembler object from the file number
  // -----------------------------------------------------------------
  void getAssemblerNodeNums(TACSAssembler *assembler, int num_nodes,
                            int *node_nums, int *num_new_nodes,
                            int **new_nodes);

  // Get the connectivity and boundary conditions
  // --------------------------------------------
  void getConnectivity(int *_num_nodes, int *_num_elements,
                       const int **_elem_node_ptr, const int **_elem_node_conn,
                       const int **_elem_component, const TacsScalar **_Xpts);
  void getBCs(int *_num_bcs, const int **_bc_nodes, const int **_bc_vars,
              const int **_bc_ptr, const TacsScalar **_bc_vals);

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
  int *file_node_nums;
  int *node_arg_sort_list;
  int *file_elem_nums;
  int *elem_arg_sort_list;

  // Reduced set of contiguous nodes
  TacsScalar *Xpts;

  // Store information about the original node numbers,
  // element numbers from the file
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
};

#endif  // TACS_MESH_LOADER_H
