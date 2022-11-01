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

#include "TACSAssembler.h"

#include "TACSElementVerification.h"
#include "TacsUtilities.h"

// Reordering implementation
#include "AMDInterface.h"

// Include the AMD package if we have it
#ifdef TACS_HAS_AMD_LIBRARY
#include "amd.h"
#endif  // TACS_HAS_AMD_LIBRARY

// The TACS-METIS header
#include "tacsmetis.h"

// TACS-BLAS/LAPACK header
#include "tacslapack.h"

/**
  Constructor for the TACSAssembler object

  @param tacs_comm the TACS communicator
  @param varsPerNode the number of degrees of freedom per node
  @param numOwnedNodes the number of locally-owned nodes
  @param numElements the number of elements in the mesh
  @param numDependentNodes the number of dependent nodes in the mesh
*/
TACSAssembler::TACSAssembler(MPI_Comm _tacs_comm, int _varsPerNode,
                             int _numOwnedNodes, int _numElements,
                             int _numDependentNodes) {
  TacsInitialize();

  // Copy the communicator for MPI communication
  tacs_comm = _tacs_comm;

  // If MPI is being used, get the rank of the current
  // process and the total number of processes
  MPI_Comm_rank(tacs_comm, &mpiRank);
  MPI_Comm_size(tacs_comm, &mpiSize);

  // Set the simulation time to 0
  time = 0.0;

  // Now set up the default pthread info with 1 thread
  thread_info = new TACSThreadInfo(1);
  thread_info->incref();
  pthread_mutex_init(&tacs_mutex, NULL);

  // Create the class that is used to
  tacsPInfo = new TACSAssemblerPthreadInfo();
  numCompletedElements = 0;

  // copy data to be used later in the program
  varsPerNode = _varsPerNode;
  numElements = _numElements;
  numOwnedNodes = _numOwnedNodes;
  numDependentNodes = _numDependentNodes;

  // These values will be computed later
  numMultiplierNodes = 0;
  numExtNodes = 0;
  numNodes = 0;
  extNodeOffset = 0;

  // Print out the number of local nodes and elements
  printf(
      "[%d] Creating TACSAssembler with numOwnedNodes = %d "
      "numElements = %d\n",
      mpiRank, numOwnedNodes, numElements);

  // Calculate the total number of nodes and elements
  int info[2], recv_info[2];
  info[0] = numOwnedNodes;
  info[1] = numElements;
  MPI_Reduce(info, recv_info, 2, MPI_INT, MPI_SUM, 0, tacs_comm);
  if (mpiSize > 1 && mpiRank == 0) {
    printf(
        "[%d] TACSAssembler: Total dof = %d Total nodes = %d "
        "Total elements = %d\n",
        mpiRank, varsPerNode * recv_info[0], recv_info[0], recv_info[1]);
  }

  // The mesh has not been initialized yet...
  meshInitializedFlag = 0;

  // Initialize the maximum element propertie
  maxElementDesignVars = 0;
  maxElementNodes = 0;
  maxElementSize = 0;
  maxElementIndepNodes = 0;

  // Set the elements array to NULL
  elements = NULL;

  // Set the auxiliary element class to NULL
  auxElements = NULL;

  // Information for setting boundary conditions and distributing variables
  nodeMap = new TACSNodeMap(tacs_comm, numOwnedNodes);
  nodeMap->incref();

  // Estimate 100 bcs at first, but this is expanded as required
  int nbc_est = 100;
  bcMap = new TACSBcMap(varsPerNode, nbc_est);
  bcInitMap = new TACSBcMap(varsPerNode, nbc_est);
  bcMap->incref();
  bcInitMap->incref();

  // Set the internal vector values to NULL
  varsVec = NULL;
  dvarsVec = NULL;
  ddvarsVec = NULL;
  xptVec = NULL;

  // Set the external node numbers to NULL
  tacsExtNodeNums = NULL;

  // Set the distribution object to NULL at first
  extDist = NULL;
  extDistIndices = NULL;

  // Reordering information old node i -> newNodeIndices[i]
  newNodeIndices = NULL;

  // parMat specific objects
  parMatIndices = NULL;

  // TACSSchurMat-specific objects
  schurBIndices = schurCIndices = NULL;
  schurBMap = schurCMap = NULL;

  // Allocate element-> node information
  elementNodeIndex = NULL;
  elementTacsNodes = NULL;

  // Null out the dependent node data
  depNodes = NULL;

  // Design variable information
  designVarsPerNode = 1;
  designNodeMap = NULL;
  designExtDist = NULL;
  designDepNodes = NULL;

  // Set the local element data to NULL
  elementData = NULL;
  elementIData = NULL;
  elementSensData = NULL;
  elementSensIData = NULL;

  // Initial condition vectors
  vars0 = NULL;
  dvars0 = NULL;
  ddvars0 = NULL;
}

/**
   Clean up the allocated memory and decref() all objects
*/
TACSAssembler::~TACSAssembler() {
  TacsFinalize();

  pthread_mutex_destroy(&tacs_mutex);
  delete tacsPInfo;

  // Go through and decref all the elements
  if (elements) {
    for (int i = 0; i < numElements; i++) {
      if (elements[i]) {
        elements[i]->decref();
      }
    }
    delete[] elements;
  }

  // Decref the variables/node vectors
  if (varsVec) {
    varsVec->decref();
  }
  if (dvarsVec) {
    dvarsVec->decref();
  }
  if (ddvarsVec) {
    ddvarsVec->decref();
  }
  if (xptVec) {
    xptVec->decref();
  }

  // Delete nodal information
  if (elementNodeIndex) {
    delete[] elementNodeIndex;
  }
  if (elementTacsNodes) {
    delete[] elementTacsNodes;
  }
  if (tacsExtNodeNums) {
    delete[] tacsExtNodeNums;
  }

  // Decreate the ref. count for the dependent node information
  if (depNodes) {
    depNodes->decref();
  }

  // Decrease the reference count to the auxiliary elements
  if (auxElements) {
    auxElements->decref();
  }

  // Decrease the reference count to objects allocated in initialize
  if (nodeMap) {
    nodeMap->decref();
  }
  if (bcMap) {
    bcMap->decref();
  }
  if (bcInitMap) {
    bcInitMap->decref();
  }
  if (extDist) {
    extDist->decref();
  }
  if (extDistIndices) {
    extDistIndices->decref();
  }

  // Free the reordering if it has been used
  if (newNodeIndices) {
    newNodeIndices->decref();
  }

  // Decrease ref. count for TACSParallelMat data
  if (parMatIndices) {
    parMatIndices->decref();
  }

  // Decrease ref. count for the TACSSchurMat data if it is allocated
  if (schurBIndices) {
    schurBIndices->decref();
  }
  if (schurCIndices) {
    schurCIndices->decref();
  }
  if (schurBMap) {
    schurBMap->decref();
  }
  if (schurCMap) {
    schurCMap->decref();
  }

  // Delete arrays allocated in initializeArrays()
  if (elementData) {
    delete[] elementData;
  }
  if (elementIData) {
    delete[] elementIData;
  }
  if (elementSensData) {
    delete[] elementSensData;
  }
  if (elementSensIData) {
    delete[] elementSensIData;
  }

  // Delete initial condition vectors
  if (vars0) {
    vars0->decref();
  }
  if (dvars0) {
    dvars0->decref();
  }
  if (ddvars0) {
    ddvars0->decref();
  }

  // Decref the thread information class
  thread_info->decref();
}

const char *TACSAssembler::tacsName = "TACSAssembler";

/**
   Return the MPI communicator for the TACSAssembler object

   @return the MPI communicator
*/
MPI_Comm TACSAssembler::getMPIComm() { return tacs_comm; }

/**
  Set the number of threads to use in the computation

  @param t The number of threads
*/
void TACSAssembler::setNumThreads(int t) { thread_info->setNumThreads(t); }

/**
   Return the thread information from the TACSAssembler object

   @return the thread info object
*/
TACSThreadInfo *TACSAssembler::getThreadInfo() { return thread_info; }

/**
   Get the number of degrees of freedom per node

   @return the number of degrees of freedom per node
*/
int TACSAssembler::getVarsPerNode() { return varsPerNode; }

/**
   Get the number of design variables per "node"

   @return the number of design variables for each design "node"
*/
int TACSAssembler::getDesignVarsPerNode() { return designVarsPerNode; }

/**
  Get the number of local nodes

  @return the number of local nodes defined on this processor
*/
int TACSAssembler::getNumNodes() { return numNodes; }

/**
  Get the number of dependent nodes

  @return the number of dependent nodes defined on this processor
*/
int TACSAssembler::getNumDependentNodes() { return numDependentNodes; }

/**
  Get the number of owned local nodes

  @return the number of locally owned nodes on this processor
*/
int TACSAssembler::getNumOwnedNodes() { return numOwnedNodes; }

/**
  Get the number of elements

  @return The number of elements defined on this processor
*/
int TACSAssembler::getNumElements() { return numElements; }

/**
  Get the node-processor assignment

  @return The map that defines the nodal assignment to each processor
*/
TACSNodeMap *TACSAssembler::getNodeMap() { return nodeMap; }

/**
  Get the node-processor design variable assignment

  @return The map that defines the design variable assignment to each processor
*/
TACSNodeMap *TACSAssembler::getDesignNodeMap() { return designNodeMap; }

/**
  Get the boundary conditions

  @return The boundary conditions (may be NULL)
*/
TACSBcMap *TACSAssembler::getBcMap() { return bcMap; }

/**
  Get the additional boundary conditions associated with the initial
  conditions.

  @return Initial boundary conditions (may be NULL)
*/
TACSBcMap *TACSAssembler::getInitBcMap() { return bcInitMap; }

/**
  Get the vector distribution object

  @return The object for distributing vector values
*/
TACSBVecDistribute *TACSAssembler::getBVecDistribute() { return extDist; }

/**
  Get the dependent node vectors

  @return The object used to define the dependent node equations
*/
TACSBVecDepNodes *TACSAssembler::getBVecDepNodes() { return depNodes; }

/**
  Get the maximum number of nodes for an element

  @return The maximum number of nodes for an element
*/
int TACSAssembler::getMaxElementNodes() { return maxElementNodes; }

/**
  Get the maximum number of element variables

  @return The maximum number of variables for an element
*/
int TACSAssembler::getMaxElementVariables() { return maxElementSize; }

/**
  Get the maximum number of design variables for an element

  @return The maximum number of design variables for an element
*/
int TACSAssembler::getMaxElementDesignVars() { return maxElementDesignVars; }

/**
  Get the array of elements from TACSAssembler

  @return The array of element objects
*/
TACSElement **TACSAssembler::getElements() { return elements; }

/**
  Get the element object and the corresponding element variables

  @param index The element index
  @param Xpts The element nodes array (may be NULL)
  @param vars The element variable array (may be NULL)
  @param dvars The element variable time derivative array (may be NULL)
  @param dvars The element variable time derivative array (may be NULL)
  @return The element object (NULL if index is out of range)
*/
TACSElement *TACSAssembler::getElement(int index, TacsScalar *Xpts,
                                       TacsScalar *vars, TacsScalar *dvars,
                                       TacsScalar *ddvars) {
  if (elements && xptVec && (index >= 0 && index < numElements)) {
    int ptr = elementNodeIndex[index];
    int len = elementNodeIndex[index + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    if (Xpts) {
      xptVec->getValues(len, nodes, Xpts);
    }
    if (vars) {
      varsVec->getValues(len, nodes, vars);
    }
    if (dvars) {
      dvarsVec->getValues(len, nodes, dvars);
    }
    if (ddvars) {
      ddvarsVec->getValues(len, nodes, ddvars);
    }

    return elements[index];
  }

  return NULL;
}

/**
  Get the element object and the pointer to the node numbers

  @param index The element index
  @param nodes The global node numbers
  @param len The number of nodes for this element
  @return The element object (NULL if index is out of range)
*/
TACSElement *TACSAssembler::getElement(int index, int *len, const int **nodes) {
  if (elements && (index >= 0 && index < numElements)) {
    int ptr = elementNodeIndex[index];
    if (len) {
      *len = elementNodeIndex[index + 1] - ptr;
    }
    if (nodes) {
      *nodes = &elementTacsNodes[ptr];
    }
    return elements[index];
  }

  if (len) {
    *len = 0;
  }
  if (nodes) {
    *nodes = NULL;
  }
  return NULL;
}

/**
  Set the element connectivity.

  Note that the number of elements that are set at this stage must be
  consistent with the number of elements passed in when TACSAssembler
  is created.  (The connectivity arrays are copied and should be freed
  by the caller.)

  @param ptr Offset into the connectivity array for this processor
  @param conn The connectivity from elements to global node index
  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::setElementConnectivity(const int *ptr, const int *conn) {
  if (meshInitializedFlag) {
    fprintf(stderr,
            "[%d] Cannot call setElementConnectivity() after initialize()\n",
            mpiRank);
    return 1;
  }
  if (tacsExtNodeNums) {
    fprintf(stderr,
            "[%d] Cannot call setElementConnectivity() after reordering\n",
            mpiRank);
    return 1;
  }

  // Free the data if it has already been set before
  if (elementTacsNodes) {
    delete[] elementTacsNodes;
  }
  if (elementNodeIndex) {
    delete[] elementNodeIndex;
  }

  int size = ptr[numElements];
  elementNodeIndex = new int[numElements + 1];
  elementTacsNodes = new int[size];
  memcpy(elementNodeIndex, ptr, (numElements + 1) * sizeof(int));

  // Copy the node numbers
  memcpy(elementTacsNodes, conn, size * sizeof(int));

  // If the elements are set, check the connectivity
  if (elements) {
    for (int i = 0; i < numElements; i++) {
      int size = elementNodeIndex[i + 1] - elementNodeIndex[i];
      if (size != elements[i]->getNumNodes()) {
        fprintf(stderr,
                "[%d] TACSAssembler: Element %s does not match "
                " connectivity\n",
                mpiRank, elements[i]->getObjectName());
        return 1;
      }
    }
  }

  // Check that the node numbers are all within range and that the
  // dependent node numbers (if any) are in range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  for (int i = 0; i < numElements; i++) {
    int jend = elementNodeIndex[i + 1];
    for (int j = elementNodeIndex[i]; j < jend; j++) {
      if (elementTacsNodes[j] >= ownerRange[mpiSize]) {
        fprintf(stderr,
                "[%d] TACSAssembler: Element %d contains node number "
                "out of range\n",
                mpiRank, i);
        return -1;
      } else if (elementTacsNodes[j] < -numDependentNodes) {
        fprintf(stderr,
                "[%d] TACSAssembler: Element %d contains dependent "
                "node number out of range\n",
                mpiRank, i);
        return -1;
      }
    }
  }
  return 0;
}

/**
  Get the element connectivity.

  @param ptr Offset into the connectivity array for this processor
  @param conn The connectivity from elements to global node index
*/
void TACSAssembler::getElementConnectivity(const int **ptr, const int **conn) {
  if (!meshInitializedFlag) {
    fprintf(stderr,
            "[%d] Cannot call getElementConnectivity() before initialize()\n",
            mpiRank);
    if (ptr) {
      *ptr = NULL;
    }
    if (conn) {
      *conn = NULL;
    }
    return;
  }

  // Free the data if it has already been set before
  if (ptr) {
    *ptr = elementNodeIndex;
  }
  if (conn) {
    *conn = elementTacsNodes;
  }
}

/**
  Set the element array within TACS.

  The number of element pointers provided should be equal to the
  number of elements specified in the TACSAssembler constructor.  Note
  that the elements themselves are not copied, just the pointers to
  them.

  @param elements The array of element pointers with length = numElements.
  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::setElements(TACSElement **_elements) {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call setElements() after initialize()\n",
            mpiRank);
    return 1;
  }

  // Check if the number of variables per node matches
  for (int i = 0; i < numElements; i++) {
    if (_elements[i]->getVarsPerNode() != varsPerNode) {
      fprintf(stderr,
              "[%d] TACSAssembler: Element %s, num displacements %d "
              "does not match variables per node %d\n",
              mpiRank, _elements[i]->getObjectName(),
              _elements[i]->getVarsPerNode(), varsPerNode);
      return 1;
    }
  }

  // Copy over the element pointers into a local array
  if (elements) {
    for (int i = 0; i < numElements; i++) {
      _elements[i]->incref();
      if (elements[i]) {
        elements[i]->decref();
      }
      elements[i] = _elements[i];
    }
  } else {
    elements = new TACSElement *[numElements];
    memset(elements, 0, numElements * sizeof(TACSElement *));
    for (int i = 0; i < numElements; i++) {
      _elements[i]->incref();
      elements[i] = _elements[i];
    }
  }

  // Keep track of the number of multiplier nodes
  numMultiplierNodes = 0;

  // Determine the maximum number of nodes per element
  maxElementDesignVars = 0;
  maxElementSize = 0;
  maxElementNodes = 0;

  for (int i = 0; i < numElements; i++) {
    // Determine if the maximum number of variables and nodes needs to
    // be changed
    int elemSize = elements[i]->getNumVariables();
    if (elemSize > maxElementSize) {
      maxElementSize = elemSize;
    }

    elemSize = elements[i]->getNumNodes();
    if (elemSize > maxElementNodes) {
      maxElementNodes = elemSize;
    }

    // Count up the multiplier nodes
    int multiplier = elements[i]->getMultiplierIndex();
    if (multiplier > 0) {
      numMultiplierNodes++;
    }

    // Keep track of the maximum number of element design variables
    int numDVs = elements[i]->getDesignVarNums(i, 0, NULL);
    if (numDVs > maxElementDesignVars) {
      maxElementDesignVars = numDVs;
    }
  }

  // If the connectivity is set, determine if it is consistent
  if (elementNodeIndex) {
    for (int i = 0; i < numElements; i++) {
      int size = elementNodeIndex[i + 1] - elementNodeIndex[i];
      if (size != elements[i]->getNumNodes()) {
        fprintf(stderr,
                "[%d] TACSAssembler: Element %s does not match "
                "connectivity\n",
                mpiRank, elements[i]->getObjectName());
        return 1;
      }
    }
  }

  return 0;
}

/**
  Set the dependent node data structure

  The dependent nodes are designated by a negative index with a 1-base
  numbering scheme. This dependent node data structure is a sparse
  matrix with constant weights that designates the weights applied to
  the global independent nodes (storred in depNodeToTacs).

  @param depNodePtr Offset into the connectivity/weight arrays
  @param depNodeToTacs The independent global nodes
  @param depNodeWeights The weight applied to each independent node
  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::setDependentNodes(const int *_depNodePtr,
                                     const int *_depNodeToTacs,
                                     const double *_depNodeWeights) {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call setDependentNodes() after initialize()\n",
            mpiRank);
    return 1;
  }
  if (tacsExtNodeNums) {
    fprintf(stderr, "[%d] Cannot call setDependentNodes() after reordering\n",
            mpiRank);
    return 1;
  }

  // Free the data if the dependent nodes have already been set
  if (depNodes) {
    depNodes->decref();
  }

  // Get the ownership range of the nodes
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  // Check that all the independent nodes are positive and are within an
  // allowable range
  for (int i = 0; i < numDependentNodes; i++) {
    for (int jp = _depNodePtr[i]; jp < _depNodePtr[i + 1]; jp++) {
      if (_depNodeToTacs[jp] >= ownerRange[mpiSize]) {
        fprintf(stderr,
                "[%d] Dependent node %d contains node number %d out of range\n",
                mpiRank, i, _depNodeToTacs[jp]);
        return 1;
      } else if (_depNodeToTacs[jp] < 0) {
        fprintf(stderr, "[%d] Dependent node %d contains dependent node %d\n",
                mpiRank, i, _depNodeToTacs[jp]);
        return 1;
      }
    }
  }

  if (numDependentNodes > 0) {
    // Allocate the new memory and copy over the data
    int *depNodePtr = new int[numDependentNodes + 1];
    memcpy(depNodePtr, _depNodePtr, (numDependentNodes + 1) * sizeof(int));

    int size = depNodePtr[numDependentNodes];
    int *depNodeToTacs = new int[size];
    memcpy(depNodeToTacs, _depNodeToTacs, size * sizeof(int));

    double *depNodeWeights = new double[size];
    memcpy(depNodeWeights, _depNodeWeights, size * sizeof(double));

    // Allocate the dependent node data structure
    depNodes = new TACSBVecDepNodes(numDependentNodes, &depNodePtr,
                                    &depNodeToTacs, &depNodeWeights);
    depNodes->incref();
  }

  return 0;
}

/**
  Add Dirichlet boundary conditions.

  These BCs are added to an object that will be associated with the
  vectors/matrices that are created using TACSAssembler.

  Fully-fixed boundary conditions can be specified by passing in
  nbcs < 0. If the vals array is NULL, then it is assumed that all values
  are fixed to zero.

  @param nnodes The number of nodes
  @param nodes The node numbers
  @param nbcs The number of nodal variables to be constrained
  @param vars The variable indices satisfying (0 <= vars[i] < varsPerNode)
  @param vals The boundary condition values
*/
void TACSAssembler::addBCs(int nnodes, const int *nodes, int nbcs,
                           const int *vars, const TacsScalar *vals) {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call addBC() after initialize()\n", mpiRank);
    return;
  }

  // Adjust the input to the addBC call. If no vars are specified,
  // set the number of boundary conditions equal to the number of
  // variables per node
  if (!vars || (nbcs < 0)) {
    nbcs = varsPerNode;
  }

  // Get the ownership range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  // Add all the boundary conditions within the specified owner range
  for (int i = 0; i < nnodes; i++) {
    if (nodes[i] >= ownerRange[mpiRank] && nodes[i] < ownerRange[mpiRank + 1]) {
      bcMap->addBC(nodes[i], nbcs, vars, vals);
    }
  }
}

/**
  Add Dirichlet boundary conditions for the initial conditions.

  These BCs are added to an object that will be associated with the
  vectors/matrices that are created using TACSAssembler.

  Fully-fixed boundary conditions can be specified by passing in
  nbcs < 0. If the vals array is NULL, then it is assumed that all values
  are fixed to zero.

  @param nnodes The number of nodes
  @param nodes The node numbers
  @param nbcs The number of nodal variables to be constrained
  @param vars The variable indices satisfying (0 <= vars[i] < varsPerNode)
  @param vals The boundary condition values
*/
void TACSAssembler::addInitBCs(int nnodes, const int *nodes, int nbcs,
                               const int *vars, const TacsScalar *vals) {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call addInitBC() after initialize()\n",
            mpiRank);
    return;
  }

  // Adjust the input to the addBC call. If no vars are specified,
  // set the number of boundary conditions equal to the number of
  // variables per node
  if (!vars || (nbcs < 0)) {
    nbcs = varsPerNode;
  }

  // Get the ownership range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  // Add all the boundary conditions within the specified owner range
  for (int i = 0; i < nnodes; i++) {
    if (nodes[i] >= ownerRange[mpiRank] && nodes[i] < ownerRange[mpiRank + 1]) {
      bcInitMap->addBC(nodes[i], nbcs, vars, vals);
    }
  }
}

/**
  Set new Dirichlet BC values at nodes where BCs are imposed

  This takes the new boundary condition values as the entries in the
  given vector where the Dirichlet boundary conditions are imposed.

  @param vec The vector containing the new boundary condition values
*/
void TACSAssembler::setBCValuesFromVec(TACSBVec *vec) {
  TacsScalar *xvals = NULL;
  vec->getArray(&xvals);

  // Get the values from the boundary condition arrays
  const int *nodes, *vars;
  TacsScalar *values;
  int nbcs = bcMap->getBCs(&nodes, &vars, &values);

  // Get the ownership range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  for (int i = 0; i < nbcs; i++) {
    if (nodes[i] >= ownerRange[mpiRank] && nodes[i] < ownerRange[mpiRank + 1]) {
      int var = varsPerNode * (nodes[i] - ownerRange[mpiRank]);
      for (int k = 0; k < varsPerNode; k++) {
        if (vars[i] & (1 << k)) {
          // Scan through the rows to be zeroed
          values[varsPerNode * i + k] = xvals[var + k];
        }
      }
    }
  }
}

/**
  Create a global vector of node locations.

  Note that the nodal coordinates are not set into the vector.

  @return A new vector for spatial coordinates
*/
TACSBVec *TACSAssembler::createNodeVec() {
  return new TACSBVec(nodeMap, TACS_SPATIAL_DIM, extDist, depNodes);
}

/**
  Set the nodal locations from the input vector

  @param X The nodal coordinate vector
*/
void TACSAssembler::setNodes(TACSBVec *X) {
  xptVec->copyValues(X);

  // Distribute the values at this point
  xptVec->beginDistributeValues();
  xptVec->endDistributeValues();
}

/**
  Get the node locations from TACS

  @param X The nodal coordinates are copied into this vector
*/
void TACSAssembler::getNodes(TACSBVec *X) {
  if (X) {
    X->copyValues(xptVec);
  }
}

/**
  Get a pointer to the internal nodal coordinate vector

  @param X A pointer to the nodal coordinate vector used by this Assembler
  object
*/
void TACSAssembler::getNodes(TACSBVec **X) {
  if (X) {
    *X = xptVec;
  }
}

/**
  Set the auxiliary elements within the TACSAssembler object

  This only needs to be done once sometime during initialization.  If
  you need to change the loads repeatedly, this can be called
  repeatedly. No check is made at this point that you haven't done
  something odd. Note that the code assumes that the elements defined
  here perfectly overlap the non-zero pattern of the elements set
  internally within TACS already.

  @param auxElems Auxiliary element object
*/
void TACSAssembler::setAuxElements(TACSAuxElements *auxElems) {
  // Increase the reference count to the input elements (Note that
  // the input may be NULL to over-write the internal object
  if (auxElems) {
    auxElems->incref();
  }

  // Decrease the reference count if the old object is not NULL
  if (auxElements) {
    auxElements->decref();
  }
  auxElements = auxElems;

  // Check whether the auxiliary elements match
  if (auxElements) {
    int naux = 0;
    TACSAuxElem *aux = NULL;
    naux = auxElements->getAuxElements(&aux);
    for (int k = 0; k < naux; k++) {
      int elem = aux[k].num;
      if (elem < 0 || elem >= numElements) {
        fprintf(stderr, "[%d] Element number %d out of range [0,%d)\n", mpiRank,
                elem, numElements);
      } else if (elements[elem]->getNumVariables() !=
                 aux[k].elem->getNumVariables()) {
        fprintf(stderr, "[%d] Auxiliary element sizes do not match\n", mpiRank);
      }
    }
  }
}

/**
  Retrieve the auxiliary element object from TACSAssembler

  Note that the auxiliary element object may be NULL

  @return The auxiliary element object
*/
TACSAuxElements *TACSAssembler::getAuxElements() { return auxElements; }

/**
  Compute the external list of nodes and sort these nodes

  This code is called automatically by TACSAssembler and is private.
  This code computes numExtNodes and allocates the array
  tacsExtNodeNums for use in other functions. This function should
  only be called once during reordering or during initialization.

  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::computeExtNodes() {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call computeExtNodes() after initialize()\n",
            mpiRank);
    return 1;
  }
  if (!elementNodeIndex) {
    fprintf(stderr,
            "[%d] Cannot call computeExtNodes() element "
            "connectivity not defined\n",
            mpiRank);
    return 1;
  }

  // Get the ownership range of the nodes
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  // Find the maximum possible size of array
  int max_size = elementNodeIndex[numElements];

  // Get the dependent node connectivity information (if any) and
  // increase the size of the array to account for any external dependent
  // nodes we may find.
  const int *depNodePtr = NULL;
  const int *depNodeConn = NULL;
  if (depNodes) {
    depNodes->getDepNodes(&depNodePtr, &depNodeConn, NULL);
    max_size += depNodePtr[numDependentNodes];
  }

  // Keep track of the external nodes that we've found
  int ext_size = 0;
  int *ext_list = new int[max_size];

  // First loop over the element connectivity
  for (int i = 0; i < elementNodeIndex[numElements]; i++) {
    int node = elementTacsNodes[i];

    // Check if the node is external
    if ((node >= 0) &&
        (node < ownerRange[mpiRank] || node >= ownerRange[mpiRank + 1])) {
      ext_list[ext_size] = node;
      ext_size++;
    }
  }

  // Loop over the dependent nodes
  if (depNodes) {
    int end = depNodePtr[numDependentNodes];
    for (int i = 0; i < end; i++) {
      int node = depNodeConn[i];

      // Check if the node is external
      if ((node >= 0) &&
          (node < ownerRange[mpiRank] || node >= ownerRange[mpiRank + 1])) {
        ext_list[ext_size] = node;
        ext_size++;
      }
    }
  }

  // Sort the list of nodes
  numExtNodes = TacsUniqueSort(ext_size, ext_list);

  // Allocate an array of the external nodes that is tight
  // to the number of external nodes
  tacsExtNodeNums = new int[numExtNodes];
  memcpy(tacsExtNodeNums, ext_list, numExtNodes * sizeof(int));

  // Free the original list of external nodes;
  delete[] ext_list;

  // Now the total number of nodes is equal to the external nodes
  // plus the locally owned nodes
  numNodes = numOwnedNodes + numExtNodes;

  // Find the offset into the external node list
  extNodeOffset = 0;
  while (extNodeOffset < numExtNodes &&
         tacsExtNodeNums[extNodeOffset] < ownerRange[mpiRank]) {
    extNodeOffset++;
  }

  return 0;
}

/**
  Compute a reordering of the nodes.

  This function should be called after the element connectivity,
  boundary conditions and optionally the dependent nodes are set, but
  before the call to initialize().

  This code computes and stores a reordering based on both the matrix
  type and the ordering type. The matrix type determines what
  constitutes a coupling node and what does not. The matrix type only
  impacts the ordering in parallel computations.

  The algorithm proceeds as follows:

  1. Compute the coupling nodes referenced by another processor
  2. Order the locally owned nodes based on the input ordering.
  3. Based on the recvied coupling nodes, set the outgoing node
  numbers back into the recieving array.
  4. Set the new values of the nodes on the requesting processes

  @param order_type The type of ordering to apply
  @param mat_type The anticipated matrix distribution type
*/
void TACSAssembler::computeReordering(OrderingType order_type,
                                      MatrixOrderingType mat_type) {
  // Return if the element connectivity not set
  if (!elementNodeIndex) {
    fprintf(stderr, "[%d] Must define element connectivity before reordering\n",
            mpiRank);
    return;
  }
  if (tacsExtNodeNums) {
    fprintf(stderr,
            "[%d] TACSAssembler::computeReordering() can only be called once\n",
            mpiRank);
    return;
  }

  // Compute the external nodes
  computeExtNodes();

  // Compute the local node numbers that correspond to the coupling
  // nodes between processors
  int *couplingNodes;
  int *extPtr, *extCount;
  int *recvPtr, *recvCount, *recvNodes;
  int numCoupling = computeCouplingNodes(&couplingNodes, &extPtr, &extCount,
                                         &recvPtr, &recvCount, &recvNodes);

  // The new node numbers to be set according to the different
  // re-ordering schemes
  int *newNodeNums = new int[numNodes];

  // Metis can't handle the non-zero diagonal in the CSR data structure
  int removeDiagonal = 0;
  if (order_type == ND_ORDER) {
    removeDiagonal = 1;
  }

  // If using only one processor, order everything. In this case
  // there is no distinction between local and global ordering.
  if (mpiSize == 1) {
    // Compute the local (global since it's serial) node to node
    // connectivity based on the element connectivity.
    int *rowp, *cols;
    computeLocalNodeToNodeCSR(&rowp, &cols, removeDiagonal);

    // Compute the reordering
    computeMatReordering(order_type, numNodes, rowp, cols, NULL, newNodeNums);

    delete[] rowp;
    delete[] cols;
  } else if (mat_type == GAUSS_SEIDEL) {
    // Compute the local node to node connectivity
    int *rowp, *cols;
    computeLocalNodeToNodeCSR(&rowp, &cols, removeDiagonal);

    // Get the owner range
    const int *ownerRange;
    nodeMap->getOwnerRange(&ownerRange);

    // Compute the node types. This is the node type relative to its
    // owner.
    int *nodeType = new int[numNodes];

    for (int i = 0; i < numNodes; i++) {
      // Set the node as a local node type
      nodeType[i] = 0;

      // Get the global node number
      int node = getGlobalNodeNum(i);

      // Get the processor that owns this node
      int owner = nodeMap->getNodeOwner(node);

      // Search the columns for the type of node
      for (int jp = rowp[i]; jp < rowp[i + 1]; jp++) {
        // Get the local column index
        int j = cols[jp];

        // Get the global node number
        int col = getGlobalNodeNum(j);

        if (col < ownerRange[owner]) {
          if (nodeType[i] == 0) {
            // Top node type
            nodeType[i] = 1;
          } else if (nodeType[i] == 3) {
            // Mid-node type
            nodeType[i] = 2;
          }
        } else if (col >= ownerRange[owner + 1]) {
          if (nodeType[i] == 0) {
            // Bottom node type
            nodeType[i] = 3;
          } else if (nodeType[i] == 1) {
            // Mid-node type
            nodeType[i] = 2;
          }
        }
      }
    }

    // This type of connectivity information is no longer required
    delete[] rowp;
    delete[] cols;

    // We've computed the local node types, but this may be affected
    // by the contributions from other processors. Now, perform some
    // collective communication to globalize the node type
    // information. This will order all the internal nodes first,
    // followed by the top nodes, mid nodes and bottom nodes,
    // respectively.

    // Now allocate space for the node types
    int *sendNodeType = new int[numExtNodes];
    for (int i = 0; i < numExtNodes; i++) {
      // Get the local node number
      int index = getLocalNodeNum(tacsExtNodeNums[i]);
      sendNodeType[i] = nodeType[index];
    }

    int *recvNodeType = new int[recvPtr[mpiSize]];
    MPI_Alltoallv(sendNodeType, extCount, extPtr, MPI_INT, recvNodeType,
                  recvCount, recvPtr, MPI_INT, tacs_comm);

    // Update the internal node types based on the recv'd values
    for (int i = 0; i < recvPtr[mpiSize]; i++) {
      int node = getLocalNodeNum(recvNodes[i]);

      if (nodeType[node] == 0) {
        // Set the node type as the recv node type
        nodeType[node] = recvNodeType[i];
      } else if (nodeType[node] != recvNodeType[i]) {
        // Mid-node type
        nodeType[node] = 2;
      }
    }

    // Free the allocated data
    delete[] sendNodeType;
    delete[] recvNodeType;

    // So now the local nodes have the correct node types,
    // compute the local ordering.
    int *reducedNodes = new int[numNodes];
    memset(reducedNodes, 0, numNodes * sizeof(int));

    // Exclude all the nodes that are external to this processor
    for (int i = 0; i < extNodeOffset; i++) {
      reducedNodes[i] = -1;
    }
    for (int i = extNodeOffset + numOwnedNodes; i < numNodes; i++) {
      reducedNodes[i] = -1;
    }

    // Loop over all of the node types. This will order all the
    // internal nodes first
    int offset = ownerRange[mpiRank];
    for (int ntype = 0; ntype < 4; ntype++) {
      int numReducedNodes = 0;
      for (int i = extNodeOffset; i < extNodeOffset + numOwnedNodes; i++) {
        if (nodeType[i] == ntype) {
          reducedNodes[i] = numReducedNodes;
          numReducedNodes++;
        } else {
          reducedNodes[i] = -1;
        }
      }

      if (numReducedNodes > 0) {
        // Compute the reordering for the reduced set of nodes
        int *newReducedNodes = new int[numReducedNodes];
        int *reduced_rowp, *reduced_cols;
        computeLocalNodeToNodeCSR(&reduced_rowp, &reduced_cols, numReducedNodes,
                                  reducedNodes, removeDiagonal);
        computeMatReordering(order_type, numReducedNodes, reduced_rowp,
                             reduced_cols, NULL, newReducedNodes);
        delete[] reduced_rowp;
        delete[] reduced_cols;

        // Set the new variable numbers for the boundary nodes
        // and include their offset
        for (int i = 0, j = 0; i < numNodes; i++) {
          if (reducedNodes[i] >= 0) {
            newNodeNums[i] = offset + newReducedNodes[j];
            j++;
          }
        }

        offset += numReducedNodes;

        // Free the new node numbers
        delete[] newReducedNodes;
      }
    }

    // Free the node type array
    delete[] nodeType;
  } else {
    // First, find the reduced nodes - the set of nodes that are only
    // referenced by this processor. These can be reordered without
    // affecting other processors.
    int *reducedNodes = new int[numNodes];
    memset(reducedNodes, 0, numNodes * sizeof(int));

    // Exclude all the nodes that are external to this processor
    for (int i = 0; i < extNodeOffset; i++) {
      reducedNodes[i] = -1;
    }
    for (int i = extNodeOffset + numOwnedNodes; i < numNodes; i++) {
      reducedNodes[i] = -1;
    }

    // Depending on the matrix type, also add the external dependent
    // nodes and nodes that also couple to those nodes
    if (mat_type == DIRECT_SCHUR) {
      for (int i = 0; i < numCoupling; i++) {
        int node = couplingNodes[i];
        reducedNodes[node] = -1;
      }
    } else if (mat_type == APPROXIMATE_SCHUR) {
      // If we want an approximate schur ordering, where the nodes
      // that couple to other processors are ordered last, we also add
      // these nodes to the reduced set.
      int *rowp, *cols;
      computeLocalNodeToNodeCSR(&rowp, &cols);

      // Order all nodes linked by an equation to an external node
      // This ordering is required for the approximate Schur method
      for (int i = 0; i < numCoupling; i++) {
        int node = couplingNodes[i];
        for (int jp = rowp[node]; jp < rowp[node + 1]; jp++) {
          reducedNodes[cols[jp]] = -1;
        }
      }

      delete[] rowp;
      delete[] cols;
    }

    // Now, order all nodes that are non-negative
    int numReducedNodes = 0;
    for (int i = 0; i < numNodes; i++) {
      if (reducedNodes[i] >= 0) {
        reducedNodes[i] = numReducedNodes;
        numReducedNodes++;
      }
    }

    // Compute the reordering for the reduced set of nodes
    int *newReducedNodes = new int[numReducedNodes];
    int *rowp, *cols;
    computeLocalNodeToNodeCSR(&rowp, &cols, numReducedNodes, reducedNodes,
                              removeDiagonal);
    computeMatReordering(order_type, numReducedNodes, rowp, cols, NULL,
                         newReducedNodes);
    delete[] rowp;
    delete[] cols;

    // Place the result back into the newNodeNums - add the
    // ownership offset
    const int *ownerRange;
    nodeMap->getOwnerRange(&ownerRange);
    int offset = ownerRange[mpiRank];
    for (int i = 0, j = 0; i < numNodes; i++) {
      if (reducedNodes[i] >= 0) {
        newNodeNums[i] = offset + newReducedNodes[j];
        j++;
      }
    }

    delete[] newReducedNodes;

    // Add the offset to the total number of reduced nodes
    offset += numReducedNodes;

    // Now, order any remaining variables that have not yet been
    // ordered. These are the coupling variables (if any) that
    // have been labeled before.
    numReducedNodes = 0;
    for (int i = extNodeOffset; i < extNodeOffset + numOwnedNodes; i++) {
      // If the node has not been ordered and is within the ownership
      // range of this process, order it now.
      if (reducedNodes[i] < 0) {
        reducedNodes[i] = numReducedNodes;
        numReducedNodes++;
      } else {
        reducedNodes[i] = -1;
      }
    }

    if (numReducedNodes > 0) {
      // Additive Schwarz ordering should number all locally owned
      // nodes first, and should not require this second ordering.
      if (mat_type == ADDITIVE_SCHWARZ) {
        fprintf(stderr, "[%d] Error in additive Schwarz reordering\n", mpiRank);
      }

      // Order any remaning variables that are locally owned
      newReducedNodes = new int[numReducedNodes];
      computeLocalNodeToNodeCSR(&rowp, &cols, numReducedNodes, reducedNodes,
                                removeDiagonal);
      computeMatReordering(order_type, numReducedNodes, rowp, cols, NULL,
                           newReducedNodes);

      // Free the allocate CSR data structure
      delete[] rowp;
      delete[] cols;

      // Set the new variable numbers for the boundary nodes
      // and include their offset
      for (int i = 0, j = 0; i < numNodes; i++) {
        if (reducedNodes[i] >= 0) {
          newNodeNums[i] = offset + newReducedNodes[j];
          j++;
        }
      }

      // Free the new node numbers
      delete[] newReducedNodes;
    }

    delete[] reducedNodes;
  }

  // So now we have new node numbers for the nodes owned by this
  // processor, but the other processors do not have these new numbers
  // yet. Find the values assigned to the nodes requested from
  // external nodes recv_nodes is now an outgoing list of nodes to
  // other processes
  for (int i = 0; i < recvPtr[mpiSize]; i++) {
    int node = getLocalNodeNum(recvNodes[i]);
    recvNodes[i] = newNodeNums[node];
  }

  // Now send the new node numbers back to the other processors
  // that reference them. This also uses all-to-all communication.
  int *newExtNodes = new int[extPtr[mpiSize]];
  MPI_Alltoallv(recvNodes, recvCount, recvPtr, MPI_INT, newExtNodes, extCount,
                extPtr, MPI_INT, tacs_comm);

  // Once the new node numbers from other processors is received
  // apply these new node numbers back to the locally owned
  // reference numbers.
  for (int i = 0; i < extNodeOffset; i++) {
    newNodeNums[i] = newExtNodes[i];
  }
  for (int i = extNodeOffset; i < numExtNodes; i++) {
    newNodeNums[i + numOwnedNodes] = newExtNodes[i];
  }

  // Free the new external node numbers
  delete[] newExtNodes;

  // Reorder the local dependent node connectivity
  int *depConn = NULL;
  const int *depPtr = NULL;
  if (depNodes) {
    depNodes->getDepNodeReorder(&depPtr, &depConn);
    int end = depPtr[numDependentNodes];
    for (int i = 0; i < end; i++) {
      int node = getLocalNodeNum(depConn[i]);
      depConn[i] = newNodeNums[node];
    }
  }

  // Reorder the element connectivity
  int end = elementNodeIndex[numElements];
  for (int i = 0; i < end; i++) {
    int node = elementTacsNodes[i];
    if (node >= 0) {
      node = getLocalNodeNum(node);
      elementTacsNodes[i] = newNodeNums[node];
    }
  }

  // If boundary conditions are set, reorder them
  if (bcMap) {
    int *nodeNums;
    int nbcs = bcMap->getBCNodeNums(&nodeNums);
    for (int i = 0; i < nbcs; i++) {
      int node = getLocalNodeNum(nodeNums[i]);
      nodeNums[i] = newNodeNums[node];
    }
  }
  if (bcInitMap) {
    int *nodeNums;
    int nbcs = bcInitMap->getBCNodeNums(&nodeNums);
    for (int i = 0; i < nbcs; i++) {
      int node = getLocalNodeNum(nodeNums[i]);
      nodeNums[i] = newNodeNums[node];
    }
  }

  // Finaly, reorder the external node numbers
  for (int i = 0; i < extNodeOffset; i++) {
    tacsExtNodeNums[i] = newNodeNums[i];
  }
  for (int i = extNodeOffset; i < numExtNodes; i++) {
    tacsExtNodeNums[i] = newNodeNums[i + numOwnedNodes];
  }

  // Resort the external node numbers - these are already unique and
  // the extNodeOffset should not change either
  TacsUniqueSort(numExtNodes, tacsExtNodeNums);

  // Save the mapping to the new numbers for later reorderings
  newNodeIndices = new TACSBVecIndices(&newNodeNums, numNodes);
  newNodeIndices->incref();

  delete[] couplingNodes;
  delete[] extPtr;
  delete[] extCount;
  delete[] recvPtr;
  delete[] recvCount;
  delete[] recvNodes;
}

/**
  Compute the reordering for the given matrix.

  This uses either Reverse Cuthill-McKee (RCM_ORDER), Approximate
  Minimum Degree (AMD) or Nested Disection (ND) to compute a
  reordering of the variables.

  The input to the function is a CSR-type data structure of the
  matrix. Note that for ND (with the external package Metis), requires
  that the diagonal be eliminated from the CSR data structure. (This
  modified data structure can be computed with the no_diagonal flag
  when calling the CSR creation routines.)

  The matrix ordering routines compute a reordering such that:

  P*A*P^{T} has fewer non-zeros.

  The function returns an array new_vars such that:

  new_vars = P^{T} old_vars
  perm = P

  @param order_type The ordering type
  @param nvars The number of variables
  @param rowp Pointer into each new row of the matrix
  @param cols Column indices for each row
  @param perm Permutation array
  @param new_vars New variable array
*/
void TACSAssembler::computeMatReordering(OrderingType order_type, int nvars,
                                         int *rowp, int *cols, int *perm,
                                         int *new_vars) {
  int *_perm = perm;
  int *_new_vars = new_vars;
  if (!perm) {
    _perm = new int[nvars];
  }
  if (!new_vars) {
    _new_vars = new int[nvars];
  }

  if (order_type == RCM_ORDER) {
    // Compute the matrix reordering using RCM TACS' version
    // of the RCM algorithm
    int root_node = 0;
    int num_rcm_iters = 1;
    TacsComputeRCMOrder(nvars, rowp, cols, _new_vars, root_node, num_rcm_iters);

    if (perm) {
      for (int k = 0; k < nvars; k++) {
        perm[_new_vars[k]] = k;
      }
    }
  } else if (order_type == AMD_ORDER) {
#if TACS_HAS_AMD_LIBRARY
    // Use the approximate minimum degree ordering
    double control[AMD_CONTROL], info[AMD_INFO];
    amd_defaults(control);  // Use the default values
    amd_order(nvars, rowp, cols, _perm, control, info);

    if (new_vars) {
      for (int k = 0; k < nvars; k++) {
        new_vars[_perm[k]] = k;
      }
    }
#else
    int use_exact_degree = 0;
    int ncoupling_nodes = 0;
    int *coupling_nodes = NULL;
    amd_order_interface(nvars, rowp, cols, _perm, coupling_nodes,
                        ncoupling_nodes, 0, NULL, NULL, NULL, use_exact_degree);

    if (new_vars) {
      for (int k = 0; k < nvars; k++) {
        new_vars[_perm[k]] = k;
      }
    }
#endif  // TACS_HAS_AMD_LIBRARY
  } else if (order_type == ND_ORDER) {
    // Set the default options in METIS
    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // Use 0-based numbering
    options[METIS_OPTION_NUMBERING] = 0;
    METIS_NodeND(&nvars, rowp, cols, NULL, options, _perm, _new_vars);
  } else if (order_type == TACS_AMD_ORDER) {
    int use_exact_degree = 0;
    int ncoupling_nodes = 0;
    int *coupling_nodes = NULL;
    amd_order_interface(nvars, rowp, cols, _perm, coupling_nodes,
                        ncoupling_nodes, 0, NULL, NULL, NULL, use_exact_degree);

    if (new_vars) {
      for (int k = 0; k < nvars; k++) {
        new_vars[_perm[k]] = k;
      }
    }
  } else if (order_type == MULTICOLOR_ORDER) {
    // Compute the matrix reordering using RCM TACS' version of the
    // RCM algorithm
    int *colors = new int[nvars];
    TacsComputeSerialMultiColor(nvars, rowp, cols, colors, new_vars);
    delete[] colors;

    if (perm) {
      for (int k = 0; k < nvars; k++) {
        perm[_new_vars[k]] = k;
      }
    }
  } else if (order_type == NATURAL_ORDER) {
    if (perm) {
      for (int k = 0; k < nvars; k++) {
        perm[k] = k;
      }
    }
    if (new_vars) {
      for (int k = 0; k < nvars; k++) {
        new_vars[k] = k;
      }
    }
  }

  if (!perm) {
    delete[] _perm;
  }
  if (!new_vars) {
    delete[] _new_vars;
  }
}

/**
  The following function returns a local node number based on the
  provided (global) TACS node number.

  If the node number is on this processor, no search is required,
  however, if the node is externally owned, then a binary search is
  needed to determine the index into the off-processor list of nodes.

  @param node The global TACS node number unique across all processors
  @return The local node number
*/
int TACSAssembler::getLocalNodeNum(int node) {
  // Get the ownership range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  if (node >= ownerRange[mpiRank] && node < ownerRange[mpiRank + 1]) {
    node = (node - ownerRange[mpiRank]) + extNodeOffset;
  } else if (node >= 0) {
    const int *ext_nodes = NULL;
    if (tacsExtNodeNums) {
      ext_nodes = tacsExtNodeNums;
    } else if (extDistIndices) {
      extDistIndices->getIndices(&ext_nodes);
    } else {
      fprintf(stderr, "[%d] External nodes not defined\n", mpiRank);
      return -1;
    }

    // Find the local index for external nodes
    int *item = TacsSearchArray(node, numExtNodes, ext_nodes);

    // Check if the item is found in the list
    if (item) {
      if (node < ownerRange[mpiRank]) {
        node = (item - ext_nodes);
      } else {
        node = numOwnedNodes + (item - ext_nodes);
      }
    } else {
      fprintf(stderr, "[%d] External node %d not in external node list\n",
              mpiRank, node);
      return -1;
    }
  } else {
    fprintf(stderr, "[%d] Cannot compute local number for dependent node %d\n",
            mpiRank, node);
    return -1;
  }

  return node;
}

/**
  Given the local node number find the corresponding global TACS node
  number

  This function is the inverse of the getLocalNodeNum() function
  defined above.

  @param node The local node number
  @return The global TACS node number
*/
int TACSAssembler::getGlobalNodeNum(int node) {
  // Get the ownership range
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  if (node < extNodeOffset) {
    const int *ext_nodes = NULL;
    if (tacsExtNodeNums) {
      ext_nodes = tacsExtNodeNums;
    } else if (extDistIndices) {
      extDistIndices->getIndices(&ext_nodes);
    } else {
      fprintf(stderr, "[%d] External nodes not defined\n", mpiRank);
      return -1;
    }

    return ext_nodes[node];
  } else if (node < extNodeOffset + numOwnedNodes) {
    return (node - extNodeOffset) + ownerRange[mpiRank];
  } else if (node < numNodes) {
    const int *ext_nodes = NULL;
    if (tacsExtNodeNums) {
      ext_nodes = tacsExtNodeNums;
    } else if (extDistIndices) {
      extDistIndices->getIndices(&ext_nodes);
    } else {
      fprintf(stderr, "[%d] External nodes not defined\n", mpiRank);
      return -1;
    }

    return ext_nodes[node - numOwnedNodes];
  } else {
    fprintf(stderr, "[%d] Local node number %d out of range\n", mpiRank, node);
    return -1;
  }

  return node;
}

/**
  The following function creates a data structure that links nodes
  to elements - this reverses the existing data structure that
  links elements to nodes but keeps the original in tact.

  The algorithm proceeds as follows:

  1. The size of the arrays are determined by finding how many
  nodes point to each element

  2. The index into the nodeElem array is determined by adding up
  the contributions from all previous entries.

  3. The original data structure is again traversed and this time
  an element number is associated with each element.

  @param _nodeElementPtr Pointer into the array nodeToElements
  @param _nodeToElements Element index associated with each node
*/
void TACSAssembler::computeNodeToElementCSR(int **_nodeElementPtr,
                                            int **_nodeToElements) {
  // Determine the node->element connectivity using local nodes
  int *nodeElementPtr = new int[numNodes + 1];
  memset(nodeElementPtr, 0, (numNodes + 1) * sizeof(int));

  // Get the dependent node connectivity information
  const int *depNodePtr = NULL;
  const int *depNodeConn = NULL;
  if (depNodes) {
    depNodes->getDepNodes(&depNodePtr, &depNodeConn, NULL);
  }

  // Loop over all the elements and count up the number of times
  // each node refers to each element
  for (int i = 0; i < numElements; i++) {
    int end = elementNodeIndex[i + 1];
    for (int jp = elementNodeIndex[i]; jp < end; jp++) {
      // Check whether this is locally owned or not
      int node = elementTacsNodes[jp];

      if (node >= 0) {
        node = getLocalNodeNum(node);
        nodeElementPtr[node + 1]++;
      } else if (node < 0) {
        // This is a dependent-node, determine which independent
        // nodes it depends on
        int dep_node = -node - 1;
        int kend = depNodePtr[dep_node + 1];
        for (int kp = depNodePtr[dep_node]; kp < kend; kp++) {
          node = getLocalNodeNum(depNodeConn[kp]);
          nodeElementPtr[node + 1]++;
        }
      }
    }
  }

  // Sum up the total size of the array
  for (int i = 0; i < numNodes; i++) {
    nodeElementPtr[i + 1] += nodeElementPtr[i];
  }

  // Allocate space for the nodeToElement connectivity
  int size = nodeElementPtr[numNodes];
  int *nodeToElements = new int[size];

  // Loop over the elements again, this time adding the nodes to the
  // connectivity
  for (int i = 0; i < numElements; i++) {
    int end = elementNodeIndex[i + 1];
    for (int jp = elementNodeIndex[i]; jp < end; jp++) {
      // Check whether this is locally owned or not
      int node = elementTacsNodes[jp];
      if (node >= 0) {
        node = getLocalNodeNum(node);
        nodeToElements[nodeElementPtr[node]] = i;
        nodeElementPtr[node]++;
      } else if (node < 0) {
        // This is a dependent-node, determine which independent
        // nodes it depends on
        int dep_node = -node - 1;
        int kend = depNodePtr[dep_node + 1];
        for (int kp = depNodePtr[dep_node]; kp < kend; kp++) {
          node = depNodeConn[kp];
          node = getLocalNodeNum(node);
          nodeToElements[nodeElementPtr[node]] = i;
          nodeElementPtr[node]++;
        }
      }
    }
  }

  // Set up the pointer array which denotes the start (and end) of each node
  for (int i = 0; i < numNodes; i++) {
    nodeElementPtr[numNodes - i] = nodeElementPtr[numNodes - i - 1];
  }
  nodeElementPtr[0] = 0;

  // Sort and unquify the CSR data structure
  TacsSortAndUniquifyCSR(numNodes, nodeElementPtr, nodeToElements);

  // Set the output pointers
  *_nodeToElements = nodeToElements;
  *_nodeElementPtr = nodeElementPtr;
}

/**
  Set up a CSR data structure pointing from local nodes to other
  local nodes.

  This function works by first estimating the number of entries in
  each row of the matrix. This information is stored temporarily in
  the array rowp. After the contributions from the elements and sparse
  constraints are added, the preceeding elements in rowp are added
  such that rowp becomes the row pointer for the matrix. Note that
  this is an upper bound because the elements and constraints may
  introduce repeated variables. Next, cols is allocated corresponding
  to the column index for each entry. This iterates back over all
  elements and constraints. At this stage, rowp is treated as an array
  of indices, that index into the i-th row of cols[:] where the next
  index should be inserted. As a result, rowp must be adjusted after
  this operation is completed.  The last step is to sort and uniquify
  each row of the matrix.

  @param _rowp The row pointer corresponding to CSR data structure
  @param cols The column indices for each row of the CSR data structure
  @param nodiag Flag to indicate whether to remove the diagonal matrix entry
*/
void TACSAssembler::computeLocalNodeToNodeCSR(int **_rowp, int **_cols,
                                              int nodiag) {
  int *cols = NULL;
  int *rowp = new int[numNodes + 1];
  memset(rowp, 0, (numNodes + 1) * sizeof(int));

  // Create the node -> element data structure
  int *nodeElementPtr = NULL;
  int *nodeToElements = NULL;
  computeNodeToElementCSR(&nodeElementPtr, &nodeToElements);

  // If we have dependent nodes, use a different algorithm
  if (depNodes) {
    const int *depNodePtr = NULL;
    const int *depNodeConn = NULL;
    depNodes->getDepNodes(&depNodePtr, &depNodeConn, NULL);

    // Count the number of nodes associated with each element
    int *nodeCount = new int[numElements];
    memset(nodeCount, 0, numElements * sizeof(int));

    for (int i = 0; i < numElements; i++) {
      int jend = elementNodeIndex[i + 1];
      for (int jp = elementNodeIndex[i]; jp < jend; jp++) {
        int node = elementTacsNodes[jp];
        if (node >= 0) {
          nodeCount[i]++;
        } else {
          // Add the number of independent nodes attached to this
          // dependent node for later use
          int dep = -node - 1;
          nodeCount[i] += depNodePtr[dep + 1] - depNodePtr[dep];
        }
      }
    }

    // First, populate rowp by finding the approximate number of
    // independent nodes per element
    for (int i = 0; i < numNodes; i++) {
      for (int jp = nodeElementPtr[i]; jp < nodeElementPtr[i + 1]; jp++) {
        int elem = nodeToElements[jp];
        rowp[i + 1] += nodeCount[elem];
      }
    }

    // Make a conservative estimate of the rowp pointer data
    for (int i = 0; i < numNodes; i++) {
      rowp[i + 1] += rowp[i];
    }

    // Set up the column indices for each row - label each one with
    // a negative index so that we know what has not been set
    int nnz = rowp[numNodes];
    cols = new int[nnz];
    for (int i = 0; i < nnz; i++) {
      cols[i] = -1;
    }

    // Add the element contribution to the column indices
    for (int i = 0; i < numNodes; i++) {
      for (int jp = nodeElementPtr[i]; jp < nodeElementPtr[i + 1]; jp++) {
        int elem = nodeToElements[jp];

        // Scan through all the nodes belonging to this element
        int kend = elementNodeIndex[elem + 1];
        int row = rowp[i];
        for (int kp = elementNodeIndex[elem]; kp < kend; kp++) {
          int node = elementTacsNodes[kp];
          if (node >= 0) {
            node = getLocalNodeNum(node);
            cols[row] = node;
            row++;
          } else {
            // This is a dependent-node, determine which independent
            // nodes it depends on
            int dep_node = -node - 1;
            int pend = depNodePtr[dep_node + 1];
            for (int p = depNodePtr[dep_node]; p < pend; p++) {
              node = depNodeConn[p];
              node = getLocalNodeNum(node);
              cols[row] = node;
              row++;
            }
          }
        }

        // Reset the pointer to this row
        rowp[i] = row;
      }
    }

    // Adjust rowp back to a zero-based index
    for (int i = numNodes; i > 0; i--) {
      rowp[i] = rowp[i - 1];
    }
    rowp[0] = 0;

    delete[] nodeCount;
  } else {
    // First, populate rowp by adding the contribution to a node from
    // all adjacent elements.
    for (int i = 0; i < numNodes; i++) {
      for (int j = nodeElementPtr[i]; j < nodeElementPtr[i + 1]; j++) {
        int elem = nodeToElements[j];
        rowp[i + 1] += elementNodeIndex[elem + 1] - elementNodeIndex[elem];
      }
    }

    // Make a conservative estimate of rowp
    for (int i = 0; i < numNodes; i++) {
      rowp[i + 1] += rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[numNodes];
    cols = new int[nnz];
    for (int i = 0; i < nnz; i++) {
      cols[i] = -1;
    }

    // Add the element contribution to the column indices
    for (int i = 0; i < numNodes; i++) {
      for (int jp = nodeElementPtr[i]; jp < nodeElementPtr[i + 1]; jp++) {
        int elem = nodeToElements[jp];

        // Add the columns to this row of the sparse matrix
        int row = rowp[i];
        int kend = elementNodeIndex[elem + 1];
        for (int kp = elementNodeIndex[elem]; kp < kend; kp++) {
          int node = elementTacsNodes[kp];
          node = getLocalNodeNum(node);
          cols[row] = node;
          row++;
        }
        rowp[i] = row;
      }
    }

    // Adjust rowp back to a zero-based index
    for (int i = numNodes; i > 0; i--) {
      rowp[i] = rowp[i - 1];
    }
    rowp[0] = 0;
  }

  // Go through and sort/uniquify each row and remove
  // the diagonal if requested
  TacsSortAndUniquifyCSR(numNodes, rowp, cols, nodiag);

  delete[] nodeElementPtr;
  delete[] nodeToElements;

  *_rowp = rowp;
  *_cols = cols;
}

/**
  Prepare a reduced CSR data structure corresponding to a matrix
  formed from a selection of the global matrix. This routine can be
  used in matrix/variable re-ordering computations.

  This function uses the same algorithm as computeLocalNodeToNodeCSR,
  but performs extra operations required to restrict the computations
  to Ar.  The rnodes array must consist of nrnodes non-negative
  integers between 0 and nrnodes-1, at any arbitrary location. All
  remaining entries of rnodes must be negative.

  @param rowp The row pointer corresponding to CSR data structure
  @param cols The column indices for each row of the CSR data structure
  @param nrnodes The number of reduced nodes
  @param rnodes The indices of the reduced nodes
  @param nodiag Flag to indicate whether to remove the diagonal matrix entry
*/
void TACSAssembler::computeLocalNodeToNodeCSR(int **_rowp, int **_cols,
                                              int nrnodes, const int *rnodes,
                                              int nodiag) {
  int *cols = NULL;
  int *rowp = new int[nrnodes + 1];
  memset(rowp, 0, (nrnodes + 1) * sizeof(int));

  // Create/get the node -> element data structure
  int *nodeElementPtr = NULL;
  int *nodeToElements = NULL;
  computeNodeToElementCSR(&nodeElementPtr, &nodeToElements);

  if (depNodes) {
    const int *depNodePtr = NULL;
    const int *depNodeConn = NULL;
    depNodes->getDepNodes(&depNodePtr, &depNodeConn, NULL);

    // Count the number of nodes associated with each element
    int *nodeCount = new int[numElements];
    memset(nodeCount, 0, numElements * sizeof(int));

    for (int i = 0; i < numElements; i++) {
      int jend = elementNodeIndex[i + 1];
      for (int j = elementNodeIndex[i]; j < jend; j++) {
        int node = elementTacsNodes[j];

        if (node >= 0) {
          // Convert to the local node number
          node = getLocalNodeNum(node);
          if (rnodes[node] >= 0) {
            nodeCount[i]++;
          }
        } else {
          // Find the dependent node
          int dep = -node - 1;
          for (int k = depNodePtr[dep]; k < depNodePtr[dep + 1]; k++) {
            node = getLocalNodeNum(depNodeConn[k]);
            if (rnodes[node] >= 0) {
              nodeCount[i]++;
            }
          }
        }
      }
    }

    // Count up the contribution to the rowp array from all elements
    // using the node->element data
    for (int i = 0; i < numNodes; i++) {
      int node = rnodes[i];
      if (node >= 0) {
        for (int j = nodeElementPtr[i]; j < nodeElementPtr[i + 1]; j++) {
          int elem = nodeToElements[j];
          rowp[node + 1] += nodeCount[elem];
        }
      }
    }

    // Make a conservative estimate of rowp
    for (int i = 0; i < nrnodes; i++) {
      rowp[i + 1] = rowp[i + 1] + rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[nrnodes];
    cols = new int[nnz];

    // Add the element contribution to the column indices
    for (int i = 0; i < numNodes; i++) {
      int node = rnodes[i];
      if (node >= 0) {
        for (int j = nodeElementPtr[i]; j < nodeElementPtr[i + 1]; j++) {
          int elem = nodeToElements[j];
          int kend = elementNodeIndex[elem + 1];

          // Scan through all the nodes belonging to this element
          int row = rowp[node];
          for (int k = elementNodeIndex[elem]; k < kend; k++) {
            int local = elementTacsNodes[k];
            if (local >= 0) {
              // Get the local node number
              local = getLocalNodeNum(local);
              int rn = rnodes[local];
              // This is an independent node
              if (rn >= 0) {
                cols[row] = rn;
                row++;
              }
            } else {
              // This is a dependent node, add the dependent node
              // variables
              int dep = -local - 1;
              int pend = depNodePtr[dep + 1];
              for (int p = depNodePtr[dep]; p < pend; p++) {
                local = depNodeConn[p];
                local = getLocalNodeNum(local);
                int rn = rnodes[local];
                if (rn >= 0) {
                  cols[row] = rn;
                  row++;
                }
              }
            }
          }
          rowp[node] = row;
        }
      }
    }

    // Adjust rowp back to a zero-based index
    for (int i = nrnodes; i > 0; i--) {
      rowp[i] = rowp[i - 1];
    }
    rowp[0] = 0;

    delete[] nodeCount;
  } else {
    // First, populate rowp by adding the contribution to a node from
    // all adjacent elements.
    for (int i = 0; i < numNodes; i++) {
      int node = rnodes[i];
      if (node >= 0) {
        for (int j = nodeElementPtr[i]; j < nodeElementPtr[i + 1]; j++) {
          int elem = nodeToElements[j];
          // Count up the number of reduced nodes that are required here.
          int count = 0;
          for (int k = elementNodeIndex[elem]; k < elementNodeIndex[elem + 1];
               k++) {
            int local = elementTacsNodes[k];
            local = getLocalNodeNum(local);
            if (rnodes[local] >= 0) {
              count++;
            }
          }

          rowp[node + 1] += count;
        }
      }
    }

    // Make a conservative estimate of rowp
    for (int i = 0; i < nrnodes; i++) {
      rowp[i + 1] = rowp[i + 1] + rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[nrnodes];
    cols = new int[nnz];

    // Add the element contribution to the column indices
    for (int i = 0; i < numNodes; i++) {
      int node = rnodes[i];
      if (node >= 0) {
        for (int j = nodeElementPtr[i]; j < nodeElementPtr[i + 1]; j++) {
          int elem = nodeToElements[j];
          int row = rowp[node];
          for (int k = elementNodeIndex[elem]; k < elementNodeIndex[elem + 1];
               k++) {
            int local = elementTacsNodes[k];
            local = getLocalNodeNum(local);
            int rn = rnodes[local];
            if (rn >= 0) {
              cols[row] = rn;
              row++;
            }
          }
          rowp[node] = row;
        }
      }
    }

    // Adjust rowp back to a zero-based index
    for (int i = nrnodes; i > 0; i--) {
      rowp[i] = rowp[i - 1];
    }
    rowp[0] = 0;
  }

  // Go through and sort/uniquify each row and remove
  // the diagonal if requested
  TacsSortAndUniquifyCSR(nrnodes, rowp, cols, nodiag);

  // Free the node -> element data structure
  delete[] nodeToElements;
  delete[] nodeElementPtr;

  *_rowp = rowp;
  *_cols = cols;
}

/**
  Compute the local node numbers that correspond to the coupling
  nodes connected to elements on other processes.

  Sort the global node numbers. Match the intervals and send them off
  to the owning process. On the owner, scan through the arrays until
  all the local coupling nodes are found.

  @param couplingNodes Local node numbers of the coupling nodes
  @param extPtr Pointer into the external node array (may be NULL)
  @param extCount External node count (may be NULL)
  @param recvPtr Incoming external ptr from other processors (may be NULL)
  @param recvCount Incoming external node count (may be NULL)
  @param recvNodes The incoming nodes from other procs (may be NULL)
*/
int TACSAssembler::computeCouplingNodes(int **_couplingNodes, int **_extPtr,
                                        int **_extCount, int **_recvPtr,
                                        int **_recvCount, int **_recvNodes) {
  // Get the ownership range and match the intervals of ownership
  const int *ownerRange;
  nodeMap->getOwnerRange(&ownerRange);

  // Get the external node numbers
  const int *extNodes = tacsExtNodeNums;
  if (extDistIndices) {
    extDistIndices->getIndices(&extNodes);
  }

  // Match the intervals for the external node numbers
  int *extPtr = new int[mpiSize + 1];
  int *extCount = new int[mpiSize];
  TacsMatchIntervals(mpiSize, ownerRange, numExtNodes, extNodes, extPtr);

  // Send the nodes owned by other processors the information. First
  // count up how many will go to each process.
  for (int i = 0; i < mpiSize; i++) {
    extCount[i] = extPtr[i + 1] - extPtr[i];
    if (i == mpiRank) {
      extCount[i] = 0;
    }
  }

  int *recvCount = new int[mpiSize];
  int *recvPtr = new int[mpiSize + 1];
  MPI_Alltoall(extCount, 1, MPI_INT, recvCount, 1, MPI_INT, tacs_comm);

  // Now prepare to send the node numbers to the other processors
  recvPtr[0] = 0;
  for (int i = 0; i < mpiSize; i++) {
    recvPtr[i + 1] = recvPtr[i] + recvCount[i];
  }

  // Number of nodes that will be received from other procs
  int *recvNodes = new int[recvPtr[mpiSize]];
  MPI_Alltoallv((void *)extNodes, extCount, extPtr, MPI_INT, recvNodes,
                recvCount, recvPtr, MPI_INT, tacs_comm);

  // Sort the recv'd nodes
  int *recvNodesSorted = NULL;
  if (_recvNodes) {
    recvNodesSorted = new int[recvPtr[mpiSize]];
    memcpy(recvNodesSorted, recvNodes, recvPtr[mpiSize] * sizeof(int));
  } else {
    recvNodesSorted = recvNodes;
  }

  // Uniquely sort the recieved nodes
  int nrecv_unique = TacsUniqueSort(recvPtr[mpiSize], recvNodesSorted);

  int num_multipliers = 0;
  int *multipliers = NULL;

  if (numMultiplierNodes > 0) {
    // Count up the number of multiplier nodes (if any). These are not
    // all of the locally owned multipliers, only those that refer to
    // external node numbers on other processors.
    multipliers = new int[numMultiplierNodes];

    // Add in the multiplier values
    for (int i = 0; i < numElements; i++) {
      // Get the multiplier index
      int multiplier = elements[i]->getMultiplierIndex();
      if (multiplier >= 0) {
        // Get the local multiplier index
        int elem_ptr = elementNodeIndex[i];

        // Check if the multiplier is locally owned. If it is not, it is
        // already in the external nodes anyway.
        int mult_node = elementTacsNodes[elem_ptr + multiplier];
        if (mult_node >= ownerRange[mpiRank] &&
            mult_node < ownerRange[mpiRank + 1]) {
          // Get the element size
          int size = elementNodeIndex[i + 1] - elem_ptr;

          // Check if any of the nodes are actually external
          for (int j = 0; j < size; j++) {
            int node = elementTacsNodes[elem_ptr + j];

            // Add the multiplier node if it's an external node or or it
            // is one of the nodes referred to on other procs.
            if (node < ownerRange[mpiRank] || node >= ownerRange[mpiRank + 1]) {
              multipliers[num_multipliers] = mult_node;
              num_multipliers++;
              break;
            } else if (TacsSearchArray(node, nrecv_unique, recvNodesSorted)) {
              multipliers[num_multipliers] = mult_node;
              num_multipliers++;
              break;
            }
          }
        }
      }
    }

    // Uniquely sort the number of multiplier nodes
    num_multipliers = TacsUniqueSort(num_multipliers, multipliers);
  }

  // Count up the number of coupling nodes
  int numCouplingNodes = nrecv_unique + numExtNodes + num_multipliers;
  if (_couplingNodes) {
    int *couplingNodes = new int[numCouplingNodes];

    // Automatically add in the external node numbers
    int index = 0;
    for (int i = 0; i < extNodeOffset; i++, index++) {
      couplingNodes[index] = i;
    }

    // Add the coupling nodes received from other processors
    if (num_multipliers > 0) {
      for (int i = 0, j = 0; (i < nrecv_unique || j < num_multipliers);
           index++) {
        if (i < nrecv_unique && j < num_multipliers) {
          if (recvNodesSorted[i] < multipliers[j]) {
            couplingNodes[index] = getLocalNodeNum(recvNodesSorted[i]);
            i++;
          } else {
            couplingNodes[index] = getLocalNodeNum(multipliers[j]);
            j++;
          }
        } else if (i < nrecv_unique) {
          couplingNodes[index] = getLocalNodeNum(recvNodesSorted[i]);
          i++;
        } else if (j < num_multipliers) {
          couplingNodes[index] = getLocalNodeNum(multipliers[j]);
          j++;
        }
      }
    } else {
      for (int i = 0; i < nrecv_unique; index++, i++) {
        couplingNodes[index] = getLocalNodeNum(recvNodesSorted[i]);
      }
    }

    // Add in the remaining external nodes
    for (int i = extNodeOffset; i < numExtNodes; i++, index++) {
      couplingNodes[index] = numOwnedNodes + i;
    }

    *_couplingNodes = couplingNodes;
  }

  // Free the multiplier index
  if (multipliers) {
    delete[] multipliers;
  }

  if (_extPtr) {
    *_extPtr = extPtr;
  } else {
    delete[] extPtr;
  }
  if (_extCount) {
    *_extCount = extCount;
  } else {
    delete[] extCount;
  }
  if (_recvPtr) {
    *_recvPtr = recvPtr;
  } else {
    delete[] recvPtr;
  }
  if (_recvCount) {
    *_recvCount = recvCount;
  } else {
    delete[] recvCount;
  }
  if (_recvNodes) {
    *_recvNodes = recvNodes;
    delete[] recvNodesSorted;
  } else {
    delete[] recvNodes;
  }

  return numCouplingNodes;
}

/**
  Compute the elements that couple with other processors.

  Compute the coupling nodes and the node to element pointer
  CSR data structure. From these, collect all elements that "own" a
  node that is referred to from another process.

  @param _couplingElems The elements coupled to other processors
*/
int TACSAssembler::computeCouplingElements(int **_couplingElems) {
  // Compute the nodes that couple to other processors
  int *couplingNodes;
  int numCouplingNodes = computeCouplingNodes(&couplingNodes);

  // Compute the node->element data structure
  int *nodeElementPtr, *nodeToElements;
  computeNodeToElementCSR(&nodeElementPtr, &nodeToElements);

  // Determine the elements that contain a coupling node
  int numCouplingElems = 0;
  int *couplingElems = new int[numElements];

  // Loop over all the coupling nodes and add all the elements
  // touched by each coupling node
  for (int i = 0; i < numCouplingNodes; i++) {
    int cnode = couplingNodes[i];

    for (int j = nodeElementPtr[cnode]; j < nodeElementPtr[cnode + 1]; j++) {
      int elem = nodeToElements[j];
      numCouplingElems =
          TacsMergeSortedArrays(numCouplingElems, couplingElems, 1, &elem);
    }
  }

  // Free the data
  delete[] nodeElementPtr;
  delete[] nodeToElements;
  delete[] couplingNodes;

  *_couplingElems = couplingElems;
  return numCouplingElems;
}

/**
  The function initialize performs a number of synchronization
  tasks that prepare the finite-element model for use.

  tacsNodeNums[i] is the global node number for the local node number i

  Two objects are required:
  1. NodeMap is constructed with the block sizes of each
  node owned by this process

  2. VecDistribute is constructed so that it takes an array and
  distributes its values to a vector or takes the vector values and
  collects them into an array This requires a sorted array of global
  node numbers.

  @return Fail flag indicating if a failure occured
*/
int TACSAssembler::initialize() {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call initialize() more than once!\n", mpiRank);
    return 1;
  }
  if (numDependentNodes > 0 && !depNodes) {
    fprintf(stderr, "[%d] Error: Dependent nodes not defined\n", mpiRank);
    return 1;
  }
  if (!elements) {
    fprintf(stderr, "[%d] Error: Elements not defined\n", mpiRank);
    return 1;
  }
  if (!elementNodeIndex) {
    fprintf(stderr, "[%d] Error: Element connectivity not defined\n", mpiRank);
    return 1;
  }

  // If the external nodes have not been computed, compute them now...
  if (!tacsExtNodeNums) {
    computeExtNodes();
  }

  // Flag to indicate that we've initialized TACSAssembler -
  // the initialization can only be done once
  meshInitializedFlag = 1;

  // Set up data for any dependent nodes. Note that the minimum
  // number of independent nodes is set as the maximum number
  // of element node by default. This is required for addMatValues()
  // to have enough memory for TRANSPOSE matrix assembly.
  maxElementIndepNodes = maxElementNodes;
  if (numDependentNodes > 0) {
    const int *depNodePtr;
    depNodes->getDepNodes(&depNodePtr, NULL, NULL);

    // Compute the maximum number of independent nodes
    for (int i = 0; i < numElements; i++) {
      int jend = elementNodeIndex[i + 1];
      int nnodes = 0;
      for (int j = elementNodeIndex[i]; j < jend; j++) {
        if (elementTacsNodes[j] >= 0) {
          nnodes++;
        } else {
          int dep = -elementTacsNodes[j] - 1;
          nnodes += depNodePtr[dep + 1] - depNodePtr[dep];
        }
      }
      if (nnodes > maxElementIndepNodes) {
        maxElementIndepNodes = nnodes;
      }
    }
  }

  // Create the distribution between the local nodes and the global ones
  extDistIndices = new TACSBVecIndices(&tacsExtNodeNums, numExtNodes);
  tacsExtNodeNums = NULL;
  extDistIndices->incref();
  extDistIndices->setUpInverse();

  // Set up the external vector distribute object
  extDist = new TACSBVecDistribute(nodeMap, extDistIndices);
  extDist->incref();

  // Scatter the boundary conditions to the external nodes
  scatterExternalBCs(bcMap);
  scatterExternalBCs(bcInitMap);

  // Allocate the vectors
  varsVec = createVec();
  varsVec->incref();
  dvarsVec = createVec();
  dvarsVec->incref();
  ddvarsVec = createVec();
  ddvarsVec->incref();
  xptVec = createNodeVec();
  xptVec->incref();

  // Allocate memory for the working array:
  // Determine the size of the data working array
  // max requirement is 4 element variable-size arrays,
  // 2 node-size arrays and either the element matrix or
  // the derivative of the residuals w.r.t. the nodes.
  int dataSize = maxElementIndepNodes + 4 * maxElementSize +
                 2 * TACS_SPATIAL_DIM * maxElementNodes;
  if (TACS_SPATIAL_DIM * maxElementNodes > maxElementSize) {
    dataSize += TACS_SPATIAL_DIM * maxElementNodes * maxElementSize;
  } else {
    dataSize += maxElementSize * maxElementSize;
  }
  elementData = new TacsScalar[dataSize];

  int idataSize = maxElementIndepNodes + maxElementNodes + 1;
  elementIData = new int[idataSize];

  // Allocate memory for the design variable data
  elementSensData = new TacsScalar[designVarsPerNode * maxElementDesignVars];
  elementSensIData = new int[maxElementDesignVars];

  // Create the design variable node mapping
  if (!designNodeMap) {
    // Get the number of design variables
    int numDVs = getNumDesignVars();

    if (mpiRank > 0) {
      numDVs = 0;
    }

    designNodeMap = new TACSNodeMap(tacs_comm, numDVs);
    designNodeMap->incref();
  }

  // Find all the design variable numbers associated
  if (!designExtDist) {
    // Get the design variables from the auxiliary elements
    int dvLen = 0;
    const int maxDVs = maxElementDesignVars;
    int *dvNums = elementSensIData;

    // Get the design variable range
    const int *range;
    designNodeMap->getOwnerRange(&range);
    int lower = range[mpiRank];
    int upper = range[mpiRank + 1];

    if (auxElements) {
      TACSAuxElem *aux = NULL;
      int naux = auxElements->getAuxElements(&aux);
      for (int i = 0; i < naux; i++) {
        int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);
        for (int j = 0; j < numDVs; j++) {
          if (dvNums[j] < lower || dvNums[j] >= upper) {
            dvLen++;
          }
        }
      }
    }
    for (int i = 0; i < numElements; i++) {
      int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
      for (int j = 0; j < numDVs; j++) {
        if (dvNums[j] < lower || dvNums[j] >= upper) {
          dvLen++;
        }
      }
    }

    if (designDepNodes) {
      const int *dep_ptr, *dep_conn;
      int num_dep = designDepNodes->getDepNodes(&dep_ptr, &dep_conn, NULL);
      for (int j = 0; j < dep_ptr[num_dep]; j++) {
        if (dep_conn[j] < lower || dep_conn[j] >= upper) {
          dvLen++;
        }
      }
    }

    // Allocate space for absolutely everything!
    int *allDVs = new int[dvLen];

    // Add all of the design variable numbers to the external node values
    dvLen = 0;
    if (auxElements) {
      TACSAuxElem *aux = NULL;
      int naux = auxElements->getAuxElements(&aux);
      for (int i = 0; i < naux; i++) {
        int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);
        for (int j = 0; j < numDVs; j++) {
          if (dvNums[j] < lower || dvNums[j] >= upper) {
            allDVs[dvLen] = dvNums[j];
            dvLen++;
          }
        }
      }
    }
    for (int i = 0; i < numElements; i++) {
      int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
      for (int j = 0; j < numDVs; j++) {
        if (dvNums[j] < lower || dvNums[j] >= upper) {
          allDVs[dvLen] = dvNums[j];
          dvLen++;
        }
      }
    }

    if (designDepNodes) {
      const int *dep_ptr, *dep_conn;
      int num_dep = designDepNodes->getDepNodes(&dep_ptr, &dep_conn, NULL);
      for (int j = 0; j < dep_ptr[num_dep]; j++) {
        if (dep_conn[j] < lower || dep_conn[j] >= upper) {
          allDVs[dvLen] = dep_conn[j];
          dvLen++;
        }
      }
    }

    dvLen = TacsUniqueSort(dvLen, allDVs);
    int *dvs = new int[dvLen];
    memcpy(dvs, allDVs, dvLen * sizeof(int));
    delete[] allDVs;

    TACSBVecIndices *dvIndices = new TACSBVecIndices(&dvs, dvLen);
    dvIndices->setUpInverse();
    designExtDist = new TACSBVecDistribute(designNodeMap, dvIndices);
    designExtDist->incref();
  }

  return 0;
}

/*
  Check the element Jacobian entries at each quadrature point to
  see if they are positive. The code prints the rank and element
  index on that rank if the check fails.
*/
void TACSAssembler::checkElementDeterminants() {
  if (elements) {
    TacsScalar *vars, *dvars, *ddvars;
    TacsScalar *elemXpts;
    getDataPointers(elementData, &vars, &dvars, &ddvars, NULL, &elemXpts, NULL,
                    NULL, NULL);

    for (int elemIndex = 0; elemIndex < numElements; elemIndex++) {
      TACSElement *element = elements[elemIndex];

      // Determine the values of the state variables for the
      // current element
      int ptr = elementNodeIndex[elemIndex];
      int len = elementNodeIndex[elemIndex + 1] - ptr;
      const int *nodes = &elementTacsNodes[ptr];
      xptVec->getValues(len, nodes, elemXpts);
      varsVec->getValues(len, nodes, vars);
      dvarsVec->getValues(len, nodes, dvars);
      ddvarsVec->getValues(len, nodes, ddvars);

      for (int n = 0; n < element->getNumQuadraturePoints(); n++) {
        double pt[3];
        element->getQuadraturePoint(n, pt);

        // Evaluate the pointwise element density
        TacsScalar density = 0.0, detXd = 0.0;
        int count = element->evalPointQuantity(elemIndex, TACS_ELEMENT_DENSITY,
                                               time, n, pt, elemXpts, vars,
                                               dvars, ddvars, &detXd, &density);

        if (count > 0) {
          if (TacsRealPart(detXd) <= 0.0) {
            printf(
                "[%d] TACS Warning: Negative determinant of the Jacobian "
                "transformation for element %d of type %s. Flipped mesh?",
                mpiRank, elemIndex, element->getObjectName());
          }
        }
      }
    }
  }
}

/*
  Scatter the boundary conditions that are shared between processors

  Note that we do not need to scatter the values along with the
  boundary condition values because the values are only used on the
  processors which own the nodes (that already have the correct
  information.)

  This function is called during initialize()
*/
void TACSAssembler::scatterExternalBCs(TACSBcMap *bcs) {
  // Get the coupling nodes shared between processors
  int *extPtr, *extCount;
  int *recvPtr, *recvCount, *recvNodes;
  computeCouplingNodes(NULL, &extPtr, &extCount, &recvPtr, &recvCount,
                       &recvNodes);

  // Get the nodes/variables
  const int *nodes, *vars;
  int nbcs = bcs->getBCs(&nodes, &vars, NULL);

  // Allocate the maximum array size that will be required
  int max_recv_size = 100;
  int *recvBCs = new int[max_recv_size];
  int index = 0;

  // Search through the coupling nodes recvd from other procs
  recvPtr[0] = 0;
  int ptr = 0;
  for (int k = 0; k < mpiSize; k++) {
    // Utilize the fact that the recvNodes are sorted
    // from each processor
    if (recvCount[k] > 0) {
      int size = recvCount[k];
      for (int i = 0; i < nbcs; i++) {
        int *item = TacsSearchArray(nodes[i], size, &recvNodes[ptr]);

        // This node is an interface node and a boundary node
        // add it to the list
        if (item) {
          // Extend the array if required
          if (2 * (index + 1) >= max_recv_size) {
            max_recv_size *= 2;
            int *temp = new int[max_recv_size];
            memcpy(temp, recvBCs, 2 * index * sizeof(int));
            delete[] recvBCs;
            recvBCs = temp;
          }

          // Set the new values into the array
          recvBCs[2 * index] = nodes[i];
          recvBCs[2 * index + 1] = vars[i];
          index++;
        }
      }

      // Update the pointer into the BC array
      ptr += size;
    }

    // Record the number of newly added nodes
    recvPtr[k + 1] = 2 * index;
    recvCount[k] = 2 * index - recvPtr[k];
  }

  // Send the counts to the other procs
  MPI_Alltoall(recvCount, 1, MPI_INT, extCount, 1, MPI_INT, tacs_comm);

  // Count up the size
  extPtr[0] = 0;
  for (int i = 0; i < mpiSize; i++) {
    extPtr[i + 1] = extPtr[i] + extCount[i];
  }

  // Allocate an array for the incoming data
  int numExtBCs = extPtr[mpiSize];
  int *extBCs = new int[numExtBCs];
  MPI_Alltoallv(recvBCs, recvCount, recvPtr, MPI_INT, extBCs, extCount, extPtr,
                MPI_INT, tacs_comm);

  // Free the data that is no longer required
  delete[] recvBCs;
  delete[] extPtr;
  delete[] extCount;
  delete[] recvPtr;
  delete[] recvCount;
  delete[] recvNodes;

  // Add the external boundary conditions
  for (int k = 0; k < numExtBCs; k += 2) {
    bcs->addBinaryFlagBC(extBCs[k], extBCs[k + 1]);
  }
  delete[] extBCs;
}

/**
  Get pointers to the element data. This code provides a way to
  automatically segment an array to avoid coding mistakes.

  Note that this is coded in such a way that you can provide NULL
  arguments
*/
void TACSAssembler::getDataPointers(TacsScalar *data, TacsScalar **v1,
                                    TacsScalar **v2, TacsScalar **v3,
                                    TacsScalar **v4, TacsScalar **x1,
                                    TacsScalar **x2, TacsScalar **weights,
                                    TacsScalar **mat) {
  int s = 0;
  if (v1) {
    *v1 = &data[s];
    s += maxElementSize;
  }
  if (v2) {
    *v2 = &data[s];
    s += maxElementSize;
  }
  if (v3) {
    *v3 = &data[s];
    s += maxElementSize;
  }
  if (v4) {
    *v4 = &data[s];
    s += maxElementSize;
  }
  if (x1) {
    *x1 = &data[s];
    s += TACS_SPATIAL_DIM * maxElementNodes;
  };
  if (x2) {
    *x2 = &data[s];
    s += TACS_SPATIAL_DIM * maxElementNodes;
  };
  if (weights) {
    *weights = &data[s];
    s += maxElementIndepNodes;
  }
  if (mat) {
    *mat = &data[s];
  }
}

/**
  Check whether a reordering has been applied to the nodes

  @return Flag indicating if the nodes have been reordered
*/
int TACSAssembler::isReordered() { return (newNodeIndices != NULL); }

/**
  Get the ordering from the old nodes to the new nodes

  @param oldToNew Array of size equal to the number of owned nodes
*/
void TACSAssembler::getReordering(int *oldToNew) {
  if (newNodeIndices) {
    // Get the new node indices
    const int *newNodes;
    newNodeIndices->getIndices(&newNodes);
    memcpy(oldToNew, &newNodes[extNodeOffset], numOwnedNodes * sizeof(int));
  } else {
    const int *ownerRange;
    nodeMap->getOwnerRange(&ownerRange);
    for (int k = 0; k < numOwnedNodes; k++) {
      oldToNew[k] = ownerRange[mpiRank] + k;
    }
  }
}

/**
  Reorder the vector using the reordering computed using the
  computeReordering() call.

  This is useful for reordering nodal vectors after the reordering has
  been applied.

  @param vec The vector to be reordered in place
*/
void TACSAssembler::reorderVec(TACSBVec *vec) {
  if (newNodeIndices) {
    // Get the ownership range
    const int *ownerRange;
    nodeMap->getOwnerRange(&ownerRange);

    // Get the vector of values from the array
    TacsScalar *x;
    int bsize = vec->getBlockSize();
    int size = vec->getArray(&x);

    // Allocate an array to store the old values and fill them in
    TacsScalar *xold = new TacsScalar[size];
    memcpy(xold, x, size * sizeof(TacsScalar));

    // Get the new node indices
    const int *newNodes;
    newNodeIndices->getIndices(&newNodes);

    // Loop through the owned nodes
    for (int i = 0; i < numOwnedNodes; i++) {
      // Get the new node value
      int node = newNodes[extNodeOffset + i];
      node = node - ownerRange[mpiRank];

      // Copy the values back to the array in the new
      // order
      for (int k = 0; k < bsize; k++) {
        x[bsize * node + k] = xold[bsize * i + k];
      }
    }

    delete[] xold;
  }
}

/**
  Reorder the nodes from the initial ordering to the final ordering

  Note that this only works for the nodes that are owned locally.
  Non-local or local non-owned nodes will not be reordered properly.

  @param num_nodes The length of the nodes array
  @param nodes The owned node numbers that will be converted
*/
void TACSAssembler::reorderNodes(int num_nodes, int *nodes) {
  if (newNodeIndices) {
    // Get the ownership range
    const int *ownerRange;
    nodeMap->getOwnerRange(&ownerRange);

    // Get the new node indices
    const int *newNodes;
    newNodeIndices->getIndices(&newNodes);

    for (int i = 0; i < num_nodes; i++) {
      if (nodes[i] >= ownerRange[mpiRank] &&
          nodes[i] < ownerRange[mpiRank + 1]) {
        int index = (nodes[i] - ownerRange[mpiRank]) + extNodeOffset;
        nodes[i] = newNodes[index];
      } else {
        nodes[i] = -1;
      }
    }
  }
}

/**
  Create a distributed design vector

  Vector classes initialized by one TACS object, cannot be used by a
  second, unless they share are exactly the parallel layout.

  @return A new design variable vector
*/
TACSBVec *TACSAssembler::createDesignVec() {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call createDeignVec() before initialize()\n",
            mpiRank);
    return NULL;
  }

  // Create the vector
  return new TACSBVec(designNodeMap, designVarsPerNode, designExtDist,
                      designDepNodes);
}

/**
  Set the design variable mapping, indicating the owners of the design vars

  @param designVarsPerNode The number of design variables for each var number
  @param designVarMap The design variable mapping
*/
void TACSAssembler::setDesignNodeMap(int _designVarsPerNode,
                                     TACSNodeMap *_designNodeMap) {
  if (meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call setDesignVarMap() after initialize()\n",
            mpiRank);
    return;
  }

  // Copy over the data
  designVarsPerNode = _designVarsPerNode;
  if (_designNodeMap) {
    _designNodeMap->incref();
  }
  if (designNodeMap) {
    designNodeMap->decref();
  }
  designNodeMap = _designNodeMap;
}

/**
  Set the dependent design variable information

  @param numDepDesignVars The number of dependent design variables
  @param depNodePtr Pointer into the depNodes array
  @param depNodes The dependent node numbers
  @param depNodeWeights The weights applied to each independent node
  @return Error code indicating failure or success
*/
int TACSAssembler::setDesignDependentNodes(int numDepDesignVars,
                                           const int *_depNodePtr,
                                           const int *_depNodeConn,
                                           const double *_depNodeWeights) {
  if (meshInitializedFlag) {
    fprintf(stderr,
            "[%d] Cannot call setDesignDependentNodes() after initialize()\n",
            mpiRank);
    return 1;
  }

  // Free the data if the dependent nodes have already been set
  if (designDepNodes) {
    designDepNodes->decref();
    designDepNodes = NULL;
  }

  // Get the ownership range of the nodes
  if (designNodeMap) {
    const int *ownerRange;
    designNodeMap->getOwnerRange(&ownerRange);

    // Check that all the independent nodes are positive and are within an
    // allowable range
    for (int i = 0; i < numDepDesignVars; i++) {
      for (int jp = _depNodePtr[i]; jp < _depNodePtr[i + 1]; jp++) {
        if (_depNodeConn[jp] >= ownerRange[mpiSize]) {
          fprintf(stderr,
                  "[%d] Dependent design node %d contains node number "
                  "%d out of range\n",
                  mpiRank, i, _depNodeConn[jp]);
          return 1;
        } else if (_depNodeConn[jp] < 0) {
          fprintf(stderr,
                  "[%d] Dependent design node %d contains dependent node %d\n",
                  mpiRank, i, _depNodeConn[jp]);
          return 1;
        }
      }
    }
  }

  if (numDepDesignVars > 0) {
    // Allocate the new memory and copy over the data
    int *depNodePtr = new int[numDepDesignVars + 1];
    memcpy(depNodePtr, _depNodePtr, (numDepDesignVars + 1) * sizeof(int));

    int size = depNodePtr[numDepDesignVars];
    int *depNodeConn = new int[size];
    memcpy(depNodeConn, _depNodeConn, size * sizeof(int));

    double *depNodeWeights = new double[size];
    memcpy(depNodeWeights, _depNodeWeights, size * sizeof(double));

    // Allocate the dependent node data structure
    designDepNodes = new TACSBVecDepNodes(numDepDesignVars, &depNodePtr,
                                          &depNodeConn, &depNodeWeights);
    designDepNodes->incref();
  }

  return 0;
}

/**
  Get the number of design variables defined by this assembler object.

  Note that when the designNodeMap is defined, the number of design variables
  are taken from the defined mapping. Otherwise, the maximum design variable
  values is computed from all elements in the TACSAssembler object.

  @return The number of design variables defined by TACSAssembler
*/
int TACSAssembler::getNumDesignVars() {
  if (designNodeMap) {
    const int *range;
    designNodeMap->getOwnerRange(&range);
    return range[mpiRank];
  } else if (elementSensIData) {
    // Get the design variables from the elements on this process
    const int maxDVs = maxElementDesignVars;
    int *dvNums = elementSensIData;

    // Set the maximum design variable number
    int maxDV = 0;

    // Get the design variables from the auxiliary elements
    if (auxElements) {
      TACSAuxElem *aux = NULL;
      int naux = auxElements->getAuxElements(&aux);
      for (int i = 0; i < naux; i++) {
        int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);

        for (int j = 0; j < numDVs; j++) {
          if (dvNums[j] > maxDV) {
            maxDV = dvNums[j];
          }
        }
      }
    }

    for (int i = 0; i < numElements; i++) {
      int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
      for (int j = 0; j < numDVs; j++) {
        if (dvNums[j] > maxDV) {
          maxDV = dvNums[j];
        }
      }
    }

    int maxDVglobal = 0;
    MPI_Allreduce(&maxDV, &maxDVglobal, 1, MPI_INT, MPI_MAX, tacs_comm);

    return maxDVglobal + 1;
  }

  return 0;
}

/**
  Collect all the design variable values assigned by this process

  This code does not ensure consistency of the design variable values
  between processes. If the values of the design variables are
  inconsistent to begin with. Call setDesignVars to force consistency.

  Each process contains objects that maintain their own design
  variable values. Ensuring the consistency of the ordering is up to
  the user. Having multiply-defined design variable numbers
  corresponding to different design variables results in undefined
  behaviour.

  @param dvs The vector of design variables
*/
void TACSAssembler::getDesignVars(TACSBVec *dvs) {
  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *dvVals = elementSensData;
  int *dvNums = elementSensIData;

  // Get the design variables from the auxiliary elements
  if (auxElements) {
    TACSAuxElem *aux = NULL;
    int naux = auxElements->getAuxElements(&aux);
    for (int i = 0; i < naux; i++) {
      int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);
      aux[i].elem->getDesignVars(aux[i].num, maxDVs, dvVals);
      dvs->setValues(numDVs, dvNums, dvVals, TACS_INSERT_NONZERO_VALUES);
    }
  }

  for (int i = 0; i < numElements; i++) {
    int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
    elements[i]->getDesignVars(i, numDVs, dvVals);
    dvs->setValues(numDVs, dvNums, dvVals, TACS_INSERT_NONZERO_VALUES);
  }

  dvs->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  dvs->endSetValues(TACS_INSERT_NONZERO_VALUES);
}

/**
  Set the design variables.

  The design variable values provided must be the same on all
  processes for consistency. This must be called by all processors.

  @param dvs The design variable values
*/
void TACSAssembler::setDesignVars(TACSBVec *dvs) {
  // Distribute the non-local design variable values
  dvs->beginDistributeValues();
  dvs->endDistributeValues();

  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *dvVals = elementSensData;
  int *dvNums = elementSensIData;

  // Get the design variables from the auxiliary elements
  if (auxElements) {
    TACSAuxElem *aux = NULL;
    int naux = auxElements->getAuxElements(&aux);
    for (int i = 0; i < naux; i++) {
      int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);
      dvs->getValues(numDVs, dvNums, dvVals);
      aux[i].elem->setDesignVars(aux[i].num, maxDVs, dvVals);
    }
  }

  for (int i = 0; i < numElements; i++) {
    int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
    dvs->getValues(numDVs, dvNums, dvVals);
    elements[i]->setDesignVars(i, numDVs, dvVals);
  }
}

/**
  Retrieve the design variable range.

  This call is collective on all TACS processes. The ranges provided
  by indivdual objects may not be consistent (if someone provided
  incorrect data they could be.)

  @param lb the lower bound on the design variables
  @param ub the upper bound on the design variables
*/
void TACSAssembler::getDesignVarRange(TACSBVec *lb, TACSBVec *ub) {
  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *dvVals = elementSensData;
  TacsScalar *ubVals = new TacsScalar[designVarsPerNode * maxDVs];
  int *dvNums = elementSensIData;

  // Get the design variables from the auxiliary elements
  if (auxElements) {
    TACSAuxElem *aux = NULL;
    int naux = auxElements->getAuxElements(&aux);
    for (int i = 0; i < naux; i++) {
      int numDVs = aux[i].elem->getDesignVarNums(aux[i].num, maxDVs, dvNums);
      aux[i].elem->getDesignVarRange(aux[i].num, maxDVs, dvVals, ubVals);
      lb->setValues(numDVs, dvNums, dvVals, TACS_INSERT_NONZERO_VALUES);
      ub->setValues(numDVs, dvNums, ubVals, TACS_INSERT_NONZERO_VALUES);
    }
  }

  for (int i = 0; i < numElements; i++) {
    int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
    elements[i]->getDesignVarRange(i, numDVs, dvVals, ubVals);
    lb->setValues(numDVs, dvNums, dvVals, TACS_INSERT_NONZERO_VALUES);
    ub->setValues(numDVs, dvNums, ubVals, TACS_INSERT_NONZERO_VALUES);
  }

  // Insert the values into the arrays
  lb->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  ub->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  lb->endSetValues(TACS_INSERT_NONZERO_VALUES);
  ub->endSetValues(TACS_INSERT_NONZERO_VALUES);

  delete[] ubVals;
}

/**
  Create a distributed solution vector.

  Vector classes initialized by one TACS object, cannot be used by a
  second, unless they share are exactly the parallel layout.

  @return A new, empty solution vector initialized to zero
*/
TACSBVec *TACSAssembler::createVec() {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call createVec() before initialize()\n",
            mpiRank);
    return NULL;
  }

  // Create the vector
  return new TACSBVec(nodeMap, varsPerNode, extDist, depNodes);
}

/**
  Apply the boundary conditions associated with the regular boundary
  conditions. This zeros all entries in the vector associated with Dirichlet
  conditions.

  @param vec Apply the boundary conditions to this vector
*/
void TACSAssembler::applyBCs(TACSVec *vec) { vec->applyBCs(bcMap); }

/**
  Apply the boundary conditions to the matrix

  @param mat Apply the boundary conditions to the rows of this matrix
*/
void TACSAssembler::applyBCs(TACSMat *mat) { mat->applyBCs(bcMap); }

/**
  Apply the boundary conditions to the tranpose of the matrix

  @param mat Apply the boundary conditions to the columns of this matrix
*/
void TACSAssembler::applyTransposeBCs(TACSMat *mat) {
  mat->applyTransposeBCs(bcMap);
}

/**
  Set the Dirichlet boundary conditions into the vector. This differs
  from applyBCs since the boundary condition values are set (not zeroed).

  @param vec Set the boundary conditions values into this vector
*/
void TACSAssembler::setBCs(TACSVec *vec) { vec->setBCs(bcMap); }

/**
  Create a distributed matrix

  This matrix is distributed in block-rows. Each processor owns a
  local part of the matrix and an off-diagonal part which couples
  between processors.

  This code creates a local array of global indices that is used to
  determine the destination for each entry in the sparse matrix.  This
  TACSBVecIndices object is reused if any subsequent parMat objects
  are created.

  @return A new parallel matrix with zeroed entries
*/
TACSParallelMat *TACSAssembler::createMat() {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call createMat() before initialize()\n",
            mpiRank);
    return NULL;
  }

  // Create the parMat indices if they do not already exist
  if (!parMatIndices) {
    // Get the global node numbering
    int *indices = new int[numNodes];
    for (int i = 0; i < numNodes; i++) {
      indices[i] = getGlobalNodeNum(i);
    }

    parMatIndices = new TACSBVecIndices(&indices, numNodes);
    parMatIndices->incref();
    parMatIndices->setUpInverse();
  }

  // Compute the local connectivity
  int *rowp, *cols;
  computeLocalNodeToNodeCSR(&rowp, &cols);

  // Create the distributed matrix class
  TACSParallelMat *dmat = new TACSParallelMat(
      thread_info, nodeMap, varsPerNode, numNodes, rowp, cols, parMatIndices);

  // Free the local connectivity
  delete[] rowp;
  delete[] cols;

  // Return the resulting matrix object
  return dmat;
}

/**
  Compute the connectivity information for the multiplier nodes

  This goes through element-by-element and finds the local node
  numbers associated with the multipliers. These multipliers must be
  ordered last. A unique set of these multipliers is determined and
  the associated element variables that must be ordered before the
  multiplier are determined. This data is stored in a CSR-like data
  structure that is sorted and uniquified before it is returned.

  @param num_multipliers The number of multipliers
  @param multipliers The multiplier node numbers
  @param indep_ptr Pointer into the array of independent node numbers
  @param indep_nodes Independent nodes associated with each multiplier
*/
void TACSAssembler::computeMultiplierConn(int *_num_multipliers,
                                          int **_multipliers, int **_indep_ptr,
                                          int **_indep_nodes) {
  *_num_multipliers = 0;
  *_multipliers = NULL;
  *_indep_ptr = NULL;
  *_indep_nodes = NULL;

  if (numMultiplierNodes == 0) {
    return;
  }

  // Determine the local multiplier index - the local node number for
  // each multiplier node.
  int *multipliers = new int[numMultiplierNodes];
  int num_multipliers = 0;
  for (int i = 0; i < numElements; i++) {
    int multiplier = elements[i]->getMultiplierIndex();
    if (multiplier >= 0) {
      // Get the local multiplier index
      int elem_ptr = elementNodeIndex[i];

      // Get the new local index for the multiplier
      int mult_node = getLocalNodeNum(elementTacsNodes[elem_ptr + multiplier]);

      // Set the multiplier node number
      multipliers[num_multipliers] = mult_node;
      num_multipliers++;
    }
  }

  // Uniquely sort the list of multipliers to determine the number of
  // true multiplier nodes
  num_multipliers = TacsUniqueSort(num_multipliers, multipliers);

  // Create the array of offsets into the pointer array
  int *ptr = new int[num_multipliers + 1];
  memset(ptr, 0, (num_multipliers + 1) * sizeof(int));

  // Now count up again to determine the total size of the non-zero
  // pattern required to store all the multiplier information
  for (int i = 0; i < numElements; i++) {
    int multiplier = elements[i]->getMultiplierIndex();
    if (multiplier >= 0) {
      // Get the local multiplier index
      int elem_ptr = elementNodeIndex[i];
      int size = elementNodeIndex[i + 1] - elem_ptr;

      // Get the new local index for the multiplier
      int mult_node = getLocalNodeNum(elementTacsNodes[elem_ptr + multiplier]);

      // Find the local index for external nodes
      int *item = TacsSearchArray(mult_node, num_multipliers, multipliers);
      if (item) {
        int index = item - multipliers;

        // Add in the element size (minus the multiplier node itself)
        ptr[index + 1] += size - 1;
      }
    }
  }

  // Set up the pointer array as an offset into each row
  for (int i = 1; i <= num_multipliers; i++) {
    ptr[i] += ptr[i - 1];
  }

  // Allocate the array of nodes
  int *nodes = new int[ptr[num_multipliers]];

  // Add in the values
  for (int i = 0; i < numElements; i++) {
    int multiplier = elements[i]->getMultiplierIndex();
    if (multiplier >= 0) {
      // Get the local multiplier index
      int elem_ptr = elementNodeIndex[i];
      int size = elementNodeIndex[i + 1] - elem_ptr;

      // Get the new local index for the multiplier
      int mult_node = getLocalNodeNum(elementTacsNodes[elem_ptr + multiplier]);

      // Find the local index for external nodes
      int *item = TacsSearchArray(mult_node, num_multipliers, multipliers);
      if (item) {
        int index = item - multipliers;
        for (int j = 0; j < size; j++) {
          if (j != multiplier) {
            // Get the local node number associated with the
            // independent degrees of freedom at this node
            nodes[ptr[index]] = getLocalNodeNum(elementTacsNodes[elem_ptr + j]);
            ptr[index]++;
          }
        }
      }
    }
  }

  // Reset the pointer array
  int offset = 0;
  for (int i = 0; i <= num_multipliers; i++) {
    int tmp = ptr[i];
    ptr[i] = offset;
    offset = tmp;
  }

  // Sort and uniquify the CSR - this compresses the CSR data
  TacsSortAndUniquifyCSR(num_multipliers, ptr, nodes);

  // Set the data and return
  *_num_multipliers = num_multipliers;
  *_multipliers = multipliers;
  *_indep_ptr = ptr;
  *_indep_nodes = nodes;
}

/**
  Create a parallel matrix for finite-element analysis.

  On the first call, this computes a reordering with the scheme
  provided. On subsequent calls, the reordering scheme is reused so
  that all TACSSchurMats, created from the same TACSAssembler object
  have the same non-zero structure. This makes adding matrices together
  easier (which is required for eigenvalue computations.)

  The first step is to determine the coupling nodes. (For a serial
  case there are no coupling nodes, so this is very simple!)  Then,
  the nodes that are not coupled to other processes are determined.
  The coupling and non-coupling nodes are ordered separately.  The
  coupling nodes must be ordered at the end of the block, while the
  local nodes must be ordered first. This type of constraint is not
  usually imposed in matrix ordering routines, so here we use a
  kludge.  First, order all the nodes and determine the ordering of
  the coupling variables within the full set. Next, order the local
  nodes. This hopefully reduces the fill-ins required, although there
  is no firm proof to back that up.

  The results from the reordering are placed in a set of objects.  The
  matrix reordering is stored in schurBIndices and schurCIndices while
  two mapping objects are created that map the variables from the
  global vector to reordered matrix.

  Mathematically this reordering can be written as follows,

  A' = (P A P^{T})

  where P^{T} is a permutation of the columns (variables), while P is
  a permutation of the rows (equations).

  @param order_type Order the Schur matrix with this type of ordeirng
  @return A new TACSSchurMat matrix with zeroed entries
*/
TACSSchurMat *TACSAssembler::createSchurMat(OrderingType order_type) {
  if (!meshInitializedFlag) {
    fprintf(stderr,
            "[%d] Cannot call createSchurMat() before "
            "initialize()\n",
            mpiRank);
    return NULL;
  }

  if (!schurBMap) {
    // The number of local nodes and the number of coupling nodes
    // that are referenced by other processors
    int nlocal_nodes = numNodes;
    int ncoupling_nodes = 0;

    // The local nodes and their
    int *perm_local_nodes = NULL;
    int *tacs_local_nodes = NULL;
    int *perm_coupling_nodes = NULL;
    int *tacs_coupling_nodes = NULL;

    if (order_type == TACS_AMD_ORDER) {
      // Use the AMD ordering scheme in TACS to compute an ordering of
      // the nodes that reduces the fill-in in the complete matrix --
      // including the off-diagonal contributions. This can reduce the
      // computational effort and the discrepancy between factorization
      // times on different processes.

      // Find the local node numbers of the coupling nodes.
      // Note that this is a sorted list
      int *coupling_nodes;
      ncoupling_nodes = computeCouplingNodes(&coupling_nodes);
      nlocal_nodes = numNodes - ncoupling_nodes;

      // Compute the CSR data structure for the node-to-node
      // connectivity without the diagonal entry
      int *rowp, *cols;
      int no_diagonal = 1;
      computeLocalNodeToNodeCSR(&rowp, &cols, no_diagonal);

      // Compute the multipliers/connectivity
      int num_multipliers = 0;
      int *multipliers = NULL;
      int *indep_ptr = NULL;
      int *indep_nodes = NULL;
      computeMultiplierConn(&num_multipliers, &multipliers, &indep_ptr,
                            &indep_nodes);

      // Here perm is the entire permutation array
      int *perm = new int[numNodes];
      int use_exact_degree = 0;  // Don't use the exact degree
      amd_order_interface(numNodes, rowp, cols, perm, coupling_nodes,
                          ncoupling_nodes, num_multipliers, multipliers,
                          indep_ptr, indep_nodes, use_exact_degree);

      // Free the rowp/cols array (which are modified by the
      // reordering anyway)
      delete[] rowp;
      delete[] cols;

      // Free the multipliers array
      delete[] multipliers;
      delete[] indep_ptr;
      delete[] indep_nodes;

      // Compute the coupling nodes based on their permutation
      tacs_local_nodes = new int[nlocal_nodes];
      tacs_coupling_nodes = new int[ncoupling_nodes];
      perm_local_nodes = new int[nlocal_nodes];
      perm_coupling_nodes = new int[ncoupling_nodes];
      for (int i = 0; i < nlocal_nodes; i++) {
        perm_local_nodes[i] = perm[i];
        tacs_local_nodes[i] = getGlobalNodeNum(perm_local_nodes[i]);
      }
      for (int i = 0; i < ncoupling_nodes; i++) {
        perm_coupling_nodes[i] = perm[i + nlocal_nodes];
        tacs_coupling_nodes[i] = getGlobalNodeNum(perm_coupling_nodes[i]);
      }

      delete[] perm;
      delete[] coupling_nodes;
    } else {
      // This scheme uses AMD or Nested Disection on the local and
      // coupling nodes independently. This ignores the off-diagonal
      // fill-ins which can be considerable!
      int no_diagonal = 0;
      if (order_type == ND_ORDER) {
        no_diagonal = 1;
      }

      // Find the local node numbers of the coupling nodes.
      // Note that this is a sorted list
      int *coupling_nodes;
      ncoupling_nodes = computeCouplingNodes(&coupling_nodes);
      nlocal_nodes = numNodes - ncoupling_nodes;

      // Set the coupling nodes for ordering
      int *all_nodes = new int[numNodes];
      for (int k = 0; k < numNodes; k++) {
        all_nodes[k] = -1;
      }

      for (int k = 0; k < ncoupling_nodes; k++) {
        all_nodes[coupling_nodes[k]] = k;
      }

      perm_coupling_nodes = new int[ncoupling_nodes];
      tacs_coupling_nodes = new int[ncoupling_nodes];

      // Now, compute the reordering for the local coupling variables
      if (ncoupling_nodes > 0) {
        int *rowp, *cols;
        computeLocalNodeToNodeCSR(&rowp, &cols, ncoupling_nodes, all_nodes,
                                  no_diagonal);

        // Compute the permutation of the coupling nodes
        computeMatReordering(order_type, ncoupling_nodes, rowp, cols,
                             perm_coupling_nodes, NULL);

        for (int i = 0; i < ncoupling_nodes; i++) {
          // Permute the coupling_nodes array - store in perm_coupling_nodes
          perm_coupling_nodes[i] = coupling_nodes[perm_coupling_nodes[i]];
          tacs_coupling_nodes[i] = getGlobalNodeNum(perm_coupling_nodes[i]);
        }

        delete[] rowp;
        delete[] cols;
      }

      // Set the remaining, local nodes for coupling
      perm_local_nodes = new int[nlocal_nodes];
      tacs_local_nodes = new int[nlocal_nodes];
      int *local_nodes = new int[nlocal_nodes];
      for (int j = 0, k = 0; k < numNodes; k++) {
        if (all_nodes[k] < 0) {
          all_nodes[k] = j;
          local_nodes[j] = k;
          j++;
        } else {
          all_nodes[k] = -1;
        }
      }

      // Now, compute the reordering for the local variables
      int *rowp, *cols;
      computeLocalNodeToNodeCSR(&rowp, &cols, nlocal_nodes, all_nodes,
                                no_diagonal);
      computeMatReordering(order_type, nlocal_nodes, rowp, cols,
                           perm_local_nodes, NULL);

      for (int i = 0; i < nlocal_nodes; i++) {
        // Permute the local nodes and record the corresponding tacs variables
        perm_local_nodes[i] = local_nodes[perm_local_nodes[i]];
        tacs_local_nodes[i] = getGlobalNodeNum(perm_local_nodes[i]);
      }

      delete[] rowp;
      delete[] cols;
      delete[] coupling_nodes;
      delete[] all_nodes;
      delete[] local_nodes;
    }

    // Create persistent objects so that all further TACSSchurMat will have
    // the same ordering.
    schurBIndices = new TACSBVecIndices(&perm_local_nodes, nlocal_nodes);
    schurCIndices = new TACSBVecIndices(&perm_coupling_nodes, ncoupling_nodes);
    schurBIndices->incref();
    schurCIndices->incref();

    TACSBVecIndices *tlocal =
        new TACSBVecIndices(&tacs_local_nodes, nlocal_nodes);
    TACSBVecIndices *tcoupling =
        new TACSBVecIndices(&tacs_coupling_nodes, ncoupling_nodes);
    schurBMap = new TACSBVecDistribute(nodeMap, tlocal);
    schurCMap = new TACSBVecDistribute(nodeMap, tcoupling);
    schurBMap->incref();
    schurCMap->incref();
  }

  // Compute he local non-zero pattern
  int *rowp, *cols;
  computeLocalNodeToNodeCSR(&rowp, &cols);

  TACSSchurMat *mat =
      new TACSSchurMat(thread_info, nodeMap, varsPerNode, numNodes, rowp, cols,
                       schurBIndices, schurBMap, schurCIndices, schurCMap);
  delete[] rowp;
  delete[] cols;

  return mat;
}

/*
  Create a serial BCSC matrix that enables sparse factorization with
  partial (row) pivoting. This matrix class is useful for debugging and
  testing the effects of roundoff errors in the solution process.  It
  does not work in parallel applications.
*/
TACSSerialPivotMat *TACSAssembler::createSerialMat() {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call createSerialMat() before initialize()\n",
            mpiRank);
    return NULL;
  }
  if (mpiSize > 1) {
    fprintf(stderr,
            "[%d] Cannot call createSerialMat() with multiple processors\n",
            mpiRank);
    return NULL;
  }

  // Compute he local non-zero pattern
  int *rowp, *cols;
  computeLocalNodeToNodeCSR(&rowp, &cols);

  // Allocate the matrix
  TACSSerialPivotMat *mat = new TACSSerialPivotMat(
      nodeMap, varsPerNode, numNodes, numNodes, rowp, cols);
  delete[] rowp;
  delete[] cols;

  return mat;
}

/**
  Retrieve the initial conditions associated with the problem

  @param vars The initial variable values (may be NULL)
  @param dvars The initial time derivative values (may be NULL)
  @param ddvars The initial second time derivative values (may be NULL)
*/
void TACSAssembler::getInitConditions(TACSBVec *vars, TACSBVec *dvars,
                                      TACSBVec *ddvars) {
  // If the user requests initial conditions of vars and vars0 has not been set,
  // zero the vector to prepare to get initial conditions from the element
  // interface
  if (vars && !vars0) {
    vars->zeroEntries();
  }
  if (dvars && !dvars0) {
    dvars->zeroEntries();
  }
  if (ddvars && !ddvars0) {
    ddvars->zeroEntries();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemDVars, *elemDDVars, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemDVars, &elemDDVars, NULL,
                  &elemXpts, NULL, NULL, NULL);

  // Retrieve the initial condition values from each element
  if ((vars && !vars0) || (dvars && !dvars0) || (ddvars && !ddvars0)) {
    for (int i = 0; i < numElements; i++) {
      int ptr = elementNodeIndex[i];
      int len = elementNodeIndex[i + 1] - ptr;
      const int *nodes = &elementTacsNodes[ptr];
      xptVec->getValues(len, nodes, elemXpts);

      // Get the initial condition values
      int nvars = elements[i]->getNumVariables();
      memset(elemVars, 0, nvars * sizeof(TacsScalar));
      memset(elemDVars, 0, nvars * sizeof(TacsScalar));
      memset(elemDDVars, 0, nvars * sizeof(TacsScalar));
      elements[i]->getInitConditions(i, elemXpts, elemVars, elemDVars,
                                     elemDDVars);

      // Set the values into the vectors
      if (vars && !vars0) {
        vars->setValues(len, nodes, elemVars, TACS_INSERT_NONZERO_VALUES);
      }
      if (dvars && !dvars0) {
        dvars->setValues(len, nodes, elemDVars, TACS_INSERT_NONZERO_VALUES);
      }
      if (ddvars && !ddvars0) {
        ddvars->setValues(len, nodes, elemDDVars, TACS_INSERT_NONZERO_VALUES);
      }
    }
  }

  // If the the initial conditions of vars is requested and the vars0 vector has
  // not been declared, then set the values of the tacs vector
  if (vars && !vars0) {
    vars->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  }
  if (dvars && !dvars0) {
    dvars->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  }
  if (ddvars && !ddvars0) {
    ddvars->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  }
  if (vars && !vars0) {
    vars->endSetValues(TACS_INSERT_NONZERO_VALUES);
  }
  if (dvars && !dvars0) {
    dvars->endSetValues(TACS_INSERT_NONZERO_VALUES);
  }
  if (ddvars && !ddvars0) {
    ddvars->beginSetValues(TACS_INSERT_NONZERO_VALUES);
  }

  // If the vars is requested and vars0 has been set, then copy the values from
  // vars0
  if (vars && vars0) {
    vars->copyValues(vars0);
  }
  if (dvars && dvars0) {
    dvars->copyValues(dvars0);
  }
  if (ddvars && ddvars0) {
    ddvars->copyValues(ddvars0);
  }
}

/**
Set the initial conditions for the problem

  @param vars The initial variable values (may be NULL)
  @param dvars The initial time derivative values (may be NULL)
  @param ddvars The initial second time derivative values (may be NULL)
*/
void TACSAssembler::setInitConditions(TACSBVec *vars, TACSBVec *dvars,
                                      TACSBVec *ddvars) {
  // Copy the values to the array.
  if (vars) {
    if (!vars0) {
      vars0 = createVec();
      vars0->incref();
    }
    vars0->copyValues(vars);
  }
  if (dvars) {
    if (!dvars0) {
      dvars0 = createVec();
      dvars0->incref();
    }
    dvars0->copyValues(dvars);
  }
  if (ddvars) {
    if (!ddvars0) {
      ddvars0 = createVec();
      ddvars0->incref();
    }
    ddvars0->copyValues(ddvars);
  }
}

/**
  Zero the entries of the local variables
*/
void TACSAssembler::zeroVariables() { varsVec->zeroEntries(); }

/**
  Zero the values of the time-derivatives of the state variables.
  This time-derivative is load-case independent.
*/
void TACSAssembler::zeroDotVariables() { dvarsVec->zeroEntries(); }

/**
  Zero the values of the time-derivatives of the state variables.
  This time-derivative is load-case independent.
*/
void TACSAssembler::zeroDDotVariables() { ddvarsVec->zeroEntries(); }

/**
  Set the value of the time/variables/time derivatives simultaneously

  @param vars The variable values (may be NULL)
  @param dvars The time derivative values (may be NULL)
  @param ddvars The second time derivative values (may be NULL)
*/
void TACSAssembler::setVariables(TACSBVec *vars, TACSBVec *dvars,
                                 TACSBVec *ddvars) {
  // Copy the values to the array.
  if (vars) {
    varsVec->copyValues(vars);
  }
  if (dvars) {
    dvarsVec->copyValues(dvars);
  }
  if (ddvars) {
    ddvarsVec->copyValues(ddvars);
  }

  // Distribute the values and evaluate the dependent nodes.
  if (vars) {
    varsVec->beginDistributeValues();
  }
  if (dvars) {
    dvarsVec->beginDistributeValues();
  }
  if (ddvars) {
    ddvarsVec->beginDistributeValues();
  }
  if (vars) {
    varsVec->endDistributeValues();
  }
  if (dvars) {
    dvarsVec->endDistributeValues();
  }
  if (ddvars) {
    ddvarsVec->endDistributeValues();
  }
}

/**
  Copy the values directly without distributing them.

  This is designed for specific instances where you want to use the
  values in the vector directly without distributing them.  This means
  that the variable values may be inconsistent across processors.
  This is not usually what you want to do.

  @param vars The variable values (may be NULL)
  @param dvars The time derivative values (may be NULL)
  @param ddvars The second time derivative values (may be NULL)
*/
void TACSAssembler::copyVariables(TACSBVec *vars, TACSBVec *dvars,
                                  TACSBVec *ddvars) {
  if (vars) {
    varsVec->copyValues(vars);
  }
  if (dvars) {
    dvarsVec->copyValues(dvars);
  }
  if (ddvars) {
    ddvarsVec->copyValues(ddvars);
  }
}

/**
  Get the variables from the vectors in TACS

  @param vars The variable values (may be NULL)
  @param dvars The time derivative values (may be NULL)
  @param ddvars The second time derivative values (may be NULL)
*/
void TACSAssembler::getVariables(TACSBVec *vars, TACSBVec *dvars,
                                 TACSBVec *ddvars) {
  // Copy the values to the array. Only local values are
  // copied, not external/dependents
  if (vars) {
    vars->copyValues(varsVec);
  }
  if (dvars) {
    dvars->copyValues(dvarsVec);
  }
  if (ddvars) {
    ddvars->copyValues(ddvarsVec);
  }
}

/**
  Get the variables from the vectors in TACS

  @param vars The variable values (may be NULL)
  @param dvars The time derivative values (may be NULL)
  @param ddvars The second time derivative values (may be NULL)
*/
void TACSAssembler::getVariables(TACSBVec **vars, TACSBVec **dvars,
                                 TACSBVec **ddvars) {
  // Copy the values to the array. Only local values are
  // copied, not external/dependents
  if (vars) {
    *vars = varsVec;
  }
  if (dvars) {
    *dvars = dvarsVec;
  }
  if (ddvars) {
    *ddvars = ddvarsVec;
  }
}

/**
  Set the simulation time internally in the TACSAssembler object

  @param _time The simulation time
*/
void TACSAssembler::setSimulationTime(double _time) { time = _time; }

/*
  Retrieve the simulation time from the TACSAssembler object

  @return The simulation time
*/
double TACSAssembler::getSimulationTime() { return time; }

/**
  Evaluates the total kinetic and potential energies of the structure

  @param Te The kinetic energy
  @param Pe The potential energy
*/
void TACSAssembler::evalEnergies(TacsScalar *Te, TacsScalar *Pe) {
  // Zero the kinetic and potential energy
  *Te = 0.0;
  *Pe = 0.0;

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, NULL, NULL, &elemXpts, NULL, NULL,
                  NULL);

  // Loop over all elements and add individual contributions to the
  // total energy
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);

    // Compute and add the element's contributions to the total
    // energy
    TacsScalar elemTe, elemPe;
    elements[i]->computeEnergies(i, time, elemXpts, vars, dvars, &elemTe,
                                 &elemPe);

    // Add up the kinetic and potential energy
    *Te += elemTe;
    *Pe += elemPe;
  }

  // Sum up the kinetic and potential energies across all processors
  TacsScalar input[2], output[2];
  input[0] = *Te;
  input[1] = *Pe;
  MPI_Allreduce(input, output, 2, TACS_MPI_TYPE, MPI_SUM, tacs_comm);

  *Te = output[0];
  *Pe = output[1];
}

/**
  Assemble the residual

  This residual includes the contributions from element tractions set
  in the auxiliary element classes. Note that the vector entries are
  zeroed first, and that the Dirichlet boundary conditions are applied
  after the assembly of the residual is complete.

  @param residual The residual vector
*/
void TACSAssembler::assembleRes(TACSBVec *residual) {
  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)
  if (auxElements) {
    auxElements->sort();
  }

  // Zero the residual
  residual->zeroEntries();

  if (thread_info->getNumThreads() > 1) {
    // Set the number of completed elements to zero
    numCompletedElements = 0;
    tacsPInfo->assembler = this;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&threads[k], &attr, TACSAssembler::assembleRes_thread,
                     (void *)tacsPInfo);
    }

    // Join all the threads
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  } else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
    getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                    NULL, NULL, NULL);

    // Get the auxiliary elements
    int naux = 0, aux_count = 0;
    TACSAuxElem *aux = NULL;
    if (auxElements) {
      naux = auxElements->getAuxElements(&aux);
    }

    // Go through and add the residuals from all the elements
    for (int i = 0; i < numElements; i++) {
      int ptr = elementNodeIndex[i];
      int len = elementNodeIndex[i + 1] - ptr;
      const int *nodes = &elementTacsNodes[ptr];
      xptVec->getValues(len, nodes, elemXpts);
      varsVec->getValues(len, nodes, vars);
      dvarsVec->getValues(len, nodes, dvars);
      ddvarsVec->getValues(len, nodes, ddvars);

      // Add the residual from the working element
      int nvars = elements[i]->getNumVariables();
      memset(elemRes, 0, nvars * sizeof(TacsScalar));
      elements[i]->addResidual(i, time, elemXpts, vars, dvars, ddvars, elemRes);

      // Add the residual from any auxiliary elements
      while (aux_count < naux && aux[aux_count].num == i) {
        aux[aux_count].elem->addResidual(i, time, elemXpts, vars, dvars, ddvars,
                                         elemRes);
        aux_count++;
      }

      // Add the residual values
      residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
    }
  }

  // Finish transmitting the residual
  residual->beginSetValues(TACS_ADD_VALUES);
  residual->endSetValues(TACS_ADD_VALUES);

  // Apply the boundary conditions for the residual
  residual->applyBCs(bcMap, varsVec);
}

/**
  Assemble the Jacobian matrix

  This function assembles the global Jacobian matrix and
  residual. This Jacobian includes the contributions from all
  elements. The Dirichlet boundary conditions are applied to the
  matrix by zeroing the rows of the matrix associated with a boundary
  condition, and setting the diagonal to unity. The matrix assembly
  also performs any communication required so that the matrix can be
  used immediately after assembly.

  @param alpha Coefficient for the variables
  @param beta Coefficient for the time-derivative terms
  @param gamma Coefficientfor the second time derivative term
  @param residual The residual of the governing equations
  @param A The Jacobian matrix
  @param matOr the matrix orientation NORMAL or TRANSPOSE
*/
void TACSAssembler::assembleJacobian(TacsScalar alpha, TacsScalar beta,
                                     TacsScalar gamma, TACSBVec *residual,
                                     TACSMat *A, MatrixOrientation matOr) {
  // Zero the residual and the matrix
  if (residual) {
    residual->zeroEntries();
  }
  A->zeroEntries();

  // Sort the list of auxiliary elements - this call only performs the
  // sort if it is required (if new elements are added)
  if (auxElements) {
    auxElements->sort();
  }

  // Run the p-threaded version of the assembly code
  if (thread_info->getNumThreads() > 1) {
    // Set the number of completed elements to zero
    numCompletedElements = 0;
    tacsPInfo->assembler = this;
    tacsPInfo->res = residual;
    tacsPInfo->mat = A;
    tacsPInfo->alpha = alpha;
    tacsPInfo->beta = beta;
    tacsPInfo->gamma = gamma;
    tacsPInfo->matOr = matOr;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&threads[k], &attr, TACSAssembler::assembleJacobian_thread,
                     (void *)tacsPInfo);
    }

    // Join all the threads
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  } else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
    TacsScalar *elemWeights, *elemMat;
    getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                    NULL, &elemWeights, &elemMat);

    // Set the data for the auxiliary elements - if there are any
    int naux = 0, aux_count = 0;
    TACSAuxElem *aux = NULL;
    if (auxElements) {
      naux = auxElements->getAuxElements(&aux);
    }

    for (int i = 0; i < numElements; i++) {
      int ptr = elementNodeIndex[i];
      int len = elementNodeIndex[i + 1] - ptr;
      const int *nodes = &elementTacsNodes[ptr];
      xptVec->getValues(len, nodes, elemXpts);
      varsVec->getValues(len, nodes, vars);
      dvarsVec->getValues(len, nodes, dvars);
      ddvarsVec->getValues(len, nodes, ddvars);

      // Get the number of variables from the element
      int nvars = elements[i]->getNumVariables();

      // Compute and add the contributions to the Jacobian
      memset(elemRes, 0, nvars * sizeof(TacsScalar));
      memset(elemMat, 0, nvars * nvars * sizeof(TacsScalar));
      elements[i]->addJacobian(i, time, alpha, beta, gamma, elemXpts, vars,
                               dvars, ddvars, elemRes, elemMat);

      // Add the contribution to the residual and the Jacobian
      // from the auxiliary elements - if any
      while (aux_count < naux && aux[aux_count].num == i) {
        aux[aux_count].elem->addJacobian(i, time, alpha, beta, gamma, elemXpts,
                                         vars, dvars, ddvars, elemRes, elemMat);
        aux_count++;
      }

      if (residual) {
        residual->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
      }
      addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
    }
  }

  // Do any matrix and residual assembly if required
  A->beginAssembly();
  if (residual) {
    residual->beginSetValues(TACS_ADD_VALUES);
  }

  A->endAssembly();
  if (residual) {
    residual->endSetValues(TACS_ADD_VALUES);
  }

  // Apply the boundary conditions
  if (residual) {
    residual->applyBCs(bcMap, varsVec);
  }

  // Apply the appropriate boundary conditions
  A->applyBCs(bcMap);
}

/**
  Assemble a matrix of a specified type. Note that all matrices
  created from the TACSAssembler object have the same non-zero pattern
  and are interchangable.

  @param matType The matrix type
  @param A The matrix to assemble (output)
  @param matOr The matrix orientation: NORMAL or TRANSPOSE
*/
void TACSAssembler::assembleMatType(ElementMatrixType matType, TACSMat *A,
                                    MatrixOrientation matOr) {
  // Zero the matrix
  A->zeroEntries();

  if (thread_info->getNumThreads() > 1) {
    // Set the number of completed elements to zero
    numCompletedElements = 0;
    tacsPInfo->assembler = this;
    tacsPInfo->mat = A;
    tacsPInfo->matType = matType;
    tacsPInfo->matOr = matOr;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_create(&threads[k], &attr, TACSAssembler::assembleMatType_thread,
                     (void *)tacsPInfo);
    }

    // Join all the threads
    for (int k = 0; k < thread_info->getNumThreads(); k++) {
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  } else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *elemXpts, *elemMat, *elemWeights;
    getDataPointers(elementData, &vars, NULL, NULL, NULL, &elemXpts, NULL,
                    &elemWeights, &elemMat);

    for (int i = 0; i < numElements; i++) {
      // Retrieve the element variables and node locations
      int ptr = elementNodeIndex[i];
      int len = elementNodeIndex[i + 1] - ptr;
      const int *nodes = &elementTacsNodes[ptr];
      xptVec->getValues(len, nodes, elemXpts);
      varsVec->getValues(len, nodes, vars);

      // Get the element matrix
      elements[i]->getMatType(matType, i, time, elemXpts, vars, elemMat);

      // Add the values into the element
      addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
    }
  }

  A->beginAssembly();
  A->endAssembly();

  A->applyBCs(bcMap);
}

/**
  Assemble a linear combination of matrices.

  This is used for some buckling/eigenvalue computations which require
  matrices that are linear combinations of specific types.

  @param matTypes The array of matrix types
  @param scale The array of scalar values
  @param nmats The number of matrices in the linear combination
  @param A The matrix to assemble (output)
  @param matOr the matrix orientation: NORMAL or TRANSPOSE
*/
void TACSAssembler::assembleMatCombo(ElementMatrixType matTypes[],
                                     TacsScalar scale[], int nmats, TACSMat *A,
                                     MatrixOrientation matOr) {
  // Zero the matrix
  A->zeroEntries();

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *elemXpts, *elemMat, *elemWeights;
  getDataPointers(elementData, &vars, NULL, NULL, NULL, &elemXpts, NULL,
                  &elemWeights, &elemMat);

  for (int i = 0; i < numElements; i++) {
    // Retrieve the element variables and node locations
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);

    for (int j = 0; j < nmats; j++) {
      // Get the element matrix
      elements[i]->getMatType(matTypes[j], i, time, elemXpts, vars, elemMat);

      // Scale the matrix
      if (scale[j] != 1.0) {
        int nvars = elements[i]->getNumVariables();
        int n = nvars * nvars, one = 1;
        BLASscal(&n, &scale[j], elemMat, &one);
      }

      // Add the values into the element
      addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
    }
  }

  A->beginAssembly();
  A->endAssembly();

  A->applyBCs(bcMap);
}

/**
  Evaluate a list of TACS functions

  First, check if the functions are initialized. Obtain the number of
  iterations over the function domain required to evaluate the
  functions.

  This function will print an error and return 0 if the underlying
  TACSAssembler object does not correspond to the TACSAssembler object.

  @param numFuncs The number of functions to evaluate
  @param funcs Array of functions to evaluate
  @param funcVals The function values
*/
void TACSAssembler::evalFunctions(int numFuncs, TACSFunction **funcs,
                                  TacsScalar *funcVals) {
  // Here we will use time-independent formulation
  TacsScalar tcoef = 1.0;

  // Check whether these are two-stage or single-stage functions
  int twoStage = 0;
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k] && funcs[k]->getStageType() == TACSFunction::TWO_STAGE) {
      twoStage = 1;
      break;
    }
  }

  // Just be lazy and call initialize on all the functions if any
  // function is two-stage
  if (twoStage) {
    for (int k = 0; k < numFuncs; k++) {
      if (funcs[k]) {
        funcs[k]->initEvaluation(TACSFunction::INITIALIZE);
      }
    }
    integrateFunctions(tcoef, TACSFunction::INITIALIZE, numFuncs, funcs);
    for (int k = 0; k < numFuncs; k++) {
      if (funcs[k]) {
        funcs[k]->finalEvaluation(TACSFunction::INITIALIZE);
      }
    }
  }

  // Perform the integration required to evaluate the function
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      funcs[k]->initEvaluation(TACSFunction::INTEGRATE);
    }
  }

  integrateFunctions(tcoef, TACSFunction::INTEGRATE, numFuncs, funcs);

  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      funcs[k]->finalEvaluation(TACSFunction::INTEGRATE);
    }
  }

  // Retrieve the function values
  for (int k = 0; k < numFuncs; k++) {
    funcVals[k] = 0.0;
    if (funcs[k]) {
      funcVals[k] = funcs[k]->getFunctionValue();
    }
  }
}

/**
  Integrate or initialize functions for a single time step of a time
  integration (or steady-state simulation).

  @param tcoef The integration coefficient
  @param ftype The type of integration to use
  @param numFuncs The number of functions
  @param funcs The array of functions
*/
void TACSAssembler::integrateFunctions(TacsScalar tcoef,
                                       TACSFunction::EvaluationType ftype,
                                       int numFuncs, TACSFunction **funcs) {
  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars;
  TacsScalar *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, NULL, &elemXpts, NULL,
                  NULL, NULL);

  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      if (funcs[k]->getDomainType() == TACSFunction::ENTIRE_DOMAIN) {
        for (int i = 0; i < numElements; i++) {
          // Determine the values of the state variables for the
          // current element
          int ptr = elementNodeIndex[i];
          int len = elementNodeIndex[i + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Evaluate the element-wise component of the function
          funcs[k]->elementWiseEval(ftype, i, elements[i], time, tcoef,
                                    elemXpts, vars, dvars, ddvars);
        }
      } else if (funcs[k]->getDomainType() == TACSFunction::SUB_DOMAIN) {
        const int *elementNums;
        int subDomainSize = funcs[k]->getElementNums(&elementNums);

        for (int i = 0; i < subDomainSize; i++) {
          int elemNum = elementNums[i];

          if (elemNum >= 0 && elemNum < numElements) {
            // Determine the values of the state variables
            // for the current element
            int ptr = elementNodeIndex[elemNum];
            int len = elementNodeIndex[elemNum + 1] - ptr;
            const int *nodes = &elementTacsNodes[ptr];
            xptVec->getValues(len, nodes, elemXpts);
            varsVec->getValues(len, nodes, vars);
            dvarsVec->getValues(len, nodes, dvars);
            ddvarsVec->getValues(len, nodes, ddvars);

            // Evaluate the element-wise component of the function
            funcs[k]->elementWiseEval(ftype, elemNum, elements[elemNum], time,
                                      tcoef, elemXpts, vars, dvars, ddvars);
          }
        }
      }
    }
  }
}

/**
  Evaluate the derivative of a list of functions w.r.t. the design
  variables.

  Note that a function should be evaluated - using evalFunction - before
  its derivatives can be evaluated.

  The design variable sensitivities are divided into two distinct
  sets: material-dependent design variables and shape design
  variables. The shape design variables are handled through the
  TACSNodeMap class. The material-dependent design variables are
  handled through the element classes themselves.

  In this code, the derivative of the function w.r.t. the
  shape-dependent design variables is handled first. The derivative of
  the function w.r.t each nodal location is determined. The
  TACSNodeMap object (if not NULL) is then used to determine the
  derivative of the nodal locations w.r.t. the design variables
  themselves.

  The material-dependent design variables are handled on an
  element-by-element and traction-by-traction dependent basis.

  Note that this function distributes the result to the processors
  through a collective communication call. No further parallel
  communication is required.

  @param coef The coefficient applied to the derivative
  @param numFuncs The number of functions - size of funcs array
  @param funcs The TACSFunction function objects
  @param dfdx The derivative vector
*/
void TACSAssembler::addDVSens(TacsScalar coef, int numFuncs,
                              TACSFunction **funcs, TACSBVec **dfdx) {
  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, NULL, &elemXpts, NULL,
                  NULL, NULL);

  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *fdvSens = elementSensData;
  int *dvNums = elementSensIData;

  // For each function, evaluate the derivative w.r.t. the
  // design variables for each element
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      if (funcs[k]->getDomainType() == TACSFunction::SUB_DOMAIN) {
        // Get the funcs[k] sub-domain
        const int *elemSubList;
        int numSubElems = funcs[k]->getElementNums(&elemSubList);

        for (int i = 0; i < numSubElems; i++) {
          int elemNum = elemSubList[i];
          // Determine the values of the state variables for subElem
          int ptr = elementNodeIndex[elemNum];
          int len = elementNodeIndex[elemNum + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Get the design variables for this element
          int numDVs =
              elements[elemNum]->getDesignVarNums(elemNum, maxDVs, dvNums);

          // Evaluate the element-wise sensitivity of the function
          memset(fdvSens, 0, numDVs * designVarsPerNode * sizeof(TacsScalar));
          funcs[k]->addElementDVSens(elemNum, elements[elemNum], time, coef,
                                     elemXpts, vars, dvars, ddvars, maxDVs,
                                     fdvSens);

          // Add the derivative values
          dfdx[k]->setValues(numDVs, dvNums, fdvSens, TACS_ADD_VALUES);
        }
      } else if (funcs[k]->getDomainType() == TACSFunction::ENTIRE_DOMAIN) {
        for (int elemNum = 0; elemNum < numElements; elemNum++) {
          // Determine the values of the state variables for elemNum
          int ptr = elementNodeIndex[elemNum];
          int len = elementNodeIndex[elemNum + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Get the design variables for this element
          int numDVs =
              elements[elemNum]->getDesignVarNums(elemNum, maxDVs, dvNums);

          // Evaluate the element-wise sensitivity of the function
          memset(fdvSens, 0, numDVs * designVarsPerNode * sizeof(TacsScalar));
          funcs[k]->addElementDVSens(elemNum, elements[elemNum], time, coef,
                                     elemXpts, vars, dvars, ddvars, maxDVs,
                                     fdvSens);

          // Add the derivative values
          dfdx[k]->setValues(numDVs, dvNums, fdvSens, TACS_ADD_VALUES);
        }
      }
    }
  }
}

/**
  Evaluate the derivative of the function w.r.t. the owned nodes.

  This code evaluates the sensitivity of the function w.r.t. the
  owned nodes for all elements in the function domain.

  Note that a function should be evaluated - using evalFunction - before
  its derivatives can be evaluated.

  This function should be preferred to the use of evalDVSens without a
  list of functions since it is more efficient!

  @param coef The coefficient applied to the derivative
  @param numFuncs The number of functions - size of funcs array
  @param funcs The TACSFunction function objects
  @param dfdXpt The derivative vector
*/
void TACSAssembler::addXptSens(TacsScalar coef, int numFuncs,
                               TACSFunction **funcs, TACSBVec **dfdXpt) {
  // First check if this is the right assembly object
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k] && this != funcs[k]->getAssembler()) {
      fprintf(stderr, "[%d] Cannot evaluate function %s, wrong TACS object\n",
              mpiRank, funcs[k]->getObjectName());
    }
  }

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars;
  TacsScalar *elemXpts, *elemXptSens;
  getDataPointers(elementData, &vars, &dvars, &ddvars, NULL, &elemXpts,
                  &elemXptSens, NULL, NULL);

  // For each function, evaluate the derivative w.r.t. the
  // nodal locations for all elements or part of the domain
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      if (funcs[k]->getDomainType() == TACSFunction::SUB_DOMAIN) {
        // Get the function sub-domain
        const int *elemSubList;
        int numSubElems = funcs[k]->getElementNums(&elemSubList);
        for (int i = 0; i < numSubElems; i++) {
          int elemNum = elemSubList[i];
          // Determine the values of the state variables for subElem
          int ptr = elementNodeIndex[elemNum];
          int len = elementNodeIndex[elemNum + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Evaluate the element-wise sensitivity of the function
          funcs[k]->getElementXptSens(elemNum, elements[elemNum], time, coef,
                                      elemXpts, vars, dvars, ddvars,
                                      elemXptSens);
          dfdXpt[k]->setValues(len, nodes, elemXptSens, TACS_ADD_VALUES);
        }
      } else if (funcs[k]->getDomainType() == TACSFunction::ENTIRE_DOMAIN) {
        for (int elemNum = 0; elemNum < numElements; elemNum++) {
          // Determine the values of the state variables for elemNum
          int ptr = elementNodeIndex[elemNum];
          int len = elementNodeIndex[elemNum + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Evaluate the element-wise sensitivity of the function
          funcs[k]->getElementXptSens(elemNum, elements[elemNum], time, coef,
                                      elemXpts, vars, dvars, ddvars,
                                      elemXptSens);
          dfdXpt[k]->setValues(len, nodes, elemXptSens, TACS_ADD_VALUES);
        }
      }
    }
  }
}

/**
  Evaluate the derivative of the function w.r.t. the state
  variables.

  This code evaluates the sensitivity of the function w.r.t. the
  state variables for all elements in the function domain. This
  code is usually much faster than the code for computing the
  derivative of the function w.r.t. the design variables.

  Note that the sensitivity vector 'vec' is assembled, and
  appropriate boundary conditions are imposed before the function
  is returned.

  @param alpha Coefficient for the variables
  @param beta Coefficient for the time-derivative terms
  @param gamma Coefficientfor the second time derivative term
  @param numFuncs The number of functions - size of funcs array
  @param funcs The TACSFunction function objects
  @param dfdu The derivative array
*/
void TACSAssembler::addSVSens(TacsScalar alpha, TacsScalar beta,
                              TacsScalar gamma, int numFuncs,
                              TACSFunction **funcs, TACSBVec **dfdu) {
  // First check if this is the right assembly object
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k] && this != funcs[k]->getAssembler()) {
      fprintf(stderr, "[%d] Cannot evaluate function %s, wrong TACS object\n",
              mpiRank, funcs[k]->getObjectName());
    }
  }

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemRes, &elemXpts,
                  NULL, NULL, NULL);

  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      if (funcs[k]->getDomainType() == TACSFunction::ENTIRE_DOMAIN) {
        for (int elemNum = 0; elemNum < numElements; elemNum++) {
          // Determine the values of the state variables for subElem
          int ptr = elementNodeIndex[elemNum];
          int len = elementNodeIndex[elemNum + 1] - ptr;
          const int *nodes = &elementTacsNodes[ptr];
          xptVec->getValues(len, nodes, elemXpts);
          varsVec->getValues(len, nodes, vars);
          dvarsVec->getValues(len, nodes, dvars);
          ddvarsVec->getValues(len, nodes, ddvars);

          // Evaluate the element-wise sensitivity of the function
          funcs[k]->getElementSVSens(elemNum, elements[elemNum], time, alpha,
                                     beta, gamma, elemXpts, vars, dvars, ddvars,
                                     elemRes);
          dfdu[k]->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
        }
      } else if (funcs[k]->getDomainType() == TACSFunction::SUB_DOMAIN) {
        const int *elementNums;
        int subDomainSize = funcs[k]->getElementNums(&elementNums);

        for (int i = 0; i < subDomainSize; i++) {
          int elemNum = elementNums[i];
          if (elemNum >= 0 && elemNum < numElements) {
            // Determine the values of the state variables for the
            // current element
            int ptr = elementNodeIndex[elemNum];
            int len = elementNodeIndex[elemNum + 1] - ptr;
            const int *nodes = &elementTacsNodes[ptr];
            xptVec->getValues(len, nodes, elemXpts);
            varsVec->getValues(len, nodes, vars);
            dvarsVec->getValues(len, nodes, dvars);
            ddvarsVec->getValues(len, nodes, ddvars);

            // Evaluate the element-wise sensitivity of the function
            funcs[k]->getElementSVSens(elemNum, elements[elemNum], time, alpha,
                                       beta, gamma, elemXpts, vars, dvars,
                                       ddvars, elemRes);
            dfdu[k]->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
          }
        }
      }

      // Add the values into the array
      dfdu[k]->beginSetValues(TACS_ADD_VALUES);
    }
  }

  // Finish adding the values
  for (int k = 0; k < numFuncs; k++) {
    if (funcs[k]) {
      dfdu[k]->endSetValues(TACS_ADD_VALUES);
      dfdu[k]->applyBCs(bcMap);
    }
  }
}

/**
  Evaluate the product of several ajdoint vectors with the derivative
  of the residual w.r.t. the design variables.

  This function is collective on all TACSAssembler processes. This
  computes the product of the derivative of the residual w.r.t. the
  design variables with several adjoint vectors simultaneously. This
  saves computational time as the derivative of the element residuals
  can be reused for each adjoint vector. This function performs the
  same task as evalAdjointResProduct, but uses more memory than
  calling it for each adjoint vector.

  @param scale Scalar factor applied to the derivative
  @param numAdjoints The number of adjoint vectors
  @param adjoint The array of adjoint vectors
  @param dfdx Product of the derivative of the residuals and the adjoint
*/
void TACSAssembler::addAdjointResProducts(TacsScalar scale, int numAdjoints,
                                          TACSBVec **adjoint, TACSBVec **dfdx) {
  // Distribute the design variable values to all processors
  for (int k = 0; k < numAdjoints; k++) {
    adjoint[k]->beginDistributeValues();
  }
  for (int k = 0; k < numAdjoints; k++) {
    adjoint[k]->endDistributeValues();
  }

  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)
  if (auxElements) {
    auxElements->sort();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars;
  TacsScalar *elemXpts, *elemAdjoint;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemAdjoint, &elemXpts,
                  NULL, NULL, NULL);

  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *fdvSens = elementSensData;
  int *dvNums = elementSensIData;

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  for (int i = 0; i < numElements; i++) {
    // Find the variables and nodes
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    // Get the design variables for this element
    int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);

    // Get the adjoint variables
    for (int k = 0; k < numAdjoints; k++) {
      memset(fdvSens, 0, numDVs * designVarsPerNode * sizeof(TacsScalar));

      // Get the element adjoint vector
      adjoint[k]->getValues(len, nodes, elemAdjoint);

      // Add the adjoint-residual product
      elements[i]->addAdjResProduct(i, time, scale, elemAdjoint, elemXpts, vars,
                                    dvars, ddvars, numDVs, fdvSens);

      dfdx[k]->setValues(numDVs, dvNums, fdvSens, TACS_ADD_VALUES);
    }

    // Add the contribution from the auxiliary elements
    if (aux_count < naux) {
      while (aux_count < naux && aux[aux_count].num == i) {
        // Get the design variables for this element
        numDVs = aux[aux_count].elem->getDesignVarNums(i, maxDVs, dvNums);

        // Get the adjoint variables
        for (int k = 0; k < numAdjoints; k++) {
          memset(fdvSens, 0, numDVs * designVarsPerNode * sizeof(TacsScalar));

          // Get the element adjoint vector
          adjoint[k]->getValues(len, nodes, elemAdjoint);

          aux[aux_count].elem->addAdjResProduct(i, time, scale, elemAdjoint,
                                                elemXpts, vars, dvars, ddvars,
                                                numDVs, fdvSens);

          dfdx[k]->setValues(numDVs, dvNums, fdvSens, TACS_ADD_VALUES);
        }
        aux_count++;
      }
    }
  }
}

/**
  Evaluate the product of several ajdoint vectors with the derivative
  of the residual w.r.t. the nodal points.

  This function is collective on all TACSAssembler processes. This
  computes the product of the derivative of the residual w.r.t. the
  nodal points with several adjoint vectors simultaneously. This
  saves computational time as the derivative of the element residuals
  can be reused for each adjoint vector.

  @param scale Scalar factor applied to the derivative
  @param numAdjoints The number of adjoint vectors
  @param adjoint The array of adjoint vectors
  @param dfdXpt Product of the derivative of the residuals and the adjoint
*/
void TACSAssembler::addAdjointResXptSensProducts(TacsScalar scale,
                                                 int numAdjoints,
                                                 TACSBVec **adjoint,
                                                 TACSBVec **dfdXpt) {
  for (int k = 0; k < numAdjoints; k++) {
    adjoint[k]->beginDistributeValues();
  }
  for (int k = 0; k < numAdjoints; k++) {
    adjoint[k]->endDistributeValues();
  }

  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)
  if (auxElements) {
    auxElements->sort();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars;
  TacsScalar *elemXpts, *elemAdjoint, *xptSens;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &elemAdjoint, &elemXpts,
                  &xptSens, NULL, NULL);

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  for (int i = 0; i < numElements; i++) {
    // Find the variables and nodes
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    // Get the adjoint variables
    for (int k = 0; k < numAdjoints; k++) {
      memset(xptSens, 0, TACS_SPATIAL_DIM * len * sizeof(TacsScalar));
      adjoint[k]->getValues(len, nodes, elemAdjoint);
      elements[i]->addAdjResXptProduct(i, time, scale, elemAdjoint, elemXpts,
                                       vars, dvars, ddvars, xptSens);

      dfdXpt[k]->setValues(len, nodes, xptSens, TACS_ADD_VALUES);
    }

    // Add the contribution from the auxiliary elements
    if (aux_count < naux) {
      while (aux_count < naux && aux[aux_count].num == i) {
        // Get the adjoint variables
        for (int k = 0; k < numAdjoints; k++) {
          memset(xptSens, 0, TACS_SPATIAL_DIM * len * sizeof(TacsScalar));

          // Get the element adjoint vector
          adjoint[k]->getValues(len, nodes, elemAdjoint);

          aux[aux_count].elem->addAdjResXptProduct(i, time, scale, elemAdjoint,
                                                   elemXpts, vars, dvars,
                                                   ddvars, xptSens);

          dfdXpt[k]->setValues(len, nodes, xptSens, TACS_ADD_VALUES);
        }
        aux_count++;
      }
    }
  }
}

/**
  Evaluate the derivative of an inner product of two vectors with a
  matrix of a given type. This code does not explicitly evaluate the
  element matrices. Instead, the inner product contribution from each
  element matrix is added to the final result. This implementation
  saves considerable computational time and memory.

  @param scale Scalar factor applied to the result
  @param matType The type of matrix
  @param psi The left-multiplying vector
  @param phi The right-multiplying vector
  @param dfdx The derivative vector
*/
void TACSAssembler::addMatDVSensInnerProduct(TacsScalar scale,
                                             ElementMatrixType matType,
                                             TACSBVec *psi, TACSBVec *phi,
                                             TACSBVec *dfdx) {
  psi->beginDistributeValues();
  if (phi != psi) {
    phi->beginDistributeValues();
  }
  psi->endDistributeValues();
  if (phi != psi) {
    phi->endDistributeValues();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemPsi, *elemPhi, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemPsi, &elemPhi, NULL, &elemXpts,
                  NULL, NULL, NULL);

  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *fdvSens = elementSensData;
  int *dvNums = elementSensIData;

  // Go through each element in the domain and compute the derivative
  // of the residuals with respect to each design variable and multiply by
  // the adjoint variables
  for (int i = 0; i < numElements; i++) {
    // Find the variables and nodes
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, elemVars);
    psi->getValues(len, nodes, elemPsi);
    phi->getValues(len, nodes, elemPhi);

    // Get the design variables for this element
    int numDVs = elements[i]->getDesignVarNums(i, maxDVs, dvNums);
    memset(fdvSens, 0, numDVs * designVarsPerNode * sizeof(TacsScalar));

    // Add the contribution to the design variable vector
    elements[i]->addMatDVSensInnerProduct(matType, i, time, scale, elemPsi,
                                          elemPhi, elemXpts, elemVars, numDVs,
                                          fdvSens);

    dfdx->setValues(numDVs, dvNums, fdvSens, TACS_ADD_VALUES);
  }
}

/**
  Evaluate the derivative of an inner product of two vectors with a
  matrix of a given type. This code does not explicitly evaluate the
  element matrices. Instead, the inner product contribution from each
  element matrix is added to the final result. This implementation
  saves considerable computational time and memory.

  @param scale Scalar factor applied to the result
  @param matType The type of matrix
  @param psi The left-multiplying vector
  @param phi The right-multiplying vector
  @param dfdXpt The derivative vector
*/
void TACSAssembler::addMatXptSensInnerProduct(TacsScalar scale,
                                              ElementMatrixType matType,
                                              TACSBVec *psi, TACSBVec *phi,
                                              TACSBVec *dfdXpt) {
  psi->beginDistributeValues();
  if (phi != psi) {
    phi->beginDistributeValues();
  }
  psi->endDistributeValues();
  if (phi != psi) {
    phi->endDistributeValues();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemPsi, *elemPhi, *elemXpts, *xptSens;
  getDataPointers(elementData, &elemVars, &elemPsi, &elemPhi, NULL, &elemXpts,
                  &xptSens, NULL, NULL);

  // Get the design variables from the elements on this process
  const int maxDVs = maxElementDesignVars;
  TacsScalar *fdvSens = elementSensData;
  int *dvNums = elementSensIData;

  // Go through each element in the domain and compute the derivative
  // of the residuals with respect to each design variable and multiply by
  // the adjoint variables
  for (int i = 0; i < numElements; i++) {
    // Find the variables and nodes
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, elemVars);
    psi->getValues(len, nodes, elemPsi);
    phi->getValues(len, nodes, elemPhi);

    memset(xptSens, 0, TACS_SPATIAL_DIM * len * sizeof(TacsScalar));

    // Add the contribution to the design variable vector
    elements[i]->addMatXptSensInnerProduct(
        matType, i, time, scale, elemPsi, elemPhi, elemXpts, elemVars, xptSens);

    dfdXpt->setValues(len, nodes, xptSens, TACS_ADD_VALUES);
  }
}

/**
  Evaluate the derivative of the inner product of two vectors with a
  matrix with respect to the state variables. This is only defined for
  nonlinear matrices, like the geometric stiffness matrix.  Instead of
  computing the derivative of the matrix for each vector component and
  then computing the inner product, this code computes the derivative
  of the inner product directly, saving computational time and memory.

  @param matType The type of matrix
  @param psi The left-multiplying vector
  @param phi The right-multiplying vector
  @param dfdu The derivative of the inner product w.r.t. the state vars
*/
void TACSAssembler::evalMatSVSensInnerProduct(ElementMatrixType matType,
                                              TACSBVec *psi, TACSBVec *phi,
                                              TACSBVec *dfdu) {
  // Zero the entries in the residual vector
  dfdu->zeroEntries();

  // Distribute the variable values
  psi->beginDistributeValues();
  if (phi != psi) {
    phi->beginDistributeValues();
  }
  psi->endDistributeValues();
  if (phi != psi) {
    phi->endDistributeValues();
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemPsi, *elemPhi, *elemRes, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemPsi, &elemPhi, &elemRes,
                  &elemXpts, NULL, NULL, NULL);

  // Go through each element in the domain and compute the derivative
  // of the residuals with respect to each design variable and multiply by
  // the adjoint variables
  for (int i = 0; i < numElements; i++) {
    // Find the variables and nodes
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, elemVars);
    psi->getValues(len, nodes, elemPsi);
    phi->getValues(len, nodes, elemPhi);

    // Add the contribution to the design variable vector
    elements[i]->getMatSVSensInnerProduct(matType, i, time, elemPsi, elemPhi,
                                          elemXpts, elemVars, elemRes);

    // Add the residual values to the local residual array
    dfdu->setValues(len, nodes, elemRes, TACS_ADD_VALUES);
  }

  dfdu->beginSetValues(TACS_ADD_VALUES);
  dfdu->endSetValues(TACS_ADD_VALUES);

  // Apply the boundary conditions to the fully assembled vector
  dfdu->applyBCs(bcMap);
}

/**
  Evaluate a Jacobian-vector product of the input vector
  x and store the result in the output vector y.

  This code does not assemble a matrix, but does compute the
  element-wise matricies. This code is not a finite-difference
  matrix-vector product implementation.

  Since the element Jacobian matrices are computed exactly, we can
  evaluate either a regular matrix-product or the transpose matrix
  product.

  @param scale The scalar coefficient
  @param alpha Coefficient on the variables
  @param beta Coefficient on the time-derivative terms
  @param gamma Coefficient on the second time derivative term
  @param x The input vector
  @param y the output vector y <- y + scale*J^{Op}*x
  @param matOr The matrix orientation
*/
void TACSAssembler::addJacobianVecProduct(TacsScalar scale, TacsScalar alpha,
                                          TacsScalar beta, TacsScalar gamma,
                                          TACSBVec *x, TACSBVec *y,
                                          MatrixOrientation matOr) {
  x->beginDistributeValues();
  x->endDistributeValues();

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *yvars, *elemXpts;
  TacsScalar *elemWeights, *elemMat;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &yvars, &elemXpts, NULL,
                  &elemWeights, &elemMat);

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  // Loop over all the elements in the model
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    // Get the number of variables from the element
    int nvars = elements[i]->getNumVariables();

    // Compute and add the contributions to the Jacobian
    memset(elemMat, 0, nvars * nvars * sizeof(TacsScalar));
    elements[i]->addJacobian(i, time, alpha, beta, gamma, elemXpts, vars, dvars,
                             ddvars, yvars, elemMat);

    // Add the contribution to the residual and the Jacobian
    // from the auxiliary elements - if any
    while (aux_count < naux && aux[aux_count].num == i) {
      aux[aux_count].elem->addJacobian(i, time, alpha, beta, gamma, elemXpts,
                                       vars, dvars, ddvars, yvars, elemMat);
      aux_count++;
    }

    // Temporarily set the variable array as the element input array
    // and get the local variable input values from the local array.
    TacsScalar *xvars = vars;
    x->getValues(len, nodes, xvars);

    // Take the matrix vector product. Note the matrix is stored in
    // row-major order and BLAS assumes column-major order. As a
    // result, the transpose arguments are reversed.
    TacsScalar zero = 0.0;
    int incx = 1;
    if (matOr == TACS_MAT_NORMAL) {
      BLASgemv("T", &nvars, &nvars, &scale, elemMat, &nvars, xvars, &incx,
               &zero, yvars, &incx);
    } else {
      BLASgemv("N", &nvars, &nvars, &scale, elemMat, &nvars, xvars, &incx,
               &zero, yvars, &incx);
    }

    // Add the residual values
    y->setValues(len, nodes, yvars, TACS_ADD_VALUES);
  }

  // Add the dependent-variable residual from the dependent nodes
  y->beginSetValues(TACS_ADD_VALUES);
  y->endSetValues(TACS_ADD_VALUES);

  // Set the boundary conditions
  y->applyBCs(bcMap);
}

/*
  Compute the sizes of the element-wise data and temporary array needed
  for a matrix-free matrix-vector product
*/
void TACSAssembler::getMatrixFreeDataSize(ElementMatrixType matType,
                                          int *_data_size, int *_temp_size) {
  int data_size = 0;
  int temp_size = 0;

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  for (int i = 0; i < numElements; i++) {
    // Compute the data for the matrix-free vector product
    int dsize, tsize;
    elements[i]->getMatVecDataSizes(matType, i, &dsize, &tsize);
    data_size += dsize;
    if (tsize > temp_size) {
      temp_size = tsize;
    }

    // Add the contribution to the residual and the Jacobian
    // from the auxiliary elements - if any
    while (aux_count < naux && aux[aux_count].num == i) {
      aux[aux_count].elem->getMatVecDataSizes(matType, i, &dsize, &tsize);
      data_size += dsize;
      if (tsize > temp_size) {
        temp_size = tsize;
      }
      aux_count++;
    }
  }

  if (_data_size) {
    *_data_size = data_size;
  }
  if (_temp_size) {
    *_temp_size = temp_size;
  }
}

/**
  Compute the element-wise data for a matrix-free matrix-vector product.

  If the data array is NULL, the code computes the size of data array required.
  This is returned

  The memory required to store the data scales with the number of quadrature
  points in the mesh.

  @param matType The type of matrix to assemble
  @param alpha Coefficient for the variables
  @param beta Coefficient for the time-derivative terms
  @param gamma Coefficientfor the second time derivative term
  @param data The data array that is used (may be NULL)
*/
void TACSAssembler::assembleMatrixFreeData(ElementMatrixType matType,
                                           TacsScalar alpha, TacsScalar beta,
                                           TacsScalar gamma,
                                           TacsScalar data[]) {
  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *yvars, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &yvars, &elemXpts, NULL,
                  NULL, NULL);

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  // Loop over all the elements in the model
  for (int i = 0; i < numElements; i++) {
    // Extract the element node locations and variable values
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    if (data) {
      xptVec->getValues(len, nodes, elemXpts);
      varsVec->getValues(len, nodes, vars);
      dvarsVec->getValues(len, nodes, dvars);
      ddvarsVec->getValues(len, nodes, ddvars);
    }

    // Get the size of the data array
    int dsize, tsize;
    elements[i]->getMatVecDataSizes(matType, i, &dsize, &tsize);

    // Compute the data for the matrix-free vector product
    elements[i]->getMatVecProductData(matType, i, time, alpha, beta, gamma,
                                      elemXpts, vars, dvars, ddvars, data);
    data += dsize;

    // Add the contribution to the residual and the Jacobian
    // from the auxiliary elements - if any
    while (aux_count < naux && aux[aux_count].num == i) {
      aux[aux_count].elem->getMatVecDataSizes(matType, i, &dsize, &tsize);
      aux[aux_count].elem->getMatVecProductData(matType, i, time, alpha, beta,
                                                gamma, elemXpts, vars, dvars,
                                                ddvars, data);
      data += dsize;
      aux_count++;
    }
  }
}

/**
  Compute a matrix-free matrix-vector product

  Note that the matType parameter must be consistent with the parameter used
  when generating the original data.

  @param matType The type of matrix to assemble
  @param data The matrix-free data computed from assembleMatrixFreeData
  @param x The input vector
  @param y The vector containing the matrix-vector product
  @param matOr The orientation of the matrix
*/
void TACSAssembler::addMatrixFreeVecProduct(ElementMatrixType matType,
                                            const TacsScalar data[],
                                            TacsScalar temp[], TACSBVec *x,
                                            TACSBVec *y,
                                            MatrixOrientation matOr) {
  x->beginDistributeValues();
  x->endDistributeValues();

  // Retrieve pointers to temporary storage
  TacsScalar *xvars, *yvars;
  getDataPointers(elementData, &xvars, &yvars, NULL, NULL, NULL, NULL, NULL,
                  NULL);

  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (auxElements) {
    naux = auxElements->getAuxElements(&aux);
  }

  // Loop over all the elements in the model
  for (int i = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];

    // Get the number of variables from the element
    int nvars = elements[i]->getNumVariables();

    // Get the values
    x->getValues(len, nodes, xvars);

    // Compute and add the contributions to the Jacobian
    memset(yvars, 0, nvars * sizeof(TacsScalar));

    // Get the size of the data array
    int dsize, tsize;
    elements[i]->getMatVecDataSizes(matType, i, &dsize, &tsize);

    elements[i]->addMatVecProduct(matType, i, data, temp, xvars, yvars);
    data += dsize;

    // Add the contribution to the residual and the Jacobian
    // from the auxiliary elements - if any
    while (aux_count < naux && aux[aux_count].num == i) {
      aux[aux_count].elem->getMatVecDataSizes(matType, i, &dsize, &tsize);
      aux[aux_count].elem->addMatVecProduct(matType, i, data, temp, xvars,
                                            yvars);
      data += dsize;
      aux_count++;
    }

    // Add the residual values
    y->setValues(len, nodes, yvars, TACS_ADD_VALUES);
  }

  // Add the dependent-variable residual from the dependent nodes
  y->beginSetValues(TACS_ADD_VALUES);
  y->endSetValues(TACS_ADD_VALUES);

  // Set the boundary conditions
  y->applyBCs(bcMap);
}

/**
  Test the implementation of the given element number.

  This tests the stiffness matrix and various parts of the
  design-sensitivities: the derivative of the determinant of the
  Jacobian, the derivative of the strain w.r.t. the nodal coordinates,
  and the state variables and the derivative of the residual w.r.t.
  the design variables and nodal coordiantes.

  @param elemNum The local element index to test
  @param print_level Print level to use
  @param dh Finite-difference (or complex-step) step length
  @param rtol Relative tolerance to apply
  @param atol Absolute tolerance to apply
*/
void TACSAssembler::testElement(int elemNum, int print_level, double dh,
                                double rtol, double atol) {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call testElement() before initialize()\n",
            mpiRank);
    return;
  } else if (elemNum < 0 || elemNum >= numElements) {
    fprintf(stderr, "[%d] Element number %d out of range [0,%d)\n", mpiRank,
            elemNum, numElements);
    return;
  }

  // Retrieve pointers to temporary storage
  TacsScalar *res;
  TacsScalar *vars, *dvars, *ddvars, *elemXpts;
  getDataPointers(elementData, &vars, &dvars, &ddvars, &res, &elemXpts, NULL,
                  NULL, NULL);

  int ptr = elementNodeIndex[elemNum];
  int len = elementNodeIndex[elemNum + 1] - ptr;
  const int *nodes = &elementTacsNodes[ptr];
  xptVec->getValues(len, nodes, elemXpts);
  varsVec->getValues(len, nodes, vars);
  dvarsVec->getValues(len, nodes, dvars);
  ddvarsVec->getValues(len, nodes, ddvars);

  int col = -1;
  TacsTestElementJacobian(elements[elemNum], elemNum, time, elemXpts, vars,
                          dvars, ddvars, col, dh, print_level, rtol, atol);

  const int maxDVs = maxElementDesignVars;
  TacsScalar *x = elementSensData;
  elements[elemNum]->getDesignVars(elemNum, maxDVs, x);
  TacsTestAdjResProduct(elements[elemNum], elemNum, time, elemXpts, vars, dvars,
                        ddvars, maxDVs, x, dh, print_level, rtol, atol);
  TacsTestAdjResXptProduct(elements[elemNum], elemNum, time, elemXpts, vars,
                           dvars, ddvars, dh, print_level, rtol, atol);

  // Test the residual computation with the Lagrange multipliers zeroed.
  // This does not work for most elements. The test requires that
  // the kinetic and potential energy functions be properly implemented. Results
  // from this test should be used with caution.
  // if (elements[elemNum]->getVarsPerNode() == 8){
  //   for ( int i = 7; i < elements[elemNum]->getNumVariables(); i++ ){
  //     vars[i] = 0.0;
  //   }
  // }
  // TacsTestElementResidual(elements[elemNum], elemNum, time, elemXpts,
  //                         vars, dvars, ddvars,
  //                         dh, print_level, rtol, atol);
  TACSElementModel *model = elements[elemNum]->getElementModel();
  if (model) {
    TacsTestElementModel(model, elemNum, time, dh, print_level, rtol, atol);
  }
  TacsTestElementMatDVSens(elements[elemNum], TACS_MASS_MATRIX, elemNum, time,
                           elemXpts, vars, maxDVs, x, dh, print_level, rtol,
                           atol);
  TacsTestElementMatDVSens(elements[elemNum], TACS_STIFFNESS_MATRIX, elemNum,
                           time, elemXpts, vars, maxDVs, x, dh, print_level,
                           rtol, atol);
  TacsTestElementMatDVSens(elements[elemNum], TACS_GEOMETRIC_STIFFNESS_MATRIX,
                           elemNum, time, elemXpts, vars, maxDVs, x, dh,
                           print_level, rtol, atol);
  TacsTestElementMatSVSens(elements[elemNum], TACS_GEOMETRIC_STIFFNESS_MATRIX,
                           elemNum, time, elemXpts, vars, dh, print_level, rtol,
                           atol);
}

/**
  Test the implementation of the function.

  This tests the state variable sensitivities and the design variable
  sensitivities of the function of interest. These sensitivities are
  computed based on a random perturbation of the input values.  Note
  that a system of equations should be solved - or the variables
  should be set randomly before calling this function, otherwise this
  function may produce unrealistic function values.

  Note that this function uses a central difference if the real code
  is compiled, and a complex step approximation if the complex version
  of the code is used.

  @param func The function to test
  @param dh The finite-difference or complex-step step size
*/
void TACSAssembler::testFunction(TACSFunction *func, double dh) {
  if (!meshInitializedFlag) {
    fprintf(stderr, "[%d] Cannot call testFunction() before initialize()\n",
            mpiRank);
    return;
  }

  // First, test the design variable values
  TACSBVec *x = createDesignVec();
  TACSBVec *xpert = createDesignVec();
  TACSBVec *xtemp = createDesignVec();
  x->incref();
  xpert->incref();
  xtemp->incref();

  xpert->setRand(-1.0, 1.0);
  getDesignVars(x);
  setDesignVars(x);

  TacsScalar fd = 0.0;
#ifdef TACS_USE_COMPLEX
  // Compute the function at the point x + dh*xpert
  xtemp->copyValues(x);
  xtemp->axpy(TacsScalar(0.0, dh), xpert);
  setDesignVars(xtemp);
  evalFunctions(1, &func, &fd);
  fd = TacsImagPart(fd) / dh;
#else
  // Compute the function at the point x + dh*xpert
  xtemp->copyValues(x);
  xtemp->axpy(dh, xpert);
  setDesignVars(xtemp);
  TacsScalar fval0;
  evalFunctions(1, &func, &fval0);

  // Compute the function at the point x - dh*xpert
  xtemp->copyValues(x);
  xtemp->axpy(-dh, xpert);
  setDesignVars(xtemp);
  TacsScalar fval1;
  evalFunctions(1, &func, &fval1);
  fd = 0.5 * (fval0 - fval1) / dh;
#endif  // TACS_USE_COMPLEX

  // Compute df/dx
  TacsScalar coef = 1.0;
  TacsScalar ftmp;
  setDesignVars(x);
  evalFunctions(1, &func, &ftmp);
  xtemp->zeroEntries();
  addDVSens(coef, 1, &func, &xtemp);

  // Compute df/dx^{T} * xpert
  TacsScalar pdf = xtemp->dot(xpert);

  if (mpiRank == 0) {
    fprintf(stderr, "Testing function %s\n", func->getObjectName());
    const char *descript = "df/dx^{T}p";
    fprintf(stderr, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
            "Analytic", "Approximate", "Rel. Error");
    if (pdf != 0.0) {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, 0,
              TacsRealPart(pdf), TacsRealPart(fd),
              fabs(TacsRealPart((pdf - fd) / pdf)));
    } else {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e\n", descript, 0, TacsRealPart(pdf),
              TacsRealPart(fd));
    }
  }

  x->decref();
  xtemp->decref();
  xpert->decref();

  TACSBVec *temp = createVec();
  TACSBVec *pert = createVec();
  TACSBVec *vars = createVec();
  temp->incref();
  pert->incref();
  vars->incref();

  getVariables(vars, NULL, NULL);

  // Set up a random perturbation
  pert->setRand(-1.0, 1.0);
  pert->applyBCs(bcMap);
  pert->scale(vars->norm() / pert->norm());

#ifdef TACS_USE_COMPLEX
  // Evaluate the function at vars + dh*pert
  temp->copyValues(vars);
  temp->axpy(TacsScalar(0.0, dh), pert);
  setVariables(temp);

  evalFunctions(1, &func, &fd);
  fd = TacsImagPart(fd) / dh;
#else
  // Evaluate the function at vars + dh*pert
  temp->copyValues(vars);
  temp->axpy(dh, pert);
  setVariables(temp);
  evalFunctions(1, &func, &fval0);

  // Evaluate the function at vars - dh*pert
  temp->copyValues(vars);
  temp->axpy(-dh, pert);
  setVariables(temp);
  evalFunctions(1, &func, &fval1);

  fd = 0.5 * (fval0 - fval1) / dh;
#endif  // TACS_USE_COMPLEX

  // Reset the variable values
  setVariables(vars);

  // Evaluate the state variable sensitivity
  evalFunctions(1, &func, &ftmp);
  TacsScalar alpha = 1.0, beta = 0.0, gamma = 0.0;
  temp->zeroEntries();
  addSVSens(alpha, beta, gamma, 1, &func, &temp);
  pdf = temp->dot(pert);

  if (mpiRank == 0) {
    const char *descript = "df/du^{T}p";
    fprintf(stderr, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
            "Analytic", "Approximate", "Rel. Error");
    if (pdf != 0.0) {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, 0,
              TacsRealPart(pdf), TacsRealPart(fd),
              fabs(TacsRealPart((pdf - fd) / pdf)));
    } else {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e\n", descript, 0, TacsRealPart(pdf),
              TacsRealPart(fd));
    }
  }

  temp->decref();
  pert->decref();
  vars->decref();

  TACSBVec *Xtemp = createNodeVec();
  TACSBVec *Xpert = createNodeVec();
  TACSBVec *Xvars = createNodeVec();
  Xtemp->incref();
  Xpert->incref();
  Xvars->incref();

  // Get the node locations
  getNodes(Xvars);

  // Set up a random perturbation
  Xpert->setRand(-1.0, 1.0);

#ifdef TACS_USE_COMPLEX
  // Evaluate the function at vars + dh*pert
  Xtemp->copyValues(Xvars);
  Xtemp->axpy(TacsScalar(0.0, dh), Xpert);
  setNodes(Xtemp);

  evalFunctions(1, &func, &fd);
  fd = TacsImagPart(fd) / dh;
#else
  // Evaluate the function at vars + dh*pert
  Xtemp->copyValues(Xvars);
  Xtemp->axpy(dh, Xpert);
  setNodes(Xtemp);
  evalFunctions(1, &func, &fval0);

  // Evaluate the function at vars - dh*pert
  Xtemp->copyValues(Xvars);
  Xtemp->axpy(-dh, Xpert);
  setNodes(Xtemp);
  evalFunctions(1, &func, &fval1);

  fd = 0.5 * (fval0 - fval1) / dh;
#endif  // TACS_USE_COMPLEX

  // Reset the variable values
  setNodes(Xvars);

  // Evaluate the state variable sensitivity
  evalFunctions(1, &func, &ftmp);
  Xtemp->zeroEntries();
  addXptSens(1.0, 1, &func, &Xtemp);
  Xtemp->beginSetValues(TACS_ADD_VALUES);
  Xtemp->endSetValues(TACS_ADD_VALUES);
  pdf = Xtemp->dot(Xpert);

  if (mpiRank == 0) {
    const char *descript = "df/dX^{T}p";
    fprintf(stderr, "%*s[   ] %15s %15s %15s\n", (int)strlen(descript), "Val",
            "Analytic", "Approximate", "Rel. Error");
    if (pdf != 0.0) {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e %15.4e\n", descript, 0,
              TacsRealPart(pdf), TacsRealPart(fd),
              fabs(TacsRealPart((pdf - fd) / pdf)));
    } else {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e\n", descript, 0, TacsRealPart(pdf),
              TacsRealPart(fd));
    }
  }

  Xtemp->decref();
  Xpert->decref();
  Xvars->decref();
}

/**
  Determine the number of components defined by elements in the
  TACSAssembler object.

  This call is collective - the number of components is obtained
  by a global reduction.

  @return The number of components
*/
int TACSAssembler::getNumComponents() {
  // Find the maximum component number on this processor
  int max_comp_num = 0;
  for (int i = 0; i < numElements; i++) {
    if (elements[i]->getComponentNum() >= max_comp_num) {
      max_comp_num = elements[i]->getComponentNum();
    }
  }
  max_comp_num++;

  int num_components = 1;
  MPI_Allreduce(&max_comp_num, &num_components, 1, MPI_INT, MPI_MAX, tacs_comm);
  return num_components;
}

/**
  Go through each element and get the output data for that element.

  The data is stored point-wise with each variable stored contiguously
  for each new point within the connectivity list. This stores the
  data at a point in memory indexed by data[node*nvals]. However,
  fewer than 'nvals' entries may be written in this call. The
  remaining data may be design variable entries that are computed
  below.

  @param elem_type The element type to match
  @param write_flag Binary flag indicating the components to write
  @param len The number of points (output)
  @param nvals The number of values at each point (output)
  @param data The data for each point for each value (output)
*/
void TACSAssembler::getElementOutputData(ElementType elem_type, int write_flag,
                                         int *_len, int *_nvals,
                                         TacsScalar **_data) {
  int nvals = TacsGetTotalOutputCount(elem_type, write_flag);
  int len = 0;
  for (int i = 0; i < numElements; i++) {
    len += elements[i]->getNumNodes();
  }
  TacsScalar *data = new TacsScalar[len * nvals];
  memset(data, 0, len * nvals * sizeof(TacsScalar));

  // Retrieve pointers to temporary storage
  TacsScalar *elemXpts, *vars, *dvars, *ddvars;
  getDataPointers(elementData, &vars, &dvars, &ddvars, NULL, &elemXpts, NULL,
                  NULL, NULL);

  for (int i = 0, offset = 0; i < numElements; i++) {
    int ptr = elementNodeIndex[i];
    int len = elementNodeIndex[i + 1] - ptr;
    const int *nodes = &elementTacsNodes[ptr];
    xptVec->getValues(len, nodes, elemXpts);
    varsVec->getValues(len, nodes, vars);
    dvarsVec->getValues(len, nodes, dvars);
    ddvarsVec->getValues(len, nodes, ddvars);

    elements[i]->getOutputData(i, elem_type, write_flag, elemXpts, vars, dvars,
                               ddvars, nvals, &data[offset]);

    offset += nvals * elements[i]->getNumNodes();
  }

  // Set the output pointers
  *_nvals = nvals;
  *_len = len;
  *_data = data;
}
