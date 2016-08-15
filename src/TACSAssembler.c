#include "TACSAssembler.h"

// Reordering implementation 
#include "FElibrary.h"
#include "MatUtils.h"
#include "AMDInterface.h"
#include "amd.h"
#include "tacsmetis.h"

// BLAS/LAPACK header
#include "tacslapack.h"

/*
  TACSAssembler implementation

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
  Not for commercial purposes.
*/

/*
  Constructor for the TACSAssembler object

  input:
  tacs_comm:           the TACS communicator 
  numOwnedNodes:       the number of locally-owned nodes
  varsPerNode:         the number of degrees of freedom per node
  numElements:         the number of elements in the mesh
  numNodes:            the number of nodes in the mesh
  nodeMaxCSRSize:      the maximum CSR size - expanded automatically
*/
TACSAssembler::TACSAssembler( MPI_Comm _tacs_comm, 
			      int numOwnedNodes, int _varsPerNode,
			      int _numElements, int _numNodes,
			      int _nodeMaxCSRsize ){
  init(_tacs_comm, numOwnedNodes, _varsPerNode,
       _numElements, _numNodes, 0,
       _nodeMaxCSRsize);
}

/*
  Constructor for the TACSAssembler object

  input:
  tacs_comm:           the TACS communicator 
  numOwnedNodes:       the number of locally-owned nodes
  varsPerNode:         the number of degrees of freedom per node
  numElements:         the number of elements in the mesh
  numNodes:            the number of nodes in the mesh
  nodeMaxCSRSize:      the maximum CSR size - expanded automatically
*/
TACSAssembler::TACSAssembler( MPI_Comm _tacs_comm, 
			      int numOwnedNodes, int _varsPerNode,
			      int _numElements, int _numNodes,
			      int _numDependentNodes,
			      int _nodeMaxCSRsize ){
  init(_tacs_comm, numOwnedNodes, _varsPerNode,
       _numElements, _numNodes, _numDependentNodes,
       _nodeMaxCSRsize);
}

/*
  Initialize the TACSAssembler object

  tacs_comm:           the TACS communicator 
  numOwnedNodes:       the number of locally-owned nodes
  varsPerNode:         the number of degrees of freedom per node
  numElements:         the number of elements in the mesh
  numNodes:            the number of nodes in the mesh
  nodeMaxCSRSize:      the maximum CSR size - expanded automatically
*/
void TACSAssembler::init( MPI_Comm _tacs_comm, 
			  int numOwnedNodes, int _varsPerNode,
			  int _numElements, int _numNodes,
			  int _numDependentNodes,
			  int _nodeMaxCSRsize ){
  TacsInitialize();
  
  // Copy the communicator for MPI communication
  tacs_comm = _tacs_comm;

  // If MPI is being used, get the rank of the current 
  // process and the total number of processes
  MPI_Comm_rank(tacs_comm, &mpiRank);
  MPI_Comm_size(tacs_comm, &mpiSize);

  // Set the simulation time to 0
  time = 0.0;

  // Now set up the default pthread info
  thread_info = new TACSThreadInfo(1); // Set up the info class with 1 thread
  thread_info->incref();
  pthread_mutex_init(&tacs_mutex, NULL);
  tacsPInfo = new TACSAssemblerPthreadInfo();
  numCompletedElements = 0;

  // copy data to be used later in the program
  varsPerNode = _varsPerNode;
  numElements = _numElements;
  numNodes = _numNodes;
  numDependentNodes = _numDependentNodes;
  nodeMaxCSRsize = _nodeMaxCSRsize;
  nodeCSRIncrement = (int)(0.5*nodeMaxCSRsize);

  // Print out the number of local nodes and elements
  printf("[%d] Creating TACSAssembler with numOwnedNodes = %d \
numElements = %d\n", mpiRank, numOwnedNodes, numElements);

  // Calculate the total number of nodes and elements
  int info[2], recv_info[2];
  info[0] = numOwnedNodes;
  info[1] = numElements;
  MPI_Reduce(info, recv_info, 2, MPI_INT, MPI_SUM, 0, tacs_comm);
  if (mpiSize > 1 && mpiRank == 0){
    printf("[%d] TACSAssembler: Total dof = %d Total nodes = %d \
Total elements = %d\n", mpiRank, varsPerNode*recv_info[0], 
           recv_info[0], recv_info[1]);
  }

  // Initialize some information about the number of items
  meshFinalizedFlag = 0;
  maxElementNodes = 0;
  maxElementSize = 0;
  maxElementIndepNodes = 0;

  // Set the auxiliary element class to NULL
  aux_elements = NULL;

  // Information for setting boundary conditions and distributing variables
  varMap = new VarMap(tacs_comm, numOwnedNodes);
  varMap->incref();

  // Estimate 100 bcs at first, but this is expanded as required
  int nbc_est = 100; 
  bcMap = new BCMap(nbc_est);
  bcMap->incref();

  // Set the distribution object to NULL at first
  vecDist = NULL;
  vecDistIndices = NULL;

  // FEMat-specific objects
  feMatBIndices = feMatCIndices = NULL;
  feMatBMap = feMatCMap = NULL;

  // Keep track of the current total nodes and elements
  currNode = 0;
  currElement = 0;
  
  // Allocate memory for the tacs node numbers
  tacsNodeNums = new int[ numNodes ];
  for ( int i = 0; i < numNodes; i++ ){
    tacsNodeNums[i] = -1;
  }

  // Allocate space for the global nodes
  int node_size = TACS_SPATIAL_DIM*(numNodes + numDependentNodes);
  Xpts = new TacsScalar[ node_size ];
  memset(Xpts, 0, node_size*sizeof(TacsScalar));

  // Allocate element-> node information
  elementNodeIndex = new int[ numElements+1 ]; 
  elementLocalNodes = new int[ nodeMaxCSRsize ];
  elementTacsNodes  = NULL;

  elementNodeIndex[0] = 0;

  // Allocate and initialize pointers to the elements
  elements = new TACSElement*[ numElements ];
  memset(elements, 0, numElements*sizeof(*elements));

  // Null out the dependent node data
  depNodeIndex = NULL;
  depNodeToLocal = NULL;
  depNodeToTacs = NULL;
  depNodeWeights = NULL;

  // Set the local residual array to NULL
  localRes = NULL;
    
  // Set the local velocity and acceleration vectors to NULL
  localVars = NULL;
  localDotVars = NULL;
  localDDotVars = NULL;

  // Set the local element data to NULL
  elementData = NULL;
  elementIData = NULL;
}

/*
  The TACSAssembler destructor. 

  Clean up the allocated memory and decref() all objects
*/
TACSAssembler::~TACSAssembler(){
  TacsFinalize();

  pthread_mutex_destroy(&tacs_mutex);
  delete tacsPInfo;

  // Go through and decref all the elements
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]){
      elements[i]->decref();
    }
  }
  delete [] elements;

  // Delete nodal information
  delete [] elementNodeIndex;
  delete [] elementLocalNodes;
  if (elementTacsNodes){ delete [] elementTacsNodes; }

  if (depNodeIndex){ delete [] depNodeIndex; }
  if (depNodeToLocal){ delete [] depNodeToLocal; }
  if (depNodeToTacs){ delete [] depNodeToTacs; }
  if (depNodeWeights){ delete [] depNodeWeights; }

  delete [] Xpts;

  // Decrease the reference count
  if (aux_elements){ aux_elements->decref(); }

  // Decrease the reference count to objects allocated in finalize
  if (varMap){ varMap->decref(); }
  if (bcMap){ bcMap->decref(); }
  if (vecDist){ vecDist->decref(); }

  if (feMatBIndices){ feMatBIndices->decref(); }
  if (feMatCIndices){ feMatCIndices->decref(); }
  if (feMatBMap){ feMatBMap->decref(); }
  if (feMatCMap){ feMatCMap->decref(); }

  if (vecDistIndices){ vecDistIndices->decref(); }
  else { delete [] tacsNodeNums; }  

  // Delete arrays allocated in initializeArrays()
  if (localRes){ delete [] localRes; }
  if (localVars){ delete [] localVars; }
  if (localDotVars){ delete [] localDotVars; }
  if (localDDotVars){ delete [] localDDotVars; }
  if (elementData){ delete [] elementData; }
  if (elementIData){ delete [] elementIData; }

  // Decref the thread information class
  thread_info->decref();
}

const char * TACSAssembler::tacsName = "TACSAssembler";

MPI_Comm TACSAssembler::getMPIComm(){
  return tacs_comm;
}

/*!
  Add a node to TACS. 

  Each node is associated with a local and a global numbering scheme.
*/
void TACSAssembler::addNode( int _localNodeNum, int _tacsNodeNum ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addNode() after finalize()\n", 
	    mpiRank);
    return;
  }

  if (_localNodeNum < 0 || _localNodeNum >= numNodes){
    fprintf(stderr, "[%d] TACSAssembler::addNode(): Local node number %d \
out of range [0,%d)\n", mpiRank, _localNodeNum, numNodes);
  }
  else if (_tacsNodeNum < 0){
    fprintf(stderr, "[%d] TACSAssembler::addNode(): Negative global \
node index %d\n", mpiRank, _tacsNodeNum);
  }
  else {
    tacsNodeNums[ _localNodeNum ] = _tacsNodeNum;    
    currNode++;
  }
}

/*!
  Add a list of nodes to TACS.
*/
void TACSAssembler::addNodes( int _localNodeNums[], int _tacsNodeNums[], 
			      int _numNodes ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addNodes() after finalize()\n", 
	    mpiRank);
    return;
  }  

  for ( int i = 0; i < _numNodes; i++ ){
    if (_localNodeNums[i] < 0 || _localNodeNums[i] >= numNodes){
      fprintf(stderr, "[%d] TACSAssembler::addNodes(): Local node number %d \
out of range [0,%d)\n", mpiRank, _localNodeNums[i], numNodes);
    }
    else if (_tacsNodeNums[i] < 0){
      fprintf(stderr, "[%d] TACSAssembler::addNode(): Negative global \
node index %d\n", mpiRank, _tacsNodeNums[i]);
    }
    else {
      tacsNodeNums[ _localNodeNums[i] ] = _tacsNodeNums[i];          
      currNode++;
    }
  }
}

/*!
  Add all the tacsNodeNums at once.
*/
void TACSAssembler::addNodes( int ** _tacsNodeNums ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addNodes() after finalize()\n", 
	    mpiRank);
    return;
  }

  if (tacsNodeNums){ delete [] tacsNodeNums; }
  tacsNodeNums = *_tacsNodeNums;
  currNode = numNodes;
  
  *_tacsNodeNums = NULL;
}

/*
  Set the dependent node data structure
*/
void TACSAssembler::setDependentNodes( int **_depNodeIndex, 
				       int **_depNodeToLocal,
				       double **_depNodeWeights ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call setDependentNodes() after finalize()\n", 
	    mpiRank);
    return;
  }

  if (depNodeIndex){ delete [] depNodeIndex; }
  if (depNodeToLocal){ delete [] depNodeToLocal; }
  if (depNodeWeights){ delete [] depNodeWeights; }

  if (_depNodeIndex){
    depNodeIndex = *_depNodeIndex;
    *_depNodeIndex = NULL;
  }
  if (_depNodeToLocal){
    depNodeToLocal = *_depNodeToLocal;
    *_depNodeToLocal = NULL;
  }
  if (_depNodeWeights){
    depNodeWeights = *_depNodeWeights;
    *_depNodeWeights = NULL;
  }
}

/*!
  Set the boundary conditions object that will be associated with the
  vectors/matrices that are created using TACSAssembler.
*/
void TACSAssembler::addBC( int nodeNum, const int bcNums[], int nbcs ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addBC() before finalize()\n", mpiRank);
  }
  
  if (nodeNum < 0 || nodeNum >= numNodes){
    fprintf(stderr, "[%d] Warning, addBC node number %d out of range\n", 
	    mpiRank, nodeNum);
  }
  else {
    bcMap->addBC(nodeNum, tacsNodeNums[nodeNum], bcNums, NULL, nbcs);
  }
}

/*!
  Set the boundary conditions into the object that will be associated
  with the vectors/matrices and with the TACSAssembler object.
*/
void TACSAssembler::addBC( int nodeNum, const int bcNums[], 
			   const TacsScalar bcVals[], int nbcs ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addBC() before finalize()\n", mpiRank);
  }
  
  if (nodeNum < 0 || nodeNum >= numNodes){
    fprintf(stderr, "[%d] Warning, addBC node number %d out of range\n", 
	    mpiRank, nodeNum);
  }
  else {   
    bcMap->addBC(nodeNum, tacsNodeNums[nodeNum], bcNums, bcVals, nbcs);
  }
}

/*
  Apply the boundary conditions to the residual vector

  r = u - lambda * u_d

  where u_d are the prescribed displacements and lambda is
  the load factor
*/
void TACSAssembler::applyBCs( BVec *residual, 
			      const TacsScalar local_vars[] ){

  const int *local, *global, *var_ptr, *vars;
  const TacsScalar *values;
  int nbcs = bcMap->getBCs(&local, &global, &var_ptr, &vars, &values);

  TacsScalar * res;
  residual->getArray(&res);

  const int * ownerRange;
  int rank, size;
  varMap->getOwnerRange(&ownerRange, &rank, &size);

  for ( int i = 0; i < nbcs; i++ ){
    int bcvar = global[i];
    int local_var = local[i];

    if (bcvar >= ownerRange[rank] && 
        bcvar <  ownerRange[rank+1]){
      bcvar = bcvar - ownerRange[rank];
      TacsScalar * r = &res[varsPerNode*bcvar];
      const TacsScalar * u = &local_vars[varsPerNode*local_var];

      for ( int j = var_ptr[i], k = 0; j < var_ptr[i+1]; j++, k++ ){
	r[vars[j]] = u[vars[j]] - values[j];
      }
    }
  }
}

/*!
  Add an element to TACS with given node numbers.
*/
int TACSAssembler::addElement( TACSElement * element, int _localNodeNums[], 
			       int numElemNodes ){
  if (currElement+1 > numElements){
    fprintf(stderr, "[%d] TACSAssembler::addElement(): Error, cannot add \
any more elements\n", mpiRank);
    return -1;
  }
  else if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addElement() after finalize()\n", 
	    mpiRank);
    return -1;
  }  

  // Check if the number of variables per node matches
  if (element->numDisplacements() != varsPerNode){
    fprintf(stderr, "[%d] Variables per node does not match element\n",
            mpiRank);
    return -1;
  }

  // Determine if the maximum number of variables and nodes needs to be changed
  int elemSize = element->numVariables();
  if (elemSize  > maxElementSize){
    maxElementSize = elemSize;
  }

  // Add the local node numbers to the current element
  int start = elementNodeIndex[currElement];  

  // Check if adding this element to the connectivity 
  // will exceed the array bounds
  if (start + numElemNodes > nodeMaxCSRsize){
    // Increase the maximum size of the array
    nodeMaxCSRsize = nodeMaxCSRsize + (nodeCSRIncrement > numElemNodes ? 
				       nodeCSRIncrement : numElemNodes);
    matutils::ExtendArray(&elementLocalNodes, start, nodeMaxCSRsize);
  }

  for ( int j = 0; j < numElemNodes; j++ ){
    // Check that the _localNodeNums are in the correct range
    if (_localNodeNums[j] >= -numDependentNodes && 
	_localNodeNums[j] < numNodes){
      elementLocalNodes[start + j] = _localNodeNums[j];
    }
    else {
      elementLocalNodes[start + j] = 0;
      fprintf(stderr, "[%d] Warning element local node %d \
out of range. Setting to zero\n", mpiRank, _localNodeNums[j]);
    }
  }
   
  // Update the pointers to the next position in the CSR 
  // data structures for the next call
  elementNodeIndex[currElement + 1] = start + numElemNodes;
  element->incref(); // Increase the reference count to the element
  elements[currElement] = element;
  
  currElement++;
  return currElement-1;
}

/*!
  Add all the elements at once
*/
void TACSAssembler::addElements( TACSElement ***_elements ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call addElements() after finalize()\n", 
	    mpiRank);
    return;
  }

  if (elements){
    if (currElement > 0){
      for ( int i = 0; i < numElements; i++ ){
	if (elements[i]){
	  elements[i]->decref();
	}
      }
    }
    delete [] elements;
  }

  elements = *_elements;
  currElement = numElements;
  *_elements = NULL;

  for ( int i = 0; i < numElements; i++ ){
    // Determine if the maximum number of 
    // variables and nodes needs to be changed  
    int elemSize = elements[i]->numVariables();
    if (elemSize > maxElementSize){
      maxElementSize = elemSize;
    }
  }  
}

/*!
  Set the element connectivity array.

  TACSAssembler then takes control of the arrays and will delete them
  after they are initialized.  
*/
void TACSAssembler::setElementConnectivity( int **elemindex, 
					    int **elemnodes ){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call setElementConnectivity() \
after finalize()\n", mpiRank);
    return;
  }  
  
  if (elementLocalNodes){ delete [] elementLocalNodes; }
  if (elementNodeIndex){ delete [] elementNodeIndex; }

  elementLocalNodes = *elemnodes;  *elemnodes = NULL;
  elementNodeIndex = *elemindex;  *elemindex = NULL;
}

/*
  Create a global vector of node locations
*/
BVec *TACSAssembler::createNodeVec(){
  return new BVec(varMap, TACS_SPATIAL_DIM);
}

/*!
  Set the nodes from the node map object
*/
void TACSAssembler::setNodes( BVec *X ){
  if (X){
    vecDist->beginForward(X, Xpts);
    vecDist->endForward(X, Xpts);
  }

  // Set the locations of the dependent nodes
  setDependentVariables(TACS_SPATIAL_DIM, Xpts);
}

/*!
  Get the node locations from TACS
*/
void TACSAssembler::getNodes( BVec *X ){
  // Set the node locations
  vecDist->beginReverse(Xpts, X, BVecDistribute::INSERT);
  vecDist->endReverse(Xpts, X, BVecDistribute::INSERT);  
}

/*!
  Get the raw nodal locations
*/
void TACSAssembler::getNodeArray( TacsScalar **_Xpts ){
  *_Xpts = Xpts;
}

/*
  Set the auxiliary elements within the TACSAssembler object
  
  This only needs to be done once sometime during initialization.  If
  you need to change the loads repeatedly, this can be called
  repeatedly. No check is made at this point that you haven't done
  something odd. Note that the code assumes that the elements defined
  here perfectly overlap the non-zero pattern of the elements set
  internally within TACS already.
*/
void TACSAssembler::setAuxElements( TACSAuxElements *_aux_elems ){
  // Increase the reference count to the input elements (Note that
  // the input may be NULL to over-write the internal object
  if (_aux_elems){
    _aux_elems->incref();
  }

  // Decrease the reference count if the old object is not NULL
  if (aux_elements){
    aux_elements->decref();
  }
  aux_elements = _aux_elems;

  // Check whether the auxiliary elements match
  int naux = 0;
  TACSAuxElem *aux = NULL;
  naux = aux_elements->getAuxElements(&aux);
  for ( int k = 0; k < naux; k++ ){
    int elem = aux[k].num;
    if (elements[elem]->numVariables() != 
        aux[k].elem->numVariables()){
      fprintf(stderr, "[%d] Auxiliary element sizes do not match\n",
              mpiRank);
    }
  }
}

/*
  Retrieve the auxiliary element object from TACSAssembler

  Warning: The auxiliary element object may be NULL
*/
TACSAuxElements* TACSAssembler::getAuxElements(){
  return aux_elements;
}

/*!  
  Compute the local node numbers that correspond to the coupling
  nodes connected to elements on other processes.

  Sort the global node numbers. Match the intervals and send them off
  to the owning process. On the owner, scan through the arrays until
  all the local coupling nodes are found.  
*/
int TACSAssembler::computeCouplingNodes( int ** _cnodes ){
  int * nodes = new int[ numNodes ];
  memcpy(nodes, tacsNodeNums, numNodes*sizeof(int));

  int nnodes = FElibrary::uniqueSort(nodes, numNodes);
  if (nnodes != numNodes){
    fprintf(stderr, "[%d] TACSAssembler error, nodes must be \
uniquely defined! \n", mpiRank);
  }

  // Get the ownership range and match the intervals of ownership
  const int * ownerRange;
  int rank, size;
  varMap->getOwnerRange(&ownerRange, &rank, &size);

  int * ext_ptr   = new int[ mpiSize+1 ];
  int * ext_count = new int[ mpiSize ];

  FElibrary::matchIntervals(mpiSize, ownerRange, nnodes, nodes, ext_ptr);

  // Send the nodes owned by other processors the information
  // First count up how many will go to each process
  for ( int i = 0; i < mpiSize; i++ ){
    ext_count[i] = ext_ptr[i+1] - ext_ptr[i];
    if (i == mpiRank){ ext_count[i] = 0; }
  }

  int * recv_count = new int[ mpiSize ];
  int * recv_ptr   = new int[ mpiSize+1 ];

  MPI_Alltoall(ext_count, 1, MPI_INT, recv_count, 1, MPI_INT, tacs_comm);

  // Now, send the node numbers to the other processors
  recv_ptr[0] = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    recv_ptr[i+1] = recv_ptr[i] + recv_count[i];
  }

  // Number of nodes that will be received from other procs
  int * recv_nodes = new int[ recv_ptr[mpiSize] ];
  MPI_Alltoallv(nodes, ext_count, ext_ptr, MPI_INT, 
		recv_nodes, recv_count, recv_ptr, MPI_INT, tacs_comm);

  // Uniquely sort the recieved nodes
  int nextern_unique = FElibrary::uniqueSort(recv_nodes, recv_ptr[mpiSize]);

  int numOwnedNodes = ownerRange[mpiRank+1] - ownerRange[mpiRank];
  int ncnodes = nextern_unique + (numNodes - numOwnedNodes);
  int * cnodes = new int[ ncnodes ];
  memset(cnodes, 0, ncnodes*sizeof(int));

  // Find recv_nodes[i] in tacsNodeNums[j]
  // assign cnodes[i] = j
  int nc = 0;
  for ( int j = 0; j < numNodes; j++ ){
    int * item = (int*)bsearch(&tacsNodeNums[j], recv_nodes, nextern_unique, 
                               sizeof(int), FElibrary::comparator);
    if (item){
      cnodes[item-recv_nodes] = j;
      nc++;
      if (nc == nextern_unique){ break; }      
    }
  }

  // Now add in the nodes that are reference locally, but out of
  // the ownership range
  for ( int k = 0, j = nextern_unique; k < numNodes; k++ ){
    if (tacsNodeNums[k] < ownerRange[mpiRank] ||
        tacsNodeNums[k] >= ownerRange[mpiRank+1]){
      cnodes[j] = k;
      j++;
    }
  }

  ncnodes = FElibrary::uniqueSort(cnodes, ncnodes);

  delete [] nodes;
  delete [] ext_ptr;
  delete [] ext_count;
  delete [] recv_ptr;
  delete [] recv_count;
  delete [] recv_nodes;  

  *_cnodes = cnodes;
  return ncnodes;
}

/*!
  Compute the elements that couple with other processors.
  
  Compute the coupling nodes and the node to element pointer
  CSR data structure. From these, collect all elements that "own" a
  node that is referred to from another process.  
*/
int TACSAssembler::computeCouplingElements( int ** _celems ){
  int * cnodes;
  int ncnodes = computeCouplingNodes(&cnodes);

  int *nodeElemIndex, *nodeElems;
  computeNodeToElementCSR(&nodeElemIndex, &nodeElems);

  int ncelems = 0;
  int * celems = new int[ numElements ];

  for ( int i = 0; i < ncnodes; i++ ){
    int cnode = cnodes[i];
    for ( int j = nodeElemIndex[cnode]; j < nodeElemIndex[cnode+1]; j++ ){
      int elem = nodeElems[j];

      ncelems = FElibrary::mergeArrays(celems, ncelems, &elem, 1);      
    }
  }

  delete [] nodeElemIndex;
  delete [] nodeElems;

  *_celems = celems;
  return ncelems;
}

/*!
  Compute a reordering of the nodes.

  1. Determine the nodes that are locally owned and those that are
  required from other processors.
  2. Distribute the unknown nodes to the owning processes. 
  These are the recv_nodes -> sort them to sorted_recv_nodes.
  3. Find the indices of the recieved nodes.
  4. Order the locally owned nodes based on the input ordering.
  5. Based on the recvied indices, set the outgoing node numbers
  back into the recieving array.
  6. Set the new values of the nodes on the requesting processes

  The nodes are transfered as follows
  extern_nodes -> sorted_extern_nodes -> recv_nodes -> sorted_recv_nodes
  -> sorted_recv_index

  Node ordering
  
  recv_nodes[i] = new_nodes[sorted_recv_index[i]]
  recv_nodes -> new_sorted_extern_nodes
  find extern_nodes[i] == sorted_extern_nodes[j]
  assign new_var = new_sorted_extern_nodes[j]
*/
void TACSAssembler::computeReordering( enum OrderingType order_type,
                                       enum MatrixOrderingType mat_type ){
  // The ownership range of the nodes
  const int * ownerRange;
  int rank, size;
  varMap->getOwnerRange(&ownerRange, &rank, &size);

  // Figure out how many of the local nodes are not owned by this
  // processor.
  int nextern_nodes = numNodes - 
    (ownerRange[mpiRank+1] - ownerRange[mpiRank]);

  // Create a sorted list of externally owned node numbers
  int * sorted_extern_nodes = new int[ nextern_nodes ];

  // Go through all the node numbers and determine which ones are
  // owned by other processors
  nextern_nodes = 0;
  for ( int i = 0; i < numNodes; i++ ){
    if (tacsNodeNums[i] < ownerRange[mpiRank] ||
        tacsNodeNums[i] >= ownerRange[mpiRank+1]){
      sorted_extern_nodes[nextern_nodes] = tacsNodeNums[i];
      nextern_nodes++;
    }
  }

  // Create a unique list of non-local nodes. Check to see that the
  // unique list is equal in length to the original array - i.e. the
  // original nodes are uniquely defined.
  if (nextern_nodes != FElibrary::uniqueSort(sorted_extern_nodes, 
					     nextern_nodes)){
    fprintf(stderr, "[%d] Error, external nodes are not unique\n", 
	    mpiRank);
  }

  // Determine the range of non-local node numbers that are owned
  // by other processors. Match the intervals of these nodes within
  // the sorted external node list in preparation for sending their
  // values to the destination processors.
  int * extern_ptr = new int[ mpiSize+1 ];
  int * extern_count = new int[ mpiSize ];
  FElibrary::matchIntervals(mpiSize, ownerRange, 
			    nextern_nodes, sorted_extern_nodes, extern_ptr);

  // Determine the number of nodes destined for each processor
  for ( int i = 0; i < mpiSize; i++ ){
    extern_count[i] = extern_ptr[i+1] - extern_ptr[i];
  }

  // Set up the receive counts/offsets from the other processors
  int * recv_count = new int[ mpiSize ];
  int * recv_ptr = new int[ mpiSize+1 ];
  MPI_Alltoall(extern_count, 1, MPI_INT, recv_count, 1, MPI_INT, tacs_comm);

  // Count up the number of received nodes and offsets from all the
  // other processors to this processor
  recv_ptr[0] = 0;
  for ( int i = 0; i < mpiSize; i++ ){
    recv_ptr[i+1] = recv_ptr[i] + recv_count[i];
  }

  // Receive the nodes from the other processors using all-to-all
  // communication with the other processors
  int nrecv_nodes = recv_ptr[mpiSize];
  int * recv_nodes = new int[ nrecv_nodes ];
  MPI_Alltoallv(sorted_extern_nodes, extern_count, extern_ptr, MPI_INT, 
		recv_nodes, recv_count, recv_ptr, MPI_INT, tacs_comm);
   
  // Sort the list of nodes that were received by other processors. 
  // These are the nodes that this processor owns that are referenced
  // by other processors.
  int * sorted_recv_nodes = new int[ nrecv_nodes ];
  memcpy(sorted_recv_nodes, recv_nodes, nrecv_nodes*sizeof(int));
  int nrecv_sorted = FElibrary::uniqueSort(sorted_recv_nodes, nrecv_nodes);

  // Find the indices of the received nodes that are referenced
  // by other processors
  int * sorted_recv_index = new int[ nrecv_sorted ];  
  int nc = 0;
  for ( int j = 0; j < numNodes; j++ ){
    int * item = (int *) bsearch(&tacsNodeNums[j], sorted_recv_nodes, 
				 nrecv_sorted, sizeof(int), 
				 FElibrary::comparator);
    if (item){
      sorted_recv_index[item-sorted_recv_nodes] = j;
      nc++;
      if (nc == nrecv_sorted){ break; }
    }
  }

  if (nc != nrecv_sorted){
    fprintf(stderr, "[%d] Error, not all received nodes \
have been set\n", mpiRank);
  }

  // The new node numbers to be set according to the different
  // re-ordering schemes
  int * new_node_nums = new int[ numNodes ];

  // Metis can't handle the non-zero diagonal in the CSR data structure
  int no_diagonal = 0;
  if (order_type == ND_ORDER){ 
    no_diagonal = 1;
  }

  // If using only one processor, order everything. In this case
  // there is no distinction between local and global ordering.
  if (mpiSize == 1){
    // The node connectivity
    int *rowp, *cols;
    computeLocalNodeToNodeCSR(&rowp, &cols, no_diagonal);
    computeMatReordering(order_type, numNodes, rowp, cols, 
                         NULL, new_node_nums);

    delete [] rowp;
    delete [] cols;
  }
  else {
    // First, find the reduced nodes - the set of nodes 
    // that are only referenced by this processor. These 
    // can be reordered without affecting other processors.
    int * reduced_nodes = new int[ numNodes ];
    memset(reduced_nodes, 0, numNodes*sizeof(int));

    // Label all the external nodes with a negative index.
    // This will eliminate them from the connectivity.
    for ( int i = 0; i < numNodes; i++ ){
      if (tacsNodeNums[i] < ownerRange[mpiRank] ||
          tacsNodeNums[i] >= ownerRange[mpiRank+1]){
        reduced_nodes[i] = -1;
      }
    }

    if (mat_type == APPROXIMATE_SCHUR ||
        mat_type == DIRECT_SCHUR){
      for ( int i = 0; i < nrecv_sorted; i++ ){
        int cnode = sorted_recv_index[i];
        reduced_nodes[cnode] = -1;
      }   

      // If we want an approximate schur ordering, where the
      // nodes that couple to other processors are ordered last,
      // we also add these nodes to the reduced set.
      if (mat_type == APPROXIMATE_SCHUR){
        // Compute the ordering for the
        int *rowp, *cols;
        computeLocalNodeToNodeCSR(&rowp, &cols);

        // Order all nodes linked by an equation to an external node
        // This ordering is better for the PMat class
        for ( int i = 0; i < nrecv_sorted; i++ ){
          int var = sorted_recv_index[i];
          for ( int j = rowp[var]; j < rowp[var+1]; j++ ){
            reduced_nodes[cols[j]] = -1;
          }
        }

        // Add to the reduced set of nodes all those nodes that
        // couple to a boundary node through the governing
        // equations. This ensures that the equations
        // that reference a boundary node are ordered last as well.
        for ( int i = 0; i < numNodes; i++ ){
          if (tacsNodeNums[i] < ownerRange[mpiRank] ||
              tacsNodeNums[i] >= ownerRange[mpiRank+1]){
            for ( int j = rowp[i]; j < rowp[i+1]; j++ ){
              reduced_nodes[cols[j]] = -1;
            }
          }
        }

        delete [] rowp;
        delete [] cols;
      }
    }   

    // Now, order all nodes that are non-negative
    int num_reduced_nodes = 0;
    for ( int i = 0; i < numNodes; i++ ){
      if (reduced_nodes[i] >= 0){
        reduced_nodes[i] = num_reduced_nodes;
        num_reduced_nodes++;
      }
    }

    // Compute the reordering for the reduced set of nodes
    int * new_reduced_vars = new int[ num_reduced_nodes ];
    int *rowp, *cols;
    computeLocalNodeToNodeCSR(&rowp, &cols, 
                              num_reduced_nodes, reduced_nodes, 
                              no_diagonal);
    computeMatReordering(order_type, num_reduced_nodes, rowp, cols,
                         NULL, new_reduced_vars);

    // Place the result back into the new_node_nums - add the
    // ownership offset
    int node_offset = ownerRange[mpiRank];
    for ( int i = 0, j = 0; i < numNodes; i++ ){
      if (reduced_nodes[i] >= 0){
        new_node_nums[i] = node_offset + new_reduced_vars[j];
        j++;
      }
    }

    // Add the offset to the total number of reduced nodes
    node_offset += num_reduced_nodes; 

    delete [] rowp;
    delete [] cols;
    delete [] new_reduced_vars;

    // Now, order any remaining variables
    num_reduced_nodes = 0;
    for ( int i = 0; i < numNodes; i++ ){
      // If the node has not been ordered and is within the ownership
      // range of this process, order it now.
      if (reduced_nodes[i] < 0 && 
          (tacsNodeNums[i] >= ownerRange[mpiRank] &&
           tacsNodeNums[i] <  ownerRange[mpiRank+1])){
        reduced_nodes[i] = num_reduced_nodes; 
        num_reduced_nodes++;
      }
      else {
        reduced_nodes[i] = -1;
      }
    }

    if (num_reduced_nodes > 0){
      // Additive Schwarz ordering should number all locally owned
      // nodes first, and should not require this second ordering.
      if (mat_type == ADDITIVE_SCHWARZ){
        fprintf(stderr, "[%d] Error in additive Schwarz reordering\n",
                mpiRank);
      }

      // Order any remaning variables that are locally owned
      new_reduced_vars = new int[ num_reduced_nodes ];
      computeLocalNodeToNodeCSR(&rowp, &cols, 
                                num_reduced_nodes, reduced_nodes, 
                                no_diagonal);
      computeMatReordering(order_type, num_reduced_nodes, rowp, cols,
                           NULL, new_reduced_vars);

      // Set the new variable numbers for the boundary nodes
      // and include their offset
      for ( int i = 0, j = 0; i < numNodes; i++ ){
        if (reduced_nodes[i] >= 0){
          new_node_nums[i] = node_offset + new_reduced_vars[j];
          j++;
        }
      }

      delete [] new_reduced_vars;
      delete [] rowp;
      delete [] cols;
    }

    delete [] reduced_nodes;
  }

  // Find the values assigned to the nodes requested from external nodes
  // recv_nodes is now an outgoing list of nodes to other processes
  for ( int i = 0; i < nrecv_nodes; i++ ){
    int * item = (int*)bsearch(&recv_nodes[i], sorted_recv_nodes, nrecv_sorted,
			       sizeof(int), FElibrary::comparator);
    int index = item - sorted_recv_nodes;
    recv_nodes[i] = new_node_nums[sorted_recv_index[index]];
  }

  // Now send the new node numbers back to the other processors 
  // that reference them. This also uses all-to-all communication.
  int * sorted_new_extern_nodes = new int[ nextern_nodes ];
  MPI_Alltoallv(recv_nodes, recv_count, recv_ptr, MPI_INT,
		sorted_new_extern_nodes, extern_count, extern_ptr, MPI_INT, 
		tacs_comm);

  // Once the new node numbers from other processors is received
  // apply these new node numbers back to the locally owned
  // reference numbers.
  for ( int i = 0; i < numNodes; i++ ){
    if (tacsNodeNums[i] < ownerRange[mpiRank] ||
        tacsNodeNums[i] >= ownerRange[mpiRank+1]){
      // Find tacsNodeNums[i] in sorted_extern_nodes
      int * item = (int *) bsearch(&tacsNodeNums[i], sorted_extern_nodes, 
				   nextern_nodes, sizeof(int), 
				   FElibrary::comparator);
      if (item){
	new_node_nums[i] = sorted_new_extern_nodes[item-sorted_extern_nodes];
      }
      else {
	fprintf(stderr, "[%d] Error, could not find recvieved node\n", 
		mpiRank);
      }
    }
  }

  // Copy over the new node numbers
  memcpy(tacsNodeNums, new_node_nums, numNodes*sizeof(int));

  // This code checks to see if the nodes have been uniquely defined.
  // This is useful for checking that the nodes have been defined
  // correctly.
  int nnodes = FElibrary::uniqueSort(new_node_nums, numNodes);
  if (nnodes != numNodes){
    fprintf(stderr, "[%d] Error, nodes are not unique %d < %d \n", 
            mpiRank, nnodes, numNodes);
  }
  
  delete [] sorted_new_extern_nodes;
  delete [] sorted_extern_nodes;
  delete [] extern_ptr;
  delete [] extern_count;
  delete [] recv_count;
  delete [] recv_ptr;
  delete [] recv_nodes;
  delete [] sorted_recv_nodes;
  delete [] sorted_recv_index;
  delete [] new_node_nums;
}

/*
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

  P * A * P^{T} has fewer non-zeros.
  
  The function returns an array new_vars such that:

  new_vars = P^{T} old_vars
  perm = P
*/
void TACSAssembler::computeMatReordering( enum OrderingType order_type, 
                                          int nvars, int * rowp, int * cols,
                                          int * perm, int * new_vars ){
  int * _perm = perm;
  int * _new_vars = new_vars;

  if (!perm){ _perm = new int[ nvars ]; }
  if (!new_vars){ _new_vars = new int[ nvars ]; }
  
  if (order_type == RCM_ORDER){
    int root_node = 0;
    int n_rcm_iters = 1;
    matutils::ComputeRCMOrder(nvars, rowp, cols,
                              _new_vars, root_node, n_rcm_iters);

    if (perm){
      for ( int k = 0; k < nvars; k++ ){
        perm[_new_vars[k]] = k;    
      }
    }
  }
  else if (order_type == AMD_ORDER){
    double control[AMD_CONTROL], info[AMD_INFO];
    amd_defaults(control); // Use the default values
    amd_order(nvars, rowp, cols, _perm, 
              control, info);
    
    if (new_vars){
      for ( int k = 0; k < nvars; k++ ){
        new_vars[_perm[k]] = k;
      }
    }
  }
  else if (order_type == ND_ORDER){
    int numflag = 0, options[8] = { 0, 0, 0, 0,  
                                    0, 0, 0, 0 };
    METIS_NodeND(&nvars, rowp, cols, &numflag, 
                 options, _perm, _new_vars);    
  }
  else if (order_type == TACS_AMD_ORDER){
    int use_exact_degree = 0;
    int ncoupling_nodes = 0;
    int * coupling_nodes = NULL;
    amd_order_interface(nvars, rowp, cols, _perm, 
                        coupling_nodes, ncoupling_nodes,
                        use_exact_degree);

    if (new_vars){
      for ( int k = 0; k < nvars; k++ ){
        new_vars[_perm[k]] = k;
      }
    }
  }
  else if (order_type == NATURAL_ORDER){
    if (perm){
      for ( int k = 0; k < nvars; k++ ){
        perm[k] = k;
      }
    }
    if (new_vars){
      for ( int k = 0; k < nvars; k++ ){
        new_vars[k] = k;
      }
    }
  }

  if (!perm){ delete [] _perm; }
  if (!new_vars){ delete [] _new_vars; }
}

/*!
  The following function creates a data structure that links nodes to
  elements - this reverses the existing data structure that links
  elements to nodes but keeps the original in tact.

  The algorithm proceeds as follows:

  1. The size of the arrays are determined by finding how many nodes
  point to each element

  2. The index into the nodeElem array is determined by adding up the
  contributions from all previous entries.

  3. The original data structure is again traversed and this time an
  element number is associated with each element.
*/
void TACSAssembler::computeNodeToElementCSR( int **_nodeElemIndex, 
					     int **_nodeElem ){
  // Determine the element adjacency
  // Go through and add up each time an element refers to a node
  int * nodeElementIndex = new int[ numNodes+1 ];
  memset(nodeElementIndex, 0, (numNodes+1)*sizeof(int));
 
  for ( int i = 0; i < numElements; i++ ){
    int end = elementNodeIndex[i+1];
    for ( int j = elementNodeIndex[i]; j < end; j++ ){
      int varNum = elementLocalNodes[j];
      if (varNum >= 0){
	// This is an independent node
	nodeElementIndex[varNum+1]++;
      }
      else {
	// This is a dependent-node, determine which independent
	// nodes it depends on
	varNum = -varNum-1;
	int kend = depNodeIndex[varNum+1];
	for ( int k = depNodeIndex[varNum]; k < kend; k++ ){
	  int node = depNodeToLocal[k];
	  nodeElementIndex[node+1]++;
	}
      }
    }
  }

  // Sum up the total size of the array
  for ( int i = 0; i < numNodes; i++ ){
    nodeElementIndex[i+1] += nodeElementIndex[i];
  }
  
  // Set up the data structure that can determine the node -> element graph
  int size = nodeElementIndex[numNodes];
  int * nodeElements = new int[ size ];
  for ( int i = 0; i < numElements; i++ ){
    int end = elementNodeIndex[i+1];
    for ( int j = elementNodeIndex[i]; j < end; j++ ){
      int varNum = elementLocalNodes[j];
      if (varNum >= 0){
	nodeElements[nodeElementIndex[varNum]] = i;
	nodeElementIndex[varNum]++;
      }
      else {
	// This is a dependent-node, determine which independent
	// nodes it depends on
	varNum = -varNum-1;
	int kend = depNodeIndex[varNum+1];
	for ( int k = depNodeIndex[varNum]; k < kend; k++ ){
	  int node = depNodeToLocal[k];
	  nodeElements[nodeElementIndex[node]] = i;
	  nodeElementIndex[node]++;
	}
      }
    }
  }

  // Set up the pointer array which denotes the start (and end) of each node
  for ( int i = 0; i < numNodes; i++ ){
    nodeElementIndex[numNodes-i] = nodeElementIndex[numNodes-i-1];
  }
  nodeElementIndex[0] = 0;

  matutils::SortAndUniquifyCSR(numNodes, nodeElementIndex, nodeElements, 0);

  *_nodeElem = nodeElements;
  *_nodeElemIndex = nodeElementIndex;
}

/*!
  Set up a CSR data structure pointing from local nodes to other
  local nodes.

  input:
  nodiag = Remove the diagonal matrix entry 

  output:
  rowp = the row pointer corresponding to CSR data structure
  cols = the column indices for each row of the CSR data structure

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
*/
void TACSAssembler::computeLocalNodeToNodeCSR( int ** _rowp, int ** _cols, 
					       int nodiag ){
  int * cols = NULL;
  int * rowp = new int[ numNodes+1 ];
  memset(rowp, 0, (numNodes+1)*sizeof(int));

  // Get the node -> element data structure
  int * nodeElemIndex = NULL;
  int * nodeElems = NULL;
  computeNodeToElementCSR(&nodeElemIndex, &nodeElems);

  if (numDependentNodes > 0){
    // Count the number of nodes associated with each element
    int * nodeCount = new int[ numElements ];
    memset(nodeCount, 0, numElements*sizeof(int));

    for ( int i = 0; i < numElements; i++ ){
      int jend = elementNodeIndex[i+1];
      for ( int j = elementNodeIndex[i]; j < jend; j++ ){
	int node = elementLocalNodes[j];
	if (node >= 0){
	  nodeCount[i]++;
	}
	else {
	  // Find the dependent node
	  int dep = -node-1;
	  nodeCount[i] += depNodeIndex[dep+1] - depNodeIndex[dep];
	}
      }
    }

    // First, populate rowp by finding the approximate number of
    // independent nodes per element
    for ( int i = 0; i < numNodes; i++ ){
      for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	int elem = nodeElems[j];
	rowp[i+1] += nodeCount[elem];
      }
    }

    // Make a conservative estimate of the rowp pointer data
    for ( int i = 0; i < numNodes; i++ ){
      rowp[i+1] += rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[numNodes];
    cols = new int[ nnz ];
    for ( int i = 0; i < nnz; i++ ){
      cols[i] = -1;
    }

    // Add the element contribution to the column indices
    for ( int i = 0; i < numNodes; i++ ){
      for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	int elem = nodeElems[j];
	int kend = elementNodeIndex[elem+1];

	// Scan through all the nodes belonging to this element
	int r = rowp[i];
	for ( int k = elementNodeIndex[elem]; k < kend; k++ ){
	  int node = elementLocalNodes[k]; 
	  if (node >= 0){
	    // This is an independent node
	    cols[r] = node;  r++;
	  }
	  else {
	    // This is a dependent node, add the dependent node
	    // variables
	    int dep = -node-1;
	    int pend = depNodeIndex[dep+1];
	    for ( int p = depNodeIndex[dep]; p < pend; p++, r++ ){
	      cols[r] = depNodeToLocal[p];
	    }
	  }
	}
	rowp[i] = r;
      }
    }

    // Adjust rowp back to a zero-based index
    for ( int i = numNodes; i > 0; i-- ){
      rowp[i] = rowp[i-1];
    }
    rowp[0] = 0;

    delete [] nodeCount;
  }
  else {
    // First, populate rowp by adding the contribution to a node from
    // all adjacent elements.
    for ( int i = 0; i < numNodes; i++ ){
      for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	int elem = nodeElems[j];
	rowp[i+1] += elementNodeIndex[elem+1] - elementNodeIndex[elem];
      }
    }

    // Make a conservative estimate of rowp
    for ( int i = 0; i < numNodes; i++ ){
      rowp[i+1] += rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[numNodes];
    cols = new int[ nnz ];
    for ( int i = 0; i < nnz; i++ ){
      cols[i] = -1;
    }

    // Add the element contribution to the column indices
    for ( int i = 0; i < numNodes; i++ ){
      for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	int elem = nodeElems[j];     
	int r = rowp[i];
	for ( int k = elementNodeIndex[elem]; 
	      k < elementNodeIndex[elem+1]; k++, r++ ){
	  cols[r] = elementLocalNodes[k];
	}
	rowp[i] = r;
      }
    }

    // Adjust rowp back to a zero-based index
    for ( int i = numNodes; i > 0; i-- ){
      rowp[i] = rowp[i-1];
    }
    rowp[0] = 0;
  }

  // Go through and sort/uniquify each row and remove the diagonal if requested
  matutils::SortAndUniquifyCSR(numNodes, rowp, cols, nodiag);

  delete [] nodeElems;
  delete [] nodeElemIndex;

  *_rowp = rowp;
  *_cols = cols;
}

/*
  Prepare a reduced CSR data structure corresponding to a matrix
  formed from a selection of the global matrix. This routine can be
  used in matrix/variable re-ordering computations.

  input:
  nrnodes = the number of reduced nodes
  rnodes = the indices of the reduced nodes
  nodiag = flag to indicate whether to remove the diagonal matrix entry

  output:
  rowp = the row pointer corresponding to CSR data structure
  cols = the column indices for each row of the CSR data structure

  This function uses the same algorithm as computeLocalNodeToNodeCSR,
  but performs extra operations required to restrict the computations
  to Ar.  The rnodes array must consist of nrnodes non-negative
  integers between 0 and nrnodes-1, at any arbitrary location. All
  remaining entries of rnodes must be negative.  
*/
void TACSAssembler::computeLocalNodeToNodeCSR( int ** _rowp, int ** _cols, 
					       int nrnodes, const int * rnodes,
					       int nodiag ){
  int * cols = NULL;
  int * rowp = new int[ nrnodes+1 ];
  memset(rowp, 0, (nrnodes+1)*sizeof(int));

  int * nodeElems = NULL;
  int * nodeElemIndex = NULL;

  computeNodeToElementCSR(&nodeElemIndex, &nodeElems);

  if (numDependentNodes > 0){
    // Count the number of nodes associated with each element
    int * nodeCount = new int[ numElements ];
    memset(nodeCount, 0, numElements*sizeof(int));

    for ( int i = 0; i < numElements; i++ ){
      int jend = elementNodeIndex[i+1];
      for ( int j = elementNodeIndex[i]; j < jend; j++ ){
	int node = elementLocalNodes[j];
	if (node >= 0){
	  if (rnodes[node] >= 0){
	    nodeCount[i]++;
	  }
	}
	else {
	  // Find the dependent node
	  int dep = -node-1;
	  for ( int k = depNodeIndex[dep]; 
		k < depNodeIndex[dep+1]; k++ ){
	    if (rnodes[depNodeToLocal[k]] >= 0){
	      nodeCount[i]++;
	    }
	  }
	}
      }
    }

    // Count up the contribution to the rowp array from all elements
    // using the node->element data
    for ( int i = 0; i < numNodes; i++ ){
      int node = rnodes[i];
      if (node >= 0){      
	for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	  int elem = nodeElems[j];
	  rowp[node+1] += nodeCount[elem];
	}
      }
    }

    // Make a conservative estimate of rowp
    for ( int i = 0; i < nrnodes; i++ ){
      rowp[i+1] = rowp[i+1] + rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[nrnodes];
    cols = new int[ nnz ];

    // Add the element contribution to the column indices
    for ( int i = 0; i < numNodes; i++ ){
      int node = rnodes[i];
      if (node >= 0){
	for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	  int elem = nodeElems[j];
	  int kend = elementNodeIndex[elem+1];

	  // Scan through all the nodes belonging to this element
	  int r = rowp[node];
	  for ( int k = elementNodeIndex[elem]; k < kend; k++ ){
	    int local = elementLocalNodes[k]; 
	    if (local >= 0){
	      int rn = rnodes[local];
	      if (rn >= 0){
		// This is an independent node
		cols[r] = rn;  
		r++;
	      }
	    }
	    else {
	      // This is a dependent node, add the dependent node
	      // variables
	      int dep = -local-1;
	      int pend = depNodeIndex[dep+1];
	      for ( int p = depNodeIndex[dep]; p < pend; p++ ){
		int rn = rnodes[depNodeToLocal[p]];
		if (rn >= 0){
		  cols[r] = rn; 
		  r++;
		}
	      }
	    }
	  }
	  rowp[node] = r;
	}
      }
    }

    // Adjust rowp back to a zero-based index
    for ( int i = nrnodes; i > 0; i-- ){
      rowp[i] = rowp[i-1];
    }
    rowp[0] = 0;

    delete [] nodeCount;
  }
  else {
    // First, populate rowp by adding the contribution to a node from 
    // all adjacent elements.
    for ( int i = 0; i < numNodes; i++ ){
      int node = rnodes[i];
      if (node >= 0){      
	for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	  int elem = nodeElems[j];
	  // Count up the number of reduced nodes that are required here.
	  int count = 0;
	  for ( int k = elementNodeIndex[elem]; 
		k < elementNodeIndex[elem+1]; k++ ){
	    if (rnodes[elementLocalNodes[k]] >= 0){
	      count++;
	    }
	  }
	
	  rowp[node+1] += count;
	}
      }
    }

    // Make a conservative estimate of rowp
    for ( int i = 0; i < nrnodes; i++ ){
      rowp[i+1] = rowp[i+1] + rowp[i];
    }

    // Set up the column indices for each row
    int nnz = rowp[nrnodes];
    cols = new int[ nnz ];

    // Add the element contribution to the column indices
    for ( int i = 0; i < numNodes; i++ ){
      int node = rnodes[i];
      if (node >= 0){
	for ( int j = nodeElemIndex[i]; j < nodeElemIndex[i+1]; j++ ){
	  int elem = nodeElems[j];     
	  int r = rowp[node];
	  for ( int k = elementNodeIndex[elem]; 
		k < elementNodeIndex[elem+1]; k++ ){
	    int rn = rnodes[elementLocalNodes[k]];
	    if (rn >= 0){
	      cols[r] = rn;
	      r++;
	    }
	  }
	  rowp[node] = r;
	}
      }
    }

    // Adjust rowp back to a zero-based index
    for ( int i = nrnodes; i > 0; i-- ){
      rowp[i] = rowp[i-1];
    }
    rowp[0] = 0;
  }

  // Go through and sort/uniquify each row and remove the diagonal if requested
  matutils::SortAndUniquifyCSR(nrnodes, rowp, cols, nodiag);

  delete [] nodeElems;
  delete [] nodeElemIndex;

  *_rowp = rowp;
  *_cols = cols;
}

/*!  
  The function finalize performs a number of synchronization tasks
  that prepare the finite-element model for use.

  tacsNodeNums[i] is the global node number for the local node number i

  Two objects are required:
  1. VarMap is constructed with the block sizes of each 
  node owned by this process

  2. VecDistribute is constructed so that it takes an array and
  distributes its values to a vector or takes the vector values and
  collects them into an array This requires a sorted array of global
  node numbers.  
*/
void TACSAssembler::finalize(){
  if (meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call finalize() more than once!\n", 
	    mpiRank);
    return;
  }
  
  // Flag to indicate that we've initialized TACSAssembler -
  // the initialization can only be done once
  meshFinalizedFlag = 1;

  if (currElement < numElements || currNode < numNodes){
    fprintf(stderr, "[%d] Error: Insufficient definition of the \
elements and/or nodes\n", mpiRank);
    return;
  }
  else if (currElement != numElements || currNode != numNodes){
    fprintf(stderr, "[%d] Warning: The elements or nodes \
were defined more than once\n", mpiRank);
  }
  if (numDependentNodes > 0 && 
      !((depNodeIndex && depNodeToLocal) && depNodeWeights)){
    fprintf(stderr, "[%d] Error: Dependent nodes not defined\n",
	    mpiRank);
    return;
  }

  // Determine the maximum number of nodes per element
  maxElementNodes = 0;
  for ( int i = 0; i < numElements; i++ ){
    int size = elementNodeIndex[i+1] - elementNodeIndex[i];
    if (size > maxElementNodes){
      maxElementNodes = size;
    }
  }

  const int * ownerRange;
  int rank, size;
  varMap->getOwnerRange(&ownerRange, &rank, &size);

  for ( int i = 0; i < numNodes; i++ ){
    if (tacsNodeNums[i] < 0 ||
        tacsNodeNums[i] >= ownerRange[mpiSize]){
      fprintf(stderr, "[%d] tacsNodeNums[%d] = %d out of range [0,%d)\n",
	      mpiRank, i, tacsNodeNums[i], ownerRange[mpiSize]);
    }
  }

  // Construct the elementTacsNodes pointer. This requires more memory,
  // but is faster than the double-indexing approach required if we
  // were to just use elementLocalNodes
  elementTacsNodes = new int[ elementNodeIndex[numElements] ];
  for ( int i = 0; i < elementNodeIndex[numElements]; i++ ){
    if (elementLocalNodes[i] >= 0 && elementLocalNodes[i] < numNodes){
      elementTacsNodes[i] = tacsNodeNums[ elementLocalNodes[i] ];
    }
    else if (elementLocalNodes[i] < 0 && 
	     elementLocalNodes[i] >= -numDependentNodes){
      elementTacsNodes[i] = elementLocalNodes[i];
    }
    else {
      fprintf(stderr, "[%d] elementLocalNodes[%d] = %d value is out of range\n",
	      mpiRank, i, elementLocalNodes[i]);
    }
  }

  // Set up data for any dependent nodes. Note that the minimum
  // number of independent nodes is set as the maximum number
  // of element node by default. This is required for addMatValues()
  // to have enough memory for TRANSPOSE matrix assembly.
  maxElementIndepNodes = maxElementNodes;
  if (numDependentNodes > 0){
    depNodeToTacs = new int[ depNodeIndex[numDependentNodes] ];
    for ( int i = 0; i < depNodeIndex[numDependentNodes]; i++ ){
      if (depNodeToLocal[i] >= 0 && depNodeToLocal[i] < numNodes){
	depNodeToTacs[i] = tacsNodeNums[ depNodeToLocal[i] ];
      }
      else {
	fprintf(stderr, "[%d] depNodeToLocal[%d] = %d values is out of range\n",
		mpiRank, i, depNodeToLocal[i]);
      }
    }

    // Compute the maximum number of independent nodes
    for ( int i = 0; i < numElements; i++ ){
      int jend = elementNodeIndex[i+1];
      int nnodes = 0;
      for ( int j = elementNodeIndex[i]; j < jend; j++ ){
	if (elementLocalNodes[j] >= 0){
	  nnodes++;
	}
	else {
	  int dep = -elementLocalNodes[j]-1;
	  nnodes += depNodeIndex[dep+1] - depNodeIndex[dep];
	}
      }
      if (nnodes > maxElementIndepNodes){
	maxElementIndepNodes = nnodes;
      }
    }
  }

  // Create the distribution between the local nodes and the global ones
  vecDistIndices = new BVecIndices(&tacsNodeNums, numNodes);
  vecDistIndices->incref();
  vecDistIndices->getIndices(&tacsNodeNums);
  
  vecDist = new BVecDistribute(varMap, vecDistIndices);
  vecDist->incref();

  // Allocate space that will be used in the analysis
  initializeArrays();
}

/*!
  Collect all the design variable values assigned by this process

  This code does not ensure consistency of the design variable values
  between processes. If the values of the design variables are
  inconsistent to begin with, the maximum design variable value is
  returned. Call setDesignVars to make them consistent.

  Each process contains objects that maintain their own design
  variable values. Ensuring the consistency of the ordering is up to
  the user. Having multiply-defined design variable numbers
  corresponding to different design variables results in undefined
  behaviour.

  dvs:    the array of design variable values (output)
  numDVs: the number of design variables
*/
void TACSAssembler::getDesignVars( TacsScalar dvs[], int numDVs ){
  TacsScalar * tempDVs = new TacsScalar[ numDVs ];
  memset(tempDVs, 0, numDVs*sizeof(TacsScalar));

  // Get the design variables from the elements on this process 
  for ( int i = 0; i < numElements; i++ ){
    elements[i]->getDesignVars(tempDVs, numDVs);
  }
  
  // Get the design variables from the auxiliary elements
  if (aux_elements){
    aux_elements->getDesignVars(tempDVs, numDVs);
  }

  MPI_Allreduce(tempDVs, dvs, numDVs, TACS_MPI_TYPE, 
		TACS_MPI_MAX, tacs_comm);
  
  // Free the allocated array
  delete [] tempDVs;
}

/*!
  Set the design variables.

  The design variable values provided must be the same on all
  processes for consistency. This call however, is not collective.

  dvs:    the array of design variable values
  numDVs: the number of design variables
*/
void TACSAssembler::setDesignVars( const TacsScalar dvs[], int numDVs ){
  for ( int i = 0; i < numElements; i++ ){
    elements[i]->setDesignVars(dvs, numDVs);
  }

  // Set the design variables in the auxiliary elements
  if (aux_elements){
    aux_elements->setDesignVars(dvs, numDVs);
  }
}

/*
  Retrieve the design variable range.

  This call is collective on all TACS processes. The ranges provided
  by indivdual objects may not be consistent (if someone provided
  incorrect data they could be.) Make a best guess; take the minimum
  upper bound and the maximum lower bound.

  lowerBound: the lower bound on the design variables (output)
  upperBound: the upper bound on the design variables (output)
  numDVs:     the number of design variables
*/
void TACSAssembler::getDesignVarRange( TacsScalar lowerBound[], 
				       TacsScalar upperBound[], 
				       int numDVs ){
  TacsScalar * tempLB = new TacsScalar[ numDVs ];
  TacsScalar * tempUB = new TacsScalar[ numDVs ];
  memset(tempLB, 0, numDVs*sizeof(TacsScalar));
  memset(tempUB, 0, numDVs*sizeof(TacsScalar));

  // Get the design variables from the elements on this process 
  for ( int i = 0; i < numElements; i++ ){
    elements[i]->getDesignVarRange(tempLB, tempUB, numDVs);
  }

  // Get the design variable range from the auxiliary elements
  if (aux_elements){
    aux_elements->getDesignVarRange(tempLB, tempUB, numDVs);
  }

  // Take the max of the lower bounds
  MPI_Allreduce(tempLB, lowerBound, numDVs, TACS_MPI_TYPE, 
		TACS_MPI_MAX, tacs_comm); 

  // Take the min of the upper bounds
  MPI_Allreduce(tempUB, upperBound, numDVs, TACS_MPI_TYPE, 
		TACS_MPI_MIN, tacs_comm); 

  delete [] tempLB;
  delete [] tempUB;
}			  
/*!
  Set the number of threads to use in the computation
*/
void TACSAssembler::setNumThreads( int t ){
  thread_info->setNumThreads(t);  
}

/*
  Retrieve the local variable arrays
*/
void TACSAssembler::getLocalArrays( const TacsScalar **_Xpts,
                                    TacsScalar **_localRes,
                                    const TacsScalar **_localVars,
                                    const TacsScalar **_localDotVars,
                                    const TacsScalar **_localDDotVars ){
  if (_Xpts){ *_Xpts = Xpts; }
  if (_localRes){ *_localRes = localRes; }
  if (_localVars){ *_localVars = localVars; }
  if (_localDotVars){ *_localDotVars = localDotVars; }
  if (_localDDotVars){ *_localDDotVars = localDDotVars; }
}

/*!
  Initialize several arrays that will be used during the course 
  of the analysis
*/
void TACSAssembler::initializeArrays(){
  // Allocate memory for the variable values
  int size = varsPerNode*(numNodes + numDependentNodes);
  localRes = new TacsScalar[ size ];
  memset(localRes, 0, size*sizeof(TacsScalar));

  localVars = new TacsScalar[ size ];
  localDotVars = new TacsScalar[ size ];
  localDDotVars = new TacsScalar[ size ];
  memset(localVars, 0, size*sizeof(TacsScalar));
  memset(localDotVars, 0, size*sizeof(TacsScalar));
  memset(localDDotVars, 0, size*sizeof(TacsScalar));

  // Allocate memory for the working array:
  // Determine the size of the data working array
  // max requirement is 4 element variable-size arrays,
  // 2 node-size arrays and either the element matrix or
  // the derivative of the residuals w.r.t. the nodes.
  int dataSize = maxElementIndepNodes + 4*maxElementSize + 
    2*TACS_SPATIAL_DIM*maxElementNodes;
  if (TACS_SPATIAL_DIM*maxElementNodes > maxElementSize){
    dataSize += TACS_SPATIAL_DIM*maxElementNodes*maxElementSize; 
  }
  else {
    dataSize += maxElementSize*maxElementSize;
  }
  elementData = new TacsScalar[ dataSize ];

  int idataSize = maxElementIndepNodes + maxElementNodes+1;
  elementIData = new int[ idataSize ];
}

/*
  Get pointers to the element data. This code provides a way to
  automatically segment an array to avoid coding mistakes.

  Note that this is coded in such a way that you can provide NULL
  arguments to
*/
void TACSAssembler::getDataPointers( TacsScalar *data,
				     TacsScalar **v1, TacsScalar **v2, 
				     TacsScalar **v3, TacsScalar **v4,
				     TacsScalar **x1, TacsScalar **x2, 
				     TacsScalar **weights,
				     TacsScalar **mat ){
  int s = 0;
  if (v1){ *v1 = &data[s];  s += maxElementSize; }
  if (v2){ *v2 = &data[s];  s += maxElementSize; }
  if (v3){ *v3 = &data[s];  s += maxElementSize; }
  if (v4){ *v4 = &data[s];  s += maxElementSize; }
  if (x1){ *x1 = &data[s];  s += TACS_SPATIAL_DIM*maxElementNodes; };
  if (x2){ *x2 = &data[s];  s += TACS_SPATIAL_DIM*maxElementNodes; };
  if (weights){ *weights = &data[s];  s += maxElementIndepNodes; }
  if (mat){
    *mat = &data[s];
  }
}

/*
  Create a distributed vector.
  
  Vector classes initialized by one TACS object, cannot be used by a
  second, unless they share are exactly the parallel layout.
*/
BVec * TACSAssembler::createVec(){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call createVec() before finalize()\n", 
	    mpiRank);
    return NULL;
  }

  return new BVec(varMap, varsPerNode, bcMap);
}

/*!
  Create a distributed matrix.
*/
DistMat * TACSAssembler::createMat(){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call createMat() before finalize()\n", 
	    mpiRank);
    return NULL;
  }
  
  // Compute the local connectivity
  int * rowp, * cols;
  computeLocalNodeToNodeCSR(&rowp, &cols);
  
  // Create the distributed matrix class
  DistMat * dmat = new DistMat(thread_info, varMap, varsPerNode,
                               numNodes, rowp, cols,
			       vecDistIndices, bcMap);

  // Free the local connectivity
  delete [] rowp;
  delete [] cols;

  // Return the resulting matrix object
  return dmat;
}

/*!  
  Create a parallel matrix specially suited for finite-element
  analysis.
  
  On the first call, this computes a reordering with the scheme
  provided. On subsequent calls, the reordering scheme is reused so
  that all FEMats, created from the same TACSAssembler object have
  the same non-zero structure. This makes adding matrices together
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
  matrix reordering is stored in feMatBIndices and feMatCIndices while
  two mapping objects are created that map the variables from the
  global vector to reordered matrix.

  Mathematically this reordering can be written as follows,

  A' = (P A P^{T})

  where P^{T} is a permutation of the columns (variables), while P is
  a permutation of the rows (equations).
*/
FEMat * TACSAssembler::createFEMat( enum OrderingType order_type ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call createFEMat() before finalize()\n", 
	    mpiRank);
    return NULL;
  }
  if (order_type == NATURAL_ORDER){
    fprintf(stderr, 
	    "[%d] Cannot call createFEMat() with \
order_type == NATURAL_ORDER\n",
	    mpiRank);
    order_type = TACS_AMD_ORDER;
  }

  if (!feMatBMap){   
    int nlocal_nodes = numNodes;
    int ncoupling_nodes = 0;
    int * perm_local_nodes = NULL;
    int * tacs_local_nodes = NULL;

    int * perm_coupling_nodes = NULL;
    int * tacs_coupling_nodes = NULL;

    if (order_type == TACS_AMD_ORDER){
      /*
        This is the new, coupled code. Use the AMD ordering scheme in
        TACS to compute an ordering of the nodes that reduces the
        fill-in in the complete matrix -- including the off-diagonal
        contributions. This can reduce the computational effort and
        the discrepancy between factorization times on different
        processes.
      */

      // Find the local node numbers of the coupling nodes.
      // Note that this is a sorted list
      int * coupling_nodes;
      ncoupling_nodes = computeCouplingNodes(&coupling_nodes);
      nlocal_nodes = numNodes - ncoupling_nodes;
      
      int *rowp, *cols;
      int no_diagonal = 1;
      computeLocalNodeToNodeCSR(&rowp, &cols, 
                                no_diagonal);
      
      // Here perm is the entire permutation array
      int * perm = new int[ numNodes ];
      int use_exact_degree = 0; // Don't use the exact degree
      amd_order_interface(numNodes, rowp, cols, perm, 
                          coupling_nodes, ncoupling_nodes, use_exact_degree);
      
      delete [] rowp;
      delete [] cols; 

      perm_coupling_nodes = new int[ ncoupling_nodes ];
      tacs_coupling_nodes = new int[ ncoupling_nodes ];
      
      perm_local_nodes = new int[ nlocal_nodes ];
      tacs_local_nodes = new int[ nlocal_nodes ];
      
      for ( int i = 0; i < nlocal_nodes; i++ ){
        perm_local_nodes[i] = perm[i];
        tacs_local_nodes[i] = tacsNodeNums[perm_local_nodes[i]];
      }

      for ( int i = 0; i < ncoupling_nodes; i++ ){
        perm_coupling_nodes[i] = perm[i + nlocal_nodes];
        tacs_coupling_nodes[i] = tacsNodeNums[perm_coupling_nodes[i]];
      }

      delete [] perm;
      delete [] coupling_nodes;
    }
    else {
      /* 
         This is the old scheme. Order the coupling nodes first, then
         order the local nodes. This ignores the off-diagonal fill-ins
         which can be considerable!
      */
      int no_diagonal = 0;
      if (order_type == ND_ORDER){
        no_diagonal = 1;
      }

      // Find the local node numbers of the coupling nodes.
      // Note that this is a sorted list
      int * coupling_nodes;
      ncoupling_nodes = computeCouplingNodes(&coupling_nodes);
      nlocal_nodes = numNodes - ncoupling_nodes;
      
      // Set the coupling nodes for ordering
      int * all_nodes = new int[ numNodes ];
      for ( int k = 0; k < numNodes; k++ ){
        all_nodes[k] = -1;
      }
      
      for ( int k = 0; k < ncoupling_nodes; k++ ){
        all_nodes[coupling_nodes[k]] = k;
      }
      
      perm_coupling_nodes = new int[ ncoupling_nodes ];    
      tacs_coupling_nodes = new int[ ncoupling_nodes ];
      
      // Now, compute the reordering for the local coupling variables
      if (ncoupling_nodes > 0){
        int *rowp, *cols;
        computeLocalNodeToNodeCSR(&rowp, &cols, 
                                  ncoupling_nodes, all_nodes,
                                  no_diagonal);
        
        // Compute the permutation of the coupling nodes
        computeMatReordering(order_type, ncoupling_nodes, rowp, cols, 
                             perm_coupling_nodes, NULL);
        
        for ( int i = 0; i < ncoupling_nodes; i++ ){
          // Permute the coupling_nodes array - store in perm_coupling_nodes
          perm_coupling_nodes[i] = coupling_nodes[perm_coupling_nodes[i]];
          tacs_coupling_nodes[i] = tacsNodeNums[perm_coupling_nodes[i]];
        }
        
        delete [] rowp;
        delete [] cols;
      }
      
      // Set the remaining, local nodes for coupling
      perm_local_nodes = new int[ nlocal_nodes ];
      tacs_local_nodes = new int[ nlocal_nodes ];
      int * local_nodes = new int[ nlocal_nodes ];
      for ( int j = 0, k = 0; k < numNodes; k++ ){
        if (all_nodes[k] < 0){
          all_nodes[k] = j;
          local_nodes[j] = k;
          j++;
        }
        else {
          all_nodes[k] = -1;
        }
      }

      // Now, compute the reordering for the local variables
      int *rowp, *cols;    
      computeLocalNodeToNodeCSR(&rowp, &cols, 
                                nlocal_nodes, all_nodes,
                                no_diagonal);
      
      computeMatReordering(order_type, nlocal_nodes, rowp, cols,
                           perm_local_nodes, NULL);
      
      for ( int i = 0; i < nlocal_nodes; i++ ){
        // Permute the local nodes and record the corresponding tacs variables
        perm_local_nodes[i] = local_nodes[perm_local_nodes[i]];
        tacs_local_nodes[i] = tacsNodeNums[perm_local_nodes[i]];
      }
      
      delete [] rowp;
      delete [] cols;
      delete [] coupling_nodes;
      delete [] all_nodes;
      delete [] local_nodes;    
    }

    // Create persistent objects so that all further FEMats will have
    // the same ordering.
    feMatBIndices = new BVecIndices(&perm_local_nodes, nlocal_nodes);
    feMatCIndices = new BVecIndices(&perm_coupling_nodes, ncoupling_nodes);
    feMatBIndices->incref();
    feMatCIndices->incref();

    BVecIndices * tlocal = new BVecIndices(&tacs_local_nodes, nlocal_nodes);
    BVecIndices * tcoupling = new BVecIndices(&tacs_coupling_nodes, 
                                              ncoupling_nodes);
    feMatBMap = new BVecDistribute(varMap, tlocal);
    feMatCMap = new BVecDistribute(varMap, tcoupling);
    feMatBMap->incref();
    feMatCMap->incref();
  }

  // Compute he local non-zero pattern
  int * rowp, * cols;
  computeLocalNodeToNodeCSR(&rowp, &cols);

  FEMat * fmat = new FEMat(thread_info, varMap, 
                           varsPerNode, numNodes, rowp, cols,
                           feMatBIndices, feMatBMap,
                           feMatCIndices, feMatCMap, bcMap);

  delete [] rowp;
  delete [] cols;

  return fmat;
}

/*
  Set the dependent variable values based on the independent variable
  values. This must be called after the local independent nodal values
  are set. Note that this is a matrix-multiply.
*/
void TACSAssembler::setDependentVariables( const int perNode,
					   TacsScalar * vars ){
  if (numDependentNodes > 0){
    int offset = perNode*numNodes;
    
    for ( int i = 0; i < numDependentNodes; i++, offset += perNode ){
      int jend = depNodeIndex[i+1];

      // Zero the variables
      for ( int k = 0; k < perNode; k++ ){
	vars[offset + k] = 0.0;
      }
      
      // Compute the weighted value of the dependent node
      for ( int j = depNodeIndex[i]; j < jend; j++ ){
	int var = perNode*depNodeToLocal[j];
	
	for ( int k = 0; k < perNode; k++ ){
	  vars[offset + k] += depNodeWeights[j]*vars[var + k];
	}
      }
    }
  }
}

/*
  Add the residual contribution that is collected from the dependent
  nodes to the independent nodes that they depend on. Note that this
  is a transpose matrix-multiply add operation.
*/
void TACSAssembler::addDependentResidual( const int perNode,
					  TacsScalar * vars ){
  if (numDependentNodes > 0){
    int offset = perNode*numNodes;
    
    for ( int i = 0; i < numDependentNodes; i++, offset += perNode ){
      int jend = depNodeIndex[i+1];

      // Compute the weighted value of the dependent node
      for ( int j = depNodeIndex[i]; j < jend; j++ ){
	int var = perNode*depNodeToLocal[j];
	
	for ( int k = 0; k < perNode; k++ ){
	  vars[var + k] += depNodeWeights[j]*vars[offset + k];
	}
      }
    }
  }
}

/*
  Retrieve the initial conditions associated with the problem
*/
void TACSAssembler::getInitConditions( BVec *vars, BVec *dvars ){
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localVars, 0, size*sizeof(TacsScalar));
  memset(localDotVars, 0, size*sizeof(TacsScalar));

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemDVars, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemDVars, NULL, NULL,
		  &elemXpts, NULL, NULL, NULL);

  // Retrieve the initial condition values from each element
  for ( int i = 0; i < numElements; i++ ){    
    getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
    elements[i]->getInitCondition(elemVars, elemDVars, elemXpts);

    // Set the values into the array
    setValues(varsPerNode, i, elemVars, localVars);
    setValues(varsPerNode, i, elemDVars, localDotVars);
  }

  // Set the dependent node values internally
  setDependentVariables(varsPerNode, localVars);
  setDependentVariables(varsPerNode, localDotVars);

  // Set the variable values
  vecDist->beginReverse(localVars, vars, 
			BVecDistribute::INSERT);
  vecDist->endReverse(localVars, vars,
		      BVecDistribute::INSERT);

  // Set the variable derivatives
  vecDist->beginReverse(localDotVars, dvars, 
			BVecDistribute::INSERT);
  vecDist->endReverse(localDotVars, dvars,
		      BVecDistribute::INSERT);
}

/*
  Zero the entries of the local variables
*/
void TACSAssembler::zeroVariables(){
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localVars, 0, size*sizeof(TacsScalar));
}

/*
  Zero the values of the time-derivatives of the state variables.
  This time-derivative is load-case independent.
*/
void TACSAssembler::zeroDotVariables(){
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localDotVars, 0, size*sizeof(TacsScalar));
}

/*
  Zero the values of the time-derivatives of the state variables.
  This time-derivative is load-case independent.
*/
void TACSAssembler::zeroDDotVariables(){
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localDDotVars, 0, size*sizeof(TacsScalar));
}

/*
  Set the value of the time/variables/time derivatives simultaneously
*/
void TACSAssembler::setVariables( double _time, 
                                  BVec *q, BVec *qdot, BVec *qddot ){
  time = _time;
  /*
  vecDist->beginForward(q, qdot, qddot, 
                        localVars, localDotVars, localDDotVars); 
  vecDist->endForward(q, qdot, qddot, 
                      localVars, localDotVars, localDDotVars); 
  */
  setDependentVariables(varsPerNode, localVars);
  setDependentVariables(varsPerNode, localDotVars);
  setDependentVariables(varsPerNode, localDDotVars);
}

/*
  Set the values of the state variables for the given load case
*/
void TACSAssembler::setVariables( BVec * stateVars ){
  vecDist->beginForward(stateVars, localVars);
  vecDist->endForward(stateVars, localVars);

  // Set the dependent node values
  setDependentVariables(varsPerNode, localVars);
}

/*
  Get the values of the state variables for the given load case
*/
void TACSAssembler::getVariables( BVec * stateVars ){
  vecDist->beginReverse(localVars, stateVars, 
			BVecDistribute::INSERT);
  vecDist->endReverse(localVars, stateVars, 
		      BVecDistribute::INSERT);
}

/*
  Set the values of the time-deriative of the state variables.

  This time-derivative is load-case independent.
*/
void TACSAssembler::setDotVariables( BVec * stateVars ){
  vecDist->beginForward(stateVars, localDotVars);
  vecDist->endForward(stateVars, localDotVars);

  // Set the dependent node values
  setDependentVariables(varsPerNode, localDotVars);
}

/*
  Set the values of the time-deriative of the state variables.

  This time-derivative is load-case independent.
*/
void TACSAssembler::setDDotVariables( BVec * stateVars ){
  vecDist->beginForward(stateVars, localDDotVars);
  vecDist->endForward(stateVars, localDDotVars);

  // Set the dependent node values
  setDependentVariables(varsPerNode, localDDotVars);
}

/*!
  Given an array of local node numbers, return the associated TACS
  numbers in the same array - over-writing the values.

  localNodes: the local node numbers to retrieve
  numNodse:   the number of nodes to retrieve
*/
void TACSAssembler::getTacsNodeNums( int localNodes[], int _numNodes ){
  for ( int i = 0; i < _numNodes; i++ ){
    if (localNodes[i] >= 0 && 
        localNodes[i] <  numNodes){
      localNodes[i] = tacsNodeNums[localNodes[i]];
    }
  }
}

/*
  Set the simulation time internally in the TACSAssembler object
*/
void TACSAssembler::setSimulationTime( double _time ){
  time = _time;
}

/*
  Retrieve the simulation time from the TACSAssembler object
*/
double TACSAssembler::getSimulationTime(){
  return time;
}

/*!
  Return the array of tacsNodeNums 

  tacsNodeNums: constant pointer to the array of all TACS variable numbers
*/
int TACSAssembler::getTacsNodeNums( const int ** _tacsNodeNums ){ 
  if (meshFinalizedFlag){
    *_tacsNodeNums = tacsNodeNums; 
    return numNodes;
  }
  
  *_tacsNodeNums = NULL;
  return 0;
}

/*!  
  Assemble the residual associated with the input load case.  
  
  This residual includes the contributions from element tractions set
  in the TACSSurfaceTraction class and any point loads. Note that the
  vector entries are zeroed first, and that the Dirichlet boundary
  conditions are applied after the assembly of the residual is
  complete.
  
  rhs:      the residual output
*/
void TACSAssembler::assembleRes( BVec *residual ){
  // Assemble the residual without applying boundary conditions
  assembleResNoBCs(residual);

  // Set the boundary conditions
  applyBCs(residual, localVars);
}

/*!
  Assemble the residuals of the finite-element problem but do not
  apply the boundary conditions. This is useful for determing the
  reaction forces.

  output:
  residual:      the residual of the governing equations
*/
void TACSAssembler::assembleResNoBCs( BVec *residual ){
  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)
  if (aux_elements){
    aux_elements->sort();
  }

  // Zero the residual
  residual->zeroEntries();

  // Zero the local residual values
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localRes, 0, size*sizeof(TacsScalar));

  if (thread_info->getNumThreads() > 1){
    // Set the number of completed elements to zero
    numCompletedElements = 0;
    tacsPInfo->tacs = this;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::assembleRes_thread,
                     (void*)tacsPInfo);
    }

    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  }
  else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
    getDataPointers(elementData, 
		    &vars, &dvars, &ddvars, &elemRes,
		    &elemXpts, NULL, NULL, NULL);

    // Get the auxiliary elements
    int naux = 0, aux_count = 0;
    TACSAuxElem *aux = NULL;
    if (aux_elements){
      naux = aux_elements->getAuxElements(&aux);
    }

    // Go through and add the residuals from all the elements
    for ( int i = 0; i < numElements; i++ ){    
      getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars, vars);
      getValues(varsPerNode, i, localDotVars, dvars);
      getValues(varsPerNode, i, localDDotVars, ddvars);

      // Add the residual from the working element
      int nvars = elements[i]->numVariables();
      memset(elemRes, 0, nvars*sizeof(TacsScalar));
      elements[i]->addResidual(time, elemRes, elemXpts, 
                               vars, dvars, ddvars);
      
      // Add the residual from any auxiliary elements
      while (aux_count < naux && aux[aux_count].num == i){
        aux[aux_count].elem->addResidual(time, elemRes, elemXpts,
                                         vars, dvars, ddvars);
        aux_count++;
      }

      // Add the residual values
      addValues(varsPerNode, i, elemRes, localRes);
    }
  }

  // Add the dependent-variable residual from the dependent nodes
  addDependentResidual(varsPerNode, localRes);

  // Assemble the full residual
  vecDist->beginReverse(localRes, residual, BVecDistribute::ADD);
  vecDist->endReverse(localRes, residual, BVecDistribute::ADD);
}

/*
  Evaluates the total kinetic and potential energies of the structure
*/
void TACSAssembler::evalEnergies( TacsScalar *Te, TacsScalar *Pe ){
  // Zero the kinetic and potential energy
  *Te = 0.0;
  *Pe = 0.0;

  if (thread_info->getNumThreads() > 1){
    fprintf(stderr, "[%d] Cannot evaluate energies: UNIMPLEMENTED \n", mpiRank);
    return;
  } 
  else {
    // Array for storing local kinetic and potential energies
    TacsScalar elem_energies[2] = {0.0, 0.0};
 
    // Retrieve pointers to temporary storage
    TacsScalar *elem_vars, *elem_dvars, *elem_xpts;
    getDataPointers(elementData, &elem_vars, &elem_dvars, 
                    NULL, NULL, &elem_xpts, NULL, NULL, NULL);

    // Loop over all elements and add individual contributions to the
    // total energy
    for ( int i = 0; i < numElements; i++ ){
      // Determine the values of the state variables for the current
      // element
      getValues(TACS_SPATIAL_DIM, i, Xpts, elem_xpts);
      getValues(varsPerNode, i, localVars, elem_vars);
      getValues(varsPerNode, i, localDotVars, elem_dvars);

      // Compute and add the element's contributions to the total
      // energy
      TacsScalar elemTe, elemPe;
      elements[i]->computeEnergies(time, &elemTe, &elemPe,
                                   elem_xpts, elem_vars, elem_dvars);

      // Add up the kinetic and potential energy
      *Te += elemTe;
      *Pe += elemPe;
    }

    // Sum up the kinetic and potential energies across all processors
    TacsScalar input[2], output[2];
    input[0] = *Te;
    input[1] = *Pe;    
    MPI_Allreduce(input, output, 2, TACS_MPI_TYPE, 
                  MPI_SUM, tacs_comm);

    *Te = output[0];
    *Pe = output[1];
  } 
} 
  
/*!
  Assemble the Jacobian matrix

  This function assembles the global Jacobian matrix and
  residual. This Jacobian includes the contributions from all
  elements. The Dirichlet boundary conditions are applied to the
  matrix by zeroing the rows of the matrix associated with a boundary
  condition, and setting the diagonal to unity. The matrix assembly
  also performs any communication required so that the matrix can be
  used immediately after assembly.

  residual:  the residual of the governing equations
  A:         the Jacobian matrix
  alpha:     coefficient on the variables
  beta:      coefficient on the time-derivative terms
  gamma:     coefficient on the second time derivative term
  matOr:     the matrix orientation NORMAL or TRANSPOSE
*/
void TACSAssembler::assembleJacobian( BVec *residual, 
				      TACSMat * A,
				      double alpha, double beta, double gamma,
				      MatrixOrientation matOr ){
  // Zero the residual and the matrix
  if (residual){ 
    int size = varsPerNode*(numNodes + numDependentNodes);
    memset(localRes, 0, size*sizeof(TacsScalar));
    residual->zeroEntries(); 
  }
  A->zeroEntries();

  // Sort the list of auxiliary elements - this call only performs the
  // sort if it is required (if new elements are added)
  if (aux_elements){
    aux_elements->sort();
  }

  // Run the p-threaded version of the assembly code
  if (thread_info->getNumThreads() > 1){
    // Set the number of completed elements to zero
    numCompletedElements = 0;
    tacsPInfo->tacs = this;
    tacsPInfo->mat = A;
    tacsPInfo->alpha = alpha;
    tacsPInfo->beta = beta;
    tacsPInfo->gamma = gamma;
    tacsPInfo->matOr = matOr;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::assembleJacobian_thread,
                     (void*)tacsPInfo);
    }

    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  }
  else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *dvars, *ddvars, *elemRes, *elemXpts;
    TacsScalar *elemWeights, *elemMat;
    getDataPointers(elementData, 
		    &vars, &dvars, &ddvars, &elemRes,
		    &elemXpts, NULL, &elemWeights, &elemMat);

    // Set the data for the auxiliary elements - if there are any
    int naux = 0, aux_count = 0;
    TACSAuxElem *aux = NULL;
    if (aux_elements){
      naux = aux_elements->getAuxElements(&aux);
    }

    for ( int i = 0; i < numElements; i++ ){
      getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars, vars);
      getValues(varsPerNode, i, localDotVars, dvars);
      getValues(varsPerNode, i, localDDotVars, ddvars);

      // Get the number of variables from the element
      int nvars = elements[i]->numVariables();

      // Compute and add the contributions to the residual
      if (residual){
        memset(elemRes, 0, nvars*sizeof(TacsScalar));
        elements[i]->addResidual(time, elemRes, elemXpts, 
                                 vars, dvars, ddvars);
      }

      // Compute and add the contributions to the Jacobian
      memset(elemMat, 0, nvars*nvars*sizeof(TacsScalar));
      elements[i]->addJacobian(time, elemMat, alpha, beta, gamma,
			       elemXpts, vars, dvars, ddvars);

      // Add the contribution to the residual and the Jacobian
      // from the auxiliary elements - if any
      while (aux_count < naux && aux[aux_count].num == i){
        if (residual){
          aux[aux_count].elem->addResidual(time, elemRes, elemXpts,
                                           vars, dvars, ddvars);
        }
        aux[aux_count].elem->addJacobian(time, elemMat, 
                                         alpha, beta, gamma,
                                         elemXpts, vars, dvars, ddvars);        
        aux_count++;
      }

      if (residual){
        addValues(varsPerNode, i, elemRes, localRes);
      }
      addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
    }
  }

  // Add the dependent-residual terms
  if (residual){ addDependentResidual(varsPerNode, localRes); }

  // Do any matrix and residual assembly if required
  A->beginAssembly();
  if (residual){
    vecDist->beginReverse(localRes, residual, BVecDistribute::ADD);
  }

  A->endAssembly();
  if (residual){
    vecDist->endReverse(localRes, residual, BVecDistribute::ADD);
  }

  // Apply the boundary conditions
  if (residual){ applyBCs(residual, localVars); }
  A->applyBCs();
}

/*!  
  Assemble a matrix of a specified type. Note that all matrices
  created from the TACSAssembler object have the same non-zero pattern
  and are interchangable.

  A:            the matrix to assemble (output)
  matType:      the matrix type defined in Element.h
  matOr:        the matrix orientation: NORMAL or TRANSPOSE
*/
void TACSAssembler::assembleMatType( ElementMatrixType matType,
                                     TACSMat *A, 
                                     MatrixOrientation matOr ){
  // Zero the matrix
  A->zeroEntries();

  if (thread_info->getNumThreads() > 1){
    // Set the number of completed elements to zero
    numCompletedElements = 0;    
    tacsPInfo->tacs = this;
    tacsPInfo->mat = A;
    tacsPInfo->matType = matType;
    tacsPInfo->matOr = matOr;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::assembleMatType_thread,
                     (void*)tacsPInfo);
    }

    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  }
  else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *elemXpts, *elemMat, *elemWeights;
    getDataPointers(elementData, &vars, NULL, NULL, NULL,
		    &elemXpts, NULL, &elemWeights, &elemMat);

    for ( int i = 0; i < numElements; i++ ){
      // Retrieve the element variables and node locations
      getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars, vars);

      // Get the element matrix
      elements[i]->getMatType(matType, elemMat, elemXpts, vars);

      // Add the values into the element
      addMatValues(A, i, elemMat, elementIData, elemWeights, matOr);
    }
  }

  A->beginAssembly();
  A->endAssembly();

  A->applyBCs();
}

/*
  Initialize a list of functions
  
  Every function must be initialized - usually just once - before it
  can be evaluated. This is handled automatically within
  TACSAssembler.

  Check whether the functions are associated with this TACSAssembler
  object.  Next, call preInitalize() for each function in the list.
  Go through the function domain and call initialize for each element
  in the domain. Finally, call post initialize for each function in
  the list.

  functions:  an array of function values
  numFuncs:   the number of functions
*/
void TACSAssembler::initializeFunctions( TACSFunction **functions, 
					 int numFuncs ){
  // First check if this is the right assembly object
  for ( int k = 0; k < numFuncs; k++ ){
    if (this != functions[k]->getTACS()){
      fprintf(stderr, "[%d] Cannot evaluate function %s, wrong \
TACSAssembler object\n", mpiRank, functions[k]->functionName());
      return;
    }
  }

  // Test which functions have been initialized
  int count = 0;
  for ( int k = 0; k < numFuncs; k++ ){
    if (!functions[k]->isInitialized()){
      count++;
    }
  }

  if (count == 0){
    return;
  }

  // Create an index-list of the functions that haven't been
  // initialized yet
  int * list = new int[ count ];
  int j = 0;
  for ( int k = 0; k < numFuncs; k++ ){
    if (!functions[k]->isInitialized()){
      list[j] = k;
      j++;
    }
  }

  // Now initialize the functions
  for ( int k = 0; k < count; k++ ){
    functions[list[k]]->preInitialize();
  }
    
  for ( int k = 0; k < count; k++ ){
    if (functions[list[k]]->getDomain() == TACSFunction::ENTIRE_DOMAIN){
      for ( int i = 0; i < numElements; i++ ){
	functions[list[k]]->elementWiseInitialize(elements[i], i);
      }
    }
    else if (functions[list[k]]->getDomain() == TACSFunction::SUB_DOMAIN){
      const int * elementNums;
      int subDomainSize = functions[list[k]]->getElements(&elementNums);
    
      for ( int i = 0; i < subDomainSize; i++ ){
	int elemNum = elementNums[i];
        if (elemNum >= 0 && elemNum < numElements){
	  functions[list[k]]->elementWiseInitialize(elements[elemNum], elemNum);
	}
      }      
    }
  }
  
  for ( int k = 0; k < count; k++ ){
    functions[list[k]]->postInitialize();
  }

  delete [] list;
}

/*
  Evaluate a list of TACS functions

  First, check if the functions are initialized. Obtain the number of
  iterations over the function domain required to evaluate the
  functions.

  This function will print an error and return 0 if the underlying
  TACSAssembler object does not correspond to the TACSAssembler object.

  input:
  functions:  array of functions to evaluate
  numFuncs:   the number of functions to evaluate

  output:
  funcVals: the values of the functions 
*/
void TACSAssembler::evalFunctions( TACSFunction **functions, 
                                   int numFuncs, 
				   TacsScalar *funcVals ){
  // Find the max. number of iterations required
  int num_iters = 0;
  for ( int k = 0; k < numFuncs; k++ ){
    int iters = functions[k]->getNumIterations();
    if (iters > num_iters){
      num_iters = iters;
    }
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemXpts;
  getDataPointers(elementData, &elemVars, NULL, NULL, NULL,
		  &elemXpts, NULL, NULL, NULL);

  // check if initialization is neccessary
  initializeFunctions(functions, numFuncs);
     
  if (thread_info->getNumThreads() > 1){
    // Initialize the threads
    tacsPInfo->tacs = this;
    tacsPInfo->functions = functions;
    tacsPInfo->numFuncs = numFuncs;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for ( int iter = 0; iter < num_iters; iter++ ){
      // Initialize the pre-evaluation iterations
      for ( int k = 0; k < numFuncs; k++ ){
        functions[k]->preEval(iter);
      }

      tacsPInfo->funcIteration = iter;
      numCompletedElements = 0;

      for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
        pthread_create(&threads[k], &attr, 
                       TACSAssembler::evalFunctions_thread,
                       (void*)tacsPInfo);
      }
      
      // Join all the threads
      for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
        pthread_join(threads[k], NULL);
      }
    
      // Initialize the pre-evaluation iterations
      for ( int k = 0; k < numFuncs; k++ ){
        functions[k]->postEval(iter);
      }
    }
  }
  else {
    for ( int iter = 0; iter < num_iters; iter++ ){
      // Initialize the pre-evaluation iterations
      for ( int k = 0; k < numFuncs; k++ ){
        functions[k]->preEval(iter);
      }
      
      // Default work arrays
      TacsScalar work_default[128];
      int iwork_default[128];

      for ( int k = 0; k < numFuncs; k++ ){
        // The work array pointers
        TacsScalar * work = NULL;
        int * iwork = NULL;
        
        // Get the size of the arrays
        int iwork_size = 0, work_size = 0;
        functions[k]->getEvalWorkSizes(&iwork_size, &work_size);
        
        // Check that the work arrays are sufficiently large
        if (work_size > 128){ work = new TacsScalar[work_size]; }
        else { work = work_default; }
        
        if (iwork_size > 128){ iwork = new int[iwork_size]; }
        else { iwork = iwork_default; }
      
        // Initialize the work arrays
        functions[k]->preEvalThread(iter, iwork, work);

        if (functions[k]->getDomain() == TACSFunction::ENTIRE_DOMAIN){
          for ( int i = 0; i < numElements; i++ ){
            // Determine the values of the state 
            // variables for the current element
	    getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
            getValues(varsPerNode, i, localVars, elemVars);
            
            // Evaluate the element-wise component of the function
            functions[k]->elementWiseEval(iter, elements[i], i, 
					  elemXpts, elemVars, iwork, work);
          }
        }
        else if (functions[k]->getDomain() == TACSFunction::SUB_DOMAIN){
          const int * elementNums;
          int subDomainSize = functions[k]->getElements(&elementNums);
          
          for ( int i = 0; i < subDomainSize; i++ ){
            int elemNum = elementNums[i];
            
            if (elemNum >= 0 && elemNum < numElements){
              // Determine the values of the state variables 
              // for the current element
	      getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
	      getValues(varsPerNode, elemNum, localVars, elemVars);
              
              // Evaluate the element-wise component of the function
              functions[k]->elementWiseEval(iter, elements[elemNum], elemNum,
					    elemXpts, elemVars, iwork, work);
            }
          }
        }
      
        functions[k]->postEvalThread(iter, iwork, work);

        // Check that the work arrays are sufficiently large
        if (work_size > 128){ delete [] work; }
        if (iwork_size > 128){ delete [] iwork; }
      }
        
      // Initialize the pre-evaluation iterations
      for ( int k = 0; k < numFuncs; k++ ){
        functions[k]->postEval(iter);
      }
    }
  }

  for ( int k = 0; k < numFuncs; k++ ){
    funcVals[k] = functions[k]->getValue();
  }

  return;
}

/*
  Evaluate the derivative of a list of functions w.r.t. the design
  variables.

  Note that a function should be evaluated - using evalFunction - before
  its derivatives can be evaluated.

  The design variable sensitivities are divided into two distinct sets:
  material-dependent design variables and shape design variables. The
  shape design variables are handled through the TACSNodeMap class. The
  material-dependent design variables are handled through the element 
  classes themselves.

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

  input:
  funcs:     the TACSFunction function objects
  numFuncs:  the number of functions - size of funcs array
  fdvSens:   the sensitivity - size numFuncs*numDVs
  numDVs:    the number of design variables
*/
void TACSAssembler::evalDVSens( TACSFunction **funcs, int numFuncs, 
                                TacsScalar *fdvSens, int numDVs ){
  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemXpts, *elemXptSens;
  getDataPointers(elementData, &elemVars, NULL, NULL, NULL,
		  &elemXpts, NULL, NULL, &elemXptSens);

  // Zero the derivative
  memset(fdvSens, 0, numDVs*numFuncs*sizeof(TacsScalar));

  // check if initialization is neccessary
  initializeFunctions(funcs, numFuncs);

  if (thread_info->getNumThreads() > 1){
    tacsPInfo->tacs = this;
    tacsPInfo->functions = funcs;
    tacsPInfo->numFuncs = numFuncs;
    tacsPInfo->numDesignVars = numDVs;
    tacsPInfo->fdvSens = fdvSens;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Now compute df/dx
    numCompletedElements = 0;
    
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::evalDVSens_thread,
                     (void*)tacsPInfo);
    }
    
    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  }
  else {
    // For each function, evaluate the derivative w.r.t. the 
    // design variables for each element
    for ( int k = 0; k < numFuncs; k++ ){
      TACSFunction * function = funcs[k];

      // Default work arrays
      TacsScalar * work = NULL;
      TacsScalar work_default[128];
      int work_size = function->getDVSensWorkSize();
      if (work_size > 128){ work = new TacsScalar[ work_size ]; }
      else { work = work_default; }

      if (function->getDomain() == TACSFunction::SUB_DOMAIN){
        // Get the function sub-domain
        const int * elemSubList;
        int numSubElems = function->getElements(&elemSubList);
      
        for ( int i = 0; i < numSubElems; i++ ){
          int elemNum = elemSubList[i];
          // Determine the values of the state variables for subElem
          getValues(varsPerNode, elemNum, localVars, elemVars);
          getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
	
          // Evaluate the element-wise sensitivity of the function
          function->elementWiseDVSens(&fdvSens[k*numDVs], numDVs,
                                      elements[elemNum], elemNum, 
                                      elemXpts, elemVars, work);	      
        }
      }
      else if (function->getDomain() == TACSFunction::ENTIRE_DOMAIN){
        for ( int elemNum = 0; elemNum < numElements; elemNum++ ){
          // Determine the values of the state variables for elemNum
          getValues(varsPerNode, elemNum, localVars, elemVars);
          getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
        
          // Evaluate the element-wise sensitivity of the function
          function->elementWiseDVSens(&fdvSens[k*numDVs], numDVs,
                                      elements[elemNum], elemNum, 
                                      elemXpts, elemVars, work);
        }
      }

      if (work_size > 128){ delete [] work; }
    }
  }

  // Reduce the gradient information on all processors
  MPI_Allreduce(MPI_IN_PLACE, fdvSens, numFuncs*numDVs, TACS_MPI_TYPE,  
		MPI_SUM, tacs_comm);
}

/*
  Evaluate the derivative of the function w.r.t. the owned nodes.

  This code evaluates the sensitivity of the function w.r.t. the 
  owned nodes for all elements in the function domain. 

  Note that a function should be evaluated - using evalFunction - before
  its derivatives can be evaluated.

  This function should be preferred to the use of evalDVSens without a 
  list of functions since it is more efficient!

  input:
  funcs:     the TACSFunction function objects
  numFuncs:  the number of functions - size of funcs array
  fXptSens:  the sensitivity - size numFuncs*numNodes*3
*/
/*
void TACSAssembler::evalXptSens( TACSFunction **funcs, int numFuncs, 
				 TacsScalar *fXptSens ){
  // First check if this is the right assembly object
  for ( int k = 0; k < numFuncs; k++ ){
    if (this != funcs[k]->getTACS()){
      fprintf(stderr, "[%d] Cannot evaluate function %s, wrong TACS object\n", 
              mpiRank, funcs[k]->functionName());
    }
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemXpts, *elemXptSens;
  getDataPointers(elementData, &elemVars, NULL, NULL, NULL,
		  &elemXpts, NULL, NULL, &elemXptSens);

  memset(fXptSens, 0, TACS_SPATIAL_DIM*numNodes*numFuncs*sizeof(TacsScalar));

  // check if initialization is neccessary
  initializeFunctions(funcs, numFuncs);
  if (thread_info->getNumThreads() > 1){
    tacsPInfo->tacs = this;
    tacsPInfo->functions = funcs;
    tacsPInfo->numFuncs = numFuncs;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Compute the derivative of the functions w.r.t. the nodes
    numCompletedElements = 0;
    tacsPInfo->fXptSens = fXptSens;

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::evalXptSens_thread,
		     (void*)tacsPInfo);
    }

    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }
  }
  else {
    // For each function, evaluate the derivative w.r.t. the 
    // nodal locations for all elements or part of the domain
    for ( int k = 0; k < numFuncs; k++ ){
      TACSFunction * function = funcs[k];
    
      // Default work arrays
      TacsScalar * work = NULL;
      TacsScalar work_default[128];
      int work_size = function->getXptSensWorkSize();
      if (work_size > 128){ work = new TacsScalar[ work_size ]; }
      else { work = work_default; }

      if (function->getDomain() == TACSFunction::SUB_DOMAIN){
	// Get the function sub-domain
	const int * elemSubList;
	int numSubElems = function->getElements(&elemSubList);
	for ( int i = 0; i < numSubElems; i++ ){
	  int elemNum = elemSubList[i];
	  // Determine the values of the state variables for subElem
	  getValues(varsPerNode, elemNum, localVars[loadCase], elemVars);
	  getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
          
	  // Evaluate the element-wise sensitivity of the function
	  function->elementWiseXptSens(elemXptSens, 
				       elements[elemNum], elemNum, 
				       elemXpts, elemVars, work);
          
	  addValues(TACS_SPATIAL_DIM, elemNum, elemXptSens,
		    &fXptSens[TACS_SPATIAL_DIM*k*numNodes]);
	}
      }
      else if (function->getDomain() == TACSFunction::ENTIRE_DOMAIN){
	for ( int elemNum = 0; elemNum < numElements; elemNum++ ){
	  // Determine the values of the state variables for elemNum
	  getValues(varsPerNode, elemNum, localVars[loadCase], elemVars);
	  getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
      
	  // Evaluate the element-wise sensitivity of the function
	  function->elementWiseXptSens(elemXptSens, 
				       elements[elemNum], elemNum, 
				       elemXpts, elemVars, work);
      
	  addValues(TACS_SPATIAL_DIM, elemNum, elemXptSens,
		    &fXptSens[TACS_SPATIAL_DIM*k*numNodes]);
	}
      }
  
      if (work_size > 128){ delete [] work; }
    }
  }
}
*/

/*
  Evaluate the derivative of the function w.r.t. the state variables.

  This code evaluates the sensitivity of the function w.r.t. the 
  state variables for all elements in the function domain. This code
  is usually much faster than the code for computing the derivative of 
  the function w.r.t. the design variables. 

  Note that the sensitivity vector 'vec' is assembled, and appropriate
  boundary conditions are imposed before the function is returned.

  function: the function pointer
  vec:      the derivative of the function w.r.t. the state variables
*/
void TACSAssembler::evalSVSens( TACSFunction *function, BVec *vec ){
  // First check if this is the right assembly object
  if (this != function->getTACS()){
    fprintf(stderr, "[%d] Cannot evaluate function, wrong TACS object\n",
	    mpiRank);
    return;
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemRes, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemRes, NULL, NULL,
		  &elemXpts, NULL, NULL, NULL);

  // Zero the vector
  vec->zeroEntries();

  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localRes, 0, size*sizeof(TacsScalar));

  // Perform the initialization if neccessary
  initializeFunctions(&function, 1);
  
  // Default work arrays
  TacsScalar * work = NULL;
  TacsScalar work_default[128];
  int work_size = function->getSVSensWorkSize();
  if (work_size > 128){ work = new TacsScalar[ work_size ]; }
  else { work = work_default; }

  if (function->getDomain() == TACSFunction::ENTIRE_DOMAIN){
    for ( int i = 0; i < numElements; i++ ){
      // Determine the values of the state variables for the current element
      getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars, elemVars);
      
      // Evaluate the element-wise sensitivity of the function
      function->elementWiseSVSens(elemRes, elements[i], i, 
				  elemXpts, elemVars, work);

      addValues(varsPerNode, i, elemRes, localRes);
    }
  }
  else if (function->getDomain() == TACSFunction::SUB_DOMAIN){
    const int * elementNums;
    int subDomainSize = function->getElements(&elementNums);
    
    for ( int i = 0; i < subDomainSize; i++ ){
      int elemNum = elementNums[i];

      if (elemNum >= 0 && elemNum < numElements){	
      	// Determine the values of the state variables for the current element
	getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);
	getValues(varsPerNode, elemNum, localVars, elemVars);

      	function->elementWiseSVSens(elemRes, elements[elemNum], elemNum, 
				    elemXpts, elemVars, work);

	addValues(varsPerNode, elemNum, elemRes, localRes);
      }
    }
  }

  if (work_size > 128){ delete [] work; }

  // Add the contribution to the derivative from the dependent-nodes
  addDependentResidual(varsPerNode, localRes);

  // Add up the derivative of the function w.r.t. the state variables
  // across all processors
  vecDist->beginReverse(localRes, vec, BVecDistribute::ADD);
  vecDist->endReverse(localRes, vec, BVecDistribute::ADD);
  vec->applyBCs();
}

/*
  Evaluate the product of several ajdoint vectors with the derivative
  of the residual w.r.t. the design variables.

  This function is collective on all TACSAssembler processes. This
  computes the product of the derivative of the residual w.r.t. the
  design variables with several adjoint vectors simultaneously. This
  saves computational time as the derivative of the element residuals
  can be reused for each adjoint vector. This function performs the
  same task as evalAdjointResProduct, but uses more memory than
  calling it for each adjoint vector.

  adjoint:     the array of adjoint vectors
  numAdjoints: the number of adjoint vectors
  dvSens:      the product of the derivative of the residuals and the adjoint
  numDVs:      the number of design variables
*/
void TACSAssembler::evalAdjointResProducts( BVec ** adjoint, int numAdjoints, 
                                            TacsScalar * dvSens, int numDVs ){
  // Allocate memory for the local components of the adjoint vector
  int nvars = varsPerNode*(numNodes + numDependentNodes);
  TacsScalar * localAdjoint = new TacsScalar[ nvars*numAdjoints ];

  // Perform the transfer of the adjoint variables 
  for ( int k = 0; k < numAdjoints; k++ ){
    // Scatter the variables to the local processes
    // now localRes == the adjoint variables
    vecDist->beginForward(adjoint[k], &localAdjoint[k*nvars]);
    vecDist->endForward(adjoint[k], &localAdjoint[k*nvars]);
    setDependentVariables(varsPerNode, &localAdjoint[k*nvars]);
  }

  // Allocate space for the design derivative
  TacsScalar * dvSensVals = new TacsScalar[ numDVs*numAdjoints ];
  memset(dvSens, 0, numDVs*numAdjoints*sizeof(TacsScalar));    
  memset(dvSensVals, 0, numDVs*numAdjoints*sizeof(TacsScalar));

  // Sort the list of auxiliary elements - this only performs the
  // sort if it is required (if new elements are added)
  if (aux_elements){
    aux_elements->sort();
  }

  if (thread_info->getNumThreads() > 1){
    tacsPInfo->tacs = this;
    tacsPInfo->numAdjoints = numAdjoints;
    tacsPInfo->adjointVars = localAdjoint;
    tacsPInfo->numDesignVars = numDVs;
    tacsPInfo->fdvSens = dvSensVals;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Now compute Phi^{T} * dR/dx 
    numCompletedElements = 0;

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
                     TACSAssembler::adjointResProduct_thread,
                     (void*)tacsPInfo);
    }
    
    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }

    // Destroy the attribute
    pthread_attr_destroy(&attr);
  }
  else {
    // Retrieve pointers to temporary storage
    TacsScalar *vars, *dvars, *ddvars;
    TacsScalar *elemXpts, *elemAdjoint;
    getDataPointers(elementData, &vars, &dvars, &ddvars, &elemAdjoint,
		    &elemXpts, NULL, NULL, NULL);

    // Set the data for the auxiliary elements - if there are any
    int naux = 0, aux_count = 0;
    TACSAuxElem *aux = NULL;
    if (aux_elements){
      naux = aux_elements->getAuxElements(&aux);
    }

    for ( int i = 0; i < numElements; i++ ){
      // Find the variables and nodes
      int nnodes = getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars, vars);
      getValues(varsPerNode, i, localDotVars, dvars);
      getValues(varsPerNode, i, localDDotVars, ddvars);
      int nevars = varsPerNode*nnodes;
      
      // Get the adjoint variables
      for ( int k = 0; k < numAdjoints; k++ ){
	double scale = 1.0;
	getValues(varsPerNode, i, &localAdjoint[nvars*k], elemAdjoint);
	
	elements[i]->addAdjResProduct(time, scale, 
                                      &dvSensVals[k*numDVs], numDVs,
				      elemAdjoint, elemXpts,
				      vars, dvars, ddvars);

        // Add the contribution from the auxiliary elements
        while (aux_count < naux && aux[aux_count].num == i){
          aux[aux_count].elem->addAdjResProduct(time, scale,
                                                &dvSensVals[k*numDVs], numDVs,
                                                elemAdjoint, elemXpts,
                                                vars, dvars, ddvars);
          aux_count++;
        }
      }
    }
  }

  // Free the local adjoint array
  delete [] localAdjoint;

  // Collect the products from all processes
  // Component wise summation
  MPI_Allreduce(dvSensVals, dvSens, numDVs*numAdjoints, TACS_MPI_TYPE, 
		MPI_SUM, tacs_comm);

  // Free the allocated design vars
  delete [] dvSensVals;
}

/*
  Evaluate the derivative of an inner product of two vectors with a
  matrix of a given type. This code does not explicitly evaluate the
  element matrices.  Instead, the inner product contribution from each
  element matrix is added to the final result. This implementation
  saves considerable computational time and memory.

  input:
  matType:   the matrix type
  psi:       the left-multiplying vector
  phi:       the right-multiplying vector
  numDVs:    the length of the design variable array

  output:
  dvSens:    the derivative of the inner product 
*/
void TACSAssembler::evalMatDVSensInnerProduct( TacsScalar scale, 
					       ElementMatrixType matType, 
					       BVec *psi, BVec *phi,
					       TacsScalar *dvSens, int numDVs ){
  // Zero the derivative of the design variables
  memset(dvSens, 0, numDVs*sizeof(TacsScalar));

  // Allocate memory for the local components of the adjoint vector
  int nvars = varsPerNode*(numNodes + numDependentNodes);
  TacsScalar * localVals = new TacsScalar[ 2*nvars ];

  // Transfer the varaibles to the local arrays
  vecDist->beginForward(psi, &localVals[0]);
  vecDist->endForward(psi, &localVals[0]);

  vecDist->beginForward(phi, &localVals[nvars]);
  vecDist->endForward(phi, &localVals[nvars]);

  setDependentVariables(varsPerNode, &localVals[0]);
  setDependentVariables(varsPerNode, &localVals[nvars]);

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemPsi, *elemPhi, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemPsi, &elemPhi, NULL,
		  &elemXpts, NULL, NULL, NULL);

  // Go through each element in the domain and compute the derivative
  // of the residuals with respect to each design variable and multiply by
  // the adjoint variables
  for ( int i = 0; i < numElements; i++ ){
    // Find the variables and nodes
    getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
    getValues(varsPerNode, i, localVars, elemVars);
    getValues(varsPerNode, i, &localVals[0], elemPsi);
    getValues(varsPerNode, i, &localVals[nvars], elemPhi);
    
    // Add the contribution to the design variable vector
    elements[i]->addMatDVSensInnerProduct(matType, scale, dvSens, numDVs,
					  elemPsi, elemPhi, elemXpts,
					  elemVars);
  }
  
  // Collect the products from all processes
  // Component wise summation
  MPI_Allreduce(MPI_IN_PLACE, dvSens, numDVs, TACS_MPI_TYPE, 
		MPI_SUM, tacs_comm);

  delete [] localVals;
}

/*
  Evaluate the derivative of the inner product of two vectors with a
  matrix with respect to the state variables. This is only defined for
  nonlinear matrices, like the geometric stiffness matrix.  Instead of
  computing the derivative of the matrix for each vector component and
  then computing the inner product, this code computes the derivative
  of the inner product directly, saving computational time and memory.

  input:
  scale:     the scaling parameter applied to the derivative
  matType:   the matrix type
  psi:       the left-multiplying vector
  phi:       the right-multiplying vector
  numDVs:    the length of the design variable array

  output:
  res:       the derivative of the inner product w.r.t. the state vars
*/
void TACSAssembler::evalMatSVSensInnerProduct( TacsScalar scale, 
					       ElementMatrixType matType, 
					       BVec *psi, BVec *phi, BVec *res ){
  // Zero the entries in the residual vector
  res->zeroEntries();

  // Zero the entries in the local copy of the residual vector
  int size = varsPerNode*(numNodes + numDependentNodes);
  memset(localRes, 0, size*sizeof(TacsScalar));

  // Allocate memory for the local components of the adjoint vector
  int nvars = varsPerNode*(numNodes + numDependentNodes);
  TacsScalar * localVals = new TacsScalar[ 2*nvars ];

  // Transfer the varaibles to the local arrays
  vecDist->beginForward(psi, &localVals[0]);
  vecDist->endForward(psi, &localVals[0]);

  vecDist->beginForward(phi, &localVals[nvars]);
  vecDist->endForward(phi, &localVals[nvars]);

  setDependentVariables(varsPerNode, &localVals[0]);
  setDependentVariables(varsPerNode, &localVals[nvars]);

  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemPsi, *elemPhi, *elemRes, *elemXpts;
  getDataPointers(elementData, &elemVars, &elemPsi, &elemPhi, &elemRes,
		  &elemXpts, NULL, NULL, NULL);

  // Go through each element in the domain and compute the derivative
  // of the residuals with respect to each design variable and multiply by
  // the adjoint variables
  for ( int i = 0; i < numElements; i++ ){
    // Find the variables and nodes
    int nnodes = getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
    getValues(varsPerNode, i, localVars, elemVars);
    getValues(varsPerNode, i, &localVals[0], elemPsi);
    getValues(varsPerNode, i, &localVals[nvars], elemPhi);
    int nevars = varsPerNode*nnodes;
    memset(elemRes, 0, nevars*sizeof(TacsScalar));

    // Add the contribution to the design variable vector
    elements[i]->getMatSVSensInnerProduct(matType, elemRes,
					  elemPsi, elemPhi, elemXpts,
					  elemVars);

    // Add the residual values to the local residual array
    addValues(varsPerNode, i, elemRes, localRes);
  }

  // Add the dependent-variable residual from the dependent nodes
  addDependentResidual(varsPerNode, localRes);

  // Assemble the full residual
  vecDist->beginReverse(localRes, res, BVecDistribute::ADD);
  vecDist->endReverse(localRes, res, BVecDistribute::ADD);

  // Apply the boundary conditions to the fully assembled vector
  res->applyBCs();
  
  // Free the allocated memory
  delete [] localVals;
}

/*
  Evaluate the matrix-free Jacobian-vector product of the input vector
  x and store the result in the output vector y.

  This code does not assemble a matrix, but does compute the
  element-wise matricies. This code is not a finite-difference
  matrix-vector product implementation.

  Since the element Jacobian matrices are computed exactly, we can
  evaluate either a regular matrix-product or the transpose matrix
  product.

  input:
  scale:     the scalar coefficient
  alpha:     coefficient on the variables
  beta:      coefficient on the time-derivative terms
  gamma:     coefficient on the second time derivative term
  x:         the input vector
  matOr:     the matrix orientation
  
  output:
  y:         the output vector y <- y + scale*J^{Op}*x
*/
void TACSAssembler::addJacobianVecProduct( TacsScalar scale, 
                                           double alpha, double beta, double gamma,
                                           BVec *x, BVec *y,
                                           MatrixOrientation matOr ){
  // Add the result to the residual
  int size = varsPerNode*(numNodes + numDependentNodes);
  TacsScalar *xlocal = new TacsScalar[ size ];
  memset(localRes, 0, size*sizeof(TacsScalar));

  // Send the variables forward
  vecDist->beginForward(x, xlocal);
  vecDist->endForward(x, xlocal);
  setDependentVariables(varsPerNode, xlocal);

  // Retrieve pointers to temporary storage
  TacsScalar *vars, *dvars, *ddvars, *yvars, *elemXpts;
  TacsScalar *elemWeights, *elemMat;
  getDataPointers(elementData, 
                  &vars, &dvars, &ddvars, &yvars,
                  &elemXpts, NULL, &elemWeights, &elemMat);
  
  // Set the data for the auxiliary elements - if there are any
  int naux = 0, aux_count = 0;
  TACSAuxElem *aux = NULL;
  if (aux_elements){
    naux = aux_elements->getAuxElements(&aux);
  }

  for ( int i = 0; i < numElements; i++ ){
    getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
    getValues(varsPerNode, i, localVars, vars);
    getValues(varsPerNode, i, localDotVars, dvars);
    getValues(varsPerNode, i, localDDotVars, ddvars);
    
    // Get the number of variables from the element
    int nvars = elements[i]->numVariables();
    
    // Compute and add the contributions to the Jacobian
    memset(elemMat, 0, nvars*nvars*sizeof(TacsScalar));
    elements[i]->addJacobian(time, elemMat, alpha, beta, gamma,
                             elemXpts, vars, dvars, ddvars);

    // Add the contribution to the residual and the Jacobian
    // from the auxiliary elements - if any
    while (aux_count < naux && aux[aux_count].num == i){
      aux[aux_count].elem->addJacobian(time, elemMat, 
                                       alpha, beta, gamma,
                                       elemXpts, vars, dvars, ddvars);      
      aux_count++;
    }

    // Temporarily set the variable array as the element input array
    // and get the local variable input values from the local array.
    TacsScalar *xvars = vars;
    getValues(varsPerNode, i, xlocal, xvars);
   
    // Take the matrix vector product. Note the matrix is stored in
    // row-major order and BLAS assumes column-major order. As a
    // result, the transpose arguments are reversed.
    TacsScalar zero = 0.0;
    int incx = 1;
    if (matOr == NORMAL){
      BLASgemv("T", &nvars, &nvars, &scale, elemMat, &nvars, 
               xvars, &incx, &zero, yvars, &incx);
    }
    else {
      BLASgemv("N", &nvars, &nvars, &scale, elemMat, &nvars, 
               xvars, &incx, &zero, yvars, &incx);
    }

    // Add the residual values
    addValues(varsPerNode, i, yvars, localRes);
  }

  // Add the dependent-variable residual from the dependent nodes
  addDependentResidual(varsPerNode, localRes);

  // Assemble the full residual
  vecDist->beginReverse(localRes, y, BVecDistribute::ADD);
  vecDist->endReverse(localRes, y, BVecDistribute::ADD);

  // Set the boundary conditions
  y->applyBCs();

  delete [] xlocal;
}

/*
  Evaluate the product of several ajdoint vectors with the derivative
  of the residual w.r.t. the nodal points.

  This function is collective on all TACSAssembler processes. This
  computes the product of the derivative of the residual w.r.t. the
  nodal points with several adjoint vectors simultaneously. This
  saves computational time as the derivative of the element residuals
  can be reused for each adjoint vector. 

  adjoint:     the array of adjoint vectors
  numAdjoints: the number of adjoint vectors
  dvSens:      the product of the derivative of the residuals and the adjoint
  numDVs:      the number of design variables
*/
/*
void TACSAssembler::evalAdjointResXptSensProducts( int loadCase, 
						   BVec ** adjoint, 
						   int numAdjoints, 
						   TacsScalar * adjXptSensProduct){
  if (loadCase < 0 || loadCase >= numLoadCases){
    fprintf(stderr, "[%d] Load case number %d out of range [0,%d)\n", 
	    mpiRank, loadCase, numLoadCases);
    return;
  }

  // Allocate memory for the local components of the adjoint vector
  int nvars = varsPerNode*(numNodes + numDependentNodes);
  TacsScalar * localAdjoint = new TacsScalar[ nvars*numAdjoints ];

  // Perform the transfer of the adjoint variables 
  for ( int k = 0; k < numAdjoints; k++ ){
    // Scatter the variables to the local processes
    // now localRes == the adjoint variables
    vecDist->beginForward(adjoint[k], &localAdjoint[k*nvars]);
    vecDist->endForward(adjoint[k], &localAdjoint[k*nvars]);
    setDependentVariables(varsPerNode, &localAdjoint[k*nvars]);
  }

  if (thread_info->getNumThreads() > 1){
    tacsPInfo->tacs = this;
    tacsPInfo->numAdjoints = numAdjoints;
    tacsPInfo->adjointVars = localAdjoint;

    // Create the joinable attribute
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // Compute adjResProduct = Phi^{T} dR/dXpts

    // Set the vector to store Phi^{T} dR/dXpts
    numCompletedElements = 0;
    tacsPInfo->adjXptSensProduct = adjXptSensProduct;

    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_create(&threads[k], &attr, 
		     TACSAssembler::adjointResXptSensProduct_thread,
		     (void*)tacsPInfo);
    }
    
    // Join all the threads
    for ( int k = 0; k < thread_info->getNumThreads(); k++ ){
      pthread_join(threads[k], NULL);
    }
  }
  else {
    // Retrieve pointers to temporary storage
    TacsScalar *elemVars, *elemResDVSens, *elemXpts, *elemResXptSens;
    getDataPointers(elementData, &elemVars, &elemResDVSens, NULL, NULL,
		    &elemXpts, NULL, NULL, &elemResXptSens);

    // Allocate memory for the element adjoint variables and 
    // elemXptSens = the product of the element adjoint variables and
    // the derivative of the residuals w.r.t. the nodes
    int s = maxElementSize;
    int sx = TACS_SPATIAL_DIM*maxElementNodes;
    TacsScalar * elemAdjoint = new TacsScalar[ s*numAdjoints ];
    TacsScalar * elemXptSens = new TacsScalar[ sx*numAdjoints ];
    
    // Get the surface traction information
    int numStElems = 0;
    const int * stElemNums = NULL;
    
    if (surfaceTractions[loadCase]){
      numStElems = surfaceTractions[loadCase]->getElementNums(&stElemNums);
    }

    TacsScalar lambda = loadFactor[loadCase];

    int stIndex = 0;      
    for ( int i = 0; i < numElements; i++ ){
      // Get the variables and nodes for this element
      int nnodes = getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);
      getValues(varsPerNode, i, localVars[loadCase], elemVars);
      int nevars = varsPerNode*nnodes;

      // Get the adjoint variables associated with this element
      for ( int k = 0; k < numAdjoints; k++ ){
	getValues(varsPerNode, i, &localAdjoint[k*nvars],
		  &elemAdjoint[k*nevars]);
      }

      // Compute the derivative of the element residuals w.r.t. the
      // element nodes
      elements[i]->getResXptSens(elemResXptSens, 
				 elemVars, elemXpts);               

      // Add the contributions from the derivative of the consistent forces
      // w.r.t. the element nodes
      while ((stIndex < numStElems) && (stElemNums[stIndex] == i)){
	TACSElementTraction * elemTraction =
	  surfaceTractions[loadCase]->getElement(stIndex);
	elemTraction->addForceXptSens(lambda, elemResXptSens,
				      elemVars, elemXpts);
	stIndex++;
      }

      // Compute the product of the derivative of the residual w.r.t. the
      // nodal coordinates and the adjoint variables
      // Need to compute: 
      // elemXptSens = elementResXptSens^{T} * elemAdjoint
      int nenodes = TACS_SPATIAL_DIM*nnodes;
      TacsScalar alpha = 1.0, beta = 0.0;
      BLASgemm("T", "N", &nenodes, &numAdjoints, &nevars,
	       &alpha, elemResXptSens, &nevars,
	       elemAdjoint, &nevars,                 
	       &beta, elemXptSens, &nenodes);
      
      // Add the values in elemXptSens into adjXptSensProduct
      for ( int k = 0; k < numAdjoints; k++ ){
	addValues(TACS_SPATIAL_DIM, i, &elemXptSens[k*nenodes],
		  &adjXptSensProduct[TACS_SPATIAL_DIM*numNodes*k]);
      }
    }
    
    delete [] elemAdjoint;
    delete [] elemXptSens;
  }
  delete [] localAdjoint;
}
*/
/*!  
  Given the load case number and the element, return the element
  object and the values of the variables and nodal locations
  associated with that element.

  This is useful for post-processing data without having to write a
  full function class.

  loadCase:    the load case number
  elemNum:     the element number
  elemVars:    the variables associated with elemNum
  elemXpts:    the nodal locations associated with elemNum

  returns:     the element pointer associated with elemNum
*/
TACSElement * TACSAssembler::getElement( int elemNum,
					 TacsScalar *elemXpts,
					 TacsScalar *vars,
					 TacsScalar *dvars,
					 TacsScalar *ddvars ){
  // Check if the element number is in the right range
  if (elemNum < 0 || elemNum >= numElements){
    fprintf(stderr, "[%d] Element number %d out of range [0,%d)\n",
	    mpiRank, elemNum, numElements);
    return NULL;
  }

  if (elemXpts){ getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts); }
  if (vars){ getValues(varsPerNode, elemNum, localVars, vars); }
  if (dvars){ getValues(varsPerNode, elemNum, localVars, dvars); }
  if (ddvars){ getValues(varsPerNode, elemNum, localVars, ddvars); }

  return elements[elemNum];
}

/*
  Test the implementation of the given element number.

  This tests the stiffness matrix and various parts of the
  design-sensitivities: the derivative of the determinant of the
  Jacobian, the derivative of the strain w.r.t. the nodal coordinates,
  and the state variables and the derivative of the residual w.r.t.
  the design variables and nodal coordiantes.

  elemNum:     the element number to test
  print_level: the print level to use   
*/
void TACSAssembler::testElement( int elemNum, int print_level ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call testElement() before finalize()\n", 
	    mpiRank);
    return;
  }
  else if (elemNum < 0 || elemNum >= numElements){
    fprintf(stderr, "[%d] Element number %d out of range [0,%d)\n",
	    mpiRank, elemNum, numElements);
    return;
  }

  // Retrieve pointers to temporary storage
  TacsScalar *elemXpts;
  getDataPointers(elementData, NULL, NULL, NULL, NULL,
		  &elemXpts, NULL, NULL, NULL);

  // Get the element node locations
  int nnodes = getValues(TACS_SPATIAL_DIM, elemNum, Xpts, elemXpts);

  // Create the element test function
  double pt[] = {0.0, 0.0, 0.0};
  TestElement * test = new TestElement(elements[elemNum],
				       elemXpts);
  test->incref();
  test->setPrintLevel(print_level);
  
  printf("Testing element %s\n", elements[elemNum]->elementName());
  if (test->testJacobian()){ printf("Stiffness matrix failed\n"); }
  else { printf("Stiffness matrix passed\n"); }
  if (test->testJacobianXptSens(pt)){ printf("Jacobian XptSens failed\n"); }
  else { printf("Jacobian XptSens passed\n"); }
  if (test->testStrainSVSens(pt)){ printf("Strain SVSens failed\n"); }
  else { printf("Strain SVSens passed\n"); }
  /*if
 (test->testResXptSens()){ printf("Res XptSens failed\n"); }
  else { printf("Res XptSens passed\n"); }
  if (test->testStrainXptSens(pt)){ printf("Strain XptSens failed\n"); }
  else { printf("Strain XptSens passed\n"); }
  if (test->testResDVSens()){ printf("Res DVSens failed\n"); }
  else { printf("Res DVSens passed\n"); }
  if (test->testMatDVSens(STIFFNESS_MATRIX)){
    printf("Stiffness matrix DVSens failed\n"); }
  else { printf("Stiffness Matrix DVSens passed\n"); }
  if (test->testMatDVSens(MASS_MATRIX)){
    printf("Mass matrix DVSens failed\n"); }
  else { printf("Mass Matrix DVSens passed\n"); }
  if (test->testMatDVSens(GEOMETRIC_STIFFNESS_MATRIX)){ 
    printf("Geometric stiffness matrix DVSens failed\n"); }
  else { printf("Geometric stiffness Matrix DVSens passed\n"); }
  */
  printf("\n");
  
  test->decref();
}

/*
  Test the implementation of the given element's constitutive class.
  
  This function tests the failure computation and the mass
  sensitivities for the given element.

  elemNum:     the element to retrieve the constitutive object from
  print_level: the print level to use for the test
*/
void TACSAssembler::testConstitutive( int elemNum, int print_level ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call testConstitutive() before finalize()\n", 
	    mpiRank);
    return;
  }
  else if (elemNum < 0 || elemNum >= numElements){
    fprintf(stderr, "[%d] Element number %d out of range [0,%d)\n",
	    mpiRank, elemNum, numElements);
    return;
  }

  TACSConstitutive * stiffness = elements[elemNum]->getConstitutive();
  
  if (stiffness){
    double pt[] = {0.0, 0.0, 0.0};
    TestConstitutive * test = new TestConstitutive(stiffness);
    test->incref();
    test->setPrintLevel(print_level);
    
    printf("Testing constitutive class %s\n", stiffness->constitutiveName());
    if (test->testFailStrainSens(pt)){ printf("Fail StrainSens failed\n"); }
    else { printf("Fail StrainSens passed\n"); }
    if (test->testBucklingStrainSens()){ printf("Buckling StrainSens failed\n"); }
    else { printf("Buckling StrainSens passed\n"); }
    /*
    if (test->testMassDVSens(pt)){ printf("Mass DVSens failed\n"); }
    else { printf("Mass DVSens passed\n"); }
    if (test->testFailDVSens(pt)){ printf("Fail DVSens failed\n"); }
    else { printf("Fail DVSens passed\n"); }
    if (test->testBucklingDVSens()){ printf("Buckling DVSens failed\n"); }
    else { printf("Buckling DVSens passed\n"); }
    */
    printf("\n");
 
    test->decref();
  }
}

/*
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

  func:            the function to test
  num_design_vars: the number of design variables to use
  dh:              the step size to use
*/
void TACSAssembler::testFunction( TACSFunction * func, 
				  int num_design_vars,
				  double dh ){
  if (!meshFinalizedFlag){
    fprintf(stderr, "[%d] Cannot call testFunction() before finalize()\n", 
	    mpiRank);
    return;
  }

  // First, test the design variable values
  TacsScalar * x = new TacsScalar[ num_design_vars ];
  TacsScalar * xpert = new TacsScalar[ num_design_vars ];
  TacsScalar * xtemp = new TacsScalar[ num_design_vars ];

  for ( int k = 0; k < num_design_vars; k++ ){
    xpert[k] = (1.0*rand())/RAND_MAX;
  }

  MPI_Bcast(xpert, num_design_vars, TACS_MPI_TYPE, 0, tacs_comm);

  getDesignVars(x, num_design_vars);
  setDesignVars(x, num_design_vars);

  TacsScalar fd = 0.0;
#ifdef TACS_USE_COMPLEX
  // Compute the function at the point x + dh*xpert
  for ( int k = 0; k < num_design_vars; k++ ){
    xtemp[k] = x[k] + TacsScalar(0.0, dh)*xpert[k];
  }
  setDesignVars(xtemp, num_design_vars);
  evalFunctions(&func, 1, &fd);
  fd = ImagPart(fd)/dh;
#else
  // Compute the function at the point x + dh*xpert
  for ( int k = 0; k < num_design_vars; k++ ){
    xtemp[k] = x[k] + dh*xpert[k];
  }
  setDesignVars(xtemp, num_design_vars);
  TacsScalar fval0;
  evalFunctions(&func, 1, &fval0);

  // Compute the function at the point x - dh*xpert
  for ( int k = 0; k < num_design_vars; k++ ){
    xtemp[k] = x[k] - dh*xpert[k];
  }
  setDesignVars(xtemp, num_design_vars);
  TacsScalar fval1;
  evalFunctions(&func, 1, &fval1);
  fd = 0.5*(fval0 - fval1)/dh;
#endif // TACS_USE_COMPLEX

  // Compute df/dx
  TacsScalar ftmp;
  setDesignVars(x, num_design_vars);
  evalFunctions(&func, 1, &ftmp);
  evalDVSens(&func, 1, xtemp, num_design_vars);
  
  // Compute df/dx^{T} * xpert
  TacsScalar pdf = 0.0;
  for ( int k = 0; k < num_design_vars; k++ ){
    pdf += xtemp[k]*xpert[k];
  }

  if (mpiRank == 0){
    fprintf(stderr, "Testing function %s\n", func->functionName());
    const char * descript = "df/dx^{T}p";
    fprintf(stderr, "%*s[   ] %15s %15s %15s\n",
	    (int)strlen(descript), "Val", "Analytic", 
            "Approximate", "Rel. Error");
    if (pdf != 0.0){
      fprintf(stderr, "%s[%3d] %15.6e %15.6e %15.4e\n", 
              descript, 0, RealPart(pdf), RealPart(fd), 
              fabs(RealPart((pdf - fd)/pdf)));
    }
    else {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e\n", 
              descript, 0, RealPart(pdf), RealPart(fd));
    }
  }

  delete [] x;
  delete [] xtemp;
  delete [] xpert;
  
  BVec * temp = createVec();
  BVec * pert = createVec();
  BVec * vars = createVec();
  temp->incref();
  pert->incref();
  vars->incref();

  getVariables(vars);

  // Set up a random perturbation 
  pert->setRand(-1.0, 1.0);
  pert->applyBCs();
  pert->scale(vars->norm()/pert->norm());

#ifdef TACS_USE_COMPLEX
  // Evaluate the function at vars + dh*pert
  temp->copyValues(vars);
  temp->axpy(TacsScalar(0.0, dh), pert);
  setVariables(temp);

  evalFunctions(&func, 1, &fd);
  fd = ImagPart(fd)/dh;
#else
  // Evaluate the function at vars + dh*pert
  temp->copyValues(vars);
  temp->axpy(dh, pert);
  setVariables(temp);
  evalFunctions(&func, 1, &fval0);

  // Evaluate the function at vars - dh*pert
  temp->copyValues(vars);
  temp->axpy(-dh, pert);
  setVariables(temp);
  evalFunctions(&func, 1, &fval1);
  
  fd = 0.5*(fval0 - fval1)/dh;
#endif // TACS_USE_COMPLEX

  // Reset the variable values
  setVariables(vars);

  // Evaluate the state variable sensitivity
  evalFunctions(&func, 1, &ftmp);
  evalSVSens(func, temp);
  pdf = temp->dot(pert);

  if (mpiRank == 0){
    const char * descript = "df/du^{T}p";
    fprintf(stderr, "%*s[   ] %15s %15s %15s\n",
	    (int)strlen(descript), "Val", "Analytic", 
            "Approximate", "Rel. Error");
    if (pdf != 0.0){
      fprintf(stderr, "%s[%3d] %15.6e %15.6e %15.4e\n", 
              descript, 0, RealPart(pdf), RealPart(fd), 
              fabs(RealPart((pdf - fd)/pdf)));
    }
    else {
      fprintf(stderr, "%s[%3d] %15.6e %15.6e\n", 
              descript, 0, RealPart(pdf), RealPart(fd));
    }
  }

  temp->decref();
  pert->decref();
  vars->decref();  
}

/*!  
  Determine the number of components defined by elements in the
  TACSAssembler object.

  This call is collective - the number of components is obtained
  by a global reduction.
*/
int TACSAssembler::getNumComponents(){
  // Find the maximum component number on this processor
  int max_comp_num = 0;
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]->getComponentNum() >= max_comp_num){
      max_comp_num = elements[i]->getComponentNum();
    }
  }
  max_comp_num++;

  int num_components = 1;
  MPI_Allreduce(&max_comp_num, &num_components, 1, MPI_INT, 
                MPI_MAX, tacs_comm);
  return num_components;
}

/*
  Return the output nodal ranges. These may be used to determine what
  range of node numbers need to be determined by this process.  
*/
void TACSAssembler::getOutputNodeRange( enum ElementType elem_type,
					int ** _node_range ){
  int nelems = 0, nodes = 0, ncsr = 0;
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]->getElementType() == elem_type){
      elements[i]->addOutputCount(&nelems, &nodes, &ncsr);
    }
  }

  int * node_range = new int[ mpiSize+1 ];
  node_range[0] = 0;
  MPI_Allgather(&nodes, 1, MPI_INT, &node_range[1], 1, MPI_INT, tacs_comm);

  for ( int i = 0; i < mpiSize; i++ ){
    node_range[i+1] += node_range[i];
  }
  *_node_range = node_range;
}

/*
  Given the element type, determine the connectivity of the global
  data structure. Record the component number for each element within
  the data structure

  input:
  elem_type: the type of element eg SHELL, EULER_BEAM, SOLID etc.
  (see all the element types in Element.h)

  output:
  component_nums: an array of size nelems of the component numbers
  csr:            the csr element->nodes information for this class
  csr_range:      the range of csr data on this processor
  node_range:     the range of nodal values on this processor
*/
void TACSAssembler::getOutputConnectivity( enum ElementType elem_type,
                                           int ** component_nums,
                                           int ** _csr, int ** _csr_range, 
					   int ** _node_range ){
  // First go through and count up the number of elements and the 
  // size of the connectivity array required
  int nelems = 0, nodes = 0, ncsr = 0;
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]->getElementType() == elem_type){
      elements[i]->addOutputCount(&nelems, &nodes, &ncsr);
    }
  }

  int *node_range = new int[ mpiSize+1 ];
  int *csr_range = new int[ mpiSize+1 ];
  node_range[0] = 0;
  MPI_Allgather(&nodes, 1, MPI_INT, &node_range[1], 1, MPI_INT, tacs_comm);
  csr_range[0] = 0;
  MPI_Allgather(&ncsr, 1, MPI_INT, &csr_range[1], 1, MPI_INT, tacs_comm);

  for ( int i = 0; i < mpiSize; i++ ){
    node_range[i+1] += node_range[i];
    csr_range[i+1]  += csr_range[i];
  }

  int *csr = new int[ ncsr ];
  int *comp_nums = new int[ nelems ];
  ncsr = 0;
  nelems = 0;
  nodes = node_range[mpiRank];
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]->getElementType() == elem_type){
      elements[i]->getOutputConnectivity(&csr[ncsr], nodes);
      int n = nelems;
      elements[i]->addOutputCount(&nelems, &nodes, &ncsr);
      int c = elements[i]->getComponentNum();
      for ( ; n < nelems; n++ ){
        comp_nums[n] = c;
      }
    }
  }

  *component_nums = comp_nums;
  *_csr = csr;
  *_csr_range = csr_range;
  *_node_range = node_range;
}

/*
  Go through each element and get the output data for that element.

  The data is stored point-wise with each variable stored contiguously
  for each new point within the connectivity list. This stores the
  data at a point in memory indexed by data[node*nvals]. However,
  fewer than 'nvals' entries may be written in this call. The
  remaining data may be design variable entries that are computed
  below.

  elem_type: the element type to match
  out_type:  the output type 
  data:      the data array - nvals x the number of elements
  nvals:     the number of values to skip at each point
*/
void TACSAssembler::getOutputData( enum ElementType elem_type,
				   unsigned int out_type,
				   double * data, int nvals ){
  // Retrieve pointers to temporary storage
  TacsScalar *elemVars, *elemXpts;
  getDataPointers(elementData, &elemVars, NULL, NULL, NULL,
		  &elemXpts, NULL, NULL, NULL);
  
  int nelems = 0, nodes = 0, ncsr = 0;
  for ( int i = 0; i < numElements; i++ ){
    if (elements[i]->getElementType() == elem_type){
      getValues(varsPerNode, i, localVars, elemVars);
      getValues(TACS_SPATIAL_DIM, i, Xpts, elemXpts);      

      elements[i]->getOutputData(out_type, &data[nvals*nodes], nvals,
				 elemXpts, elemVars);
      elements[i]->addOutputCount(&nelems, &nodes, &ncsr);
    }
  }  
}
