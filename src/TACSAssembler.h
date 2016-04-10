#ifndef TACS_ASSEMBLER_H
#define TACS_ASSEMBLER_H

/*
  TACSAssembler assembles the residuals and matrices required for
  analysis and sensitivity analysis.

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
*/

class TACSAssembler; 

#include "TACSObject.h"
#include "TACSFunction.h"
#include "TACSElement.h"

#include "BVec.h" 
#include "BVecDist.h"
#include "DistMat.h"  
#include "FEMat.h"

class TACSAssembler : public TACSObject {
 public:
  static const int TACS_SPATIAL_DIM = 3;

  enum OrderingType { NATURAL_ORDER, // Natural ordering
                      RCM_ORDER, // Reverse Cuthill Mackee ordering
                      AMD_ORDER, // Approximate minimum degree
		      ND_ORDER, // Nested disection
                      TACS_AMD_ORDER }; // Interface variables ordered last
  enum MatrixOrderingType { ADDITIVE_SCHWARZ, 
                            APPROXIMATE_SCHUR,
                            DIRECT_SCHUR };

  TACSAssembler( MPI_Comm _tacs_comm,
		 int numOwnedNodes, int _varsPerNode,
		 int _numElements, int _numNodes,
		 int _nodeMaxCSRsize );
  TACSAssembler( MPI_Comm _tacs_comm, 
		 int numOwnedNodes, int _varsPerNode,
		 int _numElements, int _numNodes,
		 int _numDependentNodes, 
		 int _nodeMaxCSRsize );
  ~TACSAssembler();

  // Return important information about the TACS object
  // --------------------------------------------------
  int getNumNodes(){ return numNodes; }
  int getNumDependentNodes(){ return numDependentNodes; }
  int getNumElements(){ return numElements; }
  VarMap * getVarMap(){ return varMap; }

  // Add nodes to TACS
  // -----------------
  void addNode( int locaNodeNum, int tacsNodeNum );
  void addNodes( int localNodeNum[], int tacsNodeNum[], int numNodes );
  void addNodes( int **tacsNodeNums );

  // Add a dependent node to TACS
  // ----------------------------
  void setDependentNodes( int **_depNodeIndex, 
			  int **_depNodeToLocal,
			  double **_depNodeWeights );

  // Add elements to TACS
  // --------------------
  int addElement( TACSElement *element, int localNodes[], int numElemNodes );
  void addElements( TACSElement ***_elements );
  void setElementConnectivity( int **elemindex, int **elemnodes );

  // Functions that may be used to perform reordering externally
  // -----------------------------------------------------------
  void computeLocalNodeToNodeCSR( int **_rowp, int **_cols, int nodiag=0 );
  int computeCouplingNodes( int **_cnodes );
  int computeCouplingElements( int **_celems );
  
  // Functions for performing ordering based on RCM
  // ----------------------------------------------
  void computeReordering( enum OrderingType order_type, 
                          enum MatrixOrderingType mat_type );

  // Finalize the mesh - no further elements or 
  // nodes may be added following this call
  // ------------------------------------------
  void finalize();

  // Return the underlying TACS node numbers
  // ---------------------------------------
  void getTacsNodeNums( int localNodes[], int numNodes );
  int getTacsNodeNums( const int ** _tacsNodeNums );

  // Set the nodes in TACS 
  // ---------------------
  void setNodes( BVec *X ); 
  void getNodes( BVec *X );
  void getNodeArray( TacsScalar **_Xpts );

  // Associate a Dirichlet boundary condition with the given variables
  // -----------------------------------------------------------------
  void addBC( int nodeNum, const int bcNums[], int nbcs );
  void addBC( int nodeNum, const int bcNums[], 
	      const TacsScalar bcVals[], int nbcs );

  // Create vectors/matrices
  // -----------------------
  BVec * createVec();
  DistMat * createMat();
  FEMat * createFEMat( enum OrderingType order_type=TACS_AMD_ORDER );

  // Set/get the simulation time
  // ---------------------------
  void setSimulationTime( double _time );
  double getSimulationTime();

  // Retrieve the initial conditions for the simulation
  // --------------------------------------------------
  void getInitConditions( BVec *vars, BVec *dvars );

  // Methods for manipulating internal variable values
  // -------------------------------------------------
  void zeroVariables();
  void zeroDotVariables();
  void zeroDDotVariables();
  
  // Methods for setting/getting variables
  // -------------------------------------
  void getVariables( BVec *stateVars );
  void setVariables( BVec *stateVars );
  void setDotVariables( BVec *stateVars );
  void setDDotVariables( BVec *stateVars );

  // Residual and Jacobian assembly
  // ------------------------------
  void assembleRes( BVec *residual );
  void assembleResNoBCs( BVec *residual );
  void assembleJacobian( BVec *residual, TACSMat *A,
			 double alpha, double beta, double gamma,
			 MatrixOrientation matOr=NORMAL );
  void assembleMatType( ElementMatrixType matType,
			TACSMat *A, MatrixOrientation matOr=NORMAL );

  // Design variable handling
  // ------------------------
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lowerBound[], 
			  TacsScalar upperBound[], int numDVs );

  // Function and sensitivity evaluation 
  // -----------------------------------
  void evalFunctions( TACSFunction **funcs, int numFuncs,
                      TacsScalar *funcVals );
  void evalDVSens( TACSFunction **funcs, int numFuncs, 
                   TacsScalar *fdvSens, int numDVs );
  void evalSVSens( TACSFunction *function, BVec *vec );
  void evalAdjointResProducts( BVec **adjoint, int numAdjoints,
                               TacsScalar * dvSens, int numDVs );
  // void evalXptSens( TACSFunction **funcs, int numFuncs,
  //                   TACSVec *fXptSens );
  // void evalAdjointResXptSensProducts( BVec ** adjoint, int numAdjoints,
  //   				      TACSVec *adjXptSensProduct );

  // Evaluate the derivative of the inner product of two vectors and a matrix
  // ------------------------------------------------------------------------
  void evalMatDVSensInnerProduct( TacsScalar scale, 
				  ElementMatrixType matType, 
				  BVec *psi, BVec *phi,
				  TacsScalar *dvSens, int numDVs );

  // Evaluate the partial derivative of the inner product with a matrix
  // ------------------------------------------------------------------
  void evalMatSVSensInnerProduct( TacsScalar scale,
				  ElementMatrixType matType, 
				  BVec *psi, BVec *phi, BVec *res );

  // Return an element and the variables associated with that element
  // ----------------------------------------------------------------
  TACSElement ** getElements(){ return elements; }
  TACSElement * getElement( int elemNum,
			    TacsScalar *elemXpts,
			    TacsScalar *vars,
			    TacsScalar *dvars,
			    TacsScalar *ddvars );

  // Get information about the output files
  // --------------------------------------
  int getNumComponents();
  void getOutputNodeRange( enum ElementType elem_type, 
			   int ** _node_range );
  void getOutputConnectivity( enum ElementType elem_type,
                              int ** _component_nums,
			      int ** _csr, int ** _csr_range, 
			      int ** _node_range );
  void getOutputData( enum ElementType elem_type,
		      unsigned int out_type,
		      double * data, int nvals );

  // Test the given element, constitutive or function class
  // ------------------------------------------------------
  void testElement( int elemNum, int print_level );
  void testConstitutive( int elemNum, int print_level );
  void testFunction( TACSFunction * func, 
                     int num_design_vars, double dh );

  // Retrieve the MPI communicator
  // -----------------------------
  MPI_Comm getMPIComm();

  // Set the number of threads to work with
  // --------------------------------------
  void setNumThreads( int t );

 private:
  // Contruct the TACSAssembler object, given the initial information
  // ----------------------------------------------------------------
  void init( MPI_Comm _tacs_comm, int numOwnedNodes, int _varsPerNode,
	     int _numElements, int _numNodes, int _numDependentNodes,
	     int _nodeMaxCSRsize );

  // Set the dependent nodal values based on the independent nodes
  // -------------------------------------------------------------
  void setDependentVariables( const int perNode, TacsScalar * vars );
  void addDependentResidual( const int perNode, TacsScalar * vars );

  // Get pointers to the start-locations within the data array
  // ---------------------------------------------------------
  void getDataPointers( TacsScalar *data, 
			TacsScalar **v1, TacsScalar **v2, 
			TacsScalar **v3, TacsScalar **v4,
			TacsScalar **x1, TacsScalar **x2,
			TacsScalar **weights, TacsScalar **mat );

  // Apply the boundary conditions to the residual vector
  // ----------------------------------------------------
  void applyBCs( BVec * res, const TacsScalar vars[] );  

  // Functions for ordering the variables
  // ------------------------------------
  void computeNodeToElementCSR( int ** _nodeElem, int ** _nodeElemIndex );
  void computeLocalNodeToNodeCSR( int ** _rowp, int ** _cols, 
				  int nrnodes, const int * rnodes,
				  int nodiag );

  // Compute the reordering for a sub-matrix
  // ---------------------------------------
  void computeMatReordering( enum OrderingType order_type, 
                             int nvars, int * rowp, int * cols,
                             int * perm, int * new_vars );

  // Initialize the internal arrays for storing load case information
  // ----------------------------------------------------------------
  void initializeArrays();

  inline int getValues( const int perNode, const int elemNum, 
			const TacsScalar *local, TacsScalar *vals );
  inline int addValues( const int perNode, const int elemNum, 
			const TacsScalar *vals, TacsScalar *local );
  inline int setValues( const int perNode, const int elemNum,
                        const TacsScalar *vals, TacsScalar *local );
  inline void addMatValues( TACSMat *A, const int elemNum, 
			    const TacsScalar *mat,
			    int *item, TacsScalar *temp );

  // Initialize the functions in the list if they have not been 
  // initialized already
  void initializeFunctions( TACSFunction ** functions, int numFuncs );

  // The static member functions that are used to p-thread TACSAssembler
  // operations... These are the most time-consuming operations.
  static void schedPthreadJob( TACSAssembler * tacs,
                               int * index, int total_size );
  static void * assembleRes_thread( void * t );
  static void * assembleJacobian_thread( void * t );
  static void * assembleMatType_thread( void * t );
  // static void * adjointResXptSensProduct_thread( void * t );
  static void * adjointResProduct_thread( void * t );
  static void * evalFunctions_thread( void * t );
  // static void * evalXptSens_thread( void * t );
  static void * evalDVSens_thread( void * t );

  VarMap * varMap; // Variable ownership map
  BCMap * bcMap; // Boundary condition data
  BVecDistribute * vecDist; // Distribute the vector
  BVecIndices * vecDistIndices; // The tacsVarNum indices

  // Additional, persistent information for the FEMat class
  // These are created once - all subsequent calls use this data.
  BVecIndices *feMatBIndices, *feMatCIndices;
  BVecDistribute *feMatBMap, *feMatCMap;

  // The global simulation time variable
  double time;

  // variables/elements have been initialized
  int meshFinalizedFlag;

  int varsPerNode; // number of variables per node
  int maxElementNodes; // maximum number of ind. and dep. nodes for any element
  int maxElementSize; // maximum number of variables for any element
  int maxElementIndepNodes; // maximum number of independent nodes 
  int numElements; // number of elements
  int numNodes; // number of nodes referenced by this process
  int numDependentNodes; // number of dependent nodes
  int *tacsNodeNums; // node numbers associated with TACS

  // Variables that define the CSR data structure to 
  // store the element -> node information
  int *elementNodeIndex;
  int *elementLocalNodes;
  int *elementTacsNodes;

  int nodeMaxCSRsize; // the maximum size of the elementLocalVars array
  int nodeCSRIncrement; // increment the size of the csr structure by this

  // Variables that define the dependent node to independent node
  // dependence
  int *depNodeIndex;
  int *depNodeToLocal;
  int *depNodeToTacs;
  double *depNodeWeights;

  // For use during set up - before call to finalize
  int currElement; // Number of elements currently added (max val. numElements)
  int currNode; // Number of nodes currently added (max val. numNodes)

  // The local list of elements
  TACSElement **elements;

  // Memory for the variables/residuals
  TacsScalar *localVars, *localDotVars, *localDDotVars;
  TacsScalar *localRes; // Local residual values being assembled

  // Memory for the element residuals and variables
  TacsScalar *elementData; // Space for element residuals/matrices
  int *elementIData; // Space for element index data

  // The x,y,z positions/sensitivities of all the local nodes
  TacsScalar *Xpts; // The nodal locations

  // The data required to perform parallel operations
  // MPI info
  int mpiRank, mpiSize;
  MPI_Comm tacs_comm;

  // Class to store specific information about the threaded
  // operations to perform. Note that assembly operations are
  // relatively easy, while design-variable dependent info is
  // much more challenging!
  class TACSAssemblerPthreadInfo {
  public:
    TACSAssemblerPthreadInfo(){
      tacs = NULL; 
      // Matrix information
      mat = NULL;
      alpha = beta = gamma = 0.0;
      matType = STIFFNESS_MATRIX;
      matOr = NORMAL;
      // Function data
      funcIteration = 0;
      numFuncs = 0;
      functions = NULL;
      // df/dx and adjoint-dR/dx product data
      numDesignVars = 0;
      numAdjoints = 0;
      fdvSens = NULL;
      fXptSens = NULL;
      adjointVars = NULL;
      adjXptSensProduct = NULL;
    }

    // The data required to perform most of the matrix
    // assembly.
    TACSAssembler *tacs;

    // Information for matrix assembly
    TACSMat *mat;
    double alpha, beta, gamma;
    ElementMatrixType matType;
    MatrixOrientation matOr;

    // Information required for the computation of f or df/dx
    unsigned int funcIteration;
    int numDesignVars;
    int numFuncs;
    TACSFunction **functions;
    TacsScalar *fdvSens; // df/dx
    TacsScalar *fXptSens; // df/dXpts

    // Information for adjoint-dR/dx products
    int numAdjoints;
    TacsScalar *adjointVars;
    TacsScalar *adjXptSensProduct;
  } *tacsPInfo;

  // The pthread data required to pthread tacs operations
  int numCompletedElements; // Keep track of how much work has been done
  TACSThreadInfo * thread_info;// The pthread object
  pthread_t threads[TACSThreadInfo::TACS_MAX_NUM_THREADS]; // The thread objects
  pthread_mutex_t tacs_mutex; // The mutex object for coordinating assembly ops.

  // The name of the TACSAssembler object
  static const char * tacsName;
};

/*!  
  Get the values associated with an element from the local array.

  These are private functions and therefore no bounds checking is
  performed.

  local:  the local values of a vector (input)
  nnodes: the number of nodes
  nodes:  the node numbers
  vals:   the values from each node (output)
*/
inline int TACSAssembler::getValues( const int perNode,
				     const int elemNum, 
				     const TacsScalar * local, 
				     TacsScalar * vals ){
  int start = elementNodeIndex[elemNum];
  int nnodes = elementNodeIndex[elemNum+1] - start;
  const int *nodes = &elementLocalNodes[start];

  if (numDependentNodes == 0){
    switch (perNode){
    case 1:
      for ( int i = 0; i < nnodes; i++ ){
	int j = nodes[i];
	vals[0] = local[j];
	vals++;
      }
      break;
    case 2:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 2*nodes[i];
	vals[0] = local[j]; 
	vals[1] = local[j+1];
	vals += 2;
      }
      break;
    case 3:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 3*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals += 3;
      }
      break;
    case 4:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 4*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals[2] = local[j+3];
	vals += 4;
      }
      break;
    case 5:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 5*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals[1] = local[j+3];
	vals[2] = local[j+4];
	vals += 5;
      }
      break;
    case 6:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 6*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals[3] = local[j+3];
	vals[4] = local[j+4];
	vals[5] = local[j+5];
	vals += 6;
      }
      break;
    case 7:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 7*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals[3] = local[j+3];
	vals[4] = local[j+4];
	vals[6] = local[j+6];
	vals += 7;
      }
      break;
    case 8:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 8*nodes[i];
	vals[0] = local[j];
	vals[1] = local[j+1];
	vals[2] = local[j+2];
	vals[3] = local[j+3];
	vals[4] = local[j+4];
	vals[5] = local[j+5];
        vals[6] = local[j+6];
        vals[7] = local[j+7];
	vals += 8;
      }
      break;
    default:
      for ( int i = 0; i < nnodes; i++ ){
	int j = perNode*nodes[i];
	for ( int n = 0; n < perNode; n++, j++ ){
	  vals[n] = local[j];
	}
	vals += perNode;
      }
      break;
    }
  }
  else {
    switch (perNode){
    case 1:
      for ( int i = 0; i < nnodes; i++ ){
	int j = nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	}
	else {
	  j = numNodes-1-j;
	  vals[0] = local[j];
	}
	vals++;
      }
      break;
    case 2:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 2*nodes[i];
	if (j >= 0){
	  vals[0] = local[j]; 
	  vals[1] = local[j+1];
	}
	else {
	  j = 2*(numNodes-1)-j;
	  vals[0] = local[j]; 
	  vals[1] = local[j+1];
	}
	vals += 2;
      }
      break;
    case 3:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 3*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	}
	else {
	  j = 3*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	}
	vals += 3;
      }
      break;
    case 4:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 4*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[2] = local[j+3];
	}
	else {
	  j = 4*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[2] = local[j+3];
	}
	vals += 4;
      }
      break;
    case 5:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 5*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[1] = local[j+3];
	  vals[2] = local[j+4];
	}
	else {
	  j = 5*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[1] = local[j+3];
	  vals[2] = local[j+4];
	}
	vals += 5;
      }
      break;
    case 6:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 6*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
	}
	else {
	  j = 6*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
	}
	vals += 6;
      }
      break;
    case 7:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 7*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
          vals[6] = local[j+6];
	}
	else {
	  j = 7*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
          vals[6] = local[j+6];
	}
	vals += 7;
      }
      break;
    case 8:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 8*nodes[i];
	if (j >= 0){
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
          vals[6] = local[j+6];
          vals[7] = local[j+7];
	}
	else {
	  j = 8*(numNodes-1)-j;
	  vals[0] = local[j];
	  vals[1] = local[j+1];
	  vals[2] = local[j+2];
	  vals[3] = local[j+3];
	  vals[4] = local[j+4];
	  vals[5] = local[j+5];
          vals[6] = local[j+6];
          vals[7] = local[j+7];
	}
	vals += 8;
      }
      break;
    default:
      for ( int i = 0; i < nnodes; i++ ){
	int j = perNode*nodes[i];
	if (j >= 0){
	  for ( int n = 0; n < perNode; n++, j++ ){
	    vals[n] = local[j];
	  }
	}
	else {
	  j = perNode*(numNodes-1)-j;
	  for ( int n = 0; n < perNode; n++, j++ ){
	    vals[n] = local[j];
	  }
	}
	vals += perNode;
      }
      break;
    }
  }

  return nnodes;
}

/*! 
  Add values associated with an element into the local array.

  These are private functions and therefore no bounds checking is
  performed.

  local:  the local values of a vector (output)
  nnodes: the number of nodes
  nodes:  the node numbers
  vals:   the values to set for each node (input)
*/
inline int TACSAssembler::addValues( const int perNode,
				     const int elemNum, 
				     const TacsScalar * vals,
				     TacsScalar * local ){
  int start = elementNodeIndex[elemNum];
  int nnodes = elementNodeIndex[elemNum+1] - start;
  const int *nodes = &elementLocalNodes[start];

  if (numDependentNodes == 0){
    switch (perNode){
    case 1:
      for ( int i = 0; i < nnodes; i++ ){
	int j = nodes[i];
	local[j] += vals[0]; 
	vals++;
      }
      break;
    case 2:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 2*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	vals += 2;
      }
      break;
    case 3:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 3*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	vals += 3;
      }
      break;
    case 4:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 4*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	local[j+3] += vals[3]; 
	vals += 4;
      }
      break;
    case 5:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 5*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	local[j+3] += vals[3]; 
	local[j+4] += vals[4]; 
	vals += 5;
      }
      break;
    case 6:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 6*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	local[j+3] += vals[3]; 
	local[j+4] += vals[4]; 
	local[j+5] += vals[5]; 
	vals += 6;
      }
      break;
    case 7:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 7*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	local[j+3] += vals[3]; 
	local[j+4] += vals[4]; 
	local[j+5] += vals[5]; 
	local[j+6] += vals[6]; 
	vals += 7;
      }
      break;
    case 8:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 8*nodes[i];
	local[j] += vals[0]; 
	local[j+1] += vals[1]; 
	local[j+2] += vals[2]; 
	local[j+3] += vals[3]; 
	local[j+4] += vals[4]; 
	local[j+5] += vals[5]; 
        local[j+6] += vals[6]; 
        local[j+7] += vals[7]; 
	vals += 8;
      }
      break;
    default:
      for ( int i = 0; i < nnodes; i++ ){
	int j = perNode*nodes[i];
	for ( int n = 0; n < perNode; n++, j++ ){
	  local[j] += vals[n];
	}
	vals += perNode;
      }
      break;
    }
  }
  else {
    switch (perNode){
    case 1:
      for ( int i = 0; i < nnodes; i++ ){
	int j = nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	}
	else {
	  j = numNodes-1-j;
	  local[j] += vals[0]; 
	}
	vals++;
      }
      break;
    case 2:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 2*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	}
	else {
	  j = 2*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	}
	vals += 2;
      }
      break;
    case 3:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 3*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	}
	else {
	  j = 3*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	}
	vals += 3;
      }
      break;
    case 4:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 4*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	}
	else {
	  j = 4*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	}
	vals += 4;
      }
      break;
    case 5:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 5*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	}
	else {
	  j = 5*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	}
	vals += 5;
      }
      break;
    case 6:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 6*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
	}
	else {
	  j = 6*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
	}
	vals += 6;
      }
      break;
    case 7:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 7*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
          local[j+6] += vals[6]; 
	}
	else {
	  j = 7*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
          local[j+6] += vals[6]; 
	}
	vals += 7;
      }
      break;
    case 8:
      for ( int i = 0; i < nnodes; i++ ){
	int j = 8*nodes[i];
	if (j >= 0){
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
          local[j+6] += vals[6]; 
          local[j+7] += vals[7]; 
	}
	else {
	  j = 8*(numNodes-1)-j;
	  local[j] += vals[0]; 
	  local[j+1] += vals[1]; 
	  local[j+2] += vals[2]; 
	  local[j+3] += vals[3]; 
	  local[j+4] += vals[4]; 
	  local[j+5] += vals[5]; 
          local[j+6] += vals[6]; 
          local[j+7] += vals[7]; 
	}
	vals += 8;
      }
      break;
    default:
      for ( int i = 0; i < nnodes; i++ ){
	int j = perNode*nodes[i];
	if (j >= 0){
	  for ( int n = 0; n < perNode; n++, j++ ){
	    local[j] += vals[n];
	  }
	}
	else {
	  j = perNode*(numNodes-1)-j;
	  for ( int n = 0; n < perNode; n++, j++ ){
	    local[j] += vals[n];
	  }
	}
	vals += perNode;
      }
      break;
    }
  }

  return nnodes;
}

/*!  
  Set the values associated with an element from the local array.

  These are private functions and therefore no bounds checking is
  performed.

  nnodes: the number of nodes
  nodes:  the node numbers
  vals:   the values from each node
  local:  the local values of a vector output
*/
inline int TACSAssembler::setValues( const int perNode,
				     const int elemNum, 
				     const TacsScalar *vals,
				     TacsScalar *local ){
  int start = elementNodeIndex[elemNum];
  int nnodes = elementNodeIndex[elemNum+1] - start;
  const int *nodes = &elementLocalNodes[start];

  switch (perNode){
  case 1:
    for ( int i = 0; i < nnodes; i++ ){
      int j = nodes[i];
      if (j >= 0){
	local[j] = vals[0];
      }
      vals++;
    }
    break;
  case 2:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 2*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
      }
      vals += 2;
    }
    break;
  case 3:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 3*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
      }
      vals += 3;
    }
    break;
  case 4:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 4*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
        local[j+3] = vals[3];
      }
      vals += 4;
    }
    break;
  case 5:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 5*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
        local[j+3] = vals[3];
        local[j+4] = vals[4];
      }
      vals += 5;
    }
    break;
  case 6:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 6*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
        local[j+3] = vals[3];
        local[j+4] = vals[4];
        local[j+5] = vals[5];
      }
      vals += 6;
    }
    break;
  case 7:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 7*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
        local[j+3] = vals[3];
        local[j+4] = vals[4];
        local[j+5] = vals[5];
        local[j+6] = vals[6];
      }
      vals += 7;
    }
    break;
  case 8:
    for ( int i = 0; i < nnodes; i++ ){
      int j = 8*nodes[i];
      if (j >= 0){
        local[j] = vals[0];
        local[j+1] = vals[1];
        local[j+2] = vals[2];
        local[j+3] = vals[3];
        local[j+4] = vals[4];
        local[j+5] = vals[5];
        local[j+6] = vals[6];
        local[j+7] = vals[7];
      }
      vals += 8;
    }
    break;
  default:
    for ( int i = 0; i < nnodes; i++ ){
      int j = perNode*nodes[i];
      if (j >= 0){
        for ( int n = 0; n < perNode; n++, j++ ){
          local[j] = vals[n];
        }
      }
      vals += perNode;
    }
    break;
  }

  return nnodes;
}

/*
  Add the values of the element matrix to the provided TACSMat. 

  This code takes into account dependent-nodes (when they exist) by
  adding the inner product of the dependent weights with the element
  matrix.  Note that the integer and scalar temporary storage should
  be allocated once for all elements for efficiency purposes. The
  maximum weight length can be determined by finding the maximum
  number of nodes + max total number of dependent->local nodes in any
  local element.

  input:
  elemNum:    the element number 
  mat:        the corresponding element matrix
  itemp:      temporary integer storage len(itemp) >= nnodes+1 + len(vars)
  temp:       temporary scalar storage len(temp) >= len(weights)

  output:
  A:          the matrix to which the element-matrix is added
*/
inline void TACSAssembler::addMatValues( TACSMat * A, 
					 const int elemNum, 
					 const TacsScalar * mat,
					 int * itemp, TacsScalar * temp ){
  int start = elementNodeIndex[elemNum];
  int end = elementNodeIndex[elemNum+1];
  int nnodes = end - start;
  int nvars = varsPerNode*nnodes;

  // Add the element values to the matrix
  const int *nodeNums = &elementTacsNodes[start];

  if (numDependentNodes == 0){
    // If we have no dependent nodes, then we don't need to do
    // anything extra here
    A->addValues(nnodes, nodeNums, nnodes, nodeNums, 
		 nvars, nvars, mat);
  }
  else {
    // If we have dependent nodes, then we have to figure out what
    // the weighting matrix is and add then add the element matrix
    
    // Set pointers to the temporary arrays
    int *varp = &itemp[0];
    int *vars = &itemp[nnodes+1];
    TacsScalar *weights = temp;

    varp[0] = 0;
    for ( int i = 0, k = 0; i < nnodes; i++ ){
      if (nodeNums[i] >= 0){
	// This is just a regular node
	weights[k] = 1.0;
	vars[k] = nodeNums[i];
	k++;
      }
      else {
	// This is a dependent node. Determine the corresponding
	// dependent node number and add the variables
	int dep = -nodeNums[i]-1;
	for ( int j = depNodeIndex[dep]; 
	      j < depNodeIndex[dep+1]; j++, k++ ){
	  weights[k] = depNodeWeights[j];
	  vars[k] = depNodeToTacs[j];
	}
      }

      varp[i+1] = k;
    }

    // Add the values to the matrix
    A->addWeightValues(nnodes, varp, vars, weights,
		       nvars, nvars, mat);
  }
}

#endif
