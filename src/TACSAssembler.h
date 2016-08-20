#ifndef TACS_ASSEMBLER_H
#define TACS_ASSEMBLER_H

/*
  TACSAssembler assembles the residuals and matrices required for
  analysis and sensitivity analysis.

  Copyright (c) 2010-2016 Graeme Kennedy. All rights reserved. 
*/

class TACSAssembler; 

// Basic analysis classes
#include "TACSObject.h"
#include "TACSFunction.h"
#include "TACSElement.h"
#include "TACSAuxElements.h"

// Linear algebra classes
#include "BVec.h"
#include "BVecDist.h"
#include "DistMat.h"
#include "FEMat.h"

/*
  TACSAssembler

  This is the main class required for structural analysis using TACS.
  Basic operations required for analysis and design optimization
  should be performed using this class, rather than with element-level
  functions. This class is used to assemble residuals and Jacobians
  required for analysis. The class also implements the routines
  necessary to compute the adjoint. These operations include parallel
  function evaluation, derivative evaluation and the implementation of
  terms required for the adjoint method.

  TACSAssembler can be instantiated and initialized directly through
  API calls, or through other TACS objects which perform the
  initialization process.  Once the TACSAssembler object is
  initialized, however, subsequent calls to the code work regardless
  of how the object was created. In addition, once TACSAssembler is
  created, the parallelism is transparent. All analysis calls are
  collective on the tacs_comm communicator.
*/
class TACSAssembler : public TACSObject {
 public:
  // There are always 3 coordinates (even for 2D problems)
  static const int TACS_SPATIAL_DIM = 3;

  enum OrderingType { NATURAL_ORDER, // Natural ordering
                      RCM_ORDER, // Reverse Cuthill Mackee ordering
                      AMD_ORDER, // Approximate minimum degree
		      ND_ORDER, // Nested disection
                      TACS_AMD_ORDER }; // Interface variables ordered last
  enum MatrixOrderingType { ADDITIVE_SCHWARZ, 
                            APPROXIMATE_SCHUR,
                            DIRECT_SCHUR };
  
  // Create the TACSAssembler object in parallel
  // -------------------------------------------
  TACSAssembler( MPI_Comm _tacs_comm, int _varsPerNode,
		 int _numOwnedNodes, int _numElements, 
		 int _numDependentNodes );
  ~TACSAssembler();

  // Set the connectivity in TACS
  // ----------------------------
  int setElementConnectivity( const int *conn, const int *ptr );
  int setElements( TACSElement **_elements );
  int setDependentNodes( const int *_depNodeIndex, 
                         const int *_depNodeToTacs,
                         const double *_depNodeWeights );

  // Associate a Dirichlet boundary condition with the given variables
  // -----------------------------------------------------------------
  void addBCs( int nnodes, const int *nodes, 
               int nbcs=-1, const int *vars=NULL, 
               const TacsScalar *vals=NULL );
  
  // Functions for performing ordering based on RCM
  // ----------------------------------------------
  // void computeReordering( enum OrderingType order_type, 
  //                         enum MatrixOrderingType mat_type );

  // Initialize the mesh
  // -------------------
  int initialize();

  // Return important information about the TACS object
  // --------------------------------------------------
  int getVarsPerNode();
  int getNumNodes();
  int getNumDependentNodes();
  int getNumElements();
  TACSVarMap *getVarMap();
  TACSBcMap *getBcMap();
  TACSBVecDistribute *getBVecDistribute();

  // Set auxiliary elements into the TACSAssembler object
  // ----------------------------------------------------
  void setAuxElements( TACSAuxElements *aux_elems );
  TACSAuxElements *getAuxElements();

  // Set the nodes in TACS 
  // ---------------------
  TACSBVec *createNodeVec();
  void setNodes( TACSBVec *X ); 
  void getNodes( TACSBVec *X );

  // Set/get the simulation time
  // ---------------------------
  void setSimulationTime( double _time );
  double getSimulationTime();

  // Create vectors
  // --------------
  TACSBVec *createVec();

  // Methods for manipulating internal variable values
  // -------------------------------------------------
  void zeroVariables();
  void zeroDotVariables();
  void zeroDDotVariables();
  
  // Methods for setting/getting variables
  // -------------------------------------
  void setVariables( TACSBVec *q, 
                     TACSBVec *qdot=NULL, TACSBVec *qddot=NULL );
  void getVariables( TACSBVec *q, 
                     TACSBVec *qdot=NULL, TACSBVec *qddot=NULL );

  // Create the matrices that can be used for analysis
  // -------------------------------------------------
  DistMat *createMat();
  FEMat *createFEMat( enum OrderingType order_type=TACS_AMD_ORDER );

  // Retrieve the initial conditions for the simulation
  // --------------------------------------------------
  void getInitConditions( TACSBVec *vars, TACSBVec *dvars );

  // Evaluate the kinetic and potential energy
  // -----------------------------------------
  void evalEnergies( TacsScalar *Te, TacsScalar *Pe );

  // Residual and Jacobian assembly
  // ------------------------------
  void assembleRes( TACSBVec *residual );
  void assembleJacobian( TACSBVec *residual, TACSMat *A,
			 double alpha, double beta, double gamma,
			 MatrixOrientation matOr=NORMAL );
  void assembleMatType( ElementMatrixType matType,
			TACSMat *A, MatrixOrientation matOr=NORMAL );
  void addJacobianVecProduct( TacsScalar scale, 
                              double alpha, double beta, double gamma,
                              TACSBVec *x, TACSBVec *y,
                              MatrixOrientation matOr=NORMAL );

  // Design variable handling
  // ------------------------
  void getDesignVars( TacsScalar dvs[], int numDVs );
  void setDesignVars( const TacsScalar dvs[], int numDVs );
  void getDesignVarRange( TacsScalar lb[], TacsScalar ub[], int numDVs );

  // Function and sensitivity evaluation
  // -----------------------------------
  void evalFunctions( TACSFunction **funcs, int numFuncs,
                      TacsScalar *funcVals );
  void addDVSens( TACSFunction **funcs, int numFuncs,
                  TacsScalar *fdvSens, int numDVs );
  void addSVSens( TACSFunction *function, TACSBVec *vec );
  void addXptSens( TACSFunction **funcs, int numFuncs,
                   TACSBVec **fXptSens );

  void addAdjointResProducts( TACSBVec **adjoint, int numAdjoints,
                              TacsScalar *dvSens, int numDVs );
  /*
  void evalAdjointResXptSensProducts( TACSBVec **adjoint, int numAdjoints,
    				      TACSBVec *adjXptSensProduct );
  */

  // Evaluate the derivative of the inner product of two vectors and a matrix
  // ------------------------------------------------------------------------
  void addMatDVSensInnerProduct( TacsScalar scale, 
                                 ElementMatrixType matType, 
                                 TACSBVec *psi, TACSBVec *phi,
                                 TacsScalar *dvSens, int numDVs );

  // Evaluate the partial derivative of the inner product with a matrix
  // ------------------------------------------------------------------
  void evalMatSVSensInnerProduct( TacsScalar scale,
				  ElementMatrixType matType, 
				  TACSBVec *psi, TACSBVec *phi, TACSBVec *res );

  // Return an element and the variables associated with that element
  // ----------------------------------------------------------------
  TACSElement **getElements(){ return elements; }

  // Get information about the output files
  // --------------------------------------
  int getNumComponents();
  void getOutputNodeRange( enum ElementType elem_type, 
			   int **_node_range );
  void getOutputConnectivity( enum ElementType elem_type,
                              int **_component_nums,
			      int **_csr, int **_csr_range, 
			      int **_node_range );
  void getOutputData( enum ElementType elem_type,
		      unsigned int out_type,
		      double *data, int nvals );

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

  // Get a constant pointer to the tacs node numbers
  // -----------------------------------------------
  int getTacsNodeNums( const int **_tacsNodeNums );

 private:
  // Get pointers to the start-locations within the data array
  // ---------------------------------------------------------
  void getDataPointers( TacsScalar *data, 
			TacsScalar **v1, TacsScalar **v2, 
			TacsScalar **v3, TacsScalar **v4,
			TacsScalar **x1, TacsScalar **x2,
			TacsScalar **weights, TacsScalar **mat );

  // Functions that are used to perform reordering
  // ---------------------------------------------
  int computeExtNodes();
  int getLocalNodeNum( int node );
  int getGlobalNodeNum( int node );
  void computeLocalNodeToNodeCSR( int **_rowp, int **_cols, int nodiag=0 );
  int computeCouplingNodes( int **_cnodes );
  int computeCouplingElements( int **_celems );

  // Functions for ordering the variables
  // ------------------------------------
  void computeNodeToElementCSR( int **_nodeElem, int **_nodeElemIndex );
  void computeLocalNodeToNodeCSR( int **_rowp, int **_cols, 
				  int nrnodes, const int *rnodes,
				  int nodiag );

  // Compute the reordering for a sub-matrix
  // ---------------------------------------
  void computeMatReordering( enum OrderingType order_type, 
                             int nvars, int *rowp, int *cols,
                             int *perm, int *new_vars );

  // Add values into the matrix
  inline void addMatValues( TACSMat *A, const int elemNum, 
			    const TacsScalar *mat,
			    int *item, TacsScalar *temp,
                            MatrixOrientation matOr );

  // Initialize the functions in the list if they have not been 
  // initialized already
  void initializeFunctions( TACSFunction **functions, int numFuncs );

  // The static member functions that are used to p-thread TACSAssembler
  // operations... These are the most time-consuming operations.
  static void schedPthreadJob( TACSAssembler *tacs,
                               int *index, int total_size );
  static void *assembleRes_thread( void *t );
  static void *assembleJacobian_thread( void *t );
  static void *assembleMatType_thread( void *t );
  // static void *adjointResXptSensProduct_thread( void *t );
  static void *adjointResProduct_thread( void *t );
  static void *evalFunctions_thread( void *t );
  // static void *evalXptSens_thread( void *t );
  static void *evalDVSens_thread( void *t );

  TACSVarMap *varMap; // Variable ownership map
  TACSBcMap *bcMap; // Boundary condition data
  TACSBVecDistribute *extDist; // Distribute the vector
  TACSBVecIndices *extDistIndices; // The tacsVarNum indices

  // Additional information information for the DistMat class
  TACSBVecIndices *distMatIndices;

  // Additional ordering information for the FEMat class
  // These are created once - all subsequent calls use this data.
  TACSBVecIndices *feMatBIndices, *feMatCIndices;
  TACSBVecDistribute *feMatBMap, *feMatCMap;

  // The global simulation time variable
  double time;

  // variables/elements have been initialized
  int meshInitializedFlag;

  // Information about the
  int varsPerNode; // number of variables per node
  int maxElementNodes; // maximum number of ind. and dep. nodes for any element
  int maxElementSize; // maximum number of variables for any element
  int maxElementIndepNodes; // maximum number of independent nodes 
  int numElements; // number of elements
  int numNodes; // number of nodes referenced by this process
  int numOwnedNodes; // number of nodes owned by this processor
  int numExtNodes; // number of extneral nodes 
  int numDependentNodes; // number of dependent nodes

  // Node numbers that are referred to from this processor
  int *tacsExtNodeNums; // node numbers associated with TACS
  int extNodeOffset; // Offset into the external nodes

  // Variables that define the CSR data structure to 
  // store the element -> node information
  int *elementNodeIndex;
  int *elementTacsNodes;

  // Variables that define the dependent node to independent node
  // dependence
  TACSBVecDepNodes *depNodes;

  // The local list of elements
  TACSElement **elements;

  // The auxiliary element class
  TACSAuxElements *auxElements;

  // The variables, velocities and accelerations
  TACSBVec *varsVec, *dvarsVec, *ddvarsVec;

  // Memory for the node locations
  TACSBVec *xptVec;

  // Memory for the element residuals and variables
  TacsScalar *elementData; // Space for element residuals/matrices
  int *elementIData; // Space for element index data

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
  TACSThreadInfo *thread_info;// The pthread object
  pthread_t threads[TACSThreadInfo::TACS_MAX_NUM_THREADS]; // The thread objects
  pthread_mutex_t tacs_mutex; // The mutex object for coordinating assembly ops.

  // The name of the TACSAssembler object
  static const char *tacsName;
};

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

  input/output:
  A:          the matrix to which the element-matrix is added
*/
inline void TACSAssembler::addMatValues( TACSMat *A, 
					 const int elemNum, 
					 const TacsScalar *mat,
					 int *itemp, TacsScalar *temp,
                                         MatrixOrientation matOr ){
  int start = elementNodeIndex[elemNum];
  int end = elementNodeIndex[elemNum+1];
  int nnodes = end - start;
  int nvars = varsPerNode*nnodes;

  // Add the element values to the matrix
  const int *nodeNums = &elementTacsNodes[start];

  if (matOr == NORMAL && numDependentNodes == 0){
    // If we have no dependent nodes, then we don't need to do
    // anything extra here
    A->addValues(nnodes, nodeNums, nnodes, nodeNums, 
		 nvars, nvars, mat);
  }
  else {
    // If we have dependent nodes, then we have to figure out what
    // the weighting matrix is and add then add the element matrix
    const int *depNodePtr = NULL;
    const int *depNodeConn = NULL;
    const double *depNodeWeights = NULL;
    if (depNodes){
      depNodes->getDepNodes(&depNodePtr, &depNodeConn, 
                            &depNodeWeights);
    }
    
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
	for ( int j = depNodePtr[dep]; j < depNodePtr[dep+1]; j++, k++ ){
	  weights[k] = depNodeWeights[j];
	  vars[k] = depNodeConn[j];
	}
      }

      varp[i+1] = k;
    }

    // Add the values to the matrix
    A->addWeightValues(nnodes, varp, vars, weights,
		       nvars, nvars, mat, matOr);
  }
}

#endif
