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

#ifndef TACS_ASSEMBLER_H
#define TACS_ASSEMBLER_H

/*
  TACSAssembler assembles the residuals and matrices required for
  analysis and sensitivity analysis.
*/

class TACSAssembler;

// Basic analysis classes
#include "TACSAuxElements.h"
#include "TACSElement.h"
#include "TACSFunction.h"
#include "TACSObject.h"

// Linear algebra classes
#include "TACSBVecDistribute.h"
#include "TACSParallelMat.h"
#include "TACSSchurMat.h"
#include "TACSSerialPivotMat.h"

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

  enum OrderingType {
    NATURAL_ORDER,   // Natural ordering
    RCM_ORDER,       // Reverse Cuthill Mackee ordering
    AMD_ORDER,       // Approximate minimum degree
    ND_ORDER,        // Nested disection
    TACS_AMD_ORDER,  // Interface variables ordered last
    MULTICOLOR_ORDER
  };  // Multicolor via greedy algorithm
  enum MatrixOrderingType {
    ADDITIVE_SCHWARZ,
    APPROXIMATE_SCHUR,
    DIRECT_SCHUR,
    GAUSS_SEIDEL
  };

  // Create the TACSAssembler object in parallel
  // -------------------------------------------
  TACSAssembler(MPI_Comm _tacs_comm, int _varsPerNode, int _numOwnedNodes,
                int _numElements, int _numDependentNodes = 0);
  ~TACSAssembler();

  // Set the connectivity in TACS
  // ----------------------------
  int setElementConnectivity(const int *ptr, const int *conn);
  void getElementConnectivity(const int **ptr, const int **conn);
  int setElements(TACSElement **_elements);
  int setDependentNodes(const int *_depNodeIndex, const int *_depNodeToTacs,
                        const double *_depNodeWeights);

  // Set additional information about the design vector
  // --------------------------------------------------
  void setDesignNodeMap(int _designVarsPerNode,
                        TACSNodeMap *_designVarMap = NULL);
  int setDesignDependentNodes(int numDepDesignVars, const int *_depNodePtr,
                              const int *_depNodes,
                              const double *_depNodeWeights);

  // Associate a Dirichlet boundary condition with the given variables
  // -----------------------------------------------------------------
  void addBCs(int nnodes, const int *nodes, int nbcs = -1,
              const int *vars = NULL, const TacsScalar *vals = NULL);
  void addInitBCs(int nnodes, const int *nodes, int nbcs = -1,
                  const int *vars = NULL, const TacsScalar *vals = NULL);

  // Set Dirichlet BC values at nodes where BCs are imposed
  // ------------------------------------------------------
  void setBCValuesFromVec(TACSBVec *vec);

  // Reorder the unknowns according to the specified reordering
  // ----------------------------------------------------------
  void computeReordering(OrderingType order_type, MatrixOrderingType mat_type);

  // Functions for retrieving the reordering
  // ---------------------------------------
  int isReordered();
  void getReordering(int *oldToNew);
  void reorderVec(TACSBVec *vec);
  void reorderNodes(int num_nodes, int *nodes);

  // Initialize the mesh
  // -------------------
  int initialize();

  // Return important information about the TACSAssembler object
  // -----------------------------------------------------------
  MPI_Comm getMPIComm();
  TACSThreadInfo *getThreadInfo();
  int getVarsPerNode();
  int getDesignVarsPerNode();
  int getNumNodes();
  int getNumDependentNodes();
  int getNumOwnedNodes();
  int getNumElements();
  TACSNodeMap *getNodeMap();
  TACSNodeMap *getDesignNodeMap();
  TACSBcMap *getBcMap();
  TACSBcMap *getInitBcMap();
  TACSBVecDistribute *getBVecDistribute();
  TACSBVecDepNodes *getBVecDepNodes();

  // Get the maximum sizes
  // ---------------------
  int getMaxElementNodes();
  int getMaxElementVariables();
  int getMaxElementDesignVars();

  // Set auxiliary elements into the TACSAssembler object
  // ----------------------------------------------------
  void setAuxElements(TACSAuxElements *aux_elems);
  TACSAuxElements *getAuxElements();

  // Set the nodes in TACS
  // ---------------------
  TACSBVec *createNodeVec();
  void setNodes(TACSBVec *X);
  void getNodes(TACSBVec *X);
  void getNodes(TACSBVec **X);

  // Check for the elements for non-positive determinants
  // ----------------------------------------------------
  void checkElementDeterminants();

  // Set/get the simulation time
  // ---------------------------
  void setSimulationTime(double _time);
  double getSimulationTime();

  // Create vectors
  // --------------
  TACSBVec *createVec();

  // Shortcut to apply boundary conditions
  void applyBCs(TACSVec *vec);
  void applyBCs(TACSMat *mat);
  void applyTransposeBCs(TACSMat *mat);

  // Set the Dirichlet boundary conditions to the state vector
  void setBCs(TACSVec *vec);

  // Methods for manipulating internal variable values
  // -------------------------------------------------
  void zeroVariables();
  void zeroDotVariables();
  void zeroDDotVariables();

  // Methods for setting/getting variables
  // -------------------------------------
  void setVariables(TACSBVec *q, TACSBVec *qdot = NULL, TACSBVec *qddot = NULL);
  void getVariables(TACSBVec *q, TACSBVec *qdot = NULL, TACSBVec *qddot = NULL);
  void getVariables(TACSBVec **q, TACSBVec **qdot = NULL,
                    TACSBVec **qddot = NULL);
  void copyVariables(TACSBVec *q, TACSBVec *qdot = NULL,
                     TACSBVec *qddot = NULL);

  // Create the matrices that can be used for analysis
  // -------------------------------------------------
  TACSParallelMat *createMat();
  TACSSchurMat *createSchurMat(OrderingType order_type = TACS_AMD_ORDER);
  TACSSerialPivotMat *createSerialMat();

  // Retrieve or set the initial conditions for the simulation
  // --------------------------------------------------
  void getInitConditions(TACSBVec *vars, TACSBVec *dvars, TACSBVec *ddvars);
  void setInitConditions(TACSBVec *vars, TACSBVec *dvars, TACSBVec *ddvars);

  // Evaluate the kinetic and potential energy
  // -----------------------------------------
  void evalEnergies(TacsScalar *Te, TacsScalar *Pe);

  // Residual and Jacobian assembly
  // ------------------------------
  void assembleRes(TACSBVec *residual);
  void assembleJacobian(TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                        TACSBVec *residual, TACSMat *A,
                        MatrixOrientation matOr = TACS_MAT_NORMAL);
  void assembleMatType(ElementMatrixType matType, TACSMat *A,
                       MatrixOrientation matOr = TACS_MAT_NORMAL);
  void assembleMatCombo(ElementMatrixType matTypes[], TacsScalar scale[],
                        int nmats, TACSMat *A,
                        MatrixOrientation matOr = TACS_MAT_NORMAL);
  void addJacobianVecProduct(TacsScalar scale, TacsScalar alpha,
                             TacsScalar beta, TacsScalar gamma, TACSBVec *x,
                             TACSBVec *y,
                             MatrixOrientation matOr = TACS_MAT_NORMAL);

  // Assemble data for and compute matrix-free matrix-vector products
  // ----------------------------------------------------------------
  void getMatrixFreeDataSize(ElementMatrixType matType, int *_data_size,
                             int *_temp_size);
  void assembleMatrixFreeData(ElementMatrixType matType, TacsScalar alpha,
                              TacsScalar beta, TacsScalar gamma,
                              TacsScalar data[]);
  void addMatrixFreeVecProduct(ElementMatrixType matType,
                               const TacsScalar data[], TacsScalar temp[],
                               TACSBVec *x, TACSBVec *y,
                               MatrixOrientation matOr = TACS_MAT_NORMAL);

  // Design variable handling
  // ------------------------
  TACSBVec *createDesignVec();
  void getDesignVars(TACSBVec *dvs);
  void setDesignVars(TACSBVec *dvs);
  void getDesignVarRange(TACSBVec *lb, TACSBVec *ub);

  // Function and sensitivity evaluation
  // -----------------------------------
  void evalFunctions(int numFuncs, TACSFunction **funcs, TacsScalar *funcVals);

  // Steady or unsteady derivative evaluation
  // ----------------------------------------
  void addDVSens(TacsScalar coef, int numFuncs, TACSFunction **funcs,
                 TACSBVec **dfdx);
  void addSVSens(TacsScalar alpha, TacsScalar beta, TacsScalar gamma,
                 int numFuncs, TACSFunction **funcs, TACSBVec **dfdu);
  void addAdjointResProducts(TacsScalar scale, int numAdjoints,
                             TACSBVec **adjoint, TACSBVec **dfdx);
  void addXptSens(TacsScalar coef, int numFuncs, TACSFunction **funcs,
                  TACSBVec **dfdXpts);
  void addAdjointResXptSensProducts(TacsScalar scale, int numAdjoints,
                                    TACSBVec **adjoint, TACSBVec **dfdXpts);

  // Advanced function interface - for time integration
  // --------------------------------------------------
  void integrateFunctions(TacsScalar tcoef, TACSFunction::EvaluationType ftype,
                          int numFuncs, TACSFunction **funcs);

  // Add the derivatives of inner products
  // -------------------------------------
  void addMatDVSensInnerProduct(TacsScalar scale, ElementMatrixType matType,
                                TACSBVec *psi, TACSBVec *phi, TACSBVec *dfdx);
  void addMatXptSensInnerProduct(TacsScalar scale, ElementMatrixType matType,
                                 TACSBVec *psi, TACSBVec *phi,
                                 TACSBVec *dfdXpts);
  void evalMatSVSensInnerProduct(ElementMatrixType matType, TACSBVec *psi,
                                 TACSBVec *phi, TACSBVec *res);

  // Return elements and node numbers
  // --------------------------------
  TACSElement **getElements();
  TACSElement *getElement(int elem, TacsScalar *Xpts = NULL,
                          TacsScalar *vars = NULL, TacsScalar *dvars = NULL,
                          TacsScalar *ddvars = NULL);
  TACSElement *getElement(int elem, int *len, const int **nodes);

  // Test the given element, constitutive or function class
  // ------------------------------------------------------
  void testElement(int elemNum, int print_level, double dh = 1e-6,
                   double rtol = 1e-8, double atol = 1e-1);
  void testFunction(TACSFunction *func, double dh);

  // Set the number of threads to work with
  // --------------------------------------
  void setNumThreads(int t);

  // Get information about the output files; For use by TACSToFH5
  // ------------------------------------------------------------
  int getNumComponents();
  void getElementOutputData(ElementType elem_type, int write_flag, int *len,
                            int *nvals, TacsScalar **data);

  // Functions for ordering the variables
  // ------------------------------------
  int getLocalNodeNum(int node);
  int getGlobalNodeNum(int node);
  void computeLocalNodeToNodeCSR(int **_rowp, int **_cols, int nodiag = 0);
  void computeNodeToElementCSR(int **_nodeElem, int **_nodeElemIndex);

 private:
  // Get the number of design variable numbers
  // -----------------------------------------
  int getNumDesignVars();

  // Get pointers to the start-locations within the data array
  // ---------------------------------------------------------
  void getDataPointers(TacsScalar *data, TacsScalar **v1, TacsScalar **v2,
                       TacsScalar **v3, TacsScalar **v4, TacsScalar **x1,
                       TacsScalar **x2, TacsScalar **weights, TacsScalar **mat);

  // Functions that are used to perform reordering
  // ---------------------------------------------
  int computeExtNodes();
  int computeCouplingNodes(int **_couplingNodes, int **_extPtr = NULL,
                           int **_extCount = NULL, int **_recvPtr = NULL,
                           int **_recvCount = NULL, int **_recvNodes = NULL);
  int computeCouplingElements(int **_celems);

  // Functions for ordering the variables
  // ------------------------------------
  void computeLocalNodeToNodeCSR(int **_rowp, int **_cols, int nrnodes,
                                 const int *rnodes, int nodiag);

  // Compute the connectivity of the multiplier information
  void computeMultiplierConn(int *_num_multipliers, int **_multipliers,
                             int **_indep_ptr, int **_indep_nodes);

  // Compute the reordering for a local matrix
  // -----------------------------------------
  void computeMatReordering(OrderingType order_type, int nvars, int *rowp,
                            int *cols, int *perm, int *new_vars);

  // Scatter the boundary conditions on external nodes
  void scatterExternalBCs(TACSBcMap *bcs);

  // Add values into the matrix
  inline void addMatValues(TACSMat *A, const int elemNum, const TacsScalar *mat,
                           int *item, TacsScalar *temp,
                           MatrixOrientation matOr);

  TACSNodeMap *nodeMap;               // Variable ownership map
  TACSBcMap *bcMap;                   // Boundary condition data
  TACSBcMap *bcInitMap;               // Initial boundary condition data
  TACSBVecDistribute *extDist;        // Distribute the vector
  TACSBVecIndices *extDistIndices;    // The tacsVarNum indices
  TACSBVecDepNodes *depNodes;         // Dependent variable information
  TACSNodeMap *designNodeMap;         // Distribution of design variables
  TACSBVecDistribute *designExtDist;  // Distribute the design variables
  TACSBVecDepNodes *designDepNodes;   // Dependent design variable information

  // Reordering information
  TACSBVecIndices *newNodeIndices;

  // Additional information information for the TACSParallel class
  TACSBVecIndices *parMatIndices;

  // Additional ordering information for the TACSSchurMat class
  // These are created once - all subsequent calls use this data.
  TACSBVecIndices *schurBIndices, *schurCIndices;
  TACSBVecDistribute *schurBMap, *schurCMap;

  // The global simulation time variable
  double time;

  // variables/elements have been initialized
  int meshInitializedFlag;

  // Information about the variables and elements
  int varsPerNode;         // number of variables per node
  int numElements;         // number of elements
  int numNodes;            // number of nodes referenced by this process
  int numOwnedNodes;       // number of nodes owned by this processor
  int numExtNodes;         // number of extneral nodes
  int numDependentNodes;   // number of dependent nodes
  int numMultiplierNodes;  // number of multiplier nodes/elements
  int designVarsPerNode;   // number of design variables at each design "node"

  // Maximum element information
  int maxElementDesignVars;  // maximum number of design variable
  int maxElementNodes;       // maximum number of ind. and dep. element nodes
  int maxElementSize;        // maximum number of variables for any element
  int maxElementIndepNodes;  // maximum number of independent nodes

  // Node numbers that are referred to from this processor
  int *tacsExtNodeNums;  // node numbers associated with TACS
  int extNodeOffset;     // Offset into the external nodes

  // Variables that define the CSR data structure to
  // store the element -> node information
  int *elementNodeIndex;
  int *elementTacsNodes;

  // The local list of elements
  TACSElement **elements;

  // The auxiliary element class
  TACSAuxElements *auxElements;

  // The variables, velocities and accelerations
  TACSBVec *varsVec, *dvarsVec, *ddvarsVec;

  // Memory for the node locations
  TACSBVec *xptVec;

  // Memory for the element residuals and variables
  TacsScalar *elementData;  // Space for element residuals/matrices
  int *elementIData;        // Space for element index data

  // Memory for the design variables and inddex data
  TacsScalar *elementSensData;
  int *elementSensIData;

  // Memory for the initial condition vectors
  TACSBVec *vars0, *dvars0, *ddvars0;

  // The data required to perform parallel operations
  // MPI info
  int mpiRank, mpiSize;
  MPI_Comm tacs_comm;

  // The static member functions that are used to p-thread TACSAssembler
  // operations... These are the most time-consuming operations.
  static void schedPthreadJob(TACSAssembler *tacs, int *index, int total_size);
  static void *assembleRes_thread(void *t);
  static void *assembleJacobian_thread(void *t);
  static void *assembleMatType_thread(void *t);

  // Class to store specific information about the threaded
  // operations to perform. Note that assembly operations are
  // relatively easy, while design-variable dependent info is
  // much more challenging!
  class TACSAssemblerPthreadInfo {
   public:
    TACSAssemblerPthreadInfo() {
      assembler = NULL;
      res = NULL;
      mat = NULL;
      alpha = beta = gamma = 0.0;
      matType = TACS_STIFFNESS_MATRIX;
      matOr = TACS_MAT_NORMAL;
      coef = 0.0;
      numFuncs = 0;
      functions = NULL;
      ftype = TACSFunction::INTEGRATE;
      numDesignVars = 0;
      numAdjoints = 0;
      fdvSens = NULL;
      fXptSens = NULL;
      adjoints = NULL;
    }

    // The data required to perform most of the matrix
    // assembly.
    TACSAssembler *assembler;

    // Information for residual assembly
    TACSBVec *res;

    // Information for matrix assembly
    TACSMat *mat;
    TacsScalar alpha, beta, gamma;
    ElementMatrixType matType;
    MatrixOrientation matOr;

    // Information required for the computation of f or df/dx
    double coef;
    int numFuncs;
    TACSFunction **functions;
    TACSFunction::EvaluationType ftype;

    int numDesignVars;
    TacsScalar *fdvSens;  // df/dx
    TACSBVec **fXptSens;

    // Information for adjoint-dR/dx products
    int numAdjoints;
    TACSBVec **adjoints;
  } * tacsPInfo;

  // The pthread data required to pthread tacs operations
  int numCompletedElements;     // Keep track of how much work has been done
  TACSThreadInfo *thread_info;  // The pthread object

  // The thread objects
  pthread_t threads[TACSThreadInfo::TACS_MAX_NUM_THREADS];
  pthread_mutex_t tacs_mutex;  // The mutex for coordinating assembly ops.

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
inline void TACSAssembler::addMatValues(TACSMat *A, const int elemNum,
                                        const TacsScalar *mat, int *itemp,
                                        TacsScalar *temp,
                                        MatrixOrientation matOr) {
  int start = elementNodeIndex[elemNum];
  int end = elementNodeIndex[elemNum + 1];
  int nnodes = end - start;
  int nvars = varsPerNode * nnodes;

  // Add the element values to the matrix
  const int *nodeNums = &elementTacsNodes[start];

  if (matOr == TACS_MAT_NORMAL && numDependentNodes == 0) {
    // If we have no dependent nodes, then we don't need to do
    // anything extra here
    A->addValues(nnodes, nodeNums, nnodes, nodeNums, nvars, nvars, mat);
  } else {
    // If we have dependent nodes, then we have to figure out what
    // the weighting matrix is and add then add the element matrix
    const int *depNodePtr = NULL;
    const int *depNodeConn = NULL;
    const double *depNodeWeights = NULL;
    if (depNodes) {
      depNodes->getDepNodes(&depNodePtr, &depNodeConn, &depNodeWeights);
    }

    // Set pointers to the temporary arrays
    int *varp = &itemp[0];
    int *vars = &itemp[nnodes + 1];
    TacsScalar *weights = temp;

    varp[0] = 0;
    for (int i = 0, k = 0; i < nnodes; i++) {
      if (nodeNums[i] >= 0) {
        // This is just a regular node
        weights[k] = 1.0;
        vars[k] = nodeNums[i];
        k++;
      } else {
        // This is a dependent node. Determine the corresponding
        // dependent node number and add the variables
        int dep = -nodeNums[i] - 1;
        for (int j = depNodePtr[dep]; j < depNodePtr[dep + 1]; j++, k++) {
          weights[k] = depNodeWeights[j];
          vars[k] = depNodeConn[j];
        }
      }

      varp[i + 1] = k;
    }

    // Add the values to the matrix
    A->addWeightValues(nnodes, varp, vars, weights, nvars, nvars, mat, matOr);
  }
}

#endif  // TACS_ASSEMBLER_H
