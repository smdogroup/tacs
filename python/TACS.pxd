# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

cdef enum OrderingType:
    NO_ORDER, RCM_ORDER, AMD_ORDER, ND_ORDER, TACS_AMD_ORDER

cdef enum MatrixOrderingType:
    ADDITIVE_SCHWARZ, APPROXIMATE_SCHUR, DIRECT_SHUR

cdef enum MatrixOrientation:
    NORMAL, TRANSPOSE

cdef extern from "TACSObject.h":
    cdef cppclass TACSObject:
       TACSObject()

cdef extern from "TACSAssembler.h":
    cdef cppclass TACSAssembler:
       TACSAssembler(MPI_Comm _tacs_comm, int numOwnedNodes,
                     int _varsPerNode, int _numElements, 
                     int _numNodes, int numDependentNodes,
                     int _nodeMaxCSRsize)
       # Return information about the TACSObject
       int getNumNodes()
       int getNumDependentNodes()
       int getNumElements()
       #VarMap getVarMap(void *_self)
       
       # Add nodes to TACS
       void addNode(int localNodeNum, int tacsNodeNum)
       void addNodes(int* localNodeNums, int *tacsNodeNums1, 
                     int num_nodes)
       void setDependentNodes(int* depNodeIndex, 
                              int* depNodeToLocal,
                              TacsScalar* depNodeWeights)
       void addBC(int nodeNum, int* bcNums, 
                  TacsScalar* bcVals, int nbcs)
       void addElement(TACSElement *element, int* localNodeNums,
                       int num_nodes)

       void setElementConnectivity()

       void setNodes(BVec* x)
       void getNodes(BVec* x)

       void computeLocalNodeToNodeCSR()
       void computeCouplingNodes()
       void computeCouplingElements()
       void computeReordering()
       void computeMatReordering()
       void computeNodeToElementCSR()

       void finalize()
       void getNumDesignVars()
       void getDesignVars(TacsScalar* dvs, int ndvs)
       void setDesignVars(TacsScalar* dvs, int ndvs)
       void getDesignVarRange(TacsScalar* lb, TacsScalar* ub,
                              int ndvs)
       void setNumThreads(int t)
       
       # Create vectors/matrices
       BVec* createVec()
       DistMat* createMat()
       FEMat* createFEMat(OrderingType order_type = TACS_AMD_ORDER)

       # Zero the variables
       void zeroVariables()
       void zeroDotVariables()
       void zeroDDotVariables()

       # Set and retrieve the state variables
       void setVariables(BVec* stateVars)
       void setDotVariables(BVec* stateVars)
       void setDDotVariables(BVec* stateVars)
       void getVariables(BVec* stateVars)

       void getTacsNodeNums(int* localNodes, int num_nodes)
       void assembleRes(BVec* residual)
       void assembleResNoBCs(BVec* residual)
       void assembleJacobian(BVec* residual,
                             TACSMat* A, double alpha,
                             double beta, double gamma,
                             MatrixOrientation matOr)
       void assembleMatType(ElementMatrixType matType, 
                            TACSMat* A, MatrixOrientation matOr) 
       void evalFunctions(TACSFunction **functions, int numFuncs,
                          TacsScalar *funcVals)
       void evalDVSens(TACSFunction **funcs, int numFuncs,
                       TacsScalar *fdvsSens, int ndvs)
       void evalXptSens(TACSFunction **funcs, int numFuncs, 
                        TacsScalar* fXptSens)
       void evalSVSens(TACSFunction **function, BVec* vec)
       void evalAdjointResProducts(BVec* adjoint, int numAdjoint,
                                   TacsScalar* dvSens, int num_dvs)
       # void evalAdjointResProductsExperimental(BVec* adjoint, 
       #                                         int numAdjoint,
       #                                         TacsScalar* dvSens, 
       #                                         int num_dvs)
       void evalMatDVSensInnerProduct(TacsScalar alpha,
                                      ElementMatrixTypes matType,
                                      BVec* psi, BVec *phi,
                                      TacsScalar* dvSens)
       void evalMatSVSensInnerProduct(TacsScalar alpha, 
                                      ElementMatrixTypes matType,
                                      BVec* psi, BVec* phi,
                                      BVec* res)
       int getElements()
       int getElement(int elemNum, TacsScalar* elemVars,
                      TacsScalar elemXpts)
       void testElement(int elemNum, int print_level)
       void testConstitutive(int elemNum, int print_level)
       void testFunction(TACSFunction *func, int num_dvs,
                         double dh)
       int getNumComponents()
       void getOutputNodeRange(ElementType elem_type, 
                               int *node_range)
       void getOutputConnectivity(ElementType elem_type)
       void getOutputData(ElementType elem_type, int out_type,
                          double *Data, int nvals)
       MPI_Comm getMPIComm()

cdef extern from "CyTACSVec.h":
    cdef cppclass CyTACSVec:
        CyTACSVec()
        void setSelfPointer(void *_self)
    
cdef extern from "CyTACSMat.h":
    cdef cppclass CyTACSMat:
        CyTACSMat()
        void setSelfPointer(void *_self)

cdef extern from "BVec.h":
    cdef cppclass VarMap:
        VarMap(MPI_Comm comm, int N, int bsize)
    cdef cppclass BCMap:
        BCMap(int num_bcs)
    cdef cppclass BVec:
        BVec(VarMap* rmap, BCMap* bcs)
        #int getSize(int *size) 
        # Basic operations for Krylov methods
        TacsScalar norm()
        void scale( TacsScalar alpha)
        TacsScalar dot(BVec* x)
        void axpy( TacsScalar alpha, TACSVec *y )
        void copyValues(TACSVec * x)
        void zeroEntries()
        void setRand(double lower, double upper)
        void applyBCs()
        
        void set( TacsScalar val)
        void restoreArray()
        void initRand()
        void placeArray(TacsScalar *x)
        void applyBCs()
        int readFromFile(const_char *filename)
        int writeToFile(const_char *filename)
cdef extern from "CyTACSMat.h":
    cdef cppclass CyTACSMat:
        CyTACSMat()
        void setSelfPointer(void *_self)

cdef extern from "PMat.h":
    cdef cppclass PMat:
        PMat()
        void zeroEntries()
        void applyBCs()
        void mult(TACSVec *x, TACSVec *y)
        void copyValues(PMat *mat)
        void scale(TacsScalar alpha)
        void axpy(TacsScalar alpha, TACSMat *mat)
        void axpby(TacsScalar alpha, TacsScalar beta,
                   PMat *mat)
cdef extern from "DistMat.h":
    cdef cppclass DistMat:
        DistMat()
        void zeroEntries()
        void applyBCs()
        void mult(TACSVec *x, TACSVec *y)
