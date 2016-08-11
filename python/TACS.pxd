# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
from libc.string cimport const_char

# Typdef 
ctypedef double TacsScalar

cdef extern from "TACSElement.h":
   enum ElementType:
       ELEMENT_NONE
       POINT_ELEMENT
       EULER_BEAM
       TIMOSHENKO_BEAM
       PLANE_STRESS
       SHELL
       SOLID
       Q3D_ELEMENT
       RIGID

   enum ElementMatrixType:
       STIFFNESS_MATRIX
       MASS_MATRIX
       GEOMETRIC_STIFFNESS_MATRIX
       STIFFNESS_PRODUCT_DERIVATIVE

cdef extern from "KSM.h":
   enum MatrixOrientation:
       NORMAL
       TRANSPOSE

# Special functions required for converting pointers
cdef extern from "":
   ScMat* _dynamicScMat "dynamic_cast<ScMat*>"(TACSMat*)
   PMat* _dynamicPMat "dynamic_cast<PMat*>"(TACSMat*)

cdef extern from "TACSObject.h":
   cdef cppclass TACSObject:
      void incref()
      void decref()

cdef extern from "KSM.h":
   cdef cppclass TACSVec(TACSObject):
      TacsScalar norm()
      void scale(TacsScalar alpha)
      TacsScalar dot(TACSVec *x)
      void axpy(TacsScalar alpha, TACSVec * x)
      void copyValues(TACSVec * x)
      void axpby(TacsScalar alpha, TacsScalar beta, TACSVec * x)
      void zeroEntries()
      void setRand(double lower, double upper)
      void applyBCs()
      
   cdef cppclass TACSMat(TACSObject):
      TACSVec *createVec()
      void zeroEntries()
      void applyBCs()
      void mult(TACSVec * x, TACSVec * y)
      void copyValues(TACSMat * mat)
      void scale(TacsScalar alpha)
      void axpy(TacsScalar alpha, TACSMat * mat)
        
   cdef cppclass TACSPc(TACSObject):
      void factor()
      void applyFactor(TACSVec *x, TACSVec *y)
        
   cdef cppclass KSMPrint(TACSObject):
      pass

   cdef cppclass KSMPrintStdout(KSMPrint):
      KSMPrintStdout(char *descript, int rank, int freq)

   cdef cppclass TACSKsm(TACSObject):
      TACSVec *createVec()
      void setOperators(TACSMat * _mat, TACSPc * _pc)
      void getOperators(TACSMat ** _mat, TACSPc ** _pc)
      void solve(TACSVec * b, TACSVec * x, int zero_guess)
      void setTolerances(double _rtol, double _atol)
      void setMonitor(KSMPrint * _monitor)

   cdef cppclass GMRES(TACSKsm):
      GMRES(TACSMat * _mat, TACSPc * _pc, int _m,
            int _nrestart, int _isFlexible )
      
cdef extern from "BVec.h":
   cdef cppclass VarMap(TACSObject):
      VarMap(MPI_Comm comm, int N, int bsize)

   cdef cppclass BCMap(TACSObject):
      BCMap(int num_bcs)
        
   cdef cppclass BVec(TACSVec):
      BVec(VarMap* rmap, BCMap* bcs)
      int getSize(int *size) 
      int getArray(TacsScalar ** array)
      void placeArray(TacsScalar * _x)
      int readFromFile(const_char *filename)
      int writeToFile(const_char *filename)
        
cdef extern from "PMat.h":
   cdef cppclass PMat(TACSMat):
      pass

   cdef cppclass AdditiveSchwarz(TACSPc):
      AdditiveSchwarz(PMat * mat, int levFill, double fill)

   cdef cppclass ApproximateSchur(TACSPc):
      ApproximateSchur(PMat * mat, int levFill, double fill, 
                       int inner_gmres_iters, double inner_rtol,
                       double inner_atol)
      
cdef extern from "DistMat.h":
   cdef cppclass DistMat(PMat):
      pass

cdef extern from "ScMat.h":
   cdef cppclass ScMat(TACSMat):
      pass

   cdef cppclass PcScMat(TACSPc):
      PcScMat(ScMat *mat, int levFill, double fill, 
              int reorder_schur_complement)
    
cdef extern from "FEMat.h":
   cdef cppclass FEMat(ScMat):
      pass
       
cdef extern from "TACSElement.h":
   cdef cppclass TACSElement(TACSObject):
      pass

cdef extern from "TACSFunction.h":
   cdef cppclass TACSFunction(TACSObject):
      pass

cdef extern from "TACSConstitutive.h":
   cdef cppclass TACSConstitutive(TACSObject):
      pass
    
# A generic wrapper class for the TACSElement object
cdef class Element:
   cdef TACSElement *ptr
   
# A generic wrapper class for the TACSFunction object
cdef class Function:
   cdef TACSFunction *ptr
   
# A generic wrapper class for the TACSConstitutive object
cdef class Constitutive:
   cdef TACSConstitutive *ptr

cdef class Vec:
   cdef BVec *ptr
   
cdef extern from "TACSAuxElements.h":
   cdef cppclass TACSAuxElements(TACSObject):
      TACSAuxElements(int)
      void addElement(int, TACSElement*)
      
cdef extern from "TACSAssembler.h":
   enum OrderingType"TACSAssembler::OrderingType":
      NATURAL_ORDER"TACSAssembler::NATURAL_ORDER"
      RCM_ORDER"TACSAssembler::RCM_ORDER"
      AMD_ORDER"TACSAssembler::AMD_ORDER"
      ND_ORDER"TACSAssembler::ND_ORDER"
      TACS_AMD_ORDER"TACSAssembler::TACS_AMD_ORDER"
      
   enum MatrixOrderingType"TACSAssembler::MatrixOrderingType":
      ADDITIVE_SCHWARZ"TACSAssembler::ADDITIVE_SCHWARZ"
      APPROXIMATE_SCHUR"TACSAssembler::APPROXIMATE_SCHUR"
      DIRECT_SCHUR"TACSAssembler::DIRECT_SCHUR"

   cdef cppclass TACSAssembler(TACSObject):
      TACSAssembler(MPI_Comm _tacs_comm, int numOwnedNodes,
                    int _varsPerNode, int _numElements, 
                    int _numNodes, int _nodeMaxCSRsize)
      TACSAssembler(MPI_Comm _tacs_comm, int numOwnedNodes,
                    int _varsPerNode, int _numElements, 
                    int _numNodes, int numDependentNodes,
                    int _nodeMaxCSRsize)
      # Return information about the TACSObject
      int getNumNodes()
      int getNumDependentNodes()
      int getNumElements()
      BVec *createNodeVec()
      void setNodes(BVec* x)
      void getNodes(BVec* x)
      void getDesignVars(TacsScalar* dvs, int ndvs)
      void setDesignVars(TacsScalar* dvs, int ndvs)
      void getDesignVarRange(TacsScalar* lb, TacsScalar* ub,
                             int ndvs)
      void setNumThreads(int t)

      # Set the auxiliary element class
      void setAuxElements(TACSAuxElements*)
      TACSAuxElements *getAuxElements()
      
      # Create vectors/matrices
      BVec* createVec()
      DistMat* createMat()
      FEMat* createFEMat(OrderingType order_type)

      # Set/get the simulation time
      void setSimulationTime(double)
      double getSimulationTime()      
      
      # Zero the variables
      void zeroVariables()
      void zeroDotVariables()
      void zeroDDotVariables()

      # Set and retrieve the state variables
      void setVariables(BVec* stateVars)
      void setDotVariables(BVec* stateVars)
      void setDDotVariables(BVec* stateVars)
      void getVariables(BVec* stateVars)

      # Add nodes to TACS
      void addNode(int locaNodeNum, int tacsNodeNum)
      void addNodes(int localNodeNum[], int tacsNodeNum[], int numNodes)

      # Add elements to TACS
      int addElement(TACSElement *element, int localNodes[], int numElemNodes)

      # Assembly routines
      void assembleRes(BVec* residual)
      void assembleJacobian(BVec* residual,
                            TACSMat* A, double alpha,
                            double beta, double gamma,
                            MatrixOrientation matOr)
      void assembleMatType(ElementMatrixType matType, 
                           TACSMat* A, MatrixOrientation matOr) 

      # Evaluation routines
      void evalFunctions(TACSFunction **functions, int numFuncs,
                         TacsScalar *funcVals)
      void evalDVSens(TACSFunction **funcs, int numFuncs,
                      TacsScalar *fdvsSens, int ndvs)
      void evalSVSens(TACSFunction *function, BVec* vec)
      void evalAdjointResProducts(BVec** adjoint, int numAdjoint,
                                  TacsScalar* dvSens, int num_dvs)

      # Test routines
      void testElement(int elemNum, int print_level)
      void testConstitutive(int elemNum, int print_level)
      void testFunction(TACSFunction *func, int num_dvs,
                        double dh)


      # MPI communicator
      MPI_Comm getMPIComm()

      # Finalize the mesh - no further elements or nodes may be added
      # following this call
      void finalize()

cdef class Assembler:
   cdef TACSAssembler *ptr
      
cdef extern from "TACSMeshLoader.h":
   cdef cppclass TACSMeshLoader(TACSObject):
      TACSMeshLoader(MPI_Comm _comm)
      int scanBDFFile(char * file_name)
      int getNumComponents()
      char* getComponentDescript(int comp_num)
      char* getElementDescript(int comp_num)
      void setElement(int component_num, TACSElement * _element)
      int getNumNodes()
      int getNumElements()
      TACSAssembler* createTACS(int vars_per_node,
                                OrderingType order_type,
                                MatrixOrderingType mat_type)
      void getConnectivity(int *_num_nodes, int *_num_elements,
                           int **_elem_node_ptr, 
                           int **_elem_node_conn,
                           double **_Xpts)

cdef extern from "TACSCreator.h":
   cdef cppclass TACSCreator(TACSObject):
      TACSCreator(MPI_Comm comm, int _vars_per_node)
      void setGlobalConnectivity(int _num_nodes, int _num_elements,
                                 int *_elem_node_ptr, 
                                 int *_elem_node_conn,
                                 int *_elem_id_nums )
      void setBoundaryConditions(int _num_bcs, int *_bc_nodes, 
                                 int *_bc_vars, int *_bc_ptr )
      void setDependentNodes(int num_dep_nodes, 
                             int *_dep_node_ptr,
                             int *_dep_node_conn, 
                             double *_dep_node_weights )
      void setElements(TACSElement **_elements, int _num_elems)
      void setNodes(TacsScalar *_Xpts)
      void setReorderingType(OrderingType _order_type,
                             MatrixOrderingType _mat_type)
      void partitionMesh(int split_size, int *part)
      int getElementPartition(const int **)
      TACSAssembler *createTACS()
      
cdef extern from "TACSToFH5.h":
   cdef cppclass TACSToFH5(TACSObject):
      TACSToFH5(TACSAssembler * _tacs, ElementType _elem_type, int _out_type)
      void setComponentName(int comp_num, char *group_name)
      void writeToFile(char *filename)

cdef extern from "TACSIntegrator.h":
   # Declare the TACSIntegrator base class
   cdef cppclass TACSIntegrator(TACSObject):      
      # Integrate forward in time
      void integrate()

      # Returns the adjoint gradient for all functions that are set into TACS
      void getFuncGrad(int num_dv, TacsScalar *x, TacsScalar *fvals,
                       TacsScalar *dfdx)

      # Returns the finite-difference gradient for all the functions that are set into TACS
      void getFDFuncGrad(int num_dv, TacsScalar *x, TacsScalar *fvals,
                         TacsScalar *dfdx, double dh)

      # Setters for class variables
      void setFunction(TACSFunction **func, int num_funcs)
      void setPrintLevel(int print_level, char *filename)
      void setRelTol(double rtol)
      void setAbsTol(double atol)
      void setMaxNewtonIters(int max_newton_iters)
      void setJacAssemblyFreq(int jac_comp_freq)
      void setUseLapack(int use_lapack)

      # Configure writing F5 files
      void configureOutput(TACSToFH5 *viewer, int write_freq, char *f5_file_fmt)

   # BDF Implementation of the integrator
   cdef cppclass TACSBDFIntegrator(TACSIntegrator):
      TACSBDFIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec,
                        int max_bdf_order)

   # DIRK Implementation of the integrator
   cdef cppclass TACSDIRKIntegrator(TACSIntegrator):
      TACSDIRKIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec,
                        int num_stages)
      
   # ABM Implementation of the integrator
   cdef cppclass TACSABMIntegrator(TACSIntegrator):
      TACSABMIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec,
                        int max_abm_order)

   # NBG Implementation of the integrator
   cdef cppclass TACSNBGIntegrator(TACSIntegrator):
      TACSNBGIntegrator(TACSAssembler *tacs,
                        double tinit, double tfinal,
                        int num_steps_per_sec)
