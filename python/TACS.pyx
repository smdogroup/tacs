# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy 
cimport numpy as np
import numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char
from libc.stdlib cimport malloc, free

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import the definitions
from TACS cimport *

cdef extern from "mpi-compat.h":
   pass

# Import the element types
PY_ELEMENT_NONE = ELEMENT_NONE
PY_POINT_ELEMENT = POINT_ELEMENT
PY_EULER_BEAM = EULER_BEAM
PY_TIMOSHENKO_BEAM = TIMOSHENKO_BEAM
PY_PLANE_STRESS = PLANE_STRESS
PY_SHELL = SHELL 
PY_SOLID = SOLID
PY_Q3D_ELEMENT = Q3D_ELEMENT

# Import the element matrix types
PY_STIFFNESS_MATRIX = STIFFNESS_MATRIX
PY_MASS_MATRIX = MASS_MATRIX
PY_GEOMETRIC_STIFFNESS_MATRIX = GEOMETRIC_STIFFNESS_MATRIX
PY_STIFFNESS_PRODUCT_DERIVATIVE = STIFFNESS_PRODUCT_DERIVATIVE

# Import the orientations
PY_NORMAL = NORMAL
PY_TRANSPOSE = TRANSPOSE
   
# A generic wrapper class for the TACSFunction object
cdef class Function:
   def __cinit__(self):
      self.ptr = NULL
      return

   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return
   
# A generic wrapper for a TACSVec class - usually BVec
cdef _init_Vec(BVec *ptr):
   vec = Vec()
   vec.ptr = ptr
   vec.ptr.incref()
   return vec

cdef class Vec:
   cdef BVec *ptr
   def __cinit__(self):
      '''
      A generic wrapper for any of the TACS vector types
      '''
      self.ptr = NULL
      return

   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return
   
   def zeroEntries(self):
      '''
      Zero the entries in the matrix
      '''
      self.ptr.zeroEntries()
      return
      
   def applyBCs(self):
      '''
      Apply the Dirichlet boundary conditions to the matrix
      '''
      self.ptr.applyBCs()
      return
      
   def norm(self):
      '''
      Vector norm
      '''
      return self.ptr.norm()
     
   def dot(self, Vec vec):
      '''
      Take the dot product with the other vector
      '''
      return self.ptr.dot(vec.ptr)
   
   def copyValues(self, Vec vec):
      '''
      Copy the values from mat         
      '''
      self.ptr.copyValues(vec.ptr)
      return
    
   def scale(self, TacsScalar alpha):
      '''
      Scale the entries in the matrix by alpha
      '''
      self.ptr.scale(alpha)
      return
    
   def axpy(self, TacsScalar alpha, Vec vec):
      '''
      Compute y <- y + alpha * x
      '''
      self.ptr.axpy(alpha, vec.ptr)
      return
    
   def axpby(self, TacsScalar alpha, TacsScalar beta, Vec vec):
      '''
      Compute y <- alpha * x + beta * y
      '''
      self.ptr.axpby(alpha, beta, vec.ptr)
      return

   def setRand(self, double lower=-1.0, double upper=1.0):
      '''
      Set random entries
      '''
      self.ptr.setRand(lower, upper)
      return
      
   def writeToFile(self, char*filename=''):
      '''
      Write the values to a file.
      
      This uses MPI file I/O. The filenames must be the same on all
      processors. The format is independent of the number of processors.
      
      The file format is as follows:
      int                       The length of the vector
      len * sizeof(TacsScalar)  The vector entries
      '''
      return self.ptr.writeToFile(&filename[0])
   
   def readFromFile(self, char*filename=''):
      '''
      Read values from a binary data file.
      
      The size of this vector must be the size of the vector
      originally stored in the file otherwise nothing is read in.
      
      The file format is as follows:
      int                       The length of the vector
      len * sizeof(TacsScalar)  The vector entries
      '''
      return self.ptr.readFromFile(&filename[0])

# Python class for corresponding instance Mat
cdef _init_Mat(TACSMat *ptr):
   mat = Mat()
   mat.ptr = ptr
   mat.ptr.incref()
   return mat

cdef class Mat:
   cdef TACSMat *ptr
   def __cinit__(self):
      '''
      A generic wrapper for any of the TACS matrix types
      '''        
      self.ptr = NULL
      return
      
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return

   def zeroEntries(self):
      '''
      Zero the entries in the matrix
      '''
      self.ptr.zeroEntries()

   def applyBCs(self):
      '''
      Apply the Dirichlet boundary conditions to the matrix
      '''
      self.ptr.applyBCs()

   def mult(self, Vec x, Vec y):
      '''
      Matrix multiplication
      '''
      self.ptr.mult(x.ptr, y.ptr)

   def copyValues(self, Mat mat):
      '''
      Copy the values from mat         
      '''
      self.ptr.copyValues(mat.ptr) 

   def scale(self, TacsScalar alpha):
      '''
      Scale the entries in the matrix by alpha
      '''
      self.ptr.scale(alpha)

   def axpy(self, TacsScalar alpha, Mat mat):
      '''
      Compute y <- y + alpha * x
      '''
      self.ptr.axpy(alpha, mat.ptr)

# Create a generic preconditioner class
cdef class Pc:
   cdef TACSPc *ptr
   def __cinit__(self, Mat mat):
      '''
      This creates a default preconditioner depending on the matrix
      type.
      '''
      cdef PMat *p_ptr = NULL #_dynamicPMat(mat.ptr)
      cdef ScMat *sc_ptr = _dynamicScMat(mat.ptr)

      # Set the defaults for the direct factorization
      cdef int lev_fill = 1000000
      cdef double fill = 10.0
      cdef int reorder = 1
      
      if sc_ptr != NULL:
         self.ptr = new PcScMat(sc_ptr, lev_fill, fill, reorder)
      else:
         self.ptr = new AdditiveSchwarz(p_ptr, 5, 10.0)
                  
      self.ptr.incref()
      return
         
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return

   def factor(self):
      '''Factor the preconditioner'''
      self.ptr.factor()

   def applyFactor(self, Vec x, Vec y):
      '''Apply the preconditioner'''
      self.ptr.applyFactor(x.ptr, y.ptr)

cdef class KSM:
   cdef TACSKsm *ptr
   def __cinit__(self, Mat mat, Pc pc, int m, int nrestart, int isFlexible=0):
      '''
      Create a GMRES object for solving a linear system with or without a preconditioner.
      
      This automatically allocates the requried Krylov subspace on
      initialization.

      input:
      mat:        the matrix operator
      pc:         the preconditioner
      m:          the size of the Krylov subspace
      nrestart:   the number of restarts before we give up
      isFlexible: is the preconditioner actually flexible? If so use FGMRES
      '''
      
      self.ptr = new GMRES(mat.ptr, pc.ptr, m, nrestart, isFlexible)
      self.ptr.incref()
      return

   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return

   def solve(self, Vec b, Vec x, int zero_guess=1):
      '''
      Try to solve the linear system using GMRES.

      The following code tries to solve the linear system using GMRES (or
      FGMRES if the preconditioner is flexible.)

      input:
      b:          the right-hand-side
      x:          the solution vector (with possibly significant entries)
      zero_guess: flag to indicate whether to zero entries of x before solution
      '''
      self.ptr.solve(b.ptr, x.ptr, zero_guess)

   def setTolerances(self, double rtol, double atol):
      '''
      Set the relative and absolute tolerances used for the stopping
      criterion.

      input:
      rtol: the relative tolerance ||r_k|| < rtol*||r_0||
      atol: the absolute tolerancne ||r_k|| < atol
      '''
      self.ptr.setTolerances(rtol, atol)

   def setMonitor(self, MPI.Comm comm, char *descript, int freq):
      '''
      Set the object to control how the convergence history is displayed
      (if at all)

      input:
      monitor: the KSMPrint monitor object
      '''
      self.ptr.setMonitor(new KSMPrintStdout(&descript[0], comm.rank, freq))
        
# Python class for corresponding instance TACSAssembler
cdef _init_Assembler(TACSAssembler *ptr):
   tacs = Assembler()
   tacs.ptr = ptr
   tacs.ptr.incref()
   return tacs

cdef class Assembler:   
   def __cinit__(self):
      '''
      Constructor for the TACSAssembler object
      '''
      self.ptr = NULL
      return
   
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return
   
   def getNumNodes(self):
      '''
      Return the number of nodes in the TACSAssembler
      '''
      return self.ptr.getNumNodes()

   def getNumDependentNodes(self):
      '''
      Return the number of dependent nodes
      '''
      return self.ptr.getNumDependentNodes()

   def getNumElements(self):
      '''
      Return the number of elements
      '''
      return self.ptr.getNumElements()
   
   def getDesignVars(self, 
                     np.ndarray[TacsScalar, ndim=1, mode='c']dvs):
      '''
      Collect all the design variable values assigned by this 
      process
      
      This code does not ensure consistency of the design variable
      values between processes. If the values of the design 
      variables are inconsistent to begin with, the maximum
      design variable value is returned. Call setDesignVars to 
      make them consistent.
    
      Each process contains objects that maintain their own design
      variable values. Ensuring the consistency of the ordering is
      up to the user. Having multiply-defined design variable 
      numbers corresponding to different design variables 
      results in undefined behaviour.
      
      dvs:    the array of design variable values (output)
      numDVs: the number of design variables
      '''
    
      # Get number of the design variables
      num_dvs = dvs.shape[0]
      self.ptr.getDesignVars(<TacsScalar*>dvs.data, num_dvs)
      return
      
   def setDesignVars(self, 
                     np.ndarray[TacsScalar, ndim=1, mode='c']dvs):
      '''
      Set the design variables.
      
      The design variable values provided must be the same on all
      processes for consistency. This call however, is not 
      collective.
      
      dvs:    the array of design variable values
      numDVs: the number of design variables
      '''
      # Get number of the design variables
      num_dvs = dvs.shape[0]
      self.ptr.setDesignVars(<TacsScalar*>dvs.data, num_dvs)
      return

   def getDesignVarRange(self, 
                         np.ndarray[TacsScalar,ndim=1,mode='c']lb,
                         np.ndarray[TacsScalar,ndim=1,mode='c']ub):
      '''
      Retrieve the design variable range.
      
      This call is collective on all TACS processes. The ranges 
      provided by indivdual objects may not be consistent (if 
      someone provided incorrect data they could be.) Make a 
      best guess; take the minimum upper bound and the maximum 
      lower bound.
      
      lowerBound: the lower bound on the design variables (output)
      upperBound: the upper bound on the design variables (output)
      numDVs:     the number of design variables
      '''
      
      # Get the number of design variables
      num_dvs = lb.shape[0]
      self.ptr.getDesignVarRange(<TacsScalar*>lb.data,
                                 <TacsScalar*>ub.data, 
                                 num_dvs)
      return
   
   def setNumThreads(self, int t):
      '''
      Set the number of threads to use in computation
      '''
      self.ptr.setNumThreads(t)
      return
      
   def createVec(self):
      '''
      Create a distributed vector.
      
      Vector classes initialized by one TACS object, cannot be 
      used by a second, unless they share are exactly the
      parallel layout.
      '''
      return _init_Vec(self.ptr.createVec())
  
   def createMat(self):
      '''
      Create a distributed matrix
      '''
      return _init_Mat(self.ptr.createMat())
   
   def createFEMat(self, OrderingType order_type=TACS_AMD_ORDER):
      '''
      Create a parallel matrix specially suited for finite-element
      analysis.
      
      On the first call, this computes a reordering with the scheme
      provided. On subsequent calls, the reordering scheme is reused s
      that all FEMats, created from the same TACSAssembler object have
      the same non-zero structure.  This makes adding matrices
      together easier (which is required for eigenvalue computations.)

      The first step is to determine the coupling nodes. (For a serial
      case there are no coupling nodes, so this is very simple!)
      Then, the nodes that are not coupled to other processes are
      determined. The coupling and non-coupling nodes are ordered
      separately.  The coupling nodes must be ordered at the end of
      the block, while the local nodes must be ordered first. This
      type of constraint is not usually imposed in matrix ordering
      routines, so here we use a kludge.  First, order all the nodes
      and determine the ordering of the coupling variables within the
      full set.  Next, order the local nodes. Tis hopefully reduces
      the fill-ins required, although there is no firm proof to back
      that up.
      
      The results from the reordering are placed in a set of
      objects. The matrix reordering is stored in feMatBIndices and
      feMatCIndices while two mapping objects are created that map the
      variables from the global vector to reordered matrix.

      Mathematically this reordering can be written as follows,

      A1 = (P A P^{T})
        
      where P^{T} is a permutation of the columns (variables), while P
      is a permutation of the rows (equations).
      '''
      return _init_Mat(self.ptr.createFEMat(order_type))
     
   def zeroVariables(self):
      '''
      Zero the entries of the local variables
      '''
      self.ptr.zeroVariables()
      return
      
   def zeroDotVariables(self):
      '''
      Zero the values of the time-derivatives of the state variables
      '''
      self.ptr.zeroDotVariables()
      return
    
   def zeroDDotVariables(self):
      '''
      Zero the values of the 2nd time-derivatives of the state
      variables
      '''
      self.ptr.zeroDDotVariables()
      return
   
   def setVariables(self, Vec stateVars):
      '''
      Set the values of the state variables
      '''
      self.ptr.setVariables(stateVars.ptr)
      return
      
   def getVariables(self, Vec stateVars):
      '''
      Get the values of the state variables
      '''
      self.ptr.getVariables(stateVars.ptr)
      return
    
   def setDotVariables(self, Vec stateVars):
      '''
      Set the values of the time-derivative of the state variables
      '''
      self.ptr.setDotVariables(stateVars.ptr)
      return
    
   def setDDotVariables(self, Vec stateVars):
      '''
      Set the values of the 2nd time derivative of the state 
      variables
      '''
      self.ptr.setDDotVariables(stateVars.ptr)
      return
    
   def assembleRes(self, Vec residual):
      '''
      Assemble the residual associated with the input load case.  
      
      This residual includes the contributions from element tractions
      set in the TACSSurfaceTraction class and any point loads. Note
      that the vector entries are zeroed first, and that the Dirichlet
      boundary conditions are applied after the assembly of the
      residual is complete.
      
      rhs:      the residual output
      '''
      self.ptr.assembleRes(residual.ptr)
      return
      
   def assembleJacobian(self, Vec residual, Mat A,
                        double alpha, double beta, double gamma,
                        MatrixOrientation matOr=NORMAL):
      '''
      Assemble the Jacobian matrix
      
      This function assembles the global Jacobian matrix and
      residual. This Jacobian includes the contributions from all
      elements. The Dirichlet boundary conditions are applied to the
      matrix by zeroing the rows of the matrix associated with a
      boundary condition, and setting the diagonal to unity. The
      matrix assembly also performs any communication required so that
      the matrix can be used immediately after assembly.

      residual:  the residual of the governing equations
      A:         the Jacobian matrix
      alpha:     coefficient on the variables
      beta:      coefficient on the time-derivative terms
      gamma:     coefficient on the second time derivative term
      matOr:     the matrix orientation NORMAL or TRANSPOSE
      '''
      self.ptr.assembleJacobian(residual.ptr, 
                                A.ptr, alpha, beta,
                                gamma, matOr)
      return

   def assembleMatType(self, ElementMatrixType matType,
                       Mat A, MatrixOrientation matOr):
      
      '''
      Assemble the Jacobian matrix
      
      This function assembles the global Jacobian matrix and
      residual. This Jacobian includes the contributions from all
      elements. The Dirichlet boundary conditions are applied to the
      matrix by zeroing the rows of the matrix associated with a
      boundary condition, and setting the diagonal to unity. The
      matrix assembly also performs any communication required so that
      the matrix can be used immediately after assembly.
        
      residual:  the residual of the governing equations
      A:         the Jacobian matrix
      alpha:     coefficient on the variables
      beta:      coefficient on the time-derivative terms
      gamma:     coefficient on the second time derivative 
      term 
      matOr:     the matrix orientation NORMAL or TRANSPOSE
      '''
      self.ptr.assembleMatType(matType, A.ptr, matOr)
      return
      
   def evalFunction(self, Function func):
      '''
      Evaluate a TACS function
      '''
      cdef int num_funcs = 1
      cdef TacsScalar val
      self.ptr.evalFunctions(&func.ptr, num_funcs, &val)
      return val

   def evalDVSens(self, Function func, 
                  np.ndarray[TacsScalar, ndim=1, mode='c'] fdvSens):
      '''
      Evaluate the derivative of a function w.r.t. the design
      variables.
      '''
      cdef int num_funcs = 1
      cdef int num_dvs = fdvSens.shape[0]
      self.ptr.evalDVSens(&func.ptr, num_funcs,
                          <TacsScalar*>fdvSens.data, num_dvs)

      return

   def evalSVSens(self, Function func, Vec vec):

      '''
      Evaluate the derivative of the function w.r.t. the state
      variables.
      
      function: the function pointer
      vec:      the derivative of the function w.r.t. the state variables 
      '''
      self.ptr.evalSVSens(func.ptr, vec.ptr)
      return

   def evalAdjointResProduct(self, Vec adjoint,
                             np.ndarray[TacsScalar,ndim=1,mode='c'] fdvSens):
      '''
      This function is collective on all TACSAssembler processes. This
      computes the product of the derivative of the residual
      w.r.t. the design variables with several adjoint vectors
      simultaneously. This saves computational time as the derivative
      of the element residuals can be reused for each adjoint
      vector. This function performs the same task as
      evalAdjointResProduct, but uses more memory than calling it for
      each adjoint vector.

      adjoint:     the array of adjoint vectors
      dvSens:      the product of the derivative of the residuals
                    and the adjoint
      num_dvs:      the number of design variables
      '''
      cdef int num_adjoints = 1
      cdef int num_dvs = fdvSens.shape[0]
      self.ptr.evalAdjointResProducts(&adjoint.ptr, num_adjoints,
                                      <TacsScalar*>fdvSens.data,
                                      num_dvs)
      return
        
   def testElement(self, int elemNum, int print_level):
      '''
      Test the implementation of the given element number.

      This tests the stiffness matrix and various parts of the
      design-sensitivities: the derivative of the determinant of the
      Jacobian, the derivative of the strain w.r.t. the nodal
      coordinates, and the state variables and the derivative of the
      residual w.r.t. the design variables and nodal coordinates.

      elemNum:     the element number to test
      print_level: the print level to use 
      '''
      self.ptr.testElement(elemNum, print_level)
      return

   def testConstitutive(self, int elemNum, print_level):
      '''
      Test the implementation of the given element constitutive
      class. 

      This function tests the failure computation and the mass
      sensitivities for the given element.

      elemNum:     the element to retrieve the constitutive object
         from 
      print_level: the print level to use for the test
      '''
      self.ptr.testConstitutive(elemNum, print_level)
      return

   def testFunction(self, Function func, int num_dvs,
                    double dh):
      '''
      Test the implementation of the function. 

      This tests the state variable sensitivities and the design
      variable sensitivities of the function of interest. These
      sensitivities are computed based on a random perturbation of
      the input values.  Note that a system of equations should be
      solved - or the variables should be set randomly before
      calling this function, otherwise this function may produce
      unrealistic function values. 

      Note that this function uses a central difference if the 
      real code is compiled, and a complex step approximation if 
      the complex version of the code is used.

      func:            the function to test
      num_dvs: the number of design variables to use
      dh:              the step size to use
      '''
      self.ptr.testFunction(func.ptr, num_dvs, dh)
      return

   def getMPIComm(self):
      '''
      Retrieve the MPI Communicator
      '''
      cdef MPI_Comm c_comm = self.ptr.getMPIComm()
      # comm = MPI.Comm()
      # comm.ob_mpi = c_comm
      return None
    
# Wrap the TACStoFH5 class
cdef class ToFH5:
   NODES = 1
   NODES = 1
   DISPLACEMENTS = 2
   STRAINS = 4
   STRESSES = 8
   EXTRAS = 16
   COORDINATES = 32

   cdef TACSToFH5 *ptr
   def __cinit__(self, Assembler tacs, ElementType elem_type,
                 int out_type):
      '''
      Create the TACSToFH5 file creation object
      
      input:
      tacs:        the instance of the TACSAssembler object
      elem_type:   the type of element to be used
      out_type:    the output type to write 
      '''
      self.ptr = new TACSToFH5(tacs.ptr, elem_type, out_type)
      self.ptr.incref()
      return
   
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return
    
   def setComponentName(self, int comp_num, char *group_name):
      '''
      Set the component name for the variable
      '''
      self.ptr.setComponentName(comp_num, &group_name[0])
      
   def writeToFile(self, char *filename):
      '''
      Write the data stored in the TACSAssembler object to filename
      '''
      self.ptr.writeToFile(&filename[0])

# Wrap the TACSCreator object
cdef class Creator:
   cdef TACSCreator *ptr
   def __cinit__(self, MPI.Comm comm, int vars_per_node):
      cdef MPI_Comm c_comm = comm.ob_mpi
      self.ptr = new TACSCreator(c_comm, vars_per_node)
      self.ptr.incref()
      return

   def __dealloc__(self):
      self.ptr.decref()
      return

   def setGlobalConnectivity(self, int num_nodes,
                             np.ndarray[int, ndim=1, mode='c'] elem_node_ptr, 
                             np.ndarray[int, ndim=1, mode='c'] elem_node_conn,
                             np.ndarray[int, ndim=1, mode='c'] elem_id_nums):
      '''Set the connectivity and the element id numbers'''
      cdef int num_elements = elem_node_ptr.shape[0]-1
      assert(num_elements == elem_id_nums.shape[0])
      self.ptr.setGlobalConnectivity(num_nodes, num_elements,
                                     <int*>elem_node_ptr.data,
                                     <int*>elem_node_conn.data,
                                     <int*>elem_id_nums.data)
      return

   def setBoundaryConditions(self, np.ndarray[int, ndim=1, mode='c'] nodes,
                             np.ndarray[int, ndim=1, mode='c'] bcvars,
                             np.ndarray[int, ndim=1, mode='c'] ptr):
      '''Set the boundary conditions'''
      cdef int num_bcs = nodes.shape[0]
      self.ptr.setBoundaryConditions(num_bcs, <int*>nodes.data,
                                     <int*>bcvars.data, <int*>ptr.data)
      return

   def setDependentNodes(self, np.ndarray[int, ndim=1, mode='c'] dep_ptr,
                         np.ndarray[int, ndim=1, mode='c'] dep_conn,
                         np.ndarray[double, ndim=1, mode='c'] dep_weights):


      return

   def setElements(self, elements):
      '''Set the elements'''

      # Allocate an array for the element pointers
      cdef TACSElement **elems
      elems = <TACSElement**>malloc(len(elements)*sizeof(TACSElement*))
      if elems is NULL:
         raise MemoryError()

      for i in xrange(len(elements)):
         elems[i] = (<Element>elements[i]).ptr

      self.ptr.setElements(elems, len(elements))

      # Free the allocated array
      free(elems)

      return
   
   def setNodes(self, np.ndarray[TacsScalar, ndim=1, mode='c'] Xpts):
      self.ptr.setNodes(<TacsScalar*>Xpts.data)
      return

   def setReorderingType(self, OrderingType order_type,
                         MatrixOrderingType mat_type):
      self.ptr.setReorderingType(order_type, mat_type)
      return

   def createTACS(self):
      return _init_Assembler(self.ptr.createTACS())
    
# Wrap the TACSMeshLoader class
cdef class MeshLoader:
   cdef TACSMeshLoader *ptr
   def __cinit__(self, MPI.Comm comm):
      '''
      This is an interface for reading the NASTRAN-style files i.e. BDF
      '''
      cdef MPI_Comm c_comm = comm.ob_mpi
      self.ptr = new TACSMeshLoader(c_comm)
      self.ptr.incref()
      
   def __dealloc__(self):
      self.ptr.decref()

   def scanBDFFile(self, char*filename):
      '''
      This scans a Nastran file - only scanning in information from the
      bulk data section
      
      The only entries scanned are the entries beginning with elem_types
      and any GRID/GRID* entries
      '''
      self.ptr.scanBDFFile(&filename[0])

   def getNumComponents(self):
      '''
      Return the number of components
      '''
      return self.ptr.getNumComponents()
   
   def getComponentDescript(self, int comp_num):
      '''
      Return the component description
      '''
      cdef bytes py_string
      py_string = self.ptr.getComponentDescript(comp_num)
      return py_string
   
   def getElementDescript(self, int comp_num):
      '''
      Retrieve the element description corresponding to the component number
      '''
      cdef bytes py_string
      py_string = self.ptr.getElementDescript(comp_num)
      return py_string
   
   def setElement(self, int comp_num, Element elem):
      '''
      Set the element associated with a given component number
      '''
      self.ptr.setElement(comp_num, elem.ptr)
      
   def getNumNodes(self):
      return self.ptr.getNumNodes()

   def getNumElement(self):
      return self.ptr.getNumElements()
    
   def createTACS(self, int varsPerNode, OrderingType order_type=NATURAL_ORDER,
                  MatrixOrderingType mat_type=DIRECT_SCHUR):
      '''
      Create a distribtued version of TACS
      '''
      return _init_Assembler(self.ptr.createTACS(varsPerNode,
                                                 order_type, mat_type))

   def getConnectivity(self):
      cdef int num_nodes
      cdef int num_elements
      cdef const int *elem_ptr
      cdef const int *elem_conn
      cdef const double *Xpts

      self.ptr.getConnectivity(&num_nodes, &num_elements,
                               &elem_ptr, &elem_conn, &Xpts)

      ptr = np.zeros(num_elements+1, dtype=np.int)
      for i in xrange(num_elements+1):
         ptr[i] = elem_ptr[i]

      conn = np.zeros(ptr[-1], dtype=np.int)
      for i in xrange(ptr[-1]):
         conn[i] = elem_conn[i]

      X = np.zeros(3*num_nodes)
      for i in xrange(3*num_nodes):
         X[i] = Xpts[i]

      return ptr, conn, X
