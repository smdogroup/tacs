# For the use of MPI
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy 
import numpy as np
cimport numpy as np

# Ensure that numpy is initialized
np.import_array()

# Import the definition required for const strings
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

cdef extern from "mpi-compat.h":
   pass

# This class wraps a C++ array with a numpy array for later useage
cdef class NpArrayWrap:
   cdef int nptype
   cdef int dim1, dim2
   cdef void *data_ptr

   cdef set_data1d(self, int nptype, int dim1, void *data_ptr):
      '''Set data in the array'''
      self.nptype = nptype
      self.dim1 = dim1
      self.dim2 = -1
      self.data_ptr = data_ptr
      return

   cdef set_data2d(self, int nptype, int dim1, int dim2, void *data_ptr):
      '''Set data in the array'''
      self.nptype = nptype
      self.dim1 = dim1
      self.dim2 = dim2
      self.data_ptr = data_ptr
      return

   cdef as_ndarray(self):
      '''Return a numpy version of the array'''
      # Set the shape of the array
      cdef int size = 1
      cdef np.npy_intp shape[2]
      cdef np.ndarray ndarray

      shape[0] = <np.npy_intp> self.dim1
      if (self.dim2 > 0):
         size = 2
         shape[1] = <np.npy_intp> self.dim2
      
      # Create the array itself
      ndarray = np.PyArray_SimpleNewFromData(size, shape,
                                             self.nptype, self.data_ptr)
      
      # Set the base class who owns the memory
      ndarray.base = <PyObject*>self
      Py_INCREF(self)
      
      return ndarray
  
# Python class for corresponding instance TACSObject
cdef class pyTACSObject:
    cdef TACSObject *this_ptr
    def __cinit__(self):
        '''
        TACSObject: The base class for all TACS objects to enable reference
        counting. In most cases this is sufficient to free any allocated 
        memory.
        '''
        self.this_ptr = new TACSObject()
        return
    def __dealloc__(self):
        del self.this_ptr

    
# Python class for corresponding instance TACSAssembler
cdef class pyTACSAssembler(pyTACSObject):
    cdef TACSAssembler *this_ptr
    def __cinit__(self, MPI.Comm comm, int numOwnedNodes,
                  int varsPerNode, int numElements, int numNodes,
                  int numDependentNodes, int nodeMaxCSRsize):
        '''
        Constructor for the TACSAssembler object

        input:
        comm:                the TACS communicator 
        numOwnedNodes:       the number of locally-owned nodes
        varsPerNode:         the number of degrees of freedom per node
        numElements:         the number of elements in the mesh
        numNodes:            the number of nodes in the mesh
        numDependentNodes:   the number of dependent nodes in the mesh
        nodeMaxCSRSize:      the maximum CSR size - expanded automatically
        '''
        # Convert the communicator
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.this_ptr = new TACSAssembler(c_comm, numOwnedNodes,
                                          varsPerNode, 
                                          numElements, numNodes,
                                          numDependentNodes,
                                          nodeMaxCSRsize)
        # Need to add incref
    def __dealloc__(self):
        del self.this_ptr
        # Need to decref


    def getNumNodes(self):
        '''
        Return the number of nodes in the TACSAssembler
        '''
        return self.this_ptr.getNumNodes()

    def getNumDependentNodes(self):
        '''
        Return the number of dependent nodes
        '''
        return self.this_ptr.getNumDependentNodes()
    
    def getNumElements(self):
        '''
        Return the number of elements
        '''
        return self.this_ptr.getNumElements()

    def getVarMap(self):
        '''
        Return VarMap object
        '''
        return self.this_ptr.getVarMap()
    
    def addNode(self, int localNodeNum, int tacsNodeNum):
        '''
        Add a node to TACS. 

        Each node is associated with a local and a global numbering scheme.
        '''
        self.this_ptr.addNode(localNodeNum, tacsNodeNum)

    
    def addNodes(self, np.ndarray[int, ndim=1, mode='c'] localNodeNums, 
                 np.ndarray[int, ndim=1, mode='c'] tacsNodeNums1 = NULL):
        '''
        Add a list of tacsNodeNums, no option to add all at once
        '''
        # Get size of localNodeNums
        num_nodes = localNodeNums.shape[0]
        self.this_ptr.addNodes(<int*>localNodeNums.data,
                               <int*>tacsNodeNums1.data,num_nodes)
        
    def setDependentNodes(self, 
                          np.ndarray[int, ndim=1, mode='c']depNodeIndex, 
                          np.ndarray[int, ndim=1, mode='c']depNodeTolLocal,
                          np.ndarray[TacsScalar, ndim=1, mode='c']
                          depNodeWeights):
        '''
        Set the dependent node data structure
        '''
        self.this_ptr.setDependentNodes(<int**>depNodeIndex, 
                                        <int**>depNodeToLocal,
                                        <TacsScalar**>depNodeWeights)
        
   
    def addBC(self, int nodeNum, 
              np.ndarray[int, ndim=1,mode='c']bcNums, 
              np.ndarray[TacsScalar, ndim=1, mode='c']bcVals= NULL):
        '''
        Set the boundary conditions into the object that will be 
        associated with the vector/matrices and with the 
        TACSAssembler object
        '''
        # Get number of bcs
        nbcs = bcNums.shape[0]
        if bcVals is NULL:
            self.this_ptr.addBC(nodeNum, <int*> bcNums.data, nbcs)
        else:
            self.this_ptr.addBC(nodeNum, <int*>bcNums.data, bcVals,
                                nbcs)
           
    def addElement(self, pyTACSElement element, 
                   np.ndarray[int, ndim=1, mode='c']localNodeNums):
        '''
        Add an element to TACS with given node numbers
        '''     
        # Get the number of nodes in element
        num_nodes = localNodeNums.shape[0]
        return self.this_ptr.addElement(element.this_ptr, 
                                        <int*>localNodeNums.data)
    
    def setElementConnectivity(self, ):
        '''
        Set the element connectivity
        '''
        self.this_ptr.setElementConnectivity()
        
    def setNodes(self, pyBVec X):
        '''
        Set the nodes from the node map object
        '''    
        self.this_ptr.setNodes(X.this_ptr)
    
    def getNodes(self, pyBVec X):
        '''
        Get the nodes from the node map object
        '''
        self.this_ptr.getNodes(X.this_ptr)
        
    def computeCouplingNodes():
        '''
        Compute the local node numbers that correspond to the 
        coupling nodes connected to elements on other processes.

        Sort the global node numbers. Match the intervals and send 
        them off to the owning process. On the owner, scan through 
        the arrays until all the local coupling nodes are found. 
        '''
        return
    
    def computeCouplingElements():
        '''
        Compute the elements that couple with other processors.
        
        Compute the coupling nodes and the node to element pointer
        CSR data structure. From these, collect all elements that 
        "own" a node that is referred to from another process.
        '''
        return

    def computeReordering(self):
        '''
        Compute a reordering of the nodes.
        
        1. Determine the nodes that are locally owned and those that
        are required from other processors.
        2. Distribute the unknown nodes to the owning processes. 
        These are the recv_nodes -> sort them to sorted_recv_nodes.
        3. Find the indices of the recieved nodes.
        4. Order the locally owned nodes based on the input ordering
        5. Based on the recvied indices, set the outgoing node 
        numbers back into the recieving array.
        6. Set the new values of the nodes on the requesting 
        processes
    
        The nodes are transfered as follows
        extern_nodes -> sorted_extern_nodes -> recv_nodes 
        ->sorted_recv_nodes-> sorted_recv_index
    
        Node ordering
    
        recv_nodes[i] = new_nodes[sorted_recv_index[i]]
        recv_nodes -> new_sorted_extern_nodes
        find extern_nodes[i] == sorted_extern_nodes[j]
        assign new_var = new_sorted_extern_nodes[j]
        '''
    
        self.this_ptr
    
    def computeMatReordering(self, ):
        '''
        Compute the reordering for the given matrix.
        
        This uses either Reverse Cuthill-McKee (RCM_ORDER),
        Approximate Minimum Degree (AMD) or Nested Disection (ND) to
        compute a reordering of the variables.
    
        The input to the function is a CSR-type data structure of 
        the matrix. Note that for ND (with the external package
        Metis), requires that the diagonal be eliminated from the 
        CSR data structure. (This modified data structure can be
        computed with the no_diagonal flag when calling the CSR
        creation routines.) 
    
        The matrix ordering routines compute a reordering such that:
        
        P * A * P^{T} has fewer non-zeros.
        
        The function returns an array new_vars such that:
        
        new_vars = P^{T} old_vars
        perm = P
        '''
        return
    
    
    def computeNodeToElementCSR(self):
        '''
        The following function creates a data structure that links
        nodes to elements - this reverses the existing data 
        structure that links elements to nodes but keeps the 
        original intact. 
    
        The algorithm proceeds as follows:
        1. The size of the arrays are determined by finding how many
        nodes point to each element
        2. The index into the nodeElem array is determined by adding
        up the contributions from all previous entries.
        3. The oirginal data structure is again traversed and this
        time an element number is associated with each element.
        '''
        return
    def computeLocalNodeToNodeCSR(self,):
        '''
        Set up a CSR data structure pointing from local nodes to 
        other local nodes.
    
        input:
        nodiag = Remove the diagonal matrix entry 
        
        output:
        rowp = the row pointer corresponding to CSR data structure
        cols = the column indices for each row of the CSR data 
        structure
    
        This function works by first estimating the number of 
        entries in each row of the matrix. This information is 
        stored temporarily in the array rowp. After the 
        contributions from the elements and sparse constraints are
        added, the preceeding elements in rowp are added such that
        rowp becomes the row pointer for the matrix. Note that this 
        is an upper bound because the elements and constraints may 
        introduce repeated variables. Next, cols is allocated corre
        sponding to the column index for each entry. This iterates
        back over all elements and constraints. At this stage, rowp
        is treated as an array of indices, that index into the i-th
        row of cols[:] where the next index should be inserted. As a
        result, rowp must be adjusted after this operation is
        completed.  The last step is to sort and uniquify 
        each row of the matrix.
        '''
        return    
    
    def finalize(self):
        '''
        The function finalize performs a number of synchronization 
        tasks that prepare the finite-element model for use.
    
        tacsNodeNums[i] is the global node number for the local 
        node number i
    
        Two objects are required:
        1. VarMap is constructed with the block sizes of each 
        node owned by this process
        
        2. VecDistribute is constructed so that it takes an array 
        and distributes its values to a vector or takes the vector
        values and collects them into an array This requires a 
        sorted array of global node numbers.
        '''
        self.this_ptr.finalize()

    def getNumDesignVars(self):
        '''
        Return the number of design variables
        '''
        return self.this_ptr.getNumDesignVars()
        
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
        self.this_ptr.getDesignVars(<TacsScalar*>dvs.data, num_dvs)

    
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
        self.this_ptr.setDesignVars(<TacsScalar*>dvs.data, num_dvs)

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
        self.this_ptr.getDesignVarRange(<TacsScalar*>lb.data,
                                        <TacsScalar*>ub.data, 
                                        num_dvs)
    
    def setNumThreads(self, int t):
        '''
        Set the number of threads to use in computation
        '''
        self.this_ptr.setNumThreads(t)
        
    def createVec(self):
        '''
        Create a distributed vector.
        
        Vector classes initialized by one TACS object, cannot be 
        used by a second, unless they share are exactly the
        parallel layout.
        '''
        return self.this_ptr.createVec()
    
    def createMat(self):
        '''
        Create a distributed matrix
        '''
        return self.this_ptr.createMat()

    def createFEMat(self):
        '''
        Create a parallel matrix specially suited for finite-element
        analysis.
        
        On the first call, this computes a reordering with the 
        scheme provided. On subsequent calls, the reordering scheme
        is reused s that all FEMats, created from the same 
        TACSAssembler object have the same non-zero structure. 
        This makes adding matrices together easier (which is 
        required for eigenvalue computations.)

        The first step is to determine the coupling nodes. (For a 
        serial case there are no coupling nodes, so this is very 
        simple!)  Then, the nodes that are not coupled to other 
        processes are determined. The coupling and non-coupling 
        nodes are ordered separately.  The coupling nodes must be 
        ordered at the end of the block, while the local nodes must
        be ordered first. This type of constraint is not usually 
        imposed in matrix ordering routines, so here we use a 
        kludge.  First, order all the nodes and determine the 
        ordering of the coupling variables within the full set. 
        Next, order the local nodes. Tis hopefully reduces the 
        fill-ins required, although there is no firm proof to back 
        that up.

        The results from the reordering are placed in a set of 
        objects. The matrix reordering is stored in feMatBIndices 
        and feMatCIndices while two mapping objects are created 
        that map the variables from the global vector to 
        reordered matrix.

        Mathematically this reordering can be written as follows,

        A' = (P A P^{T})
        
        where P^{T} is a permutation of the columns (variables), 
        while 
        P is a permutation of the rows (equations).
        '''
        return self.this_ptr.createFEMat()
     
    def zeroVariables(self):
        '''
        Zero the entries of the local variables
        '''
        self.this_ptr.zeroVariables()
    
    def zeroDotVariables(self):
        '''
        Zero the values of the time-derivatives of the state 
        variables
        '''
        self.this_ptr.zeroDotVariables()
    
    def zeroDDotVariables(self):
        '''
        Zero the values of the 2nd time-derivatives of the state 
        variables 
        '''
        self.this_ptr.zeroDDotVariables()
    
    def setVariables(self, pyBVec stateVars):
        '''
        Set the values of the state variables
        '''
        self.this_ptr.setVariables(stateVars.this_ptr)

    def getVariables(self, pyBVec stateVars):
        '''
        Get the values of the state variables
        '''
        self.this_ptr.getVariables(stateVars.this_ptr)
    
    def setDotVariables(self, pyBVec stateVars):
        '''
        Set the values of the time-derivative of the state variables
        '''
        self.this_ptr.setDotVariables(stateVars.this_ptr)
    
    def setDDotVariables(self, pyBVec stateVars):
        '''
        Set the values of the 2nd time derivative of the state 
        variables
        '''
        self.this_ptr.setDDotVariables(stateVars.this_ptr)

    def getTacsNodeNums(self, 
                        np.ndarray[int,ndim=1,mode='c']localNodes):
        
        '''
        Given an array of local node numbers, return the associated
        TACS numbers in the same array - over-writing the values. 
        
        localNodes: the local node numbers to retrieve
        '''
        # Get the number of local nodes 
        num_nodes = localNodes.shape[0]
        self.this_ptr.getTacsNodeNums(<int*>localNodes.data,
                                      num_nodes)
    
    def assembleRes(self, pyBVec residual):
        '''
        Assemble the residual associated with the input load case.  
  
        This residual includes the contributions from element
        tractions set in the TACSSurfaceTraction class and any point
        loads. Note that the vector entries are zeroed first, and 
        that the Dirichlet boundary conditions are applied after the
        assembly of the residual is complete.

        rhs:      the residual output
        '''
        self.this_ptr.assembleRes(residual.this_ptr)

    def assembleResNoBCs(self, pyBVec residual):
        '''
        Assemble the residuals of the finite-element problem but do
        not apply the boundary conditions. This is useful for
        determing the reaction forces.
    
        residual:      the residual (output)
        '''
        self.this_ptr.assembleResNoBCs(residual.this_ptr)

    def assembleJacobian(self, pyBVec residual, pyTACSMat A,
                         double alpha, double beta, double gamma,
                         MatrixOrientation matOr):
        '''
        Assemble the Jacobian matrix

        This function assembles the global Jacobian matrix and
        residual. This Jacobian includes the contributions from all
        elements. The Dirichlet boundary conditions are applied to 
        the matrix by zeroing the rows of the matrix associated 
        with a boundary condition, and setting the diagonal to
        unity. The matrix assembly also performs any communication
        required so that the matrix can be used immediately after
        assembly. 

        residual:  the residual of the governing equations
        A:         the Jacobian matrix
        alpha:     coefficient on the variables
        beta:      coefficient on the time-derivative terms
        gamma:     coefficient on the second time derivative term
        matOr:     the matrix orientation NORMAL or TRANSPOSE
        '''
        self.this_ptr.assembleJacobian(residual.this_ptr, 
                                       A.this_ptr, alpha, beta,
                                       gamma, matOr)

    def assembleMatType(self, ElementMatrixTypes matType,
                        pyTACSMat A,
                        MatrixOrientation matOr):
        
        '''
        Assemble the Jacobian matrix
        
        This function assembles the global Jacobian matrix and
        residual. This Jacobian includes the contributions 
        from all elements. The Dirichlet boundary conditions are
        applied to the matrix by zeroing the rows of the matrix
        associated with a boundary condition, and setting the
        diagonal to unity. The matrix assembly also performs any
        communication required so that the matrix can be 
        used immediately after assembly.
        
        residual:  the residual of the governing equations
        A:         the Jacobian matrix
        alpha:     coefficient on the variables
        beta:      coefficient on the time-derivative terms
        gamma:     coefficient on the second time derivative 
                       term 
        matOr:     the matrix orientation NORMAL or TRANSPOSE
        '''
        self.this_ptr.assembleMatType(matType, A.this_ptr, 
                                      matOr)
   
    def evalFunctions(self, pyTACSFunction functions, int numFuncs,
                      np.ndarray[TacsScalar, ndim=1,mode='c']fVals):
        '''
        Evaluate a list of TACS functions

        First, check if the functions are initialized. Obtain the
        number of iterations over the function domain required to
        evaluate the functions.

        This function will print an error and return 0 if the
        underlying TACSAssembler object does not correspond to the
        TACSAssembler object .
    
        functions:  array of functions to evaluate
        numFuncs:   the number of functions to evaluate

        output:
        funcVals: the values of the functions 
        '''
        self.this_ptr.evalFunctions(functions.this_ptr, numFuncs,
                                    <TacsScalar*>fVals.data)

    def evalDVSens(self, pyTACSFunction funcs, int numFuncs, 
                   np.ndarray[TacsScalar, ndim=1, mode='c']fdvSens,
                   int num_dvs):
        '''
        Evaluate the derivative of a list of functions w.r.t. the
        design variables.

        Note that a function should be evaluated - using 
        evalFunction - before its derivatives can be evaluated.

        The design variable sensitivities are divided into two
        distinct sets: material-dependent design variables and shape
        design variables. The shape design variables are handled
        through the TACSNodeMap class. The material-dependent design
        variables are handled through the element classes 
        themselves.

        In this code, the derivative of the function w.r.t. the
        shape-dependent design variables is handled first. The
        derivative of the function w.r.t each nodal location is
        determined. The TACSNodeMap object (if not NULL) is then 
        used to determine the derivative of the nodal locations 
        w.r.t. the design variables themselves. 
  
        The material-dependent design variables are handled on an
        element-by-element and traction-by-traction dependent basis.

        Note that this function distributes the result to the
        processors through a collective communication call. No 
        further parallel communication is required.

        input:
        funcs:     the TACSFunction function objects
        numFuncs:  the number of functions - size of funcs array
        fdvSens:   the sensitivity - size numFuncs*numDVs
        numDVs:    the number of design variables
        '''
        self.this_ptr.evalDVSens(funcs.this_ptr, numFuncs,
                                 <TacsScalar*> fdvSens.data, numDVs)

    def evalXptSens(self, pyTACSFunction funcs, numFuncs, 
                    np.ndarray[TacsScalar,ndim=1,mode='c']fXptSens):
        '''
        Evaluate the derivative of the function w.r.t. the owned
        nodes. 
        
        This code evaluates the sensitivity of the function w.r.t. 
        the owned nodes for all elements in the function domain. 
        
        Note that a function should be evaluated - using 
        evalFunction - before its derivatives can be evaluated.
        
        This function should be preferred to the use of evalDVSens
        without a list of functions since it is more efficient!

        input:
        funcs:     the TACSFunction function objects
        numFuncs:  the number of functions - size of funcs array
        fXptSens:  the sensitivity - size numFuncs*numNodes*3
        '''
        self.this_ptr.evalXptSens(funcs.this_ptr, numFuncs, 
                                  <TacsScalar*>fXptSens.data)

    def evalSVSens(self, pyTACSFunction function, pyBVec vec):
        '''
        Evaluate the derivative of the function w.r.t. the state
        variables. 

        This code evaluates the sensitivity of the function w.r.t. 
        the state variables for all elements in the function
        domain. This code is usually much faster than the code for
        computing the derivative of the function w.r.t. the design
        variables.  

        Note that the sensitivity vector 'vec' is assembled, and
        appropriate boundary conditions are imposed before the
        function is returned. 

        function: the function pointer
        vec:      the derivative of the function w.r.t. the state
                  variables 
        '''
        self.this_ptr.evalSVSens(function.this_ptr, vec.this_ptr)

    def evalAdjointResProducts(self, pyBVec adjoint,int numAdjoint, 
                               np.ndarray[TacsScalar,ndim=1,mode='c']dvSens,
                               int num_dvs):
        '''
        This function is collective on all TACSAssembler
        processes. This computes the product of the derivative of 
        the residual w.r.t. the design variables with several 
        adjoint vectors simultaneously. This saves computational 
        time as the derivative of the element residuals can be 
        reused for each adjoint vector. This function performs the 
        same task as evalAdjointResProduct, but uses more memory 
        than calling it for each adjoint vector.

        adjoint:     the array of adjoint vectors
        numAdjoints: the number of adjoint vectors
        dvSens:      the product of the derivative of the residuals
                     and the adjoint
        num_dvs:      the number of design variables
        '''
        self.this_ptr.evalAdjointResProducts(adjoint.this_ptr, 
                                             numAdjoints, 
                                             <TacsScalar*>dvSens.data,
        int num_dvs)

    # def evalAdjointResProductsExperimental(self, pyBVec adjoint, 
    #                                        int numAdjoints, 
    #                                        np.ndarray[TacsScalar,ndim=1,mode='c']dvSens,
    #                                        int num_dvs):
    #     '''
    #     Evaluate the product of several ajdoint vectors with the
    #     derivative of the residual w.r.t. the design variables.

    #     This function is collective on all TACSAssembler
    #     processes. This computes the product of the derivative of 
    #     the residual w.r.t. the design variables with several 
    #     adjoint vectors simultaneously. This saves computational 
    #     time as the derivative of the element residuals 
    #     can be reused for each adjoint vector. This function 
    #     performs the same task as evalAdjointResProduct, but uses 
    #     more memory than calling it for each adjoint vector.
        
    #     This implementation is experimental because it is not yet
    #     fully implemented for all element types, does not handle
    #     geometric variables and does not yet cover design variables
    #     in tractions. Nevertheless, the short-term plan is to finish
    #     this implementation b/c it has the potential to be faster.

    #     adjoint:     the array of adjoint vectors
    #     numAdjoints: the number of adjoint vectors
    #     dvSens:      the product of the derivative of the residuals
    #                  and the adjoint   
    #     numDVs:      the number of design variables
    #     '''
    #     self.this_ptr.evalAdjointResProductsExperiemental(adjoint.this_ptr,
    #     numAdjoint, <TacScalar*>dvSens.data, num_dvs)

    def evalMatDVSensInnerProduct(self, TacsScalar alpha, 
                                  ElementMatrixTypes matType,
                                  pyBVec psi, pyBVec phi,
                                  np.ndarray[TacsScalar, ndim=1,
                                  mode='c']dvSens):
        '''
        Evaluate the derivative of an inner product of two vectors
        with a matrix of a given type. This code does not explicitly
        evaluate the element matrices.  Instead, the inner product
        contribution from each element matrix is added to the final
        result. This implementation saves considerable computational
        time and memory. 

        input:
        alpha:     the scaling parameter applied to the derivative
        matType:   the matrix type
        psi:       the left-multiplying vector
        phi:       the right-multiplying vector
        numDVs:    the length of the design variable array
        
        output:
        dvSens:    the derivative of the inner product 
        '''
        # Obtain number of design variables
        num_dvs = dvSens.shape[0]
        self.this_ptr.evalMatDVSensInnerProduct(alpha, matType,
                                                psi.this_ptr,
                                                phi.this_ptr,
                                                <TacsScalar*>dvSens.data, num_dvs)

    def evalMatSVSensInnerProduct(self, TacsScalar alpha, 
                                  ElementMatrixTypes matType,
                                  pyBVec psi, pyBVec phi,
                                  pyBVec res):
        '''
        Evaluate the derivative of the inner product of two vectors
        with a matrix with respect to the state variables. This is
        only defined for  nonlinear matrices, like the geometric
        stiffness matrix.  Instead of computing the derivative of 
        the matrix for each vector component and then computing the
        inner product, this code computes the derivative of the 
        inner product directly, saving computational time and memory
        
        input:
        
        alpha:     the scaling parameter applied to the derivative
        matType:   the matrix type
        psi:       the left-multiplying vector
        phi:       the right-multiplying vector
        numDVs:    the length of the design variable array
        
        output:
        res:       the derivative of the inner product w.r.t. the
                   state vars 
        '''
        self.this_ptr.evalMatSVSensInnerProduct(alpha, matType, 
                                                psi.this_ptr, 
                                                phi.this_ptr,
                                                res.this_ptr)
        
    def getElements(self):
        '''
        Return the TACSElement object
        '''
        return self.this_ptr.getElements()
        
    def getElement(self, int elemNum, 
                   np.ndarray[TacsScalar,ndim=1,mode='c']elemVars,
                   np.ndarray[TacsScalar,ndim=1,mode='c']elemXpts):
        '''
        Given the load case number and the element, return the 
        element object and the values of the variables and nodal
        locations associated with that element.

        This is useful for post-processing data without having to
        write a full function class.

        elemNum:     the element number
        elemVars:    the variables associated with elemNum
        elemXpts:    the nodal locations associated with elemNum
        
        returns:     the element pointer associated with elemNum
        '''
        return self.this_ptr.getElement(elemNum, 
                                        <TacsScalar*>elemVars.data,
                                        <TacsScalar*>elemXpts.data)

    def testElement(self, int elemNum, int print_level):
        '''
        Test the implementation of the given element number.

        This tests the stiffness matrix and various parts of the
        design-sensitivities: the derivative of the determinant of 
        the Jacobian, the derivative of the strain w.r.t. the nodal
        coordinates, and the state variables and the derivative of 
        the residual w.r.t. the design variables and nodal
        coordinates. 
 
        elemNum:     the element number to test
        print_level: the print level to use 
        '''
        self.this_ptr.testElement(elemNum, print_level)

    def testConstitutive(self, int elemNum, print_level):
        '''
        Test the implementation of the given element's constitutive
        class. 
  
        This function tests the failure computation and the mass
        sensitivities for the given element.
    
        elemNum:     the element to retrieve the constitutive object
                     from 
        print_level: the print level to use for the test
        '''
        self.this_ptr.testConstitutive(elemNum, print_level)
        
    def testFunction(self, pyTACSFunction func, int num_dvs,
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
        self.this_ptr.testFunction(func.this_ptr, num_dvs, dh)

    def getNumComponents(self):
        '''
        Determine the number of components defined by elements in 
        the TACSAssembler object.
        
        This call is collective - the number of components is 
        obtained by a global reduction.
        '''
        return this_ptr.getNumComponents()

    def getOutputNodeRange(self, ElementType elem_type,
                           np.ndarray[int, ndim=1,
                           mode='c']node_range):
        '''
        Return the output nodal ranges. These may be used to 
        determine what range of node numbers need to be determined 
        by this process.  
        '''
        self.this_ptr.getOutputNodeRange(elem_type, 
                                         <int**>node_range.data)

    def getOutputConnectivity(self, ElementType elem_type):
        '''
        Given the element type, determine the connectivity of the
        global data structure. Record the component number for each
        element within the data structure

        input:
        elem_type: the type of element eg SHELL, EULER_BEAM, SOLID
        etc. 
        (see all the element types in Element.h)
        
        output:
        component_nums: an array of size nelems of the component
                        numbers 
        csr:            the csr element->nodes information for this
                        class
        csr_range:      the range of csr data on this processor
        node_range:     the range of nodal values on this processor

        '''
        
        return
    def getOutputData(self, ElementType elem_type, int out_type, 
                      np.ndarray[double, ndim=1, mode='c']Data,
                      int nvals):
        '''
        Go through each element and get the output data for that
        element. 
        The data is stored point-wise with each variable stored
        contiguously for each new point within the connectivity
        list. This stores the data at a point in memory indexed by
        data[node*nvals]. However, fewer than 'nvals' entries may be
        written in this call. The remaining data may be design
        variable entries that are computed below.
    
        elem_type: the element type to match
        out_type:  the output type 
        data:      the data array - nvals x the number of elements
        nvals:     the number of values to skip at each point
        '''
        self.this_ptr.getOutputData(elem_type, out_type,
                                    <double*>Data.data, nvals)

    def getMPIComm(self):
        '''
        Retrieve the MPI Communicator
        '''
        return self.this_ptr.getMPIComm()

# Wrap the abstract vector base class
cdef class pyTACSVec(pyTACSObject):
    cdef CyTACSVec *this_ptr
    def __cinit__(self):
        '''
        The abstract vector class. 

        All TACS KSM objects use this interface for matrix-vector products
        and preconditioning operations. In practice, within TACS itself,
        most of the matricies require the BVec base class. However, when
        coupling to other disciplines it is frequently convenient to define
        a hybrid vector consisting of multiple components. In these cases,
        it is handy to have a base vector class. 
        '''
        self.this_ptr = new CyTACSVec()
        self.this_ptr.setSelfPointer(<void*>self)

    def __dealloc__(self):
        del this_ptr

cdef class pyVarMap(pyTACSObject):
    cdef VarMap *this_ptr
    def __cinit__(self, MPI.Comm comm, int N, int bsize):
        '''
        Variable map for the parallel distribution of a vector
        '''
        # Convert the communicator
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.this_ptr = new VarMap(c_comm, N, bsize)

    def __dealloc(self):
        del this_ptr

cdef class pyBCMap(pyTACSObject):
    cdef BCMap *this_ptr
    def __cinit__(self, int num_bcs):
        '''
        Define boundary conditions that are applied after all 
        the matrix/vector values have been set.
        '''
        self.this_ptr = new BCMap(num_bcs)
    def __dealloc__(self):
        del this_ptr

cdef class pyBVec(pyTACSVec):
    cdef BVec *this_ptr
    def __cinit__(self, pyVarMap rmap, pyBCMap bcs):
        '''
        The basic parallel block-based vector class, distributed
        based on the variable map provided. The boundary conditions
        zero the rows.
        '''   
        self.this_ptr = new BVec(rmap.this_ptr, bcs.this_ptr)

    def __dealloc__(self):
        del this_ptr

    # Basic vector operations
    # def getSize(self, int *size):
    #     '''
    #     Get the local size of the vector on this processor
    #     '''
    #     return
    
    def norm(self):
        '''
        Compute the norm of the vector
        '''
        return self.this_ptr.norm()

    def scale(self, TacsScalar alpha):
        '''
        Scale the vector by alpha
        '''
        self.this_ptr.scale(alpha)
    
    def dot(self, pyTACSVec tvec):
        '''
        Compute the dot product of two vectors
        '''
        return self.this_ptr.dot(tvec.this_ptr)
        
    def mdot(self, np.ndarray[pyTACSVec,ndim=2, mode='c'] tvec,
             np.ndarray[TacsScalar, ndim=1,mode='c']ans):
        '''
        Compute multiple dot products. This is more efficient for parallel
        computations since there are fewer gather operations for the same
        number of dot products.
        '''
        # Obtain the number of rows in tvecs
        nvec = tvec.shape[1]
        self.this_ptr.mdot(tvec.this_ptr, <TacsScalar*>ans.data, nvecs)

    def axpy(self, TacsScalar alpha, pyTACSVec tvec):
        '''
        Compute y = alpha*x + y
        '''
        self.this_ptr.axpy(alpha, tvec.this_ptr)

    def axpby(self, TacsScalar alpha, TacsScalar beta,
              pyTACSVec tvec):
        '''
        Compute x <- alpha * vec + beta * x
        '''
        self.this_ptr.axpby(alpha, beta, tvec.this_ptr)

    def copyValues(self, pyTACSVec tvec):
        '''
        Copy the values x <- vec->x
        '''
        self.this_ptr.copyValues(tvec.this_ptr)

    def zeroEntries(self):
        '''
        Zero all entries in the vector
        '''
        self.this_ptr.zeroEntries()

    def set(self, TacsScalar val):
        '''
        Set all the entries in the vector to val
        '''
        self.this_ptr.set(val)

    def initRand(self):
        '''
        Initialize the random value generator
        '''
        self.this_ptr.initRand()

    def setRand(self, double lower, double upper):
        '''
        Set all the values in the vector using a uniform pseudo-random
        distribution over an interval between lower/upper.
        '''
        self.this_ptr.setRand(lower, upper)

    def placeArray(self, np.ndarray[TacsScalar, ndim=1, mode='c']x):
        '''
        Place an array into the place of the local storage for this vector
        '''
        self.this_ptr.placeArray(<TacsScalar*>x.data)

    def restoreArray(self):
        '''
        Restore the array displaced by the call to placeArray() with the
        original arrayremove the 'placed' array with the original array 
        '''
        self.this_ptr.restoreArray()

    def applyBCs(self):
        '''
        Apply the Dirichlet boundary conditions to the vector
        '''
        self.this_ptr.applyBCs()

    def writeToFile(self, char*filename=''):
        '''
        Write the values to a file.

        This uses MPI file I/O. The filenames must be the same on all
        processors. The format is independent of the number of processors.
    
        The file format is as follows:
        int                       The length of the vector
        len * sizeof(TacsScalar)  The vector entries
        '''
        return self.this_ptr.writeToFile(&filename[0])

    def readFromFile(self, char*filename=''):
        '''
        Read values from a binary data file.

        The size of this vector must be the size of the vector
        originally stored in the file otherwise nothing is read in.

        The file format is as follows:
        int                       The length of the vector
        len * sizeof(TacsScalar)  The vector entries
        '''
        return self.this_ptr.readFromFile(&filename[0])

# "Wrap" the abstract matrix base class TACSMat
cdef class pyTACSMat:
    cdef CyTACSMat *this_ptr
    def __init__(self):
        '''
        The abstract matrix base class for all TACS matrices. All of
        these operations must be defined to be utilized by TACSAssembler.
        '''
        # Create the pointer to the underlying C++ object
        self.this_ptr = new CyTACSMat()
        self.this_ptr.setSelfPointer(<void*>self)
        return

    def __dealloc__(self):
        del self.this_ptr
        return

# Python class for corresponding instance PMat
cdef class pyPMat(pyTACSMat):
    cdef PMat *this_ptr
    def __cinit__(self):
        '''
        The set up for the parallel block-CSR matrix.

        The parallel matrix is split into two parts that are identified
        in the initialization. The diagonal matrix and the off-diagonal 
        matrix. The off-diagonal matrix corresponds to the coupling 
        terms to the external-interface unknowns. The internal-interface
        unknowns must be ordered last on each process. External-interface
        unknowns can only be coupled to other interface-unknowns (either
        external or internal). Thus the global matrix can be represented as
        follows
        
        A_i = [ B_i, F_i ; G_i, C_i ]
        u_i = [ x_i, y_i ]^{T}

        On each process the unknowns are divided into internal-variables x_i, 
        and internal-iterface variables y_i.
  
        Each domain is coupled to other domains only through the interface
        variables y_i.

        A_i u_i + P * E_{ij} y_j = b_i

        where P = [ 0, I_{size(y_i)} ]^{T}

        The matrix structure outlined above can be exploited to achieve 
        efficient and effective parallel preconditioning.
        '''
        self.this_ptr = new PMat()
        
    def __dealloc__(self):
        del self.this_ptr
        
    def zeroEntries(self):
        '''
        Zero the entries in the matrix
        '''
        self.this_ptr.zeroEntries()
    
    def applyBCs(self):
        '''
        Apply the Dirichlet boundary conditions to the matrix
        '''
        self.this_ptr.applyBCs()
    
    def mult(self, pyTACSVec x, pyTACSVec y):
        '''
        Matrix multiplication
        '''
        self.this_ptr.mult(x.this_ptr, y.this_ptr)
    
    def copyValues(self, pyPMat mat):
        '''
        Copy the values from mat         
        '''
        self.this_ptr.copyValues(mat.this_ptr) 
    
    def scale(self, TacsScalar alpha):
        '''
        Scale the entries in the matrix by alpha
        '''
        self.this_ptr.scale(alpha)
    
    def axpy(self, TacsScalar alpha, pyTACSMat mat):
        '''
        Compute y <- y + alpha * x
        '''
        self.this_ptr.axpy(alpha, mat.this_ptr)
    
    def axpby(self, TacsScalar alpha, TacsScalar beta, pyPMat mat):
        '''
        Compute y <- alpha * x + beta * y
        '''
        self.this_ptr.axbpy(alpha, beta, mat.this_ptr)
        
cdef class pyDistMat(pyTACSMat):
    cdef DistMat *this_ptr
    def __cinit__(self):
        self.this_ptr = new DistMat()
    
    def __dealloc__(self):
        del self.this_ptr
        
    def zeroEntries(self):
        '''
        Zero the entries in the matrix
        '''
        self.this_ptr.zeroEntries()

    def applyBCs(self):
        '''
        Apply the Dirichlet boundary conditions to the matrix
        '''
        self.this_ptr.applyBCs()
    
    def mult(self, pyTACSVec x, pyTACSVec y):
        '''
        Matrix multiplication
        '''
        self.this_ptr.mult(x.this_ptr, y.this_ptr)
    
cdef class pyScMat(pyTACSMat):
    cdef ScMat *this_ptr
    def __cinit__(self):
        self.this_ptr = new ScMat()
    
    def __dealloc__(self):
        del self.this_ptr
    
    def zeroEntries(self):
        '''
        Zero the entries in the matrix
        '''
        self.this_ptr.zeroEntries()
    
    def applyBCs(self):
        '''
        Apply the Dirichlet boundary conditions to the matrix
        '''
        self.this_ptr.applyBCs()
    
    def mult(self, pyTACSVec x, pyTACSVec y):
        '''
        Matrix multiplication
        '''
        self.this_ptr.mult(x.this_ptr, y.this_ptr)
    
    def copyValues(self, pyScMat mat):
        '''
        Copy the values from mat         
        '''
        self.this_ptr.copyValues(mat.this_ptr) 
    
    def scale(self, TacsScalar alpha):
        '''
        Scale the entries in the matrix by alpha
        '''
        self.this_ptr.scale(alpha)
    
    def axpy(self, TacsScalar alpha, pyTACSMat mat):
        '''
        Compute y <- y + alpha * x
        '''
        self.this_ptr.axpy(alpha, mat.this_ptr)
    
    def axpby(self, TacsScalar alpha, TacsScalar beta, pyScMat mat):
        '''
        Compute y <- alpha * x + beta * y
        '''
        self.this_ptr.axbpy(alpha, beta, mat.this_ptr)
        
cdef class pyFEMat(pyTACSMat):
    cdef FEMat *this_ptr
    def __cinit__(self):
        '''
        
        '''
        self.this_ptr = new FEMat()

    def __dealloc__(self):
        del self.this_ptr
    
    def zeroEntries(self):
        '''
        Zero the entries in the matrix
        '''
        self.this_ptr.zeroEntries()
    
    def applyBCs(self):
        '''
        Apply the Dirichlet boundary conditions to the matrix
        '''
        self.this_ptr.applyBCs()
    
    def mult(self, pyTACSVec x, pyTACSVec y):
        '''
        Matrix multiplication
        '''
        self.this_ptr.mult(x.this_ptr, y.this_ptr)
        
# Wrap the abstract base class KSMPrint
cdef class pyKSMPrint(pyTACSObject):
    cdef CyKSMPrint *this_ptr
    def __init__(self):
        '''
        KSMPrint base class
        '''
        self.this_ptr = new CyKSMPrint()
        self.this_ptr.setSelfPointer(<void*>self)
            
    def __dealloc__(self):
        del self.this_ptr

cdef class pyKSMPrintStdout(pyKSMPrint):
    cdef pyKSMPrintStdout *this_ptr
    def __cinit__(self, char*descript = '', int rank, int freq):
        '''
        KSMPrint standard print class
        '''
        self.this_ptr = new KSMPrintStdout(&descript[0], rank, freq)

    def __dealloc__(self):
        del self.this_ptr
    
    def printResidual(self, int iterate, TacsScalar res):
        '''
        Print the residual norm for the given iteration on the root processor
    
        input:
        iter:  the iteration count
        res:   the residual normal
        '''
        self.this_ptr.printResidual(iterate, res)
    
# Wrap the abstract pre conditioner class
cdef class pyTACSPc(pyTACSObjec):
    cdef CyTACSPc *this_ptr
    def __init__(self):
        '''
        Abstract Preconditioner base class TACSPc
        '''
        self.this_ptr = new CyTACSPc()
        self.this_ptr.setSelfPointer(<void*>self)
            
    def __dealloc__(self):
        del self.this_ptr
    
# Wrap interpolation class BVecInterp
cdef class pyBVecInterp(pyTACSObject):
    cdef BVecInterp *this_ptr
    def __cinit__(self):
        '''
        Interpolation class BVecInterp
        '''
        self.this_ptr = new BVecInterp()

    def __dealloc__(self):
        del self.this_ptr
    
    def mult(self, pyBVec inVec, pyBvec outVec):
        '''
        Perform the interpolation from inVec to outVec:
    
        Interp*inVec -> outVec
        
        input:
        inVec:  the input vector
    
        output:
        outVec: the interpolated output vector
        '''
        self.this_ptr.mult(inVec.this_ptr, outVec.this_ptr)
        
    def multAdd(self, pyBVec inVec, pyBVec addVec, pyBVec outVec):
        '''
        Perform the interpolation from inVec to outVec:
    
        addVec + Interp*inVec -> outVec
        
        input:
        inVec:  the input vector
        addVec: the vector to add to the output
    
        output:
        outVec: the interpolated output vector
        '''
        self.this_ptr.multAdd(inVec.this_ptr, addVec.this_ptr, outVec.this_ptr)
    
    def multTranspose(self, pyBvec inVec, pyBVec outVec):
        '''
        Perform the interpolation from inVec to outVec:
    
        Interp*inVec -> outVec
    
        input:
        inVec:  the input vector
    
        output:
        outVec: the interpolated output vector
        '''
        self.this_ptr.multTranspose(inVec.this_ptr, outVec.this_ptr)
    
    def multTransposeAdd(self, pyBVec inVec, pyBVec addVec, pyBVec outVec):
        '''
        Perform the interpolation from inVec to outVec:
    
        addVec + Interp*inVec -> outVec
    
        input:
        inVec:  the input vector
        addVec: the vector to add to the output
    
        output:
        outVec: the interpolated output vector
        '''
        self.this_ptr.multTransposeAdd(inVec.this_ptr, addVec.this_ptr, outVec.this_ptr)

cdef class pyAdditiveSchwarz(pyTACSPc):
    cdef AdditivieSchwarz *this_ptr
    def __cinit__(self, pyPMat mat, int levFill, double fill):
        '''
        Additive Schwarz preconditioner
        '''
        self.this_ptr = new AdditiveSchwarz(mat.this_ptr, levFill, fill)
    def __dealloc__(self):
        del self.this_ptr
    
    def factor(self):
        self.this_ptr.factor()
    
    def applyFactor(self, pyTACSVec xvec, pyTACSVec yvec):
        '''
        Apply the preconditioner to the input vector
    
        For the additive Schwarz method that simply involves apply the ILU 
        factorization of the diagonal to the input vector:
    
        y = U^{-1} L^{-1} x
        '''
        self.this_ptr.applyFactor(xvec.this_ptr,yvec.this_ptr)
    
    def applyFactor(self, pyTACSVec xvec):
        '''
        Apply the preconditioner to the input vector
    
        For the additive Schwarz method that simply involves apply the ILU 
        factorization of the diagonal to the input vector:
    
        y = U^{-1} L^{-1} y
        '''
        self.this_ptr.applyFactor(xvec.this_ptr)
    
cdef class pyApproximateSchur(pyTACSPc):
    cdef ApproximateSchur *this_ptr
    def __cinit__(self, pyPMat mat, int levFill, double fill,
                    int inner_gmres_iters, double inner_rtol = 1e-3,
                    double inner_atol = 1e-30):
        '''
        Approximate Schur preconditioner class
        '''
        self.this_ptr = new ApproximateSchur(mat.this_ptr, levFill, fill,
                                                inner_gmres_iters, inner_rtol,
                                                inner_atol) 
    
    def __dealloc__(self):
        del self.this_ptr
    
    def setMonitor(self, pyKSMPrint ksm_print):
        '''
        Monitor the output to the file/screen
        '''
        self.this_ptr.setMonitor(ksm_print.this_ptr)
        
    def factor(self):
       '''
        Factor preconditioner based on the values in the matrix.
        '''
        self.this_ptr.factor()
    
    def applyFactor(self, pyTACSVec xvec, pyTACSVec yvec):
       '''
        Apply the factorization for the approximate Schur method
        '''
        self.this_ptr.applyFactor(xvec.this_ptr, yvec.this_ptr)
    
cdef class pyPcScMat(pyTACSPc):
    cdef PcScMat *this_ptr
    def __cinit__(self, pyScMat smat, int levFill, double fill,
                    int reorder_schur_complement):
        '''
        Global Schur preconditioner
        '''        
        self.this_ptr = new PcScMat(smat.this_ptr, levFill, fill,
                                    reorder_schur_complement)
    def __deallloc__(self):
        del self.this_ptr
        
    def factor(self):
        '''
        Factor preconditioner based on the values in the matrix.
        '''
        self.this_ptr.factor()
    
    def applyFactor(self, pyTACSVec xvec, pyTACSVec yvec):
        '''
        Apply the factorization for the approximate Global Schur method
        '''
        self.this_ptr.applyFactor(xvec.this_ptr, yvec.this_ptr)
    
    def setMonitorFactorFlag(self, int flag):
        '''
        Set the flag that prints out the factorization time
    
        input:
        flag: the flag value for the factor-time monitor
        '''
        self.this_ptr.setMonitorFactorFlag(flag)
    
    def setMonitorBackSolveFlag(self, int flag):
        '''
        Set the flag that prints out the factorization time
    
        input:
        flag: the flag value for the factor-time monitor
        '''
        self.this_ptr.setMonitorBackSolveFlage(flag)
    
    def setAlltoallAssemblyFlag(self, int flag):
        '''
        Set the flag that controls which matrix assembly code to use.
    
        This flag controls whether the Alltoall version of the PDMat matrix
        assembly is used. When true, this uses a faster, but more
        memory-intensive version of the code. Be aware that this can cause
        the code to run out of memory, but can be very beneficial in terms
        of CPU time.  
    
        input:
        flag:  the flag value to use for the Alltoall flag
        '''
       self.this_ptr.setAlltoallAssemblyFlag(flag)
    
cdef class pyTACSMg(pyTACSPc):
    cdef TACSMg *this_ptr
    def __cinit__(self, MPI.Comm comm, int nlevels, double sor_omega
                    int sor_iters, int sor_symmetric):
        '''
        Set up the TACS multi-grid object with a given number of multi-grid
        levels.
    
        Each level is initialized with the given successive over relaxation
        (SOR) factor (sor_omega), and the number of relaxations to perform
        at each level. The flag sor_symmetric indicates whether or not to
        use a symmetric SOR iteration.
    
        input:
        comm:          MPI communicator
        nlevels:       number of multi-grid levels
        sor_omega:     SOR factor
        sor_iters:     number of SOR iterations to perform at each level
        sor_symmetric: symmetric SOR flag
    
        '''
        # Convert the communicator
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.this_ptr = new TACSMg(c_comm, nlevels, sor_omega,
                                    sor_iters, sor_symmetric)
        
    def __dealloc__(self):
        del self.this_ptr
    
    def setLevel(self, int level, pyTACSAssembler tacs, pyBVecInterp restrct,
                 pyBVecInterp interp, int iters = 1):
        '''
        Set the data for the given multi-grid level.
    
        This consists of setting the TACS finite-element model, and the
        restriction and interpolation operators for this level. (These are
        not requried for the lowest level)
    
        input:
        level:     the multigrid level
        tacs:      the TACSAssembler object
        restrict:  the restriction operator
        interp:    the interpolation operator
        iters:     the number of iterations to take at this level
        '''
        self.this_ptr.setLevel(level, tacs.this_ptr, restrct.this_ptr,
                                interp.this_ptr, iters)
    def setVariables(self, pyBVec vec):
        '''
        Set the state variables for this, and all subsequent levels.
    
        input:
        vec:      the input vector of state variables
        '''
        self.this_ptr.setVariables(vec.this_ptr)
    
    def setDesignVars(self, np.ndarray[TacsScalar, ndim=1,mode='c']dvs):
        '''
        Set the design variables for all multi-grid levels.
    
        This call ensures that all the TACSAssembler objects, and the
        objects they reference, share the same set of design variables.
    
        input:
        dvs:     the design variable values
        '''
        # Number of design variables
        ndvs = dvs.shape[0]
        self.setDesignVars(<TacsScalar*>dvs.data, ndvs)
    
    def assembleMatType(self, ElementMatrixType matType,
                        MatrixOrientation matOr):
        '''
        Set up the multi-grid data by computing the matrices at each
        multi-grid level within the problem.
        '''
        self.this_ptr.assembleMatType(matType, matOr)

    def applyFactor(self, pyTACSVec x, pyTACSVec y):
        '''
        Apply the multi-grid preconditioner to try and solve the problem
        x = A^{-1} b

        Assume an initial guess of zero.
        '''
        self.this_ptr.applyFactor(x.this_ptr, y.this_ptr)
    
    def factor(self):
        '''
        Factor the SOR objects for each level of multi-grid, and 
        the direct solver for the lowest level.
        '''
        self.this_ptr.factor()
        
    def solve(self, pyBVec bvec, pyBVec xvec, int max_iters,
                double rtol, double atol):
        '''
        Solve the problem by repeatedly applying the multi-grid method
        '''
        self.this_ptr.solve(bvec.this_ptr, xvec.this_ptr,
                            max_iters, rtol, atol)       
        
    def getMat(self, int level):
        '''
        Retreive the matrix at the specified multi-grid level
        '''
        return self.this_ptr.getMat(level)
    
    def setMonitor(self, pyKSMPrint monitor):
        '''
        Set the monitor to use internally for printing out convergence data.
        '''
        self.this_ptr.setMonitor(monitor.this_ptr)
    
# Wrap the abstract base KSM class
cdef class pyTACSKsm(pyTACSObject):
    cdef CyTACSKsm *this_ptr
    def __init__(self):
        '''
        The abstract Krylov-subspace method class
    
        Solve the linear system A*x = b with a Krylov subspace method,
        possibly using some preconditioner.
        '''
        self.this_ptr = new CyTACSKsm()
        self.this_ptr.setSelfPointer(<void*>self)

    def __dealloc__(self):
        del self.this_ptr

cdef class pyGMRES(pyKsm):
    cdef GMRES *this_ptr
    def __cinit__(self, pyTACSMat mat, pyTACSPc pc, int m, int nrestart, int isFlexible = NULL ):
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
        if isFlexible:
            self.this_ptr.GMRES(mat.this_ptr, pc.this_ptr, m, nrestart, isFlexible)

        else:
            self.this_ptr.GMRES(mat.this_ptr, m, nrestart)

    def solve(self, pyTACSVec b, pyTACSVec x, int zero_guess = 1):
        '''
        Try to solve the linear system using GMRES.
    
        The following code tries to solve the linear system using GMRES (or
        FGMRES if the preconditioner is flexible.)
    
        input:
        b:          the right-hand-side
        x:          the solution vector (with possibly significant entries)
        zero_guess: flag to indicate whether to zero entries of x before solution
        '''
        self.this_ptr.solve(b.this_ptr, x.this_ptr, zero_guess)
    
    def setTolerances(self, double rtol, double atol):
        '''
        Set the relative and absolute tolerances used for the stopping
        criterion.
    
        input:
        rtol: the relative tolerance ||r_k|| < rtol*||r_0||
        atol: the absolute tolerancne ||r_k|| < atol
        '''
        self.this_ptr.setTolerances(rtol, atol)

    def setMonitor(self, pyKSMPrint monitor):
        '''
        Set the object to control how the convergence history is displayed
        (if at all)
    
        input:
        monitor: the KSMPrint monitor object
        '''
        self.this_ptr.setMonitor(monitor.this_ptr)

cdef class pyGCROT(pyTACSKsm):
    cdef GCROT* this_ptr
    def __cinit__(self, pyTACSMat mat, pyTACSPc pc, int outer_iters,
                    int max_outer, int inner_iters, int isFlexible = NULL):
        '''
        Create the GCROT linear system solver
    
        This constructor creates an object that uses a simplified variant of
        GCROT described by Hicken and Zingg with or without a preconditioner.
        
        input:
        mat:        the matrix operator
        pc:         the preconditioner
        outer:      the number of outer vectors
        max_outer:  the maximum number of outer iterations before we give up
        msub:       the size of the underlying GMRES (FGMRES) subspace
        isFlexible: flag to indicate required use of flexible GCROT       
        '''
        # With preconditioner
        if isFlexible:
            self.this_ptr = new GCROT(mat.this_ptr, pc.this_ptr, outer_iters,
                                      max_outer, inner_iters, isFlexible)
        # Without preconditioner
        else:
            self.this_ptr = new GCROT(mat.this_ptr, NULL, outer_iters,
                                      max_outer, inner_iters, 0)
    
    def __dealloc__(self):
        del self.this_ptr

    def solve(self, pyTACSVec b, pyTACSVec x, int zero_guess = 1):
        '''
        Try to solve the linear system using GMRES.
    
        The following code tries to solve the linear system using the
        given (possibly flexible) preconditioner
    
        input:
        b:          the right-hand-side
        x:          the solution vector (with possibly significant entries)
        zero_guess: flag to indicate whether to zero entries of x before solution
        '''
        self.this_ptr.solve(b.this_ptr, x.this_ptr, zero_guess)
    
    def setTolerances(self, double rtol, double atol):
        '''
        Set the relative and absolute tolerances used for the stopping
        criterion.
    
        input:
        rtol: the relative tolerance ||r_k|| < rtol*||r_0||
        atol: the absolute tolerancne ||r_k|| < atol
        '''
        self.this_ptr.setTolerances(rtol, atol)

    def setMonitor(self, pyKSMPrint monitor):
        '''
        Set the object to control how the convergence history is displayed
        (if at all)
    
        input:
        monitor: the KSMPrint monitor object
        '''
        self.this_ptr.setMonitor(monitor.this_ptr)    

cdef class pyPCG(pyTACSKsm):
    cdef PCG *this_ptr
    def __cinit__(self, pyTACSMat mat, pyTACSPc pc, int reset, int nouter):
        '''
        The preconditioned conjugate gradient method
    
        This object is used to perform the preconditioned conjugate gradient
        method for a linear system. The preconditioner cannot be flexible.
    
        input:
        mat:    the matrix operator
        pc:     the preconditioner operator
        reset:  reset the CG iterations every 'reset' iterations
        nouter: the number of resets to try before giving up
        '''
        self.this_ptr = new PCG(mat.this_ptr, pc.this_ptr, reset, nouter)
    
    def __dealloc__(self):
        del self.this_ptr
    
    def solve(self, pyTACSVec b, pyTACSVec x, int zero_guess = 1):
        '''
        Solve the linear system with the preconditioned conjugate gradient
        method
    
        input:
        b:          the right-hand-side
        x:          the solution vector
        zero_guess: flag to indicate whether to start with x = 0
        '''
        self.this_ptr(b.this_ptr, x.this_ptr, zero_guess)

    def setTolerances(self, double rtol, double atol):
        '''
        Set the relative and absolute tolerances used for the stopping
        criterion.
    
        input:
        rtol: the relative tolerance ||r_k|| < rtol*||r_0||
        atol: the absolute tolerancne ||r_k|| < atol
        '''
        self.this_ptr.setTolerances(rtol, atol)

    def setMonitor(self, pyKSMPrint monitor):
        '''
        Set the object to control how the convergence history is displayed
        (if at all)
    
        input:
        monitor: the KSMPrint monitor object
        '''
        self.this_ptr.setMonitor(monitor.this_ptr)
    
# Wrap the TACStoFH5 class
cdef class pyTACSToFH5(pyTACSObject):
    cdef TACSToFH5 *this_ptr
    def __cinit__(self, pyTACSAssembler tacs, ElementType elem_type,
                  int out_type):
        '''
        Create the TACSToFH5 file creation object
        
        input:
        tacs:        the instance of the TACSAssembler object
        elem_type:   the type of element to be used
        out_type:    the output type to write 
        '''
        self.this_ptr = new TACSToFH5(tacs.this_ptr, elem_type, out_type)
    
    def __dealloc__(self):
        del self.this_ptr
    
    def setComponentName( int comp_num, char *group_name):
        '''
        Set the component name for the variable
        '''
        self.this_ptr.setComponentName(comp_num, &group_name[0])
    
    def writeToFile( char*filename):
        '''
        Write the data stored in the TACSAssembler object to filename
        '''
        self.this_ptr.writeToFile(&filename[0])

# Wrap the TACSMeshLoader class
cdef class pyTACSMeshLoader(pyTACSObject):
    cdef TACSMeshLoader *this_ptr
    def __cinit__(self, MPI.Comm comm):
        '''
        This is an interface for reading the NASTRAN-style files i.e. BDF
        '''
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.this_ptr = new TACSMeshLoader(c_comm)
    
    def __dealloc__(self):
        del self.this_ptr

    def scanBdfFile(self, char*filename):
        '''
        This scans a Nastran file - only scanning in information from the
        bulk data section

        The only entries scanned are the entries beginning with elem_types
        and any GRID/GRID* entries
        '''
        self.this_ptr.scanBdfFile(&filename[0])

    def getNumComponents(self):
        '''
        Return the number of components
        '''
        return self.this_ptr.getNumComponents()

    def getComponentDescript(self, int comp_num):
        '''
        Return the component description
        '''
        return self.this_ptr.getComponentDescript(comp_num)

    def getElementDescript(self, int comp_num):
        '''
        Retrieve the element description corresponding to the component number
        '''
        return self.getElementDescript(comp_num)

    def setElement(self, int comp_num, pyTACSElement elem):
        '''
        Set the element associated with a given component number
        '''
        self.this_ptr.setElement(comp_num, elem)

    def getComponentNums(self, np.ndarray[int, ndim=1, mode='c']comp_num)
        '''
        Retrieve the array of elements corresponding to the component
        '''
        # Retrieve the number of elements
        num_elems = comp_num.shape[0]
        self.this_ptr.getComponentNums(<int*>comp_num.data, num_elems)

    def setFunctionDomain(self, pyTACSFunction function,
                          np.ndarray[int, ndim=1, mode='c']comp_nums):
        '''
        Set the function domain

        Given the function, and the set of component numbers that define
        the domain of interest, set the element numbers in the function that
        '''
        self.this_ptr.setFunctionDomain(function.this_ptr, <int*>comp_nums.data)

    def getNumElementsForComps(self, np.ndarray[int, ndim=1, mode='c']numElem):
        '''
        Determine the number of elements that are in each component
        '''
        num_elem_comp = numElem.shape[0]
        self.this_ptr.getNumElementsForComps(<int*>numElem.data, num_elem_comp)

    def getElementComponents(self, np.ndarray[int, ndim=1,mode='c']compIDs):
        '''
        Return the component number for each element
        '''
        num_comp_ids = compIDs.shape[0]
        self.this_ptr.getElementComponents(<int*>compIDs.data, num_comp_ids)
    
    def getTotalNumElements(self):
        '''
        Return the total number of elements
        '''
        return self.this_ptr.getTotalNumElements()
        
    def getOrigNums(self, np.ndarray[double, ndim=1, mode='c']xOrig):
        '''
        Returns the original nodes in the BDF ordering
        '''
        num_xorig = xOrig.shape[0]
        self.this_ptr.getOrigNodes(<double*>xOrig.data, num_xorig)

    def getOrigNodeNums(self, np.ndarray[int, ndim=1, mode='c']nodeNumsOrig):
        '''
        Returns the original node numbers in the BDF ordering
        '''
        n = nodeNumsOrig.shape[0]
        self.this_ptr.getOrigNodeNums(<int*>nodeNumOrig.data, n)

    def getNumNodes(self):
        '''
        Retrieve the number of nodes in the model
        '''
        return self.this_ptr.getNumNodes()

    def createTACSToFH5(self, pyTACSAssembler tacs, ElementType elem_type,
                        int write_flag):
        '''
        Create a TACSToFH5 file creation object
        '''
        return self.this_ptr.createTACSToFH5(tacs.this_ptr, elem_type, write_flag)
    
    def createTACS(self, int varsPerNode, pyTACSAssembler.OrderingType order_type,
                   pyTACSAssembler.MatrixOrdering mat_type):
        '''
        Create a distribtued version of TACS
        '''
        return self.this_ptr.createTACS(varsPerNode, order_type, mat_type)

    def createSerialTACS(self, int split_size, int varsPerNode):
        '''
        Given the split size, create a serial version of TACS
        '''
        return self.this_ptr.createSerialTACS(split_size, varsPerNode)

    def getConnectivityForComp(self, int compID, np.ndarray[int, ndim=1,mode='c']conn):
        '''
        Return the second-order connectivities for elements in component
        'compID'
        '''
        size_conn = conn.shape[0]
        self.this_ptr.getConnectivityForComp(compID, <int*>conn.data, size_conn)
        
    def getNodes(self, np.ndarray[int, ndim=1, mode='c']nodeList,
                 np.ndarray[double, ndim=1,mode='c']pts):
        '''
        This function returns the nodes specified by the indices in nodeList
        '''
        numNodes = nodeList.shape[0]
        npts = pts.shape[0]
        self.this_ptr.getNodes(<int*>nodeList.data, numNodes, <double*>pts.data, npts)
        
    def writeBDF(self, char *filename, np.ndarray[TacsScalar, ndim=1, mode='c']bdf_nodes,
                 np.ndarray[int, ndim=1, mode='c']node_nums):
        '''
        Write a BDF file with the material properties taken from the
        TACSAssembler object
        '''
        num_bdf_nodes = bdf_nodes.shape[0]
        self.this_ptr.writeBDF(&filename[0], <TacsScalar*>bdf_nodes.data,
                               <int*>node_nums.data, num_bdf_nodes)
                               
# TACSMesh class and associated objects
cdef class pyTACSNode(pyTACSObject):
    cdef TACSNode *this_ptr
    def __cinit__(self):
        '''
        '''
        self.this_ptr = new TACSNode()

    def __dealloc__(self):
        del self.this_ptr

    def getLocalVar(self):
        '''
        '''
        return self.this_ptr.getLocalVar()

    def isMPIOwner(self):
        '''
        '''
        return self.this_ptr.isMPIOwner()

    def setNode(self, pyTACSAssembler tacs, np.ndarray[TacsScalar, ndim=1, mode='c']):
        '''
        '''
        self.this_ptr.setNode()

# Linear buckling and frequency analysis in TACS
cdef class pyTACSLinearBuckling(pyTACSObject):
    cdef TACSLinearBuckling *this_ptr
    def __cinit__(self, pyTACSAssembler tacs, TacsScalar sigma, pyTACSMat gmat,
                  pyTACSMat kmat, pyTACSMat aux_mat, pyTACSKsm solver,
                  int max_lanczos_vecs, int num_eigvals, double eig_tol):
        '''
        Linear buckling analysis object.

        This object uses a shift and invert strategy in conjunction with a
        full-orthogonalization Lanczos method to compute the buckling
        eigenvalues and eigenvectors. The method requires objects for the
        stiffness and geometric stiffness matrices as well as an auxiliary
        matrix that is used to store a linear combination of the two. The
        solver must be associated with the auxiliary matrix.
        
        Note: all the matrices supplied must be of the same type and support
        copy/axpy/axpby etc. operations.

        input:
        tacs:         The TACS model corresponding to the analysis problem
        sigma:        The spectral shift
        gmat:         The geometric stiffness matrix
        kmat:         The stiffness matrix
        aux_mat:      The auxiliary matrix associated with the solver
        solver:       Whatever KSM object you create  
        max_lanczos:  Maximum size of the projected subspace
        num_eigvals:  Number of converged eigenvalues required
        eig_tol:      Tolerance of the eigenvalues
        '''
        self.this_ptr = new TACSLinearBuckling(tacs.this_ptr, sigma,
                                               gmat.this_ptr,
                                               kmat.this_ptr,
                                               aux_mat.this_ptr,
                                               solver.this_ptr,
                                               max_lanczos_vecs,
                                               num_eigvals, eig_tol)
    
    def __dealloc__(self):
        del self.this_ptr

    def setSgima(self, TacsScalar sigma):
        '''
        Set the shift value
        '''
        self.this_ptr.setSigma(sigma)

    def solve(self, pyKSMPrint ksm_print = NULL):
        '''
        Solve the linearized buckling problem about x = 0.
  
        This code determine the lowest magnitude eigenvalue such that, 
  
        K x = - lambda * G x

        where K is the stiffness matrix and G is the geometric stiffness
        matrix.  The difficulty is that the geometric stiffness matrix is
        not positive definite necessarily. As a a result, a modified
        eigenvalue problem must be solved instead. In order to attain better
        convergence behaviour from the eigensolver, it is necessary to
        employ a shift and invert strategy to isolate a segment of the
        eigen-spectrum that is of primary interest in the analysis.  For
        buckling problems, the lowest buckling modes must be determined.
        The shift-invert strategy is employed such that the eigenvector
        remains unmodified, while the eigenvalue is shifted.
    
        (K + sigma G) x = - (lambda - sigma) G x = (lambda - sigma)/lambda K x

        The modified problem is then,

        (K + sigma G)^{-1} K x = lambda/(lambda - sigma) x 
        '''
        self.this_ptr.solve(ksm_print.this_ptr)

    def evalEigenDVSens(self, int n,
                        np.ndarray[TacsScalar, ndim=1, mode='c']fdvSens):
        '''
        The function computes the derivatives of the buckling eigenvalues.

        Compute the derivative of the eignevalues w.r.t. the design
        variables. This function must be called after the solve function has
        been called. The stiffness matrix and geometric stiffness matrix
        cannot be modified from the previous call to solve.

        The original eigenvalue problem is
        
        K*u + lambda*G*u = 0

        The derivative of the eigenvalue problem is given as follows:

        d(lambda)/dx = - u^{T}*(dK/dx + lambda*dG/dx)*u/(u^{T}*G*u)

        The difficulty is that the load path is determined by solving an
        auxiliary linear system:

        K*path = f

        Since the geometric stiffness matrix is a function of the path, we
        must compute the total derivative of the inner product of the
        geometric stiffness matrix as follows:

        d(u^{T}*G*u)/dx = [ p(u^{T}*G*u)/px - psi*d(K*path)/dx ]

        where the adjoint variables psi are found by solving the linear
        system:

        K*psi = d(u^{T}*G*u)/d(path)
        '''
        # Retrieve the number of design variables 
        numDVs = fdvSens.shape[0]
        self.this_ptr.evalEigenDVSens(n, <TacsScalar*>fdvSens.data, numDVs)

    def extractEigenvalue(self, int n,
                          np.ndarray[TacsScalar, ndim=1, mode='c']error):
        '''
        Extract the eigenvalue from the analysis.
        '''
        return self.this_ptr.extractEigenvalue(n, <TacsScalar*>error.data)

    def extractEigenvector(self, int n, pyBVec ans,
                           np.ndarray[TacsScalar, ndim=1,mode='c']error):
        '''
        Extract the eigenvector and eigenvalue from the eigenvalue analysis
        
        '''
        return self.this_ptr.extractEigenvector(n, ans.this_ptr,<TacsScalar*>error.data)

    def checkOrthogonality(self):
        '''
        Return ||I - Q^{T}Q ||_{F}
        '''
        return self.this_ptr.checkOrthogonality()
        
cdef class pyTACSFrequencyAnalysis(pyTACSObject):
    cdef TACSFrequencyAnalysis *this_ptr
    def __cinit__(self, pyTACSAssembler tacs, TacsScalar sigma,
                  pyTACSMat mmat, pyTACSMat kmat, pyTACSKsm solver,
                  int max_lanczos, int num_eigvals, double eig_tol):
        '''
        The following code computes the eigenvalues and eigenvectors
        for the natural frequency eigenproblem:

        K u = lambda M u

        The method uses a shift and invert strategy in conjunction with the
        Lanczos method with full orthogonalization.

        Input:
        tacs:        the TACS assembler object
        sigma:       the initial value of the shift-invert
        mmat:        the mass matrix object
        kmat:        the stiffness matrix object
        solver:      the Krylov subspace method associated with the kmat
        max_lanczos: the maximum number of Lanczos vectors to use
        num_eigvals: the number of eigenvalues to use
        eig_tol:     the eigenproblem tolerance
        '''
        self.this_ptr = new TACSFrequencyAnalysis(tacs.this_ptr,
                                                  sigma,
                                                  mmat.this_ptr,
                                                  kmat.this_ptr,
                                                  solver.this_ptr,
                                                  max_lanczos,
                                                  num_eigvals,
                                                  eig_tols)
        
    def __dealloc__(self):
        del self.this_ptr

    def setSgima(self, TacsScalar sigma):
        '''
        Set the shift value
        '''
        self.this_ptr.setSigma(sigma)

    def solve(self, pyKSMPrint ksm_print = NULL):
        '''
        Solve the eigenvalue problem
        '''
        self.this_ptr.solve(ksm_print.this_ptr)

    def evalEigenDVSens(self, int n,
                        np.ndarray[TacsScalar, ndim=1, mode='c']fdvSens):
        '''
        The function computes the derivatives of the buckling eigenvalues.

        Compute the derivative of the eignevalues w.r.t. the design
        variables. This function must be called after the solve function has
        been called. The stiffness matrix and geometric stiffness matrix
        cannot be modified from the previous call to solve.

        The original eigenvalue problem is
        
        K*u = lambda*M*u

        The derivative of the eigenvalue problem is given as follows,

        dK/dx * u + K * du/dx = 
        d lambda/dx M u + lambda * dM/dx * u + lambda M * du/dx

        Since M = M^{T} and K = K^{T}, pre-multiplying by u^{T} gives,

        u^{T} * dK/dx * u = d lambda/dx ( u^{T} * M * u ) + lambda * u^{T} * dM/dx * u

        Rearranging gives, 

        ( u^{T} * M * u ) [ d lambda/dx ] = u^{T} * ( dK/dx - lambda * dM/dx ) * u       
        
        '''
        # Retrieve the number of design variables 
        numDVs = fdvSens.shape[0]
        self.this_ptr.evalEigenDVSens(n, <TacsScalar*>fdvSens.data, numDVs)

    def extractEigenvalue(self, int n,
                          np.ndarray[TacsScalar, ndim=1, mode='c']error):
        '''
        Extract the eigenvalue from the analysis.
        '''
        return self.this_ptr.extractEigenvalue(n, <TacsScalar*>error.data)

    def extractEigenvector(self, int n, pyBVec ans,
                           np.ndarray[TacsScalar, ndim=1,mode='c']error):
        '''
        Extract the eigenvector and eigenvalue from the eigenvalue analysis
        
        '''
        return self.this_ptr.extractEigenvector(n, ans.this_ptr,<TacsScalar*>error.data)

    def checkOrthogonality(self):
        '''
        Return ||I - Q^{T}Q ||_{F}
        '''
        return self.this_ptr.checkOrthogonality()
    
