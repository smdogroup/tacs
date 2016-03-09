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
       self.this_ptr = new TACSObject()
       return
    def __dealloc__(self):
        del self.this_ptr
        return

# Python class for corresponding instance TACSAssembler
cdef class pyTACSAssembler:
    cdef TACSAssembler *this_ptr
    def __cinit__(self, MPI.Comm comm, int numOwnedNodes,
                  int varsPerNode, int numElements, int numNodes,
                  int numDependentNodes, int nodeMaxCSRsize):
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
    
    def addNode(self, int localNodeNum, int tacsNodeNum):
        self.this_ptr.addNode(localNodeNum, tacsNodeNum)

    
    def addNodes(self, np.ndarray[int, ndim=1, mode='c'] localNodeNums, 
                 np.ndarray[int, ndim=1, mode='c'] tacsNodeNums1 =
    NULL):
        '''
        Add a list of tacsNodeNums, no option to add all at once
        '''
        # Get size of localNodeNums
        num_nodes = localNodeNums.shape[0]
        self.this_ptr.addNodes(<int*>localNodeNums.data,
                               <int*>tacsNodeNums1.data,num_nodes)
        
    def setDependentNodes(self, 
                          np.ndarray[int, ndim=1, mode='c']
    depNodeIndex, 
                          np.ndarray[int, ndim=1, mode='c']
    depNodeTolLocal,
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
    
    def applyBCs(self, pyBVec residual, 
                 np.ndarray[TacsScalar, ndim=1, mode='c']local_vars):
        '''
        Apply the boundary conditions to the residual vector
    
        r = u - lambda * u_d
        
        where u_d are the prescribed displacements and lambda is
        the load factor
        '''
        self.this_ptr.applyBCs(residual.this_ptr, 
                               <TacsScalar*>local_vars.data)
        
    def addElement(self, pyTACSElement element, 
                   np.ndarray[int, ndim=1, mode='c']localNodeNums):
        '''
        Add an element to TACS with given node numbers
        '''     
        # Get the number of nodes in element
        num_nodes = localNodeNums.shape[0]
        return self.this_ptr.addElement(element.this_ptr, 
                                        <int*>localNodeNums.data)

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
    def initializeArrays(self):
        '''
        Initialize several arrays that will be used during the 
        course of the analysis
        '''    
        self.this_ptr.initializeArrays()
    
    
    def getDataPointers(self, ):
        '''
        Get pointers to the element data. This code provides a way 
        to automatically segment an array to avoid coding mistakes.
    
        Note that this is coded in such a way that you can provide 
        NULL arguments to
        '''

    
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
    
    
    def setDependentVariables(self, int perNode, 
                              np.ndarray[TacsScalar, ndim=1, mode='c']Vars):
        self.this_ptr.setDependentVariables(perNode, 
                                            <TacsScalar*>Vars.data):
        '''
        Set the dependent variable values based on the independent 
        variable values. This must be called after the local 
        independent nodal values are set. Note that this is a 
        matrix-multiply.
        '''

    
    def addDependentResidual(self, int perNode, 
                             np.ndarray[TacsScalar, ndim=1, mode='c']Vars):
        self.this_ptr.addDependentResidual(perNode, 
                                           <TacsScalar*>Vars.data):
        '''
        Add the residual contribution that is collected from the 
        dependent nodes to the independent nodes that they depend 
        on. 
        Note that this is a transpose matrix-multiply add operation.
        '''
        
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

    def assembleMatType(self, pyTACSMat A, 
                        TacsScalar scaleFactor, 
                        ElementMatrixTypes matType,
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
        self.this_ptr.assembleMatType(A.this_ptr, scaleFactor,
                                      matType, matOr)
        
    def initializeFunctions(self, pyTACSFunction functions,
                            int numFuncs):
        '''
        Initialize a list of functions
        
        Every function must be initialized - usually just once -
        before it can be evaluated. This is handled automatically
        within TACSAssembler.
 
        Check whether the functions are associated with this
        TACSAssembler object.  Next, call preInitalize() for each
        function in the list. Go through the function domain and 
        call initialize for each element in the domain. Finally, 
        call post initialize for each function in the list.

        functions:  an array of function values
        numFuncs:   the number of functions
        '''
        self.this_ptr.initializeFunctions(functions.this_ptr, 
                                          numFuncs)

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

    def evalAdjointResProductsExperimental(self, pyBVec adjoint, 
                                           int numAdjoints, 
                                           np.ndarray[TacsScalar,ndim=1,mode='c']dvSens,
                                           int num_dvs):
        '''
        Evaluate the product of several ajdoint vectors with the
        derivative of the residual w.r.t. the design variables.

        This function is collective on all TACSAssembler
        processes. This computes the product of the derivative of 
        the residual w.r.t. the design variables with several 
        adjoint vectors simultaneously. This saves computational 
        time as the derivative of the element residuals 
        can be reused for each adjoint vector. This function 
        performs the same task as evalAdjointResProduct, but uses 
        more memory than calling it for each adjoint vector.
        
        This implementation is experimental because it is not yet
        fully implemented for all element types, does not handle
        geometric variables and does not yet cover design variables
        in tractions. Nevertheless, the short-term plan is to finish
        this implementation b/c it has the potential to be faster.

        adjoint:     the array of adjoint vectors
        numAdjoints: the number of adjoint vectors
        dvSens:      the product of the derivative of the residuals
                     and the adjoint   
        numDVs:      the number of design variables
        '''
        self.this_ptr.evalAdjointResProductsExperiemental(adjoint.this_ptr,
        numAdjoint, <TacScalar*>dvSens.data, num_dvs)

    def evalMatDVSensInnerProduct(self, TacsScalar alpha, 
                                  ElementMatrixTypes, matType,
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
        
    def evalAdjointResXptSensProducts(self, pyBVec adjoint, 
                                      int numAdjoints,
                                      np.ndarray[TacsScalar,
                                                 ndim=1,mode='c']adjXptSensProduct):
        '''
        Evaluate the product of several ajdoint vectors with the
        derivative of the residual w.r.t. the nodal points.

        This function is collective on all TACSAssembler
        processes. This computes the product of the derivative of 
        the residual w.r.t. the nodal points with several adjoint
        vectors simultaneously. This saves computational time as the
        derivative of the element residuals can be reused for each
        adjoint vector.  

        loadCase:    the load case number
        adjoint:     the array of adjoint vectors
        numAdjoints: the number of adjoint vectors
        dvSens:      the product of the derivative of the residuals
                     and the adjoint 
        numDVs:      the number of design variables
        '''
        self.evalAdjointResXptSensProducts(adjoint.this_ptr,
                                           numAdjoints, 
                                           np.ndarray[TacsScalar,ndim=1,mode='c']
        adjXptSensProduct)
        
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
        self.testElement(elemNum, print_level)

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
        self.testConstitutive(elemNum, print_level)
        
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
    def getOutputDesignVarData(self, np.ndarray[TacsScalar,
                                                ndim=1,mode='c']x, 
                               ElementType elem_type,
                               np.ndarray[double,ndim=1,mode='c']Data,
                               int nvals):
        '''
        Store the design variable data from the file point-wise

        This can produce a large number of entries if there are many
        design variables stored for each element.

        intput:
        x:         the design variables
        numDVs:    the number of design variables/length of array x
        elem_type: the type of element we're going to try to match
        data:      the data array
        nvals:     the number of entries to skip at each node
        '''
        # Get the number of design variables
        num_dvs = x.shape[0]
        
        self.this_ptr.getOutputDesignVarData(<TacsScalar*>x.data,
                                             num_dvs,
                                             elem_type,
                                             <double*>Data.data,
                                             nvals)
cdef class pyTACSVec(pyTACSObject):
    cdef CyTACSVec *this_ptr

    def __cinit__(self):
        self.this_ptr = new CyTACSVec()

    def __dealloc__(self):
        del this_ptr

cdef class pyVarMap(pyTACSObject):
    '''
    Variable map for the parallel distribution of a vector
    '''
    cdef VarMap *this_ptr
    def __cinit__(self, MPI.Comm comm, int N, int bsize):
        # Convert the communicator
        cdef MPI_Comm c_comm = comm.ob_mpi
        self.this_ptr = new VarMap(c_comm, N, bsize)

    def __dealloc(self):
        del this_ptr

cdef class pyBCMap(pyTACSObject):
    
    '''
    Define boundary conditions that are applied after all 
    the matrix/vector values have been set.
    '''
    cdef BCMap *this_ptr
    def __cinit__(self, int num_bcs):
        self.this_ptr = new BCMap(num_bcs)
    def __dealloc__(self):
        del this_ptr

cdef class pyBVec(pyTACSVec):
    '''
    The basic parallel block-based vector class, distributed
    based on the variable map provided. The boundary conditions
    zero the rows.
    '''
    cdef BVec *this_ptr

    def __cinit__(self, pyVarMap rmap, pyBCMap bcs):
        self.this_ptr = new BVec(rmap.this_ptr, bcs.this_ptr)

    def __dealloc__(self):
        del this_ptr

    # Basic vector operations
    def getSize(self, int *size):
        '''
        Get the local size of the vector on this processor
        '''
        return

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
        return self.this_ptr.writeToFile(&checkpoint[0])

    def readFromFile(self, char*filename=''):
        '''
        Read values from a binary data file.

        The size of this vector must be the size of the vector
        originally stored in the file otherwise nothing is read in.

        The file format is as follows:
        int                       The length of the vector
        len * sizeof(TacsScalar)  The vector entries
        '''
        return self.this_ptr.writeToFile(&checkpoint[0])

# "Wrap" the abstract base class TACSMat
cdef class pyTACSMat:
    cdef CyTACSMat *this_ptr

    def __init__(self):
        # Create the pointer to the underlying C++ object
        self.this_ptr = new CyTACSMat()
        self.this_ptr.setSelfPointer(<void*>self)
        self.this_ptr.zeroEntries(_zeroentries)
        self.this_ptr.applyBCs(_applybcs)
        self.this_ptr.mult(_mult)
        return

    def __dealloc__(self):
        del self.this_ptr
        return

# Python class for corresponding instance PMat
cdef class pyPMat:
   cdef PMat *this_ptr

   def __cinit__(self):
       self.this_ptr = new PMat()

   def __dealloc 
