#  This file is part of TACS: The Toolkit for the Analysis of Composite
#  Structures, a parallel finite-element code for structural and
#  multidisciplinary design optimization.
#
#  Copyright (C) 2014 Georgia Tech Research Corporation
#
#  TACS is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this software except in compliance with
#  the License.  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0

# distutils: language=c++

from __future__ import print_function, division

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

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

# Import the element types
PY_ELEMENT_NONE = TACS_ELEMENT_NONE
PY_POINT_ELEMENT = TACS_POINT_ELEMENT
PY_EULER_BEAM = TACS_EULER_BEAM
PY_TIMOSHENKO_BEAM = TACS_TIMOSHENKO_BEAM
PY_PLANE_STRESS = TACS_PLANE_STRESS
PY_SHELL = TACS_SHELL
PY_SOLID = TACS_SOLID
PY_Q3D_ELEMENT = TACS_Q3D_ELEMENT
PY_RIGID = TACS_RIGID
PY_POISSON_2D_ELEMENT = TACS_POISSON_2D_ELEMENT
PY_POISSON_3D_ELEMENT = TACS_POISSON_3D_ELEMENT

# Import the element matrix types
PY_STIFFNESS_MATRIX = STIFFNESS_MATRIX
PY_MASS_MATRIX = MASS_MATRIX
PY_GEOMETRIC_STIFFNESS_MATRIX = GEOMETRIC_STIFFNESS_MATRIX
PY_STIFFNESS_PRODUCT_DERIVATIVE = STIFFNESS_PRODUCT_DERIVATIVE

# Import the orientations
PY_NORMAL = NORMAL
PY_TRANSPOSE = TRANSPOSE

# Import the ordering types
PY_NATURAL_ORDER = NATURAL_ORDER
PY_RCM_ORDER = RCM_ORDER
PY_AMD_ORDER = AMD_ORDER
PY_ND_ORDER = ND_ORDER
PY_TACS_AMD_ORDER = TACS_AMD_ORDER
PY_MULTICOLOR_ORDER = MULTICOLOR_ORDER

# Import the matrix ordering types
PY_ADDITIVE_SCHWARZ = ADDITIVE_SCHWARZ
PY_APPROXIMATE_SCHUR = APPROXIMATE_SCHUR
PY_DIRECT_SCHUR = DIRECT_SCHUR
PY_GAUSS_SEIDEL = GAUSS_SEIDEL

# JDRecycleType
SUM_TWO = JD_SUM_TWO
NUM_RECYCLE = JD_NUM_RECYCLE

# The vector sum/set operations
INSERT_VALUES = TACS_INSERT_VALUES
ADD_VALUES = TACS_ADD_VALUES
INSERT_NONZERO_VALUES = TACS_INSERT_NONZERO_VALUES

# A generic wrapper class for the TACSFunction object
cdef class Function:
    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return
    def setDomain(self, list elem_index):
        cdef int num_elems = len(elem_index)
        cdef int *elem_ind = NULL

        elem_ind = <int*>malloc(num_elems*sizeof(int));
        for i in range(num_elems):
            elem_ind[i] = <int>elem_index[i]

        self.ptr.setDomain(elem_ind, num_elems)
        free(elem_ind)
        return
# A generic wrapper class for the TACSElement object
cdef class Element:
    '''Base element class'''
    def __cinit__(self):
        self.ptr = NULL
        return
    def setComponentNum(self, int comp_num):
        self.ptr.setComponentNum(comp_num)
        return
    def numNodes(self):
        return self.ptr.numNodes()
    def getConstitutive(self):
        return _init_Constitutive(self.ptr.getConstitutive())

# A generic wrapper class for the TACSConstitutive object
cdef class Constitutive:
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL
        return

    def getDVOutputValue(self, int dvIndex,
                         np.ndarray[double, ndim=1, mode='c'] pt):
        return self.ptr.getDVOutputValue(dvIndex, <double*>pt.data)

# This wraps a C++ array with a numpy array for later useage
cdef inplace_array_1d(int nptype, int dim1, void *data_ptr,
                      PyObject *ptr):
    '''Return a numpy version of the array'''
    # Set the shape of the array
    cdef int size = 1
    cdef np.npy_intp shape[1]
    cdef np.ndarray ndarray

    # Set the first entry of the shape array
    shape[0] = <np.npy_intp>dim1

    # Create the array itself - Note that this function will not
    # delete the data once the ndarray goes out of scope
    ndarray = np.PyArray_SimpleNewFromData(size, shape,
                                           nptype, data_ptr)

    # Set the base class who owns the memory
    if ptr != NULL:
        ndarray.base = ptr

    return ndarray

# A generic wrapper for a TACSVec class - usually TACSBVec
cdef class Vec:
    def __cinit__(self, VarMap vmap=None, int bsize=1):
        '''
        A generic wrapper for any of the TACS vector types
        '''
        self.ptr = NULL
        if vmap is not None:
            self.ptr = new TACSBVec(vmap.ptr, bsize)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def getVarsPerNode(self):
        '''
        Return the number of variables per node
        '''
        return self.ptr.getBlockSize()

    def zeroEntries(self):
        '''
        Zero the entries in the matrix
        '''
        self.ptr.zeroEntries()
        return

    def getArray(self):
        '''
        Get the local values
        '''
        cdef TacsScalar *array
        cdef int size = self.ptr.getArray(&array)
        arry = inplace_array_1d(TACS_NPY_SCALAR, size, <void*>array,
                                <PyObject*>self)
        Py_INCREF(self)
        return arry

    def getSize(self):
        '''
        Length of the array
        '''
        cdef int size
        self.ptr.getSize(&size)
        return size

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

    def getValues(self, np.ndarray[int, ndim=1] var):
        '''
        Get the values from the given global indices
        '''
        cdef int fail = 0
        cdef int length = 0
        cdef int bsize = 0
        cdef np.ndarray values
        bsize = self.ptr.getBlockSize()
        length = var.shape[0]
        values = np.zeros(bsize*length)
        fail = self.ptr.getValues(length, <int*>var.data, <TacsScalar*>values.data)
        if fail:
            errmsg = 'Vec: Failed on get values. Incorrect indices'
            raise RuntimeError(errmsg)
        return values

    def setValues(self, np.ndarray[int, ndim=1] var,
                  np.ndarray[TacsScalar, ndim=1] values,
                  TACSBVecOperation op=ADD_VALUES):
        '''
        Set the values into the given vector components.

        Note: Vector indices are global
        '''
        cdef int bsize = 0
        cdef int length = 0
        bsize = self.ptr.getBlockSize()
        if bsize*var.shape[0] != values.shape[0]:
            errmsg = 'Vec: Inconsistent arrays. Must be of size (%d) and (%d)'%(
                var.shape[0], bsize*var.shape[0])
            raise ValueError(errmsg)
        length = var.shape[0]
        self.ptr.setValues(length, <int*>var.data, <TacsScalar*>values.data, op)
        return

    def beginSetValues(self, TACSBVecOperation op=ADD_VALUES):
        '''Begin setting the values: Collective on the TACS communicator'''
        self.ptr.beginSetValues(op)
        return

    def endSetValues(self, TACSBVecOperation op=ADD_VALUES):
        '''Finish setting the values: Collective on the TACS communicator'''
        self.ptr.endSetValues(op)
        return

    def distributeValues(self):
        '''
        Distribute values: Collective on the TACS communicator
        '''
        self.ptr.beginDistributeValues()
        self.ptr.endDistributeValues()
        return

    def beginDistributeValues(self):
        '''
        Begin distributing the values to attain consistent
        local/global entries.  This function is collective on the TACS
        communicator.
        '''
        self.ptr.beginDistributeValues()
        return

    def endDistributeValues(self):
        '''
        Finishe distributing values to attain consistent local/global
        entries. This function is collective on the TACS communicator.
        '''
        self.ptr.endDistributeValues()
        return

    def getVarMap(self):
        '''
        Get the TACSVarMap object from the class
        '''
        return _init_VarMap(self.ptr.getVarMap())

    def writeToFile(self, fname):
        '''
        Write the values to a file.

        This uses MPI file I/O. The filenames must be the same on all
        processors. The format is independent of the number of processors.
        '''
        cdef char *filename = convert_to_chars(fname)
        return self.ptr.writeToFile(filename)

    def readFromFile(self, fname):
        '''
        Read values from a binary data file.

        The size of this vector must be the size of the vector
        originally stored in the file otherwise nothing is read in.
        '''
        cdef char *filename = convert_to_chars(fname)
        return self.ptr.readFromFile(filename)

cdef class VarMap:
    def __cinit__(self, MPI.Comm comm=None, int owned_size=0):
        if comm is None:
            self.ptr = NULL
        else:
            self.ptr = new TACSVarMap(comm.ob_mpi, owned_size)
            self.ptr.incref()

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

cdef class VecInterp:
    def __cinit__(self, inobj=None, outobj=None,
                  int vars_per_node=1):
        self.ptr = NULL
        if isinstance(inobj, VarMap) and isinstance(outobj, VarMap):
            self.ptr = new TACSBVecInterp((<VarMap>inobj).ptr,
                                          (<VarMap>outobj).ptr, vars_per_node)
            self.ptr.incref()
        elif isinstance(inobj, Assembler) and isinstance(outobj, Assembler):
            self.ptr = new TACSBVecInterp((<Assembler>inobj).ptr,
                                          (<Assembler>outobj).ptr)
            self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def initialize(self):
        self.ptr.initialize()
        return

    def mult(self, Vec input_vec, Vec output_vec):
        self.ptr.mult(input_vec.ptr, output_vec.ptr)
        return

    def multAdd(self, Vec input_vec, Vec add_vec, Vec output_vec):
        self.ptr.multAdd(input_vec.ptr, add_vec.ptr, output_vec.ptr)
        return

    def multTranspose(self, Vec input_vec, Vec output_vec):
        self.ptr.multTranspose(input_vec.ptr, output_vec.ptr)
        return

    def multTransposeAdd(self, Vec input_vec, Vec add_vec, Vec output_vec):
        self.ptr.multTransposeAdd(input_vec.ptr, add_vec.ptr, output_vec.ptr)
        return

    def multWeightTranspose(self, Vec input_vec, Vec output_vec):
        self.ptr.multWeightTranspose(input_vec.ptr, output_vec.ptr)
        return

cdef class AuxElements:
    cdef TACSAuxElements *ptr
    def __cinit__(self):
        self.ptr = new TACSAuxElements(100)
        self.ptr.incref()
        return

    def __dealloc__(self):
        self.ptr.decref()
        return

    def addElement(self, int num, Element elem):
        self.ptr.addElement(num, elem.ptr)
        return

cdef class Mat:
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

    def getDenseMatrix(self):
        '''
        Get the dense, column-major order of the matrix.  This only
        works for serial cases.
        '''
        cdef int n = 0
        cdef int m = 0
        cdef int bs = 0
        cdef np.ndarray A = None
        cdef ScMat *sc_ptr = NULL
        cdef BCSRMat *bcsr = NULL
        cdef TACSBVecDistribute *dist = NULL
        cdef TACSBVecIndices *indices = NULL
        cdef const int *indx = NULL

        sc_ptr = _dynamicScMat(self.ptr)
        if sc_ptr != NULL:
            sc_ptr.getBCSRMat(&bcsr, NULL, NULL, NULL)
            bs = bcsr.getBlockSize()
            dist = sc_ptr.getLocalMap()
            indices = dist.getIndices()
            if bcsr != NULL:
                n = bcsr.getRowDim()*bs
                m = bcsr.getColDim()*bs
                A = np.zeros((m, n), dtype=dtype)
                bcsr.getDenseColumnMajor(<TacsScalar*>A.data)

                # Reorder the matrix
                if indices != NULL:
                    indices.getIndices(&indx)
                    P = np.zeros((m, n), dtype=dtype)
                    for j in range(bcsr.getColDim()):
                        for i in range(bcsr.getRowDim()):
                            P[bs*indx[i]:bs*(indx[i]+1), bs*indx[j]:bs*(indx[j]+1)] =\
                                                          A[bs*i:bs*(i+1), bs*j:bs*(j+1)]
                    return P.T
                else:
                    return A.T

        return None

# Create a generic preconditioner class
cdef class Pc:
    def __cinit__(self, Mat mat=None, *args, **kwargs):
        '''
        This creates a default preconditioner depending on the matrix
        type.
        '''
        # Set the defaults for the direct factorization
        cdef int lev_fill = 1000000
        cdef double fill = 10.0
        cdef int reorder = 1
        cdef TACSPMat *p_ptr = NULL
        cdef ScMat *sc_ptr = NULL

        if mat is not None:
            p_ptr = _dynamicPMat(mat.ptr)
            sc_ptr = _dynamicScMat(mat.ptr)

        self.ptr = NULL
        if sc_ptr != NULL:
            self.ptr = new PcScMat(sc_ptr, lev_fill, fill, reorder)
            self.ptr.incref()
        elif p_ptr != NULL:
            self.ptr = new TACSAdditiveSchwarz(p_ptr, 5, 10.0)
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

    def getMat(self):
        '''Retrieve the associated matrix'''
        cdef TACSMat *mat = NULL
        self.ptr.getMat(&mat)
        if mat:
            return _init_Mat(mat)
        return None

    def setMonitorFlags(self, int flag=1):
        '''
        Monitor the time taken in the back-solve
        '''
        cdef PcScMat *pc_ptr = NULL
        pc_ptr = _dynamicPcScMat(self.ptr)
        if pc_ptr is not NULL:
            pc_ptr.setMonitorFactorFlag(flag)
            pc_ptr.setMonitorBackSolveFlag(flag)
        return

cdef class Mg(Pc):
    def __cinit__(self, MPI.Comm comm=None, int num_levs=-1, double omega=0.5,
                  int num_smooth=1, int mg_symm=0):
        if comm is not None and num_levs >= 2:
            self.mg = new TACSMg(comm.ob_mpi, num_levs, omega, num_smooth, mg_symm)
            self.mg.incref()
        else:
            self.mg = NULL
        self.ptr = self.mg

    def setLevel(self, int lev, Assembler assembler, VecInterp interp=None,
                 int num_iters=0, Mat mat=None, Pc pc=None):
        cdef TACSBVecInterp *_interp = NULL
        cdef TACSMat *_mat = NULL
        cdef TACSPc *_pc = NULL
        if interp is not None:
            _interp = interp.ptr
        if mat is not None:
            _mat = mat.ptr
        if pc is not None:
            _pc = pc.ptr
        self.mg.setLevel(lev, assembler.ptr, _interp, num_iters, _mat, _pc)
        return

    def setVariables(self, Vec vec):
        self.mg.setVariables(vec.ptr)

    def assembleJacobian(self, double alpha, double beta, double gamma,
                         Vec residual=None,
                         MatrixOrientation matOr=NORMAL):
        '''Assemble the Jacobian for all levels'''
        cdef TACSBVec *res = NULL
        if residual is not None:
            res = residual.ptr
        self.mg.assembleJacobian(alpha, beta, gamma, res, matOr)
        return

    def assembleMatType(self, ElementMatrixType matType,
                        MatrixOrientation matOr):
        self.mg.assembleMatType(matType, matOr)
        return

    def setMonitor(self, MPI.Comm comm,
                   _descript='GMRES', int freq=10):
        '''
        Set the object to control how the convergence history is displayed
        (if at all)

        input:
        monitor: the KSMPrint monitor object
        '''
        cdef char *descript = convert_to_chars(_descript)
        self.mg.setMonitor(new KSMPrintStdout(descript, comm.rank, freq))

cdef class KSM:
    def __cinit__(self, Mat mat, Pc pc, int m,
                  int nrestart=1, int isFlexible=0):
        '''
        Create a GMRES object for solving a linear system with or
        without a preconditioner.

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

        The following code tries to solve the linear system using GMRES
        (or FGMRES if the preconditioner is flexible.)

        input:
        b:          the right-hand-side
        x:          the solution vector (with possibly significant entries)
        zero_guess:  indicate whether to zero entries of x before solution
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

    def setMonitor(self, MPI.Comm comm,
                   _descript='GMRES', int freq=10):
        '''
        Set the object to control how the convergence history is displayed
        (if at all)

        input:
        monitor: the KSMPrint monitor object
        '''
        cdef char *descript = convert_to_chars(_descript)
        self.ptr.setMonitor(new KSMPrintStdout(descript, comm.rank, freq))

    def setTimeMonitor(self):
        cdef GMRES *gmres_ptr = NULL
        gmres_ptr = _dynamicGMRES(self.ptr)
        if gmres_ptr != NULL:
            gmres_ptr.setTimeMonitor()


cdef class Assembler:
    def __cinit__(self):
        '''
        Constructor for the TACSAssembler object
        '''
        self.ptr = NULL
        return

    @staticmethod
    def create(MPI.Comm comm, int varsPerNode,
               int numOwnedNodes, int numElements,
               int numDependentNodes=0):
        '''
        Static factory method for creating an instance of Assembler
        '''
        cdef MPI_Comm c_comm = comm.ob_mpi
        cdef TACSAssembler *tacs = NULL
        tacs = new TACSAssembler(c_comm, varsPerNode,
                                        numOwnedNodes, numElements,
                                        numDependentNodes)
        return _init_Assembler(tacs)

    def __dealloc__(self):
        '''
        Destructor for Assembler
        '''
        if self.ptr:
            self.ptr.decref()
        return

    def setElementConnectivity(self,
                               np.ndarray[int, ndim=1, mode='c'] conn,
                               np.ndarray[int, ndim=1, mode='c'] ptr):
        '''Set the connectivity'''
        cdef int num_elements = ptr.shape[0]-1
        if num_elements != self.getNumElements():
            raise ValueError('Connectivity must match number of elements')

        # Set the connectivity into TACSAssembler
        self.ptr.setElementConnectivity(<int*>conn.data, <int*>ptr.data)

        return

    def setDependentNodes(self,
                          np.ndarray[int, ndim=1, mode='c'] ptr,
                          np.ndarray[int, ndim=1, mode='c'] conn,
                          np.ndarray[double, ndim=1, mode='c'] weights):
        '''Set the dependent node connectivity'''
        self.ptr.setDependentNodes(<int*>ptr.data, <int*>conn.data,
                                   <double*>weights.data)
        return

    def setElements(self, elements):
        '''Set the elements in to TACSAssembler'''
        if len(elements) != self.getNumElements():
            raise ValueError('Element list must match number of elements')

        # Allocate an array for the element pointers
        cdef TACSElement **elems
        elems = <TACSElement**>malloc(len(elements)*sizeof(TACSElement*))
        if elems is NULL:
            raise MemoryError()

        for i in range(len(elements)):
            elems[i] = (<Element>elements[i]).ptr

        self.ptr.setElements(elems)

        # Free the allocated array
        free(elems)

        return

    def addBCs(self, np.ndarray[int, ndim=1, mode='c'] nodes,
               np.ndarray[int, ndim=1, mode='c'] _vars=None,
               np.ndarray[int, ndim=1, mode='c'] values=None):
        cdef int nnodes = nodes.shape[0]
        cdef int *node_nums = <int*>nodes.data
        cdef int vars_dim = -1
        cdef int *vars_data = NULL
        cdef TacsScalar *values_data = NULL

        # Unwrap the boundary condition information
        if _vars is not None:
            vars_dim = _vars.shape[0]
            vars_data = <int*>_vars.data
        if values is not None:
            values_data = <TacsScalar*>values.data
        self.ptr.addBCs(nnodes, node_nums, vars_dim, vars_data, values_data)
        return

    def addInitBCs(self, np.ndarray[int, ndim=1, mode='c'] nodes,
                   np.ndarray[int, ndim=1, mode='c'] _vars=None,
                   np.ndarray[int, ndim=1, mode='c'] values=None):
        cdef int nnodes = nodes.shape[0]
        cdef int *node_nums = <int*>nodes.data
        cdef int vars_dim = -1
        cdef int *vars_data = NULL
        cdef TacsScalar *values_data = NULL

        # Unwrap the boundary condition information
        if _vars is not None:
            vars_dim = _vars.shape[0]
            vars_data = <int*>_vars.data
        if values is not None:
            values_data = <TacsScalar*>values.data
        self.ptr.addInitBCs(nnodes, node_nums, vars_dim,
                            vars_data, values_data)
        return

    def computeReordering(self, OrderingType order_type,
                          MatrixOrderingType mat_type):
        '''Compute a reordering of the unknowns before initialize()'''
        self.ptr.computeReordering(order_type, mat_type)
        return

    def initialize(self):
        '''
        Function to call after all the nodes and elements have been
        added into the created instance of TACS. This function need not
        be called when tacs is created using TACSCreator class.
        '''
        self.ptr.initialize()
        return

    def getNumNodes(self):
        '''
        Return the number of nodes in the TACSAssembler
        '''
        return self.ptr.getNumNodes()

    def getNumOwnedNodes(self):
        '''
        Return the number of owned nodes
        '''
        return self.ptr.getNumOwnedNodes()

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

    def getVarsPerNode(self):
        '''
        Return the number of variables per node
        '''
        return self.ptr.getVarsPerNode()

    def getOwnerRange(self):
        '''
        Get the ranges of global node numbers owned by each processor
        '''
        cdef MPI_Comm c_comm
        cdef TACSVarMap *varmap = NULL
        cdef const int *owner_range = NULL
        cdef int size = 0
        varmap = self.ptr.getVarMap()
        c_comm = varmap.getMPIComm()
        varmap.getOwnerRange(&owner_range)
        MPI_Comm_size(c_comm, &size)
        rng = np.zeros(size+1, dtype=np.intc)
        for i in range(size+1):
            rng[i] = owner_range[i]
        return rng

    def getElements(self):
        '''Get the elements'''
        # Allocate an array for the element pointers
        cdef int num_elems = 0
        cdef TACSElement **elements
        num_elems = self.ptr.getNumElements()
        elements = self.ptr.getElements()
        e = []
        for i in range(num_elems):
            e.append(_init_Element(elements[i]))
        return e

    def getElementData(self, int num):
        '''Return the element data associated with the element'''
        cdef TACSElement *element = NULL
        cdef int nnodes = 0
        cdef int nvars = 0
        cdef np.ndarray Xpt
        cdef np.ndarray vars0
        cdef np.ndarray dvars
        cdef np.ndarray ddvars

        if num >= 0 and num < self.ptr.getNumElements():
            # Get the element and query the size of the entries
            element = self.ptr.getElement(num, NULL, NULL, NULL, NULL)
            nnodes = element.numNodes()
            nvars = element.numVariables()
            compID = element.getComponentNum()

            # Allocate the numpy array and retrieve the internal data
            Xpt = np.zeros(3*nnodes, dtype=dtype)
            vars0 = np.zeros(nvars, dtype=dtype)
            dvars = np.zeros(nvars, dtype=dtype)
            ddvars = np.zeros(nvars, dtype=dtype)
            self.ptr.getElement(num, <TacsScalar*>Xpt.data,
                                <TacsScalar*>vars0.data, <TacsScalar*>dvars.data,
                                <TacsScalar*>ddvars.data)
        else:
            raise ValueError('Element index out of range')

        return _init_Element(element), Xpt, vars0, dvars, ddvars, compID

    def getElementNodes(self, int num):
        '''Get the node numbers associated with the given element'''
        cdef int num_nodes = 0
        cdef const int *node_nums = NULL
        cdef np.ndarray nodes

        # Get the node numbers
        if num >= 0 and num < self.ptr.getNumElements():
            self.ptr.getElement(num, &node_nums, &num_nodes)
            nodes = np.zeros(num_nodes, dtype=np.int)
            for i in range(num_nodes):
                nodes[i] = node_nums[i]

        return nodes

    def getDesignVars(self,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] dvs):
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

        dvs:     the array of design variable values (output)
        numDVs: the number of design variables
        '''

        # Get number of the design variables
        num_dvs = dvs.shape[0]
        self.ptr.getDesignVars(<TacsScalar*>dvs.data, num_dvs)
        return

    def setDesignVars(self,
                      np.ndarray[TacsScalar, ndim=1, mode='c'] dvs):
        '''
        Set the design variables.

        The design variable values provided must be the same on all
        processes for consistency. This call however, is not
        collective.

        dvs:     the array of design variable values
        numDVs: the number of design variables
        '''
        # Get number of the design variables
        num_dvs = dvs.shape[0]
        self.ptr.setDesignVars(<TacsScalar*>dvs.data, num_dvs)
        return

    def getDesignVarRange(self,
                          np.ndarray[TacsScalar, ndim=1, mode='c'] lb,
                          np.ndarray[TacsScalar, ndim=1, mode='c'] ub):
        '''
        Retrieve the design variable range.

        This call is collective on all TACS processes. The ranges
        provided by indivdual objects may not be consistent (if
        someone provided incorrect data they could be.) Make a
        best guess; take the minimum upper bound and the maximum
        lower bound.

        lowerBound: the lower bound on the design variables (output)
        upperBound: the upper bound on the design variables (output)
        numDVs:      the number of design variables
        '''

        # Get the number of design variables
        num_dvs = lb.shape[0]
        self.ptr.getDesignVarRange(<TacsScalar*>lb.data,
                                   <TacsScalar*>ub.data,
                                   num_dvs)
        return

    def getMPIComm(self):
        '''
        Get the MPI communicator
        '''
        cdef MPI.Comm comm = MPI.Comm()
        comm.ob_mpi = self.ptr.getMPIComm()
        return comm

    def setNumThreads(self, int t):
        '''
        Set the number of threads to use in computation
        '''
        self.ptr.setNumThreads(t)
        return

    def setAuxElements(self, AuxElements elems=None):
        '''Set the auxiliary elements'''
        cdef TACSAuxElements *ptr = NULL
        if elems is not None:
            ptr = elems.ptr
        self.ptr.setAuxElements(ptr)
        return

    def createNodeVec(self):
        '''
        Create a distributed node vector
        '''
        return _init_Vec(self.ptr.createNodeVec())

    def setNodes(self, Vec X):
        '''
        Set the node locations
        '''
        self.ptr.setNodes(X.ptr)
        return

    def getNodes(self, Vec X):
        '''
        Get the node locations
        '''
        self.ptr.getNodes(X.ptr)
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

    def reorderVec(self, Vec vec):
        '''Reorder the vector based on the TACSAssembler reordering'''
        self.ptr.reorderVec(vec.ptr)
        return

    def getReordering(self):
        '''Get the reordering'''
        cdef int size = 0
        cdef np.ndarray oldtonew
        size = self.ptr.getNumOwnedNodes()
        oldtonew = np.zeros(size, dtype=np.intc)
        self.ptr.getReordering(<int*>oldtonew.data)
        return oldtonew

    def applyBCs(self, Vec vec):
        '''Apply boundary conditions to the vector'''
        self.ptr.applyBCs(vec.ptr)
        return

    def applyMatBCs(self, Mat mat):
        '''Apply boundary conditions to the matrix'''
        self.ptr.applyBCs(mat.ptr)
        return

    def setBCs(self, Vec vec):
        '''Apply the Dirichlet boundary conditions to the state vector'''
        self.ptr.setBCs(vec.ptr)
        return

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

    def setSimulationTime(self, double time):
        '''Set the simulation time within TACS'''
        self.ptr.setSimulationTime(time)
        return

    def getSimulationTime(self):
        '''Retrieve the simulation time from TACS'''
        return self.ptr.getSimulationTime()

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

    def setVariables(self, Vec vec=None,
                     Vec dvec=None, Vec ddvec=None):
        '''
        Set the values of the state variables
        '''
        cdef TACSBVec *cvec = NULL
        cdef TACSBVec *cdvec = NULL
        cdef TACSBVec *cddvec = NULL

        if vec is not None:
            cvec = vec.ptr
        if dvec is not None:
            cdvec = dvec.ptr
        if ddvec is not None:
            cddvec = ddvec.ptr

        self.ptr.setVariables(cvec, cdvec, cddvec)
        return

    def getVariables(self, Vec vec=None,
                     Vec dvec=None, Vec ddvec=None):
        '''
        Set the values of the state variables
        '''
        cdef TACSBVec *cvec = NULL
        cdef TACSBVec *cdvec = NULL
        cdef TACSBVec *cddvec = NULL

        if vec is not None:
            cvec = vec.ptr
        if dvec is not None:
            cdvec = dvec.ptr
        if ddvec is not None:
            cddvec = ddvec.ptr

        self.ptr.getVariables(cvec, cdvec, cddvec)
        return

    def copyVariables(self, Vec vec=None,
                      Vec dvec=None, Vec ddvec=None):
        '''
        Set the values of the state variables
        '''
        cdef TACSBVec *cvec = NULL
        cdef TACSBVec *cdvec = NULL
        cdef TACSBVec *cddvec = NULL

        if vec is not None:
            cvec = vec.ptr
        if dvec is not None:
            cdvec = dvec.ptr
        if ddvec is not None:
            cddvec = ddvec.ptr

        self.ptr.copyVariables(cvec, cdvec, cddvec)
        return

    def getInitConditions(self, Vec vec=None,
                          Vec dvec=None, Vec ddvec=None):
        '''
        Retrieve the initial conditions
        '''
        cdef TACSBVec *cvec = NULL
        cdef TACSBVec *cdvec = NULL
        cdef TACSBVec *cddvec = NULL

        if vec is not None:
            cvec = vec.ptr
        if dvec is not None:
            cdvec = dvec.ptr
        if ddvec is not None:
            cddvec = ddvec.ptr

        self.ptr.getInitConditions(cvec, cdvec, cddvec)
        return

    def evalEnergies(self):
        '''Evaluate the kinetic and potential energies'''
        cdef TacsScalar Te, Pe
        self.ptr.evalEnergies(&Te, &Pe)
        return Te, Pe

    def assembleRes(self, Vec residual):
        '''
        Assemble the residual associated with the input load case.

        This residual includes the contributions from element tractions
        set in the TACSSurfaceTraction class and any point loads. Note
        that the vector entries are zeroed first, and that the Dirichlet
        boundary conditions are applied after the assembly of the
        residual is complete.

        rhs:        the residual output
        '''
        self.ptr.assembleRes(residual.ptr)
        return

    def assembleJacobian(self, double alpha, double beta, double gamma,
                         Vec residual, Mat A,
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

        alpha:      coefficient on the variables
        beta:        coefficient on the time-derivative terms
        gamma:      coefficient on the second time derivative term
        residual:  the residual of the governing equations
        A:            the Jacobian matrix
        matOr:      the matrix orientation NORMAL or TRANSPOSE
        '''
        cdef TACSBVec *res = NULL
        if residual is not None:
            res = residual.ptr

        self.ptr.assembleJacobian(alpha, beta, gamma,
                                  res, A.ptr, matOr)
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
        A:            the Jacobian matrix
        alpha:      coefficient on the variables
        beta:        coefficient on the time-derivative terms
        gamma:      coefficient on the second time derivative
        term
        matOr:      the matrix orientation NORMAL or TRANSPOSE
        '''
        self.ptr.assembleMatType(matType, A.ptr, matOr)
        return

    def evalFunctions(self, funclist):
        '''
        Evaluate a list of TACS function
        '''

        # Allocate the array of TACSFunction pointers
        cdef TACSFunction **funcs
        funcs = <TACSFunction**>malloc(len(funclist)*sizeof(TACSFunction*))
        if funcs is NULL:
            raise MemoryError()

        for i in range(len(funclist)):
            funcs[i] = (<Function>funclist[i]).ptr

        # Allocate the numpy array of function values
        cdef np.ndarray fvals = np.zeros(len(funclist), dtype)

        self.ptr.evalFunctions(funcs, len(funclist), <TacsScalar*>fvals.data)

        # Free the allocated array
        free(funcs)

        return fvals

    def evalDVSens(self, func, np.ndarray[TacsScalar, mode='c'] A):
        '''
        Evaluate the derivative of a list of functions w.r.t. the design
        variables.
        '''
        cdef int num_funcs = 1
        cdef TACSFunction **funcs = &((<Function>func).ptr)
        cdef int num_design_vars = A.shape[0]
        cdef TacsScalar *Avals = <TacsScalar*>A.data
        cdef MPI_Comm comm
        cdef int size = 0

        # The the communicator
        comm = self.ptr.getMPIComm()

        # Evaluate the derivative of the functions
        A[:] = 0.0
        self.ptr.addDVSens(1.0, funcs, num_funcs, Avals, num_design_vars)

        # The total size
        size = num_design_vars*num_funcs

        # Allreduce the contributions in-place across all processors
        MPI_Allreduce(MPI_IN_PLACE, Avals, size,
                          TACS_MPI_TYPE, MPI_SUM, comm)
        return

    def evalSVSens(self, Function func, Vec vec):

        '''
        Evaluate the derivative of the function w.r.t. the state
        variables.

        function: the function pointer

        vec:        the derivative of the function w.r.t. the state variables
        '''
        cdef int num_funcs = 1
        cdef TACSFunction **funcs = &((<Function>func).ptr)
        cdef TACSBVec **vecs = &((<Vec>vec).ptr)

        # Compute the derivative
        vec.zeroEntries()
        self.ptr.addSVSens(1.0, 0.0, 0.0, funcs, num_funcs, vecs)
        return

    def evalAdjointResProduct(self, adjoint,
                              np.ndarray[TacsScalar, mode='c'] A):
        '''
        This function is collective on all TACSAssembler processes. This
        computes the product of the derivative of the residual
        w.r.t. the design variables with several adjoint vectors
        simultaneously. This saves computational time as the derivative
        of the element residuals can be reused for each adjoint
        vector. This function performs the same task as
        evalAdjointResProduct, but uses more memory than calling it for
        each adjoint vector.

        adjoint: the array of adjoint vectors
        dvSens: the product of the derivative of the residuals and the adjoint
        num_dvs: the number of design variables
        '''
        cdef int num_adj = 1
        cdef TACSBVec **adj = &((<Vec>adjoint).ptr)
        cdef TacsScalar *Avals = <TacsScalar*>A.data
        cdef int num_design_vars = A.shape[0]
        cdef MPI_Comm comm

        # Get the communicator
        comm = self.ptr.getMPIComm()

        # Add the derivative of the product of the adjoint and residual
        A[:] = 0.0
        self.ptr.addAdjointResProducts(1.0, adj, num_adj,
                                       Avals, num_design_vars)

        # Allreduce the contributions in-place across all processors
        MPI_Allreduce(MPI_IN_PLACE, Avals, num_design_vars*num_adj,
                      TACS_MPI_TYPE, MPI_SUM, comm)
        return

    def evalXptSens(self, Function func, Vec vec):
        '''
        Evaluate the derivative of a list of functions w.r.t.
        the node locations
        '''
        cdef int num_funcs = 1
        cdef TACSFunction **funcs = &((<Function>func).ptr)
        cdef TACSBVec **vecs = &((<Vec>vec).ptr)

        # Zero the entries of the vector
        vec.zeroEntries()

        # Compute the derivative
        self.ptr.addXptSens(1.0, funcs, num_funcs, vecs)

        # Complete the parallel assembly of the sensitivity
        vec.ptr.beginSetValues(TACS_ADD_VALUES)
        vec.ptr.endSetValues(TACS_ADD_VALUES)

        return

    def evalAdjointResXptSensProduct(self, Vec adjoint, Vec prod):
        '''
        This function is collective on all TACSAssembler processes. This
        computes the product of the derivative of the residual
        w.r.t. the node locations with several adjoint vectors
        simultaneously.
        '''
        cdef int num_adj = 1
        cdef TACSBVec **adj = &((<Vec>adjoint).ptr)
        cdef TACSBVec **prods = &((<Vec>prod).ptr)

        # Zero the entries in the product
        prod.zeroEntries()

        # Add the derivative of the product of the adjoint and residual
        self.ptr.addAdjointResXptSensProducts(1.0, adj, num_adj, prods)

        # Add the contributions across all processors
        prod.ptr.beginSetValues(TACS_ADD_VALUES)
        prod.ptr.endSetValues(TACS_ADD_VALUES)

        return

    def addMatDVSensInnerProduct(self, double scale,
                                 ElementMatrixType matType,
                                 Vec psi, Vec phi,
                                 np.ndarray[TacsScalar, mode='c'] A):
        '''
        Add the derivative of the inner product of the specified
        matrix with the input vectors to the design variable
        sensitivity vector A.
        '''
        cdef int num_design_vars = A.shape[0]
        cdef TacsScalar *Avals = <TacsScalar*>A.data
        self.ptr.addMatDVSensInnerProduct(scale, matType, psi.ptr, phi.ptr,
                                          Avals, num_design_vars)
        return

    def evalMatSVSensInnerProduct(self, ElementMatrixType matType,
                                  Vec psi, Vec phi, Vec res):
        self.ptr.evalMatSVSensInnerProduct(matType,
                                           psi.ptr, phi.ptr, res.ptr)
        return

    def addJacobianVecProduct(self, TacsScalar scale,
                              double alpha, double beta, double gamma,
                              Vec x, Vec y, MatrixOrientation matOr=NORMAL):
        '''
        Compute the Jacobian-vector product
        '''
        self.ptr.addJacobianVecProduct(scale, alpha, beta, gamma,
                                       x.ptr, y.ptr, matOr)
        return

    def testElement(self, int elemNum, int print_level,
                    double dh=1e-6, double rtol=1e-8, double atol=1e-1):
        '''
        Test the implementation of the given element number.

        This tests the stiffness matrix and various parts of the
        design-sensitivities: the derivative of the determinant of the
        Jacobian, the derivative of the strain w.r.t. the nodal
        coordinates, and the state variables and the derivative of the
        residual w.r.t. the design variables and nodal coordinates.

        elemNum:      the element number to test
        print_level:  the print level to use
        '''
        self.ptr.testElement(elemNum, print_level, dh, rtol, atol)
        return

    def testConstitutive(self, int elemNum, int print_level):
        '''
        Test the implementation of the given element constitutive
        class.

        This function tests the failure computation and the mass
        sensitivities for the given element.

        elemNum:     the element to retrieve the constitutive object
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

        func:    the function to test
        num_dvs: the number of design variables to use
        dh:      the step size to use
        '''
        self.ptr.testFunction(func.ptr, num_dvs, dh)
        return

# Wrap the TACStoFH5 class
cdef class ToFH5:
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
        tacs:         the instance of the TACSAssembler object
        elem_type:    the type of element to be used
        out_type:     the output type to write
        '''
        self.ptr = new TACSToFH5(tacs.ptr, elem_type, out_type)
        self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def setComponentName(self, int comp_num, _group_name):
        '''
        Set the component name for the variable
        '''
        cdef char *group_name = convert_to_chars(_group_name)
        self.ptr.setComponentName(comp_num, group_name)

    def writeToFile(self, fname):
        '''
        Write the data stored in the TACSAssembler object to filename
        '''
        cdef char *filename = convert_to_chars(fname)
        self.ptr.writeToFile(filename)

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
                              np.ndarray[int, ndim=1, mode='c'] node_ptr,
                              np.ndarray[int, ndim=1, mode='c'] node_conn,
                              np.ndarray[int, ndim=1, mode='c'] id_nums):
        '''Set the connectivity and the element id numbers'''
        cdef int num_elements = node_ptr.shape[0]-1
        if num_elements != id_nums.shape[0]:
            raise ValueError('Connectivity must match number of element ids')
        self.ptr.setGlobalConnectivity(num_nodes, num_elements,
                                       <int*>node_ptr.data,
                                       <int*>node_conn.data,
                                       <int*>id_nums.data)
        return

    def setBoundaryConditions(self,
                              np.ndarray[int, ndim=1, mode='c'] nodes,
                              np.ndarray[int, ndim=1, mode='c'] ptr=None,
                              np.ndarray[int, ndim=1, mode='c'] bcvars=None):
        '''Set the boundary conditions'''
        cdef int num_bcs = nodes.shape[0]
        if ptr is not None and bcvars is not None:
            self.ptr.setBoundaryConditions(num_bcs,
                                           <int*>nodes.data,
                                           <int*>ptr.data,
                                           <int*>bcvars.data)
        else:
            self.ptr.setBoundaryConditions(num_bcs,
                                           <int*>nodes.data,
                                           NULL, NULL)
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

        for i in range(len(elements)):
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

    def getElementPartition(self):
        '''Retrieve the element partition'''
        cdef const int *part = NULL
        cdef int nelems = 0
        nelems = self.ptr.getElementPartition(&part)

        # Create the partition array and return it
        partition = np.zeros(nelems, dtype=np.intc)
        for i in range(nelems):
            partition[i] = part[i]

        # Retrun the copy of the partition array
        return partition

    def createTACS(self):
        return _init_Assembler(self.ptr.createTACS())

    def getTacsNodeNums(self, Assembler assembler,
                        np.ndarray[int, ndim=1, mode='c'] nodes=None):
        cdef int num_orig_nodes = 0
        cdef int *orig_nodes = NULL
        cdef int *new_nodes = NULL
        cdef int num_new_nodes = 0

        # If the array of nodes exists, set the correct pointers/shapes
        if nodes is not None:
            num_orig_nodes = nodes.shape[0]
            orig_nodes = <int*>nodes.data
        self.ptr.getTacsNodeNums(assembler.ptr, orig_nodes, num_orig_nodes,
                                 &new_nodes, &num_new_nodes)

        cdef np.ndarray array = np.zeros(num_new_nodes, dtype=np.int)
        for i in range(num_new_nodes):
            array[i] = new_nodes[i]

        # Free the array from C++
        deleteArray(new_nodes)
        return array

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

    def scanBDFFile(self, fname):
        '''
        This scans a Nastran file - only scanning in information from the
        bulk data section

        The only entries scanned are the entries beginning with elem_types
        and any GRID/GRID* entries
        '''
        cdef char *filename = convert_to_chars(fname)
        self.ptr.scanBDFFile(filename)

    def getNumComponents(self):
        '''
        Return the number of components
        '''
        return self.ptr.getNumComponents()

    def getComponentDescript(self, int comp_num):
        '''
        Return the component description
        '''
        cdef const char *descript
        cdef bytes py_string
        py_string = self.ptr.getComponentDescript(comp_num)
        return convert_bytes_to_str(py_string)

    def getElementDescript(self, int comp_num):
        '''
        Retrieve the element description corresponding to
        the component number
        '''
        cdef const char *descript
        cdef bytes py_string
        py_string = self.ptr.getElementDescript(comp_num)
        return convert_bytes_to_str(py_string)

    def setElement(self, int comp_num, Element elem):
        '''
        Set the element associated with a given component number
        '''
        self.ptr.setElement(comp_num, elem.ptr)

    def setConvertToCoordinate(self, int flag):
        self.ptr.setConvertToCoordinate(flag)
        return

    def getNumNodes(self):
        return self.ptr.getNumNodes()

    def getNumElements(self):
        return self.ptr.getNumElements()

    def createTACS(self, int varsPerNode,
                   OrderingType order_type=NATURAL_ORDER,
                   MatrixOrderingType mat_type=DIRECT_SCHUR):
        '''
        Create a distribtued version of TACS
        '''
        return _init_Assembler(self.ptr.createTACS(varsPerNode,
                                                   order_type, mat_type))

    def addAuxElement(self, AuxElements aux, int comp_num, Element elem):
        '''
        Add the auxiliary element to the given component
        '''
        self.ptr.addAuxElement(aux.ptr, comp_num, elem.ptr)

    def addFunctionDomain(self, Function func, list comp_list):
        '''Add the specified components to the domain of the function'''
        cdef int num_comps = len(comp_list)
        cdef int *comps = NULL
        comps = <int*>malloc(num_comps*sizeof(int))
        for i in range(num_comps):
            comps[i] = <int>comp_list[i]
        self.ptr.addFunctionDomain(func.ptr, comps, num_comps)
        free(comps)
        return

    def getConnectivity(self):
        '''
        Return the connectivity of the mesh
        '''
        cdef int num_nodes
        cdef int num_elements
        cdef const int *elem_ptr
        cdef const int *elem_conn
        cdef const int *elem_comps
        cdef const TacsScalar *Xpts

        self.ptr.getConnectivity(&num_nodes, &num_elements,
                                 &elem_ptr, &elem_conn, &elem_comps, &Xpts)

        cdef np.ndarray ptr = np.zeros(num_elements+1, dtype=np.int)
        if elem_ptr is not NULL:
            for i in range(num_elements+1):
                ptr[i] = elem_ptr[i]

        cdef np.ndarray conn = np.zeros(ptr[-1], dtype=np.int)
        if elem_conn is not NULL:
            for i in range(ptr[-1]):
                conn[i] = elem_conn[i]

        cdef np.ndarray comps = np.zeros(num_elements, dtype=np.int)
        if elem_comps is not NULL:
            for i in range(num_elements):
                comps[i] = elem_comps[i]

        cdef np.ndarray X = np.zeros(3*num_nodes, dtype)
        if Xpts is not NULL:
            for i in range(3*num_nodes):
                X[i] = Xpts[i]

        return ptr, conn, comps, X

    def getBCs(self):
        '''
        Return the boundary conditions associated with the file
        '''
        cdef int num_bcs
        cdef const int *bc_nodes
        cdef const int *bc_vars
        cdef const int *bc_ptr
        cdef const TacsScalar *bc_vals

        self.ptr.getBCs(&num_bcs, &bc_nodes, &bc_vars, &bc_ptr, &bc_vals)

        cdef np.ndarray nodes = np.zeros(num_bcs, dtype=np.int)
        if bc_nodes is not NULL:
            for i in range(num_bcs):
                nodes[i] = bc_nodes[i]

        cdef np.ndarray ptr = np.zeros(num_bcs+1, dtype=np.int)
        if bc_ptr is not NULL:
            for i in range(num_bcs+1):
                ptr[i] = bc_ptr[i]

        cdef np.ndarray bvars = np.zeros(ptr[-1], dtype=np.int)
        cdef np.ndarray vals = np.zeros(ptr[-1], dtype)
        if bc_vars is not NULL and bc_vals is not NULL:
            for i in range(ptr[-1]):
                bvars[i] = bc_vars[i]
                vals[i] = bc_vals[i]

        return nodes, ptr, bvars, vals

cdef class FrequencyAnalysis:
    cdef TACSFrequencyAnalysis *ptr
    def __cinit__(self, Assembler assembler, TacsScalar sigma,
                  Mat M, Mat K, KSM solver, int max_lanczos=100,
                  int num_eigs=5, double eig_tol=1e-6, double eig_rtol=1e-9,
                  Mat PC=None, Pc pc=None, int fgmres_size=5,
                  double eig_atol=1e-30, int num_recycle=0,
                  JDRecycleType recycle_type=JD_NUM_RECYCLE):
        if solver is None:
            self.ptr = new TACSFrequencyAnalysis(assembler.ptr, sigma, M.ptr,
                                                 K.ptr, PC.ptr, pc.ptr, max_lanczos,
                                                 fgmres_size, num_eigs, eig_tol,
                                                 eig_rtol, eig_atol,
                                                 num_recycle, recycle_type)
        else:
            self.ptr = new TACSFrequencyAnalysis(assembler.ptr, sigma, M.ptr,
                                                 K.ptr, solver.ptr, max_lanczos,
                                                 num_eigs, eig_tol)
        self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getSigma(self):
        return self.ptr.getSigma()

    def setSigma(self, TacsScalar sigma):
        self.ptr.setSigma(sigma)

    def solve(self, print_flag=True, int freq=10):
        cdef MPI_Comm comm
        cdef int rank
        cdef TACSAssembler *assembler = NULL
        cdef KSMPrint *ksm_print = NULL

        if print_flag:
            assembler = self.ptr.getTACS()
            comm = assembler.getMPIComm()
            MPI_Comm_rank(comm, &rank)
            ksm_print = new KSMPrintStdout("FrequencyAnalysis", rank, freq)

        self.ptr.solve(ksm_print)
        return

    def extractEigenvalue(self, int eig):
        cdef TacsScalar err = 0.0
        cdef TacsScalar eigval = 0.0
        eigval = self.ptr.extractEigenvalue(eig, &err)
        return eigval, err

    def extractEigenvector(self, int eig, Vec vec):
        cdef TacsScalar err = 0.0
        cdef TacsScalar eigval = 0.0
        eigval = self.ptr.extractEigenvector(eig, vec.ptr, &err)
        return eigval, err

cdef class BucklingAnalysis:
    cdef TACSLinearBuckling *ptr
    def __cinit__(self, Assembler assembler, TacsScalar sigma,
                  Mat G, Mat K, KSM solver, int max_lanczos=100,
                  int num_eigs=5, double eig_tol=1e-6):
        # Get the auxiliary matrix from the solver
        cdef TACSMat *aux_mat
        solver.ptr.getOperators(&aux_mat, NULL)

        # Create the linear buckling class
        self.ptr = new TACSLinearBuckling(assembler.ptr, sigma, G.ptr,
                                          K.ptr, aux_mat, solver.ptr, max_lanczos,
                                          num_eigs, eig_tol)
        self.ptr.incref()
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()

    def getSigma(self):
        return self.ptr.getSigma()

    def setSigma(self, TacsScalar sigma):
        self.ptr.setSigma(sigma)

    def solve(self, Vec force=None, print_flag=True, int freq=10):
        cdef TACSBVec *f = NULL
        cdef MPI_Comm comm
        cdef int rank
        cdef TACSAssembler *assembler = NULL
        cdef KSMPrint *ksm_print = NULL

        if force is not None:
            f = force.ptr

        if print_flag:
            assembler = self.ptr.getTACS()
            comm = assembler.getMPIComm()
            MPI_Comm_rank(comm, &rank)
            ksm_print = new KSMPrintStdout("BucklingAnalysis", rank, freq)

        self.ptr.solve(f, ksm_print)
        return

    def extractEigenvalue(self, int eig):
        cdef TacsScalar err = 0.0
        cdef TacsScalar eigval = 0.0
        eigval = self.ptr.extractEigenvalue(eig, &err)
        return eigval, err

    def extractEigenvector(self, int eig, Vec vec):
        cdef TacsScalar err = 0.0
        cdef TacsScalar eigval = 0.0
        eigval = self.ptr.extractEigenvector(eig, vec.ptr, &err)
        return eigval, err

# A generic abstract class for all integrators implemented in TACS
cdef class Integrator:
    '''
    Class containing functions for solving the equations forward in
    time and adjoint.
    '''
    cdef TACSIntegrator *ptr

    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr:
            self.ptr.decref()
        return

    def setRelTol(self, double rtol):
        '''
        Relative tolerance of Newton solver
        '''
        self.ptr.setRelTol(rtol)
        return

    def setAbsTol(self, double atol):
        '''
        Absolute tolerance of Newton solver
        '''
        self.ptr.setAbsTol(atol)
        return

    def setMaxNewtonIters(self, int max_newton_iters):
        '''
        Maximum iteration in Newton solver
        '''
        self.ptr.setMaxNewtonIters(max_newton_iters)
        return

    def setPrintLevel(self, int print_level, fname=None):
        '''
        Level of print from TACSIntegrator
        0: off
        1: summary each step
        2: summary each newton iteration
        '''
        cdef char *filename = NULL
        if fname is not None:
            filename = convert_to_chars(fname)
        self.ptr.setPrintLevel(print_level, filename)
        return

    def setJacAssemblyFreq(self, int freq):
        '''
        How frequent to assemble the Jacobian for nonlinear solve
        '''
        self.ptr.setJacAssemblyFreq(freq)
        return

    def setUseLapack(self, int use_lapack):
        '''
        Should TACSIntegrator use lapack for linear solve. This will
        be slow and need to be in serial mode only
        '''
        self.ptr.setUseLapack(use_lapack)
        return

    def setUseFEMat(self, int use_femat, OrderingType order_type):
        '''
        Set FEMAT = 1 for parallel execution
        '''
        self.ptr.setUseFEMat(use_femat, order_type)
        return

    def setInitNewtonDeltaFraction(self, double frac):
        '''
        Parameter for globalization in Newton solver
        '''
        self.ptr.setInitNewtonDeltaFraction(frac)
        return

    def setKrylovSubspaceMethod(self, KSM ksm):
        '''
        Make TACS use this linear solver
        '''
        self.ptr.setKrylovSubspaceMethod(ksm.ptr)
        return

    def setTimeInterval(self, double tinit, double tfinal):
        '''
        Set the time interval for the simulation
        '''
        self.ptr.setTimeInterval(tinit, tfinal)
        return

    def setFunctions(self, list funcs, int num_dvs,
                     int start_plane=-1, int end_plane=-1):
        '''
        Sets the functions for obtaining the derivatives.
        '''
        # Allocate the array of TACSFunction pointers
        cdef TACSFunction **fn = NULL
        cdef int nfuncs = 0
        nfuncs = len(funcs)
        fn = <TACSFunction**>malloc(nfuncs*sizeof(TACSFunction*))
        for i in range(nfuncs):
            if funcs[i] is None:
                fn[i] = NULL
            else:
                fn[i] = (<Function>funcs[i]).ptr
        self.ptr.setFunctions(fn, nfuncs, num_dvs, start_plane, end_plane)
        free(fn)
        return

    def lapackNaturalFrequencies(self, Vec q, Vec qdot, Vec qddot,
                                 int use_gyroscopic=1):
        cdef int size = 0
        cdef np.ndarray eigvals
        size = q.ptr.getArray(NULL)
        eigvals = np.zeros(size)
        self.ptr.lapackNaturalFrequencies(use_gyroscopic,
                                          q.ptr, qdot.ptr, qddot.ptr,
                                          <TacsScalar*>eigvals.data,
                                          NULL)
        return eigvals

    def lapackNaturalModes(self, Vec q, Vec qdot, Vec qddot,
                           int use_gyroscopic=1):
        cdef int size = 0
        cdef np.ndarray eigvals
        size = q.ptr.getArray(NULL)
        eigvals = np.zeros(size)
        cdef np.ndarray modes
        modes = np.zeros((size,size))
        self.ptr.lapackNaturalFrequencies(use_gyroscopic,
                                          q.ptr, qdot.ptr, qddot.ptr,
                                          <TacsScalar*>eigvals.data,
                                          <TacsScalar*>modes.data)
        return eigvals, modes

    def iterate(self, int step_num, Vec forces=None):
        '''
        Solve the nonlinear system at current time step
        '''
        cdef TACSBVec *fvec = NULL
        if forces is not None:
            fvec = forces.ptr
        return self.ptr.iterate(step_num, fvec)

    def integrate(self):
        '''
        Integrates the governing equations forward in time
        '''
        return self.ptr.integrate()

    def evalFunctions(self, funclist):
        '''
        Evaluate a list of TACS function in time
        '''

        # Allocate the array of TACSFunction pointers
        cdef TACSFunction **funcs
        funcs = <TACSFunction**>malloc(len(funclist)*sizeof(TACSFunction*))
        for i in range(len(funclist)):
            if funclist[i] is None:
                funcs[i] = NULL
            else:
                funcs[i] = (<Function>funclist[i]).ptr

        # Allocate the numpy array of function values
        cdef np.ndarray fvals = np.zeros(len(funclist), dtype)
        self.ptr.evalFunctions(<TacsScalar*>fvals.data)

        # Free the allocated array
        free(funcs)

        return fvals

    def initAdjoint(self, int step_num):
        '''
        Initialize adjoint step
        '''
        self.ptr.initAdjoint(step_num)
        return

    def iterateAdjoint(self, int step_num, list adjlist=None):
        '''
        Perform one iteration in reverse mode
        '''
        cdef TACSBVec **adjoint = NULL
        cdef int nadj = 0
        if adjlist is not None:
            nadj = len(adjlist)
            adjoint = <TACSBVec**>malloc(nadj*sizeof(TACSBVec*))

            for i in range(nadj):
                adjoint[i] = (<Vec>adjlist[i]).ptr
        self.ptr.iterateAdjoint(step_num, adjoint)
        if adjlist is not None:
            free(adjoint)
        return

    def postAdjoint(self, int step_num):
        '''
        Terminate adjoint step
        '''
        self.ptr.postAdjoint(step_num)
        return

    def integrateAdjoint(self):
        '''
        Integrates the adjoint backwards in time
        '''
        self.ptr.integrateAdjoint()
        return

    def getAdjoint(self, int step_num, int func_num):
        '''
        Get the adjoint vector at the given step
        '''
        cdef TACSBVec *adjoint = NULL
        self.ptr.getAdjoint(step_num, func_num, &adjoint)
        return _init_Vec(adjoint)

    def getGradient(self, np.ndarray[TacsScalar, ndim=1, mode='c'] dfdx):
        '''
        Get the time-dependent derivative of functionals
        '''
        self.ptr.getGradient(<TacsScalar*>dfdx.data)
        return

    def getXptGradient(self, int func_num):
        '''
        Get the time-dependent nodal derivatives of the functional
        '''
        cdef TACSBVec *dfdXpt = NULL
        self.ptr.getXptGradient(func_num, &dfdXpt)
        return _init_Vec(dfdXpt)

    def getStates(self, int time_step):
        '''
        TACS state vectors are returned at the given time step
        '''
        cdef double time
        cdef TACSBVec *cq = NULL
        cdef TACSBVec *cqdot = NULL
        cdef TACSBVec *cqddot = NULL
        time = self.ptr.getStates(time_step, &cq, &cqdot, &cqddot)
        return time, _init_Vec(cq), _init_Vec(cqdot), _init_Vec(cqddot)

    def checkGradients(self, double dh):
        '''
        Performs a FD/CSD verification of the gradients
        '''
        self.ptr.checkGradients(dh)
        return

    def setOutputPrefix(self, _prefix):
        '''
        Output directory to use for f5 files
        '''
        cdef char *prefix = convert_to_chars(_prefix)
        self.ptr.setOutputPrefix(prefix)
        return

    def setOutputFrequency(self, int write_freq=0):
        '''
        Configure how frequent to write f5 files
        '''
        self.ptr.setOutputFrequency(write_freq)
        return

    def setRigidOutput(self, ToFH5 f5):
        '''
        Configure the export of rigid bodies
        '''
        self.ptr.setRigidOutput(f5.ptr)
        return

    def setShellOutput(self, ToFH5 f5):
        '''
        Configure the export of shell bodies
        '''
        self.ptr.setShellOutput(f5.ptr)
        return

    def setBeamOutput(self, ToFH5 f5):
        '''
        Configure the export of beam elements
        '''
        self.ptr.setBeamOutput(f5.ptr)
        return

    def setSolidOutput(self, ToFH5 f5):
        '''
        Configure the export of solid elements
        '''
        self.ptr.setSolidOutput(f5.ptr)
        return

    def getNumTimeSteps(self):
        '''
        Get the number of time steps
        '''
        return self.ptr.getNumTimeSteps()

    def writeRawSolution(self, fname, int format_flag=2):
        cdef char *filename = convert_to_chars(fname)
        self.ptr.writeRawSolution(filename, format_flag)
        return

    def persistStates(self, int step_num, prefix=''):
        '''
        Writes the states variables to disk. The string argument
        prefix can be used to put the binaries in a separate
        directory.
        '''
        # Get the current state varibles from integrator
        t, q, qd, qdd = self.getStates(step_num)

        # Make filenames for each state vector
        qfnametmp = '%sq-%d.bin' % (prefix, step_num)
        qdfnametmp = '%sqd-%d.bin' % (prefix, step_num)
        qddfnametmp = '%sqdd-%d.bin' % (prefix, step_num)
        cdef char *qfname = convert_to_chars(qfnametmp)
        cdef char *qdfname = convert_to_chars(qdfnametmp)
        cdef char *qddfname = convert_to_chars(qddfnametmp)

        # Write states to disk
        flag1 = q.writeToFile(qfname)
        flag2 = qd.writeToFile(qdfname)
        flag3 = qdd.writeToFile(qddfname)
        flag = max(flag1, flag2, flag3)

        return flag

    def loadStates(self, int step_num, prefix=''):
        '''
        Loads the states variables to disk. The string argument
        prefix can be used to put the binaries in a separate
        directory.
        '''
        # Make filenames for each state vector
        qfnametmp = '%sq-%d.bin' % (prefix, step_num)
        qdfnametmp = '%sqd-%d.bin' % (prefix, step_num)
        qddfnametmp = '%sqdd-%d.bin' % (prefix, step_num)
        cdef char *qfname = convert_to_chars(qfnametmp)
        cdef char *qdfname = convert_to_chars(qdfnametmp)
        cdef char *qddfname = convert_to_chars(qddfnametmp)

        # Get the current state varibles from integrator
        t, q, qd, qdd = self.getStates(step_num)

        # Store values read from file
        return max(q.readFromFile(qfname),
                   qd.readFromFile(qdfname),
                   qdd.readFromFile(qddfname))

cdef class BDFIntegrator(Integrator):
    '''
    Backward-Difference method for integration. This currently
    supports upto third order accuracy in time integration.
    '''
    def __cinit__(self, Assembler tacs,
                  double tinit, double tfinal,
                  double num_steps,
                  int max_bdf_order):
        '''
        Constructor for BDF Integrators of order 1, 2 and 3
        '''
        self.ptr = new TACSBDFIntegrator(tacs.ptr, tinit, tfinal,
                                         num_steps, max_bdf_order)
        self.ptr.incref()
        return

cdef class DIRKIntegrator(Integrator):
    '''
    Diagonally-Implicit-Runge-Kutta integration class. This supports
    upto fourth order accuracy in time and domain. One stage DIRK is
    second order accurate, two stage DIRK is third order accurate and
    '''
    def __cinit__(self, Assembler tacs,
                  double tinit, double tfinal,
                  double num_steps,
                  int stages):
        self.ptr = new TACSDIRKIntegrator(tacs.ptr, tinit, tfinal,
                                          num_steps, stages)
        self.ptr.incref()
        return

cdef class ABMIntegrator(Integrator):
    '''
    Adams-Bashforth-Moulton method for integration. This currently
    supports upto sixth order accuracy in time integration.
    '''
    def __cinit__(self, Assembler tacs,
                  double tinit, double tfinal,
                  double num_steps,
                  int max_abm_order):
        '''
        Constructor for ABM Integrators of order 1, 2, 3, 4, 5 and 6
        '''
        self.ptr = new TACSBDFIntegrator(tacs.ptr, tinit, tfinal,
                                         num_steps, max_abm_order)
        self.ptr.incref()
        return

cdef class NBGIntegrator(Integrator):
    '''
    Newmark-Beta-Gamma method for integration.
    '''
    def __cinit__(self, Assembler tacs,
                      double tinit, double tfinal,
                      double num_steps,
                      int order):
        '''
        Constructor for Newmark-Beta-Gamma method of integration
        '''
        self.ptr = new TACSBDFIntegrator(tacs.ptr, tinit, tfinal,
                                         num_steps, order)
        self.ptr.incref()
        return
