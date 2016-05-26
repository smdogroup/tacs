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

# Import the definitions
from TACS cimport *
from functions cimport *

cdef extern from "mpi-compat.h":
   pass

# A generic wrapper class for the TACSFunction object
cdef class Function:
   def __cinit__(self):
      self.ptr = NULL
      return
   def __dealloc__(self):
      if self.ptr:
         self.ptr.decref()
      return

cdef class compliance(Function):
    def __cinit__(self, Assembler tacs):
        '''
        Wrap the function Compliance
        '''
        self.ptr = new Compliance(tacs.ptr)
        self.ptr.incref()
        return
    
cdef class mass(Function):
    def __cinit__(self, Assembler tacs):
        '''
        Wrap the function StructuralMass
        '''
        self.ptr = new StructuralMass(tacs.ptr)
        self.ptr.incref()
        return

cdef class ksfailure(Function):
    def __cinit__(self, Assembler tacs, double ksWeight, double alpha=1.0):
        '''
        Wrap the function KSFailure
        '''
        self.ptr = new KSFailure(tacs.ptr, ksWeight, alpha)
        self.ptr.incref()
        return
    
cdef class ksbuckling(Function):
    def __cinit__(self, Assembler tacs, double ksWeight):
        '''
        Wrap the function KSBuckling
        '''
        self.ptr = new KSBuckling(tacs.ptr, ksWeight)
        self.ptr.incref()
        return

cdef class ksdisplacement(Function):
    def __cinit__(self, Assembler tacs, np.ndarray[TacsScalar, ndim=1, mode='c']d,
                  double ksWeight, double alpha=1.0):
        '''
        Wrap the function KSDisplacement
        '''
        self.ptr = new KSDisplacement(tacs.ptr, <TacsScalar*>d.data, ksWeight, alpha)
        self.ptr.incref()
        return

cdef class inducedbuckling(Function):
    def __cinit__(self, Assembler tacs, np.ndarray[int, ndim=1, mode='c']elemNums, double P):
        '''
        Wrap the function InducedBuckling
        '''
        self.ptr = new InducedBuckling(tacs.ptr, <int*>elemNums.data,
                                       elemNums.shape[0], P)
        self.ptr.incref()
        return

cdef class inducedfailure(Function):
    def __cinit__(self, Assembler tacs, np.ndarray[int, ndim=1, mode='c']elemNums, double P):
        '''
        Wrap the function InducedFailure
        '''
        self.ptr = new InducedFailure(tacs.ptr, <int*>elemNums.data,
                                      elemNums.shape[0], P)
        self.ptr.incref()
        return
