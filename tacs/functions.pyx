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

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
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

cdef class Compliance(Function):
    def __cinit__(self, Assembler tacs):
        '''
        Wrap the function Compliance
        '''
        self.ptr = new TACSCompliance(tacs.ptr)
        self.ptr.incref()
        return
    
cdef class StructuralMass(Function):
    def __cinit__(self, Assembler tacs):
        '''
        Wrap the function StructuralMass
        '''
        self.ptr = new TACSStructuralMass(tacs.ptr)
        self.ptr.incref()
        return

cdef class KSFailure(Function):
    def __cinit__(self, Assembler tacs, double ksWeight, double alpha=1.0):
        '''
        Wrap the function KSFailure
        '''
        self.ptr = new TACSKSFailure(tacs.ptr, ksWeight,
                                     KS_FAILURE, alpha)
        self.ptr.incref()
        return

cdef class InducedFailure(Function):
    def __cinit__(self, Assembler tacs, double P):
        '''
        Wrap the function InducedFailure
        '''
        self.ptr = new TACSInducedFailure(tacs.ptr, P,
                                          INDUCED_FAILURE)
        self.ptr.incref()
        return
