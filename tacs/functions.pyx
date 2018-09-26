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
    cdef TACSKSFailure *ksptr
    def __cinit__(self, Assembler tacs, double ksWeight, double alpha=1.0):
        '''
        Wrap the function KSFailure
        '''
        self.ksptr = new TACSKSFailure(tacs.ptr, ksWeight,
                                       KS_FAILURE, alpha)
        self.ptr = self.ksptr
        self.ptr.incref()
        return

    def setKSFailureType(self, ftype='discrete'):
        if ftype == 'discrete':
            self.ksptr.setKSFailureType(KS_FAILURE_DISCRETE)
        elif ftype == 'continuous':
            self.ksptr.setKSFailureType(KS_FAILURE_CONTINUOUS)
        elif ftype == 'pnorm-discrete':
            self.ksptr.setKSFailureType(PNORM_FAILURE_DISCRETE)
        elif ftype == 'pnorm-continuous':
            self.ksptr.setKSFailureType(PNORM_FAILURE_CONTINUOUS)
        return

cdef class KSDisplacement(Function):
    cdef TACSKSDisplacement *ksptr
    def __cinit__(self, Assembler tacs, double ksWeight, dirs):
        '''
        Wrap the function KSFailure
        '''
        cdef TacsScalar _dirs[3]
        _dirs[0] = dirs[0]
        _dirs[1] = dirs[1]
        _dirs[2] = dirs[2]
        self.ksptr = new TACSKSDisplacement(tacs.ptr, ksWeight, _dirs)
        self.ptr = self.ksptr
        self.ptr.incref()
        return

    def setKSDispType(self, ftype='discrete'):
        if ftype == 'discrete':
            self.ksptr.setKSDispType(KS_DISP_DISCRETE)
        elif ftype == 'continuous':
            self.ksptr.setKSDispType(KS_DISP_CONTINUOUS)
        elif ftype == 'pnorm-discrete':
            self.ksptr.setKSDispType(PNORM_DISP_DISCRETE)
        elif ftype == 'pnorm-continuous':
            self.ksptr.setKSDispType(PNORM_DISP_CONTINUOUS)
        return

cdef class DisplacementIntegral(Function):
    cdef TACSDisplacementIntegral *dptr
    def __cinit__(self, Assembler tacs, dirs):
        '''
        Wrap the function KSFailure
        '''
        cdef TacsScalar _dirs[3]
        _dirs[0] = dirs[0]
        _dirs[1] = dirs[1]
        _dirs[2] = dirs[2]
        self.dptr = new TACSDisplacementIntegral(tacs.ptr, _dirs)
        self.ptr = self.dptr
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
