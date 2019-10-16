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
from constitutive cimport *

# Include the definitions
include "TacsDefs.pxi"

# Include the mpi4py header
cdef extern from "mpi-compat.h":
    pass

cdef class Constitutive:
    def __cinit__(self, *args, **kwargs):
        self.ptr = NULL

    def __dealloc__(self):
        if self.ptr != NULL:
            self.ptr.decref()

    def getNumStresses(self):
        if self.ptr != NULL:
            return self.ptr.getNumStresses()
        return 0


cdef class MaterialProperties:
    cdef TACSMaterialProperties *ptr
    def __cinit__(self, **kwargs):
        cdef TacsScalar rho = 2700.0
        cdef TacsScalar specific_heat = 921.0
        cdef TacsScalar E = 70e9
        cdef TacsScalar nu = 0.3
        cdef TacsScalar ys = 270e6
        cdef TacsScalar alpha = 24e-6
        cdef TacsScalar kappa = 230.0

        if 'rho' in kwargs:
            rho = kwargs['rho']
        if 'specific_heat' in kwargs:
            specific_heat = kwargs['specific_heat']
        if 'E' in kwargs:
            E = kwargs['E']
        if 'nu' in kwargs:
            nu = kwargs['nu']
        if 'ys' in kwargs:
            ys = kwargs['ys']
        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
        if 'kappa' in kwargs:
            kappa = kwargs['kappa']

        self.ptr = new TACSMaterialProperties(rho, specific_heat,
                                              E, nu, ys, alpha, kappa)
        self.ptr.incref()

    def __dealloc__(self):
        self.ptr.decref()


cdef class PlaneStressConstitutive(Constitutive):
    def __cinit__(self, MaterialProperties props, **kwargs):
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        self.cptr = new TACSPlaneStressConstitutive(props.ptr, t, tNum,
                                                    tlb, tub)
        self.ptr = self.cptr
        self.ptr.incref()


cdef class SolidConstitutive(Constitutive):
    def __cinit__(self, MaterialProperties props, **kwargs):
        cdef TacsScalar t = 1.0
        cdef int tNum = -1
        cdef TacsScalar tlb = 0.0
        cdef TacsScalar tub = 10.0

        if 't' in kwargs:
            t = kwargs['t']
        if 'tNum' in kwargs:
            tNum = kwargs['tNum']
        if 'tlb' in kwargs:
            tlb = kwargs['tlb']
        if 'tub' in kwargs:
            tub = kwargs['tub']

        self.cptr = new TACSSolidConstitutive(props.ptr, t, tNum,
                                              tlb, tub)
        self.ptr = self.cptr
        self.ptr.incref()
