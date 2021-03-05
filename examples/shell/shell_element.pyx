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

# Import the TACS module
from tacs.TACS cimport *
from tacs.constitutive cimport *
from tacs.elements cimport *

cdef extern from "mpi-compat.h":
    pass

cdef extern from "TACSShellElementTransform.h":
    cdef cppclass TACSShellTransform(TACSObject):
        TACSShellTransform()

    cdef cppclass TACSShellNaturalTransform(TACSShellTransform):
        TACSShellNaturalTransform()

    cdef cppclass TACSShellRefAxisTransform(TACSShellTransform):
        TACSShellRefAxisTransform(TacsScalar*)

cdef extern from "TACSShellElement.h":
    cdef cppclass TACSQuadLinearShell"TACSShellElement<TACSQuadLinearQuadrature, TACSShellQuadLinearBasis, TACSLinearizedRotation, TACSShellLinearModel>"(TACSElement):
        TACSQuadLinearShell(TACSShellTransform*,
                            TACSShellConstitutive*)

    cdef cppclass TACSQuadQuadraticShell"TACSShellElement<TACSQuadQuadraticQuadrature, TACSShellQuadQuadraticBasis, TACSLinearizedRotation, TACSShellLinearModel>"(TACSElement):
        TACSQuadQuadraticShell(TACSShellTransform*,
                               TACSShellConstitutive*)

cdef class ShellTransform:
    cdef TACSShellTransform *ptr
    def __cinit__(self):
        self.ptr = NULL
        return

    def __dealloc__(self):
        if self.ptr != NULL:
            self.ptr.decref()

cdef class ShellNaturalTransform(ShellTransform):
    def __cinit__(self):
        self.ptr = new TACSShellNaturalTransform()
        self.ptr.incref()

cdef class ShellRefAxisTransform(ShellTransform):
    def __cinit__(self, axis):
        cdef TacsScalar a[3]
        a[0] = axis[0]
        a[1] = axis[1]
        a[2] = axis[2]
        self.ptr = new TACSShellRefAxisTransform(a)
        self.ptr.incref()

cdef class QuadLinearShell(Element):
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        self.ptr = new TACSQuadLinearShell(transform.ptr, con.cptr)
        self.ptr.incref()

cdef class QuadQuadraticShell(Element):
    def __cinit__(self, ShellTransform transform, ShellConstitutive con):
        self.ptr = new TACSQuadQuadraticShell(transform.ptr, con.cptr)
        self.ptr.incref()
