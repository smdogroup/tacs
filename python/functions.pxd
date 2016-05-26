# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

# Import from TACS for definitions
from TACS cimport *

cdef extern from "TACSFunction.h":
    enum FunctionDomain:
        ENTIRE_DOMAIN
        SUB_DOMAIN
        NO_DOMAIN        
    
cdef extern from "Compliance.h":
    cdef cppclass Compliance(TACSFunction):
        Compliance(TACSAssembler *tacs)

cdef extern from "StructuralMass.h":
    cdef cppclass StructuralMass(TACSFunction):
        StructuralMass(TACSAssembler *tacs)

cdef extern from "KSFailure.h":
    cdef cppclass KSFailure(TACSFunction):
        KSFailure(TACSAssembler *tacs, double ksWeight, double alpha)

cdef extern from "KSBuckling.h":
    cdef cppclass KSBuckling(TACSFunction):
        KSBuckling(TACSAssembler* tacs, TacsScalar ksWeight)

cdef extern from "KSDisplacement.h":
    cdef cppclass KSDisplacement(TACSFunction):
        KSDisplacement(TACSAssembler* tacs, TacsScalar d[],
                       double ksWeight, double alpha)

cdef extern from "InducedBuckling.h":
    cdef cppclass InducedBuckling(TACSFunction):
        InducedBuckling(TACSAssembler* tacs, int elementNums[],
                        int numElements, double P)

cdef extern from "InducedFailure.h":
    cdef cppclass InducedFailure(TACSFunction):
        InducedFailure(TACSAssembler* tacs, int elementNums[],
                       int numElements, double P)
    
