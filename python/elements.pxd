# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

# Import from constitutive for definitions
from constitutive cimport *
from TACS cimport *

cdef extern from "TACSElement.h":
    enum ElementBehaviorType:
        LINEAR
        NONLINEAR
        LARGE_ROTATION

cdef extern from "TACSGibbsVector.h":
    cdef cppclass TACSGibbsVector(TACSObject):
        TACSGibbsVector(TacsScalar x[])

# Template
cdef extern from "TACSElementTemplates.h":
    # Declare the PlaneStressQuad elements
    cdef cppclass PlaneStressQuad2(TACSElement):
        PlaneStressQuad2(PlaneStressStiffness *stiff, 
                         ElementBehaviorType type, int)
        
    cdef cppclass PlaneStressQuad3(TACSElement):
        PlaneStressQuad3(PlaneStressStiffness *stiff, 
                         ElementBehaviorType type, int)
        
    cdef cppclass PlaneStressQuad4(TACSElement):
        PlaneStressQuad4(PlaneStressStiffness *stiff, 
                         ElementBehaviorType type, int)

    # Declare the MITCShell elements
    cdef cppclass MITCShell2(TACSElement):
        MITCShell2(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell3(TACSElement):
        MITCShell3(FSDTStiffness *stiff, ElementBehaviorType type, int)
        
    cdef cppclass MITCShell4(TACSElement):
        MITCShell4(FSDTStiffness *stiff, ElementBehaviorType type, int)

    # Declare the Solid elements
    cdef cppclass Solid2(TACSElement):
        Solid2(SolidStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass Solid3(TACSElement):
        Solid3(SolidStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass Solid4(TACSElement):
        Solid4(SolidStiffness *stiff, ElementBehaviorType type, int)


cdef extern from "PlaneStressTri6.h":
    cdef cppclass PlaneStressTri6(TACSElement):
        PlaneStressTri6(PlaneStressStiffness *stiff, 
                        ElementBehaviorType type, int)

cdef extern from  "MITC9.h":
    cdef cppclass MITC9(TACSElement):
        MITC9(FSDTStiffness *_stiff, TACSGibbsVector *_gravity=NULL,
              TACSGibbsVector *_vInit=NULL, TACSGibbsVector *_omegaInit=NULL)
