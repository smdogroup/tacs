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

# cdef extern from "TACSGibbsVector.h":
#     cdef cppcalss TACSGibbsVector(TACSObject):
#         TACSGibbsVector(TacsScalar _x[]);

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

cdef extern from "PlaneStressTri6.h":
    cdef cppclass PlaneStressTri6(TACSElement):
        PlaneStressTri6(PlaneStressStiffness *stiff, 
                        ElementBehaviorType type, int)

# cdef extern from "MITCShell.h":
#     cdef cppclass MITCShell(TACSShell)[int]:
#         MITCShell(FSDTStiffness *_stiff, ElementBehaviorType _type, int)

# cdef extern from  "MITC9.h":
#     cdef cppclass MITC9(TACSElement):
#         MITC9(FSDTStiffness *_stiff, TACSGibbsVector *_gravity=NULL,
#               TACSGibbsVector *_vInit=NULL, TACSGibbsVector *_omegaInit=NULL)

# cdef extern from "Solid.h":
#     enum SolidElementType"Solid::SolidElementType":
#         LINEAR"Solid::LINEAR"
#         NONLINEAR"Solid::NONLINEAR"
#     cdef cppclass Solid(TACS3DElement)[int]:
#         Solid(SolidStiffness *_stiff, SolidElementType elem_type=LINEAR, int)
