# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

# Import from TACS for definitions
from TACS cimport *

cdef extern from "FSDTStiffness.h":
    cdef cppclass FSDTStiffness(TACSConstitutive):
        FSDTStiffness()

cdef extern from "isoFSDTStiffness.h":
    cdef cppclass isoFSDTStiffness(FSDTStiffness):
        isoFSDTStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu, 
                         TacsScalar kcorr, TacsScalar yieldStress, 
                         TacsScalar thickness, int tNum,
                         TacsScalar minThickness, 
                         TacsScalar maxThickness)

cdef extern from "PlaneStressStiffness.h":
    cdef cppclass PlaneStressStiffness(TACSConstitutive):
        PlaneStressStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu)
        
cdef extern from "SolidStiffness.h":
    cdef cppclass SolidStiffness(TACSConstitutive):
        SolidStiffness(TacsScalar rho, TacsScalar E, TacsScalar nu)

cdef class FSDT(Constitutive):
    pass

cdef class PlaneStress(Constitutive):
    pass
   
cdef class Solid(Constitutive):
    pass

# Special functions required for converting pointers
cdef extern from "":
    PlaneStressStiffness* _dynamicPlaneStress"dynamic_cast<PlaneStressStiffness*>"(TACSConstitutive*) except NULL
    FSDTStiffness* _dynamicFSDT"dynamic_cast<FSDTStiffness*>"(TACSConstitutive*) except NULL
