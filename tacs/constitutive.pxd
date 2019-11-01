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

# For MPI capabilities
from mpi4py.libmpi cimport *
cimport mpi4py.MPI as MPI

# Import numpy
import numpy as np
cimport numpy as np
from libc.string cimport const_char

# Import C methods for python
from cpython cimport PyObject, Py_INCREF

# Import from TACS for definitions
from TACS cimport *

cdef extern from "TACSMaterialProperties.h":
    cdef cppclass TACSMaterialProperties(TACSObject):
        TACSMaterialProperties(TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar)

cdef class MaterialProperties:
    cdef TACSMaterialProperties *ptr

cdef inline _init_MaterialProperties(TACSMaterialProperties *ptr):
    props = MaterialProperties()
    props.ptr = ptr
    props.ptr.incref()
    return props

cdef extern from "TACSPlaneStressConstitutive.h":
    cdef cppclass TACSPlaneStressConstitutive(TACSConstitutive):
        TACSPlaneStressConstitutive(TACSMaterialProperties*,
                                    TacsScalar, int, TacsScalar, TacsScalar)

cdef class PlaneStressConstitutive(Constitutive):
    cdef TACSPlaneStressConstitutive *cptr

cdef extern from "TACSSolidConstitutive.h":
    cdef cppclass TACSSolidConstitutive(TACSConstitutive):
        TACSSolidConstitutive(TACSMaterialProperties*,
                              TacsScalar, int, TacsScalar, TacsScalar)

cdef class SolidConstitutive(Constitutive):
    cdef TACSSolidConstitutive *cptr

cdef extern from "TACSConstitutiveVerification.h":
    int TacsTestConstitutive(TACSConstitutive*, int, double, int, double, double)
