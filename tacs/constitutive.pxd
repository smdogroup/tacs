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

#distutils: language=c++

# Import from cpp headers for definitions
from tacs.cpp_headers.TACS cimport *
from tacs.cpp_headers.constitutive cimport *

# Import cython headers
from tacs.TACS cimport *

cdef class MaterialProperties:
    cdef TACSMaterialProperties *ptr
    cdef int nastranID

cdef inline _init_MaterialProperties(TACSMaterialProperties *ptr):
    props = MaterialProperties()
    props.ptr = ptr
    props.ptr.incref()
    return props

cdef class PlaneStressConstitutive(Constitutive):
    cdef TACSPlaneStressConstitutive *cptr

cdef class PhaseChangeMaterialConstitutive(Constitutive):
    cdef TACSPhaseChangeMaterialConstitutive *cptr

cdef class SolidConstitutive(Constitutive):
    cdef TACSSolidConstitutive *cptr

cdef class ShellConstitutive(Constitutive):
    cdef TACSShellConstitutive *cptr

cdef class BladeStiffenedShellConstitutive(ShellConstitutive):
    cdef TACSBladeStiffenedShellConstitutive *blade_ptr

cdef class BeamConstitutive(Constitutive):
    cdef TACSBeamConstitutive *cptr

cdef class GeneralMassConstitutive(Constitutive):
    cdef TACSGeneralMassConstitutive *cptr

cdef class GeneralSpringConstitutive(Constitutive):
    cdef TACSGeneralSpringConstitutive *cptr
