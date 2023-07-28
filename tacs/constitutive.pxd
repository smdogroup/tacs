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

cdef extern from "TACSIsoShellConstitutive.h":
    cdef cppclass TACSIsoShellConstitutive(TACSShellConstitutive):
        TACSIsoShellConstitutive(TACSMaterialProperties*, TacsScalar, int,
                                 TacsScalar, TacsScalar)

cdef extern from "TACSCompositeShellConstitutive.h":
    cdef cppclass TACSCompositeShellConstitutive(TACSShellConstitutive):
        TACSCompositeShellConstitutive(int, TACSOrthotropicPly**, const TacsScalar*,
                                       const TacsScalar*, TacsScalar)
        void getPlyThicknesses(TacsScalar*);
        void getPlyAngles(TacsScalar*);

cdef extern from "TACSSmearedCompositeShellConstitutive.h":
    cdef cppclass TACSSmearedCompositeShellConstitutive(TACSShellConstitutive):
        TACSSmearedCompositeShellConstitutive(int, TACSOrthotropicPly**, TacsScalar,
                                       const TacsScalar*, const TacsScalar*,
                                       int, const int*,
                                       TacsScalar, TacsScalar,
                                       const TacsScalar*, const TacsScalar*)
        TacsScalar getLaminateThickness();
        void getPlyAngles(TacsScalar*);
        void getPlyFractions(TacsScalar *);

cdef extern from "TACSLamParamShellConstitutive.h":
    cdef cppclass TACSLamParamShellConstitutive(TACSShellConstitutive):
        TACSLamParamShellConstitutive(TACSOrthotropicPly*, TacsScalar, int, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, TacsScalar, int, int, int,
                                      TacsScalar, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, int, int, TacsScalar, TacsScalar)

cdef extern from "TACSBladeStiffenedShellConstitutive.h":
    cdef cppclass TACSBladeStiffenedShellConstitutive(TACSShellConstitutive):
        TACSBladeStiffenedShellConstitutive(
            TACSOrthotropicPly*, # panelPly
            TACSOrthotropicPly*, # stiffenerPly
            TacsScalar, # kcorr
            TacsScalar, # panelLength
            int, # panelLengthNum
            TacsScalar, # stiffenerPitch
            int, # stiffenerPitchNum
            TacsScalar, # panelThick
            int, # panelThickNum
            int, # numPanelPlies
            TacsScalar[], # panelPlyAngles
            TacsScalar[], # panelPlyFracs
            int[], # panelPlyFracNums
            TacsScalar, # stiffenerHeight
            int, # stiffenerHeightNum
            TacsScalar, # stiffenerThick
            int, # stiffenerThickNum
            int, # numStiffenerPlies
            TacsScalar[], # stiffenerPlyAngles
            TacsScalar[], # stiffenerPlyFracs
            int[], # stiffenerPlyFracNums
            TacsScalar # flangeFraction
        )
        int getNumPanelPlies()
        int getNumStiffenerPlies()
        void setKSWeight(double ksWeight)
        void setStiffenerPitchBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerHeightBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerThicknessBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setPanelThicknessBounds(TacsScalar lowerBound, TacsScalar upperBound)
        void setStiffenerPlyFractionBounds(TacsScalar[] lowerBound, TacsScalar[] upperBound)
        void setPanelPlyFractionBounds(TacsScalar[] lowerBound, TacsScalar[] upperBound)

cdef class BladeStiffenedShellConstitutive(ShellConstitutive):
    cdef TACSBladeStiffenedShellConstitutive *blade_ptr

cdef class BeamConstitutive(Constitutive):
    cdef TACSBeamConstitutive *cptr

cdef class GeneralMassConstitutive(Constitutive):
    cdef TACSGeneralMassConstitutive *cptr

cdef class GeneralSpringConstitutive(Constitutive):
    cdef TACSGeneralSpringConstitutive *cptr
