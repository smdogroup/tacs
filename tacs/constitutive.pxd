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
cimport TACS
from TACS cimport *

cdef extern from "TACSMaterialProperties.h":
    cdef cppclass TACSMaterialProperties(TACSObject):
        TACSMaterialProperties(TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar)
        TACSMaterialProperties(TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar,
                               TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar,
                               TacsScalar, TacsScalar, TacsScalar)
        void setDensity(TacsScalar)
        void setSpecificHeat(TacsScalar)
        void getOrthotropicProperties(TacsScalar*, TacsScalar*, TacsScalar*,
                                      TacsScalar*, TacsScalar*, TacsScalar*,
                                      TacsScalar*, TacsScalar*, TacsScalar*)
        void getStrengthProperties(TacsScalar*, TacsScalar*, TacsScalar*,
                                   TacsScalar*, TacsScalar*, TacsScalar*,
                                   TacsScalar*, TacsScalar*, TacsScalar*)
        void getCoefThermalExpansion(TacsScalar*, TacsScalar*, TacsScalar*)
        void getThermalConductivity(TacsScalar*, TacsScalar*, TacsScalar*)

    cdef cppclass TACSOrthotropicPly(TACSObject):
        TACSOrthotropicPly(TacsScalar, TACSMaterialProperties*)
        void setKSWeight(TacsScalar)
        void setUseMaxStrainCriterion()
        void setUseTsaiWuCriterion()

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

cdef extern from "TACSPhaseChangeMaterialConstitutive.h":
    cdef cppclass TACSPhaseChangeMaterialConstitutive(TACSConstitutive):
        TACSPhaseChangeMaterialConstitutive(TACSMaterialProperties*,
                                            TACSMaterialProperties*,
                                            TacsScalar, TacsScalar,
                                            TacsScalar, TacsScalar,
                                            int, TacsScalar, TacsScalar)

cdef class PhaseChangeMaterialConstitutive(Constitutive):
    cdef TACSPhaseChangeMaterialConstitutive *cptr

cdef extern from "TACSSolidConstitutive.h":
    cdef cppclass TACSSolidConstitutive(TACSConstitutive):
        TACSSolidConstitutive(TACSMaterialProperties*,
                              TacsScalar, int, TacsScalar, TacsScalar)
        TACSMaterialProperties* getMaterialProperties()

cdef class SolidConstitutive(Constitutive):
    cdef TACSSolidConstitutive *cptr

cdef extern from "TACSShellConstitutive.h":
    cdef cppclass TACSShellConstitutive(TACSConstitutive):
        void setDrillingRegularization(double)

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

cdef extern from "TACSLamParamShellConstitutive.h":
    cdef cppclass TACSLamParamShellConstitutive(TACSShellConstitutive):
        TACSLamParamShellConstitutive(TACSOrthotropicPly*, TacsScalar, int, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, TacsScalar, int, int, int,
                                      TacsScalar, TacsScalar, TacsScalar,
                                      TacsScalar, TacsScalar, int, int, TacsScalar, TacsScalar)

cdef extern from "TACSBeamConstitutive.h":
    cdef cppclass TACSBeamConstitutive(TACSConstitutive):
        pass

cdef class BeamConstitutive(Constitutive):
    cdef TACSBeamConstitutive *cptr

cdef extern from "TACSBasicBeamConstitutive.h":
    cdef cppclass TACSBasicBeamConstitutive(TACSBeamConstitutive):
        TACSBasicBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar, TacsScalar,
                                  TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSIsoTubeBeamConstitutive.h":
    cdef cppclass TACSIsoTubeBeamConstitutive(TACSBeamConstitutive):
        TACSIsoTubeBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar,
                                    int, int, TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSIsoRectangleBeamConstitutive.h":
    cdef cppclass TACSIsoRectangleBeamConstitutive(TACSBeamConstitutive):
        TACSIsoRectangleBeamConstitutive(TACSMaterialProperties*, TacsScalar, TacsScalar,
                                         int, int, TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSGeneralMassConstitutive.h":
    cdef cppclass TACSGeneralMassConstitutive(TACSConstitutive):
        TACSGeneralMassConstitutive(const TacsScalar*)

cdef class GeneralMassConstitutive(Constitutive):
    cdef TACSGeneralMassConstitutive *cptr

cdef extern from "TACSPointMassConstitutive.h":
    cdef cppclass TACSPointMassConstitutive(TACSGeneralMassConstitutive):
        TACSPointMassConstitutive(TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSGeneralSpringConstitutive.h":
    cdef cppclass TACSGeneralSpringConstitutive(TACSConstitutive):
        TACSGeneralSpringConstitutive(TacsScalar*)

cdef class GeneralSpringConstitutive(Constitutive):
    cdef TACSGeneralSpringConstitutive *cptr

cdef extern from "TACSDOFSpringConstitutive.h":
    cdef cppclass TACSDOFSpringConstitutive(TACSGeneralSpringConstitutive):
        TACSDOFSpringConstitutive(TacsScalar*)

# Special functions required for converting pointers
cdef extern from "":
    TACSBeamConstitutive* _dynamicBeamConstitutive"dynamic_cast<TACSBeamConstitutive*>"(TACSConstitutive*)

cdef extern from "TACSConstitutiveVerification.h":
    int TacsTestConstitutiveDensity(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                    const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveSpecificHeat(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                         const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveHeatFlux(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                     const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveStress(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                   const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveThermalStrain(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                          const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveFailure(TACSConstitutive*, int, const double*, const TacsScalar*, int,
                                    const TacsScalar*, double, int, double, double)
    int TacsTestConstitutiveFailureStrainSens(TACSConstitutive*, int, const double*, const TacsScalar*,
                                              double, int, double, double)
    int TacsTestConstitutive(TACSConstitutive*, int, double, int, double, double)
