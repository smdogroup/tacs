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

# Import from TACS for definitions
from tacs.cpp_headers.TACS cimport *

cdef extern from "TACSFunction.h":
    enum FunctionDomain:
        ENTIRE_DOMAIN
        SUB_DOMAIN
        NO_DOMAIN

    # Declare the C++ enum under a private Cython name to avoid conflicting with
    # the Python IntEnum class also named KSAggregationType.
    enum _CKSAggregationType "KSAggregationType":
        _CKSAGG_DISCRETE "KS_DISCRETE"
        _CKSAGG_CONTINUOUS "KS_CONTINUOUS"
        _CKSAGG_PNORM_DISCRETE "PNORM_DISCRETE"
        _CKSAGG_PNORM_CONTINUOUS "PNORM_CONTINUOUS"
        _CKSAGG_DISCRETE_AVERAGE "KS_DISCRETE_AVERAGE"

    # Integer aliases for populating the Python IntEnum
    int _KSAGG_DISCRETE "KS_DISCRETE"
    int _KSAGG_CONTINUOUS "KS_CONTINUOUS"
    int _KSAGG_PNORM_DISCRETE "PNORM_DISCRETE"
    int _KSAGG_PNORM_CONTINUOUS "PNORM_CONTINUOUS"
    int _KSAGG_DISCRETE_AVERAGE "KS_DISCRETE_AVERAGE"

cdef extern from "TACSStructuralMass.h":
    cdef cppclass TACSStructuralMass(TACSFunction):
        TACSStructuralMass(TACSAssembler*)

cdef extern from "TACSEnclosedVolume.h":
    cdef cppclass TACSEnclosedVolume(TACSFunction):
        TACSEnclosedVolume(TACSAssembler*)

cdef extern from "TACSCenterOfMass.h":
    cdef cppclass TACSCenterOfMass(TACSFunction):
        TACSCenterOfMass(TACSAssembler*, const double*)

cdef extern from "TACSMomentOfInertia.h":
    cdef cppclass TACSMomentOfInertia(TACSFunction):
        TACSMomentOfInertia(TACSAssembler*, const double*,  const double*, int)

cdef extern from "TACSCompliance.h":
    cdef cppclass TACSCompliance(TACSFunction):
        TACSCompliance(TACSAssembler*)
        void setComplianceType(int)

cdef extern from "TACSAverageTemperature.h":
    cdef cppclass TACSAverageTemperature(TACSFunction):
        TACSAverageTemperature(TACSAssembler*, TacsScalar)

cdef extern from "TACSKSTemperature.h":
    cdef cppclass TACSKSTemperature(TACSFunction):
        TACSKSTemperature(TACSAssembler*, double, double)
        void setKSAggregationType(_CKSAggregationType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxFailOffset(TacsScalar)

cdef extern from "TACSKSFailure.h":
    cdef cppclass TACSKSFailure(TACSFunction):
        TACSKSFailure(TACSAssembler*, double, double, double)
        void setKSAggregationType(_CKSAggregationType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxFailOffset(TacsScalar)

cdef extern from "TACSKSDisplacement.h":
    cdef cppclass TACSKSDisplacement(TACSFunction):
        TACSKSDisplacement(TACSAssembler*, double, const double*, double)
        void setKSAggregationType(_CKSAggregationType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxDispOffset(TacsScalar)

cdef extern from "TACSInducedFailure.h":
    enum InducedNormType"InducedNormType":
        INDUCED_EXPONENTIAL"EXPONENTIAL"
        INDUCED_POWER"POWER"
        INDUCED_EXPONENTIAL_SQUARED"EXPONENTIAL_SQUARED"
        INDUCED_POWER_SQUARED"POWER_SQUARED"
        INDUCED_DISCRETE_EXPONENTIAL"DISCRETE_EXPONENTIAL"
        INDUCED_DISCRETE_POWER"DISCRETE_POWER"
        INDUCED_DISCRETE_EXPONENTIAL_SQUARED"DISCRETE_EXPONENTIAL_SQUARED"
        INDUCED_DISCRETE_POWER_SQUARED"DISCRETE_POWER_SQUARED"

    cdef cppclass TACSInducedFailure(TACSFunction):
        TACSInducedFailure(TACSAssembler*, double)
        void setInducedType(InducedNormType)
        void setParameter(double)
        double getParameter()
        void setMaxFailOffset(TacsScalar _max_fail)

cdef extern from "TACSHeatFlux.h":
    cdef cppclass TACSHeatFlux(TACSFunction):
        TACSHeatFlux(TACSAssembler*, int*, int*, int)

# cdef extern from "TACSDisplacementIntegral.h":
#     cdef cppclass TACSDisplacementIntegral(TACSFunction):
#         TACSDisplacementIntegral(TACSAssembler *tacs, TacsScalar dir[])
