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

# Import from TACS for definitions
from TACS cimport *

cdef extern from "TACSFunction.h":
    enum FunctionDomain:
        ENTIRE_DOMAIN
        SUB_DOMAIN
        NO_DOMAIN

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
    enum KSTemperatureType"TACSKTemperature::KSTemperatureType":
        KS_TEMPERATURE_DISCRETE"TACSKSTemperature::DISCRETE"
        KS_TEMPERATURE_CONTINUOUS"TACSKSTemperature::CONTINUOUS"
        PNORM_TEMPERATURE_DISCRETE"TACSKSTemperature::PNORM_DISCRETE"
        PNORM_TEMPERATURE_CONTINUOUS"TACSKSTemperature::PNORM_CONTINUOUS"

    cdef cppclass TACSKSTemperature(TACSFunction):
        TACSKSTemperature(TACSAssembler*, double, double)
        void setKSTemperatureType(KSTemperatureType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxFailOffset(TacsScalar)

cdef extern from "TACSKSFailure.h":
    enum KSFailureType"TACSKFailure::KSFailureType":
        KS_FAILURE_DISCRETE"TACSKSFailure::DISCRETE"
        KS_FAILURE_CONTINUOUS"TACSKSFailure::CONTINUOUS"
        PNORM_FAILURE_DISCRETE"TACSKSFailure::PNORM_DISCRETE"
        PNORM_FAILURE_CONTINUOUS"TACSKSFailure::PNORM_CONTINUOUS"

    cdef cppclass TACSKSFailure(TACSFunction):
        TACSKSFailure(TACSAssembler*, double, double, double)
        void setKSFailureType(KSFailureType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxFailOffset(TacsScalar)

cdef extern from "TACSKSDisplacement.h":
    enum KSDisplacementType"TACSKDisplacement::KSDisplacementType":
        KS_DISPLACEMENT_DISCRETE"TACSKSDisplacement::DISCRETE"
        KS_DISPLACEMENT_CONTINUOUS"TACSKSDisplacement::CONTINUOUS"
        PNORM_DISPLACEMENT_DISCRETE"TACSKSDisplacement::PNORM_DISCRETE"
        PNORM_DISPLACEMENT_CONTINUOUS"TACSKSDisplacement::PNORM_CONTINUOUS"

    cdef cppclass TACSKSDisplacement(TACSFunction):
        TACSKSDisplacement(TACSAssembler*, double, const double*, double)
        void setKSDisplacementType(KSDisplacementType ftype)
        double getParameter()
        void setParameter(double)
        void setMaxDispOffset(TacsScalar)

cdef extern from "TACSInducedFailure.h":
    enum InducedNormType"TACSInducedFailure::InducedNormType":
        INDUCED_EXPONENTIAL"TACSInducedFailure::EXPONENTIAL"
        INDUCED_POWER"TACSInducedFailure::POWER"
        INDUCED_EXPONENTIAL_SQUARED"TACSInducedFailure::EXPONENTIAL_SQUARED"
        INDUCED_POWER_SQUARED"TACSInducedFailure::POWER_SQUARED"
        INDUCED_DISCRETE_EXPONENTIAL"TACSInducedFailure::DISCRETE_EXPONENTIAL"
        INDUCED_DISCRETE_POWER"TACSInducedFailure::DISCRETE_POWER"
        INDUCED_DISCRETE_EXPONENTIAL_SQUARED"TACSInducedFailure::DISCRETE_EXPONENTIAL_SQUARED"
        INDUCED_DISCRETE_POWER_SQUARED"TACSInducedFailure::DISCRETE_POWER_SQUARED"

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
