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

# Import from TACS for definitions
from TACS cimport *

cdef extern from "TACSFunction.h":
    enum FunctionDomain:
        ENTIRE_DOMAIN
        SUB_DOMAIN
        NO_DOMAIN

cdef extern from "TACSStructuralMass.h":
    cdef cppclass TACSStructuralMass(TACSFunction):
        TACSStructuralMass(TACSAssembler *tacs)

cdef extern from "TACSKSFailure.h":
    enum KSFailureType"TACSKFailure::KSFailureType":
        KS_FAILURE_DISCRETE"TACSKSFailure::DISCRETE"
        KS_FAILURE_CONTINUOUS"TACSKSFailure::CONTINUOUS"
        PNORM_FAILURE_DISCRETE"TACSKSFailure::PNORM_DISCRETE"
        PNORM_FAILURE_CONTINUOUS"TACSKSFailure::PNORM_CONTINUOUS"

    enum KSConstitutiveFunction"TACSKSFailure::KSConstitutiveFunction":
        KS_FAILURE"TACSKSFailure::FAILURE"
        KS_BUCKLING"TACSKSFailure::BUCKLING"

    cdef cppclass TACSKSFailure(TACSFunction):
        TACSKSFailure(TACSAssembler *tacs, double ksWeight,
                      KSConstitutiveFunction func,
                      double alpha)
        void setKSFailureType(KSFailureType ftype)
        double getParameter()
        void setParameter(double _ksWeight)
        void setLoadFactor(TacsScalar _loadFactor)
        void setMaxFailOffset(TacsScalar _maxFail)

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

    enum InducedConstitutiveFunction"TACSInducedFailure::InducedConstitutiveFunction":
        INDUCED_FAILURE"TACSInducedFailure::FAILURE"
        INDUCED_BUCKLING"TACSInducedFailure::BUCKLING"

    cdef cppclass TACSInducedFailure(TACSFunction):
        TACSInducedFailure(TACSAssembler *tacs, double P,
                           InducedConstitutiveFunction func)
        void setInducedType(InducedNormType _norm_type)
        void setParameter(double _P)
        double getParameter()
        void setMaxFailOffset(TacsScalar _max_fail)

cdef extern from "TACSHeatFlux.h":
    cdef cppclass HeatFluxIntegral( TACSFunction ):
        HeatFluxIntegral(TACSAssembler*, int*, int*, int)

# cdef extern from "TACSDisplacementIntegral.h":
#     cdef cppclass TACSDisplacementIntegral(TACSFunction):
#         TACSDisplacementIntegral(TACSAssembler *tacs, TacsScalar dir[])
