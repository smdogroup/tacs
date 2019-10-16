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
        TACSGibbsVector(TacsScalar, TacsScalar, TacsScalar)

cdef extern from "TACSRigidBody.h":
    cdef cppclass TACSRefFrame(TACSObject):
        TACSRefFrame(TACSGibbsVector*, TACSGibbsVector*, TACSGibbsVector*)

    cdef cppclass TACSRigidBodyViz(TACSObject):
        TACSRigidBodyViz(int, int, TacsScalar*, int*, TACSGibbsVector*)
        TACSRigidBodyViz(TacsScalar, TacsScalar, TacsScalar)

    cdef cppclass TACSRigidBody(TACSElement):
        TACSRigidBody(TACSRefFrame*, const TacsScalar, const TacsScalar*,
                      const TacsScalar*, TACSGibbsVector*,
                      TACSGibbsVector*, TACSGibbsVector*, TACSGibbsVector*)
        void setDesignVarNums(int, const int*, const int*)
        void setVisualization(TACSRigidBodyViz*)

cdef extern from "KinematicConstraints.h":
    cdef cppclass TACSFixedConstraint(TACSElement):
        TACSFixedConstraint(TACSRigidBody *bodyA,
                            TACSGibbsVector *point)

    cdef cppclass TACSSphericalConstraint(TACSElement):
        TACSSphericalConstraint(TACSRigidBody *bodyA, TACSRigidBody *bodyB,
                                TACSGibbsVector *point)
        TACSSphericalConstraint(TACSRigidBody *bodyA,
                                TACSGibbsVector *point)

    cdef cppclass TACSRevoluteConstraint(TACSElement):
        TACSRevoluteConstraint(TACSRigidBody *bodyA, TACSRigidBody *bodyB,
                               TACSGibbsVector *point, TACSGibbsVector *eA)
        TACSRevoluteConstraint(TACSRigidBody *bodyA,
                               TACSGibbsVector *point, TACSGibbsVector *eA)

    cdef cppclass TACSRigidLink(TACSElement):
        TACSRigidLink(TACSRigidBody*)

    cdef cppclass TACSRevoluteDriver(TACSElement):
        TACSRevoluteDriver(TACSGibbsVector*, TacsScalar)

    cdef cppclass TACSMotionDriver(TACSElement):
        TACSMotionDriver(TACSGibbsVector*, TacsScalar, int)

    cdef cppclass TACSAverageConstraint(TACSElement):
        TACSAverageConstraint(TACSRigidBody*, TACSGibbsVector*,
                              TACSRefFrame*, int)

cdef extern from "TACSElementWrapper.h":
    cdef cppclass TACSElementWrapper(TACSElement):
        TACSElementWrapper(PyObject *obj, int, int)

        # Member functions
        void (*getinitconditions)( void *, int, int,
                                   TacsScalar *, TacsScalar *,
                                   TacsScalar *, const TacsScalar * )

        void (*addresidual)( void *, int, int, double time, TacsScalar *,
                             const TacsScalar *,
                             const TacsScalar *,
                             const TacsScalar *,
                             const TacsScalar * )

        void (*addjacobian)( void *, int, int, double time, TacsScalar *,
                             double alpha, double beta, double gamma,
                             const TacsScalar *,
                             const TacsScalar *,
                             const TacsScalar *,
                             const TacsScalar * )

        void addadjresproduct(void * , int, int, double,
                           TacsScalar, TacsScalar *, int,
                           const TacsScalar *,
                           const TacsScalar *,
                           const TacsScalar *,
                           const TacsScalar *,
                           const TacsScalar * )
