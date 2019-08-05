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

cdef extern from "RigidBody.h":
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

    cdef cppclass TACSLinearizedMotionDriver(TACSElement):
        TACSLinearizedMotionDriver(TACSGibbsVector*, TacsScalar, int)

    cdef cppclass TACSCylindricalConstraint(TACSElement):
        TACSCylindricalConstraint(TACSRigidBody *bodyA, TACSRigidBody *bodyB,
                                  TACSGibbsVector *point, TACSGibbsVector *eA)
        TACSCylindricalConstraint(TACSRigidBody *bodyA,
                                  TACSGibbsVector *point, TACSGibbsVector *eA)

    cdef cppclass TACSPrismaticConstraint(TACSElement):
        TACSPrismaticConstraint(TACSRigidBody *bodyA, TACSRigidBody *bodyB,
                                TACSGibbsVector *point, TACSGibbsVector *eA)
        TACSPrismaticConstraint(TACSRigidBody *bodyA,
                                TACSGibbsVector *point, TACSGibbsVector *eA)

    cdef cppclass TACSSlidingPivotConstraint(TACSElement):
        TACSSlidingPivotConstraint(TACSRigidBody *bodyA, TACSRigidBody *bodyB,
                                   TACSGibbsVector *point, TACSGibbsVector *eA)
        TACSSlidingPivotConstraint(TACSRigidBody *bodyA,
                                   TACSGibbsVector *point, TACSGibbsVector *eA)

    cdef cppclass TACSAverageConstraint(TACSElement):
        TACSAverageConstraint(TACSRigidBody*, TACSGibbsVector*,
                              TACSRefFrame*, int)

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

    cdef cppclass PlaneStressQuad5(TACSElement):
        PlaneStressQuad5(PlaneStressStiffness *stiff,
                         ElementBehaviorType type, int)

    # Declare theee PSQuadTraction elements
    cdef cppclass PSQuadTraction2(TACSElement):
        PSQuadTraction2(int, TacsScalar, TacsScalar)
        PSQuadTraction2(int, TacsScalar*, TacsScalar*)

    cdef cppclass PSQuadTraction3(TACSElement):
        PSQuadTraction3(int, TacsScalar, TacsScalar)
        PSQuadTraction3(int, TacsScalar*, TacsScalar*)

    cdef cppclass PSQuadTraction4(TACSElement):
        PSQuadTraction4(int, TacsScalar, TacsScalar)
        PSQuadTraction4(int, TacsScalar*, TacsScalar*)

    cdef cppclass PSQuadTraction5(TACSElement):
        PSQuadTraction5(int, TacsScalar, TacsScalar)
        PSQuadTraction5(int, TacsScalar*, TacsScalar*)

    # Declare the Poisson Quadrilateral elements
    ctypedef void (*shell_evalf)(void*, const TacsScalar*, TacsScalar*)

    # Declare shell traction
    cdef cppclass TACSShellTraction2(TACSElement):
        TACSShellTraction2(TacsScalar, TacsScalar, TacsScalar)
        TACSShellTraction2(TacsScalar*, TacsScalar*, TacsScalar*)
        TACSShellTraction2(void*, shell_evalf)

    cdef cppclass TACSShellTraction3(TACSElement):
        TACSShellTraction3(TacsScalar, TacsScalar, TacsScalar)
        TACSShellTraction3(TacsScalar*, TacsScalar*, TacsScalar*)
        TACSShellTraction3(void*, shell_evalf)

    cdef cppclass TACSShellTraction4(TACSElement):
        TACSShellTraction4(TacsScalar, TacsScalar, TacsScalar)
        TACSShellTraction4(TacsScalar*, TacsScalar*, TacsScalar*)
        TACSShellTraction4(void*, shell_evalf)

    cdef cppclass TACSShellTraction5(TACSElement):
        TACSShellTraction5(TacsScalar, TacsScalar, TacsScalar)
        TACSShellTraction5(TacsScalar*, TacsScalar*, TacsScalar*)
        TACSShellTraction5(void*, shell_evalf)

    # Declare the 3D traction classes
    cdef cppclass TACS3DTraction2(TACSElement):
        TACS3DTraction2(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DTraction2(int, TacsScalar[], TacsScalar[], TacsScalar[])

    cdef cppclass TACS3DTraction3(TACSElement):
        TACS3DTraction3(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DTraction3(int, TacsScalar[], TacsScalar[], TacsScalar[])

    cdef cppclass TACS3DTraction4(TACSElement):
        TACS3DTraction4(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DTraction4(int, TacsScalar[], TacsScalar[], TacsScalar[])

    cdef cppclass TACS3DTraction5(TACSElement):
        TACS3DTraction5(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DTraction5(int, TacsScalar[], TacsScalar[], TacsScalar[])

    cdef cppclass TACS3DPressureTraction2(TACSElement):
        TACS3DPressureTraction2(int, TacsScalar)

    cdef cppclass TACS3DPressureTraction3(TACSElement):
        TACS3DPressureTraction3(int, TacsScalar)

    cdef cppclass TACS3DPressureTraction4(TACSElement):
        TACS3DPressureTraction4(int, TacsScalar)

    cdef cppclass TACS3DPressureTraction5(TACSElement):
        TACS3DPressureTraction5(int, TacsScalar)

    cdef cppclass TACS3DBoundingTraction2(TACSElement):
        TACS3DBoundingTraction2(int, TacsScalar, TacsScalar, TacsScalar, TacsScalar[])
        TACS3DBoundingTraction2(int, TacsScalar[], TacsScalar[], TacsScalar[], TacsScalar[])

    # Declare the MITCShell elements
    cdef cppclass MITCShell2(TACSElement):
        MITCShell2(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell32(TACSElement):
        MITCShell32(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell3(TACSElement):
        MITCShell3(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell4(TACSElement):
        MITCShell4(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell43(TACSElement):
        MITCShell43(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell5(TACSElement):
        MITCShell5(FSDTStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass MITCShell54(TACSElement):
        MITCShell54(FSDTStiffness *stiff, ElementBehaviorType type, int)

    # Declare the Solid elements
    cdef cppclass Solid2(TACSElement):
        Solid2(SolidStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass Solid3(TACSElement):
        Solid3(SolidStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass Solid4(TACSElement):
        Solid4(SolidStiffness *stiff, ElementBehaviorType type, int)

    cdef cppclass Solid5(TACSElement):
        Solid5(SolidStiffness *stiff, ElementBehaviorType type, int)

    # Declare the Poisson Quadrilateral elements
    ctypedef TacsScalar (*poisson_evalf)(void*, const TacsScalar*)

    cdef cppclass PoissonQuad2(TACSElement):
        PoissonQuad2(TacsScalar*)
        PoissonQuad2(void*, poisson_evalf)

    cdef cppclass PoissonQuad3(TACSElement):
        PoissonQuad3(TacsScalar*)
        PoissonQuad3(void*, poisson_evalf)

    cdef cppclass PoissonQuad4(TACSElement):
        PoissonQuad4(TacsScalar*)
        PoissonQuad4(void*, poisson_evalf)

    cdef cppclass PoissonQuad5(TACSElement):
        PoissonQuad5(TacsScalar*)
        PoissonQuad5(void*, poisson_evalf)

    # Declare the Plane Stress Thermoelastic Quad elements
    cdef cppclass PSThermoQuad2(TACSElement):
        PSThermoQuad2(CoupledThermoPlaneStressStiffness *stiff,
                      ElementBehaviorType type, int)

    cdef cppclass PSThermoQuad3(TACSElement):
        PSThermoQuad3(CoupledThermoPlaneStressStiffness *stiff,
                      ElementBehaviorType type, int)

    cdef cppclass PSThermoQuad4(TACSElement):
        PSThermoQuad4(CoupledThermoPlaneStressStiffness *stiff,
                      ElementBehaviorType type, int)

    cdef cppclass PSThermoQuad5(TACSElement):
        PSThermoQuad5(CoupledThermoPlaneStressStiffness *stiff,
                      ElementBehaviorType type, int)

    cdef cppclass PSThermoQuad6(TACSElement):
        PSThermoQuad6(CoupledThermoPlaneStressStiffness *stiff,
                      ElementBehaviorType type, int)

    # Declare the Plane Stress Thermoelastic Traction Quad elements
    cdef cppclass PSThermoQuadTraction2(TACSElement):
        PSThermoQuadTraction2(int, TacsScalar, TacsScalar)
        PSThermoQuadTraction2(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadTraction3(TACSElement):
        PSThermoQuadTraction3(int, TacsScalar, TacsScalar)
        PSThermoQuadTraction3(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadTraction4(TACSElement):
        PSThermoQuadTraction4(int, TacsScalar, TacsScalar)
        PSThermoQuadTraction4(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadTraction5(TACSElement):
        PSThermoQuadTraction5(int, TacsScalar, TacsScalar)
        PSThermoQuadTraction5(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadTraction6(TACSElement):
        PSThermoQuadTraction6(int, TacsScalar, TacsScalar)
        PSThermoQuadTraction6(int, TacsScalar*, TacsScalar*)

    # Declare Plane Stress Thermoelastic Heat Flux elements
    cdef cppclass PSThermoQuadHF2(TACSElement):
        PSThermoQuadHF2(int, TacsScalar, TacsScalar)
        PSThermoQuadHF2(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadHF3(TACSElement):
        PSThermoQuadHF3(int, TacsScalar, TacsScalar)
        PSThermoQuadHF3(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadHF4(TACSElement):
        PSThermoQuadHF4(int, TacsScalar, TacsScalar)
        PSThermoQuadHF4(int, TacsScalar*, TacsScalar*)
        
    cdef cppclass PSThermoQuadHF5(TACSElement):
        PSThermoQuadHF5(int, TacsScalar, TacsScalar)
        PSThermoQuadHF5(int, TacsScalar*, TacsScalar*)

    # Declare the Solid Thermoelastic elements
    cdef cppclass SolidThermo2(TACSElement):
        SolidThermo2(CoupledThermoSolidStiffness *stiff,
                     ElementBehaviorType type, int)
        
    cdef cppclass SolidThermo3(TACSElement):
        SolidThermo3(CoupledThermoSolidStiffness *stiff,
                     ElementBehaviorType type, int)
        
    cdef cppclass SolidThermo4(TACSElement):
        SolidThermo4(CoupledThermoSolidStiffness *stiff,
                     ElementBehaviorType type, int)
        
    cdef cppclass SolidThermo5(TACSElement):
        SolidThermo5(CoupledThermoSolidStiffness *stiff,
                     ElementBehaviorType type, int)
        
    cdef cppclass SolidThermo6(TACSElement):
        SolidThermo6(CoupledThermoSolidStiffness *stiff,
                     ElementBehaviorType type, int)

    # Declare the Solid Thermoelastic Traction elements
    cdef cppclass TACS3DThermoTraction2(TACSElement):
        TACS3DThermoTraction2(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoTraction2(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoTraction3(TACSElement):
        TACS3DThermoTraction3(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoTraction3(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoTraction4(TACSElement):
        TACS3DThermoTraction4(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoTraction4(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoTraction5(TACSElement):
        TACS3DThermoTraction5(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoTraction5(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoTraction6(TACSElement):
        TACS3DThermoTraction6(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoTraction6(int, TacsScalar*, TacsScalar*, TacsScalar*)

    # Declare the Solid Thermoelastic Pressure Traction elements
    cdef cppclass TACS3DThermoPressureTraction2(TACSElement):
        TACS3DThermoPressureTraction2(int, TacsScalar)

    cdef cppclass TACS3DThermoPressureTraction3(TACSElement):
        TACS3DThermoPressureTraction3(int, TacsScalar)

    cdef cppclass TACS3DThermoPressureTraction4(TACSElement):
        TACS3DThermoPressureTraction4(int, TacsScalar)

    cdef cppclass TACS3DThermoPressureTraction5(TACSElement):
        TACS3DThermoPressureTraction5(int, TacsScalar)

    cdef cppclass TACS3DThermoPressureTraction6(TACSElement):
        TACS3DThermoPressureTraction6(int, TacsScalar)

    # Declare the Solid Thermoelastic Heat Flux elements
    cdef cppclass TACS3DThermoHF2(TACSElement):
        TACS3DThermoHF2(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoHF2(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoHF3(TACSElement):
        TACS3DThermoHF3(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoHF3(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoHF4(TACSElement):
        TACS3DThermoHF4(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoHF4(int, TacsScalar*, TacsScalar*, TacsScalar*)
        
    cdef cppclass TACS3DThermoHF5(TACSElement):
        TACS3DThermoHF5(int, TacsScalar, TacsScalar, TacsScalar)
        TACS3DThermoHF5(int, TacsScalar*, TacsScalar*, TacsScalar*)

    # Declare the Solid Thermoelastic Normal Heat Flux elements
    cdef cppclass TACS3DThermoNormalHF2(TACSElement):
        TACS3DThermoNormalHF2(int, TacsScalar)

    cdef cppclass TACS3DThermoNormalHF3(TACSElement):
        TACS3DThermoNormalHF3(int, TacsScalar)

    cdef cppclass TACS3DThermoNormalHF4(TACSElement):
        TACS3DThermoNormalHF4(int, TacsScalar)

    cdef cppclass TACS3DThermoNormalHF5(TACSElement):
        TACS3DThermoNormalHF5(int, TacsScalar)

cdef extern from "PlaneStressTri6.h":
    cdef cppclass PlaneStressTri6(TACSElement):
        PlaneStressTri6(PlaneStressStiffness *stiff,
                        ElementBehaviorType type, int)

cdef extern from  "MITC9.h":
    cdef cppclass MITC9(TACSElement):
        MITC9(FSDTStiffness *_stiff, TACSGibbsVector *_gravity,
              TACSGibbsVector *_vInit, TACSGibbsVector *_omegaInit)

cdef extern from  "MITC3.h":
    cdef cppclass MITC3(TACSElement):
        MITC3(TimoshenkoStiffness *_stiff, TACSGibbsVector *_gravity,
              TACSGibbsVector *_vInit, TACSGibbsVector *_omegaInit)

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
