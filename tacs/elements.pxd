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

# Import numpy
from libc.string cimport const_char

# Import the major python version
from cpython.version cimport PY_MAJOR_VERSION

# Import numpy
cimport numpy as np
import numpy as np

# Import from constitutive for definitions
from constitutive cimport *
from TACS cimport *

cdef extern from "TACSElementTypes.h":
    const int TACS_ELEMENT_DENSITY
    const int TACS_STRAIN_ENERGY_DENSITY
    const int TACS_FAILURE_INDEX
    const int TACS_HEAT_FLUX
    const int TACS_TEMPERATURE
    const int TACS_TOTAL_STRAIN_ENERGY_DENSITY
    const int TACS_ELEMENT_DISPLACEMENT
    const int TACS_ELEMENT_STRAIN
    const int TACS_ELEMENT_STRESS

cdef extern from "TACSElementVerification.h":
    int TacsTestElementBasisFunctions(TACSElementBasis*, double, int, double, double)
    int TacsTestElementBasisFaceNormals(TACSElementBasis*, double, int, double, double)
    int TacsTestElementBasisJacobianTransform(TACSElementBasis*, double, int, double, double)
    int TacsTestElementBasis(TACSElementBasis*, double, int, double, double)
    int TacsTestElementModelJacobian(TACSElementModel*, int, const double, double, int, double, double)
    int TacsTestElementModelAdjXptSensProduct(TACSElementModel*, int, const double, double, int, double, double)
    int TacsTestElementResidual(TACSElement*, int, double, const TacsScalar*,
                                const TacsScalar*, const TacsScalar*,
                                const TacsScalar*, double, int, double,
                                double)
    int TacsTestElementJacobian(TACSElement*, int, double, const TacsScalar*,
                                const TacsScalar*, const TacsScalar*,
                                const TacsScalar*, int, double, int, double,
                                double)
    int TacsTestAdjResProduct(TACSElement*, int, double, const TacsScalar*,
                              const TacsScalar*, const TacsScalar*, const TacsScalar*,
                              int, const TacsScalar*, double, int, double, double)
    int TacsTestAdjResXptProduct(TACSElement*, int, double, const TacsScalar*,
                                const TacsScalar*, const TacsScalar*,
                                const TacsScalar*, double, int, double,
                                double)
    int TacsTestElementMatDVSens(TACSElement*, ElementMatrixType, int, double, const TacsScalar*, const TacsScalar*,
                                 int, const TacsScalar*, double, int, double, double)
    int TacsTestElementMatSVSens(TACSElement*, ElementMatrixType, int, double, const TacsScalar*, const TacsScalar*,
                                 double, int, double, double)

cdef extern from "TACSTetrahedralBasis.h":
    cdef cppclass TACSLinearTetrahedralBasis(TACSElementBasis):
        TACSLinearTetrahedralBasis()

    cdef cppclass TACSQuadraticTetrahedralBasis(TACSElementBasis):
        TACSQuadraticTetrahedralBasis()

    cdef cppclass TACSCubicTetrahedralBasis(TACSElementBasis):
        TACSCubicTetrahedralBasis()

cdef extern from "TACSHexaBasis.h":
    cdef cppclass TACSLinearHexaBasis(TACSElementBasis):
        TACSLinearHexaBasis()

    cdef cppclass TACSQuadraticHexaBasis(TACSElementBasis):
        TACSQuadraticHexaBasis()

    cdef cppclass TACSCubicHexaBasis(TACSElementBasis):
        TACSCubicHexaBasis()

cdef extern from "TACSQuadBasis.h":
    cdef cppclass TACSLinearQuadBasis(TACSElementBasis):
        TACSLinearQuadBasis()

    cdef cppclass TACSQuadraticQuadBasis(TACSElementBasis):
        TACSQuadraticQuadBasis()

    cdef cppclass TACSCubicQuadBasis(TACSElementBasis):
        TACSCubicQuadBasis()

    cdef cppclass TACSQuarticQuadBasis(TACSElementBasis):
        TACSQuarticQuadBasis()

    cdef cppclass TACSQuinticQuadBasis(TACSElementBasis):
        TACSQuinticQuadBasis()

cdef extern from "TACSTriangularBasis.h":
    cdef cppclass TACSLinearTriangleBasis(TACSElementBasis):
        TACSLinearTriangleBasis()

    cdef cppclass TACSQuadraticTriangleBasis(TACSElementBasis):
        TACSQuadraticTriangleBasis()

    cdef cppclass TACSCubicTriangleBasis(TACSElementBasis):
        TACSCubicTriangleBasis()

cdef extern from "TACSHeatConduction.h":
    cdef cppclass TACSHeatConduction2D(TACSElementModel):
        TACSHeatConduction2D(TACSPlaneStressConstitutive*)

    cdef cppclass TACSHeatConduction3D(TACSElementModel):
        TACSHeatConduction3D(TACSSolidConstitutive*)

cdef extern from "TACSLinearElasticity.h":
    enum ElementStrainType:
        TACS_LINEAR_STRAIN
        TACS_NONLINEAR_STRAIN

    cdef cppclass TACSLinearElasticity2D(TACSElementModel):
        TACSLinearElasticity2D(TACSPlaneStressConstitutive*,
                               ElementStrainType)

    cdef cppclass TACSLinearElasticity3D(TACSElementModel):
        TACSLinearElasticity3D(TACSSolidConstitutive*,
                               ElementStrainType)
        TACSSolidConstitutive* getConstitutive()

cdef extern from "TACSThermoelasticity.h":
    const int TACS_STEADY_STATE_MECHANICAL
    const int TACS_STEADY_STATE_THERMAL

    cdef cppclass TACSLinearThermoelasticity2D(TACSElementModel):
        TACSLinearThermoelasticity2D(TACSPlaneStressConstitutive*,
                                     ElementStrainType, int)

    cdef cppclass TACSLinearThermoelasticity3D(TACSElementModel):
        TACSLinearThermoelasticity3D(TACSSolidConstitutive*,
                                     ElementStrainType, int)

cdef extern from "TACSPlateModel.h":
    cdef cppclass TACSPlateModel(TACSElementModel):
        TACSPlateModel(TACSShellConstitutive*)

cdef extern from "TACSThermoelasticPlateModel.h":
    cdef cppclass TACSThermoelasticPlateModel(TACSElementModel):
        TACSThermoelasticPlateModel(TACSShellConstitutive*)

cdef extern from "TACSElement2D.h":
    cdef cppclass TACSElement2D(TACSElement):
        TACSElement2D(TACSElementModel*, TACSElementBasis*)

cdef extern from "TACSElement3D.h":
    cdef cppclass TACSElement3D(TACSElement):
        TACSElement3D(TACSElementModel*, TACSElementBasis*)

cdef extern from "TACSTraction2D.h":
    cdef cppclass TACSTraction2D(TACSElement):
        TACSTraction2D(int, int, TACSElementBasis*, TacsScalar*, int)

cdef extern from "TACSTraction3D.h":
    cdef cppclass TACSTraction3D(TACSElement):
        TACSTraction3D(int, int, TACSElementBasis*, TacsScalar*, int)

cdef extern from "TACSConvectiveTraction2D.h":
    cdef cppclass TACSConvectiveTraction2D(TACSElement):
        TACSConvectiveTraction2D(int, int, int, TacsScalar,
                                 TacsScalar,
                                 TACSElementBasis*)

cdef extern from "TACSConvectiveTraction3D.h":
    cdef cppclass TACSConvectiveTraction3D(TACSElement):
        TACSConvectiveTraction3D(int, int, int, TacsScalar,
                                 TacsScalar,
                                 TACSElementBasis*)

cdef extern from "TACSShellElementTransform.h":
    cdef cppclass TACSShellTransform(TACSObject):
        pass

    cdef cppclass TACSShellNaturalTransform(TACSShellTransform):
        TACSShellNaturalTransform()

    cdef cppclass TACSShellRefAxisTransform(TACSShellTransform):
        TACSShellRefAxisTransform(TacsScalar*)

cdef extern from "TACSShellElementDefs.h":
    cdef cppclass TACSQuad4Shell(TACSElement):
        TACSQuad4Shell(TACSShellTransform*,
                       TACSShellConstitutive*)

    cdef cppclass TACSQuad9Shell(TACSElement):
        TACSQuad9Shell(TACSShellTransform*,
                       TACSShellConstitutive*)

    cdef cppclass TACSQuad16Shell(TACSElement):
        TACSQuad16Shell(TACSShellTransform*,
                       TACSShellConstitutive*)

    cdef cppclass TACSTri3Shell(TACSElement):
        TACSTri3Shell(TACSShellTransform*,
                      TACSShellConstitutive*)

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

cdef extern from "TACSKinematicConstraints.h":
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
        TACSRevoluteConstraint( int fixed_ref_point,
                                TACSGibbsVector *point,
                                TACSGibbsVector *eAVec,
                                int inertial_rev_axis )

    cdef cppclass TACSRigidLink(TACSElement):
        TACSRigidLink(TACSRigidBody*)

    cdef cppclass TACSRevoluteDriver(TACSElement):
        TACSRevoluteDriver(TACSGibbsVector*, TacsScalar)

    cdef cppclass TACSMotionDriver(TACSElement):
        TACSMotionDriver(TACSGibbsVector*, TacsScalar, int)

    cdef cppclass TACSAverageConstraint(TACSElement):
        TACSAverageConstraint(TACSRigidBody*, TACSGibbsVector*,
                              TACSRefFrame*, int)

cdef extern from "RBE2.h":
    cdef cppclass RBE2(TACSElement):
        RBE2(int, int*)
        void setScalingParameters(double, double)

cdef extern from "RBE3.h":
    cdef cppclass RBE3(TACSElement):
        RBE3(int, int*, double*, int*)
        void setScalingParameters(double, double)

cdef extern from  "MITC3.h":
    cdef cppclass MITC3(TACSElement):
        MITC3(TACSTimoshenkoConstitutive *_stiff,
              TACSGibbsVector *_gravity,
              TACSGibbsVector *_vInit,
              TACSGibbsVector *_omegaInit)

cdef extern from "TACSElementWrapper.h":
    cdef cppclass TACSElementWrapper(TACSElement):
        TACSElementWrapper(PyObject*, int, int)

        int (*getmultiplierindex)(void*);
        void (*getinitconditions)(void*, int, int, const TacsScalar*, int,
                                  TacsScalar*, TacsScalar*, TacsScalar*)
        void (*addresidual)(void*, int, double, int, const TacsScalar*,
                            int, const TacsScalar*, const TacsScalar*,
                            const TacsScalar*, TacsScalar*)
        void (*addjacobian)(void*, int, double, TacsScalar, TacsScalar, TacsScalar,
                            int, const TacsScalar*, int, const TacsScalar*,
                            const TacsScalar*, const TacsScalar*,
                            TacsScalar*, TacsScalar*)
